import type { LocalBundle } from "./classify-bundle";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { gzip, ungzip } from "pako";

import { partFileMap } from "./classify-bundle";

import { uploadChunkedToS3 } from "@/lib/upload/uploadChunkedToS3";
import { uploadImageToS3 } from "@/lib/upload/uploadImageToS3";
import { uploadZarrFolderToS3 } from "@/lib/upload/uploadZarrFolderToS3";

/**
 * Save a folder that was dropped for a browser preview: every part goes to
 * S3 exactly as it is on disk — nothing is converted — each as its own
 * dataset, and the cells are told to draw the molecules over themselves and
 * the image under, so the saved link opens the way the preview looked.
 */

export interface SavedBundle {
  cellsId?: string;
  moleculesId?: string;
  imageId?: string;
  /** Where to open what was saved. */
  viewerUrl: string;
}

type Progress = (fraction: number, message: string) => void;

async function sha256Hex(text: string): Promise<string> {
  const digest = await crypto.subtle.digest(
    "SHA-256",
    new TextEncoder().encode(text),
  );

  return Array.from(new Uint8Array(digest), (b) =>
    b.toString(16).padStart(2, "0"),
  ).join("");
}

async function postJson(url: string, body: unknown, method = "POST") {
  const res = await fetch(url, {
    method,
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  const json = await res.json().catch(() => ({}));

  if (!res.ok) {
    throw new Error(json.error || json.message || `${url} failed`);
  }

  return json;
}

const bytesOf = (files: Map<string, File>) =>
  Array.from(files.values()).reduce((sum, f) => sum + f.size, 0);

/** A merfisheyes cell folder, file for file. */
async function saveChunkedCells(
  files: Map<string, File>,
  title: string,
  onProgress: Progress,
): Promise<string> {
  const manifestText = await files.get("manifest.json")!.text();
  const manifest = JSON.parse(manifestText);
  // Only the name it was saved under changes
  const named = new Blob([JSON.stringify({ ...manifest, name: title })], {
    type: "application/json",
  });
  const { datasetId, uploadId } = await uploadChunkedToS3({
    fingerprint: await sha256Hex(manifestText),
    metadata: {
      title,
      numCells: manifest.statistics.total_cells,
      numGenes: manifest.statistics.total_genes,
      platform: manifest.type,
    },
    files: Array.from(files, ([key, file]) =>
      key === "manifest.json"
        ? {
            key,
            blob: named,
            size: named.size,
            contentType: "application/json",
          }
        : {
            key,
            blob: file,
            size: file.size,
            contentType: file.type || "application/octet-stream",
          },
    ),
    onProgress: (percent, message) => onProgress(percent / 100, message),
  });

  await postJson(`/api/datasets/${datasetId}/complete`, { uploadId });

  return datasetId;
}

/** A SpatialData store, with the annotations folder beside it if there is one. */
async function saveSpatialDataStore(
  bundle: LocalBundle,
  root: string,
  cells: StandardizedDataset,
  title: string,
  onProgress: Progress,
): Promise<string> {
  // Keys under the dataset: the store keeps its folder name ("store.zarr/…")
  // and its sidecar sits beside it, the same layout the viewer reads from S3.
  const folder = root || `${bundle.name.replace(/\.zarr$/, "")}.zarr`;
  const files = Array.from(partFileMap(bundle, root), ([key, file]) => ({
    key: `${folder}/${key}`,
    file,
  }));

  if (root) {
    for (const [key, file] of bundle.files) {
      if (key.startsWith("annotations/")) files.push({ key, file });
    }
  }

  const { datasetId } = await uploadZarrFolderToS3({
    title,
    platform: "spatialdata",
    numCells: cells.getPointCount(),
    numGenes: cells.genes.length,
    files,
    lastKey: `${folder}/zarr.json`,
    onProgress: (p) =>
      onProgress(
        p.loaded / Math.max(p.total, 1),
        `Uploaded ${p.fileIndex}/${p.fileCount} files`,
      ),
  });

  return datasetId;
}

/** A merfisheyes molecule folder, file for file. */
async function saveMolecules(
  files: Map<string, File>,
  title: string,
  onProgress: Progress,
): Promise<string> {
  const upload = new Map(files);
  const plain = upload.get("manifest.json");
  const zipped = upload.get("manifest.json.gz");
  const manifestText = zipped
    ? ungzip(new Uint8Array(await zipped.arrayBuffer()), { to: "string" })
    : await plain!.text();
  const manifest = JSON.parse(manifestText);

  // The viewer fetches the gzipped manifest; a plain one is stored that way
  if (plain && !zipped) {
    upload.delete("manifest.json");
    upload.set(
      "manifest.json.gz",
      new File([gzip(manifestText)], "manifest.json.gz", {
        type: "application/gzip",
      }),
    );
  }

  const contentType = (key: string) =>
    key.endsWith(".gz") ? "application/gzip" : "application/octet-stream";
  const { datasetId, uploadId, uploadUrls } = await postJson(
    "/api/single-molecule/initiate",
    {
      fingerprint: await sha256Hex(manifestText),
      metadata: {
        title,
        numMolecules: manifest?.statistics?.total_molecules,
        numGenes:
          manifest?.statistics?.unique_genes ??
          manifest?.genes?.unique_gene_names?.length,
        platform: manifest?.type,
        description: "",
      },
      manifest,
      files: Array.from(upload, ([key, file]) => ({
        key,
        size: file.size,
        contentType: contentType(key),
      })),
    },
  );
  const total = bytesOf(upload);
  let loaded = 0;
  let done = 0;

  for (const [key, file] of upload) {
    const put = await fetch(uploadUrls[key], {
      method: "PUT",
      body: file,
      headers: { "Content-Type": contentType(key) },
    });

    if (!put.ok) throw new Error(`Upload of ${key} failed (${put.status})`);
    await postJson(
      `/api/single-molecule/${datasetId}/files/${encodeURIComponent(key)}/complete`,
      { uploadId },
    );
    loaded += file.size;
    onProgress(
      loaded / Math.max(total, 1),
      `Uploaded ${++done}/${upload.size} files`,
    );
  }
  await postJson(`/api/single-molecule/${datasetId}/complete`, { uploadId });

  return datasetId;
}

export async function saveLocalBundle({
  bundle,
  name,
  cells,
  onProgress,
}: {
  bundle: LocalBundle;
  name: string;
  /** The opened cells, when the folder has any (for a store's counts). */
  cells?: StandardizedDataset | null;
  onProgress: Progress;
}): Promise<SavedBundle> {
  const imageFiles =
    bundle.image !== undefined ? partFileMap(bundle, bundle.image) : null;
  const moleculeFiles =
    bundle.molecules !== undefined
      ? partFileMap(bundle, bundle.molecules)
      : null;
  const cellFiles = bundle.cells
    ? partFileMap(bundle, bundle.cells.root)
    : null;
  // One bar across the parts, each taking its share of the bytes
  const sizes = [imageFiles, moleculeFiles, cellFiles].map((f) =>
    f ? bytesOf(f) : 0,
  );
  const total = Math.max(sizes[0] + sizes[1] + sizes[2], 1);
  const stage = (index: number, label: string): Progress => {
    const before = sizes.slice(0, index).reduce((a, b) => a + b, 0);

    return (fraction, message) =>
      onProgress(
        (before + fraction * sizes[index]) / total,
        `${label}: ${message}`,
      );
  };
  const parts =
    (imageFiles ? 1 : 0) + (moleculeFiles ? 1 : 0) + (cellFiles ? 1 : 0);
  const titled = (kind: string) => (parts > 1 ? `${name} (${kind})` : name);
  const saved: SavedBundle = { viewerUrl: "" };

  if (imageFiles) {
    saved.imageId = (
      await uploadImageToS3({
        title: titled("image"),
        files: Array.from(imageFiles, ([key, file]) => ({ key, file })),
        onProgress: (p) =>
          stage(0, "Image")(
            p.loaded / Math.max(p.total, 1),
            `uploaded ${p.fileIndex}/${p.fileCount} files`,
          ),
      })
    ).datasetId;
  }
  if (moleculeFiles) {
    saved.moleculesId = await saveMolecules(
      moleculeFiles,
      titled("molecules"),
      stage(1, "Molecules"),
    );
  }
  if (bundle.cells && cellFiles) {
    saved.cellsId =
      bundle.cells.kind === "chunked"
        ? await saveChunkedCells(cellFiles, titled("cells"), stage(2, "Cells"))
        : await saveSpatialDataStore(
            bundle,
            bundle.cells.root,
            cells!,
            titled("cells"),
            stage(2, "Cells"),
          );
  }

  // Link them, and group them in one project
  onProgress(1, "Linking the datasets...");
  let projectId: string | undefined;

  if (saved.cellsId && saved.moleculesId) {
    projectId = (
      await postJson(`/api/datasets/${saved.cellsId}/overlay`, {
        smDatasetId: saved.moleculesId,
        force: true,
      })
    ).projectId;
  }
  if (saved.cellsId && saved.imageId) {
    await postJson(
      `/api/datasets/${saved.cellsId}/overlay`,
      { imageDatasetId: saved.imageId },
      "PATCH",
    );
  }
  const ids = [saved.cellsId, saved.moleculesId, saved.imageId].filter(
    (id) => id !== undefined,
  );

  if (ids.length > 1) {
    projectId ??= (await postJson("/api/projects", { title: name })).id;
    // Adding one that is already in the project is refused; that is fine
    await Promise.all(
      ids.map((datasetId) =>
        postJson(`/api/projects/${projectId}/datasets`, { datasetId }).catch(
          () => {},
        ),
      ),
    );
  }

  saved.viewerUrl = saved.cellsId
    ? `${bundle.cells!.kind === "spatialdata" ? "/spatialdata-viewer" : "/viewer"}/${saved.cellsId}`
    : saved.moleculesId
      ? `/sm-viewer/${saved.moleculesId}`
      : `/image-viewer/${saved.imageId}`;

  return saved;
}
