import { ungzip } from "pako";

import { resolveZarrUrl } from "@/lib/spatialdata/resolveUrl";
import { parseS3Url, testManifestAccess } from "@/lib/utils/s3-url-parser";

/**
 * Which viewer opens what a URL holds, and its path: a merfisheyes folder
 * (manifest.json[.gz]) in the cell or molecule viewer by what the manifest
 * says; otherwise a SpatialData store (or a folder holding one) or an
 * OME-Zarr image. The one implementation behind the "Load from URL" dialog
 * and `/open?url=` (which `merfisheyes-mcp view` opens for a local folder).
 * Throws an Error whose message says what to check.
 */
export async function viewerPathForUrl(
  input: string,
  opts: {
    /** The homepage's mode toggle: decides only an ambiguous manifest. */
    datasetType?: "single_cell" | "single_molecule";
    /** Progress, for a dialog to show. */
    onStep?: (step: string) => void;
  } = {},
): Promise<string> {
  const step = opts.onStep ?? (() => {});

  step("Validating URL...");
  const parsed = parseS3Url(input);

  step("Testing manifest access...");
  // Test if manifest exists (tries .gz first, then .json)
  const manifestTest = await testManifestAccess(parsed.baseUrl);

  if (!manifestTest.exists) {
    // Not a merfisheyes folder: a SpatialData store, an OME-Zarr image,
    // or a folder holding a store. Each is read as it is by its viewer.
    step("Looking for a SpatialData store or an image...");
    const target = await resolveZarrUrl(parsed.baseUrl);

    if (!target) {
      throw new Error(
        `Nothing to open at that URL. Please ensure:\n` +
          `1. The URL is a merfisheyes folder (manifest.json or manifest.json.gz), a SpatialData store or an OME-Zarr image (zarr.json), or a folder holding a store\n` +
          `2. A folder names its store in mapping.json ("cells"), or the store is named after the folder (sample/sample.zarr)\n` +
          `3. The S3 bucket/folder has public read access\n` +
          `4. CORS is configured to allow requests from ${window.location.origin}`,
      );
    }

    // A store's link keeps the URL as entered: the viewer finds the
    // store in a folder the same way.
    return target.kind === "image"
      ? `/image-viewer/from-s3?url=${encodeURIComponent(target.url)}`
      : `/spatialdata-viewer/from-s3?url=${encodeURIComponent(parsed.baseUrl)}`;
  }

  step("Downloading and parsing manifest...");
  const response = await fetch(manifestTest.url);

  if (!response.ok) {
    throw new Error(`Failed to download manifest: ${response.statusText}`);
  }

  const manifest = JSON.parse(
    manifestTest.url.endsWith(".gz")
      ? ungzip(new Uint8Array(await response.arrayBuffer()), { to: "string" })
      : await response.text(),
  );

  // Always detect the type from the manifest's own fields — that is the
  // reliable signal. The `datasetType` toggle must NOT override a manifest
  // that clearly says otherwise: pasting a single-molecule URL while the
  // toggle sits on "single cell" used to force the cell viewer, which then
  // crashed on a molecule manifest. It decides only an ambiguous manifest.
  //
  // Reliable, mutually-exclusive fields (see scripts/process_*.py):
  //   SM manifest → genes.unique_gene_names, statistics.total_molecules
  //   SC manifest → statistics.total_cells, files.expression_chunks
  // Do NOT key off `manifest.type`: both pipelines write the *format*
  // there ("h5ad" | "xenium" | "merscope"), so "merscope"/"xenium" appear
  // on both kinds and can't discriminate.
  const looksSingleMolecule =
    Boolean(manifest.genes?.unique_gene_names) ||
    manifest.statistics?.total_molecules != null;
  const looksSingleCell =
    manifest.statistics?.total_cells != null ||
    manifest.files?.expression_chunks != null;

  let detected: "single_cell" | "single_molecule";

  if (looksSingleMolecule && !looksSingleCell) detected = "single_molecule";
  else if (looksSingleCell && !looksSingleMolecule) detected = "single_cell";
  else if (opts.datasetType) detected = opts.datasetType;
  else {
    throw new Error(
      "Could not determine from the manifest whether this is a single-cell " +
        "or single-molecule dataset.",
    );
  }

  // The toggle said one thing but the manifest is another: route by the
  // manifest, and say so, so the other viewer is not a surprise.
  if (opts.datasetType && detected !== opts.datasetType) {
    step(
      `This is a ${
        detected === "single_molecule" ? "single-molecule" : "single-cell"
      } dataset — opening the matching viewer…`,
    );
  }

  step("Loading dataset...");
  const encoded = encodeURIComponent(parsed.baseUrl);

  return detected === "single_molecule"
    ? `/sm-viewer/from-s3?url=${encoded}`
    : `/viewer/from-s3?url=${encoded}`;
}
