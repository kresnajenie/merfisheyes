"use client";

import type { ReactNode } from "react";

import { useCallback, useRef, useState } from "react";
import { useRouter } from "next/navigation";
import { useSession } from "next-auth/react";
import { toast } from "react-toastify";
import * as Comlink from "comlink";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import { Input } from "@heroui/input";
import { Button } from "@heroui/button";
import { Checkbox } from "@heroui/checkbox";

import { StandardizedDataset } from "@/lib/StandardizedDataset";
import { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import {
  MoleculeDatasetType,
  MoleculeColumnMapping,
  MOLECULE_COLUMN_MAPPINGS,
} from "@/lib/config/moleculeColumnMappings";
import { createPortal } from "react-dom";
import { Spinner } from "@heroui/react";
import { readMoleculePreview } from "@/lib/services/molecule-preview";
import { SingleMoleculeHeaderModal } from "@/components/single-molecule-header-modal";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { getSingleMoleculeWorker } from "@/lib/workers/singleMoleculeWorkerManager";
import { resetHooks, markDatasetLoaded } from "@/lib/utils/test-hooks";
import {
  classifyFolder,
  fallbackSingleMoleculeFile,
} from "@/lib/ingest/classify-folder";
import { useUploadStore } from "@/lib/stores/uploadStore";
import { AuthModal } from "@/components/auth-modal";

type RawFile = { key: string; file: File; contentType: string };

/** A server upload prepared and waiting for the user to name it. */
interface PendingServerUpload {
  kind: "single_cell" | "single_molecule";
  files: RawFile[];
  processingParams: Record<string, unknown>;
  /**
   * Set when a single-cell folder also contained single-molecule data. Offers
   * to upload it as a separate SM dataset and auto-create a project linking the
   * two (with a mapping.json overlay).
   */
  sm?: { files: RawFile[]; fileName: string };
}

/**
 * Reserved raw key for the user-supplied cell-type CSV.
 *
 * The worker moves this file out of the raw input tree before running the
 * processor and passes it as --mmc-csv. It is named explicitly in
 * processingParams rather than detected by extension, because a MERSCOPE
 * upload consists entirely of CSVs.
 */
const ANNOTATION_CSV_KEY = "__annotations__/mapping.csv";

// "folder" is a meta-type that auto-detects the dropped folder shape
// (zarr / chunked / xenium / merscope / h5ad-inside) and dispatches to the
// appropriate handler.
//
// "file" (single-molecule only) is the parquet/csv equivalent: sniff the
// dropped file's column names and pick xenium / merscope / custom on the fly.
type UploadType =
  | "h5ad"
  | "h5ad-zarr"
  | "xenium"
  | "merscope"
  | "chunked"
  | "folder"
  | "file";

// Map UploadType to MoleculeDatasetType for single molecule datasets
const UPLOAD_TYPE_TO_PARQUET_TYPE: Record<UploadType, MoleculeDatasetType> = {
  h5ad: "custom",
  "h5ad-zarr": "custom",
  xenium: "xenium",
  merscope: "merscope",
  chunked: "custom",
  folder: "custom",
  file: "custom",
};

interface FileUploadProps {
  type: UploadType;
  title: string;
  description: string;
  singleMolecule?: boolean;
  /**
   * When true, this card uploads the raw bytes to the server-side ingestion
   * pipeline (no browser parsing) instead of parsing locally and navigating to
   * the viewer. Requires an authenticated session.
   */
  serverUpload?: boolean;
  /**
   * Optional per-cell label CSV uploaded alongside the dataset (single-cell
   * only). Merged during chunking, so its columns get palettes and DE stats
   * like any other cluster column.
   */
  annotationCsv?: File | null;
  /**
   * Optional row of small icons rendered under the description. Used by the
   * unified "Folder" card to advertise which formats it accepts.
   */
  icons?: ReactNode;
}

export function FileUpload({
  type,
  title,
  description,
  singleMolecule = false,
  serverUpload = false,
  annotationCsv,
  icons,
}: FileUploadProps) {
  const [isDragging, setIsDragging] = useState(false);
  const [progress, setProgress] = useState(0);
  const [progressMessage, setProgressMessage] = useState("");

  const { data: session } = useSession();

  // Server-side raw upload runs in a global store so the transfer survives the
  // navigation to /explore, where the top bar shows its progress.
  const startUpload = useUploadStore((s) => s.start);

  // Auth-on-drop: a signed-out user can drop a file; we stash it and pop the
  // sign-in modal, then resume the upload on success (no sign-in wall first).
  const [authModalOpen, setAuthModalOpen] = useState(false);
  const pendingFilesRef = useRef<File[] | null>(null);

  // Single-molecule "confirm your columns" step. When a raw parquet/CSV is
  // dropped we read a preview and pause here until the user confirms/remaps the
  // gene / x / y / z / cell-id columns.
  const [headerModal, setHeaderModal] = useState<{
    fileName: string;
    columns: string[];
    rows: Record<string, unknown>[];
    autoMapping: MoleculeColumnMapping;
  } | null>(null);
  const mappingResolverRef = useRef<
    ((m: MoleculeColumnMapping | null) => void) | null
  >(null);
  // Reading a parquet/CSV preview to detect columns can take a few seconds on a
  // large file; show an explicit indicator so it does not look frozen.
  const [columnPreviewLoading, setColumnPreviewLoading] = useState(false);
  // Files whose in-browser processing failed — offer to process them on the
  // server instead (same flow as Upload & Save), reusing the dropped file.
  const [serverFallbackFiles, setServerFallbackFiles] = useState<File[] | null>(null);

  // Read a preview, open the confirm modal, and resolve with the user's chosen
  // mapping (or null if they cancel).
  const awaitColumnMapping = useCallback(
    async (file: File): Promise<MoleculeColumnMapping | null> => {
      setColumnPreviewLoading(true);
      let preview;
      try {
        preview = await readMoleculePreview(file);
      } finally {
        setColumnPreviewLoading(false);
      }
      const autoMapping = MOLECULE_COLUMN_MAPPINGS[preview.autoType];

      return new Promise<MoleculeColumnMapping | null>((resolve) => {
        mappingResolverRef.current = resolve;
        setHeaderModal({
          fileName: file.name,
          columns: preview.columns,
          rows: preview.rows,
          autoMapping,
        });
      });
    },
    [],
  );

  const resolveMapping = (mapping: MoleculeColumnMapping | null) => {
    setHeaderModal(null);
    mappingResolverRef.current?.(mapping);
    mappingResolverRef.current = null;
  };

  // A prepared server upload awaiting a title. Auto-derived names collide (every
  // "transcripts.csv" would share one title), so the user confirms/renames here
  // before the transfer starts.
  const [pendingUpload, setPendingUpload] =
    useState<PendingServerUpload | null>(null);
  const [titleDraft, setTitleDraft] = useState("");
  // Whether to also upload the detected single-molecule data (combined flow).
  const [includeSm, setIncludeSm] = useState(true);

  // Use appropriate store based on singleMolecule mode
  const cellStore = useDatasetStore();
  const smStore = useSingleMoleculeStore();

  const { setLoading, setError } = singleMolecule ? smStore : cellStore;
  const router = useRouter();

  const inputId = `${singleMolecule ? "sm" : "sc"}-file-input-${type}`;

  const handleDragOver = useCallback((e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(true);
  }, []);

  const handleDragLeave = useCallback((e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
  }, []);

  // NOTE: plain functions (not useCallback) so they always call the current
  // handleFiles closure. Memoizing them would capture a stale handleFiles and
  // miss prop changes like `serverUpload` toggling on.
  const handleDrop = async (e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);

    // dataTransfer.files does NOT recurse into dropped folders — for that
    // we have to walk dataTransfer.items via webkitGetAsEntry. We patch
    // webkitRelativePath on each yielded File so the detection helpers and
    // downstream upload code (which all key off webkitRelativePath) work
    // identically to the click-to-pick path.
    let files: File[] = [];

    if (e.dataTransfer.items && e.dataTransfer.items.length > 0) {
      // Capture entries synchronously; FileSystemEntry handles become
      // invalid after the drop event finishes propagating.
      const entries = Array.from(e.dataTransfer.items)
        .filter((item) => item.kind === "file")
        .map((item) => (item as any).webkitGetAsEntry?.())
        .filter(Boolean);

      try {
        files = await collectFilesFromEntries(entries);
      } catch (err) {
        console.error("[FileUpload] folder walk failed:", err);
        files = Array.from(e.dataTransfer.files);
      }
    } else {
      files = Array.from(e.dataTransfer.files);
    }

    handleFiles(files);
  };

  const handleFileInput = (e: React.ChangeEvent<HTMLInputElement>) => {
    const files = e.target.files ? Array.from(e.target.files) : [];

    handleFiles(files);
  };

  /**
   * Detect if files represent a chunked dataset folder (from Python script)
   */
  const isChunkedDatasetFolder = (files: File[]): boolean => {
    const fileNames = files.map((f) => {
      const path = f.webkitRelativePath || f.name;
      // Remove root folder name
      const parts = path.split("/");

      return parts.length > 1 ? parts.slice(1).join("/") : path;
    });

    // Check for required chunked dataset files
    const hasManifest = fileNames.some((name) => name === "manifest.json");
    const hasExprIndex = fileNames.some((name) => name === "expr/index.json");
    const hasChunks = fileNames.some((name) =>
      /^expr\/chunk_\d{5}\.bin\.gz$/.test(name),
    );
    const hasSpatial = fileNames.some(
      (name) => name === "coords/spatial.bin.gz",
    );

    return hasManifest && hasExprIndex && hasChunks && hasSpatial;
  };

  /**
   * Detect if a folder looks like an h5ad-zarr store (i.e. contains a top-level
   * .zgroup / .zarray / zarr.json marker).
   */
  const isZarrFolder = (files: File[]): boolean => {
    return files.some((f) => {
      const path = f.webkitRelativePath || f.name;
      const parts = path.split("/");
      const rel = parts.length > 1 ? parts.slice(1).join("/") : path;

      return rel === ".zgroup" || rel === ".zarray" || rel === "zarr.json";
    });
  };

  /** Detect Xenium folder by canonical files (cells.csv[.gz]). */
  const isXeniumFolder = (files: File[]): boolean => {
    return files.some((f) => {
      const name = (f.name || "").toLowerCase();

      return name === "cells.csv" || name === "cells.csv.gz";
    });
  };

  /** Detect MERSCOPE folder by canonical files (cell_metadata.csv). */
  const isMerscopeFolder = (files: File[]): boolean => {
    return files.some((f) => {
      const name = (f.name || "").toLowerCase();

      return name === "cell_metadata.csv";
    });
  };

  /** Detect a single .h5ad file living inside a dropped folder. */
  const findH5adFile = (files: File[]): File | null => {
    return (
      files.find((f) => {
        const name = (f.name || "").toLowerCase();
        const rel = (f.webkitRelativePath || "").toLowerCase();

        return name.endsWith(".h5ad") || rel.endsWith(".h5ad");
      }) ?? null
    );
  };

  /**
   * Auto-detect the format of a dropped folder. Order matters — zarr's
   * .zgroup is unambiguous; chunked's manifest cluster is unambiguous;
   * xenium / merscope have format-specific files. h5ad-inside is the catch
   * for "folder containing an .h5ad" (matches existing MERSCOPE fallback).
   */
  const detectFolderFormat = (
    files: File[],
  ): "h5ad-zarr" | "chunked" | "xenium" | "merscope" | "h5ad" | null => {
    if (isZarrFolder(files)) return "h5ad-zarr";
    if (isChunkedDatasetFolder(files)) return "chunked";
    if (isXeniumFolder(files)) return "xenium";
    if (isMerscopeFolder(files)) return "merscope";
    if (findH5adFile(files)) return "h5ad";

    return null;
  };

  /**
   * Detect if files represent a chunked single molecule dataset folder (from Python script)
   */
  const isChunkedSingleMoleculeFolder = (files: File[]): boolean => {
    const fileNames = files.map((f) => {
      const path = f.webkitRelativePath || f.name;
      // Remove root folder name
      const parts = path.split("/");

      return parts.length > 1 ? parts.slice(1).join("/") : path;
    });

    // Check for required chunked single molecule files
    // Accept both compressed (manifest.json.gz) and uncompressed (manifest.json) manifests
    const hasManifest = fileNames.some(
      (name) => name === "manifest.json.gz" || name === "manifest.json",
    );
    const hasGenesFolder = fileNames.some((name) => name.startsWith("genes/"));
    const hasGeneBinFiles = fileNames.some((name) =>
      /^genes\/[^\/]+\.bin\.gz$/.test(name),
    );

    return hasManifest && hasGenesFolder && hasGeneBinFiles;
  };

  // Server-side ingestion entry point. If the user isn't signed in yet, stash
  // the dropped files and open the auth modal — the upload resumes on success,
  // so there's no sign-in wall before the drop.
  const handleServerUpload = async (files: File[]) => {
    if (files.length === 0) return;

    if (!session?.user) {
      pendingFilesRef.current = files;
      setAuthModalOpen(true);

      return;
    }

    await doServerUpload(files);
  };

  // The actual raw-bytes upload. Assumes an authenticated session — the entry
  // point above guarantees it, and the ingest API re-checks via cookie.
  const doServerUpload = async (files: File[]) => {
    const kind = singleMolecule ? "single_molecule" : "single_cell";
    const title = files[0].webkitRelativePath?.split("/")[0] || files[0].name;

    // A picked folder is filtered down to the files the processor actually
    // reads — a Xenium export is mostly images and QC output, and uploading it
    // whole means sending gigabytes nothing will open.
    const classification = classifyFolder(files);
    let selected =
      kind === "single_molecule"
        ? classification.singleMolecule
        : classification.singleCell;

    // Explicit single-molecule upload with no transcripts-style name (e.g.
    // "spots_set4.parquet"): the worker accepts any lone .csv/.parquet, so a
    // single candidate is unambiguous — mirror that instead of erroring.
    if (!selected && kind === "single_molecule") {
      selected = fallbackSingleMoleculeFile(files) ?? undefined;
    }

    if (!selected) {
      toast.error(
        kind === "single_molecule"
          ? "No single-molecule data found — pick one .parquet or .csv transcripts file."
          : "No single-cell data found in that folder.",
      );

      return;
    }

    const selectedKeys = new Set(selected.files.map((f) => f.key));
    const ignored = classification.ignored.filter(
      (e) => !selectedKeys.has(e.key),
    );

    if (ignored.length > 0) {
      const skippedBytes = ignored.reduce((sum, f) => sum + f.file.size, 0);

      toast(
        `Uploading ${selected.files.length} of ${files.length} files — skipping ${ignored.length} the processor doesn't read (${(skippedBytes / 1024 ** 2).toFixed(0)} MB).`,
      );
    }

    const rawFiles = selected.files.map(({ key, file }) => ({
      key,
      file,
      contentType: "application/octet-stream",
    }));

    // The user's cell-type CSV rides along under a reserved key. Named
    // explicitly in processingParams rather than sniffed for, because a
    // MERSCOPE upload is nothing but CSVs.
    if (annotationCsv && kind === "single_cell") {
      rawFiles.push({
        key: ANNOTATION_CSV_KEY,
        file: annotationCsv,
        contentType: "text/csv",
      });
    }

    // Hand the transfer to the global upload store (so it survives navigation)
    // and go straight to /explore, where the top bar shows byte progress and
    // "Your uploads" shows the dataset processing once the bytes are up.
    //
    // v1 pipeline: chunk only. chunkSize 1 = one file per gene, so a gene
    // selection in the viewer fetches only that gene instead of a whole chunk.
    const processingParams = {
      kind,
      stages: {
        chunk: {
          chunkSize: 1,
          ...(annotationCsv && kind === "single_cell"
            ? { mmcCsv: ANNOTATION_CSV_KEY }
            : {}),
        },
      },
    };

    // A single-cell folder that also contains single-molecule data (a
    // transcripts file) offers a combined upload: a separate SM dataset + an
    // auto-created project linking the two.
    const sm =
      kind === "single_cell" && classification.singleMolecule
        ? {
            files: classification.singleMolecule.files.map(({ key, file }) => ({
              key,
              file,
              contentType: "application/octet-stream",
            })),
            fileName: classification.singleMolecule.files[0].key,
          }
        : undefined;

    // Stage the upload and let the user confirm/rename the title before the
    // transfer begins (see startPreparedUpload). Navigation to /explore happens
    // on confirm, once the transfer is actually kicked off.
    setPendingUpload({ kind, files: rawFiles, processingParams, sm });
    setTitleDraft(title);
    setIncludeSm(true);
  };

  // Best-effort: create a project grouping the SC + SM datasets from a combined
  // upload. Runs once both uploads are initiated (their rows exist).
  const createCombinedProject = async (
    title: string,
    scId: string,
    smId: string,
  ) => {
    try {
      const res = await fetch("/api/projects", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ title }),
      });

      if (!res.ok) return;
      const project = await res.json();

      await Promise.all(
        [scId, smId].map((datasetId) =>
          fetch(`/api/projects/${project.id}/datasets`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ datasetId }),
          }),
        ),
      );
    } catch {
      // Non-fatal — the datasets still upload; only the auto-grouping is missed.
    }
  };

  // Kick off the staged upload with the (possibly edited) title, then go to
  // /explore where the top bar tracks byte + processing progress.
  const startPreparedUpload = async () => {
    if (!pendingUpload) return;
    const pending = pendingUpload;
    const base = titleDraft.trim() || "Untitled Dataset";

    // A single-molecule file is being uploaded for server processing — either a
    // standalone SM upload, or the SM half of a combined SC+SM upload. Confirm
    // its columns first so the worker gets an explicit gene/x/y/z/cell-id
    // mapping instead of auto-detecting.
    const smFile =
      pending.kind === "single_molecule"
        ? pending.files[0]?.file
        : pending.sm && includeSm
          ? pending.sm.files[0]?.file
          : undefined;

    // Close the name dialog first so only the column modal is shown.
    setPendingUpload(null);

    let smMapping: MoleculeColumnMapping | undefined;

    if (smFile) {
      try {
        const mapping = await awaitColumnMapping(smFile);

        if (!mapping) return; // user cancelled the column confirm — abort
        smMapping = mapping;
      } catch (err) {
        // Reading the file's columns failed (unreadable parquet encoding, etc.).
        // Don't strand the upload — fall back to server-side auto-detection and
        // tell the user what happened instead of silently aborting.
        console.error("[FileUpload] Column preview failed:", err);
        toast.error(
          `Couldn't read columns from ${smFile.name} — uploading with automatic column detection.`,
        );
      }
    }

    // Combined SC + SM upload: two separate datasets (each tracked as its own
    // bar) + an auto-created project. The SC upload carries linkedSmDatasetId so
    // the server writes a mapping.json (overlay) when it finishes processing.
    if (pending.sm && includeSm) {
      // Start the SM upload; once it's initiated (id known), start the SC upload
      // linked to it, and once the SC is initiated too, group both in a project.
      void startUpload({
        kind: "single_molecule",
        title: `${base} (molecules)`,
        files: pending.sm.files,
        processingParams: {
          kind: "single_molecule",
          ...(smMapping ? { columnMapping: smMapping } : {}),
          stages: { chunk: { chunkSize: 1 } },
        },
        onInitiated: (smId) => {
          void startUpload({
            kind: pending.kind,
            title: `${base} (cells)`,
            files: pending.files,
            processingParams: {
              ...pending.processingParams,
              linkedSmDatasetId: smId,
            },
            onInitiated: (scId) => {
              void createCombinedProject(base, scId, smId);
            },
          });
        },
      });
      router.push("/explore");

      return;
    }

    // Single dataset (single-cell only, or a standalone single-molecule upload).
    startUpload({
      kind: pending.kind,
      title: base,
      files: pending.files,
      processingParams: smMapping
        ? { ...pending.processingParams, columnMapping: smMapping }
        : pending.processingParams,
    });
    router.push("/explore");
  };

  const handleFiles = async (files: File[]) => {
    if (files.length === 0) return;

    if (serverUpload) {
      await handleServerUpload(files);

      return;
    }

    try {
      resetHooks();
      setLoading(true);
      setError(null);
      setProgress(0);
      setProgressMessage("Starting...");

      const onProgress = (prog: number, msg: string) => {
        console.log(`[FileUpload] Progress: ${prog}% - ${msg}`);
        setProgress(prog);
        setProgressMessage(msg);
      };

      let dataset: StandardizedDataset | SingleMoleculeDataset;

      if (singleMolecule) {
        // Single molecule mode
        if (type === "chunked") {
          // Chunked folder upload - verify it's a valid chunked single molecule dataset
          if (!isChunkedSingleMoleculeFolder(files)) {
            throw new Error(
              "Invalid chunked single molecule dataset folder. Make sure it contains manifest.json.gz and genes/*.bin.gz files",
            );
          }

          toast.info(
            `Loading pre-chunked single molecule dataset (${files.length} files)...`,
          );
          console.log(
            "=== Loading pre-chunked single molecule dataset (ready for upload) ===",
          );
          console.log("Files:", files.length);

          onProgress(10, "Creating dataset from chunked files...");

          // Load dataset using ProcessedSingleMoleculeAdapter
          dataset = await SingleMoleculeDataset.fromLocalChunked(
            files,
            onProgress,
          );

          // Mark dataset as pre-chunked
          (dataset as any).isPreChunked = true;

          console.log("=== Pre-chunked Single Molecule Dataset loaded ===");
          console.log("Dataset ID:", dataset.id);
          console.log("Molecule count:", dataset.getMoleculeCount());
          console.log("Gene count:", dataset.genes.length);
          console.log("Spatial dimensions:", dataset.dimensions);
          console.log("Ready for upload!");
        } else {
          // Single parquet/csv file upload
          const file = files[0];
          const fileExtension = file.name.split(".").pop()?.toLowerCase();

          toast.info(`Processing ${file.name}...`);
          console.log(
            "=== Starting Single Molecule file processing (in worker) ===",
          );
          console.log("File:", file.name, "Size:", file.size, "bytes");
          console.log("File extension:", fileExtension);

          // For the unified "file" upload, sniff column names to pick the
          // right schema (xenium / merscope / custom). For the legacy
          // type="xenium" / type="merscope" cards, keep the explicit choice.
          let parquetDatasetType: MoleculeDatasetType =
            UPLOAD_TYPE_TO_PARQUET_TYPE[type];
          let mappingOverride: MoleculeColumnMapping | undefined;

          if (type === "file") {
            // Show the "confirm your columns" step: preview the file and let the
            // user check/remap gene / x / y / z / cell-id before processing.
            onProgress(2, "Inspecting file schema...");
            const mapping = await awaitColumnMapping(file);

            if (!mapping) {
              // User cancelled — abort cleanly (no error toast).
              setLoading(false);
              setProgress(0);
              setProgressMessage("");

              return;
            }
            mappingOverride = mapping;
            parquetDatasetType = "custom";
            console.log("[FileUpload] Confirmed column mapping:", mapping);
          } else {
            console.log("Dataset type:", type, "→", parquetDatasetType);
          }

          // Get singleton worker instance
          const workerApi = await getSingleMoleculeWorker();

          // Wrap progress callback with Comlink.proxy for cross-thread communication
          const proxiedProgress = Comlink.proxy(onProgress);

          // Parse file in web worker
          let serializedData;

          if (fileExtension === "parquet") {
            serializedData = await workerApi.parseParquet(
              file,
              parquetDatasetType,
              proxiedProgress,
              mappingOverride,
            );
          } else if (fileExtension === "csv") {
            serializedData = await workerApi.parseCSV(
              file,
              parquetDatasetType,
              proxiedProgress,
              mappingOverride,
            );
          } else {
            throw new Error(
              `Unsupported file type: .${fileExtension}. Only .parquet and .csv files are supported.`,
            );
          }

          // Reconstruct dataset from serialized data
          dataset = SingleMoleculeDataset.fromSerializedData(serializedData);

          console.log(
            "=== Single Molecule Dataset created successfully (from worker) ===",
          );
          console.log("Dataset ID:", dataset.id);
          console.log("Dataset name:", dataset.name);
          console.log("Dataset type:", dataset.type);
          console.log("Molecule count:", dataset.getMoleculeCount());
          console.log("Gene count:", dataset.genes.length);
          console.log("Spatial dimensions:", dataset.dimensions);
          console.log("Scaling factor:", dataset.scalingFactor);
          console.log("Summary:", dataset.getSummary());
        }
      } else {
        // Single cell mode.
        //
        // The unified "Folder" card uses type="folder" — auto-detect which of
        // zarr / chunked / xenium / merscope / h5ad-inside the dropped folder
        // is, then dispatch to the existing per-format branch below.
        let effectiveType: UploadType = type;

        if (type === "folder") {
          const detected = detectFolderFormat(files);

          if (!detected) {
            throw new Error(
              "Could not recognize this folder. Supported formats: " +
                "zarr (with .zgroup), pre-chunked (manifest.json + expr/), " +
                "Xenium (cells.csv), MERSCOPE (cell_metadata.csv), " +
                "or any folder containing a .h5ad file.",
            );
          }
          effectiveType = detected;
          console.log(`[FileUpload] Auto-detected folder format: ${detected}`);
          toast.info(`Detected ${detected.toUpperCase()} folder`);
        }

        if (effectiveType === "chunked") {
          // Chunked folder upload - verify it's a valid chunked dataset
          if (!isChunkedDatasetFolder(files)) {
            throw new Error(
              "Invalid chunked dataset folder. Make sure it contains manifest.json, expr/index.json, expr/chunk_*.bin.gz, and coords/spatial.bin.gz",
            );
          }
          toast.info(`Loading pre-chunked dataset (${files.length} files)...`);
          console.log("=== Loading pre-chunked dataset (ready for upload) ===");
          console.log("Files:", files.length);

          // Read manifest to get metadata
          onProgress(10, "Reading manifest...");
          const manifestFile = files.find((f) => {
            const path = f.webkitRelativePath || f.name;
            const parts = path.split("/");
            const relativePath =
              parts.length > 1 ? parts.slice(1).join("/") : path;

            return relativePath === "manifest.json";
          });

          if (!manifestFile) {
            throw new Error("manifest.json not found in chunked folder");
          }

          const manifestText = await manifestFile.text();
          const manifest = JSON.parse(manifestText);

          console.log("Manifest loaded:", manifest);

          // Create file map for later upload
          const fileMap = new Map<string, File>();

          for (const file of files) {
            const relativePath = file.webkitRelativePath;
            const parts = relativePath.split("/");
            const fileKey = parts.slice(1).join("/"); // Remove root folder

            fileMap.set(fileKey, file);
          }

          onProgress(30, "Creating dataset from manifest...");

          // Create a minimal StandardizedDataset from manifest
          // We don't actually load the chunked data, just create metadata
          dataset = await StandardizedDataset.fromLocalChunked(
            files,
            onProgress,
          );

          // Mark dataset as pre-chunked and attach file map
          (dataset as any).isPreChunked = true;
          (dataset as any).chunkedFiles = fileMap;
          (dataset as any).manifestData = manifest;

          onProgress(100, "Pre-chunked dataset ready for upload!");
          console.log(
            "Pre-chunked dataset ready. File map size:",
            fileMap.size,
          );
        } else if (effectiveType === "h5ad") {
          // h5ad can come from the .h5ad file card OR an auto-detected
          // folder containing a .h5ad (e.g., a MERSCOPE export with one).
          const h5adFile = findH5adFile(files) ?? files[0];

          toast.info(`Processing ${h5adFile.name}...`);
          console.log("=== Starting H5AD file processing ===");
          console.log("File:", h5adFile.name, "Size:", h5adFile.size, "bytes");
          dataset = await StandardizedDataset.fromH5ad(h5adFile, onProgress);
        } else if (effectiveType === "h5ad-zarr") {
          if (!isZarrFolder(files)) {
            throw new Error(
              "Dropped folder is not a zarr store (expected a .zgroup, .zarray, or zarr.json marker at the root).",
            );
          }
          toast.info(`Processing zarr folder (${files.length} files)...`);
          console.log("=== Starting h5ad-zarr folder processing ===");
          console.log("Files:", files.length);
          dataset = await StandardizedDataset.fromH5adZarr(files, onProgress);
        } else if (effectiveType === "xenium") {
          toast.info(`Processing Xenium folder (${files.length} files)...`);
          console.log("=== Starting Xenium folder processing ===");
          console.log("Files:", files.length);
          dataset = await StandardizedDataset.fromXenium(files, onProgress);
        } else {
          toast.info(`Processing MERSCOPE folder (${files.length} files)...`);
          console.log("=== Starting MERSCOPE folder processing ===");
          console.log("Files:", files.length);
          dataset = await StandardizedDataset.fromMerscope(files, onProgress);
        }

        console.log("=== Dataset created successfully ===");
        console.log("Dataset ID:", dataset.id);
        console.log("Dataset name:", dataset.name);
        console.log("Dataset type:", dataset.type);
        console.log(
          "Point count:",
          (dataset as StandardizedDataset).getPointCount(),
        );
        console.log("Gene count:", dataset.genes.length);
        console.log(
          "Spatial dimensions:",
          (dataset as StandardizedDataset).spatial.dimensions,
        );
        console.log(
          "Available embeddings:",
          Object.keys((dataset as StandardizedDataset).embeddings),
        );
        console.log("Cluster info:", (dataset as StandardizedDataset).clusters);
        console.log("Summary:", dataset.getSummary());
      }

      // Add to appropriate Zustand store
      setProgress(100);
      setProgressMessage("Complete!");

      if (singleMolecule) {
        const sm = dataset as SingleMoleculeDataset;

        smStore.addDataset(sm);
        markDatasetLoaded({
          geneCount: sm.uniqueGenes.length,
          dimensions: sm.dimensions,
          dataType: "single_molecule",
        });
      } else {
        const sc = dataset as StandardizedDataset;

        cellStore.addDataset(sc);
        markDatasetLoaded({
          pointCount: sc.getPointCount(),
          geneCount: sc.genes.length,
          dimensions: sc.spatial.dimensions,
          dataType: "single_cell",
        });
      }

      toast.success(`Dataset loaded successfully!`);
      setLoading(false);

      // Navigate to from-local viewer (not /{id} which is for already-uploaded datasets)
      router.push(
        singleMolecule ? `/sm-viewer/from-local` : `/viewer/from-local`,
      );
    } catch (error) {
      console.error(`Error processing ${type} data:`, error);
      const errorMessage =
        error instanceof Error ? error.message : "Unknown error";

      setError(errorMessage);
      toast.error(`Failed to process data: ${errorMessage}`);
      setLoading(false);
      setProgress(0);
      setProgressMessage("");
      // Some datasets that fail in the browser (memory limits, int64 matrices,
      // unusual encodings) process fine on the server — offer that path with
      // the file the user already dropped.
      setServerFallbackFiles(files);
    }
  };

  const handleClick = () => {
    document.getElementById(inputId)?.click();
  };

  const isFolder =
    (type === "xenium" ||
      type === "merscope" ||
      type === "chunked" ||
      type === "h5ad-zarr" ||
      type === "folder") &&
    (singleMolecule ? type === "chunked" : true);

  // Get isLoading from appropriate store
  const cellIsLoading = useDatasetStore((state) => state.isLoading);
  const smIsLoading = useSingleMoleculeStore((state) => state.isLoading);
  const isLoading = singleMolecule ? smIsLoading : cellIsLoading;

  // Determine accepted file types
  const getAcceptedFileTypes = () => {
    if (singleMolecule) {
      return ".parquet,.csv";
    }
    if (type === "file") {
      return ".parquet,.csv,.tsv,.txt";
    }

    return type === "h5ad" ? ".h5ad" : ".csv,.tsv,.txt";
  };

  return (
    <div className="w-full">
      <div
        aria-label={`${title} — ${description}. Click or drop files to upload.`}
        className={`
          relative border-2 border-dashed rounded-lg p-3 text-center cursor-pointer
          transition-all duration-200 ease-in-out flex items-center justify-center
          min-h-[6rem]
          ${
            isDragging
              ? "border-primary bg-primary/10 scale-[1.02]"
              : "border-default-300 hover:border-primary/50 hover:bg-default-100/50"
          }
          ${isLoading ? "pointer-events-none opacity-60" : ""}
        `}
        role="button"
        tabIndex={isLoading ? -1 : 0}
        onClick={handleClick}
        onDragLeave={handleDragLeave}
        onDragOver={handleDragOver}
        onDrop={handleDrop}
        onKeyDown={(e) => {
          if (e.key === "Enter" || e.key === " ") {
            e.preventDefault();
            handleClick();
          }
        }}
      >
        <input
          accept={getAcceptedFileTypes()}
          className="hidden"
          id={inputId}
          multiple={isFolder}
          type="file"
          onChange={handleFileInput}
          {...(isFolder ? { webkitdirectory: "", directory: "" } : {})}
        />

        <div className="flex flex-col items-center gap-1.5 w-full">
          <p className="text-lg font-semibold text-foreground">{title}</p>
          <p className="text-sm text-default-500">{description}</p>

          {icons ? (
            <div className="flex items-center justify-center gap-2 text-default-500">
              {icons}
            </div>
          ) : null}

          {isLoading && progress > 0 && (
            <div className="w-full mt-4 px-4">
              <div className="w-full bg-default-200 rounded-full h-2 overflow-hidden">
                <div
                  className="bg-primary h-full transition-all duration-300 ease-out"
                  data-progress={progress}
                  data-testid="upload-progress"
                  style={{ width: `${progress}%` }}
                />
              </div>
              <p className="text-sm text-default-500 mt-2">{progressMessage}</p>
            </div>
          )}
        </div>
      </div>

      {serverFallbackFiles &&
        createPortal(
          <div
            className="fixed inset-0 z-[var(--z-modal-top)] flex items-center justify-center bg-black/70 p-6 backdrop-blur-md"
            role="presentation"
            onClick={(e) => {
              if (e.target === e.currentTarget) setServerFallbackFiles(null);
            }}
          >
            <div className="glass-panel w-full max-w-md rounded-2xl p-7" role="dialog">
              <h2 className="text-xl font-semibold text-foreground">
                Process on our server instead?
              </h2>
              <p className="mt-2 text-sm text-default-500">
                This dataset could not be processed in your browser. Our server
                handles larger datasets and formats the browser can\u2019t \u2014 same
                zero-setup flow, and you\u2019ll get a shareable link when it\u2019s ready.
              </p>
              <div className="mt-6 flex justify-end gap-2">
                <Button
                  variant="light"
                  onPress={() => setServerFallbackFiles(null)}
                >
                  Cancel
                </Button>
                <Button
                  color="primary"
                  onPress={() => {
                    const files = serverFallbackFiles;

                    setServerFallbackFiles(null);
                    if (files) void handleServerUpload(files);
                  }}
                >
                  Process on server
                </Button>
              </div>
            </div>
          </div>,
          document.body,
        )}

      {columnPreviewLoading && (
        <div className="fixed inset-0 z-[var(--z-modal-top)] flex flex-col items-center justify-center gap-3 bg-black/70 backdrop-blur-md">
          <Spinner size="lg" />
          <p className="text-sm text-default-300">Reading column headers…</p>
        </div>
      )}

      {headerModal && (
        <SingleMoleculeHeaderModal
          autoMapping={headerModal.autoMapping}
          columns={headerModal.columns}
          fileName={headerModal.fileName}
          isOpen={true}
          rows={headerModal.rows}
          onCancel={() => resolveMapping(null)}
          onConfirm={(m) => resolveMapping(m)}
        />
      )}

      {serverUpload && (
        <AuthModal
          isOpen={authModalOpen}
          onAuthenticated={() => {
            setAuthModalOpen(false);
            const files = pendingFilesRef.current;

            pendingFilesRef.current = null;
            if (files) doServerUpload(files);
          }}
          onClose={() => {
            setAuthModalOpen(false);
            pendingFilesRef.current = null;
          }}
        />
      )}

      {/* Name-your-dataset step for the server upload — pre-filled from the
          file/folder name, editable so multiple same-named files don't collide. */}
      <Modal
        disableAnimation
        isOpen={pendingUpload !== null}
        size="md"
        onClose={() => setPendingUpload(null)}
      >
        <ModalContent>
          <ModalHeader className="flex flex-col gap-1">
            <span>Name your dataset</span>
            <span className="text-xs font-normal text-default-400">
              {pendingUpload?.files.length ?? 0}{" "}
              {(pendingUpload?.files.length ?? 0) === 1 ? "file" : "files"} ·{" "}
              {pendingUpload?.kind === "single_molecule"
                ? "single molecule"
                : "single cell"}
            </span>
          </ModalHeader>
          <ModalBody className="gap-4">
            <Input
              description={
                pendingUpload?.sm && includeSm
                  ? "Used for the project. Datasets become “… (cells)” and “… (molecules)”."
                  : "You can reuse a name — this is just how it's labeled."
              }
              label={
                pendingUpload?.sm && includeSm ? "Project name" : "Dataset name"
              }
              value={titleDraft}
              onKeyDown={(e) => {
                if (e.key === "Enter" && titleDraft.trim()) {
                  startPreparedUpload();
                }
              }}
              onValueChange={setTitleDraft}
            />

            {pendingUpload?.sm && (
              <div className="rounded-xl border border-default-200 bg-default-100/50 p-3">
                <Checkbox
                  isSelected={includeSm}
                  size="sm"
                  onValueChange={setIncludeSm}
                >
                  <span className="text-sm">
                    Also upload the single-molecule data
                  </span>
                </Checkbox>
                <p className="mt-1 ml-7 text-xs text-default-500">
                  We detected single-molecule data in{" "}
                  <span className="font-medium text-default-600">
                    {pendingUpload.sm.fileName}
                  </span>
                  . It&apos;ll upload as a separate dataset, and we&apos;ll
                  group both into one project and overlay the molecules on the
                  cells.
                </p>
              </div>
            )}
          </ModalBody>
          <ModalFooter>
            <Button variant="light" onPress={() => setPendingUpload(null)}>
              Cancel
            </Button>
            <Button
              color="primary"
              isDisabled={!titleDraft.trim()}
              onPress={startPreparedUpload}
            >
              {pendingUpload?.sm && includeSm ? "Upload both" : "Start upload"}
            </Button>
          </ModalFooter>
        </ModalContent>
      </Modal>
    </div>
  );
}

/**
 * Recursively walk a list of FileSystemEntry trees (from dataTransfer.items
 * via webkitGetAsEntry) and return all leaf File objects.
 *
 * Each yielded File gets its `webkitRelativePath` patched to the dropped-tree
 * path so downstream detection (`isXeniumFolder`, etc.) and file-map building
 * (`fromH5adZarr`, `fromLocalChunked`) behave the same as the click-to-pick
 * `<input webkitdirectory>` path.
 */
async function collectFilesFromEntries(entries: any[]): Promise<File[]> {
  const out: File[] = [];

  await Promise.all(entries.map((entry) => walkEntry(entry, "", out)));

  return out;
}

async function walkEntry(
  entry: any,
  parentPath: string,
  out: File[],
): Promise<void> {
  if (!entry) return;

  const namePath = parentPath ? `${parentPath}/${entry.name}` : entry.name;

  if (entry.isFile) {
    const file: File = await new Promise((resolve, reject) => {
      entry.file(resolve, reject);
    });

    // webkitRelativePath is a getter on File.prototype; defining an own
    // property on the instance shadows it. This is the standard trick for
    // drag-and-drop folder uploads — supported in Chromium, Firefox, Safari.
    try {
      Object.defineProperty(file, "webkitRelativePath", {
        value: namePath,
        configurable: true,
        writable: false,
      });
    } catch {
      // Some engines lock down the prototype getter; fall back silently.
      // The detection helpers will still see file.name; folder routing may
      // misclassify but the user gets a friendly error.
    }
    out.push(file);

    return;
  }

  if (entry.isDirectory) {
    const reader = entry.createReader();

    // readEntries returns up to ~100 entries per call; loop until empty.
    while (true) {
      const batch: any[] = await new Promise((resolve, reject) => {
        reader.readEntries(resolve, reject);
      });

      if (batch.length === 0) break;
      await Promise.all(batch.map((child) => walkEntry(child, namePath, out)));
    }
  }
}
