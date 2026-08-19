/**
 * Decide which files in a picked folder are worth uploading, and whether the
 * folder holds single-cell data, single-molecule data, or both.
 *
 * Server-side ingestion uploads raw bytes, so an unfiltered folder means
 * sending gigabytes the processor will never open — a Xenium export is mostly
 * images and QC artefacts. This narrows the upload to the files the Python
 * actually reads.
 *
 * IMPORTANT — these rules mirror the processing scripts and must be kept in
 * step with them:
 *   - `detect_input_format()` in scripts/process_spatial_data.py (Xenium/MERSCOPE)
 *   - `load_xenium_data()` in the same file (which expression files it tries)
 *   - `match_transcript_file()` in scripts/process_single_molecule.py
 * The tests in tests/unit/classify-folder.test.ts pin the behaviour against real
 * folder listings.
 */

export type SingleCellFormat = "h5ad" | "xenium" | "merscope";

export interface ClassifiedFile {
  /** Key relative to the dataset root, e.g. "cell_metadata.csv". */
  key: string;
  file: File;
}

export interface FolderClassification {
  singleCell?: { format: SingleCellFormat; files: ClassifiedFile[] };
  singleMolecule?: { files: ClassifiedFile[] };
  /** Everything not selected, so the UI can say what it is leaving behind. */
  ignored: ClassifiedFile[];
}

/**
 * Strip the picked directory's own name from a file's relative path.
 *
 * The worker unpacks `raw/{id}/` and hands that directory straight to the
 * processor, whose format detection only looks at files *directly* inside it.
 * Keeping `webkitRelativePath` verbatim would bury everything one level down
 * under the folder's name and detection would find nothing. Deeper structure
 * (e.g. `cell_feature_matrix/features.tsv.gz`) is preserved, because the Xenium
 * loader looks for exactly that.
 */
export function datasetRelativeKey(file: File): string {
  const path = file.webkitRelativePath || file.name;
  const parts = path.split("/");

  return parts.length > 1 ? parts.slice(1).join("/") : path;
}

/** macOS and archive-tool droppings; never part of a dataset. */
function isJunk(key: string): boolean {
  return key
    .split("/")
    .some(
      (segment) =>
        segment.startsWith(".") ||
        segment.startsWith("._") ||
        segment === "__MACOSX",
    );
}

function baseName(key: string): string {
  return key.slice(key.lastIndexOf("/") + 1);
}

/** Top level of the dataset, i.e. no "/" left after the root was stripped. */
function isTopLevel(key: string): boolean {
  return !key.includes("/");
}

/**
 * Mirrors normalize_filename() in process_single_molecule.py: drop the
 * extension, lowercase, and turn separators into spaces.
 */
function normalizeStem(name: string): string {
  const stem = name.replace(/\.[^.]+$/, "");

  return stem.toLowerCase().replace(/[-_.]/g, " ");
}

/**
 * Mirrors match_transcript_file(). Returns 0 for a non-match.
 *
 * `.parquet` is accepted here even though the Python pattern lists only `.csv`:
 * that restriction applies to its *directory scan*, while read_parquet_file()
 * handles parquet fine when the path is passed explicitly — which is how the
 * worker invokes it.
 */
export function scoreTranscriptFile(name: string): number {
  const lower = name.toLowerCase();

  if (!lower.endsWith(".csv") && !lower.endsWith(".parquet")) return 0;

  const normalized = normalizeStem(name);

  for (const excluded of ["metadata", "cell by gene"]) {
    if (normalized.includes(excluded)) return 0;
  }

  const hasKeyword = (kw: string) =>
    normalized.includes(kw) || normalized.includes(`${kw}s`);

  if (hasKeyword("detected") && hasKeyword("transcript")) return 100;
  if (hasKeyword("transcript")) return 80;

  return 0;
}

function isMerscopeMetadata(name: string): boolean {
  const lower = name.toLowerCase();

  return (
    lower.endsWith(".csv") &&
    lower.includes("metadata") &&
    !lower.includes("gene")
  );
}

function isMerscopeCellByGene(name: string): boolean {
  const lower = name.toLowerCase();

  return (
    lower.endsWith(".csv") && lower.includes("cell") && lower.includes("gene")
  );
}

/** Files load_xenium_data() will look for, beyond cells.csv itself. */
const XENIUM_EXPRESSION = [
  "cell_by_gene.csv",
  "cell_feature_matrix.h5",
  "features.tsv",
  "features.tsv.gz",
  "barcodes.tsv",
  "barcodes.tsv.gz",
  "matrix.mtx",
  "matrix.mtx.gz",
];

function isXeniumCells(key: string): boolean {
  const name = baseName(key).toLowerCase();

  return isTopLevel(key) && (name === "cells.csv" || name === "cells.csv.gz");
}

function isXeniumExpression(key: string): boolean {
  const name = baseName(key).toLowerCase();

  if (!XENIUM_EXPRESSION.includes(name)) return false;

  // Either directly in the dataset root, or in the cell_feature_matrix/ folder
  // the Xenium loader also checks.
  return (
    isTopLevel(key) || key.toLowerCase().startsWith("cell_feature_matrix/")
  );
}

/**
 * Fallback for an EXPLICIT single-molecule upload whose file matches no
 * transcripts-style name (e.g. "spots_set4.parquet"). The worker's
 * find_single_molecule_input() accepts any lone .csv/.parquet under the raw
 * dir regardless of its name — the keyword scoring only exists to pick one
 * file out of many. So when the user already said "this is single-molecule
 * data" and there is exactly one candidate, take it. Multiple unnamed
 * candidates stay ambiguous (null) so we never guess which file to upload.
 *
 * Only for the single-molecule upload mode — never for the combined-upload
 * detection inside single-cell folders, where a loose match would misread
 * arbitrary CSVs as transcripts.
 */
export function fallbackSingleMoleculeFile(
  files: File[],
): { files: ClassifiedFile[] } | null {
  const candidates = files
    .map((file) => ({ key: datasetRelativeKey(file), file }))
    .filter((e) => !isJunk(e.key))
    .filter((e) => /\.(csv|parquet)$/i.test(baseName(e.key)));

  return candidates.length === 1 ? { files: candidates } : null;
}

export function classifyFolder(files: File[]): FolderClassification {
  const entries: ClassifiedFile[] = files.map((file) => ({
    key: datasetRelativeKey(file),
    file,
  }));

  const usable = entries.filter((e) => !isJunk(e.key));
  const claimed = new Set<ClassifiedFile>();

  const claim = (matches: ClassifiedFile[]) => {
    matches.forEach((m) => claimed.add(m));

    return matches;
  };

  let singleCell: FolderClassification["singleCell"];

  // A folder of .h5ad files is a single-cell dataset per file; the caller picks
  // one. Checked first because an .h5ad is unambiguous.
  const h5ads = usable.filter(
    (e) => isTopLevel(e.key) && e.key.toLowerCase().endsWith(".h5ad"),
  );

  const xeniumCells = usable.filter((e) => isXeniumCells(e.key));
  const merscopeFiles = usable.filter(
    (e) =>
      isTopLevel(e.key) &&
      (isMerscopeMetadata(baseName(e.key)) ||
        isMerscopeCellByGene(baseName(e.key))),
  );

  if (h5ads.length > 0) {
    singleCell = { format: "h5ad", files: claim(h5ads) };
  } else if (xeniumCells.length > 0) {
    // Xenium wins over MERSCOPE when both could match, matching the order
    // detect_input_format() checks them in.
    const expression = usable.filter((e) => isXeniumExpression(e.key));

    singleCell = {
      format: "xenium",
      files: claim([...xeniumCells, ...expression]),
    };
  } else if (merscopeFiles.length > 0) {
    singleCell = { format: "merscope", files: claim(merscopeFiles) };
  }

  // Single molecule: the single best-scoring transcripts file, matching the
  // Python's "best match wins" directory scan.
  let singleMolecule: FolderClassification["singleMolecule"];
  let best: ClassifiedFile | undefined;
  let bestScore = 0;

  for (const entry of usable) {
    if (claimed.has(entry)) continue;
    const score = scoreTranscriptFile(baseName(entry.key));

    if (score > bestScore) {
      best = entry;
      bestScore = score;
    }
  }

  if (best) {
    singleMolecule = { files: claim([best]) };
  }

  return {
    singleCell,
    singleMolecule,
    ignored: entries.filter((e) => !claimed.has(e)),
  };
}
