import { parquetMetadataAsync } from "hyparquet";

import {
  MOLECULE_COLUMN_MAPPINGS,
  MoleculeDatasetType,
} from "../config/moleculeColumnMappings";

/**
 * Read a parquet file's column names by parsing only the footer/header
 * metadata — no row data. Works on Files of any size.
 */
async function readParquetColumns(file: File): Promise<string[]> {
  const asyncBuffer = {
    byteLength: file.size,
    slice: async (start: number, end?: number): Promise<ArrayBuffer> => {
      return await file.slice(start, end).arrayBuffer();
    },
  };
  const meta = await parquetMetadataAsync(asyncBuffer);

  // meta.schema is a depth-first list of schema nodes with name + path.
  // Leaf names are the column names. We collect every node name minus the
  // root (typically "schema") to be permissive across writers.
  const names = new Set<string>();

  for (const node of meta.schema) {
    if (!node || !node.name) continue;
    if (node.name === "schema" || node.name === "root") continue;
    names.add(node.name);
  }

  return Array.from(names);
}

/** Read a CSV's header row by slicing the first chunk of the file. */
async function readCsvHeader(file: File): Promise<string[]> {
  const blob = file.slice(0, 16 * 1024);
  const text = await blob.text();
  const firstLine = text.split(/\r?\n/)[0] ?? "";

  if (!firstLine) return [];

  // Tab- or comma-separated; we don't try to parse quoted commas because
  // headers rarely have them.
  const sep = firstLine.includes("\t") && !firstLine.includes(",")
    ? "\t"
    : ",";

  return firstLine.split(sep).map((s) => s.trim().replace(/^"|"$/g, ""));
}

/**
 * Pick the right MoleculeDatasetType from a list of column names.
 *
 * Heuristic: each known schema declares a (gene, x) column pair. We pick the
 * first schema whose pair is fully present. Falls back to "custom" so the
 * existing default mapping kicks in (and the user can adjust later).
 */
function pickSchema(columnNames: string[]): MoleculeDatasetType {
  const cols = new Set(columnNames);

  for (const schema of ["xenium", "merscope"] as const) {
    const map = MOLECULE_COLUMN_MAPPINGS[schema];

    if (cols.has(map.gene) && cols.has(map.x)) {
      return schema;
    }
  }

  return "custom";
}

/**
 * Detect a single-molecule file's schema (xenium / merscope / custom)
 * by inspecting its column names.
 */
export async function detectMoleculeFileType(
  file: File,
): Promise<{ type: MoleculeDatasetType; columns: string[] }> {
  const ext = file.name.toLowerCase().split(".").pop();
  let columns: string[] = [];

  if (ext === "parquet") {
    columns = await readParquetColumns(file);
  } else if (ext === "csv" || ext === "tsv" || ext === "txt") {
    columns = await readCsvHeader(file);
  } else {
    throw new Error(
      `Unsupported file type: .${ext}. Expected .parquet, .csv, .tsv, or .txt`,
    );
  }

  return { type: pickSchema(columns), columns };
}
