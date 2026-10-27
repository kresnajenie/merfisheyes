import { columnPalette, compareLabels } from "@/lib/utils/color-palette";

/**
 * One obs column built from a row table (cells.csv / cell_metadata.csv), in
 * the adapters' ClusterData shape: indexed representation (sorted unique
 * labels + per-cell codes) and a default palette dealt in compareLabels
 * order — the order the SpatialData viewer deals the saved store's colours
 * in, so a column keeps its colours when it is saved.
 */
export interface TableObsColumn {
  column: string;
  values: string[];
  valueIndices: Uint16Array | Uint32Array;
  palette: Record<string, string>;
  uniqueValues: string[];
}

const ID_KEYS = ["cell_id", "id", "barcode", "barcodes", "cell", "cellid", "EntityID"];
const Z_KEYS = ["cell_centroid_z", "z_centroid", "centroid_z", "center_z", "centerz", "z"];

/**
 * Every metadata column of a row table worth exporting — mirrors the server
 * pipeline's rule (scripts/process_spatial_data.py, Xenium + MERSCOPE
 * loaders): skip the coordinate columns, id columns, gene-name columns and
 * columns that are entirely empty. Missing values become the label "".
 */
export function buildObsColumnsFromRows(
  rows: Record<string, any>[],
  opts: {
    coordinateKeys: (string | null | undefined)[];
    extraExcludeKeys?: (string | null | undefined)[];
    geneNames?: Iterable<string>;
  },
): TableObsColumn[] {
  if (rows.length === 0) return [];

  const exclude = new Set<string>();

  for (const k of opts.coordinateKeys) if (k) exclude.add(k);
  for (const k of opts.extraExcludeKeys ?? []) if (k) exclude.add(k);
  const lowerKeys = new Map(Object.keys(rows[0]).map((k) => [k.toLowerCase(), k]));

  for (const cand of [...ID_KEYS, ...Z_KEYS]) {
    const actual = lowerKeys.get(cand.toLowerCase());

    if (actual) exclude.add(actual);
  }
  const genes = new Set(opts.geneNames ?? []);
  const columns: TableObsColumn[] = [];

  for (const key of Object.keys(rows[0])) {
    if (exclude.has(key) || genes.has(key)) continue;
    // A column with no name is the file's row index (MERSCOPE's
    // cell_metadata.csv starts with one), never a property of the cells;
    // nor is an adapter's own "_" field (a MERSCOPE row's _umap_x / _umap_y)
    if (key.trim() === "" || /^Unnamed: \d+$/.test(key) || key.startsWith("_")) {
      continue;
    }

    const vals = new Array<string>(rows.length);
    let anyValue = false;

    for (let i = 0; i < rows.length; i++) {
      const raw = rows[i]?.[key];
      const str = raw == null ? "" : String(raw).trim() === "" ? "" : String(raw);

      if (str !== "") anyValue = true;
      vals[i] = str;
    }
    if (!anyValue) continue; // entirely empty column

    const uniq = Array.from(new Set(vals)).sort(compareLabels);
    const index = new Map(uniq.map((u, i) => [u, i]));
    const IndexArray = uniq.length <= 65535 ? Uint16Array : Uint32Array;
    const valueIndices = new IndexArray(vals.length);

    for (let i = 0; i < vals.length; i++) valueIndices[i] = index.get(vals[i])!;
    columns.push({
      column: key,
      values: [],
      valueIndices,
      palette: columnPalette(uniq),
      uniqueValues: uniq,
    });
  }

  return columns;
}
