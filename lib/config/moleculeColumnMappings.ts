export interface MoleculeColumnMapping {
  gene: string;
  x: string;
  y: string;
  z: string;
  cellId: string | null;
}

export type MoleculeDatasetType = "xenium" | "merscope" | "custom";

/**
 * Column name mappings for different single molecule dataset types
 */
export const MOLECULE_COLUMN_MAPPINGS: Record<
  MoleculeDatasetType,
  MoleculeColumnMapping
> = {
  xenium: {
    gene: "feature_name",
    x: "x_location",
    y: "y_location",
    z: "z_location",
    cellId: null,
  },
  merscope: {
    gene: "gene",
    x: "global_x",
    y: "global_y",
    z: "global_z",
    cellId: "cell_id",
  },
  custom: {
    gene: "feature_name",
    x: "x_location",
    y: "y_location",
    z: "z_location",
    cellId: null,
  },
};

// Cell-id values that mean "this molecule wasn't assigned to a cell". Each
// platform encodes it differently:
//   MERSCOPE (Vizgen)  → integer  -1
//   Xenium   (10x)     → string   "UNASSIGNED"
//   CosMx    (NanoString) → integer 0  (assigned cell_ID is >= 1, per FOV)
// So any numeric <= 0 is unassigned (covers -1 and 0), plus NaN/blank and the
// common string tokens. A non-token, non-numeric string (e.g. a Xenium cell id
// like "aaabbbcc-1") is an ASSIGNED cell.
const UNASSIGNED_TOKENS = new Set([
  "unassigned",
  "none",
  "nan",
  "na",
  "null",
  "",
]);

/** Whether a cell-id value denotes an unassigned molecule. */
export function isUnassignedCellId(value: unknown): boolean {
  if (value === null || value === undefined) return true;
  if (typeof value === "number") return Number.isNaN(value) || value <= 0;

  const s = String(value).trim().toLowerCase();

  if (UNASSIGNED_TOKENS.has(s)) return true;
  // Numeric-looking string (e.g. CosMx "0", MERSCOPE "-1"): unassigned if <= 0.
  const n = Number(s);

  return !Number.isNaN(n) && n <= 0;
}
