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

// Cell-id values that mean "this molecule wasn't assigned to a cell". Different
// platforms encode it differently: MERSCOPE uses numeric -1, Xenium uses the
// STRING "UNASSIGNED", and various exports use 0 / "" / "None". The old code
// only caught numeric -1, so Xenium's unassigned transcripts were silently
// counted as assigned.
const UNASSIGNED_TOKENS = new Set([
  "unassigned",
  "none",
  "nan",
  "na",
  "null",
  "",
  "-1",
]);

/** Whether a cell-id value denotes an unassigned molecule. */
export function isUnassignedCellId(value: unknown): boolean {
  if (value === null || value === undefined) return true;
  if (typeof value === "number") return value === -1 || Number.isNaN(value);

  return UNASSIGNED_TOKENS.has(String(value).trim().toLowerCase());
}
