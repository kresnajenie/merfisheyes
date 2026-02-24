export interface MoleculeColumnMapping {
  gene: string;
  x: string;
  y: string;
  z: string;
  cellId?: string;
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
    cellId: "cell_id",
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
  },
};
