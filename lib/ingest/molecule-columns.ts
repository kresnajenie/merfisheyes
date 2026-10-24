/**
 * Whether a table's columns are molecules: a gene column and x / y
 * coordinates, by name in any case. The same names, in the same order, as
 * `MOLECULE_COLUMN_NAMES` in worker/entrypoint.py, which reads the columns
 * of whatever is uploaded as molecules — keep the two equal.
 */

const NAMES = {
  gene: [
    "gene",
    "feature_name",
    "feature",
    "target",
    "gene_name",
    "target_name",
    "genes",
  ],
  x: ["x", "x_location", "global_x", "x_um", "x_global", "x_position"],
  y: ["y", "y_location", "global_y", "y_um", "y_global", "y_position"],
  z: ["z", "z_location", "global_z", "z_um", "z_global", "z_position"],
} as const;

export type MoleculeColumns = {
  gene: string;
  x: string;
  y: string;
  z?: string;
};

/** gene / x / y (and z) among `columns`; null without all three. */
export function matchMoleculeColumns(
  columns: string[],
): MoleculeColumns | null {
  const byName = new Map<string, string>();

  for (const column of columns) {
    const name = column.trim().toLowerCase();

    if (!byName.has(name)) byName.set(name, column);
  }
  const pick = (role: keyof typeof NAMES) =>
    NAMES[role].map((n) => byName.get(n)).find((c) => c !== undefined);
  const [gene, x, y, z] = [pick("gene"), pick("x"), pick("y"), pick("z")];

  return gene && x && y ? { gene, x, y, ...(z && { z }) } : null;
}
