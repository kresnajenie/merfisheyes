import type { MoleculeDatasetType } from "../config/moleculeColumnMappings";

import { pickSchema } from "./molecule-file-sniffer";

export interface MoleculePreview {
  columns: string[];
  /** First few rows, as { column: value } objects, for the confirm UI. */
  rows: Record<string, unknown>[];
  /** Auto-detected schema (xenium / merscope / custom) from the columns. */
  autoType: MoleculeDatasetType;
}

/**
 * Read a single-molecule file's columns AND a handful of sample rows, cheaply
 * and on the main thread — enough for the "confirm your columns" step. Parquet
 * reads only the first `nRows`; CSV uses PapaParse's `preview`.
 */
export async function readMoleculePreview(
  file: File,
  nRows = 8,
): Promise<MoleculePreview> {
  const ext = file.name.toLowerCase().split(".").pop();

  if (ext === "parquet") {
    const { parquetReadObjects, parquetMetadataAsync } = await import(
      "hyparquet"
    );
    const { compressors } = await import("hyparquet-compressors");
    const asyncBuffer = {
      byteLength: file.size,
      slice: async (start: number, end?: number): Promise<ArrayBuffer> =>
        file.slice(start, end).arrayBuffer(),
    };

    const meta = await parquetMetadataAsync(asyncBuffer);
    const columns: string[] = [];

    for (const node of meta.schema) {
      if (!node || !node.name) continue;
      if (node.name === "schema" || node.name === "root") continue;
      columns.push(node.name);
    }

    const rows = (await parquetReadObjects({
      file: asyncBuffer,
      compressors,
      rowStart: 0,
      rowEnd: nRows,
    })) as Record<string, unknown>[];

    return { columns, rows, autoType: pickSchema(columns) };
  }

  if (ext === "csv" || ext === "tsv" || ext === "txt") {
    const Papa = (await import("papaparse")).default;

    return new Promise<MoleculePreview>((resolve, reject) => {
      // `preview` on a File input isn't in PapaParse's local-config types, so
      // cast like the streaming parse in SingleMoleculeDataset does.

      (Papa.parse as any)(file, {
        header: true,
        preview: nRows,
        skipEmptyLines: true,
        complete: (res: {
          data: Record<string, unknown>[];
          meta: { fields?: string[] };
        }) => {
          const columns = (res.meta.fields ?? []).map((c) => String(c).trim());

          resolve({
            columns,
            rows: res.data,
            autoType: pickSchema(columns),
          });
        },
        error: (err: unknown) => reject(err),
      });
    });
  }

  throw new Error(
    `Unsupported file type: .${ext}. Expected .parquet, .csv, .tsv, or .txt`,
  );
}
