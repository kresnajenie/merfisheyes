import type { MoleculeDatasetType } from "../config/moleculeColumnMappings";

import { pickSchema, readCsvHeader } from "./molecule-file-sniffer";

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

    // PapaParse streams a File in chunks. On very large files it can finish
    // without ever surfacing the header (meta.fields empty), which left the
    // confirm-columns modal blank with nothing to map. Treat anything that
    // does not yield columns — an error, a hang, or an empty field list — as
    // "no preview" and fall back to reading just the header bytes, which is
    // O(16 KB) regardless of file size.
    const PREVIEW_TIMEOUT_MS = 15_000;

    const viaPapa = await new Promise<{
      columns: string[];
      rows: Record<string, unknown>[];
    } | null>((resolve) => {
      let settled = false;
      const done = (
        v: { columns: string[]; rows: Record<string, unknown>[] } | null,
      ) => {
        if (settled) return;
        settled = true;
        resolve(v);
      };
      const timer = setTimeout(() => done(null), PREVIEW_TIMEOUT_MS);

      try {
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
            clearTimeout(timer);
            done({
              columns: (res.meta?.fields ?? []).map((c) => String(c).trim()),
              rows: res.data ?? [],
            });
          },
          error: () => {
            clearTimeout(timer);
            done(null);
          },
        });
      } catch {
        clearTimeout(timer);
        done(null);
      }
    });

    if (viaPapa && viaPapa.columns.length > 0) {
      return {
        columns: viaPapa.columns,
        rows: viaPapa.rows,
        autoType: pickSchema(viaPapa.columns),
      };
    }

    // Fallback: header bytes only. No sample rows, but the dropdowns are
    // populated and auto-detection still works.
    const columns = await readCsvHeader(file);

    return { columns, rows: [], autoType: pickSchema(columns) };
  }

  throw new Error(
    `Unsupported file type: .${ext}. Expected .parquet, .csv, .tsv, or .txt`,
  );
}
