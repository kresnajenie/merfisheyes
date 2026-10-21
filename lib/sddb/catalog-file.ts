/**
 * Read a `catalog.parquet` release asset (Node only: hyparquet's file
 * buffer). Checks the file against the vendored contract before handing
 * rows back, so a changed upstream schema fails here, not in the database.
 */
import {
  READ_COLUMNS,
  REQUIRED_COLUMNS,
  type CatalogRow,
} from "./catalog-schema";

export interface CatalogFile {
  rows: CatalogRow[];
  /** `sddb_generated_at` from the file's metadata (a date), when present. */
  generatedAt: string | null;
  /** Known columns the file no longer carries (a warning, not an error). */
  missingColumns: string[];
}

export class CatalogSchemaError extends Error {}

export async function readCatalogFile(path: string): Promise<CatalogFile> {
  // hyparquet is ESM-only; a dynamic import loads its Node build (which adds
  // asyncBufferFromFile, declared in ./hyparquet-node.d.ts) whether the
  // caller is a tsx script (CommonJS) or vitest (ESM).
  const { asyncBufferFromFile, parquetMetadata, parquetReadObjects } =
    await import("hyparquet");
  const file = await asyncBufferFromFile(path);
  const meta = parquetMetadata(await file.slice(0, file.byteLength));
  const present = new Set(
    meta.schema.slice(1).map((el: { name: string }) => el.name),
  );
  const missingRequired = REQUIRED_COLUMNS.filter((c) => !present.has(c));

  if (missingRequired.length) {
    throw new CatalogSchemaError(
      `catalog file lacks required column(s): ${missingRequired.join(", ")}`,
    );
  }
  const missingColumns = READ_COLUMNS.filter((c) => !present.has(c));
  const columns = READ_COLUMNS.filter((c) => present.has(c));
  const generatedAt =
    meta.key_value_metadata?.find(
      (kv: { key: string }) => kv.key === "sddb_generated_at",
    )?.value ?? null;

  const rows = (await parquetReadObjects({
    file,
    columns: [...columns],
  })) as CatalogRow[];

  for (const c of REQUIRED_COLUMNS) {
    if (rows.length && rows.every((r) => r[c] === null || r[c] === undefined)) {
      throw new CatalogSchemaError(`required column ${c} is all-null`);
    }
  }

  return { rows, generatedAt, missingColumns: [...missingColumns] };
}
