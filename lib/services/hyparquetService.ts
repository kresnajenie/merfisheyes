// lib/services/hyparquetService.ts
import { parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";

/**
 * Column data structure from hyparquet onPage callback
 * Note: We define this locally to ensure type compatibility across environments
 * as hyparquet's exported types can vary between TypeScript versions
 */
interface ColumnData {
  pathInSchema: string[];
  columnData: ArrayLike<any>;
  rowStart: number;
  rowEnd: number;
}

/**
 * Service for reading parquet files using hyparquet
 */
class HyparquetService {
  private static instance: HyparquetService;

  private constructor() {}

  static getInstance(): HyparquetService {
    if (!HyparquetService.instance) {
      HyparquetService.instance = new HyparquetService();
    }

    return HyparquetService.instance;
  }

  /**
   * Read parquet file and collect column data via onPage callback
   * Returns a Map of column name to concatenated data arrays
   */
  async readParquetColumns(
    file: File,
    columnNames: string[],
    onProgress?: (progress: number, message: string) => Promise<void> | void,
  ): Promise<Map<string, any[]>> {
    if (typeof self === "undefined") {
      throw new Error("Hyparquet service can only be used in browser");
    }

    // Warn for very large files
    const fileGB = file.size / 1_000_000_000;

    if (fileGB > 2) {
      console.warn(
        `⚠️ Large file detected: ${fileGB.toFixed(1)}GB - This may cause memory issues`,
      );
    }

    await onProgress?.(5, "Reading file into memory...");
    console.log(columnNames);

    // Convert File to ArrayBuffer (hyparquet accepts ArrayBuffer as AsyncBuffer in browser)
    // Note: This loads entire file into memory - may fail for very large files (>2-3GB)
    const arrayBuffer = await file.arrayBuffer();

    await onProgress?.(10, "Parsing parquet structure...");

    // Store raw chunk references per column — no copying until the end
    const columnChunks = new Map<string, ArrayLike<any>[]>();
    const columnLengths = new Map<string, number>();

    columnNames.forEach((name) => {
      columnChunks.set(name, []);
      columnLengths.set(name, 0);
    });

    let relevantPages = 0;
    let lastReportedProgress = 10;

    // Read parquet file with columns filter + onPage callback
    // The `columns` option tells hyparquet to skip decompressing/parsing unwanted columns
    await parquetRead({
      file: arrayBuffer,
      columns: columnNames,
      compressors,
      onPage: ((page: ColumnData) => {
        const columnName = page.pathInSchema[0];
        const chunks = columnChunks.get(columnName);

        if (chunks) {
          relevantPages++;
          // Store reference to raw data — avoid Array.from() copy
          chunks.push(page.columnData);
          columnLengths.set(
            columnName,
            columnLengths.get(columnName)! + page.columnData.length,
          );

          const progress =
            10 + Math.floor((relevantPages / (columnNames.length * 10)) * 20);

          if (progress >= lastReportedProgress + 5) {
            onProgress?.(progress, `Reading parquet data...`);
            lastReportedProgress = progress;
          }
        }
      }) as any,
    });

    await onProgress?.(25, "Combining column data...");

    // Concatenate chunks into final arrays (single allocation + copy)
    const columnData = new Map<string, any[]>();

    for (const [columnName, chunks] of columnChunks.entries()) {
      const totalLen = columnLengths.get(columnName)!;
      const result = new Array(totalLen);
      let offset = 0;

      for (const chunk of chunks) {
        for (let i = 0; i < chunk.length; i++) {
          result[offset++] = chunk[i];
        }
      }
      columnData.set(columnName, result);
    }

    await onProgress?.(30, "Column extraction complete");

    // Validate all requested columns were found
    for (const columnName of columnNames) {
      const data = columnData.get(columnName);

      if (!data || data.length === 0) {
        throw new Error(
          `Column '${columnName}' not found or empty in parquet file`,
        );
      }
    }

    return columnData;
  }
}

export const hyparquetService = HyparquetService.getInstance();
