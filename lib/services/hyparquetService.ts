// lib/services/hyparquetService.ts
import { parquetRead, type ColumnData } from "hyparquet";
import { compressors } from "hyparquet-compressors";

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
    if (typeof window === "undefined") {
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

    // Map to store accumulated column data as arrays of chunks
    // This avoids repeated array concatenation which creates many copies
    const columnChunks = new Map<string, any[][]>();

    columnNames.forEach((name) => columnChunks.set(name, []));

    let totalPages = 0;
    let relevantPages = 0;
    let lastReportedProgress = 10;

    // Read parquet file with onPage callback
    await parquetRead({
      file: arrayBuffer,
      compressors,
      onPage: (page: ColumnData) => {
        const columnName = page.columnName;

        // Only process requested columns
        if (columnNames.includes(columnName)) {
          totalPages++;
          relevantPages++;

          // Store chunks instead of concatenating - more memory efficient
          const chunks = columnChunks.get(columnName)!;

          chunks.push(Array.from(page.columnData));

          // Only report progress every 5% to reduce spam
          const progress =
            10 + Math.floor((relevantPages / (columnNames.length * 10)) * 20); // 10-30% estimate

          if (progress >= lastReportedProgress + 5) {
            onProgress?.(progress, `Reading parquet data...`);
            lastReportedProgress = progress;
          }
        }
        // Silently skip unwanted columns - no processing, no progress reporting
      },
    });

    await onProgress?.(25, "Combining column data...");

    // Flatten chunks into final arrays (done once at the end)
    const columnData = new Map<string, any[]>();

    for (const [columnName, chunks] of columnChunks.entries()) {
      columnData.set(columnName, chunks.flat());
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
