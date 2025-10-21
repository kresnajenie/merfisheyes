// lib/services/hyparquetService.ts
import { parquetRead } from "hyparquet";
import { compressors } from "hyparquet-compressors";

/**
 * Column data structure returned by onPage callback from hyparquet
 */
interface ColumnData {
  columnName: string;
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
    onProgress?: (progress: number, message: string) => Promise<void> | void
  ): Promise<Map<string, any[]>> {
    if (typeof window === "undefined") {
      throw new Error("Hyparquet service can only be used in browser");
    }

    await onProgress?.(5, "Reading file into memory...");
    console.log(columnNames);

    // Convert File to ArrayBuffer (hyparquet accepts ArrayBuffer as AsyncBuffer in browser)
    const arrayBuffer = await file.arrayBuffer();

    await onProgress?.(10, "Parsing parquet structure...");

    // Map to store accumulated column data
    const columnData = new Map<string, any[]>();
    columnNames.forEach((name) => columnData.set(name, []));

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

          const existingData = columnData.get(columnName) || [];

          // Concatenate page data (columnData is ArrayLike)
          // Use concat instead of spread to avoid "Maximum call stack size exceeded"
          const newData = Array.from(page.columnData);
          const combined = existingData.concat(newData);

          columnData.set(columnName, combined);

          // Only report progress every 5% to reduce spam
          const progress = 10 + Math.floor((relevantPages / (columnNames.length * 10)) * 20); // 10-30% estimate
          if (progress >= lastReportedProgress + 5) {
            onProgress?.(progress, `Reading parquet data...`);
            lastReportedProgress = progress;
          }
        }
        // Silently skip unwanted columns - no processing, no progress reporting
      },
    });

    await onProgress?.(30, "Column extraction complete");

    // Validate all requested columns were found
    for (const columnName of columnNames) {
      const data = columnData.get(columnName);
      if (!data || data.length === 0) {
        throw new Error(
          `Column '${columnName}' not found or empty in parquet file`
        );
      }
    }

    return columnData;
  }
}

export const hyparquetService = HyparquetService.getInstance();
