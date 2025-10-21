"use client";

import { parquetService } from "@/lib/services/parquetService";

interface ParquetColumnMapping {
  gene: string;
  x: string;
  y: string;
  z: string;
}

interface MoleculeData {
  genes: string[];
  coordinates: number[];
  dimensions: 2 | 3;
}

export type ParquetDatasetType = "xenium" | "merscope" | "custom";

/**
 * Predefined column mappings for different dataset types
 */
const DATASET_COLUMN_MAPPINGS: Record<
  ParquetDatasetType,
  ParquetColumnMapping
> = {
  xenium: {
    gene: "feature_name",
    x: "x_location",
    y: "y_location",
    z: "z_location",
  },
  merscope: {
    gene: "gene",
    x: "global_x",
    y: "global_y",
    z: "global_z",
  },
  custom: {
    gene: "feature_name",
    x: "x_location",
    y: "y_location",
    z: "z_location",
  },
};

/**
 * Adapter for reading single molecule data from Parquet files
 */
export class ParquetAdapter {
  private data: any[] = [];
  private columnMapping: ParquetColumnMapping;
  private datasetType: ParquetDatasetType;

  constructor(datasetType: ParquetDatasetType = "xenium") {
    this.datasetType = datasetType;
    this.columnMapping = { ...DATASET_COLUMN_MAPPINGS[datasetType] };
  }

  /**
   * Initialize the adapter with a parquet or CSV file
   */
  async initialize(
    file: File,
    onProgress?: (progress: number, message: string) => Promise<void> | void
  ): Promise<void> {
    const fileExtension = file.name.split(".").pop()?.toLowerCase();

    if (fileExtension === "csv") {
      await this.initializeFromCSV(file, onProgress);
    } else {
      await this.initializeFromParquet(file, onProgress);
    }
  }

  /**
   * Initialize from a CSV file
   */
  private async initializeFromCSV(
    file: File,
    onProgress?: (progress: number, message: string) => Promise<void> | void
  ): Promise<void> {
    await onProgress?.(10, "Reading CSV file...");

    try {
      const text = await file.text();
      await onProgress?.(30, "Parsing CSV data...");

      const lines = text.trim().split("\n");
      if (lines.length === 0) {
        throw new Error("CSV file is empty");
      }

      // Parse header
      const header = lines[0].split(",").map((col) => col.trim());

      // Parse rows
      const rows: any[] = [];
      for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(",").map((val) => val.trim());
        const row: any = {};
        header.forEach((col, idx) => {
          row[col] = values[idx];
        });
        rows.push(row);
      }

      this.data = rows;

      await onProgress?.(45, `Loaded ${rows.length} molecules`);
    } catch (error) {
      throw new Error(`Failed to parse CSV file: ${error}`);
    }
  }

  /**
   * Initialize from a Parquet file using parquet-wasm
   */
  private async initializeFromParquet(
    file: File,
    onProgress?: (progress: number, message: string) => Promise<void> | void
  ): Promise<void> {
    await onProgress?.(10, "Reading parquet file...");

    try {
      await onProgress?.(30, "Parsing parquet data...");

      // Parse parquet file using parquet-wasm (returns Apache Arrow Table)
      const table = await parquetService.readParquet(file);

      // Convert Arrow Table to rows (objects)
      const rows: any[] = [];
      const numRows = table.numRows;

      // Get column names from schema
      const columnNames = table.schema.fields.map((field) => field.name);

      // Extract data column by column using Arrow API
      const columnData: Record<string, any[]> = {};
      for (const columnName of columnNames) {
        const column = table.getChild(columnName);
        if (column) {
          // Convert Arrow Vector to JS array
          columnData[columnName] = column.toArray();
        }
      }

      // Convert columnar data to row-based objects
      for (let i = 0; i < numRows; i++) {
        const row: any = {};
        for (const columnName of columnNames) {
          row[columnName] = columnData[columnName]?.[i];
        }
        rows.push(row);
      }

      this.data = rows;

      await onProgress?.(45, `Loaded ${rows.length} molecules`);
    } catch (error) {
      throw new Error(`Failed to parse parquet file: ${error}`);
    }
  }

  /**
   * Set custom column mapping for different data formats
   */
  setColumnMapping(mapping: Partial<ParquetColumnMapping>): void {
    this.columnMapping = { ...this.columnMapping, ...mapping };
  }

  /**
   * Load molecule data from the parquet file
   * Returns genes array and flattened coordinates array
   */
  async loadMoleculeData(): Promise<MoleculeData> {
    if (!this.data || this.data.length === 0) {
      throw new Error("No data loaded. Call initialize() first.");
    }

    const genes: string[] = [];
    const coordinates: number[] = [];
    let hasZCoordinate = false;

    // Check if z coordinate exists in the first row
    if (this.data.length > 0) {
      const firstRow = this.data[0];
      hasZCoordinate = this.columnMapping.z in firstRow;
    }

    const dimensions: 2 | 3 = hasZCoordinate ? 3 : 2;

    // Process each row (molecule)
    for (const row of this.data) {
      // Extract gene name
      const gene = row[this.columnMapping.gene];
      if (gene === undefined || gene === null) {
        console.warn("Row missing gene column:", row);
        continue;
      }
      genes.push(String(gene));

      // Extract coordinates
      const x = Number(row[this.columnMapping.x]);
      const y = Number(row[this.columnMapping.y]);
      const z = hasZCoordinate ? Number(row[this.columnMapping.z]) : 0;

      if (isNaN(x) || isNaN(y) || (hasZCoordinate && isNaN(z))) {
        console.warn("Row has invalid coordinates:", row);
        continue;
      }

      // Add to flat array: [x1,y1,z1, x2,y2,z2, ...]
      coordinates.push(x, y, z);
    }

    return {
      genes,
      coordinates,
      dimensions,
    };
  }

  /**
   * Get dataset information
   */
  getDatasetInfo() {
    return {
      moleculeCount: this.data.length,
      datasetType: this.datasetType,
      columnMapping: this.columnMapping,
    };
  }

  /**
   * Get the column mapping being used
   */
  getColumnMapping(): ParquetColumnMapping {
    return { ...this.columnMapping };
  }

  /**
   * Get the dataset type
   */
  getDatasetType(): ParquetDatasetType {
    return this.datasetType;
  }
}
