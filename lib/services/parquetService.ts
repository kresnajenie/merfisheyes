// lib/services/parquetService.ts
import * as arrow from "apache-arrow";
import { parseTable } from "arrow-js-ffi";

type ParquetModule = typeof import("parquet-wasm");

class ParquetService {
  private static instance: ParquetService;
  private parquetModule: ParquetModule | null = null;
  private initPromise: Promise<ParquetModule> | null = null;
  private wasmInitialized = false;

  private constructor() {}

  static getInstance(): ParquetService {
    if (!ParquetService.instance) {
      ParquetService.instance = new ParquetService();
    }
    return ParquetService.instance;
  }

  private async initialize(): Promise<ParquetModule> {
    if (typeof window === "undefined") {
      throw new Error("Parquet service can only be used in browser");
    }

    if (this.parquetModule && this.wasmInitialized) {
      return this.parquetModule;
    }

    if (this.initPromise) {
      return this.initPromise;
    }

    this.initPromise = import("parquet-wasm").then(async (module) => {
      // Initialize the WebAssembly module
      await module.default();
      this.wasmInitialized = true;
      this.parquetModule = module;
      return module;
    });

    return this.initPromise;
  }

  // Read parquet and return Arrow Table
  async readParquet(input: File | ArrayBuffer | Uint8Array): Promise<arrow.Table> {
    const module = await this.initialize();

    // Convert input to Uint8Array
    let uint8Array: Uint8Array;
    if (input instanceof File) {
      const buffer = await input.arrayBuffer();
      uint8Array = new Uint8Array(buffer);
    } else if (input instanceof ArrayBuffer) {
      uint8Array = new Uint8Array(input);
    } else {
      uint8Array = input;
    }

    // Read parquet file using parquet-wasm
    const wasmArrowTable = module.readParquet(uint8Array).intoFFI();

    // Get reference to WASM memory
    const wasmMemory = module.wasmMemory();

    // Parse the Arrow table from WASM memory
    const table: arrow.Table = parseTable(
      wasmMemory.buffer,
      wasmArrowTable.arrayAddrs(),
      wasmArrowTable.schemaAddr()
    );

    // IMPORTANT: Release WASM memory
    wasmArrowTable.drop();

    return table;
  }
}

export const parquetService = ParquetService.getInstance();
