"use client";

import { useCallback, useState } from "react";
import { useRouter } from "next/navigation";
import { toast } from "react-toastify";
import * as Comlink from "comlink";

import { StandardizedDataset } from "@/lib/StandardizedDataset";
import { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import { MoleculeDatasetType } from "@/lib/config/moleculeColumnMappings";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { getSingleMoleculeWorker } from "@/lib/workers/singleMoleculeWorkerManager";

type UploadType = "h5ad" | "xenium" | "merscope" | "chunked";

// Map UploadType to MoleculeDatasetType for single molecule datasets
const UPLOAD_TYPE_TO_PARQUET_TYPE: Record<UploadType, MoleculeDatasetType> = {
  h5ad: "custom",
  xenium: "xenium",
  merscope: "merscope",
  chunked: "custom",
};

interface FileUploadProps {
  type: UploadType;
  title: string;
  description: string;
  singleMolecule?: boolean;
}

export function FileUpload({
  type,
  title,
  description,
  singleMolecule = false,
}: FileUploadProps) {
  const [isDragging, setIsDragging] = useState(false);
  const [progress, setProgress] = useState(0);
  const [progressMessage, setProgressMessage] = useState("");

  // Use appropriate store based on singleMolecule mode
  const cellStore = useDatasetStore();
  const smStore = useSingleMoleculeStore();

  const { setLoading, setError } = singleMolecule ? smStore : cellStore;
  const router = useRouter();

  const inputId = `${singleMolecule ? "sm" : "sc"}-file-input-${type}`;

  const handleDragOver = useCallback((e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(true);
  }, []);

  const handleDragLeave = useCallback((e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
  }, []);

  const handleDrop = useCallback(
    (e: React.DragEvent<HTMLDivElement>) => {
      e.preventDefault();
      e.stopPropagation();
      setIsDragging(false);

      const files = Array.from(e.dataTransfer.files);

      handleFiles(files);
    },
    [type],
  );

  const handleFileInput = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const files = e.target.files ? Array.from(e.target.files) : [];

      handleFiles(files);
    },
    [type],
  );

  /**
   * Detect if files represent a chunked dataset folder (from Python script)
   */
  const isChunkedDatasetFolder = (files: File[]): boolean => {
    const fileNames = files.map((f) => {
      const path = f.webkitRelativePath || f.name;
      // Remove root folder name
      const parts = path.split("/");

      return parts.length > 1 ? parts.slice(1).join("/") : path;
    });

    // Check for required chunked dataset files
    const hasManifest = fileNames.some((name) => name === "manifest.json");
    const hasExprIndex = fileNames.some((name) => name === "expr/index.json");
    const hasChunks = fileNames.some((name) =>
      /^expr\/chunk_\d{5}\.bin\.gz$/.test(name),
    );
    const hasSpatial = fileNames.some(
      (name) => name === "coords/spatial.bin.gz",
    );

    return hasManifest && hasExprIndex && hasChunks && hasSpatial;
  };

  /**
   * Detect if files represent a chunked single molecule dataset folder (from Python script)
   */
  const isChunkedSingleMoleculeFolder = (files: File[]): boolean => {
    const fileNames = files.map((f) => {
      const path = f.webkitRelativePath || f.name;
      // Remove root folder name
      const parts = path.split("/");

      return parts.length > 1 ? parts.slice(1).join("/") : path;
    });

    // Check for required chunked single molecule files
    const hasManifest = fileNames.some((name) => name === "manifest.json.gz");
    const hasGenesFolder = fileNames.some((name) => name.startsWith("genes/"));
    const hasGeneBinFiles = fileNames.some((name) =>
      /^genes\/[^\/]+\.bin\.gz$/.test(name),
    );

    return hasManifest && hasGenesFolder && hasGeneBinFiles;
  };

  const handleFiles = async (files: File[]) => {
    if (files.length === 0) return;

    try {
      setLoading(true);
      setError(null);
      setProgress(0);
      setProgressMessage("Starting...");

      const onProgress = (prog: number, msg: string) => {
        console.log(`[FileUpload] Progress: ${prog}% - ${msg}`);
        setProgress(prog);
        setProgressMessage(msg);
      };

      let dataset: StandardizedDataset | SingleMoleculeDataset;

      if (singleMolecule) {
        // Single molecule mode
        if (type === "chunked") {
          // Chunked folder upload - verify it's a valid chunked single molecule dataset
          if (!isChunkedSingleMoleculeFolder(files)) {
            throw new Error(
              "Invalid chunked single molecule dataset folder. Make sure it contains manifest.json.gz and genes/*.bin.gz files",
            );
          }

          toast.info(
            `Loading pre-chunked single molecule dataset (${files.length} files)...`,
          );
          console.log(
            "=== Loading pre-chunked single molecule dataset (ready for upload) ===",
          );
          console.log("Files:", files.length);

          onProgress(10, "Creating dataset from chunked files...");

          // Load dataset using ProcessedSingleMoleculeAdapter
          dataset = await SingleMoleculeDataset.fromLocalChunked(
            files,
            onProgress,
          );

          // Mark dataset as pre-chunked
          (dataset as any).isPreChunked = true;

          console.log("=== Pre-chunked Single Molecule Dataset loaded ===");
          console.log("Dataset ID:", dataset.id);
          console.log("Molecule count:", dataset.getMoleculeCount());
          console.log("Gene count:", dataset.genes.length);
          console.log("Spatial dimensions:", dataset.dimensions);
          console.log("Ready for upload!");
        } else {
          // Single parquet/csv file upload
          const file = files[0];
          const fileExtension = file.name.split(".").pop()?.toLowerCase();

          toast.info(`Processing ${file.name}...`);
          console.log(
            "=== Starting Single Molecule file processing (in worker) ===",
          );
          console.log("File:", file.name, "Size:", file.size, "bytes");
          console.log("File extension:", fileExtension);
          console.log(
            "Dataset type:",
            type,
            "→",
            UPLOAD_TYPE_TO_PARQUET_TYPE[type],
          );

          const parquetDatasetType = UPLOAD_TYPE_TO_PARQUET_TYPE[type];

          // Get singleton worker instance
          const workerApi = await getSingleMoleculeWorker();

          // Wrap progress callback with Comlink.proxy for cross-thread communication
          const proxiedProgress = Comlink.proxy(onProgress);

          // Parse file in web worker
          let serializedData;

          if (fileExtension === "parquet") {
            serializedData = await workerApi.parseParquet(
              file,
              parquetDatasetType,
              proxiedProgress,
            );
          } else if (fileExtension === "csv") {
            serializedData = await workerApi.parseCSV(
              file,
              parquetDatasetType,
              proxiedProgress,
            );
          } else {
            throw new Error(
              `Unsupported file type: .${fileExtension}. Only .parquet and .csv files are supported.`,
            );
          }

          // Reconstruct dataset from serialized data
          dataset = SingleMoleculeDataset.fromSerializedData(serializedData);

          console.log(
            "=== Single Molecule Dataset created successfully (from worker) ===",
          );
          console.log("Dataset ID:", dataset.id);
          console.log("Dataset name:", dataset.name);
          console.log("Dataset type:", dataset.type);
          console.log("Molecule count:", dataset.getMoleculeCount());
          console.log("Gene count:", dataset.genes.length);
          console.log("Spatial dimensions:", dataset.dimensions);
          console.log("Scaling factor:", dataset.scalingFactor);
          console.log("Summary:", dataset.getSummary());
        }
      } else {
        // Single cell mode
        if (type === "chunked") {
          // Chunked folder upload - verify it's a valid chunked dataset
          if (!isChunkedDatasetFolder(files)) {
            throw new Error(
              "Invalid chunked dataset folder. Make sure it contains manifest.json, expr/index.json, expr/chunk_*.bin.gz, and coords/spatial.bin.gz",
            );
          }
          toast.info(`Loading pre-chunked dataset (${files.length} files)...`);
          console.log("=== Loading pre-chunked dataset (ready for upload) ===");
          console.log("Files:", files.length);

          // Read manifest to get metadata
          onProgress(10, "Reading manifest...");
          const manifestFile = files.find((f) => {
            const path = f.webkitRelativePath || f.name;
            const parts = path.split("/");
            const relativePath =
              parts.length > 1 ? parts.slice(1).join("/") : path;

            return relativePath === "manifest.json";
          });

          if (!manifestFile) {
            throw new Error("manifest.json not found in chunked folder");
          }

          const manifestText = await manifestFile.text();
          const manifest = JSON.parse(manifestText);

          console.log("Manifest loaded:", manifest);

          // Create file map for later upload
          const fileMap = new Map<string, File>();

          for (const file of files) {
            const relativePath = file.webkitRelativePath;
            const parts = relativePath.split("/");
            const fileKey = parts.slice(1).join("/"); // Remove root folder

            fileMap.set(fileKey, file);
          }

          onProgress(30, "Creating dataset from manifest...");

          // Create a minimal StandardizedDataset from manifest
          // We don't actually load the chunked data, just create metadata
          dataset = await StandardizedDataset.fromLocalChunked(
            files,
            onProgress,
          );

          // Mark dataset as pre-chunked and attach file map
          (dataset as any).isPreChunked = true;
          (dataset as any).chunkedFiles = fileMap;
          (dataset as any).manifestData = manifest;

          onProgress(100, "Pre-chunked dataset ready for upload!");
          console.log(
            "Pre-chunked dataset ready. File map size:",
            fileMap.size,
          );
        } else if (type === "h5ad") {
          const file = files[0];

          toast.info(`Processing ${file.name}...`);
          console.log("=== Starting H5AD file processing ===");
          console.log("File:", file.name, "Size:", file.size, "bytes");
          dataset = await StandardizedDataset.fromH5ad(file, onProgress);
        } else if (type === "xenium") {
          toast.info(`Processing Xenium folder (${files.length} files)...`);
          console.log("=== Starting Xenium folder processing ===");
          console.log("Files:", files.length);
          dataset = await StandardizedDataset.fromXenium(files, onProgress);
        } else {
          // MERSCOPE uploads may contain a preprocessed H5AD
          const h5adFile = files.find((f) => {
            const name = f.name.toLowerCase();
            const rel = (f.webkitRelativePath || "").toLowerCase();

            return name.endsWith(".h5ad") || rel.endsWith(".h5ad");
          });

          if (h5adFile) {
            toast.info(
              `Detected H5AD file (${h5adFile.name}); processing as H5AD...`,
            );
            console.log(
              "=== H5AD detected inside MERSCOPE upload; routing to H5AD parser ===",
            );
            console.log(
              "File:",
              h5adFile.name,
              "Size:",
              h5adFile.size,
              "bytes",
            );
            dataset = await StandardizedDataset.fromH5ad(h5adFile, onProgress);
          } else {
            toast.info(`Processing MERSCOPE folder (${files.length} files)...`);
            console.log("=== Starting MERSCOPE folder processing ===");
            console.log("Files:", files.length);
            dataset = await StandardizedDataset.fromMerscope(files, onProgress);
          }
        }

        console.log("=== Dataset created successfully ===");
        console.log("Dataset ID:", dataset.id);
        console.log("Dataset name:", dataset.name);
        console.log("Dataset type:", dataset.type);
        console.log(
          "Point count:",
          (dataset as StandardizedDataset).getPointCount(),
        );
        console.log("Gene count:", dataset.genes.length);
        console.log(
          "Spatial dimensions:",
          (dataset as StandardizedDataset).spatial.dimensions,
        );
        console.log(
          "Available embeddings:",
          Object.keys((dataset as StandardizedDataset).embeddings),
        );
        console.log("Cluster info:", (dataset as StandardizedDataset).clusters);
        console.log("Summary:", dataset.getSummary());
      }

      // Add to appropriate Zustand store
      setProgress(100);
      setProgressMessage("Complete!");

      if (singleMolecule) {
        smStore.addDataset(dataset as SingleMoleculeDataset);
      } else {
        cellStore.addDataset(dataset as StandardizedDataset);
      }

      toast.success(`Dataset loaded successfully!`);
      setLoading(false);

      // Navigate to appropriate viewer page
      router.push(singleMolecule ? "/sm-viewer" : "/viewer");
    } catch (error) {
      console.error(`Error processing ${type} data:`, error);
      const errorMessage =
        error instanceof Error ? error.message : "Unknown error";

      setError(errorMessage);
      toast.error(`Failed to process data: ${errorMessage}`);
      setLoading(false);
      setProgress(0);
      setProgressMessage("");
    }
  };

  const handleClick = () => {
    document.getElementById(inputId)?.click();
  };

  const isFolder =
    (type === "xenium" || type === "merscope" || type === "chunked") &&
    (singleMolecule ? type === "chunked" : true);

  // Get isLoading from appropriate store
  const cellIsLoading = useDatasetStore((state) => state.isLoading);
  const smIsLoading = useSingleMoleculeStore((state) => state.isLoading);
  const isLoading = singleMolecule ? smIsLoading : cellIsLoading;

  // Determine accepted file types
  const getAcceptedFileTypes = () => {
    if (singleMolecule) {
      return ".parquet,.csv";
    }

    return type === "h5ad" ? ".h5ad" : ".csv,.tsv,.txt";
  };

  return (
    <div className="w-full">
      <div
        className={`
          relative border-2 border-dashed rounded-lg p-6 text-center cursor-pointer
          transition-all duration-200 ease-in-out aspect-square flex items-center justify-center
          ${
            isDragging
              ? "border-primary bg-primary/10 scale-105"
              : "border-default-300 hover:border-primary/50 hover:bg-default-100/50"
          }
          ${isLoading ? "pointer-events-none opacity-60" : ""}
        `}
        onClick={handleClick}
        onDragLeave={handleDragLeave}
        onDragOver={handleDragOver}
        onDrop={handleDrop}
      >
        <input
          accept={getAcceptedFileTypes()}
          className="hidden"
          id={inputId}
          multiple={isFolder}
          type="file"
          onChange={handleFileInput}
          {...(isFolder ? { webkitdirectory: "", directory: "" } : {})}
        />

        <div className="flex flex-col items-center gap-2 w-full">
          <p className="text-lg font-semibold text-foreground">{title}</p>
          <p className="text-xs text-default-500">{description}</p>

          {isLoading && progress > 0 && (
            <div className="w-full mt-4 px-4">
              <div className="w-full bg-default-200 rounded-full h-2 overflow-hidden">
                <div
                  className="bg-primary h-full transition-all duration-300 ease-out"
                  style={{ width: `${progress}%` }}
                />
              </div>
              <p className="text-xs text-default-500 mt-2">{progressMessage}</p>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
