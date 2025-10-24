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

type UploadType = "h5ad" | "xenium" | "merscope";

// Map UploadType to MoleculeDatasetType for single molecule datasets
const UPLOAD_TYPE_TO_PARQUET_TYPE: Record<UploadType, MoleculeDatasetType> = {
  h5ad: "custom",
  xenium: "xenium",
  merscope: "merscope",
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
        // Single molecule mode - expects single parquet/csv file
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
      } else {
        // Single cell mode - original behavior
        if (type === "h5ad") {
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
            console.log("File:", h5adFile.name, "Size:", h5adFile.size, "bytes");
            dataset = await StandardizedDataset.fromH5ad(
              h5adFile,
              onProgress,
            );
          } else {
            toast.info(
              `Processing MERSCOPE folder (${files.length} files)...`,
            );
            console.log("=== Starting MERSCOPE folder processing ===");
            console.log("Files:", files.length);
            dataset = await StandardizedDataset.fromMerscope(
              files,
              onProgress,
            );
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
    !singleMolecule && (type === "xenium" || type === "merscope");

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
