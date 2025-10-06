"use client";

import { useCallback, useState } from "react";
import { useRouter } from "next/navigation";
import { toast } from "react-toastify";
import { StandardizedDataset } from "@/lib/StandardizedDataset";
import { useDatasetStore } from "@/lib/stores/datasetStore";

type UploadType = "h5ad" | "xenium" | "merscope";

interface FileUploadProps {
  type: UploadType;
  title: string;
  description: string;
}

export function FileUpload({ type, title, description }: FileUploadProps) {
  const [isDragging, setIsDragging] = useState(false);
  const [progress, setProgress] = useState(0);
  const [progressMessage, setProgressMessage] = useState("");
  const { addDataset, setLoading, setError } = useDatasetStore();
  const router = useRouter();

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

  const handleDrop = useCallback((e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);

    const files = Array.from(e.dataTransfer.files);
    handleFiles(files);
  }, [type]);

  const handleFileInput = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const files = e.target.files ? Array.from(e.target.files) : [];
      handleFiles(files);
    },
    [type]
  );

  const handleFiles = async (files: File[]) => {
    if (files.length === 0) return;

    try {
      setLoading(true);
      setError(null);
      setProgress(0);
      setProgressMessage("Starting...");

      const onProgress = async (prog: number, msg: string) => {
        setProgress(prog);
        setProgressMessage(msg);
        // Yield to allow React to re-render
        await new Promise(resolve => setTimeout(resolve, 0));
      };

      let dataset: StandardizedDataset;

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
        // merscope
        toast.info(`Processing MERSCOPE folder (${files.length} files)...`);
        console.log("=== Starting MERSCOPE folder processing ===");
        console.log("Files:", files.length);
        dataset = await StandardizedDataset.fromMerscope(files, onProgress);
      }

      console.log("=== Dataset created successfully ===");
      console.log("Dataset ID:", dataset.id);
      console.log("Dataset name:", dataset.name);
      console.log("Dataset type:", dataset.type);
      console.log("Point count:", dataset.getPointCount());
      console.log("Gene count:", dataset.genes.length);
      console.log("Spatial dimensions:", dataset.spatial.dimensions);
      console.log("Available embeddings:", Object.keys(dataset.embeddings));
      console.log("Cluster info:", dataset.clusters);
      console.log("Summary:", dataset.getSummary());

      // Add to Zustand store
      setProgress(100);
      setProgressMessage("Complete!");
      addDataset(dataset);

      toast.success(`Dataset loaded successfully!`);
      setLoading(false);

      // Navigate to viewer page
      router.push("/viewer");
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
    document.getElementById(`file-input-${type}`)?.click();
  };

  const isFolder = type === "xenium" || type === "merscope";
  const { isLoading } = useDatasetStore();

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
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onDrop={handleDrop}
        onClick={handleClick}
      >
        <input
          id={`file-input-${type}`}
          type="file"
          multiple={isFolder}
          className="hidden"
          onChange={handleFileInput}
          accept={type === "h5ad" ? ".h5ad" : ".csv,.tsv,.txt"}
          {...(isFolder
            ? { webkitdirectory: "", directory: "" }
            : {})}
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
