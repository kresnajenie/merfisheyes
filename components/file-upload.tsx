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

      let dataset: StandardizedDataset;

      if (type === "h5ad") {
        const file = files[0];
        toast.info(`Processing ${file.name}...`);
        console.log("=== Starting H5AD file processing ===");
        console.log("File:", file.name, "Size:", file.size, "bytes");
        dataset = await StandardizedDataset.fromH5ad(file);
      } else if (type === "xenium") {
        toast.info(`Processing Xenium folder (${files.length} files)...`);
        console.log("=== Starting Xenium folder processing ===");
        console.log("Files:", files.length);
        dataset = await StandardizedDataset.fromXenium(files);
      } else {
        // merscope
        toast.info(`Processing MERSCOPE folder (${files.length} files)...`);
        console.log("=== Starting MERSCOPE folder processing ===");
        console.log("Files:", files.length);
        dataset = await StandardizedDataset.fromMerscope(files);
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
    }
  };

  const handleClick = () => {
    document.getElementById(`file-input-${type}`)?.click();
  };

  const isFolder = type === "xenium" || type === "merscope";

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

        <div className="flex flex-col items-center gap-2">
          <p className="text-lg font-semibold text-foreground">{title}</p>
          <p className="text-xs text-default-500">{description}</p>
        </div>
      </div>
    </div>
  );
}
