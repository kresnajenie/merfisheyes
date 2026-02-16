"use client";

import type { LocalDatasetMetadata } from "@/lib/services/localDatasetDB";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";

import { useCallback, useState } from "react";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
  Button,
  Chip,
  Progress,
} from "@heroui/react";
import { toast } from "react-toastify";
import * as Comlink from "comlink";

interface LocalDatasetReuploadModalProps {
  isOpen: boolean;
  onClose: () => void;
  metadata: LocalDatasetMetadata;
  expectedDatasetId: string;
  onDatasetLoaded: (
    dataset: StandardizedDataset | SingleMoleculeDataset,
  ) => void;
}

const FORMAT_LABELS: Record<string, string> = {
  h5ad: "H5AD",
  xenium: "Xenium",
  merscope: "MERSCOPE",
  parquet: "Parquet",
  csv: "CSV",
  processed_chunked: "Pre-chunked",
};

export function LocalDatasetReuploadModal({
  isOpen,
  onClose,
  metadata,
  expectedDatasetId,
  onDatasetLoaded,
}: LocalDatasetReuploadModalProps) {
  const [progress, setProgress] = useState(0);
  const [progressMessage, setProgressMessage] = useState("");
  const [isProcessing, setIsProcessing] = useState(false);

  const isFolder =
    metadata.type === "xenium" ||
    metadata.type === "merscope" ||
    metadata.type === "processed_chunked";

  const isMolecule = metadata.datasetCategory === "molecule";
  const inputId = "reupload-file-input";

  const onProgress = useCallback((prog: number, msg: string) => {
    setProgress(prog);
    setProgressMessage(msg);
  }, []);

  const processFiles = async (files: File[]) => {
    if (files.length === 0) return;

    try {
      setIsProcessing(true);
      setProgress(0);
      setProgressMessage("Starting...");

      let dataset: StandardizedDataset | SingleMoleculeDataset;

      if (isMolecule) {
        const { SingleMoleculeDataset } = await import(
          "@/lib/SingleMoleculeDataset"
        );

        if (metadata.type === "processed_chunked") {
          dataset = await SingleMoleculeDataset.fromLocalChunked(
            files,
            onProgress,
          );
        } else {
          const file = files[0];
          const ext = file.name.split(".").pop()?.toLowerCase();
          const { getSingleMoleculeWorker } = await import(
            "@/lib/workers/singleMoleculeWorkerManager"
          );
          const workerApi = await getSingleMoleculeWorker();
          const proxiedProgress = Comlink.proxy(onProgress);

          let serializedData;

          if (ext === "parquet") {
            serializedData = await workerApi.parseParquet(
              file,
              "custom",
              proxiedProgress,
            );
          } else if (ext === "csv") {
            serializedData = await workerApi.parseCSV(
              file,
              "custom",
              proxiedProgress,
            );
          } else {
            throw new Error(`Unsupported file type: .${ext}`);
          }

          dataset = SingleMoleculeDataset.fromSerializedData(serializedData);
        }
      } else {
        const { StandardizedDataset } = await import(
          "@/lib/StandardizedDataset"
        );

        if (metadata.type === "processed_chunked") {
          dataset = await StandardizedDataset.fromLocalChunked(
            files,
            onProgress,
          );
        } else if (metadata.type === "h5ad") {
          dataset = await StandardizedDataset.fromH5ad(files[0], onProgress);
        } else if (metadata.type === "xenium") {
          dataset = await StandardizedDataset.fromXenium(files, onProgress);
        } else if (metadata.type === "merscope") {
          // Check for H5AD inside MERSCOPE folder
          const h5adFile = files.find((f) => {
            const name = f.name.toLowerCase();
            const rel = (f.webkitRelativePath || "").toLowerCase();

            return name.endsWith(".h5ad") || rel.endsWith(".h5ad");
          });

          if (h5adFile) {
            dataset = await StandardizedDataset.fromH5ad(h5adFile, onProgress);
          } else {
            dataset = await StandardizedDataset.fromMerscope(files, onProgress);
          }
        } else {
          throw new Error(`Unknown dataset type: ${metadata.type}`);
        }
      }

      // Override dataset ID to match the URL
      (dataset as any).id = expectedDatasetId;

      toast.success("Dataset loaded successfully!");
      onDatasetLoaded(dataset);
    } catch (error) {
      console.error("Error re-processing dataset:", error);
      const msg = error instanceof Error ? error.message : "Unknown error";

      toast.error(`Failed to process data: ${msg}`);
    } finally {
      setIsProcessing(false);
      setProgress(0);
      setProgressMessage("");
    }
  };

  const handleFileInput = (e: React.ChangeEvent<HTMLInputElement>) => {
    const files = e.target.files ? Array.from(e.target.files) : [];

    processFiles(files);
  };

  const handleClick = () => {
    document.getElementById(inputId)?.click();
  };

  const stats: string[] = [];

  if (metadata.pointCount)
    stats.push(`${metadata.pointCount.toLocaleString()} cells`);
  if (metadata.moleculeCount)
    stats.push(`${metadata.moleculeCount.toLocaleString()} molecules`);
  if (metadata.geneCount) stats.push(`${metadata.geneCount} genes`);
  if (metadata.uniqueGeneCount) stats.push(`${metadata.uniqueGeneCount} genes`);
  if (metadata.spatialDimensions) stats.push(`${metadata.spatialDimensions}D`);

  const getAcceptedFileTypes = () => {
    if (isMolecule && !isFolder) return ".parquet,.csv";
    if (metadata.type === "h5ad") return ".h5ad";

    return ".csv,.tsv,.txt";
  };

  return (
    <Modal
      backdrop="blur"
      isDismissable={false}
      isOpen={isOpen}
      size="2xl"
      onClose={onClose}
    >
      <ModalContent>
        <ModalHeader className="flex flex-col gap-1">
          Dataset Not in Memory
        </ModalHeader>
        <ModalBody>
          <div className="flex flex-col gap-4">
            <p className="text-sm text-default-600">
              This dataset was previously loaded locally but is no longer in
              memory. Please re-select the original file to continue.
            </p>

            <div className="flex flex-col gap-2 p-4 bg-default-100 rounded-lg">
              <div className="flex flex-wrap items-center gap-2">
                <span className="text-sm font-semibold break-all">
                  {metadata.name}
                </span>
                <Chip color="primary" size="sm" variant="flat">
                  {FORMAT_LABELS[metadata.type] || metadata.type}
                </Chip>
                <Chip
                  color={isMolecule ? "secondary" : "success"}
                  size="sm"
                  variant="flat"
                >
                  {isMolecule ? "Molecule" : "Cell"}
                </Chip>
              </div>
              {stats.length > 0 && (
                <p className="text-xs text-default-500">{stats.join(" | ")}</p>
              )}
            </div>

            <div
              className={`
                border-2 border-dashed rounded-lg p-6 text-center cursor-pointer
                transition-all duration-200
                ${isProcessing ? "pointer-events-none opacity-60" : "border-default-300 hover:border-primary/50 hover:bg-default-100/50"}
              `}
              role="button"
              tabIndex={0}
              onClick={handleClick}
              onKeyDown={(e) => {
                if (e.key === "Enter" || e.key === " ") handleClick();
              }}
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
              <p className="text-sm text-default-600">
                {isFolder ? "Click to select folder" : "Click to select file"}
              </p>
              <p className="text-xs text-default-400 mt-1">
                {isFolder
                  ? `Select the original ${FORMAT_LABELS[metadata.type] || metadata.type} folder`
                  : `Select the original .${metadata.type} file`}
              </p>
            </div>

            {isProcessing && progress > 0 && (
              <div className="w-full">
                <Progress
                  aria-label="Processing progress"
                  className="w-full"
                  color="primary"
                  size="sm"
                  value={progress}
                />
                <p className="text-xs text-default-500 mt-1">
                  {progressMessage}
                </p>
              </div>
            )}
          </div>
        </ModalBody>
        <ModalFooter>
          <Button color="danger" variant="light" onPress={onClose}>
            Go Home
          </Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}
