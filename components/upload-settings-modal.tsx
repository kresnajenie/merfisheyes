"use client";

import { useState } from "react";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
  Button,
  Input,
  Select,
  SelectItem,
} from "@heroui/react";
import { StandardizedDataset } from "@/lib/StandardizedDataset";
import { GeneChunkProcessor } from "@/lib/utils/GeneChunkProcessor";
import { toast } from "react-toastify";

interface UploadSettingsModalProps {
  isOpen: boolean;
  onClose: () => void;
  dataset: StandardizedDataset | null;
}

export function UploadSettingsModal({
  isOpen,
  onClose,
  dataset,
}: UploadSettingsModalProps) {
  const [chunkSize, setChunkSize] = useState<string>("auto");
  const [customChunkSize, setCustomChunkSize] = useState<string>("100");
  const [datasetName, setDatasetName] = useState<string>(
    dataset?.name || "dataset"
  );
  const [isProcessing, setIsProcessing] = useState(false);
  const [progress, setProgress] = useState(0);
  const [progressMessage, setProgressMessage] = useState("");

  const handleUpload = async () => {
    if (!dataset) {
      toast.error("No dataset loaded");
      return;
    }

    setIsProcessing(true);
    setProgress(0);
    setProgressMessage("Initializing...");

    try {
      // Determine chunk size
      const actualChunkSize =
        chunkSize === "auto" ? null : parseInt(customChunkSize);
      const processor = new GeneChunkProcessor(actualChunkSize);

      // Process genes
      setProgressMessage("Processing genes into chunks...");
      const { chunks, index } = await processor.processGenes(
        dataset,
        (prog, msg) => {
          setProgress(prog * 0.7); // Genes take 70% of progress
          setProgressMessage(msg);
        }
      );

      // Process coordinates
      setProgressMessage("Processing coordinates...");
      setProgress(70);
      const coordinates = await processor.processCoordinates(dataset);

      // Process metadata
      setProgressMessage("Processing metadata...");
      setProgress(85);
      const metadata = await processor.processMetadata(dataset);

      // Save all files
      setProgressMessage("Saving files...");
      setProgress(90);
      await saveFiles(datasetName, chunks, index, coordinates, metadata);

      setProgress(100);
      setProgressMessage("Complete!");
      toast.success("Dataset processed and saved successfully!");

      setTimeout(() => {
        onClose();
        setIsProcessing(false);
        setProgress(0);
        setProgressMessage("");
      }, 1000);
    } catch (error) {
      console.error("Upload processing error:", error);
      toast.error(
        `Failed to process dataset: ${error instanceof Error ? error.message : "Unknown error"}`
      );
      setIsProcessing(false);
      setProgress(0);
      setProgressMessage("");
    }
  };

  return (
    <Modal isOpen={isOpen} onClose={onClose} size="2xl">
      <ModalContent>
        <ModalHeader className="flex flex-col gap-1">
          Upload Dataset Settings
        </ModalHeader>
        <ModalBody>
          <div className="flex flex-col gap-4">
            <Input
              label="Dataset Name"
              placeholder="Enter dataset name"
              value={datasetName}
              onValueChange={setDatasetName}
              description="Name for the saved files"
              isDisabled={isProcessing}
            />

            <Select
              label="Chunk Size"
              placeholder="Select chunk size"
              selectedKeys={[chunkSize]}
              onSelectionChange={(keys) => {
                const value = Array.from(keys)[0] as string;
                setChunkSize(value);
              }}
              description="Number of genes per chunk (auto-determines based on total genes)"
              isDisabled={isProcessing}
            >
              <SelectItem key="auto">Auto (Recommended)</SelectItem>
              <SelectItem key="50">50 genes/chunk</SelectItem>
              <SelectItem key="100">100 genes/chunk</SelectItem>
              <SelectItem key="150">150 genes/chunk</SelectItem>
              <SelectItem key="custom">Custom</SelectItem>
            </Select>

            {chunkSize === "custom" && (
              <Input
                type="number"
                label="Custom Chunk Size"
                placeholder="Enter chunk size"
                value={customChunkSize}
                onValueChange={setCustomChunkSize}
                min="10"
                max="500"
                isDisabled={isProcessing}
              />
            )}

            {dataset && (
              <div className="bg-default-100 p-4 rounded-lg">
                <h4 className="text-sm font-semibold mb-2">Dataset Info</h4>
                <div className="text-xs space-y-1">
                  <p>
                    <span className="font-medium">Genes:</span>{" "}
                    {dataset.genes.length.toLocaleString()}
                  </p>
                  <p>
                    <span className="font-medium">Cells:</span>{" "}
                    {dataset.getPointCount().toLocaleString()}
                  </p>
                  <p>
                    <span className="font-medium">Type:</span> {dataset.type}
                  </p>
                  {chunkSize !== "custom" && (
                    <p>
                      <span className="font-medium">Estimated chunks:</span>{" "}
                      {Math.ceil(
                        dataset.genes.length /
                          new GeneChunkProcessor(
                            chunkSize === "auto"
                              ? null
                              : parseInt(chunkSize)
                          ).determineChunkSize(dataset.genes.length)
                      )}
                    </p>
                  )}
                </div>
              </div>
            )}

            {isProcessing && (
              <div className="space-y-2">
                <div className="w-full bg-default-200 rounded-full h-2 overflow-hidden">
                  <div
                    className="bg-primary h-full transition-all duration-300"
                    style={{ width: `${progress}%` }}
                  />
                </div>
                <p className="text-xs text-default-500 text-center">
                  {progressMessage}
                </p>
              </div>
            )}
          </div>
        </ModalBody>
        <ModalFooter>
          <Button
            color="danger"
            variant="light"
            onPress={onClose}
            isDisabled={isProcessing}
          >
            Cancel
          </Button>
          <Button
            color="primary"
            onPress={handleUpload}
            isLoading={isProcessing}
            isDisabled={!dataset}
          >
            {isProcessing ? "Processing..." : "Process & Save"}
          </Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}

/**
 * Save all processed files to disk
 */
async function saveFiles(
  datasetName: string,
  chunks: any[],
  index: any,
  coordinates: Record<string, Blob>,
  metadata: Blob
) {
  // Create a folder name
  const folderName = `${datasetName}_${Date.now()}`;

  // Save gene chunks
  for (const chunk of chunks) {
    await saveBlob(chunk.data, `${folderName}/genes/${chunk.filename}`);
  }

  // Save gene index
  const indexBlob = new Blob([JSON.stringify(index, null, 2)], {
    type: "application/json",
  });
  await saveBlob(indexBlob, `${folderName}/genes/index.json`);

  // Save coordinates
  for (const [name, blob] of Object.entries(coordinates)) {
    await saveBlob(blob, `${folderName}/coordinates/${name}.bin.gz`);
  }

  // Save metadata
  await saveBlob(metadata, `${folderName}/metadata.json.gz`);

  toast.info(
    `Files saved to downloads folder: ${folderName} (${chunks.length + Object.keys(coordinates).length + 2} files)`
  );
}

/**
 * Save a blob to disk using the browser download API
 */
async function saveBlob(blob: Blob, filepath: string) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filepath.replace(/\//g, "_"); // Flatten path for browser download
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);

  // Small delay to prevent browser from blocking multiple downloads
  await new Promise((resolve) => setTimeout(resolve, 100));
}
