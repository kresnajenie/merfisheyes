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
import { toast } from "react-toastify";

import { StandardizedDataset } from "@/lib/StandardizedDataset";
import { GeneChunkProcessor } from "@/lib/utils/GeneChunkProcessor";
import { generateDatasetFingerprint } from "@/lib/utils/fingerprint";

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
    dataset?.name || "dataset",
  );
  const [email, setEmail] = useState<string>("");
  const [isProcessing, setIsProcessing] = useState(false);
  const [progress, setProgress] = useState(0);
  const [progressMessage, setProgressMessage] = useState("");
  const [uploadProgress, setUploadProgress] = useState(0);
  const [uploadMessage, setUploadMessage] = useState("");
  const [uploadComplete, setUploadComplete] = useState(false);
  const [uploadedDatasetId, setUploadedDatasetId] = useState<string>("");

  // Check if dataset is pre-chunked
  const isPreChunked = (dataset as any)?.isPreChunked || false;

  // Email validation
  const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
  const isEmailValid = emailRegex.test(email);

  const handleUpload = async () => {
    if (!dataset) {
      toast.error("No dataset loaded");

      return;
    }

    if (!email || !isEmailValid) {
      toast.error("Please provide a valid email address");

      return;
    }

    setIsProcessing(true);
    setProgress(0);
    setProgressMessage("Generating fingerprint...");

    try {
      let filesToUpload: Array<{ key: string; blob: Blob; size: number; contentType: string }>;
      let datasetId: string;
      let fingerprint: string;

      if (isPreChunked) {
        // Pre-chunked dataset - skip processing, use existing files
        console.log("=== Uploading pre-chunked dataset ===");

        const chunkedFiles = (dataset as any).chunkedFiles as Map<string, File>;
        const manifest = (dataset as any).manifestData;

        if (!chunkedFiles || !manifest) {
          throw new Error("Pre-chunked dataset is missing file map or manifest data");
        }

        // Generate fingerprint from manifest (simpler for pre-chunked)
        // We don't want to load all expression data just for fingerprinting
        const manifestString = JSON.stringify({
          cells: manifest.statistics.total_cells,
          genes: manifest.statistics.total_genes,
          type: dataset.type,
          genes_list: dataset.genes.sort().join(","),
          clusters: dataset.clusters?.map((c) => c.column).sort().join(","),
        });

        // Hash the manifest string to create a proper fingerprint
        const encoder = new TextEncoder();
        const data = encoder.encode(manifestString);
        const hashBuffer = await crypto.subtle.digest("SHA-256", data);
        const hashArray = Array.from(new Uint8Array(hashBuffer));
        fingerprint = hashArray.map((b) => b.toString(16).padStart(2, "0")).join("");

        console.log("Pre-chunked dataset fingerprint:", fingerprint.substring(0, 50));

        // Check for duplicates
        setProgressMessage("Checking for duplicates...");
        setProgress(5);
        const duplicateExists = await checkDuplicate(fingerprint);

        if (duplicateExists) {
          const confirmed = window.confirm(
            `A dataset with this content already exists.\n\nDataset: ${duplicateExists.title}\nUploaded: ${new Date(duplicateExists.createdAt).toLocaleString()}\n\nDo you want to upload anyway?`,
          );

          if (!confirmed) {
            setIsProcessing(false);
            setProgress(0);
            setProgressMessage("");

            return;
          }
        }

        // Generate dataset ID
        datasetId = `${dataset.type}_${datasetName}_${Date.now()}_${fingerprint.substring(0, 9)}`;

        // Update manifest with new dataset ID and name
        const updatedManifest = {
          ...manifest,
          dataset_id: datasetId,
          name: datasetName,
        };

        // Prepare files from the chunked folder
        setProgressMessage("Preparing pre-chunked files for upload...");
        setProgress(20);
        filesToUpload = [];

        // Add manifest
        const manifestBlob = new Blob([JSON.stringify(updatedManifest, null, 2)], { type: "application/json" });
        filesToUpload.push({
          key: "manifest.json",
          blob: manifestBlob,
          size: manifestBlob.size,
          contentType: "application/json",
        });

        // Add all other files from the chunked folder
        for (const [fileKey, file] of chunkedFiles.entries()) {
          if (fileKey !== "manifest.json") {
            filesToUpload.push({
              key: fileKey,
              blob: file,
              size: file.size,
              contentType: file.type || "application/octet-stream",
            });
          }
        }

        console.log(`Prepared ${filesToUpload.length} pre-chunked files for upload`);
        setProgress(55);
      } else {
        // Regular dataset - process and chunk
        // Generate fingerprint
        fingerprint = await generateDatasetFingerprint(dataset);

        console.log("Dataset fingerprint:", fingerprint);

        // Check for duplicates
        setProgressMessage("Checking for duplicates...");
        setProgress(5);
        const duplicateExists = await checkDuplicate(fingerprint);

        if (duplicateExists) {
          const confirmed = window.confirm(
            `A dataset with this content already exists.\n\nDataset: ${duplicateExists.title}\nUploaded: ${new Date(duplicateExists.createdAt).toLocaleString()}\n\nDo you want to upload anyway?`,
          );

          if (!confirmed) {
            setIsProcessing(false);
            setProgress(0);
            setProgressMessage("");

            return;
          }
        }

        // Determine chunk size
        const actualChunkSize =
          chunkSize === "auto" ? null : parseInt(customChunkSize);
        const processor = new GeneChunkProcessor(actualChunkSize);

        // Process genes
        setProgressMessage("Processing genes into chunks...");
        const { chunks, index } = await processor.processGenes(
          dataset,
          (prog, msg) => {
            setProgress(5 + prog * 0.35); // Genes take 35% of progress (5-40%)
            setProgressMessage(msg);
          },
        );

        // Process coordinates
        setProgressMessage("Processing coordinates...");
        setProgress(40);
        const coordinates = await processor.processCoordinates(dataset);

        // Process observations
        setProgressMessage("Processing observations...");
        setProgress(45);
        const observations = await processor.processObservations(dataset);

        // Process palettes
        setProgressMessage("Processing color palettes...");
        setProgress(47);
        const palettes = await processor.processPalettes(dataset);

        // Generate dataset ID
        datasetId = `${dataset.type}_${datasetName}_${Date.now()}_${fingerprint.substring(0, 9)}`;

        // Create manifest
        setProgressMessage("Creating manifest...");
        setProgress(50);
        const manifestJson = await createManifest(
          dataset,
          datasetId,
          datasetName,
          chunks,
          index,
          coordinates,
          observations.metadata,
        );

        // Prepare files for upload
        setProgressMessage("Preparing files for upload...");
        setProgress(55);
        filesToUpload = await prepareFilesForUpload(
          chunks,
          index,
          coordinates,
          observations.files,
          observations.metadata,
          palettes,
          manifestJson,
        );
      }

      // Initiate upload
      setProgressMessage("Initiating upload...");
      setProgress(60);
      const uploadSession = await initiateUpload(
        fingerprint,
        datasetName,
        dataset,
        filesToUpload,
      );

      // Upload files to S3
      setUploadMessage("Uploading files to S3...");
      setUploadProgress(0);
      await uploadFilesToS3(
        uploadSession.uploadUrls,
        filesToUpload,
        uploadSession.datasetId,
        uploadSession.uploadId,
        (prog, msg) => {
          setUploadProgress(prog);
          setUploadMessage(msg);
        },
      );

      // Complete upload
      setProgressMessage("Completing upload...");
      setProgress(95);
      await completeUpload(
        uploadSession.datasetId,
        uploadSession.uploadId,
        email,
        datasetName,
        dataset,
      );

      setProgress(100);
      setUploadProgress(100);
      setProgressMessage("Complete!");
      setUploadMessage("Upload successful!");
      setUploadedDatasetId(uploadSession.datasetId);
      setUploadComplete(true);
      setIsProcessing(false);
    } catch (error) {
      console.error("Upload processing error:", error);
      toast.error(
        `Failed to process dataset: ${error instanceof Error ? error.message : "Unknown error"}`,
      );
      setIsProcessing(false);
      setProgress(0);
      setProgressMessage("");
      setUploadProgress(0);
      setUploadMessage("");
    }
  };

  const baseUrl = process.env.NEXT_PUBLIC_APP_URL || "http://localhost:3000";

  return (
    <Modal isOpen={isOpen} size="2xl" onClose={onClose}>
      <ModalContent>
        <ModalHeader className="flex flex-col gap-1">
          {uploadComplete ? "Upload Successful!" : "Upload Dataset Settings"}
        </ModalHeader>
        <ModalBody>
          {uploadComplete ? (
            <div className="flex flex-col gap-4 text-center py-4">
              <div className="text-success text-6xl">✓</div>
              <h3 className="text-xl font-semibold">Upload Successful</h3>
              <div className="bg-default-100 p-4 rounded-lg">
                <p className="text-sm mb-2">Here is your dataset link:</p>
                <a
                  className="text-primary hover:underline break-all"
                  href={`${baseUrl}/viewer/${uploadedDatasetId}`}
                  rel="noopener noreferrer"
                  target="_blank"
                >
                  {baseUrl}/viewer/{uploadedDatasetId}
                </a>
              </div>
              <p className="text-sm text-default-500">
                The link has also been sent to your email:{" "}
                <strong>{email}</strong>
              </p>
            </div>
          ) : (
            <div className="flex flex-col gap-4">
              <Input
                description="Name for the saved files"
                isDisabled={isProcessing}
                label="Dataset Name"
                placeholder="Enter dataset name"
                value={datasetName}
                onValueChange={setDatasetName}
              />

              <Input
                isRequired
                description="Get notified when processing is complete"
                errorMessage={
                  email.length > 0 && !isEmailValid
                    ? "Please enter a valid email address"
                    : ""
                }
                isDisabled={isProcessing}
                isInvalid={email.length > 0 && !isEmailValid}
                label="Email"
                placeholder="Enter your email"
                type="email"
                value={email}
                onValueChange={setEmail}
              />

              {!isPreChunked && (
                <>
                  <Select
                    description="Number of genes per chunk (auto-determines based on total genes)"
                    isDisabled={isProcessing}
                    label="Chunk Size"
                    placeholder="Select chunk size"
                    selectedKeys={[chunkSize]}
                    onSelectionChange={(keys) => {
                      const value = Array.from(keys)[0] as string;

                      setChunkSize(value);
                    }}
                  >
                    <SelectItem key="auto">Auto (Recommended)</SelectItem>
                    <SelectItem key="50">50 genes/chunk</SelectItem>
                    <SelectItem key="100">100 genes/chunk</SelectItem>
                    <SelectItem key="150">150 genes/chunk</SelectItem>
                    <SelectItem key="custom">Custom</SelectItem>
                  </Select>

                  {chunkSize === "custom" && (
                    <Input
                      isDisabled={isProcessing}
                      label="Custom Chunk Size"
                      max="500"
                      min="10"
                      placeholder="Enter chunk size"
                      type="number"
                      value={customChunkSize}
                      onValueChange={setCustomChunkSize}
                    />
                  )}
                </>
              )}

              {isPreChunked && (
                <div className="bg-primary-50 dark:bg-primary-950 p-3 rounded-lg">
                  <p className="text-sm text-primary-600 dark:text-primary-400">
                    ✓ This dataset is pre-chunked and ready for upload. No processing needed.
                  </p>
                </div>
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
                    {!isPreChunked && chunkSize !== "custom" && (
                      <p>
                        <span className="font-medium">Estimated chunks:</span>{" "}
                        {Math.ceil(
                          dataset.genes.length /
                            new GeneChunkProcessor(
                              chunkSize === "auto" ? null : parseInt(chunkSize),
                            ).determineChunkSize(dataset.genes.length),
                        )}
                      </p>
                    )}
                    {isPreChunked && (
                      <p>
                        <span className="font-medium">Status:</span> Pre-chunked (ready for upload)
                      </p>
                    )}
                  </div>
                </div>
              )}

              {isProcessing && (
                <div className="space-y-4">
                  <div className="space-y-2">
                    <p className="text-xs font-medium">Processing</p>
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

                  {uploadMessage && (
                    <div className="space-y-2">
                      <p className="text-xs font-medium">Upload</p>
                      <div className="w-full bg-default-200 rounded-full h-2 overflow-hidden">
                        <div
                          className="bg-success h-full transition-all duration-300"
                          style={{ width: `${uploadProgress}%` }}
                        />
                      </div>
                      <p className="text-xs text-default-500 text-center">
                        {uploadMessage}
                      </p>
                    </div>
                  )}
                </div>
              )}
            </div>
          )}
        </ModalBody>
        <ModalFooter>
          {uploadComplete ? (
            <Button
              color="primary"
              onPress={() => {
                setUploadComplete(false);
                setUploadedDatasetId("");
                setProgress(0);
                setProgressMessage("");
                setUploadProgress(0);
                setUploadMessage("");
                onClose();
              }}
            >
              Close
            </Button>
          ) : (
            <>
              <Button
                color="danger"
                isDisabled={isProcessing}
                variant="light"
                onPress={onClose}
              >
                Cancel
              </Button>
              <Button
                color="primary"
                isDisabled={!dataset || !email || !isEmailValid}
                isLoading={isProcessing}
                onPress={handleUpload}
              >
                {isProcessing ? "Processing..." : "Process & Save"}
              </Button>
            </>
          )}
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}

/**
 * Create manifest.json
 */
async function createManifest(
  dataset: StandardizedDataset,
  datasetId: string,
  datasetName: string,
  chunks: any[],
  index: any,
  coordinates: Record<string, Blob>,
  obsMetadata: Record<string, any>,
): Promise<string> {
  // Determine spatial dimensions
  let spatialDimensions = 2;

  if (
    dataset.spatial &&
    dataset.spatial.coordinates &&
    dataset.spatial.coordinates.length > 0
  ) {
    spatialDimensions = dataset.spatial.coordinates[0].length;
  }

  // Get available embeddings
  const availableEmbeddings: string[] = [];

  if (dataset.embeddings) {
    availableEmbeddings.push(...Object.keys(dataset.embeddings));
  }

  // Count clusters
  let clusterCount = 0;

  if (dataset.clusters && dataset.clusters.length > 0) {
    clusterCount = new Set(dataset.clusters).size;
  }

  const manifest = {
    version: "2.0",
    created_at: new Date().toISOString(),
    dataset_id: datasetId,
    name: datasetName,
    type: dataset.type,
    statistics: {
      total_cells: dataset.getPointCount(),
      total_genes: dataset.genes.length,
      spatial_dimensions: spatialDimensions,
      available_embeddings: availableEmbeddings,
      cluster_count: clusterCount,
    },
    files: {
      coordinates: Object.keys(coordinates),
      expression: {
        num_chunks: chunks.length,
        chunk_size: index.chunk_size,
        total_genes: dataset.genes.length,
      },
      observations: Object.keys(obsMetadata),
      palettes: [], // TODO: Add palette support
    },
    processing: {
      chunk_strategy: "adaptive",
      compression: "gzip",
      sparse_format: true,
      created_by: "MERFISH Visualizer",
      source_file: dataset.name || "unknown",
    },
  };

  return JSON.stringify(manifest, null, 2);
}

/**
 * Save all processed files with proper S3 folder structure
 */
async function saveFilesWithStructure(
  datasetName: string,
  chunks: any[],
  index: any,
  coordinates: Record<string, Blob>,
  observationFiles: Record<string, Blob>,
  observationMetadata: Record<string, any>,
  paletteFiles: Record<string, Blob>,
  manifestJson: string,
) {
  const folderName = `${datasetName}_${Date.now()}`;
  let fileCount = 0;

  // Save manifest
  const manifestBlob = new Blob([manifestJson], { type: "application/json" });

  await saveBlob(manifestBlob, `${folderName}/manifest.json`);
  fileCount++;

  // Save gene chunks to expr/
  for (const chunk of chunks) {
    await saveBlob(chunk.data, `${folderName}/expr/${chunk.filename}`);
    fileCount++;
  }

  // Save gene index to expr/
  const indexBlob = new Blob([JSON.stringify(index, null, 2)], {
    type: "application/json",
  });

  await saveBlob(indexBlob, `${folderName}/expr/index.json`);
  fileCount++;

  // Save coordinates to coords/
  for (const [name, blob] of Object.entries(coordinates)) {
    await saveBlob(blob, `${folderName}/coords/${name}.bin.gz`);
    fileCount++;
  }

  // Save observations to obs/
  for (const [name, blob] of Object.entries(observationFiles)) {
    await saveBlob(blob, `${folderName}/obs/${name}.json.gz`);
    fileCount++;
  }

  // Save observation metadata to obs/
  const obsMetadataBlob = new Blob(
    [JSON.stringify(observationMetadata, null, 2)],
    { type: "application/json" },
  );

  await saveBlob(obsMetadataBlob, `${folderName}/obs/metadata.json`);
  fileCount++;

  // Save palettes to palettes/
  for (const [name, blob] of Object.entries(paletteFiles)) {
    await saveBlob(blob, `${folderName}/palettes/${name}.json`);
    fileCount++;
  }

  toast.info(
    `Files saved to downloads folder: ${folderName} (${fileCount} files)`,
  );
}

/**
 * Check for duplicate dataset
 */
async function checkDuplicate(
  fingerprint: string,
): Promise<{ title: string; createdAt: string } | null> {
  try {
    const response = await fetch(
      `/api/datasets/check-duplicate/${fingerprint}`,
    );

    if (response.status === 404) {
      return null; // No duplicate
    }

    if (response.ok) {
      const data = await response.json();

      return data.dataset;
    }

    throw new Error("Failed to check for duplicates");
  } catch (error) {
    console.error("Duplicate check error:", error);

    return null; // Continue on error
  }
}

/**
 * Prepare files for upload with proper structure
 */
async function prepareFilesForUpload(
  chunks: any[],
  index: any,
  coordinates: Record<string, Blob>,
  observationFiles: Record<string, Blob>,
  observationMetadata: Record<string, any>,
  paletteFiles: Record<string, Blob>,
  manifestJson: string,
): Promise<{ key: string; blob: Blob; size: number; contentType: string }[]> {
  const files: {
    key: string;
    blob: Blob;
    size: number;
    contentType: string;
  }[] = [];

  // Manifest
  const manifestBlob = new Blob([manifestJson], { type: "application/json" });

  files.push({
    key: "manifest.json",
    blob: manifestBlob,
    size: manifestBlob.size,
    contentType: "application/json",
  });

  // Gene chunks
  for (const chunk of chunks) {
    files.push({
      key: `expr/${chunk.filename}`,
      blob: chunk.data,
      size: chunk.data.size,
      contentType: "application/gzip",
    });
  }

  // Gene index
  const indexBlob = new Blob([JSON.stringify(index, null, 2)], {
    type: "application/json",
  });

  files.push({
    key: "expr/index.json",
    blob: indexBlob,
    size: indexBlob.size,
    contentType: "application/json",
  });

  // Coordinates
  for (const [name, blob] of Object.entries(coordinates)) {
    files.push({
      key: `coords/${name}.bin.gz`,
      blob: blob,
      size: blob.size,
      contentType: "application/gzip",
    });
  }

  // Observations
  for (const [name, blob] of Object.entries(observationFiles)) {
    files.push({
      key: `obs/${name}.json.gz`,
      blob: blob,
      size: blob.size,
      contentType: "application/gzip",
    });
  }

  // Observation metadata
  const obsMetadataBlob = new Blob(
    [JSON.stringify(observationMetadata, null, 2)],
    { type: "application/json" },
  );

  files.push({
    key: "obs/metadata.json",
    blob: obsMetadataBlob,
    size: obsMetadataBlob.size,
    contentType: "application/json",
  });

  // Palettes
  for (const [name, blob] of Object.entries(paletteFiles)) {
    files.push({
      key: `palettes/${name}.json`,
      blob: blob,
      size: blob.size,
      contentType: "application/json",
    });
  }

  return files;
}

/**
 * Initiate upload with backend
 */
async function initiateUpload(
  fingerprint: string,
  title: string,
  dataset: any,
  files: { key: string; size: number; contentType: string }[],
): Promise<{
  uploadId: string;
  datasetId: string;
  uploadUrls: Record<string, string>;
}> {
  const response = await fetch("/api/datasets/initiate", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      fingerprint,
      metadata: {
        title,
        numCells: dataset.getPointCount(),
        numGenes: dataset.genes.length,
        platform: dataset.type,
        description: "",
      },
      files: files.map((f) => ({
        key: f.key,
        size: f.size,
        contentType: f.contentType,
      })),
    }),
  });

  if (!response.ok) {
    const error = await response.json();

    throw new Error(error.error || "Failed to initiate upload");
  }

  const result = await response.json();

  // Map uploadUrls to presignedUrls for consistency
  return {
    uploadId: result.uploadId,
    datasetId: result.datasetId,
    uploadUrls: result.uploadUrls,
  };
}

/**
 * Upload files to S3 using presigned URLs
 */
async function uploadFilesToS3(
  uploadUrls: Record<string, string>,
  files: { key: string; blob: Blob; contentType: string }[],
  datasetId: string,
  uploadId: string,
  onProgress: (progress: number, message: string) => void,
) {
  const MAX_RETRIES = 3;
  let completed = 0;

  for (const file of files) {
    const url = uploadUrls[file.key];

    if (!url) {
      console.warn(`No presigned URL for ${file.key}, skipping`);
      continue;
    }

    let retries = 0;
    let success = false;

    while (retries < MAX_RETRIES && !success) {
      try {
        // Upload to S3
        const response = await fetch(url, {
          method: "PUT",
          body: file.blob,
          headers: {
            "Content-Type": file.contentType || "application/octet-stream",
          },
          mode: "cors", // Explicitly set CORS mode
        });

        if (!response.ok) {
          throw new Error(`Upload failed: ${response.statusText}`);
        }

        // Mark file as complete in database
        await markFileComplete(datasetId, file.key, uploadId);

        success = true;
        completed++;
        const progress = (completed / files.length) * 100;

        onProgress(progress, `Uploaded ${completed}/${files.length} files`);
      } catch (error) {
        retries++;
        if (retries >= MAX_RETRIES) {
          throw new Error(
            `Failed to upload ${file.key} after ${MAX_RETRIES} retries: ${error}`,
          );
        }
        console.warn(`Retry ${retries}/${MAX_RETRIES} for ${file.key}:`, error);
        await new Promise((resolve) => setTimeout(resolve, 1000 * retries));
      }
    }
  }
}

/**
 * Mark a file as complete in the database
 */
async function markFileComplete(
  datasetId: string,
  fileKey: string,
  uploadId: string,
): Promise<void> {
  const response = await fetch(
    `/api/datasets/${datasetId}/files/${encodeURIComponent(fileKey)}/complete`,
    {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ uploadId }),
    },
  );

  if (!response.ok) {
    console.warn(`Failed to mark ${fileKey} as complete in database`);
    // Don't throw - the file is uploaded to S3, this is just metadata
  }
}

/**
 * Complete upload
 */
async function completeUpload(
  datasetId: string,
  uploadId: string,
  email: string,
  datasetName: string,
  dataset: StandardizedDataset,
): Promise<void> {
  const response = await fetch(`/api/datasets/${datasetId}/complete`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ uploadId }),
  });

  if (!response.ok) {
    const error = await response.json();

    throw new Error(error.error || "Failed to complete upload");
  }

  // Prepare metadata for email
  const metadata = {
    numCells: dataset.getPointCount(),
    numGenes: dataset.genes.length,
    platform: dataset.type,
    clusterCount: dataset.clusters?.length || 0,
  };

  // Send email notification
  try {
    await fetch("/api/send-email", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ email, datasetId, datasetName, metadata }),
    });
  } catch (error) {
    console.error("Failed to send email notification:", error);
    // Don't throw - email notification failure shouldn't fail the upload
  }
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
