"use client";

import { useEffect, useRef, useState } from "react";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
  Button,
  Input,
} from "@heroui/react";
import { gzip } from "pako";
import { toast } from "react-toastify";

import { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import { SingleMoleculeProcessor } from "@/lib/utils/SingleMoleculeProcessor";
import { generateSingleMoleculeFingerprint } from "@/lib/utils/fingerprint";

interface SingleMoleculeUploadModalProps {
  isOpen: boolean;
  onClose: () => void;
  dataset: SingleMoleculeDataset | null;
}

export function SingleMoleculeUploadModal({
  isOpen,
  onClose,
  dataset,
}: SingleMoleculeUploadModalProps) {
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

  // Email validation
  const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
  const isEmailValid = emailRegex.test(email);

  const resetState = (targetDataset: SingleMoleculeDataset | null) => {
    setDatasetName(targetDataset?.name || "dataset");
    setEmail("");
    setIsProcessing(false);
    setProgress(0);
    setProgressMessage("");
    setUploadProgress(0);
    setUploadMessage("");
    setUploadComplete(false);
    setUploadedDatasetId("");
  };

  const previousDatasetIdRef = useRef<string | null>(null);

  useEffect(() => {
    const currentId = dataset?.id ?? null;

    if (previousDatasetIdRef.current !== currentId) {
      resetState(dataset);
      previousDatasetIdRef.current = currentId;
    }
  }, [dataset?.id]);

  useEffect(() => {
    if (isOpen) {
      setDatasetName(dataset?.name || "dataset");
    }
  }, [dataset, isOpen]);

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
      // Generate fingerprint
      const fingerprint = await generateSingleMoleculeFingerprint(dataset);

      console.log("Single molecule dataset fingerprint:", fingerprint);

      // Check for duplicates
      setProgressMessage("Checking for duplicates...");
      setProgress(2);

      const duplicateResponse = await fetch(
        `/api/single-molecule/check-duplicate/${fingerprint}`,
      );

      if (duplicateResponse.ok) {
        const duplicateData = await duplicateResponse.json();

        if (duplicateData.exists) {
          const confirmed = window.confirm(
            `A dataset with this content already exists.\n\nDataset: ${duplicateData.dataset.title}\nUploaded: ${new Date(duplicateData.dataset.createdAt).toLocaleString()}\nMolecules: ${duplicateData.dataset.numMolecules.toLocaleString()}\n\nDo you want to upload anyway?`,
          );

          if (!confirmed) {
            setIsProcessing(false);
            setProgress(0);
            setProgressMessage("");

            return;
          }
        }
      }

      // Check if this is a pre-chunked dataset
      const isPreChunked = (dataset as any).isPreChunked === true;
      let datasetId: string;
      let uploadId: string;

      if (isPreChunked) {
        // Pre-chunked dataset - upload directly without processing
        setProgressMessage("Uploading pre-chunked dataset...");
        setProgress(5);

        const fileMap = (dataset as any).fileMap as Map<string, File>;

        if (!fileMap) {
          throw new Error("Pre-chunked dataset missing file map");
        }

        const adapter = (dataset as any).adapter;
        const manifest = adapter.getManifest();

        // Normalize manifest format: ensure manifest.json.gz exists for S3 consistency
        if (fileMap.has("manifest.json") && !fileMap.has("manifest.json.gz")) {
          const plainManifestFile = fileMap.get("manifest.json")!;
          const manifestText = await plainManifestFile.text();
          const compressed = gzip(manifestText);
          const gzippedBlob = new Blob([compressed], { type: "application/gzip" });
          const gzippedFile = new File([gzippedBlob], "manifest.json.gz", {
            type: "application/gzip",
          });

          fileMap.delete("manifest.json");
          fileMap.set("manifest.json.gz", gzippedFile);
        }

        // Build files list from fileMap for the API
        const filesList = Array.from(fileMap.entries()).map(([key, file]) => ({
          key,
          size: file.size,
          contentType: key.endsWith(".gz")
            ? "application/gzip"
            : "application/octet-stream",
        }));

        // Initiate upload with correct request body
        const initiateResponse = await fetch("/api/single-molecule/initiate", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            fingerprint,
            metadata: {
              title: datasetName,
              numMolecules: dataset.getMoleculeCount(),
              numGenes: dataset.uniqueGenes.length,
              platform: dataset.type,
              description: "",
            },
            manifest,
            files: filesList,
          }),
        });

        if (!initiateResponse.ok) {
          const error = await initiateResponse.json();

          throw new Error(error.error || "Failed to initiate upload");
        }

        const { datasetId: newDatasetId, uploadId: newUploadId, uploadUrls } =
          await initiateResponse.json();

        datasetId = newDatasetId;
        uploadId = newUploadId;

        // Upload files directly to S3
        const files = Array.from(fileMap.entries());
        const totalFiles = files.length;

        setProgressMessage(`Uploading ${totalFiles} files...`);

        for (let i = 0; i < files.length; i++) {
          const [key, file] = files[i];
          const presignedUrl = uploadUrls[key];

          if (!presignedUrl) {
            console.warn(`No presigned URL for ${key}, skipping...`);
            continue;
          }

          // Upload to S3
          const uploadResponse = await fetch(presignedUrl, {
            method: "PUT",
            body: file,
            headers: {
              "Content-Type": key.endsWith(".gz")
                ? "application/gzip"
                : "application/octet-stream",
            },
          });

          if (!uploadResponse.ok) {
            throw new Error(`Failed to upload ${key} to S3`);
          }

          // Mark file as complete
          await fetch(`/api/single-molecule/${datasetId}/files/${encodeURIComponent(key)}/complete`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ uploadId }),
          });

          const progress = 5 + ((i + 1) / totalFiles) * 85;

          setProgress(progress);
          setProgressMessage(
            `Uploaded ${i + 1}/${totalFiles} files...`,
          );
        }

        setProgress(90);
      } else {
        // Regular dataset - process and upload
        setProgressMessage("Processing dataset...");
        setProgress(5);

        const result = await SingleMoleculeProcessor.processAndUpload(
          dataset,
          datasetName,
          fingerprint,
          (prog, msg) => {
            setProgress(prog);
            setProgressMessage(msg);
          },
          (prog, msg) => {
            setUploadProgress(prog);
            setUploadMessage(msg);
          },
        );

        datasetId = result.datasetId;
        uploadId = result.uploadId;
      }

      // Complete upload
      setProgressMessage("Completing upload...");
      setProgress(95);

      const completeResponse = await fetch(
        `/api/single-molecule/${datasetId}/complete`,
        {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ uploadId, email }),
        },
      );

      if (!completeResponse.ok) {
        const error = await completeResponse.json();

        throw new Error(error.error || "Failed to complete upload");
      }

      setProgress(100);
      setUploadProgress(100);
      setProgressMessage("Complete!");
      setUploadMessage("Upload successful!");
      setUploadedDatasetId(datasetId);
      setUploadComplete(true);
      setIsProcessing(false);

      toast.success("Dataset uploaded successfully!");
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

  return (
    <Modal isOpen={isOpen} size="2xl" onClose={onClose}>
      <ModalContent>
        <ModalHeader className="flex flex-col gap-1">
          Upload Single Molecule Dataset
        </ModalHeader>
        <ModalBody>
          <div className="flex flex-col gap-4">
            <Input
              description="Name for the saved files"
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
              isInvalid={email.length > 0 && !isEmailValid}
              label="Email"
              placeholder="Enter your email"
              type="email"
              value={email}
              onValueChange={setEmail}
            />

            {dataset && (
              <div className="bg-default-100 p-4 rounded-lg">
                <h4 className="text-sm font-semibold mb-2">
                  Dataset Information
                </h4>
                <div className="text-xs space-y-1">
                  <p>
                    <span className="font-medium">Dataset Name:</span>{" "}
                    {dataset.name}
                  </p>
                  <p>
                    <span className="font-medium">Type:</span> {dataset.type}
                  </p>
                  <p>
                    <span className="font-medium">Total Molecules:</span>{" "}
                    {dataset.getMoleculeCount().toLocaleString()}
                  </p>
                  <p>
                    <span className="font-medium">Unique Genes:</span>{" "}
                    {dataset.uniqueGenes.length.toLocaleString()}
                  </p>
                  <p>
                    <span className="font-medium">Dimensions:</span>{" "}
                    {dataset.dimensions}D
                  </p>
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

                {uploadProgress > 0 && (
                  <div className="space-y-2">
                    <p className="text-xs font-medium">Uploading to S3</p>
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

            {uploadComplete && (
              <div className="bg-success-50 border border-success-200 p-4 rounded-lg">
                <h4 className="text-sm font-semibold text-success-700 mb-2">
                  Upload Complete!
                </h4>
                <p className="text-xs text-success-600 mb-3">
                  Your single molecule dataset has been uploaded successfully
                  and is ready to view.
                </p>
                <a
                  className="inline-block text-xs font-medium text-success-700 hover:text-success-800 underline"
                  href={`/sm-viewer/${uploadedDatasetId}`}
                  rel="noopener noreferrer"
                  target="_blank"
                >
                  Open in Single Molecule Viewer →
                </a>
              </div>
            )}
          </div>
        </ModalBody>
        <ModalFooter>
          {uploadComplete ? (
            <Button color="primary" onPress={onClose}>
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
                {isProcessing ? "Processing..." : "Process & Upload"}
              </Button>
            </>
          )}
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}
