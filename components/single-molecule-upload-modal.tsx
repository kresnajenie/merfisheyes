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
} from "@heroui/react";
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

      // Process and upload to S3
      setProgressMessage("Processing dataset...");
      setProgress(5);

      const { datasetId, uploadId } =
        await SingleMoleculeProcessor.processAndUpload(
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
