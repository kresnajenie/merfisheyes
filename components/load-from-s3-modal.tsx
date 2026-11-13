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
  Spinner,
} from "@heroui/react";
import { useRouter } from "next/navigation";

import { viewerPathForUrl } from "@/lib/viewer/viewer-for-url";

interface LoadFromS3ModalProps {
  isOpen: boolean;
  onClose: () => void;
  /** Dataset type hint - if provided, skips auto-detection */
  datasetType?: "single_cell" | "single_molecule";
}

export function LoadFromS3Modal({
  isOpen,
  onClose,
  datasetType,
}: LoadFromS3ModalProps) {
  const router = useRouter();
  const [s3Url, setS3Url] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [validationStep, setValidationStep] = useState<string>("");

  const handleLoad = async () => {
    try {
      setIsLoading(true);
      setError(null);
      const path = await viewerPathForUrl(s3Url, {
        datasetType,
        onStep: setValidationStep,
      });

      router.push(path);
      onClose();
    } catch (err) {
      setError(
        err instanceof Error ? err.message : "Failed to load dataset from S3",
      );
      setIsLoading(false);
    }
  };

  const handleClose = () => {
    if (!isLoading) {
      setS3Url("");
      setError(null);
      setValidationStep("");
      onClose();
    }
  };

  return (
    <Modal isOpen={isOpen} size="2xl" onClose={handleClose}>
      <ModalContent>
        <ModalHeader className="flex flex-col gap-1">
          Load Dataset from S3
          {datasetType && (
            <span className="text-sm font-normal text-default-500">
              {datasetType === "single_cell"
                ? "Single Cell"
                : "Single Molecule"}{" "}
              Dataset
            </span>
          )}
        </ModalHeader>
        <ModalBody>
          <div className="flex flex-col gap-4">
            <p className="text-sm text-default-600">
              Enter the S3 URL to your dataset folder: a merfisheyes folder
              (with a{" "}
              <code className="text-xs bg-default-100 px-1 py-0.5 rounded">
                manifest.json.gz
              </code>
              ), a SpatialData{" "}
              <code className="text-xs bg-default-100 px-1 py-0.5 rounded">
                .zarr
              </code>{" "}
              store, an OME-Zarr image, or a folder holding a store with its
              molecules beside it — those open overlaid.
            </p>

            <Input
              label="Dataset Folder URL"
              placeholder="https://my-bucket.s3.us-east-1.amazonaws.com/datasets/my-data"
              value={s3Url}
              variant="bordered"
              onValueChange={setS3Url}
              description="S3 URL or any public HTTP URL of a dataset folder, store or image"
            />

            {isLoading && (
              <div className="flex items-center gap-3 p-4 bg-default-50 rounded-lg">
                <Spinner size="sm" />
                <p className="text-sm text-default-600">{validationStep}</p>
              </div>
            )}

            {error && (
              <p className="whitespace-pre-line rounded-lg bg-danger-50 p-3 text-sm text-danger">
                {error}
              </p>
            )}

            <div className="text-xs text-default-500 space-y-3">
              <div>
                <p className="font-semibold mb-2">Requirements:</p>
                <ul className="list-disc list-inside space-y-1 ml-2">
                  <li>
                    S3 bucket/folder must have{" "}
                    <strong>public read access</strong>
                  </li>
                  <li>CORS must be configured (see below)</li>
                  <li>
                    A merfisheyes folder (from the MERFISHeyes Python scripts),
                    a SpatialData store or an OME-Zarr image, on Zarr v3
                  </li>
                </ul>
              </div>

              <div className="bg-default-100 p-3 rounded-lg">
                <p className="font-semibold mb-2">CORS Configuration:</p>
                <p className="mb-2">
                  Go to:{" "}
                  <strong>
                    S3 Bucket → Permissions → Cross-origin resource sharing
                    (CORS)
                  </strong>
                </p>
                <p className="mb-1">Paste this configuration:</p>
                <pre className="bg-default-50 p-2 rounded text-[10px] overflow-x-auto">
                  {`[
  {
    "AllowedHeaders": ["*"],
    "AllowedMethods": ["GET", "HEAD"],
    "AllowedOrigins": [
      "https://www.merfisheyes.com",
      "https://merfisheyes.com"
    ],
    "ExposeHeaders": [],
    "MaxAgeSeconds": 3600
  }
]`}
                </pre>
              </div>
            </div>

            {!datasetType && (
              <div className="text-xs text-default-500 p-3 bg-primary-50 rounded-lg">
                <p className="font-semibold text-primary">
                  Auto-detection enabled
                </p>
                <p className="mt-1">
                  Dataset type (Single Cell vs Single Molecule) will be
                  automatically detected from the manifest file.
                </p>
              </div>
            )}
          </div>
        </ModalBody>
        <ModalFooter>
          <Button
            color="danger"
            isDisabled={isLoading}
            variant="light"
            onPress={handleClose}
          >
            Cancel
          </Button>
          <Button
            color="primary"
            isDisabled={!s3Url.trim() || isLoading}
            isLoading={isLoading}
            onPress={handleLoad}
          >
            Load Dataset
          </Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}
