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

import { parseS3Url, getS3FileUrl, testManifestAccess } from "@/lib/utils/s3-url-parser";
import { ungzip } from "@/lib/utils/gzip";

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
      setValidationStep("Validating S3 URL...");

      // Parse and normalize S3 URL
      const parsed = parseS3Url(s3Url);
      console.log("Parsed S3 URL:", parsed);

      setValidationStep("Testing manifest access...");

      // Test if manifest exists
      const manifestUrl = getS3FileUrl(parsed.baseUrl, "manifest.json.gz");
      const manifestExists = await testManifestAccess(parsed.baseUrl);

      if (!manifestExists) {
        throw new Error(
          `Cannot access manifest at ${manifestUrl}. Please ensure:\n` +
            `1. The folder URL is correct\n` +
            `2. manifest.json.gz exists in the folder\n` +
            `3. The S3 bucket/folder has public read access\n` +
            `4. CORS is configured to allow requests from ${window.location.origin}`
        );
      }

      setValidationStep("Downloading and parsing manifest...");

      // Download and parse manifest to detect dataset type
      const response = await fetch(manifestUrl);
      if (!response.ok) {
        throw new Error(`Failed to download manifest: ${response.statusText}`);
      }

      const manifestCompressed = await response.arrayBuffer();
      const manifestJson = ungzip(new Uint8Array(manifestCompressed), {
        to: "string",
      });
      const manifest = JSON.parse(manifestJson);

      console.log("Manifest loaded:", manifest);

      // Auto-detect dataset type if not provided
      let detectedType: "single_cell" | "single_molecule";

      if (datasetType) {
        detectedType = datasetType;
      } else {
        // Detect based on manifest structure
        if (manifest.type === "single_molecule" || manifest.genes?.unique_gene_names) {
          detectedType = "single_molecule";
        } else if (manifest.type === "single_cell" || manifest.cells || manifest.spatial_dimensions) {
          detectedType = "single_cell";
        } else {
          throw new Error(
            "Could not auto-detect dataset type from manifest. Please specify dataset type explicitly."
          );
        }
      }

      console.log("Detected dataset type:", detectedType);

      setValidationStep("Loading dataset...");

      // Encode base URL to pass in URL parameter
      const encodedBaseUrl = encodeURIComponent(parsed.baseUrl);

      // Navigate to appropriate viewer with custom S3 URL
      if (detectedType === "single_molecule") {
        router.push(`/sm-viewer/from-s3?url=${encodedBaseUrl}`);
      } else {
        router.push(`/viewer/from-s3?url=${encodedBaseUrl}`);
      }

      // Close modal
      onClose();
    } catch (err) {
      console.error("Error loading from S3:", err);
      setError(err instanceof Error ? err.message : "Failed to load dataset from S3");
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
    <Modal
      isOpen={isOpen}
      size="2xl"
      onClose={handleClose}
    >
      <ModalContent>
        <ModalHeader className="flex flex-col gap-1">
          Load Dataset from S3
          {datasetType && (
            <span className="text-sm font-normal text-default-500">
              {datasetType === "single_cell" ? "Single Cell" : "Single Molecule"} Dataset
            </span>
          )}
        </ModalHeader>
        <ModalBody>
          <div className="flex flex-col gap-4">
            <p className="text-sm text-default-600">
              Enter the S3 URL to your dataset folder. The folder should contain a{" "}
              <code className="text-xs bg-default-100 px-1 py-0.5 rounded">
                manifest.json.gz
              </code>{" "}
              file and associated data files.
            </p>

            <Input
              label="S3 Folder URL"
              placeholder="https://my-bucket.s3.us-east-1.amazonaws.com/datasets/my-data"
              value={s3Url}
              variant="bordered"
              isDisabled={isLoading}
              isInvalid={!!error}
              errorMessage={error}
              onValueChange={setS3Url}
              description="Can be a folder URL or direct manifest URL"
            />

            {isLoading && (
              <div className="flex items-center gap-3 p-4 bg-default-50 rounded-lg">
                <Spinner size="sm" />
                <p className="text-sm text-default-600">{validationStep}</p>
              </div>
            )}

            <div className="text-xs text-default-500 space-y-2">
              <p className="font-semibold">Requirements:</p>
              <ul className="list-disc list-inside space-y-1 ml-2">
                <li>S3 bucket/folder must have <strong>public read access</strong></li>
                <li>
                  CORS must be configured to allow requests from{" "}
                  <code className="bg-default-100 px-1 py-0.5 rounded">
                    {typeof window !== "undefined" ? window.location.origin : "this domain"}
                  </code>
                </li>
                <li>Dataset must be processed using the MERFISHeyes Python scripts</li>
              </ul>
            </div>

            {!datasetType && (
              <div className="text-xs text-default-500 p-3 bg-primary-50 rounded-lg">
                <p className="font-semibold text-primary">Auto-detection enabled</p>
                <p className="mt-1">
                  Dataset type (Single Cell vs Single Molecule) will be automatically detected
                  from the manifest file.
                </p>
              </div>
            )}
          </div>
        </ModalBody>
        <ModalFooter>
          <Button
            color="danger"
            variant="light"
            onPress={handleClose}
            isDisabled={isLoading}
          >
            Cancel
          </Button>
          <Button
            color="primary"
            onPress={handleLoad}
            isDisabled={!s3Url.trim() || isLoading}
            isLoading={isLoading}
          >
            Load Dataset
          </Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}
