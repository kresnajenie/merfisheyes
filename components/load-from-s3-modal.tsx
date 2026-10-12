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
import { ungzip } from "pako";

import { resolveZarrUrl } from "@/lib/spatialdata/resolveUrl";
import { parseS3Url, testManifestAccess } from "@/lib/utils/s3-url-parser";

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

      // Test if manifest exists (tries .gz first, then .json)
      const manifestTest = await testManifestAccess(parsed.baseUrl);

      if (!manifestTest.exists) {
        // Not a merfisheyes folder: a SpatialData store, an OME-Zarr image,
        // or a folder holding a store. Each is read as it is by its viewer.
        setValidationStep("Looking for a SpatialData store or an image...");
        const target = await resolveZarrUrl(parsed.baseUrl);

        if (!target) {
          throw new Error(
            `Nothing to open at that URL. Please ensure:\n` +
              `1. The URL is a merfisheyes folder (manifest.json or manifest.json.gz), a SpatialData store or an OME-Zarr image (zarr.json), or a folder holding a store\n` +
              `2. A folder names its store in mapping.json ("cells"), or the store is named after the folder (sample/sample.zarr)\n` +
              `3. The S3 bucket/folder has public read access\n` +
              `4. CORS is configured to allow requests from ${window.location.origin}`,
          );
        }
        // A store's link keeps the URL as entered: the viewer finds the
        // store in a folder the same way.
        router.push(
          target.kind === "image"
            ? `/image-viewer/from-s3?url=${encodeURIComponent(target.url)}`
            : `/spatialdata-viewer/from-s3?url=${encodeURIComponent(parsed.baseUrl)}`,
        );
        onClose();

        return;
      }

      setValidationStep("Downloading and parsing manifest...");

      // Download and parse manifest to detect dataset type
      const response = await fetch(manifestTest.url);

      if (!response.ok) {
        throw new Error(`Failed to download manifest: ${response.statusText}`);
      }

      let manifestJson: string;

      // Check if it's gzipped based on URL
      if (manifestTest.url.endsWith(".gz")) {
        const manifestCompressed = await response.arrayBuffer();

        manifestJson = ungzip(new Uint8Array(manifestCompressed), {
          to: "string",
        });
      } else {
        manifestJson = await response.text();
      }

      const manifest = JSON.parse(manifestJson);

      console.log("Manifest loaded:", manifest);

      // Always detect the type from the manifest's own fields — that is the
      // reliable signal. The `datasetType` prop is only the homepage's mode
      // toggle, and it must NOT override a manifest that clearly says
      // otherwise: pasting a single-molecule URL while the toggle sits on
      // "single cell" used to force the cell viewer, which then crashed on a
      // molecule manifest. We only fall back to the toggle when the manifest
      // is genuinely ambiguous.
      //
      // Reliable, mutually-exclusive fields (see scripts/process_*.py):
      //   SM manifest → genes.unique_gene_names, statistics.total_molecules
      //   SC manifest → statistics.total_cells, files.expression_chunks
      // Do NOT key off `manifest.type`: both pipelines write the *format*
      // there ("h5ad" | "xenium" | "merscope"), so "merscope"/"xenium" appear
      // on both kinds and can't discriminate.
      const looksSingleMolecule =
        Boolean(manifest.genes?.unique_gene_names) ||
        manifest.statistics?.total_molecules != null;
      const looksSingleCell =
        manifest.statistics?.total_cells != null ||
        manifest.files?.expression_chunks != null;

      let detectedType: "single_cell" | "single_molecule";

      if (looksSingleMolecule && !looksSingleCell) {
        detectedType = "single_molecule";
      } else if (looksSingleCell && !looksSingleMolecule) {
        detectedType = "single_cell";
      } else if (datasetType) {
        // Ambiguous or self-contradictory manifest — defer to the mode toggle.
        detectedType = datasetType;
      } else {
        throw new Error(
          "Could not determine from the manifest whether this is a single-cell " +
            "or single-molecule dataset.",
        );
      }

      console.log("Detected dataset type:", detectedType);

      // The toggle said one thing but the manifest is another — we route by the
      // manifest (correct), but tell the user so the different viewer isn't a
      // surprise.
      if (datasetType && detectedType !== datasetType) {
        setValidationStep(
          `This is a ${
            detectedType === "single_molecule"
              ? "single-molecule"
              : "single-cell"
          } dataset — opening the matching viewer…`,
        );
      }

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
