import { gzip } from "pako";

import { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";

/**
 * Manifest structure for single molecule datasets
 */
interface SingleMoleculeManifest {
  version: string;
  created_at: string;
  dataset_id: string;
  name: string;
  type: string;
  statistics: {
    total_molecules: number;
    unique_genes: number;
    spatial_dimensions: number;
    max_raw_coordinate?: number;
  };
  genes: {
    unique_gene_names: string[];
  };
  processing: {
    compression: string;
    coordinate_format: string;
    coordinate_range: string;
    scaling_factor: number;
    created_by: string;
    source_file: string;
  };
}

/**
 * Single Molecule Dataset Processor
 * Processes single molecule datasets for download
 */
export class SingleMoleculeProcessor {
  /**
   * Create manifest for single molecule dataset
   */
  static createManifest(
    dataset: SingleMoleculeDataset,
    datasetId: string,
    datasetName: string,
  ): SingleMoleculeManifest {
    const manifest: SingleMoleculeManifest = {
      version: "1.0",
      created_at: new Date().toISOString(),
      dataset_id: datasetId,
      name: datasetName,
      type: dataset.type,
      statistics: {
        total_molecules: dataset.getMoleculeCount(),
        unique_genes: dataset.uniqueGenes.length,
        spatial_dimensions: dataset.dimensions,
        max_raw_coordinate: dataset.metadata?.maxRawCoordinate || 0,
      },
      genes: {
        unique_gene_names: [...dataset.uniqueGenes],
      },
      processing: {
        compression: "gzip",
        coordinate_format: "float32_flat_array",
        coordinate_range: "normalized_[-1,1]",
        scaling_factor: dataset.scalingFactor,
        created_by: "MERFISH Eyes - Single Molecule Viewer",
        source_file: dataset.name || "unknown",
      },
    };

    return manifest;
  }

  /**
   * Process all genes and create compressed binary files
   */
  static async processGenes(
    dataset: SingleMoleculeDataset,
    onProgress?: (progress: number, message: string) => void,
  ): Promise<Record<string, Blob>> {
    const geneFiles: Record<string, Blob> = {};
    const totalGenes = dataset.uniqueGenes.length;

    for (let i = 0; i < totalGenes; i++) {
      const gene = dataset.uniqueGenes[i];
      const progress = ((i + 1) / totalGenes) * 100;

      await onProgress?.(
        progress,
        `Processing gene ${i + 1}/${totalGenes}: ${gene}`,
      );

      // Get coordinates for this gene (already normalized in flat array format)
      const coordinates = dataset.getCoordinatesByGene(gene);

      if (coordinates && coordinates.length > 0) {
        // Convert number array to Float32Array for binary storage
        const float32Array = new Float32Array(coordinates);

        // Get the underlying buffer
        const buffer = float32Array.buffer;

        // Compress with gzip
        const compressed = gzip(new Uint8Array(buffer));
        const blob = new Blob([compressed], { type: "application/gzip" });

        // Sanitize gene name for filename
        const sanitizedName = this.sanitizeGeneName(gene);

        geneFiles[sanitizedName] = blob;
      }

      // Yield to browser every 10 genes
      if (i % 10 === 0) {
        await new Promise((resolve) => setTimeout(resolve, 0));
      }
    }

    return geneFiles;
  }

  /**
   * Sanitize gene name for use as filename
   * Replace spaces and special characters with underscores
   */
  static sanitizeGeneName(geneName: string): string {
    return geneName
      .replace(/[^a-zA-Z0-9]/g, "_")
      .replace(/_+/g, "_")
      .replace(/^_|_$/g, "");
  }

  /**
   * Create compressed manifest blob
   */
  static async createManifestBlob(
    manifest: SingleMoleculeManifest,
  ): Promise<Blob> {
    const manifestJson = JSON.stringify(manifest, null, 2);
    const compressed = gzip(manifestJson);

    return new Blob([compressed], { type: "application/gzip" });
  }

  /**
   * Save a blob to disk using the browser download API
   */
  static async saveBlob(blob: Blob, filepath: string): Promise<void> {
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

  /**
   * Process and download all files with folder structure
   */
  static async processAndDownload(
    dataset: SingleMoleculeDataset,
    datasetName: string,
    onProgress?: (progress: number, message: string) => void,
  ): Promise<void> {
    const datasetId = `${dataset.type}_${datasetName}_${Date.now()}`;
    const folderName = `${datasetName}_${Date.now()}`;

    // Step 1: Create manifest (5%)
    await onProgress?.(5, "Creating manifest...");
    const manifest = this.createManifest(dataset, datasetId, datasetName);
    const manifestBlob = await this.createManifestBlob(manifest);

    // Step 2: Process genes (5% - 90%)
    await onProgress?.(10, "Processing genes...");
    const geneFiles = await this.processGenes(
      dataset,
      (geneProgress, geneMsg) => {
        // Map gene processing progress to 10-90% range
        const overallProgress = 10 + geneProgress * 0.8;

        onProgress?.(overallProgress, geneMsg);
      },
    );

    // Step 3: Download manifest (92%)
    await onProgress?.(92, "Downloading manifest...");
    await this.saveBlob(manifestBlob, `${folderName}/manifest.json.gz`);

    // Step 4: Download gene files (92% - 100%)
    const geneFileEntries = Object.entries(geneFiles);
    const totalFiles = geneFileEntries.length;

    for (let i = 0; i < totalFiles; i++) {
      const [sanitizedName, blob] = geneFileEntries[i];
      const progress = 92 + ((i + 1) / totalFiles) * 8;

      await onProgress?.(
        progress,
        `Downloading gene file ${i + 1}/${totalFiles}...`,
      );

      await this.saveBlob(blob, `${folderName}/genes/${sanitizedName}.bin.gz`);
    }

    await onProgress?.(100, "Download complete!");
  }

  /**
   * Process and upload all files to S3
   */
  static async processAndUpload(
    dataset: SingleMoleculeDataset,
    datasetName: string,
    fingerprint: string,
    onProgress?: (progress: number, message: string) => void,
    onUploadProgress?: (progress: number, message: string) => void,
  ): Promise<{ datasetId: string; uploadId: string }> {
    // Step 1: Create manifest (5%) - use temporary ID, will be replaced by API
    await onProgress?.(5, "Creating manifest...");
    const tempDatasetId = `temp_${Date.now()}`;
    const manifest = this.createManifest(dataset, tempDatasetId, datasetName);
    const manifestBlob = await this.createManifestBlob(manifest);

    // Step 2: Process genes (5% - 45%)
    await onProgress?.(10, "Processing genes...");
    const geneFiles = await this.processGenes(
      dataset,
      (geneProgress, geneMsg) => {
        const overallProgress = 10 + geneProgress * 0.35;

        onProgress?.(overallProgress, geneMsg);
      },
    );

    // Step 3: Prepare files for upload (45%)
    await onProgress?.(45, "Preparing files for upload...");
    const filesToUpload: Array<{
      key: string;
      blob: Blob;
      size: number;
      contentType: string;
    }> = [];

    // Add manifest
    filesToUpload.push({
      key: "manifest.json.gz",
      blob: manifestBlob,
      size: manifestBlob.size,
      contentType: "application/gzip",
    });

    // Add gene files
    for (const [sanitizedName, blob] of Object.entries(geneFiles)) {
      filesToUpload.push({
        key: `genes/${sanitizedName}.bin.gz`,
        blob: blob,
        size: blob.size,
        contentType: "application/gzip",
      });
    }

    // Step 4: Initiate upload (50%)
    await onProgress?.(50, "Initiating upload...");
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
        manifest: manifest,
        files: filesToUpload.map((f) => ({
          key: f.key,
          size: f.size,
          contentType: f.contentType,
        })),
      }),
    });

    if (!initiateResponse.ok) {
      const error = await initiateResponse.json();

      throw new Error(error.error || "Failed to initiate upload");
    }

    const { datasetId, uploadId, uploadUrls } = await initiateResponse.json();

    // Step 5: Upload files to S3 (55% - 95%)
    await onProgress?.(55, "Uploading files to S3...");
    await onUploadProgress?.(0, "Starting upload...");

    let completed = 0;
    const MAX_RETRIES = 3;

    for (const file of filesToUpload) {
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
            mode: "cors",
          });

          if (!response.ok) {
            throw new Error(`Upload failed: ${response.statusText}`);
          }

          // Mark file as complete in database
          await fetch(
            `/api/single-molecule/${datasetId}/files/${encodeURIComponent(file.key)}/complete`,
            {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ uploadId }),
            },
          );

          success = true;
          completed++;
          const uploadProgress = (completed / filesToUpload.length) * 100;
          const overallProgress = 55 + uploadProgress * 0.4;

          await onProgress?.(
            overallProgress,
            `Uploading files (${completed}/${filesToUpload.length})...`,
          );
          await onUploadProgress?.(
            uploadProgress,
            `Uploaded ${completed}/${filesToUpload.length} files`,
          );
        } catch (error) {
          retries++;
          if (retries >= MAX_RETRIES) {
            throw new Error(
              `Failed to upload ${file.key} after ${MAX_RETRIES} retries: ${error}`,
            );
          }
          console.warn(
            `Retry ${retries}/${MAX_RETRIES} for ${file.key}:`,
            error,
          );
          await new Promise((resolve) => setTimeout(resolve, 1000 * retries));
        }
      }
    }

    return { datasetId, uploadId };
  }
}
