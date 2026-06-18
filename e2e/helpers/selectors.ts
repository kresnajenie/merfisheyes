import type { DatasetFormat } from "./databank";

/**
 * Central selector registry. Prefers stable ids / data-testid / roles so the
 * tests survive styling changes. Hooks added to the app live in
 * lib/utils/test-hooks.ts and the data-testid attributes referenced here.
 */

/**
 * File-input id for a given format (see components/file-upload.tsx).
 *
 * The homepage consolidated to two cards per pipeline:
 *  - single cell:    `h5ad` (single .h5ad file) + `folder` (auto-detects
 *    xenium / merscope / chunked / zarr / h5ad-in-folder from the dropped dir)
 *  - single molecule: `file` (parquet/csv, schema sniffed) + `chunked` (folder)
 * So every folder-based single-cell format goes through the one `folder` input,
 * and single-molecule parquet goes through the `file` input.
 */
export function fileInputSelector(format: DatasetFormat): string {
  switch (format) {
    case "h5ad":
      return "#sc-file-input-h5ad";
    case "xenium":
    case "merscope":
    case "chunked":
      return "#sc-file-input-folder";
    case "single-molecule":
      return "#sm-file-input-file";
  }
}

export const SELECTORS = {
  // Canvas containers (data-testid added in the scene components).
  scCanvas: '[data-testid="sc-scene-canvas"] canvas',
  smCanvas: '[data-testid="sm-scene-canvas"] canvas',
  uploadProgress: '[data-testid="upload-progress"]',
  selectedGeneBadge: '[data-testid="selected-gene-badge"]',

  // Single-cell controls.
  geneModeButton: 'button:has-text("Gene")',
  celltypeModeButton: 'button:has-text("Celltype")',
  geneSearchInput: 'input[placeholder="Search gene"]',

  // Single-molecule controls.
  smGenesButton: 'button:has-text("Genes")',
  smGeneSearchInput: 'input[placeholder="Search genes..."]',

  // Upload flow.
  uploadSaveButton: 'button:has-text("Upload & Save")',
  processSaveButton: 'button:has-text("Process & Save")',
  processUploadButton: 'button:has-text("Process & Upload")',
  uploadSuccessHeading: 'text=/Upload Success|Upload Complete/i',

  // Toasts / errors.
  toastError: ".Toastify__toast--error",
  errorBoundary: 'text="Something went wrong!"',
} as const;
