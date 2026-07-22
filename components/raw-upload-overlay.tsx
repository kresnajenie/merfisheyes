"use client";

import { Button } from "@heroui/button";

export type RawUploadStatus =
  | "preparing"
  | "uploading"
  | "processing"
  | "done"
  | "error";

interface RawUploadOverlayProps {
  status: RawUploadStatus;
  loaded: number;
  total: number;
  fileIndex: number;
  fileCount: number;
  error?: string;
  /** Live processing stage from /status, e.g. "Expression chunks". */
  stage?: string;
  stagePercent?: number;
  /** Set when the upload succeeded but no processing job was submitted. */
  submitWarning?: string;
  viewerUrl?: string | null;
  onCancel: () => void;
  onClose: () => void;
}

function mb(bytes: number): string {
  return (bytes / (1024 * 1024)).toFixed(1);
}

/**
 * Full-screen blocking overlay for the raw server-side upload (design §10) plus
 * live processing progress (§5.5). Card tone tracks status: blue while
 * preparing/uploading/processing, green on success, danger on error.
 */
export function RawUploadOverlay({
  status,
  loaded,
  total,
  fileIndex,
  fileCount,
  error,
  stage,
  stagePercent,
  submitWarning,
  viewerUrl,
  onCancel,
  onClose,
}: RawUploadOverlayProps) {
  const pct = total > 0 ? Math.min(100, (loaded / total) * 100) : 0;

  // Green means "done and good". An upload whose processing never started is
  // NOT that, so it gets amber rather than a success-coloured card.
  const tone =
    status === "done"
      ? submitWarning
        ? "bg-amber-600/35 border-amber-400/60"
        : "bg-emerald-600/40 border-emerald-400/60"
      : status === "error"
        ? "bg-danger-600/35 border-danger-400/60"
        : "bg-blue-600/40 border-blue-400/60";

  return (
    <div
      data-ui-overlay
      className="fixed inset-0 z-[10001] flex items-center justify-center p-6 bg-black/60 backdrop-blur-md"
    >
      <style>{`
        @keyframes rawIndeterminate {
          0% { transform: translateX(-100%); }
          100% { transform: translateX(400%); }
        }
      `}</style>

      <div
        data-testid="raw-upload-card"
        data-status={status}
        className={`w-full max-w-lg rounded-2xl border p-8 shadow-2xl backdrop-blur-2xl transition-colors duration-500 ${tone}`}
      >
        {status === "preparing" && (
          <>
            <h2 className="text-lg font-semibold text-foreground">
              Uploading your dataset
            </h2>
            <p className="mt-1 text-sm text-default-400">
              Please don&apos;t close this tab.
            </p>

            <div className="mt-6 w-full bg-default-200/30 rounded-full h-3 overflow-hidden">
              <div
                className="h-full w-1/4 rounded-full bg-blue-400"
                style={{ animation: "rawIndeterminate 1.2s ease-in-out infinite" }}
              />
            </div>

            <p className="mt-3 text-xs text-default-400">
              Fetching upload links… preparing {fileCount}{" "}
              {fileCount === 1 ? "file" : "files"}.
            </p>

            <div className="mt-6 flex justify-end">
              <Button size="sm" variant="flat" color="danger" onPress={onCancel}>
                Cancel
              </Button>
            </div>
          </>
        )}

        {status === "uploading" && (
          <>
            <h2 className="text-lg font-semibold text-foreground">
              Uploading your dataset
            </h2>
            <p className="mt-1 text-sm text-default-400">
              Please don&apos;t close this tab.
            </p>

            <div className="mt-6 w-full bg-default-200/30 rounded-full h-3 overflow-hidden">
              <div
                data-testid="raw-upload-progress"
                data-progress={pct.toFixed(1)}
                className="bg-blue-400 h-full transition-all duration-200 ease-out"
                style={{ width: `${pct}%` }}
              />
            </div>

            <div className="mt-3 flex items-center justify-between text-xs text-default-400">
              <span>
                {mb(loaded)} MB / {mb(total)} MB
              </span>
              <span>
                File {Math.min(fileIndex + 1, fileCount)} of {fileCount} ·{" "}
                {pct.toFixed(0)}%
              </span>
            </div>

            <div className="mt-6 flex justify-end">
              <Button size="sm" variant="flat" color="danger" onPress={onCancel}>
                Cancel
              </Button>
            </div>
          </>
        )}

        {status === "processing" && (
          <>
            <h2 className="text-lg font-semibold text-foreground">
              ✅ Uploaded — processing on the server
            </h2>
            <p className="mt-1 text-sm text-default-300">
              You can close this tab; we&apos;ll email you when it&apos;s done.
            </p>

            <div className="mt-6 w-full bg-default-200/30 rounded-full h-3 overflow-hidden">
              {typeof stagePercent === "number" ? (
                <div
                  data-testid="raw-process-progress"
                  data-progress={stagePercent}
                  className="bg-blue-400 h-full transition-all duration-300 ease-out"
                  style={{ width: `${Math.min(100, stagePercent)}%` }}
                />
              ) : (
                <div
                  className="h-full w-1/4 rounded-full bg-blue-400"
                  style={{ animation: "rawIndeterminate 1.2s ease-in-out infinite" }}
                />
              )}
            </div>

            <div className="mt-3 flex items-center justify-between text-xs text-default-300">
              <span data-testid="raw-process-stage">{stage || "Queued…"}</span>
              {typeof stagePercent === "number" && <span>{stagePercent}%</span>}
            </div>

            <div className="mt-6 flex justify-end">
              <Button size="sm" variant="flat" onPress={onClose}>
                Close
              </Button>
            </div>
          </>
        )}

        {status === "done" && (
          <>
            <h2 className="text-lg font-semibold text-foreground">
              ✅ {submitWarning ? "Uploaded" : "Dataset ready"}
            </h2>

            {submitWarning ? (
              <div className="mt-3 rounded-lg border border-amber-300/50 bg-black/25 p-3">
                <p className="text-sm font-medium text-amber-100">
                  Processing hasn&apos;t started
                </p>
                <p className="mt-1 break-words text-xs text-amber-100/80">
                  Your files uploaded successfully, but the processing job could
                  not be submitted: {submitWarning}
                </p>
              </div>
            ) : (
              <p className="mt-2 text-sm text-default-300">
                Processing finished. We&apos;ve also emailed you the link.
              </p>
            )}

            <div className="mt-6 flex justify-end gap-2">
              {viewerUrl && !submitWarning && (
                <Button as="a" href={viewerUrl} size="sm" color="success">
                  Open dataset
                </Button>
              )}
              <Button size="sm" variant="flat" onPress={onClose}>
                Close
              </Button>
            </div>
          </>
        )}

        {status === "error" && (
          <>
            <h2 className="text-lg font-semibold text-danger">
              {stage ? "Processing failed" : "Upload failed"}
            </h2>
            <p className="mt-2 break-words text-sm text-default-300">
              {error || "Something went wrong."}
            </p>
            <div className="mt-6 flex justify-end">
              <Button size="sm" variant="flat" onPress={onClose}>
                Close
              </Button>
            </div>
          </>
        )}
      </div>
    </div>
  );
}
