"use client";

import { Button } from "@heroui/button";

export type RawUploadStatus = "preparing" | "uploading" | "done" | "error";

interface RawUploadOverlayProps {
  status: RawUploadStatus;
  loaded: number;
  total: number;
  fileIndex: number;
  fileCount: number;
  error?: string;
  onCancel: () => void;
  onClose: () => void;
}

function mb(bytes: number): string {
  return (bytes / (1024 * 1024)).toFixed(1);
}

/**
 * Full-screen blocking overlay for the raw server-side upload (design §10).
 * Card tone tracks status: blue while preparing/uploading, green on success,
 * danger on error. "preparing" (fetching presigned URLs) shows an indeterminate
 * bar so the UI never sits at a stalled-looking 0 MB.
 */
export function RawUploadOverlay({
  status,
  loaded,
  total,
  fileIndex,
  fileCount,
  error,
  onCancel,
  onClose,
}: RawUploadOverlayProps) {
  const pct = total > 0 ? Math.min(100, (loaded / total) * 100) : 0;

  // Clearly colored translucent card (not the muted global .glass bg): blue
  // while preparing/uploading, green on success, danger on error.
  const tone =
    status === "done"
      ? "bg-emerald-600/40 border-emerald-400/60"
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

        {status === "done" && (
          <>
            <h2 className="text-lg font-semibold text-foreground">
              ✅ Uploaded
            </h2>
            <p className="mt-2 text-sm text-default-300">
              We&apos;ll email you when processing finishes. You can safely close
              this tab or keep working.
            </p>
            <div className="mt-6 flex justify-end">
              <Button size="sm" color="success" onPress={onClose}>
                Done
              </Button>
            </div>
          </>
        )}

        {status === "error" && (
          <>
            <h2 className="text-lg font-semibold text-danger">Upload failed</h2>
            <p className="mt-2 break-words text-sm text-default-300">
              {error || "Something went wrong during the upload."}
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
