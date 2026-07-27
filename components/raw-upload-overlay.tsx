"use client";

import { Button } from "@heroui/button";
import { AlertTriangle, CheckCircle2, Loader2, XCircle } from "lucide-react";
import { useEffect, useRef } from "react";
import { createPortal } from "react-dom";

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

/** Human bytes, rolling over to GB so a multi-GB upload doesn't read "8123.4 MB". */
function formatBytes(bytes: number): string {
  const mb = bytes / 1024 ** 2;

  return mb >= 1024
    ? `${(bytes / 1024 ** 3).toFixed(2)} GB`
    : `${mb.toFixed(1)} MB`;
}

/**
 * Full-screen blocking overlay for the raw server-side upload (design §10) plus
 * live processing progress (§5.5). Status is carried by the icon; the card
 * itself is the shared `.glass-panel` so it matches every other panel.
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
  const cardRef = useRef<HTMLDivElement>(null);

  // Escape dismisses only the non-destructive states — during an active upload
  // the sole dismiss is Cancel, which aborts, and a stray Escape shouldn't
  // throw away a multi-GB transfer.
  const dismissable =
    status === "processing" || status === "done" || status === "error";

  useEffect(() => {
    cardRef.current?.focus();
  }, []);

  useEffect(() => {
    if (!dismissable) return;
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") onClose();
    };

    window.addEventListener("keydown", onKey);

    return () => window.removeEventListener("keydown", onKey);
  }, [dismissable, onClose]);

  const uploading = status === "preparing" || status === "uploading";
  const doneWarning = status === "done" && Boolean(submitWarning);

  const heading =
    status === "processing"
      ? "Uploaded — processing on the server"
      : status === "done"
        ? doneWarning
          ? "Uploaded"
          : "Dataset ready"
        : status === "error"
          ? stage
            ? "Processing failed"
            : "Upload failed"
          : "Uploading your dataset";

  const Icon = () => {
    if (uploading || status === "processing") {
      return <Loader2 className="w-6 h-6 text-blue-300 animate-spin" />;
    }
    if (status === "error") return <XCircle className="w-6 h-6 text-danger" />;
    if (doneWarning) {
      return <AlertTriangle className="w-6 h-6 text-amber-300" />;
    }

    return <CheckCircle2 className="w-6 h-6 text-emerald-300" />;
  };

  const indeterminateBar = (
    <div className="mt-6 w-full bg-white/10 rounded-full h-2.5 overflow-hidden">
      <div
        className="h-full w-1/4 rounded-full bg-blue-400"
        style={{ animation: "glassIndeterminate 1.2s ease-in-out infinite" }}
      />
    </div>
  );

  if (typeof document === "undefined") return null;

  // Portal to <body>: the upload card lives inside a transformed ancestor (the
  // crossfade layout's translate/scale), which would otherwise anchor this
  // `fixed` overlay to that ancestor instead of the viewport.
  return createPortal(
    <div
      data-ui-overlay
      className="fixed inset-0 z-[var(--z-modal)] flex items-center justify-center p-6 bg-black/60 backdrop-blur-md"
    >
      <div
        ref={cardRef}
        aria-label={heading}
        aria-modal="true"
        className="glass-panel w-full max-w-lg rounded-2xl p-8 outline-none"
        data-status={status}
        data-testid="raw-upload-card"
        role="dialog"
        tabIndex={-1}
      >
        <div className="flex items-center gap-3">
          <Icon />
          <h2
            className={`text-xl font-semibold ${status === "error" ? "text-danger" : "text-foreground"}`}
          >
            {heading}
          </h2>
        </div>

        {/* Progress + status text — announced to assistive tech as it changes. */}
        <div aria-live="polite">
          {uploading && (
            <>
              <p className="mt-1 text-base text-default-600">
                Please don&apos;t close this tab.
              </p>

              {status === "preparing" ? (
                indeterminateBar
              ) : (
                <div className="mt-6 w-full bg-white/10 rounded-full h-2.5 overflow-hidden">
                  <div
                    className="bg-blue-400 h-full transition-all duration-200 ease-out"
                    data-progress={pct.toFixed(1)}
                    data-testid="raw-upload-progress"
                    style={{ width: `${pct}%` }}
                  />
                </div>
              )}

              <div className="mt-3 flex items-center justify-between text-sm text-default-500">
                <span>
                  {status === "preparing"
                    ? `Fetching upload links… preparing ${fileCount} ${fileCount === 1 ? "file" : "files"}.`
                    : `${formatBytes(loaded)} / ${formatBytes(total)}`}
                </span>
                {status === "uploading" && (
                  <span>
                    File {Math.min(fileIndex + 1, fileCount)} of {fileCount} ·{" "}
                    {pct.toFixed(0)}%
                  </span>
                )}
              </div>
            </>
          )}

          {status === "processing" && (
            <>
              <p className="mt-1 text-base text-default-600">
                You can close this tab; we&apos;ll email you when it&apos;s
                done.
              </p>

              {typeof stagePercent === "number" ? (
                <div className="mt-6 w-full bg-white/10 rounded-full h-2.5 overflow-hidden">
                  <div
                    className="bg-blue-400 h-full transition-all duration-300 ease-out"
                    data-progress={stagePercent}
                    data-testid="raw-process-progress"
                    style={{ width: `${Math.min(100, stagePercent)}%` }}
                  />
                </div>
              ) : (
                indeterminateBar
              )}

              <div className="mt-3 flex items-center justify-between text-sm text-default-500">
                <span data-testid="raw-process-stage">
                  {stage || "Queued…"}
                </span>
                {typeof stagePercent === "number" && (
                  <span>{stagePercent}%</span>
                )}
              </div>
            </>
          )}

          {doneWarning && (
            <div className="mt-3 rounded-lg border border-amber-300/50 bg-black/25 p-3">
              <p className="text-sm font-medium text-amber-100">
                Processing hasn&apos;t started
              </p>
              <p className="mt-1 break-words text-sm text-amber-100/90">
                Your files uploaded successfully, but the processing job could
                not be submitted: {submitWarning}
              </p>
            </div>
          )}

          {status === "done" && !doneWarning && (
            <p className="mt-2 text-base text-default-600">
              Processing finished. We&apos;ve also emailed you the link.
            </p>
          )}

          {status === "error" && (
            <p className="mt-2 break-words text-base text-default-600">
              {error || "Something went wrong."}
            </p>
          )}
        </div>

        <div className="mt-6 flex justify-end gap-2">
          {status === "done" && viewerUrl && !submitWarning && (
            <Button as="a" color="success" href={viewerUrl} size="sm">
              Open dataset
            </Button>
          )}
          {uploading ? (
            <Button color="danger" size="sm" variant="flat" onPress={onCancel}>
              Cancel
            </Button>
          ) : (
            <Button size="sm" variant="flat" onPress={onClose}>
              Close
            </Button>
          )}
        </div>
      </div>
    </div>,
    document.body,
  );
}
