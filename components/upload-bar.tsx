"use client";

import type { UploadJob } from "@/lib/stores/uploadStore";

import { AlertTriangle, CheckCircle2, Loader2, XCircle } from "lucide-react";
import Link from "next/link";
import { useEffect, useRef, useState } from "react";

import { useUploadStore } from "@/lib/stores/uploadStore";

/** Human bytes, rolling to GB so a multi-GB upload doesn't read "8123.4 MB". */
function formatBytes(bytes: number): string {
  const mb = bytes / 1024 ** 2;

  return mb >= 1024
    ? `${(bytes / 1024 ** 3).toFixed(2)} GB`
    : `${mb.toFixed(1)} MB`;
}

function formatRate(bps: number): string {
  if (bps <= 0) return "—";
  const mbps = bps / 1024 ** 2;

  return mbps >= 1
    ? `${mbps.toFixed(1)} MB/s`
    : `${(bps / 1024).toFixed(0)} KB/s`;
}

function formatEta(sec: number | null): string {
  if (sec == null || !isFinite(sec) || sec < 0) return "Calculating…";
  if (sec < 60) {
    const s = Math.ceil(sec);

    return `${s} second${s === 1 ? "" : "s"}`;
  }
  const m = Math.round(sec / 60);

  if (m < 60) return `${m} minute${m === 1 ? "" : "s"}`;
  const h = Math.floor(m / 60);
  const mm = m % 60;

  return `${h} hr ${mm} min`;
}

/** One colored banner row for a single upload job. */
function UploadJobRow({ job }: { job: UploadJob }) {
  const cancel = useUploadStore((s) => s.cancel);
  const dismiss = useUploadStore((s) => s.dismiss);

  const {
    id,
    status,
    loaded,
    total,
    fileCount,
    title,
    error,
    submitWarning,
    rateBps,
    etaSeconds,
    stage,
    percent,
    viewerUrl,
  } = job;

  const uploading = status === "preparing" || status === "uploading";
  const processingPhase = status === "queued" || status === "processing";
  const complete = status === "complete";
  const busy = uploading || processingPhase;

  const pct = total > 0 ? Math.min(100, (loaded / total) * 100) : 0;
  const remaining = Math.max(0, total - loaded);
  const remainingPct = Math.max(0, 100 - pct);
  const doneWarning = status === "done" && Boolean(submitWarning);
  const preparing = status === "preparing";
  const hasPercent = processingPhase && percent != null;

  const bg =
    status === "error"
      ? "bg-danger"
      : doneWarning
        ? "bg-amber-600"
        : complete
          ? "bg-emerald-700"
          : "bg-[#0b6bcb]";

  const heading =
    status === "error"
      ? "Upload failed"
      : complete
        ? "Ready to view"
        : status === "queued"
          ? "Queued"
          : status === "processing"
            ? "Processing"
            : status === "done"
              ? "Uploaded"
              : preparing
                ? "Preparing upload"
                : "Uploading";

  const Icon = () => {
    if (busy) return <Loader2 className="h-5 w-5 animate-spin" />;
    if (status === "error") return <XCircle className="h-5 w-5" />;
    if (doneWarning) return <AlertTriangle className="h-5 w-5" />;

    return <CheckCircle2 className="h-5 w-5" />;
  };

  return (
    <div className={`w-full text-white ${bg}`} data-status={status}>
      <div className="mx-auto max-w-6xl px-6 py-4">
        <div className="flex items-center gap-4">
          <Icon />
          <span className="min-w-0 max-w-[14rem] shrink-0 truncate text-lg font-semibold">
            {title}
          </span>
          <span className="shrink-0 text-sm text-white/80">{heading}</span>

          {/* Progress bar. Indeterminate while preparing or while the server
              works without a reported percent; determinate otherwise. */}
          <div className="h-2 flex-1 overflow-hidden rounded-full bg-white/25">
            {preparing || (processingPhase && !hasPercent) ? (
              <div
                className="h-full w-1/4 rounded-full bg-white"
                style={{
                  animation: "glassIndeterminate 1.2s ease-in-out infinite",
                }}
              />
            ) : (
              <div
                className="h-full rounded-full bg-white transition-all duration-300"
                data-progress={(complete
                  ? 100
                  : hasPercent
                    ? (percent ?? 0)
                    : pct
                ).toFixed(1)}
                data-testid="upload-bar-progress"
                style={{
                  width: `${
                    status === "error" || complete
                      ? 100
                      : hasPercent
                        ? (percent ?? 0)
                        : pct
                  }%`,
                }}
              />
            )}
          </div>

          {(status === "uploading" || complete || hasPercent) && (
            <span className="w-12 shrink-0 text-right text-sm tabular-nums">
              {(complete ? 100 : hasPercent ? (percent ?? 0) : pct).toFixed(0)}%
            </span>
          )}

          {complete && viewerUrl && (
            <Link
              className="shrink-0 rounded-md bg-white px-4 py-1.5 text-sm font-semibold text-emerald-700 transition-colors hover:bg-white/90"
              href={viewerUrl}
            >
              View dataset
            </Link>
          )}

          <button
            className="shrink-0 rounded-md border border-white/70 px-4 py-1.5 text-sm font-medium transition-colors hover:bg-white/15"
            type="button"
            onClick={() => (uploading ? cancel(id) : dismiss(id))}
          >
            {uploading ? "Cancel" : "Close"}
          </button>
        </div>

        {/* Detail lines (AWS-style) while transferring. */}
        {status === "uploading" && (
          <div className="mt-2 space-y-0.5 pl-9 text-sm text-white/90">
            <div>
              Total remaining: {fileCount} file{fileCount === 1 ? "" : "s"} ·{" "}
              {formatBytes(remaining)} ({remainingPct.toFixed(2)}%)
            </div>
            <div>Estimated time remaining: {formatEta(etaSeconds)}</div>
            <div>Transfer rate: {formatRate(rateBps)}</div>
            <div className="font-bold">
              Do not close this tab — closing cancels the upload.
            </div>
          </div>
        )}

        {preparing && (
          <div className="mt-2 pl-9 text-sm text-white/90">
            Preparing {title} — {formatBytes(total)}…
          </div>
        )}

        {status === "done" && (
          <div
            className={`mt-2 pl-9 text-sm text-white/90 ${doneWarning ? "font-bold" : ""}`}
          >
            {doneWarning
              ? `Your files uploaded, but processing didn't start: ${submitWarning}`
              : `${title} uploaded — now processing on the server.`}
          </div>
        )}

        {status === "queued" && (
          <div className="mt-2 space-y-0.5 pl-9 text-sm text-white/90">
            <div>{title} uploaded — queued for processing…</div>
            <div className="font-bold">
              You can close this tab — we&rsquo;ll email you a link when
              it&rsquo;s ready.
            </div>
          </div>
        )}

        {status === "processing" && (
          <div className="mt-2 space-y-0.5 pl-9 text-sm text-white/90">
            <div>
              Processing {title}
              {stage ? ` — ${stage}` : "…"}
            </div>
            <div className="font-bold">
              You can close this tab — we&rsquo;ll email you a link when
              it&rsquo;s ready.
            </div>
          </div>
        )}

        {complete && (
          <div className="mt-2 pl-9 text-sm text-white/90">
            {title} is ready. Open it with “View dataset”. We’ve also emailed
            you the link.
          </div>
        )}

        {status === "error" && (
          <div className="mt-2 pl-9 text-sm font-bold text-white/90">
            {error || "Something went wrong during the upload."}
          </div>
        )}
      </div>
    </div>
  );
}

/**
 * Persistent top banner for global raw uploads. Mounted once in the layout so
 * it survives navigation — a transfer starts on the homepage and this banner
 * follows the user to /explore. Fixed at the top with an in-flow spacer that
 * pushes page content down, AWS-style. Renders one row per in-flight upload, so
 * a combined single-cell + single-molecule upload shows two stacked bars.
 */
export function UploadBar() {
  const jobs = useUploadStore((s) => s.jobs);

  const anyUploading = jobs.some(
    (j) => j.status === "preparing" || j.status === "uploading",
  );

  const barRef = useRef<HTMLDivElement>(null);
  const [barHeight, setBarHeight] = useState(0);

  useEffect(() => {
    const el = barRef.current;

    if (!el) return;
    const ro = new ResizeObserver(() => setBarHeight(el.offsetHeight));

    ro.observe(el);
    setBarHeight(el.offsetHeight);

    return () => ro.disconnect();
  }, [jobs.length]);

  // Warn before the tab is closed mid-transfer — closing aborts the in-flight
  // upload(s) (same guard as any large web upload).
  useEffect(() => {
    if (!anyUploading) return;
    const handler = (e: BeforeUnloadEvent) => {
      e.preventDefault();
      e.returnValue = "";
    };

    window.addEventListener("beforeunload", handler);

    return () => window.removeEventListener("beforeunload", handler);
  }, [anyUploading]);

  if (jobs.length === 0) return null;

  return (
    <>
      {/* In-flow spacer: reserves the fixed bar's height so content is pushed
          down rather than hidden underneath it. */}
      <div aria-hidden className="shrink-0" style={{ height: barHeight }} />

      {/* z-modal sits above the fixed z-0 light-source backgrounds; fixed keeps
          it visible through the whole scroll. Rows stack for multiple uploads. */}
      <div
        ref={barRef}
        className="fixed inset-x-0 top-0 z-[var(--z-modal)] flex w-full flex-col"
        data-testid="upload-bar"
      >
        {jobs.map((job) => (
          <UploadJobRow key={job.id} job={job} />
        ))}
      </div>
    </>
  );
}
