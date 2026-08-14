"use client";

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

/**
 * Persistent top banner for the global raw upload. Mounted once in the layout
 * so it survives navigation — the transfer starts on the homepage and this
 * banner follows the user to /explore. Fixed at the top so it stays visible
 * through the whole page scroll, with an in-flow spacer that pushes the page
 * content down, AWS-style.
 */
export function UploadBar() {
  const active = useUploadStore((s) => s.active);
  const status = useUploadStore((s) => s.status);
  const loaded = useUploadStore((s) => s.loaded);
  const total = useUploadStore((s) => s.total);
  const fileCount = useUploadStore((s) => s.fileCount);
  const title = useUploadStore((s) => s.title);
  const error = useUploadStore((s) => s.error);
  const submitWarning = useUploadStore((s) => s.submitWarning);
  const rateBps = useUploadStore((s) => s.rateBps);
  const etaSeconds = useUploadStore((s) => s.etaSeconds);
  const stage = useUploadStore((s) => s.stage);
  const percent = useUploadStore((s) => s.percent);
  const viewerUrl = useUploadStore((s) => s.viewerUrl);
  const cancel = useUploadStore((s) => s.cancel);
  const dismiss = useUploadStore((s) => s.dismiss);

  const uploading = status === "preparing" || status === "uploading";
  // Server-side processing phase: bytes are up, the worker is running.
  const processingPhase = status === "queued" || status === "processing";
  const complete = status === "complete";
  // Anything with a spinning indicator (bytes moving or the server working).
  const busy = uploading || processingPhase;

  // The bar is `fixed` (stays visible through the whole page scroll), so an
  // in-flow spacer of the same height reserves its space and pushes the page
  // content down — matching the AWS layout without hiding anything under it.
  // The layout wrapper is `h-screen`, which would clip a `sticky` bar after one
  // viewport of scroll, so `fixed` + spacer is used instead.
  const barRef = useRef<HTMLDivElement>(null);
  const [barHeight, setBarHeight] = useState(0);

  useEffect(() => {
    const el = barRef.current;

    if (!el) return;
    const ro = new ResizeObserver(() => setBarHeight(el.offsetHeight));

    ro.observe(el);
    setBarHeight(el.offsetHeight);

    return () => ro.disconnect();
  }, [active, status]);

  // Warn before the tab is closed mid-transfer — closing aborts the in-flight
  // upload (same guard as any large web upload).
  useEffect(() => {
    if (!uploading) return;
    const handler = (e: BeforeUnloadEvent) => {
      e.preventDefault();
      e.returnValue = "";
    };

    window.addEventListener("beforeunload", handler);

    return () => window.removeEventListener("beforeunload", handler);
  }, [uploading]);

  if (!active) return null;

  const pct = total > 0 ? Math.min(100, (loaded / total) * 100) : 0;
  const remaining = Math.max(0, total - loaded);
  const remainingPct = Math.max(0, 100 - pct);
  const doneWarning = status === "done" && Boolean(submitWarning);
  const preparing = status === "preparing";
  // The processing bar is indeterminate unless the worker reports a percent.
  const hasPercent = processingPhase && percent != null;

  const bg =
    status === "error"
      ? "bg-danger"
      : doneWarning
        ? "bg-amber-600"
        : complete
          ? // Fixed Tailwind green, not the HeroUI `success-*` scale — that scale
            // inverts in dark mode (higher number → lighter), which made the bar
            // paler instead of darker.
            "bg-emerald-700"
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
    <>
      {/* In-flow spacer: reserves the fixed bar's height so content is pushed
          down rather than hidden underneath it. `shrink-0` because it's an
          empty flex item in the layout's flex-col — flex would otherwise
          collapse it to 0 and the bar would overlap the navbar. */}
      <div aria-hidden className="shrink-0" style={{ height: barHeight }} />

      {/* z-modal sits above the fixed z-0 light-source backgrounds (LightRays /
          ExploreBackground); fixed keeps it visible through the whole scroll. */}
      <div
        ref={barRef}
        className={`fixed inset-x-0 top-0 z-[var(--z-modal)] w-full text-white ${bg}`}
        data-testid="upload-bar"
      >
        <div className="mx-auto max-w-6xl px-6 py-4">
          <div className="flex items-center gap-4">
            <Icon />
            <span className="shrink-0 text-lg font-semibold">{heading}</span>

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
                {(complete ? 100 : hasPercent ? (percent ?? 0) : pct).toFixed(
                  0,
                )}
                %
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
              onClick={uploading ? cancel : dismiss}
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
    </>
  );
}
