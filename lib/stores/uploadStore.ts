import { create } from "zustand";

import {
  uploadRawToS3,
  abortRawUpload,
  pollIngestStatus,
  type RawUploadFile,
  type RawUploadKind,
} from "@/lib/upload/uploadRawToS3";

/**
 * Global raw-upload state.
 *
 * The server-side upload transfers bytes straight from the browser, so it must
 * outlive the page that started it: the user drops a file on the homepage and
 * is sent to /explore, where a top bar shows progress. Keeping the transfer in
 * a component would abort it on that navigation — so it lives here instead, and
 * a single `<UploadBar>` mounted in the layout renders it from anywhere.
 */

export type UploadStatus =
  | "idle"
  | "preparing"
  | "uploading"
  // Bytes are up; the server-side job is now running. The bar keeps tracking it
  // via pollIngestStatus instead of vanishing.
  | "queued"
  | "processing"
  | "complete"
  // "done" is the transient state right after bytes finish but before the first
  // status poll (or the terminal state when the job failed to even submit).
  | "done"
  | "error";

interface StartOpts {
  kind: RawUploadKind;
  title: string;
  files: RawUploadFile[];
  processingParams: Record<string, unknown>;
}

interface UploadStore {
  active: boolean;
  status: UploadStatus;
  loaded: number;
  total: number;
  fileIndex: number;
  fileCount: number;
  title: string;
  datasetId: string | null;
  error: string;
  submitWarning: string;
  /** Bytes/sec, averaged over the transfer so far (0 until measurable). */
  rateBps: number;
  /** Seconds left at the current rate (null until measurable). */
  etaSeconds: number | null;
  /** Server-side processing stage, e.g. "Expression chunks" (blank until known). */
  stage: string;
  /** Processing percent when the worker reports one (null = indeterminate). */
  percent: number | null;
  /** Link into the viewer once processing completes. */
  viewerUrl: string | null;
  _abort: AbortController | null;
  _uploadStartAt: number | null;
  /** Stops the in-flight status poll; set while a job is being tracked. */
  _stopPoll: (() => void) | null;
  start: (opts: StartOpts) => Promise<void>;
  cancel: () => void;
  dismiss: () => void;
}

export const useUploadStore = create<UploadStore>((set, get) => ({
  active: false,
  status: "idle",
  loaded: 0,
  total: 0,
  fileIndex: 0,
  fileCount: 0,
  title: "",
  datasetId: null,
  error: "",
  submitWarning: "",
  rateBps: 0,
  etaSeconds: null,
  stage: "",
  percent: null,
  viewerUrl: null,
  _abort: null,
  _uploadStartAt: null,
  _stopPoll: null,

  async start({ kind, title, files, processingParams }) {
    const controller = new AbortController();

    // A previous job's poll must not outlive the bar being reused.
    get()._stopPoll?.();

    set({
      active: true,
      status: "preparing",
      loaded: 0,
      total: files.reduce((s, f) => s + f.file.size, 0),
      fileIndex: 0,
      fileCount: files.length,
      title,
      datasetId: null,
      error: "",
      submitWarning: "",
      rateBps: 0,
      etaSeconds: null,
      stage: "",
      percent: null,
      viewerUrl: null,
      _abort: controller,
      _uploadStartAt: null,
      _stopPoll: null,
    });

    try {
      const { datasetId, submitError } = await uploadRawToS3({
        kind,
        title,
        processingParams,
        files,
        signal: controller.signal,
        // Stamp the clock when bytes actually start flowing, so the rate isn't
        // skewed by the initiate round-trip that precedes it.
        onPhase: (phase) =>
          set((s) =>
            phase === "uploading" && s._uploadStartAt == null
              ? { status: phase, _uploadStartAt: Date.now() }
              : { status: phase },
          ),
        onProgress: (p) =>
          set((s) => {
            let { rateBps, etaSeconds } = s;

            if (s._uploadStartAt != null && p.loaded > 0) {
              const elapsed = (Date.now() - s._uploadStartAt) / 1000;

              if (elapsed > 0.4) {
                rateBps = p.loaded / elapsed;
                etaSeconds =
                  rateBps > 0 ? (p.total - p.loaded) / rateBps : null;
              }
            }

            return {
              loaded: p.loaded,
              total: p.total,
              fileIndex: p.fileIndex,
              fileCount: p.fileCount,
              rateBps,
              etaSeconds,
            };
          }),
      });

      set({ datasetId, status: "done", submitWarning: submitError ?? "" });

      // The job never made it to the queue — leave the bar up (amber) so the
      // failure is noticed; there's nothing to poll.
      if (submitError) return;

      // Bytes are up and the job is queued. Track it to completion in-place so
      // the bar shows Processing → Ready instead of vanishing. The store is a
      // singleton, so this keeps polling as the user navigates to /explore.
      set({ status: "queued", stage: "", percent: null, viewerUrl: null });

      const stop = pollIngestStatus(datasetId, (s) => {
        // Ignore stale updates if the bar has moved on to another upload.
        if (get().datasetId !== datasetId) return;

        if (s.status === "PROCESSING") {
          set({
            status: "processing",
            stage: s.progress?.stage ?? "",
            percent: s.progress?.percent ?? null,
          });
        } else if (s.status === "COMPLETE") {
          set({
            status: "complete",
            stage: "",
            percent: 100,
            viewerUrl: s.viewerUrl ?? null,
          });
        } else if (s.status === "FAILED") {
          set({
            status: "error",
            error: s.errorMessage || "Processing failed on the server.",
          });
        } else if (s.status === "QUEUED") {
          set({ status: "queued" });
        }
      });

      set({ _stopPoll: stop });
    } catch (err: unknown) {
      if (err instanceof DOMException && err.name === "AbortError") {
        set({ active: false, status: "idle" });

        return;
      }
      set({
        status: "error",
        error: err instanceof Error ? err.message : "Upload failed",
      });
    } finally {
      set({ _abort: null });
    }
  },

  cancel() {
    const { _abort, _stopPoll, datasetId } = get();

    _stopPoll?.();
    _abort?.abort();
    if (datasetId) abortRawUpload(datasetId);
    set({
      active: false,
      status: "idle",
      _abort: null,
      _stopPoll: null,
      datasetId: null,
    });
  },

  dismiss() {
    get()._stopPoll?.();
    set({ active: false, status: "idle", _stopPoll: null });
  },
}));
