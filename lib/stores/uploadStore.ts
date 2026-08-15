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
  /** Called with the dataset id as soon as the upload is initiated. */
  onInitiated?: (datasetId: string) => void;
}

/**
 * One tracked upload. Several can run at once — e.g. a combined single-cell +
 * single-molecule folder starts two, each rendered as its own bar. Every job
 * owns its own abort controller and status poll so they're fully independent.
 */
export interface UploadJob {
  /** Local id used to key + target this job (not the dataset id). */
  id: string;
  kind: RawUploadKind;
  title: string;
  status: UploadStatus;
  loaded: number;
  total: number;
  fileIndex: number;
  fileCount: number;
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
}

interface UploadStore {
  jobs: UploadJob[];
  start: (opts: StartOpts) => Promise<void>;
  cancel: (id: string) => void;
  dismiss: (id: string) => void;
}

let jobSeq = 0;

export const useUploadStore = create<UploadStore>((set, get) => ({
  jobs: [],

  async start({ kind, title, files, processingParams, onInitiated }) {
    const id = `up_local_${++jobSeq}`;
    const controller = new AbortController();

    // Update just this job, leaving any other in-flight jobs untouched.
    const patch = (partial: Partial<UploadJob>) =>
      set((s) => ({
        jobs: s.jobs.map((j) => (j.id === id ? { ...j, ...partial } : j)),
      }));
    const getJob = () => get().jobs.find((j) => j.id === id);

    set((s) => ({
      jobs: [
        ...s.jobs,
        {
          id,
          kind,
          title,
          status: "preparing",
          loaded: 0,
          total: files.reduce((a, f) => a + f.file.size, 0),
          fileIndex: 0,
          fileCount: files.length,
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
        },
      ],
    }));

    try {
      const { datasetId, submitError } = await uploadRawToS3({
        kind,
        title,
        processingParams,
        files,
        signal: controller.signal,
        // Surface the id early (before bytes finish) and record it on the job.
        onInitiated: (dsId) => {
          patch({ datasetId: dsId });
          onInitiated?.(dsId);
        },
        // Stamp the clock when bytes actually start flowing, so the rate isn't
        // skewed by the initiate round-trip that precedes it.
        onPhase: (phase) => {
          const j = getJob();

          if (phase === "uploading" && j && j._uploadStartAt == null) {
            patch({ status: phase, _uploadStartAt: Date.now() });
          } else {
            patch({ status: phase });
          }
        },
        onProgress: (p) => {
          const j = getJob();
          let rateBps = j?.rateBps ?? 0;
          let etaSeconds = j?.etaSeconds ?? null;

          if (j?._uploadStartAt != null && p.loaded > 0) {
            const elapsed = (Date.now() - j._uploadStartAt) / 1000;

            if (elapsed > 0.4) {
              rateBps = p.loaded / elapsed;
              etaSeconds = rateBps > 0 ? (p.total - p.loaded) / rateBps : null;
            }
          }

          patch({
            loaded: p.loaded,
            total: p.total,
            fileIndex: p.fileIndex,
            fileCount: p.fileCount,
            rateBps,
            etaSeconds,
          });
        },
      });

      patch({ datasetId, status: "done", submitWarning: submitError ?? "" });

      // The job never made it to the queue — leave the bar up (amber) so the
      // failure is noticed; there's nothing to poll.
      if (submitError) return;

      // Bytes are up and the job is queued. Track it to completion so the bar
      // shows Processing → Ready instead of vanishing. Lives in the store, so
      // polling continues as the user navigates to /explore.
      patch({ status: "queued", stage: "", percent: null, viewerUrl: null });

      const stop = pollIngestStatus(datasetId, (s) => {
        // Job was cancelled/dismissed — stop touching it.
        if (!getJob()) return;

        if (s.status === "PROCESSING") {
          patch({
            status: "processing",
            stage: s.progress?.stage ?? "",
            percent: s.progress?.percent ?? null,
          });
        } else if (s.status === "COMPLETE") {
          patch({
            status: "complete",
            stage: "",
            percent: 100,
            viewerUrl: s.viewerUrl ?? null,
          });
        } else if (s.status === "FAILED") {
          patch({
            status: "error",
            error: s.errorMessage || "Processing failed on the server.",
          });
        } else if (s.status === "QUEUED") {
          patch({ status: "queued" });
        }
      });

      patch({ _stopPoll: stop });
    } catch (err: unknown) {
      if (err instanceof DOMException && err.name === "AbortError") {
        set((s) => ({ jobs: s.jobs.filter((j) => j.id !== id) }));

        return;
      }
      patch({
        status: "error",
        error: err instanceof Error ? err.message : "Upload failed",
      });
    } finally {
      patch({ _abort: null });
    }
  },

  cancel(id) {
    const job = get().jobs.find((j) => j.id === id);

    if (!job) return;
    job._stopPoll?.();
    job._abort?.abort();
    if (job.datasetId) abortRawUpload(job.datasetId);
    set((s) => ({ jobs: s.jobs.filter((j) => j.id !== id) }));
  },

  dismiss(id) {
    get()
      .jobs.find((j) => j.id === id)
      ?._stopPoll?.();
    set((s) => ({ jobs: s.jobs.filter((j) => j.id !== id) }));
  },
}));
