import { create } from "zustand";

import {
  uploadRawToS3,
  abortRawUpload,
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
  _abort: AbortController | null;
  _uploadStartAt: number | null;
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
  _abort: null,
  _uploadStartAt: null,

  async start({ kind, title, files, processingParams }) {
    const controller = new AbortController();

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
      _abort: controller,
      _uploadStartAt: null,
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

      // The bytes are up; the dataset now processes server-side and shows in
      // "Your uploads" on /explore. Auto-clear the bar shortly after a clean
      // finish (keep it up if the job failed to submit, so it's noticed).
      if (!submitError) {
        setTimeout(() => {
          if (get().status === "done" && get().datasetId === datasetId) {
            set({ active: false, status: "idle" });
          }
        }, 5000);
      }
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
    const { _abort, datasetId } = get();

    _abort?.abort();
    if (datasetId) abortRawUpload(datasetId);
    set({ active: false, status: "idle", _abort: null, datasetId: null });
  },

  dismiss() {
    set({ active: false, status: "idle" });
  },
}));
