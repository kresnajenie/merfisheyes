// lib/workers/standardizedDatasetWorkerManager.ts
import type { StandardizedDatasetWorkerApi } from "./standardized-dataset.worker";

import * as Comlink from "comlink";

let worker: Worker | null = null;
let workerProxy: Comlink.Remote<StandardizedDatasetWorkerApi> | null = null;

/**
 * Get the singleton StandardizedDataset worker instance.
 * Creates the worker on first call, returns the same instance on subsequent calls.
 */
export async function getStandardizedDatasetWorker(): Promise<
  Comlink.Remote<StandardizedDatasetWorkerApi>
> {
  // Only create the worker and proxy if they don't exist yet
  if (!workerProxy) {
    console.log(
      "[StandardizedDatasetWorkerManager] Creating singleton worker...",
    );

    worker = new Worker(
      new URL("./standardized-dataset.worker.ts", import.meta.url),
      { type: "module" },
    );

    // Add error handling for worker creation/loading
    worker.onerror = (event) => {
      console.error("[StandardizedDatasetWorkerManager] Worker error:", event);
      // Nullify proxy/worker so retry is possible
      workerProxy = null;
      worker = null;
    };

    workerProxy = Comlink.wrap<StandardizedDatasetWorkerApi>(worker);

    console.log("[StandardizedDatasetWorkerManager] Singleton worker created");
  }

  // Add a check in case worker creation failed asynchronously
  if (!workerProxy) {
    throw new Error(
      "Failed to initialize singleton StandardizedDataset worker.",
    );
  }

  // Always return the same proxy instance
  return workerProxy;
}

/**
 * Terminate the singleton worker.
 * Should only be called when the entire application is shutting down.
 */
export function terminateStandardizedDatasetWorker() {
  if (worker) {
    console.log(
      "[StandardizedDatasetWorkerManager] Terminating singleton worker...",
    );
    worker.terminate();
    worker = null;
    workerProxy = null;
  }
}

/**
 * Check if the worker is currently initialized
 */
export function isStandardizedDatasetWorkerInitialized(): boolean {
  return workerProxy !== null;
}
