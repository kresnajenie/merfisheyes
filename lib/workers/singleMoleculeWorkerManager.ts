// lib/workers/singleMoleculeWorkerManager.ts
import type { SingleMoleculeWorkerApi } from "./single-molecule.worker";

import * as Comlink from "comlink";

let worker: Worker | null = null;
let workerProxy: Comlink.Remote<SingleMoleculeWorkerApi> | null = null;

/**
 * Get the singleton SingleMolecule worker instance.
 * Creates the worker on first call, returns the same instance on subsequent calls.
 */
export async function getSingleMoleculeWorker(): Promise<
  Comlink.Remote<SingleMoleculeWorkerApi>
> {
  // Only create the worker and proxy if they don't exist yet
  if (!workerProxy) {
    console.log("[SingleMoleculeWorkerManager] Creating singleton worker...");

    worker = new Worker(
      new URL("./single-molecule.worker.ts", import.meta.url),
      { type: "module" },
    );

    // Add error handling for worker creation/loading
    worker.onerror = (event) => {
      console.error("[SingleMoleculeWorkerManager] Worker error:", event);
      // Nullify proxy/worker so retry is possible
      workerProxy = null;
      worker = null;
    };

    workerProxy = Comlink.wrap<SingleMoleculeWorkerApi>(worker);

    console.log("[SingleMoleculeWorkerManager] Singleton worker created");
  }

  // Add a check in case worker creation failed asynchronously
  if (!workerProxy) {
    throw new Error("Failed to initialize singleton SingleMolecule worker.");
  }

  // Always return the same proxy instance
  return workerProxy;
}

/**
 * Terminate the singleton worker.
 * Should only be called when the entire application is shutting down.
 */
export function terminateSingleMoleculeWorker() {
  if (worker) {
    console.log(
      "[SingleMoleculeWorkerManager] Terminating singleton worker...",
    );
    worker.terminate();
    worker = null;
    workerProxy = null;
  }
}

/**
 * Check if the worker is currently initialized
 */
export function isSingleMoleculeWorkerInitialized(): boolean {
  return workerProxy !== null;
}
