/**
 * Debounced URL writer that coordinates updates from both panels.
 * Uses window.history.replaceState to avoid triggering React re-renders.
 */

let pendingLeft: string | null | undefined = undefined; // undefined = no pending change
let pendingRight: string | null | undefined = undefined;
let debounceTimer: ReturnType<typeof setTimeout> | null = null;

const DEBOUNCE_MS = 400;

function flushToUrl() {
  if (typeof window === "undefined") return;

  const url = new URL(window.location.href);

  if (pendingLeft !== undefined) {
    if (pendingLeft === null) {
      url.searchParams.delete("v");
    } else {
      url.searchParams.set("v", pendingLeft);
    }
    pendingLeft = undefined;
  }

  if (pendingRight !== undefined) {
    if (pendingRight === null) {
      url.searchParams.delete("rv");
    } else {
      url.searchParams.set("rv", pendingRight);
    }
    pendingRight = undefined;
  }

  window.history.replaceState(null, "", url.toString());
}

export function scheduleUrlUpdate(
  panel: "left" | "right",
  encoded: string | null,
) {
  if (panel === "left") {
    pendingLeft = encoded;
  } else {
    pendingRight = encoded;
  }

  if (debounceTimer) clearTimeout(debounceTimer);
  debounceTimer = setTimeout(flushToUrl, DEBOUNCE_MS);
}

export function readUrlVizState(): {
  left: string | null;
  right: string | null;
} {
  if (typeof window === "undefined") return { left: null, right: null };

  const params = new URLSearchParams(window.location.search);

  return {
    left: params.get("v"),
    right: params.get("rv"),
  };
}
