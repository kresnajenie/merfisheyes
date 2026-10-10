/**
 * Debounced URL writer that coordinates updates from the left panel (v=),
 * the right panel (rv=), and the SM overlay on the SC viewer of either panel
 * (ov= on the left, rov= on the right).
 * Uses window.history.replaceState to avoid triggering React re-renders.
 */

let pendingLeft: string | null | undefined = undefined; // undefined = no pending change
let pendingRight: string | null | undefined = undefined;
let pendingOverlay: string | null | undefined = undefined;
let pendingRightOverlay: string | null | undefined = undefined;
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

  if (pendingOverlay !== undefined) {
    if (pendingOverlay === null) {
      url.searchParams.delete("ov");
    } else {
      url.searchParams.set("ov", pendingOverlay);
    }
    pendingOverlay = undefined;
  }

  if (pendingRightOverlay !== undefined) {
    if (pendingRightOverlay === null) {
      url.searchParams.delete("rov");
    } else {
      url.searchParams.set("rov", pendingRightOverlay);
    }
    pendingRightOverlay = undefined;
  }

  window.history.replaceState(null, "", url.toString());
}

export function scheduleUrlUpdate(
  panel: "left" | "right" | "overlay" | "rightOverlay",
  encoded: string | null,
) {
  if (panel === "left") {
    pendingLeft = encoded;
  } else if (panel === "right") {
    pendingRight = encoded;
  } else if (panel === "rightOverlay") {
    pendingRightOverlay = encoded;
  } else {
    pendingOverlay = encoded;
  }

  if (debounceTimer) clearTimeout(debounceTimer);
  debounceTimer = setTimeout(flushToUrl, DEBOUNCE_MS);
}

export function readUrlVizState(): {
  left: string | null;
  right: string | null;
  overlay: string | null;
  rightOverlay: string | null;
} {
  if (typeof window === "undefined")
    return { left: null, right: null, overlay: null, rightOverlay: null };

  const params = new URLSearchParams(window.location.search);

  return {
    left: params.get("v"),
    right: params.get("rv"),
    overlay: params.get("ov"),
    rightOverlay: params.get("rov"),
  };
}
