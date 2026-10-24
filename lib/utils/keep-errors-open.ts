import type { toast as Toast } from "react-toastify";

/**
 * An error stays on screen until it is closed, whatever a call site asked
 * for: whoever has to act on it must get to read it (a failed upload's
 * reason used to vanish after a few seconds). Covers both ways an error
 * toast is made — `toast.error(…)`, and `toast.update(id, { type: "error" })`
 * on a loading toast that failed.
 */
export function keepErrorsOpen(toast: typeof Toast): void {
  const error = toast.error;
  const update = toast.update;

  toast.error = ((content, options) =>
    error(content, { ...options, autoClose: false })) as typeof toast.error;
  toast.update = ((id, options) =>
    update(
      id,
      options?.type === "error" ? { ...options, autoClose: false } : options,
    )) as typeof toast.update;
}
