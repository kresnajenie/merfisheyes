import { nextComputeTier, type ComputeTier } from "@/lib/batch";
import { OUT_OF_MEMORY_LABEL } from "@/lib/ingest/error-classification";

/**
 * Decide whether a failed processing run should be retried automatically on
 * a bigger machine instead of being marked FAILED.
 *
 * Only out-of-memory failures escalate (a malformed file won't get better
 * with more RAM), only one step up the tier ladder at a time, and only while
 * there is a bigger tier left — standard → large → xlarge, so at most two
 * automatic retries per dataset.
 */
export const MAX_AUTO_ATTEMPTS = 3;

export function planOomEscalation(input: {
  /** classifyFailure(...).label for the failure being recorded. */
  label: string;
  /** Tier the failed job ran on (stored on the dataset; null = legacy row). */
  currentTier: string | null | undefined;
  /** Jobs submitted so far for this dataset (the failed one included). */
  attempts: number;
}): { nextTier: ComputeTier } | null {
  if (input.label !== OUT_OF_MEMORY_LABEL) return null;
  if (input.attempts >= MAX_AUTO_ATTEMPTS) return null;
  const nextTier = nextComputeTier(input.currentTier);

  return nextTier ? { nextTier } : null;
}
