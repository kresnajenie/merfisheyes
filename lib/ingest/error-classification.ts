import { DatasetFaultCategory } from "@prisma/client";

/**
 * Best-effort classification of a processing failure from its (free-text) error
 * message. The pipeline stores no structured error code, so we pattern-match the
 * message to guess whose fault it was and produce a plain-language hint for the
 * owner's email. An admin can always override the result from the dashboard.
 */
export interface FailureClassification {
  category: DatasetFaultCategory;
  /** Short label for the dashboard, e.g. "Out of memory". */
  label: string;
  /** One or two sentences shown to the owner in the failure email. */
  userHint: string;
}

/** Dashboard label of the out-of-memory rule — the trigger for automatic
 *  compute-tier escalation on retry. */
export const OUT_OF_MEMORY_LABEL = "Out of memory";

interface Rule {
  category: DatasetFaultCategory;
  label: string;
  userHint: string;
  test: RegExp;
}

// Order matters: the first matching rule wins. Platform/infra signals are
// checked before generic data errors so an OOM that surfaces as a KeyError-ish
// message isn't misread as a user problem.
const RULES: Rule[] = [
  // ── Platform / our side ──────────────────────────────────────────────
  {
    category: "PLATFORM",
    label: OUT_OF_MEMORY_LABEL,
    userHint:
      "The dataset ran out of memory while processing on our servers. This is on our side — we've been notified and will look into giving it more resources. You don't need to do anything.",
    test: /out of memory|outofmemory|memoryerror|cannot allocate|exit code 137|oom[\s-]?kill|killed due to memory|memory usage/i,
  },
  {
    category: "PLATFORM",
    label: "Timed out",
    userHint:
      "Processing took longer than allowed and timed out on our side. We've been notified and will look into it — no action needed from you.",
    test: /timed out|timeout|dockertimeouterror|time limit|deadline exceeded/i,
  },
  {
    category: "PLATFORM",
    label: "Job lost / expired",
    userHint:
      "We lost track of this processing job on our infrastructure, so its outcome couldn't be recovered. This is on our side. Please try uploading again, and let us know if it happens twice.",
    test: /no longer known to aws batch|records expire|job record expired|failed to start|cannotpullcontainer|resourceinitialization|no reason reported/i,
  },
  {
    category: "PLATFORM",
    label: "No output produced",
    userHint:
      "Processing finished but didn't produce the expected output. This is on our side — we're investigating. Please try uploading again.",
    test: /produced no manifest|no manifest|no output/i,
  },
  {
    category: "PLATFORM",
    label: "Infrastructure error",
    userHint:
      "An infrastructure error interrupted processing. This is on our side — we've been notified. Please try again shortly.",
    test: /accessdenied|s3 error|network|connection reset|econn|internal error|503|unavailable/i,
  },

  // ── User / the uploaded data ─────────────────────────────────────────
  {
    category: "USER",
    label: "Missing/unknown column",
    userHint:
      "A required column was missing or named unexpectedly in your file. Please check that the gene and X/Y(/Z) coordinate columns are present and correctly named, then re-upload.",
    test: /missing (required )?column|required column|column .* not found|no such column|keyerror|expected column|unknown column|could not find column/i,
  },
  {
    category: "USER",
    label: "Unreadable / corrupt file",
    userHint:
      "We couldn't read or parse your file — it may be corrupt, truncated, or not the format we expected (.h5ad / .parquet / .csv). Please re-export it and upload again.",
    test: /could not parse|failed to parse|not a valid|unsupported format|corrupt|malformed|unreadable|bad magic|unable to read|invalid file|not a parquet|not an? h5ad|decode error|unexpected end of/i,
  },
  {
    category: "USER",
    label: "No valid coordinates",
    userHint:
      "We couldn't find valid spatial coordinates in your file. Please check the coordinate columns contain numeric values, then re-upload.",
    test: /no valid coordinates|missing coordinates|invalid coordinates|no spatial|insufficient valid/i,
  },
  {
    category: "USER",
    label: "Empty / no data",
    userHint:
      "Your file appears to be empty or contains no usable rows. Please check the file and upload again.",
    test: /empty file|no data|0 rows|no rows|no genes found|no cells found|file is empty/i,
  },
  {
    category: "USER",
    label: "Duplicate upload",
    userHint:
      "This dataset is identical to one already in your library, so we didn't create a duplicate. Open the existing one instead.",
    test: /duplicate of already-ingested|already in your library/i,
  },
];

export function classifyFailure(
  errorMessage: string | null | undefined,
): FailureClassification {
  const msg = errorMessage ?? "";

  for (const rule of RULES) {
    if (rule.test.test(msg)) {
      return {
        category: rule.category,
        label: rule.label,
        userHint: rule.userHint,
      };
    }
  }

  return {
    category: "UNKNOWN",
    label: "Needs review",
    userHint:
      "Processing failed for a reason we couldn't automatically categorize. Our team has been notified and will take a look. Feel free to try uploading again, or reply to this email with details.",
  };
}
