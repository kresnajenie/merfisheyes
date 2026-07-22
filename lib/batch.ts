import { BatchClient, SubmitJobCommand } from "@aws-sdk/client-batch";

// Submits the server-side ingestion job to AWS Batch (Fargate). The app never
// blocks on processing — SubmitJob returns immediately and all status flows back
// through /api/ingest/[id]/callback (design §5.4, and Vercel functions are
// short-lived anyway).

const batch = new BatchClient({
  region: process.env.AWS_REGION || "us-west-2",
  ...(process.env.AWS_ACCESS_KEY_ID
    ? {
        credentials: {
          accessKeyId: process.env.AWS_ACCESS_KEY_ID,
          secretAccessKey: process.env.AWS_SECRET_ACCESS_KEY || "",
        },
      }
    : {}),
});

const JOB_QUEUE = process.env.AWS_BATCH_JOB_QUEUE || "merfisheyes-ingest";
const JOB_DEFINITION =
  process.env.AWS_BATCH_JOB_DEFINITION || "merfisheyes-ingest";

/**
 * Per-user compute sizing (design §5.6). Fargate only accepts specific
 * vCPU↔memory pairs (4 vCPU → 8–30 GB, 8 → 16–60 GB, 16 → 32–120 GB), so tiers
 * must snap to valid combos.
 */
const COMPUTE_TIERS: Record<string, { vcpu: string; memoryMiB: string }> = {
  standard: { vcpu: "4", memoryMiB: "16384" },
  large: { vcpu: "8", memoryMiB: "32768" },
  xlarge: { vcpu: "16", memoryMiB: "65536" },
};

export interface SubmitIngestJobParams {
  datasetId: string;
  kind: "single_cell" | "single_molecule";
  processingParams: unknown;
  computeTier?: string;
}

/** True when the app has enough config to submit jobs at all. */
export function isBatchConfigured(): boolean {
  return Boolean(process.env.CALLBACK_SECRET && callbackBaseUrl());
}

/**
 * Vercel injects this once "Protection Bypass for Automation" is configured.
 * VERCEL_BYPASS_SECRET is accepted as a manual override (e.g. local testing).
 */
function bypassSecret(): string | undefined {
  return (
    process.env.VERCEL_AUTOMATION_BYPASS_SECRET ||
    process.env.VERCEL_BYPASS_SECRET ||
    undefined
  );
}

/**
 * Public origin Fargate POSTs callbacks to.
 *
 * Falls back to Vercel's own system vars so **preview deployments work without
 * any config** — VERCEL_BRANCH_URL is the stable per-branch alias, VERCEL_URL
 * the per-deployment one. Set INGEST_CALLBACK_BASE_URL to override (e.g. to
 * pin callbacks at a specific domain, or a tunnel during local testing).
 */
function callbackBaseUrl(): string {
  const explicit =
    process.env.INGEST_CALLBACK_BASE_URL ||
    process.env.NEXT_PUBLIC_APP_URL ||
    process.env.NEXT_PUBLIC_BASE_URL;

  if (explicit) return explicit.replace(/\/$/, "");

  const vercelHost = process.env.VERCEL_BRANCH_URL || process.env.VERCEL_URL;

  return vercelHost ? `https://${vercelHost.replace(/\/$/, "")}` : "";
}

/**
 * Submit the processing job. Returns the Batch job id to store on the Dataset.
 */
export async function submitIngestJob({
  datasetId,
  kind,
  processingParams,
  computeTier = "standard",
}: SubmitIngestJobParams): Promise<string> {
  const base = callbackBaseUrl();

  if (!base) {
    throw new Error(
      "No callback base URL configured (INGEST_CALLBACK_BASE_URL / NEXT_PUBLIC_APP_URL)",
    );
  }
  if (!process.env.CALLBACK_SECRET) {
    throw new Error("CALLBACK_SECRET is not configured");
  }

  const tier = COMPUTE_TIERS[computeTier] ?? COMPUTE_TIERS.standard;

  const command = new SubmitJobCommand({
    jobName: `ingest-${datasetId}`.slice(0, 128),
    jobQueue: JOB_QUEUE,
    jobDefinition: JOB_DEFINITION,
    containerOverrides: {
      resourceRequirements: [
        { type: "VCPU", value: tier.vcpu },
        { type: "MEMORY", value: tier.memoryMiB },
      ],
      environment: [
        { name: "DATASET_ID", value: datasetId },
        { name: "DATASET_KIND", value: kind },
        { name: "PROCESSING_PARAMS", value: JSON.stringify(processingParams ?? {}) },
        { name: "CALLBACK_URL", value: `${base}/api/ingest/${datasetId}/callback` },
        { name: "CALLBACK_SECRET", value: process.env.CALLBACK_SECRET },
        // Lets the worker through Vercel Deployment Protection. Vercel injects
        // VERCEL_AUTOMATION_BYPASS_SECRET once a bypass secret is configured;
        // absent it, the callback would 401 on a protected deployment.
        ...(bypassSecret()
          ? [{ name: "VERCEL_BYPASS_SECRET", value: bypassSecret() }]
          : []),
      ],
    },
  });

  const result = await batch.send(command);

  if (!result.jobId) throw new Error("Batch SubmitJob returned no jobId");

  return result.jobId;
}
