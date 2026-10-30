import { SESClient, SendEmailCommand } from "@aws-sdk/client-ses";

import { viewerHref } from "@/lib/datasets/kind";

const ses = new SESClient({
  region: process.env.AWS_SES_REGION || "us-west-2",
  credentials: {
    accessKeyId: process.env.AWS_SES_ACCESS_KEY_ID!,
    secretAccessKey: process.env.AWS_SES_SECRET_ACCESS_KEY!,
  },
});

const FROM_EMAIL = process.env.SES_FROM_EMAIL || "noreply@merfisheyes.com";

export interface DatasetReadyEmailParams {
  email: string;
  datasetId: string;
  /** Its `datasetType`: which viewer the link opens (lib/datasets/kind.ts). */
  datasetType?: string | null;
  /** Whether its count is of molecules (a molecule dataset or store). */
  molecules?: boolean;
  datasetName?: string;
  metadata?: {
    numCells?: number;
    numGenes?: number;
    platform?: string;
    clusterCount?: number;
  };
  /**
   * The other datasets saved with this one (a folder's molecules and image
   * beside its cells), linked below it, so one save sends one email.
   */
  alsoSaved?: {
    datasetId: string;
    datasetType?: string | null;
    datasetName?: string;
  }[];
}

/**
 * Send the "your dataset is ready" email.
 *
 * Extracted from app/api/send-email so the server-side ingestion callback can
 * reuse it directly — a server-to-server HTTP call to our own route gets
 * blocked by Vercel deployment protection.
 */
export async function sendDatasetReadyEmail({
  email,
  datasetId,
  datasetType,
  molecules,
  datasetName,
  metadata,
  alsoSaved = [],
}: DatasetReadyEmailParams): Promise<void> {
  const baseUrl = process.env.NEXT_PUBLIC_APP_URL || "http://localhost:3000";
  const others = alsoSaved.map((d) => ({
    name: d.datasetName || "Untitled",
    link: `${baseUrl}${viewerHref({ id: d.datasetId, datasetType: d.datasetType ?? null })}`,
  }));
  // The viewer of what it is (a SpatialData store, an image, cells or
  // molecules); without a type, an sm_ id is still a molecule dataset
  const isSingleMolecule = molecules ?? datasetId.startsWith("sm_");
  const link = `${baseUrl}${viewerHref({
    id: datasetId,
    datasetType:
      datasetType ?? (datasetId.startsWith("sm_") ? "single_molecule" : null),
  })}`;
  // numCells stores the molecule count for molecule datasets, so label the
  // row accordingly.
  const countLabel = isSingleMolecule ? "Molecules" : "Cells";

  const subject = datasetName
    ? `${datasetName} - Dataset Ready - MERFISHEYES`
    : "Your Dataset is Ready - MERFISHEYES";

  const command = new SendEmailCommand({
    Source: FROM_EMAIL,
    Destination: {
      ToAddresses: [email],
    },
    Message: {
      Subject: { Data: subject },
      Body: {
        Text: {
          Data: `Your dataset is ready.\n\nDataset: ${datasetName || "Untitled"}\n${countLabel}: ${metadata?.numCells?.toLocaleString() || "N/A"}\nGenes: ${metadata?.numGenes?.toLocaleString() || "N/A"}\nPlatform: ${metadata?.platform || "N/A"}\n\nOpen your dataset: ${link}${others.length ? `\n\nAlso in this save:\n${others.map((o) => `${o.name}: ${o.link}`).join("\n")}` : ""}`,
        },
        Html: {
          Data: `<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
  </head>
  <body style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; line-height: 1.6; color: #333; max-width: 520px; margin: 0 auto; padding: 20px;">
    <div style="padding: 24px 0 16px 0;">
      <span style="font-size: 14px; font-weight: 600; color: #111; letter-spacing: 0.5px;">MERFISHEYES</span>
    </div>
    <div style="border-top: 1px solid #e0e0e0; padding-top: 24px;">
      <p style="margin: 0 0 4px 0; font-size: 14px; color: #666;">Your dataset is ready</p>
      <h2 style="margin: 0 0 24px 0; font-size: 22px; font-weight: 600; color: #111;">${datasetName || "Untitled Dataset"}</h2>
      <a href="${link}" style="display: inline-block; background: #2563eb; color: #fff; text-decoration: none; padding: 12px 28px; border-radius: 6px; font-weight: 500; font-size: 15px;">Open Dataset &rarr;</a>${
        others.length
          ? `
      <p style="margin: 24px 0 4px 0; font-size: 13px; color: #666;">Also in this save</p>${others
        .map(
          (o) => `
      <p style="margin: 0; font-size: 14px;"><a href="${o.link}" style="color: #2563eb; text-decoration: none;">${escapeHtml(o.name)}</a></p>`,
        )
        .join("")}`
          : ""
      }
    </div>
    <div style="margin-top: 32px; padding: 16px; background: #f8f8f8; border-radius: 6px;">
      <table style="width: 100%; border-collapse: collapse; font-size: 13px;">
        <tr>
          <td style="padding: 4px 0; color: #888;">${countLabel}</td>
          <td style="padding: 4px 0; color: #333; text-align: right;">${metadata?.numCells?.toLocaleString() || "N/A"}</td>
        </tr>
        <tr>
          <td style="padding: 4px 0; color: #888;">Genes</td>
          <td style="padding: 4px 0; color: #333; text-align: right;">${metadata?.numGenes?.toLocaleString() || "N/A"}</td>
        </tr>
        <tr>
          <td style="padding: 4px 0; color: #888;">Platform</td>
          <td style="padding: 4px 0; color: #333; text-align: right; text-transform: uppercase;">${metadata?.platform || "N/A"}</td>
        </tr>${
          metadata?.clusterCount
            ? `
        <tr>
          <td style="padding: 4px 0; color: #888;">Clusters</td>
          <td style="padding: 4px 0; color: #333; text-align: right;">${metadata.clusterCount}</td>
        </tr>`
            : ""
        }
      </table>
    </div>
    <div style="margin-top: 24px; font-size: 12px; color: #999; word-break: break-all;">${link}</div>
    <div style="margin-top: 32px; padding-top: 16px; border-top: 1px solid #e0e0e0;">
      <table style="width: 100%;">
        <tr>
          <td style="font-size: 11px; color: #bbb;">&copy; ${new Date().getFullYear()} MERFISHEYES</td>
          <td style="text-align: right; font-size: 12px;">
            <a href="https://github.com/kresnajenie/merfisheyes" style="color: #999; text-decoration: none; margin-left: 12px;">GitHub</a>
            <a href="https://x.com/merfisheyes" style="color: #999; text-decoration: none; margin-left: 12px;">X</a>
            <a href="https://discord.gg/BRp6C2EVHU" style="color: #999; text-decoration: none; margin-left: 12px;">Discord</a>
            <a href="https://www.patreon.com/cw/MERFISHEYES" style="color: #999; text-decoration: none; margin-left: 12px;">Patreon</a>
          </td>
        </tr>
      </table>
    </div>
  </body>
</html>`,
        },
      },
    },
  });

  await ses.send(command);
}

export interface DatasetFailedEmailParams {
  email: string;
  datasetId: string;
  datasetName?: string;
  /** The raw processing error, shown verbatim so the user can act / forward it. */
  errorMessage: string;
  /** Plain-language hint about whose fault it is and what to do next. */
  userHint: string;
}

/**
 * Tell the owner their dataset failed to process. Shows the raw error text
 * (so they — or we — can act on the specifics) plus a friendly hint about
 * whether it's likely their file or our side.
 */
export async function sendDatasetFailedEmail({
  email,
  datasetId,
  datasetName,
  errorMessage,
  userHint,
}: DatasetFailedEmailParams): Promise<void> {
  const baseUrl = process.env.NEXT_PUBLIC_APP_URL || "http://localhost:3000";
  const uploadLink = `${baseUrl}/`;
  const name = datasetName || "Untitled Dataset";
  // Keep the shown error bounded even though the DB column already caps it.
  const shownError = (errorMessage || "Processing failed").slice(0, 1500);

  const command = new SendEmailCommand({
    Source: FROM_EMAIL,
    Destination: { ToAddresses: [email] },
    Message: {
      Subject: { Data: `${name} — processing failed — MERFISHEYES` },
      Body: {
        Text: {
          Data: `Your dataset "${name}" failed to process.\n\nWhat happened:\n${shownError}\n\n${userHint}\n\nReference: ${datasetId}\nTry again: ${uploadLink}`,
        },
        Html: {
          Data: `<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
  </head>
  <body style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; line-height: 1.6; color: #333; max-width: 520px; margin: 0 auto; padding: 20px;">
    <div style="padding: 24px 0 16px 0;">
      <span style="font-size: 14px; font-weight: 600; color: #111; letter-spacing: 0.5px;">MERFISHEYES</span>
    </div>
    <div style="border-top: 1px solid #e0e0e0; padding-top: 24px;">
      <p style="margin: 0 0 4px 0; font-size: 14px; color: #b91c1c;">Processing failed</p>
      <h2 style="margin: 0 0 16px 0; font-size: 22px; font-weight: 600; color: #111;">${name}</h2>
      <p style="margin: 0 0 20px 0; font-size: 14px; color: #555;">${userHint}</p>
      <div style="margin: 0 0 20px 0; padding: 12px 14px; background: #fef2f2; border: 1px solid #fecaca; border-radius: 6px;">
        <p style="margin: 0 0 6px 0; font-size: 12px; color: #991b1b; font-weight: 600;">Error details</p>
        <pre style="margin: 0; font-size: 12px; color: #7f1d1d; white-space: pre-wrap; word-break: break-word; font-family: ui-monospace, SFMono-Regular, Menlo, monospace;">${escapeHtml(shownError)}</pre>
      </div>
      <a href="${uploadLink}" style="display: inline-block; background: #2563eb; color: #fff; text-decoration: none; padding: 12px 28px; border-radius: 6px; font-weight: 500; font-size: 15px;">Try uploading again &rarr;</a>
    </div>
    <div style="margin-top: 24px; font-size: 12px; color: #999;">Reference: ${datasetId}</div>
    <div style="margin-top: 32px; padding-top: 16px; border-top: 1px solid #e0e0e0; font-size: 12px; color: #999;">
      Reply to this email if you'd like a hand — we're happy to take a look.
    </div>
  </body>
</html>`,
        },
      },
    },
  });

  await ses.send(command);
}

/** Minimal HTML-escape for embedding untrusted error text in the email body. */
function escapeHtml(s: string): string {
  return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

export interface VerificationCodeEmailParams {
  email: string;
  code: string;
  /** Minutes until the code stops working, quoted in the email. */
  expiresInMinutes: number;
}

/**
 * Send a sign-in code.
 *
 * Deliberately does NOT include a clickable sign-in link. The code is short
 * enough to retype, and a link in an email gets prefetched by scanners and
 * corporate link-checkers — which would consume the single-use token before
 * the person ever clicks it.
 */
export async function sendVerificationCodeEmail({
  email,
  code,
  expiresInMinutes,
}: VerificationCodeEmailParams): Promise<void> {
  const command = new SendEmailCommand({
    Source: FROM_EMAIL,
    Destination: {
      ToAddresses: [email],
    },
    Message: {
      Subject: { Data: `${code} is your MERFISHEYES sign-in code` },
      Body: {
        Text: {
          Data: `Your MERFISHEYES sign-in code is ${code}\n\nIt expires in ${expiresInMinutes} minutes and can only be used once.\n\nIf you didn't try to sign in, you can ignore this email.`,
        },
        Html: {
          Data: `<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
  </head>
  <body style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; line-height: 1.6; color: #333; max-width: 520px; margin: 0 auto; padding: 20px;">
    <div style="padding: 24px 0 16px 0;">
      <span style="font-size: 14px; font-weight: 600; color: #111; letter-spacing: 0.5px;">MERFISHEYES</span>
    </div>
    <div style="border-top: 1px solid #e0e0e0; padding-top: 24px;">
      <p style="margin: 0 0 16px 0; font-size: 14px; color: #666;">Your sign-in code</p>
      <div style="font-size: 34px; font-weight: 600; letter-spacing: 8px; color: #111; padding: 16px 0;">${code}</div>
      <p style="margin: 16px 0 0 0; font-size: 13px; color: #888;">Expires in ${expiresInMinutes} minutes. Can only be used once.</p>
    </div>
    <div style="margin-top: 32px; padding-top: 16px; border-top: 1px solid #e0e0e0; font-size: 12px; color: #999;">
      If you didn't try to sign in, you can ignore this email.
    </div>
  </body>
</html>`,
        },
      },
    },
  });

  await ses.send(command);
}

export interface DatasetRetryingEmailParams {
  email: string;
  datasetId: string;
  datasetName?: string;
  /** Why the job is being re-run. */
  reason: "out_of_memory" | "admin";
  /** Memory (GB) of the tier that failed and of the tier being tried now. */
  previousMemoryGb: number;
  nextMemoryGb: number;
}

/**
 * Tell the owner their dataset is being processed again — either because it
 * ran out of memory and we're automatically retrying on a larger machine, or
 * because an admin restarted it. No action needed from them; the usual
 * ready/failed email follows when the new run finishes.
 */
export async function sendDatasetRetryingEmail({
  email,
  datasetId,
  datasetName,
  reason,
  previousMemoryGb,
  nextMemoryGb,
}: DatasetRetryingEmailParams): Promise<void> {
  const name = datasetName || "Untitled Dataset";
  const what =
    reason === "out_of_memory"
      ? `Your dataset ran out of memory while processing on our ${previousMemoryGb} GB tier. We're automatically retrying it on a ${nextMemoryGb} GB machine.`
      : previousMemoryGb === nextMemoryGb
        ? `We've restarted processing for your dataset.`
        : `We've restarted processing for your dataset on a larger ${nextMemoryGb} GB machine (it previously ran on ${previousMemoryGb} GB).`;
  const followUp =
    "No action is needed — you'll get another email when it's ready, or if it fails again.";

  const command = new SendEmailCommand({
    Source: FROM_EMAIL,
    Destination: { ToAddresses: [email] },
    Message: {
      Subject: { Data: `${name} — retrying processing — MERFISHEYES` },
      Body: {
        Text: {
          Data: `${what}\n\n${followUp}\n\nReference: ${datasetId}`,
        },
        Html: {
          Data: `<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
  </head>
  <body style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; line-height: 1.6; color: #333; max-width: 520px; margin: 0 auto; padding: 20px;">
    <div style="padding: 24px 0 16px 0;">
      <span style="font-size: 14px; font-weight: 600; color: #111; letter-spacing: 0.5px;">MERFISHEYES</span>
    </div>
    <div style="border-top: 1px solid #e0e0e0; padding-top: 24px;">
      <p style="margin: 0 0 4px 0; font-size: 14px; color: #b45309;">Retrying processing</p>
      <h2 style="margin: 0 0 16px 0; font-size: 22px; font-weight: 600; color: #111;">${escapeHtml(name)}</h2>
      <p style="margin: 0 0 12px 0; font-size: 14px; color: #555;">${escapeHtml(what)}</p>
      <p style="margin: 0 0 20px 0; font-size: 14px; color: #555;">${followUp}</p>
    </div>
    <div style="margin-top: 24px; font-size: 12px; color: #999;">Reference: ${datasetId}</div>
  </body>
</html>`,
        },
      },
    },
  });

  await ses.send(command);
}
