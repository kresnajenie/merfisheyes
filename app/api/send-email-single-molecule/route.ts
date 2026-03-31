import { SESClient, SendEmailCommand } from "@aws-sdk/client-ses";
import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";

const ses = new SESClient({
  region: process.env.AWS_SES_REGION || "us-west-2",
  credentials: {
    accessKeyId: process.env.AWS_SES_ACCESS_KEY_ID!,
    secretAccessKey: process.env.AWS_SES_SECRET_ACCESS_KEY!,
  },
});

const FROM_EMAIL = process.env.SES_FROM_EMAIL || "noreply@merfisheyes.com";

export async function POST(request: NextRequest) {
  try {
    const { email, datasetId } = await request.json();

    if (!email || !datasetId) {
      return NextResponse.json(
        { error: "Email and datasetId are required" },
        { status: 400 },
      );
    }

    const dataset = await prisma.dataset.findUnique({
      where: { id: datasetId },
      select: {
        title: true,
        numCells: true,
        numGenes: true,
      },
    });

    if (!dataset) {
      return NextResponse.json({ error: "Dataset not found" }, { status: 404 });
    }

    const baseUrl = process.env.NEXT_PUBLIC_BASE_URL || "http://localhost:3000";
    const datasetUrl = `${baseUrl}/sm-viewer/${datasetId}`;

    const command = new SendEmailCommand({
      Source: FROM_EMAIL,
      Destination: {
        ToAddresses: [email],
      },
      Message: {
        Subject: { Data: "Your Single Molecule Dataset is Ready - MERFISHEYES" },
        Body: {
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
      <h2 style="margin: 0 0 24px 0; font-size: 22px; font-weight: 600; color: #111;">${dataset.title || "Untitled Dataset"}</h2>
      <a href="${datasetUrl}" style="display: inline-block; background: #111; color: #fff; text-decoration: none; padding: 12px 28px; border-radius: 6px; font-weight: 500; font-size: 15px;">Open Dataset &rarr;</a>
    </div>
    <div style="margin-top: 32px; padding: 16px; background: #f8f8f8; border-radius: 6px;">
      <table style="width: 100%; border-collapse: collapse; font-size: 13px;">
        <tr>
          <td style="padding: 4px 0; color: #888;">Molecules</td>
          <td style="padding: 4px 0; color: #333; text-align: right;">${dataset.numCells.toLocaleString()}</td>
        </tr>
        <tr>
          <td style="padding: 4px 0; color: #888;">Genes</td>
          <td style="padding: 4px 0; color: #333; text-align: right;">${dataset.numGenes.toLocaleString()}</td>
        </tr>
      </table>
    </div>
    <div style="margin-top: 24px; font-size: 12px; color: #999; word-break: break-all;">${datasetUrl}</div>
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

    return NextResponse.json({ success: true });
  } catch (error: any) {
    console.error("SES error:", error);

    return NextResponse.json(
      {
        error: "Failed to send email",
        message: error.message,
      },
      { status: 500 },
    );
  }
}
