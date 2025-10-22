import sgMail from "@sendgrid/mail";
import { NextRequest } from "next/server";

sgMail.setApiKey(process.env.SENDGRID_API_KEY!);

export async function POST(req: NextRequest) {
  try {
    const { email, datasetId } = await req.json();

    const baseUrl = process.env.NEXT_PUBLIC_APP_URL || "http://localhost:3000";
    const link = `${baseUrl}/viewer/${datasetId}`;

    const msg = {
      to: email,
      from: "noreply@merfisheyes.com",
      subject: "Your Dataset is Ready - MERFISHeyes",
      text: `Your dataset has been processed and is ready to view.\n\nDataset ID: ${datasetId}\n\nView your dataset here: ${link}`,
      html: `
        <div style="font-family: Arial, sans-serif; max-width: 600px; margin: 0 auto;">
          <h2 style="color: #333;">Your Dataset is Ready!</h2>
          <p>Your dataset has been processed and is ready to view.</p>
          <p><strong>Dataset ID:</strong> ${datasetId}</p>
          <p style="margin: 30px 0;">
            <a href="${link}" style="background-color: #0070f3; color: white; padding: 12px 24px; text-decoration: none; border-radius: 5px; display: inline-block;">
              View Dataset
            </a>
          </p>
          <p style="color: #666; font-size: 14px;">
            Or copy this link: <a href="${link}">${link}</a>
          </p>
        </div>
      `,
    };

    await sgMail.send(msg);

    return Response.json({ success: true });
  } catch (error) {
    console.error("SendGrid Error:", error);

    return Response.json(
      {
        error: error instanceof Error ? error.message : "Failed to send email",
      },
      { status: 500 },
    );
  }
}
