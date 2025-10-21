import sgMail from "@sendgrid/mail";
import { NextRequest, NextResponse } from "next/server";
import { prisma } from "@/lib/prisma";

sgMail.setApiKey(process.env.SENDGRID_API_KEY!);

export async function POST(request: NextRequest) {
  try {
    const { email, datasetId } = await request.json();

    if (!email || !datasetId) {
      return NextResponse.json(
        { error: "Email and datasetId are required" },
        { status: 400 }
      );
    }

    // Get dataset info
    const dataset = await prisma.dataset.findUnique({
      where: { id: datasetId },
      select: {
        title: true,
        numCells: true, // molecule count
        numGenes: true,
      },
    });

    if (!dataset) {
      return NextResponse.json(
        { error: "Dataset not found" },
        { status: 404 }
      );
    }

    const baseUrl = process.env.NEXT_PUBLIC_BASE_URL || "http://localhost:3000";
    const datasetUrl = `${baseUrl}/sm-viewer/${datasetId}`;

    // Send email using SendGrid
    const msg = {
      to: email,
      from: process.env.SENDGRID_FROM_EMAIL || "noreply@merfisheyes.com",
      subject: "Your Single Molecule Dataset is Ready!",
      html: `
        <!DOCTYPE html>
        <html>
          <head>
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Single Molecule Dataset Ready</title>
          </head>
          <body style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; line-height: 1.6; color: #333; max-width: 600px; margin: 0 auto; padding: 20px;">
            <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 30px; border-radius: 10px 10px 0 0; text-align: center;">
              <h1 style="color: white; margin: 0; font-size: 28px;">MERFISH Eyes</h1>
              <p style="color: rgba(255,255,255,0.9); margin: 10px 0 0 0; font-size: 16px;">Single Molecule Viewer</p>
            </div>

            <div style="background: #ffffff; padding: 30px; border: 1px solid #e0e0e0; border-top: none; border-radius: 0 0 10px 10px;">
              <h2 style="color: #667eea; margin-top: 0;">Your Dataset is Ready!</h2>

              <p>Your single molecule dataset has been successfully processed and is now ready to explore.</p>

              <div style="background: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0;">
                <h3 style="margin-top: 0; color: #333; font-size: 16px;">Dataset Information</h3>
                <table style="width: 100%; border-collapse: collapse;">
                  <tr>
                    <td style="padding: 8px 0; color: #666; font-weight: 500;">Name:</td>
                    <td style="padding: 8px 0; color: #333;">${dataset.title || "Untitled Dataset"}</td>
                  </tr>
                  <tr>
                    <td style="padding: 8px 0; color: #666; font-weight: 500;">Total Molecules:</td>
                    <td style="padding: 8px 0; color: #333;">${dataset.numCells.toLocaleString()}</td>
                  </tr>
                  <tr>
                    <td style="padding: 8px 0; color: #666; font-weight: 500;">Unique Genes:</td>
                    <td style="padding: 8px 0; color: #333;">${dataset.numGenes.toLocaleString()}</td>
                  </tr>
                  <tr>
                    <td style="padding: 8px 0; color: #666; font-weight: 500;">Dataset ID:</td>
                    <td style="padding: 8px 0; color: #333; font-family: monospace; font-size: 12px;">${datasetId}</td>
                  </tr>
                </table>
              </div>

              <div style="text-align: center; margin: 30px 0;">
                <a href="${datasetUrl}" style="display: inline-block; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; text-decoration: none; padding: 14px 32px; border-radius: 6px; font-weight: 600; font-size: 16px; box-shadow: 0 4px 6px rgba(102, 126, 234, 0.3);">
                  View Dataset
                </a>
              </div>

              <div style="background: #fff3cd; border: 1px solid #ffeaa7; padding: 15px; border-radius: 6px; margin: 20px 0;">
                <p style="margin: 0; color: #856404; font-size: 14px;">
                  <strong>📋 Note:</strong> Save this link to access your dataset anytime:
                </p>
                <p style="margin: 10px 0 0 0; word-break: break-all;">
                  <a href="${datasetUrl}" style="color: #667eea; text-decoration: none; font-family: monospace; font-size: 12px;">${datasetUrl}</a>
                </p>
              </div>

              <hr style="border: none; border-top: 1px solid #e0e0e0; margin: 30px 0;">

              <p style="color: #666; font-size: 14px; margin-bottom: 0;">
                Need help? Visit our <a href="${baseUrl}" style="color: #667eea; text-decoration: none;">documentation</a> or contact support.
              </p>
            </div>

            <div style="text-align: center; padding: 20px; color: #999; font-size: 12px;">
              <p style="margin: 0;">© ${new Date().getFullYear()} MERFISH Eyes - Single Molecule Viewer</p>
              <p style="margin: 10px 0 0 0;">Powered by spatial transcriptomics visualization technology</p>
            </div>
          </body>
        </html>
      `,
    };

    await sgMail.send(msg);

    return NextResponse.json({ success: true });
  } catch (error: any) {
    console.error("SendGrid error:", error);
    return NextResponse.json(
      {
        error: "Failed to send email",
        message: error.message,
      },
      { status: 500 }
    );
  }
}
