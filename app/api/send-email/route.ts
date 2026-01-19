import sgMail from "@sendgrid/mail";
import { NextRequest } from "next/server";

sgMail.setApiKey(process.env.SENDGRID_API_KEY!);

export async function POST(req: NextRequest) {
  try {
    const { email, datasetId, datasetName, metadata } = await req.json();

    const baseUrl = process.env.NEXT_PUBLIC_APP_URL || "http://localhost:3000";
    const link = `${baseUrl}/viewer/${datasetId}`;

    // Use dataset name in subject, fallback to generic message
    const subject = datasetName
      ? `${datasetName} - Dataset Ready - MERFISHeyes`
      : "Your Dataset is Ready - MERFISHeyes";

    const msg = {
      to: email,
      from: "noreply@merfisheyes.com",
      subject: subject,
      text: `Your dataset has been processed and is ready to view.\n\nDataset: ${datasetName || "Untitled"}\nCells: ${metadata?.numCells?.toLocaleString() || "N/A"}\nGenes: ${metadata?.numGenes?.toLocaleString() || "N/A"}\nPlatform: ${metadata?.platform || "N/A"}\n\nView your dataset here: ${link}`,
      html: `
        <!DOCTYPE html>
        <html>
          <head>
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Dataset Ready</title>
          </head>
          <body style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; line-height: 1.6; color: #333; max-width: 600px; margin: 0 auto; padding: 20px;">
            <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 30px; border-radius: 10px 10px 0 0; text-align: center;">
              <h1 style="color: white; margin: 0; font-size: 28px;">MERFISH Eyes</h1>
              <p style="color: rgba(255,255,255,0.9); margin: 10px 0 0 0; font-size: 16px;">Single Cell Viewer</p>
            </div>

            <div style="background: #ffffff; padding: 30px; border: 1px solid #e0e0e0; border-top: none; border-radius: 0 0 10px 10px;">
              <h2 style="color: #667eea; margin-top: 0;">Your Dataset is Ready!</h2>

              <p>Your single cell dataset has been successfully processed and is now ready to explore.</p>

              <div style="background: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0;">
                <h3 style="margin-top: 0; color: #333; font-size: 16px;">Dataset Information</h3>
                <table style="width: 100%; border-collapse: collapse;">
                  <tr>
                    <td style="padding: 8px 0; color: #666; font-weight: 500;">Name:</td>
                    <td style="padding: 8px 0; color: #333;">${datasetName || "Untitled Dataset"}</td>
                  </tr>
                  <tr>
                    <td style="padding: 8px 0; color: #666; font-weight: 500;">Total Cells:</td>
                    <td style="padding: 8px 0; color: #333;">${metadata?.numCells?.toLocaleString() || "N/A"}</td>
                  </tr>
                  <tr>
                    <td style="padding: 8px 0; color: #666; font-weight: 500;">Unique Genes:</td>
                    <td style="padding: 8px 0; color: #333;">${metadata?.numGenes?.toLocaleString() || "N/A"}</td>
                  </tr>
                  <tr>
                    <td style="padding: 8px 0; color: #666; font-weight: 500;">Platform:</td>
                    <td style="padding: 8px 0; color: #333; text-transform: uppercase;">${metadata?.platform || "N/A"}</td>
                  </tr>
                  ${
                    metadata?.clusterCount
                      ? `
                  <tr>
                    <td style="padding: 8px 0; color: #666; font-weight: 500;">Cluster Columns:</td>
                    <td style="padding: 8px 0; color: #333;">${metadata.clusterCount}</td>
                  </tr>
                  `
                      : ""
                  }
                  <tr>
                    <td style="padding: 8px 0; color: #666; font-weight: 500;">Dataset ID:</td>
                    <td style="padding: 8px 0; color: #333; font-family: monospace; font-size: 12px;">${datasetId}</td>
                  </tr>
                </table>
              </div>

              <div style="text-align: center; margin: 30px 0;">
                <a href="${link}" style="display: inline-block; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; text-decoration: none; padding: 14px 32px; border-radius: 6px; font-weight: 600; font-size: 16px; box-shadow: 0 4px 6px rgba(102, 126, 234, 0.3);">
                  View Dataset
                </a>
              </div>

              <div style="background: #fff3cd; border: 1px solid #ffeaa7; padding: 15px; border-radius: 6px; margin: 20px 0;">
                <p style="margin: 0; color: #856404; font-size: 14px;">
                  <strong>📋 Note:</strong> Save this link to access your dataset anytime:
                </p>
                <p style="margin: 10px 0 0 0; word-break: break-all;">
                  <a href="${link}" style="color: #667eea; text-decoration: none; font-family: monospace; font-size: 12px;">${link}</a>
                </p>
              </div>

              <hr style="border: none; border-top: 1px solid #e0e0e0; margin: 30px 0;">

              <p style="color: #666; font-size: 14px; margin-bottom: 0;">
                Need help? Visit our <a href="${baseUrl}" style="color: #667eea; text-decoration: none;">documentation</a> or contact support.
              </p>
            </div>

            <div style="text-align: center; padding: 20px; color: #999; font-size: 12px;">
              <p style="margin: 0;">© ${new Date().getFullYear()} MERFISH Eyes - Single Cell Viewer</p>
              <p style="margin: 10px 0 0 0;">Powered by spatial transcriptomics visualization technology</p>
            </div>
          </body>
        </html>
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
