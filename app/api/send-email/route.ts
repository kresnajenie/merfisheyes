import { NextRequest } from "next/server";

import { sendDatasetReadyEmail } from "@/lib/ses";

export async function POST(req: NextRequest) {
  try {
    const { email, datasetId, datasetName, metadata } = await req.json();

    await sendDatasetReadyEmail({ email, datasetId, datasetName, metadata });

    return Response.json({ success: true });
  } catch (error) {
    console.error("SES Error:", error);

    return Response.json(
      {
        error: error instanceof Error ? error.message : "Failed to send email",
      },
      { status: 500 },
    );
  }
}
