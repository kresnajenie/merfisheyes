import { NextRequest } from "next/server";

import { sendDatasetReadyEmail } from "@/lib/ses";
import { prisma } from "@/lib/prisma";

export async function POST(req: NextRequest) {
  try {
    const { email, datasetId, datasetName, metadata } = await req.json();

    // The link opens the viewer of what the dataset is
    const dataset = await prisma.dataset.findUnique({
      where: { id: String(datasetId) },
      select: { datasetType: true },
    });

    await sendDatasetReadyEmail({
      email,
      datasetId,
      datasetType: dataset?.datasetType,
      datasetName,
      metadata,
    });

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
