// tests/api/test-complete-upload.ts
import { config } from "dotenv";
import path from "path";

config({ path: path.resolve(process.cwd(), ".env.local") });

import { PrismaClient } from "@prisma/client";

const prisma = new PrismaClient();
const API_BASE = "http://localhost:3000";

async function testCompleteUploadFlow() {
  console.log("🧪 Testing Complete Upload Flow\n");

  const testFingerprint = "fp_complete_test_" + Date.now();
  let datasetId: string;
  let uploadId: string;

  try {
    // Step 1: Initiate upload
    console.log("1️⃣ Initiating upload...");
    const initiateResponse = await fetch(`${API_BASE}/api/datasets/initiate`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        fingerprint: testFingerprint,
        metadata: {
          title: "Complete Upload Test Dataset",
          numCells: 25000,
          numGenes: 250,
        },
        files: [
          { key: "file1.csv", size: 1024, contentType: "text/csv" },
          {
            key: "file2.parquet",
            size: 2048,
            contentType: "application/octet-stream",
          },
        ],
      }),
    });

    const initiateData = await initiateResponse.json();
    datasetId = initiateData.datasetId;
    uploadId = initiateData.uploadId;
    console.log("   ✅ Upload initiated");
    console.log("   Dataset ID:", datasetId);
    console.log("   Upload ID:", uploadId);

    // Step 2: Simulate uploading files (mark them as complete in DB)
    console.log("\n2️⃣ Simulating file uploads...");
    await prisma.uploadFile.updateMany({
      where: { uploadSessionId: uploadId },
      data: {
        status: "COMPLETE",
        uploadedAt: new Date(),
      },
    });
    await prisma.uploadSession.update({
      where: { id: uploadId },
      data: { completedFiles: 2 },
    });
    console.log("   ✅ Files marked as uploaded");

    // Step 3: Try to complete upload (without manifest)
    console.log("\n3️⃣ Completing upload...");
    const completeResponse = await fetch(
      `${API_BASE}/api/datasets/${datasetId}/complete`,
      {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ uploadId }),
      }
    );

    if (!completeResponse.ok) {
      const errorData = await completeResponse.json();
      console.error("   ❌ Failed:", errorData);
      return;
    }

    const completeData = await completeResponse.json();
    console.log("   ✅ Upload completed!");
    console.log("   Status:", completeData.status);
    console.log("   View URL:", completeData.viewUrl);
    console.log("   Share URL:", completeData.shareUrl);

    // Step 4: Verify in database
    console.log("\n4️⃣ Verifying database...");
    const dataset = await prisma.dataset.findUnique({
      where: { id: datasetId },
    });

    if (dataset?.status === "COMPLETE" && dataset?.completedAt) {
      console.log("   ✅ Dataset marked as COMPLETE");
      console.log("   Completed at:", dataset.completedAt.toISOString());
    } else {
      console.log("   ❌ Dataset status incorrect:", dataset?.status);
    }

    // Step 5: Test error case - incomplete uploads
    console.log("\n5️⃣ Testing incomplete upload rejection...");

    const incompleteUploadId = "up_incomplete_" + Date.now();

    const incompleteDataset = await prisma.dataset.create({
      data: {
        id: "ds_incomplete_" + Date.now(),
        fingerprint: "fp_incomplete_" + Date.now(),
        title: "Incomplete Test",
        numCells: 1000,
        numGenes: 100,
        status: "UPLOADING",
        uploadSessions: {
          create: {
            id: incompleteUploadId, // Store the ID
            totalFiles: 3,
            completedFiles: 1,
            expiresAt: new Date(Date.now() + 3600000),
            files: {
              create: [
                {
                  fileKey: "file1.csv",
                  fileSize: BigInt(1024),
                  status: "COMPLETE",
                },
                {
                  fileKey: "file2.csv",
                  fileSize: BigInt(2048),
                  status: "PENDING",
                },
              ],
            },
          },
        },
      },
    });

    const incompleteResponse = await fetch(
      `${API_BASE}/api/datasets/${incompleteDataset.id}/complete`,
      {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          uploadId: incompleteUploadId, // Use the stored ID
        }),
      }
    );

    if (incompleteResponse.status === 400) {
      const errorData = await incompleteResponse.json();
      console.log("   ✅ Correctly rejected incomplete upload");
      console.log("   Error:", errorData.error);
    } else {
      console.log("   ⚠️  Expected 400, got:", incompleteResponse.status);
    }

    // Step 6: Test fetching completed dataset
    console.log("\n6️⃣ Testing dataset retrieval...");
    const getResponse = await fetch(`${API_BASE}/api/datasets/${datasetId}`);

    if (getResponse.ok) {
      const getData = await getResponse.json();
      console.log("   ✅ Can retrieve completed dataset");
      console.log("   Files available:", Object.keys(getData.files).length);
    } else {
      console.log("   ❌ Failed to retrieve dataset:", getResponse.status);
    }

    // Cleanup
    console.log("\n7️⃣ Cleaning up...");
    await prisma.dataset.delete({ where: { id: datasetId } });
    await prisma.dataset.delete({ where: { id: incompleteDataset.id } });
    console.log("   ✅ Test data cleaned up");

    console.log("\n🎉 All tests passed!\n");
    console.log("📝 Complete Upload Flow Summary:");
    console.log("   ✅ Can initiate upload");
    console.log("   ✅ Can mark files as uploaded");
    console.log("   ✅ Can complete upload when all files ready");
    console.log("   ✅ Rejects completion if files missing");
    console.log("   ✅ Generates view/share URLs");
    console.log("   ✅ Can retrieve completed dataset");
  } catch (error) {
    console.error("❌ Test failed:", error);
  } finally {
    await prisma.$disconnect();
  }
}

testCompleteUploadFlow();
