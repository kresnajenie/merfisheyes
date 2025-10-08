// test-prisma.ts
import { PrismaClient } from "@prisma/client";

const prisma = new PrismaClient({
  log: ["query", "info", "warn", "error"],
});

async function main() {
  console.log("🔍 Testing Prisma connection...\n");

  try {
    // Test 1: Create a dataset
    console.log("1️⃣ Creating test dataset...");
    const dataset = await prisma.dataset.create({
      data: {
        id: "ds_test_" + Date.now(),
        fingerprint: "fp_" + Math.random().toString(36).substring(7),
        title: "Test MERFISH Dataset",
        numCells: 100000,
        numGenes: 500,
        status: "UPLOADING",
        manifestUrl:
          "https://s3.amazonaws.com/bucket/datasets/ds_test/manifest.json.gz",
      },
    });
    console.log("✅ Created dataset:", {
      id: dataset.id,
      title: dataset.title,
      numCells: dataset.numCells,
      status: dataset.status,
    });

    // Test 2: Create upload session
    console.log("\n2️⃣ Creating upload session...");
    const uploadSession = await prisma.uploadSession.create({
      data: {
        id: "up_test_" + Date.now(),
        datasetId: dataset.id,
        totalFiles: 3,
        completedFiles: 0,
        expiresAt: new Date(Date.now() + 3600000), // 1 hour
      },
    });
    console.log("✅ Created upload session:", {
      id: uploadSession.id,
      totalFiles: uploadSession.totalFiles,
      expiresAt: uploadSession.expiresAt,
    });

    // Test 3: Create upload files
    console.log("\n3️⃣ Creating upload files...");
    const files = await prisma.uploadFile.createMany({
      data: [
        {
          uploadSessionId: uploadSession.id,
          fileKey: "cell_boundaries.parquet",
          fileSize: BigInt(1024 * 1024 * 50), // 50MB
          status: "PENDING",
        },
        {
          uploadSessionId: uploadSession.id,
          fileKey: "cell_metadata.csv",
          fileSize: BigInt(1024 * 100), // 100KB
          status: "PENDING",
        },
        {
          uploadSessionId: uploadSession.id,
          fileKey: "detected_transcripts.csv",
          fileSize: BigInt(1024 * 1024 * 200), // 200MB
          status: "PENDING",
        },
      ],
    });
    console.log("✅ Created", files.count, "upload files");

    // Test 4: Query with relations
    console.log("\n4️⃣ Querying dataset with relations...");
    const datasetWithRelations = await prisma.dataset.findUnique({
      where: { id: dataset.id },
      include: {
        uploadSessions: {
          include: {
            files: true,
          },
        },
      },
    });
    console.log("✅ Dataset with relations:", {
      dataset: datasetWithRelations?.title,
      sessions: datasetWithRelations?.uploadSessions.length,
      totalFiles: datasetWithRelations?.uploadSessions[0]?.files.length,
    });

    // Test 5: Update file status
    console.log("\n5️⃣ Updating file status...");
    const updatedFile = await prisma.uploadFile.updateMany({
      where: {
        uploadSessionId: uploadSession.id,
        fileKey: "cell_metadata.csv",
      },
      data: {
        status: "COMPLETE",
        uploadedAt: new Date(),
      },
    });
    console.log("✅ Updated", updatedFile.count, "file(s)");

    // Test 6: Check duplicate by fingerprint
    console.log("\n6️⃣ Testing duplicate check...");
    const duplicate = await prisma.dataset.findUnique({
      where: { fingerprint: dataset.fingerprint },
    });
    console.log("✅ Found duplicate:", duplicate ? "Yes" : "No");

    // Test 7: Count datasets by status
    console.log("\n7️⃣ Counting datasets by status...");
    const uploadingCount = await prisma.dataset.count({
      where: { status: "UPLOADING" },
    });
    console.log("✅ Uploading datasets:", uploadingCount);

    // Test 8: Cleanup
    console.log("\n8️⃣ Cleaning up test data...");
    await prisma.dataset.delete({
      where: { id: dataset.id },
    });
    console.log("✅ Cleanup complete (cascade deleted session and files)");

    console.log("\n🎉 All tests passed!");
  } catch (error) {
    console.error("❌ Test failed:", error);
    process.exit(1);
  }
}

main()
  .catch((e) => {
    console.error("Fatal error:", e);
    process.exit(1);
  })
  .finally(async () => {
    await prisma.$disconnect();
  });
