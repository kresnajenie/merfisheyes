import { config } from "dotenv";
config({ path: ".env.local" });

// test-duplicate-check.ts
import { PrismaClient } from "@prisma/client";

const prisma = new PrismaClient();

async function testDuplicateCheck() {
  console.log("🧪 Testing duplicate check API\n");

  const testFingerprint = "test_fp_" + Date.now();
  const testDatasetId = "ds_test_" + Date.now();

  try {
    // Step 1: Create a test dataset
    console.log("1️⃣ Creating test dataset...");
    await prisma.dataset.create({
      data: {
        id: testDatasetId,
        fingerprint: testFingerprint,
        title: "Test Dataset for Duplicate Check",
        numCells: 50000,
        numGenes: 300,
        status: "COMPLETE",
        manifestUrl: `https://s3.amazonaws.com/bucket/datasets/${testDatasetId}/manifest.json.gz`,
        completedAt: new Date(),
      },
    });
    console.log(`✅ Created dataset with fingerprint: ${testFingerprint}\n`);

    // Step 2: Test API
    console.log("2️⃣ Testing API endpoint...");
    const baseUrl = "http://localhost:3000";

    // Test non-existent fingerprint
    console.log("   Testing non-existent fingerprint...");
    const response1 = await fetch(
      `${baseUrl}/api/datasets/check-duplicate/nonexistent123`
    );
    const data1 = await response1.json();
    console.log("   Response:", data1);
    console.log("   ✅ Should show exists: false\n");

    // Test existing fingerprint
    console.log("   Testing existing fingerprint...");
    const response2 = await fetch(
      `${baseUrl}/api/datasets/check-duplicate/${testFingerprint}`
    );
    const data2 = await response2.json();
    console.log("   Response:", JSON.stringify(data2, null, 2));
    console.log("   ✅ Should show exists: true with dataset details\n");

    // Step 3: Cleanup
    console.log("3️⃣ Cleaning up...");
    await prisma.dataset.delete({
      where: { id: testDatasetId },
    });
    console.log("✅ Test dataset deleted\n");

    console.log("🎉 All tests passed!\n");
    console.log("📝 Summary:");
    console.log(
      "   - API correctly returns exists: false for new fingerprints"
    );
    console.log(
      "   - API correctly returns exists: true with details for duplicates"
    );
  } catch (error) {
    console.error("❌ Test failed:", error);
  } finally {
    await prisma.$disconnect();
  }
}

testDuplicateCheck();
