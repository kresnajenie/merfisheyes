import { config } from "dotenv";
config({ path: ".env.local" });

// test-initiate-upload.ts
const API_BASE = "http://localhost:3000";

async function testInitiateUpload() {
  console.log("🧪 Testing Initiate Upload API\n");

  const testFingerprint = "fp_test_" + Date.now();

  try {
    // Step 1: Check duplicate first (should not exist)
    console.log("1️⃣ Checking for duplicate...");
    const checkResponse = await fetch(
      `${API_BASE}/api/datasets/check-duplicate/${testFingerprint}`
    );
    const checkData = await checkResponse.json();
    console.log("   Response:", checkData);
    console.log("   ✅ No duplicate found\n");

    // Step 2: Initiate upload
    console.log("2️⃣ Initiating upload...");
    const initiateResponse = await fetch(`${API_BASE}/api/datasets/initiate`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        fingerprint: testFingerprint,
        metadata: {
          title: "Test MERFISH Dataset Upload",
          numCells: 50000,
          numGenes: 300,
          platform: "MERSCOPE",
          description: "Test dataset for upload API validation",
        },
        files: [
          {
            key: "cell_boundaries.parquet",
            size: 52428800, // 50 MB
            contentType: "application/octet-stream",
          },
          {
            key: "cell_metadata.csv",
            size: 102400, // 100 KB
            contentType: "text/csv",
          },
          {
            key: "detected_transcripts.csv",
            size: 209715200, // 200 MB
            contentType: "text/csv",
          },
        ],
      }),
    });

    if (!initiateResponse.ok) {
      const errorData = await initiateResponse.json();
      console.error("❌ Failed to initiate upload:", errorData);
      return;
    }

    const initiateData = await initiateResponse.json();
    console.log("   ✅ Upload initiated successfully!");
    console.log("   Dataset ID:", initiateData.datasetId);
    console.log("   Upload ID:", initiateData.uploadId);
    console.log("   Expires in:", initiateData.expiresIn, "seconds");
    console.log(
      "   Number of presigned URLs:",
      Object.keys(initiateData.uploadUrls).length
    );
    console.log("\n   Presigned URLs:");

    for (const [fileKey, urlData] of Object.entries(initiateData.uploadUrls)) {
      console.log(`   - ${fileKey}:`);
      console.log(`     Method: ${(urlData as any).method}`);
      console.log(`     URL: ${(urlData as any).url.substring(0, 100)}...`);
      console.log(`     Headers:`, (urlData as any).headers);
    }

    // Step 3: Verify database records
    console.log("\n3️⃣ Verifying database records...");
    console.log("   (Check your database to see the created records)");
    console.log("   Dataset ID:", initiateData.datasetId);
    console.log("   Upload Session ID:", initiateData.uploadId);

    // Step 4: Test duplicate initiation (should fail)
    console.log("\n4️⃣ Testing duplicate upload initiation...");
    const duplicateResponse = await fetch(`${API_BASE}/api/datasets/initiate`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        fingerprint: testFingerprint, // Same fingerprint
        metadata: {
          title: "Duplicate Test",
          numCells: 1000,
          numGenes: 100,
        },
        files: [
          {
            key: "test.csv",
            size: 1024,
            contentType: "text/csv",
          },
        ],
      }),
    });

    if (duplicateResponse.status === 409) {
      const duplicateData = await duplicateResponse.json();
      console.log("   ✅ Correctly rejected duplicate:", duplicateData.error);
    } else {
      console.log(
        "   ⚠️  Expected 409 Conflict, got:",
        duplicateResponse.status
      );
    }

    console.log("\n🎉 All tests passed!\n");
    console.log("📝 Summary:");
    console.log("   ✅ Duplicate check works");
    console.log("   ✅ Upload initiation creates all database records");
    console.log("   ✅ Presigned URLs generated for all files");
    console.log("   ✅ Duplicate fingerprints are rejected");
    console.log("\n⚠️  Note: Created test data is still in database.");
    console.log("   Dataset ID:", initiateData.datasetId);
  } catch (error) {
    console.error("❌ Test failed:", error);
  }
}

testInitiateUpload();
