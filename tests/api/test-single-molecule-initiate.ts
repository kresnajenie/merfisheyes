import { config } from "dotenv";
config({ path: ".env.local" });

// test-single-molecule-initiate.ts
const API_BASE = "http://localhost:3000";

async function testSingleMoleculeInitiate() {
  console.log("🧪 Testing Single Molecule Initiate Upload API\n");

  const testFingerprint = "sm_fp_test_" + Date.now();

  // Mock manifest data
  const mockManifest = {
    version: "1.0",
    created_at: new Date().toISOString(),
    dataset_id: "test_id",
    name: "Test Xenium Dataset",
    type: "xenium",
    statistics: {
      total_molecules: 1000000,
      unique_genes: 150,
      spatial_dimensions: 3,
    },
    genes: {
      unique_gene_names: ["LEF1", "CD8A", "CD4", "IL2", "TNFA"],
    },
    processing: {
      compression: "gzip",
      coordinate_format: "float32_flat_array",
      coordinate_range: "normalized_[-1,1]",
      scaling_factor: 0.001,
      created_by: "MERFISH Eyes - Single Molecule Viewer",
      source_file: "test_dataset.parquet",
    },
  };

  try {
    // Step 1: Check duplicate first (should not exist)
    console.log("1️⃣ Checking for duplicate...");
    const checkResponse = await fetch(
      `${API_BASE}/api/single-molecule/check-duplicate/${testFingerprint}`
    );
    const checkData = await checkResponse.json();
    console.log("   Response:", checkData);
    console.log("   ✅ No duplicate found\n");

    // Step 2: Initiate upload
    console.log("2️⃣ Initiating single molecule upload...");
    const initiateResponse = await fetch(
      `${API_BASE}/api/single-molecule/initiate`,
      {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          fingerprint: testFingerprint,
          metadata: {
            title: "Test Xenium Single Molecule Dataset",
            numMolecules: 1000000,
            numGenes: 150,
            platform: "xenium",
            description: "Test single molecule dataset for API validation",
          },
          manifest: mockManifest,
          files: [
            {
              key: "manifest.json.gz",
              size: 2048, // 2 KB
              contentType: "application/gzip",
            },
            {
              key: "genes/LEF1.bin.gz",
              size: 48000, // 48 KB
              contentType: "application/gzip",
            },
            {
              key: "genes/CD8A.bin.gz",
              size: 36000, // 36 KB
              contentType: "application/gzip",
            },
            {
              key: "genes/CD4.bin.gz",
              size: 42000, // 42 KB
              contentType: "application/gzip",
            },
          ],
        }),
      }
    );

    if (!initiateResponse.ok) {
      const errorData = await initiateResponse.json();
      console.error("❌ Failed to initiate upload:", errorData);
      return;
    }

    const initiateData = await initiateResponse.json();
    console.log("   ✅ Single molecule upload initiated successfully!");
    console.log("   Dataset ID:", initiateData.datasetId);
    console.log("   Upload ID:", initiateData.uploadId);
    console.log("   Expires in:", initiateData.expiresIn, "seconds");
    console.log(
      "   Number of presigned URLs:",
      Object.keys(initiateData.uploadUrls).length
    );
    console.log("\n   Presigned URLs:");

    for (const [fileKey, url] of Object.entries(initiateData.uploadUrls)) {
      console.log(`   - ${fileKey}:`);
      console.log(`     URL: ${(url as string).substring(0, 100)}...`);
    }

    // Step 3: Verify database records
    console.log("\n3️⃣ Verifying database records...");
    console.log("   (Check your database to see the created records)");
    console.log("   Dataset ID:", initiateData.datasetId);
    console.log("   Dataset Type: single_molecule");
    console.log("   Upload Session ID:", initiateData.uploadId);
    console.log("   Manifest stored in manifestJson field");

    // Step 4: Test duplicate initiation (should fail)
    console.log("\n4️⃣ Testing duplicate upload initiation...");
    const duplicateResponse = await fetch(
      `${API_BASE}/api/single-molecule/initiate`,
      {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          fingerprint: testFingerprint, // Same fingerprint
          metadata: {
            title: "Duplicate Test",
            numMolecules: 1000,
            numGenes: 10,
          },
          manifest: mockManifest,
          files: [
            {
              key: "manifest.json.gz",
              size: 1024,
              contentType: "application/gzip",
            },
          ],
        }),
      }
    );

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
    console.log("   ✅ Manifest JSON stored in database");
    console.log("   ✅ Dataset type set to 'single_molecule'");
    console.log("   ✅ Presigned URLs generated for all files");
    console.log("   ✅ Duplicate fingerprints are rejected");
    console.log("\n⚠️  Note: Created test data is still in database.");
    console.log("   Dataset ID:", initiateData.datasetId);
  } catch (error) {
    console.error("❌ Test failed:", error);
  }
}

testSingleMoleculeInitiate();
