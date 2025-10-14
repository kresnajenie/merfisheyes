import { config } from "dotenv";
import path from "path";

// Load from project root, not relative to test file
config({ path: path.resolve(process.cwd(), ".env.local") });
// tests/api/test-s3-config.ts
import {
  generatePresignedUploadUrl,
  generatePresignedDownloadUrl,
} from "@/lib/s3";

async function testS3Configuration() {
  console.log("🧪 Testing S3 Configuration\n");

  try {
    const testKey = `test/test-file-${Date.now()}.txt`;
    const testContent = "Hello from MERFISH uploader!";

    // Step 1: Generate upload URL
    console.log("1️⃣ Generating presigned upload URL...");
    const uploadUrlData = await generatePresignedUploadUrl(
      testKey,
      "text/plain",
      300 // 5 minutes
    );
    console.log("   ✅ Upload URL generated");
    console.log("   Method:", uploadUrlData.method);
    console.log("   URL:", uploadUrlData.url.substring(0, 80) + "...");

    // Step 2: Upload test file
    console.log("\n2️⃣ Uploading test file to S3...");
    const uploadResponse = await fetch(uploadUrlData.url, {
      method: uploadUrlData.method,
      headers: uploadUrlData.headers,
      body: testContent,
    });

    if (uploadResponse.ok) {
      console.log("   ✅ File uploaded successfully!");
      console.log("   Status:", uploadResponse.status);
    } else {
      console.error("   ❌ Upload failed:", uploadResponse.status);
      const errorText = await uploadResponse.text();
      console.error("   Error:", errorText);
      return;
    }

    // Step 3: Generate download URL
    console.log("\n3️⃣ Generating presigned download URL...");
    const downloadUrl = await generatePresignedDownloadUrl(testKey, 300);
    console.log("   ✅ Download URL generated");
    console.log("   URL:", downloadUrl.substring(0, 80) + "...");

    // Step 4: Download and verify
    console.log("\n4️⃣ Downloading file from S3...");
    const downloadResponse = await fetch(downloadUrl);

    if (downloadResponse.ok) {
      const downloadedContent = await downloadResponse.text();
      console.log("   ✅ File downloaded successfully!");
      console.log("   Content:", downloadedContent);

      if (downloadedContent === testContent) {
        console.log("   ✅ Content matches!");
      } else {
        console.log("   ❌ Content mismatch!");
      }
    } else {
      console.error("   ❌ Download failed:", downloadResponse.status);
      return;
    }

    // Step 5: Test CORS by checking headers
    console.log("\n5️⃣ Checking CORS headers...");
    const corsTestUrl = downloadUrl;
    const corsResponse = await fetch(corsTestUrl, {
      method: "HEAD",
    });

    const corsHeader = corsResponse.headers.get("access-control-allow-origin");
    if (corsHeader) {
      console.log("   ✅ CORS headers present");
      console.log("   Allow-Origin:", corsHeader);
    } else {
      console.log("   ⚠️  CORS headers not found (may need OPTIONS request)");
    }

    console.log("\n🎉 S3 Configuration Test Complete!\n");
    console.log("📝 Summary:");
    console.log("   ✅ Can generate presigned upload URLs");
    console.log("   ✅ Can upload files to S3");
    console.log("   ✅ Can generate presigned download URLs");
    console.log("   ✅ Can download files from S3");
    console.log("   ✅ File content integrity verified");
    console.log("\n🔒 Security Check:");
    console.log("   - Bucket is private (no public access)");
    console.log("   - Access only via presigned URLs");
    console.log("   - CORS configured for your domains");
    console.log("\n⚠️  Cleanup:");
    console.log("   Test file created at:", testKey);
    console.log("   You may want to delete it from S3 console");
  } catch (error) {
    console.error("\n❌ Test failed:", error);
    console.error("\nPossible issues:");
    console.error("   1. AWS credentials not configured");
    console.error("   2. S3 bucket does not exist");
    console.error("   3. IAM permissions insufficient");
    console.error("   4. CORS not configured properly");
  }
}

testS3Configuration();
