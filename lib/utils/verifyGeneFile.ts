import { ungzip } from "pako";

/**
 * Verify a gene .bin.gz file by decompressing and inspecting contents
 * Use this for debugging/testing downloaded gene files
 */
export async function verifyGeneFile(file: File): Promise<{
  isValid: boolean;
  byteLength: number;
  coordinateCount: number;
  firstFewCoordinates: number[];
  stats: {
    min: number;
    max: number;
    avgX: number;
    avgY: number;
    avgZ?: number;
  };
}> {
  try {
    // Read file as array buffer
    const compressedBuffer = await file.arrayBuffer();

    // Decompress with pako
    const decompressed = ungzip(new Uint8Array(compressedBuffer));

    // Convert to Float32Array
    const float32Array = new Float32Array(decompressed.buffer);

    // Calculate stats
    const coordinateCount = float32Array.length;
    const moleculeCount = coordinateCount / 3; // Assuming 3D coordinates

    // Get first few coordinates for inspection
    const firstFewCoordinates = Array.from(float32Array.slice(0, 15)); // First 5 molecules (x,y,z each)

    // Calculate min/max/avg
    let minVal = Infinity;
    let maxVal = -Infinity;
    let sumX = 0,
      sumY = 0,
      sumZ = 0;

    for (let i = 0; i < coordinateCount; i += 3) {
      const x = float32Array[i];
      const y = float32Array[i + 1];
      const z = float32Array[i + 2];

      minVal = Math.min(minVal, x, y, z);
      maxVal = Math.max(maxVal, x, y, z);

      sumX += x;
      sumY += y;
      sumZ += z;
    }

    return {
      isValid: true,
      byteLength: decompressed.byteLength,
      coordinateCount,
      firstFewCoordinates,
      stats: {
        min: minVal,
        max: maxVal,
        avgX: sumX / moleculeCount,
        avgY: sumY / moleculeCount,
        avgZ: sumZ / moleculeCount,
      },
    };
  } catch (error) {
    console.error("Failed to verify gene file:", error);

    return {
      isValid: false,
      byteLength: 0,
      coordinateCount: 0,
      firstFewCoordinates: [],
      stats: { min: 0, max: 0, avgX: 0, avgY: 0 },
    };
  }
}

/**
 * Quick test function you can call from browser console
 *
 * Usage in browser console:
 * 1. Upload a .bin.gz file via file input
 * 2. const input = document.createElement('input'); input.type = 'file'; input.onchange = async (e) => { const file = e.target.files[0]; const result = await verifyGeneFile(file); console.log(result); }; input.click();
 */
export function createGeneFileVerifier() {
  const input = document.createElement("input");

  input.type = "file";
  input.accept = ".bin.gz,.gz";

  input.onchange = async (e: Event) => {
    const target = e.target as HTMLInputElement;
    const file = target.files?.[0];

    if (!file) {
      console.error("No file selected");

      return;
    }

    console.log(`\n📁 Verifying file: ${file.name} (${file.size} bytes)`);
    const result = await verifyGeneFile(file);

    if (result.isValid) {
      console.log("✅ File is valid!");
      console.log(`📊 Decompressed size: ${result.byteLength} bytes`);
      console.log(`🧬 Total coordinates: ${result.coordinateCount}`);
      console.log(`🧬 Total molecules: ${result.coordinateCount / 3}`);
      console.log(`📍 First few coordinates:`, result.firstFewCoordinates);
      console.log(`📈 Stats:`, result.stats);
      console.log(
        `   - Range: [${result.stats.min.toFixed(4)}, ${result.stats.max.toFixed(4)}]`,
      );
      console.log(`   - Should be normalized to [-1, 1] range`);
    } else {
      console.error("❌ File is invalid or corrupt");
    }
  };

  input.click();
}

// Make it available globally for testing
if (typeof window !== "undefined") {
  (window as any).verifyGeneFile = verifyGeneFile;
  (window as any).testGeneFile = createGeneFileVerifier;
}
