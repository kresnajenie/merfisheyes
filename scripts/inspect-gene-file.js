#!/usr/bin/env node

/**
 * Inspect a .bin.gz gene file and display coordinates
 *
 * Usage:
 *   node scripts/inspect-gene-file.js path/to/LEF1.bin.gz
 *   node scripts/inspect-gene-file.js path/to/LEF1.bin.gz --limit 10
 *   node scripts/inspect-gene-file.js path/to/LEF1.bin.gz --stats
 */

const fs = require('fs');
const zlib = require('zlib');
const path = require('path');

// Get filename from command line
const args = process.argv.slice(2);
const filename = args.find(arg => !arg.startsWith('--'));
const showStats = args.includes('--stats');
const limitArg = args.find(arg => arg.startsWith('--limit'));
const limit = limitArg ? parseInt(limitArg.split('=')[1]) : null;

if (!filename) {
  console.error('Usage: node inspect-gene-file.js <file.bin.gz> [--limit=N] [--stats]');
  process.exit(1);
}

if (!fs.existsSync(filename)) {
  console.error(`File not found: ${filename}`);
  process.exit(1);
}

// Read and decompress file
const compressed = fs.readFileSync(filename);
const decompressed = zlib.gunzipSync(compressed);

// Convert to Float32Array
const float32Array = new Float32Array(decompressed.buffer, decompressed.byteOffset, decompressed.byteLength / 4);

const totalCoords = float32Array.length;
const moleculeCount = totalCoords / 3;

console.log(`\n📁 File: ${path.basename(filename)}`);
console.log(`📊 Compressed size: ${compressed.length.toLocaleString()} bytes`);
console.log(`📊 Decompressed size: ${decompressed.length.toLocaleString()} bytes`);
console.log(`🧬 Total coordinates: ${totalCoords.toLocaleString()}`);
console.log(`🧬 Total molecules: ${moleculeCount.toLocaleString()}`);
console.log(`📦 Compression ratio: ${(compressed.length / decompressed.length * 100).toFixed(1)}%\n`);

// Calculate statistics
if (showStats) {
  let minVal = Infinity;
  let maxVal = -Infinity;
  let sumX = 0, sumY = 0, sumZ = 0;

  for (let i = 0; i < totalCoords; i += 3) {
    const x = float32Array[i];
    const y = float32Array[i + 1];
    const z = float32Array[i + 2];

    minVal = Math.min(minVal, x, y, z);
    maxVal = Math.max(maxVal, x, y, z);

    sumX += x;
    sumY += y;
    sumZ += z;
  }

  console.log('📈 Statistics:');
  console.log(`   Range: [${minVal.toFixed(6)}, ${maxVal.toFixed(6)}]`);
  console.log(`   Average X: ${(sumX / moleculeCount).toFixed(6)}`);
  console.log(`   Average Y: ${(sumY / moleculeCount).toFixed(6)}`);
  console.log(`   Average Z: ${(sumZ / moleculeCount).toFixed(6)}`);
  console.log('');
}

// Display coordinates
const displayCount = limit || Math.min(moleculeCount, 10);
console.log(`📍 First ${displayCount} molecules (x, y, z):\n`);

for (let i = 0; i < displayCount * 3; i += 3) {
  const x = float32Array[i];
  const y = float32Array[i + 1];
  const z = float32Array[i + 2];

  const moleculeNum = (i / 3) + 1;
  console.log(`  ${moleculeNum.toString().padStart(5)}: [${x.toFixed(6)}, ${y.toFixed(6)}, ${z.toFixed(6)}]`);
}

if (moleculeCount > displayCount) {
  console.log(`\n  ... and ${(moleculeCount - displayCount).toLocaleString()} more molecules`);
}

console.log('');
