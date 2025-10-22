// src/utils/gzip.ts
import { ungzip } from "pako";

/**
 * Return file contents as UTF-8 text.
 * Supports plain text and .gz transparently.
 */
export async function fileToTextMaybeGz(file: File): Promise<string> {
  const name = (file.name || "").toLowerCase();
  const isGz = name.endsWith(".gz");

  // Fast path: plain file
  if (!isGz) return file.text();

  // Prefer built-in DecompressionStream when available (Chrome/Edge, Safari 17+)
  if ("DecompressionStream" in window) {
    const ds = new DecompressionStream("gzip");
    const decompressed = file.stream().pipeThrough(ds);
    const resp = new Response(decompressed);

    return await resp.text();
  }

  // Fallback: pako
  const buf = await file.arrayBuffer();
  const text = ungzip(new Uint8Array(buf), { to: "string" });

  return text;
}
