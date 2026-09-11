/**
 * Per-cell segmentation surfaces for the labelled-molecule viewer.
 *
 * Two kinds share this format and this loader:
 *   cells   written by scripts/spiralia/export_meshes.py
 *   nuclei  written by scripts/spiralia/export_nuclei.py
 *
 * Both are keyed by the same `cell` obs labels, so a nucleus follows the same
 * selection as the cell that contains it.
 *
 * Vertices are in the molecules' µm frame — cells use `vertices_xyz`, not
 * `vertices_reoriented_xyz` which is a different frame; nuclei are converted
 * from segm voxels at 1.0 µm in z and 0.708 µm in x/y — so no transform is
 * applied here.
 *
 *   <index>.json   per-cell label + vertex/index slices
 *   <bin>.bin.gz   u32 version | u32 cells | u32 nVerts | u32 nIndices
 *                  Float32 vertices[nVerts * 3]
 *                  Uint32  indices[nIndices]   (local to each cell)
 */

export interface CellMesh {
  /** Matches a value in the `cell` obs column, so it can follow a selection. */
  label: string;
  cellId: number;
  positions: Float32Array;
  indices: Uint32Array;
}

interface MeshIndexEntry {
  label: string;
  cell_id: number;
  vertex_offset: number;
  vertex_count: number;
  index_offset: number;
  index_count: number;
}

async function gunzip(response: Response): Promise<ArrayBuffer> {
  const stream = response.body!.pipeThrough(new DecompressionStream("gzip"));

  return await new Response(stream).arrayBuffer();
}

export type SurfaceKind = "cells" | "nuclei";

const FILES: Record<SurfaceKind, { index: string; bin: string }> = {
  cells: { index: "meshes/index.json", bin: "meshes/cells.bin.gz" },
  nuclei: { index: "meshes/nuclei_index.json", bin: "meshes/nuclei.bin.gz" },
};

/**
 * Fetch one kind of surface for a dataset, or null when it has none.
 *
 * Returns null rather than throwing on a missing file: a dataset without
 * meshes is a normal state, not a failure.
 */
export async function loadSurfaceMeshes(
  baseUrl: string,
  kind: SurfaceKind,
): Promise<CellMesh[] | null> {
  const base = baseUrl.replace(/\/+$/, "");
  const files = FILES[kind];

  try {
    const [indexRes, binRes] = await Promise.all([
      fetch(`${base}/${files.index}`),
      fetch(`${base}/${files.bin}`),
    ]);

    if (!indexRes.ok || !binRes.ok) return null;

    const index = (await indexRes.json()) as { cells: MeshIndexEntry[] };
    const buf = await gunzip(binRes);
    const header = new Uint32Array(buf, 0, 4);
    const nVerts = header[2];
    const nIndices = header[3];

    const vertBytes = 16;
    const idxBytes = vertBytes + nVerts * 3 * 4;
    const allPositions = new Float32Array(buf, vertBytes, nVerts * 3);
    const allIndices = new Uint32Array(buf, idxBytes, nIndices);

    return index.cells.map((c) => ({
      label: c.label,
      cellId: c.cell_id,
      // slice() rather than subarray(): each geometry owns its buffer, so
      // disposing one cell can't invalidate another.
      positions: allPositions.slice(
        c.vertex_offset * 3,
        (c.vertex_offset + c.vertex_count) * 3,
      ),
      indices: allIndices.slice(c.index_offset, c.index_offset + c.index_count),
    }));
  } catch {
    return null;
  }
}
