/**
 * Per-cell segmentation meshes for the labelled-molecule viewer.
 *
 * Written by scripts/spiralia/export_meshes.py — the two must stay in step.
 * Vertices are in the same µm frame as the molecules (the exporter uses
 * `vertices_xyz`, not `vertices_reoriented_xyz`, which is a different frame),
 * so no transform is applied here.
 *
 *   meshes/index.json    per-cell label + vertex/index slices
 *   meshes/cells.bin.gz  u32 version | u32 cells | u32 nVerts | u32 nIndices
 *                        Float32 vertices[nVerts * 3]
 *                        Uint32  indices[nIndices]   (local to each cell)
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

/**
 * Fetch the cell meshes for a dataset, or null when it has none.
 *
 * Returns null rather than throwing on a missing file: most datasets have no
 * meshes, and that is a normal state, not a failure.
 */
export async function loadCellMeshes(
  baseUrl: string,
): Promise<CellMesh[] | null> {
  const base = baseUrl.replace(/\/+$/, "");

  try {
    const [indexRes, binRes] = await Promise.all([
      fetch(`${base}/meshes/index.json`),
      fetch(`${base}/meshes/cells.bin.gz`),
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
