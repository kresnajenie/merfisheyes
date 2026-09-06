#!/usr/bin/env python3
"""Export per-cell segmentation meshes alongside a labelled-molecule dataset.

    python export_meshes.py MER6-2_E3_1 /path/to/MER6-2_E3_1_lm

Writes into the dataset folder:

    meshes/index.json     per-cell label + vertex/index slices
    meshes/cells.bin.gz   Float32 vertices, then Uint32 indices

and sets `has_meshes: true` in the manifest.

Uses `vertices_xyz`, NOT `vertices_reoriented_xyz`. Verified against
MER6-2_E3_1: per-cell centroids sit 1.5–8 µm from the molecule centroids for
every small cell and 97–100% of each cell's molecules fall inside its own mesh
bounding box, whereas the reoriented vertices are ~40 µm out. The reoriented
set is a different frame and would render plausibly but wrong.

Cell labels must match the dataset's `cell` obs column exactly, including the
"(id)" suffix applied when two cells share a name, or the mesh cannot be tied
to a selection. That is asserted rather than assumed.
"""

import argparse
import gzip
import json
import struct
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

BOGDAN = Path("/home/data/yiqun-spiralia/Sep2026/Segmentation/Bogdan")
MASK_MESHES = BOGDAN / "annotated_domain_viz_cell_meshes"
PC_MESHES = BOGDAN / "annotated_domain_viz_pointcloud_nonoverlap_cell_meshes"


def dataset_cell_labels(out: Path) -> list[str]:
    """The `cell` obs dictionary — the labels the viewer will filter on."""
    raw = gzip.open(out / "obs/cell.bin.gz", "rb").read()
    _, _, _, _, _, dict_len = struct.unpack_from("<IIBBHI", raw, 0)

    return json.loads(raw[16 : 16 + dict_len].decode())


def display_labels(embryo: str) -> dict[str, str]:
    """cell_id -> display label, mirroring export_parquet's disambiguation.

    Identities come from the pointcloud mesh file, the only variant that
    records `cell_identity`; the mask-derived file keys groups by id alone.
    """
    with h5py.File(PC_MESHES / f"{embryo}_cell_meshes.h5", "r") as f:
        ident = {cid: str(f[f"meshes/{cid}"].attrs["cell_identity"]) for cid in f["meshes"]}

    counts = Counter(ident.values())

    return {cid: (v if counts[v] == 1 else f"{v} ({cid})") for cid, v in ident.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("embryo")
    ap.add_argument("dataset_dir")
    args = ap.parse_args()

    out = Path(args.dataset_dir)
    labels = display_labels(args.embryo)
    known = set(dataset_cell_labels(out))

    verts: list[np.ndarray] = []
    idx: list[np.ndarray] = []
    cells = []
    v_off = i_off = 0
    skipped = []

    with h5py.File(MASK_MESHES / f"{args.embryo}_cell_meshes.h5", "r") as f:
        for cid in sorted(f["meshes"], key=int):
            label = labels.get(cid)

            # A mesh we can't tie to a selectable cell is worse than no mesh.
            if label not in known:
                skipped.append((cid, label))
                continue

            v = np.ascontiguousarray(f[f"meshes/{cid}/vertices_xyz"][:], dtype=np.float32)
            i = np.ascontiguousarray(f[f"meshes/{cid}/faces"][:], dtype=np.uint32).ravel()

            verts.append(v)
            idx.append(i)
            cells.append({
                "label": label, "cell_id": int(cid),
                "vertex_offset": v_off, "vertex_count": len(v),
                "index_offset": i_off, "index_count": len(i),
            })
            v_off += len(v)
            i_off += len(i)

    if skipped:
        print(f"  skipped {len(skipped)} mesh(es) with no matching cell: {skipped}")
    if not cells:
        raise SystemExit("no meshes matched the dataset's cell labels")

    V = np.concatenate(verts)
    I = np.concatenate(idx)

    (out / "meshes").mkdir(parents=True, exist_ok=True)
    with gzip.open(out / "meshes/cells.bin.gz", "wb") as fh:
        # Indices are LOCAL to each cell, so a slice builds a geometry directly.
        fh.write(struct.pack("<IIII", 1, len(cells), len(V), len(I)))
        fh.write(V.tobytes())
        fh.write(I.tobytes())

    (out / "meshes/index.json").write_text(json.dumps({
        "version": 1,
        "source": "annotated_domain_viz_cell_meshes",
        "vertex_key": "vertices_xyz",
        "units": "um",
        "cells": cells,
    }, indent=2))

    manifest_path = out / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["has_meshes"] = True
    manifest_path.write_text(json.dumps(manifest, indent=2))

    size = (out / "meshes/cells.bin.gz").stat().st_size
    lo, hi = V.min(axis=0), V.max(axis=0)
    print(f"  {len(cells)} cells  {len(V):,} vertices  {len(I) // 3:,} triangles"
          f"  {size / 1e3:.0f} KB")
    print(f"  extent {' × '.join(f'{a:.1f}..{b:.1f}' for a, b in zip(lo, hi))} µm")


if __name__ == "__main__":
    main()
