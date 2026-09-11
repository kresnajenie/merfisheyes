#!/usr/bin/env python3
"""Export per-cell nuclei surfaces alongside a labelled-molecule dataset.

    python export_nuclei.py MER6-2_E3_1 /path/to/MER6-2_E3_1_lm

Writes into the dataset folder:

    meshes/nuclei_index.json   per-nucleus label + vertex/index slices
    meshes/nuclei.bin.gz       Float32 vertices, then Uint32 indices

and sets `has_nuclei: true` in the manifest. Same binary layout as
export_meshes.py writes for cells, so the viewer reads both the same way.

Source is `Segmentation/Bogdan/Old/{embryo}_segm.npz`, the volume the analysis
notebook loads under "# Load corresponding nuclei masks". Despite the napari
layer being named "nuclei" and its labels being *cell* ids, the geometry really
is nuclear: it fills 0.12% of the volume where the DS4 cell masks fill ~30%,
a nucleus:cell ratio of about 0.4%, which is right for blastomeres this yolky.

Frame, verified rather than assumed:

    segm is 2,4,4 downsampled in (z,x,y) from Xh, and Xh is pixels at 0.5 µm
    in z and 0.177 µm in xy, so one segm voxel is 1.0 µm in z and 0.708 µm in
    x and y. On MER6-2_E3_1 that reproduces the object's own `cms_nuc` centres
    for all 26 nuclei to within one voxel, and every nucleus lands inside its
    own cell's molecule cloud.

Decimated to step 2 by default: full resolution is 1.76 MB gzipped against
294 KB for the same embryo's cell meshes, and stepping to 2 costs under 0.5 µm
of extent on features tens of µm across while cutting it to 410 KB.

Labels must match the dataset's `cell` obs dictionary exactly, including the
"name (id)" suffix that process_labelled_molecules.py applies when two cells
share a name. That is checked, not assumed.
"""

import argparse
import gzip
import json
import struct
from pathlib import Path

import dill
import numpy as np
from skimage import measure

ROOT = Path("/home/data/yiqun-spiralia/Sep2026")
SEGM_DIR = ROOT / "Segmentation/Bogdan/Old"

# One segm voxel, in µm. The 2/4 are the z/xy downsample of segm relative to
# Xh; the 0.5/0.177 are Xh's own pixel size, as used by export_parquet.py.
UM_Z = 2 * 0.5
UM_XY = 4 * 0.177


def dataset_cell_labels(out: Path) -> list[str]:
    """The `cell` obs dictionary — the labels the viewer filters on."""
    raw = gzip.open(out / "obs/cell.bin.gz", "rb").read()
    _, _, _, _, _, dict_len = struct.unpack_from("<IIBBHI", raw, 0)

    return json.loads(raw[16 : 16 + dict_len].decode())


def display_labels(embryo: str, present: set[int]) -> dict[int, str]:
    """cell_id -> display label, reproducing the ingest's disambiguation.

    process_labelled_molecules.disambiguate() keys on the cell id and suffixes
    only labels that collide, so two polar bodies become "pb (1)" and "pb (3)"
    while every unique name is left alone.
    """
    obj = ROOT / f"Final_analyzed_objects/{embryo}/{embryo}_RNA_domain_annotated"
    names = dill.load(open(obj, "rb"))["cell_names"]
    plain = {
        int(i): str(v).split("::")[1]
        for i, v in names.items()
        if int(i) in present
    }
    counts: dict[str, int] = {}

    for v in plain.values():
        counts[v] = counts.get(v, 0) + 1

    return {
        i: (v if counts[v] == 1 else f"{v} ({i})") for i, v in plain.items()
    }


def surface(mask: np.ndarray, step: int):
    """Marching cubes on one nucleus, in µm, as (x, y, z) vertices.

    The volume is padded first: a nucleus touching the array edge would
    otherwise come out as an open shell with a hole where it was clipped.
    """
    padded = np.pad(mask, 1)
    verts, faces, _, _ = measure.marching_cubes(
        padded.astype(np.float32),
        level=0.5,
        spacing=(UM_Z, UM_XY, UM_XY),
        step_size=step,
    )
    # marching_cubes returns vertices in array order, which here is (z, x, y).
    return verts[:, [1, 2, 0]].astype(np.float32), faces.astype(np.uint32)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("embryo")
    ap.add_argument("dataset_dir")
    ap.add_argument(
        "--step",
        type=int,
        default=2,
        help="marching-cubes step size; >1 decimates (default 2)",
    )
    args = ap.parse_args()

    out = Path(args.dataset_dir)
    segm = np.load(SEGM_DIR / f"{args.embryo}_segm.npz")["segm"]

    ids = [int(v) for v in np.unique(segm) if v != 0]
    labels = display_labels(args.embryo, set(ids))
    known = set(dataset_cell_labels(out))

    verts: list[np.ndarray] = []
    idx: list[np.ndarray] = []
    nuclei = []
    v_off = i_off = 0
    skipped = []

    for cid in ids:
        label = labels.get(cid)

        # A surface we cannot tie to a selectable cell is worse than none: it
        # would float in the scene with no way to hide it.
        if label not in known:
            skipped.append((cid, label))
            continue

        mask = segm == cid
        # Crop to the nucleus before marching: the full volume is 194M voxels
        # and all but a rounding error of it is background for any one cell.
        z, x, y = mask.nonzero()
        sub = mask[z.min() : z.max() + 1, x.min() : x.max() + 1, y.min() : y.max() + 1]
        v, f = surface(sub, args.step)

        # Undo the crop, in µm.
        v += np.array(
            [(x.min() - 1) * UM_XY, (y.min() - 1) * UM_XY, (z.min() - 1) * UM_Z],
            dtype=np.float32,
        )

        f = f.ravel()
        verts.append(v)
        idx.append(f)
        nuclei.append({
            "label": label, "cell_id": cid,
            "vertex_offset": v_off, "vertex_count": len(v),
            "index_offset": i_off, "index_count": len(f),
        })
        v_off += len(v)
        i_off += len(f)

    if skipped:
        print(f"  skipped {len(skipped)} nucleus/nuclei with no matching cell: {skipped}")
    if not nuclei:
        raise SystemExit("no nuclei matched the dataset's cell labels")

    V = np.concatenate(verts)
    I = np.concatenate(idx)

    (out / "meshes").mkdir(parents=True, exist_ok=True)
    with gzip.open(out / "meshes/nuclei.bin.gz", "wb") as fh:
        # Indices are LOCAL to each nucleus, so a slice builds a geometry.
        fh.write(struct.pack("<IIII", 1, len(nuclei), len(V), len(I)))
        fh.write(V.tobytes())
        fh.write(I.tobytes())

    (out / "meshes/nuclei_index.json").write_text(json.dumps({
        "version": 1,
        "source": "Segmentation/Bogdan/Old/*_segm.npz",
        "units": "um",
        "voxel_um": [UM_XY, UM_XY, UM_Z],
        "step_size": args.step,
        "cells": nuclei,
    }, indent=2))

    manifest_path = out / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["has_nuclei"] = True
    manifest_path.write_text(json.dumps(manifest, indent=2))

    size = (out / "meshes/nuclei.bin.gz").stat().st_size
    lo, hi = V.min(axis=0), V.max(axis=0)
    print(f"  {len(nuclei)} nuclei  {len(V):,} vertices  {len(I) // 3:,} triangles"
          f"  {size / 1e3:.0f} KB")
    print(f"  extent {' × '.join(f'{a:.1f}..{b:.1f}' for a, b in zip(lo, hi))} µm")


if __name__ == "__main__":
    main()
