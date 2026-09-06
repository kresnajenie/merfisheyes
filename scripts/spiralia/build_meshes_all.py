#!/usr/bin/env python3
"""Export cell meshes for every embryo that has them, and verify each one.

    python build_meshes_all.py              # all embryos with meshes
    python build_meshes_all.py MER1_E00 ... # just these

Only 31 of the 45 analysed embryos have prebuilt meshes; the 14 without are the
MER2 set.

Every embryo is checked, not just the first: for each cell, the mesh centroid is
compared with the centroid of that cell's molecules, and the fraction of the
cell's molecules falling inside its own mesh bounding box is measured. That
check is what established `vertices_xyz` is the molecule frame and
`vertices_reoriented_xyz` is not — assuming the frames agree everywhere would be
exactly the mistake it caught.

Containment is the reliable signal, but it has to be weighted by molecule count.
Polar bodies are a few hundred molecules wrapped in a 15–30 vertex mesh, and
they score 30–70% on their own while every real cell scores 95–100%; gating on
the raw minimum flags healthy embryos. Centroid offset is likewise expected to
be large for yolk-heavy macromeres, whose molecules sit far from the surface
centroid, so it is reported but not gated.
"""

import gzip
import json
import struct
import subprocess
import sys
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("/home/data/yiqun-spiralia/Sep2026")
OUT = ROOT / "merfisheyes_export"
REPO = Path("/home/kjenie/merfisheyes")
MASK = ROOT / "Segmentation/Bogdan/annotated_domain_viz_cell_meshes"
PC = ROOT / "Segmentation/Bogdan/annotated_domain_viz_pointcloud_nonoverlap_cell_meshes"

# Share of an embryo's molecules that must fall inside their own cell's mesh.
# A frame mismatch drags this far below the threshold; a coarse polar-body mesh
# cannot, because it holds a negligible share of the molecules.
MIN_CONTAINMENT = 95.0
# Cells below this share of the embryo are too small for their containment to
# mean much on its own.
SMALL_CELL_SHARE = 0.005


def read_dataset(dir_: Path):
    """Molecule coordinates (µm) and per-molecule cell label index."""
    raw = gzip.open(dir_ / "coords/spatial.q16.bin.gz", "rb").read()
    n, dims = struct.unpack_from("<II", raw, 0)
    off = 8
    mn = np.frombuffer(raw, np.float32, dims, off)
    off += 4 * dims
    sc = np.frombuffer(raw, np.float32, dims, off)
    off += 4 * dims
    codes = np.frombuffer(raw, np.uint16, n * dims, off).reshape(n, dims)
    xyz = mn + codes.astype(np.float32) * sc

    b = gzip.open(dir_ / "obs/cell.bin.gz", "rb").read()
    _, ncell, _, width, _, dlen = struct.unpack_from("<IIBBHI", b, 0)
    dic = json.loads(b[16 : 16 + dlen].decode())
    o = 16 + dlen + ((-dlen) % 4)
    cc = np.frombuffer(b, {1: np.uint8, 2: np.uint16, 4: np.uint32}[width], ncell, o)

    return xyz, dic, cc


def verify(embryo: str, dir_: Path):
    xyz, dic, cc = read_dataset(dir_)

    with h5py.File(PC / f"{embryo}_cell_meshes.h5", "r") as f:
        ident = {c: str(f[f"meshes/{c}"].attrs["cell_identity"]) for c in f["meshes"]}
    counts = Counter(ident.values())
    disp = {c: (v if counts[v] == 1 else f"{v} ({c})") for c, v in ident.items()}

    offsets, unmatched = [], 0
    inside_total = counted_total = 0
    worst_label, worst_pct = None, 100.0
    total_molecules = len(cc)

    with h5py.File(MASK / f"{embryo}_cell_meshes.h5", "r") as f:
        for cid in f["meshes"]:
            label = disp.get(cid)

            if label not in dic:
                unmatched += 1
                continue
            P = xyz[cc == dic.index(label)]

            if len(P) == 0:
                continue
            v = f[f"meshes/{cid}/vertices_xyz"][:]
            lo, hi = v.min(0), v.max(0)
            inside = ((P >= lo - 1) & (P <= hi + 1)).all(1)

            offsets.append(float(np.linalg.norm(v.mean(0) - P.mean(0))))
            inside_total += int(inside.sum())
            counted_total += len(P)

            pct = float(inside.mean() * 100)

            # Only cells big enough for the number to mean something.
            if len(P) / total_molecules >= SMALL_CELL_SHARE and pct < worst_pct:
                worst_label, worst_pct = label, pct

    return {
        "cells": len(offsets),
        "unmatched": unmatched,
        "median_offset": float(np.median(offsets)) if offsets else float("nan"),
        "contained": 100.0 * inside_total / counted_total if counted_total else 0.0,
        "worst_cell": worst_label,
        "worst_pct": worst_pct,
    }


def main():
    embryos = sys.argv[1:] or sorted(
        p.name.replace("_cell_meshes.h5", "") for p in MASK.glob("*_cell_meshes.h5")
    )
    rows, problems = [], []

    print(f"{'embryo':<16}{'cells':>6}{'median Δ':>11}{'in mesh':>9}"
          f"{'worst real cell':>22}  status")
    for emb in embryos:
        dir_ = OUT / f"{emb}_lm"

        if not (dir_ / "manifest.json").exists():
            print(f"{emb:<16}{'':>6}{'':>11}{'':>10}  no dataset built — skipped")
            continue

        r = subprocess.run(
            [sys.executable, str(REPO / "scripts/spiralia/export_meshes.py"), emb, str(dir_)],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            print(f"{emb:<16}{'':>6}{'':>11}{'':>10}  EXPORT FAILED")
            problems.append((emb, r.stderr.strip().splitlines()[-1:] or ["?"]))
            continue

        v = verify(emb, dir_)
        bad = v["contained"] < MIN_CONTAINMENT
        status = "MISALIGNED" if bad else "ok"
        worst = (
            f"{v['worst_cell']} {v['worst_pct']:.0f}%" if v["worst_cell"] else "-"
        )

        print(f"{emb:<16}{v['cells']:>6}{v['median_offset']:>9.1f} µm"
              f"{v['contained']:>8.1f}%{worst:>22}  {status}")
        rows.append({"embryo": emb, **v})
        if bad:
            problems.append((emb, f"only {v['contained']:.1f}% of molecules in mesh"))

    unmatched = [(r["embryo"], r["unmatched"]) for r in rows if r["unmatched"]]

    print(f"\nexported {len(rows)} embryos")
    if unmatched:
        # Meshes for cells the ingest dropped (unnamed ids, unassigned). They
        # hold no molecules in the shipped dataset, so they are simply omitted.
        print("meshes with no matching cell (dropped at ingest, expected): "
              + ", ".join(f"{e}×{n}" for e, n in unmatched))
    if problems:
        print("needs a look:")
        for emb, why in problems:
            print(f"  {emb}: {why}")
    else:
        print("every embryo's meshes align with its molecules")


if __name__ == "__main__":
    main()
