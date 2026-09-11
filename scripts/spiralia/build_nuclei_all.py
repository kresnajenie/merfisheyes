#!/usr/bin/env python3
"""Export nuclei surfaces for every built spiralia dataset.

    python build_nuclei_all.py [--out DIR] [--step 2] [--only MER6-2_E3_1]

Reads build_report.csv for the list, runs export_nuclei.py against each built
dataset folder, and verifies the result before moving on. Unlike the cell
meshes — which exist as prebuilt files for only 31 of the 45 — every embryo has
a `_segm.npz`, so this covers all of them.

Verification per embryo, because a surface in the wrong frame renders
plausibly and is wrong:

  * every nucleus label is in the dataset's `cell` obs dictionary
  * the nuclei extent sits inside the molecules' own bounding box

Failures are recorded and the run continues; the summary lists them.
"""

import argparse
import csv
import gzip
import json
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = Path("/home/data/yiqun-spiralia/Sep2026/merfisheyes_export")


def cell_labels(ds: Path) -> set[str]:
    raw = gzip.open(ds / "obs/cell.bin.gz", "rb").read()
    _, _, _, _, _, dict_len = struct.unpack_from("<IIBBHI", raw, 0)

    return set(json.loads(raw[16 : 16 + dict_len].decode()))


def molecule_bounds(ds: Path) -> tuple[np.ndarray, np.ndarray]:
    """Molecule bounding box in µm, dequantised from the q16 coordinates."""
    raw = gzip.open(ds / "coords/spatial.q16.bin.gz", "rb").read()
    n, dims = struct.unpack_from("<II", raw, 0)
    off = 8
    mn = np.frombuffer(raw, np.float32, dims, off)
    off += 4 * dims
    sc = np.frombuffer(raw, np.float32, dims, off)
    off += 4 * dims
    codes = np.frombuffer(raw, np.uint16, n * dims, off).reshape(n, dims)

    return mn, mn + codes.max(axis=0).astype(np.float32) * sc


def verify(ds: Path) -> str | None:
    """None when the export looks right, else why it doesn't."""
    idx = json.loads((ds / "meshes/nuclei_index.json").read_text())
    known = cell_labels(ds)
    unknown = [c["label"] for c in idx["cells"] if c["label"] not in known]

    if unknown:
        return f"{len(unknown)} label(s) not in the cell column: {unknown[:3]}"

    raw = gzip.open(ds / "meshes/nuclei.bin.gz", "rb").read()
    _, _, n_verts, _ = struct.unpack_from("<IIII", raw, 0)
    V = np.frombuffer(raw, np.float32, n_verts * 3, 16).reshape(-1, 3)

    lo, hi = molecule_bounds(ds)
    # A tolerance, not equality: a nucleus can bulge a little past the outermost
    # molecule, and the q16 bounds are themselves quantised.
    pad = 0.05 * (hi - lo)

    if (V.min(axis=0) < lo - pad).any() or (V.max(axis=0) > hi + pad).any():
        return (
            f"nuclei extent {V.min(axis=0).round(1)}..{V.max(axis=0).round(1)} "
            f"outside molecules {lo.round(1)}..{hi.round(1)}"
        )

    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--step", type=int, default=2)
    ap.add_argument("--only", help="single embryo, for a quick check")
    args = ap.parse_args()

    rows = [
        r for r in csv.DictReader(open(HERE / "build_report.csv"))
        if r["status"] == "ok" and (not args.only or r["embryo"] == args.only)
    ]

    ok, failed = [], []

    for i, r in enumerate(rows, 1):
        emb = r["embryo"]
        ds = args.out / f"{emb}_lm"

        print(f"[{i}/{len(rows)}] {emb}", flush=True)
        if not ds.exists():
            failed.append((emb, "no built dataset folder"))
            print("  no built dataset folder")
            continue

        proc = subprocess.run(
            [sys.executable, str(HERE / "export_nuclei.py"), emb, str(ds),
             "--step", str(args.step)],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            failed.append((emb, proc.stderr.strip().splitlines()[-1:] or "failed"))
            print("  FAILED:", proc.stderr.strip().splitlines()[-1:])
            continue

        print(proc.stdout.rstrip())

        why = verify(ds)

        if why:
            failed.append((emb, why))
            print(f"  VERIFY FAILED: {why}")
        else:
            ok.append(emb)

    print(f"\n{len(ok)}/{len(rows)} exported and verified")
    for emb, why in failed:
        print(f"  FAILED {emb}: {why}")

    total = sum(
        (args.out / f"{e}_lm/meshes/nuclei.bin.gz").stat().st_size for e in ok
    )

    print(f"  {total / 1e6:.1f} MB of nuclei across {len(ok)} datasets")

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
