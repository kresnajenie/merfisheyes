#!/usr/bin/env python3
"""Build chunked labelled-molecule datasets for the analysed spiralia embryos.

    python build_all.py                  # every analysed embryo
    python build_all.py MER1_E00 ...     # just these

Per embryo: dill object -> parquet -> chunked folder. The parquet is an
intermediate and is removed once the chunked folder exists.

Each result is checked against `embryo_count_noblank` in the embryo metadata
CSV, which is an independently produced molecule count — a mismatch means the
filtering diverged from the source analysis and the dataset should not ship.
"""

import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path("/home/data/yiqun-spiralia/Sep2026")
OUT = ROOT / "merfisheyes_export"
REPO = Path("/home/kjenie/merfisheyes")
PY_DILL = "/home/kjenie/miniconda3/envs/napari/bin/python"  # has dill + pyarrow
PY = "/home/kjenie/miniconda3/bin/python"

META = pd.read_csv(ROOT / "embryo_metadata_deDup_May12_2026.csv").set_index("embryo")


def build(emb: str) -> dict:
    pq = OUT / f"{emb}_molecules.parquet"
    out = OUT / f"{emb}_lm"
    log = OUT / f"{emb}.log"
    stage = str(META.loc[emb, "stage"]) if emb in META.index else None
    t0 = time.perf_counter()

    with open(log, "w") as lf:
        r = subprocess.run(
            [PY_DILL, str(REPO / "scripts/spiralia/export_parquet.py"), emb, str(OUT)],
            stdout=lf, stderr=subprocess.STDOUT,
        )
        if r.returncode != 0:
            return {"embryo": emb, "status": "export failed"}

        cmd = [
            PY, str(REPO / "scripts/process_labelled_molecules.py"), str(pq), str(out),
            "--category", "gene=feature_name",
            "--category", "domain=rna_domain_anno",
            "--category", "cell=cell_id", "--label", "cell=cell_name",
            "--drop-unassigned", "cell_id",
            "--drop-control-genes", "feature_name",
            "--name", emb,
        ]
        if stage:
            cmd += ["--stage", stage]
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            return {"embryo": emb, "status": "ingest failed"}

    manifest = json.loads((out / "manifest.json").read_text())
    obs = json.loads((out / "obs/metadata.json").read_text())
    n = manifest["statistics"]["total_cells"]
    expected = META.loc[emb, "embryo_count_noblank"] if emb in META.index else None
    matches = expected is not None and int(expected) == n

    pq.unlink(missing_ok=True)

    return {
        "embryo": emb, "status": "ok", "stage": stage, "molecules": n,
        "expected": int(expected) if expected is not None else None,
        "count_ok": matches,
        "genes": obs.get("gene", {}).get("unique_values"),
        "domains": obs.get("domain", {}).get("unique_values"),
        "cells": obs.get("cell", {}).get("unique_values"),
        "mb": round(sum(f.stat().st_size for f in out.rglob("*") if f.is_file()) / 1e6, 1),
        "secs": round(time.perf_counter() - t0),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    embryos = sys.argv[1:] or sorted(p.name for p in (ROOT / "Final_analyzed_objects").iterdir() if p.is_dir())
    rows = []

    for i, emb in enumerate(embryos, 1):
        r = build(emb)
        rows.append(r)
        flag = "" if r.get("count_ok", True) else "  <-- COUNT MISMATCH"
        print(
            f"[{i}/{len(embryos)}] {emb:<15} {r['status']:<14} "
            f"{str(r.get('stage','')):<12} {r.get('molecules','') or '':>10} "
            f"g{r.get('genes','-')} d{r.get('domains','-')} c{r.get('cells','-')} "
            f"{r.get('mb','')}MB {r.get('secs','')}s{flag}",
            flush=True,
        )

    pd.DataFrame(rows).to_csv(OUT / "build_report.csv", index=False)
    ok = [r for r in rows if r["status"] == "ok"]
    bad = [r for r in ok if not r["count_ok"]]

    print(f"\nbuilt {len(ok)}/{len(rows)}; count mismatches: {len(bad)}")
    for r in rows:
        if r["status"] != "ok":
            print(f"  FAILED {r['embryo']}: {r['status']}")
    for r in bad:
        print(f"  MISMATCH {r['embryo']}: got {r['molecules']:,} expected {r['expected']:,}")


if __name__ == "__main__":
    main()
