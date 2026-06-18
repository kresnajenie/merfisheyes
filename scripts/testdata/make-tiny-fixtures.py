#!/usr/bin/env python3
"""Build the committed *tiny* test fixtures for every supported format.

These are intentionally small (default 300 cells x 40 genes) so they can be
committed to the repo and used in the fast `@quick` CI suite. Run once after
checkout if the fixtures are missing or you change the tiny spec; the output is
deterministic for a given seed, so re-running yields identical data.

  python scripts/testdata/make-tiny-fixtures.py

Outputs into tests/data/tiny/:
  h5ad/tiny.h5ad
  xenium/tiny/...
  merscope/tiny/...
  single-molecule/tiny.parquet
  chunked/tiny/...            (produced from the h5ad via process_spatial_data.py)
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TINY = ROOT / "tests" / "data" / "tiny"
GEN = ROOT / "scripts" / "testdata" / "generate.py"
PROCESS = ROOT / "scripts" / "process_spatial_data.py"

CELLS = 300
GENES = 40
MOLECULES = 5000
SEED = 42


def run(cmd):
    print("  $", " ".join(str(c) for c in cmd))
    subprocess.run([sys.executable, *map(str, cmd)], check=True)


def main():
    print(f"Generating tiny fixtures ({CELLS} cells x {GENES} genes) into {TINY}")

    run([GEN, "--format", "h5ad", "--cells", CELLS, "--genes", GENES,
         "--seed", SEED, "--out", TINY / "h5ad" / "tiny.h5ad"])
    run([GEN, "--format", "xenium", "--cells", CELLS, "--genes", GENES,
         "--seed", SEED, "--out", TINY / "xenium" / "tiny"])
    run([GEN, "--format", "merscope", "--cells", CELLS, "--genes", GENES,
         "--seed", SEED, "--out", TINY / "merscope" / "tiny"])
    # 3D: the Xenium single-molecule column mapping requires a z_location column.
    run([GEN, "--format", "single-molecule", "--genes", GENES,
         "--molecules", MOLECULES, "--seed", SEED, "--dims", "3",
         "--out", TINY / "single-molecule" / "tiny.parquet"])

    # Chunked folder: reuse the production preprocessing pipeline on the h5ad.
    chunked_out = TINY / "chunked" / "tiny"
    if PROCESS.exists():
        try:
            run([PROCESS, TINY / "h5ad" / "tiny.h5ad", chunked_out, "--format", "h5ad"])
        except subprocess.CalledProcessError as e:
            print(f"  WARN: chunked generation failed ({e}); skipping chunked fixture.")
    else:
        print("  WARN: scripts/process_spatial_data.py not found; skipping chunked fixture.")

    print("Done. Remember to commit tests/data/tiny/ and tests/data/manifest.json.")


if __name__ == "__main__":
    main()
