#!/usr/bin/env python3
"""Run scripts/prepare.py on each tiny source and keep what it writes.

  python scripts/testdata/make-tiny-parity-stores.py

tests/unit/preview-store-parity.test.ts opens each source the way the browser
preview does and the store prepare.py made from it the way the SpatialData
viewer does, and compares what a person would see (docs/plans/PREVIEW-SAVE.md,
D5). The unit tests run without Python, so the stores are committed here;
tests/python/test_parity_stores.py fails when they no longer match a fresh
prepare.py run, and this script is how to regenerate them.

Outputs tests/data/tiny/parity/<case>/: the store (<case>.zarr) and
mapping.json. stats.json is left out: it records the machine's package
versions and the git commit, which change on every run.

Run with the packages of scripts/requirements.txt.
"""

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import prepare  # noqa: E402

TINY = ROOT / "tests" / "data" / "tiny"
OUT = TINY / "parity"

# case -> source, as a person would drop it on the preview
CASES = {
    "tiny_h5ad": TINY / "h5ad" / "tiny.h5ad",
    "nullable": TINY / "h5ad" / "nullable.h5ad",
    "named_index": TINY / "h5ad" / "named_index.h5ad",
    "missing_xy": TINY / "h5ad" / "missing_xy.h5ad",
    "merscope": TINY / "merscope" / "tiny",
}


def build(case: str, source: Path, out: Path) -> None:
    shutil.rmtree(out, ignore_errors=True)
    code = prepare.main([str(source), str(out), "--name", case])
    if code != 0:
        raise SystemExit(f"{case}: prepare.py exited {code}")
    (out / "stats.json").unlink(missing_ok=True)


if __name__ == "__main__":
    for case, source in CASES.items():
        build(case, source, OUT / case)
        print(f"wrote {OUT / case}")
