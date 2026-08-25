#!/usr/bin/env python3
"""Inventory the benchmark drop folder (bench/datasets/) → bench/datasets.json.

Drop real-world datasets into bench/datasets/<format>/ and run this. It records
what the benchmark needs to interpret a result: file size, the shape that
actually drives memory (cells × genes, or molecules), and any structure known
to blow the peak up (a full-width matrix parked in obsm, raw/layers, dense X).

Format detection reuses the pipeline's own `detect_input_format`, so a dataset
that lands in the wrong folder is reported rather than silently benchmarked as
something else.

    python3 scripts/bench/discover.py            # inventory + table
    python3 scripts/bench/discover.py --json     # just the JSON path

Sizes are read from metadata only — no dataset is loaded, so this stays fast on
a 50 GB drop folder. Row counts for very large CSVs are sampled (flagged
"estimated": true).
"""
import argparse
import contextlib
import gzip
import importlib.util
import io
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DROP = ROOT / "bench" / "datasets"
DEFAULT_OUT = ROOT / "bench" / "datasets.json"

# Folders users drop into. "single-molecule" holds .parquet / .csv files.
FORMAT_DIRS = ["h5ad", "xenium", "merscope", "single-molecule"]
# Count rows exactly below this; sample-estimate above it.
EXACT_COUNT_LIMIT = 200 * 1024**2
GB = 1024**3
# How deep to look for the real dataset root inside a dropped folder.
MAX_NEST = 3

# The only files each folder format's processor actually opens. A Xenium export
# is mostly images/zarr/boundaries the pipeline never reads, so "relevant" size
# — not folder size — is what gets uploaded and what the benchmark should quote.
RELEVANT = {
    "xenium": ["cells.csv", "cells.csv.gz", "cell_feature_matrix",
               "cell_feature_matrix.tar.gz", "cell_feature_matrix.tar", "cell_feature_matrix.zip",
               "features.tsv", "features.tsv.gz",
               "matrix.mtx", "matrix.mtx.gz", "barcodes.tsv", "barcodes.tsv.gz",
               "analysis", "analysis.tar.gz", "analysis.zip"],
    "merscope": ["cell_metadata.csv", "cell_metadata.csv.gz",
                 "cell_by_gene.csv", "cell_by_gene.csv.gz",
                 "cell_categories.csv", "cell_categories.csv.gz"],
}
# A folder is a dataset root when it holds one of these.
MARKERS = {
    "xenium": ["cells.csv", "cells.csv.gz"],
    "merscope": ["cell_metadata.csv", "cell_metadata.csv.gz"],
}


def _load_pipeline():
    """`detect_input_format` from the processor (worker/entrypoint.py needs boto3)."""
    spec = importlib.util.spec_from_file_location(
        "psd", ROOT / "scripts" / "process_spatial_data.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def dir_bytes(p: Path) -> int:
    if p.is_file():
        return p.stat().st_size
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())


def _open_text(p: Path):
    if p.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def count_rows(p: Path):
    """(rows_excluding_header, estimated?) — exact for small files, sampled for big ones."""
    size = p.stat().st_size
    if size == 0:
        return 0, False
    if size <= EXACT_COUNT_LIMIT or p.suffix == ".gz":
        with _open_text(p) as f:
            n = sum(1 for _ in f)
        return max(0, n - 1), False
    with open(p, "rb") as f:
        head = f.read(1 * 1024**2)
    lines = head.split(b"\n")[1:-1]  # drop header + partial last line
    if not lines:
        return None, True
    avg = sum(len(x) + 1 for x in lines) / len(lines)
    return int(size / avg) - 1, True


def header_cols(p: Path):
    with _open_text(p) as f:
        line = f.readline()
    return [c.strip().strip('"') for c in line.rstrip("\n").split(",")] if line else []


def find_one(base: Path, names):
    for n in names:
        c = base / n
        if c.exists():
            return c
    return None


# --- per-format probes --------------------------------------------------------

def probe_h5ad(path: Path) -> dict:
    import h5py

    info, notes = {}, []
    with h5py.File(path, "r") as h:
        x = h.get("X")
        cells = genes = None
        if isinstance(x, h5py.Dataset):
            cells, genes = int(x.shape[0]), int(x.shape[1])
            info["xEncoding"] = f"dense {x.dtype}"
        elif x is not None:  # sparse group
            shape = x.attrs.get("shape")
            if shape is not None:
                cells, genes = int(shape[0]), int(shape[1])
            enc = x.attrs.get("encoding-type", b"sparse")
            enc = enc.decode() if isinstance(enc, bytes) else str(enc)
            nnz = int(x["data"].shape[0]) if "data" in x else None
            info["xEncoding"] = enc
            info["nnz"] = nnz
            if nnz and cells and genes:
                info["density"] = round(nnz / (cells * genes), 4)

        obs = h.get("obs")
        if obs is not None:
            cols = [k for k in obs.keys() if not k.startswith("_")]
            info["obsColumns"] = len(cols)

        obsm = {}
        if "obsm" in h:
            for k, v in h["obsm"].items():
                if isinstance(v, h5py.Dataset) and v.ndim == 2:
                    obsm[k] = [int(v.shape[0]), int(v.shape[1])]
                    # A full-width matrix in obsm is read like X would be.
                    if genes and v.shape[1] >= genes:
                        notes.append(
                            f"obsm['{k}'] is full-width ({v.shape[0]:,}×{v.shape[1]:,}"
                            f" ≈ {v.shape[0] * v.shape[1] * 4 / GB:.1f} GB dense f32)"
                        )
                elif hasattr(v, "attrs") and "column-order" in getattr(v, "attrs", {}):
                    ncol = len(v.attrs["column-order"])
                    obsm[k] = ["dataframe", ncol]
                    if ncol > 3:
                        notes.append(f"obsm['{k}'] is a {ncol}-column dataframe")
        info["obsm"] = obsm

        for grp in ("raw", "layers"):
            if grp in h:
                keys = list(h[grp].keys())
                if keys:
                    info[grp] = keys
                    notes.append(f"has {grp} ({', '.join(keys[:3])}) — extra copies of X")

    if cells and genes:
        notes.insert(0, f"dense f32 matrix ≈ {cells * genes * 4 / GB:.1f} GB")
    return {"cells": cells, "genes": genes, "molecules": None}, info, notes


def probe_single_molecule(path: Path) -> dict:
    info, notes = {}, []
    molecules = None
    if path.suffix.lower() == ".parquet":
        import pyarrow.parquet as pq

        md = pq.read_metadata(str(path))  # metadata only: free
        molecules = int(md.num_rows)
        cols = list(pq.read_schema(str(path)).names)
        info["rowGroups"] = md.num_row_groups
    else:
        molecules, est = count_rows(path)
        cols = header_cols(path)
        if est:
            info["rowsEstimated"] = True
            notes.append("molecule count is a sample-based estimate")
    info["columns"] = cols
    # Mirrors detect_molecule_type() in worker/entrypoint.py (not importable
    # here: that module needs boto3).
    if "global_x" in cols:
        info["preset"] = "merscope"
        info["cellIdCol"] = "cell_id" if "cell_id" in cols else None
    elif "x_location" in cols:
        info["preset"] = "xenium"
        info["cellIdCol"] = None
    else:
        info["preset"] = None
        notes.append("no 'global_x' (MERSCOPE) or 'x_location' (Xenium) column — the worker would reject this")
    if info.get("cellIdCol"):
        notes.append(f"has cell-id column '{info['cellIdCol']}' → assigned/unassigned split")
    return {"cells": None, "genes": None, "molecules": molecules}, info, notes


def probe_xenium(path: Path) -> dict:
    info, notes = {}, []
    cells_file = find_one(path, ["cells.csv", "cells.csv.gz"])
    cells = None
    if cells_file:
        cells, est = count_rows(cells_file)
        if est:
            info["rowsEstimated"] = True
        info["obsColumns"] = max(0, len(header_cols(cells_file)) - 3)
    else:
        notes.append("no cells.csv(.gz) — the processor would reject this folder")
    feat = find_one(path, ["cell_feature_matrix/features.tsv", "cell_feature_matrix/features.tsv.gz",
                           "features.tsv", "features.tsv.gz"])
    genes = None
    if feat:
        with _open_text(feat) as f:
            lines = [ln for ln in (l.strip() for l in f) if ln]
        # A real Xenium features.tsv is headerless (id \t name \t type); some
        # exports carry a one-column header row — don't count it as a gene.
        if lines and len(lines[0].split("\t")) == 1 and lines[0].lower() in (
            "gene", "genes", "feature", "features", "name", "id"
        ):
            lines = lines[1:]
        genes = len(lines)
    if genes is None:
        # No loose features.tsv (10x ships it inside cell_feature_matrix.tar.gz)
        # — the .h5 carries the matrix shape in its metadata, read for free.
        h5 = path / "cell_feature_matrix.h5"
        if h5.exists():
            try:
                import h5py
                with h5py.File(h5, "r") as f:
                    shp = f["matrix"]["shape"][()]
                    genes = int(shp[0])
                    info["h5Cells"] = int(shp[1])
            except Exception:  # noqa: BLE001
                pass
    # process_spatial_data.py reads matrix.mtx(.gz), or extracts
    # cell_feature_matrix.{zip,tar,tar.gz} to find one.
    mtx = find_one(path, ["cell_feature_matrix/matrix.mtx", "cell_feature_matrix/matrix.mtx.gz",
                          "matrix.mtx", "matrix.mtx.gz"])
    archive = find_one(path, ["cell_feature_matrix.zip", "cell_feature_matrix.tar",
                              "cell_feature_matrix.tar.gz"])
    info["expressionSource"] = ("matrix.mtx" if mtx else archive.name if archive else None)
    if not mtx and not archive:
        if (path / "cell_feature_matrix.h5").exists():
            notes.append("only cell_feature_matrix.h5 — the loader reads matrix.mtx or a "
                         "cell_feature_matrix .zip/.tar/.tar.gz, so expression would be empty")
        else:
            notes.append("no expression matrix (matrix.mtx or cell_feature_matrix archive) — "
                         "gene expression would be empty")
    if find_one(path, ["transcripts.parquet", "transcripts.csv", "transcripts.csv.gz"]):
        notes.append("contains a transcripts file → also usable as a single-molecule input")
    return {"cells": cells, "genes": genes, "molecules": None}, info, notes


def probe_merscope(path: Path) -> dict:
    info, notes = {}, []
    meta = find_one(path, ["cell_metadata.csv", "cell_metadata.csv.gz"])
    cells = None
    est = False
    if meta:
        cells, est = count_rows(meta)
        if est:
            info["rowsEstimated"] = True
        info["obsColumns"] = max(0, len(header_cols(meta)) - 3)
    else:
        notes.append("no cell_metadata.csv at any depth — raw/unsegmented data, "
                     "not a MERSCOPE analysis output the processor can read")
    cbg = find_one(path, ["cell_by_gene.csv", "cell_by_gene.csv.gz"])
    genes = None
    if cbg:
        genes = max(0, len(header_cols(cbg)) - 1)
        # The loader requires both tables to describe the same cells in the same
        # order and aborts otherwise, so a mismatch here means the folder is not
        # a usable dataset — cheaper to learn now than after an upload + a run.
        cbg_rows, cbg_est = count_rows(cbg)
        info["cellByGeneRows"] = cbg_rows
        if cells and cbg_rows:
            ratio = cbg_rows / cells
            if ratio < 0.95 or ratio > 1.05:
                notes.append(
                    f"ROW MISMATCH: cell_metadata.csv ~{cells:,} rows vs "
                    f"cell_by_gene.csv ~{cbg_rows:,}"
                    + (" (both sampled)" if (est or cbg_est) else "")
                    + " — the processor rejects this folder"
                )
    else:
        notes.append("no cell_by_gene.csv — gene expression would be empty")
    if find_one(path, ["detected_transcripts.csv", "detected_transcripts.csv.gz"]):
        notes.append("contains detected_transcripts → also usable as a single-molecule input")
    return {"cells": cells, "genes": genes, "molecules": None}, info, notes


PROBES = {
    "h5ad": probe_h5ad,
    "xenium": probe_xenium,
    "merscope": probe_merscope,
    "single-molecule": probe_single_molecule,
}


def slug(name: str) -> str:
    keep = [c if (c.isalnum() or c in "-_") else "-" for c in name]
    return "".join(keep).strip("-")[:60]


def dataset_slug(path: Path, fmt: str) -> str:
    """Stable id fragment. Single-molecule keeps its extension — the same export
    can be present as both .parquet and .csv — and keeps it through truncation."""
    if fmt == "single-molecule" and path.is_file():
        return f"{slug(path.stem)[:50]}-{path.suffix.lstrip('.').lower()}"
    return slug(path.stem if path.is_file() else path.name)


def relevant_bytes(path: Path, fmt: str):
    """Bytes of the files this format's processor actually reads (None if n/a)."""
    names = RELEVANT.get(fmt)
    if not names or not path.is_dir():
        return None
    total = 0
    for n in names:
        p = path / n
        if p.exists():
            total += dir_bytes(p)
    return total


def resolve_root(path: Path, fmt: str) -> Path:
    """The folder that actually holds the dataset.

    Real exports are often wrapped in a download folder (logs + a nested region
    directory), so a dropped folder may be one or two levels above the files.
    """
    markers = MARKERS.get(fmt)
    if not markers or not path.is_dir():
        return path
    if any((path / m).exists() for m in markers):
        return path
    best = None
    stack = [(path, 0)]
    while stack:
        d, depth = stack.pop()
        if depth >= MAX_NEST:
            continue
        try:
            children = sorted(c for c in d.iterdir() if c.is_dir() and not c.name.startswith("."))
        except OSError:
            continue
        for c in children:
            if any((c / m).exists() for m in markers):
                # Prefer the largest match: a multi-region export has several.
                size = sum((c / m).stat().st_size for m in markers if (c / m).exists())
                if best is None or size > best[1]:
                    best = (c, size)
            stack.append((c, depth + 1))
    return best[0] if best else path


def entries_in(dirpath: Path):
    """One dataset per top-level entry: a file (h5ad/parquet/csv) or a folder."""
    if not dirpath.is_dir():
        return []
    out = []
    for p in sorted(dirpath.iterdir()):
        # "._x" are macOS AppleDouble sidecars, not datasets.
        if p.name.startswith(".") or p.name.startswith("._"):
            continue
        if p.name.lower() in ("readme.md", "readme"):
            continue
        out.append(p)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=DEFAULT_DROP,
                    help=f"databank root (default: {DEFAULT_DROP.relative_to(ROOT)})")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help="where to write the inventory JSON")
    ap.add_argument("--json", action="store_true", help="print only the output path")
    args = ap.parse_args()

    drop, out_path = args.root.resolve(), args.out
    if not drop.is_dir():
        print(f"No such databank: {drop}", file=sys.stderr)
        return 2

    psd = _load_pipeline()
    datasets, problems = [], []

    for fmt in FORMAT_DIRS:
        for path in entries_in(drop / fmt):
            rel = path.relative_to(drop).as_posix()

            # A single-molecule dataset is one file; folders there are the
            # source exports the flat symlinks point into (see link-transcripts.py).
            if fmt == "single-molecule" and path.is_dir():
                continue

            # Real exports are often wrapped in a download folder.
            root_path = resolve_root(path, fmt)
            size = dir_bytes(root_path)
            rec = {
                "id": f"{fmt}__{dataset_slug(path, fmt)}",
                "declaredFormat": fmt,
                "path": rel,
                "bytes": size,
                "sizeGB": round(size / GB, 3),
            }
            if root_path != path:
                rec["datasetRoot"] = root_path.relative_to(drop).as_posix()

            # Only some of a folder export is ever read/uploaded.
            rb = relevant_bytes(root_path, fmt)
            if rb is not None:
                rec["relevantBytes"] = rb
                rec["relevantGB"] = round(rb / GB, 3)

            # Cross-check the folder against the pipeline's own detector.
            if fmt == "single-molecule":
                detected = "single-molecule" if path.suffix.lower() in (".parquet", ".csv") else None
                if detected is None:
                    problems.append(f"{rel}: expected a .parquet or .csv single-molecule file")
            else:
                try:
                    # detect_input_format() narrates to stdout; keep our table clean.
                    with contextlib.redirect_stdout(io.StringIO()):
                        detected = psd.detect_input_format(root_path)
                except Exception as e:  # noqa: BLE001
                    detected = None
                    problems.append(f"{rel}: format not detected ({e})")
            rec["detectedFormat"] = detected
            if detected and detected != fmt:
                problems.append(f"{rel}: sits in {fmt}/ but detects as '{detected}'")

            try:
                shape, info, notes = PROBES[fmt](root_path)
            except Exception as e:  # noqa: BLE001
                shape, info, notes = ({"cells": None, "genes": None, "molecules": None},
                                      {}, [f"probe failed: {e}"])
                problems.append(f"{rel}: probe failed ({e})")
            rec["shape"] = shape
            rec["detail"] = info
            rec["notes"] = notes
            if rb is not None and size and rb < size * 0.5:
                rec["notes"].append(
                    f"only {rb / GB:.2f} GB of {size / GB:.1f} GB is read by the pipeline "
                    f"(images/zarr/boundaries are skipped)"
                )
            datasets.append(rec)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps({
        "generatedAt": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "root": str(drop),
        "datasets": datasets,
    }, indent=1) + "\n")

    if args.json:
        print(out_path)
        return 0

    if not datasets:
        print(f"No datasets found. Drop files into {drop}/<format>/ "
              f"({', '.join(FORMAT_DIRS)}) and re-run — see bench/README.md.")
        return 0

    def num(v):
        return format(v, ",") if v else "-"

    print("{:<46} {:>8} {:>8}  {:>11} {:>7} {:>13}".format(
        "id", "GB", "used", "cells", "genes", "molecules"))
    print("-" * 100)
    total = 0
    for d in datasets:
        s = d["shape"]
        total += d["bytes"]
        used = "{:.2f}".format(d["relevantGB"]) if "relevantGB" in d else "-"
        print("{:<46} {:>8.2f} {:>8}  {:>11} {:>7} {:>13}".format(
            d["id"], d["sizeGB"], used, num(s["cells"]), num(s["genes"]), num(s["molecules"])))
        for n in d["notes"]:
            print("{:<46}   · {}".format("", n))
    print("-" * 100)
    print(f"{len(datasets)} dataset(s), {total / GB:.2f} GB on disk → {out_path}")
    if problems:
        print("\nProblems:")
        for p in problems:
            print(f"  ! {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
