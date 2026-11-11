#!/usr/bin/env python3
"""Rebuild a BIL sample's pipeline inputs from what the legacy pipeline put on S3.

  python reconstruct_inputs.py SAMPLE_DIR OUT [--molecules-csv] [--fit-sample N]

While the BIL cluster is down the raw slices cannot be read, but every sample's
legacy output is still in s3://merfisheyes-bil/bil-psc-data2/<sample>/:

  meyes_output/   the cells, as process_spatial_data.py wrote them
                  (manifest.json, coords/, expr/, obs/ — docs/PIPELINE-SPEC.md)
  sm_output/      the molecules, as process_single_molecule.py wrote them:
                  one folder per slice (sm_output/<slice>/) for a multi-slice
                  sample, or the folder itself for a single-slice one

Download those two folders (they are all that is publicly readable):

  aws s3 sync s3://merfisheyes-bil/bil-psc-data2/<sample>/ <sample>/ \\
      --exclude "*" --include "meyes_output/*" --include "sm_output/*"

and point this script at <sample>/. It writes, never touching SAMPLE_DIR:

  multi-slice (the cells have a `_sample_id` column):
    OUT/input/<slice>/cell_metadata.csv          EntityID, center_x/y, min/max, volume…
    OUT/input/<slice>/cell_by_gene.csv           cell, <gene>…   (integer counts)
    OUT/input/<slice>/detected_transcripts.parquet   gene, global_x, global_y[, global_z][, cell_id]
    → the shape combine_slices_v3.py reads:  tests/bil/run_dry_run.sh OUT/input …

  single slice (an h5ad or Xenium sample):
    OUT/input/cell_by_gene.h5ad                  X counts, obsm['spatial']
    OUT/input/segmented_spot_table.parquet       gene_names, x, y[, z][, cell_ids]
    → the shape prepare.py reads:  prepare.py OUT/input STORE_DIR

  both:
    OUT/labels/mapping_output.csv   the MapMyCells columns the old run added, keyed
                                    as the new pipeline keys cells, for --labels
    OUT/reconstruction.json         counts, each slice's recovered shift and how
                                    well it fits, and everything that was lost

What cannot come back, and what stands in for it
------------------------------------------------
* The cells of a multi-slice sample were grid-shifted into one frame by the
  combine, and the per-slice shift (sample_offsets.json) was never uploaded.
  The molecules, though, are still in each slice's own frame. So each slice's
  shift is **measured**: a coarse cross-correlation of cell density against
  assigned-molecule density, then refined until each assigned molecule's
  nearest cell is, on average, centred on it. The cells are moved back by it;
  the molecules are written exactly as stored. reconstruction.json records the
  shift and the fit (median distance from a molecule to its nearest cell).
* Cell ids: the slice's EntityIDs were stored as float32 and collapsed, so each
  slice's cells are numbered 1…n in their stored order.
* A molecule's cell: at most "assigned or not" was kept. Where the old run
  split them, `cell_id` is 1 for a molecule that was in a cell (which cell is
  unknown) and -1 for one that was not — all the new pipeline reads from it
  (meyes_common.unassigned_mask). The MERSCOPE runs never split them, so
  those tables have no cell column at all.
* Cells the old run dropped (its artifact mask) are gone. Masking the rebuilt
  sample again removes more; run the new pipeline with --mask-percentile 0.
* Control probes (Blank-*, NegControl*) were filtered out and do not return.
* Coordinates are the stored ones: raw microns rounded to 2 dp, float32.

Exit codes: 0 written; 2 SAMPLE_DIR is not a legacy sample; 3 written but
incomplete (files missing from the download); 4 written but a slice's shift
did not fit (see reconstruction.json before using it).
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import struct
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# The MapMyCells columns process_spatial_data.py merged into obs (--mmc-csv).
# They are a result, not an input, so they go to labels/mapping_output.csv.
LABEL_COLUMN = re.compile(
    r"^((class|subclass|supertype|cluster)_)?"
    r"(name|label|alias|bootstrapping_probability|correlation_coefficient)$"
)
LABEL_LEVELS = ("class", "subclass", "supertype", "cluster")
# What the combine added to each cell; it adds them again.
COMBINE_COLUMNS = {"_source_file", "_sample_id"}
# The cell id the combine read, collapsed by float32. Replaced by EntityID.
COLLAPSED_ID_COLUMNS = {"id", "EntityID", "cell_id"}
UNASSIGNED_SUFFIX = "_uuuuuuuuuu"

# A slice's shift is accepted when the median distance from an assigned
# molecule to its nearest cell is under this many microns (a cell's radius is
# ~5 µm; a wrong shift puts molecules tens of microns from any cell).
FIT_MEDIAN_UM = 6.0

T0 = time.perf_counter()


def log(msg: str) -> None:
    print(f"[{time.perf_counter() - T0:7.1f}s] {msg}", flush=True)


class NotLegacy(Exception):
    pass


# ── reading the legacy cell format ──────────────────────────────────────

def read_gz(path: Path) -> bytes:
    with gzip.open(path, "rb") as f:
        return f.read()


def read_coords(path: Path) -> np.ndarray:
    """coords/<name>.bin.gz: [u32 n][u32 dims][f32 n*dims]."""
    buf = read_gz(path)
    n, dims = struct.unpack_from("<II", buf, 0)
    return np.frombuffer(buf, np.float32, n * dims, 8).reshape(n, dims).copy()


def read_obs_column(obs_dir: Path, name: str, info: dict, n: int):
    """One obs column, binary v1 (lib/utils/obs-binary.ts) or a legacy JSON
    array. Returns (values, kind) or None when the file is missing."""
    if info.get("format") == "bin-v1":
        path = obs_dir / f"{name}.bin.gz"
        if not path.exists():
            return None
        buf = read_gz(path)
        version, count, kind, width, _, dict_bytes = struct.unpack_from("<IIBBHI", buf, 0)
        if version != 1 or count != n:
            raise ValueError(f"obs/{path.name}: version {version}, {count} cells (expected 1, {n})")
        offset = 16 + ((dict_bytes + 3) & ~3)
        if kind == 0:
            labels = np.array(json.loads(buf[16:16 + dict_bytes].decode("utf-8")), dtype=object)
            codes = np.frombuffer(buf, {1: np.uint8, 2: np.uint16, 4: np.uint32}[width], n, offset)
            return pd.Categorical.from_codes(codes.astype(np.int64), categories=pd.Index(labels).astype(str)), "categorical"
        return np.frombuffer(buf, np.float32, n, offset).copy(), "numerical"

    path = obs_dir / f"{name}.json.gz"
    if not path.exists():
        return None
    values = json.loads(read_gz(path).decode("utf-8"))
    if len(values) != n:
        raise ValueError(f"obs/{path.name}: {len(values)} values, expected {n}")
    if info.get("type") == "categorical":
        return pd.Categorical(["" if v is None else str(v) for v in values]), "categorical"
    return pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(np.float32), "numerical"


def read_expression(expr_dir: Path, n_cells: int, missing: list):
    """expr/index.json + chunk_*.bin.gz → (CSC float32 cells×genes, genes)."""
    from scipy import sparse

    index = json.loads((expr_dir / "index.json").read_text())
    genes = [None] * index["total_genes"]
    for g in index["genes"]:
        genes[g["chunk_id"] * index["chunk_size"] + g["position_in_chunk"]] = g["name"]

    rows: list = [np.empty(0, np.int32)] * len(genes)
    vals: list = [np.empty(0, np.float32)] * len(genes)
    for chunk in range(index["num_chunks"]):
        path = expr_dir / f"chunk_{chunk:05d}.bin.gz"
        if not path.exists():
            missing.append(str(path.relative_to(expr_dir.parent.parent)))
            continue
        buf = read_gz(path)
        _, num_genes, _, total = struct.unpack_from("<IIII", buf, 0)
        if total != n_cells:
            raise ValueError(f"{path.name}: written for {total} cells, the manifest has {n_cells}")
        for i in range(num_genes):
            gene_idx, offset, _, _, nnz, _ = struct.unpack_from("<IIIIII", buf, 16 + i * 24)
            rows[gene_idx] = np.frombuffer(buf, np.uint32, nnz, offset + 8).astype(np.int32)
            vals[gene_idx] = np.frombuffer(buf, np.float32, nnz, offset + 8 + 4 * nnz).copy()

    indptr = np.concatenate([[0], np.cumsum([r.size for r in rows])]).astype(np.int64)
    X = sparse.csc_matrix((np.concatenate(vals), np.concatenate(rows), indptr), shape=(n_cells, len(genes)))
    return X, genes


def read_cells(meyes: Path, missing: list) -> dict:
    manifest = json.loads((meyes / "manifest.json").read_text())
    if manifest.get("normalized", True):
        raise NotLegacy(f"{meyes}: normalized coordinates; only raw-micron outputs can be rebuilt")
    n = manifest["statistics"]["total_cells"]

    spatial = read_coords(meyes / "coords" / "spatial.bin.gz")
    if spatial.shape[0] != n:
        raise ValueError(f"coords/spatial.bin.gz has {spatial.shape[0]} cells, the manifest {n}")
    embeddings = {}
    for name in manifest["statistics"].get("available_embeddings", []):
        path = meyes / "coords" / f"{name}.bin.gz"
        if path.exists():
            embeddings[name] = read_coords(path)

    meta = json.loads((meyes / "obs" / "metadata.json").read_text())
    obs, kinds = {}, {}
    for name, info in meta.items():
        got = read_obs_column(meyes / "obs", name, info, n)
        if got is None:
            missing.append(f"meyes_output/obs/{name}")
            continue
        obs[name], kinds[name] = got
    obs = pd.DataFrame(obs)

    X, genes = read_expression(meyes / "expr", n, missing)
    integer = bool(np.all(X.data == np.round(X.data)))
    log(f"cells: {n:,} × {len(genes)} genes, {X.nnz:,} non-zeros "
        f"({'integer counts' if integer else 'NOT integer — not raw counts'}), "
        f"{spatial.shape[1]}-D, {len(obs.columns)} obs columns")
    return {"manifest": manifest, "spatial": spatial, "embeddings": embeddings,
            "obs": obs, "kinds": kinds, "X": X, "genes": genes, "integer": integer}


# ── reading the legacy molecule format ──────────────────────────────────

def sanitize_gene_name(gene: str) -> str:
    """The gene's file name — process_single_molecule.py's sanitize_gene_name."""
    name = re.sub(r"[^a-zA-Z0-9]", "_", gene)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")


def molecule_manifest(folder: Path) -> dict | None:
    for name in ("manifest.json.gz", "manifest.json"):
        path = folder / name
        if path.exists():
            raw = read_gz(path) if name.endswith(".gz") else path.read_bytes()
            return json.loads(raw.decode("utf-8"))
    return None


def iter_gene_files(folder: Path, manifest: dict, missing: list):
    """(gene, assigned xyz f32 (n,3), unassigned xyz f32 (m,3)) per gene.
    A file the manifest counts molecules for but the download lacks is
    recorded in `missing` and read as empty."""
    genes = manifest["genes"]["unique_gene_names"]
    counts = manifest["genes"].get("molecule_counts", {})
    names = [sanitize_gene_name(g) for g in genes]
    if len(set(names)) != len(names):
        raise ValueError(f"{folder}: two genes share a file name after sanitising")
    for gene, name in zip(genes, names):
        parts = []
        for suffix, key in (("", "assigned"), (UNASSIGNED_SUFFIX, "unassigned")):
            path = folder / "genes" / f"{name}{suffix}.bin.gz"
            count = counts.get(gene)
            expected = count.get(key) if isinstance(count, dict) else (count if not suffix else None)
            if path.exists():
                xyz = np.frombuffer(read_gz(path), np.float32).reshape(-1, 3)
                if expected is not None and xyz.shape[0] != expected:
                    raise ValueError(f"{path}: {xyz.shape[0]} molecules, the manifest says {expected}")
            else:
                # an unassigned file is only written when non-empty
                if (expected or 0) > 0 or (not suffix and expected is None):
                    missing.append(str(path))
                xyz = np.empty((0, 3), np.float32)
            parts.append(xyz)
        yield gene, parts[0], parts[1]


def write_molecules(folder: Path, out: Path, names: dict, fit_sample: int, missing: list) -> dict:
    """Stream one molecule folder into one parquet (or CSV) table, keeping a
    spread-out sample of the assigned molecules' x, y for fitting the shift."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    manifest = molecule_manifest(folder)
    dims = int(manifest["statistics"].get("spatial_dimensions", 3))
    total = int(manifest["statistics"]["total_molecules"])
    # Only a run given a cell column split assigned from unassigned (MERSCOPE's
    # default mapping had none). Without the split there is no cell column:
    # writing one would claim every molecule was in a cell.
    split = bool(manifest.get("has_unassigned"))
    counts = manifest["genes"].get("molecule_counts", {})
    total_assigned = sum(c.get("assigned", 0) for c in counts.values() if isinstance(c, dict)) or total
    keep = min(1.0, fit_sample / max(total_assigned, 1))
    rng = np.random.default_rng(0)
    coords = [names["x"], names["y"]] + ([names["z"]] if dims == 3 else [])
    cols = [names["gene"]] + coords + ([names["cell"]] if split else [])
    schema = pa.schema([(cols[0], pa.dictionary(pa.int32(), pa.string()))]
                       + [(c, pa.float32()) for c in coords] + ([(names["cell"], pa.int32())] if split else []))

    out.parent.mkdir(parents=True, exist_ok=True)
    csv = out.suffix == ".csv"
    writer = None if csv else pq.ParquetWriter(out, schema, compression="zstd")
    if csv:
        out.write_text(",".join(cols) + "\n")
    written = assigned = 0
    sample = []
    try:
        for gene, a, u in iter_gene_files(folder, manifest, missing):
            xyz = np.concatenate([a, u])
            if xyz.shape[0] == 0:
                continue
            cell = np.concatenate([np.ones(len(a), np.int32), np.full(len(u), -1, np.int32)])
            if len(a):
                sample.append(a[rng.random(len(a)) < keep, :2])
            if csv:
                frame = pd.DataFrame({cols[0]: gene, **{c: xyz[:, i] for i, c in enumerate(coords)},
                                      **({names["cell"]: cell} if split else {})})
                frame.to_csv(out, mode="a", header=False, index=False, float_format="%.2f")
            else:
                gene_col = pa.DictionaryArray.from_arrays(pa.array(np.zeros(len(xyz), np.int32)), pa.array([gene]))
                arrays = [gene_col] + [pa.array(xyz[:, i]) for i in range(dims)] + ([pa.array(cell)] if split else [])
                writer.write_table(pa.Table.from_arrays(arrays, schema=schema))
            written += len(xyz)
            assigned += len(a)
    finally:
        if writer is not None:
            writer.close()
    # What the gene files should hold: the per-gene counts. Older runs' total
    # also counted the control probes they filtered out, so it is only a fallback.
    per_gene = sum(c.get("assigned", 0) + c.get("unassigned", 0) for c in counts.values() if isinstance(c, dict))
    expected = per_gene or total
    if written != expected:
        log(f"  ⚠ {folder.name}: wrote {written:,} molecules, the manifest's gene counts say {expected:,}")
    elif per_gene and per_gene != total:
        log(f"  ({total - per_gene:,} control-probe molecules were filtered by the old run)")
    log(f"  molecules: {written:,} " + (f"({assigned:,} in a cell" if split else "(never split by cell")
        + f", {dims}-D) → {out.name}")
    fit = np.concatenate(sample) if sample else np.empty((0, 2), np.float32)
    return {"written": written, "manifest_total": total, "split_by_cell": split,
            "assigned": assigned if split else None,
            "dims": dims, "fit": fit.astype(np.float64)}


# ── measuring a slice's shift ──────────────────────────────────────────

def coarse_shift(cells: np.ndarray, molecules: np.ndarray, bin_um: float) -> np.ndarray:
    """d such that cells ≈ molecules + d, from the peak of the cross-correlation
    of the two density images."""
    a0, m0 = cells.min(0), molecules.min(0)
    shape_a = np.ceil((cells.max(0) - a0) / bin_um).astype(int) + 1
    shape_m = np.ceil((molecules.max(0) - m0) / bin_um).astype(int) + 1
    size = shape_a + shape_m

    def image(points, origin, shape):
        h, _, _ = np.histogram2d(points[:, 0], points[:, 1], bins=shape,
                                 range=[[origin[0], origin[0] + shape[0] * bin_um],
                                        [origin[1], origin[1] + shape[1] * bin_um]])
        h = np.log1p(h).astype(np.float32)  # tissue shape, not a few dense spots
        return h - h.mean()

    A = np.fft.rfft2(image(cells, a0, shape_a), s=size)
    M = np.fft.rfft2(image(molecules, m0, shape_m), s=size)
    corr = np.fft.irfft2(A * np.conj(M), s=size)
    k = np.array(np.unravel_index(np.argmax(corr), corr.shape))
    k = np.where(k > size // 2, k - size, k)
    return a0 - m0 + k * bin_um


def fit_shift(cells: np.ndarray, molecules: np.ndarray) -> dict:
    """Coarse cross-correlation, then refine: with cells moved back by d, each
    assigned molecule's nearest cell is offset from it by -(d - d_true) on
    average, so d ← d + median(nearest cell − molecule) until it settles."""
    from scipy.spatial import cKDTree

    d = coarse_shift(cells, molecules, 10.0)
    coarse = d.copy()
    steps = []
    for _ in range(60):
        tree = cKDTree(cells - d)
        dist, idx = tree.query(molecules, workers=-1)
        near = dist < 15.0
        if near.sum() < 100:
            break
        step = np.median((cells[idx[near]] - d) - molecules[near], axis=0)
        d = d + step
        steps.append(float(np.hypot(*step)))
        if steps[-1] < 0.005:
            break
    dist, _ = cKDTree(cells - d).query(molecules, workers=-1)
    # the same fit at a far-off shift, so the number above means something
    baseline, _ = cKDTree(cells - d - 37.0).query(molecules, workers=-1)
    return {"dx": float(d[0]), "dy": float(d[1]),
            "coarse_dx": float(coarse[0]), "coarse_dy": float(coarse[1]), "iterations": len(steps),
            "last_step_um": steps[-1] if steps else None,
            "median_nearest_cell_um": float(np.median(dist)),
            "within_10um": float(np.mean(dist < 10)),
            "median_nearest_cell_um_if_wrong": float(np.median(baseline)),
            "molecules_used": int(len(molecules))}


def fit_zero(cells: np.ndarray, molecules: np.ndarray) -> dict:
    """How well cells and molecules agree as they stand (no shift)."""
    from scipy.spatial import cKDTree

    dist, _ = cKDTree(cells).query(molecules, workers=-1)
    return {"median_nearest_cell_um": float(np.median(dist)), "within_10um": float(np.mean(dist < 10)),
            "molecules_used": int(len(molecules))}


# ── writing ────────────────────────────────────────────────────────────

def split_labels(obs: pd.DataFrame):
    label_cols = [c for c in obs.columns if LABEL_COLUMN.match(c)
                  and (c.split("_")[0] in LABEL_LEVELS or c in ("correlation_coefficient", "bootstrapping_probability"))]
    return obs.drop(columns=label_cols), obs[label_cols]


def write_cell_by_gene(path: Path, X, genes, ids, integer: bool) -> None:
    X = X.tocsr()
    header = True
    for start in range(0, X.shape[0], 20000):
        block = X[start:start + 20000].toarray()
        frame = pd.DataFrame(block.astype(np.int64) if integer else block, columns=genes)
        frame.insert(0, "cell", ids[start:start + 20000])
        frame.to_csv(path, mode="w" if header else "a", header=header, index=False)
        header = False


def plain(series):
    """An obs column as CSV-ready values: categories as text ("" = missing)."""
    if isinstance(series.dtype, pd.CategoricalDtype):
        return series.astype(str)
    return series


def rebuild_slices(cells: dict, sm: Path, out: Path, args, report: dict, missing: list) -> bool:
    obs = cells["obs"]
    sids = obs["_sample_id"].astype(str).to_numpy()
    rest, labels = split_labels(obs.drop(columns=[c for c in obs.columns if c in COMBINE_COLUMNS | COLLAPSED_ID_COLUMNS]))
    bbox = {c: c[-1] for c in ("min_x", "max_x", "min_y", "max_y") if c in rest.columns}
    names = {"gene": "gene", "x": "global_x", "y": "global_y", "z": "global_z", "cell": "cell_id"}
    ok = True
    label_rows = []
    for sid in pd.unique(sids):
        log(f"slice {sid}")
        rows = np.flatnonzero(sids == sid)
        folder = sm / sid
        slice_dir = out / "input" / sid
        entry = {"cells": int(len(rows))}
        if molecule_manifest(folder) is None:
            log(f"  ⚠ no molecules at {folder}: cells written in the combined frame, unshifted")
            missing.append(str(folder))
            d = np.zeros(2)
            entry["shift"] = None
        else:
            suffix = ".csv" if args.molecules_csv else ".parquet"
            mol = write_molecules(folder, slice_dir / f"detected_transcripts{suffix}", names, args.fit_sample, missing)
            fit = fit_shift(cells["spatial"][rows, :2].astype(np.float64), mol["fit"])
            mol.pop("fit")
            entry["molecules"] = mol
            entry["shift"] = fit
            d = np.array([fit["dx"], fit["dy"]])
            good = fit["median_nearest_cell_um"] < FIT_MEDIAN_UM
            ok &= good
            log(f"  shift: combined = raw + ({fit['dx']:.2f}, {fit['dy']:.2f}) µm; median molecule→cell "
                f"{fit['median_nearest_cell_um']:.2f} µm (wrong shift: {fit['median_nearest_cell_um_if_wrong']:.1f}), "
                f"{fit['within_10um']:.0%} within 10 µm {'✓' if good else '✗ DID NOT FIT'}")

        ids = np.arange(1, len(rows) + 1)
        meta = pd.DataFrame({"EntityID": ids})
        xy = cells["spatial"][rows].astype(np.float64)
        meta["center_x"] = np.round(xy[:, 0] - d[0], 2)
        meta["center_y"] = np.round(xy[:, 1] - d[1], 2)
        if xy.shape[1] > 2:
            meta["center_z"] = xy[:, 2]
        for col in rest.columns:
            if col in ("center_x", "center_y", "center_z"):
                continue
            values = plain(rest[col].iloc[rows]).to_numpy()
            if col in bbox:
                values = np.round(values.astype(np.float64) - d[0 if bbox[col] == "x" else 1], 2)
            meta[col] = values
        slice_dir.mkdir(parents=True, exist_ok=True)
        meta.to_csv(slice_dir / "cell_metadata.csv", index=False)
        write_cell_by_gene(slice_dir / "cell_by_gene.csv", cells["X"][rows], cells["genes"], ids, cells["integer"])
        log(f"  cells: {len(rows):,} → cell_metadata.csv, cell_by_gene.csv")
        if len(labels.columns):
            part = labels.iloc[rows].apply(plain)
            part.insert(0, "cell_id", [f"{sid}_{i}" for i in ids])
            label_rows.append(part)
        report["slices"][sid] = entry
    if label_rows:
        write_labels(pd.concat(label_rows), out, report)
    return ok


def rebuild_single(cells: dict, sm: Path, out: Path, args, report: dict, missing: list) -> bool:
    import anndata as ad

    rest, labels = split_labels(cells["obs"])
    ids = np.array([str(i) for i in range(1, len(rest) + 1)])
    rest = rest.copy()
    rest.index = ids
    adata = ad.AnnData(X=cells["X"].tocsr(), obs=rest, var=pd.DataFrame(index=pd.Index(cells["genes"])))
    adata.obsm["spatial"] = cells["spatial"].astype(np.float64)
    for name, emb in cells["embeddings"].items():
        adata.obsm[f"X_{name}"] = emb
    (out / "input").mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out / "input" / "cell_by_gene.h5ad", compression="gzip")
    log(f"cells: {adata.n_obs:,} → cell_by_gene.h5ad")

    entry = {"cells": int(adata.n_obs)}
    ok = True
    if molecule_manifest(sm) is None:
        log(f"  ⚠ no molecules at {sm}")
        missing.append(str(sm))
    else:
        names = {"gene": "gene_names", "x": "x", "y": "y", "z": "z", "cell": "cell_ids"}
        suffix = ".csv" if args.molecules_csv else ".parquet"
        mol = write_molecules(sm, out / "input" / f"segmented_spot_table{suffix}", names, args.fit_sample, missing)
        fit = fit_zero(cells["spatial"][:, :2].astype(np.float64), mol.pop("fit"))
        entry["molecules"] = mol
        entry["alignment"] = fit
        ok = fit["median_nearest_cell_um"] < FIT_MEDIAN_UM
        log(f"  cells and molecules as stored: median molecule→cell {fit['median_nearest_cell_um']:.2f} µm "
            f"{'✓' if ok else '✗ they do not line up'}")
    report["slices"]["(single)"] = entry
    if len(labels.columns):
        part = labels.apply(plain)
        part.insert(0, "cell_id", ids)
        write_labels(part, out, report)
    return ok


def write_labels(labels: pd.DataFrame, out: Path, report: dict) -> None:
    (out / "labels").mkdir(parents=True, exist_ok=True)
    labels.to_csv(out / "labels" / "mapping_output.csv", index=False)
    report["labels"] = [c for c in labels.columns if c != "cell_id"]
    log(f"labels: {len(report['labels'])} MapMyCells columns → labels/mapping_output.csv")


# ── main ───────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("sample_dir", type=Path, help="the downloaded <sample>/ holding meyes_output/ and sm_output/")
    ap.add_argument("out", type=Path, help="where to write (must not exist or be empty)")
    ap.add_argument("--molecules-csv", action="store_true", help="write molecules as CSV instead of parquet (much larger)")
    ap.add_argument("--fit-sample", type=int, default=2_000_000, help="assigned molecules per slice used to measure its shift")
    args = ap.parse_args(argv)

    meyes, sm = args.sample_dir / "meyes_output", args.sample_dir / "sm_output"
    if not (meyes / "manifest.json").exists():
        print(f"{meyes}/manifest.json not found: not a legacy sample folder", file=sys.stderr)
        return 2
    if args.out.exists() and any(args.out.iterdir()):
        print(f"{args.out} is not empty", file=sys.stderr)
        return 2
    try:
        args.out.resolve().relative_to(args.sample_dir.resolve())
        print("OUT must not be inside SAMPLE_DIR", file=sys.stderr)
        return 2
    except ValueError:
        pass

    missing: list = []
    try:
        cells = read_cells(meyes, missing)
    except NotLegacy as e:
        print(e, file=sys.stderr)
        return 2
    multi = "_sample_id" in cells["obs"].columns and molecule_manifest(sm) is None
    layout = "slices" if multi else "single"
    log(f"layout: {layout}")
    report = {"sample": args.sample_dir.resolve().name, "layout": layout,
              "cells": int(cells["X"].shape[0]), "genes": len(cells["genes"]),
              "integer_counts": cells["integer"], "slices": {},
              "source": {"meyes_output": cells["manifest"]}}

    ok = (rebuild_slices if multi else rebuild_single)(cells, sm, args.out, args, report, missing)
    report["missing"] = missing
    report["shift_fit_threshold_um"] = FIT_MEDIAN_UM
    report["lost"] = [
        "original cell ids (EntityID stored as float32): cells renumbered 1..n per slice",
        "molecule→cell assignment: cell_id is 1 (in a cell) / -1 (in none), "
        "or no cell column when the old run never split them",
        "cells removed by the old artifact mask",
        "control probes and blanks",
    ] + (["per-slice combine shift: measured from the molecules, see slices.*.shift"] if multi else [])
    (args.out / "reconstruction.json").write_text(json.dumps(report, indent=2, default=str))

    name = report["sample"]
    print()
    if multi:
        print("Next — the new pipeline on the rebuilt slices:")
        print(f"  tests/bil/run_dry_run.sh {args.out}/input {args.out}/run {name.replace('-', '_')} \\")
        print(f"      --mask-percentile 0" + (f" --labels {args.out}/labels/mapping_output.csv" if report.get("labels") else " --reference …"))
    else:
        print("Next — prepare.py on the rebuilt table and molecules:")
        print(f"  python scripts/prepare.py {args.out}/input {args.out}/store --name {name.replace('-', '_')}"
              + (f" --labels {args.out}/labels/mapping_output.csv" if report.get("labels") else ""))
    if missing:
        log(f"⚠ INCOMPLETE: {len(missing)} files missing from the download (reconstruction.json → missing)")
        return 3
    if not ok:
        log("⚠ a slice's shift did not fit — check reconstruction.json before using it")
        return 4
    log("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
