#!/usr/bin/env python3
"""Turn a raw MERFISH run into one SpatialData store MERFISHEYES opens as it is.

  python scripts/prepare_merfish_run.py RUN FOLDER [--name NAME] [--sets 1,2,9]
                                        [--labels-from H5AD] [--keep-z]
                                        [--workers N] [--resume]

RUN is a run in the baseline shape (docs/specs/raw-merfish-profile.md):
final_spots/ (per set: XF, genes, …; de-duplicated and score-filtered) and
final_mats/ or final_matsWithGFP/ (per field: Mcells, Mgenes, Mcounts, Xum).
Writes

  FOLDER/NAME.zarr/   tables/table: one row per cell after the field-overlap
                      rule (counts × the panel's genes, X as CSC; obs section,
                      field, n_molecules, n_blank, gfp_count, centroid_z;
                      obsm spatial — the sections side by side, y
                      flipped so each reads as in the lab's figures);
                      points/molecules: every molecule but the blanks, tiled,
                      with its section and field; DE stats for the table's
                      categorical columns
  FOLDER/mapping.json names the store

Each set is a section, on a grid read left to right, top to bottom (as
scripts/combine_slices_v3.py lays out slices). A set's stitched frame has y
going down; the store flips it (the viewer draws +y upwards) so a section
looks as in the lab's own plots: x = x_set + dx, y = dy - y_set, with each
set's dx, dy in the table's uns['merfisheyes_sections']. No rotation.

--labels-from adds an .h5ad's obs columns and embeddings (obsm X_umap, …) to
the cells with the same names (her MapMyCells labels, leiden, UMAP). With
--resume, sets whose molecules an interrupted run already wrote are kept.

scripts/prepare.py runs this for a run folder. RUN is never changed.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import points_zarr as pz  # noqa: E402
import prepare_folder  # noqa: E402
import raw_merfish as rm  # noqa: E402
from meyes_common import should_filter_gene  # noqa: E402

PADDING_UM = 250.0  # between sections on the grid
CHUNK_ROWS = 4_000_000  # molecules read at a time: ~300 MB with their names
PROGRESS = "progress.json"  # in the molecules' scratch folder: the sets written
# The obsm tables taken from an .h5ad (--labels-from): embeddings by name
EMBEDDING = re.compile(r"^X_(umap|tsne|phate|mde|draw_graph|diffmap|fa)(?![a-z])", re.I)


# ── sections ─────────────────────────────────────────────────────────────

def place(x, y, sh: dict):
    """A set's own (stitched, y-down) µm into the store's frame: y flipped,
    so a section reads as in the lab's own figures (the viewer draws +y
    upwards), and shifted onto its grid slot. The one rule cells, molecules
    and extents share; back again: x - dx, dy - y."""
    return x + sh["dx"], sh["dy"] - y


def layout(sets: list[rm.SpotSet], padding: float = PADDING_UM) -> dict[str, dict]:
    """Each set's slot on a grid: ceil(√n) columns, each as wide as its
    widest set plus `padding`, rows likewise, the first row at the top and
    rows going down (combine_slices_v3.py's grid: rows along -y). Returns
    the shifts `place` applies."""
    n = len(sets)
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    boxes = [s.extent_um() for s in sets]
    w = [b[2] - b[0] for b in boxes]
    h = [b[3] - b[1] for b in boxes]
    col_w = [max(w[i] for i in range(n) if i % cols == c) + padding for c in range(cols)]
    row_h = [max(h[i] for i in range(n) if i // cols == r) + padding for r in range(rows)]
    out = {}
    for i, (s, b) in enumerate(zip(sets, boxes)):
        left = sum(col_w[: i % cols])
        top = -sum(row_h[: i // cols])
        # the set's left edge (raw x0) at `left`, its top (raw y0, flipped) at `top`
        out[s.id] = {"section": s.id, "dx": left - b[0], "dy": top + b[1], "flip_y": True,
                     "extent_um": [round(b[2] - b[0], 2), round(b[3] - b[1], 2)],
                     "grid": [i // cols, i % cols]}
    return out


# ── cells ────────────────────────────────────────────────────────────────

def build_table(run: rm.Run, sets: list[rm.SpotSet], shifts: dict, labels_from: Path | None, log):
    """The cells as an AnnData: every field's cells after the overlap rule."""
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    run_sets = rm.Run(run.root, sets, run.mats, run.gfp_fields, run.notes)
    panel, fields, problems = rm.read_cells(run_sets, progress=log)
    for p in problems:
        log(f"WARN {p}")
    if not fields:
        return None
    blank = np.array([should_filter_gene(g) for g in panel])
    names, section, field, counts, xy, z, gfp = [], [], [], [], [], [], []
    any_gfp = any(f.gfp is not None for f in fields)
    for f in fields:
        k = f.keep
        if not k.any():
            continue
        sh = shifts[f.set]
        names.append(f.cells[k])
        section.append(np.full(int(k.sum()), f.set, dtype=object))
        field.append(np.full(int(k.sum()), f.name, dtype=object))
        counts.append(sp.csr_matrix(f.counts[k].astype(np.float32)))
        px, py = place(f.xum[k][:, 1], f.xum[k][:, 2], sh)
        xy.append(np.c_[px, py])
        z.append(f.xum[k][:, 0])
        if any_gfp:
            gfp.append(f.gfp[k].astype(np.float64) if f.gfp is not None else np.full(int(k.sum()), np.nan))
    X = sp.vstack(counts).tocsr()
    obs = pd.DataFrame(index=pd.Index(np.concatenate(names).astype(str), name="cell_id"))
    obs["section"] = pd.Categorical(np.concatenate(section), categories=[s.id for s in sets])
    obs["field"] = pd.Categorical(np.concatenate(field))
    obs["n_molecules"] = np.asarray(X[:, ~blank].sum(1)).ravel().astype(np.int64)
    obs["n_blank"] = np.asarray(X[:, blank].sum(1)).ravel().astype(np.int64)
    if any_gfp:
        obs["gfp_count"] = np.concatenate(gfp)
    obs["centroid_z"] = np.concatenate(z)
    if not obs.index.is_unique:
        dup = obs.index[obs.index.duplicated()][:5].tolist()
        raise SystemExit(f"cell names repeat across fields ({dup}): Mcells must be unique in a run")
    adata = ad.AnnData(
        X=sp.csc_matrix(X), obs=obs, var=pd.DataFrame(index=pd.Index(panel, name="gene")),
        obsm={"spatial": np.concatenate(xy)},
    )
    adata.uns["merfisheyes_sections"] = json.dumps(list(shifts.values()))
    adata.uns["merfisheyes_source"] = json.dumps({
        "kind": "merfish_run", "run": str(run.root),
        "cells_in_mats": int(sum(len(f.cells) for f in fields)),
        "rule": "a cell is kept in the field whose centre is nearest (Chebyshev); exact ties to the lower index",
    })
    log(f"cells: {adata.n_obs:,} after the overlap rule (of {sum(len(f.cells) for f in fields):,} in the mats),"
        f" {int((~blank).sum())} genes + {int(blank.sum())} blanks counted as n_blank"
        + (f", gfp_count from {run.gfp_fields} fields" if any_gfp else ""))
    if labels_from:
        add_labels(adata, labels_from, log)
    return adata


def add_labels(adata, h5ad: Path, log) -> None:
    """obs columns and embeddings of an .h5ad, matched by cell name. A cell
    the file lacks gets "" (a label) or NaN (a number, an embedding)."""
    import anndata as ad
    import pandas as pd

    src = ad.read_h5ad(h5ad, backed="r")
    try:
        found = adata.obs_names.isin(src.obs_names)
        log(f"labels from {h5ad.name}: {int(found.sum()):,} of {adata.n_obs:,} cells matched by name"
            f" ({src.n_obs:,} in the file)")
        if not found.any():
            log(f"WARN {h5ad.name}: no cell names in common (its first: {list(src.obs_names[:3])})")
            return
        obs = src.obs.reindex(adata.obs_names)
        for col in src.obs.columns:
            name = col if col not in adata.obs.columns else f"{col}_labels"
            v = obs[col]
            if isinstance(src.obs[col].dtype, pd.CategoricalDtype) or v.dtype == object:
                cats = [str(c) for c in (src.obs[col].cat.categories if hasattr(src.obs[col], "cat")
                                         else pd.unique(src.obs[col].astype(str)))]
                adata.obs[name] = pd.Categorical(v.astype(object).where(v.notna(), "").astype(str),
                                                 categories=sorted(set(cats) | {""}))
            elif np.issubdtype(v.dtype, np.number) or v.dtype == bool:
                adata.obs[name] = v.astype(np.float64)
        rows = src.obs_names.get_indexer(adata.obs_names)
        for key in src.obsm.keys():
            # embeddings only: a table with another name may be a spatial
            # layout of its own (her X_shifted lays the sets out differently
            # from this store, and the viewer would offer it beside the cells)
            if not EMBEDDING.search(key):
                log(f"  obsm {key} left out (not an embedding)")
                continue
            m = np.asarray(src.obsm[key])
            if m.ndim != 2 or m.shape[1] > 3:
                continue
            out = np.full((adata.n_obs, m.shape[1]), np.nan, np.float32)
            out[rows >= 0] = m[rows[rows >= 0]]
            adata.obsm[key] = out
            log(f"  obsm {key} added")
    finally:
        src.file.close()


# ── molecules ────────────────────────────────────────────────────────────

def _spill_set(job: dict) -> dict:
    """One set's molecules into the writer's band files (a worker process)."""
    s: rm.SpotSet = job["set"]
    w: pz.PointsWriter = job["writer"]
    codes = rm.GeneCodes(job["panel"])
    kept_of = job["kept_of"]  # panel code -> feature code, -1 = left out
    sh, keep_z = job["shift"], job["keep_z"]
    field0, sec = job["field0"], job["section_code"]
    xf, g = s.xf(), s.genes()
    n_kept = n_left = n_unknown = 0
    lo, hi = [np.inf] * 2, [-np.inf] * 2
    for ci, a in enumerate(range(0, s.n, CHUNK_ROWS)):
        c = codes.codes(g.rows(a, CHUNK_ROWS))
        known = c < len(kept_of)
        f = np.full(len(c), -1, np.int64)
        f[known] = kept_of[c[known]]
        n_unknown += int((~known).sum())
        keep = f >= 0
        n_left += int(known.sum() - keep.sum())
        rows = xf.rows(a, CHUNK_ROWS)[keep]
        if len(rows):
            x, y = place(rows[:, rm.XF_X].astype(np.float64) * s.xy_um,
                         rows[:, rm.XF_Y].astype(np.float64) * s.xy_um, sh)
            spill = {"x": x, "y": y, "f": f[keep].astype(w.fdtype),
                     "section": np.full(len(rows), sec, np.uint8),
                     "field": (rows[:, rm.XF_FIELD].astype(np.int64) + field0).astype(np.uint16)}
            if keep_z:
                spill["z"] = rows[:, rm.XF_Z].astype(np.float64) * float(s.pix[0])
            w.spill(spill, seed=job["seed"] * 100_000 + ci, tag=f"set{s.id}c{ci}")
            n_kept += len(rows)
            lo = [min(lo[0], float(x.min())), min(lo[1], float(y.min()))]
            hi = [max(hi[0], float(x.max())), max(hi[1], float(y.max()))]
        print(f"set {s.id}: {min(a + CHUNK_ROWS, s.n):,} / {s.n:,} molecules", flush=True)
    return {"set": s.id, "molecules": s.n, "written": n_kept, "left_out": n_left,
            "unknown": n_unknown, "min": lo, "max": hi}


def write_molecules(store: Path, sets: list[rm.SpotSet], shifts: dict, panel: list[str], *,
                    keep_z: bool, workers: int, resume: bool, log) -> dict:
    """Every set's molecules (blanks left out) into points/molecules."""
    import zarr

    kept = [g for g in panel if not should_filter_gene(g)]
    kept_of = np.array([kept.index(g) if g in kept else -1 for g in panel], np.int64)
    fields = [f for s in sets for f in s.fovs]
    if len(fields) >= 2**16:
        raise SystemExit(f"{len(fields)} fields: more than the field column holds")
    field0, k = {}, 0
    for s in sets:
        field0[s.id] = k
        k += len(s.fovs)
    boxes = [(s.extent_um(), shifts[s.id]) for s in sets]
    corners = [place(np.array([b[0], b[2]]), np.array([b[1], b[3]]), sh) for b, sh in boxes]
    x0 = min(float(cx.min()) for cx, _ in corners)
    x1 = max(float(cx.max()) for cx, _ in corners)
    y0 = min(float(cy.min()) for _, cy in corners)
    y1 = max(float(cy.max()) for _, cy in corners)
    pad = 1.0
    root = zarr.open_group(store, mode="a")
    root.attrs.setdefault("spatialdata_attrs", {"version": "0.2"})
    root.require_group("points")  # its own zarr.json, or consolidation never finds the element
    out = store / "points" / "molecules"
    tmp = store.parent / ".molecules.tmp"
    done_file = tmp / PROGRESS
    done: dict = json.loads(done_file.read_text()) if resume and done_file.exists() else {}
    w = pz.PointsWriter(
        out, tmp, features=kept, dims=3 if keep_z else 2,
        origin=(x0 - pad, y0 - pad), extent=(x1 - x0 + 2 * pad, y1 - y0 + 2 * pad),
        columns=[("section", "u1", [s.id for s in sets]), ("field", "u2", fields)],
        copied=["section"], keep_tmp=bool(done),
        band_um=pz.band_for(sum(s.n for s in sets), y1 - y0 + 2 * pad),
    )
    # a set an interrupted run left half-written is written again
    if done:
        for f in (tmp / "L0").glob("*.bin") if (tmp / "L0").exists() else []:
            tag = f.name.split(".")[1]
            if tag.split("c")[0].removeprefix("set") not in done:
                f.unlink()
        log(f"molecules: resuming — sets {', '.join(done)} were written before")
    todo = [s for s in sets if s.id not in done]
    total = sum(s.n for s in sets)
    log(f"molecules: {total:,} in {len(sets)} sets ({len(todo)} to read), {len(kept)} genes"
        f" ({len(panel) - len(kept)} blanks left out), {'3' if keep_z else '2'}D,"
        f" extent {x1 - x0:,.0f} × {y1 - y0:,.0f} µm, bands of {w.band_um:.0f} µm")
    tmp.mkdir(parents=True, exist_ok=True)
    jobs = [{"set": s, "writer": w, "panel": panel, "kept_of": kept_of, "shift": shifts[s.id],
             "keep_z": keep_z, "field0": field0[s.id], "section_code": sets.index(s), "seed": sets.index(s) + 1}
            for s in todo]
    t0 = time.time()

    def record(r: dict) -> None:
        done[r["set"]] = r
        done_file.write_text(json.dumps(done))
        log(f"set {r['set']} written: {r['written']:,} molecules ({len(done)} / {len(sets)} sets,"
            f" {time.time() - t0:,.0f} s)")

    if workers <= 1:  # in this process: one set at a time
        for job in jobs:
            record(_spill_set(job))
    else:
        with ProcessPoolExecutor(min(workers, len(jobs) or 1)) as ex:
            for r in ex.map(_spill_set, jobs):
                record(r)
    unknown = sum(r["unknown"] for r in done.values())
    if unknown:
        log(f"WARN {unknown:,} molecules name a gene outside the panel (Mgenes): left out")
    log("molecules: building the levels ...")
    w.finish()
    shutil.rmtree(tmp, ignore_errors=True)
    zarr.consolidate_metadata(str(store))
    return {"molecules": total, "written": sum(r["written"] for r in done.values()),
            "blanks_left_out": sum(r["left_out"] for r in done.values()), "unknown": unknown,
            "seconds": round(time.time() - t0, 1)}


# ── the whole run ────────────────────────────────────────────────────────

def prepare_merfish_run(source: Path, folder: Path, *, name: str | None = None,
                        sets: list[str] | None = None, labels_from: Path | None = None,
                        keep_z: bool = False, workers: int = 4, cluster: bool = False,
                        label_csv=None, resume: bool = False, log=print):
    try:
        run = rm.detect(source)
    except rm.NotARun as e:
        raise SystemExit(str(e))
    chosen = run.sets
    if sets:
        want = [str(s) for s in sets]
        missing = [s for s in want if s not in {x.id for x in run.sets}]
        if missing:
            raise SystemExit(f"no set {', '.join(missing)} in {run.root.name} (sets: "
                             f"{', '.join(x.id for x in run.sets)})")
        chosen = [s for s in run.sets if s.id in want]
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    name = name or run.root.name
    store = folder / f"{name}.zarr"
    log(f"{run.root.name}: {len(chosen)} of {len(run.sets)} sets, "
        f"{sum(s.n for s in chosen):,} molecules → {store.name}")
    shifts = layout(chosen)
    for line in rm.undecoded_fields(run, chosen):
        log(f"WARN {line}")

    table_done = store.exists() and (store / "tables" / "table").exists() and resume
    if not table_done:
        if store.exists():
            shutil.rmtree(store)
        adata = build_table(run, chosen, shifts, labels_from, log)
        if adata is None:
            raise SystemExit(f"{run.root.name}: no cell mats for these sets — molecules alone are"
                             " not converted yet")
        prepare_folder.write_store_from_adata(adata, store, log)
        del adata
    else:
        log("cells: written before (resuming)")
    panel = rm.panel_of(run)
    info = write_molecules(store, chosen, shifts, panel, keep_z=keep_z, workers=workers,
                           resume=resume, log=log)
    log(f"molecules: {info['written']:,} written ({info['blanks_left_out']:,} blanks left out) in"
        f" {info['seconds']:,.0f} s")
    return prepare_folder.prepare(folder, cells=store, cluster=cluster, labels=label_csv, log=log)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run", type=Path, help="the run folder (holding final_spots/)")
    ap.add_argument("folder", type=Path, help="the folder to write")
    ap.add_argument("--name", help="the store's name (default: the run folder's)")
    ap.add_argument("--sets", help="only these sets, e.g. 1,2,9")
    ap.add_argument("--labels-from", type=Path, help="an .h5ad of the same run: its labels and embeddings")
    ap.add_argument("--keep-z", action="store_true", help="keep the molecules' z (3-D)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--resume", action="store_true", help="keep what an interrupted run wrote")
    a = ap.parse_args(argv)
    r = prepare_merfish_run(a.run, a.folder, name=a.name, sets=a.sets.split(",") if a.sets else None,
                            labels_from=a.labels_from, keep_z=a.keep_z, workers=a.workers, resume=a.resume)
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
