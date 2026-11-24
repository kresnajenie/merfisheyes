#!/usr/bin/env python3
"""
Analysis steps an agent runs on a dataset open in MERFISHEYES (MCP.md
pillar D, tier 2): run by the @merfisheyes/mcp package in its own Python, one
step per process, as a background job.

    python meyes_analysis.py STEP --store S.zarr --work WORKDIR --args ARGS.json --out RESULT.json

STEP is one of:
  profile    what the table holds: columns and their values (with counts, per
             section when asked), embeddings, the spacing of the cells (their
             units), which marker genes the panel measures, stored DE stats
  niche      the cells within a radius of anchor cells (by label), per
             section; a label column `niche_<name>`
  recluster  for a subset of cells: normalise, PCA, neighbours, UMAP and
             Leiden at several resolutions; columns `<name>_leiden_<res>`, an
             embedding `X_umap_<name>`, and each cluster's top genes
  markers    top genes per value of a label column (the store's DE stats when
             it has them for that column, else computed the same way), and a
             heatmap of chosen genes
  save       a copy of the dataset with the workspace's results: .zarr (a
             SpatialData store, with DE stats for the new labels) or .h5ad

Results of niche and recluster go into the workspace (WORKDIR) in the table's
cell order — `obs__<name>.parquet` (one column), `obsm__<name>.npy` (NaN for
cells not covered) — and are listed in WORKDIR/workspace.json. The package
shows them in the open view at once; `save` writes them into a copy. The
store is never changed.

Expression is log1p(normalize_total(X)) for everything (de_stats.log_normalize),
the numbers the viewer's DEG panel shows.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

T0 = time.time()


def say(msg: str) -> None:
    """A progress line the package relays to the agent."""
    print(f"[{time.time() - T0:5.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------- the table


class Table:
    """The default table of a SpatialData store, read piece by piece (obs,
    var, obsm, X only when needed): a store of millions of cells is not
    loaded whole for a profile."""

    def __init__(self, store: Path):
        import zarr
        from validate_store import default_table

        self.store = Path(store)
        root = zarr.open_group(str(self.store), mode="r", use_consolidated=False)
        tables = sorted(root["tables"].group_keys()) if "tables" in root else []
        if not tables:
            raise SystemExit(f"{self.store} has no table.")
        self.name = default_table(dict(root.attrs), tables) or tables[0]
        self.group = root[f"tables/{self.name}"]
        self._obs = None
        self._x = None

    @property
    def obs(self) -> pd.DataFrame:
        import anndata as ad

        if self._obs is None:
            self._obs = ad.io.read_elem(self.group["obs"])
        return self._obs

    @property
    def genes(self) -> list[str]:
        import anndata as ad

        return [str(g) for g in ad.io.read_elem(self.group["var"]).index]

    def obsm_keys(self) -> list[str]:
        return sorted(self.group["obsm"].array_keys()) if "obsm" in self.group else []

    def coords(self) -> np.ndarray:
        keys = self.obsm_keys()
        key = "spatial" if "spatial" in keys else "X_spatial" if "X_spatial" in keys else None
        if key is None:
            raise SystemExit("The table has no obsm/spatial coordinates.")
        return np.asarray(self.group["obsm"][key][:])[:, :2].astype(np.float64)

    def counts(self):
        """X as CSR float32 (counts as stored)."""
        import anndata as ad
        import scipy.sparse as sp

        if self._x is None:
            self._x = sp.csr_matrix(ad.io.read_elem(self.group["X"])).astype(np.float32)
        return self._x

    def stored_de(self) -> list[str]:
        if "uns" not in self.group or "merfisheyes_de" not in self.group["uns"]:
            return []
        return sorted(self.group["uns"]["merfisheyes_de"].group_keys())


class Workspace:
    """An agent's results for one dataset, in the table's cell order."""

    def __init__(self, work: Path, n_cells: int):
        self.dir = Path(work)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.n = n_cells
        self.index_file = self.dir / "workspace.json"
        self.index = json.loads(self.index_file.read_text()) if self.index_file.exists() else {"obs": {}, "obsm": {}}

    def put_obs(self, name: str, values, how: dict) -> None:
        col = pd.Categorical(values)
        pd.DataFrame({name: col}).to_parquet(self.dir / f"obs__{name}.parquet")
        self.index["obs"][name] = {"how": how, "values": [str(c) for c in col.categories][:200]}
        self._write()

    def put_obsm(self, name: str, array: np.ndarray, how: dict) -> None:
        np.save(self.dir / f"obsm__{name}.npy", array.astype(np.float32))
        self.index["obsm"][name] = {"how": how, "shape": list(array.shape)}
        self._write()

    def obs(self, name: str) -> pd.Series:
        return pd.read_parquet(self.dir / f"obs__{name}.parquet")[name]

    def obsm(self, name: str) -> np.ndarray:
        return np.load(self.dir / f"obsm__{name}.npy")

    def _write(self) -> None:
        self.index_file.write_text(json.dumps(self.index, indent=2))


def column(table: Table, ws: Workspace, name: str) -> pd.Series:
    """A label column: the workspace's own first (an agent's results), else
    the table's."""
    if name in ws.index["obs"]:
        return ws.obs(name).astype(str)
    if name in table.obs.columns:
        return table.obs[name].astype(str).reset_index(drop=True)
    have = sorted(set(table.obs.columns) | set(ws.index["obs"]))
    raise SystemExit(f'No column "{name}". Columns: {", ".join(have)}.')


def pick(table: Table, ws: Workspace, spec: dict | None) -> np.ndarray:
    """A boolean mask from {column, values} (all cells when absent)."""
    if not spec:
        return np.ones(len(table.obs), dtype=bool)
    col = column(table, ws, spec["column"])
    values = [str(v) for v in spec.get("values", [])]
    missing = [v for v in values if v not in set(col.unique())]
    if missing:
        near = [v for v in col.unique() if any(w.lower() in v.lower() for m in missing for w in m.split() if len(w) >= 3)]
        raise SystemExit(
            f'"{spec["column"]}" has no value {", ".join(map(repr, missing))}.'
            + (f" Values containing those words: {', '.join(map(repr, near[:12]))}." if near else "")
        )
    return np.array(col.isin(values).to_numpy(), dtype=bool)  # writable


def log_normalized(counts):
    from de_stats import log_normalize

    return log_normalize(counts)


def group_stats(lognorm, counts, codes: np.ndarray, n_groups: int):
    """Mean log-normalised expression and fraction of cells expressing, per
    group (genes × groups), as de_stats writes them."""
    import scipy.sparse as sp

    ok = codes >= 0
    membership = sp.csr_matrix(
        (np.ones(int(ok.sum()), dtype=np.float32), (np.flatnonzero(ok), codes[ok])),
        shape=(len(codes), n_groups),
    )
    sizes = np.bincount(codes[ok], minlength=n_groups).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.nan_to_num(np.asarray((lognorm.T @ membership).todense()) / sizes)
        pct = np.nan_to_num(np.asarray(((counts > 0).astype(np.float32).T @ membership).todense()) / sizes)
    return mean, pct, sizes


def top_genes(mean: np.ndarray, pct: np.ndarray, sizes: np.ndarray, genes: list[str], groups: list[str], n: int) -> dict:
    """Each group's top genes: mean in the group minus the mean over the other
    cells (weighted by group size), then the fraction expressing."""
    total = (mean * sizes).sum(axis=1)
    out = {}
    for j, g in enumerate(groups):
        rest_n = sizes.sum() - sizes[j]
        rest = (total - mean[:, j] * sizes[j]) / rest_n if rest_n > 0 else np.zeros(len(genes))
        diff = mean[:, j] - rest
        order = np.argsort(-diff)[:n]
        out[g] = [
            {"gene": genes[i], "mean": round(float(mean[i, j]), 3), "rest": round(float(rest[i]), 3),
             "pct": round(float(pct[i, j]), 3)}
            for i in order
        ]
    return out


# -------------------------------------------------------------------- steps


def profile(table: Table, ws: Workspace, a: dict) -> dict:
    obs = table.obs
    max_values = int(a.get("max_values", 25))
    say(f"{len(obs):,} cells, {len(table.genes):,} genes")
    columns = []
    for name in obs.columns:
        s = obs[name]
        numeric = pd.api.types.is_numeric_dtype(s) and not isinstance(s.dtype, pd.CategoricalDtype)
        n_unique = int(s.nunique(dropna=True))
        entry = {"name": name, "kind": "numeric" if numeric else "categorical", "n_values": n_unique}
        if numeric:
            entry.update(min=float(np.nanmin(s)), max=float(np.nanmax(s)))
        if not numeric or n_unique <= 50:
            vc = s.astype(str).value_counts()
            entry["values"] = [{"value": v, "cells": int(c)} for v, c in vc.head(max_values).items()]
        columns.append(entry)
    out = {"cells": len(obs), "genes": len(table.genes), "table": table.name, "columns": columns,
           "obsm": table.obsm_keys(), "stored_de_columns": table.stored_de(),
           "workspace": {"obs": list(ws.index["obs"]), "obsm": list(ws.index["obsm"])}}

    # The cells' spacing: what units the coordinates are in
    xy = table.coords()
    from scipy.spatial import cKDTree

    rng = np.random.default_rng(0)
    sample = xy[rng.choice(len(xy), size=min(5000, len(xy)), replace=False)]
    d, _ = cKDTree(xy).query(sample, k=2)
    spacing = float(np.median(d[:, 1]))
    out["coordinates"] = {
        "median_nearest_cell": round(spacing, 3),
        "extent": [round(float(v), 1) for v in (xy[:, 0].min(), xy[:, 0].max(), xy[:, 1].min(), xy[:, 1].max())],
        "units_guess": "microns" if 3 <= spacing <= 40 else "not microns (normalised or pixels?)",
    }

    if a.get("breakdown"):
        b = a["breakdown"]
        rows = column(table, ws, b["column"])
        by = column(table, ws, b["by"])
        ct = pd.crosstab(rows, by)
        if b.get("values"):
            ct = ct.loc[[v for v in b["values"] if v in ct.index]]
        out["breakdown"] = {"column": b["column"], "by": b["by"], "counts": {r: {c: int(v) for c, v in row.items()} for r, row in ct.head(200).iterrows()}}

    if a.get("genes"):
        genes = table.genes
        lower = {g.lower(): i for i, g in enumerate(genes)}
        want = [g for g in a["genes"]]
        found = [(g, lower.get(g.lower())) for g in want]
        idx = [i for _, i in found if i is not None]
        detect = {}
        if idx:
            x = table.counts()[:, idx]
            frac = np.asarray((x > 0).mean(axis=0)).ravel()
            detect = {genes[i]: round(float(f), 4) for i, f in zip(idx, frac)}
        out["markers"] = {
            "measured": [{"gene": genes[i], "fraction_of_cells": detect.get(genes[i])} for g, i in found if i is not None],
            "not_in_panel": [g for g, i in found if i is None],
        }
    return out


def niche(table: Table, ws: Workspace, a: dict) -> dict:
    from scipy.spatial import cKDTree

    name = a["name"]
    radius = float(a["radius_um"])
    min_anchors = int(a.get("min_anchors", 1))
    xy = table.coords()
    anchors = pick(table, ws, a["anchors"])
    excluded = pick(table, ws, a["exclude"]) if a.get("exclude") else np.zeros(len(xy), bool)
    anchors &= ~excluded
    sections = column(table, ws, a["within"]).to_numpy() if a.get("within") else np.zeros(len(xy), dtype=object)
    inside = np.zeros(len(xy), dtype=bool)
    per_section = {}
    for sec in pd.unique(sections):
        here = sections == sec
        a_idx = np.flatnonzero(here & anchors)
        if len(a_idx) == 0:
            per_section[str(sec)] = {"anchors": 0, "cells": 0}
            continue
        tree = cKDTree(xy[a_idx])
        cells = np.flatnonzero(here & ~excluded)
        if min_anchors <= 1:
            d, _ = tree.query(xy[cells], distance_upper_bound=radius)
            hit = np.isfinite(d)
        else:
            hit = np.array([len(n) >= min_anchors for n in tree.query_ball_point(xy[cells], radius)], dtype=bool)
        inside[cells[hit]] = True
        per_section[str(sec)] = {"anchors": int(len(a_idx)), "cells": int(hit.sum())}
    say(f"niche {name}: {int(inside.sum()):,} cells within {radius:g} µm of {int(anchors.sum()):,} anchors")
    col = f"niche_{name}"
    ws.put_obs(col, np.where(inside, name, ""), {"step": "niche", **a})
    out = {"column": col, "cells": int(inside.sum()), "anchors": int(anchors.sum()), "excluded": int(excluded.sum()),
           "per_section": per_section if a.get("within") else None}
    if a.get("describe_by"):
        lab = column(table, ws, a["describe_by"])[inside].value_counts()
        out["composition"] = {"column": a["describe_by"], "top": [{"value": v, "cells": int(c)} for v, c in lab.head(25).items()]}
    return out


def recluster(table: Table, ws: Workspace, a: dict) -> dict:
    import anndata as ad
    import scanpy as sc

    name = a["name"]
    mask = pick(table, ws, a["cells"])
    n = int(mask.sum())
    if n < 20:
        raise SystemExit(f"Only {n} cells selected: too few to cluster.")
    say(f"{n:,} cells selected")
    counts = table.counts()[mask]
    genes = table.genes
    lognorm = log_normalized(counts)
    sub = ad.AnnData(lognorm.copy())
    n_comps = max(1, min(int(a.get("n_pcs", 30)), n - 1, len(genes) - 1))
    sc.pp.pca(sub, n_comps=n_comps)
    sc.pp.neighbors(sub, n_neighbors=min(int(a.get("n_neighbors", 15)), n - 1))
    say("PCA and neighbours")
    sc.tl.umap(sub, random_state=0)
    say("UMAP")
    umap = np.full((len(mask), 2), np.nan, dtype=np.float32)
    umap[mask] = sub.obsm["X_umap"]
    ws.put_obsm(f"X_umap_{name}", umap, {"step": "recluster", **a})
    resolutions = [float(r) for r in a.get("resolutions", [0.5])]
    primary = float(a.get("primary_resolution", resolutions[0]))
    label_col = a.get("label_column")
    labels = column(table, ws, label_col)[mask].to_numpy() if label_col else None
    clusterings = []
    for res in resolutions:
        key = f"leiden_{res:g}"
        sc.tl.leiden(sub, resolution=res, key_added=key, flavor="igraph", n_iterations=2, random_state=0)
        codes = sub.obs[key].cat.codes.to_numpy()
        names = [f"C{c}" for c in sub.obs[key].cat.categories]
        full = np.full(len(mask), "", dtype=object)
        full[mask] = np.array(names, dtype=object)[codes]
        col = f"{name}_leiden_{res:g}"
        ws.put_obs(col, full, {"step": "recluster", "resolution": res, **a})
        say(f"Leiden {res:g}: {len(names)} clusters")
        entry = {"column": col, "clusters": len(names)}
        if res == primary:
            mean, pct, sizes = group_stats(lognorm, counts, codes, len(names))
            top = top_genes(mean, pct, sizes, genes, names, int(a.get("top_n", 8)))
            clusters = []
            for j, cname in enumerate(names):
                c = {"cluster": cname, "cells": int(sizes[j]), "top_genes": [g["gene"] for g in top[cname]]}
                if labels is not None:
                    vc = pd.Series(labels[codes == j]).value_counts()
                    c[f"main_{label_col}"] = [{"value": v, "share": round(float(k) / sizes[j], 3)} for v, k in vc.head(3).items()]
                if a.get("within"):
                    sec = column(table, ws, a["within"])[mask].to_numpy()[codes == j]
                    vc = pd.Series(sec).value_counts()
                    c["by_section"] = {str(v): int(k) for v, k in vc.items()}
                clusters.append(c)
            entry["detail"] = clusters
        clusterings.append({"resolution": res, **entry})
    return {"cells": n, "embedding": f"X_umap_{name}", "clusterings": clusterings, "primary_resolution": primary}


def markers(table: Table, ws: Workspace, a: dict) -> dict:
    col_name = a["column"]
    labels = column(table, ws, col_name)
    mask = pick(table, ws, a.get("cells")) & (labels != "").to_numpy()
    genes = table.genes
    top_n = int(a.get("top_n", 10))
    groups = sorted(pd.unique(labels[mask]))
    used = "computed"
    if col_name in table.stored_de() and not a.get("cells"):
        # the store's own DE stats (de_stats.py): no X read
        import anndata as ad

        rec = ad.io.read_elem(table.group["uns"]["merfisheyes_de"][col_name])
        stored_groups = [str(g) for g in rec["categories"]]
        mean = np.asarray(table.group["varm"][f"de_mean__{col_name}"][:])
        pct = np.asarray(table.group["varm"][f"de_pct__{col_name}"][:])
        sizes = np.asarray(rec["cell_counts"], dtype=np.float64)
        keep = [i for i, g in enumerate(stored_groups) if g in set(groups)]
        mean, pct, sizes = mean[:, keep], pct[:, keep], sizes[keep]
        groups = [stored_groups[i] for i in keep]
        used = "stored"
    else:
        counts = table.counts()[mask]
        lognorm = log_normalized(counts)
        codes = pd.Categorical(labels[mask], categories=groups).codes
        mean, pct, sizes = group_stats(lognorm, counts, codes, len(groups))
    say(f"stats for {len(groups)} values of {col_name} ({used})")
    if a.get("values"):
        groups_shown = [g for g in a["values"] if g in groups]
    else:
        groups_shown = groups
    top = top_genes(mean, pct, sizes, genes, groups, top_n)
    out = {"column": col_name, "stats": used, "method": "mean of log1p(normalize_total(X)) per value; ranked by mean minus the mean over the other cells",
           "values": {g: {"cells": int(sizes[groups.index(g)]), "top": top[g]} for g in groups_shown[:60]}}
    if a.get("genes"):
        lower = {g.lower(): i for i, g in enumerate(genes)}
        rows = [(g, lower.get(g.lower())) for g in a["genes"]]
        have = [(genes[i], i) for g, i in rows if i is not None]
        out["genes_not_in_panel"] = [g for g, i in rows if i is None]
        if have and a.get("heatmap"):
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            cols = [groups.index(g) for g in groups_shown[:60]]
            m = mean[[i for _, i in have]][:, cols]
            if a.get("scale", "zscore") == "zscore":
                sd = m.std(axis=1, keepdims=True)
                sd[sd == 0] = 1
                m = (m - m.mean(axis=1, keepdims=True)) / sd
                label, cmap, kw = "row z-score of mean expression", "RdBu_r", {"vmin": -2.5, "vmax": 2.5}
            else:
                lo, hi = m.min(axis=1, keepdims=True), m.max(axis=1, keepdims=True)
                m = (m - lo) / np.where(hi > lo, hi - lo, 1)
                label, cmap, kw = "min–max scaled mean", "viridis", {}
            fig, ax = plt.subplots(figsize=(max(4, 0.45 * len(cols) + 2), max(3, 0.3 * len(have) + 1.5)))
            im = ax.imshow(m, aspect="auto", cmap=cmap, **kw)
            ax.set_yticks(range(len(have)), [g for g, _ in have])
            ax.set_xticks(range(len(cols)), [groups_shown[k][:24] for k in range(len(cols))], rotation=60, ha="right")
            ax.set_title(a.get("title", f"{col_name}: chosen genes"))
            fig.colorbar(im, ax=ax, label=label)
            fig.tight_layout()
            fig.savefig(a["heatmap"], dpi=130)
            plt.close(fig)
            out["heatmap"] = a["heatmap"]
    return out


def save(table: Table, ws: Workspace, a: dict) -> dict:
    import anndata as ad
    import zarr

    fmt = a["format"]
    dest = Path(a["dest"])
    names_obs = a.get("columns") or list(ws.index["obs"])
    names_obsm = a.get("embeddings") or list(ws.index["obsm"])
    if dest.exists() and not a.get("overwrite"):
        raise SystemExit(f"{dest} exists; choose another name or pass overwrite.")
    if fmt == "zarr":
        from spatialdata_cluster import write_column

        if dest.exists():
            shutil.rmtree(dest)
        say(f"copying the store to {dest}")
        shutil.copytree(table.store, dest)
        for n in names_obs:
            write_column(dest, table.name, n, pd.Categorical(ws.obs(n).astype(str)))
        group = zarr.open_group(str(dest), mode="r+", use_consolidated=False)[f"tables/{table.name}"]
        obsm = group.require_group("obsm")
        for n in names_obsm:
            if n in obsm:
                del obsm[n]
            ad.io.write_elem(obsm, n, ws.obsm(n))
        de_cols = [n for n in names_obs if len(set(ws.obs(n))) <= 500]
        if de_cols:
            import de_stats

            de_stats.write_de_stats(group, say, columns=de_cols)
        zarr.consolidate_metadata(str(dest))
        say("written")
    elif fmt == "h5ad":
        source = a.get("original")
        if source and str(source).lower().endswith(".h5ad"):
            adata = ad.read_h5ad(source)
            order = table.obs.index.astype(str)
            if list(adata.obs_names.astype(str)) == list(order):
                take = np.arange(len(order))
            else:
                at = pd.Series(np.arange(len(order)), index=order)
                take = at.reindex(adata.obs_names.astype(str)).to_numpy()
                if np.isnan(take).all():
                    raise SystemExit("The original .h5ad's cells do not match the dataset shown (no shared cell ids).")
        else:
            adata = ad.io.read_elem(table.group)
            take = np.arange(adata.n_obs)
        ok = ~pd.isna(take)
        idx = np.where(ok, take, 0).astype(int)
        for n in names_obs:
            v = ws.obs(n).astype(str).to_numpy()[idx]
            adata.obs[n] = pd.Categorical(np.where(ok, v, ""))
        for n in names_obsm:
            m = ws.obsm(n)[idx]
            m[~ok] = np.nan
            adata.obsm[n] = m
        dest.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(dest)
        say("written")
    else:
        raise SystemExit(f"Unknown format {fmt!r}: zarr or h5ad.")
    return {"path": str(dest), "format": fmt, "columns": names_obs, "embeddings": names_obsm}


STEPS = {"profile": profile, "niche": niche, "recluster": recluster, "markers": markers, "save": save}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("step", choices=sorted(STEPS))
    p.add_argument("--store", required=True)
    p.add_argument("--work", required=True)
    p.add_argument("--args", required=True, help="a JSON file of the step's arguments")
    p.add_argument("--out", required=True, help="where to write the result (JSON)")
    o = p.parse_args()
    args = json.loads(Path(o.args).read_text())
    table = Table(Path(o.store))
    ws = Workspace(Path(o.work), len(table.obs))
    result = STEPS[o.step](table, ws, args)
    Path(o.out).write_text(json.dumps({"ok": True, "seconds": round(time.time() - T0, 1), **result}, default=str))
    say("done")


if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        if e.code not in (0, None) and not isinstance(e.code, int):
            print(f"ERROR: {e.code}", flush=True)
            sys.exit(2)
        raise
