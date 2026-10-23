"""Per-category expression summaries, kept inside the table itself.

The viewer's DEG panel ranks genes for a category from two numbers per gene
and category: the mean of the log-normalised expression (normalize_total,
then log1p, as scanpy does) and the fraction of cells expressing the gene
(count > 0). They are written for every categorical obs column into the
AnnData table, so the store carries them and no folder beside it is needed:

  tables/<table>/varm/de_mean__<column>   float32, genes × categories
  tables/<table>/varm/de_pct__<column>    float32, genes × categories
  tables/<table>/uns/merfisheyes_de/<column>/
      categories    the category of each varm column, in order
      cell_counts   cells in each category
      method        how the numbers were computed

varm and uns are ordinary AnnData elements: anndata and spatialdata read and
write them with the table. lib/spatialdata/SpatialDataTableAdapter.ts reads
them (and falls back to an annotations/ folder for stores written before).
"""
from __future__ import annotations

from pathlib import Path

DE_KEY = "merfisheyes_de"
METHOD = "mean of log1p(normalize_total(X)) (scanpy defaults) and fraction of cells with X > 0, per category"
# A column with more categories than this is an identifier, not a grouping
MAX_CATEGORIES = 500


def categorical_columns(table) -> list[str]:
    """The obs columns to summarise: categorical, with 2..MAX_CATEGORIES
    categories, not the index nor the SpatialData join keys."""
    obs = table["obs"]
    skip = {obs.attrs.get("_index", "_index"), table.attrs.get("region_key"), table.attrs.get("instance_key")}
    found = []
    for name in obs.attrs.get("column-order", []):
        if name in skip or name not in obs:
            continue
        node = obs[name]
        if not hasattr(node, "attrs") or node.attrs.get("encoding-type") != "categorical":
            continue
        n = node["categories"].shape[0]
        if 2 <= n <= MAX_CATEGORIES:
            found.append(name)
    return found


def log_normalize(counts):
    """log1p(normalize_total(counts)) as scanpy computes them (sc.pp.normalize_total
    with its defaults, then sc.pp.log1p): every cell scaled to the median of
    the cells' totals (cells with none left at zero), then natural log1p.
    Here without scanpy, which the server does not need for this."""
    import numpy as np
    import scipy.sparse as sp

    counts = sp.csr_matrix(counts, dtype=np.float32)
    totals = np.ravel(counts.sum(axis=1))
    nonzero = totals > 0
    target = float(np.median(totals[nonzero])) if nonzero.any() else 1.0
    scale = np.zeros_like(totals)
    scale[nonzero] = target / totals[nonzero]
    lognorm = sp.diags(scale.astype(np.float32)) @ counts
    lognorm.data = np.log1p(lognorm.data)
    return sp.csr_matrix(lognorm)


def write_de_stats(table, log=print, columns: list[str] | None = None) -> list[str]:
    """Write the summaries for `columns` (default: every categorical column)
    into the open table group `table`. Returns the columns written."""
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    columns = categorical_columns(table) if columns is None else columns
    if not columns or "X" not in table:
        return []

    counts = sp.csr_matrix(ad.io.read_elem(table["X"])).astype(np.float32)
    expressing = (counts > 0).astype(np.float32)
    lognorm = log_normalize(counts)
    obs = ad.io.read_elem(table["obs"])

    varm = table.require_group("varm")
    uns = table.require_group("uns")
    record = {}
    for column in columns:
        values = obs[column]
        categories = [str(c) for c in values.cat.categories]
        codes = values.cat.codes.to_numpy()
        assigned = codes >= 0
        membership = sp.csr_matrix(
            (np.ones(int(assigned.sum()), dtype=np.float32), (np.flatnonzero(assigned), codes[assigned])),
            shape=(len(codes), len(categories)),
        )
        cell_counts = np.bincount(codes[assigned], minlength=len(categories))
        with np.errstate(divide="ignore", invalid="ignore"):
            means = np.nan_to_num(np.asarray((lognorm.T @ membership).todense()) / cell_counts)
            pct = np.nan_to_num(np.asarray((expressing.T @ membership).todense()) / cell_counts)
        for key, array in ((f"de_mean__{column}", means), (f"de_pct__{column}", pct)):
            if key in varm:
                del varm[key]
            ad.io.write_elem(varm, key, np.ascontiguousarray(array, dtype=np.float32))
        record[column] = {
            "categories": np.array(categories, dtype=object),
            "cell_counts": cell_counts.astype(np.int64),
            "method": METHOD,
        }
        log(f"  DE stats for obs/{column}: {len(categories)} categories × {means.shape[0]} genes")
    if DE_KEY in uns:
        del uns[DE_KEY]
    ad.io.write_elem(uns, DE_KEY, record)
    return columns


def write_store_de_stats(store: Path, table_name: str = "table", log=print) -> list[str]:
    """write_de_stats on `tables/<table_name>` of the store at `store`."""
    import prepare_folder

    return write_de_stats(prepare_folder.table_group(Path(store), table_name), log)
