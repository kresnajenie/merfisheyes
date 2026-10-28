"""Cell ids repeat across slices, and the pipeline must survive it.

A MERSCOPE instrument numbers cells per slice, so two slices of one dataset
hand out the same id to different cells: 54,553 of ace-dip-vat's 474,986 are
duplicated. What identifies a cell in a combined dataset is
`<_sample_id>_<id>` — which is what `combine_slices_v3.py` tags, what
`map_my_cell.py` writes in `mapping_output.csv`, and what
`process_spatial_data.py` indexes by.

Before this was handled, `counts.loc[ids]` on the duplicated index multiplied
those rows together and `prepare.py` died on 474,986 cells against 584,092
rows of expression — after the combine and MapMyCells had already run.

  python -m pytest tests/bil/test_duplicate_cell_ids.py -q
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from conftest import GENES, ROOT, SCRIPTS

N = 60


@pytest.fixture(scope="module")
def combined(tmp_path_factory) -> Path:
    """A combined_output whose two slices share every cell id."""
    out = tmp_path_factory.mktemp("dupes") / "combined_output"
    out.mkdir(parents=True)
    rng = np.random.default_rng(3)
    ids = [f"{rng.integers(10**17, 10**18)}{i:04d}" for i in range(N)]

    meta, counts = [], []
    for s, (sample, dx) in enumerate([("1231122828", 0.0), ("1231122831", 5000.0)]):
        cx = rng.uniform(0, 100, N).round(3) + dx
        cy = rng.uniform(0, 80, N).round(3)
        meta.append(pd.DataFrame({
            "id": ids, "fov": 0, "volume": rng.uniform(100, 900, N).round(3),
            "center_x": cx, "center_y": cy,
            "min_x": cx - 3, "min_y": cy - 3, "max_x": cx + 3, "max_y": cy + 3,
            "_source_file": f"/in/{sample}/cellpose_metadata.csv", "_sample_id": sample,
        }))
        c = pd.DataFrame(rng.integers(0, 12, (N, len(GENES))), columns=GENES)
        c.insert(0, "cell", ids)
        counts.append(c)

    pd.concat(meta).to_csv(out / "cell_metadata.csv", index=False)
    pd.concat(counts).to_csv(out / "cell_by_gene.csv", index=False)
    # the mask the combine writes: keyed by the id the slice gave the cell
    pd.DataFrame({"cell_id": ids * 2,
                  "is_artifact": ["True"] * 5 + ["False"] * (2 * N - 5)}).to_csv(
        out / "artifact_mask_p25.csv", index=False)
    # what MapMyCells writes: keyed by <sample>_<id>
    pd.DataFrame({
        "cell_id": [f"1231122828_{i}" for i in ids] + [f"1231122831_{i}" for i in ids],
        "class_name": ["Glut"] * N + ["GABA"] * N,
    }).to_csv(out / "mapping_output.csv", index=False)
    return out


def test_a_store_is_written_and_every_cell_is_kept_once(combined, tmp_path):
    """The whole point: 2 x N cells in, 2 x N cells out."""
    out = tmp_path / "store"
    r = subprocess.run(
        [sys.executable, "-I", str(SCRIPTS / "prepare.py"), str(combined), str(out), "--name", "dupes"],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert r.returncode == 0, r.stdout[-2500:] + r.stderr[-2500:]

    import anndata as ad

    table = ad.read_zarr(out / "dupes.zarr" / "tables" / "table")
    assert table.n_obs == 2 * N, "rows were multiplied together, or dropped"
    assert table.obs_names.is_unique, "the index does not identify a cell"
    assert list(table.obs_names[:1])[0].startswith("1231122828_")
    assert "slice_cell_id" in table.obs.columns


def test_labels_keyed_either_way_land_on_the_right_cells(combined, tmp_path):
    """MapMyCells keys by <sample>_<id> and joins by it. The artifact mask
    keys by the slice's own id, which two slices share, so it cannot join by
    id and must fall back to the order it was written in — never to a join
    that hands one cell's answer to another."""
    out = tmp_path / "store_labels"
    r = subprocess.run(
        [sys.executable, "-I", str(SCRIPTS / "prepare.py"), str(combined), str(out), "--name", "dupes",
         "--labels", str(combined / "mapping_output.csv"),
         "--labels", str(combined / "artifact_mask_p25.csv")],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert r.returncode == 0, r.stdout[-2500:] + r.stderr[-2500:]
    # MapMyCells keys by <sample>_<id>, which identifies a cell, so it joins
    assert "mapping_output.csv: matched by cell_id" in r.stdout, r.stdout[-2000:]
    # The mask keys by the id the slice gave the cell, which does NOT identify
    # one — two slices use it. Joining on it would hand one cell's answer to
    # another's, so the rows are taken in order, which is the order the
    # combine wrote them in.
    assert "artifact_mask_p25.csv: no cell_id the table knows" in r.stdout, r.stdout[-2000:]

    import anndata as ad

    obs = ad.read_zarr(out / "dupes.zarr" / "tables" / "table").obs
    # sample ids are digits, so pandas reads them as integers unless told
    by_sample = obs.groupby(obs["_sample_id"].astype(str), observed=True)["class_name"].agg(
        lambda s: set(s.astype(str)))
    assert by_sample["1231122828"] == {"Glut"}
    assert by_sample["1231122831"] == {"GABA"}
    # taken in order: the five flagged rows are the first five cells
    assert (obs["is_artifact"].astype(str) == "True").sum() == 5
    assert list(obs["is_artifact"].astype(str)[:5]) == ["True"] * 5
