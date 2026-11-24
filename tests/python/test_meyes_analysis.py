"""scripts/meyes_analysis.py: the analysis steps an agent runs on an open
dataset (MCP.md pillar D, tier 2), on the committed tiny store."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import zarr

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
FIXTURE = ROOT / "tests" / "data" / "tiny" / "zarr" / "spatialdata_v3_de.zarr"
sys.path.insert(0, str(SCRIPTS))

# the analysis runs in the pipeline's Python (scripts/requirements.txt)
pytest.importorskip("scanpy")
pytest.importorskip("leidenalg")

import meyes_analysis as m  # noqa: E402


@pytest.fixture()
def env(tmp_path):
    store = tmp_path / "s.zarr"
    shutil.copytree(FIXTURE, store)
    table = m.Table(store)
    return store, table, m.Workspace(tmp_path / "work", len(table.obs))


def test_profile(env):
    store, table, ws = env
    out = m.profile(table, ws, {"genes": [table.genes[0], "NotAGene"], "breakdown": {"column": "celltype", "by": "region"}})
    assert out["cells"] == 30 and out["genes"] == 8
    kinds = {c["name"]: c["kind"] for c in out["columns"]}
    assert kinds["celltype"] == "categorical" and kinds["area"] == "numeric"
    assert {v["value"]: v["cells"] for v in next(c for c in out["columns"] if c["name"] == "celltype")["values"]} == {
        "astro": 10, "microglia": 10, "neuron": 10}
    assert out["markers"]["not_in_panel"] == ["NotAGene"]
    assert out["markers"]["measured"][0]["gene"] == table.genes[0]
    assert out["stored_de_columns"] == ["celltype"]
    assert out["breakdown"]["counts"]["neuron"]["cell_labels"] == 10
    assert out["coordinates"]["median_nearest_cell"] > 0


def test_niche_is_within_the_radius(env):
    store, table, ws = env
    out = m.niche(table, ws, {"name": "N", "anchors": {"column": "celltype", "values": ["neuron"]}, "radius_um": 10})
    inside = ws.obs("niche_N").to_numpy() == "N"
    assert out["cells"] == int(inside.sum()) and out["anchors"] == 10
    # by hand: every cell within 10 of a neuron, the neurons included
    xy = table.coords()
    anchors = xy[table.obs["celltype"].to_numpy() == "neuron"]
    d = np.sqrt(((xy[:, None, :] - anchors[None]) ** 2).sum(-1)).min(1)
    np.testing.assert_array_equal(inside, d <= 10)
    # an excluded label is in neither the anchors nor the niche
    m.niche(table, ws, {"name": "M", "anchors": {"column": "celltype", "values": ["neuron", "astro"]},
                        "radius_um": 10, "exclude": {"column": "celltype", "values": ["astro"]}})
    inside_m = ws.obs("niche_M").to_numpy() == "M"
    assert not inside_m[table.obs["celltype"].to_numpy() == "astro"].any()


def test_unknown_values_name_the_close_ones(env):
    store, table, ws = env
    with pytest.raises(SystemExit, match="Values containing those words: 'neuron'"):
        m.niche(table, ws, {"name": "X", "anchors": {"column": "celltype", "values": ["neuron cells"]}, "radius_um": 5})


def test_recluster_markers_and_save(env, tmp_path):
    store, table, ws = env
    out = m.recluster(table, ws, {"name": "all", "cells": {"column": "celltype", "values": ["astro", "microglia", "neuron"]},
                                  "resolutions": [0.5, 1.0], "label_column": "celltype", "n_pcs": 5, "n_neighbors": 5})
    assert [c["resolution"] for c in out["clusterings"]] == [0.5, 1.0]
    primary = out["clusterings"][0]
    assert sum(c["cells"] for c in primary["detail"]) == 30
    assert ws.obsm("X_umap_all").shape == (30, 2)
    labels = ws.obs(primary["column"])
    assert labels.str.startswith("C").all()

    stored = m.markers(table, ws, {"column": "celltype", "top_n": 3})
    assert stored["stats"] == "stored" and set(stored["values"]) == {"astro", "microglia", "neuron"}
    heat = tmp_path / "h.png"
    computed = m.markers(table, ws, {"column": primary["column"], "genes": table.genes[:4] + ["Nope"], "heatmap": str(heat)})
    assert computed["stats"] == "computed" and computed["genes_not_in_panel"] == ["Nope"] and heat.exists()

    dest = tmp_path / "copy.zarr"
    saved = m.save(table, ws, {"format": "zarr", "dest": str(dest)})
    assert saved["embeddings"] == ["X_umap_all"]
    copy = ad.read_zarr(dest / "tables" / "table")
    assert primary["column"] in copy.obs and "X_umap_all" in copy.obsm
    # the new labels carry DE stats, and the store's own are kept
    assert {"celltype", primary["column"]} <= set(copy.uns["merfisheyes_de"])
    check = subprocess.run([sys.executable, str(SCRIPTS / "validate_store.py"), str(dest)], capture_output=True, text=True)
    assert check.returncode == 0, check.stdout + check.stderr
    # the original store is untouched
    assert primary["column"] not in ad.read_zarr(store / "tables" / "table").obs

    h5 = tmp_path / "copy.h5ad"
    m.save(table, ws, {"format": "h5ad", "dest": str(h5)})
    back = ad.read_h5ad(h5)
    assert back.n_obs == 30 and primary["column"] in back.obs and back.obsm["X_umap_all"].shape == (30, 2)


def test_the_command_line(env, tmp_path):
    store, _, _ = env
    args = tmp_path / "a.json"
    out = tmp_path / "o.json"
    args.write_text(json.dumps({"genes": ["NotAGene"]}))
    run = subprocess.run([sys.executable, str(SCRIPTS / "meyes_analysis.py"), "profile", "--store", str(store),
                          "--work", str(tmp_path / "w"), "--args", str(args), "--out", str(out)], capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    assert json.loads(out.read_text())["ok"] is True
    args.write_text(json.dumps({"name": "X", "anchors": {"column": "nope", "values": ["a"]}, "radius_um": 5}))
    run = subprocess.run([sys.executable, str(SCRIPTS / "meyes_analysis.py"), "niche", "--store", str(store),
                          "--work", str(tmp_path / "w"), "--args", str(args), "--out", str(out)], capture_output=True, text=True)
    assert run.returncode != 0 and 'ERROR: No column "nope"' in run.stdout
