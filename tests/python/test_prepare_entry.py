"""scripts/prepare.py, the one entry point a person runs and the server runs:
what it recognises, exactly what it writes, that it leaves INPUT alone, that
it writes the same bytes twice, and the rules it applies on the way (control
probes out, uns entries the writer refuses left out, coordinates from obs,
per-cell labels in)."""
from __future__ import annotations

import filecmp
import json
import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import zarr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import meyes_common  # noqa: E402
import prepare  # noqa: E402

TINY = ROOT / "tests" / "data" / "tiny"


def cells(path: Path, *, obsm=True, genes=("Gad1", "Sst", "Pvalb", "NegControlProbe_00042"), n=30, uns=None) -> Path:
    rng = np.random.default_rng(0)
    obs = pd.DataFrame(
        {"celltype": pd.Categorical(["a", "b", "c"] * (n // 3)), "center_x": rng.uniform(0, 100, n), "center_y": rng.uniform(0, 100, n)},
        index=pd.Index([f"cell{i}" for i in range(n)]),
    )
    adata = ad.AnnData(
        X=rng.poisson(1.0, (n, len(genes))).astype(np.float32),
        obs=obs,
        var=pd.DataFrame(index=list(genes)),
        obsm={"spatial": obs[["center_x", "center_y"]].to_numpy()} if obsm else {},
        uns=uns or {},
    )
    adata.write_h5ad(path)
    return path


def molecules(path: Path, n=2000) -> Path:
    """A lab's own molecule table: Gene, X, Y, Z (and columns that are not)."""
    rng = np.random.default_rng(1)
    genes = np.array(["Gad1", "Sst", "Pvalb", "Blank-3"])[rng.integers(0, 4, n)]
    pq.write_table(
        pa.table({"iGene": rng.integers(0, 4, n), "Gene": genes, "X": rng.uniform(0, 100, n), "Y": rng.uniform(0, 100, n), "Z": rng.uniform(0, 5, n), "Cell ID": rng.integers(0, 30, n)}),
        path,
    )
    return path


def run(input_path: Path, out: Path, *args) -> int:
    return prepare.main([str(input_path), str(out), *args])


def tree(folder: Path) -> list[str]:
    return sorted(str(p.relative_to(folder)) for p in folder.rglob("*"))


def test_recognises_each_kind_by_content(tmp_path):
    h5ad = cells(tmp_path / "anything.h5ad")
    lab = molecules(tmp_path / "set1__single_molecules.parquet")
    parts = tmp_path / "parts"
    parts.mkdir()
    shutil.copy(h5ad, parts / "c.h5ad")
    shutil.copy(lab, parts / "m.parquet")

    assert prepare.recognise(h5ad) == "h5ad"
    assert prepare.recognise(lab) == "transcripts"
    assert prepare.recognise(TINY / "single-molecule" / "tiny.parquet") == "transcripts"
    assert prepare.recognise(parts) == "folder"
    with pytest.raises(prepare.Unrecognised):
        prepare.recognise(ROOT / "README.md")


def test_unrecognised_input_exits_2(tmp_path):
    assert run(ROOT / "README.md", tmp_path / "out") == 2


def test_a_folder_of_parts_becomes_exactly_the_contract(tmp_path):
    parts = tmp_path / "sample"
    parts.mkdir()
    cells(parts / "cells.h5ad")
    molecules(parts / "set1__single_molecules.parquet")
    before = tree(parts)

    assert run(parts, tmp_path / "out") == 0
    out = tmp_path / "out"
    assert sorted(p.name for p in out.iterdir()) == ["mapping.json", "sample.zarr", "stats.json"]
    assert tree(parts) == before, "INPUT is never changed"
    assert json.loads((out / "mapping.json").read_text()) == {"cells": "sample.zarr"}

    stats = json.loads((out / "stats.json").read_text())
    assert stats["input"] == "folder"
    assert (stats["cells"], stats["genes"]) == (30, 3), "the control probe is left out"
    assert stats["molecule_genes"] == 3, "Blank-3 is left out"
    assert 0 < stats["molecules"] < 2000
    assert stats["pipeline"] and "anndata" in stats["packages"]

    # Readable by the libraries everyone else uses
    table = ad.read_zarr(out / "sample.zarr" / "tables" / "table")
    assert list(table.var_names) == ["Gad1", "Sst", "Pvalb"]
    assert "X_csc" in table.layers


def test_same_input_same_bytes(tmp_path):
    parts = tmp_path / "sample"
    parts.mkdir()
    cells(parts / "cells.h5ad")
    molecules(parts / "m.parquet")
    assert run(parts, tmp_path / "a") == 0
    assert run(parts, tmp_path / "b") == 0

    a, b = tmp_path / "a" / "sample.zarr", tmp_path / "b" / "sample.zarr"
    assert tree(a) == tree(b)
    differ = [p for p in tree(a) if (a / p).is_file() and not filecmp.cmp(a / p, b / p, shallow=False)]
    assert differ == []


def test_molecules_only(tmp_path):
    assert run(molecules(tmp_path / "lab.parquet"), tmp_path / "out") == 0
    stats = json.loads((tmp_path / "out" / "stats.json").read_text())
    assert (stats["input"], stats["cells"], stats["molecule_genes"]) == ("transcripts", 0, 3)


def test_an_uns_entry_the_writer_refuses_is_left_out(tmp_path, capsys):
    h5ad = cells(tmp_path / "c.h5ad")
    adata = ad.read_h5ad(h5ad)
    out = zarr.open_group(tmp_path / "t.zarr", mode="w", zarr_format=3)

    adata.uns = {"kept": "fine", "refused": object()}
    import prepare_folder

    prepare_folder.write_table(out, adata, print)
    written = ad.read_zarr(tmp_path / "t.zarr" / "table")
    assert written.uns == {"kept": "fine"}
    assert "uns['refused'] left out" in capsys.readouterr().out


def test_coordinates_from_obs_when_obsm_has_none(tmp_path):
    h5ad = cells(tmp_path / "c.h5ad", obsm=False)
    assert run(h5ad, tmp_path / "out") == 0
    table = ad.read_zarr(tmp_path / "out" / "c.zarr" / "tables" / "table")
    np.testing.assert_allclose(table.obsm["spatial"], table.obs[["center_x", "center_y"]].to_numpy())


def test_per_cell_labels_are_merged(tmp_path):
    h5ad = cells(tmp_path / "c.h5ad")
    labels = tmp_path / "mmc.csv"
    labels.write_text(
        "# MapMyCells output\ncell_id,class_name,score\n"
        + "".join(f"cell{i},{'neuron' if i % 2 else 'glia'},{i / 30:.3f}\n" for i in reversed(range(30)))
    )
    assert run(h5ad, tmp_path / "out", "--labels", str(labels)) == 0
    table = ad.read_zarr(tmp_path / "out" / "c.zarr" / "tables" / "table")
    assert table.obs["class_name"].dtype == "category"
    assert table.obs.loc["cell3", "class_name"] == "neuron", "matched by cell_id, not row order"
    assert table.obs.loc["cell3", "score"] == pytest.approx(0.1)


def test_the_column_names_agree_everywhere():
    # worker/entrypoint.py and the browser read the same names as the scripts
    sys.path.insert(0, str(ROOT / "worker"))
    import types

    sys.modules.setdefault("boto3", types.ModuleType("boto3"))
    import entrypoint

    assert entrypoint.MOLECULE_COLUMN_NAMES == meyes_common.MOLECULE_COLUMN_NAMES
