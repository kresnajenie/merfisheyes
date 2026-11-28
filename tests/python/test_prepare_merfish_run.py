"""scripts/prepare_merfish_run.py (through scripts/prepare.py): a raw MERFISH
run becomes one store — cells after the overlap rule, molecules without
blanks, sets as sections side by side — on a small run made here."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import prepare  # noqa: E402
import prepare_merfish_run as pmr  # noqa: E402
import zarr  # noqa: E402
from test_raw_merfish import PANEL, make_run  # noqa: E402


def quiet(*_):
    pass


def table(store: Path):
    import anndata as ad

    return ad.read_zarr(store / "tables" / "table")


def molecules(store: Path) -> int:
    return int(zarr.open_group(store / "points" / "molecules", mode="r").attrs["points:n_points"])


def labels_file(path: Path, names: list[str]) -> Path:
    import anndata as ad
    import pandas as pd

    obs = pd.DataFrame(index=pd.Index(names))
    obs["subclass_name"] = pd.Categorical(["Astro" if i % 2 else "Oligo" for i in range(len(names))])
    obs["score"] = np.arange(len(names), dtype=float)
    a = ad.AnnData(obs=obs, obsm={"X_umap": np.c_[np.arange(len(names)), -np.arange(len(names))].astype(float),
                                  "X_spatial": np.zeros((len(names), 2)),
                                  "X_shifted": np.zeros((len(names), 2))})
    a.write_h5ad(path)
    return path


def test_a_run_becomes_one_store(tmp_path):
    truth = make_run(tmp_path / "run", grid=3)
    assert prepare.prepare(tmp_path / "run", tmp_path / "out", name="run", log=quiet) == 0
    store = tmp_path / "out" / "run.zarr"
    t = table(store)
    assert t.n_obs == truth["cells"]
    assert list(t.var_names) == PANEL[:5]  # blanks out of X, counted per cell
    assert t.obs["n_blank"].sum() > 0 and "gfp_count" not in t.obs
    assert list(t.obs["section"].cat.categories) == ["1", "2"]
    # side by side: set 2's cells start right of set 1's, each in its own frame too
    s1 = t.obsm["spatial"][t.obs["section"] == "1"]
    s2 = t.obsm["spatial"][t.obs["section"] == "2"]
    assert s2[:, 0].min() > s1[:, 0].max() + pmr.PADDING_UM
    assert "spatial_section" not in t.obsm  # one spatial frame only
    non_blank = sum(v for sid in truth["per_gene"] for g, v in truth["per_gene"][sid].items() if not g.startswith("blank"))
    assert molecules(store) == non_blank
    sec = json.loads(t.uns["merfisheyes_sections"])
    assert [s["section"] for s in sec] == ["1", "2"]
    stats = json.loads((tmp_path / "out" / "stats.json").read_text())
    assert stats["input"] == "merfish_run" and stats["molecules"] == non_blank and stats["cells"] == truth["cells"]
    assert "section" in stats["de_stats"]


def test_sections_read_top_down_each_flipped(tmp_path):
    """y is flipped (the viewer draws +y upwards) so a section reads as in
    the lab's plots; rows of sections go down, the first set top-left."""
    import raw_merfish as rm

    make_run(tmp_path / "run", sets=("1", "2", "3"), grid=3)
    assert prepare.prepare(tmp_path / "run", tmp_path / "out", name="run", log=quiet) == 0
    t = table(tmp_path / "out" / "run.zarr")
    shifts = {s["section"]: s for s in json.loads(t.uns["merfisheyes_sections"])}
    ys = {k: t.obsm["spatial"][(t.obs["section"] == k).to_numpy(), 1] for k in "123"}
    assert ys["3"].max() < ys["1"].min() - pmr.PADDING_UM + 5  # set 3 starts the second row, below set 1
    # every cell: x = x_set + dx, y = dy - y_set, from its own mat
    run = rm.detect(tmp_path / "run")
    _, fields, _ = rm.read_cells(run, progress=quiet)
    raw = {n: xum for f in fields for n, xum, k in zip(f.cells, f.xum, f.keep) if k}
    for name, (x, y), sec in zip(t.obs_names[:200], t.obsm["spatial"][:200], t.obs["section"][:200]):
        sh = shifts[str(sec)]
        assert np.isclose(x, raw[name][1] + sh["dx"]) and np.isclose(y, sh["dy"] - raw[name][2])
    # the molecules follow the same rule: their extent sits inside the cells' slots
    import points_zarr as pz

    m = pz.read_dataframe(tmp_path / "out" / "run.zarr" / "points" / "molecules", level=0, columns=["x", "y"])
    assert m["y"].max() <= max(v.max() for v in ys.values()) + 5
    assert m["y"].min() >= min(v.min() for v in ys.values()) - 5


def test_some_sets_and_labels_by_name(tmp_path):
    make_run(tmp_path / "run", grid=3)
    one = make_run(tmp_path / "ref", sets=("2",), grid=3)  # the same set 2, to count it
    names = [f"Conv_zscan2_0{i}--{k}" for i in range(3) for k in range(5)]
    h5ad = labels_file(tmp_path / "labels.h5ad", names)
    code = prepare.prepare(tmp_path / "run", tmp_path / "out", name="two", sets=["2"],
                           labels_from=h5ad, log=quiet)
    assert code == 0
    t = table(tmp_path / "out" / "two.zarr")
    assert t.n_obs == one["cells"] and list(t.obs["section"].cat.categories) == ["2"]
    have = t.obs["subclass_name"].astype(str)
    assert set(have[have != ""]) == {"Astro", "Oligo"} and (have == "").sum() > 0
    assert "X_umap" in t.obsm and "X_spatial" not in t.obsm
    assert "X_shifted" not in t.obsm  # a layout of its own, not an embedding
    assert np.isnan(t.obsm["X_umap"][have.to_numpy() == ""]).all()


def test_an_interrupted_conversion_resumes(tmp_path, monkeypatch):
    truth = make_run(tmp_path / "run", grid=3)
    real = pmr._spill_set

    def fail_on_set_2(job):
        if job["set"].id == "2":
            raise RuntimeError("the computer went to sleep")
        return real(job)

    monkeypatch.setattr(pmr, "_spill_set", fail_on_set_2)
    with pytest.raises(RuntimeError):
        pmr.prepare_merfish_run(tmp_path / "run", tmp_path / "out", name="run", workers=1, log=quiet)
    progress = json.loads((tmp_path / "out" / ".molecules.tmp" / pmr.PROGRESS).read_text())
    assert list(progress) == ["1"]

    calls = []
    monkeypatch.setattr(pmr, "_spill_set", lambda job: calls.append(job["set"].id) or real(job))
    lines = []
    r = pmr.prepare_merfish_run(tmp_path / "run", tmp_path / "out", name="run", workers=1,
                                resume=True, log=lines.append)
    assert not r.failed and calls == ["2"]  # set 1 was not read again
    assert any("resuming" in str(m) for m in lines)
    store = tmp_path / "out" / "run.zarr"
    non_blank = sum(v for sid in truth["per_gene"] for g, v in truth["per_gene"][sid].items() if not g.startswith("blank"))
    assert molecules(store) == non_blank
    assert not (tmp_path / "out" / ".molecules.tmp").exists()


def test_bands_shrink_for_dense_wide_elements():
    import points_zarr as pz

    assert pz.band_for(1_000_000, 2_000) == 800.0  # sparse: the default
    assert pz.band_for(112_843_931, 2_331) == 500.0  # two sets
    assert pz.band_for(1_709_331_981, 11_243) == 100.0  # 18 side by side: one tile row
    assert pz.band_for(10**11, 1_000) == 100.0  # never below a tile


def test_prepare_recognises_a_run_beside_h5ad_files(tmp_path):
    make_run(tmp_path / "run", sets=("1",))
    labels_file(tmp_path / "run" / "combined.h5ad", ["Conv_zscan1_00--0"])
    assert prepare.recognise(tmp_path / "run") == "merfish_run"
    (tmp_path / "bad" / "final_spots").mkdir(parents=True)
    with pytest.raises(prepare.Unrecognised, match="holds no .npz spot files"):
        prepare.recognise(tmp_path / "bad")
