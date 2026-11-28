"""scripts/raw_merfish.py: a raw MERFISH run recognised, read in place and
checked (docs/specs/raw-merfish-profile.md), on a small run made here."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import raw_merfish as rm  # noqa: E402

PANEL = ["Slc17a7", "Gad2", "Mbp", "Sst", "Vip", "blank0001", "blank0002"]
PIX = np.array([0.4, 0.5, 0.5])  # µm per z plane, x pixel, y pixel
FOV = np.array([10, 100, 100])  # a field is 50 µm wide
STEP = 90  # pixels: neighbours overlap by 10 pixels (5 µm)


def fields_of(grid: int) -> np.ndarray:
    return np.array([(i * STEP, j * STEP) for j in range(grid) for i in range(grid)], np.float32)


def make_run(root: Path, sets=("1", "2"), grid=2, dedupe=True, bad_field=None,
             compressed=False, seed=0) -> dict:
    """A run: per set, molecules kept in their nearest field (unless
    dedupe=False: every field keeps its own), and per field a mat whose
    cells near an edge appear again in the neighbour's mat."""
    rng = np.random.default_rng(seed)
    (root / "final_spots").mkdir(parents=True)
    (root / "final_mats").mkdir()
    truth = {"cells": 0, "molecules": {}, "per_gene": {}}
    for sid in sets:
        pos = fields_of(grid)
        centres_px = pos + FOV[1:] / 2
        fovs = [f"Conv_zscan{sid}_{i:02d}" for i in range(len(pos))]
        xs, genes = [], []
        for i, p in enumerate(pos):
            n = 4000
            xy = p + rng.uniform(0, FOV[1], (n, 2))
            if dedupe:
                own = rm.nearest_field(xy * PIX[1], centres_px * PIX[1]) == i
                xy = xy[own]
            g = rng.choice(PANEL[:5], len(xy), p=[0.4, 0.3, 0.2, 0.06, 0.04]).astype("<U13")
            if i == bad_field:  # a field that decoded badly: a third blanks
                g[: len(g) // 3] = "blank0001"
            else:
                g[rng.random(len(g)) < 0.005] = "blank0002"
            z = rng.uniform(0, FOV[0], len(xy))
            xs.append(np.c_[z, xy, rng.uniform(-0.5, 0, len(xy)), np.full(len(xy), i)].astype(np.float32))
            genes.append(g)
        XF, G = np.concatenate(xs), np.concatenate(genes)
        save = np.savez_compressed if compressed else np.savez
        save(root / "final_spots" / f"spots_allth3_final_set{sid}.npz", XF=XF, genes=G,
             pix_size=PIX, fov_size=FOV, pos_fovs_pix=pos, fovs=np.array(fovs))
        truth["molecules"][sid] = len(XF)
        truth["per_gene"][sid] = {g: int((G == g).sum()) for g in PANEL}

        # cells: a grid of true cells over the set; each field's mat holds
        # every cell inside it (so overlap cells are in two or four mats)
        span = STEP * (grid - 1) + FOV[1]  # from 5 px: no cell on a midpoint between centres
        cx, cy = np.meshgrid(np.arange(5, span - 2, 7.0), np.arange(5, span - 2, 7.0))
        cells_px = np.c_[cx.ravel(), cy.ravel()]
        owner = rm.nearest_field(cells_px * PIX[1], centres_px * PIX[1])
        truth["cells"] += len(cells_px)
        for i, p in enumerate(pos):
            inside = np.all((cells_px >= p) & (cells_px < p + FOV[1]), axis=1)
            idx = np.flatnonzero(inside)
            counts = rng.poisson(20, (len(idx), len(PANEL))).astype(np.uint32)
            counts[:, 5:] = rng.poisson(0.1, (len(idx), 2))
            if i == bad_field:
                counts[:, 5] = 40
            jitter = np.where(owner[idx] == i, 0.0, 0.05)[:, None]  # a duplicate's centroid differs a little
            xum = np.c_[rng.uniform(0, 4, len(idx)), (cells_px[idx] + jitter) * PIX[1]]
            np.savez(root / "final_mats" / f"{fovs[i]}--mat.npz",
                     Mcells=np.array([f"{fovs[i]}--{k}" for k in range(len(idx))]),
                     Mgenes=np.array(PANEL), Mcounts=counts, Xum=xum)
    return truth


def test_detect_reads_the_run(tmp_path):
    make_run(tmp_path / "run", sets=("10", "2", "1"))
    run = rm.detect(tmp_path / "run")
    assert [s.id for s in run.sets] == ["1", "2", "10"]  # by number, not by name
    assert len(run.mats) == 12 and run.gfp_fields == 0
    s = run.sets[0]
    assert s.stored and s.fovs[0] == "Conv_zscan1_00"
    assert s.extent_um() == (0.0, 0.0, 95.0, 95.0)


@pytest.mark.parametrize("compressed", [False, True])
def test_members_read_in_place_equal_numpy(tmp_path, compressed):
    make_run(tmp_path / "run", sets=("1",), compressed=compressed)
    run = rm.detect(tmp_path / "run")
    s = run.sets[0]
    assert s.stored is (not compressed)
    with np.load(s.path) as z:
        XF, G = z["XF"], z["genes"]
    assert np.array_equal(s.xf().rows(100, 250), XF[100:350])
    assert np.array_equal(s.genes().rows(s.n - 5, 50), G[-5:])  # past the end: cut
    assert np.array_equal(s.xf().rows(0, s.n), XF)


def test_gene_codes_ignore_string_width():
    codes = rm.GeneCodes(["Gad2", "blank0001"])
    a = codes.codes(np.array(["Gad2", "Mbp", "blank0001", "Mbp"], dtype="<U13"))
    b = codes.codes(np.array(["Gad2", "Mbp"], dtype="<U4"))
    assert codes.names == ["Gad2", "blank0001", "Mbp"]
    assert list(a) == [0, 2, 1, 2] and list(b) == [0, 2]
    assert list(codes.is_blank()) == [False, True, False]


def test_an_exact_tie_goes_to_the_lowest_field():
    centres = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 0.0]])  # field 2 re-imaged field 0
    xy = np.array([[1.0, 0.0], [9.0, 0.0], [5.0, 0.0]])
    assert list(rm.nearest_field(xy, centres)) == [0, 1, 0]


def test_cells_keep_one_copy_per_overlap(tmp_path):
    truth = make_run(tmp_path / "run")
    run = rm.detect(tmp_path / "run")
    panel, fields, problems = rm.read_cells(run, progress=lambda m: None)
    assert panel == PANEL and problems == []
    in_mats = sum(len(f.cells) for f in fields)
    kept = sum(int(f.keep.sum()) for f in fields)
    assert in_mats > truth["cells"]  # the overlaps are in two mats or more
    assert kept == truth["cells"]


def test_check_end_to_end(tmp_path):
    truth = make_run(tmp_path / "run", grid=3, bad_field=3)  # outliers are judged in sets of 5+ fields
    run = rm.detect(tmp_path / "run")
    out = rm.check(run, tmp_path / "work", ["slc17a7", "Gad1"], "2", "molecules_per_cell",
                   progress=lambda m: None)
    assert out["sets"] == 2 and out["fields"] == 18
    assert out["molecules"] == sum(truth["molecules"].values())
    assert out["cells"] == truth["cells"]
    assert out["panel"] == {"genes": 5, "blanks": 2}
    for k in ("overview", "fields", "genes", "molecules"):
        assert Path(out["plots"][k]).stat().st_size > 1000
    flags = " | ".join(out["flags"])
    assert "Conv_zscan1_03" in flags and "blank fraction" in flags  # the bad field
    assert "Conv_zscan2_03" in flags
    assert "Not in the panel: Gad1" in flags
    assert f"Set 2, molecules of the chosen genes: Slc17a7 {truth['per_gene']['2']['Slc17a7']:,}" in flags
    assert "not de-duplicated" not in flags
    assert all(p["outside_own_field"] == 0 for p in out["per_set"])
    assert json.loads((tmp_path / "work" / "summary.json").read_text())["cells"] == truth["cells"]


def test_spots_not_deduplicated_are_flagged(tmp_path):
    make_run(tmp_path / "run", sets=("1",), dedupe=False)
    out = rm.check(rm.detect(tmp_path / "run"), tmp_path / "w", [], None, "cells",
                   progress=lambda m: None)
    assert any("not de-duplicated" in f for f in out["flags"])


@pytest.mark.parametrize("again,flagged", [(0, False), (2, True)])
def test_a_field_imaged_twice(tmp_path, again, flagged):
    """The first fields imaged again at the end of a set is the lab's way
    (a note); any other field imaged twice is a flag."""
    make_run(tmp_path / "run", sets=("1",), grid=3)
    f = tmp_path / "run" / "final_spots" / "spots_allth3_final_set1.npz"
    with np.load(f) as z:
        d = dict(z)
    d["pos_fovs_pix"] = np.r_[d["pos_fovs_pix"], d["pos_fovs_pix"][again:again + 1] + [0, 3]]
    d["fovs"] = np.r_[d["fovs"], ["Conv_zscan1_09"]]
    np.savez(f, **d)
    out = rm.check(rm.detect(tmp_path / "run"), tmp_path / "w", [], None, "cells",
                   progress=lambda m: None)
    pair = f"Conv_zscan1_0{again} and Conv_zscan1_09"
    assert any(pair in x for x in out["flags"]) is flagged
    assert any("imaged again at the end" in x for x in out["notes"]) is not flagged


def test_fields_imaged_but_never_decoded_are_flagged(tmp_path):
    make_run(tmp_path / "run", sets=("1",), grid=3)
    (tmp_path / "run" / "POS").mkdir()
    (tmp_path / "run" / "POS" / "pos_1.txt").write_text("".join(f"{i}.0,0.0\n" for i in range(13)))
    run = rm.detect(tmp_path / "run")
    assert rm.undecoded_fields(run) == [
        "Set 1: 13 positions imaged (POS/pos_1.txt), 9 fields decoded — 4 never decoded (9–12): holes in the section"]
    out = rm.check(run, tmp_path / "w", [], None, "cells", progress=lambda m: None)
    assert any("never decoded (9–12)" in f for f in out["flags"])
    assert rm.field_index("Conv_zscan10_090") == 90


def test_missing_mats_are_flagged(tmp_path):
    make_run(tmp_path / "run", sets=("1",))
    (tmp_path / "run" / "final_mats" / "Conv_zscan1_02--mat.npz").unlink()
    out = rm.check(rm.detect(tmp_path / "run"), tmp_path / "w", [], None, "cells",
                   progress=lambda m: None)
    assert any("1 of 4 fields have no cell mat (Conv_zscan1_02)" in f for f in out["flags"])


def test_other_layouts_are_refused_with_a_reason(tmp_path):
    with pytest.raises(rm.NotARun, match="No final_spots/"):
        rm.detect(tmp_path)
    per_field = tmp_path / "roy" / "FinalDecodedSpots"
    per_field.mkdir(parents=True)
    np.savez(per_field / "saveDecodedFinalV1--MERFISH_Analysis--Conv_zscan__001--_set1.npz", XF=np.zeros((1, 17)))
    with pytest.raises(rm.NotARun, match="saveDecodedFinalV1"):
        rm.detect(tmp_path / "roy")
    run = tmp_path / "half"
    (run / "final_spots").mkdir(parents=True)
    np.savez(run / "final_spots" / "spots_set1.npz", XF=np.zeros((3, 5), np.float32),
             genes=np.array(["a", "b", "c"]))
    with pytest.raises(rm.NotARun, match="lacks pix_size, fov_size, pos_fovs_pix, fovs"):
        rm.detect(run)


def test_cli_refusal_is_an_error_line(tmp_path, capsys):
    assert rm.main(["check", "--run", str(tmp_path), "--work", str(tmp_path / "w")]) == 2
    assert capsys.readouterr().out.startswith("ERROR: No final_spots/")
