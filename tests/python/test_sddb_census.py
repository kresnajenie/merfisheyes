"""scripts/sddb_census.py on the SpatialData-DB roots in tests/data/sddb/roots
(no network: the parquet footers are not read)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ROOTS = ROOT / "tests" / "data" / "sddb" / "roots"
sys.path.insert(0, str(ROOT / "scripts"))

import sddb_census as census  # noqa: E402


def survey(uid: str, key: str) -> dict:
    a = {"uid": uid, "key": key, "size": 0, "n_files": 0, "description": "tier silver", "ulabels": ["CC BY 4.0"]}
    root = json.loads((ROOTS / f"{uid}.json").read_text())
    return census.survey(a, {"root": root, "status": 200, "parquet": {}})


def test_visium_targets_the_full_image_system():
    r = survey("fJSOO03VzCTYtD1o", "visium/x.zarr")
    c = r["coords"]
    sample = "Targeted_Visium_Human_Cerebellum_Neuroscience"
    assert r["technology"] == "Visium"
    assert not c["has_global"]
    assert (c["target"], c["rule"]) == (sample, "region identity")
    assert c["full_images"] == [f"images/{sample}_full_image"]
    assert c["downscaled_images"] == [f"images/{sample}_hires_image", f"images/{sample}_lowres_image"]
    assert c["unreachable"] == [] and c["samples"] == 1 and c["multi_entry"]
    assert r["tables"][0]["csc"] == "csc_counts"
    assert r["license_class"] == "commercial" and r["tier"] == "silver"


def test_visium_hd_default_table_and_its_csc():
    r = survey("haOzFbJ7Bc9x7NW6", "visiumhd/x.zarr")
    assert r["default_table"] == "square_008um"
    csc = {t["name"]: t["csc"] for t in r["tables"]}
    assert csc == {"square_002um": None, "square_008um": None, "square_016um": "csc_counts", "square_016um_hvg": None}
    assert r["coords"]["target"] == "Visium_HD_3prime_Human_Ovarian_Cancer_FF"
    assert any("square_008um" in f for f in r["validator"]["fail"])


def test_xenium_and_atera_use_global():
    for uid in ("0W4MBvLWKZXauWhZ", "jCLxdY1tWKmbfWOB", "gZwZDQTAM8vvm8sE"):
        c = survey(uid, "xenium/x.zarr")["coords"]
        assert (c["target"], c["rule"]) == ("global", "global"), uid
        assert c["unsupported"] == []


def test_shapes_and_points_unread_without_footers():
    r = survey("0W4MBvLWKZXauWhZ", "xenium/x.zarr")
    assert {s["kind"] for s in r["shapes"]} == {"unread"}
    assert [p["layout"] for p in r["points"]] == ["unread"]
    assert "unknown" not in r["profile"]
