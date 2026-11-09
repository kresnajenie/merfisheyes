"""The JSON files written beside a store (stats.json, qc.json, catalog.json,
mapping.json), checked by scripts/validate_sidecars.py against
docs/specs/schemas/ — the data the agent API reads (docs/specs/MCP-PROTOCOL.md).

The fixtures are two BIL datasets as the pipeline wrote them on 2026-11-05:
ace-dip-vat (MERSCOPE, four taxonomy levels) and ace-owl-cup (Xenium, no
cell types, most QC values null).

  python -m pytest tests/python/test_sidecars.py -q
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures" / "sidecars"
sys.path.insert(0, str(ROOT / "scripts"))

import validate_sidecars  # noqa: E402


def levels(report, level):
    return [(where, msg) for lv, where, msg in report.lines if lv == level]


@pytest.mark.parametrize("code", ["ace-dip-vat", "ace-owl-cup"])
def test_fixtures_pass(code):
    r = validate_sidecars.check(str(FIXTURES / code))
    assert not levels(r, "FAIL"), levels(r, "FAIL")
    assert not levels(r, "WARN"), levels(r, "WARN")


def test_missing_files_are_not_failures(tmp_path):
    r = validate_sidecars.check(str(tmp_path))
    assert not levels(r, "FAIL")
    assert {w for w, _ in levels(r, "WARN")} == set(validate_sidecars.FILES)


def test_catalog_with_only_a_code_passes(tmp_path):
    (tmp_path / "catalog.json").write_text(json.dumps({"bilCode": "ace-abc-def", "source": {}, "written_at": "2026-11-05T23:15:00+00:00"}))
    assert not levels(validate_sidecars.check(str(tmp_path)), "FAIL")


def mutate(tmp_path, name, change):
    folder = tmp_path / "ds"
    shutil.copytree(FIXTURES / "ace-dip-vat", folder)
    doc = json.loads((folder / name).read_text())
    change(doc)
    (folder / name).write_text(json.dumps(doc))
    return validate_sidecars.check(str(folder))


@pytest.mark.parametrize("name, change, where", [
    ("stats.json", lambda d: d.pop("store"), "stats.json"),
    ("stats.json", lambda d: d.update(cells=-1), "stats.json"),
    ("stats.json", lambda d: d.update(dims=4), "stats.json"),
    ("qc.json", lambda d: d.update(n_class="34"), "qc.json"),            # per-level column, wrong type
    ("qc.json", lambda d: d.update(mean_boot_class=1.5), "qc.json"),     # a bootstrap score is in [0, 1]
    ("qc.json", lambda d: d.update(mmc_levels="cluster; class"), "qc.json"),
    ("catalog.json", lambda d: d.pop("bilCode"), "catalog.json"),
    ("catalog.json", lambda d: d["source"].update(title="somewhere"), "catalog.json"),
    ("catalog.json", lambda d: d.update(bilCode="ace-dip-van"), "consistency"),  # disagrees with qc.json
    ("mapping.json", lambda d: d.update(cells="other.zarr"), "consistency"),     # not the store stats names
])
def test_violations_fail(tmp_path, name, change, where):
    fails = levels(mutate(tmp_path, name, change), "FAIL")
    assert any(w == where for w, _ in fails), fails


def test_level_without_count_warns(tmp_path):
    r = mutate(tmp_path, "qc.json", lambda d: d.pop("n_supertype"))
    assert any(w == "qc.json" and "supertype" in m for w, m in levels(r, "WARN"))
