"""scripts/metadata/validate_metadata.py: the Lamin side's bionty check (docs/plans/METADATA.md
step 2). Needs bionty, so it runs in the meyes-lamin environment and is skipped elsewhere:

  conda activate meyes-lamin && python -m pytest tests/python/test_metadata_bionty.py -q

The first run downloads the pinned ontology tables (bionty caches them).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

pytest.importorskip("bionty.base")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "metadata"))
import validate_metadata as vm  # noqa: E402


def meta(**kw):
    m = {"organism": "mouse", "organism_id": "NCBITaxon:10090", "tissue": "brain", "tissue_id": "UBERON:0000955",
         "disease": "unknown", "disease_id": None, "meyes_sex": "male", "meyes_sex_id": "PATO:0000384",
         "timepoint": "91.0 Days", "timepoint_kind": "age", "timepoint_days": 91.0}
    m.update(kw)
    return m


def test_a_good_store_passes_and_gets_its_stage():
    m, block = vm.validate(meta())
    assert block["_outcome"] == "pass" and not block["errors"]
    assert (m["organism"], m["development_stage"], m["meyes_development_stage_id"]) == \
        ("mouse", "13-week-old stage", "MmusDv:0000159")   # not a mouse strain sharing the taxon id
    assert block["sources"]["tissue"] == "uberon 2024-08-07"


def test_names_become_terms_and_ids():
    m, block = vm.validate(meta(organism="Macaque", organism_id=None, tissue="Hippocampal Formation",
                                tissue_id=None, disease="Alzheimer disease", meyes_sex="Female",
                                meyes_sex_id=None, timepoint=None, timepoint_days=None))
    assert block["_outcome"] == "pass", block["errors"]
    assert (m["organism"], m["organism_id"]) == ("macaque", "NCBITaxon:9544")
    assert (m["tissue"], m["tissue_id"]) == ("hippocampal formation", "UBERON:0002421")
    assert m["disease_id"] == "MONDO:0004975" and m["meyes_sex_id"] == "PATO:0000383"


def test_prenatal_human_stage():
    m, _ = vm.validate(meta(organism="human", organism_id="NCBITaxon:9606", timepoint="20.0 Gestational weeks",
                            timepoint_kind="stage", timepoint_days=126.0))
    assert m["development_stage"] == "18th week post-fertilization stage"


def test_wrong_terms_fail_with_suggestions():
    _, block = vm.validate(meta(organism_id="NCBITaxon:9606", tissue="hipocampus", tissue_id=None))
    assert block["_outcome"] == "fail"
    assert any("NCBITaxon:10090, not NCBITaxon:9606" in e for e in block["errors"])
    assert any("hipocampus" in e and "closest" in e for e in block["errors"])


def test_unknown_tissue_is_incomplete_not_invalid():
    _, block = vm.validate(meta(tissue="unknown", tissue_id=None))
    assert block["_outcome"] == "incomplete" and not block["errors"]
