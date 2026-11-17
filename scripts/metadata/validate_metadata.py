#!/usr/bin/env python3
"""Check a store's metadata against the ontologies with bionty, before it is registered in Lamin.

    python scripts/metadata/validate_metadata.py STORE.zarr [--write] [--json]

The Lamin side's gate (docs/plans/METADATA.md M3, step 2): runs where lamindb and bionty are
installed (the meyes-lamin environment), never in the pipeline. Reads the store's
meyes_metadata (scripts/metadata/build_metadata.py wrote it) and, against the ontology sources
and versions SpatialData-DB pins (docs/specs/conventions.json metadata.ontologies):

  organism            Ensembl names (release-112), NCBITaxon ids
  tissue              UBERON
  disease             MONDO, or `normal` (PATO:0000461), or `unknown`
  meyes_sex           PATO
  development_stage   HsapDv / MmusDv, filled from timepoint_days when empty (metadata.md §7.4)

A name is matched by name or synonym, case-insensitively, and replaced by the term's own name;
an id must exist and agree with the name. With --write the standardised values go back into
the store's root zarr.json (meyes_metadata, metadata_general's labels) with
meyes_provenance.validation = {status, at, bionty, sources, errors, notes}.

Exit 0 pass; 3 incomplete (organism or tissue unknown or missing: not publishable); 4 a term
is invalid (the message names the closest terms).
"""
from __future__ import annotations

import argparse
import datetime as dt
import difflib
import json
import sys
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONVENTIONS = json.loads((HERE.parent.parent / "docs" / "specs" / "conventions.json").read_text())
PINS = CONVENTIONS["metadata"]["ontologies"]
NORMAL = CONVENTIONS["metadata"]["terms"]["disease"]["normal"]   # ["normal", "PATO:0000461"]


@lru_cache(maxsize=None)
def table(kind: str, organism: str | None = None):
    """{lower name or synonym → (name, id)}, {id → name} for one pinned source."""
    import bionty.base as bb

    if kind == "organism":
        pin = PINS["organismNames"]
        df = bb.Organism(taxa="vertebrates", source=pin["source"], version=pin["version"]).to_dataframe()
    elif kind == "tissue":
        df = bb.Tissue(source=PINS["tissue"]["source"], version=PINS["tissue"]["version"]).to_dataframe()
    elif kind == "disease":
        df = bb.Disease(source=PINS["disease"]["source"], version=PINS["disease"]["version"]).to_dataframe()
    elif kind == "phenotype":
        df = bb.Phenotype(source=PINS["phenotype"]["source"], version=PINS["phenotype"]["version"]).to_dataframe()
    elif kind == "stage":
        pin = PINS["developmentStageHuman" if organism == "human" else "developmentStageMouse"]
        df = bb.DevelopmentalStage(organism=organism, source=pin["source"], version=pin["version"]).to_dataframe()
    else:
        raise ValueError(kind)
    df = df.reset_index()
    by_name: dict = {}
    by_id: dict = {}
    for row in df.itertuples():
        oid, name = getattr(row, "ontology_id"), getattr(row, "name")
        if not isinstance(oid, str) or not isinstance(name, str):
            continue
        by_id.setdefault(oid, name)
        by_name.setdefault(name.lower(), (name, oid))
    for row in df.itertuples():
        syn = getattr(row, "synonyms", None)
        if isinstance(syn, str):
            for s in syn.split("|"):
                by_name.setdefault(s.strip().lower(), (row.name, row.ontology_id))
    return by_name, by_id


def source_label(kind: str, organism: str | None = None) -> str:
    key = {"organism": "organismNames", "tissue": "tissue", "disease": "disease", "phenotype": "phenotype",
           "stage": "developmentStageHuman" if organism == "human" else "developmentStageMouse"}[kind]
    return f"{PINS[key]['source']} {PINS[key]['version']}"


def check(kind: str, name, oid, organism: str | None = None) -> tuple[str | None, str | None, str | None]:
    """(standard name, id, error). An id wins when it exists; it must agree with the name."""
    by_name, by_id = table(kind, organism)
    if oid:
        if oid not in by_id:
            return name, oid, f"{oid} is not a term of {source_label(kind, organism)}"
        std = by_id[oid]
        hit = by_name.get(str(name).strip().lower()) if name else None
        if hit and hit[1] != oid:
            return name, oid, f"{name!r} is {hit[1]}, not {oid} ({std})"
        # Several Ensembl entries share one NCBITaxon id (mouse strains): the name that
        # agrees with the id is kept, not whichever entry the id points to first
        return (hit[0] if hit else std), oid, None
    if not name:
        return None, None, None
    hit = by_name.get(str(name).strip().lower())
    if hit:
        return hit[0], hit[1], None
    close = difflib.get_close_matches(str(name).lower(), list(by_name), n=3, cutoff=0.75)
    hint = f"; closest: {', '.join(repr(by_name[c][0]) + ' ' + by_name[c][1] for c in close)}" if close else ""
    return name, None, f"{name!r} is not a term of {source_label(kind, organism)}{hint}"


def ordinal(n: int) -> str:
    return f"{n}{'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')}"


def stage_for(organism: str, kind: str | None, days: float | None) -> str | None:
    """The most specific HsapDv / MmusDv term for an age, or None (metadata.md §7.4)."""
    if organism not in ("human", "mouse") or days is None:
        return None
    by_name, _ = table("stage", organism)
    names = []
    if kind == "stage":
        if organism == "human":
            names.append(f"{ordinal(int(days // 7))} week post-fertilization stage")
    elif kind == "age":
        d = round(days, 4)
        if d == int(d):
            names.append(f"{int(d)}-day-old stage")
        if d % 7 == 0:
            names.append(f"{int(d // 7)}-week-old stage")
        months = d / 30.4375
        if abs(months - round(months)) < 1e-3:
            names.append(f"{int(round(months))}-month-old stage")
        years = d / 365.25
        if abs(years - round(years)) < 1e-3:
            names.append(f"{int(round(years))}-year-old stage")
    for n in names:
        if n.lower() in by_name:
            return by_name[n.lower()][0]
    return None


def validate(meta: dict) -> tuple[dict, dict]:
    """(standardised meyes_metadata, validation block)."""
    import bionty

    m = dict(meta)
    errors, notes, used = [], [], {}

    def run(kind, ncol, icol, organism=None):
        name, oid, err = check(kind, m.get(ncol), m.get(icol), organism)
        if err:
            errors.append(f"{ncol}: {err}")
        else:
            if name is not None:
                m[ncol], m[icol] = name, oid
                used[ncol] = source_label(kind, organism)

    run("organism", "organism", "organism_id")
    if m.get("tissue") not in (None, "unknown") or m.get("tissue_id"):
        run("tissue", "tissue", "tissue_id")
    disease = (m.get("disease") or "unknown")
    if disease == "unknown":
        m["disease"], m["disease_id"] = "unknown", None
    elif disease.lower() == NORMAL[0] or m.get("disease_id") == NORMAL[1]:
        m["disease"], m["disease_id"] = NORMAL
        used["disease"] = source_label("phenotype")
    else:
        run("disease", "disease", "disease_id")
    if m.get("meyes_sex") or m.get("meyes_sex_id"):
        run("phenotype", "meyes_sex", "meyes_sex_id")
    org = m.get("organism")
    if org in ("human", "mouse"):
        if m.get("development_stage") or m.get("meyes_development_stage_id"):
            run("stage", "development_stage", "meyes_development_stage_id", org)
        else:
            stage = stage_for(org, m.get("timepoint_kind"), m.get("timepoint_days"))
            if stage:
                m["development_stage"] = stage
                m["meyes_development_stage_id"] = table("stage", org)[0][stage.lower()][1]
                used["development_stage"] = source_label("stage", org)
            elif m.get("timepoint"):
                notes.append(f"development_stage: no {source_label('stage', org)} term for {m['timepoint']!r}")
    elif m.get("timepoint"):
        notes.append(f"development_stage: no development-stage ontology for {org or 'this organism'}")

    missing = [c for c in ("organism_id", "tissue_id") if not m.get(c)]
    status = "fail" if errors else ("incomplete" if missing else "pass")
    if missing and not errors:
        notes.append(f"not publishable until {', '.join(missing)} {'is' if len(missing) == 1 else 'are'} known")
    block = {
        "status": "pass" if status == "pass" else "fail",
        "at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "bionty": bionty.__version__,
        "sources": used,
        "errors": errors,
        "notes": notes,
    }
    block["_outcome"] = status
    return m, block


def write(store: Path, meta: dict, block: dict) -> None:
    """Put the standardised values and the verdict back into the root zarr.json."""
    path = store / "zarr.json"
    root = json.loads(path.read_text())
    attrs = root["attributes"]
    attrs["meyes_metadata"] = meta
    g = attrs.setdefault("metadata_general", {})
    for key, col in (("Organism", "organism"), ("Tissue", "tissue"), ("Disease", "disease"),
                     ("Development Stage", "development_stage")):
        if meta.get(col) is not None:
            g[key] = meta[col]
    attrs.setdefault("meyes_provenance", {})["validation"] = {k: v for k, v in block.items() if not k.startswith("_")}
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(root, indent=2, ensure_ascii=False) + "\n")
    tmp.replace(path)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("store", type=Path)
    ap.add_argument("--write", action="store_true", help="write the standardised values and the verdict into the store")
    ap.add_argument("--json", action="store_true", help="print the result as JSON")
    a = ap.parse_args(argv)
    try:
        attrs = json.loads((a.store / "zarr.json").read_text()).get("attributes", {})
    except (OSError, ValueError) as e:
        print(f"error: {a.store}/zarr.json: {e}", file=sys.stderr)
        return 2
    if "meyes_metadata" not in attrs:
        print(f"error: {a.store.name} has no meyes_metadata: run scripts/metadata/build_metadata.py first", file=sys.stderr)
        return 3
    meta, block = validate(attrs["meyes_metadata"])
    outcome = block.pop("_outcome")
    if a.write:
        write(a.store, meta, block)
    if a.json:
        print(json.dumps({"outcome": outcome, "validation": block, "meyes_metadata": meta}, indent=2, ensure_ascii=False))
    else:
        print(f"{a.store.name}: {outcome}")
        for k in ("organism", "tissue", "disease", "meyes_sex", "development_stage"):
            ident = {"organism": "organism_id", "tissue": "tissue_id", "disease": "disease_id",
                     "meyes_sex": "meyes_sex_id", "development_stage": "meyes_development_stage_id"}[k]
            print(f"  {k:18} {meta.get(k)!s:28} {meta.get(ident) or ''}")
        for e in block["errors"]:
            print(f"  ERROR {e}")
        for n in block["notes"]:
            print(f"  note  {n}")
    return {"pass": 0, "incomplete": 3, "fail": 4}[outcome]


if __name__ == "__main__":
    sys.exit(main())
