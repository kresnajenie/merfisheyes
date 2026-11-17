#!/usr/bin/env python3
"""Count, for every public MERFISHEYES catalog row, the profile columns it can fill today.

    python scripts/metadata/inventory.py --rows DIR [--stores DIR] [--registered FILE] > out.md

The gate of step 1 of docs/plans/METADATA.md. Read-only: it writes nothing but stdout.

  --rows        one JSON file per catalog row, as GET /api/agent/v1/catalog/<id> returns it
                (curated, bil and community rows; SpatialData-DB rows are theirs already).
  --stores      folders <bil_code>/store/ holding stats.json, qc.json and catalog.json beside a
                reprocessed store (the BIL pipeline's output). A row with one is also counted
                as it will be once that store is what it opens.
  --registered  JSON {"<code>/<name>.zarr": [uid, size_bytes]} of store artifacts in Lamin
                eyes-catalog.

Each required and recommended column of docs/specs/schemas/metadata.schema.json ($defs/indexRow,
x-meyes-level) gets one state per row:

  have     the value is in the row or its files as it is
  derive   a fixed rule makes it from a value that is there (a spelling table, a parse, an
           ontology lookup of an unambiguous name) — no person needed
  fetch    the source's API has it (the BIL API's submission date), not kept today
  person   a value exists but turning it into the column needs a person (which instrument a
           "MERFISH" row ran on; a free-text tissue)
  missing  nothing to make it from

"Fills today" is have + derive. Columns the index build derives (x-meyes-level "derived") and
n/a columns are not counted; nor is meyes_thumbnail_url, which the agent API does not return.
"""
from __future__ import annotations

import argparse
import collections
import json
import re
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCHEMA = REPO / "docs" / "specs" / "schemas" / "metadata.schema.json"
NOT_SEEN = {"meyes_thumbnail_url"}

# Today's spellings -> what a rule can make of them. Unambiguous names only; anything else is
# "person". The ids are checked with bionty in step 2, not here.
ORGANISMS = {"mouse", "mice", "human", "rat", "macaque", "rhesus macaque", "marmoset",
             "drosophila melanogaster", "schmidtea mediterranea", "zebrafish"}
TISSUES = {"brain", "brainstem", "breast", "skin", "whole embryo", "embryo"}
PLATFORMS = {"xenium": "10x_xenium", "merscope": "vizgen_merscope", "stereo-seq": "stomics_stereo_seq",
             "cosmx": "bruker_cosmx", "visium": "10x_visium", "visium hd": "10x_visium_hd"}
INPUTS = {"merscope": "vizgen_merscope", "xenium": "10x_xenium", "visium": "10x_visium"}  # stats.json input
LICENSES = {"CC-BY-4.0", "CC-BY-NC-4.0", "CC-BY-SA-4.0", "CC-BY-NC-SA-4.0", "CC0-1.0"}
AGE = re.compile(r"^\s*([0-9.]+)\s*(days?|weeks?|months?|years?|gestational weeks?)\s*$", re.I)
DOI = re.compile(r"10\.[0-9]{4,9}/[^\s?#]+")


def levels() -> dict[str, str]:
    props = json.loads(SCHEMA.read_text())["$defs"]["indexRow"]["properties"]
    return {c: p["x-meyes-level"] for c, p in props.items()
            if p["x-meyes-level"] in ("required", "recommended") and c not in NOT_SEEN}


def states(row: dict, store: dict | None, reg: tuple | None) -> dict[str, str]:
    """store: {"stats", "qc", "catalog"} of the row's reprocessed store, or None."""
    m = row.get("metadata") or {}
    cat = (store or {}).get("catalog") or {}
    cm = cat.get("metadata") or {}
    st = (store or {}).get("stats") or {}
    qc = (store or {}).get("qc") or {}
    sd_entry = any(e.get("datasetType") == "spatialdata" for e in row.get("entries", []))
    has_store = bool(store) or sd_entry
    meta = {**m, **cm}

    def v(*vals):
        return next((x for x in vals if x not in (None, "", [], {})), None)

    species, tissue = v(row.get("species"), cat.get("species")), v(row.get("tissue"), cat.get("tissue"))
    platform = v(row.get("platform"), cat.get("platform"))
    # The instrument's own output names the technology; "MERFISH" alone does not (academic
    # MERFISH or Vizgen MERSCOPE?).
    tech = INPUTS.get(st.get("input")) or PLATFORMS.get((platform or "").lower())
    age, lic = v(meta.get("age")), v(row.get("license"), meta.get("license"))
    doi = DOI.search(v(row.get("publicationLink"), row.get("externalLink")) or "")
    code = v(row.get("bilCode"), cat.get("bilCode"))
    from_store = "derive" if sd_entry and not store else "have"   # a registered URL store is read, not filed

    def known(x, table):
        return "missing" if not x else ("derive" if x.lower() in table else "person")

    s = {}
    s["uid"] = "have" if reg else "missing"
    s["canonical_key"] = s["study_id"] = "derive"
    s["organism"] = s["organism_id"] = known(species, ORGANISMS)
    s["tissue"] = s["tissue_id"] = known(tissue, TISSUES)
    s["technology"] = s["vendor"] = s["assay"] = s["product"] = (
        "derive" if tech else ("person" if platform else "missing"))
    s["modality"] = s["modalities"] = "derive" if (tech or platform) else "missing"
    s["validation_status"] = "derive" if s["organism_id"] == s["tissue_id"] == "derive" else "missing"
    s["zarr_url"] = "have" if has_store else "missing"
    s["n_obs"] = "have" if (st.get("cells") or row.get("cells")) else ("derive" if sd_entry else "missing")
    s["n_features"] = "have" if (st.get("genes") or row.get("genes")) else ("derive" if sd_entry else "missing")
    s["published_at"] = "derive"
    s["title"] = "have" if row.get("title") else "missing"
    d = row.get("disease")
    s["disease"] = s["disease_id"] = "derive" if d else "missing"
    if age and AGE.match(age):
        s["timepoint"] = "have"
        s["timepoint_kind"] = s["timepoint_days"] = "derive"
        s["development_stage"] = s["meyes_development_stage_id"] = (
            "derive" if s["organism"] == "derive" else "person")
    else:
        for c in ("timepoint", "timepoint_kind", "timepoint_days", "development_stage", "meyes_development_stage_id"):
            s[c] = "person" if age else "missing"
    s["panel_size"] = "have" if (st.get("genes") or qc.get("n_genes_imaged") or row.get("genes")) else "missing"
    s["license"] = "have" if lic in LICENSES else ("person" if lic else "missing")
    s["total_counts"] = "have" if qc.get("total_counts") is not None else ("derive" if sd_entry else "missing")
    s["card_url"] = s["size_bytes"] = "have" if reg else "missing"
    s["biomaterial_type"] = "derive"
    s["hne_image"] = s["if_image"] = s["default_table"] = from_store if has_store else "missing"
    s["n_transcripts"] = ("have" if "molecules" in st else from_store) if has_store else "missing"
    s["publication_date"] = "fetch" if code else "missing"
    s["doi"] = "derive" if doi else "missing"
    s["publication_year"] = "fetch" if doi else "missing"
    s["publication_verify_status"] = "derive" if doi else "person"
    s["sample_id"] = "derive" if code else "person"
    s["dataset_url"] = "have" if row.get("externalLink") or cat.get("externalLink") else "missing"
    s["genotype"] = "derive" if meta.get("genotype") else "missing"
    s["obs_unit"] = s["observation_unit"] = s["observation_unit_source"] = (
        "derive" if any(e.get("kind") == "single_cell" for e in row.get("entries", [])) or store else "missing")
    s["segmentation_source"] = "derive" if st.get("input") else "person"
    s["meyes_id"] = s["meyes_source"] = s["meyes_page_url"] = "have"
    s["meyes_kind"] = "derive"
    s["meyes_viewer_url"] = "derive" if has_store else "missing"
    s["meyes_description"] = "have" if v(row.get("description"), cat.get("description")) else "missing"
    s["meyes_sex"] = s["meyes_sex_id"] = "derive" if meta.get("sex") else "missing"
    s["meyes_n_sections"] = s["meyes_sections"] = "have" if qc.get("n_samples") else from_store if has_store else "missing"
    s["meyes_dims"] = "have" if st.get("dims") else from_store if has_store else "missing"
    s["meyes_de_columns"] = "have" if "de_stats" in st else from_store if has_store else "missing"
    s["meyes_pipeline"] = "have" if st.get("pipeline") else "missing"
    # No pipeline writes it yet; any store can have it computed (matrix_hash.py)
    s["meyes_matrix_hash"] = "derive" if has_store else "missing"
    return s


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--rows", required=True)
    ap.add_argument("--stores")
    ap.add_argument("--registered")
    a = ap.parse_args()

    lv = levels()
    rows = [json.loads(p.read_text()) for p in sorted(Path(a.rows).glob("*.json"))]
    rows = [r for r in rows if r.get("source") in ("curated", "bil", "community")]
    stores = {}
    if a.stores:
        for d in Path(a.stores).glob("*/store"):
            stores[d.parent.name] = {k: json.loads((d / f"{k}.json").read_text())
                                     for k in ("stats", "qc", "catalog") if (d / f"{k}.json").exists()}
    reg = {}
    if a.registered:
        for key, val in json.loads(Path(a.registered).read_text()).items():
            reg[key.split("/")[0]] = val
    missing_cols = set(lv) - set(states({}, None, None))
    if missing_cols:
        print(f"no rule for: {sorted(missing_cols)}", file=sys.stderr)
        return 1

    req = [c for c, l in lv.items() if l == "required"]
    rec = [c for c, l in lv.items() if l == "recommended"]
    out, per_col = [], collections.defaultdict(collections.Counter)
    for r in rows:
        code = r.get("bilCode")
        today = states(r, None, None)
        later = states(r, stores[code], reg.get(code)) if code in stores else None
        for c, s in today.items():
            if c in lv:
                per_col[(r["source"], c)][s] += 1
        n = lambda s, cols: sum(s[c] in ("have", "derive") for c in cols)
        out.append({"id": r["id"], "source": r["source"], "code": code or "", "title": r["title"],
                    "req": n(today, req), "rec": n(today, rec),
                    "fetch": sum(today[c] == "fetch" for c in lv),
                    "person": sum(today[c] == "person" for c in lv),
                    "req_store": n(later, req) if later else None, "rec_store": n(later, rec) if later else None})

    print(f"# Inventory: profile columns each row can fill\n")
    print(f"{len(rows)} rows; profile: {len(req)} required, {len(rec)} recommended columns counted "
          f"({SCHEMA.relative_to(REPO)}). 'with store' = once the row opens its reprocessed store "
          f"({len(stores)} stores, {len(reg)} registered in Lamin).\n")
    print("| source | rows | required filled (median, min–max) | recommended filled (median, min–max) | rows with a store | required with store (median) | recommended with store (median) |")
    print("|---|---|---|---|---|---|---|")
    for src in ("bil", "curated", "community"):
        R = [o for o in out if o["source"] == src]
        if not R:
            continue
        rq, rc = [o["req"] for o in R], [o["rec"] for o in R]
        W = [o for o in R if o["req_store"] is not None]
        ws = (f"{statistics.median(o['req_store'] for o in W):g}", f"{statistics.median(o['rec_store'] for o in W):g}") if W else ("—", "—")
        print(f"| {src} | {len(R)} | {statistics.median(rq):g} ({min(rq)}–{max(rq)}) of {len(req)} | "
              f"{statistics.median(rc):g} ({min(rc)}–{max(rc)}) of {len(rec)} | {len(W)} | {ws[0]} | {ws[1]} |")
    print("\n## By column (today)\n")
    print("have / derive / fetch / person / missing, per source\n")
    print("| column | level | bil | curated | community |")
    print("|---|---|---|---|---|")
    for c in req + rec:
        cells = []
        for src in ("bil", "curated", "community"):
            k = per_col[(src, c)]
            cells.append(" / ".join(str(k[s]) for s in ("have", "derive", "fetch", "person", "missing")))
        print(f"| `{c}` | {lv[c]} | " + " | ".join(cells) + " |")
    print("\n## By row\n")
    print("| source | BIL code | title | required | recommended | with store (req / rec) | fetch | person |")
    print("|---|---|---|---|---|---|---|---|")
    for o in sorted(out, key=lambda o: (o["source"], o["code"], o["title"])):
        ws = f"{o['req_store']} / {o['rec_store']}" if o["req_store"] is not None else "—"
        title = o["title"].replace("|", "/")[:70]
        print(f"| {o['source']} | {o['code']} | {title} | {o['req']} | {o['rec']} | {ws} | {o['fetch']} | {o['person']} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
