#!/usr/bin/env python3
"""Write a store's metadata into its root zarr.json (docs/specs/metadata.md §4).

    python scripts/metadata/build_metadata.py OUT [--metadata FILE] [--print]

OUT is a prepared folder: the store mapping.json names (else the one *.zarr in it), with the
files the pipeline writes beside it — stats.json, and for BIL datasets qc.json and
catalog.json. --metadata FILE adds what only a person knows (organism, tissue, title, …), as
index column names (docs/specs/schemas/metadata.schema.json $defs/indexRow); it wins over
everything else.

Writes four root attributes and leaves every other one alone:

  metadata_general       SpatialData-DB's block: labels
  metadata_<technology>  their per-technology block, when the technology is known
  meyes_metadata         the facts with ontology ids, under the index's column names
  meyes_provenance       what wrote it, where each value came from, the matrix hash's recipe

No ontology check runs here: names map to ids only through the spelling tables of
docs/specs/conventions.json (metadata.terms), and validation.status is "not_validated".
The check with bionty is scripts/metadata/validate_metadata.py, on the Lamin side.

With --print the blocks are printed and nothing is written. Exit 0 written (or printed), 2
no store found, 3 the store's table cannot be read.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parent
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(HERE))

import matrix_hash  # noqa: E402

CONVENTIONS = json.loads((SCRIPTS.parent / "docs" / "specs" / "conventions.json").read_text())
META = CONVENTIONS["metadata"]
VERSION = "1"
SECTION_COLUMN = "_sample_id"   # combine_slices tags each cell with its section
OURS = ("metadata_general", "meyes_metadata", "meyes_provenance")

AGE = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?)\s*(gestational\s+weeks?|days?|weeks?|months?|years?)\s*$", re.I)
EDAY = re.compile(r"^\s*E\s*([0-9]+(?:\.[0-9]+)?)\s*$", re.I)
PDAY = re.compile(r"^\s*P\s*([0-9]+(?:\.[0-9]+)?)\s*$", re.I)
DAYS = {"day": 1.0, "week": 7.0, "month": 30.4375, "year": 365.25}


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


def find_store(out: Path) -> Path | None:
    name = read_json(out / "mapping.json").get("cells")
    if name and (out / name).is_dir():
        return out / name
    stores = sorted(p for p in out.glob("*.zarr") if p.is_dir())
    return stores[0] if len(stores) == 1 else None


def term(table: str, value) -> tuple[str | None, str | None]:
    """Our spelling → (term name, id) through conventions.json; unknown → (value, None)."""
    if not value:
        return None, None
    hit = META["terms"][table].get(str(value).strip().lower())
    return (hit[0], hit[1]) if hit else (str(value).strip(), None)


def timepoint(text) -> dict:
    """The source's age text → timepoint, timepoint_kind, timepoint_days (metadata.md §7.4)."""
    if not text:
        return {}
    text = str(text).strip()
    out = {"timepoint": text}
    if m := EDAY.match(text):
        return {**out, "timepoint_kind": "stage", "timepoint_days": float(m.group(1))}
    if m := PDAY.match(text):
        return {**out, "timepoint_kind": "age", "timepoint_days": float(m.group(1))}
    m = AGE.match(text)
    if not m:
        return {**out, "timepoint_kind": "other"}
    n, unit = float(m.group(1)), m.group(2).lower()
    if unit.startswith("gestational"):
        # Gestational weeks count from the last menstrual period, two weeks before fertilization
        return {**out, "timepoint_kind": "stage", "timepoint_days": round((n - 2) * 7, 4)}
    return {**out, "timepoint_kind": "age", "timepoint_days": round(n * DAYS[unit.rstrip("s")], 4)}


def genotype(text) -> str | None:
    if not text:
        return None
    m = re.fullmatch(r"\s*Strain_Name\s*=\s*'([^']*)'\s*", str(text))
    g = m.group(1) if m else str(text).strip()
    return re.sub(r"(?i)^c57bl/?6", "C57BL/6", g) or None


def doi(*texts) -> str | None:
    for t in texts:
        if t and (m := re.search(r"10\.[0-9]{4,9}/[^\s?#]+", str(t))):
            return m.group(0).rstrip(".")
    return None


def from_catalog(cat: dict) -> tuple[dict, dict]:
    """A catalog.json (one bil-completed.json item) → columns, and where each came from."""
    m = cat.get("metadata") or {}
    src = cat.get("source") or {}
    cols: dict = {}
    where: dict = {}

    def put(col, value, field):
        if value not in (None, "", []):
            cols[col] = value
            where[col] = src.get(field) or src.get("metadata") or "catalog.json"

    code = cat.get("bilCode")
    if code:
        key = f"{META['keyPrefixes']['bil']}:{code}"
        put("canonical_key", key, "bilCode")
        put("study_id", key, "bilCode")
        put("sample_id", code, "bilCode")
        put("meyes_bil_code", code, "bilCode")
    put("title", cat.get("title"), "title")
    put("meyes_description", cat.get("description"), "description")
    name, oid = term("organism", cat.get("species"))
    put("organism", name, "species")
    put("organism_id", oid, "species")
    name, tid = term("tissue", cat.get("tissue"))
    put("tissue", name, "tissue")
    put("tissue_id", tid, "tissue")
    tech = META["technologyFrom"]["platform"].get(str(cat.get("platform") or "").lower())
    put("technology", tech, "platform")
    put("dataset_url", cat.get("externalLink"), "externalLink")
    put("doi", doi(cat.get("publicationLink")), "publicationLink")
    for k, v in timepoint(m.get("age")).items():
        put(k, v, "metadata")
    sex, sid = term("sex", m.get("sex"))
    put("meyes_sex", sex.lower() if sex else None, "metadata")
    put("meyes_sex_id", sid, "metadata")
    put("genotype", genotype(m.get("genotype")), "metadata")
    lic = m.get("license")
    put("license", lic if lic in META["licenses"] else None, "metadata")
    put("meyes_investigator", m.get("investigator"), "metadata")
    put("meyes_co_investigators", m.get("coInvestigators"), "metadata")
    put("meyes_institution", m.get("institution"), "metadata")
    put("meyes_funding", m.get("funding"), "metadata")
    return cols, where


def from_store(store: Path, stats: dict, qc: dict) -> tuple[dict, dict]:
    """What the store itself says: counts, sections, the matrix hash."""
    import anndata as ad
    import prepare_folder

    cols: dict = {}
    hashinfo: dict = {}
    tables = store / "tables"
    if tables.is_dir() and any(tables.iterdir()):
        name = prepare_folder.default_table_name(store)
        adata = ad.read_zarr(str(tables / name))
        cols["default_table"] = name
        cols["n_obs"], cols["n_features"] = int(adata.n_obs), int(adata.n_vars)
        cols["panel_size"] = int(adata.n_vars)
        cols["meyes_kind"] = "single_cell"
        cols["obs_unit"], cols["observation_unit"], cols["observation_unit_source"] = "cells", "cell", "technology_map"
        digests = matrix_hash.cell_digests(adata.X, list(adata.var_names))
        cols["meyes_matrix_hash"], sketch = matrix_hash.combine(digests)
        import scipy.sparse as sp

        hashinfo = {"version": matrix_hash.VERSION, "matrix": "X", "genes": "symbol",
                    "n_cells": len(digests), "n_nonzero": int(sp.csr_matrix(adata.X).count_nonzero()),
                    "sketch": sketch}
        if SECTION_COLUMN in adata.obs:
            col = adata.obs[SECTION_COLUMN].astype(str).to_numpy()
            sections = []
            for s in sorted(set(col), key=lambda x: (len(x), x)):
                mine = [d for d, c in zip(digests, col) if c == s]
                sections.append({"id": s, "n_obs": len(mine), "meyes_matrix_hash": matrix_hash.combine(mine)[0]})
            cols["meyes_sections"] = sections
            cols["meyes_n_sections"] = len(sections)
        if "total_counts" not in qc:
            import numpy as np

            cols["total_counts"] = float(np.asarray(adata.X.sum()))
    else:
        cols["meyes_kind"] = "single_molecule"
    if stats.get("molecules"):
        cols["n_transcripts"] = int(stats["molecules"])
    if stats.get("molecule_genes"):
        cols["meyes_n_molecule_genes"] = int(stats["molecule_genes"])
    if stats.get("dims"):
        cols["meyes_dims"] = int(stats["dims"])
    if "de_stats" in stats:
        cols["meyes_de_columns"] = list(stats["de_stats"])
    images = [i.lower() for i in stats.get("images") or []]
    cols["hne_image"] = any(re.search(r"(^|_)(he|hne|h&e)(_|$)", i) for i in images)
    cols["if_image"] = bool(images) and not all(re.search(r"(^|_)(he|hne|h&e)(_|$)", i) for i in images)
    cols["ftu_annotation"] = False
    if qc.get("total_counts") is not None:
        cols["total_counts"] = float(qc["total_counts"])
    if qc.get("n_samples") and "meyes_n_sections" not in cols:
        cols["meyes_n_sections"] = int(qc["n_samples"])
    if qc.get("mmc_methods"):
        cols["meyes_cell_type_method"] = f"MapMyCells {qc['mmc_methods']}"
    if qc.get("mmc_levels"):
        cols["meyes_cell_type_levels"] = qc["mmc_levels"].split(",")
    if stats.get("pipeline"):
        cols["meyes_pipeline"] = stats["pipeline"]
    if stats.get("input"):
        cols["meyes_pipeline_input"] = stats["input"]
        tech = META["technologyFrom"]["input"].get(stats["input"])
        if tech:
            cols["technology"] = tech
            cols["segmentation_source"] = "vendor"
        elif stats["input"] in ("folder", "h5ad", "table"):
            cols["segmentation_source"] = "depositor"
    return cols, hashinfo


def complete(cols: dict) -> dict:
    """Fill what follows from what is there: technology → vendor/assay/product/modality,
    licence label, disease default, sex id, the label/id pairs of input values."""
    tech = META["technologies"].get(cols.get("technology") or "")
    if tech:
        cols.setdefault("vendor", tech["vendor"])
        cols.setdefault("assay", tech["assay"])
        cols.setdefault("product", tech["product"])
        cols.setdefault("modality", tech["modality"])
    cols.setdefault("modality", "RNA")
    cols.setdefault("modalities", [cols["modality"]])
    cols.setdefault("biomaterial_type", "Specimen from Organism")
    if cols.get("disease") in (None, ""):
        cols["disease"], cols["disease_id"] = "unknown", None
    elif not cols.get("disease_id"):
        cols["disease"], cols["disease_id"] = term("disease", cols["disease"])
    for col, table in (("organism", "organism"), ("tissue", "tissue")):
        if cols.get(col) and not cols.get(f"{col}_id"):
            cols[col], cols[f"{col}_id"] = term(table, cols[col])
    if cols.get("meyes_sex") and not cols.get("meyes_sex_id"):
        name, sid = term("sex", cols["meyes_sex"])
        cols["meyes_sex"], cols["meyes_sex_id"] = name, sid
    if cols.get("tissue") is None:
        cols["tissue"], cols["tissue_id"] = "unknown", None
    pub = cols.get("publication_verify_status")
    if not pub and cols.get("doi"):
        cols["publication_verify_status"] = "provisional"
    return cols


def general(cols: dict) -> dict:
    """metadata_general: SpatialData-DB's keys (conventions.json metadata.sddb.generalKeys),
    labels only (metadata.md §4)."""
    g = {}
    for col, key in META["sddb"]["generalKeys"].items():
        v = cols.get(col)
        if col == "license":
            v = (META["licenses"].get(v or "") or {}).get("label")
        elif col in ("hne_image", "if_image", "ftu_annotation"):
            v = bool(v)
        if v is not None or col == "disease":
            g[key] = v
    g["General Metadata Version"] = META["sddb"]["generalMetadataVersion"]
    return g


def build(out: Path, metadata: dict | None = None, log=print) -> tuple[Path, dict]:
    """The four blocks for OUT's store (not written). Deterministic: the same folder gives the
    same bytes (no timestamp — the pipeline's outputs are compared byte for byte)."""
    store = find_store(out)
    if not store:
        raise FileNotFoundError(f"no store in {out} (mapping.json names none, and not exactly one *.zarr)")
    stats, qc, cat = (read_json(out / f) for f in ("stats.json", "qc.json", "catalog.json"))
    cols, where = from_catalog(cat) if cat else ({}, {})
    facts, hashinfo = from_store(store, stats, qc)
    for k, v in facts.items():
        if k == "technology" and k in cols:
            continue   # what the platform said is replaced only by the input file
        cols[k] = v
        where.setdefault(k, "store")
    for k, v in (metadata or {}).items():
        cols[k] = v
        where[k] = "input"
    if not cols.get("canonical_key") and cols.get("meyes_id"):
        key = f"{META['keyPrefixes']['meyes']}:{cols['meyes_id']}"
        cols["canonical_key"] = cols.get("canonical_key") or key
    cols.setdefault("title", store.stem)
    cols = complete(cols)
    # The index's own columns are not facts about the data (metadata.md §2)
    for k in ("uid", "zarr_url", "card_url", "size_bytes", "published_at", "validation_status",
              "meyes_id", "meyes_source", "meyes_page_url", "meyes_viewer_url"):
        cols.pop(k, None)
    meyes = {"version": VERSION, **{k: cols[k] for k in sorted(cols)}}
    provenance = {
        "canonical_key": cols.get("canonical_key"),
        "study_id": cols.get("study_id"),
        "pipeline": stats.get("pipeline"),
        "input": stats.get("input"),
        "sources": {k: where[k] for k in sorted(where) if k in cols},
        "validation": {"status": "not_validated", "at": None, "sources": {}},
        "packages": stats.get("packages") or {},
    }
    if hashinfo:
        provenance["matrix_hash"] = hashinfo
    blocks = {"metadata_general": general(cols), "meyes_metadata": meyes, "meyes_provenance": provenance}
    block = META["sddb"]["technologyBlocks"].get(cols.get("technology") or "")
    if block:
        blocks[block["attribute"]] = {block["versionKey"]: block["version"]}
    return store, blocks


def write(store: Path, blocks: dict) -> None:
    """Replace our blocks in the root zarr.json's attributes; touch nothing else."""
    path = store / "zarr.json"
    root = json.loads(path.read_text())
    attrs = root.setdefault("attributes", {})
    technology_blocks = {b["attribute"] for b in META["sddb"]["technologyBlocks"].values()}
    if "sddb_provenance" not in attrs:   # never rewrite a store SpatialData-DB converted
        for k in list(attrs):
            if k in OURS or k in technology_blocks:
                del attrs[k]
    attrs.update(blocks)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(root, indent=2, ensure_ascii=False) + "\n")
    tmp.replace(path)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out", type=Path, help="the prepared folder (the store and its sidecars)")
    ap.add_argument("--metadata", type=Path, help="JSON of index columns a person supplies")
    ap.add_argument("--print", action="store_true", help="print the blocks, write nothing")
    a = ap.parse_args(argv)
    extra = None
    if a.metadata:
        try:
            extra = json.loads(a.metadata.read_text())
        except (OSError, ValueError) as e:
            print(f"error: cannot read {a.metadata}: {e}", file=sys.stderr)
            return 2
    try:
        store, blocks = build(a.out, extra)
    except FileNotFoundError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    except Exception as e:  # noqa: BLE001 — an unreadable table
        print(f"error: {type(e).__name__}: {e}", file=sys.stderr)
        return 3
    if "sddb_provenance" in json.loads((store / "zarr.json").read_text()).get("attributes", {}):
        print(f"error: {store.name} was converted by SpatialData-DB; its metadata is theirs", file=sys.stderr)
        return 2
    if a.print:
        print(json.dumps(blocks, indent=2, ensure_ascii=False, default=str))
        return 0
    write(store, blocks)
    m = blocks["meyes_metadata"]
    print(f"{store}: metadata written — {m.get('organism')} / {m.get('tissue')} / "
          f"{m.get('technology') or 'technology unknown'}, {m.get('n_obs')} cells, "
          f"{m.get('meyes_n_sections') or 1} section(s), {m.get('meyes_matrix_hash', '')[:16]}…")
    return 0


if __name__ == "__main__":
    sys.exit(main())
