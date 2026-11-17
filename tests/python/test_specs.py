"""The format specifications, checked from outside the viewer:

- scripts/validate_store.py passes the tiny fixtures composed into one
  store and fails each profile violation it is meant to catch;
- scripts/legacy_to_points.py turns a legacy molecule folder into a points
  element the validator and the reference reader accept, molecule for
  molecule;
- the spatialdata library reads our table and image, and a store it
  writes passes our validator (skipped when spatialdata is not installed).

  python -m pytest tests/python -q
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
import zarr

ROOT = Path(__file__).resolve().parents[2]
TINY = ROOT / "tests" / "data" / "tiny"
ZARR = TINY / "zarr"
sys.path.insert(0, str(ROOT / "scripts"))

import legacy_to_points  # noqa: E402
import points_zarr as pz  # noqa: E402
import validate_store  # noqa: E402


def compose_store(dst: Path, *, points=True, shapes=True) -> Path:
    """The tiny fixtures as one store, the way the E2E server serves them."""
    store = dst / "store.zarr"
    shutil.copytree(ZARR / "spatialdata_v3.zarr", store)
    if points:
        shutil.copytree(ZARR / "points_v1.zarr" / "points", store / "points")
    if shapes:
        shutil.copytree(ZARR / "shapes_v1", store / "shapes")
        # The fixture holds the element only; the store needs the group
        (store / "shapes" / "zarr.json").write_text(json.dumps({"zarr_format": 3, "node_type": "group"}))
    zarr.consolidate_metadata(str(store))
    return store


def levels(report: validate_store.Report):
    return {(lv, w) for lv, w, _ in report.lines}


def fails(report: validate_store.Report):
    return [m for lv, _, m in report.lines if lv == "FAIL"]


# ── validate_store ────────────────────────────────────────────────────────

def test_composed_tiny_store_passes(tmp_path):
    store = compose_store(tmp_path)
    r = validate_store.validate(store, sample=1000)
    assert not r.failed, fails(r)
    assert ("OK", "tables/table") in levels(r)
    assert ("OK", "points/transcripts") in levels(r)
    assert ("OK", "shapes/cell_boundaries") in levels(r)
    assert ("OK", "images/dapi") in levels(r)
    assert validate_store.main([str(store), "--sample", "1000"]) == 0


def test_conventions_file_matches_its_schema():
    """docs/specs/conventions.json, shared by the viewer and this validator."""
    import jsonschema

    schema = json.loads((ROOT / "docs" / "specs" / "schemas" / "conventions.schema.json").read_text())
    jsonschema.Draft202012Validator(schema).validate(validate_store.CONVENTIONS)
    # A broken file is refused
    bad = json.loads(json.dumps(validate_store.CONVENTIONS))
    bad["shapes"]["pointSiblingSuffix"] = ""
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema).validate(bad)
    # And the validator reads it: the default table follows the file
    named = validate_store.CONVENTIONS["tables"]["curatedDefaultTable"]
    attrs = {named["attribute"]: {named["key"]: "b"}}
    assert validate_store.default_table(attrs, ["a", "b"]) == "b"


# SpatialData-DB's REQUIRED_COLUMNS (src/sddb/_catalog_schema.py, schema 3).
SDDB_REQUIRED = ["uid", "technology", "assay", "organism", "validation_status", "zarr_url"]


def metadata_example() -> dict:
    """Root attributes of the reprocessed ace-dip-van store, as docs/specs/metadata.md §4 lays them out
    (values from its catalog.json and stats.json; section ids abbreviated)."""
    return {
        "metadata_general": {
            "Product": "MERSCOPE", "Assay": "Vizgen MERSCOPE", "Biomaterial Type": "Specimen from Organism",
            "Organism": "mouse", "Tissue": "brain", "Modality": "RNA", "Disease": "unknown",
            "Default Table": "table", "Dataset Url": "https://api.brainimagelibrary.org/web/view?bildid=ace-dip-van",
            "Sample ID at Source": "ace-dip-van", "License": "CC BY-NC 4.0", "Development Stage": "13-week-old stage",
            "HnE Image Available": False, "IF Image Available": False, "FTU Annotation Available": False,
            "Number of Transcripts": 472686422, "General Metadata Version": "0.3.2",
        },
        "metadata_merscope": {"Merscope Metadata Version": "0.1.0"},
        "meyes_metadata": {
            "version": "1", "canonical_key": "bil:ace-dip-van", "study_id": "bil:ace-dip-van",
            "title": "MERFISH Imaging of mesocorticolimbic regions in the mouse brain (infected with EcoHIV virus; systemic)",
            "organism": "mouse", "organism_id": "NCBITaxon:10090",
            "tissue": "brain", "tissue_id": "UBERON:0000955",
            "disease": "unknown", "disease_id": None,
            "development_stage": "13-week-old stage", "meyes_development_stage_id": "MmusDv:0000159",
            "technology": "vizgen_merscope", "vendor": "vizgen", "assay": "Vizgen MERSCOPE",
            "n_obs": 243224, "n_features": 300, "n_transcripts": 472686422,
            "timepoint": "91.0 Days", "timepoint_kind": "age", "timepoint_days": 91.0,
            "genotype": "C57BL/6J", "condition": "EcoHIV infection, systemic",
            "meyes_sex": "male", "meyes_sex_id": "PATO:0000384", "license": "CC-BY-NC-4.0",
            "meyes_kind": "single_cell", "meyes_n_sections": 3,
            "meyes_dims": 2, "meyes_sections": [{"id": "s1"}, {"id": "s2"}, {"id": "s3"}],
            # computed by scripts/metadata/matrix_hash.py on the store, 2026-11-16
            "meyes_matrix_hash": "mh1:d0ff37d6072dd0d8ef11bae317b5fa6e573f59364d40c61962dd183432dfb025",
        },
        "meyes_provenance": {
            "canonical_key": "bil:ace-dip-van", "pipeline": "93ab1c0", "input": "merscope",
            "written_at": "2026-11-16T00:00:00Z",
            "validation": {"status": "pass", "at": "2026-11-16T00:00:00Z",
                           "sources": {"organism": "ncbitaxon 2023-06-20", "tissue": "uberon 2024-08-07"}},
        },
    }


def test_metadata_schema_profiles_sddb_schema_3():
    """docs/specs/schemas/metadata.schema.json: every index column has a level, SpatialData-DB's
    required columns stay required, the pinned ontologies are the ones the spec names, and a
    store's root attributes laid out as the spec says validate — a wrong id does not."""
    import jsonschema

    schema = json.loads((ROOT / "docs" / "specs" / "schemas" / "metadata.schema.json").read_text())
    jsonschema.Draft202012Validator.check_schema(schema)
    row = schema["$defs"]["indexRow"]
    cols = row["properties"]
    assert all(p["x-meyes-level"] in ("required", "recommended", "optional", "derived", "n/a")
               for p in cols.values())
    sddb = [c for c in cols if not c.startswith("meyes_")]
    assert len(sddb) == 93 and all("x-sddb-dtype" in cols[c] for c in sddb)
    # Their required columns are ours too: now, or once the fleet is reprocessed (METADATA.md M10)
    assert set(SDDB_REQUIRED) <= set(row["required"]) | set(row["x-meyes-required-after-reprocessing"])
    pins = validate_store.CONVENTIONS["metadata"]
    assert pins["sddb"]["catalogSchemaVersion"] == 3
    assert schema["$defs"]["metadataGeneral"]["properties"]["General Metadata Version"]["const"] == \
        pins["sddb"]["generalMetadataVersion"]

    v = jsonschema.Draft202012Validator(schema)
    good = metadata_example()
    v.validate(good)
    bad = json.loads(json.dumps(good))
    bad["meyes_metadata"]["tissue_id"] = "brain"
    with pytest.raises(jsonschema.ValidationError):
        v.validate(bad)
    no_id = json.loads(json.dumps(good))
    del no_id["meyes_metadata"]["organism_id"]
    with pytest.raises(jsonschema.ValidationError):
        v.validate(no_id)


def test_matrix_hash_names_the_same_data_the_same():
    """scripts/metadata/matrix_hash.py (mh1, docs/specs/metadata.md §7.7): storage, dtype,
    cell order and gene order do not change the hash; one count does; the sketch sees
    shared cells; the value of a fixed matrix is pinned, so the algorithm cannot drift."""
    import scipy.sparse as sp
    sys.path.insert(0, str(ROOT / "scripts" / "metadata"))
    import matrix_hash as mh

    genes = ["Gad1", "Slc17a7", "Sst"]
    m = np.array([[0, 3, 1], [2, 0, 0], [0, 0, 0], [0.5, 1, 0]])
    h = mh.matrix_hash(m, genes)
    assert h == "mh1:359b1643cfa569e53c2ce4f02614079bfbb88c7aea83a68045b78772cd25b9b9"
    assert mh.matrix_hash(sp.csc_matrix(m.astype(np.float32)), genes) == h
    assert mh.matrix_hash(sp.csr_matrix(m), [" GAD1", "slc17a7", "SST "]) == h
    cells, cols = [3, 0, 2, 1], [2, 0, 1]
    assert mh.matrix_hash(m[cells][:, cols], [genes[i] for i in cols]) == h
    with_empty_gene = np.hstack([m, np.zeros((4, 1))])
    assert mh.matrix_hash(with_empty_gene, genes + ["Vip"]) == h
    changed = m.copy()
    changed[1, 0] = 3
    assert mh.matrix_hash(changed, genes) != h
    with pytest.raises(ValueError):
        mh.matrix_hash(m, ["Gad1", "GAD1", "Sst"])

    rng = np.random.default_rng(0)
    big = sp.random(400, 50, density=0.2, random_state=1, data_rvs=lambda n: rng.integers(1, 20, n))
    names = [f"g{i}" for i in range(50)]
    d = mh.cell_digests(big, names)
    _, whole = mh.combine(d)
    _, half = mh.combine(d[:200] + mh.cell_digests(sp.random(200, 50, density=0.2, random_state=2), names))
    assert mh.sketch_overlap(whole, whole) == 1.0
    assert 0.1 < mh.sketch_overlap(whole, half) < 0.6   # true Jaccard 200/600 = 0.33


def test_build_metadata_writes_the_four_blocks(tmp_path):
    """scripts/metadata/build_metadata.py (docs/specs/metadata.md §4): from a BIL-shaped
    catalog.json and the store, the root attributes gain metadata_general (labels),
    metadata_merscope, meyes_metadata (ids, counts, hash) and meyes_provenance; they pass the
    schema; spatialdata_attrs and the consolidated metadata are untouched; validate_store
    warns until bionty has checked the terms, fails with --require-metadata, and fails when
    the metadata no longer matches the table."""
    sys.path.insert(0, str(ROOT / "scripts" / "metadata"))
    import build_metadata

    store = compose_store(tmp_path)
    (tmp_path / "mapping.json").write_text(json.dumps({"cells": store.name}))
    (tmp_path / "stats.json").write_text(json.dumps({
        "input": "merscope", "store": store.name, "molecules": 1234, "dims": 2, "de_stats": [],
        "pipeline": "abc1234", "packages": {"anndata": "0.13.4"}}))
    (tmp_path / "catalog.json").write_text(json.dumps({
        "bilCode": "ace-xyz-abc", "title": "A tiny mouse brain", "species": "Mouse", "tissue": "Brain",
        "platform": "MERFISH", "externalLink": "https://api.brainimagelibrary.org/web/view?bildid=ace-xyz-abc",
        "metadata": {"age": "20.0 Gestational weeks", "sex": "Female", "genotype": "Strain_Name='C57Bl/6J'",
                     "license": "CC-BY-SA-4.0", "investigator": "Doe, Jane"},
        "source": {"title": "bil-api", "species": "bil-api", "metadata": "bil-api"}}))
    before = json.loads((store / "zarr.json").read_text())
    r = validate_store.validate(store)
    assert ("WARN", "metadata") in levels(r) and not r.failed

    s, blocks = build_metadata.build(tmp_path)
    build_metadata.write(s, blocks)
    root = json.loads((store / "zarr.json").read_text())
    attrs = root["attributes"]
    assert attrs["spatialdata_attrs"] == before["attributes"]["spatialdata_attrs"]
    assert root["consolidated_metadata"] == before["consolidated_metadata"]
    g, m, p = attrs["metadata_general"], attrs["meyes_metadata"], attrs["meyes_provenance"]
    assert attrs["metadata_merscope"] == {"Merscope Metadata Version": "0.1.0"}
    assert (g["Organism"], g["Tissue"], g["Assay"], g["License"]) == ("mouse", "brain", "Vizgen MERSCOPE", "CC BY-SA 4.0")
    assert (m["organism_id"], m["tissue_id"], m["technology"], m["vendor"]) == \
        ("NCBITaxon:10090", "UBERON:0000955", "vizgen_merscope", "vizgen")
    assert (m["canonical_key"], m["sample_id"], m["genotype"]) == ("bil:ace-xyz-abc", "ace-xyz-abc", "C57BL/6J")
    assert (m["timepoint_kind"], m["timepoint_days"]) == ("stage", 126.0)   # GW20 = 18 weeks post-fertilization
    assert (m["meyes_sex"], m["meyes_sex_id"], m["disease"]) == ("female", "PATO:0000383", "unknown")
    assert m["meyes_matrix_hash"].startswith("mh1:") and m["n_transcripts"] == 1234
    assert p["validation"]["status"] == "not_validated" and p["sources"]["organism"] == "bil-api"
    assert not validate_store.schema_errors(
        {k: v for k, v in attrs.items() if k.startswith(("metadata_", "meyes_"))}, "metadata.schema.json")

    r = validate_store.validate(store)
    assert not r.failed and ("OK", "metadata") in levels(r)
    assert validate_store.validate(store, require_metadata=True).failed   # not validated yet
    attrs["meyes_metadata"]["n_obs"] += 1
    (store / "zarr.json").write_text(json.dumps(root))
    assert any("stale" in f for f in fails(validate_store.validate(store)))


def test_modality_fixtures_pass():
    """The stores shaped like SpatialData-DB's CODEX, MALDI-MSI and CosMx ones
    (scripts/testdata/make-tiny-modalities.py): a dense X, a feature index
    column that is not `_index`, an image placed by a rotating affine, and a
    table whose cells' positions are only in its shapes."""
    codex = validate_store.validate(ZARR / "codex_v3.zarr")
    assert not codex.failed, fails(codex)
    assert any(w == "tables/table" and "X is array" in m for lv, w, m in codex.lines)

    maldi = validate_store.validate(ZARR / "maldi_v3.zarr")
    assert not maldi.failed, fails(maldi)
    assert any("placed by an affine" in m and "rotated 30.0°" in m for _, _, m in maldi.lines)

    cosmx = validate_store.validate(ZARR / "cosmx_v3.zarr")
    assert not cosmx.failed, fails(cosmx)
    assert any("positions come from shapes/cells_polygons_for_vitessce" in m for _, _, m in cosmx.lines)


def test_table_without_positions_anywhere_fails(tmp_path):
    """No obsm/spatial and no region shapes element: the validator says so."""
    store = tmp_path / "store.zarr"
    shutil.copytree(ZARR / "cosmx_v3.zarr", store)
    shutil.rmtree(store / "shapes" / "cells_polygons")
    shutil.rmtree(store / "shapes" / "cells_polygons_for_vitessce")
    zarr.consolidate_metadata(str(store))
    r = validate_store.validate(store)
    assert any("no obsm/spatial, and no region shapes element" in m for m in fails(r))


def test_mapping_beside_the_store(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=False)
    (tmp_path / "mapping.json").write_text(json.dumps({"cells": "store.zarr", "image": "https://x/dapi.zarr"}))
    r = validate_store.validate(store)
    assert ("OK", "mapping.json") in levels(r)
    (tmp_path / "mapping.json").write_text(json.dumps({"links": {"__all__": "u"}}))
    r = validate_store.validate(store)
    assert any("linkColumn" in m for m in fails(r))


def test_missing_consolidated_metadata_fails(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=False)
    root = json.loads((store / "zarr.json").read_text())
    del root["consolidated_metadata"]
    (store / "zarr.json").write_text(json.dumps(root))
    r = validate_store.validate(store)
    assert any("consolidated" in m for m in fails(r))


def test_stale_consolidated_metadata_fails(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=True)
    # Written after the last consolidation: the viewer would never see it
    root = json.loads((store / "zarr.json").read_text())
    root["consolidated_metadata"]["metadata"] = {
        k: v for k, v in root["consolidated_metadata"]["metadata"].items() if not k.startswith("shapes")
    }
    (store / "zarr.json").write_text(json.dumps(root))
    r = validate_store.validate(store)
    assert any("misses" in m for m in fails(r))


def test_shapes_as_a_directory_fails(tmp_path):
    store = compose_store(tmp_path, points=False)
    f = store / "shapes" / "cell_boundaries" / "shapes.parquet"
    part = f.read_bytes()
    f.unlink()
    f.mkdir()
    (f / "part.0.parquet").write_bytes(part)
    r = validate_store.validate(store)
    assert any("directory of parts" in m for m in fails(r))


def test_csr_without_csc_copy_fails(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=False)
    x = store / "tables" / "table" / "X" / "zarr.json"
    meta = json.loads(x.read_text())
    meta.setdefault("attributes", {})["encoding-type"] = "csr_matrix"
    x.write_text(json.dumps(meta))
    # The tiny table is CSR with the copy the profile asks for: accepted
    assert not validate_store.validate(store).failed
    shutil.rmtree(store / "tables" / "table" / "layers" / "X_csc")
    zarr.consolidate_metadata(str(store))
    r = validate_store.validate(store)
    assert any("CSR" in m for m in fails(r))


SDDB_ROOTS = ROOT / "tests" / "data" / "sddb" / "roots"


def test_root_json_reads_a_store_from_its_root_alone():
    # Visium (cerebellum): CSR X with a csc_counts layer, as the reader accepts it
    root = validate_store.read_json(SDDB_ROOTS / "fJSOO03VzCTYtD1o.json")
    r = validate_store.validate(validate_store.RootStore(root))
    assert not r.failed, fails(r)
    assert any("layers/csc_counts" in m for lv, w, m in r.lines if w == "tables/table")
    assert any("root only" in m for lv, w, m in r.lines if w.startswith("shapes/"))
    assert validate_store.main(["--root-json", str(SDDB_ROOTS / "fJSOO03VzCTYtD1o.json")]) == 0


def test_root_json_fails_a_csr_table_with_no_csc_layer():
    # Visium HD (ovary): only the 16 um table has a CSC layer
    root = validate_store.read_json(SDDB_ROOTS / "haOzFbJ7Bc9x7NW6.json")
    r = validate_store.validate(validate_store.RootStore(root))
    failed = {w for lv, w, _ in r.lines if lv == "FAIL"}
    assert failed == {"tables/square_002um", "tables/square_008um", "tables/square_016um_hvg"}
    assert ("OK", "tables/square_016um") in levels(r)


def _root_with(nodes: dict) -> dict:
    return {
        "zarr_format": 3,
        "node_type": "group",
        "attributes": {"spatialdata_attrs": {"version": "0.2"}},
        "consolidated_metadata": {"kind": "inline", "must_understand": False, "metadata": nodes},
    }


def _scale(s, out="global", axes=("x", "y")):
    return {"type": "scale", "scale": list(s), "input": {"name": "xy", "axes": [{"name": a} for a in axes]}, "output": {"name": out}}


def _identity(out="global"):
    return {"type": "identity", "output": {"name": out}}


def test_coordinates_points_must_share_the_cells_frame():
    # As our converted Atera store was written: cells scaled into global,
    # the molecules declaring identity (their data is in the cells' frame)
    cells = {"attributes": {"coordinateTransformations": [_scale([4.70588, 4.70588])], "axes": ["x", "y"]}}
    table = {"attributes": {"region": "cell_circles", "instance_key": "cell_id"}}
    root = _root_with({
        "tables/table": table,
        "shapes/cell_circles": cells,
        "points/transcripts": {"attributes": {"coordinateTransformations": [_identity()], "axes": ["x", "y"]}},
    })
    r = validate_store.validate(validate_store.RootStore(root))
    assert any(w == "points/transcripts" and "frame other than the cells'" in m for lv, w, m in r.lines if lv == "FAIL")
    # The same transform as the cells: fine
    root["consolidated_metadata"]["metadata"]["points/transcripts"]["attributes"]["coordinateTransformations"] = [
        _scale([4.70588, 4.70588])
    ]
    r = validate_store.validate(validate_store.RootStore(root))
    assert not any(w == "points/transcripts" for lv, w, _ in r.lines if lv == "FAIL")


def test_coordinates_visium_targets_the_spots_and_places_an_affine():
    root = validate_store.read_json(SDDB_ROOTS / "fJSOO03VzCTYtD1o.json")
    cs = validate_store.coordinate_systems(validate_store.RootStore(root))
    sample = "Targeted_Visium_Human_Cerebellum_Neuroscience"
    assert (cs["target"], cs["rule"], cs["region"]) == (sample, "region identity", f"shapes/{sample}")
    assert cs["place"](f"images/{sample}_full_image")[0] == validate_store.IDENTITY
    m, _ = cs["place"](f"images/{sample}_hires_image")
    assert round(m[0], 3) == 6.666
    # An affine image is placed, rotation and all, and the report says so
    nodes = root["consolidated_metadata"]["metadata"]
    ms = nodes[f"images/{sample}_hires_image"]["attributes"]["ome"]["multiscales"][0]
    ms["coordinateTransformations"] = [{
        "type": "affine",
        "affine": [[1, 0, 0, 0], [0, 0, 2, 5], [0, -2, 0, 7]],
        "input": {"name": "cyx", "axes": [{"name": "c"}, {"name": "y"}, {"name": "x"}]},
        "output": {"name": sample, "axes": [{"name": "c"}, {"name": "y"}, {"name": "x"}]},
    }]
    cs = validate_store.coordinate_systems(validate_store.RootStore(root))
    m, _ = cs["place"](f"images/{sample}_hires_image")
    assert [round(v, 6) for v in m] == [0.0, -2.0, 2.0, 0.0, 7.0, 5.0]
    r = validate_store.validate(validate_store.RootStore(root))
    assert any("rotated" in m for lv, w, m in r.lines if lv == "OK" and w.endswith("hires_image"))
    # A broken affine (no rows) still reads as the identity, not a crash
    ms["coordinateTransformations"] = [{"type": "affine", "affine": [], "output": {"name": sample}}]
    cs = validate_store.coordinate_systems(validate_store.RootStore(root))
    assert cs["place"](f"images/{sample}_hires_image")[0] == validate_store.IDENTITY


def test_maldi_and_cosmx_roots_pass():
    """The real roots (tests/data/sddb/roots): MALDI's affine image is placed,
    CosMx's cells take their positions from the shapes."""
    maldi = validate_store.validate(validate_store.RootStore(validate_store.read_json(SDDB_ROOTS / "fF2M7WUbkFnMMhpC.json")))
    assert not maldi.failed, fails(maldi)
    assert any("rotated -90.9°" in m for _, _, m in maldi.lines)
    cosmx = validate_store.validate(validate_store.RootStore(validate_store.read_json(SDDB_ROOTS / "571yTEm4hqlkXRPf.json")))
    assert not cosmx.failed, fails(cosmx)
    assert any("positions come from shapes/cells_polygons_for_vitessce" in m for _, _, m in cosmx.lines)


def test_points_writer_keeps_the_source_transforms(tmp_path):
    out = tmp_path / "store.zarr" / "points" / "transcripts"
    w = pz.PointsWriter(
        out, tmp_path / "tmp", features=["a"], dims=2, origin=(0, 0), extent=(10, 10), n_levels=1,
        transformations=[_scale([4.70588, 4.70588])],
    )
    w.spill({"x": np.array([1.0, 2.0], np.float32), "y": np.array([1.0, 2.0], np.float32), "f": np.zeros(2, np.uint16)}, seed=0, tag="t")
    w.finish()
    _, attrs = pz.open_element(out)
    [t] = attrs["coordinateTransformations"]
    assert (t["type"], t["scale"], t["output"]["name"]) == ("scale", [4.70588, 4.70588], "global")
    with pytest.raises(ValueError):
        pz.PointsWriter(out, tmp_path / "t2", features=["a"], dims=3, origin=(0, 0), extent=(1, 1),
                        transformations=[_scale([1, 1])])


def test_image_axes_and_dtype_are_checked(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=False)
    g = store / "images" / "dapi" / "zarr.json"
    meta = json.loads(g.read_text())
    meta["attributes"]["ome"]["multiscales"][0]["axes"] = [{"name": "y"}, {"name": "x"}]
    g.write_text(json.dumps(meta))
    r = validate_store.validate(store)
    assert any("axes" in m for m in fails(r))


def test_ngff_04_attributes_in_a_v3_group_are_read(tmp_path):
    """The layout older ome-zarr-py wrote (and still reads): `multiscales`
    and `omero` at the top of a Zarr v3 group's attributes, no `ome` wrapper.
    Read, warned about, never written — image-profile.md, and
    lib/spatialdata/ome.ts for the viewer's copy of the rule."""
    store = compose_store(tmp_path, points=False, shapes=False)
    g = store / "images" / "dapi" / "zarr.json"
    meta = json.loads(g.read_text())
    ome = meta["attributes"].pop("ome")
    ome["multiscales"][0]["version"] = "0.4"
    meta["attributes"] = {k: ome[k] for k in ("multiscales", "omero") if k in ome}
    g.write_text(json.dumps(meta))

    r = validate_store.validate(store)

    assert ("OK", "images/dapi") in levels(r), "the image should still validate"
    assert not any("ome.multiscales" in m for m in fails(r))
    assert any("not under `ome`" in m for _, _, m in r.lines), "it should be flagged"
    # The same metadata, through the helper the readers share
    assert validate_store.image_multiscale(meta["attributes"])["axes"]
    assert validate_store.ome_attrs(meta["attributes"]).get("omero")


def test_points_metadata_schema_is_enforced(tmp_path):
    store = compose_store(tmp_path, shapes=False)
    g = store / "points" / "transcripts" / "zarr.json"
    meta = json.loads(g.read_text())
    del meta["attributes"]["points:n_points"]
    g.write_text(json.dumps(meta))
    r = validate_store.validate(store, sample=100)
    assert any("points:n_points" in m for m in fails(r))


# ── legacy_to_points ──────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def legacy_folder(tmp_path_factory) -> Path:
    """A legacy molecule folder from the tiny parquet, as the processor
    writes it (5,000 molecules, 3D, raw coordinates)."""
    import process_single_molecule as psm

    out = tmp_path_factory.mktemp("legacy") / "molecules"
    psm.process_single_molecule_data(
        input_file=str(TINY / "single-molecule" / "tiny.parquet"),
        output_folder=str(out),
        dataset_type="xenium",
    )
    return out


def test_legacy_folder_becomes_a_valid_points_element(tmp_path, legacy_folder):
    store = tmp_path / "store.zarr"
    out = legacy_to_points.convert(legacy_folder, store, levels=3, tile_size=200, log=lambda *_: None)
    assert out == store / "points" / "transcripts"
    r = validate_store.validate(store, sample=2000)
    assert not r.failed, fails(r)

    # Molecule for molecule: every gene's x, y, z in the folder are in the element
    manifest = legacy_to_points.read_manifest(legacy_folder)
    genes = manifest["genes"]["unique_gene_names"]
    _, attrs = pz.open_element(out)
    assert attrs["spatialdata_attrs"]["feature_key"] == "gene"
    assert attrs["points:n_points"] == manifest["statistics"]["total_molecules"]
    for i, gene in enumerate(genes[:3]):
        expected = legacy_to_points.read_gene(legacy_folder, gene)
        got = pz.feature_points(out, gene, level=0)
        assert len(got["x"]) == len(expected)
        for k, axis in enumerate(("x", "y", "z")):
            assert np.array_equal(np.sort(np.asarray(got[axis], dtype=np.float32)), np.sort(expected[:, k]))


def test_legacy_folder_can_drop_z(tmp_path, legacy_folder):
    store = tmp_path / "flat.zarr"
    out = legacy_to_points.convert(legacy_folder, store, levels=2, tile_size=200, drop_z=True, log=lambda *_: None)
    _, attrs = pz.open_element(out)
    assert attrs["axes"] == ["x", "y"]
    assert not validate_store.validate(store, sample=500).failed


# ── spatialdata round-trip ────────────────────────────────────────────────

spatialdata = pytest.importorskip("spatialdata", reason="the spatialdata package is the external reader")


@pytest.mark.xfail(
    reason="spatialdata 0.8 rejects the tiny fixture's element version ('0.5'): "
    "the fixture is written by anndata and our scripts, not by spatialdata. "
    "The other direction (a store spatialdata writes passes our validator) holds.",
    strict=True,
)
def test_spatialdata_reads_our_table_and_image(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=True)
    sdata = spatialdata.read_zarr(store, selection=("tables", "images", "shapes"))
    table = sdata.tables["table"]
    assert table.n_obs == 30
    assert "spatial" in table.obsm
    assert "dapi" in sdata.images
    assert "cell_boundaries" in sdata.shapes
    assert len(sdata.shapes["cell_boundaries"]) == 30


def test_a_store_spatialdata_writes_passes_the_validator(tmp_path):
    import anndata as ad
    import pandas as pd
    from spatialdata.models import Image2DModel, TableModel

    n = 12
    rng = np.random.default_rng(1)
    table = ad.AnnData(
        X=rng.poisson(1.0, size=(n, 5)).astype(np.float32),
        obs=pd.DataFrame({"cell_type": pd.Categorical(["a", "b"] * (n // 2)), "region": "img", "cell_id": np.arange(n)}, index=[f"c{i}" for i in range(n)]),
    )
    table.obsm["spatial"] = rng.uniform(0, 64, size=(n, 2)).astype(np.float32)
    table.var_names = [f"g{i}" for i in range(5)]
    table = TableModel.parse(table, region="img", region_key="region", instance_key="cell_id")
    image = Image2DModel.parse(rng.integers(0, 255, size=(1, 64, 64), dtype=np.uint8), dims=("c", "y", "x"), scale_factors=[2])
    sdata = spatialdata.SpatialData(images={"img": image}, tables={"table": table})
    store = tmp_path / "written.zarr"
    sdata.write(store)
    r = validate_store.validate(store)
    # What spatialdata leaves out is reported, never silently accepted
    assert not r.failed, fails(r)
