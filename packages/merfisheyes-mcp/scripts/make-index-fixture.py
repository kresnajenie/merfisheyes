"""Cut the committed test index (test/fixtures/index/) from a real one.

    python scripts/make-index-fixture.py INDEX_DIR

INDEX_DIR holds the published index (.index/latest/ downloaded): catalog,
rollup, technology, gene index and gene lookup. The fixture keeps ~100 rows
chosen to exercise each search rule — ontology descendants, a re-release
pair, rows from Explore only, ages, sexes, a few genes — and remaps the gene
bitsets to its own row order, as the index build would.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

SRC = Path(sys.argv[1])
OUT = Path(__file__).resolve().parent.parent / "test" / "fixtures" / "index"
GENES = {"gad1", "rab3b", "slc17a7"}
COLUMNS = [  # src/library/data.ts CATALOG_COLUMNS
    "uid", "canonical_key", "study_id", "collection_name", "organism",
    "organism_id", "tissue", "tissue_id", "tissue_system", "tissue_system_id",
    "tissue_organ", "tissue_organ_id", "disease", "disease_id",
    "disease_class", "disease_group", "development_stage", "technology",
    "vendor", "assay", "panel_size", "validation_status", "tier", "license",
    "n_obs", "n_features", "n_transcripts", "published_at",
    "publication_date", "title", "sample_title", "sample_id", "donor_id",
    "doi", "pmid", "modality", "obs_unit", "timepoint", "timepoint_kind",
    "timepoint_days", "treatment", "genotype", "condition",
    "same_section_group", "same_section_relation", "same_section_release",
    "meyes_id", "meyes_source", "meyes_kind", "meyes_page_url",
    "meyes_viewer_url", "meyes_description", "meyes_sex", "meyes_n_sections",
    "meyes_dims", "meyes_format", "meyes_metadata_from", "meyes_bil_code",
    "meyes_investigator", "meyes_institution", "meyes_cell_type_levels",
]

cat = pq.read_table(SRC / "catalog.parquet", columns=COLUMNS)
df = cat.to_pandas()
pick = set()


def take(mask, n):
    pick.update(df.index[mask][:n])


take(df.meyes_metadata_from == "store", 20)
take(df.meyes_metadata_from == "explore", 6)
take(df.same_section_relation == "rerelease", 4)
take(df.tissue.str.contains("hippocamp|dentate", na=False) & (df.organism == "mouse"), 5)
take((df.organism == "human") & (df.disease_class == "cancer"), 10)
take((df.organism == "mouse") & (df.tissue_organ == "spinal cord"), 5)
take(df.timepoint_kind.eq("age") & df.timepoint_days.notna() & (df.organism == "mouse"), 8)
take(df.disease.eq("Alzheimer disease"), 4)
take(df.meyes_sex.eq("female"), 3)
# Re-release pairs whole, so the fold has two members.
groups = set(df.loc[list(pick)].same_section_group.dropna())
take(df.same_section_group.isin(groups) & (df.same_section_relation == "rerelease"), 50)
rows = sorted(pick)
sub = cat.take(pa.array(rows))
OUT.mkdir(parents=True, exist_ok=True)
pq.write_table(sub, OUT / "catalog.parquet", compression="snappy")

ids = sub.column("meyes_id").to_pylist()
(OUT / "gene_search_meta.json").write_text(json.dumps({
    "meyes_ids": ids,
    "order_hash": hashlib.sha256(",".join(ids).encode()).hexdigest(),
}))

lk = pq.read_table(SRC / "gene_lookup.parquet").to_pandas()
lk = lk[lk.symbol.isin(GENES)].copy()
n = len(df)


def remap(b):
    bits = np.unpackbits(np.frombuffer(b, np.uint8), bitorder="little")[:n]
    return np.packbits(bits[rows].astype(np.uint8), bitorder="little").tobytes()


lk["bitset"] = lk.bitset.map(remap)
lk["n_datasets"] = lk.bitset.map(
    lambda b: int(np.unpackbits(np.frombuffer(b, np.uint8)).sum())).astype("int32")
lk = lk.sort_values(["symbol", "entity_id"])
pq.write_table(
    pa.Table.from_pandas(lk, preserve_index=False), OUT / "gene_lookup.parquet",
    row_group_size=8, compression="zstd", write_statistics=["symbol"])

gi = pq.read_table(SRC / "gene_index.parquet").to_pandas()
gi = gi[gi.meyes_id.isin(ids) & gi.symbol.str.lower().isin(GENES)]
pq.write_table(pa.Table.from_pandas(gi, preserve_index=False), OUT / "gene_index.parquet")

rollup = json.loads((SRC / "ontology_rollup.json").read_text())
used = set(df.loc[rows, "tissue_id"].dropna()) | set(df.loc[rows, "disease_id"].dropna())
keep = set(used)
for t in used:
    keep.update(rollup.get(t, {}).get("ancestors_present", []))
(OUT / "ontology_rollup.json").write_text(
    json.dumps({k: v for k, v in rollup.items() if k in keep}, indent=1))
(OUT / "technology.json").write_text((SRC / "technology.json").read_text())
(OUT / "release_info.json").write_text(json.dumps({
    "release": "fixture", "built_at": "2026-11-18T00:00:00+00:00",
    "cut_from": json.loads((SRC / "release_info.json").read_text())["release"],
    "rows": len(rows)}, indent=1))
print(len(rows), "rows,", len(lk), "gene lookup rows,", len(gi), "gene index rows ->", OUT)
