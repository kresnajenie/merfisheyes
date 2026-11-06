# Plan: MERFISHEYES datasets in a Lamin catalog

Record every dataset MERFISHEYES serves in the Lamin instance
**`merfisheyes/eyes-catalog`**: what it is, where its bytes are, who added
it, which pipeline made it. S3 keeps the bytes and the viewer keeps reading
them by URL; Lamin is the record of them. This is the reverse of
[SDDB.md](SDDB.md), where we *read* someone else's Lamin instance.

Keep this file current: a step that changes a fact or a decision edits it
in the same commit.

## 1. Facts (2026-11-05)

- **Instance** `merfisheyes/eyes-catalog`, Postgres, lamindb 2.12.2,
  module `bionty`. Readable anonymously; writing needs `lamin login` with an
  account that is a member.
- **Storage locations** registered with it:

  | root | what |
  |---|---|
  | `s3://lamin-us-west-2/RebS9q9gmXSZ` | the default; Lamin's own bucket, private, reached through Lamin |
  | `s3://merfisheyes-us-west-2?endpoint_url=…r2.cloudflarestorage.com` | Cloudflare R2, ours, unused so far |
  | `s3://merfisheyes-bil/bil-psc-zarr` | new SpatialData stores, registered 2026-11-05 |
  | `s3://merfisheyes-bil/bil-psc-data2` | legacy chunked folders (step 2) |

- **Registration is in place.** `aws s3 sync` puts the bytes where the
  viewer reads them; `ln.Artifact("s3://…").save()` records that path. The
  key is the path below the storage root and is *not* virtual, so the URL
  the viewer opens and the path the record names are the same files.
  Uploading a local folder through Lamin instead would write it to
  `<root>/.lamindb/<uid>…` — a path the viewer does not know.
- **Never `Artifact.from_spatialdata()` and never `.load()`** on our
  stores: both go through `spatialdata`, which (0.8) cannot read our points
  element, and the first rewrites the store.
- Registering a prefix writes one marker, `<root>/.lamindb/storage_uid.txt`.
- **Local env** `meyes-lamin` (conda, a clone of `meyes-prepare` plus
  lamindb). Not the pinned `meyes-prepare` itself: lamindb moves `anndata`
  0.13.4 → 0.13.2, off the worker's pins.
- **AWS** from a workstation: IAM user `merfisheyes-data-sync`, scoped to
  the `merfisheyes-bil` prefixes above. Never a root key.

## 2. The tool

[`scripts/bil-scripts/register_lamin.py`](../../scripts/bil-scripts/register_lamin.py)
registers samples and **checks each artifact**, exiting non-zero if any check
fails:

| check | passes when |
|---|---|
| `record` | the Lamin path is the S3 path and the key is not virtual |
| `listing` | a fresh S3 listing has the record's file count and bytes |
| `public` | the entry file (`zarr.json`, `manifest.json`, `manifest.json.gz`) answers on the public URL with CORS for www.merfisheyes.com; a store's `zarr.json` has consolidated metadata |
| `lamin` | the entry file reads back through `Artifact.get(key=…).path` |

Each artifact gets features `bil_code`, `data_format`, `dataset_kind`,
`species`, `tissue`, `platform`, `investigator`, `n_cells`, `n_genes`,
`n_molecules`, `pipeline_commit`, `viewer_url` (+ the QC ones below on a
store), and the project `BIL`. Re-running is safe: same bytes return the same
record, values are replaced.

**One row per dataset in the record sheet "BIL datasets"** (`ln.Record`
of that type, named by BIL code; checks `row` and `catalog`). Lamin fills
nothing in by itself — not artifacts' features, not records: every value
comes from a file the pipeline wrote on that run.

| from | file (beside the store, synced with it) | written by | fields |
|---|---|---|---|
| the store | `stats.json` | `prepare.py`, `slices_to_points.py` | `n_cells`, `n_genes`, `n_molecules`, `pipeline_commit` |
| QC | `qc.json` | `qc_store.sbatch` (`collect_dataset_stats.py --json`) | `n_sections`, `n_genes_detected`, `median_counts_per_cell`, `mask_percentile`, `pct_artifact`, `mmc_levels` |
| catalog | `catalog.json` | `catalog_store.sbatch` (`catalog_entry.py`, BIL API) | `title`, `species`, `tissue`, `platform`, `investigator`, `institution`, `publication_year`, `license`, `bil_url` |
| this script | — | `register_lamin.py` | `data_format`, `viewer_url`, `store` (the artifact the viewer opens), `artifacts` (all of them) |

A sheet accepts only the features its schema lists, and a schema is found by
its content: adding a column to `SHEET_COLUMNS` makes a new schema and the
sheet is moved onto it, rows keep their values. `--catalog` is the fallback
for a store with no `catalog.json`, and the only source for legacy folders.

`build_bil_catalog.py` turns the `catalog.json` files back into
`bil-completed.json` for the website, so Explore and Lamin are made from the
same file.

```bash
python register_lamin.py ace-dip-vat --catalog bil-completed.json          # SpatialData
python register_lamin.py --format legacy --catalog bil-completed.json ace-can-ill …
python register_lamin.py --check-only ace-dip-vat                           # verify only
```

On the cluster, `launch_spatialdata_pipeline.sh --sync --lamin` runs
`register_lamin.sbatch` after the sync.

## 3. Steps

1. **One SpatialData store, end to end** — done 2026-11-05. `ace-dip-vat`
   (474,986 cells, 646 M molecules) synced to
   `s3://merfisheyes-bil/bil-psc-zarr/ace-dip-vat/` (3,034 files,
   9,383,391,282 B, identical to the local store) and registered; all four
   checks pass. Store artifact `RIeFAgjpb8WpQqpW0000`.
   Backfilled 2026-11-05 with `qc.json`, `catalog.json` and the first sheet
   row.
2. **10–20 datasets** — done 2026-11-05. The BIL servers being down, 13 samples rebuilt from
   their legacy `bil-psc-data2` output (`reconstruct_inputs.py`:
   raw-shaped inputs, the old MapMyCells labels, masked cells gone so
   `--mask-percentile 0`) run through the new pipeline locally, synced to
   `bil-psc-zarr` and registered — with `ace-dip-vat`, 14 datasets. The
   rebuilt `ace-dip-vat` was not synced over the store made from its raw
   data. Gate: every artifact and sheet row passes every check. Registering
   the legacy folders in place stays possible (`--format legacy`) and needs
   read access to `bil-psc-data2`. Result: 70 artifacts (5 per dataset) and
   14 sheet rows, every check passing; ace-irk-shy is cells only (its legacy
   output had no molecules) and the three ace-owl-* Xenium sets have no
   cell-type labels.
3. **The cluster env** — `~/merfisheyes_lamin_env` with lamindb on PSC,
   `lamin login`, then the first `--sync --lamin` run.
4. **Open questions** — whether the app's `Dataset` table links to Lamin
   records or Explore queries Lamin directly; species/tissue as `bionty`
   ontology terms instead of free text; whether R2 replaces S3 for serving;
   whether sessions are tracked as Lamin runs (`ln.track()`).
