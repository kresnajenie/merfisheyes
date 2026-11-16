# Inventory: what each MERFISHEYES row can fill

Status: measured 2026-11-16 for step 1 of [`docs/plans/METADATA.md`](../plans/METADATA.md);
the profile is [`metadata.md`](metadata.md), with the interim required set
(M10). Every public curated, BIL and community row of production Explore
(149), read through `GET /api/agent/v1/catalog/<id>`, and the 62 reprocessed
BIL stores' `stats.json`, `qc.json` and `catalog.json` on the workstation.
Lamin is not read: the stores are registered again when the fleet is
reprocessed (METADATA.md step 5). Produced by:

```bash
python scripts/metadata/inventory.py --rows <rows dir> --stores ~/data/bil-data/processed
```

The states (have / derive / fetch / person / missing) are defined in the script's
docstring; "filled" is have + derive. SpatialData-DB rows are not here: they are
their own rows, already in schema 3.

149 rows; profile: 16 required, 44 recommended columns counted (docs/specs/schemas/metadata.schema.json). 'with store' = once the row opens its reprocessed store (62 stores, 0 registered in Lamin).

| source | rows | required filled (median, min–max) | recommended filled (median, min–max) | rows with a store | required with store (median) | recommended with store (median) |
|---|---|---|---|---|---|---|
| bil | 92 | 13 (12–13) of 16 | 20 (9–20) of 44 | 59 | 16 | 31 |
| curated | 56 | 13 (11–16) of 16 | 20 (14–22) of 44 | 3 | 16 | 34 |
| community | 1 | 12 (12–12) of 16 | 11 (11–11) of 44 | 0 | — | — |

## By column (today)

have / derive / fetch / person / missing, per source

| column | level | bil | curated | community |
|---|---|---|---|---|
| `canonical_key` | required | 0 / 92 / 0 / 0 / 0 | 0 / 56 / 0 / 0 / 0 | 0 / 1 / 0 / 0 / 0 |
| `organism` | required | 0 / 92 / 0 / 0 / 0 | 0 / 56 / 0 / 0 / 0 | 0 / 1 / 0 / 0 / 0 |
| `organism_id` | required | 0 / 92 / 0 / 0 / 0 | 0 / 56 / 0 / 0 / 0 | 0 / 1 / 0 / 0 / 0 |
| `tissue` | required | 0 / 92 / 0 / 0 / 0 | 0 / 51 / 0 / 3 / 2 | 0 / 1 / 0 / 0 / 0 |
| `tissue_id` | required | 0 / 92 / 0 / 0 / 0 | 0 / 51 / 0 / 3 / 2 | 0 / 1 / 0 / 0 / 0 |
| `validation_status` | required | 0 / 92 / 0 / 0 / 0 | 0 / 51 / 0 / 0 / 5 | 0 / 1 / 0 / 0 / 0 |
| `n_obs` | required | 0 / 0 / 0 / 0 / 92 | 8 / 1 / 0 / 0 / 47 | 0 / 0 / 0 / 0 / 1 |
| `n_features` | required | 87 / 0 / 0 / 0 / 5 | 55 / 1 / 0 / 0 / 0 | 0 / 0 / 0 / 0 / 1 |
| `zarr_url` | required | 0 / 0 / 0 / 0 / 92 | 2 / 0 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `published_at` | required | 0 / 92 / 0 / 0 / 0 | 0 / 56 / 0 / 0 / 0 | 0 / 1 / 0 / 0 / 0 |
| `title` | required | 92 / 0 / 0 / 0 / 0 | 56 / 0 / 0 / 0 / 0 | 1 / 0 / 0 / 0 / 0 |
| `meyes_id` | required | 92 / 0 / 0 / 0 / 0 | 56 / 0 / 0 / 0 / 0 | 1 / 0 / 0 / 0 / 0 |
| `meyes_source` | required | 92 / 0 / 0 / 0 / 0 | 56 / 0 / 0 / 0 / 0 | 1 / 0 / 0 / 0 / 0 |
| `meyes_kind` | required | 0 / 92 / 0 / 0 / 0 | 0 / 56 / 0 / 0 / 0 | 0 / 1 / 0 / 0 / 0 |
| `meyes_page_url` | required | 92 / 0 / 0 / 0 / 0 | 56 / 0 / 0 / 0 / 0 | 1 / 0 / 0 / 0 / 0 |
| `meyes_viewer_url` | required | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `uid` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 0 / 0 / 0 / 56 | 0 / 0 / 0 / 0 / 1 |
| `study_id` | recommended | 0 / 92 / 0 / 0 / 0 | 0 / 56 / 0 / 0 / 0 | 0 / 1 / 0 / 0 / 0 |
| `disease` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 1 / 0 / 0 / 55 | 0 / 0 / 0 / 0 / 1 |
| `disease_id` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 1 / 0 / 0 / 55 | 0 / 0 / 0 / 0 / 1 |
| `development_stage` | recommended | 0 / 87 / 0 / 0 / 5 | 0 / 47 / 0 / 0 / 9 | 0 / 0 / 0 / 0 / 1 |
| `technology` | recommended | 0 / 0 / 0 / 92 / 0 | 0 / 55 / 0 / 1 / 0 | 0 / 0 / 0 / 1 / 0 |
| `vendor` | recommended | 0 / 0 / 0 / 92 / 0 | 0 / 55 / 0 / 1 / 0 | 0 / 0 / 0 / 1 / 0 |
| `assay` | recommended | 0 / 0 / 0 / 92 / 0 | 0 / 55 / 0 / 1 / 0 | 0 / 0 / 0 / 1 / 0 |
| `panel_size` | recommended | 87 / 0 / 0 / 0 / 5 | 55 / 0 / 0 / 0 / 1 | 0 / 0 / 0 / 0 / 1 |
| `license` | recommended | 87 / 0 / 0 / 0 / 5 | 47 / 0 / 0 / 0 / 9 | 0 / 0 / 0 / 0 / 1 |
| `total_counts` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `card_url` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 0 / 0 / 0 / 56 | 0 / 0 / 0 / 0 / 1 |
| `size_bytes` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 0 / 0 / 0 / 56 | 0 / 0 / 0 / 0 / 1 |
| `product` | recommended | 0 / 0 / 0 / 92 / 0 | 0 / 55 / 0 / 1 / 0 | 0 / 0 / 0 / 1 / 0 |
| `modality` | recommended | 0 / 92 / 0 / 0 / 0 | 0 / 56 / 0 / 0 / 0 | 0 / 1 / 0 / 0 / 0 |
| `modalities` | recommended | 0 / 92 / 0 / 0 / 0 | 0 / 56 / 0 / 0 / 0 | 0 / 1 / 0 / 0 / 0 |
| `biomaterial_type` | recommended | 0 / 92 / 0 / 0 / 0 | 0 / 56 / 0 / 0 / 0 | 0 / 1 / 0 / 0 / 0 |
| `hne_image` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `if_image` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `n_transcripts` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `publication_date` | recommended | 0 / 0 / 87 / 0 / 5 | 0 / 0 / 47 / 0 / 9 | 0 / 0 / 0 / 0 / 1 |
| `doi` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 1 / 0 / 0 / 0 |
| `publication_year` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 0 / 2 / 0 / 54 | 0 / 0 / 1 / 0 / 0 |
| `publication_verify_status` | recommended | 0 / 0 / 0 / 92 / 0 | 0 / 2 / 0 / 54 / 0 | 0 / 1 / 0 / 0 / 0 |
| `sample_id` | recommended | 0 / 87 / 0 / 5 / 0 | 0 / 47 / 0 / 9 / 0 | 0 / 0 / 0 / 1 / 0 |
| `default_table` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `dataset_url` | recommended | 92 / 0 / 0 / 0 / 0 | 54 / 0 / 0 / 0 / 2 | 1 / 0 / 0 / 0 / 0 |
| `timepoint` | recommended | 87 / 0 / 0 / 0 / 5 | 47 / 0 / 0 / 0 / 9 | 0 / 0 / 0 / 0 / 1 |
| `timepoint_kind` | recommended | 0 / 87 / 0 / 0 / 5 | 0 / 47 / 0 / 0 / 9 | 0 / 0 / 0 / 0 / 1 |
| `timepoint_days` | recommended | 0 / 87 / 0 / 0 / 5 | 0 / 47 / 0 / 0 / 9 | 0 / 0 / 0 / 0 / 1 |
| `genotype` | recommended | 0 / 68 / 0 / 0 / 24 | 0 / 0 / 0 / 0 / 56 | 0 / 0 / 0 / 0 / 1 |
| `obs_unit` | recommended | 0 / 90 / 0 / 0 / 2 | 0 / 9 / 0 / 0 / 47 | 0 / 1 / 0 / 0 / 0 |
| `observation_unit` | recommended | 0 / 90 / 0 / 0 / 2 | 0 / 9 / 0 / 0 / 47 | 0 / 1 / 0 / 0 / 0 |
| `observation_unit_source` | recommended | 0 / 90 / 0 / 0 / 2 | 0 / 9 / 0 / 0 / 47 | 0 / 1 / 0 / 0 / 0 |
| `segmentation_source` | recommended | 0 / 0 / 0 / 92 / 0 | 0 / 0 / 0 / 56 / 0 | 0 / 0 / 0 / 1 / 0 |
| `meyes_description` | recommended | 92 / 0 / 0 / 0 / 0 | 55 / 0 / 0 / 0 / 1 | 1 / 0 / 0 / 0 / 0 |
| `meyes_sex` | recommended | 0 / 84 / 0 / 0 / 8 | 0 / 47 / 0 / 0 / 9 | 0 / 0 / 0 / 0 / 1 |
| `meyes_sex_id` | recommended | 0 / 84 / 0 / 0 / 8 | 0 / 47 / 0 / 0 / 9 | 0 / 0 / 0 / 0 / 1 |
| `meyes_development_stage_id` | recommended | 0 / 87 / 0 / 0 / 5 | 0 / 47 / 0 / 0 / 9 | 0 / 0 / 0 / 0 / 1 |
| `meyes_n_sections` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `meyes_sections` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `meyes_dims` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `meyes_de_columns` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 2 / 0 / 0 / 54 | 0 / 0 / 0 / 0 / 1 |
| `meyes_pipeline` | recommended | 0 / 0 / 0 / 0 / 92 | 0 / 0 / 0 / 0 / 56 | 0 / 0 / 0 / 0 / 1 |

## By row

| source | BIL code | title | required | recommended | with store (req / rec) | fetch | person |
|---|---|---|---|---|---|---|---|
| bil |  | MERFISH Imaging of mesocorticolimbic regions in the mouse brain | 12 | 9 | — | 0 | 7 |
| bil |  | MERSCOPE Imaging of WT, Trem2-R47H, 5xFAD and Trem2-R47H; 5xFAD sample | 12 | 9 | — | 0 | 7 |
| bil |  | Spatial Transcriptomic Atlas of Adult Human Basal Ganglia: H22.30.001. | 12 | 9 | — | 0 | 7 |
| bil |  | Subcortical Spatial Transcriptomic Atlas of the Adult Rhesus Macaque | 12 | 9 | — | 0 | 7 |
| bil |  | VA00265 300-gene MERFISH results of human visual cortex (GW24 & GW34) | 12 | 9 | — | 0 | 7 |
| bil | ace-can-ill | VA00265 300-gene MERFISH results of a GW20 human prefrontal cortex | 13 | 17 | — | 1 | 6 |
| bil | ace-can-imp | VA00265 300-gene MERFISH results of a GW34 human prefrontal cortex | 13 | 19 | 16 / 29 | 1 | 6 |
| bil | ace-can-ink | VA00265 300-gene MERFISH results of a GW24 human prefrontal cortex | 13 | 17 | — | 1 | 6 |
| bil | ace-den-fix | Mouse whole-brain transcriptomic cell type atlas - MERSCOPE v1: 638850 | 13 | 20 | — | 1 | 6 |
| bil | ace-den-fog | Mouse whole-brain transcriptomic cell type atlas - MERSCOPE v1: 609889 | 13 | 20 | — | 1 | 6 |
| bil | ace-den-hen | P20 Mouse whole-brain transcriptomic cell type atlas - MERSCOPE v1 o - | 13 | 20 | — | 1 | 6 |
| bil | ace-dim-nut | Mouse whole-brain transcriptomic cell type atlas - MERSCOPE v1 | 13 | 20 | — | 1 | 6 |
| bil | ace-dip-tub | MERFISH Imaging of mesocorticolimbic regions in the mouse brain (Chron | 13 | 20 | — | 1 | 6 |
| bil | ace-dip-tug | MERFISH Imaging of mesocorticolimbic regions in the mouse brain (Chron | 13 | 20 | — | 1 | 6 |
| bil | ace-dip-use | MERFISH Imaging of mesocorticolimbic regions in the mouse brain (infec | 13 | 20 | — | 1 | 6 |
| bil | ace-dip-van | MERFISH Imaging of mesocorticolimbic regions in the mouse brain (infec | 13 | 20 | 16 / 35 | 1 | 6 |
| bil | ace-dip-vat | MERFISH Imaging of mesocorticolimbic regions in the mouse brain (infec | 13 | 20 | — | 1 | 6 |
| bil | ace-dry-dip | VA00265 300-gene MERFISH results of an infancy stage human prefrontal  | 13 | 19 | — | 1 | 6 |
| bil | ace-dry-dog | VA00265 300-gene MERFISH results of a GW24 human visual cortex | 13 | 17 | — | 1 | 6 |
| bil | ace-dry-dry | VA00265 300-gene MERFISH results of a GW34 human visual cortex | 13 | 19 | — | 1 | 6 |
| bil | ace-dry-dub | VA00265 300-gene MERFISH results of an infancy stage human visual cort | 13 | 19 | 16 / 34 | 1 | 6 |
| bil | ace-dud-vex | Subcortical Spatial Transcriptomic Atlas of the Adult Rhesus Macaque:  | 13 | 19 | 16 / 34 | 1 | 6 |
| bil | ace-dud-vow | Subcortical Spatial Transcriptomic Atlas of the Adult Rhesus Macaque:  | 13 | 19 | 16 / 34 | 1 | 6 |
| bil | ace-dud-wag | Subcortical Spatial Transcriptomic Atlas of the Adult Rhesus Macaque:  | 13 | 19 | — | 1 | 6 |
| bil | ace-dud-war | Subcortical Spatial Transcriptomic Atlas of the Adult Rhesus Macaque:  | 13 | 19 | — | 1 | 6 |
| bil | ace-dud-was | Subcortical Spatial Transcriptomic Atlas of the Adult Rhesus Macaque:  | 13 | 19 | — | 1 | 6 |
| bil | ace-dug-fin | Mouse whole-brain transcriptomic cell type atlas - MERSCOPE v1 | 13 | 20 | — | 1 | 6 |
| bil | ace-ear-dud | Mouse whole-brain transcriptomic cell type atlas - MERSCOPE v1 | 13 | 20 | — | 1 | 6 |
| bil | ace-irk-sag | Spatial Transcriptomic Atlas of Adult Human Basal Ganglia: H22.30.001. | 13 | 19 | — | 1 | 6 |
| bil | ace-irk-sap | Spatial Transcriptomic Atlas of Adult Human Basal Ganglia: H22.30.001. | 13 | 19 | — | 1 | 6 |
| bil | ace-irk-sat | Spatial Transcriptomic Atlas of Adult Human Basal Ganglia: H22.30.001. | 13 | 19 | — | 1 | 6 |
| bil | ace-irk-saw | Spatial Transcriptomic Atlas of Adult Human Basal Ganglia: H22.30.001. | 13 | 19 | — | 1 | 6 |
| bil | ace-irk-say | Spatial Transcriptomic Atlas of Adult Human Basal Ganglia: H22.30.001. | 13 | 19 | — | 1 | 6 |
| bil | ace-irk-set | Spatial Transcriptomic Atlas of Adult Human Basal Ganglia: H22.30.001. | 13 | 19 | — | 1 | 6 |
| bil | ace-irk-shy | Spatial Transcriptomic Atlas of Adult Human Basal Ganglia: H22.30.001. | 13 | 19 | 16 / 34 | 1 | 6 |
| bil | ace-let-law | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | — | 1 | 6 |
| bil | ace-let-lax | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | — | 1 | 6 |
| bil | ace-let-lay | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | — | 1 | 6 |
| bil | ace-lot-all | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-lot-web | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-lot-wet | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-lot-who | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-lot-wig | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-lot-win | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-lot-wit | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-lot-yes | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-lot-yet | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-lot-zoo | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-ace | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-act | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-add | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-age | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-aim | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-air | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-and | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-ant | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-ape | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-arm | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-art | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-ash | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-ask | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bad | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bag | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-ban | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bar | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bat | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bay | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bed | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 17 | 16 / 31 | 1 | 6 |
| bil | ace-low-bet | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bid | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-big | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bin | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bit | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bog | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-boo | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-box | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bud | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bug | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bun | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-bus | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-cab | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-can | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-cap | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-car | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-cat | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-cop | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-cot | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-cow | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-cry | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-cub | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | 16 / 31 | 1 | 6 |
| bil | ace-low-zip | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 17 | 16 / 31 | 1 | 6 |
| bil | ace-met-war | Developing mouse whole-brain transcriptomic cell type atlas - MERSCOPE | 13 | 20 | — | 1 | 6 |
| community |  | Re-activation of neurogenic niches | 12 | 11 | — | 1 | 6 |
| curated |  | Allen Brain Mouse Whole-Brain Transcriptomic Cell Type Atlas | 11 | 14 | — | 0 | 3 |
| curated |  | Drosophila Embryo and Larva in 3D (Flysta3D) | 11 | 14 | — | 0 | 5 |
| curated |  | Human Breast Cancer, Cells and Molecules (Xenium) | 14 | 16 | — | 0 | 3 |
| curated |  | Human Embryo Atlas, CS12–CS23 (HESTA) | 14 | 14 | — | 0 | 3 |
| curated |  | Human Skin, Cells, Polygons and Molecules (Xenium) | 16 | 22 | — | 0 | 3 |
| curated |  | Mouse Brain Saggital 10 Slices | 13 | 16 | — | 0 | 7 |
| curated |  | Mouse Organogenesis Sections (MOSTA) | 14 | 14 | — | 0 | 3 |
| curated |  | Mouse Organogenesis in 3D (MOSTA-3D) | 11 | 16 | — | 1 | 4 |
| curated |  | Planarian Whole-Body Regeneration in 3D (PRISTA4D) | 11 | 16 | — | 1 | 4 |
| curated | ace-owl-cry | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-cub | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-cup | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | 16 / 34 | 1 | 2 |
| curated | ace-owl-cut | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-day | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-den | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-did | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-die | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-dig | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-dim | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-dip | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-dog | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-dry | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-dub | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-dud | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-dug | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-ear | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-eat | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | 16 / 34 | 1 | 2 |
| curated | ace-owl-eel | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-egg | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-elf | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-elk | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-elm | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-end | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fan | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-far | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fat | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fed | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-few | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | 16 / 34 | 1 | 2 |
| curated | ace-owl-fib | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fig | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fin | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fit | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fix | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fly | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fog | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-foo | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fox | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fry | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-fun | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-gab | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-gag | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-gap | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-gas | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-gel | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-gem | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
| curated | ace-owl-get | Spatial transcriptomic atlas of the marmoset basal ganglia: CJ23.56.00 | 13 | 20 | — | 1 | 2 |
