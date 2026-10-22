# Table profile

What a SpatialData table must be for `lib/spatialdata/SpatialDataTableAdapter.ts`
to read it as cells: coordinates, cluster columns, genes, and one gene's
expression on demand. A table is an AnnData group in Zarr v3
(`anndata` ≥ 0.11 writes it; `tables/<name>` in the store).

| Requirement | Level |
|---|---|
| The store's root has `spatialdata_attrs` and exactly one `tables/<name>` is read: the curated `metadata_general["Default Table"]` when the root names one, else `tables/table`, else the first by name (`defaultTableName` in `lib/spatialdata/store.ts`) | MUST |
| `obsm/spatial`: a 2-D numeric array (float, or int64 as Visium writes it) of `n_obs × 2` or `× 3`, the cells' coordinates in their intrinsic system (microns, or full-image pixels for Visium). Without it, the region shapes element the table annotates gives each cell its position: a Point element (SpatialData-DB's `<region>_for_vitessce`, with the cell's `radius` and id) is read first, else the polygons' centres; cells are matched by `instance_key`, and the element must be in the cells' frame (a CosMx table has no `obsm`) | MUST, one or the other |
| `obs/_index` (or the `_index` named in `obs` attrs): cell ids, as a string array or a `nullable-string` encoding | MUST |
| `var`'s index column (the `_index` named in `var` attrs, `_index` by default; `formula_adduct` in a MALDI table, `gene` in a CosMx one): feature names, likewise. The viewer calls every feature a gene, whether it is a gene, an antibody, a peak or an ion | MUST for expression |
| `X`: `csc_matrix` encoding (`indptr`, `indices`, `data`), so one gene is a contiguous slice; or `layers/X_csc` holding that copy when `X` is CSR; or a dense `array` (a CODEX or PhenoCycler table: every image channel a feature, `marker_role` telling markers from blanks), where one feature reads every row chunk | MUST for gene expression; a CSR `X` with no CSC copy reads the whole matrix per gene |
| `X`'s `indices` and `data` chunked small enough that one gene reads a few chunks (the ten-section store: sharded, 17 reads a gene) | SHOULD |
| Root attributes `instance_key` and `region_key` naming the obs columns that join the table to a shapes element and name the region; these columns are not offered as cluster columns | MUST when a shapes element is to be matched to the cells |
| `obs/<column>`: `categorical` (`codes` + `categories`), `string` / `nullable-string`, or numeric; columns with ≤ 100 unique values are categorical in the viewer, others numerical | SHOULD |
| `annotations/annotations.json` beside the store (the folder, not the store): `{ columns: { <column>: { categories: [...] } }, de_stats: [...] }`, with a file per annotated column, for cell type names written after the table | MAY |
| Everything else in the group (`uns`, `obsp`, `varm`, `raw`) | ignored |

Reading (`docs/specs/LOADING.md`): the root `zarr.json`'s consolidated
metadata answers every node's `zarr.json`; coordinates and one obs column
are read at open; other columns on demand; a gene through `indptr` once and
the overlapping chunks of `indices` and `data`; the shapes' instance ids
through `readInstanceIds()` from the `instance_key` column.

Known stores against this profile: the ten-section store and set1 pass
(CSC, sharded); atera passes (CSC); `skin_adapt_v2` fails — its root
attributes name a `layers/X_csc` that does not exist, so every gene scans
the CSR matrix. The fix is in the writer: `scripts/validate_store.py`
reports it.
