# SpatialData-DB root metadata

The root `zarr.json` of one artifact per technology in
`scverse/spatialdata-db`, fetched on 2026-10-18 from
`https://lamin.ai/storage/s3/scverse-spatial-eu-central-1//.lamindb/<uid>.zarr/zarr.json`,
as the raw material for `docs/plans/SDDB.md`. Re-fetch with curl (the
proxy answers 403 to Python's default user agent).

| File | Technology | Dataset |
|---|---|---|
| `jCLxdY1tWKmbfWOB.json` | 10x Atera | human breast (breast cancer) |
| `0W4MBvLWKZXauWhZ.json` | 10x Xenium | mouse bone |
| `fJSOO03VzCTYtD1o.json` | Visium | human cerebellum, targeted neuroscience panel |
| `haOzFbJ7Bc9x7NW6.json` | Visium HD | human ovary (ovarian serous adenocarcinoma) |
| `gZwZDQTAM8vvm8sE.json` | spatial transcriptomics | mouse spinal cord (ALS), NCBI137 |
| `YRIWXhfhZHe30Vf6.json` | CODEX (2026-10-21, from the bucket) | human pancreas: dense `X`, an 84-channel uint16 image chunked one channel at a time |
| `fF2M7WUbkFnMMhpC.json` | MALDI-MSI (2026-10-21) | human kidney: `var` index `formula_adduct`, the optical image placed by an affine rotated 91° |
| `571yTEm4hqlkXRPf.json` | CosMx, protein (2026-10-21) | human tonsil: 1,464,516 cells, no `obsm/spatial`, positions in `cells_polygons(_for_vitessce)` |
