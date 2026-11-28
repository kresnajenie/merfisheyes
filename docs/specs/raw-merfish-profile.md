# Raw MERFISH run profile

What a raw MERFISH run folder must hold for `scripts/raw_merfish.py` to
check it (`meyes_raw_check`) and convert it into one SpatialData store
([RAW-MERFISH.md](../plans/RAW-MERFISH.md)). The baseline is the lab
pipeline as Jess's runs leave it (decode → stitch per set → drop overlap
duplicates → count per cell): "what the data needs to look like before
conversion". Runs processed another way are refused with what they lack
until an adapter brings them to this shape.

The names below (folders, members, columns) are defined once, at the top of
`scripts/raw_merfish.py`.

## Members

| Requirement | Level |
|---|---|
| `final_spots/` holding one `.npz` per **set** (the set's id is the number after `set` in the file name, else the file's stem) | MUST |
| Each spot file's `XF`: float32 or float64, N × ≥ 5, columns `z` (plane), `x`, `y` (stitched pixels), `score`, field (an index into `fovs`); rows grouped by field | MUST |
| Each spot file's `genes`: N strings, one per molecule — panel genes and `blank…` codewords | MUST |
| Each spot file's `pix_size` (µm per z plane, x pixel, y pixel — x and y equal) and `fov_size` (planes, pixels, pixels) | MUST |
| Each spot file's `fovs` (field names) and `pos_fovs_pix` (each field's x, y corner in the set's stitched pixels; a field spans `pos + [0, fov_size)`) | MUST |
| `.npz` members stored, not deflated, so they are read in place and in pieces | SHOULD (a deflated member is read whole: slower, more memory) |
| `final_matsWithGFP/` or `final_mats/` holding `<field>--mat.npz` for every field of every set, matched by name to `fovs` | MUST for cells |
| Each mat's `Mcells` (cell names, unique across the run), `Mgenes` (the panel, blanks included, the same in every mat), `Mcounts` (cells × genes, molecules per cell), `Xum` (cells × 3: z, x, y in µm, the set's frame) | MUST for cells |
| A mat's `ctsGFP` (GFP spots per cell) | MAY — read as `gfp_count` |
| `POS/pos_<K>.txt`: the stage positions imaged, one line a field, a field's number its line's | MAY — the check flags fields imaged but never decoded (a decode stopped early leaves holes) |
| `GFP_xfit/XGFPF_set<K>.npz`, `Segmentation2d/`, mosaics, an h5ad of the run | MAY — later steps; ignored by the check |
| Pickles (`*.pkl`) | never opened: loading one runs code |

## Rules the pipeline has already applied

- **Score**: molecules are kept at `score > -0.5` (`allth3`). No threshold is
  applied again; `score` is kept on each molecule.
- **Overlaps**: neighbouring fields overlap (~14 µm on a 3,000-pixel field).
  A molecule is kept only in the field whose centre is nearest, by Chebyshev
  distance — `KDTree(centres).query(xy, p=inf)` keeps a molecule when the
  nearest centre is its own field's. The **spots are already
  de-duplicated**; the **mats are not**: the same rule is applied to each
  cell's `Xum` x, y against the field centres (`pos_fovs_pix + fov_size/2`,
  in µm) when cells are read.
- **Fields imaged twice**: the lab's runs image a set's first fields (00,
  01) again as its last two, a few pixels apart. A cell goes to the nearer
  copy's centre; an exact tie (Chebyshev distances tie whenever x
  dominates) goes to the lower index. The lab's own runs split such ties by
  sub-pixel stage positions that the saved float32 `pos_fovs_pix` no longer
  hold, so a run's cell count can differ from the lab's h5ad by a cell or
  two: 313,075 vs 313,076 (`20260326_WT_BRBB_Mer`), 296,057 vs 296,058
  (`MouseA56`). Any other pair of fields within 5 % of a field's width is
  flagged.

## Units and frames

- µm = `XF[:, :3] × pix_size` (z, x, y); `Xum` is already in µm, in the
  same frame.
- Each set is its own frame, its origin at the set's lowest field corner.
  Sets do not share coordinates; in a store each set is a `section`,
  shifted onto a grid.
- A set's frame has y going down (stitched image pixels); the viewer draws
  +y upwards. A store flips y so each section reads as in the lab's own
  plots, and lays the sections on a grid read left to right, top to
  bottom: x = x_set + dx, y = dy − y_set, each set's dx, dy in the table's
  `uns['merfisheyes_sections']` (back: x − dx, dy − y). No rotation: the
  fields' own orientation is right as written (cells imaged twice in an
  overlap sit within a micron of their twin only as written).

## Identity and blanks

- A field is `fovs[i]` of its set; a cell is its `Mcells` name.
- Molecules carry no cell id.
- A codeword whose name starts with `blank` is a blank: counted for QC (the
  misidentification estimate is the blank fraction), never drawn as a gene.
