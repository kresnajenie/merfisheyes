# Questions on the Sep2026 spiralia dataset

Context: porting `Final_analyzed_objects/*_RNA_domain_annotated` + the cell
segmentation into MERFISHeyes as a 3D single-molecule viewer. Starting with
**MER6-2_E3_1**. Everything below is either a real ambiguity in the files or a
fact I derived and would like confirmed before it gets baked in.

---

## A. Cell segmentation — which one is canonical?

**A1.** There are three voxel masks per embryo, all with an identical label set
(`1,3,4…24, 5000–5006` for E3_1):

| file | shape (E3_1) | voxel size |
|---|---|---|
| `Segmentation/Bogdan/Old/MER6-2_E3_1_segm.npz['segm']` | 370 × 725 × 725 | 1.0 / 0.708 / 0.708 µm |
| `Segmentation/Bogdan/AutoRefined_Masks/..._cell_mask_DS4.npz` | 93 × 182 × 182 | 4.0 / 2.832 / 2.832 µm |
| `Segmentation/Bogdan/ManualCurate_Masks/..._curated_cell_mask_DS4.npz` | 93 × 182 × 182 | 4.0 / 2.832 / 2.832 µm |

Which is the one to render? What is the difference between AutoRefined and
ManualCurate, and is the `Old/*_segm.npz` one superseded by them or just a
higher-resolution version of the same thing?

**A2.** There are also two prebuilt mesh sets, 26 cells each for E3_1:

- `annotated_domain_viz_cell_meshes/` — marching cubes on the manually curated
  mask (`mask_source: manual_curated`).
- `annotated_domain_viz_pointcloud_nonoverlap_cell_meshes/` — built from the
  molecule point cloud on a shared grid, guaranteed non-overlapping. This one
  additionally carries a `cell_identity` attribute per mesh (`3B`, `1a1`, `pb`…).

Which should the viewer show? The pointcloud one is the only one with cell
identities attached, but the mask-derived one is presumably closer to the actual
segmentation.

**A3.** Each mesh group has both `vertices_xyz` and `vertices_reoriented_xyz`.
**What is the reorientation?** If it is a per-embryo alignment (animal–vegetal
axis? a common reference frame?), then using the reoriented vertices requires
applying the same transform to the molecule coordinates, and that transform is
not stored in the h5. Is it recoverable, or should we just use `vertices_xyz`?

**A4.** Only 31 of the 45 embryos have mesh files. Is the remainder still coming,
or are those 31 the ones considered good?

---

## B. Coordinates and units

**B1.** I derived that `FISH['Xh']` is in **pixels** with voxel size
**z = 0.5 µm, x = y = 0.177 µm** — by dividing `domain_meta.csv`'s
`center_z/x/y_um` by the mean `Xh` position of each domain. The ratio is exact
to 1e-10 across all 18 domains. Please confirm, and confirm whether it is the
same for all 45 embryos / all five experiments.

**B2.** The notebook applies `layer.scale = (1.2, -1, 1)` in napari. The `-1`
flips one axis. Is that a real handedness correction (i.e. the data is
left-handed and needs flipping for correct anatomy), or just a viewing
convenience? And where does `1.2` come from — the true z:xy aspect after the
2/4/4 downsample is 1.41, not 1.2.

---

## C. Molecule filtering

**C1.** The panel has 306 entries. I find these break down as:

- **213 real genes**
- **80 blank control probes** (`blank0000::all::MERFISH` …)
- **13 duplicate/alternate probe entries** (`::alt`, `::dup-0`, `::dup-1`,
  or a second `tar_region` for the same node) — these are exactly the entries
  with `gns_names_keep == False`, and **no** molecule survives `molecules_keep`
  on any of them, so that flag is redundant once `molecules_keep` is applied.

Correct?

**C2.** I found that **`analysis_mask` = `molecules_keep` AND (not a blank probe)**
— verified exactly, and it accounts for all 8,108 molecules of difference. It
also means the `RNA_domain_anno == "unassigned"` category is *precisely the blank
control probes*, not unassigned tissue. Confirm? And which mask do you consider
the right one for display — I'm currently using `molecules_keep` and letting the
viewer drop the blanks, which lands on `analysis_mask`.

**C3.** 14,430 molecules under `molecules_keep` have `cell_id == 0`. What does 0
mean here — outside the embryo, in an inter-cellular gap, or unsegmented? We
plan to show them as their own "unassigned" category rather than dropping them.

**C4.** `nuc_ids_eroded` is non-zero for only 62,385 of 3.82M molecules (1.6%),
and never disagrees with `cell_ids` where it is non-zero. Is this the intended
"molecule is in the nucleus" flag, and what erosion was applied? (Asking because
1.6% nuclear seems low, and we may want a nuclear/cytoplasmic toggle.)

---

## D. Naming

**D1.** For the viewer we're displaying genes as `SYMBOL (full probe string)`,
e.g. `ACTC_2 (NODE168517::ACTC_2|P12716::all::MERFISH)`. Is `ACTC_2` the right
symbol to pull out, or is there a preferred short name per node?

**D2.** `NODE29150::CTNB|P35224` and `NODE23626::CTNB_2|P35224` are separate
panel entries mapping to the same UniProt (P35224). Different paralogs, or two
probe sets against the same gene? Should they be merged in the UI?

**D3.** Cell ids for E3_1 are `1,3,4…24` and then `5000–5006`. Why the jump — is
5000+ meaningful (a later curation batch? a different segmentation pass?) or an
artifact we should ignore?

**D4.** Two cells in E3_1 are both named `pb` (ids 1 and 3). Polar bodies, I
assume — should they be shown, and should they be distinguishable?

---

## E. Scope

**E1.** Is the 306-probe panel identical across MER1 / MER2 / MER5-1 / MER6-1 /
MER6-2, or does it differ per experiment?

**E2.** `domain_meta.csv` has `stage` (`16-cells`, etc.), `annotation`,
`annotation_cluster`, `annotation_cluster_fine`, `annotation-macro`. Which of
these is the one to surface in the UI as the compartment label? We're currently
using the per-molecule `RNA_domain_anno`.
