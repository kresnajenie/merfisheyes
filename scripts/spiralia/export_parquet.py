"""Export one spiralia embryo to a MERFISHeyes-ready single-molecule parquet.

Reads the dill `*_RNA_domain_annotated` object and writes Xenium-schema columns
(so the uploader's sniffer auto-detects the mapping) plus the extra per-molecule
annotations this dataset carries.

  python export_parquet.py [EMBRYO_ID] [OUT_DIR]
"""
import re
import sys

import dill
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = "/home/data/yiqun-spiralia/Sep2026"
EMB = sys.argv[1] if len(sys.argv) > 1 else "MER6-2_E3_1"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else f"{ROOT}/merfisheyes_export"

# Verified against domain_meta.csv µm columns across every domain (ratio is
# exact to 1e-10): Xh is in pixels, z 0.5 µm, x/y 0.177 µm.
UM_PER_PX_Z = 0.5
UM_PER_PX_XY = 0.177

# Mirrors shouldFilterGene() in lib/utils/gene-filters.ts — a gene name hitting
# any of these would be silently dropped by the viewer.
VIEWER_GENE_FILTERS = [
    r"negative[\s_-]*control", r"neg[\s_-]*ctrl", r"unassigned",
    r"deprecated", r"codeword", r"blank", r"negcontrol",
]


def sanitize(name):
    """Mirrors the gene -> S3 filename rule in SingleMoleculeDataset.ts."""
    return re.sub(r"_+", "_", re.sub(r"[^a-zA-Z0-9]", "_", name)).strip("_")


def gene_label(anno, full):
    """Display name: '<symbol> (<full probe string>)'.

    anno is 'NODE168517::ACTC_2|P12716' for a real gene and plain 'blank0000'
    for a control probe, so the symbol is only split out when it's there.
    """
    symbol = anno.split("::")[1].split("|")[0] if "::" in anno else anno
    return f"{symbol} ({full})"


print(f"loading {EMB} ...", flush=True)
F = dill.load(open(f"{ROOT}/Final_analyzed_objects/{EMB}/{EMB}_RNA_domain_annotated", "rb"))

# Molecule filter: molecules_keep only. gns_names_keep (which drops 13
# duplicate/alternate panel entries) is redundant — molecules_keep already
# excludes every molecule on those entries.
icodes = F["icodesN"]
sel = F["molecules_keep"]
print(f"  {len(sel)} molecules -> {sel.sum()} after molecules_keep")

sc = F["gns_names_sc"]
labels = np.array([gene_label(a, s) for a, s in zip(F["gns_names_anno"], sc)], dtype=object)
present = np.unique(icodes[sel])

# Sanity checks against the viewer's own rules.
clashes = pd.Series([sanitize(s) for s in labels[present]]).duplicated().sum()
assert clashes == 0, f"{clashes} gene names collide after filename sanitization"

# The panel carries 80 blank control probes. The viewer's shouldFilterGene()
# drops anything matching /blank/, so they never reach the scene — which lands
# the rendered set exactly on the object's own `analysis_mask` (verified).
auto_filtered = [s for s in labels[present]
                 if any(re.search(p, s.lower()) for p in VIEWER_GENE_FILTERS)]
n_filtered_mols = int(np.isin(icodes[sel], [i for i in present if any(
    re.search(p, labels[i].lower()) for p in VIEWER_GENE_FILTERS)]).sum())
print(f"  panel entries present: {len(present)} "
      f"({len(present) - len(auto_filtered)} real + {len(auto_filtered)} blank controls)")
print(f"  the viewer will auto-filter the blanks: {n_filtered_mols:,} molecules "
      f"({100 * n_filtered_mols / sel.sum():.2f}%)")

Xh = F["Xh"]
cell_ids = F["cell_ids"][sel].astype(np.int32)

# cell_names is a Series indexed by cell id; 0 = not in any cell.
names = F["cell_names"]
name_map = {int(i): str(v).split("::")[1] for i, v in names.items()}

# cell_ids can carry labels that have no cell_names entry — segmentation
# regions that didn't survive curation. They are not cells, so fold them into
# 0 ("no cell") and let the unassigned filter drop them. Without this they were
# kept as a phantom cell literally named "unassigned"; on MER5-1_E1_01 that was
# 374 molecules and one extra cell versus the counts in the embryo metadata.
unnamed = np.array([c not in name_map for c in cell_ids])
if unnamed.any():
    ids = sorted({int(c) for c in cell_ids[unnamed]})
    print(f"  dropping {int(unnamed.sum()):,} molecules on {len(ids)} unnamed "
          f"cell id(s): {ids}")
    cell_ids = np.where(unnamed, 0, cell_ids).astype(np.int32)

cell_name = pd.Series(cell_ids).map(name_map).fillna("unassigned").to_numpy()

domain_anno = np.asarray(F["RNA_domain_anno"], dtype=object)[sel]

table = pa.table({
    # Xenium schema -> pickSchema() auto-detects, no manual remap needed.
    "feature_name": pa.array(labels[icodes[sel]].astype(str)).dictionary_encode(),
    "x_location": pa.array((Xh[sel, 1] * UM_PER_PX_XY).astype(np.float32)),
    "y_location": pa.array((Xh[sel, 2] * UM_PER_PX_XY).astype(np.float32)),
    "z_location": pa.array((Xh[sel, 0] * UM_PER_PX_Z).astype(np.float32)),
    # Extra per-molecule annotations. cell_id is read by the current pipeline
    # (assigned vs unassigned only); the rest are carried for the compartment
    # work and ignored for now.
    "cell_id": pa.array(cell_ids),
    "cell_name": pa.array(cell_name.astype(str)).dictionary_encode(),
    "rna_domain_id": pa.array(F["RNA_domain_id"][sel].astype(np.int32)),
    "rna_domain_anno": pa.array(domain_anno.astype(str)).dictionary_encode(),
})

import os
os.makedirs(OUT_DIR, exist_ok=True)
out = f"{OUT_DIR}/{EMB}_molecules.parquet"
pq.write_table(table, out, compression="snappy", row_group_size=500_000)

print(f"\nwrote {out}")
print(f"  rows        {table.num_rows:,}")
print(f"  genes       {len(np.unique(labels[icodes[sel]]))}")
print(f"  cells       {len(np.unique(cell_ids))} (incl. 0 = unassigned: {(cell_ids == 0).sum():,} molecules)")
print(f"  domains     {len(np.unique(domain_anno))} -> {sorted(np.unique(domain_anno))}")
print(f"  x µm        {table['x_location'].to_numpy().min():.1f} .. {table['x_location'].to_numpy().max():.1f}")
print(f"  y µm        {table['y_location'].to_numpy().min():.1f} .. {table['y_location'].to_numpy().max():.1f}")
print(f"  z µm        {table['z_location'].to_numpy().min():.1f} .. {table['z_location'].to_numpy().max():.1f}")
print(f"  size        {os.path.getsize(out) / 1e6:.1f} MB")
