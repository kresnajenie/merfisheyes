"""Load one spiralia embryo in napari, following final_analyzed_obj_explore.ipynb.

Usage: python view_napari.py [EMBRYO_ID] [SAMPLE_FACTOR]
"""
import sys
import dill
import numpy as np
import pandas as pd
import plotly.express as px
import napari

ROOT = "/home/data/yiqun-spiralia/Sep2026"
EMB = sys.argv[1] if len(sys.argv) > 1 else "MER6-2_E3_1"
SAMPLE = int(sys.argv[2]) if len(sys.argv) > 2 else 40
RESIZE = np.array([2, 4, 4])          # segm is 2,4,4 downsampled vs Xh
LAYER_SCALE = (1.2, -1, 1)


def load_object(filename):
    with open(filename, "rb") as f:
        return dill.load(f)


print(f"loading {EMB} ...", flush=True)
FISH = load_object(f"{ROOT}/Final_analyzed_objects/{EMB}/{EMB}_RNA_domain_annotated")
nuc = np.load(f"{ROOT}/Segmentation/Bogdan/Old/{EMB}_segm.npz")["segm"]
print(f"  {len(FISH['Xh'])} molecules, nuc mask {nuc.shape}", flush=True)

colors_seq = np.array(px.colors.qualitative.Light24 + px.colors.qualitative.Dark24)

# nucleus centres of mass, for the text labels
cell_sheet = pd.DataFrame(FISH["cell_names"])
cell_sheet["cell_name"] = [i.split("::")[1] for i in np.array(cell_sheet.iloc[:, 0])]
for c in "xyz":
    cell_sheet[c] = 0.0
for i in np.unique(nuc):
    if i in cell_sheet.index:
        ind = np.where(nuc == i)
        cell_sheet.loc[i, "z"] = np.mean(ind[0])
        cell_sheet.loc[i, "x"] = np.mean(ind[1])
        cell_sheet.loc[i, "y"] = np.mean(ind[2])

keep = FISH["molecules_keep"][::SAMPLE]
X_plt = (FISH["Xh"][::SAMPLE, :3] / RESIZE)[keep]

V = napari.Viewer(title=f"{EMB}")
for colorby, key in [("domain", "RNA_domain_id"), ("cell", "cell_ids")]:
    vals = FISH[key][::SAMPLE][keep]
    # cell_ids are non-contiguous (…24, 5000–5006), so compact them before
    # indexing the palette; the notebook indexes raw and overflows here.
    _, compact = np.unique(vals, return_inverse=True)
    face = colors_seq[compact % len(colors_seq)]
    V.add_points(X_plt, face_color=face, size=1.5, opacity=0.6, border_width=0,
                 border_color="white", name=f"{EMB}::{colorby}",
                 visible=(colorby == "domain"))

V.add_points(cell_sheet[["z", "x", "y"]].to_numpy(),
             text={"string": list(cell_sheet["cell_name"]), "color": "red"},
             size=12, face_color="white", border_width=0, name="cell_name")
V.add_labels(nuc, opacity=0.6, name="nuclei")

for layer in V.layers:
    layer.scale = LAYER_SCALE

print("viewer ready", flush=True)
napari.run()
