#!/usr/bin/env python3
"""Build the tiny .h5ad files the browser preview used to misread.

  python scripts/testdata/make-tiny-h5ad-edge.py

Outputs, in tests/data/tiny/h5ad/:

  nullable.h5ad     the encodings anndata 0.13 + pandas 3 write: the var and
                    obs indexes and a string obs column as
                    nullable-string-array groups (values + mask), a nullable
                    integer and a nullable boolean column, one of each with a
                    missing entry. The preview read none of them (0 genes).
  named_index.h5ad  a var index named `gene_ids` (Ensembl ids) beside a `gene`
                    symbol column: var's `_index` attribute names the index.
  missing_xy.h5ad   24 cells, no obsm: positions in obs center_x / center_y, with rows
                    3 and 7 missing an x. The preview dropped those rows and
                    shifted every later position onto the wrong cell.

Needs anndata >= 0.11 (scripts/requirements.txt pins 0.13).
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "tests" / "data" / "tiny" / "h5ad"
N, GENES = 12, ["Gad1", "Sst", "Pvalb", "Vip"]
rng = np.random.default_rng(7)
X = rng.poisson(2.0, (N, len(GENES))).astype(np.float32)
xy = np.c_[np.arange(N) * 10.0, np.arange(N) * 5.0]

# Write the nullable encodings whatever pandas' default string dtype is
if hasattr(ad.settings, "allow_write_nullable_strings"):
    ad.settings.allow_write_nullable_strings = True

strings = lambda values: pd.array(values, dtype="string")  # noqa: E731

nullable = ad.AnnData(
    X=X,
    obs=pd.DataFrame(
        {
            "celltype": pd.Categorical(["a", "b", "c"] * (N // 3)),
            "region": strings(["cortex", "striatum", None] + ["cortex"] * (N - 3)),
            "count": pd.array([5, None] + list(range(N - 2)), dtype="Int64"),
            "flag": pd.array([True, False, None] + [True] * (N - 3), dtype="boolean"),
        },
        index=pd.Index(strings([f"cell{i}" for i in range(N)])),
    ),
    var=pd.DataFrame(index=pd.Index(strings(GENES))),
    obsm={"spatial": xy},
)
nullable.write_h5ad(OUT / "nullable.h5ad")

named = ad.AnnData(
    X=X,
    obs=pd.DataFrame(
        {"celltype": pd.Categorical(["a", "b", "c"] * (N // 3))},
        index=pd.Index([f"cell{i}" for i in range(N)], dtype=object),
    ),
    var=pd.DataFrame(
        {"gene": pd.array(GENES, dtype=object)},
        index=pd.Index([f"ENSMUSG{i:011d}" for i in range(len(GENES))], dtype=object, name="gene_ids"),
    ),
    obsm={"spatial": xy},
)
named.write_h5ad(OUT / "named_index.h5ad")

# 24 cells, so two missing still clear the preview's 90%-valid rule
M = 24
x = np.arange(M) * 10.0
x[[3, 7]] = np.nan
missing = ad.AnnData(
    X=rng.poisson(2.0, (M, len(GENES))).astype(np.float32),
    obs=pd.DataFrame(
        {"celltype": pd.Categorical(["a", "b", "c"] * (M // 3)), "center_x": x, "center_y": np.arange(M) * 5.0},
        index=pd.Index([f"cell{i}" for i in range(M)], dtype=object),
    ),
    var=pd.DataFrame(index=pd.Index(GENES, dtype=object)),
)
missing.write_h5ad(OUT / "missing_xy.h5ad")

for name in ("nullable.h5ad", "named_index.h5ad", "missing_xy.h5ad"):
    print(f"wrote {OUT / name} ({(OUT / name).stat().st_size:,} bytes)")
