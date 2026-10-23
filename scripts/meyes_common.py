"""Rules every script that prepares data shares, so they cannot drift apart:
which features are not genes, and which columns of a table are molecules.

Used by scripts/prepare.py (and the prepare_* scripts it runs), the server
worker (worker/entrypoint.py) and the legacy processors.
"""
from __future__ import annotations

import re

# Control probes / unassigned codewords are not genes. Mirrors
# lib/utils/gene-filters.ts shouldFilterGene() — the two lists must stay equal.
_GENE_FILTER_PATTERNS = [
    re.compile(p, re.IGNORECASE)
    for p in (
        r"negative[\s_-]*control",
        r"neg[\s_-]*ctrl",
        r"unassigned",
        r"deprecated",
        r"codeword",
        r"blank",
        r"negcontrol",
    )
]


def should_filter_gene(gene: str) -> bool:
    g = str(gene)
    return any(p.search(g) for p in _GENE_FILTER_PATTERNS)


# The names a molecule table's columns may have, each list tried in order and
# without regard to case: a platform's own (Xenium, MERSCOPE) and a lab's own
# (Gene, X, Y, Z). worker/entrypoint.py and lib/ingest/molecule-columns.ts
# recognise the same names — tests keep the three equal.
MOLECULE_COLUMN_NAMES = {
    "gene": ["gene", "feature_name", "feature", "target", "gene_name", "target_name", "genes"],
    "x": ["x", "x_location", "global_x", "x_um", "x_global", "x_position"],
    "y": ["y", "y_location", "global_y", "y_um", "y_global", "y_position"],
    "z": ["z", "z_location", "global_z", "z_um", "z_global", "z_position"],
}


def match_molecule_columns(columns) -> dict | None:
    """gene / x / y (and z when there is one) among `columns`, by name in any
    case; None without a gene, an x and a y column."""
    by_name: dict[str, str] = {}
    for column in columns:
        by_name.setdefault(str(column).strip().lower(), column)

    found = {}
    for role, names in MOLECULE_COLUMN_NAMES.items():
        match = next((by_name[n] for n in names if n in by_name), None)
        if match is not None:
            found[role] = match

    return found if {"gene", "x", "y"} <= found.keys() else None


# Where a table without obsm['spatial'] keeps its cells' coordinates: obs
# columns, each list tried in order — the same as process_spatial_data.py's
# fallback and the browser's H5adAdapter.
OBS_COORDINATE_COLUMNS = {
    "x": ["center_x", "centerX", "centerx", "x", "X", "x_centroid", "centroid_x", "x_location", "center_X", "spatial_x", "spatial_X"],
    "y": ["center_y", "centerY", "centery", "y", "Y", "y_centroid", "centroid_y", "y_location", "center_Y", "spatial_y", "spatial_Y"],
    "z": ["center_z", "centerZ", "centerz", "z", "Z", "z_centroid", "centroid_z", "z_location", "center_Z", "spatial_z", "spatial_Z"],
}


def obs_coordinate_columns(columns) -> list[str] | None:
    """[x, y] or [x, y, z] among an obs frame's columns; None without x and y."""
    names = set(columns)
    found = [next((c for c in OBS_COORDINATE_COLUMNS[a] if c in names), None) for a in ("x", "y", "z")]
    if found[0] is None or found[1] is None:
        return None
    return found[:2] + ([found[2]] if found[2] else [])
