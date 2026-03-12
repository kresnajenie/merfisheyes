"""
combine_slices.py

Discovers all cell metadata CSVs under a target folder using fuzzy matching,
separates their coordinates into a dynamic grid, and saves concatenated CSVs
back into the target folder.

Also concatenates the corresponding cell_by_gene.csv (purely positional, no
coordinate modification) and validates that row counts match cell_metadata.

Usage:
    python combine_slices.py /path/to/dataset_folder/

Output:
    <target_folder>/combined_slices.csv
    <target_folder>/combined_cell_by_gene.csv
"""

import argparse
import math
import sys
import re
from pathlib import Path
import pandas as pd
import numpy as np


# ─────────────────────────────────────────────
# FUZZY FILE MATCHING
# ─────────────────────────────────────────────
FILE_PATTERNS = {
    'cell_metadata': {
        'keywords': ['metadata'],
        'alternative_keywords': [['cell', 'meta']],
        'exclude': ['gene', 'detected', 'transcript'],
        'extensions': ['.csv'],
        'description': 'Cell metadata'
    },
    'cell_by_gene': {
        'keywords': ['cell', 'gene'],
        'exclude': ['metadata', 'detected', 'transcript'],
        'extensions': ['.csv'],
        'description': 'Cell by gene matrix'
    },
}


def normalize_filename(filename):
    """Lowercase, strip extension, replace separators with spaces."""
    name = Path(filename).stem.lower()
    return re.sub(r'[-_.]', ' ', name)


def match_file_pattern(filename, pattern_config):
    """
    Returns (score, reason). score > 0 means a match.
    """
    normalized = normalize_filename(filename)
    ext = Path(filename).suffix.lower()

    if ext not in pattern_config['extensions']:
        return 0, f"Extension {ext} not in {pattern_config['extensions']}"

    for excl in pattern_config.get('exclude', []):
        if excl in normalized:
            return 0, f"Contains excluded keyword '{excl}'"

    # Primary keywords — all must be present
    keywords = pattern_config['keywords']
    hits = [kw for kw in keywords
            if kw in normalized or kw + 's' in normalized
            or (kw.endswith('s') and kw[:-1] in normalized)]

    if len(hits) == len(keywords):
        return 100, f"Primary match: all keywords {keywords} found"

    # Alternative keyword combos
    for alt_combo in pattern_config.get('alternative_keywords', []):
        alt_hits = [kw for kw in alt_combo
                    if kw in normalized or kw + 's' in normalized
                    or (kw.endswith('s') and kw[:-1] in normalized)]
        if len(alt_hits) == len(alt_combo):
            return 80, f"Alternative match: keywords {alt_combo} found"

    return 0, f"Partial match only: found {hits}, missing {set(keywords) - set(hits)}"


def categorize_file(filename):
    """
    Returns (category, score, reason) or (None, 0, reason) if ambiguous/no match.
    """
    matches = []
    for pattern_name, pattern_config in FILE_PATTERNS.items():
        score, reason = match_file_pattern(filename, pattern_config)
        if score > 0:
            matches.append((pattern_name, score, reason))

    if not matches:
        return None, 0, "No pattern match"

    if len(matches) > 1:
        scores = [s for _, s, _ in matches]
        if max(scores) - min(scores) < 30:
            desc = ", ".join(f"{n} (score: {s})" for n, s, _ in matches)
            return None, 0, f"AMBIGUOUS: {desc}"

    best = max(matches, key=lambda x: x[1])
    return best


def find_file_in_dir(directory, target_category):
    """
    Scan a directory for a file matching target_category using fuzzy matching.
    Returns (Path, reason) or (None, reason) if not found.
    Exits on ambiguity (multiple matches).
    """
    candidates = []
    for f in directory.iterdir():
        if not f.is_file():
            continue
        cat, score, reason = categorize_file(f.name)
        if cat == target_category:
            candidates.append((f, score, reason))

    if len(candidates) == 1:
        f, score, reason = candidates[0]
        return f, reason

    if len(candidates) > 1:
        print(f"ERROR: Multiple files matched '{target_category}' in {directory}:")
        for f, score, reason in candidates:
            print(f"  {f.name}  (score: {score}, {reason})")
        sys.exit(1)

    return None, "Not found"


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Separate overlapping tissue sample coordinates.")
parser.add_argument("target_folder", help="Path to the folder containing sample sub-folders.")
# changed
parser.add_argument("output_dir", 
                    help="Where to create combined_output/. Defaults to target folder.")
parser.add_argument("--padding", type=float, default=1000,
                    help="Gap between samples in coordinate units (default: 1000)")


args = parser.parse_args()

target = Path.cwd() / args.target_folder

if not target.is_dir():
    print(f"ERROR: '{target}' is not a valid directory.")
    sys.exit(1)


# ─────────────────────────────────────────────
# DISCOVER DIRECTORIES CONTAINING cell_metadata
# ─────────────────────────────────────────────
candidate_dirs = set()
for f in target.rglob("*"):
    if not f.is_file():
        continue
    cat, score, reason = categorize_file(f.name)
    if cat == "cell_metadata":
        candidate_dirs.add(f.parent)

candidate_dirs = sorted(candidate_dirs)

if not candidate_dirs:
    print(f"ERROR: No cell metadata files found under '{target}'.")
    sys.exit(1)

print(f"Found {len(candidate_dirs)} sample directory(ies):")
for d in candidate_dirs:
    print(f"  {d}")


# ─────────────────────────────────────────────
# LOAD cell_metadata AND cell_by_gene PER DIR
# ─────────────────────────────────────────────
dfs     = []   # cell_metadata DataFrames
cbg_dfs = []   # cell_by_gene DataFrames

for d in candidate_dirs:
    # --- cell_metadata ---
    meta_path, meta_reason = find_file_in_dir(d, "cell_metadata")
    if meta_path is None:
        print(f"ERROR: No cell_metadata file found in {d}")
        sys.exit(1)

    df = pd.read_csv(meta_path)
    df["_source_file"] = str(meta_path)
    df["_sample_id"] = meta_path.relative_to(target).parts[0]
    dfs.append(df)
    print(f"  cell_metadata  : {meta_path.name}  ({meta_reason})")

    # --- cell_by_gene ---
    cbg_path, cbg_reason = find_file_in_dir(d, "cell_by_gene")
    if cbg_path is None:
        print(f"ERROR: No cell_by_gene file found in {d}")
        sys.exit(1)

    cbg_df = pd.read_csv(cbg_path)

    # Row count validation
    if len(cbg_df) != len(df):
        print(
            f"ERROR: Row count mismatch in {d}:\n"
            f"  {meta_path.name}  -> {len(df):,} rows\n"
            f"  {cbg_path.name}   -> {len(cbg_df):,} rows"
        )
        sys.exit(1)

    cbg_dfs.append(cbg_df)
    print(f"  cell_by_gene   : {cbg_path.name}  ({cbg_reason})")

print(f"\nAll files found and row counts validated.")


# ─────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────
x_col = "center_x"
y_col = "center_y"

bbox_cols = {
    "min_x": "x",
    "max_x": "x",
    "min_y": "y",
    "max_y": "y",
}

padding = args.padding


# ─────────────────────────────────────────────
# BOUNDING BOXES
# ─────────────────────────────────────────────
def get_bbox(df):
    return (
        df[x_col].min(), df[y_col].min(),
        df[x_col].max(), df[y_col].max(),
    )

bboxes  = [get_bbox(df) for df in dfs]
widths  = [b[2] - b[0] for b in bboxes]
heights = [b[3] - b[1] for b in bboxes]


# ─────────────────────────────────────────────
# DYNAMIC GRID
# ─────────────────────────────────────────────
n      = len(dfs)
n_cols = math.ceil(math.sqrt(n))
n_rows = math.ceil(n / n_cols)

col_widths = []
for c in range(n_cols):
    indices = [i for i in range(n) if i % n_cols == c]
    col_widths.append(max(widths[i] for i in indices) + padding)

row_heights = []
for r in range(n_rows):
    indices = [i for i in range(n) if i // n_cols == r]
    row_heights.append(max(heights[i] for i in indices) + padding)

col_offsets = [sum(col_widths[:c])  for c in range(n_cols)]
row_offsets = [sum(row_heights[:r]) for r in range(n_rows)]

grid_origins = []
for i in range(n):
    row = i // n_cols
    col = i %  n_cols
    grid_origins.append((col_offsets[col], -row_offsets[row]))


# ─────────────────────────────────────────────
# SHIFT COORDINATES
# ─────────────────────────────────────────────
def shift_df(df, origin, old_bbox):
    df = df.copy()
    dx = origin[0] - old_bbox[0]
    dy = origin[1] - old_bbox[1]

    df[x_col] += dx
    df[y_col] += dy

    for col, axis in bbox_cols.items():
        if col in df.columns:
            df[col] += dx if axis == "x" else dy

    return df

separated_dfs = [
    shift_df(df, origin, bb)
    for df, origin, bb in zip(dfs, grid_origins, bboxes)
]


# ─────────────────────────────────────────────
# CONCATENATE & SAVE cell_metadata
# ─────────────────────────────────────────────
combined = pd.concat(separated_dfs, ignore_index=True)

out_dir = Path.cwd() / args.output_dir / "combined_output"
print("OUTIDIR", out_dir)
out_dir.mkdir(parents = True, exist_ok=True)

out_path = out_dir / "cell_metadata.csv"
combined.to_csv(out_path, index=False)

print(f"\nDone. {n} samples | {n_rows}x{n_cols} grid | {len(combined):,} total rows")
print(f"Output -> {out_path}")

# SANITY CHECK PLOTS
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

cmap   = plt.colormaps.get_cmap("tab20")
colors = [cmap(i / max(n - 1, 1)) for i in range(n)]

# Plot 1: spatial layout
fig, ax = plt.subplots(figsize=(10, 10))
for i, df in enumerate(separated_dfs):
    ax.scatter(df[x_col], df[y_col], c=[colors[i]], s=0.5, alpha=0.4, rasterized=True)
ax.set_title(f"Cell metadata spatial layout — {len(combined):,} cells")
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig(out_dir / "check_spatial.png", dpi=150)
plt.close()

# ─────────────────────────────────────────────
# CONCATENATE & SAVE cell_by_gene
# ─────────────────────────────────────────────
combined_cbg = pd.concat(cbg_dfs, ignore_index=True)

# alignment is positional — both lists were built in the same sample order
# just verify row counts match, don't sort by ID since IDs repeat across samples
assert len(combined) == len(combined_cbg), \
    f"ERROR: Row count mismatch: metadata={len(combined):,}, cell_by_gene={len(combined_cbg):,}"

out_cbg_path = out_dir / "cell_by_gene.csv"
combined_cbg.to_csv(out_cbg_path, index=False)

print(f"Output -> {out_cbg_path} ({len(combined_cbg):,} total rows)")

# Plot 2: raw cell_by_gene total counts per cell
gene_cols = [c for c in combined_cbg.columns if c != "cell"]
gene_sums = combined_cbg[gene_cols].sum(axis=1).values
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(gene_sums, bins=100, color="steelblue", edgecolor="none")
ax.axvline(np.median(gene_sums), color="red", linestyle="--", label=f"median={np.median(gene_sums):.1f}")
ax.set_title("Distribution of total gene counts per cell")
ax.set_xlabel("total gene counts")
ax.set_ylabel("# cells")
ax.legend()
plt.tight_layout()
plt.savefig(out_dir / "check_expression.png", dpi=150)
plt.close()
