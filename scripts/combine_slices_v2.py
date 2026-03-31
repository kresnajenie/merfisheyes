"""
combine_slices_v2.py

Same as combine_slices.py but uses BFS with os.scandir() for fast file
discovery. Avoids traversing directories with thousands of non-CSV files
(e.g. imaging data).

For each child of the target folder, does a BFS:
  - Single scandir pass per directory: collect CSVs + subdirectories
  - If both cell_metadata and cell_by_gene CSVs found → record, skip subtree
  - Otherwise queue subdirectories for next level

Usage:
    python combine_slices_v2.py /path/to/dataset_folder/ /path/to/output_dir/

Output:
    <output_dir>/combined_output/cell_metadata.csv
    <output_dir>/combined_output/cell_by_gene.csv
"""

import argparse
import math
import os
import sys
import re
import time
from collections import deque
from pathlib import Path
import pandas as pd
import numpy as np


def fmt_elapsed(seconds):
    """Format elapsed time as human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(seconds, 60)
    return f"{int(m)}m {s:.1f}s"


def log(msg, t0=None):
    """Print a timestamped log message. If t0 given, also prints elapsed."""
    elapsed = ""
    if t0 is not None:
        elapsed = f" [{fmt_elapsed(time.perf_counter() - t0)}]"
    print(f"[{time.strftime('%H:%M:%S')}]{elapsed} {msg}", flush=True)


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
    Only checks .csv files via scandir (skips non-CSV entirely).
    Returns (Path, reason) or (None, reason) if not found.
    Exits on ambiguity (multiple matches).
    """
    candidates = []
    for entry in os.scandir(directory):
        if not entry.is_file():
            continue
        if not entry.name.lower().endswith('.csv'):
            continue
        cat, score, reason = categorize_file(entry.name)
        if cat == target_category:
            candidates.append((Path(entry.path), score, reason))

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
# BFS DIRECTORY DISCOVERY
# ─────────────────────────────────────────────
def scan_dir_once(directory):
    """
    Single os.scandir pass. Returns (csv_names, subdirs).
    Only collects .csv file names and subdirectory paths.
    Skips all non-CSV files without further inspection.
    """
    csv_names = []
    subdirs = []
    try:
        for entry in os.scandir(directory):
            if entry.is_dir(follow_symlinks=False):
                subdirs.append(entry.path)
            elif entry.is_file(follow_symlinks=False) and entry.name.lower().endswith('.csv'):
                csv_names.append(entry.name)
    except PermissionError:
        pass
    return csv_names, subdirs


def has_both_targets(csv_names):
    """Check if a list of CSV filenames contains both cell_metadata and cell_by_gene."""
    found_categories = set()
    for name in csv_names:
        cat, score, reason = categorize_file(name)
        if cat in ('cell_metadata', 'cell_by_gene'):
            found_categories.add(cat)
        if len(found_categories) == 2:
            return True
    return False


def discover_sample_dirs_bfs(target, t0):
    """
    BFS from each child of target. For each child branch:
    - Scan directory: collect CSVs + subdirs in one pass
    - If both targets found → record directory, skip subtree
    - Otherwise → queue subdirs
    """
    candidate_dirs = []

    # Get top-level children
    top_children = []
    for entry in os.scandir(target):
        if entry.is_dir(follow_symlinks=False):
            top_children.append((entry.name, entry.path))
    top_children.sort()

    log(f"Scanning {len(top_children)} top-level directories...", t0)

    for idx, (child_name, child_path) in enumerate(top_children, 1):
        log(f"  [{idx}/{len(top_children)}] Searching '{child_name}'...", t0)

        # BFS within this child branch
        queue = deque([child_path])
        dirs_scanned = 0

        while queue:
            current = queue.popleft()
            csv_names, subdirs = scan_dir_once(current)
            dirs_scanned += 1

            if has_both_targets(csv_names):
                candidate_dirs.append(Path(current))
                log(f"    FOUND at {current} (scanned {dirs_scanned} dirs)", t0)
                break

            # queue subdirs for next level
            subdirs.sort()
            queue.extend(subdirs)
        else:
            log(f"    not found (scanned {dirs_scanned} dirs)", t0)

    return candidate_dirs


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Separate overlapping tissue sample coordinates.")
parser.add_argument("target_folder", help="Path to the folder containing sample sub-folders.")
parser.add_argument("output_dir",
                    help="Where to create combined_output/. Defaults to target folder.")
parser.add_argument("--padding", type=float, default=1000,
                    help="Gap between samples in coordinate units (default: 1000)")


args = parser.parse_args()
t_start = time.perf_counter()

target = Path.cwd() / args.target_folder

if not target.is_dir():
    print(f"ERROR: '{target}' is not a valid directory.")
    sys.exit(1)

log(f"Target folder: {target}", t_start)


# ─────────────────────────────────────────────
# DISCOVER DIRECTORIES CONTAINING cell_metadata
# ─────────────────────────────────────────────
log("=== STEP 1/6: Discovering sample directories (BFS) ===", t_start)
candidate_dirs = discover_sample_dirs_bfs(target, t_start)

if not candidate_dirs:
    log(f"ERROR: No sample directories found under '{target}'.", t_start)
    sys.exit(1)

log(f"Discovery complete. Found {len(candidate_dirs)} sample(s):", t_start)
for d in candidate_dirs:
    log(f"  {d}")


# ─────────────────────────────────────────────
# LOAD cell_metadata AND cell_by_gene PER DIR
# ─────────────────────────────────────────────
log(f"=== STEP 2/6: Loading CSV files ({len(candidate_dirs)} samples) ===", t_start)
dfs     = []   # cell_metadata DataFrames
cbg_dfs = []   # cell_by_gene DataFrames

for i, d in enumerate(candidate_dirs, 1):
    log(f"  [{i}/{len(candidate_dirs)}] Loading from {d.name}...", t_start)

    # --- cell_metadata ---
    meta_path, meta_reason = find_file_in_dir(d, "cell_metadata")
    if meta_path is None:
        log(f"ERROR: No cell_metadata file found in {d}", t_start)
        sys.exit(1)

    t_read = time.perf_counter()
    df = pd.read_csv(meta_path)
    df["_source_file"] = str(meta_path)
    df["_sample_id"] = meta_path.relative_to(target).parts[0]
    dfs.append(df)
    log(f"    cell_metadata  : {meta_path.name} ({len(df):,} rows, read in {fmt_elapsed(time.perf_counter() - t_read)})", t_start)

    # --- cell_by_gene ---
    cbg_path, cbg_reason = find_file_in_dir(d, "cell_by_gene")
    if cbg_path is None:
        log(f"ERROR: No cell_by_gene file found in {d}", t_start)
        sys.exit(1)

    t_read = time.perf_counter()
    cbg_df = pd.read_csv(cbg_path)
    log(f"    cell_by_gene   : {cbg_path.name} ({len(cbg_df):,} rows, read in {fmt_elapsed(time.perf_counter() - t_read)})", t_start)

    # Row count validation
    if len(cbg_df) != len(df):
        log(
            f"ERROR: Row count mismatch in {d}:\n"
            f"  {meta_path.name}  -> {len(df):,} rows\n"
            f"  {cbg_path.name}   -> {len(cbg_df):,} rows",
            t_start
        )
        sys.exit(1)

    cbg_dfs.append(cbg_df)

total_rows = sum(len(df) for df in dfs)
log(f"All files loaded and validated. Total rows: {total_rows:,}", t_start)


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
# BOUNDING BOXES & GRID
# ─────────────────────────────────────────────
log("=== STEP 3/6: Computing grid layout ===", t_start)

def get_bbox(df):
    return (
        df[x_col].min(), df[y_col].min(),
        df[x_col].max(), df[y_col].max(),
    )

bboxes  = [get_bbox(df) for df in dfs]
widths  = [b[2] - b[0] for b in bboxes]
heights = [b[3] - b[1] for b in bboxes]

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

log(f"Grid: {n_rows} rows x {n_cols} cols, padding={padding}", t_start)


# ─────────────────────────────────────────────
# SHIFT COORDINATES
# ─────────────────────────────────────────────
log("=== STEP 4/6: Shifting coordinates ===", t_start)

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

log("Coordinates shifted.", t_start)


# ─────────────────────────────────────────────
# CONCATENATE & SAVE cell_metadata
# ─────────────────────────────────────────────
log("=== STEP 5/6: Saving cell_metadata.csv ===", t_start)
combined = pd.concat(separated_dfs, ignore_index=True)

out_dir = Path.cwd() / args.output_dir / "combined_output"
out_dir.mkdir(parents = True, exist_ok=True)

out_path = out_dir / "cell_metadata.csv"
t_write = time.perf_counter()
combined.to_csv(out_path, index=False)
log(f"Saved {out_path} ({len(combined):,} rows, wrote in {fmt_elapsed(time.perf_counter() - t_write)})", t_start)


# ─────────────────────────────────────────────
# CONCATENATE & SAVE cell_by_gene
# ─────────────────────────────────────────────
log("=== STEP 6/6: Saving cell_by_gene.csv ===", t_start)
combined_cbg = pd.concat(cbg_dfs, ignore_index=True)

assert len(combined) == len(combined_cbg), \
    f"ERROR: Row count mismatch: metadata={len(combined):,}, cell_by_gene={len(combined_cbg):,}"

out_cbg_path = out_dir / "cell_by_gene.csv"
t_write = time.perf_counter()
combined_cbg.to_csv(out_cbg_path, index=False)
log(f"Saved {out_cbg_path} ({len(combined_cbg):,} rows, wrote in {fmt_elapsed(time.perf_counter() - t_write)})", t_start)


# ─────────────────────────────────────────────
# SANITY CHECK PLOTS
# ─────────────────────────────────────────────
log("Generating sanity check plots...", t_start)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
log(f"Saved {out_dir / 'check_spatial.png'}", t_start)

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
log(f"Saved {out_dir / 'check_expression.png'}", t_start)

# ─────────────────────────────────────────────
# DONE
# ─────────────────────────────────────────────
log(f"=== COMPLETE ===", t_start)
log(f"  Samples: {n}", t_start)
log(f"  Grid: {n_rows}x{n_cols}", t_start)
log(f"  Total rows: {len(combined):,}", t_start)
log(f"  Output dir: {out_dir}", t_start)
log(f"  Total time: {fmt_elapsed(time.perf_counter() - t_start)}", t_start)

