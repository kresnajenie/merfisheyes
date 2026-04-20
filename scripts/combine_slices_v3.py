"""
combine_slices_v3.py

Memory-optimized version. Writes CSVs incrementally (one sample at a time)
instead of concatenating everything in memory. Peak memory ~ 2 samples
instead of all samples x2.

Key changes from v2:
  - Pass 1: Read each cell_metadata, compute bounding boxes, record file paths
            (only keep bboxes + paths, not DataFrames)
  - Pass 2: Compute grid layout from bboxes
  - Pass 3: Re-read each cell_metadata one at a time, shift coords, append to output CSV
  - Pass 4: Append each cell_by_gene one at a time to output CSV (no coord changes)
  - Sanity plots use the already-written combined CSV (chunked read for histogram)

Usage:
    python combine_slices_v3.py /path/to/dataset_folder/ /path/to/output_dir/

Output:
    <output_dir>/combined_output/cell_metadata.csv
    <output_dir>/combined_output/cell_by_gene.csv
"""

import argparse
import gc
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
                    help="Where to create combined_output/.")
parser.add_argument("--padding", type=float, default=1000,
                    help="Gap between samples in coordinate units (default: 1000)")
parser.add_argument("--mask-only", action="store_true",
                    help="Skip combining. target_folder should point to an existing "
                         "combined_output/ directory with cell_metadata.csv and "
                         "cell_by_gene.csv already present. Runs only step 6b onward.")


args = parser.parse_args()
t_start = time.perf_counter()

target = Path.cwd() / args.target_folder

if not target.is_dir():
    print(f"ERROR: '{target}' is not a valid directory.")
    sys.exit(1)

log(f"Target folder: {target}", t_start)

x_col = "center_x"
y_col = "center_y"

if args.mask_only:
    # ── MASK-ONLY MODE: skip combining, jump straight to 6b ──
    out_dir = Path(target)
    out_meta_path = out_dir / "cell_metadata.csv"
    out_cbg_path = out_dir / "cell_by_gene.csv"

    if not out_meta_path.exists() or not out_cbg_path.exists():
        log(f"ERROR: --mask-only requires cell_metadata.csv and cell_by_gene.csv in {out_dir}", t_start)
        sys.exit(1)

    # Count rows for downstream logging
    rows_written = sum(len(chunk) for chunk in pd.read_csv(out_meta_path, usecols=[0], chunksize=100_000))
    log(f"Mask-only mode: found {rows_written:,} rows in existing combined CSVs", t_start)

else:

    # ─────────────────────────────────────────────
    # STEP 1: DISCOVER DIRECTORIES
    # ─────────────────────────────────────────────
    log("=== STEP 1/8: Discovering sample directories (BFS) ===", t_start)
    candidate_dirs = discover_sample_dirs_bfs(target, t_start)
    
    if not candidate_dirs:
        log(f"ERROR: No sample directories found under '{target}'.", t_start)
        sys.exit(1)
    
    log(f"Discovery complete. Found {len(candidate_dirs)} sample(s):", t_start)
    for d in candidate_dirs:
        log(f"  {d}")
    
    
    # ─────────────────────────────────────────────
    # STEP 2: FIRST PASS — read cell_metadata for
    #   bounding boxes only, record file paths,
    #   validate row counts, then FREE the DataFrame
    # ─────────────────────────────────────────────
    log(f"=== STEP 2/8: First pass — bounding boxes & validation ({len(candidate_dirs)} samples) ===", t_start)

    bbox_cols = {
        "min_x": "x",
        "max_x": "x",
        "min_y": "y",
        "max_y": "y",
    }
    
    padding = args.padding
    
    # Store lightweight info per sample (no DataFrames kept)
    sample_info = []  # list of dicts: {meta_path, cbg_path, bbox, num_rows, sample_id}
    
    for i, d in enumerate(candidate_dirs, 1):
        log(f"  [{i}/{len(candidate_dirs)}] Scanning {d.name}...", t_start)
    
        # --- find files ---
        meta_path, meta_reason = find_file_in_dir(d, "cell_metadata")
        if meta_path is None:
            log(f"ERROR: No cell_metadata file found in {d}", t_start)
            sys.exit(1)
    
        cbg_path, cbg_reason = find_file_in_dir(d, "cell_by_gene")
        if cbg_path is None:
            log(f"ERROR: No cell_by_gene file found in {d}", t_start)
            sys.exit(1)
    
        # --- read cell_metadata ---
        t_read = time.perf_counter()
        df = pd.read_csv(meta_path)
        meta_rows_raw = len(df)
    
        # Find metadata ID column (EntityID or id)
        metadata_id_col = None
        for candidate in ['EntityID', 'id', 'cell_id']:
            if candidate in df.columns:
                metadata_id_col = candidate
                break
        if metadata_id_col is None:
            first_col = df.columns[0]
            log(f"WARNING: No named ID column in {meta_path.name}. "
                f"Falling back to first column: '{first_col}'. Renaming to 'cell_id'.", t_start)
            df.rename(columns={first_col: 'cell_id'}, inplace=True)
            metadata_id_col = 'cell_id'
    
        # Validate metadata IDs are numeric
        df[metadata_id_col] = df[metadata_id_col].astype(str)
    
        meta_id_set = set(df[metadata_id_col].values)
        meta_columns = list(df.columns)
        log(f"    cell_metadata  : {meta_path.name} ({meta_rows_raw:,} rows, "
            f"ID col='{metadata_id_col}', read in {fmt_elapsed(time.perf_counter() - t_read)})", t_start)
    
        # --- read cell_by_gene cell IDs ---
        t_read = time.perf_counter()
        cbg_all_cols = pd.read_csv(cbg_path, nrows=0).columns.tolist()
        if 'cell' in cbg_all_cols:
            cbg_id_col = 'cell'
        else:
            cbg_id_col = cbg_all_cols[0]
            log(f"WARNING: No named ID column in {cbg_path.name}. "
                f"Falling back to first column: '{cbg_id_col}'. Renaming to 'cell'.", t_start)
            
        id_chunks = []
        for chunk in pd.read_csv(cbg_path, usecols=[cbg_id_col], chunksize=100_000):
            if cbg_id_col != 'cell':
                chunk = chunk.rename(columns={cbg_id_col: 'cell'})
            id_chunks.append(chunk['cell'])
        expr_ids = pd.concat(id_chunks, ignore_index=True)
        del id_chunks
    
        expr_ids = expr_ids.astype(str)
    
        cbg_rows_raw = len(expr_ids)
        expr_id_set = set(expr_ids.values)
        log(f"    cell_by_gene   : {cbg_path.name} ({cbg_rows_raw:,} rows, "
            f"IDs read in {fmt_elapsed(time.perf_counter() - t_read)})", t_start)
    
        # --- Inner join on cell IDs ---
        common_ids = expr_id_set & meta_id_set
        if len(common_ids) == 0:
            log(f"ERROR: No matching cell IDs between {meta_path.name} (col '{metadata_id_col}') "
                f"and {cbg_path.name} (col 'cell'). "
                f"Meta IDs (first 5): {sorted(meta_id_set)[:5]}, "
                f"Expr IDs (first 5): {sorted(expr_id_set)[:5]}", t_start)
            sys.exit(1)
    
        expr_only = expr_id_set - common_ids
        meta_only = meta_id_set - common_ids
        log(f"    ID alignment   : {len(common_ids):,} common, "
            f"{len(expr_only):,} expr-only, {len(meta_only):,} meta-only", t_start)
    
        if expr_only:
            sample_ids_list = sorted(expr_only)[:20]
            log(f"      ⚠ Discarding {len(expr_only):,} from {cbg_path.name}: "
                f"{sample_ids_list}{'...' if len(expr_only) > 20 else ''}", t_start)
        if meta_only:
            sample_ids_list = sorted(meta_only)[:20]
            log(f"      ⚠ Discarding {len(meta_only):,} from {meta_path.name}: "
                f"{sample_ids_list}{'...' if len(meta_only) > 20 else ''}", t_start)
    
        # Valid cell IDs in cell_by_gene order
        mask = expr_ids.isin(common_ids)
        valid_cell_ids = expr_ids[mask].values
        del expr_ids
        gc.collect()
    
        # Filter metadata to aligned cells, compute bbox
        df = df.set_index(metadata_id_col).loc[valid_cell_ids].reset_index()
        bbox = (
            df[x_col].min(), df[y_col].min(),
            df[x_col].max(), df[y_col].max(),
        )
        log(f"    bbox (aligned) : x=[{bbox[0]:.1f}, {bbox[2]:.1f}] y=[{bbox[1]:.1f}, {bbox[3]:.1f}]", t_start)
    
        del df
        gc.collect()
    
        sample_info.append({
            'meta_path': meta_path,
            'cbg_path': cbg_path,
            'bbox': bbox,
            'num_rows': len(valid_cell_ids),
            'sample_id': meta_path.relative_to(target).parts[0],
            'meta_columns': meta_columns,
            'metadata_id_col': metadata_id_col,
            'valid_cell_ids': valid_cell_ids,
            'cbg_id_col': cbg_id_col
        })
    
    total_rows = sum(s['num_rows'] for s in sample_info)
    log(f"All samples validated. Total rows: {total_rows:,}", t_start)
    
    
    # ─────────────────────────────────────────────
    # STEP 3: COMPUTE GRID LAYOUT (from bboxes)
    # ─────────────────────────────────────────────
    log("=== STEP 3/8: Computing grid layout ===", t_start)
    
    n = len(sample_info)
    bboxes = [s['bbox'] for s in sample_info]
    widths = [b[2] - b[0] for b in bboxes]
    heights = [b[3] - b[1] for b in bboxes]
    
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
    
    col_offsets = [sum(col_widths[:c]) for c in range(n_cols)]
    row_offsets = [0]
    for r in range(1, n_rows):
        row_offsets.append(row_offsets[-1] + row_heights[r])
    
    grid_origins = []
    for i in range(n):
        row = i // n_cols
        col = i % n_cols
        grid_origins.append((col_offsets[col], -row_offsets[row]))
    
    log(f"Grid: {n_rows} rows x {n_cols} cols, padding={padding}", t_start)
    
    
    # ─────────────────────────────────────────────
    # STEP 4: WRITE cell_metadata.csv INCREMENTALLY
    #   Re-read each sample, shift coords, append to file
    # ─────────────────────────────────────────────
    log("=== STEP 4/8: Writing cell_metadata.csv (incremental) ===", t_start)
    
    out_dir = Path.cwd() / args.output_dir / "combined_output"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_meta_path = out_dir / "cell_metadata.csv"
    
    rows_written = 0
    
    for i, info in enumerate(sample_info):
        t_sample = time.perf_counter()
        log(f"  [{i + 1}/{n}] Processing {info['meta_path'].name}...", t_start)
    
        df = pd.read_csv(info['meta_path'])
    
        # Align rows to cell_by_gene order (inner join by cell ID)
        id_col = info['metadata_id_col']
        # Re-apply rename if original column was unnamed (fallback case)
        if id_col not in df.columns:
            df.rename(columns={df.columns[0]: id_col}, inplace=True)
        df[id_col] = df[id_col].astype(str)
        df = df.set_index(id_col).loc[info['valid_cell_ids']].reset_index()
    
        df["_source_file"] = str(info['meta_path'])
        df["_sample_id"] = info['sample_id']
    
        # Shift coordinates
        origin = grid_origins[i]
        old_bbox = info['bbox']
        dx = origin[0] - old_bbox[0]
        dy = origin[1] - old_bbox[1]
    
        df[x_col] += dx
        df[y_col] += dy
    
        for col, axis in bbox_cols.items():
            if col in df.columns:
                df[col] += dx if axis == "x" else dy
    
        # Append to CSV
        write_header = (i == 0)
        df.to_csv(out_meta_path, mode='a' if i > 0 else 'w', header=write_header, index=False)
        rows_written += len(df)
    
        del df
        gc.collect()
        log(f"    Appended {info['num_rows']:,} rows ({rows_written:,} total, {fmt_elapsed(time.perf_counter() - t_sample)})", t_start)
    
    log(f"Saved {out_meta_path} ({rows_written:,} rows)", t_start)
    
    
    # ─────────────────────────────────────────────
    # STEP 5: SPATIAL SANITY CHECK PLOT
    #   Re-read combined cell_metadata per sample for scatter
    # ─────────────────────────────────────────────
    log("=== STEP 5/8: Generating spatial layout plot ===", t_start)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    cmap = plt.colormaps.get_cmap("tab20")
    colors = [cmap(i / max(n - 1, 1)) for i in range(n)]
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    row_offset = 0
    for i, info in enumerate(sample_info):
        nrows = info['num_rows']
        chunk = pd.read_csv(out_meta_path, skiprows=range(1, row_offset + 1), nrows=nrows)
        ax.scatter(chunk[x_col], chunk[y_col], c=[colors[i]], s=0.5, alpha=0.4, rasterized=True)
        row_offset += nrows
        del chunk
        gc.collect()
    
    ax.set_title(f"Cell metadata spatial layout -- {rows_written:,} cells")
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig(out_dir / "check_spatial.png", dpi=150)
    plt.close()
    log(f"  Saved {out_dir / 'check_spatial.png'}", t_start)
    
    
    # ─────────────────────────────────────────────
    # STEP 6: WRITE cell_by_gene.csv INCREMENTALLY
    #   Read each sample's file and append directly
    # ─────────────────────────────────────────────
    log("=== STEP 6/8: Writing cell_by_gene.csv (incremental) ===", t_start)
    
    out_cbg_path = out_dir / "cell_by_gene.csv"
    cbg_rows_written = 0
    canonical_columns = None  # Column order from first sample — all others must match
    
    for i, info in enumerate(sample_info):
        t_sample = time.perf_counter()
        log(f"  [{i + 1}/{n}] Appending {info['cbg_path'].name}...", t_start)
    
        # Read and write in chunks to limit memory
        write_header = (i == 0)
        first_chunk = True
        sample_rows = 0
    
        valid_id_set = set(info['valid_cell_ids'])
        for chunk in pd.read_csv(info['cbg_path'], chunksize=50_000):
            # Filter to aligned cell IDs
            if info['cbg_id_col'] != 'cell':
                chunk = chunk.rename(columns={info['cbg_id_col']: 'cell'})
            chunk['cell'] = chunk['cell'].astype(str)
            chunk = chunk[chunk['cell'].isin(valid_id_set)]
    
            if len(chunk) == 0:
                continue
    
            if canonical_columns is None:
                canonical_columns = list(chunk.columns)
            else:
                # Reorder columns to match first sample's order
                # (different MERSCOPE slices may have genes in different column orders)
                chunk = chunk[canonical_columns]
    
            if write_header and first_chunk:
                chunk.to_csv(out_cbg_path, mode='w', header=True, index=False)
                first_chunk = False
            else:
                chunk.to_csv(out_cbg_path, mode='a', header=False, index=False)
            sample_rows += len(chunk)
    
        cbg_rows_written += sample_rows
        log(f"    Appended {sample_rows:,} rows ({cbg_rows_written:,} total, {fmt_elapsed(time.perf_counter() - t_sample)})", t_start)
    
        gc.collect()
    
    # Validate total
    if cbg_rows_written != rows_written:
        log(f"ERROR: Row count mismatch: metadata={rows_written:,}, cell_by_gene={cbg_rows_written:,}", t_start)
        sys.exit(1)
    
    log(f"Saved {out_cbg_path} ({cbg_rows_written:,} rows)", t_start)


# ─────────────────────────────────────────────
# STEP 6b: BORDER ARTIFACT DETECTION
#   Compute per-cell row sums, generate threshold
#   comparison plot, and write artifact mask CSVs
# ─────────────────────────────────────────────
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

log("=== STEP 6b: Border artifact threshold comparison ===", t_start)

# Compute per-cell row sums from combined cell_by_gene
cell_sums_chunks = []
n_genes = None
for chunk in pd.read_csv(out_cbg_path, chunksize=50_000):
    gene_cols = [c for c in chunk.columns if c != "cell"]
    if n_genes is None:
        n_genes = len(gene_cols)
    cell_sums_chunks.append(chunk[gene_cols].sum(axis=1).values)
cell_sums = np.concatenate(cell_sums_chunks)
del cell_sums_chunks

# Read cell ID column from cell_metadata.csv in chunks
meta_header = pd.read_csv(out_meta_path, nrows=0).columns.tolist()
meta_id_col = None
for candidate in ("EntityID", "id", "cell_id"):
    if candidate in meta_header:
        meta_id_col = candidate
        break
if meta_id_col is None:
    meta_id_col = meta_header[0]

cell_ids_chunks = []
for chunk in pd.read_csv(out_meta_path, usecols=[meta_id_col], dtype={meta_id_col: str}, chunksize=50_000):
    cell_ids_chunks.append(chunk[meta_id_col].values)
cell_ids = np.concatenate(cell_ids_chunks)
del cell_ids_chunks

log(f"  Computed row sums for {len(cell_sums):,} cells across {n_genes} genes", t_start)

# Read spatial coords for plotting
meta_all = pd.read_csv(out_meta_path, usecols=[x_col, y_col])
cx_all = meta_all[x_col].values
cy_all = meta_all[y_col].values
del meta_all

# Determine percentile thresholds to compare
compare_percentiles = list(range(5, 85, 5))
compare_thresholds = [(p, np.percentile(cell_sums, p)) for p in compare_percentiles]

# Generate threshold comparison grid
n_thresh = len(compare_thresholds)
n_plot_cols = min(n_thresh, 4)
n_plot_rows = math.ceil(n_thresh / n_plot_cols)

plt.style.use("dark_background")
fig, axes = plt.subplots(n_plot_rows, n_plot_cols, figsize=(7 * n_plot_cols, 7 * n_plot_rows),
                         squeeze=False)

for idx, (pctl, rs_thresh) in enumerate(compare_thresholds):
    ax = axes[idx // n_plot_cols][idx % n_plot_cols]
    mask = cell_sums > rs_thresh
    n_kept = mask.sum()
    n_removed = (~mask).sum()
    pct_kept = n_kept / len(cell_sums) * 100

    # Grey = removed (artifacts), lime = kept (real cells)
    ax.scatter(cx_all[~mask], cy_all[~mask], s=0.1, c='dimgrey', alpha=0.3, rasterized=True)
    ax.scatter(cx_all[mask], cy_all[mask], s=0.1, c='lime', alpha=0.4, rasterized=True)
    ax.set_title(f"percentile={pctl:.0f}  (row_sum={rs_thresh:.1f})\n"
                 f"kept={n_kept:,} ({pct_kept:.1f}%)  removed={n_removed:,}",
                 fontsize=10)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=6)

# Hide empty subplots
for idx in range(n_thresh, n_plot_rows * n_plot_cols):
    axes[idx // n_plot_cols][idx % n_plot_cols].set_visible(False)

fig.suptitle("Border artifact filtering — grey=removed, green=kept", fontsize=14, y=1.01)
plt.tight_layout()
thresh_plot_path = out_dir / "check_artifact_thresholds.png"
plt.savefig(thresh_plot_path, dpi=150, bbox_inches='tight')
plt.close()
log(f"  Saved {thresh_plot_path}", t_start)

# Write artifact mask CSVs for each percentile
for percentile in compare_percentiles:
    threshold = np.percentile(cell_sums, percentile)
    is_artifact = cell_sums <= threshold
    n_flagged = is_artifact.sum()
    pct_flagged = n_flagged / len(cell_sums) * 100

    mask_df = pd.DataFrame({"cell_id": cell_ids, "is_artifact": is_artifact})
    mask_path = out_dir / f"artifact_mask_p{percentile}.csv"
    mask_df.to_csv(mask_path, index=False)

    log(f"  p{percentile}: threshold={threshold:.4f}, flagged {n_flagged:,} / {len(cell_sums):,} ({pct_flagged:.1f}%)", t_start)
    del mask_df

del cell_sums, cell_ids, cx_all, cy_all
gc.collect()


if not args.mask_only:
    # ─────────────────────────────────────────────
    # STEP 7: PER-GENE SPATIAL SANITY CHECK PLOTS
    #   Plot expression of known marker genes on spatial coords
    # ─────────────────────────────────────────────
    log("=== STEP 7/8: Generating per-gene expression plots ===", t_start)

    MARKER_GENES = ["Slc17a7", "Gfap", "Gad2", "Drd1", "VIM", "KLHL1"]

    # Read available columns from combined cell_by_gene header
    cbg_all_cols = pd.read_csv(out_cbg_path, nrows=0).columns.tolist()
    genes_to_plot = [g for g in MARKER_GENES if g in cbg_all_cols]

    if not genes_to_plot:
        log(f"  No marker genes found in cell_by_gene columns. Checked: {MARKER_GENES}", t_start)
    else:
        log(f"  Found {len(genes_to_plot)} marker gene(s): {genes_to_plot}", t_start)

        # Read spatial coordinates from combined cell_metadata
        t_plot = time.perf_counter()
        log(f"  Reading spatial coordinates from cell_metadata...", t_start)
        meta_coords = pd.read_csv(out_meta_path, usecols=[x_col, y_col])
        cx = meta_coords[x_col].values
        cy = meta_coords[y_col].values
        del meta_coords
        gc.collect()

        for gene in genes_to_plot:
            t_gene = time.perf_counter()
            log(f"  Plotting {gene}...", t_start)

            # Read only this gene's column
            expr = pd.read_csv(out_cbg_path, usecols=[gene])[gene].values
            expr_norm = expr / expr.max() if expr.max() > 0 else expr
            nmax = np.percentile(expr_norm, 80)
            if nmax > 0:
                expr_norm = np.clip(expr_norm, 0, nmax) / nmax

            cmap = plt.cm.coolwarm(expr_norm)

            plt.style.use("dark_background")
            plt.figure(figsize=(20, 20))
            plt.scatter(cx, cy, s=1, c=cmap)
            plt.axis("equal")
            plt.title(f"{gene} expression -- {rows_written:,} cells")
            plt.tight_layout()
            plt.savefig(out_dir / f"check_gene_{gene}.png", dpi=150)
            plt.close()

            del expr, expr_norm, cmap
            gc.collect()
            log(f"    Saved check_gene_{gene}.png ({fmt_elapsed(time.perf_counter() - t_gene)})", t_start)

        del cx, cy
        gc.collect()
        log(f"  Per-gene plots done ({fmt_elapsed(time.perf_counter() - t_plot)})", t_start)


    # ─────────────────────────────────────────────
    # STEP 8: EXPRESSION SANITY CHECK PLOT
    #   Read combined cell_by_gene in chunks for histogram
    # ─────────────────────────────────────────────
    log("=== STEP 8/8: Generating expression distribution plot ===", t_start)

    all_gene_sums = []
    header_cols = None

    for chunk in pd.read_csv(out_cbg_path, chunksize=50_000):
        if header_cols is None:
            header_cols = [c for c in chunk.columns if c != "cell"]
        gene_sums = chunk[header_cols].sum(axis=1).values
        all_gene_sums.append(gene_sums)

    all_gene_sums = np.concatenate(all_gene_sums)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(all_gene_sums, bins=100, color="steelblue", edgecolor="none")
    ax.axvline(np.median(all_gene_sums), color="red", linestyle="--", label=f"median={np.median(all_gene_sums):.1f}")
    ax.set_title("Distribution of total gene counts per cell")
    ax.set_xlabel("total gene counts")
    ax.set_ylabel("# cells")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "check_expression.png", dpi=150)
    plt.close()

    del all_gene_sums
    gc.collect()
    log(f"  Saved {out_dir / 'check_expression.png'}", t_start)

# ─────────────────────────────────────────────
# DONE
# ─────────────────────────────────────────────
log(f"=== COMPLETE ===", t_start)
if not args.mask_only:
    log(f"  Samples: {n}", t_start)
    log(f"  Grid: {n_rows}x{n_cols}", t_start)
log(f"  Total rows: {rows_written:,}", t_start)
log(f"  Output dir: {out_dir}", t_start)
log(f"  Total time: {fmt_elapsed(time.perf_counter() - t_start)}", t_start)

