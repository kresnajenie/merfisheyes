# MERFISHeyes Benchmark Suite

Automated benchmarks for measuring load time, memory usage, and rendering performance across all supported data formats and sizes.

## Prerequisites

```bash
# Python dependencies
pip install playwright numpy matplotlib anndata scipy pandas

# Install Chromium for Playwright
python -m playwright install chromium

# Start the MERFISHeyes dev server (must be running for all browser benchmarks)
npm run dev
```

## Quick Start

Run a single cell scaling benchmark with synthetic data:

```bash
# 1. Generate test data (1K to 1M cells, 500 genes, all 3 formats)
python scripts/benchmark/generate_single_cell.py \
  --cells 1000 10000 50000 100000 300000 500000 1000000 \
  --genes 500

# 2. Run browser benchmarks
python scripts/benchmark/run_browser_benchmarks.py \
  --types single_cell --headless

# 3. Generate charts
python scripts/benchmark/generate_sc_scaling_chart.py
```

---

## Benchmark Scripts

### 1. Synthetic Data Generators

These create test datasets at configurable sizes.

#### `generate_single_cell.py` — H5AD, Xenium, MERSCOPE

```bash
# Generate all combinations (dry run to preview)
python scripts/benchmark/generate_single_cell.py --dry-run

# Generate specific sizes
python scripts/benchmark/generate_single_cell.py \
  --cells 1000 10000 100000 500000 1000000 \
  --genes 500 \
  --filetypes h5ad xenium merscope

# Limit max file size
python scripts/benchmark/generate_single_cell.py --max-file-size-gb 2
```

Output: `benchmark_data/synthetic/single_cell/{h5ad,xenium,merscope}/`

#### `generate_single_molecule.py` — Parquet, CSV

```bash
python scripts/benchmark/generate_single_molecule.py \
  --molecules 1000000 10000000 50000000 \
  --genes 1000 \
  --filetypes parquet csv
```

Output: `benchmark_data/synthetic/single_molecule/{parquet,csv}/`

### 2. Browser Benchmarks (Playwright)

These automate the browser to upload files and measure performance. **The dev server must be running.**

#### `run_browser_benchmarks.py` — Synthetic data, all formats

Loads synthetic datasets through the browser UI and measures load time, gene query time, memory, and FPS.

```bash
# Run all single cell benchmarks
python scripts/benchmark/run_browser_benchmarks.py \
  --types single_cell \
  --headless

# Filter by format and size
python scripts/benchmark/run_browser_benchmarks.py \
  --types single_cell \
  --sc-filetypes h5ad merscope \
  --max-cells 500000 \
  --headless

# Run single molecule benchmarks
python scripts/benchmark/run_browser_benchmarks.py \
  --types single_molecule \
  --sm-filetypes parquet \
  --headless

# Resume after interruption (skips already-tested datasets)
python scripts/benchmark/run_browser_benchmarks.py --resume --headless

# Custom output path
python scripts/benchmark/run_browser_benchmarks.py \
  --output benchmark_results/my_run.csv
```

Output: `benchmark_results/browser_benchmarks.csv` (or custom path)

**Options:**
| Flag | Description |
|------|-------------|
| `--headless` | Run browser without GUI |
| `--resume` | Skip datasets already in the output CSV |
| `--types` | `single_cell`, `single_molecule`, or both |
| `--sc-filetypes` | `h5ad`, `xenium`, `merscope` |
| `--sm-filetypes` | `parquet`, `csv` |
| `--max-cells` | Skip SC datasets above this cell count |
| `--max-molecules` | Skip SM datasets above this molecule count |
| `--max-genes` | Skip datasets above this gene count |
| `--timeout` | Seconds per test (default: 600) |
| `--url` | Dev server URL (default: `http://localhost:3000`) |

#### `run_real_data_benchmarks.py` — Real datasets from disk

Benchmarks real MERSCOPE/Xenium datasets from `/data/merfisheyes-test-data/test-data-sizes/`.

```bash
# Run everything
python scripts/benchmark/run_real_data_benchmarks.py --headless

# Single molecule only
python scripts/benchmark/run_real_data_benchmarks.py --only sm --headless

# Single cell MERSCOPE only
python scripts/benchmark/run_real_data_benchmarks.py --only sc-merscope --headless

# Single cell Xenium only
python scripts/benchmark/run_real_data_benchmarks.py --only sc-xenium --headless

# Skip files larger than 3 GB
python scripts/benchmark/run_real_data_benchmarks.py --max-size 3000 --headless
```

Output: `benchmark_results/real_data_benchmark_<timestamp>.json`

#### `run_s3_sm_benchmarks.py` — S3 lazy-loading

Benchmarks loading single molecule datasets from S3 via the "Load from S3" modal. Tests manifest download and initial rendering time.

```bash
python scripts/benchmark/run_s3_sm_benchmarks.py --headless
```

The S3 URLs are hardcoded in the script (BIL bucket). Edit `S3_DATASETS` in the script to add/change URLs.

**Requirements:**
- S3 bucket must be publicly readable
- CORS must allow requests from `http://localhost:3000`
- Each dataset folder needs `manifest.json.gz` and `genes/*.bin.gz`

Output: `benchmark_results/s3_sm_benchmark.csv`

### 3. Chart Generators

Generate charts from benchmark results. All charts are saved to `benchmark_results/charts/`.

```bash
# Single cell scaling chart (cells vs load time, by format)
python scripts/benchmark/generate_sc_scaling_chart.py

# S3 lazy-loading chart (molecules vs initial load time)
python scripts/benchmark/generate_s3_chart.py

# Real data dashboard (load time, memory, FPS summary)
python scripts/benchmark/generate_charts.py
```

### 4. Configuration

`config.py` contains all default parameters:
- Cell/molecule count ranges
- Gene count ranges
- File type lists
- Sparsity settings (adaptive based on matrix size)
- File size estimators
- Result schema (`BenchmarkResult` dataclass)

---

## Output Files

| File | Description |
|------|-------------|
| `benchmark_results/sc_load_benchmark.csv` | Synthetic single cell results |
| `benchmark_results/s3_sm_benchmark.csv` | S3 lazy-loading results |
| `benchmark_results/real_data_consolidated.json` | Real data results (all formats) |
| `benchmark_results/real_data_benchmark_*.json` | Individual real data runs |
| `benchmark_results/charts/sc_load_scaling.png` | SC load time vs cell count |
| `benchmark_results/charts/sc_filesize_scaling.png` | SC load time vs file size |
| `benchmark_results/charts/s3_sm_load_scaling.png` | S3 SM lazy-load scaling |
| `benchmark_results/charts/dashboard.png` | Real data summary dashboard |
| `benchmark_results/charts/load_times.png` | Real data load times |
| `benchmark_results/charts/memory_usage.png` | Real data memory |
| `benchmark_results/charts/fps.png` | Real data FPS |
| `benchmark_results/charts/sm_scaling.png` | SM molecule scaling |

---

## What Gets Measured

Each browser benchmark records:

| Metric | Description |
|--------|-------------|
| `load_time_s` | Time from file selection to canvas render |
| `gene_query_time_s` | Time to select and visualize a gene |
| `memory_peak_mb` | JS heap usage after loading |
| `fps_after_load` | Frames per second (2-second sample) |
| `ui_responsive` | Whether the UI remained interactive |
| `crashed` | Whether the browser tab crashed (OOM/timeout) |

---

## Typical Results (as of April 2026)

### Single Cell Local Load (500 genes)

| Cells | H5AD | Xenium | MERSCOPE |
|-------|------|--------|----------|
| 1K | 4.0s | 4.1s | 4.1s |
| 10K | 4.1s | 4.3s | 4.2s |
| 100K | 4.8s | 8.1s | 6.0s |
| 500K | 8.1s | 16.0s | 14.4s |
| 1M | 12.7s | -- | 16.1s |

### Single Molecule S3 Lazy-Loading

| Molecules | Load Time |
|-----------|-----------|
| 11.5M | 8.2s |
| 29M | 6.6s |
| 67M | 6.2s |
| 102M | 5.9s |
| 1.2B | 6.2s |

Load time is constant (~6s) because only the manifest is downloaded initially. Gene data is fetched on demand.
