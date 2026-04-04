"""Central configuration for the MERFISHeyes benchmark pipeline.

All benchmark parameters are customizable via CLI arguments in each script.
These defaults match the testing checklist specification.
"""

from dataclasses import dataclass, field, asdict
from datetime import datetime
from typing import Optional
import subprocess

# ---------------------------------------------------------------------------
# Single Cell Defaults
# ---------------------------------------------------------------------------
SINGLE_CELL_SIZES = [
    1_000, 10_000, 300_000, 500_000,
    1_000_000, 2_000_000, 3_000_000, 4_000_000,
    5_000_000, 6_000_000, 7_000_000, 10_000_000,
]
SINGLE_CELL_GENES = [100, 500, 1_000, 2_000, 5_000, 10_000, 15_000, 20_000]
SINGLE_CELL_FILETYPES = ["h5ad", "xenium", "merscope"]
SINGLE_CELL_OBS_COLUMNS = 20

# ---------------------------------------------------------------------------
# Single Molecule Defaults
# ---------------------------------------------------------------------------
SINGLE_MOLECULE_SIZES = [
    10_000_000, 30_000_000, 50_000_000, 100_000_000,
    200_000_000, 400_000_000, 500_000_000, 600_000_000,
    700_000_000, 800_000_000, 900_000_000, 1_000_000_000,
]
SINGLE_MOLECULE_GENES = [1_000, 2_000]
SINGLE_MOLECULE_FILETYPES = ["parquet", "csv"]
MAL_FRACTION = 0.5          # 50% of all molecules → "Mal"
ASSIGNED_FRACTION = 0.4     # Of the remaining 50%, 40% assigned, 60% unassigned

# ---------------------------------------------------------------------------
# Browser Benchmark
# ---------------------------------------------------------------------------
DEV_SERVER_URL = "http://localhost:3000"
BENCHMARK_TIMEOUT_S = 600   # 10 min max per test
GENE_TO_QUERY = "gene_0"    # Gene to measure query time (single cell)
SM_GENE_TO_QUERY = "Mal"    # Gene to measure query time (single molecule)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DEFAULT_SYNTH_DIR = "benchmark_data/synthetic"
DEFAULT_RESULTS_DIR = "benchmark_results"

# ---------------------------------------------------------------------------
# Max file size guard (GB) – skip generation if estimated size exceeds this
# ---------------------------------------------------------------------------
MAX_FILE_SIZE_GB = 10.0

# ---------------------------------------------------------------------------
# Cluster / cell-type label pools
# ---------------------------------------------------------------------------
CELL_TYPES = [
    "Neuron", "Astrocyte", "Oligodendrocyte", "Microglia", "Endothelial",
    "Pericyte", "OPC", "Ependymal", "Fibroblast", "Macrophage",
    "T_Cell", "B_Cell", "NK_Cell", "Monocyte", "Schwann",
]
LEIDEN_CLUSTERS = [f"cluster_{i}" for i in range(20)]

# Control-probe names used for unassigned molecules
CONTROL_PROBES = [
    "NegControlProbe_00001", "NegControlProbe_00002", "NegControlProbe_00003",
    "NegControlProbe_00004", "NegControlProbe_00005", "NegControlCodeword_0001",
    "NegControlCodeword_0002", "NegControlCodeword_0003", "antisense_00001",
    "antisense_00002", "BLANK_0001", "BLANK_0002", "BLANK_0003", "BLANK_0004",
    "BLANK_0005", "DeprecatedCodeword_0001", "DeprecatedCodeword_0002",
    "Unassigned", "unassigned", "UnassignedCodeword_0001",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def get_git_version() -> str:
    """Return the short git commit hash of the current repo."""
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        return "unknown"


def get_size_tier(n_items: int, data_type: str) -> str:
    if data_type == "single_cell":
        if n_items <= 10_000:
            return "tiny"
        if n_items <= 100_000:
            return "small"
        if n_items <= 1_000_000:
            return "medium"
        if n_items <= 5_000_000:
            return "large"
        return "xl"
    else:
        if n_items <= 10_000_000:
            return "tiny"
        if n_items <= 50_000_000:
            return "small"
        if n_items <= 200_000_000:
            return "medium"
        if n_items <= 500_000_000:
            return "large"
        return "xl"


def get_sparsity(n_cells: int, n_genes: int) -> float:
    """Adaptive sparsity to keep memory / file size practical."""
    total = n_cells * n_genes
    if total < 1_000_000:
        return 0.10
    if total < 100_000_000:
        return 0.05
    if total < 1_000_000_000:
        return 0.01
    if total < 10_000_000_000:
        return 0.005
    return 0.001


# ---------------------------------------------------------------------------
# File-size estimators (MB)
# ---------------------------------------------------------------------------
def estimate_h5ad_size_mb(n_cells: int, n_genes: int) -> float:
    density = get_sparsity(n_cells, n_genes)
    nnz = int(n_cells * n_genes * density)
    sparse_bytes = nnz * 8 + n_cells * 4
    obs_bytes = n_cells * 100
    obsm_bytes = n_cells * 16
    return (sparse_bytes + obs_bytes + obsm_bytes) * 0.5 / (1024 * 1024)


def estimate_xenium_size_mb(n_cells: int, n_genes: int) -> float:
    cells_csv = n_cells * 250          # id + coords + 20 cluster cols
    features_tsv = n_genes * 20
    density = get_sparsity(n_cells, n_genes)
    nnz = int(n_cells * n_genes * density)
    matrix_mtx = nnz * 20             # text "row col val\n"
    return (cells_csv + features_tsv + matrix_mtx) / (1024 * 1024)


def estimate_merscope_size_mb(n_cells: int, n_genes: int) -> float:
    metadata_csv = n_cells * 250
    categories_csv = n_cells * 50
    density = get_sparsity(n_cells, n_genes)
    # Wide-format cell_by_gene.csv
    header = n_genes * 15
    avg_row = 10 + n_genes * 2 + int(n_genes * density) * 5
    data_csv = n_cells * avg_row
    return (metadata_csv + categories_csv + header + data_csv) / (1024 * 1024)


def estimate_parquet_size_mb(n_molecules: int) -> float:
    return n_molecules * 7 / (1024 * 1024)


def estimate_csv_size_mb(n_molecules: int) -> float:
    return n_molecules * 40 / (1024 * 1024)


ESTIMATORS = {
    "h5ad": estimate_h5ad_size_mb,
    "xenium": estimate_xenium_size_mb,
    "merscope": estimate_merscope_size_mb,
    "parquet": lambda n, **_: estimate_parquet_size_mb(n),
    "csv": lambda n, **_: estimate_csv_size_mb(n),
}


# ---------------------------------------------------------------------------
# Result schema
# ---------------------------------------------------------------------------
@dataclass
class BenchmarkResult:
    """One row of benchmark output."""
    platform: str = "MERFISHeyes"
    data_type: str = ""           # single_cell | single_molecule
    size_tier: str = ""
    n_cells: Optional[int] = None
    n_molecules: Optional[int] = None
    n_genes: int = 0
    path: str = ""
    file_size_mb: float = 0.0
    filetype: str = ""

    load_time_s: Optional[float] = None
    gene_query_time_s: Optional[float] = None
    memory_peak_mb: Optional[float] = None
    ui_responsive: Optional[bool] = None
    fps_after_load: Optional[float] = None
    crashed: bool = False
    error: Optional[str] = None
    notes: str = ""
    timestamp: str = field(
        default_factory=lambda: datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    version: str = field(default_factory=get_git_version)

    def to_dict(self) -> dict:
        return asdict(self)

    # CSV column order matching the user's spec
    CSV_COLUMNS = [
        "platform", "data_type", "size_tier",
        "n_cells", "n_molecules", "n_genes",
        "path", "file_size_mb", "filetype",
        "load_time_s", "gene_query_time_s", "memory_peak_mb",
        "ui_responsive", "fps_after_load",
        "crashed", "error", "notes",
        "timestamp", "version",
    ]
