"""The worker's reading of a molecule table: which columns it takes for gene,
x, y, z (worker/entrypoint.py), and gene names stored as bytes
(scripts/process_single_molecule.py)."""
from __future__ import annotations

import sys
import types
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "worker"))
sys.path.insert(0, str(ROOT / "scripts"))
# Only the S3 client needs boto3; the column logic does not
sys.modules.setdefault("boto3", types.ModuleType("boto3"))

import entrypoint  # noqa: E402
import process_single_molecule  # noqa: E402


def table(path: Path, columns: dict) -> Path:
    pq.write_table(pa.table(columns), path)
    return path


def test_platform_layouts_keep_their_presets(tmp_path):
    merscope = table(tmp_path / "m.parquet", {"gene": ["A"], "global_x": [1.0], "global_y": [1.0], "cell_id": [1]})
    xenium = table(tmp_path / "x.parquet", {"feature_name": ["A"], "x_location": [1.0], "y_location": [1.0]})

    assert entrypoint.detect_molecule_type(merscope)[0] == ["--dataset-type", "merscope", "--cell-id-col", "cell_id"]
    assert entrypoint.detect_molecule_type(xenium)[0] == ["--dataset-type", "xenium"]


def test_a_labs_own_columns_are_read_by_name(tmp_path):
    # Beaumont4GBM's merged parquet: Gene, X, Y, Z beside iGene, Cell ID, ...
    path = table(tmp_path / "s.parquet", {
        "iGene": [1], "X": [1.0], "Y": [2.0], "Z": [0.0], "Cell ID": [0],
        "Gene": ["GFAP"], "RNA_iset": [0], "iset": [2], "ifov": [0],
    })
    args, _ = entrypoint.detect_molecule_type(path)

    assert args == ["--dataset-type", "custom", "--gene-col", "Gene", "--x-col", "X", "--y-col", "Y", "--z-col", "Z"]


def test_csv_header_and_case(tmp_path):
    path = tmp_path / "t.csv"
    path.write_text("TARGET,x_um,y_um\nA,1,2\n")

    assert entrypoint.detect_molecule_type(path)[0] == [
        "--dataset-type", "custom", "--gene-col", "TARGET", "--x-col", "x_um", "--y-col", "y_um",
    ]


def test_no_molecule_columns_names_what_is_there(tmp_path):
    path = table(tmp_path / "n.parquet", {"name": ["A"], "u": [1.0], "v": [1.0]})

    with pytest.raises(RuntimeError, match="Columns in the file: name, u, v"):
        entrypoint.detect_molecule_type(path)


def test_the_browser_recognises_the_same_names():
    # lib/ingest/molecule-columns.ts mirrors MOLECULE_COLUMN_NAMES
    ts = (ROOT / "lib" / "ingest" / "molecule-columns.ts").read_text()
    for role, names in entrypoint.MOLECULE_COLUMN_NAMES.items():
        listed = ts.split(f"{role}: [", 1)[1].split("]", 1)[0]
        assert [n.strip().strip('"') for n in listed.split(",") if n.strip()] == names, role


def test_gene_names_stored_as_bytes_are_read_as_text(tmp_path):
    # Xenium Onboard Analysis 1.x writes feature_name as parquet BINARY
    path = table(tmp_path / "b.parquet", {
        "feature_name": pa.array([b"Gad1", b"NegControlProbe_00042", b"Gad1"], type=pa.binary()),
        "x_location": [1.0, 2.0, 3.0],
        "y_location": [1.0, 2.0, 3.0],
    })
    genes, coords, dims, _ = process_single_molecule.read_parquet_file(
        str(path), "feature_name", "x_location", "y_location", None
    )

    assert list(genes) == ["Gad1", "NegControlProbe_00042", "Gad1"]
    assert not process_single_molecule.should_filter_gene("Gad1")
    assert process_single_molecule.should_filter_gene("NegControlProbe_00042")
