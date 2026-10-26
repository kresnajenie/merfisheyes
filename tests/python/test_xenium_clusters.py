"""prepare_xenium.graph_clusters: 10x's own clustering of a Xenium export's
cells — the column the browser preview shows as "Cluster" — read from
analysis.tar.gz or an unpacked analysis/ folder, so it survives a save
(docs/plans/PREVIEW-SAVE.md, D5). A real export runs in
test_prepare_xenium.py (MERFISHEYES_XENIUM); these need none."""
from __future__ import annotations

import io
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import prepare_xenium  # noqa: E402

CSV = "Barcode,Cluster\naaafiiei-1,11\naaanbaof-1,10\naabdiein-1,2\n"


def test_reads_the_archive(tmp_path):
    with tarfile.open(tmp_path / "analysis.tar.gz", "w:gz") as archive:
        data = CSV.encode()
        info = tarfile.TarInfo(f"analysis/{prepare_xenium.GRAPHCLUST.removeprefix('analysis/')}")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
        # The other clusterings beside it are not the preview's column
        other = b"Barcode,Cluster\naaafiiei-1,1\n"
        info = tarfile.TarInfo("analysis/clustering/gene_expression_kmeans_2_clusters/clusters.csv")
        info.size = len(other)
        archive.addfile(info, io.BytesIO(other))

    column, labels = prepare_xenium.graph_clusters(tmp_path)

    assert column == "Cluster"
    assert labels.to_dict() == {"aaafiiei-1": "11", "aaanbaof-1": "10", "aabdiein-1": "2"}


def test_reads_an_unpacked_folder(tmp_path):
    path = tmp_path / prepare_xenium.GRAPHCLUST
    path.parent.mkdir(parents=True)
    path.write_text(CSV)

    column, labels = prepare_xenium.graph_clusters(tmp_path)

    assert column == "Cluster" and labels["aabdiein-1"] == "2"


def test_an_export_without_analysis_has_none(tmp_path):
    assert prepare_xenium.graph_clusters(tmp_path) is None
