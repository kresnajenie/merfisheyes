"""Cut the unit-test fixtures for the SpatialData-DB catalog sync from a release.

    python scripts/testdata/make-sddb-catalog-fixture.py <release-dir>

<release-dir> holds catalog.parquet and collections.parquet downloaded from a
`catalog-<date>` release of theislab/spatialdata-db-client. Writes the first
three catalog rows, schema and metadata intact, to
tests/data/sddb/catalog-fixture.parquet, and the first three collection rows
to tests/data/sddb/collections-fixture.parquet (a file in the *other*
published schema, which the reader must refuse).
"""
import sys
from pathlib import Path

import pyarrow.parquet as pq

src = Path(sys.argv[1])
out = Path(__file__).resolve().parents[2] / "tests" / "data" / "sddb"
out.mkdir(parents=True, exist_ok=True)

for name in ("catalog", "collections"):
    table = pq.read_table(src / f"{name}.parquet")
    pq.write_table(table.slice(0, 3), out / f"{name}-fixture.parquet", compression="snappy")
    print(name, "→", out / f"{name}-fixture.parquet", table.slice(0, 3).num_rows, "rows")
