"""Scrape BIL Xenium dataset listing into a CSV of (bilCode, bilPath).

Output columns: bilCode, bilPath
- bilCode: BIL data ID (e.g. ace-owl-get)
- bilPath: filesystem path on BIL (e.g. /bil/data/f2/59/.../CJ23.56.004.CX/1388176147)
"""
import csv
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

LIST_URL = "https://api.brainimagelibrary.org/web/list?fields=FullText&searchtxt=xenium"
VIEW_URL = "https://api.brainimagelibrary.org/web/view?bildid={bildid}"
OUTPUT_CSV = Path(__file__).parent / "bil-xenium-datasets.csv"

BILDID_RE = re.compile(r'view\?bildid=([a-z0-9-]+)')
PATH_RE = re.compile(r'(/bil/data/[A-Za-z0-9._/-]+?)(?=<|\s|$)')


def fetch(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="replace")


def get_bil_path(bildid: str) -> tuple[str, str]:
    try:
        page = fetch(VIEW_URL.format(bildid=bildid))
    except Exception as e:
        return bildid, f"ERROR: {e}"
    m = PATH_RE.search(page)
    return bildid, m.group(1).rstrip("/") if m else ""


def main() -> int:
    page = fetch(LIST_URL)
    bildids = sorted(set(BILDID_RE.findall(page)))
    if not bildids:
        print("No BIL IDs parsed — HTML structure may have changed.", file=sys.stderr)
        return 1

    print(f"Found {len(bildids)} Xenium datasets, fetching paths…", file=sys.stderr)

    results: list[tuple[str, str]] = []
    with ThreadPoolExecutor(max_workers=8) as pool:
        for bildid, path in pool.map(get_bil_path, bildids):
            results.append((bildid, path))
            print(f"  {bildid}\t{path}", file=sys.stderr)

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["bilCode", "bilPath"])
        writer.writerows(results)

    missing = sum(1 for _, p in results if not p or p.startswith("ERROR"))
    print(f"Wrote {len(results)} rows to {OUTPUT_CSV} ({missing} without path)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
