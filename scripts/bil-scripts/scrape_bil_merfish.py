"""Scrape BIL MERFISH dataset listing into a CSV for bil-examples.json population.

Columns: title, description, species, tissue, platform, externalLink, bilCode
Skips: genes (per-dataset gene files loaded later) and entries (populated from S3).
"""
import csv
import html
import re
import sys
import urllib.request
from pathlib import Path

LIST_URL = "https://api.brainimagelibrary.org/web/list?fields=FullText&searchtxt=merfish"
VIEW_URL = "https://api.brainimagelibrary.org/web/view?bildid={bildid}"
OUTPUT_CSV = Path(__file__).parent / "bil-merfish-datasets.csv"

ROW_RE = re.compile(
    r'<a href="view\?bildid=(?P<bildid>[^"]+)">(?P=bildid)</a>.*?'
    r'<b>Title:\s*</b>\s*(?P<title>.*?)\s*<br\s*/?>'
    r'.*?<b>Abstract:\s*</b>\s*(?P<abstract>.*?)\s*<br\s*/?>'
    r'.*?<b>Species:\s*</b>\s*(?P<species>.*?)\s*<br\s*/?>',
    re.DOTALL,
)


def clean(text: str) -> str:
    text = html.unescape(text)
    text = re.sub(r"<[^>]+>", "", text)
    return re.sub(r"\s+", " ", text).strip()


def fetch(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="replace")


def scrape() -> list[dict]:
    page = fetch(LIST_URL)
    rows = []
    for match in ROW_RE.finditer(page):
        bildid = match.group("bildid").strip()
        rows.append({
            "title": clean(match.group("title")),
            "description": clean(match.group("abstract")),
            "species": clean(match.group("species")).title(),
            "tissue": "Brain",
            "platform": "MERFISH",
            "externalLink": VIEW_URL.format(bildid=bildid),
            "bilCode": bildid,
        })
    return rows


def main() -> int:
    rows = scrape()
    if not rows:
        print("No datasets parsed — HTML structure may have changed.", file=sys.stderr)
        return 1

    fieldnames = ["title", "description", "species", "tissue", "platform", "externalLink", "bilCode"]
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} datasets to {OUTPUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
