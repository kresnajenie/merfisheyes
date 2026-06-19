"""Merge per-dataset gene lists from bil-genes-long.tsv into bil-examples.json.

TSV format: one line per gene — `bilCode<TAB>geneName`.
For each bilCode found in the TSV, the matching entry's `genes` field is
replaced with the full ordered list. Entries without gene data keep their
existing value (usually the ["..."] placeholder).
"""
import csv
import json
from collections import OrderedDict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
TSV_PATH = Path(__file__).parent / "bil-genes-long.tsv"
JSON_PATH = REPO_ROOT / "prisma" / "bil-examples.json"


def load_genes() -> dict[str, list[str]]:
    genes: dict[str, list[str]] = OrderedDict()
    with TSV_PATH.open(newline="", encoding="utf-8") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 2:
                continue
            bilcode, gene = row[0].strip(), row[1].strip()
            if not bilcode or not gene:
                continue
            genes.setdefault(bilcode, []).append(gene)
    return genes


def main() -> None:
    genes_by_code = load_genes()
    entries = json.loads(JSON_PATH.read_text())
    known_codes = {e.get("bilCode") for e in entries if e.get("bilCode")}

    matched = 0
    for entry in entries:
        code = entry.get("bilCode")
        if code and code in genes_by_code:
            entry["genes"] = genes_by_code[code]
            matched += 1

    JSON_PATH.write_text(json.dumps(entries, indent=2) + "\n")

    unmatched = sorted(set(genes_by_code) - known_codes)
    print(f"Updated {matched} entries with gene lists ({len(entries)} total).")
    if unmatched:
        print(f"Skipped {len(unmatched)} TSV bilCodes not found in JSON: {unmatched}")


if __name__ == "__main__":
    main()
