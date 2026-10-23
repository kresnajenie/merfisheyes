"""Populate prisma/bil-examples.json from the scraped BIL MERFISH CSV.

Preserves existing fully-populated entries (matched by bilCode) and appends
placeholder entries for the rest with genes/entries set to ["..."].
"""
import csv
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CSV_PATH = Path(__file__).parent / "bil-merfish-datasets.csv"
JSON_PATH = REPO_ROOT / "prisma" / "bil-examples.json"


def main() -> None:
    existing = json.loads(JSON_PATH.read_text())
    by_code = {entry.get("bilCode"): entry for entry in existing if entry.get("bilCode")}

    with CSV_PATH.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    out = []
    for row in rows:
        code = row["bilCode"]
        if code in by_code:
            out.append(by_code[code])
            continue
        out.append({
            "title": row["title"],
            "description": row["description"],
            "species": row["species"],
            "tissue": row["tissue"],
            "platform": row["platform"],
            "genes": ["..."],
            "externalLink": row["externalLink"],
            "bilCode": code,
            "entries": ["..."],
        })

    JSON_PATH.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Wrote {len(out)} entries to {JSON_PATH} (preserved {sum(1 for r in rows if r['bilCode'] in by_code)} existing)")


if __name__ == "__main__":
    main()
