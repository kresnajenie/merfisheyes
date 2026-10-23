"""Populate `entries[]` in bil-examples.json from an S3 recursive listing.

Input: /tmp/bil_recursive.txt — output of
    aws s3 ls s3://merfisheyes-bil/bil-psc-data2/ --recursive

Rules (match the existing ace-dip-use example):
- Each bilCode with any `meyes_output/` key →
    {"label": "All Single Cell", "datasetType": "single_cell",
     "s3BaseUrl": ".../bil-psc-data2/{bilCode}/meyes_output/"}        (trailing /)
- Each (bilCode, sampleId) under `sm_output/{sampleId}/` →
    {"label": "{sampleId}", "datasetType": "single_molecule",
     "s3BaseUrl": ".../bil-psc-data2/{bilCode}/sm_output/{sampleId}"}  (no trailing /)

Only the canonical `meyes_output/` (not numbered variants) and canonical
single-slash paths are considered. test-* bilCodes are skipped.
"""
import json
import re
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
LISTING_PATH = Path("/tmp/bil_recursive.txt")
JSON_PATH = REPO_ROOT / "prisma" / "bil-examples.json"
BUCKET_BASE = "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-data2"

# Matches: bil-psc-data2/{bilCode}/meyes_output/...  (single slash, no numbered variants)
MEYES_RE = re.compile(r"bil-psc-data2/([^/]+)/meyes_output/")
# Matches: bil-psc-data2/{bilCode}/sm_output/{sampleId}/...
SM_RE = re.compile(r"bil-psc-data2/([^/]+)/sm_output/([^/]+)/")


def parse_listing() -> tuple[set[str], dict[str, set[str]]]:
    has_meyes: set[str] = set()
    sm_samples: dict[str, set[str]] = defaultdict(set)
    with LISTING_PATH.open() as f:
        for line in f:
            parts = line.rstrip().split(None, 3)
            if len(parts) < 4:
                continue
            key = parts[3]
            if "bil-psc-data2//" in key:  # skip legacy double-slash residuals
                continue
            m = MEYES_RE.search(key)
            if m:
                has_meyes.add(m.group(1))
                continue
            m = SM_RE.search(key)
            if m:
                bilcode, sample = m.group(1), m.group(2)
                sm_samples[bilcode].add(sample)
    return has_meyes, sm_samples


def build_entries(bilcode: str, has_meyes: set[str], sm: dict[str, set[str]]) -> list[dict]:
    entries: list[dict] = []
    if bilcode in has_meyes:
        entries.append({
            "label": "All Single Cell",
            "datasetType": "single_cell",
            "s3BaseUrl": f"{BUCKET_BASE}/{bilcode}/meyes_output/",
        })
    for sample in sorted(sm.get(bilcode, []), reverse=True):
        entries.append({
            "label": sample,
            "datasetType": "single_molecule",
            "s3BaseUrl": f"{BUCKET_BASE}/{bilcode}/sm_output/{sample}",
        })
    return entries


def main() -> None:
    has_meyes, sm_samples = parse_listing()
    # Drop test-* bilcodes
    has_meyes = {b for b in has_meyes if not b.startswith("test-")}
    sm_samples = {b: v for b, v in sm_samples.items() if not b.startswith("test-")}

    entries = json.loads(JSON_PATH.read_text())
    known_codes = {e.get("bilCode") for e in entries if e.get("bilCode")}

    populated = 0
    for entry in entries:
        code = entry.get("bilCode")
        built = build_entries(code, has_meyes, sm_samples) if code else []
        if built:
            entry["entries"] = built
            populated += 1
        # else: leave existing ["..."] placeholder untouched

    JSON_PATH.write_text(json.dumps(entries, indent=2) + "\n")

    s3_only = sorted((has_meyes | set(sm_samples)) - known_codes)
    print(f"Populated entries for {populated} / {len(entries)} datasets")
    if s3_only:
        print(f"S3 bilCodes not in JSON (skipped): {s3_only}")


if __name__ == "__main__":
    main()
