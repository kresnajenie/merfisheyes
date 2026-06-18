"""
Scrape metadata from the BIL API for all datasets in bil-examples.json.
Populates the 'metadata' field with investigator, institution, funding, etc.

Usage:
    python prisma/scrape-bil-metadata.py
"""

import json
import time
import urllib.request
import urllib.error
import sys

BIL_API = "https://api.brainimagelibrary.org/retrieve?bildid="
INPUT_FILE = "prisma/bil-examples.json"
OUTPUT_FILE = "prisma/bil-examples.json"


def fetch_bil_metadata(bil_code: str) -> dict:
    """Fetch metadata from BIL API and extract relevant fields."""
    url = f"{BIL_API}{bil_code}"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            raw = json.loads(resp.read().decode())
    except (urllib.error.URLError, json.JSONDecodeError) as e:
        print(f"  ERROR fetching {bil_code}: {e}")
        return {}

    # Data is nested under retjson[0]
    retjson = raw.get("retjson", [])
    if not retjson:
        return {}
    data = retjson[0]

    metadata = {}

    # Submission date
    submission = data.get("Submission", {})
    if submission.get("bildate"):
        metadata["publicationYear"] = submission["bildate"][:4]

    # Contributors — filter out organizational entries, find project leaders
    contributors = data.get("Contributors", [])
    personal = [c for c in contributors if c.get("nametype") != "Organizational"]
    if personal:
        leaders = [c for c in personal if c.get("contributortype") == "ProjectLeader"]
        if leaders:
            metadata["investigator"] = leaders[0].get("contributorname", "")
            metadata["institution"] = leaders[0].get("affiliation", "")
            # Deduplicate co-investigators by name
            seen = {leaders[0].get("contributorname", "")}
            co = []
            for c in leaders[1:]:
                name = c.get("contributorname", "")
                if name and name not in seen:
                    co.append(name)
                    seen.add(name)
            if co:
                metadata["coInvestigators"] = co
        else:
            contacts = [c for c in personal if c.get("contributortype") == "ContactPerson"]
            pick = contacts[0] if contacts else personal[0]
            metadata["investigator"] = pick.get("contributorname", "")
            metadata["institution"] = pick.get("affiliation", "")

    # Funding
    funders = data.get("Funders", [])
    if funders:
        funder = funders[0]
        parts = []
        if funder.get("fundername"):
            parts.append(funder["fundername"])
        if funder.get("award_number"):
            parts.append(funder["award_number"])
        if parts:
            metadata["funding"] = " — ".join(parts)

    # Specimen
    specimens = data.get("Specimen", [])
    if specimens:
        spec = specimens[0]
        age = spec.get("age", "")
        age_unit = spec.get("ageunit", "")
        if age and age != "Unknown":
            metadata["age"] = f"{age} {age_unit}".strip()
        sex = spec.get("sex", "")
        if sex and sex != "Unknown":
            metadata["sex"] = sex
        genotype = spec.get("genotype", "")
        if genotype and genotype != "Unknown":
            metadata["genotype"] = genotype

    # Dataset info
    datasets = data.get("Dataset", [])
    if datasets:
        ds = datasets[0]
        rights = ds.get("rightsidentifier", "") or ds.get("rights", "")
        if rights:
            metadata["license"] = rights
        technique = ds.get("technique", "")
        if technique:
            metadata["technique"] = technique

    return metadata


def main():
    with open(INPUT_FILE) as f:
        datasets = json.load(f)

    total = len(datasets)
    print(f"Scraping metadata for {total} datasets...")

    for i, ds in enumerate(datasets):
        bil_code = ds.get("bilCode", "")
        if not bil_code:
            print(f"  [{i+1}/{total}] SKIP — no bilCode")
            continue

        print(f"  [{i+1}/{total}] {bil_code}...", end=" ", flush=True)
        meta = fetch_bil_metadata(bil_code)
        if meta:
            ds["metadata"] = meta
            print(f"OK — {meta.get('investigator', '?')}")
        else:
            print("EMPTY")

        # Be polite to the API
        time.sleep(0.3)

    with open(OUTPUT_FILE, "w") as f:
        json.dump(datasets, f, indent=2)

    print(f"\nDone. Updated {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
