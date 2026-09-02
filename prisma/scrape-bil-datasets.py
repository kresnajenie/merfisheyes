"""
Discover and fill in BIL catalog datasets for prisma/bil-examples.json.

Two things this does, in one pass over https://download.brainimagelibrary.org/meyes/:

1. New datasets: bilCodes present on the BIL server but not yet in bil-examples.json.
   Full record is built from scratch (title/description/species/tissue/platform/
   metadata from the BIL API, genes + entries scraped from the server).

2. Backfill: existing bil-examples.json rows whose "genes" and/or "entries" are the
   literal placeholder ["..."] instead of real data. Only those two fields are
   touched — title/description/species/tissue/platform/metadata/thumbnailUrl are
   left exactly as curated, since only genes/entries were ever stubbed out.

For both, dataset "shape" (single cell vs single molecule, and gene list) comes
from the live folder layout under /meyes/{bilCode}/:
  - meyes_output/          -> single_cell entry, genes from meyes_output/expr/index.json
  - sm_output/genes/...    -> single_molecule entry (flat), genes from *.bin.gz filenames
                              (each gene has a base file and a "_uuuuuuuuuu" duplicate —
                              dedup to the base name)
  - sm_output/<slice-id>/genes/... -> one single_molecule entry per slice folder, gene
                              list is the union across slices
  - anything else (only combined_output/ and/or mmc_output/, no meyes_output/sm_output)
                              -> pipeline hasn't reached the chunking step yet, nothing to
                              scrape; the row is left alone and reported as skipped

A bilCode only gets a full record if the BIL API (api.brainimagelibrary.org) has a
submission for it — that's what excludes internal test/scratch folders (test-*,
mapmycells-reference, *-combine-v1, ...) without needing a hardcoded denylist.

Usage:
    python prisma/scrape-bil-datasets.py            # dry run — prints a summary
    python prisma/scrape-bil-datasets.py --apply    # writes prisma/bil-examples.json
"""

import json
import re
import sys
import time
import urllib.error
import urllib.request

DOWNLOAD_BASE = "https://download.brainimagelibrary.org/meyes"
BIL_API = "https://api.brainimagelibrary.org/retrieve?bildid="
S3_BASE = "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-data2"
INPUT_FILE = "prisma/bil-examples.json"
OUTPUT_FILE = "prisma/bil-examples.json"
REQUEST_DELAY = 0.15


def get(url: str, timeout: int = 20):
    """GET url, return (status, body_text). Never raises."""
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json, text/html"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, resp.read().decode(errors="replace")
    except urllib.error.HTTPError as e:
        return e.code, ""
    except Exception as e:
        return None, str(e)


def head_ok(url: str, timeout: int = 10) -> bool:
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status == 200
    except Exception:
        return False


def list_folders(url: str):
    status, body = get(url)
    if status != 200:
        return None
    names = re.findall(r'<a href="([^"]+)/">', body)
    return [n for n in names if n != ".."]


def list_files(url: str, suffix: str):
    status, body = get(url)
    if status != 200:
        return []
    return re.findall(r'<a href="([^"]+' + re.escape(suffix) + r')">', body)


def fetch_bil_record(bil_code: str):
    """Full BIL API submission record, or None if unregistered."""
    status, body = get(f"{BIL_API}{bil_code}")
    if status != 200:
        return None
    try:
        raw = json.loads(body)
    except json.JSONDecodeError:
        return None
    retjson = raw.get("retjson", [])
    return retjson[0] if retjson else None


def build_title_fields(record: dict):
    ds = (record.get("Dataset") or [{}])[0]
    spec = (record.get("Specimen") or [{}])[0]
    technique = ds.get("technique", "")
    platform = ds.get("other") if technique == "other" and ds.get("other") else technique
    return {
        "title": ds.get("title", ""),
        "description": ds.get("abstract", ""),
        "species": spec.get("species", ""),
        "tissue": spec.get("organname", ""),
        "platform": platform,
    }


def build_metadata(record: dict) -> dict:
    """Same extraction as scrape-bil-metadata.py, kept in sync intentionally."""
    metadata = {}

    submission = record.get("Submission", {})
    if submission.get("bildate"):
        metadata["publicationYear"] = submission["bildate"][:4]

    contributors = record.get("Contributors", [])
    personal = [c for c in contributors if c.get("nametype") != "Organizational"]
    if personal:
        leaders = [c for c in personal if c.get("contributortype") == "ProjectLeader"]
        if leaders:
            metadata["investigator"] = leaders[0].get("contributorname", "")
            metadata["institution"] = leaders[0].get("affiliation", "")
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

    funders = record.get("Funders", [])
    if funders:
        funder = funders[0]
        parts = [p for p in (funder.get("fundername"), funder.get("award_number")) if p]
        if parts:
            metadata["funding"] = " — ".join(parts)

    specimens = record.get("Specimen", [])
    if specimens:
        spec = specimens[0]
        age, age_unit = spec.get("age", ""), spec.get("ageunit", "")
        if age and age != "Unknown":
            metadata["age"] = f"{age} {age_unit}".strip()
        if spec.get("sex", "") not in ("", "Unknown"):
            metadata["sex"] = spec["sex"]
        if spec.get("genotype", "") not in ("", "Unknown"):
            metadata["genotype"] = spec["genotype"]

    datasets = record.get("Dataset", [])
    if datasets:
        ds = datasets[0]
        rights = ds.get("rightsidentifier", "") or ds.get("rights", "")
        if rights:
            metadata["license"] = rights
        if ds.get("technique"):
            metadata["technique"] = ds["technique"]

    return metadata


def scrape_genes_single_cell(bil_code: str):
    status, body = get(f"{DOWNLOAD_BASE}/{bil_code}/meyes_output/expr/index.json")
    if status != 200:
        return None
    try:
        data = json.loads(body)
    except json.JSONDecodeError:
        return None
    return [g["name"] for g in data.get("genes", [])]


def scrape_genes_from_dir(url: str):
    files = list_files(url, ".bin.gz")
    genes = set()
    for f in files:
        name = f[: -len(".bin.gz")]
        name = re.sub(r"_uuuuuuuuuu$", "", name)
        genes.add(name)
    return sorted(genes)


def thumbnail_or_none(url: str):
    return url if head_ok(url) else None


def build_record(bil_code: str, structure: list):
    """Returns (genes, entries) for a bilCode given its top-level folder listing,
    or (None, None) if there's nothing scrapable yet."""
    if "meyes_output" in structure:
        genes = scrape_genes_single_cell(bil_code)
        if genes is None:
            return None, None
        thumb_url = f"{S3_BASE}/{bil_code}/meyes_output/thumbnail.jpg"
        entries = [
            {
                "label": "All Single Cell",
                "datasetType": "single_cell",
                "s3BaseUrl": f"{S3_BASE}/{bil_code}/meyes_output/",
                "thumbnailUrl": thumbnail_or_none(thumb_url),
            }
        ]
        return genes, entries

    if "sm_output" in structure:
        sub = list_folders(f"{DOWNLOAD_BASE}/{bil_code}/sm_output/")
        if sub is None:
            return None, None

        if sub == ["genes"] or (len(sub) == 1 and sub[0] == "genes"):
            genes = scrape_genes_from_dir(f"{DOWNLOAD_BASE}/{bil_code}/sm_output/genes/")
            if not genes:
                return None, None
            thumb_url = f"{S3_BASE}/{bil_code}/sm_output/genes/thumbnail.jpg"
            entries = [
                {
                    "label": "genes",
                    "datasetType": "single_molecule",
                    "s3BaseUrl": f"{S3_BASE}/{bil_code}/sm_output",
                    "thumbnailUrl": thumbnail_or_none(thumb_url),
                }
            ]
            return genes, entries

        # Multi-slice: one entry per slice sub-folder, gene list is the union.
        all_genes = set()
        entries = []
        for slice_id in sub:
            slice_genes = scrape_genes_from_dir(
                f"{DOWNLOAD_BASE}/{bil_code}/sm_output/{slice_id}/genes/"
            )
            all_genes.update(slice_genes)
            thumb_url = f"{S3_BASE}/{bil_code}/sm_output/{slice_id}/thumbnail.jpg"
            entries.append(
                {
                    "label": slice_id,
                    "datasetType": "single_molecule",
                    "s3BaseUrl": f"{S3_BASE}/{bil_code}/sm_output/{slice_id}",
                    "thumbnailUrl": thumbnail_or_none(thumb_url),
                }
            )
        if not all_genes:
            return None, None
        return sorted(all_genes), entries

    return None, None


def main():
    apply = "--apply" in sys.argv

    with open(INPUT_FILE) as f:
        datasets = json.load(f)
    existing_codes = {d["bilCode"] for d in datasets}

    print("Listing https://download.brainimagelibrary.org/meyes/ ...")
    server_folders = list_folders(f"{DOWNLOAD_BASE}/")
    if server_folders is None:
        print("ERROR: could not list the BIL download server.")
        sys.exit(1)
    server_folders = [f for f in server_folders if not f.startswith("_")]

    new_codes = sorted(set(server_folders) - existing_codes)

    placeholder_gene_codes = {d["bilCode"] for d in datasets if d.get("genes") == ["..."]}
    placeholder_entry_codes = {
        d["bilCode"]
        for d in datasets
        if any(not isinstance(e, dict) for e in d.get("entries", []))
    }
    backfill_codes = sorted(placeholder_gene_codes | placeholder_entry_codes)

    print(f"Candidate new bilCodes on server: {len(new_codes)}")
    print(f"Existing rows needing backfill (placeholder genes/entries): {len(backfill_codes)}")
    print()

    added, backfilled, skipped_no_bil, skipped_no_data = [], [], [], []

    # --- New datasets ---
    for i, code in enumerate(new_codes):
        print(f"[new {i + 1}/{len(new_codes)}] {code} ...", end=" ", flush=True)
        record = fetch_bil_record(code)
        if record is None:
            print("SKIP (no BIL API record)")
            skipped_no_bil.append(code)
            time.sleep(REQUEST_DELAY)
            continue

        structure = list_folders(f"{DOWNLOAD_BASE}/{code}/") or []
        genes, entries = build_record(code, structure)
        if genes is None:
            print(f"SKIP (no meyes_output/sm_output yet — has {structure})")
            skipped_no_data.append(code)
            time.sleep(REQUEST_DELAY)
            continue

        fields = build_title_fields(record)
        new_row = {
            "title": fields["title"],
            "description": fields["description"],
            "species": fields["species"],
            "tissue": fields["tissue"],
            "platform": fields["platform"],
            "genes": genes,
            "externalLink": f"https://api.brainimagelibrary.org/web/view?bildid={code}",
            "bilCode": code,
            "entries": entries,
            "metadata": build_metadata(record),
            "thumbnailUrl": None,
        }
        datasets.append(new_row)
        added.append(code)
        print(f"OK — {len(genes)} genes, {len(entries)} entr{'y' if len(entries)==1 else 'ies'}")
        time.sleep(REQUEST_DELAY)

    print()

    # --- Backfill existing placeholder rows ---
    by_code = {d["bilCode"]: d for d in datasets}
    for i, code in enumerate(backfill_codes):
        print(f"[backfill {i + 1}/{len(backfill_codes)}] {code} ...", end=" ", flush=True)
        structure = list_folders(f"{DOWNLOAD_BASE}/{code}/")
        if structure is None:
            print("SKIP (not found on server)")
            skipped_no_data.append(code)
            time.sleep(REQUEST_DELAY)
            continue

        genes, entries = build_record(code, structure)
        if genes is None:
            print(f"SKIP (no meyes_output/sm_output yet — has {structure})")
            skipped_no_data.append(code)
            time.sleep(REQUEST_DELAY)
            continue

        row = by_code[code]
        touched = []
        if row.get("genes") == ["..."]:
            row["genes"] = genes
            touched.append("genes")
        if any(not isinstance(e, dict) for e in row.get("entries", [])):
            row["entries"] = entries
            touched.append("entries")
        backfilled.append(code)
        print(f"OK — updated {', '.join(touched)} ({len(genes)} genes, {len(entries)} entries)")
        time.sleep(REQUEST_DELAY)

    print()
    print(f"New rows added:        {len(added)}")
    print(f"Existing rows backfilled: {len(backfilled)}")
    print(f"Skipped — no BIL record:  {len(skipped_no_bil)}  {skipped_no_bil}")
    print(f"Skipped — no data yet:    {len(skipped_no_data)}  {skipped_no_data}")

    if not apply:
        print("\nDry run — pass --apply to write prisma/bil-examples.json")
        return

    with open(OUTPUT_FILE, "w") as f:
        json.dump(datasets, f, indent=2)
    print(f"\nWrote {OUTPUT_FILE} ({len(datasets)} total datasets).")


if __name__ == "__main__":
    main()
