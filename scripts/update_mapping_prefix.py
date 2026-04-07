#!/usr/bin/env python3
"""
Update S3 prefix in an existing mapping.json file.

Usage:
    python update_mapping_prefix.py mapping.json --new-prefix https://new-bucket.s3.region.amazonaws.com/new_prefix

This replaces the S3 prefix in all link values while preserving the sample directory names.
For example:
    "https://old-bucket.s3.../old_prefix/sample_A"
  becomes:
    "https://new-bucket.s3.../new_prefix/sample_A"
"""

import argparse
import json
import sys
from pathlib import Path


def update_prefix(mapping_file: str, new_prefix: str):
    """
    Update the S3 prefix in all links of a mapping.json file.
    Preserves the last path segment (sample directory name) of each link.
    """
    mapping_path = Path(mapping_file)
    if not mapping_path.exists():
        print(f"Error: File not found: {mapping_file}", file=sys.stderr)
        sys.exit(1)

    with open(mapping_path, "r") as f:
        mapping = json.load(f)

    if "links" not in mapping:
        print("Error: No 'links' key found in mapping file.", file=sys.stderr)
        sys.exit(1)

    new_prefix = new_prefix.rstrip("/")

    old_links = mapping["links"]
    new_links = {}

    print(f"Updating {len(old_links)} link(s)...")
    for key, old_url in old_links.items():
        # Extract the last path segment (sample directory name)
        sample_name = old_url.rstrip("/").rsplit("/", 1)[-1]
        new_url = f"{new_prefix}/{sample_name}"
        new_links[key] = new_url
        print(f"  {key}: {old_url}")
        print(f"       → {new_url}")

    mapping["links"] = new_links

    with open(mapping_path, "w") as f:
        json.dump(mapping, f, indent=2)

    print(f"\nUpdated mapping written to: {mapping_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Update S3 prefix in a mapping.json file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "mapping_file",
        help="Path to the mapping.json file to update",
    )
    parser.add_argument(
        "--new-prefix",
        required=True,
        help="New S3 base URL prefix to replace the existing one",
    )

    args = parser.parse_args()
    update_prefix(args.mapping_file, args.new_prefix)


if __name__ == "__main__":
    main()
