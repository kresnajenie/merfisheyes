#!/usr/bin/env python3
"""Symlink every single-molecule file in a databank into its single-molecule/ folder.

A Xenium export carries `transcripts.parquet` and a MERSCOPE export carries
`detected_transcripts.csv` — the same download is both a single-cell dataset and
a single-molecule dataset. This collects them (and any already nested under
single-molecule/) into one flat folder so the benchmark sees one entry per
dataset, without copying a byte.

    python3 scripts/bench/link-transcripts.py --root <databank> [--dry-run]

Idempotent: an existing correct link is left alone, a stale one is repointed.
Broken sources and unsupported types (.csv.gz — the processor takes only
.parquet and .csv) are reported and skipped.
"""
import argparse
import os
import sys
from pathlib import Path

# What each platform's export calls its molecule table.
XENIUM_NAMES = ["transcripts.parquet", "transcripts.csv"]
MERSCOPE_GLOB = "detected_transcripts*.csv"
# Real but unusable: the processor accepts only .parquet / .csv.
UNSUPPORTED = ["transcripts.csv.gz", "detected_transcripts.csv.gz"]
MAX_DEPTH = 3


def walk_dirs(base: Path, max_depth: int):
    """`base` and its subdirectories, up to max_depth levels down."""
    if not base.is_dir():
        return
    yield base
    if max_depth <= 0:
        return
    for child in sorted(base.iterdir()):
        if child.is_dir() and not child.name.startswith("."):
            yield from walk_dirs(child, max_depth - 1)


def clean(name: str) -> str:
    return "".join(c if (c.isalnum() or c in "-_.") else "-" for c in name).strip("-.")


def find_sources(root: Path):
    """[(platform, label, source_path)] for every molecule file in the databank."""
    found, skipped = [], []

    for platform, names in (("xenium", XENIUM_NAMES), ("merscope", None)):
        base = root / platform
        for d in walk_dirs(base, MAX_DEPTH):
            candidates = (
                [d / n for n in names] if names else sorted(d.glob(MERSCOPE_GLOB))
            )
            for src in candidates:
                if not src.exists():
                    if src.is_symlink():
                        skipped.append((str(src.relative_to(root)), "broken symlink"))
                    continue
                # Label = the dataset folder, plus the subfolder when nested,
                # plus the file stem when a folder holds several tables.
                rel_dir = d.relative_to(base)
                parts = [p for p in rel_dir.parts]
                stem = src.name.split(".")[0]
                if stem not in ("transcripts", "detected_transcripts"):
                    parts.append(stem)
                label = clean("-".join(parts)) or platform
                found.append((platform, label, src))

    # Files already sitting under single-molecule/, just nested in folders.
    sm = root / "single-molecule"
    for d in walk_dirs(sm, MAX_DEPTH):
        if d == sm:
            continue
        for src in sorted(list(d.glob("*.parquet")) + list(d.glob("*.csv"))):
            if src.is_symlink() and not src.exists():
                skipped.append((str(src.relative_to(root)), "broken symlink"))
                continue
            rel = d.relative_to(sm)
            platform = rel.parts[0] if rel.parts else "sm"
            parts = list(rel.parts[1:])
            stem = src.name.split(".")[0]
            if stem not in ("transcripts", "detected_transcripts"):
                parts.append(stem)
            label = clean("-".join(parts)) or clean(d.name)
            found.append((platform, label, src))

    for d in walk_dirs(root / "xenium", MAX_DEPTH):
        for bad in UNSUPPORTED:
            if (d / bad).exists():
                skipped.append((str((d / bad).relative_to(root)), "unsupported (.gz)"))

    return found, skipped


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f}{unit}" if unit == "B" else f"{n:.1f}{unit}"
        n /= 1024.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True, help="databank root (holds h5ad/ xenium/ merscope/ single-molecule/)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = args.root.resolve()
    sm_dir = root / "single-molecule"
    if not root.is_dir():
        print(f"No such databank: {root}", file=sys.stderr)
        return 2
    sm_dir.mkdir(exist_ok=True)

    found, skipped = find_sources(root)
    created = relinked = kept = 0

    for platform, label, src in sorted(found, key=lambda t: (t[0], t[1])):
        ext = ".parquet" if src.suffix.lower() == ".parquet" else ".csv"
        link = sm_dir / f"{platform}-{label}{ext}"
        target = os.path.realpath(src)

        if link.is_symlink() or link.exists():
            if link.is_symlink() and os.path.realpath(link) == target:
                kept += 1
                continue
            if link.is_symlink():
                action, relinked = "relink", relinked + 1
                if not args.dry_run:
                    link.unlink()
            else:
                print(f"  ! {link.name}: a real file is already there — leaving it")
                continue
        else:
            action, created = "link", created + 1

        size = human(os.path.getsize(target))
        print(f"  {action:<7} {link.name:<58} -> {src.relative_to(root)}  ({size})")
        if not args.dry_run:
            link.symlink_to(target)

    print(f"\n{created} created, {relinked} repointed, {kept} already correct"
          + ("  [dry run]" if args.dry_run else ""))
    if skipped:
        print("\nSkipped:")
        for path, why in skipped:
            print(f"  - {path}  ({why})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
