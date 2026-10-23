#!/usr/bin/env bash
# Find the pipeline scripts that are referenced by the docs but exist in no
# branch of the repository (docs/plans/BIL-CLUSTER.md §1.3, "Tier B").
#
#   align2.py       per-set translation of molecules onto cells
#   build_store.py  the 32-worker parallel spill that built the 11.24 B-molecule
#                   12_16_2025_all_sets store (docs/molecule-pyramid.md:281)
#
# build_store.py matters most: it is the parallel-spill wrapper that
# scripts/transcripts_to_points.py lacks (BIL-CLUSTER.md §1.4). Finding it may
# remove most of Step 4.
#
# Run it wherever the work was done — the BIL cluster over SSH, a laptop, a
# scratch volume:
#
#   bash scripts/bil-scripts/find_missing_scripts.sh
#   bash scripts/bil-scripts/find_missing_scripts.sh /path/to/search ...
#
# Searches by file name and, because the files may have been renamed, by the
# identifiers only they would contain.

set -uo pipefail

ROOTS=("$@")
if [ ${#ROOTS[@]} -eq 0 ]; then
    ROOTS=("$HOME" /bil/users /bil/data/meyes /scratch /tmp)
fi

NAMES=(align2.py build_store.py)
# Identifiers distinctive enough that a hit is almost certainly the real file.
NEEDLES=(PointsWriter 'spill(' feature_major tile_major band_um)

echo "Searching for the missing pipeline scripts"
echo "roots: ${ROOTS[*]}"
echo

echo "== by file name =="
for root in "${ROOTS[@]}"; do
    [ -d "$root" ] || continue
    for n in "${NAMES[@]}"; do
        find "$root" -maxdepth 8 -type f -name "$n" 2>/dev/null
    done
done | sort -u | sed 's/^/  /' || true
echo "  (nothing above = not found by name)"
echo

echo "== by content: a .py mentioning PointsWriter outside this repo =="
for root in "${ROOTS[@]}"; do
    [ -d "$root" ] || continue
    grep -rl --include='*.py' -e 'PointsWriter' "$root" 2>/dev/null
done | grep -v '/merfisheyes/scripts/points_zarr\.py$' \
     | grep -v '/merfisheyes/scripts/\(transcripts_to_points\|legacy_to_points\)\.py$' \
     | sort -u | sed 's/^/  /' || true
echo "  (nothing above = no other caller of the writer on this machine)"
echo

echo "== by content: a .py that spills in parallel =="
for root in "${ROOTS[@]}"; do
    [ -d "$root" ] || continue
    grep -rl --include='*.py' -e 'ProcessPoolExecutor' -e 'multiprocessing' "$root" 2>/dev/null
done | sort -u | while read -r f; do
    grep -lq -e 'spill' -e 'feature_major' -e 'tile_major' "$f" 2>/dev/null && echo "  $f"
done || true
echo

echo "== notebooks and archives worth opening by hand =="
for root in "${ROOTS[@]}"; do
    [ -d "$root" ] || continue
    find "$root" -maxdepth 6 -type f \
        \( -name '*.ipynb' -o -name '*align*' -o -name '*build_store*' \) \
        2>/dev/null | grep -iE 'align|store|point|molecule' | head -40
done | sort -u | sed 's/^/  /' || true

echo
echo "Found something? Put it in scripts/ and record it in docs/plans/BIL-CLUSTER.md §1.3."
echo "Found nothing anywhere? Say so in §1.3 — Step 4 then rebuilds the wrapper."
