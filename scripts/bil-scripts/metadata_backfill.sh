#!/bin/bash
# Give stores that are already in S3 their metadata, without reprocessing them
# (docs/plans/METADATA.md step 5; docs/specs/metadata.md §4–5).
#
#   metadata_backfill.sh prepare <folder>...   build + bionty check + validate, locally
#   metadata_backfill.sh publish <folder>...   then upload each changed root zarr.json and register
#
# A <folder> holds one store and its sidecars (stats.json, qc.json, catalog.json,
# mapping.json), e.g. ~/data/bil-data/processed/<code>/store; the BIL code is the name of
# its parent folder. Per store:
#
#   1. build     scripts/metadata/build_metadata.py (meyes-prepare), with this code's entry
#                of metadata-overrides.json when it has one
#   2. check     scripts/metadata/validate_metadata.py --write (meyes-lamin: bionty)
#   3. validate  scripts/validate_store.py --require-metadata
#   4. upload    aws s3 cp of the root zarr.json only — the sync is --size-only and would
#                skip a same-size rewrite (METADATA.md, handoff); skipped when S3 already
#                has these bytes
#   5. register  scripts/bil-scripts/register_lamin.py, all passing stores in one call
#
# Only the root zarr.json of each store changes. Steps 4–5 run with `publish` only, and only
# for stores that passed 1–3. Logs and a summary go to $META_RUNS/<timestamp>/.
# Env: MEYES_PREPARE_PY, MEYES_LAMIN_PY (interpreters), META_JOBS (parallel builds, default 6).

set -uo pipefail
MODE="${1:-}"; shift || true
[[ "$MODE" == prepare || "$MODE" == publish ]] && [ $# -gt 0 ] || { sed -n '2,25p' "$0"; exit 2; }

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PREP_PY="${MEYES_PREPARE_PY:-$HOME/miniconda3/envs/meyes-prepare/bin/python}"
LAMIN_PY="${MEYES_LAMIN_PY:-$HOME/miniconda3/envs/meyes-lamin/bin/python}"
BUCKET="merfisheyes-bil"; PREFIX="bil-psc-zarr"
OVERRIDES="$REPO/scripts/bil-scripts/metadata-overrides.json"
RUN="${META_RUNS:-$HOME/data/bil-data/metadata-runs}/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$RUN"
echo "run: $RUN  mode: $MODE  stores: $#"

one () {   # prepare one folder; prints "<code>\t<status>\t<note>"
    local dir="${1%/}" code store name extra=""
    code="$(basename "$(dirname "$dir")")"
    store="$(ls -d "$dir"/*.zarr 2>/dev/null | head -1)"
    local log="$RUN/$code.log"
    [ -n "$store" ] || { printf '%s\tfail\tno store in %s\n' "$code" "$dir"; return; }
    if "$PREP_PY" -I -c "import json,sys; o=json.load(open('$OVERRIDES')); sys.exit(0 if '$code' in o else 1)"; then
        extra="$RUN/$code.metadata.json"
        "$PREP_PY" -I -c "import json; o=json.load(open('$OVERRIDES'))['$code']; json.dump({k: v for k, v in o.items() if not k.startswith('_')}, open('$extra', 'w'))"
    fi
    { "$PREP_PY" -I "$REPO/scripts/metadata/build_metadata.py" "$dir" ${extra:+--metadata "$extra"}; } >>"$log" 2>&1 \
        || { printf '%s\tfail\tbuild (see %s)\n' "$code" "$log"; return; }
    "$LAMIN_PY" -I "$REPO/scripts/metadata/validate_metadata.py" "$store" --write >>"$log" 2>&1
    local v=$?
    [ $v -eq 0 ] || { printf '%s\tfail\tbionty exit %s (see %s)\n' "$code" "$v" "$log"; return; }
    "$PREP_PY" -I "$REPO/scripts/validate_store.py" "$store" --require-metadata --sample 1000 >>"$log" 2>&1 \
        || { printf '%s\tfail\tvalidate_store (see %s)\n' "$code" "$log"; return; }
    local tech
    tech="$("$PREP_PY" -I -c "import json; m=json.load(open('$store/zarr.json'))['attributes']['meyes_metadata']; print(m.get('organism'), m.get('tissue'), m.get('technology'), m.get('development_stage'), sep=' / ')")"
    printf '%s\tpass\t%s\n' "$code" "$tech"
}
export -f one; export RUN PREP_PY LAMIN_PY REPO OVERRIDES

printf '%s\n' "$@" | xargs -P "${META_JOBS:-6}" -I{} bash -c 'one "$@"' _ {} | tee "$RUN/summary.tsv"
PASS=$(awk -F'\t' '$2=="pass"' "$RUN/summary.tsv" | wc -l)
echo "prepared: $PASS pass, $(( $# - PASS )) fail"
[ "$MODE" == publish ] || exit $(( $# - PASS > 0 ))

uploaded=(); codes=()
for dir in "$@"; do
    dir="${dir%/}"; code="$(basename "$(dirname "$dir")")"
    awk -F'\t' -v c="$code" '$1==c && $2=="pass"{f=1} END{exit !f}' "$RUN/summary.tsv" || continue
    store="$(ls -d "$dir"/*.zarr | head -1)"; key="$PREFIX/$code/$(basename "$store")/zarr.json"
    remote="$(curl -sf "https://$BUCKET.s3.us-west-2.amazonaws.com/$key" | md5sum | cut -c1-32)"
    if [ "$remote" != "$(md5sum < "$store/zarr.json" | cut -c1-32)" ]; then
        aws s3 cp "$store/zarr.json" "s3://$BUCKET/$key" --content-type application/json --only-show-errors \
            || { echo "$code: upload FAILED"; continue; }
        uploaded+=("$code")
    fi
    codes+=("$code")
done
echo "uploaded ${#uploaded[@]} root zarr.json: ${uploaded[*]:-none}"
[ ${#codes[@]} -gt 0 ] || exit 1
"$LAMIN_PY" -I "$REPO/scripts/bil-scripts/register_lamin.py" "${codes[@]}" 2>&1 | grep -v -i futurewarning | tee "$RUN/register.log"
exit "${PIPESTATUS[0]}"
