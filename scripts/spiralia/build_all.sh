#!/usr/bin/env bash
# Build chunked labelled-molecule datasets for every analysed spiralia embryo.
#
#   ./build_all.sh                 # all embryos
#   ./build_all.sh MER1_E00 ...    # just these
#
# Each embryo: dill object -> parquet -> chunked folder. The parquet is an
# intermediate and is deleted once the chunked folder exists.
set -uo pipefail

ROOT=/home/data/yiqun-spiralia/Sep2026
OUT=$ROOT/merfisheyes_export
PY=/home/kjenie/miniconda3/envs/napari/bin/python
PY_INGEST=/home/kjenie/miniconda3/bin/python
REPO=/home/kjenie/merfisheyes

mkdir -p "$OUT"
if [ $# -gt 0 ]; then EMBRYOS=("$@"); else
  mapfile -t EMBRYOS < <(ls "$ROOT/Final_analyzed_objects")
fi

ok=0; fail=0
printf '%-16s %-10s %10s %8s %8s %8s %8s\n' EMBRYO STATUS MOLECULES GENES DOMAINS CELLS MB
for e in "${EMBRYOS[@]}"; do
  pq="$OUT/${e}_molecules.parquet"
  dir="$OUT/${e}_lm"
  log="$OUT/${e}.log"

  {
    "$PY" "$REPO/scripts/spiralia/export_parquet.py" "$e" "$OUT" &&
    "$PY_INGEST" "$REPO/scripts/process_labelled_molecules.py" "$pq" "$dir" \
        --category gene=feature_name \
        --category domain=rna_domain_anno \
        --category cell=cell_id --label cell=cell_name \
        --drop-unassigned cell_id \
        --drop-control-genes feature_name \
        --name "$e"
  } >"$log" 2>&1

  if [ -f "$dir/manifest.json" ]; then
    read -r n g d c < <("$PY_INGEST" - "$dir" <<'PYEOF'
import json,sys
from pathlib import Path
d=Path(sys.argv[1])
m=json.load(open(d/'manifest.json')); o=json.load(open(d/'obs/metadata.json'))
print(m['statistics']['total_cells'],
      o.get('gene',{}).get('unique_values','-'),
      o.get('domain',{}).get('unique_values','-'),
      o.get('cell',{}).get('unique_values','-'))
PYEOF
)
    mb=$(du -sm "$dir" | cut -f1)
    printf '%-16s %-10s %10s %8s %8s %8s %8s\n' "$e" OK "$n" "$g" "$d" "$c" "$mb"
    ok=$((ok+1))
    rm -f "$pq"
  else
    printf '%-16s %-10s %s\n' "$e" FAILED "$(tail -2 "$log" | tr '\n' ' ' | cut -c1-90)"
    fail=$((fail+1))
  fi
done
echo
echo "built $ok, failed $fail"
