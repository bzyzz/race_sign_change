#!/bin/bash
# Heuristic team-race scans (team_scan.py) for the rho >= 8 races with q <= 200 that have no lead below 1e12.
# UMAX is about 3x the random-phase prediction (ld_preds_T300.txt); q = 24 was scanned to 30000 earlier.
cd "$(dirname "$0")"
PY=${PYTHON:-python}; export OMP_NUM_THREADS=${NT:-12}
for job in "156 1000" "165 1000" "132 1000" "105 1500" "56 1500" "84 2000" "40 2000" "60 12000" "168 20000" "120 100000"; do
  set -- $job
  echo "=== q=$1 UMAX=$2 $(date +%T)"
  $PY team_scan.py $1 $2 2000 10000 0.3 > scans/team_q$1.out 2>&1
  $PY - "$1" "$2" <<'P'
import json, sys
r = json.load(open(f"scans/team_q{sys.argv[1]}_U{sys.argv[2]}.json"))
print("  candidates", r["n_candidates"], "best_jp", r["best_jp"])
print("  first", r["first"])
P
done
echo "SCANS DONE $(date)"
