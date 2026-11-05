#!/bin/bash
# Full-weight scans (zeros <= 1e4, margin 0.15) from u = 27 (the 1e12 count) up to each first region found by team_scan.
cd "$(dirname "$0")"
PY=${PYTHON:-python}; export OMP_NUM_THREADS=${NT:-12}
while pgrep -f run_scans_rho8.sh > /dev/null; do sleep 60; done
for job in "156 45.80" "132 55.94" "56 76.80" "40 97.84" "165 108.47" "105 180.55" "84 227.24"; do
  set -- $job
  echo "=== q=$1 to $2 $(date +%T)"
  $PY fullscan.py $1 27 $2 10000 0.15 > scans/full_q$1.json
  $PY -c "
import json; d = json.load(open('scans/full_q$1.json'))
print('  max_E %.4f at %.5f; near:' % (d['max_E'], d['u_at_max']), [(round(r['u'], 4), round(r['E_max'], 3)) for r in d['near']][:15])"
done
echo "FULLSCANS DONE $(date)"
