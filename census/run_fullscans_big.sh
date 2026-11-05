#!/bin/bash
# Full-weight scans below the first regions of the four largest races (zeros <= T, margin 0.15).
cd "$(dirname "$0")"
PY=${PYTHON:-python}; export OMP_NUM_THREADS=${NT:-12}
for job in "60 1602.72 10000" "168 5345.40 5000" "24 13619.12 5000" "120 9814.28 5000"; do
  set -- $job
  echo "=== q=$1 to $2 T=$3 $(date +%T)"
  $PY fullscan.py $1 27 $2 $3 0.15 > scans/full_q$1.json
  $PY -c "
import json; d = json.load(open('scans/full_q$1.json'))
print('  max_E %.4f at %.5f; near:' % (d['max_E'], d['u_at_max']), [(round(r['u'], 4), round(r['E_max'], 3)) for r in d['near']][:15])"
done
echo "FULLSCANS DONE $(date)"
