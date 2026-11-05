#!/bin/bash
# Item 4: truncated explicit formula vs exact lead regions (zones_1e12.txt), all rho >= 4 races, T = 5e3 then 1e4.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=${NT:-3}
for T in 5000 10000; do
  ${PYTHON:-python} near_miss.py zones_1e12.txt 27.6310211 $T $(cat reps_rho4plus.txt)
done
echo "NEAR_MISS DONE $(date)"
