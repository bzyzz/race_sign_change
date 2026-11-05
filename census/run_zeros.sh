#!/bin/bash
# Float zeros (heuristic) up to height T for the characters listed in LIST (lines "d n"; default need_200.txt,
# written by races.py), or for every primitive character of conductor <= DMAX if LIST=none.
cd "$(dirname "$0")"
PY=${PYTHON:-python}; T=${T:-5000}; NP=${NP:-20}; DMAX=${DMAX:-30}; LIST=${LIST:-need_200.txt}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
if [ "$LIST" = none ]; then
  $PY list_chars.py $DMAX > chars_$DMAX.txt
  [ "${REAL_ONLY:-1}" = 1 ] && LINES=$(awk '$4 == 1 {print $1, $2}' chars_$DMAX.txt) || LINES=$(awk '{print $1, $2}' chars_$DMAX.txt)
else
  LINES=$(cat $LIST)
fi
echo "$LINES" | while read d n; do
  out=zeros/z_${d}_${n}_T${T}.npy
  [ -f "$out" ] && continue
  echo "=== $d.$n $(date +%T)"
  $PY float_zeros_any.py $d $n $T $NP $out | head -1
done
echo "ALL DONE $(date)"
