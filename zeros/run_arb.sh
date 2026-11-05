#!/bin/bash
# Certified zeros (float location, Arb certification, Turing) up to T for the characters given as arguments "d.n".
# usage: run_arb.sh T NP d.n [d.n ...]      (run from anywhere; writes ../data/float_d_n_T.npy and ../data/arb_d_n)
# TAG (env, optional) is appended to the output names, e.g. TAG=_T1e5 -> data/arb_d_n_T1e5
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python}; T=$1; NP=$2; shift 2
TS=$(seq -s, 5000 5000 $T)
for c in "$@"; do
  d=${c%.*}; n=${c#*.}
  echo "=== $c start $(date)"
  $PY zeros/float_zeros.py $d $n $((T + 100)) $NP data/float_${d}_${n}${TAG}.npy
  $PY zeros/arb_certify.py $d $n data/float_${d}_${n}${TAG}.npy data/arb_${d}_${n}${TAG} $NP
  $PY zeros/turing.py data/arb_${d}_${n}${TAG}/certified.npz 40 $TS
  echo "=== $c done $(date)"
done
