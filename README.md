# When do the quadratic residues first lead?

Code, data and certificates for a census of residue races modulo q ≤ 200.

For q ≢ 2 (mod 4), let R be the squares in (ℤ/qℤ)^×. Let ρ = ρ(q) be the index of R, which is the number of real characters mod q. The residues *lead* at x when the primes up to x in the residue classes outnumber, per class, those in the non-residue classes:

G_q(x) = ρ·#{p ≤ x : p mod q ∈ R} − #{p ≤ x : p ∤ q} = Σ_{χ real, χ ≠ χ₀} π(x, χ) > 0.

Chebyshev's bias keeps the non-residues ahead most of the time. The first lead x_q is Fiorilli's generalized Skewes number (Algebra Number Theory 8 (2014), Conj. 1.19). Determining it is Problem 9 of the comparative prime number theory problem list (PIMS 2024). Moduli with the same real characters (e.g. 5 and 25) give the same race. There are 111 distinct races with q ≤ 200: 46 with ρ = 2, 47 with ρ = 4, 16 with ρ = 8 and 2 with ρ = 16.

## Results

**ρ = 4.** All 47 races are settled: 44 exact values and 3 proven intervals (`census/team_count_*.txt`, `census/check_primesieve_*.txt`). The largest exact values:

| q | x_q |
|---|---|
| 20, 100 | 64,039,880,506,961 |
| 21, 63, 147, 189 | 4,243,973,906,701 |
| 28, 196 | 3,843,952,272,217 |
| 33, 99 | 296,644,194,181 |
| 55 | 4,347,138,491 |

The three proven intervals:
- **q = 8:** 10¹⁵ < x₈ < 10^28.1270. The upper bound is certified in the companion repository [`skewes_bound`](https://github.com/bzyzz/skewes_bound).
- **q = 12:** 10¹⁵ < x₁₂ < 10^147.5007, also from `skewes_bound`.
- **q = 15:** 10¹⁵ < x₁₅ < 10^47.9376, certified here.

Many races lead at a tiny prime, because the first few primes happen to be residues: 26 of the 47 lead below 100.

**ρ ≥ 8.** Seven races lead below 10¹² (88, 104, 136, 140, 152, 184, 195). The other eleven have no lead below 10¹⁵. For eight of these, certificates give proven upper bounds:

| q | 156 | 132 | 56 | 40 | 165 | 105 | 84 | 60 |
|---|---|---|---|---|---|---|---|---|
| x_q < 10^… | 19.8929 | 24.2940 | 28.1310 | 42.4953 | 47.0413 | 78.4148 | 98.6909 | 696.0560 |

For the last three, q = 24, 120 and 168, the first positive region of the explicit formula (near 10^5915, 10^4262 and 10^2321) is narrower than any mollifier allowed by the known GRH heights. They are out of reach of this unconditional method.

**Near-misses.** There are 672 exact lead regions with e²⁰ ≤ x ≤ 10¹² in the races with ρ ≥ 4. Against these, the truncated explicit formula behaves as follows (`census/near_miss_T*.json`):
- a positive sum (zeros ≤ 10⁴) always came with a lead;
- 25% of the true leads have a negative sum;
- none of them is below −0.15.

So a heuristic "first sign change" is reliable only after every near-miss above about −0.15 has been checked. The exact values for q = 20 and 28 lie at such near-misses.

## Method

- **Exact counts** (`census/team_count.c`). A segmented sieve of Eratosthenes over blocks of 10⁹ integers, run in parallel. For each block it records the sum of the increments of G_q and their maximal prefix sum. A serial pass finds the first block that crosses, and that block is sieved again to locate the prime. Everything is integer arithmetic.
- **Independent checks.** `check_primesieve.c` uses the same block scheme with the primesieve library, and `count_races.c` is single-threaded with separate bookkeeping. They reproduce every exact value above. `check_primesieve` gives the 10¹⁵ lower bounds.
- **Certificates** (`certify/`). A two-scale compactly supported kernel (a cos^{2m} core convolved with a high-order cos^{2M} mollifier) in a Lehman-type explicit-formula criterion. The code is the residue-race version of the reviewed race code in [`skewes_bound/race`](https://github.com/bzyzz/skewes_bound/tree/main/race).
  - The bias is (ρ − 1).
  - Each real character enters with coefficient 1.
  - Below Dusart's range, an exact count of π(√x) replaces his bound.
  - The core can be cos⁸ (m = 4), with the condition η₁T₀ > max(3, m)π.
  - Every quantity is computed in Arb ball arithmetic. Zeros are certified by Arb sign changes and proved complete by Turing's method (Trudgian and Rumely bounds). GRH heights are Platt's (2016).

## Layout

```text
census/races.py              the 111 distinct races with q <= 200 (races_200.json) and their characters
census/team_count.c          exact first leads (OpenMP segmented sieve)
census/check_primesieve.c    independent count with primesieve; used for the 1e15 bounds
census/count_races.c         single-threaded check (two-way and team races, q <= 30)
census/zones.c               every lead region (zones_1e12.txt)
census/float_zeros_any.py    float zeros of L(1/2+it, chi), any primitive chi; run_zeros.sh
census/team_scan.py          smoothed + full-weight scan of the truncated explicit formula
census/fullscan.py           full-weight scan on an interval; recheck.py: pointwise recheck with the certified zeros
census/near_miss.py          calibration against the exact lead regions (near_miss_T*.json)
census/ld_team.py            random-phase (large-deviation + Rice) prediction; model_calibration.json
census/*.txt, *.json         outputs of the runs (counts, scans, calibration, Fiorilli data)
certify/team_criterion.py    the criterion: characters, Platt heights, error terms R3, R4, R5, eps0, bias
certify/team_certify.py      Arb certificate -> team_cert_q<q>_T<T0>_w<omega>.json
certify/team_search.py       floating-point parameter search
certify/validate_team.py     check of the explicit formula against actual primes
certify/prime_pi.py          exact pi(y) for the bias term below Dusart's range
certify/race_criterion.py, certify.py   shared race code (from skewes_bound/race/criterion)
zeros/                       float location, Arb certification (arb_certify.py) and Turing check (turing.py)
data/arb_<d>_<n>[_T...]/     certified zero intervals (certified.npz) and Turing certificates (turing_*.json)
LITERATURE.md                literature sweep on first sign changes of prime races
```

## Reproducing

Requirements:
- Python 3 with `python-flint` 0.9 (Arb), `numpy` and `scipy`;
- a C compiler with OpenMP;
- for `check_primesieve` and `zones`, primesieve ≥ 12 (`conda install -c conda-forge primesieve`).

```bash
make PRIMESIEVE=$CONDA_PREFIX

# exact counts
census/team_count 1e12 $(cat census/race_reps_200.txt)          # all 111 races, ~1 h on 4 threads
census/check_primesieve 1e15 8 12 15 24 40 56 60 84 105 120 132 156 165 168   # ~33 h on 12 threads

# certificates (seconds each, from the stored zeros); output JSON in certify/
cd certify
python team_certify.py 156 20000 45.804 6e-4 3.3e-4 40 0,1,0 16
python team_certify.py 132 20000 55.938 6e-4 3e-4 40 0,1,0 16
python team_certify.py 56 20000 64.7731 6e-4 1.41e-4 40 0,0,1 16
python team_certify.py 40 20000 97.8481 6e-4 1.2e-4 40 0.04,0.52,0.44 16
python team_certify.py 165 20000 108.3154 6e-4 5.2e-4 40 0,0,1 16
python team_certify.py 15 100000 110.3802155 1.46e-4 4.8e-5 40 0,0,0,1 16
python team_certify.py 105 20000 180.5552 1e-3 4.92e-4 60 0,0,0,1 16
python team_certify.py 84 20000 227.2428 7e-4 5.2e-4 60 0,0,0,1 16
python team_certify.py 60 20000 1602.724 1e-3 3e-3 300 0,0,0,1 16
```

The arguments are q, T₀, ω, η₁, δ, M, the core weights (λ₁, …, λ₄) and the number of processes. A certificate states that the residues lead for some log x ∈ [ω − η₁ − δ, ω + η₁ + δ] when L > 0. Each command regenerates the corresponding JSON in `certify/` byte for byte.

To recompute the zeros from scratch, run `zeros/run_arb.sh T NP d.n ...`. It does float location, then Arb certification, then the Turing check, and writes `data/arb_d_n/`. At 2·10⁴ it takes minutes to an hour per character on 16 processes; at 10⁵, several hours.

## Status

These results are computational and have not been peer reviewed. The heuristic statements (first positive regions, near-miss rechecks, the random-phase model) are labelled as such and prove nothing. A manuscript is in preparation.
