// First lead of the residues in the team race mod q, by direct counting (exact; segmented sieve, OpenMP).
// G_q(x) = rho #{p <= x : p mod q is a square of a unit} - #{p <= x : p does not divide q}, rho = [units : squares].
// The residues lead (average count per class) iff G_q(x) > 0. For each modulus q given, report the first prime p
// with G_q(p) > 0, or 0 if there is none <= X, together with G_q(X) and max_{p <= X} G_q(p).
// Pass 1 computes, per block of BLK integers and per q, the sum of the increments and their maximal prefix sum.
// Pass 2 re-sieves the first block in which the running maximum becomes positive and locates the prime.
// usage: team_count X q1 [q2 ...]          cc -O3 -march=native -fopenmp team_count.c -o team_count -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

#define QMAXN 512
static int nq, qs[QMAXN], *wt[QMAXN];          // wt[i][a] = rho - 1 (a residue), -1 (other unit), 0 (non-unit)
static int64_t *bp; static int64_t nbp;
static const int64_t BLK = 1000000000LL;       // integers per block
static const int64_t SUB = 1LL << 21;          // odd numbers per sieve segment

static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

// Sieve the odd numbers in [lo, hi) (lo, hi even or odd; handled by index), call f on each prime.
// The callback is expressed through a mode: 0 = accumulate sum/max-prefix per q, 1 = locate first positive.
static void sieve_block(int64_t lo, int64_t hi, int mode, int64_t *sum, int64_t *mx, int64_t *start, int64_t *firstp) {
    char *seg = malloc(SUB);
    int64_t i0 = lo / 2, i1 = (hi + 1) / 2;           // odd n = 2i + 1 with lo <= n < hi  <=> i in [lo/2, (hi-1)/2]
    if ((hi - 1) % 2 == 0) i1 = (hi - 1) / 2; else i1 = (hi - 1) / 2 + 1;
    int64_t pref[QMAXN];
    for (int k = 0; k < nq; k++) { pref[k] = mode ? start[k] : 0; if (!mode) { sum[k] = 0; mx[k] = INT64_MIN; } }
    if (lo <= 2 && 2 < hi) {                          // the prime 2
        for (int k = 0; k < nq; k++) {
            pref[k] += wt[k][2 % qs[k]];
            if (!mode) { if (pref[k] > mx[k]) mx[k] = pref[k]; }
            else if (pref[k] > 0 && !firstp[k]) firstp[k] = 2;
        }
    }
    for (int64_t s = i0; s < i1; s += SUB) {
        int64_t e = s + SUB < i1 ? s + SUB : i1;
        memset(seg, 0, e - s);
        if (s == 0) seg[0] = 1;                      // n = 1
        for (int64_t k = 1; k < nbp; k++) {
            int64_t p = bp[k];
            int64_t j = (p * p - 1) / 2;             // index of p^2
            if (j >= e) break;
            if (j < s) { int64_t n0 = 2 * s + 1; int64_t m = (n0 + p - 1) / p * p; if (m % 2 == 0) m += p; j = (m - 1) / 2; }
            for (; j < e; j += p) seg[j - s] = 1;
        }
        for (int64_t j = s; j < e; j++) if (!seg[j - s]) {
            int64_t n = 2 * j + 1;
            for (int k = 0; k < nq; k++) {
                int w = wt[k][n % qs[k]];
                if (!w) continue;
                pref[k] += w;
                if (!mode) { if (pref[k] > mx[k]) mx[k] = pref[k]; }
                else if (pref[k] > 0 && !firstp[k]) firstp[k] = n;
            }
        }
    }
    if (!mode) for (int k = 0; k < nq; k++) sum[k] = pref[k];
    free(seg);
}

int main(int argc, char **argv) {
    int64_t X = (int64_t)atof(argv[1]);
    for (int a = 2; a < argc; a++) {
        int q = atoi(argv[a]), i = nq++; qs[i] = q; wt[i] = calloc(q, sizeof(int));
        int *isres = calloc(q, sizeof(int)), nu = 0, nr = 0;
        for (int y = 1; y < q; y++) if (gcd(y, q) == 1) { isres[(y * y) % q] = 1; nu++; }
        for (int b = 1; b < q; b++) if (gcd(b, q) == 1 && isres[b]) nr++;
        int rho = nu / nr;
        for (int b = 1; b < q; b++) if (gcd(b, q) == 1) wt[i][b] = isres[b] ? rho - 1 : -1;
        if (q == 1) wt[i][0] = 0;
        free(isres);
    }
    int64_t R = (int64_t)sqrt((double)X) + 2;
    char *bs = calloc(R + 1, 1); bp = malloc(sizeof(int64_t) * (R / 2 + 10));
    for (int64_t i = 2; i <= R; i++) if (!bs[i]) { bp[nbp++] = i; for (int64_t j = i * i; j <= R; j += i) bs[j] = 1; }
    int64_t nb = (X + BLK) / BLK;                    // blocks [b BLK, (b+1) BLK) covering [0, X]
    int64_t *sum = malloc(sizeof(int64_t) * nb * nq), *mx = malloc(sizeof(int64_t) * nb * nq);
    int64_t done = 0;
    #pragma omp parallel for schedule(dynamic, 1)
    for (int64_t b = 0; b < nb; b++) {
        int64_t lo = b * BLK, hi = (b + 1) * BLK; if (hi > X + 1) hi = X + 1;
        sieve_block(lo, hi, 0, sum + b * nq, mx + b * nq, NULL, NULL);
        #pragma omp atomic
        done++;
        if (done % 100 == 0) { fprintf(stderr, "%lld/%lld blocks\n", (long long)done, (long long)nb); }
    }
    printf("# X=%lld\n", (long long)X);
    for (int k = 0; k < nq; k++) {
        int64_t off = 0, gmax = INT64_MIN, first = 0;
        for (int64_t b = 0; b < nb; b++) {
            int64_t m = off + mx[b * nq + k];
            if (mx[b * nq + k] != INT64_MIN && m > gmax) gmax = m;
            if (!first && mx[b * nq + k] != INT64_MIN && m > 0) {   // locate in this block
                int64_t lo = b * BLK, hi = (b + 1) * BLK; if (hi > X + 1) hi = X + 1;
                int64_t st[QMAXN] = {0}, fp[QMAXN] = {0};
                int save = nq; int qsave = qs[0]; int *wsave = wt[0];
                qs[0] = qs[k]; wt[0] = wt[k]; nq = 1; st[0] = off;
                sieve_block(lo, hi, 1, NULL, NULL, st, fp);
                nq = save; qs[0] = qsave; wt[0] = wsave;
                first = fp[0];
            }
            off += sum[b * nq + k];
        }
        printf("team %d first %lld G(X) %lld maxG %lld\n", qs[k], (long long)first, (long long)off, (long long)gmax);
    }
    return 0;
}
