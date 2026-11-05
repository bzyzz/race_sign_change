// All lead zones of the team races mod q up to X: maximal runs of consecutive primes p with G_q(p) > 0, where
// G_q = rho #{p in R} - #{p not dividing q} (as in team_count.c). Primes from primesieve. Prints one line per zone with
// start prime >= ZMIN: "zone q start end maxG nprimes" (end = last prime of the run; the lead holds on [start, next prime)).
// At most CAP zones per modulus are printed (a final "capped q" line says so).
// Zones separated by less than a factor (1 + GAP) (i.e. GAP in log x) are merged into one region.
// usage: zones X ZMIN CAP GAP q1 [q2 ...]
// cc -O2 zones.c -I$ENV/include -L$ENV/lib -Wl,-rpath,$ENV/lib -lprimesieve -o zones
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <primesieve.h>

#define NQ 128
static int nq, qs[NQ], w[NQ][256];

int main(int argc, char **argv) {
    uint64_t X = (uint64_t)atof(argv[1]), ZMIN = (uint64_t)atof(argv[2]);
    long CAP = atol(argv[3]); double GAP = atof(argv[4]);
    for (int a = 5; a < argc; a++) {
        int q = atoi(argv[a]), k = nq++; qs[k] = q;
        int unit[256] = {0}, sq[256] = {0}, nu = 0, nr = 0;
        for (int b = 1; b < q; b++) { int g = q, h = b; while (h) { int t = g % h; g = h; h = t; } if (g == 1) { unit[b] = 1; nu++; } }
        for (int b = 1; b < q; b++) if (unit[b]) sq[(b * b) % q] = 1;
        for (int b = 1; b < q; b++) if (unit[b] && sq[b]) nr++;
        for (int b = 0; b < q; b++) w[k][b] = unit[b] ? (sq[b] ? nu / nr - 1 : -1) : 0;
    }
    int64_t G[NQ] = {0}, zmax[NQ] = {0}; uint64_t zs[NQ] = {0}, zlast[NQ] = {0}, zn[NQ] = {0}; long cnt[NQ] = {0};
    primesieve_iterator it; primesieve_init(&it);
    uint64_t p;
    while ((p = primesieve_next_prime(&it)) <= X) {
        for (int k = 0; k < nq; k++) {
            int c = w[k][p % qs[k]];
            if (!c) continue;
            G[k] += c;
            if (G[k] > 0) {
                if (!zs[k]) { zs[k] = p; zmax[k] = G[k]; zn[k] = 0; }
                if (G[k] > zmax[k]) zmax[k] = G[k];
                zlast[k] = p; zn[k]++;
            } else if (zs[k] && (double)p > (double)zlast[k]*(1.0 + GAP)) {   // region ended more than GAP ago
                if (zs[k] >= ZMIN && cnt[k] < CAP) {
                    printf("zone %d %llu %llu %lld %llu\n", qs[k], (unsigned long long)zs[k], (unsigned long long)zlast[k],
                           (long long)zmax[k], (unsigned long long)zn[k]);
                    if (++cnt[k] == CAP) printf("capped %d\n", qs[k]);
                }
                zs[k] = 0;
            }
        }
    }
    for (int k = 0; k < nq; k++)
        if (zs[k] && zs[k] >= ZMIN && cnt[k] < CAP)
            printf("zone %d %llu %llu %lld %llu open\n", qs[k], (unsigned long long)zs[k], (unsigned long long)zlast[k],
                   (long long)zmax[k], (unsigned long long)zn[k]);
    primesieve_free_iterator(&it);
    return 0;
}
