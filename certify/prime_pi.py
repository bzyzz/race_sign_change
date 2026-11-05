#!/usr/bin/env python3
"""Exact prime counting pi(y) for y up to a few 1e9 by a segmented sieve of Eratosthenes (numpy, integers only).
Checked against pi(10^k) = 4, 25, 168, 1229, 9592, 78498, 664579, 5761455, 50847534 (k = 1..9) by `prime_pi.py test`.
"""
import math, sys
import numpy as np


def pi_exact(y, seg=1 << 24):
    y = int(y)
    if y < 2:
        return 0
    r = math.isqrt(y)
    base = np.ones(r + 1, dtype=bool); base[:2] = False
    for i in range(2, math.isqrt(r) + 1):
        if base[i]:
            base[i*i::i] = False
    bp = np.nonzero(base)[0]
    count = 0
    for lo in range(0, y + 1, seg):
        hi = min(lo + seg, y + 1)
        s = np.ones(hi - lo, dtype=bool)
        if lo == 0:
            s[:min(2, hi)] = False
        for p in bp:
            p = int(p)
            if p*p >= hi:
                break
            start = max(p*p, (lo + p - 1)//p*p)
            s[start - lo::p] = False
        count += int(s.sum())
    return count


if __name__ == "__main__":
    if sys.argv[1:] == ["test"]:
        known = [4, 25, 168, 1229, 9592, 78498, 664579, 5761455, 50847534]
        for k, v in enumerate(known, 1):
            got = pi_exact(10**k)
            assert got == v, (k, got, v)
            print(f"pi(1e{k}) = {got} ok", flush=True)
    else:
        print(pi_exact(int(float(sys.argv[1]))))
