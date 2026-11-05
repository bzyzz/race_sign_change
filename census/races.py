#!/usr/bin/env python3
"""Distinct residue vs non-residue ("team") races modulo q <= QMAX, q != 2 (mod 4), grouped by their set of real
characters (q and q' give the same race iff they have the same real characters, e.g. 5 and 25).

Writes races_<QMAX>.json: one entry per race with the smallest modulus q, the other moduli, rho = number of real
characters (the index of the squares), and the labels d.n of the real primitive characters (conductor > 1).
Also writes need_<QMAX>.txt: the characters in the order their zeros should be computed (rho = 4 races first).
usage: races.py [QMAX]   (needs python-flint)
"""
import json, math, os, sys
from flint import dirichlet_char

HERE = os.path.dirname(os.path.abspath(__file__))
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200


def real_prims(q):
    out = []
    for d in range(3, q + 1):
        if q % d:
            continue
        for n in range(1, d + 1):
            if math.gcd(n, d) == 1:
                c = dirichlet_char(d, n)
                if c.is_primitive() and c.is_real():
                    out.append((d, n, int(c.parity())))
    return out


races = {}
for q in range(3, QMAX + 1):
    if q % 4 == 2:
        continue
    rp = real_prims(q)
    key = tuple((d, n) for d, n, _ in rp)
    if key in races:
        races[key]["moduli"].append(q)
    else:
        races[key] = dict(q=q, moduli=[q], rho=len(rp) + 1, chars=[f"{d}.{n}" for d, n, _ in rp],
                          parity={f"{d}.{n}": p for d, n, p in rp})
out = sorted(races.values(), key=lambda r: (r["rho"], r["q"]))
json.dump(out, open(os.path.join(HERE, f"races_{QMAX}.json"), "w"), indent=1)
order = []
for rho in (4, 8, 2, 16):
    for r in out:
        if r["rho"] == rho:
            order += [c for c in r["chars"] if c not in order]
with open(os.path.join(HERE, f"need_{QMAX}.txt"), "w") as f:
    for c in order:
        d, n = c.split(".")
        f.write(f"{d} {n}\n")
for rho in sorted({r["rho"] for r in out}):
    print(f"rho={rho}: {sum(r['rho'] == rho for r in out)} races")
print(len(order), "characters")
