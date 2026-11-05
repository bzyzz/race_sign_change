#!/usr/bin/env python3
"""How negative can the truncated explicit formula be at a true lead? (item 4)

For each team race mod q (rho >= 4), compare the full-weight truncated sum
  E_T(u) = -(rho-1)(1 + 2/u + 8/u^2) - sum_chi 2 Re sum_{0<gamma<=T} e^{i gamma u}(1/rho + 1/(rho^2 u))
on 20 <= u <= U1 (= log of the counting limit) with the exact lead regions from zones.c (merged at GAP = 1e-3).
Outputs, per q and T:
  - for every true region: max E_T over [log start, log end] (a lead the scan "sees" iff this is > -c);
  - alarm segments {E_T > -c} (merged at GAP): hits (contain a true region) and false alarms, for several c.
Floating point, heuristic; the exact data are the regions.
usage: near_miss.py zones_file U1 T q1 [q2 ...]     -> near_miss_T<T>.json (appends per q)
"""
import json, math, os, sys
import numpy as np
import team_scan as ts

HERE = os.path.dirname(os.path.abspath(__file__))
GAP = 1e-3
CS = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)


def regions(fn, q):
    out = []
    for l in open(fn):
        t = l.split()
        if t[0] == "zone" and int(t[1]) == q:
            out.append(dict(start=int(t[2]), end=int(t[3]), maxG=int(t[4]), n=int(t[5]), us=math.log(int(t[2])), ue=math.log(int(t[3]))))
    return out


def segments(mask, u):
    """Merged (gap < GAP) runs of True in mask, as [u_lo, u_hi]."""
    idx = np.nonzero(mask)[0]
    if len(idx) == 0:
        return []
    br = np.nonzero(np.diff(u[idx]) > GAP)[0]
    starts = np.concatenate([[idx[0]], idx[br + 1]]); ends = np.concatenate([idx[br], [idx[-1]]])
    return [(float(u[a]), float(u[b])) for a, b in zip(starts, ends)]


def main():
    fn, U1, T = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    out_fn = os.path.join(HERE, f"near_miss_T{int(T)}.json")
    res = json.load(open(out_fn)) if os.path.exists(out_fn) else {}
    for q in map(int, sys.argv[4:]):
        if str(q) in res:
            continue
        h = 0.02/T; U0 = 20.0
        _, _, negE = ts.scan_negE(q, "full", T, U0, U1, h, full=True)
        E = -negE; u = U0 + h*np.arange(len(E))
        regs = regions(fn, q)
        for r in regs:
            a = int(np.searchsorted(u, r["us"])); b = int(np.searchsorted(u, r["ue"], side="right"))
            a = max(0, min(a, len(E) - 1)); b = max(a + 1, min(b, len(E)))
            r["maxE"] = float(E[a:b].max())
        rows = {}
        for c in CS:
            segs = segments(E > -c, u)
            hit = [s for s in segs if any(r["us"] <= s[1] + GAP and r["ue"] >= s[0] - GAP for r in regs)]
            miss = [r for r in regs if r["maxE"] <= -c]
            rows[str(c)] = dict(alarms=len(segs), hits=len(hit), false_alarms=len(segs) - len(hit), missed_regions=len(miss),
                                alarm_measure=float((E > -c).sum()*h))
        res[str(q)] = dict(rho=len(ts.real_prims(q)) + 1, T=T, U1=U1, n_regions=len(regs), regions=regs, by_c=rows,
                           min_maxE_true=min((r["maxE"] for r in regs), default=None))
        print(q, "regions", len(regs), "min max E_T over true regions", res[str(q)]["min_maxE_true"],
              {c: (v["hits"], v["false_alarms"], v["missed_regions"]) for c, v in rows.items()}, flush=True)
        json.dump(res, open(out_fn, "w"))


if __name__ == "__main__":
    main()
