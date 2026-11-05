#!/usr/bin/env python3
"""Pointwise full-weight E(u) for the team race mod q near u0 with the certified zeros (data/arb_*), truncated at several T.
E(u) = -(rho-1)(1 + 2/u + 8/u^2) - sum_chi 2 Re sum_{0<gamma<=T} e^{i gamma u}(1/rho + 1/(rho^2 u)).  Heuristic only.
usage: recheck.py q u0 halfwidth step T1,T2,..."""
import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certify"))
import team_criterion as tc
q, u0, hw, st = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
Ts = [float(x) for x in sys.argv[5].split(",")]
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
chars = tc.characters(q); r1 = len(chars)
gs = []
for (d, n, *_) in chars:
    p = os.path.join(D, f"arb_{d}_{n}_T1e5")       # the 1e5 runs, where present
    z = np.load(os.path.join(p if os.path.isdir(p) else os.path.join(D, f"arb_{d}_{n}"), "certified.npz")); g = 0.5*(z["lo"] + z["hi"]); gs.append(g[g <= max(Ts)])
g = np.concatenate(gs); rho = 0.5 + 1j*g
u = u0 + np.arange(-int(hw/st), int(hw/st) + 1)*st
for T in Ts:
    m = g <= T; r = rho[m]; gg = g[m]
    E = np.empty(len(u))
    for k in range(0, len(u), 100):
        uu = u[k:k + 100, None]
        E[k:k + 100] = -r1*(1 + 2/u[k:k + 100] + 8/u[k:k + 100]**2) - 2*np.sum((np.exp(1j*gg*uu)*(1/r + 1/(r**2*uu))).real, axis=1)
    pos = u[E > 0]
    print(f"q={q} T={T:.0f}: max E {E.max():+.4f} at {u[E.argmax()]:.7f}; positive width {len(pos)*st:.1e}"
          + (f" [{pos[0]:.7f}, {pos[-1]:.7f}]" if len(pos) else ""), flush=True)
