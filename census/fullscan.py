#!/usr/bin/env python3
"""Full-weight (sharp cutoff) heuristic scan of E(u) for the team race mod q on [U0, U1], zeros <= T, step 0.02/T.
Reports every block whose maximum of E is above -MARGIN, in increasing u. Catches narrow regions a smoothed scan can miss.
usage: fullscan.py q U0 U1 [T MARGIN]"""
import json, sys, numpy as np
import team_scan as ts
q, U0, U1 = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
T = float(sys.argv[4]) if len(sys.argv) > 4 else 1e4
MARGIN = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1
h = 0.02/T; blk = 4096
bmin, bu, _ = ts.scan_negE(q, "full", T, U0, U1, h, blk=blk)
idx = np.nonzero(bmin < MARGIN)[0]
rows = [dict(u=float(bu[i]), E_max=float(-bmin[i])) for i in idx]
print(json.dumps(dict(q=q, U0=U0, U1=U1, T=T, n_blocks=len(bmin), max_E=float(-bmin.min()), u_at_max=float(bu[np.argmin(bmin)]),
                      near=rows[:60]), indent=1))
