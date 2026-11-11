#!/usr/bin/env python3
"""Floating-point feasibility test: can the races mod 24, 120, 168 be certified, with the cos^{2M} mollifier at large M
or with a Logan mollifier (the "hybrid" kernel K = core * Logan_delta)?

Nothing here is a certificate. The signal uses float zeros (census/zeros, T <= 1e4) and the error terms of
team_criterion.py (R3, R4, eps0 and the bias in Arb at the scan centre), with R5 replaced as follows.

Mollifiers (both nonnegative with mass 1 on [-delta, delta], so |psi^(t)| <= 1 for real t and R4 is unchanged):
  cos:   hat kappa_{2M}(xi) = Gamma(M+1)^2 / (Gamma(M+1+xi/pi) Gamma(M+1-xi/pi))   (exact; O(1) per zero for any M).
         R5 is the existing rigorous bound team_criterion.R5_bound (2M-fold integration by parts).
  logan: K_c(y) = c/(2 sinh c) I_0(c sqrt(1-y^2)) on |y| < 1 (Buethe 2015), scaled to [-delta, delta];
         hat = l_c(delta t),  l_c(t) = (c/sinh c) sin(sqrt(t^2-c^2))/sqrt(t^2-c^2).
         R5 here is a MODEL: the leading term of the zeros above A_chi, possibly off the line,
           2 e^{u_hi/2} e^{eta0/2} (c/sinh c) (1.01/(delta sqrt(1 - (c/(delta A))^2))) sum_{|gamma|>A} gamma^{-2},
         using |l_c(w)| <= 1.01 (c/sinh c)/|sqrt(w^2-c^2)| for |Im w| <= delta/2 and Re w >= 1.05 c, and
         sum_{|gamma|>A} gamma^{-2} <= (2/pi)(log(qA/2pi) + 1)/A. A rigorous version must also expand the slowly varying
         factor u e^{-u/2} li(e^{rho u}) e^{-(rho-1/2)u} (as in Buethe 2015, (7.2)-(7.8)); that changes lower-order terms.

usage: hybrid_search.py q uc [T0=10000] [half=0.006]
"""
import ctypes, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from flint import arb
import race_criterion as rc, team_criterion as tc

lib = ctypes.CDLL(os.path.join(HERE, "libscan.so")); dp = ctypes.POINTER(ctypes.c_double)
lib.scan.argtypes = [dp]*5 + [ctypes.c_long, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_long, ctypes.c_long, dp, dp, dp]
lgam = np.vectorize(math.lgamma)


def moll_cos(xi, M):
    x = np.abs(xi)/np.pi
    if np.any(x >= M + 1):
        raise ValueError("xi beyond the analytic range of the Gamma formula")
    return np.exp(2*math.lgamma(M + 1) - lgam(M + 1 + x) - lgam(M + 1 - x))


def moll_logan(t, c):
    t = np.abs(t)
    if np.any(t >= c):
        raise ValueError("delta * T0 >= c: not needed here")
    s = np.sqrt(c*c - t*t)
    return (c/s)*np.exp(s - c)*(-np.expm1(-2*s))/(-math.expm1(-2*c))


def tail_inv2(q, A):
    return 2/math.pi*(math.log(q*A/(2*math.pi)) + 1)/A


def R5_logan(chars, u_hi, eta0, delta, c):
    """log of the model R5 (see the docstring)."""
    tot = 0.0
    for (q, n, d, A, _) in chars:
        A = float(A.mid()); r = c/(delta*A)
        if r > 1/1.05:
            return math.inf
        lg = (u_hi/2 + eta0/2 + math.log(2*c) - c + math.log(1.01/(delta*math.sqrt(1 - r*r)))
              + math.log(tail_inv2(q, A)))
        tot = max(tot, lg) if tot else lg
    return tot + math.log(len(chars))       # crude: (number of characters) x the largest term


def load(d, n, T0):
    z = np.load(os.path.join(ROOT, "census", "zeros", f"z_{d}_{n}_T10000.npy"))
    return z[z <= T0]


def main():
    q, uc = int(sys.argv[1]), float(sys.argv[2])
    T0 = float(sys.argv[3]) if len(sys.argv) > 3 else 1e4
    half = float(sys.argv[4]) if len(sys.argv) > 4 else 6e-3
    chars = tc.characters(q); r = len(chars) + 1
    zs = [load(cq, cn, T0) for (cq, cn, *_ ) in chars]
    g_all = np.concatenate(zs); rho = 0.5 + 1j*g_all
    s2 = [arb(2*np.sum(g**-2.0)) + rc.sum_inv_pow_tail(arb(T0), cq, 2, 2*len(g)) for (cq, *_), g in zip(chars, zs)]
    s3 = [arb(2*np.sum(g**-3.0)) + rc.sum_inv_pow_tail(arb(T0), cq, 3, 2*len(g)) for (cq, *_), g in zip(chars, zs)]
    Amin = min(float(A.mid()) for (*_, A, _) in chars)
    print(f"q={q} rho={r} chars={len(chars)} zeros={len(g_all)} T0={T0:g} A_min={Amin:.3g} uc={uc}", flush=True)
    u_hi0 = uc + 0.05
    # cos^{2M}: R5 ~ exp(u/2 - 2M) at M = A delta/(e pi); floor delta = e pi u/(4 A)
    floor_cos = math.e*math.pi*uc/(4*Amin)
    # logan: need c > u/2 + log(...), delta A >= 1.05 c
    designs = []
    for f in (1.15, 1.3, 1.6):
        delta = f*floor_cos
        best = None
        for M in sorted({int(Amin*delta/(math.e*math.pi)*k) for k in (0.8, 0.9, 1.0, 1.1, 1.25)}):
            try:
                l5 = float((tc.R5_bound(chars, arb(uc - 0.05), arb(u_hi0), arb(delta), M)).log().mid())/math.log(10)
            except ValueError:
                continue
            if best is None or l5 < best[1]:
                best = (M, l5)
        if best:
            designs.append(("cos", delta, best[0], best[1]))
    for extra in (20, 40, 80):
        c = u_hi0/2 + extra
        for f in (1.06, 1.15, 1.3):
            delta = f*c/Amin
            l5 = R5_logan(chars, u_hi0, 0.01, delta, c)/math.log(10)
            designs.append(("logan", delta, c, l5))
    print(f"cos floor delta = {floor_cos:.3e};  logan floor delta ~ {u_hi0/2/Amin:.3e}")
    results = []
    for kind, delta, par, l5 in designs:
        if l5 > -6:
            print(f"  {kind:5s} delta={delta:.3e} par={par:.0f}: log10 R5 = {l5:.1f} (too large)"); continue
        R5 = 10**l5
        for cc in (10, 12, 14, 16, 20):
            eta1 = cc/T0
            for lam in [(0, 0, 1), (0.04, 0.52, 0.44), (0, 1, 0), (0, 0, 0, 1)]:
                m_max = len(lam)
                if eta1*T0 <= max(3, m_max)*math.pi:
                    continue
                core = sum(l*rc.khat_f(eta1*g_all, m) for m, l in enumerate(lam, 1) if l)
                try:
                    mol = moll_cos(delta*g_all, par) if kind == "cos" else moll_logan(delta*g_all, par)
                except ValueError as e:
                    print(f"  {kind} delta={delta:.2e}: {e}"); continue
                w = core*mol
                A_ = 2*w/rho; B_ = 2*w/rho**2
                h = eta1/20; npts = int(2*half/h) + 1; blk = 256; nb = (npts + blk - 1)//blk
                arrs = [np.ascontiguousarray(v, dtype=np.float64) for v in (g_all, A_.real, A_.imag, B_.real, B_.imag)]
                bmin, bu, full = np.empty(nb), np.empty(nb), np.empty(npts); p = lambda v: v.ctypes.data_as(dp)
                lib.scan(*[p(v) for v in arrs], len(g_all), 0.0, float(uc - half), float(h), npts, blk, p(bmin), p(bu), p(full))
                eta0 = eta1 + delta; om = arb(uc); u_lo, u_hi = om - eta0, om + eta0
                lamA = [arb(l) for l in lam]
                bias = float(((r - 1)*tc.P_max(u_lo, u_hi)).mid())
                R3 = float(tc.R3_bound(u_lo, arb(eta0), s2, s3).mid())
                R4 = float(tc.R4_bound(chars, om, arb(eta1), arb(T0), lamA, n_grid=1500).mid())
                e0 = float(tc.eps0_bound(chars, tc.phi(q), u_lo, u_hi).mid())
                Lcurve = -full - bias - R3 - R4 - R5 - e0
                k = int(Lcurve.argmax()); u = uc - half + h*np.arange(npts)
                sig = -full[k] - bias
                results.append((Lcurve[k], kind, delta, par, cc, lam, u[k], sig, R4))
                print(f"  {kind:5s} delta={delta:.3e} par={par:7.0f} eta1T0={cc} lam={lam}: max L={Lcurve[k]:+.4f} at u={u[k]:.6f} "
                      f"(signal {sig:+.4f}, R4 {R4:.4f}, log10 R5 {l5:.1f})", flush=True)
    results.sort(key=lambda t: -t[0])
    print("BEST:", results[:3])


if __name__ == "__main__":
    main()
