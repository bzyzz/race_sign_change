#!/usr/bin/env python3
"""Two-scale criterion for the residue-vs-nonresidue ("team") race modulo q.

Let R be the squares in (Z/qZ)^x and rho = [G : R], the number of real characters mod q. Put
G(x) := rho #{p <= x : p mod q in R} - #{p <= x : p not dividing q}. The residues lead (their average count exceeds
that of the non-residues) iff G(x) > 0. Since 1_R = (1/rho) sum_{chi real} chi,

  G(x) = sum_{chi in X_q} pi(x, chi),   pi(x, chi) = sum_{p <= x} chi(p),

with X_q the rho - 1 real nonprincipal characters mod q. For q in {8, 12} this is the G of ../skewes_bound/race/team
(phi = rho = 4). Write pi(x, chi) = S(chi; x) - (1/2) #{p <= sqrt x : p not dividing q} - sum_{m >= 3} ..., use
chi(p)^2 = 1, replace S(chi) by S(chi*) for the primitive chi* inducing chi (FH (2.13) corrected) and use the corrected
explicit formula (../skewes_bound/race/THEOREM.md section 2). Then

  u e^{-u/2} G(e^u) >= -(rho-1) P(u) - sum_{chi in X_q} sum_rho u e^{-u/2} li(e^{rho u}) - eps0(u),

where P(u) = (1/2) u e^{-u/2} pi(e^{u/2}) (the +(rho-1) omega(q)/2 term is dropped) and eps0 collects
u e^{-u/2} [e^{-u} + (1-delta) log u + |C_chi| + n_chi (log u + 1.37)] per character and the prime powers m >= 3,
|sum_chi sum_{m >= 3} chi(p)^m/m| <= (rho-1) 1.3 x^{1/3}/log x (FH (2.12), x >= e^30; we use phi(q) >= rho-1).
Integrating against the nonnegative kernel K (mass 1, support [omega-eta0, omega+eta0]) gives

  I(omega) >= L := -(rho-1) P_max - sum_chi S_chi(omega, T0) - R3 - R4 - R5 - eps0,

and L > 0 implies G(e^u) > 0 for some u in [omega - eta0, omega + eta0]. P_max bounds P on the window:
Dusart's Phi_D(u_lo) if u_lo >= 45 (e^{u/2} >= 4e9), otherwise (1/2) u_hi e^{-u_lo/2} pi(e^{u_hi/2}) with pi counted
exactly (prime_pi.py). R3, R4, R5 are the race bounds (race_criterion.py) with coefficient 1 per character.
"""
import math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from flint import arb, dirichlet_char
import race_criterion as rc
from race_criterion import Phi_D, Phi2m, C_M, bmor_E, bmor_dE, C_const, require

COEF = 1     # each character enters with coefficient 1


def platt_height(d):
    """GRH height for conductor d <= 400000 (Platt 2016, Thm 7.1): max(1e8/d, c/d + 200), c = 7.5e7 (even d) or
    3.75e7 (odd d); rounded down to two significant digits."""
    h = max(1e8/d, (7.5e7 if d % 2 == 0 else 3.75e7)/d + 200)
    e = 10**(int(math.floor(math.log10(h))) - 1)
    return arb(int(h // e)*e)


def characters(q):
    """[(d, n, parity, A_chi, n_imprimitive)] for the real primitive characters d.n (d > 1, d | q) inducing the real
    nonprincipal characters mod q; n_imprimitive = #{p | q : p does not divide d}."""
    ps = [p for p in range(2, q + 1) if q % p == 0 and all(p % r for r in range(2, int(p**0.5) + 1))]
    out = []
    for d in range(3, q + 1):
        if q % d:
            continue
        for n in range(1, d + 1):
            if math.gcd(n, d) == 1:
                c = dirichlet_char(d, n)
                if c.is_primitive() and c.is_real():
                    out.append((d, n, int(c.parity()), platt_height(d), sum(1 for p in ps if d % p)))
    return out


def phi(q):
    return sum(1 for a in range(1, q + 1) if math.gcd(a, q) == 1)


def R3_bound(u_lo, eta0, sums2, sums3):
    tot = arb(0)
    for s2, s3 in zip(sums2, sums3):
        tot += 2*s3 + eta0/2*s2
    return COEF*tot/(u_lo*u_lo)


def R4_bound(chars, omega, eta1, T0, lam, n_grid=6000):
    pi = arb.pi()
    m_max = max(m for m, l in enumerate(lam, 1) if l != 0)
    require(eta1*T0 > max(3, m_max)*pi, f"eta1 T0 must exceed {max(3, m_max)} pi")
    H = lambda xi: sum(l*Phi2m(xi, m) for m, l in enumerate(lam, 1) if l != 0)
    ph = lambda t: H(eta1*t)*(1/t + 1/(t*t*omega))
    tot = arb(0)
    for (q, n, d, A, _) in chars:
        require(T0 < A, f"T0 must be below the GRH height A_chi for {q}.{n}")
        r = (A/T0)**(arb(1)/n_grid); t = T0; integ = arb(0)
        for k in range(n_grid):
            t2 = t*r
            integ += ph(t)*((q*t2/(2*pi)).log()/(2*pi) + bmor_dE(t, q)/2)*(t2 - t)
            t = t2
        tot += 2*COEF*(integ + bmor_E(T0, q)*ph(T0))
    return tot.upper()


def R5_bound(chars, u_lo, u_hi, delta, M):
    pi = arb.pi(); N = 2*M; R = u_lo/2; tot = arb(0)
    for (q, n, d, A, _) in chars:
        rate = (M*pi/delta + N/R)/A
        require(rate < 1, "R5 rate must be < 1")
        tail = (N + 1)/pi*((q*A/(2*pi)).log()/N + 1/arb(N)**2)
        tot += COEF*(u_hi/2).exp()*4*C_M(M)*rate**N*tail
    return tot.upper()


def eps0_bound(chars, phi_q, u_lo, u_hi):
    tot = arb(0)
    for (q, n, d, A, nim) in chars:
        tot += COEF*((-u_lo).exp() + (1 - d)*u_hi.log() + C_const(q, n)) + COEF*nim*(u_hi.log() + arb("1.37"))
    # prime powers m >= 3: (rho - 1) 1.3 x^{1/3}/log x <= phi 1.3 x^{1/3}/log x (FH (2.12)), times u e^{-u/2}
    return (u_lo*(-u_lo/2).exp()*tot + phi_q*arb("1.3")*(-u_lo/6).exp()).upper()


def P_max(u_lo, u_hi):
    """Upper bound for (1/2) u e^{-u/2} pi(e^{u/2}) on [u_lo, u_hi]."""
    if u_lo >= 45:
        return Phi_D(u_lo).upper()      # Dusart: pi(y) <= (y/log y)(1 + 1/log y + 2/log^2 y + 7.32/log^3 y), y >= 4e9
    import prime_pi
    y = int(math.floor(float((u_hi/2).exp().upper()))) + 1
    return (u_hi*(-u_lo/2).exp()*prime_pi.pi_exact(y)/2).upper()


def errors(q, omega, eta1, delta, M, T0, lam, sums2, sums3):
    chars = characters(q); rho = len(chars) + 1
    omega, eta1, delta, T0 = arb(omega), arb(eta1), arb(delta), arb(T0)
    eta0 = eta1 + delta; u_lo, u_hi = omega - eta0, omega + eta0
    # u_lo >= 30: Ford-Hudson (2.12) needs x >= e^30; (2.10) needs x > 10; P_max switches from Dusart (u >= 44.22)
    # to an exact prime count below 45; the R5 disk argument and the li expansions hold for any u_lo > 0.
    require(u_lo >= 30, "need u_lo >= 30")
    out = dict(bias=((rho - 1)*P_max(u_lo, u_hi)).upper(), R3=R3_bound(u_lo, eta0, sums2, sums3).upper(),
               R4=R4_bound(chars, omega, eta1, T0, lam), R5=R5_bound(chars, u_lo, u_hi, delta, M),
               eps0=eps0_bound(chars, phi(q), u_lo, u_hi))
    out["total"] = (out["bias"] + out["R3"] + out["R4"] + out["R5"] + out["eps0"]).upper()
    return out


def min_delta(chars, u_lo, u_hi, M, target):
    lo, hi = 1e-7, 1e-1
    for _ in range(60):
        mid = math.sqrt(lo*hi)
        try:
            ok = R5_bound(chars, arb(u_lo), arb(u_hi), arb(mid), M) <= target
        except ValueError:
            ok = False
        lo, hi = (lo, mid) if ok else (mid, hi)
    return hi
