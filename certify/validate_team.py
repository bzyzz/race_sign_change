# Validate u e^{-u/2} G(e^u) ~ -(rho-1) Phi(u) - sum_chi 2 Re sum_{gamma>0} K(gamma) e^{i gamma u}(1/rho + 1/(rho^2 u))
# by smoothing, where G = rho #{p in R} - #{p not dividing q} (team_criterion.py). Zeros: certified ones in ../data if
# present, else the float zeros of ../census/zeros (T = 1e4 is plenty: the kernel transform at eta*gamma = 200 is ~1e-14).
# usage: validate_team.py q [q ...]
import glob, math, os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import team_criterion as tc


def zeros(d, n):
    p = os.path.join(ROOT, "data", f"arb_{d}_{n}", "certified.npz")
    if os.path.exists(p):
        z = np.load(p); return 0.5*(z["lo"] + z["hi"])
    return np.load(sorted(glob.glob(os.path.join(ROOT, "census", "zeros", f"z_{d}_{n}_T*.npy")))[-1])


for q in map(int, sys.argv[1:]):
    chars = tc.characters(q); rho = len(chars) + 1
    R = {(a*a) % q for a in range(1, q) if math.gcd(a, q) == 1}
    for omega in (16.5, 17.3, 18.1):
        eta = 0.02; lo, hi = omega - eta, omega + eta; X = int(math.exp(hi)) + 2
        sv = np.ones(X + 1, bool); sv[:2] = False
        for i in range(2, int(X**0.5) + 1):
            if sv[i]: sv[i*i::i] = False
        P = np.nonzero(sv)[0]; P = P[np.gcd(P, q) == 1]
        w = np.where(np.isin(P % q, list(R)), rho - 1.0, -1.0)
        lp = np.log(P.astype(float))
        G0 = w[lp <= lo].sum(); jl = lp[lp > lo]; jc = w[lp > lo]
        C3 = 4**3*math.factorial(3)**2/(2*math.factorial(6))
        grid = np.linspace(lo, hi, 400001); K = C3*np.cos(np.pi*(grid - omega)/(2*eta))**6/eta
        cs = np.concatenate([[0], np.cumsum(jc)]); G = G0 + cs[np.searchsorted(jl, grid, side="right")]
        lhs = np.trapezoid(K*grid*np.exp(-grid/2)*G, grid)
        rhs = -(rho - 1)*np.trapezoid(K*(1 + 2/grid + 8/grid**2 + 48/grid**3), grid)
        for (d, n, *_ ) in chars:
            g = zeros(d, n); g = g[g <= 2e4]
            r = 0.5 + 1j*g; xi = eta*g; kh = np.sinc(xi/np.pi)
            for k in (1, 2, 3): kh = kh*(k*np.pi)**2/((k*np.pi)**2 - xi*xi)
            rhs -= 2*np.sum((kh*np.exp(1j*g*omega)*(1/r + 1/(r**2*omega))).real)
        print(f"q={q} rho={rho} omega={omega}: smoothed u e^(-u/2) G from primes = {lhs:+.4f}; bias + zero formula = {rhs:+.4f}; "
              f"diff {lhs - rhs:+.4f}", flush=True)
