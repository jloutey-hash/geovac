"""IS1 remainder: the Mellin-route operator of F_n, n = 6, 7, 8, checked on F_n by quadrature.
F_n(t) = int_0^oo K0(x)^(n-1) K0(t x) dx; theta = t d/dt acts only on K0(t x).
Reports the relative residual |sum| / sum|terms| at several t, and a perturbed-coefficient control."""
import sys, json, time
sys.path.insert(0, "tests")
import mpmath as mp
import sympy as sp
from test_paper59_four_density_rung import _t_ode, u

DPS = 45
mp.mp.dps = DPS
TS = [mp.mpf("0.5"), mp.mpf("0.7"), mp.mpf("1.3"), mp.mpf("1.7")]

def theta_tables(J):
    al, be = [sp.Integer(1)], [sp.Integer(0)]   # theta^j K0(u) = al K0(u) + be u K1(u)
    for _ in range(J):
        a_, b_ = al[-1], be[-1]
        al.append(sp.expand(u * sp.diff(a_, u) - u**2 * b_))
        be.append(sp.expand(u * sp.diff(b_, u) - a_))
    return ([sp.lambdify(u, a_, "mpmath") for a_ in al], [sp.lambdify(u, b_, "mpmath") for b_ in be])

def residual(Q, n, T, alf, bef, J, pert=None):
    pts = [0, mp.mpf(1)/4, mp.mpf(1)/2, 1, 2, 4, 8, 16, mp.inf]
    cache = CACHE.setdefault(T, {})       # Bessel values memoised per node: all J+1 integrals share nodes
    def bes(x):
        if x not in cache:
            cache[x] = (mp.besselk(0, x), mp.besselk(0, T * x), mp.besselk(1, T * x))
        return cache[x]
    def integrand(x, j):
        k0, k0t, k1t = bes(x)
        return k0**(n - 1) * (alf[j](T * x) * k0t + bef[j](T * x) * T * x * k1t)
    thF = [mp.quad(lambda x, j=j: integrand(x, j), pts) for j in range(J + 1)]
    tot = mp.mpf(0); scale = mp.mpf(0)
    for i, q in enumerate(Q):
        for (deg,), cf in sp.Poly(q, u).terms():
            c = int(cf)
            if pert is not None and pert == (i, deg):
                c += 1
            term = T**(2 * i) * mp.mpf(c) * thF[deg]
            tot += term; scale += abs(term)
    return tot / scale, tot

CACHE = {}
out = {}
NS = [int(a) for a in sys.argv[1:]] or [6, 7, 8]
for n in NS:
    t0 = time.time()
    Q = _t_ode(n)
    J = max(sp.Poly(q, u).degree() for q in Q)
    alf, bef = theta_tables(J)
    # perturb the constant-in-theta term of the top t-power polynomial (any nonzero entry)
    top = len(Q) - 1
    terms = sp.Poly(Q[top], u).terms()
    pdeg = terms[-1][0][0]
    rec = {"order": J, "t_degree": 2 * top, "points": []}
    for T in TS:
        r, raw = residual(Q, n, T, alf, bef, J)
        rp, _ = residual(Q, n, T, alf, bef, J, pert=(top, pdeg))
        rec["points"].append({"t": str(T), "rel_residual": mp.nstr(r, 5), "raw": mp.nstr(raw, 5),
                              "perturbed_rel_residual": mp.nstr(rp, 5)})
        print(n, T, "rel", mp.nstr(r, 5), "perturbed", mp.nstr(rp, 5), flush=True)
    rec["seconds"] = round(time.time() - t0, 1)
    out[n] = rec
json.dump({"dps": DPS, "results": out}, open(f"debug/data/ising_class_is1_quad_{'_'.join(map(str, NS))}.json", "w"), indent=1)
