"""Cold review of debug/tree_attain_memo.md: independent checks (reviewer's own code except the rebuild of a).

P0  toy counterexample: d/dt (t-p) d/dt kills log(t-p); and 1/(t-p) lies in C(t).1.
P1  (*): phi(b) has order weight(b), leading coefficient t^-w sigma(b)(1,-t), on random b (exact).
P2  star 4: a rebuilt with the memo's analyse(); EXACT pointwise check in the function ring
    Q[x,1/x][chi,K0,K1] that a (chi^3/x) contains no chi and equals the memo's g.
P3  star 4: phi(a) T(H) = T(aH) in mpmath at 60 digits, phi applied under the integral sign.
P4  top constant C for branched trees by nested Gauss-Laguerre (independent of the memo's tail script).
"""
from __future__ import annotations

import json
import random
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import mpmath as mp
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "debug"))
T0 = time.time()
OUT: Dict = {}


def log(*a) -> None:
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


# ---------------------------------------------------------------- P0
def p0() -> None:
    t, p = sp.symbols("t p")
    f = sp.log(t - p)
    r = sp.simplify(sp.diff((t - p) * sp.diff(f, t), t))
    OUT["P0_toy_log_killed"] = bool(r == 0)
    log("P0: d/dt (t-p) d/dt log(t-p) =", r)


# ---------------------------------------------------------------- P1
# operators in t: dict {(m, j): c} = sum c t^m theta_t^j ; theta t^m = t^m (theta + m)
def op_mul(A: Dict, B: Dict) -> Dict:
    out: Dict = {}
    for (m1, j1), c1 in A.items():
        for (m2, j2), c2 in B.items():
            # theta^j1 t^m2 = t^m2 (theta+m2)^j1
            for i in range(j1 + 1):
                c = c1 * c2 * sp.binomial(j1, i) * Fraction(m2) ** (j1 - i)
                if c != 0:
                    key = (m1 + m2, i + j2)
                    out[key] = out.get(key, 0) + Fraction(int(c.numerator), int(c.denominator)) if isinstance(c, Fraction) else out.get(key, 0) + Fraction(c)
    return {k: v for k, v in out.items() if v != 0}


def op_add(A: Dict, B: Dict, s=1) -> Dict:
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v != 0}


PHI_Z = {(-2, 2): Fraction(1)}
PHI_TH = {(0, 0): Fraction(-1), (0, 1): Fraction(-1)}


def binom(n, k):
    from math import comb
    return comb(n, k)


def op_mul_fast(A: Dict, B: Dict) -> Dict:
    out: Dict = {}
    for (m1, j1), c1 in A.items():
        for (m2, j2), c2 in B.items():
            for i in range(j1 + 1):
                c = c1 * c2 * binom(j1, i) * Fraction(m2) ** (j1 - i)
                if c != 0:
                    key = (m1 + m2, i + j2)
                    out[key] = out.get(key, 0) + c
    return {k: v for k, v in out.items() if v != 0}


def phi(b: Dict) -> Dict:
    """b = {(k, j): c} = sum c z^k theta^j  (z to the left)."""
    out: Dict = {}
    for (k, j), c in b.items():
        term = {(0, 0): Fraction(c)}
        for _ in range(k):
            term = op_mul_fast(term, PHI_Z)
        for _ in range(j):
            term = op_mul_fast(term, PHI_TH)
        out = op_add(out, term)
    return out


def p1() -> None:
    random.seed(7)
    ok_all = True
    tt = sp.symbols("t")
    for trial in range(12):
        w = random.randint(3, 9)
        b = {}
        for k in range(0, w // 2 + 1):
            for j in range(0, w - 2 * k + 1):
                if random.random() < 0.6:
                    b[(k, j)] = random.randint(-5, 5)
        # make sure there is a top-weight term
        k0_ = random.randint(0, w // 2)
        b[(k0_, w - 2 * k0_)] = random.choice([-3, -2, -1, 1, 2, 3])
        b = {k: v for k, v in b.items() if v != 0}
        P = phi(b)
        order = max(j for (_, j) in P)
        lead = sum(sp.Rational(c.numerator, c.denominator) * tt ** m for (m, j), c in P.items() if j == order)
        sig = sum(c * (-tt) ** j for (k, j), c in b.items() if 2 * k + j == w)   # sigma(1, -t)
        ok = (order == w) and sp.expand(lead - tt ** (-w) * sig) == 0
        ok_all &= bool(ok)
    # the homomorphism property itself: phi(theta z) = phi(z theta) + 2 phi(z)
    lhs = op_mul_fast(PHI_TH, PHI_Z)
    rhs = op_add(op_mul_fast(PHI_Z, PHI_TH), {k: 2 * v for k, v in PHI_Z.items()})
    hom = (lhs == rhs)
    OUT["P1_star_rule_12_random_b"] = bool(ok_all)
    OUT["P1_phi_respects_[theta,z]=2z"] = bool(hom)
    log("P1: (*) on 12 random b:", ok_all, "; phi respects [theta,z]=2z:", hom)


# ---------------------------------------------------------------- P2 exact pointwise check for stars
def theta_ring(e: Dict) -> Dict:
    """theta on sum c x^ex chi^ec K0^e0 K1^e1; theta chi = x^2 K1, theta K0 = -x K1, theta K1 = -x K0 - K1."""
    out: Dict = {}

    def add(k, c):
        out[k] = out.get(k, 0) + c

    for (ex, ec, e0, e1), c in e.items():
        if ex:
            add((ex, ec, e0, e1), ex * c)
        if ec:
            add((ex + 2, ec - 1, e0, e1 + 1), ec * c)
        if e0:
            add((ex + 1, ec, e0 - 1, e1 + 1), -e0 * c)
        if e1:
            add((ex + 1, ec, e0 + 1, e1 - 1), -e1 * c)
            add((ex, ec, e0, e1), -e1 * c)
    return {k: v for k, v in out.items() if v != 0}


def star_exact(n: int, a: Dict, g: Dict) -> Dict:
    """a = {(k,j): c}; H = chi^(n-1) x^(3-n) (hat normalisation, P_leaf = chi/x)."""
    cur = {(3 - n, n - 1, 0, 0): Fraction(1)}
    jmax = max(j for (_, j) in a)
    tot: Dict = {}
    for j in range(jmax + 1):
        for (k, jj), c in a.items():
            if jj == j:
                for (ex, ec, e0, e1), v in cur.items():
                    key = (ex + 2 * k, ec, e0, e1)
                    tot[key] = tot.get(key, 0) + Fraction(int(c)) * v
        cur = theta_ring(cur)
    tot = {k: v for k, v in tot.items() if v != 0}
    chi_free = all(k[1] == 0 for k in tot)
    # memo's g: sum_i g_i(z) K0^(n-1-i) D^i, D = -x K1
    gd: Dict = {}
    zz = sp.symbols("z")
    for i, p in g.items():
        for (kk,), cf in sp.Poly(p, zz).terms():
            key = (2 * kk + i, 0, n - 1 - i, i)
            gd[key] = gd.get(key, 0) + Fraction(int(cf)) * (-1) ** i
    gd = {k: v for k, v in gd.items() if v != 0}
    return {"aH_has_no_chi": bool(chi_free), "aH_equals_memo_g": bool(tot == gd), "n_terms_aH": len(tot),
            "aH_nonzero": bool(len(tot) > 0)}


# ---------------------------------------------------------------- P3 high-precision transfer for star 4
def chi_mp(x):
    return -x * mp.besselk(0, x) + (mp.pi * x / 2) * (mp.besselk(0, x) * mp.struvel(-1, x) + mp.besselk(1, x) * mp.struvel(0, x))


def theta_K_tables(mmax: int):
    """theta^m K0(y) = A_m(y) K0 + B_m(y) K1, A, B integer polynomials in y (lists of coefficients)."""
    y = sp.symbols("y")
    A, B = [sp.Integer(1)], [sp.Integer(0)]
    for _ in range(mmax):
        al, be = A[-1], B[-1]
        A.append(sp.expand(y * sp.diff(al, y) - y * be))
        B.append(sp.expand(y * sp.diff(be, y) - y * al - be))
    return y, A, B


def p3(n: int, a: Dict, g: Dict, tq: Fraction, dps: int = 60) -> Dict:
    mp.mp.dps = dps
    th = sp.symbols("theta")
    # phi(a) on a function of y = t x :  (t^-2 theta^2)^k = t^(-2k) prod_{i<k} (theta-2i)^2, times (-1-theta)^j
    tS = sp.Rational(tq.numerator, tq.denominator)
    pol = sp.Poly(sp.expand(sum(sp.Integer(int(c)) * tS ** (-2 * k) * sp.prod([(th - 2 * i) ** 2 for i in range(k)]) * (-1 - th) ** j
                                for (k, j), c in a.items())), th)
    coeffs = pol.all_coeffs()[::-1]
    y, A, B = theta_K_tables(len(coeffs) - 1)
    al = sp.Poly(sp.expand(sum(c * A[m] for m, c in enumerate(coeffs))), y)
    be = sp.Poly(sp.expand(sum(c * B[m] for m, c in enumerate(coeffs))), y)
    alc = [mp.mpf(sp.Rational(c).p) / mp.mpf(sp.Rational(c).q) for c in al.all_coeffs()]
    bec = [mp.mpf(sp.Rational(c).p) / mp.mpf(sp.Rational(c).q) for c in be.all_coeffs()]
    t = mp.mpf(tq.numerator) / tq.denominator
    zz = sp.symbols("z")
    gp = {i: [mp.mpf(int(c)) for c in sp.Poly(p, zz).all_coeffs()] for i, p in g.items()}

    def H(x):
        return chi_mp(x) ** (n - 1) * x ** (3 - n)

    def lhs_int(x):
        yv = t * x
        return (mp.polyval(alc, yv) * mp.besselk(0, yv) + mp.polyval(bec, yv) * mp.besselk(1, yv)) * H(x)

    def rhs_int(x):
        k0v, D = mp.besselk(0, x), -x * mp.besselk(1, x)
        s = mp.mpf(0)
        for i, cs in gp.items():
            s += mp.polyval(cs, x * x) * k0v ** (n - 1 - i) * D ** i
        return mp.besselk(0, t * x) * s

    pts = [0, mp.mpf(1) / 8, mp.mpf(1) / 2, 1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 400]
    lhs = mp.quad(lhs_int, pts, maxdegree=9)
    rhs = mp.quad(rhs_int, pts, maxdegree=9)
    big = max(abs(lhs_int(mp.mpf(v))) for v in (1, 4, 16, 32, 64))
    res = {"t": str(tq), "dps": dps, "lhs_phi_a_T_H": mp.nstr(lhs, 25), "rhs_T_aH": mp.nstr(rhs, 25),
           "rel_diff": mp.nstr(abs(lhs - rhs) / abs(rhs), 5), "max_integrand_sampled": mp.nstr(big, 5)}
    log(f"P3 n={n} t={tq}: lhs {res['lhs_phi_a_T_H']}  rhs {res['rhs_T_aH']}  rel {res['rel_diff']}")
    return res


# ---------------------------------------------------------------- P4 top constant by nested Gauss-Laguerre
from scipy.special import k0e  # noqa: E402

GL_X, GL_W = np.polynomial.laguerre.laggauss(48)


def top_constant(n: int, edges: List[Tuple[int, int]], root: int, xs=(100.0, 200.0, 400.0, 800.0, 1600.0)) -> Dict:
    adj = {v: [] for v in range(n)}
    for a_, b_ in edges:
        adj[a_].append(b_)
        adj[b_].append(a_)
    children, size = {}, {}

    def rec(v, p):
        children[v] = [w for w in adj[v] if w != p]
        size[v] = 1 + sum(rec(w, v) for w in children[v])
        return size[v]
    rec(root, -1)

    def wsc(u, s):
        """e^(k_u s) W_u^full(s): K0 times the full tails of the children."""
        out = k0e(s)
        for c in children[u]:
            out = out * esc(c, s)
        return out

    def esc(u, x):
        """e^(k_u x) E_u^full(x),  E = -(1/x) int_0^oo W(x+v) (x+v) v dv."""
        k = size[u]
        x = np.asarray(x, dtype=float)
        v = GL_X / k                                   # e^(-k v) dv -> nodes / k, weights / k
        s = x[..., None] + v
        return -(1.0 / x) * np.sum(GL_W / k * wsc(u, s) * s * v, axis=-1)

    vals = []
    for x in xs:
        h = x ** 2
        for w in children[root]:
            h = h * float(esc(w, np.array(x)))
        vals.append(h * x ** (-(5 - n) / 2))            # -> C (1 + c1/x + ...)
    # Richardson in 1/x with doubling
    tab = list(vals)
    lvl = 1
    while len(tab) > 1:
        tab = [(2 ** lvl * tab[i + 1] - tab[i]) / (2 ** lvl - 1) for i in range(len(tab) - 1)]
        lvl += 1
    prod = 1
    for u in range(n):
        if u != root:
            prod *= size[u] ** 2
    pred = (-1) ** (n - 1) * (np.pi / 2) ** ((n - 1) / 2) / prod
    return {"prod_k2": prod, "predicted_C": pred, "raw_at_x": dict(zip(map(str, xs), vals)), "richardson": tab[0],
            "rel_err": abs(tab[0] / pred - 1)}


def p4() -> None:
    cases = {
        "star5": (5, [(0, 1), (0, 2), (0, 3), (0, 4)], 0),
        "fork5_root0": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 0),
        "fork5_root_branch": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 1),
        "fork5_root3": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 3),
        "fork5_root4": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 4),
        "path5_centre": (5, [(0, 1), (1, 2), (0, 3), (3, 4)], 0),
        "spider6": (6, [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5)], 0),
        "tree7_branched": (7, [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (5, 6)], 0),
        "tree8_branched": (8, [(0, 1), (0, 2), (2, 3), (2, 4), (4, 5), (4, 6), (6, 7)], 0),
        # wrong-answer control: the rule with k_u = 1 everywhere (the star constant) must be far off
    }
    OUT["P4"] = {}
    for name, (n, edges, root) in cases.items():
        r = top_constant(n, edges, root)
        OUT["P4"][name] = r
        log(f"P4 {name}: prod k^2 = {r['prod_k2']}, predicted C = {r['predicted_C']:.10g}, "
            f"Richardson = {r['richardson']:.10g}, rel err {r['rel_err']:.2e}")


if __name__ == "__main__":
    p0()
    p1()
    p4()
    import tree_attain_algebra as TA  # noqa: E402  (memo's module model, used only to rebuild a and g)
    OUT["P2"], OUT["P3"] = {}, {}
    for name in ("star3_root_centre", "star4"):
        n, edges = TA.TREES[name]
        res, a, g = TA.analyse(TA.Tree(n, edges), name)
        ex = star_exact(n, a, g)
        OUT["P2"][name] = ex
        log(f"P2 {name}: {ex}")
        # P3 (the transfer in high precision) lives in tree_attain_review_p3.py
        (ROOT / "debug" / "data" / "tree_attain_review_main.json").write_text(json.dumps(OUT, indent=1, default=str))
    (ROOT / "debug" / "data" / "tree_attain_review_main.json").write_text(json.dumps(OUT, indent=1, default=str))
    log("done")
