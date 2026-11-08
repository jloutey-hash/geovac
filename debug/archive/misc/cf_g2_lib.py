"""CF5-CF7 helpers: the record-g32 libraries (correlated_triangle_g32_lib*, imported, not modified) with the
singular-point lists re-pointed to gamma = 2, plus local-data and right-factor survey utilities.

At gamma = 2 the candidate true singular points of Y's annihilator are t = |1-gamma|, 2, 1+gamma, 2 gamma
= 1, 2, 3, 4 (x = 1, 4, 9, 16), plus t = 0 and infinity.  lib3/lib4 hard-wire the gamma = 3/2 points;
set_points() overrides their module globals in this process only."""
from __future__ import annotations
import json
from fractions import Fraction
from pathlib import Path
import flint
import sympy as sp
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import correlated_triangle_g32_lib3 as Lc
import correlated_triangle_g32_lib4 as S

DATA = Path(__file__).parent / "data"
u = G.u_sym
xs = sp.Symbol("x")


def set_points(tpos):
    """tpos: positive t-values (Fractions) of the true singular points other than 0"""
    tpos = [Fraction(v) for v in tpos]
    Lc.XPTS = [v * v for v in tpos]
    S.TPTS = [s * v for v in tpos for s in (1, -1)]
    S.XOF = {a: a * a for a in S.TPTS}


def load_Q(path):
    d = json.loads(Path(path).read_text())
    return G.from_Q([sp.sympify(q, locals={"u": u}) for q in d["Q"]])


def lead_factors(L):
    return [(f, int(e)) for f, e in L[-1].factor()[1]]


def rational_roots(L):
    out = []
    for f, e in lead_factors(L):
        if f.degree() == 1:
            c = [int(v) for v in f.coeffs()]
            out.append((Fraction(-c[0], c[1]), e))
    return out


def apparent_check(L, f, nprimes=4, N=60):
    """for an irreducible factor f (fmpz_poly) of the leading coefficient of degree >= 2: at roots of f mod p,
    are all exponents non-negative integers < some bound with no log?  Returns list of (p, root, exps, logfree)."""
    res = []
    r = len(L) - 1
    for p in H.PRIMES[50:400]:
        fp = flint.nmod_poly([int(v) % p for v in f.coeffs()], p)
        if fp.degree() < f.degree():
            continue
        lin = [c for c, m in fp.factor()[1] if c.degree() == 1]
        if not lin:
            continue
        cc = [int(v) for v in lin[0].coeffs()]
        x0 = (-cc[0]) * G.inv(cc[1], p) % p
        Lp = G.op_mod(L, p)
        B, _ = G.hform(Lp, x0)
        ind = [int(c) % p for c in B[0]]
        # integer roots in [-60, 60]
        def ev(v):
            acc = 0
            for c in reversed(ind):
                acc = (acc * v + c) % p
            return acc
        ints = [n for n in range(-60, 61) if ev(n % p) == 0]
        logfree = []
        for e in ints:
            try:
                G.frob_mod(Lp, x0, Fraction(e), N, p)
                logfree.append(True)
            except ZeroDivisionError:
                logfree.append(False)
            except AssertionError:
                logfree.append(None)
        # multiplicity: count distinct integer roots; an apparent point of order r has r distinct ones
        res.append({"p": p, "int_exps": ints, "logfree": logfree})
        if len(res) >= nprimes:
            break
    return res


def local_report(L, xpts, verbose=True, name="L"):
    out = {"order": len(L) - 1, "x_degree": G.op_degree(L)}
    out["t0"] = [str(e) for e in Lc.exps_t0(L)[0]]
    out["inf_mu"] = [str(e) for e in Lc.exps_inf_growth(L)[0]]
    out["points"] = {}
    for x0 in xpts:
        out["points"][str(x0)] = {"exps": [str(e) for e in Lc.exps_finite_t(L, x0)[0]],
                                   "mono": Lc.monodromy_type(L, x0)}
    out["mono_t0"] = Lc.monodromy_type(L, Fraction(0))
    out["mono_inf"] = Lc.monodromy_type(L, "inf")
    out["lead"] = [(str(f), e) for f, e in lead_factors(L)]
    if verbose:
        print(f"== {name}: order {out['order']}, x-degree {out['x_degree']}", flush=True)
        print("   t=0 exponents:", out["t0"], flush=True)
        print("   infinity growth mu:", out["inf_mu"], flush=True)
        print("   t=0 classes:", [c["exponents"] for c in out["mono_t0"]], flush=True)
        print("   inf classes:", [c["exponents"] for c in out["mono_inf"]], flush=True)
        for k, v in out["points"].items():
            print(f"   x={k}: exps {v['exps']}; classes", [c["exponents"] for c in v["mono"]], flush=True)
        print("   leading factors:", [(f[:70], e) for f, e in out["lead"]], flush=True)
    return out


def min_annihilator_of_solution(L, x0, e, rmax, dmax, p=H.P0, rmin=1):
    fn = H.frob_series_fn(L, Fraction(x0), Fraction(e))
    return H.find_order_degree(fn, Fraction(x0), rmax, dmax, p=p, rmin=rmin, verbose=False)


def divides(L, R):
    m, Qt, rem = G.pseudo_rdiv(L, R)
    return G.is_zero_op(rem), (G.content_free(Qt) if G.is_zero_op(rem) else None)


def _null_at(series_fn, x0, r, d, p, margin=30):
    U = (r + 1) * (d + 1)
    N = U + margin + r + 4
    sl = series_fn(p, N)
    rows, labels = G.columns(sl, x0, p, r, d)
    return len(G.nullspace_mod(rows, U, p))


def fast_order_degree(series_fn, x0, rmax, dmax, p=H.P0, rmin=1):
    """smallest order r <= rmax with an annihilator of x-degree <= dmax (one test per r at d = dmax),
    then the minimal degree for that r by bisection.  Returns (r, d, nullity) or None."""
    x0 = Fraction(x0)
    for r in range(rmin, rmax + 1):
        if _null_at(series_fn, x0, r, dmax, p):
            lo, hi = -1, dmax
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if _null_at(series_fn, x0, r, mid, p):
                    hi = mid
                else:
                    lo = mid
            return r, hi, _null_at(series_fn, x0, r, hi, p)
    return None


def fast_min_annihilator(L, x0, e, rmax, dmax, p=H.P0, rmin=1):
    return fast_order_degree(H.frob_series_fn(L, Fraction(x0), Fraction(e)), Fraction(x0), rmax, dmax, p=p, rmin=rmin)
