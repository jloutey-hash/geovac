"""Parametric hybrid closed forms: geovac.two_center_eri's own shell algebra, run with SYMBOLIC exponents.

The evaluator (`hybrid_closed_form_shell`) is used unchanged; only its three radial helpers are swapped for
versions that accept sympy expressions for Z and drop the per-orbital normalisation (a nonzero constant factor,
irrelevant to the order in R).  Every sign decision the evaluator makes (`is_negative`, `> 0`, `== 0`) is then
taken by sympy on the symbolic rates.  We make those decisions decidable by parametrising each OPEN region of the
exponent ratio with one positive symbol u, so that every rate combination the algebra meets is +- a positive
combination of {1, u}:

  rates  a_c = Z_A/n_c,  a_p = Z_A/n_a + Z_A/n_b  (the same-centre pair),  a_d = Z_B/n_d.  By scaling (R -> s R)
  only the ratio a_d/a_c matters; the coincidences are a_d = a_c and a_d = a_c + a_p.

    region "lt"  (a_d < a_c)            : a_d = 1,             a_c = 1 + u
    region "mid" (a_c < a_d < a_c + a_p): a_p = 1 + u,         a_c = (a_c/a_p)(1 + u),  a_d = a_c + u
    region "gt"  (a_d > a_c + a_p)      : a_c = 1,             a_d = 1 + a_p/a_c + u
    point  "eq1" (a_d = a_c)            : a_c = a_d = 1
    point  "eq2" (a_d = a_c + a_p)      : a_c = 1,             a_d = 1 + a_p/a_c

u ranges over (0, oo) and sweeps the whole open region exactly once (up to the scaling).
Nothing here edits a paper, test or doc.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path
from math import comb, factorial

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import geovac.two_center_eri as T  # noqa: E402

u = sp.Symbol("u", positive=True)
R = T.R_s


def _sym(Z):
    if isinstance(Z, Fraction):
        return sp.Rational(Z.numerator, Z.denominator)
    return sp.sympify(Z)


def radial_poly_sym(Z, n: int, l: int):
    Z = _sym(Z)
    lam = 2 * Z / n
    lag = [sp.Rational((-1) ** j * comb(n + l, n - l - 1 - j), factorial(j)) for j in range(n - l)]
    return {l + j: sp.expand(c * lam ** j) for j, c in enumerate(lag)}, Z / n


def radial_norm_one(Z, n: int, l: int):
    return sp.Integer(1)


def radial_product_sym(Z1, n1, l1, Z2, n2, l2):
    c1, a1 = radial_poly_sym(Z1, n1, l1)
    c2, a2 = radial_poly_sym(Z2, n2, l2)
    out = {}
    for k1, v1 in c1.items():
        for k2, v2 in c2.items():
            out[k1 + k2] = sp.expand(out.get(k1 + k2, 0) + v1 * v2)
    return out, sp.expand(a1 + a2)


T.radial_poly = radial_poly_sym
T.radial_norm = radial_norm_one
T.radial_product = radial_product_sym


# sympy's expand multiplies out DENOMINATORS too: c/(2u+1)^2 * e^{-...} becomes c/(4u^2 e^{..} + 4u e^{..} + e^{..}),
# which the evaluator's term separator cannot read.  The proxy below re-factors any expanded denominator that is
# a sum (factor treats the exp's as generators, so it returns monomial * (polynomial in u)^k).  Same function,
# different spelling; the evaluator's algebra is untouched.
import types  # noqa: E402

_real_expand = sp.expand


_FCACHE: dict = {}


def _fix_term(t):
    n, d = t.as_numer_denom()
    if not d.has(sp.Add):
        return t
    if d not in _FCACHE:
        _FCACHE[d] = sp.factor(d)
    return n / _FCACHE[d]


def expand_fixed(e, *args, **kw):
    ee = _real_expand(e, *args, **kw)
    if not isinstance(ee, sp.Basic):
        return ee
    terms = sp.Add.make_args(ee)
    if not any(tt.as_numer_denom()[1].has(sp.Add) for tt in terms):
        return ee
    return sp.Add(*[_fix_term(tt) for tt in terms])


class _SpProxy(types.ModuleType):
    def __getattr__(self, name):
        return getattr(sp, name)


_proxy = _SpProxy("sympy_proxy")
_proxy.expand = expand_fixed
T.sp = _proxy

ORB = {"1s": (1, 0, 0), "2s": (2, 0, 0), "3s": (3, 0, 0), "2p0": (2, 1, 0), "2p1": (2, 1, 1),
       "2p-1": (2, 1, -1), "3p0": (3, 1, 0), "3p1": (3, 1, 1), "3d0": (3, 2, 0), "3d1": (3, 2, 1),
       "3d-1": (3, 2, -1), "3d2": (3, 2, 2), "3d-2": (3, 2, -2)}


def rates_of(oa, ob, oc, od):
    """(a_p, a_c, a_d) per unit Z_A, Z_B: a_p = Z_A*kp, a_c = Z_A*kc, a_d = Z_B*kd."""
    kp = sp.Rational(1, oa[0]) + sp.Rational(1, ob[0])
    return kp, sp.Rational(1, oc[0]), sp.Rational(1, od[0])


def region_Z(region, oa, ob, oc, od):
    """(Z_A, Z_B, a_p, a_c, a_d) as expressions in u for the region."""
    kp, kc, kd = rates_of(oa, ob, oc, od)
    rho = kp / kc                                    # a_p / a_c
    if region == "lt":
        ac, ad = 1 + u, sp.Integer(1)
    elif region == "mid":
        ac = (1 + u) / rho
        ad = ac + u
    elif region == "gt":
        ac, ad = sp.Integer(1), 1 + rho + u
    elif region == "eq1":
        ac, ad = sp.Integer(1), sp.Integer(1)
    elif region == "eq2":
        ac, ad = sp.Integer(1), 1 + rho
    else:
        raise ValueError(region)
    ZA, ZB = ac / kc, ad / kd
    return sp.expand(ZA), sp.expand(ZB), sp.expand(rho * ac), sp.expand(ac), sp.expand(ad)


def shell_parts(ZA, oa, ob, oc, ZB, od, log=print):
    """`T.hybrid_closed_form_shell` step for step (same calls, same order), but returning the list of
    unexpanded pieces whose sum is the closed form (the final global sp.expand is the slow step), with timings."""
    import time
    t0 = time.time()
    w_s = sp.Symbol("w", positive=True)
    x_s = T.x_s
    parts = []
    for coeff, L, Lp, Mp, ld, md in T.hybrid_terms(ZA, oa, ob, oc, ZB, od):
        rad1, b1 = None, None
        for LL, _M, _g, rad, b in T.multipole_decomposition(ZA, *oa, ZA, *ob):
            if LL == L:
                rad1, b1 = rad, b
                break
        weight = (4 * sp.pi / (2 * L + 1) * sum(c * x_s ** (k + 2) for k, c in rad1.items()) * sp.exp(-b1 * x_s))
        e_lt = T.shell_inner_double(ZA, oa, ob, oc, ZB, od, L, Lp, Mp, ld, md, x_s, R, branch="lt")
        e_gt = T.shell_inner_double(ZA, oa, ob, oc, ZB, od, L, Lp, Mp, ld, md, x_s, R, branch="gt")
        assert not e_gt.atoms(sp.Ei)
        log(f"  term L={L} L'={Lp}: inner doubles {time.time() - t0:.1f}s")
        lt = expand_fixed(weight * e_lt)
        gt = expand_fixed(weight * e_gt)
        for e in (lt, gt):
            for term in sp.Add.make_args(e):
                pw = sum(int(f.exp) if (f.is_Pow and f.base == x_s) else (1 if f == x_s else 0)
                         for f in sp.Mul.make_args(term))
                assert pw >= 0
        p1 = T._integrate_with_e1(lt, x_s, R)
        log(f"  term L={L}: x<R integral {time.time() - t0:.1f}s")
        p2 = T._integrate_with_e1(expand_fixed(gt.subs(x_s, R + w_s)), w_s, sp.oo)
        log(f"  term L={L}: x>R integral {time.time() - t0:.1f}s")
        parts.append(2 * sp.pi / R * coeff * p1)
        parts.append(2 * sp.pi / R * coeff * p2)
    return parts


def hybrid_param(case, region, log=print):
    a, b, c, d = (ORB[x] for x in case)
    ZA, ZB, ap, ac, ad = region_Z(region, a, b, c, d)
    parts = shell_parts(ZA, a, b, c, ZB, d, log=log)
    return parts, dict(ZA=ZA, ZB=ZB, a_p=ap, a_c=ac, a_d=ad)


if __name__ == "__main__":
    import time
    t0 = time.time()
    parts, info = hybrid_param(("2p0", "2p0", "1s", "1s"), sys.argv[1] if len(sys.argv) > 1 else "mid",
                               log=lambda m: print(m, flush=True))
    print(info, flush=True)
    print("seconds", round(time.time() - t0, 1), flush=True)
