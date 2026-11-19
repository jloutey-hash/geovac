"""Helpers for tests/test_paper58_exchange_reverse.py (the four mixed exchange integrals at a > b).

Convention: a = exponent of the l > 0 orbital on A, b = exponent of 1s on B, t = b/a.  The frozen closed forms
(tests/data/paper58_pachucki_ext_2026_11_18.json for 2psig, 3dsig; tests/data/paper58_exchange_m_2026_11_18.json
for 2ppi, 3ddel) are K(R; 1, t) = sum_i c_i(t, R) g_i(t, R) on the nine functions of Paper 58; here they are read
at t in (0, 1).

Two numerical evaluators, written separately:
  * `value_real`: the stored coordinates specialised at t = b/a, rescaled to (a, b), turned into a Func over Q with
    the helper's own nine functions (X.nine_functions), evaluated on the REAL branch: LC = ln((a+b)/|b-a|), Ei of
    a positive argument = principal value (mpmath ei).  Every base value is asserted to be a real mpf before it is
    used, so an expression that only becomes real after taking a real part is refused, not silently repaired.
    (X.evaluate takes ln((a+b)/(b-a)), complex at a > b, and is not used here.)
  * `K1_typed`: the nine functions typed from their definitions, coefficients evaluated from integer numerator /
    denominator polynomials of the stored strings, so t may be complex.  Branch 'real' (t real) or 'cont': the
    holomorphic continuation FROM t > 1: LC = Log(1+t) - Log(t-1); a key whose Ei argument w is positive at t > 1
    uses Ei(w) = gamma + Log w - Ein(-w) = Log w - Log(-w) - E1(-w); one whose argument is negative at t > 1 uses
    -E1(-w).  f0 = int_R^oo F^(-1)(s)/s ds in either branch.
Exact pieces: `checked_solves` wraps every solve_L0 / solve_D2 of the reviewed helpers so that each output is
re-checked as an exact identity and its free unknowns and keys are recorded.
Nothing here reads or imports debug/.
"""
from __future__ import annotations

import contextlib
import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Tuple

import mpmath as mp
import sympy as sp

import _exchange_m_nonzero as M
import _pachucki_extension as X

M0_CASES = ("2psig", "3dsig")
MN_CASES = ("2ppi", "3ddel")
CASES = M0_CASES + MN_CASES

_DATA = Path(__file__).resolve().parent / "data"
EXT = json.loads((_DATA / "paper58_pachucki_ext_2026_11_18.json").read_text())
MNZ = json.loads((_DATA / "paper58_exchange_m_2026_11_18.json").read_text())


def stored(case: str) -> Dict[str, Dict[str, str]]:
    """The frozen closed form {name: {power of R: rational function of t}} (a = 1, b = t)."""
    return EXT["sym"][case] if case in M0_CASES else MNZ["sym"][case]


def audit(case: str) -> dict:
    """The frozen audit over Q(t) of every intermediate jet (denominator factors, function keys)."""
    return EXT["denoms"][case] if case in M0_CASES else MNZ["audit"][case]


def mpq(v) -> mp.mpf:
    v = Fraction(v)
    return mp.mpf(v.numerator) / v.denominator


# ------------------------------------------------------------------ the stored form at (a, b)
def stored_coords(case: str, a: Fraction, b: Fraction) -> List[Dict[int, Fraction]]:
    """Coordinates of K(R; a, b) on the nine functions at (a, b): K(R; a, b) = a K(aR; 1, b/a) and
    f0(R; a, b) = a^-2 f0(aR; 1, b/a), so the coefficient of R^p g_i is a^(p+1) c_ip(b/a) for i >= 1 and
    a^(p+3) c_0p(b/a) for f0."""
    a, b = Fraction(a), Fraction(b)
    co = X.coeff_strings_at(stored(case), b / a)
    return [{p: (a ** 3 if i == 0 else a) * c * a ** p for p, c in d.items()} for i, d in enumerate(co)]


def coords_to_func(a: Fraction, b: Fraction, co) -> X.Func:
    ctx = X.Ctx("Q", Fraction(a), Fraction(b))
    K: X.Func = {}
    for (g, _), d in zip(X.nine_functions(ctx), co):
        for p, c in d.items():
            K = X.fadd(K, X.fscale(g, c, p))
    return K


# ------------------------------------------------------------------ real-branch evaluation of a Func over Q
def lc_real(a: Fraction, b: Fraction) -> mp.mpf:
    """LC = ln((a+b)/|b-a|), the real branch for either ordering."""
    a, b = Fraction(a), Fraction(b)
    return mp.log(mpq((a + b) / abs(b - a)))


def _real(v, what):
    if not isinstance(v, mp.mpf):
        raise AssertionError(f"{what} is not real: {v}")
    return v


def evaluate_real(a: Fraction, b: Fraction, A: X.Func, x, f0=None) -> mp.mpf:
    lc = _real(lc_real(a, b), "LC")
    s = mp.mpf(0)
    for k, d in A.items():
        if k[0] == "f":
            base = f0
        elif k[0] == "E":
            base = mp.exp(mpq(k[1]) * x) * lc ** k[2]
        else:
            assert k[2] != 0, k
            base = mp.exp(mpq(k[1]) * x) * mp.ei(mpq(k[2]) * x)
        _real(base, f"base value of key {k}")
        for p, c in d.items():
            s += mpq(c) * x ** p * base
    return _real(s, "sum")


def f0_real(a: Fraction, b: Fraction, x) -> mp.mpf:
    fm1 = X.Calc(X.Ctx("Q", Fraction(a), Fraction(b))).FM1
    return mp.quad(lambda s: evaluate_real(a, b, fm1, s) / s, [x, x + 1, x + 2, x + 4, x + 10, x + 30, mp.inf])


def value_func(a, b, K: X.Func, R, dps: int = 60) -> mp.mpf:
    a, b = Fraction(a), Fraction(b)
    with mp.workdps(dps):
        x = mpq(R)
        f0 = f0_real(a, b, x) if ("f",) in K else None
        return +evaluate_real(a, b, K, x, f0)


def value_real(case: str, a, b, R, dps: int = 60) -> mp.mpf:
    """The stored closed form of `case` at (a, b), real branch, at R."""
    return value_func(a, b, coords_to_func(a, b, stored_coords(case, a, b)), R, dps)


# ------------------------------------------------------------------ typed evaluator, complex t allowed
@lru_cache(maxsize=None)
def coeff_polys(case: str) -> Tuple[Dict[int, Tuple[Tuple[int, ...], Tuple[int, ...]]], ...]:
    """Per function (NAMES9 order): {power of R: (numerator, denominator) integer coefficient lists in t}."""
    t = sp.Symbol("t")
    out = []
    for nm in X.NAMES9:
        d = {}
        for p, s in stored(case).get(nm, {}).items():
            n, dd = sp.fraction(sp.together(sp.sympify(s, locals={"t": t})))
            d[int(p)] = (tuple(int(c) for c in sp.Poly(n, t).all_coeffs()),
                         tuple(int(c) for c in sp.Poly(dd, t).all_coeffs()))
        out.append(d)
    return tuple(out)


def ei_from_positive(w):
    """The continuation of Ei from w > 0: gamma + Log w - Ein(-w), Ein(z) = E1(z) + Log z + gamma."""
    return mp.log(w) - mp.log(-w) - mp.e1(-w)


def ei_from_negative(w):
    """The continuation of Ei from w < 0: -E1(-w)."""
    return -mp.e1(-w)


def lc_cont(t):
    return mp.log(1 + t) - mp.log(t - 1)


def nine_typed(t, R, branch: str):
    """[None, g2, ..., g9] at (a, b) = (1, t), and F^(-1) = c (psi_a - psi_b + g7 - g6), c = 1/(2(1 - t^2))."""
    if branch == "real":
        LC = mp.log((1 + t) / abs(t - 1))
        Ep = Em = mp.ei
    elif branch == "cont":
        LC = lc_cont(t)
        Ep, Em = ei_from_positive, ei_from_negative
    else:
        raise ValueError(branch)
    g = [None] * 9
    g[1] = Em(-(1 + t) * R)
    g[2] = mp.exp(-(1 + t) * R)
    g[3] = mp.exp(-2 * R)
    g[4] = mp.exp(-2 * t * R)
    E2 = Em(-2 * (1 + t) * R)
    g[5] = mp.exp(2 * R) * E2
    g[6] = mp.exp(2 * t * R) * E2
    # at t > 1 the psi_a arguments (1 - t) R are negative and the psi_b arguments (t - 1) R positive
    g[7] = mp.exp(-2 * R) * (Em(2 * (1 - t) * R) - 2 * Em((1 - t) * R) - LC)
    g[8] = mp.exp(-2 * t * R) * (Ep(2 * (t - 1) * R) - 2 * Ep((t - 1) * R) - LC)
    c = 1 / (2 * (1 - t * t))
    return g, c * (g[7] - g[8] + g[6] - g[5])


def f0_typed(t, R, branch: str):
    return mp.quad(lambda s: nine_typed(t, s, branch)[1] / s, [R, R + 1, R + 2, R + 4, R + 10, R + 30, mp.inf])


def K1_typed(case: str, t, R, branch: str = "real"):
    """K(R; 1, t) from the stored coefficients and the typed nine functions."""
    g, _ = nine_typed(t, R, branch)
    g[0] = f0_typed(t, R, branch)
    s = 0
    for i, d in enumerate(coeff_polys(case)):
        for p, (n, dd) in d.items():
            s += mp.polyval(list(n), t) / mp.polyval(list(dd), t) * R ** p * g[i]
    return s


# ------------------------------------------------------------------ fresh exact solves, every solve re-checked
@contextlib.contextmanager
def checked_solves():
    """Wrap X.Calc.solve_L0 and M.MCalc.solve_D2: re-check each output as an exact identity (L0 u = src, resp.
    u'' = src), record its free unknowns and its keys.  Restores the originals on exit."""
    stats = {"solves": 0, "free": [], "residual_failures": 0, "keys": set()}
    o_l0, o_d2 = X.Calc.solve_L0, M.MCalc.solve_D2

    def l0(self, src):
        sol = o_l0(self, src)
        stats["solves"] += 1
        stats["free"].append(self.last_free)
        if X.fadd(self.L0(sol), src, -1) != {}:
            stats["residual_failures"] += 1
        stats["keys"].update(sol)
        return sol

    def d2(self, src):
        sol = o_d2(self, src)
        stats["solves"] += 1
        stats["free"].append(self.last_free)
        if X.fadd(self.D2(sol), src, -1) != {}:
            stats["residual_failures"] += 1
        stats["keys"].update(sol)
        return sol

    X.Calc.solve_L0, M.MCalc.solve_D2 = l0, d2
    try:
        yield stats
    finally:
        X.Calc.solve_L0, M.MCalc.solve_D2 = o_l0, o_d2


def fresh_solve(case: str, a: Fraction, b: Fraction):
    """(K over Q, its coordinates on the nine functions or None, stats) from a fresh exact solve at (a, b)."""
    a, b = Fraction(a), Fraction(b)
    ctx = X.Ctx("Q", a, b)
    with checked_solves() as st:
        if case in M0_CASES:
            calc, fS, Fs, P2, nfree = X.solve_all(ctx, X.case_downset(case))
            K = X.assemble_K(ctx, case, fS)
            st["returned_free"] = nfree
        else:
            calc, fl, nf0, free = M.solve_case(ctx, case)
            K = M.assemble_K(ctx, case, fl)
            st["returned_free"] = nf0 + sum(free)
    return K, X.in_nine(ctx, K), st


def key_rate(k):
    """Growth rate at R -> oo of a key over Q (None for f0).  Ei(mu R) ~ e^{mu R}/(mu R) for mu of either sign
    (principal value for mu > 0), so the rate of e^{lam R} Ei(mu R) is lam + mu."""
    if k[0] == "f":
        return None
    if k[0] == "E":
        return Fraction(k[1])
    assert k[2] != 0, k
    return Fraction(k[1]) + Fraction(k[2])


# ------------------------------------------------------------------ root counting for leg (D)
def roots_in(g, lo, hi) -> int:
    """Distinct real roots of the sympy Poly g in the open interval (lo, hi); hi may be sp.oo."""
    if g.is_zero:
        raise AssertionError("zero polynomial")
    if g.degree() <= 0:
        return 0
    n = int(g.count_roots(lo, hi))
    for e in (lo, hi):
        if e != sp.oo and g.eval(e) == 0:
            n -= 1
    return n
