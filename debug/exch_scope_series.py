"""Scoping (row 3c): exact log-Taylor series at R=0 of the 1s exchange integral, and an
operator guesser mod p.  Control only at equal zeta (OP_EXCH_1 reproduction); at unequal
zeta this script ONLY times the series build and prints nothing about its content.

Representation: f(R) = sum_{const monomial m, j in {0,1}} m * (ln R)^j * F_{m,j}(R),
F_{m,j} Laurent series with Fraction coefficients.  Constants (gamma_E, ln of primes) are
kept as independent sympy symbols, as in tests/test_paper58_r_operators.py.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction
from math import factorial
from typing import Dict, List, Tuple

import sympy as sp

sys.path.insert(0, "."); sys.path.insert(0, "tests")
from geovac import qfd_core as Q  # noqa: E402
from geovac.two_center_eri import R_s  # noqa: E402
from test_paper58_r_operators import OP_EXCH_1, decompose  # noqa: E402

R = R_s
S1 = (1, 0, 0)
G = sp.Symbol("gE")

Series = Dict[Tuple[sp.Expr, int], Dict[int, Fraction]]


def _add(out: Series, key, k: int, v: Fraction) -> None:
    if v == 0:
        return
    d = out.setdefault(key, {})
    d[k] = d.get(k, Fraction(0)) + v


def _poly_coeffs(p: sp.Expr) -> Dict[int, Fraction]:
    P = sp.Poly(sp.expand(p), R)
    return {m[0]: Fraction(int(c.p), int(c.q)) for m, c in P.terms()}


def _inv_series(den: Dict[int, Fraction], N: int) -> Tuple[int, List[Fraction]]:
    v = min(den)
    d = [den.get(v + i, Fraction(0)) for i in range(N + 1)]
    inv = [Fraction(0)] * (N + 1)
    inv[0] = 1 / d[0]
    for n in range(1, N + 1):
        s = sum(d[i] * inv[n - i] for i in range(1, n + 1))
        inv[n] = -s / d[0]
    return v, inv


def _ln_number(k: sp.Expr) -> sp.Expr:
    return sp.expand_log(sp.log(k), force=True)


def expand_atom(c, kind, k, coeff: sp.Expr, reps, N: int) -> Series:
    """Series of coeff(R, consts) * e^{cR} * kind(kR) to R^N (relative order N)."""
    inv_reps = {v: kk for kk, v in reps.items()}
    coeff = sp.together(coeff)
    num, den = sp.fraction(coeff)
    csyms = sorted(set(reps.values()), key=str)
    numP = sp.Poly(sp.expand(num), *csyms) if csyms else None
    num_parts = ({m: cf for m, cf in numP.terms()} if csyms else {(): num})
    dshift, dinv = _inv_series(_poly_coeffs(den), N + 20)
    ec = [Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)) ** n / factorial(n)
          for n in range(N + 30)]
    kq = Fraction(int(sp.Rational(k).p), int(sp.Rational(k).q)) if k != 0 else Fraction(0)
    # kind series: list of (const expr, logpow, {power: coeff})
    if kind == "1":
        ks = [(sp.Integer(1), 0, {0: Fraction(1)})]
    elif kind == "log":
        ks = [(sp.Integer(1), 1, {0: Fraction(1)})]
    elif kind == "E1":
        ent = {n: -(-kq) ** n / (n * factorial(n)) for n in range(1, N + 30)}
        ks = [(G, 0, {0: Fraction(-1)}), (_ln_number(sp.Rational(k)), 0, {0: Fraction(-1)}),
              (sp.Integer(1), 1, {0: Fraction(-1)}), (sp.Integer(1), 0, ent)]
    elif kind == "Ei":
        ent = {n: kq ** n / (n * factorial(n)) for n in range(1, N + 30)}
        ks = [(G, 0, {0: Fraction(1)}), (_ln_number(abs(sp.Rational(k))), 0, {0: Fraction(1)}),
              (sp.Integer(1), 1, {0: Fraction(1)}), (sp.Integer(1), 0, ent)]
    else:
        raise ValueError(kind)
    out: Series = {}
    for mono, cf in num_parts.items():
        cexpr = sp.Mul(*[s ** e for s, e in zip(csyms, mono)]) if csyms else sp.Integer(1)
        cexpr = sp.expand(cexpr.xreplace(inv_reps).xreplace({sp.EulerGamma: G}))
        npoly = _poly_coeffs(cf)
        # rational part: num * 1/den, Laurent starting at -dshift
        rat: Dict[int, Fraction] = {}
        for a, va in npoly.items():
            for i, vi in enumerate(dinv):
                if vi:
                    rat[a + i - dshift] = rat.get(a + i - dshift, Fraction(0)) + va * vi
        lo = min(rat)
        # times e^{cR}
        re: Dict[int, Fraction] = {}
        for a, va in rat.items():
            for n in range(0, N + 25):
                if a + n > lo + N + 10:
                    break
                re[a + n] = re.get(a + n, Fraction(0)) + va * ec[n]
        for kc, lp, kser in ks:
            for ce, cv in sp.expand(cexpr * kc).as_coefficients_dict().items():
                key = (ce, lp)
                for a, va in re.items():
                    for b, vb in kser.items():
                        if a + b <= lo + N + 10:
                            _add(out, key, a + b, va * vb * Fraction(int(cv.p), int(cv.q)))
    return out


def series_of(expr: sp.Expr, N: int) -> Series:
    atoms, reps = decompose(expr)
    tot: Series = {}
    for (c, kind, k), v in atoms.items():
        for key, d in expand_atom(c, kind, k, v, reps, N).items():
            for a, val in d.items():
                _add(tot, key, a, val)
    return {k: {a: x for a, x in d.items() if x != 0} for k, d in tot.items()
            if any(x != 0 for x in d.values())}


# ------------------------------------------------------------------ guesser mod p
def guess_mod_p(ser: Series, r: int, d: int, p: int, mmax: int) -> List[List[List[int]]]:
    """Nullspace of L = sum c_ij R^j D^i acting on every component, mod p.
    Equations: coefficient of m * lnR^lp * R^t in L f, for t <= mmax."""
    comps: Dict[sp.Expr, Dict[int, Dict[int, Fraction]]] = {}
    for (ce, lp), dd in ser.items():
        comps.setdefault(ce, {})[lp] = dd
    unk = [(i, j) for i in range(r + 1) for j in range(d + 1)]
    rows = []
    for ce, parts in comps.items():
        F0, F1 = parts.get(0, {}), parts.get(1, {})
        # D^i of R^k ln R = R^{k-i}( ff(k,i) ln R + dff(k,i) )
        eq: Dict[Tuple[int, int], Dict[int, Fraction]] = {}
        for idx, (i, j) in enumerate(unk):
            for k, v in F0.items():
                ff = 1
                for s in range(i):
                    ff *= (k - s)
                if ff:
                    eq.setdefault((0, k - i + j), {})
                    eq[(0, k - i + j)][idx] = eq[(0, k - i + j)].get(idx, 0) + v * ff
            for k, v in F1.items():
                # d^i/dR^i (R^k lnR) = R^{k-i}[ P_i(k) lnR + P_i'(k) ]
                ff = 1
                for s in range(i):
                    ff *= (k - s)
                dff = sp.diff(sp.prod([sp.Symbol("x") - s for s in range(i)]), sp.Symbol("x")).subs(sp.Symbol("x"), k) if i else 0
                dff = int(dff)
                if ff:
                    eq.setdefault((1, k - i + j), {})
                    eq[(1, k - i + j)][idx] = eq[(1, k - i + j)].get(idx, 0) + v * ff
                if dff:
                    eq.setdefault((0, k - i + j), {})
                    eq[(0, k - i + j)][idx] = eq[(0, k - i + j)].get(idx, 0) + v * dff
        for (lp, t), row in eq.items():
            if t <= mmax:
                rows.append([int((row.get(u, 0)).numerator * pow((row.get(u, Fraction(0))).denominator, -1, p)) % p
                             if isinstance(row.get(u, 0), Fraction) else int(row.get(u, 0)) % p
                             for u in range(len(unk))])
    ns = nullspace_mod_p(rows, len(unk), p)
    return [[[vec[unk.index((i, j))] for j in range(d + 1)] for i in range(r + 1)] for vec in ns]


def nullspace_mod_p(rows: List[List[int]], n: int, p: int) -> List[List[int]]:
    M = [row[:] for row in rows]
    piv_cols, rr = [], 0
    for c in range(n):
        pr = next((i for i in range(rr, len(M)) if M[i][c] % p), None)
        if pr is None:
            continue
        M[rr], M[pr] = M[pr], M[rr]
        inv = pow(M[rr][c], -1, p)
        M[rr] = [x * inv % p for x in M[rr]]
        for i in range(len(M)):
            if i != rr and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rr])]
        piv_cols.append(c)
        rr += 1
    free = [c for c in range(n) if c not in piv_cols]
    basis = []
    for fc in free:
        v = [0] * n
        v[fc] = 1
        for i, pc in enumerate(piv_cols):
            v[pc] = (-M[i][fc]) % p
        basis.append(v)
    return basis


def op_mod_p(strs, p, d):
    out = []
    for s in strs:
        P = sp.Poly(sp.expand(sp.sympify(s, locals={"R": R})), R)
        cs = [0] * (d + 1)
        for (m,), c in P.terms():
            cs[m] = int(c) % p
        out.append(cs)
    return out


def proportional(a, b, p) -> bool:
    fa = [x for row in a for x in row]
    fb = [x for row in b for x in row]
    i = next(k for k, x in enumerate(fb) if x)
    lam = fa[i] * pow(fb[i], -1, p) % p
    return all((x - lam * y) % p == 0 for x, y in zip(fa, fb))


def main() -> None:
    p = 1000003
    # ---- control: equal zeta = 1
    t0 = time.time()
    f = Q.exchange_closed_form(Fraction(1), S1, S1, Fraction(1), S1, S1, R, tau_max=6)
    t1 = time.time()
    N = 70
    ser = series_of(f, N)
    t2 = time.time()
    print(f"[control zeta=1] closed form {t1-t0:.1f}s, series to R^{N} {t2-t1:.1f}s, "
          f"components {len(ser)}")
    for r, d in [(3, 8), (3, 12), (4, 7), (4, 8), (4, 10)]:
        ns = guess_mod_p(ser, r, d, p, mmax=N - 15)
        msg = f"  guess r={r} d={d}: nullity {len(ns)}"
        if r == 4 and d == 8 and len(ns) == 1:
            msg += f"; proportional to OP_EXCH_1 mod p: {proportional(ns[0], op_mod_p(OP_EXCH_1, p, 8), p)}"
        print(msg)
    # ---- cost only at unequal zeta (content NOT printed or inspected)
    for tm in (4, 10):
        t0 = time.time()
        terms = Q.exchange_closed_form(Fraction(1), S1, S1, Fraction(2), S1, S1, R,
                                       tau_max=tm, return_terms=True)
        t1 = time.time()
        _ = [series_of(T, 20) for T in terms if T != 0]
        t2 = time.time()
        print(f"[cost zeta=(1,2)] tau_max={tm}: closed forms {t1-t0:.1f}s, "
              f"series to rel. order 20 {t2-t1:.1f}s")


if __name__ == "__main__":
    main()
