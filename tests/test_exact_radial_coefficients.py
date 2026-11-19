"""Regression guard: exact radial coefficients (no nsimplify) in
geovac.two_center_eri.radial_product and geovac.qfd_core.s_radial.

sp.nsimplify applied to an exact surd sometimes returns a different surd (relative error 1e-18 .. 1.7e-17), e.g.
s_radial at Z = 19/2, 17/4, 19/3 and radial_product for pairs involving Z = 19/10 or 99/100.  Both now use
sp.sympify.  The same fix in radial_norm is guarded in tests/test_paper58_exchange_reverse.py (F1-F3); not repeated.

Reference construction: the textbook N^2 (a Rational, independent of the code under test) times the exact
polynomial coefficients.  A value x is accepted iff sign(x) == sign(ref_signed) and x**2 == ref**2 exactly
(both sides Rational), so a wrong surd (even at 1e-18) is rejected.  No float tolerance anywhere.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import factorial

import pytest
import sympy as sp

from geovac import qfd_core as Q
from geovac import two_center_eri as TCE


def _norm2(Z: F, n: int, l: int):
    """Textbook N^2 for radial_poly's convention (lam^(2l+3) (n-l-1)!/(2n (n+l)!))."""
    lam = sp.Rational(2 * Z.numerator, Z.denominator * n)
    return lam ** (2 * l + 3) * sp.Rational(factorial(n - l - 1), 2 * n * factorial(n + l))


def _same(x, ref2, ref_sign):
    """x exactly equals the number with square ref2 and sign ref_sign (ref2 Rational)."""
    x = sp.sympify(x)
    if ref2 == 0:
        return x == 0
    return (x ** 2).is_Rational and x ** 2 == ref2 and sp.sign(x) == ref_sign


SWEEP_S = sorted({F(p, q) for q in range(1, 21) for p in range(1, 21)} | {F(99, 100), F(199, 200), F(97, 100)})
PAIR_Z = [F(19, 10), F(99, 100), F(19, 2), F(17, 4), F(19, 3), F(1), F(2), F(3, 2), F(31, 15), F(40, 11)]
PAIR_NL = [(1, 0), (2, 0), (2, 1)]
PAIR_NL_SLOW = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]


def _old_s_radial(Z, n):
    Zf = Z if isinstance(Z, F) else F(Z)
    c, a = TCE.radial_poly(Zf, n, 0)
    N = TCE.radial_norm(Zf, n, 0)
    return tuple(sorted((k, sp.nsimplify(N * v)) for k, v in c.items())), a


def _old_radial_product(Z1, n1, l1, Z2, n2, l2):
    c1, a1 = TCE.radial_poly(Z1, n1, l1)
    c2, a2 = TCE.radial_poly(Z2, n2, l2)
    N1, N2 = TCE.radial_norm(Z1, n1, l1), TCE.radial_norm(Z2, n2, l2)
    out = {}
    for k1, v1 in c1.items():
        for k2, v2 in c2.items():
            out[k1 + k2] = out.get(k1 + k2, 0) + N1 * N2 * v1 * v2
    return {k: sp.nsimplify(v) for k, v in out.items()}, a1 + a2


def _bad_s_radial(fn, zs):
    bad = []
    for Z in zs:
        for n in (1, 2, 3):
            coeffs, a = fn(Z, n)
            c, a_ref = TCE.radial_poly(Z, n, 0)
            if a != a_ref or [k for k, _ in coeffs] != sorted(c):
                bad.append((str(Z), n, "structure"))
                continue
            N2 = _norm2(Z, n, 0)
            for k, x in coeffs:
                v = c[k]
                if not _same(x, N2 * v ** 2, sp.sign(v)):
                    bad.append((str(Z), n, k))
    return bad


def _bad_product(fn, zs, nls):
    bad = []
    for Z1 in zs:
        for Z2 in zs:
            for n1, l1 in nls:
                for n2, l2 in nls:
                    out, a = fn(Z1, n1, l1, Z2, n2, l2)
                    c1, a1 = TCE.radial_poly(Z1, n1, l1)
                    c2, a2 = TCE.radial_poly(Z2, n2, l2)
                    ref = {}
                    for k1, v1 in c1.items():
                        for k2, v2 in c2.items():
                            ref[k1 + k2] = ref.get(k1 + k2, 0) + v1 * v2
                    n12 = _norm2(Z1, n1, l1) * _norm2(Z2, n2, l2)
                    if a != a1 + a2 or sorted(out) != sorted(ref):
                        bad.append((str(Z1), n1, l1, str(Z2), n2, l2, "structure"))
                        continue
                    for k, s in ref.items():
                        if not _same(out[k], n12 * s ** 2, sp.sign(s)):
                            bad.append((str(Z1), n1, l1, str(Z2), n2, l2, k))
    return bad


def test_s_radial_sweep_is_exact():
    """G1.  s_radial(Z, n) coefficients are exactly N * (radial_poly coefficient): coefficient^2 is the Rational
    N^2 v^2 (textbook N^2) with the right sign, for 258 distinct rational Z (p/q, p, q <= 20, plus 99/100, 199/200,
    97/100) and n <= 3.

    Rejects: s_radial built with sp.nsimplify (wrong surd at Z = 19/2, 17/4, 19/3, ...); a dropped or doubled
    normalisation; a wrong sign or exponent set.
    """
    assert len(SWEEP_S) > 250
    assert _bad_s_radial(Q.s_radial, SWEEP_S) == []


def test_radial_product_pairs_are_exact():
    """G2.  radial_product coefficients equal the exact double sum N1 N2 v1 v2 (square compared as a Rational)
    for all ordered pairs of 10 Z values (incl. 19/10, 99/100) and (n, l) in {1s, 2s, 2p}: 900 cases.

    Rejects: radial_product with sp.nsimplify on the coefficients (wrong surd for pairs with Z = 19/10 or 99/100);
    a dropped normalisation; a wrong combined rate or power set.
    """
    assert _bad_product(TCE.radial_product, PAIR_Z, PAIR_NL) == []


@pytest.mark.slow
def test_radial_product_pairs_are_exact_large():
    """G2-slow.  As G2 over n <= 3 (all l) and the wider Z list."""
    assert _bad_product(TCE.radial_product, PAIR_Z + [F(37, 3), F(7, 5)], PAIR_NL_SLOW) == []


def test_nsimplify_versions_would_be_caught():
    """G3.  Positive control: the pre-fix nsimplify versions, rebuilt inline, fail G1 and G2's checks on at least
    one sweep case each.  Skipped (not failed) if a future sympy makes nsimplify exact on all of them.

    Rejects: a G1/G2 that cannot tell the old and new versions apart.
    """
    old_s = _bad_s_radial(_old_s_radial, [F(19, 2), F(17, 4), F(19, 3)])
    old_p = _bad_product(_old_radial_product, [F(19, 10), F(99, 100)], PAIR_NL)
    if not old_s and not old_p:
        pytest.skip("this sympy's nsimplify is exact on all of these inputs")
    assert old_s, "nsimplify s_radial not caught at Z = 19/2, 17/4, 19/3"
    assert old_p, "nsimplify radial_product not caught at Z = 19/10, 99/100"


# FIRE-TEST (2026-11-19; planted by monkeypatching the old nsimplify versions, geovac/ untouched)
#   plant                                              G1 s_radial   G2 radial_product   G3 control
#   qfd_core.s_radial   -> nsimplify(N*v)              FIRED         (n/a)               passes
#   two_center_eri.radial_product -> nsimplify(v)      (n/a)         FIRED               passes
# Clean: fast 3 passed, 1 skipped (slow G2) in ~8 s; --slow 4 passed in ~20 s.
