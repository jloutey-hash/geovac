"""Paper 58 sec:census, one-centre classes: each integral is ONE rational times ONE square root.

The [SYMBOLIC PROOF] there: (ab|cd) = sum_k R^k c^k(a,b) c^k(d,c); in Racah's formula the only k-dependent root
of a Gaunt factor is sqrt((k-q)! (k+q)!), q = m_a - m_b, and the partner factor carries the same q, so it squares
out and every k-term lies in one square class.  The wrong answer this guard rejects: an integral whose k-terms
carry two different square classes (a sum like q1 + q2 sqrt(s), not q sqrt(s)).  Checked exactly (sympy 3j, no
floats) over every complex-harmonic pattern with l <= 3, and for real harmonics at l <= 2.

Fire-tested 2026-11-16 (debug/qa/fire_test.py): the partner factor taken as a Gaunt integral, sqrt(2k+1) off c^k,
FIRED; a classifier that sees one class FIRED.  A sqrt(k!) on every c^k did not fire, correctly: it squares out
in the product, which is the mechanism, so it is not a wrong answer.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

import sympy as sp
from sympy.physics.wigner import wigner_3j

LMAX_COMPLEX = 3
LMAX_REAL = 2


@lru_cache(maxsize=None)
def ck(l1: int, m1: int, l2: int, m2: int, k: int) -> sp.Expr:
    """c^k(l1 m1, l2 m2) = (-1)^m1 sqrt((2l1+1)(2l2+1)) (l1 k l2; 0 0 0)(l1 k l2; -m1 q m2), q = m1 - m2."""
    q = m1 - m2
    if abs(q) > k or (l1 + l2 + k) % 2:
        return sp.Integer(0)
    return sp.radsimp(sp.Integer(-1) ** m1 * sp.sqrt((2 * l1 + 1) * (2 * l2 + 1))
                      * wigner_3j(l1, k, l2, 0, 0, 0) * wigner_3j(l1, k, l2, -m1, q, m2))


def square_class(v: sp.Expr) -> int:
    """square-free part of the rational v^2 (v a rational times one square root)."""
    r = sp.Rational(sp.radsimp(sp.expand(v ** 2)))
    s = 1
    for pr, e in sp.factorint(r.p * r.q).items():
        if e % 2:
            s *= pr
    return s


def classes_of(terms: dict) -> set:
    return {square_class(t) for t in terms.values() if t != 0}


def k_terms(a, b, c, d) -> dict:
    if a[1] + c[1] != b[1] + d[1]:
        return {}
    kmax = min(a[0] + b[0], c[0] + d[0])
    return {k: ck(a[0], a[1], b[0], b[1], k) * ck(d[0], d[1], c[0], c[1], k) for k in range(kmax + 1)}


def _orbs(lmax: int):
    return [(l, m) for l in range(lmax + 1) for m in range(-l, l + 1)]


def test_classifier_rejects_a_two_class_sum():
    """negative control: the instrument sees a planted q1 sqrt(2) + q2 sqrt(3) as two classes."""
    assert classes_of({0: sp.Rational(1, 3) * sp.sqrt(2), 2: sp.Rational(2, 5) * sp.sqrt(3)}) == {2, 3}
    assert classes_of({0: sp.Rational(1, 3) * sp.sqrt(8), 2: sp.Rational(2, 5) * sp.sqrt(18)}) == {2}


def test_gaunt_factor_is_irrational_with_k_dependent_root():
    """the cancellation is not vacuous: single Gaunt factors carry a k-dependent root that only the product kills."""
    single = {k: ck(2, 0, 2, 1, k) for k in (2, 4)}
    assert all(v != 0 for v in single.values())
    assert len(classes_of(single)) == 2
    assert len(classes_of(k_terms((2, 0), (2, 1), (2, 1), (2, 0)))) == 1


def test_complex_harmonics_single_square_class_l_le_3():
    n_nonzero, mixed = 0, []
    for a, b, c, d in product(_orbs(LMAX_COMPLEX), repeat=4):
        t = k_terms(a, b, c, d)
        cl = classes_of(t)
        if cl:
            n_nonzero += 1
            if len(cl) > 1:
                mixed.append((a, b, c, d, cl))
    assert n_nonzero == 3972       # every nonzero complex pattern, l <= 3 (counted 2026-11-16)
    assert mixed == [], mixed[:5]


def _real(l: int, m: int):
    if m == 0:
        return [((l, 0), sp.Integer(1))]
    a = abs(m)
    if m > 0:
        return [((l, a), (-1) ** a / sp.sqrt(2)), ((l, -a), 1 / sp.sqrt(2))]
    return [((l, -a), sp.I / sp.sqrt(2)), ((l, a), -sp.I * (-1) ** a / sp.sqrt(2))]


def test_real_harmonics_single_square_class_l_le_2():
    n_nonzero, mixed = 0, []
    for A, B, C, D in product(_orbs(LMAX_REAL), repeat=4):
        acc: dict = {}
        for (a, ca), (b, cb), (c, cc), (d, cd) in product(_real(*A), _real(*B), _real(*C), _real(*D)):
            for k, t in k_terms(a, b, c, d).items():
                if t != 0:
                    acc[k] = acc.get(k, 0) + sp.conjugate(ca) * cb * sp.conjugate(cc) * cd * t
        acc = {k: sp.radsimp(sp.expand(v)) for k, v in acc.items()}
        acc = {k: v for k, v in acc.items() if v != 0}
        if acc:
            n_nonzero += 1
            assert all(sp.im(v) == 0 for v in acc.values()), (A, B, C, D, acc)
            cl = classes_of(acc)
            if len(cl) > 1:
                mixed.append((A, B, C, D, acc))
    assert n_nonzero == 753        # every nonzero real pattern, l <= 2 (counted 2026-11-16)
    assert mixed == [], mixed[:5]
