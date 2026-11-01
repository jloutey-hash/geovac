"""Second route for record HG's occupation weights: Mardoyan-Pogosyan-Sissakian-Ter-Antonyan
(arXiv:2310.17336), eq. (2.2.29), the Tarter 3F2 form of the spherical<->parabolic coefficients.

  W^{nl}_{n1 n2 m} = (-1)^{l+(m-|m|)/2}
       sqrt[ (2l+1)(n1+|m|)!(n2+|m|)!(l+|m|)! / (n1! n2! (l-|m|)! (n-l-1)! (n+l)!) ]
       * (n-|m|-1)!/|m|! * 3F2(-n2, -l+|m|, l+|m|+1 ; |m|+1, -n+|m|+1 ; 1),
  n = n1 + n2 + |m| + 1.

No Clebsch-Gordan routine is used for W.  W^2 is exact rational.  Label map (from eq. 2.2.25, whose CG has
projections (n1-n2+|m|)/2 on the first j and (n2-n1+|m|)/2 on the second):  mu = (m - n1 + n2)/2, i.e. the
parabolic state (n1, n2, m) is the J-_z = mu state of record HG (A_z eigenvalue n1 - n2 in the CG ordering).

Comparison target: the record's occupation |<j, m-mu; j, mu | l, m>|^2 (sympy CG, as in
tests/test_paper25_hopf_grouping.py::hopf_M).

Usage: python debug/hg_mardoyan_second_route.py
"""
from fractions import Fraction
from math import factorial as f
from typing import Callable, Dict, Tuple

import sympy as sp
from sympy.physics.quantum.cg import CG


def poch(a: int, k: int) -> int:
    p = 1
    for i in range(k):
        p *= a + i
    return p


def hyp3f2(a: Tuple[int, int, int], b: Tuple[int, int], kmax: int) -> Fraction:
    s = Fraction(0)
    for k in range(kmax + 1):
        den = poch(b[0], k) * poch(b[1], k) * f(k)
        if den == 0:
            break
        s += Fraction(poch(a[0], k) * poch(a[1], k) * poch(a[2], k), den)
    return s


def W2_book(n1: int, n2: int, m: int, l: int, lower_shift: int = 0) -> Fraction:
    """Square of eq. (2.2.29). lower_shift plants a typo in the second lower 3F2 parameter (control)."""
    am = abs(m)
    n = n1 + n2 + am + 1
    if not (am <= l <= n - 1):
        return Fraction(0)
    pref = Fraction((2 * l + 1) * f(n1 + am) * f(n2 + am) * f(l + am),
                    f(n1) * f(n2) * f(l - am) * f(n - l - 1) * f(n + l))
    F = hyp3f2((-n2, -l + am, l + am + 1), (am + 1, -n + am + 1 + lower_shift), n2)
    return pref * Fraction(f(n - am - 1), f(am)) ** 2 * F * F


def W2_book_n1(n1: int, n2: int, m: int, l: int) -> Fraction:
    """Control: the 3F2's first upper parameter -n2 replaced by -n1, prefactor kept (a transcription slip)."""
    am = abs(m)
    n = n1 + n2 + am + 1
    pref = Fraction((2 * l + 1) * f(n1 + am) * f(n2 + am) * f(l + am),
                    f(n1) * f(n2) * f(l - am) * f(n - l - 1) * f(n + l))
    F = hyp3f2((-n1, -l + am, l + am + 1), (am + 1, -n + am + 1), n1)
    return pref * Fraction(f(n - am - 1), f(am)) ** 2 * F * F


def mu_true(n1: int, n2: int, m: int) -> Fraction:
    return Fraction(m - n1 + n2, 2)


def book_occupations(nmax: int, mu_map: Callable[[int, int, int], Fraction] = mu_true,
                     lower_shift: int = 0, use_n1: bool = False) -> Dict[Tuple[int, int, int, Fraction], Fraction]:
    occ: Dict[Tuple[int, int, int, Fraction], Fraction] = {}
    for n in range(1, nmax + 1):
        for m in range(-(n - 1), n):
            for n1 in range(0, n - abs(m)):
                n2 = n - 1 - abs(m) - n1
                mu = mu_map(n1, n2, m)
                for l in range(abs(m), n):
                    key = (n, l, m, mu)
                    occ[key] = occ.get(key, Fraction(0)) + (
                        W2_book_n1(n1, n2, m, l) if use_n1 else W2_book(n1, n2, m, l, lower_shift))
    return occ


def record_occupations(nmax: int) -> Dict[Tuple[int, int, int, Fraction], Fraction]:
    occ = {}
    for n in range(1, nmax + 1):
        j = sp.Rational(n - 1, 2)
        for l in range(n):
            for m in range(-l, l + 1):
                for k in range(n):
                    mu = sp.Rational(2 * k - (n - 1), 2)
                    if abs(m - mu) > j:
                        continue
                    c2 = sp.nsimplify(CG(j, m - mu, j, mu, l, m).doit() ** 2)
                    occ[(n, l, m, Fraction(int(mu * 2), 2))] = Fraction(int(c2.p), int(c2.q))
    return occ


def compare(a: Dict, b: Dict) -> Tuple[Fraction, object]:
    worst, where = Fraction(0), None
    for key in set(a) | set(b):
        d = abs(a.get(key, Fraction(0)) - b.get(key, Fraction(0)))
        if d > worst:
            worst, where = d, key
    return worst, where


def main() -> None:
    # self-checks of the transcribed formula: orthonormality of the book matrix within each (n, m)
    for n in range(1, 5):
        for m in range(-(n - 1), n):
            for n1 in range(0, n - abs(m)):
                n2 = n - 1 - abs(m) - n1
                assert sum(W2_book(n1, n2, m, l) for l in range(abs(m), n)) == 1, (n, m, n1)
            for l in range(abs(m), n):
                assert sum(W2_book(n1, n - 1 - abs(m) - n1, m, l) for n1 in range(0, n - abs(m))) == 1
    print("book W^2: rows and columns sum to 1 exactly for every (n, m), n <= 4")

    rec4 = record_occupations(4)
    for nmax in range(1, 5):
        rec = {k: v for k, v in rec4.items() if k[0] <= nmax}
        bk = book_occupations(nmax)
        d, w = compare(bk, rec)
        nz = sum(1 for v in rec.values() if v != 0)
        print(f"n_max={nmax}: entries {len(rec)} (nonzero {nz}), max |book - record| = {d} at {w}")

    print("\ncontrols (n_max = 4):")
    swap = lambda n1, n2, m: Fraction(m - n2 + n1, 2)  # n1 <-> n2 at fixed m
    ctrls = {
        "C-a n1<->n2 swapped at fixed m": dict(mu_map=swap),
        "C-b mu = (n2 - n1)/2 (m dropped)": dict(mu_map=lambda n1, n2, m: Fraction(n2 - n1, 2)),
        "C-c mu = (m - n2 + n1)/2 for m < 0 only (|m| map)": dict(
            mu_map=lambda n1, n2, m: mu_true(n1, n2, m) if m >= 0 else swap(n1, n2, m)),
        "C-e mu -> -mu at fixed m (L_z sign error)": dict(mu_map=lambda n1, n2, m: -mu_true(n1, n2, m)),
        "C-f n1 <-> n2 in the 3F2 only (-n2 -> -n1)": dict(mu_map=swap, lower_shift=0, use_n1=True),
        "C-d 3F2 lower parameter -n+|m|+1 -> -n+|m|+2": dict(lower_shift=1),
    }
    for name, kw in ctrls.items():
        try:
            d, w = compare(book_occupations(4, **kw), rec4)
            print(f"  {name}: max diff {d} at {w}")
        except ZeroDivisionError:
            print(f"  {name}: ZeroDivisionError")
    # why C-a cannot discriminate: CG^2 is symmetric under m1 <-> m2 when j1 = j2
    sym = all(CG(sp.Rational(n - 1, 2), a, sp.Rational(n - 1, 2), sp.Rational(m, 1) - a, l, m).doit() ** 2
              == CG(sp.Rational(n - 1, 2), sp.Rational(m, 1) - a, sp.Rational(n - 1, 2), a, l, m).doit() ** 2
              for n in range(1, 5) for l in range(n) for m in range(-l, l + 1)
              for a in [sp.Rational(2 * k - (n - 1), 2) for k in range(n)]
              if abs(m - a) <= sp.Rational(n - 1, 2))
    print(f"  CG^2(j,a; j,m-a | l,m) == CG^2(j,m-a; j,a | l,m) for all n <= 4: {sym}")


if __name__ == "__main__":
    main()
