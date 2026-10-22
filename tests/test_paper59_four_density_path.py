"""Paper 59 sec:four_density, the real three-dimensional path of four densities.

Record docs/predictions/2026-10-22_four_density_path.md (FP1-FP3, committed 66af9822 before computing).
Written as a separate step from the result; no debug/ import.  The exact pieces found by search are
loaded from tests/data/paper59_ladder_2026_10_22.json and re-verified here.  Shared operator helpers
come from tests/test_paper59_generic_slice.py and tests/test_paper59_four_density_star.py.

On the slice (1,t,1,1) the path is P(t) = (16/pi^4) t^3 int_0^oo K0(t x) chi(x) H(x) dx with
H(x) = int_0^oo min(x,s) K0(s) chi(s) ds.  W12 = Y4 o R8 annihilates P, with R8 = X4 o R3 o (theta - 2),
R3 = L1 o R2 and X4 = Z2 o A14^-2 o Z2' (the stored Q-forms clear denominators on the left, so the
stored Z2' carries A14^2).  Y4 is equivalent to F_4's operator (eq:four_density_op); R2, Z2' and Z2 are
equivalent to Legendre's.  Not re-run here: the order-74 module operator and its exact division by W12
(python-flint, in the record), and the momentum-space check of the normalisation.
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

import test_paper59_four_density_star as ST
import test_paper59_generic_slice as GS

t, u = GS.t, GS.u
DATA = json.loads((Path(__file__).parent / "data" / "paper59_ladder_2026_10_22.json").read_text())["path"]
LEG = [-t**2, -2 * t**2, 1 - t**2]                  # theta^2 - t^2 (theta + 1)^2, Legendre's operator
TH2 = [sp.Integer(-2), sp.Integer(1)]
A14 = 44*t**14 + 2865*t**12 + 12733*t**10 - 2360*t**8 - 26328*t**6 + 29360*t**4 - 23424*t**2 + 5760


def _Q(key):
    return GS._coeff_list([sp.sympify(s, locals={"u": u}) for s in DATA[key]])


def _map(name):
    keys = [k for k in DATA["maps"] if k == name] or [k for k in DATA["maps"] if k.startswith(name)]
    assert len(keys) == 1, keys
    return [sp.sympify(s, locals={"t": t}) for s in DATA["maps"][keys[0]]]


def _zero(rem):
    return all(sp.cancel(r) == 0 for r in rem)


def _prop(a, b):
    r = sp.cancel(a[-1] / b[-1])
    return len(a) == len(b) and all(sp.cancel(x - r * y) == 0 for x, y in zip(a, b))


def _is_equivalence(target, A, source):
    """target o A = 0 mod source, A nonzero of order below the irreducible source: A maps Sol(source)
    injectively into Sol(target); equal orders make it an equivalence.  A 1/1000 change must break it."""
    ok = _zero(GS._rdiv(GS._compose(target, A), source)[1])
    A2 = list(A)
    A2[0] = A2[0] * sp.Rational(1001, 1000)
    broken = not _zero(GS._rdiv(GS._compose(target, A2), source)[1])
    return ok and broken and len(A) < len(source) and any(sp.cancel(a) != 0 for a in A)


def test_path_bottom_pieces():
    """Rejects a wrong bottom of the path's factor: R2 is an exact right factor of R3 with an order-one
    quotient, and R2 is equivalent to Legendre's operator by the stored map."""
    R2, R3 = _Q("R2_Q"), _Q("R3_Q")
    q, rem = GS._rdiv(R3, R2)
    assert _zero(rem) and len(q) - 1 == 1
    assert len(R2) - 1 == 2 and _is_equivalence(R2, _map("Legendre theta^2"), LEG)


@pytest.mark.slow
def test_path_middle_pieces_are_legendre():
    """Rejects 'the order-four middle piece X4 holds a K3 or F_4 piece': X4 right-divides by Z2' with
    quotient Z2 o A14^-2 (Z2 up to a rational gauge), and both Z2' and Z2 are equivalent to Legendre's
    operator by the stored maps.  Also R8 = X4 o R3 o (theta - 2) up to a left factor."""
    X4, Z2, Z2p, R3, R8 = _Q("X4_Q"), _Q("Z2_Q"), _Q("Z2p_Q"), _Q("R3_Q"), _Q("R8_Q")
    q, rem = GS._rdiv(X4, Z2p)
    assert _zero(rem) and _prop(GS._compose(q, [A14**2]), Z2)
    assert _is_equivalence(Z2p, _map("Legendre -> Z2'"), LEG)
    assert _is_equivalence(Z2, _map("Legendre -> Z2"), LEG)
    q, rem = GS._rdiv(R8, GS._compose(R3, TH2))
    assert _zero(rem) and _prop(q, X4)


@pytest.mark.slow
def test_path_f4_piece_is_exact():
    """Rejects 'no piece of the path's factor is equivalent to F_4's' and a wrong stored piece: W12 is
    Y4 o R8 up to a left factor, and Y4 is equivalent to F_4's operator by the stored map."""
    W12, Y4, R8 = _Q("W12_Q"), _Q("Y4_Q"), _Q("R8_Q")
    q, rem = GS._rdiv(W12, R8)
    assert _zero(rem) and _prop(q, Y4)
    assert len(Y4) - 1 == 4 and _is_equivalence(Y4, _map("F_4"), GS._coeff_list(ST.Q4))


# ------------------------------------------------------------------ numerics on P
def _grids(dps=32, hinv=24):
    """DE nodes for the path (weight chi H) and the star (weight chi^3 / x); H by quadrature between
    consecutive nodes, H(x) = int_0^x s h + x int_x^oo h, h = K0 chi."""
    mp.mp.dps = dps
    hh = mp.mpf(1) / hinv
    xs, ws = [], []
    for k in range(int(-5 * hinv), int(5.2 * hinv) + 1):
        sv = k * hh
        xs.append(mp.exp(sv - mp.exp(-sv)))
        ws.append(hh * xs[-1] * (1 + mp.exp(-sv)))
    f = lambda s: mp.besselk(0, s) * ST._chi(s)
    lo, a, acc = [], mp.mpf(0), mp.mpf(0)
    for x in xs:
        acc += mp.quad(lambda s: s * f(s), [a, x])
        lo.append(acc)
        a = x
    hi = [mp.mpf(0)] * len(xs)
    acc = mp.quad(f, [xs[-1], mp.inf])
    hi[-1] = acc
    for k in range(len(xs) - 2, -1, -1):
        acc += mp.quad(f, [xs[k], xs[k + 1]])
        hi[k] = acc
    chis = [ST._chi(x) for x in xs]
    path = [(x, w, c * (p + x * q)) for x, w, c, p, q in zip(xs, ws, chis, lo, hi)]
    star = [(x, w, c**3 / x) for x, w, c in zip(xs, ws, chis)]
    return path, star


@pytest.mark.slow
def test_w12_annihilates_path_and_needs_f4():
    """Rejects 'the path sits below the F_4 piece' (R8 P = 0), a wrong order-four piece (a 1/200 change in
    Y4's theta^2 coefficient must spoil the annihilation) and a wrong integrand (W12 must not kill the
    star): W12 kills P to 1e-20 of its largest term at 32 digits, at t = 1.3 and 2.5."""
    Y4, R8 = _Q("Y4_Q"), _Q("R8_Q")
    W12 = GS._compose(Y4, R8)
    Yp = list(Y4)
    Yp[2] = Y4[2] * sp.Rational(201, 200)
    W12p = GS._compose(Yp, R8)
    path, star = _grids()
    for tv in ("1.3", "2.5"):
        assert ST._apply(path, W12, tv) < mp.mpf(10) ** -20
        assert ST._apply(path, W12p, tv) > mp.mpf(10) ** -4
        assert ST._apply(path, R8, tv) > mp.mpf(10) ** -3
        assert ST._apply(star, W12, tv) > mp.mpf(10) ** -6
