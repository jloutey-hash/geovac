"""Guards for Paper 59 sec:density_ring, the kite at a blind second point and its ninth singular point
(records docs/predictions/2026-10-26_parametric_kite_operator.md, 2026-10-26_kite_ninth_point.md,
2026-10-27_kite_ninth_point_origin.md).  Data: tests/data/paper59_kite_p2_2026_10_27.json, built by
debug/build_kite_guard_data_2026_10_27.py (a copy of stored results; no debug/ import here).

Claims guarded, as the paper states them:
  * at (a,b,g1,g2) = (1,3/2,7/3,8/3) the kite's operator has eighteen masters and order eleven (degree 120);
  * the eight cut points, minus the squared sums and differences of the rates on each cut, are singular with the
    predicted exponents (-1/2 at the crossed three-line cuts; doubled exponents at the two-line cuts); the
    pairings of adjacent lines give none;
  * a ninth genuine point at z9 = -111/1120 with exponents -7/2, 0, 1, ..., 9;
  * eq:kite_ninth_point, z9 = (A G2 - B G1)(A - B - G1 + G2) / ((A - G1)(B - G2)) in squared rates, reproduces the
    measured ninth point (53 rate sets mod P; the family (1,2,3/2,g2) at two primes), is zero at the original rates
    (a g2 = b g1), has the two stated symmetries and poles at adjacent equal rates;
  * mechanism: with all five lines on shell and the loop momenta in one plane, elimination gives exactly that factor
    and no other; the multipliers are nonzero and U vanishes on them.
"""
import json
from pathlib import Path

import sympy as sp
from flint import fmpq, fmpq_poly

DATA = json.loads((Path(__file__).parent / "data" / "paper59_kite_p2_2026_10_27.json").read_text())
R = sp.Rational
A_, B_, G1_, G2_ = sp.symbols("A B G1 G2")
P2 = {A_: R(1), B_: R(9, 4), G1_: R(49, 9), G2_: R(64, 9)}
Z9 = (A_ * G2_ - B_ * G1_) * (A_ - B_ - G1_ + G2_) / ((A_ - G1_) * (B_ - G2_))


def _op():
    return [fmpq_poly([int(c) for c in p]) for p in DATA["p2_op"]["polys"]]


def local_exponents(P, v):
    """(multiplicity of v in the leading coefficient, roots of the indicial polynomial at z = v)."""
    sh = [p(fmpq_poly([fmpq(int(v.p), int(v.q)), 1])) for p in P]
    n = len(P) - 1
    vals = [next((i for i, c in enumerate(p.coeffs()) if c != 0), None) for p in sh]
    delta = min(j - k for k, j in enumerate(vals) if j is not None)
    r = sp.Symbol("r")
    ind = 0
    for k, p in enumerate(sh):
        cs = p.coeffs()
        j = k + delta
        if 0 <= j < len(cs) and cs[j] != 0:
            ind += R(int(cs[j].p), int(cs[j].q)) * sp.ff(r, k)
    mult = next(i for i, c in enumerate(sh[n].coeffs()) if c != 0)
    return mult, sp.roots(sp.Poly(sp.expand(ind), r))


def test_blind_operator_has_order_eleven_and_eighteen_masters():
    P = _op()
    assert len(DATA["p2_op"]["masters"]) == 18
    assert len(P) - 1 == 11
    assert P[-1].degree() == 120


def test_cut_points_and_adjacent_pairings():
    P = _op()
    a, b, g1, g2 = R(1), R(3, 2), R(7, 3), R(8, 3)
    three_line = [-(a + g2) ** 2, -(a - g2) ** 2, -(b + g1) ** 2, -(b - g1) ** 2]
    two_line = [-(a + b) ** 2, -(a - b) ** 2, -(g1 + g2) ** 2, -(g1 - g2) ** 2]
    adjacent = [-(a + g1) ** 2, -(a - g1) ** 2, -(b + g2) ** 2, -(b - g2) ** 2]
    for v in three_line:
        mult, ex = local_exponents(P, v)
        assert mult == 1 and ex == {R(-1, 2): 1, **{k: 1 for k in range(10)}}, (v, ex)
    for v in two_line:
        mult, ex = local_exponents(P, v)
        assert mult == 3 and ex.get(0) == 2 and sum(ex.values()) == 11, (v, ex)
        assert all(sp.Rational(e).q == 1 for e in ex), (v, ex)
    for v in adjacent:
        assert local_exponents(P, v)[0] == 0, v


def test_ninth_point_is_genuine_with_exponent_minus_seven_halves():
    assert Z9.subs(P2) == R(-111, 1120)
    mult, ex = local_exponents(_op(), R(-111, 1120))
    assert mult == 1 and ex == {R(-7, 2): 1, **{k: 1 for k in range(10)}}, ex


def test_closed_form_reproduces_the_measured_ninth_point():
    num, den = sp.fraction(sp.together(Z9.subs(A_, 1)))
    fnum, fden = sp.lambdify((B_, G1_, G2_), num), sp.lambdify((B_, G1_, G2_), den)
    C = DATA["stageC"]
    P = C["P"]
    for i, (b2, g12, g22) in enumerate(C["pts"]):
        n_, d_ = int(fnum(b2, g12, g22)) % P, int(fden(b2, g12, g22)) % P
        assert n_ == C["vals"][str(i)] * d_ % P, i
    # the family (1, 2, 3/2, g2): s = g2^2, z9 = (s - 9)(4 s - 21) / (5 (s - 4))
    s = sp.Symbol("s")
    fam = sp.cancel(Z9.subs({A_: 1, B_: 4, G1_: R(9, 4), G2_: s}))
    assert sp.simplify(fam - (s - 9) * (4 * s - 21) / (5 * (s - 4))) == 0
    for st in DATA["stageA"]:
        P = st["P"]
        assert len(st["z9"]) >= 7
        for sv, zv in zip(st["s"], st["z9"]):
            assert (sv - 9) * (4 * sv - 21) % P == 5 * (sv - 4) * zv % P


def test_closed_form_structure():
    base = {A_: 1, B_: 4, G1_: R(9, 4), G2_: 9}                    # (1, 2; 3/2, 3): a g2 = b g1
    assert Z9.subs(base) == 0
    assert sp.simplify(Z9 - Z9.subs({A_: B_, B_: A_, G1_: G2_, G2_: G1_}, simultaneous=True)) == 0
    assert sp.simplify(Z9 - Z9.subs({A_: G1_, G1_: A_, B_: G2_, G2_: B_}, simultaneous=True)) == 0
    den = sp.fraction(sp.cancel(Z9))[1]
    assert sp.factor(den).has(A_ - G1_) and sp.factor(den).has(B_ - G2_)


def test_on_shell_elimination_gives_the_model_factor_only():
    k1, k2, q1, q2, om, sat, z = sp.symbols("k1 k2 q1 q2 om sat z")
    pts = [(R(1), R(9, 4), R(49, 9), R(64, 9)),
           (R(100), R(81, 16), R(1521, 16), R(729, 4)),
           (R(225), R(1, 16), R(729, 49), R(1296))]
    for A, B, G1, G2 in pts:
        eqs = [k1**2 + k2**2 + A, (om - k1)**2 + k2**2 + B, (k1 - q1)**2 + (k2 - q2)**2 + G1,
               (om - k1 + q1)**2 + (q2 - k2)**2 + G2, q1**2 + q2**2, sat * k2 - 1]
        gb = sp.groebner(eqs, sat, k1, k2, q1, q2, om, order="lex")
        elim = [g for g in gb.exprs if g.free_symbols <= {om}]
        assert len(elim) == 1
        ez = sp.Poly(elim[0], om)
        assert all(d % 2 == 0 for (d,) in ez.monoms())
        pz = sp.Poly(sum(c * z ** (d // 2) for (d,), c in ez.terms()), z)
        z9 = Z9.subs({A_: A, B_: B, G1_: G1, G2_: G2})
        assert pz.degree() == 1 and pz.eval(z9) == 0, (A, B, G1, G2, pz)


def test_landau_multipliers_nonzero_and_U_vanishes():
    A, B, G1, G2 = R(1), R(9, 4), R(49, 9), R(64, 9)
    w0 = sp.sqrt(R(-111, 1120))
    K1 = sp.simplify((w0**2 - A + B) / (2 * w0))
    kv = sp.Matrix([K1, sp.sqrt(-A - K1**2)])
    pv = sp.Matrix([w0, 0]) - kv
    Q1, Q2 = sp.symbols("Q1 Q2")
    sol = sp.solve([kv.dot(sp.Matrix([Q1, Q2])) - (G1 - A) / 2, pv.dot(sp.Matrix([Q1, Q2])) - (B - G2) / 2],
                   [Q1, Q2], dict=True)[0]
    qv = sp.Matrix([sp.simplify(sol[Q1]), sp.simplify(sol[Q2])])
    for v in (kv.dot(kv) + A, pv.dot(pv) + B, qv.dot(qv), (kv - qv).dot(kv - qv) + G1, (pv + qv).dot(pv + qv) + G2):
        assert sp.simplify(sp.expand(v)) == 0
    al = sp.symbols("al1:6")
    L1 = al[0] * kv - al[1] * pv + al[2] * (kv - qv) - al[3] * (pv + qv)
    L2 = -al[2] * (kv - qv) + al[3] * (pv + qv) + al[4] * qv
    M = sp.Matrix([[sp.diff(e, a) for a in al] for e in list(L1) + list(L2)]).applyfunc(
        lambda e: sp.radsimp(sp.simplify(e)))
    ns = M.nullspace(simplify=True)
    assert len(ns) == 1
    v = (ns[0] / ns[0][4]).applyfunc(lambda e: sp.radsimp(sp.simplify(e)))
    assert list(v) == [R(-49000, 9711), R(44800, 9711), R(49000, 5511), R(-44800, 5511), 1]
    a1, a2, a3, a4, a5 = v
    assert sp.simplify((a1 + a2) * (a3 + a4) + a5 * (a1 + a2 + a3 + a4)) == 0
