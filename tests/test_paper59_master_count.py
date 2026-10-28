"""Guard for Paper 59 sec:density_ring, the master count of the one-mass slice and the two-offset fibre (record
docs/predictions/2026-10-27_master_count_tabulation.md, MC0-MC2b, MC6).

Masters are counted as critical points of the twisted one-form (Euler characteristic), excluding removed points.
On the slice, L(D, rho) = int_1^oo e^{-D x} x^{1-n1} (x^2-1)^{-1/2} (sigma x^2 + 1 - sigma)^{-n2/2} dx with the
exponential twist (Route X), or its Mellin dual (Route M).  Claims guarded, as the paper states them:
  * the slice needs four masters (n1 = 1), by both routes, with generic and with the actual exponents;
  * the production index range needs a fifth: n1 >= 2 gives five (the pole at x = 0), Route X;
  * the physical fibre, with two offsets, needs six (conic quotient, Route F);
  * controls of the counting method: Bessel 2, Gauss 2F1 2, the one-loop bubble 3 (Lee-Pomeransky, mod p).
Code: a copy of the relevant parts of debug/mct_counts.py (no debug/ import).
"""
import itertools

import sympy as sp

NU = [sp.Rational(3, 7), sp.Rational(5, 11), sp.Rational(2, 13), sp.Rational(7, 17), sp.Rational(4, 19)]
x, tt = sp.symbols("x t")
PTS = [(sp.Rational(13, 7), sp.Rational(5, 9)), (sp.Rational(41, 3), sp.Rational(17, 29)),
       (sp.Rational(2, 37), sp.Rational(88, 13))]


def count_1d(omega, var, removed):
    num, _ = sp.fraction(sp.cancel(sp.together(omega)))
    N = sp.Poly(sp.Poly(num, var).clear_denoms()[1].as_expr(), var)
    Nsq = sp.Poly(sp.quo(N, sp.gcd(N, N.diff(var))), var)
    for g in removed:
        G = sp.Poly(sp.Poly(g, var).clear_denoms()[1].as_expr(), var)
        if G.degree() <= 0:
            continue
        while True:
            h = sp.gcd(Nsq, G)
            if h.degree() <= 0:
                break
            Nsq = sp.Poly(sp.quo(Nsq, h), var)
    return Nsq.degree()


def _quad(sig):
    return sig * x**2 + 1 - sig


def slice_count(dhat, sig, n1, n2, generic, route):
    if route == "X":
        ex = 1 - n1
        xe = NU[3] if (generic and ex != 0) else ex
    else:
        xe = NU[4] + 1 - n1
    if generic:
        facs = [(x - 1, NU[0]), (x + 1, NU[1]), (_quad(sig), NU[2]), (x, xe)]
    else:
        facs = [(x - 1, sp.Rational(-1, 2)), (x + 1, sp.Rational(-1, 2)), (_quad(sig), sp.Rational(-n2, 2)), (x, xe)]
    facs = [(g, e) for g, e in facs if e != 0]
    om = (-dhat if route == "X" else 0) + sum(e * sp.diff(g, x) / g for g, e in facs)
    return count_1d(om, x, [g for g, _ in facs])


def test_slice_needs_four_masters():
    for dh, sg in PTS:
        for route, generic in itertools.product("XM", (True, False)):
            assert slice_count(dh, sg, 1, 1, generic, route) == 4, (dh, sg, route, generic)
        for n2 in (2, 3, 5):
            assert slice_count(dh, sg, 1, n2, False, "X") == 4, (dh, sg, n2)


def test_production_index_range_needs_a_fifth():
    for dh, sg in PTS:
        for n1, n2 in itertools.product((2, 3, 5), (1, 3, 5)):
            assert slice_count(dh, sg, n1, n2, False, "X") == 5, (dh, sg, n1, n2)
            assert slice_count(dh, sg, n1, n2, True, "X") == 5, (dh, sg, n1, n2)


def test_two_offset_fibre_needs_six():
    for c1, r, D1, D2 in [(sp.Rational(7, 53), sp.Rational(11, 37), sp.Rational(9, 4), sp.Rational(5, 13)),
                          (sp.Rational(30, 43), sp.Rational(3, 41), sp.Rational(1, 7), sp.Rational(22, 3))]:
        c2 = c1 * r**2
        X = (2 * c1 * tt - c1 * tt**2 - c2) / (c2 - c1 * tt**2)
        W = 1 + tt * (X - 1)
        assert sp.simplify(c2 * (X**2 - 1) - c1 * (W**2 - 1)) == 0
        branch = [tt, tt - 1, tt - c2 / c1]
        f = -D1 * X - D2 * W
        rem = branch + [c2 - c1 * tt**2]
        assert count_1d(sp.diff(f, tt) + sum(NU[j] / branch[j] for j in range(3)), tt, rem) == 6
        assert count_1d(sp.diff(f, tt) - sum(sp.Rational(1, 2) / b for b in branch), tt, rem) == 6


def test_counting_method_controls():
    D = sp.Rational(13, 7)
    assert count_1d(-D + NU[0] / (x - 1) + NU[1] / (x + 1), x, [x - 1, x + 1]) == 2          # Bessel
    zz = sp.Rational(5, 17)
    assert count_1d(NU[0] / x - NU[1] / (1 - x) - NU[2] * zz / (1 - zz * x), x, [x, 1 - x, 1 - zz * x]) == 2
    # bubble by Lee-Pomeransky critical points, mod p
    x1, x2, ts = sp.symbols("x1 x2 tsat")
    U = x1 + x2
    G = sp.expand(U + x1 * x2 * sp.Rational(37, 101) + U * (x1 + sp.Rational(9, 4) * x2))
    s = sp.Rational(9, 23)
    eqs = [sp.expand(NU[0] * G - s * x1 * sp.diff(G, x1)), sp.expand(NU[1] * G - s * x2 * sp.diff(G, x2)),
           ts * x1 * x2 * G - 1]
    gb = sp.groebner([sp.together(e).as_numer_denom()[0] for e in eqs], ts, x1, x2, order="grevlex", modulus=32003)
    lead = [sp.Poly(e, ts, x1, x2, modulus=32003).monoms(order="grevlex")[0] for e in gb.exprs]
    bound = [min(m[i] for m in lead if all(m[j] == 0 for j in range(3) if j != i)) for i in range(3)]
    cnt = sum(1 for mono in itertools.product(*[range(b) for b in bound])
              if not any(all(mono[j] >= m[j] for j in range(3)) for m in lead))
    assert cnt == 3
