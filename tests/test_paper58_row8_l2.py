"""Paper 58 sec:census, record R8 (docs/predictions/2026-11-14_row8_two_centre_l2.md):
two-centre integrals over hydrogenic p and d orbitals (l <= 2, n <= 3) as functions of R.

Written 2026-11-14 as a SEPARATE pass after the scoring run (CLAUDE.md sec. 9).  Closed forms
are rebuilt live from geovac (`two_center_one_electron`, `two_center_eri`); the R-operators are
frozen literals from the scored panel; the operator machinery (annihilation by sympy's own
derivatives, rank minimality, Frobenius count at 0, apparentness over Q[x]/(factor), Newton
polygon at infinity) is the one already guarded in tests/test_paper58_r_operators.py and is
loaded from there, not re-implemented.

What each test rejects:

* ``test_s_type_matches_qfd_core``: an l>0 engine that is wrong in a way invisible at l = 0
  is not caught here, but one whose conventions (normalisation, centre placement, kinetic
  sign) differ from the independent s-type engine is.
* ``test_matches_independent_quadrature``: any wrong closed form at l > 0, m != 0 -- the
  referee is cylindrical quadrature built in this file from scipy's genlaguerre and lpmv, so
  it shares no code with the prolate-spheroidal route (normalisations, Legendre signs, the
  1/r cancellation).
* ``test_kinetic_routes_agree``: a radial-Laplacian bug (the two routes share only S and V).
* ``test_operator_*``: a stored operator that does not kill the live closed form; a claimed
  order with a lower-order operator; an order that is not (rows 2/3a) the number of distinct
  rates, or (row 3b) that number plus one exactly when E1 is present.
* ``test_singularities``: R = 0 not regular singular, a log solution at 0, a non-apparent
  finite point, or characteristic roots at infinity other than the atom rates.
* ``test_function_content``: the s-type-pair rule (an s-type same-centre pair (ab) gives an
  elementary hybrid even when a_d > a_c, H9) against the l > 0 pairs that carry Ei (H5) and
  the gamma_E + ln R pair at equal rates (H2); and no ln R in f at R = 0.

NOT guarded: the full 48-case panel (only the cases below), exchange at l > 0, (AA|BC) at
l > 0, the attainment questions of Paper 63.
"""
from __future__ import annotations

import importlib.util
from fractions import Fraction
from functools import lru_cache
from math import factorial, sqrt
from pathlib import Path

import pytest
import sympy as sp

from geovac import qfd_core as Q
from geovac import two_center_one_electron as E
from geovac.two_center_eri import aabb_closed_form, hybrid_closed_form

_spec = importlib.util.spec_from_file_location(
    "_r_ops", Path(__file__).with_name("test_paper58_r_operators.py"))
RO = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(RO)

R = E.R
F = Fraction


def o(s: str) -> tuple:
    """'3d-1' -> (3, 2, -1)."""
    return (int(s[0]), "spd".index(s[1]), int(s[2:]) if len(s) > 2 else 0)


def A(s: str, Z) -> tuple:
    return ("A", F(Z)) + o(s)


def B(s: str, Z) -> tuple:
    return ("B", F(Z)) + o(s)


# ---------------------------------------------------------------- (1) s-type limit


@pytest.mark.parametrize("za,zb", [(1, 1), (2, 1), (F(3, 2), 3)])
def test_s_type_matches_qfd_core(za, zb):
    a, b = A("1s", za), B("1s", zb)
    qa, qb = ("A", F(za), 1), ("B", F(zb), 1)
    assert sp.expand(E.S(a, b) - Q.overlap(qa, qb, R)) == 0
    assert sp.expand(E.T_lap(a, b) - Q.kinetic(qa, qb, R)) == 0
    assert sp.expand(E.VA(a, b) - Q._inv_r(qa, qb, "A", R)) == 0
    assert sp.expand(E.VB(a, b) - Q._inv_r(qa, qb, "B", R)) == 0
    assert sp.expand(E.VB(a, A("1s", za)) - Q._inv_r(qa, ("A", F(za), 1), "B", R)) == 0


# ---------------------------------------------------------------- (2) independent referee


def _psi(orb, Rv, x_rho, x_z):
    import numpy as np
    from scipy.special import genlaguerre, lpmv

    cen, Z, n, l, m = orb
    Z = float(Z)
    zc = x_z if cen == "A" else x_z - Rv
    r = np.hypot(x_rho, zc)
    ct = zc / r
    rho_ = 2 * Z * r / n
    Nn = sqrt((2 * Z / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    Rnl = Nn * np.exp(-rho_ / 2) * rho_ ** l * genlaguerre(n - l - 1, 2 * l + 1)(rho_)
    Ylm = sqrt((2 * l + 1) / (4 * np.pi) * factorial(l - m) / factorial(l + m)) * lpmv(m, l, ct)
    return Rnl * Ylm


def _quad(oa, ob, Rv, which):
    import numpy as np
    from scipy import integrate

    def f(x_z, x_rho):
        k = {"S": 1.0, "VA": 1.0 / np.hypot(x_rho, x_z),
             "VB": 1.0 / np.hypot(x_rho, x_z - Rv)}[which]
        return 2 * np.pi * x_rho * _psi(oa, Rv, x_rho, x_z) * _psi(ob, Rv, x_rho, x_z) * k

    val, _ = integrate.dblquad(f, 0, 40, lambda _: -40, lambda _: 40 + Rv,
                               epsabs=1e-13, epsrel=1e-11)
    return val


SLOW = pytest.mark.slow   # dblquad over 15 s each
QUAD = [  # (a, b, kind); m >= 0 so lpmv's Condon-Shortley phase applies to both factors
    pytest.param(A("2p0", 3), B("1s", 1), "S", marks=SLOW, id="2p0-1sB-S"),
    pytest.param(A("3p1", 3), B("3d1", 2), "S", id="3p1-3d1B-S"),
    pytest.param(A("3p1", 3), B("3d1", 2), "VA", marks=SLOW, id="3p1-3d1B-VA"),
    pytest.param(A("3p1", 3), B("3d1", 2), "VB", id="3p1-3d1B-VB"),
    pytest.param(A("3d2", 3), B("3d2", F(3, 2)), "VB", id="3d2-3d2B-VB"),
    pytest.param(A("3d1", 3), A("2p1", 3), "VB", marks=SLOW, id="3d1-2p1A-VB"),
    pytest.param(A("3d0", 3), A("3s", 3), "VB", marks=SLOW, id="3d0-3sA-VB"),
]


@pytest.mark.parametrize("oa,ob,kind", QUAD)
def test_matches_independent_quadrature(oa, ob, kind):
    Rv = sp.Rational(7, 5)
    exact = float(sp.re(sp.N({"S": E.S, "VA": E.VA, "VB": E.VB}[kind](oa, ob).subs(R, Rv), 30)))
    assert abs(exact) > 1e-3                      # a vanishing pair would prove nothing
    assert abs(exact - _quad(oa, ob, float(Rv), kind)) < 1e-10


PAIRS = [(A("2p0", 3), B("1s", 1)), (A("2p1", 2), B("2p1", 2)), (A("3d0", 3), B("2p0", 2)),
         (A("3d2", 3), B("3d2", F(3, 2))), (A("3p1", 3), B("3d1", 2)), (A("2s", 2), B("3d0", 3))]


@pytest.mark.parametrize("oa,ob", PAIRS)
def test_kinetic_routes_agree(oa, ob):
    t = E.T_lap(oa, ob)
    assert t != 0
    assert sp.simplify(t - E.T_eig(oa, ob)) == 0


# ---------------------------------------------------------------- (3) frozen R-operators

OP_E5VB = [
    "2*(1625*R**5 + 7075*R**4 + 180*R**3 - 105192*R**2 - 311040*R - 186624)",
    "R*(8125*R**4 + 38625*R**3 + 86400*R**2 + 103680*R + 62208)",
    "3*R**2*(1625*R**3 + 8700*R**2 + 17280*R + 10368)",
]
OP_E4T = [
    "R*(81*R**6 - 162*R**5 - 9450*R**4 - 68544*R**3 - 228096*R**2 - 364800*R - 243200)",
    "3*(81*R**7 - 162*R**6 - 7128*R**5 - 51984*R**4 - 209664*R**3 - 523776*R**2 - 768000*R - 512000)",
    "2*R*(81*R**6 - 162*R**5 - 7020*R**4 - 42624*R**3 - 124416*R**2 - 192000*R - 128000)",
]
OP_F5 = [
    "4*R*(R**4 + 8*R**3 + 24*R**2 + 30*R + 15)",
    "2*(2*R**6 + 18*R**5 + 81*R**4 + 228*R**3 + 414*R**2 + 450*R + 225)",
    "R*(R**2 + 4*R + 5)*(2*R**3 + 9*R**2 + 18*R + 15)",
]
OP_C4 = [
    "2025*R*(736305040625*R**12 + 736305040625*R**11 + 14481910894000*R**10 - 5859181790000*R**9 - 46953912582720*R**8 - 75370560511488*R**7 - 88178222045184*R**6 - 86806907117568*R**5 - 64091836071936*R**4 - 33698355609600*R**3 - 12807084441600*R**2 - 3643696742400*R - 662490316800)",
    "497005902421875*R**14 + 1954889882859375*R**13 + 9130286637862500*R**12 + 28406062461287500*R**11 + 10707845696493000*R**10 - 45686067858063600*R**9 - 109266967860855360*R**8 - 155670299042696064*R**7 - 172382962045215744*R**6 - 158669798367320064*R**5 - 130278063245230080*R**4 - 95402505393930240*R**3 - 52574347371479040*R**2 - 19215511254466560*R - 3493729318993920",
    "R*(364470995109375*R**13 + 1124337797034375*R**12 + 6627950731703125*R**11 + 13276975114865000*R**10 + 12767175873135200*R**9 - 1577148496344960*R**8 - 27834426015648384*R**7 - 54091817332792320*R**6 - 65981564444657664*R**5 - 59212810302898176*R**4 - 42188963742351360*R**3 - 23817087776194560*R**2 - 8807109324963840*R - 1601292604538880)",
    "2*R**2*(33133726828125*R**12 + 108236840971875*R**11 + 602540975609375*R**10 + 1097444833058750*R**9 + 761581503718200*R**8 - 558785228968560*R**7 - 2124013107627264*R**6 - 2986516932410880*R**5 - 2740305568702464*R**4 - 1895321121914880*R**3 - 1078529617428480*R**2 - 400323151134720*R - 72786027479040)",
]
OP_H2 = [
    "-R*(98304*R**13 + 2322432*R**12 + 3089408*R**11 - 121239552*R**10 - 1762550208*R**9 - 12392834736*R**8 - 56243894784*R**7 - 171565839216*R**6 - 364950313380*R**5 - 557753553225*R**4 - 621502210152*R**3 - 493069581018*R**2 - 253270521360*R - 63317630340)",
    "-4*(32768*R**14 + 798720*R**13 + 1272832*R**12 - 30509312*R**11 - 438172992*R**10 - 3342346176*R**9 - 17900499792*R**8 - 66521709648*R**7 - 178591383252*R**6 - 358769115864*R**5 - 546512850630*R**4 - 620703676188*R**3 - 501197390343*R**2 - 259516103472*R - 64879025868)",
    "2*R*(32768*R**13 + 757760*R**12 - 1451008*R**11 - 111757824*R**10 - 1111132992*R**9 - 5912333232*R**8 - 19424173824*R**7 - 43628274576*R**6 - 70486389540*R**5 - 82324114857*R**4 - 66447760968*R**3 - 35835355314*R**2 - 14417561304*R - 3604390326)",
    "4*R**2*(32768*R**12 + 798720*R**11 - 300032*R**10 - 77891840*R**9 - 755301696*R**8 - 3890218944*R**7 - 13321825104*R**6 - 32718801744*R**5 - 58806872532*R**4 - 76064134104*R**3 - 67074571206*R**2 - 36043903260*R - 9010975815)",
    "R**3*(32768*R**11 + 806912*R**10 - 300032*R**9 - 77816832*R**8 - 716393280*R**7 - 3352923984*R**6 - 9968901120*R**5 - 20257675344*R**4 - 28420359516*R**3 - 26328504951*R**2 - 14417561304*R - 3604390326)",
]
OP_H9 = [
    "4*R*(17419049611490912000000*R**19 + 220267932257836170400000*R**18 + 1812715566206497460640000*R**17 + 6535866442604325218865000*R**16 - 21368869984233509426095000*R**15 - 465198734650211653894522500*R**14 - 3299379514745754635567286750*R**13 - 13871862990614710624745750625*R**12 - 34487324562236869802824146765*R**11 - 31040867569865773303893750714*R**10 + 120665768194330276018069532862*R**9 + 743945293995275605489071328818*R**8 + 2350314111622672725592182030288*R**7 + 5459071275960571755789899479488*R**6 + 10115694634232259299331594457344*R**5 + 14864915516574750161826472505856*R**4 + 16469961535823140821637536171264*R**3 + 12794971591829482413129368020224*R**2 + 6226761625898370651786524152320*R + 1436944990591931688873813265920)",
    "2*(100159535266072744000000*R**20 + 1434198962993158007800000*R**19 + 13306000818759676739855000*R**18 + 71402629728486118896542500*R**17 + 178368580910448334440375000*R**16 - 648692979365403655261012500*R**15 - 8946619184179592307764732250*R**14 - 46363315965928129619456593500*R**13 - 144602826755550875689791889305*R**12 - 270714443226898612892440224438*R**11 - 137454115967856650350414797525*R**10 + 1118159696870217960901508951721*R**9 + 4953118081342980116575617032880*R**8 + 12997109815258063346633638585728*R**7 + 26110180239991547024989034381232*R**6 + 42454440735370398261087083492688*R**5 + 55423809573969396174390960406272*R**4 + 55782361485089064320398133907456*R**3 + 40634203752948510589471489916928*R**2 + 19219396917145390546256891744256*R + 4435245442418167049136205787136)",
    "3*R*(56611911237345464000000*R**19 + 907480325564367585800000*R**18 + 8928132049756616735005000*R**17 + 53012192972253555507005000*R**16 + 182607532818055710446152500*R**15 + 67591437392581019092848750*R**14 - 3075169444704973416614368500*R**13 - 19180357648030380614654246625*R**12 - 67250865593070875354386522905*R**11 - 156865397749628610565591281288*R**10 - 208276324711158711803026075866*R**9 + 77461707366040855829280246324*R**8 + 1280098997028680309479854746928*R**7 + 4044205227585352872792958085088*R**6 + 8434516244962466544279806918304*R**5 + 13051837387150874808269744633472*R**4 + 14950074017682246395982758556672*R**3 + 11932652140114060902750902913024*R**2 + 5872593502461091555800716921856*R + 1355213885183328820569396212736)",
    "9*R**2*(4354762402872728000000*R**18 + 76840795078822682600000*R**17 + 793174335475032693985000*R**16 + 4993062154183003532060000*R**15 + 18655713743711673294682500*R**14 + 18114835468152775288076250*R**13 - 228190802453542828068147000*R**12 - 1528071542722385315296822125*R**11 - 5173143507159298104183578685*R**10 - 10872766709857747404080041941*R**9 - 11003055727077863644503524940*R**8 + 13576092993057048356139382368*R**7 + 85937375470153517631244781424*R**6 + 211934199656388469871098951056*R**5 + 355361588549805698969232753792*R**4 + 429941848277304240762758215680*R**3 + 356266896848226509052043183104*R**2 + 177957378862457319872748997632*R + 41067087429797843047557460992)",
]


@lru_cache(maxsize=None)
def live(cid: str) -> sp.Expr:
    return {
        "E5VB": lambda: E.VB(A("3p1", 3), B("3d1", 2)),
        "E4T": lambda: E.T_lap(A("3d2", 3), B("3d2", F(3, 2))),
        "F5": lambda: E.VB(A("3d2", 3), A("3d2", 3)),
        "C4": lambda: aabb_closed_form(F(3), o("2p0"), o("3d0"), F(2), o("1s"), o("2p0")),
        "H2": lambda: hybrid_closed_form(F(2), o("2p1"), o("2p1"), o("2p1"), F(2), o("2p1")),
        "H5": lambda: hybrid_closed_form(F(3), o("3d1"), o("3d1"), o("3d1"), F(3), o("2p1")),
        "H9": lambda: hybrid_closed_form(F(2), o("2s"), o("2s"), o("3d0"), F(3), o("3d0")),
    }[cid]()


# id: (operator, order, R = 0 exponents, extra order from E1)
FROZEN = {
    "E5VB": (OP_E5VB, 2, {3, -4}, 0),
    "E4T": (OP_E4T, 2, {-5, 0}, 0),
    "F5": (OP_F5, 2, {-5, 0}, 0),
    "C4": (OP_C4, 3, {-3, -5, 0}, 0),
    "H2": (OP_H2, 4, {4, -3, -5, 0}, 1),
    "H9": (OP_H9, 3, {-3, -5, 0}, 0),
}


@pytest.mark.parametrize("cid", list(FROZEN))
def test_operator_annihilates_and_is_minimal(cid):
    ops, r, _, extra = FROZEN[cid]
    P = RO._ops(ops)
    f = live(cid)
    assert len(P) - 1 == r
    assert RO.annihilates(P, f)
    assert RO.no_lower_order(f, r)
    assert r == len(RO.atom_rates(f)) + extra


@pytest.mark.slow
def test_operator_H5_order_at_least_five():
    """H5's operator is too large to freeze; its lower bound (no order-4 operator) is guarded."""
    f = live("H5")
    assert len(RO.atom_rates(f)) == 4
    assert RO.no_lower_order(f, 5)


@pytest.mark.parametrize("cid", list(FROZEN))
def test_singularities(cid):
    ops, r, exps, _ = FROZEN[cid]
    P = RO._ops(ops)
    ex = RO.exponents_at_zero(P)
    assert set(ex) == exps and all(m == 1 for m in ex.values())
    # Local solution space at R = 0: log-free, except that an operator whose f carries E1
    # has exactly one log solution (the E1 atom's own; measured in this guard pass, and the
    # same for the 1s-era hybrid operators HYB33/HYB34).  f itself is log-free:
    # test_function_content.
    has_e1 = "E1" in _kinds(live(cid))[0]
    dim, _ = RO.laurent_dim_at_zero(P)
    assert dim == r - (1 if has_e1 else 0)
    rank, rates, deg = RO.infinity_rank_and_rates(P)
    assert rank == 1 and deg == r
    assert rates == RO.atom_rates(live(cid))
    assert RO.finite_points_apparent(P)


# ---------------------------------------------------------------- (4) function content


def _kinds(f):
    atoms, reps = RO.decompose(f)
    return {k for (_, k, _) in atoms}, set(reps)


def test_function_content():
    k9, c9 = _kinds(live("H9"))          # s-type pair, a_d > a_c
    assert k9 == {"1"} and sp.EulerGamma not in c9
    k5, c5 = _kinds(live("H5"))          # l > 0 pair, a_d > a_c
    assert "Ei" in k5 and "log" not in k5 and sp.EulerGamma not in c5
    k2, c2 = _kinds(live("H2"))          # l > 0 pair, a_d = a_c
    assert {"E1", "log"} <= k2 and "Ei" not in k2 and sp.EulerGamma in c2
    for cid in ("H2", "H5"):
        assert RO.lnR_coefficient(live(cid)) == 0
