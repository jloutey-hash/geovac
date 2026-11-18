"""Paper 58 -- the order of the unequal-exponent 1s exchange operator, derived from Pachucki's ODE.

Guard for the sentence proposed for Paper 58 (papers/group2_quantum_chemistry/paper_58_abelian_residue.tex):

  "At (a,b) = (1,3/2) and (1,2) ... differentiating Pachucki's second-order equation for the two-centre master
  integral once in each of the four exponents at the exchange point, where p^2 vanishes to second order, gives
  inhomogeneous equations whose unique decaying solutions lie, by exact linear algebra, in the span of nine
  functions ... The resulting closed form agrees with the independently computed integral ... The order is
  therefore at most 9, and the operator the closed form yields coincides with the guessed one modulo both primes."

Written as a SEPARATE pass from the derivation (CLAUDE.md Sec. 9).  Prior art: K. Pachucki, Phys. Rev. A 80, 032520
(2009), arXiv:0905.1635; its formulas are transcribed in tests/_pachucki_exchange.py with the arXiv equation label
at each (p^2 label 20, the ODE 28, the source 30-34, f0 and F^(-1) 59-60).  Nothing here imports debug/.

K(R) = (1s_A^a 1s_B^b | 1s_A^a 1s_B^b), normalised Slater 1s, = 16 a^3 b^3 / R * d_u3 d_u2 d_w2 d_w3 f at the
exchange point u3 = w2 = a, u2 = w3 = b.

WHAT IS PROVED HERE, AND WHAT IS NOT.
  * Exact (Fractions / python-flint, LC = ln((a+b)/(b-a)) carried as an indeterminate): the 15 derivative equations
    have a solution in the span of 11 decaying functions, with ODE residual identically zero and no free ansatz
    column; K and its derivatives lie in the Q[R, 1/R]-span of NINE functions; an operator of order 9, degree 16
    annihilates the closed form, and it is the committed guessed operator modulo each prime (and the CRT-reconstructed
    operator over Q of tests/test_paper58_exchange_unequal.py, exactly).  Hence order <= 9.
  * Numerical: the closed form is the exchange integral (geovac.qfd_core.exchange_hp, an independent route).
  * NOT proved: that the solution found is f_S rests on uniqueness of the decaying solution (the kernel of (r d)'
    is {1, ln r}); the ansatz functions all decay (checked), the decay of f_S itself is an argument, not a test.
    Order >= 9 for the true function needs linear independence of the nine functions over C(R); the rank 9 computed
    here is the rank in the formal span.

THE NINE FUNCTIONS (a correction to the sentence as drafted).  Over Q(R) the two brackets must carry Pachucki's
constant (label 60):   psi_a = e^{-2aR} [Ei(2(a-b)R) - 2 Ei((a-b)R) + ln((b-a)/(b+a))],  and its a <-> b twin.
With the bare bracket e^{-2aR}[Ei(2(a-b)R) - 2 Ei((a-b)R)] the coefficient of e^{-2aR} contains ln((a+b)/(b-a)),
so K is then in the span over Q(R, ln((a+b)/(b-a))), not over Q(R) (test_K_lies_in_the_nine_function_span).

What each test rejects is named in its docstring; fire-test plants are listed at the end of the file.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

import _pachucki_exchange as PX

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "debug" / "data"

PAIRS = {"1_2": (Fraction(1), Fraction(2)), "1_3half": (Fraction(1), Fraction(3, 2))}

# The committed guessed operators (records EX and EXb; `minimal_operator.coeffs_mod_p`, ascending powers of R).
OPS = {
    "1_2": ("exch_run_1_2_67108859_EXS.json", "exch_run_1_2_33554393_EXS.json"),
    "1_3half": ("exch_run_1_3_2_67108859_EXB.json", "exch_run_1_3_2_33554393_EXB.json"),
}
# sha256 of the committed files, CRLF folded to LF (git may check text out with either)
SHA = {
    "exch_run_1_2_67108859_EXS.json": "e84e0e89e423bf7fd708912764a447d37d450f6af198abc99185421b3a60414a",
    "exch_run_1_2_33554393_EXS.json": "8a21085017375f2957d2bf700c9286c492b9ec9929bc47c08052218377d1c96b",
    "exch_run_1_3_2_67108859_EXB.json": "931f996379642a88a33867926c9a7ae927de4f4abfe5d91bcf9911996f8352f7",
    "exch_run_1_3_2_33554393_EXB.json": "72735d32f8ce2041a4df67d6b13b876a1b42ed9b767a8ff938e122698267a8a4",
}
PRIMES = {"exch_run_1_2_67108859_EXS.json": 67108859, "exch_run_1_2_33554393_EXS.json": 33554393,
          "exch_run_1_3_2_67108859_EXB.json": 67108859, "exch_run_1_3_2_33554393_EXB.json": 33554393}
FILE_CASES = [(key, name) for key in sorted(OPS) for name in OPS[key]]

ORDER, DEGREE = 9, 16
COEFF_DEGREES = [14, 15] + [16] * 8
TIE = -2                                  # coordinate of Ei((a-b)R) / coordinate of Ei(2(a-b)R), under e^{-2aR}
R0_EXPONENTS = {-1: 2, 0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 6: 1, 9: 1}

# exchange_hp at dps 80, tau_max = 32 (records EX / EXb numeric step), typed from the records
RECORD_DPS80 = {
    ("1_3half", "1.4"): "0.272600564598402583271880658867844616019656149590167032857415",
    ("1_3half", "3.015"): "0.0283454185592348529587015928859179181830153685961323648344686",
    ("1_2", "1.4"): "0.20632030202807079004021890670600820962904556376580789324225",
    ("1_2", "3.015"): "0.0152431491774871955879056964903847142468863512279614803377193",
}


# ---------------------------------------------------------------- shared exact objects
@lru_cache(maxsize=None)
def solved(key: str):
    a, b = PAIRS[key]
    alg = PX.Algebra(a, b)
    fS, free, Fs, P2 = alg.solve_all()
    K = alg.K(fS)
    ders = alg.derivatives(K, 10)
    return alg, fS, free, Fs, P2, K, ders


@lru_cache(maxsize=None)
def derived_operator(key: str):
    ders = solved(key)[6]
    ops = PX.operators_annihilating(ders, ORDER, DEGREE)
    assert len(ops) == 1, f"{len(ops)} independent operators of order <= {ORDER}, degree <= {DEGREE}"
    op = [PX.trim(c) for c in ops[0]]
    if op[-1][-1] < 0:
        op = [[-x for x in c] for c in op]
    return op


def _bytes(name: str) -> bytes:
    return (DATA / name).read_bytes().replace(b"\r\n", b"\n")


def committed(name: str):
    J = json.loads(_bytes(name).decode())
    assert J["prime"] == PRIMES[name]
    return [[int(x) for x in c] for c in J["minimal_operator"]["coeffs_mod_p"]], J["prime"]


def nine_functions(alg):
    """The nine functions as elements of the function algebra (LC = ln((a+b)/(b-a)) an indeterminate)."""
    a, b = alg.a, alg.b
    one = Fraction(1)

    def psi(x, y):   # e^{-2xR}[Ei(2(x-y)R) + TIE Ei((x-y)R) - LC]
        return {("Ei", -2 * x, 2 * (x - y)): {0: one}, ("Ei", -2 * x, (x - y)): {0: one * TIE},
                ("E", -2 * x, 1): {0: -one}}

    return [
        ({("f",): {0: one}}, ("f",)),
        ({("Ei", Fraction(0), -(a + b)): {0: one}}, ("Ei", Fraction(0), -(a + b))),        # -E1((a+b)R)
        ({("E", -(a + b), 0): {0: one}}, ("E", -(a + b), 0)),
        ({("E", -2 * a, 0): {0: one}}, ("E", -2 * a, 0)),
        ({("E", -2 * b, 0): {0: one}}, ("E", -2 * b, 0)),
        ({("Ei", 2 * a, -2 * (a + b)): {0: one}}, ("Ei", 2 * a, -2 * (a + b))),            # -e^{2aR} E1(2(a+b)R)
        ({("Ei", 2 * b, -2 * (a + b)): {0: one}}, ("Ei", 2 * b, -2 * (a + b))),
        (psi(a, b), ("Ei", -2 * a, 2 * (a - b))),
        (psi(b, a), ("Ei", -2 * b, 2 * (b - a))),
    ]


# ---------------------------------------------------------------- (1) the exchange point
def test_p2_vanishes_to_second_order_at_the_exchange_point():
    """Symbolic in (a, b): p^2 (label 20) and its gradient vanish at u3 = w2 = a, u2 = w3 = b, and the Hessian is
    the stated one, with entry d_{w2} d_{w3} p^2 = 4ab(a^2+b^2)/(a^2-b^2)^2.

    Rejects: a p^2 that vanishes only to first order (then first derivatives of the ODE keep a p^2 r f term and the
    equations are not the ones solved), a mistranscribed p^2, a wrong Hessian (e.g. a^2 - b^2 in the numerator).
    """
    A, B = sp.symbols("A B", positive=True)
    p2 = PX.p2_expr()
    at = PX.exchange_point(A, B)
    assert sp.simplify(p2.subs(at)) == 0
    for v in PX.THETA:
        assert sp.simplify(sp.diff(p2, v).subs(at)) == 0
    den = (A**2 - B**2) ** 2
    h1 = 4 * A * B * (A**2 + B**2) / den
    h2 = 8 * A**2 * B**2 / den
    # variable order (u3, u2, w2, w3) = exponents of r1A, r1B, r2A, r2B
    expected = sp.Matrix([[-h2, h1, h2, -h1], [h1, -h2, -h1, h2], [h2, -h1, -h2, h1], [-h1, h2, h1, -h2]])
    H = sp.Matrix(4, 4, lambda i, j: sp.diff(p2, PX.THETA[i], PX.THETA[j]).subs(at))
    assert (H - expected).applyfunc(sp.simplify) == sp.zeros(4, 4)
    assert sp.simplify(sp.diff(p2, PX.w2, PX.w3).subs(at) - h1) == 0
    assert h1.subs({A: 1, B: 2}) == sp.Rational(40, 9)
    # ln|x| of labels 33-34 is entered as ln(-x): x = -(a+b)/(b-a) at the exchange point
    assert sp.simplify(PX.ln_argument().subs(at) + (A + B) / (B - A)) == 0


@pytest.mark.parametrize("key", sorted(PAIRS))
def test_label60_is_the_antiderivative_of_label30(key):
    """d/dr F^(-1) (label 60, typed directly into the algebra) = F (labels 30-34, differentiated symbolically and
    parsed) at the exchange point, exactly; and every function of the ansatz decays.

    Rejects: a sign or factor error in either transcription (e.g. +2 for -2 on one Ei of label 60, a wrong a+-/b+-),
    and an ansatz containing a non-decaying function (which would void uniqueness).
    """
    alg, _, _, Fs, _, _, _ = solved(key)
    assert Fs[(0, 0, 0, 0)], "empty source"
    assert alg.diff(alg.FM1) == Fs[(0, 0, 0, 0)]
    for k in alg.keys11() + alg.keys_lc():
        if k[0] == "E":
            assert k[1] < 0
        elif k[0] == "Ei":
            assert k[1] + k[2] < 0 and k[2] != 0     # e^{lam r} Ei(mu r) ~ e^{(lam+mu) r} / (mu r)
    assert len(alg.keys11()) == len(set(alg.keys11())) == 11


# ---------------------------------------------------------------- (2) the fifteen equations
def _residual(alg, al, fS, Fs, P2):
    """(r f_S')' + F_S - r sum_{0 != U subset S} (d_U p^2) f_{S \\ U}, written independently of the helper's rhs."""
    res = PX.fadd(alg.L0(fS[al]), Fs[al])
    for U in PX.SUBSETS:
        if sum(U) == 0 or any(u > s for u, s in zip(U, al)):
            continue
        T = tuple(s - u for s, u in zip(al, U))
        res = PX.fadd(res, PX.fscale(fS[T], -P2[U], 1))
    return res


@pytest.mark.parametrize("key", sorted(PAIRS))
def test_fifteen_derivative_equations_solve_exactly_in_the_span(key):
    """All 15 equations (*) are solved in the 11-function span (Laurent window r^-6..r^6): ODE residual identically
    zero, no free ansatz column (the span holds no kernel element, so the solution in it is unique), and the
    second-order coupling d_U p^2 is the Hessian of test (1).

    Rejects: a solution of the wrong equation (p^2 coupling dropped or mis-signed, wrong source), a system solved
    only approximately, an ansatz with a free direction.
    """
    alg, fS, free, Fs, P2, _, _ = solved(key)
    a, b = alg.a, alg.b
    assert sorted(free) == sorted(al for al in PX.SUBSETS if sum(al)) and len(free) == 15
    assert set(free.values()) == {0}
    h1 = 4 * a * b * (a * a + b * b) / (a * a - b * b) ** 2
    h2 = 8 * a * a * b * b / (a * a - b * b) ** 2
    assert all(P2[U] == 0 for U in PX.SUBSETS if sum(U) <= 1)
    assert P2[(1, 1, 0, 0)] == P2[(0, 0, 1, 1)] == h1 and P2[(0, 1, 1, 0)] == P2[(1, 0, 0, 1)] == -h1
    assert P2[(1, 0, 1, 0)] == P2[(0, 1, 0, 1)] == h2
    assert all(P2[U] != 0 for U in PX.SUBSETS if sum(U) == 3)        # the coupling reaches f_T with |T| = 1
    allowed = set(alg.keys11() + alg.keys_lc())
    n_coupled = 0
    for al in PX.SUBSETS:
        if sum(al) == 0:
            continue
        assert fS[al] and set(fS[al]) <= allowed
        assert _residual(alg, al, fS, Fs, P2) == {}, al
        if sum(al) >= 2:
            # the coupling is not vacuous: without it the residual is not zero
            assert PX.fadd(alg.L0(fS[al]), Fs[al]) != {}
            n_coupled += 1
    assert n_coupled == 11


@pytest.mark.parametrize("key", sorted(PAIRS))
def test_no_function_can_be_removed_from_the_span(key):
    """Removing any one of the 11 functions from the ansatz makes the top equation (S = all four exponents, with
    the true lower f_T) unsolvable.

    Rejects: a span stated larger than needed, and a solver that reports success whatever the ansatz.
    """
    alg, fS, _, Fs, P2, _, _ = solved(key)
    top = (1, 1, 1, 1)
    src = alg.rhs(top, Fs, P2, fS)
    full = alg.keys11() + alg.keys_lc()
    sol, nfree = alg.solve_L0(src, full)
    assert sol == fS[top] and nfree == 0
    for k in alg.keys11():
        with pytest.raises(RuntimeError):
            alg.solve_L0(src, [x for x in full if x != k])


# ---------------------------------------------------------------- (2') nine functions, rank nine
@pytest.mark.parametrize("key", sorted(PAIRS))
def test_K_lies_in_the_nine_function_span(key):
    """K, K', ..., K^(10) are Q[R, 1/R]-combinations of the nine functions: under e^{-2aR} the coordinates of
    Ei((a-b)R), Ei(2(a-b)R) and of the constant LC are tied as -2 : 1 : -1 (Pachucki's bracket of label 60), and
    likewise under e^{-2bR}; every one of the nine occurs in K.

    Rejects: a tie with the wrong ratio (+2 for -2), a tenth function, and the bare bracket without the constant
    (then the coefficient of e^{-2aR} is not in Q(R)).
    """
    alg, _, _, _, _, K, ders = solved(key)
    nine = nine_functions(alg)
    for i, A in enumerate(ders):
        rebuilt = {}
        for fn, lead in nine:
            for p, c in A.get(lead, {}).items():
                rebuilt = PX.fadd(rebuilt, PX.fscale(fn, c, p))
        assert rebuilt == A, f"derivative {i} leaves the nine-function span"
    for fn, lead in nine:
        assert K.get(lead), f"{lead} does not occur in K"
    # the constant is needed: LC e^{-2aR} occurs, with the coordinate of the bracket
    for x, y in ((alg.a, alg.b), (alg.b, alg.a)):
        assert K[("E", -2 * x, 1)] == {p: -c for p, c in K[("Ei", -2 * x, 2 * (x - y))].items()}


@pytest.mark.parametrize("key", sorted(PAIRS))
def test_module_generated_by_K_has_rank_nine(key):
    """rank{K, ..., K^(8)} = 9 at a rational point (so no operator of order <= 8 annihilates the closed form in the
    formal span), rank{K, ..., K^(10)} = 9, exactly one operator of order <= 9 and degree <= 16 up to scale, none
    of degree <= 15.

    Rejects: rank 8 or 10/11 (a missed or a spurious tie), a non-unique operator, a degree other than 16.
    """
    ders = solved(key)[6]
    keys = PX.all_keys(ders)
    assert len(keys) == 13                    # 11 functions + LC e^{-2aR}, LC e^{-2bR}
    for x in (Fraction(7, 5), Fraction(-3, 11)):
        M = PX.coordinate_matrix_at(ders, keys, x)
        assert PX.rank_q(M[:8]) == 8 and PX.rank_q(M[:9]) == 9 and PX.rank_q(M) == 9
    op = derived_operator(key)                # asserts uniqueness at (9, 16)
    assert len(op) == ORDER + 1 and [len(c) - 1 for c in op] == COEFF_DEGREES
    assert PX.operators_annihilating(ders, ORDER, DEGREE - 1) == []
    assert PX.operators_annihilating(ders, ORDER - 1, 40) == []


# ---------------------------------------------------------------- (3) the closed form is the exchange integral
def closed_form_value(key: str, Rv: str, dps: int):
    alg, K = solved(key)[0], solved(key)[5]
    with mp.workdps(dps):
        x = mp.mpf(Rv)
        return +PX.evaluate(alg, K, x, PX.f0_numeric(alg, x))


def _exchange_hp(key: str, Rv: str, tau_max: int, dps: int):
    from geovac import qfd_core as Q
    a, b = PAIRS[key]
    S1 = (1, 0, 0)
    return Q.exchange_hp(sp.Rational(str(a)), S1, S1, sp.Rational(str(b)), S1, S1, sp.Rational(Rv),
                         tau_max=tau_max, dps=dps)


@pytest.mark.parametrize("key,Rv", sorted(RECORD_DPS80))
def test_closed_form_equals_the_recorded_integral(key, Rv):
    """K = 16 a^3 b^3 / R * d^4 f against the records' exchange_hp values (dps 80, tau <= 32; 60 printed digits),
    to 55 digits.  The only quadrature is the one-dimensional tail integral for f0.

    Rejects: a wrong prefactor (8 for 16, a^2 b^4 for a^3 b^3), a closed form that solves the equations but is not
    the integral, a wrong LC.
    """
    with mp.workdps(70):
        v = closed_form_value(key, Rv, 70)
        ref = mp.mpf(RECORD_DPS80[(key, Rv)])
        assert abs(v - ref) < mp.mpf(10) ** -55 * ref, (mp.nstr(v, 62), RECORD_DPS80[(key, Rv)])


@pytest.mark.parametrize("key", sorted(PAIRS))
def test_closed_form_equals_live_exchange_hp_fast(key):
    """Closed form vs the live production route geovac.qfd_core.exchange_hp at R = 1.4, tau <= 8, dps 40: agreement
    to the size of the omitted tau tail (below 1e-15), which ties the typed record values above to the live code.

    Rejects: the same as the record test, at the 1e-15 level, without relying on typed digits.
    """
    with mp.workdps(40):
        v = closed_form_value(key, "1.4", 45)
        vh, per = _exchange_hp(key, "1.4", tau_max=8, dps=40)
        assert abs(per[-1]) < mp.mpf(10) ** -15
        assert abs(v - vh) < 100 * abs(per[-1]), (mp.nstr(v, 30), mp.nstr(vh, 30))


# tau_max chosen so that the omitted tau tail is below 1e-46 of the value (measured; asserted below)
LIVE = [("1_2", "1.4", 18), ("1_2", "2", 20), ("1_2", "3.015", 22),
        ("1_3half", "1.4", 15), ("1_3half", "2", 17), ("1_3half", "3.015", 18)]


@pytest.mark.slow
@pytest.mark.parametrize("key,Rv,tau", LIVE)
def test_closed_form_equals_live_exchange_hp_to_45_digits(key, Rv, tau):
    """Closed form vs live exchange_hp at R = 1.4, 2 and 3.015, dps 60: relative agreement to 45 digits or better
    (measured: 46 to 49 digits, limited by the tau truncation of the reference, not by the closed form).

    Rejects: a closed form right at one distance only (a wrong R-dependence), at 45 digits.
    """
    with mp.workdps(60):
        v = closed_form_value(key, Rv, 60)
        vh, per = _exchange_hp(key, Rv, tau_max=tau, dps=60)
        tail = abs(per[-1]) ** 2 / abs(per[-2])          # the tau terms fall geometrically; next omitted term
        assert abs(per[-1]) < abs(per[-2]) / 100 and tail < mp.mpf(10) ** -46 * abs(vh)
        assert abs(v - vh) < mp.mpf(10) ** -45 * abs(vh), (mp.nstr(v, 50), mp.nstr(vh, 50))
        assert abs(v - vh) < 100 * tail                  # the whole difference is the omitted tau tail


# ---------------------------------------------------------------- (4) the committed guessed operators
@pytest.mark.parametrize("name", sorted(SHA))
def test_committed_operator_files_are_the_scored_ones(name):
    assert hashlib.sha256(_bytes(name)).hexdigest() == SHA[name]


@pytest.mark.parametrize("key,name", FILE_CASES)
def test_committed_operator_annihilates_the_closed_form_mod_p(key, name):
    """The committed operator, applied to the nine-function representation of K with all arithmetic modulo its
    prime, has every coordinate (13 keys x all powers of R) equal to zero.

    Rejects: a guessed operator that does not annihilate the derived closed form (one coefficient off, the other
    pair's operator), and an application routine that returns "zero" for anything (control: one coefficient + 1).
    """
    cm, prime = committed(name)
    ders = solved(key)[6]
    assert len(cm) == ORDER + 1 and [len(c) - 1 for c in cm] == COEFF_DEGREES
    assert all(c[-1] % prime for c in cm)
    assert PX.apply_operator_mod(cm, ders, prime) == {}
    for i, j in ((0, 0), (5, 7), (9, 16)):
        bad = [list(c) for c in cm]
        bad[i][j] = (bad[i][j] + 1) % prime
        assert PX.apply_operator_mod(bad, ders, prime) != {}


@pytest.mark.parametrize("key,name", FILE_CASES)
def test_derived_operator_equals_the_committed_one_mod_p(key, name):
    """The operator built from the closed form (order 9, degree 16, coefficient degrees 14, 15, 16 x 8) equals the
    committed one modulo the prime, every coefficient, up to ONE scalar.

    Rejects: an operator that merely shares order and degree, or agrees up to a different scalar per coefficient.
    """
    cm, prime = committed(name)
    op = derived_operator(key)
    assert len(op) == ORDER + 1 and [len(c) - 1 for c in op] == COEFF_DEGREES
    assert op[-1][-1] % prime != 0
    s = cm[-1][-1] * pow(op[-1][-1], -1, prime) % prime
    assert [[x * s % prime for x in c] for c in op] == cm


@pytest.mark.parametrize("key", sorted(PAIRS))
def test_derived_operator_equals_the_reconstructed_operator_over_Q(key):
    """Stronger than the mod-p statement: the derived operator IS the CRT-reconstructed operator over Q frozen in
    tests/test_paper58_exchange_unequal.py (records EX, EXb), coefficient for coefficient (both primitive).

    Rejects: a rational reconstruction that is right modulo the two primes compared above but wrong over Q.
    """
    import test_paper58_exchange_unequal as U
    frozen = {"1_2": U.OP_EXCH_1_2, "1_3half": U.OP_EXCH_1_3HALF}[key]
    Rs = sp.Symbol("R")
    q = [[int(c) for c in reversed(sp.Poly(sp.sympify(s), Rs).all_coeffs())] for s in frozen]
    if q[-1][-1] < 0:
        q = [[-x for x in c] for c in q]
    assert derived_operator(key) == q


# ---------------------------------------------------------------- (5) singular data of the derived operator
@pytest.mark.parametrize("key", sorted(PAIRS))
def test_rates_at_infinity_and_exponents_at_zero(key):
    """Read from the DERIVED operator: characteristic roots at infinity {0 (x2), 2a, 2b, -(a+b), -2a (x2), -2b (x2)}
    (no rate +-(a-b)); R = 0 Fuchsian with exponents {-1 (x2), 0, 1, 2, 3, 4, 6, 9}.

    Rejects: a rate a - b, unsplit rates, a single 0 or -2a; exponent 5 for 6, a simple -1.
    """
    a, b = (sp.Rational(str(x)) for x in PAIRS[key])
    op = derived_operator(key)
    deg, rates = PX.rates_at_infinity(op)
    expected = {}
    for x, m in ((sp.Integer(0), 2), (2 * a, 1), (2 * b, 1), (-(a + b), 1), (-2 * a, 2), (-2 * b, 2)):
        expected[x] = expected.get(x, 0) + m
    assert deg == ORDER and rates == expected
    deg0, roots = PX.indicial_roots_at_zero(op)
    assert deg0 == ORDER and roots == {sp.Integer(k): m for k, m in R0_EXPONENTS.items()}


# ---------------------------------------------------------------------------------------------------------------
# TIMINGS (2026-11-18).  Fast tier 33 tests, about 66 s (live exchange_hp tau <= 8: 16-18 s per pair; the symbolic
# exponent derivatives 5-8 s per pair, once; everything else under 2.5 s).  Slow tier 6 tests, about 18 min (--slow;
# live exchange_hp at dps 60, 80-350 s per point).
#
# DIGITS REACHED (closed form at dps 60 vs live exchange_hp; relative difference, which equals the omitted tau tail
# of the reference in every case):  (1,2): R = 1.4 2.4e-49, R = 2 3.3e-50, R = 3.015 4.2e-49;
# (1,3/2): R = 1.4 below 1e-45 at tau <= 15 (6.2e-44 at tau <= 14), R = 2 below 1e-45 at tau <= 17 (9.1e-47 at
# tau <= 16), R = 3.015 8.3e-48.  Against the records' dps-80 values (R = 1.4, 3.015): 1e-59 to 1e-62, the limit of
# the 60 printed digits.
#
# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper58_pachucki_exchange.py --plant "OLD=>NEW"
# [--plant-in tests/_pachucki_exchange.py] -k "<sel> and not 45_digits"), written and fired as a separate step after
# the guard.  H = planted in the helper, T = in this file.  41 runs, 39 fired; the 2 that did not are explained.
#   drop the p^2 coupling term (H)          "src = fadd(src, fscale(fS[be], P2[U], 1))" -> "src = src"
#                                           -k fifteen / recorded / annihilates / equals_the_committed    FIRED x4
#   coupling sign (H)                       "P2[U], 1))" -> "-P2[U], 1))"            -k fifteen           FIRED
#   prefactor 8 (H)                         "PREFACTOR = 16" -> "= 8"    -k recorded / live_exchange_hp_fast   FIRED x2
#   prefactor a^2 b^4 (H)                   "self.a**3 * self.b**3" -> "self.a**2 * self.b**4"   same     FIRED x2
#   label 60: -2 Ei -> +2 Ei (H)            '("Ei", -2 * a, (a - b)): {0: -2 * c}' -> "2 * c"
#                                           -k label60 / nine_function                                    FIRED x2
#   tie -2 -> +2 (T)                        "TIE = -2" -> "TIE = 2"                  -k nine_function     FIRED
#   remove e^{-(a+b)R} from the span (H)    keys11                                   -k fifteen           FIRED
#   remove E1((a+b)R) from the span (H)     keys11                                   -k fifteen           FIRED
#   remove LC e^{-2aR} from the span (H)    keys_lc                                  -k fifteen           FIRED
#   swap the (1,2) and (1,3/2) operators (T)  "for name in OPS[key]]" -> swapped key
#                                           -k annihilates / equals_the_committed                         FIRED x2
#   pinned hash altered (T)                 "e84e0e89" -> "e84e0e88"                 -k scored_ones       FIRED
#   Hessian a^2+b^2 -> a^2-b^2 (T)          h1                                       -k p2_vanishes       FIRED
#   Hessian value 40/9 -> 32/9 (T)                                                   -k p2_vanishes       FIRED
#   p^2 mistranscribed (H, label 20)        second factor u2^2 w3^2 - u3^2 w2^2      -k p2_vanishes       FIRED
#   p^2 first-factor sign (H, label 20)                                 -k p2_vanishes / fifteen          FIRED x2
#   label 33 Ei sign (H)                    "+ Ei(r * (w2 - w3))" -> "-"             -k label60           FIRED
#   a+ sign (H, label 22)                   "+ u3" -> "- u3" in ap      -k recorded / annihilates / equals_the_committed   FIRED x3
#   same                                                                             -k label60           DID NOT FIRE
#   same                                                                             -k fifteen           DID NOT FIRE
#       Explained: a+ = 0 and F_1^- = 0 at the exchange point, so label 30 there does not see a+ (label60 test), and
#       a wrong source is still solved exactly in the same span (fifteen test: it checks that the equations as
#       transcribed are solved, not that they are Pachucki's).  The transcription of the parts of F that vanish at
#       the exchange point is guarded only end to end: by the numeric tests and by the operator comparisons.
#   exponent 6 -> 5 (T)                     R0_EXPONENTS                             -k rates             FIRED
#   rate -(a+b) -> a-b (T)                                                           -k rates             FIRED
#   degrees [14, 15] -> [15, 15] (T)        COEFF_DEGREES                            -k rank_nine         FIRED
#   order 9 -> 8 (T)                        "ORDER, DEGREE = 9, 16" -> "8, 16"       -k rank_nine         FIRED
#   record value, digit 50 altered (T)      "...4556376580789..." -> "...4556376581789..."   -k recorded  FIRED
#   Q operator p_0 coefficient + 1          in tests/test_paper58_exchange_unequal.py   -k reconstructed_operator_over_Q   FIRED
#   prefactor 8 (H), slow tier              -k "45_digits and 1_3half and 1.4"                            FIRED
#   committed operator, one coefficient + 1, on a TEMP COPY of the four files (DATA repointed by a pytest plugin;
#   debug/data untouched):  hash pin as is  -k scored_ones                                                FIRED
#                           hash re-pinned  -k annihilates / equals_the_committed                         FIRED x2
#                           (control: unmodified copy, all pass)
