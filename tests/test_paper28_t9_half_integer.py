"""Paper 28 Theorem T9 scope: integer s is pure-Tate, half-integer s carries odd zeta.

T9: zeta_{D^2}(s) = 2^{2s-1} [lambda(2s-2) - lambda(2s)],
lambda(x) = (1 - 2^{-x}) zeta(x), for the Camporesi-Higuchi spectrum on unit S^3
(|lambda_n| = n + 3/2, total multiplicity 2(n+1)(n+2)).

(a) Integer s: zeta_{D^2}(2) = pi^2 - pi^4/12, zeta_{D^2}(3) = pi^4/3 - pi^6/30,
    each matching the direct spectral sum to 1e-20.
(b) Half-integer s: zeta_{D^2}(5/2) = 14 zeta(3) - (31/2) zeta(5) against the direct
    spectral sum to 1e-20. This is the claim that retired "no odd-zeta content at any s"
    (DELTA #20). Wrong answer it rejects: any other rational coefficient on zeta(5)
    (fire-tested with 31/2 -> 15/2), and any sqrt(pi)-multiple reading.
"""
import mpmath
import pytest



@pytest.fixture
def hp():
    old = mpmath.mp.dps
    mpmath.mp.dps = 40
    try:
        yield
    finally:
        mpmath.mp.dps = old


def direct_sum(s):
    """Direct spectral sum sum_n 2(n+1)(n+2) (n+3/2)^{-2s}, via Hurwitz zeta.

    2(n+1)(n+2) = 2 (n+3/2)^2 - 1/2, so the sum equals
    2 zeta(2s-2, 3/2) - (1/2) zeta(2s, 3/2) exactly; cross-checked against a
    brute partial sum + Euler-Maclaurin tail by mpmath.nsum below.
    """
    s = mpmath.mpf(s)
    return 2 * mpmath.zeta(2 * s - 2, mpmath.mpf(3) / 2) - mpmath.zeta(2 * s, mpmath.mpf(3) / 2) / 2


def brute_sum(s):
    s = mpmath.mpf(s)
    return mpmath.nsum(lambda n: 2 * (n + 1) * (n + 2) * (n + mpmath.mpf(3) / 2) ** (-2 * s), [0, mpmath.inf])


def t9_formula(s):
    s = mpmath.mpf(s)
    lam = lambda x: (1 - mpmath.power(2, -x)) * mpmath.zeta(x)
    return mpmath.power(2, 2 * s - 1) * (lam(2 * s - 2) - lam(2 * s))


def test_direct_sum_routes_agree(hp):
    for s in (3, mpmath.mpf(5) / 2):
        assert abs(direct_sum(s) - brute_sum(s)) < mpmath.mpf(10) ** -20


def test_integer_s_pure_tate(hp):
    pi = mpmath.pi
    closed = {2: pi ** 2 - pi ** 4 / 12, 3: pi ** 4 / 3 - pi ** 6 / 30}
    for s, val in closed.items():
        assert abs(t9_formula(s) - val) < mpmath.mpf(10) ** -20
        assert abs(direct_sum(s) - val) < mpmath.mpf(10) ** -20


def test_half_integer_odd_zeta(hp):
    s = mpmath.mpf(5) / 2
    closed = 14 * mpmath.zeta(3) - mpmath.mpf(31) / 2 * mpmath.zeta(5)
    assert abs(direct_sum(s) - closed) < mpmath.mpf(10) ** -20
    assert abs(t9_formula(s) - closed) < mpmath.mpf(10) ** -20
    assert abs(closed - mpmath.mpf("0.7564164395")) < mpmath.mpf(10) ** -10


def test_poles_and_heat_trace_coefficients(hp):
    """Paper 28 remark "Where sqrt(pi) lives": the heat trace of D^2 (|lambda| = n+3/2,
    deg 2(n+1)(n+2)) is sqrt(pi)/2 t^{-3/2} - sqrt(pi)/4 t^{-1/2} up to exponentially small
    terms, so zeta_{D^2} has simple poles at s = 3/2, 1/2 with RATIONAL residues 1 and -1/4;
    volume-normalised (x (4 pi)^{3/2}) the coefficients are a_0 = 4 pi^2, a_1 = -2 pi^2.
    Rejects a sqrt(pi)-valued residue and a wrong heat-trace coefficient."""
    t = mpmath.mpf(1) / 50
    tr = mpmath.nsum(lambda n: 2 * (n + 1) * (n + 2) * mpmath.e ** (-t * (n + mpmath.mpf(3) / 2) ** 2), [0, mpmath.inf])
    c0, c1 = mpmath.sqrt(mpmath.pi) / 2, -mpmath.sqrt(mpmath.pi) / 4
    assert abs(tr - (c0 * t ** mpmath.mpf(-1.5) + c1 * t ** mpmath.mpf(-0.5))) < mpmath.mpf(10) ** -30
    norm = (4 * mpmath.pi) ** mpmath.mpf(1.5)
    assert abs(c0 * norm - 4 * mpmath.pi ** 2) < mpmath.mpf(10) ** -30
    assert abs(c1 * norm + 2 * mpmath.pi ** 2) < mpmath.mpf(10) ** -30
    lam = lambda x: (1 - mpmath.mpf(2) ** (-x)) * mpmath.zeta(x)
    f = lambda s: mpmath.mpf(2) ** (2 * s - 1) * (lam(2 * s - 2) - lam(2 * s))
    eps = mpmath.mpf(10) ** -20
    assert abs(eps * f(mpmath.mpf(3) / 2 + eps) - 1) < mpmath.mpf(10) ** -15
    assert abs(eps * f(mpmath.mpf(1) / 2 + eps) + mpmath.mpf(1) / 4) < mpmath.mpf(10) ** -15
    assert abs(c0 / mpmath.gamma(mpmath.mpf(3) / 2) - 1) < mpmath.mpf(10) ** -30
    assert abs(c1 / mpmath.gamma(mpmath.mpf(1) / 2) + mpmath.mpf(1) / 4) < mpmath.mpf(10) ** -30
