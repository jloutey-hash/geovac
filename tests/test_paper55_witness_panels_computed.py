"""Paper 55 M1 / M2 witness panels and the S^(3) anchor -- COMPUTED, then compared.

Guard (group3 guard pass, 2026-11-03).  The earlier panels in
test_paper55_m1_pure_tate.py / test_paper55_M2_mixed_tate.py typed the panel
value in and PSLQ'd it (plants 17 pi^8 -> 16 pi^8 and 4/pi -> 3/pi did not
fire).  Here every witness is computed from its spectral source, by a route
that does not use the closed form, and only then compared with the panel
entry as printed in Paper 55 (the PANEL dicts below are the paper's table
cells, the only literals in the file).

M1 panel (Paper 55 sec:m1 Witnesses)
  Vol(S^2)/4 = pi                    -- quadrature of the area element
  4/pi Fejer-moment constant         -- lim n gamma_n / log n, fitted from the
                                        exact sum rule at n = 40..448; the
                                        quadrature route gamma_rate agrees
  2/pi on the unit S^3               -- same moment with distance chi/2
  1/pi^2 F-theorem cross-product     -- (F_D - 2 F_s) / (3 zeta(3)/4) from the
                                        Hurwitz-form production F values
  (1/pi^2 Pythagorean prefactor is posited, "class only": not computable.)

M2 panel (Paper 55 sec:m2 Witnesses)
  a_0, a_1 (Dirac, vol-norm), a_k = 0 (k >= 2) -- small-t fit of the CH heat
      trace sum_n 2(n+1)(n+2) e^{-t (n+3/2)^2}, times (4 pi)^{3/2}
  a_k scalar = 2 pi^2/k!             -- small-t fit of sum_n n^2 e^{-t(n^2-1)}
  zeta_{D^2}(s), s = 2..5            -- direct spectral sum (mpmath nsum)

S^(3) (Paper 55 tab. depth-k; Paper 28)
  The factorised-chain definition S^(3) = sum_m phi(m) L(m)^2 summed directly
  to N = 48000 and extrapolated in the Euler-Maclaurin basis
  {ln^k N / N^j}; pins the canonical 31.5725612075120... to < 1e-11, where the
  previous definitional tie was the rigorous bracket [31.57063, 31.57300]
  (about 4 digits).  Independent of the t-value decomposition.
"""
from __future__ import annotations

import mpmath as mp
import pytest
import sympy as sp

from geovac.central_fejer_su2 import gamma_n_via_sum_rule, gamma_rate
from geovac.qed_two_loop import dirac_F_theorem, scalar_F_theorem

# ---- Paper 55 table cells (sympy strings, exactly as printed) -------------
PANEL_M1 = {
    "Vol_S2_over_4": "pi",
    "fejer_dual_coxeter": "4/pi",
    "fejer_unit_S3": "2/pi",
    "F_theorem_cross_product": "1/pi**2",
}
PANEL_M2_DIRAC = {0: "4*pi**2", 1: "-2*pi**2", 2: "0", 3: "0"}
PANEL_M2_SCALAR = "2*pi**2/factorial(k)"
PANEL_M2_ZETA = {
    2: "pi**2 - pi**4/12",
    3: "pi**4/3 - pi**6/30",
    4: "2*pi**6/15 - 17*pi**8/1260",
    5: "pi**8*(306 - 31*pi**2)/5670",
}
# Paper 28 / Paper 55 canonical S^(3) (first 30 digits)
S3_CANONICAL = "31.5725612075120227547619295450"


def _val(expr: str, **subs) -> mp.mpf:
    e = sp.sympify(expr).subs(subs)
    return mp.mpf(sp.N(e, mp.mp.dps + 10))


# ---------------------------------------------------------------------------
# M1
# ---------------------------------------------------------------------------

def test_m1_volume_slot_from_area_quadrature():
    mp.mp.dps = 40
    area = mp.quad(lambda th, ph: mp.sin(th), [0, mp.pi], [0, 2 * mp.pi])
    assert abs(area / 4 - _val(PANEL_M1["Vol_S2_over_4"])) < mp.mpf(10) ** -35


@pytest.fixture(scope="module")
def fejer_fit():
    ns = [40, 56, 80, 112, 160, 224, 320, 448]
    g = {n: gamma_n_via_sum_rule(n, prec=30) for n in ns}
    mp.mp.dps = 30

    def basis(n):
        n = mp.mpf(n)
        l = mp.log(n)
        return [l, 1, l / n, 1 / n, l / n ** 2, 1 / n ** 2, l / n ** 3, 1 / n ** 3]

    sol = mp.lu_solve(mp.matrix([basis(n) for n in ns]), mp.matrix([n * g[n] for n in ns]))
    return sol[0], sol[1]


def test_m1_fejer_sum_rule_matches_quadrature_route():
    """The exact sum rule used for the fit agrees with direct quadrature of
    the Fejer kernel moment (two routes to gamma_n)."""
    for n in (3, 6):
        a = gamma_n_via_sum_rule(n, prec=30)
        b = gamma_rate(n, prec=30)
        assert abs(a - b) < mp.mpf(10) ** -20


def test_m1_fejer_moment_constant_dual_coxeter(fejer_fit):
    slope, intercept = fejer_fit
    mp.mp.dps = 30
    assert abs(slope - _val(PANEL_M1["fejer_dual_coxeter"])) < 1e-6
    # Paper 18 MR-C next-order constant c = 4.10932146748987...
    assert abs(intercept - mp.mpf("4.1093214674898756204229")) < 1e-5


def test_m1_fejer_moment_constant_unit_S3(fejer_fit):
    """On the unit S^3 the geodesic distance is chi/2: the moment halves."""
    slope, _ = fejer_fit
    assert abs(slope / 2 - _val(PANEL_M1["fejer_unit_S3"])) < 1e-6


def test_m1_F_theorem_cross_product_factor():
    mp.mp.dps = 50
    m1 = (dirac_F_theorem() - 2 * scalar_F_theorem()) / (3 * mp.zeta(3) / 4)
    assert abs(m1 - _val(PANEL_M1["F_theorem_cross_product"])) < mp.mpf(10) ** -40


# ---------------------------------------------------------------------------
# M2
# ---------------------------------------------------------------------------

def _small_t_coeffs(K, n_coef: int, t_lo: str, t_hi: str):
    """Solve K(t) t^{3/2} = sum_{k<n_coef} c_k t^k at n_coef small t-points.

    The exponentially small modular remainder (e^{-pi^2/t}) is far below the
    working precision on the chosen t-range.
    """
    lo, hi = mp.mpf(t_lo), mp.mpf(t_hi)
    ts = [lo + (hi - lo) * i / (n_coef - 1) for i in range(n_coef)]
    M = mp.matrix([[t ** k for k in range(n_coef)] for t in ts])
    b = mp.matrix([K(t) * t ** mp.mpf(1.5) for t in ts])
    return mp.lu_solve(M, b)


def _K_dirac(t):
    tot = mp.mpf(0)
    n = 0
    while True:
        term = 2 * (n + 1) * (n + 2) * mp.exp(-t * (n + mp.mpf(3) / 2) ** 2)
        tot += term
        if n > 10 and term < mp.mpf(10) ** (-mp.mp.dps - 5):
            return tot
        n += 1


def _K_scalar(t):
    tot = mp.mpf(0)
    n = 1
    while True:
        term = n * n * mp.exp(-t * (n * n - 1))
        tot += term
        if n > 10 and term < mp.mpf(10) ** (-mp.mp.dps - 5) * tot:
            return tot
        n += 1


def test_m2_dirac_seeley_dewitt_from_heat_trace():
    mp.mp.dps = 60
    c = _small_t_coeffs(_K_dirac, 6, "0.01", "0.03")
    vol = (4 * mp.pi) ** mp.mpf(1.5)
    for k in range(4):
        assert abs(c[k] * vol - _val(PANEL_M2_DIRAC[k])) < mp.mpf(10) ** -15, (k, c[k] * vol)


def test_m2_scalar_seeley_dewitt_from_heat_trace():
    mp.mp.dps = 80
    c = _small_t_coeffs(_K_scalar, 14, "0.01", "0.05")
    vol = (4 * mp.pi) ** mp.mpf(1.5)
    for k in range(5):
        target = _val(PANEL_M2_SCALAR.replace("factorial(k)", f"factorial({k})"))
        assert abs(c[k] * vol - target) < mp.mpf(10) ** -12, (k, c[k] * vol, target)


@pytest.mark.parametrize("s", [2, 3, 4, 5])
def test_m2_spectral_zeta_direct_sum(s):
    mp.mp.dps = 40
    direct = mp.nsum(lambda n: 2 * (n + 1) * (n + 2) / (n + mp.mpf(3) / 2) ** (2 * s),
                     [0, mp.inf])
    assert abs(direct - _val(PANEL_M2_ZETA[s])) < mp.mpf(10) ** -30, (s, direct)


# ---------------------------------------------------------------------------
# S^(3) anchor
# ---------------------------------------------------------------------------

def test_s3_definition_pins_canonical_value():
    mp.mp.dps = 40
    pi = mp.pi
    # P = sum_{n>=1} phi(n) = zeta_{D^2}(2) - phi(0), phi(0) = 64/81
    P = pi ** 2 - pi ** 4 / 12 - mp.mpf(64) / 81
    checkpoints = [2000, 3000, 4000, 6000, 8000, 12000, 16000, 24000, 32000, 48000]
    head = mp.mpf(0)
    cum_phi = mp.mpf(0)
    cum_nphi = mp.mpf(0)
    S = {}
    for m in range(1, checkpoints[-1] + 1):
        x2 = mp.mpf(2 * m + 3) ** 2
        pm = 8 * (x2 - 1) / (x2 * x2)
        A = cum_nphi + m * (P - cum_phi)       # sum_n min(n, m) phi(n), full range
        L = 2 * A - P - pm
        head += pm * L * L
        cum_phi += pm
        cum_nphi += m * pm
        if m in checkpoints:
            S[m] = +head

    def basis(N):
        N = mp.mpf(N)
        l = mp.log(N)
        out = [mp.mpf(1)]
        for j in (1, 2, 3):
            out += [l ** 2 / N ** j, l / N ** j, 1 / N ** j]
        return out

    sol = mp.lu_solve(mp.matrix([basis(N) for N in checkpoints]),
                      mp.matrix([S[N] for N in checkpoints]))
    extrap = sol[0]
    assert abs(extrap - mp.mpf(S3_CANONICAL)) < mp.mpf(10) ** -11, mp.nstr(extrap, 20)
    # S_N = S_inf - tail: the leading tail term is the analytic 32 ln^2 N / N
    # (phi ~ 2/m^2, L ~ 4 ln m), so its fitted coefficient is -32.
    assert abs(sol[1] + 32) < 1e-3
