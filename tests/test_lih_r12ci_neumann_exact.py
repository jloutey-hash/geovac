"""Fast guard for the exact ordered-integral prolate-Neumann operator (Phase 0b of the LiH marriage
build, 2026-09-23; geovac/lih_r12ci/neumann_exact.py, switch geovac.lih_r12ci.kernels.USE_EXACT_NEUMANN).

What each guard EXCLUDES (the wrong answer it was fire-tested against, debug/firetest_lih_r12ci_neumann_exact.py):

  test_exact_operator_reproduces_5zeta_over_8
      The legacy cumulative-Gauss-Legendre radial step, which integrates the kinked kernel
      P_l(xi_<) Q_l(xi_>) across its diagonal and is wrong by +5.0e-3 (zeta = 2.6875) and +8.5e-3
      (zeta = 4.5) relative on the Li 1s self-Coulomb -- the +8.09 mHa defect Phase 0 found in
      <V_ee> of the CI reference.  Also fires if the switch is dead (exact=True silently falling
      through to the legacy branch), if P and Q are interchanged in the operator, or if the sub-rule
      is too coarse to integrate the degree <= 109 P-side polynomial exactly.
  test_legacy_operator_is_wrong_at_the_pinned_level
      The built-in fire test of the guard above: the SAME closed form applied to the legacy path must
      FAIL by more than 1e-3.  Excludes a rewrite that would make the two paths agree by degrading the
      exact one (or a "fix" of the default path that would silently move the banked -7.9168 / -7.9420).
  test_two_centre_aabb_matches_the_hartree_closed_form
      Wrong two-centre normalisation / prefactor (2/R)(2 pi) a^3 of the potential: (aa|bb) is anchored
      to the closed-form Hartree potential of a 1s density (energy._hartree_1s), in BOTH dressing
      directions.  NOTE energy.V_aabb is the f-GEMINAL integral (aa|f|bb), not this Coulomb integral.
  test_general_m_mode_potential_matches_the_solid_harmonic_closed_form
      The general-m path of triangle.coul_mode_potential (the triangle reduction's Coulomb): the
      m = 1 mode density rho_cyl e^{-2 zeta r} is the cos(phi) component of a degree-1 solid-harmonic
      density whose potential is closed-form; the legacy m = 1 path is off by 8e-3 sup-norm and must
      fail the same anchor.  Excludes a wrong (x^2-1)^{m/2} handling or a wrong P_l^m / Q_l^m table.
  test_default_switch_is_off_and_legacy_path_is_untouched
      A default of USE_EXACT_NEUMANN = True (or a consumer reading the flag by from-import at import
      time) would change the legacy anchors of tests/test_lih_r12ci.py; the default path must be the
      explicit exact=False path bit-for-bit.

All anchors are closed forms; no debug/ artifact is read.  Wall < 10 s (imports hVee/triangle only;
gVee's phi-kernels are NOT built).
"""
import math

import numpy as np
import pytest
from scipy.special import gammainc

import geovac.lih_r12ci.kernels as KG
from geovac.lih_r12ci.hVee import neumann_potential, d_aa, d_bb, NXI, NETA
from geovac.lih_r12ci.triangle import coul_mode_potential
from geovac.lih_r12ci.energy import _hartree_1s, ZA, ZB, a

SHAPE = (NXI, NETA)
RA = KG.rA.ravel()
RB = KG.rB.ravel()
RHO_CYL = KG.RHO_CYL.ravel()


def _rho_1s(zeta: float, r: np.ndarray) -> np.ndarray:
    return (zeta ** 3 / np.pi) * np.exp(-2.0 * zeta * r)


def _self_coulomb(zeta: float, r: np.ndarray, exact: bool) -> float:
    rho = _rho_1s(zeta, r)
    V = neumann_potential(rho.reshape(SHAPE), exact=exact).reshape(-1)
    return float(KG.grid_int(rho * V))


def _solid_harmonic_potential(r: np.ndarray, zeta: float, m: int) -> np.ndarray:
    """cos(m phi) coefficient of the Coulomb potential of rho_cyl^m cos(m phi) e^{-2 zeta r}, times
    2 pi a^3 (the normalisation triangle.coul_mode_potential carries)."""
    s = 2 * m + 3
    lower = gammainc(s, 2.0 * zeta * r) * math.factorial(s - 1) / (2.0 * zeta) ** s
    upper = np.exp(-2.0 * zeta * r) * (r / (2.0 * zeta) + 1.0 / (2.0 * zeta) ** 2)
    rr = np.maximum(r, 1e-300)
    return 2 * np.pi * a ** 3 * (4.0 * np.pi / (2 * m + 1)) * RHO_CYL ** m * (lower / rr ** (2 * m + 1) + upper)


@pytest.mark.parametrize("zeta", [4.5, 2.6875])
def test_exact_operator_reproduces_5zeta_over_8(zeta):
    """Exact ordered-integral operator: Neumann self-Coulomb of a 1s(zeta) STO on centre A == 5 zeta/8."""
    J = _self_coulomb(zeta, RA, exact=True)
    ref = 5.0 * zeta / 8.0
    assert abs(J - ref) / ref <= 1e-7, f"zeta={zeta}: exact operator {J:.10f} vs 5zeta/8 {ref:.10f}"


@pytest.mark.parametrize("zeta", [4.5, 2.6875])
def test_legacy_operator_is_wrong_at_the_pinned_level(zeta):
    """The wrong answer the guard above excludes: the legacy cumulative-GL path misses 5 zeta/8 by > 1e-3."""
    J = _self_coulomb(zeta, RA, exact=False)
    ref = 5.0 * zeta / 8.0
    err = abs(J - ref) / ref
    assert err > 1e-3, f"zeta={zeta}: legacy path is no longer wrong ({err:.2e}) -- the pinned defect moved"
    assert err < 2e-2, f"zeta={zeta}: legacy error {err:.2e} is not the documented O(1e-3..1e-2) level"


def test_two_centre_aabb_matches_the_hartree_closed_form():
    """(aa|bb) Coulomb integral, both dressing directions, vs the closed-form Hartree-potential anchor."""
    anchor = float(KG.grid_int(d_bb * _hartree_1s(RA, ZA)))
    anchor_alt = float(KG.grid_int(d_aa * _hartree_1s(RB, ZB)))
    assert abs(anchor - anchor_alt) / anchor < 1e-12
    V_aa = neumann_potential(d_aa.reshape(SHAPE), exact=True).reshape(-1)
    V_bb = neumann_potential(d_bb.reshape(SHAPE), exact=True).reshape(-1)
    dress_iso = float(KG.grid_int(d_bb * V_aa))
    dress_other = float(KG.grid_int(d_aa * V_bb))
    assert abs(dress_iso - anchor) / anchor <= 1e-7, f"(aa|bb) dress rho_aa: {dress_iso:.10f} vs {anchor:.10f}"
    assert abs(dress_other - anchor) / anchor <= 1e-7, f"(aa|bb) dress rho_bb: {dress_other:.10f} vs {anchor:.10f}"


def test_general_m_mode_potential_matches_the_solid_harmonic_closed_form():
    """m = 1 mode potential (triangle path) on the exact operator == closed form; the legacy m = 1 path does not."""
    zeta, m = 2.6875, 1
    Bm = RHO_CYL ** m * np.exp(-2.0 * zeta * RA)
    Va = _solid_harmonic_potential(RA, zeta, m)
    Ve = coul_mode_potential(Bm.reshape(SHAPE), m, exact=True).reshape(-1)
    Vl = coul_mode_potential(Bm.reshape(SHAPE), m, exact=False).reshape(-1)
    sup = np.max(np.abs(Va))
    err_exact = np.max(np.abs(Ve - Va)) / sup
    err_legacy = np.max(np.abs(Vl - Va)) / sup
    assert err_exact <= 1e-8, f"m=1 exact mode potential sup-rel {err_exact:.2e}"
    assert err_legacy > 1e-3, f"m=1 legacy mode potential is no longer wrong ({err_legacy:.2e})"


def test_default_switch_is_off_and_legacy_path_is_untouched():
    """Default = legacy, bit-for-bit: the banked -7.9168 / -7.9420 (tests/test_lih_r12ci.py) rest on it."""
    assert KG.USE_EXACT_NEUMANN is False
    rho = _rho_1s(4.5, RA).reshape(SHAPE)
    V_default = neumann_potential(rho)
    V_legacy = neumann_potential(rho, exact=False)
    V_exact = neumann_potential(rho, exact=True)
    assert np.array_equal(V_default, V_legacy), "default path is not the legacy path"
    assert not np.array_equal(V_default, V_exact), "exact and legacy paths coincide -- the switch does nothing"
    Bm = RHO_CYL * np.exp(-2.0 * 2.6875 * RA)
    assert np.array_equal(coul_mode_potential(Bm.reshape(SHAPE), 1), coul_mode_potential(Bm.reshape(SHAPE), 1, exact=False))
