"""Regression backing for the migrated LiH R12-CI engine (geovac/lih_r12ci).

Pins the two-center four-electron explicit-r12 CI energy for both geminals against the
validated debug-tree / VMC ground truths (Paper 12 Sec. "Explicit correlation"; CHANGELOG
v5.15.10-.14):

  * linexp (cusp-correct f = r e^{-g r}): E_R12 = -7.9168 Ha, dE = -29.0 mHa, matching the
    independent VMC correlation lowering -29.5 mHa. This is the item-(1) result -- the
    Cov[F,Y_sum] "grid-limited" miss closed by exact-isotropic Yukawa dressings (Cov[F,Y]
    moved -0.079 -> -0.1025, MC -0.1018); the guard here would FAIL on the pre-fix -0.079.
  * exp (f = e^{-g r}): E_R12 = -7.9420 Ha, matching the same-geminal VMC 2x2 -7.94121.

PoC (ionic single-zeta reference, not spectroscopic LiH): the deliverable is the RI-free
reduction reproduced from the production package, not the correlation energy recovered.

Slow: the first energy() call builds the shared prolate grids/kernels (~2 min); both geminals
share the build. Run with:  pytest tests/test_lih_r12ci.py --slow
"""
import pytest

from geovac.lih_r12ci import energy, VMC_TARGETS
from geovac.lih_r12ci.assembly import ANALYTIC_REF

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def results():
    return {"exp": energy("exp"), "linexp": energy("linexp")}


def test_linexp_lands_the_cusp_correct_headline(results):
    """linexp E_R12 = -7.9168 (dE -29.0 mHa), variational, matching the VMC dE -29.5."""
    r = results["linexp"]
    assert r.variational, "linexp not variational"
    assert r.E_R12 < r.E0, "no correlation lowering"
    # absolute analytic value (deterministic, RI-free)
    assert abs(r.E_R12 - (-7.9168)) < 5e-4, f"linexp E_R12 {r.E_R12} != -7.9168"
    # the apples-to-apples comparison is the correlation lowering vs the VMC dE
    assert abs(r.dE_mHa - VMC_TARGETS["linexp"]["dE_mHa"]) < 1.0, (
        f"linexp dE {r.dE_mHa:.2f} mHa not within 1 mHa of VMC {VMC_TARGETS['linexp']['dE_mHa']}")


def test_linexp_covfy_is_the_fixed_value_not_the_grid_limited_one(results):
    """The item-(1) fix: exact-isotropic Yukawa dressings give Cov[F,Y] ~ -0.1025 (MC -0.1018).
    This guard FAILS on the pre-fix psi_yuk value -0.079 -- it pins the fix, not just the sign."""
    covfy = results["linexp"].pieces["CovFY"]
    assert covfy < -0.098, f"Cov[F,Y] {covfy:.5f} not at the fixed value (pre-fix was -0.079)"
    assert covfy > -0.106, f"Cov[F,Y] {covfy:.5f} overshot the MC target -0.1018"


def test_exp_reproduces_the_validated_2x2(results):
    """exp E_R12 = -7.9420 (matches the same-geminal VMC 2x2 -7.94121), variational.
    The exp path keeps psi_yuk dressings by design (validated cancellation, NOT the item-1 fix)."""
    r = results["exp"]
    assert r.variational, "exp not variational"
    assert abs(r.E_R12 - (-7.9420)) < 5e-4, f"exp E_R12 {r.E_R12} != -7.9420"
    assert abs(r.E_R12 - VMC_TARGETS["exp"]["E_R12"]) < 1.0e-3, (
        f"exp E_R12 {r.E_R12} not within 1 mHa of VMC 2x2 {VMC_TARGETS['exp']['E_R12']}")


@pytest.mark.parametrize("gem", ["exp", "linexp"])
def test_pieces_match_analytic_reference(results, gem):
    """Every matrix-element component matches the debug-tree analytic reference."""
    r = results[gem]; ref = ANALYTIC_REF[gem]
    assert abs(r.sigma2 - ref["sigma2"]) < 1e-4, f"{gem} sigma2 {r.sigma2} vs {ref['sigma2']}"
    # h and g reproduced to sub-mHa (exp allows the deterministic-ab vs MC-ab drift)
    assert abs(r.h - ref["h"]) < 1e-3, f"{gem} h {r.h} vs {ref['h']}"
    assert abs(r.g - ref["g"]) < 2e-3, f"{gem} g {r.g} vs {ref['g']}"
    for key in ("h_T", "h_Vne", "h_Vee", "g_T", "g_Vee"):
        assert abs(r.pieces[key] - ref[key]) < 2e-3, (
            f"{gem} {key} {r.pieces[key]:.5f} vs ref {ref[key]:.5f}")


def test_overlap_is_diagonal_and_reference_is_variational(results):
    """S01 = <F-Fbar> = 0 exactly (diagonal overlap), and E0 > exact -8.070."""
    for gem in ("exp", "linexp"):
        r = results[gem]
        assert r.E0 > -8.070, f"{gem} E0 {r.E0} below exact (bug)"
        assert -8.070 < r.E_R12 < r.E0, f"{gem} E_R12 {r.E_R12} out of (exact, E0)"
