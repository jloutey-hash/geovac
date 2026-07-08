"""Backing test for Paper 19, "Depth versus shape: the curvature read-out".

Claim: the balanced-coupled LiH residual is a well-SHAPE error, not a depth
error. Even though the energy converges (0.20% at n_max=3), the well is too
stiff: a harmonic fit to the n_max=2 balanced PES about its minimum gives
omega_e ~ 2040 cm^-1, ~45% above the experimental 1406 cm^-1 (force constant
~0.14 vs the experimental ~0.066 Ha/bohr^2), and R_eq drifts outward to ~3.22.

Slow (~1-2 min, live solver at n_grid_vne=8000). Run with `pytest --slow`.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import lih_spec

# LiH X^1Sigma+ experimental constants (Huber-Herzberg)
CM1_TO_HA = 1.0 / 219474.6313702
MU_ME = (7.016004 * 1.007825) / (7.016004 + 1.007825) * 1822.888486  # electron masses
W_EXP_CM1 = 1405.65
K_TRUE = MU_ME * (W_EXP_CM1 * CM1_TO_HA) ** 2  # ~0.0659 Ha/bohr^2


def _balanced_energy(R: float) -> float:
    """Full balanced-coupled LiH total energy (E_coupled) at bond length R."""
    spec = lih_spec(R=R, max_n=2)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=8000, L_max=4,
        screened_cross_center=False, verbose=False,
    )
    return float(coupled_fci_energy(ham, n_electrons=n_e, verbose=False)['E_coupled'])


def _balanced_energy_knobs(R: float, n_grid_vne: int, L_max: int) -> float:
    """Balanced LiH E_coupled at bond length R with explicit quadrature/angular knobs."""
    spec = lih_spec(R=R, max_n=2)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=n_grid_vne, L_max=L_max,
        screened_cross_center=False, verbose=False,
    )
    return float(coupled_fci_energy(ham, n_electrons=n_e, verbose=False)['E_coupled'])


@pytest.mark.slow
def test_paper19_well_shape_localized_to_orbital_basis():
    """The well-SHAPE residual (tilt, curvature) is bit-invariant to the cross-V_ne
    angular multipole order L_max and the radial quadrature n_grid_vne, and the
    electronic gradient at R_true is ~9% too weak.  => the shape defect lives in the
    max_n orbital basis, not in integral evaluation or angular truncation
    (debug/sprint_abc_connections_test_memo.md)."""
    R_TRUE = 3.015
    h = 0.05
    Rs = [R_TRUE - h, R_TRUE, R_TRUE + h]

    def tilt_curv(n_grid, L):
        E = np.array([_balanced_energy_knobs(r, n_grid, L) for r in Rs])
        tilt = (E[2] - E[0]) / (2 * h)               # central 1st derivative
        curv = (E[2] - 2 * E[1] + E[0]) / h ** 2     # central 2nd derivative
        return tilt, curv

    tilt_L2, curv_L2 = tilt_curv(8000, 2)
    tilt_L6, _ = tilt_curv(8000, 6)                  # angular refinement
    tilt_g2k, _ = tilt_curv(2000, 4)                 # coarse quadrature
    base_tilt, _ = tilt_curv(8000, 4)

    # (1) angular multipole order is bit-irrelevant (Gaunt-terminated at L=2)
    assert abs(tilt_L2 - tilt_L6) < 1e-9, f"L_max moved the tilt: {tilt_L2} vs {tilt_L6}"
    # (2) radial quadrature is bit-irrelevant (converged by n_grid=2000)
    assert abs(base_tilt - tilt_g2k) < 1e-9, f"n_grid moved the tilt: {base_tilt} vs {tilt_g2k}"
    # (3) the defect is real and outward (negative tilt ~ -0.030 Ha/bohr)
    assert -0.045 < base_tilt < -0.020, f"tilt={base_tilt:.4f}"
    # (4) electronic gradient 9% too weak: dE/dR = dV_NN/dR + d<el>/dR, V_NN=3/R exact,
    #     required d<el>/dR = +3/R^2 for zero net tilt at the true minimum
    req = 3.0 / R_TRUE ** 2
    elec_grad = base_tilt + req                       # = base_tilt - dV_NN/dR
    deficit = (req - elec_grad) / req
    assert 0.05 < deficit < 0.13, f"electronic-gradient deficit={deficit:.3f}"


def _omega_e(curv_ha_bohr2: float) -> float:
    """Harmonic frequency (cm^-1) from a force constant in Ha/bohr^2."""
    return np.sqrt(max(curv_ha_bohr2, 0.0) / MU_ME) / CM1_TO_HA


@pytest.mark.slow
def test_paper19_balanced_lih_well_too_stiff():
    # bracket the balanced minimum (~3.22 bohr) with a window wide enough that
    # the well rise (~few mHa) dominates the ~1 mHa grid noise at n_grid=8000
    R = np.array([3.02, 3.10, 3.18, 3.26, 3.34, 3.42])
    E = np.array([_balanced_energy(r) for r in R])

    # quartic fit; curvature and minimum are offset-immune (the R-independent
    # core-energy convention constant drops out of the derivatives)
    c = np.polyfit(R - 3.22, E, 4)
    p = np.poly1d(c)
    dp, ddp = p.deriv(1), p.deriv(2)
    roots = dp.r[np.isreal(dp.r)].real
    mins = [r for r in roots if ddp(r) > 0 and (R.min() - 3.22) < r < (R.max() - 3.22)]
    assert mins, "no interior PES minimum found (monotonic PES?) -- well-shape claim not testable"
    r_eq = min(mins, key=lambda r: p(r)) + 3.22
    curv = float(ddp(r_eq - 3.22))
    omega = _omega_e(curv)

    # (1) R_eq drifts outward to ~3.22 (energy is separately 0.20% at n_max=3)
    assert 3.18 < r_eq < 3.28, f"R_eq={r_eq:.3f}"
    # (2) the well is ~2x too stiff in the force constant
    assert 1.9 < curv / K_TRUE < 2.6, f"curv/k_true={curv / K_TRUE:.2f}"
    # (3) => omega_e ~ 2040 cm^-1, ~45% above experiment (the paper's number)
    assert 1900.0 < omega < 2200.0, f"omega_e={omega:.0f} cm^-1"
    assert 1.35 < omega / W_EXP_CM1 < 1.55, f"omega/exp={omega / W_EXP_CM1:.2f}"
