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
