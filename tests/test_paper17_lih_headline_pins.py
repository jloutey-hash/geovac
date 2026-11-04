"""Paper 17 LiH headline pins (owed guards, /qa group2 FULL 2026-11-03).

Claims (paper_17_composed_geometries.tex, abstract + Table tab:lih +
"l-dependent pseudopotential" paragraph):

  * composed LiH, ab initio PK (B from the <1/r^2>-weighted core radius),
    channel-blind, l_max = 2: R_eq = 3.22 bohr (6.9 % vs 3.015),
    omega_e = 1465 cm^-1 (4.2 % vs 1406);
  * l-dependent PK (delta_{l,0}) at l_max = 2: R_eq error 5.3 %
    -> R_eq = 3.176 bohr.

tests/test_ab_initio_pk_v2.py only asserts 2.5 <= R_eq <= 4.0, so swapping the
barrier-width formula (B_method 'inv2' -> 'median') stayed green. These pins
use windows set by the paper's printed precision: 6.9 % and 5.3 % pin R_eq to
+-0.0015 bohr (half a unit in the last printed digit of the percentage) and
1465 cm^-1 pins omega_e to +-1 cm^-1. (Tightened 2026-11-04 after an
independent review showed a 2 % change of the PK width passed the old
+-0.01 / +-10 windows.) Same R grid as tests/test_ab_initio_pk_v2.py.

Each PES scan is ~4.5 min, hence slow.
"""
import numpy as np
import pytest

from geovac.composed_diatomic import ComposedDiatomicSolver

PAPER_REQ_CHANNEL_BLIND = 3.222   # bohr  (paper: 3.22, 6.9 % above 3.015)
PAPER_OMEGA_CHANNEL_BLIND = 1465.0  # cm^-1
PAPER_REQ_L_DEPENDENT = 3.176     # bohr  (paper: 5.3 % above 3.015)


def _lih_r_grid() -> np.ndarray:
    return np.concatenate([
        np.linspace(2.0, 2.5, 3),
        np.linspace(2.7, 4.0, 10),
        np.linspace(4.5, 7.0, 5),
    ])


def _scan(mode: str) -> ComposedDiatomicSolver:
    s = ComposedDiatomicSolver.LiH_ab_initio(
        l_max=2, verbose=False, pk_channel_mode=mode,
    )
    s.solve_core()
    s.scan_pes(R_grid=_lih_r_grid(), n_Re=300)
    s.fit_spectroscopic_constants(fit_window=1.5)
    return s


@pytest.mark.slow
def test_paper17_lih_channel_blind_headline() -> None:
    """Channel-blind ab initio PK, l_max=2: R_eq 3.222 +- 0.01, omega 1465 +- 10."""
    s = _scan('channel_blind')
    assert s.pk_mode == 'ab_initio'
    R_eq = s.spectro['R_eq']
    omega = s.spectro['omega_e']
    pct = (R_eq / 3.015 - 1.0) * 100.0
    assert 6.85 <= pct < 6.95, f"R_eq = {R_eq:.4f} bohr = {pct:.3f} % (paper 6.9 %)"
    assert abs(omega - PAPER_OMEGA_CHANNEL_BLIND) < 1.0, (
        f"omega_e = {omega:.1f} cm^-1, paper 1465 (+-1)")


@pytest.mark.slow
def test_paper17_lih_l_dependent_req() -> None:
    """l-dependent PK, l_max=2: R_eq 3.176 +- 0.01 (5.3 % vs expt)."""
    s = _scan('l_dependent')
    R_eq = s.spectro['R_eq']
    pct = (R_eq / 3.015 - 1.0) * 100.0
    assert 5.25 <= pct < 5.35, f"R_eq = {R_eq:.4f} bohr = {pct:.3f} % (paper 5.3 %)"
