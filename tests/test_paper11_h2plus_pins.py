"""Paper 11 headline pins (owed guards, /qa group2 FULL 2026-11-03).

Claims (paper_11_prolate_spheroidal.tex, abstract-adjacent results paragraph
and Table tab:spectral_convergence):

  * spectral Laguerre radial solver, n_basis = 20, R = 2.0 bohr:
    E_total reproduces Peek's -0.6026342144949 Ha to 3.6e-14 Ha
    ("reference-limited").
  * the equilibrium bond length agrees with the exact R_eq = 1.99719 bohr.

The earlier tests (tests/test_paper11_spectral_radial.py) accepted 0.01 % in
energy and 1 % in R_eq, so a +1.21e-6 Ha offset planted in the solver stayed
green. These pins exclude that.

Expected values are EXTERNAL literature numbers (Peek 1965 for E(R=2);
Bates-Ledsham-Stewart / Peek for R_eq = 1.99719 bohr), never read from the
code under test. R_eq is located by a direct 1-D minimization of the total
energy (Brent), not by the module's own quadratic-fit helper, so the guard
does not share the fit's coarse-grid bias.
"""
import os
import sys

import pytest
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_spheroidal_lattice import ProlateSpheroidalLattice

PEEK_E_R2 = -0.6026342144949   # Ha, Peek (1965), quoted to 13 decimals
EXACT_R_EQ = 1.99719           # bohr, exact H2+ equilibrium (BLS 1953 / Peek)


def _e_total(R: float) -> float:
    lat = ProlateSpheroidalLattice(R=R, radial_method='spectral', n_basis=20)
    return lat.total_energy()


def test_paper11_energy_R2_reference_limited() -> None:
    """|E(R=2.0) - Peek| < 1e-12 Ha at n_basis = 20 (paper: 3.6e-14)."""
    E = _e_total(2.0)
    err = abs(E - PEEK_E_R2)
    assert err < 1e-12, f"E(R=2)={E!r}, |E - Peek| = {err:.3e} Ha >= 1e-12"


def test_paper11_req_matches_exact() -> None:
    """Direct minimizer of the n_basis=20 PES within 1e-5 relative of 1.99719.

    The owed-guard spec asked for 1e-4; the measured minimizer agrees to
    ~2e-6 relative (the reference itself is quoted to 2.5e-6), so 1e-5 is used:
    it still admits the reference rounding but fires on an R-dependent plant
    of slope 1.2e-5 Ha/bohr that 1e-4 would accept.
    """
    res = minimize_scalar(_e_total, bracket=(1.9, 2.0, 2.1), tol=1e-10)
    R_eq = res.x
    rel = abs(R_eq - EXACT_R_EQ) / EXACT_R_EQ
    assert 1.9 < R_eq < 2.1, f"minimizer left the bracket: R_eq={R_eq}"
    assert rel < 1e-5, f"R_eq={R_eq:.7f}, rel. dev. {rel:.2e} from 1.99719"
