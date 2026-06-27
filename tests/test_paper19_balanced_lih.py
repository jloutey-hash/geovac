"""
Paper 19 (Balanced Coupled Composition) — permanent backing for the
balanced-LiH headline numbers.

Paper 19 §IV reports, for LiH at n_max=2 with the *balanced* coupled
Hamiltonian (all one-body + two-body + cross-center V_ne terms, no PK):

    E (n_max=2)  = -7.924 .. -7.929 Ha   (1.7-1.8% vs exact -8.071)
    R_eq         =  3.226 .. 3.227 bohr  (7.0% vs exact 3.015)
    D_e          =  0.037 Ha             (the only bound coupling scheme)

These numbers had no permanent regression backing — they reproduced only
under an *archived* driver (debug/archive/.../chemistry_solver_retest_lih*.py),
which is pruned by design (CLAUDE.md §9 clean-room rule).  This module
replaces that backing with a recompute-from-framework test.

The "TC convention" (the accounting the paper's energies use):
    E_min(TC) = E_min_fit - E_core_off
where
    E_min_fit  = quadratic-fit minimum of the balanced-coupled FCI energy
                 coupled_fci_energy(build_balanced_hamiltonian(
                     lih_spec(), nuclei=None, R=R))['E_coupled']
                 over a fine R-grid, and
    E_core_off = geovac.molecular_spec._FIRST_ROW_CORE_ENERGY[3] (= -7.2799).

The raw E_coupled (~-15.21 Ha) carries E_core + V_NN baked into the
Hamiltonian constant; subtracting E_core_off recovers the -8.071-comparable
total.  (This is NOT a fitted parameter — it is the He-like Li(2+) core
energy already used inside molecular_spec; the test reads it from the
module rather than hardcoding it.)

Verification type (CLAUDE.md §13.4a): numerical cross-check against the
exact LiH energy (-8.071 Ha) and experimental R_eq (3.015 bohr), with a
genuine variational-type sanity bound (E above exact).

The build is ~40 s per R-point, so the recompute is marked @slow.
The n_max=3 result (0.20% energy error) is larger and much slower
(Q=84, FCI dim grows); it is left as a deferred coverage item.
"""

import numpy as np
import pytest

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import lih_spec, _FIRST_ROW_CORE_ENERGY


# Reference values (NOT framework outputs — exact/experimental anchors).
E_EXACT_LIH = -8.071      # Ha, exact LiH total energy (Paper 19 Table)
R_EQ_EXPT = 3.015         # bohr, experimental LiH bond length


def _balanced_lih_energy(spec, R: float) -> float:
    """Balanced-coupled 4-electron FCI energy of LiH at bond length R.

    Mirrors exactly the paper's production path:
    coupled_fci_energy(build_balanced_hamiltonian(spec, nuclei=None, R=R)).
    """
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(spec, nuclei=None, R=float(R))
    fci = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
    return float(fci['E_coupled'])


def _parabolic_min(R: np.ndarray, E: np.ndarray):
    """Quadratic-fit minimum using the 3 points centred on the discrete min."""
    i = int(np.argmin(E))
    i = min(max(i, 1), len(R) - 2)  # keep an interior 3-window
    R3, E3 = R[i - 1:i + 2], E[i - 1:i + 2]
    a, b, c = np.polyfit(R3, E3, 2)
    R_eq = -b / (2.0 * a)
    E_min = a * R_eq * R_eq + b * R_eq + c
    return float(R_eq), float(E_min)


@pytest.fixture(scope="module")
def balanced_lih_pes():
    """Balanced-coupled LiH PES on a fine R-grid bracketing the minimum.

    Builds the spec once (n_max=2 default) and recomputes the balanced
    Hamiltonian + FCI at each R, exactly as the paper's driver did.
    """
    spec = lih_spec()  # LiH, max_n=2 (Paper 19 n_max=2 row)
    R_grid = np.array([3.0, 3.1, 3.2, 3.3, 3.4])
    E = np.array([_balanced_lih_energy(spec, R) for R in R_grid])
    R_eq_fit, E_min_fit = _parabolic_min(R_grid, E)
    E_core_off = _FIRST_ROW_CORE_ENERGY[3]
    return {
        'R_grid': R_grid,
        'E_coupled': E,
        'n_electrons': sum(b.n_electrons for b in spec.blocks),
        'R_eq_fit': R_eq_fit,
        'E_min_fit': E_min_fit,
        'E_core_off': E_core_off,
        'E_min_TC': E_min_fit - E_core_off,
    }


@pytest.mark.slow
class TestPaper19BalancedLiH:
    """Recompute the Paper 19 balanced-LiH n_max=2 headline numbers."""

    def test_four_electrons(self, balanced_lih_pes):
        """LiH default spec is the 4-electron system (2 core + 2 valence)."""
        assert balanced_lih_pes['n_electrons'] == 4

    def test_core_offset_is_module_value(self, balanced_lih_pes):
        """E_core_off is the He-like Li(2+) core energy from molecular_spec."""
        # Read from the module (not hardcoded); pinned to the documented value.
        assert balanced_lih_pes['E_core_off'] == pytest.approx(-7.2799, abs=1e-4)

    def test_pes_is_bound(self, balanced_lih_pes):
        """The balanced scheme binds: a minimum exists in the grid interior."""
        E = balanced_lih_pes['E_coupled']
        i_min = int(np.argmin(E))
        assert 0 < i_min < len(E) - 1, (
            f"min at edge (i={i_min}); PES not bracketed: {E}")

    def test_E_min_TC_matches_exact(self, balanced_lih_pes):
        """E_min(TC) ≈ -7.93 Ha: within 2.5% of exact -8.071, above exact."""
        E_TC = balanced_lih_pes['E_min_TC']
        err = abs(E_TC - E_EXACT_LIH) / abs(E_EXACT_LIH)
        print(f"\n  E_min(TC) = {E_TC:.4f} Ha  "
              f"(exact {E_EXACT_LIH}, err {err*100:.2f}%)")
        # Paper: 1.7-1.8% error. Variational-type: TC energy above exact.
        assert err < 0.025, f"E_min(TC) error {err*100:.2f}% exceeds 2.5%"
        assert E_TC > E_EXACT_LIH, "E_min(TC) below exact (non-physical)"
        # Pin the magnitude near the paper's -7.93 (not just <2.5%).
        assert -7.99 < E_TC < -7.87, f"E_min(TC)={E_TC:.4f} not near -7.93"

    def test_R_eq_seven_percent(self, balanced_lih_pes):
        """R_eq ≈ 3.23 bohr: 6-8% error vs experimental 3.015 (Paper 7.0%)."""
        R_eq = balanced_lih_pes['R_eq_fit']
        err = abs(R_eq - R_EQ_EXPT) / R_EQ_EXPT
        print(f"  R_eq = {R_eq:.4f} bohr  "
              f"(expt {R_EQ_EXPT}, err {err*100:.2f}%)")
        assert 0.06 <= err <= 0.08, (
            f"R_eq error {err*100:.2f}% outside paper's 6-8% band "
            f"(R_eq={R_eq:.4f})")
