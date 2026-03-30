"""
Tests for the Berry phase module (geovac/berry_phase.py).

Validates plaquette Berry phases on the geometric lattice using topological
weights w = 1/(n1*n2). The plaquette loop traverses:

    |n,l,m> -> T+ -> |n+1,l,m> -> L+ -> |n+1,l,m+1> -> T- -> |n,l,m+1> -> L- -> |n,l,m>

The log-plaquette phase is:
    ln(w_T+) + ln(w_L+) - ln(w_T-) - ln(w_L-)
    = ln(1/(n(n+1))) + ln(1/(n+1)^2) - ln(1/(n(n+1))) - ln(1/n^2)
    = -2 * ln((n+1)/n)
"""

import numpy as np
import pytest

from geovac._archive.dead_ends.berry_phase import (
    compute_plaquette_phase,
    compute_berry_phase_at_n,
    compute_phase_convergence,
    fit_power_law,
)
from geovac.lattice import GeometricLattice


# ---------------------------------------------------------------------------
# Reference values
# ---------------------------------------------------------------------------
# Plaquette phase at n=2: -2 * ln(3/2)
PHASE_N2_EXACT = -2.0 * np.log(3.0 / 2.0)  # ~ -0.8109

# At n=2: l=0 -> 0 plaquettes, l=1 -> m=-1,0 -> 2 plaquettes. Total = 2.
PLAQUETTES_AT_N2 = 2

# Total Berry phase at n=2: 2 * (-2*ln(3/2))
TOTAL_PHASE_N2_EXACT = PLAQUETTES_AT_N2 * PHASE_N2_EXACT

TOL = 1e-12


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def lattice_topo() -> GeometricLattice:
    """Lattice with max_n=4, topological weights."""
    return GeometricLattice(max_n=4, nuclear_charge=1, topological_weights=True)


@pytest.fixture
def lattice_binary() -> GeometricLattice:
    """Lattice with max_n=4, binary (uniform) weights."""
    return GeometricLattice(max_n=4, nuclear_charge=1, topological_weights=False)


# ---------------------------------------------------------------------------
# TestBerryPhaseComputation
# ---------------------------------------------------------------------------
class TestBerryPhaseComputation:
    """Tests for individual plaquette and shell-level Berry phase values."""

    def test_plaquette_phase_analytical(self, lattice_topo: GeometricLattice) -> None:
        """Plaquette phase at (n=2, l=1, m=0) equals -2*ln(3/2) exactly.

        The plaquette loop:
            (2,1,0) -T+-> (3,1,0) -L+-> (3,1,1) -T--> (2,1,1) -L--> (2,1,0)

        Weights:
            T+ : 1/(2*3) = 1/6
            L+ : 1/(3*3) = 1/9
            T- : 1/(2*3) = 1/6
            L- : 1/(2*2) = 1/4

        Phase = ln(1/6) + ln(1/9) - ln(1/6) - ln(1/4)
              = ln(1/9) - ln(1/4) = ln(4/9) = -2*ln(3/2)
        """
        phase = compute_plaquette_phase(lattice_topo, n=2, l=1, m=0)
        np.testing.assert_allclose(
            phase, PHASE_N2_EXACT, atol=TOL,
            err_msg=f"Plaquette phase at (2,1,0): got {phase}, expected {PHASE_N2_EXACT}"
        )

    def test_plaquette_phase_binary_zero(self, lattice_binary: GeometricLattice) -> None:
        """With uniform (binary) weights, all plaquette phases vanish.

        If every edge weight is 1.0, then ln(1)+ln(1)-ln(1)-ln(1) = 0.
        """
        phase = compute_plaquette_phase(lattice_binary, n=2, l=1, m=0)
        np.testing.assert_allclose(
            phase, 0.0, atol=TOL,
            err_msg="Binary (uniform) weights should give zero plaquette phase"
        )

    def test_berry_phase_at_n_count(self, lattice_topo: GeometricLattice) -> None:
        """Total Berry phase at n=2 equals 2 plaquettes * (-2*ln(3/2)).

        At n=2, valid plaquettes are:
            l=1, m=-1 -> plaquette (2,1,-1)
            l=1, m=0  -> plaquette (2,1,0)
        Total = 2 * phase_per_plaquette.
        """
        total = compute_berry_phase_at_n(lattice_topo, n=2)
        np.testing.assert_allclose(
            total, TOTAL_PHASE_N2_EXACT, atol=TOL,
            err_msg=f"Total phase at n=2: got {total}, expected {TOTAL_PHASE_N2_EXACT}"
        )

    def test_invalid_n_raises(self, lattice_topo: GeometricLattice) -> None:
        """n >= max_n should raise ValueError (no n+1 shell available)."""
        with pytest.raises(ValueError, match="n"):
            compute_berry_phase_at_n(lattice_topo, n=4)
        with pytest.raises(ValueError, match="n"):
            compute_berry_phase_at_n(lattice_topo, n=5)

    def test_plaquette_missing_states_returns_zero(
        self, lattice_topo: GeometricLattice
    ) -> None:
        """Plaquette at (n=2, l=0, m=0) returns 0: m+1=1 > l=0, state doesn't exist."""
        phase = compute_plaquette_phase(lattice_topo, n=2, l=0, m=0)
        assert phase == 0.0, "Plaquette with missing corner states should return 0.0"


# ---------------------------------------------------------------------------
# TestExponentConvergence
# ---------------------------------------------------------------------------
class TestExponentConvergence:
    """Tests for power-law scaling of Berry phase with shell number."""

    def test_power_law_fit_quality(self) -> None:
        """Power-law fit to theta(n) vs n should have R^2 > 0.95."""
        max_n = 30
        lattice = GeometricLattice(
            max_n=max_n, nuclear_charge=1, topological_weights=True
        )
        n_vals = np.arange(2, max_n)  # skip n=1 (zero phase)
        phases = np.array([
            compute_berry_phase_at_n(lattice, n) for n in n_vals
        ])

        abs_phases = np.abs(phases)
        mask = abs_phases > 1e-15
        assert np.sum(mask) >= 5, "Need at least 5 nonzero phase values"

        k, A = fit_power_law(n_vals[mask], abs_phases[mask])

        # Compute R^2
        log_n = np.log(n_vals[mask].astype(float))
        log_phase = np.log(abs_phases[mask])
        log_fit = np.log(A) - k * log_n
        ss_res = np.sum((log_phase - log_fit) ** 2)
        ss_tot = np.sum((log_phase - np.mean(log_phase)) ** 2)
        r_squared = 1.0 - ss_res / ss_tot

        assert r_squared > 0.95, (
            f"Power-law fit R^2 = {r_squared:.4f}, expected > 0.95. "
            f"Exponent k = {k:.4f}, amplitude A = {A:.6f}"
        )

    def test_exponent_reasonable(self) -> None:
        """Per-plaquette phase decays; total phase grows (negative k).

        The total Berry phase at shell n sums n*(n-1) plaquettes, each
        of magnitude ~2/n. So total ~ 2*(n-1) grows linearly, giving
        negative k (~ -1.3). The per-plaquette phase decays as ~n^(-1).
        We check both behaviors are consistent.
        """
        max_n = 20
        lattice = GeometricLattice(
            max_n=max_n, nuclear_charge=1, topological_weights=True
        )
        n_vals = np.arange(2, max_n)
        phases = np.array([
            compute_berry_phase_at_n(lattice, n) for n in n_vals
        ])

        # Total phases should grow in magnitude
        abs_phases = np.abs(phases)
        assert abs_phases[-1] > abs_phases[0], (
            "Total Berry phase should grow with n (more plaquettes)"
        )

        # Per-plaquette phase should decay
        n_plaquettes = np.array([n * (n - 1) for n in n_vals], dtype=float)
        per_plaquette = abs_phases / np.maximum(n_plaquettes, 1.0)
        mask = per_plaquette > 1e-15
        k_per, _ = fit_power_law(n_vals[mask], per_plaquette[mask])
        assert 0.5 <= k_per <= 3.0, (
            f"Per-plaquette exponent k = {k_per:.4f} outside [0.5, 3.0]"
        )
