"""
Test: Origin of the matched scale kappa = -1/16
================================================

Precision audit of the claim that -1/16 is a "topological invariant."

This test answers the question: does -1/16 emerge from the topology alone,
or is it the scale factor that sends the graph's band edge to -1/2?

Method:
  1. Construct the graph Laplacian at several max_n values
  2. Find the extremal eigenvalue lambda_max at each resolution
  3. Compute the scale factor needed to match E_1 = -0.5 Ha exactly
  4. Plot convergence and fit the limiting value
  5. Report the honest answer

Author: Precision Audit Suite
Date: February 2026
"""

import numpy as np
import sys
import os
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
from typing import Tuple, List, Dict

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac import GeometricLattice


# =============================================================================
# Core Analysis Functions
# =============================================================================

def build_graph_laplacian(max_n: int) -> Tuple[np.ndarray, int]:
    """
    Build the pure graph Laplacian L = D - A for hydrogen at given max_n.

    Returns the Laplacian as a sparse matrix and the number of states.
    """
    lattice = GeometricLattice(max_n=max_n)
    adjacency = lattice.adjacency
    n_states = lattice.num_states

    degree = np.array(adjacency.sum(axis=1)).flatten()
    D = diags(degree, 0, shape=(n_states, n_states), format='csr')
    laplacian = D - adjacency

    return laplacian, n_states


def get_lowest_eigenvalues(max_n: int, n_eigs: int = 5) -> np.ndarray:
    """
    Get the lowest eigenvalues of the pure (unscaled) graph Laplacian.

    The graph Laplacian L = D - A is positive semi-definite.
    The lowest eigenvalue is 0 (constant mode on each connected component).
    We want the lowest non-trivial eigenvalue(s).
    """
    laplacian, n_states = build_graph_laplacian(max_n)
    k = min(n_eigs, n_states - 1)
    eigenvalues, _ = eigsh(laplacian, k=k, which='SA')
    return np.sort(eigenvalues)


def compute_required_scale_factor(max_n: int, target_energy: float = -0.5) -> Dict:
    """
    Compute the scale factor kappa needed so that kappa * lambda_min = target_energy.

    For the hydrogen ground state: target_energy = -0.5 Ha.
    """
    eigenvalues = get_lowest_eigenvalues(max_n, n_eigs=5)

    # The graph Laplacian eigenvalues are non-negative.
    # The Hamiltonian H = kappa * L, where kappa < 0 gives negative energies.
    # The ground state energy = kappa * lambda_max (largest eigenvalue of L
    # becomes the most negative energy when multiplied by negative kappa).
    #
    # Actually, we need to think carefully:
    # H = kappa * L where kappa = -1/16
    # Eigenvalues of H = kappa * eigenvalues of L
    # Since kappa < 0, the LARGEST eigenvalue of L gives the LOWEST energy.
    # So we need the largest eigenvalue of L, not the smallest.

    # Get largest eigenvalue
    laplacian, n_states = build_graph_laplacian(max_n)
    k = min(10, n_states - 1)
    largest_eigs, _ = eigsh(laplacian, k=k, which='LA')
    lambda_max = np.max(largest_eigs)

    # Required: kappa * lambda_max = target_energy
    # So kappa = target_energy / lambda_max
    kappa_required = target_energy / lambda_max

    return {
        'max_n': max_n,
        'n_states': n_states,
        'lambda_max': lambda_max,
        'kappa_required': kappa_required,
        'kappa_ratio_to_1_16': kappa_required / (-1/16),
        'smallest_eigs': eigenvalues[:5],
        'largest_eigs': np.sort(largest_eigs)[-5:],
    }


# =============================================================================
# Pytest Tests
# =============================================================================

class TestUniversalConstantOrigin:
    """
    Precision audit: what is -1/16 and where does it come from?
    """

    def test_graph_laplacian_is_dimensionless(self):
        """
        Verify the graph Laplacian has no intrinsic physical scale.

        The eigenvalues of L = D - A are pure numbers determined by
        graph connectivity alone. No physical constants appear.
        """
        for max_n in [5, 10, 15]:
            eigs = get_lowest_eigenvalues(max_n, n_eigs=3)
            # All eigenvalues should be real, non-negative (PSD property)
            assert np.all(eigs >= -1e-10), (
                f"Graph Laplacian should be PSD, got negative eigenvalue at max_n={max_n}"
            )

    def test_scale_factor_convergence(self):
        """
        The scale factor kappa that maps graph eigenvalues to -0.5 Ha
        converges to -1/16 = -0.0625 as max_n increases.

        This is THE key test: does -1/16 emerge, and from what?
        """
        results = []
        max_n_values = [5, 10, 15, 20, 25, 30]

        for max_n in max_n_values:
            result = compute_required_scale_factor(max_n)
            results.append(result)

        # The required kappa should converge toward -1/16
        kappas = [r['kappa_required'] for r in results]

        # At large max_n, kappa should be close to -0.0625
        assert abs(kappas[-1] - (-0.0625)) / 0.0625 < 0.10, (
            f"At max_n={max_n_values[-1]}, kappa={kappas[-1]:.6f}, "
            f"expected ~-0.0625 (error: {abs(kappas[-1] - (-0.0625)) / 0.0625 * 100:.1f}%)"
        )

        # Kappa should be monotonically converging (getting closer to -0.0625)
        errors = [abs(k - (-0.0625)) for k in kappas]
        for i in range(1, len(errors)):
            # Allow small non-monotonicity due to finite-size effects
            assert errors[i] < errors[0] * 1.5, (
                f"Convergence not improving: error at max_n={max_n_values[i]} "
                f"({errors[i]:.6f}) worse than at max_n={max_n_values[0]} ({errors[0]:.6f})"
            )

    def test_kappa_is_not_independently_derived(self):
        """
        Verify that -1/16 is determined BY the requirement to match physics,
        not derived from topology alone.

        The graph Laplacian eigenvalues are dimensionless numbers.
        Without knowing the target energy (-0.5 Ha), there is no way to
        determine kappa from the graph alone.  kappa is a matched scale on
        the graph's extremal eigenvalue (a 2D-grid band edge); the graph
        spectrum is NOT the hydrogen spectrum (the edge Laplacian does not
        converge to Delta_{S^3}, Paper 7).

        Honest framing: kappa = -1/16 matches lambda_max to -1/2 by
        construction; the bridge from the dimensionless labels to
        dimensionful energies is Fock's energy shell p0^2 = -2E, not kappa.
        kappa's coincidence with the Fock-Jacobian 1/16 is an Observation.
        """
        # Compute kappa at two different max_n values
        r1 = compute_required_scale_factor(max_n=10)
        r2 = compute_required_scale_factor(max_n=20)

        # Key observation: kappa changes with max_n because the graph
        # eigenvalues change. As max_n -> infinity (lambda_max -> 8, the band
        # edge), kappa -> -1/16.
        # kappa is a matched scale whose required value tends to -1/16 as lambda_max -> 8.
        kappa_10 = r1['kappa_required']
        kappa_20 = r2['kappa_required']

        # They should be different (proving kappa depends on resolution)
        # but converging toward the same limit
        assert kappa_10 != kappa_20, (
            "Kappa should vary with max_n (it's resolution-dependent)"
        )

        # Both should be in the neighborhood of -0.0625
        for max_n, kappa in [(10, kappa_10), (20, kappa_20)]:
            assert abs(kappa - (-0.0625)) / 0.0625 < 0.15, (
                f"Kappa at max_n={max_n} should be near -0.0625, got {kappa:.6f}"
            )

    def test_pure_laplacian_does_not_resolve_shells(self):
        """
        IMPORTANT AUDIT FINDING: The pure graph Laplacian (D-A) alone
        does NOT resolve the Rydberg shell structure (E_n = -1/2n^2).

        The eigenvalues of kappa*(D-A) all cluster near -0.5 Ha,
        without clear separation between n=1, n=2, n=3 shells.

        This finding is consistent with the honest framing: the graph
        spectrum is not the hydrogen spectrum, and kappa only matches the
        band edge to -1/2.
        """
        from geovac import AtomicSolver

        solver = AtomicSolver(max_n=20, Z=1, kinetic_scale=-1/16,
                              hamiltonian="graph")
        energies, _ = solver.compute_ground_state(n_states=10)
        energies = np.sort(energies)

        # All eigenvalues should be negative
        assert np.all(energies < 0), "All eigenvalues should be negative"

        # The ground state should be near -0.5 Ha
        assert abs(energies[0] - (-0.5)) / 0.5 < 0.05, (
            f"Ground state {energies[0]:.6f} should be near -0.5 Ha"
        )

        # Key finding: eigenvalues are CLUSTERED, not well-separated
        # The spread (max-min) should be small relative to the ground state
        spread = abs(energies[-1] - energies[0])
        relative_spread = spread / abs(energies[0])
        assert relative_spread < 0.20, (
            f"Eigenvalue spread {relative_spread:.4f} should be small, "
            f"confirming the pure Laplacian clusters eigenvalues"
        )

    def test_z_scaling_preserves_kappa(self):
        """
        kappa * Z^2 * lambda_max matches -Z^2/2 by construction for every Z.

        This checks only the legacy graph model's Z^2 scaling of its
        matched ground-state scale (kappa * Z^2 * lambda_max); kappa is a
        matched scale, not a unique bridge constant (the bridge to
        dimensionful energies is Fock's energy shell).
        """
        from geovac import AtomicSolver

        max_n = 20
        exact_energies = {1: -0.5, 2: -2.0, 3: -4.5}

        for Z, E_exact in exact_energies.items():
            solver = AtomicSolver(max_n=max_n, Z=Z, kinetic_scale=-1/16,
                              hamiltonian="graph")
            energies, _ = solver.compute_ground_state(n_states=1)
            E_computed = energies[0]
            error_pct = abs((E_computed - E_exact) / E_exact) * 100

            assert error_pct < 3.0, (
                f"Z={Z}: E={E_computed:.6f}, exact={E_exact}, error={error_pct:.2f}%"
            )


# =============================================================================
# Standalone Report Generator
# =============================================================================

def generate_full_report() -> str:
    """
    Generate the complete precision audit report for the -1/16 constant.
    Returns the report as a string.
    """
    lines = []
    lines.append("=" * 74)
    lines.append("PRECISION AUDIT: Origin of the matched scale kappa = -1/16")
    lines.append("=" * 74)
    lines.append("")

    # === Part 1: Convergence Analysis ===
    lines.append("PART 1: Scale Factor Convergence")
    lines.append("-" * 74)
    lines.append(f"{'max_n':>6s} {'N_states':>10s} {'lambda_max':>12s} "
                 f"{'kappa_req':>12s} {'ratio':>8s} {'error%':>8s}")
    lines.append("-" * 74)

    max_n_values = [5, 10, 15, 20, 25, 30]
    results = []

    for max_n in max_n_values:
        r = compute_required_scale_factor(max_n)
        results.append(r)
        error_pct = abs(r['kappa_required'] - (-0.0625)) / 0.0625 * 100
        lines.append(
            f"{r['max_n']:>6d} {r['n_states']:>10d} {r['lambda_max']:>12.6f} "
            f"{r['kappa_required']:>12.8f} {r['kappa_ratio_to_1_16']:>8.4f} "
            f"{error_pct:>8.2f}"
        )

    lines.append("")

    # === Part 2: Eigenvalue Ratios (topological content) ===
    lines.append("PART 2: Eigenvalue Ratios (Pure Topology)")
    lines.append("-" * 74)
    lines.append("These ratios are determined by graph structure alone,")
    lines.append("independent of any physical scale factor.")
    lines.append("")

    for max_n in [10, 20, 30]:
        laplacian, n_states = build_graph_laplacian(max_n)
        k = min(15, n_states - 1)
        largest_eigs, _ = eigsh(laplacian, k=k, which='LA')
        largest_eigs = np.sort(largest_eigs)[::-1]

        lines.append(f"  max_n = {max_n} (N = {n_states} states):")
        lines.append(f"    Top 5 eigenvalues of L: "
                     f"{', '.join(f'{e:.4f}' for e in largest_eigs[:5])}")

        if len(largest_eigs) >= 2:
            lines.append(f"    lambda_2 / lambda_1 = {largest_eigs[1]/largest_eigs[0]:.6f}")
        lines.append("")

    # === Part 3: Z-scaling verification ===
    lines.append("PART 3: Z-Scaling Verification")
    lines.append("-" * 74)
    lines.append("kappa*Z^2*lambda_max matches -Z^2/2 by construction for every Z.")
    lines.append("")

    from geovac import AtomicSolver
    exact_energies = {1: -0.5, 2: -2.0, 3: -4.5, 4: -8.0}
    for Z, E_exact in exact_energies.items():
        solver = AtomicSolver(max_n=20, Z=Z, kinetic_scale=-1/16,
                              hamiltonian="graph")
        energies, _ = solver.compute_ground_state(n_states=1)
        E_computed = energies[0]
        error_pct = abs((E_computed - E_exact) / E_exact) * 100
        lines.append(f"  Z={Z}: E_computed={E_computed:>12.6f} Ha, "
                     f"E_exact={E_exact:>8.2f} Ha, error={error_pct:.4f}%")

    lines.append("")

    # === Part 4: Honest Assessment ===
    lines.append("=" * 74)
    lines.append("HONEST ASSESSMENT")
    lines.append("=" * 74)
    lines.append("")
    lines.append("Q: Does -1/16 emerge from the topology alone?")
    lines.append("")
    lines.append("A: No. Here is the precise chain of reasoning:")
    lines.append("")
    lines.append("  1. The graph Laplacian L = D - A is a dimensionless matrix")
    lines.append("     whose eigenvalues are determined purely by graph connectivity.")
    lines.append("     This is genuinely topological.")
    lines.append("")
    lines.append("  2. The graph spectrum is NOT the hydrogen spectrum: D - A has")
    lines.append("     no l-changing edge, so it splits into max_n disconnected")
    lines.append("     l-blocks with spectrum in [0, 8] (Paper 7, WALL). The")
    lines.append("     n^2 - 1 ladder is the S3 Casimir on the labels, not a graph limit.")
    lines.append("")
    lines.append("  3. kappa = -1/16 sends the band edge (8) to E_1 = -1/2 Ha;")
    lines.append("     it does not map the graph spectrum onto the Rydberg ladder.")
    lines.append("     This match requires EXTERNAL INPUT: the physical energy")
    lines.append("     scale (the Hartree).")
    lines.append("")
    lines.append("  4. -1/16 maps the extremal eigenvalue (band edge 8) to -1/2;")
    lines.append("     its equality with the Fock-Jacobian 1/16 is an Observation")
    lines.append("     (Papers 0/7), a matching condition, not a derivation.")
    lines.append("")
    lines.append("RECOMMENDED FRAMING:")
    lines.append("  '-1/16 is the matched scale that sends the graph's band")
    lines.append("  edge to the hydrogen ground state (an Observation). The")
    lines.append("  bridge to dimensionful energies is Fock's energy shell")
    lines.append("  p0^2 = -2E, not the graph spectrum.'")
    lines.append("")

    # Convergence data for potential plotting
    lines.append("=" * 74)
    lines.append("CONVERGENCE DATA (for plotting)")
    lines.append("=" * 74)
    lines.append("max_n, n_states, lambda_max, kappa_required, error_pct")
    for r in results:
        error_pct = abs(r['kappa_required'] - (-0.0625)) / 0.0625 * 100
        lines.append(f"{r['max_n']}, {r['n_states']}, {r['lambda_max']:.6f}, "
                     f"{r['kappa_required']:.8f}, {error_pct:.4f}")

    return "\n".join(lines)


if __name__ == "__main__":
    report = generate_full_report()
    print(report)

    # Save to debug/data/
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'debug', 'data')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'constant_origin_audit.txt')
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"\nReport saved to: {output_path}")

    # Generate convergence plot if matplotlib available
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        max_n_values = [5, 10, 15, 20, 25, 30]
        kappas = []
        for max_n in max_n_values:
            r = compute_required_scale_factor(max_n)
            kappas.append(r['kappa_required'])

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Plot 1: kappa vs max_n
        ax1.plot(max_n_values, kappas, 'bo-', markersize=8, linewidth=2)
        ax1.axhline(y=-0.0625, color='r', linestyle='--', linewidth=1.5,
                     label='$\\kappa = -1/16$')
        ax1.set_xlabel('$n_{max}$ (lattice truncation)', fontsize=12)
        ax1.set_ylabel('Required $\\kappa$', fontsize=12)
        ax1.set_title('Scale Factor Convergence', fontsize=14)
        ax1.legend(fontsize=11)
        ax1.grid(True, alpha=0.3)

        # Plot 2: % error vs max_n
        errors = [abs(k - (-0.0625)) / 0.0625 * 100 for k in kappas]
        ax2.semilogy(max_n_values, errors, 'rs-', markersize=8, linewidth=2)
        ax2.set_xlabel('$n_{max}$ (lattice truncation)', fontsize=12)
        ax2.set_ylabel('Error vs $-1/16$ (%)', fontsize=12)
        ax2.set_title('Convergence Rate', fontsize=14)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plot_dir = os.path.join(os.path.dirname(__file__), '..', 'debug', 'plots')
        os.makedirs(plot_dir, exist_ok=True)
        plot_path = os.path.join(plot_dir, 'constant_origin_convergence.png')
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {plot_path}")

    except ImportError:
        print("matplotlib not available, skipping plot generation")
