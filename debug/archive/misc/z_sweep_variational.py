"""
Z-sweep investigation: Where does the graph-native CI variational bound break?

For He-like ions (2 electrons, nuclear charge Z), the graph-native CI uses:
  h1 = exact diagonal (-Z^2/2n^2) + graph off-diagonal (kappa=-1/16)
  V_ee = Slater integrals at orbital exponent k=Z

This is a MODEL Hamiltonian. At Z=2 (He), E_CI > E_exact (variational).
At Z=1 (H-), E_CI < E_exact (non-variational, over-binding by 21%).
This script maps the crossover as a function of Z.

Track DI investigation, April 2026.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import time
import numpy as np
from geovac.casimir_ci import two_electron_integral
from geovac.lattice import GeometricLattice


def exact_energy_he_like(Z: float) -> float:
    """Approximate exact ground-state energy for 2-electron ion with charge Z.

    Uses known values where available, otherwise 1/Z perturbation expansion.
    For Z < Z_critical ~ 0.9112, no bound state exists (returns NaN).
    """
    known = {
        1.0: -0.52775,    # H-
        2.0: -2.90372,    # He
        3.0: -7.27991,    # Li+
        5.0: -22.0310,    # B3+
        10.0: -93.907,    # Ne8+
    }
    if Z in known:
        return known[Z]

    Z_crit = 0.9112
    if Z < Z_crit:
        return float('nan')

    # 1/Z expansion: E(Z) ~ -Z^2 + (5/8)Z - 0.15767 + O(1/Z)
    return -Z**2 + (5.0/8.0) * Z - 0.15767


# Cache lattice data since it's Z-independent
_lattice_cache = {}

def _get_lattice_data(n_max: int):
    """Get cached lattice adjacency and orbital list."""
    if n_max not in _lattice_cache:
        lattice = GeometricLattice(max_n=n_max)
        orbitals = list(lattice.states)
        A = lattice.adjacency
        A_dense = A.toarray() if hasattr(A, 'toarray') else np.array(A)
        # Precompute singlet configs
        n_spatial = len(orbitals)
        configs = []
        for i in range(n_spatial):
            for j in range(i, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == 0:
                    configs.append((i, j))
        _lattice_cache[n_max] = (orbitals, A_dense, configs, n_spatial)
    return _lattice_cache[n_max]


# Cache V_ee integrals (Z-dependent only through k_orb linear scaling)
_vee_cache = {}

def _get_vee_matrix(n_max: int, orbitals, configs):
    """Get cached V_ee FCI matrix at k_orb=1.0 (scale by Z later)."""
    if n_max not in _vee_cache:
        n_configs = len(configs)
        H_vee = np.zeros((n_configs, n_configs))
        parity = 1.0  # singlet

        for I in range(n_configs):
            i, j = configs[I]
            for J in range(I, n_configs):
                p, q = configs[J]

                bra_perms = [(i, j, 1.0)]
                if i != j:
                    bra_perms.append((j, i, parity))
                ket_perms = [(p, q, 1.0)]
                if p != q:
                    ket_perms.append((q, p, parity))

                N_IJ = np.sqrt(float(len(bra_perms)))
                N_PQ = np.sqrt(float(len(ket_perms)))

                me = 0.0
                for a, b, s_bra in bra_perms:
                    for c, d, s_ket in ket_perms:
                        sign = s_bra * s_ket
                        na, la, ma = orbitals[a]
                        nb, lb, mb = orbitals[b]
                        nc, lc, mc = orbitals[c]
                        nd, ld, md = orbitals[d]
                        me += sign * two_electron_integral(
                            na, la, ma, nb, lb, mb,
                            nc, lc, mc, nd, ld, md, 1.0)

                me /= (N_IJ * N_PQ)
                H_vee[I, J] = me
                H_vee[J, I] = me

        _vee_cache[n_max] = H_vee
    return _vee_cache[n_max]


def _get_h1_fci_matrix(n_max: int, orbitals, A_dense, configs, n_spatial):
    """Get FCI matrix of h1 at Z=1 (diagonal=-1/2n^2, off-diag=kappa).

    Since diagonal scales as Z^2 and off-diag is Z-independent,
    we precompute two matrices: h1_diag_fci (Z^2-scaled) and h1_offdiag_fci (fixed).
    """
    cache_key = ('h1_parts', n_max)
    if cache_key not in _lattice_cache:
        n_configs = len(configs)
        parity = 1.0

        # Build h1_diag and h1_offdiag as separate n_spatial x n_spatial matrices
        h1_diag = np.zeros((n_spatial, n_spatial))
        h1_offdiag = np.zeros((n_spatial, n_spatial))

        for i, (n, l, m) in enumerate(orbitals):
            h1_diag[i, i] = -1.0 / (2.0 * n * n)  # at Z=1

        kappa = -1.0 / 16.0
        for i in range(n_spatial):
            for j in range(n_spatial):
                if i != j and abs(A_dense[i, j]) > 1e-15:
                    h1_offdiag[i, j] = kappa * (-A_dense[i, j])

        # Now lift these to FCI space
        H_diag = np.zeros((n_configs, n_configs))
        H_offdiag = np.zeros((n_configs, n_configs))

        for I in range(n_configs):
            i, j = configs[I]
            for J in range(I, n_configs):
                p, q = configs[J]

                bra_perms = [(i, j, 1.0)]
                if i != j:
                    bra_perms.append((j, i, parity))
                ket_perms = [(p, q, 1.0)]
                if p != q:
                    ket_perms.append((q, p, parity))

                N_IJ = np.sqrt(float(len(bra_perms)))
                N_PQ = np.sqrt(float(len(ket_perms)))

                me_d = 0.0
                me_o = 0.0
                for a, b, s_bra in bra_perms:
                    for c, d, s_ket in ket_perms:
                        sign = s_bra * s_ket
                        if b == d:
                            me_d += sign * h1_diag[a, c]
                            me_o += sign * h1_offdiag[a, c]
                        if a == c:
                            me_d += sign * h1_diag[b, d]
                            me_o += sign * h1_offdiag[b, d]

                me_d /= (N_IJ * N_PQ)
                me_o /= (N_IJ * N_PQ)
                H_diag[I, J] = me_d
                H_diag[J, I] = me_d
                H_offdiag[I, J] = me_o
                H_offdiag[J, I] = me_o

        _lattice_cache[cache_key] = (H_diag, H_offdiag)
    return _lattice_cache[cache_key]


def build_fci_fast(Z: float, n_max: int):
    """Build FCI Hamiltonian efficiently using cached components.

    H_FCI(Z) = Z^2 * H_diag + H_offdiag + Z * H_vee

    Returns H_FCI, configs, H_diag_fci, H_offdiag_fci, H_vee_fci
    """
    orbitals, A_dense, configs, n_spatial = _get_lattice_data(n_max)
    H_vee = _get_vee_matrix(n_max, orbitals, configs)
    H_diag, H_offdiag = _get_h1_fci_matrix(n_max, orbitals, A_dense, configs, n_spatial)

    # H(Z) = Z^2 * H_diag + H_offdiag + Z * H_vee
    H = Z * Z * H_diag + H_offdiag + Z * H_vee
    return H, configs, H_diag, H_offdiag, H_vee


def main():
    Z_values = [0.5, 0.75, 0.90, 0.95, 1.0, 1.1, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0]
    n_max_range = list(range(1, 8))  # 1 through 7

    print("=" * 100)
    print("Z-SWEEP: Graph-Native CI Variational Bound Investigation")
    print("=" * 100)

    # Precompute lattice and integrals for each n_max
    print("\nPrecomputing lattice structures and integrals...")
    for n_max in n_max_range:
        t0 = time.time()
        orbitals, A_dense, configs, n_spatial = _get_lattice_data(n_max)
        _ = _get_vee_matrix(n_max, orbitals, configs)
        _ = _get_h1_fci_matrix(n_max, orbitals, A_dense, configs, n_spatial)
        t1 = time.time()
        print(f"  n_max={n_max}: {len(configs)} configs, precomputed in {t1-t0:.1f}s")

    # Main sweep
    results = {}

    for Z in Z_values:
        E_exact = exact_energy_he_like(Z)
        results[Z] = {'E_exact': E_exact, 'E_CI': {}, 'n_configs': {}}

        print(f"\n{'='*80}")
        print(f"Z = {Z:.2f}   E_exact = {E_exact:.6f} Ha")
        print(f"{'='*80}")
        print(f"{'n_max':>6} {'n_cfg':>6} {'E_CI':>14} {'E_CI-E_ex':>14} {'error%':>10} {'bound?':>8}")
        print("-" * 70)

        for n_max in n_max_range:
            H, configs, _, _, _ = build_fci_fast(Z, n_max)
            n_configs = len(configs)

            if n_configs == 0:
                print(f"{n_max:>6} {'--':>6}")
                continue

            evals, evecs = np.linalg.eigh(H)
            E_CI = evals[0]

            results[Z]['E_CI'][n_max] = E_CI
            results[Z]['n_configs'][n_max] = n_configs

            if np.isnan(E_exact):
                delta = E_CI  # compare to 0 or just show raw
                err_pct = float('nan')
                bound_str = "N/A"
            else:
                delta = E_CI - E_exact
                err_pct = 100.0 * (E_CI - E_exact) / abs(E_exact)
                bound_str = "YES" if delta > 0 else "NO"

            print(f"{n_max:>6} {n_configs:>6} {E_CI:>14.6f} {delta:>14.6f} {err_pct:>10.4f} {bound_str:>8}")

    # Summary table
    print("\n\n" + "=" * 100)
    print("SUMMARY: Signed Error at n_max=7")
    print("=" * 100)
    print(f"{'Z':>8} {'E_CI(7)':>14} {'E_exact':>14} {'Delta':>14} {'err%':>10} {'Variational?':>14}")
    print("-" * 80)

    for Z in Z_values:
        E_exact = results[Z]['E_exact']
        E_CI_7 = results[Z]['E_CI'].get(7, None)
        if E_CI_7 is None:
            print(f"{Z:>8.2f} {'N/A':>14}")
            continue

        if np.isnan(E_exact):
            print(f"{Z:>8.2f} {E_CI_7:>14.6f} {'unbound':>14} {'N/A':>14} {'N/A':>10} {'N/A':>14}")
        else:
            delta = E_CI_7 - E_exact
            err_pct = 100.0 * delta / abs(E_exact)
            bound_str = "YES" if delta > 0 else "NO"
            print(f"{Z:>8.2f} {E_CI_7:>14.6f} {E_exact:>14.6f} {delta:>14.6f} {err_pct:>10.4f} {bound_str:>14}")

    # Find crossover from main sweep
    print("\n\n" + "=" * 100)
    print("CROSSOVER ANALYSIS (from main sweep)")
    print("=" * 100)

    for n_max in [3, 5, 7]:
        print(f"\n--- n_max = {n_max} ---")
        prev_Z = None
        prev_delta = None
        for Z in Z_values:
            E_exact = results[Z]['E_exact']
            E_CI = results[Z]['E_CI'].get(n_max, None)
            if E_CI is None or np.isnan(E_exact):
                prev_Z = None
                prev_delta = None
                continue
            delta = E_CI - E_exact
            if prev_delta is not None and prev_delta * delta < 0:
                Z_cross = prev_Z + (Z - prev_Z) * (-prev_delta) / (delta - prev_delta)
                print(f"  Sign change between Z={prev_Z:.2f} (delta={prev_delta:.6f}) "
                      f"and Z={Z:.2f} (delta={delta:.6f})")
                print(f"  Linear interpolation: Z_crossover ~ {Z_cross:.4f}")
            prev_Z = Z
            prev_delta = delta

    # Refined bisection at n_max=5 (fast) and n_max=7
    for n_max_bisect in [5, 7]:
        print(f"\n\n{'='*100}")
        print(f"REFINED CROSSOVER SEARCH (n_max={n_max_bisect})")
        print(f"{'='*100}")

        Z_lo, Z_hi = None, None
        for Z in Z_values:
            E_exact = exact_energy_he_like(Z)
            E_CI = results[Z]['E_CI'].get(n_max_bisect, None)
            if E_CI is None or np.isnan(E_exact):
                continue
            delta = E_CI - E_exact
            if delta < 0:
                Z_lo = Z
            elif delta > 0 and Z_lo is not None and Z_hi is None:
                Z_hi = Z

        if Z_lo is not None and Z_hi is not None:
            print(f"Bracket: Z in [{Z_lo:.2f}, {Z_hi:.2f}]")
            for _ in range(30):
                Z_mid = (Z_lo + Z_hi) / 2.0
                E_exact_mid = exact_energy_he_like(Z_mid)
                H, _, _, _, _ = build_fci_fast(Z_mid, n_max_bisect)
                evals, _ = np.linalg.eigh(H)
                E_CI_mid = evals[0]
                delta_mid = E_CI_mid - E_exact_mid

                if delta_mid < 0:
                    Z_lo = Z_mid
                else:
                    Z_hi = Z_mid

                if abs(Z_hi - Z_lo) < 1e-6:
                    break

            Z_cross = (Z_lo + Z_hi) / 2.0
            E_exact_c = exact_energy_he_like(Z_cross)
            H_c, _, _, _, _ = build_fci_fast(Z_cross, n_max_bisect)
            evals_c, _ = np.linalg.eigh(H_c)
            E_CI_c = evals_c[0]
            print(f"Converged: Z_crossover = {Z_cross:.6f}")
            print(f"  E_CI = {E_CI_c:.6f}, E_exact = {E_exact_c:.6f}, delta = {E_CI_c - E_exact_c:.8f}")
        else:
            print("Could not bracket crossover.")

    # Energy decomposition at n_max=7
    print("\n\n" + "=" * 100)
    print("ENERGY DECOMPOSITION at n_max=7")
    print("=" * 100)
    print(f"{'Z':>8} {'<h1_diag>':>14} {'<h1_offdiag>':>14} {'<V_ee>':>14} "
          f"{'E_CI':>14} {'<Vee>/|<h1>|':>14} {'1/(8Z^2)':>12}")
    print("-" * 110)

    for Z in Z_values:
        E_CI = results[Z]['E_CI'].get(7, None)
        if E_CI is None:
            continue

        H, configs, H_diag, H_offdiag, H_vee = build_fci_fast(Z, 7)
        evals, evecs = np.linalg.eigh(H)
        evec = evecs[:, 0]

        exp_diag = Z * Z * (evec @ H_diag @ evec)
        exp_offdiag = evec @ H_offdiag @ evec
        exp_vee = Z * (evec @ H_vee @ evec)
        exp_h1 = exp_diag + exp_offdiag
        E_check = exp_diag + exp_offdiag + exp_vee

        ratio_vee_h1 = abs(exp_vee) / abs(exp_h1) if abs(exp_h1) > 1e-15 else float('nan')
        offdiag_scale = 1.0 / (8.0 * Z * Z)

        E_exact = results[Z]['E_exact']
        print(f"{Z:>8.2f} {exp_diag:>14.6f} {exp_offdiag:>14.6f} {exp_vee:>14.6f} "
              f"{E_check:>14.6f} {ratio_vee_h1:>14.6f} {offdiag_scale:>12.6f}")

    # Off-diagonal analysis
    print("\n\n" + "=" * 100)
    print("GRAPH OFF-DIAGONAL ANALYSIS")
    print("=" * 100)
    print("Graph h1 off-diagonal (kappa=-1/16) is Z-INDEPENDENT.")
    print("Diagonal scales as Z^2. Relative importance: |kappa|/(Z^2/2) = 1/(8Z^2)")
    print()
    print(f"{'Z':>8} {'diag(1s)':>12} {'|off-diag|':>12} {'ratio':>12} {'comment':>30}")
    print("-" * 80)
    for Z in Z_values:
        diag_1s = -Z * Z / 2.0
        off_diag = 1.0 / 16.0
        ratio = off_diag / abs(diag_1s)
        if ratio > 0.5:
            comment = "OFF-DIAG DOMINATES"
        elif ratio > 0.1:
            comment = "comparable"
        elif ratio > 0.01:
            comment = "perturbative"
        else:
            comment = "negligible"
        print(f"{Z:>8.2f} {diag_1s:>12.6f} {off_diag:>12.6f} {ratio:>12.6f} {comment:>30}")

    # Convergence behavior
    print("\n\n" + "=" * 100)
    print("CONVERGENCE: E_CI(n_max) - E_CI(n_max-1) for key Z values")
    print("=" * 100)
    print("Negative = energy lowering with more basis functions")
    print(f"{'Z':>8}", end='')
    for n_max in range(2, 8):
        print(f" {'dE('+str(n_max)+')':>12}", end='')
    print()
    print("-" * 85)
    for Z in [0.95, 1.0, 1.25, 1.5, 2.0, 5.0]:
        print(f"{Z:>8.2f}", end='')
        for n_max in range(2, 8):
            E_prev = results[Z]['E_CI'].get(n_max-1, None)
            E_curr = results[Z]['E_CI'].get(n_max, None)
            if E_prev is not None and E_curr is not None:
                dE = E_curr - E_prev
                print(f" {dE:>12.6f}", end='')
            else:
                print(f" {'N/A':>12}", end='')
        print()

    # Physical interpretation
    print("\n\n" + "=" * 100)
    print("PHYSICAL INTERPRETATION")
    print("=" * 100)
    print("""
The graph-native CI Hamiltonian H = Z^2 * H_diag + H_offdiag + Z * H_vee has three terms:

1. H_diag (scales as Z^2): exact hydrogenic diagonal energies -Z^2/(2n^2)
2. H_offdiag (Z-independent): graph Laplacian inter-shell coupling, kappa = -1/16
3. H_vee (scales as Z): Slater electron-electron repulsion at k=Z

The model is NOT the exact projected Hamiltonian, so there is no guaranteed variational
bound. The crossover occurs because:

- At large Z: H_diag dominates. The graph off-diagonal is a tiny perturbation (~1/(8Z^2)
  of the diagonal). V_ee/Z^2 -> 0 so correlation is small. The model closely
  approximates the exact Hamiltonian. E_CI > E_exact (variational-like).

- At Z ~ 1: H_offdiag is comparable to H_diag. The kappa=-1/16 couplings create
  effective kinetic energy channels that don't exist in the exact Hamiltonian. These
  extra channels allow the CI to find states that would not exist in the exact basis,
  leading to artificial lowering of the energy below exact. E_CI < E_exact.

- At Z < 1: H_offdiag DOMINATES. The graph topology imposes a structure that is
  qualitatively wrong for weakly-bound systems. The artificial inter-shell coupling
  creates binding energy that doesn't physically exist.

Key insight: The graph off-diagonal kappa=-1/16 is the graph Laplacian's representation
of kinetic energy. In the exact Hamiltonian, kinetic energy is diagonal in the hydrogen
eigenbasis (it equals -E_n = Z^2/(2n^2)). The graph adds EXTRA off-diagonal kinetic
couplings that are artifacts of the discrete topology. These artifacts are negligible
at high Z but catastrophic at low Z.
""")

    print("\nDone.")


if __name__ == '__main__':
    main()
