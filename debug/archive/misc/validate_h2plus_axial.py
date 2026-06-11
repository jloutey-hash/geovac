"""
Minimal H2+ test: 1D axial graph Laplacian between two protons.

Tests whether a simple nearest-neighbor chain along the internuclear axis
produces a PES with a minimum near R_eq = 2 bohr.

H2+ exact: R_eq = 1.997 bohr, E = -0.6026 Ha, D_e = 0.1026 Ha

Approach:
  - Vertices on the z-axis (protons at z=0 and z=R)
  - Edges: nearest-neighbor with kinetic weight
  - V = -1/|z - R_A| - 1/|z - R_B| (softened Coulomb)
  - H = T + V, solve for ground state, add V_NN = 1/R
"""
import numpy as np
from typing import List, Tuple
import os

# Exact H2+ values
EXACT_R_EQ = 1.997   # bohr
EXACT_E = -0.6026    # Ha
EXACT_DE = 0.1026    # Ha (relative to H + p)
E_H = -0.5           # Ha (hydrogen atom)


def build_h2plus_hamiltonian(
    R: float,
    N: int = 5,
    padding: float = 0.0,
    softening: float = 0.1,
    kinetic: str = 'fd',
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Build H2+ Hamiltonian on a 1D axial grid.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    N : int
        Number of grid points.
    padding : float
        Extra space beyond nuclei on each side.
    softening : float
        Coulomb softening parameter (bohr).
    kinetic : str
        'fd' = finite-difference -1/(2h²), 'topo' = topological -1/16.

    Returns
    -------
    H : ndarray (N, N)
    z : ndarray (N,)
    V_NN : float
    """
    if padding > 0:
        z = np.linspace(-padding, R + padding, N)
    else:
        z = np.linspace(0, R, N)
    h = z[1] - z[0]

    # Kinetic energy: graph Laplacian of a chain
    # L_{ij} = degree(i) delta_{ij} - A_{ij}
    # For a chain: L = tridiag(-1, 2, -1) at interior, (1,-1) at ends
    L = np.zeros((N, N))
    for i in range(N):
        neighbors = 0
        if i > 0:
            L[i, i - 1] = -1.0
            neighbors += 1
        if i < N - 1:
            L[i, i + 1] = -1.0
            neighbors += 1
        L[i, i] = float(neighbors)

    if kinetic == 'fd':
        # Finite-difference: T = -1/(2h²) * (-L) = 1/(2h²) * L
        scale = 1.0 / (2.0 * h**2)
    elif kinetic == 'topo':
        # Topological: T = (-1/16) * L  (GeoVac universal constant)
        scale = -1.0 / 16.0
    else:
        raise ValueError(f"Unknown kinetic type: {kinetic}")

    T = scale * L

    # Nuclear potential
    V = np.zeros(N)
    for k in range(N):
        r_A = np.sqrt(z[k]**2 + softening**2)
        r_B = np.sqrt((z[k] - R)**2 + softening**2)
        V[k] = -1.0 / r_A - 1.0 / r_B

    H = T + np.diag(V)
    V_NN = 1.0 / R

    return H, z, V_NN


def scan_pes(
    R_values: np.ndarray,
    N: int = 5,
    padding: float = 0.0,
    softening: float = 0.1,
    kinetic: str = 'fd',
) -> List[Tuple[float, float, float]]:
    """Scan PES: returns list of (R, E_elec, E_total)."""
    results = []
    for R in R_values:
        H, z, V_NN = build_h2plus_hamiltonian(R, N, padding, softening, kinetic)
        evals = np.linalg.eigvalsh(H)
        E_elec = evals[0]
        E_total = E_elec + V_NN
        results.append((R, E_elec, E_total))
    return results


def fit_pes(results: List[Tuple[float, float, float]]) -> dict:
    """Fit PES to extract R_eq, E_min, D_e, k."""
    Rs = np.array([r[0] for r in results])
    Es = np.array([r[2] for r in results])

    idx_min = np.argmin(Es)

    # Check if minimum is at boundary
    if idx_min == 0 or idx_min == len(Rs) - 1:
        return {
            'R_eq': Rs[idx_min], 'E_min': Es[idx_min],
            'D_e': E_H - Es[idx_min], 'k': 0.0,
            'bound': Es[idx_min] < E_H, 'boundary': True,
        }

    # Parabolic fit around minimum
    lo = max(0, idx_min - 2)
    hi = min(len(Rs), idx_min + 3)
    R_fit = Rs[lo:hi]
    E_fit = Es[lo:hi]
    coeffs = np.polyfit(R_fit, E_fit, 2)
    a, b, c = coeffs

    R_eq = -b / (2 * a)
    E_min = np.polyval(coeffs, R_eq)
    k = 2 * a
    D_e = E_H - E_min

    return {
        'R_eq': R_eq, 'E_min': E_min, 'D_e': D_e, 'k': k,
        'bound': E_min < E_H, 'boundary': False,
    }


def main():
    print("=" * 72)
    print("  H2+ Axial Lattice Test")
    print("  Exact: R_eq=1.997, E=-0.6026 Ha, D_e=0.1026 Ha")
    print("=" * 72)

    R_values = np.arange(0.5, 8.1, 0.25)

    # ===================================================================
    # Test 1: User's exact proposal — 5 vertices on [0, R], no padding
    # ===================================================================
    print("\n--- Test 1: 5 vertices on [0, R], FD kinetic, softening=0.1 ---")
    results = scan_pes(R_values, N=5, padding=0.0, softening=0.1, kinetic='fd')
    fit = fit_pes(results)
    print(f"  Bound: {fit['bound']}, Boundary: {fit['boundary']}")
    print(f"  R_eq = {fit['R_eq']:.3f} (exact 1.997, err {abs(fit['R_eq']-EXACT_R_EQ)/EXACT_R_EQ*100:.1f}%)")
    print(f"  E_min = {fit['E_min']:.4f} (exact -0.6026)")
    print(f"  D_e = {fit['D_e']:.4f} (exact 0.1026)")
    print(f"  k = {fit['k']:.4f}")

    # Print PES
    print(f"\n  {'R':>6s}  {'E_elec':>10s}  {'E_total':>10s}")
    for R, E_el, E_tot in results:
        marker = " <-- min" if abs(E_tot - fit['E_min']) < 0.001 else ""
        print(f"  {R:6.2f}  {E_el:10.6f}  {E_tot:10.6f}{marker}")

    # ===================================================================
    # Test 2: More vertices, with padding (convergence study)
    # ===================================================================
    print("\n--- Test 2: Convergence with N (padding=4.0, softening=0.1) ---")
    print(f"  {'N':>5s}  {'R_eq':>7s}  {'E_min':>10s}  {'D_e':>8s}  {'R_err%':>7s}  {'E_err%':>7s}")
    for N in [5, 11, 21, 41, 81, 161]:
        results_N = scan_pes(R_values, N=N, padding=4.0, softening=0.1, kinetic='fd')
        fit_N = fit_pes(results_N)
        if fit_N['bound'] and not fit_N['boundary']:
            R_err = abs(fit_N['R_eq'] - EXACT_R_EQ) / EXACT_R_EQ * 100
            E_err = abs(fit_N['E_min'] - EXACT_E) / abs(EXACT_E) * 100
            print(f"  {N:5d}  {fit_N['R_eq']:7.3f}  {fit_N['E_min']:10.6f}  {fit_N['D_e']:8.4f}  {R_err:7.1f}  {E_err:7.1f}")
        else:
            print(f"  {N:5d}  {'UNBOUND' if not fit_N['bound'] else 'BNDRY':>7s}  {fit_N['E_min']:10.6f}  {fit_N['D_e']:8.4f}")

    # ===================================================================
    # Test 3: Softening parameter study
    # ===================================================================
    print("\n--- Test 3: Softening parameter (N=81, padding=4.0) ---")
    print(f"  {'eps':>6s}  {'R_eq':>7s}  {'E_min':>10s}  {'D_e':>8s}  {'R_err%':>7s}  {'E_err%':>7s}")
    for eps in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
        results_e = scan_pes(R_values, N=81, padding=4.0, softening=eps, kinetic='fd')
        fit_e = fit_pes(results_e)
        if fit_e['bound'] and not fit_e['boundary']:
            R_err = abs(fit_e['R_eq'] - EXACT_R_EQ) / EXACT_R_EQ * 100
            E_err = abs(fit_e['E_min'] - EXACT_E) / abs(EXACT_E) * 100
            print(f"  {eps:6.3f}  {fit_e['R_eq']:7.3f}  {fit_e['E_min']:10.6f}  {fit_e['D_e']:8.4f}  {R_err:7.1f}  {E_err:7.1f}")
        else:
            print(f"  {eps:6.3f}  {'UNBOUND' if not fit_e['bound'] else 'BNDRY':>7s}  {fit_e['E_min']:10.6f}")

    # ===================================================================
    # Test 4: Topological kinetic (-1/16) vs finite-difference
    # ===================================================================
    print("\n--- Test 4: Topological kinetic (-1/16) vs FD (N=41, padding=4.0) ---")
    for kin_label, kin in [('fd', 'fd'), ('topo', 'topo')]:
        results_k = scan_pes(R_values, N=41, padding=4.0, softening=0.1, kinetic=kin)
        fit_k = fit_pes(results_k)
        if fit_k['bound'] and not fit_k['boundary']:
            R_err = abs(fit_k['R_eq'] - EXACT_R_EQ) / EXACT_R_EQ * 100
            print(f"  {kin_label:>5s}: R_eq={fit_k['R_eq']:.3f} ({R_err:.1f}%), E={fit_k['E_min']:.6f}, D_e={fit_k['D_e']:.4f}")
        else:
            print(f"  {kin_label:>5s}: {'UNBOUND' if not fit_k['bound'] else 'BNDRY'}, E={fit_k['E_min']:.6f}")

    # ===================================================================
    # Test 5: Best converged PES for comparison
    # ===================================================================
    print("\n--- Test 5: Best PES (N=161, padding=6.0, eps=0.01) ---")
    R_fine = np.arange(0.5, 8.01, 0.1)
    results_best = scan_pes(R_fine, N=161, padding=6.0, softening=0.01, kinetic='fd')
    fit_best = fit_pes(results_best)
    if fit_best['bound'] and not fit_best['boundary']:
        R_err = abs(fit_best['R_eq'] - EXACT_R_EQ) / EXACT_R_EQ * 100
        E_err = abs(fit_best['E_min'] - EXACT_E) / abs(EXACT_E) * 100
        print(f"  R_eq = {fit_best['R_eq']:.4f} (exact 1.997, err {R_err:.2f}%)")
        print(f"  E_min = {fit_best['E_min']:.6f} (exact -0.6026, err {E_err:.2f}%)")
        print(f"  D_e = {fit_best['D_e']:.6f} (exact 0.1026)")
        print(f"  k = {fit_best['k']:.4f}")
    else:
        print(f"  {'UNBOUND' if not fit_best['bound'] else 'BNDRY'}")

    # Save data
    outdir = os.path.join(os.path.dirname(__file__), 'data')
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, 'h2plus_axial.txt')
    with open(outfile, 'w') as f:
        f.write("# H2+ Axial Lattice PES (N=161, padding=6.0, eps=0.01)\n")
        f.write(f"# R_eq={fit_best['R_eq']:.4f} E_min={fit_best['E_min']:.6f} "
                f"D_e={fit_best['D_e']:.6f}\n")
        f.write(f"# {'R':>8s} {'E_elec':>12s} {'E_total':>12s}\n")
        for R, E_el, E_tot in results_best:
            f.write(f"  {R:8.3f} {E_el:12.6f} {E_tot:12.6f}\n")
    print(f"\n  Data saved to {outfile}")


if __name__ == '__main__':
    main()
