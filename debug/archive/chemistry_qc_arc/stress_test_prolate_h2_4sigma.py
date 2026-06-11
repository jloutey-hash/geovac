"""
Phase 3: H2 4σ-Orbital CI on Prolate Spheroidal Lattice
========================================================
Extend the minimal-basis (1σ_g, 1σ_u) CI to include 2σ_g and 2σ_u orbitals
from the H2+ solver. This gives a 6×6 symmetry-adapted CI in the ¹Σ_g⁺ sector.

The 6 CSFs (configuration state functions):
  1. |1σ_g²⟩         closed-shell
  2. |1σ_u²⟩         closed-shell
  3. |2σ_g²⟩         closed-shell
  4. |2σ_u²⟩         closed-shell
  5. ¹(1σ_g, 2σ_g)   open-shell singlet (gerade × gerade = gerade)
  6. ¹(1σ_u, 2σ_u)   open-shell singlet (ungerade × ungerade = gerade)

All orbitals have m=0, so the existing azimuthal-averaging V_ee code works
unchanged.

Date: 2026-03-13
"""

import warnings
warnings.filterwarnings('ignore')

import sys
import os
# Fix Windows console encoding for unicode characters (sigma, etc.)
if sys.stdout.encoding and sys.stdout.encoding.lower().startswith('cp'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
from scipy.linalg import eigh
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_spheroidal_lattice import (
    ProlateSpheroidalLattice,
    fit_spectroscopic_constants,
)

# Reuse the orbital/integral infrastructure from Phase 2
from debug.stress_test_prolate_h2 import (
    get_orbital_on_grid,
    compute_vee_integral,
)


# ============================================================
# Integral cache — avoid recomputing symmetric integrals
# ============================================================
class IntegralCache:
    """Cache for two-electron integrals (ab|cd) with 8-fold symmetry.

    Symmetries exploited (real orbitals):
      (ab|cd) = (ba|cd) = (ab|dc) = (ba|dc)    [particle 1,2 internal]
      (ab|cd) = (cd|ab)                          [particle exchange]
    """

    def __init__(self, orbitals: dict):
        self._orbs = orbitals
        self._cache = {}
        self._n_computed = 0

    def _canonical_key(self, a: str, b: str, c: str, d: str) -> tuple:
        """Return canonical key exploiting (ab|cd) = (cd|ab) and particle symmetry."""
        # Our compute_vee_integral(a,b,c,d) computes ∫ a(1)c(1)/r12 × b(2)d(2)
        # = ⟨ab|cd⟩_phys.
        # Symmetries of ⟨ab|cd⟩: ⟨ab|cd⟩ = ⟨ba|dc⟩ = ⟨cd|ab⟩ = ⟨dc|ba⟩
        # Also for real orbitals: ⟨ab|cd⟩ = ⟨ab|cd⟩* (trivially real)
        keys = [
            (a, b, c, d),
            (b, a, d, c),
            (c, d, a, b),
            (d, c, b, a),
        ]
        return min(keys)

    def get(self, a: str, b: str, c: str, d: str) -> float:
        """Get ⟨ab|cd⟩_phys = their(a,b,c,d), using cache."""
        key = self._canonical_key(a, b, c, d)
        if key not in self._cache:
            p, q, r, s = key
            val = compute_vee_integral(
                self._orbs[p], self._orbs[q], self._orbs[r], self._orbs[s]
            )
            self._cache[key] = val
            self._n_computed += 1
        return self._cache[key]


# ============================================================
# Build the 6×6 CI Hamiltonian
# ============================================================
def build_ci_matrix_6x6(
    eps: dict,
    integrals: IntegralCache,
    orb_labels: list,
    verbose: bool = True,
) -> np.ndarray:
    """Build the 6×6 ¹Σ_g⁺ CI matrix.

    CSFs:
      0: |g²⟩   1: |u²⟩   2: |G²⟩   3: |U²⟩
      4: ¹(g,G)  5: ¹(u,U)

    where g=1σ_g, u=1σ_u, G=2σ_g, U=2σ_u.

    Parameters
    ----------
    eps : dict
        One-electron energies {label: E_elec}.
    integrals : IntegralCache
        Two-electron integral cache.
    orb_labels : list
        ['g', 'u', 'G', 'U'] orbital labels.
    """
    g, u, G, U = orb_labels
    I = integrals.get  # shorthand: I(a,b,c,d) = ⟨ab|cd⟩_phys

    H = np.zeros((6, 6))

    # --- Diagonal elements ---

    # Closed-shell: H_ii = 2*eps_i + J_ii  where J_ii = ⟨ii|ii⟩
    H[0, 0] = 2*eps[g] + I(g, g, g, g)
    H[1, 1] = 2*eps[u] + I(u, u, u, u)
    H[2, 2] = 2*eps[G] + I(G, G, G, G)
    H[3, 3] = 2*eps[U] + I(U, U, U, U)

    # Open-shell singlet ¹(a,b): H = eps_a + eps_b + J_ab + K_ab
    # J_ab = ⟨ab|ab⟩, K_ab = ⟨ab|ba⟩
    H[4, 4] = eps[g] + eps[G] + I(g, G, g, G) + I(g, G, G, g)
    H[5, 5] = eps[u] + eps[U] + I(u, U, u, U) + I(u, U, U, u)

    # --- Off-diagonal: closed ↔ closed (double excitation) ---
    # ⟨i²|H|j²⟩ = ⟨ii|jj⟩  (for i≠j)
    H[0, 1] = H[1, 0] = I(g, g, u, u)
    H[0, 2] = H[2, 0] = I(g, g, G, G)
    H[0, 3] = H[3, 0] = I(g, g, U, U)
    H[1, 2] = H[2, 1] = I(u, u, G, G)
    H[1, 3] = H[3, 1] = I(u, u, U, U)
    H[2, 3] = H[3, 2] = I(G, G, U, U)

    # --- Off-diagonal: closed ↔ open-shell singlet ---
    # ⟨i²|H|¹(a,b)⟩:
    #   If i=a: √2 * ⟨ii|ib⟩
    #   If i=b: √2 * ⟨ii|ai⟩
    #   Otherwise: 0 (>2 orbital diff)

    # CSF4 = ¹(g,G):  couples to |g²⟩ (i=g=a) and |G²⟩ (i=G=b)
    H[0, 4] = H[4, 0] = np.sqrt(2) * I(g, g, g, G)   # i=g, a=g, b=G
    H[2, 4] = H[4, 2] = np.sqrt(2) * I(G, G, g, G)   # i=G, a=g, b=G → √2⟨GG|gG⟩
    # |u²⟩ and |U²⟩ don't couple to ¹(g,G): differ in >2 orbitals
    # H[1,4] = H[3,4] = 0

    # CSF5 = ¹(u,U):  couples to |u²⟩ (i=u=a) and |U²⟩ (i=U=b)
    H[1, 5] = H[5, 1] = np.sqrt(2) * I(u, u, u, U)   # i=u, a=u, b=U
    H[3, 5] = H[5, 3] = np.sqrt(2) * I(U, U, u, U)   # i=U, a=u, b=U → √2⟨UU|uU⟩
    # H[0,5] = H[2,5] = 0

    # --- Off-diagonal: open ↔ open (double excitation) ---
    # ¹(g,G) ↔ ¹(u,U): all 4 orbitals differ → double excitation
    # For singlet CSFs ¹(a,b) and ¹(c,d):
    # ⟨¹(a,b)|H|¹(c,d)⟩ = ⟨ac|bd⟩ + ⟨ad|bc⟩  (for the singlet coupling)
    # Wait — need to be careful. Let me derive this.
    #
    # ¹(a,b) = (1/√2)(|a↑b↓⟩ - |a↓b↑⟩)
    # ¹(c,d) = (1/√2)(|c↑d↓⟩ - |c↓d↑⟩)
    #
    # ⟨¹(a,b)|V|¹(c,d)⟩ = (1/2)[
    #   ⟨a↑b↓|V|c↑d↓⟩ - ⟨a↑b↓|V|c↓d↑⟩
    #   - ⟨a↓b↑|V|c↑d↓⟩ + ⟨a↓b↑|V|c↓d↑⟩
    # ]
    #
    # Each det-det matrix element via Slater-Condon (double excitation):
    # ⟨a↑b↓|V|c↑d↓⟩: a↑→c↑, b↓→d↓
    #   = ⟨ab|cd⟩×δ(↑↑)δ(↓↓) - ⟨ab|dc⟩×δ(↑↓)δ(↓↑) = ⟨ab|cd⟩
    #
    # ⟨a↑b↓|V|c↓d↑⟩: a↑→d↑, b↓→c↓
    #   = ⟨ab|dc⟩×δ(↑↑)δ(↓↓) - ⟨ab|cd⟩×δ(↑↓)δ(↓↑) = ⟨ab|dc⟩
    #
    # ⟨a↓b↑|V|c↑d↓⟩: a↓→d↓, b↑→c↑
    #   = ⟨ab|dc⟩×δ(↓↓)δ(↑↑) - ⟨ab|cd⟩×δ(↓↑)δ(↑↓) = ⟨ab|dc⟩
    #
    # ⟨a↓b↑|V|c↓d↑⟩: a↓→c↓, b↑→d↑
    #   = ⟨ab|cd⟩×δ(↓↓)δ(↑↑) - ⟨ab|dc⟩×δ(↓↑)δ(↑↓) = ⟨ab|cd⟩
    #
    # Total: (1/2)[⟨ab|cd⟩ - ⟨ab|dc⟩ - ⟨ab|dc⟩ + ⟨ab|cd⟩]
    #      = ⟨ab|cd⟩ - ⟨ab|dc⟩

    H[4, 5] = H[5, 4] = I(g, G, u, U) - I(g, G, U, u)

    if verbose:
        csf_names = ['|g²⟩', '|u²⟩', '|G²⟩', '|U²⟩', '¹(g,G)', '¹(u,U)']
        print("\n  6×6 CI matrix (electronic, no V_NN):")
        print(f"  {'':8s}", end='')
        for name in csf_names:
            print(f"{name:>10s}", end='')
        print()
        for i in range(6):
            print(f"  {csf_names[i]:8s}", end='')
            for j in range(6):
                print(f"{H[i,j]:10.4f}", end='')
            print()

    return H


# ============================================================
# 4σ CI driver
# ============================================================
def h2_4sigma_ci(
    R: float,
    N_xi_solve: int = 5000,
    N_grid: int = 60,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> dict:
    """Build 4σ-orbital (6×6) CI for H2.

    Generates 1σ_g, 1σ_u, 2σ_g, 2σ_u from the H2+ solver, computes all
    two-electron integrals via azimuthal averaging, builds and diagonalizes
    the 6×6 ¹Σ_g⁺ CI matrix.
    """
    t0 = time.time()

    if verbose:
        print(f"\n  H2 4σ CI at R={R:.3f} bohr")
        print(f"  Grid: {N_grid}×{N_grid}, N_xi_solve={N_xi_solve}")

    # --- Step 1: Generate all 4 orbitals ---
    orb_info = [
        ('g', 0, '1σ_g'),
        ('u', 1, '1σ_u'),
        ('G', 2, '2σ_g'),
        ('U', 3, '2σ_u'),
    ]

    orbitals = {}
    eps = {}

    for label, n_ang, name in orb_info:
        if verbose:
            print(f"  Solving for {name} (n_angular={n_ang})...")
        try:
            orb = get_orbital_on_grid(
                R=R, n_angular=n_ang, N_xi_solve=N_xi_solve,
                N_xi_grid=N_grid, N_eta_grid=N_grid,
                xi_max_grid=xi_max_grid, m=0,
            )
            orbitals[label] = orb
            eps[label] = orb['E_elec']
            if verbose:
                print(f"    E_elec({name}) = {orb['E_elec']:.6f} Ha")
        except Exception as e:
            if verbose:
                print(f"    FAILED: {e}")
            return {'R': R, 'E_total': np.nan, 'error': str(e)}

    t_orb = time.time() - t0
    if verbose:
        print(f"  Orbital generation: {t_orb:.1f}s")

    # --- Step 2: Compute all two-electron integrals ---
    if verbose:
        print("  Computing V_ee integrals...")

    t1 = time.time()
    cache = IntegralCache(orbitals)
    labels = ['g', 'u', 'G', 'U']

    # Pre-compute all needed integrals (the cache handles deduplication)
    # Diagonal J_ii
    for a in labels:
        cache.get(a, a, a, a)

    # Off-diagonal ⟨ii|jj⟩ for closed-shell coupling
    for i, a in enumerate(labels):
        for b in labels[i+1:]:
            cache.get(a, a, b, b)

    # J_ab and K_ab for open-shell singlets
    for a, b in [('g', 'G'), ('u', 'U')]:
        cache.get(a, b, a, b)  # J
        cache.get(a, b, b, a)  # K

    # Closed ↔ open-shell: ⟨ii|ib⟩ type
    for a, b in [('g', 'G'), ('u', 'U')]:
        cache.get(a, a, a, b)  # |a²⟩ ↔ ¹(a,b)
        cache.get(b, b, a, b)  # |b²⟩ ↔ ¹(a,b)

    # Open ↔ open: ⟨gG|uU⟩ and ⟨gG|Uu⟩
    cache.get('g', 'G', 'u', 'U')
    cache.get('g', 'G', 'U', 'u')

    t_vee = time.time() - t1
    if verbose:
        print(f"    {cache._n_computed} unique integrals computed in {t_vee:.1f}s")

    # --- Step 3: Build CI matrix ---
    H_CI = build_ci_matrix_6x6(eps, cache, labels, verbose=verbose)

    # --- Step 4: Diagonalize ---
    evals, evecs = eigh(H_CI)
    E_CI_elec = evals[0]
    V_NN = 1.0 / R
    E_total = E_CI_elec + V_NN

    # Also extract the 2×2 result for comparison
    H_2x2 = H_CI[:2, :2]
    evals_2x2, _ = eigh(H_2x2)
    E_2x2 = evals_2x2[0] + V_NN

    # Hartree-Fock (no correlation)
    E_HF = 2*eps['g'] + cache.get('g', 'g', 'g', 'g') + V_NN

    dt = time.time() - t0

    if verbose:
        csf_names = ['|g²⟩', '|u²⟩', '|G²⟩', '|U²⟩', '¹(g,G)', '¹(u,U)']
        c = evecs[:, 0]
        print(f"\n  CI eigenvalues: {evals}")
        print(f"\n  Ground state decomposition:")
        for i, name in enumerate(csf_names):
            print(f"    {name:10s}: c={c[i]:+.4f}  |c|²={c[i]**2:.4f}")

        print(f"\n  Energies:")
        print(f"    E_HF      = {E_HF:.6f} Ha")
        print(f"    E_2×2_CI  = {E_2x2:.6f} Ha")
        print(f"    E_6×6_CI  = {E_total:.6f} Ha")
        print(f"    E_corr(2×2) = {E_2x2 - E_HF:.6f} Ha")
        print(f"    E_corr(6×6) = {E_total - E_HF:.6f} Ha")
        print(f"    Improvement = {E_2x2 - E_total:.6f} Ha")
        print(f"  Total time: {dt:.1f}s")

    return {
        'R': R,
        'E_total': E_total,
        'E_2x2': E_2x2,
        'E_HF': E_HF,
        'E_CI_elec': E_CI_elec,
        'V_NN': V_NN,
        'eps': dict(eps),
        'evals': evals,
        'evecs': evecs[:, 0],
        'n_integrals': cache._n_computed,
        'time': dt,
    }


# ============================================================
# PES scan
# ============================================================
def scan_h2_4sigma_pes(
    R_values: np.ndarray,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> dict:
    """Scan H2 PES with 4σ CI."""
    results = []

    for R in R_values:
        try:
            res = h2_4sigma_ci(
                R=R, N_xi_solve=N_xi_solve, N_grid=N_grid,
                xi_max_grid=xi_max_grid, verbose=verbose,
            )
            results.append(res)
        except Exception as e:
            print(f"  R={R:.3f}: FAILED - {e}")
            results.append({
                'R': R, 'E_total': np.nan, 'E_2x2': np.nan,
                'E_HF': np.nan, 'time': 0, 'error': str(e),
            })

    R_arr = np.array([r['R'] for r in results])
    E_arr = np.array([r['E_total'] for r in results])
    E_2x2_arr = np.array([r.get('E_2x2', np.nan) for r in results])
    E_HF_arr = np.array([r.get('E_HF', np.nan) for r in results])

    return {
        'R': R_arr,
        'E_total': E_arr,
        'E_2x2': E_2x2_arr,
        'E_HF': E_HF_arr,
        'results': results,
    }


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("  H2 4σ-ORBITAL CI ON PROLATE SPHEROIDAL LATTICE")
    print("  Orbitals: 1σ_g, 1σ_u, 2σ_g, 2σ_u  |  CI: 6×6 ¹Σ_g⁺")
    print("  Date: 2026-03-13")
    print("=" * 70)

    # --- Reference values ---
    E_H = -0.5  # H atom energy
    E_H2_exact = -1.1745  # Exact H2 total energy at R_eq
    R_eq_exact = 1.401    # bohr
    D_e_exact = 0.1745    # Ha

    # --- Test 1: Single point at R=1.4 (near exact R_eq) ---
    print("\n" + "="*60)
    print("  TEST 1: Single point at R=1.4 bohr")
    print("="*60)
    res = h2_4sigma_ci(R=1.4, N_xi_solve=5000, N_grid=50, verbose=True)

    if not np.isnan(res['E_total']):
        D_e_6x6 = 2*E_H - res['E_total']
        D_e_2x2 = 2*E_H - res['E_2x2']
        print(f"\n  Summary at R=1.4:")
        print(f"    D_e(2×2) = {D_e_2x2:.6f} Ha ({D_e_2x2/D_e_exact*100:.1f}% of exact)")
        print(f"    D_e(6×6) = {D_e_6x6:.6f} Ha ({D_e_6x6/D_e_exact*100:.1f}% of exact)")
        print(f"    Improvement: {D_e_6x6 - D_e_2x2:.6f} Ha")

    # --- Test 2: Grid convergence at R=1.4 ---
    print("\n" + "="*60)
    print("  TEST 2: Grid convergence at R=1.4")
    print("="*60)
    for N_grid in [30, 40, 50, 60]:
        res_conv = h2_4sigma_ci(R=1.4, N_xi_solve=5000, N_grid=N_grid, verbose=False)
        if not np.isnan(res_conv['E_total']):
            De = 2*E_H - res_conv['E_total']
            print(f"  N_grid={N_grid:3d}: E_6×6={res_conv['E_total']:.6f}  "
                  f"D_e={De:.6f}  [{res_conv['time']:.1f}s]")

    # --- Test 3: PES scan ---
    print("\n" + "="*60)
    print("  TEST 3: H2 PES scan (4σ CI)")
    print("="*60)
    R_vals = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0])
    pes = scan_h2_4sigma_pes(
        R_vals, N_xi_solve=5000, N_grid=50,
        xi_max_grid=15.0, verbose=False,
    )

    print(f"\n  {'R':>6s} {'E_6×6':>12s} {'E_2×2':>12s} {'E_HF':>12s} "
          f"{'D_e(6×6)':>10s} {'D_e(2×2)':>10s}")
    print("-" * 70)
    for i, R in enumerate(pes['R']):
        E6 = pes['E_total'][i]
        E2 = pes['E_2x2'][i]
        Ehf = pes['E_HF'][i]
        De6 = 2*E_H - E6
        De2 = 2*E_H - E2
        print(f"  {R:6.3f} {E6:12.6f} {E2:12.6f} {Ehf:12.6f} "
              f"{De6:10.6f} {De2:10.6f}")

    # Fit spectroscopic constants
    valid = ~np.isnan(pes['E_total'])
    if np.sum(valid) >= 5:
        fit_6x6 = fit_spectroscopic_constants(pes['R'][valid], pes['E_total'][valid])
        fit_2x2 = fit_spectroscopic_constants(pes['R'][valid], pes['E_2x2'][valid])
        fit_hf = fit_spectroscopic_constants(pes['R'][valid], pes['E_HF'][valid])

        print(f"\n  Spectroscopic constants:")
        print(f"  {'':20s} {'6×6 CI':>12s} {'2×2 CI':>12s} {'HF':>12s} {'Exact':>12s}")
        R6 = fit_6x6['R_eq']
        R2 = fit_2x2['R_eq']
        Rhf = fit_hf['R_eq']
        E6 = fit_6x6['E_min']
        E2 = fit_2x2['E_min']
        Ehf = fit_hf['E_min']
        D6 = 2*E_H - E6
        D2 = 2*E_H - E2
        Dhf = 2*E_H - Ehf

        print(f"  {'R_eq (bohr)':20s} {R6:12.4f} {R2:12.4f} {Rhf:12.4f} {R_eq_exact:12.4f}")
        print(f"  {'E_min (Ha)':20s} {E6:12.6f} {E2:12.6f} {Ehf:12.6f} {E_H2_exact:12.6f}")
        print(f"  {'D_e (Ha)':20s} {D6:12.6f} {D2:12.6f} {Dhf:12.6f} {D_e_exact:12.6f}")
        print(f"  {'D_e (% exact)':20s} {D6/D_e_exact*100:11.1f}% {D2/D_e_exact*100:11.1f}% "
              f"{Dhf/D_e_exact*100:11.1f}%")

        print(f"\n  Errors vs exact:")
        print(f"    R_eq: 6×6={abs(R6-R_eq_exact)/R_eq_exact*100:.1f}%  "
              f"2×2={abs(R2-R_eq_exact)/R_eq_exact*100:.1f}%")
        print(f"    D_e:  6×6={abs(D6-D_e_exact)/D_e_exact*100:.1f}%  "
              f"2×2={abs(D2-D_e_exact)/D_e_exact*100:.1f}%")

        # Verdict
        print(f"\n  {'='*55}")
        if D6 > D2:
            print(f"  4σ CI IMPROVES over minimal basis!")
            print(f"    D_e improvement: {D6-D2:.4f} Ha ({(D6-D2)/D2*100:.1f}%)")
        else:
            print(f"  4σ CI shows NO improvement (unexpected)")
        if D6 > 0:
            print(f"  H2 IS BOUND (D_e = {D6:.4f} Ha, {D6/D_e_exact*100:.1f}% of exact)")
        print(f"  {'='*55}")

    # --- Save results ---
    outpath = os.path.join(os.path.dirname(__file__), 'data', 'prolate_h2_4sigma_ci.txt')
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write("# H2 4σ-Orbital CI on Prolate Spheroidal Lattice\n")
        f.write(f"# Date: 2026-03-13\n")
        f.write(f"# Orbitals: 1σ_g, 1σ_u, 2σ_g, 2σ_u  |  CI: 6×6 ¹Σ_g⁺\n")
        f.write(f"# N_xi_solve=5000, N_grid=50\n")
        f.write(f"# E(H) = {E_H:.6f} Ha\n")
        f.write(f"#\n")
        f.write(f"# {'R':>8s} {'E_6x6':>14s} {'E_2x2':>14s} {'E_HF':>14s} "
                f"{'D_e_6x6':>14s} {'D_e_2x2':>14s}\n")
        for i, R in enumerate(pes['R']):
            E6 = pes['E_total'][i]
            E2 = pes['E_2x2'][i]
            Ehf = pes['E_HF'][i]
            f.write(f"  {R:8.4f} {E6:14.8f} {E2:14.8f} {Ehf:14.8f} "
                    f"{2*E_H-E6:14.8f} {2*E_H-E2:14.8f}\n")

        if np.sum(valid) >= 5:
            f.write(f"#\n# Spectroscopic constants:\n")
            f.write(f"#   6×6 CI: R_eq={R6:.4f} bohr, E_min={E6:.6f} Ha, D_e={D6:.6f} Ha\n")
            f.write(f"#   2×2 CI: R_eq={R2:.4f} bohr, E_min={E2:.6f} Ha, D_e={D2:.6f} Ha\n")
            f.write(f"#   Exact:  R_eq={R_eq_exact:.4f} bohr, D_e={D_e_exact:.6f} Ha\n")

    print(f"\n  Results saved to {outpath}")


if __name__ == '__main__':
    main()
