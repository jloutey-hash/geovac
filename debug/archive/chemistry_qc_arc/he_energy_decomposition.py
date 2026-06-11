"""
He ground-state energy decomposition: graph-native CI.

Decomposes E_total = <h1> + <V_ee> at each n_max (1..7) for the
graph-native FCI solver in casimir_ci.py.

- h1 = graph Laplacian one-body operator (kinetic + nuclear combined)
- V_ee = two-electron Slater integrals at k=Z=2

Also examines the h1 eigenvalue spectrum vs exact -Z^2/(2n^2).
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from geovac.casimir_ci import (
    _build_graph_h1,
    two_electron_integral,
)


def build_separate_fci_matrices(Z: int, n_max: int):
    """Build separate h1-only and V_ee-only FCI matrices for singlet M_L=0.

    Returns H_h1, H_vee, configs, orbitals such that
    H_total = H_h1 + H_vee (verified to machine precision).
    """
    h1_mat, orbitals = _build_graph_h1(Z, n_max)
    n_spatial = len(orbitals)
    k_orb = float(Z)

    # Build singlet M_L=0 configurations (same as build_graph_native_fci)
    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if orbitals[i][2] + orbitals[j][2] == 0:
                configs.append((i, j))

    n_configs = len(configs)
    H_h1 = np.zeros((n_configs, n_configs))
    H_vee = np.zeros((n_configs, n_configs))

    parity = 1.0  # singlet

    def g_int(a, b, c, d):
        na, la, ma = orbitals[a]
        nb, lb, mb = orbitals[b]
        nc, lc, mc = orbitals[c]
        nd, ld, md = orbitals[d]
        return two_electron_integral(na, la, ma, nb, lb, mb,
                                     nc, lc, mc, nd, ld, md, k_orb)

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

            me_h1 = 0.0
            me_vee = 0.0
            for a, b, s_bra in bra_perms:
                for c, d, s_ket in ket_perms:
                    sign = s_bra * s_ket
                    if b == d:
                        me_h1 += sign * h1_mat[a, c]
                    if a == c:
                        me_h1 += sign * h1_mat[b, d]
                    me_vee += sign * g_int(a, b, c, d)

            me_h1 /= (N_IJ * N_PQ)
            me_vee /= (N_IJ * N_PQ)

            H_h1[I, J] = me_h1
            H_h1[J, I] = me_h1
            H_vee[I, J] = me_vee
            H_vee[J, I] = me_vee

    return H_h1, H_vee, configs, orbitals


def analyze_h1_spectrum(Z: int, n_max: int):
    """Analyze the h1 eigenvalues vs exact -Z^2/(2n^2)."""
    h1_mat, orbitals = _build_graph_h1(Z, n_max)

    # Block-diagonalize by (l, m) sector
    sectors = {}
    for idx, (n, l, m) in enumerate(orbitals):
        key = (l, m)
        if key not in sectors:
            sectors[key] = []
        sectors[key].append(idx)

    all_evals = []
    for (l, m), indices in sorted(sectors.items()):
        block = h1_mat[np.ix_(indices, indices)]
        evals = np.linalg.eigvalsh(block)
        for ev in evals:
            all_evals.append((l, m, ev))

    # Exact hydrogenic eigenvalues
    exact_evals = []
    for n, l, m in orbitals:
        exact_evals.append(-Z * Z / (2.0 * n * n))

    return all_evals, exact_evals, orbitals


def main():
    Z = 2
    E_exact = -2.903724  # Ha (exact non-relativistic He)

    # Exact component values for reference
    T_exact = 2.9037  # Ha
    Vne_exact = -6.7530  # Ha
    Vee_exact = 0.9456  # Ha
    vee_ratio_exact = Vee_exact / abs(E_exact)

    print("=" * 90)
    print("He Ground-State Energy Decomposition: Graph-Native CI")
    print("=" * 90)
    print()
    print(f"Exact He: E = {E_exact:.6f} Ha, T = {T_exact:.4f}, "
          f"V_ne = {Vne_exact:.4f}, V_ee = {Vee_exact:.4f}")
    print(f"Exact <V_ee>/|E| = {vee_ratio_exact:.4f}")
    print()

    # --- Part 1: Energy decomposition table ---
    print("-" * 90)
    print(f"{'n_max':>5} {'n_cfg':>6} {'E_total':>12} {'<h1>':>12} "
          f"{'<V_ee>':>12} {'<V_ee>/|E|':>12} {'exact_ratio':>12} "
          f"{'E_err%':>8}")
    print("-" * 90)

    results = []
    for n_max in range(1, 8):
        print(f"  Building n_max={n_max}...", end="", flush=True)
        H_h1, H_vee, configs, orbitals = build_separate_fci_matrices(Z, n_max)
        H_total = H_h1 + H_vee

        # Verify against build_graph_native_fci
        from geovac.casimir_ci import build_graph_native_fci
        H_ref = build_graph_native_fci(Z, n_max)
        max_diff = np.max(np.abs(H_total - H_ref))
        assert max_diff < 1e-12, f"Decomposition error: {max_diff}"

        # Solve for ground state
        evals, evecs = np.linalg.eigh(H_total)
        E_gs = evals[0]
        psi = evecs[:, 0]

        # Expectation values
        h1_exp = psi @ H_h1 @ psi
        vee_exp = psi @ H_vee @ psi
        E_check = h1_exp + vee_exp

        # Verification
        assert abs(E_check - E_gs) < 1e-12, \
            f"Decomposition check failed: {E_check} vs {E_gs}"

        n_configs = len(configs)
        vee_ratio = vee_exp / abs(E_gs)
        err_pct = (E_gs - E_exact) / abs(E_exact) * 100

        results.append({
            'n_max': n_max,
            'n_configs': n_configs,
            'E_total': E_gs,
            'h1_exp': h1_exp,
            'vee_exp': vee_exp,
            'vee_ratio': vee_ratio,
            'err_pct': err_pct,
        })

        print(f"\r{n_max:>5} {n_configs:>6} {E_gs:>12.6f} {h1_exp:>12.6f} "
              f"{vee_exp:>12.6f} {vee_ratio:>12.4f} {vee_ratio_exact:>12.4f} "
              f"{err_pct:>+8.3f}")

    print("-" * 90)
    print()

    # --- Part 2: Convergence analysis ---
    print("=" * 70)
    print("Convergence Analysis")
    print("=" * 70)
    print()
    print(f"{'n_max':>5} {'delta_E':>12} {'delta_h1':>12} {'delta_Vee':>12} "
          f"{'Vee_conv':>10}")
    print("-" * 70)
    for i in range(1, len(results)):
        r = results[i]
        rp = results[i - 1]
        dE = r['E_total'] - rp['E_total']
        dh1 = r['h1_exp'] - rp['h1_exp']
        dvee = r['vee_exp'] - rp['vee_exp']
        # Ratio: how much of delta_E comes from V_ee?
        if abs(dE) > 1e-15:
            vee_frac = dvee / dE * 100
        else:
            vee_frac = 0.0
        print(f"{r['n_max']:>5} {dE:>12.6f} {dh1:>12.6f} {dvee:>12.6f} "
              f"{vee_frac:>9.1f}%")
    print("-" * 70)
    print("  Vee_conv = fraction of delta_E from V_ee change")
    print()

    # --- Part 3: Distance from exact component values ---
    print("=" * 70)
    print("Distance from Exact Component Values")
    print("=" * 70)
    print()
    print("On the graph, h1 encodes T+V_ne combined. We cannot separate them.")
    print("But we can check the virial-like ratio <V_ee>/|E|.")
    print()
    print(f"  Exact <V_ee>/|E| = {vee_ratio_exact:.6f}")
    for r in results:
        print(f"  n_max={r['n_max']}: <V_ee>/|E| = {r['vee_ratio']:.6f}  "
              f"(delta = {r['vee_ratio'] - vee_ratio_exact:+.6f})")
    print()

    # Check: <h1> + 2*<V_ee> ratio (virial-related)
    print("=" * 70)
    print("Virial-Related Check: <h1> + 2*<V_ee>")
    print("=" * 70)
    print()
    print("For exact He: <h1> + 2*<V_ee> = (T+V_ne) + 2*V_ee")
    print(f"             = {T_exact + Vne_exact:.4f} + 2*{Vee_exact:.4f}")
    print(f"             = {T_exact + Vne_exact + 2 * Vee_exact:.4f}")
    print(f"  Since virial: 2T + V = 0, i.e. V = -2T = 2E, so")
    print(f"  V_ne + V_ee = 2E => h1 - T + V_ee = 2E => h1 + V_ee = 2E + T = 2E - E = E")
    print(f"  That's just E_total = h1 + V_ee (trivial).")
    print()
    print("  Better check: virial says T = -E, so <h1> = T + V_ne = -E + V_ne")
    print(f"  Exact: <h1> = {T_exact + Vne_exact:.4f} Ha")
    print(f"  Exact: <V_ee> = {Vee_exact:.4f} Ha")
    print()

    # --- Part 4: h1 eigenvalue analysis at n_max=3 ---
    print("=" * 70)
    print("h1 Eigenvalue Analysis at n_max=3 (Z=2)")
    print("=" * 70)
    print()

    h1_evals, exact_evals, orbs = analyze_h1_spectrum(Z, 3)
    print(f"{'Orbital':>12} {'h1_eval':>12} {'exact -Z^2/2n^2':>16} {'diff':>12}")
    print("-" * 56)

    # Sort by eigenvalue for comparison
    # Group by (l,m) and sort within each group
    sectors = {}
    for i, (n, l, m) in enumerate(orbs):
        key = (l, m)
        if key not in sectors:
            sectors[key] = []
        sectors[key].append(i)

    h1_mat, _ = _build_graph_h1(Z, 3)
    for (l, m), indices in sorted(sectors.items()):
        block = h1_mat[np.ix_(indices, indices)]
        evals = np.linalg.eigvalsh(block)
        for ev_idx, ev in enumerate(evals):
            n_val = min(i_orb for i_orb in indices
                        if orbs[i_orb][1] == l)  # not great, use index
            # The n values in this sector
            n_vals = sorted([orbs[i][0] for i in indices])
            n_here = n_vals[ev_idx]
            exact_e = -Z * Z / (2.0 * n_here * n_here)
            diff = ev - exact_e
            print(f"  ({n_here},{l},{m:+d}){'':<4} {ev:>12.6f} {exact_e:>16.6f} "
                  f"{diff:>12.6f}")
    print()

    # Also print the full h1 matrix at n_max=2 for inspection
    print("=" * 70)
    print("Full h1 matrix at n_max=2 (Z=2)")
    print("=" * 70)
    print()
    h1_2, orbs_2 = _build_graph_h1(Z, 2)
    print("Orbitals:", [(n, l, m) for n, l, m in orbs_2])
    print()
    print("h1 matrix:")
    np.set_printoptions(precision=6, linewidth=120, suppress=True)
    print(h1_2)
    print()
    print("Eigenvalues:", np.sort(np.linalg.eigvalsh(h1_2)))
    exact_2 = sorted([-Z * Z / (2.0 * n * n) for n, l, m in orbs_2])
    print("Exact -Z^2/2n^2:", exact_2)
    print()

    # --- Part 5: Summary table for paper ---
    print("=" * 90)
    print("SUMMARY TABLE (for potential paper use)")
    print("=" * 90)
    print()
    print(f"{'n_max':>5} | {'N_cfg':>5} | {'E (Ha)':>11} | "
          f"{'<h1> (Ha)':>11} | {'<V_ee> (Ha)':>11} | "
          f"{'<V_ee>/|E|':>10} | {'err (%)':>8}")
    print("-" * 80)
    for r in results:
        print(f"{r['n_max']:>5} | {r['n_configs']:>5} | "
              f"{r['E_total']:>11.6f} | {r['h1_exp']:>11.6f} | "
              f"{r['vee_exp']:>11.6f} | {r['vee_ratio']:>10.6f} | "
              f"{r['err_pct']:>+8.3f}")
    print("-" * 80)
    print(f"{'exact':>5} | {'':>5} | {E_exact:>11.6f} | "
          f"{T_exact + Vne_exact:>11.4f} | {Vee_exact:>11.4f} | "
          f"{vee_ratio_exact:>10.6f} |")
    print()


if __name__ == '__main__':
    main()
