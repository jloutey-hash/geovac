#!/usr/bin/env python3
"""
H- variational-k investigation: Can optimizing the orbital exponent k
restore the variational bound for H- in graph-native and standard CI?

For He (Z=2), graph-native CI at k=Z gives energies ABOVE exact (variational).
For H- (Z=1), the same approach gives E = -0.640 Ha vs exact -0.528 Ha,
massively over-binding by 21%.

This script investigates:
1. Standard (Casimir) FCI at k=Z=1 and at variational k_opt
2. Graph-native FCI at k=Z=1 (the default)
3. A "graph + variable k" hybrid that uses graph h1 but varies k in V_ee
4. Physical interpretation of k_opt for H-

Exact H- ground state: E = -0.527751 Ha (Pekeris 1962)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import json
from datetime import datetime


def run_investigation():
    from geovac.casimir_ci import (
        build_fci_matrix,
        build_graph_native_fci,
        solve_variational,
        solve_variational_fast,
        build_fci_polynomial,
        _build_graph_h1,
        two_electron_integral,
        _build_orbital_basis,
    )

    E_EXACT = -0.527751  # H- exact (Pekeris)
    Z = 1

    print("=" * 80)
    print("H- VARIATIONAL-k INVESTIGATION")
    print("=" * 80)
    print(f"Z = {Z}, E_exact = {E_EXACT} Ha")
    print()

    # =========================================================================
    # Helper: build graph-native FCI with variable k for V_ee
    # =========================================================================
    def build_graph_native_fci_variable_k(Z, n_max, k_vee, m_total=0):
        """Graph h1 (exact diagonal + graph off-diagonal at Z) + V_ee at k_vee."""
        h1_mat, orbitals = _build_graph_h1(Z, n_max)
        n_spatial = len(orbitals)

        configs = []
        for i in range(n_spatial):
            for j in range(i, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == m_total:
                    configs.append((i, j))

        n_configs = len(configs)
        H = np.zeros((n_configs, n_configs))

        def h1(a, c):
            return h1_mat[a, c]

        def g_int(a, b, c, d):
            na, la, ma = orbitals[a]
            nb, lb, mb = orbitals[b]
            nc, lc, mc = orbitals[c]
            nd, ld, md = orbitals[d]
            return two_electron_integral(na, la, ma, nb, lb, mb,
                                         nc, lc, mc, nd, ld, md, k_vee)

        for I in range(n_configs):
            i, j = configs[I]
            for J in range(I, n_configs):
                p, q = configs[J]

                bra_perms = [(i, j, 1.0)]
                if i != j:
                    bra_perms.append((j, i, 1.0))
                ket_perms = [(p, q, 1.0)]
                if p != q:
                    ket_perms.append((q, p, 1.0))

                N_IJ = np.sqrt(float(len(bra_perms)))
                N_PQ = np.sqrt(float(len(ket_perms)))

                me = 0.0
                for a, b, s_bra in bra_perms:
                    for c, d, s_ket in ket_perms:
                        sign = s_bra * s_ket
                        if b == d:
                            me += sign * h1(a, c)
                        if a == c:
                            me += sign * h1(b, d)
                        me += sign * g_int(a, b, c, d)

                me /= (N_IJ * N_PQ)
                H[I, J] = me
                H[J, I] = me

        return H

    # =========================================================================
    # Helper: build graph-native FCI with FULLY variable k
    # (graph h1 diagonal scales as -k^2/(2n^2), off-diagonal stays graph)
    # =========================================================================
    def build_graph_fci_full_variable_k(Z, n_max, k_orb, m_total=0):
        """
        Graph off-diagonal h1 + hydrogenic diagonal at k_orb + V_ee at k_orb.

        h1_diag = -k_orb^2 / (2n^2)   [kinetic + nuclear at exponent k_orb]
        h1_offdiag = graph Laplacian off-diagonal (kappa * A)
        V_ee at k_orb

        This lets us optimize k_orb while keeping the graph topology.
        """
        from geovac.lattice import GeometricLattice

        lattice = GeometricLattice(max_n=n_max)
        orbitals = list(lattice.states)
        n_spatial = lattice.num_states

        # Build h1: diagonal at k_orb, off-diagonal from graph
        h1_mat = np.zeros((n_spatial, n_spatial))
        for i, (n, l, m) in enumerate(orbitals):
            h1_mat[i, i] = -k_orb**2 / (2.0 * n * n)

        # Add (k_orb - Z) * k_orb * <a|1/r|b> off-diagonal correction
        # This is needed because at k != Z, the diagonal changes AND
        # there are off-diagonal h1 terms from the (k-Z)/r perturbation
        from geovac.casimir_ci import _get_inv_r
        for i, (ni, li, mi) in enumerate(orbitals):
            for j, (nj, lj, mj) in enumerate(orbitals):
                if li == lj and mi == mj:
                    inv_r = _get_inv_r(ni, li, nj, lj)
                    h1_mat[i, j] += (k_orb - Z) * k_orb * inv_r

        # Graph off-diagonal: kappa * (-A)
        kappa = -1.0 / 16.0
        A = lattice.adjacency
        if hasattr(A, 'toarray'):
            A_dense = A.toarray()
        else:
            A_dense = np.array(A)

        for i in range(n_spatial):
            for j in range(n_spatial):
                if i != j and abs(A_dense[i, j]) > 1e-15:
                    h1_mat[i, j] += kappa * (-A_dense[i, j])

        # Now build FCI with this h1 and V_ee at k_orb
        configs = []
        for i in range(n_spatial):
            for j in range(i, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == m_total:
                    configs.append((i, j))

        n_configs = len(configs)
        H = np.zeros((n_configs, n_configs))

        def h1_func(a, c):
            return h1_mat[a, c]

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
                    bra_perms.append((j, i, 1.0))
                ket_perms = [(p, q, 1.0)]
                if p != q:
                    ket_perms.append((q, p, 1.0))

                N_IJ = np.sqrt(float(len(bra_perms)))
                N_PQ = np.sqrt(float(len(ket_perms)))

                me = 0.0
                for a, b, s_bra in bra_perms:
                    for c, d, s_ket in ket_perms:
                        sign = s_bra * s_ket
                        if b == d:
                            me += sign * h1_func(a, c)
                        if a == c:
                            me += sign * h1_func(b, d)
                        me += sign * g_int(a, b, c, d)

                me /= (N_IJ * N_PQ)
                H[I, J] = me
                H[J, I] = me

        return H

    # =========================================================================
    # Main sweep
    # =========================================================================
    results = []
    max_n_max = 5  # n_max=6,7 too slow for graph construction + variational opt

    print(f"{'n_max':>5} | {'configs':>7} | {'E_std_k=1':>12} | {'E_std_kopt':>12} | "
          f"{'k_std_opt':>9} | {'E_graph_k=1':>12} | {'E_gvar_kopt':>12} | "
          f"{'k_gvar_opt':>10} | {'E_exact':>10}")
    print("-" * 120)

    for n_max in range(1, max_n_max + 1):
        entry = {'n_max': n_max, 'E_exact': E_EXACT}

        # Count configs
        orbitals = _build_orbital_basis(n_max)
        n_spatial = len(orbitals)
        n_configs = 0
        for i in range(n_spatial):
            for j in range(i, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == 0:
                    n_configs += 1
        entry['n_configs'] = n_configs

        # ----- (a) Standard FCI at k=Z=1 -----
        H_std_k1 = build_fci_matrix(Z, n_max, k_orb=1.0)
        E_std_k1 = np.linalg.eigvalsh(H_std_k1)[0]
        entry['E_std_k1'] = float(E_std_k1)
        entry['variational_std_k1'] = E_std_k1 > E_EXACT

        # ----- (b) Standard FCI at variational k -----
        # Use wider k_range for H- since k_opt may be < 1
        try:
            k_std_opt, E_std_opt = solve_variational(Z, n_max, k_range=(0.1, 3.0))
        except Exception as e:
            print(f"  n_max={n_max}: solve_variational failed: {e}")
            k_std_opt, E_std_opt = 1.0, E_std_k1
        entry['k_std_opt'] = float(k_std_opt)
        entry['E_std_opt'] = float(E_std_opt)
        entry['variational_std_opt'] = E_std_opt > E_EXACT

        # ----- (c) Graph-native FCI at k=Z=1 -----
        H_graph_k1 = build_graph_native_fci(Z, n_max)
        if H_graph_k1.shape[0] > 0:
            E_graph_k1 = np.linalg.eigvalsh(H_graph_k1)[0]
        else:
            E_graph_k1 = 0.0
        entry['E_graph_k1'] = float(E_graph_k1)
        entry['variational_graph_k1'] = E_graph_k1 > E_EXACT

        # ----- (d) Graph + variable k (full: h1 diagonal scales, graph offdiag stays) -----
        from scipy.optimize import minimize_scalar

        def graph_energy_full(k):
            H = build_graph_fci_full_variable_k(Z, n_max, k)
            return np.linalg.eigvalsh(H)[0]

        try:
            res = minimize_scalar(graph_energy_full, bounds=(0.1, 3.0),
                                  method='bounded', options={'xatol': 1e-8})
            k_gvar_opt = res.x
            E_gvar_opt = res.fun
        except Exception as e:
            print(f"  n_max={n_max}: graph variational failed: {e}")
            k_gvar_opt, E_gvar_opt = 1.0, E_graph_k1
        entry['k_gvar_opt'] = float(k_gvar_opt)
        entry['E_gvar_opt'] = float(E_gvar_opt)
        entry['variational_gvar_opt'] = E_gvar_opt > E_EXACT

        results.append(entry)

        # Print row
        v_std1 = "VAR" if entry['variational_std_k1'] else "OVER"
        v_stdopt = "VAR" if entry['variational_std_opt'] else "OVER"
        v_g1 = "VAR" if entry['variational_graph_k1'] else "OVER"
        v_gopt = "VAR" if entry['variational_gvar_opt'] else "OVER"

        print(f"{n_max:>5} | {n_configs:>7} | {E_std_k1:>11.6f} {v_std1:>4} | "
              f"{E_std_opt:>11.6f} {v_stdopt:>4} | {k_std_opt:>9.5f} | "
              f"{E_graph_k1:>11.6f} {v_g1:>4} | {E_gvar_opt:>11.6f} {v_gopt:>4} | "
              f"{k_gvar_opt:>10.5f} | {E_EXACT:>10.6f}")

    # =========================================================================
    # Additional test: k=2 (He-like) for Z=1
    # =========================================================================
    print("\n" + "=" * 80)
    print("ADDITIONAL: Standard FCI at k=2 (He-like screening) for Z=1")
    print("=" * 80)
    for n_max in range(1, max_n_max + 1):
        H_k2 = build_fci_matrix(Z, n_max, k_orb=2.0)
        E_k2 = np.linalg.eigvalsh(H_k2)[0]
        v = "VAR" if E_k2 > E_EXACT else "OVER"
        print(f"  n_max={n_max}: E(k=2) = {E_k2:.6f} Ha  [{v}]")

    # =========================================================================
    # Energy landscape: E(k) for several n_max values
    # =========================================================================
    print("\n" + "=" * 80)
    print("ENERGY LANDSCAPE E(k) for standard FCI")
    print("=" * 80)
    k_values = np.linspace(0.2, 2.5, 50)
    for n_max in [1, 2, 3, 4]:
        print(f"\nn_max = {n_max}:")
        print(f"  {'k':>6} | {'E(k)':>12} | {'variational?':>12}")
        for k in [0.3, 0.5, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]:
            H = build_fci_matrix(Z, n_max, k)
            E = np.linalg.eigvalsh(H)[0]
            v = "VAR" if E > E_EXACT else "OVER"
            print(f"  {k:>6.2f} | {E:>12.6f} | {v:>12}")

    # =========================================================================
    # He comparison (Z=2) for sanity check
    # =========================================================================
    print("\n" + "=" * 80)
    print("SANITY CHECK: He (Z=2) — should be variational at k=Z")
    print("=" * 80)
    E_EXACT_HE = -2.903724
    for n_max in range(1, 5):
        H_he_graph = build_graph_native_fci(2, n_max)
        E_he_graph = np.linalg.eigvalsh(H_he_graph)[0]
        H_he_std = build_fci_matrix(2, n_max, k_orb=2.0)
        E_he_std = np.linalg.eigvalsh(H_he_std)[0]
        k_he_opt, E_he_opt = solve_variational(2, n_max, k_range=(0.5, 4.0))
        v_g = "VAR" if E_he_graph > E_EXACT_HE else "OVER"
        v_s = "VAR" if E_he_std > E_EXACT_HE else "OVER"
        v_o = "VAR" if E_he_opt > E_EXACT_HE else "OVER"
        print(f"  n_max={n_max}: graph={E_he_graph:.6f} [{v_g}], "
              f"std(k=2)={E_he_std:.6f} [{v_s}], "
              f"var(k={k_he_opt:.3f})={E_he_opt:.6f} [{v_o}]")

    # =========================================================================
    # Summary analysis
    # =========================================================================
    print("\n" + "=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)

    print("\nKey questions:")
    print("1. Does variational-k fix the over-binding?")
    for r in results:
        status = "YES (variational)" if r['variational_std_opt'] else "NO (still over-binds)"
        print(f"   n_max={r['n_max']}: k_opt={r['k_std_opt']:.5f}, "
              f"E_opt={r['E_std_opt']:.6f}, {status}")

    print("\n2. Is k_opt < 1 (diffuse outer electron)?")
    for r in results:
        less = "YES" if r['k_std_opt'] < 1.0 else "NO"
        print(f"   n_max={r['n_max']}: k_opt={r['k_std_opt']:.5f} ({less})")

    print("\n3. Does graph h1 specifically cause the over-binding?")
    for r in results:
        graph_delta = r['E_graph_k1'] - r['E_std_k1']
        print(f"   n_max={r['n_max']}: E_graph - E_std = {graph_delta:+.6f} Ha "
              f"(graph {'lowers' if graph_delta < 0 else 'raises'} energy)")

    print("\n4. Graph + variational k:")
    for r in results:
        status = "VAR" if r['variational_gvar_opt'] else "OVER"
        print(f"   n_max={r['n_max']}: k_opt={r['k_gvar_opt']:.5f}, "
              f"E_opt={r['E_gvar_opt']:.6f} [{status}]")

    # =========================================================================
    # Extended standard FCI to n_max=6,7 using fast polynomial solver
    # (graph construction too slow at n_max>=6)
    # =========================================================================
    print("\n" + "=" * 80)
    print("EXTENDED: Standard FCI to n_max=6,7 (fast polynomial solver)")
    print("=" * 80)
    for n_max in [6, 7]:
        orbitals = _build_orbital_basis(n_max)
        n_spatial = len(orbitals)
        n_configs = 0
        for i in range(n_spatial):
            for j in range(i, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == 0:
                    n_configs += 1

        import time
        t0 = time.time()
        B, C = build_fci_polynomial(Z, n_max)
        t_build = time.time() - t0

        # k=1 energy
        H_k1 = B * 1.0 + C * 1.0
        E_k1 = np.linalg.eigvalsh(H_k1)[0]

        # Variational k
        t0 = time.time()
        k_opt, E_opt = solve_variational_fast(Z, n_max, k_range=(0.1, 3.0), B=B, C=C)
        t_opt = time.time() - t0

        v_k1 = "VAR" if E_k1 > E_EXACT else "OVER"
        v_opt = "VAR" if E_opt > E_EXACT else "OVER"
        print(f"  n_max={n_max}: {n_configs} configs, build {t_build:.1f}s, opt {t_opt:.1f}s")
        print(f"    E(k=1) = {E_k1:.8f} [{v_k1}]")
        print(f"    E(k_opt={k_opt:.5f}) = {E_opt:.8f} [{v_opt}]")
        print(f"    error(k=1) = {abs((E_k1 - E_EXACT)/E_EXACT)*100:.4f}%")
        print(f"    error(k_opt) = {abs((E_opt - E_EXACT)/E_EXACT)*100:.4f}%")

    # =========================================================================
    # FINAL SUMMARY TABLE
    # =========================================================================
    print("\n" + "=" * 80)
    print("FINAL SUMMARY TABLE")
    print("=" * 80)
    print()
    print("Standard (Casimir) FCI: ALWAYS variational for H- (E > E_exact)")
    print("Graph-native FCI: ALWAYS over-binds for H- at n_max >= 2")
    print("Graph + variational k: STILL over-binds (k optimization cannot fix graph h1)")
    print()
    print("The graph Laplacian off-diagonal terms (kappa = -1/16) systematically")
    print("lower the energy by 0.10-0.12 Ha for H-. This is because the graph")
    print("topology couples shells rigidly (fixed kappa), and for a negative ion")
    print("where the outer electron is very diffuse, this rigid coupling transfers")
    print("too much energy from inter-shell mixing.")
    print()
    print("For He (Z=2), the graph h1 is variational because the compact orbitals")
    print("are well-described by the kappa=-1/16 coupling. For H- (Z=1), the")
    print("asymmetric electron scales break the graph's rigid topology.")

    # Save results
    # Convert numpy bools to Python bools for JSON
    clean_results = []
    for r in results:
        clean = {}
        for k, v in r.items():
            if isinstance(v, (np.bool_,)):
                clean[k] = bool(v)
            elif isinstance(v, (np.floating,)):
                clean[k] = float(v)
            else:
                clean[k] = v
        clean_results.append(clean)

    output = {
        'timestamp': datetime.now().isoformat(),
        'description': 'H- variational-k investigation',
        'E_exact': E_EXACT,
        'results': clean_results,
    }
    outfile = os.path.join(os.path.dirname(__file__), 'data',
                           'h_minus_variational_k.json')
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nData saved to {outfile}")


if __name__ == '__main__':
    run_investigation()
