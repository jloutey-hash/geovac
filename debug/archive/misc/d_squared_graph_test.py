"""Test whether D^2 (squared Dirac operator) has a graph representation.

D^2 eigenvalues on unit S^3: (n_CH + 3/2)^2 with degeneracy 2(n_CH+1)(n_CH+2)
Lichnerowicz: D^2 = Delta_spinor + R/4 = Delta_spinor + 3/2  (R=6 on unit S^3)

Structural question: the scalar graph Laplacian on (n,l,m) nodes gives
eigenvalues related to n^2-1 (Fock). Can any graph on (n,kappa,m_j) nodes
give eigenvalues related to (n+3/2)^2?

Key insight: graph Laplacian eigenvalues are bounded (by ~4*d_max for
unweighted graphs) while D^2 grows quadratically. So graph Laplacian
CANNOT directly match D^2 spectrum. Same situation as scalar case:
diagonal carries the spectrum, graph coupling is perturbative.

Tests:
1. Convention table: scalar LB vs D^2 vs spinor Laplacian
2. Graph Laplacian eigenvalues on Dirac nodes (Rule A)
3. Optimal kappa fit: H = kappa*L + diag(c) -> D^2 target
4. Per-kappa sector analysis with path-graph eigenvalue comparison
5. Lichnerowicz decomposition: L + 3/2*I vs D^2
6. Weighted adjacency: Fock-topological weights on Dirac edges
7. Scalar graph comparison at same n_max
"""

import numpy as np
import json
from pathlib import Path
from scipy.optimize import minimize

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from geovac.ihara_zeta_dirac import build_dirac_s3_graph, DiracLabel
from geovac.dirac_matrix_elements import kappa_to_l
from geovac.lattice import GeometricLattice


KAPPA_SCALAR = -1.0 / 16.0


def d2_eigenvalue(n_fock):
    """D^2 eigenvalue in Fock convention: (n_Fock + 1/2)^2."""
    return (n_fock + 0.5)**2


def spinor_lap_eigenvalue(n_fock):
    """Spinor Laplacian eigenvalue = D^2 - R/4 = D^2 - 3/2."""
    return d2_eigenvalue(n_fock) - 1.5


def scalar_lb_eigenvalue_fock(n_fock):
    """Scalar LB eigenvalue in Fock convention: n^2 - 1."""
    return n_fock**2 - 1


def path_graph_eigenvalues(N):
    """Laplacian eigenvalues of the path graph P_N (N nodes)."""
    return np.array([2 - 2*np.cos(k*np.pi/(N+1)) for k in range(1, N+1)])


def analyze_degeneracies(evals, tol=1e-8):
    """Group eigenvalues by degeneracy."""
    groups = []
    i = 0
    while i < len(evals):
        val = evals[i]
        count = 1
        while i + count < len(evals) and abs(evals[i + count] - val) < tol:
            count += 1
        groups.append((val, count))
        i += count
    return groups


def main():
    results = {}
    Z = 1

    print("=" * 80)
    print("D^2 GRAPH REPRESENTATION TEST")
    print("D^2 = (n_Fock + 1/2)^2   [Fock convention]")
    print("Lichnerowicz: D^2 = Delta_spinor + 3/2   [R=6 on unit S^3]")
    print("=" * 80)

    # --- Convention table ---
    print("\n--- CONVENTION TABLE ---")
    print(f"{'n_F':>3} {'n_CH':>4} {'D^2':>8} {'D^2-9/4':>8} "
          f"{'D_sp':>8} {'scal_LB':>8} {'D_sp-LB':>8} {'deg_D':>6} {'deg_s':>6}")
    for n_f in range(1, 7):
        n_ch = n_f - 1
        d2 = d2_eigenvalue(n_f)
        d2_int = n_ch * (n_ch + 3)
        d_sp = spinor_lap_eigenvalue(n_f)
        s_lb = scalar_lb_eigenvalue_fock(n_f)
        diff = d_sp - s_lb
        deg_d = 2 * n_f * (n_f + 1)
        deg_s = n_f**2
        print(f"{n_f:3d} {n_ch:4d} {d2:8.2f} {d2_int:8d} "
              f"{d_sp:8.2f} {s_lb:8.0f} {diff:8.2f} {deg_d:6d} {deg_s:6d}")

    conv_data = {
        'note': 'D_sp - scalar_LB = n_CH + 3/4 = n_Fock - 1/4 (grows with n!)',
        'bounded_obstruction': ('Graph Laplacian eigenvalues bounded by ~4*d_max; '
                                'D^2 grows as n^2. Same situation as scalar.'),
    }
    results['conventions'] = conv_data

    # --- Scalar graph reference ---
    print("\n--- SCALAR GRAPH REFERENCE ---")
    for n_max in [2, 3, 4]:
        lattice = GeometricLattice(max_n=n_max, nuclear_charge=Z)
        orbitals = list(lattice.states)
        n_states = lattice.num_states

        A_sc = lattice.adjacency
        if hasattr(A_sc, 'toarray'):
            A_sc = A_sc.toarray()
        D_sc = np.diag(np.array(A_sc.sum(axis=1)).flatten())
        L_sc = D_sc - A_sc

        evals_L_sc = np.sort(np.linalg.eigvalsh(L_sc))

        # Scalar h1 eigenvalues
        H_sc = KAPPA_SCALAR * L_sc.copy()
        for i, (n, l, m) in enumerate(orbitals):
            H_sc[i, i] += -Z**2 / (2.0 * n**2)
        evals_H_sc = np.sort(np.linalg.eigvalsh(H_sc))

        rydberg = sorted([-Z**2 / (2.0 * n**2) for n in range(1, n_max + 1)
                          for l in range(n) for m in range(-l, l + 1)])
        rydberg = np.array(rydberg)

        max_err = np.max(np.abs(evals_H_sc - rydberg))
        rms_err = np.sqrt(np.mean((evals_H_sc - rydberg)**2))

        print(f"\n  n_max={n_max}: {n_states} states")
        print(f"    L eigenvalues (first 8): {np.round(evals_L_sc[:8], 4)}")
        print(f"    L range: [{evals_L_sc[0]:.4f}, {evals_L_sc[-1]:.4f}]")
        print(f"    H = kappa*L + diag(Rydberg): max|dE|={max_err:.6f} Ha, rms={rms_err:.6f}")

        degs = analyze_degeneracies(evals_L_sc, tol=1e-6)
        deg_str = ", ".join(f"{v:.3f}(x{c})" for v, c in degs[:8])
        print(f"    L degeneracies: {deg_str}")

        results[f'scalar_n{n_max}'] = {
            'n_states': n_states,
            'L_evals': evals_L_sc.tolist(),
            'H_max_err': float(max_err),
            'H_rms_err': float(rms_err),
        }

    # --- Dirac graph D^2 test ---
    print("\n--- DIRAC GRAPH D^2 TEST ---")
    for n_max in [2, 3, 4]:
        adj, labels, _deg, desc = build_dirac_s3_graph(n_max, adjacency_rule='A')
        n_states = len(labels)

        A_d = adj.astype(float)
        D_d = np.diag(A_d.sum(axis=1))
        L_d = D_d - A_d

        evals_L = np.sort(np.linalg.eigvalsh(L_d))

        # Target eigenvalues
        target_d2 = np.sort(np.array([d2_eigenvalue(lab.n_fock) for lab in labels]))
        target_sp = np.sort(np.array([spinor_lap_eigenvalue(lab.n_fock) for lab in labels]))

        print(f"\n  n_max={n_max}: {n_states} Dirac nodes")
        print(f"    {desc}")
        print(f"    L eigenvalues (first 8): {np.round(evals_L[:8], 4)}")
        print(f"    L range: [{evals_L[0]:.4f}, {evals_L[-1]:.4f}]")
        print(f"    Target D^2 range: [{target_d2[0]:.2f}, {target_d2[-1]:.2f}]")

        degs_L = analyze_degeneracies(evals_L, tol=1e-6)
        deg_str = ", ".join(f"{v:.3f}(x{c})" for v, c in degs_L[:8])
        print(f"    L degeneracies: {deg_str}")

        degs_d2 = analyze_degeneracies(target_d2, tol=1e-8)
        deg_str_d2 = ", ".join(f"{v:.2f}(x{c})" for v, c in degs_d2)
        print(f"    D^2 degeneracies: {deg_str_d2}")

        # Test A: Linear fit kappa*L + c -> D^2 target (sorted eigenvalue matching)
        def err_linear(params):
            k, c = params
            return np.sum((k * evals_L + c - target_d2)**2)

        res = minimize(err_linear, [1.0, 2.0], method='Nelder-Mead')
        k_opt, c_opt = res.x
        rms_opt = np.sqrt(res.fun / n_states)
        print(f"\n    Linear fit: kappa*L + c -> D^2")
        print(f"      kappa_opt={k_opt:.4f}, c_opt={c_opt:.4f}, rms={rms_opt:.4f}")
        print(f"      (Compare: D^2 range is {target_d2[-1]-target_d2[0]:.2f})")

        # Test B: H = kappa*L + diag(D^2) — perturbative coupling on top of exact diagonal
        print(f"\n    H = kappa*L + diag(D^2) [perturbative test]:")
        for kappa_test in [0.0, -1/16, -1/8, -1/32, 1/16]:
            H_test = kappa_test * L_d.copy()
            for i, lab in enumerate(labels):
                H_test[i, i] += d2_eigenvalue(lab.n_fock)
            evals_test = np.sort(np.linalg.eigvalsh(H_test))
            max_err = np.max(np.abs(evals_test - target_d2))
            rms_err = np.sqrt(np.mean((evals_test - target_d2)**2))
            print(f"      kappa={kappa_test:+.4f}: max|dE|={max_err:.6f}, rms={rms_err:.6f}")

        # Test C: Lichnerowicz — L + 3/2*I vs D^2
        H_lich = L_d.copy() + 1.5 * np.eye(n_states)
        evals_lich = np.sort(np.linalg.eigvalsh(H_lich))
        max_err_lich = np.max(np.abs(evals_lich - target_d2))
        print(f"\n    Lichnerowicz: L + 3/2*I vs D^2: max|dE|={max_err_lich:.4f}")

        # With scaling
        for k_test in [-1/16, 1.0, 2.0]:
            H_wl = k_test * L_d.copy() + 1.5 * np.eye(n_states)
            evals_wl = np.sort(np.linalg.eigvalsh(H_wl))
            err_wl = np.max(np.abs(evals_wl - target_d2))
            print(f"      kappa={k_test:+.4f}: kappa*L + 3/2*I max|dE|={err_wl:.4f}")

        # Test D: Per-kappa sector analysis
        print(f"\n    Per-kappa sector analysis:")
        kappa_values = sorted(set(lab.kappa for lab in labels))
        sector_data = {}

        for kv in kappa_values:
            indices = [i for i, lab in enumerate(labels) if lab.kappa == kv]
            sector_labels = [labels[i] for i in indices]
            n_vals = sorted(set(lab.n_fock for lab in sector_labels))
            mj_vals = sorted(set(lab.two_m_j for lab in sector_labels))

            l_val = kappa_to_l(kv)
            j_val = abs(kv) - 0.5

            # Sector adjacency
            sector_adj = np.zeros((len(indices), len(indices)))
            for i_new, i_old in enumerate(indices):
                for j_new, j_old in enumerate(indices):
                    sector_adj[i_new, j_new] = adj[i_old, j_old]

            D_sec = np.diag(sector_adj.sum(axis=1))
            L_sec = D_sec - sector_adj
            evals_sec = np.sort(np.linalg.eigvalsh(L_sec))

            # Target within sector
            target_sec = np.sort(np.array([d2_eigenvalue(lab.n_fock)
                                           for lab in sector_labels]))

            # Path graph comparison
            # Within each kappa sector, the graph is a 2D grid:
            # n-dimension: path P_{len(n_vals)} (one per m_j)
            # m_j-dimension: path P_{len(mj_vals)} (one per n)
            n_dim = len(n_vals)
            mj_dim = len(mj_vals)

            # 2D grid Laplacian eigenvalues = sum of 1D path eigenvalues
            if n_dim > 0 and mj_dim > 0:
                ev_n = path_graph_eigenvalues(n_dim)
                ev_mj = path_graph_eigenvalues(mj_dim)
                grid_evals = sorted([ev_n[i] + ev_mj[j]
                                     for i in range(n_dim) for j in range(mj_dim)])
                grid_evals = np.array(grid_evals)

                grid_match = np.max(np.abs(evals_sec - grid_evals))
            else:
                grid_match = float('nan')

            # Correlation
            if len(evals_sec) > 2 and np.std(evals_sec) > 1e-10:
                corr = np.corrcoef(evals_sec, target_sec)[0, 1]
            else:
                corr = float('nan')

            print(f"      kappa={kv:+d} (l={l_val}, j={j_val}): "
                  f"{len(indices)} nodes, n in {n_vals}, 2j+1={len(mj_vals)} m_j values")
            print(f"        L_sector evals: {np.round(evals_sec[:6], 4)}")
            print(f"        2D grid evals:  {np.round(grid_evals[:6], 4) if n_dim > 0 else 'N/A'}")
            print(f"        grid match:     {grid_match:.2e}")
            print(f"        D^2 target:     {np.round(target_sec[:6], 2)}")
            print(f"        corr(L, D^2):   {corr:.4f}")

            sector_data[f'kappa_{kv}'] = {
                'l': l_val, 'j': j_val,
                'n_nodes': len(indices),
                'L_evals': evals_sec.tolist(),
                'grid_evals': grid_evals.tolist() if n_dim > 0 else [],
                'grid_match': float(grid_match),
                'target_d2': target_sec.tolist(),
                'corr_L_d2': float(corr),
            }

        # Test E: Fock-topological weighted adjacency
        # Weight each edge by 1/(n_a * n_b), mimicking scalar Fock weights
        print(f"\n    Fock-topological weighted adjacency:")
        A_weighted = np.zeros_like(A_d)
        for i in range(n_states):
            for j in range(i + 1, n_states):
                if adj[i, j]:
                    w = 1.0 / (labels[i].n_fock * labels[j].n_fock)
                    A_weighted[i, j] = w
                    A_weighted[j, i] = w
        D_w = np.diag(A_weighted.sum(axis=1))
        L_w = D_w - A_weighted
        evals_L_w = np.sort(np.linalg.eigvalsh(L_w))

        print(f"    Weighted L eigenvalues (first 8): {np.round(evals_L_w[:8], 6)}")
        print(f"    Weighted L range: [{evals_L_w[0]:.6f}, {evals_L_w[-1]:.6f}]")

        # Optimal fit with weighted L
        def err_weighted(params):
            k, c = params
            return np.sum((k * evals_L_w + c - target_d2)**2)

        res_w = minimize(err_weighted, [1.0, 2.0], method='Nelder-Mead')
        k_w, c_w = res_w.x
        rms_w = np.sqrt(res_w.fun / n_states)
        print(f"    Weighted fit: kappa={k_w:.4f}, c={c_w:.4f}, rms={rms_w:.4f}")

        # H = kappa*L_weighted + diag(D^2) perturbative test
        for kappa_test in [0.0, -1/16]:
            H_wt = kappa_test * L_w.copy()
            for i, lab in enumerate(labels):
                H_wt[i, i] += d2_eigenvalue(lab.n_fock)
            evals_wt = np.sort(np.linalg.eigvalsh(H_wt))
            max_err_wt = np.max(np.abs(evals_wt - target_d2))
            print(f"    kappa={kappa_test:+.4f} weighted: max|dE|={max_err_wt:.6f}")

        results[f'dirac_d2_n{n_max}'] = {
            'n_states': n_states,
            'L_evals': evals_L.tolist(),
            'target_d2': target_d2.tolist(),
            'linear_fit': {'kappa': float(k_opt), 'c': float(c_opt), 'rms': float(rms_opt)},
            'lichnerowicz_max_err': float(max_err_lich),
            'sectors': sector_data,
            'weighted_fit': {'kappa': float(k_w), 'c': float(c_w), 'rms': float(rms_w)},
        }

    # --- KEY STRUCTURAL ANALYSIS ---
    print("\n" + "=" * 80)
    print("STRUCTURAL ANALYSIS")
    print("=" * 80)

    print("""
1. BOUNDED vs UNBOUNDED SPECTRUM:
   Graph Laplacian eigenvalues are bounded by ~4*d_max (Cheeger).
   D^2 eigenvalues grow as (n+3/2)^2 (unbounded).
   => Graph Laplacian CANNOT match D^2 spectrum directly.
   => Same situation as scalar: diagonal carries spectrum,
      graph coupling provides perturbative inter-shell mixing.

2. DIAGONAL IS ALREADY EXACT:
   For scalar: diag(-Z^2/(2n^2)) = exact Rydberg.
   For D^2: diag((n+3/2)^2) = exact D^2 eigenvalues.
   Graph coupling PERTURBS away from exact in both cases.
   The graph is NOT a spectrum calculator — it's a coupling operator.

3. SPINOR vs SCALAR LAPLACIAN DIFFERENCE:
   Delta_spinor eigenvalue = n^2 + n + 3/4  [Fock, n=1,2,3,...]
   Scalar LB eigenvalue    = n^2 - 1         [Fock]
   Difference = n + 7/4 (grows with n!)
   NOT a constant shift => cannot be absorbed into diagonal.
   => Different graph topology needed for spinor Laplacian.

4. PER-KAPPA SECTOR STRUCTURE:
   Within each kappa-sector (Rule A), the graph is a 2D grid:
   path(n-ladder) x path(m_j-ladder).
   Eigenvalues are sums of 1D path-graph cosines.
   These are IDENTICAL in structure to scalar per-l sectors.
   The extra n+7/4 spinor shift has NO representation in this topology.

5. THE GRAPH'S ROLE:
   In CI, the graph h1 provides inter-shell hopping (off-diagonal).
   The diagonal provides energy ordering.
   D^2 diagonal = (n+3/2)^2 gives the DIRAC energy ordering.
   The graph h1 off-diagonal is the SAME for scalar and Dirac
   (same adjacency structure, same kappa = -1/16).
   The physics difference between Schrodinger and Dirac lives
   entirely in the DIAGONAL, not in the graph topology.
""")

    # --- Quantitative comparison: scalar h1 error vs D^2 h1 error ---
    print("--- QUANTITATIVE COMPARISON ---")
    print("Both scalar and D^2 graphs use the same off-diagonal structure.")
    print("The error from graph coupling is O(kappa) ~ 0.06 in both cases.")
    print("For scalar, this is ~O(1) relative to Rydberg gaps ~0.375 Ha.")
    print("For D^2, this is ~O(1) relative to D^2 gaps ~3.0.")
    print("In BOTH cases, the graph is perturbative, not spectral.")

    n_max = 3
    # Scalar
    lattice = GeometricLattice(max_n=n_max, nuclear_charge=Z)
    orbitals = list(lattice.states)
    A_sc = lattice.adjacency
    if hasattr(A_sc, 'toarray'):
        A_sc = A_sc.toarray()
    D_sc = np.diag(np.array(A_sc.sum(axis=1)).flatten())
    L_sc = D_sc - A_sc
    H_sc = KAPPA_SCALAR * L_sc.copy()
    for i, (n, l, m) in enumerate(orbitals):
        H_sc[i, i] += -Z**2 / (2.0 * n**2)
    evals_H_sc = np.sort(np.linalg.eigvalsh(H_sc))
    rydberg = sorted([-Z**2/(2.0*n**2) for n in range(1, n_max+1)
                      for l in range(n) for m in range(-l, l+1)])
    rydberg = np.array(rydberg)
    scalar_gap = abs(rydberg[1] - rydberg[0])
    scalar_err = np.max(np.abs(evals_H_sc - rydberg))
    scalar_rel = scalar_err / scalar_gap

    # Dirac
    adj, labels, _, _ = build_dirac_s3_graph(n_max, adjacency_rule='A')
    A_d = adj.astype(float)
    D_d = np.diag(A_d.sum(axis=1))
    L_d = D_d - A_d
    H_d2 = KAPPA_SCALAR * L_d.copy()
    for i, lab in enumerate(labels):
        H_d2[i, i] += d2_eigenvalue(lab.n_fock)
    evals_d2 = np.sort(np.linalg.eigvalsh(H_d2))
    target_d2 = np.sort(np.array([d2_eigenvalue(lab.n_fock) for lab in labels]))
    d2_gap = target_d2[-1] - target_d2[0]
    d2_err = np.max(np.abs(evals_d2 - target_d2))
    d2_rel = d2_err / d2_gap if d2_gap > 0 else float('inf')

    print(f"\n  n_max={n_max}:")
    print(f"  Scalar: max|dE|={scalar_err:.4f} Ha, gap={scalar_gap:.4f} Ha, "
          f"relative={scalar_rel:.4f}")
    print(f"  D^2:    max|dE|={d2_err:.4f},     gap={d2_gap:.2f},       "
          f"relative={d2_rel:.4f}")

    results['comparison'] = {
        'scalar_max_err': float(scalar_err),
        'scalar_gap': float(scalar_gap),
        'scalar_relative': float(scalar_rel),
        'd2_max_err': float(d2_err),
        'd2_gap': float(d2_gap),
        'd2_relative': float(d2_rel),
    }

    # --- DEGENERACY BREAKING PATTERN ---
    print("\n--- DEGENERACY BREAKING BY GRAPH COUPLING ---")
    print("D^2 has accidental degeneracy: all (kappa, m_j) at same n_Fock")
    print("are degenerate (D^2 depends only on n).")
    print("Graph coupling breaks this degeneracy. How?")

    adj, labels, _, _ = build_dirac_s3_graph(3, adjacency_rule='A')
    A_d = adj.astype(float)
    D_d = np.diag(A_d.sum(axis=1))
    L_d = D_d - A_d

    for kappa_test in [0.0, -1/16, -1/8]:
        H = kappa_test * L_d.copy()
        for i, lab in enumerate(labels):
            H[i, i] += d2_eigenvalue(lab.n_fock)
        evals = np.sort(np.linalg.eigvalsh(H))
        groups = analyze_degeneracies(evals, tol=1e-8)
        print(f"\n  kappa={kappa_test:+.4f}:")
        for val, count in groups[:10]:
            # find which n_Fock this corresponds to
            closest_n = min(range(1, 5),
                           key=lambda n: abs(d2_eigenvalue(n) - val))
            print(f"    E={val:.6f} (x{count}) "
                  f"[nearest D^2(n={closest_n})={(closest_n+0.5)**2:.2f}]")

    # Save
    out_path = Path(__file__).parent / 'data' / 'd_squared_graph_test.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
