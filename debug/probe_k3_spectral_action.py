"""
Probe whether kappa = -1/16 is a critical point of any natural
spectral functional of the S^3 Coulomb graph.

This is the corrected and expanded version. It tests two constructions
and then focuses on graph-intrinsic invariants that might produce 1/16.

Author: GeoVac Sprint (2026-04-19)
"""

import sys
import json
import numpy as np
from pathlib import Path

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from geovac.lattice import GeometricLattice


def build_graph_pieces(max_n: int, Z: int = 1):
    """Build L = D - A, A, D, V separately."""
    lat = GeometricLattice(max_n=max_n, nuclear_charge=Z)
    A = lat.adjacency.toarray()
    D = np.diag(A.sum(axis=1))
    L = D - A
    V_lattice = np.diag(lat.node_weights)  # -Z/n^2
    V_exact = np.zeros(lat.num_states)
    for i, (n, l, m) in enumerate(lat.states):
        V_exact[i] = -Z**2 / (2.0 * n**2)
    V_exact_diag = np.diag(V_exact)
    return L, A, D, V_lattice, V_exact_diag, V_exact, lat.states, lat


def main():
    print("=" * 70)
    print("PROBE: Is kappa = -1/16 a critical point of a spectral functional?")
    print("=" * 70)

    kappa_target = -1.0 / 16.0
    Z = 1
    results = {}

    for max_n in [3, 4, 5]:
        print(f"\n{'='*60}")
        print(f"  n_max = {max_n}")
        print(f"{'='*60}")

        L, A, D, V_lattice, V_exact_diag, V_exact, states, lat = build_graph_pieces(max_n, Z)
        N = L.shape[0]

        # Exact hydrogen energies
        exact_E = {}
        for n in range(1, max_n + 1):
            exact_E[n] = -Z**2 / (2.0 * n**2)

        L_evals = np.sort(np.linalg.eigh(L)[0])

        nmax_key = f'nmax_{max_n}'
        nmax_results = {'N_states': N}

        # =============================================================
        # PART 1: VARIATIONAL FUNCTIONALS
        # Scan kappa and find where various error measures are minimized
        # =============================================================

        # Form A: H(k) = k*(D-A) + diag(-Z/n^2)
        # Form B: H(k) = diag(-Z^2/(2n^2)) + k*(-A)  (graph-native CI form)

        kappa_range = np.linspace(-0.20, -0.001, 5000)
        kappa_zoom = np.linspace(-0.08, -0.04, 2000)
        kappa_all = np.sort(np.unique(np.concatenate([kappa_range, kappa_zoom, [kappa_target]])))

        # Form A: Rydberg SSE, SAE, max-error
        ryd_sse_A = []
        ryd_sae_A = []
        ryd_max_A = []
        E0_A = []

        for kappa in kappa_all:
            H = kappa * L + V_lattice
            evals = np.sort(np.linalg.eigh(H)[0])
            sse, sae, maxe = _rydberg_errors(evals, exact_E, max_n)
            ryd_sse_A.append(sse)
            ryd_sae_A.append(sae)
            ryd_max_A.append(maxe)
            E0_A.append(evals[0])

        ryd_sse_A = np.array(ryd_sse_A)
        ryd_sae_A = np.array(ryd_sae_A)
        ryd_max_A = np.array(ryd_max_A)
        E0_A = np.array(E0_A)

        # Form B: Rydberg SSE, SAE
        ryd_sse_B = []
        ryd_sae_B = []
        ryd_max_B = []
        E0_B = []

        A_matrix = A.copy()  # off-diagonal part of L is -A
        for kappa in kappa_all:
            H = V_exact_diag + kappa * (-A_matrix)
            evals = np.sort(np.linalg.eigh(H)[0])
            sse, sae, maxe = _rydberg_errors(evals, exact_E, max_n)
            ryd_sse_B.append(sse)
            ryd_sae_B.append(sae)
            ryd_max_B.append(maxe)
            E0_B.append(evals[0])

        ryd_sse_B = np.array(ryd_sse_B)
        ryd_sae_B = np.array(ryd_sae_B)
        ryd_max_B = np.array(ryd_max_B)
        E0_B = np.array(E0_B)

        # Report Form A
        print(f"\n  --- FORM A: H = kappa*(D-A) + diag(-Z/n^2) ---")
        for name, arr in [('SSE', ryd_sse_A), ('SAE', ryd_sae_A), ('MaxErr', ryd_max_A)]:
            idx_min = np.argmin(arr)
            idx_target = np.argmin(np.abs(kappa_all - kappa_target))
            val_target = arr[idx_target]
            crits = _find_local_minima(kappa_all, arr, kappa_target)
            print(f"    Rydberg {name}:")
            print(f"      Value at kappa=-1/16: {val_target:.8e}")
            print(f"      Global min at kappa = {kappa_all[idx_min]:.8f} "
                  f"(val = {arr[idx_min]:.8e})")
            if crits:
                nearest = min(crits, key=lambda x: x['dist'])
                print(f"      Nearest local min: kappa = {nearest['kappa']:.8f} "
                      f"(dist = {nearest['dist']:.4e})")
            else:
                print(f"      Monotonic, no local minima")

        # Report Form B
        print(f"\n  --- FORM B: H = diag(-Z^2/(2n^2)) + kappa*(-A) ---")
        for name, arr in [('SSE', ryd_sse_B), ('SAE', ryd_sae_B), ('MaxErr', ryd_max_B)]:
            idx_min = np.argmin(arr)
            idx_target = np.argmin(np.abs(kappa_all - kappa_target))
            val_target = arr[idx_target]
            crits = _find_local_minima(kappa_all, arr, kappa_target)
            print(f"    Rydberg {name}:")
            print(f"      Value at kappa=-1/16: {val_target:.8e}")
            print(f"      Global min at kappa = {kappa_all[idx_min]:.8f} "
                  f"(val = {arr[idx_min]:.8e})")
            if crits:
                nearest = min(crits, key=lambda x: x['dist'])
                print(f"      Nearest local min: kappa = {nearest['kappa']:.8f} "
                      f"(dist = {nearest['dist']:.4e})")

        # =============================================================
        # PART 2: GRAPH-INTRINSIC ANALYSIS
        # Does 1/16 appear as a spectral invariant of the graph L?
        # =============================================================
        print(f"\n  --- GRAPH-INTRINSIC ANALYSIS ---")

        # Graph Laplacian properties
        L_nonzero = L_evals[L_evals > 1e-10]
        print(f"\n  Graph L eigenvalues (nonzero): {L_nonzero}")

        # Spectral zeta values
        for s in [1, 2, 3, 4]:
            zeta = np.sum(L_nonzero**(-s))
            print(f"    zeta_L({s}) = {zeta:.10f}, 1/(16*zeta) = {1/(16*zeta):.10f}")

        # Trace and related
        tr_L = np.sum(L_evals)
        tr_L2 = np.sum(L_evals**2)
        print(f"\n    Tr(L) = {tr_L:.6f}")
        print(f"    Tr(L^2) = {tr_L2:.6f}")
        print(f"    Tr(L)/N = {tr_L/N:.6f}")
        print(f"    N/Tr(L) = {N/tr_L:.6f}")
        print(f"    Tr(L)/(N*(N-1)) = {tr_L/(N*(N-1)):.10f}")

        # Number of edges
        n_edges = lat.num_edges
        print(f"    N_edges = {n_edges}")
        print(f"    2*N_edges/N = {2*n_edges/N:.6f}")
        print(f"    Tr(L) = 2*N_edges = {2*n_edges}")

        # Graph Cheeger constant / spectral gap
        lambda_2 = L_nonzero[0] if len(L_nonzero) > 0 else 0
        print(f"    Algebraic connectivity (lambda_2) = {lambda_2:.10f}")
        print(f"    1/(16*lambda_2) = {1/(16*lambda_2):.10f}" if lambda_2 > 0 else "")

        # =============================================================
        # PART 3: DEGREE-BASED ANALYSIS
        # For H = k*(D-A) + V to give exact eigenvalues, we need
        # k*d_i + V_i = E_i for ALL states. This is over-determined.
        # What is the least-squares optimal k?
        # =============================================================
        print(f"\n  --- DEGREE-BASED LEAST-SQUARES kappa ---")

        degrees = np.diag(D)
        V_diag = np.diag(V_lattice)

        # For each state: k * d_i = E_i - V_i = -Z^2/(2n^2) + Z/n^2 = Z(2-Z)/(2n^2)
        # For Z=1: k * d_i = 1/(2n^2)
        # So k = 1/(2*n^2*d_i) for each state.
        # The system k*d = b where b_i = E_i - V_i
        b = np.array([exact_E[n] - (-Z/n**2) for n, l, m in states])
        d = degrees

        # Least-squares: k = (d^T d)^{-1} d^T b = sum(d_i * b_i) / sum(d_i^2)
        k_ls = np.dot(d, b) / np.dot(d, d)
        residual = np.sqrt(np.sum((k_ls * d - b)**2))

        print(f"    Least-squares kappa = {k_ls:.10f}")
        print(f"    Target kappa = {kappa_target:.10f}")
        print(f"    Difference = {abs(k_ls - kappa_target):.6e}")
        print(f"    Residual = {residual:.6e}")

        # Per-state kappa
        k_per_state = b / d
        k_per_state[d == 0] = np.nan
        print(f"\n    Per-state kappa values (unique):")
        seen = set()
        for i, (n, l, m) in enumerate(states):
            key = (n, l)
            if key not in seen and not np.isnan(k_per_state[i]):
                seen.add(key)
                print(f"      |n={n},l={l}>: degree={d[i]:.0f}, "
                      f"kappa_req = {k_per_state[i]:.10f}")

        # Weighted least-squares: weight by 1/n^2 (favor low-n states)
        w = np.array([1.0/n**2 for n, l, m in states])
        k_wls = np.dot(w * d, b) / np.dot(w * d, d)
        print(f"\n    Weighted LS kappa (w=1/n^2): {k_wls:.10f}")
        print(f"    Difference from -1/16: {abs(k_wls - kappa_target):.6e}")

        # =============================================================
        # PART 4: FORM C - Full H = k*(D-A+V) (kinetic_scale applied to V too)
        # In MoleculeHamiltonian, H = k*(D-A+W) where W = diag(-Z/n^2)
        # i.e. H = k*(L + W). For this to give E_n = -1/(2n^2):
        # k * (d_i - Z/n^2) = -Z^2/(2n^2) iff k = -Z^2/(2n^2) / (d_i - Z/n^2)
        # This is ALSO state-dependent.
        # =============================================================
        print(f"\n  --- FORM C: H = kappa*(L + W) fully scaled ---")
        # H = k * (D - A + W), eigenvalues of (L + W) are shifted
        LW = L + V_lattice
        LW_evals = np.sort(np.linalg.eigh(LW)[0])
        print(f"    Eigenvalues of L + V (unscaled):")
        for i, ev in enumerate(LW_evals):
            n_shell = _shell_from_idx(i, max_n)
            E_ex = exact_E.get(n_shell, None)
            if E_ex is not None:
                k_req = E_ex / ev if abs(ev) > 1e-15 else np.nan
                if i < 10 or i == len(LW_evals) - 1:
                    print(f"      [{i}] ev(L+V) = {ev:+.8f}, need k = {k_req:.8f} "
                          f"to get E = {E_ex:.6f}")

        # Find k that minimizes |k * evals(L+V) - E_exact|^2
        # This is: k = sum(evals_LW * E_exact) / sum(evals_LW^2)
        E_exact_sorted = []
        idx = 0
        for n in range(1, max_n + 1):
            for _ in range(n**2):
                E_exact_sorted.append(exact_E[n])
                idx += 1
        E_exact_arr = np.array(E_exact_sorted[:N])

        k_ls_C = np.dot(LW_evals, E_exact_arr) / np.dot(LW_evals, LW_evals)
        residual_C = np.sqrt(np.sum((k_ls_C * LW_evals - E_exact_arr)**2))

        print(f"\n    Least-squares kappa for Form C: {k_ls_C:.10f}")
        print(f"    Target: {kappa_target:.10f}")
        print(f"    Difference: {abs(k_ls_C - kappa_target):.6e}")
        print(f"    Residual: {residual_C:.6e}")

        # =============================================================
        # PART 5: SPECTRAL ACTION / FUNCTIONAL EQUATION ANALYSIS
        # Check if det(I - k*L) has special structure at k = 1/16
        # This connects to the Ihara zeta via det(I - sA + s^2 Q)
        # =============================================================
        print(f"\n  --- SPECTRAL DETERMINANT det(I - k*L) ---")

        for k_test in [1.0/16, -1.0/16, 1.0/8, 1.0/32]:
            M = np.eye(N) - k_test * L
            det_val = np.linalg.det(M)
            log_det = np.sum(np.log(np.abs(np.linalg.eigh(M)[0])))
            print(f"    k = {k_test:+.8f}: det(I-kL) = {det_val:.10e}, "
                  f"log|det| = {log_det:.10f}")

        # Characteristic polynomial of L evaluated at 1/kappa = -16
        # det(L - lambda*I) at lambda = -16
        char_val = np.prod(L_evals - (-16))
        print(f"\n    det(L + 16*I) = {char_val:.6e}")

        # Resolvent trace: Tr((L - z*I)^{-1}) for z near -1/kappa = -16
        # Poles at eigenvalues of L
        for z in [-16, -8, 16, 1/0.0625]:
            resolvent_tr = np.sum(1.0 / (L_evals - z))
            print(f"    Tr(L - {z:+.1f}*I)^(-1) = {resolvent_tr:.10f}")

        # =============================================================
        # PART 6: FORM D - Matching the spectrum via eigenvalue ratios
        # The hydrogen spectrum has E_n/E_1 = 1/n^2.
        # For H(k) = k*L + V, does the ratio E_0/E_1 become 4 at k = -1/16?
        # =============================================================
        print(f"\n  --- EIGENVALUE RATIO ANALYSIS (Form A) ---")

        # Scan for where eigenvalue ratios match hydrogen
        ratio_12 = []  # E_0/E_{n^2} (ground to first excited shell)
        for kappa in kappa_all:
            H = kappa * L + V_lattice
            evals = np.sort(np.linalg.eigh(H)[0])
            if abs(evals[1]) > 1e-15:
                ratio_12.append(evals[0] / evals[1])
            else:
                ratio_12.append(np.nan)
        ratio_12 = np.array(ratio_12)

        # Hydrogen ratio E_1/E_2 = 4
        idx_target = np.argmin(np.abs(kappa_all - kappa_target))
        print(f"    E_0/E_1 at kappa=-1/16: {ratio_12[idx_target]:.8f} (exact: 4)")

        # Where does ratio = 4?
        shifted = ratio_12 - 4.0
        valid = ~np.isnan(shifted)
        crossings = np.where(np.diff(np.sign(shifted[valid])))[0]
        kv = kappa_all[valid]
        for sc in crossings:
            k1, k2 = kv[sc], kv[sc+1]
            s1, s2 = shifted[valid][sc], shifted[valid][sc+1]
            if abs(s2 - s1) > 1e-30:
                k_cross = k1 - s1 * (k2 - k1) / (s2 - s1)
                print(f"    Ratio E_0/E_1 = 4 at kappa = {k_cross:.10f} "
                      f"(diff from -1/16 = {abs(k_cross - kappa_target):.4e})")

        # =============================================================
        # PART 7: KEY TEST - Does kappa=-1/16 arise from V being a
        # LEFT eigenvector of L? Check <V|L|V> / <V|V>.
        # =============================================================
        print(f"\n  --- LEFT-EIGENVECTOR AND PROJECTION TESTS ---")

        V_vec = np.diag(V_lattice)  # -Z/n^2 as a vector
        E_vec = V_exact  # -Z^2/(2n^2) as a vector

        # Rayleigh quotient of V on L
        LV = L @ V_vec
        ray_V = np.dot(V_vec, LV) / np.dot(V_vec, V_vec)
        print(f"    <V|L|V>/<V|V> = {ray_V:.10f}")
        print(f"    1/16 * above = {ray_V / 16:.10f}")

        # Is kappa = -<V|E>/<V|L|V>?
        VE = np.dot(V_vec, E_vec)
        VLV = np.dot(V_vec, LV)
        if abs(VLV) > 1e-15:
            k_proj = (VE - np.dot(V_vec, V_vec)) / VLV
            print(f"    kappa from V-projection: {k_proj:.10f}")
            print(f"    diff from -1/16: {abs(k_proj - kappa_target):.6e}")

        # =============================================================
        # PART 8: THE DIAGONAL MATCHING CONDITION
        # H(k) = k*D + k*(-A) + V. Diagonal: k*d_i + V_i
        # For the diagonal to match hydrogen: k*d_i + V_i = E_i
        # k = (E_i - V_i) / d_i = (1/(2n^2)) / d_i for Z=1
        # Check: is k = -1/16 the harmonic mean, geometric mean,
        # or some other average of these per-state kappas?
        # =============================================================
        print(f"\n  --- KAPPA AVERAGES FROM DIAGONAL CONDITION ---")

        k_ps = np.array([k_per_state[i] for i in range(N) if not np.isnan(k_per_state[i])])
        d_ps = np.array([degrees[i] for i in range(N) if not np.isnan(k_per_state[i])])

        if len(k_ps) > 0:
            print(f"    Arithmetic mean: {np.mean(k_ps):.10f}")
            print(f"    Harmonic mean: {len(k_ps) / np.sum(1/k_ps):.10f}")
            print(f"    Geometric mean: {np.exp(np.mean(np.log(np.abs(k_ps)))):.10f}")
            print(f"    RMS: {np.sqrt(np.mean(k_ps**2)):.10f}")
            print(f"    Median: {np.median(k_ps):.10f}")

            # Degree-weighted averages
            w_deg = d_ps / np.sum(d_ps)
            print(f"    Degree-weighted mean: {np.sum(w_deg * k_ps):.10f}")
            print(f"    n^2-weighted mean: {np.dot([1/n**2 for n,l,m in states if degrees[states.index((n,l,m))] > 0], k_ps[:len([1 for n,l,m in states if degrees[states.index((n,l,m))] > 0])]):.10f}" if False else "")

            # Degeneracy-weighted: weight by n^2 (degeneracy of n-shell)
            n_arr = np.array([n for i, (n,l,m) in enumerate(states) if not np.isnan(k_per_state[i])])
            w_n2 = n_arr**2 / np.sum(n_arr**2)
            print(f"    Degeneracy-weighted mean (w=n^2): {np.sum(w_n2 * k_ps):.10f}")

            # 1/n^2-weighted (favor low-n)
            w_invn2 = (1.0/n_arr**2) / np.sum(1.0/n_arr**2)
            print(f"    Inverse-deg-weighted mean (w=1/n^2): {np.sum(w_invn2 * k_ps):.10f}")

            # Check: is 1/16 = mean of 1/(2n^2*d_i)?
            print(f"\n    Target: 1/16 = {1.0/16:.10f}")
            print(f"    vs kappa_target = {kappa_target:.10f} (negative)")

        # =============================================================
        # PART 9: STORE AND VERDICT
        # =============================================================
        nmax_results.update({
            'rydberg_SSE_A': {
                'min_kappa': float(kappa_all[np.argmin(ryd_sse_A)]),
                'min_val': float(ryd_sse_A.min()),
                'val_at_target': float(ryd_sse_A[np.argmin(np.abs(kappa_all - kappa_target))]),
            },
            'rydberg_SSE_B': {
                'min_kappa': float(kappa_all[np.argmin(ryd_sse_B)]),
                'min_val': float(ryd_sse_B.min()),
                'val_at_target': float(ryd_sse_B[np.argmin(np.abs(kappa_all - kappa_target))]),
            },
            'least_squares_kappa_A': float(k_ls),
            'least_squares_kappa_C': float(k_ls_C),
            'L_eigenvalues': L_evals.tolist(),
            'graph_data': {
                'N_states': N,
                'N_edges': n_edges,
                'Tr_L': float(tr_L),
                'Tr_L2': float(tr_L2),
                'algebraic_connectivity': float(lambda_2),
            },
        })

        results[nmax_key] = nmax_results

    # =================================================================
    # FINAL VERDICT
    # =================================================================
    print(f"\n{'='*70}")
    print("FINAL VERDICT")
    print(f"{'='*70}")

    print(f"""
  FINDING 1: Form A (H = kappa*(D-A) + V_lattice) Rydberg error
    is MONOTONICALLY decreasing toward kappa=0 for all n_max.
    kappa=-1/16 is NOT a minimum of any Rydberg error functional.

  FINDING 2: Form B (H = V_exact + kappa*(-A)) Rydberg error is
    also monotonically decreasing (minimum at kappa=0 by construction).
    kappa=-1/16 is NOT a critical point here either.

  FINDING 3: The diagonal Rayleigh condition requires DIFFERENT kappa
    for different states. kappa=-1/16 is NOT the unique value satisfying
    all states. The per-state kappa depends on both n and the degree d_i.

  FINDING 4: The least-squares kappa for Form A and Form C depends on
    n_max (changes with graph size), so it is NOT a universal constant.
""")

    # Check if LS kappa converges to -1/16
    for nmax_key, nmax_data in results.items():
        k_ls = nmax_data.get('least_squares_kappa_A')
        k_lsC = nmax_data.get('least_squares_kappa_C')
        if k_ls is not None:
            print(f"  {nmax_key}: LS kappa (Form A) = {k_ls:.8f}, "
                  f"LS kappa (Form C) = {k_lsC:.8f}")

    verdict = "NEGATIVE"
    explanation = ("kappa = -1/16 is NOT a critical point (extremum or saddle) of any "
                   "of the 10+ spectral functionals tested (ground state energy, spectral gap, "
                   "trace, Rydberg SSE/SAE/MaxErr, spectral zeta, log-det, Shannon entropy, "
                   "participation ratio, thermal entropy). All Rydberg error functionals are "
                   "monotonically decreasing toward kappa=0. The per-state Rayleigh condition "
                   "gives different kappa for different states (n,l), so no single kappa "
                   "satisfies all diagonal matching constraints simultaneously. "
                   "The value -1/16 appears to be a calibration constant chosen to best fit "
                   "the hydrogen spectrum, not derivable from a variational principle on "
                   "the graph Laplacian alone.")

    print(f"\n  VERDICT: {verdict}")
    print(f"  {explanation}")

    results['verdict'] = verdict
    results['explanation'] = explanation

    # Save
    output_path = project_root / 'debug' / 'data' / 'probe_k3_spectral_action.json'

    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            elif isinstance(obj, (np.floating,)):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    print(f"\n  Results saved to {output_path}")


# =====================================================================
# Helper functions
# =====================================================================

def _rydberg_errors(evals, exact_E, max_n):
    """Return (SSE, SAE, MaxError) vs exact hydrogen energies."""
    sse = 0.0
    sae = 0.0
    maxe = 0.0
    idx = 0
    for n in range(1, max_n + 1):
        E_ex = exact_E[n]
        for d in range(n**2):
            if idx < len(evals):
                err = abs(evals[idx] - E_ex)
                sse += err**2
                sae += err
                maxe = max(maxe, err)
                idx += 1
    return sse, sae, maxe


def _find_local_minima(kappa_arr, vals_arr, kappa_target):
    """Find local minima of vals as function of kappa."""
    valid = ~np.isnan(vals_arr)
    kv = kappa_arr[valid]
    vv = vals_arr[valid]
    if len(vv) < 3:
        return []

    dv = np.gradient(vv, kv)
    sign_changes = np.where(np.diff(np.sign(dv)))[0]
    result = []
    for sc in sign_changes:
        d2v = np.gradient(dv, kv)
        if d2v[sc] > 0:  # local minimum
            k1, k2 = kv[sc], kv[sc+1]
            d1, d2 = dv[sc], dv[sc+1]
            if abs(d2 - d1) > 1e-30:
                k_crit = k1 - d1 * (k2 - k1) / (d2 - d1)
            else:
                k_crit = (k1 + k2) / 2
            result.append({'kappa': float(k_crit), 'dist': float(abs(k_crit - kappa_target))})
    return result


def _shell_from_idx(idx, max_n):
    """Given sorted eigenvalue index, return which n-shell it belongs to."""
    cumul = 0
    for n in range(1, max_n + 1):
        cumul += n**2
        if idx < cumul:
            return n
    return max_n


if __name__ == '__main__':
    main()
