"""
α-SA Probe: Spectral Action on the Hopf Graph

Goal: Evaluate spectral action candidates on the Fock-projected S³ graph
at cutoffs m=1..5 and determine whether any naturally produces
K(m)/π = B(m) + F − Δ(m) at m=3.

Theory under test: 1/α is the total spectral content of the Hopf gauge
spectral triple, decomposed by the bundle structure into
matter (base) + gauge (fiber) − fermion boundary.
"""

import numpy as np
import json
from collections import defaultdict

# ============================================================
# Part A: Build the S³ Fock graph
# ============================================================

def build_fock_graph(m_max, edge_rule='dipole'):
    """Build the S³ Fock graph at cutoff m_max.

    Nodes: (n, l, m_l) with n=1..m_max, l=0..n-1, m_l=-l..l

    Edge rules:
    - 'fock_raise': Δn=+1, Δl=+1, |Δm|≤1 (Fock raising only)
    - 'dipole': Δn=±1, Δl=±1, |Δm|≤1 (full dipole transitions)
    - 'fock_raise_m0': Δn=+1, Δl=+1, Δm=0 (m-preserving Fock)
    """
    nodes = []
    node_index = {}
    for n in range(1, m_max + 1):
        for l in range(n):
            for ml in range(-l, l + 1):
                idx = len(nodes)
                nodes.append((n, l, ml))
                node_index[(n, l, ml)] = idx

    V = len(nodes)
    edges = []

    for i, (n1, l1, m1) in enumerate(nodes):
        for j, (n2, l2, m2) in enumerate(nodes):
            if j <= i:
                continue
            dn = n2 - n1
            dl = l2 - l1
            dm = m2 - m1

            if edge_rule == 'fock_raise':
                if dn == 1 and dl == 1 and abs(dm) <= 1:
                    edges.append((i, j))
            elif edge_rule == 'dipole':
                if abs(dn) == 1 and abs(dl) == 1 and abs(dm) <= 1:
                    edges.append((i, j))
            elif edge_rule == 'fock_raise_m0':
                if dn == 1 and dl == 1 and dm == 0:
                    edges.append((i, j))

    return nodes, edges, V


def build_incidence_matrix(V, edges):
    E = len(edges)
    B = np.zeros((V, E))
    for e, (i, j) in enumerate(edges):
        B[i, e] = 1.0
        B[j, e] = -1.0
    return B


# ============================================================
# Part B: Reference values K(m)
# ============================================================

def K_ingredients(m):
    B = sum(n**2 * (n**2 - 1) // 2 for n in range(1, m + 1))
    V = sum(n**2 for n in range(1, m + 1))
    F = np.pi**2 / 6
    Delta = 1.0 / (2 * (m + 1) * (m + 2))
    K_over_pi = B + F - Delta
    K = np.pi * K_over_pi
    return {'m': m, 'B': B, 'V': V, 'F': F, 'Delta': Delta,
            'K_over_pi': K_over_pi, 'K': K, 'B_over_V': B / V if V > 0 else 0}


# ============================================================
# Part C: Hodge decomposition and spectral invariants
# ============================================================

def hodge_spectra(B_inc):
    L0 = B_inc @ B_inc.T
    L1 = B_inc.T @ B_inc
    evals_L0 = np.sort(np.linalg.eigvalsh(L0))
    evals_L1 = np.sort(np.linalg.eigvalsh(L1))
    return L0, L1, evals_L0, evals_L1


def spectral_invariants(evals, label=""):
    nz = [x for x in evals if x > 1e-10]
    inv = {}
    inv[f'{label}_trace'] = sum(evals)
    inv[f'{label}_n_zero'] = len(evals) - len(nz)
    inv[f'{label}_n_nonzero'] = len(nz)
    if nz:
        inv[f'{label}_logdet'] = sum(np.log(x) for x in nz)
        inv[f'{label}_det'] = np.exp(inv[f'{label}_logdet'])
        for s in [0.5, 1.0, 1.5, 2.0]:
            inv[f'{label}_zeta({s})'] = sum(x**(-s) for x in nz)
        for t in [0.5, 1.0, np.pi/2, np.pi]:
            inv[f'{label}_heat(t={t:.3f})'] = sum(np.exp(-t*x) for x in evals)
    return inv


# ============================================================
# Part D: Hopf quotient (S² base graph)
# ============================================================

def build_hopf_quotient(m_max, nodes, edges):
    sectors = []
    sector_map = {}
    for n in range(1, m_max + 1):
        for l in range(n):
            idx = len(sectors)
            sectors.append((n, l))
            sector_map[(n, l)] = idx

    V_base = len(sectors)
    fiber_sizes = [2*l+1 for (n, l) in sectors]
    casimir = [l*(l+1) for (n, l) in sectors]

    q_edges_set = set()
    for i, j in edges:
        n1, l1, _ = nodes[i]
        n2, l2, _ = nodes[j]
        s1 = sector_map[(n1, l1)]
        s2 = sector_map[(n2, l2)]
        if s1 != s2:
            q_edges_set.add((min(s1, s2), max(s1, s2)))

    q_edges = sorted(q_edges_set)
    E_base = len(q_edges)

    B_base = np.zeros((V_base, E_base)) if E_base > 0 else np.zeros((V_base, 0))
    for e, (i, j) in enumerate(q_edges):
        B_base[i, e] = 1.0
        B_base[j, e] = -1.0

    if E_base > 0:
        L0_b = B_base @ B_base.T
        L1_b = B_base.T @ B_base
        ev_L0 = np.sort(np.linalg.eigvalsh(L0_b))
        ev_L1 = np.sort(np.linalg.eigvalsh(L1_b))
    else:
        L0_b = np.zeros((V_base, V_base))
        ev_L0 = np.zeros(V_base)
        L1_b = np.array([])
        ev_L1 = np.array([])

    return {
        'sectors': sectors, 'V_base': V_base, 'E_base': E_base,
        'fiber_sizes': fiber_sizes, 'casimir': casimir,
        'q_edges': q_edges,
        'evals_L0_base': ev_L0.tolist(),
        'evals_L1_base': ev_L1.tolist(),
    }


# ============================================================
# Part E: Physical Dirac / scalar spectral data
# ============================================================

def physical_spectral_data(m):
    """Camporesi-Higuchi Dirac and scalar Laplacian on S³ at cutoff m."""

    # Scalar Laplacian: eigenvalues n²-1, degeneracy n², n=1..m (Fock)
    scalar_evals = []
    for n in range(1, m + 1):
        scalar_evals.extend([(n**2 - 1)] * n**2)

    # Dirac: eigenvalues ±(n+3/2), degeneracy 2(n+1)(n+2), n_CH=0..m-1
    # (m Fock shells = m CH levels: 0..m-1)
    dirac_evals_sq = []  # D² eigenvalues
    for n_ch in range(m):
        lam = n_ch + 1.5  # |λ| = n_CH + 3/2
        g = 2 * (n_ch + 1) * (n_ch + 2)  # single chirality
        dirac_evals_sq.extend([lam**2] * (2 * g))  # both chiralities

    # Dirac boundary (the NEXT level, not included)
    g_boundary = 2 * (m + 1) * (m + 2)

    return {
        'scalar_evals': sorted(scalar_evals),
        'n_scalar': len(scalar_evals),
        'dirac_evals_sq': sorted(dirac_evals_sq),
        'n_dirac': len(dirac_evals_sq),
        'dirac_boundary_g': g_boundary,
    }


# ============================================================
# Part F: Candidate action functionals
# ============================================================

def test_candidates(m, graph_data, quotient, phys, K_data):
    """Test spectral action candidates against K(m)/π."""
    target = K_data['K_over_pi']
    candidates = {}

    evals_L0 = graph_data['evals_L0']
    evals_L1 = graph_data['evals_L1']
    V = graph_data['V']
    E = graph_data['E']
    nodes = graph_data['nodes']

    nz_L0 = [x for x in evals_L0 if x > 1e-10]
    nz_L1 = [x for x in evals_L1 if x > 1e-10]

    # --- Graph Laplacian candidates ---

    # C1: Basic traces
    candidates['Tr(L0)'] = sum(evals_L0)
    candidates['Tr(L1)'] = sum(evals_L1)
    candidates['Str = Tr(L0)-Tr(L1)'] = sum(evals_L0) - sum(evals_L1)

    # C2: Spectral zetas (supertrace)
    if nz_L0 and nz_L1:
        for s in [0.5, 1.0, 1.5, 2.0]:
            z0 = sum(x**(-s) for x in nz_L0)
            z1 = sum(x**(-s) for x in nz_L1)
            candidates[f'ζ_L0({s})'] = z0
            candidates[f'ζ_L1({s})'] = z1
            candidates[f'Str_ζ({s})'] = z0 - z1

    # C3: Heat kernel supertrace
    if len(evals_L0) > 0:
        for t in [0.5, 1.0, np.pi/2, np.pi, 2*np.pi]:
            hk0 = sum(np.exp(-t*x) for x in evals_L0)
            hk1 = sum(np.exp(-t*x) for x in evals_L1) if len(evals_L1) > 0 else 0
            candidates[f'Str_heat(t={t:.3f})'] = hk0 - hk1

    # C4: Log determinant supertrace
    if nz_L0 and nz_L1:
        ld0 = sum(np.log(x) for x in nz_L0)
        ld1 = sum(np.log(x) for x in nz_L1)
        candidates['Str_logdet'] = ld0 - ld1

    # C5: Betti numbers and Euler char
    beta0 = len(evals_L0) - len(nz_L0)
    beta1 = len(evals_L1) - len(nz_L1) if len(evals_L1) > 0 else 0
    candidates['β₀'] = beta0
    candidates['β₁'] = beta1
    candidates['χ = β₀-β₁'] = beta0 - beta1

    # --- Casimir-weighted candidates ---
    casimir_node = np.array([l*(l+1) for (n,l,ml) in nodes])
    fiber_node = np.array([2*l+1 for (n,l,ml) in nodes])

    # C6: Casimir sum = B(m)
    candidates['Σ(2l+1)l(l+1) = B'] = float(sum(fiber_node * casimir_node))

    # C7: Casimir-weighted graph zeta Tr(C·L0^{-s})
    if nz_L0:
        L0_mat = graph_data['L0']
        eigvals, eigvecs = np.linalg.eigh(L0_mat)
        C_diag = np.diag(casimir_node.astype(float))
        FC_diag = np.diag((fiber_node * casimir_node).astype(float))

        for s in [0.5, 1.0, 1.5, 2.0]:
            val_c = 0.0
            val_fc = 0.0
            for k in range(len(eigvals)):
                if eigvals[k] > 1e-10:
                    wc = eigvecs[:, k] @ C_diag @ eigvecs[:, k]
                    wfc = eigvecs[:, k] @ FC_diag @ eigvecs[:, k]
                    val_c += wc * eigvals[k]**(-s)
                    val_fc += wfc * eigvals[k]**(-s)
            candidates[f'Tr(C·L0^(-{s}))'] = val_c
            candidates[f'Tr(FC·L0^(-{s}))'] = val_fc

    # --- Physical Dirac candidates ---

    # C8: Scalar spectral action (counting)
    candidates['N_scalar'] = phys['n_scalar']
    candidates['N_dirac'] = phys['n_dirac']
    candidates['Str_count = N_s - N_d/4'] = phys['n_scalar'] - phys['n_dirac']/4

    # C9: Scalar/Dirac spectral zeta
    sc_nz = [x for x in phys['scalar_evals'] if x > 0]
    di_nz = [x for x in phys['dirac_evals_sq'] if x > 0]

    if sc_nz:
        for s in [0.5, 1.0, 1.5, 2.0]:
            zs = sum(x**(-s) for x in sc_nz)
            candidates[f'ζ_scalar({s})'] = zs

    if di_nz:
        for s in [0.5, 1.0, 1.5, 2.0]:
            zd = sum(x**(-s) for x in di_nz)
            candidates[f'ζ_D²({s})'] = zd

    if sc_nz and di_nz:
        for s in [0.5, 1.0, 1.5, 2.0]:
            zs = sum(x**(-s) for x in sc_nz)
            zd = sum(x**(-s) for x in di_nz)
            candidates[f'Str_phys({s}) = ζ_s-ζ_D/4'] = zs - zd/4

    # C10: Mixed: Casimir sum + graph zeta + Dirac boundary
    if nz_L0:
        for s in [0.5, 1.0, 1.5, 2.0]:
            z0 = sum(x**(-s) for x in nz_L0)
            candidates[f'B + ζ_L0({s})'] = K_data['B'] + z0
            candidates[f'B + ζ_L0({s}) - Δ'] = K_data['B'] + z0 - K_data['Delta']

    # C11: Mixed with fiber spectral content
    # F_partial(m) = Σ_{n=1}^m 1/n² (partial sum of ζ(2))
    F_partial = sum(1.0/n**2 for n in range(1, m+1))
    candidates['F_partial'] = F_partial
    candidates['B + F_partial'] = K_data['B'] + F_partial
    candidates['B + F_partial - Δ'] = K_data['B'] + F_partial - K_data['Delta']
    candidates['B + F - Δ (=K/π target)'] = K_data['K_over_pi']

    # C12: Quotient Laplacian
    if quotient['E_base'] > 0:
        qev = [x for x in quotient['evals_L0_base'] if x > 1e-10]
        if qev:
            candidates['Tr(L0_base)'] = sum(quotient['evals_L0_base'])
            for s in [0.5, 1.0, 1.5, 2.0]:
                candidates[f'ζ_base({s})'] = sum(x**(-s) for x in qev)

    # C13: Per-shell Casimir c(n) = n²(n²-1)/2 as "action density"
    # Action = Σ c(n) weighted by something from the Dirac spectrum
    for n_fock in range(1, m+1):
        n_ch = n_fock - 1
        lam = n_ch + 1.5
        g_dirac = 2 * (n_ch + 1) * (n_ch + 2)
        cn = n_fock**2 * (n_fock**2 - 1) / 2

    # C14: Wilson-type: ratio of Hodge spectral content
    if nz_L0 and nz_L1:
        candidates['ζ_L0(1)/ζ_L1(1)'] = sum(1/x for x in nz_L0) / sum(1/x for x in nz_L1)
        candidates['det(L0_nz)/det(L1_nz)'] = np.exp(
            sum(np.log(x) for x in nz_L0) - sum(np.log(x) for x in nz_L1))

    # C15: Physical Dirac mode-count decomposition via EM
    # S_D(m) = 4(m+1)(m+2)(m+3)/3 (total Dirac modes up to n_CH=m)
    S_D_total = 4 * (m) * (m + 1) * (m + 2) / 3  # using n_CH = 0..m-1
    candidates['S_D_total'] = S_D_total
    candidates['S_D_total / Δ⁻¹'] = S_D_total / (2*(m+1)*(m+2)) if m > 0 else 0

    return candidates


# ============================================================
# Main
# ============================================================

def main():
    alpha_inv = 137.035999084
    K_target = alpha_inv
    Kpi_target = K_target / np.pi

    print("=" * 78)
    print("α-SA PROBE: Spectral Action on the Hopf Graph")
    print("=" * 78)
    print(f"Target: K = 1/α = {K_target:.10f}")
    print(f"Target: K/π = {Kpi_target:.10f}")

    results = {}

    for edge_rule in ['fock_raise', 'dipole']:
        print(f"\n{'#'*78}")
        print(f"# Edge rule: {edge_rule}")
        print(f"{'#'*78}")

        results[edge_rule] = {}

        for m in range(1, 6):
            print(f"\n{'='*60}")
            print(f"  m = {m}  (edge rule: {edge_rule})")
            print(f"{'='*60}")

            nodes, edges, V = build_fock_graph(m, edge_rule)
            E = len(edges)
            K_data = K_ingredients(m)

            print(f"  V={V}, E={E}, B={K_data['B']}, B/V={K_data['B_over_V']:.4f}")
            print(f"  K(m)/π = {K_data['K_over_pi']:.6f}")

            graph_data = {'nodes': nodes, 'V': V, 'E': E}

            if E > 0:
                B_inc = build_incidence_matrix(V, edges)
                L0, L1, ev0, ev1 = hodge_spectra(B_inc)
                graph_data['L0'] = L0
                graph_data['evals_L0'] = ev0
                graph_data['evals_L1'] = ev1

                beta0 = sum(1 for x in ev0 if x < 1e-10)
                beta1 = sum(1 for x in ev1 if x < 1e-10)
                print(f"  L0 spec: {np.round(ev0, 3).tolist()}")
                print(f"  L1 spec: {np.round(ev1, 3).tolist()}")
                print(f"  β₀={beta0}, β₁={beta1}, χ={beta0-beta1}")
            else:
                graph_data['L0'] = np.zeros((V,V))
                graph_data['evals_L0'] = np.zeros(V)
                graph_data['evals_L1'] = np.array([])

            quotient = build_hopf_quotient(m, nodes, edges)
            print(f"  S² quotient: V={quotient['V_base']}, E={quotient['E_base']}")
            print(f"  Sectors: {quotient['sectors']}")
            print(f"  Fibers: {quotient['fiber_sizes']}, Casimir: {quotient['casimir']}")
            if quotient['evals_L0_base']:
                print(f"  L0_base spec: {np.round(quotient['evals_L0_base'], 3).tolist()}")

            phys = physical_spectral_data(m)

            candidates = test_candidates(m, graph_data, quotient, phys, K_data)

            # Sort by proximity to K(m)/π
            target = K_data['K_over_pi']
            scored = []
            for name, val in candidates.items():
                if isinstance(val, (int, float)) and np.isfinite(val) and abs(val) > 1e-15:
                    rel = abs(val / target - 1.0) if abs(target) > 1e-15 else float('inf')
                    scored.append((name, val, rel))
            scored.sort(key=lambda x: x[2])

            print(f"\n  Top 20 closest to K(m)/π = {target:.6f}:")
            for name, val, rel in scored[:20]:
                hit = " <<<" if rel < 0.001 else (" <~" if rel < 0.05 else "")
                print(f"    {name:40s} = {val:14.6f}  rel_err={rel:.6f}{hit}")

            results[edge_rule][m] = {
                'V': V, 'E': E,
                'K_data': {k: float(v) for k, v in K_data.items()},
                'evals_L0': graph_data['evals_L0'].tolist() if isinstance(graph_data['evals_L0'], np.ndarray) else [],
                'evals_L1': graph_data['evals_L1'].tolist() if isinstance(graph_data['evals_L1'], np.ndarray) else [],
                'quotient_summary': {
                    'V_base': quotient['V_base'], 'E_base': quotient['E_base'],
                    'sectors': quotient['sectors'],
                    'evals_L0_base': quotient['evals_L0_base'],
                },
                'top10': [(n, float(v), float(r)) for n, v, r in scored[:10]],
            }

    # ============================================================
    # SYNTHESIS: What works at m=3?
    # ============================================================
    print(f"\n{'='*78}")
    print("SYNTHESIS: Which candidates hit K(3)/π = 43.619900 across edge rules?")
    print(f"{'='*78}")

    for rule in ['fock_raise', 'dipole']:
        if 3 in results[rule]:
            print(f"\n  Edge rule: {rule}")
            for name, val, rel in results[rule][3]['top10']:
                hit = " <<<" if rel < 0.001 else (" <~" if rel < 0.05 else "")
                print(f"    {name:40s} = {val:14.6f}  rel_err={rel:.6f}{hit}")

    # ============================================================
    # CRITICAL TEST: Does any graph invariant track K(m)/π across m?
    # ============================================================
    print(f"\n{'='*78}")
    print("TRACKING TEST: Do any candidates track K(m)/π across m=2..5?")
    print(f"{'='*78}")

    for rule in ['fock_raise', 'dipole']:
        print(f"\n  Edge rule: {rule}")
        # Collect candidate names that appear at all m
        all_names = set()
        for m in range(2, 6):
            if m in results[rule]:
                for name, _, _ in results[rule][m]['top10']:
                    all_names.add(name)

        # For each candidate, check if it tracks K(m)/π
        tracking = []
        for name in all_names:
            ratios = []
            for m in range(2, 6):
                if m in results[rule]:
                    for n, v, r in results[rule][m]['top10']:
                        if n == name:
                            ratios.append((m, v / results[rule][m]['K_data']['K_over_pi']))
            if len(ratios) >= 3:
                ratio_vals = [r for _, r in ratios]
                cv = np.std(ratio_vals) / np.mean(ratio_vals) if np.mean(ratio_vals) != 0 else float('inf')
                tracking.append((name, ratios, cv))

        tracking.sort(key=lambda x: x[2])
        for name, ratios, cv in tracking[:5]:
            print(f"    {name:40s} CV={cv:.4f}  ratios={[(m, f'{r:.4f}') for m, r in ratios]}")

    out_path = 'debug/data/alpha_sa_hopf_action.json'
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved to {out_path}")


if __name__ == '__main__':
    main()
