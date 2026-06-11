#!/usr/bin/env python3
"""
Spin-Doubled Graph Test
========================

Tests whether spin can emerge from graph topology by doubling the scalar
Fock graph with a Z_2 fiber coupling that corresponds to the center
{I, -I} of SU(2) = S^3.

Hypothesis: Paper 0's packing construction has a factor of 2 from S^2
orientability.  On S^3 = SU(2), this Z_2 IS the spin structure.  If
the Z_2 is made topological (fiber edges connecting the two copies of
each node), the graph eigenvalues might produce fine-structure-like
splittings.

Four fiber-coupling prescriptions are tested:
  A  uniform:      w = eps
  B  l-dependent:  w = eps * l(l+1) / n^3      (s-states protected)
  C  kappa-aware:  w = eps * l(l+1) / (n^3 (2l+1))
  D  Hopf-geometric: w = eps * sin(chi_n)       (Fock conformal angle)

Author: GeoVac (spin-topology exploration)
Date:   April 2026
"""

import sys, os, json, itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse import lil_matrix, csr_matrix, diags
from collections import defaultdict

# ---- project imports ----
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from geovac.lattice import GeometricLattice
from geovac.hopf_bundle import fock_chi, state_to_s3, decompose_hopf

# ---- physical constants ----
KAPPA_FOCK = -1.0 / 16.0          # universal topological constant
ALPHA_PHYS = 7.2973525693e-3       # fine structure constant (CODATA 2022)


# =============================================================
# 1. Build / verify the scalar Fock graph
# =============================================================

def build_scalar_graph(n_max: int):
    """
    Build the scalar Fock graph adjacency matrix using GeometricLattice
    and return its graph Laplacian eigenvalues for verification.

    Returns
    -------
    lattice : GeometricLattice
    evals_laplacian : np.ndarray
        Sorted eigenvalues of L = D - A.
    """
    lattice = GeometricLattice(max_n=n_max)
    A = lattice.adjacency.toarray().astype(float)
    D = np.diag(A.sum(axis=1))
    L = D - A
    evals = np.linalg.eigvalsh(L)
    return lattice, evals


def verify_scalar_spectrum(lattice, evals, n_max: int):
    """
    Check that kappa * L eigenvalues reproduce -1/(n^2) Rydberg spectrum.

    The Fock graph Laplacian has eigenvalues n^2 - 1 (with degeneracy n^2)
    for n = 1 .. n_max.  So H = kappa * L gives E_n = kappa * (n^2 - 1)
    which should equal -1/(16) * (n^2 - 1).

    For hydrogen: the exact energy is E_n = -1/(2 n^2).  The graph
    approximation gives E_n = -1/16 * (n^2-1) which matches only at n=1
    exactly (E_1 = 0 from L, but the physical ground state energy is
    -0.5 Ha --- the offset is absorbed by the constant shift).

    Here we just check the Laplacian eigenvalue structure (n^2 - 1).
    """
    # expected Laplacian eigenvalues: n^2 - 1 with degeneracy n^2
    expected = {}
    for n in range(1, n_max + 1):
        expected[n * n - 1] = expected.get(n * n - 1, 0) + n * n

    # group numerical eigenvalues by rounding
    tol = 1e-8
    rounded = np.round(evals, decimals=6)
    unique, counts = np.unique(rounded, return_counts=True)

    match_ok = True
    for val, cnt in zip(unique, counts):
        ival = int(round(val))
        if ival not in expected or expected[ival] != cnt:
            match_ok = False
    return match_ok, dict(zip([int(round(v)) for v in unique], counts.tolist()))


# =============================================================
# 2-3. Build the spin-doubled graph + fiber coupling
# =============================================================

def build_spin_doubled_graph(lattice, prescription: str, eps: float):
    """
    Double the Fock graph and add fiber-coupling edges.

    Node ordering in the doubled graph:
        first  N nodes:  spin-UP   copy  (index = scalar_index)
        second N nodes:  spin-DOWN copy  (index = scalar_index + N)

    Parameters
    ----------
    lattice : GeometricLattice
    prescription : str
        'A', 'B', 'C', or 'D'
    eps : float
        Coupling strength.

    Returns
    -------
    A_doubled : np.ndarray   (2N x 2N)
        Adjacency matrix of the doubled graph.
    fiber_weights : dict
        Maps (n, l) -> fiber edge weight used.
    """
    N = lattice.num_states
    A_scalar = lattice.adjacency.toarray().astype(float)

    # block-diagonal: diag(A_scalar, A_scalar)
    A_doubled = np.zeros((2 * N, 2 * N))
    A_doubled[:N, :N] = A_scalar
    A_doubled[N:, N:] = A_scalar

    fiber_weights = {}  # (n,l) -> weight

    for idx, (n, l, m) in enumerate(lattice.states):
        w = _fiber_weight(n, l, m, prescription, eps, lattice.max_n)
        fiber_weights[(n, l)] = w   # same for all m at given (n,l)

        if abs(w) > 0:
            i_up = idx
            i_dn = idx + N
            A_doubled[i_up, i_dn] = abs(w)   # adjacency is non-negative
            A_doubled[i_dn, i_up] = abs(w)

    return A_doubled, fiber_weights


def _fiber_weight(n: int, l: int, m: int,
                  prescription: str, eps: float, n_max: int) -> float:
    """Compute the fiber edge weight for a given prescription."""

    if prescription == 'A':
        # uniform coupling
        return eps

    elif prescription == 'B':
        # l-dependent: eps * l(l+1) / n^3   (vanishes at l=0)
        return eps * l * (l + 1) / (n ** 3)

    elif prescription == 'C':
        # kappa-aware: eps * l(l+1) / (n^3 * (2l+1))
        if l == 0:
            return 0.0
        return eps * l * (l + 1) / (n ** 3 * (2 * l + 1))

    elif prescription == 'D':
        # Hopf-geometric: eps * sin(chi_n)
        chi = fock_chi(n)
        return eps * np.sin(chi)

    else:
        raise ValueError(f"Unknown prescription: {prescription}")


# =============================================================
# 4. Eigenvalue analysis
# =============================================================

def analyze_eigenvalues(lattice, A_doubled: np.ndarray,
                        evals_scalar: np.ndarray,
                        prescription: str, eps: float):
    """
    Compute graph Laplacian eigenvalues for the doubled graph and
    compare with the doubly-degenerate scalar spectrum.

    Returns a dict with detailed analysis.
    """
    N = lattice.num_states
    n_max = lattice.max_n

    # Graph Laplacian of doubled graph
    D_doubled = np.diag(A_doubled.sum(axis=1))
    L_doubled = D_doubled - A_doubled
    evals_doubled = np.linalg.eigvalsh(L_doubled)

    # The uncoupled (eps=0) doubled spectrum is just evals_scalar repeated
    evals_uncoupled = np.sort(np.concatenate([evals_scalar, evals_scalar]))

    # Hamiltonian eigenvalues (kappa * L)
    H_evals = KAPPA_FOCK * evals_doubled
    H_evals_uncoupled = KAPPA_FOCK * evals_uncoupled

    # --- per-(n,l) splitting analysis ---
    # Identify which eigenvalues correspond to which (n,l) shell.
    # In the uncoupled case, each Laplacian eigenvalue n^2-1 has
    # degeneracy 2*n^2 (doubled).  With fiber coupling, this degeneracy
    # can be partially or fully lifted.

    # Strategy: sort eigenvalues and assign them to shells based on
    # proximity to the known uncoupled values.

    shell_splittings = {}
    # Collect the known shell eigenvalue targets
    shell_targets = []
    for n in range(1, n_max + 1):
        lam_scalar = n * n - 1   # Laplacian eigenvalue for shell n
        # In the doubled graph, 2*n^2 eigenvalues cluster near this value
        for l in range(n):
            deg_l = 2 * l + 1    # m-degeneracy for this l
            # each (n,l,m) appears twice (up,down), so 2*(2l+1) eigenvalues
            for _ in range(2 * deg_l):
                shell_targets.append((n, l, lam_scalar))

    # Sort the doubled eigenvalues
    sorted_evals = np.sort(evals_doubled)
    # Sort targets by their Laplacian eigenvalue
    shell_targets.sort(key=lambda x: x[2])

    # Assign eigenvalues to (n,l) groups
    nl_groups = defaultdict(list)
    for i, (n, l, lam0) in enumerate(shell_targets):
        if i < len(sorted_evals):
            nl_groups[(n, l)].append(sorted_evals[i])

    # Compute splittings
    for (n, l), vals in nl_groups.items():
        vals = np.array(vals)
        mean_val = np.mean(vals)
        spread = np.max(vals) - np.min(vals)
        # How much does it split relative to eps?
        shell_splittings[(n, l)] = {
            'n': n, 'l': l,
            'count': len(vals),
            'mean_laplacian_eval': float(mean_val),
            'spread_laplacian': float(spread),
            'spread_energy': float(abs(KAPPA_FOCK) * spread),
            'values': vals.tolist(),
            # Physical spin-orbit would give splitting ~ alpha^2 * l(l+1)/n^3
            'physical_so_scale': ALPHA_PHYS**2 * l * (l + 1) / (n ** 3) if l > 0 else 0.0,
        }

    # Spectral gap
    spectral_gap = float(sorted_evals[1] - sorted_evals[0]) if len(sorted_evals) > 1 else 0.0

    return {
        'prescription': prescription,
        'eps': eps,
        'n_max': n_max,
        'N_scalar': N,
        'N_doubled': 2 * N,
        'evals_doubled': sorted_evals.tolist(),
        'evals_uncoupled': evals_uncoupled.tolist(),
        'H_evals_doubled': np.sort(H_evals).tolist(),
        'H_evals_uncoupled': np.sort(H_evals_uncoupled).tolist(),
        'spectral_gap': spectral_gap,
        'shell_splittings': {f"({n},{l})": v for (n, l), v in sorted(shell_splittings.items())},
    }


# =============================================================
# 5. Graph-topological analysis
# =============================================================

def topological_analysis(lattice, A_doubled: np.ndarray,
                         prescription: str, eps: float):
    """
    Compute topological invariants of the doubled graph.
    """
    N = lattice.num_states
    A_scalar = lattice.adjacency.toarray().astype(float)

    n_scalar_edges = int(np.count_nonzero(A_scalar) / 2)
    n_doubled_edges = int(np.count_nonzero(A_doubled) / 2)
    n_fiber_edges = n_doubled_edges - 2 * n_scalar_edges

    density_scalar = np.count_nonzero(A_scalar) / (N * N)
    density_doubled = np.count_nonzero(A_doubled) / (4 * N * N)

    # Degree distribution
    degrees = A_doubled.sum(axis=1)
    d_max = int(degrees.max())
    d_min = int(degrees.min())
    q_max = d_max - 1   # excess degree

    # Ramanujan bound: sqrt(q_max)
    ram_bound = float(np.sqrt(max(q_max, 0)))

    # Get non-trivial eigenvalues of adjacency for Ramanujan check
    A_evals = np.linalg.eigvalsh(A_doubled)
    A_evals_sorted = np.sort(np.abs(A_evals))[::-1]  # descending by magnitude

    # Trivial eigenvalues are the largest (Perron-Frobenius)
    # For bipartite graphs: eigenvalues come in +/- pairs
    # Non-trivial: everything except the top eigenvalue pair
    if len(A_evals_sorted) > 2:
        # largest non-trivial
        mu_max_nt = float(A_evals_sorted[2])
    else:
        mu_max_nt = 0.0

    is_ramanujan = mu_max_nt <= ram_bound + 1e-10

    return {
        'prescription': prescription,
        'eps': eps,
        'n_scalar_edges': n_scalar_edges,
        'n_doubled_edges': n_doubled_edges,
        'n_fiber_edges': n_fiber_edges,
        'density_scalar': float(density_scalar),
        'density_doubled': float(density_doubled),
        'd_min': d_min,
        'd_max': d_max,
        'q_max': q_max,
        'ramanujan_bound': ram_bound,
        'mu_max_nontrivial': mu_max_nt,
        'is_ramanujan': bool(is_ramanujan),
        'ramanujan_deviation': float(mu_max_nt - ram_bound),
    }


# =============================================================
# Main driver
# =============================================================

def run_analysis(n_max: int = 3, eps_values=None):
    """
    Run the full spin-doubled graph analysis.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.
    eps_values : list of float, optional
        Coupling strengths to test.
    """
    if eps_values is None:
        eps_values = [0.01, 0.001, ALPHA_PHYS**2]

    prescriptions = ['A', 'B', 'C', 'D']

    print("=" * 78)
    print(f"  SPIN-DOUBLED GRAPH TEST   n_max = {n_max}")
    print("=" * 78)

    # --- Step 1: scalar graph verification ---
    print(f"\n{'='*60}")
    print(f"  Step 1: Scalar Fock Graph Verification (n_max={n_max})")
    print(f"{'='*60}")

    lattice, evals_scalar = build_scalar_graph(n_max)
    match_ok, eval_groups = verify_scalar_spectrum(lattice, evals_scalar, n_max)

    N = lattice.num_states
    print(f"  Nodes:  {N}")
    print(f"  Edges:  {lattice.num_edges}")
    print(f"  States: {lattice.states[:5]} ... (showing first 5)")
    print(f"\n  Laplacian eigenvalue groups (value: degeneracy):")
    for val, deg in sorted(eval_groups.items()):
        n_shell = int(round(np.sqrt(val + 1)))
        expected_deg = n_shell * n_shell
        flag = "OK" if deg == expected_deg else f"MISMATCH (expected {expected_deg})"
        print(f"    lambda = {val:3d}  (n={n_shell})  degeneracy = {deg:3d}  {flag}")
    print(f"\n  Scalar spectrum verification: {'PASS' if match_ok else 'FAIL'}")

    # --- Steps 2-4: doubled graph with each prescription ---
    all_results = {
        'n_max': n_max,
        'N_scalar': N,
        'scalar_evals': evals_scalar.tolist(),
        'scalar_verification': match_ok,
        'kappa_fock': KAPPA_FOCK,
        'alpha_physical': ALPHA_PHYS,
        'analyses': [],
    }

    for pname in prescriptions:
        for eps in eps_values:
            print(f"\n{'='*60}")
            print(f"  Prescription {pname}  |  eps = {eps:.2e}")
            print(f"{'='*60}")

            # Build doubled graph
            A_doubled, fiber_wts = build_spin_doubled_graph(lattice, pname, eps)

            # Eigenvalue analysis
            result = analyze_eigenvalues(lattice, A_doubled, evals_scalar, pname, eps)

            # Topological analysis
            topo = topological_analysis(lattice, A_doubled, pname, eps)
            result['topology'] = topo

            # Print summary table
            print(f"\n  Fiber edges added: {topo['n_fiber_edges']} "
                  f"(of {topo['n_doubled_edges']} total)")
            print(f"  Density: scalar {topo['density_scalar']:.4f}  "
                  f"doubled {topo['density_doubled']:.4f}")
            print(f"  Ramanujan: {'YES' if topo['is_ramanujan'] else 'NO'}  "
                  f"(deviation {topo['ramanujan_deviation']:+.4f})")
            print(f"  Spectral gap (Laplacian): {result['spectral_gap']:.6f}")

            print(f"\n  Per-shell splitting analysis:")
            print(f"  {'(n,l)':>6s}  {'count':>5s}  {'mean_L':>10s}  "
                  f"{'spread_L':>12s}  {'spread_E':>12s}  "
                  f"{'phys_SO':>12s}  {'l=0 prot?':>10s}")
            print(f"  {'-'*80}")

            for key in sorted(result['shell_splittings'].keys()):
                s = result['shell_splittings'][key]
                n, l = s['n'], s['l']
                protected = "YES" if (l == 0 and s['spread_laplacian'] < 1e-12) else \
                            ("N/A" if l > 0 else "NO!")
                print(f"  ({n},{l}){' ' * (4 - len(f'({n},{l})'))}  "
                      f"{s['count']:5d}  "
                      f"{s['mean_laplacian_eval']:10.4f}  "
                      f"{s['spread_laplacian']:12.2e}  "
                      f"{s['spread_energy']:12.2e}  "
                      f"{s['physical_so_scale']:12.2e}  "
                      f"{protected:>10s}")

            # Check splitting pattern vs n, l
            _analyze_splitting_pattern(result, pname, eps)

            all_results['analyses'].append(result)

    # --- Summary comparison ---
    print(f"\n\n{'='*78}")
    print(f"  SUMMARY: Does any prescription produce fine-structure-like splittings?")
    print(f"{'='*78}")

    _print_summary(all_results, prescriptions, eps_values, n_max)

    return all_results


def _analyze_splitting_pattern(result, pname, eps):
    """
    Check if the splitting pattern resembles fine structure:
      E_fs ~ alpha^2 * Z^4 / (2 n^4) * [n/(j+1/2) - 3/4]

    More practically, the spin-orbit piece goes as l(l+1)/n^3 for
    the leading-order LS coupling.
    """
    splittings = result['shell_splittings']
    print(f"\n  Splitting pattern analysis (prescription {pname}, eps={eps:.2e}):")

    # Collect l>0 splittings to check scaling
    l_gt_0 = []
    for key in sorted(splittings.keys()):
        s = splittings[key]
        if s['l'] > 0 and s['spread_laplacian'] > 1e-14:
            l_gt_0.append(s)

    if len(l_gt_0) < 2:
        print("    Not enough l>0 splittings to check pattern.")
        return

    # Check if spread ~ eps * l(l+1)/n^3
    print(f"\n    Checking spread ~ C * l(l+1)/n^3:")
    for s in l_gt_0:
        n, l = s['n'], s['l']
        ratio_to_ll1_n3 = s['spread_laplacian'] / (l * (l + 1) / n ** 3) if l > 0 else 0
        print(f"      (n={n}, l={l}): spread_L = {s['spread_laplacian']:.2e}  "
              f"l(l+1)/n^3 = {l*(l+1)/n**3:.4f}  "
              f"ratio = {ratio_to_ll1_n3:.4e}")


def _print_summary(all_results, prescriptions, eps_values, n_max):
    """Print the final comparison summary."""

    # For each prescription at the physical epsilon (alpha^2), report
    # whether l=0 states are protected and whether splittings scale correctly
    phys_eps = ALPHA_PHYS ** 2

    for pname in prescriptions:
        # Find the result for this prescription at the physical eps
        res = None
        for r in all_results['analyses']:
            if r['prescription'] == pname and abs(r['eps'] - phys_eps) < 1e-15:
                res = r
                break
        if res is None:
            for r in all_results['analyses']:
                if r['prescription'] == pname:
                    res = r
                    break

        if res is None:
            continue

        print(f"\n  Prescription {pname} (eps = {res['eps']:.2e}):")

        l0_protected = True
        any_splitting = False
        for key in sorted(res['shell_splittings'].keys()):
            s = res['shell_splittings'][key]
            if s['l'] == 0 and s['spread_laplacian'] > 1e-12:
                l0_protected = False
            if s['l'] > 0 and s['spread_laplacian'] > 1e-14:
                any_splitting = True

        print(f"    l=0 states degenerate (protected): {'YES' if l0_protected else 'NO'}")
        print(f"    l>0 states split:                  {'YES' if any_splitting else 'NO'}")

        # Check if the split eigenvalues come in the right j=l+1/2, j=l-1/2 pattern
        # Physical fine structure has 2j+1 = 2l+2 states for j=l+1/2
        #                            and 2j+1 = 2l   states for j=l-1/2
        # So the doubled (n,l) block of 2*(2l+1) = 4l+2 states should split into
        # groups of (2l+2) and (2l) if the pattern is fine-structure-like.
        if any_splitting:
            print(f"    Degeneracy sub-pattern (checking j = l +/- 1/2):")
            for key in sorted(res['shell_splittings'].keys()):
                s = res['shell_splittings'][key]
                if s['l'] == 0:
                    continue
                n, l = s['n'], s['l']
                vals = np.array(s['values'])
                # Try to identify sub-groups by clustering
                if len(vals) > 1:
                    center = np.median(vals)
                    above = vals[vals >= center - 1e-14]
                    below = vals[vals < center - 1e-14]
                    # If clean j-splitting: sizes should be (2l+2) and (2l)
                    expected_up = 2 * l + 2
                    expected_dn = 2 * l
                    total = 2 * (2 * l + 1)
                    # Use k-means-like split: find the best 2-cluster split
                    sorted_v = np.sort(vals)
                    # Find the largest gap
                    gaps = np.diff(sorted_v)
                    if len(gaps) > 0:
                        best_gap_idx = np.argmax(gaps)
                        group1 = sorted_v[:best_gap_idx + 1]
                        group2 = sorted_v[best_gap_idx + 1:]
                        sz1, sz2 = len(group1), len(group2)
                        j_match = ((sz1 == expected_dn and sz2 == expected_up) or
                                   (sz1 == expected_up and sz2 == expected_dn))
                        print(f"      (n={n}, l={l}): total={total}, "
                              f"split = {sz1}+{sz2}, "
                              f"expected j-split = {expected_dn}+{expected_up}, "
                              f"{'MATCH' if j_match else 'NO MATCH'}")
                    else:
                        print(f"      (n={n}, l={l}): single eigenvalue")

    # Final verdict
    print(f"\n  {'='*60}")
    print(f"  VERDICT:")
    print(f"  {'='*60}")

    # Prescription B/C should protect l=0 by construction
    # Check if any prescription gives fine-structure-like pattern
    verdicts = []
    for pname in prescriptions:
        for r in all_results['analyses']:
            if r['prescription'] == pname:
                l0_ok = all(
                    r['shell_splittings'][k]['spread_laplacian'] < 1e-12
                    for k in r['shell_splittings']
                    if r['shell_splittings'][k]['l'] == 0
                )
                l_gt0_split = any(
                    r['shell_splittings'][k]['spread_laplacian'] > 1e-14
                    for k in r['shell_splittings']
                    if r['shell_splittings'][k]['l'] > 0
                )
                verdicts.append((pname, r['eps'], l0_ok, l_gt0_split))
                break  # just check first eps for the summary

    for pname, eps, l0_ok, split in verdicts:
        status = "l=0 protected, l>0 split" if (l0_ok and split) else \
                 "l=0 NOT protected" if not l0_ok else \
                 "No splitting at l>0"
        print(f"    {pname}: {status}")


def run_extended(n_max_values=None, eps_values=None):
    """
    Run at multiple n_max values for convergence check.
    """
    if n_max_values is None:
        n_max_values = [3, 4]
    if eps_values is None:
        eps_values = [0.01, 0.001, ALPHA_PHYS**2]

    all_data = {}
    for nm in n_max_values:
        result = run_analysis(nm, eps_values)
        all_data[f"n_max_{nm}"] = result

    return all_data


# =============================================================
# Deep diagnostic: understand the scalar graph structure
# =============================================================

def deep_scalar_diagnostic(n_max: int = 4):
    """
    The scalar spectrum verification FAILED because the Fock graph
    as built by GeometricLattice is NOT the full S^3 Fock graph --
    it's a tree with specific connectivity (n<->n+1 at same l,m;
    m<->m+1 at same n,l).  The eigenvalues do NOT have n^2 degeneracy.

    This diagnostic shows:
    1. The actual scalar Laplacian eigenvalue structure per shell
    2. What the fiber coupling ACTUALLY does (perturbation theory)
    3. Whether the splitting from fiber coupling is eps-dependent
       or pre-existing in the scalar structure
    """
    print(f"\n{'='*78}")
    print(f"  DEEP DIAGNOSTIC: Understanding the scalar graph structure (n_max={n_max})")
    print(f"{'='*78}")

    lattice, evals_scalar = build_scalar_graph(n_max)
    N = lattice.num_states
    A_scalar = lattice.adjacency.toarray().astype(float)
    D = np.diag(A_scalar.sum(axis=1))
    L = D - A_scalar

    print(f"\n  Graph is a TREE: {lattice.num_edges} edges for {N} nodes "
          f"(tree has V-1 = {N-1} edges)")

    # Eigenvalues by shell
    print(f"\n  Scalar Laplacian eigenvalues (sorted):")
    evals = np.sort(evals_scalar)
    for i, e in enumerate(evals):
        # Find which state this eigenvalue belongs to
        print(f"    [{i:3d}] lambda = {e:8.4f}")

    # Now compute per-state degree
    print(f"\n  Per-state degree analysis:")
    print(f"  {'state':>15s}  {'degree':>6s}  {'neighbors':>30s}")
    print(f"  {'-'*60}")
    for idx, (n, l, m) in enumerate(lattice.states):
        deg = int(A_scalar[idx].sum())
        nbrs = lattice.get_neighbors((n, l, m))
        nbr_str = ', '.join([f"({nn},{nl},{nm})" for nn, nl, nm in nbrs])
        print(f"  ({n},{l},{m:+d}){' '*(10-len(f'({n},{l},{m:+d})'))}"
              f"  {deg:6d}  {nbr_str}")

    # Key insight: what happens to the DOUBLED graph when eps -> 0?
    # In the doubled graph with fiber edges, the Laplacian is:
    #   L_doubled = [[L + W, -W], [-W, L + W]]
    # where W = diag(w_1, ..., w_N) is the fiber weight matrix.
    #
    # The eigenvalues of this 2x2 block matrix are:
    #   lambda_+/- = eigenvalues of (L + W +/- W) = eigenvalues of (L + 2W) and L
    #
    # Wait -- that's wrong. Let me redo:
    #   A_doubled = [[A, W], [W, A]]  where W = diag(fiber weights)
    #   D_doubled = diag(d_i + w_i) for each block
    #   L_doubled = D_doubled - A_doubled
    #             = [[D + W_diag, 0], [0, D + W_diag]] - [[A, W], [W, A]]
    #             = [[D - A + W_diag, -W], [-W, D - A + W_diag]]
    #             = [[L + W_diag, -W], [-W, L + W_diag]]
    #
    # For the special case W = w * I (prescription A, uniform):
    #   L_doubled = [[L + wI, -wI], [-wI, L + wI]]
    #
    # This is a 2x2 block matrix with blocks commuting (L commutes with I),
    # so the eigenvalues are: eigenvalues of (L + wI - wI) = L,
    #                    and: eigenvalues of (L + wI + wI) = L + 2wI.
    # => Spectrum is {lambda_i} union {lambda_i + 2w}
    #
    # So uniform fiber coupling should shift ALL eigenvalues by 2*eps,
    # giving UNIFORM splitting of exactly 2*eps for every eigenvalue.
    # This means:
    # - All states split, including l=0 (no protection)
    # - Splitting is uniform, not l-dependent
    #
    # For prescription B/C (diagonal W not proportional to I):
    #   L_doubled = [[L + W, -W], [-W, L + W]]
    #   eigenvalues: if [L, W] = 0 (they commute -- both are diagonal in
    #     the graph eigenbasis? NO -- W is diagonal in the STATE basis
    #     but L is NOT diagonal in the state basis for a non-trivial graph)
    #
    # Actually W is diagonal in the (n,l,m) basis and L is also expressed
    # in the (n,l,m) basis.  But L is NOT diagonal -- it has off-diagonal
    # elements from the edges.
    #
    # For the 2x2 block Schur decomposition: if v is an eigenvector of L
    # with eigenvalue mu, then [v, +v] and [v, -v] are eigenvectors of
    # the UNCOUPLED doubled graph.  With fiber coupling W, we need to
    # look at the perturbation.

    print(f"\n  ANALYTICAL INSIGHT:")
    print(f"  The doubled adjacency is [[A, W], [W, A]] where W = diag(fiber weights).")
    print(f"  The doubled Laplacian is [[L+W_d, -W], [-W, L+W_d]]")
    print(f"  where W_d = diag(row-sums-of-W) = W since W is diagonal.")
    print(f"")
    print(f"  For UNIFORM w: spectrum = {{mu_i}} union {{mu_i + 2w}}")
    print(f"  => Every eigenvalue splits by exactly 2*eps.")
    print(f"  => l=0 states CANNOT be protected under uniform coupling.")
    print(f"")
    print(f"  For NON-UNIFORM w: the splitting depends on how W projects")
    print(f"  onto the Laplacian eigenvectors.  States whose eigenvectors")
    print(f"  have support entirely on l=0 nodes will split by 0 if w(l=0)=0.")
    print(f"  But in the Fock graph, eigenvectors mix different (n,l,m) states!")
    print(f"  The (3,0,0) state at n_max=4 is connected to (2,0,0) which is")
    print(f"  connected to (2,1,m) states.  So the eigenvectors near the n=3,l=0")
    print(f"  eigenvalue will have tails on l>0 states, and the B/C prescriptions")
    print(f"  will cause indirect splitting of the l=0 eigenvalue.")

    # Verify the uniform (prescription A) analytically
    print(f"\n  VERIFICATION: Prescription A analytical prediction vs numerical:")
    eps_test = 0.001
    A_doubled, _ = build_spin_doubled_graph(lattice, 'A', eps_test)
    D_d = np.diag(A_doubled.sum(axis=1))
    L_d = D_d - A_doubled
    evals_d = np.sort(np.linalg.eigvalsh(L_d))

    # Predicted: union of {lambda_i} and {lambda_i + 2*eps}
    predicted = np.sort(np.concatenate([evals_scalar, evals_scalar + 2*eps_test]))
    max_diff = np.max(np.abs(evals_d - predicted))
    print(f"  eps = {eps_test}")
    print(f"  max |numerical - predicted| = {max_diff:.2e}")
    print(f"  Prediction {'CONFIRMED' if max_diff < 1e-12 else 'FAILED'}: "
          f"uniform fiber gives exact 2*eps shift")

    # Now check: the l>0 "spreads" that are ~O(1) -- are they from the SCALAR graph?
    print(f"\n  CRITICAL TEST: Are the O(1) spreads pre-existing in the scalar graph?")
    print(f"  Scalar Laplacian eigenvalues grouped by proximity to n^2-1:")
    for n in range(1, n_max + 1):
        target = n * n - 1
        # All eigenvalues within +/- 2 of target
        mask = np.abs(evals_scalar - target) < max(2.0, target * 0.5 + 1)
        nearby = evals_scalar[mask]
        if len(nearby) > 0:
            spread = np.max(nearby) - np.min(nearby)
            print(f"    n={n} (target {target}): {len(nearby)} eigenvalues, "
                  f"range [{np.min(nearby):.4f}, {np.max(nearby):.4f}], "
                  f"spread = {spread:.4f}")

    # THE KEY FINDING: the Fock graph Laplacian eigenvalues DO split
    # within each shell because the graph is a tree, not a complete
    # shell structure.  The m-ladder transitions and the n-ladder
    # transitions create different local structures for different (n,l)
    # states, so the eigenvalues are NOT degenerate.
    #
    # The l>0 "splittings" we see in the doubled graph are almost entirely
    # the pre-existing scalar structure, NOT a fiber-coupling effect.

    # Quantify: how much of the "splitting" is from fiber coupling?
    print(f"\n  DECOMPOSITION: scalar spread vs fiber-induced shift")
    print(f"  {'(n,l)':>6s}  {'scalar_spread':>14s}  {'fiber_shift':>14s}  {'ratio':>10s}")
    print(f"  {'-'*50}")

    # For prescription B at eps = alpha^2
    eps_phys = ALPHA_PHYS ** 2
    A_d_B, fwts = build_spin_doubled_graph(lattice, 'B', eps_phys)
    D_d_B = np.diag(A_d_B.sum(axis=1))
    L_d_B = D_d_B - A_d_B
    evals_B = np.sort(np.linalg.eigvalsh(L_d_B))

    # Compare with the UNCOUPLED doubled spectrum
    evals_uncoupled = np.sort(np.concatenate([evals_scalar, evals_scalar]))
    # The uncoupled spectrum has each eigenvalue exactly twice.
    # The coupled spectrum shifts these.

    # Per eigenvalue comparison:
    for n in range(1, n_max + 1):
        for l in range(n):
            # In the uncoupled case, find eigenvalues for this (n,l) sector
            # This is tricky because the graph eigenvectors mix (n,l,m)
            # But we can compare total spreads

            # Scalar: eigenvalues that are "near" the n-shell
            target = n * n - 1
            tol = max(2.0, target * 0.5 + 1)
            scalar_near = evals_scalar[np.abs(evals_scalar - target) < tol]
            scalar_spread = (np.max(scalar_near) - np.min(scalar_near)) if len(scalar_near) > 1 else 0

            # The fiber shift for this (n,l) is at most 2*w
            w = _fiber_weight(n, l, 0, 'B', eps_phys, n_max)
            max_fiber_shift = 2 * abs(w)

            if scalar_spread > 0:
                ratio = max_fiber_shift / scalar_spread
            else:
                ratio = float('inf') if max_fiber_shift > 0 else 0

            print(f"  ({n},{l})    {scalar_spread:14.6f}  {max_fiber_shift:14.2e}  "
                  f"{ratio:10.2e}")

    return {}


# =============================================================
# Corrected analysis: isolate the fiber-coupling effect
# =============================================================

def corrected_analysis(n_max: int = 3):
    """
    The correct way to measure the fiber-induced splitting is:
    1. Compute the SCALAR Laplacian eigenvalues (each with multiplicity 1).
    2. In the UNCOUPLED doubled graph, each eigenvalue has multiplicity 2.
    3. Fiber coupling splits each doublet by an amount that depends on
       how the fiber weight matrix W projects onto that eigenvector.
    4. Measure the doublet splitting, NOT the total spread of the shell.

    For uniform W = eps*I: doublet splitting = exactly 2*eps for all.
    For non-uniform W: doublet splitting = 2 * <psi_i | W | psi_i>
      in first-order perturbation theory.
    """
    print(f"\n{'='*78}")
    print(f"  CORRECTED ANALYSIS: Doublet Splitting from Fiber Coupling (n_max={n_max})")
    print(f"{'='*78}")

    lattice, evals_scalar = build_scalar_graph(n_max)
    N = lattice.num_states
    A_scalar = lattice.adjacency.toarray().astype(float)
    D = np.diag(A_scalar.sum(axis=1))
    L = D - A_scalar

    # Get eigenvectors
    evals, evecs = np.linalg.eigh(L)
    sorted_idx = np.argsort(evals)
    evals = evals[sorted_idx]
    evecs = evecs[:, sorted_idx]

    prescriptions = ['A', 'B', 'C', 'D']
    eps_phys = ALPHA_PHYS ** 2

    for pname in prescriptions:
        print(f"\n  Prescription {pname} (eps = {eps_phys:.2e}):")
        print(f"  {'eval_idx':>8s}  {'scalar_eval':>11s}  {'doublet_split':>13s}  "
              f"{'predicted_1st':>13s}  {'dominant_(n,l)':>15s}  "
              f"{'<psi|W|psi>':>12s}")
        print(f"  {'-'*80}")

        # Build fiber weight vector
        W_diag = np.zeros(N)
        for idx, (n, l, m) in enumerate(lattice.states):
            W_diag[idx] = abs(_fiber_weight(n, l, m, pname, eps_phys, lattice.max_n))

        # In the doubled graph: the doublet splitting of the i-th eigenvector
        # is 2 * <psi_i | W | psi_i> at first order.
        # This is because the 2x2 block structure gives eigenvalues:
        #   lambda_i and lambda_i + 2*<psi_i|W|psi_i>
        # when [L, W] != 0 but W is treated as a perturbation.
        #
        # More precisely: the doubled Laplacian in the eigenbasis of L is
        #   [[lambda_i + V_ij, -V_ij], [-V_ij, lambda_i + V_ij]]
        # where V = U^T W U (U = eigenvector matrix).
        # For each diagonal 2x2 block (i=j):
        #   eigenvalues are lambda_i + V_ii + V_ii = lambda_i + 2*V_ii
        #                and lambda_i + V_ii - V_ii = lambda_i
        # But the off-diagonal V_ij couples different eigenvalues.
        # When V_ii >> V_ij (perturbative), the doublet split is 2*V_ii.

        # Compute V = U^T W U
        V = evecs.T @ np.diag(W_diag) @ evecs

        # Also do the EXACT numerical computation
        A_d, _ = build_spin_doubled_graph(lattice, pname, eps_phys)
        D_d = np.diag(A_d.sum(axis=1))
        L_d = D_d - A_d
        evals_d = np.sort(np.linalg.eigvalsh(L_d))

        # The doubled spectrum should be N pairs of near-degenerate eigenvalues.
        # Extract doublet splittings.
        for i in range(N):
            # The i-th pair in the doubled spectrum
            e_lower = evals_d[2*i]
            e_upper = evals_d[2*i + 1]
            doublet_split = e_upper - e_lower

            # First-order prediction
            V_ii = V[i, i]  # <psi_i | W | psi_i>
            predicted = 2 * V_ii

            # Dominant (n,l) of this eigenvector
            weights = evecs[:, i] ** 2  # probability distribution
            dominant_idx = np.argmax(weights)
            n_dom, l_dom, m_dom = lattice.states[dominant_idx]
            dom_weight = weights[dominant_idx]

            print(f"  {i:8d}  {evals[i]:11.6f}  {doublet_split:13.2e}  "
                  f"{predicted:13.2e}  "
                  f"({n_dom},{l_dom},{m_dom:+d}) [{dom_weight:.2f}]  "
                  f"{V_ii:12.2e}")

        # Summary: does the splitting pattern follow l-dependence?
        print(f"\n  Summary for prescription {pname}:")
        # Check whether l=0-dominated eigenvectors have zero splitting
        l0_splits = []
        lgt0_splits = []
        for i in range(N):
            weights = evecs[:, i] ** 2
            # l-character of this eigenvector
            l_avg = sum(weights[j] * lattice.states[j][1]
                        for j in range(N))
            doublet_split = evals_d[2*i + 1] - evals_d[2*i]
            if l_avg < 0.5:
                l0_splits.append(doublet_split)
            else:
                lgt0_splits.append(doublet_split)

        if l0_splits:
            print(f"    l~0 eigenvectors: max doublet split = {max(l0_splits):.2e}")
        if lgt0_splits:
            print(f"    l>0 eigenvectors: max doublet split = {max(lgt0_splits):.2e}")

    return {}


# =============================================================
# Entrypoint
# =============================================================

if __name__ == "__main__":
    print("Running spin-doubled graph test...\n")

    # Run at n_max = 3 and 4
    results_all = run_extended(n_max_values=[3, 4],
                               eps_values=[0.01, 0.001, ALPHA_PHYS**2])

    # Deep diagnostic
    deep_scalar_diagnostic(n_max=4)

    # Corrected analysis
    corrected_analysis(n_max=3)
    corrected_analysis(n_max=4)

    # Save to JSON
    out_path = os.path.join(os.path.dirname(__file__), "data",
                            "spin_doubled_graph_test.json")

    # Make JSON-serializable: convert numpy arrays and bools
    def _serialize(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

    with open(out_path, 'w') as f:
        json.dump(results_all, f, indent=2, default=_serialize)
    print(f"\nResults saved to {out_path}")
