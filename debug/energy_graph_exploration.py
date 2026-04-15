"""
Characterization sprint: the V_ee "energy graph" on S^3 pair-states.

Nodes  = symmetric (singlet, spin-adapted) spatial pair-states |(ab)> for
         spatial orbitals (n,l,m) with n <= n_max, restricted to total
         M_L = 0 so the graph is finite and physically canonical.
Edges  = the two-electron Coulomb matrix element <(ab)|1/r12|(cd)>, which
         evaluates through the GeoVac hypergeometric Slater integrals
         (exact float, machine precision).

This is a STRUCTURAL CHARACTERIZATION, not a solver. No energies. No fits.

Output:
    debug/data/energy_graph_nmax3.json
    debug/data/energy_graph_nmax4.json  (if tractable)
    debug/energy_graph_exploration.md   (written separately)

PM sprint: "V_ee as its own graph" (April 2026).
"""

from __future__ import annotations

import json
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import List, Tuple

import numpy as np

from geovac.casimir_ci import (
    _build_orbital_basis,
    _build_graph_h1,
    two_electron_integral,
)
from geovac.lattice import GeometricLattice


# ---------------------------------------------------------------------------
# Node set: singlet (symmetric) pair states at total M_L = 0
# ---------------------------------------------------------------------------

def build_pair_nodes(n_max: int, m_total: int = 0):
    """Enumerate singlet pair-state nodes: (i,j) with i<=j and m_i+m_j=m_total."""
    orbitals = _build_orbital_basis(n_max)
    pairs = []
    for i in range(len(orbitals)):
        for j in range(i, len(orbitals)):
            if orbitals[i][2] + orbitals[j][2] == m_total:
                pairs.append((i, j))
    return pairs, orbitals


def pair_label(pair, orbitals):
    (i, j) = pair
    return f"({orbitals[i]},{orbitals[j]})"


# ---------------------------------------------------------------------------
# Build V_ee adjacency matrix in the singlet pair-state basis
# ---------------------------------------------------------------------------

def build_vee_matrix(n_max: int, m_total: int = 0):
    """V_ee only (no h1) in singlet pair basis, via casimir_ci's integrator."""
    pairs, orbitals = build_pair_nodes(n_max, m_total=m_total)
    N = len(pairs)
    V = np.zeros((N, N))

    def g(a, b, c, d):
        na, la, ma = orbitals[a]
        nb, lb, mb = orbitals[b]
        nc, lc, mc = orbitals[c]
        nd, ld, md = orbitals[d]
        return two_electron_integral(
            na, la, ma, nb, lb, mb,
            nc, lc, mc, nd, ld, md,
            k_orb=1.0,
        )

    for I, (i, j) in enumerate(pairs):
        bra = [(i, j)] + ([(j, i)] if i != j else [])
        NI = np.sqrt(len(bra))
        for J in range(I, N):
            p, q = pairs[J]
            ket = [(p, q)] + ([(q, p)] if p != q else [])
            NJ = np.sqrt(len(ket))
            me = 0.0
            for (a, b) in bra:
                for (c, d) in ket:
                    me += g(a, b, c, d)
            me /= (NI * NJ)
            V[I, J] = me
            V[J, I] = me
    return V, pairs, orbitals


# ---------------------------------------------------------------------------
# Wavefunction-graph Laplacian projected to the same pair basis
# ---------------------------------------------------------------------------

def build_h1_pair_matrix(n_max: int, Z: int, m_total: int = 0):
    """The one-body GeoVac graph-h1 operator in the same singlet pair basis.

    h1 is a one-body operator; its action on |(ij)> gives the symmetric
    pair sum h1|i> |j> + |i> h1|j>.  We project onto the pair basis to
    make the one-body graph Laplacian live on the same node set as V.
    """
    h1_spatial, orbitals = _build_graph_h1(Z, n_max)
    pairs, _ = build_pair_nodes(n_max, m_total=m_total)
    N = len(pairs)
    H1 = np.zeros((N, N))

    for I, (i, j) in enumerate(pairs):
        bra = [(i, j)] + ([(j, i)] if i != j else [])
        NI = np.sqrt(len(bra))
        for J in range(I, N):
            p, q = pairs[J]
            ket = [(p, q)] + ([(q, p)] if p != q else [])
            NJ = np.sqrt(len(ket))
            me = 0.0
            for (a, b) in bra:
                for (c, d) in ket:
                    # h1(1)*delta(b,d) + delta(a,c)*h1(b,d)
                    if b == d:
                        me += h1_spatial[a, c]
                    if a == c:
                        me += h1_spatial[b, d]
            me /= (NI * NJ)
            H1[I, J] = me
            H1[J, I] = me
    return H1, pairs, orbitals


# ---------------------------------------------------------------------------
# Symmetry labels: total L, parity, m_total (already pinned)
# ---------------------------------------------------------------------------

def symmetry_labels(pairs, orbitals):
    """For each node, tabulate (l1+l2 parity, |l1-l2| min_L, l1+l2 max_L)."""
    lab = []
    for (i, j) in pairs:
        n1, l1, m1 = orbitals[i]
        n2, l2, m2 = orbitals[j]
        lab.append(dict(
            l1=l1, l2=l2,
            parity=(l1 + l2) % 2,
            L_min=abs(l1 - l2),
            L_max=l1 + l2,
            n1=n1, n2=n2,
            m1=m1, m2=m2,
        ))
    return lab


# ---------------------------------------------------------------------------
# Characterization driver
# ---------------------------------------------------------------------------

def characterize(n_max: int, tol: float = 1e-12):
    t0 = time.time()
    V, pairs, orbitals = build_vee_matrix(n_max, m_total=0)
    H1, _, _ = build_h1_pair_matrix(n_max, Z=2, m_total=0)
    build_time = time.time() - t0

    N = V.shape[0]
    # ------------- Section 2: sparsity -----------------
    nz_mask = np.abs(V) > tol
    n_nz_total = int(nz_mask.sum())
    n_nz_offdiag = int((nz_mask & ~np.eye(N, dtype=bool)).sum())
    total_offdiag = N * N - N
    density_offdiag = n_nz_offdiag / total_offdiag if total_offdiag else 0.0
    density_total = n_nz_total / (N * N)

    # Selection-rule check: zero V edges should mostly come from Gaunt.
    # Count edges killed by (parity mismatch) of (l1+l2) across pair-states.
    lab = symmetry_labels(pairs, orbitals)
    parity_block_pairs = Counter((lab[i]['parity'], lab[j]['parity'])
                                 for i in range(N) for j in range(N))
    parity_killed = 0
    for i in range(N):
        for j in range(N):
            if lab[i]['parity'] != lab[j]['parity'] and abs(V[i, j]) < tol:
                parity_killed += 1

    # ------------- Section 3: spectrum -----------------
    eigvals_V = np.linalg.eigvalsh(V)
    eigvals_V = np.sort(eigvals_V)
    # degeneracy counts (round to 10 decimals)
    rounded = np.round(eigvals_V, 10)
    degens = Counter(rounded.tolist())
    n_distinct = len(degens)

    # Try integer-rational detection: see if eigvals cluster on small rationals
    # (we don't expect this on V since 1/r integrals are not integer, but check)
    near_rational = []
    for ev in eigvals_V[:20]:
        # rational approx with denominator <= 256
        from fractions import Fraction
        f = Fraction(float(ev)).limit_denominator(256)
        if abs(float(f) - ev) < 1e-6:
            near_rational.append((float(ev), str(f)))

    # Wavefunction-graph spectrum on same pair basis (H1)
    eigvals_H1 = np.sort(np.linalg.eigvalsh(H1))

    # ------------- Section 4: degree distribution -----------------
    # weighted degree = sum |V_ij| over j != i
    W = np.abs(V).copy()
    np.fill_diagonal(W, 0.0)
    deg = W.sum(axis=1)
    diag = np.diag(V).copy()

    # unweighted degree: count nonzero off-diagonal edges
    udeg = (W > tol).sum(axis=1)

    # sort to find hubs
    hub_idx = np.argsort(-deg)
    top_hubs = []
    for r in hub_idx[:min(10, N)]:
        top_hubs.append(dict(
            node=pair_label(pairs[r], orbitals),
            weighted_degree=float(deg[r]),
            unweighted_degree=int(udeg[r]),
            diagonal=float(diag[r]),
        ))

    # ------------- Section 5: recurrence search -----------------
    # Test whether V[(n1,l,m),(n2,l,m)|(n1',l,m),(n2',l,m)]
    # satisfies a three-term rule in (n1, n2).
    # Simplest check: diagonal <(1s,ns)|1/r|(1s,ns)> as function of n.
    # Also check <(1s,ns)|1/r|(1s,(n+1)s)> versus <(1s,n+1)s|...|...>
    recurrence_data = {}
    # Find indices of (1s 1s), (1s 2s), (1s 3s), ..., and (2s 2s), (2s 3s), ...
    def find_pair(nA, lA, nB, lB):
        try:
            ia = orbitals.index((nA, lA, 0))
            ib = orbitals.index((nB, lB, 0))
        except ValueError:
            return None
        a, b = sorted((ia, ib))
        if (a, b) in pairs:
            return pairs.index((a, b))
        return None

    # (1s ns) diagonal V_ee for n = 1..n_max
    vee_1s_ns_diag = {}
    for n in range(1, n_max + 1):
        idx = find_pair(1, 0, n, 0)
        if idx is not None:
            vee_1s_ns_diag[n] = float(V[idx, idx])

    # (1s ns) -> (1s (n+1)s) coupling: off-diagonal V for fixed first electron 1s
    vee_1s_ns_offdiag = {}
    for n in range(1, n_max):
        i1 = find_pair(1, 0, n, 0)
        i2 = find_pair(1, 0, n + 1, 0)
        if i1 is not None and i2 is not None:
            vee_1s_ns_offdiag[n] = float(V[i1, i2])

    # (ns,ns) diagonal V_ee for n = 1..n_max
    vee_nsns_diag = {}
    for n in range(1, n_max + 1):
        idx = find_pair(n, 0, n, 0)
        if idx is not None:
            vee_nsns_diag[n] = float(V[idx, idx])

    # Ratios and differences for recurrence hunt
    ratios = {}
    keys = sorted(vee_1s_ns_diag)
    for i in range(len(keys) - 1):
        n = keys[i]; m = keys[i + 1]
        if abs(vee_1s_ns_diag[m]) > 1e-15:
            ratios[f"V(1s{n}s)/V(1s{m}s)"] = vee_1s_ns_diag[n] / vee_1s_ns_diag[m]

    recurrence_data = dict(
        vee_1s_ns_diag=vee_1s_ns_diag,
        vee_1s_ns_offdiag=vee_1s_ns_offdiag,
        vee_nsns_diag=vee_nsns_diag,
        ratios=ratios,
    )

    # ------------- Section 6: compare with H1 wave graph -----------------
    # Commutator norm [H1, V]
    C = H1 @ V - V @ H1
    commutator_norm = float(np.linalg.norm(C))
    relative_comm_norm = commutator_norm / max(np.linalg.norm(H1), 1e-12)
    # Shared eigenvectors?  If they commute, they share a basis; measure how
    # much the spectrum of V splits within each H1 eigenspace.
    eH1, VH1 = np.linalg.eigh(H1)
    V_in_H1basis = VH1.T @ V @ VH1
    V_offdiag_in_H1 = np.sqrt(max(0.0, float(np.sum(V_in_H1basis ** 2)
                                         - np.sum(np.diag(V_in_H1basis) ** 2))))
    V_frob = float(np.linalg.norm(V))

    # ------------- Section 7: cusp signature -----------------
    # Which diagonal V entries are largest?  The "cusp" concentrates on tight
    # pair states (both electrons at small n).  Report ordering of V_ii
    diag_sorted_idx = np.argsort(-np.diag(V))
    diag_hot = []
    for r in diag_sorted_idx[:min(10, N)]:
        diag_hot.append(dict(
            node=pair_label(pairs[r], orbitals),
            Vii=float(V[r, r]),
        ))

    # Also: largest off-diagonal entry
    absV = np.abs(V - np.diag(np.diag(V)))
    i_hot, j_hot = np.unravel_index(np.argmax(absV), absV.shape)
    hottest_offdiag = dict(
        a=pair_label(pairs[i_hot], orbitals),
        b=pair_label(pairs[j_hot], orbitals),
        value=float(V[i_hot, j_hot]),
    )

    summary = dict(
        n_max=n_max,
        m_total=0,
        n_nodes=N,
        build_seconds=build_time,
        sparsity=dict(
            total_entries=N * N,
            nonzero_total=n_nz_total,
            nonzero_offdiag=n_nz_offdiag,
            density_total=density_total,
            density_offdiag=density_offdiag,
            parity_killed_offdiag=parity_killed,
        ),
        spectrum=dict(
            n_distinct_eigenvalues=n_distinct,
            eigenvalues=[float(x) for x in eigvals_V],
            eigenvalues_H1=[float(x) for x in eigvals_H1],
            near_rational_first_20=near_rational,
            degeneracy_counts=[(float(k), int(v)) for k, v in sorted(degens.items())],
        ),
        degree_distribution=dict(
            top_hubs=top_hubs,
            diag_hot=diag_hot,
            unweighted_degree_histogram=dict(Counter(udeg.tolist())),
        ),
        recurrence=recurrence_data,
        h1_vs_v=dict(
            commutator_norm=commutator_norm,
            relative_commutator_norm=relative_comm_norm,
            V_frobenius=V_frob,
            V_offdiagonal_in_H1_eigenbasis=V_offdiag_in_H1,
            V_diagonal_fraction_in_H1_eigenbasis=(
                float(np.sqrt(np.sum(np.diag(V_in_H1basis) ** 2))) / V_frob
                if V_frob > 0 else 0.0
            ),
        ),
        cusp_signature=dict(
            diagonal_top=diag_hot,
            hottest_offdiag=hottest_offdiag,
        ),
    )
    return summary, V, H1, pairs, orbitals


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    outdir = Path("debug/data")
    outdir.mkdir(parents=True, exist_ok=True)

    for nmax, tag in [(3, "nmax3"), (4, "nmax4")]:
        print(f"\n=== n_max = {nmax} ===")
        t0 = time.time()
        summary, V, H1, pairs, orbitals = characterize(nmax)
        elapsed = time.time() - t0
        summary['total_seconds'] = elapsed
        print(f"  N_nodes = {summary['n_nodes']}")
        print(f"  density_offdiag = {summary['sparsity']['density_offdiag']:.3f}")
        print(f"  n_distinct_eigs = {summary['spectrum']['n_distinct_eigenvalues']}")
        print(f"  [H1,V] rel norm = {summary['h1_vs_v']['relative_commutator_norm']:.3e}")
        print(f"  elapsed = {elapsed:.1f}s")

        outpath = outdir / f"energy_graph_{tag}.json"
        with open(outpath, "w") as f:
            json.dump(summary, f, indent=2, default=str)
        print(f"  saved -> {outpath}")

        if elapsed > 300 and nmax == 3:
            print("  (skipping n_max=4: n_max=3 was already slow)")
            break


if __name__ == "__main__":
    main()
