"""
Bargmann–Segal discretization of the 3D isotropic harmonic oscillator.

This module implements the Sprint 2 construction from Track NK: the HO graph
obtained from the holomorphic (Hardy) sector of the Bargmann–Segal space on
ℂ³, restricted to the SU(3) symmetric irreps (N, 0).

Key facts
---------
Nodes
    Quantum numbers (N, l, m_l) with N ≥ 0, l = N, N-2, ..., (0 or 1),
    m_l = -l,...,+l.  Total count per shell: (N+1)(N+2)/2 = dim (N,0) SU(3).

Diagonal ("Hamiltonian")
    The Euler / number operator N̂ = Σ zᵢ ∂/∂zᵢ acts as N on monomials of
    degree N.  The HO Hamiltonian is ℏω (N̂ + 3/2), so the diagonal entry
    at node (N, l, m_l) is exactly ℏω (N + 3/2).  Rational in units of ℏω.

Adjacency
    Dipole-like transition operators (creation / annihilation of one
    quantum in direction i, corresponding to zᵢ or ∂ᵢ on the Bargmann
    space) connect (N, l, m_l) ↔ (N+1, l±1, m_l').  The SQUARED matrix
    elements of the ``r`` operator (isotropic position) between HO radial
    states (n_r, l) and (n_r', l±1) involve integer-coefficient polynomials
    in N, l divided by 2 (with factors √2 that cancel on squaring).

Rationality
    All squared edge weights ⟨N+1, l±1, m_l' | r² | N, l, m_l⟩-style
    reduced elements are rational.  All diagonal entries are rational
    multiples of ℏω.  No π, no transcendentals.

This module uses ``fractions.Fraction`` throughout for the π-free
verification and provides a ``build_bargmann_graph`` entry point that
produces an exact-arithmetic graph, plus a numerical helper for
eigenvalue checks.

References
----------
- Bargmann (1961), Comm. Pure Appl. Math.  14, 187.
- Segal (1963), in "Applications of Mathematics to Problems in
  Theoretical Physics", Gordon and Breach.
- Jauch & Hill (1940), Phys. Rev. 57, 641 (SU(3) uniqueness for 3D HO).
- Folland & Stein (1974), Comm. Pure Appl. Math. 27, 429 (Kohn Laplacian).
- Hall (1994), J. Funct. Anal. 122, 103 (Segal–Bargmann for compact groups).
- Rowe (2010), "Nuclear Collective Motion" (SU(3) shell model).
- GeoVac Paper 0, Paper 7, Paper 18.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, List, Tuple

import numpy as np
from scipy import sparse


# ── node enumeration ────────────────────────────────────────────────────────

Node = Tuple[int, int, int]   # (N, l, m_l)


def _l_values(N: int) -> List[int]:
    """l = N, N-2, ..., 0 or 1 (same parity as N)."""
    return list(range(N, -1, -2))


def enumerate_nodes(N_max: int) -> List[Node]:
    """
    List of (N, l, m_l) triples for N = 0..N_max, in a canonical order:
    by N ascending, then l ascending, then m_l ascending.

    Total count = Σ_{N=0}^{N_max} (N+1)(N+2)/2.
    """
    nodes: List[Node] = []
    for N in range(N_max + 1):
        for l in sorted(_l_values(N)):
            for m in range(-l, l + 1):
                nodes.append((N, l, m))
    return nodes


def shell_degeneracy(N: int) -> int:
    """Degeneracy of HO shell N (spinless): (N+1)(N+2)/2."""
    return (N + 1) * (N + 2) // 2


def total_nodes(N_max: int) -> int:
    """Σ_{N=0}^{N_max} (N+1)(N+2)/2 = (N_max+1)(N_max+2)(N_max+3)/6."""
    return (N_max + 1) * (N_max + 2) * (N_max + 3) // 6


# ── radial reduced matrix elements of the dipole operator r ────────────────
#
# On the (n_r, l) radial basis of the 3D isotropic harmonic oscillator,
# the operator r connects (n_r, l) to (n_r', l ± 1) with the well-known
# selection rule ΔN = ±1 in the principal quantum number N = 2 n_r + l.
#
# The nonzero reduced radial integrals (in units where ℏ = m = ω = 1) are:
#
#     ⟨n_r, l+1 | r | n_r,   l⟩ = √((n_r + l + 3/2))
#     ⟨n_r-1, l+1 | r | n_r, l⟩ = -√(n_r)
#
# Equivalently, going from (n_r, l) at shell N = 2 n_r + l to shell N+1:
#     - (n_r, l+1), ΔN = +1, matrix element = √(n_r + l + 3/2) = √((N + l + 3)/2)
#        Wait — 2 n_r + (l+1) = 2 n_r + l + 1 = N + 1.  OK.
#     - (n_r+1, l-1), ΔN = +1 (since 2(n_r+1) + (l-1) = N + 1), element = √(n_r + 1)
#        ...via ⟨n_r+1, l-1| r | n_r, l⟩ = √(n_r + 1) · (sign).
#
# The SQUARED matrix elements are rational.  Specifically:
#
#     |⟨n_r, l+1 | r | n_r, l⟩|²       = n_r + l + 3/2 = (2 n_r + 2 l + 3)/2
#     |⟨n_r+1, l-1 | r | n_r, l⟩|²     = n_r + 1
#
# We store these as Fraction.  The Wigner–Eckart angular part for a rank-1
# tensor (r is a rank-1 spherical tensor) is given by a Wigner 3j symbol.
# The SQUARED 3j symbols for the (l, 1, l±1) coupling are also rational
# (in fact, very simple rationals — no surds at all in the squared form).


def radial_r_squared_up_lplus(n_r: int, l: int) -> Fraction:
    """
    |⟨n_r, l+1 | r | n_r, l⟩|² (shell N = 2 n_r + l → shell N+1 = 2 n_r + l + 1)

    = n_r + l + 3/2 = (2 n_r + 2 l + 3)/2
    """
    return Fraction(2 * n_r + 2 * l + 3, 2)


def radial_r_squared_up_lminus(n_r: int, l: int) -> Fraction:
    """
    |⟨n_r+1, l-1 | r | n_r, l⟩|² (shell N → shell N+1; l must be ≥ 1)

    = n_r + 1
    """
    if l < 1:
        return Fraction(0)
    return Fraction(n_r + 1, 1)


# ── angular (Wigner–Eckart) squared 3j ─────────────────────────────────────
#
# For a rank-1 spherical tensor T^1_q acting between |l, m⟩ and |l', m'⟩,
# the Wigner 3j symbol that appears is
#
#     ( l  1  l' )
#     ( m  q  -m' )
#
# and the sum over m, q, m' of |3j|² gives (2l+1)^{-1} · (l' selection).
# The SQUARED 3j for ΔL = ±1 rank-1 transitions between SO(3) states are
# rational numbers.  We compute them directly with exact arithmetic using
# the closed-form formulae.
#
# We only need the m-summed reduced squared matrix element for the graph
# edge (since the graph is naturally built on (N, l, m_l) with edges summed
# over intermediate m_l', or alternatively on (N, l) shell-reduced nodes).
#
# For the simple graph built here, we connect (N, l, m) ↔ (N+1, l', m')
# with m' - m = q ∈ {-1, 0, +1} whenever the 3j is nonzero, and take the
# SQUARED matrix element as the edge weight.  These are rational by the
# explicit formulas below.

def wigner3j_rank1_squared(l: int, l_prime: int, m: int, m_prime: int) -> Fraction:
    """
    |( l   1   l' )|²
    |( m   q   -m')|      with q = m' - m (implicitly determined).

    Nonzero only for |l - l'| = 1 (or 0, but we skip l=l' here since we
    want ΔN = ±1 with Δl = ±1 selection for the simple dipole graph).

    Closed-form values (cf. Varshalovich §8.5, Edmonds §5.4).
    """
    q = m_prime - m
    if abs(q) > 1:
        return Fraction(0)

    # Triangle condition
    if abs(l - l_prime) != 1:
        return Fraction(0)

    # Case l' = l + 1
    if l_prime == l + 1:
        if q == 0:
            # (l 1 l+1; m 0 -m)² = (l+1-m)(l+1+m) / [(2l+1)(2l+2)(2l+3)]
            num = (l + 1 - m) * (l + 1 + m)
        elif q == 1:
            # (l 1 l+1; m 1 -(m+1))² = (l+m+1)(l+m+2) / [2(2l+1)(2l+2)(2l+3)]
            num_full = Fraction((l + m + 1) * (l + m + 2),
                                2 * (2 * l + 1) * (2 * l + 2) * (2 * l + 3))
            return num_full
        else:  # q == -1
            # (l 1 l+1; m -1 -(m-1))² = (l-m+1)(l-m+2) / [2(2l+1)(2l+2)(2l+3)]
            num_full = Fraction((l - m + 1) * (l - m + 2),
                                2 * (2 * l + 1) * (2 * l + 2) * (2 * l + 3))
            return num_full
        den = (2 * l + 1) * (2 * l + 2) * (2 * l + 3)
        return Fraction(num, den)

    # Case l' = l - 1
    if l_prime == l - 1:
        if q == 0:
            # (l 1 l-1; m 0 -m)² = (l-m)(l+m) / [(2l-1)(2l)(2l+1)]
            num = (l - m) * (l + m)
            den = (2 * l - 1) * (2 * l) * (2 * l + 1)
            return Fraction(num, den)
        if q == 1:
            # (l 1 l-1; m 1 -(m+1))² = (l-m)(l-m-1) / [2(2l-1)(2l)(2l+1)]
            return Fraction((l - m) * (l - m - 1),
                            2 * (2 * l - 1) * (2 * l) * (2 * l + 1))
        # q == -1
        return Fraction((l + m) * (l + m - 1),
                        2 * (2 * l - 1) * (2 * l) * (2 * l + 1))

    return Fraction(0)


# ── graph construction ────────────────────────────────────────────────────

@dataclass
class BargmannGraph:
    """
    Bargmann–Segal HO graph on the SU(3) (N, 0) sector.

    Attributes
    ----------
    N_max
        Maximum shell index.
    nodes
        List of (N, l, m) triples in canonical order.
    index
        Dict mapping (N, l, m) → integer index 0..n_nodes-1.
    diagonal
        List of Fraction values: HO energy (N + 3/2) at each node
        (in units of ℏω).
    adjacency
        Dict ((i, j) → Fraction) with i < j storing squared dipole matrix
        elements connecting ΔN = +1, Δl = ±1 nodes.  Symmetric.
    """

    N_max: int
    nodes: List[Node]
    index: Dict[Node, int]
    diagonal: List[Fraction]
    adjacency: Dict[Tuple[int, int], Fraction] = field(default_factory=dict)

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    def hamiltonian_dense(self, hw: float = 1.0) -> np.ndarray:
        """Return the diagonal HO Hamiltonian as a dense numpy matrix (units of ℏω)."""
        n = self.n_nodes
        H = np.zeros((n, n), dtype=float)
        for i, d in enumerate(self.diagonal):
            H[i, i] = float(d) * hw
        return H

    def adjacency_dense(self) -> np.ndarray:
        """Return the adjacency (squared dipole matrix elements) as a dense numpy matrix."""
        n = self.n_nodes
        A = np.zeros((n, n), dtype=float)
        for (i, j), w in self.adjacency.items():
            wf = float(w)
            A[i, j] = wf
            A[j, i] = wf
        return A

    def graph_laplacian_dense(self) -> np.ndarray:
        """Return L = D - A."""
        A = self.adjacency_dense()
        D = np.diag(A.sum(axis=1))
        return D - A


def build_bargmann_graph(N_max: int) -> BargmannGraph:
    """
    Build the Bargmann–Segal HO graph up to shell N_max using pure
    rational arithmetic.

    Nodes: (N, l, m_l) with l of same parity as N, 0 ≤ l ≤ N.
    Diagonal: Fraction(2N + 3, 2)  (= N + 3/2, units of ℏω).
    Adjacency: squared dipole matrix elements between (N, l, m) and
               (N+1, l', m') with l' = l ± 1, q = m' - m ∈ {-1, 0, +1}.
               The squared element is (radial r² factor) × (squared 3j).

    All arithmetic is exact (Fraction).
    """
    if N_max < 0:
        raise ValueError("N_max must be ≥ 0")

    nodes = enumerate_nodes(N_max)
    idx: Dict[Node, int] = {nd: i for i, nd in enumerate(nodes)}

    # Diagonal: HO energy = N + 3/2
    diagonal = [Fraction(2 * N + 3, 2) for (N, _, _) in nodes]

    adjacency: Dict[Tuple[int, int], Fraction] = {}

    # Populate edges: ΔN = +1, Δl = ±1
    for (N, l, m) in nodes:
        if N + 1 > N_max:
            continue
        n_r = (N - l) // 2  # radial quantum number for this (N, l) shell

        # (a) l → l + 1:  (n_r, l) → (n_r, l+1), radial r² = n_r + l + 3/2
        l_up = l + 1
        if l_up in _l_values(N + 1):
            radial_sq = radial_r_squared_up_lplus(n_r, l)  # Fraction
            for m_up in range(-l_up, l_up + 1):
                q = m_up - m
                if abs(q) > 1:
                    continue
                ang_sq = wigner3j_rank1_squared(l, l_up, m, m_up)
                if ang_sq == 0:
                    continue
                w = radial_sq * ang_sq
                i, j = idx[(N, l, m)], idx[(N + 1, l_up, m_up)]
                a, b = (i, j) if i < j else (j, i)
                adjacency[(a, b)] = adjacency.get((a, b), Fraction(0)) + w

        # (b) l → l - 1:  (n_r, l) → (n_r+1, l-1), radial r² = n_r + 1
        if l >= 1:
            l_dn = l - 1
            if l_dn in _l_values(N + 1):
                radial_sq = radial_r_squared_up_lminus(n_r, l)
                for m_dn in range(-l_dn, l_dn + 1):
                    q = m_dn - m
                    if abs(q) > 1:
                        continue
                    ang_sq = wigner3j_rank1_squared(l, l_dn, m, m_dn)
                    if ang_sq == 0:
                        continue
                    w = radial_sq * ang_sq
                    i, j = idx[(N, l, m)], idx[(N + 1, l_dn, m_dn)]
                    a, b = (i, j) if i < j else (j, i)
                    adjacency[(a, b)] = adjacency.get((a, b), Fraction(0)) + w

    return BargmannGraph(
        N_max=N_max,
        nodes=nodes,
        index=idx,
        diagonal=diagonal,
        adjacency=adjacency,
    )


# ── π-free verification ────────────────────────────────────────────────────

def verify_pi_free(N_max: int) -> Dict:
    """
    Build the Bargmann graph with rational arithmetic and certify that
    every entry of the diagonal and adjacency is rational.

    Returns a certificate dict with keys:
        N_max
        n_nodes, n_edges
        all_diagonal_rational
        all_adjacency_rational
        irrationals_encountered  (list, should be empty)
    """
    g = build_bargmann_graph(N_max)
    irrationals: List = []
    all_diag = all(isinstance(d, Fraction) for d in g.diagonal)
    all_adj = all(isinstance(w, Fraction) for w in g.adjacency.values())
    for d in g.diagonal:
        if not isinstance(d, Fraction):
            irrationals.append(("diagonal", d))
    for key, w in g.adjacency.items():
        if not isinstance(w, Fraction):
            irrationals.append(("adjacency", key, w))
    return {
        "N_max": N_max,
        "n_nodes": g.n_nodes,
        "n_edges": len(g.adjacency),
        "all_diagonal_rational": all_diag,
        "all_adjacency_rational": all_adj,
        "irrationals_encountered": irrationals,
        "pi_free": all_diag and all_adj and not irrationals,
    }


# ── eigenvalue helpers (numerical, for tests and comparison) ───────────────

def bargmann_ho_spectrum(N_max: int, hw: float = 1.0) -> np.ndarray:
    """
    Return the HO spectrum produced by the diagonal of the Bargmann graph
    (which is exact, not an approximation).

    Sorted ascending; each energy ℏω(N + 3/2) appears with multiplicity
    (N+1)(N+2)/2.
    """
    g = build_bargmann_graph(N_max)
    H = g.hamiltonian_dense(hw=hw)
    return np.sort(np.diag(H))


def bargmann_graph_laplacian_spectrum(N_max: int) -> np.ndarray:
    """
    Return the eigenvalues of the graph Laplacian L = D - A of the
    Bargmann adjacency.  These do NOT equal the HO spectrum — the
    HO spectrum is encoded in the diagonal, not in (D - A).  The
    Laplacian spectrum reflects the SU(3) recoupling structure.
    """
    g = build_bargmann_graph(N_max)
    L = g.graph_laplacian_dense()
    return np.sort(np.linalg.eigvalsh(L))


# ── public API ─────────────────────────────────────────────────────────────

__all__ = [
    "Node",
    "BargmannGraph",
    "enumerate_nodes",
    "shell_degeneracy",
    "total_nodes",
    "radial_r_squared_up_lplus",
    "radial_r_squared_up_lminus",
    "wigner3j_rank1_squared",
    "build_bargmann_graph",
    "verify_pi_free",
    "bargmann_ho_spectrum",
    "bargmann_graph_laplacian_spectrum",
]
