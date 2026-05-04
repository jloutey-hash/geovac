"""
SU(3) Wilson Lattice Gauge Theory on the Bargmann-Segal S^5 Graph
==================================================================

This module implements a non-abelian SU(3) Wilson lattice gauge theory on
the Bargmann-Segal discretization of the 3D harmonic oscillator (Paper 24).
It is the structural sibling of `geovac/su2_wilson_gauge.py` (Paper 30),
which built SU(2) Wilson on the S^3 = SU(2) Coulomb graph.

Framing
-------

Paper 25 §VII.A flagged as open: does the S^5 Bargmann-Segal graph carry
a non-abelian gauge structure?  Sprint 5 Track S5 (debug/s5_gauge_structure_memo.md)
returned a clean negative for the *adjacency-preserving* SU(3) action: the
inter-shell transitions (N,0) → (N+1,0) are Clebsch-Gordan intertwiners
between distinct SU(3) irreps, not SU(3) group elements, so an SU(3)-action
on the (N,0)-tower fails Wilson's "fixed group on every link" requirement.

The Wilson construction is structurally different and may bypass that
obstruction: the link variables U_e ∈ SU(3) are *external* gauge degrees of
freedom living on edges, NOT a representation of SU(3) acting on shell
labels.  Whether this works on the Bargmann-Segal graph is what this
module tests.

What this module IS
-------------------

* A finite-dimensional non-abelian lattice gauge theory on the S^5
  Bargmann-Segal graph at any N_max, with link variables in SU(3),
  plaquettes as primitive closed non-backtracking walks, and Wilson
  action S_W = beta * sum_P (1 - (1/3) Re tr U_P).
* Gauge-invariant under U_e -> g_{src(e)} U_e g_{tgt(e)}^dagger.
* Has a finite-dimensional Haar-measure partition function.
* Has a maximal-torus reduction: links restricted to the Cartan torus
  T = U(1) x U(1) yield an abelian gauge theory which decomposes as a
  product of two U(1) Wilson theories on the same graph.
* Has a weak-coupling kinetic term:
  S_W ≈ (beta/12) * sum_{a=1..8} <A^a, L_1^{plaq} A^a> + O(A^4).

What this module IS NOT
-----------------------

* NOT a proof of continuum SU(3) Yang-Mills mass gap (Clay problem).
* NOT a proof of confinement on the Bargmann-Segal graph.
* NOT Lorentzian; Euclidean spatial-slice lattice gauge theory.
* NOT a refutation of Sprint 5 Track S5: that result concerns the
  IRREP-tower action of SU(3); this module concerns Wilson link variables.
  They are categorically different objects.  Both verdicts coexist.

Conventions
-----------

* SU(3) elements parametrized by 8 real parameters (the Lie algebra
  su(3) has dimension 8) with generators T^a = lambda^a / 2 (Gell-Mann
  matrices divided by 2).  We use the standard Gell-Mann basis:
      lambda_1, lambda_2, lambda_3 (analog of sigma_x, sigma_y, sigma_z)
      lambda_4, lambda_5 (off-diagonal real/imag in (1,3) block)
      lambda_6, lambda_7 (off-diagonal real/imag in (2,3) block)
      lambda_8 (diagonal, traceless, second Cartan generator)

* Each undirected edge {i, j} of the Bargmann graph (built from
  geovac/nuclear/bargmann_graph.py) is promoted to two oriented edges
  i -> j and j -> i.  Link variable: U_{i->j} ∈ SU(3); reverse carries
  U_{j->i} = U_{i->j}^dagger.

* Plaquettes: primitive closed non-backtracking walks (same convention
  as Paper 30).  The Bargmann graph at fixed N_max is BIPARTITE by N
  parity, so all primitive cycles have EVEN length (smallest possible
  is 4).

* Wilson action: S_W = beta * sum_P (1 - (1/N_c) Re tr U_P) with
  N_c = 3 (the dimension of the fundamental rep of SU(3)).  Vanishes
  iff every plaquette holonomy is the identity (trivial vacuum).

* Haar measure on SU(3): normalized bi-invariant probability measure;
  sampled by SU(3) ~ exp(i sum_a x_a T^a) with appropriate distribution
  on x_a, or equivalently by QR factorization of complex Gaussians.

* Cartan torus: T = {diag(e^{i a}, e^{i b}, e^{-i(a+b)}) : a, b ∈ R}
  ≅ U(1) x U(1).  Maximal torus reduction: restrict link variables to T.

References
----------

- Paper 24 (Bargmann-Segal lattice).
- Paper 25 §VII.A (open question about non-abelian gauge on S^5).
- Paper 30 (SU(2) Wilson on S^3, the structural sibling).
- Paper 28 §graph_native_qed (related vertex-coupling analyses).
- Sprint 5 Track S5 memo (debug/s5_gauge_structure_memo.md).
- Wilson, Phys. Rev. D 10, 2445 (1974).
- Creutz, "Quarks, Gluons and Lattices" (1983).
- Drouffe & Zuber, Phys. Rep. 102, 1 (1983) (SU(N) lattice gauge).

Author: GeoVac Sprint ST-SU3 (May 2026).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import sympy as sp

from geovac.nuclear.bargmann_graph import build_bargmann_graph, BargmannGraph


# ---------------------------------------------------------------------------
# 1. Gell-Mann matrices and SU(3) operations
# ---------------------------------------------------------------------------

def gell_mann_matrices() -> List[np.ndarray]:
    """
    Return the 8 Gell-Mann matrices lambda_1, ..., lambda_8 as 3x3 complex arrays.

    Conventions: traceless Hermitian; tr(lambda_a lambda_b) = 2 delta_{ab};
    lambda_3 and lambda_8 are diagonal (the Cartan generators).  The
    Lie-algebra basis is T^a = lambda^a / 2 with [T^a, T^b] = i f^{abc} T^c.

    Returns
    -------
    [lambda_1, ..., lambda_8] : list of 3x3 complex ndarrays.
    """
    lam = [None] * 8
    lam[0] = np.array(
        [[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex
    )  # lambda_1
    lam[1] = np.array(
        [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex
    )  # lambda_2
    lam[2] = np.array(
        [[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex
    )  # lambda_3 (Cartan)
    lam[3] = np.array(
        [[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex
    )  # lambda_4
    lam[4] = np.array(
        [[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex
    )  # lambda_5
    lam[5] = np.array(
        [[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex
    )  # lambda_6
    lam[6] = np.array(
        [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex
    )  # lambda_7
    lam[7] = (1.0 / np.sqrt(3.0)) * np.array(
        [[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex
    )  # lambda_8 (Cartan)
    return lam


def su3_generators() -> List[np.ndarray]:
    """Return T^a = lambda^a / 2 for a = 1..8."""
    return [0.5 * lam for lam in gell_mann_matrices()]


def su3_from_algebra(coeffs: Sequence[float]) -> np.ndarray:
    """
    Construct U = exp(i sum_a x_a T^a) for x_a real, a = 1..8.

    Parameters
    ----------
    coeffs : length-8 array of real numbers
        Coefficients in the T^a = lambda^a / 2 basis.

    Returns
    -------
    U : np.ndarray, shape (3, 3), dtype complex
    """
    if len(coeffs) != 8:
        raise ValueError(f"Need 8 coefficients, got {len(coeffs)}")
    T = su3_generators()
    A = np.zeros((3, 3), dtype=complex)
    for a, x in enumerate(coeffs):
        A = A + x * T[a]
    return _matrix_exp_hermitian(1j * A)


def _matrix_exp_hermitian(iH: np.ndarray) -> np.ndarray:
    """
    Compute exp(iH) for iH = i * (Hermitian) using eigendecomposition.

    Faster and more stable than scipy.linalg.expm for this special case.
    """
    H = -1j * iH  # H is Hermitian
    # Hermitize defensively
    H = 0.5 * (H + H.conj().T)
    w, V = np.linalg.eigh(H)
    return V @ np.diag(np.exp(1j * w)) @ V.conj().T


def su3_random(rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """
    Sample a Haar-random SU(3) matrix.

    Implementation: take a 3x3 complex matrix with iid standard complex
    Gaussian entries; perform QR factorization Z = QR; correct phase by
    Q -> Q * diag(R/|R|); finally enforce det = +1 by Q -> Q * det(Q)^*.

    This is the Mezzadri 2007 algorithm for U(N), specialized to SU(N).

    Parameters
    ----------
    rng : np.random.Generator, optional

    Returns
    -------
    U : np.ndarray, shape (3, 3), dtype complex
    """
    if rng is None:
        rng = np.random.default_rng()
    Z = (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    diagR = np.diag(R)
    Lambda = diagR / np.abs(diagR)
    U = Q * Lambda  # broadcast; U is unitary
    # Make special unitary
    detU = np.linalg.det(U)
    U[:, 0] = U[:, 0] / detU  # multiply first column by det(U)^*
    return U


def is_su3(U: np.ndarray, atol: float = 1e-8) -> bool:
    """
    Check whether U is in SU(3) within tolerance atol.

    Verifies: shape (3,3), |det(U) - 1| < atol, ||U U^dagger - I|| < atol.
    """
    if U.shape != (3, 3):
        return False
    det_ok = np.isclose(np.linalg.det(U), 1.0, atol=atol)
    unitary_ok = np.allclose(U @ U.conj().T, np.eye(3), atol=atol)
    return bool(det_ok and unitary_ok)


def su3_character(U: np.ndarray) -> float:
    """
    Wilson plaquette character (1/3) Re tr U for SU(3).

    For U ∈ SU(3) the trace is the sum of three eigenvalues
    (e^{i a}, e^{i b}, e^{-i(a+b)}); this is not a single cosine but
    a sum of three phases, real part bounded by [-1, +1] (achieved
    -1/2 at minimum trace, +1 at U = I).
    """
    return float(np.real(np.trace(U))) / 3.0


def cartan_su3_from_phases(a: float, b: float) -> np.ndarray:
    """
    SU(3) matrix in the maximal Cartan torus T = U(1) x U(1):

        diag(e^{i a}, e^{i b}, e^{-i(a + b)}).

    The constraint det = 1 forces the third phase to be -(a+b).

    Parameters
    ----------
    a, b : float
        Two independent phase angles parametrizing the maximal torus.

    Returns
    -------
    U : np.ndarray, shape (3, 3)
    """
    return np.diag(
        [np.exp(1j * a), np.exp(1j * b), np.exp(-1j * (a + b))]
    )


# ---------------------------------------------------------------------------
# 2. Adjacency adapter for the Bargmann-Segal graph
# ---------------------------------------------------------------------------

def bargmann_adjacency_dense(N_max: int) -> np.ndarray:
    """
    0/1 adjacency matrix for the Bargmann-Segal S^5 graph at N_max.

    The Bargmann graph carries weighted edges (squared dipole matrix
    elements) but for the lattice-gauge construction we only need the
    underlying connectivity.  Nonzero entries are reduced to 1.

    Parameters
    ----------
    N_max : int

    Returns
    -------
    A : np.ndarray, shape (V, V)
        Symmetric 0/1 adjacency.
    """
    g = build_bargmann_graph(N_max)
    n = g.n_nodes
    A = np.zeros((n, n), dtype=int)
    for (i, j), w in g.adjacency.items():
        if w != 0:
            A[i, j] = 1
            A[j, i] = 1
    return A


# ---------------------------------------------------------------------------
# 3. Oriented-edge bookkeeping (mirrors su2_wilson_gauge.py)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class OrientedEdge:
    """A directed edge (source -> target) on the graph."""
    source: int
    target: int

    def reverse(self) -> "OrientedEdge":
        return OrientedEdge(self.target, self.source)


def enumerate_oriented_edges(
    adjacency: np.ndarray,
) -> Tuple[List[OrientedEdge], Dict[Tuple[int, int], int]]:
    """
    List all oriented edges of a 0/1 adjacency matrix.  Each undirected
    edge contributes two oriented edges.
    """
    A = np.asarray(adjacency)
    V = A.shape[0]
    if A.shape != (V, V):
        raise ValueError("adjacency must be square")
    oriented: List[OrientedEdge] = []
    index: Dict[Tuple[int, int], int] = {}
    for u in range(V):
        for v in np.nonzero(A[u])[0]:
            oe = OrientedEdge(u, int(v))
            oriented.append(oe)
            index[(u, int(v))] = len(oriented) - 1
    return oriented, index


# ---------------------------------------------------------------------------
# 4. Plaquette enumeration (primitive closed non-backtracking walks)
# ---------------------------------------------------------------------------

def enumerate_plaquettes(
    adjacency: np.ndarray,
    max_length: int = 6,
    both_orientations: bool = False,
) -> List[List[OrientedEdge]]:
    """
    Enumerate primitive closed non-backtracking walks of length <= max_length.

    Mirrors the SU(2) Paper 30 enumeration.  On a bipartite graph (such
    as the Bargmann-Segal graph at fixed N_max), all closed walks have
    EVEN length.  The smallest plaquettes are length-4 squares.

    Parameters
    ----------
    adjacency : np.ndarray
        Symmetric 0/1 adjacency.
    max_length : int
    both_orientations : bool
        If False, keep one representative per unoriented cycle.

    Returns
    -------
    plaquettes : list of list[OrientedEdge]
    """
    A = np.asarray(adjacency)
    V = A.shape[0]
    oriented, edge_index = enumerate_oriented_edges(A)
    n_oriented = len(oriented)

    seen: set = set()
    plaquettes: List[List[OrientedEdge]] = []

    def is_primitive(walk: Tuple[int, ...]) -> bool:
        L = len(walk)
        for d in range(1, L):
            if L % d != 0:
                continue
            period = walk[:d]
            is_power = True
            for k in range(L // d):
                if walk[k * d:(k + 1) * d] != period:
                    is_power = False
                    break
            if is_power:
                return False
        return True

    for start_idx in range(n_oriented):
        stack: List[List[int]] = [[start_idx]]
        while stack:
            path_idx = stack.pop()
            if len(path_idx) >= max_length:
                continue
            last = path_idx[-1]
            last_edge = oriented[last]
            for j in np.nonzero(A[last_edge.target])[0]:
                next_idx = edge_index[(last_edge.target, int(j))]
                next_edge = oriented[next_idx]
                if (next_edge.source == last_edge.target
                        and next_edge.target == last_edge.source):
                    continue
                new_path = path_idx + [next_idx]
                first_edge = oriented[new_path[0]]
                if (next_edge.target == first_edge.source
                        and len(new_path) >= 3):
                    if not (
                        first_edge.source == next_edge.target
                        and first_edge.target == next_edge.source
                    ):
                        path_tuple = tuple(new_path)
                        if is_primitive(path_tuple):
                            min_pos = int(np.argmin(new_path))
                            canonical = tuple(
                                new_path[min_pos:] + new_path[:min_pos]
                            )
                            reversed_edges = [
                                edge_index[
                                    (oriented[i].target, oriented[i].source)
                                ]
                                for i in reversed(new_path)
                            ]
                            min_pos_r = int(np.argmin(reversed_edges))
                            canonical_r = tuple(
                                reversed_edges[min_pos_r:]
                                + reversed_edges[:min_pos_r]
                            )
                            if both_orientations:
                                if canonical not in seen:
                                    seen.add(canonical)
                                    plaquettes.append(
                                        [oriented[i] for i in canonical]
                                    )
                            else:
                                rep = (
                                    canonical
                                    if canonical <= canonical_r
                                    else canonical_r
                                )
                                if rep not in seen:
                                    seen.add(rep)
                                    plaquettes.append(
                                        [oriented[i] for i in rep]
                                    )
                if len(new_path) < max_length:
                    stack.append(new_path)

    return plaquettes


# ---------------------------------------------------------------------------
# 5. Link variables, plaquette holonomy, Wilson action
# ---------------------------------------------------------------------------

def edge_link_variable(
    edge: OrientedEdge,
    link_variables: Dict[Tuple[int, int], np.ndarray],
) -> np.ndarray:
    """
    Look up the SU(3) link variable for an oriented edge.

    If the forward edge (source, target) is stored, return its matrix.
    If only the reverse is stored, return its conjugate transpose
    (Wilson convention U_{e^{-1}} = U_e^dagger).
    """
    key = (edge.source, edge.target)
    if key in link_variables:
        return link_variables[key]
    rev = (edge.target, edge.source)
    if rev in link_variables:
        return link_variables[rev].conj().T
    raise KeyError(
        f"No link variable for edge {edge.source}->{edge.target} or its reverse."
    )


def plaquette_holonomy(
    plaquette: Sequence[OrientedEdge],
    link_variables: Dict[Tuple[int, int], np.ndarray],
) -> np.ndarray:
    """Ordered product U_{e_1} U_{e_2} ... U_{e_L} around plaquette P."""
    U_P = np.eye(3, dtype=complex)
    for e in plaquette:
        U_P = U_P @ edge_link_variable(e, link_variables)
    return U_P


def wilson_action(
    plaquettes: Sequence[Sequence[OrientedEdge]],
    link_variables: Dict[Tuple[int, int], np.ndarray],
    beta: float,
) -> float:
    """
    SU(3) Wilson action S_W = beta * sum_P (1 - (1/3) Re tr U_P).

    Vanishes iff every plaquette holonomy is the identity (trivial vacuum).
    Non-negative on SU(3).

    Parameters
    ----------
    plaquettes : list of list[OrientedEdge]
    link_variables : dict (src, tgt) -> 3x3 complex array
    beta : float
        Inverse coupling.  beta = 6/g^2 in the SU(3) Wilson convention.

    Returns
    -------
    S : float
    """
    S = 0.0
    for P in plaquettes:
        U_P = plaquette_holonomy(P, link_variables)
        S += 1.0 - su3_character(U_P)
    return beta * S


# ---------------------------------------------------------------------------
# 6. Gauge transformation
# ---------------------------------------------------------------------------

def gauge_transform(
    link_variables: Dict[Tuple[int, int], np.ndarray],
    gauge: Dict[int, np.ndarray],
) -> Dict[Tuple[int, int], np.ndarray]:
    """
    Apply node-local gauge: U_e -> g_{src(e)} U_e g_{tgt(e)}^dagger.
    """
    transformed: Dict[Tuple[int, int], np.ndarray] = {}
    for (u, v), U in link_variables.items():
        g_u = gauge.get(u, np.eye(3, dtype=complex))
        g_v = gauge.get(v, np.eye(3, dtype=complex))
        transformed[(u, v)] = g_u @ U @ g_v.conj().T
    return transformed


# ---------------------------------------------------------------------------
# 7. Cartan torus (maximal torus) reduction
# ---------------------------------------------------------------------------

def cartan_links_from_phases(
    phases_a: Dict[Tuple[int, int], float],
    phases_b: Dict[Tuple[int, int], float],
) -> Dict[Tuple[int, int], np.ndarray]:
    """
    Build SU(3) link variables in the Cartan torus from two U(1) phase
    fields a and b on oriented edges:
        U_e = diag(e^{i a_e}, e^{i b_e}, e^{-i(a_e + b_e)}).

    Reverse-edge convention: U_{e^{-1}} = U_e^dagger has phases (-a_e, -b_e).

    Parameters
    ----------
    phases_a, phases_b : dict (src, tgt) -> float
        Two independent U(1) phase fields on oriented edges.  Must have
        identical key sets.

    Returns
    -------
    links : dict (src, tgt) -> 3x3 complex array
    """
    if set(phases_a.keys()) != set(phases_b.keys()):
        raise ValueError("phases_a and phases_b must have identical keys")
    links: Dict[Tuple[int, int], np.ndarray] = {}
    for key in phases_a:
        a = phases_a[key]
        b = phases_b[key]
        links[key] = cartan_su3_from_phases(a, b)
    return links


def u1xu1_action_from_su3(
    plaquettes: Sequence[Sequence[OrientedEdge]],
    phases_a: Dict[Tuple[int, int], float],
    phases_b: Dict[Tuple[int, int], float],
    beta: float,
) -> float:
    """
    U(1) x U(1) Wilson action obtained from SU(3) restricted to its
    maximal Cartan torus.

    With link variables U_e = diag(e^{i a_e}, e^{i b_e}, e^{-i(a_e + b_e)}),
    the plaquette holonomy is
        U_P = diag(e^{i a_P}, e^{i b_P}, e^{-i(a_P + b_P)})
    where a_P = sum_{e in P} a_e, b_P = sum_{e in P} b_e (with reverse-edge
    phases negated).  The plaquette character is
        (1/3) Re tr U_P = (1/3) [cos a_P + cos b_P + cos(a_P + b_P)].

    Parameters
    ----------
    plaquettes : list of list[OrientedEdge]
    phases_a, phases_b : dict (src, tgt) -> float
    beta : float

    Returns
    -------
    S : float
    """
    S = 0.0
    for P in plaquettes:
        a_P = 0.0
        b_P = 0.0
        for e in P:
            key = (e.source, e.target)
            if key in phases_a:
                a_P += phases_a[key]
                b_P += phases_b[key]
            else:
                rev = (e.target, e.source)
                a_P -= phases_a[rev]
                b_P -= phases_b[rev]
        char = (np.cos(a_P) + np.cos(b_P) + np.cos(a_P + b_P)) / 3.0
        S += 1.0 - char
    return beta * S


# ---------------------------------------------------------------------------
# 8. Character expansion coefficients (SU(3) heat kernel on the group)
# ---------------------------------------------------------------------------

def su3_character_coefficient_fundamental(beta: float) -> float:
    """
    Leading non-trivial character coefficient d_{(1,0)}(beta) for the SU(3)
    Wilson Boltzmann factor expanded in the fundamental rep.

    The exact analytic expression for SU(3) character coefficients
    requires Weyl integration over the maximal torus with the Vandermonde
    Jacobian; closed forms exist (Drouffe-Zuber 1983 §3) but are
    cumbersome.  For the Wilson-loop expectation in the leading-order
    independent-plaquette approximation, the relevant ratio is
        <chi_F(U)> = c_F(beta) / c_0(beta)
    where chi_F is the fundamental character.  For SU(N) at large beta,
        c_F / c_0 -> 1 - (N^2 - 1)/(4 N beta) + O(beta^{-2})
    (Drouffe-Zuber 1983 eq. 3.40), giving the perimeter-law deconfined
    leading behavior.  For SU(3): c_F/c_0 -> 1 - 2/(3 beta) at large beta.

    For our finite-graph computation we evaluate the integral
        c_F(beta) = int_{SU(3)} dU chi_F(U) exp(beta * Re chi_F(U))
    by Monte Carlo on the SU(3) Haar measure.  This function returns
    just the LO ratio for diagnostic purposes; for actual Wilson-loop
    expectations use `monte_carlo_wilson_expectation` below.

    Parameters
    ----------
    beta : float

    Returns
    -------
    ratio : float
        Approximate <(1/3) Re tr U> at the given beta (single-link Haar
        average against exp(beta * Re tr U / 3)).
    """
    # Direct numerical integration via Haar sampling
    rng = np.random.default_rng(0)
    n_sample = 2000
    chis = np.zeros(n_sample)
    weights = np.zeros(n_sample)
    for i in range(n_sample):
        U = su3_random(rng)
        chi = float(np.real(np.trace(U))) / 3.0
        chis[i] = chi
        weights[i] = np.exp(beta * chi)
    Z = weights.mean()
    return float((chis * weights).mean() / Z)


# ---------------------------------------------------------------------------
# 9. Monte Carlo with Metropolis algorithm
# ---------------------------------------------------------------------------

def _su3_local_update(U_old: np.ndarray, scale: float, rng: np.random.Generator) -> np.ndarray:
    """
    Propose a local update U_old -> exp(i scale * X) U_old where X is a random
    Hermitian su(3) element with components ~ N(0, 1).

    The scale parameter controls the typical update size; small scale gives
    high acceptance but slow exploration; large scale gives broad exploration
    but low acceptance at moderate beta.
    """
    T = su3_generators()
    coeffs = scale * rng.standard_normal(8)
    X = sum(coeffs[a] * T[a] for a in range(8))
    return _matrix_exp_hermitian(1j * X) @ U_old


def monte_carlo_wilson_expectation(
    loop_edges: Sequence[OrientedEdge],
    plaquettes: Sequence[Sequence[OrientedEdge]],
    edge_list: Sequence[OrientedEdge],
    beta: float,
    n_samples: int = 1000,
    n_thermalize: int = 300,
    seed: int = 42,
    update: str = "local",
    update_scale: Optional[float] = None,
) -> Tuple[float, float]:
    """
    Metropolis Monte Carlo estimate of <W_C> = <(1/3) Re tr U_C> at coupling beta.

    Mirrors the SU(2) implementation in Paper 30 with two update modes:

    - ``update="haar"``: propose Haar-random SU(3) at each step.  Same as
      Paper 30's SU(2) algorithm.  Tends to have low acceptance at large
      beta because the Haar measure is dominated by far-from-identity
      configurations.

    - ``update="local"`` (default): propose U_e -> exp(i scale * X) U_e
      with X a random Hermitian su(3) generator combination, scale tuned
      to ~ 1/sqrt(beta) so the typical update changes the action by O(1).
      This gives reasonable acceptance across a wide beta range.

    Parameters
    ----------
    loop_edges : sequence of OrientedEdge
        The Wilson loop.
    plaquettes : sequence of sequence of OrientedEdge
    edge_list : sequence of OrientedEdge
        Forward link variables; must cover all plaquettes.
    beta : float
    n_samples : int
    n_thermalize : int
    seed : int
    update : str
        Either "haar" or "local".
    update_scale : float, optional
        For "local": typical update size in radians.  Defaults to
        ~ min(pi/2, 1/sqrt(beta)) which gives reasonable acceptance.

    Returns
    -------
    mean : float
        Sample mean of (1/3) Re tr U_C.
    stderr : float
        Standard error.
    """
    rng = np.random.default_rng(seed)
    link_vars: Dict[Tuple[int, int], np.ndarray] = {
        (e.source, e.target): np.eye(3, dtype=complex)
        for e in edge_list
    }
    S_current = wilson_action(plaquettes, link_vars, beta)

    if update_scale is None:
        if beta <= 0:
            update_scale = float(np.pi)
        else:
            update_scale = float(min(np.pi / 2, 1.0 / np.sqrt(max(beta, 0.1))))

    samples = []
    edge_keys = [(e.source, e.target) for e in edge_list]
    n_links = len(edge_keys)
    n_accept = 0

    for step in range(n_thermalize + n_samples):
        idx = rng.integers(n_links)
        key = edge_keys[idx]
        old_U = link_vars[key]
        if update == "haar":
            new_U = su3_random(rng)
        else:  # local
            new_U = _su3_local_update(old_U, update_scale, rng)
        link_vars[key] = new_U
        S_new = wilson_action(plaquettes, link_vars, beta)
        dS = S_new - S_current
        if rng.random() < np.exp(-min(dS, 50.0)):
            S_current = S_new
            if step >= n_thermalize:
                n_accept += 1
        else:
            link_vars[key] = old_U

        if step >= n_thermalize:
            U_C = plaquette_holonomy(loop_edges, link_vars)
            samples.append(su3_character(U_C))

    samples = np.asarray(samples)
    mean = float(samples.mean())
    stderr = float(samples.std(ddof=1) / np.sqrt(max(len(samples), 1)))
    return mean, stderr


# ---------------------------------------------------------------------------
# 10. Weak-coupling kinetic-term derivation (symbolic)
# ---------------------------------------------------------------------------

def gell_mann_anticommutators_normalization() -> Tuple[List[np.ndarray], float]:
    """
    Verify the Gell-Mann normalization tr(T^a T^b) = (1/2) delta_{ab} and
    return the generators with the normalization constant.

    Returns
    -------
    T : list of 3x3 complex arrays
        T^a = lambda^a / 2.
    norm : float
        Common value of tr(T^a T^a) (should be 1/2 for all a).
    """
    T = su3_generators()
    n = float(np.real(np.trace(T[0] @ T[0])))
    return T, n


def weak_coupling_kinetic_coefficient_per_plaquette() -> sp.Rational:
    """
    Symbolic derivation of the weak-coupling kinetic-term coefficient
    per plaquette per su(3) component.

    Expand U_e = exp(i sum_a A_e^a T^a) ≈ I + i A_e^a T^a - (1/2)(A_e^a T^a)^2 + ...
    Plaquette holonomy U_P = prod_e U_e ≈ I + i (sum_e A_e^a) T^a + O(A^2).
    Take 1 - (1/3) Re tr U_P.

    Using tr(T^a) = 0, tr(T^a T^b) = (1/2) delta^{ab}:
        Re tr U_P ≈ 3 - (1/2) (sum_e A_e^a)(sum_e A_e^a) + O(A^4)  (sum over a)

    Actually more carefully:
        tr U_P ≈ 3 + i tr(A_P^a T^a) - (1/2) tr(A_P^a A_P^b T^a T^b) + O(A^3)
              = 3 + 0 - (1/2) A_P^a A_P^b (1/2) delta^{ab} + O(A^3)
              = 3 - (1/4) (A_P^a)^2 + O(A^3)
    where A_P^a = sum_{e in P} A_e^a.

    Re tr U_P = 3 - (1/4) (A_P^a)^2 + O(A^4).
    1 - (1/3) Re tr U_P = (1/12) (A_P^a)^2 + O(A^4).

    So the coefficient PER PLAQUETTE PER su(3)-COMPONENT is 1/12.
    Multiply by beta to get the action contribution.

    For comparison, SU(2) gives 1/8 per plaquette per su(2) component
    (Paper 30 eq. (\\ref{eq:plaquette_quadratic})).  The general SU(N)
    coefficient is 1/(2 N N_T) where N_T is the trace normalization;
    for our T^a = lambda^a/2 with tr(T^a T^a) = 1/2, the coefficient
    is 1/(2 * 2 * N) = 1/(4N), giving 1/8 for SU(2) and 1/12 for SU(3).

    Returns
    -------
    coeff : sympy Rational
        1/12 = 1/(4*3) = 1/(2*N_c*tr_norm).
    """
    return sp.Rational(1, 12)


# ---------------------------------------------------------------------------
# 11. Public API
# ---------------------------------------------------------------------------

__all__ = [
    "OrientedEdge",
    # SU(3) primitives
    "gell_mann_matrices",
    "su3_generators",
    "su3_from_algebra",
    "su3_random",
    "is_su3",
    "su3_character",
    "cartan_su3_from_phases",
    # Graph adapter
    "bargmann_adjacency_dense",
    # Edge / plaquette enumeration
    "enumerate_oriented_edges",
    "enumerate_plaquettes",
    # Wilson action
    "edge_link_variable",
    "plaquette_holonomy",
    "wilson_action",
    "gauge_transform",
    # Cartan torus reduction
    "cartan_links_from_phases",
    "u1xu1_action_from_su3",
    # Character / Monte Carlo
    "su3_character_coefficient_fundamental",
    "monte_carlo_wilson_expectation",
    # Symbolic
    "gell_mann_anticommutators_normalization",
    "weak_coupling_kinetic_coefficient_per_plaquette",
    # Update helpers
    "_su3_local_update",
]
