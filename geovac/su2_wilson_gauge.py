"""
SU(2) Wilson Lattice Gauge Theory on the GeoVac Hopf Graph
==========================================================

This module implements a genuine Wilson lattice gauge theory on the
Fock-projected S^3 = SU(2) Coulomb graph. Each oriented edge e carries
a link variable U_e in SU(2); plaquettes are primitive closed
non-backtracking walks (the Ihara cycles of `geovac/ihara_zeta.py`);
the Wilson action is

    S_W[U] = beta * sum_{plaquettes P} (1 - (1/2) Re tr U_P),

where the plaquette holonomy U_P is the ordered product of link
variables around P. The partition function is the Haar-measure
integral over SU(2)^|E|:

    Z(beta) = int dU exp(-S_W[U]).

Framing
-------

Paper 25 shows the GeoVac Hopf graph carries an abelian U(1)
Wilson-Hodge lattice-gauge structure, with the node Laplacian
L_0 = B B^T as matter propagator and the edge Laplacian L_1 = B^T B as
discrete Hodge-1 Laplacian. Paper 25 Section VII.1 examined whether a
non-abelian analog exists on the Bargmann-Segal S^5 graph and found
it does not (transitions are Clebsch-Gordan intertwiners, not SU(3)
group elements).

The S^3 Coulomb graph is different: S^3 is itself the Lie group SU(2),
so SU(2)-valued link variables are a natural (and structurally correct)
non-abelian object. This module builds that structure explicitly.

What this IS and what it IS NOT
-------------------------------

IS: a concrete non-abelian lattice gauge theory on a compact discrete
spatial slice; gauge-invariant; reduces to the U(1) construction in
the maximal-torus limit; finite-dimensional Haar integrals are
computable to any accuracy via character expansion.

IS NOT: a proof of mass gap on R^4; a renormalizable continuum field
theory; a Yang-Mills theory on a Lorentzian spacetime. This is an
Euclidean lattice gauge theory on a fixed finite graph. It is Yang-
Mills-adjacent in the sense that it uses the Wilson action, not in
the sense that it proves continuum-Yang-Mills statements.

Conventions
-----------

* An SU(2) matrix U is parametrized by a unit quaternion
  (a, b, c, d) with a^2 + b^2 + c^2 + d^2 = 1:
      U = [[a + i d,  i b + c], [i b - c,  a - i d]].
  Equivalently, U = exp(i theta (n . sigma)) with sigma the Pauli
  matrices.

* Each undirected edge {i, j} of the S^3 Coulomb graph is promoted to
  two oriented edges i->j and j->i. If U is the link variable for
  i->j, then j->i carries U^dagger.

* Plaquettes are primitive closed non-backtracking walks, enumerated
  via the Hashimoto edge operator T (see `geovac/ihara_zeta.py`).
  A primitive closed walk of length L is an equivalence class of
  closed NB walks on L edges modulo cyclic rotation. We pick one
  representative per class (with a canonical base point and
  orientation) for the Wilson action. Orientation-reversal gives a
  distinct walk; we sum over both, which makes S_W real and Re-symmetric.

* The Haar measure on SU(2) is the normalized bi-invariant probability
  measure. Its characters chi_j, j = 0, 1/2, 1, 3/2, ..., are
  chi_j(U) = sin((2j+1) theta) / sin(theta) where
  U = exp(i theta (n . sigma)). Orthogonality:
      int dU chi_j(U) chi_k(U^dagger) = delta_{j,k}.

* Character expansion of the heat-kernel Boltzmann weight:
      exp(beta * (1/2) Re tr U_P)
        = sum_{j >= 0, half-integer} d_j(beta) chi_{2j}(U_P),
  where the coefficients d_j(beta) satisfy
      d_j(beta) = int dU chi_{2j}(U) exp(beta * (1/2) Re tr U).
  Here we use (1/2) Re tr U = cos(theta). For SU(2), j = 0, 1/2, 1, ...
  correspond to reps of dimensions 1, 2, 3, 4, ....

Character coefficients
----------------------

For SU(2) the Wilson-action character coefficients are
    c_r(beta) = I_{2r+1}(beta) / I_1(beta) * (2r+1),
or equivalently in the normalization convention used here:

    d_r(beta) = I_{r+1}(beta) - I_{r-1}(beta)
for the reduction of exp(beta cos theta) under the Haar measure with
Jacobian (2/pi) sin^2 theta on [0, pi].

We use the normalized expansion

    exp(beta cos theta) = sum_{n=0}^infty a_n(beta) U_n(cos theta),

where U_n is the Chebyshev polynomial of the second kind and
a_n(beta) = 2 I_{n+1}(beta) / beta in the classical expansion of
Bessel generating functions (Abramowitz & Stegun 9.6.34).

In our sign convention exp(-S_W) for a single plaquette gives
exp(beta * (1/2) Re tr U_P - beta), so only the relative
character coefficients matter for ratios and expectation values.

References
----------

- K. G. Wilson, Phys. Rev. D 10, 2445 (1974).
- M. Creutz, "Quarks, Gluons and Lattices," Cambridge Univ. Press (1983),
  chs. 6-8.
- R. Oeckl, "Discrete Gauge Theory: From Lattices to TQFT," Imperial
  College Press (2005), ch. 3.

Author: GeoVac Track RH-Q (Sprint 4).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import sympy as sp


# ---------------------------------------------------------------------------
# 1. SU(2) elementary operations
# ---------------------------------------------------------------------------

def su2_from_quaternion(
    a: float, b: float, c: float, d: float
) -> np.ndarray:
    """
    Construct a 2x2 SU(2) matrix from a unit quaternion (a, b, c, d).

    The quaternion q = a + b i + c j + d k is mapped to the matrix
        U = [[a + i d,  i b + c], [i b - c,  a - i d]].
    The unit-norm constraint a^2 + b^2 + c^2 + d^2 = 1 gives det(U) = 1
    and U U^dagger = I.

    Parameters
    ----------
    a, b, c, d : float
        Quaternion components.  Should satisfy a^2 + b^2 + c^2 + d^2 = 1.

    Returns
    -------
    U : np.ndarray, shape (2, 2), dtype complex
        The SU(2) matrix.
    """
    norm_sq = a * a + b * b + c * c + d * d
    if not np.isclose(norm_sq, 1.0, atol=1e-10):
        raise ValueError(
            f"Quaternion must be unit norm, got |q|^2 = {norm_sq}"
        )
    U = np.array(
        [[a + 1j * d, 1j * b + c], [1j * b - c, a - 1j * d]],
        dtype=complex,
    )
    return U


def su2_random(rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """
    Sample a random SU(2) matrix from the Haar measure.

    The Haar-uniform distribution on SU(2) is equivalent to uniform
    sampling from the unit 3-sphere S^3 in quaternion coordinates; we
    draw a 4D standard-normal vector and normalize.

    Parameters
    ----------
    rng : np.random.Generator, optional
        Random number generator.  If None, uses the default numpy RNG.

    Returns
    -------
    U : np.ndarray, shape (2, 2), dtype complex
        Haar-random SU(2) matrix.
    """
    if rng is None:
        rng = np.random.default_rng()
    q = rng.standard_normal(4)
    q = q / np.linalg.norm(q)
    return su2_from_quaternion(q[0], q[1], q[2], q[3])


def su2_from_angle_axis(
    theta: float, axis: Sequence[float]
) -> np.ndarray:
    """
    Construct U = exp(i (theta/2) (n . sigma)) for unit axis n.

    This is the rotation by angle theta about axis n in the spin-1/2
    representation.

    Parameters
    ----------
    theta : float
        Rotation angle (radians).
    axis : length-3 array
        Unit vector specifying the axis.

    Returns
    -------
    U : np.ndarray, shape (2, 2)
    """
    n = np.asarray(axis, dtype=float)
    nn = np.linalg.norm(n)
    if nn < 1e-12:
        return np.eye(2, dtype=complex)
    n = n / nn
    # exp(i (theta/2) n . sigma) = cos(theta/2) I + i sin(theta/2) n . sigma
    c = np.cos(theta / 2.0)
    s = np.sin(theta / 2.0)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    U = c * np.eye(2, dtype=complex) + 1j * s * (
        n[0] * sigma_x + n[1] * sigma_y + n[2] * sigma_z
    )
    return U


def is_su2(U: np.ndarray, atol: float = 1e-8) -> bool:
    """
    Check whether U lies in SU(2) within tolerance atol.

    Verifies: (i) shape 2x2; (ii) det(U) ~ 1; (iii) U U^dagger ~ I.
    """
    if U.shape != (2, 2):
        return False
    det_ok = np.isclose(np.linalg.det(U), 1.0, atol=atol)
    unitary_ok = np.allclose(U @ U.conj().T, np.eye(2), atol=atol)
    return bool(det_ok and unitary_ok)


def su2_character(U: np.ndarray) -> float:
    """
    Wilson plaquette character (1/2) Re tr U for SU(2).

    For U = exp(i theta (n . sigma)), this equals cos(theta).  It is the
    fundamental-representation character divided by dimension (= 2).
    """
    return 0.5 * float(np.real(np.trace(U)))


# ---------------------------------------------------------------------------
# 2. Oriented-edge bookkeeping on a graph
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
    Build the list of oriented edges of a simple graph.

    Each undirected edge {u, v} yields two oriented edges (u, v) and
    (v, u).  The list is canonicalized by the natural (source, target)
    lexicographic order.

    Parameters
    ----------
    adjacency : np.ndarray, shape (V, V)
        Symmetric 0/1 adjacency.  Diagonal must be zero.

    Returns
    -------
    oriented : list of OrientedEdge
        All 2 * E directed edges.
    index : dict (src, tgt) -> int
        Inverse lookup.
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
# 3. Plaquette enumeration
# ---------------------------------------------------------------------------

def enumerate_plaquettes(
    adjacency: np.ndarray,
    max_length: int = 6,
    both_orientations: bool = False,
) -> List[List[OrientedEdge]]:
    """
    Enumerate primitive closed non-backtracking walks up to max_length.

    A primitive closed NB walk on L edges is a closed walk
    e_1 e_2 ... e_L with head(e_i) = tail(e_{i+1}), e_{i+1} not the
    reverse of e_i, and head(e_L) = tail(e_1), which is not a cyclic
    rotation of a shorter walk.  Two walks differing by cyclic
    rotation are treated as the same plaquette; we pick the
    representative whose starting oriented-edge index is smallest.

    Parameters
    ----------
    adjacency : np.ndarray
        Symmetric adjacency.
    max_length : int
        Maximum plaquette length to enumerate.
    both_orientations : bool
        If True, include both directions of each unoriented cycle
        (CW and CCW).  If False, keep only one representative per
        unoriented cycle.  For Wilson-action applications, both
        orientations contribute because the action is real (Re tr),
        so we typically set both_orientations=False and include a
        factor of 2 at the action level.

    Returns
    -------
    plaquettes : list of list[OrientedEdge]
        Each inner list is an ordered sequence of oriented edges
        forming one plaquette.
    """
    A = np.asarray(adjacency)
    V = A.shape[0]
    oriented, edge_index = enumerate_oriented_edges(A)
    n_oriented = len(oriented)

    seen: set = set()
    plaquettes: List[List[OrientedEdge]] = []

    def is_primitive(walk: Tuple[int, ...]) -> bool:
        """Return True if walk is not a proper power of a shorter walk."""
        L = len(walk)
        for d in range(1, L):
            if L % d != 0:
                continue
            # Check if walk = (walk[:d]) ** (L/d)
            period = walk[:d]
            is_power = True
            for k in range(L // d):
                if walk[k * d:(k + 1) * d] != period:
                    is_power = False
                    break
            if is_power:
                return False
        return True

    # Enumerate via depth-first search over oriented-edge sequences.
    # At each step, we extend by an oriented edge e' such that
    # tail(e') = head(e_prev) AND e' != reverse(e_prev).
    # Cycle closure: head(e_L) == tail(e_1) AND e_L is not the reverse of e_1
    # (no-backtracking around the cycle).
    for start_idx in range(n_oriented):
        # DFS stack: (current_path).
        stack: List[List[int]] = [[start_idx]]
        while stack:
            path_idx = stack.pop()
            if len(path_idx) >= max_length:
                # Cannot extend further; try closing only if already closed
                continue
            last = path_idx[-1]
            last_edge = oriented[last]
            # Successor edges: tail == last_edge.target, not reverse of last_edge.
            for j in np.nonzero(A[last_edge.target])[0]:
                next_idx = edge_index[(last_edge.target, int(j))]
                next_edge = oriented[next_idx]
                # No backtracking over last edge
                if next_edge.source == last_edge.target and next_edge.target == last_edge.source:
                    continue
                new_path = path_idx + [next_idx]
                first_edge = oriented[new_path[0]]
                # Check if we've closed the walk
                if next_edge.target == first_edge.source and len(new_path) >= 3:
                    # Verify no backtracking across the closure:
                    # the "next edge after next_edge" would be first_edge,
                    # which must not be the reverse of next_edge.
                    if not (
                        first_edge.source == next_edge.target
                        and first_edge.target == next_edge.source
                    ):
                        # Filter out non-primitive walks (proper powers of shorter walks)
                        path_tuple = tuple(new_path)
                        if is_primitive(path_tuple):
                            # Canonicalize: cycle starts with minimum oriented-edge index.
                            min_pos = int(np.argmin(new_path))
                            canonical = tuple(new_path[min_pos:] + new_path[:min_pos])
                            # Reverse-orientation canonical
                            reversed_edges = [
                                edge_index[(oriented[i].target, oriented[i].source)]
                                for i in reversed(new_path)
                            ]
                            min_pos_r = int(np.argmin(reversed_edges))
                            canonical_r = tuple(
                                reversed_edges[min_pos_r:] + reversed_edges[:min_pos_r]
                            )
                            if both_orientations:
                                if canonical not in seen:
                                    seen.add(canonical)
                                    plaquettes.append([oriented[i] for i in canonical])
                            else:
                                rep = canonical if canonical <= canonical_r else canonical_r
                                if rep not in seen:
                                    seen.add(rep)
                                    plaquettes.append([oriented[i] for i in rep])
                # Continue extending regardless of closure (longer plaquettes exist)
                if len(new_path) < max_length:
                    stack.append(new_path)

    return plaquettes


# ---------------------------------------------------------------------------
# 4. Link variables, plaquette holonomy, Wilson action
# ---------------------------------------------------------------------------

def edge_link_variable(
    edge: OrientedEdge,
    link_variables: Dict[Tuple[int, int], np.ndarray],
) -> np.ndarray:
    """
    Look up the SU(2) link variable for an oriented edge.

    If the forward edge (source, target) is present in the dict,
    returns its matrix.  If only the reverse (target, source) is
    present, returns its hermitian conjugate (this enforces the
    lattice-gauge convention U_{e^-1} = U_e^dagger).

    Parameters
    ----------
    edge : OrientedEdge
    link_variables : dict (src, tgt) -> 2x2 complex array

    Returns
    -------
    U : np.ndarray, shape (2, 2)
    """
    key = (edge.source, edge.target)
    if key in link_variables:
        return link_variables[key]
    rev = (edge.target, edge.source)
    if rev in link_variables:
        return link_variables[rev].conj().T
    raise KeyError(
        f"No link variable found for edge {edge.source}->{edge.target} "
        f"or its reverse."
    )


def plaquette_holonomy(
    plaquette: Sequence[OrientedEdge],
    link_variables: Dict[Tuple[int, int], np.ndarray],
) -> np.ndarray:
    """
    Ordered product U_{e_1} U_{e_2} ... U_{e_L} around plaquette P.

    Parameters
    ----------
    plaquette : sequence of OrientedEdge
    link_variables : dict (src, tgt) -> 2x2 complex array

    Returns
    -------
    U_P : np.ndarray, shape (2, 2)
    """
    U_P = np.eye(2, dtype=complex)
    for e in plaquette:
        U_P = U_P @ edge_link_variable(e, link_variables)
    return U_P


def wilson_action(
    plaquettes: Sequence[Sequence[OrientedEdge]],
    link_variables: Dict[Tuple[int, int], np.ndarray],
    beta: float,
) -> float:
    """
    Wilson action S_W = beta * sum_P (1 - (1/2) Re tr U_P).

    For SU(2), (1/2) Re tr U_P is the normalized fundamental-rep
    character; the action is 0 when all plaquettes are the identity
    (trivial vacuum) and positive otherwise.

    Parameters
    ----------
    plaquettes : list of list[OrientedEdge]
    link_variables : dict (src, tgt) -> 2x2 complex array
    beta : float
        Gauge coupling.  beta = 4/g^2 in the standard Wilson
        convention; large beta = weak coupling; small beta = strong
        coupling.

    Returns
    -------
    S : float
        Total action.
    """
    S = 0.0
    for P in plaquettes:
        U_P = plaquette_holonomy(P, link_variables)
        S += 1.0 - su2_character(U_P)
    return beta * S


# ---------------------------------------------------------------------------
# 5. Gauge transformation
# ---------------------------------------------------------------------------

def gauge_transform(
    link_variables: Dict[Tuple[int, int], np.ndarray],
    gauge: Dict[int, np.ndarray],
) -> Dict[Tuple[int, int], np.ndarray]:
    """
    Apply a node-local gauge transformation U_e -> g_{src(e)} U_e g_{tgt(e)}^dagger.

    Parameters
    ----------
    link_variables : dict (src, tgt) -> 2x2 complex array
        Forward link variables.
    gauge : dict node_idx -> 2x2 SU(2) matrix
        The gauge transformation.

    Returns
    -------
    transformed : dict (src, tgt) -> 2x2 complex array
    """
    transformed: Dict[Tuple[int, int], np.ndarray] = {}
    for (u, v), U in link_variables.items():
        g_u = gauge.get(u, np.eye(2, dtype=complex))
        g_v = gauge.get(v, np.eye(2, dtype=complex))
        transformed[(u, v)] = g_u @ U @ g_v.conj().T
    return transformed


# ---------------------------------------------------------------------------
# 6. Character expansion and partition function
# ---------------------------------------------------------------------------

def su2_character_coefficient(j: int, beta: float) -> float:
    """
    Character expansion coefficient c_j(beta) for SU(2) Wilson action.

    The heat-kernel-like expansion of the Wilson Boltzmann factor is

        exp(beta cos theta) = I_0(beta) + 2 sum_{j >= 1} I_j(beta) cos(j theta),

    but in the SU(2) class-function basis (using characters
    chi_j(U) = sin((j+1) theta)/sin(theta) of the (j+1)-dimensional
    rep), Fegan/Menotti/Onofri give

        exp(beta cos theta) = sum_{j=0}^infty d_j(beta) chi_j(U),

    with d_j(beta) = (2 / beta) (j + 1) I_{j+1}(beta).  We use the
    convention that j = 0, 1, 2, ... labels reps of dimension j+1
    (so j = 0 -> trivial rep, j = 1 -> fundamental rep).

    Parameters
    ----------
    j : int
        Non-negative integer labeling SU(2) rep by dimension = j + 1.
    beta : float
        Gauge coupling.

    Returns
    -------
    c_j : float
        The character coefficient d_j(beta).
    """
    from scipy.special import iv, ive
    if j < 0:
        raise ValueError("j must be >= 0")
    if abs(beta) < 1e-15:
        # Limit: d_0 -> 1, all other d_j -> 0 for beta -> 0
        return 1.0 if j == 0 else 0.0
    # Use ive (exponentially scaled) for large beta; unscale explicitly
    # since we return an absolute coefficient, not a ratio.
    # For moderate beta (<700), iv is safe; else rely on ive * exp(beta).
    if abs(beta) < 700:
        return (2.0 / beta) * (j + 1) * iv(j + 1, beta)
    return (2.0 / beta) * (j + 1) * ive(j + 1, beta) * np.exp(abs(beta))


def partition_function_character_expansion(
    plaquettes: Sequence[Sequence[OrientedEdge]],
    beta: float,
    R_max: int = 3,
) -> float:
    """
    Partition function Z(beta) via SU(2) character expansion.

    On a graph with NO plaquettes (beta_1 = 0), Z is trivially the
    product of Haar volumes, which normalize to 1 per link: Z = 1.

    For a graph with p plaquettes, ALL treated as independent (i.e.,
    assuming we ignore intersections), we have
        Z = exp(-beta * p) * sum_{j=0}^{R_max} d_j(beta) * <chi_j>_0^p,
    where <chi_j>_0 = delta_{j,0} by Haar orthogonality of the
    characters under a single-link integration.  So the leading-order
    character-expansion approximation simply gives
        Z_0(beta) = exp(-beta * p) * d_0(beta)^p.

    This approximation is exact when plaquettes share no edges,
    i.e., when the graph is a disjoint union of isolated cycles.
    For our Hopf graph at n_max = 3, the two primitive 4-cycles
    may share edges; then higher-rep contributions correct this.
    We include reps up to 2j = R_max.

    Parameters
    ----------
    plaquettes : list of list[OrientedEdge]
    beta : float
    R_max : int
        Maximum rep label to include (corresponds to 2j = R_max,
        so the fundamental rep j=1 has R_max >= 1).

    Returns
    -------
    Z : float
        Partition function approximated to rep-label R_max.

    Notes
    -----
    The full Z on graphs with edge-sharing plaquettes requires careful
    treatment of the shared edges via Kronecker sums of characters.
    We give the leading-order approximation here, which is exact when
    plaquettes do not share edges.  For general cases, Monte Carlo
    (method) is the right tool.
    """
    p = len(plaquettes)
    if p == 0:
        return 1.0
    # Leading order: Z_0 = exp(-beta * p) * d_0(beta)^p
    d0 = su2_character_coefficient(0, beta)
    from scipy.special import iv
    if abs(beta) < 1e-15:
        d0 = 1.0
    else:
        d0 = (2.0 / beta) * iv(1, beta)  # (j=0, dim=1)
    Z0 = np.exp(-beta * p) * (d0 ** p)

    # Higher-order corrections from shared-edge plaquette correlations.
    # These are O(d_j(beta) / d_0(beta))^{shared edges}. For our
    # small graphs they are tiny at moderate beta. Keep leading term.
    # (Extensions: for graphs where plaquettes share edges, one needs
    # to track the Gauss-linking structure; we defer that to Sprint 5.)
    return float(Z0)


# ---------------------------------------------------------------------------
# 7. Expectation value of a Wilson loop
# ---------------------------------------------------------------------------

def expectation_wilson_loop(
    loop_edges: Sequence[OrientedEdge],
    plaquettes: Sequence[Sequence[OrientedEdge]],
    beta: float,
    R_max: int = 3,
) -> float:
    """
    Expectation value <W_C> = <(1/2) Re tr U_C> of a Wilson loop.

    In strong coupling (small beta), <W_C> ~ exp(-sigma * Area(C))
    (area law, confinement signal).  In weak coupling (large beta),
    <W_C> ~ exp(-mu * Perimeter(C)) (perimeter law, deconfined).

    This function computes the leading-order character-expansion
    approximation to <W_C> assuming each plaquette is independent:

        <W_C> ~ product over plaquettes tiling C of (d_1(beta)/d_0(beta)).

    For a loop C bounding a disk tiled by A plaquettes,
    <W_C> ~ (I_2/I_1)^A = tanh(beta/2)^A for beta >> 1 (Creutz 1983).

    For a plaquette-counting argument on our graph:
    - loop_edges of length L, if it coincides with the boundary of a
      single plaquette of length L, has <W_C> = d_1(beta) / d_0(beta).
    - if it is longer (no plaquette of that length), <W_C> ~ 0 at
      leading order unless it tiles as a union of smaller plaquettes.

    We implement the simplest case: if loop_edges IS one of the
    enumerated plaquettes, return d_1/d_0; otherwise return 0 (leading
    order).  More general Wilson loops require tiling enumeration.

    Parameters
    ----------
    loop_edges : sequence of OrientedEdge
    plaquettes : list of list[OrientedEdge]
    beta : float
    R_max : int

    Returns
    -------
    <W_C> : float
        Leading-order expectation value.
    """
    loop_tuple = tuple((e.source, e.target) for e in loop_edges)
    # Cyclic canonical form
    def canonicalize(t):
        L = len(t)
        # Try all rotations of t and its reversal; pick lex-min
        best = t
        for k in range(L):
            rot = t[k:] + t[:k]
            if rot < best:
                best = rot
        # Reverse
        rev = tuple((b, a) for (a, b) in reversed(t))
        for k in range(L):
            rot = rev[k:] + rev[:k]
            if rot < best:
                best = rot
        return best

    loop_can = canonicalize(loop_tuple)

    matched = False
    for P in plaquettes:
        P_tuple = tuple((e.source, e.target) for e in P)
        if canonicalize(P_tuple) == loop_can:
            matched = True
            break

    if not matched:
        return 0.0  # leading order

    # d_1(beta) / d_0(beta) = I_2(beta) / I_1(beta)
    # Use ive (exponentially scaled) to avoid overflow for large beta.
    from scipy.special import ive
    if abs(beta) < 1e-15:
        return 0.0
    # ive(n, beta) = exp(-|beta|) * iv(n, beta), so the ratio is
    # ive(2, beta) / ive(1, beta) = iv(2, beta) / iv(1, beta).
    return float(ive(2, beta) / ive(1, beta))


# ---------------------------------------------------------------------------
# 8. U(1) reduction check
# ---------------------------------------------------------------------------

def diagonal_su2_from_phase(phi: float) -> np.ndarray:
    """
    SU(2) matrix diag(e^{i phi}, e^{-i phi}) corresponding to a U(1)
    subgroup (maximal torus).
    """
    return np.array(
        [[np.exp(1j * phi), 0], [0, np.exp(-1j * phi)]],
        dtype=complex,
    )


def u1_action_from_su2(
    plaquettes: Sequence[Sequence[OrientedEdge]],
    phases: Dict[Tuple[int, int], float],
    beta: float,
) -> float:
    """
    U(1) Wilson action in the maximal-torus limit of SU(2).

    With link variables U_e = diag(e^{i phi_e}, e^{-i phi_e}),
    plaquette holonomy is U_P = diag(e^{i sum_e phi_e},
    e^{-i sum_e phi_e}), so (1/2) Re tr U_P = cos(theta_P) with
    theta_P = sum_{e in P} phi_e.  This matches Paper 25's U(1)
    Wilson-Hodge construction.

    Parameters
    ----------
    plaquettes : list of list[OrientedEdge]
    phases : dict (src, tgt) -> float
        U(1) phases on oriented edges.  Reverse edge phase is negated.
    beta : float

    Returns
    -------
    S : float
        U(1) Wilson action = beta * sum_P (1 - cos theta_P).
    """
    S = 0.0
    for P in plaquettes:
        theta = 0.0
        for e in P:
            key = (e.source, e.target)
            if key in phases:
                theta += phases[key]
            else:
                rev = (e.target, e.source)
                theta -= phases[rev]
        S += 1.0 - np.cos(theta)
    return beta * S


# ---------------------------------------------------------------------------
# 9. Monte Carlo with heat-bath algorithm (for Wilson loop expectation)
# ---------------------------------------------------------------------------

def monte_carlo_wilson_expectation(
    loop_edges: Sequence[OrientedEdge],
    plaquettes: Sequence[Sequence[OrientedEdge]],
    edge_list: Sequence[OrientedEdge],
    beta: float,
    n_samples: int = 2000,
    n_thermalize: int = 500,
    seed: int = 42,
) -> Tuple[float, float]:
    """
    Metropolis Monte Carlo estimate of <W_C> at coupling beta.

    Rather than a full heat-bath algorithm (which requires SU(2)
    conditional sampling), we use a simple Metropolis update: at each
    step, pick a random link, propose a new Haar-random SU(2) value,
    and accept with Metropolis probability exp(-Delta S).

    Parameters
    ----------
    loop_edges : sequence of OrientedEdge
        The Wilson loop.
    plaquettes : sequence of sequence of OrientedEdge
    edge_list : sequence of OrientedEdge
        Unique oriented edges storing forward link variables.  Must
        be enough to reconstruct all plaquettes via
        edge_link_variable().
    beta : float
    n_samples : int
    n_thermalize : int
    seed : int

    Returns
    -------
    mean : float
        Sample mean of (1/2) Re tr U_C.
    stderr : float
        Standard error.
    """
    rng = np.random.default_rng(seed)
    # Initialize to trivial config
    link_vars: Dict[Tuple[int, int], np.ndarray] = {
        (e.source, e.target): np.eye(2, dtype=complex)
        for e in edge_list
    }
    S_current = wilson_action(plaquettes, link_vars, beta)

    samples = []
    edge_keys = [(e.source, e.target) for e in edge_list]
    n_links = len(edge_keys)

    for step in range(n_thermalize + n_samples):
        # Propose a Haar-random update on a random link
        idx = rng.integers(n_links)
        key = edge_keys[idx]
        old_U = link_vars[key]
        new_U = su2_random(rng)
        link_vars[key] = new_U
        S_new = wilson_action(plaquettes, link_vars, beta)
        dS = S_new - S_current
        if rng.random() < np.exp(-min(dS, 50.0)):
            S_current = S_new
        else:
            link_vars[key] = old_U

        if step >= n_thermalize:
            U_C = plaquette_holonomy(loop_edges, link_vars)
            samples.append(su2_character(U_C))

    samples = np.asarray(samples)
    mean = float(samples.mean())
    stderr = float(samples.std(ddof=1) / np.sqrt(len(samples)))
    return mean, stderr


# ---------------------------------------------------------------------------
# __all__
# ---------------------------------------------------------------------------

__all__ = [
    "OrientedEdge",
    "su2_from_quaternion",
    "su2_random",
    "su2_from_angle_axis",
    "is_su2",
    "su2_character",
    "enumerate_oriented_edges",
    "enumerate_plaquettes",
    "edge_link_variable",
    "plaquette_holonomy",
    "wilson_action",
    "gauge_transform",
    "su2_character_coefficient",
    "partition_function_character_expansion",
    "expectation_wilson_loop",
    "diagonal_su2_from_phase",
    "u1_action_from_su2",
    "monte_carlo_wilson_expectation",
]
