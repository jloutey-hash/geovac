"""
Weighted Ihara zeta function (Track RH-K).
===========================================

This module implements the **edge-weighted** Ihara-Bass determinantal
formula (Mizuno-Sato 2004; Setyadi-Storm 2022) and the associated
weighted Hashimoto edge operator, then applies them to the Dirac-S^3
graph of Track RH-C with edge weights derived from the Breit-Pauli
spin-orbit coupling (proportional to alpha^2).

Weighted Ihara-Bass identity
-----------------------------
For a finite graph G = (V, E) with positive edge weights w : E -> R_{>0},
define the **weighted adjacency** A_w (symmetric, A_w[i,j] = w(edge ij))
and the **weighted generalized degree diagonal**

    Q_w[i,i] = (sum over j ~ i) w(edge ij)  -  1.

The weighted Ihara zeta is defined by the Euler product over primitive
closed walks C,

    zeta_G^w(s) = prod_{[C]} (1 - W(C) s^{length(C)})^{-1},

where W(C) = prod over edges e in C of w(e). The Mizuno-Sato
generalization of the Bass formula gives

    zeta_G^w(s)^{-1} = (1 - s^2)^{r - c}  *  det(I - s A_w + s^2 Q_w).

For w identically 1, A_w = A and Q_w = Q = diag(d_v - 1), and we recover
the ordinary Ihara-Bass formula (Track RH-A).

Weighted Hashimoto edge operator
---------------------------------
Let T_w be the (2E) x (2E) non-backtracking operator on oriented edges,

    (T_w)[(u->v), (x->y)] = sqrt(w(uv) * w(xy))   if v = x and y != u,
                         = 0                     otherwise.

Then as in the unweighted case (Setyadi-Storm 2022, Eq. 4.6)

    zeta_G^w(s)^{-1} = det(I - s T_w).

Edge-weight convention chosen for Track RH-K
---------------------------------------------
Edges (a, b) on the Dirac-S^3 graph (Rule A or B) are weighted by

    w(a, b) = 1  +  alpha^2 * omega(a, b),

where omega(a, b) is a positive, dimensionless, symmetric pair functional:

    omega(a, b) = | f_SO(a) + f_SO(b) |  / 2,

and f_SO(lab) is the Breit-Pauli spin-orbit diagonal for the Dirac
orbital, taken at the symbolic Z -> 1 (hydrogen-like) and with the
alpha^2 factor STRIPPED so that alpha^2 appears exactly once as the
perturbation parameter in w(a,b):

    f_SO(n, kappa)  =  -(kappa + 1) / [4 * n^3 * l * (l+1/2) * (l+1)]
                       * Z_eff^4   (= 1 here)

    (zero at l = 0 by Kramers cancellation)

Motivation:
 (i) w is symmetric and strictly positive for any real alpha (the
     symmetric mean of two signed SO diagonals, wrapped in | . |,
     is non-negative; adding the constant 1 keeps it > 0 even when
     omega = 0).
 (ii) alpha -> 0 gives w = 1 uniformly -> the unweighted zeta.
 (iii) alpha^2 enters linearly (not quadratically) in w, so the first
       non-trivial correction to the zeta polynomial is O(alpha^2).
 (iv) omega is "local" to the endpoints of each edge and does NOT
      require computing an off-diagonal SO matrix element (which
      vanishes anyway: H_SO is diagonal in the kappa-basis). We use
      the symmetric AVERAGE of the endpoint SO diagonals as a proxy
      for the pair's spin-orbit environment.
 (v) The factor of 1/4 in f_SO matches the denominator in
     Eq. (T2) of `geovac/spin_orbit.py`; we omit the Z-prefactor
     and treat Z_eff = 1 uniformly so that alpha^2 is the ONLY small
     parameter in the weights.

Alternative conventions considered and explicitly rejected:
 * Edge weight = |<H_SO>_a * <H_SO>_b|^(1/2): zero for any l=0
   endpoint, kills about half of the edges in Rule A and is
   artificially drastic.
 * Edge weight = alpha^2 * (additive SO): breaks the alpha -> 0
   sanity check (edges with omega = 0 would have w = 0, i.e. be
   SEVERED from the graph in the unweighted limit).
 * Off-diagonal <H_SO(a,b)>: H_SO is diagonal in (kappa, m_j), so
   this would be identically zero and carry no information.

Transcendental taxonomy (Paper 18)
-----------------------------------
The unweighted Ihara zeta is pi-free (Paper 29, integer-coefficient
polynomials). Adding alpha^2 edge weights introduces alpha as a new
transcendental seed, so the weighted zeta lives in the ring
Q(alpha^2)[s] — exactly the "spinor-intrinsic" tier of Paper 18
(ring R_sp, without the gamma content). This is the FIRST appearance
of a genuinely alpha^2-charged graph invariant in the framework;
previous alpha^2 content was confined to the H_SO diagonal.

References
----------
- Y. Mizuno and I. Sato, "Weighted zeta functions of graphs", J.
  Combin. Theory Ser. B 91 (2004) 169-183.
- H. M. Stark and A. A. Terras, "Zeta functions of finite graphs and
  coverings, III", Adv. Math. 208 (2007) 467-489.
- B. Setyadi and C. K. Storm, "Edge zeta functions", preprint (2022).
- A. Terras, "Zeta Functions of Graphs: A Stroll through the Garden"
  (Cambridge, 2011), Chap. 18.
- GeoVac Paper 18 (exchange constants taxonomy).
- GeoVac Paper 29 (unweighted Ihara Ramanujan observation).
- ``geovac.ihara_zeta`` — Track RH-A unweighted machinery.
- ``geovac.ihara_zeta_dirac`` — Track RH-C Dirac-S^3 graph builder.

Author: GeoVac Track RH-K (April 2026, Sprint 3 of the RH series).
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import sympy as sp

from geovac.dirac_matrix_elements import DiracLabel, kappa_to_l
from geovac.ihara_zeta import _count_components, _degree_sequence


# ---------------------------------------------------------------------------
# Module-level symbolic parameters
# ---------------------------------------------------------------------------

alpha_sym = sp.Symbol("alpha", positive=True)
s_sym = sp.Symbol("s")


# ---------------------------------------------------------------------------
# Generic weighted-graph machinery (symbolic sympy)
# ---------------------------------------------------------------------------

def _require_symmetric(A_w) -> sp.Matrix:
    """Coerce to a sympy square symmetric matrix, zero diagonal."""
    if not isinstance(A_w, sp.Matrix):
        M = sp.Matrix(A_w)
    else:
        M = A_w.copy()
    n, m = M.shape
    if n != m:
        raise ValueError("weighted adjacency must be square")
    for i in range(n):
        if M[i, i] != 0:
            raise ValueError("weighted adjacency must be loopless (zero diag)")
        for j in range(i + 1, n):
            if sp.simplify(M[i, j] - M[j, i]) != 0:
                raise ValueError("weighted adjacency must be symmetric")
    return M


def weighted_adjacency(
    edges: Sequence[Tuple[int, int]],
    weights: Sequence[sp.Expr],
    n_nodes: int,
) -> sp.Matrix:
    """Build the sympy-symbolic weighted adjacency from an edge list.

    Parameters
    ----------
    edges : sequence of (i, j) with i < j
        Undirected edges.
    weights : sequence of sympy Expr
        Positive weight for each edge; len must match ``edges``.
    n_nodes : int
        Total number of nodes (determines matrix size).

    Returns
    -------
    sympy.Matrix of shape (n_nodes, n_nodes), symmetric, zero diagonal.
    """
    if len(edges) != len(weights):
        raise ValueError("len(edges) must equal len(weights)")
    A = sp.zeros(n_nodes, n_nodes)
    for (i, j), w in zip(edges, weights):
        if i == j:
            raise ValueError(f"self-loop detected at node {i}")
        if i >= n_nodes or j >= n_nodes or i < 0 or j < 0:
            raise ValueError(f"edge ({i}, {j}) out of bounds")
        A[i, j] = w
        A[j, i] = w
    return A


def weighted_q_matrix(A_w: sp.Matrix) -> sp.Matrix:
    """Weighted excess-degree diagonal Q_w = diag(row-sum(A_w) - 1).

    This is the Mizuno-Sato generalization of the unweighted
    Q = diag(d_v - 1): for each node v, the diagonal entry is
        Q_w[v, v] = (sum of edge weights incident at v)  -  1.

    Parameters
    ----------
    A_w : sympy.Matrix
        Symmetric weighted adjacency.

    Returns
    -------
    sympy.Matrix, diagonal.
    """
    M = _require_symmetric(A_w)
    n = M.shape[0]
    Q = sp.zeros(n, n)
    for v in range(n):
        row_sum = sum(M[v, j] for j in range(n))
        Q[v, v] = sp.simplify(row_sum - 1)
    return Q


def connectivity_pattern(A_w: sp.Matrix) -> np.ndarray:
    """Reduce a symbolic weighted adjacency to its 0/1 connectivity pattern.

    An edge is present iff the weight is a nonzero sympy expression.
    Used to compute V, E, r = Betti-1 via the existing unweighted
    helpers (these are combinatorial quantities and do not depend on
    the weights).
    """
    M = _require_symmetric(A_w)
    n = M.shape[0]
    B = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if M[i, j] != 0:
                B[i, j] = 1
                B[j, i] = 1
    return B


def weighted_ihara_zeta_bass(
    A_w: sp.Matrix,
    num_components: Optional[int] = None,
    s: sp.Symbol = s_sym,
) -> sp.Expr:
    """Weighted Ihara-Bass zeta as a symbolic sympy expression.

    Implements the Mizuno-Sato identity

        zeta_G^w(s)^{-1} = (1 - s^2)^{r - c}  det(I - s A_w + s^2 Q_w),

    returning zeta^{-1} as an EXPANDED sympy expression in s and the
    symbolic parameters that appear in A_w (typically alpha^2).

    Parameters
    ----------
    A_w : sympy.Matrix
        Symmetric weighted adjacency with symbolic (or numeric) weights.
    num_components : int, optional
        Override for c. If None, inferred from the 0/1 connectivity pattern.
    s : sympy.Symbol
        Variable of the zeta function.

    Returns
    -------
    sympy.Expr
        zeta_G^w(s)^{-1}, an expanded polynomial in s (coefficients
        lie in Q(alpha^2) or whatever ring the weights live in).
    """
    M = _require_symmetric(A_w)
    V = M.shape[0]
    pattern = connectivity_pattern(M)
    E = int(pattern.sum()) // 2
    c = _count_components(pattern) if num_components is None else int(num_components)
    r = E - V + c  # Betti-1

    Q = weighted_q_matrix(M)
    I = sp.eye(V)
    mat = I - s * M + s * s * Q
    det = sp.expand(mat.det())
    # (1 - s^2)^{r - c}
    exp_pref = r - c
    prefactor = (1 - s * s) ** exp_pref
    zeta_inv = sp.expand(prefactor * det)
    return zeta_inv


def weighted_hashimoto_matrix(A_w: sp.Matrix) -> np.ndarray:
    """Build the (2E) x (2E) weighted Hashimoto edge operator as a
    numeric matrix. ``A_w`` must have numeric (non-symbolic) weights.

    (T_w)[(u->v), (x->y)] = sqrt(w(uv) * w(xy))  if v = x and y != u,
                         = 0                    otherwise.

    The sqrt() is needed so that det(I - s T_w) reproduces the
    Mizuno-Sato weighted Ihara-Bass identity (equivalently, the
    directed-edge operator uses sqrt(w) on each endpoint leg).

    Parameters
    ----------
    A_w : sympy.Matrix
        Weighted adjacency. If the weights are symbolic, substitute
        numeric values BEFORE calling this function (or use
        weighted_ihara_zeta_bass which handles symbolic natively).

    Returns
    -------
    numpy.ndarray of shape (2E, 2E)
    """
    M = _require_symmetric(A_w)
    n = M.shape[0]
    # Numericise; require purely numeric at this point
    A_num = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            val = M[i, j]
            val_num = float(sp.nsimplify(val).evalf()) if val != 0 else 0.0
            A_num[i, j] = val_num
    # Enumerate oriented edges.
    oriented: List[Tuple[int, int]] = []
    for u in range(n):
        for v in range(n):
            if A_num[u, v] != 0.0:
                oriented.append((u, v))
    m2 = len(oriented)
    T = np.zeros((m2, m2), dtype=float)
    for i, (u, v) in enumerate(oriented):
        w_uv = A_num[u, v]
        for j, (x, y) in enumerate(oriented):
            if v == x and y != u:
                w_xy = A_num[x, y]
                T[i, j] = float(np.sqrt(w_uv * w_xy))
    return T


def is_weighted_ramanujan(
    A_w: sp.Matrix,
    tol: float = 1e-9,
) -> Tuple[bool, float, str]:
    """Decide whether the weighted graph satisfies the weighted Kotani-
    Sunada Ramanujan bound.

    For a WEIGHTED graph, the natural generalization (Setyadi-Storm
    2022; Terras 2011 Chap. 18) uses the weighted q_max,

        q_max^w  =  max_v (sum over j~v of w(vj))  -  1
                =  max_v Q_w[v, v],

    and the graph is (weakly) weighted-Ramanujan iff every
    non-trivial eigenvalue mu of T_w satisfies

        |mu| <= sqrt(q_max^w).

    The Perron eigenvalue rho(T_w) = q_max^w for a regular weighted
    graph; for irregular weighted graphs we identify the trivial
    eigenvalues with the spectral radius and its conjugates.

    Parameters
    ----------
    A_w : sympy.Matrix
        Numeric-valued weighted adjacency. Substitute alpha to a
        number (e.g. 1/137.036 or 1) BEFORE calling this function.

    Returns
    -------
    (is_ram, deviation, explanation)
        deviation = max|mu_nontrivial| - sqrt(q_max^w).
        Negative => Ramanujan.
    """
    M = _require_symmetric(A_w)
    V = M.shape[0]
    # Weighted q_max = max row sum of A_w, minus 1.
    row_sums = [float(sum(M[i, j] for j in range(V)).evalf()) for i in range(V)]
    q_max_w = max(row_sums) - 1.0

    T = weighted_hashimoto_matrix(M)
    if T.size == 0:
        return True, 0.0, "trivial: zero-edge graph"
    ev = np.linalg.eigvals(T)
    mags = np.abs(ev)
    rho = float(mags.max())
    non_trivial = mags < rho - tol
    mu_nt_max = float(mags[non_trivial].max()) if non_trivial.any() else 0.0

    bound = float(np.sqrt(max(q_max_w, 0.0)))
    deviation = mu_nt_max - bound
    is_ram = deviation <= tol

    text = (
        f"V={V}, 2E={T.shape[0]}, "
        f"q_max_w={q_max_w:.6f}, sqrt(q_max_w)={bound:.6f}, "
        f"rho(T_w)={rho:.6f}, max|mu_nontrivial|={mu_nt_max:.6f}, "
        f"deviation={deviation:+.6f} "
        f"({'Ramanujan' if is_ram else 'NOT Ramanujan'})"
    )
    return is_ram, float(deviation), text


# ---------------------------------------------------------------------------
# Dirac-S^3-specific edge weighting (alpha^2 spin-orbit)
# ---------------------------------------------------------------------------

def _f_so(lab: DiracLabel) -> sp.Expr:
    """Dimensionless Breit-Pauli SO diagonal (alpha^2 stripped, Z=1).

        f_SO(n, kappa) = -(kappa+1) / [ 4 * n^3 * l * (l+1/2) * (l+1) ]
                       = 0 if l = 0 (Kramers).

    Returns an exact sympy Rational.
    """
    n = lab.n_fock
    kappa = lab.kappa
    l = kappa_to_l(kappa)
    if l == 0:
        return sp.Integer(0)
    # Note: the denominator l*(l+1/2)*(l+1) matches spin_orbit.py.
    numer = -sp.Integer(kappa + 1)
    denom = sp.Integer(4) * sp.Integer(n)**3 * sp.Integer(l) \
            * sp.Rational(2 * l + 1, 2) * sp.Integer(l + 1)
    return sp.Rational(numer, denom) if denom != 0 else sp.Integer(0)


def dirac_s3_edge_weight(
    a: DiracLabel,
    b: DiracLabel,
    alpha: sp.Expr = alpha_sym,
) -> sp.Expr:
    """alpha^2-perturbed SO edge weight on the Dirac-S^3 graph.

        w(a, b) = 1 + alpha^2 * | f_SO(a) + f_SO(b) | / 2

    This is > 0 for all real alpha and = 1 when alpha = 0, so the
    unweighted Ihara zeta is recovered as alpha -> 0.

    Parameters
    ----------
    a, b : DiracLabel
        Endpoints of the edge.
    alpha : sympy Expr, default symbolic alpha_sym.

    Returns
    -------
    sympy Expr.
    """
    fa = _f_so(a)
    fb = _f_so(b)
    mean = (fa + fb) / sp.Integer(2)
    # sp.Abs is idempotent for nonneg, and correctly symbolic for signed.
    return sp.Integer(1) + alpha ** 2 * sp.Abs(mean)


def build_weighted_dirac_adjacency(
    n_max: int,
    adjacency_rule: str,
    alpha: sp.Expr = alpha_sym,
) -> Tuple[sp.Matrix, List[DiracLabel]]:
    """Build the symbolic weighted adjacency for the Dirac-S^3 graph.

    Delegates the 0/1 topology construction to
    ``geovac.ihara_zeta_dirac.build_dirac_s3_graph`` and then lifts
    each edge to its alpha^2-SO weight via ``dirac_s3_edge_weight``.

    Parameters
    ----------
    n_max : int
    adjacency_rule : {"A", "B"}
    alpha : sympy.Expr, symbolic alpha by default.

    Returns
    -------
    (A_w, labels) with A_w a sympy Matrix.
    """
    # Import here to avoid circular dependency at module load.
    from geovac.ihara_zeta_dirac import build_dirac_s3_graph
    A01, labels, _, _ = build_dirac_s3_graph(n_max, adjacency_rule)
    V = A01.shape[0]
    A_w = sp.zeros(V, V)
    for i in range(V):
        for j in range(i + 1, V):
            if A01[i, j]:
                w = dirac_s3_edge_weight(labels[i], labels[j], alpha=alpha)
                A_w[i, j] = w
                A_w[j, i] = w
    return A_w, labels


# ---------------------------------------------------------------------------
# Convenience: pair-correlation diagnostic for Hashimoto spectrum
# ---------------------------------------------------------------------------

def hashimoto_pair_correlation_cv(A_w: sp.Matrix) -> Dict[str, float]:
    """Coefficient-of-variation diagnostic for nearest-neighbor
    spacings of the Hashimoto spectrum. Used as a crude surrogate
    for "how close to GUE" the spectrum sits.

    For a Poisson (uncorrelated) spectrum, CV -> 1.
    For a GUE (level-repelled) spectrum, CV -> 0.522.
    For a GOE spectrum, CV -> 0.523.

    This is NOT a rigorous RMT test — it is a first-moment diagnostic
    that distinguishes Poisson-like from Wigner-like.

    Returns
    -------
    dict with keys 'cv_spacing', 'n_spacings', 'mean_spacing'.
    """
    T = weighted_hashimoto_matrix(A_w)
    if T.shape[0] < 3:
        return {"cv_spacing": float("nan"), "n_spacings": 0,
                "mean_spacing": 0.0}
    ev = np.linalg.eigvals(T)
    mags = np.sort(np.abs(ev))
    # Unfold crudely by normalizing to mean spacing 1.
    spacings = np.diff(mags)
    # Drop zero spacings (degeneracies):
    nz = spacings[spacings > 1e-9]
    if len(nz) < 2:
        return {"cv_spacing": float("nan"), "n_spacings": int(len(nz)),
                "mean_spacing": 0.0}
    mean = float(nz.mean())
    std = float(nz.std())
    cv = std / mean if mean > 0 else float("nan")
    return {
        "cv_spacing": cv,
        "n_spacings": int(len(nz)),
        "mean_spacing": mean,
    }


__all__ = [
    "alpha_sym",
    "weighted_adjacency",
    "weighted_q_matrix",
    "weighted_ihara_zeta_bass",
    "weighted_hashimoto_matrix",
    "is_weighted_ramanujan",
    "connectivity_pattern",
    "dirac_s3_edge_weight",
    "build_weighted_dirac_adjacency",
    "hashimoto_pair_correlation_cv",
]
