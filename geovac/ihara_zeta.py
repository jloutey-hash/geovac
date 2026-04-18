"""
Ihara zeta function of a finite graph.

This module implements the Ihara zeta function of a finite, undirected,
loopless graph and the associated graph-theoretic analogs of the
Riemann Hypothesis.

The three canonical objects
---------------------------

1. **Euler product over primitive closed walks.**  A closed walk is
   primitive if it has no backtracking (no immediate reversal of an
   edge) and no tail (it is a genuine closed path, not a path with a
   repeated prefix/suffix).  Two closed walks that differ by a cyclic
   rotation represent the same equivalence class.  The Ihara zeta is
   defined by

       zeta_G(s) = prod over [C] primitive, (1 - s^{len(C)})^{-1}

   where the product runs over equivalence classes of primitive closed
   walks and len(C) is the edge length.  The variable s plays the role
   of e^{-s} in the Selberg-zeta analogy.

2. **Ihara-Bass determinantal formula (Bass 1992).**  For a finite
   connected graph with V vertices and E edges,

       zeta_G(s)^{-1} = (1 - s^2)^{r-1} det(I - s A + s^2 Q)

   where A is the V x V adjacency matrix, Q = diag(d_v - 1) is the
   "excess-degree" diagonal, and r = E - V + 1 is the first Betti
   number (number of independent cycles).  For a graph with c
   connected components, r = E - V + c.

3. **Ihara-Hashimoto edge operator.**  Let T be the 2E x 2E
   "non-backtracking" matrix, indexed by oriented edges:
   T_{(u->v),(x->y)} = 1 iff v = x and y != u, and 0 otherwise.
   Then

       zeta_G(s)^{-1} = det(I - s T)

   so the zeros of zeta_G in s are the reciprocals of the nonzero
   eigenvalues of T.  This formulation gives the correct notion of
   the "graph Riemann Hypothesis" for irregular graphs.

Graph Riemann Hypothesis (Ramanujan bound)
------------------------------------------

For a (q+1)-regular graph, the non-trivial zeros of zeta_G(s) lie on
|s| = 1/sqrt(q) iff the graph is Ramanujan, i.e. all non-trivial
eigenvalues of A satisfy |lambda| <= 2 sqrt(q).  The analog for
irregular graphs uses the Hashimoto edge operator T: the graph is
"weakly Ramanujan" iff every eigenvalue mu of T satisfies either
|mu| = 1 (trivial), |mu| = sqrt(q_max) (sharp), or |mu| <= sqrt(q_max)
(bulk), where q_max = max_v (d_v - 1) is the largest excess degree.
The "graph RH" for irregular graphs is understood in the sense of
Kotani-Sunada (2000): the non-trivial spectrum of T (excluding the
Perron eigenvalue and its conjugates) should satisfy
sqrt(q_min) <= |mu| <= sqrt(q_max) with the "critical line"
|mu| = sqrt(q_avg).  A cleaner RH statement for non-regular graphs
uses the bipartite double cover and measures deviation from the
Ramanujan bound sqrt(q_max).

In this module the check `is_ramanujan` uses the Hashimoto-spectrum
definition.  It reports the deviation

    deviation = max_{mu nontrivial} |mu| - sqrt(q_max)

where q_max = max_v(d_v - 1).  Negative deviation means the graph
respects the Ramanujan bound.

References
----------
- Y. Ihara, "On discrete subgroups of the two by two projective linear
  group over p-adic fields," J. Math. Soc. Japan 18 (1966), 219-235.
- H. Bass, "The Ihara-Selberg zeta function of a tree lattice," Int.
  J. Math. 3 (1992), 717-797.
- M. Kotani and T. Sunada, "Zeta functions of finite graphs,"
  J. Math. Sci. Univ. Tokyo 7 (2000), 7-25.
- A. Terras, "Zeta Functions of Graphs: A Stroll through the Garden,"
  Cambridge Univ. Press, 2011.  (Primary reference for this module.)
- H. Stark and A. Terras, "Zeta functions of finite graphs and
  coverings," Adv. Math. 121 (1996), 124-165.

Author: GeoVac Track RH-A (Sprint: "hey-buddy-we-need-crystalline-sprout")
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import sympy as sp


# ----------------------------------------------------------------------------
# Helper: adjacency normalization
# ----------------------------------------------------------------------------

def _as_integer_adjacency(adjacency) -> np.ndarray:
    """
    Coerce an adjacency-like object to a square, symmetric, zero-diagonal
    numpy array of nonnegative integers.  If the input is weighted (e.g.
    rational squared dipole elements from the Bargmann graph), we reduce
    it to its unweighted 0/1 connectivity pattern: an edge is present iff
    the weight is nonzero.  The Ihara zeta is a *combinatorial* invariant,
    it does not depend on edge weights.
    """
    if hasattr(adjacency, "toarray"):
        A = np.asarray(adjacency.toarray())
    else:
        A = np.asarray(adjacency)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("adjacency must be square, 2-D")
    if not np.allclose(A, A.T):
        raise ValueError("adjacency must be symmetric")
    # Coerce to 0/1 connectivity pattern (Ihara zeta is combinatorial).
    B = (np.abs(A) > 1e-12).astype(int)
    np.fill_diagonal(B, 0)  # strip self-loops (Ihara zeta is for simple graphs)
    return B


def _degree_sequence(A: np.ndarray) -> np.ndarray:
    """Return the degree vector of the adjacency A."""
    return A.sum(axis=1).astype(int)


def _count_components(A: np.ndarray) -> int:
    """
    Count connected components of the graph whose adjacency is A via BFS.
    """
    n = A.shape[0]
    seen = [False] * n
    comps = 0
    for s in range(n):
        if seen[s]:
            continue
        comps += 1
        stack = [s]
        seen[s] = True
        while stack:
            u = stack.pop()
            neigh = np.nonzero(A[u])[0]
            for v in neigh:
                if not seen[v]:
                    seen[v] = True
                    stack.append(int(v))
    return comps


# ----------------------------------------------------------------------------
# Ihara-Bass determinantal zeta
# ----------------------------------------------------------------------------

def ihara_zeta_bass(
    adjacency,
    num_components: Optional[int] = None,
    symbolic: bool = True,
):
    """
    Compute the Ihara zeta function via the Ihara-Bass determinantal
    formula:

        zeta_G(s)^{-1} = (1 - s^2)^{r-1} det(I - s A + s^2 Q)

    where Q = diag(d_v - 1), r = E - V + c.  This is the most efficient
    closed-form route to zeta_G on a general graph.

    Parameters
    ----------
    adjacency
        V x V adjacency (scipy sparse or numpy array).  Weighted inputs
        are reduced to their 0/1 connectivity pattern — the Ihara zeta
        is a combinatorial invariant.
    num_components
        Optional override for the number of connected components c.
        If None, c is computed by BFS.
    symbolic
        If True (default), return a sympy rational function zeta_G(s)
        (as an unevaluated expression).  If False, return a tuple
        (numerator_coeffs, denominator_coeffs) of numpy polynomial
        coefficients, highest-degree-first.

    Returns
    -------
    sympy.Expr  (if symbolic=True)
        The rational function zeta_G(s) = 1 / [(1-s^2)^{r-1} det(...)].
    (num, den)  (if symbolic=False)
        Two numpy arrays of polynomial coefficients.

    References
    ----------
    Bass (1992), Terras (2011) Theorem 2.5.
    """
    A = _as_integer_adjacency(adjacency)
    V = A.shape[0]
    E = int(A.sum()) // 2
    c = _count_components(A) if num_components is None else int(num_components)
    r = E - V + c  # total Betti-1 (sum over components)

    d = _degree_sequence(A)
    Q = np.diag(d - 1)

    s = sp.symbols("s")

    # Sympy determinant of (I - s A + s^2 Q).  The matrix is V x V with
    # integer entries in s, so sympy handles it in pure rational arithmetic.
    I = sp.eye(V)
    # Build the polynomial-in-s matrix.
    M = I - s * sp.Matrix(A.tolist()) + s * s * sp.Matrix(Q.tolist())
    det_poly = sp.expand(M.det())

    # zeta^{-1} = (1 - s^2)^{r - c} * det(...)
    # Terras's connected formula has (1 - s^2)^{r - 1}, which assumes
    # c = 1.  On c > 1 components, ζ_G factors as the product of
    # ζ_{G_i} over components, giving a combined prefactor
    # ∏_i (1 - s^2)^{r_i - 1} = (1 - s^2)^{(r - c)} because
    # sum_i r_i = r and sum_i 1 = c.  This is the correct convention
    # for graphs that may be disconnected.
    prefactor = (1 - s * s) ** (r - c)
    zeta_inv = sp.expand(prefactor * det_poly)

    if symbolic:
        return 1 / zeta_inv

    # Convert to polynomial coefficient arrays.  Note that zeta^{-1}
    # may in general be a Laurent polynomial when r < 1 (because of
    # the (1 - s^2)^{r-1} factor with negative exponent).  In that
    # case we return numerator and denominator separately.
    num, den = sp.fraction(sp.together(zeta_inv))
    num_poly = sp.Poly(num, s) if num != 0 else sp.Poly(1, s)
    den_poly = sp.Poly(den, s) if den != 0 else sp.Poly(1, s)
    num_coeffs = np.array([sp.nsimplify(c).evalf() for c in num_poly.all_coeffs()],
                          dtype=complex)
    den_coeffs = np.array([sp.nsimplify(c).evalf() for c in den_poly.all_coeffs()],
                          dtype=complex)
    # Since zeta = 1 / zeta_inv, swap:
    return (den_coeffs, num_coeffs)


# ----------------------------------------------------------------------------
# Ihara-Hashimoto edge operator
# ----------------------------------------------------------------------------

def hashimoto_matrix(adjacency) -> np.ndarray:
    """
    Build the Ihara-Hashimoto non-backtracking edge matrix T of size
    (2E) x (2E).

    T is indexed by oriented edges e = (u -> v).  We set
    T[(u->v), (x->y)] = 1 iff v = x AND y != u (no backtracking),
    else 0.

    Parameters
    ----------
    adjacency
        V x V adjacency (scipy sparse or numpy array).  Reduced to
        connectivity pattern.

    Returns
    -------
    T : numpy.ndarray of shape (2E, 2E)
        The Hashimoto edge matrix.  Rows/columns index a canonical
        orientation of the edges: for each undirected edge {u, v},
        both (u -> v) and (v -> u) appear.

    References
    ----------
    Hashimoto (1989); Bass (1992); Terras (2011) §4.
    """
    A = _as_integer_adjacency(adjacency)
    V = A.shape[0]
    # Enumerate oriented edges.
    oriented: List[Tuple[int, int]] = []
    edge_index: Dict[Tuple[int, int], int] = {}
    for u in range(V):
        for v in np.nonzero(A[u])[0]:
            oriented.append((u, int(v)))
            edge_index[(u, int(v))] = len(oriented) - 1

    m2 = len(oriented)
    T = np.zeros((m2, m2), dtype=int)
    for i, (u, v) in enumerate(oriented):
        for j, (x, y) in enumerate(oriented):
            if v == x and y != u:
                T[i, j] = 1
    return T


def is_ramanujan(
    adjacency,
    tol: float = 1e-9,
) -> Tuple[bool, float, str]:
    """
    Decide whether the graph satisfies the (weak) Ramanujan bound via
    the spectrum of the Hashimoto edge operator.

    For a (q+1)-regular graph, the graph is Ramanujan iff
    every non-trivial eigenvalue mu of T satisfies |mu| = sqrt(q) or
    |mu| <= sqrt(q) (Bass; Lubotzky-Phillips-Sarnak).  For irregular
    graphs, Kotani-Sunada (2000) gives the appropriate generalization
    using q_max = max_v(d_v - 1).

    This function computes:
        q_max          = max_v (d_v - 1)
        q_min          = min_v (d_v - 1)
        mu_max_nt      = max |mu| over non-trivial T-eigenvalues
                         (where "trivial" means real |mu| = 1 associated
                         with the Perron eigenvalue and degree-1 legs).

    A graph is called **(q_max-)Ramanujan** if mu_max_nt <= sqrt(q_max).

    Parameters
    ----------
    adjacency
        V x V adjacency.
    tol
        Numerical tolerance for "trivial" vs "non-trivial" eigenvalue
        separation.

    Returns
    -------
    (is_ram, deviation, explanation)
        - is_ram: True iff mu_max_nt <= sqrt(q_max) + tol.
        - deviation = mu_max_nt - sqrt(q_max).  Negative iff Ramanujan.
        - explanation: short textual description of the verdict.
    """
    A = _as_integer_adjacency(adjacency)
    d = _degree_sequence(A)
    q_max = int(d.max()) - 1
    q_min = int(d.min()) - 1
    regular = (d.min() == d.max())

    T = hashimoto_matrix(A)
    ev = np.linalg.eigvals(T.astype(float))
    mags = np.abs(ev)

    # The "trivial" eigenvalues of T lie on |mu| = 1 (from dead-ends /
    # leaves) and at mu = +/-Perron eigenvalue.  For a connected graph
    # of Perron spectral radius rho(T), the two sharpest trivial values
    # are +rho and -rho (for bipartite graphs, both; for non-bipartite,
    # only +rho is certain).  We identify "trivial" as: eigenvalues
    # whose |mu| matches max |mu| within tol.  All other eigenvalues
    # are "non-trivial".
    rho_T = float(mags.max())
    non_trivial_mask = mags < rho_T - tol
    if non_trivial_mask.any():
        mu_nt_max = float(mags[non_trivial_mask].max())
    else:
        mu_nt_max = 0.0

    bound = float(np.sqrt(max(q_max, 0)))
    deviation = mu_nt_max - bound
    is_ram = deviation <= tol

    text = (
        f"V={A.shape[0]}, E={A.sum()//2}, regular={regular}, "
        f"q_max={q_max}, q_min={q_min}, "
        f"rho(T)={rho_T:.6f} (Perron), "
        f"max|mu_nontrivial|={mu_nt_max:.6f}, "
        f"sqrt(q_max)={bound:.6f}, "
        f"deviation={deviation:+.6f} "
        f"({'Ramanujan' if is_ram else 'NOT Ramanujan'})"
    )
    return is_ram, float(deviation), text


# ----------------------------------------------------------------------------
# Euler-product cross-check (truncated to walks of length <= L_max)
# ----------------------------------------------------------------------------

def _count_closed_nonbacktracking_walks(A: np.ndarray, max_length: int) -> Dict[int, int]:
    """
    Count, for each 1 <= L <= max_length, the number of CLOSED non-
    backtracking walks of edge length L.

    A closed non-backtracking walk of length L is a sequence of oriented
    edges (e_1, ..., e_L) with head(e_i) = tail(e_{i+1}), e_{i+1} not
    the reverse of e_i (no backtracking), and head(e_L) = tail(e_1).
    Such walks are counted with orientation (starting edge is part of
    the data).  They include all cyclic rotations of a given primitive
    closed walk.

    Uses T^L and tr(T^L) = number of closed NB walks of length L.
    """
    T = hashimoto_matrix(A).astype(float)
    T_power = np.eye(T.shape[0])
    counts = {}
    for L in range(1, max_length + 1):
        T_power = T_power @ T
        counts[L] = int(round(T_power.trace()))
    return counts


def _mobius_to_primitive(tr_counts: Dict[int, int]) -> Dict[int, int]:
    """
    From the number of closed non-backtracking walks N_L = tr(T^L)
    (counted with orientation and starting point), extract the number
    N_prim(L) of primitive closed walks of length L (up to cyclic
    rotation).

    A classical Mobius-inversion identity (Terras 2011 §4) gives:

        tr(T^L) = sum_{d | L} d * N_prim(d)

    (each primitive walk of length d accounts for d rotations in length
    L if d divides L, and 0 otherwise).  Therefore:

        L * N_prim(L) = sum_{d | L} mu(L/d) * tr(T^d)

    Returns {L: N_prim(L)}.
    """
    def mobius(n: int) -> int:
        if n == 1:
            return 1
        # factorize
        factors = {}
        m = n
        p = 2
        while p * p <= m:
            while m % p == 0:
                factors[p] = factors.get(p, 0) + 1
                m //= p
            p += 1
        if m > 1:
            factors[m] = factors.get(m, 0) + 1
        for exp in factors.values():
            if exp >= 2:
                return 0
        return (-1) ** len(factors)

    prim = {}
    for L in tr_counts:
        s = 0
        for d in range(1, L + 1):
            if L % d == 0:
                s += mobius(L // d) * tr_counts[d]
        prim[L] = s // L if s % L == 0 else None
    return prim


def ihara_zeta_euler(
    adjacency,
    max_length: int = 6,
    symbolic: bool = True,
):
    """
    Compute the truncated Euler product for the Ihara zeta:

        zeta_G(s) ~ prod_{L=1}^{max_length} (1 - s^L)^{-N_prim(L)}

    valid to order s^{max_length}.  This is used as a consistency
    check against `ihara_zeta_bass`.

    Parameters
    ----------
    adjacency
        V x V adjacency.
    max_length
        Maximum primitive-walk length to retain.
    symbolic
        If True, return a sympy rational function; else return the
        dictionary {L: N_prim(L)}.

    Returns
    -------
    sympy.Expr or dict
        Truncated zeta as a product of factors (1 - s^L)^{-N_prim(L)},
        or the raw primitive-walk multiplicities.
    """
    A = _as_integer_adjacency(adjacency)
    tr = _count_closed_nonbacktracking_walks(A, max_length)
    prim = _mobius_to_primitive(tr)
    if not symbolic:
        return prim
    s = sp.symbols("s")
    expr = sp.Integer(1)
    for L, n in prim.items():
        if n is None or n == 0:
            continue
        expr = expr * (1 - s ** L) ** (-n)
    return sp.simplify(expr)


# ----------------------------------------------------------------------------
# Functional equation report
# ----------------------------------------------------------------------------

def functional_equation_report(
    adjacency,
    num_components: Optional[int] = None,
) -> Dict:
    """
    Identify the functional equation satisfied by zeta_G.

    For a (q+1)-regular graph, the Ihara zeta satisfies a clean
    Riemann-type functional equation in the variable u = q s^2:

        Xi_G(s) := (1 - s^2)^{r-1} (1 - q^2 s^2)^{r-1} zeta_G(s)^{-1}
                  = Xi_G(1 / (q s))         (Stark-Terras 1996, Thm 3).

    Equivalently, the zeros of zeta_G are symmetric under s -> 1/(q s)
    up to the trivial poles at s = +/- 1 and s = +/- 1/q.

    For an irregular graph, a cleanly stated functional equation uses
    the Hashimoto edge operator T and its spectrum; the completed
    Xi-function has a less compact closed form but the *zero set* of
    zeta_G is still symmetric under s -> 1/s_conj where s_conj is the
    reciprocal of any T-eigenvalue conjugate.  We report:

      - regularity: bool, with the effective q if regular.
      - if regular: the functional-equation exponent and critical
        radius 1/sqrt(q).
      - if irregular: the spectrum of T and the observed symmetry
        (closure of the spectrum under mu -> q_max/mu, which for
        Ramanujan graphs gives the standard critical circle
        |s| = 1/sqrt(q_max)).

    Returns
    -------
    dict with keys:
        'regular', 'q_effective', 'critical_radius',
        'functional_equation', 'zero_symmetry_observed'.
    """
    A = _as_integer_adjacency(adjacency)
    d = _degree_sequence(A)
    regular = bool(d.min() == d.max())
    V = A.shape[0]
    E = int(A.sum()) // 2
    c = _count_components(A) if num_components is None else int(num_components)
    r = E - V + c

    info: Dict = {
        "V": V,
        "E": E,
        "c": c,
        "r_betti1": r,
        "regular": regular,
    }

    if regular:
        q = int(d[0]) - 1
        info.update({
            "q_effective": q,
            "critical_radius": float(1.0 / np.sqrt(q)) if q > 0 else None,
            "functional_equation": (
                f"Xi_G(s) := (1 - s^2)^{r-1} (1 - ({q}*s)^2)^{r-1} "
                f"zeta_G(s)^(-1) satisfies Xi_G(s) = Xi_G(1/({q}*s))"
            ),
            "zero_symmetry": f"s -> 1/({q}*s)",
        })
    else:
        q_max = int(d.max()) - 1
        q_min = int(d.min()) - 1
        T = hashimoto_matrix(A)
        ev = np.linalg.eigvals(T.astype(float))
        info.update({
            "q_max": q_max,
            "q_min": q_min,
            "q_avg_nb": float((d - 1).mean()),
            "critical_radius_qmax": float(1.0 / np.sqrt(q_max)) if q_max > 0 else None,
            "critical_radius_qmin": float(1.0 / np.sqrt(max(q_min, 1))) if q_min >= 1 else None,
            "functional_equation": (
                "Irregular graph: the symmetric completed xi has no "
                "uniform closed form in s alone (Terras 2011 §3, §4). "
                "The Hashimoto T-spectrum {mu} determines the zeros of "
                "zeta_G via s = 1/mu.  For zeros on a single critical "
                "circle |s| = 1/sqrt(q), the graph must be regular or "
                "satisfy the Kotani-Sunada uniform bound."
            ),
            "zero_symmetry": "s_j = 1 / mu_j (mu_j eigenvalues of Hashimoto T)",
            "hashimoto_spectral_radius": float(np.max(np.abs(ev))),
        })
    return info


# ----------------------------------------------------------------------------
# Zero extraction (for reports)
# ----------------------------------------------------------------------------

def zeta_zeros_from_hashimoto(
    adjacency,
    return_eigs: bool = False,
):
    """
    Return the zeros of zeta_G(s) as the reciprocals of the nonzero
    eigenvalues of the Hashimoto edge operator T.  (Bass 1992:
    zeta^{-1} = det(I - sT).)

    Parameters
    ----------
    adjacency
        V x V adjacency.
    return_eigs
        If True, also return the T-eigenvalues themselves.

    Returns
    -------
    zeros : numpy array of complex numbers
        zeros[k] = 1 / mu_k for each nonzero eigenvalue mu_k of T.
    (optional) eigs : numpy array
        The Hashimoto eigenvalues.
    """
    A = _as_integer_adjacency(adjacency)
    T = hashimoto_matrix(A)
    ev = np.linalg.eigvals(T.astype(float))
    nonzero = np.abs(ev) > 1e-12
    zeros = 1.0 / ev[nonzero]
    if return_eigs:
        return zeros, ev
    return zeros


__all__ = [
    "ihara_zeta_bass",
    "ihara_zeta_euler",
    "hashimoto_matrix",
    "is_ramanujan",
    "functional_equation_report",
    "zeta_zeros_from_hashimoto",
]
