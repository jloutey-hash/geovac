"""
Electron propagator for graph-native QED on the GeoVac Dirac graph.
====================================================================

The GeoVac Dirac graph has nodes labeled |n, κ, m_j⟩ (DiracLabel),
diagonal weights λ_n = χ·(n + 3/2) from the Camporesi-Higuchi spectrum,
and off-diagonal adjacency from E1 dipole selection rules (Rule B from
the Ihara zeta Dirac analysis).

The Dirac graph operator is:

    D_graph[i, j] = λ_i  δ_{ij}  +  t · A[i,j]

where
  - λ_i = χ_i · (n_i + 3/2)  is the diagonal Camporesi-Higuchi eigenvalue
    (exact rational, with sign from chirality)
  - A[i,j] ∈ {0, 1}  is the E1 dipole adjacency (Rule B)
  - t  is the off-diagonal hopping strength

The electron propagator is:

    G_e = (D_graph)^{-1}       (massless)
    G_e(m) = (D_graph - m·I)^{-1}   (massive, mass parameter m)

This module computes D_graph and G_e in exact sympy rational arithmetic
at n_max=2 (10 states in atomic mode), and in numpy at n_max=3 (28 states).

Two exact inversion methods are available:
  - Direct: M.inv() via sympy (Gauss elimination / adjugate).
  - Neumann series: G_e = [sum_k (-t)^k (Lambda^{-1} A)^k] Lambda^{-1},
    which terminates at order N-1 by Cayley-Hamilton.  Algebraically
    cleaner at t != 0 because expression trees stay factored.
    See neumann_electron_propagator().

Off-diagonal coupling t:
------------------------
Three natural choices are explored:
  t=0   : diagonal-only; D_graph = diag(λ); G_e = diag(1/λ). π-free, rational.
  t=1   : unit-weight adjacency (graph topology only)
  t=κ_graph = -1/16 : the scalar graph topological coupling constant

The key finding: the scalar κ = -1/16 arises from the Fock projection
(chordal-distance weighting on S³). The Dirac graph adjacency (Rule B,
E1 dipole) is a *different* graph on the same node set — it encodes
transition selection rules, not nearest-neighbor Fock coupling. Thus
t = κ_graph is not the natural Dirac hopping. Instead:
  - t=0 gives the free Dirac propagator (exact rational, clean baseline)
  - t=1 gives the adjacency-perturbed propagator (checks the topology)
  - A physically motivated coupling would come from computing the
    Szmytkowski σ·r̂ reduced matrix element between adjacent states.

The diagonal propagator (t=0) is the most natural graph-QED object:
G_e[n,κ,m_j] = 1 / (χ·(n+3/2)), exact rational in ℚ.

π-free certificate:
-------------------
For t=0: all entries of D_graph and G_e are in ℚ (exact rationals).
For t=1: D_graph ∈ ℚ (adjacency is 0/1); G_e is a matrix of rationals
provided det(D_graph) in ℚ \\ {0}. Since all entries are in ℚ, Cramer's rule
gives G_e entries in ℚ. Verified by explicit sympy computation.

Transcendental taxonomy (Paper 18):
------------------------------------
- D_graph diagonal: ℚ (half-integer eigenvalues, exact Rational)
- D_graph off-diagonal: ℚ (0 or t, adjacency integer weights)
- G_e (t=0): ℚ (diagonal, 2/(2n+3) per state)
- G_e (t≠0): ℚ (rational matrix inverse of a rational matrix)
All quantities are in the intrinsic tier — no π, no ζ, no Hurwitz.

References
----------
- Camporesi & Higuchi (1996), J. Geom. Phys. 20 (1996) 1-18.
- GeoVac dirac_lattice.py (DiracLattice infrastructure)
- GeoVac dirac_matrix_elements.py (DiracLabel, iter_dirac_labels)
- GeoVac dirac_s3.py (dirac_eigenvalue_abs, Camporesi-Higuchi spectrum)
- CLAUDE.md §2 Dirac-on-S³ infrastructure (Track D1, Tier 1)
"""

from __future__ import annotations

from fractions import Fraction
from typing import Dict, List, Optional, Tuple

import numpy as np
import sympy as sp
from sympy import Matrix, Rational, zeros as sp_zeros

from geovac.dirac_lattice import DiracLattice
from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l

__all__ = [
    "DiracGraphOperator",
    "build_dirac_graph_operator",
    "electron_propagator",
    "neumann_electron_propagator",
    "massive_propagator",
    "propagator_spectrum",
    "propagator_pi_free_certificate",
    "propagator_summary",
]

# The scalar-graph topological coupling constant (Paper 0 / Paper 7)
KAPPA_SCALAR = Rational(-1, 16)


# ---------------------------------------------------------------------------
# Core class
# ---------------------------------------------------------------------------

class DiracGraphOperator:
    """
    The Dirac graph operator D_graph on the finite GeoVac Dirac lattice.

    D_graph = Λ + t·A

    where Λ = diag(λ_0, ..., λ_{N-1}) contains the Camporesi-Higuchi
    eigenvalues (with chirality sign) and A is the E1 dipole adjacency.

    Parameters
    ----------
    n_max : int
        Maximum Fock principal quantum number (≥ 1).
    t : sympy.Rational or int
        Off-diagonal hopping strength. Default 0 (free diagonal propagator).
    mode : str
        'atomic' (l < n_fock, hydrogen convention) or 's3' (full S³ spinor
        harmonics including l = n_fock boundary states).
    """

    def __init__(self, n_max: int, t=Rational(0), mode: str = 'atomic'):
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        self.n_max = n_max
        self.t = Rational(t) if not isinstance(t, sp.Expr) else t
        self.mode = mode

        # Build lattice for adjacency and eigenvalues
        self._lattice = DiracLattice(n_max=n_max, mode=mode)
        self.labels: List[DiracLabel] = self._lattice.labels
        self.N = len(self.labels)

        # Camporesi-Higuchi eigenvalues as exact rationals
        # λ_i = χ_i · (n_fock_i + 1/2)  [Fock convention: n+1/2 = (n_CH)+3/2]
        # In Fock convention: |λ| = n_fock + 1/2 = n_CH + 3/2
        self._eigenvalues: List[Rational] = []
        for i, lab in enumerate(self.labels):
            chi = int(self._lattice.chirality[i])
            # |λ| in Fock convention: n_fock + 1/2
            lam_abs = Rational(2 * lab.n_fock + 1, 2)
            self._eigenvalues.append(Rational(chi) * lam_abs)

        # Build sympy matrix (exact rational)
        self._matrix_sym: Optional[Matrix] = None
        self._adj_csr = self._lattice.adjacency  # scipy sparse

        # numpy float matrix for large n_max
        self._matrix_np: Optional[np.ndarray] = None

    @property
    def eigenvalues(self) -> List[Rational]:
        """Diagonal Camporesi-Higuchi eigenvalues (with chirality sign)."""
        return self._eigenvalues

    def matrix_sympy(self) -> Matrix:
        """Build D_graph as an exact sympy Matrix.

        D_graph[i,i] = λ_i  (exact Rational)
        D_graph[i,j] = t    if A[i,j] = 1
        D_graph[i,j] = 0    otherwise

        Returns a sympy Matrix. Suitable for n_max ≤ 2 (10×10 or smaller).
        """
        if self._matrix_sym is not None:
            return self._matrix_sym

        N = self.N
        M = sp_zeros(N, N)

        # Diagonal: Camporesi-Higuchi eigenvalues
        for i in range(N):
            M[i, i] = self._eigenvalues[i]

        # Off-diagonal: E1 dipole adjacency × t
        if self.t != Rational(0):
            cx = self._adj_csr.tocoo()
            for r, c in zip(cx.row, cx.col):
                if r != c:
                    M[r, c] = self.t

        self._matrix_sym = M
        return M

    def matrix_numpy(self) -> np.ndarray:
        """Build D_graph as a float64 numpy array.

        Suitable for all n_max. Off-diagonal entries are float(t).
        """
        if self._matrix_np is not None:
            return self._matrix_np

        N = self.N
        M = np.zeros((N, N), dtype=np.float64)

        # Diagonal
        for i in range(N):
            M[i, i] = float(self._eigenvalues[i])

        # Off-diagonal
        t_float = float(self.t)
        if t_float != 0.0:
            cx = self._adj_csr.tocoo()
            for r, c in zip(cx.row, cx.col):
                if r != c:
                    M[r, c] = t_float

        self._matrix_np = M
        return M

    def spectrum_numpy(self) -> np.ndarray:
        """Eigenvalues of D_graph (numpy, sorted)."""
        M = self.matrix_numpy()
        eigs = np.linalg.eigvalsh(M)
        return np.sort(eigs)

    def det_sympy(self) -> sp.Expr:
        """Exact determinant of D_graph (sympy)."""
        return self.matrix_sympy().det()

    def is_invertible(self) -> bool:
        """True if D_graph is invertible (det ≠ 0, checked numerically)."""
        M = self.matrix_numpy()
        return abs(np.linalg.det(M)) > 1e-12

    def condition_number(self) -> float:
        """Condition number of D_graph (numpy)."""
        M = self.matrix_numpy()
        return np.linalg.cond(M)

    def spectral_gap(self) -> float:
        """Smallest |eigenvalue| of D_graph."""
        eigs = self.spectrum_numpy()
        return float(np.min(np.abs(eigs)))

    def __repr__(self) -> str:
        return (
            f"DiracGraphOperator(n_max={self.n_max}, t={self.t}, "
            f"mode={self.mode!r}, N={self.N})"
        )


# ---------------------------------------------------------------------------
# Propagator computation
# ---------------------------------------------------------------------------

def build_dirac_graph_operator(
    n_max: int,
    t=Rational(0),
    mode: str = 'atomic'
) -> DiracGraphOperator:
    """Construct the Dirac graph operator for the given parameters.

    Parameters
    ----------
    n_max : int
        Maximum Fock quantum number.
    t : rational or int
        Off-diagonal hopping (default 0 = free diagonal propagator).
    mode : str
        'atomic' or 's3'.

    Returns
    -------
    DiracGraphOperator
    """
    return DiracGraphOperator(n_max=n_max, t=t, mode=mode)


def electron_propagator(
    op: DiracGraphOperator,
    exact: bool = True,
    method: str = 'auto',
) -> Tuple[object, bool]:
    """Compute the massless electron propagator G_e = D_graph^{-1}.

    Parameters
    ----------
    op : DiracGraphOperator
        The Dirac graph operator.
    exact : bool
        If True, use exact sympy Matrix inverse (suitable for N <= ~30).
        If False, use numpy (fast, float64).
    method : str
        'auto' : use Neumann series when exact=True and t != 0;
                 use direct M.inv() when exact=True and t == 0;
                 use numpy when exact=False.
        'direct' : always use M.inv() (original behavior).
        'neumann' : always use the Neumann/Born series (exact sympy only).

    Returns
    -------
    (propagator, is_rational)
        propagator: sympy.Matrix or np.ndarray depending on `exact`
        is_rational: True if all entries are rational (pi-free certificate)
    """
    if method not in ('auto', 'direct', 'neumann'):
        raise ValueError(f"method must be 'auto', 'direct', or 'neumann', got {method!r}")

    # Decide which path to use
    use_neumann = False
    if method == 'neumann':
        if not exact:
            raise ValueError("method='neumann' requires exact=True")
        use_neumann = True
    elif method == 'auto':
        use_neumann = exact and (op.t != Rational(0))

    if use_neumann:
        G, is_rational, _info = neumann_electron_propagator(op)
        return G, is_rational

    if exact:
        M = op.matrix_sympy()
        G = M.inv()
        is_rational = _check_rational_matrix(G)
        return G, is_rational
    else:
        M = op.matrix_numpy()
        G = np.linalg.inv(M)
        return G, None  # Cannot certify pi-free from float


def neumann_electron_propagator(
    op: DiracGraphOperator,
    max_order: Optional[int] = None,
) -> Tuple[Matrix, bool, Dict]:
    """Compute the electron propagator via a Cayley-Hamilton polynomial.

    Uses the factorization:

        G_e = (Lambda + t*A)^{-1}
            = (I + t*Lambda^{-1}*A)^{-1} * Lambda^{-1}

    Let M = I + t*B where B = Lambda^{-1}*A.  The characteristic polynomial
    of M is p(x) = det(x*I - M) = x^N + c_{N-1}*x^{N-1} + ... + c_0.

    By Cayley-Hamilton, p(M) = 0, so:

        M^N = -(c_{N-1}*M^{N-1} + ... + c_1*M + c_0*I)

    Since det(M) = c_0 * (-1)^N, and M^{-1} = -(1/c_0)*(M^{N-1} + c_{N-1}*M^{N-2}
    + ... + c_1*I), this gives M^{-1} as an EXACT polynomial of degree N-1 in M
    (equivalently, in t*B).

    When max_order < N-1, this falls back to the truncated Neumann (geometric)
    series sum_{k=0}^{max_order} (-t*B)^k, which is an approximation.

    Parameters
    ----------
    op : DiracGraphOperator
        The Dirac graph operator. Must have all eigenvalues nonzero.
    max_order : int, optional
        Maximum series order K.  Defaults to None which triggers the exact
        Cayley-Hamilton inversion.  Set to an integer for a truncated
        Neumann (geometric) series approximation of that order.

    Returns
    -------
    (G_e, is_rational, info)
        G_e : sympy.Matrix
            The electron propagator matrix (exact rational).
        is_rational : bool
            True if all entries are rational.
        info : dict
            Diagnostic information including:
            - 'orders_used': number of matrix powers computed
            - 'max_order': the K value used (or 'cayley-hamilton')
            - 'N': matrix dimension
            - 'method': 'cayley-hamilton' or 'truncated-neumann'
            - 'converged': True if the result is exact
    """
    N = op.N

    # Build Lambda^{-1} as diagonal sympy matrix
    Lambda_inv = sp_zeros(N, N)
    for i in range(N):
        lam = op._eigenvalues[i]
        if lam == 0:
            raise ValueError(
                f"Eigenvalue at index {i} is zero; Lambda is singular. "
                "Neumann series requires invertible diagonal."
            )
        Lambda_inv[i, i] = Rational(1) / lam

    # Build adjacency A as sympy matrix (sparse, 0/1 entries)
    A_sym = sp_zeros(N, N)
    cx = op._adj_csr.tocoo()
    for r, c in zip(cx.row, cx.col):
        if r != c:
            A_sym[r, c] = sp.Integer(1)

    # B = Lambda^{-1} * A  (exact rational)
    B = Lambda_inv * A_sym

    t_val = op.t

    if max_order is not None:
        # Truncated Neumann (geometric) series:
        #   S = sum_{k=0}^{K} (-t)^k B^k
        neg_t_B = (-t_val) * B
        S = sp.eye(N)
        P = sp.eye(N)
        for k in range(1, max_order + 1):
            P = P * neg_t_B
            S = S + P

        G_raw = S * Lambda_inv

        # Simplify with cancel() (exact rational simplification)
        G_e = sp_zeros(N, N)
        for i in range(N):
            for j in range(N):
                G_e[i, j] = sp.cancel(G_raw[i, j])

        is_rational = _check_rational_matrix(G_e)

        info = {
            'orders_used': max_order + 1,
            'max_order': max_order,
            'N': N,
            'method': 'truncated-neumann',
            'converged': False,
            't': str(t_val),
        }
        return G_e, is_rational, info

    # Exact Cayley-Hamilton inversion of M = I + t*B.
    #
    # The characteristic polynomial is p(x) = det(xI - M).
    # By Cayley-Hamilton, p(M) = 0.
    # M^{-1} = -(1/c_0)*(M^{N-1} + c_{N-1}*M^{N-2} + ... + c_1*I)
    # where c_0 = det(-M) * (-1)^N = det(M) (up to sign conventions).
    #
    # We compute the charpoly coefficients and build M^{-1} directly.

    M = sp.eye(N) + t_val * B

    # Get characteristic polynomial: det(x*I - M) = x^N + a_{N-1}x^{N-1} + ... + a_0
    # sympy's charpoly returns coefficients [1, a_{N-1}, ..., a_0] (leading coeff first)
    x = sp.Symbol('x')
    charpoly_coeffs = M.charpoly(x).all_coeffs()
    # charpoly_coeffs[0] = 1 (leading), charpoly_coeffs[N] = a_0 = (-1)^N det(M)
    a_0 = charpoly_coeffs[N]  # constant term

    if a_0 == 0:
        raise ValueError("Matrix M = I + t*B is singular (det = 0).")

    # M^{-1} = -(1/a_0) * (M^{N-1} + a_1*M^{N-2} + ... + a_{N-1}*I)
    # where a_k = charpoly_coeffs[k] (k=0 is leading 1, k=1 is a_1, etc.)
    #
    # Build M_powers: M^0 = I, M^1, ..., M^{N-1}
    M_powers = [sp.eye(N)]
    for k in range(1, N):
        M_powers.append(M_powers[k - 1] * M)

    # Accumulate: sum = M^{N-1} + a_1*M^{N-2} + ... + a_{N-1}*I
    # In terms of charpoly_coeffs: a_k = charpoly_coeffs[k], sum over k=0..N-1
    # of charpoly_coeffs[k] * M^{N-1-k}
    accumulator = sp_zeros(N, N)
    for k in range(N):
        coeff = charpoly_coeffs[k]
        power_idx = N - 1 - k
        if coeff != 0:
            accumulator = accumulator + coeff * M_powers[power_idx]

    # M^{-1} = -(1/a_0) * accumulator
    M_inv_raw = (Rational(-1) / a_0) * accumulator

    # G_e = M^{-1} * Lambda^{-1}
    G_raw = M_inv_raw * Lambda_inv

    # Simplify with cancel() (exact rational simplification)
    G_e = sp_zeros(N, N)
    for i in range(N):
        for j in range(N):
            G_e[i, j] = sp.cancel(G_raw[i, j])

    is_rational = _check_rational_matrix(G_e)

    info = {
        'orders_used': N,
        'max_order': 'cayley-hamilton',
        'N': N,
        'method': 'cayley-hamilton',
        'converged': True,
        't': str(t_val),
    }

    return G_e, is_rational, info


def massive_propagator(
    op: DiracGraphOperator,
    mass: object,
    exact: bool = True
) -> Tuple[object, bool]:
    """Massive electron propagator G_e(m) = (D_graph - m·I)^{-1}.

    Parameters
    ----------
    op : DiracGraphOperator
        The Dirac graph operator.
    mass : rational, int, or float
        Mass parameter m. m=0 gives the massless propagator.
    exact : bool
        Whether to use exact sympy or numpy.

    Returns
    -------
    (propagator, is_rational)
    """
    if exact:
        m = Rational(mass) if not isinstance(mass, sp.Expr) else mass
        M = op.matrix_sympy()
        N = op.N
        M_shifted = M - m * sp.eye(N)
        G = M_shifted.inv()
        is_rational = _check_rational_matrix(G)
        return G, is_rational
    else:
        m_float = float(mass)
        M = op.matrix_numpy()
        N = op.N
        M_shifted = M - m_float * np.eye(N)
        G = np.linalg.inv(M_shifted)
        return G, None


def propagator_spectrum(op: DiracGraphOperator) -> Dict:
    """Compute spectral properties of D_graph.

    Returns
    -------
    dict with keys:
        eigenvalues_np : np.ndarray
            All eigenvalues of D_graph (sorted).
        spectral_gap : float
            Smallest |eigenvalue|.
        condition_number : float
            Condition number κ(D_graph) = σ_max / σ_min.
        bare_eigenvalues : list of Rational
            The diagonal Camporesi-Higuchi eigenvalues (no hopping).
        min_bare : Rational
            Minimum |bare eigenvalue|.
        max_bare : Rational
            Maximum |bare eigenvalue|.
    """
    eigs = op.spectrum_numpy()
    bare = op.eigenvalues
    bare_abs = [abs(e) for e in bare]
    return {
        "eigenvalues_np": eigs,
        "spectral_gap": float(np.min(np.abs(eigs))),
        "condition_number": float(op.condition_number()),
        "bare_eigenvalues": bare,
        "min_bare": min(bare_abs),
        "max_bare": max(bare_abs),
        "n_positive_chi": sum(1 for e in bare if e > 0),
        "n_negative_chi": sum(1 for e in bare if e < 0),
    }


def propagator_pi_free_certificate(
    op: DiracGraphOperator,
    test_mass: Optional[Rational] = None
) -> Dict:
    """Verify that D_graph and G_e are π-free (all entries in ℚ).

    Checks:
    1. All diagonal entries are rational half-integers.
    2. All off-diagonal entries are rational (0 or t).
    3. det(D_graph) is a nonzero rational.
    4. All entries of G_e are rational.

    Parameters
    ----------
    op : DiracGraphOperator
        Operator to certify.
    test_mass : Rational, optional
        If given, also certify G_e(test_mass).

    Returns
    -------
    dict with boolean fields:
        diag_rational, offdiag_rational, det_rational, det_nonzero,
        propagator_rational, all_pass
    """
    result = {}

    # 1. Check diagonal entries
    diag_ok = all(isinstance(e, sp.Rational) for e in op.eigenvalues)
    result["diag_rational"] = diag_ok

    # 2. Check off-diagonal (t must be rational, adjacency is 0/1)
    t_ok = isinstance(op.t, (sp.Rational, int)) or op.t == Rational(0)
    result["offdiag_rational"] = t_ok

    # 3. Exact determinant
    M = op.matrix_sympy()
    det = M.det()
    det_rat = isinstance(det, (sp.Rational, sp.Integer, sp.core.numbers.Zero))
    result["det_rational"] = det_rat
    result["det_value"] = det
    result["det_nonzero"] = (det != 0)

    # 4. Propagator entries
    if det != 0:
        G, is_rat = electron_propagator(op, exact=True)
        result["propagator_rational"] = is_rat
        result["propagator_shape"] = (op.N, op.N)

        # Verify D_graph · G_e = I to machine precision
        M_np = op.matrix_numpy()
        G_np = np.array(G.tolist(), dtype=np.float64)
        identity_residual = np.max(np.abs(M_np @ G_np - np.eye(op.N)))
        result["identity_residual"] = float(identity_residual)
    else:
        result["propagator_rational"] = None
        result["identity_residual"] = None

    # 5. Optional mass test
    if test_mass is not None:
        Gm, is_rat_m = massive_propagator(op, test_mass, exact=True)
        result["massive_propagator_rational"] = is_rat_m
        result["mass_tested"] = float(test_mass)

    result["all_pass"] = all([
        result["diag_rational"],
        result["offdiag_rational"],
        result["det_rational"],
        result["det_nonzero"],
        result.get("propagator_rational", False),
    ])

    return result


def propagator_summary(n_max_list: List[int], t_values: List = None) -> Dict:
    """Comprehensive summary of D_graph and G_e at multiple n_max and t values.

    For each (n_max, t) pair:
    - Matrix size
    - Rank
    - Spectral gap
    - Condition number
    - π-free certificate (exact, for N ≤ 30)
    - Structural observations

    Parameters
    ----------
    n_max_list : list of int
        Which n_max values to analyze.
    t_values : list, optional
        Which hopping values to try. Default: [0, 1, Rational(-1, 16)].

    Returns
    -------
    dict: keyed by (n_max, t_str)
    """
    if t_values is None:
        t_values = [Rational(0), Rational(1), KAPPA_SCALAR]

    results = {}

    for n_max in n_max_list:
        for t in t_values:
            key = f"n_max={n_max}_t={t}"
            op = DiracGraphOperator(n_max=n_max, t=t)
            spec = propagator_spectrum(op)

            entry = {
                "n_max": n_max,
                "t": str(t),
                "t_float": float(t),
                "N": op.N,
                "n_edges": op.n_edges,
                "sparsity": float(op.sparsity),
                "spectral_gap": spec["spectral_gap"],
                "condition_number": spec["condition_number"],
                "n_positive_chi": spec["n_positive_chi"],
                "n_negative_chi": spec["n_negative_chi"],
                "bare_eigenvalues_str": [str(e) for e in spec["bare_eigenvalues"]],
            }

            # Exact certificate only for small matrices
            if op.N <= 30:
                cert = propagator_pi_free_certificate(op)
                entry["det_value"] = str(cert["det_value"])
                entry["det_nonzero"] = bool(cert["det_nonzero"])
                entry["propagator_rational"] = cert.get("propagator_rational")
                entry["identity_residual"] = cert.get("identity_residual")
                entry["pi_free_all_pass"] = cert.get("all_pass")

            results[key] = entry

    return results


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _check_rational_matrix(G: Matrix) -> bool:
    """Return True if all entries of the sympy matrix G are rational."""
    for entry in G:
        if not isinstance(entry, (sp.Rational, sp.Integer,
                                  sp.core.numbers.Zero, sp.core.numbers.One,
                                  sp.core.numbers.NegativeOne)):
            # Try to simplify and check again
            s = sp.simplify(entry)
            if not isinstance(s, (sp.Rational, sp.Integer,
                                  sp.core.numbers.Zero,
                                  sp.core.numbers.NegativeOne,
                                  sp.core.numbers.One)):
                return False
    return True


# ---------------------------------------------------------------------------
# Convenience: node-indexed propagator access
# ---------------------------------------------------------------------------

def propagator_entry(
    G: Matrix,
    labels: List[DiracLabel],
    state_a: DiracLabel,
    state_b: DiracLabel
) -> sp.Expr:
    """Return G[a, b] by DiracLabel lookup.

    Parameters
    ----------
    G : sympy Matrix
        The propagator matrix.
    labels : list of DiracLabel
        Node labels (from DiracGraphOperator.labels).
    state_a, state_b : DiracLabel
        The bra and ket states.

    Returns
    -------
    sympy expression G[i_a, i_b].
    """
    label_index = {lab: i for i, lab in enumerate(labels)}
    if state_a not in label_index:
        raise KeyError(f"state_a {state_a} not in lattice")
    if state_b not in label_index:
        raise KeyError(f"state_b {state_b} not in lattice")
    return G[label_index[state_a], label_index[state_b]]


# ---------------------------------------------------------------------------
# Sparsity helpers (attach to DiracGraphOperator)
# ---------------------------------------------------------------------------

@property
def _n_edges(self) -> int:
    return self._lattice.num_edges


@property
def _sparsity(self) -> float:
    return self._lattice.sparsity()


DiracGraphOperator.n_edges = property(lambda self: self._lattice.num_edges)
DiracGraphOperator.sparsity = property(lambda self: self._lattice.sparsity())
