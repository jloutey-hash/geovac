"""Connes distance d_{D_{n_max}}(phi, psi) on the GeoVac truncated operator
system O_{n_max} = P_{n_max} C^infty(S^3) P_{n_max}.

This module implements WH1 round-3: the metric step of the Connes-vS
spectral truncation framework on the Fock-projected S^3, computing the
Connes distance between pure node-evaluation states phi_v(x) := <v|x|v>
on the operator system O_{n_max} of round 2 (geovac/operator_system.py).

Reference: A. Connes & W. D. van Suijlekom, "Spectral Truncations in
Noncommutative Geometry and Operator Systems," CMP 383 (2021),
arXiv:2004.14115. Equation (2), line 1090 of debug/connes_vs_2004_14115.txt:

    d(phi, psi) = sup { |phi(x) - psi(x)| : x in E, ||[D, x]||_op <= 1 }.

For pure node-evaluation states phi_v(x) = <v|x|v> = Tr(E_v x), with
E_v = |v><v| the rank-1 projector onto the basis vector v, the objective
becomes a real linear functional on x:

    phi_v(x) - phi_w(x) = Tr((E_v - E_w) x).

The SDP formulation in cvxpy (with x parameterized as x = sum_k c_k M_k):

    maximize  | Tr((E_v - E_w) x) |
    subject to  x in O (linear membership constraint, automatic via parameterization)
                x = x^*  (Hermitian)
                ||[D, x]||_op <= 1  (operator-norm constraint, SDP form)

The operator-norm constraint ||[D,x]||_op <= 1 is equivalent to the LMI

    [[I, [D,x]], [[D,x]^*, I]] >> 0,

which is what cvxpy lowers to under the hood when you write
``cp.norm(commutator, 2) <= 1``.

The Dirac operator
==================

For the round-3 computation we need a Dirac-like self-adjoint operator on
H_{n_max}. We use the simplest faithful realization:

    D = diag(|lambda_n|)  where  |lambda_n| = n + 3/2  (Camporesi-Higuchi)

acting on H_{n_max} by promotion to the basis {|n,l,m>}: each basis vector
gets eigenvalue n + 3/2. This is the SCALAR analog of the Dirac operator
on S^3 (the spinor lift would double the Hilbert space and act non-trivially
on the spin index, but for the algebraic structure relevant to the metric
it suffices to use the eigenvalue magnitudes on the scalar Fock basis).

The convention is documented in the dataclass DiracProxy below. Different
conventions (D = n shifted by a constant, D scaled) change numerical values
of d_Connes by a fixed prefactor / shift but do not change the structural
pattern (monotonicity in graph distance, ratio of distances, etc.).

Placeholder caveat
==================

The TruncatedOperatorSystem of round 2 uses a placeholder for the
SO(4) radial overlap I_{nNn'}^{lLl'}. The propagation number of round 2 is
proven invariant under this placeholder (it is a structural invariant of
the SO(4) selection rules). The Connes distance, by contrast, is
value-sensitive: replacing the placeholder by the genuine Avery-Wen-Avery
3-Y integral on S^3 will change the numerical distance values.

What is invariant under reasonable placeholders:

  - d(v, v) = 0 (positivity at the diagonal)
  - d(v, w) = d(w, v) (symmetry)
  - Triangle inequality d(u, w) <= d(u, v) + d(v, w)
  - Whether d_Connes is monotone in d_graph (the qualitative metric structure)

What is convention-dependent:

  - The absolute numerical value of d_Connes(v, w)
  - The exact ratio d_Connes / d_graph

The memo wh1_round3_connes_distance_memo.md frames each result honestly
as either "placeholder-robust structural" or "placeholder-dependent numerical".

Public API
==========

The main entry points are:

  - compute_connes_distance(op_sys, v, w, D=None) -> float
  - compute_distance_matrix(op_sys, D=None) -> ndarray
  - fock_graph_distance(v, w) -> int        (combinatorial baseline)

with helpers DiracProxy and a few SDP-internal utilities.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple, Union

import numpy as np

try:
    import cvxpy as cp
    HAS_CVXPY = True
except ImportError:  # pragma: no cover - optional dep
    HAS_CVXPY = False
    cp = None  # type: ignore

from geovac.operator_system import HyperLabel, TruncatedOperatorSystem


# ---------------------------------------------------------------------------
# Dirac operator on H_{n_max}
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DiracProxy:
    """Dirac-proxy operator D on H_{n_max}.

    The Camporesi-Higuchi spectrum of the Dirac operator on the round S^3
    is |lambda_n| = n + 3/2 with multiplicity g_n^Dirac = 2(n+1)(n+2). The
    natural scalar promotion to the Fock basis {|n,l,m>} would set
    D|n,l,m> = (n + 3/2)|n,l,m>; but this gives a *diagonal* D in the
    multiplier basis, and ker([D, .]) cap O is then very large (it
    contains every multiplier M_{N, L=0, M=0} that is also diagonal-by-
    shell). The resulting Connes distance is +infinity on most pairs.

    To get a finite, well-posed metric on every pair (v, w), the Dirac
    needs to break the (n, l, m)-diagonal structure of the multipliers --
    i.e. it must have OFF-diagonal couplings in the (n, l, m) basis.

    Three presets are provided:

      - mode = "shell_scalar": D = diag(n + shift).
          Diagonal, n-shell-degenerate. Gives +infinity on most pairs.
          This is the "truthful" Camporesi-Higuchi scalar Dirac and the
          honest diagnostic of the degeneracy obstruction.

      - mode = "nondegenerate": D = diag(n + l_weight*l + m_weight*m).
          Diagonal, non-degenerate spectrum. Still has very large
          ker([D, .]) cap O because the diagonal-by-shell multipliers
          all commute with any diagonal D. Also gives +infinity on
          most pairs in the present truncation.

      - mode = "offdiag" (DEFAULT for round 3): D = D_diag + alpha * H,
          where H is the Hermitian part of a small symmetric perturbation
          breaking the (n, l, m) diagonal structure. We use H built from
          the n-ladder couplings of the GeoVac Fock graph: H[i, j] = 1
          if |n_i - n_j| = 1 and (l_i, m_i) = (l_j, m_j), zero otherwise.
          This is the natural "raising/lowering" coupling on the SO(4)
          principal QN, which is the closest off-diagonal analog to the
          d/dx Dirac on S^1 in Connes-vS's Toeplitz example.

    The "offdiag" choice is non-degenerate AND breaks commutation with
    the M_{N, 0, 0} multipliers. With it, ker([D, .]) cap O is reduced to
    the identity span only, and d_D(v, w) becomes finite for (almost)
    every pair (v, w). Numerical values depend on alpha; structural
    pattern (positivity, triangle inequality, monotonicity in graph
    distance) does not.
    """

    n_max: int
    shift: float = 0.0
    mode: str = "offdiag"
    l_weight: float = 0.1
    m_weight: float = 0.01
    alpha: float = 1.0

    def matrix(self, basis: List[HyperLabel]) -> np.ndarray:
        """Build the Dirac matrix on the given basis."""
        N = len(basis)
        if self.mode == "shell_scalar":
            diag = np.array(
                [b.n + self.shift for b in basis], dtype=np.complex128,
            )
            return np.diag(diag)
        if self.mode == "nondegenerate":
            diag = np.array(
                [
                    b.n + self.shift + self.l_weight * b.l + self.m_weight * b.m
                    for b in basis
                ],
                dtype=np.complex128,
            )
            return np.diag(diag)
        if self.mode == "offdiag":
            diag = np.array(
                [
                    b.n + self.shift + self.l_weight * b.l + self.m_weight * b.m
                    for b in basis
                ],
                dtype=np.complex128,
            )
            D = np.diag(diag)
            # Off-diagonal SO(4) raising/lowering couplings: hop between
            # (n, l, m) and (n', l', m') whenever the SO(4) ladder
            # selection rules of an E1-style "Dirac-like" derivative
            # operator are satisfied. We use the rule
            #    |Delta n| = 1, |Delta l| = 1, |Delta m| <= 1.
            # This is the Camporesi-Higuchi spinor-Dirac off-diagonal
            # support pattern projected onto the scalar Fock basis.
            #
            # The resulting D is non-diagonal in m, l, and n, so its
            # commutator [D, .] kills only the identity multiplier
            # M_{1, 0, 0} (and any operator proportional to I); ker([D, .])
            # cap O = C * I, the smallest possible kernel.
            for i, bi in enumerate(basis):
                for j, bj in enumerate(basis):
                    if i == j:
                        continue
                    if (
                        abs(bi.n - bj.n) == 1
                        and abs(bi.l - bj.l) == 1
                        and abs(bi.m - bj.m) <= 1
                    ):
                        D[i, j] = self.alpha
            return D
        raise ValueError(f"unknown DiracProxy mode {self.mode!r}")


def default_dirac(op_sys: TruncatedOperatorSystem) -> np.ndarray:
    """Default Dirac proxy: ``DiracProxy(n_max, mode='offdiag')``.

    Diagonal n + 0.1*l + 0.01*m plus n-ladder unit off-diagonal couplings on
    (l, m)-conserving pairs. This breaks the commutation of the diagonal-
    Dirac with the axisymmetric M_{N,0,0} multipliers, giving a finite
    Connes distance on (almost) every pair (v, w).
    """
    return DiracProxy(n_max=op_sys.n_max).matrix(op_sys.basis)


# ---------------------------------------------------------------------------
# Combinatorial Fock-graph distance (baseline)
# ---------------------------------------------------------------------------


def fock_graph_distance(v: HyperLabel, w: HyperLabel) -> int:
    """Combinatorial graph distance between two Fock-graph nodes.

    We use the simplest "shell ladder" graph metric on the (n, l, m) labels:

        d_graph(v, w) = |n_v - n_w| + |l_v - l_w| + |m_v - m_w|

    which is the L^1 distance in the (n, l, m) coordinate. This is the
    natural Cayley-graph distance under the SO(4) ladder generators
    Delta_n = +/- 1, Delta_l = +/- 1, Delta_m = +/- 1.

    Other natural choices include:

      - d = max(|Dn|, |Dl|, |Dm|)  (Chebyshev / shell-ladder)
      - d = |n_v - n_w|              (only the principal-QN ladder)

    They are all monotone in n-distance for nodes in different shells, and
    the qualitative comparisons in the memo do not depend on the specific
    choice.
    """
    return abs(v.n - w.n) + abs(v.l - w.l) + abs(v.m - w.m)


# ---------------------------------------------------------------------------
# SDP for one Connes distance
# ---------------------------------------------------------------------------


def _build_basis_real_imag(
    op_sys: TruncatedOperatorSystem,
) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Build a Hermitian-spanning real basis for O.

    The multiplier matrices {M_k} are complex Hermitian-spanning generators
    of O over the COMPLEX field. To parameterize Hermitian elements x in O
    over the REALS, we expand each M_k into its Hermitian/anti-Hermitian
    parts:

        M_k = H_k + i K_k,  H_k = (M_k + M_k^*)/2,  K_k = (M_k - M_k^*)/(2i)

    Then any Hermitian x in O can be written as x = sum_k (a_k H_k + b_k (i K_k * i))
    with real a_k, b_k. Wait, more carefully: x in O means x = sum c_k M_k for
    complex c_k. With c_k = a_k + i b_k:

        x = sum_k (a_k H_k - b_k K_k) + i sum_k (a_k K_k + b_k H_k)

    For x to be Hermitian: imag part = 0, i.e.
        sum_k (a_k K_k + b_k H_k) = 0.

    A simpler equivalent route: the Hermitian elements of O are
    O_h := { x in O : x = x^* }. Since O is *-closed (round 2), O_h is a
    real-linear subspace of dim_R = 2 dim_C(O) - dim_C(O cap O^*) ... but
    actually since O is *-closed, dim_R(O_h) = dim_C(O).

    The cleanest parameterization: for each generator M_k, the Hermitization
    H_k = (M_k + M_k^*)/2 is in O (since O is *-closed) and Hermitian, and
    (i K_k) = (M_k - M_k^*)/2 / i = (M_k - M_k^*) / (2i)... hmm, this is
    getting tangled.

    ALTERNATIVE (used here): just parameterize x by its complex coefficients
    c_k in front of M_k directly, then add a Hermitian constraint x = x^* in
    cvxpy as a linear constraint. This is what we do.

    Return: (real_basis, imag_basis) where real_basis[k] = M_k and
    imag_basis[k] = M_k as numpy arrays. (Trivial pass-through; the
    Hermiticity is handled in the SDP layer.)
    """
    return list(op_sys.multiplier_matrices), list(op_sys.multiplier_matrices)


def compute_connes_distance(
    op_sys: TruncatedOperatorSystem,
    v: Union[HyperLabel, int],
    w: Union[HyperLabel, int],
    D: Optional[np.ndarray] = None,
    *,
    solver: str = "SCS",
    solver_kwargs: Optional[dict] = None,
    verbose: bool = False,
) -> float:
    """Compute d_D(phi_v, phi_w) on the truncated operator system.

    Parameters
    ----------
    op_sys : TruncatedOperatorSystem
        The truncated operator system O_{n_max} from round 2.
    v, w : HyperLabel or int
        The two basis vectors (states phi_v(x) = <v|x|v>). Can be passed
        as HyperLabel(n, l, m) or as an integer index into op_sys.basis.
    D : np.ndarray, optional
        The Dirac operator on H_{n_max} as a Hermitian (dim_H, dim_H) matrix.
        Default: DiracProxy(n_max).matrix(basis) = diag(n + 3/2).
    solver : str
        cvxpy SDP solver (default "SCS"; CLARABEL also works).
    solver_kwargs : dict, optional
        Extra keyword arguments to pass to ``problem.solve``.
    verbose : bool
        Print solver progress.

    Returns
    -------
    distance : float
        The Connes distance d_D(phi_v, phi_w) >= 0. Returns 0.0 exactly if
        v == w.

    Raises
    ------
    ImportError if cvxpy is not installed.
    """
    if not HAS_CVXPY:
        raise ImportError(
            "cvxpy is required for Connes distance computation. "
            "Install with `pip install cvxpy`."
        )
    # Resolve integer indices to HyperLabels.
    if isinstance(v, int):
        v_idx = v
    else:
        v_idx = op_sys.basis.index(v)
    if isinstance(w, int):
        w_idx = w
    else:
        w_idx = op_sys.basis.index(w)
    if v_idx == w_idx:
        return 0.0

    if D is None:
        D = default_dirac(op_sys)
    D = np.asarray(D, dtype=np.complex128)
    N = op_sys.dim_H
    if D.shape != (N, N):
        raise ValueError(f"D shape {D.shape} != ({N}, {N})")

    K = len(op_sys.multiplier_matrices)
    # Stack the K multipliers as a (K, N, N) complex tensor for vectorized assembly.
    M_stack = np.stack(op_sys.multiplier_matrices, axis=0)

    # Cvxpy variables: complex coefficients c_k = a_k + i b_k.
    # We use a complex variable directly (cvxpy handles this).
    c = cp.Variable(K, complex=True)
    # x = sum_k c_k M_k. Vectorize via reshape: x.flatten() = M_flat @ c
    # where M_flat is (N*N, K) the stacked vec(M_k)'s.
    M_flat = M_stack.transpose(1, 2, 0).reshape(N * N, K)
    x_vec = cp.Constant(M_flat) @ c
    x_expr = cp.reshape(x_vec, (N, N), order="C")
    # x in O is automatic. Hermitian constraint:
    constraints = [x_expr == cp.conj(x_expr.T)]

    # Commutator [D, x] = D @ x - x @ D.
    D_const = cp.Constant(D)
    commutator = D_const @ x_expr - x_expr @ D_const

    # Operator norm constraint ||commutator||_2 <= 1.
    constraints.append(cp.norm(commutator, 2) <= 1.0)

    # Gauge fixing: the SDP is invariant under x -> x + y for any
    # y in ker([D, .]) (any matrix that commutes with D), and the
    # ambient operator system O typically contains MORE than the identity
    # in this kernel -- e.g. M_{N,0,0} multipliers (axisymmetric f's
    # depending only on the polar angle chi) are diagonal-by-shell and
    # therefore commute with any diagonal D.
    #
    # Two cases:
    #   (a) The objective functional E_v - E_w has nonzero inner product
    #       with some y in ker([D, .]) cap O_h: the SDP is unbounded along
    #       this direction (free shifts in objective at zero commutator
    #       cost), so d_D(phi_v, phi_w) = +infinity. Detect and return.
    #   (b) E_v - E_w is HS-orthogonal to ker([D, .]) cap O_h: the SDP is
    #       well-posed once we pin x orthogonally to the kernel.
    kernel_basis = _kernel_of_commutator_in_O(op_sys, D, tol=1e-8)
    E_diff = np.zeros((N, N), dtype=np.complex128)
    E_diff[v_idx, v_idx] = 1.0
    E_diff[w_idx, w_idx] = -1.0
    for K_mat in kernel_basis:
        proj = np.real(np.trace(K_mat.conj().T @ E_diff))
        if abs(proj) > 1e-9:
            return float("inf")
    # All kernel elements are orthogonal to E_diff -> add gauge-fix
    # constraints pinning x to ker^perp.
    for K_mat in kernel_basis:
        K_const = cp.Constant(K_mat)
        constraints.append(
            cp.real(cp.trace(cp.conj(K_const.T) @ x_expr)) == 0.0
        )

    # Objective: maximize |<E_v - E_w, x>_HS| = |x[v,v] - x[w,w]|.
    # x[v_idx, v_idx] is the diagonal, real for Hermitian x.
    diff = cp.real(x_expr[v_idx, v_idx]) - cp.real(x_expr[w_idx, w_idx])
    # cvxpy doesn't directly maximize an absolute value (non-DCP).
    # Standard trick: solve two LPs (max diff, max -diff), take the larger.
    best = -np.inf
    for sense in (+1, -1):
        prob = cp.Problem(cp.Maximize(sense * diff), constraints)
        kwargs = dict(solver=solver, verbose=verbose)
        if solver_kwargs:
            kwargs.update(solver_kwargs)
        try:
            prob.solve(**kwargs)
        except cp.SolverError as e:  # pragma: no cover
            raise RuntimeError(f"SDP solver {solver} failed: {e}")
        if prob.status not in ("optimal", "optimal_inaccurate"):
            raise RuntimeError(
                f"SDP did not converge: status={prob.status}, sense={sense}, "
                f"v={v_idx}, w={w_idx}"
            )
        val = float(prob.value)
        if val > best:
            best = val
    return max(best, 0.0)


def compute_distance_matrix(
    op_sys: TruncatedOperatorSystem,
    D: Optional[np.ndarray] = None,
    *,
    solver: str = "SCS",
    solver_kwargs: Optional[dict] = None,
    verbose: bool = False,
    progress: bool = False,
) -> np.ndarray:
    """Compute the full N x N Connes distance matrix on op_sys.

    Returns a real symmetric numpy array of shape (N, N) with zero diagonal.

    Parameters
    ----------
    op_sys : TruncatedOperatorSystem
    D : np.ndarray, optional
        Dirac operator (default: diag(n + 3/2)).
    solver, solver_kwargs, verbose : as in compute_connes_distance.
    progress : bool
        Print which (v, w) pair is being solved.

    Notes
    -----
    For an N x N matrix only the strict upper triangle (N(N-1)/2 pairs) is
    computed; the diagonal is zero by definition and the lower triangle is
    filled by symmetry. Each pair requires two SDPs (max(diff), max(-diff)),
    so the total cost is N(N-1) SDPs.

    At n_max=2 this is 5*4 = 20 SDPs (~5 s with SCS). At n_max=3 this is
    14*13 = 182 SDPs (~minutes).
    """
    N = op_sys.dim_H
    if D is None:
        D = default_dirac(op_sys)

    dist = np.zeros((N, N), dtype=float)
    for v_idx in range(N):
        for w_idx in range(v_idx + 1, N):
            d = compute_connes_distance(
                op_sys, v_idx, w_idx, D=D,
                solver=solver, solver_kwargs=solver_kwargs, verbose=verbose,
            )
            dist[v_idx, w_idx] = d
            dist[w_idx, v_idx] = d
            if progress:
                print(
                    f"  d({v_idx}, {w_idx}) = "
                    f"d({op_sys.basis[v_idx]}, {op_sys.basis[w_idx]}) = "
                    f"{d:.6f}"
                )
    return dist


# ---------------------------------------------------------------------------
# Convenience: graph-distance matrix and pretty labels
# ---------------------------------------------------------------------------


def _kernel_of_commutator_in_O(
    op_sys: TruncatedOperatorSystem,
    D: np.ndarray,
    *,
    tol: float = 1e-8,
) -> List[np.ndarray]:
    """Find an orthonormal HS-basis of { x in O_h : [D, x] = 0 }.

    Strategy:
      1. Build the Hermitian-spanning real basis of O. Since O is *-closed
         (round 2), {(M_k + M_k^*)/2, (M_k - M_k^*)/(2i)} for k=1..K span
         the real Hermitian subspace O_h. Many of these are linearly
         dependent; we extract a basis via SVD on the Hermitian-vec stack.
      2. For each basis element B_j of O_h, compute [D, B_j].
      3. The kernel is the null space of the linear map B -> [D, B] on O_h.
    """
    M_list = op_sys.multiplier_matrices
    K = len(M_list)
    N = op_sys.dim_H
    if K == 0:
        return []

    # Build real Hermitian-spanning generators: H_k = (M_k + M_k^*)/2 and
    # K_k = (M_k - M_k^*)/(2i). Both are Hermitian.
    herm_generators: List[np.ndarray] = []
    for M in M_list:
        H = (M + M.conj().T) / 2
        K_part = (M - M.conj().T) / (2j)
        herm_generators.append(H)
        herm_generators.append(K_part)

    # Filter zero / nearly-zero generators.
    herm_generators = [
        H for H in herm_generators if np.linalg.norm(H) > 1e-14
    ]

    # Stack as columns of a (N^2, num_gen) complex matrix; extract a
    # complex-rank basis via SVD.
    if not herm_generators:
        return []
    cols = np.zeros(
        (N * N, len(herm_generators)), dtype=np.complex128,
    )
    for j, H in enumerate(herm_generators):
        cols[:, j] = H.reshape(-1)
    U, S, _ = np.linalg.svd(cols, full_matrices=False)
    if S.size == 0:
        return []
    rank = int(np.sum(S > tol * max(S.max(), 1.0)))
    # The first `rank` columns of U are an orthonormal HS-basis of O_h
    # (within complex M_N(C), but the basis vectors will be Hermitian
    # because the generators were Hermitian and the column space is the
    # Hermitian subspace).
    O_h_basis = [U[:, j].reshape((N, N)) for j in range(rank)]
    # Re-Hermitize to kill any numerical anti-Hermitian residual.
    O_h_basis = [(B + B.conj().T) / 2 for B in O_h_basis]

    # Now compute the kernel of B -> [D, B] on O_h.
    # Build the matrix L of shape (N^2, rank) where column j is vec([D, B_j]).
    L_cols = np.zeros((N * N, rank), dtype=np.complex128)
    for j, B in enumerate(O_h_basis):
        L_cols[:, j] = (D @ B - B @ D).reshape(-1)
    # Null space of L (treating columns as the rank-many B-coefficients).
    # We want { c in R^rank : L c = 0 }, then kernel basis vectors are
    # B(c) = sum_j c_j B_j. Use SVD null space.
    U_L, S_L, Vh_L = np.linalg.svd(L_cols, full_matrices=True)
    if S_L.size == 0:
        kernel_dim = rank
        null_basis = np.eye(rank, dtype=np.complex128)
    else:
        cutoff = tol * max(S_L.max(), 1.0)
        nonzero_count = int(np.sum(S_L > cutoff))
        kernel_dim = rank - nonzero_count
        # Right-singular vectors corresponding to zero singular values.
        null_basis = Vh_L[nonzero_count:, :].conj().T  # shape (rank, kernel_dim)

    if kernel_dim == 0:
        return []

    # Reconstruct the kernel matrices in M_N(C).
    kernel_matrices = []
    for k in range(kernel_dim):
        coeffs = null_basis[:, k]
        K_mat = sum(coeffs[j] * O_h_basis[j] for j in range(rank))
        # Hermitize and normalize to unit HS norm.
        K_mat = (K_mat + K_mat.conj().T) / 2
        nrm = np.linalg.norm(K_mat)
        if nrm > 1e-14:
            K_mat = K_mat / nrm
            kernel_matrices.append(K_mat)
    return kernel_matrices


def graph_distance_matrix(op_sys: TruncatedOperatorSystem) -> np.ndarray:
    """Combinatorial Fock-graph distance matrix.

    d_graph(v, w) = |Dn| + |Dl| + |Dm| (L^1 distance on (n, l, m) labels).
    """
    N = op_sys.dim_H
    dist = np.zeros((N, N), dtype=int)
    for i in range(N):
        for j in range(N):
            dist[i, j] = fock_graph_distance(op_sys.basis[i], op_sys.basis[j])
    return dist


def basis_label_strings(op_sys: TruncatedOperatorSystem) -> List[str]:
    """Pretty label strings ``|n,l,m>`` for each basis vector."""
    return [f"|{b.n},{b.l},{b.m:+d}>" for b in op_sys.basis]
