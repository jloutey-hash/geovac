"""Krein-positive state space — the state-side complement of the L3a-1 substrate.

Sprint L3b first move (2026-05-17): construct the state-space-level
Krein-positive cone on a Lorentzian truncated operator system on a
Krein space.  This is the GENUINELY NEW mathematical object motivated by
L3a-1's structural finding that operator-multiplier-level Krein-positivity
is trivial (every scalar multiplier preserves K^+ = J-positive cone in
the chirality-doubled spinor lift).

The L3a-1 finding shifted the Krein-positivity question from the operator
side to the STATE side: a Krein-positive state on the operator system O^L
is a linear functional omega : O^L -> C with omega(a^* J a) >= 0 for all
a in O^L (and omega(I) = 1).  The set of such states is the state-space-
level Krein-positive cone, the natural substrate for the propinquity /
Wasserstein-Kantorovich question downstream of L3a-1.

This module's role
==================

This module provides:

  (i) A computable representation of the state space on O^L: a state
      omega corresponds to a Hermitian positive trace-1 matrix rho on
      the truncated Krein space such that omega(a) = Tr(rho J a) for
      Krein-positive states (where J is the Krein fundamental symmetry).

  (ii) Krein-positivity checks: omega(a^* J a) >= 0 verified pointwise
       on a finite sample of generators a in O^L.

  (iii) A Wasserstein-Kantorovich-style distance via cvxpy SDP between
        two pure node-evaluation states under the operator system's
        Lipschitz seminorm.  This is the foundation analog of Paper 38's
        L4 (Berezin reconstruction) -- it gives a distance on the state
        space that the propinquity bound will eventually compare against.

  (iv) The Wick-rotated Hilbert-space inner product on K^+: at the K^+
       cone, the Krein product <psi, J phi> reduces to a positive-definite
       inner product, so K^+ carries a genuine Hilbert-space structure
       on which standard propinquity machinery applies.

Honest scope
============

This is a FOUNDATION module.  It does not:
  - Prove the propinquity bound (that's the L1'-L5 program).
  - Compute the full Wasserstein-Kantorovich completion of the state
    space (Paper 38 §VIII open problem on full state space).
  - Handle Krein-indefinite states (non-positive functionals); we
    restrict to the positive cone by construction.

What it provides is the computable substrate on which propinquity-
related distances can be evaluated, with cvxpy SDP infrastructure for
pure-state distance computation analogous to compute_connes_distance
on Riemannian truncated operator systems.

API
===

  KreinPositiveStateSpace(op_sys)
      .J                         : fundamental symmetry (dim_K, dim_K)
      .pure_state_density(idx)  : projector |v><v| for v = basis[idx]
      .is_krein_positive_state(rho) : check rho satisfies the
                                       Krein-positivity inequality
                                       on the multipliers of op_sys
      .wasserstein_distance_pure(idx_v, idx_w, D)
                                  : SDP for d_D(omega_v, omega_w)
      .restrict_to_K_plus(rho)   : project rho to the K^+ block
      .K_plus_dim                : dim K^+ = dim_K / 2 (chirality-block
                                    convention)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple, Union

import numpy as np

try:
    import cvxpy as cp
    HAS_CVXPY = True
except ImportError:  # pragma: no cover
    HAS_CVXPY = False


# ---------------------------------------------------------------------------
# State space wrapper
# ---------------------------------------------------------------------------


@dataclass
class KreinPositiveStateSpace:
    """State space of a Lorentzian truncated operator system on a Krein space.

    Parameters
    ----------
    op_sys : LorentzianTruncatedOperatorSystem or
             CompactTemporalTruncatedOperatorSystem
        Either grid- or compact-temporal version.  Must have:
          - .krein with .J, .dim
          - .basis_spatial
          - .multiplier_matrices

    Attributes
    ----------
    op_sys : the operator system
    J : (dim_K, dim_K) Krein fundamental symmetry
    dim_K : int
    K_plus_eigvals, K_plus_eigvecs : J = U diag(lam) U^dagger
        lam = +/- 1; K^+ basis = eigenvectors for +1 eigenvalue
    K_plus_dim : int
    K_minus_dim : int
    P_plus, P_minus : projectors onto K^+, K^-
    """

    op_sys: object  # duck type
    J: np.ndarray = None
    dim_K: int = 0
    K_plus_eigvecs: np.ndarray = None
    K_minus_eigvecs: np.ndarray = None
    K_plus_dim: int = 0
    K_minus_dim: int = 0
    P_plus: np.ndarray = None
    P_minus: np.ndarray = None

    def __post_init__(self) -> None:
        self.J = self.op_sys.krein.J
        self.dim_K = self.op_sys.krein.dim
        # Diagonalize J: it should be Hermitian with eigenvalues +/- 1
        # (or +1 in the chirality-swap convention from L2-B).
        eigvals, eigvecs = np.linalg.eigh(self.J)
        # Round eigenvalues to +/- 1 (they should be exactly +1 or -1)
        rounded = np.round(eigvals)
        # K^+ = +1 eigenspace; K^- = -1 eigenspace
        plus_mask = np.isclose(rounded, +1.0)
        minus_mask = np.isclose(rounded, -1.0)
        self.K_plus_eigvecs = eigvecs[:, plus_mask]
        self.K_minus_eigvecs = eigvecs[:, minus_mask]
        self.K_plus_dim = int(np.sum(plus_mask))
        self.K_minus_dim = int(np.sum(minus_mask))

        # Projectors
        self.P_plus = self.K_plus_eigvecs @ self.K_plus_eigvecs.conj().T
        self.P_minus = self.K_minus_eigvecs @ self.K_minus_eigvecs.conj().T

    # -----------------------------------------------------------------
    # Pure-state densities
    # -----------------------------------------------------------------

    def pure_state_density(
        self, vec_or_idx: Union[int, np.ndarray]
    ) -> np.ndarray:
        """Return rho = |v><v| for a pure state.

        Parameters
        ----------
        vec_or_idx : int or array
            If int: use the canonical basis vector e_i of dim_K.
            If array: use the given (normalized) vector.

        Returns
        -------
        rho : (dim_K, dim_K) Hermitian, positive, trace-1.
        """
        if isinstance(vec_or_idx, (int, np.integer)):
            if not (0 <= int(vec_or_idx) < self.dim_K):
                raise ValueError(
                    f"index {vec_or_idx} out of range [0, {self.dim_K})"
                )
            v = np.zeros(self.dim_K, dtype=np.complex128)
            v[int(vec_or_idx)] = 1.0
        else:
            v = np.asarray(vec_or_idx, dtype=np.complex128).reshape(-1)
            norm = np.linalg.norm(v)
            if norm < 1e-30:
                raise ValueError("vector has zero norm")
            v = v / norm
        return np.outer(v, v.conj())

    # -----------------------------------------------------------------
    # Krein-positivity checks
    # -----------------------------------------------------------------

    def krein_expectation(
        self, rho: np.ndarray, a: np.ndarray
    ) -> complex:
        """Compute omega_rho(a) = Tr(rho J a) (Krein-twisted expectation).

        Convention: a Krein-positive state on O^L is one for which
        omega(a^* J a) >= 0 for every a.  We use omega(a) = Tr(rho J a)
        because in the Krein-positive cone (rho supported on K^+)
        this reduces to the standard Hilbert-space expectation Tr(rho a).
        """
        return complex(np.trace(rho @ self.J @ a))

    def hilbert_expectation(
        self, rho: np.ndarray, a: np.ndarray
    ) -> complex:
        """Standard expectation omega_rho(a) = Tr(rho a)."""
        return complex(np.trace(rho @ a))

    def is_krein_positive_state(
        self, rho: np.ndarray, *, sample_size: Optional[int] = None,
        tol: float = -1e-10
    ) -> Tuple[bool, dict]:
        """Check rho satisfies the Krein-positivity inequality.

        For every generator M_i in op_sys.multiplier_matrices, check that

            omega(M_i^* J M_i) = Tr(rho J M_i^* J M_i) >= 0.

        Note: this is necessary but not sufficient (the full state-
        positivity requires the inequality on every a in O^L, not just
        on generators).  For a finite-dimensional operator system, the
        generators span O^L linearly, and checking at the generator
        level is the natural finite-cutoff version.

        Parameters
        ----------
        rho : (dim_K, dim_K) Hermitian
        sample_size : if not None, only check first `sample_size` generators
        tol : threshold for treating real(value) < 0 as a violation
            (default -1e-10 allows machine-precision noise)

        Returns
        -------
        (ok, details_dict)
        """
        if rho.shape != (self.dim_K, self.dim_K):
            raise ValueError(
                f"rho shape {rho.shape} != ({self.dim_K}, {self.dim_K})"
            )

        gens = self.op_sys.multiplier_matrices
        if sample_size is not None:
            gens = gens[:sample_size]

        n_total = len(gens)
        n_violated = 0
        min_real = float("inf")
        max_imag = 0.0
        for M in gens:
            val = np.trace(rho @ self.J @ M.conj().T @ self.J @ M)
            re_v = float(np.real(val))
            im_v = float(np.abs(np.imag(val)))
            if re_v < tol:
                n_violated += 1
            min_real = min(min_real, re_v)
            max_imag = max(max_imag, im_v)

        ok = (n_violated == 0)
        details = {
            "n_total": n_total,
            "n_violated": n_violated,
            "min_real_value": min_real,
            "max_imag_value": max_imag,
            "tol": tol,
        }
        return ok, details

    # -----------------------------------------------------------------
    # Restrict to K^+
    # -----------------------------------------------------------------

    def restrict_to_K_plus(self, rho: np.ndarray) -> np.ndarray:
        """Project rho onto K^+ via P_plus rho P_plus.

        Returns the K^+ block as a (K_plus_dim, K_plus_dim) matrix in
        the K^+ eigenbasis (i.e., expressed in the K_plus_eigvecs basis).
        """
        if rho.shape != (self.dim_K, self.dim_K):
            raise ValueError(
                f"rho shape {rho.shape} != ({self.dim_K}, {self.dim_K})"
            )
        # Expressed in the K^+ basis:
        U_plus = self.K_plus_eigvecs
        rho_plus = U_plus.conj().T @ rho @ U_plus
        return rho_plus

    def K_plus_inner_product(
        self, psi: np.ndarray, phi: np.ndarray
    ) -> complex:
        """The K^+ Hilbert-space inner product on K^+ vectors.

        On K^+ (J psi = psi), the Krein product <psi, J phi> equals
        <psi, phi>, so K^+ inherits a genuine Hilbert space inner
        product.  This routine projects psi, phi onto K^+ and returns
        the standard inner product on the projected vectors.
        """
        if psi.shape != (self.dim_K,) or phi.shape != (self.dim_K,):
            raise ValueError(
                f"shapes must be ({self.dim_K},); "
                f"got {psi.shape}, {phi.shape}"
            )
        psi_plus = self.P_plus @ psi
        phi_plus = self.P_plus @ phi
        return complex(np.vdot(psi_plus, phi_plus))

    # -----------------------------------------------------------------
    # SDP-based Wasserstein-Kantorovich distance on pure node states
    # -----------------------------------------------------------------

    def wasserstein_distance_pure(
        self,
        idx_v: int,
        idx_w: int,
        D: np.ndarray,
        *,
        solver: str = "SCS",
        solver_kwargs: Optional[dict] = None,
        verbose: bool = False,
    ) -> float:
        """Compute the SDP-based distance d_D(omega_v, omega_w) on
        pure node states.

        This is the Lorentzian / Krein-positive analog of
        geovac.connes_distance.compute_connes_distance.  The distance is

            d_D(omega_v, omega_w) := sup { |omega_v(a) - omega_w(a)|
                                            : a in O^L,
                                              ||[D, a]||_op <= 1 },

        with omega_v(a) = <e_v, a e_v> for canonical basis e_v.

        This is a foundation analog: it provides a computable distance
        on the state space that the propinquity bound will eventually
        compare against.  It is NOT the full propinquity / Latremoliere
        construction; that requires the L4 / L5 lemmas.

        Honest scope
        ------------
        The distance returned is positive-definite for v != w in the
        generic case (when the operator system is "rich enough" to
        separate states).  Returns +infinity if the operator system
        contains an element commuting with D but distinguishing
        omega_v from omega_w (gauge-fixing failure).

        Parameters
        ----------
        idx_v, idx_w : int
            Indices into the standard basis of K (0 <= idx < dim_K).
        D : np.ndarray
            (dim_K, dim_K) Hermitian Dirac proxy.
        solver : str
        solver_kwargs : dict
        verbose : bool

        Returns
        -------
        distance : float >= 0 or +inf.
        """
        if not HAS_CVXPY:
            raise ImportError(
                "cvxpy is required.  Install with `pip install cvxpy`."
            )
        if not (0 <= idx_v < self.dim_K and 0 <= idx_w < self.dim_K):
            raise ValueError(
                f"indices {idx_v}, {idx_w} out of range [0, {self.dim_K})"
            )
        if idx_v == idx_w:
            return 0.0
        if D.shape != (self.dim_K, self.dim_K):
            raise ValueError(
                f"D shape {D.shape} != ({self.dim_K}, {self.dim_K})"
            )

        gens = self.op_sys.multiplier_matrices
        K = len(gens)
        N = self.dim_K
        if K == 0:
            return 0.0

        # Stack generators
        M_stack = np.stack(gens, axis=0)  # (K, N, N)
        M_flat = M_stack.transpose(1, 2, 0).reshape(N * N, K)

        c = cp.Variable(K, complex=True)
        x_vec = cp.Constant(M_flat) @ c
        x_expr = cp.reshape(x_vec, (N, N), order="C")
        constraints = [x_expr == cp.conj(x_expr.T)]

        D_const = cp.Constant(D)
        commutator = D_const @ x_expr - x_expr @ D_const
        constraints.append(cp.norm(commutator, 2) <= 1.0)

        # Gauge fixing: check if E_v - E_w is HS-orthogonal to ker([D, .]) in O
        E_diff = np.zeros((N, N), dtype=np.complex128)
        E_diff[idx_v, idx_v] = 1.0
        E_diff[idx_w, idx_w] = -1.0

        # Find kernel of [D, .] on the span of generators (numerical SVD)
        commutator_basis_cols = []
        for M in gens:
            commutator_basis_cols.append(
                (D @ M - M @ D).reshape(-1)
            )
        cb_mat = np.array(commutator_basis_cols).T  # (N^2, K)
        # Generator coefficients producing zero commutator:
        # null space of cb_mat is the kernel of [D, .] on span O.
        # We do SVD and identify the right singular vectors with small s.
        U, S, Vh = np.linalg.svd(cb_mat, full_matrices=False)
        tol_svd = 1e-8 * max(S.max() if S.size > 0 else 1.0, 1.0)
        null_mask = S < tol_svd
        # Construct kernel basis matrices from corresponding columns of M_stack
        kernel_basis = []
        for k_idx in range(len(null_mask)):
            if null_mask[k_idx]:
                coeffs = Vh[k_idx, :].conj()
                K_mat = np.zeros((N, N), dtype=np.complex128)
                for j in range(K):
                    K_mat = K_mat + coeffs[j] * gens[j]
                kernel_basis.append(K_mat)

        # Detect free directions for the objective
        for K_mat in kernel_basis:
            proj = np.real(np.trace(K_mat.conj().T @ E_diff))
            if abs(proj) > 1e-9:
                return float("inf")

        for K_mat in kernel_basis:
            K_const = cp.Constant(K_mat)
            constraints.append(
                cp.real(cp.trace(cp.conj(K_const.T) @ x_expr)) == 0.0
            )

        diff = cp.real(x_expr[idx_v, idx_v]) - cp.real(x_expr[idx_w, idx_w])
        prob = cp.Problem(cp.Maximize(diff), constraints)
        kwargs = dict(solver=solver, verbose=verbose)
        if solver_kwargs:
            kwargs.update(solver_kwargs)
        try:
            prob.solve(**kwargs)
        except cp.SolverError as e:
            raise RuntimeError(f"SDP solver {solver} failed: {e}")
        if prob.status not in ("optimal", "optimal_inaccurate"):
            raise RuntimeError(
                f"SDP did not converge: status={prob.status}"
            )
        return max(float(prob.value), 0.0)

    # -----------------------------------------------------------------
    # Sanity / structural checks
    # -----------------------------------------------------------------

    def verify_J_eigendecomp(self, tol: float = 1e-12) -> Tuple[bool, dict]:
        """Verify J = U diag(lam) U^dagger with lam in {+1, -1}.

        The Krein fundamental symmetry on the chirality-doubled space
        has eigenvalues +1 (chirality + sector) and -1 (chirality -
        sector), each with multiplicity dim_K / 2.
        """
        eigvals, _ = np.linalg.eigh(self.J)
        within_pm1 = np.all(
            np.isclose(eigvals, +1, atol=tol)
            | np.isclose(eigvals, -1, atol=tol)
        )
        details = {
            "K_plus_dim": self.K_plus_dim,
            "K_minus_dim": self.K_minus_dim,
            "total": self.K_plus_dim + self.K_minus_dim,
            "dim_K": self.dim_K,
            "eigenvalues_pm1": bool(within_pm1),
            "min_eig": float(eigvals.min()),
            "max_eig": float(eigvals.max()),
        }
        ok = (
            bool(within_pm1)
            and self.K_plus_dim + self.K_minus_dim == self.dim_K
        )
        return ok, details

    def __repr__(self) -> str:
        return (
            f"KreinPositiveStateSpace(dim_K={self.dim_K}, "
            f"K_plus_dim={self.K_plus_dim}, "
            f"K_minus_dim={self.K_minus_dim})"
        )
