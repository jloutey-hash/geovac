"""Joint Berezin reconstruction map B_{n_max, N_t}: C(S^3 x S^1_T) -> O^L (L3b-2 Sub-sprint C).

This module realizes Lemma L4 of the L3b-2 propinquity-convergence proof
on the compact-temporal Lorentzian truncated operator system.  It is the
tensor-product joint extension of `geovac.berezin_reconstruction.BerezinReconstruction`,
combining the SU(2) Peter-Weyl Berezin convolution by the L2 central
spectral Fejer kernel with a U(1) Fourier convolution by the standard
Fejer kernel on S^1_T.

Mathematical setup
==================

Per Sub-sprint A `debug/l3b_2_sub_sprint_A_lichnerowicz_memo.md`, the
compact-temporal operator system is a pure tensor product:

    O^L_{n_max, N_t, T} = span_C { M^{spat}_{N, L, M} (x) g_p(omega_k) }

where:

  - M^{spat}_{N, L, M} is the chirality-doubled scalar 3-Y multiplier
    on H_GV (dim_spatial = 2/3 n_max(n_max+1)(n_max+2));
  - g_p(omega_k) = diag(omega_k^p) is the momentum-polynomial diagonal
    multiplier on C^{N_t}, p = 0, ..., N_t - 1.

Joint Berezin map
-----------------

For f in C^infty(S^3 x S^1_T) with the joint Peter-Weyl x Fourier
expansion

    f(omega, t) = sum_{N, L, M} sum_{q in Z} c_{N L M q} Y^{(3)}_{N L M}(omega) e^{i q t / R_T}

where R_T = T / (2 pi) is the natural circle radius (so e^{i q t / R_T}
has period T), the JOINT Berezin map is

    B_{n_max, N_t, T}(f) := sum_{N <= n_max} sum_{L=0}^{N-1} sum_{|M| <= L}
                              sum_{|q| <= K_max}
                                hat{K}^{SU(2)}_{n_max}(N) *
                                hat{K}^{U(1)}_{N_t}(q) *
                                c_{N L M q} *
                                (M^{spat}_{N L M} (x) M^{temp}_q),

where:

  - hat{K}^{SU(2)}_{n_max}(N) = N / Z^{SU(2)}_{n_max} is the SU(2) Plancherel
    weight from Paper 38 L2 (`geovac.central_fejer_su2`).
  - hat{K}^{U(1)}_{N_t}(q) = max(0, (N_t + 1 - 2|q|) / (N_t + 1)) is the
    standard Fejer-on-S^1 weight (`geovac.central_fejer_compact_temporal`).
  - K_max = floor((N_t - 1) / 2) is the temporal-mode truncation.
  - M^{temp}_q is the momentum-mode-q multiplier: in the natural
    (Fourier-mode) representation, M^{temp}_q is the diagonal indicator
    of momentum eigenvalue q -- specifically diag(delta_{k, q}) in the
    momentum basis (k in {-K_max, ..., +K_max}).  The Berezin map then
    weights mode q by hat{K}^{U(1)}(q) and SUMS over all q in the kept
    range to build the Fourier-projected temporal multiplier.

For pure-tensor symbols f = f_s(omega) * f_t(t) the map factorizes
exactly:

    B^{joint}(f_s * f_t) = B^{SU(2)}(f_s) (x) B^{U(1)}(f_t).

Tensor-product factorization of the four L4 properties
======================================================

Per Sub-sprints A + B, the four L4 properties factor across the SU(2)
and U(1) factors:

(1) Positivity.  If f >= 0 in C(S^3 x S^1_T), then B^{joint}(f) is
    Hermitian PSD.  Mechanism: B^{joint}(f) = P^{joint} (K^{joint} * f)
    P^{joint}, the joint kernel K^{joint} = K^{SU(2)} (x) K^{U(1)} is
    pointwise non-negative (each factor is), and convolution preserves
    non-negativity on a compact group product.  Compression-of-PSD is
    PSD by the standard 〈psi | P g P | psi〉 = ∫ g |P psi|^2 >= 0.

(2) Contractivity.  ||B^{joint}(f)||_op <= ||f||_infty.  Direct corollary
    of Sub-sprint B (joint cb-norm = 2/(n_max+1) <= 1) and the standard
    factorization ||P(K*f)P||_op <= ||K*f||_infty <= ||K||_{L^1} ||f||_infty
    = ||f||_infty (K^{joint} is normalized on the product Haar measure).

(3) Approximate identity.  ||B^{joint}(f) - P^{joint} M_f P^{joint}||_op
    <= gamma^{joint}_{n_max, N_t, T} ||nabla^{joint} f||_infty, with the
    joint gamma rate gamma^{joint} = O(log n_max / n_max + 1 / N_t)
    from L3b first move §5 (factor-additive in either L^1 or L^2 joint
    metric, both -> 0 in the joint limit).

(4) L3 compatibility.  ||[D_L, B^{joint}(f)]||_op <= C_3 * ||nabla^{joint}
    f||_infty with C_3 inherited from Sub-sprint A (the cross term
    [gamma^0 (x) ∂_t, a_s (x) a_t] VANISHES identically under the
    momentum-polynomial convention; the joint commutator reduces to the
    spatial commutator i [D_GV, a_s] (x) a_t).  Sub-sprint A established
    C_3^{joint} <= C_3^{SU(2)} <= 1.

K^+ positivity (the L3b foundation structural finding)
======================================================

Per `debug/l3b_first_move_memo.md` §3.2 and `geovac/krein_positive_state_space.py`,
the operator-multiplier-level Krein-positive restriction is TRIVIAL on the
compact-temporal substrate: every chirality-doubled scalar multiplier
M^{spat} = blkdiag(W, W) commutes with J = J_spatial (x) I_{N_t}, and
every momentum-polynomial temporal multiplier commutes with the identity
in the temporal slot.  Hence every generator preserves K^+, and so does
every B^{joint}(f) (Berezin maps weighted linear combinations of these).
The non-trivial K^+ structure lives at the STATE level, not the operator
level.  In this module we expose a `verify_krein_positive_preservation`
that confirms this finding for any B^{joint}(f) we construct, and the
positivity property restricts cleanly to the K^+ subspace via the
standard P_+ (.) P_+ compression.

Honest scope
============

(i) The mapping uses the FACTOR-WISE convolution-and-compression form:
    B^{joint}(f) = P^{joint} (K^{joint} * f) P^{joint}, equivalent to the
    spectral form on the central subalgebra via Plancherel.

(ii) Test functions are PURE TENSOR PRODUCT f = f_s(omega) * f_t(t) by
     default (the construction is bilinear and extends to general f by
     linearity).  We provide both pure-tensor and non-separable test
     functions in the panel.

(iii) The Fourier-mode q is signed (q in {-K_max, ..., +K_max}); the
      Fejer weight hat{K}(q) depends only on |q|.

(iv) The "natural Lipschitz norm" of f on S^3 x S^1_T uses the joint
     gradient ||nabla^{joint, L^2} f||_infty := sup sqrt(|nabla_x f|^2
     + |R_T^{-1} partial_t f|^2) per Sub-sprint A §1.4 (Pythagorean form).

API
===

  JointBerezinReconstruction(n_max, N_t, T)
      Class wrapping the joint Berezin map.

      .apply(f) -> matrix
      .verify_positivity(f)
      .verify_contractivity(f, f_infty)
      .verify_l3_compatibility(f, lipschitz_inf)
      .approximate_identity_residual(f)
      .verify_krein_positive_preservation(f)
      .factorize_pure_tensor(f_s, f_t) -> tensor product

  JointTestFunction
      Dataclass: spatial + temporal coefficient dictionaries.

  joint_panel(n_max, N_t)
      Standard verification panel of joint test functions.

Verification
============

Per CLAUDE.md §13.4a, every equation has a corresponding unit test in
`tests/test_joint_berezin_compact_temporal.py`.  The driver
`debug/l3b_2_sub_sprint_C_compute.py` regenerates the panel data.

References
==========

  Loutey, J. "GeoVac Paper 38: SU(2) propinquity convergence on the
  Camporesi-Higuchi spectral triple," Zenodo 2026; §L4 Berezin
  reconstruction.

  Sub-sprint A memo: debug/l3b_2_sub_sprint_A_lichnerowicz_memo.md
  Sub-sprint B memo: debug/l3b_2_sub_sprint_B_cb_norm_memo.md
  L3b foundation:    debug/l3b_first_move_memo.md
  Paper 38 L4 memo:  debug/r25_l4_proof_memo.md
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import sympy as sp

from geovac.berezin_reconstruction import (
    BerezinReconstruction,
    PlancherelSymbol,
)
from geovac.central_fejer_compact_temporal import (
    _u1_K_max,
    cb_norm_circle,
    plancherel_symbol_circle,
)
from geovac.central_fejer_su2 import normalization_constant
from geovac.krein_space_compact_temporal import (
    CompactTemporalKreinSpace,
    fourier_momentum_grid,
)
from geovac.lorentzian_dirac_compact import lorentzian_dirac_compact_matrix
from geovac.operator_system import TruncatedOperatorSystem
from geovac.operator_system_compact_temporal import (
    CompactTemporalTruncatedOperatorSystem,
    compact_temporal_multiplier_matrices,
)
from geovac.r25_l3_lipschitz_bound import TestFunction, make_test_function


# ---------------------------------------------------------------------------
# Joint Plancherel weights (Sub-sprint B factorization)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class JointPlancherelSymbol:
    """Joint Plancherel symbol on SU(2) x U(1).

    Factorizes by construction (Sub-sprint B Theorem 7.1):

        hat{K}^{joint}_{n_max, N_t}(N, q)
            = hat{K}^{SU(2)}_{n_max}(N) * hat{K}^{U(1)}_{N_t}(q)
            = (N / Z^{SU(2)}_{n_max}) * max(0, (N_t + 1 - 2|q|) / (N_t + 1)).

    The L^infty norm equals 2/(n_max + 1) (factor of 1 from the U(1)
    factor, attained at q = 0).
    """

    n_max: int
    N_t: int

    def __post_init__(self) -> None:  # pragma: no cover - sanity
        if self.n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {self.n_max}")
        if self.N_t < 1:
            raise ValueError(f"N_t must be >= 1, got {self.N_t}")

    @property
    def Z_su2(self) -> int:
        return normalization_constant(self.n_max)

    @property
    def K_max_u1(self) -> int:
        return _u1_K_max(self.N_t)

    def weight_su2(self, N: int) -> sp.Rational:
        """SU(2) factor weight hat{K}^{SU(2)}(N) = N/Z for N<=n_max, else 0."""
        if N < 1:
            raise ValueError(f"N must be >= 1, got {N}")
        if N > self.n_max:
            return sp.Rational(0)
        return sp.Rational(N, self.Z_su2)

    def weight_u1(self, q: int) -> sp.Rational:
        """U(1) factor weight hat{K}^{U(1)}(q) (standard Fejer)."""
        return plancherel_symbol_circle(self.N_t, q)

    def weight(self, N: int, q: int) -> sp.Rational:
        """Joint weight hat{K}^{joint}(N, q) = weight_su2(N) * weight_u1(q)."""
        return self.weight_su2(N) * self.weight_u1(q)

    def linfty_norm(self) -> sp.Rational:
        """||hat{K}^{joint}||_{l^infty} = 2 / (n_max + 1) (Sub-sprint B)."""
        return sp.Rational(2, self.n_max + 1) * cb_norm_circle(self.N_t)


# ---------------------------------------------------------------------------
# Joint Test Function (pure-tensor or non-separable)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class JointTestFunction:
    """A test function f on S^3 x S^1_T as a sum of pure-tensor terms.

    f(omega, t) = sum_alpha c_alpha * Y^{(3)}_{N_alpha, L_alpha, M_alpha}(omega)
                              * exp(i q_alpha t / R_T),

    where R_T = T / (2 pi) so the temporal mode q has natural period T.

    Fields
    ------
    name : str
    spatial_temporal_coeffs : dict ((N, L, M), q) -> complex
        Joint coefficients in the (spatial Peter-Weyl) x (temporal Fourier)
        basis.

    Empty dict = identically zero.
    """

    name: str
    coeffs: Tuple[Tuple[Tuple[Tuple[int, int, int], int], complex], ...]

    @property
    def coeff_dict(self) -> Dict[Tuple[Tuple[int, int, int], int], complex]:
        return dict(self.coeffs)

    def is_pure_tensor(self) -> bool:
        """Return True iff f = f_s * f_t (separable spatial x temporal)."""
        if not self.coeffs:
            return True
        spatial_part: Dict[Tuple[int, int, int], complex] = {}
        temporal_part: Dict[int, complex] = {}
        for ((N, L, M), q), c in self.coeff_dict.items():
            spatial_part[(N, L, M)] = spatial_part.get((N, L, M), 0) + c
            temporal_part[q] = temporal_part.get(q, 0) + c
        # Check rank-1 condition: c_{nm} = a_n * b_m for a outer product.
        # Equivalently: the (s, t) coeff matrix has rank 1.
        spatial_keys = sorted(spatial_part.keys())
        temporal_keys = sorted(temporal_part.keys())
        mat = np.zeros((len(spatial_keys), len(temporal_keys)), dtype=np.complex128)
        sp_idx = {k: i for i, k in enumerate(spatial_keys)}
        tp_idx = {k: i for i, k in enumerate(temporal_keys)}
        for ((N, L, M), q), c in self.coeff_dict.items():
            mat[sp_idx[(N, L, M)], tp_idx[q]] = c
        if mat.size == 0:
            return True
        rank = int(np.linalg.matrix_rank(mat, tol=1e-12))
        return rank <= 1

    def spatial_marginal_coeffs(self) -> Dict[Tuple[int, int, int], complex]:
        """For a pure tensor f = f_s * f_t with f_t(0) != 0, return f_s."""
        out: Dict[Tuple[int, int, int], complex] = {}
        for ((N, L, M), q), c in self.coeff_dict.items():
            if q == 0:
                out[(N, L, M)] = out.get((N, L, M), 0) + c
        return out

    def temporal_marginal_coeffs(self) -> Dict[int, complex]:
        """Sum over spatial labels at each q -- only meaningful for separable f."""
        out: Dict[int, complex] = {}
        for ((N, L, M), q), c in self.coeff_dict.items():
            out[q] = out.get(q, 0) + c
        return out


def make_joint_test_function(
    name: str,
    coeffs: Dict[Tuple[Tuple[int, int, int], int], complex],
) -> JointTestFunction:
    """Construct a JointTestFunction from a coefficient dict."""
    coeffs_tuple = tuple(sorted(coeffs.items()))
    return JointTestFunction(name=name, coeffs=coeffs_tuple)


def pure_tensor_function(
    name: str,
    spatial_coeffs: Dict[Tuple[int, int, int], complex],
    temporal_coeffs: Dict[int, complex],
) -> JointTestFunction:
    """Construct a pure-tensor f = (sum c_s Y) * (sum c_t e^{i q t/R_T})."""
    joint: Dict[Tuple[Tuple[int, int, int], int], complex] = {}
    for (N, L, M), c_s in spatial_coeffs.items():
        for q, c_t in temporal_coeffs.items():
            joint[((N, L, M), q)] = c_s * c_t
    return make_joint_test_function(name, joint)


# ---------------------------------------------------------------------------
# Joint Berezin reconstruction map
# ---------------------------------------------------------------------------


class JointBerezinReconstruction:
    """Joint Berezin map B^{joint}: C(S^3 x S^1_T) -> O^L (L3b-2 Sub-sprint C).

    For a JointTestFunction with coefficients c_{(N, L, M), q}, the map
    returns

        B^{joint}(f) = sum_{N <= n_max, |q| <= K_max}
                          hat{K}^{SU(2)}_{n_max}(N) *
                          hat{K}^{U(1)}_{N_t}(q) *
                          c_{(N, L, M), q} *
                          (M^{spat}_{N L M} (x) M^{temp}_q),

    where M^{temp}_q is the diagonal indicator-of-momentum matrix on
    C^{N_t} (entry 1 at momentum index k = q, 0 elsewhere).  For pure
    tensor input the construction factorizes as

        B^{joint}(f_s * f_t) = B^{SU(2)}(f_s) (x) B^{U(1)}(f_t).

    Parameters
    ----------
    n_max : int
    N_t : int
    T : float
    op_sys : CompactTemporalTruncatedOperatorSystem, optional

    Attributes
    ----------
    op_sys : the compact-temporal operator system
    n_max, N_t, T : truncation parameters
    plancherel : JointPlancherelSymbol
    spatial_berezin : BerezinReconstruction (the SU(2) factor)
    momentum_grid : ndarray of momentum indices
    """

    def __init__(
        self,
        n_max: int,
        N_t: int,
        T: float = 2.0 * np.pi,
        *,
        op_sys: Optional[CompactTemporalTruncatedOperatorSystem] = None,
    ) -> None:
        self.n_max = n_max
        self.N_t = N_t
        self.T = T
        if op_sys is None:
            self.op_sys = CompactTemporalTruncatedOperatorSystem(
                n_max=n_max, N_t=N_t, T=T
            )
        else:
            if op_sys.n_max != n_max or op_sys.N_t != N_t:
                raise ValueError(
                    f"op_sys ({op_sys.n_max}, {op_sys.N_t}) != "
                    f"({n_max}, {N_t})"
                )
            self.op_sys = op_sys

        self.plancherel = JointPlancherelSymbol(n_max=n_max, N_t=N_t)
        self.momentum_grid = fourier_momentum_grid(N_t)
        self._k_to_idx = {int(k): i for i, k in enumerate(self.momentum_grid)}

        # Spatial Berezin map (SU(2) factor) -- we use the scalar operator
        # system here so the Berezin weights are computed correctly; the
        # actual matrices we build are the chirality-doubled lifts via
        # op_sys._spat_matrices.
        self._scalar_op_sys = TruncatedOperatorSystem(n_max=n_max)
        self.spatial_berezin = BerezinReconstruction(
            n_max=n_max, op_sys=self._scalar_op_sys
        )

        # Cache the chirality-doubled spatial multipliers indexed by (N, L, M)
        self._spat_label_to_full = dict(
            zip(self.op_sys.spat_labels, self.op_sys._spat_matrices)
        )

        # Identity in temporal slot for marker-of-momentum-q diagonal
        self._I_t = np.eye(N_t, dtype=np.complex128)

    # -----------------------------------------------------------------
    # Plancherel weights
    # -----------------------------------------------------------------

    @property
    def dim_K(self) -> int:
        return self.op_sys.dim_K

    @property
    def dim_spatial(self) -> int:
        return self.op_sys.dim_spatial

    def all_weights(self) -> Dict[Tuple[int, int], sp.Rational]:
        """Joint weight grid {(N, q): hat{K}^{joint}(N, q)} on the kept range."""
        out: Dict[Tuple[int, int], sp.Rational] = {}
        K_max = self.plancherel.K_max_u1
        for N in range(1, self.n_max + 1):
            for q in range(-K_max, K_max + 1):
                out[(N, q)] = self.plancherel.weight(N, q)
        return out

    # -----------------------------------------------------------------
    # Temporal-mode-q indicator multiplier
    # -----------------------------------------------------------------

    def momentum_mode_matrix(self, q: int) -> np.ndarray:
        """Diagonal indicator matrix M^{temp}_q on C^{N_t}.

        Returns diag(delta_{k, q}) for k in the momentum grid; zero matrix
        if q is outside the kept range.
        """
        out = np.zeros((self.N_t, self.N_t), dtype=np.complex128)
        if q in self._k_to_idx:
            i = self._k_to_idx[q]
            out[i, i] = 1.0
        return out

    # -----------------------------------------------------------------
    # Apply the joint Berezin map
    # -----------------------------------------------------------------

    def apply(self, f: JointTestFunction) -> np.ndarray:
        """Compute B^{joint}(f) as a (dim_K, dim_K) complex matrix.

        Spectral form:

            B^{joint}(f) = sum_{(N, L, M, q): N <= n_max, |q| <= K_max}
                              hat{K}^{joint}(N, q) *
                              c_{(N L M), q} *
                              (M^{spat}_{N L M} (x) diag-q).
        """
        out = np.zeros((self.dim_K, self.dim_K), dtype=np.complex128)
        K_max = self.plancherel.K_max_u1
        for ((N, L, M), q), c in f.coeff_dict.items():
            if N > self.n_max:
                continue
            if abs(q) > K_max:
                continue
            if (N, L, M) not in self._spat_label_to_full:
                continue
            weight = self.plancherel.weight(N, q)
            w_float = complex(float(weight))
            if abs(w_float) < 1e-30:
                continue
            s_mat = self._spat_label_to_full[(N, L, M)]
            t_mat = self.momentum_mode_matrix(q)
            out += w_float * complex(c) * np.kron(s_mat, t_mat)
        return out

    def apply_unweighted(self, f: JointTestFunction) -> np.ndarray:
        """Reference: P^{joint} M_f P^{joint} (unweighted joint compression).

        Same construction but with weights = 1 instead of hat{K}^{joint}.
        Used in the L4 (c) approximate-identity residual.
        """
        out = np.zeros((self.dim_K, self.dim_K), dtype=np.complex128)
        K_max = self.plancherel.K_max_u1
        for ((N, L, M), q), c in f.coeff_dict.items():
            if N > self.n_max:
                continue
            if abs(q) > K_max:
                continue
            if (N, L, M) not in self._spat_label_to_full:
                continue
            s_mat = self._spat_label_to_full[(N, L, M)]
            t_mat = self.momentum_mode_matrix(q)
            out += complex(c) * np.kron(s_mat, t_mat)
        return out

    # -----------------------------------------------------------------
    # Pure-tensor factorization
    # -----------------------------------------------------------------

    def apply_pure_tensor(
        self,
        f_s: TestFunction,
        f_t: Dict[int, complex],
    ) -> np.ndarray:
        """For f = f_s * f_t separable, compute B^{joint}(f).

        Uses the factorization

            B^{joint}(f_s * f_t)
              = B^{SU(2)}(f_s)_chirality_doubled (x) B^{U(1)}(f_t)

        where B^{U(1)}(f_t) = sum_q hat{K}^{U(1)}(q) c_t(q) diag-q.

        Parameters
        ----------
        f_s : TestFunction (SU(2) spatial)
        f_t : dict q -> coefficient (temporal Fourier mode)
        """
        # Spatial: build full chirality-doubled spatial Berezin
        # NOTE: spatial_berezin is on the SCALAR Fock op-sys; we use it
        # only to compute the weighted scalar action, then lift to
        # chirality-doubled.  But the simplest path is to directly
        # apply the joint formula by constructing the equivalent
        # joint test function.
        joint_coeffs: Dict[Tuple[Tuple[int, int, int], int], complex] = {}
        for (N, L, M), c_s in f_s.coeff_dict.items():
            for q, c_t in f_t.items():
                joint_coeffs[((N, L, M), q)] = c_s * c_t
        joint_f = make_joint_test_function(
            f"pure_tensor_{f_s.name}", joint_coeffs
        )
        return self.apply(joint_f)

    def factor_check(
        self,
        f_s: TestFunction,
        f_t: Dict[int, complex],
        tol: float = 1e-12,
    ) -> Tuple[bool, float]:
        """Verify B^{joint}(f_s * f_t) = B^{spatial-lift}(f_s) (x) B^{U(1)}(f_t).

        The "spatial Berezin" here is the chirality-doubled-spinor-lift
        version, with weights hat{K}^{SU(2)}_{n_max}(N) applied to the
        spinor multipliers M^{spat}_{N L M} = blkdiag(W, W) from
        op_sys._spat_matrices.  This is the SAME content as Paper 38 §L4
        but represented on the spinor bundle instead of the scalar Fock
        basis.

        Returns (factored_ok, residual_F_norm).
        """
        joint_apply = self.apply_pure_tensor(f_s, f_t)

        # Spatial Berezin lift directly on the chirality-doubled spinor
        # multipliers
        spat_doubled = np.zeros(
            (self.dim_spatial, self.dim_spatial), dtype=np.complex128
        )
        for (N, L, M), c_s in f_s.coeff_dict.items():
            if N > self.n_max:
                continue
            if (N, L, M) not in self._spat_label_to_full:
                continue
            weight = float(self.plancherel.weight_su2(N))
            spat_doubled += weight * complex(c_s) * \
                self._spat_label_to_full[(N, L, M)]

        # Temporal Berezin: sum_q hat{K}(q) c_q diag-q
        temp = np.zeros((self.N_t, self.N_t), dtype=np.complex128)
        K_max = self.plancherel.K_max_u1
        for q, c_t in f_t.items():
            if abs(q) > K_max:
                continue
            weight = float(self.plancherel.weight_u1(q))
            temp += weight * complex(c_t) * self.momentum_mode_matrix(q)

        factored = np.kron(spat_doubled, temp)
        residual = float(np.linalg.norm(joint_apply - factored))
        return residual < tol, residual

    # -----------------------------------------------------------------
    # Property verifications
    # -----------------------------------------------------------------

    def verify_positivity(
        self, f: JointTestFunction, *, tol: float = 1e-9
    ) -> Tuple[bool, float]:
        """L4 (a): if f >= 0 pointwise, B^{joint}(f) is PSD.

        Returns (is_PSD, min_eigenvalue).  Caller is responsible for
        ensuring f >= 0; see joint_panel for canonical positive entries.
        """
        B = self.apply(f)
        B_herm = (B + B.conj().T) / 2
        eigvals = np.linalg.eigvalsh(B_herm)
        min_eig = float(np.min(eigvals.real)) if eigvals.size > 0 else 0.0
        return (min_eig >= -tol, min_eig)

    def operator_norm(self, f: JointTestFunction) -> float:
        """||B^{joint}(f)||_op (largest singular value)."""
        B = self.apply(f)
        if B.size == 0:
            return 0.0
        return float(np.linalg.norm(B, ord=2))

    def verify_contractivity(
        self,
        f: JointTestFunction,
        f_infty_norm: float,
        *,
        tol: float = 1e-9,
    ) -> Tuple[bool, float, float]:
        """L4 (b): ||B^{joint}(f)||_op <= ||f||_infty.

        Returns (is_contractive, op_norm, ratio = op_norm / ||f||_infty).
        """
        op_norm = self.operator_norm(f)
        if f_infty_norm < tol:
            return (op_norm <= tol, op_norm, 0.0)
        ratio = op_norm / f_infty_norm
        return (ratio <= 1.0 + tol, op_norm, ratio)

    def approximate_identity_residual(
        self, f: JointTestFunction
    ) -> Tuple[float, np.ndarray]:
        """L4 (c): ||B^{joint}(f) - P^{joint} M_f P^{joint}||_op.

        Returns (residual_op_norm, residual_matrix).
        """
        B = self.apply(f)
        P_M_P = self.apply_unweighted(f)
        delta = B - P_M_P
        if delta.size == 0:
            return 0.0, delta
        return float(np.linalg.norm(delta, ord=2)), delta

    def commutator_with_lorentzian_dirac(
        self, f: JointTestFunction
    ) -> np.ndarray:
        """Compute [D_L, B^{joint}(f)] for the L4 (d) compatibility check."""
        D_L = lorentzian_dirac_compact_matrix(self.op_sys.krein)
        B = self.apply(f)
        return D_L @ B - B @ D_L

    def verify_l3_compatibility(
        self,
        f: JointTestFunction,
        lipschitz_inf: float,
        *,
        tol: float = 1e-9,
    ) -> Tuple[bool, float, float]:
        """L4 (d): ||[D_L, B^{joint}(f)]||_op <= C_3 * ||nabla^{joint} f||_infty.

        Returns (is_compatible, commutator_op_norm, ratio).  C_3 = 1 from
        Sub-sprint A inherits to the joint setting; we report the ratio
        and pass if ratio <= 1 + tol.
        """
        comm = self.commutator_with_lorentzian_dirac(f)
        comm_norm = float(np.linalg.norm(comm, ord=2))
        if lipschitz_inf < tol:
            return (comm_norm <= tol, comm_norm, 0.0)
        ratio = comm_norm / lipschitz_inf
        return (ratio <= 1.0 + tol, comm_norm, ratio)

    # -----------------------------------------------------------------
    # K^+ positivity preservation (structural finding from L3a-1)
    # -----------------------------------------------------------------

    def verify_krein_positive_preservation(
        self, f: JointTestFunction, *, tol: float = 1e-10
    ) -> Tuple[bool, float]:
        """Verify B^{joint}(f) commutes with J (preserves K^+).

        Per L3a-1 + L3b foundation: chirality-doubled scalar multipliers
        commute with J = J_spatial (x) I_{N_t} structurally.  Their linear
        combinations (Berezin weights are scalar) inherit the same
        property.  This is the operator-multiplier-level finding that
        K^+-positivity is TRIVIAL at the operator level (the non-trivial
        program shifts to the STATE level via KreinPositiveStateSpace).

        Returns (preserves_K_plus, residual ||[J, B]||_F).
        """
        B = self.apply(f)
        J = self.op_sys.krein.J
        comm = J @ B - B @ J
        residual = float(np.linalg.norm(comm))
        return residual < tol, residual

    # -----------------------------------------------------------------
    # Riemannian limit (load-bearing check)
    # -----------------------------------------------------------------

    def reduce_to_paper38_at_N_t_1(
        self, f_s: TestFunction, *, tol: float = 1e-12
    ) -> Tuple[bool, dict]:
        """At N_t = 1 with f = f_s * 1 (constant temporal), B^{joint}
        equals the chirality-doubled spinor-lift Berezin map.

        Mechanism: the joint apply at q = 0 sums over (N, L, M) with
        spatial multipliers M^{spat}_{N L M} = blkdiag(W_{N L M}, W_{N L M})
        (chirality-doubled Weyl spinor lift, NOT the scalar Fock-basis
        matrix from the SU(2) Paper 38 Berezin).  The temporal factor is
        the trivial 1x1 identity I_1.

        The "Paper 38 reference" therefore is: the SU(2) Plancherel
        weights hat{K}^{SU(2)}(N) applied to the chirality-doubled
        Weyl spinor multipliers, summed over the support of f_s.

        This is the SAME content as Paper 38 §L4 (which works on the
        scalar Fock basis), but acting on a DIFFERENT representation
        (the spinor lift), with weights inherited identically.

        Constructs the reference manually from the op_sys spatial
        multipliers and the SU(2) Plancherel weights.
        """
        if self.N_t != 1:
            jb_1 = JointBerezinReconstruction(
                n_max=self.n_max, N_t=1, T=self.T
            )
            return jb_1.reduce_to_paper38_at_N_t_1(f_s, tol=tol)

        # At N_t = 1: joint apply with f_t = {0: 1} (constant temporal)
        joint_coeffs = {
            ((N, L, M), 0): c
            for (N, L, M), c in f_s.coeff_dict.items()
        }
        joint_f = make_joint_test_function(
            f"reduction_{f_s.name}", joint_coeffs
        )
        joint_apply = self.apply(joint_f)

        # Reference: weighted sum of chirality-doubled spatial multipliers,
        # tensored with I_1 (which is just a scalar multiplication).
        ref = np.zeros(
            (self.dim_spatial, self.dim_spatial), dtype=np.complex128
        )
        for (N, L, M), c in f_s.coeff_dict.items():
            if N > self.n_max:
                continue
            if (N, L, M) not in self._spat_label_to_full:
                continue
            weight = float(self.plancherel.weight_su2(N))
            spat_mat = self._spat_label_to_full[(N, L, M)]
            ref += weight * complex(c) * spat_mat
        # N_t = 1: temporal factor is the trivial 1x1 identity
        ref_full = np.kron(ref, self._I_t)

        residual = float(np.linalg.norm(joint_apply - ref_full))
        details = {
            "N_t": self.N_t,
            "n_max": self.n_max,
            "f_s_name": f_s.name,
            "residual_F_norm": residual,
            "ref_dim": ref_full.shape[0],
            "joint_dim": joint_apply.shape[0],
        }
        return residual < tol, details


# ---------------------------------------------------------------------------
# Joint test panel
# ---------------------------------------------------------------------------


def joint_constant_function() -> JointTestFunction:
    """The constant function 1 = Y^{(3)}_{1,0,0} * e^{i 0 t/R_T}.

    The unique non-negative pure-tensor function with both factors
    constant.  B^{joint}(constant) is a scalar multiple of the
    chirality-doubled identity, manifestly PSD.
    """
    return make_joint_test_function(
        "joint_constant_Y3_(1,0,0)_q0",
        {((1, 0, 0), 0): 1.0},
    )


def joint_axisymmetric_positive(
    n_max: int, N_t: int, eps_s: float = 0.01, eps_t: float = 0.05
) -> JointTestFunction:
    """A positive axisymmetric joint function f = f_s * f_t.

    f_s = Y_{1,0,0} + eps_s * Y_{2,0,0} (spatial positivity from Paper 38 L4)
    f_t = 1 + 2 * eps_t * cos(t / R_T) (temporal positivity from Fejer-on-S^1
                                        if eps_t < 1/2; here also a sum
                                        e^{i q} + e^{-i q} = 2 cos)

    The product is pointwise non-negative when both factors are; suitable
    for L4 (a) positivity test.
    """
    spatial = {(1, 0, 0): 1.0, (2, 0, 0): eps_s}
    temporal = {0: 1.0, 1: eps_t, -1: eps_t}
    return pure_tensor_function(
        f"joint_axisymmetric_positive_eps{eps_s}_eps{eps_t}",
        spatial,
        temporal,
    )


def joint_separable_single_mode(
    N: int, L: int, M: int, q: int
) -> JointTestFunction:
    """Pure-tensor with single spatial Y_{N,L,M} and single temporal mode q.

    Used for L4 (b)-(d) testing (NOT in general positive on S^3 x S^1_T
    since single Y for N >= 2 has nodes).
    """
    return make_joint_test_function(
        f"Y3_({N},{L},{M})_q{q}",
        {((N, L, M), q): 1.0},
    )


def joint_non_separable(
    N1: int, L1: int, M1: int, q1: int,
    N2: int, L2: int, M2: int, q2: int,
) -> JointTestFunction:
    """Sum of two pure-tensor terms with no shared structure (non-separable).

    f = Y_{N1,L1,M1} * e^{i q1 t/R_T} + Y_{N2,L2,M2} * e^{i q2 t/R_T},
    with (N1,L1,M1) != (N2,L2,M2) and q1 != q2 in general so the
    spatial-temporal coefficient matrix has rank 2.
    """
    return make_joint_test_function(
        f"Y3_({N1},{L1},{M1})_q{q1} + Y3_({N2},{L2},{M2})_q{q2}",
        {
            ((N1, L1, M1), q1): 1.0,
            ((N2, L2, M2), q2): 1.0,
        },
    )


def joint_panel(n_max: int, N_t: int) -> List[JointTestFunction]:
    """L4 joint verification panel.

    Returns a list of JointTestFunction covering:
      - constant (PSD-applicable)
      - axisymmetric positive perturbation (PSD-applicable)
      - single pure-tensor Y_{N,L,M} * e^{i q t} for low (N, L, M, q)
      - non-separable sum with rank-2 coefficient matrix

    Only entries with N <= n_max and |q| <= K_max contribute non-trivially.
    """
    panel: List[JointTestFunction] = []

    # Constant (PSD-applicable)
    panel.append(joint_constant_function())

    # Axisymmetric positive (PSD-applicable when eps_s, eps_t small enough)
    panel.append(joint_axisymmetric_positive(n_max, N_t))

    # Single pure-tensor modes (general, not necessarily positive)
    K_max = _u1_K_max(N_t)
    for N in range(2, min(n_max, 3) + 1):
        for L in range(min(N, 2)):
            for M in range(-min(L, 1), min(L, 1) + 1):
                for q in [0, 1, -1]:
                    if abs(q) > K_max:
                        continue
                    panel.append(joint_separable_single_mode(N, L, M, q))

    # Non-separable sum (rank-2 coefficient matrix)
    if n_max >= 2 and K_max >= 1:
        panel.append(
            joint_non_separable(2, 0, 0, 0, 2, 1, 0, 1)
        )
    if n_max >= 3 and K_max >= 1:
        panel.append(
            joint_non_separable(2, 0, 0, 1, 3, 0, 0, -1)
        )

    return panel


# ---------------------------------------------------------------------------
# Helpers for joint Lipschitz norm
# ---------------------------------------------------------------------------


def temporal_lipschitz_inf(
    f_t_coeffs: Dict[int, complex], T: float = 2.0 * np.pi
) -> float:
    """||partial_t f_t||_infty for f_t = sum c_q e^{i q t / R_T}.

    partial_t f_t = (1/R_T) sum c_q (i q) e^{i q t/R_T},
    R_T = T / (2 pi).  Bounded by (2 pi / T) sum |c_q| * |q|.
    """
    R_T = T / (2.0 * np.pi)
    s = 0.0
    for q, c in f_t_coeffs.items():
        s += float(abs(c)) * abs(q)
    return s / R_T


def joint_lipschitz_inf_pure_tensor(
    f_s_lip_inf: float,
    f_s_inf: float,
    f_t_lip_inf: float,
    f_t_inf: float,
    metric: str = "L2",
) -> float:
    """Joint Lipschitz norm for a pure-tensor f = f_s * f_t.

    Per Sub-sprint A §1.4:
      L^1 form: ||nabla_x f_s||_inf * ||f_t||_inf + ||f_s||_inf * ||partial_t f_t||_inf
      L^2 form: sup sqrt(|nabla_x f|^2 + |partial_t f|^2), bounded above by
                sqrt((||nabla_x f_s||_inf * ||f_t||_inf)^2 +
                     (||f_s||_inf * ||partial_t f_t||_inf)^2).

    For pure-tensor symbols, both forms are upper bounds on the joint
    Lipschitz norm; the L^2 form is sharper.
    """
    a = f_s_lip_inf * f_t_inf
    b = f_s_inf * f_t_lip_inf
    if metric == "L1":
        return a + b
    if metric == "L2":
        return float(np.sqrt(a * a + b * b))
    raise ValueError(f"metric must be 'L1' or 'L2', got {metric}")


def joint_lipschitz_inf_approx(
    f: JointTestFunction,
    f_s_inf_lookup: Optional[Dict[Tuple[int, int, int], float]] = None,
    f_s_lip_lookup: Optional[Dict[Tuple[int, int, int], float]] = None,
    T: float = 2.0 * np.pi,
    metric: str = "L2",
) -> float:
    """Compute the joint Lipschitz norm via the term-by-term bound.

    For f = sum c_{NLMq} Y_{NLM} * e^{i q t/R_T}, the triangle inequality
    gives

        ||nabla^{joint} f||_inf
            <= sum |c| * (||nabla Y||_inf * 1 + ||Y||_inf * |q|/R_T)   [L^1]
            <= sum |c| * sqrt(||nabla Y||_inf^2 + (||Y||_inf * q / R_T)^2)  [L^2 per-term]

    Note this is an UPPER BOUND; the true joint Lipschitz norm may be
    smaller (the triangle inequality is generally not tight).

    Args:
        f : JointTestFunction
        f_s_inf_lookup : optional {(N, L, M): ||Y||_inf} (defaults to 1)
        f_s_lip_lookup : optional {(N, L, M): ||nabla Y||_inf} (defaults to 1)
        T : circumference
        metric : 'L1' or 'L2'

    Returns:
        Upper bound on ||nabla^{joint} f||_inf.
    """
    R_T = T / (2.0 * np.pi)
    total = 0.0
    for ((N, L, M), q), c in f.coeff_dict.items():
        y_inf = 1.0 if f_s_inf_lookup is None else f_s_inf_lookup.get((N, L, M), 1.0)
        y_lip = 1.0 if f_s_lip_lookup is None else f_s_lip_lookup.get((N, L, M), 1.0)
        spatial_contrib = y_lip
        temporal_contrib = y_inf * abs(q) / R_T
        if metric == "L1":
            total += abs(c) * (spatial_contrib + temporal_contrib)
        elif metric == "L2":
            total += abs(c) * float(np.sqrt(
                spatial_contrib ** 2 + temporal_contrib ** 2
            ))
        else:
            raise ValueError(f"metric must be 'L1' or 'L2', got {metric}")
    return total
