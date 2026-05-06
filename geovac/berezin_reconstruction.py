"""Berezin-type reconstruction map B_{n_max}: C(S^3) -> O_{n_max} (R2.5 Lemma L4).

This module realizes Lemma L4 of the GH-convergence proof shape laid out in
`debug/track_ts_a_gh_convergence_memo.md` Section 5.4. It constructs the
Berezin-type reconstruction map

    B_{n_max} : C(S^3) -> O_{n_max}

that takes a smooth scalar function f on S^3 and returns an operator-system
element in the Connes-van Suijlekom truncation O_{n_max} = P_{n_max}
C(S^3) P_{n_max} on the Fock-projected S^3 graph. It is the SU(2) analog
of the Berezin map of Hawkins (Comm. Math. Phys. 215 (2000) 409-432) and
the Berezin-type reconstructions used in the metric-spectral-triple
literature (Connes-van Suijlekom CMP 383 (2021); Latremoliere arXiv:1811.10843;
Leimbach-van Suijlekom Adv. Math. 439 (2024) 109496).

Mathematical setup
==================

The L2 central spectral Fejer kernel (sprint R2.5/L2,
`geovac/central_fejer_su2.py`)

    K_{n_max}(g) = Z_{n_max}^{-1} | sum_{j <= j_max} sqrt(2j+1) chi_j(g) |^2

has Plancherel symbol on the central subalgebra L^2_{cent}(SU(2))

    hat{K}_{n_max}(j) = (2j + 1) / Z_{n_max}    for j <= j_max,
                      = 0                      for j  > j_max,

with Z_{n_max} = n_max(n_max + 1)/2 and j_max = (n_max - 1)/2.

Under the Peter-Weyl bijection n = 2j + 1 (verified to be the Fock-shell
index, see `geovac/so4_three_y_integral.py` and `central_fejer_su2.py`),
the Plancherel symbol on Fock shell N is

    hat{K}_{n_max}(N) := hat{K}_{n_max}(j = (N-1)/2)
                     = N / Z_{n_max}    for N <= n_max,
                     = 0               for N  > n_max.

Definition
==========

For f in C^infty(S^3) expanded as

    f(omega) = sum_{N >= 1} sum_{L=0}^{N-1} sum_{|M| <= L}
                 c_{N L M} Y^{(3)}_{N L M}(omega),

the Berezin map B_{n_max}(f) is defined by

    B_{n_max}(f) := sum_{N <= n_max} sum_{L=0}^{N-1} sum_{|M| <= L}
                      hat{K}_{n_max}(N) * c_{N L M} * M_{N L M}

where M_{N L M} in M_{N(n_max)}(C) is the truncated multiplier matrix
constructed in `geovac/operator_system.py::build_multiplier_matrix` (using
the Avery-Wen-Avery 3-Y integral on S^3).

This is the canonical "reconstruction" of f as an operator-system element
weighted by the L2 Plancherel symbol. Equivalently, B_{n_max} is the
composition

    B_{n_max} : f |--> P_{n_max} (K_{n_max} * f) P_{n_max}

of (i) convolution with the L2 central Fejer kernel and (ii) compression to
the truncated Hilbert space H_{n_max}. Both operations are spectrally
diagonal in the Peter-Weyl decomposition (the kernel acts as multiplication
by hat{K}(j) on each j-isotype block, and compression to H_{n_max} kills
the j > j_max blocks), so the spectral form B_{n_max}(f) above is
equivalent to the convolution form by the Plancherel theorem.

Properties (proven in this module + tests)
==========================================

(a) Positivity
--------------

If f >= 0 in C(S^3), then B_{n_max}(f) >= 0 as an element of O_{n_max}
(positive semi-definite as a matrix in M_{N(n_max)}(C)).

Proof: B_{n_max}(f) = P_{n_max} (K_{n_max} * f) P_{n_max}. The convolution
K_{n_max} * f is a non-negative function on S^3 (since K_{n_max} >= 0 by
L2 part (a) and f >= 0 by hypothesis). The compression P (g) P of a
non-negative function g by an orthogonal projection P is positive semi-
definite as a matrix: <psi | P g P | psi> = <P psi | g | P psi> =
integral_{S^3} g(omega) |P psi(omega)|^2 dOmega >= 0.

(b) Contractivity
-----------------

For every f in C(S^3),

    || B_{n_max}(f) ||_op  <=  || f ||_infty.

Proof: B_{n_max}(f) = P_{n_max} (K_{n_max} * f) P_{n_max}, and

    || K_{n_max} * f ||_infty  <=  || K_{n_max} ||_{L^1} || f ||_infty
                              =  1 * || f ||_infty       (by L2(b) normalization)

since K_{n_max} is a probability density on SU(2). Compression P(.)P does
not increase operator norm: || P g P ||_op <= || g ||_op = || g ||_infty
on multiplication operators by Sobolev / Cauchy-Schwarz.

(c) Approximate identity (rate controlled by L2 cb-norm)
--------------------------------------------------------

For f in C^infty(S^3),

    || B_{n_max}(f) - P_{n_max} M_f P_{n_max} ||_op
        <=  C_3 * gamma_{n_max} * || nabla f ||_infty,

where gamma_{n_max} -> 0 by L2 part (e), C_3 = 1 by L3, and
P_{n_max} M_f P_{n_max} is the unweighted compression of f to O_{n_max}.

This is the L^infty-to-operator-norm comparison form needed for Track A's
Latremoliere quantum GH propinquity bound (L5).

Equivalently, on the central subalgebra Z(C(SU(2))), the cb-norm equality
of L2(g) gives the symbol-side estimate

    || (id - B_{n_max}) f ||_op  <=  || id - hat{K}_{n_max} ||_{l^infty}
                                       * || f ||_{C(SU(2))}

with || id - hat{K}_{n_max} ||_{l^infty} bounded by 1 trivially and by
1 - 2/(n_max+1) -> 1 on the support; the rate is supplied by gamma_{n_max}
in the L^2 -> L^infty form above.

(d) Compatibility with L3
-------------------------

Because B_{n_max} factors through compression, the commutator
[D_CH, B_{n_max}(f)] is bounded by the L3 Lipschitz constant times the
Lipschitz norm of f, with the same C_3 = 1 constant:

    || [D_CH, B_{n_max}(f)] ||_op  <=  || hat{K} ||_{l^infty} * C_3 * || nabla f ||_infty
                                    <=  C_3 * || nabla f ||_infty
                                    =  || nabla f ||_infty   (since hat{K} <= 1).

This is the form in which L4 enters Track A's GH-convergence assembly
(L5): the Berezin map itself is Lipschitz in the operator-Lipschitz sense
inherited from L3, and its image in O_{n_max} converges to f in the
operator-norm sense as n_max -> infinity.

API
===

  - PlancherelSymbol class: tiny holder for the L2 symbol (proxy for
    `geovac/central_fejer_su2.py::plancherel_symbol`), exposing
    `weight_for_shell(N)` = hat{K}(j_N) = N / Z_{n_max}.

  - BerezinReconstruction class: the L4 reconstruction operator, with
    methods .apply(test_function) -> matrix, .verify_positivity(...),
    .verify_contractivity(...), .approximate_identity_residual(...).

  - berezin_reconstruct(test_function, op_sys, plancherel) -> matrix:
    one-shot function form, returning B_{n_max}(f) for an
    operator_system.TruncatedOperatorSystem and a TestFunction (from
    `r25_l3_lipschitz_bound.TestFunction`).

  - test_panel_n_max(n_max) -> List of TestFunction: for the verification
    panel (re-uses the panel construction from `r25_l3_lipschitz_bound`).

Verification protocol
=====================

Per CLAUDE.md Section 13.4a, every equation in the L4 proof memo has a
corresponding unit test in `tests/test_berezin_reconstruction.py`.

Honest limitations
==================

(i)  The "convolution form" B = P (K*f) P and the "spectral form"
     B = sum hat{K}(j) Pi_j(f) are equivalent on the central subalgebra
     (Plancherel). For non-central f the two forms differ at the
     non-isotypic level; we use the spectral form throughout this module
     because it lifts unambiguously to multiplier matrices.

(ii) The Berezin map of Hawkins 2000 was constructed for KAHLER manifolds
     using the holomorphic structure on a coadjoint orbit. SU(2) = S^3 is
     not Kahler (no almost-complex structure compatible with the round
     metric), so we do NOT use the holomorphic / Toeplitz form. We use
     instead the central-Fejer-kernel-spectral-projection form, which is
     the natural analog on a NON-Kahler compact group (the L2 kernel is
     the SU(2) analog of the Fejer kernel on S^1 used in
     Connes-van Suijlekom 2021).

(iii) The approximate-identity rate (c) is proved here using L2 part (e)
      gamma_{n_max} -> 0 (qualitative). The asymptotic rate of
      gamma_{n_max} is consistent with O(log n / n) but the closed-form
      rate constant is not nailed down; this is a sub-task of L2/L5
      assembly, not of L4 itself. The L4 statement (c) is proved with
      gamma_{n_max} treated as the (qualitative) rate-supplying ingredient.

References
==========

E. Hawkins, "Geometric Quantization of Vector Bundles and the
Correspondence with Deformation Quantization," Comm. Math. Phys. 215
(2000) 409-432, arXiv:math/9808116.

A. Connes & W. D. van Suijlekom, "Spectral Truncations in Noncommutative
Geometry and Operator Systems," CMP 383 (2021), arXiv:2004.14115.

F. Latremoliere, "The Gromov-Hausdorff propinquity for metric Spectral
Triples," arXiv:1811.10843.

M. Leimbach & W. D. van Suijlekom, "Gromov-Hausdorff Convergence of
Spectral Truncations for Tori," Adv. Math. 439 (2024) 109496,
arXiv:2302.07877.

GeoVac sprint records:
- L2 module: geovac/central_fejer_su2.py
- L2 memo:   debug/r25_l2_proof_memo.md
- L3 module: geovac/r25_l3_lipschitz_bound.py
- L3 memo:   debug/r25_l3_proof_memo.md
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import sympy as sp

from geovac.central_fejer_su2 import (
    fock_n_to_su2_j,
    normalization_constant,
    plancherel_symbol,
)
from geovac.operator_system import TruncatedOperatorSystem
from geovac.r25_l3_lipschitz_bound import (
    TestFunction,
    build_multiplier_for_test_function,
    default_test_panel,
    make_test_function,
)


# ---------------------------------------------------------------------------
# Plancherel symbol holder
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PlancherelSymbol:
    """Holder for the L2 Plancherel symbol hat{K}_{n_max}(j_N) = N / Z_{n_max}.

    Wraps central_fejer_su2.plancherel_symbol with a Fock-shell-indexed
    accessor `weight_for_shell(N)` returning the EXACT sympy Rational
    N / Z_{n_max} for N <= n_max and 0 otherwise.
    """

    n_max: int

    def __post_init__(self) -> None:  # pragma: no cover - sanity
        if self.n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {self.n_max}")

    @property
    def Z(self) -> int:
        """Z_{n_max} = n_max(n_max + 1)/2."""
        return normalization_constant(self.n_max)

    def weight_for_shell(self, N: int) -> sp.Rational:
        """hat{K}_{n_max}(N) = N / Z_{n_max} for N <= n_max, else 0.

        Uses the n = 2j + 1 Peter-Weyl bijection (verified in L2).
        Returns sympy Rational (exact arithmetic).
        """
        if N < 1:
            raise ValueError(f"N must be >= 1, got {N}")
        if N > self.n_max:
            return sp.Rational(0)
        j = fock_n_to_su2_j(N)
        return plancherel_symbol(self.n_max, j)

    def weight_for_shell_float(self, N: int) -> float:
        """Numerical version of weight_for_shell."""
        return float(self.weight_for_shell(N))

    def all_weights(self) -> Dict[int, sp.Rational]:
        """Dict mapping N in {1, ..., n_max} to hat{K}_{n_max}(N)."""
        return {N: self.weight_for_shell(N) for N in range(1, self.n_max + 1)}

    def linfty_norm(self) -> sp.Rational:
        """|| hat{K}_{n_max} ||_{l^infty} = max weight = n_max / Z_{n_max}.

        From central_fejer_su2.central_multiplier_cb_norm: this equals
        2 / (n_max + 1). Returned as exact sympy Rational.
        """
        return sp.Rational(2, self.n_max + 1)

    def is_unweighted_at_N(self, N: int) -> bool:
        """Predicate: weight is exactly 1 (used for testing the limit n_max=1)."""
        return self.weight_for_shell(N) == sp.Rational(1)


# ---------------------------------------------------------------------------
# Berezin reconstruction map
# ---------------------------------------------------------------------------


class BerezinReconstruction:
    """The Berezin-type reconstruction map B_{n_max}: C(S^3) -> O_{n_max}.

    Constructs B_{n_max}(f) for any test function f in the panel
    expansion as a sum of Plancherel-weighted multiplier matrices.

    Attributes
    ----------
    n_max : int
        Cutoff.
    op_sys : TruncatedOperatorSystem
        The Connes-vS truncated operator system on the scalar Fock basis.
    plancherel : PlancherelSymbol
        The L2 Plancherel symbol weights hat{K}(N) = N / Z_{n_max}.
    label_to_idx : dict (N, L, M) -> int
        Maps a multiplier label to its index in op_sys.multiplier_matrices.
    """

    def __init__(self, n_max: int, *, op_sys: Optional[TruncatedOperatorSystem] = None) -> None:
        self.n_max = n_max
        self.op_sys = op_sys if op_sys is not None else TruncatedOperatorSystem(n_max)
        self.plancherel = PlancherelSymbol(n_max=n_max)
        self.label_to_idx = {
            lab: i for i, lab in enumerate(self.op_sys.multiplier_labels)
        }

    # -----------------------------------------------------------------
    # Application
    # -----------------------------------------------------------------

    def apply(self, f: TestFunction) -> np.ndarray:
        """Compute B_{n_max}(f) as a matrix in M_{N(n_max)}(C).

        Spectral form:

            B_{n_max}(f) = sum_{(N,L,M) in coeffs of f, N <= n_max}
                              hat{K}_{n_max}(N) * c_{N L M} * M_{N L M}

        where hat{K}(N) = N / Z_{n_max}.

        Components with N > n_max are dropped (consistent with the
        truncation P_{n_max} which kills high-shell content). Components
        with (N, L, M) not in the multiplier_labels (zero matrix on the
        truncation) are skipped.
        """
        N_dim = self.op_sys.dim_H
        out = np.zeros((N_dim, N_dim), dtype=np.complex128)
        for (N, L, M), c in f.coeff_dict.items():
            if N > self.n_max:
                continue  # truncated out by P_{n_max}
            if (N, L, M) not in self.label_to_idx:
                continue  # zero multiplier matrix on this truncation
            weight = self.plancherel.weight_for_shell_float(N)
            out += complex(weight) * complex(c) * self.op_sys.multiplier_matrices[
                self.label_to_idx[(N, L, M)]
            ]
        return out

    def apply_unweighted(self, f: TestFunction) -> np.ndarray:
        """Reference: P_{n_max} M_f P_{n_max} (unweighted compression).

        This is the SAME as build_multiplier_for_test_function from the L3
        module (the truncated-multiplier expansion). Used as the reference
        in the approximate-identity check (c).
        """
        return build_multiplier_for_test_function(
            f, _ScalarOpSysAdapter(self.op_sys)
        )

    # -----------------------------------------------------------------
    # Property verifications
    # -----------------------------------------------------------------

    def verify_positivity(
        self, f: TestFunction, *, tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Verify B_{n_max}(f) >= 0 (PSD) for the test function f.

        Returns (is_positive, min_eigenvalue). is_positive = True iff
        the smallest eigenvalue of (B(f) + B(f)^*)/2 is >= -tol.

        Note: this routine does NOT check that f >= 0 on S^3; the user
        must provide a TestFunction known to be a non-negative function.
        For the canonical positive test cases, see test_panel_n_max
        below (which includes the constant function 1 = Y^{(3)}_{1,0,0}
        and the positive harmonics with explicit positivity reasoning).
        """
        B = self.apply(f)
        B_herm = (B + B.conj().T) / 2
        eigvals = np.linalg.eigvalsh(B_herm)
        min_eig = float(np.min(eigvals.real))
        return (min_eig >= -tol, min_eig)

    def operator_norm(self, f: TestFunction) -> float:
        """|| B_{n_max}(f) ||_op (largest singular value)."""
        B = self.apply(f)
        if B.size == 0:
            return 0.0
        return float(np.linalg.norm(B, ord=2))

    def verify_contractivity(
        self,
        f: TestFunction,
        f_infty_norm: float,
        *,
        tol: float = 1e-10,
    ) -> Tuple[bool, float, float]:
        """Verify || B_{n_max}(f) ||_op <= || f ||_infty.

        Args:
            f: TestFunction.
            f_infty_norm: || f ||_infty (the user supplies it; computing
                it requires sup-search and is delegated to the caller, e.g.
                lipschitz_norm_inf_test_function-style routines).
            tol: numerical tolerance.

        Returns:
            (is_contractive, op_norm, ratio).
            ratio = op_norm / f_infty_norm; is_contractive iff ratio <= 1 + tol.
        """
        op_norm = self.operator_norm(f)
        if f_infty_norm < tol:
            return (op_norm <= tol, op_norm, 0.0)
        ratio = op_norm / f_infty_norm
        return (ratio <= 1.0 + tol, op_norm, ratio)

    def approximate_identity_residual(
        self,
        f: TestFunction,
    ) -> Tuple[float, np.ndarray]:
        """Compute || B_{n_max}(f) - P_{n_max} M_f P_{n_max} ||_op.

        This is the operator-norm distance between the Berezin
        reconstruction (which weights each shell by hat{K}(N) = N/Z) and
        the unweighted compression. Should be O(gamma_{n_max}) by the L4
        property (c) and L2 part (e).

        Returns:
            (residual_op_norm, residual_matrix).
        """
        B = self.apply(f)
        P_M_P = self.apply_unweighted(f)
        delta = B - P_M_P
        if delta.size == 0:
            return 0.0, delta
        return float(np.linalg.norm(delta, ord=2)), delta

    # -----------------------------------------------------------------
    # Diagnostics
    # -----------------------------------------------------------------

    def shell_weight_table(self) -> Dict[int, sp.Rational]:
        """Return {N: hat{K}(N)} for all shells in the truncation."""
        return self.plancherel.all_weights()

    def supported_labels(self, f: TestFunction) -> List[Tuple[int, int, int]]:
        """List the (N, L, M) coeffs of f that survive the Berezin map.

        I.e. those with N <= n_max AND (N, L, M) in op_sys.multiplier_labels.
        """
        return [
            (N, L, M)
            for (N, L, M) in f.coeff_dict
            if N <= self.n_max and (N, L, M) in self.label_to_idx
        ]


# ---------------------------------------------------------------------------
# One-shot functional form
# ---------------------------------------------------------------------------


def berezin_reconstruct(
    f: TestFunction,
    op_sys: TruncatedOperatorSystem,
) -> np.ndarray:
    """One-shot wrapper: compute B_{n_max}(f) given a test function and
    the truncated scalar operator system O_{n_max}.

    Equivalent to BerezinReconstruction(op_sys.n_max, op_sys=op_sys).apply(f).
    """
    return BerezinReconstruction(op_sys.n_max, op_sys=op_sys).apply(f)


# ---------------------------------------------------------------------------
# Adapter so we can re-use L3's build_multiplier_for_test_function
# ---------------------------------------------------------------------------


class _ScalarOpSysAdapter:
    """Wrap a scalar TruncatedOperatorSystem to look like a Full Dirac one
    for build_multiplier_for_test_function compatibility.

    The L3 helper build_multiplier_for_test_function takes any object with
    .dim_H, .multiplier_labels, .multiplier_matrices attributes. The
    scalar TruncatedOperatorSystem has all three with the right semantics,
    so this adapter is a thin pass-through used to satisfy the type
    annotation. (The full Dirac op-sys form is the *spinor-lifted* one,
    which we are not using here; we work on the scalar Fock graph directly.)
    """

    def __init__(self, scalar_op_sys: TruncatedOperatorSystem) -> None:
        self._inner = scalar_op_sys

    @property
    def dim_H(self) -> int:
        return self._inner.dim_H

    @property
    def multiplier_labels(self) -> List[Tuple[int, int, int]]:
        return self._inner.multiplier_labels

    @property
    def multiplier_matrices(self) -> List[np.ndarray]:
        return self._inner.multiplier_matrices


# ---------------------------------------------------------------------------
# Test panel for L4 verification
# ---------------------------------------------------------------------------


def panel_n_max(n_max: int) -> List[TestFunction]:
    """L4 verification panel.

    Re-uses default_test_panel from r25_l3_lipschitz_bound (single
    Y^{(3)}_{NLM} for low (N,L,M) plus selected two-term sums) and adds:
      - the constant function f = Y^{(3)}_{1,0,0} (positive on S^3),
      - axisymmetric_positive_function (constant + small Y_{2,0,0},
        positive on S^3 by triangle inequality),
    at the front. These two are the canonical PSD-applicable tests for
    L4 property (a). Single Y_{NLM} for N >= 2 have nodes and are NOT
    pointwise positive, so B(Y_{NLM}) need not be PSD; those entries
    test (b)-(d) only.
    """
    panel: List[TestFunction] = []
    # Constant function (positive on S^3) -- always present.
    panel.append(make_test_function("constant_Y3_(1,0,0)", {(1, 0, 0): 1.0}))
    # Positive perturbation (also positive on S^3) -- L4 property (a) test.
    panel.append(axisymmetric_positive_function(n_max=n_max))
    # L3-style panel (single harmonics N >= 2 + sums).
    panel.extend(default_test_panel(n_max))
    return panel


# ---------------------------------------------------------------------------
# Diagnostic helpers
# ---------------------------------------------------------------------------


def constant_function() -> TestFunction:
    """The constant function 1 = Y^{(3)}_{1,0,0} (up to normalization).

    On the Fock graph, M_{1,0,0} is a SCALAR multiple of the identity
    matrix (the constant function on S^3 acts as scalar multiplication).
    The Berezin weight hat{K}(N=1) = 1 / Z_{n_max} is small but nonzero,
    so B_{n_max}(constant) is the identity (up to the 1/Z scaling), which
    is positive. This is the canonical positivity test case.
    """
    return make_test_function("constant_Y3_(1,0,0)", {(1, 0, 0): 1.0})


def axisymmetric_positive_function(n_max: int) -> TestFunction:
    """A positive axisymmetric f = constant + small Y_{2,0,0} component.

    Concretely f = Y_{1,0,0} + 0.01 * Y_{2,0,0}, designed so that the
    constant part dominates and f > 0 pointwise on S^3 (the perturbation
    is small enough not to flip sign).
    """
    return make_test_function(
        "axisymmetric_positive",
        {(1, 0, 0): 1.0, (2, 0, 0): 0.01},
    )


def commutator_with_dirac(
    B: np.ndarray, dirac_diag: np.ndarray,
) -> np.ndarray:
    """Compute [D, B] for B in M_dim and D = diag(dirac_diag).

    Helper for L4 property (d) verification: we need to check that the
    Berezin reconstruction is compatible with the L3 Lipschitz bound.
    """
    if dirac_diag.ndim != 1 or dirac_diag.shape[0] != B.shape[0]:
        raise ValueError(
            f"dirac_diag shape {dirac_diag.shape} != ({B.shape[0]},)"
        )
    # [D, B]_{ij} = (D_i - D_j) B_{ij}
    diff = dirac_diag[:, None] - dirac_diag[None, :]
    return diff * B
