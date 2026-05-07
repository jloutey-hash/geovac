"""Chirality grading γ_GV on the GeoVac full-Dirac truncated operator system.

Sprint G3 (Electroweak chirality co-location on S^3): this module is
**G3-A**, the GeoVac-side construction. It provides the Z_2 grading
γ_GV that is required for the inner-fluctuation Higgs mechanism on the
almost-commutative extension A_GV ⊗ (C ⊕ H) (Sprint H1, Paper 32 §VIII.C).
The companion sprint G3-B audits the AC-factor side γ_F. Sprint G3-C
then tests the identification (γ_GV ⊗ I_F) ≟ (I_GV ⊗ γ_F) on the
tensor-product spectral triple.

Mathematical setup
==================

On the full-Dirac sector basis (geovac/full_dirac_operator_system.py)
the Hilbert space H_GV(n_max) decomposes by chirality as

    H_GV = H_+ ⊕ H_-,        dim H_± = spinor_dim(n_max),

with the explicit basis ordering [chi=+1 block, chi=-1 block]. The
truthful Camporesi-Higuchi Dirac in this basis is

    D_truthful = diag(chi · (n_fock + 1/2))
               = sigma_z ⊗ D_+    (block form),

where D_+ = diag(n_fock + 1/2) is the magnitude on the Weyl basis. The
scalar-function multipliers act identically on both chirality blocks,

    M_full = I_2 ⊗ M_Weyl       (block form).

Choice of γ_GV
--------------

The NCG chirality grading must satisfy three axioms (Connes 1995,
van Suijlekom 2015 Ch. 3):

  (i)  γ² = I.
  (ii) {γ, D} = 0.
  (iii) [γ, a] = 0 for every a in the algebra.

In the [chi=+1, chi=-1] block ordering, the unique candidate (up to
overall sign) that satisfies all three with the truthful CH Dirac is

    γ_GV = sigma_x ⊗ I_{dim_weyl},

i.e. the operator that SWAPS the two chirality blocks pointwise on the
Weyl basis. Verification:

  (i)  (sigma_x ⊗ I)^2 = sigma_x^2 ⊗ I = I.
  (ii) {sigma_x ⊗ I, sigma_z ⊗ D_+}
       = (sigma_x sigma_z + sigma_z sigma_x) ⊗ D_+
       = 0 ⊗ D_+ = 0    (Pauli matrices anticommute).
  (iii) [sigma_x ⊗ I, I_2 ⊗ M_Weyl] = sigma_x ⊗ M_Weyl − sigma_x ⊗ M_Weyl = 0.

In particular γ M γ = M holds *exactly* for every scalar multiplier
(γ acts trivially on the operator system), so γ preserves O in the
strongest possible sense.

Comparison with the diagonal Z_2 candidate
-------------------------------------------

The naive "diagonal in chirality" choice γ_diag = sigma_z ⊗ I (i.e.
+I on the chi=+1 block, -I on the chi=-1 block) satisfies (i) and (iii)
but FAILS (ii):

    {sigma_z ⊗ I, sigma_z ⊗ D_+} = 2 sigma_z^2 ⊗ D_+ = 2 I ⊗ D_+ ≠ 0.

So γ_diag is NOT the chirality grading; it is a separate Z_2 sign
(the "D-eigenvalue sign") that *commutes* with D rather than
anticommuting. We provide γ_diag for diagnostic use only.

Behaviour on the offdiag CH Dirac
---------------------------------

For the Sprint-R3.5 offdiag CH Dirac
(camporesi_higuchi_offdiag_dirac_matrix), the within-chirality E1-style
ladder couplings are *equal* on the two blocks (offdiag_alpha is the
same constant on each chirality), so the within-chirality offdiag part
has the form [[X, 0], [0, X]] = I ⊗ X. This part *commutes* with
γ_GV = sigma_x ⊗ I, contributing 2 I ⊗ X to {γ_GV, D_offdiag} rather
than 0. Hence

    {γ_GV, D_offdiag} ≠ 0.

This is consistent with the WH1 round-3 finding that the offdiag CH
Dirac is an SDP-bounding device for the Connes distance, NOT a
spectral-action-foundational object (Sprint TS-C, debug/real_structure
_finite_nmax_memo.md). The chirality grading and the offdiag
construction are *different* perturbations of the truthful CH Dirac
and are not designed to be mutually consistent.

For the G3-C tensor-product test we therefore use γ_GV with the
TRUTHFUL CH Dirac, where {γ, D} = 0 holds exactly.

References
==========

A. Connes, "Noncommutative geometry and reality," J. Math. Phys. 36
(1995) 6194-6231.

W. D. van Suijlekom, "Noncommutative Geometry and Particle Physics,"
Springer 2015. Ch. 3 (Z_2 grading on a spectral triple).

R. Camporesi & A. Higuchi, "On the eigenfunctions of the Dirac operator
on spheres and real hyperbolic spaces," J. Geom. Phys. 20 (1996) 1-18.

GeoVac sprint records:
- Paper 32 §VIII.B (electroweak co-location target)
- Sprint H1 memo: debug/h1_ac_extension_memo.md
- R3.5 memo: debug/wh1_r35_full_dirac_memo.md
- Real structure J: geovac/real_structure.py
- Full Dirac operator system: geovac/full_dirac_operator_system.py
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    camporesi_higuchi_offdiag_dirac_matrix,
    full_dirac_basis,
    full_dirac_dim,
)
from geovac.spinor_operator_system import spinor_dim


# ---------------------------------------------------------------------------
# ChiralityGrading dataclass
# ---------------------------------------------------------------------------


@dataclass
class ChiralityGrading:
    """Z_2 chirality grading γ_GV on the GeoVac full-Dirac sector at n_max.

    Stored as a Hermitian unitary matrix with γ² = I.

    Attributes
    ----------
    n_max : int
        Fock cutoff.
    basis : list of FullDiracLabel
        Full-Dirac basis from `full_dirac_basis(n_max)`.
    matrix : complex ndarray of shape (dim_H, dim_H)
        γ_GV in the same basis ordering as
        `FullDiracTruncatedOperatorSystem(n_max).basis`.
    convention : str
        "sigma_x" (default, the NCG chirality grading; γMγ = M, {γ,D}=0)
        or "sigma_z" (diagnostic, the D-eigenvalue sign; commutes with D).

    Conventions
    -----------
    Basis ordering follows `full_dirac_basis`: chi = +1 block first
    (length spinor_dim(n_max)), chi = -1 block second (same length).
    Pauli matrices act on the chirality index with this ordering.

    γ_GV = sigma_x ⊗ I  (default, "sigma_x")
        Swaps the two chirality blocks. Satisfies γ²=I, {γ,D_truthful}=0,
        γ M γ = M for every scalar multiplier.

    γ_GV = sigma_z ⊗ I  ("sigma_z", diagnostic)
        +I on chi=+1, -I on chi=-1. Satisfies γ²=I and γMγ = M but
        COMMUTES with D_truthful instead of anticommuting. Provided
        only for negative-control tests.
    """

    n_max: int
    basis: Sequence[FullDiracLabel]
    matrix: np.ndarray
    convention: str

    @property
    def dim(self) -> int:
        return self.matrix.shape[0]

    def gamma_matrix(self) -> np.ndarray:
        """Return γ_GV as a complex ndarray (Hermitian unitary, γ² = I)."""
        return self.matrix

    # -- Axiom verifications ------------------------------------------------

    def verify_grading_squared(self) -> float:
        """Return max|γ² − I|. Should be ~0 for any Z_2 grading."""
        I = np.eye(self.dim, dtype=np.complex128)
        return float(np.max(np.abs(self.matrix @ self.matrix - I)))

    def verify_anticommutes_dirac(
        self,
        dirac_mode: str = "truthful",
        *,
        D: Optional[np.ndarray] = None,
        offdiag_kwargs: Optional[dict] = None,
    ) -> float:
        """Return max|{γ, D}| = max|γD + Dγ|.

        Parameters
        ----------
        dirac_mode : str
            "truthful": camporesi_higuchi_full_dirac_matrix (block-diagonal
                in chirality, eigenvalue chi*(n_fock + 1/2)). Expected
                anticommutator: ~0 for sigma_x convention, ≠ 0 for sigma_z.
            "offdiag": camporesi_higuchi_offdiag_dirac_matrix. Expected
                anticommutator: NONZERO for both conventions because the
                within-chirality offdiag part is symmetric across blocks
                (commutes with sigma_x).
            "custom": use the provided `D` argument directly.
        D : ndarray, optional
            Required when dirac_mode="custom".
        offdiag_kwargs : dict, optional
            Forwarded to `camporesi_higuchi_offdiag_dirac_matrix` when
            dirac_mode="offdiag".
        """
        if dirac_mode == "truthful":
            D_mat = camporesi_higuchi_full_dirac_matrix(self.basis)
        elif dirac_mode == "offdiag":
            kwargs = offdiag_kwargs or {}
            D_mat = camporesi_higuchi_offdiag_dirac_matrix(self.basis, **kwargs)
        elif dirac_mode == "custom":
            if D is None:
                raise ValueError("dirac_mode='custom' requires D argument")
            D_mat = np.asarray(D, dtype=np.complex128)
        else:
            raise ValueError(f"unknown dirac_mode {dirac_mode!r}")
        anticom = self.matrix @ D_mat + D_mat @ self.matrix
        return float(np.max(np.abs(anticom)))

    def verify_commutes_with_multiplier(
        self, M: np.ndarray
    ) -> float:
        """Return max|[γ, M]| = max|γM − Mγ| for a single multiplier M."""
        if M.shape != (self.dim, self.dim):
            raise ValueError(
                f"M shape {M.shape} != ({self.dim}, {self.dim})"
            )
        return float(np.max(np.abs(self.matrix @ M - M @ self.matrix)))

    def verify_preserves_operator_system(
        self,
        op_sys: FullDiracTruncatedOperatorSystem,
        *,
        tol: float = 1e-10,
    ) -> Tuple[bool, float, List[Tuple[int, float]]]:
        """For each multiplier M ∈ O, check γMγ ∈ O.

        Returns
        -------
        all_in_O : bool
        worst_residual : float
            Max residual over all generators.
        failures : list of (i, residual) for any generator that fails.
        """
        failures: List[Tuple[int, float]] = []
        worst = 0.0
        for i, M in enumerate(op_sys.multiplier_matrices):
            transformed = self.matrix @ M @ self.matrix
            in_O, residual = op_sys.contains(transformed, tol=tol)
            if residual > worst:
                worst = residual
            if not in_O:
                failures.append((i, residual))
        return (len(failures) == 0, worst, failures)

    def verify_unitary(self) -> float:
        """Return max|γγ^† − I|. Should be ~0 for any unitary."""
        I = np.eye(self.dim, dtype=np.complex128)
        return float(np.max(np.abs(self.matrix @ self.matrix.conj().T - I)))

    def verify_hermitian(self) -> float:
        """Return max|γ − γ^†|. Should be ~0 for a Hermitian Z_2 grading."""
        return float(np.max(np.abs(self.matrix - self.matrix.conj().T)))


# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------


def build_gamma_GV(
    n_max: int, *, convention: str = "sigma_x"
) -> ChiralityGrading:
    """Build the chirality grading γ_GV on the full-Dirac sector at n_max.

    Parameters
    ----------
    n_max : int
        Fock cutoff (>= 1).
    convention : str
        "sigma_x": natural NCG chirality grading. Anticommutes with the
        truthful CH Dirac, commutes with all scalar multipliers, and
        γ M γ = M exactly. **This is the production choice for G3-C.**
        "sigma_z": diagonal D-eigenvalue sign. Commutes with D rather
        than anticommuting; provided as a diagnostic / negative control.

    Returns
    -------
    ChiralityGrading

    Notes
    -----
    The basis ordering is [chi = +1 block, chi = -1 block] from
    `full_dirac_basis(n_max)`. Pauli matrices below act on the chirality
    label in this ordering; the inner factor is the identity on the
    Weyl-sector subspace of dim spinor_dim(n_max).
    """
    if n_max < 1:
        raise ValueError(f"n_max must be >= 1, got {n_max}")
    if convention not in ("sigma_x", "sigma_z"):
        raise ValueError(
            f"convention must be 'sigma_x' or 'sigma_z', got {convention!r}"
        )

    basis = full_dirac_basis(n_max)
    dim_full = len(basis)
    dim_weyl = spinor_dim(n_max)
    assert 2 * dim_weyl == dim_full

    I_w = np.eye(dim_weyl, dtype=np.complex128)
    if convention == "sigma_x":
        sigma = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    else:  # sigma_z
        sigma = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
    gamma = np.kron(sigma, I_w)

    return ChiralityGrading(
        n_max=n_max, basis=basis, matrix=gamma, convention=convention,
    )


# ---------------------------------------------------------------------------
# High-level audit driver
# ---------------------------------------------------------------------------


def audit_chirality_grading(
    n_max: int,
    *,
    convention: str = "sigma_x",
    verbose: bool = False,
) -> dict:
    """Run the full G3-A axiom audit at cutoff n_max.

    Parameters
    ----------
    n_max : int
    convention : str
        "sigma_x" (default, NCG chirality) or "sigma_z" (diagnostic).
    verbose : bool

    Returns
    -------
    dict with structured residuals for all axioms across both Dirac
    modes (truthful, offdiag).
    """
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    gamma = build_gamma_GV(n_max, convention=convention)

    # Axiom (i): γ² = I
    sq_residual = gamma.verify_grading_squared()
    # Hermiticity / unitarity sanity
    herm_residual = gamma.verify_hermitian()
    unit_residual = gamma.verify_unitary()
    # Axiom (ii) for two Dirac modes
    anticom_truthful = gamma.verify_anticommutes_dirac("truthful")
    anticom_offdiag = gamma.verify_anticommutes_dirac("offdiag")
    # Axiom (iii) preservation of O
    pres_ok, pres_worst, pres_failures = gamma.verify_preserves_operator_system(op_sys)
    # Per-generator commutator [γ, M]
    comm_residuals = [
        gamma.verify_commutes_with_multiplier(M)
        for M in op_sys.multiplier_matrices
    ]
    worst_comm = max(comm_residuals) if comm_residuals else 0.0

    result = {
        "n_max": n_max,
        "convention": convention,
        "dim_H": op_sys.dim_H,
        "dim_O": op_sys.dim,
        "n_generators": len(op_sys.multiplier_matrices),
        "gamma_squared_minus_I": sq_residual,
        "hermitian_residual": herm_residual,
        "unitary_residual": unit_residual,
        "anticommutator_truthful": anticom_truthful,
        "anticommutator_offdiag": anticom_offdiag,
        "preserves_O": pres_ok,
        "preserves_O_worst_residual": pres_worst,
        "preserves_O_failures": len(pres_failures),
        "max_commutator_with_M": worst_comm,
    }

    if verbose:
        print(f"audit_chirality_grading(n_max={n_max}, convention={convention!r})")
        print(f"  dim_H = {op_sys.dim_H}, dim_O = {op_sys.dim}, "
              f"|generators| = {len(op_sys.multiplier_matrices)}")
        print(f"  γ² − I       max residual: {sq_residual:.3e}")
        print(f"  γ Hermitian  max residual: {herm_residual:.3e}")
        print(f"  γ unitary    max residual: {unit_residual:.3e}")
        print(f"  {{γ, D_truthful}} max:     {anticom_truthful:.3e}")
        print(f"  {{γ, D_offdiag}}  max:     {anticom_offdiag:.3e}  "
              f"(expected nonzero — within-chirality offdiag is symmetric)")
        print(f"  γ preserves O: {pres_ok}  worst residual: {pres_worst:.3e}  "
              f"failures: {len(pres_failures)}")
        print(f"  max [γ, M]:   {worst_comm:.3e}")

    return result
