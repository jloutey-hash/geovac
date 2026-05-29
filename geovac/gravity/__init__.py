"""GeoVac gravity subpackage.

Discrete-substrate gravity-arc infrastructure (sprint sequence G4-3
through G4-6). The continuum gravity-arc results (G1-G8) are
captured in Paper 28 ; the discrete-substrate
program is the multi-month commitment opened by G4-3.

Modules
-------
warped_dirac : G4-4a first move (2026-05-28). Discrete spinor bundle
    on the cigar D^2 x S^2 at constant warp, with the Dirac operator
    factorization at the operator level, the constant-warp tensor
    product, and the load-bearing falsifiers F1 (factorization) and
    F2 (chirality grading).
"""

from geovac.gravity.warped_dirac import (
    GAMMA_2D_1,
    GAMMA_2D_2,
    GAMMA_2D_5,
    I_2,
    DiscreteDirac2D,
    DiscreteDiskDirac,
    DiscreteDiskScalar,
    DiscreteWedgeDirac,
    S2DiracSpectrum,
    VariableWarpDirac,
    WarpedDiracConstant,
    verify_F1_factorization,
    verify_F2_chirality,
    verify_F3_continuum_recovery_rough,
    verify_F4_tip_regular,
    verify_F6_riemannian_limit,
    verify_F7_factorization_loss,
    verify_gamma_algebra_2d,
)

__all__ = [
    "GAMMA_2D_1",
    "GAMMA_2D_2",
    "GAMMA_2D_5",
    "I_2",
    "DiscreteDirac2D",
    "DiscreteDiskDirac",
    "DiscreteDiskScalar",
    "DiscreteWedgeDirac",
    "S2DiracSpectrum",
    "VariableWarpDirac",
    "WarpedDiracConstant",
    "verify_F1_factorization",
    "verify_F2_chirality",
    "verify_F3_continuum_recovery_rough",
    "verify_F4_tip_regular",
    "verify_F6_riemannian_limit",
    "verify_F7_factorization_loss",
    "verify_gamma_algebra_2d",
]
