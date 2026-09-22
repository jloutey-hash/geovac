"""Two-center four-electron explicit-r12 CI on a minimal ionic LiH reference.

Production home of the LiH R12-CI proof-of-concept engine migrated from the sprint
tree (debug/lih_r12ci_*). It assembles the fully analytic, resolution-of-identity-free
{Phi0, (F-Fbar)Phi0} 2x2 for a James-Coolidge geminal on the ionic single-zeta reference
Phi0 = |1s_A^2 1s_B^2|, and diagonalizes to E_R12. Backing for Paper 12 Sec. "Explicit
correlation" (the two-center 4e RI-free reduction: 4-body bridge + 3-body triangle).

This is a PROOF OF CONCEPT (ionic single-zeta reference, not a spectroscopic LiH); the
deliverable is the RI-free reduction, not the correlation energy recovered. See
CHANGELOG v5.15.10-.14 and debug/lih_r12_build_plan.md for the full chronicle.

Public API:
    energy(geminal='exp') -> R12Result   # 'exp': f=e^{-g r} (-7.942); 'linexp': f=r e^{-g r} (-7.917)

The submodules (fourbody, energy, basis, kernels, hT, hVee, gVne, triangle, gVee, gT)
carry the validated primitives; assembly.py orchestrates them into the 2x2.
"""
from .assembly import energy, R12Result, VMC_TARGETS

__all__ = ["energy", "R12Result", "VMC_TARGETS"]
