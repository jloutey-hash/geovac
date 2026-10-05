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
    configure_geminal(gamma) -> float    # set the geminal range g (default 0.5); returns the old g
    geminal_gamma() -> float             # the current g

Geminal range.  g is ONE process-wide setting held in ``kernels.GAM`` (default 0.5) and set only via
``configure_geminal``; the geminal KIND ('exp' / 'linexp') is the per-call argument of energy(),
enumerators.Geminal(name) and pi_m16.GeminalM(name).  Every consumer reads g at call / construction
time; gVee's gamma tables (Kf, Kf2, _KS, Fm) are built at import for the current g and rebuilt on a
change, and the gamma-blind derived caches are cleared.  Calling it BEFORE importing the heavy
submodules avoids building the default-g tables first (and is the order the golden oracle uses).
Build Geminal objects after configuring.  ``<module>.GAM`` attributes remain readable (they return
the current g).  The former source override (source_override.install_kernels) is now a thin
wrapper over this; only the kernel GRID and the triangle cutoffs still need it.

The submodules (fourbody, energy, basis, kernels, hT, hVee, gVne, triangle, gVee, gT)
carry the validated primitives; assembly.py orchestrates them into the 2x2.

Many-determinant ("marriage") stack, promoted from debug/ in the LiH productionization
step 2c (2026-10-04; bit-identical against tests/data/lih_golden/golden.json; the debug/
originals are alias shims):

    debug/lih_marriage_phase0.py     -> ref_tensor      (74-det sigma reference, pair tensor T)
    debug/lih_marriage_phase1.py     -> bridge_gate     (gate G-leaf, 4-body bridge on real leaves)
    debug/lih_marriage_phase2.py     -> enumerators     (seven T-generalised enumerators, gate G1)
    debug/lih_marriage_phase3.py     -> sigma_marriage  (74-det sigma E_R12, mixing c, G2 VMC)
    debug/lih_pi_m16.py              -> pi_m16          (M=16 pi reference + m-aware enumerators)
    debug/lih_gvee_sparse.py         -> gvee_sparse     (O(nnz) g_Vee; the trunc-1e-3 pi marriage)
    debug/lih_exp_geminal_run.py     -> configure_geminal (+ source_override for grid / triangle cutoffs)
    debug/lih_pi_m16_vmc.py          -> pi_m16.build_wf_from_no (only that helper)
    debug/lih_r12_4body_integral.py  -> fourbody        (was already a verbatim copy here)
    (debug/lih_vmc.py                -> geovac.lih_vmc)
    assembly.energy's kernel closures + the phase2 copy -> kernels.build_kernel_gen / _m_gen

Hazards (import ``assembly.energy`` explicitly): after any submodule import of
``geovac.lih_r12ci.energy`` the package attribute ``energy`` is that SUBMODULE, so
``from geovac.lih_r12ci import energy`` is order-dependent -- use
``from geovac.lih_r12ci.assembly import energy``.  ``ref_tensor.grid_setup(exact)`` sets
``kernels.USE_EXACT_NEUMANN``; module RNG streams ``triangle.rng`` / ``pi_m16._RSVD_RNG`` are
consumed by the randomized eigh/SVD (bit-identity needs the same call order).
"""
from .assembly import energy, R12Result, VMC_TARGETS


def configure_geminal(gamma: float) -> float:
    """Set the geminal range gamma for the whole package; returns the previous gamma.
    Thin lazy wrapper over ``kernels.configure_geminal`` (kept lazy so that importing the package
    does not import kernels -- the grid source override needs kernels un-imported)."""
    from .kernels import configure_geminal as _configure
    return _configure(gamma)


def geminal_gamma() -> float:
    """The current geminal range (``kernels.GAM``)."""
    from . import kernels as _KG
    return _KG.GAM


__all__ = ["energy", "R12Result", "VMC_TARGETS", "configure_geminal", "geminal_gamma"]
