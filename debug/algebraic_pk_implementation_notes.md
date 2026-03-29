# Algebraic PK Projector — Implementation Notes

**Date:** 2026-03-27
**Version:** v2.0.6-dev

---

## Summary

Replaced the Gaussian PK barrier `V_PK(r) = A·exp(-Br²)/r²` with the exact
Phillips-Kleinman rank-1 projector `V_PK = E_shift × |core⟩⟨core|` expressed
in the Level 4 channel basis.  The new mode coexists with existing Gaussian
modes (`channel_blind`, `l_dependent`, `r_dependent`) as `pk_mode='algebraic'`.

## Files Modified

| File | Changes |
|:-----|:--------|
| `geovac/core_screening.py` | Added `_extract_core_eigenvector()` to persist Level 3 ground-state angular eigenvector during `solve()`. Stores `core_eigenvector`, `core_l0_wavefunction`, `core_n_alpha`, `core_h_alpha`, `core_R_representative`. |
| `geovac/ab_initio_pk.py` | Added `algebraic_projector()` method returning a dict with the core l=0 wavefunction, grid parameters, and energy shift. |
| `geovac/composed_diatomic.py` | Added `pk_mode='algebraic'` support in constructor, `solve_core()`, and `_solve_valence_at_R()`. Added `LiH_algebraic_pk()` class method. Added `pk_projector` attribute. |
| `geovac/level4_multichannel.py` | Added `pk_projector` parameter threaded through `solve_level4_h2_multichannel()` → `compute_adiabatic_curve_mc()` → `solve_angular_multichannel()` → `build_angular_hamiltonian()`. Added algebraic PK injection block after existing Gaussian PK block. |

## Files Created

| File | Purpose |
|:-----|:--------|
| `tests/test_algebraic_pk.py` | 13 tests covering eigenvector persistence, projector shape, backward compatibility, smoke test, rank-1 verification, and channel coupling. |

## Design Decisions

1. **Atomic-limit approximation:** The core 1s² state is >98% in the l=0 channel
   (Paper 13, Section XII), so the Level 3 → Level 4 mapping is
   `|core⟩ ≈ δ_{(l1,l2),(0,0)} × u_core(α)`. No Clebsch-Gordan decomposition needed.

2. **Representative R:** The core eigenvector is extracted at the peak of the
   hyperradial wavefunction |F(R)|². This is the R value that contributes most
   to the density — a natural choice for a single-R representative.

3. **α-grid interpolation:** The Level 3 core uses n_alpha=200 (CoreScreening default)
   while the Level 4 solver uses n_alpha=100 (ComposedDiatomicSolver default).
   Cubic spline interpolation bridges the grids when they differ.

4. **R_e scaling:** The projector is injected as `H += R_e * E_shift * |c⟩⟨c|`,
   matching the charge-function units of the existing Hamiltonian (`R_e * V` terms).

5. **Two-electron factor:** The projector `|core⟩⟨core|` represents the
   two-electron core state (both electrons in the same angular configuration).
   The energy_shift = |E_core/N_core - E_val| × N_core already includes the
   N_core multiplier, so the projector is applied once (not per-electron).

## Test Results

All 13 new tests pass. Existing `test_composed_diatomic.py` regression tests
remain green (Gaussian PK behavior unchanged).

## Issues / Future Work

- **R_e dependence:** The current implementation uses a fixed core eigenvector
  (at the peak-R of the hyperradial wavefunction). A more rigorous approach would
  evaluate the core eigenvector at each R_e during the adiabatic sweep. This is
  straightforward but would require caching (similar to the existing rho-collapse
  cache infrastructure).

- **Beyond atomic limit:** For higher accuracy, the (l=0) → (0,0) mapping could
  be improved with Clebsch-Gordan corrections that populate (l,l) channels.
  The current implementation is correct to ~2% (the l>0 channel weight of the
  1s² core).
