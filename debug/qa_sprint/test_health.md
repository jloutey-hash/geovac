# Test Suite Health Check

**Date:** 2026-03-15
**Status:** Complete

---

## 1. Overall Results

| Metric | Count |
|--------|:-----:|
| Total tests | ~579 |
| Passed | 569+ |
| Failed | 0 |
| Skipped (runtime) | varies |
| xfailed | 1 |
| Slow (>10 min) | 10 (prolate_heteronuclear_scf) |

The full suite passes with zero failures. One xfail is expected (Sturmian structural theorem). The 10 prolate heteronuclear SCF tests are extremely slow (>10 min total) and are a CI bottleneck.

---

## 2. Tests With No Assertions (ISSUE)

**`tests/test_lih_validation.py`** -- 2 tests are structurally broken:

- `test_lih_bridge_r_dependence()`
- `test_lih_conformal_factors()`

Both compute results and return a dict with `passed` boolean, but **never call `assert`**. Pytest counts them as PASSED regardless of actual result. Pytest warns: `PytestReturnNotNoneWarning: Test functions should return None`.

**Fix:** Add `assert result['passed']` before the return statement in both tests.

---

## 3. Conditionally Gated Test (WEAK)

**`tests/test_hylleraas.py::test_energy_below_atoms`** (line 272):

```python
if result['E_total'] < -1.0:
    assert result['D_e'] > 0
```

If the solver degrades and returns E >= -1.0, the test passes with zero assertions. The conditional is documented in the docstring ("may fail with very small grids") but means solver regression could go undetected.

**Fix:** Add `assert result['E_total'] < -1.0` as a separate assertion.

---

## 4. Loose Tolerances (>20%)

| File | Test | Tolerance | Notes |
|------|------|:---------:|:------|
| test_mo_fci.py:115 | test_lih_dissociation | 30% | xfail anyway |
| test_mo_fci.py:187 | test_j00_magnitude | 30% | J_00 vs atomic Li |
| test_prolate_grid_scf.py:94 | test_h2plus_eigenvalue | 40% | 2D FD inherently less accurate |
| test_molecular_sturmian.py:132 | beta spread | 30% | Spread/mean ratio |
| test_universal_constant_origin.py:234 | scale factor spread | 20% | Cross-Z relative spread |
| test_lih_fci.py:149 | binding energy diff | 20% | Energy separation |

These are all for methods with known limited accuracy. None are for core benchmarks.

---

## 5. No `assert True` Found

Zero instances of `assert True` in any test file. Good.

---

## 6. xfail and Skip Markers

### xfail (1 in test_ files)
- `tests/test_mo_fci.py::test_lih_dissociation` -- "Sturmian structural theorem: shared-p0 MO Sturmian FCI fails for heteronuclear molecules." Known theoretical limitation, documented.

### Runtime skips (~30 instances across suite)
- **Numba-dependent:** 6 in test_hylleraas.py (pass when Numba installed)
- **Complex scaling:** 5 in test_complex_scaling.py (no resonances / insufficient continuum)
- **Hyperspherical coupled:** 3 (correlation/width data insufficiency)
- **Prolate SCF:** ~8 (SCF/optimization failures -- guard against numerical instability)
- **MO FCI:** 2 (insufficient sigma/pi MOs)

**Concern:** The prolate SCF skips mean a broken solver could produce silent passes. Consider adding a minimum-pass-count assertion.

---

## 7. Solver Verification

### Key tests that exercise real solvers (not cached data)

| Test | Solver called | Tolerance |
|------|:-------------|:---------:|
| test_hyperspherical_he::test_energy_close_to_exact | `solve_helium()` | < 0.1% |
| test_complex_scaling::test_bound_state_real | `solve_ecs_single_channel()` | Im < 0.001 |
| test_complex_scaling::test_coupled_ground_state_real | `solve_ecs_coupled()` | Im < 0.001 |
| test_hylleraas::test_convergence_with_grid | Hylleraas solver | variational |
| test_berry_phase::test_plaquette_phase_analytical | `compute_plaquette_phase()` | 1e-12 |
| test_direct_ci::test_he_consistency | `DirectCISolver` | < 1e-8 Ha |
| test_h2_energy_decomposition::test_total_energy | Full CI pipeline | < 2.5% |

All critical benchmarks call real solvers. No mock data or cached results for core tests.

---

## 8. Recommended Fixes

1. **`test_lih_validation.py`**: Add assertions to 2 test functions
2. **`test_hylleraas.py::test_energy_below_atoms`**: Make the E < -1.0 check unconditional
3. **`test_prolate_heteronuclear_scf.py`**: Consider marking slow tests with `@pytest.mark.slow` and excluding from default runs
4. Consider adding a CI-friendly test subset marker for the ~570 fast tests (< 5 min total)
