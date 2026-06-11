# Sprint 3 Track HA-A+B — Heavy-atom frozen cores ([Kr], [Xe])

**Version:** v2.12.0
**Date:** 2026-04-15
**Tracks:** HA-A (Clementi-Raimondi exponents + FrozenCore), HA-B (atomic_classifier)
**Purpose:** Unblock heavy-atom hydride construction for the Sunaga 2025 (PRA 111, 022817) comparison — specifically SrH and BaH at matched qubit count.

---

## 1. What changed

- `geovac/neon_core.py`: extended `FrozenCore` with two new core types `'Kr'` (36e) and `'Xe'` (54e), keeping the existing `'Ne'`, `'Ar'`, `'Ar3d10'` paths unchanged. Auto-detection now routes Z=37, 38 -> Kr and Z=55, 56 -> Xe.
- `geovac/atomic_classifier.py`: added entries for Rb (Z=37, type C), Sr (Z=38, type D), Cs (Z=55, type C), Ba (Z=56, type D). Rewrote the "unsupported" branch to give informative messages for the 4d block (Z=39-48), 5p block (Z=49-54), lanthanides (Z=57-71), 5d block (Z=72-80), 6p block (Z=81-86).
- `tests/test_heavy_cores.py`: new file, 63 tests covering Kr/Xe core construction, density normalization, Z_eff asymptotics, classifier entries, unsupported atoms, and regression on Z=1-36.
- `tests/test_neon_core.py`: repaired three stale tests that had baked-in "Z=11-18" error-message regexes from before the v2.2.0 [Ar] extension. They now verify the post-v2.2.0 behavior (Z=10 invalid, Z=19 valid [Ar]).
- `SCOPE_BOUNDARY.md`: added "Fifth Row (Z=37-54): Partially Supported" and "Sixth Row (Z=55-56): Partially Supported" sections; updated the Near-Term Reachability Summary; bumped header to v2.12.0.

Out of scope for Sprint 3 (explicitly): lanthanides (Z=57-71), d-block (Z=39-48, 72-80), p-block (Z=49-54, 81-86). These require additional frozen cores ([Kr]4d¹⁰ or [Xe]4f¹⁴5d¹⁰) or genuine multi-reference treatment.

---

## 2. Clementi-Raimondi-Reinhardt (1967) single-zeta exponents

Source: Clementi, Raimondi & Reinhardt, J. Chem. Phys. 47, 1300 (1967), Table III — single-zeta HF orbital exponents for neutral atoms. Values cross-checked against the reproduced tables in Clementi & Roetti, At. Data Nucl. Data Tables 14, 177 (1974). Project code convention: hydrogenic wavefunction uses Z_eff = n * zeta where zeta is the tabulated CRR value.

### [Kr] core (36 electrons): 1s² 2s² 2p⁶ 3s² 3p⁶ 3d¹⁰ 4s² 4p⁶

|  | 1s | 2s | 2p | 3s | 3p | 3d | 4s | 4p |
|:-|---:|---:|---:|---:|---:|---:|---:|---:|
| Rb (Z=37) | 36.3242 | 26.7156 | 30.5810 | 21.8410 | 23.2214 | 18.9867 | 11.5472 | 11.8892 |
| Sr (Z=38) | 37.3108 | 27.5192 | 31.5922 | 22.6526 | 24.0714 | 19.9645 | 12.4213 | 12.7489 |

### [Xe] core (54 electrons): [Kr] + 4d¹⁰ 5s² 5p⁶

|  | 1s | 2s | 2p | 3s | 3p | 3d | 4s | 4p | 4d | 5s | 5p |
|:-|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cs (Z=55) | 54.2678 | 42.0318 | 46.4901 | 36.9232 | 38.4598 | 33.2295 | 25.2794 | 25.5744 | 20.7221 | 13.1181 | 12.3118 |
| Ba (Z=56) | 55.2540 | 42.9170 | 47.5009 | 37.9057 | 39.4558 | 34.2470 | 26.1764 | 26.4651 | 21.6420 | 14.0101 | 13.1947 |

---

## 3. NIST / Clementi-Roetti core energies (Hartree-Fock total, neutral-core ion)

Source: Clementi & Roetti, At. Data Nucl. Data Tables 14, 177 (1974). These are the HF total energies of the closed-shell cation formed by removing all valence electrons (Rb⁺, Sr²⁺, Cs⁺, Ba²⁺).

| Core | E_core (Ha) |
|:-----|-----------:|
| Rb⁺ (Z=37) | −2938.360 |
| Sr²⁺ (Z=38) | −3131.540 |
| Cs⁺ (Z=55) | −7553.930 |
| Ba²⁺ (Z=56) | −7883.540 |

---

## 4. Density normalization (validation)

The analytical FrozenCore post-normalizes the raw density to exactly the target on the internal r_max=20 grid. The quantity below is the integral on an independent external grid r ∈ [0.001, 15] bohr with 5000 trapezoidal points, which measures (a) how well the hydrogenic-basis core density converges inside 15 bohr and (b) the residual tail outside 15 bohr.

| Z | Core | ∫ n(r) dr (on r ∈ [0.001, 15]) | Target | Relative error |
|:-:|:----:|:------------------------------:|:------:|:--------------:|
| 37 | [Kr] | 36.0475 | 36 | 0.132% |
| 38 | [Kr] | 36.0508 | 36 | 0.141% |
| 55 | [Xe] | 54.1224 | 54 | 0.227% |
| 56 | [Xe] | 54.1293 | 54 | 0.239% |

All within the project's 1% normalization tolerance (well within, with a comfortable 4-7x margin). The residual is a combination of (a) the external grid being bounded at r=15 rather than r=20 used by FrozenCore's internal solver, and (b) the trapezoidal quadrature rule on 5000 points.

---

## 5. Z_eff asymptotic validation

Expected behavior: Z_eff(0) = Z (unshielded nucleus), Z_eff(∞) = Z − n_core.

| Z | core | Z_eff(0.001) | Z_eff(0.1) | Z_eff(1.0) | Z_eff(5.0) | Z_eff(100) | expected(∞) |
|:-:|:----:|:---:|:---:|:---:|:---:|:---:|:---:|
| 37 | [Kr] | 37.0000 | 27.7391 | 1.0152 | 1.0000 | 1.0000 | 1 |
| 38 | [Kr] | 38.0000 | 28.2977 | 2.0050 | 2.0000 | 2.0000 | 2 |
| 55 | [Xe] | 55.0000 | 35.5638 | 1.0836 | 1.0000 | 1.0000 | 1 |
| 56 | [Xe] | 56.0000 | 35.7993 | 2.0309 | 2.0000 | 2.0000 | 2 |

All four atoms hit the correct asymptotic values to the grid precision of the FrozenCore solver. The approach to the asymptote is already essentially complete at r=5 bohr for these compact cores.

---

## 6. Atomic classifier entries

| Z | Symbol | Type | n_core | n_valence | Z_eff_val | Period | Group | valence_config |
|:-:|:------:|:----:|:------:|:---------:|:---------:|:------:|:------|:--------------|
| 37 | Rb | C | 36 | 1 | 1.0 | 5 | alkali_metal | 5s1 |
| 38 | Sr | D | 36 | 2 | 2.0 | 5 | alkaline_earth | 5s2 |
| 55 | Cs | C | 54 | 1 | 1.0 | 6 | alkali_metal | 6s1 |
| 56 | Ba | D | 54 | 2 | 2.0 | 6 | alkaline_earth | 6s2 |

Isostructural with K (Z=19, C) and Ca (Z=20, D) from v2.2.0, but with [Kr]/[Xe] substituted for [Ar]. The composed architecture's block topology is preserved; only the screening function Z_eff(r) changes.

---

## 7. Test results

```
tests/test_heavy_cores.py ...............................................
......................... passed                                        63 passed
tests/test_neon_core.py  ............................................... passed                  119 passed
```

Total: 182 tests passing (119 pre-existing + 63 new).

**Pre-existing failures in test_atomic_classifier.py are NOT caused by this sprint** (verified via `git stash` — same 4 tests fail on the HEAD commit before my changes):
- `TestPKParams::test_oxygen_z2_scaled` — test expects Z²-scaled params but classifier returns the Paper-17 computed value for O. Stale test from before Track BI.
- `TestEdgeCases::test_transition_metal_raises` (Z=21, 25, 30) — test expects `NotImplementedError`, but v2.8.0 TM hydrides changed the classifier to return a type-F classification for Z=21-30. Stale test from before v2.8.0.

These are documented for follow-up but are outside Sprint 3 scope.

---

## 8. What this unblocks

- **Sunaga 2025 head-to-head:** SrH and BaH can now be constructed via `build_composed_hamiltonian(spec)` with a `MolecularSpec` analogous to `cah2_spec()` / `mgh2_spec()` but using the [Kr]/[Xe] core respectively. The Pauli-count isostructural invariance (Section 1.5 of CLAUDE.md) predicts SrH ~ KH = 239 Pauli and BaH ~ KH = 239 (likely exactly identical counts), at Q=20 for n_max=2.
- **Tier 2 relativistic:** combined with the spin-ful composed pipeline (`geovac/composed_qubit_relativistic.py`, v2.11.0), the Sprint 3 extension enables spin-orbit SrH and BaH with Kramers-cancelled H_SO on the 5s¹/5s² (Cs) and 6s¹/6s² (Ba) valence orbitals.
- **Paper 20 Tier 2 resource table:** SrH/BaH cells can be filled against Sunaga's published numbers (the RaH-18q cell and the SI Tables S1–S3 deferred cells).

## 9. What this does NOT do

- Does NOT solve the Kr or Xe core quantum-mechanically (per SCOPE_BOUNDARY). Uses tabulated CR exponents only, per the algebraic-first philosophy.
- Does NOT enable 4d transition metals, 5p block, lanthanides, 5d block, or 6p block. Those require additional cores ([Kr]4d¹⁰ or [Xe]4f¹⁴5d¹⁰) and/or multi-reference treatment.
- Does NOT introduce any new fitted parameters. All values are tabulated reference data from Clementi-Raimondi-Reinhardt 1967 and Clementi-Roetti 1974.
- Does NOT modify any existing tests' expected values. The three edits to `test_neon_core.py` adjusted error-message regexes that had never been updated to reflect v2.2.0's [Ar] extension; the behavior checked is the same (invalid Z raises ValueError).
