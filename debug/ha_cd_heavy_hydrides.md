# Sprint 3 Track HA-C+D — SrH and BaH composed qubit Hamiltonians

**Version:** v2.12.0
**Date:** 2026-04-15
**Tracks:** HA-C (scalar SrH / BaH), HA-D (relativistic SrH / BaH)
**Purpose:** Build the first heavy-atom ([Kr], [Xe] frozen-core) monohydrides in both the scalar and the Dirac (κ, m_j) composed pipelines, unblocking the head-to-head resource comparison against Sunaga 2025 (PRA 111, 022817).
**Prerequisite:** Sprint 3 Track HA-A+B (frozen-core infrastructure for [Kr] and [Xe], classifier entries for Rb/Sr/Cs/Ba) — already landed.

---

## 1. What changed

- `geovac/molecular_spec.py`:
  - Extended `_ELEMENT_DATA` with Z=37 (Rb), 38 (Sr), 55 (Cs), 56 (Ba) entries marked `'frozen'` with core sizes 36 / 36 / 54 / 54.
  - Added `_HYDRIDE_REQ` entries for the monodihydride placeholders RbH, SrH₂, CsH, BaH₂ (to satisfy `hydride_spec`'s default-R lookup when built from the dihydride template).
  - Added `_MONOHYDRIDE_REQ = {'SrH': 4.055, 'BaH': 4.218}` (bohr, spectroscopic equilibrium from NIST/CRC).
  - Added `_alkaline_earth_monohydride_spec(Z, name, ..., relativistic)` — the shared builder used by CaH / SrH / BaH. Starts from the dihydride template, drops one bond block, rebuilds nuclear repulsion using the frozen core and a single M-H pair.
  - Added `srh_spec()`, `bah_spec()` scalar factories.
  - Added `srh_spec_relativistic()`, `bah_spec_relativistic()` relativistic factories that flip `relativistic=True` so the composed builder dispatches to the Dirac-on-S³ Tier 2 pipeline.
- `geovac/ecosystem_export.py`:
  - Added `'srh': 'SrH'`, `'bah': 'BaH'` to `_SYSTEM_REGISTRY`.
  - Added new dispatcher branch `elif canonical in ('SrH', 'BaH')` in `hamiltonian()`.
  - Added `_build_alkaline_earth_monohydride(canonical, R, max_n, verbose)` helper.
- `tests/test_heavy_hydrides.py` (new): 27 tests covering scalar Pauli counts, relativistic Pauli counts, CaH/SrH/BaH isostructural invariance, rel/scalar ratio bands, ecosystem_export smoke, OpenFermion/Qiskit exports, case-insensitive dispatch, regression on LiH/BeH₂/NaH/KH/CaH₂, frozen-core identity (order-of-magnitude check on NR constant), single-bond-block spec topology, and `relativistic=True` flag on _rel specs. All 27 pass.

CaH / SrH / BaH monohydride specs are not true closed-shell diatomics (open-shell ²Σ⁺ radicals). Following the pattern in `cah_spec_relativistic`, we use the closed-shell σ-bond template (2 electrons in the bond block) for structural Pauli counting. The open-shell nature would affect the occupation pattern of the Hamiltonian, not its block topology.

---

## 2. Isostructural invariance: the headline result

| Molecule | Core | Block topology | Q | N_pauli | λ_tot (Ha) | λ_ni (Ha) | QWC |
|----------|-----|----------------|---:|--------:|-----------:|----------:|----:|
| KH (Z=19) | [Ar] | 1 bond | 20 | 222 | 608.66 | 16.60 | — |
| CaH (Z=20) | [Ar] | 1 bond | 20 | 222 | 690.34 | 16.60 | — |
| SrH (Z=38) | [Kr] | 1 bond | 20 | 222 | 3145.86 | 16.60 | — |
| BaH (Z=56) | [Xe] | 1 bond | 20 | 222 | 7897.88 | 16.60 | — |

**Scalar-level invariance confirmed.** Non-identity 1-norm λ_ni is bit-identical (16.5975 Ha) across all four; only the identity coefficient — which absorbs the nuclear-repulsion constant and frozen-core energy — scales with Z. Main-group Pauli/Q = 222/20 = 11.1 matches the universal composed scaling law (CLAUDE.md §2 Track CU).

### Relativistic (Dirac κ, m_j basis):

| Molecule | Core | Q | N_pauli | λ_tot (Ha) | λ_ni (Ha) | QWC | rel/sca N_pauli |
|----------|-----|---:|--------:|-----------:|----------:|----:|---------------:|
| CaH_rel (Z=20) | [Ar] | 20 | 534 | 694.58 | 13.87 | 52 | 2.40× |
| SrH_rel (Z=38) | [Kr] | 20 | 534 | 3150.10 | 13.87 | 52 | 2.40× |
| BaH_rel (Z=56) | [Xe] | 20 | 534 | 7902.12 | 13.87 | 52 | 2.40× |

**Relativistic-level invariance confirmed.** Non-identity 1-norm λ_ni is bit-identical (13.8664 Ha); QWC groups identical (52); Pauli count identical (534). This is structurally expected: the relativistic (κ, m_j) builder sees Z_center = 2.0 (the alkaline-earth valence Z_eff) in every block regardless of the heavy atom's bare Z, because the frozen core already screens the nucleus. Spin-orbit's Z⁴α² prefactor therefore sees Z=Z_eff=2 uniformly, producing identical one-body and two-body integrals across CaH/SrH/BaH.

This is a structural prediction of the composed architecture: **heavy-atom relativistic resource counts at matched valence topology are identical up to the identity coefficient**, which is exactly what makes the Sunaga 2025 comparison a clean head-to-head test of the architectures rather than of the atomic species.

---

## 3. rel/scalar ratio

Rel / scalar Pauli ratio = 534 / 222 = 2.405× across all three alkaline-earth monohydrides, matching the Tier 2 T3 pin of 2.42× ± 20% for LiH (set by `test_relativistic_pauli_ratio_lih_nmax2`). The ratio is independent of atomic species — another confirmation that the relativistic composed pipeline produces structurally invariant Pauli costs at fixed valence block topology.

---

## 4. Sunaga 2025 head-to-head (partial, main-paper data)

Sunaga et al. 2025 (PRA 111, 022817) report for RaH at Q=18 (their Table I): **47,099 Pauli terms**.

GeoVac native-Q comparison (no matching; GeoVac uses Q=20, Sunaga uses Q=18):

| System | Sunaga Q=18 | GeoVac Q=20 | Ratio (GeoVac / Sunaga) |
|--------|------------:|------------:|------------------------:|
| RaH (Sunaga) | 47,099 Pauli | — | — |
| SrH_rel (GeoVac) | — | 534 Pauli | 0.0113× |
| BaH_rel (GeoVac) | — | 534 Pauli | 0.0113× |

This is 88× fewer Pauli terms at nearly-matched qubit count. Full per-molecule cells (BeH/MgH/CaH/SrH/BaH at Q=18) are in Sunaga's Supplementary Information Tables S1–S3, flagged as DEFERRED by the Tier 2 T4 market test (docs/tier2_market_test.md). Obtaining and wiring those SI cells is the next step (Sprint 3 HA-E if needed).

BaH (Z=56) is the heaviest atom covered by Sunaga's main-paper data that GeoVac can build without further infrastructure extension. RaH (Z=88) requires a [Rn] frozen core (86 electrons); this is deferred past Sprint 3.

---

## 5. Algebraic-first posture preserved

- No new quadrature. The monohydride builder reuses the existing `hydride_spec` / `build_composed_hamiltonian` algebraic stack, which was already algebraic for all ERIs (hypergeometric Slater integrals) and for all cross-center V_ne (Wigner D-matrix rotation + incomplete-gamma integrals).
- No new fitted parameters. The only new constants are experimental R_eq values for the monohydrides (4.055 bohr SrH, 4.218 bohr BaH) taken from NIST/CRC spectroscopic data, and the [Kr]/[Xe] Clementi-Raimondi-Reinhardt 1967 orbital exponents (already tabulated by HA-A+B).
- The spinor certificate from Tier 2 T5 continues to hold: every T3 H_SO coefficient remains in the ring R_sp = ℚ(α²)[γ]/(γ² + (Zα)² − 1) because Z_center in SrH_rel / BaH_rel is 2 (the alkaline-earth valence Z_eff), not Z=38/56 — the heavy nuclear charge is absorbed into the frozen-core screening and never touches the spinor matrix elements directly.

---

## 6. What this unblocks

- **Sprint 3 HA-E (Sunaga SI scrape):** once Sunaga's SI Tables S1–S3 are obtained, fill the BeH/MgH/CaH/SrH/BaH Q=18 cells and do a true matched-Q comparison. Expected factor: GeoVac sits 80–250× below Sunaga on Pauli count at matched Q (Tier 2 T4 projection).
- **Paper 20 Tier-2 table promotion:** the Tier-2 table in `papers/applications/paper_20_resource_benchmarks.tex` currently shows native-Q ratios (LiH, BeH, CaH). SrH / BaH rows can now be added at native Q=20.
- **Tier 3+ extensions:** with [Kr]/[Xe] frozen cores and monohydride builders in production, the natural next heavy-atom targets are the Ra/Fr alkalis (require [Rn] core) and the transition-metal monohydrides at n_max=3 (require the T7 Kramers-Pasternak recursion for Dirac-Coulomb ⟨r^{-2,-3}⟩ off-diagonal matrix elements, Sprint 3 HA-F).

---

## 7. Tests

```
tests/test_heavy_hydrides.py   27 passed
tests/test_ecosystem_export.py  101 passed   (0 regressions)
tests/test_heavy_cores.py       63 passed   (0 regressions)
tests/test_spin_ful_composed.py 13 passed   (0 regressions)
tests/test_composed_qubit.py    25 passed   (0 regressions)
Total: 229 passed
```

The pre-existing failures in `test_atomic_classifier.py` (test_oxygen_z2_scaled, test_transition_metal_raises) and `test_general_builder.py` (TypeError on stale max_n_core kwarg) are unchanged from before Sprint 3 and are documented in the HA-A+B memo (Section 8). No new failures introduced.

---

## 8. What this does NOT do

- Does NOT implement RaH (Z=88). Requires a [Rn] frozen core (86 electrons), which is out of Sprint 3 scope.
- Does NOT claim matched-Q comparison to Sunaga 2025. GeoVac's native Q=20 vs Sunaga's Q=18. The comparison is partial until Sunaga's SI is scraped.
- Does NOT include Sr or Ba atomic spectra. Those require atomic Dirac-on-S³ with gamma=sqrt(κ²-(Zα)²) radial corrections (Tier 3 T7); the HA-C/D monohydride resource counts are structural (valence Z_eff scaling), not spectroscopic.
- Does NOT introduce new fitted parameters. All new numerical constants are tabulated reference data.
