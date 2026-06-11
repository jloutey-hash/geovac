# Screening kernel upgrade memo — multi-zeta extension and Cs HFS retest

**Date:** 2026-05-09 (post-Sprint Cs-HFS-v2; companion to `debug/cs_hfs_v2_compute_memo.md`).
**Sprint:** Multi-observable focal-length decomposition program (CLAUDE.md §1.8) — kernel upgrade probe at heavy-atom regime (Z=55).
**Verdict:** **POSITIVE on the engineering** (multi-zeta machinery installed, backward compat preserved, 25 new tests + 208/208 regression all green); **CLEAN NEGATIVE on the hypothesis** that a heuristic two-zeta upgrade closes the Cs HFS framework-native gap. The Z=55 framework cliff is **deeper than "single-zeta vs multi-zeta"**: the Clementi-Raimondi 1967 single-zeta exponents themselves are non-faithful for heavy-atom outer shells, and a heuristic two-zeta split inherits the same problem. Cs HFS framework-native residual goes from −47% (single-zeta) to −90% (two-zeta heuristic) — wrong direction. The diagnostic is unambiguous: closing the heavy-atom kernel gap requires either a published full multi-zeta tabulation (Bunge-Barrientos-Bunge 1993 / Koga-Tatewaki-Thakkar 2000 for Z>54) or a self-consistent HF iteration in `geovac/neon_core.py` (multi-week scope). The multi-zeta machinery added in this sprint is the foundation for either path.

---

## 1. U1 — Tabulation source choice

**Choice: hybrid path.** For light cores (Ne, Z=10), use Bunge-Barrientos-Bunge 1993 Table I directly (BBB93 is the canonical published tabulation, covers Z=2–54, 5-zeta expansions per orbital). For the [Xe] core (Z=54, used at Cs Z=55 and Ba Z=56), CR74 / BBB93 stop at Z=54 and a heuristic two-zeta extrapolation from CR67 single-zeta is the only sprint-window-feasible approach. Koga-Tatewaki-Thakkar 1993 (Theor. Chim. Acta 86, 425) extends RHF to Z=55–92 but the tabulated coefficients are not as widely accessible as BBB93's.

**Rationale.** The full BBB93 Table I for Xe is ~200 (c_i, ζ_i) entries at 7-digit precision across 11 orbitals × 5–9 STOs each. Hand-tabulating those would consume the entire sprint window before any compute can run. The two-zeta heuristic is a **scoping probe**: if it shifts E_6s and the framework HFS in the right direction proportionally, the kernel-gap diagnosis is confirmed and the path forward is full BBB93 tabulation. If it shifts in the wrong direction (or by less than expected), the diagnosis sharpens to "the limit is deeper than single-vs-multi-zeta count" and a different engineering response is required.

**Specific values used.** For Ne (regression case), three orbitals (1s, 2s, 2p) each as 4–5 STOs from BBB93 Table I leading digits:

| Orbital | STO primitives (n_i, ζ_i) | Coefficients c_i |
|---------|---------------------------|-------------------|
| Ne 1s   | (1, 9.485), (1, 15.566), (2, 2.06), (2, 4.00), (2, 6.43)   | (0.393, 0.629, −0.0001, 0.0023, −0.0006) |
| Ne 2s   | same primitives                                              | (−0.108, −0.213, 0.699, 0.402, −0.079)   |
| Ne 2p   | (2, 1.45), (2, 2.38), (2, 4.485), (2, 9.135)                | (0.218, 0.533, 0.329, 0.019)             |

These reproduce neutral Ne to within 0.1% of total electron count under direct numerical integration.

For the Xe core (used at Cs Z=55), the two-zeta expansion is constructed from the CR67 single-zeta exponent ζ_CR per orbital with inner/outer ratios `(zeta_inner/zeta_CR, zeta_outer/zeta_CR, c_in, c_out)`:

```
(n=1, l=0): (1.20, 0.83, 0.65, 0.40)   # 1s
(n=2, l=0): (1.20, 0.83, 0.60, 0.45)   # 2s
(n=2, l=1): (1.20, 0.83, 0.55, 0.50)   # 2p
... (similar for n=3,4,5)
(n=5, l=0): (1.20, 0.83, 0.60, 0.45)   # 5s
(n=5, l=1): (1.20, 0.83, 0.55, 0.50)   # 5p
```

The coefficients `c_in, c_out` are renormalized at build time so that integral |R_nl|² r² dr = 1 on a fine geometric grid; the renormalization scale is typically within 5–15% of the input values.

**Honest scope:** these ratios are *derived from* BBB93 Kr (Z=36) tabulations applied to Xe assuming the inner/outer split structure transfers. After diagnostic iteration (see §3), they were tightened to be more symmetric (mean ratio ≈ 1) but the diagnostic outcome is unchanged: **the heuristic two-zeta is not the right object for Cs**.

---

## 2. U2 — FrozenCore extension

**Implementation.**

New module `geovac/multi_zeta_orbitals.py` (~340 lines):

- `STO` dataclass: normalized Slater-type orbital primitive `χ_i(r) = N(n,ζ) · r^{n-1} · exp(-ζr)` with `N = (2ζ)^{n+1/2} / √((2n)!)`.
- `MultiZetaOrbital` dataclass: `R_{nl}(r) = Σ c_i χ_i(r)` with validation (l < n, primitive n ≥ l+1, length match).
- `_build_ne_orbitals_neutral()`: returns the BBB93 Ne tabulation as a list of three MultiZetaOrbital objects.
- `_build_two_zeta_orbital(n, l, occ, ζ_cr)`: constructs a two-zeta orbital from a single CR67 ζ value using the lookup table `_TWO_ZETA_SPLITS`.
- `build_two_zeta_xe_orbitals_from_cr(zetas)`: builds the 11-orbital [Xe] core for given CR67 ζ tuple.
- `density_from_orbitals(orbs, r)`: returns Σ occ_i · |R_i|² · r² as the radial density weight.
- `core_electron_count(orbs)`: sanity check (sums occupancies).
- `warn_multi_zeta_unavailable(core_type)`: emits UserWarning for cores without tabulation.

`FrozenCore` extension in `geovac/neon_core.py`:

- New constructor kwarg `screening: str = 'single_zeta'` (defaults preserve legacy behavior).
- New class attribute `_MULTI_ZETA_AVAILABLE = {'Xe'}` listing which cores have a multi-zeta path; others fall back to single-zeta with a UserWarning.
- New method `_build_density_multi_zeta(r)` dispatches to the multi-zeta builder for the registered core.
- `solve()` checks `self.screening` and chooses the appropriate density builder.
- `_solve_screened_radial`, `_solve_screened_radial_log`, `screened_psi_origin_squared` all gain a `screening='single_zeta'` parameter that threads through to the FrozenCore constructor.

**Backward compatibility.**
- Default `screening='single_zeta'` produces bit-identical Z_eff values to the legacy path, verified by `TestFrozenCoreMultiZetaIntegration.test_frozen_core_bit_identical_default_vs_explicit_single`.
- All 149 pre-existing `test_neon_core.py` tests pass without modification.
- All 17 `test_hyperfine_a_constant.py` and 17 `test_nuclear_electronic.py` (Track-NI) tests pass.
- The Sprint 7b Ca/Sr/Ba screened SO splitting computation (`screened_so_splitting`) uses the default screening path and is unaffected.

**New test coverage** in `tests/test_multi_zeta_orbitals.py` (26 tests, 1 slow-skipped):

- `TestSTOPrimitives` (3): normalization, formula correctness, compactness scaling.
- `TestMultiZetaOrbital` (5): construction, validation paths (l<n, length match, n≥l+1), density contribution.
- `TestBBB93Neon` (3): orbital count, electron count, density normalization to 10±1%.
- `TestXeTwoZetaBuilder` (5): orbital count (11), electron count for Cs/Ba (54), density normalization to 0.5%, error path for wrong tuple length.
- `TestFrozenCoreMultiZetaIntegration` (8): default screening, backward-compat bit-identical, invalid-screening raises, multi_zeta works for Xe, fallback warning for unsupported cores, density normalization, Z_eff at origin/infinity, existing callers unaffected.
- `TestCsValenceWithMultiZeta` (1, slow): Cs 6s with multi-zeta gives finite, bound eigenvalue.

Combined regression: **208 passed, 3 slow-skipped, 0 failures** across `test_neon_core.py`, `test_multi_zeta_orbitals.py`, `test_hyperfine_a_constant.py`, `test_nuclear_electronic.py`.

---

## 3. U3 — Cs HFS re-compute with multi-zeta

**Compute driver:** `debug/calc_track_cs_hfs_v3.py` runs the v2 five-component Roothaan autopsy at `screening='single_zeta'` and `screening='multi_zeta'` side-by-side using the same physical-constants block, the same `_solve_screened_radial_log` solver, and the same Layer-2 conventions.

**Headline result:**

| Component                                  | Single-zeta (v2) | Multi-zeta (v3) | Δ |
|--------------------------------------------|------------------|-----------------|----|
| C1: BF strict (Dirac g)                    | 782.96 MHz       | 147.92 MHz      | −81% |
| C3: BF + Schwinger a_e                     | 783.87 MHz       | 148.09 MHz      | −81% |
| C4: BF + Schwinger + Casimir F_R           | 1219.12 MHz      | 230.32 MHz      | −81% |
| **Framework-native total**                 | **1219.12 MHz**  | **230.32 MHz**  | **−81%** |
| C5 Zemach + BW corrections                 | +14.51 MHz       | +2.74 MHz       |    |
| L Multi-loop QED (Layer-2)                 | +0.12 MHz        | +0.02 MHz       |    |
| **Framework + Layer-2 total**              | **1233.75 MHz**  | **233.09 MHz**  | **−81%** |
| Experimental A                             | 2298.16 MHz      | 2298.16 MHz     |    |
| **Residual (framework-native)**            | **−46.95%**      | **−89.98%**     | **WORSE** |
| **Residual (framework + Layer-2)**         | **−46.32%**      | **−89.86%**     | **WORSE** |
| E_6s [eV]                                  | −1.475           | −1.170          |    |
| E_6s relerr vs NIST                        | −62.13%          | −69.95%         | WORSE |
| ψ_6s(0)² [bohr⁻³]                          | 1.328            | 0.251           | −81% |

**The diagnostic is unambiguous: the heuristic two-zeta upgrade SHIFTS THE FRAMEWORK IN THE WRONG DIRECTION.** E_6s gets less bound, |ψ(0)|² gets smaller, and the residual on A_HFS doubles (from −47% to −90%).

**Mechanism analysis.** The Z_eff(r) profiles at intermediate r (where the Cs 5p shell density peaks):

| r [bohr] | Z_eff (single) | Z_eff (multi)  | Δ%      |
|----------|----------------|----------------|---------|
| 0.001    | 55.0000        | 55.0000        | 0%      |
| 0.05     | 47.93          | 47.75          | −0.4%   |
| 0.10     | 35.56          | 34.08          | −4%     |
| 0.30     | 12.42          | 9.36           | **−25%** |
| 0.50     | 7.17           | 3.36           | **−53%** |
| 1.00     | 1.08           | 1.02           | −6%     |
| > 2.0    | 1.00           | 1.00           | 0%      |

The two-zeta expansion under-screens at intermediate r (0.3–0.5 bohr), making the inner core of Cs effectively too compact and reducing the residual nuclear charge that the 6s orbital sees in the screening-transition region. This makes the 6s orbital experience a more diffuse screening profile and pulls its eigenvalue UP (less bound).

**Root cause.** The Clementi-Raimondi 1967 single-zeta exponents for Z=55 outer shells are physically non-faithful: ζ_CR(5p, Cs) = 12.31 gives Z_eff(5p) = 5·12.31 = 61.56, peaking at r ≈ 0.63 bohr — but the real Cs 5p orbital physically peaks at r ≈ 2 bohr with Z_eff(5p)_physical ≈ 3–5. The CR67 hydrogenic mapping `R_{nl}(r; Z_eff)` is a Slater-orbital fit constrained to reproduce the total energy of the atom, not the actual orbital extent. A two-zeta expansion *built from these CR67 exponents* inherits the problem and amplifies it.

**Implication.** The screening-kernel limit at Z=55 is **NOT just a matter of "single-zeta vs multi-zeta count"**. Real closure requires either:

(a) **Published full multi-zeta tabulation** (Bunge-Barrientos-Bunge 1993 Table I for Xe; or Koga-Tatewaki-Thakkar 1993/2000 for Z>54). Each orbital is 5–9 STOs with proper coefficients (some negative, encoding orthogonalization to inner shells, which is what produces the radial nodes). This is a multi-day data-entry sprint; the multi-zeta machinery in this commit is the foundation.

(b) **Self-consistent Hartree-Fock iteration** in `geovac/neon_core.py` with a basis of STO primitives and SCF convergence. Multi-week scope; would close not only the Cs gap but extend the framework to general atoms with arbitrary screening.

The two-zeta heuristic between (a) and (b) is **not viable** as a kernel upgrade for the heavy-atom precision regime.

---

## 4. U4 — Light-atom regression and Z-cliff verdict

**Light-atom regression check.** None of the existing precision-catalogue computations are affected by the multi-zeta extension because:

1. **H 21 cm hyperfine, Mu HFS, Mu 1S-2S, Ps 1S HFS, Ps 1S-2S** all use Z=1 (or equal-mass leptonic systems) where there is no FrozenCore.
2. **He 2³P, He 2¹S-2³S** use Z=2 with explicit graph-native CI; no FrozenCore.
3. **D 1S HFS** uses Z=1 nucleus + electron; no FrozenCore.
4. **μH Lamb / μH HFS** use Z=1 nucleus + lepton; no FrozenCore.
5. **NaH/MgH₂/HCl/H₂S/SiH₄/PH₃ chemistry** uses [Ne] FrozenCore at Z=11–18; the multi-zeta path is NOT registered for [Ne], so a `screening='multi_zeta'` request falls back to single-zeta with a UserWarning. The chemistry computations all use the default `screening='single_zeta'` and are bit-identical to legacy.
6. **CaH/SrH/BaH SO splittings** (Sprint 7b) use `screened_so_splitting` with the default screening path; bit-identical to legacy.

**Bit-identical regression verified by `test_frozen_core_bit_identical_default_vs_explicit_single`.** All 208 existing tests pass with the new code.

**Z-cliff verdict (the load-bearing finding):**

| Z range | Production status | Kernel finding | Action |
|---------|-------------------|----------------|--------|
| Z=1     | sub-100 ppm       | Hydrogenic exact | — |
| Z=2     | 0.014% (He fs)    | Z_eff fits FCI well | — |
| Z=11–18 | sub-percent (chemistry) | CR67 [Ne] adequate | — |
| Z=19–20 | sub-percent (chemistry) | CR67 [Ar] adequate | — |
| Z=21–30 | engineering frontier | CR67 [Ar] qualitative for d-block | (open) |
| Z=31–36 | engineering frontier | CR67 [Ar]3d10 qualitative | (open) |
| Z=37–38 | qualitative (SO splittings −70%) | CR67 [Kr] qualitative | (open) |
| Z=55–56 | **−47% (Cs HFS)**    | **CR67 [Xe] non-faithful for outer shells** | **Z-cliff confirmed; heuristic two-zeta does NOT close** |

The Z-cliff sits at Z ≈ 30 for SO splittings (Sprint 7b 60–70% errors for Ca/Sr/Ba) and at Z ≈ 55 for HFS (Cs 47% gap). **The cliff is structural to the Clementi-Raimondi 1967 hydrogenic single-zeta orbital fit**, not just a matter of zeta count. The cliff cannot be closed within the GeoVac framework using the existing screening machinery; it requires either a full RHF tabulation upgrade or self-consistent HF.

---

## 5. Pattern-finding: §1.8 directive value

This sprint executes the §1.8 multi-observable focal-length directive: it tested whether a kernel upgrade closes a previously identified framework gap (the −47% Cs HFS residual from Sprint Cs-HFS-v2). The diagnostic outcome is:

- **Class 2 (framework kernel approximation gap):** SHARPENED. The kernel limit is not single-vs-multi-zeta; it is the underlying CR67 orbital fit being non-faithful at high Z for outer shells. Sprint Cs-HFS-v2 already named "Clementi-Raimondi single-zeta is qualitative for the Cs valence regime"; this sprint refines that to "the entire CR67 hydrogenic orbital structure is non-faithful at Z=55 outer shells, not just the single-zeta count."

- **Class 1 (literature convention mismatch):** Out of reach until the kernel gap closes. The framework-native residual −47% (or now worse −90% under heuristic two-zeta) is two orders of magnitude larger than the Roberts-Ginges 2015 vs Porsev-Beloy-Derevianko 2009 atomic-structure spread (~0.5–1%).

- **Class 3 (focal-length decomposition cataloguing):** The five-component decomposition (BF strict / Schwinger / Casimir F_R / Zemach / multi-loop QED) is structurally complete and sound. The components are correctly attributed; only the kernel-derived |ψ(0)|² is the load-bearing wrong piece. The decomposition discipline successfully isolates which projection step contains the limiting error.

**The diagnostic-before-engineering memory file (`feedback_diagnostic_before_engineering.md`) is validated again.** When PI's intuition flagged "missing something" after Sprint Cs-HFS-v2's −47% residual, this sprint took ONE additional engineering step to test the obvious next hypothesis (multi-zeta upgrade) and produced a cleanly informative negative. Without this diagnostic step, the next sprint would have been a multi-week BBB93 tabulation effort under the unverified assumption that more zetas closes the gap. We now know: more zetas alone does not close the gap; the underlying orbital model needs replacement.

---

## 6. Files

**Modified:**
- `geovac/neon_core.py` (~70 lines added): new `screening` constructor kwarg on FrozenCore; new `_MULTI_ZETA_AVAILABLE` registry; new `_build_density_multi_zeta` dispatch; threaded through `_solve_screened_radial` (+1 kwarg), `_solve_screened_radial_log` (+1 kwarg), `screened_psi_origin_squared` (+1 kwarg). Backward compat: default `screening='single_zeta'` is bit-identical to legacy.

**Added:**
- `geovac/multi_zeta_orbitals.py` (~340 lines): new module with STO/MultiZetaOrbital primitives, BBB93 Ne tabulation, Xe two-zeta builder from CR67, density helpers.
- `tests/test_multi_zeta_orbitals.py` (~250 lines, 26 tests + 1 slow): coverage of STO normalization, MultiZetaOrbital construction/validation, BBB93 Ne reference, Xe builder, FrozenCore integration, backward-compat bit-identity.
- `debug/calc_track_cs_hfs_v3.py` (~280 lines): kernel-comparison driver that reruns the v2 five-component autopsy at both screening levels.
- `debug/data/cs_hfs_v3.json`: full v3 result set with both single and multi-zeta decompositions.
- `debug/data/screening_kernel_comparison.json`: side-by-side comparison + verdict text.
- `debug/screening_kernel_upgrade_memo.md`: this memo.

**Test summary:**
| Test file                                    | Pass | Slow-skip | Status |
|----------------------------------------------|------|-----------|--------|
| `tests/test_neon_core.py` (existing)         | 149  | 2         | green  |
| `tests/test_multi_zeta_orbitals.py` (new)    | 25   | 1         | green  |
| `tests/test_hyperfine_a_constant.py`         | 17   | 0         | green  |
| `tests/test_nuclear_electronic.py` (Track NI)| 17   | 0         | green  |
| **Total**                                    | **208** | **3** | **green** |

No production GeoVac code modified outside `neon_core.py`. Backward compat verified by direct bit-identity test, and indirectly by all existing tests remaining green.

---

## 7. Open follow-ups (for the next sprint window)

**Tooling-addressable, named:**

1. **Full BBB93 Table I tabulation for Xe** (multi-day data entry, then ~1 day for Cs HFS retest). Hand-tabulate ~200 (c_i, ζ_i) values from BBB93 for Xe Z=54, plug into the existing `_build_xe_orbitals_bbb93()` slot in `multi_zeta_orbitals.py` (function not yet written; placeholder named in U2 implementation). Predicted outcome at full BBB93: E_6s should land near NIST −3.89 eV at ~5%, and ψ_6s(0)² should rise correspondingly. If the framework HFS residual then closes from −47% to ~−10%, the Z-cliff is identified and closed for the Xe core; if it stays at ~−30%, the residual is in the relativistic-enhancement Casimir vs full Bohr-Weisskopf gap.

2. **Self-consistent HF iteration in geovac/neon_core.py** (multi-week scope; the principled framework path). A SCF loop on a Slater-type basis would close all heavy-atom screening kernels in one engineering pass. This is the multi-week Sprint W1c-residual continuation already flagged in CLAUDE.md.

3. **Full Bohr-Weisskopf relativistic enhancement** (1–2 weeks; Tier 3 / Sprint TR-class). Even with a perfect screening kernel, the Casimir leading-order F_R = 1.555 vs full Bohr-Weisskopf F_R ≈ 2.6 leaves a ~+30% headroom. This is INDEPENDENT of the screening upgrade.

**Diagnostic-instrument extension** (long-term, §1.8 directive value): Cs PNC sprint, blocked behind both (1)+(3) above.

**Out of scope for this sprint:** Ne, Ar multi-zeta tabulations (would extend `_MULTI_ZETA_AVAILABLE` to {'Ne', 'Ar', 'Xe'}; light-atom kernel is NOT a current bottleneck per §4 above).

---

## 8. Verdict

**Engineering closure: SUCCESS.** Multi-zeta machinery installed, backward-compatible, tested, and ready for the proper full-tabulation upgrade.

**Z-cliff resolution: STRUCTURAL GAP REMAINS, NEW DIAGNOSIS.** The Z=55 framework cliff is not closed by a heuristic two-zeta; the residual goes from −47% to −90% in the wrong direction. The diagnostic refines the gap to "CR67 hydrogenic orbital fit is non-faithful for outer shells at high Z, not just single vs two zeta count." Closing requires full BBB93/KTT tabulation OR self-consistent HF.

**Cs PNC diagnostic-instrument target: STILL BLOCKED.** Until the kernel upgrade closes the framework HFS residual to sub-percent (currently −47% with single-zeta, would need to reach < 1% to expose any literature convention mismatch in Cs PNC).

**§1.8 directive value: VALIDATED.** This sprint demonstrates the discipline of taking one additional engineering step to test the obvious hypothesis BEFORE committing to a multi-week effort. The negative is cleanly informative: more zetas is not the gap-closer.
