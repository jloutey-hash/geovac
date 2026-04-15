# Dirac-on-S³ Tier 2 Sprint — Verdict (Track T6)

**Date:** 2026-04-15
**Status:** Tier 2 complete.
**Input tracks:** T0 (spinor ERI density), T1 (spinor matrix elements), T2 (Breit-Pauli SO), T3 (spin-ful composed pipeline), T4 (market test), T5 (π-free spinor certificate + Paper 18 subtier).
**Method:** Pure synthesis. No new computation. All numbers quoted exactly from T0–T5 memos.

---

## 1. Sprint outcome (one paragraph)

The Tier 2 Dirac-on-S³ sprint produced three deliverables, all positive. **(a) Infrastructure:** `geovac/dirac_matrix_elements.py` (T1, 117 tests), `geovac/spin_orbit.py` (T2, 22 tests), and `geovac/composed_qubit_relativistic.py` (T3, 13 tests + 164 pre-existing regression tests preserved bit-exactly) constitute a full algebraic-first pipeline from (κ, m_j) quantum numbers to JW qubit operators for three relativistic molecules (LiH, BeH, CaH). Every matrix element is closed form in exact sympy arithmetic or machine-precision float; no numerical quadrature enters the spinor path. **(b) Paper 22 extension:** T0 computed the spinor analog d_spinor(l_max) of Paper 22's potential-independent ERI density in both the Paper-22-native pair-diagonal convention and the full-Gaunt convention, finding d_spinor ≤ d_scalar at every l_max with ratio 1/4 (l_max=0, pure spin-dilution) → 0.92 (l_max=5, asymptotically equal). The sparsity-exponent result of Paper 22 extends verbatim; the prefactor is fixed by jj-coupling. **(c) Paper 18 completion:** T5 landed `verify_spinor_pi_free` and a drop-in for Paper 18 §IV adding the "spinor-intrinsic" subtier covering α² (Breit-Pauli coupling) and γ = √(1−(Zα)²) (Martínez-y-Romero radial factor, reserved). The coefficient ring is R_sp := ℚ(α²)[γ]/(γ²+(Zα)²−1); every T3 coefficient lies in R_sp under the certifier, with zero π, ζ, log, or E₁ contamination.

## 2. Headline result table

### 2.1 Three relativistic composed molecules (T3)

| n_max | Molecule | Q | N_Pauli (scalar) | N_Pauli (rel) | rel/scalar | λ_ni (Ha, rel) | QWC (rel) |
|:-----:|:---------|:-:|:----------------:|:-------------:|:----------:|:--------------:|:---------:|
| 1 | LiH | 6  |     9   |     9   | 1.00× |  10.15 |   1 |
| 1 | BeH | 6  |     9   |     9   | 1.00× |  63.44 |   1 |
| 1 | CaH | 4  |     6   |     6   | 1.00× |   2.03 |   1 |
| 2 | LiH | 30 |   333   |   805   | 2.42× |  35.90 |  55 |
| 2 | BeH | 30 |   333   |   805   | 2.42× | 141.32 |  52 |
| 2 | CaH | 20 |   222   |   534   | 2.41× |  13.87 |  52 |
| 3 | LiH | 84 |  7 878  |  46 434 | 5.89× | 126.75 | 6 571 |
| 3 | BeH | 84 |  7 878  |  46 434 | 5.89× | 297.03 | 6 571 |
| 3 | CaH | 56 |  5 252  |  30 940 | 5.89× |  65.81 | 6 571 |

Rel/scalar Pauli ratio is isostructural across the three molecules (1.00×/2.42×/5.89× at n_max=1/2/3). 1-norm rel vs scalar stays flat or decreases at n_max=3 (QPE-favorable; see §4 below).

### 2.2 Sunaga head-to-head (T4)

Sunaga et al. 2025 (PRA 111, 022817) publishes one calibrated Pauli-count cell: **RaH-18q = 47 099 Pauli strings**. Per-molecule comparisons for BeH/MgH/CaH/SrH/BaH at 18q require Supplemental Material Tables S1–S3 which are **deferred** (PDF not accessible to this sprint).

| Molecule | Q (GeoVac rel) | N_Pauli (GeoVac rel) | Sunaga RaH-18q | Ratio |
|:---------|:--------------:|:--------------------:|:--------------:|:-----:|
| LiH  (n_max=2) | 30 |   805 | 47 099 | **0.017×** |
| BeH  (n_max=2) | 30 |   805 | 47 099 | **0.017×** |
| CaH  (n_max=2) | 20 |   534 | 47 099 | **0.011×** |

At matched Q=18 (extrapolation via Paper 14 §IV.B's O(Q^2.5) law × T3 rel/scalar 2.4×), the projected advantage is **150×–250×**. Per-molecule Sunaga numbers would sharpen this.

### 2.3 Fine-structure sanity (T4)

| System | Reference (MHz) | GeoVac Z_eff (MHz) | Sign? | OoM? | Rel error |
|:-------|:---------------:|:------------------:|:-----:|:----:|:---------:|
| Li 2²P_{3/2}−2²P_{1/2} | 1.005×10⁴ | 3.127×10⁴ | ✓ | ✓ | +211 % |
| He 2³P total span      | 3.191×10⁴ | 1.095×10⁴ | ✓ | ✓ | −66 %  |
| Be 2s2p ³P span        | 7.260×10⁵ | 1.583×10⁵ | ✓ | ✓ | −78 %  |

All three are correct sign + correct order of magnitude (the Tier 2 Explorer T2-3 sanity criterion). The 20–50 % absolute-accuracy target from the sprint plan is **not** met for any of the three atoms — see §4 flag.

## 3. New ingredients for Paper 18

T5 delivered the spinor-intrinsic subtier. The updated operator-order × bundle structural grid now reads:

|                           | **2nd-order operator**               | **1st-order operator**                            |
|:--------------------------|:------------------------------------:|:-------------------------------------------------:|
| **scalar bundle**         | calibration π (even ζ)               | Tier-1 odd-zeta (ζ_R(3), ζ_R(5), …)               |
| **spinor bundle**         | *(not yet encountered — open slot)*  | **Tier-2 spinor-intrinsic (α², γ)** — NEW         |

Three of four cells populated. The empty cell (2nd-order on spinor bundle) would describe e.g. a squared Dirac operator or a Pauli–Villars regularized 2-spinor action; **flagged as a conjectural slot for future QED-on-S³ work** (Tier 3+ or beyond).

## 4. Honest framing reminders for T6 downstream drafts

The paper drop-ins (Paper 14 §V, Paper 20 Tier-2 table, Paper 22 spinor section, Paper 18 subtier) must respect the following constraints, which Tier 2's honest framing has surfaced:

1. **Not a new Dirac Fock projection.** Tier 1 Explorer Finding 2.1 (BJL Z₂-supersymmetric algebra) obstructs lifting the Schrödinger Fock map to the Dirac sector. Tier 2 builds spin-orbit *corrections* on the existing scalar S³ graph via Szmytkowski angular tables and (reserved) Martínez-y-Romero radial recursions. Papers must state this explicitly when describing the physical setup.

2. **Resource advantage is real; accuracy is not spectroscopic.** GeoVac composed architecture has known R_eq errors of 5–26 % for the molecules involved (Paper 17). Fine-structure splittings match sign + OoM only, 66–211 % relative error. The claim is structural resource sparsity at single-point fixed geometries, not a replacement for Dirac-Coulomb CASSCF or relativistic FCI.

3. **Sunaga SI Tables S1–S3 deferred.** Only the RaH-18q Pauli count is available from the main paper. Per-molecule (BeH, MgH, CaH, SrH, BaH) comparison is flagged as future work requiring SI extraction.

4. **Breit-Pauli scope is narrow.** T2 implements the leading-order α²/(n³·l(l+½)(l+1)) × L·S spin-orbit term only. γ radial corrections, Darwin term, and mass-velocity are explicitly out of scope. Reserved symbols exist (α, γ in T1) but are not bound to any T2/T3 matrix element.

5. **No TC combinations, no S⁵/S⁷, no graph Dirac operator.** All standing Tier-1/Tier-2 guardrails carry through. CUSP-3 (TC dead in second quantization) remains enforced.

## 5. Recommendation for Tier 3

Tier 2 closes the spin-orbit engineering upgrade. Four natural Tier 3 extensions, in priority order:

1. **Frozen-core [Kr] infrastructure for SrH / RaH.** Required to push into the Sunaga comparison range at native (not extrapolated) Q=18. Same frozen-core template as v2.2.0 ([Ne], [Ar]) but at Z=36 cutoff. Would enable point-by-point Sunaga comparison at matched Q.

2. **γ radial corrections via Martínez-y-Romero recursions.** T1 reserved α and γ; T2 binds α² only. Tier 3 would bind γ by importing Martínez-y-Romero's Dirac-Coulomb three-term radial recurrence and wiring it into `dirac_matrix_elements.py`. R_sp ring already accommodates γ (T5); the certifier would accept. First-row impact is ~1 % (Zα ≪ 1); third-row and heavier impact scales as (Zα)².

3. **Darwin + mass-velocity + SS / SOO** for sub-MHz fine structure. To convert the 66–211 % OoM-correct splittings into 1–5 % spectroscopically meaningful numbers. Multi-electron SS/SOO requires 6j-recoupling-level infrastructure in the two-body block; Darwin and mass-velocity are single-body diagonal additions compatible with T2's structure.

4. **QED-on-S³.** Populating the empty lower-left cell (2nd-order on spinor bundle) would require a squared Dirac or Pauli–Villars construction. No immediate motivation, but the Paper 18 grid is explicit about this being the open slot — any future QED extension would naturally slot here and complete the taxonomy.

**Do not open Tier 3 immediately.** Tier 2 is publication-ready; the PI's decision on commit timing for the T6 proposals determines whether Tier 3 follows in the same release or later.

---

## 6. Files (T6 scope)

| File | Purpose |
|:-----|:--------|
| `docs/tier2_verdict.md` | this memo |
| `docs/paper14_section5_proposal.tex` | Paper 14 §V drop-in (spinor composed encoding) |
| `docs/paper22_spinor_section_proposal.tex` | Paper 22 new section (d_spinor(l_max)) |
| `docs/paper20_tier2_table_proposal.tex` | Paper 20 Tier-2 resource table |
| `docs/paper18_spinor_subtier_proposal.tex` | Paper 18 §IV subtier (T5 output, unchanged; cross-referenced) |
| `docs/claude_md_tier2_updates.md` | CLAUDE.md v2.10.0 → v2.11.0 edits |

No papers or CLAUDE.md modified in-place. All six files are drop-in proposals for PI review.
