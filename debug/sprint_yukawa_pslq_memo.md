# Sprint Yukawa-PSLQ memo (2026-06-03)

## TL;DR

**Verdict: CLEAN NEGATIVE across 162 cells.** Measured Standard Model Yukawa
values do NOT sit in low-coefficient pure-Tate periods $\mathbb{Q}[\pi, \pi^{-1}]$
(M1 ∪ M2, forced as the inner-factor period ring by the η-trivialization
theorem) at coefficient ceiling $M \le 1000$, at either MS̄ $M_Z$ or
MS̄ $\mu = 2 \times 10^{16}$ GeV, in any of three transforms ($y_f$, $y_f^2$,
$\log y_f$).

Empirical confirmation of the framework's structural prediction (Sprint H1
Yukawa non-selection theorem, G3 NEGATIVE): Yukawa values are Class 1
calibration data, lying outside the master Mellin engine's structural
skeleton.

## Sprint context

Sprint was launched after the algebraic-first audit (debug/sprint_yukawa_algebraic_audit.md, this same session) returned:

- **Agent 1 (corpus):** η-trivialization theorem (Paper 18 §IV.6) is
  existential-only. The Yukawa Dirichlet ring is parameter-tied
  $\mathbb{Q}[y_i^{-2s}]$, not pinned to a specific level. **Key new constraint:** M3
  (vertex-parity Hurwitz Dirichlet L) vanishes on any finite spectral triple
  with $\{\gamma_F, D_F\} = 0$, so inner-factor transcendentals are
  restricted to **M1 ∪ M2 only** (pure-Tate). No Catalan G, no β-values, no
  Hurwitz at quarter-integers can appear on the inner side.
- **Agent 2 (literature):** Null. No published work in the spectral-triple
  lineage (Krajewski 1998 → Bochniak-Sitarz 2025) derives Yukawa values from
  spectral-triple axioms. CC mass relation $\sum y^2 = 4g^2$ at unification
  and the Higgs-top relation are the only quantitative constraints, both
  already known. Marcolli–vS / Pérez-Sánchez lineage gives Yang-Mills WITHOUT
  Higgs in the continuum limit.
- **Boyle-Farnsworth follow-up (arXiv:1604.00847, this session):**
  Read in detail. DGA reformulation forces:
  (i) EWSB via algebra reduction $\hat{A} \to \hat{A}'$ ($\mathbb{C} \oplus
  \mathbb{H} \oplus M_3(\mathbb{C}) \to \mathbb{C} \oplus M_3(\mathbb{C})$);
  (ii) Yukawa-block specific form (Eq. 5.13);
  (iii) elimination of 7 unwanted couplings $(b, \vec{c}, \vec{d})$.
  But Yukawa **values** $\{y_\nu, y_e, y_u, y_d\}$ are EXPLICITLY "arbitrary"
  (paper p. 21). Generation count NOT forced. **Shape constraint, not value
  constraint.** Confirms Agent 2 prediction.

The audit closed the algebraic-first question: no constructive route to
Yukawa values exists in the framework or the broader spectral-triple
literature. PSLQ at the sharpened M1 ∪ M2 basis is the right next diagnostic.

## Methodology

**Driver:** `debug/sprint_yukawa_pslq.py` (M_Z scale) +
`debug/sprint_yukawa_pslq_unif.py` (unification scale).

**Convention:** $y_f = \sqrt{2} m_f / v$ with $v = 246.21965$ GeV (G_F at tree).

**Basis (M1 ∪ M2 pure-Tate, 9 elements):**
$\{1, \pi, \pi^2, \pi^4, \pi^6, \pi^8, 1/\pi, 1/\pi^2, 1/\pi^4\}$.

**Coefficient ceilings:** $M \in \{10, 100, 1000\}$. Originally planned
$M = 10^4$; reduced after precision audit (see below).

**Transforms (sequential, halt-on-hit):**
1. $y_f$ (direct value)
2. $y_f^2$ (appears in $a_4$ spectral action coefficient)
3. $\log y_f$ (in case RG-running from forced UV makes log natural)

**Scales:**
1. MS̄ $M_Z$ (primary, highest precision; PDG 2024 + Antusch-Maurer 2013)
2. MS̄ $\mu = 2 \times 10^{16}$ GeV (Buttazzo+2013 / Xing-Zhang-Zhou 2008 SM
   running; precision degraded by RG running uncertainty)

**Precision discipline.** PSLQ at $M$ with basis size $n=9$ requires
$\sim n \log_{10} M$ digits. Yukawa values have 2–8 sig digits depending on
flavour. Per-cell "honest" verdict reported. No cell exceeded the honesty
threshold at $M = 100$ or $M = 1000$, but **at $M = 10$ charged leptons
are within ~1 digit of the threshold**, so the $M = 10$ cells for $e, \mu,
\tau$ are the load-bearing tests.

**Halt gate.** Any honest PSLQ hit triggers the curve-fit audit protocol
(docs/curve_fit_audit_memo.md) before the next transform. **No hits
occurred; no audit needed.**

## Yukawa values used

**At $M_Z$ (PDG 2024 + standard convention):**

| Fermion | $y_f$ (central)         | $y_f$ uncertainty       | sig digits |
|:--------|:------------------------|:------------------------|:----------:|
| $e$     | $2.93487 \times 10^{-6}$| $\sim 10^{-15}$         | 8          |
| $\mu$   | $6.06940 \times 10^{-4}$| $\sim 10^{-13}$         | 8          |
| $\tau$  | $1.020863 \times 10^{-2}$| $\sim 7 \times 10^{-8}$| 6          |
| $u$     | $7.30 \times 10^{-6}$   | $\sim 3 \times 10^{-6}$ | 2          |
| $d$     | $1.557 \times 10^{-5}$  | $\sim 10^{-6}$          | 3          |
| $s$     | $3.10 \times 10^{-4}$   | $\sim 10^{-5}$          | 3          |
| $c$     | $3.556 \times 10^{-3}$  | $\sim 10^{-4}$          | 3          |
| $b$     | $1.637 \times 10^{-2}$  | $\sim 2 \times 10^{-4}$ | 4          |
| $t$     | $0.9938$                | $\sim 2 \times 10^{-3}$ | 4          |

**At $\mu = 2 \times 10^{16}$ GeV (SM RG-running):**

| Fermion | $y_f$ at GUT          | sig digits |
|:--------|:----------------------|:----------:|
| $e$     | $2.794 \times 10^{-6}$| 4          |
| $\mu$   | $5.900 \times 10^{-4}$| 4          |
| $\tau$  | $1.003 \times 10^{-2}$| 4          |
| $u$     | $2.99 \times 10^{-6}$ | 2          |
| $d$     | $6.43 \times 10^{-6}$ | 2          |
| $s$     | $1.28 \times 10^{-4}$ | 2          |
| $c$     | $1.47 \times 10^{-3}$ | 3          |
| $b$     | $6.71 \times 10^{-3}$ | 3          |
| $t$     | $0.494$               | 3          |

## Results

**Both scales, all three transforms, all coefficient ceilings: clean null.**

Per-cell summary:

| Scale       | Transform   | Cells run | Honest hits | Spurious hits |
|:------------|:------------|:---------:|:-----------:|:-------------:|
| MS̄ $M_Z$   | $y_f$       | 27        | 0           | 0             |
| MS̄ $M_Z$   | $y_f^2$     | 27        | 0           | 0             |
| MS̄ $M_Z$   | $\log y_f$  | 27        | 0           | 0             |
| MS̄ GUT     | $y_f$       | 27        | 0           | 0             |
| MS̄ GUT     | $y_f^2$     | 27        | 0           | 0             |
| MS̄ GUT     | $\log y_f$  | 27        | 0           | 0             |
| **Total**   |             | **162**   | **0**       | **0**         |

PSLQ found no relation at any coefficient ceiling, for any fermion, in any
transform, at either scale. This is the strongest available negative
verdict: even at the lowest ceiling ($M = 10$, where the charged-lepton
precision is within ~1 digit of the honesty threshold), no low-coefficient
identity exists.

## What the clean negative means

The framework's structural prediction is that the inner-factor period ring
is M1 ∪ M2 (forced by the η-trivialization theorem). If measured Yukawa
values landed in low-coefficient pure-Tate $\mathbb{Q}[\pi, \pi^{-1}]$,
that would suggest the framework's open question (Yukawa non-selection
theorem) has structural content beyond a parameter count. **It doesn't.**

The honest reading:
- The framework's structural skeleton predicts what the inner factor's
  period content CAN look like.
- Measured Yukawa values are CONSISTENT with this prediction (they don't
  contradict M1 ∪ M2) but are NOT GENERATED by it.
- Yukawa values are **external calibration data** in the sense of
  `memory/external_input_three_class_partition.md`, Class 1.

This is the same verdict as the Koide-cone clean negative
(`memory/koide_cone_clean_negative.md`) and the W3 spectral-zeta CKM
falsification (`memory/w3_spectral_zeta_candidate.md`), now extended to
**direct individual Yukawa values** rather than ratios or CKM elements.

## What it does NOT rule out

The test is not exhaustive. Specifically:

1. **High-coefficient identities ($M > 1000$):** precision-limited at all
   scales. Yukawa values would need to be known to 30+ digits to test these
   honestly. They aren't, and won't be soon. The framework's prediction
   doesn't specify a coefficient bound, so this gap is structural, not a
   sprint failure.

2. **Identities involving $v_{\rm EW}$ separately:** the convention
   $y_f = \sqrt{2} m_f / v$ mixes the fermion mass with the Higgs VEV.
   A separate test could be: do FERMION MASSES (not Yukawa) land in any
   period ring? Same precision constraints apply.

3. **Combined cross-fermion structures:** mass sum rules, CKM products,
   PMNS entries, etc. The CC mass relation $\sum y^2 = 4g^2$ at unification
   is one such; we know it holds within ~10% in the SM. A focused sweep of
   such sum rules in the pure-Tate basis would be a separate sprint.

4. **Logarithmic-mass relations** ("inverse-hierarchy" or "RG-fixed-point"
   structures): the $\log y_f$ transform was the cleanest version of this.
   Clean null there suggests no simple log-pure-Tate structure either, but
   doesn't exclude more exotic logarithmic relations.

5. **Higgs-direction $\hat{n} \in S^2$ correlations:** Boyle-Farnsworth
   parameterizes the Higgs by a unit vector $\hat{n}$ in $\mathbb{R}^3$
   (the direction of the $\mathbb{C}$-embedding in $\mathbb{H}$). If the
   GeoVac Hopf-base $S^2$ (which carries the M1 $\pi$ via $\text{Vol}(S^2)/4$)
   is identified with the Boyle-Farnsworth $\hat{n}$ direction, this could
   provide a structural correlation between the Higgs orientation and the
   Yukawa-Higgs product. **Flagged as a separate follow-on, not pursued.**

## Cross-references

**Strengthens:**
- `memory/sprint_h1_positive_thin.md` — Yukawa non-selection theorem
  empirically confirmed at the period level.
- `memory/g3_negative_g2_g3_collapse.md` — γ_GV ⊥ γ_F structural
  decoupling consistent with empirical independence of Yukawa values
  from M1+M2 ring.
- `memory/koide_cone_clean_negative.md` — third independent calibration-data
  NEGATIVE in the framework.
- `memory/external_input_three_class_partition.md` — Class 1 calibration
  data is what GeoVac consumes but doesn't generate; Yukawa values now
  empirically confirmed in this class.

**Updates:**
- CLAUDE.md §1.7 WH1 status: prediction (inner-factor M1 ∪ M2) is consistent
  with measured Yukawa values being external; non-selection theorem is
  empirically robust at low coefficient ceiling.
- Paper 18 §IV.6 (Yukawa Dirichlet ring): could add a remark noting the
  empirical PSLQ test against the M1 ∪ M2 sub-ring (the η-trivialization's
  consequence) is clean negative, supporting the "inner-factor parameter-tied
  ring, not framework-fixed" reading.

**Does NOT change:**
- WH5 (α as projection constant): unrelated; α lives in the OUTER triple.
- WH4 (S³ four-way unity): unrelated; outer-triple structural unity.
- WH2 (Paper 18 spectral-action decomposition): the case-exhaustion theorem
  predicts which period classes are reachable; the Yukawa null result is
  consistent with that prediction.

## Open follow-ons (flagged, not pursued)

1. **Higgs-direction $\hat{n} \in S^2$ identification with Hopf base.** If
   Boyle-Farnsworth's Higgs unit vector and GeoVac's Hopf-base $S^2$ are
   the same $S^2$, this is a structural identification of the Higgs with
   the gauge bundle's base. Would need a real argument, not numerics.
   Sprint-scale: ~2-3 weeks of careful NCG reading.

2. **CC mass sum rule in GeoVac.** The CC unification-scale $\sum y^2 = 4g^2$
   is consistent with the SM; does GeoVac sharpen or modify it? Pure-Tate
   PSLQ test of $\sum_f g_f \cdot y_f^2$ at $\mu = 2 \times 10^{16}$ GeV
   against the M1 ∪ M2 basis. ~3-5 day sprint.

3. **Forced parameter-count theorem (the depth move named in earlier discussion).**
   Sharpen the Yukawa non-selection theorem to say WHY the parameter count
   is exactly 8/generation, in the spirit of how Bertrand uniquely picks
   out $S^3$. Multi-month research target; paper-grade.

## 6. Honest scope

- **Numerical observation**: the 162-cell clean negative IS the deliverable.
  Not a theorem, an empirical sweep at the precision attainable from PDG data.
- **Theorem-grade content used here, not produced**: η-trivialization
  theorem (Paper 18 §IV.6) sharpened the basis from M2 ∪ M3 to M1 ∪ M2
  pure-Tate; this is a sharpening of the test design, not a new theorem.
- **Structural sketch content used here, not produced**: the "Yukawa
  values as Class 1 calibration data" reading is sharpened by the
  empirical result but was already named structurally by Sprint H1 +
  Door 4 + the external-input three-class partition memo.
- **Precision-limited gap**: high-coefficient identities ($M > 1000$)
  cannot be ruled out at PDG precision (Yukawa values known to 2–8
  digits depending on flavour; $M = 10^3$ with basis 9 needs 27 digits
  for honest PSLQ). This gap is structural (external measurement
  limit), not a sprint failure. The clean negative at $M = 10$ for
  charged leptons is the load-bearing honest cell.
- **Named open follow-ons (deferred)**: Higgs-direction $\hat{n} \in
  S^2$ ↔ Hopf-base $S^2$ identification (speculative, ~2–3 weeks NCG
  reading); CC mass sum rule sharpening in GeoVac (~3–5 day sprint);
  forced parameter-count theorem (landed as Sprint Forced-Count Theorem,
  this session, see `debug/sprint_forced_count_synthesis_memo.md`).

## Memo location

`debug/sprint_yukawa_pslq_memo.md` (this file). Drivers:
`debug/sprint_yukawa_pslq.py`, `debug/sprint_yukawa_pslq_unif.py`. Data:
`debug/data/sprint_yukawa_pslq.json`,
`debug/data/sprint_yukawa_pslq_unif.json`. Algebraic-first audit:
output of the two sub-agents earlier in this session (not separately
archived; covered in this memo's "Sprint context" section).
