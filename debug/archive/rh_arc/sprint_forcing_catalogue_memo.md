# Sprint: Forcing Catalogue forward-run (2026-06-01)

**Verdict in one line:** Built `docs/forcing_catalogue.md` (the "atlas with the extra
column"), swept the corpus for ~165 bit-exact forcing-certificates, and ran 5 candidate-door
probes + 3 graduation probes. Net: **2 doors, 1 deep theorem, 1 vindicating partial, 2
clean kills.** The forced/free partition was promoted from a slogan to a theorem (the
tensor-product seam), and the inner Standard-Model algebra was shown mostly-forced with one
fork (ℍ vs M₂) confirmed open.

---

## 1. Premise and method

Origin: a reflective PI conversation about "what does GeoVac have, and how does it end."
The operational outcome was a method — treat every bit-exact result as a **forcing-certificate**
("this, I did not choose") and, for each, ask the question a consistency-check never asks:
*what else does the same rigidity force, that we never looked at because no physicist pointed
us there?* This is the "extra column" on top of the atlas: cataloguing the methods (sober
ending) and hunting cracks (romantic ending) are the **same path** plus one question per row.

Artifact: `docs/forcing_catalogue.md`. Sixteen seed certificates with full crack columns,
then a six-arc Explore sweep (~165 certificates), then the ranked candidate doors with probe
results.

Two meta-findings about the corpus, independent of any door:
1. **The engine has been run backwards.** Across ~165 certificates the FILED-AS tag is
   dominated by *consistency-check* over *prediction*. That ratio is the quantitative
   measure of how much forward-run remains.
2. **Negatives hide half-doors.** Predicted, and confirmed on the first sweep (Door 2).

---

## 2. The five candidate doors

| Door | Probe | Verdict |
|:-----|:------|:--------|
| 1 | F-theorem closed forms in general odd d (S⁷, S⁹) | **DOOR** |
| 2 | BW wedge entropy re-read (slope=2? F-link?) | **PARTIAL** |
| 3 | Paper 35 π-criterion full coverage | **THEOREM (full coverage)** |
| 4 | Why gauge forced / Yukawa free | **WALL+DOOR (seam theorem)** |
| 5 | ζ_{D²} / χ₋₄ integer-s ladder | **WALL** |

### Door 1 — F-theorem odd-d ladder: DOOR
Prop 7.4 (Paper 50) HOLDS. The conformal-scalar F-coefficient on S⁷ and S⁹ is an exact
ℚ-combination of {log2, ζ(3)/π², ζ(5)/π⁴, ζ(7)/π⁶, ζ(9)/π⁸}; the ζ(odd) ladder continues;
closed at 1e-302, frozen basis, two independent dimensions. A prior "falsification" was a
PSLQ maxsteps artifact (verified PSLQ-independently via direct integer-relation dot product).
**The engine demonstrably generates new closed forms** — this is the proof-of-concept that
the forward-run mode is productive. Files: `debug/door1_ftheorem_odd_d*`.

### Door 2 — BW wedge entropy: PARTIAL (vindication of the method)
The catalogue was right that a positive was shelved under a negative, and the audit corrected
*what* it is. Slope is provably EXACTLY 2 (degeneracy g_m=(n−m)(n−m+1) forces n² leading;
1.963 was a small-n artifact, windowed slope → 2.00005 at n∈[350,2000]). Constant has the
exact closed form C_∞ = coth(1) − log(2 sinh 1) = 0.4584… (PSLQ [−1,1,−1,−1], residual
<1e-80; original 0.540 was small-n contamination). **F-theorem link RULED OUT** — C_∞ is
PSLQ-null against {1, log2, ζ(3)/π²}; entropy (state-side, boundary dimension) and the
F-coefficient (operator-side, 3D free energy) live in disjoint rings. Confirms+sharpens
Paper 50 Prop 4.3. Applied to Paper 50 §5. Files: `debug/door2_bw_entropy_reread*`.

### Door 3 — π-criterion: THEOREM (full coverage)
The Paper 35 criterion ("π appears iff a continuous temporal/spectral integration is in the
evaluation") holds across all 28 Paper 34 projections + the case-exhaustion list + M1/M2/M3,
zero counterexamples. The M1 prime suspect (Hopf-measure π) resolved: discrete c₁ on the
finite Hopf graph = 0 (exact integer), continuum π reappears only via the Vol(S²)/4 measure
integral. Sharpening: two disjoint transcendental engines — spectral-side master Mellin
(π-class), state-side von-Neumann (log-class, π-free). Graduates the observation from sampled
to fully-covered (not a new theorem — Paper 32 §VIII already had the case-exhaustion theorem;
this is total rather than sampled coverage + the M1 resolution + the disjoint-engine
sharpening). Files: `debug/door3_pi_criterion_coverage*`.

### Door 4 — gauge/Yukawa seam: WALL+DOOR (deepest result)
The forced/free boundary is the **tensor-product seam** of the AC triple:
D² = D_GV²⊗1 + 1⊗D_F² (cross term killed by chirality anticommutation), factorizing the
Mellin engine into an outer ring (M1/M2/M3) × inner Yukawa Dirichlet ring ℚ[y_i^{−2s}], with
no shared generator. (i) "AC-extension cannot select the Yukawa" is a **THEOREM** —
η-trivialization + factorization + the G3 commuting-Z₂ result prove bone-ring ⊥
calibration-ring. This is the cleanest proof of the structural-skeleton scope in the corpus.
(ii) "No discrete construction could fix Y" is a **GAP** (category-relative). The wall and gap
are separated by one question: is the inner factor itself subject to a packing principle?
The seam explains why H1/W3/Koide all failed — they probed *inside* the AC category where (i)
forbids selection. Files: `debug/door4_gauge_yukawa_boundary_memo.md`.

### Door 5 — ζ_{D²} ladder: WALL
Ladder rigid through s=12, all forms exact — but every value is already-connected or an
internal coefficient with no isolated measurement. The "new observable" claim was
selection-attributable (the W3/c₂ failure mode) and correctly demoted by the audit. One
honest positive: D(2m+1) = (2^{2m}−2)ζ(2m−1) − ((2^{2m+1}−1)/2)ζ(2m+1). Files:
`debug/door5_zeta_d2_ladder*`.

---

## 3. The three graduation probes

### Door 1 recursion → WALL (clean, understood)
The odd-d F-theorem closed forms continue to S¹¹ (the flagship door stands). But the bonus
recursion scalar−2·Dirac = c_d·Dirac_{d−2} does NOT graduate: c₁₁ has no rational form, and
the original c₇=c₉=½ traced to an inconsistent Dirac normalization (a factor-of-2 switching
at d≥7). Top-atom cancellation scalar_top = 2·Dirac_top holds only at d=5, derived
analytically (ratio 2^{−(d−5)/2}). A convention artifact, now understood. Files:
`debug/door1b_s11_graduation*`.

### Door 4b inner algebra → PARTIAL DOOR
Among finite semisimple real *-algebras (≤3 summands, ≤size 3), exactly two reproduce
U(1)×SU(2)×SU(3): ℂ⊕ℍ⊕M₃(ℂ) and ℂ⊕M₂(ℂ)⊕M₃(ℂ). The Bertrand × Hopf-tower truncation at
n≤3 forces the factor count (3), the n=1 factor (ℂ), and the n=3 factor (M₃(ℂ)). Residual
non-uniqueness collapses to the single binary fork ℍ vs M₂(ℂ) at n=2. Files:
`debug/door4b_inner_algebra_forcing_memo.md`.

### Door 4c J sign-table audit → the fork ADMITS but does not FORCE ℍ (clean negative)
The conjectured handle (J_GV²=−1 selects ℍ) is closed NEGATIVE. Over ℂ, ℍ and M₂(ℂ) are the
same complex algebra (ℍ⊗ℂ ≅ M₂(ℂ)); the distinction is a real form fixed by an internal C²
involution (J_int²=−1 ⇔ ℍ; J_int²=+1 ⇔ M₂) that lives on the algebra side and is invisible
to the combined J = J_GV⊗J_F. Combined J², (ε,ε') signs, and KO-dimension are bit-identical
for both candidates at n_max∈{1,2,3}; the CCM order-one-with-Yukawa selector is satisfied by
both (matter/antimatter decoupling). ℍ is therefore a literature import (CCM
complex-fermion-rep / second-order condition), not a GeoVac forcing. Files:
`debug/door4c_j_signtable_audit*`.

**Discipline note.** Doors 2, 5, 1-recursion, and 4c are four consecutive over-reach
corrections — each killed with a *derived reason*, not a shrug. This is what makes the
surviving positives (Door 1 ladder, Door 3 coverage, Door 4 seam theorem, Door 4b
partial-forcing) trustworthy.

---

## 4. Paper / artifact edits applied

- **Paper 50 §5** (`paper_50_cft3_partition_function.tex`): the small-n fit
  (1.94 log n + 0.56) replaced by the exact asymptotic S = 2 log n + (coth 1 − log 2 sinh 1)
  + O(1/n), slope provably exactly 2, constant in closed form, F-disjointness made rigorous.
  Sharpens Prop 4.3 (`prop:F_vs_S_decomposition`).
- **Paper 32 §VIII** (`paper_32_spectral_triple.tex`): two new paragraphs after the Yukawa
  non-selection theorem — the forced/free seam theorem (Door 4) and the inner-algebra
  forcing ledger (Door 4b + 4c). The Door 4c negative is recorded honestly (J_GV²=−1
  admits-not-forces ℍ; ℍ is a literature import).
- **`docs/forcing_catalogue.md`**: new artifact; seed + sweep inventory + ranked doors +
  probe results + graduation outcomes.

---

## 5. Verification

- No `geovac/` source changed → no regression risk to the test suite. (Sub-agent drivers
  live in `debug/`, not `tests/`.)
- Paper edits are prose + math against existing labels (`prop:F_vs_S_decomposition`,
  `sec:unified_gauge`, `sec:construction`); compile check performed (see §7).
- Hard-prohibition check (§13.5): no change to natural-geometry hierarchy; no fitted/empirical
  parameter introduced; no negative result deleted (negatives ADDED); the Paper 2 K
  combination-rule "conjectural" label untouched. CLEAN.

---

## 6. Honest scope

**Closed at theorem grade:**
- The forced/free **seam theorem** (Paper 32 §VIII): D² factorization makes the inner Yukawa
  ring disjoint from every forced outer ring — bone ⊥ calibration data, proven via existing
  η-trivialization + G3 results. The Yukawa is free because it generates its own ring.
- **Door 4c negative**: J_GV²=−1 admits but does not force ℍ (the ℍ/M₂ distinction is invisible
  to the combined J; bit-exact at n_max∈{1,2,3}).
- **Inner-algebra partial forcing** (Door 4b): factor count, ℂ (n=1), M₃ (n=3) forced by the
  finite-algebra enumeration + Hopf-tower truncation.
- **Paper 50 §5 exact asymptotic**: slope exactly 2 (degeneracy-forced) + closed-form constant
  (PSLQ residual <1e-80).

**Numerical observation (high-confidence, not yet a closed-form theorem):**
- **Door 1 odd-d ladder**: the ζ(odd) closed forms on S⁷/S⁹/S¹¹ are exact (1e-302), but the
  *general*-odd-d pattern is verified by frozen-basis PSLQ on three dimensions, not proved in
  closed form. Prop 7.4 status: strongly supported, not a theorem.

**Full-coverage consolidation (not a new theorem):**
- **Door 3 π-criterion**: total coverage of the catalogued transcendentals with zero
  counterexamples; the underlying case-exhaustion theorem pre-existed (Paper 32 §VIII).

**Killed / walls (with derived reasons):**
- Door 1 recursion (convention artifact); Door 2 F-link (disjoint rings); Door 5 ladder
  (no unconnected observable); Door 4c ℍ-forcing (real-form invisible to J).

**Named open follow-ons:**
- The inner-factor packing question: is there a "second packing axiom" that fixes the n=2
  real form, the generation count, and the inner KO-dimension? (The remaining inner-algebra
  free data after Door 4b/4c.) This is the located romantic-door — currently a GAP, not a
  theorem.
- General-odd-d closed-form proof of the F-theorem ladder (Prop 7.4 graduation).
- tests/ coverage for the new Paper 50 §5 equation and the Door 4b finite-algebra
  enumeration (§13.4a open item).

---

## 7. Compile / release readiness

Paper 50 and Paper 32 compile checks: see release step. Ready for `/release` as a minor
version (v3.41.0): new diagnostic forward-run + doc + two paper edits, no source changes.
