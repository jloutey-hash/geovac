# Sprint: /aha cross-corpus scan + 7-track parallel investigation (2026-08-21)

**Origin:** PI-invoked /aha generative pass (plan: `debug/aha_cross_corpus_scan_2026_08_21.md`), then PI direction
"all in parallel" — 7 subagent tracks (T1, T2a, T2b, T3, T4a, T4b, T4c). This is the ONE canonical memo; per-track
detail lives in `debug/aha_track*_findings.md` + drivers `debug/aha_t*_*.py` + `debug/data/aha_*`.
**No paper or CLAUDE.md edits were made** — recommended edits are listed in §9 for PI decision.

## 1. Verdict table

| Track | Question | Verdict |
|:--|:--|:--|
| T1 | Paper 60 conditioning ↔ v4.73.0 commutator = one σ-spectrum? | **BORDERLINE gate / MAJOR positive**: theorem PASS 5.7e-15; exponents DERIVED = exactly 2; Woodbury leg NEGATIVE |
| T2a | Seam involution/double-cover mechanism | **BORDERLINE**: hypothesis's ρ-leg wrong (s↔t ⇒ ρ↦1/ρ, not 1−ρ); T-1 map tautology-blocked; **T2 lives on X₀(2)** (promotable); seam surprise quantified ~1 bit |
| T2b | Conductor-4 base-rate census | Blind trials **0/30**; forced mechanisms **3/3** ⇒ constrained landscape, CONVERGENT quantified |
| T3 | Recursive taxonomy: algebraic resurgent skeleton | **GO 4/4**: E₁ seed + hybrid ERI both algebraic-up-to-π-power; forced→free boundary-period transition localized between weight 1 and genus 1 |
| T4a | Defect↔polarizability; banked n_max=4 decider | UNDERPOWERED (n=4 real systems); decider **cost-walled** (16.0M dets, 21.7× jump) — needs matrix-free Davidson CI |
| T4b | Bargmann/S⁵ circle CM check | **SPLIT**: CM Hodge structure ABSENT (7th-asymmetry-layer candidate) BUT conductor-4 period content PRESENT via metaplectic 3/2: Z₋(4)=2(β(4)−G) = ¼ of P28's value |
| T4c | Where does T2's residual live (principled vs inherited)? | **MIXED**: ~14-digit outer plateau = inherited (s,t)-chart artifact; second-cusp Gevrey-1/Stokes structure = principled, chart-invariant |

## 2. T1 — principal-angle unification (drivers `aha_t1_*.py`, findings `aha_track1_findings.md`)

- **Theorem (now verified 5.7e-15, production matrices):** SW intra-block exactly I ⇒ spec(S)={1±σ_k},
  **cond(S)=(1+σ_max)/(1−σ_max)**, ‖[P_A,P_B]‖=max σ_k√(1−σ_k²). Both walls are moments of one
  canonical-correlation spectrum.
- **Exponents derived, = 2 exactly.** Cross-block = finite section of multiplication by W(χ)=j₀(kR·cot(χ/2)) ⇒
  band-limited concentration: **1−σ_max = c_sym·π²/n²**, c_sym=(kR)²/24 (SW), R_OH²/24−R_HH²/96 (water A₁).
  Published N^1.85 / N^1.97 are pre-asymptotic windows (local exponent →1.986/1.985 at n=160; prefactor match
  0.986–0.991). Collapse: (1−σ_max)(n/kR)² → π²/24 = 0.41123 across all R. Reconciles the memo-1.81 vs paper-1.85
  discrepancy (window effect).
- **Gerade lever = exact constant:** cond(gerade) → 2/(1+min_x j₀(x)) = **2.555041…** for every R, N
  (measured 2.553–2.554 at n=160). Paper 60's "flat κ≈2" is this constant at small n.
- **Woodbury/hot-key NEGATIVE with mechanism:** O↔H σ-spectrum not low-rank (participation ratio 7.40; Szegő ⇒
  rank grows ∝ n). Since 1−σ_k ≈ c_sym(kπ)²/n², fixed rank-r divides κ by (r+1)² only; flattening needs r∝n
  (= the Löwdin dead-end). Constant-factor d_eff wins 39–74× at fixed N, but d_eff ~ N^1.81 vs d_base ~ N^2.10.
- **Commutator-norm deflation:** over 216 sweep points, once cond>20 the commutator pins at **0.484±0.025**
  (algebraic ceiling 0.5 at σ=1/√2). The v4.73.0 "‖[P_A,P_B]‖=0.50" was the saturation value — the
  non-commuting-projections *diagnosis* stands; the *norm* carries no severity information.
- SW metric has **no Z-dependence** (one-parameter family in kR); Z enters only the L²/Goscinskian metric
  (measured α≈3.6→3.1, worse, different mechanism; large-R Z=(1,8) rows R²=0.84, not quotable as exponents).

## 3. T2a — seam involution (driver `aha_t2a_involution.py`, findings `aha_track2a_findings.md`; 44/44 checks)

- **Correction to the /aha hypothesis:** c₁=s(1−s), c₂=t(1−t), ρ=c₂/c₁ ⇒ s↔t ⟺ **ρ↦1/ρ** (λ↦λ/(λ−1),
  T-coset fixing the cusp), NOT ρ↦1−ρ. Full Feynman-square symmetry group induces exactly {id, ρ↦1/ρ};
  τ=i / ρ=1/2 is not a fixed point of anything physical (it is half of the orbit {1/2, 2}).
- Elliptic lift real: M=[[1,−1],[2,−1]], M²=−I, fixed τ*=(1+i)/2 (disc −4, j=1728, λ=2) — **off the physical
  domain**. Physical fixed locus (t=s, t=1−s) = ρ=1 = the cusp = pure-Tate nodal degeneration.
- **Promotable positive:** ⟨Γ(2), M⟩ = Γ₀(2); the co-area fold Φ(ρ)=Φ(1/ρ)/ρ² IS the s↔t statement ⇒
  **T2 = (16/π)∫₀¹Φ dρ integrates over a Γ₀(2) fundamental domain — T2's modular home is X₀(2), not X(2)**.
  Symmetry-forced, fibre-independent.
- **T-1 map tautology-blocked structurally:** route (a)'s object is a complex structure on a compact-group rep
  (not H¹ of anything, no Galois action); no functor to H₁(curve). Also the audit's own T-1 polarization
  discriminator is **void**: h(ℚ(i))=1 ⇒ the polarization leg passes automatically even for tautological maps.
  T-1 must be restated as functoriality (exhibit V_fund as H₁ of a GeoVac-constructed abelian variety).
- **Deflation, full strength:** PSL₂(ℤ)≅ℤ₂∗ℤ₃ ⇒ any order-2 symmetry of a modular family pins a ℤ[i] point
  (all transposition fixed points have j=1728), order-3 pins ℤ[ζ₃]. "Two ℤ₂'s land on ℚ(i)" ≈ two coin flips
  agreeing (~1 bit). Sharpest CONVERGENT statement to date. Stabilizer leg verified exactly: S acts on H₁(E_i)
  as mult-by-i with matrix = P56's J, Ω = P56's Q, ΩJ=I.

## 4. T2b — denominator census (`debug/data/aha_period_census.md`, findings `aha_track2b_findings.md`)

- **Population A (blind PSLQ/scan trials, N=30): 0 confirmed disc-4 hits** — including the two campaigns that
  targeted disc-4 on purpose (W1e basis; S_min χ₋₄-depth-2 sweep, 35 attempts).
- **Population B (forced constructions, N=9):** conditioned on reaching beyond pure-Tate at all, disc-4 landing
  = **3/3** (vertex-parity M3; P56 Kramers-CM; P59 elliptic CM fibers). Paper 34's 28-row dictionary
  cross-checks: exactly 1/28 rows carries disc-4 (the same vertex-parity route).
- **Verdict: constrained landscape** — ~3 known roads out of pure-Tate, each landing on the first CM field its
  own rational arithmetic reaches. Quantitative version of CONVERGENT. The one live upgrade test remains
  β(2)-in-T2 at ≥32 digits (does disc-4 survive to the integrated observable).

## 5. T3 — resurgent Stokes battery (drivers `aha_t3_*.py`, findings `aha_track3_findings.md`) — GO

- **Object 1, E₁ seed:** Borel transform 1/(1+ζ), one simple pole ζ=−1, Stokes constant **2πi·1** — derived and
  verified 3 ways (9.7e-63 / 4.4e-61 / 5.0e-46). Trans-series terminates.
- **Object 2, hybrid two-centre ERI {E₁, ln} (9 quartets, Z_A 2–5, Z_B 1–3, l 1–2):** every Borel singularity at
  ζ=±a_d (far-centre exponent), 36/36; Stokes constants **2πi×ℚ(√d)** — blind large-order extraction matches
  closed form to 35–81 digits, PSLQ exact 8/8, branch-cut cross-check 16/16 (≤2.3e-46). **Logs do not enter the
  Stokes data**: entire log content obeys L(R) = −P(R)·ln Λ with Λ = the multiplicative cross-ratio of the four
  Borel actions ∈ ℚ (9/9, decoy-controlled, mechanism derived). Bonus: sector Stokes constants cancel pairwise ⇒
  physical ERI **cut-free** in R (|Disc F|/|F| ≤ 2.6e-49) — median resummation IS the closed form.
- **Pattern 4/4** (E₁ seed, hybrid ERI, N(D), T2 cusp): Stokes data algebraic up to a π-power, and **the π-power
  is fixed by the Borel-singularity type** (pole→2πi; √-branch→π⁰; 3/2-branch→1/√π). The boundary period is the
  non-uniform part: absent (E₁) / forced ln Λ (hybrid, weight 1) / free K(1−ρ) (N(D), genus 1) / open Γ(2) MMV
  (T2). **The forced→free transition sits between weight 1 and genus 1** — the sharp form of the proposed new
  taxonomy axis: grade each Layer-2 transcendental by (field of Stokes data, weight/genus of boundary period,
  whether the latter is forced by the former).
- Honest cap: Object 2's resurgent data derives from the validated closed form (large-order route = internal
  cross-check, not a second determination). **Named next falsifier:** the exchange-class γ — scoping shows
  coeff(γ)=coeff(ln R) identically (5/5; γ never multiplies E₁) but the R-dependent log ⇒ logarithmic Borel
  branch point ⇒ separate computation, left unclassified.

## 6. T4a — probes (drivers `aha_t4a_*.py`, findings `aha_track4a_findings.md`)

- **Defect↔polarizability: UNDERPOWERED.** Only 4 genuinely distinct systems carry an R_eq defect; at n=4 even
  ρ=1 fails significance (exact permutation p=0.083). LiH variants are pseudoreplication (defect 1.5%→63.5% at
  one polarizability). Observed ρ≈+0.13/+0.20 (n.s.). Polarizability values agent-knowledge, confidence-tagged.
  WH7 de-compactification reading neither supported nor killed — untestable on this atlas; would need per-system
  defects on a common solver.
- **Banked n_max=4 A-vs-C decider: cost-walled.** Sector dim 741,321 → 16,040,025 (21.7×); already ~2.3 h/point
  at n_max=3 with the pure-Python sparse CI builder. Run attempted, guard invoked, not restarted. **Unblock =
  matrix-free (Davidson) CI**, not a bigger probe. Decider stays open/banked.

## 7. T4b — Bargmann/S⁵ circle (driver `aha_t4b_bargmann_circle.py`, findings `aha_track4b_findings.md`) — SPLIT

- **Leg A, CM/Hodge structure: ABSENT** — three independent exact obstructions (positive generator ⇒ no ±-paired
  spectrum / traceless lemma; i^N ℚ-rational exactly where it squares to +1; antiholomorphic half discarded by
  the Hardy construction). ⇒ **7th Coulomb/HO asymmetry layer candidate** (draft in findings §4, incl. an
  explicit companion NON-layer for χ₋₄).
- **Leg B, conductor-4 period content: PRESENT** (unanticipated; ~5e-51, three routes): ℤ₂-graded Hardy zeta
  **Z₋(s) = 2^{s−3}(β(s)−β(s−2))**, so Z₋(4) = 2(β(4)−G) = **exactly ¼ of P28 thm:chi4's 8(β(4)−G)**.
  Ablation computed: integer shift → pure Tate; any half-integer shift → χ₋₄ (not χ₋₃). Source = metaplectic
  zero-point 3/2 × bipartite grading — **no spin required**.
- Ungraded partition function = 1/(1−q)³ (Hilbert series of ℂ[z₁,z₂,z₃], 𝔸³, pure Tate, conductor 1) —
  separated from θ₃² by recurrence/multiplicativity/conductor tests.
- **Sharpened discriminator (three-track coherent, see §8):** conductor-4 *periods* need only an order-2 datum
  lifted to order 4 (spin OR metaplectic OR modular torsion); a CM *Hodge structure* additionally needs the
  spin-specific ±-paired ℚ-rational carrier. "Not with compactness" confirmed and strengthened (third compact
  circle pure-Tate on the ungraded ledger).
- Honest fences: leg B is mechanism transfer, not a measured observable; normal-ordering the 3/2 away kills it
  (the Dirac chirality-forced shift cannot be normal-ordered away — the residual asymmetry a reviewer would push).

## 8. Cross-track synthesis (the actual /aha yield)

1. **The conductor-4 story is now one coherent statement across T2a+T2b+T4b:** ℚ(i)/conductor-4 content appears
   wherever an order-2 datum acts through its order-4 lift — three independent instances (spin double cover;
   metaplectic half-integer shift; PSL₂(ℤ) torsion) — and this is a *low-information* event (~1 bit; the
   landscape's first exit from pure-Tate). CM *Hodge* structure is strictly harder: it needs the ±-paired
   ℚ-rational (quaternionic) carrier, which only the spin side has. The QI-SEAM verdict CONVERGENT is not just
   sustained but *explained and quantified*. Untraveled-road prediction: an order-3 physical symmetry on a
   modular family should pin disc −3/ℚ(ζ₃) (μ₃ grading/counterfactual computed in T4b: sixths, conductor | 12;
   the substrate has no ℤ₃ grading — so absence is itself the expected outcome).
2. **T1's derivation converts Paper 60's two empirical exponents into one theorem + one constant** and retires
   the commutator norm as a severity metric while confirming the diagnosis it stood for.
3. **T3's forced/free boundary-period split** gives the recursive-taxonomy claim its precise shape and locates
   the skeleton/projection boundary *within* resurgence at weight 1 → genus 1 — consistent with the corpus-wide
   pattern (decidability ends where genus begins, cf. P58/P59).
4. **T4c legitimizes the digit plan:** the ~14-digit plateau was chart-inherited (already relocated); only the
   second-cusp Stokes structure is invariant — so ≥32 digits (the β(2) seam test) is a resummation problem, not
   a quadrature problem, confirming v4.100.0's call and striking the K₀/Parseval chart from the candidate list.

## 9. Recommended edits (PI-gated; NOT applied)

1. **Paper 60** (certified — edits would compound a re-review): replace fitted N^1.85/N^1.97 with the derived
   1−σ_max = c_sym π²/n² law (+ collapse constant π²/24); state cond(S)=(1+σ_max)/(1−σ_max) via the exact
   intra-identity; replace "flat κ≈2" with the exact gerade constant 2.555041…; add the Woodbury negative
   (fixed-rank ⇒ (r+1)² only) closing the hot-key idea; new backing tests `test_paper60_sigma_law` etc.
2. **Paper 59 sec:modular:** add the X₀(2) statement (s↔t = Fricke-type ρ↦1/ρ; co-area integral over a Γ₀(2)
   fundamental domain) + backing test; optionally the T4c chart-audit sentence (outer plateau inherited,
   second cusp principled).
3. **Seam audit memo / P56 rem:paper59_cm:** restate T-1 (polarization leg void by h(ℚ(i))=1 ⇒ functoriality
   test instead); weaken the τ=i leg per T2a (τ=i is not a physical fixed point); record the census numbers
   (0/30 blind, 3/3 forced) as the quantitative CONVERGENT.
4. **Paper 56 rem:kms_hodge_circle:** qualification-2 discriminant is necessary-not-sufficient (HO μ₄ faithful,
   still no CM) — add the ±-paired-carrier condition.
5. **Paper 24:** 7th asymmetry layer (CM Hodge structure) + the companion statement that conductor-4 *period*
   content transfers (Z₋(4)=2(β(4)−G), ¼ of P28) via the metaplectic shift; fix Layer-6 wording that invites
   conflation with the zero-point half-integer. Transcendental tags: β(4), G already catalogued (P28 χ₋₄ family);
   the new appearance is the same M-class vertex-parity/χ₋₄ projection, now via metaplectic shift — Paper 34
   chain note recommended.
6. **Paper 18/34 (bigger, separate decision):** the resurgent-skeleton axis (Stokes field, boundary-period
   weight/genus, forced-vs-free). New axis = structural change; suggest an [OBSERVATION]-tier subsection first.

## 10. Named follow-ons (not started)

- **Matrix-free Davidson CI** for the balanced builder → unblocks the banked n_max=4 A-vs-C decider (T4a).
- **Exchange-class γ resurgence** (logarithmic Borel branch point) — T3's sharpest next falsifier.
- **β(2)-in-T2 at ≥32 digits** via near-cusp spectral evaluator + median resummation (T4c says this is the only
  invariant wall; T3's cut-free result on the hybrid class is the proof-of-concept that median resummation can
  BE the closed form).
- Per-system common-solver defect atlas if the WH7 defect↔continuum question is ever to be powered (T4a).
