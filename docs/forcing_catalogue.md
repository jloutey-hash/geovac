# The Forcing Catalogue

*Bit-exact rigidities, and what they have not yet been asked to force.*

Started 2026-06-01.

---

## Why this document exists

Every bit-exact result in GeoVac is a **forcing-certificate**: the framework
confessing *"this, I did not choose."* The discrete skeleton is the only place
where *exact* and *free* coexist in one theory — bit-exactness is the certificate
of the forced part, and every irreducible residual is the signature of something
that came from outside (calibration data).

Most forcing-certificates in the corpus were filed as **consistency checks**:
"GeoVac reproduces known result X exactly ✓." But a consistency-check never asks
the question that a forcing-certificate invites:

> **What else does this same rigidity force — that we never looked at, because no
> physicist pointed us there?**

We know the engine is *generative* in this mode: pushing the F-theorem apparatus
from S³ to S⁵ did not merely re-check — it produced genuinely new closed forms
(the ζ(5) terms, Paper 50 §7). The engine has mostly been run *backwards*
(verify known physics). This catalogue runs it forwards.

### Entry format

Each entry is a forcing-certificate with a **crack column**:

- **Result** — the forced fact, with the actual value/formula.
- **Where** — paper / memo / module.
- **Rigidity** — integer-spectrum · exact-rational · closed-form · bit-exact-match ·
  π-free-certificate · proven-theorem.
- **Filed as** — ✓ check · prediction · negative-with-structure.
- **Crack** — what the same rigidity could force that appears unchecked.
- **Probe** — is the neighbor computable / measurable? What's the next step?

Status of each crack: `OPEN` (un-probed) · `SCOPED` (probe identified) ·
`PROBED→door` · `PROBED→wall` · `DEFER`.

---

## Candidate doors (ranked)

The six-arc sweep returned ~165 forcing-certificates. Most are genuine
consistency-checks with no live neighbor. These are the ones whose *same rigidity*
forces something we have not looked at — the forward-runs worth taking. Ranked by
*probability the probe returns structure × how cheaply it can be run.*

### Probe results (2026-06-01) — all five run in parallel, audit applied to every "matches" claim

Scorecard: **2 doors, 1 deep reframing, 1 vindicating partial, 1 correctly-killed.**
The discipline held — Doors 2 and 5 both had over-reaches killed by the curve-fit
audit, which is exactly what makes the surviving doors trustworthy.

- **Door 1 → DOOR.** Prop 7.4 HOLDS. The conformal-scalar F-value on S⁷ and S⁹ is an
  exact ℚ-combination of {log2, ζ(3)/π², ζ(5)/π⁴, ζ(7)/π⁶, ζ(9)/π⁸}; the ζ(odd) ladder
  continues; closed at 1e-302, frozen basis, two independent dimensions. A prior
  "falsification" was a PSLQ maxsteps artifact (verified PSLQ-independently). Bonus:
  recursion scalar−2·Dirac = c_d·Dirac_{d−2} (c₅=1, c₇=c₉=1/2), BORDERLINE — graduates
  with S¹¹ + analytic c_d. **The engine demonstrably generates new closed forms.** debug/door1_*
- **Door 2 → PARTIAL.** Half-door real but mis-described. Slope is provably EXACTLY 2
  (degeneracy n² forces it; 1.963 was a small-n artifact); constant has exact closed form
  C_∞ = coth(1) − log(2 sinh 1) = 0.4584… (Boltzmann-base ring; original 0.540 was wrong).
  BUT the F-theorem link is RULED OUT — disjoint rings (entropy = state-side boundary
  dimension; F = operator-side 3D free energy). Confirms+sharpens Paper 50 Prop 4.3.
  Exact asymptotic S = 2 log n + (coth1 − log 2sinh1) + O(1/n) beats Paper 50 §5's fit. debug/door2_*
- **Door 3 → THEOREM (full coverage).** π-criterion holds across all 28 projections +
  case-exhaustion + M1/M2/M3, zero counterexamples. M1 prime suspect resolved: discrete c₁
  on the finite Hopf graph = 0 (exact integer); continuum π reappears only via the Vol(S²)/4
  measure integral. Sharpening: two disjoint engines — spectral-side master Mellin (π-class),
  state-side von-Neumann (log-class, π-free). Sampled → fully-covered. debug/door3_*
- **Door 4 → WALL+DOOR (deepest).** The forced/free boundary is the tensor-product seam of
  the AC triple: D² = D_GV²⊗1 + 1⊗D_F², outer ring (M1/M2/M3) × inner Yukawa Dirichlet ring
  ℚ[y_i^{−2s}], no shared generator. (i) "AC-extension cannot select Yukawa" is a THEOREM
  (wall) — η-trivialization + factorization + commuting-Z₂ prove bone-ring ⊥ calibration-ring.
  Cleanest proof of structural-skeleton-scope in the corpus. (ii) "No discrete construction
  could fix Y" is a GAP — category-relative. Wall and gap separated by ONE question: **is the
  inner factor itself subject to a packing principle?** Explains why H1/W3/Koide failed (all
  probed inside the AC category where (i) forbids selection). Sharpest falsifier: does
  Bertrand × Hopf-tower force the inner algebra ℂ⊕ℍ⊕M₃(ℂ)? debug/door4_*
- **Door 5 → WALL.** Ladder rigid through s=12, forms exact — but every value is
  already-connected or an internal coefficient with no isolated measurement. "New observable"
  claim was selection-attributable (W3/c₂ failure mode), correctly demoted. One honest
  positive: D(2m+1) = (2^{2m}−2)ζ(2m−1) − ((2^{2m+1}−1)/2)ζ(2m+1). debug/door5_*

### Graduation tests (2026-06-01)

Doors 1 and 4 each got a second probe to graduate or break.

- **Door 1 recursion → WALL (clean, understood).** The odd-d F-theorem closed forms
  continue to S¹¹ — the flagship door stands. But the bonus recursion
  scalar−2·Dirac = c_d·Dirac_{d−2} does NOT graduate: c₁₁ has no rational form, and the
  original c₇=c₉=½ traced to an inconsistent Dirac normalization (a factor-of-2 switching
  at d≥7). Top-atom cancellation holds only at d=5, derived analytically. A convention
  artifact, now understood — the audit caught one a third time. debug/door1b_*
- **Door 4 inner algebra → PARTIAL DOOR (the session's deepest forward result).** The inner
  algebra STRUCTURE is mostly FORCED: among finite real *-algebras reproducing
  U(1)×SU(2)×SU(3), only ℂ⊕ℍ⊕M₃(ℂ) and ℂ⊕M₂(ℂ)⊕M₃(ℂ) survive — factor count, ℂ (n=1) and
  M₃ (n=3) all pinned by the Hopf tower. Residual non-uniqueness = one binary fork
  (ℍ vs M₂ at n=2). The conjectured handle (J_GV²=−1 selects ℍ) was tested — **Door 4c
  (J sign-table audit) closed it NEGATIVE:** over ℂ, ℍ and M₂(ℂ) are the same algebra
  (ℍ⊗ℂ≅M₂); the real-form distinction is an internal involution invisible to the combined
  J (combined J²/signs/KO-dim bit-identical for both). J_GV²=−1 ADMITS but does not FORCE
  ℍ — ℍ is a literature import (CCM), not GeoVac-forced. Honest ledger: factor count + ℂ
  + M₃ FORCED; n=2 real form ADMITTED-not-forced; generation count + Yukawa FREE. Written
  up in Paper 32 §VIII. debug/door4b_*, door4c_*

- **Door 4f (Falsifier A'' on Upgrade B, 2026-06-02) → PARTIAL DOOR FINAL (Outcome 3 NEUTRAL).**
  12 sub-tests on whether adopting Upgrade B (sphere-Lie-group axiom from Door 4e) introduces
  new downstream constraints on existing GeoVac observables: **0 new constraints**, **0
  contradictions**, 5 silent (T1 inner KO-dim, T2 γ_F, T3 Yukawa, T4 N_gen, T10 cross-rung
  CKM), 3 compatible (T5/T6/T7 Connes order-0/1/2), 2 identical to CCM (T8 spectral action
  gauge coefficient, T9 Higgs vacuum S²), 1 clean extensibility rule via Adams 1958 (T11),
  1 consistent strengthening of Bertrand × Hopf-tower reading (T12). Honest axiom-savings
  accounting: Upgrade B replaces 2 of CCM's 3 axioms (complex chiral fermion rep + 2N²=32
  dimension count) but the **CCM second-order condition stays in both formulations** because
  it constrains D_F's off-diagonal Yukawa structure, not A_F. Net axiom savings: **1 axiom**
  (not 2). The **Door 4 series is structurally complete at PARTIAL-DOOR FINAL.** Remaining
  open questions are PI/community-side judgment (adopt Upgrade B vs CCM-3-axiom) and the
  deeper wall (force N_gen / inner KO-dim from a packing principle — no known handle,
  deferred indefinitely). debug/door4f_*

- **Door 4e (Falsifier A' on DAS, 2026-06-02) → PARTIAL DOOR (CONFIRMED).**
  Literal Falsifier A' tested three natural AC tensor-product compatibility conditions for
  distinguishing ℍ from M₂(ℂ) at n=2: (C1) SU(2)-equivariance under Ad action — both pass
  (ℍ minimal, M₂(ℂ) maximal SU(2)-Ad-invariant *-subalgebra; no distinguishing); (C2) Hopf
  U(1)-equivariance under scalar phase — vacuous (every matrix commutes with e^{iθ}I);
  (C3) principal-bundle compatibility with Hopf U(1) fiber — both pass at machine precision.
  **No standard AC condition forces ℍ over M₂(ℂ).** Substantive content: CCM ℍ-selection
  requires 3 axioms (second-order condition + complex chiral fermion rep + dimension count);
  DAS requires 1 axiom (sphere = unit-norm subset of division algebra). Both are imports;
  DAS is **leaner**. FULL-DOOR upgrade path identified: adopt **Upgrade B (sphere-Lie-group
  axiom)** as foundational — "the inner algebra at rung n is the *-algebra whose unit group
  IS the rung-n sphere as a Lie group, when one exists; else fallback M_n(ℂ)." Reads the
  Hopf-rung sphere FULLY (Lie structure + topology) where the existing construction reads it
  PARTIALLY (gauge group only). Paper-level, sprint-reachable. Falsifier A'' (does Upgrade B
  introduce new downstream constraints?) is the next sharpest thread. debug/door4e_*

- **Door 4d (DAS, 2026-06-02) → PARTIAL DOOR (a NEW handle, structurally orthogonal to 4c).**
  Division-Algebra-of-the-Sphere criterion: the inner algebra at Hopf rung n is the natural
  associative real normed division algebra whose unit-norm sphere IS the rung's Hopf bundle
  total space S^(2n−1). Bit-exact verification at all three rungs: S¹ = unit ℂ (closure
  1.1e−16); S³ = unit ℍ = Sp(1) = SU(2) (closure 2.2e−16, SU(2)-iso 3.3e−16); S⁵ has NO
  associative-division-algebra realization (Hurwitz: only ℝ, ℂ, ℍ at dim 1, 2, 4 →
  fallback M₃(ℂ)). DAS **agrees** with Door 4b at n=1 (ℂ) and n=3 (M₃(ℂ)) and **closes**
  the n=2 fork by selecting ℍ — the direct (no-quotient) realization where U(ℍ)=Sp(1)=SU(2)
  matches S³ directly, while M₂(ℂ)'s SU(2) is reached only via the U(2)/U(1) unimodularity
  quotient. Structurally orthogonal to Door 4c: no J or sign-table data is invoked. Verdict
  PARTIAL-DOOR because Paper 32 §VIII.B's existing argument extracts only the gauge group
  from the Hopf tower, not the algebra realization; adoption as a forcing requires one
  structural theorem (Falsifier A', sprint-reachable, paper-level): the AC tensor-product
  H = H_GV ⊗ H_F transfers the rung-n sphere's division-algebra structure to H_F. If that
  closes, PARTIAL-DOOR upgrades to DOOR and the inner-algebra status becomes "fully
  forced except for N_gen, inner KO-dim, and Yukawa values." debug/door4d_*

### Door 1 — F-theorem closed forms in general odd d  `PROBED→DOOR`  ★ flagship
- **Forced:** the spectral-zeta F-coefficient is an exact closed form, bit-exact to
  Klebanov–Pufu–Safdi on S³ (61+ digits, PSLQ [8,2,−3]). Pushed to S⁵ the *same
  machinery already produced genuinely new closed forms* with ζ(5) terms (Paper 50 §7).
- **Crack:** Prop 7.4 *conjectures* the general-odd-d pattern and leaves it OPEN.
  This is the one place we have **direct empirical proof the engine runs forward** —
  it didn't re-check, it generated. Run it on S⁷, S⁹: does the ζ(3),ζ(5),ζ(7),…
  ladder continue, and does the scalar/Dirac dual-basis projection (Thm 3.4/3.5)
  survive in general odd d?
- **Probe:** compute ζ′_Δ(0) for scalar + Weyl-Dirac on S⁷, S⁹; PSLQ the forms;
  test the Prop 7.4 conjecture. Pure computation, sprint-scale.
- **Why a door not a wall:** the generativity is already demonstrated one dimension up.

### Door 2 — the BW wedge entropy was mis-filed as a failed area-law  `PROBED→PARTIAL`  ★ vindication
- **Forced:** S = 1.963·log(n_max) + 0.540, R² = 0.99992 over n_max ∈ [2,12].
- **Crack:** filed under §3 dead-ends as *"area-law REJECTED (R²=0.83)."* But the
  data it rejected area-law *in favor of* is a clean **Cardy–Calabrese boundary-log**:
  slope ≈ 2 = the dimension of the S² Hopf base, the same boundary the F-theorem
  (Door 1) lives on. We hunted area-law, didn't find it, logged the negative, and
  shelved the positive result sitting underneath it. This is exactly the failure mode
  we predicted — a negative tested against one hypothesis hiding structure — and the
  *first* sweep found it.
- **Probe:** (a) extend the panel; test whether the slope is *exactly* 2 and whether
  0.540 has a closed form — **AUDIT the "≈2" before claiming it** (fitted slope, not
  yet an identity); (b) connect to Paper 50: is wedge entropy the entanglement-side
  complement of the F-theorem F-coefficient on the same S² boundary?

### Door 3 — settle the Paper 35 π-criterion (full coverage)  `PROBED→THEOREM`  ★ sharpest unique-to-discreteness test
- **Forced:** "a GeoVac observable contains π iff its evaluation includes a continuous
  integration over a temporal/spectral parameter promoted from the discrete spectrum."
  This is the sharpest *discreteness-unique* claim in the corpus.
- **Crack:** it's stated as a falsifiable prediction but tested on a sample, not the
  whole transcendental inventory. Cross-check it against **every** π in Paper 34's 28
  projections + the case-exhaustion list. One counterexample = wall (informative);
  full coverage = the observation **graduates to a theorem.** Either outcome is a result.
- **Probe:** systematic sweep, Paper 34 inventory × the criterion. Sprint-scale.

### Door 4 — why does discreteness force the gauge group but free the Yukawas?  `PROBED→WALL+DOOR`  ★ deepest
- **Forced:** U(1)×SU(2)×SU(3) is *forced* (Bertrand × Hopf-tower truncating at n≤3).
  The Yukawa non-selection theorem *proves* the framework cannot force the matter
  sector at the current AC-extension level.
- **Crack:** the framework forces *half* the Standard Model and confesses it can't
  force the other half — and the **location of that line** is itself a structural fact
  we have described but never interrogated. If there's a "why" to where forced meets
  free, the why tells us whether the Yukawa wall is real or merely unexplored.
- **Probe:** characterize the structural reason (Hopf-tower truncation vs inner-factor
  calibration data). Analysis, not a single compute; medium-term. Low odds of a Yukawa
  derivation; high odds of sharpening the structural-skeleton boundary — itself a result.

### Door 5 — the ζ_{D²} / χ₋₄ integer-s ladder  `PROBED→WALL`  (narrow)
- **Forced:** ζ_{D²}(2) = π²−π⁴/12 (closed form); D_even−D_odd = 2^(s−1)(β(s)−β(s−2))
  for all integer s ≥ 2.
- **Crack:** the s-ladder is rigid. Has every integer s been swept for a spectral-zeta
  value that predicts an *unmeasured* QED / spectral observable, or only s = 2, 4?
- **Probe:** enumerate the ladder, map each value to a physical observable. Narrow compute.

### The overarching hunt — the discreteness-residue scan  `OPEN`  ★ the romantic program
Doors 1–5 are instances of one program. **Magic numbers (F-12) are the existence
proof** that integer rigidity from the skeleton can survive projection into a
*measurable.* The un-run question is whether there are others — a systematic scan,
across every forcing-certificate below, for *"where does an integer/rational
constraint of the bone poke through the transcendental skin into something you can
measure?"* That scan is the difference between the sober ending (atlas) and the
romantic one (a new door). It has never been run as a scan; it has only been
stumbled into.

---

## What the sweep revealed (meta)

Two findings about the corpus as a whole, independent of any single door:

1. **The engine has been run backwards.** Across ~165 certificates the FILED-AS tag
   is dominated by *consistency-check* (reproduces known physics) over *prediction*.
   That ratio **is** the quantitative measure of how much forward-run is left: every
   "✓ reproduces X exactly" is a forcing-certificate that was never asked for its
   neighbor. The catalogue is long because the backward-run was thorough; the doors
   are few because the forward-run has barely started.

2. **Negatives hide half-doors.** We predicted that results logged as negatives
   *against a specific hypothesis* would contain walked-past structure. Door 2 (BW
   entropy = Cardy–Calabrese, filed as "area-law rejected") confirms it on the first
   sweep. **Action:** re-read every `negative-with-structure` entry below asking not
   "what failed" but "what exact thing did the failure leave on the floor."

---

## A. Foundations / skeleton

### F-1. Graph Laplacian spectrum is exactly n²−1
- **Where:** Paper 0, Paper 7
- **Rigidity:** integer-spectrum (the root certificate)
- **Filed as:** prediction (foundational)
- **Crack:** this is the root rigidity every downstream certificate inherits. The
  question isn't "what does it force" (it forces everything) but "is the packing
  construction, followed rigorously, the *unique* generator of this spectrum?"
  (WH3 falsifier territory.)
- **Probe:** `OPEN` — re-derive the spectrum from the packing axiom alone, no
  physics input, and check for hidden choices.

### F-2. Eighteen symbolic proofs of S³ conformal equivalence
- **Where:** Paper 7; `tests/test_fock_projection.py`, `tests/test_fock_laplacian.py`
- **Rigidity:** bit-exact symbolic
- **Filed as:** ✓ check (topological integrity gate)
- **Crack:** all 18 are run as *checks*. Does any one of them, read as a forcing
  statement, have an unchecked neighbor? (e.g. the chordal-distance identity that
  makes 1/r the projection distortion — does the same identity force a second
  observable?)
- **Probe:** `OPEN` — audit the 18, tag each check vs. invertible.

---

## B. Spectral triple / math.OA

### F-3. F-theorem match to Klebanov–Pufu–Safdi, 61+ digits
- **Where:** Paper 50 §3 (scalar F_s, PSLQ relation [8,2,−3]; Dirac F_D symbolic-exact)
- **Rigidity:** bit-exact-match
- **Filed as:** ✓ check (reproduces KPS 2011)
- **Crack:** the **generativity proof of the whole catalogue.** Pushed to S⁵ it
  already produced NEW closed forms with ζ(5) (Paper 50 §7). Prop 7.4 states the
  dual-basis projection does NOT extend to S⁵ and *conjectures* the general-odd-d
  pattern — that conjecture is OPEN and is exactly a "what else does it force."
- **Probe:** `SCOPED` — Paper 50 §7 / Prop 7.4: settle general odd-d. This is the
  single clearest live forward-run.

### F-4. Master Mellin engine: M1/M2/M3 = 𝓜[Tr(D^k·e^{−tD²})], k∈{0,1,2}
- **Where:** Paper 18 §III.7, Paper 32 §VIII; `mellin_engine_domain_partition.md`
- **Rigidity:** proven-theorem (case-exhaustion)
- **Filed as:** organizing principle
- **Crack:** memory says the engine is "predictive *within* each domain." Has each
  domain's prediction actually been *extracted*, or only the classification stated?
  M1↔propinquity rates, M2↔heat-kernel, M3↔vertex-parity — is there an un-computed
  observable in any one domain the engine already determines?
- **Probe:** `OPEN` — per-domain: list what the engine predicts but we haven't computed.

### F-5. WH1 propinquity bound Λ ≤ C₃·γ → 0, C₃=1, rate 4/π
- **Where:** Paper 38; `wh1_proven.md`, `l2_quantitative_rate_4_over_pi.md`
- **Rigidity:** proven-theorem; 4/π = Vol(S²)/π² is the M1 signature
- **Filed as:** prediction (the "something new")
- **Crack:** 4/π is universal across all compact Lie groups (Paper 40). Universality
  is a forcing. Does the rate constant carry a *physical* prediction (a convergence
  speed that maps to an observable), or is it purely internal math?
- **Probe:** `OPEN`.

---

## C. QED / gauge

### F-6. Gauge group U(1)×SU(2)×SU(3) forced by Bertrand × Hopf-tower (n≤3)
- **Where:** `bertrand_sm_gauge_truncation.md`; Paper 32 §VIII.B; Sprint G4a
- **Rigidity:** proven-theorem (structural forcing)
- **Filed as:** semi-prediction
- **Crack:** **the single most principled candidate.** The framework forces *half*
  the Standard Model (gauge) and confesses (Yukawa non-selection theorem, 8 free
  params/gen) it cannot force the other half (matter). *Why does the line fall
  exactly there?* The location of the forced/free boundary is itself a structural
  fact, described but never interrogated.
- **Probe:** `SCOPED` — characterise *why* discreteness forces gauge and frees
  Yukawa. If there's a "why," the why is the door.

### F-7. ζ_{D²}(2) = π² − π⁴/12 (and the χ₋₄ identity)
- **Where:** Paper 28; D_even−D_odd = 2^(s−1)(β(s)−β(s−2))
- **Rigidity:** closed-form
- **Filed as:** new identity
- **Crack:** integer-s closed forms exist; have all integer s been swept, and does
  the s-ladder predict an unmeasured spectral-zeta value?
- **Probe:** `OPEN`.

### F-8. F₂ = 5√2/3 (graph-native anomalous moment, π-free)
- **Where:** Paper 33 (GN arc)
- **Rigidity:** π-free-certificate
- **Filed as:** prediction
- **Crack:** the π-free certificate sits next to the calibration 1/(4π) per loop.
  Does the certificate force the *next* coefficient's algebraic class?
- **Probe:** `OPEN`.

---

## D. Gravity

### F-9. Spectral action on S³ is two-term-exact (pure Einstein–Hilbert + Λ), uniquely
- **Where:** gravity arc (Paper 51); `fifth_asymmetry_gravity_termination.md`
- **Rigidity:** proven-theorem (Bernoulli identity; S³ unique, S⁵ has R²)
- **Filed as:** semi-prediction
- **Crack:** the *uniqueness* (only S³ terminates at 2 terms) is a forcing. Is there
  an observable that distinguishes "lives on the two-term manifold" from a generic
  higher-curvature theory?
- **Probe:** `OPEN`.

### F-10. Graviton kinetic spectrum is exactly (2k)²
- **Where:** Sprint G6-Full; `g6_full_graviton_integer_kinetic.md`
- **Rigidity:** integer-spectrum (bit-exact to 5 digits)
- **Filed as:** ✓ check
- **Crack:** an integer kinetic spectrum is rigid. Neighbor: does it force a
  degeneracy / mode-count that maps to a physical graviton DOF prediction?
- **Probe:** `OPEN`.

### F-11. BW wedge entropy S ~ 2·log(n_max) — NOT area law
- **Where:** Sprint BH-Phase0; `bh_phase0_entanglement_entropy_memo.md`
- **Rigidity:** scaling result (area-law rejected, R²=0.83)
- **Filed as:** **negative-with-structure**
- **Crack:** filed as a negative because it wasn't area-law. But 2·log(count) is a
  *dimension* — it counts the degeneracy of the lowest K_α shell. We hunted area-law,
  didn't find it, and may have walked past *what the number actually is.*
- **Probe:** `PROBED→PARTIAL` (Door 2, 2026-06-01) — half-door real but mis-described.
  Slope provably EXACTLY 2 (degeneracy n² forces it; 1.963 was a small-n artifact);
  constant has exact closed form C_∞ = coth(1) − log(2 sinh 1) = 0.4584…; F-theorem link
  RULED OUT (disjoint rings). Applied to Paper 50 §5. The negative correctly rejected
  area-law; the exact replacement is strictly better than the fit.

---

## E. Chemistry / quantum computing / nuclear

### F-12. Nuclear magic numbers from the packing skeleton
- **Where:** Paper 23
- **Rigidity:** integer (shell closures)
- **Filed as:** prediction
- **Crack:** **the existence proof that discreteness leaves a measurable residue** —
  an integer constraint from the skeleton that survives projection and that you can
  *measure.* The un-run question: are there OTHER residues? A systematic scan for
  "where does integer rigidity poke through into a measurable" has never been done.
- **Probe:** `SCOPED` — this is the romantic-ending hunt. Systematic residue scan.

### F-13. Pauli count = 11.11 × Q, universal across 28 molecules; isostructural invariance
- **Where:** Paper 14; CLAUDE.md §1.5
- **Rigidity:** exact-rational / integer (coefficient stable to ±0.1)
- **Filed as:** prediction
- **Crack:** universality + isostructural invariance (CO=N₂=1111) is a strong
  forcing. Does the invariance class predict anything beyond resource count — e.g.
  a conserved structural quantity shared by isostructural molecules?
- **Probe:** `OPEN`.

### F-14. Angular sparsity 1.44% depends only on l_max, not V(r)
- **Where:** Paper 22 (potential-independent sparsity theorem)
- **Rigidity:** proven-theorem
- **Filed as:** prediction (universal)
- **Crack:** potential-independence is a forcing across ALL potentials. Has it been
  used to predict sparsity for a potential class nobody has encoded yet?
- **Probe:** `OPEN`.

---

## F. Precision / entropy

### F-15. HO ground state has exactly zero entanglement entropy (rigidity)
- **Where:** Paper 27; `ep2b_ho_zero_entropy.md` ([H_HO,V_central]=0 by Moshinsky–Talmi)
- **Rigidity:** proven-theorem (exact zero)
- **Filed as:** prediction
- **Crack:** an exact zero is maximal rigidity. What perturbation is the *minimal*
  one that breaks it, and does that minimal breaker predict a structure?
- **Probe:** `OPEN`.

### F-16. Spatial Casimir on S³ = 1/240 (scalar), 17/480 (Dirac) — exact rational
- **Where:** Paper 35
- **Rigidity:** exact-rational (π enters ONLY via Matsubara temporal compactification)
- **Filed as:** prediction (falsifiable: π iff temporal integration)
- **Crack:** the falsifiable prediction (π appears iff continuous temporal/spectral
  integration) is itself the sharpest discreteness-unique claim in the corpus. Has
  it been tested against EVERY catalogued transcendental, or only a sample?
- **Probe:** `SCOPED` — cross-check Paper 34's full transcendental list against the
  Paper 35 criterion. A single counterexample is a wall; full coverage is a theorem.

---

## Sweep inventory (condensed)

Six arc-sweeps complete (2026-06-01); ~165 certificates total. The seeded entries
(F-1…F-16) above carry full crack columns. The remainder are preserved here in
one-line form — `result — rigidity — filed-as`. Cracks on these are surfaced
collectively via the **discreteness-residue scan** (Door, above): re-read each as a
forcing statement and ask for its measurable neighbor. `negative-with-structure`
rows are the priority re-read (see meta-finding 2).

### Gravity (32 found)
- ζ_unit(−k)=0 ∀k≥0 → two-term-exact spectral action (Bernoulli) — proven — check
- S(R,Λ)=φ(3/2)(ΛR)³ − ¼φ(1/2)(ΛR) + O(e^…) — closed-form — prediction
- u_crit²=φ(1/2)/(12φ(3/2))=1/6 (Gaussian) — exact-rational — check
- CH spinor spectrum |λ_n|=n+3/2, g_n=2(n+1)(n+2) — integer-spectrum — prediction
- heat-kernel factorization S³×S¹; ζ_{4D}(−k)=0 — proven — check
- 4D action βR³Λ⁴/4 − βRΛ²/8 (EH + Λ, exact) — closed-form — prediction
- a_k^Δ=2π²/k! (scalar Laplacian, infinite series); two-term is SPINOR-specific — proven — negative-with-structure
- S² Dirac |λ_n|=n+1, mult 4(n+1), a_1=−4π/3 — integer-spectrum — check
- S_BH=r_h²Λ²/3=AΛ²/(12π); G_N=3π/Λ² — closed-form — prediction
- G7/G4-2 factor-2 forced by Wald (not calibration) — proven — check
- s_min=−Λ/(12√6) dS vacuum density; R_crit·Λ=1/√6 across G1/G2/G5 — closed-form — prediction
- (1,1)=J0⊕J1⊕J2; S^(2) J-blind (Schur) — proven — negative-with-structure
- ∫(4u²−2)e^{−u²}=0 (background extremality) — exact-rational — check
- G_eff=6π/(φ(1)Λ²), Λ_cc=6φ(2)/φ(1)·Λ² — closed-form — prediction
- Δ_K^{Dirac,tip}(α)=−1/12(1/α−α) — bit-exact (5 sig) — prediction
- dΔ_K/dα|_{α=1}=+1/6 (replica derivative) — closed-form — check
- sector-wise Mellin map: tip↔φ(0), EH↔φ(1), Λ_cc↔φ(2) — proven — prediction

### QED / gauge (30 found)
- α³−Kα+1=0, K=π(42+π²/6−1/40) — exact-rational — prediction *(combination rule stays conjectural — §13.5)*
- B=42 / F=π²/6 / Δ=1/40 three spectral homes; no common generator (12 mechanisms) — proven — negative-with-structure
- Σ(n_ext=0)=0 self-energy structural zero — proven — check
- ζ_{D²}(s)=2^{2s−1}[λ(2s−2)−λ(2s)] (π^even only, T9) — proven — check
- D_even(4)−D_odd(4)=8(β(4)−G) (χ₋₄) — proven — check
- F₂/(α/2π)=1.084 on S³; Parker–Toms 0.4% — bit-exact — prediction
- S_min = 8π²ln2 − (2/3)π⁴ln2 − 3π²ζ(3) + (1/2)π⁴ζ(3) − (5/2)π²ζ(5) + π⁶/4 − (3/2)π⁴ − π⁸/96 (identified 2026-06-11; prior irreducibility = basis-coverage artifact) — closed-form — check
- pendant-edge Σ_scalar(GS)=2(n−1)/n — exact-rational — negative-with-structure
- SU(2) maximal-torus → U(1) exact; L₁ kinetic 1/(4N_c) universal — proven — check
- K π-prefactor = Vol(S²)/4 = Vol(S³)/Vol(S¹) (M1) — closed-form — check
- plaquette counts β₁=0,2,8 at n_max=2,3,4 — exact-rational — check
- Dirac D(s)=2(2^{s−2}−1)ζ(s−2)−½(2^s−1)ζ(s) — exact-rational — check

### math.OA / spectral triple (33 found)
- J²=−I, JD=+DJ, JOJ⁻¹=O bit-exact at n_max=1,2,3 (truthful CH) — bit-exact — prediction
- Spec(K_α^W)⊂ℤ; σ_{2π}(O)=O bit-exact (BW-α + BW-γ period closure) — proven — prediction
- Krein axioms J²=+I, D_L^×=D_L Frob=0; four BBB signs bit-exact at (m,n)=(4,6) — bit-exact — prediction
- C₃^{(2)}(n,n)=√(n/(n+2)); C₃^op=√(1−1/n_max) — closed-form — proven
- joint cb-norm = 2/(n_max+1) (saturation) — exact-rational — proven
- Riemannian-limit recovery at N_t=1 bit-exact; panel 2.0746/1.6101/1.3223 — bit-exact — prediction *(load-bearing falsifier)*
- universal rate 4/π rank-invariant (M1 signature) — closed-form — proven
- [D_L, a_s⊗a_t]=i[D_GV,a_s]⊗a_t (time Lipschitz-invisible) — bit-exact — prediction
- KO-dim 3 mod 8; order-zero fails 5–8% (finite-res artifact); offdiag breaks JD=+DJ (resid 2) — proven — negative-with-structure

### Foundations / exchange constants (31 found)
- κ=−1/16 = 1/Ω⁴(0) from Fock projection — exact-rational — prediction
- g_n=n² shell degeneracy; SO(6) Casimir ν(ν+4)/2 — exact-rational — prediction
- Fock coupling c²(n,l)=(1/16)[1−l(l+1)/(n(n+1))]; c²(4,3)=1/40=Δ — closed-form — prediction
- Bargmann π-free at N_max=5 (56 nodes, 165 edges, 0 irrationals) — π-free — prediction
- HO spectrum E_N=ℏω(N+3/2) exact, π-free (linear projection) — π-free — prediction
- Seeley–DeWitt a₀=a₁=√π, a₂=√π/8 on unit S³ — closed-form — prediction
- HO rigidity theorem (3D iso HO unique with SU(3) on S⁵) — proven — prediction
- Ollivier–Ricci κ_OR≡0 on S³ graph (κ lives in Fock interface, not discrete curvature) — proven — negative-with-structure
- case-exhaustion: every π ∈ {M1, M2, M3} — proven — prediction
- CH modular residual ε(t): all coeffs in √π·ℚ ⊕ π²·ℚ — closed-form — prediction

### Chemistry / QC / nuclear (22 found)
- O(Q^3.15) atomic Pauli scaling — proven — prediction
- N_Pauli=11.10×Q exact (38 molecules); 9.23×Q for d-blocks — exact-rational — prediction
- isostructural invariance NaH=KH=239; CO=N₂=F₂=1111 — proven — check
- split-region Legendre terminates exactly at L_max=2·l_max — proven — prediction
- Sturmian structural theorem H∝S → R-independent eigenvalues — proven — (guardrail)
- Fock projection rigidity (S³ unique to −Z/r; HO/Woods-Saxon excluded) — proven — prediction
- deuteron 592 Pauli/16q; He-4 712 Pauli/16q (12.25× Hilbert, 1.2× Pauli) — exact-integer — prediction
- HO shell closures 2,8,20,40,70,112 exact from graph counting — proven — check
- E*(He, n_max=1) = −729/256 Ha exact rational — exact-rational — prediction
- FCI basis invariance < 3×10⁻¹⁵ Ha — bit-exact — check
- Pauli ratio balanced/composed ~ O(B²): 2.63→4.77→7.45 — proven — prediction

### Precision / entropy (17 found)
- energy–entanglement decoupling (off-diag h1 <0.2% of S; all S from V_ee) — proven — structural
- S~Z^{−2.56} (n_max=4), Z^{−2.613} (n_max=3) — exact-rational — prediction
- ERI-density step function {0.424, 0.992} at θ=0 vs θ>0 — proven — structural
- core–valence I_cv<10⁻³ by Z=4 — integer-spectrum — check
- composed block entropy S_bond=0.303 nats R-independent — proven — structural
- one-body entanglement-inert S_kin/S_full~10⁻¹⁴ — proven — structural
- (1s,1s) Slater F⁰=5/8 exact (cusp hot-node) — exact-rational — check
- universal γ_∞≈1.96 (Richardson) — exact-rational — prediction
- area-law as pair counting: log A_n = 4 log n — exact-rational — structural
- correlation entropy PSLQ ring-orphan (∉ master Mellin, 12,312-form basis) — proven — negative-with-structure
- thermo S=k_B·S_info bit-exact (resid 4.8×10⁻¹⁴) — proven — check

---

*Source: six-arc Explore sweep, 2026-06-01. Full per-entry detail (WHERE / section /
eq) is recoverable from the papers and debug memos named in each domain. Next action:
pick a door and run its probe — Door 1 (F-theorem odd-d) and Door 2 (BW-entropy
re-read) are the cheapest forward-runs.*
