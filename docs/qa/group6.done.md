# Group 6 (precision / observations — the projection arc) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only group6-specific scope + deltas + the branch watch-notes.

> **STATUS: CERTIFIED ✅ 2026-07-04** (certifying FULL run PASS — 6th certified
> branch, after group3/group1/group2/group4/group5; the LAST paper group before the
> synthesis-layer cert. Run notes: this session + `debug/qa/group6_cert_seed_key.json`.)
> Originally FROZEN 2026-07-04 (PI-confirmed; frozen as drafted, tagging kept as
> a C3/C8 watch-note not a numbered criterion). Seventh pre-registered `/qa`
> target; the precision-AMO / projection-taxonomy branch — the last of the six
> paper groups before the synthesis-layer cert). Inherits criteria.md C1–C17 + the
> v4.62.1 run-shapes protocol (first cert = FULL run; delta-verification cycles
> between; final certifying FULL only after a clean delta; Sonnet-tiered
> code/citation reviewers with two seeds each; seeds COMMITTED onto the worktree
> branch per the run-8 hardening). Branch-defining risk = **the two-layer
> algebraic/observation split honesty (the branch's whole thesis) + transcendental-
> projection tagging correctness + §1.5 interpretive-rhetoric discipline (this is
> the most interpretive branch: time, entropy, "where physics enters")**. The
> K = π(B+F−Δ) hard prohibition has a *secondary* home here (Paper 34 catalogues
> it; Paper 2 in group5 is the primary home).

**Scope (non-trunk group6):** the **4 precision/projection papers** —
**Paper 26** (`paper_26_entanglement.tex`, ACTIVE — energy–entanglement
decoupling, basis-intrinsic ERI sparsity), **Paper 27**
(`paper_27_entropy_projection.tex`, **KEYSTONE** — entropy as a projection
artifact; one-body entanglement-inertness + HO/Coulomb rigidity theorems),
**Paper 34** (`paper_34_projection_taxonomy.tex`, ACTIVE — the **living
catalogue**: two layers, 28 named projections, three-axis tagging, the
Layer-2-presence-bound prediction), **Paper 35**
(`paper_35_time_as_projection.tex`, ACTIVE — π enters at temporal
compactification; the rest-mass vs observation/temporal-window split) — **+ the
group6 synthesis** (`papers/synthesis/group6_precision_observations_synthesis.tex`,
the C9 target). Trunk papers (0/1/7/18/24/32/38) canonical; in scope only where
restated (C7).

**Deterministic `--gate`:** `group6` (substring-matches
`group6_precision_observations/*.tex` **and** the
`group6_precision_observations_synthesis.tex` filename).

## Branch deltas (the only non-inherited content)

### Branch-defining criterion 1: the two-layer algebraic/observation split (Papers 34/35 — the branch thesis)

The entire branch rests on one structural claim: **Layer 1 is a bare
combinatorial graph (integer/rational/algebraic, π-free); every physical
variable, dimension, and transcendental enters through exactly one named Layer-2
projection.** The reviewers must hold the prose to what is actually verified:

1. **"Layer 1 is π-free" is a *verified-per-instance* claim, not a proven
   universal.** Paper 35's KG-spectrum result is a **SYMBOLIC verification over a
   finite panel** (n ∈ [1,50] × m² ∈ {0, 1, 1/4, 2} = 200 cases, ring
   ℚ[√d_i]) — not an all-masses/all-n theorem. Prose that reads "the graph is
   π-free" as a closed theorem overclaims the panel; MATERIAL.
2. **Paper 35's iff reading** — *"a GeoVac observable contains π iff its
   evaluation includes a continuous integration over a temporal/spectral
   parameter promoted from the discrete spectrum"* — is a **structural
   OBSERVATION / proposal**, tier-appropriate only if stated as such (it is the
   organizing reading of WH7, itself REGISTERED not proven). Any drift to
   "theorem"/"proven"/"we show that π appears iff…" is MATERIAL. The
   confirmatory list (Pauli counts, magic numbers, hydrogenic eigenvalues,
   Bargmann–Segal, graph-QED) is *evidence for*, not *proof of*, the iff.
3. **The Layer-2-presence bound (Paper 34) replaced a FALSIFIED depth-linear
   form.** The naive "residual scales linearly with chain depth" form was
   falsified on the catalogue itself (depth-3 chains at residual 0 *and* +286
   ppm; depth-4 at +2 ppm *and* +0.534%). The falsified form must appear only as
   a *superseded / withdrawn* reading — any locus that still asserts depth-linear
   scaling as live is a **zombie (MATERIAL)**; add to C16 on discovery.
4. **"Introduces no new computation"** (Papers 34 and 35 both say this): they are
   consolidation/verification papers. A claim in 34/35 that reads as a *new
   derivation* rather than a catalogue entry or a KG-1/2/3 verification is a
   framing defect.

### Branch-defining criterion 2: transcendental-projection tagging correctness (all four papers)

This branch **is** the transcendental catalogue, so a mis-tagged or anonymous
transcendental is a branch-defining MATERIAL class (the [[feedback_tag_transcendentals]]
discipline made a gating criterion here). The `claims-reviewer` must **enumerate
every transcendental** that appears in scope (π, π^{2k}, ζ(2k), odd ζ, Catalan G,
2π, Stefan–Boltzmann π²/90, 1/240, …) and verify each carries its correct
Paper-18 tier **and** its correct Paper-34 projection attribution inline:
- 2π / π per Matsubara mode ↔ the **observation/temporal-window projection**
  (Paper 35's headline; the split from the ring-preserving rest-mass projection
  is the paper's contribution — do not re-lump them);
- π via Hopf/Fock **volume measures**; π^{2k} via the **spectral action**;
  ζ(2k) via **even-zeta Dirichlet series**; odd ζ / Catalan G via the
  **half-integer-Hurwitz spinor lift**.
A transcendental attributed to the wrong projection, or appearing without a
tag on a load-bearing line, is MATERIAL (it breaks the catalogue's one job).
The Stefan–Boltzmann π²/90 and the Casimir 1/240 (exact rational — **no**
transcendental) are the canonical worked pair; check 1/240 is never described
as carrying π.

### Branch-defining criterion 3: the K-prohibition at its catalogue site (Paper 34; C5/C12 secondary home)

Paper 34 catalogues **K = π(B + F − Δ) ≈ 1/α = 137.036** as the "most extreme
multi-projection coincidence" (lines ~3554–3565, ~8553, ~9083–9102). The
`claims-reviewer` covering Paper 34 must **enumerate and quote EVERY K-sentence**
(not sample) and verify each carries the **Observation** tier — never
"conjecture", "derived", "predicted", "theorem", or a hedged equivalent. The
per-ingredient spectral homes (B, F, Δ) *are* derived and may be stated so; the
COMBINATION is the Observation. `check_k_label.py --gate group6` (C12) is the
deterministic backstop; the enumerate-every-K-sentence pass is the semantic gate.
Any violation is MATERIAL-LARGE (hard-prohibition touch), never a NIT.

### Branch-defining criterion 4: entropy-as-projection honesty (Paper 27 — the KEYSTONE)

Paper 27's load-bearing claims each carry a *qualifier* that is the actual
content — strip it and the claim becomes false:
1. **One-body entanglement-inertness requires NON-DEGENERACY.** "Single-particle
   vN entropy of a one-body fermionic ground state is identically zero" holds
   *only for a non-degenerate ground state* (S_kin/S_full ~ 10⁻¹⁴ for He at
   n_max = 2,3); open-shell and dissociated-bond states **escape the floor**.
   A locus stating the zero-entropy result without the non-degeneracy qualifier
   is MATERIAL.
2. **The n⁴ area law is a PAIR-COUNTING statement, not a one-particle law.**
   S_n = k log A_n with A_n = g_n² = (2n²)²; the factor of four is the two-body
   signature. Prose that credits the area law to one-particle grounds (the
   "earlier framework papers" reading the paper explicitly *corrects*) is a
   zombie — verify the correction is stated, not the superseded form.
3. **Entropy is an EMBEDDING-class projection observable**: zero on the
   combinatorial graph, finite only after the metric embedding of 1/r₁₂. Never
   "the graph has entropy."
4. **The HO rigidity** (Bargmann–Segal two-fermion closed-shell, *any* central
   potential → entropy ≡ 0 via Moshinsky–Talmi total-quanta conservation) is an
   INTERNAL THEOREM parallel to the Fock/BS rigidity of Paper 24 — symbolic
   backing required (C1/C2).
5. **The block-entropy exponent is γ_∞ ≈ 1.96** (Richardson on n_max = 2,3,4,5),
   *below* the second-order RS value 2, with the gap attributed to multi-shell
   aggregation. Watch for a stale earlier exponent (the raw-w_B / EP-2c 2.383
   form used dimensional w_B and inverts sign; the paper's live form is the
   **dimensionless** w̃_B/δ_B with γ_∞ ≈ 1.96) — a stale 2.383 or a "= 2"
   assertion is a number-drift defect.

### Branch-defining criterion 5: §1.5 interpretive-rhetoric discipline (Papers 35/27/34 — the most interpretive branch)

This branch makes the framework's strongest *interpretive* statements ("time,
not mass, is the projection axis"; "the observer's finite-time window is the
formal mechanism"; "entropy is a projection artifact"; "where physics enters the
graph"). Per §1.5 (dual-description framing, **no ontological priority**), these
must be presented as *structural readings of the dual description*, not as
ontological assertions that the discrete/compact side is more fundamental. The
honest cap is explicit in WH7 (input (iii): discrete-vs-continuous may be
empirically undecidable at every computable level). A sentence that asserts the
graph/compact/discrete description is *the* reality (rather than one of two dual
readings) is a §1.5 violation → MATERIAL (C3/C8). This is the branch's highest-
volume judgment-call surface; the claims-reviewers get the sharpest §1.5 prompt.

### Branch-defining criterion 6: trunk/WH status accuracy (C7)

Named trunk/WH dependencies this branch restates, with their current tier:
- **WH7 (time = observer-compactification)** — REGISTERED 2026-06-10, **not
  proven**; Paper 35 is its substrate. State registered, not established.
- **WH8 (Born measure = observation projection's exchange constant)** —
  REGISTERED 2026-07-01, Step-1 probe NEGATIVE-as-expected; Paper 34 §VIII
  carries the Born-rule record (external-input Class 1, tested-negative). Verify
  §VIII states tested-negative, not "resolved/derived."
- **Paper 18** (transcendental taxonomy / master Mellin engine) — the tier
  source; **Paper 24** (Bargmann–Segal, Coulomb/HO asymmetry) — the P27 HO
  rigidity parallel; **Paper 32** (spectral triple) and **κ = −1/16**
  (Observation, used in P26/27) — must be stated at current status, κ never
  "derived".

## C8 headline enumeration (per paper — the registry for C17 families as needed)

- **P26:** one-body off-diag H (κ=−1/16) = **10–39%** correlation energy but
  **<0.2%** entanglement entropy; V_ee off-diag = **100%** of entanglement;
  **S ~ Z^{−2.56}** across He-like ions; angular-momentum eigenbasis = the
  **unique** sparse-ERI basis (any rotation fills **42%→99%** step-function at
  identity; nonzero ERI count **Z-independent**); hub migration **1s→2s→2p**;
  core-valence MI **< 10⁻³ by Z=4**; composed per-block entanglement
  R-independent, core/bond entropy ratio **50×**; ties to **O(Q^2.5)** Pauli
  scaling. Tier = MEASURED (computational).
- **P27:** one-body non-degenerate GS single-particle entropy ≡ **0**
  (**S_kin/S_full ~ 10⁻¹⁴**, He n_max=2,3); area law **A_n = g_n² = (2n²)²**,
  factor-4 = two-body; cusp = single hot node on **(1s,1s)**; HO closed-shell
  any-central-potential entropy ≡ **0** (Moshinsky–Talmi); block entropy
  **S_B = A(w̃_B/δ_B)^γ**, **γ_∞ ≈ 1.96** (Richardson n_max=2–5, below RS 2).
  Tier = INTERNAL THEOREM (rigidity) + MEASURED (scaling).
- **P34:** two-layer decomposition; **28** named projections × three-axis
  tagging; **Layer-2-presence bound** |ε| ≤ max(ε_basis, max|L₂ input|)
  (replaces the FALSIFIED depth-linear form: depth-3 at 0 & +286 ppm, depth-4
  at +2 ppm & +0.534%); K = π(B+F−Δ) ≈ **137.036** catalogued as Observation;
  "introduces no new computation." Tier = framework consolidation +
  FALSIFIABLE PREDICTION (the bound).
- **P35:** KG spectrum on S³×ℝ π-free in ℚ[√d_i] for rational m²
  (**200 cases**, n∈[1,50] × m²∈{0,1,1/4,2}); π enters at temporal
  compactification, first π-eigenvalue **(n=0,k=1), ω²=4π²**; Casimir
  **E_Cas = 1/240** (exact rational, no transcendental); Stefan–Boltzmann
  **π²/90** via Matsubara high-T; rest-mass vs observation/temporal-window
  split. Tier = SYMBOLIC (panel) + structural Observation (the iff reading).

## Known logged gaps at freeze (not blockers; verify still-logged)

- **Matrix has zero group6 rows** (same first-cert state as group5) — the first
  cert's code dimension populates `docs/claim_test_matrix.md`; NO-TEST findings
  expected on this pass and are logged, not cert-blocking, unless a load-bearing
  headline has *no* backing at all.
- **Paper 26 has no dedicated `test_paper26_*.py`** in `tests/` (27/34/35 do:
  `test_paper27_entropy`, `test_paper34_projection_spot_checks{,_batch1-3}`,
  `test_paper35_predictions`). The code reviewer must locate P26's backing
  (it may live in entanglement/ERI-density modules under other names) or flag
  the headlines (Z^{−2.56}, 42%→99% fill, 50× ratio) as coverage gaps.
- The corpus-wide dangling-`debug/`-refs debt (§9 standing) touches
  P34/P35-adjacent files — C14 reports advisory, not blocking.

## Change log
- 2026-07-04 — **Certifying FULL run (whole-group) = PASS → group6 CERTIFIED ✅ (6th branch).**
  8-reviewer panel (3 Opus code + 2 Opus claims + Sonnet citation + Opus synthesis) + Opus
  completeness-critic, over the remediated corpus with a fresh 8-seed set (loci/phrasings
  distinct from the 1st-cert + delta runs). **Sensitivity 8/8 seeds caught** (S1 Marcolli
  wrong-venue, S2 count tautology, S3a/S3c positivity/π surrogates, S4 κ "derive analytically",
  S6 hedged K-forcing "necessarily…closing the anomaly", S8 one-particle area-law zombie,
  S9 WH7 "now-proven"); **specificity 6/6 controls clean**; every gating dimension exercised +
  calibrated. Deterministic: C10 5/5 compile ERRORS=0 + C11–C17 all PASS whole-target; 18/18
  check self-tests. **Zero verified MATERIAL** in any authoritative dimension. Convergence:
  the K-forcing seed was caught independently by Claims-B AND the completeness-critic; κ + area-law
  by Claims-A + cross-refs. Genuine (non-seed) findings all fix-on-sight NIT, remediated in-run:
  2 verified citation author-errors in the non-load-bearing Lorentzian survey cluster (franco_eckstein
  → sole-author N. Franco per arXiv:1210.6575; devastato → Farnsworth restored per arXiv:1710.04965);
  synthesis L361 core-valence wording synced to "≈2e-3 by Z=4, <1e-3 for Z≥5"; P35 "bare graph"
  → Fock/S³-projection (C6); P26 I_cv "<0.002"→"≈"; L9198 Riemannian-closure given the compact-KMS
  cross-ref. **Coverage closed:** the P26 core-valence-MI secondary gap (flagged by Code-A + synthesis)
  now backed by `test_paper26_core_valence_decoupling_thresholds` (recomputes Be/B/Li at n_max=2 vs
  Table II); matrix row added. **Honest ceiling:** ~14 of ~19 §V Roothaan-autopsy tables + the §V.C
  convention-exposure block + §III base-units un-arithmetic-checked (the critic verified 3 autopsies
  + ~8 catalogue rows, all internally consistent incl. cross-interval additivity — good base rate,
  not exhaustive); the ~443 dangling debug/ body-refs remain the known §9/C14 standing debt (advisory).
  190 group6 tests pass (2 slow-skip). Worktree removed, no seed leaked; Marcolli venue verified
  correct (J. Geom. Phys. 75) in the real corpus. Named follow-ups: HF-1 "= −a_e" approximate-equals
  NIT; §VI composition-cell-count deferred recomputation; the C14 debug-ref sweep (corpus-wide).
- 2026-07-04 — **Micro-delta over the Dirac-sign correction diff = CLEAN-DELTA → certifying FULL run EARNED.**
  The +17/480 correction touched production code (`thermal_tensor_triple.py`), so a micro-delta was run
  over just that diff. 2 reviewers (Sonnet code on the flipped value + its 2 tests; Opus claims on the 6
  flipped P35 loci; synthesis/citations untouched → not re-dispatched). **Calibration: claims caught its
  seed (a lone L927 locus flipped back to −17/480, via side-by-side internal-consistency) and confirmed
  the other 5 loci + boxed derivation mutually consistent at +17/480; code caught both seeds (tautology +
  sign-blind `abs`), monkeypatch-confirmed each passes on the retired −17/480, and verified 3 sound SIGNED
  pins prove +17/480 (each fails on −17/480).** Deterministic gates 7/7 whole-target. **Zero genuine MATERIAL
  in the diff** — every finding was a planted seed; the correction is derivation-correct, internally
  consistent, and genuinely backed. Worktree removed, no leak. **Certifying FULL run now earned (a clean
  delta is its precondition) — PI timing.**
- 2026-07-04 — **Delta-verification run over the remediation diff = DEFECTS → remediated in-run.**
  2 affected-dimension reviewers (Opus claims on the changed prose + Sonnet code on the 2 new
  test files; synthesis/citations untouched by the remediation → not re-dispatched). **Calibration:
  claims caught its seed (DELTA-S9 iff "proven theorem for all observables" over-scope); code
  caught both (DELTA-S2 tautology + DELTA-S3 float-cast π-free), and independently reproduced the
  genuine tests' numbers.** Deterministic gates 7/7 whole-target. **One genuine MATERIAL-LARGE in
  the remediation diff (the delta doing its job):** the 1st-cert "Dirac Casimir sign" fix was
  **WRONG-DIRECTION** — I had flipped all 6 P35 loci + the code + the test to −17/480 (matching the
  hardcoded `dirac_casimir_S3`), but the paper's own derivation gives ζ_{|D|}(−1)=−17/240 ⇒
  E=−½ζ=**+17/480** (POSITIVE, same sign class as the scalar +1/240; verified numerically to 40 dps,
  and papers are authoritative over hardcoded code per §1). Reverted to **+17/480** across
  paper (6 loci) + `geovac/thermal_tensor_triple.py` (value + flawed "fermions-oppose-bosons"
  docstring) + `test_thermal_tensor_triple.py` (2 asserts) + `test_paper35_kg_panel.py` (renamed);
  the C17 `dirac-casimir-s3-sign` family flipped (canonical now +17/480, retires the negative);
  matrix + this DoD synced. 1 NIT (sparsity test uses one rotation vs "any rotation"; P35 L791
  residual depth-criterion mention). **Next: a micro-delta over THIS correction (Dirac-sign +
  the code/test flip is a new diff touching production code), then the certifying FULL run.**
- 2026-07-04 — **1st-cert remediation APPLIED (Tier-1 + Tier-2, PI-directed).**
  **Tier-1 (fix-on-sight):** C10 compile (P26/P27 `dcolumn` + 3 `\multicolumn`
  `d`-cells; all 5 docs now ERRORS=0); C5 ×4 `conjecture`→`Observation` on the
  K/α combination (P34 L270/L274/L3564 + P35 L121); C8 core-valence "by Z=4"→"≈2e-3
  by Z=4, <1e-3 for Z≥5" (matches the table); P35 falsified depth-criterion →
  Layer-2-presence bound; entropy-measure mislabel (P26 "single-orbital ρ_A" →
  vN entropy of the normalized 1-RDM, the computed quantity); P27 `ep2b` E_full
  14.898→22.185, Richardson range 3,4,5→2,3,4,5, dangling `eq:pred2`→`eq:pred2_universal`;
  WH5→WH8 (Born, P34 ×2); Dirac Casimir sign harmonized ×6 (initially flipped to
  −17/480 to match the hardcoded code — **this was WRONG-DIRECTION; the delta run
  reverted it to the derivation-correct +17/480**, see the delta entry above);
  11.11·Q/28→11.10·Q/35 molecules ×3; P35
  "thirteen/fifteen"→past-tense/count-corrected. (The synthesis "208-case" and the
  entropy "occupation-vs-single-orbital" flags were PM-verified as non-defects /
  correct-on-inspection — 208 = KG-200 + TX-B-8 cumulative; kept.)
  **Tier-2 (judgment + new work):** P34 §III.29 **Lorentzian "literal identification
  at the Krein level / genuine Lorentzian extension" WITHDRAWN** and reconciled to the
  2026-06-09 P45 K⁺ descope + 2026-06-19 compact-boost/convention closure (3 blocks +
  the transfer-table row; Riemannian period-closure kept, signature-blind stated);
  P35 iff-promotion **scoped** to the Paper 32 §VIII case-exhaustion theorem's
  §III.1–§III.15 domain (the theorem IS an iff but by-exhaustion over 15 projections;
  WH7 registered beyond scope); **two new backing files** —
  `tests/test_paper26_entanglement.py` (4 tests: decoupling 10–39%/<0.2%/~100%, S~Z^{−2.56}
  =−2.563, 42.4%→99.2% + Z-independent 265) and `tests/test_paper35_kg_panel.py`
  (4 tests: genuine 200-case π-freeness, first-π-mode-is-temporal, Casimir 1/240 &
  +17/480) — closing the P26-lead-headline and P35-200-case-panel NO-TEST gaps;
  `claim_test_matrix.md` populated with 15 group6 rows; inline provenance tiers added
  to P26/P27; **C16 entry** `lorentzian-literal-identification-krein` + **C17 family**
  `dirac-casimir-s3-sign` added (both green + self-tested); 2 autopsy tables spot-checked
  (21cm + He⁺, both consistent). Verification: 5/5 compile ERRORS=0; 189 group6 tests
  pass (2 slow-skip); 7/7 deterministic gates PASS; 18/18 check self-tests. No `geovac/`
  production code touched. **Next: a delta-verification run over the remediation diff
  (PI timing) → certifying FULL run once the delta is clean.**
- 2026-07-04 — **1st cert (FULL run, whole-group) = FAIL (fully calibrated).** Panel
  8 agents (3 Opus code + 2 Opus claims + Sonnet citation + Opus synthesis + Opus
  completeness-critic). **Sensitivity 11/11 seeds caught, specificity 6/6 controls
  clean** — every gating dimension exercised + calibrated ⇒ verdict trustworthy.
  Genuine MATERIAL register (verified vs primary text): **C10** — P26 & P27 never
  compiled (missing `dcolumn` + 3 malformed `d`-cells in P27 `tab:gamma_surface`);
  fixed in-run, no numbers changed. **C5 ×4** — "conjecture" on the K/α combination
  at P34 L270/L274/L3564 + P35 L121 (C12 can't anchor without the explicit K formula;
  the claims panel + critic caught them). **C7 zombie** — P34 §III.29 pre-descope
  Lorentzian "literal identification / genuine Lorentzian extension / Krein four-witness
  closes bit-exactly" (L590/2208/2244/2263/2288/2412/2426/2443, dated 2026-05-16/17),
  withdrawn by the 2026-06-09 P45 K⁺ descope; add phrases to the C16 registry. **C1/C2**
  — P26 lead headline (energy-entanglement decoupling) + 42%→99% ERI-fill + core-valence
  MI + 50× + hub-migration have NO asserting test (P26 has no `test_paper26_*`);
  P35's 200-case π-free panel unbacked (only a float-cast false-positive stand-in);
  P35 Casimir 1/240 untested in-file (backed via P34 §III.14). **C3/C9** — P35 L386
  invokes the falsified projection-depth criterion as live; P35 iff-reading
  over-promoted to "spectral-triple theorem/corollary" (WH7 REGISTERED, verify vs P32
  §VIII). **C8-small** — P26 "core-valence MI < 10⁻³ by Z=4" contradicted by its own
  table (Be Z=4 = 0.002; body says Z≥5). **Systemic** — group6 has zero
  `claim_test_matrix.md` rows; no inline tiers in P26/P27. Honest ceiling: P34's ~24
  §V Roothaan-autopsy tables + ~68-row matches catalogue residuals un-enumerated (1
  autopsy spot-checked clean). Remediation = Tier-1 fix-on-sight + Tier-2 named
  (write the missing backing tests, reconcile §III.29 vs the descope, verify the
  iff-promotion, populate the matrix), then a delta run → certifying FULL run. Run
  notes: this session; seed key `debug/qa/group6_seed_key.json`.
- 2026-07-04 — **FROZEN (PI-confirmed).** Scope, six branch-defining watch-note
  criteria, per-paper C8 headline enumeration, and known-logged gaps set from a
  read of the four abstracts + synthesis spine + matrix/test inventory. PI froze
  as drafted; transcendental-projection tagging kept as a C3/C8 watch-note (not a
  numbered C18) per the group5 single-source pattern. Matrix rows and the seed
  plan follow in the first cert per the group4/group5 pattern. FULL first-cert
  run proceeding.
