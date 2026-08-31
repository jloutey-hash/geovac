<!-- CERT-STALENESS-BANNER -->
> ### ⚠ RE-CERTIFICATION OWED
> This record certifies the state as of **2026-08-30**. Since then **5 `.tex` changed**: group6_precision_observations_synthesis.tex, paper_26_entanglement.tex, paper_27_entropy_projection.tex, paper_34_projection_taxonomy.tex, paper_35_time_as_projection.tex.
>
> **The CERTIFIED verdict below is therefore historical, not current.** Do not cite it as present-tense status.
>
> Re-measure rather than trusting this banner — it is itself a snapshot and will go stale the same way:
> `python debug/qa/check_cert_staleness.py --detail`
<!-- /CERT-STALENESS-BANNER -->

# Group 6 (precision / observations — the projection arc) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only group6-specific scope + deltas + the branch watch-notes.

> **STATUS: NOT CERTIFIED — superseded 2026-08-22 (v5.0.0 retraction re-review)
> and again by the FULL certifying runs of 2026-08-28 (FAIL) and 2026-08-29
> (FAIL). See the change log and the run record before treating any part of
> this profile as a passed gate.**
>
> *Historical:* **CERTIFIED 2026-07-04** (certifying FULL run PASS — 6th certified
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
   falsified on the catalogue itself (depth-3 chains at residual 0 *and* -211
   ppm; depth-4 at +2 ppm *and* -0.534%). The falsified form must appear only as
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
5. **The block-entropy exponent is γ_∞ ≈ 1.93–1.96** (extrapolated from n_max = 2,3,4,5;
   the precise value is extrapolation-order dependent, the sub-2 conclusion is not),
   *below* the second-order RS value 2, with the gap attributed to multi-shell
   aggregation. Watch for a stale earlier exponent (the raw-w_B / EP-2c 2.383
   form used dimensional w_B and inverts sign; the paper's live form is the
   **dimensionless** w̃_B/δ_B with γ_∞ ≈ 1.93–1.96) — a stale 2.383 or a "= 2"
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
  **S ~ Z^{−2.56}** (measured −2.563) across He-like ions; angular-momentum
  eigenbasis = the **essentially unique** sparse-ERI basis (up to
  (l,m)-label-preserving rotations; any generic rotation fills **42%→100%**
  in a rapid but graded transition — intermediate densities are generator-
  and cutoff-dependent, the endpoints are not; nonzero ERI count
  **Z-independent**); hub migration **1s→2s→2p** (N/O/F values are
  multiplet-member-dependent — ground states **4/7/6**-fold degenerate
  for N/O/F (MEASURED; corrected 2026-08-30 from a stale "4/6/4 =
  C(4, n_p−2)" — the combinatorial formula does not reproduce the
  exact-rule degeneracies, and `tests/test_paper26_entanglement.py` pins
  {7:4, 8:7, 9:6}); range 0 to the attained ceilings
  no analytic ceiling — the former "ln 8 ≈ 2.08 N/F, ln 16 ≈ 2.77 O"
  ceilings are **RETRACTED** (artifacts of a diagonal-only single-orbital
  entropy routine that discarded spin coherence and so violated
  subadditivity; the maxima are MEASURED, not bounded). Core-valence MI
  **≈4×10⁻³ at Z=4** (Be, small-but-nonzero; `test_paper26_entanglement.py`
  pins 0.00365), decaying through **1.7×10⁻³ at Z=5** and vanishing
  **exactly for Z≥7** — NOT "exactly 0 for Z≥5", which the Z=5 measurement
  contradicts. Composed per-block entanglement R-independent, core/bond
  entropy ratio **40×**. Tier = MEASURED (computational); the core-closure
  statement is a derivation from a measured premise.
  *(Corrected 2026-08-30: this block previously carried the retracted
  ceilings, ≈2×10⁻³, "exactly 0 for Z≥5", a 50× ratio, and a tie to the
  retired O(Q^2.5) Pauli exponent. Criteria that assert withdrawn claims
  actively mis-direct reviewers.)*
  *(C8 block updated 2026-08-28; the pre-update block protected the
  withdrawn 99.2%→"step-function" characterization, the unqualified
  "unique", and a false <10⁻³ Be threshold.)*
- **P27:** one-body non-degenerate GS single-particle entropy ≡ **0**
  (**S_kin/S_full ~ 10⁻¹⁴**, He n_max=2,3 — the VALID one-body theorem); area law
  **A_n = g_n² = (2n²)²**, factor-4 = two-body; cusp = single hot node on **(1s,1s)**;
  **EP-2b RETRACTED (v5.0.0, moshinsky N_tot guard):** the "HO closed-shell
  any-central-potential entropy ≡ 0 (Moshinsky–Talmi)" claim + its mechanism + its
  INTERNAL THEOREM tier are **WITHDRAWN** — corrected S_full = **0.0671/0.0716/0.0833
  nats** (N_max=2/3/4), grows with basis, ground state NOT a single Slater determinant,
  ‖[H_HO,V]‖ = O(1) (0.74/0.63/0.67); block entropy **S_B = A(w̃_B/δ_B)^γ**, **γ_∞ ≈ 1.93–1.96**
  (Richardson n_max=2–5, below RS 2). Tier = INTERNAL THEOREM (one-body inertness — valid)
  + MEASURED (scaling); the two-fermion HO rigidity corollary is retracted.
  *(Criteria-currency correction 2026-08-24 to the PI-approved retraction — a live
  re-assertion of the ≡0 claim is now MATERIAL, not a headline to protect.)*
- **P34:** two-layer decomposition; **28** named projections × three-axis
  tagging; **Layer-2-presence bound** |ε| ≤ max(ε_basis, max|L₂ input|)
  (replaces the FALSIFIED depth-linear form: depth-3 at 0 & -211 ppm, depth-4
  at +2 ppm & -0.534%); K = π(B+F−Δ) ≈ **137.036** catalogued as Observation;
  "introduces no new computation." Tier = framework consolidation +
  FALSIFIABLE PREDICTION (the bound).

  **C8 EXTENSION — 2026-08-28 (PI-directed, BEFORE the freeze of this delta run).**
  Five remarks added to P34 after the 2026-07-04 certification (v5.1.5–v5.1.7) carry
  new headline numbers. They were unenumerated, and per the completeness-critic
  doctrine an unenumerated headline is an *unmeasured* criterion, not a passed one.
  Enumerated here so the delta can gate them:
  - `rem:minimal_presentation` — the **two-primitive** packaging (rational skeleton +
    projection family) with a **four-instance** operator-collapse table; the atomic
    one-electron identity **h₁ = k²(𝟙 − S/2) − Zk·diag(1/n)**, verified against an
    independent quadrature route at **6×10⁻¹⁰**; the boundary control — the
    two-electron tensor is NOT metric-generated (**39%** residual vs ~10⁻¹² for a
    metric-generated control). Tier: the four collapse instances individually
    proven/measured at their citations; **minimality is an OBSERVATION** with a stated
    falsifier. *Gate: the minimality claim must not read as proven.*
  - `rem:ee_partial_split` (base) — exact min/max split per multipole channel;
    **rank-two** head (labels×metric at L=0; A = ⟨r^L⟩ bandwidth **L+1**; B dense at
    L≥1); head supports **no correlation** (bit-exact vs dressed one-body, N=2,3);
    **g(k) = k·g(1)** with g(1) rational, **⟨1s²|g|1s²⟩ = 5k/8**;
    **⟨1s²|W⟩/⟨1s²|g⟩ = 3/5** exact; W compressible at basis-independent rank
    (**1 mHa @ rank 3–4**, **0.01 mHa @ rank 5–7**, n_s = 3–6); Frobenius tails after
    four terms **1.6/0.7/0.1/0.03%** for L = 0–3. Tier: exact (machine-verified) for
    the identities; MEASURED for the rank-flatness. *Gate: the lineage disclaimer
    (Beebe–Linderberg 1977 Cholesky; low-rank ERI structure is CLASSICAL and not
    claimed as new) must remain inline and unhedged.*
  - `rem:ee_partial_split` (sharpenings) — the **2n−1 exact-rank law** (pair densities
    span exactly 2n−1 dims ⇒ every radial-kernel matricization has rank ≤ 2n−1;
    integer dependence **−ρ₁₁ + 2ρ₁₂ − 3ρ₁₃ + 2ρ₂₂ ≡ 0** at n=3; nullity **(n−1)²**);
    the **closed-form NEGATIVE** — W(1) exactly rational (**⟨1s²|W|1s²⟩ = 3/8**) but
    active characteristic polynomials **irreducible over ℚ** (n_s=2 cubic
    **131072λ³ − 50688λ² − 25200λ − 675**), Löwdin spectrum **grows** with n_s;
    surviving skeleton form is the ODE **λf″(u) = w(u)f(u)** in u = 1/r.
    Tier: exact over ℚ; MEASURED for the spectral growth. *Gate: the negative must
    stay a negative — no drift toward "the eigenvectors have closed forms".*
  - `rem:ee_partial_split` (l>0 clause) — Gaunt-coupled s+p full-CI, L = 0,1,2 active:
    W-rank **3–4 → <1 mHa**, flat as the radial-pair space grows **10 → 28**.
    Tier: MEASURED.
  - `rem:ee_partial_split` (resource clause) — JW LCU 1-norm down **15–17%** at
    bit-identical energy; beats a density-fitting control by **13–17%** at matched
    accuracy; **λ_split/λ_full = 0.83, 0.85, 0.89, 0.92, 0.94** at n_s = 4,6,8,10,12,
    monotone toward unity. Tier: MEASURED. *Gate: the DEGRADATION is part of the
    claim — a reading that keeps the 15% and drops the decay is MATERIAL, and the
    clause must not read as a scaling result.*

  Backing: `tests/test_paper34_minimal_presentation.py` (5 legs),
  `tests/test_paper34_ee_split.py` (13 legs). Provenance:
  `debug/sprint_minimal_presentation_memo.md` §§7–11.
- **P35:** KG spectrum on S³×ℝ π-free in ℚ[√d_i] for rational m²
  (**200 cases**, n∈[1,50] × m²∈{0,1,1/4,2}); π enters at temporal
  compactification, first π-eigenvalue **(n=0,k=1), ω²=4π²**; Casimir
  **E_Cas = 1/240** (exact rational, no transcendental); Stefan–Boltzmann
  **π²/90** via Matsubara high-T; rest-mass vs observation/temporal-window
  split. Tier = SYMBOLIC (panel) + structural Observation (the iff reading).

## C8 EXTENSION -- 2026-08-29 (follow-on sprint, BEFORE the freeze of the FULL certifying run)

The 2026-08-29 follow-on sprint added new headline claims to P34, P27 and P26.
They are enumerated here so the certifying run gates them.  **This extension adds
criteria; it relaxes none.**  Disclosure: the criteria were edited by the PM
between the delta run and this certifying run (the completeness-critic's G15
class) -- the additions are listed in full so a reviewer can see exactly what
was added and when.

- **P34, deuteron polarizability channel (new, `sec:autopsy_d_hfs` +
  `sec:conv_d_polarizability` + `sec:prediction`).** Sourced TPE decomposition
  (Bonilla et al., arXiv:2508.18776): elastic **-41.50**, polarizability
  **+110.16**, single-proton **-33.14**, single-neutron **+8.95** kHz, total
  **+44.5(1.1) kHz = +135.9 ppm**; inelastic **+336.5 ppm** outweighs all
  elastic-class pieces **-200.6 ppm**, flipping the total positive.  Framework
  requires **+112.9 ppm**; difference **+23.0 ppm** attributed (as an
  INFERENCE from the H parallel, not an itemization) to missing higher-order
  QED; H 21cm comparator **+18.4 ppm**.  RETIRED: the "+44 ppm polarizability"
  (a kHz-read-as-ppm unit slip of the TOTAL) and the "~+200 ppm" entry.
  Prediction status: the D row is a *diagnosed missing-channel* exception, and
  a chain consuming the TPE is predicted to land in the tens-of-ppm class.
  Tier = MEASURED (framework side) against a cited Layer-2 evaluation.

- **P34, alkali cliff closed form (new, `eq:alkali_cliff`).**
  cliff = (Z/Z_eff^3)(n/nu)^3 F_rel(Z alpha), nu = sqrt(Ry/E_ion), F_rel =
  [gamma(2gamma-1)]^-1.  Predicts **2.9 / 4.4 / 22.3 / 41.5 / 56.9** against
  tabulated **2.8 / 4.0 / 22.0 / 43.1 / 63.0** (errors 2/10/1/4/10%).
  RETIRED: the **~Z^1.2 growth law** (misses Na by 126%; the fitted slope 1.16
  survives only as a data pin discriminating the retired ~Z^2.5).  Structural
  claim: nu is near-constant (**1.59-1.87**) while n runs 2..6; the
  quantum-defect factor spans 17x and is monotone while the charge factor
  spans <3x and is NOT; the K "jump" is a CR67-denominator artifact (the
  NEEDED density is smooth).  *Honest scope:* a diagnosis consuming nu from
  experiment, not a framework prediction.

- **P27, entropy-cusp co-location RESTORED to MEASURED (`sec:cusp`).**
  Edge-resolved leave-one-out decomposition on the pair-state graph:
  **42.5%** of S on the hottest V_ee edge (1s,1s)<->(1s,2s) at n_max=3,
  **40.0%** at n_max=4 (same edge); a comparably large non-reference edge
  costs **-1.7e-15** nats; edge ordering tracks V_ee ordering at Spearman
  **rho = +0.94** (+0.92 at n_max=4).  Explicitly NOT proportionality (the
  (2s,2s) edge carries 29.8% on a 6x smaller coupling).  Tier = MEASURED;
  the attribution is leave-one-edge-out, well-defined but not unique.

- **P26, molecular headline now BACKED (`eq:rindep`, `eq:ratio`).**
  S_bond = **0.3033139** nats R-independent on [0.5, 10] bohr with the
  mechanism pinned (composed electronic blocks **bit-identical** across R);
  S_bond/S_core = **49.6** at LiH R_eq; dissociation does NOT reach ln 2.
  (Previously a C8 headline with no test through three certifying runs.)

- **P34, citation layer (new `sec:citation_status`).**  13 references verified
  against primary text and added to the bibliography; "Eides 2024" corrected
  to Eides--Grotch--Shelyuto 2007 at 9 loci (no 2024 edition exists);
  two load-bearing attributions named as UNRESOLVED rather than fixed
  (the "Pachucki--Yerokhin 2010" D-HFS itemization, and the -(11/48)
  coefficient's textbook source).  A reviewer must check that these are
  presented as open, not as resolved.

**New deterministic gate in force for this run: C20** (inline-attribution
resolvability, ratchet semantics, `debug/qa/check_inline_attributions.py`,
baseline `debug/qa/inline_attribution_baseline.json`).  Registered in
`docs/qa/criteria.md` 2026-08-29; discrimination proven the same day.

## Known logged gaps at freeze (not blockers; verify still-logged)

- **Matrix has zero group6 rows** (same first-cert state as group5) — the first
  cert's code dimension populates `docs/claim_test_matrix.md`; NO-TEST findings
  expected on this pass and are logged, not cert-blocking, unless a load-bearing
  headline has *no* backing at all.
- **Paper 26 has no dedicated `test_paper26_*.py`** in `tests/` (27/34/35 do:
  `test_paper27_entropy`, `test_paper34_projection_spot_checks{,_batch1-3}`,
  `test_paper35_predictions`). The code reviewer must locate P26's backing
  (it may live in entanglement/ERI-density modules under other names) or flag
  the headlines (Z^{−2.56}, 42%→100% fill, 50× ratio) as coverage gaps.
- The corpus-wide dangling-`debug/`-refs debt (§9 standing) touches
  P34/P35-adjacent files — C14 reports advisory, not blocking.

## FULL certifying run 2026-08-28 (v5.1.8) — **FAIL**, remediated same day

Status: **NOT CERTIFIED.** Discharges the "Re-review OWED (2026-08-22, v5.0.0)"
item below by *running* the pass; the pass failed, was remediated, and the branch
now needs a fresh delta-verification run before any certifying attempt.

**Calibration.** 9 seeds planted, **2 VOID by PM seeding error** (planted on LaTeX
comment lines — see the new seed-placement rule in `docs/qa/seed_defects.md`). Of
the 7 valid: **7/7 caught**, four fire-tested by their finders. **7/7 controls
clean.** Per-dimension: code P26/P27/P34/P35, claims-A, claims-B **calibrated**;
citation **thin** (1 valid seed where its tier requires 2); synthesis
**UNCALIBRATED** — its only seed was void, so its clean verdict carries no measured
discriminating power and that dimension cannot be certified from this run.

**Cert-blocker.** D 1S HFS Bohr–Fermi baseline high by **496.5 ppm** (the
`g_atomic = 2μ_I/μ_N` convention needs division by **2I**; the driver divided by
**m_d/m_p = 1.99900750** and used `m_e/m_d` where the mass factor is `m_e/m_p` for
every nucleus). Two errors that nearly cancel; H and T unaffected (I=1/2 ⇒ both
factors unity), which is why the H sanity check hid it through three
certifications. Found by the **completeness-critic**, not the panel — §V.C has no
test file and its defects are column arithmetic, so no dimension reads it.
Corrected to 327.234993 / 327.315305; residual **+285.6 → −210.9 ppm**; propagated
to 26 loci including **certified group4**. The "+285.6 matches −286 ppm PY budget"
closure is **withdrawn**, with no replacement invented.

**Also fixed:** both Paper 26 LARGEs (§IV computed in the wrong basis; sparsity step
reproducing as neither claimed value), the incomplete-propagation family across
P27/P34/P35/synthesis/**this document**, three vacuous `sp.pi in free_symbols`
guards, one null test, and the missing C17 EP-2b family. Two upgrades applied
(core-valence at the 1e-15 noise floor; exponent agreement 0.11% not 0.6%).

**Verification after remediation:** C11/C13/C14/C16/C17/C18/C19 PASS on group6 *and*
group4; 6/6 papers compile three-pass clean with 0 undefined references; 219 tests
pass, 2 skipped.

**Open, not remediated (PI calls):** Li-7 HFS baseline (critic-flagged, **not**
PM-verified — depends on an assumed Z_eff, not a clean ratio test); the
Kennedy–Critchley–Dowker attribution in P35; ~50 informal "Author Year"
attributions in P34 §V; the fact that the corrected D itemization no longer closes;
the remaining 17 §V.C autopsies (never cross-anchored); P34 spot-check batches
(critic estimates 11–13 of ~35 functions lack discriminating power; 4 fixed).

Run notes: `debug/qa/group6_full_run_2026_08_28_notes.md`.

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
  =−2.563, 42.4%→100% + Z-independent 265 — **both corrected 2026-08-29 to 17.1%→100% and 107**; the same exact-rule correction later moved the P26 molecular cells (S_bond 0.303→0.330, still exactly R-independent; S_core 0.006→0.008; the 50×/49.6 ratio → 40.1) and the Minnesota adjudication was OVERTURNED by cert-3 item 3 (the ~31× is a zero-crossing artifact; production convention 0.48×), see debug/sprint_eri_evaluator_defects_memo.md) and `tests/test_paper35_kg_panel.py`
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
  — P26 lead headline (energy-entanglement decoupling) + 42%→100% ERI-fill + core-valence
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

---

## Re-review OWED (2026-08-22, v5.0.0) — group6

**Papers in this target changed after certification.** Logged here, at the
owning source, so the certified status is not read as covering text that
post-dates it.

- **Paper 27** (`paper_27_entropy_projection.tex`): EP-2b, one of the
  paper's **two headline "Rigidity Results"** and carried at **INTERNAL THEOREM**
  tier, is **RETRACTED**. Abstract, introduction, provenance-tier paragraph,
  taxonomy cross-reference, section lead-in, §sec:pred1 (full rewrite with the
  corrected `tab:ep2b`), and the conclusion all brought into line. Only ONE of
  the paper's two original predictions now survives.
- Pre-existing C10 defect fixed: missing `GeoVac_Paper19` bibitem.
- Backing tests rewritten (`test_paper27_entropy.py`, both EP-2b tests).
- **Note for the re-review:** the provenance-tier paragraph changed, so C3 is
  directly in scope for this target.

**Discharge condition:** a CLEAN DELTA over the changed loci (diff-scoped, per
the qa.md run-shape rule), which is also the standing precondition for the next
FULL certifying pass. The delta must re-test *these specific defects* rather
than trust this entry (qa.md hard rule, added the same day).

**Deterministic layer already re-run and GREEN on this target:** C10 (compiles,
with the aux-clean fix), C13, C14, C16, C17, C18, C19, C5/C12 (now corpus-wide
after the scope fix). What is owed is the LLM-judgment layer: claims, and
synthesis where the target has one.

**DISCHARGED — DEFECTS found + REMEDIATED, 2026-08-24 (v5.0.9).** Delta scope = the
changed loci (P27 `sec:pred1` EP-2b retraction + the group6 synthesis). Two dimensions
dispatched, both **exercised + calibrated**: claims on P27 (opus, 1/1 seed — a "stays
pinned near zero" re-assertion contradicting the growing 0.0671→0.0716→0.0833 table,
caught) and synthesis on group6 (opus, 1/1 seed — a factor-8 vs (2n²)²=4n⁴ slip, caught).

**Genuine MATERIAL found (the retraction re-review working as intended):** the group6
synthesis §"A second rigidity statement" was a **LIVE ZOMBIE** — it carried the *retracted*
claim forward in full at INTERNAL THEOREM tier ("entropy identically zero for any central
two-body potential", Moshinsky–Talmi N_tot conservation, single determinant |(0s)²⟩,
commutator ~10⁻¹⁶). The v5.0.0 retraction fixed P24/P27/P23/group4-synth but never
propagated to the group6 synthesis, and C16 `p24-entanglement-rigidity` was scoped
group3+group4 (+its file list omitted group6), so the deterministic gate never checked it —
the "deterministic GREEN" note above was NOT accurate for this target.

**Remediation (all verified):** (1) the synthesis paragraph rewritten as a retraction
(corrected 0.0671/0.0716/0.0833 nats, commutator O(1) 0.74/0.63/0.67, not a single Slater
determinant, tier withdrawn, no rigidity dual — parent HO rigidity theorem preserved);
(2) the minor P27 "S_HO ≪ S_Coulomb" wording tightened; (3) **C16 gate fixed** — scope +
file list widened to group6, two patterns added to match the synthesis phrasing, two
pre-existing false-positive patterns tightened (`= 0\b` matched "= 0.902"; `identically zero
for` matched an unrelated P23 kinetic-term line), discrimination proven BOTH directions
(fired on the live zombie: 2 hits; silent after remediation: PASS); (4) pre-existing C10
broken `\ref`s in paper_34 (×7) and paper_35 fixed on sight. Deterministic layer green
**now** (C10/C16/C19 PASS). A remediation-delta re-scan of the changed loci is CLEAN. Seed
key `debug/qa/group36_delta_seed_key.json`. Branch stays CERTIFIED; the standing precondition
for the next FULL certifying pass (a clean delta) is met after this remediation.

**FULL certifying pass = FAIL → REMEDIATED, 2026-08-24 (v5.0.9).** Fired after the clean delta
(PI-directed re-cert sequence group4→group6→group3). Seeded scratch copy of all 5 group6 papers
(synced from the working tree). Panel: 7 reviewers (claims-A P34 opus, claims-B P26+P27+P35 opus,
citations all-4 OPUS, synthesis opus, code-P27 + code-P34 sonnet, completeness-critic opus).
**Calibration: 100% sensitivity** — every planted seed caught (K "derivation" S5, P27 "does commute"
S8, synthesis "single determinant" S8, two same-cite-two-arXiv S1, code b1/b2 ×2); specificity clean.
Criteria-currency fix applied pre-freeze (§C8 P27 headline still asserted the retracted ≡0 at
INTERNAL THEOREM tier → corrected to the v5.0.0 values).

**The FULL pass FAILED — it found genuine material the delta + deterministic layer + the v5.0.0
retraction itself all missed (the retraction was INCOMPLETE):**
- **Three un-propagated retraction zombies** (v5.0.0 fixed §sec:pred1 but not these): synthesis
  §"What is robust" L465 ("the HO zero-entropy rigidity are INTERNAL THEOREMS"); Paper 27 L1069–77
  Equation-Verification bullet (marked the retracted S_HO=0 / (2,0,0,0) / noise-floor commutator
  "Verified" by now-nonexistent `..._is_single_determinant` tests); Paper 27 L1013–16 conclusion
  ("sharp contrast between Coulomb and HO … falsifiable"). All rewritten to the corrected facts.
- **Two genuine citation defects:** Krauth muonic-He-4 `arXiv:2102.05728` (both loci) → an unrelated
  football-biomechanics paper (dropped the wrong eprint, kept Nature 589,527); `arXiv:2508.18776`
  misattributed to "Patkos–Pachucki–Yerokhin" → resolves to Bonilla et al. 2025 (corrected).
- **Stale `claim_test_matrix.md` P27 row** (retracted claim as BACKED-SOUND/INTERNAL THEOREM, dead
  test name) → rewritten to the retraction + live tests.
- SMALLs fixed: P34 abstract Lamb sign (+0.534→−0.534), P34 conclusion stale −3.10% (LS-6a supersedes
  with −0.534%), P26 abstract "unique basis" overclaim (→ "essentially unique, up to (l,m)-preserving
  rotations"). P34 Layer-2-presence bound NO-TEST → logged in the matrix + **raised to PI** (hedged
  Prediction, not an overclaim). **CLOSED 2026-08-24 (PI-directed):** new
  `tests/test_paper34_l2_presence_bound.py` (6 tests) — the three exact L2-count=0 anchors COMPUTED
  (S³ Casimir 1/240, Stefan–Boltzmann π²/90, H 1S polarizability 9/2 a₀³ via a Dalgarno–Lewis
  symbolic solve → residual 0), the depth-linear form FALSIFIED, the bound checked on every catalogue
  row, L2-count-separates-where-depth-cannot, + a violation guard; cited inline in P34 §prediction;
  matrix row 118 NO-TEST → BACKED-SOUND. C10/C13/C19 green.

**Gate hardening this pass:** C16 gained a `zero-entropy rigidity` / `entropy-side rigidity|dual`
pattern; the too-broad bare `artifact` exempt marker was removed (a "projection artifact" 4 lines
away had false-exempted the L465 zombie under WINDOW=5); an `S_HO=0` generalization was tried and
REVERTED after it false-positived on a legit "…small but nonzero rather than S_HO=0…" comparison.
Discrimination re-proven; C16 group3/group6 + selftest PASS.

**Post-remediation state:** deterministic layer green (C10/C16/C18/C19 PASS); every confirmed
genuine defect remediated; completeness-critic's exhaustive zombie re-sweep found NO further live
zombie (only the L465 it re-saw in the pre-remediation scratch). **Honest ceilings (logged, not
blockers):** P34 §V.C autopsy region (~24 subsections, L3376–7293) spot-checked not line-by-line;
P34 C3 inline register tiers near-absent (adjective-tiering; 1 bracketed tag in ~10k lines);
Paper 27 §pred2's four (A,γ) framings to confirm mutually consistent; trunk restatements
(28-count, M1/M2/M3 map) out-of-path. **VERDICT: FAIL(genuine)→REMEDIATED.** Per the run-shapes
rule, a clean delta-verification over the remediated loci is the precondition for the FULL
certifying PASS; that re-verify is OWED (PI-gated). Seed key `debug/qa/group6_fullcert_seed_key.json`.

**RE-VERIFY DELTA = CLEAN → FULL certifying PASS EARNED → group6 RE-CERTIFIED ✅ (2026-08-24, v5.0.9).**
PI-directed immediately after the remediation. Diff-scoped pasted-hunks delta over the remediated loci,
2 Opus reviewers, both **calibrated**: claims-delta caught a planted re-introduced-zombie seed
(HUNK A "...is an INTERNAL THEOREM alongside it") and verified the real hunks B–E CLEAN (P27 Verified
bullet, P27 conclusion, P26 "unique", P34 Lamb ceiling all state the post-retraction facts);
citation-delta caught a planted wrong-arXiv seed (`2058.18776`) and CONFIRMED the real Krauth Nature
589,527 fix + the `2508.18776`/Bonilla attribution. Seed-leak clean (real corpus carries the correct
text). Deterministic layer green whole-target (C10/C16/C17/C18/C19/C11). Per the run-shapes rule a clean
delta over the remediated surface unlocks the PASS — **group6 is re-certified** (its retraction cleanup
is now complete: the v5.0.0 retraction's three missed loci are closed and C16 is hardened to catch the
class). Re-verify seed key `debug/qa/group6_reverify_seed_key.json`. Standing honest ceilings (P34 §V.C
autopsy line-by-line, P34 C3 tiers, Layer-2-bound NO-TEST → PI) carry forward, non-blocking.

<!-- Change log: 2026-08-29 -- FULL certifying run #2 = FAIL (calibrated: 9/9 dims, 14/15 first-pass seed catches all recovered on re-dispatch, 8/8 controls); full remediation applied same day (record: debug/qa/group6_full_run_2026_08_28_notes.md). Canonical corrections rippling into C8-adjacent facts: Minnesota contrast -0.55/+17.3 (~31x; the -0.81/+17.2 reproduced under no convention), Friar profile factor 3pi/8 (not 4/pi), K residual 8.8e-8 = Paper 2's self-consistency-root residual (plain sum 4.8e-7). Re-cert DELTA owed; plant a transcendental-tag-class seed there (critic G6). -->
