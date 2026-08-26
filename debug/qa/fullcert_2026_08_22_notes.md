# /qa FULL certifying run — 2026-08-22

**Target:** papers 58, 59, 60, 19 (group2) + 56, 24 (group3) + the group2 and
group3 syntheses.
**Shape:** FULL (all dimensions, whole-paper enumeration, completeness-critic).
**Seed key:** `debug/qa/fullcert_seed_key.json` — 14 seeds, 6 controls.
**Criteria:** `docs/qa/criteria.md` C1–C19, frozen before dispatch, plus the
per-paper profiles `docs/qa/paper_{58,59,60}.done.md` and the group2/group3 +
synthesis profiles. Nothing was relaxed for this run; three items were *added*
since the previous run (C9-currency ratified, C3-boundary adopted, C19 new) and
additions tighten rather than relax.

---

## Step 1 — deterministic layer (run whole-target, before any LLM dispatch)

| Gate | Result |
|:--|:--|
| C5/C12 K-label | PASS |
| C10 compiles + undefined refs | **found defects → remediated → PASS** (below) |
| C11 internal titles | PASS |
| C13 paper↔test refs | PASS |
| C14 paper↔file refs | PASS |
| C15 inline arXiv IDs | PASS |
| C16 retracted terms / zombies | PASS |
| C17 headline-number registry | PASS |
| C18 duration language | PASS |
| C19 eaten-escape corruption | PASS |

### C10 was not actually being measured

C10 had been certified by running pdflatex and reading its exit code. That is
not sufficient:

> **`pdflatex -halt-on-error` exits 0 on undefined references and citations.**
> They are warnings, not errors.

So a dangling `\ref{sec:does_not_exist}` — which renders as a bare `??` in the
PDF — could pass every "C10 green" ever reported. This run built the deterministic
implementation (`debug/qa/check_compiles.py`) and it immediately found live
defects, all pre-existing:

**In scope, remediated before dispatch:**
- **Paper 19 — three dangling `\ref`s**, all *cross-document* (`sec:relationship`,
  `eq:Vpk` → labels in Paper 17; `sec:conv_w1e_cross_domain_wall` → a Paper 34
  section that does not exist). LaTeX cannot resolve a label across documents, so
  these could never have bound. They carried load-bearing pointers (the PK
  definition; the "sixth instance of the multi-focal wall pattern" attribution),
  which under C10's own clause — *"not MATERIAL unless they break a load-bearing
  reference"* — makes them MATERIAL. Replaced with prose pointers.
- **Papers 19, 24 and the group3 synthesis — undefined `Note1` endnote citations**,
  each from a lone `\footnote` in a document whose class turns footnotes into
  endnote citations requiring bibliography machinery the papers do not run.
  Under C10's carve-out these are NIT (undefined *citations*, not load-bearing
  references); fixed on sight anyway by inlining the note content verbatim.

**A second flaw, in the new gate itself.** Its first version cleaned aux files
only *after* compiling, so a leftover `.aux` from an earlier build could make an
endnote resolve that would not resolve from a clean checkout — the target-scoped
run reported P19/P24 ok while the corpus sweep, which started clean, reported
`Note1` undefined in both. A gate whose verdict depends on what ran before it is
untrustworthy. Fixed to clean *before* compiling; the target was then re-measured
from a clean state, which is where the P19/P24 `Note1` defects actually surfaced.

**Post-remediation, all 8 in-scope documents PASS** the deterministic C10 from a
clean start (zero undefined references, zero undefined citations). Paper 56's
three warnings are font-substitution only (`T1/cmr/m/scit`), which is a local
font-availability matter, not a document defect.

### Out-of-scope corpus debt (logged, NOT part of this target's verdict)

The corpus-wide sweep found **12 further documents** carrying the same classes,
none in this target:

| Document | Undefined |
|:--|:--|
| paper_11_prolate_spheroidal | ref `fig:pes` |
| paper_13_hyperspherical | refs `sec:coupled_channel`, `sec:dboc`, `sec:liouville`; cite `paper15` |
| paper_15_level4_geometry | cites `loutey_paper11`, `loutey_paper13` |
| paper_17_composed_geometries | 5 refs; cite `Note1` |
| paper_18_exchange_constants | refs `sec:tensor_product_verification`, `thm:T9`; cite `paper7` |
| paper_22_angular_sparsity | 7 refs |
| paper_25_hopf_gauge_structure | cite `Paper16` |
| paper_27_entropy_projection | cite `GeoVac_Paper19` |
| paper_33_qed_selection_rules | ref `prop:reality` |
| paper_34_projection_taxonomy | 10 refs |
| paper_35_time_as_projection | ref `sec:projection_inventory` |
| group6_precision_observations_synthesis | cite `Note1` |

These are **standing debt for a separate sweep**, not findings against this
target. Several touch CERTIFIED branches, which means the branches were certified
under the exit-code interpretation of C10; the honest statement is that this
class was never measured before today, not that the earlier certifications were
mis-scored against the criteria as then implemented.

---

## Step 3 — seeds planted (worktree only)

14 seeds across every gating dimension, plus 6 known-good controls.

| Dimension | Seeds | Classes |
|:--|:--:|:--|
| code | 6 | S3 gutted tolerance ×2, S3 gutted band ×1, S2 tautological ×3 |
| claims | 2 | S4 overclaim + C3-boundary; S9 keystone status inversion |
| citation | 4 | S1 wrong metadata ×4 |
| synthesis | 2 | S8 zombie (non-existence from a bounded search); S9 scope creep |

Controls cover the results this session established or deliberately hedged: the
census counts and the 79.5/20.5 split, the eq:K0 tier upgrade, "apparently
novel", the exact gerade constant, and the deliberately *capped* 60-digit H2
figure (an underclaim — flagging it as wrong would be a false positive).

## Step 4 — panel dispatched

11 reviewers, all path-pinned to `C:/Users/jlout/Desktop/geovac-qa-fullcert` and
explicitly forbidden the real corpus:

- **code** (Opus, 1 seed each): one per paper — 58, 59, 60, 19, 56, 24
- **claims** (Opus): A = 58/59/60, B = 19/56/24
- **citation** (Sonnet, 2 seeds each): A = 58/59/60/19, B = 56/24 + both syntheses
- **synthesis** (Opus, 2 seeds): both syntheses, with the C9-currency leg exercised explicitly

Completeness-critic follows the panel.

---

*(Results, calibration scorecard, and verdict appended below as the panel returns.)*

---

# VERDICT: **FAIL**

Every gating dimension was **exercised** and **calibrated**, and multiple verified
MATERIAL defects were found. That is the definition of FAIL — *target not done* —
as distinct from INCONCLUSIVE, which the previous run returned. The panel is
trustworthy this time; the corpus is what did not pass.

## Per-dimension scorecard

| Dimension | Exercised | Calibrated | Clean |
|:--|:--:|:--:|:--:|
| Deterministic (C5, C10–C19; 10 gates, whole target) | YES | n/a | **NO** → remediated → YES |
| Claims P58/59/60 | YES | YES (1/1, 0 FP on 6 controls) | NO |
| Claims P19/56/24 | YES | YES (1/1) | NO |
| Citations P58/59/60/19 | YES | YES (2/2) | YES (seeds only) |
| Citations P56/24 + syntheses | YES (×2) | Sonnet 1/2 → **Opus 2/2** | NO |
| Synthesis (both) | YES | YES (2/2) | NO |
| Code P58 | YES | YES (1/1) | NO |
| Code P59 | YES | YES (1/1) | NO |
| Code P60 | YES | YES (1/1) | NO |
| Code P19 | YES | YES (1/1) | NO |
| Code P56 | YES | YES (1/1) | NO |
| Code P24 | YES | YES (1/1) | NO |

**Sensitivity 14/14** (after the protocol-mandated Opus re-dispatch of the one
de-calibrated agent). **Specificity 6/6 — zero false positives on controls.**

## Cert-blocking MATERIAL defects

**Wrong numbers / wrong equations**
1. **P58 §C8 headline #2 is wrong.** Cross-center magnitudes stated as intervals
   0.08–0.52 (overlap) and 0.19–0.80 Ha (one-body). The maxima reproduce exactly;
   the minima do not — the smallest overlap is **0.011** on the (2p0,2p0) pair,
   and the two stated minima came from *inconsistent subsets* of the census
   pairs. Confirmed by three independent computations (the reviewer's exact
   prolate machinery, its quadrature cross-check, and my own recomputation, which
   reproduced 0.0109 exactly). Restated on the verified maxima.
2. **P56 `thm:endo_rigidity` closed form is wrong.** n(n+1)(n+2)(3n+**5**)/24
   gives 2,11,35,85,175; the sum it claims to equal — and the values printed
   beside it — are 4,19,55,125,245. Correct factor is (3n+**13**). It survived
   because PS-4 (872 of the headline 5,864 residuals) had **no test at all**.
   Fixed, and `tests/test_paper56_endo_rigidity.py` added with a tripwire that
   rejects the retired polynomial.
3. **P24 `eq:edge-weight` does not describe the code that produced its
   certificate** — it writes a Clebsch–Gordan coefficient squared where the code
   computes a Wigner 3j squared, differing by (2l'+1) (verified: 1/2 vs 3/2 on
   the (0,0,0)→(1,1,0) edge). The prose two lines below already said "3j".

**False statements about our own work**
4. **P19's "Reproducibility note" was false.** It claimed "no number in this
   paper changes" after the sign fix, reconciling −8.055224 against −8.055 — but
   those are *different observables*: the former is the energy at the well's own
   minimum with a different offset, the latter at fixed R=3.015. From the banked
   data, the corrected values at R=3.015 are −7.9256/−8.0509/−8.0982 vs the
   published faithful −7.9297/−8.0548/−8.1020, so the n_max=3 error moves
   0.19% → 0.24%. Rewritten to state both conventions and what actually does not
   move (the geometry).
5. **P58 `tab:backing` asserted LIVE backing for the LiH 30-digit energy** that
   no test computes. The value is corroborated (reviewer reproduced 17 digits
   independently); the *label* was false → OWED.

**Citation misattributions (all load-bearing)**
6. **P56 presented a CONJECTURE as a PROOF in four places** and used it as the
   forcing argument for its comparison target: Eskandari–Murty–Nemoto 2025 do not
   prove Catalan G is not a period of MT(Q). Corrected everywhere; added to the
   C16 registry, **which immediately found a fifth locus in the group3 synthesis
   that no reviewer had flagged**.
7. **P56's Perez-Sanchez characterization contradicts both sources** —
   arXiv:2401.03705 derives Yang–Mills–*Higgs* (the opposite of the "without a
   Higgs" claim) and arXiv:2508.17338 is Yang–Mills, not the Standard Model.
8. **P56's `brown2017` bibitem is the wrong Brown paper** for the surjection
   claim it underwrites (ICM 2014 survey does not discuss U/U* at all).
9. **Bargmann 1961 mis-cited as CPAM 20, 1 (1967)** — that is Part II, a
   different paper — while the corpus's own group3 synthesis has it right.
10. Wybourne Ch. 16 / Iachello Ch. 4 pointers for the `thm:rigidity` converse are
    UNVERIFIABLE (Iachello Ch. 4 appears to be the wrong chapter);
    `follandstein1974` is the Heisenberg-group paper, not the sphere □_b
    spectrum; `neumann1878` resolves to no such work.

**Backing that does not prove what it is cited for**
11. **P59's PSLQ guard is precision-fragile and inverts** at the paper's own
    certified precision: the relative decoy rule passes at dps 19 and FAILS at
    dps 40 and 64. The sibling file states the correct methodology and explicitly
    refuses to do this. The *claim* survives (absolute height 19489 ≫ 10).
12. **P59's cited diagonal-A negative covers 2 of its 5 named generators** — and
    omits K_1, precisely the generator a closed form would need. The reviewer
    built the full ring and the negative survives; the backing did not establish
    it.
13. **P19's quadrature-invariance leg was a no-op** — it varied `n_grid_vne`,
    which the code explicitly ignores ("kept for API compatibility"), so the
    assertion could not fail. The framework's one live radial quadrature (the
    cross-block ERI trapezoid) is *not* invariant (~0.9% tilt spread). Rewritten
    as an honest no-op guard that says so.
14. **P60's abstract tagged the σ-law `[SYMBOLIC PROOF]`** where its own body
    says "the remainder is measured, not bounded". (This was a planted seed —
    but the reviewer also supplied a constructive fix: the leading order *is*
    derivable via a dual-Dirichlet argument, so the honest cure is to write the
    derivation, not to cut the claim.)

**Synthesis (C9 / C9-currency)**
15. group2 synthesis restored verbatim the reading Paper 58 exists to deny
    ("counted *pending an independence argument*" — the "ordinary work" reading;
    the paper says weight one is not the same as decidable, E1-value independence
    is open transcendence, and γ is not known even to be irrational). It also
    carried the stale 80% (decided 79.5% this run) and omitted §sec:qfd entirely.
16. group3 synthesis called Paper 31 "theorem-grade" 184 lines after saying it is
    "a structural reformulation, not a new theorem"; carried the retired
    "propinquity convergence rate" label at two loci — outside every registry
    file list, which is how it survived. New group3-scoped C16 entry **found a
    further live locus in Paper 18**.

## What the deterministic layer found — and a gate that was not measuring itself

**C10 was never actually being measured.** It was certified by reading pdflatex's
exit code, but `pdflatex -halt-on-error` **exits 0 on undefined references** —
they are warnings. A dangling `\ref` renders as `??` and passes silently. Built
`debug/qa/check_compiles.py`; it immediately found three dangling *cross-document*
`\ref`s in Paper 19 (labels belonging to Papers 17 and 34, unresolvable by
construction) plus undefined `Note1` endnotes in Papers 19/24 and the group3
synthesis. Then the gate itself proved non-deterministic — it cleaned aux files
only *after* compiling, so leftover state changed its verdict; fixed to clean
before, which is how the P19/P24 endnotes actually surfaced.

**Corpus-wide, 12 further documents carry the same classes** (Paper 34 has 10
dangling refs, Paper 22 has 7, Paper 17 has 5), several in already-CERTIFIED
branches. Logged as standing debt, out of scope for this target. The honest
statement is that this class was never measured before today — not that those
certifications were mis-scored against C10 as then implemented.

**New gate C19** (eaten-escape corruption) also earned its place: it caught a
live corrupted `\ref` in Paper 59 introduced by my own earlier remediation, and
bit again *while the gate was being written*.

## Two-way discipline

Applied upgrades: P60 compound-matrix [MEASURED] → [SYMBOLIC + MEASURED]; P60
abstract now names the exact gerade constant 2.555041…; P58 abstract now carries
the quadrature-free section it had omitted entirely.

Dissolved against primary text: the claim that P24 understates Paper 43 (the
reviewer preferred a `debug/` memo; Paper 43 itself states the N_t>1 leg as
*verified*, so P24's hedge is correct as written), and a P60 ratio NIT computed
from the wrong Gaussian convention.

## Honest ceiling

FAIL here means the panel was trustworthy and found real defects — not that the
corpus is now clean. The listed items are remediated, but this run did **not**
re-verify the remediated text with a fresh calibrated panel. Per the run-shapes
rule, the next step is a DELTA-verification pass over this diff; only a
subsequent FULL run on clean text can emit PASS. A substantial tail of NITs and
coverage gaps (notably P59's untested second-cusp Stokes block, P60's untested
v4.105 probe paragraph, and P56's ~49% non-discriminating residual panel) is
logged and not yet closed.

---

## Step 5 — calibration scorecard (as reviewers return)

| Agent | Tier | Seeds | Caught | Calibration |
|:--|:--|:--|:--:|:--|
| claims-A (P58/59/60) | Opus | K7 | **1/1** | CALIBRATED — 0 false positives on 6 controls |
| claims-B (P19/56/24) | Opus | K8 | **1/1** | CALIBRATED |
| citation-A (P58/59/60/19) | Sonnet | K9, K10 | **2/2** | CALIBRATED |
| citation-B (P56/24 + syntheses) | Sonnet | K11, K12 | **1/2** | **DE-CALIBRATED** → re-dispatched on Opus |
| synthesis (both) | Opus | K13, K14 | **2/2** | CALIBRATED |
| code-P24 | Opus | K6 | **1/1** | CALIBRATED |
| code-P58/59/60/19/56 | Opus | K1–K5 | pending | running |

**Specificity: zero false positives on controls so far.** The strongest single
signal is control M6 (the deliberately *capped* 60-digit H2 figure, an
underclaim): claims-A did not "correct" it — it explicitly praised the cap as
correct discipline and asked only that the achieved 6e-85 agreement be given
equal prominence. That is the control behaving exactly as designed.

### The one de-calibration, and what it teaches

citation-B (Sonnet) caught K12 (Bander–Itzykson triple-field drift) but missed
K11 (Bargmann 1961 CPAM 14,187 → 20,1,1967), and said why in its own report: it
verified the modern NCG/motives entries individually but treated the mid-century
classical entries as *"unambiguous, non-fabricable canonical works"* and did not
digit-audit them. That is precisely the class where quiet drift survives longest
— everyone recognises the title, nobody re-checks the digits — and the prompt had
warned about it explicitly. Re-dispatched on Opus with a mandate to verify every
bibliographic field individually and to report how many classical-era entries
were digit-audited.

## Genuine (seed-free) defects verified against primary text and remediated

**Citation dimension**
- **Paper 56 overstated a CONJECTURE as a PROOF, load-bearingly, in four places.**
  It asserted that Eskandari–Murty–Nemoto 2025 (arXiv:2510.20648) *prove* Catalan
  G is not a period of MT(Z)/MT(Q), and used that negative as the FORCING argument
  for adopting G_4 over Brown's G_MT(Z). The source establishes only the positive
  half; its abstract claims no negative result at all. Corrected at every locus:
  the level-4 choice is now *motivated* by where G is known to live, not *forced*
  by a proven exclusion. Added to the C16 registry — **which immediately found a
  fifth locus, in the group3 synthesis, that no reviewer had flagged.**
- Latrémolière ×2 (titles missing "Gromov–Hausdorff" — a documented recurring
  risk class for this project), Sugiura page range, "Brown 2017" for a 2014 paper.

**Claims dimension**
- P59 abstract stated the WRONG degeneration condition for K0(D) ("when the two
  orbital scales coincide" — that is rho→1, the genus-zero case the same abstract
  assigns two lines above; the K0 limit is rho→0, one scale decoupling).
- P19 "sub-linear in the block count (~B^1.14)" — B^1.14 is *super*-linear; the
  true claim is sub-QUADRATIC. **Introduced by my own previous remediation**
  (absent at HEAD) when the measured exponent replaced "O(B^2)" without fixing
  the adjective.
- P19 two contradictory W1e aggregate figures ("together at most 10–25%" vs
  "stack to ~35–45%"); 10.2 + 25.7 = 35.9, so only the aggregate was misstated.
- P56 rank remark still carried the pre-correction "3N+1 independent" against its
  own corrected 9N+1 rank — **also introduced by my previous remediation**.
- P56 used the monomorphism arrow ↪ at three loci for a map its own theorem
  proves is not injective.
- P58 "no fifth currency" after enumerating five instances → sixth.

**Synthesis dimension**
- group2 synthesis restored verbatim the reading Paper 58 exists to deny
  ("remain counted *pending an independence argument*" = the "ordinary work"
  reading; Paper 58: *"weight one is not the same as decidable"*, E1-value
  independence is open transcendence and gamma is not known even to be
  irrational). Also carried the stale 80% (decided 79.5% this run).
- group3 synthesis called Paper 31 "theorem-grade" 184 lines after saying it is
  "a structural reformulation, not a new theorem."
- group3 synthesis carried the retired "propinquity convergence rate" label at
  two loci — outside every existing registry file list, which is how it survived.
  Added a group3-scoped C16 entry; **it immediately found a further live locus in
  Paper 18**, also corrected.

**Code dimension (P24)**
- **Eq. `eq:edge-weight` does not describe the code that produced its
  certificate.** The equation writes a Clebsch–Gordan coefficient squared; the
  code computes a Wigner 3j squared — differing by (2l'+1), verified on the
  (0,0,0)→(1,1,0) edge (code 1/2 vs equation 3/2). The prose two lines below
  already says "3j", so the equation was the outlier. This is an equation-gate
  (§13.4a) defect: a reader reimplementing it verbatim builds a different graph.

## Two-way upgrades applied (prose weaker than backing)

- P60 compound-matrix result [MEASURED] → [SYMBOLIC + MEASURED] (an algebraic
  compound-matrix identity plus by-construction orthonormality; only the
  1×/8×/9× figures are measured).
- P60 abstract now names the exact gerade constant 2/(1+min j_0) = 2.555041…,
  basis- and R-independent, rather than only its small-n reading.
- P58 abstract now carries the quadrature-free section, which it omitted
  entirely — a substantive result invisible above the fold.

## Findings DISSOLVED against primary text (two-way discipline)

- **P24 "understates Paper 43"** (claims-B, upgrade-direction): the reviewer
  relied on a `debug/` memo claiming the HS-orthogonality is closed at general
  (n_max, N_t). Paper 43 — the owning paper — states it as a Corollary at N_t = 1
  and says the N_t > 1 extension is *"verified bit-exactly across the panel"*
  plus a structural argument. Paper 24's hedge ("a formal proof at general
  (n_max, N_t) is a named follow-on") is therefore CORRECT AS WRITTEN. The
  current-state rule decides: the owning paper outranks the memo. No edit.
- **P60 "10^3–10^4× understated"** (claims-A NIT): computed from the ratio-1.6
  Gaussian variant; the sentence is explicitly scoped to ratio 2, with the ratio
  dependence disclosed in the very next clause. No edit.

---

# DELTA-VERIFICATION RUN — 2026-08-22 — VERDICT: **DEFECTS**

Fired after the FULL run's FAIL→remediation, per the run-shapes rule. Scope =
the 42 loci changed **after** the FULL panel was dispatched, i.e. exactly the
text no calibrated reviewer had seen. Seed key `debug/qa/delta_seed_key.json`;
worktree destroyed, zero leakage verified against all five seeds.

## Calibration: 5/5 sensitivity, 4/4 specificity

| Agent | Seeds | Caught | False positives |
|:--|:--|:--:|:--:|
| claims (P58/59/60/19/56/24) | D1, D2 | **2/2** | 0 of 4 controls |
| synthesis (both) | D4 | **1/1** | 0 |
| citation (P19/58/56 + synth) | D3 | **1/1** | 0 |
| code (2 new + 1 modified + 2 gates) | D5 | **1/1** | 0 |

Every seed was a **regression of the exact defect the remediation had just
fixed**, so each agent was asked "did this repair hold?" rather than "is
anything wrong here?". Three catches went beyond the plant: the claims reviewer
caught D1 on an internal contradiction I had not anticipated (the inserted
clause referred to "the quoted lower end" when the corrected text quotes none),
and the code reviewer ran a mutation returning wrong census counts and showed it
passed all six tests in the gutted file.

## The delta's real value: my own remediation was incomplete in eleven places

**Stale siblings the FULL run's fixes did not reach**
1. Paper 58 still carried a stale `80\%` at l.433 — in the same paragraph as
   "only the $20.5\%$", summing to 100.5. (My grep escaping had missed it.)
2. The Eskandari–Murty–Nemoto correction had not propagated to two further loci,
   one **in the abstract**, inside a theorem's stated proof.
3. Paper 24 still said "the footnote at the head of Section…" after I had
   converted that footnote to an inline note.

**Corrections that were themselves wrong or lossy**
4. The W1e aggregate remained arithmetically self-contradicting: joint 35–45%
   implies a 55–65% residual, not the canonical 65–90%.
5. My rank-clause repair made the sentence vacuous — after 3N+1 → 9N+1 the "of
   which" clause stated the same number twice.
6. My group2-synthesis repair fixed the *direction* but dropped the
   **unconditional Liouville–Rosenlicht upgrade** Paper 58 records; the paper
   says in as many words that the position is "materially better than counted".

**Defects the FULL run had not reached at all**
7. **The "Brown establishes a surjection" attribution is in neither candidate
   Brown paper.** The cited ICM survey never mentions Connes, Marcolli, cosmic
   Galois, surjection or renormalisation; and Brown's actual cosmic-Galois paper
   explicitly declines the connection — *"It is not clear if it is at all
   related to the groups defined here."* The paper was flagging a "common
   framing slip" with a claim its own source disowns.
8. **Paper 19's `[MEASURED]` "bit-invariant to the radial quadrature" claim is
   false.** The test varied `n_grid_vne`, which `shibuya_wulfman` explicitly
   ignores (V_ne is analytic). The one live radial quadrature — the cross-block
   ERI trapezoid — moves the tilt 0.15% over 2000–8000 and 0.87% over 1000–8000.
   The FULL run repaired the test comment; the paper claim stood until now.
9. Paper 19's gradient-deficit headline reads 9.0% and measures **8.80%**; the
   test band (0.05–0.13, ±44%) could not discriminate it. Band tightened to
   0.084–0.093.
10. **Paper 56's `thm:inverse_limit` contradicts both the code and its own
    sibling theorem**: it states the sector set as `0 <= l < n`, while
    `geovac/pro_system.py` enumerates `l <= n` and `thm:endo_rigidity`'s
    (i+1)x(j+1) blocks *require* `l <= n`.
11. The class split is an n_max=2 figure (79.5%); at n_max=3 it is 80.0%. Tagged.

## Two gates I wrote were themselves defective

- **C19 missed most of its own defect class.** `\r` was caught only via a
  two-entry hardcoded tail list, so `\right`, `\rho`, `\rangle`, `\rule`,
  `\raggedright` all passed — forms far more common in this corpus's math than
  `\ref`. Its selftest contained no case for any of them, so it certified the
  covered half. Fixed: `\r` is now caught generically (CRLF is normalised before
  scanning, so a surviving lone CR is never legitimate), TAB handling is
  line-aware (indentation benign, mid-line suspicious — the previous lookbehind
  missed a TAB preceded by a space), and the selftest went 4→10 positives. The
  docstring's `\'e` claim was **withdrawn**: a swallowed `\'` leaves no control
  character, so it is structurally undetectable here and must not be advertised.
- **C10's cleanup list had a dead entry.** `".Notes.bib"` never matched, because
  revtex4-2 writes `<base>Notes.bib` with no dot; the gate left strays behind.
  Fixed and the strays removed.

**And one registry entry I added was a dead guard.** The new
`brown-surjection-attribution` C16 pattern was written with doubled escapes
(`\\\\emph` matches two literal backslashes; LaTeX has one) and fired on
**nothing** — C16 reported PASS and the class would have been believed guarded.
Caught only because I tested it, which I did only because a *pre-existing*
backspace-corrupted C17 regex had turned up earlier in this session. Two
incidents ⇒ `qa.md` now carries a hard **REGISTRY DISCRIMINATION RULE**: every
registry add/edit must prove it fires on the retired wording and stays silent on
the corrected one, with both results reported.

## Two-way upgrades earned

- **Paper 56's endomorphism closed form is a symbolic identity, not a panel
  result.** The reviewer derived it: `3n^2+19n+26 = (n+2)(3n+13)`. The paper's
  "Empirically across n_max in {1,2,3,4}" framing was dropped, a sympy proof of
  the full summation identity added to the test, and the citation raised to the
  symbolic tier. The backing now *proves* what the paper asserts.
- Paper 58's census-provenance paragraph states a general lemma (a same-centre
  Gaunt range is never empty), not a two-cutoff coincidence — the 176 removals
  at n_max=2 and 7,640 at n_max=3 follow at every cutoff.

## Verification

All nine deterministic gates PASS; all 8 target documents compile clean from a
cold start with zero undefined references; C19 selftest 10 pos / 5 neg, 0
failures; 40 passed on the census + endo-rigidity + topological-integrity
suites; the P19 slow leg passes with the tightened band; zero stray build
artifacts; worktree destroyed with zero seed leakage.

## Verdict

**DEFECTS** — not clean-delta. The panel was fully calibrated (5/5, 4/4) and
found eleven genuine incompletions plus three defective gates/registry entries,
all now remediated. Per the run-shapes rule a delta can never emit PASS, and
this one cannot emit CLEAN-DELTA either: too much was found.

**The honest read is that this delta earned its cost.** Four of the eleven
findings were defects *introduced or left behind by the FULL run's own
remediation* — which is precisely the failure mode a delta exists to catch, and
precisely why remediated text must not be assumed clean. A second delta over
*this* diff is the next step before any full certifying run.
