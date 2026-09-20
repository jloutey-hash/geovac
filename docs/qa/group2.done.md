<!-- CERT-STALENESS-BANNER -->
> ### ⚠ RE-CERTIFICATION OWED
> This record certifies the state as of **2026-06-28**. Since then **7 `.tex` changed** (plus 3 same-day, ambiguous): Paper_8_Bond_Sphere_Sturmian.tex, group2_quantum_chemistry_synthesis.tex, paper_11_prolate_spheroidal.tex, paper_12_algebraic_vee.tex, paper_17_composed_geometries.tex, paper_19_coupled_composition.tex, ….
>
> **The CERTIFIED verdict below is therefore historical, not current.** Do not cite it as present-tense status.
>
> Re-measure rather than trusting this banner — it is itself a snapshot and will go stale the same way:
> `python debug/qa/check_cert_staleness.py --detail`
<!-- /CERT-STALENESS-BANNER -->

# Group 2 (Quantum chemistry) — `/qa` profile

> **SCOPE EXPANSION 2026-09-13 (PI direction) — Papers 58, 59, 60 added.** These
> three papers are physically in `papers/group2_quantum_chemistry/` but post-dated
> the original group2 scope, so every prior group2 run silently excluded them. The
> scope in `debug/qa/qa_scopes.py` now resolves to **13 documents** (was 10). Their
> criteria are their own pre-registered per-paper definition-of-done files, which
> this group DoD incorporates by reference:
> - **Paper 58** (Abelian residue): `docs/qa/paper_58.done.md`
> - **Paper 59** (Elliptic Bessel moment): `docs/qa/paper_59.done.md` — note the
>   59/61 seam (Paper 61 is group3; they describe one object).
> - **Paper 60** (Sturmian secular / QC reading): `docs/qa/paper_60.done.md` —
>   NOT certified; carries a declared literal-registration debt (partly discharged
>   2026-09-13) and has never had a clean delta.
>
> **This is a baseline-establishing FULL run (PI direction), fired ahead of clean
> deltas on 58/59/60 deliberately** — like the trunk FULL runs #1–#8, its value is
> to re-measure the whole group's current surface, not to certify. The 2026-06-28
> CERTIFIED verdict below is historical and covers only the original 9 papers.

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only group2-specific scope + deltas + the branch watch-notes.

> **STATUS: CERTIFIED ✅ 2026-06-28 (v4.51.0, PI direction).** Five whole-group `/qa group2`
> runs (v4.50.0→v4.50.4); runs #4 and #5 returned PERFECT calibration (11/11 sens, 6/6 spec)
> with only thin secondary/provenance residuals (the thin-residual asymptote), all fixed.
> Every §C8 authoritative headline is verified-correct + soundly-backed; no headline was ever
> wrong. Certified under the 2026-06-28 secondary-number/provenance NIT carve-out
> (criteria.md §"Material vs nit"). Third certified branch (after group3, group1). See the
> Change log below + CHANGELOG v4.50.0–v4.51.0.

**Scope (non-trunk group2):** the **9 quantum-chemistry papers** —
**Paper 8** (Bond Sphere / Sturmian), **Paper 11** (prolate spheroidal H$_2^+$),
**Paper 12** (algebraic $V_{ee}$ / Neumann), **Paper 13** (hyperspherical He),
**Paper 15** (Level-4 molecule-frame geometry, H$_2$), **Paper 17** (composed
geometries, LiH/BeH$_2$/H$_2$O), **Paper 19** (coupled composition / balanced),
**FCI-atoms** (`paper_fci_atoms.tex`), **FCI-molecules** (`paper_fci_molecules.tex`)
— **+ the group2 quantum-chemistry synthesis**
(`papers/synthesis/group2_quantum_chemistry_synthesis.tex`, drafted this sprint;
this is the C9 target). Trunk/foundation papers (0, 1, 7) taken as already-canonical;
in scope only where a group2 paper restates them (C7).

**Deterministic `--gate`:** `group2` (path-substring match → `group2_quantum_chemistry/*.tex`).

## Branch deltas (the only non-inherited content)

### Branch-defining criterion: benchmarking-rule + guardrail-negative honesty

This is the group2 analog of group1's descope-accuracy criterion — the highest-risk
class for a *quantum-chemistry* branch, and it is a **sharpening of the shared C8
(headline honesty) + C3 (prose ≤ tier) + C5 (no §3 negative suppressed)**, not a new
numbered criterion. (Encoded as a watch-note, not a new number, to (a) respect the
criteria.md thin-profile / no-proliferation design and (b) avoid the existing
"C14" label collision — criteria.md's shared **C14 = paper↔file ref integrity**, while
group1's profile reused "C14" for its branch criterion; that pre-existing infra
inconsistency is flagged to the PI separately, not resolved here. If a second
chemistry-style branch needs this, promote a shared criterion.)

The reviewers (claims-reviewer, per paper, enumeration-forced) must verify ALL of:

1. **Benchmarking rule (§1.5).** Every accuracy headline *names the baseline it is
   measured against* and uses the **strongest relevant baseline** (cc-pVTZ or better
   for atoms; explicitly-correlated / near-exact for small molecules — Pekeris/Hylleraas
   He, Kołos–Wolniewicz H$_2$, exact H$_2^+$). Where a comparison is **unfavorable**, the
   paper says so honestly and identifies what the framework offers *instead* (O($V$)
   sparsity, baked-in angular selection rules, zero-parameter construction, scaling /
   structural insight). A favorable-looking number measured only against a weak baseline
   (STO-3G alone) without the strong-baseline context is MATERIAL.
2. **Guardrail negatives stay negative.** The two load-bearing NEGATIVE theorems are
   presented *as negatives* and not re-asserted anywhere as a working method:
   - **Sturmian Structural Theorem (Paper 8 / Papers 8–9 guardrail, §3.5):**
     single-center / shared-exponent / unified-basis molecular encodings give
     $H_{ij}\propto S_{ij}$ ⇒ **R-independent eigenvalues ⇒ no equilibrium**.
   - **FCI-M graph-concatenation (paper_fci_molecules, guardrail, §3.5):** LCAO graph
     concatenation gives R-independent kinetic energy ⇒ **monotonically attractive PES,
     no minimum.**
3. **PK / l_max ceilings honest.** The PK pseudopotential is the composed-accuracy
   bottleneck; **l_max = 2 is the optimal composed operating point** (l_max divergence is
   *structural*, not a convergence artifact — §3). No paper claims composed accuracy
   improves monotonically with l_max, and no §3 dead-end (the 6 PK modifications, the
   cusp treatments, the nested/Löwdin molecular encodings) is re-asserted as working.
4. **Positioning (§1.5).** The framework is a **research instrument demonstrating
   structural results**, NOT a production-chemistry replacement; the zero-parameter
   construction is stated as exactly that (no fitted/empirical parameter dressed as
   derived — C5).

### Per-criterion watch-notes

- **C8 (headline honesty), per-paper — the enumerated headlines + tiers.**
  - **Paper 11 — REWRITTEN 2026-09-19 (PI direction), mid-run:** H$_2^+$ energy is
    reproduced to **machine precision** (spectral Laguerre, $n_{\rm basis}=20$):
    measured $|E - E_{\rm ref}| = 3.6\times10^{-14}$~Ha against
    $E_{\rm ref} = -0.6026342144949$~Ha, i.e. at or below the precision to which
    that reference is conventionally quoted, so **no percentage is a valid
    criterion here** and none is stated. The FD 1.01% is an *artifact* (must be
    flagged, not a competing result).
    [retracted 2026-09-19: p11-h2plus-0002pct-retired] The criterion this line carried until 2026-09-19 was
    **0.0002%** ($=1.21\times10^{-6}$~Ha), which understates the method by ~7.6
    orders: even $n_{\rm basis}=5$ (3.5e-6 %) is 57x better, and nine
    quantity-x-reference combinations failed to reproduce it (nearest 4.64e-4 %,
    a coarse-grid PES fit, off by 2.3x). Also retired with it: R_eq 2.005 bohr /
    0.38% (a coarse-grid *fit* artifact; fine-grid fit gives 1.99726, +0.013%)
    and the "5000x accuracy improvement" derived as 1.01%/0.0002%.
    Registry: `p11_h2plus_err_ha`, `p11_h2plus_req_bohr`; C17 family
    `p11-h2plus-0002pct-retired`.
    **Grading note:** the CODE dimension of the 2026-09-19 run was graded against
    the OLD criterion and reported its falsification; any later dimension is
    graded against this one.
  - **Paper 12 — watch-note refreshed 2026-09-20 (/qa DELTA):** the live
    headline is **99.97% of $D_e$ / 0.053 mHa via the explicit-$r_{12}$ CI**
    (James--Coolidge, exact algebraic integrals, no quadrature; v5.14.10;
    registry `p12_r12_de_pct` / `p12_r12_err_mha`; backing
    `tests/test_paper12_r12.py`). *Coverage note:* the $(j_{\max},l_{\max})=(3,4)$,
    $n=416$ headline point is variationally **bracketed** by the backing test's
    $(3,2)$, $n=160$ rung ($>99.94\%$, monotone-nested, above exact) but is not
    directly regression-pinned at $n=416$ (the test docstring, which claimed it was
    "@slow-verified", was corrected 2026-09-20). The re-based CI **99.81% / 0.32 mHa**
    at the variational optimum $\alpha=1.40$, and **99.767% / 0.41 mHa** the fixed-$\alpha=1.0$ ladder
    endpoint (registry `p12_rebased_de_pct_aopt`, `p12_rebased_de_pct`; both real,
    neither superseding the other). The 92.4% below is the $\sigma$-only monomial
    ceiling and 99.1% the azimuthal-restored monomial value — both historical
    rungs, not the current criterion. **$\mu>0$ is quadrature-FREE** since
    v5.13.4 (closed-form $\{E_1,\gamma,\ln\}$ B-seeds, enforced by
    `test_B_table_uses_closed_form_not_quadrature`); the "spectral quadrature"
    clause below is superseded. **The 80.1% numerical comparator has NO backing
    test** and its only guard tolerates a ~31 mHa wrong-direction swing — raised
    to the PI 2026-09-19, unresolved.
    H$_2$ Neumann $V_{ee}$ recovers **92.4%** of $D_e$ vs 80.1% numerical.
    **AMENDED 2026-09-14:** the 92.4% is the **σ-only** ceiling, and the 7.6% gap is
    the absent $m \neq 0$ configurations, **not** the cusp — restoring them in the same
    basis reaches **99.1%**. The old watch-note ("the 7.6% gap is the cusp") is
    retired; grading a corrected paper against it would mark the correction as a
    defect. Algebraic recurrence = *exact* **for σ**; μ>0 uses spectral quadrature.
    The surviving transcendental seed is named.
  - **Paper 13:** He **0.022% raw** ($l_{\max}=7$) / **0.004% cusp-corrected**
    ($l_{\max}=4$), 2D variational, *properly variational upper bounds*; the 0.05%
    single-channel adiabatic is **non-variational** (lucky cancellation — must be
    flagged); graph-native CI **0.216%** at $n_{\max}=7$ (MEASURED 2026-09-19; the 0.19% this line carried was stale), zero parameters. (Note:
    CLAUDE.md §5 still drifts to "0.20% / $n_{\max}=9$" — papers win; the synthesis +
    §2 table use the paper numbers.)
  - **Paper 15:** H$_2$ **96.0%** of $D_e$ ($l_{\max}=6$, 61 channels, ~97% CBS);
    the adiabatic over-estimate (~11%) is a flagged artifact.
  - **Paper 17:** LiH composed $R_{\rm eq}$ **5.3%** ($l$-dependent PK); BeH$_2$ **11.7%**;
    H$_2$O **19.4%** (*uncoupled* five-block, $R_{\rm eq}=1.459$ bohr); LiH 4N $R_{\rm eq}
    \approx 64\%$ (**unbound** $D_e$, no PK — the equilibrium-without-PK control); 144×
    angular compression ($l_{\max}^2$ vs $l_{\max}^{11}$, cost of PK).
  - **Paper 19:** balanced coupled LiH **0.20%** *energy* ($n_{\max}=3$) with **structural
    $R_{\rm eq}$ drift (~8.8%)** — energy converges, geometry drifts; the 29% unbalanced
    figure is the negative control.
  - **FCI-atoms:** graph-native CI He 0.216% (2026-09-19; was 0.19%) / Be 0.71% / Li 1.03% (zero-parameter, exact
    rational Slater integrals); H$^-$ bound but over-binds 21% ($Z_c\approx1.84$ boundary).
  - **FCI-molecules:** the graph-concatenation **negative** (no minimum) — the headline
    *is* the negative result.
- **C4 (citations) — chemistry-method + reference-value surface.** Verify the reference
  energies / baselines actually say what is attributed: Pekeris/Hylleraas He ($-2.9037$ Ha),
  Kołos–Wolniewicz H$_2$ $D_e$, NIST atomic data, the cc-pVTZ/cc-pVDZ/STO-3G
  Gaussian-basis comparators, Phillips–Kleinman, Shibuya–Wulfman (Coulomb Sturmian),
  Neumann-expansion / prolate-spheroidal literature. Wrong reference value or
  misattributed method = MATERIAL.
- **C6 (discrete-vs-continuum precision) — the *mechanism* of the negatives.** The graph
  Laplacian is **dimensionless / scale-invariant / R-independent**; this R-independence is
  precisely *why* the Sturmian and graph-concatenation encodings fail (no equilibrium).
  The paper must state R-independence as the mechanism, and must not claim the *discrete
  graph* "has/produces" a continuum ($-(n^2-1)$ / Rydberg) spectrum as a bare graph
  property — it *converges* to it under the Fock projection (energy-shell constraint
  $p_0^2=-2E$, $\kappa=-1/16$).
- **C3 (prose ≤ tier).** Algebraic results (Gaunt/3j angular couplings, Neumann $V_{ee}$,
  split-region Legendre 3j-termination) = *exact / algebraic*; accuracy results =
  *matched / converges*, never *derived*. The natural-geometry construction is
  *zero-parameter* — assert exactly that and no more.
- **C7 (trunk-dependent status).** Where a group2 paper restates the Fock S³ equivalence,
  $\kappa=-1/16$, or the natural-geometry hierarchy, it states the current tier
  (Paper 7 = 18 symbolic proofs; $\kappa$ = Observation, not derived).
- **C9 — the new group2 synthesis** is the synthesis-faithfulness target (separate
  claims-reviewer dispatch): every claim traces to one of the 9 papers; the hierarchy
  spine + the two guardrail negatives + the honest accuracy ceilings are all faithful;
  no §3 dead-end re-asserted as working.

## First bite (PI-confirmed at FREEZE)

- **Whole-group** — all 9 papers + the new synthesis in one `/qa group2` invocation
  (PI direction, this sprint; same granularity as the certified group1 whole-group runs).
  Per `qa.md` step 4: deterministic layer (C10–C16) whole-group first; **code-reviewer
  1 paper/agent** (9 agents); **claims-reviewer chunked ~3–4 papers/agent** (≈3 chunks);
  **citation-reviewer chunked ~5–6 papers/agent** (≈2 chunks); **C9 synthesis** one
  dispatch; **completeness-critic** one agent. Per-chunk seeding: ≥1 seed catchable by
  each (dimension × chunk), planted in a non-first paper of each chunk, in the throwaway
  worktree only.

## Change log
- 2026-09-19 — **CRITERIA REWRITTEN MID-RUN (PI direction).** The `/qa group2 full`
  run of 2026-09-19 (PI-scoped to Papers 11/12/13) verified this file frozen at
  commit `4bd5a36`, ran the **CODE** dimension, and found the Paper-11 C8 headline
  **0.0002%** wrong by ~7.6 orders (measured 3.6e-14 Ha at $n_{\rm basis}=20$).
  The Paper-11 and Paper-12 watch-notes above were rewritten at PI direction
  **after** that dimension and **before** claims / citations / synthesis /
  completeness-critic, which had not run (session rate limit). **Consequence, stated
  so a later verdict stays auditable:** the CODE dimension was graded against the
  OLD criteria and reported their falsification; any subsequent dimension is graded
  against the NEW ones. The 2026-09-19 run's status is **INCONCLUSIVE** (four gating
  dimensions unexercised) and is not a certification under either set.
- 2026-09-13 — **BASELINE FULL run COMPLETE (v5.11.14–17), 4 batches, 13 docs.** A baseline
  re-measure (PI direction), **NOT a certification** — 58/59/60 are in by reference; clean
  deltas on them plus the owed in-paper extrapolation footnotes remain before any group2 cert.
  Per-document dimension results (**PASS** = reviewer returned clean, no defect; **fixed** =
  defect found → remediated this run):
  - **Paper 58 — code PASS, citations PASS, claims fixed** (1 SMALL: polyatomic FCI-vs-RHF
    label). The cleanest paper of the run: code (all 25 slow legs pass, mutation-guarded
    certificate + fire-tested decider, no LARGE) and citations (25+ externally grounded, no
    WRONG-ID/orphans) both **passed outright**.
  - **Paper 11 — code PASS** (216 tests). **Paper 12 — code PASS** (32 tests; citations fixed: KW literal).
  - **Paper 13 — code PASS** (citations fixed: Mitnik WRONG-ID, Madden misattribution).
  - **Paper 15 — code PASS** (citation fixed: Kolos digit; the run-#1 deleted-code regression is CLOSED).
  - **Paper 17 — code PASS** (claims fixed: §VIII.B adiabatic-vs-2D self-contradiction).
  - **FCI-molecules — code PASS** (guardrail-negative clean, 48 tests; M3/M4 overclaim softenings fixed).
  - **Paper 8 — code fixed** (LARGE: Sturmian-guardrail scope over-strong at summary surfaces
    + a NaN'd backing test; PI-approved scope propagation) + citations fixed (KW literal).
  - **Paper 19 — fixed** (MATERIAL: pair-diagonal zombie in prose vs its own re-measured tables).
  - **FCI-atoms — code fixed** (LARGE: 3 RED Li/Be backing tests + stale Table I; the 0.027-Ha
    Be shift independently cross-checked before re-pinning; propagation completed in Batch 4
    after a locus-incomplete first pass).
  - **group2 synthesis (C9) — fixed** (LARGE: carried the retired FCI-atoms energies; this is
    what surfaced the incomplete Batch-3 propagation).
  - **completeness-critic:** GAPS → all closed. Deterministic layer **8/8 PASS**
    (C10 compile 13 papers, C13/C14/C16/C17/C19/C21/C22).
  All findings remediated; every paper compiles and now passes its reviewed dimensions
  post-remediation. Full chronicle: CHANGELOG v5.11.14–17.
- 2026-06-26 — **DRAFTED** by PM for PI freeze (fourth pre-registered `/qa` target;
  first quantum-chemistry branch). Inherits criteria.md C1–C16. Branch-defining risk =
  benchmarking-rule + guardrail-negative honesty (encoded as a C8/C3/C5 sharpening, not
  a new number — avoids the criteria.md "C14" label collision, which is flagged to PI
  separately). C9 = the new group2 synthesis (drafted + audited vs primary text this
  sprint: H$_2$O Table-I label corrected to *uncoupled*; He/4N/compression numbers
  confirmed against Papers 13/17). First bite = whole-group (PI direction).
- 2026-06-26 — **FROZEN + whole-group run #1 = FAIL** (PI-confirmed freeze). Panel
  FULLY CALIBRATED (sensitivity **10/10**, specificity clean). Deterministic
  C10–C16 PASS (fixed a genuine pre-existing paper_15 compile bug; 6 missing figures
  flagged). Extensive verified MATERIAL defects, several STRUCTURAL — most notably
  **headlines backed by code deleted in v2.7.0** (Paper 11 spectral solver; Paper 15
  spectral 16×/20×/269×) and **retired `MolecularLatticeIndex`** (FCI-atoms LiH;
  entire FCI-molecules pipeline), plus stale/non-reproducible numbers (P17 BeH₂
  11.7%→19.7%; P19 balanced-LiH 0.20% non-reproducible post V_NN fix; P13 H₂
  rovib +10.5%), a confirmed l_max-divergence §3-suppression + P17 self-contradiction
  (1068 vs 1173), Paper 12 "no integrals" contradicted by `scipy.quad` B_l, and ~5
  genuine citation splices. "Code decides" (PI): the recompute CONFIRMS the paper He
  numbers (0.022%/0.004%/0.19%) and the §2 fix; §5 + claims_register still stale.
  Full inventory: `debug/qa/group2_whole_run_notes.md`, seed key
  `debug/qa/group2_seed_key.json`. **Disposition: TIER-1 structural items need PI
  strategic direction (restore deleted code vs descope/reframe headlines) before
  remediation.**
- 2026-06-27/28 — **Re-cert arc #2–#5 → CERTIFIED ✅ (2026-06-28, v4.51.0).** Five whole-group
  runs total. The FAILs peeled monotonically thinning layers: #1 restorations (v2.7.0-dropped
  code) → #2 framing/§3-TC-reassertion + P8 theorem re-derivation → #3 cross-paper second-locus
  propagation (~25 loci) → #4 two thin framing-syncs → #5 two thin secondary/provenance residuals.
  **Runs #4 and #5 both returned PERFECT calibration (sensitivity 11/11, specificity 6/6)** with
  only thin DIFFERENT secondary/provenance residuals (never a §C8 headline) — the thin-residual
  asymptote. **NO group2 paper headline was ever wrong** across the arc (the run-#1 FAIL resolved
  into restorations + one regression fix + one improvement; the genuine corrections were a single
  false claim [FCI-A "Li beats HF", run #3] + framing/secondary-number/provenance fixes). Every
  §C8 authoritative headline is verified-correct + soundly-backed, twice over. **Certified per the
  2026-06-28 secondary-number/provenance NIT carve-out** (criteria.md §"Material vs nit"): the
  remaining long-tail residuals are fix-on-sight NITs, not cert blockers. Per-run detail:
  CHANGELOG v4.50.0–v4.51.0 + `debug/qa/group2_{whole,recert,recert3,recert4,recert5}_run_notes.md`.
