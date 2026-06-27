# Group 2 (Quantum chemistry) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only group2-specific scope + deltas + the branch watch-notes.

> **STATUS: FROZEN 2026-06-26 (PI-confirmed).** Pre-registration complete; the
> whole-group `/qa group2` run is authorized. Criteria are locked for this run —
> no goalpost-moving in either direction. Corrections are applied in place.

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
  - **Paper 11:** H$_2^+$ **0.0002%** energy (spectral Laguerre, $n_{\rm basis}=20$);
    the FD 1.01% is an *artifact* (must be flagged, not a competing result).
  - **Paper 12:** H$_2$ Neumann $V_{ee}$ recovers **92.4%** of $D_e$ vs 80.1% numerical;
    the **7.6%** gap is the cusp (a diagnosed limitation, stated as such). Algebraic
    recurrence = *exact*; the surviving transcendental seed is named.
  - **Paper 13:** He **0.022% raw** ($l_{\max}=7$) / **0.004% cusp-corrected**
    ($l_{\max}=4$), 2D variational, *properly variational upper bounds*; the 0.05%
    single-channel adiabatic is **non-variational** (lucky cancellation — must be
    flagged); graph-native CI **0.19%** at $n_{\max}=7$, zero parameters. (Note:
    CLAUDE.md §5 still drifts to "0.20% / $n_{\max}=9$" — papers win; the synthesis +
    §2 table use the paper numbers.)
  - **Paper 15:** H$_2$ **96.0%** of $D_e$ ($l_{\max}=6$, 61 channels, ~97% CBS);
    the adiabatic over-estimate (~11%) is a flagged artifact.
  - **Paper 17:** LiH composed $R_{\rm eq}$ **5.3%** ($l$-dependent PK); BeH$_2$ **11.7%**;
    H$_2$O **26%** (*uncoupled* five-block, $R_{\rm eq}=1.34$ bohr); LiH 4N $R_{\rm eq}
    \approx 64\%$ (**unbound** $D_e$, no PK — the equilibrium-without-PK control); 144×
    angular compression ($l_{\max}^2$ vs $l_{\max}^{11}$, cost of PK).
  - **Paper 19:** balanced coupled LiH **0.20%** *energy* ($n_{\max}=3$) with **structural
    $R_{\rm eq}$ drift (~8.8%)** — energy converges, geometry drifts; the 29% unbalanced
    figure is the negative control.
  - **FCI-atoms:** graph-native CI He 0.19% / Be 0.90% / Li 1.07% (zero-parameter, exact
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
