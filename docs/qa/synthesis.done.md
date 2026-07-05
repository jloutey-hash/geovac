# Synthesis layer (the field guide + cross-group narrative) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only synthesis-specific scope + deltas + watch-notes.

> **STATUS: CERTIFIED ✅ 2026-07-04** (7th and FINAL certified target — the whole-corpus
> QA sweep is complete: trunk + group1–6 + synthesis). Path: 1st FULL run INCONCLUSIVE
> (PM worktree setup error) → field-guide re-run CLEAN → citation delta DEFECTS (3 LARGE
> group1-bib wrong-cites the branch cert missed) → remediated (v4.67.1/.2) → citation
> **re-delta CLEAN + calibrated (2/2 seeds, specificity clean, zero MATERIAL; 1 NIT
> minguzzi-title fixed on sight)**. Per PI direction a clean citation delta = the
> certifying pass; all gating dims (claims re-run, citations re-delta, synthesis-
> faithfulness, deterministic) calibrated + clean; code/test-backing WAIVED (inherited).
> **Honest ceiling:** the group1 bibliography needed two remediation cycles — a full
> external-citation audit of all group syntheses remains a standing candidate (PI-deferred);
> groups 2–6 synthesis bibs got light spot-checks this cert (they passed their branch certs).
> Run notes: this session + `debug/qa/synthesis_{seed,delta,redelta}_seed_key.json` (gitignored).
> Originally FROZEN 2026-07-04 (PI-confirmed; frozen as drafted. **Code/test-backing
> dimension WAIVED by PI direction** — the synthesis layer inherits its backing from
> the already-certified group-synthesis papers, so C1/C2 is **out of scope** for this
> target, NOT an unexercised gating dimension. The gating dimensions are therefore
> **claims (C3/C5/C6/C8) + citations (C4) + synthesis-faithfulness (C9) +
> deterministic (C10–C17)**; the verdict is the AND across those four.) Eighth and
> FINAL pre-registered
> `/qa` target — the top-of-tree readiness gate, fired after all six paper groups
> certified: trunk + group1–6). Inherits criteria.md C1–C17 + the v4.62.1 run-shapes
> protocol (first cert = FULL run; delta cycles between; final certifying FULL only
> after a clean delta; Sonnet-tiered citation reviewers with two seeds each). The
> synthesis layer's **backing is the already-certified corpus** — its job is
> *faithful aggregation*, so the branch-defining risk is **cross-group fidelity: the
> field guide must restate each certified result at its certified TIER and each
> headline NUMBER at its certified value, resurrect no withdrawn/descoped result,
> and never promote the α-combination past Observation.**

**Scope (the synthesis target):**
- **PRIMARY — the field guide** `papers/synthesis/geovac_field_guide.tex` (the
  cross-group "what is GeoVac" narrative identity document; ~931 lines, 27 external
  bibitems, inline `\begin{thebibliography}`). This is the one **uncertified**
  synthesis document — it gets the full review across every dimension.
- **SECONDARY — the six group syntheses** (`group{1..6}_*_synthesis.tex`). Each was
  already the **C9 dimension of its own branch cert** (all six branches CERTIFIED),
  so they are **not** re-reviewed in full here; they are checked for (a) **mutual
  cross-group consistency** and (b) **faithfulness to the now-frozen certified
  papers** — carried by the deterministic gates (whole-`synthesis` scope) + one
  cross-group claims pass, not six fresh per-synthesis dispatches.

**Deterministic `--gate`:** `synthesis` (substring-matches every
`papers/synthesis/*.tex` — the field guide **and** the six group syntheses).

## Branch deltas (the only non-inherited content)

### Branch-defining criterion 1: cross-group tier fidelity (the field guide restates certified results — C3/C7/C9)

The field guide narrates results from all six groups. Every load-bearing restatement
must carry the **certified tier** of its source paper — no promotion. The
`claims-reviewer` (field guide) + the cross-group `claims-reviewer` must enumerate
and verify each against the certified corpus:
- **WH1 = PROVEN unconditional** (Paper 38, state-space GH); **WH5 = Observation**
  (α is a projection constant, not derivable); **WH7 = REGISTERED, not proven**;
  **WH8 = tested-negative** (Born).
- **Paper 2 K = π(B+F−Δ) = Observation** — never derived/conjecture/theorem. The
  field guide has a bullet literally headed **"The α-derivation"** (≈L610) whose body
  correctly frames it as an unresolved coincidence (WH5, twelve mechanisms); verify
  the *heading + body together* assert no more than Observation (a heading that reads
  as a derivation claim, uncured by the body, is MATERIAL under the enumerate-and-quote
  rule). C12 backstops the formula; the semantic pass owns the framing.
- **Descoped/withdrawn results must stay withdrawn:** Paper 45 DESCOPED (K⁺ theorem);
  the Lorentzian "literal identification at the Krein level" (withdrawn 2026-06); the
  depth-linear residual form (falsified → Layer-2-presence bound); the withdrawn
  Pythagorean C₃<1 form. Any resurrection is a MATERIAL zombie (C16 backstops the
  known phrases; the reviewer owns new ones).
- **κ = −1/16 = Observation** (matches the geometric 1/16; no derivation bridge),
  never "derived".

### Branch-defining criterion 2: headline-number fidelity (C8/C17)

Every headline NUMBER the field guide quotes must equal its certified-paper value.
The C17 registry currently gates the group syntheses but **not** the field guide
(`GROUP4_FILES` omits it) — so the `claims-reviewer` must **enumerate every number**
in the field guide and match it to its source, and any correction ADDS the field
guide to the relevant C17 family (maintenance rule). The enumerated headlines:
- Natural-geometry table (≈L225): He **0.004%** (cusp) / the per-level best results —
  must match the CLAUDE.md §5 / Paper 13 certified table.
- Resource headline (≈L243): **O(Q^2.5)** composed Pauli scaling, **O(Q^1.69)**
  sub-quadratic 1-norm (R²=0.997), **37-system** library — must match certified
  group4/Paper 14 (verified GROUNDED at draft: paper_14 L37/635/758 = 1.694/0.997).
  The Pauli-multiplier "near parity with STO-3G" framing must carry the v4.59.0 M-A
  honest read (not the retired raw-2.7× headline).
- α block (≈L479–615): **α⁻¹ ≈ 137.036**, K = π(B+F−Δ), **twelve** eliminated
  mechanisms, agreement 8.8×10⁻⁸ — must match Paper 2 / the group5 cert.

### Branch-defining criterion 3: aggregation faithfulness (C9 — the core dimension)

The field guide is a *synthesis of syntheses*. The cross-group `claims-reviewer`
verifies it faithfully aggregates the six certified group syntheses + the papers:
no claim that overstates convergence between groups, no result attributed to the
wrong group, no "the framework proves X" where X is a working hypothesis, and the
§1.5 dual-description / no-ontological-priority rhetoric held throughout (the field
guide is the most narrative document in the corpus — highest interpretive-drift
surface). The "forced/free seam" and "packing → periodic table" origin-story framing
(WH3) may be stated as the project's internal reading but not as established physics.

### Branch-defining criterion 4: code-backing is OUT OF SCOPE (C1/C2 waived — PI direction 2026-07-04)

The synthesis layer has **no own backing tests** — every field-guide headline traces
to an already-certified paper claim, and **all six group syntheses already passed
their branch certs (C9 + the papers' C1/C2)**. Per PI direction at freeze, the
code/test-backing dimension is therefore **WAIVED / out of scope for the synthesis
target** — no `code-reviewer` is dispatched and its absence does **not** force
INCONCLUSIVE (it is a pre-registered scoping decision, not an unexercised gating
dimension). Headline-number correctness is still guarded — by **C17** (deterministic)
+ the **claims-reviewer's** enumerate-every-number pass against the certified corpus
(criterion 2) — so a wrong *quoted* number is still MATERIAL; it is only the
*re-running of the underlying tests* that is waived (they were run at branch-cert
time). A field-guide number with no certified-paper source is still a coverage gap → raised.

## C8 headline enumeration (the field guide — the C17 families to extend as needed)

- **Identity:** discrete spectral triple = Marcolli–vS 2014 graph-network instance at
  the data level + forced-graph claim + master Mellin engine; transcendentals =
  cyclotomic mixed-Tate periods. Tier = the corpus's certified structural claims.
- **κ = −1/16** (Observation); **λ_n = −(n²−1)** = Fock-projected continuum spectrum
  (NOT a bare-graph property — C6).
- **Natural-geometry hierarchy** table: He 0.004% cusp / 0.022% raw / 0.19% CI; H₂⁺
  0.0002%; H₂ 96.0% D_e; LiH R_eq 5.3%; etc. (Paper 13/11/15/17 certified values).
- **O(Q^2.5)** Pauli, **O(Q^1.69)** 1-norm (R²=0.997), **37-system** library (Paper 14).
- **α⁻¹ ≈ 137.036**, K = π(B+F−Δ) Observation, twelve mechanisms, 8.8×10⁻⁸ (Paper 2).

## Known logged gaps at freeze (not blockers; verify still-logged)

- **C17 does not yet gate the field guide** (`GROUP4_FILES` etc. omit it) — expected
  NO-REGISTRY-COVERAGE on the field-guide numbers; the claims pass carries them this
  run, and corrections extend the families.
- The corpus-wide dangling-`debug/`-refs debt (§9 standing) — C14 advisory.
- The six group syntheses are **certified**; a genuine NEW defect in one is possible
  (post-freeze drift) but the expected yield is cross-group-consistency NITs, not
  fresh MATERIAL.

## Change log
- 2026-07-04 — **Citation RE-DELTA = CLEAN + calibrated → synthesis layer CERTIFIED ✅ (v4.67.3; PI-directed clean-citation-delta = pass).**
  1 Sonnet citation reviewer, 2 fresh seeds (latremoliere_metric_st vol 404→402, deligne2010 IHES 112→114).
  **Sensitivity 2/2; specificity clean** — all 7 v4.67.2-remediated cites re-verified CORRECT (hekkelman2022=2111.13865, hekkelman_mcdonald honest-flag, latremoliere_metric_st 404/108393, hawkins2000 title, latremoliere2025 title, mondino_samann no-"Synthetic", connes1995 label 1994). Deterministic 7/7 whole-target. **Zero MATERIAL.** One genuine non-seed NIT: `minguzzi_suhr2024` title "Bounded Lorentzian metric spaces" → real "Lorentzian metric spaces and their Gromov–Hausdorff convergence" (arXiv:2209.14384 / LMP 114/73 correct — referent right; non-load-bearing) — fixed on sight (v4.67.3). Worktree removed, no leak (seeds 402/114 absent; correct 404/112 present; both wrong arXiv IDs gone). **Verdict per PI direction: clean calibrated citation delta = the certifying pass.**
- 2026-07-04 — **Citation-dimension delta (PI-directed, "clean = pass") = DEFECTS → LARGEs remediated in-run; NOT a pass.**
  1 Sonnet citation reviewer, 2 seeds (re-corrupted the two v4.67.1-fixed cites: latremoliere2018 368→366, glanois 160→162).
  **Calibration: 2/2 seeds caught; specificity clean** (the 5 remediated cites verified: connes_vs2021/connes1995-text/brown2012 correct; the two "still-wrong" were the seeds). Deterministic 7/7 whole-target.
  **But the full group1-synthesis bibliography scan surfaced GENUINE non-seed MATERIAL (verified vs primary sources) — in a CERTIFIED branch the group1 cert missed:**
  (i) `hekkelman2022` arXiv **2206.13744 = a Kerr–Melvin black-hole imaging paper** (Hou et al.) — real: "Truncated geometry on the circle," Lett. Math. Phys. 112 (2022) 20, **arXiv:2111.13865**. FIXED.
  (ii) `hekkelman_mcdonald2024` arXiv **2403.18619 = an OpenMP shortest-path CS paper** (Calderón et al.) — the real "Spectral truncations of Tᵈ" ID is unconfirmable via search; the wrong ID REMOVED + flagged pending author confirmation.
  (iii) `latremoliere_metric_st_2017` "Adv. Math. **415** (2023), 108876, 88pp" → real **404 (2022), 108393, 56pp** (arXiv:1811.10843 journal-ref). FIXED.
  SMALLs fixed: `connes1995` in-text label 1995→1994; `mondino_samann2025` title dropped extraneous "Synthetic". SMALL/NIT remaining (PI direction): `hawkins2000` title-mashup, `latremoliere2025_hypertopology` title paraphrase, `farsi_latremoliere2024` (no arXiv ID, UNVERIFIABLE).
  All three wrong-arXiv/venue defects sat in group1's load-bearing "position against the literature" paragraph (L297–303). **Implication (raised to PI):** the group1 branch cert's citation coverage had a real gap; a full citation re-sweep of the group syntheses may be warranted. **Next:** finish the remaining SMALLs, then re-run the citation delta (clean = pass); the fixes are uncommitted pending PI review.
- 2026-07-04 — **1st cert (FULL run) = INCONCLUSIVE (setup error) → genuine MATERIAL found + remediated → field-guide re-run CLEAN.**
  Panel: 3 LLM agents (2 Opus claims [field guide + cross-group/6-syntheses] + Sonnet citation; code dim WAIVED per freeze) + a re-run.
  **Calibration: sensitivity 6/6 planted seeds (run-1) + 2/2 (re-run); specificity clean.** Deterministic C10 (7/7 compile) + C11–C17 PASS.
  **SETUP ERROR (PM):** the seed worktree was built from HEAD, but `geovac_field_guide.tex` carried *uncommitted* working-tree edits (pre-session PI fixes: He 0.019%→0.004%, 28→37 molecules, two→six syntheses). So the claims dimensions reviewed STALE field-guide text → their verdict on the real corpus is untrustworthy ⇒ **INCONCLUSIVE** (primary target not validly exercised). Lesson: for any target with uncommitted in-scope edits, SYNC the working-tree files into the worktree before seeding (as done for the group6 untracked synthesis).
  **Genuine MATERIAL found by the VALID dimensions (bib unchanged HEAD↔worktree; verified vs primary sources) + remediated:**
  (i) group1 synthesis `latremoliere2018` — "The Gromov–Hausdorff propinquity, TAMS 370 (2018)" → real "**The Quantum** Gromov–Hausdorff Propinquity, TAMS **368** (2016)" (arXiv:1302.4058); load-bearing at group1 L1418; a genuine drift in a *certified* branch (the cross-cutting synthesis cert's value).
  (ii) field guide `glanois` journal — "J. Number Theory 182, 36–90 (2018)" → real "**160, 334–384 (2016)**".
  (iii) field guide D6 — Paper 49 "retracted vs closed" tension scoped to *metric-level* retraction (OSLPLS algebra-level bridge survives, Q2′).
  NITs fixed: `connes_vs2021` p.2059→2067; group1 `connes1995` y.1995→1994; CLAUDE.md §2 κ "Derivable from Fock" → Observation.
  **Re-run** (worktree correctly synced to the current corpus; fresh 2 seeds He-0.4%-drift + Lorentzian-re-zombie): **CLEAN + calibrated** — both defects were the planted seeds; controls (κ=Observation, α-block, WH1, retraction paragraph, correct Lorentzian sibling loci, 37-system, C6, §1.5) verified faithful; 2 secondary items PM-verified faithful (5,864 residuals = Paper 56 L1720; CC φ(2)/φ(1) dimensionful = consistent with Paper 51's dimensionless φ(2)/φ(1)²).
  **Net:** the corrected corpus is clean, but this INVOCATION cannot PASS (setup error + genuine MATERIAL). Next: a fresh certifying FULL run on the corrected corpus (correctly-synced worktree) — expected clean. Seeds: `debug/qa/synthesis_seed_key.json` (gitignored). Worktrees removed, no leak.
- 2026-07-04 — **FROZEN (PI-confirmed).** Scope (field guide full + six syntheses
  cross-consistency), four branch-defining criteria, C8 field-guide headline
  enumeration. PI froze as drafted and **WAIVED the code/test-backing dimension** for
  the synthesis target ("skip code tests on the synthesis layer — the group synthesis
  papers already passed"): C1/C2 is out of scope, gating dimensions = claims +
  citations + synthesis-faithfulness + deterministic, verdict = AND across those four.
  FULL first-cert run proceeding.
