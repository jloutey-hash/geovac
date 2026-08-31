<!-- CERT-STALENESS-BANNER -->
> ### ⚠ RE-CERTIFICATION OWED
> This record certifies the state as of **2026-08-30**. Since then **4 `.tex` changed**: group4_quantum_computing_synthesis.tex, paper_14_qubit_encoding.tex, paper_20_resource_benchmarks.tex, paper_23_nuclear_shell.tex.
>
> **The CERTIFIED verdict below is therefore historical, not current.** Do not cite it as present-tense status.
>
> Re-measure rather than trusting this banner — it is itself a snapshot and will go stale the same way:
> `python debug/qa/check_cert_staleness.py --detail`
<!-- /CERT-STALENESS-BANNER -->

# Group 4 (Quantum computing) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only group4-specific scope + deltas + the branch watch-notes.

> **STATUS: CERTIFIED ✅ 2026-07-02 (v4.63.0)** — the 8th `/qa group4` run (the FULL
> certifying pass, fired after delta-1 CLEAN under the v4.62.1 run-shapes protocol)
> **PASSED**: 13/13 sensitivity (incl. all 10 Sonnet-tier seeds) / 6/6 specificity, zero
> remaining verified MATERIAL after three in-run loop-until-dry cycles, deterministic ×7
> green, all new pins green, zero seed leakage. Fourth certified branch (after group3
> v4.21.1, group1 v4.49.0, group2 v4.51.0). Chronicle in the change log below; canonical
> memo `debug/sprint_group4_8thcert_memo.md`.
>
> *(Original freeze note, retained for the record:)* DRAFTED 2026-06-28 for PI freeze
> (fifth pre-registered `/qa` target; the QC/NISQ/VQE branch). Inherits criteria.md
> C1–C17. Branch-defining risk = **QC-resource-claim honesty + the CF-1
> pair-diagonal-ERI disposition** (encoded as a C8/C3/C5/§1.5 sharpening, not a new
> number). CF-1 resolved **A=disclose** (v4.53.0).

**Scope (non-trunk group4):** the **4 quantum-computing papers** —
**Paper 14** (qubit encoding / Pauli-sparsity — the KEYSTONE), **Paper 16**
(chemical periodicity as $S_N$ representation theory — the group-theoretic foundation),
**Paper 20** (resource benchmarks / 38-molecule library), **Paper 23** (nuclear shell-model
qubit Hamiltonians) — **+ the group4 quantum-computing synthesis**
(`papers/synthesis/group4_quantum_computing_synthesis.tex`, drafted this pre-work; the C9
target). Trunk/foundation papers (0, 1, 7, 22) taken as already-canonical; in scope only
where a group4 paper restates them (C7). Paper 16 also sits on the foundations boundary
(it is the $S_N$ basis for `atomic_classifier.py`); it is reviewed here as the group4
periodicity foundation.

**Deterministic `--gate`:** `group4` (path-substring match → `group4_quantum_computing/*.tex`).

## Branch deltas (the only non-inherited content)

### Branch-defining criterion: QC-resource-claim honesty + CF-1 disposition

This is the group4 analog of group2's benchmarking-rule/guardrail-negative criterion and
group1's descope-accuracy criterion — the highest-risk class for a *quantum-computing*
branch. It is a **sharpening of the shared C8 (headline honesty) + C3 (prose ≤ tier) + C5
(no fitted parameter / no §3 negative suppressed) + §1.5 (benchmarking + positioning)**,
not a new numbered criterion.

The reviewers (claims-reviewer, per paper, enumeration-forced) must verify ALL of:

1. **Resource-claim baseline + matched-axis honesty (§1.5 benchmarking rule).** Every
   sparsity / Pauli-count / 1-norm advantage *names the Gaussian baseline it is measured
   against* (STO-3G / cc-pVDZ / cc-pVTZ / the published `trenev2025` exponents) and states
   whether the comparison is at **matched qubit count** or **matched accuracy** — never
   conflated. The standing caveat "*these are matched-qubit, not matched-accuracy*
   comparisons; the composed basis carries 5–26% classical $R_{\rm eq}$ errors" must accompany
   every composed Pauli-advantage headline. A favorable multiplier stated without the
   matched-axis caveat, or measured only against the weakest baseline, is MATERIAL.

2. **Everything is rule-sensitive (the CF-1 partition, INVERTED 2026-08-30).**
   ⚠ **This item previously said the opposite and was wrong.** It classed scaling exponents as
   *ROBUST* under CF-1 — on the premise that "the pair-diagonal vs global-$M_L$ re-pricing is a
   constant factor within valence class, so the log-log slope and the linearity are unchanged" —
   and licensed them as MEASURED. **That premise is falsified.** The "constant $2.51\times$" was a
   single-point artifact, and *every* exponent in the exempted list moved. A reviewer following the
   old text would have skipped exactly the class that turned out to be wrong. There is now **one
   tier, not two**: both exponents and multipliers are rule-sensitive and both must be checked.
   - **Scaling exponents — MOVED, all of them.** Canonical: atomic $O(Q^{3.8})$ (retired: $3.15$);
     composed within-molecule $3.17$ small-basis fit, local slope $\sim\!3.8$ by $n_{\max}=4$
     (retired: $O(Q^{2.5})$); Paper 20 two-point ($n_{\max}=1,2$) $\alpha = 2.816$ **identically across
     all six molecules** (retired: $2.21 \pm 0.02$), forced by $1 + \log_5 18.6$; $N_{\rm Pauli} =
     27.90\,Q$ main-group / $30.03\,Q$ $d$-block (retired: $11.10$ / $9.23$). The 1-norm $O(Q^{1.69})$
     and QWC $O(Q^{3.36})$ exponents are the two that were **not** re-derived in the correction —
     treat them as UNVERIFIED under the exact rule, not as robust.
   - **Absolute multipliers / market-test lines** — sensitive as before. LiH re-prices 333→837;
     equal-qubit H$_2$O $54\times$–$317\times$ (retired: $51\times$–$1712\times$); cc-pVDZ $76\times$
     (retired: $190\times$); the $d$-block coefficient $9.23 \to 30.03$, becoming **denser** than
     main-group's $27.90$ — the "$d$ sparser" claim REVERSES.
   - **Reviewer instruction.** Any exponent or multiplier stated as MEASURED without being
     traceable to a post-2026-08-29 measurement is MATERIAL. Do not treat a clean log-log fit as
     evidence: this whole class read $O(Q^{2.5})$ for months *because* the fit was clean.
     See `debug/qa/delta4_run_notes.md` + `debug/sprint_eri_evaluator_defects_memo.md`.

3. **CF-1 is DISSOLVED (2026-08-29, PI direction) — supersedes the former "DECIDED A".**
   ⚠ **This item previously recorded a decision that no longer holds.** It said the PI chose
   **option A** — keep the pair-diagonal rule ($m_a=m_c \wedge m_b=m_d$) as a *quality QC
   sparsity approximation*, disclose it, and leave **B (global-$M_L$) deliberately on the table**.
   There was never a choice to make: **the pair-diagonal "convention" was a wrong-sign Gaunt
   argument** (`q = m_c - m_a` where the Wigner-3j bottom row requires `q = m_a - m_c`), found in
   seven production modules. It is a bug, not a rule. The PI's call is the **exact global-$M_L$
   rule everywhere**, with the corpus re-priced honestly.
   **The cert criterion, re-aimed.** What is MATERIAL is no longer "failing to name which of two
   live rules is in use" — there is only one rule. It is now: **any pair-diagonal-era number
   presented as a live value**, and **any surviving framing that treats the two as alternative
   conventions** (a *framing zombie*: the arithmetic may be right while the framing implies a
   choice that does not exist). Enforced by the `claims-reviewer` (enumerate every
   sparsity/density/Pauli/exponent claim and trace it to a post-correction measurement), the
   deterministic **C16 entry `pair-diagonal-as-exact-sparsity`** and **C17 entries
   `composed-rule-a-retired-figures` + `composed-retired-scaling-and-counts`**, and the witness
   test `tests/test_paper14_eri_rule.py` (sweeps all seven $c^k$ implementations against an
   independent reference value). Historical mentions remain legitimate **when explicitly marked
   retired/corrected** — that is what the C17 exemption windows are for.
   **Known unremediated instance, disclosed not hidden:** `composed_qubit_relativistic` still
   carries the factor-order half of the defect; both orderings are recorded in
   `tests/test_paper14_spinor_eri_ordering.py` and `tab:spinor_resource` carries a vintage note.
   Reviewers should confirm the note is present and honest, **not** flag the table's numbers as
   defects. See `debug/qa/delta4_run_notes.md` ("SEVENTEENTH SITE").

4. **Nuclear honesty (Paper 23).** The nuclear binding energies at $N_{\rm shells}=2$ are
   **encoding-validation benchmarks, far from experiment** — stated as such, never as a
   nuclear-structure calculation. All nuclear inputs (Minnesota, $v_{ls}/\hbar\omega$,
   $\hbar\omega$) are standard shell-model literature, not GeoVac-derived. The GeoVac
   contribution is the **angular machinery + qubit encoding only**, NOT the conformal
   projection (Fock rigidity theorem: $S^3$ is unique to Coulomb, does not extend to HO/
   Woods–Saxon). Any claim that GeoVac *predicts* nuclear structure is MATERIAL.

5. **Periodicity honesty (Paper 16).** The four-type ($A/B/C/D$) classification and
   $\nu=N-2$ **map onto, rather than predict, the known periodic table** — stated as such
   (no claim that GeoVac *derives* the periodic law). The Dirac instability is a **metric**
   singularity on $S^{3N-1}$, not a topological one; the topology is smooth through
   $Z=1/\alpha$. No fitted parameter (the $v_{ls}$, $d_{\ell\ell}$ in Paper 23 are literature
   inputs, labeled as such — C5).

### Per-criterion watch-notes

- **C8 (headline honesty), per-paper — the enumerated headlines + tiers.**
  *(⚑ marked the old CF-1 "pending FREEZE disposition" tier. CF-1 is DISSOLVED — see item 3 —
  so ⚑ now means simply: re-priced 2026-08-29/30, verify against a post-correction measurement.)*
  - **Paper 14 (keystone)** — ⚑ **entire list re-priced 2026-08-30**; the values below are the
    canonical ones, and the retired figures are given so a reviewer can recognise a zombie on
    sight. Atomic **$O(Q^{3.8})$** Pauli (retired: $O(Q^{3.15})$) vs Gaussian $Q^{4.25}$ LiH /
    $Q^{3.92}$ H$_2$O (`trenev2025`); composed within-molecule **$3.17$** small-basis fit, local
    slope $\sim\!3.8$ by $n_{\max}=4$ (retired: $O(Q^{2.5})$ "universal, exponent spread 0.02" —
    the *spread* claim is also retired: the Paper 20 two-point exponent is now $2.816$
    **identically**, spread $0.000$); $N_{\rm Pauli}=27.90\,Q$ exact, $30.03\,Q$ $d$-block
    (retired: $11.10$ / $9.23$); per-block $279$ non-identity Pauli from $107$ ERIs
    (retired: $111$ / $65$); ERI density $\sim 1/M$ (corrected v4.54.0).
    **QWC groups $O(Q^{3.36})$ and 1-norm $\lambda\ O(Q^{1.69})$, $R^2=0.997$ were NOT re-derived
    under the exact rule — treat as UNVERIFIED, not as canonical.** ⚑ equal-qubit advantage
    $54\times$–$317\times$ (retired: "two-or-more orders of magnitude / 51×–1712×"); the
    "$d$-block sparser" claim **REVERSES** ($30.03 > 27.90$); matched-qubit-not-accuracy caveat
    MANDATORY.
  - **Paper 16:** $\mu_{\rm free}=\nu(\nu+3N-2)/2$ (SO(3N) Casimir); $\nu=N-2$ universal for
    $S<N/2$; **5 atom types A/B/C/D/E** (FIXED 2026-06-28 — abstract+§IV synced to the Table+
    conclusion+code; was a stale "4"); **"maps onto, rather than predicts"**; Dirac instability =
    metric (not topological) singularity, smooth through $Z=1/\alpha$.
  - **Paper 20:** ⚑ LiH composed **838 Pauli @ 30q vs STO-3G 907 @ 12q** — near parity on raw
    counts, and $3.0\times$ MORE than the qubit-reduced 276 (retired: "334 Pauli, 13× fewer QWC",
    which read as a decisive win; under the exact rule the decisive wins are instead the exact
    $Q$-linearity and the $76\times$ cc-pVDZ reduction); balanced coupled (PK-free)
    binds LiH at **$R_{\rm eq}=3.227$ bohr computed (7.0% above the experimental 3.015), 878 Pauli @ 30q, 0.20%** single-point energy at
    the minimum; **row-conditional** chemistry-accuracy (first-row binds; second-row NaH↓
    monotone overattraction — the honest §scope_boundary); library **37 systems** (35 composed
    + He + H2, $Z=1$–56 H–Ba; decided 2026-06-28, this watch-note synced 2026-07-01); $O(Q^{2.5})$ (retired-rule vintage 2026-08-29) universal vs
    Gaussian $O(Q^{3.9-4.3})$; ⚑ $11.10\,Q$ / $9.23$ $d$-block; frozen cores enter via
    identity only.
  - **Paper 23:** HO closures **2,8,20,40,70,112** from graph state counting; magic
    **2,8,20,28,50,82,126** at $v_{ls}/\hbar\omega\approx0.17$, $d_{\ell\ell}/\hbar\omega
    \approx0.02$; deuteron **16q / 688 Pauli (80 Z-only + 608 XY) / 1-norm 383.7 MeV**;
    He-4 **16q / 828 Pauli / 1-norm 511.8/507.2 MeV** (no-Coulomb/Coulomb)
    (12.25× Hilbert, 1.20× Pauli); composed nuc-electronic deuterium **26q / 710 Pauli**,
    ~$10^{13}$ scale ratio; **Fock rigidity theorem**; binding energies = encoding-validation
    benchmarks (far from experiment, stated as such).
    *(Numbers corrected 2026-08-24 to the v5.0.0 retraction values — the moshinsky N_tot
    guard removal moved deuteron 592→688 / He-4 712→828 / 1-norms 342→383.7, 467/462→511.8/507.2,
    composed 614→710. This is a criteria-currency correction to a PI-approved retraction, NOT a
    relaxation; the +20.3% structural claim survives, 828/688 = 712/592. The stale pre-retraction
    numbers were the ones this DoD had frozen at the v4.63 cert.)*
- **C4 (citations) — QC-resource + nuclear surface.** Verify the comparators say what is
  attributed: `trenev2025` ($Q^{4.25}$ LiH / $Q^{3.92}$ H$_2$O Gaussian exponents), Szabo–
  Ostlund (H$_2$ integrals), the STO-3G/cc-pVDZ/cc-pVTZ Pauli/1-norm comparators, the
  OpenFermion/Qiskit/PennyLane export targets, Jordan–Wigner; for Paper 23: Minnesota NN
  potential, Moshinsky–Talmi brackets, Mayer–Jensen spin–orbit, the magic-number/shell-model
  literature. Wrong exponent, misattributed baseline, or fabricated arXiv ID = MATERIAL.
- **C3 (prose ≤ tier).** Scaling exponents / linearity / ERI-density-decay = MEASURED
  (log-log fits, report $R^2$ + residual per §13.4a scaling-law rule); angular selection
  rules (Gaunt/3j) = exact/algebraic; the composed sparsity *advantage magnitude* = MEASURED
  under a named ERI rule, never "exact". Zero-parameter construction asserted as exactly that.
- **C6 (discrete-vs-continuum precision).** Pauli counts / ERI density are properties of the
  *graph-derived* Hamiltonian after JW; where a paper invokes the Fock $S^3$ equivalence or
  the angular-sparsity theorem (Paper 22) it states the current tier. Paper 23's Fock
  rigidity theorem (Coulomb-unique) is stated as the theorem it is.
- **C7 (trunk-dependent status).** Where a group4 paper restates Fock $S^3$, $\kappa=-1/16$,
  the angular-sparsity theorem (Paper 22), or the natural-geometry hierarchy, it states the
  current tier (Paper 7 = 18 symbolic proofs; $\kappa$ = Observation; Paper 22 angular
  sparsity = potential-independent theorem, density 6.06% universal / 1.44% pair-diagonal).
- **C9 — the new group4 synthesis** is the synthesis-faithfulness target (separate
  claims-reviewer dispatch): every claim traces to one of the 4 papers; the QC-resource
  spine (encoding 14 → resources 20 → periodicity foundation 16 → nuclear extension 23) is
  faithful; the matched-qubit caveat + the CF-1 disposition + the nuclear/periodicity honesty
  ceilings are all carried; no §3 dead-end (the 5 relativistic-Z₂-tapering negatives, the
  cross-block-h1 16× over-binding, the non-abelian-gauge Pauli-reduction NO) re-asserted as working.

## Known drifts surfaced in pre-work (cert-time disposition)

The pre-work claim→test mapping (`docs/claim_test_matrix.md` group4 section, 2026-06-28;
all cited tests RUN GREEN: 254 + 105 + 78 passed) surfaced four cross-corpus inconsistencies
and a coverage profile the cert reviewers should treat as enumeration targets. **Backfill +
inconsistency pass done 2026-06-28** (PI-directed "address the backfill/inconsistencies first"):

1. ✅ **CF-1 pair-diagonal ERI** — **DISSOLVED 2026-08-29 (PI direction); supersedes the former
   "DECIDED A (disclose)".** The pair-diagonal rule was a wrong-sign Gaunt argument, not a
   convention, so there was no disposition to choose (see item 3). Its quantification is retired
   with it: the "constant 2.51× main-group / 3.25× $d$-block, scaling robust" line was a
   **single-point artifact** — the re-pricing is not a constant factor and the exponents moved.
   Canonical: exact global-$M_L$ everywhere; LiH 333→838 @ Q30 (near parity vs raw STO-3G 907).
   Applied: the `sec:eri_rule` disclosure subsections in Papers 14/20 + the qualified d-block/
   market-test/TM-caption loci; the **shared [criteria.md "Dual-rule ERI framing"](criteria.md)**
   rule + the **C16 `pair-diagonal-as-exact-sparsity`** deterministic backstop + the
   **`tests/test_paper14_eri_rule.py`** characterization test (product realizes A; B re-prices
   2.51×/3.25×) + the `angular_zero_count` docstring fix. **Cross-branch flag (NOT auto-fixed):**
   group2 Paper 13/fci_atoms "exact rational Slater integrals" doesn't disclose the pair-diagonal
   angular subset (radial R^k exact, angular = A; energy-negligible) — for a future group2 re-touch.
2. ✅ **Library size** — **DECIDED 37 (ship as-is, PI direction 2026-06-28) + applied corpus-wide.**
   `_SYSTEM_REGISTRY` = **37 systems** (35 composed molecules + He + H2). Synced everywhere:
   Paper 14 (28/30→37 systems / 35 composed), Paper 20 (38→37; the **3 non-buildable organics
   CH$_2$O/C$_2$H$_2$/C$_2$H$_6$ removed** from the abstract + multi-center table, 8→5 multi-center),
   the synthesis, ecosystem docstring (28→37), test comment (40→37), claims_register (40→37),
   CLAUDE.md §2 + §1.1 + §1.5 (40→37; the §1.1/§1.5 count edits applied under explicit PI "37"
   direction). All edited papers + synthesis compile clean. *(Pre-existing, NOT from this pass:
   Paper 20 has 4 broken section-`\ref`s — sec:spinor\_composed/composed, subsec:spinor\_scope —
   and a missing `Childs2021` bibitem; surfaced incidentally, flag for the cert run.)*
   *(**RESOLVED 2026-08-24**, retraction re-review delta — it was a **C10-gate procedure artifact,
   not a paper defect.** `Childs2021` IS in `paper_20_refs.bib` (l.181); paper_20 is the only active
   paper using an external `\bibliography{}`, and `check_compiles.py` ran pdflatex-only, so every
   citation came back undefined. Under the correct revtex build paper_20 compiles with **zero**
   undefined refs/cites — but it needs a **double** bibtex cycle because `\cite{Childs2021}` sits
   inside a `\bibnote`. Fixed the gate to run bibtex when the `.aux` carries `\bibdata`;
   discrimination re-proven (still catches a bogus cite AND a bogus ref). C10 `--gate group4` now
   PASSES all 5 documents.)*
3. ✅ **Paper 16 structure-type count** — **FIXED**: synced to **5 types (A/B/C/D/E)**
   consistently (abstract + §IV intro + new Type E subsection now match the Table + conclusion +
   code + CLAUDE.md §1.5); code's E/F split noted as an implementation refinement. P16 compiles.
4. ✅ **Paper 23 1-norm drift** — **FIXED + PINNED**: deuteron 227→**342 MeV**, He-4 557/552→
   **467/462 MeV** (code stable since v2.7.0 — drafting-era staleness, not a regression; term
   structures always correct). `test_paper23_resource_counts.py` pins all three + the counts.
5. ✅ **LiH composed Pauli 333 vs 334** — **RESOLVED (not a defect)**: 334 is the shipping
   `hamiltonian()` API value (= paper + CLAUDE.md); 333 is the raw-builder count (excludes the
   identity term). No fix needed.

**Coverage profile — backfill done: 7 of 14 NO-TEST gaps CLOSED.** The *structural* spine was
already BACKED-SOUND (JW=FCI, QWC correctness, antisymmetrisation, magic-number recovery, qubit
counts, linear-in-Q ratios). **The four keystone scaling exponents are now PINNED from GeoVac
data** (`test_paper14_scaling.py`, slow): atomic $O(Q^{3.15})$ (3-pt 3.10), 1-norm $O(Q^{1.69})$
(3-pt 1.67, sub-quadratic), QWC $O(Q^{3.36})$ (3-pt 3.356), composed $O(Q^{2.5})$ (retired-rule vintage 2026-08-29) (~2.5,
CF-1-robust) — replacing the synthetic-only `test_fit_scaling`. The nuclear counts + 1-norms are
pinned (`test_paper23_resource_counts.py`). **7 gaps remain OPEN** (none a §C8 headline blocker):
ERI $1/M^2$ (decays slower — BACKED-WEAK, consider rewording to ~$1/M$); double-factorization
rank; 13× QWC factor; balanced-LiH 878/0.20%@n=3 (cross-ref group2 P19); per-pair-4 asymptotic;
$Z=137$ metric-not-topological; Fock rigidity theorem; the He-4 1.20×/12.25× ratios.
`test_paper14_revision.py`
is a file-string sanity check, NOT physics backing — do not count it as coverage.

## First bite (proposed; PI confirms at FREEZE)

- **Whole-group** — all 4 papers + the new synthesis in one `/qa group4` invocation (same
  granularity as the certified group1/group2 whole-group runs; group4 is small — 4 papers).
  Per `qa.md` step 4: deterministic layer (C10–C16) whole-group first; **code-reviewer
  1 paper/agent** (the papers with backing tests — 14/16/20/23 all have tests); **claims-
  reviewer chunked** (≈2 chunks: {14, 20} resource-claims together since they share the
  Pauli/CF-1 surface; {16, 23} foundation+nuclear); **citation-reviewer** one chunk (4
  papers); **C9 synthesis** one dispatch; **completeness-critic** one agent. Per-chunk
  seeding: ≥1 seed catchable by each (dimension × chunk), planted in a non-first paper of
  each chunk, in the throwaway worktree only. **Recommend batching the FIRST cert** (qa.md:
  fresh group → smaller remediation cycles), though at 4 papers whole-group is also tractable.

## Change log
- 2026-06-28 — **DRAFTED** by PM for PI freeze (fifth pre-registered `/qa` target; first
  QC branch). Inherits criteria.md C1–C16. Branch-defining risk = QC-resource-claim honesty.
  The CF-1 "disposition" question is CLOSED — dissolved 2026-08-29, the pair-diagonal rule being
  a bug rather than an alternative convention. ⚠ The quantification once recorded here
  (`debug/qa/group4_cf1_library_sweep_memo.md`: "constant 2.51× main-group / 3.25× $d$-block;
  scaling robust") is **retired as a single-point artifact** — superseded by
  `debug/sprint_eri_evaluator_defects_memo.md` + `debug/qa/delta4_run_notes.md`. **The former
  FREEZE requirement to resolve A-vs-B is void — the §C8
  ⚑ numbers depend on it.** Also flagged: library-size inconsistency (14:"30" / 20:"38" /
  CLAUDE.md:"40" / 35 shipping). C9 = the new group4 synthesis (drafted this pre-work;
  `papers/synthesis/group4_quantum_computing_synthesis.tex`, 4pp, three-pass clean,
  CF-1-honest from the start). **Pre-work complete (4/4 items):** (1) CF-1 library sweep
  quantified; (2) this DoD drafted; (3) synthesis drafted; (4) `docs/claim_test_matrix.md`
  group4 rows populated (was zero; all cited tests RUN GREEN). **Deterministic layer
  pre-validated GREEN** incl. the new synthesis — C11/C13/C14/C15/C16 PASS, synthesis compiles
  clean (C11 caught + fixed two wrong bibitem titles in the draft). The 5 known drifts + 14
  NO-TEST gaps (above) folded in as enumeration targets. First bite = whole-group (proposed; PI confirms).
- 2026-06-28 — **FROZEN + first cert (whole-group) = FAIL → REMEDIATED (v4.54.0).** Panel
  FULLY CALIBRATED (sensitivity **8/8**, specificity **5/5** — every gating dimension's plant
  caught, zero false positives; the disclosed `sec:eri_rule` correctly NOT flagged → the new
  framing dimension is calibrated). Verified MATERIAL defects (seeds excluded), all remediated:
  (1) ~10 C16-dodging **framing zombies** (P14/P20 undisclosed d-block-"cheaper/sparser/
  economical" + "2.7×" market-test loci) — disclosed; **C16 `pair-diagonal-as-exact-sparsity`
  pattern broadened** to catch the dodgers; (2) **μ_free code bug** (`atomic_classifier` `2ν²`
  → Casimir `ν(ν+3N−2)/2`; test + 3 P16 loci re-pinned; diagnostic-only, no Hamiltonian
  affected); (3) **`test_balanced_row2` RED** (12 fails → 39 pass; spec-factory module +
  R-alias kwargs); (4) **trenev2025 misattribution** (15 P14 loci + P20 caption → GeoVac
  OpenFermion recompute, Trenev = methodology + range); (5) **"~1/M²"→"~1/M"** (9 P14 loci);
  (6) **P23 "conjecture"→"observation"** (C5); (7) **Paper-38 "Latrémolière propinquity"→
  "state-space GH"** (P14/P20/P23/synthesis/bibitem); (8) **P20 "38 molecules"→37**. All
  deterministic gates PASS; affected tests green (atomic_classifier 197, balanced_row2 39,
  eri_rule 3); all 5 papers compile clean. **(9) placeholder cites RESOLVED** (PI-directed
  explorer search): `Sunaga2025`→**Swain et al. arXiv:2211.06907** (was misattributed; RaH
  18q/47,099 numbers exact, re-keyed across P14/P20/.bib), `caesura2025`→**PRX Quantum 6,
  030337** (+BLISS-THC), `ChildsBerry`→**arXiv:1501.01715**, `MartinezYRomero2004`→
  **physics/0402061**, `BJL`→**withdrawn** (unverifiable). C4 clean; P14/P20 recompile clean.
  **All 9 material findings remediated** → ready to re-run `/qa group4` for the certified PASS.
  Per-run detail: CHANGELOG v4.54.0.
- 2026-06-28 — **re-cert (whole-group) = FAIL → ALL 4 REMEDIATED; recert HELD per PI (v4.55.0).**
  Panel FULLY CALIBRATED again (sensitivity **8/8**, specificity **5/5** — the v4.54.0-disclosed
  `sec:eri_rule` and the now-correct Swain cite both correctly NOT flagged). 4 verified MATERIAL
  findings (seeds excluded): **(A)** the v4.54.0 `Sunaga2025`→`Swain2022` fill was ITSELF a
  misattribution — the RaH-18q benchmark is **Chawla et al. arXiv:2406.04992 = PRA 111, 022817**;
  **47,099 = two-electron integrals, NOT Pauli; rel Pauli = 12,556** (non-rel 2,740) — re-keyed
  `Swain2022`→`Chawla2024` (P14 bibitem + P20 `@article` + all `\cite`), recomputed the P20
  ratio column on 12,556 (×3.751); **(B)** trenev2025 (vibrational-spectra) still credited with
  the Q^3.9–4.3 range → **methodology cite ONLY** (range/exponents are GeoVac's recompute; P14 7
  loci + bibitem, P20 caption + bib, synthesis); **(C)** P14 §origin framing zombie ("advantage
  that grows with angular complexity") → disclosed as a pair-diagonal artifact (exact rule = d-block
  denser); **(D)** P20 STO-3G market test convention mix — LiH "907@Q12" is RAW while every other
  Gaussian count is 2-qubit-reduced (matches `GAUSSIAN_LIH_PUBLISHED` 276@Q10) → caption + prose
  DISCLOSE that the 2.7× is raw-vs-raw and narrows to parity under uniform reduction (not a number
  reversal). NITs swept: 30→**37 systems** (P14), "Propinquity-derived"→"GH-convergence-derived"
  §heading + 382 (P20). Deterministic gates C11/C13/C14/C15/C16 PASS; all 5 papers compile clean
  (P14 26 / P16 7 / P20 12 / P23 11 / synth 4 pp); no production code edited. **The certified-PASS
  confirmation run is the next `/qa group4` — HELD for PI timing, NOT auto-fired.** Per-run detail:
  CHANGELOG v4.55.0 + `debug/sprint_group4_recert_remediation_memo.md`.
- 2026-06-29 — **confirmation cert (3rd run, whole-group) = FAIL → ALL remediated; recert HELD per PI (v4.56.0).**
  Panel FULLY CALIBRATED (sensitivity **7/7** valid seeds — c1 excluded as inert/dud; specificity **6/6**,
  incl. the citation-reviewer *confirming* the now-correct Chawla2024 cite). The run peeled a deeper layer
  (a full relativistic-table internal-consistency cross-check prior runs' chunking skipped). Findings:
  **(M1)** P14 `tab:sunaga` stale GeoVac native 805/534 vs the paper's own `tab:spinor_resource` 1413/942
  → fixed (advantage 16–24×→**9–13×**, obs-3 QWC 6571→11865; v4.55.0's A-fix had left P14↔P20 inconsistent);
  **(M2)** P20 abstract "binds at R_eq=3.015" vs body's computed **3.227** → reworded (n=2 resource vs n=3
  accuracy separated); **(M3)** 0.20% n=3 headline has no test → logged coverage gap; **(M4)** synth
  1/M²→1/M; **(M5)** P23 §4 title→"First **Two-Species**…"; **(projected)** P14↔P20→honest **17–32×**.
  **(rel λ_ni conflict — "diagnose first" per PI):** code's first-row LiH/BeH rel 1-norm matched neither
  table (BeH n=2 code 143.96 vs table 40.26) while frozen-core CaH matched exactly → diagnosed as a
  **stale table, NOT a code regression** (the table's BeH rel λ 40.26 was physically impossible: 3.5×
  BELOW its own scalar 139.12; the first-row λ path drifted after the table's v2.15.0 vintage, Pauli
  counts unchanged) → tables + obs-2 + P20 prose corrected to code values + **new pinning test
  `tests/test_paper14_rel_lambda.py`** (4 passed; λ was never test-guarded — the gap that let it drift).
  Deterministic gates PASS; all 5 papers compile clean. Deferred NITs (citation source-check):
  Pachucki 2023-vs-2018, rocca/caesura authors, ScH 277/278. **Certified-PASS run is the next `/qa group4`,
  HELD for PI.** Per-run detail: CHANGELOG v4.56.0 + `debug/sprint_group4_confirm_cert_remediation_memo.md`.
- 2026-06-29 — **4th cert (PASS-confirmation, whole-group) = FAIL → ALL remediated; recert HELD per PI (v4.57.0).**
  Panel FULLY CALIBRATED (sensitivity **8/8**, specificity **6/6** — confirmed every v4.55/v4.56 fix). A thin
  converging layer. **HEADLINE — Trenev attribution REVERSED (corrects the v4.54.0 finding 4 + v4.55.0 finding B
  recorded above):** a "diagnose-first" web-verification of Trenev (arXiv:2311.03719) found its **Appendix B /
  Table 5 ("Electronic structure vs Vibrational structure") DOES tabulate the electronic Gaussian JW Pauli counts
  for LiH/H2O — all six GeoVac values (276/5851/63519; 551/8921/107382) match Table 5 exactly.** So the earlier
  "Trenev is vibrational-only / no electronic counts / counts are GeoVac's own recompute" was WRONG (it missed
  Appendix B); the **code (`composed_qubit.py` "Source: Trenev Table 5") was right all along.** Reverted corpus-wide
  to *counts = Trenev Table 5 (App. B); exponents = GeoVac's log-log fit of those published counts* (P14 11 loci +
  bibitem, P20 caption + refs.bib, synthesis body + bibitem, code docstrings); **no numeric value changed.** Other
  findings: **(F2)** P20:957 "/balanced" zombie dropped (balanced LiH binds — flagship result); **(S1)** synthesis
  R_eq=3.015→3.227 (the v4.56.0 M2 fix hadn't propagated to the synthesis); **(F3/F5)** P14 matched-qubit wording +
  P20 RaH species-mismatch caveat; **NITs** ScH 277 non-identity / rocca P.J.~Ollitrault / synth "Pauli counts".
  **New tests:** `test_paper20_library.py` (37 systems) + `test_paper16_dirac_metric.py` (§VI formula/smoothness) —
  6/6 green, closing two recurring no-test flags. Deterministic gates PASS; edited papers compile clean. Deferred
  NITs (non-blocking carve-out): Pachucki year/FW source-check, Navrátil title-venue, duplicate caesura, dangling
  P20 .bib, "M=n_max²" text, rel n=3 + 0.20% n=3 coverage. Per-run detail: CHANGELOG v4.57.0 +
  `debug/sprint_group4_4thcert_remediation_memo.md`.
- 2026-06-29 — **5th cert (whole-group) = FAIL → remediated; CLOSEST to PASS; recert HELD per PI (v4.58.0).**
  Panel FULLY CALIBRATED (sensitivity **8/8**, specificity **6/6**); **6 of 8 dimensions clean-except-their-seed**,
  and **every accumulated v4.55–4.57 fix confirmed accepted** — the Trenev reversal **triple-confirmed** (two code
  reviewers + the citation reviewer reading the PDF Appendix B). A thin framing/attribution-precision layer.
  **(#2)** P20 §VI GH bound: γ is a proven state-space-GH **convergence rate, NOT a direct energy bound** (was
  oversold via "this energy sits within γ" + "error bound", contra the own footnote) → reworded. **(M-A)** the raw
  STO-3G **907** (the 2.7×/13× denominator) was over-blanketed under "Trenev Table 5 (2-qubit reduction)" by my own
  v4.57.0 reversal — but Table 5's reduced STO-3G LiH is **276**; caption now carves out 907 as the raw-JW
  matched-raw baseline. **(M-C)** "Z=1–36" → **"Z=1–56 (H through Ba)"** (registry has SrH Z=38 + BaH Z=56,
  probe-confirmed). **Claims SMALL:** §1.5 "comparable accuracy" (P14:2782) / "comparable basis quality" (P20:475)
  → matched-qubit; CF-1 lead +pair-diagonal qualifier (P14:1580); P14:980 "composed 7–8.8% [Paper17]" → **balanced
  [Paper19]** (composed=5.3%). **NITs:** NaH worked example (showed balanced 239 + non-existent attr) → 223/171.46
  via `H.one_norm`; atomic_classifier docstrings +'F'/'d_block'; Pachucki prose 2023→2018 (bibitem verified to
  support claim); P14 "recomputed"→"published Trenev counts". Deterministic gates PASS; papers compile clean.
  **PI-flagged (not auto-decided):** M-A market-direction (raw-907 2.7× vs Trenev-reduced 276 parity — caption
  brackets both); M-B (0.20% n=3 no CI test, heavy 84q FCI) deferred; "three rows" framing sweep; secondary-number
  reconciliations; citation NITs. Per-run detail: CHANGELOG v4.58.0 + `debug/sprint_group4_5thcert_remediation_memo.md`.
- 2026-07-01 — **6th cert (whole-group, PI-fired) = FAIL → ALL remediated; recert HELD per PI (v4.60.0).**
  Panel FULLY CALIBRATED (sensitivity **8/8**, specificity **6/6**); one calibration wrinkle: the P16 code
  reviewer #1 MISSED its seed (ordering-only divergence asserts blessed as SOUND) → clean discarded, fresh
  strength-matching re-dispatch CAUGHT it (capped-δ counterfactual) — the per-dimension fix-and-re-run path.
  Seed s1 caught twice independently. 3 verified genuine MATERIALs, all remediated: **(M1)** the v4.58.0 M-C
  Z=1–36→Z=1–56 fix had missed P20's conclusion + 4 P14 loci (second-locus propagation class); **(M2)** P14
  ℓ-parity "verified bit-exact spectrum preservation (test test_extended_tapering.py)" — the cited file had
  NO eigenvalue comparison (`_spectrum_lowest` defined, never called) → new
  `test_extended_hopf_ell_spectrum_preserved_h2` PASSES <1e-10, sentence repointed; **(M3)** P14
  tab:multi_center still carried the v4.52.0-de-shipped organics CH₂O/C₂H₂/C₂H₆ (8 rows vs prose "these
  five") → removed + **C16 `organics-in-library` registry entry**. NITs: LiH λ 33.3→32.6 / 0.97×→0.95×
  (closes the v4.59.0-deferred 33.3/32.59 reconciliation; §1.5 echo PI-flagged), FriarPayne→PRA 56,5173(1997),
  goings2022 bibitem added, Migdalek–Bylicki wrong-ID fixed, Pachucki2018 re-key, duplicate caesura merged,
  all six tab:metric rows pinned, 10¹³ guard →1e12, matrix upgrades (Z=137 gap #10 CLOSED; He-4 ratios
  →BACKED-SOUND), DoD watch-notes synced (~1/M; 37 systems Z=1–56). Coverage closure: completeness-critic
  surfaced P14 §V.G (never enumerated by any prior panel) → focused claims re-dispatch CLEAN across all 8
  regions (zero MATERIAL); focused citation re-dispatch found the Migdalek–Bylicki wrong-ID. Post-remediation:
  6 gates PASS, 4 papers 0-errors, 144 affected tests green, zero seed leakage. **Certifying run (7th) HELD
  for PI.** Per-run detail: CHANGELOG v4.60.0 + `debug/sprint_group4_6thcert_remediation_memo.md`.
- 2026-07-02 — **7th cert (whole-group, PI-fired) = FAIL → ALL remediated; certifying (8th) run HELD per PI (v4.62.0).**
  Panel FULLY CALIBRATED (sensitivity **8/8**, specificity **6/6**; fresh seed classes — every seed caught by its
  own agent, three caught cross-dimension; two detected-but-severity-downgraded with sound reasoning, one seed
  genuinely inert in situ). 3 verified genuine MATERIALs, all remediated: **(M1)** P16:292 printed a FALSE identity
  (μ_free/N² = 2−8/N+8/N²; correct 2−6/N+4/N², matching the paper's own table) → fixed; **(M2)** P14:2781 floor
  "190×–1,712×" vs the cited table's own 51× floor → fixed; **(M3)** the magic-number presence scan (min_gap=1e-10)
  is necessary-not-sufficient on a §C8 headline — PM probe: six lower magics ARE the six dominant gaps ≤126 but
  **126 is non-dominant** (0.107ℏω, ten sub-shell boundaries larger; >126 = truncation-edge artifacts) → new
  `tests/test_paper23_magic_gaps.py` (gap-structure pins + dominance + non-dominance bands) + honest P23 disclosure
  sentence; production scan untouched. **Critic follow-through:** live sweep of ALL 12 balanced tab:molecules cells —
  counts exact, **λ_ni column stale 9/12** (0.3–7%) → re-synced + `tests/test_paper20_balanced_lambda.py` pins every
  cell; the run-6 KH fix had never landed (PM propagation miss) AND its 28.15 proved unreproducible → **31.6**
  (factory experimental R, pinned; reproduce-before-syncing reaffirmed). NITs: GH "truncation-error bounds"→
  "truncation-convergence rates" (5 loci; reviewer severity-split recorded), Navratil title↔venue, 13×-QWC qualifier,
  P14 symmetry-adapted sentence honesty, deuteron round-trip docstring, dirac-metric magnitude assert, matrix
  upgrades (878→SOUND, Z=1–56, O(Q^2.5) LiH-fit note, magic→SOUND, +2 §VII rows). Post-remediation: 6 gates PASS,
  5 papers 0-errors, 12/12+3/3+7/7 new/affected pins green, zero seed leakage. **Certifying run (8th) HELD for PI.**
  Per-run detail: CHANGELOG v4.62.0 + `debug/sprint_group4_7thcert_remediation_memo.md`.
- 2026-07-02 — **8th cert (FULL certifying pass, PI-fired post-clean-delta) = PASS → group4 CERTIFIED ✅ (v4.63.0).**
  First full run under the v4.62.1 cost package. Calibration **13/13 sens / 6/6 spec** (2 seeds per
  Sonnet-tiered agent ×5 + 1 per Opus ×3; every seed caught by its own agent; P23-code proved the planted
  "hw-grid interpolation" rationale fabricated by reading the production source; P20-code proved the planted
  justifying comment false by checking the claimed alternate coverage). Three in-run loop-until-dry cycles:
  **(1) panel** → 4 genuine MATERIALs remediated (P14 stale TM-automation ×2 loci vs the registry-verified
  all-ten state; §hopf 6-system/254 claims un-tested → He spectrum test + per_block 254-sum pin + honest
  scoping; P20 FCIDUMP seven-system list → the actually-tested set; matched-qubit→matched-raw-JW-convention
  ×2) + the §BeH₂–H₂O λ-convention re-verification (identity-INCLUDED; composed-BeH₂ 354.9 was
  legacy-builder vintage → 373.4 live, QPE direction survives at 18%; 6 cells pinned + new C17 family);
  **(2) critic→gap-closure** (12 never-enumerated regions) → 1 MATERIAL: **eq:dirac_fs printed (Zα)⁴ vs the
  Z⁴α² its own code/test/sibling-term carry** (the run-7 P16 wrong-printed-equation class) → fixed, 43/43;
  **(3) re-scan of the remediation diff** (18 hunks, seeded) → 1 tail defect my own sync introduced
  (P14:1036 sibling 355/361 → 66/359 live) → fixed, dry. Backfills: P16 heavy-row pins (Kr/Xe/Rn/Og),
  ℓ=2 multipole termination (live-verified tight), balanced-MPO drift guard, 32.6 pin, P23 3×10⁵→5×10⁵
  body locus, citations (Tung Wei-Cheng, Burkat J./N., Chawla v2 title, 5 P16 orphans removed).
  Protocol hardening adopted: commit seeds onto the worktree branch (two reviewers legitimately enumerated
  uncommitted seeds via git diff). Cost ≈2.49M subagent tokens, ~half billed at Sonnet tier, 13 seeds,
  3 cycles. Canonical memo `debug/sprint_group4_8thcert_memo.md`.
- 2026-07-02 — **Delta-verification #1 (PI-fired, first run under the v4.62.1 run-shapes protocol) = CLEAN-DELTA → the full certifying (8th) run is EARNED.**
  Scope = diff 51a8b1a..HEAD (the 7th-cert remediation: 206 insertions across 9 files). Deterministic gates ×7
  (incl. the new C17) whole-target PASS. Three delta reviewers, hunks pasted (no file reads, no worktree, no critic):
  **sensitivity 4/4** (code-Sonnet caught BOTH its seeds incl. an independent top-6 recompute; synthesis caught its
  QWC-qualifier flip; claims #1 MISSED its "126-dominant" self-contradiction seed by constructing a charitable
  reading → de-calibrated, discarded; sharpened internal-consistency re-dispatch caught it immediately — the lesson
  baked into qa.md), **specificity clean** (no genuine hunk falsely flagged). All genuine remediation hunks verified
  across dimensions; citation delta (Navratil title) = the calibrated run-7 reviewer's own verified correction,
  PM-applied, disclosed. **Cost: ~334k subagent tokens (~7.5x below a full run), including one recalibration
  re-dispatch.** Seed key `debug/qa/group4_delta1_seed_key.json`.

---

## Re-review OWED (2026-08-22, v5.0.0) — group4

**Papers in this target changed after certification.** Logged here, at the
owning source, so the certified status is not read as covering text that
post-dates it.

- **Paper 23** (`paper_23_nuclear_shell.tex`): every nuclear resource
  number changed after the `moshinsky.py` N_tot guard was removed — deuteron
  592→**688** non-I Pauli (XY 512→**608**), He-4 712→**828**, 1-norms
  342.2→**383.7**, 466.9→**511.8**, 462.4→**507.2** MeV, composed
  nuclear-electronic total 614→**710** (measured, and now pinned). Qubit counts
  unchanged. **The +20.3% structural claim survives exactly** (828/688 =
  712/592).
- **group4 synthesis**: the same two counts.
- C17 registry entry `p23-nuclear-resource-counts` added; it caught synthesis
  propagation that manual sweeping had missed.
- New backing test `test_paper23_composed_nuclear_electronic_counts` — the
  composed decomposition had **no test at all**, which is how a blanket
  substitution was able to leave a component (688) larger than the stated
  total (614).

**Discharge condition:** a CLEAN DELTA over the changed loci (diff-scoped, per
the qa.md run-shape rule), which is also the standing precondition for the next
FULL certifying pass. The delta must re-test *these specific defects* rather
than trust this entry (qa.md hard rule, added the same day).

**Deterministic layer already re-run and GREEN on this target:** C10 (compiles,
with the aux-clean fix), C13, C14, C16, C17, C18, C19, C5/C12 (now corpus-wide
after the scope fix). What is owed is the LLM-judgment layer: claims, and
synthesis where the target has one.

**DISCHARGED — CLEAN-DELTA, 2026-08-24 (v5.0.9).** Delta scope = the changed loci
(P23 retraction counts + the group4 synthesis's two counts) **plus** the newly-banked
Paper 14 `sec:tc_atomic_sparsity` (v5.0.9). Pasted-hunks delta (HEAD=v4.85.0 with the
whole arc uncommitted, so no worktree; seeds lived only in the reviewer prompts, never
the corpus). Three dimensions dispatched, all **exercised + calibrated**: claims (opus,
2/2 seeds — inverted-honest-scope + a 712≠710 component-sum, both caught), synthesis
(opus, 1/1 — a 1.45× ratio contradicting 828/688=1.20, caught), code (sonnet, 2/2 — a
`n_l3==n_l3` tautology + a `n_nonid>0` weak guard, both caught); specificity clean.
**Zero verified MATERIAL** once seeds removed (real corpus carries `12`/710 and the
honest "not absolute accuracy claims"). Genuine NITs found + fixed on sight: P14 prose
"density stays single-digit" (s+p is 14.8%) -> stated range; the P23 test module
docstring's pre-retraction 592/512/712 "already correct" -> corrected to 688/608/828 +
history; the P14 non-tautology guard extended to lmax=3. Also this run: **C10 gate fixed**
(it ran pdflatex-only, spuriously failing paper_20 — the only external-`\bibliography{}`
paper; now runs bibtex + a double cycle for note-nested `\cite`; discrimination re-proven).
Seed key `debug/qa/group4_delta2_seed_key.json`. The branch stays CERTIFIED; a clean
delta is the standing precondition for the next FULL certifying pass.

**RE-CERTIFIED — FULL certifying pass = PASS, 2026-08-24 (v5.0.9).** Fired after the clean
delta (PI-directed). Seeded isolated copy of all 5 group4 papers (synced from the working
tree, not HEAD=v4.85.0); code dimension = real repo + pasted seeds. Panel: 7 reviewers over
3 waves — claims-A (P14+P16, opus), claims-B (P20+P23, opus), citations (all 4, sonnet→**opus
recovery**), synthesis (opus), code-P14 + code-P23 (sonnet), completeness-critic (opus).
**Calibration: 100% sensitivity** — every planted seed caught (P16 maps→predicts, P23 far→within-few-percent,
McArdle↔Bauer arXiv swap, synthesis cross-block-h₁ dead-end flip, code b1/b2 ×2); **specificity
clean** (zero controls flagged). **Zero verified MATERIAL in the real corpus** — every reviewer
MATERIAL was its planted seed (scratch/prompt-only; real text confirmed: P16 "maps onto not
predicts", P23 "far from experiment", synthesis "16× over-binding", bib IDs correct).

- **Citation recovery (run-6 pattern):** the Sonnet citation reviewer MISSED the McArdle↔Bauer
  swap (didn't fetch + didn't cross-check doi↔eprint) → de-calibrated → re-dispatched on Opus,
  which caught it via the internal-consistency pass. Citation dimension calibrated via recovery.
- **Criteria-currency fix (pre-freeze):** §C8 P23 headline updated 592/712/342/467-462 → the
  v5.0.0 retraction values 688/828/383.7/511.8-507.2 (+ composed 710). Correction, not relaxation.
- **Suspiciously-clean audit:** the completeness-critic's #1 flag — "exactly 21 QWC groups for
  ALL composed molecules" (LiH Q30 … F₂ Q100) — VERIFIED real, not a bug: reproduced live
  (21 for LiH/H₂O/BeH₂/N₂/F₂/CH₄) and structurally explained (block-diagonal composed ⇒
  different-block terms are qubitwise-commuting ⇒ QWC count is per-block-*type*, size-independent).
- **Remediation applied this pass:** `claim_test_matrix.md` P23 rows synced to 688/828/383.7/
  511.8/507.2 + composed-710 (code-P23 finding); **C18 gate-miss fixed** — it returned green while
  "one to several weeks of focused work" (P23) was live, because the `units-of-project-work`
  pattern required "of" *immediately* before "work"; pattern widened to allow an intervening
  adjective, discrimination proven both ways (fires on the live phrase, silent after the text
  → "a scoped follow-up sprint"), selftest 17/12; lowsu2024 SIAM venue-year 2024→2026;
  `eq:dirac_fs` = Z⁴α² confirmed (former error site intact).

**Honest ceilings + NIT log (non-cert-blocking, fix-on-sight-deferred):** C3 inline provenance
tiers are conveyed by adjectives ("measured"/"exact"/"observation"), not bracketed [MEASURED]/
[SYMBOLIC] tags — the *substance* (prose ≤ backing) was verified but the FORM is thin (standing
characteristic accepted at the v4.63 cert; report as thin-surface). P14 `M = n_{\max}^2` is
wrong (correct: M = Σ_{n=1}^{n_max} n² = 5/14/30/55, which the tables use) and cascades into the
ρ_ERI order-of-magnitude estimate — deferred to a careful focused edit (headline exponents are
fit from data, unaffected). Wording NITs (retired-rule vintage): P20 O(Q^2.5)-headline vs table-2.2 band, "exactly
linear" vs affine 11.10Q+1, 37-systems vs 33-tabulated, He 0.55%→0.544%, H₂O 1-norm 359 vs 361,
LiH 1-norm convention variants, cross-paper 1-norm/Pauli convention consistency; citation:
Chawla/rocca truncated author lists, PachuckiYerokhin2010 key label, NIST version drift; the
flat-21 QWC deserves a one-line structural note in the paper. **PI-scoped process finding:** the
`/regression` workflow never passes `--slow`, so ALL of Paper 14's headline claims (3.15/1.69/
3.36/2.5/2.51×/MPO/TC — all `@pytest.mark.slow`) have no day-to-day CI protection; they pass at
cert time (via `/qa` + a manual `--slow` run, 35/35) but a Gaunt-sparsity regression would slip
`/regression full`. Seed key `debug/qa/group4_fullcert_seed_key.json`; run notes in this record.

---

## Post-certification touch: the group6 D-HFS convention arc (2026-08-28/29)

**Status: certification stands; the touched content is re-verified, not re-certified.**

The `/qa group6` FULL runs corrected a Fermi-contact convention defect whose root
cause (`buggy = correct x 2(m_p/m_N)`) also lived in **Paper 23**.  `eq:d-hfs-bf`
and its four derived values were corrected here as a consequence:

| quantity | was | now |
|---|---|---|
| the equation | `g_d^atomic * (m_e/m_d)` | `(g_d^atomic / 2) * (m_e/m_p)` |
| D 1S BF baseline | 327.3975 MHz | **327.2350 MHz** |
| BF-strict residual | +40 ppm | **-456 ppm** |
| full-chain value | 327.4779 MHz | **327.3153 MHz** |
| full-chain residual | +286 ppm | **-211 ppm** |

Mechanism, in the paper's own correction note: the divisor restoring
`mu_I/(I mu_N)` is **2, not 2I** (they coincide only at I=1), and the mass factor
is `m_e/m_p` for **every** nucleus because `mu_N = e hbar / 2m_p` is defined with
the proton mass.  The two errors nearly cancelled at I=1 (`m_d/m_p = 1.99900750`),
leaving the baseline high by 496.5 ppm; I=1/2 nuclei were unaffected, which is why
the H and T cross-checks never exposed it.

**Re-verification performed 2026-08-29** (not a re-certification -- no seeded panel):
Paper 23 compiles three-pass clean with zero undefined references; no
pre-correction value survives anywhere in the paper (swept for 327.3975 / 327.4779
/ +40 ppm / +286 ppm; the only `2I` occurrences are the legitimate CG factor and
the correction note itself); the corrected values agree with Paper 34's
independently corrected autopsy (-456 BF strict, -211 full chain); C13/C16/C19 and
the inline-arXiv gate all PASS at `--gate group4`; `test_paper23_magic_gaps.py` and
`test_paper23_resource_counts.py` pass; and no test in the suite pins a retired
value (the five corpus-wide matches are correction notes describing the
retirement).

**Not re-opened:** the group4 headline claims (Pauli counts, scaling exponents,
tapering) are untouched by this arc.  The standing `--slow` CI gap recorded above
is unchanged.

## Post-certification touch — 2026-08-29 (ERI evaluator correction, PI-directed)

The wrong-sign Gaunt argument (`q = mc - ma`) behind the "pair-diagonal rule
A" convention (CF-1) was identified as a sign error and removed from all
production evaluators by PI direction ("exact rule everywhere, re-price
everything honestly").  **CF-1 is DISSOLVED.**  Every certified group4
resource number was re-measured, not scaled:

- N_Pauli = 27.90 x Q (main-group) / 30.03 x Q (d-block), exact, replacing
  11.10/9.23; the d-block sparsity-ordering claim REVERSES.
- Within-molecule exponent 2.5 -> 3.17 universal (exact c(n_max) x Q
  factorization; local slope ~3.8 by n_max=4); the "constant-factor A->B"
  claim was a single-point artifact.
- LiH 334 -> 838; 190x -> 76x (cc-pVDZ); 51-1712x -> 54-317x (H2O
  equal-qubit); atomic He advantage INVERTS at Q=10 (288 vs 156).
- Papers 14 + 20 + the group4 synthesis re-priced in place with correction
  records; tables regenerated from live builds; figures regenerated.
- C17 family `composed-rule-a-retired-figures` (discrimination proven both
  ways) guards every retired figure.
- Vintage-marked pending re-measurement: TC-composed table, VQE experiment
  tables, O(Q^1.69) simulation-cost exponent, balanced QWC counts, He
  n_max=5 atomic Pauli count, Chawla matched-Q projection.

**Certification status: the 2026-07-02 certification predates these
corrections; the branch is OWED a re-review** (the same status the v5.0.0
retraction arc created for group3/4/6).  Full record:
debug/sprint_eri_evaluator_defects_memo.md; CHANGELOG v5.1.12.
