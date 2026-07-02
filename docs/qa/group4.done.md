# Group 4 (Quantum computing) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only group4-specific scope + deltas + the branch watch-notes.

> **STATUS: DRAFTED 2026-06-28 for PI freeze** (fifth pre-registered `/qa` target; the
> QC/NISQ/VQE branch). Inherits criteria.md C1–C16. Branch-defining risk = **QC-resource-
> claim honesty + the CF-1 pair-diagonal-ERI disposition** (encoded as a C8/C3/C5/§1.5
> sharpening, not a new number). The §C8 headline numbers below are pending the **CF-1
> decision** (disclose-and-keep vs switch-to-global; see `docs/qa/group4.carryforward.md`
> + the 2026-06-28 sweep) — FREEZE must resolve which set is authoritative.

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

2. **Scaling vs prefactor (the CF-1 partition).** The reviewers must distinguish two tiers
   of claim and check each at its own tier:
   - **Scaling exponents** (atomic $O(Q^{3.15})$, composed $O(Q^{2.5})$, 1-norm $O(Q^{1.69})$,
     QWC $O(Q^{3.36})$, $N_{\rm Pauli}=11.10\,Q$ linearity) — ROBUST under CF-1 (the 2026-06-28
     sweep: the pair-diagonal vs global-$M_L$ re-pricing is a constant factor within valence
     class, so the log-log slope and the linearity are unchanged). These may be stated as MEASURED.
   - **Absolute multipliers / market-test lines** (the "51×–1712× vs Gaussian", the
     "LiH 334 vs STO-3G 907 ≈ 2.7× fewer" market test, the "$d$-orbital blocks are *sparser*
     (9.23 < 11.10)" claim) — CF-1-SENSITIVE. Under the physical global-$M_L$ rule LiH
     re-prices 333→837 (**2.51×**, parity vs STO-3G), and the $d$-block coefficient re-prices
     9.23→**30.0**, becoming DENSER than main-group's 27.9 (the "$d$ sparser" claim REVERSES).
     See `group4.carryforward.md` + `debug/qa/group4_cf1_library_sweep_memo.md`.

3. **CF-1 disclosure — DECIDED A (disclose), applied 2026-06-28; now governed by the shared
   [criteria.md "Dual-rule ERI framing"](criteria.md) rule.** The PI chose **option A**: keep
   the pair-diagonal rule ($m_a=m_c \wedge m_b=m_d$) as the *quality QC sparsity approximation*
   and **disclose** it; **B (global-$M_L$) was deliberately left on the table** (it is the
   physics-accuracy rule, used in the precision-physics paths). Applied: Papers 14 + 20 each
   carry a `\label{sec:eri_rule}` disclosure subsection (the rule, the constant 2.51×/3.25×
   re-pricing, the STO-3G→parity and $d$-block-reversal consequences, B-left-on-table); the
   abstract d-block claim + the 2.7×-market-test line + the TM-table caption now carry the
   pair-diagonal qualifier. **The cert criterion** (now a *framing-zombie* check): every
   sparsity/density/Pauli claim names its rule (A disclosed for sparsity, B for accuracy); an
   **undisclosed pair-diagonal multiplier presented as the exact selection-rule value is
   MATERIAL** — enforced by the `claims-reviewer` (enumerate every such claim) + the
   deterministic **C16 entry `pair-diagonal-as-exact-sparsity`** + the characterization test
   `tests/test_paper14_eri_rule.py` (pins that the product realizes A). The repo study
   confirmed the QC product is **uniformly A** (atomic `lattice_index` + composed
   `composed_qubit`); see `debug/sprint_group4_prework_memo.md`.

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
  *(Numbers marked ⚑ are CF-1-sensitive and pending the FREEZE disposition.)*
  - **Paper 14 (keystone):** atomic **$O(Q^{3.15})$** Pauli (vs Gaussian $Q^{4.25}$ LiH /
    $Q^{3.92}$ H$_2$O, `trenev2025`); QWC groups **$O(Q^{3.36})$**; 1-norm $\lambda$
    **$O(Q^{1.69})$**, $R^2=0.997$ (the key FT result); composed **$O(Q^{2.5})$** universal
    (exponent spread 0.02); $N_{\rm Pauli}=11.10\,Q$ exact (9.23 $d$-block); ERI density
    $\sim 1/M$ (corrected v4.54.0; this watch-note synced 2026-07-01). ⚑ "two-or-more orders of magnitude / 51×–1712× vs Gaussian" + "$d$-block
    sparser"; matched-qubit-not-accuracy caveat MANDATORY.
  - **Paper 16:** $\mu_{\rm free}=\nu(\nu+3N-2)/2$ (SO(3N) Casimir); $\nu=N-2$ universal for
    $S<N/2$; **5 atom types A/B/C/D/E** (FIXED 2026-06-28 — abstract+§IV synced to the Table+
    conclusion+code; was a stale "4"); **"maps onto, rather than predicts"**; Dirac instability =
    metric (not topological) singularity, smooth through $Z=1/\alpha$.
  - **Paper 20:** ⚑ LiH composed **334 Pauli @ 30q vs STO-3G 907 @ 12q, 13× fewer QWC**
    (the market test — re-prices to parity under global-$M_L$); balanced coupled (PK-free)
    binds LiH at **$R_{\rm eq}=3.015$ bohr, 878 Pauli @ 30q, 0.20%** single-point energy at
    the minimum; **row-conditional** chemistry-accuracy (first-row binds; second-row NaH↓
    monotone overattraction — the honest §scope_boundary); library **37 systems** (35 composed
    + He + H2, $Z=1$–56 H–Ba; decided 2026-06-28, this watch-note synced 2026-07-01); $O(Q^{2.5})$ universal vs
    Gaussian $O(Q^{3.9-4.3})$; ⚑ $11.10\,Q$ / $9.23$ $d$-block; frozen cores enter via
    identity only.
  - **Paper 23:** HO closures **2,8,20,40,70,112** from graph state counting; magic
    **2,8,20,28,50,82,126** at $v_{ls}/\hbar\omega\approx0.17$, $d_{\ell\ell}/\hbar\omega
    \approx0.02$; deuteron **16q / 592 Pauli / 1-norm 342 MeV** (FIXED 2026-06-28 from a stale
    227.3, pinned); He-4 **16q / 712 Pauli / 1-norm 467/462 MeV** (FIXED from stale 557/552, pinned)
    (12.25× Hilbert, 1.20× Pauli); composed nuc-electronic deuterium **26q**, ~$10^{13}$
    scale ratio; **Fock rigidity theorem**; binding energies = encoding-validation
    benchmarks (far from experiment, stated as such).
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

1. ✅ **CF-1 pair-diagonal ERI** — **DECIDED A (disclose, PI direction 2026-06-28) + applied;
   codified as a shared QA rule.** Quantified: constant 2.51× main-group / 3.25× $d$-block;
   LiH→parity vs STO-3G; scaling robust. **The repo study confirmed the QC product is uniformly
   A** (atomic + composed pair-diagonal; global-$M_L$ B lives only in the precision-physics paths).
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
(3-pt 1.67, sub-quadratic), QWC $O(Q^{3.36})$ (3-pt 3.356), composed $O(Q^{2.5})$ (~2.5,
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
  QC branch). Inherits criteria.md C1–C16. Branch-defining risk = QC-resource-claim honesty
  + CF-1 pair-diagonal-ERI disposition (C8/C3/C5/§1.5 sharpening, no new number). The CF-1
  re-pricing is quantified (`debug/qa/group4_cf1_library_sweep_memo.md`: constant 2.51×
  main-group / 3.25× $d$-block; scaling robust, prefactor + STO-3G market test + "$d$ sparser"
  re-price). **FREEZE must resolve the CF-1 disposition (A disclose vs B switch) — the §C8
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
