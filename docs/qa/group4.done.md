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
    $\sim 1/M^2$. ⚑ "two-or-more orders of magnitude / 51×–1712× vs Gaussian" + "$d$-block
    sparser"; matched-qubit-not-accuracy caveat MANDATORY.
  - **Paper 16:** $\mu_{\rm free}=\nu(\nu+3N-2)/2$ (SO(3N) Casimir); $\nu=N-2$ universal for
    $S<N/2$; **5 atom types A/B/C/D/E** (FIXED 2026-06-28 — abstract+§IV synced to the Table+
    conclusion+code; was a stale "4"); **"maps onto, rather than predicts"**; Dirac instability =
    metric (not topological) singularity, smooth through $Z=1/\alpha$.
  - **Paper 20:** ⚑ LiH composed **334 Pauli @ 30q vs STO-3G 907 @ 12q, 13× fewer QWC**
    (the market test — re-prices to parity under global-$M_L$); balanced coupled (PK-free)
    binds LiH at **$R_{\rm eq}=3.015$ bohr, 878 Pauli @ 30q, 0.20%** single-point energy at
    the minimum; **row-conditional** chemistry-accuracy (first-row binds; second-row NaH↓
    monotone overattraction — the honest §scope_boundary); library **38 molecules** (Z=1–36;
    ⚑ reconcile vs Paper 14 "30" / CLAUDE.md "40" / 35 shipping); $O(Q^{2.5})$ universal vs
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
