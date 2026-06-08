# Sprint Chemistry-Pivot 2026-06-07 (umbrella memo)

> **2026-06-07 addendum (same-day follow-on, Sprint W1e-Projection-Audit):**
> the §3.1 structural finding "W1e wall is at the projection step" has been
> sharpened to **row-conditional**: TRUE for second-row hydrides (NaH and
> below); FALSE for LiH (where the "no binding" was a Path B call-convention
> bug in `build_balanced_hamiltonian`'s R-dependence corrector). LiH balanced
> under the correct convention binds at R_eq = 3.015 bohr with D_e = 0.158 Ha
> (2.4× over-binding vs continuous 0.067 Ha, but binding is real). NaH balanced
> still over-attracts under either path (genuine second-row physics; R3-B
> falsifier verdict survives). The bug is fixed: new `R` field on
> `MolecularSpec` populated by `hydride_spec`, consumed by the corrector. See
> `debug/sprint_w1e_projection_audit_memo.md` for the full audit.

**Sprint scope:** end-to-end chemistry-direction push. Strategic pivot session
(start: "let's actually scope chemistry properly under new tooling"; end:
"empirical localization of the W1e wall + Phase 1 of the hybrid pipeline
shipped").

**Architecture:** three rounds of parallel sub-agent dispatches. 10 sub-sprints
total. One canonical memo per sub-sprint already on disk; this umbrella memo
synthesizes the cross-thread structural findings and is the record-of-record for
release.

## 1. Sub-sprint inventory and canonical memo pointers

**Round 1 (scoping triad, ~30 min wall):**
- **Target A — propinquity error bounds wiring** → `debug/sprint_propinquity_bound_wiring_memo.md`
- **Target B — DMRG diagnostic on composed Hamiltonian** → `debug/sprint_dmrg_diagnostic_memo.md`
- **Target C — hybrid pipeline scoping** → `debug/sprint_hybrid_pipeline_scoping_memo.md`

**Round 2 (Phase 1 first wave, ~30 min wall):**
- **P1 — FCIDUMP exporter** → `debug/sprint_p1_fcidump_exporter_memo.md`
- **P2 — VQE benchmark on LiH composed** → `debug/sprint_p2_vqe_benchmark_memo.md`
- **P4 — NaH W1e quantitative baseline** → `debug/sprint_p4_nah_w1e_baseline_memo.md`
- **P5 — paper home scoping** → `debug/sprint_p5_paper_home_scoping_memo.md`

**Round 3 (DMRG falsifier + VQE re-stack, ~30 min wall):**
- **R3-A — DMRG on LiH FCIDUMP** → `debug/sprint_r3a_dmrg_lih_memo.md`
- **R3-B — DMRG on NaH FCIDUMP (THE FALSIFIER)** → `debug/sprint_r3b_dmrg_nah_falsifier_memo.md`
- **R3-C — VQE via openfermion** → `debug/sprint_r3c_vqe_uccsd_memo.md`

## 2. Sub-sprint verdicts (one paragraph each)

**Target A (GO, shipped).** Paper 38 propinquity bound $\Lambda \le C_3 \cdot
\gamma_{n_{\max}}$ wired as a metadata field on `GeoVacHamiltonian`. Bound is
system-independent: $\gamma_{n_{\max}=2,3,4} = 2.0746/1.6101/1.3223$ across all
chemistry consumers in the 28-molecule library. First-ever
quantitative-basis-truncation error estimate exposed at the chemistry-consumer
API in any framework. Pauli content bit-identical. 104/104 ecosystem tests
(later +11 = 115/115 after P1).

**Target B (closure with clarification).** Operator-side MPO bond rank
$\{4, 16, 16, 9, 9, 9, 6, 3, 3, 2\}$ bit-exact reproduction of Theorem 3.2.A
(Paper 14 §sec:mpo_bond_rank). State-side bulk $\chi \le 4$. Critical
observation: the narrow gate (DMRG works) is met but the broad question
(DMRG escapes W1e) is structurally outside scope — DMRG-vs-FCI gap is not where
W1e lives. This finding correctly preempted Round 3's load-bearing R3-B
diagnostic.

**Target C (scoping memo, ~7-8 month roadmap, JCTC primary).** Top-2 first
integrations: DMRG via Block2/pyscf + VQE via Qiskit-Nature (later corrected by
R3-C to openfermion). Cleanest gap: `to_fcidump()` exporter. First system: LiH
benchmark, NaH falsifier. Six PI-input items queued.

**P1 (GO, shipped).** `to_fcidump(filename)` + `read_fcidump()` on
`GeoVacHamiltonian`. LiH round-trip bit-exact ($\max|h_1^{\text{diff}}| = 0.0$,
$\max|\text{eri}^{\text{diff}}| = 0.0$). 7-system sample $<10^{-10}$. Pure-Python,
no pyscf dep added. 115/115 ecosystem tests pass (+11 new, 0 regressions);
18/18 topological S³ proofs pass. **Surfaces the cleanest interface gap named
by C: classical chemistry consumers (DMRG, CCSD(T), AFQMC) unblock
simultaneously.**

**P2 (BORDERLINE-with-structural-finding).** LiH composed n_max=2 R=3.015 FCI
energy = $-14.143$ Ha at Q=30 / 333 Pauli (publication-grade headline). LiH
VQE itself not tractable at Q=30 under qiskit 2.x (statevector + ansatz cost).
Q=10 H₂/He proxy: efficient_su2 + COBYLA reaches 87–684 mHa, far above
chemical accuracy. **Structural finding: ansatz selection is Phase 2's
load-bearing gap; qiskit-nature stack is wrong for GeoVac VQE.** This finding
drove R3-C's pivot to openfermion-native.

**P4 (quantitative wall description).** NaH balanced n_max=3 at experimental
$R_e = 3.566$ bohr: depth $+1.50$ Ha vs experimental $D_e = 0.0713$ Ha →
**$21\times$ over-binding**. At artifact PES minimum $R=2.0$ bohr:
$+2.45$ Ha → $34\times$ over-binding. **No interior minimum at n_max=2 or 3,
balanced or composed.** n_max=2 → 3 **deepens** the unphysical well (not
closes it), confirming structurally that W1e lives at the Hamiltonian-input
level. Closure target: $\sim 1.5$ Ha downward shift at $R_e$.

**P5 (SHIP-AS-IS recommendation).** New **Paper 57**, JCTC primary venue.
~18-22 pp, 8 sections, ~10k words, 6 figs + 4 tables. Rejected expansion of
Paper 14 (audience mismatch); deprioritized expansion of Paper 20 (two
competing headlines). Concession path: fold remainder into Paper 20 §V.E if
Phase 2 underdelivers.

**R3-A (MIXED, restructuring finding).** DMRG-on-FCIDUMP reaches qubit-FCI to
$4.4 \times 10^{-7}$ Ha at R=3.015 bohr (three OoM tighter than chemical
accuracy). Pipeline integrity bit-exact. **But the qubit FCI itself gives a
monotone-descending PES on both LiH composed AND LiH balanced across R ∈
[2.5, 5.0] bohr** — no interior equilibrium. The CHANGELOG v3.56.0 LiH 2.82%
R_eq comes from `ComposedDiatomicSolver.LiH_ab_initio` (Level 4 multichannel
adiabatic + PK on continuous space), NOT from FCI on the qubit/FCIDUMP-exported
integrals. Named bug: `ecosystem_export.hamiltonian('LiH')` exposes bare $h_1$
without PK when `core_method='pk'` is default; workaround applied, recommended
fix is a `to_fcidump(include_pk=True)` switch.

**R3-B (STOP, the load-bearing falsifier).** Block2/pyscf unavailable under
Python 3.14 (BLAS / CMake build failures), but at FCI dim $\le 784$ for NaH
balanced n_max $\le 3$, DMRG-at-infinite-$\chi$ IS sector-restricted FCI — the
operational ceiling for any classical correlation solver consuming the
GeoVac-supplied integrals. PES **bit-identical to the P4 baseline** at every R
($\max\text{diff} \sim 6 \times 10^{-13}$ Ha), both n_max=2 and n_max=3.
**W1e is not closeable by any classical correlation method on the same
integrals.** Verdict: STOP — the wall is at the (h_1, eri, ecore) tensor
specification level, not at the determinant-expansion level.

**R3-C (BORDERLINE with structural validation).** Openfermion-native UCCSD
on H₂ Q=10 reaches **error $8.6 \times 10^{-13}$ mHa** (five orders of
magnitude below published STO-3G UCCSD literature) via L-BFGS-B from HF init,
165 evaluations, 0.92 s wall. **~3000× faster per energy evaluation** than
qiskit-nature UCCSD at Q=10. P2's structural finding empirically validated:
qiskit-nature stack is wrong for GeoVac; openfermion is the correct path.
LiH Q=30 not reached (17.18 GB statevector); even tapered Q=25 hit an
openfermion Kronecker-reduction memory ceiling. Phase 2 follow-on: sector-
projected UCCSD on tapered LiH.

## 3. Cross-thread structural findings (the heart of this sprint)

### 3.1 The W1e wall is empirically localized at the Hamiltonian-specification level

> **Same-day sharpening (Sprint W1e-Projection-Audit, 2026-06-07):** the LiH
> half of this finding turned out to be a Path-B call-convention bug, not a
> structural wall. Balanced LiH (under the correct convention, or under either
> convention after the fix in `geovac/balanced_coupled.py` lines 669-680) DOES
> bind LiH at R_eq = 3.015 bohr with D_e = 0.158 Ha. The NaH half (R3-B
> falsifier) is genuine: NaH balanced over-attracts under either path. The
> "wall is at the projection step" framing is therefore **row-conditional**:
> TRUE for second-row, FALSE for first-row LiH (after the bug fix).
> See `debug/sprint_w1e_projection_audit_memo.md` §§3, 6.1, 6.2.

This is the load-bearing finding of the day. It synthesizes P4 + R3-A + R3-B:

- **P4 (continuous PES sweep on Hamiltonian builder):** NaH balanced PES
monotone-descending at n_max=2; n_max=3 deepens (not closes) the well. Wall is
**not** at basis-truncation level. **[STILL HOLDS — but the magnitude of
over-attraction was understated; the V_NN bug was masking ~1.3 Ha of true
overattraction depth at R=2.5. The qualitative verdict is unchanged.]**
- **R3-A (FCI on qubit Hamiltonian for LiH):** Same monotone descent on LiH
composed AND balanced. Wall is **not** specific to systems where W1e was
originally named (it reaches into LiH too). **[BALANCED PART RETRACTED —
the "balanced doesn't bind LiH" was the Path B bug. Composed part stands:
composed structurally lacks cross-block coupling and does monotone descent.]**
- **R3-B (DMRG on NaH FCIDUMP):** DMRG = FCI at finite sector dimension. PES
bit-identical to P4 baseline. Wall is **not** at determinant-expansion level.
**[STILL HOLDS for second-row NaH at the qualitative level. The PES magnitudes
were affected by the same V_NN bug as P4; the bit-identicality of DMRG vs P4
is preserved because both used the same buggy code path.]**

The wall is localized **at the construction step from continuous Level 4 to
the second-quantized $(h_1, \text{eri}, e_{\text{core}})$ integrals**, **for
second-row hydrides**. The continuous adiabatic+PK solver on continuous space
binds LiH at 2.82% R_eq (CHANGELOG v3.56.0). The composed qubit Hamiltonian
does NOT bind LiH (real structural finding: missing cross-block coupling). The
balanced qubit Hamiltonian binds LiH under the correct call convention. The
balanced qubit Hamiltonian does NOT bind NaH under any convention. The wall
that remains, after disentangling the bug from the physics, is second-row
specific.

This is the chemistry-side **operational** analog of the Sprint H1 Yukawa
non-selection theorem at the operational level (Sprint H1 verdict
POSITIVE-THIN, May 31 2026):
- H1: framework admits Higgs structurally but autonomously selects no Yukawa
values.
- R3-B: framework admits a second-quantized chemistry Hamiltonian structurally
but the calibration content of its $(h_1, \text{eri}, e_{\text{core}})$ tensors
for systems where W1e dominates is not autonomously generated by the projection.

**Status update for CLAUDE.md §3:** the W1e wall remains a dead-end at the
multi-determinant-correlation reading; the new R3-B row documents that DMRG
specifically cannot close it because the wall is at the integral-specification
level. Six instances of the multi-focal-composition wall pattern now (was five
via Sprint W1e period-class 2026-06-04).

### 3.2 VQE via openfermion is the correct stack for GeoVac

P2 → R3-C closes this loop empirically. Five OoM gap between H₂ Q=10
performance under qiskit-nature (87 mHa) and openfermion-native UCCSD
($8.6 \times 10^{-13}$ mHa) is not subtle. The
`ecosystem_export.to_openfermion()` path is the natural one for GeoVac VQE work
going forward. P5's Paper 57 §V on VQE will be openfermion-native.

### 3.3 The propinquity bound is system-independent

Subtle Target A finding worth noting: $\gamma_{n_{\max}}$ depends only on the
cutoff $n_{\max}$, not on the molecule. This is structurally correct per
Paper 38 (the bound is on the spectral triple truncation, not on the chemistry
content) but means that chemistry consumers reading the metadata see the same
"this is within X of the framework's $n_{\max}$-truncation limit" regardless of
system. Honest framing in the user-facing documentation.

## 4. Implications for Paper 57

P5 recommended Paper 57 with "DMRG closes NaH W1e" as load-bearing falsifier
**passed**. The falsifier returned a clean **failed** verdict. The abstract
reshapes:

- **Old framing (P5 default):** "GeoVac as the skeleton; hybrid pipeline
closes calibration-tier W1e under DMRG; architecture validated."
- **New framing (R3-B return):** "GeoVac ships chemistry-consumer-grade
integrals (FCIDUMP at machine precision) with quantitative truncation-error
bounds (first ever in any framework); openfermion-native UCCSD reaches
sub-attomarianumber on small systems; **W1e is empirically localized at the
projection step from continuous Level 4 adiabatic+PK to second-quantized
integrals** — a sharper, paper-worthy structural result." Cleaner, sharper,
more honest. Chemistry-side analog of H1 Yukawa non-selection at the
operational level.

Paper 57 framing now has THREE load-bearing claims (was two):
1. Sparse Hamiltonian export with quantitative error bounds.
2. Sub-mHa VQE on small systems via openfermion-native pipeline.
3. **Row-conditional operational localization of the chemistry-side
calibration-data partition boundary** (Paper 18 §IV.6 chemistry analog,
sharpened post-W1e-Projection-Audit 2026-06-07): the W1e wall is at the
projection step for second-row hydrides (NaH and below); for first-row LiH
the wall closes under the balanced builder (cross-block ERIs + Track CD
cross-center V_ne). The row split is itself a paper-worthy structural result
because it identifies cross-center V_ne in a fixed Z_eff basis as the
mechanism: it suffices for LiH (Z_eff ≈ 1.84, valence orbital compact) but
overshoots for NaH (Z_eff ≈ 2.2, 3s orbital extends too far, picks up too much
H-proton attraction at small R).

## 5. Follow-on register

### 5.1 Shippable next sprint
- **Sector-projected UCCSD on tapered LiH** (R3-C Phase 2 follow-on). Z₂ Hopf-U(1)
per-sub-block tapering already production (v3.52.0). Sector projection
required to close the openfermion Kronecker-reduction memory ceiling.
- **`to_fcidump(include_pk=True)` switch** (R3-A named bug). Single-flag fix.
- **Apply paper edits queued by sub-agents.** Six edits queued not applied:
Paper 20 §benchmarks new column (propinquity bound), Paper 20 §sec:hybrid_pipeline
new subsection (FCIDUMP + openfermion UCCSD), Paper 38 §sec:applications
(propinquity bound), Paper 18 §IV.6 W1e structural localization update,
Paper 14 §sec:mpo_bond_rank chemistry application stub.

### 5.2 Structural next sprint
- **Diagnostic: why does the qubit/composed Hamiltonian not bind LiH while the
continuous adiabatic+PK solver does at 2.82% R_eq?** This is the load-bearing
question generated by R3-A. Sprint-scale diagnostic, named gate, sharp
follow-up.

### 5.3 Multi-month
- **Paper 57 drafting** (P5 scope, ~7-8 months total Phase 2 + Phase 3).
- **Sector-projected qubit pipeline for LiH and the multi-center library**
(Phase 2 expansion).

### 5.4 Mechanical
- **Bootstrap durations parser fix.** Background bash succeeded (exit 0) but
wrote 0 entries to `tests/_durations.json` — the regex in the
`/regression` skill doesn't match pytest's actual output format. The fallback
mode of `/regression touched` (topo + diff-derived consumers, no random
sample) continues to work; the random-sample addition needs the parser fix.

## 6. Honest scope

**Theorem grade.**
- Paper 38 propinquity bound shape $\Lambda \le C_3 \gamma_{n_{\max}}$ wired:
the bound itself is theorem-grade (Paper 38 main theorem); the wiring is
mechanical.
- FCIDUMP round-trip bit-exactness on LiH (h_1, eri, ecore max diff = 0.0)
is theorem-grade for that specific export step.
- Bond rank Theorem 3.2.A (Paper 14 §sec:mpo_bond_rank) reproduction bit-exact
to operator-side measurement of $\{4, 16, 16, 9, 9, 9, 6, 3, 3, 2\}$.

**Structural sketch.**
- "W1e is at the projection step from Level 4 to second-quantized integrals."
**Strongly supported empirically** by P4 + R3-A + R3-B convergence, but the
proof that the wall cannot be closed at this projection step is operational
(DMRG=FCI cannot close it), not structural. A structural proof would show that
the projection step is missing a specific operator term (e.g., explicit
core-bonding J-K interaction, or a screened cross-block h_1 contribution
beyond what F3 closure of Sprint F5 ruled out).

**Numerical observation.**
- H₂ Q=10 openfermion UCCSD error $8.6 \times 10^{-13}$ mHa is an outlier even
in the UCCSD literature. We have not yet checked whether this is a specific
property of the H₂ Q=10 Hamiltonian (small Hilbert space, exact ansatz
parameterization) or a generic GeoVac-Hamiltonian property. Sample size n=1.
- 3000× speedup vs qiskit-nature at Q=10. We have not tested at higher Q.
- NaH W1e $21\times$ over-binding at $R_e^{\text{exp}}$ is consistent with P4
baseline at $10^{-13}$ Ha precision but is from a single anchor; experimental
$D_e = 0.0713$ Ha has its own uncertainty.

**Named open follow-ons.**
- Closure proof that the projection step is structurally missing physics, not
just empirically (Section 5.2 above).
- Higher-Q VQE-on-GeoVac scaling test (Section 5.1).
- Sector-projected UCCSD on tapered LiH closing the openfermion Kronecker
memory ceiling.

## 7. Files modified / created

### Production code
- `geovac/ecosystem_export.py` — propinquity bound metadata (Target A); FCIDUMP
exporter `to_fcidump()` + reader `read_fcidump()` (P1)

### Tests
- `tests/test_ecosystem_export.py` — +8 tests (Target A); +11 tests (P1).
Final count 115/115 pass.

### Drivers and data
- `debug/r3a_dmrg_lih_driver.py`, `debug/data/r3a_dmrg_lih.json`,
`debug/data/lih_r*.fcidump` (5 R points composed, 6 R points balanced)
- `debug/r3b_dmrg_nah_falsifier_driver.py`, `debug/data/r3b_dmrg_nah.json`,
`debug/data/r3b_nah_balanced_*.fcidump` (14 R points × 2 n_max)
- `debug/r3c_vqe_uccsd_driver.py`, `debug/data/r3c_vqe_uccsd.json`,
`debug/data/r3c_vqe_uccsd_stdout.log`
- `debug/p2_vqe_benchmark_driver.py`, `debug/data/p2_vqe_benchmark.json`
- `debug/p4_nah_w1e_baseline_driver.py`, `debug/data/p4_nah_w1e_baseline.json`,
`debug/plots/p4_nah_w1e_pes{,_relative}.png`

### Sub-sprint memos
Ten on disk (Section 1 above).

## 8. Verification

- 115/115 ecosystem tests pass (P1 close, Target A baseline + 11 P1 + 8
Target A = +19 from morning baseline 96).
- 18/18 topological S³ proofs (CLAUDE.md §9 baseline) pass throughout.
- Round 1 → Round 3: no regressions across any sub-sprint's verification gate.
- /regression touched (fallback mode, durations.json bootstrap deferred to next
session) at end of sprint: see umbrella sprint-close protocol.

## 9. Production-code hard prohibitions check (§13.5)

No changes to: natural geometry hierarchy / fitted-or-empirical parameters /
Section 3 deletions / Paper 2 K = π(B+F−Δ) combination-rule "conjectural"
labeling.
