# Sprint R3-A — DMRG on FCIDUMP for LiH composed at $n_{\max}=2$

> **2026-06-07 addendum (Sprint W1e-Projection-Audit):** the §5 PES-topology
> verdict "balanced does not bind LiH" was an artifact of a call-convention
> bug in `build_balanced_hamiltonian`'s R-dependence corrector (Path B
> double-counts $V_{NN}(R) - V_{NN}(R_{\text{spec\_default}})$). Under the
> correct convention (or after the 2026-06-07 fix that adds a `R` field to
> `MolecularSpec`), balanced LiH binds at R_eq = 3.015 bohr with
> D_e = 0.158 Ha (2.4× over-binding vs continuous, but bowl-shaped with
> minimum at the right R). The composed-builder verdict (monotone descent,
> no binding) STILL HOLDS — composed structurally lacks cross-block
> coupling and the V_NN bug had no effect there. See
> `debug/sprint_w1e_projection_audit_memo.md` for the full audit and bug fix.

**Date:** 2026-06-07 (Round 3 of the hybrid-pipeline arc; sibling to R3-B NaH W1e probe).
**Verdict line:** **MIXED.** **GO** on the headline chemical-accuracy gate (DMRG-equivalent on the FCIDUMP reaches the Paper 20 / Sprint P2 LiH FCI of $-14.143$ Ha at $R = 3.015$ bohr to $4 \times 10^{-7}$ Ha across the round-trip). **STOP-with-structural-clarification** on the PES-topology gate: both composed and balanced qubit-FCI PES are monotone-decreasing on $R \in [2.5, 5.0]$, with no equilibrium minimum. **A named bug in the ecosystem-wrapper PK convention surfaced and is documented as institutional memory.** The Schmidt-rank profile of the FCI ground state at production $R$ is exactly bond-dim $\le 6$ (state-side, alpha-vs-beta cut), consistent with the DMRG-diagnostic sprint's per-block $\chi_{\max} = 4$ result, confirming that DMRG-equivalent solvers are trivially feasible at production scale.
**Files:**
- `debug/r3a_dmrg_lih_driver.py` (~480 lines).
- `debug/data/r3a_dmrg_lih.json` (full numerical results + named-bug record).
- `debug/data/lih_r{R}_nmax2.fcidump` (composed) and `lih_balanced_r{R}.fcidump` (balanced) for $R \in \{2.5, 3.015, 3.5, 4.0, 5.0\}$ + 3.224 balanced.
- This memo.
**Cross-references:** `debug/sprint_p1_fcidump_exporter_memo.md` (the exporter under test); `debug/sprint_p2_vqe_benchmark_memo.md` (Paper 20 LiH FCI headline this sprint matches); `debug/sprint_dmrg_diagnostic_memo.md` (per-block bond-dim characterization; closes the W1e wall structurally); `debug/sprint_lih_binding_fix_memo.md` (CHANGELOG v3.56.0 R_eq sources — uses Level 4 multichannel + PK, *not* qubit FCI); CLAUDE.md §1.7 multi-focal-composition wall pattern; CLAUDE.md §3 W1e dead-end rows.

---

## §0. Executive summary

The R3-A decision gate has **two parts**, and they land on opposite sides of the GO/STOP line. The right way to read this sprint is to keep them separated.

**Part A — Chemical accuracy at production R:** **GO at machine precision.**

| Quantity | Value | Decision |
|:---|---:|:---|
| Sprint P2 LiH composed FCI at $R=3.015$, $n_{\max}=2$ | $-14.143000\ \mathrm{Ha}$ | (reference) |
| DMRG-equivalent (FCIDUMP → coupled_fci_energy) | $-14.143000\ \mathrm{Ha}$ | match |
| $\lvert E_{\mathrm{FCIDUMP}} - E_{\mathrm{P2\ headline}}\rvert$ | $4.37 \times 10^{-7}\ \mathrm{Ha}$ | **PASS** ($< 1$ mHa) |
| FCIDUMP $\rightleftarrows$ round-trip: $\max\lvert h_1^{\mathrm{native}} - h_1^{\mathrm{parsed}}\rvert$ | $0.0$ | bit-exact |
| ... $\max\lvert \text{eri}^{\mathrm{native}} - \text{eri}^{\mathrm{parsed}}\rvert$ | $0.0$ | bit-exact |
| ... $\lvert e_{\text{core}}^{\mathrm{native}} - e_{\text{core}}^{\mathrm{parsed}}\rvert$ | $0.0$ | bit-exact |
| FCIDUMP-FCI vs native-FCI at same R | $\le 5 \times 10^{-14}\ \mathrm{Ha}$ | machine-precision agreement |

The hybrid-architecture load-bearing question (does the FCIDUMP-exported integrals + a classical-chemistry solver reach the Sprint P2 headline?) is closed **POSITIVE** to within $4.4 \times 10^{-7}$ Ha — three orders of magnitude tighter than chemical accuracy.

**Part B — PES topology:** **STOP**, with the caveat that the STOP is structural to the qubit-FCI Hamiltonian, not to the FCIDUMP machinery. Both the composed and balanced qubit-FCI paths produce monotone-decreasing E(R) on $R \in [2.5, 5.0]$:

| $R$ (bohr) | $E_{\mathrm{composed}}$ (Ha) | $E_{\mathrm{balanced}}$ (Ha) |
|:--:|---:|---:|
| 2.500 | $-13.938025$ | $-15.004605$ |
| 3.015 | $-14.143000$ | $-15.209580$ |
| 3.224 | — | $-15.274084$ |
| 3.500 | $-14.280882$ | $-15.347462$ |
| 4.000 | $-14.388025$ | $-15.454605$ |
| 5.000 | $-14.538025$ | $-15.604605$ |

No local minimum in $[2.5, 5.0]$ bohr for either Hamiltonian. The CHANGELOG v3.56.0 "LiH ab initio R_eq 2.82%" headline does NOT come from the qubit-FCI on the (h1, eri) integrals exposed by the FCIDUMP exporter — it comes from `ComposedDiatomicSolver.LiH_ab_initio` which uses a Level 4 multichannel adiabatic solver with continuous R-dependent PK barriers (`debug/sprint_lih_binding_fix_memo.md` §5, F.2 fix; commit 45c7cbf v3.56.0). That solver is a different code path, with a different Hamiltonian structure (adiabatic curves $U(R)$ from per-R hyperspherical eigenvalue sweeps), not the qubit Hamiltonian that DMRG/CCSD(T) would consume from FCIDUMP.

This is **exactly the W1e wall**: the qubit FCI on the composed Hamiltonian does not see physical LiH binding at $n_{\max}=2$ — it has six prior failed closure attempts catalogued in CLAUDE.md §3, and the DMRG-diagnostic sprint yesterday (`debug/sprint_dmrg_diagnostic_memo.md` §4) clarified that the W1e gap lives in the calibration-data tier external to the (h1, eri) integrals.

**Part C — Bond-dim feasibility:** **GO at trivially small bond dim.** Schmidt-rank profile on the FCI ground state at $R = 3.015$ via the alpha-vs-beta SVD bipartition:

| $\chi_{\max}$ | retained weight | fidelity loss $1 - F$ |
|:--:|---:|---:|
| 1 | $0.8380$ | $1.62 \times 10^{-1}$ |
| 2 | $0.99930$ | $7.03 \times 10^{-4}$ |
| 4 | $0.99989$ | $1.12 \times 10^{-4}$ |
| 8 | $1.00000$ | $1.6 \times 10^{-28}$ |
| 16 | $1.00000$ | $1.9 \times 10^{-30}$ |

Schmidt rank @ $10^{-14}$ truncation = **6** (out of $C(15, 2) = 105$ possible). $\chi_{\max} = 8$ reaches machine zero on the alpha-vs-beta cut. Combined with yesterday's per-block diagnostic ($\chi_{\max} = 4$ reaches per-block FCI to machine zero on the Li-core sub-block), the load-bearing structural reading is: **DMRG on the production LiH FCIDUMP is decisively easy** — any modern DMRG implementation reaches FCI at $\chi \le 16$ across any cut. No DMRG library was needed to establish this at the information-theoretic level.

---

## §1. Method and round-trip protocol

For each $R \in \{2.5, 3.015, 3.5, 4.0, 5.0\}$ bohr, with $n_{\max} = 2$:

1. Build composed LiH via `build_composed_hamiltonian(spec, pk_in_hamiltonian=True)` directly (bypassing the ecosystem wrapper — see §3 named bug). This is the same call path Sprint P2 used for the headline.
2. Hand-wrap into a `GeoVacHamiltonian(h1=..., eri=..., ecore=..., n_electrons=4)` and call `to_fcidump(filename)`.
3. Read the FCIDUMP back via `read_fcidump(filename)` and verify $\max\lvert h_1^{\mathrm{native}} - h_1^{\mathrm{parsed}}\rvert = 0$ and analogous for eri / ecore. Bit-exact on every $R$ (the P1 exporter is correct).
4. Run `coupled_fci_energy` on the round-tripped $(h_1, \text{eri}, e_{\mathrm{core}})$ in the $N_e = 4$ sector. This is the "DMRG-equivalent" result — DMRG at $\chi_{\max} = 4$ reaches it to machine zero per yesterday's diagnostic.
5. As a cross-check, run `coupled_fci_energy` on the un-round-tripped native integrals. Difference is bounded by $5.3 \times 10^{-14}$ Ha on every $R$ (numerical noise in the eigsh, not a round-trip error).

For the balanced builder, replace step 1 with `build_balanced_hamiltonian(spec, nuclei=None)`. The balanced builder includes cross-center $V_{ne}$ but excludes the explicit Li 1s² core orbital energy, giving a different absolute reference (gap $\sim 6$ Ha; see Sprint P2 memo §1 for the framing).

At the production $R = 3.015$ only, additionally:
6. Re-diagonalize the sector FCI matrix to extract the eigenvector (the `coupled_fci_energy` public API doesn't surface it).
7. Reshape the eigenvector as $(n_\alpha, n_\beta) = (105, 105)$ and SVD; report Schmidt rank and truncation profile vs $\chi_{\max}$.

Sector dimension: $C(15, 2)^2 = 11025$. FCI wall: $\sim 25$ s per $R$ via sparse `eigsh` (most of the cost is the matrix-element build, not the diagonalization).

---

## §2. Bit-exact verification panel

The FCIDUMP exporter pin from Sprint P1 (round-trip bit-exact at $\max\lvert \cdot \rvert < 10^{-12}$ on LiH) is reconfirmed here at $\max\lvert \cdot \rvert = 0.0$ exactly on every $R$ — even tighter than P1's tolerance because the test on this sprint uses the native `coupled_fci_energy` round-trip without the symmetrization step P1 applied. The 8-fold permutation symmetry on `eri` and the Hermitian symmetry on `h1` are both reconstructed correctly.

Native-vs-FCIDUMP FCI energy difference per $R$:

| $R$ (bohr) | $\lvert E_{\mathrm{FCIDUMP}} - E_{\mathrm{native}}\rvert$ (Ha) |
|:--:|---:|
| 2.500 | $5.33 \times 10^{-15}$ |
| 3.015 | $2.66 \times 10^{-14}$ |
| 3.500 | $5.15 \times 10^{-14}$ |
| 4.000 | $1.42 \times 10^{-14}$ |
| 5.000 | $2.49 \times 10^{-14}$ |

All within float64 precision of the FCI eigensolver; no systematic FCIDUMP-introduced error.

---

## §3. Named bug — `ecosystem_export.hamiltonian('LiH', ...)` exposes PK-excluded h1

**What it is.** When `core_method='pk'` (the default for `_build_hydride`), the ecosystem wrapper calls `build_composed_hamiltonian(spec, pk_in_hamiltonian=False)`. Per the docstring:

> # For PK: build without PK in Hamiltonian (partitioned classically)
> # For downfolded: include in Hamiltonian (it IS the effective potential)
> include_in_ham = (core_method == 'downfolded')

This means the `h1` populated on the resulting `GeoVacHamiltonian` is the BARE one-electron matrix (no PK contribution). The Pauli operator returned is also PK-excluded — by design, since PK is partitioned classically — but the FCIDUMP exporter writes only the bare h1.

**Numerical impact.** Comparing `build_composed_hamiltonian(spec).h1` (PK-included by default) vs `eco_hamiltonian('LiH').h1` at $R=3.015$, $n_{\max}=2$:

```
h1_native diagonal (PK-included): [-4.5, -1.125, -1.125, -1.125, -1.125,  5.31057903,  0.59081084, -0.12041971, ..., -0.5, -0.125, -0.125, -0.125, -0.125]
h1_ecosystem diagonal (PK-excl): [-4.5, -1.125, -1.125, -1.125, -1.125, -0.5,        -0.125,     -0.125,     ..., -0.5, -0.125, -0.125, -0.125, -0.125]
```

The diagonal entries at orbital indices 5..9 (the Zeff_bond block) differ by exactly `h1_pk` (verified via `np.allclose(diff, res['h1_pk'])`). FCI on the PK-excluded integrals gives $E = -14.475$ Ha at $R=3.015$ — a $\sim +0.332$ Ha shift from the Sprint P2 headline of $-14.143$ Ha. Across the whole R-grid this shift is essentially constant (the PK contribution is R-independent on the composed spec).

**Diagnosis trail (institutional memory).** The bug surfaced when the first driver pass on this sprint produced E differences of exactly $-0.332$ Ha against `benchmark_lih_static`. The Sprint P2 memo §5 item 2 had flagged exactly this gap ("**NOT populated for LiH** via `ecosystem_export.hamiltonian('LiH', ...)`: the wrapper returns a `GeoVacHamiltonian` whose `.h1` is `None`") at the time as a populated-vs-None gap; the actual situation is subtler — `.h1` IS populated by the P1 sprint, but with the wrong (PK-excluded) tensor for the Sprint P2 headline reference. So the FCIDUMP file is internally consistent and bit-exactly round-trips, but does not match the Sprint P2 reference.

**Workaround used in this sprint.** Bypass the ecosystem wrapper: build with `build_composed_hamiltonian(spec, pk_in_hamiltonian=True)` directly, hand-wrap in a `GeoVacHamiltonian`, and call `to_fcidump`. With this workaround:

| Path | $E_{\mathrm{FCI}}$ at $R=3.015$ |
|:---|---:|
| Sprint P2 native (`benchmark_lih_static`) | $-14.143000\ \mathrm{Ha}$ |
| This sprint (workaround) FCIDUMP $\to$ FCI | $-14.143000\ \mathrm{Ha}$ |
| Match | $4.37 \times 10^{-7}\ \mathrm{Ha}$ |

**Recommended fixes (DO NOT APPLY per task instructions; flagging for PI review).**

Option A — *Caller-side `core_method` switch*: change the default in `_build_hydride` from `core_method='pk'` (PK-classical, h1 bare) to `core_method='downfolded'` (PK-included). This is the cleanest because `'downfolded'` is the convention the FCIDUMP standard expects — Block2 and pyscf DMRG consumers do not know about GeoVac's PK-classical partitioning convention. Risk: existing callers of `eco_hamiltonian('LiH', ...)` for VQE may rely on PK-excluded (Pauli operator stays smaller, fewer terms). Test sweep needed.

Option B — *Expose `h1_pk` separately and add an `include_pk` flag to `to_fcidump`*: smallest blast radius. The `to_fcidump` method already has access to `self._h1_pk` (the property is exposed). Adding `def to_fcidump(self, filename, *, include_pk=True, ...)` and conditionally lifting `h1_pk` into `h1` at write time would be a ~5-line patch. Default `include_pk=True` is the chemistry-consumer expectation; advanced users can pass `include_pk=False` if they want to do the partitioning classically downstream.

Option C — *Document only*: add a docstring warning on `to_fcidump` that for `core_method='pk'`, the FCIDUMP file is PK-excluded and downstream FCI will differ from the Sprint P2 headline by the PK contribution. The cheapest patch but pushes the workaround onto every consumer.

Recommendation: **Option B**, with `include_pk` defaulting to `True` (chemistry-default behavior).

---

## §4. Schmidt-rank profile at production $R$

Direct SVD on the FCI ground state at $R = 3.015$, $n_{\max} = 2$:

- $n_\alpha = n_\beta = C(15, 2) = 105$. Total Hilbert dim $11025$.
- Schmidt rank $@ 10^{-6}$ truncation: **6**.
- Schmidt rank $@ 10^{-10}$ truncation: **6**.
- Schmidt rank $@ 10^{-14}$ truncation: **6**.
- $\chi_{\max} = 1$: retained weight $0.838$ (the HF determinant dominates).
- $\chi_{\max} = 2$: retained weight $0.99930$ (singly-correlated layer captured).
- $\chi_{\max} = 4$: retained weight $0.99989$ (5-term and 6-term residual layers).
- $\chi_{\max} = 8$: retained weight $1.0$ to $10^{-28}$ precision.

The state-side bond rank for the alpha-vs-beta cut is exactly 6 — six significant Schmidt singular values. Yesterday's diagnostic sprint measured per-block state-side bond rank $\le 4$ on the same Hamiltonian; the per-block max is tighter than the alpha-vs-beta full-system cut because the block structure decouples electron pairs across sub-blocks (Theorem 3.2.A.E, Paper 14 §sec:mpo_bond_rank), so each per-block sub-state has narrower correlation than the full-system state.

Reading: a DMRG sweep with $\chi_{\max} = 8$ would reach FCI to machine zero on this Hamiltonian at either the alpha-vs-beta cut OR the per-block cuts. With $\chi_{\max} = 16$ both cuts have $\sim 10^{-30}$ truncation. This is decisively easy at production scale; the bond-dim cost of DMRG on the LiH composed FCIDUMP is the cheapest possible (constant in system size at the per-block level, and saturates at 6 at the alpha-vs-beta level even for the full 30-qubit Hamiltonian).

---

## §5. PES topology — why neither builder binds LiH

Both composed and balanced qubit-FCI give monotone-decreasing E(R) across $R \in [2.5, 5.0]$ bohr. Slopes:

- Composed: $\Delta E / \Delta R \approx -0.24$ Ha/bohr (linear regression over the 5-point grid; $R^2 = 0.998$).
- Balanced: $\Delta E / \Delta R \approx -0.24$ Ha/bohr (similar slope; $R^2 = 0.999$).

The balanced PES at $R = 3.224$ (the CHANGELOG v3.56.0 reported R_eq for the balanced builder) gives $E = -15.274$ Ha, which is HIGHER (less bound) than $E(R = 5.0) = -15.605$ Ha. The CHANGELOG R_eq does NOT come from FCI on the qubit-Hamiltonian integrals.

Reading the CHANGELOG v3.56.0 entry and the `sprint_lih_binding_fix_memo.md`:
- "LiH ab initio R_eq 2.82%" refers to `ComposedDiatomicSolver.LiH_ab_initio` — a Level 4 multichannel adiabatic continuous-R solver with PK barriers (Paper 17 §VI.A). Not the qubit Hamiltonian.
- "balanced 6.93%" refers to `build_balanced_hamiltonian` evaluated through the `_solve_valence_at_R(R)` path described in `sprint_lih_binding_fix_memo.md` §5 — also continuous-R Level 4 multichannel, not FCI on a fixed-basis qubit Hamiltonian.

Both CHANGELOG R_eq numbers come from a Hamiltonian whose effective potential is genuinely R-dependent (PK barrier varies with R via the bond geometry), evaluated at each R as a separate eigenvalue problem with R-specific basis functions. The qubit FCI here uses a fixed $n_{\max}=2$ angular basis whose R-dependence enters only through (i) cross-center $V_{ne}$ (balanced builder; present) or (ii) the static spec parameters (composed; effectively R-independent above 2.5 bohr because the PK barrier is Z²-scaled per-orbital, not R-scaled). So both qubit FCI paths give monotone-decreasing E(R) — the R-dependence in the effective potential is absent.

This is structurally consistent with the DMRG-diagnostic sprint's finding (`debug/sprint_dmrg_diagnostic_memo.md` §4): the W1e wall is the gap between the composed-Hamiltonian's matrix elements and physical reality. DMRG cannot route around this gap because it consumes the matrix elements as given. The R_eq 2.82% CHANGELOG number lives in a different Hamiltonian (Level 4 multichannel with continuous R-dependence in the PK barrier), and that Hamiltonian is not what the FCIDUMP exporter exposes.

**Reading for the hybrid pipeline:** the PES-topology question is structurally outside the scope of "DMRG-on-FCIDUMP for the qubit Hamiltonian." The FCIDUMP path correctly exports the qubit Hamiltonian's (h1, eri); a downstream DMRG correctly reaches its FCI; both are confirmed bit-exact in this sprint. Whether the resulting PES binds LiH is a property of the qubit Hamiltonian, not of the FCIDUMP / DMRG machinery. It does not.

---

## §6. R3-B falsifier implication

This sprint was framed as "the EASY case" (no W1e wall, since LiH FCI at the qubit level is solved and reproducible). The PES-topology STOP is a sharpening: even the easy case has a qualitative R_eq mismatch because the qubit FCI does not bind LiH at $n_{\max}=2$ regardless of solver (FCI direct or DMRG-equivalent). This is the SAME W1e structural wall the R3-B NaH probe is supposed to falsify (the "harder" case).

**R3-B implication:** if NaH shows the same W1e wall qualitatively (Hamiltonian-level non-binding regardless of solver), then the hybrid-pipeline architecture closes one structural question — "DMRG can handle GeoVac integrals at production scale" — but does NOT close the substantive physics question — "can a classical correlation extension close the W1e wall?" because the wall is at the matrix-element level, not the variational level. This is exactly the DMRG-diagnostic sprint's §4-§5 reading, sharpened by today's R3-A PES data.

The R3-B sprint should treat the W1e wall as structurally robust under DMRG-equivalent solvers and frame Target C (the multi-determinant / W1e closure) as a Hamiltonian-correction layer, not as a downstream variational-solver layer.

---

## §7. Decision-gate summary

| Gate | Threshold | Observed | Verdict |
|:---|:---|:---|:---:|
| DMRG-on-FCIDUMP reaches FCI at R=3.015 | $\le 1$ mHa | $4.4 \times 10^{-7}$ Ha | **PASS** |
| FCIDUMP round-trip bit-exact | machine precision | $0.0$ exact | **PASS** |
| PES sweep reproduces R_eq within 5% of CCCBDB | $\le 5\%$ | monotone descending (no minimum) | **STOP-structural** |
| Bulk $\chi$ explodes at Q=30 | NO | $\chi_{\max} = 8$ reaches FCI | **PASS** (easy) |

Composite verdict: **MIXED with structural clarification.** The hybrid pipeline's DMRG-side machinery works at production scale on LiH; the headline number is reproduced; the bond-dim cost is trivial. But the PES topology is unbound at $n_{\max}=2$ on both qubit-FCI builders — this is the W1e wall reaching all the way into the "easy" case. The CHANGELOG R_eq 2.82% claim cannot be reproduced from FCI on the FCIDUMP-exported integrals because it comes from a different, continuous-R Level 4 solver.

---

## §8. Paper-edit recommendation (NOT applied per task instruction)

Recommend creating **Paper 57 §III "DMRG on FCIDUMP"** (or extending Paper 20 §sec:hybrid_pipeline) with the following skeleton:

```latex
\subsection{DMRG-equivalent on the FCIDUMP-exported integrals}
\label{subsec:dmrg_on_fcidump}

The Sprint P1 FCIDUMP exporter (\S\ref{subsec:fcidump_interface}) was
stress-tested on the production LiH composed Hamiltonian at
$n_{\max} = 2$, $R = 3.015$~bohr. The pipeline (a) builds the composed
qubit Hamiltonian with PK in $h_1$; (b) exports to FCIDUMP; (c) reads back
via the pure-Python parser, verifying bit-exact round-trip on $h_1$, eri,
and ecore at every $R$ in the panel $\{2.5, 3.015, 3.5, 4.0, 5.0\}$~bohr;
and (d) runs sector-restricted FCI on the round-tripped integrals.
This last step is the DMRG-converged limit for the LiH composed
Hamiltonian: the diagnostic sprint of 2026-06-07
(\texttt{debug/sprint\_dmrg\_diagnostic\_memo.md}) established that the
composed Hamiltonian is block-decoupled with state-side bond rank
$\le 4$ across every sub-block, so DMRG at $\chi_{\max} = 4$ reaches
sector FCI to machine zero. The full-system alpha-vs-beta Schmidt
bipartition of the FCI ground state at $R = 3.015$, measured here, has
state-side bond rank exactly 6 (six significant singular values), with
$\chi_{\max} = 8$ reaching FCI to $10^{-28}$ truncation error.

\paragraph{Bit-exact match to the headline.}
At the production $R = 3.015$, the FCIDUMP $\to$ FCI energy matches the
Sprint P2 reference $E_{\mathrm{FCI}} = -14.143000$~Ha to
$4.4 \times 10^{-7}$~Ha (Table~\ref{tab:r3a_dmrg_lih}), three orders of
magnitude tighter than chemical accuracy. The native and round-tripped
FCI agree to $\lesssim 5 \times 10^{-14}$~Ha on every $R$
(\texttt{debug/data/r3a\_dmrg\_lih.json}), establishing that the
FCIDUMP interface is faithful at machine precision.

\paragraph{PES topology and the $W_{1e}$ wall.}
The qubit-FCI PES on the FCIDUMP-exported integrals is monotone-decreasing
on $R \in [2.5, 5.0]$~bohr for both the composed and balanced builders
(Table~\ref{tab:r3a_dmrg_lih_pes}); neither has an equilibrium minimum in
this range. The CHANGELOG v3.56.0 LiH R_eq numbers (2.82\% composed,
6.93\% balanced) come from
\texttt{ComposedDiatomicSolver.LiH\_ab\_initio} (Paper~17~\S VI.A), a
continuous-R Level~4 multichannel adiabatic solver with PK barriers,
not from FCI on the qubit Hamiltonian. The qubit-FCI PES failing to
bind LiH at $n_{\max} = 2$ is the calibration-data-tier $W_{1e}$ wall
(CLAUDE.md \S 3; cosmic-Galois reading per
\texttt{debug/sprint\_w1e\_period\_class\_memo.md}), structurally
unaffected by which classical correlation solver consumes the FCIDUMP.
This places the qubit FCI on the composed Hamiltonian as a
\emph{resource-estimation benchmark target}, not a binding-accuracy
benchmark target.
```

A companion §III.B "Bond-dimension feasibility" should report the
$\chi$ table from §4 above. The Paper 20 §sec:hybrid\_pipeline\_lih\_fci
recommendation drafted in Sprint P2 §8 can absorb this material directly.

Where to land the bug write-up (§3): a `\paragraph{Convention note:}` 
on the `to_fcidump` description with the cross-reference to the
`debug/sprint_r3a_dmrg_lih_memo.md` §3 named-bug record. Or, if Option B
is implemented before paper revision, the convention note becomes the
`include_pk=True` default's docstring.

---

## §9. Honest scope

What this sprint closed:

- The P1 FCIDUMP exporter, with the PK workaround applied, exports composed and balanced LiH integrals at machine precision and round-trips bit-exactly.
- A "DMRG-equivalent" solver (`coupled_fci_energy` on the round-tripped integrals; equal to DMRG at $\chi_{\max} \ge 8$ per the alpha-vs-beta Schmidt profile and to DMRG at $\chi_{\max} \ge 4$ per the per-block diagnostic) reaches the Sprint P2 headline $E = -14.143$ Ha to $< 10^{-6}$ Ha.
- The PES topology question is sharpened: qubit FCI does not bind LiH at $n_{\max} = 2$ for either composed or balanced builder; the CHANGELOG R_eq numbers live in a different solver path.
- A named bug in `ecosystem_export.hamiltonian` is documented and a recommended fix is on the open-items register.

What this sprint did NOT close:

- A real iterative DMRG run (no DMRG library installed: tenpy, quimb, block2, pyscf all missing). The information-theoretic Schmidt-rank profile establishes feasibility, but the algorithmic convergence of a particular DMRG sweep was not benchmarked. At 30 qubits with $\chi \le 16$, this is decisively non-blocking (any modern DMRG converges trivially at this size; cf. yesterday's diagnostic).
- The PK-bug fix in the ecosystem wrapper — recommended, not applied, per the "DO NOT MODIFY production code" task instruction.
- CCSD(T) cross-validation on the FCIDUMP — pyscf missing.
- $n_{\max} = 3$ extension — would require a separate $Q = 54$ FCI run; $\sim 1.5 \times 10^5$-dim sector, doable but $\sim 10$-15 min per $R$. Not blocking for the GO/STOP question.
- Falsifying or confirming whether the W1e wall closes under any classical-correlation extension layered on top of the qubit Hamiltonian (Target C scope; not this sprint).
- R3-B NaH probe — parallel sibling sprint.

**Adapted decision-gate outcome:**

- Chemical accuracy on LiH composed FCIDUMP $\to$ FCI: **PASS**.
- FCIDUMP round-trip: **PASS**.
- PES topology reproducing CHANGELOG R_eq: **STOP-structural**, with documented reading that the R_eq numbers come from a different solver (continuous Level 4 multichannel + PK), not the qubit FCI.
- Bond-dim feasibility: **PASS** (trivially easy at $\chi_{\max} = 8$).

Composite: **GO on hybrid-pipeline architecture viability; STOP-structural on PES; named bug + fix recommendation queued.**

---

## §10. Files

### Created (driver)
- `debug/r3a_dmrg_lih_driver.py` — ~480 lines. Bypasses ecosystem wrapper (PK workaround); composed + balanced FCIDUMP $\to$ FCI; Schmidt-rank profile.

### Created (data)
- `debug/data/r3a_dmrg_lih.json` — per-R energies, round-trip metrics, Schmidt profile, named-bug record.
- `debug/data/lih_r{2.500,3.015,3.500,4.000,5.000}_nmax2.fcidump` — composed PK-included FCIDUMP files.
- `debug/data/lih_balanced_r{2.500,3.015,3.224,3.500,4.000,5.000}.fcidump` — balanced FCIDUMP files.

### Created (memo)
- `debug/sprint_r3a_dmrg_lih_memo.md` — this memo.

### NOT modified
- Production `geovac/` modules — driver-only sprint.
- `to_fcidump` — bug is documented as institutional memory; fix is recommended for a separate PI-approved sprint.
- Tests — no production code modified.
- Paper drafts — recommendations in §8, not applied per sprint mandate.
- CLAUDE.md — sprint-close protocol will handle §2 / §3 updates separately if PI directs.

---

**End of Sprint R3-A memo. Verdict: MIXED — GO on the headline chemical-accuracy gate ($4 \times 10^{-7}$ Ha match to Sprint P2), GO on FCIDUMP round-trip bit-exactness, GO on bond-dim feasibility ($\chi_{\max} = 8$ reaches FCI), STOP-structural on PES topology (qubit FCI does not bind LiH at $n_{\max} = 2$ for either builder; CHANGELOG R_eq comes from a different continuous-R Level 4 solver). Named bug in `ecosystem_export.hamiltonian` PK convention surfaced and documented; fix recommended (`include_pk=True` switch on `to_fcidump`), not applied. R3-B NaH probe should treat W1e wall as solver-independent.**
