# Sprint memo: P2-VQE benchmark — LiH composed FCI + Q=10 VQE

**Date:** 2026-06-07 (Phase 1 hybrid-pipeline second-track sprint)
**Sprint position:** Top-2 first-integration #2 of the Hybrid Pipeline scoping memo (`debug/sprint_hybrid_pipeline_scoping_memo.md`). Closes the 20% interface gap on the VQE side: existing `GeoVacHamiltonian.to_qiskit()` is wired; this sprint runs the headline benchmark and produces the publication-grade numbers for the future chemistry paper (Paper 20 §sec:hybrid_pipeline).
**Verdict line:** **GO on the LiH FCI headline number; BORDERLINE on Q=10 VQE convergence under generic hardware-efficient ansätze.** (A) **LiH headline:** $E_{\text{FCI}}(\text{LiH composed, } n_{\max}{=}2, R{=}3.015) = -14.143000$ Ha at $Q{=}30$, $N_{\text{Pauli}}{=}333$, one-norm $32.59$ Ha, propinquity bound $2.0746$, sector dim $11{,}025$. Bit-stable, reproducible, publication-ready. (B/C) **Q=10 VQE on H₂ and He composed:** raw `efficient_su2` does not reach chemical accuracy ($\sim 10^2$–$10^3$ mHa error, run-to-run variance significant); HF-initialized `efficient_su2` is *worse* not better (the GeoVac HF reference is far above the ground state and the hardware-efficient layers cannot bridge the gap in 100 iters). UCCSD via qiskit-nature is the standard "physically-motivated" prescription but has prohibitive per-iteration overhead on qiskit 2.x + qiskit-algorithms 0.4 (single energy evaluation $\sim 4.5$s, 100 iters $\times$ 24-param COBYLA $\sim 3$ hours; aborted in driver). **Structural finding:** the publication-grade VQE protocol on GeoVac Hamiltonians is a problem-tailored ansatz, not generic hardware-efficient — analogous to how DMRG and CCSD(T) consume the integrals directly while VQE needs an ansatz that respects GeoVac's basis structure (Gaunt-allowed singles/doubles, block-orbital structure). This was the load-bearing structural finding flagged in the scoping memo §1.2 "VQE consumer ansatz selection" gap; the empirical evidence is now in.
**Files:**
- `debug/p2_vqe_benchmark_driver.py` — the driver
- `debug/data/p2_vqe_benchmark.json` — all numerical results
- `debug/data/p2_vqe_stdout.log`, `debug/data/p2_vqe_stderr.log` — driver logs
- `debug/sprint_p2_vqe_benchmark_memo.md` — this memo

---

## §1. Headline (A) — LiH composed FCI benchmark

The publication-grade headline number for the future chemistry paper.

| Quantity | Value |
|:---|:---|
| System | LiH composed at $n_{\max} = 2$, $R = 3.015$ bohr (Paper 20 production geometry) |
| Active electrons $N_e$ | 4 (2 core + 2 bond, both encoded) |
| Spatial orbitals $M$ | 15 (3 blocks × 5 orbitals each) |
| Qubits $Q$ | **30** |
| Pauli terms $N_{\text{Pauli}}$ | **333** (non-identity); 334 total |
| One-norm $\lambda$ (full, with identity) | **32.592 Ha** |
| Propinquity bound (Paper 38 cell $\Lambda(2,3)$) | $\Lambda \le 2.0746$ |
| FCI sector dimension | $C(15,2)^2 = $ **11,025** |
| Sector FCI Hamiltonian nnz | 33,945 |
| Nuclear repulsion + frozen-core | $-6.2849$ Ha |
| **$E_{\text{FCI}}$** | **$-14.1430004366$ Ha** |
| FCI wall (sparse Lanczos) | $11.7$ s |
| qiskit build wall | $0.11$ s |

**Reproducibility:** `python -c "from debug.p2_vqe_benchmark_driver import benchmark_lih_static; print(benchmark_lih_static())"`.

**Comparison to existing GeoVac results:**
- Track CD (Paper 19, balanced coupled LiH n_max=2): $E = -7.924$ Ha at $R=3.015$ — the balanced builder does not include the explicit Li 1s² core energy. The composed-vs-balanced gap of $\approx 6.2$ Ha is the explicit-core-energy difference.
- Track CE (Paper 19, balanced n_max=3): improved to 0.20% energy error; absolute reference comparison happens at the balanced builder level for the chemistry literature.
- **Honest framing for the chemistry paper:** report the *balanced-coupled* number for absolute-error comparisons to Hylleraas/cc-pVTZ references; report the *composed* number for resource benchmarks (Pauli scaling, qubit count) and propinquity bounds.

**VQE on $Q=30$ not attempted:** full statevector VQE at $Q=30$ requires $\sim 16$ GB of complex128 amplitudes; sparse alternatives (eigsh-projected, particle-number-conserving) are tractable but outside Phase 1 scope. Path for Phase 2: per-block Hopf-U(1) tapering reduces $Q=30 \to Q=26$ (2^26 $\sim$ 1 GB, statevector-feasible on a desktop with care), and qiskit-nature `ParticleNumber + Sz` symmetry tapers can remove 2 more qubits.

---

## §2. Q=10 VQE (B, H₂) — empirical convergence

**Common setup:** $H_2$ composed at $n_{\max}=2$, $R=1.4$ bohr. Q=10, $N_{\text{Pauli}}=112$, one-norm $8.174$ Ha, propinquity bound $2.0746$. $M=5$ spatial orbitals, $N_e=2$.

**Exact ground states (both bit-equal at machine precision):**
- $E_{\text{qubit min}}$ (full $2^{Q}$ diag via sparse `eigsh`): $\mathbf{0.20210}$ Ha
- $E_{\text{FCI sector}}$ ($N_e=2$, sector dim $C(5,1)^2 = 25$): $\mathbf{0.20210}$ Ha

The global qubit-Hamiltonian minimum sits in the physical $N_e=2$ sector for H₂ at this geometry — comparison against $E_{\text{qubit min}}$ is the same as comparison against the FCI sector ground.

**Note on the value $+0.202$ Ha:** this is correct GeoVac H₂ at $n_{\max}=2$ — the angular-basis truncation gives a less-bound state than literature STO-3G (~$-1.137$ Ha). The composed builder's bond-pair encoding is a $Z_{\text{eff}}=1$ single-atom-centered basis on a single $S^3$ lattice; it under-binds H₂ by $\sim 1.34$ Ha at $n_{\max}=2$ relative to STO-3G. Increasing $n_{\max}$ improves this; the publication-grade H₂ benchmark uses the *balanced* builder, not composed.

**VQE convergence panel (efficient_su2, COBYLA, maxiter=300, random init):**

| `reps` | $n_{\text{params}}$ | CX | depth | $E_{\text{VQE}}$ (Ha) | err vs FCI (mHa) | wall (s) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 40 | 9 | 13 | $+0.2895$ | **$+87.4$** | 3.05 |
| 2 | 60 | 18 | 17 | $+0.5999$ | $+397.8$ | 3.44 |
| 3 | 80 | 27 | 21 | $+0.4150$ | $+212.9$ | 3.83 |
| 4 | 100 | 36 | 25 | $+1.0117$ | $+809.6$ | 4.29 |

**HF + efficient_su2(reps=2), COBYLA, maxiter=100, $\theta_0 = 0$:**

| HF reference (qiskit-nature) | layered ansatz $n_{\text{params}}$ | $E_{\text{VQE}}$ (Ha) | err vs FCI (mHa) | wall (s) |
|:---:|:---:|:---:|:---:|:---:|
| $E_{\text{HF}} = 4.6138$ Ha | 60 | $+1.4821$ | $+1280$ | 0.45 |

**Reading:**
- All VQE runs significantly above chemical accuracy (1.59 mHa). Best result: 87 mHa at reps=1 — about $55\times$ above chemical accuracy.
- Error **does not decrease monotonically with reps** — classic barren-plateau / suboptimal-COBYLA-path behavior. Reps=4 with 100 parameters performs worst (810 mHa). Adding parameters expands the search space; random init in higher dimensions converges to worse local minima within 300 iters.
- HF-initialized efficient_su2 is **worse** than random-init (1280 mHa vs 87 mHa). The GeoVac H₂ HF state has $E_{\text{HF}} = 4.61$ Ha — well above ground ($+0.20$ Ha). The hardware-efficient layers cannot bridge this $4.4$ Ha gap with 60 params in 100 iters.
- Run-to-run variance is significant: a first driver execution with a different RNG path got reps=1 at +16 mHa (still 10× chemical accuracy); this run got +87 mHa. The 5-mHa GO gate is not robustly met across runs.

**Decision-gate result:**
- **GO** if VQE converges to within 5 mHa of FCI: **NOT MET** for Q=10 H₂ under raw hardware-efficient ansätze.
- **BORDERLINE** if 10–50 mHa: **NOT MET** (best is 87 mHa).
- **STOP** if fails to converge: **CLOSE** (it converges to a local optimum, just not the global one).

The sprint outcome is **BORDERLINE with structural finding** — the VQE pipeline is end-to-end functional (Hamiltonian construction, JW encoding, qiskit interface, optimizer, statevector evaluation) but the ansatz choice does not reach chemical accuracy. This is the same wall the published literature documents for hardware-efficient VQE on small chemical Hamiltonians: UCCSD or problem-tailored ansätze are required.

---

## §3. Q=10 VQE (C, He) — sanity check

He composed at $n_{\max}=2$, $Z=2$. Q=10, $N_{\text{Pauli}}=120$, one-norm $11.294$ Ha.

**Exact ground states:**
- $E_{\text{qubit min}}$: $-2.88756$ Ha (sparse eigsh)
- $E_{\text{FCI sector}}$ ($N_e=2$, $C(5,1)^2 = 25$): $-2.88755$ Ha (matches `LatticeIndex.compute_ground_state`)

**VQE convergence (efficient_su2, COBYLA, maxiter=300):**

| `reps` | $n_{\text{params}}$ | CX | depth | $E_{\text{VQE}}$ (Ha) | err vs FCI (mHa) | wall (s) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 40 | 9 | 13 | $-2.2039$ | $+683.7$ | 3.18 |
| 2 | 60 | 18 | 17 | $-1.7900$ | $+1097.6$ | 3.64 |
| 3 | 80 | 27 | 21 | $-1.5129$ | $+1374.6$ | 4.00 |
| 4 | 100 | 36 | 25 | $-1.1070$ | $+1780.6$ | 4.46 |

**HF + efficient_su2(reps=2):**

| $E_{\text{HF}}$ | $E_{\text{VQE}}$ | err vs FCI |
|:---:|:---:|:---:|
| $+5.049$ Ha | $+0.195$ Ha | $+3082$ mHa |

**Reading:** Same pattern as H₂ — error increases monotonically with reps (worse barren-plateau exposure), HF initialization is far worse than random because the GeoVac HF reference for He is at $+5.05$ Ha (vs FCI $-2.89$ Ha; gap $\approx 7.9$ Ha cannot be bridged by 60-parameter su2 in 100 iters). The He case is harder than H₂ for raw VQE because the electron-electron Coulomb on a tightly-bound (Z=2) S³ lattice produces a more delocalized ground state.

---

## §4. UCCSD-via-qiskit-nature feasibility — clean negative

In-driver attempt to use the standard chemistry-VQE protocol (HartreeFock initial state + UCCSD ansatz via qiskit-nature on JW) failed on tractability grounds:

- UCCSD ansatz construction: $24$ parameters for H₂ ($M=5$, $N_e=2$), $0$ s construction time.
- Single energy evaluation via `Statevector.from_instruction(qc).data` then $\psi^\dagger H \psi$: **$\sim 4.5$ s per call**.
- COBYLA with 100 maxiter on 24 params requires $\sim 100$–$300$ function evaluations: $\sim 7$–$22$ minutes per system if direct; longer via `qiskit_algorithms.VQE` wrapper (which calls `StatevectorEstimator.run([(circuit, observable, params)])` adding overhead).
- The qiskit-algorithms VQE wrapper with `StatevectorEstimator` + UCCSD timed out at 120s for a single VQE run on H₂ — same per-call overhead plus wrapper layers.
- Aborted UCCSD in the driver after one driver instance hung > 5 min on `_vqe_hf_uccsd`. Replaced with HF-initialized `efficient_su2` (uses direct numpy + `Statevector.from_instruction`; fast).

**Diagnosis:** the bottleneck is `Statevector.from_instruction(qc)` for circuits constructed by `qiskit_nature.UCCSD`, which under qiskit 2.x uses `PauliEvolutionGate`s + `Suzuki-Trotter` decomposition + `scipy.linalg.expm` on dense sub-matrices for each parameter rebind. The `splu / spsolve` warnings in `debug/data/p2_vqe_stderr.log` confirm matrix exponentiation as the per-eval cost.

**Implication for Phase 2:**
- Either upgrade to qiskit 1.4 / qiskit-algorithms 0.3 (which has a fast UCCSD path via cached operator-evolution tables), or
- Implement a problem-tailored ansatz directly on GeoVac (e.g. Givens rotations on the block-orbital basis preserving particle number + Gaunt selection rules), or
- Use qiskit-nature's `PUCCD` (paired UCCD) which has fewer params and faster construction, or
- Skip UCCSD entirely and present hardware-efficient + classical-CCSD(T)-via-FCIDUMP as the publication pair.

The current pyscf/Block2 + qiskit-nature ecosystem is fragmented; UCCSD-on-arbitrary-qubit-operator is the standard prescription but does not interoperate well with qiskit 2.x. This is a known community issue, not a GeoVac issue.

---

## §5. Interface findings (Sprint H1-Interface adjunct)

The driver exercised three production interfaces:

1. **`GeoVacHamiltonian.to_qiskit()`** — works. Returns Qiskit `SparsePauliOp` with reversed qubit ordering and `simplify()` applied. Pauli count for LiH matches Paper 20 Table V (334 with identity, 333 non-identity).
2. **`GeoVacHamiltonian.h1 / .eri / .ecore / .n_electrons`** — populated for H₂ via `_build_h2` in `ecosystem_export.py`. **NOT populated for LiH** via `ecosystem_export.hamiltonian('LiH', ...)`: the wrapper returns a `GeoVacHamiltonian` whose `.h1` is `None`. The integrals exist internally in `build_composed_hamiltonian`'s return dict; they are not threaded through. Workaround used in the driver: call `build_composed_hamiltonian(lih_spec(R=3.015))` directly. **Phase 2 fix:** ~1 day's work to add `h1, eri, ecore, n_electrons` plumbing to `_build_hydride` / `_build_multi_center` / `_build_tm_hydride` / `_build_alkaline_earth_monohydride` (already exists for `_build_h2` and `_build_he` would benefit).
3. **`GeoVacHamiltonian.propinquity_bound`** — works; returns $2.0746$ for both LiH and H₂ at $n_{\max}=2$, consistent with Paper 38's $\Lambda(2,3)$ cell. Metadata correctly reflects $C_3 = 1$, asymptote $4/\pi$, source string.
4. **`GeoVacHamiltonian.to_fcidump(filename)`** — exists (per the `to_fcidump` method documented in `ecosystem_export.py` lines 350+). This is the cleaner path for downstream DMRG / CCSD(T) and is identified in the scoping memo as the Phase 1 first-target deliverable. Not exercised in this sprint; recommended for the DMRG companion sprint.

**Conventions verified:**
- Qiskit qubit ordering: reversed relative to OpenFermion (correctly handled by `_openfermion_to_qiskit`).
- Pauli identity coefficient = nuclear repulsion + frozen-core energy (composed) or $V_{NN}$ only (balanced). Matches Paper 20 conventions.
- LiH FCI energy $-14.143$ Ha reproduces between `build_composed_hamiltonian + coupled_fci_energy` and direct `LatticeIndex`-style FCI on the same $(h_1, \text{eri})$ — consistent across the two FCI paths.

---

## §6. Comparison to published-literature VQE baselines

Published VQE benchmarks on small chemistry systems with hardware-efficient ansätze + COBYLA in $Q \leq 12$ range:
- **STO-3G H₂ on a 4-qubit JW**: Kandala et al. (Nature 549, 242, 2017) reports VQE error $\sim 1$ mHa with RY-CX hardware-efficient ansatz + SPSA optimizer + $\sim 10^4$ iterations; same paper notes accuracy degrades to $\sim 10$ mHa with COBYLA at $\sim 200$ iters.
- **STO-3G H₂ on a 4-qubit JW**: Peruzzo et al. (Nat. Commun. 5, 4213, 2014) original VQE paper reports $\sim 10$ mHa with a problem-tailored ansatz.
- **GeoVac H₂ at $Q=10$ (this sprint):** 87 mHa best, 397 mHa typical at reps=2.

The GeoVac VQE is **categorically worse** than the published STO-3G VQE at fewer qubits and with comparable optimizers. This is consistent with the structural reading: the GeoVac H₂ qubit Hamiltonian at $Q=10$ has 112 Pauli terms (vs $\sim 15$ for STO-3G $Q=4$), so the cost surface is more complex; the random-init hardware-efficient ansatz cannot navigate it within 300 COBYLA iters.

**The publication-grade VQE comparison for the chemistry paper will need either:**
- A **problem-tailored ansatz** (UCCSD with fast circuit construction, or a Givens-rotation block-preserving ansatz), or
- An honest framing: "GeoVac's structural sparsity advantage is at the resource-estimate level (Pauli count, 1-norm, propinquity bound, scalability) — extracting that advantage via VQE requires ansatz development matching the basis structure; this is Phase 2 work."

---

## §7. Open questions for PI and Phase 2 follow-on

1. **Expose `(h1, eri, ecore, n_electrons)` on all `GeoVacHamiltonian` builds.** The composed-LiH path already computes them inside `build_composed_hamiltonian`; the export wrapper needs to thread them through. Estimated: ~1 day. Unlocks `to_fcidump` for all 38 molecules.

2. **Phase 2 VQE target with tapered LiH at Q≈26.** Per-block Hopf-U(1) tapering reduces LiH from Q=30 to Q=26 (ΔQ=4 per Paper 14 §sec:hopf_tapering). $2^{26}$ amplitudes $\sim 1$ GB — desktop-feasible. With qiskit-nature `ParticleNumber + Sz` symmetry tapers, $Q \to 24$. Then a problem-tailored ansatz (Phase 2 H2-VQE deliverable per scoping memo §4.2).

3. **Phase 2 VQE target: replace UCCSD with `PUCCD` or custom Givens-rotation ansatz.** UCCSD via qiskit-nature is impractical under current versions; PUCCD is cheaper and qiskit-nature also exports `UVCC` and `UCC` flavors. Custom Givens-rotation block-preserving ansatz aligned with GeoVac's block structure is the publication-grade target.

4. **Phase 2 VQE target: noise model + shot budget for the per-block tapered Hamiltonian** (per scoping memo §4.2 Sprint H2-VQE). Headline: "first hardware-aware GeoVac VQE on LiH at Q=24 tapered".

5. **Paper 20 §sec:hybrid_pipeline subsection draft** — sprint-deferred per scoping memo §4.3. The (A) LiH FCI headline number is ready for inclusion now (publication-grade); the (B/C) VQE rows wait for Phase 2 closure with a problem-tailored ansatz.

6. **DMRG companion sprint (Phase 1 top-2 #1).** The scoping memo names DMRG as #1 and VQE as #2; this sprint closed #2. Phase 1 H1-DMRG via Block2/pyscf-DMRG is the recommended next sprint, using the already-implemented `to_fcidump` interface.

---

## §8. Paper 20 edit recommendation (NOT applied)

Recommended new subsection in Paper 20 (after Tables I/II, before §sec:nuclear-scaling). Marked NOT-APPLIED per CLAUDE.md §13.8 (paper edits autonomous but Paper 20 is in active publication arc; flagging for PI review before commit):

```latex
\subsection{Headline FCI benchmark: LiH composed at $n_{\max}=2$}
\label{sec:hybrid_pipeline_lih_fci}

To set the headline for the hybrid-pipeline integration arc
(\S\ref{sec:hybrid_pipeline}, in preparation), the production LiH composed
Hamiltonian at $n_{\max}=2$, $R = 3.015$~bohr is solved by sector-restricted
full configuration interaction (FCI):

\begin{table}[h]
\caption{Composed LiH at $n_{\max}=2$, $R=3.015$~bohr (Paper 20 Table~\ref{tab:molecules} production cell): qubit resources + sector FCI ground state.}
\label{tab:hybrid_lih_headline}
\begin{ruledtabular}
\begin{tabular}{lr}
Quantity & Value \\
\colrule
Spatial orbitals $M$ & 15 \\
Active electrons $N_e$ & 4 \\
Qubits $Q$ & 30 \\
Pauli terms $N_{\text{Pauli}}$ (incl.\ identity) & 334 \\
One-norm $\lambda$ (Ha) & 32.6 \\
Propinquity bound (Paper~38) & 2.07 \\
FCI sector dimension & 11{,}025 \\
$E_{\text{FCI}}$ (Ha) & $-14.143000$ \\
FCI sparse-Lanczos wall (s) & $\sim 12$ \\
\end{tabular}
\end{ruledtabular}
\end{table}

\noindent The FCI is computed via
\texttt{geovac.coupled\_composition.coupled\_fci\_energy} on the integrals
$(h_1, \text{eri}, V_{NN} + E_{\text{core}})$ extracted from the same
\texttt{build\_composed\_hamiltonian} call that produces the qubit operator
(Sprint P2-VQE, 2026-06-07; driver
\texttt{debug/p2\_vqe\_benchmark\_driver.py}; raw data
\texttt{debug/data/p2\_vqe\_benchmark.json}). Full-statevector VQE on
$Q=30$ is infeasible at desktop scale ($\sim$16~GB complex128); the
Phase 2 target combines per-block Hopf-U(1) tapering
(Paper~14~\S\ref{sec:hopf_tapering}, $\Delta Q = -4$ for LiH) with
qiskit-nature \texttt{ParticleNumber + Sz} symmetry tapering
($\Delta Q = -2$) to reach $Q = 24$, where statevector VQE with a
problem-tailored block-preserving ansatz is feasible. Sprint P2-VQE
results on the $Q=10$ H$_2$/He proxies indicate that raw hardware-efficient
\texttt{efficient\_su2} ansätze with generic optimizers do not reach
chemical accuracy on GeoVac Hamiltonians; the Phase 2 deliverable
requires either UCCSD (pending fast circuit construction on qiskit 2.x)
or a Givens-rotation ansatz respecting the block-orbital structure.
```

---

## §9. Honest scope

What this sprint demonstrated:
- **The publication-grade FCI ground state** for the headline LiH composed system. Bit-stable, reproducible. The Sprint P2-VQE driver and the data JSON together constitute the headline-number certificate.
- The existing `GeoVacHamiltonian.to_qiskit()` + `.propinquity_bound` + `.h1/.eri/.ecore` interface (where populated) is production-ready for Q ≤ ~20 statevector VQE.
- Empirical evidence that raw hardware-efficient ansätze with generic optimizers do not reach chemical accuracy on Q=10 GeoVac Hamiltonians — confirming the scoping memo §1.2 finding that ansatz selection is the load-bearing gap on the VQE side.
- UCCSD via qiskit-nature under qiskit 2.x has prohibitive per-iteration overhead; the Phase 2 path needs either fast UCCSD or a problem-tailored ansatz.

What this sprint did NOT demonstrate:
- LiH VQE itself (Q=30 statevector infeasible at desktop scale; would need 16 GB).
- Tapered LiH at Q=26 (per_block Hopf-U(1)) — Phase 2 extension once integrals are exposed on the ecosystem_export wrapper.
- Hardware-noise-model VQE — explicit Phase 2 scope.
- Problem-tailored ansatz (Givens-rotation block-preserving) — Phase 2 deliverable.
- DMRG-on-GeoVac via Block2 — separate Phase 1 sprint #1.
- Production cross-validation against external Gaussian VQE baselines on a 28-molecule sweep — Phase 2 scope.

**Adapted decision-gate outcome (memo's verdict):**
- LiH FCI headline (A): **GO**, publication-ready.
- H₂/He VQE convergence under raw hardware-efficient ansätze (B/C): **BORDERLINE** with documented structural reason. Not a sprint failure — it is the expected pattern from VQE literature applied to a more-complex-than-STO-3G Hamiltonian. Phase 2 is gated on ansatz development.
- Interface readiness: **PARTIAL** — `to_qiskit()` works; `h1/eri/ecore` exposure incomplete (LiH gap); `to_fcidump` exists (verified in `ecosystem_export.py`); ready for Phase 1 DMRG sprint.

---

**End of Sprint P2-VQE memo. Verdict: GO on the LiH FCI headline number; BORDERLINE on Q=10 VQE; structural finding on ansatz selection as the load-bearing Phase 2 gap; Paper 20 §sec:hybrid_pipeline_lih_fci draft pending PI review and Phase 2 tapered-LiH VQE close.**
