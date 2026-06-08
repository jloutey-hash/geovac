# Sprint R3-B — NaH DMRG-on-FCIDUMP falsifier (Phase 1 hybrid pipeline)

> **2026-06-07 addendum (Sprint W1e-Projection-Audit, same day):** the
> qualitative verdict ("balanced over-attracts NaH; DMRG = FCI at finite
> sector dim does not close W1e") SURVIVES. The PES *magnitudes* reported in
> this memo were affected by the same V_NN call-convention bug that affected
> R3-A: pre-fix Path B was masking ~1.3 Ha of true NaH over-attraction depth
> at R=2.5 bohr. With the bug fixed (`MolecularSpec.R` field, 2026-06-07
> patch in `geovac/balanced_coupled.py`), the corrected NaH PES at R=2.5 is
> $-171.097$ Ha (was $-169.781$ Ha in this memo's Path B). The well is even
> *more* deeply over-attractive than reported; the qualitative "no interior
> minimum" / "DMRG cannot close it" verdict is structurally unchanged. See
> `debug/sprint_w1e_projection_audit_memo.md` §§6.1 and 6.2 for the full
> Path-A vs Path-B comparison and bug fix.

**Date:** 2026-06-07.
**Sprint position:** Round 3 of the hybrid-pipeline Phase 1. Load-bearing falsifier for the cosmic-Galois Class 1 / Class 3 reading of W1e. Companion to Sprint P1 (FCIDUMP exporter) and Sprint P4 (NaH balanced FCI baseline). No new physics — applies the existing pipeline to the falsifier question.
**Verdict line:** **STOP. DMRG-on-GeoVac-balanced-FCIDUMP DOES NOT close the W1e wall on NaH. The PES is monotone-ascending toward large R at both $n_{\max}=2$ and $n_{\max}=3$, with no interior minimum anywhere in $R \in [2.5, 6.0]$ bohr. At the experimental $R_e = 3.566$ bohr the framework over-binds by $-1.691$ Ha (n_max=2, $23.7\times$ exp $D_e$) and $-1.502$ Ha (n_max=3, $21.1\times$ exp $D_e$) — bit-identical to the P4 baseline.** W1e is NOT pure Class-1 / Class-3 multi-determinant correlation; it sits at the Hamiltonian specification level and is not solver-closable.

**Cross-references:** `debug/sprint_p4_nah_w1e_baseline_memo.md` (anchor numbers), `debug/sprint_p1_fcidump_exporter_memo.md` (FCIDUMP exporter), `debug/sprint_w1e_period_class_memo.md` (structural classification under test), `papers/group3_foundations/paper_18_exchange_constants.tex` §IV.6 (inner-factor input-data tier — chemistry-side analog), CLAUDE.md §1.7 multi-focal-composition wall pattern (W1e is the 6th instance), CLAUDE.md §3 W1e dead-end rows.

---

## §0. Executive summary

### §0.1 The falsifier and the structural answer

The hybrid-pipeline arc rests on the working hypothesis that the W1e chemistry wall is multi-determinant correlation outside GeoVac's outer-factor period engine but accessible to classical correlation solvers (DMRG, CCSD(T), AFQMC). The Class 1 reading splits the wall into (a) calibration/input-data tier (Paper 18 §IV.6) and (b) multi-determinant correlation (Class 3); the GO scenario is that classical correlation solvers close (b), leaving only (a) as residual.

The structural answer is mathematically prior to any DMRG run: **the framework's `coupled_fci_energy` IS the exact FCI of the GeoVac Hamiltonian**. At the production NaH balanced parameters the FCI Hilbert space is small enough ($\dim_{\text{FCI}} = 100$ at $n_{\max}=2$, $784$ at $n_{\max}=3$) that the $\chi \to \infty$ DMRG limit is just sector-restricted FCI. Any DMRG run on the same FCIDUMP at sufficient bond dimension MUST reproduce the framework FCI; it cannot produce a different answer.

This sprint runs the falsifier honestly anyway, with the operational test: (Path A) export FCIDUMP from `build_balanced_hamiltonian`, parse it back, run an independent FCI on the reparsed integrals; (Path B) run the framework's own `coupled_fci_energy` on the in-memory result; and (C) compare both against the P4 baseline. The three numbers must agree, and the resulting PES is the *χ → ∞ DMRG-on-the-balanced-FCIDUMP* — the empirical ceiling of any classical correlation solver against the same Hamiltonian.

### §0.2 Numerical verdict (bit-identical to P4 at both $n_{\max}$)

| $n_{\max}$ | $R$ (bohr) | $E_{\text{R3B direct}}$ (Ha) | $E_{\text{R3B via FCIDUMP}}$ (Ha) | $E_{\text{P4}}$ (Ha) | direct vs FCIDUMP |
|:---:|:---:|---:|---:|---:|---:|
| 2 | 2.500 | $-169.781376$ | $-169.781376$ | $-169.7814$ | $5.7 \times 10^{-13}$ |
| 2 | 3.000 | $-169.391423$ | $-169.391423$ | $-169.3914$ | $2.3 \times 10^{-13}$ |
| 2 | **3.566** | $-169.114635$ | $-169.114635$ | $-169.1146$ | $5.1 \times 10^{-13}$ |
| 2 | 4.000 | $-168.966136$ | $-168.966136$ | $-168.9661$ | $5.1 \times 10^{-13}$ |
| 2 | 4.500 | $-168.808556$ | $-168.808556$ | $-168.8086$ | $0$ |
| 2 | 5.000 | $-168.636527$ | $-168.636527$ | $-168.6365$ | $1.7 \times 10^{-13}$ |
| 2 | 6.000 | $-168.235596$ | $-168.235596$ | $-168.2356$ | $2.8 \times 10^{-13}$ |
| 3 | 2.500 | $-170.052963$ | $-170.052963$ | $-170.0530$ | $2.0 \times 10^{-13}$ |
| 3 | 3.000 | $-169.754125$ | $-169.754125$ | $-169.7541$ | $5.1 \times 10^{-13}$ |
| 3 | **3.566** | $-169.568414$ | $-169.568414$ | $-169.5684$ | $6.3 \times 10^{-13}$ |
| 3 | 4.000 | $-169.460487$ | $-169.460487$ | $-169.4605$ | $3.4 \times 10^{-13}$ |
| 3 | 4.500 | $-169.319704$ | $-169.319704$ | $-169.3197$ | $5.7 \times 10^{-14}$ |
| 3 | 5.000 | $-169.147213$ | $-169.147213$ | $-169.1472$ | $6.3 \times 10^{-13}$ |
| 3 | 6.000 | $-168.747304$ | $-168.747304$ | $-168.7473$ | $1.4 \times 10^{-13}$ |

The two R3-B paths (direct FCI vs FCIDUMP round-trip + independent FCI) agree to $\sim 10^{-13}$ Ha at every grid point — the FCIDUMP exporter preserves the Hamiltonian content losslessly through the bit-exact text round-trip. R3-B agrees with the P4 baseline to the rounding precision of the published P4 table ($\sim 10^{-5}$ Ha = the 4-digit display in `sprint_p4_nah_w1e_baseline_memo.md`).

FCIDUMP round-trip diagnostics: $\max |h_1 - \text{parsed}(h_1)| = 0$ at every point; $\max |\text{eri} - \text{parsed}(\text{eri})| \in \{0, 1.7 \times 10^{-16}\}$; $|\text{ecore} - \text{parsed}(\text{ecore})| = 0$. The single ULP discrepancy on $\text{eri}$ at $n_{\max}=3$ is at machine epsilon — not a real difference.

### §0.3 The load-bearing PES shape

Falsifier gate: GO iff interior minimum at $R_{eq} \in [3.0, 4.5]$ bohr with $D_e \in [0.04, 0.15]$ Ha (within factor 2 of experimental $D_e = 0.0713$ Ha).

**Outcome at $n_{\max} = 2$:** No interior minimum. PES monotone ascending across $R \in [2.5, 6.0]$ bohr. Slope $\approx +0.30$ Ha/bohr in the $R \in [3.566, 5.0]$ region (ascending toward dissociation, i.e. the absolute minimum at $R = 2.5$ bohr or smaller — same artifact-well as P4).

**Outcome at $n_{\max} = 3$:** No interior minimum. PES monotone ascending. Slope $\approx +0.26$ Ha/bohr in the $R \in [3.566, 5.0]$ region. Same artifact-well shape; $n_{\max} = 2 \to 3$ DEEPENS the well at every $R$ by $\sim 0.2-0.5$ Ha (more basis flexibility ⇒ more overattraction), exactly matching the P4 §2.1 observation.

**Wall quantification at the experimental $R_e^{\text{exp}} = 3.566$ bohr** (anchored against the P4 $R = 8$ bohr dissociation reference, which R3-B's grid does not extend to but which P4 establishes as $-167.4236$ Ha at $n_{\max}=2$ and $-168.0667$ Ha at $n_{\max}=3$):

| $n_{\max}$ | $E(R_e^{\text{exp}}) - E(R=8)$ | as multiple of exp $D_e = 0.0713$ Ha |
|:---:|---:|:---:|
| 2 | $-1.691$ Ha (over-binds) | $23.7\times$ |
| 3 | $-1.502$ Ha (over-binds) | $21.1\times$ |

These are bit-identical to P4 §0.2 — the W1e wall on the balanced NaH FCIDUMP is unchanged by the DMRG-at-infinite-χ route. The Path-A vs Path-B internal consistency check ($\sim 10^{-13}$ Ha) confirms this is not a pipeline plumbing artifact.

### §0.4 The load-bearing structural verdict

**The W1e wall lives at the Hamiltonian-specification level, NOT at the solver level.** Any classical correlation engine (DMRG, CCSD(T), AFQMC, MRCI) consuming the same balanced FCIDUMP can at best match the framework FCI; it cannot fix an overattraction artifact baked into the (h1, eri, ecore) integrals before any solver is run.

The R3-B falsifier is decisive: the cosmic-Galois reading that W1e is pure Class-1 + Class-3 multi-determinant correlation accessible to classical correlation solvers is **falsified for NaH at balanced architecture**. W1e survives the DMRG limit because the wall's content is in the Hamiltonian, not in the solver's correlation depth.

---

## §1. Method

### §1.1 The two-path test

The falsifier compares **two independent FCI calculations on the same Hamiltonian** at each $R$:

* **Path A — FCIDUMP loop.** `build_balanced_hamiltonian(spec, R, L_max=4)` → wrap into `GeoVacHamiltonian` → `to_fcidump(file)` (Sprint P1 exporter, Knowles-Handy format) → `read_fcidump(file)` (Sprint P1 pure-Python reader) → independent `coupled_fci_energy` on the parsed $(h_1, \text{eri}, \text{ecore})$. This is the "DMRG-at-infinite-$\chi$" benchmark: it exercises everything that any classical correlation solver would touch — the standardised text exchange format and a particle-number-projected FCI on its output.

* **Path B — direct FCI.** `coupled_fci_energy` called on the in-memory `build_balanced_hamiltonian` result. This is the P4 baseline path.

Path A = Path B mathematically (the FCIDUMP round-trip is lossless and the FCI is the same projection on the same Hilbert sector). Observing equality empirically confirms (i) the FCIDUMP exporter is intact, (ii) the read-back tensors are bit-identical, (iii) any future DMRG/CCSD(T) on the same FCIDUMP is bound from above by this energy. Observing equality between Path B and the P4 published numbers confirms reproducibility of the falsifier baseline.

### §1.2 Why DMRG-at-infinite-χ is the correct benchmark for the load-bearing question

The hybrid-pipeline scoping memo named DMRG (via block2 / pyscf-DMRG) as the natural Phase-1 solver. Three observations make exact FCI the correct gold-standard at these sizes:

1. **FCI dimension is tiny.** At $n_{\max} = 2$, $M = 10$ spatial orbitals, $N_{\text{el}} = 2$, the sector-restricted FCI is $\binom{10}{1}^2 = 100$ determinants. At $n_{\max} = 3$, $M = 28$, $N_{\text{el}} = 2$, it is $\binom{28}{1}^2 = 784$ determinants. Both are tractable in $< 200$ ms via sparse diagonalisation.

2. **Any DMRG result is bounded above by FCI.** A DMRG calculation at bond dimension $\chi$ is variational over the rank-$\chi$ MPS manifold; at $\chi \to \infty$ it equals the FCI (in the Wilks–Vidal sense — every state on a finite spin chain has an MPS representation of finite bond dimension). For closed-shell singlets on $\leq 56$ qubits the required $\chi$ to reach FCI is bounded by $\chi_{\max} = \sqrt{\dim_{\text{FCI}}}$ which is $10$ at $n_{\max} = 2$ and $\approx 28$ at $n_{\max} = 3$ — well below any production-DMRG bond-dimension regime.

3. **The PES shape is the load-bearing observable.** Whether DMRG closes W1e is a question about the location of the minimum and its depth, not about the absolute energy. The FCIDUMP-at-infinite-χ test gives the *bestcase*: if even the exact FCI of the Hamiltonian does not produce an interior minimum near $R_e^{\text{exp}}$, no finite-$\chi$ DMRG can.

The structural reading P4 §0.3 stated already gives the answer; the run reported below verifies it empirically and confirms pipeline integrity at machine precision.

### §1.3 R-grid

$R \in \{2.5, 3.0, 3.566, 4.0, 4.5, 5.0, 6.0\}$ bohr; the task spec called for $\{2.5, 3.0, 3.566, 4.0, 4.5, 5.0, 6.0\}$, which is what was run. P4 covers $R \in [2.0, 8.0]$ on a 9-point grid; the R3-B grid is the inner 7-point subset, sufficient to characterise the PES shape and the wall, and chosen to keep the seven-point gate identical to the falsifier-task scope. The $R = 8$ bohr P4 dissociation anchor is referenced as fixed external input; the wall-quantification ratios use that value verbatim.

### §1.4 Falsifier gate (verbatim from task)

* **GO**: interior minimum at $R_{eq} \in [3.0, 4.5]$ bohr with $D_e \in [0.04, 0.15]$ Ha (within factor 2 of experimental $0.0713$ Ha).
* **BORDERLINE**: interior minimum exists but $R_{eq}$ off by $>15\%$ or $D_e$ off by $>$ factor 3.
* **STOP**: monotone descent / no interior minimum / overbind by $\geq 10\times$.

---

## §2. Results

### §2.1 PES tables (already given in §0.2 above)

### §2.2 Pipeline-integrity certificates

For each $(n_{\max}, R)$:

| certificate | observed | gate |
|:---|---:|:---:|
| max $|h_1^{\text{export}} - h_1^{\text{parsed}}|$ | $0.0$ everywhere | $< 10^{-12}$ |
| max $|\text{eri}^{\text{export}} - \text{eri}^{\text{parsed}}|$ | $0$ (n_max=2), $1.7 \times 10^{-16}$ (n_max=3, single ULP) | $< 10^{-12}$ |
| $|\text{ecore}^{\text{export}} - \text{ecore}^{\text{parsed}}|$ | $0$ everywhere | $< 10^{-14}$ |
| Path A vs Path B FCI energy | $\leq 6.3 \times 10^{-13}$ Ha (max over 14 points) | $< 10^{-10}$ Ha |
| Path B vs P4 published | $\leq 4.5 \times 10^{-5}$ Ha (rounding precision of P4 table) | matches P4 rounding |

All certificates pass at machine precision. The R3-B Path-A loop is the formal "FCIDUMP-out, FCIDUMP-in, FCI-solver" round-trip and constitutes a non-trivial pipeline-validation result for the Sprint P1 exporter, independent of the falsifier question: it confirms that the exporter is consumable by an independent FCI solver and that the resulting energy matches the framework's own.

### §2.3 Interior-minimum search

At $n_{\max} = 2$: the 7-point PES is monotone ascending — $E$ increases at every consecutive $R_{i} \to R_{i+1}$ step. The minimum on the grid sits at the smallest $R$ ($R = 2.5$ bohr) at $E = -169.781$ Ha, which is the absolute artifact-well of P4 §0.2 (the well at $R = 2.0$ bohr is even deeper at P4's $E = -170.31$ Ha but R3-B does not extend to $R = 2.0$). No turn-over inside the grid.

At $n_{\max} = 3$: same shape. Minimum on grid at $R = 2.5$ bohr at $E = -170.053$ Ha.

The slope analysis:

| $R$ interval | $n_{\max}=2$ slope (Ha/bohr) | $n_{\max}=3$ slope (Ha/bohr) |
|:---|---:|---:|
| 2.5 → 3.0 | $+0.78$ | $+0.60$ |
| 3.0 → 3.566 | $+0.49$ | $+0.33$ |
| 3.566 → 4.0 | $+0.34$ | $+0.25$ |
| 4.0 → 4.5 | $+0.32$ | $+0.28$ |
| 4.5 → 5.0 | $+0.34$ | $+0.34$ |
| 5.0 → 6.0 | $+0.40$ | $+0.40$ |

All slopes are positive everywhere — the PES is ascending toward the dissociation limit at every step. No sign change of $dE/dR$ anywhere in the test region. This is the operational definition of "no interior minimum" and triggers the STOP verdict.

### §2.4 W1e wall direction

P4 §2.2 quantified the W1e wall at the experimental $R_e$ as $+1.691$ Ha over-binding ($n_{\max}=2$) and $+1.502$ Ha ($n_{\max}=3$) relative to the framework's $R = 8$ dissociation anchor. The R3-B numbers are bit-identical (the same FCI on the same Hamiltonian); the closure target of $\sim 1.5$ Ha downward shift identified in P4 §0.3 is therefore equally unachieved at the DMRG-at-infinite-$\chi$ limit. The framework over-binds by $20\text{--}24 \times$ experimental $D_e$ at the same Hamiltonian, regardless of the classical solver applied to it.

---

## §3. Verdict and structural implications

### §3.1 Net verdict

**Net verdict: STOP. The cosmic-Galois Class-1 reading of W1e — that classical correlation solvers consuming the balanced FCIDUMP would close the wall — is empirically falsified for NaH at the production architectures.** The DMRG-at-infinite-$\chi$ limit reproduces the P4 baseline bit-exactly and does not shift the PES shape away from the artifact-well configuration. W1e survives the DMRG limit because the wall is not in the determinant expansion; it is in the Hamiltonian specification.

This sharpens (does not refute) the Sprint W1e period-class reading: W1e is structurally in the **inner-factor input-data tier** (Paper 18 §IV.6 chemistry-side analog). The current sprint adds the operational consequence: the wall is not solver-closable on the balanced FCIDUMP, no matter which classical correlation solver is used.

### §3.2 Reclassification of W1e

The earlier reading (CLAUDE.md §1.7 multi-focal-composition wall pattern, with W1e as the 6th instance) split the wall into three contribution classes:
* **Class 1**: external atomic-physics calibration data (Clementi-Raimondi exponents for the [Ne] frozen core, hydrogenic baseline, FrozenCore profile).
* **Class 2**: multi-focal composition (the structural frame — outer-factor machinery cannot autonomously generate external inputs).
* **Class 3**: multi-determinant correlation (the partial contribution visible only in explicit [Ne] active-space promotion).

R3-B's bit-identical reproduction of the P4 baseline confirms: the full sprint's-worth of "going beyond GeoVac's FCI to a more powerful classical solver" closes Class 3 *completely* (DMRG-at-infinite-χ ≡ FCI on the same Hamiltonian; the FCI dimension at $n_{\max} = 2,3$ is small enough that the classical solver has perfect access to all determinants). Yet **the W1e wall is essentially unchanged**.

Therefore the Class-3 contribution to W1e on the *balanced* FCIDUMP is structurally bounded: anything Class-3 that lives in this 2-electron / $M=10\text{ or }28$ Hilbert sector has already been included in the framework FCI. The remaining $\sim 1.5$ Ha gap is **almost entirely Class 1 + Class 2**, with vanishing Class-3 contribution from this Hamiltonian.

The Class-3 reading is not abandoned globally — it would still apply to (a) basis enlargement beyond max_n=3 (where the F6 sprint already saw 10.2% closure ceiling at max_n=4), or (b) explicit promotion of the [Ne] frozen-core electrons into the active space (10 electron Hilbert sector, $\dim_{\text{FCI}} \sim 10^9$ — outside this sprint's scope and the load-bearing claim). But at the production NaH balanced FCIDUMP that the hybrid pipeline was naturally going to consume, Class 3 is closed and Class 1 dominates the residual.

This is the **chemistry-side analog of the Sprint H1 Yukawa non-selection theorem with the falsifier closed**: not just structurally (we cannot derive the Yukawa from outer-factor machinery), but operationally — the standard external-solver pipeline cannot derive missing chemistry calibration content from the same outer-factor data either.

### §3.3 Implications for Paper 57 abstract

The downstream implications for the Paper 57 abstract:

* **The "W1e closure under DMRG" structure is now ruled out** as a content vector. Paper 57 cannot frame W1e as solver-closable on the balanced FCIDUMP.
* **The natural Paper 57 framing is**: "GeoVac as Hamiltonian producer + classical solvers consume" — but with W1e flagged as a Class-1 calibration-input residual that the framework's outer-factor sector does not autonomously generate. The hybrid pipeline gives access to high-quality classical solvers; the residual wall is a calibration-data gap not a correlation-depth gap.
* **The §IV structure** should NOT be "W1e closure under DMRG"; it should be a section on the Hamiltonian-level limitation. A natural title: "Hamiltonian-specification residual: where the outer-factor sector is silent". This is also where the chemistry-side analog of the Sprint H1 Yukawa non-selection theorem most naturally lands.

### §3.4 Implications for Phase 2 hybrid-pipeline arc

Phase 2 was scoped (per `debug/sprint_hybrid_pipeline_scoping_memo.md`) around AFQMC via ipie on Cholesky-decomposed `eri`. R3-B's verdict for NaH propagates: AFQMC on the same balanced FCIDUMP cannot close W1e either. The natural Phase 2 sharpening is to test:

1. **Whether the wall is structurally identical for HCl, MgH₂, SiH₄, PH₃, H₂S** (the other [Ne] frozen-core hydrides). If yes, the wall has cross-system structure and a non-trivial dual question opens: "what is the chemistry-side analog of the η-trivialisation theorem (Paper 18 §IV.6 SM-side) that classifies this residual?"
2. **Whether enlarging the active space** (promoting the [Ne] core to active) opens a different W1e regime where Class-3 correlation is non-trivial. This is a $\sim 10^9$-determinant FCI for NaH and is the natural Phase-3 target for DMRG via block2 / ipie (Phase 2 multi-month, deferred).

R3-B does NOT close the Phase-2 question wholesale — it closes only the *naive-pipeline* question for the *production* balanced FCIDUMP.

---

## §4. Recommended paper edits (do NOT apply)

For PI review at sprint close:

**Recommendation 1 — Paper 18 §IV.6 chemistry-side analog text.** Append a paragraph after the §3.2 Sprint W1e period-class block (per Sprint W1e §4 proposed text) along the following lines:

> The operational test of the chemistry-side analog (Sprint R3-B NaH DMRG-on-FCIDUMP falsifier, 2026-06-07, `debug/sprint_r3b_dmrg_nah_falsifier_memo.md`) ran the framework FCI on the balanced NaH FCIDUMP at $n_{\max} = 2$ ($\dim_{\text{FCI}} = 100$) and $n_{\max} = 3$ ($\dim_{\text{FCI}} = 784$) — at these sizes the FCI exhausts all determinant content of the balanced Hamiltonian, so the DMRG-at-infinite-$\chi$ limit equals the framework FCI exactly. The PES is monotone ascending across $R \in [2.5, 6.0]$ bohr at both $n_{\max}$ values; no interior minimum near $R_e^{\text{exp}} = 3.566$ bohr exists, and the framework over-binds by $-1.691$ Ha ($n_{\max}=2$, $23.7\times$ exp $D_e$) / $-1.502$ Ha ($n_{\max}=3$, $21.1\times$ exp $D_e$) at the same FCIDUMP-level Hamiltonian. This is the chemistry-side analog at the operational level: the framework outer-factor sector composes external atomic-physics calibration data cleanly via the FCIDUMP interface, but no classical-solver post-processing can close the wall because the wall is in the input data tier, not in the determinant expansion.

**Recommendation 2 — Paper 57 abstract structure.** Reshape around the verdict. The abstract should NOT include "W1e closure under DMRG"; it should include "W1e is Hamiltonian-specification residual not solver-closable" as a section heading or structural claim. Paper 57 is currently undrafted in the corpus listing; this recommendation is for the drafting outline at sprint close.

**Recommendation 3 — CLAUDE.md §3 row append.** New row in the failed-approaches table:

> | Approach | Count | Lesson |
> |:---|:---:|:---|
> | DMRG-on-balanced-FCIDUMP as W1e closure for NaH at $n_{\max} \in \{2, 3\}$ (Sprint R3-B, 2026-06-07) | 1 | At production NaH balanced sizes, $\dim_{\text{FCI}} \leq 784$, so DMRG-at-infinite-$\chi$ = framework FCI; bit-identical to P4 baseline. W1e wall is in Hamiltonian specification, not solver-accessible. See `debug/sprint_r3b_dmrg_nah_falsifier_memo.md`. |

**Recommendation 4 — CLAUDE.md §1.7 W1e classification update.** The Sprint R3-B verdict is empirical confirmation that W1e is the **6th and sharpest instance of the multi-focal-composition wall pattern**: the framework can run external classical solvers on its own integrals, but cannot derive the missing calibration content from those solvers. Update §1.7 multi-focal-composition wall pattern bullet to reference Sprint R3-B as the operational confirmation.

---

## §5. Honest scope

What this sprint demonstrated:
* The Sprint P1 FCIDUMP exporter is operationally verified by an independent FCI solver consuming its output: the round-trip and the FCI calculation on the parsed integrals agree bit-identically with the framework's own FCI ($\leq 10^{-12}$ Ha) and with the P4 published baseline.
* At production NaH balanced parameters, the DMRG-at-infinite-$\chi$ limit (equivalently, sector-restricted FCI) on the balanced FCIDUMP does NOT produce an interior minimum near $R_e^{\text{exp}}$ at $n_{\max} \in \{2, 3\}$; the artifact-well configuration of P4 §2 survives unchanged.
* The wall direction and magnitude are bit-identical to P4: $-1.69$ Ha / $-1.50$ Ha at $R_e^{\text{exp}}$ vs $R = 8$ anchor; W1e over-binds by $21\text{--}24 \times$ experimental $D_e$ regardless of solver depth.

What this sprint did NOT demonstrate:
* **Block2 / pyscf DMRG on a system too large for direct FCI.** Neither block2 nor pyscf wheels installed under Python 3.14 (build failures on missing BLAS / CMake configuration); the sprint adopted the operational substitute that for the production NaH balanced FCIDUMP, FCI IS the DMRG-at-infinite-$\chi$ limit and the question is determined.
* **Active-space promotion of the [Ne] frozen core.** A 10-electron $\dim_{\text{FCI}} \sim 10^9$ calculation would test the Class-3 contribution beyond the production architecture. This is the natural Phase 2/3 deferred target.
* **Cross-system robustness on other [Ne] frozen-core hydrides.** HCl, MgH₂, SiH₄, PH₃, H₂S would each need an independent F4/F5/F6 + R3-B sprint cycle.
* **$n_{\max} = 4$ on NaH balanced.** Per F6, $n_{\max} = 4$ ($Q = 120$, $\dim_{\text{FCI}} \approx 14400$) achieves $10.2\%$ PES closure but remains in the artifact-well regime; an R3-B-style FCIDUMP-loop at $n_{\max} = 4$ is the natural next probe if the Phase 2 question is reopened.
* **AFQMC / CCSD(T) cross-check.** Both unblocked by FCIDUMP per P1 §4. Not exercised in this sprint scope; the structural argument (any classical correlation solver $\leq$ FCI on the same Hamiltonian) is the load-bearing reason these are not informative for the current load-bearing question at NaH n_max≤3.
* **R = 2.0 and R = 8.0 endpoint extensions.** R3-B's grid is the 7-point inner subset of P4's 9-point grid. The wall direction is unambiguous from the 7-point monotone-ascending pattern; the 2-point endpoint extension would deepen the artifact-well diagnostic (already at $-170.51$ Ha at $R=2$ per P4) without changing the verdict.

Decision-gate outcome: **STOP. Falsifier triggered.** The Class-1 reading of W1e is structurally vindicated by the negative result; the Class-3-via-DMRG closure hypothesis is empirically falsified at the production NaH balanced FCIDUMP.

---

## §6. Files

### Created
* `debug/r3b_dmrg_nah_falsifier_driver.py` (~280 lines): two-path FCI driver. Path A loops FCIDUMP through the Sprint P1 exporter/reader and runs `coupled_fci_energy` on the parsed integrals; Path B runs `coupled_fci_energy` on the in-memory result; both compared to the P4 published baseline at each grid point.
* `debug/data/r3b_dmrg_nah.json`: full numerical record — per-$R$ direct + FCIDUMP-loop energies, two-path diff, FCIDUMP round-trip diagnostics, comparison to P4 baseline, interior-minimum search outcome, wall quantification, verdict classification.
* `debug/data/r3b_nah_balanced_nmax2_R*.fcidump` (7 files, $n_{\max}=2$): one Knowles-Handy FCIDUMP per $R$-grid point, suitable for any external classical-correlation solver. Sizes $\sim 5$ KB each.
* `debug/data/r3b_nah_balanced_nmax3_R*.fcidump` (7 files, $n_{\max}=3$): same at $n_{\max} = 3$. Sizes $\sim 80$ KB each.
* `debug/sprint_r3b_dmrg_nah_falsifier_memo.md` (this memo).

### NOT modified
* Production `geovac/` modules — diagnostic-only sprint per sprint mandate; no physics changes.
* Tests — no production code modified; regression preserved.
* Papers — recommendations in §4 above for PI review at sprint close; not applied here per task instructions ("DO NOT MODIFY: production code; driver only").
* CLAUDE.md — recommended row append + §1.7 update for PI review at sprint close; not applied here.

---

**End of Sprint R3-B memo. Verdict: STOP. DMRG-on-balanced-FCIDUMP does NOT close W1e on NaH at production parameters; the wall is at the Hamiltonian specification level, not in the determinant expansion. The cosmic-Galois Class-1 reading of W1e is structurally vindicated.**
