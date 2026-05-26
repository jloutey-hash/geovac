# Track AdS-B — Diagnostic-before-engineering audit on Ryu-Takayanagi-style holographic entanglement entropy on the discrete S³ graph

**Date:** 2026-05-25
**Sprint:** Track AdS-B (read-only diagnostic; no production code modified)
**Builds on:** Sprint BH-Phase0 (2026-05-22, `debug/bh_phase0_diagnostic_memo.md`); Sprint TD Track 5 (2026-05-09, `debug/sprint_td_track5_memo.md`); Bekenstein–Hawking sprint scoping memo (2026-05-16, `debug/bekenstein_hawking_sprint_scoping_memo.md`); Papers 27, 42–46.

---

## Verdict

**BLOCKED-AS-RT, but GO-FAST as Cardy–Calabrese-style discrete boundary entropy (3–6 weeks).**

The naive "discrete Ryu–Takayanagi on the S³ Fock graph" reading is structurally blocked by three independent obstructions documented below. The Cardy–Calabrese-flavor closure — log-scaling entropy with slope = boundary dimension, computed bit-exactly on the operator-system wedge — **has already been delivered in Sprint BH-Phase0** at the operator-system level. Extending it to a sub-region-as-graph-cut diagnostic (the only remaining substantive question that's a genuine RT-flavored computation) is 3–6 weeks of work, NOT 3–4 months. The 3–4 month figure in the original scope is dominated by infrastructure that already exists.

This is the structural-skeleton-scope pattern (CLAUDE.md §1.7 / multi-focal wall) repeating in the holography direction: the framework delivers a clean discrete entropy formula, but it is NOT the area-law statement the original RT framing wanted.

---

## §1. What the framework already has for this question

### §1.1 Graph-cut / minimum-cut infrastructure

`geovac/fock_graph_hodge.py` is the natural home, but with a critical structural constraint that limits the RT analog:

> The Fock graph has $\beta_0 = n_{\max}$ connected components, one per angular momentum sector $l = 0, 1, \ldots, n_{\max}-1$. T± connects same-$l$ states across shells; L± connects same-$(n, l)$ states across $m$ values; no transition changes $l$, so different $l$-sectors are disconnected. (`fock_graph_hodge.py` docstring §Betti number)

**This is a hard structural blocker for a vanilla RT formula.** RT entropy = (minimum cut between sub-region and complement). On a graph with $n_{\max}$ disconnected components, the minimum cut between any spatial sub-region $A$ and its complement $A^c$ is **identically zero** if $A$ and $A^c$ live in different connected components, and reduces to a single-$l$-sector min-cut otherwise. There is no graph-theoretic "boundary" between sub-regions on the S³ Fock graph in the RT sense.

What `fock_graph_hodge.py` does provide cleanly:
- Signed incidence matrix $B$ (V × E) and node/edge Laplacians $L_0 = BB^T$, $L_1 = B^T B$.
- Per-$l$-sector decomposition — substantially useful for per-sector entropy.
- π-free certificate for both $L_0$ and $L_1$ spectra (all eigenvalues algebraic integers).

What it does NOT provide:
- A weighted-edge version where T± and L± have different weights and an "effective boundary" could carry RT-like weight.
- A bipartite graph cut function (no `min_cut(A)` method).
- Any notion of bulk vs boundary partition.

The connected-Dirac-graph Rule B variant (Paper 29 §RH-C) does NOT have the per-$l$ disconnection (E1 dipole couples $\Delta l = \pm 1$). It is the natural substrate IF a graph-cut RT version were to be attempted — but see §3.4 below for why even Rule B is the wrong object.

### §1.2 Reduced-density-matrix machinery

`geovac/casimir_ci.py` produces FCI ground-state wavefunctions for He/Li⁺ at $n_{\max} = 2{-}9$. The spatial 1-RDM and its von Neumann entropy `S_full(GS)` are computed in `debug/entanglement_geometry.py` and (high-precision) in `debug/sprint_td_track5.py`. Paper 27 documents the entropy decomposition (one-body = entanglement-inert; V_ee carries all correlation entropy).

The wedge KMS state $\rho_W = e^{-K_\alpha^W} / Z$ on the truncated Camporesi–Higuchi spectral triple lives in `geovac/modular_hamiltonian.py` (Paper 42 / 43 / 45 substrate). Sprint BH-Phase0 used it for the entropy diagnostic.

### §1.3 What's missing

- No discretized sub-region structure on $S^3$ (e.g., spherical-cap discretization with explicit interior/exterior labels).
- No discretized minimum-cut function over a weighted Fock graph.
- No bulk-vs-boundary partition function.
- No definition of a "discrete CFT₃ on the Fock graph" against which RT would be tested.

---

## §2. Standard RT in CFT₃ on round S³ — the continuum baseline

The Casini–Huerta–Myers (CHM) construction (arXiv:1102.0440) gives the exact continuum RT for a spherical sub-region of radius $R$ in a CFT on round $S^d$ at vacuum: the reduced density matrix is a thermal state on hyperbolic space $H^{d-1}$ at temperature $T = 1/(2\pi R)$, and

$$S_{\text{RT}}(A) = \frac{c}{6} \log\left(\frac{2R}{\epsilon}\right) \quad (d = 2) \qquad \text{or} \qquad S_{\text{RT}}(A) = \frac{\text{Area}(\gamma_A)}{4 G_N} \quad (d \geq 3)$$

For a hemisphere of $S^3$ in CFT₃, the leading divergence is the AREA-law $S \sim A_{\partial A}/\epsilon$ (proportional to the equator $S^2$ area at the UV cutoff scale). The renormalized RT entropy of a hemisphere is finite and computable from CHM.

**What does this predict for the discrete Fock graph?** It predicts $S$ should scale with the dimension of the equator $S^2$, which on the truncated S³ at level $n_{\max}$ has node count $n_{\max}(n_{\max} + 1)$. The continuum CFT₃ wants $S \sim n_{\max}^2 / \epsilon$ at the UV cutoff $\epsilon \sim 1/n_{\max}$, giving $S \sim n_{\max}^3$ — volume-law on the discrete graph. **This is precisely the cubic scaling Sprint BH-Phase0 falsified at $R^2 = 0.99996$ versus the log-scaling fit's $R^2 = 0.99991$ with one fewer parameter.**

Conclusion: the continuum CFT₃-RT expectation IS the BH-Phase0 falsified hypothesis (F1 in the BH memo). The naive RT-on-the-Fock-graph reading is dead at this same scope.

---

## §3. Is chemistry entropy the same KIND of object as CFT entanglement entropy? (The critical question)

**NO, at three independent levels.**

### §3.1 Master Mellin engine partition (Sprint TD Track 5)

Sprint TD Track 5 (`debug/sprint_td_track5_memo.md`) closed this in May 2026 via a PSLQ probe at 12,312-element basis × 100 dps. GeoVac correlation entropy `S_full(GS)` for He and Li⁺ at $n_{\max} = 3, 4$ is **PSLQ-disjoint from the master Mellin engine ring** (M1 ∪ M2 ∪ M3) at the tested precision. CFT entanglement entropy in the CHM construction lives squarely in M2 (Seeley–DeWitt heat-kernel mechanism); GeoVac correlation entropy does not.

This is decisive at the algebraic level: the two objects are not in the same transcendental class. Whatever GeoVac correlation entropy is, it is not the Seeley–DeWitt CFT entanglement entropy of a sub-region.

### §3.2 State substrate

CFT entanglement entropy of a sub-region $A$ is $S(\text{Tr}_{A^c} |\Omega\rangle\langle\Omega|)$, where $|\Omega\rangle$ is the vacuum and $A^c$ is the spatial complement. The wedge KMS state of Papers 42–46 corresponds to this in the BW limit ($|\Omega\rangle$ = vacuum, $A^c$ = the wedge complement on the Cauchy slice).

GeoVac correlation entropy of Paper 27 is $S(\text{Tr}_{\text{electron 2}} |\Psi_{\text{GS}}\rangle\langle\Psi_{\text{GS}}|)$, where $|\Psi_{\text{GS}}\rangle$ is the FCI ground state of a Coulomb 2-electron system and the trace is over *one electron*, not a spatial sub-region. This is the per-particle reduced density matrix, not a spatial-bipartition state.

These are different traces. The first is geometric (spatial-region partition); the second is particle-statistical (one-electron-out-of-two partition). Sprint TD Track 5's PSLQ negative is the algebraic confirmation that they live in categorically different rings.

### §3.3 Sprint BH-Phase0's actual finding

Sprint BH-Phase0 computed $S(\rho_W)$ for the spatial wedge KMS state — the right object for RT analogy — and found $S \approx 2 \log(n_{\max})$ at $R^2 = 0.99991$. This is **Cardy–Calabrese**-style log scaling with slope = 2 (the dimension of the wedge boundary $\partial W = S^2$). It is NOT RT area-law.

The slope = boundary dimension is structurally meaningful (it's the discrete operator-system analog of $S = (c/3) \log L$ for a 1+1D CFT on a region of length $L$), but the scaling is logarithmic, not power-law, because the wedge KMS state is effectively concentrated on the lowest $K_\alpha$ shell (the equator). The Hilbert-space dimension grows as $n_{\max}^3$, the wedge dimension as $n_{\max}^3 / 2$, but the entropy only sees the boundary shell.

**This is a positive structural finding, but it is NOT RT.** It is closer to Cardy–Calabrese boundary entropy or to the Page-curve plateau of black hole evaporation than to the RT minimal-surface entropy of a sub-region.

### §3.4 Why the Dirac Rule B graph doesn't save the situation

Rule B is the connected graph option (E1 dipole, $\Delta l = \pm 1$, no per-$l$ disconnection). One might hope to compute min-cut entropy on Rule B as a discrete RT analog. Three problems:

1. The Rule B graph is the Dirac-spin-graph (Paper 29 §RH-C), not a discretization of the spatial geometry. Edges represent dipole transitions in the (n, l, m) eigenbasis, not bulk-vs-boundary spatial neighborhood. A min-cut on Rule B partitions transitions, not spatial points.
2. Min-cut RT requires a notion of "boundary" — which sub-region on Rule B corresponds to a hemispheric cap on $S^3$? The (n, l, m) eigenstates are spatially delocalized; the spectral basis has no spatial-locality structure that would let us identify "sub-region A".
3. Even if a sub-region were defined, the min-cut would be Hodge-theoretic ($\beta_1 = E - V + c$ for the cut graph), which gives integer counts — combinatorial, no Newton's constant calibration available.

This is the same structural reason CFT₃-RT does not transport to GeoVac: the framework operates in the spectral basis (energy/angular momentum eigenstates), not the spatial basis (points on $S^3$). The Fourier transform between the two is the Casimir spectrum, and the wedge KMS state's entropy reflects spectral structure, not spatial geometry.

---

## §4. HaPPY-code / tensor network framing: is GeoVac a candidate discrete holographic TN?

**Tempting analogy, but the structural answer is NO at the load-bearing level.**

HaPPY codes (Pastawski–Yoshida–Harlow–Preskill 2015, arXiv:1503.06237) build a tensor network on a hyperbolic tiling with perfect tensors at each node. The RT formula is **exact on a HaPPY code** because the perfect-tensor property forces the minimum cut to equal the entanglement entropy of the boundary sub-region. Discrete RT on HaPPY is a theorem.

The Hopf bundle $S^3 \to S^2$ and the Bargmann–Segal $S^5 \to \mathbb{CP}^2$ have explicit fiber structure, but:

1. **Wrong manifold.** HaPPY codes live on hyperbolic tilings (e.g., {5, 4}); GeoVac graphs live on $S^3$ (positively curved). The hyperbolic-tiling structure is what makes the RT formula go through — the negative curvature gives the bulk geometry that minimal surfaces sit in. $S^3$ is the boundary, not the bulk, in the AdS₄/CFT₃ correspondence. Bringing in $H^4$ or $\mathbb{H}^3$ would require an entirely new construction.
2. **No perfect-tensor structure at GeoVac nodes.** The Wigner 3j coefficients at Fock-graph vertices are NOT isometries from input legs to output legs (they're CG-projection tensors with specific selection rules). They do not satisfy the "any half of the legs is an isometry" property that makes HaPPY codes' RT formula a theorem.
3. **Fiber structure ≠ bulk geometry.** The Hopf fibration's $S^1$ fibers index a U(1) gauge structure (Paper 25), not a bulk radial direction. A holographic bulk would need a manifold-with-boundary structure where the bulk is one-dimension-higher than the boundary; GeoVac has fibration over the same-dim base, not over a lower-dim boundary.

A genuine HaPPY-style construction on the GeoVac geometry would require: (i) lifting to $H^4 \times S^3$ or AdS₄ × S³; (ii) putting perfect tensors on a hyperbolic-tile bulk; (iii) reading the $S^3$ boundary as the CFT side. This is a 6–12 month construction-from-scratch sprint, not a 3–4 month diagnostic.

**The honest answer:** GeoVac has spectral-triple machinery on $S^3$, not holographic-bulk machinery on AdS₄. The two are different objects.

---

## §5. Cleanest single deliverable for a 3–6 week sprint

Given the BH-Phase0 result already in hand and the structural blockers documented above, the cleanest sprint deliverable is:

**Discrete sub-region-on-S³-Fock-graph Cardy–Calabrese check.** Use the per-$l$-sector decomposition of `fock_graph_hodge.py` to partition the Fock graph into (i) "boundary" states (low-$l$, e.g., $l = 0, 1$) and (ii) "bulk" states (high-$l$). Compute the per-sector reduced density matrix of the FCI ground state from `casimir_ci.py`, and verify that $S(\rho_{\text{boundary}})$ scales as $\log(n_{\text{boundary states}})$ at the same slope (≈ 2) as Sprint BH-Phase0.

This is concretely:
- ~1 week to wire `geovac/casimir_ci.py` 1-RDM machinery to accept a sector projection (e.g., $P_{l \leq L_0}$).
- ~1 week to compute the per-sector entropy across $n_{\max} = 2{-}6$ and $L_0 = 0, 1, \ldots, n_{\max} - 1$ panel.
- ~1 week to compare against Sprint BH-Phase0's modular-Hamiltonian entropy and check whether the two log-scaling slopes agree.
- ~1 week for paper edits (Paper 27 §VIII or new subsection, Paper 34 §V.B Cardy-style row) and writeup.

**Expected positive outcome:** the FCI-state per-$l$-sector entropy scales as $\log(n_{\max})$ with slope ≈ 2, matching the wedge KMS state's BH-Phase0 result. This would be the **first cross-check that the framework's two natural entropies (state-side FCI 1-RDM and spectral-side wedge KMS) agree on the boundary-dimension scaling**.

**Expected negative outcome:** the scalings disagree (different slopes), confirming Sprint TD Track 5's PSLQ structural finding that chemistry entropy and spectral-side entropy are categorically different objects at all levels (not just at the master Mellin engine ring level). Either outcome is informative.

**Decision gate for the sprint:** at week 2, run the per-sector entropy on $n_{\max} = 2, 3$. If the slope is in $[1.5, 2.5]$ (matching BH-Phase0 within 25%), proceed to full panel + paper edits. If the slope is < 1 or > 3, write a clean-negative memo and close.

This is **NOT** a Ryu–Takayanagi closure. It is a discrete boundary-entropy closure with one structurally meaningful number (slope = boundary dim).

---

## §6. Recent prior art (2024–2026): scoop risk and citation targets

The web search across "2024 2025 discrete RT", "noncommutative geometry spectral triple holographic", "Casini–Huerta CFT3 lattice", and "modular Hamiltonian JLMS bulk boundary spectral truncation" returned:

- **No published work on discrete RT directly on a spectral-triple operator-system truncation** as of May 2026. The Connes–van Suijlekom 2021 / Hekkelman–McDonald 2024 / Leimbach–van Suijlekom 2024 lineage (extensively cited in Papers 38–47) explicitly does NOT address modular Hamiltonian or entanglement entropy.
- **The 2026 JLMS-with-error-correction paper** (arXiv:2601.00442) is the closest active research direction but lives in the OAQEC framework of AdS/CFT, not in spectral-triple truncations.
- **The Sept 2025 random-tensor-network statistical-mechanics paper** (arXiv:2508.16570) extends the RT-from-random-TN paradigm; not directly relevant to GeoVac's deterministic spectral structure but useful as a citation if a Cardy-flavor entropy paper is drafted.
- **Casini–Huerta 2009 and 2011** (CHM) remain the standard references for the continuum CFT₃ entanglement entropy targets.
- **Zhu–Casini 2020** (`debug/l1_modular_hamiltonian_architecture_memo.md` line 476) on lattice BW modular Hamiltonian on critical chains is the closest precedent for finite-lattice BW reading; relevant comparison for any GeoVac Cardy paper.

**Scoop risk:** LOW. The intersection of (spectral-triple / Connes–vS truncation framework) + (modular-Hamiltonian-entropy or RT-style entropy) is not occupied. The Mondino–Sämann synthetic Lorentzian GH program (extensively covered in Papers 47–49) is the closest adjacent direction; it does not extend to entanglement entropy.

**Citation targets for any sprint output:** CHM 2011 (arXiv:1102.0440) as the continuum baseline; HaPPY 2015 (arXiv:1503.06237) as the discrete-holography contrast; Zhu–Casini 2020 as finite-lattice BW precedent; JLMS 2016 (arXiv:1512.06431) and the 2026 follow-up; the Hekkelman–McDonald 2024 (arXiv:2412.00628) and Leimbach–van Suijlekom 2024 papers for the spectral-truncation Riemannian lineage; Sprint BH-Phase0 and Sprint TD Track 5 as internal precedents.

---

## §7. Honest risk: is GeoVac's structure the right substrate for discrete RT?

**No, not for a faithful RT analog. Yes, for a Cardy–Calabrese-style boundary-dim log-scaling diagnostic.**

The RT formula in CFT₃/AdS₄ requires (i) a CFT on the boundary, (ii) a bulk classical-gravity dual, (iii) a minimal surface in the bulk. GeoVac has none of these:

- GeoVac on $S^3$ is the boundary, with no bulk dual constructed. The Bargmann–Segal $S^5$ is a different system (HO basis, Paper 24), not the bulk of $S^3$ in any AdS sense.
- No graph-theoretic minimal-surface structure exists on the disconnected Fock graph.
- The "CFT" reading of the GeoVac Hopf graph is structurally weak — Paper 32's spectral-triple framing is operator-algebraic, not CFT-stress-tensor-algebra-based.

What DOES survive cleanly is the **operator-system entropy** $S(\rho_W)$ of the wedge KMS state (Sprint BH-Phase0), which is structurally the discrete analog of CFT thermal entropy on the BW Rindler wedge, not RT entropy of a spatial sub-region.

**The honest framing of any sprint deliverable would be:** "Discrete boundary-dimension log-scaling entropy on the truncated S³ Fock graph, computed from both the spectral-side wedge KMS state (BH-Phase0) and the state-side FCI 1-RDM per-$l$-sector reduction (this sprint), with cross-check that the slopes agree at ≈ 2." This is a paper, but it is NOT an RT paper. It is a CC-flavor / spectral-triple-entropy paper.

If the PI wants a literal discrete RT formula, the answer is **BLOCKED** at the level of (i) constructing AdS₄ × S³ in spectral-triple language, (ii) wiring perfect-tensor structure into Wigner-3j vertices, (iii) defining bulk vs boundary on a single-manifold spectral triple. Each is a 6–12 month construction sprint, total 18–36 months. The 3–4 month original scope is not sufficient.

---

## §8. Overlap with Tracks A (partition function) and C (modular JLMS)

The wedge KMS state $\rho_W = e^{-K_\alpha^W}/Z$ from Paper 42 is the **common substrate** for all three tracks:

- **Track A (partition function $Z$):** $Z = \text{Tr}(e^{-K_\alpha^W}) = \sum_{m_j > 0} e^{-2m_j} \cdot g(2m_j) = e^{-1} n_{\max}(n_{\max}+1) + O(e^{-3})$. Closed-form expression at every $n_{\max}$; already implicit in Sprint BH-Phase0 (`debug/data/bh_phase0_entanglement_entropy.json`).
- **Track B (entanglement entropy):** $S(\rho_W) = -\text{Tr}(\rho_W \log \rho_W) = \langle K_\alpha^W \rangle - \log Z + \log Z = \log Z + \beta \langle K_\alpha^W \rangle$ at BW canonical $\beta = 2\pi$. Computed in Sprint BH-Phase0.
- **Track C (JLMS bulk identification of $K_\alpha$):** the modular Hamiltonian $K_\alpha^W = J_{\text{polar}}/(2\pi)$ has integer spectrum two_m_j; $K_\alpha^W$ is already the operator that Paper 42 identifies as the BW boost generator. The JLMS step is to identify $K_\alpha^W$ as a bulk modular Hamiltonian of an enveloping bulk theory — which would require the bulk-dual construction the framework does not have (§7 above).

**Shared infrastructure:** all three tracks reuse `geovac/modular_hamiltonian.py::HemisphericWedge` and `restrict_K_alpha_to_wedge()`. Track A is essentially a 1-line addition to BH-Phase0's compute (extract $Z$ values per panel cell). Track C is the open structural problem (no bulk dual) — multi-month frontier.

**Recommendation for sprint organization:** A combined Track A+B sprint is a clean 4–6 week deliverable. Track A is a 1-week side-task of Track B (compute $Z$, log it, paper-edit). Track C should be separated as a multi-month frontier item, blocked at the bulk-dual-construction level.

---

## §9. Decision summary

| Sub-question | Verdict | Effort |
|:---|:---:|:---:|
| Discrete graph-cut RT on Fock graph | **BLOCKED** | Multi-month (no bulk, no perfect tensors, disconnected graph) |
| Cardy–Calabrese log-scaling on wedge | **DONE (BH-Phase0)** | 0 weeks |
| Per-$l$-sector FCI entropy cross-check | **GO-FAST** | 3–6 weeks (one cross-check deliverable) |
| HaPPY-style holographic TN on Hopf/Bargmann | **BLOCKED** | 6–12 month construction |
| JLMS-style modular bulk identification | **BLOCKED** (no bulk dual) | Multi-month frontier |
| Combined Track A+B partition-function-and-entropy paper | **GO-FAST** | 4–6 weeks |

**Verdict on the original 3–4 month scope:** **NEITHER GO-AS-PROJECTED NOR BLOCKED.** The original scope conflates two different deliverables:

1. **A literal discrete RT formula** — BLOCKED at the bulk-dual-construction level.
2. **A discrete CC-style log-scaling diagnostic** — already 80% done in Sprint BH-Phase0; cross-check sprint is 3–6 weeks.

The 3–4 month estimate would be honest if the goal were to build a HaPPY-style TN on a hyperbolic-tile bulk + S³ boundary, with perfect-tensor structure at vertices — a genuine new construction. But that's structurally a different sprint from the question as framed ("can the framework's existing entropy / graph-cut machinery give an RT formula").

The honest deliverable in 3–6 weeks: a state-side / spectral-side cross-check on the slope-2 log scaling, plus a clean closure memo on what the framework cannot do (literal RT). This sharpens the structural-skeleton-scope statement on the holography axis.

---

## §10. Recommendation

**Do not run an "RT-on-the-graph" sprint as such.** The diagnostic-before-engineering rule says: if ≥ 2 honest negatives accumulate on a question (here: TD Track 5 PSLQ negative, BH-Phase0 area-law negative), the next sprint should be a sharpened diagnostic, not another engineering attempt at the same target. Both negatives point at the same structural fact: the framework's natural entropies (FCI 1-RDM correlation entropy AND wedge KMS modular entropy) are NOT in the same class as continuum CFT RT entropy.

**Do consider running a combined Track A+B "discrete boundary entropy" sprint** (4–6 weeks, single deliverable). This would:
- Compute the partition function $Z(n_{\max})$ closed form (Track A, ~1 week).
- Reframe BH-Phase0's $S(\rho_W)$ result as a Cardy–Calabrese-flavor boundary-dim entropy (Track B core, mostly already done).
- Add a state-side cross-check via per-$l$-sector FCI 1-RDM entropy from `casimir_ci.py` (the new content, ~2 weeks).
- Produce a Paper 27 §VIII or §IX extension and possibly a small standalone math.OA-adjacent note (4–6 page).

**Do flag the literal RT question as multi-month frontier** with the specific blocker being construction of AdS₄ × S³ in spectral-triple language. Pin this to the same "structural-skeleton-scope" list as W1c-residual, W3 calibration data, multi-loop QED — all are real open questions, all are structurally distinct from the framework's current scope.

This is not the answer the original sprint scope was looking for, but it is the honest one: the framework has already delivered the discrete-entropy result that the RT framing wanted; what it cannot do is bulk-dual construction or perfect-tensor TN structure. Those are 18–36 month sprints, not 3–4 month sprints, and they require structurally different infrastructure than the framework currently has.

---

## §11. Files referenced (no modifications)

**Production code (read-only):**
- `geovac/fock_graph_hodge.py` — Hodge decomposition, per-$l$-sector connectedness
- `geovac/casimir_ci.py` — FCI matrix builder, 1-RDM extraction substrate
- `geovac/modular_hamiltonian.py` — `HemisphericWedge`, `for_bisognano_wichmann()`, wedge KMS state machinery
- `geovac/composed_qubit.py` — qubit Hamiltonian substrate

**Prior sprint outputs (read-only):**
- `debug/bh_phase0_diagnostic_memo.md` — 2026-05-22 BH Phase 0 closure; log-scaling slope 2
- `debug/bh_phase0_entanglement_entropy.py` — driver
- `debug/data/bh_phase0_entanglement_entropy.json` — raw data
- `debug/sprint_td_track5_memo.md` — 2026-05-09 PSLQ negative on master Mellin ring
- `debug/data/sprint_td_track5.json` — raw S_full(GS) at 150 dps
- `debug/bekenstein_hawking_sprint_scoping_memo.md` — 2026-05-16 original BH scoping

**Memory files (read-only):**
- `memory/pythagorean_orthogonality.md` — Paper 42 §7.2 / §8 / Paper 43 §10.2 substrate
- `memory/paper27_entropy_thesis.md` — Paper 27 thesis (entropy ≡ V_ee, one-body inert)

**Output:**
- `debug/sprint_ads_track_b_entanglement_audit_memo.md` — this memo (~2400 words)

No production `geovac/`, `papers/`, or `tests/` modifications.

---

## §12. Citation context for the verdict

The structural pattern (framework delivers discrete operator-system result that is structurally distinct from the continuum-CFT version it was hoped to reproduce) repeats across:

- **Sprint H1 (Higgs, 2026-05-06):** AC extension admits Higgs structurally, but Yukawa is calibration input.
- **Sprint LS-8a (2026-05-07):** Iterated CC spectral action reproduces UV-divergent integrand, but Z_2/δm renormalization is external.
- **Sprint HF-3/4/5 (2026-05-07):** Bohr-Fermi + Schwinger reach +18 ppm on 21cm, but recoil/Zemach/multi-loop a_e require external input.
- **W3 / Wolfenstein (2026-05-08):** Mechanical-basis search returned z = −0.52 (null) on master Mellin ring identification.
- **Sprint TD Track 5 (2026-05-09):** GeoVac correlation entropy NOT in master Mellin ring.
- **Sprint BH-Phase0 (2026-05-22):** Wedge KMS entropy log-scales, not power-law-area.
- **(this sprint, 2026-05-25):** RT-style sub-region entropy on Fock graph is BLOCKED at construction level; CC-style boundary-dim log scaling is GO-FAST.

This is the structural-skeleton-scope pattern. The framework's master Mellin engine + spectral-triple machinery + Wick-rotation-functor bridge (Papers 38–49) gives clean discrete operator-algebraic structure, but does not autonomously generate the calibration data (Newton's constant, Yukawa couplings, renormalization counterterms, area-law area coefficient) that links the operator-system structure to continuum geometric expressions.

The "second packing axiom" question — what generates the calibration data — remains the framework's principal open structural problem. RT-style holographic entanglement entropy on a discrete graph would be a calibration question of exactly this kind: the area-law coefficient is precisely the discrete analog of $1/(4G_N)$, which has no autonomous spectral-triple origin in the framework's current structure.
