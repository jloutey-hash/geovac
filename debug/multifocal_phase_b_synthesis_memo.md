# Multi-Focal Composition Phase B Synthesis

**Date:** 2026-05-07
**Author:** PM (Phase B synthesis; no production files modified, no paper edits applied)
**Sources:** All seven Phase B memos read in full —
- `debug/multifocal_b_w1a_diag_memo.md` (~5400 words; verdict (b) HYBRID, algebraic backbone closed; Roothaan 1951 cross-register ERI is the headline)
- `debug/multifocal_b_w1b_diag_memo.md` (~3400 words; verdict (b) downstream of W1a; augmented sprint closes both)
- `debug/multifocal_b_w1c_diag_memo.md` (~3500 words; verdict (b) tooling-addressable with non-trivial extension; 10× overattraction at NaH R=3.5 verified)
- `debug/multifocal_b_w2a_diag_memo.md` (~5400 words; verdict (b)+(c) hybrid + frontier-of-field; HM 2024 is Tauberian residue of master Mellin engine; no published framework gives counterterms; §VIII edit drafted)
- `debug/multifocal_b_w2b_diag_memo.md` (~4100 words; verdict (a) keystone reachable for W2b-easy; W2b-medium/hard structurally blocked; strategic restructure proposed)
- `debug/multifocal_b_w3_diag_memo.md` (~3700 words; verdict (b) own wall but no concrete proposal; WH4 deflated to one Fock-projection statement plus three forced consequences)
- `debug/multifocal_b_position_memo.md` (~4900 words; T1 Track NI gap confirmed; T2 Zenodo memo not Paper 39; T3 §VIII.D LaTeX paste-ready ~970 words)

Also re-loaded: `debug/multifocal_phase_a_synthesis_memo.md` for continuity.

---

## Section 1: Refined wall taxonomy after diagnostics

**[Established by Phase B.]** All seven diagnostics returned. The wall taxonomy that Phase A sketched at the verdict-pending level is now classified at the verdict-conditional level. The key finding is that the panel reduced more cleanly than Phase A predicted, and one wall (W2b) bifurcated into a tractable easy case and structurally-obstructed medium/hard cases that need to be tracked separately.

| Wall | B-diagnostic verdict | Frontier classification | Phase C effort | Dependencies |
|:-----|:---------------------|:------------------------|:---------------|:-------------|
| **W1a** (cross-register coordinate operator) | **(b) HYBRID** — algebraic backbone closed (Roothaan 1951 + multipole termination at $L_{\max} = l_a + l_b$ trivially preserved across mismatched exponents). Engineering-only. | GeoVac-internal, tooling-addressable. The "most likely failure mode" Phase A flagged (multipole termination) was a phantom; the real risk is engineering and validation. | 6–10 weeks standalone; **4–6 weeks streamlined** if downstream of W2b-easy keystone. | Independent; *enables* Phase C-W1a-physics if executed second. |
| **W1b** (magnetization-distribution) | **(b) reduces to W1a.** Operator infrastructure shared (same multipole expansion across operator-valued $\hat{\mathbf{R}}_p$); $r_Z$ adds as one Layer-2 calibration scalar from Eides 2024. | GeoVac-internal, tooling-addressable. Different inner-fluctuation component on the same composed triple, analogous to $\omega_\text{gauge}$ vs $\omega_\text{Higgs}$ in the Connes SM construction. | **Folded into Phase C-W1a (no separate sprint).** Incremental cost: one inner-fluctuation component + one Layer-2 scalar. ~100–200 additional cross-register Pauli terms. | Folded into Phase C-W1a-physics. |
| **W1c** (frozen-core cross-center screening) | **(b) tooling-addressable with non-trivial extension.** Quantitatively confirmed: 10× overattraction at NaH R=3.5 bohr matches Sprint 7 PES failure mechanism. Multipole termination preserved; mechanical extension of `_split_integral_analytical`. | GeoVac-internal, narrowest in scope of W1 sub-walls. **Independent of W2b/W1a** — it is a single-integral-evaluation extension, no architectural change. | 5–7 weeks first pass. Risk profile: low. | **Independent.** Can run in parallel with W2b-easy and W1a. |
| **W2a** (multi-loop UV/IR composition) | **(b) hybrid + (c) frontier-of-field.** HM 2024 is the Tauberian residue of the master Mellin engine on Dirac-$S^3$ at $s = d/2$ ($\sqrt{\pi}$ M2 cancels $\Gamma$ normalization, recovers $2/3$). Marcolli–vS rationality constrains the ring. **Neither generates counterterms.** | Frontier-of-field. No published spectral-action framework solves multi-cutoff renormalization autonomously. | **Not a Phase C sprint.** Paper 32 §VIII frontier-of-field framing edit only (~250 words drafted). | None — paper edit is independent positioning. |
| **W2b-easy** ($\mathcal{T}_{S^3}^{\lambda_a} \otimes \mathcal{T}_{S^3}^{\lambda_b}$) | **(a) keystone reachable.** Each of Paper 38's five lemmas has a structurally clean tensor extension (L1' chirality $\gamma_a \otimes \gamma_b$; L2 product Plancherel; L3 Connes–Marcolli joint Lipschitz $C_3 \le 2$; L4 $B_a \otimes B_b$; L5 tunneling-pair joint reach). | NCG-frontier (genuinely open in published literature; Track 3 surprise S2 confirmed via WebFetch on five papers). | 4–8 weeks symbolic + 2–4 weeks numerical + 2 weeks paper draft = **8–14 weeks total**. | Independent. **Strategic keystone if executed.** |
| **W2b-medium** ($\mathcal{T}_{S^3} \otimes \mathcal{T}_\text{Hardy}(S^5)$) | **(d) STRUCTURALLY BLOCKED.** Coulomb/HO category mismatch at the operator-type level. Bargmann transform and Kustaanheimo–Stiefel do NOT bridge spectral triples. | NCG-framework-extension required. Multi-paper / multi-year if at all. | Not a sprint target. | Open structural question; record in Paper 32 §VIII.C (already done). |
| **W2b-hard** ($\mathcal{T}_{S^3} \otimes \mathcal{T}_{S^5}^\text{Riem}$) | **(d) prerequisite-blocked.** Riemannian $S^5$ propinquity convergence is unpublished (Paper 38 covers only $\mathrm{SU}(2) = S^3$ rank-1 case). | NCG-frontier; lower payoff than W2b-easy. | Not a sprint target. | Would need a separate "Paper 38b for $S^5$" before tensor product is even attempted. |
| **W3** (inner-factor parameter selection) | **(b) own wall, no concrete proposal.** 14 catalogued speculations; 0 concrete; 2 candidates (Phase 4H Tracks SM-B, SM-C) tested and falsified. | Calibration-side. May be permanent. | **Not a closure sprint.** Paper 18 §IV.6 framework's-open-question addition (~150 words) only. | None. |

**Two compressed observations from the table:**

1. **The W1 cluster compressed.** Phase A treated W1a, W1b, W1c as three independent walls. Phase B reduced this to **one architectural sprint (W1a-physics) absorbing W1b, plus one independent mechanical fix (W1c)**. The infrastructure shared between W1a and W1b is not a coincidence — both walls require promoting the Track NI proton register to a spectral triple in its own right and building cross-register operators as multipliers on the joint algebra. Whether that promotion is the recoil $V_{eN}$ (W1a) or the magnetization convolution $\rho_E \star \rho_M$ (W1b) is a choice of inner-fluctuation component, not a separate construction.

2. **W2b bifurcated.** Phase A treated W2b as a single frontier-of-NCG wall; Phase B revealed three sub-cases with categorically different verdicts. **W2b-easy is reachable**; W2b-medium is structurally blocked by the Paper 24 §V Coulomb/HO asymmetry; W2b-hard would need a separate $S^5$ propinquity paper before being attempted. Treating W2b as one wall obscures that the sprint-tractable case sits inside the same wall as the structurally-obstructed cases.

---

## Section 2: The strategic restructure

**[Recommendation, not yet established.]** W2b-diag's headline strategic claim is that the multi-focal-composition sprint should restructure around **W2b-easy as the keystone**, with W1a-physics as a streamlined application built on the W2b-easy NCG foundation. This replaces Phase A's Section 4 ranking (W1a as Candidate 1, 6–10 weeks standalone).

### Phase A's ordering (atomic-physics-first)

- **Phase C-W1a Pachucki port** as standalone Candidate 1, 6–10 weeks, novel-math sprint at NCG-atomic-physics intersection.
- Spectral-triple framing of the cross-register architecture is a "polish" deferred to Paper 32 §VIII.C addendum or a candidate Paper 39.
- Single publishable artifact: GeoVac's first explicit multi-focal composition theorem at the spatial-coordinate level, calibrated against Pachucki–Patkóš–Yerokhin 2023.

### Phase B-W2b-diag's restructure (keystone-first)

- **Phase C-W2b-easy as keystone first**, 4–8 weeks symbolic + 2–4 weeks numerical, NCG theorem extending Paper 38 from one $\mathrm{SU}(2)$ factor to a tensor product of two factors.
- **Phase C-W1a-physics second**, streamlined to 4–6 weeks because the algebraic foundation is already laid by the W2b-easy theorem (cross-register operator is a multiplier on the joint algebra; multipole-expansion termination on $S^3 \times S^3$ inherits the joint Gaunt rule).
- Two publishable artifacts: (a) NCG-frontier theorem closing one of the two-infinite-metric tensor-product cases that Track 3's S2 surprise documents as genuinely open; (b) atomic-physics application validated against Pachucki et al. 2023.
- Total effort: 8–14 weeks. Comparable to Phase A's 6–10 week standalone but with two artifacts instead of one.

### Pros / cons

| Dimension | Atomic-physics-first (Phase A) | Keystone-first (Phase B) |
|:----------|:--------------------------------|:--------------------------|
| Total effort | 6–10 weeks | 8–14 weeks |
| Published artifacts | 1 | 2 |
| First-deliverable scope | Atomic-physics calibration | NCG theorem |
| Risk of W2b-easy hidden obstacle | N/A (not attempted) | Medium-low (joint L3/L5 unpublished but reachable) |
| W1a-physics deferred risk | Direct attempt; medium risk in engineering only | Streamlined attempt after foundation laid; lower engineering risk |
| Time-to-first-publication | Faster (6–10 weeks) | Slower (4–8 weeks for W2b-easy + paper draft) |
| Lineage payoff | Pachucki port is a calibration win | NCG theorem closes a published-open-question (Track 3 S2) |
| WH1 register | Stable | New sibling theorem (WH1-class headline) |

### Recommendation

**[Recommendation, PI decision.]** **Keystone-first is well-defended** and is the right move for the project's positioning. Three reasons:

1. **W2b-easy is genuinely open in the published NCG literature** (Track 3 S2; Phase B-W2b-diag §2 verified via WebFetch on Latrémolière 2026, Farsi–Latrémolière 2024 + 2025, Aguilar 2019, Hekkelman–McDonald 2024). Closing it is a publishable contribution at the NCG / spectral-action frontier, parallel to Paper 38's WH1 closure.
2. **The engineering of W1a-physics is non-trivial but the algebra is closed.** Phase B-W1a-diag verified that Roothaan 1951's $J_0(\lambda_e, \lambda_n) = \lambda_e \lambda_n (\lambda_e^2 + 3\lambda_e \lambda_n + \lambda_n^2)/(\lambda_e + \lambda_n)^3$ is the relevant closed form, and Phase B-W2b-diag's roadmap shows that the joint operator system on $S^3 \times S^3$ provides exactly the algebraic foundation W1a-physics needs. Splitting the deliverable means each sprint stays inside its natural specialty.
3. **Phase A's worry about W1a's algebraic risk is now refuted.** Phase B-W1a-diag explicitly retracted Phase A's "medium-high risk" assessment to "medium risk; the algebraic backbone is closed; the engineering is the bottleneck." This makes the keystone-first ordering safer, not riskier — the sequencing protects the smaller second sprint from any unanticipated mathematical obstruction at the foundation.

**Weakness check.** The strategic restructure does have one defensible weakness: if Phase C-W2b-easy hits an unanticipated joint L3 or L5 obstacle (e.g., Bożejko–Fendler product cb-norm fails, or Latrémolière 2017/2023 propinquity height bookkeeping does not lift cleanly to the tensor product), the W1a-physics sprint is delayed without a foundation. The mitigation is that the W1a-physics sprint can be rescoped back to its standalone Phase A form on short notice — Phase B-W1a-diag's task list (sub-tracks C-W1a.1 through C-W1a.4) is intact and runnable as a 6–10 week standalone if W2b-easy stalls. **The restructure is recoverable.**

---

## Section 3: Phase C dispatch plan

**[Recommendation.]** Phase C decomposes into six sprints / artifacts. Three are dispatchable in parallel; two are blocking-sequential after the keystone; one is a small positioning move.

### Sprint table

| Sprint | Target wall | Deliverable | Effort | Dependencies | Zenodo deposit value |
|:-------|:------------|:------------|:-------|:-------------|:---------------------|
| **Phase C-W2b-easy** | W2b-easy (NCG keystone) | Theorem: $\Lambda(\mathcal{T}_{n_a, n_b}, \mathcal{T}_{S^3} \otimes \mathcal{T}_{S^3}) \le C_3 \max(\gamma_{n_a}, \gamma_{n_b})$ via tensor extension of Paper 38's five lemmas. Numerical verification at $(n_a, n_b) \in \{2,3,4\}^2$. Paper draft (standalone "Paper 38b" or §VI of Paper 38 v2). | 8–14 weeks (4–8 symbolic + 2–4 numerical + 2 paper draft) | None | **High.** Closes one of the two-infinite-metric tensor-product cases that Track 3's S2 surprise documents as genuinely open in NCG. WH1-class headline. |
| **Phase C-W1a-physics** | W1a (+W1b folded) | Module `geovac/cross_register_vne.py` (bilinear ERI builder, HO + Sturmian); promote `R_PROTON_BOHR` to operator-valued $\hat{\mathbf{R}}_n$; magnetization-density inner-fluctuation component for W1b; validate against Pachucki–Patkóš–Yerokhin 2023 PRL recoil + Eides 2024 hyperfine. | 4–6 weeks (streamlined) or 6–10 weeks (standalone if W2b-easy stalls) | **Sequential after W2b-easy** for streamlined version; standalone fallback ready. | Medium-high. Atomic-physics application of NCG framework; calibrated against published literature. Resource benchmarks for hydrogen 1S Lamb shift recoil, 21 cm hyperfine. |
| **Phase C-W1c** | W1c (frozen-core cross-center) | Module `geovac/cross_center_screened_vne.py` (extend `_split_integral_analytical` to multi-shell exponentials); wire into `balanced_coupled.py` behind `screened_cross_center=True` kwarg. NaH/MgH$_2$/HCl/H$_2$S/PH$_3$/SiH$_4$ PES regression. Extend to 14 frozen-core species in the library. | 5–7 weeks | **Independent.** Can run in parallel with W2b-easy and W1a-physics. | Medium. Recovers PES capability for second-row balanced builders; closes one of three currently-empirical Sprint 7 negative results. |
| **Track NI Zenodo deposit** | None — positioning | ~1500–2500 word standalone Zenodo memo at `papers/observations/track_ni_atom_spectral_triple.tex` claiming "an explicit … construction" not "the first." Add cross-reference paragraphs to Paper 23 §VI.5 and Paper 32 §V. | ~1 week (memo) + ~1 day (paragraph additions) | None — small positioning sprint. | Medium-high. T1 confirmed published NCG ecosystem does not contain a real-space multi-particle Connes-style spectral triple. Earns a DOI for a structural observation that surrounds the entire Phase C plan. |
| **Paper 32 §VIII.D + citation update** | W2a + framing | Insert ~970 words drafted in B-position T3 (paste-ready LaTeX) as new §VIII.D after §VIII.C. Add three bibitems (`hekkelman_mcdonald2024`, `latremoliere2026`, `paper36`). Flag existing `hekkelman2024` placeholder for cleanup. | 1–2 days (paste + sign-off) | None — positioning, no compute. | Low (positioning only); but tightens Paper 32's lineage placement and names W2a as frontier-of-field rather than GeoVac-internal failure. |
| **Paper 18 §IV.6 W3 framework's-open-question addition** | W3 | Insert ~150 words drafted in B-W3-diag §8 after the existing inner-factor input data tier. Names the second-packing-axiom question as open without committing to a closure attempt. | 1–2 days (paste + sign-off) | None | Low (positioning only); but completes the six-tier taxonomy with an honest scope statement. |

### Parallelization map

```
                                         time →
 ┌────────────────────┐
 │ Track NI Zenodo    │  (~1 week, parallel)
 └────────────────────┘
 ┌────────────────────────────────┐
 │ Paper 32 §VIII.D paste         │  (~1–2 days, parallel)
 └────────────────────────────────┘
 ┌────────────────────────────────┐
 │ Paper 18 §IV.6 paste           │  (~1–2 days, parallel)
 └────────────────────────────────┘
 ┌─────────────────────────────────────────────────────────┐
 │ Phase C-W1c (independent, mechanical)                   │  (5–7 weeks, parallel)
 └─────────────────────────────────────────────────────────┘
 ┌──────────────────────────────────────────────────────────────────────────┐
 │ Phase C-W2b-easy (keystone, NCG theorem)                                 │  (8–14 weeks, blocking-sequential)
 └──────────────────────────────────────────────────────────────────────────┘
                                                              ┌────────────────────────────┐
                                                              │ Phase C-W1a-physics        │  (4–6 weeks streamlined,
                                                              │ (after W2b-easy)           │   sequential)
                                                              └────────────────────────────┘
```

The four small-paste / small-sprint deliverables (Track NI Zenodo, §VIII.D, §IV.6, W1c) can run in parallel with the keystone. W1a-physics is sequential after the keystone in the streamlined version; standalone fallback (6–10 weeks) is available if the keystone stalls.

**Total Phase C duration:** ~14–20 weeks if executed cleanly with parallelization. The keystone (W2b-easy) plus its application (W1a-physics) is the critical path at 12–20 weeks combined.

---

## Section 4: Pattern across all returns

**[Established by Phase B.]** Reading the seven memos as a panel, the structural pattern is sharper than Phase A's "the framework couples discrete labels but not space" framing.

### What the framework genuinely cannot do

- **W2a** (autonomous multi-loop counterterm generation): no published spectral-action framework solves this. Marcolli–vS 2014 constrains the ring; Hekkelman–McDonald 2024 extracts the leading Tauberian residue; Connes–Chamseddine 1996 treats renormalization as imposed externally. Closing W2a would require an extension of the spectral-action framework to multi-cutoff matching — a multi-paper / multi-year direction at the program's frontier, not a sprint target.
- **W2b-medium / W2b-hard** (cross-manifold spectral-triple composition with Hardy-sector or dimension mismatch): structurally blocked by the Paper 24 §V Coulomb/HO category mismatch. The Bargmann transform and Kustaanheimo–Stiefel do not bridge spectral triples. Closing this would require either a generalized propinquity framework for Hardy-sector × Riemannian factors, or a unifying ambient-manifold construction. Neither exists in the published literature.
- **W3** (autonomous inner-factor parameter selection): the framework admits Yukawa, magnetization, generation count as inputs but does not select them. Two specific second-axiom candidates (Phase 4H Tracks SM-B, SM-C) were tested and falsified. Fourteen catalogued speculations contain zero concrete proposals. The Bertrand × Hopf-tower upper bound has no matching lower-bound forcing argument. This may be a permanent frontier of the spectral-action program.

### What the framework can do but had not yet built

- **W1a / W1b / W1c**: all three are tooling-addressable. The mathematics is closed-form (Roothaan 1951 + incomplete-gamma + erfc identities for W1a/b; multi-shell exponential extension of `_split_integral_analytical` for W1c). The architecture supports the operators (Track NI cross-register tensor product; Sprint H1 inner-fluctuation classification). The engineering work is non-trivial but bounded (4–10 weeks per sprint), and each closure has a calibrated atomic-physics target (Pachucki et al. 2023; Eides 2024; Sprint 7 NaH/MgH$_2$ overattraction).
- **W2b-easy**: structurally a 4–8 week extension of Paper 38's machinery. Each of Paper 38's five lemmas has a clean factor-by-factor extension. The closure is a publishable NCG theorem at the program's frontier, not a routine application.

### The honest synthesis

The "structural-skeleton scope" framing is **partially correct and partially over-stated**. The framework genuinely:
- Maps selection rules cleanly (Paper 14, Paper 22).
- Classifies transcendentals via the master Mellin engine (Paper 32 §VIII; Paper 18 §III.7).
- Forces the upper bound on gauge content via Bertrand × Hopf-tower (May 2026 memo).
- Cannot autonomously generate counterterms, cannot select Yukawas, cannot bridge Coulomb and HO at the spectral-triple level.

But Phase B reveals the framework also:
- Has all the algebraic infrastructure for cross-register coordinate operators (W1a) and magnetization densities (W1b); it had simply not been built.
- Has closed-form integrals for cross-center frozen-core screening (W1c); it had simply not been wired through `balanced_coupled.py`.
- Has Paper 38 as the single-factor case of a tensor-product theorem (W2b-easy) that is reachable in 4–8 weeks of focused extension.

**The "multi-focal-composition wall" of CLAUDE.md §2 was correct empirically but the framing under-stated the framework's capacity.** Five of the six refined walls (W1a, W1b, W1c, W2a-as-Tauberian-side, W2b-easy) reduce, port, or are tooling-addressable. Only W2a's renormalization side, W2b-medium/hard, and W3 are genuine structural endings. This sharpens the structural-skeleton scope statement: the framework's structural ending is at *autonomous calibration data generation* and *multi-cutoff renormalization*, not at *multi-focal composition* per se.

---

## Section 5: WH register implications

**[Recommendation, PI decision territory per CLAUDE.md §1.7 + §13.5 access control.]** B-W3-diag §5 argued that **WH4 (the four-way S³ coincidence)** is structurally a one-way Fock-projection statement plus three forced consequences:

1. (Coulomb $\to S^3$) is forced by Bertrand's theorem + SO(4) symmetry of the Coulomb Hamiltonian (Paper 7).
2. ($S^3 = \mathrm{SU}(2)$ parallelizability $\to$ Dirac spinor bundle): once $S^3$ is established, Camporesi–Higuchi spinor bundle is the canonical object.
3. ($S^3 \to S^2$ Hopf base): the maximal-torus reduction of $\mathrm{SU}(2)$ on itself is automatic.
4. (Wilson SU(2) on $S^3 = \mathrm{SU}(2)$): strict-natural per Paper 30.

Under this reading, the four-way coincidence is **a structural unfolding**, not a "deep coincidence" that motivates the spectral-triple framing more than any single computation. The deflation argument is technically clean (Sprint TS-D's master-Mellin-engine reading + the Bertrand × Hopf-tower memo + Paper 30's strict-naturalness combine to make it a forced consequence) and consistent with WH1 PROVEN's actual content.

### Recommended language for §1.7 update if PI accepts

The current WH4 entry reads:
> **WH4 — The four-way S³ coincidence is one structure expressing itself four times.** S³ appears as (a) the Fock projection image of the hydrogenic spectrum, (b) the base of the Hopf bundle carrying Paper 2's α construction, (c) the spin-carrier for the Dirac sector, (d) the manifold of SU(2) gauge. ...

A deflated version could read:
> **WH4 (deflated, 2026-05-07) — The four-way $S^3$ structural unity is a Fock-projection statement plus three forced consequences.** Once Bertrand's theorem (closed-orbit selection of $-Z/r$) plus SO(4) symmetry forces the outer manifold to be $S^3$, three consequences follow without independent input: (b) $S^3 \to S^2$ Hopf base is automatic from $\mathrm{SU}(2)$'s maximal-torus reduction; (c) Camporesi–Higuchi spinor bundle is canonical via parallelizability; (d) Wilson SU(2) on $S^3 = \mathrm{SU}(2)$ is strict-natural. The "coincidence" reading was the chronological order of discovery, not a structural fact — the four roles are one role plus three consequences.
> *Status:* deflated to a **single-input forcing statement** plus three forced consequences. Does NOT extend to inner-factor selection (Yukawa, magnetization, generation count remain unforced). The *outer* spectral-triple's structural unity is essentially closed by Sprint TS-D + Bertrand × Hopf-tower + master Mellin engine; the *inner* factor's data remains W3.

### Why this matters

The deflation does not weaken the framework's positioning — Bertrand's forcing is the strongest single structural argument in the project record, and naming WH4 as its consequence rather than as an independent "deep" observation is **more honest, not less**. It also clarifies what the framework's structural unity does and does not do: it forces the outer triple but is silent on the inner factor (W3 / Yukawa / generation count).

### Three handling options

- **Accept the deflation** and update §1.7 with the language above. Clean and immediate.
- **Reject the deflation** if the PI views the four-way coincidence as load-bearing in some way Phase B did not capture. The original framing remains.
- **Scope a verification sub-sprint** to test the deflation's three forced consequences explicitly. (Recommended only if there is reasonable doubt; B-W3-diag's argument is structurally sharp enough that a verification sprint is unlikely to surface new content.)

**Recommendation:** accept the deflation. The argument is technically sound, the rewording is honest, and it tightens the WH register without weakening any active research direction.

---

## Section 6: Open questions for PI sign-off

The Phase B output produces the following decision points before any Phase C dispatch:

1. **Confirm strategic restructure (keystone-first vs atomic-physics-first)?** B-W2b-diag's §7 makes the keystone-first case; B-W1a-diag's §10 supports it via the algebraic-backbone-closed finding. Recommendation: keystone-first. PI: yes / no / modify.

2. **Approve Phase C-W2b-easy dispatch as the next major sprint?** 8–14 weeks, NCG keystone, 4–8 weeks symbolic + 2–4 weeks numerical + 2 weeks paper draft. Roadmap is the tensor extension of Paper 38's five lemmas. Risk profile: low-medium (joint L3 / L5 unpublished but reachable). PI: yes / no / scope-modify.

3. **Approve Track NI Zenodo deposit (option a) as a small parallel sprint?** ~1 week memo + ~1 day paragraph additions to Papers 23 and 32. Title: "The composed nuclear–electronic deuterium Hamiltonian as an explicit Connes-style real-space multi-particle spectral triple." Skeleton in B-position T5. PI: yes / no / change scope.

4. **Approve §VIII.D LaTeX paste into Paper 32?** ~970 words drafted in B-position T3, paste-ready. Names W1a/b/c as GeoVac-internal architectural gaps and W2a/W2b as broader spectral-action open problems. Adds three bibitems (`hekkelman_mcdonald2024`, `latremoliere2026`, `paper36`). Flags existing `hekkelman2024` placeholder for cleanup. PI: paste as-is / edit / defer.

5. **Approve Paper 18 §IV.6 W3 framework's-open-question addition?** ~150 words drafted in B-W3-diag §8. Names the second-packing-axiom question as open without committing to closure. PI: paste as-is / edit / defer.

6. **WH4 deflation: accept, reject, or scope a verification sub-sprint?** Recommendation: accept. Updated §1.7 language drafted in Section 5 above. PI: accept / reject / verify.

7. **Bertrand × Hopf-tower lower-bound gap: scope a Phase C-W3 attempt or leave as open question?** B-W3-diag §6 confirmed no structural argument fixes the lower bound; the closest target ("any GeoVac packing produces only spheres of dim 3 and 5") is a desired theorem with no proof attempt. Recommendation: leave as open question and record in Paper 18 §IV.6 (option 5 above). A Phase C attempt would be high-risk and exploratory. PI: scope / leave open.

**Two minor questions** flagged by B-W2a-diag and B-position that do not require PI sign-off but should be visible:
- (a) The existing Paper 32 `hekkelman2024` bibitem at line 2693–2697 has a placeholder arXiv number (`arXiv:2403.xxxxx`). B-position T4 flagged this for cleanup. Recommended action: PM repairs this when next touching Paper 32's bibliography (no separate PI sign-off needed).
- (b) Numerical compatibility check between HM 2024 and the master Mellin engine on a panel of 5–10 Dirac-$S^3$ observables (~1–2 weeks). Optional follow-up to B-W2a-diag; would tighten the Paper 32 §VIII.D framing edit. Recommendation: defer unless §VIII.D is later challenged on rigor grounds.

---

## Section 7: Phase C structure recommendation

**[Recommendation.]** Given Sections 1–6, the cleanest Phase C structure is **one large sprint (W2b-easy) as the headline, three small parallel positioning artifacts, and one independent mechanical sprint (W1c)**, with W1a-physics as the streamlined sequential follow-up.

### The single most important question Phase C must answer

**Does Phase 38's five-lemma machinery extend to the tensor product of two Camporesi–Higuchi spectral triples on $S^3$ at distinct focal lengths?**

This question is load-bearing for three reasons:

1. **It is genuinely open in the published NCG literature** (Track 3 surprise S2; B-W2b-diag §2 verified). Its closure would be a publishable NCG-frontier theorem parallel to Paper 38's WH1 closure.
2. **It is the algebraic foundation for W1a-physics.** Closing it converts W1a-physics from "build a new infrastructure from atomic-physics-side" to "apply a published NCG theorem to atomic spectroscopy" — a much more efficient sprint.
3. **It tests the framework's reach.** Paper 38 closed the single-$\mathrm{SU}(2)$-factor case. Closing the tensor-product case demonstrates that the GeoVac spectral-triple framework extends to multi-particle configurations at the GH-convergence level — a structural strengthening of WH1 PROVEN.

### Sprint ordering

The blocking-sequential order is:

1. **Week 0:** Three small positioning artifacts dispatched in parallel (~1–2 days each):
   - Track NI Zenodo deposit memo + paragraph additions to Papers 23 and 32.
   - Paper 32 §VIII.D paste (B-position T3 LaTeX).
   - Paper 18 §IV.6 W3 framework's-open-question paste (B-W3-diag §8 LaTeX).
   - WH4 §1.7 deflation update (Section 5 above) if PI accepts.
   - All four are independent; can land in a single PM session.

2. **Weeks 1–8 (parallel tracks):**
   - **Phase C-W2b-easy** (keystone, weeks 1–8 symbolic): tensor extension of Paper 38's five lemmas. L1' / L2 / L3 / L4 lift factor-by-factor; joint L3 ($C_3 \le 2$ general, $C_3 = 1$ on factorized panel) and joint L5 (Latrémolière 2017/2023 propinquity tunneling-pair height bookkeeping) are the two technical risks.
   - **Phase C-W1c** (independent mechanical, weeks 1–7): module `cross_center_screened_vne.py`; wire into `balanced_coupled.py`; NaH/MgH$_2$/HCl/etc. PES regression. Parallel and independent.

3. **Weeks 8–12:** Phase C-W2b-easy numerical verification (`(n_a, n_b) \in \{2,3,4\}^2$ panel via existing `gh_convergence.py` and `connes_distance.py`) + paper draft.

4. **Weeks 12–18:** **Phase C-W1a-physics** (streamlined sequential after W2b-easy): module `cross_register_vne.py` (HO + Sturmian); promote `R_PROTON_BOHR` to operator-valued; magnetization-density component for W1b; Pachucki–Patkóš–Yerokhin 2023 + Eides 2024 calibration. 4–6 weeks streamlined.

### Expected Phase C duration

**14–20 weeks total** for the full plan if W2b-easy executes cleanly and W1a-physics streamlined applies. Phase C produces:

- **One NCG-frontier theorem** (W2b-easy; "Paper 38b" or §VI of Paper 38 v2).
- **One atomic-physics calibration** (W1a-physics; Lamb shift recoil + 21 cm hyperfine validation against Pachucki + Eides).
- **One tooling fix** (W1c; second-row PES capability restored for NaH/MgH$_2$/HCl/H$_2$S/PH$_3$/SiH$_4$ + 8 third-row species).
- **One Zenodo memo** (Track NI architectural observation).
- **Two paper edits** (Paper 32 §VIII.D + Paper 18 §IV.6).
- **One WH register update** (WH4 deflation, if PI accepts).

If W2b-easy stalls on joint L3 or L5: rescope to standalone W1a-physics (6–10 weeks atomic-physics-first), keep W1c and the positioning artifacts intact. **Total Phase C remains 12–20 weeks; the headline reverts from NCG-frontier to atomic-physics calibration but the deliverable count is preserved.**

### Phase C exit

Phase C exits when:
- W2b-easy theorem is proved or cleanly stalled (with diagnostic memo).
- W1a-physics sprint has produced cross-register operator + Pachucki/Eides validation, or has been rescoped to standalone form.
- W1c has restored second-row PES capability or has documented residual structural error.
- Track NI Zenodo memo is deposited.
- §VIII.D and §IV.6 are pasted.

The Phase C synthesis memo at exit will close out the multi-focal-composition sprint and either (a) declare the wall taxonomy resolved at the level Phase B's diagnostics permit, or (b) identify the remaining open structural questions with verdict-conditional Phase D recommendations. CLAUDE.md §2 multi-focal-wall entry should be revised at Phase C exit to reflect the post-keystone state.

---

**End of Phase B synthesis memo.**
