# Sprint P5 — Hybrid-Pipeline Paper-Home Scoping Memo

**Date:** 2026-06-07.
**Sprint position:** P5 (paper-home decision) within the hybrid-pipeline scoping arc launched in `debug/sprint_hybrid_pipeline_scoping_memo.md` (2026-06-07). Phase 1 is 0–2 months; Phase 3 (publication) lands at month 7–8. This sprint resolves the paper-home question NOW so that Phase 1 writes against a clear target.
**Verdict line:** **(b) — write a new Paper 57 dedicated to the hybrid-pipeline chemistry result.** JCTC as primary venue; J. Chem. Phys. as fallback; Quantum / npj QI for an optional VQE-focused spinoff. Reasoning below.
**Cross-references:** `debug/sprint_hybrid_pipeline_scoping_memo.md` (parent scoping); `papers/group4_quantum_computing/paper_14_qubit_encoding.tex` (3,039 lines, qubit-encoding theory paper, large); `papers/group4_quantum_computing/paper_20_resource_benchmarks.tex` (841 lines, chemistry-consumer resource benchmarks); CLAUDE.md §6 group folder structure; CLAUDE.md §1.5 (positioning, chemistry-consumer rhetoric); `memory/external_input_three_class_partition.md` (Class 1/2/3 partition).

---

## §0. TL;DR

Parent scoping §5 item 4 defaulted to **(a) Paper 20 extension**. This sprint **revises that default to (b) new Paper 57** after sizing Papers 14 (3,039 lines, ~50 pp) and 20 (841 lines, ~22 pp) and reading their headlines.

Driving observation: the hybrid-pipeline story has a **different load-bearing claim** than either. Paper 14: "GeoVac integrals are structurally sparser than Gaussian." Paper 20: "GeoVac integrals reach reference accuracy at low resource cost." Hybrid: "GeoVac integrals + external correlation engine = production-grade quantum chemistry that respects the structural-skeleton-scope boundary (Paper 18 §IV.6)." A *third* claim with its own audience (computational chemists already using DMRG/CCSD(T)) and its own falsifiability test (NaH W1e closure under DMRG correlation = Class 3 empirically validated).

A dedicated Paper 57 gives the cosmic-Galois Class 1/2/3 partition its first empirical operationalization in a single venue — hard to do as a §V subsection inside a paper with a different headline.

Decision-gate verdict: **SHIP-AS-IS**, with one PI-input flag in §7 (optional VQE-focused Paper 58 if Phase 2 H2-VQE lands clean).

---

## §1. Candidate analysis

### §1.1 Candidate (a): expand Paper 20 §benchmarks

**Current audience and venue.** Resource benchmarks for the QC-applied-to-chemistry community. 841 lines (~22 pp). Already includes "Installation and usage," "Available Systems" inventory, accuracy benchmarks (LiH/BeH₂/H₂O), and Pauli/QWC/1-norm comparisons vs Gaussian. Implicitly JCTC-shaped.

**Audience fit.** Good — overlapping community, but a non-trivial narrative shift: Paper 20's headline is sparsity-and-accuracy; the hybrid headline is sparsity-plus-external-correlation-equals-production-grade.

**Sacrificed.** Adding a 28-molecule DMRG sweep + VQE hardware-aware benchmark + W1e classification adds 8–12 pages (→ 30–34 total). Limitations + Conclusion need rewriting (W1e flips from "framework limitation" to "external-tier consumed by hybrid"). Narrative tightness is lost.

**Gained.** Zero new-paper overhead; cross-ref network already correct (Papers 14/17/19/22 cited); ~6 wk Phase-3 writing vs ~2 mo for new paper; existing arxiv mirror.

**Verdict (a).** Workable but uncomfortable. Paper 20 becomes a 30+ pp hybrid with two competing headlines.

### §1.2 Candidate (b): new Paper 57 dedicated to hybrid chemistry

**Audience.** Computational chemists who already use DMRG/CCSD(T)/AFQMC/VQE — **categorically different** from Paper 14's theory readers and **broader** than Paper 20's resource-benchmark readers. The pitch "we provide integrals, you bring your correlation engine" is the natural JCTC/JCP framing.

**Headline claim.** "Structurally sparse GeoVac Hamiltonians, consumed by DMRG/CCSD(T)/VQE, reach production-grade quantum chemistry while preserving a structurally-derived skeleton/calibration boundary (Class 1 / Class 3) that maps each method onto its native tier." Own arc; not a subsection.

**Audience-shift requirement.** None. Writes directly at JCTC audience; CLAUDE.md §1.5 chemistry-consumer rhetoric applies cleanly.

**Sacrificed.** ~2 mo extra writing; new arxiv ID; new bib (~50 entries, ~30 overlap with Paper 20); +1 paper in group 4.

**Gained.** (i) clean headline; (ii) clean audience; (iii) NaH W1e closure foregrounded as load-bearing structural test, not buried; (iv) first chemistry-side empirical operationalization of Paper 18 §IV.6 in a dedicated venue (matters for parent scoping §5 item 6); (v) "see Paper 57" cites cleanly; (vi) Paper 20 stays crisp.

**Verdict (b).** Clean fit. Recommendation.

### §1.3 Candidate (c): expand Paper 14 with hybrid-application section

**Current audience.** NCG / QC-theory community. Headline: "qubit Hamiltonians from spectral graph theory are structurally sparse." Already 3,039 lines (~50 pp), covering Pauli/QWC/1-norm scaling, composed/second-row/multi-center/TM-hydride/nested/spinor extensions, DF comparison, and long Discussion + Limitations.

**Audience fit.** Poor. Paper 14 readers want algebraic / spectral-theoretic arguments; hybrid-pipeline readers want benchmark tables and FCIDUMP docs. Hard audience mismatch.

**Sacrificed.** Adding ~10–15 pp pushes Paper 14 past 60 pp. Narrative splits across "theory of structural sparsity" + "applied chemistry pipeline." Paper 14 currently frames N-electron correlation as out-of-scope; the hybrid section would reverse that. Major rewrites of Discussion/Limitations/Conclusion required.

**Gained.** Nothing that (a) or (b) don't get more cleanly.

**Verdict (c).** Reject. Audience mismatch + length blowout.

---

## §2. Recommended choice with reasoning

**Recommendation: (b) — new Paper 57 dedicated to hybrid chemistry.**

Decisive factor: **headline claim distinctness**. The hybrid result is neither a benchmark refinement of Paper 20 nor an application extension of Paper 14. It is a *new structural claim*: the Class 1 / Class 3 partition (Paper 18 §IV.6) is the correct architectural reading, empirically validated by NaH W1e closure under DMRG. Three reasons it deserves its own home:

1. **Falsifier visibility.** Sprint H1-NaH (parent scoping §4.1) IS the Paper 18 §IV.6 falsifier. Buried in Paper 20 §V it reads as a benchmark line; foregrounded as Paper 57's main result it reads as a structural test.
2. **Audience.** Paper 20 talks to chemists who want a black-box Hamiltonian; Paper 57 talks to chemists who want to integrate GeoVac into their pipeline. Different conversations; JCTC routes them to different referee pools.
3. **Corpus hygiene.** Paper 20 stays the crisp resource-benchmarks reference; Paper 57 stays the crisp hybrid-pipeline reference. Forcing both into Paper 20 produces a 30+ pp paper with two competing headlines.

Cost (~2 mo extra writing + 1 more paper) is small vs the ~7–8 mo total roadmap and is paid only at Phase 3.

**Concession to (a):** if Phase 2 produces fewer headline-grade results than projected (e.g. DMRG closes NaH cleanly but the sweep is uninteresting), fold what remains into Paper 20 §V.E. Phase-3-entry reassessment, not a now-decision.

---

## §3. Journal/venue targeting

Parent scoping memo defaulted to **JCTC**. Validating below.

| Venue | Audience | Hybrid-pipeline fit | Verdict |
|:------|:---------|:--------------------|:--------|
| **JCTC** (*J. Chem. Theory and Comput.*) | Comp chemists; method-development focus; high-volume DMRG/CCSD(T)/AFQMC/VQE methods papers | **Best fit.** JCTC routinely publishes "method X coupled to method Y for system Z" papers. Hybrid-pipeline architecture maps directly to JCTC's editorial taste. Standard 12–20 page format. Strong open-access pathway via ACS AuthorChoice. | **Primary venue.** |
| **JCP** (*J. Chem. Phys.*) | Comp chemists + AMO + spectroscopy; broader audience than JCTC | Workable. JCP would take this but the venue rewards methodology-development with quantitative benchmark contributions more than structural-framework arguments. JCP referees may push back on the Class 1/2/3 partition language as too NCG-flavored. | **Primary fallback.** |
| **Quantum** | Quantum-computing community; theoretical + applied | Good fit IF the VQE result is the headline. Quantum's audience reads Pauli scaling + hardware-aware results as the main content. If Phase 2 H2-VQE lands clean (sub-percent tapered VQE on LiH at N_qubits=27), a Quantum spinoff makes sense. | **Spinoff venue (conditional on Phase 2 H2-VQE).** |
| **npj Quantum Inf.** | Quantum-computing community; Nature group venue, broader reach | Same as Quantum but with higher visibility / higher rejection rate. Spinoff option. | **Alternative spinoff.** |
| **PRX Quantum** | Quantum-computing community; rigor + breadth | Possible if VQE result has structural-novelty hook (e.g. per-block Hopf-U(1) tapering as a new sparsity primitive). Hard sell for the chemistry-pipeline content alone. | **Less likely; defer.** |
| **PR Applied** | Applied physics; less chemistry-aware | Wrong audience. | **Reject.** |

**Net targeting recommendation:**
- **Paper 57 → JCTC** as primary submission. Fallback to JCP if JCTC rejects (~25% rejection rate at JCTC for methods papers in this class).
- **Optional Paper 58 (VQE-focused) → Quantum** if Phase 2 H2-VQE produces a publishable headline. This is a Phase-3-entry decision, not a now-decision.
- Both manuscripts mirrored to arxiv (`quant-ph` for both) per CLAUDE.md §1.

---

## §4. Phase 3 paper outline (Paper 57)

**Title (working):** *Structurally Sparse Quantum-Chemistry Pipeline: GeoVac Integrals Consumed by DMRG and VQE for Production-Grade Molecular Simulation*

**Page budget:** 18–22 pages (REVTeX twocolumn, JCTC format). Bibliography ~50 entries (overlap with Paper 20 ~30, new chemistry-pipeline entries ~20).

**Headline result statement:** "GeoVac-derived molecular integrals, consumed by external correlation engines (DMRG, CCSD(T), VQE), reach sub-percent reference accuracy across 28 molecules in the GeoVac library while preserving O(Q^2.5) Pauli scaling and per-block Hopf-U(1) qubit tapering. The architecture cleanly separates the structural skeleton (sparse integrals + selection rules + qubit encoding, computed by GeoVac) from the multi-determinant correlation tier (Class 3 per Paper 18 §IV.6, handled by the external engine). The NaH wall closure under DMRG correlation validates the Class 1 / Class 3 cosmic-Galois partition at production scale."

**Section structure:**

1. **Introduction** (~2 pp). Pipeline-architecture problem. Sparse-integral approaches vs Gaussian. Where multi-determinant correlation lives. Brief Paper 18 §IV.6 motivation. Headline numbers up front.
2. **GeoVac integrals: brief recap** (~2 pp). Pointer to Papers 14/17/22. Two figures (Pauli scaling; O(Q^2.5) across 28 molecules).
3. **Interface design** (~3 pp). `to_fcidump`, the 6 convention-pinning items, RHF cross-check methodology, per-block Hopf-U(1) tapering for non-VQE consumers.
4. **DMRG benchmarks** (~4 pp). Block2/pyscf-DMRG on LiH at increasing bond dimension; convergence to GeoVac-FCI; cc-pVTZ comparison. NaH W1e closure (load-bearing). 28-molecule sweep table.
5. **VQE benchmarks** (~3 pp). Tapered VQE on LiH at N_qubits=27; UCCSD on Gaunt-allowed singles/doubles; noise-aware shot budget; OpenFermion/Qiskit/PennyLane platform-agnostic.
6. **Class 1 / Class 3 partition empirically** (~3 pp). Per-molecule decomposition into Class-1 calibration mismatch, Class-3 multi-determinant correlation handled by DMRG, residual. NaH + TM-hydride case studies.
7. **Discussion + limitations** (~2 pp). Position within DF/THC/compressed-DF landscape. Honest scope: pipeline closes Class 3, NOT Class 1 (calibration external) or Class 2 (composition, framework-internal).
8. **Conclusion** (~1 pp).

**Figures (~6) + tables (~4).** GeoVac vs Gaussian Pauli scaling; architecture schematic; LiH bond-dim convergence; NaH PES (GeoVac-FCI vs DMRG vs cc-pVTZ); tapered-VQE LiH under noise; Class 1/3 decomposition. Tables: convention-pinning checklist, 28-molecule benchmark, VQE shot-budget, DF baseline comparison.

**Word budget:** ~18–22 pp × 350 w/2col-page ≈ 6,500–8,000 w main + ~2,000 w appendices ≈ 10,000 w total.

---

## §5. Cross-references to GeoVac corpus

**Papers Paper 57 cites:**

| Cited paper | Role in Paper 57 |
|:-----------|:------------------|
| Paper 7 | Foundational S³ Fock projection (one-line cite in §2). |
| Paper 14 | O(Q^2.5) Pauli scaling; per-block Hopf-U(1) tapering; ERI density decay. Heavy cite in §2 (recap section). |
| Paper 17 | Composed-geometry architecture (LiH, BeH₂, H₂O). Cite in §2 and §6. |
| Paper 18 | §IV.6 chemistry-side analog and the Class 1 / Class 3 partition. Load-bearing cite in §6. §III.7 master Mellin engine (one-line cite at structural-skeleton scope). |
| Paper 19 | Balanced-coupled architecture (cross-center V_ne). Cite in §2. |
| Paper 20 | Resource benchmarks + Available Systems inventory. Cite in §1 (motivation) and §4 (benchmark reference). |
| Paper 22 | Universal angular sparsity theorem; ERI density independent of V(r). Cite in §2 and §6. |
| Paper 23 | Composed nuclear-electronic architecture (one-line cite in §6 for the multi-particle skeleton extension as forward-looking material). |
| Paper 31 | Universal/Coulomb partition; A/D split. Cite in §1 (motivation: framework-side reading of why the hybrid pipeline works). |

**Papers that cite Paper 57 back:**

| Paper | Citation context |
|:------|:-----------------|
| Paper 14 (next revision) | "Validated empirically by Paper 57's hybrid pipeline." Single-line forward-ref in Discussion. |
| Paper 18 (next revision) | "§IV.6 chemistry-side analog operationalized empirically in Paper 57." This is the load-bearing back-cite — Paper 57 becomes Paper 18's chemistry-side empirical anchor. |
| Paper 20 (next revision) | "See Paper 57 for the hybrid-pipeline integration of these benchmarks with external correlation engines." Single-line forward-ref in Limitations. |
| Future Paper 58 (VQE spinoff, conditional) | Paper 58 would cite Paper 57 as the parent hybrid-pipeline paper if it's split off. |

**Net corpus impact:** +1 paper in group 4 (now 5 papers + Track NI memo); ~7 new inbound cite edges; cross-group structural reading: Paper 57 is the chemistry-empirical operationalization of Paper 18 §IV.6 (cross-group 3↔4 anchor; presently no such anchor exists at empirical level).

---

## §6. Risk register

| Risk | Severity | Mitigation |
|:-----|:---------|:-----------|
| **JCTC referees push back on Class 1/2/3 partition language as too NCG-flavored.** | Medium. JCTC audience is computational-chemistry-applied; the cosmic-Galois reading may read as gratuitous theory to them. | Section 6 lead-paragraph reframes Class 1/2/3 in pure-chemistry vocabulary first (calibration / composition / multi-determinant correlation) before citing Paper 18 §IV.6 in a footnote. Paper-internal vocabulary is chemistry-native; corpus cross-ref is the technical citation. |
| **Scope creep: Paper 57 becomes 30+ pages and overlaps Paper 20.** | Medium. Three benchmarks (DMRG, CCSD(T), VQE) × 28 molecules × multiple references → table sprawl easy. | Phase 3 writing sprint enforces 20-page hard cap. Move table appendices to supplementary information (JCTC accepts SI freely). Only 3 figures + 2 tables in main text; rest to SI. |
| **VQE Phase 2 underperforms (shot noise + hardware error make sub-1% out of reach on LiH at N_qubits=27).** | Medium-low. NISQ-era VQE on first-row systems is well-characterized; sub-1% is achievable under noise-mitigated protocols but not guaranteed. | Paper 57 §5 hedges on VQE numbers (report shot budget, error bars, noise model); the load-bearing claim is DMRG-based. If VQE underperforms, §5 stays as a hardware-applicability demonstration not a hardware-accuracy claim. Spinoff Paper 58 only happens if VQE lands clean. |
| **NaH W1e wall does NOT close cleanly under DMRG (parent scoping §5 item 6 falsifier).** | Medium. The Class 1 / Class 3 partition predicts DMRG closes W1e because W1e is Class 3; if DMRG does NOT close it, the partition is partially falsified and Paper 57's structural reading needs to revise. | Per parent scoping §5 item 6, treat as finding-not-exit. Paper 57 §6 reports the actual W1e closure result; if it's partial, the paper documents the partial closure with the structural implication (W1e has Class 3 + something else). This is publishable either way — failed falsifier with structural implication is also strong content. |
| **Reviewers want comparison with Gaussian+DMRG (the existing baseline) rather than GeoVac-FCI.** | Low-medium. Pyscf + Block2 on cc-pVTZ + LiH is the standard external comparison. Phase 2 sweep already plans to run this. | Build the comparison into the Phase 2 sweep architecture; report both energies in the 28-molecule benchmark table. Phase 1 H1-Reference sprint covers cc-pVTZ-FCI on LiH directly. |
| **Project corpus growing toward 60 papers.** | Low (philosophical, not technical). | Paper 57 fills a real audience gap (chemistry-pipeline, no current paper covers this). Splinter-control discipline lives in CLAUDE.md §13.8; Paper 57 has its own headline and cross-reference network, so it's not a splinter. |
| **Paper 14 / Paper 20 narrative drift if the hybrid story changes them.** | Low. Both papers cited cleanly; single-line forward-refs added in their next revisions; no structural changes needed. | Restrict cross-paper edits to single-line forward-refs ("see Paper 57"). No major rewrites to Paper 14 or Paper 20. |

---

## §7. Decision-gate verdict

**SHIP-AS-IS.** Recommendation is clear: write a new Paper 57 dedicated to the hybrid-pipeline chemistry result, target JCTC primary / JCP fallback. Phase 3 writing sprint begins month 7–8 per parent scoping roadmap; this memo sets the writing target so Phase 1 + Phase 2 can be executed against a clear destination.

**One PI-input flag for later (not blocking Phase 1 launch):** the Phase-3-entry decision of whether to split off a VQE-focused Paper 58 (target: *Quantum*) depends on Phase 2 H2-VQE performance under realistic noise. If H2-VQE lands sub-1% with the tapered Hamiltonian, Paper 58 is worth ~6 weeks extra writing time for a high-visibility Quantum publication. If H2-VQE comes in at 1–5% (likely), keep VQE as §5 of Paper 57 and skip Paper 58. **PI input needed at Phase 3 entry (~month 6.5–7), not now.**

**Concession path:** if Phase 2 sweep produces fewer headline-grade results than projected (DMRG-on-GeoVac matches DMRG-on-Gaussian within a few percent and NaH wall closes but uninterestingly), fold the result into Paper 20 §V.E and skip Paper 57. This is a Phase-3-entry reassessment trigger, not a now-decision.

---

## §8. Honest scope

What this scoping demonstrated:

- The hybrid-pipeline result has a structurally distinct headline claim from both Paper 14 and Paper 20.
- Paper 14 is the wrong venue (audience mismatch + length blowout).
- Paper 20 expansion is workable but uncomfortable (length blowout + two competing headlines).
- New Paper 57 is the cleanest fit on audience, headline, and falsifier visibility.
- JCTC is the right primary venue; JCP is the fallback; Quantum is the conditional spinoff venue if Phase 2 VQE lands clean.

What this scoping did NOT demonstrate:

- **No code, no paper edits.** Pure scoping document per sprint mandate. Paper 57 outline is rough.
- **No bibliography construction.** ~50 bibitems projected; actual list built at Phase 3 entry.
- **No assessment of competing hybrid-pipeline efforts.** A literature scan for "DMRG + sparse-integral generators" / "VQE + custom-basis pipelines" is a Phase 2 prerequisite to scope reviewer-expectation alignment. Out of scope here.
- **No VQE-spinoff Paper 58 outline.** Paper 58 is conditional on Phase 2 results; no outline drafted.
- **No assessment of arxiv-listing strategy.** Cross-listing decisions (`quant-ph` primary; `physics.chem-ph` secondary?) deferred to Phase 3.

---

## §9. Files

### Created (this sprint)
- `debug/sprint_p5_paper_home_scoping_memo.md` (this memo).

### NOT modified
- Production `geovac/` modules — pure scoping sprint per sprint mandate.
- Papers — Paper 14, Paper 20, and Paper 18 NOT edited; cross-references to a future Paper 57 are not added in this sprint (they would be added at Phase 3 entry, only if Paper 57 is greenlit by PI sign-off on this memo).
- CLAUDE.md — no edits; scoping memo only. §6 paper inventory would gain Paper 57 only at Phase 3 entry.
- Tests — no code changes.

---

**End of Sprint P5 paper-home scoping memo. Verdict: (b) new Paper 57, JCTC primary, JCP fallback, optional Quantum spinoff at Phase 3 entry. Phase 1 launches against this writing target.**
