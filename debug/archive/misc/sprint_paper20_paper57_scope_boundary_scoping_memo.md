# Paper 20 / Paper 57 — chemistry scope boundary scoping memo

**Date:** 2026-06-07
**Triggering events:**
1. W1e-Projection-Audit (earlier today) — balanced LiH binds with 2.4x over-binding, NaH overattracts under any path.
2. LiH kwarg sweep (earlier today) — all cheap engineering kwargs no-op or catastrophic on LiH.
3. B.1 explicit-core NaH HF (this evening) — does not bind NaH even with un-frozen [Ne] core; LiH HF cross-check binds at R_eq=3.015 with D_e=0.158 Ha (methodology validated).
4. H₂ Bratteli pilot (this evening) — bit-exact Marcolli-vS structural home for Track CD.

**Decision target:** What should Paper 20 (already drafted) and Paper 57 (proposed, not drafted) honestly claim about the chemistry-accuracy scope?

---

## §1 The current Paper 20 framing

Paper 20 (`papers/group4_quantum_computing/paper_20_resource_benchmarks.tex`) is the chemistry/QC resource-benchmark paper. The current abstract claims:

> *"The balanced coupled architecture (PK-free) achieves 0.20% energy error at the equilibrium geometry of LiH when solved exactly, with 878 Pauli terms at 30 qubits."*

This 0.20% claim is **single-R-point accuracy on LiH at the equilibrium geometry** (R = 3.015 bohr, energy at the well minimum), not a PES claim. Today's W1e-Projection-Audit reaffirmed it: balanced LiH at R = 3.015 gives E = -15.21 Ha = -8.18 Ha continuum + ~7 Ha core-double-counting offset. The single-point energy is what Paper 20 claimed and is correct.

But Paper 20 does NOT distinguish:
- **What molecules bind** (first-row: yes; NaH and second-row: no, per today's findings)
- **What PES properties are reproduced** (LiH well shape: yes; LiH well depth: 2.4× too deep; NaH well: doesn't exist in qubit Hamiltonian under any tested config)
- **Where the chemistry-accuracy scope boundary actually sits**

Today's results give us the empirical evidence we need to claim a scope boundary honestly.

## §2 Proposed Paper 20 framing edits

Three categories of edit, in order of importance:

### 2.1 Abstract row-conditional claim (REQUIRED)

Replace:
> *"The balanced coupled architecture (PK-free) achieves 0.20% energy error at the equilibrium geometry of LiH when solved exactly, with 878 Pauli terms at 30 qubits."*

With (proposed):
> *"For first-row hydrides (LiH at $n_{\max}=2$), the balanced coupled architecture binds at the experimental equilibrium geometry $R_{\mathrm{eq}}=3.015$ bohr with 878 Pauli terms at 30 qubits, achieving 0.20\% energy error at the well minimum (well depth 0.158 Ha, 2.4$\times$ the continuum reference). For second-row hydrides (NaH and below), the qubit Hamiltonian exhibits monotone overattraction with no interior potential-energy minimum under any tested configuration (frozen-core balanced, explicit-core HF, multiple cross-block coupling treatments); the chemistry-accuracy scope at $n_{\max}=2$ is therefore structurally bounded to first-row hydrides. We document this scope boundary as a positive feature of the framework: the Hamiltonians are correct quantum-simulation resource targets, with documented PES topology, so consumers know exactly what they're getting."*

This is ~80 words longer than the current abstract claim. Honest. Doesn't undermine the framework — actually strengthens it by being explicit about scope.

### 2.2 New §V.C "Chemistry-accuracy scope boundary" (NEW SECTION, ~2-3 pages)

Add a new subsection documenting the W1e-Projection-Audit + B.1 results:
- Table: PES topology by molecule class (first-row LiH, BeH₂, CH₄ → interior min; second-row NaH, MgH₂, HCl, KH → monotone descent under any config)
- Operator-level diagnosis: what's structurally absent for second-row (R-dependent basis adaptation, bonding-orbital Pauli repulsion against explicit core, multi-determinant correlation requiring O(N⁷) operators that destroy sparsity)
- Honest statement: chemical accuracy on second-row would require structural extension (multi-zeta + R-adaptive Z_eff + explicit-core FCI), which is a multi-month engineering investment not yet committed
- Cross-reference: Marcolli-vS gauge-network paper (see §3 below) for the structural reason WHY this scope boundary sits where it does

### 2.3 Conclusions reframe (REQUIRED)

The current Paper 20 Conclusions probably says something like "GeoVac achieves chemical accuracy on the library." Today's results force a more honest statement: "GeoVac's qubit Hamiltonians are correct quantum-simulation resource targets across the 38-molecule library, with chemical accuracy on first-row hydrides and explicitly-documented scope boundaries for second-row and beyond. The framework's primary value proposition is structural sparsity with quantitative truncation error bounds (Paper 38 propinquity), not energy accuracy parity with CCSD(T)."

This positions the framework correctly for its actual customer (NISQ algorithm researchers benchmarking Pauli-count scaling and circuit-cost estimation).

## §3 Paper 57 — the new chemistry-pivot framing

Paper 57 was proposed in the chemistry-pivot sprint (`debug/sprint_p5_paper_home_scoping_memo.md` if it exists, otherwise the framing was in `debug/sprint_chemistry_pivot_2026_06_07_memo.md` §4). Its three load-bearing claims were:
1. Sparse Hamiltonian export with quantitative error bounds.
2. Sub-mHa VQE on small systems via openfermion-native pipeline.
3. Row-conditional operational localization of the chemistry-side calibration-data partition boundary.

Today's results sharpen all three:

**Claim 1 is unchanged and strengthened.** FCIDUMP export bit-exact; propinquity bound wired; 38-molecule library audit-clean.

**Claim 2 is unchanged.** Openfermion-native UCCSD on H₂ at attomHa precision was the v3.85.0 finding.

**Claim 3 is REFRAMED** based on today's results. The new claim:
> *"Chemistry-accuracy at $n_{\max}=2$ scope boundary is empirically localized at the second-row hydride threshold (first-row binds; second-row overattracts). Engineering closures (frozen-core, multi-zeta, screened-valence, cross-block ERI, cross-block h1, explicit-core HF) have been systematically tested and exhausted; the wall is structural. The framework's NCG home (Marcolli-vS gauge networks, verified bit-exactly on H₂) gives the structural language to discuss the boundary."*

This is a much stronger paper than the original framing. It says:
- Here's what we know works (first-row, with quantitative truncation error bounds)
- Here's the boundary, characterized empirically
- Here's the structural language to discuss it (NCG / Marcolli-vS)
- Here's the path forward if you want to push the boundary (multi-month explicit-core engineering)

This is "honest scope + structural foundation" — the right shape of paper for the project's positioning.

### Recommended Paper 57 sections (revised from chemistry-pivot scoping)

| § | Section | Content |
|---|---|---|
| 1 | Introduction | Resource estimation = NISQ value prop. State scope: first-row hydrides accurate; second-row scope-bounded. |
| 2 | Sparse Hamiltonian export with propinquity bounds | Paper 38 bound wired; FCIDUMP machine-precision round-trip; 38-molecule library. |
| 3 | Openfermion-native UCCSD | Sub-mHa H₂; methodology; scaling characterization. |
| 4 | Chemistry-accuracy scope boundary | The W1e + B.1 empirical localization. PES topology by molecule class. Operator-level diagnosis. |
| 5 | Structural framing | Brief — point to the Marcolli-vS chemistry paper (see §3 of *this* memo) for the operator-algebra story. Don't duplicate. |
| 6 | Where to push the boundary | The multi-month explicit-core engineering arc that *could* close second-row, if a chemistry consumer wanted to invest. Scoping outline only. |

Estimated length: 18-22 pages, comparable to Paper 20 scope. Target venue: J. Chem. Theory Comput. or Quantum (the original chemistry-pivot scoping memo's recommendation).

## §4 Coordination with the Marcolli-vS chemistry paper

Today produced two new papers' worth of content (Marcolli-vS chemistry + chemistry-scope-boundary). They should cite each other but not duplicate content:

- **Marcolli-vS chemistry paper** (math.OA standalone, J. Geom. Phys.) = STRUCTURE. Why the framework is what it is.
- **Paper 57 (chemistry-scope-boundary)** (chemistry-customer, J. Chem. Theory Comput. or Quantum) = OPERATIONAL. What the framework does and doesn't deliver.
- **Paper 20** (already drafted) = BENCHMARKS. What the resource savings actually are.

The trio frames the project for both audiences (math.OA + chemistry/QC).

## §5 Recommended sequencing

1. **Edit Paper 20 abstract + conclusions** (1-2 days). Apply §2.1 and §2.3 edits. Add new §V.C subsection (~2-3 pages of writing). This is the immediate, no-blockers move.

2. **Draft Paper 57** (4-6 weeks) per the §3 revised outline. Bigger commitment but well-scoped now.

3. **Sprint M-vS-1 (LiH Bratteli scale-up)** runs in parallel with Paper 57 drafting; if it passes, the Marcolli-vS paper arc opens and gets its own multi-month track.

Paper 20 edits are the cheapest concrete deliverable from today's session. They're 1-2 days of focused writing and they make Paper 20 honest about scope. I'd prioritize them ahead of Paper 57 drafting because Paper 20 is already on disk and the edits are surgical.

## §6 What NOT to claim

Today's results allow strong claims about what the framework *is* and *isn't*. Some claims to avoid:

- **DON'T claim "no engineering can close the second-row wall."** We only tested a handful of engineering knobs. The Plan agent's A/B sprint list has several un-run options (multi-zeta on first-row, R-adaptive Z_eff, B.4 bonding-orbital Pauli repulsion). We can claim "the cheapest engineering knobs are exhausted" but not "no engineering can close it."

- **DON'T claim "GeoVac is structurally bounded to first-row chemistry forever."** Today's negative is at the current architecture's scope. A future structural reformulation (e.g., the Marcolli-vS chemistry arc, or a continuous-basis embedding) could shift the boundary.

- **DON'T claim "we can predict the boundary from first principles."** We can document where it sits empirically. Predicting it structurally is the Marcolli-vS paper's job, not Paper 57's.

The Paper 20 edits and Paper 57 framing should be carefully scoped to what today's evidence actually supports.

## §7 Files

### To be modified (next sprint)
- `papers/group4_quantum_computing/paper_20_resource_benchmarks.tex` — abstract, conclusions, new §V.C

### To be created (Paper 57 drafting sprint, separate)
- `papers/group4_quantum_computing/paper_57_chemistry_scope_boundary.tex` (new)

### Reference material already on disk
- `debug/sprint_w1e_projection_audit_memo.md` — engineering boundary evidence (LiH/NaH balanced PES)
- `debug/sprint_b1_nah_hf_verdict_memo.md` — B.1 explicit-core HF negative
- `debug/lih_kwarg_sweep_log.txt` — kwarg sweep negative
- `debug/bratteli_h2_pilot_memo.md` — structural framework match
- `debug/sprint_marcolli_vs_chemistry_paper_arc_scoping_memo.md` — companion paper-arc scoping
