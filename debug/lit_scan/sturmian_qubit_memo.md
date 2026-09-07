# Lit scan — Sturmian/isoenergetic basis as a qubit-encoding / block-encoding target

**Date:** 2026-09-06
**Scope:** adversarial prior-art / novelty scan for Paper 60's distinctive claims. Refreshes and
extends the prior internal scout `debug/sturmian_quantum_priorart.md` (2026-08-17), which was framed
around the QPE/root-search *algorithm*; this pass is framed around the **1-norm / block-encoding
resource** claims. Web-verified citations only; anything not re-confirmed this pass is marked.

## The three claims under test (verbatim from task)

1. **DISTINCTIVE / load-bearing.** Reading Avery's isoenergetic Coulomb-Sturmian secular equation as
   a QC target: for **atoms** the metric drops out (ordinary eigenproblem) and the block-encoding LCU
   1-norm scales **SUBLINEARLY** (~K^0.84) in configuration/basis number K — structurally better than
   Gaussian-basis encodings whose 1-norm grows faster.
2. For **molecules** the sublinearity is lost (collective-scale diagonal T^0 becomes non-diagonal);
   molecules pay ordinary polynomial 1-norm (~n_orb^2.2); the distinctive value there is metric
   "levers" (flat gerade sector; large-R conditioning), not a sublinear matrix.
3. Composed qubit Hamiltonians with **Pauli-term count exactly linear** in qubit number across
   molecules at fixed basis (N_Pauli = 27.9·Q), from angular (Gaunt/3j) selection-rule sparsity.

---

## Ranked papers (most dangerous first)

### A. Babbush, Berry, McClean, Neven — "Quantum Simulation of Chemistry with Sublinear Scaling in Basis Size"
- npj Quantum Information **5**, 92 (2019). arXiv:1807.09802. **[VERIFIED this pass]**
- **What it says:** first-quantized **plane-wave** dynamics; total **gate complexity** Õ(N^{1/3} η^{8/3}),
  i.e. *sublinear in the number of basis functions N*. The sublinear quantity is the **gate count**
  (interaction-picture Trotter), NOT an LCU 1-norm, and NOT for atoms specifically.
- **Touches:** claim 1 — **on WORDING ONLY**. Same headline phrase ("sublinear scaling in basis size"),
  entirely different quantity (per-evolution gate complexity vs. LCU 1-norm of a fixed-atom secular
  matrix) and different basis (plane waves vs. Coulomb-Sturmians). **This is the most dangerous paper:**
  a referee who reads only the abstract will file Paper 60 as derivative of this famous result.
  Paper 60 must disambiguate explicitly (quantity = λ = ‖M‖₁, not gate count; system = atom at fixed
  size, growing config count, not first-quantized dynamics).
- **Verdict contribution:** NOT a substantive collision. Wording collision → mandatory disambiguation.

### B. Koridon, Yalouz, Senjean, Buda, O'Brien, Visscher — "Orbital transformations to reduce the 1-norm of the electronic structure Hamiltonian for quantum computing applications"
- Phys. Rev. Research **3**, 033127 (2021). arXiv:2103.14753. DOI 10.1103/PhysRevResearch.3.033127. **[VERIFIED this pass]**
- **What it says:** the electronic-structure **1-norm λ** and its **scaling with basis/orbital size**
  is basis-dependent; localization lowers it. Reported exponents (per search summary of the full
  text — treat exponent digits as reported-not-independently-recomputed): canonical MOs ~N^2.2–2.3;
  **localized orbitals ~N^1.34–1.38** (H-chains, alkanes). Gaussian orbitals.
- **Touches:** claim 1 (directly) and claim 2. This is the **correct comparator** for "basis choice
  changes the 1-norm exponent." Its best is **N^1.34 — still SUPERLINEAR** and for a *growing molecule*
  (orbital count), not a fixed atom with growing config count. Our K^0.84 (exponent < 1) is a
  qualitatively stronger regime than anything here.
- **Verdict contribution:** **ADJACENT.** Establishes that "basis lowers λ-exponent" is a known lever,
  but **no sublinear (exponent < 1) λ is reported by anyone.** Paper 60 should cite this as the
  best-known prior lever and state the K^0.84 is below the N^1.34 floor.

### C. Rajchel-Mieldzioć, Pliś, Żak — "Quantum algorithm for solving generalized eigenvalue problems with application to the Schrödinger equation"
- Phys. Rev. Research (2025/26). arXiv:2506.13534. **[carried from prior scout 2026-08-17; not re-fetched this pass]**
- **What it says:** solves Hψ=ESψ (pseudospectral/distributed-Gaussian) by an outer scan for the
  parameter α where the block-encoded residue M(α) becomes singular; explicitly avoids S⁻¹/κ(S).
  Combined cost Õ(√(NK)) via fused amplitude amplification.
- **Touches:** claim 1(b) — closest structural analog to "energy-as-eigenvalue." **Key difference:**
  the scanned α *is the eigenvalue E itself* (a null-space/singularity test), NOT the isoenergetic
  metric-free eigenproblem where the eigenvalues are pκ=√(−2E) *directly* and there is no metric.
- **Verdict contribution:** **ADJACENT**, not collision. Does not exploit the isoenergetic inversion.

### D. Georges, Bothe, Sünderhauf, Berntson, Izsák, Ivanov — "Quantum simulations of chemistry in first quantization with any basis set"
- npj Quantum Information **11**, 55 (2025). arXiv:2408.03145. **[VERIFIED this pass]**
- **What it says:** first-quantized ground-state chemistry with an **arbitrary** basis (molecular
  orbitals / dual plane waves); asymptotic Toffoli speedup vs. second quantization.
- **Touches:** claim 1/2 tangentially. "Any basis" is **generic (MO / dual-plane-wave), not Sturmian
  or hyperspherical or momentum-space Fock.** No sublinear-λ-for-atoms claim.
- **Verdict contribution:** ADJACENT (shows "basis freedom in first quantization" is live), no collision.

### E. Babbush, Somma, Rubin — "Fast quantum simulation of electronic structure by spectrum amplification"
- Phys. Rev. X (2025). arXiv:2502.15882. (research.google; PRX 10.1103/pb2g-j9cw). **[VERIFIED this pass]**
- **What it says:** sum-of-squares Hamiltonians let ground-state estimation cost scale as √(2Λ·E_gap)
  in the block-encoding normalization Λ — an *effective* reduction of the 1-norm penalty, but by a
  new **algorithmic** device (spectrum amplification + a DF/THC-interpolating factorization), not by
  a **basis choice** that makes ‖M‖₁ itself sublinear.
- **Touches:** claim 1. Different mechanism (amplify spectrum, keep basis) vs. ours (change the posing
  so the matrix norm is small). Not a collision; worth citing as "the other way people beat λ."
- **Verdict contribution:** ADJACENT.

### F. Berry, Gidney, Motta, McClean, Babbush — "Qubitization of Arbitrary Basis Quantum Chemistry Leveraging Sparsity and Low Rank Factorization"
- Quantum **3**, 208 (2019). (arXiv:1902.02134 — **arXiv ID UNVERIFIED this pass**; venue verified via quantum-journal.org.)
- **What it says:** arbitrary-basis qubitization; λ reduced via sparsity + low-rank factorization.
  Baseline for "structure in the Coulomb operator lowers λ."
- **Touches:** claim 1/3 (term-count & 1-norm sparsity). Gaussian/arbitrary basis, no Sturmian.
- **Verdict contribution:** ADJACENT context; the general "sparsity lowers cost" lineage.

### G. Loaiza, Izmaylov (BLISS) — "Global Minimization of Electronic Hamiltonian 1-Norm via Linear Programming in the Block Invariant Symmetry Shift"
- J. Chem. Theory Comput. **21**, 703 (2025). arXiv:2409.18277. **[VERIFIED this pass — search]**
- **What it says:** lowers λ by symmetry-shifting the Hamiltonian off unphysical particle-number
  sectors via LP. Another λ-reduction lever, orthogonal to basis choice.
- **Touches:** claim 1/2 context. Not a collision.

### H. Herbst, Avery, Dreuw — "Quantum chemistry with Coulomb Sturmians …"
- Phys. Rev. A **99**, 012512 (2019). arXiv:1811.05777. **[VERIFIED this pass]**
- **What it says:** classical HF-level Coulomb-Sturmian basis construction/convergence. Same Avery
  lineage as the isoenergetic method. **Zero quantum-computing content.**
- **Touches:** confirms the *basis* side has no QC prior art either.

### Supporting (from prior scout, carried; GEVP-conditioning line, not re-fetched this pass)
- Liang, Shen, Li, Fei — Quantum Inf. Process. **21**, 23 (2022), arXiv:2112.02554 **[seen in search this pass]**:
  quantum GEVP with explicit κ(S) multiplicative cost — the metric penalty our atomic case *removes*.
- Baek et al. (NOQE) — PRX Quantum **4**, 030307 (2023), arXiv:2205.09039: κ(S) reappears as amplified
  measurement error in the non-orthogonal near-term route.

---

## Verdicts

**Claim 1 (atomic, metric-free, sublinear λ ~ K^0.84): OPEN.**
Decision-gate answers: (a) NO prior work uses a Sturmian/hyperspherical basis for qubit simulation —
the QC electronic-structure literature is Gaussian / plane-wave / dual-plane-wave / MO only (A, D, F, G);
the only Coulomb-Sturmian QC-adjacent paper (H) is purely classical. (b) NO prior work exploits the
isoenergetic *metric-free* (energy-as-eigenvalue) form as an encoding target — the closest (C) scans
the eigenvalue as a singularity test and keeps the GEVP framing; the fault-tolerant GEVP solvers
(Liang) pay κ(S) explicitly. (c) NO prior work reports a **sublinear (exponent < 1)** 1-norm from
basis choice: the best documented basis lever is Koridon's **N^1.34** (B), still superlinear and for a
*growing molecule*, not a fixed atom. The famous "sublinear scaling in basis size" (A) is a **gate**
count in plane waves — a wording collision, not a substantive one. → **OPEN, no collision.** The one
required action is to disambiguate against A so the headline is not read as derivative.

**Claim 2 (molecular ~n_orb^2.2, levers not a sublinear matrix): ADJACENT / consistent — honest, not novel-as-such.**
n_orb^2.2 sits squarely in the known band (Koridon canonical ~N^2.2; localized ~N^1.34). Paper 60
already frames molecules as paying ordinary polynomial cost, which is exactly right. The gerade-sector
/ large-R metric-conditioning "levers" are adjacent to the GEVP-conditioning (Liang) and
λ-reduction (F, G, BLISS) literature; distinctive in detail (symmetry-flat metric block) but not a
standalone novelty claim. No collision.

**Claim 3 (N_Pauli = 27.9·Q exactly linear across molecules, from Gaunt/3j sparsity): ADJACENT.**
Term-count sparsity from structure is a known genus — plane-wave dual basis gives O(N^2) term count,
and sparsity/low-rank qubitization (F) exploits Coulomb structure. **Exact Q-linearity across molecules**
from angular selection rules is a distinctive empirical regularity and I found no direct prior
statement of it. Caveat (already in CLAUDE.md §1.5): Pauli **term count** is a weaker metric than λ —
it does not set qubitization cost — so this claim is real but should not be over-weighted against the
1-norm literature. No collision.

## Most dangerous paper
**A — Babbush, Berry, McClean, Neven, npj QI 5:92 (2019), "…Sublinear Scaling in Basis Size."**
Not because it does what Paper 60 does (it does not — gate count in plane waves, not λ for atoms),
but because it **owns the exact headline phrase**. Any reader/referee pattern-matches on it first.
Paper 60's abstract/intro must state, in one sentence, that the sublinear quantity here is the LCU
1-norm ‖M‖₁ of the fixed-size atomic isoenergetic secular matrix (growing configuration count), which
is a different object from the first-quantized per-evolution gate complexity of Babbush 2019.

## Net
The distinctive claim (1) survives the scan: **genuinely OPEN, no substantive prior-art collision.**
The gap the prior scout flagged (Sturmian/isoenergetic × quantum) is confirmed still open on the
1-norm axis specifically, and sharpened: **nobody reports a sub-*linear* (exponent < 1) 1-norm from
basis choice for atoms; the best prior basis lever is N^1.34.** The only exposure is the wording
collision with the 2019 "sublinear scaling in basis size" title, which is a presentation fix, not a
novelty defect.
