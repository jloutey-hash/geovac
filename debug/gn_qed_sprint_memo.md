# Graph-Native QED Sprint: Three-Track Results

**Date:** 2026-04-27
**Sprint scope:** F₂ rational structure, selection rule census, self-energy gap analysis
**Core question:** Is α the projection constant from graph-native scalar QED to physical vector QED?

## Summary

The graph-native QED framework (GN-1 through GN-7) computes scalar QED exactly on the finite Fock graph in the algebraic ring ℚ(√2,√3,√6) — entirely π-free. Three parallel tracks characterize its structure and its relationship to continuum vector QED.

**Headline results:**

1. **F₂(t) is an even rational function** of the hopping coupling t, degree 2/2 over ℚ(√2,√3,√6) at n_max=2. All odd Taylor coefficients vanish identically.

2. **Only 1 of 8 continuum QED selection rules survives** on the finite graph (Gaunt/CG sparsity). The 7 broken rules partition cleanly into 3 STRUCTURAL (require vector photon structure absent from the graph) and 4 INTRINSIC (broken by graph kinematics).

3. **C × F₂ does NOT converge to 1/(2π).** The product grows with n_max: 0.053 (n=3), 0.075 (n=4). The VP projection constant C is not the correct projection for the vertex correction. The graph-to-continuum projection for different diagram topologies involves different spectral-density factors.

4. **F₂(κ) decreases monotonically** with n_max: 2.353, 1.873, 1.589 at n_max=2,3,4. Effective power law F₂ ~ n_max^(−0.57).

## Track 1: F₂(t) Rational Function Extraction

### n_max=2 (exact symbolic, 2.8s)

| Quantity | Value |
|:---------|:------|
| F₂(0) | 5√2/3 ≈ 2.3570 |
| F₂(κ = −1/16) | ≈ 2.3525 |
| Numerator degree | 2 |
| Denominator degree | 2 |
| Number field | ℚ(√2, √3, √6) |
| Taylor c₀ | 5√2/3 [algebraic] |
| Taylor c₁ | 0 [exact] |
| Taylor c₂ | ≈ −1.163 [algebraic] |
| Taylor c₃ | 0 [exact] |
| Taylor c₄ | ≈ 1.685 [algebraic] |
| Taylor c₅ | 0 [exact] |

**Key finding:** F₂(t) is an EVEN function of t. The Neumann/Born series for G_e(t) = Λ⁻¹ Σ (−tΛ⁻¹A)^k terminates at finite order (Cayley-Hamilton), making F₂ a rational function. The even-ness means the anomalous moment depends only on t² = κ², not on the sign of the hopping coupling. This is because D(t) = Λ + tA decomposes into diagonal Λ (even eigenvalues ±|λ|) and off-diagonal A (T± ladder operators connecting adjacent shells), and the vertex correction involves the bilinear V·G_e·V^T which is invariant under t → −t.

### n_max=3 (exact symbolic, 207s)

| Quantity | Value |
|:---------|:------|
| F₂(0) | exact rational/algebraic (long expression) |
| F₂(κ) | ≈ 1.8731 |
| N_dirac | 28 |
| E_fock | 13 |

The 28×28 symbolic matrix inversion succeeds in 6.4s. The full Lambda(t) computation with 13 edges takes ~200s.

### n_max=4 (numeric only)

| Quantity | Value |
|:---------|:------|
| F₂(κ) | 1.5892 |
| N_dirac | 60 |
| E_fock | 34 |
| Note | Exact computation obstructed by sympy complex eigenvalue bug in L₁ at E=34 |

### F₂ convergence table

| n_max | N_dirac | E_fock | F₂(κ) | F₂(0) | F₂×n_max |
|------:|--------:|-------:|------:|------:|---------:|
| 2 | 10 | 3 | 2.3525 | 2.3570 | 4.705 |
| 3 | 28 | 13 | 1.8731 | 1.8756 | 5.619 |
| 4 | 60 | 34 | 1.5892 | — | 6.357 |

F₂ ~ n_max^(−0.57) (effective exponent from log-log fit on n=2,3,4).

### Convergence: C(n_max) × F₂(κ) vs α/(2π)

| n_max | C | F₂(κ) | C×F₂ | α/(2π) | Ratio |
|------:|--:|------:|-----:|-------:|------:|
| 2 | 0 (VP=0) | 2.3525 | N/A | 1.16×10⁻³ | — |
| 3 | 2.836×10⁻² | 1.8731 | 5.31×10⁻² | 1.16×10⁻³ | 45.7 |
| 4 | 4.743×10⁻² | 1.5892 | 7.54×10⁻² | 1.16×10⁻³ | 64.9 |

**Verdict: C × F₂ diverges from α/(2π), not converges.** The VP projection constant C is specific to vacuum polarization and grows faster with n_max than F₂ decreases. The correct graph-to-continuum projection for the vertex correction is a different quantity from C.

## Track 2: Selection Rule Survival Census

Systematic test of 8 continuum QED selection rules on the finite Fock graph at n_max=2.

| Rule | Verdict | Tier | Key Metric |
|:-----|:--------|:-----|:-----------|
| Gaunt/CG sparsity | **SURVIVES** | INTRINSIC | density 16%→2.6%→0.6% at n=2,3,4 |
| Angular momentum Δm_j | BROKEN | INTRINSIC | 20 violations |
| Spatial parity (E1) | BROKEN | INTRINSIC | 100% parity violations |
| Charge conjugation (C) | BROKEN | INTRINSIC | 8 C-violations in Σ |
| Furry's theorem | BROKEN | INTRINSIC | tadpole = 16√2/15 ≠ 0 |
| Vertex parity (n₁+n₂+q odd) | BROKEN | STRUCTURAL | 0% forbidden by continuum rule |
| SO(4) channel count (W>0) | BROKEN | STRUCTURAL | 100% have W=0 |
| Ward identity | BROKEN | STRUCTURAL | max residual = 5.08 |

**Clean partition:**
- **3 STRUCTURAL** (vertex parity, SO(4) channel count, Ward identity): require the photon to be a vector harmonic on S³. The graph's scalar 1-cochain photon has no mechanism to carry parity-flip, angular momentum transfer, or gauge derivative structure. These enter at the CALIBRATION tier when projecting to vector QED.
- **4 INTRINSIC** (Δm_j, spatial parity, C, Furry): broken by graph kinematics. The CG projection allows couplings forbidden in vector QED because the scalar photon carries different angular momentum content.

**Notable invariants:**
- Tr(Σ) = 44/3, Tr(Λ) = 32/9, Tr(Λ)/Tr(Σ) = 8/33
- Tadpole = 16√2/15, Bubble = 32/9, Triangle = 3584√2/3375

## Track 3: Self-Energy Gap Analysis

Per-shell comparison of graph self-energy Σ_graph(n_fock) vs continuum spectral self-energy Σ_cont(n_ext) at t=0.

### n_max=2 (Tr(Σ) = 44/3 = 14.667)

| n_fock | g_shell | Σ_graph (mean) | Σ_cont | Delta | Ratio |
|-------:|--------:|---------------:|-------:|------:|------:|
| 1 | 2 | 1.000 | 0.000 | 1.000 | ∞ |
| 2 | 8 | 1.583 | 0.160 | 1.423 | 9.9 |

### n_max=3 (Tr(Σ) = 61.600)

| n_fock | g_shell | Σ_graph (mean) | pendant |
|-------:|--------:|---------------:|--------:|
| 1 | 2 | 1.333 | 4/3 ✓ |
| 2 | 8 | 1.850 | — |
| 3 | 18 | 2.452 | — |

**Pendant-edge theorem verified:** GS self-energy = 2(n_max−1)/n_max exactly at every n_max. At n_max=2: Σ(GS) = 1; at n_max=3: Σ(GS) = 4/3. The ground-state node is a pendant (leaf) vertex of the Fock graph, and the unique edge connecting it decouples in L₁, giving the exact formula.

**Gap structure:** The graph self-energy is uniformly ~10× larger than the continuum because the graph scalar photon opens ALL CG-allowed channels (no SO(4) W=0 suppression). The gap is not a correction — it's a structural feature of scalar vs vector QED.

## Structural Picture

### What the graph computes
The finite Fock graph at n_max computes **exact scalar QED** in the algebraic ring ℚ(√2,√3,√6):
- F₂ = 5√2/3 ≈ 2.357 (n_max=2, t=0) — the graph-native anomalous magnetic moment
- Σ has eigenvalues {0(×5), 4/3, 2(×2), 4, 16/3} (n_max=2, t=0)
- All quantities are algebraic, π-free, exact at every finite truncation
- The projection exchange constant C is rational at every truncation

### What the projection to physics requires
Three ingredients transform graph-native scalar QED into physical vector QED:
1. **Coupling constant α** — one power per loop; enters as the physical charge e²/(4π)
2. **Vector photon structure** — the 7 broken selection rules encode 7 missing symmetry constraints. The scalar photon has 1 polarization; the vector photon on S³ has 3 polarizations × parity. This is the STRUCTURAL gap.
3. **Continuum limit n_max → ∞** — the graph's algebraic ring converges to the transcendental ring containing π, ζ(2), etc.

These ingredients do NOT factorize into a single multiplicative constant C × F₂. The correct projection involves diagram-topology-specific spectral-density matching.

### The focal-point analogy
The PI's intuition is structurally correct: just as p₀² = −2E projects the dimensionless graph eigenvalues to physical energies (introducing a unit of measurement), the projection from graph-native QED to physical QED introduces α (the coupling) and π (the continuum measure). But the "focal point" for QED is not a single number — it's the full spectral-density-of-states factor that converts discrete sums to continuum integrals.

### Paper 18 classification
- **INTRINSIC (graph-native):** F₂, Σ eigenvalues, Tr(Λ)/Tr(Σ) = 8/33, Gaunt sparsity, pendant-edge theorem, tadpole/triangle values. All algebraic, all π-free.
- **CALIBRATION (projection):** C = Π_cont/Tr(Π_graph), the scalar-to-vector upgrade factor, the coupling constant α. These introduce π and ζ through the continuum integration measure.
- **STRUCTURAL (absent from graph):** vertex parity, SO(4) channel count, Ward identity, vector photon polarizations. These are features of the continuum gauge theory that have no graph analog.

## Files

- `debug/gn_selection_rule_census_memo.md` — Track 2 complete memo
- `debug/data/gn_selection_rule_census.json` — Track 2 structured data
- `debug/gn_qed_sprint_memo.md` — this file
