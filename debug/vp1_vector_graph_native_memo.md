# VP-1 Sprint Memo: Vector-photon promotion of graph-native QED

**Date:** 2026-05-02
**Sprint:** VP-1 (vector-photon graph-native vertex correction)
**Status:** **NEGATIVE** with three concrete structural findings

## Goal

Test whether replacing the SCALAR (Fock-edge 1-cochain) photon in the
graph-native QED vertex correction with a VECTOR photon ((q, m_q) modes
with Wigner 3j vertex coupling) closes the `C × F2` negative result
(CLAUDE.md `single-constant graph-to-continuum QED projection` failed
mode).

Hypothesis: with vector photons, `C × F2_vector` should converge to
`α/(2π) ≈ 1.16140973×10⁻³` as n_max grows.

## Construction

Hybrid graph-native + vector-photon pipeline:

  - **Electron propagator**: graph-native Dirac graph at `t=0`,
    `G_e[i,i] = 1/λ_i` with Camporesi-Higuchi `λ_n = ±(n+1/2)`. Exact
    rational, π-free.
  - **Photon propagator**: vector-photon `G_γ(q) = 1/[q(q+2)]` from
    Hodge-1 spectrum on S³.
  - **Vertex**: `dirac_vertex_coupling` from `geovac/vector_qed.py` with
    Wigner 3j, E1 parity, j-triangle. Contains `1/√(4π)`.
  - **Lambda_total**: `Λ[a,c] = Σ_{b,k} V[a,b,k]·G_e[b,b]·G_γ(q_k)·V[c,b,k]`
  - **F2 extraction**: graph-native style, `F2 = Tr(Λ)/Tr(V_bare·G_e)`.
    But `V_bare` diagonal is **identically zero** by the Dirac spinor
    phase constraint (Furry), so the original normalization is degenerate.
    Three fallback normalizations are reported (TrGe, ||V||²_F, N_e).

## Results

### Exact sympy values (confirms π-content of vector promotion)

| n_max | Tr(Λ) (exact) | Tr(G_e) (exact) | F2_vector = Tr(Λ)/Tr(G_e) (exact) | Float |
|:-----:|:--------------|:----------------|:----------------------------------|:------|
| 2     | `4/(5π)`      | `44/15`         | `3/(11π)`                         | 0.08681179 |
| 3     | `319/(42π)`   | `488/105`       | `1595/(976π)`                     | 0.52018880 |

**Key observation**: F2_vector contains an explicit `1/π` from the
spherical-harmonic normalization `1/√(4π)` in the vertex coupling.
Vector promotion injects π immediately — this is the calibration
tier of Paper 18 by construction. The graph-native scalar F2 at n_max=2
was `5√2/3` (algebraic, π-free); the vector version at n_max=2 is
`3/(11π)` (single-π, transcendental). Vector promotion is exactly
the calibration step that turns graph-intrinsic algebraic F2 into
continuum-style transcendental F2.

### F2 convergence table (numpy, three normalization choices)

| n_max | F2_scalar (CLAUDE) | F2_vec / Tr(G_e) | F2_vec / ‖V‖²_F | F2_vec / N_e |
|:-----:|:------------------:|:----------------:|:-----------------:|:------------:|
| 2     | 2.357023           | 8.6812e-02       | 3.6364e-02        | 2.5465e-02   |
| 3     | 1.873000           | 5.2019e-01       | 1.6511e-02        | 8.6344e-02   |
| 4     | 1.589000           | 1.5555e+00       | 7.6620e-03        | 1.6658e-01   |

Power-law fits (n_max=2,3,4):
  - `F2_vec / Tr(G_e)`:  α = **+4.18**, R² = 0.998 (INCREASING with n_max)
  - `F2_vec / ‖V‖²_F`:   α = **−2.23**, R² = 0.992 (DECREASING with n_max)
  - F2_scalar reference: α = **−0.573** (CLAUDE.md, n_max=2..6)

The TrGe-normalized F2_vector grows fast with n_max — opposite sign
to the scalar's −0.573. The ‖V‖²_F-normalized F2_vector decreases
faster than the scalar's −0.573. **Neither vector normalization matches
the scalar convergence trend, and neither converges toward `α/(2π)`.**

### C × F2_vector vs α/(2π) at n_max=3

The graph-native projection constant `C(n_max=3) = 50471424/1779441125
≈ 2.836×10⁻²` was extracted from the SCALAR-VP / continuum-VP ratio.

| Quantity                  | Value      |
|:--------------------------|:-----------|
| C(n_max=3)                | 2.836e-02 |
| F2_vec/TrGe at n_max=3    | 5.202e-01 |
| C × F2_vec (TrGe)         | **1.475e-02** |
| α/(2π)                    | **1.161e-03** |
| Ratio = (C × F2_vec) / (α/(2π)) | **12.70** |

C × F2_vector overshoots α/(2π) by **12.7×** at n_max=3. Compare to
the scalar baseline `C × F2_scalar = 0.0284 × 1.873 = 0.0531` which
overshoots α/(2π) by 45.7× — so vector promotion gets closer (3.6×
better) but still nowhere near 1.0.

## Verdict: NEGATIVE

The hypothesis "vector promotion closes the C × F2 negative result"
is **NOT confirmed**. Three concrete structural reasons:

### 1. Wrong-sign convergence (TrGe normalization)

F2_vec(TrGe) at n_max=2,3,4 INCREASES with n_max (slope +4.18 in log-log
fit). The scalar F2(n_max) DECREASES (slope −0.573). Vector promotion
flips the sign of the convergence direction. The scalar pendant-edge
mechanism that drove F2_scalar → 0 is absent in the vector version because:

  - In the scalar graph, GS is a pendant vertex, the unique edge e_0
    connecting it has `G_γ[e_0,e_0] = (n_max−1)/n_max` from the path-graph
    Laplacian, giving `Σ(GS) = 2(n_max−1)/n_max`.
  - In the vector vertex, GS structural zero is from the **angular
    momentum** triangle (l_a + l_b + q parity), not from edge-graph
    topology. The pendant-edge mechanism does not transfer.

The two GS protection mechanisms produce different convergence
behaviors. The scalar F2 → 0 fast; the vector F2 has no such
suppression and grows with the size of the angular momentum sum
(more allowed (q, m_q) channels).

### 2. C is the WRONG projection constant for vector vertex

The projection constant `C(n_max=3) = 50471424/1779441125` was extracted
from `C = Π_continuum_truncated / Tr(Π_graph_scalar)`. This compares
SCALAR graph-native VP (Tr(Π_graph) at n_max=3 = 3872/735) to
truncated SO(4)-channel-counted continuum VP (538361856/3603000625).

The vector vertex produces a DIFFERENT graph-native quantity (a vector
QED vertex correction Λ, not a vacuum polarization Π). Mixing
`C_VP_scalar × F2_vector` is dimensionally and structurally
inconsistent — confirming the CLAUDE.md "Single-constant graph-to-
continuum QED projection" failed mode where different diagram
topologies require different projection constants.

The correct vector-aware projection would be a NEW constant
`C_F2_vector` extracted from the ratio of the truncated continuum
F2 (Schwinger spectral mode sum) to the graph-vector-vertex F2.
This reduces to a tautology: since F2_vector contains `1/π` already,
any rescaling that produces α/(2π) at one n_max won't generalize
without becoming C(n_max).

### 3. Vector promotion injects π already (per Paper 18)

The exact sympy values at n_max=2 (`F2 = 3/(11π)`) and n_max=3
(`F2 = 1595/(976π)`) confirm that the vector vertex **directly
injects π** through the spherical-harmonic normalization `1/√(4π)`.
This is the **calibration tier** of Paper 18: continuum-style
transcendental content is paid for at the vertex level by importing
angular-momentum quantum numbers `(q, m_q)`.

Reaching `α/(2π)` would require:

  (a) A vector vertex (this sprint provides it: 8/8 selection rules,
      π-content, GS structural zero from angular-momentum parity);
  (b) An α-dependent coupling at the action level (NOT in the graph
      vertex itself, but in the Wilson action coupling β = 1/g² that
      appears in Paper 30's SU(2) gauge construction);
  (c) A correct projection from finite-graph mode sums to infinite
      Hurwitz-zeta-regularized continuum sums — this is where π² and
      the Schwinger result α/(2π) live.

VP-1 supplies (a) but does not address (b) or (c). The gap between
`C × F2_vec = 1.5e-2` and `α/(2π) = 1.2e-3` is the (b)+(c) deficit.

## Side findings (positive)

  - **8/8 selection rules verified at the vertex level** (Lambda_GS_block_max = 0
    exactly at all n_max=2,3,4; V_bare_diag = 0 exactly by Furry/Dirac spinor
    phase). The vector vertex is correct.
  - **Sparsity is high**: vertex tensor 93.3% / 94.8% / 95.8% sparse at
    n_max=2/3/4 — Wigner 3j enforces angular momentum conservation
    aggressively.
  - **Exact sympy results at n_max=2,3** confirm π-injection at the
    one-loop vertex level: `F2_vec = 3/(11π) and 1595/(976π)`. These
    are clean rational multiples of `1/π`, which is the simplest
    possible π-content (one factor of π for one loop, matching the
    Schwinger continuum result `α/(2π)` having one factor of π).

## Recommendation

The hypothesis "vector promotion closes the C × F2 negative result"
is closed as a clean negative. Two follow-ups would be productive:

  - **VP-2 (recommended)**: derive the n_max-dependent
    `C_F2_vector(n_max)` by computing the truncated continuum F2
    spectral mode sum (Schwinger evaluation on S³) at matching n_max,
    and check whether `C_F2_vector × F2_vec` converges to `α/(2π)`.
    This requires defining the matching truncation, but the data
    (F2_vec at n_max=2,3,4) is in this sprint.
  - **VP-3 (longer term)**: VP-1's exact result `F2_vec = 3/(11π)` at
    n_max=2 hints that the natural normalization is the bare
    `Tr(Λ)·(4π)/(some ℚ)` extraction, where the `4π` from the vertex
    normalization cancels and gives a clean rational. Worth one PSLQ
    diagnostic if a rational pattern emerges across `4π·F2_vec(n_max)`.

## Files

  - `debug/vp1_vector_graph_native.py` — implementation
  - `debug/data/vp1_vector_graph_native.json` — full numerical data
  - `debug/vp1_vector_graph_native_memo.md` — this memo
