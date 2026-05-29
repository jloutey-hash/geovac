# Sprint G6-Diag first-pass — graviton diagnostic on the CH Dirac substrate

**Date:** 2026-05-28
**Path:** Gravity arc, G6-Diag first-pass (scoping-then-implementation). Tests whether the second variation of the spectral action around the CH Dirac background carries eigenmodes with $(j_L, j_R) = (1, 1)$ "graviton-irrep" content of $SO(4) = SU(2)_L \times SU(2)_R$.
**Verdict:** **POSITIVE-FIRST-PASS.** At $n_{\max} = 1$, four (1,1)-irrep candidate off-blocks identified with NONZERO $S^{(2)}$ eigenvalues — graviton-irrep modes exist and propagate at the substrate level. Necessary condition for graviton dynamics is met. Full sufficient analysis (Fierz-Pauli kinetic, gauge invariance, propagator structure) is the multi-month G6-Diag-Full / G6 sprint.

## 1. Setup

**Background:** truncated CH Dirac $D_0$ on $\mathcal{H}_{n_{\max} = 1}$, dim 16. Four sectors:

| Sector | $\lambda$ | $(j_L, j_R)$ | dim |
|---|---|---|---|
| $S_1$ | $+3/2$ | $(1/2, 0)$ | 2 |
| $S_2$ | $-3/2$ | $(0, 1/2)$ | 2 |
| $S_3$ | $+5/2$ | $(1, 1/2)$ | 6 |
| $S_4$ | $-5/2$ | $(1/2, 1)$ | 6 |

**Perturbation:** $D = D_0 + \epsilon V$ where $V$ is Hermitian. Spectral action expands as
$$S[D_0 + \epsilon V] = S[D_0] + \epsilon S^{(1)}[V] + \tfrac{\epsilon^2}{2} S^{(2)}[V, V] + O(\epsilon^3)$$

For Gaussian cutoff $f(x) = e^{-x}$ at $\Lambda^2 = 6$.

**Question:** does the quadratic form $S^{(2)}$ on the 256-dim real space $\text{Herm}(\mathcal{H})$ have eigenmodes with $(j_L, j_R) = (1, 1)$ content?

## 2. Analytical structure of $S^{(2)}$

In the $D_0$ eigenbasis, the quadratic form is DIAGONAL: each Hermitian basis matrix is an eigenvector. The eigenvalues split into two classes:

**Within-sector (degenerate background eigenvalue):** for any Hermitian mode $E$ supported entirely within sector $S_i$ (eigenvalue $\lambda_i$),
$$S^{(2)}[E, E] = a_{\lambda_i}\left(\frac{4\lambda_i^2}{\Lambda^4} - \frac{2}{\Lambda^2}\right) =: A_{\lambda_i}$$
where $a_\lambda = e^{-\lambda^2/\Lambda^2}$.

**Cross-sector (different background eigenvalues $\lambda_a \neq \lambda_b$):** for a real-Hermitian off-diagonal mode connecting sectors $S_a$ and $S_b$,
$$S^{(2)} = -\frac{2}{\Lambda^2}\,\frac{\lambda_a a_{\lambda_a} - \lambda_b a_{\lambda_b}}{\lambda_a - \lambda_b} =: B_{ab}$$

(The imaginary-Hermitian modes have the same eigenvalue.)

## 3. Numerical verification

Implemented in `debug/g6_diag_quadratic_form.py` and verified via 5-point finite-difference stencil on $g(\epsilon) = \mathrm{Tr}\,e^{-(D_0 + \epsilon V)^2/\Lambda^2}$ for representative test modes. **All analytical eigenvalues match finite-difference to rel diff $\sim 10^{-9}$** (limited by float64, not by the formulas).

**Eigenvalue table** at $n_{\max} = 1$, $\Lambda^2 = 6$ (all 10 sector pairs):

| Sector pair | $\lambda_a$ | $\lambda_b$ | Eigenvalue |
|---|---|---|---|
| $S_1 \times S_1$ | $+1.50$ | $+1.50$ | $-0.0573$ |
| $S_1 \times S_2$ | $+1.50$ | $-1.50$ | $-0.2291$ |
| $S_1 \times S_3$ | $+1.50$ | $+2.50$ | $+0.0496$ |
| $S_1 \times S_4$ | $+1.50$ | $-2.50$ | $-0.1594$ |
| $S_2 \times S_2$ | $-1.50$ | $-1.50$ | $-0.0573$ |
| $S_2 \times S_3$ | $-1.50$ | $+2.50$ | $-0.1594$ |
| $S_2 \times S_4$ | $-1.50$ | $-2.50$ | $+0.0496$ |
| $S_3 \times S_3$ | $+2.50$ | $+2.50$ | $+0.1274$ |
| $S_3 \times S_4$ | $+2.50$ | $-2.50$ | $-0.1176$ |
| $S_4 \times S_4$ | $-2.50$ | $-2.50$ | $+0.1274$ |

## 4. $SO(4)$ classification: which off-blocks carry $(1, 1)$ graviton irrep?

For each off-block $S_i \otimes S_j^*$, decompose into $SU(2)_L \times SU(2)_R$ irreps via tensor product:

| Sector pair | $SU(2)_L$ decomp | $SU(2)_R$ decomp | $(1, 1)$ mult | Eigenvalue |
|---|---|---|---|---|
| $S_1 \times S_1$ | $\{0, 1\}$ | $\{0\}$ | 0 | — |
| $S_1 \times S_2$ | $\{1/2\}$ | $\{1/2\}$ | 0 | — |
| $S_1 \times S_3$ | $\{1/2, 3/2\}$ | $\{1/2\}$ | 0 | — |
| **$S_1 \times S_4$** | $\{0, 1\}$ | $\{1\}$ | **1** | **$-0.1594$** |
| $S_2 \times S_2$ | $\{0\}$ | $\{0, 1\}$ | 0 | — |
| **$S_2 \times S_3$** | $\{1\}$ | $\{0, 1\}$ | **1** | **$-0.1594$** |
| $S_2 \times S_4$ | $\{1/2\}$ | $\{1/2, 3/2\}$ | 0 | — |
| **$S_3 \times S_3$** | $\{0, 1, 2\}$ | $\{0, 1\}$ | **1** | **$+0.1274$** |
| $S_3 \times S_4$ | $\{1/2, 3/2\}$ | $\{1/2, 3/2\}$ | 0 | — |
| **$S_4 \times S_4$** | $\{0, 1\}$ | $\{0, 1, 2\}$ | **1** | **$+0.1274$** |

**Four off-blocks contain the $(j_L, j_R) = (1, 1)$ graviton irrep**, each with multiplicity 1 (so 9 real Hermitian modes per block in the (1,1) sector). All four blocks have nonzero $S^{(2)}$ eigenvalue.

## 5. Verdict: POSITIVE-FIRST-PASS

**Named falsifier from scoping memo:**

> POSITIVE (G6 viable): At least one eigenmode $V_\alpha$ of $S^{(2)}$ has nonzero $(j_L, j_R) = (1, 1)$ content under $SO(4)$ decomposition, with nonzero eigenvalue $\kappa_\alpha$.

**Result:** four such off-blocks identified, all with nonzero eigenvalue. **POSITIVE** at the first-pass level.

**Implication:** GeoVac's discrete substrate at $n_{\max} = 1$ structurally supports spin-2 perturbations with nonzero quadratic-form eigenvalues. Gravitons are POSSIBLE on the substrate.

This invalidates the scoping memo's deepest worry: I had argued that "spin-1/2 $\otimes$ spin-1/2 only reaches spin-1," but at level $n = 1$ the spinor harmonics are already in $(1, 1/2)$ and $(1/2, 1)$ representations (not pure spin-1/2), so tensor products naturally reach $(1, 1)$. The structural obstacle to gravitons on a spin-1/2 Hilbert space is less severe than the scoping memo anticipated.

## 6. Sign structure observation

The four (1,1)-graviton blocks split into two classes:

- **Within-sector ($S_3 \times S_3$, $S_4 \times S_4$): eigenvalue $+0.127$** — positive, consistent with propagating physical modes
- **Cross-sector ($S_1 \times S_4$, $S_2 \times S_3$): eigenvalue $-0.159$** — negative, potentially indicating instability OR non-physical modes

For graviton dynamics, we want positive kinetic energy (positive $S^{(2)}$ eigenvalue on the propagating modes). The two within-sector blocks are candidates for this. The two cross-sector blocks with negative eigenvalues need physical interpretation:

- If they're gauge modes (non-propagating, projected out by gauge fixing), the negative eigenvalue is harmless
- If they're physical modes with negative kinetic term, the framework has graviton-mode instability at this background
- If they're "fermion-doubling-like" artifacts of the truncated discretization, they may go away in the continuum limit

Distinguishing these requires the multi-month full analysis.

## 7. Honest scope of this first-pass

**Reached:**
- Existence of $(1, 1)$ irrep content in the second-variation quadratic form ✓
- Nonzero eigenvalues for all four graviton-candidate blocks ✓
- Analytical eigenvalues verified by finite difference to rel diff $10^{-9}$ ✓

**Not reached:**
- Verification of Fierz-Pauli kinetic structure (the (1,1) irrep contains 9 modes; physical graviton has 2; the 7-dim difference is gauge structure)
- Gauge invariance / longitudinal mode identification
- Convergence in $n_{\max}$ (test only at $n_{\max} = 1$; need $n_{\max} = 2, 3$ for convergence)
- Sign-structure interpretation (positive vs negative eigenvalue blocks)
- Propagator structure for the $(1, 1)$ modes in the discrete substrate
- Connection to continuum graviton via Paper 38-style propinquity convergence

Full sufficient analysis is the multi-month G6-Diag-Full or G6 sprint per the scoping memo (Path P1 = explicit gamma-matrix re-derivation most likely route).

## 8. Recommendation

The first-pass diagnostic is **POSITIVE** for the named falsifier (graviton-irrep content exists with nonzero eigenvalue). The gravity-arc forward direction is clarified:

- G6 in some form is justified (not structurally blocked)
- The natural next sprint is **G6-Diag-Full** (2-3 weeks): extend to $n_{\max} = 2, 3$, verify convergence, identify gauge modes via Killing-vector analysis, classify the (1, 1) modes into "physical graviton-like" vs "gauge"
- After G6-Diag-Full, if still positive, the full G6 (multi-month, Path P1) would build the explicit Fierz-Pauli derivation

**Alternative:** G4 (Bekenstein-Hawking on cigar) remains a multi-month gravity-side target that doesn't require graviton dynamics. Decision between G6-Diag-Full and G4 is a PI call about which gravity-arc thread to pursue first.

## 9. Cross-references

- **Sprint G3** (`sprint_g3_scalar_TT_S3.md`) — identified that the graviton (TT-tensor) sector has standard CC continuum expansion; G6-Diag now tests whether GeoVac's substrate hosts the discrete-graviton-mode analog
- **Scoping memo** (`debug/g6_scoping_memo.md`) — predicted POSITIVE/NEGATIVE/AMBIGUOUS; this result is POSITIVE-FIRST-PASS with caveats
- **Paper 23 Fock rigidity theorem** — Fock projection fixes the $S^3$ metric; the graviton modes here are perturbations of $D$ at the operator level, not direct metric perturbations
- **External_input_three_class_partition** memory — calibration data lives outside framework scope; if G6 ultimately succeeds, gravitons move INTO framework scope; if NEGATIVE results emerge at later stages, gravitons remain Class 1 external

## 10. Files produced

- `debug/g6_diag_quadratic_form.py` (~330 lines) — driver with analytical formulas, FD verification, $SO(4)$ tensor product decomposition
- `debug/data/g6_diag_quadratic_form.json` — structured numerical results
- `debug/g6_diag_memo.md` — this memo (first-pass canonical)

No CLAUDE.md or Paper updates from this first-pass — substantive paper-level claims require the full G6-Diag-Full sprint to verify Fierz-Pauli kinetic structure.
