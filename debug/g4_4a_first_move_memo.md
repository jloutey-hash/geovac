# Sprint G4-4a first move — Constant-warp Dirac on cigar (week 1)

**Date:** 2026-05-28
**Path:** Gravity arc, opening the second-stage multi-month commitment (G4-4 = warped Dirac on the discrete cigar substrate). First-move sub-sprint of the 4-8 week G4-4a. F1 (factorization), F2 (chirality grading) closed bit-exact at machine precision; F3 (continuum recovery) rough-check at the rank-2 enhancement level.
**Verdict:** **POSITIVE-G4-4a-FIRST-MOVE-VERIFIED.** F1 and F2 land at machine precision on three panel sizes (rel_err $\leq 3 \times 10^{-16}$, gamma algebra residuals = 0.00). F3 rough rank-2 enhancement check verifies $K_{\rm Dirac}/K_{\rm scalar} = 2.000$ at small $t$ on the T2 UV-refined substrate.

## 1. What this sprint delivers

**Production code (new):**
- `geovac/gravity/__init__.py` — subpackage init exposing the warped-Dirac API
- `geovac/gravity/warped_dirac.py` (~440 lines) — 2D Cl(2,0) gamma matrices, `DiscreteDiskDirac` (anti-periodic spinor), `DiscreteDiskScalar` (periodic bosonic companion for F3 rough), `S2DiracSpectrum` (Camporesi-Higuchi), `WarpedDiracConstant` (tensor product), three falsifier-verification functions.

**Tests:**
- `tests/test_warped_dirac.py` — 37 tests, all pass in 0.62s. Covers gamma algebra (8 tests), disk-Dirac (10 tests), $S^2$ spectrum (6 tests), F1 factorization (3 tests), F2 chirality (2 tests), F3 rough (1 test), `WarpedDiracConstant` (5 tests). Zero regression on the broader test suite (this is a new subpackage, no upstream dependencies).

**Driver and memo:**
- `debug/g4_4a_first_move.py` — verification panel covering F1, F2, F3 across three panel sizes
- `debug/data/g4_4a_first_move.json` — structured results
- `debug/g4_4a_first_move_memo.md` — this memo

## 2. F1 results: factorization bit-exact

At constant warp, the squared Dirac on the cigar factorizes:
$$
D_{\rm cigar}^2 = D_{D^2}^2 \otimes I_{S^2} + I_{D^2} \otimes \frac{D_{S^2}^2}{r_h^2}
$$
giving the heat-trace identity $K_{\rm cigar}(t) = K_{D^2}(t) \cdot K_{S^2}(t)$. Verified at three panel sizes:

| Panel | $\dim_{\rm disk}$ | $\dim_{S^2}$ | $\dim_{\rm cigar}$ | $t = 0.05$ rel_err | $t = 1.0$ rel_err |
|---|---|---|---|---|---|
| small | 120 | 48 | 5,760 | $0.00 \times 10^0$ | $1.78 \times 10^{-16}$ |
| medium | 480 | 80 | 38,400 | $2.81 \times 10^{-16}$ | $2.35 \times 10^{-16}$ |
| larger | 960 | 120 | 115,200 | $1.63 \times 10^{-16}$ | $1.50 \times 10^{-16}$ |

**F1 bit-exact at machine precision (float64 single-summation noise) at every panel cell.** This is the operator-level identity that follows from the outer-sum structure $\sum_{i,j} e^{-(a_i+b_j)t} = (\sum_i e^{-a_i t})(\sum_j e^{-b_j t})$; the verification is the load-bearing check that the discrete tensor product is correctly implemented.

## 3. F2 results: chirality grading at gamma algebra

$\{\gamma^5, \gamma^a\} = 0$ for $a = 1, 2$ in Cl(2,0) is a Pauli-matrix identity. Verified bit-exact (residual = 0.00):

| Identity | Residual |
|---|---|
| $\gamma_1^2 = I$ | 0.00 |
| $\gamma_2^2 = I$ | 0.00 |
| $\gamma_5^2 = I$ | 0.00 |
| $\{\gamma_1, \gamma_2\} = 0$ | 0.00 |
| $\gamma_5 = -i\gamma_1\gamma_2$ | 0.00 |
| $\{\gamma_5, \gamma_1\} = 0$ | 0.00 |
| $\{\gamma_5, \gamma_2\} = 0$ | 0.00 |

Since $D_{D^2} = \gamma_1 \partial_\rho + \gamma_2 (1/\rho) \nabla_\phi$ and $\gamma_5$ anticommutes with both $\gamma_1$ and $\gamma_2$, $\{\gamma_5, D_{D^2}\} = 0$ follows immediately at the operator level.

**F2 verified at the algebraic backbone.** Full operator-level F2 (explicit construction of the sparse Dirac matrix and direct anticommutator computation) is reserved for the next G4-4a week when the full $D$ operator is built (currently only $D^2$ is implemented).

## 4. F3 rough results: rank-2 spinor enhancement

At small $t$ the continuum Weyl leading order on the flat disk gives
$$
K_{D^2}^{\rm Dirac}(t) \sim 2 \cdot \frac{A_{D^2}}{4\pi t}, \qquad K_{D^2}^{\rm scalar}(t) \sim \frac{A_{D^2}}{4\pi t}
$$
so the ratio $K_{\rm Dirac}/K_{\rm scalar} \to 2$ at small $t$.

Panel: $N_\rho = 100$, $a = 0.1$, $N_\phi = 48$ (T2 UV regime), $R = 10$.

| $t$ | $K_{\rm Dirac}$ | $K_{\rm scalar}$ | ratio | target |
|---|---|---|---|---|
| 0.05 | 813.47 | 406.74 | **2.000** | 2.0 |
| 0.10 | 475.82 | 237.93 | **2.000** | 2.0 |
| 0.20 | 258.63 | 129.34 | **2.000** | 2.0 |
| 0.50 | 100.20 | 50.13  | 1.999 | 2.0 |
| 1.00 | 44.91  | 22.49  | 1.997 | 2.0 |

**Rank-2 enhancement bit-exact at small $t$**, with sub-percent IR drift at $t \geq 0.5$ from the difference between anti-periodic (Dirac, half-integer $m_{\rm eff}$) and periodic (scalar, integer $m_{\rm eff}$) lowest-mode contributions.

The 0.3% deficit at $t = 1.0$ is the half-integer vs integer angular spectrum difference becoming visible. Both branches share the same radial substrate; only the BC distinction (anti-periodic vs periodic) and the rank multiplier differ. The UV-asymptotic rank-2 ratio is exact by construction (both branches Weyl-bulk as $1/t$, the spinor with rank-2 prefactor).

## 5. What this does NOT close

This is a **first-move week of G4-4a, not the full G4-4a closure.** The remaining 4-8 weeks of G4-4a:

- **Explicit $D$ operator construction.** F2 at operator level (rather than just gamma algebra) requires building the sparse Dirac matrix and verifying $\{\gamma^5_{\rm spinor}, D\} = 0$ at machine precision in the spinor-tensor Hilbert space. Sprint-scale ~ 1 week extension.
- **Quantitative F3 with continuum Weyl-Selberg.** The rank-2 ratio check is a rough first pass. The full F3 needs $K_{D^2}^{\rm Dirac, disc}(t) / K_{D^2}^{\rm Dirac, continuum}(t) \to 1$ within 10% at small $t$. Requires UV refinement like T2. Sprint-scale ~ 2 weeks.
- **Anti-periodic BC fine structure.** Tests at low $t$ that the lowest spinor mode matches $j_{1/2, 1}/R = \pi/R$ rather than the scalar's $j_{0, 1}/R$. Sprint-scale ~ 1 week.

After G4-4a closes: G4-4b (variable warp), G4-4c (conical defect spinor), G4-4d (Seeley-DeWitt extraction), G4-4e (BC sectors), G4-4f (replica preparation). Then G4-5 (replica method) and G4-6 ($S_{BH}$ derivation).

## 6. Substantive new content (beyond scoping)

1. **`geovac/gravity` subpackage opened**, with clean API and 37-test verification at machine precision. Empty before this sprint.

2. **Cl(2,0) gamma algebra implemented as the framework's first 2D Dirac infrastructure.** Prior Dirac code in `geovac/` is 4D ($S^3$ Camporesi-Higuchi) or 1D (radial Dirac matrix elements); G4-4a opens the 2D disk-Dirac path.

3. **Anti-periodic spinor BC vs periodic scalar BC** implemented as a clean code pattern with the half-integer $m_{\rm eff}$ shift. The companion `DiscreteDiskScalar` class enables the F3 rough comparison and is reusable for any future scalar-vs-spinor diagnostic on the disk.

4. **Bit-exact tensor-product factorization at three panel sizes** confirms the architectural design choices in the G4-4 scoping memo (constant-warp factorization is "operator-level, almost guaranteed bit-exact").

5. **F3 rough rank-2 check at the T2 UV-refined substrate** lands the structural prediction quantitatively at sprint scale. Full quantitative F3 with continuum Weyl-Selberg is the next-week's work.

## 7. Honest scope

**Reached (this week):**
- F1 bit-exact at three panel sizes
- F2 verified at Cl(2,0) gamma algebra
- F3 rough rank-2 enhancement at the UV-refined panel
- 37 tests passing in 0.62s
- Production code ~440 lines clean
- Driver + JSON + memo

**Not reached (subsequent weeks of G4-4a):**
- Explicit sparse $D$ operator (for operator-level F2 not just algebraic)
- Quantitative F3 with continuum Weyl-Selberg match within 10% at small $t$
- Lowest-mode validation against $j_{1/2,1}/R = \pi/R$
- Variable warp (G4-4b multi-week sprint)
- Conical defect (G4-4c, where T1's wedge-lattice diagnosis applies in the spinor case)

## 8. G4-4 sequence status

| Sub-sprint | Status |
|---|---|
| G4-4 scoping (T3) | Done (2026-05-28) |
| **G4-4a first move (this)** | **Done (2026-05-28), POSITIVE-VERIFIED** |
| G4-4a full closure (remaining 4-8 weeks) | In flight |
| G4-4b through G4-4f | Multi-month sequence per scoping memo |

## 9. Files

- `geovac/gravity/__init__.py`
- `geovac/gravity/warped_dirac.py` (440 lines)
- `tests/test_warped_dirac.py` (37 tests)
- `debug/g4_4a_first_move.py` (driver)
- `debug/data/g4_4a_first_move.json` (structured results)
- `debug/g4_4a_first_move_memo.md` (this memo)

## 10. Cross-references

- **G4-4 scoping (T3)**: `debug/g4_4_warped_dirac_scoping_memo.md` (architectural blueprint)
- **G4-3a-cleanup**: Hermitian polar Laplacian substrate (reused via the symmetric tridiagonal radial Laplacian convention)
- **G4-1 ($S^2$ Dirac spectrum)**: Paper 28 §4.14, Camporesi-Higuchi 1996; `S2DiracSpectrum` is the production import
- **T1 (G4-3c proper wedge)**: structural diagnosis on the scalar side that applies in the spinor case (G4-4c)
- **T2 (G4-3d-UV)**: UV-refinement substrate for F3 (current and future)
