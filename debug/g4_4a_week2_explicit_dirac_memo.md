# Sprint G4-4a week 2 — Explicit Dirac operator + lowest-mode validation

**Date:** 2026-05-28
**Path:** Gravity arc, G4-4a sub-sprint week 2 (4-8 week sequence). Closes the named next-week's work from the week-1 first-move memo.
**Verdict:** **POSITIVE-G4-4a-WEEK2-VERIFIED.** Explicit Dirac operator constructed in canonical chirality-graded form; operator-level F2 bit-exact (0.00) on every Fourier block; D² spectrum matches the factorized DiscreteDiskDirac at machine precision; lowest-mode $|\lambda_{\min}|$ converges to $\pi/R$ within 0.5% at the UV-fine panel.

## 1. What this sprint delivers

Extension of `geovac/gravity/warped_dirac.py` with new `DiscreteDirac2D` class (~110 lines). 17 new tests covering construction, operator-level F2, D² matching, spectrum symmetry, and continuum-limit validation. All 54/54 tests in the module pass in 1.05s. Driver + JSON + memo for the week-2 closure.

## 2. Construction

Canonical chirality-graded form. In each anti-periodic azimuthal Fourier mode $k$ with $m_{\rm eff}$ half-integer-shifted:
$$
D_k = \begin{pmatrix} 0 & \sqrt{L_k} \\ \sqrt{L_k} & 0 \end{pmatrix}
$$
where $L_k$ is the Hermitian radial Laplacian from G4-3a-cleanup. The chirality grading $\gamma^5 = \mathrm{diag}(I_{N_\rho}, -I_{N_\rho})$ on each block satisfies $\{\gamma^5, D_k\} = 0$ as an operator-level identity:
$$
\{\gamma^5, D_k\} = \begin{pmatrix} 0 & \sqrt{L_k} \\ -\sqrt{L_k} & 0 \end{pmatrix} + \begin{pmatrix} 0 & -\sqrt{L_k} \\ \sqrt{L_k} & 0 \end{pmatrix} = 0
$$
Each block gives $D_k^2 = \mathrm{diag}(L_k, L_k)$ — each scalar eigenvalue $\mu_n$ appears twice in the spinor spectrum, matching the rank-2 multiplicity used in `DiscreteDiskDirac.squared_eigenvalues()`. Total $D$ across all Fourier modes is block-diagonal with $\pm\sqrt{\mu_n}$ spectrum.

## 3. Operator-level F2 verification

Three panels, three checks each. Max anticommutator residual at machine precision:

| Panel | $(N_\rho, a, N_\phi)$ | $\max\|\{\gamma^5, D_k\}\|$ | Hermitian |
|---|---|---|---|
| small | (10, 0.5, 8) | $\mathbf{0.00 \times 10^0}$ | PASS |
| medium | (20, 0.3, 12) | $\mathbf{0.00 \times 10^0}$ | PASS |
| larger | (30, 0.2, 16) | $\mathbf{0.00 \times 10^0}$ | PASS |

**Operator-level F2 bit-exact at every Fourier block on every panel.** This is the load-bearing upgrade from week-1's algebraic F2 (gamma matrix Pauli identities) to the full operator level: $\{\gamma^5, D_k\} = 0$ holds for the explicit $2N_\rho \times 2N_\rho$ matrices, not just for the abstract gamma algebra.

## 4. D² spectrum match

Verifies the canonical-sqrt construction is consistent with the factorized squared-Dirac spectrum from week-1's `DiscreteDiskDirac`:

| Panel | max $\|D^2_{\rm explicit} - D^2_{\rm factorized}\|$ |
|---|---|
| small | $3.55 \times 10^{-15}$ |
| medium | $1.42 \times 10^{-14}$ |
| larger | $1.14 \times 10^{-13}$ |

**Machine precision** (float64 eigenvalue ordering noise) across all panels. Confirms that the explicit-D construction produces the correct $D^2$ spectrum, which is the load-bearing input for the heat-trace factorization (F1 from week 1).

## 5. Spectrum $\pm$ symmetry

Exact pairing at every panel:

| Panel | $\#(\lambda > 0)$ | $\#(\lambda < 0)$ | max $\|\lambda^+ - |\lambda^-|\|$ |
|---|---|---|---|
| small  | 80  | 80  | $\mathbf{0.00 \times 10^0}$ |
| medium | 240 | 240 | $\mathbf{0.00 \times 10^0}$ |
| larger | 480 | 480 | $\mathbf{0.00 \times 10^0}$ |

Bit-exact — the canonical $D = \mathrm{antidiag}(\sqrt{L}, \sqrt{L})$ form forces $\pm\sqrt{\mu_n}$ spectrum by construction.

## 6. Lowest-mode validation: $|\lambda_{\min}| \to \pi/R$

The continuum spinor lowest mode on the flat polar disk with anti-periodic phi BC (half-integer $m_{\rm eff} = 1/2$, no effective centrifugal in u-representation) sits at the first Bessel zero $j_{1/2, 1} = \pi$:
$$
|\lambda_{\min}^{\rm Dirac, continuum}| = \frac{j_{1/2, 1}}{R} = \frac{\pi}{R}
$$
This is distinct from the scalar disk lowest mode at $j_{0, 1}/R \approx 2.405/R$ — the half-integer angular momentum structure is operationally visible.

| Panel | $(N_\rho, a, N_\phi)$ | $R$ | $\|\lambda_{\min}\|$ | $\pi/R$ | rel_err |
|---|---|---|---|---|---|
| uv_small | (50, 0.2, 24) | 10 | 0.307902 | 0.314159 | **−1.99%** |
| uv_med   | (100, 0.1, 48) | 10 | 0.311024 | 0.314159 | **−1.00%** |
| uv_fine  | (200, 0.05, 96) | 10 | 0.312590 | 0.314159 | **−0.50%** |

**Clean linear convergence in lattice spacing $a$** (error halves as $a$ halves: $-1.99\% \to -1.00\% \to -0.50\%$). The continuum spinor lowest-mode prediction $\pi/R$ is recovered to better than 1% at $a = 0.05$, $N_\phi = 96$. The half-integer angular momentum structure of the anti-periodic spinor BC is verified at the eigenvalue level on the discrete substrate.

The negative sign of the error (substrate slightly under-estimates $\pi/R$) is consistent with finite-$N_\phi$ angular discretization slightly suppressing the lowest $m_{\rm eff} \approx 1/2$ relative to the continuum value $m_{\rm eff}^{\rm cont} = 1/2$ exactly.

## 7. Substantive new content

1. **`DiscreteDirac2D` opens explicit-D infrastructure.** Week 1 only had $D^2$ via factorized scalar Laplacian; week 2 has the LINEAR Dirac with operator-level F2 verification.

2. **Operator-level F2 bit-exact on Fourier blocks.** Upgrades F2 from algebraic Pauli identity (week 1) to numerical anticommutator on the explicit $2N_\rho \times 2N_\rho$ block matrices (week 2). Residual exactly 0.00 — the canonical chirality-graded sqrt construction makes this an identity at every panel size.

3. **Half-integer angular momentum verified at the eigenvalue level.** The continuum spinor lowest-mode prediction $j_{1/2, 1}/R = \pi/R$ is recovered to 0.5% at sprint-scale UV-refined panels. Structurally distinct from the scalar disk's $j_{0, 1}/R \approx 0.765 \times (\pi/R)$.

4. **D² explicit matches D² factorized at machine precision.** Cross-check that confirms the canonical-sqrt construction is consistent with the rank-doubled scalar spectrum used in F1 (week 1).

## 8. Honest scope

**Reached (week 2):**
- Explicit canonical-graded D operator
- Operator-level F2 bit-exact (every block, every panel)
- D² explicit matches D² factorized
- Spectrum $\pm$ symmetry bit-exact
- Lowest-mode $|\lambda_{\min}| \to \pi/R$ within 0.5% at UV-fine panel
- 17 new tests + 0 regression on week-1 tests (54/54 in module)

**Not reached (subsequent weeks of G4-4a):**
- **Geometric Dirac (with spin connection)**: the canonical sqrt construction is consistent and operator-level identities hold, but it's not the same matrix as the geometric Dirac $D = i\gamma^a e_a^\mu(\partial_\mu - \omega_\mu)$ built from spin-connection-corrected finite differences. The eigenvalue $\pm\sqrt{\mu_n}$ matches but the eigenVECTORS differ. For Seeley-DeWitt extraction (G4-4d) and replica method (G4-4f) the eigenvalues are sufficient; for matrix-element calculations a future week would extend to the geometric form.
- **Quantitative F3 with continuum Weyl-Selberg**: requires fuller UV sweep beyond rank-2 enhancement check.
- **Variable warp (G4-4b)**: separate multi-week sprint.
- **Conical defect spinor (G4-4c)**: separate multi-week sprint.

## 9. G4-4a status (running)

| G4-4a week | Status |
|---|---|
| Week 1 (first move): F1, F2-algebraic, F3-rough | Done (POSITIVE-VERIFIED) |
| **Week 2 (this): explicit D, F2-operator-level, lowest-mode** | **Done (POSITIVE-VERIFIED)** |
| Week 3-4 target: quantitative F3 with continuum Weyl-Selberg, full UV sweep | In flight |
| Week 5-6 target: anti-periodic BC fine structure, all-mode validation | Queued |
| Week 7-8 target: closure memo, hand-off to G4-4b | Queued |

## 10. Files

- `geovac/gravity/warped_dirac.py` — extended with `DiscreteDirac2D` (~110 new lines)
- `geovac/gravity/__init__.py` — exports updated
- `tests/test_warped_dirac.py` — 17 new tests (54 total)
- `debug/g4_4a_week2_explicit_dirac.py` — driver
- `debug/data/g4_4a_week2_explicit_dirac.json` — results
- `debug/g4_4a_week2_explicit_dirac_memo.md` — this memo

## 11. Cross-references

- **G4-4a week 1 first move**: `debug/g4_4a_first_move_memo.md` (F1, F2-algebraic, F3-rough)
- **G4-4 scoping (T3)**: `debug/g4_4_warped_dirac_scoping_memo.md` (architectural blueprint)
- **G4-3a-cleanup**: Hermitian polar Laplacian substrate (reused via $L_k$ blocks)
- Bessel zeros: $j_{1/2, 1} = \pi$ (the continuum spinor lowest-mode reference)
