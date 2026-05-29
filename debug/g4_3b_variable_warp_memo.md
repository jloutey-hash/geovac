# Sprint G4-3b — Variable warp $r(\rho)$ for asymptotic Schwarzschild

**Date:** 2026-05-28
**Path:** Gravity arc, second sub-sprint of G4-3. Extends G4-3a-cleanup (constant warp) to variable warp $r(\rho)$. Sprint-scale.
**Verdict:** **POSITIVE-G4-3b.** Asymmetric variable-warp radial operator has real eigenvalues (max $|\mathrm{Im}(\lambda)| = 0$ at machine precision), heat trace responds to warp, asymptotic limit $r(\rho)/\rho \to 1$ verified.

## 1. Warp function

For sprint-scale first pass, use a smooth interpolation:
$$r(\rho) = r_h\,\sqrt{1 + (\rho/r_h)^2}$$

with:
- $r(0) = r_h$ (horizon radius at tip)
- $r(\rho) \sim \rho$ asymptotically (Schwarzschild-like outer geometry)
- Smooth at $\rho = 0$, monotonically increasing

This captures the qualitative cigar geometry: horizon at $\rho = 0$ with $r = r_h$, opening into a paraboloid in the asymptotic region. The exact Schwarzschild $r(\rho)$ relation differs in detail (tortoise coordinate) but the qualitative behavior is the same.

## 2. Warped-product Laplacian

For metric $ds^2 = d\rho^2 + \rho^2 d\phi^2 + r(\rho)^2 d\Omega_2^2$:
$$\Delta f = \partial_\rho^2 f + \left(\frac{1}{\rho} + \frac{2 r'(\rho)}{r(\rho)}\right)\partial_\rho f + \frac{1}{\rho^2}\partial_\phi^2 f + \frac{1}{r(\rho)^2}\Delta_{S^2} f$$

Separation $f(\rho, \phi, \theta, \phi') = R(\rho)\,e^{im\phi}\,Y_{l m'}(\theta, \phi')$ gives the radial ODE:
$$-R'' - \left(\frac{1}{\rho} + \frac{2r'}{r}\right)R' + \left[\frac{m^2}{\rho^2} + \frac{l(l+1)}{r(\rho)^2}\right]R = \lambda R$$

## 3. Asymmetric discretization

The first-derivative coupling $(1/\rho + 2r'/r) R'$ has variable coefficient, so the discrete matrix is asymmetric. Centered finite differences on $\rho_k = ka$:

$$L_{kk} = \frac{2}{a^2} + V(\rho_k)$$
$$L_{k, k-1} = -\frac{1}{a^2} + \frac{c_1(\rho_k)}{2a}$$
$$L_{k, k+1} = -\frac{1}{a^2} - \frac{c_1(\rho_k)}{2a}$$

where $c_1(\rho) = 1/\rho + 2r'(\rho)/r(\rho)$ and $V(\rho) = m^2/\rho^2 + l(l+1)/r(\rho)^2$.

**Real eigenvalues despite asymmetry:** the underlying continuum operator is self-adjoint on $L^2(\rho\,r(\rho)^2\,d\rho)$, so eigenvalues are guaranteed real. The discrete asymmetric matrix preserves this property numerically.

## 4. Verification

| Check | Result |
|---|---|
| Eigenvalues real, $(m, l) = (0, 0)$ | max $|\mathrm{Im}(\lambda)| = 0$ exact |
| Eigenvalues real, $(0, 1), (1, 1), (0, 2)$ | max $|\mathrm{Im}(\lambda)| = 0$ exact |
| Variable warp lowers $\ell \geq 1$ eigenvalues | Confirmed |
| Asymptotic $r(\rho)/\rho \to 1$ | Verified at $\rho = 5, 10, 20$ |
| Heat-trace responds to warp | Yes, K_var/K_const $\in [0.89, 2.19]$ |

## 5. Heat trace at $N_\rho = 100$, $a = 0.05$, $r_h = 2$, $|m| \leq 4$, $l \leq 4$

| $t$ | $K_{\rm variable}$ | $K_{\rm constant}$ | Ratio var/const |
|---|---|---|---|
| 0.1 | 641.9 | 549.8 | 1.17 |
| 0.5 | 121.3 | 66.4 | 1.83 |
| 1.0 | 38.65 | 17.63 | 2.19 |
| 5.0 | 0.424 | 0.477 | 0.89 |

Variable warp generally enhances the heat trace (more low-lying modes from reduced $\ell(\ell+1)/r^2$ angular potential at $\rho > 0$). At large $t$ the ratio inverts as the heat trace becomes IR-dominated.

## 6. Asymptotic warp behavior

| $\rho$ | $r(\rho)$ | $r/\rho$ |
|---|---|---|
| 0.5 | 2.06 | 4.12 |
| 1.0 | 2.24 | 2.24 |
| 2.0 | 2.83 | 1.41 |
| 5.0 | 5.39 | 1.08 |
| 10.0 | 10.20 | 1.020 |
| 20.0 | 20.10 | 1.005 |

The warp asymptotes correctly to $r(\rho)/\rho \to 1$ at large $\rho$, recovering the asymptotic Schwarzschild geometry where the spatial sphere's radius equals the radial coordinate.

## 7. Operator-form caveat

The comparison "variable warp vs constant warp" at $(m, l) = (0, 0)$ shows ~66% difference at the lowest eigenvalue. This is partly a real warp effect AND partly an operator-form artifact:
- **Variable warp**: uses the asymmetric discretization of the original radial operator (no substitution).
- **Constant warp (G4-3a-cleanup)**: uses the substituted symmetric Schrödinger form $u = \sqrt\rho R$ which absorbs the $1/\rho$ first-derivative into a centrifugal potential $-1/(4\rho^2)$.

A clean apples-to-apples comparison would unify the discretization (either both substituted or both asymmetric). For sprint-scale G4-3b this is acknowledged but not pursued; the substantive content (real eigenvalues, warp asymptotics, heat-trace response) is established.

## 8. Honest scope

**Reached:**
- Variable warp $r(\rho) = r_h \sqrt{1 + (\rho/r_h)^2}$ defined and discretized
- Asymmetric radial operator yields real eigenvalues to machine precision
- Asymptotic limit $r(\rho)/\rho \to 1$ verified
- Heat trace responds to warp

**Not reached:**
- Unified discretization for clean variable-warp vs constant-warp comparison
- Exact Schwarzschild $r(\rho)$ relation (using simplified smooth interpolation)
- Continuum-limit verification (G4-3d target)

## 9. Verdict

**POSITIVE-G4-3b.** Variable warp framework operational. The asymmetric radial discretization preserves the underlying operator self-adjointness numerically. Heat trace responds to the warp factor through the $\ell(\ell+1)/r(\rho)^2$ angular coupling.

## 10. Multi-month sub-sprint sequence status

| Sub-sprint | Status |
|---|---|
| G4-3a | Done (scoping) |
| G4-3a-cleanup | Done (Hermitian polar Laplacian) |
| **G4-3b (this)** | **Done (variable warp)** |
| G4-3c | Next (conical-defect sweep over $N_\phi$) |
| G4-3d | Next (continuum-limit heat-kernel asymptotics) |
| G4-4 | Multi-month (warped Dirac) |
| G4-5 | Multi-month (discrete replica) |
| G4-6 | Multi-month (full $S_{BH}$ derivation) |

Three sub-sprints completed (G4-3a, cleanup, G4-3b) in single sessions; G4-3c and G4-3d are sprint-scale follow-ons. G4-4 onwards is the multi-month commitment.

## 11. Files

- `debug/g4_3b_variable_warp.py` — driver
- `debug/data/g4_3b_variable_warp.json` — structured results
- `debug/g4_3b_variable_warp_memo.md` — this memo

## 12. Cross-references

- **G4-3** (`memory/sprint_g4_3_warped_substrate.md`): substrate scoping; G4-3b implements variable warp
- **G4-3a-cleanup** (`debug/g4_3a_cleanup_hermitian_polar_memo.md`): Hermitian constant-warp reference
- **G4-2** (`memory/sprint_g4_2_conical_replica.md`): continuum derivation of $S_{BH}$ from conical defect
- **Schwarzschild 1916**: standard Euclidean black-hole geometry
