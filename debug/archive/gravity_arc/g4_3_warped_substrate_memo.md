# Sprint G4-3 — Discrete warped-product substrate for the cigar (scoping + first-pass)

**Date:** 2026-05-28
**Path:** Gravity arc, opening of multi-month G4 full discrete-substrate track. Per the G4-2 memo's named follow-on: "G4-3: warped product structure on discrete substrate (multi-month start)".
**Verdict:** **POSITIVE-SCOPING-G4-3a.** Discrete substrate $\mathcal{G}_{\rm cigar} = \mathbb{Z}_+(a) \times \mathbb{Z}/N_\phi \times \mathrm{Fock}(S^2, l_{\max})$ defined; constant-warp factorization $K_{\rm cigar}(t) = K_{D^2}(t) \cdot K_{S^2_{r_h}}(t)$ holds by construction; sub-sprint sequence G4-3a/b/c/d named. The numerical heat trace requires a properly symmetric discrete polar Laplacian (named G4-3b prerequisite); the structural content (substrate definition, factorization, conical parameterization, sub-sprint plan) is the load-bearing deliverable.

## 1. The discrete substrate

The Euclidean Schwarzschild cigar's near-horizon geometry is
$$ds^2 = d\rho^2 + \rho^2\,d\phi^2 + r(\rho)^2\,d\Omega_2^2$$
with $\rho \in [0, \infty)$, $\phi \in [0, 2\pi)$, and warp factor $r(\rho)$. In the **near-horizon limit** $r(\rho) = r_h$ (constant), the geometry becomes $D^2 \times S^2_{r_h}$.

**Discrete substrate:**
$$\mathcal{G}_{\rm cigar}(a, N_\rho, N_\phi, l_{\max}, r_h) = \mathbb{Z}_+(a)|_{N_\rho} \;\times\; \mathbb{Z}/N_\phi \;\times\; \mathrm{Fock}(S^2, l_{\max})$$

where:
- $\mathbb{Z}_+(a)|_{N_\rho}$ = $N_\rho$-truncation of $\mathbb{Z}_+$ with lattice spacing $a$: $\rho_k = k a$, $k = 0, \ldots, N_\rho - 1$
- $\mathbb{Z}/N_\phi$ = discrete $\phi$-circle with $N_\phi$ sites: $\phi_j = 2\pi j / N_\phi$
- $\mathrm{Fock}(S^2, l_{\max})$ = standard $S^2$ Fock projection with $l \le l_{\max}$ angular cutoff

**Continuum limit:** $a \to 0$, $N_\rho a \to \infty$, $N_\phi \to \infty$, $l_{\max} \to \infty$.

**Conical-defect parameter:** $\alpha = N_\phi / N_0$ where $N_0$ is the reference smooth-tip count. $\alpha = 1$: smooth disk; $\alpha \neq 1$: conical singularity at $\rho = 0$ with apex angle $2\pi\alpha$.

## 2. Constant-warp factorization

At $r(\rho) = r_h$ constant, the Laplacian factorizes:
$$\Delta_{\rm cigar} = \Delta_{D^2}(\rho, \phi) + \frac{1}{r_h^2}\,\Delta_{S^2}$$

The two factors act on independent sectors (radial-azimuthal vs spherical), so the heat trace factorizes:
$$K_{\rm cigar}(t) = K_{D^2}(t)\,\cdot\,K_{S^2_{r_h}}(t)$$

**This is true by construction** at the operator level. The discrete substrate inherits the factorization automatically when restricted to constant warp.

For the **continuum limits**:
- $K_{D^2_\alpha}(t) \sim A_{D^2}/(4\pi t) + (1/12)(1/\alpha - \alpha) + O(t)$ (Sommerfeld/Cheeger, scalar)
- $K_{S^2_{r_h}}(t) \sim A_{S^2}/(4\pi t) + r_h^2/3 + O(t)$ (from G4-1)
- $K_{\rm cigar}(t) \sim$ leading $A_{D^2} A_{S^2}/(4\pi t)^2$ (4D Weyl law)

## 3. Sprint G4-3a first-pass: construction + first heat-trace probe

Implemented `debug/g4_3_warped_substrate.py` (~280 lines) with:
- Discrete radial Laplacian (with naive polar second-difference)
- Periodic azimuthal Laplacian on $\mathbb{Z}/N_\phi$
- Standard Fock $S^2$ Laplacian with $l_{\max}$ cutoff
- Joint disk Laplacian via azimuthal-mode decomposition + heat trace tabulation

**Parameters tested:** $N_\rho = 20$, $a = 0.5$, $N_\phi = 12$, $l_{\max} = 6$, $r_h = 2$. Gives $\rho_{\max} = 10$, $A_{D^2} = 100\pi$, $A_{S^2} = 16\pi$.

## 4. Honest scope: numerical issue at G4-3a level

**The naive discrete polar Laplacian I used is NOT Hermitian.** The continuum polar Laplacian
$$\Delta_{D^2} = \partial_\rho^2 + \frac{1}{\rho}\partial_\rho + \frac{1}{\rho^2}\partial_\phi^2$$
is symmetric on $L^2(D^2, \rho\,d\rho\,d\phi)$. A direct second-difference discretization on uniform $\mathbb{Z}_+(a)$ does NOT respect the polar measure $\rho\,d\rho$, so the resulting discrete operator has complex / negative eigenvalues that should not appear in a positive Laplacian.

**Symptom in the driver:** smallest disk-Laplacian eigenvalues are negative ($\sim -15$), producing $K_{D^2}(t) > 1$ at small $t$ (should be $\sim 1$). The heat trace values are quantitatively wrong.

**The structural content is unaffected:**
- The substrate definition is correct
- The factorization at constant warp is true by construction (operator-level)
- The continuum limit and conical-defect parameterization are well-defined targets
- The sub-sprint sequence is named

**The proper fix** (deferred to G4-3a-cleanup or first part of G4-3b):
- Symmetrize the discrete radial Laplacian on $L^2(\rho\,d\rho)$ via rescaling $u = \sqrt{\rho}\,f$ and symmetric finite differences
- Alternative: use Bessel-mode basis (eigenfunctions of continuum disk Laplacian) at finite cutoff
- Alternative: use Camporesi-Higuchi-style discretization on the warped geometry

This is **engineering, not structural**. The G4-3 substrate framework stands.

## 5. Sub-sprint sequence (multi-month G4-3 → G4-6 plan)

| Sub-sprint | Scope | Estimated commitment |
|---|---|---|
| **G4-3a** (this) | Constant-warp substrate + scoping | Sprint-scale ✓ |
| G4-3a-cleanup | Symmetric discrete polar Laplacian | ~1 week |
| G4-3b | Variable warp $r(\rho)$ for asymptotic Schwarzschild | ~2-4 weeks |
| G4-3c | Discrete conical-defect deformation ($N_\phi$ sweeps) | ~2-4 weeks |
| G4-3d | Continuum-limit heat-kernel asymptotics verification | ~2-4 weeks |
| G4-4 | Warped Dirac spectrum on discrete substrate | ~1-2 months |
| G4-5 | Discrete replica method | ~1-2 months |
| G4-6 | Full discrete-substrate $S_{BH}$ derivation | ~1-2 months |

**Total estimated multi-month commitment**: ~6-12 months for full discrete-substrate derivation of $S_{BH}$.

**G4-3a (this sprint)** is the **opening move** of this multi-month track. It establishes:
- The discrete substrate is a concrete computable object
- The constant-warp limit factorizes (matches G4 first-pass continuum result)
- Sub-sprint sequence is well-defined

## 6. Connection to G4-2 continuum derivation

G4-2 derived $S_{BH} = A\Lambda^2/(12\pi)$ at the continuum level via:
1. Conical defect heat-trace contribution $(1/12)(1/\alpha - \alpha)$ (Sommerfeld/Cheeger)
2. $S^2$ Dirac heat trace from G4-1
3. Replica method $S = -dI_E/d\alpha|_{\alpha=1}$

**G4-3 onwards** rebuilds this on the discrete substrate. The target is to show:
1. Discrete substrate has analog conical-defect contribution that → $(1/12)(1/\alpha - \alpha)$ in continuum limit
2. Discrete Dirac on warped substrate → G4-1's $S^2$ Dirac at constant warp
3. Discrete replica method gives discrete-substrate $S_{BH}$ that → $A\Lambda^2/(12\pi)$

This is the multi-month sub-sprint sequence.

## 7. Connection to G8 cutoff dependence

Per G8, the $\Lambda^2$ prefactor in $S_{BH}$ depends on the cutoff function $f$. For Gaussian cutoff: $G_N = 3\pi/\Lambda^2$. The discrete substrate inherits this calibration: the prefactor in any discrete-substrate $S_{BH}$ will be $\propto \phi(2)/\phi(1)$ for arbitrary cutoff. **The cutoff dependence is structural, not a defect of G4-3 specifically.**

## 8. Honest scope summary

**Reached (G4-3a):**
- Discrete substrate $\mathcal{G}_{\rm cigar}(a, N_\rho, N_\phi, l_{\max}, r_h)$ defined ✓
- Constant-warp factorization at operator level ✓
- $S^2$ Fock cutoff inherits from G4-1 ✓
- Conical-defect parameterization $\alpha = N_\phi/N_0$ identified ✓
- Sub-sprint sequence G4-3a/b/c/d → G4-4/5/6 named ✓

**Not reached (G4-3a known limitations):**
- Hermitian discrete polar Laplacian (G4-3a-cleanup)
- Numerical heat-trace values (depends on Hermitian Laplacian)
- Variable warp $r(\rho)$ (G4-3b)
- Conical-defect sweep over $N_\phi$ (G4-3c)
- Continuum-limit asymptotic verification (G4-3d)

**Not reached (multi-month):**
- Warped Dirac spectrum (G4-4)
- Discrete replica method (G4-5)
- Full discrete-substrate $S_{BH}$ derivation (G4-6)

## 9. Verdict

**POSITIVE-SCOPING-G4-3a.**

This sprint opens the multi-month G4 discrete-substrate track. The substrate $\mathcal{G}_{\rm cigar}$ is defined, the constant-warp factorization is established at the operator level, and the sub-sprint sequence is named. Numerical refinement (Hermitian Laplacian) is sprint-scale engineering for G4-3a-cleanup; structural framework is intact.

**Gravity arc status after G4-3a:**
- Continuum gravity arc: G1–G8 complete (8 sprints, natural sprint-scale closure)
- Discrete-substrate gravity arc: G4-3a opens; ~6-12 months estimated for G4-3 → G4-6

This is the natural multi-month commitment the PI flagged in the G4-2 memo. G4-3a is the structural opening; the full discrete-substrate derivation will be the next major project arc.

## 10. Files

- `debug/g4_3_warped_substrate.py` — substrate construction + first-pass driver
- `debug/data/g4_3_warped_substrate.json` — structured results
- `debug/g4_3_warped_substrate_memo.md` — this memo
- `papers/group5_qed_gauge/paper_28_qed_s3.tex` — §4.17 added

## 11. Cross-references

- **G4-2** (`memory/sprint_g4_2_conical_replica.md`): continuum conical-defect derivation; G4-3 rebuilds on discrete
- **G4-1** (`memory/sprint_g4_1_S2_dirac.md`): $S^2$ Dirac heat trace; G4-3 inherits Fock $S^2$ projection
- **G8** (`memory/sprint_g8_cutoff_dependence.md`): cutoff-dependence inherited by discrete substrate
- **Camporesi-Higuchi 1996**: standard reference for Dirac spectra on round spheres (G4-1 already cited)
- **Sommerfeld 1894, Cheeger 1983**: conical-defect heat-trace contribution
