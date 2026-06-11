# Sprint G4-4 scoping — Warped Dirac operator on the discrete cigar substrate

**Date:** 2026-05-28
**Path:** Gravity arc, opening of the **second-stage multi-month commitment** (G4-3 closed the substrate; G4-4 builds the dynamics that produces $S_{BH}$).
**Verdict:** **POSITIVE-SCOPING-G4-4.** Discrete spinor bundle and warped Dirac operator are architecturally well-defined on the G4-3 substrate $\mathcal{G}_{\rm cigar}$. Three load-bearing falsifiers identified at theorem-grade rigor. Sub-sprint sequence G4-4a/b/c named, with G4-4a (spinor bundle + Riemannian-limit recovery) sized as a 4-8 week sprint. Multi-month commitment $\sim$ 4-7 months end-to-end to first $S_{BH}$ on discrete substrate.

---

## §1. Context

**G4 strategic position.** Sprint TD Track 4 + G4 sequence reproduced Hawking $T_H = 1/(8\pi M)$ and Bekenstein-Hawking $S_{BH} = A/(4G)$ at the **continuum** Connes–Chamseddine spectral-action level (Paper 28 §4.11, §4.15). The structural-skeleton-scope reading per [[geovac_structural_skeleton_scope_pattern]] is that the framework reproduces the **standard** CC derivation of the BH entropy; the next question is whether the framework's distinguishing feature — that it is a **discrete** spectral triple per Paper 38's WH1 PROVEN closure — has any new content for BH entropy.

The G4-3 sub-sprint sequence (2026-05-28 close-of-day) closed five of seven sub-sprints opening the discrete-substrate program: G4-3a, G4-3a-cleanup (Hermitian polar Laplacian), G4-3b (variable warp), G4-3c (naive sweep + proper wedge), G4-3d (continuum Weyl law + UV extension). The G4-3 substrate $\mathcal{G}_{\rm cigar}$ is now a **scalar Laplacian** infrastructure. **G4-4 builds the Dirac operator on this substrate.**

The Dirac operator is the load-bearing object for $S_{BH}$ via two independent paths:
1. **Spectral action**: $S_{\rm spec} = \text{Tr}\,f(D/\Lambda)$; per G7 (Paper 28 §4.13), matching to Einstein–Hilbert gives $G_N = 6\pi/\Lambda^2$.
2. **Replica / conical defect**: $S_{BH} = -dI_E/d\alpha|_{\alpha=1}$ where $I_E$ is the Euclidean Dirac heat-trace action on a cone (G4-2, Paper 28 §4.15).

Both paths require the warped Dirac, NOT the scalar Laplacian.

## §2. Target physics

### §2.1 Continuum warped Dirac on the cigar

The Euclidean Schwarzschild cigar's geometry, near horizon, is
$$ds^2 = d\rho^2 + \rho^2\,d\phi^2 + r(\rho)^2\,d\Omega_2^2$$
with $r(\rho) = r_h\sqrt{1 + (\rho/r_h)^2}$ (G4-3b smooth-tip warp) and $\phi$ periodic with apex angle $2\pi\alpha$ ($\alpha = 1$: smooth disk; $\alpha \ne 1$: conical singularity at $\rho = 0$).

The 4D Dirac operator on this warped product splits as
$$D_{\rm cigar} = D_{D^2_\alpha} \otimes I_{S^2} + \gamma^5_{D^2} \otimes \frac{D_{S^2}}{r(\rho)}$$
(in any chirality-grading convention compatible with the four spacetime dimensions; $\gamma^5_{D^2}$ is the 2D chirality matrix on the disk sector).

The squared Dirac:
$$D_{\rm cigar}^2 = D_{D^2_\alpha}^2 \otimes I_{S^2} + I_{D^2} \otimes \frac{D_{S^2}^2}{r(\rho)^2} + \text{warp-coupling terms involving } r'(\rho)$$

### §2.2 Heat trace and $S_{BH}$

The heat trace
$$K_{\rm cigar}^{\rm Dirac}(t; \alpha) := \text{Tr}_{\rm spinor}\,e^{-t D_{\rm cigar}^2}$$
gives the Euclidean action $I_E(\alpha) = -\int_0^\infty (dt/t) K^{\rm Dirac}(t; \alpha) f(t\Lambda^2)$ (cutoff-regularized) and $S_{BH} = -dI_E/d\alpha|_{\alpha=1}$ per replica method.

**G4-4 target:** compute $K_{\rm cigar}^{\rm Dirac}(t; \alpha)$ on the discrete substrate, extract $S_{BH}^{\rm discrete}$, compare to continuum $S_{BH} = A/(4G)$.

## §3. Discrete substrate inputs (from G4-3)

**Already built (G4-3a/cleanup/b/c/d/UV-T2):**

- **Radial sector** $\mathbb{Z}_+(a)|_{N_\rho}$ with sites $\rho_k = k a$ for $k = 1, \ldots, N_\rho$ (the apex $\rho = 0$ is excluded as a singular point).
- **Hermitian polar Laplacian** (G4-3a-cleanup): symmetric tridiagonal $H_{\rm rad}$ via the substitution $u = \sqrt{\rho}\,f$, Bessel-zero asymptotic convergence verified.
- **Azimuthal sector** $\mathbb{Z}/N_\phi$ for the disk; periodic BC; conical-defect $\alpha = N_\phi/N_0$ (T1 wedge lattice).
- **Variable warp** $r(\rho)$ (G4-3b): smooth-tip with asymptotic Schwarzschild $r(\rho)/\rho \to 1$. Heat-trace ratios $K_{\rm var}/K_{\rm const} \in [0.89, 2.19]$.
- **$S^2$ sector** Fock projection at $l_{\max}$: integer spectrum $l(l+1)$ with multiplicity $2l+1$ (Paper 22 angular sparsity substrate).
- **UV regime verified** (T2, G4-3d-UV): Weyl law recovers within 1% at $t = 0.1$, $N_\phi = 192$.

**Need to add for G4-4:**

- **Discrete spinor bundle** over $\mathcal{G}_{\rm cigar}$
- **Discrete warped Dirac operator** $D_{\rm cigar}^{\rm disc}$
- **Spinor heat trace** $\text{Tr}_{\rm spinor}\,e^{-t (D_{\rm cigar}^{\rm disc})^2}$

## §4. Discrete spinor bundle on $\mathcal{G}_{\rm cigar}$

### §4.1 Bundle structure

The continuum 4D spinor bundle on the cigar factorizes (at constant warp) as
$$S(D^2 \times S^2) = S(D^2) \otimes S(S^2)$$
with $S(D^2)$ a rank-2 complex bundle (2D Dirac), $S(S^2)$ a rank-2 complex bundle (Weyl). Total rank: 4 = $4 \times \binom{4}{2}/2 = 4$ (4D Cl(4,0) γ-algebra has $2^{[4/2]} = 4$ minimal irrep).

### §4.2 Discrete spinor space

$$\mathcal{H}^{\rm spin}_{\rm cigar} = \mathcal{H}^{\rm spin}_{D^2} \otimes \mathcal{H}^{\rm spin}_{S^2}$$

where:

**$\mathcal{H}^{\rm spin}_{D^2}$** is the spinor space over the discrete disk substrate. At each radial site $\rho_k$, a 2-component complex vector lives (the 2D Dirac spinor); coupled to the azimuthal $N_\phi$ sites via the angular spin connection (1/2-twist under $\phi \to \phi + 2\pi$ for fermions, opposite of bosons).
$$\dim \mathcal{H}^{\rm spin}_{D^2} = 2 \cdot N_\rho \cdot N_\phi$$

**$\mathcal{H}^{\rm spin}_{S^2}$** is the Camporesi–Higuchi spinor space on $S^2$ (Paper 28 §4.14 / Paper 23 §V): integer Dirac spectrum $|\lambda_n| = n + 1$ with multiplicity $4(n+1)$ at half-integer $j = n + 1/2$. At $l_{\max}$ cutoff:
$$\dim \mathcal{H}^{\rm spin}_{S^2} = \sum_{n=0}^{l_{\max}} 4(n+1) = 2(l_{\max}+1)(l_{\max}+2)$$

**Total spinor Hilbert space:**
$$\dim \mathcal{H}^{\rm spin}_{\rm cigar} = 2 \cdot N_\rho \cdot N_\phi \cdot 2(l_{\max}+1)(l_{\max}+2)$$
At sprint-scale parameters $(N_\rho, N_\phi, l_{\max}) = (200, 120, 4)$: $\dim \sim 2.9 \cdot 10^6$. Eigenvalues solvable via Lanczos / shift-invert at this scale.

### §4.3 Fermionic anti-periodic BC

At Hawking temperature $T_H = 1/(8\pi M)$, the Euclidean $\phi$-circle has fermionic anti-periodic BC: $\psi(\phi + 2\pi) = -\psi(\phi)$. Discretely: the azimuthal Laplacian eigenvalues shift from $(2/h_\phi)^2 \sin^2(\pi k / N_\phi)$ for bosons to $(2/h_\phi)^2 \sin^2(\pi (k + 1/2)/N_\phi)$ for fermions. The fermionic spectrum has NO zero mode, consistent with the absence of fermion thermal IR divergence.

## §5. Discrete warped Dirac operator

### §5.1 Constant-warp (G4-4a target)

At $r(\rho) = r_h$ constant, the discrete warped Dirac factorizes:
$$D_{\rm cigar}^{\rm disc, const} = D_{D^2}^{\rm disc} \otimes I_{S^2}^{\rm spin} + \gamma^5_{D^2} \otimes \frac{D_{S^2}^{\rm spin}}{r_h}$$

**$D_{D^2}^{\rm disc}$ (2D disk Dirac):**
$$D_{D^2}^{\rm disc} = \gamma^1 \partial_\rho^{\rm disc} + \gamma^2 \frac{1}{\rho_k}\,\partial_\phi^{\rm disc}$$
with Cl(2,0) γ-matrices $\gamma^1 = \sigma_1$, $\gamma^2 = \sigma_2$ (Pauli). The finite-difference $\partial_\rho^{\rm disc}$ uses centered FD with Dirichlet zero at $k = 0, N_\rho + 1$ (consistent with G4-3a-cleanup substitution). The $\partial_\phi^{\rm disc}$ uses centered FD on $\mathbb{Z}/N_\phi$ with anti-periodic BC.

**$D_{S^2}^{\rm spin}$ (Camporesi-Higuchi spinor Dirac on $S^2$):**
Diagonal in $(n, m_j)$ basis with eigenvalues $\pm(n+1)$ and multiplicity $2(n+1)$ at each sign.

**Squared Dirac:**
$$(D_{\rm cigar}^{\rm disc, const})^2 = (D_{D^2}^{\rm disc})^2 \otimes I_{S^2} + I_{D^2} \otimes \frac{(D_{S^2}^{\rm spin})^2}{r_h^2}$$
(cross term cancels because $\{D_{D^2}, \gamma^5_{D^2}\} = 0$).

### §5.2 Variable-warp (G4-4b target)

At variable $r(\rho)$:
$$D_{\rm cigar}^{\rm disc, var} = D_{D^2}^{\rm disc} \otimes I_{S^2} + \gamma^5_{D^2}(\rho) \otimes \frac{D_{S^2}^{\rm spin}}{r(\rho)} + \text{warp-derivative terms}$$

The warp-derivative term involves $r'(\rho)$ and arises from the spin connection on the warped product. Standard form (Camporesi 1996):
$$D_{\rm cigar}^{\rm var} = D_{D^2} \otimes I + \gamma^5_{D^2} \otimes \frac{D_{S^2}}{r(\rho)} + \frac{r'(\rho)}{r(\rho)} \gamma^\rho \otimes I_{S^2}$$

For smooth-tip $r(\rho) = r_h\sqrt{1 + (\rho/r_h)^2}$: $r'(\rho)/r(\rho) = \rho/(\rho^2 + r_h^2)$. Tip-regular ($\to \rho/r_h^2$ as $\rho \to 0$), asymptotic-free ($\to 1/\rho$ as $\rho \to \infty$).

### §5.3 Heat trace

$$K_{\rm cigar}^{\rm Dirac}(t; \alpha) = \text{Tr}_{\rm spinor}\,e^{-t (D_{\rm cigar}^{\rm disc})^2}$$

Constant-warp factorization at the heat-trace level:
$$K_{\rm cigar}^{\rm Dirac, const}(t; \alpha) = K_{D^2}^{\rm Dirac}(t; \alpha) \cdot K_{S^2_{r_h}}^{\rm Dirac}(t)$$

The disk-Dirac heat trace at $\alpha = 1$ is the new computational object G4-4 must produce. The $S^2$-Dirac heat trace is already characterized (G4-1, Paper 28 §4.14): full Seeley-DeWitt series, $a_1^{S^2,{\rm Dirac}} = -4\pi/3$.

## §6. Load-bearing falsifiers

Following the L3a-1 / L3b-2 first-move discipline, G4-4 lives at finite cutoff and depends on three load-bearing falsifiers passing at every tested $(N_\rho, N_\phi, l_{\max})$ panel cell:

### (F1) Riemannian-limit factorization at constant warp

At $r(\rho) = r_h$, $\alpha = 1$:
$$\text{Tr}\,e^{-t (D_{\rm cigar}^{\rm disc, const})^2} \stackrel{?}{=} \text{Tr}\,e^{-t (D_{D^2}^{\rm disc})^2} \cdot \text{Tr}\,e^{-t (D_{S^2}^{\rm spin})^2/r_h^2}$$
**Bit-exact** at every test point (residual = 0 in float64). The operator-level factorization is preserved at the heat-trace level by trace cyclicity + tensor decomposition.

### (F2) Chirality grading

$\{\gamma^5_{\rm cigar}, D_{\rm cigar}^{\rm disc}\} = 0$ at every cell. $\gamma^5_{\rm cigar} = \gamma^5_{D^2} \otimes \gamma^5_{S^2}$ in the chiral basis convention; the 4D chirality is the product of disk and sphere chiralities.

### (F3) Continuum heat-trace recovery

$K_{D^2}^{\rm Dirac, disc}(t = 0.5; N_\rho = 200, N_\phi = 192, \alpha = 1)$ matches continuum Weyl-Selberg prediction within 10% (matching T2 G4-3d-UV scalar verification). This is the falsifier that pins the spinor-substrate operationally.

**Quick scope check.** Continuum disk Dirac small-$t$ asymptotics:
$$K_{D^2}^{\rm Dirac}(t) \sim \frac{2 A_{D^2}}{4\pi t} - \frac{L_{D^2}}{4\sqrt{\pi t}} + O(t^0)$$
The leading $1/t$ coefficient is **twice** the scalar (factor of 2 = rank of 2D spinor bundle). Discrete recovery to 10% is plausibly reachable at the same UV refinement that worked for the scalar (T2).

### Open falsifier (deferred to G4-4b)

(F4) Variable-warp asymptotic-free behavior: $K_{\rm cigar}^{\rm Dirac, var}(t \to 0)$ matches the warped Weyl law with appropriate boundary corrections. Tests whether the discrete substrate handles the warp-derivative spin connection consistently.

## §7. G4-4 sub-sprint sequence

Six sub-sprints, sized per the L3a-1 / L3b-2 cadence (1 sub-sprint $\approx$ 1 week of focused PM time + sub-agent dispatch for parallel verification):

| Sub-sprint | Scope | Effort |
|---|---|---|
| **G4-4a** | Constant-warp Dirac on disk; load-bearing falsifiers F1, F2, F3 | 4-8 weeks |
| **G4-4b** | Variable-warp Dirac with smooth-tip; falsifier F4; tip-regularity | 4-8 weeks |
| **G4-4c** | Wedge-Dirac (conical-defect $\alpha \ne 1$); reciprocal cancellation test | 4 weeks |
| **G4-4d** | Spinor heat trace at small $t$; Seeley-DeWitt $a_0, a_1, a_2$ extraction | 2-4 weeks |
| **G4-4e** | Anti-periodic vs periodic BC: fermion vs boson sectors | 2 weeks |
| **G4-4f** | Replica-method preparation: heat-trace differential $dK/d\alpha\big|_{\alpha=1}$ | 4-6 weeks |

**Total G4-4 commitment: 5-8 months.**

After G4-4: G4-5 (discrete replica method, 2-4 months) and G4-6 (full $S_{BH}$ derivation from G4-4 + G4-5 fusion, 2-4 months). End-to-end **G4-4 → G4-6: 9-16 months**.

## §8. First-move plan for G4-4a

### §8.1 Module structure

New production module `geovac/gravity/warped_dirac.py` (estimated ~600 lines):

```
- DiscreteSpinorBundle class
  - radial spinor sites
  - azimuthal spinor sites (anti-periodic BC)
  - S^2 spinor Fock projection (Camporesi-Higuchi import)

- DiscreteDiskDirac function
  - gamma matrices (sigma_1, sigma_2)
  - finite-difference rho-derivative
  - finite-difference phi-derivative
  - Hermiticity verification

- WarpedDirac function (constant warp)
  - tensor product assembly D_cigar = D_disk x I + gamma5 x D_S2/r_h
  - chirality grading verification
  - heat trace computation (full diagonalization or Lanczos)
```

### §8.2 Test architecture

New test file `tests/test_warped_dirac.py` (estimated ~25 tests + 3 slow):
- Riemannian-limit factorization at constant warp (F1, bit-exact)
- Chirality anticommutation (F2)
- Discrete vs continuum disk-Dirac at small $t$ (F3, within 10%)
- Anti-periodic BC: no zero mode for fermions
- gamma matrix conventions (Hermitian, square to I, anticommute)
- Bessel-zero asymptotic convergence at large $\rho$
- $S^2$ spinor Fock sector matches Paper 23 §V analytical

### §8.3 Driver and memo

- `debug/g4_4a_constant_warp_dirac.py`: driver
- `debug/g4_4a_constant_warp_dirac_memo.md`: closure memo (target verdict POSITIVE-G4-4a-VERIFIED with F1/F2/F3 bit-exact)
- `debug/data/g4_4a_constant_warp_dirac.json`: numerical panel

### §8.4 Sequencing for G4-4a (4-8 week breakdown)

- **Week 1-2**: spinor bundle construction; γ-matrix conventions; Riemannian-limit factorization (F1)
- **Week 3-4**: disk-Dirac at $\alpha = 1$; chirality grading (F2); small-$t$ heat trace
- **Week 5-6**: UV-regime $N_\phi$ refinement on disk-Dirac; continuum recovery (F3)
- **Week 7-8**: anti-periodic BC sector; documentation; closure memo; sub-sprint review

## §9. Multi-month commitment scoping

**The G4-4 commitment is structurally additive on top of G4-3.** The G4-3 substrate is intact; G4-4 adds the spinor bundle and Dirac dynamics. No re-derivation of G4-3 results.

**Risk tier 1 (low):**
- Constant-warp factorization (F1): operator-level identity, almost-guaranteed bit-exact
- Chirality grading (F2): convention choice, verifiable in closed form
- Anti-periodic BC: standard FD modification

**Risk tier 2 (medium):**
- Continuum recovery at small $t$ (F3): inherits T2's UV regime; should work with $N_\phi \geq 144$
- Tip-regularity of warp-derivative term in G4-4b: requires careful $\rho \to 0$ behavior

**Risk tier 3 (higher):**
- Spinor heat trace at conical defect $\alpha \ne 1$ (G4-4c): G4-3c proper wedge negative-with-diagnosis on the SCALAR side suggests this won't be bit-exact for the spinor either; the SC tip term lives below discretization floor at sprint scale
- Replica-method $dK/d\alpha$ extraction (G4-4f): the differential at $\alpha = 1$ is well-defined in continuum but discrete substrate has $\alpha$ in integer intervals; need analytic continuation in $\alpha$ (the actual content of "replica method")

**The dominant risk is G4-4c + G4-4f**, which together provide the load-bearing $S_{BH}$ extraction. Mitigation:
- Pursue both (a) G4-4c on discrete substrate AND (b) "continuum Sommerfeld-Cheeger on the disk-Dirac side analytically" to give a continuum-vs-discrete bridge that may close the literal extraction gap.
- Defer literal $1/12 \cdot (1/\alpha - \alpha)$ extraction to G4-5 (discrete replica method) if G4-4c + G4-4f surface obstructions.

## §10. Honest scope (sprint vs multi-month)

This memo is a **scoping document**. It does NOT contain new computational results; it:
- Confirms G4-3 substrate is sufficient input
- Names the architectural extensions G4-4 needs (spinor bundle + Dirac operator)
- Names load-bearing falsifiers F1, F2, F3 (operationally testable at sprint scale)
- Names sub-sprint sequence G4-4a/b/c/d/e/f with effort estimates
- Names risk-tier classification and mitigation

The PI may launch G4-4a as the next focused sprint or defer to a future cycle. The sub-sprint sequence is designed for parallelizable PM sessions: G4-4a/b can run sequentially (architectural), G4-4c is downstream, G4-4d/e are spectral-extraction work that can parallel to G4-4b.

## §11. Cross-references

- G4-3 substrate scoping: `debug/g4_3_warped_substrate_memo.md`
- G4-3a-cleanup Hermitian polar Laplacian: `debug/g4_3a_cleanup_hermitian_polar_memo.md`
- G4-3b variable warp: `debug/g4_3b_variable_warp_memo.md`
- G4-3c proper wedge (T1, same sprint): `debug/g4_3c_proper_wedge_memo.md`
- G4-3d UV extension (T2, same sprint): `debug/g4_3d_uv_extension_memo.md`
- G4-1 $S^2$ Dirac: Paper 28 §4.14
- G4-2 conical replica derivation: Paper 28 §4.15
- L3a-1 first-move template: `debug/l3a_1_lorentzian_operator_system_memo.md`
- Camporesi 1996, "Harmonic analysis and propagators on homogeneous spaces": warped-product spin connection
- Paper 23 §V: Camporesi-Higuchi $S^2$ spinor spectrum
- Paper 22: angular sparsity substrate

## §12. Files

- `debug/g4_4_warped_dirac_scoping_memo.md` — this memo
- (Future, when G4-4a launches): `geovac/gravity/warped_dirac.py`, `tests/test_warped_dirac.py`, `debug/g4_4a_constant_warp_dirac.py`, `debug/data/g4_4a_constant_warp_dirac.json`, `debug/g4_4a_constant_warp_dirac_memo.md`
