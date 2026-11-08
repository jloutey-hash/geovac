# Sprint G3 — Scalar Laplacian, TT-tensor Lichnerowicz, and the spinor-bundle specificity of two-term exactness

**Date:** 2026-05-28
**Path:** Gravity Path 1, third sprint. Tests whether Paper 28's two-term exactness on the Dirac sector of $S^3$ propagates to other natural operators — the scalar Laplacian (spin-0) and the Lichnerowicz Laplacian on transverse-traceless symmetric rank-2 tensors (the graviton sector, spin-2).
**Verdict:** POSITIVE-WITH-STRUCTURAL-DISTINCTION. The answer is NO — two-term exactness is **specific to the half-integer-shifted Dirac spectrum**. The scalar and TT-tensor sectors have closed-form heat traces with infinitely many Seeley-DeWitt coefficients, the standard CC continuum expansion. Substantive new content: closed-form identity $a_k = 2\pi^2 / k!$ for the scalar Laplacian on unit $S^3$, plus a more intricate closed form for the TT-tensor sector with mixed half-integer and integer power-of-$t$ contributions.

## 1. The question

Sprint G1 found that the spectral action on $S^3_R$ is EXACTLY two-term modulo exp-small via the structural identity $\zeta_{\rm unit}(-k) = 0$ for all $k \geq 0$ (Paper 28 §4.7, $\zeta$-side of two-term exactness). G2 showed this propagates to $S^3 \times S^1_\beta$ as exact 4D CC Einstein-Hilbert + cosmological constant.

The natural Path-1 follow-on: does this clean structure extend to OTHER operators on $S^3$? Specifically, the operators relevant to gravitons (TT-tensor Lichnerowicz) and to scalar matter (scalar Laplacian).

If YES (extends): GeoVac substrate has uniform two-term exactness across all spin sectors. The framework's "clean gravity" is a property of the manifold $S^3$, not the Dirac bundle.

If NO (specific to Dirac): two-term exactness is a feature of the spinor bundle (half-integer-shifted spectrum). The gravity sector inherits the standard CC continuum expansion.

## 2. Three closed forms

For all three operators on unit $S^3$ we derive closed-form heat traces.

### 2.1 Dirac (Paper 28 two-term exactness)

Spectrum: $|\lambda_n^{D^2}| = (n+3/2)^2$ with multiplicity $g_n = 2(n+1)(n+2)$, $n \geq 0$. Heat trace:
$$\boxed{\;K_{D^2}(t) = \mathrm{Tr}\,e^{-tD^2}\bigm|_{S^3,{\rm unit}} = \frac{\sqrt{\pi}}{2}\,t^{-3/2} - \frac{\sqrt{\pi}}{4}\,t^{-1/2} + O(e^{-\pi^2/t})\;}$$

Exactly two power-law terms. Mechanism: half-integer shift in the spectrum, Jacobi $\theta_2$ inversion. Paper 28 verified at 80 dps.

### 2.2 Scalar Laplacian (new closed form)

Spectrum: $\lambda_n^\Delta = n(n+2)$ with multiplicity $g_n^{\rm scalar} = (n+1)^2$, $n \geq 0$. Heat trace:
$$\boxed{\;K_\Delta(t) = \mathrm{Tr}\,e^{-t\Delta}\bigm|_{S^3,{\rm unit}} = \frac{\sqrt{\pi}}{4} \cdot \frac{e^t}{t^{3/2}} + O(e^{-\pi^2/t})\;}$$

**Substantive new closed form**, verified bit-exactly (rel diff $10^{-50}$) at $t \leq 0.1$. Derivation: let $u = n+1$, then $K_\Delta(t) = e^t \sum_{u \geq 1} u^2 e^{-u^2 t}$. The inner sum is $-\frac{1}{2}\frac{d}{dt}\theta_3(0, e^{-t}) = \frac{\sqrt{\pi}}{4} t^{-3/2}$ + exp-small via Jacobi $\theta_3$ inversion.

### 2.3 TT-tensor Lichnerowicz (new closed form)

Spectrum: $\lambda_n^{TT} = n(n+2) - 2$ for $n \geq 2$ (graviton TT modes, standard de Sitter result). Multiplicity: $g_n^{TT} = 2(n-1)(n+3)$. Heat trace:
$$\boxed{\;K_{\Delta_L^{TT}}(t) = e^{3t}\!\left[\frac{\sqrt{\pi}}{2}\,t^{-3/2} - 4\sqrt{\pi}\,t^{-1/2} + 4\right] + 6 e^{2t} + O(e^{-\pi^2/t})\;}$$

**Also substantive new closed form**, verified bit-exactly (rel diff $5 \times 10^{-51}$) at $t \leq 0.05$. Derivation: let $u = n+1 \geq 3$, then $g_n^{TT} = 2(u-2)(u+2) = 2u^2 - 8$, and $\lambda^{TT} = u^2 - 3$. The sum splits into a smooth Weyl part (Jacobi-asymptotic terms) plus discrete corrections from the missing $u=1, u=2$ levels:
$$K_{TT}(t) = e^{3t}\left[2\sum_{u \geq 1} u^2 e^{-u^2 t} - 8\sum_{u \geq 1} e^{-u^2 t} - 2(1 e^{-t} + 4 e^{-4t}) + 8(e^{-t} + e^{-4t})\right]$$

The Jacobi-asymptotic terms give the half-integer-power content $(\sqrt{\pi}/2) t^{-3/2} - 4\sqrt{\pi}\,t^{-1/2}$. The constant $-8 \cdot (-1/2) = 4$ comes from the $-1/2$ in the $\theta_3$ asymptotic. The $6 e^{-t} \cdot e^{3t} = 6 e^{2t}$ is the discrete correction from the missing $u=1$ level.

## 3. Seeley-DeWitt coefficient comparison

Expanding into standard form $K(t) = (4\pi t)^{-3/2} \sum_k a_k t^k$ at small $t$:

### 3.1 Dirac (Paper 28)

$$a_0^{D^2} = 4\pi^2, \quad a_1^{D^2} = -2\pi^2, \quad a_k^{D^2} = 0\text{ for }k \geq 2.$$

The vanishing of all higher coefficients is the two-term exactness theorem.

### 3.2 Scalar (new closed form)

Expanding $K_\Delta(t) = (\sqrt{\pi}/4) \sum_k t^{k-3/2}/k!$ gives
$$\boxed{\;a_k^\Delta = \frac{2\pi^2}{k!} \quad \forall k \geq 0\;}$$

Beautiful explicit form. Verified numerically (rel diff $10^{-8}$ to $10^{-6}$ for $k = 0, 1$, accuracy degrading for higher $k$ due to fit conditioning).

Values: $a_0 = a_1 = 2\pi^2 \approx 19.74$, $a_2 = \pi^2 \approx 9.87$, $a_3 = \pi^2/3 \approx 3.29$, $a_4 = \pi^2/12 \approx 0.82$, $a_5 = \pi^2/60 \approx 0.16$, $\ldots$. All nonzero, factorially decreasing.

### 3.3 TT-tensor

The TT closed form mixes half-integer powers (from the smooth Weyl part with $e^{3t}$ prefactor) AND integer powers (from the constant $4$ and discrete $6 e^{2t}$ correction). The standard SD expansion $\sum_k a_k t^k$ captures ONLY the half-integer powers; the integer-power content is "non-SD."

Half-integer-power content from $e^{3t}[(\sqrt{\pi}/2) t^{-3/2} - 4\sqrt{\pi}\,t^{-1/2}]$:
- coefficient of $t^{-3/2}$: $\sqrt{\pi}/2 \Rightarrow a_0^{TT} = 4\pi^2$
- coefficient of $t^{-1/2}$: $(\sqrt{\pi}/2)(3) - 4\sqrt{\pi} = -5\sqrt{\pi}/2 \Rightarrow a_1^{TT} = -20\pi^2$
- coefficient of $t^{1/2}$: $(\sqrt{\pi}/2)(9/2) - 4\sqrt{\pi}(3) = -39\sqrt{\pi}/4 \Rightarrow a_2^{TT} = -78\pi^2$

So $a_k^{TT}$ are all nonzero (with explicit factors of $3^k$ from the $e^{3t}$ expansion).

Integer-power content from $4 \cdot e^{3t} + 6 e^{2t}$: contributes Taylor-series corrections $t^0, t^1, t^2, \ldots$ which do not fit the standard SD form. These are "non-CC" structural pieces specific to the truncated $n \geq 2$ spectrum.

## 4. The structural distinction

| Operator | Spectrum | Two-term exact? | Closed-form $a_k$ | Higher curvature |
|---|---|---|---|---|
| Dirac $D^2$ | half-integer $(n+3/2)^2$ | YES (Paper 28) | $a_0 = 4\pi^2$, $a_1 = -2\pi^2$ | $a_k = 0$ for $k \geq 2$ |
| Scalar Laplacian $\Delta$ | integer $n(n+2)$ | NO | $a_k = 2\pi^2/k!$ | all nonzero, factorial decay |
| TT Lichnerowicz $\Delta_L^{TT}$ | integer $n(n+2) - 2$, $n \geq 2$ | NO | mixed half-integer + integer | all nonzero + non-SD pieces |

**Mechanism.** Two-term exactness comes from Jacobi theta inversion. The Dirac uses $\theta_2$ (half-integer-shifted sum), which after differentiation gives exactly $\sqrt{\pi}/2 \cdot t^{-3/2}$ — no other power-law terms. The scalar Laplacian uses $\theta_3$ (integer sum), which gives $\sqrt{\pi}/2 \cdot t^{-1/2}$ for the sum (not derivative), plus the $-1/2$ constant offset, and after multiplication by $e^t$ (from $u^2 - 1$ shift) gets an infinite series of higher powers.

## 5. Implications for GeoVac's gravity sector

**The two-term exactness theorem (Paper 28 / Sprint G1) is a feature of the spinor bundle on $S^3$, not of $S^3$ as a manifold.**

Specifically:
- Dirac's clean two-term form is what makes GeoVac's spectral action on $S^3_R$ EXACTLY Einstein-Hilbert + cosmological constant with no higher-curvature corrections at any order (G1, G2).
- For SCALAR matter coupled to gravity, GeoVac on $S^3$ would inherit the standard CC infinite SD expansion via the scalar Laplacian sector — same as continuum CC.
- For GRAVITONS (spin-2, TT-tensor modes), the Lichnerowicz Laplacian gives the full SD expansion plus additional non-SD discrete corrections from the $n \geq 2$ cutoff of the graviton spectrum. The graviton sector is MORE complex than CC's continuum result, not less.

**Reading.** GeoVac's "clean gravity" claim from G1/G2 is sector-specific. The clean two-term form lives in the spin-1/2 sector. The spin-0 and spin-2 sectors behave like (or are MORE complex than) standard CC continuum.

This sharpens the structural-skeleton-scope statement: the framework's clean structure is specific to the Dirac sector that the Fock projection sources. Gravity, as the spin-2 sector, isn't part of this clean structure — it lives in a separate (and standard CC-like) sector.

## 6. Honest scope

**Reached:**
- Closed-form $K_\Delta(t) = (\sqrt{\pi}/4)\,e^t/t^{3/2}$ for scalar Laplacian on $S^3$, verified bit-exactly.
- Closed-form $a_k^\Delta = 2\pi^2/k!$ for all $k \geq 0$, substantive new structural identity.
- Closed-form for TT-tensor heat trace with mixed half-integer + integer + discrete contributions.
- Structural distinction: Dirac's two-term exactness is half-integer-spectrum-specific.

**Not reached:**
- Construction of GeoVac's native graviton sector. The framework's natural Hilbert space is built on the Dirac operator (spin-1/2). Spin-2 modes don't naturally live in the framework's data; they'd need to be ADDED via CC-style metric variation. The TT-tensor analysis here is an EXTERNAL computation testing what would happen IF gravitons were present.
- Linearized graviton equations from spectral-action variation. To get the Fierz-Pauli graviton kinetic term, one would need to vary the spectral action $\mathrm{Tr}\,f(D^2/\Lambda^2)$ over the metric $g$, get Einstein-Hilbert in the asymptotic, then linearize. This is a multi-month construction beyond a sprint.
- Tests of whether the framework's discrete substrate constrains the graviton sector. The Fock projection is rigid (Paper 23 rigidity theorem); whether this transfers any structural constraints to the spin-2 sector is unknown.

## 7. Three closed forms in summary

$$\boxed{
\begin{aligned}
K_{D^2}(t)\Big|_{S^3} &= \frac{\sqrt{\pi}}{2}\,t^{-3/2} - \frac{\sqrt{\pi}}{4}\,t^{-1/2} + O(e^{-\pi^2/t}) \\[2pt]
K_\Delta(t)\Big|_{S^3} &= \frac{\sqrt{\pi}}{4} \cdot \frac{e^t}{t^{3/2}} + O(e^{-\pi^2/t}) \\[2pt]
K_{\Delta_L^{TT}}(t)\Big|_{S^3} &= e^{3t}\!\left[\frac{\sqrt{\pi}}{2}\,t^{-3/2} - 4\sqrt{\pi}\,t^{-1/2} + 4\right] + 6 e^{2t} + O(e^{-\pi^2/t})
\end{aligned}
}$$

Each tells a structurally distinct story:
- $D^2$: pure two-term, no Taylor corrections (half-integer Dirac spectrum)
- $\Delta$: pure exponential prefactor, all $a_k = 2\pi^2/k!$ nonzero (integer scalar spectrum from $n=0$)
- $\Delta_L^{TT}$: mixed structure with non-SD discrete corrections (integer TT spectrum starting at $n=2$)

## 8. Cross-references

- **Sprint G1** (`debug/g1_spectral_action_S3_radius_memo.md`) — one-parameter spectral action; identifies the $\zeta_{\rm unit}(-k) = 0$ identity equivalent to Paper 28's two-term exactness.
- **Sprint G2** (`debug/g2_spectral_action_S3_x_S1_memo.md`) — extends to $S^3 \times S^1_\beta$ with exact 4D CC form. Same two-term exactness, propagated via heat-kernel factorization.
- **Paper 28 §4.7 and §4.8** (G1, G2 additions) — Dirac two-term exactness and its $\zeta$-zero corollary.
- **Paper 18 §III.7** master Mellin engine — M2 (Seeley-DeWitt) signature. G3 shows that the M2 expansion is BUNDLE-DEPENDENT: half-integer (Dirac) sector is two-term, integer (scalar/tensor) sectors are full series.
- **Paper 23 Fock rigidity theorem** — the $S^3$ Dirac spectrum is forced by the $-Z/r$ Coulomb projection. Sprint G3's spinor-bundle specificity reading is structurally consistent: the framework's clean structure lives in the projection's image, which is the spinor sector.
- **External_input_three_class_partition** memory — the graviton sector inherits Class 1 calibration data (Newton constant, gauge fixing, etc.) from external CC inputs, NOT from GeoVac substrate.

## 9. Forward in gravity arc

G3 closes the "structural-test" phase of Path 1 with three substantive closed forms. Next steps:

- **G4** — cigar geometry parameterized by black-hole mass $M$. Tests gravity-side spectral action on a non-spherical 2-parameter family. Stefan-Boltzmann and $T_H = 1/(8\pi M)$ already in Sprint TD Track 4; G4 would chase Bekenstein-Hawking $S = A/4$ via spectral-action computation. Multi-month; Paper 49 §11 BCFM is the bulk-dual anchor.
- **G5** — spectral action on $S^3 \times \mathbb{R}_\tau$ (decompactified). G2 identified the $\beta \to \infty$ minimum as zero-T de Sitter; G5 would explicit this with non-compact temporal direction.
- **G6 (new candidate after G3)** — graviton-on-GeoVac from CC metric variation. The substantive question is: would the inner-fluctuation prescription of $D$ extended to include metric variation yield a sensible graviton kinetic term on the discrete $S^3$ substrate? Multi-month conceptual sprint.

The gravity arc has natural depth: structural skeleton (G1/G2/G3 sprint-scale) → 4D black-hole spectroscopy (G4) → genuine graviton dynamics (G6). Path 1 done up to the structural-test phase.

## 10. Files produced

- `debug/g3_graviton_lichnerowicz_S3.py` (~330 lines) — driver with 7 steps (scalar closed-form, scalar SD coefficients, numerical SD fit, TT heat trace, TT SD fit, Dirac recap, structural summary)
- `debug/data/g3_graviton_lichnerowicz_S3.json` — structured numerical results
- `debug/g3_graviton_lichnerowicz_S3_memo.md` — this memo (canonical)
