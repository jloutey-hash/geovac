# Sprint G2 — Spectral action on $S^3_R \times S^1_\beta$

**Date:** 2026-05-28
**Path:** Gravity Path 1 follow-on — extension of G1 from one-parameter $\{S^3_R\}$ to the two-parameter family $\{S^3_R \times S^1_\beta\}$ (thermal compactification with antiperiodic boundary conditions).
**Verdict:** POSITIVE. Two-term exactness propagates from $S^3$ to $S^3 \times S^1$ via factorization. The spectral action on the 4-manifold $S^3 \times S^1$ has EXACTLY two terms in its UV asymptotic — $\Lambda^4$ cosmological constant and $\Lambda^2$ Einstein-Hilbert — with NO higher-curvature corrections at any order. Stefan-Boltzmann thermal physics lives in the exp-small (IR) corrections, structurally distinct from the asymptotic. Formal extremum in $(R, \beta)$ at $R\Lambda = 1/\sqrt{6}$ (same as G1), with $\beta \to \infty$ minimum — corresponding to the zero-temperature de Sitter vacuum.

## 1. Setup

Dirac operator on $S^3_R \times S^1_\beta$:
$$D = \gamma^0 \otimes \partial_\tau + I \otimes D_{S^3}$$

with antiperiodic boundary conditions in the temporal direction (fermions). Since $\{\gamma^0, D_{S^3}\} = 0$, the squared Dirac operator is
$$D^2 = -\partial_\tau^2 \otimes I + I \otimes D_{S^3}^2.$$

Squared eigenvalues:
$$\lambda^2_{n,k} = \frac{(n+3/2)^2}{R^2} + \omega_k^2, \qquad \omega_k = \frac{2\pi(k+1/2)}{\beta}, \quad k \in \mathbb{Z}$$

Multiplicity: $g_n = 2(n+1)(n+2)$ (Camporesi-Higuchi on $S^3$, carried through the tensor product).

## 2. Heat kernel factorization

Since $D^2 = D_{S^3}^2 \otimes I + I \otimes (-\partial_\tau^2)$ acts on a tensor-product Hilbert space,
$$\boxed{\;\mathrm{Tr}\,e^{-tD^2}\bigm|_{S^3_R \times S^1_\beta} = K_{S^3}(t, R) \cdot K_{S^1}(t, \beta)\;}$$

with
$$K_{S^3}(t, R) = K_{\rm unit}(t/R^2) = \frac{\sqrt{\pi}}{2} R^3 t^{-3/2} - \frac{\sqrt{\pi}}{4} R t^{-1/2} + O(e^{-\pi^2 R^2 / t})$$

(Paper 28 two-term exactness theorem) and

$$K_{S^1}(t, \beta) = \sum_k e^{-t\omega_k^2} = \frac{\beta}{2\sqrt{\pi t}}\,\theta_4\!\left(0,\, e^{-\beta^2/(4t)}\right) = \frac{\beta}{2\sqrt{\pi t}}\left[1 + 2\sum_{m \geq 1} (-1)^m e^{-m^2\beta^2/(4t)}\right]$$

(Poisson resummation on the antiperiodic Matsubara sum). **Verified numerically** (`debug/g2_spectral_action_S3_x_S1.py`, Steps 1-2): the factorization holds to bit-exact precision; the leading Poisson term $\beta/(2\sqrt{\pi t})$ matches the direct sum to rel diff $< 10^{-50}$ at small $t$.

## 3. Two-term exactness propagates

**Substantive structural finding:** the UV asymptotic of $\mathrm{Tr}\,e^{-tD^2}$ on $S^3_R \times S^1_\beta$ has EXACTLY two power-law terms:

$$\boxed{\;\mathrm{Tr}\,e^{-tD^2}\bigm|_{S^3_R \times S^1_\beta} = \frac{\beta R^3}{4}\,t^{-2} - \frac{\beta R}{8}\,t^{-1} + O(\text{exp small})\;}$$

with NO Taylor terms in $t$. The exp-small corrections are exponential in $\min(R^2/t, \beta^2/t)$.

**Proof.** $K_{S^3}(t, R)$ has only two power-law terms (Paper 28 two-term exactness, Sprint G1's $\zeta_{\rm unit}(-k) = 0$ identity). $K_{S^1}(t, \beta)$ has ONE power-law term ($\beta/(2\sqrt{\pi t})$) plus exp-small. The product has $2 \times 1 = 2$ power-law terms; cross-products of "exp small" times "power law" remain exp small. QED.

**Numerical verification.** From `debug/data/g2_spectral_action_S3_x_S1.json`:

| $(R, \beta, t)$ | $K_{\rm direct}$ | $K_{\rm asymp}$ | rel diff |
|---|---|---|---|
| $(1, 1, 10^{-2})$ | $2.4875 \times 10^3$ | $2.4875 \times 10^3$ | $3 \times 10^{-11}$ |
| $(1, 1, 10^{-3})$ | $2.4988 \times 10^5$ | $2.4988 \times 10^5$ | $1 \times 10^{-17}$ |
| $(0.5, 1, 10^{-3})$ | $3.1188 \times 10^4$ | $3.1188 \times 10^4$ | $6 \times 10^{-51}$ |
| $(2, 2, 10^{-2})$ | $3.995 \times 10^4$ | $3.995 \times 10^4$ | $1 \times 10^{-43}$ |

Bit-exact agreement across all $(R, \beta)$ at $t$ in the UV regime.

## 4. Spectral-zeta identity $\zeta_{S^3 \times S^1}(-k) = 0$

**Corollary.** The spectral zeta of $D^2$ on $S^3_R \times S^1_\beta$ vanishes at all non-positive integers:
$$\zeta_{S^3 \times S^1}(-k) = 0 \qquad \forall k = 0, 1, 2, \dots$$

This propagates G1's identity from $S^3$ to $S^3 \times S^1$. **Proof:** the Mellin transform of $\mathrm{Tr}\,e^{-tD^2}$ relates Taylor coefficients of $K(t)$ at small $t$ to $\zeta(-k)$ values via $\Gamma$-residue at $s = -k$. Since $K(t)$ has no Taylor terms (only power-law $t^{-2}, t^{-1}$ and exp-small), all $\zeta(-k) = 0$.

## 5. Spectral action: CC Einstein-Hilbert + cosmological constant

For Gaussian cutoff $f(x) = e^{-x}$:
$$S_{\rm Gauss}(R, \beta, \Lambda) = \mathrm{Tr}\,e^{-D^2/\Lambda^2} = K(1/\Lambda^2)$$

(factorizes because $f = e^{-x}$ preserves the additive structure of $D^2$). Substituting the two-term asymptotic:
$$\boxed{\;S_{\rm Gauss}(R, \beta, \Lambda) = \frac{\beta R^3}{4}\,\Lambda^4 - \frac{\beta R}{8}\,\Lambda^2 + O(\text{exp small})\;}$$

This is the **Connes-Chamseddine 4D Einstein-Hilbert + cosmological constant form** on the 4-manifold $S^3 \times S^1$:

- $\Lambda^4$ term: cosmological constant. Coefficient $= (\beta R^3 / 4) = \mathrm{Vol}(S^3_R \times S^1_\beta) / (8\pi^2)$ where $\mathrm{Vol}(S^3 \times S^1) = 2\pi^2 R^3 \beta$.
- $\Lambda^2$ term: Einstein-Hilbert. Coefficient $= -(\beta R / 8) = -\frac{1}{96\pi^2} \int R_{\rm scalar} \sqrt{g} d^4x$ where $\int R_{\rm scalar} d^4x = 12\pi^2 R \beta$ on $S^3 \times S^1$ (only $S^3$ contributes; $S^1$ is flat).
- Higher $\Lambda$ powers: **zero** identically (two-term exactness extends).

**Verified numerically** (`g2_spectral_action_S3_x_S1.py` Step 5). Example values at $\Lambda = 10$:

| $(R, \beta)$ | $S_{\rm QM}$ | $S_{\rm asymp}$ | rel diff |
|---|---|---|---|
| $(0.5, 2)$ | $612.5$ | $612.5$ | $7 \times 10^{-44}$ |
| $(1, 1)$ | $2487.5$ | $2487.5$ | $3 \times 10^{-11}$ |
| $(1, 2)$ | $4975.0$ | $4975.0$ | $7 \times 10^{-44}$ |
| $(2, 2)$ | $39950$ | $39950$ | $7 \times 10^{-44}$ |

Bit-exact at large $\Lambda$ across all tested geometries.

## 6. Formal extremum and de Sitter vacuum reading

The asymptotic spectral action $S(R, \beta, \Lambda) = A \beta R^3 - B \beta R$ with $A = \Lambda^4/4 > 0$, $B = \Lambda^2/8 > 0$ has:

**At fixed $\beta, \Lambda$:** $\partial S / \partial R = 0 \Rightarrow R_{\rm crit}^2 = \frac{1}{6\Lambda^2}$, i.e., $R_{\rm crit} \Lambda = 1/\sqrt{6} \approx 0.408$.

**Same as G1.** Adding the $S^1$ factor doesn't change the formal extremum location in $R$ — the $\beta$-dependence factors out and the extremum condition is identical to the one-parameter family.

**At $R = R_{\rm crit}$, $\partial S / \partial \beta < 0$:** the action decreases linearly with $\beta$. So the joint $(R, \beta)$-minimum sits at $R = R_{\rm crit}$ and $\beta \to \infty$.

**Physical reading.** $\beta \to \infty$ is the decompactified $S^1 = \mathbb{R}_\tau$ limit — zero temperature. The spectral action prefers the **zero-temperature de Sitter vacuum** with radius $R_{\rm crit} \sim 1/\Lambda \sim$ Planck length. This is the standard CC prediction: spectral action picks out Planck-scale de Sitter as the preferred vacuum.

The value of $S$ at the extremum: $S(R_{\rm crit}, \beta, \Lambda) = -\beta\Lambda/(12\sqrt{6})$. Negative, linear in $\beta$. Standard interpretation: cosmological-constant problem (the predicted $\Lambda_{cc}$ scale is the cutoff $\Lambda$, not the observed $10^{-60}$ smaller value). Same as in CC's continuum framework.

**Honest scope.** As in G1, the formal extremum sits at $u_{\rm crit} = R_{\rm crit} \Lambda = 1/\sqrt{6}$, which is the boundary of the asymptotic regime of validity (numerical agreement breaks down for $\Lambda R \lesssim 1$). The literal QM spectral action does not have a literal extremum — it's the asymptotic (Reading B in G1) that has one. CC adopts Reading B; under Reading A, the framework does not select a metric.

## 7. Stefan-Boltzmann lives in exp-small corrections

The Stefan-Boltzmann free energy for free fermions is
$$F_{\rm SB}(T) = -\frac{7\pi^2}{180} V T^4 \quad \text{(per Weyl)}, \qquad T = 1/\beta$$

This is an IR / thermodynamic quantity. It does NOT appear in the UV asymptotic of the spectral action. Rather, it lives in the exp-small corrections from the higher Matsubara modes ($m \geq 1$ in the Poisson resummation).

To extract $F_{\rm SB}$, one computes $-T \ln Z = -(1/\beta) \ln \mathrm{Tr}\,e^{-\beta H}$ where $H$ is the Hamiltonian, NOT $S = \mathrm{Tr}\,f(D^2/\Lambda^2)$. These are two different objects: $S$ is the CC effective action (UV regime), while $F$ is the thermal free energy (IR regime).

**Structural finding (G2 corollary):** The CC spectral action expansion on $S^3 \times S^1$ contains the local gravitational invariants (Einstein-Hilbert + cosmological constant) AT THE UV LEVEL, while thermal physics (Stefan-Boltzmann, anomalies) appears at the IR / convergent level. The two regimes are STRUCTURALLY DISTINCT.

This is already known in CC's continuum spectral-action literature, but G2 verifies it explicitly on the GeoVac discrete substrate.

## 8. What this DOES and DOES NOT establish

**Reached:**
- Factorization theorem: $K_{S^3 \times S^1}(t) = K_{S^3}(t) \cdot K_{S^1}(t)$ — exact, structural.
- Two-term exactness propagates: only $t^{-2}$ and $t^{-1}$ terms at small $t$, no Taylor terms.
- $\zeta_{S^3 \times S^1}(-k) = 0$ for all $k \geq 0$ — propagates G1's identity to the 4-manifold.
- Spectral action on the 4-manifold has EXACT 4D CC Einstein-Hilbert + cosmological constant structure with no higher-curvature corrections.
- Formal extremum at $R\Lambda = 1/\sqrt{6}$ (same as G1) with $\beta \to \infty$ — the zero-temperature de Sitter vacuum.

**Not reached:**
- Einstein equations from variation over the FULL metric (not just the scale). Same gap as G1.
- Stefan-Boltzmann or thermal physics from the spectral action UV expansion (lives in IR).
- Cosmological-constant fix to observed scale — same gap as CC's continuum spectral action.
- Gravitational coupling (Newton constant) from GeoVac first principles — set by cutoff $\Lambda$, calibration data.

## 9. Cross-references

- **Sprint G1** (`debug/g1_spectral_action_S3_radius_memo.md`) — one-parameter analog. Same formal extremum at $u_{\rm crit} = 1/\sqrt{6}$.
- **Paper 28 §4.7** (`sec:parametric_spectral_action`, added in v3.4.0) — $\zeta_{\rm unit}(-k) = 0$ theorem on $S^3$. G2 extends to $S^3 \times S^1$.
- **Sprint TD Track 4** (Hawking $T_H = 1/(8\pi M)$ from cigar) — gravity-adjacent Stefan-Boltzmann content. G2 places this in the IR sector, distinct from UV spectral action.
- **Paper 18 §III.7** master Mellin engine — M2 (Seeley-DeWitt) signature. G2 is M2 on $S^3 \times S^1$.
- **CLAUDE.md memory** `external_input_three_class_partition.md` — gravitational coupling and cosmological constant scale are Class 1 calibration data, same as Yukawas and $\alpha$.

## 10. Forward paths (gravity arc, after G2)

- **G3** — linearized graviton spectrum on $S^3$ via Lichnerowicz Laplacian. Tests whether the discrete substrate supports a graviton mode. 2-4 weeks.
- **G4** — cigar geometry parameterized by black-hole mass $M$. Reached via Wick rotation of the Schwarzschild metric. Sprint TD Track 4 already did $T_H = 1/(8\pi M)$; G4 would chase Bekenstein-Hawking $S = A/4$ via spectral-action computation on the cigar. Multi-month; Paper 49 §11 BCFM bulk-dual is the anchor.
- **G5 (new candidate after G2)** — repeat G1 on $S^3 \times \mathbb{R}_\tau$ (decompactified time). The $\beta \to \infty$ limit explored explicitly. Tests whether the formal-extremum reading at zero-temperature de Sitter is the correct interpretation.

## 11. Files produced

- `debug/g2_spectral_action_S3_x_S1.py` (~340 lines) — driver covering 6 steps (factorization, Poisson resummation, two-term asymptotic, SD coefficients, spectral action panel, extremum analysis)
- `debug/data/g2_spectral_action_S3_x_S1.json` — structured numerical results
- `debug/g2_spectral_action_S3_x_S1_memo.md` — this memo
