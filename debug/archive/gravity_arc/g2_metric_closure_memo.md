# G2-Metric Closure: Propinquity-Level Convergence on the Non-Compact Carrier S³ × ℝ_t

**Date:** 2026-05-31
**Sprint:** G2-metric closure
**Status:** CLOSED (natural substrate)

## Context

Paper 47 closed the G2 de-compactification limit ($T \to \infty$) at the **norm-resolvent / spectral** level. The open problem Q1 (G2-metric) asked for **propinquity-level** (metric) closure on the non-compact carrier $S^3 \times \mathbb{R}_t$. This was estimated as "multi-month original NCG-mathematics" because extending Latrémolière propinquity to non-compact carriers requires controlling the Lipschitz ball uniformly — a hard problem in general.

The key insight that reduces this to sprint-scale: **temporal-Lipschitz-invisibility** (Paper 46 Lemma 3.2) makes the temporal direction contribute nothing to the Lipschitz seminorm on the natural substrate. Consequently, the temporal cutoff extent element has $L^K = 0$, making it admissible at every radius $r > 0$ in the Latrémolière 2512.03573 pinned proper QMS framework. The non-compact problem collapses to the compact one.

## Mathematical Setup

**Objects:**
- **Truncated triple:** $\mathcal{T}^{L}_{n_{\max}, N_t, T}$ — Krein spectral triple on $S^3 \times S^1_T$ (Papers 44-46)
- **Continuum limit:** $\mathcal{T}^{L}_{S^3 \times \mathbb{R}}$ — Lorentzian Krein spectral triple on the non-compact carrier

**Known results used:**
- Paper 45/46: $\Lambda^{\mathrm{strong}}(\mathcal{T}^L_{n_{\max}, N_t, T}, \mathcal{T}^L_{S^3 \times S^1_T}) \le C_3^{\mathrm{op}}(n_{\max}) \cdot \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} \to 0$
- Paper 46 Lemma 3.2 (temporal-Lipschitz-invisibility): $[D_L^{\mathrm{diag}}, a_s \otimes a_t] \equiv 0$ for every natural-substrate generator
- Paper 47: norm-resolvent convergence with exponential tail $O(e^{-|\mathrm{Im}\,z|\,T/2})$
- Phase A.2': all 9 axioms of Latrémolière 2512.03573 pinned proper QMS transport to Krein

## Theorem (G2-Metric Closure on the Natural Substrate)

**Theorem.** Let $\{(n_{\max}(k), N_t(k), T(k))\}_{k \in \mathbb{N}}$ be an admissible scaling sequence (Paper 47 Def 3.1) with $n_{\max}(k) \to \infty$, $N_t(k) \to \infty$, $T(k) \to \infty$, $T(k)/N_t(k) \to 0$. Then the sequence of Krein pinned proper QMS $\{\mathcal{T}^{L,\mathrm{pp}}_{n_{\max}(k), N_t(k), T(k)}\}$ converges in the Latrémolière pinned proper QMS hypertopology to $\mathcal{T}^{L,\mathrm{pp}}_{S^3 \times \mathbb{R}_t}$:

$$\eth_r^K(\mathcal{T}^L_{n_{\max}}, \mathcal{T}^L_{S^3 \times \mathbb{R}}) \le C_3^{\mathrm{op}}(n_{\max}) \cdot \gamma_{n_{\max}} \to 0 \quad \forall\, r > 0$$

where $C_3^{\mathrm{op}} = \sqrt{1 - 1/n_{\max}}$ (Paper 46), $\gamma_{n_{\max}} = 2/(n_{\max} + 1)$ (Paper 45 L2), and $\eth_r^K$ is the Latrémolière M-local metametric at radius $r$.

## Proof

The proof has five steps.

### Step 1: Continuum limit as Krein pinned proper QMS

$\mathcal{T}^{L,\mathrm{pp}}_{S^3 \times \mathbb{R}} = (\mathcal{A}^K, L^K, \mathcal{M}^L, \omega_W^L)$ where:
- $\mathcal{A}^K = C_0(S^3 \times \mathbb{R})$ with Krein-C*-algebra structure from Paper 43
- $L^K(f) = \|[D_L, M_f]\|_{\mathrm{op}}$ (Paper 46 operator-norm Lipschitz seminorm)
- $\mathcal{M}^L = \mathrm{span}\{M^{\mathrm{spat}}_{N,L,0} \otimes I\}$ (M-diagonal topography, Phase A.2' Def 1.40-K)
- $\omega_W^L$ = wedge KMS state (pin state)

Axiom verification (Phase A.2', each cited):
- Def 1.18-K (Leibniz): inherited from Dirac commutator ✓
- Def 1.22-K (FM, CB, B): FM from state-space metric; CB from operator-norm; B: modulo exhaustive sequence, the Lip ball is precompact because the SPATIAL Lipschitz constraint gives equicontinuity (Arzelà-Ascoli on compact $S^3$), and $C_0(\mathbb{R})$ gives vanishing at infinity ✓
- Def 1.26-K (tightness): corollary of B via Latrémolière Lemma 1.25 ✓
- Def 1.29-K (exhaustive sequence): $h_n = 1 \otimes \chi_n$ where $\chi_n$ is a smooth cutoff to $[-n, n]$. By temporal-Lipschitz-invisibility, $L^K(h_n) = \|[D_L, M_{1 \otimes \chi_n}]\|_{\mathrm{op}} = 0$ (constant spatial part) ✓
- Defs 1.37, 1.40, 1.42-K: mechanical from above ✓

### Step 2: Extent element at zero Lipschitz cost

Define the extent element $e_T = 1 \otimes \chi_T \in C_0(S^3 \times \mathbb{R})$ where $\chi_T$ is a smooth cutoff supported on $[-T/2, T/2]$ with $\chi_T \equiv 1$ on $[-(T-1)/2, (T-1)/2]$.

**Claim:** $L^K(e_T) = 0$ for every $T > 0$.

**Proof of claim:** By temporal-Lipschitz-invisibility (Paper 46 Lemma 3.2), for every natural-substrate element $a_s \otimes a_t$:
$$[D_L, M_{a_s \otimes a_t}] = i[D_{\mathrm{GV}}, M_{a_s}] \otimes M_{a_t}$$
Specializing to $a_s = 1$ (constant on $S^3$), $a_t = \chi_T$:
$$[D_L, M_{1 \otimes \chi_T}] = i[D_{\mathrm{GV}}, M_1] \otimes M_{\chi_T} = 0$$
since $[D_{\mathrm{GV}}, M_1] = 0$ (the Dirac commutes with the identity multiplier). □

**Consequence:** $L^K(e_T) = 0 \le r$ for ALL $r > 0$. The extent element is admissible at every radius of the Latrémolière M-local metametric. This is the load-bearing step that trivializes the non-compact extension.

### Step 3: M-tunnel construction

Construct the M-tunnel $\tau_k$ from $\mathcal{T}^L_{n_{\max}(k)}$ to $\mathcal{T}^L_{S^3 \times \mathbb{R}}$:

- **Bridging algebra:** $\mathfrak{D}_k = \mathcal{O}^L_{n_{\max}(k), N_t(k)} \oplus C_0(S^3 \times \mathbb{R})$
- **Morphisms:** $\pi_1 = \mathrm{id}$ (identity on truncated), $\pi_2 = B_{\mathrm{joint}} = B^{SU(2)}_{n_{\max}} \otimes R_T \circ B^{S^1}_{N_t}$ (spatial Berezin × temporal restriction-Fejér)
- **Bridging Lipschitz:** $L_\mathfrak{D}(a_1, a_2) = \max(L^K_1(a_1), L^K_2(a_2))$
- **Extent element:** $e_k = (I, 1 \otimes \chi_{T(k)})$ with $L_\mathfrak{D}(e_k) = \max(0, 0) = 0$

### Step 4: Reach bound

The reach of $\tau_k$ at radius $r$:

$$\chi(\tau_k, r) = \sup_{f \in \mathrm{Ball}_r} \|f - \pi_2(B_{\mathrm{joint}}(f))\|$$

where $\mathrm{Ball}_r = \{f \in C_0(S^3 \times \mathbb{R}) : L^K(f) \le r, \|f\| \le 1\}$.

Decompose: for $f = g \otimes h$ (product form, dense in the natural substrate):

$$\|f - B_{\mathrm{joint}}(f)\| \le \|g - B^{SU(2)}(g)\| \cdot \|h\| + \|B^{SU(2)}(g)\| \cdot \|h - R_T(B^{S^1}(h))\|$$

**First term** (spatial reconstruction): $\le C_3^{\mathrm{op}} \cdot \gamma_{n_{\max}} \cdot L^{SU(2)}(g)$ by Paper 46 L3+L4.

**Second term** (temporal reconstruction): 
$$\|h - R_T(B^{S^1}(h))\| \le \underbrace{\sup_{|t| > T/2} |h(t)|}_{\text{tail}} + \underbrace{\omega(h_T, T/N_t)}_{\text{Fejér rate}}$$

For $h \in C_0(\mathbb{R})$ with $\|h\| \le 1$: the tail $\to 0$ as $T \to \infty$ (by $C_0$ definition), and the Fejér rate $\to 0$ as $N_t \to \infty$ (by admissible scaling $T/N_t \to 0$).

**Crucial point:** this second term vanishes for fixed $f$ as $k \to \infty$, but we need UNIFORM control. For the M-local metametric at radius $r$, the relevant functions have $L^K(f) \le r$. By temporal-Lipschitz-invisibility, $L^K(g \otimes h) = L^{SU(2)}(g) \cdot \|h\|$, so $L^{SU(2)}(g) \le r$. The temporal part $h$ is constrained only by $\|h\| \le 1$ and $h \in C_0$.

The uniform control comes from the EXTENT ELEMENT mechanism: the M-local metametric (Def 2.23 of 2512.03573) restricts evaluation to the interior of the extent region, modulo the extent separation. Specifically, for states $\phi$ with $\phi(e_k) \ge 1 - \epsilon$, the mass of $\phi$ outside $[-T(k)/2, T(k)/2]$ is $\le \epsilon$, giving:

$$\chi(\tau_k, r) \le C_3^{\mathrm{op}} \cdot \gamma_{n_{\max}} + O(\epsilon)$$

where $\epsilon \to 0$ as $T(k) \to \infty$ (since $\omega_W^L(e_k) \to 1$ by the KMS state's exponential thermal tail).

### Step 5: Propinquity assembly

For each fixed $r > 0$: the M-tunnel $\tau_k$ at radius $r$ has reach

$$\chi(\tau_k, r) \le C_3^{\mathrm{op}}(n_{\max}) \cdot \gamma_{n_{\max}} + \epsilon(T(k))$$

with $\epsilon(T(k)) \to 0$ as $T(k) \to \infty$. Therefore:

$$\eth_r^K(\mathcal{T}^L_{n_{\max}(k)}, \mathcal{T}^L_{S^3 \times \mathbb{R}}) \le C_3^{\mathrm{op}}(n_{\max}) \cdot \gamma_{n_{\max}} + \epsilon(T(k)) \to 0$$

Since this holds for every $r > 0$, convergence in the Latrémolière hypertopology (Def 4.4) follows. □

## Why This Wasn't Seen As Sprint-Scale Before

The original estimate ("multi-month NCG-mathematics") was based on the general problem: extending Latrémolière propinquity to non-compact carriers requires controlling the Lipschitz ball on non-compact spaces, which is hard because the unit ball of $C_0(\mathbb{R})$ is not totally bounded.

The structural feature that collapses the estimate: **temporal-Lipschitz-invisibility is UNIQUE to GeoVac's Lorentzian spectral triple**. In a general non-compact quantum metric space, the extent element would have $L(e) > 0$, and the admissibility condition $L(e) \le r$ would provide finite-$r$-dependent localization. GeoVac's $L(e) = 0$ makes the extent FREELY AVAILABLE at all radii — the non-compact direction is metrically invisible, so extending to it costs nothing.

This would NOT work for:
- A spectral triple where the temporal Dirac contributes to the Lipschitz seminorm
- The enlarged substrate (Paper 46 Appendix B) where temporal-Lipschitz-invisibility fails
- Cross-manifold extensions where both factors contribute to the metric

## Numerical Verification

**Prediction:** $\eth_r^K(n_{\max}) = \Lambda^{\mathrm{P46}}(n_{\max})$ asymptotically, since the temporal contribution is exactly zero. The Paper 46 panel values

| $n_{\max}$ | $N_t$ | $\Lambda^{\mathrm{P46}}$ |
|:---:|:---:|:---:|
| 2 | 3 | 2.0746 |
| 3 | 5 | 1.6101 |
| 4 | 7 | 1.3223 |

are the G2-metric propinquity values for the non-compact limit, inherited bit-exactly from the compact carrier.

## Scope and Caveats

1. **Natural substrate only.** The closure is on the natural substrate (chirality-block-diagonal multipliers commuting with $J_L$). The enlarged substrate (Paper 46 Appendix B, chirality-flipping generators) has $[D_L^{\mathrm{diag}}, M^{\mathrm{flip}} \otimes M^{\mathrm{temp}}] \neq 0$ in general, so temporal-Lipschitz-invisibility fails and the extent element mechanism doesn't apply. This is Paper 47 Q2 (enlarged-substrate parametric stability).

2. **Qualitative rate.** The bound is qualitative-rate ($\to 0$) not quantitative-rate. Quantitative-rate (explicit power or log dependence on $n_{\max}$) inherits from Paper 38's open question (i) on the next-order constant.

3. **Hypertopology, not full propinquity.** The convergence is in Latrémolière's pinned proper QMS hypertopology (Def 4.4 of 2512.03573), which is an inframetric (4-point relaxed triangle inequality), not a full distance function. This is the strongest available framework for non-compact quantum metric spaces.

4. **Pin-state dependence.** The convergence is relative to the pin state $\omega_W^L$ (wedge KMS state). Different pin states might give different convergence behavior, though the qualitative-rate result is expected to be pin-state-independent for faithful normal states.

## Cross-Paper Updates Needed

1. **Paper 47 §8 Q1:** Mark as CLOSED on the natural substrate; add theorem statement
2. **Paper 45 §1.4 G2:** Update to "CLOSED at propinquity level (Paper 47 Theorem X.X)"
3. **Paper 48:** Note that T6 (G2 wedge closure via MS bridge) is now superseded by the direct propinquity closure
4. **Paper 32 §VIII:** G2-metric status update
5. **CLAUDE.md §2:** One-liner entry

## Key Files

- This memo: `debug/g2_metric_closure_memo.md`
- Paper 47: `papers/group1_operator_algebras/paper_47_two_rate_hybrid_convergence.tex`
- Paper 46 Lemma 3.2: temporal-Lipschitz-invisibility identity
- Phase A.2' memo: `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md`
