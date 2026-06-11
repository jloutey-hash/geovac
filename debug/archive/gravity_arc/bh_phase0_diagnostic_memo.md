# BH Phase 0 Diagnostic — Closure Memo

**Date:** 2026-05-22
**Sprint:** Bekenstein-Hawking probe via modular Hamiltonian, Phase 0 (diagnostic only)
**Status:** CLOSED — clean negative on naive area-law; substantive structural finding documented.
**Builds on:** Papers 42–46 (modular-Hamiltonian construction at finite cutoff), Sprint TD Track 4 (Hawking temperature reproduction via M1).
**Decision:** Do NOT proceed to Phase 1 (area-operator definition + area-coefficient extraction).
**Sprint duration:** 1 session.

---

## §1. The question and the verdict

**Question.** Does the BW canonical wedge KMS state $\rho_W = e^{-K_\alpha^W}/Z$ on the truncated Camporesi–Higuchi spectral triple of Papers 42–46 carry a Bekenstein–Hawking-style area-law signature at finite cutoff? Concretely, does $S(\rho_W) = -\mathrm{Tr}(\rho_W \log \rho_W)$ scale as $n_{\max}^2$ (matching the equator $S^2$ "area" in the truncated $S^3$)?

**Verdict.** NO. $S(\rho_W)$ scales as $\mathbf{2 \log(n_{\max})}$ at the BW canonical normalization, with log-linear fit $S \approx 1.94 \cdot \log(n_{\max}) + 0.56$ at $R^2 = 0.99991$ across $n_{\max} \in \{2, 3, 4, 5, 6, 7\}$. There is no polynomial scaling.

**Structural reading.** $S(\rho_W) \approx \log(n_{\mathrm{equator}})$ where $n_{\mathrm{equator}} = n_{\max}(n_{\max}+1)$ is the count of states with smallest $|m_j| = 1/2$. The wedge KMS state is effectively the **maximally mixed state on the lowest-$K_\alpha$ shell**, with exponentially-suppressed contributions from higher shells. The entropy is the logarithm of the boundary-shell dimension, NOT the boundary-shell dimension itself.

---

## §2. Data

### §2.1 BW canonical entropy at $n_{\max} = 2, \dots, 7$

| $n_{\max}$ | $\dim H$ | $\dim W$ | $n_{\mathrm{equator}}$ | $S(\rho_W)$ | $\log n_{\mathrm{equator}}$ | $\log \dim W$ |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 2 | 16  | 8   | 6  | 1.9222 | 1.7918 | 2.0794 |
| 3 | 40  | 20  | 12 | 2.6910 | 2.4849 | 2.9957 |
| 4 | 80  | 40  | 20 | 3.2501 | 2.9957 | 3.6889 |
| 5 | 140 | 70  | 30 | 3.6886 | 3.4012 | 4.2485 |
| 6 | 224 | 112 | 42 | 4.0491 | 3.7377 | 4.7185 |
| 7 | 336 | 168 | 56 | 4.3549 | 4.0254 | 5.1240 |

Key inequalities verified at every cell:

$$\log(n_{\mathrm{equator}}) \;\leq\; S(\rho_W) \;\leq\; \log(\dim W).$$

The two bounds are sharply realized in the $\beta$-deformation limits.

### §2.2 Asymptotic limits ($s$-deformation $\rho(s) = e^{-s K_\alpha^W}/Z$, BW canonical is $s = 1$)

| Limit | Physical reading | Predicted $S$ | Measured $S$ at $n_{\max} = 7$ |
|:--|:--|:--:|:--:|
| $s \to 0^+$ (hot) | $\rho \to I/\dim W$ maximally mixed | $\log \dim W = 5.124$ | $5.1235$ ($s = 0.01$) |
| $s = 1$ (BW canonical) | between the two limits | — | $4.3549$ |
| $s \to \infty$ (cold) | $\rho \to I_{n_{\mathrm{eq}}}/n_{\mathrm{eq}}$ on lowest shell | $\log n_{\mathrm{equator}} = 4.025$ | $4.0254$ ($s = 6.28$) |

**Both limits are saturated bit-exactly.** This confirms the structural reading: the entropy lives between two Hilbert-space-dimensional logarithms, no polynomial scaling appears anywhere in the panel.

### §2.3 Scaling fits at BW canonical

| Fit form | Parameter(s) | $R^2$ |
|:--|:--|:--:|
| $S = a \log n_{\max} + b$ | $a = 1.9435$, $b = 0.5646$ | $0.99991$ |
| $S = c \cdot n_{\max}^\alpha$ | $\alpha = 0.6480$, $c = 1.278$ | $0.98644$ |
| $S = a_2 n_{\max}^2 + a_1 n_{\max} + a_0$ | $a_2 = -0.0555$ (negative!) | $0.99886$ |
| $S = a_3 n_{\max}^3 + \dots$ | $a_3 = +0.0084$ (overfit with 4 params, 6 points) | $0.99996$ |

**The log-linear fit is the cleanest and most parsimonious.** The $n_{\max}^2$ poly fit has a *negative* leading coefficient (the function is concave, not convex), which is the polynomial-fit signature of underlying logarithmic growth. The $n_{\max}^3$ cubic R² gain is the $\Delta = 5 \times 10^{-5}$ that two extra parameters always buy on 6 data points — does not survive Occam.

### §2.4 Equator boundary fit

$n_{\mathrm{equator}}$ measured: $\{6, 12, 20, 30, 42, 56\} = \{n_{\max}(n_{\max}+1)\}$ exactly. Power-law fit $n_{\mathrm{equator}} \sim 1.71 \cdot n_{\max}^{1.78}$ at small $n_{\max}$ asymptotes to $n_{\max}^2$ for large $n_{\max}$ (the $1.78$ exponent is a finite-size effect from fitting $n(n+1)$ as a pure power).

---

## §3. Mechanism

The BW canonical wedge KMS state has eigenvalues

$$\lambda_i = \frac{e^{-\,\mathrm{two\_m}_j(i)}}{Z}, \qquad \mathrm{two\_m}_j(i) \in \{1, 3, 5, \dots, 2n_{\max} - 1\},$$

with multiplicities $g(2k+1) = (n_{\max} - k)(n_{\max} - k + 1)$ for $k = 0, 1, \dots, n_{\max} - 1$.

At $\mathrm{two\_m}_j = 1$ (the equator shell), the per-state weight is $e^{-1}/Z \approx 0.368/Z$. At $\mathrm{two\_m}_j = 3$, it drops to $e^{-3}/Z \approx 0.050/Z$ — a factor of $e^{-2} \approx 0.135$ suppression per shell.

The partition function is dominated by the equator shell: $Z \approx g(1) \cdot e^{-1} = n_{\max}(n_{\max}+1) \cdot e^{-1} \approx 0.368 \, n_{\max}^2$.

Effective probability per equator state: $p_{\mathrm{eq}} \approx 1/n_{\mathrm{equator}}$.

Hence $S \approx -n_{\mathrm{equator}} \cdot p_{\mathrm{eq}} \log p_{\mathrm{eq}} \approx \log n_{\mathrm{equator}} \approx 2 \log n_{\max}$.

**The slope $\approx 2$ is exactly $\dim(S^2_{\mathrm{equator}})$**, the dimension of the wedge boundary on $S^3$. This is a *CFT-style* log scaling along a $d=2$ boundary, not the QFT-BH area-law $A \cdot \Lambda^2 / 4$.

---

## §4. Why the naive BH expectation fails

The Phase 0 scope sketched $S(\rho_W) \sim c \cdot n_{\max}^2$ as the area-law signature, by analogy with continuum QFT where the unrenormalized entanglement entropy diverges as $A \cdot \Lambda_{\mathrm{UV}}^2 / 4$ for a 3+1D field on a 3-volume with UV cutoff $\Lambda_{\mathrm{UV}}$.

The diagnostic shows this analogy **does not transfer** to the operator-system truncation:

1. **The $K_\alpha$ spectrum is bounded** by $2n_{\max} - 1$ — there is no "UV tower" of high-$K$ modes contributing per-cell entropy.
2. **At BW canonical $\beta = 2\pi$, only the lowest few $K$-shells have significant weight.** Specifically, the equator shell ($K=1$) carries $\sim 1 - e^{-2} \approx 86\%$ of the partition function for any $n_{\max}$.
3. **The wedge state is effectively confined to a $\sim n_{\max}^2$-dimensional subspace** (the equator shell), giving entropy $\log n_{\max}^2 = 2 \log n_{\max}$.

The structural conclusion: the Paper 42 modular-Hamiltonian construction preserves the **four-witness Wick-rotation theorem** (period closure, KMS condition, $J^2$ signatures) but does **NOT preserve the BH area-law entanglement** of the continuum vacuum.

This is consistent with — and sharpens — the structural-skeleton scope statement (CLAUDE.md §2 pattern crystallization, 2026-05-07): the framework determines the algebraic skeleton of modular flow but does not autonomously generate the area-law calibration data that links wedge-entropy to area in physical units.

---

## §5. Falsification of two adjacent hypotheses

The diagnostic also closes two adjacent readings that an optimist might have hoped for:

**(F1)** "Maybe $S(\rho_W)$ scales as $n_{\max}^3$ (volume-law UV divergence)." — FALSE. The cubic poly fit has leading coefficient $0.0084 \ll 1$ and the log-linear fit is one parameter cheaper at essentially the same $R^2$. Pure log scaling.

**(F2)** "Maybe a non-BW normalization gives area-law." — FALSE. The $s$-scan across $s \in \{0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 2\pi\}$ shows every cell sits between $\log(n_{\mathrm{equator}})$ and $\log(\dim W)$. The power-law exponent $\alpha$ for $S \sim n_{\max}^\alpha$ ranges from $0.64$ (cold) to $0.72$ (hot) across all $s$, never close to $2$ or $3$. **No regime of the $s$-parameter gives polynomial scaling.**

---

## §6. What this DOES say (the substantive positive)

1. **The framework's wedge KMS state is a $\log(\text{boundary-shell dim})$-entropy state** at any finite $n_{\max}$. This is structurally cleaner than the QFT area-law would be — the entropy is bounded by a pure Hilbert-space-dimensional quantity with no calibration parameter.

2. **The entropy scaling slope ($\approx 2$) matches the dimension of the wedge boundary** ($\dim S^2 = 2$). This is a discrete operator-system analog of the Cardy-Calabrese $S = (c/3) \log L$ entropy at criticality.

3. **The boundary-state count formula $n_{\mathrm{equator}} = n_{\max}(n_{\max} + 1)$** is exact at every tested $n_{\max}$, providing a clean closed-form "boundary area" that the framework computes from first principles (no calibration input).

4. **The $s \to 0$ and $s \to \infty$ limits saturate bit-exactly to $\log(\dim W)$ and $\log(n_{\mathrm{equator}})$.** The framework respects information-theoretic bounds at both temperature extremes. This is a structural consistency check that passes cleanly.

---

## §7. Next-step options (NOT executed this sprint)

The clean Phase 0 negative on naive area-law cleanly closes the BH probe **in its current form** (wedge KMS state $\rho_W = e^{-K_\alpha^W}/Z$ on the operator-system truncation). Three follow-on routes could be considered, ordered by how much new architecture each requires:

| Route | What to test | Expected outcome | Cost |
|:--|:--|:--|:--:|
| (A) | Reduced density matrix of a global pure state (e.g., the GeoVac Dirac ground state $\ket{\Psi_0}$) on $H_W$ — different state, same wedge | Possibly different scaling; tests whether the BH analog lives at the "wrong" state | 2–4 days |
| (B) | Cigar-construction Euclidean spectral action with $S_2$ Seeley–DeWitt module — the original May-16 scoping route | $S_{BH} = A/(4G)$ with calibration $\Lambda^2 \leftrightarrow 1/G$ | 8–14 weeks (per May-16 memo) |
| (C) | Lorentzian variant of Phase 0 at $N_t > 1$ — adds temporal-grid contribution to the wedge | Same log-scaling pattern, with possible $N_t$-dependent prefactor | 1–2 days |

**Recommendation:** (A) is the cheapest concrete diagnostic that could distinguish "wedge KMS state is the wrong object" from "the operator-system truncation structurally cuts BH area-law." If (A) also returns log scaling, the Phase 0 negative extends to all wedge-reduced states and the BH-area-law-from-modular-Hamiltonian program is closed for this truncation. If (A) returns polynomial scaling, the BW canonical was the "wrong" state and (B) becomes a more attractive sprint.

**Not recommending immediate (A) execution.** The Phase 0 result is already crisp and informative; the closure is on the books. Reopen if the PI wants to invest the additional 2–4 days for the cleaner two-sided closure.

---

## §8. Files

- `debug/bh_phase0_entanglement_entropy.py` — diagnostic driver (~310 lines)
- `debug/data/bh_phase0_entanglement_entropy.json` — raw eigenvalues, entropies, fits per $(n_{\max}, s)$ cell
- `debug/data/bh_phase0_entanglement_entropy.log` — full computation log
- `debug/bh_phase0_diagnostic_memo.md` — this memo

No production `geovac/` code modified. No tests added (diagnostic-only sprint; ~310 lines of debug code).

---

## §9. Paper update

**Paper 34 §V.B** gains one row documenting the Phase 0 result as a machinery-witness off-precision entry — same flavor as the Hawking-$T$ Sprint TD Track 4 row, but on the entropy side rather than the temperature side.

Proposed row (to be applied):

| Observable | Reference | Framework | Residual | Error class | Source |
|:--|:--|:--|:--:|:--:|:--|
| Wedge KMS entanglement entropy $S(\rho_W)$ scaling at finite $n_{\max}$ | Naive QFT area-law expectation $S \sim n_{\max}^2$ | log-linear $S \approx 2\log(n_{\max})$ (R² = 0.99991, $n_{\max} = 2{-}7$) | clean negative on area-law | **C** (calibration; framework's operator-system truncation does not preserve continuum QFT vacuum entanglement) | Sprint BH-Phase0 |

**Paper 32 §VIII** (`rem:bisognano_wichmann_reading` or new `rem:bh_phase0_negative`): gains a short remark noting that the four-witness Wick-rotation theorem closure at finite $n_{\max}$ (Papers 42–46) does **not** extend to area-law BH entanglement; the structural reading is that modular-flow period closure and BH-style vacuum entanglement are independent structural features, with the framework capturing the former cleanly and not the latter.

---

## §10. Reading

The result is a clean negative on the BH-style area-law expectation **via this specific construction**. It is a positive structural finding about the entanglement content of the operator-system truncation: the wedge KMS state has $\log$-of-boundary-dim entropy, slope $\approx \dim(\partial W) = 2$. This is the Cardy-Calabrese-flavor entropy of a $d=2$ boundary, NOT the BH area-law.

The framework's modular-Hamiltonian skeleton (Papers 42–46) and the Bekenstein-Hawking area-law are *structurally distinct objects*. The four-witness Wick-rotation theorem identifies the modular-flow period; it does not autonomously generate the BH area-law calibration. This is consistent with the structural-skeleton-scope pattern documented across H1, LS-8a, HF-3/4/5, multi-focal-composition wall, and W3 calibration-data falsifications.

If the BH area-law is to be reproduced from GeoVac structure, the modular-Hamiltonian wedge KMS state at finite cutoff is **not** the right object. Route (B) — the May-16-scoping Euclidean cigar Connes-Chamseddine spectral action with one-parameter Newton-constant calibration — remains the more promising sprint (8–14 weeks), but with the honest scope statement that Newton's constant is calibration input, not derivation.
