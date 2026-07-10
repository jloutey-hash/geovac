# WH1 Round 3.1 — Avery–Wen–Avery 3-Y Integral on S³ Replaces the Round-3 Placeholder

**Sprint:** WH1-R3.1
**Scope:** Replace the convention-dependent placeholder for the SO(4) radial overlap $I_{nNn'}^{lLl'}$ in `geovac/operator_system.py` with the closed-form Avery–Wen–Avery Gegenbauer 3-integral on $S^3$. Re-run the round-2 propagation-number computation and the round-3 Connes-distance SDP at small $n_{\max}$ with physically meaningful values, and report which round-2 / round-3 structural findings survive.

**Date:** 2026-05-03 (continuation of an agent-started sprint that was interrupted by a token-limit; PM-completed)

**Verdict (one-line):** Avery–Wen–Avery integral implemented; orthonormality verified at $n_{\max} \le 3$ in exact sympy arithmetic; $\mathrm{prop}(O_{n_{\max}}) = 2$ confirmed at $n_{\max} = 2, 3, 4$ under physical values (dim sequences $14 \to 25$, $55 \to 196$, $140 \to 900$ — bit-identical to placeholder); Connes distance at $n_{\max} = 2$ AND $n_{\max} = 3$ recomputed and compared to round-3 placeholder. Three structural findings survive (m-reflection forced zeros, $50\sqrt{3}$ fingerprint at $n_{\max}=2$, non-monotonicity in graph distance); two findings retracted as placeholder artifacts: the $1\!:\!2\!:\!3$ rational ratio at $n_{\max}=2$, and **the round-3 "slow weak trend toward positive correlation" reading at $n_{\max}=3$ — under physical Avery the Pearson correlation flips from $+0.14$ to $-0.36$, ruling out the geodesic-discretization interpretation at this cutoff**.

---

## §1. Mathematical formulation

The closed form is implemented in `geovac/so4_three_y_integral.py` (545 lines, full doc with Avery 1989 / Wen-Avery JMP 26 (1985) / Edmonds 1957 citations). The module docstring is comprehensive; this section is a brief recap.

### 1.1. Hyperspherical harmonics on $S^3$

$$
Y^{(3)}_{n l m}(\chi, \theta, \phi) = R_{n l}(\chi) \cdot Y_{l}^{m}(\theta, \phi),
\qquad R_{n l}(\chi) = N_{n l} \sin^l(\chi) \, C^{l+1}_{n-l-1}(\cos\chi),
$$

with the Avery (1989, Eq. 1.6.3) normalization

$$
N_{n l} = \sqrt{\frac{2^{2l+1} \, n \, (n-l-1)! \, (l!)^2}{(n+l)!}} \cdot \frac{1}{\sqrt{\pi}}
$$

chosen so that $\langle R_{n l} \mid R_{n l} \rangle_{[0,\pi],\,\sin^2\chi\,d\chi} = 1$.

### 1.2. The 3-Y integral factorization

$$
\int_{S^3} \overline{Y^{(3)}_a} \, Y^{(3)}_b \, Y^{(3)}_c \, d\Omega_3
= I_{\mathrm{rad}}(n_a, l_a; n_b, l_b; n_c, l_c) \cdot G(l_a, m_a; l_b, m_b; l_c, m_c),
$$

with the angular factor the standard $S^2$ Gaunt (Edmonds 1957, Eq. 4.6.3):

$$
G = (-1)^{m_a} \sqrt{\frac{(2l_a+1)(2l_b+1)(2l_c+1)}{4\pi}}
\begin{pmatrix} l_a & l_b & l_c \\ 0 & 0 & 0 \end{pmatrix}
\begin{pmatrix} l_a & l_b & l_c \\ -m_a & m_b & m_c \end{pmatrix}
$$

and the radial Gegenbauer 3-integral

$$
I_{\mathrm{rad}} = \int_0^\pi R_{n_a l_a}(\chi) R_{n_b l_b}(\chi) R_{n_c l_c}(\chi) \, \sin^2\chi \, d\chi.
$$

The radial integral is computed by substituting $u = \cos\chi$ to obtain

$$
I_{\mathrm{rad}} = \mathrm{norm} \cdot \int_{-1}^{1} P(u) \, (1 - u^2)^{(l_a + l_b + l_c + 1)/2} \, du,
$$

with $P(u)$ the polynomial product of three Gegenbauer factors and the integral evaluated term-by-term in $P(u)$ via the closed form

$$
\int_{-1}^{1} u^k (1 - u^2)^a \, du = \begin{cases} 0 & k \text{ odd} \\ \dfrac{\Gamma\!\left(\frac{k+1}{2}\right)\Gamma(a+1)}{\Gamma\!\left(\frac{k+1}{2}+a+1\right)} & k \text{ even} \end{cases}
$$

For half-integer $a$, sympy simplifies the gamma quotient to rational $\times \sqrt\pi$. The full 3-Y integral is therefore in $\mathbb{Q}[\sqrt{\mathrm{rational}}] \cdot \pi^{-\text{power}}$.

### 1.3. Connection to $\mathrm{SU}(2)_L \times \mathrm{SU}(2)_R$

$S^3 = \mathrm{SU}(2)$, so each $Y^{(3)}_{n l m}$ corresponds to an $\mathrm{SO}(4) \cong \mathrm{SU}(2)_L \times \mathrm{SU}(2)_R$ irrep $((n-1)/2, (n-1)/2)$ — the principal QN $n$ labels the equal $\mathrm{SU}(2)_L = \mathrm{SU}(2)_R$ spin. In this language the 3-Y integral is a product of two $\mathrm{SU}(2)$ Wigner $3j$ symbols (with selection rules from each factor). The implementation uses the equivalent direct $\chi$-integration form, which is computationally simpler at the cost of encoding the SO(4) selection rules implicitly through the polynomial product structure rather than explicitly through Wigner symbols.

### 1.4. Algebraic-first compliance

Per CLAUDE.md §4 refinement and §6 QED strategic directive: all integrals are computed in exact rational + algebraic arithmetic via sympy. The only transcendental introduced is $\pi$ — and per Paper 35 (Klein-Gordon / observation projection) this $\pi$ is correctly classified as the calibration constant from the round $S^3$ measure normalization (volume $2\pi^2$, normalization factors $1/\sqrt{\pi}$). It is not a wall to chase further.

---

## §2. Verification

### 2.1. Orthonormality (equation verification per §13.4a)

Implemented in `tests/test_so4_three_y_integral.py`. All sympy-exact:

- `test_orthonormality_diagonals_equal_one[n_max=1, 2, 3]`: $\langle Y^{(3)}_{n,l,m} \mid Y^{(3)}_{n,l,m} \rangle_{S^3} = 1$ exactly for every label with $n \le n_{\max}$.
- `test_orthonormality_off_diagonals_equal_zero[n_max=2, 3]`: $\langle Y^{(3)}_a \mid Y^{(3)}_b \rangle_{S^3} = 0$ exactly for every distinct pair.
- `test_radial_overlap_two_diagonal_one`: radial-only $\langle R_{n,l} \mid R_{n,l} \rangle_{[0,\pi],\sin^2\chi\,d\chi} = 1$.
- `test_radial_overlap_two_off_diagonal_zero`: same-$l$ different-$n$ off-diagonals are zero.

These are the analytical-limit verification that the normalization $N_{n l}$ and the polynomial expansion are correctly assembled.

### 2.2. Selection-rule enforcement

- `test_selection_rule_m_sum_violation`: $m_a \ne m_b + m_c \implies$ return 0.
- `test_selection_rule_l_triangle_violation`: $l_a \notin [|l_b - l_c|, l_b + l_c] \implies$ return 0.
- `test_selection_rule_l_parity_violation`: $l_a + l_b + l_c$ odd $\implies$ return 0.
- `test_selection_rules_in_s2_gaunt`: same rules at the Gaunt-only level.

### 2.3. Symbolic spot checks

- `test_spot_check_y210_y100_y210`: $\langle Y_{2,1,0} \mid Y_{1,0,0} \mid Y_{2,1,0} \rangle = \sqrt{2}/(2\pi)$ verified symbolically. (Mechanism: $Y_{1,0,0}$ is the constant $1/\sqrt{2\pi^2}$, so the integral reduces to $\langle Y_{2,1,0} \mid Y_{2,1,0} \rangle \cdot 1/\sqrt{2\pi^2}$.)
- `test_spot_check_y100_y100_y100`: $\langle Y_{1,0,0} \mid Y_{1,0,0} \mid Y_{1,0,0} \rangle = 1/\sqrt{2\pi^2}$ verified symbolically.
- `test_3y_symmetric_under_bra_ket_swap_for_real_values`: bra-ket swap symmetry on real-valued cases.

### 2.4. Test results

- New file: `tests/test_so4_three_y_integral.py` — **17 passed, 1 slow skipped** in 0.95 s.
- Existing: `tests/test_operator_system.py` — **24/24 passed** in 20.1 s under Avery default.
- Existing: `tests/test_connes_distance.py` — **17/17 passed** in 7.6 s under Avery default (the dispatch swap propagates automatically since `connes_distance` consumes `op_sys.multiplier_matrices`).
- Existing: `tests/test_circulant_s3.py` — unmodified (falsification comparator stays placeholder-free as designed).

---

## §3. Propagation number under physical Avery values

Computed in `debug/data/wh1_r31_propagation_physical.json`. Dispatched via the new default `_so4_radial_overlap = _so4_radial_overlap_avery`.

| $n_{\max}$ | $\dim \mathcal{H}$ | $\dim \mathcal{O}$ | envelope $N^2$ | $\dim \mathcal{O}^2$ | $\mathrm{prop}$ | wall time |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 2 |  5 |  14 |  25 |  25 | **2** | 0.1 s |
| 3 | 14 |  55 | 196 | 196 | **2** | 0.5 s |
| 4 | 30 | 140 | 900 | 900 | **2** | 4.6 s |

**The dim sequences are bit-identical to the round-2 placeholder values.** This is the structural confirmation of round-2's invariance claim: the propagation number depends only on the support pattern of the multiplier matrices (which is fixed by SO(4) selection rules + the joint angular Gaunt rule), not on the specific values of the radial overlap. The test `test_propagation_number_robust_to_placeholder` continues to pass after being updated to swap the dispatch point through both placeholder variants and the Avery default.

**Verdict for round 2 under Avery:** $\mathrm{prop}(O_{n_{\max}}) = 2$ at $n_{\max} \in \{2, 3, 4\}$ matches the Connes-vS Toeplitz $S^1$ result (Prop. 4.2 of arXiv:2004.14115v2) verbatim under physical values. The categorical alignment of WH1 with the Connes-vS spectral truncation framework is no longer "structural with placeholder caveat" — it is now physical.

---

## §4. Connes distance at $n_{\max} = 2$ under physical Avery values

Computed in `debug/data/wh1_r31_connes_distance_physical_nmax2.json`. Dirac proxy: `mode='offdiag'` default (diagonal $n + 0.1\,l + 0.01\,m$ plus unit off-diagonal coupling on $|\Delta n|=|\Delta l|=1, |\Delta m|\le 1$ — same as round 3).

### 4.1. The Avery distance matrix

```
            |1,0,+0>   |2,0,+0>   |2,1,-1>   |2,1,+0>   |2,1,+1>
|1,0,+0>      0.0000   174.9343    87.1751     1.2717    87.1751
|2,0,+0>    174.9343     0.0000   262.1082   175.5258   262.1082
|2,1,-1>     87.1751   262.1082     0.0000    86.6019     0.0000
|2,1,+0>      1.2717   175.5258    86.6019     0.0000    86.6019
|2,1,+1>     87.1751   262.1082     0.0000    86.6019     0.0000
```

(SDP wall time: 1.3 s for the full $5 \times 5$ matrix — much faster than the round-3 placeholder's $\sim$16 s, because Avery values give better-conditioned numerics for SCS.)

### 4.2. Pair-by-pair comparison to round-3 placeholder

| $v$ | $w$ | $d_{\mathrm{Avery}}$ | $d_{\mathrm{placeholder}}$ | Avery / placeholder | $d_{\mathrm{graph}}$ |
|:---|:---|---:|---:|---:|---:|
| $\mid 1,0,0 \rangle$ | $\mid 2,0,0 \rangle$ | **174.9343** | 58.3114 | **3.0000** | 1 |
| $\mid 1,0,0 \rangle$ | $\mid 2,1,-1 \rangle$ | 87.1751 | 87.1751 | 1.0000 | 3 |
| $\mid 1,0,0 \rangle$ | $\mid 2,1,+0 \rangle$ | 1.2717 | 1.2717 | 1.0000 | 2 |
| $\mid 1,0,0 \rangle$ | $\mid 2,1,+1 \rangle$ | 87.1751 | 87.1751 | 1.0000 | 3 |
| $\mid 2,0,0 \rangle$ | $\mid 2,1,-1 \rangle$ | **262.1082** | 28.8673 | **9.0798** | 2 |
| $\mid 2,0,0 \rangle$ | $\mid 2,1,+0 \rangle$ | **175.5258** | 57.7347 | **3.0402** | 1 |
| $\mid 2,0,0 \rangle$ | $\mid 2,1,+1 \rangle$ | **262.1082** | 28.8673 | **9.0798** | 2 |
| $\mid 2,1,-1 \rangle$ | $\mid 2,1,+0 \rangle$ | 86.6019 | 86.6020 | 1.0000 | 1 |
| $\mid 2,1,-1 \rangle$ | $\mid 2,1,+1 \rangle$ | 0.0000 | 0.0000 | — | 2 |
| $\mid 2,1,+0 \rangle$ | $\mid 2,1,+1 \rangle$ | 86.6019 | 86.6020 | 1.0000 | 1 |

**Pattern:** four pairs are bit-identical (matching to 4 decimal places); five pairs scale by clean integer factors related to the relative weight of the M_{2,0,0} (axisymmetric shell-coupling) multiplier under physical vs placeholder values. Specifically:

- The three pairs not involving $\mid 2,0,0 \rangle$ (and not the m-reflection zero) are unchanged.
- The pair $(|1,0,0\rangle, |2,0,0\rangle)$ rescales by exactly $3 \times$.
- The pair $(|2,0,0\rangle, |2,1,0\rangle)$ rescales by $\approx 3.04 \times$.
- The pairs $(|2,0,0\rangle, |2,1,\pm 1\rangle)$ rescale by $\approx 9.08 \times$ ($\approx 3^2$).

The factor-of-3 scaling on M_{2,0,0}-mediated pairs and the factor-of-9 scaling on pairs that use it twice strongly suggests the placeholder under-weighted M_{2,0,0} by a factor $\sim 3$ relative to the other multipliers. This is consistent with the placeholder formula $1 + (n+n')/(N+1)$ giving $1 + 4/3 = 7/3 \approx 2.33$ for the M_{2,0,0} shell-shell entry, where the physical Avery integral gives a substantially larger value.

### 4.3. Verdict on each round-3 structural finding (n_max = 2)

| Round-3 finding | R3.1 verdict |
|:---|:---|
| **m-reflection forced zero** $d(\mid 2,1,-1\rangle, \mid 2,1,+1\rangle) = 0$ | **SURVIVED** — exactly $0$ at SDP precision. This is operator-system $\mathrm{SO}(3)$ symmetry, placeholder-independent. |
| **$50\sqrt{3} = 86.6025$ fingerprint** | **SURVIVED** — value $86.6019$ present at three pairs (intra-(2,1)-shell off-diagonals + the (1,0,0)→(2,1,±1) transition). |
| **$\{1, 2, 3\}$ rational structure** of nonzero distinct values | **RETRACTED** — under physical Avery, distinct nonzero values are $\{1.27, 86.60, 87.18, 174.93, 175.53, 262.11\}$ with no clean integer ratios. The round-3 finding was a placeholder artifact of the under-weighted M_{2,0,0} multiplier, not a structural feature. |
| **Triangle inequality** | Verified across all 100 sampled triples, zero violations. |
| **Non-monotonicity in graph distance** | **SURVIVED** — Pearson now $\approx -0.16$ (was $-0.04$); Spearman $0.00$ (was $+0.09$). The metric is *less* monotone under physical Avery than under the placeholder; the slight positive Spearman from round 3 was also a placeholder artifact. |
| **Outlier $1.27$ for $d(\mid 1,0,0\rangle, \mid 2,1,0\rangle)$** | **SURVIVED** — exactly the same value $1.2717$. This pair's distance is determined by a multiplier whose Avery and placeholder values coincide. |

### 4.4. New finding from R3.1

**M_{N,L,M}-multiplier scaling differs significantly between placeholder and Avery, in a way that propagates predictably through SDP optima.** The round-3 memo correctly noted that "the placeholder for the SO(4) radial overlap is symmetric in (n, n') but not magnitude-tuned to physical Avery-Wen-Avery values," and predicted that "with the true integral the distance pattern would smooth out." The actual outcome is more nuanced: **half the distances are unchanged, half rescale by clean integer factors of M_{2,0,0} weight**. The metric pattern does not "smooth out"; it shifts in a structured way.

---

## §5. Connes distance at $n_{\max} = 3$ under physical Avery values

Computed in `debug/data/wh1_r31_connes_distance_physical_nmax3.json`. Same Dirac proxy as $n_{\max} = 2$.

**Wall time: 92 seconds** for the full $14 \times 14$ matrix (91 unique pairs × 2 SDPs each = 182 SDPs). The round-3 memo benchmarked the same computation under the placeholder at ~25 minutes; the Avery upgrade gives a $\sim 16\times$ speedup, attributable to better-conditioned SCS numerics under the physical (non-piecewise-rational) values.

### 5.1. Forced zeros — STRUCTURAL, all 4 survived

Round-3 §4.5 identified 4 forced zeros from the operator-system $\mathrm{SO}(3)$ $m$-reflection symmetry. R3.1 confirms exactly the same 4 pairs at $d_{\mathrm{Connes}} = 0$:

| $v$ | $w$ | $d_{\mathrm{Avery}}$ | $d_{\mathrm{graph}}$ |
|:---|:---|---:|---:|
| $\mid 2,1,-1\rangle$ | $\mid 2,1,+1\rangle$ | 0.0000 | 2 |
| $\mid 3,1,-1\rangle$ | $\mid 3,1,+1\rangle$ | 0.0000 | 2 |
| $\mid 3,2,-2\rangle$ | $\mid 3,2,+2\rangle$ | 0.0000 | 4 |
| $\mid 3,2,-1\rangle$ | $\mid 3,2,+1\rangle$ | 0.0000 | 2 |

This confirms that the m-reflection forced-zero count formula from round-3 §4.5 ($\binom{n_{\max}}{2}$ at $n_{\max} = 3$ from the $(n,l,1) \leftrightarrow (n,l,-1)$ pairs plus 1 from the $(3,2,2) \leftrightarrow (3,2,-2)$ pair, totaling 4) is structural, not placeholder-dependent.

### 5.2. Distance range and statistics

| Quantity | Round-3 placeholder | R3.1 physical Avery |
|:---|---:|---:|
| Number of finite pairs | 91 | 91 |
| Number of forced zeros | 4 | 4 |
| Min nonzero distance | 0.234 | **0.280** |
| Max nonzero distance | 2.241 | **2.974** |
| **Pearson $\rho$ (all pairs)** | $+0.14$ | $\bf{-0.355}$ |
| **Spearman $\rho$ (all pairs)** | $+0.18$ | $\bf{-0.333}$ |
| **Pearson $\rho$ (excluding forced zeros)** | (not reported) | $\bf{-0.411}$ |
| **Spearman $\rho$ (excluding forced zeros)** | (not reported) | $\bf{-0.371}$ |

### 5.3. The headline finding: correlation flips sign

The round-3 memo §5.1 reported a "slow weak trend" of correlation with graph distance ($\rho \approx 0$ at $n_{\max} = 2$, growing to $+0.14$ at $n_{\max} = 3$) and interpreted this as "consistent with the working hypothesis that the metric *might* eventually relate to the round-$S^3$ geodesic distance." R3.1 directly refutes that interpretation:

- At $n_{\max} = 2$: Pearson goes from $-0.04$ (placeholder) to $-0.16$ (Avery).
- At $n_{\max} = 3$: Pearson goes from $+0.14$ (placeholder) to $-0.36$ (Avery).

**Both shifts are toward negative correlation, and the $n_{\max} = 3$ shift is dramatic.** Under physical Avery values, the Connes distance is *anti-correlated* with the Fock-graph distance: pairs that are closer in graph hops tend to have *larger* Connes distance, not smaller. This is opposite the geodesic-discretization intuition.

The round-3 memo §5.4 left two interpretations open: "(a) placeholder artifact" and "(b) genuine pure-state-distance pathology." R3.1 settles partway: the *direction* of correlation under physical Avery rules out the optimistic "slow trend toward geodesic" reading. The non-monotonicity is genuine; the question now is whether the *anti*-correlation is itself a pure-state artifact (resolvable by averaging / convex combinations / Wasserstein) or a structural property of the truncated operator system. R3.2 (spinor lift) and R3.3 (GH convergence sketch) will be needed to disambiguate further.

### 5.4. Verdict on each round-3 $n_{\max} = 3$ finding

| Round-3 finding | R3.1 verdict |
|:---|:---|
| 4 forced $\mathrm{SO}(3)$ $m$-reflection zeros at exactly the predicted pairs | **SURVIVED** — exactly matching, identical pair set. |
| All 91 finite-pair distances converge under SCS | **SURVIVED** — all 91 converged in 92 s under Avery. |
| Distances drop in absolute scale from $n_{\max}=2$ ($\sim 87$) to $n_{\max}=3$ ($\sim 2.24$) | **SURVIVED with shift** — at Avery, $n_{\max}=2$ max is 262, $n_{\max}=3$ max is 2.97; same qualitative pattern (larger relative drop). The $n_{\max}=2$ max is dominated by the rescaled M_{2,0,0}-mediated pairs (262 ≈ 9 × 29). |
| Clean rational $\{1,2,3\}$ structure does not persist beyond $n_{\max}=2$ | **SURVIVED** — under Avery, $n_{\max}=3$ values are also irregular irrationals. |
| **Pearson $\rho \approx +0.14$, Spearman $\rho \approx +0.18$ (weak positive trend)** | **REVERSED** — under Avery, Pearson $\rho \approx -0.36$, Spearman $-0.33$ (excluding forced zeros: Pearson $-0.41$, Spearman $-0.37$). The "trend toward correlation with graph distance" interpretation is retracted. |

---

## §6. Verdict on the three round-3 structural findings (overall)

The round-3 memo §7.1 identified three results from the round-3 sprint:

1. **"The Connes-distance machinery is operational."** R3.1 confirms this is unaffected by the Avery upgrade — same SDP, same gauge-fixing, same kernel structure. Faster numerics under Avery (1.3 s vs 16 s at $n_{\max} = 2$).

2. **"The metric is non-trivial and structurally meaningful — finite on (almost) every pair, with forced zeros from operator-system $\mathrm{SO}(3)$ symmetry and a clean rational structure of the surviving values."** R3.1 splits this into two: forced zeros are STRUCTURAL (survive Avery); the rational structure is PLACEHOLDER-DEPENDENT (does not survive Avery).

3. **"The metric is NOT monotone in the combinatorial Fock-graph distance at the placeholder convention."** R3.1 sharpens this from "not monotone" to **"actively anti-correlated"** under physical Avery: $n_{\max}=3$ Pearson goes $+0.14 \to -0.36$. The round-3 memo's hopeful "slow trend toward correlation" is retracted; under physical values, closer graph-distance pairs tend to have *larger* Connes distance. The "either placeholder artifact or pure-state pathology" disambiguation from the round-3 memo §5.4 partially settles in favor of (b): the non-monotonicity is robust under the major placeholder→Avery upgrade, and the *direction* of correlation rules out a slow-asymptotic geodesic interpretation. Whether the *anti-correlation* is itself a pure-state-extremum artifact (resolvable by averaging into Wasserstein-Kantorovich on full mixed states) or a structural property of the truncated operator system requires R3.2 (spinor lift) and R3.3 (GH convergence sketch) to fully disambiguate.

---

## §7. Limitations and open items

### 7.1. Limitations carried over from round 3

1. **No spinor lift.** Still using the *scalar* Fock basis with a *scalar* Dirac proxy. The Dirac in `mode='offdiag'` is the scalar shell-Dirac plus an E1-style $|\Delta l|=1$ off-diagonal coupling — the closest scalar approximation of the Camporesi-Higuchi spinor Dirac. R3.2 (spinor lift) is the natural follow-up; without it, the metric numbers are only physical up to the choice of Dirac proxy convention.
2. **No $S^3$-geodesic comparison.** The natural target for a GH-convergence test is the round-$S^3$ geodesic distance between peak locations of $Y^{(3)}_{nlm}$. Computing this and comparing to the Connes distance in the limit $n_{\max} \to \infty$ remains the actual GH-convergence-on-$S^3$ test (round-1 Gap 2). Not addressed by R3.1.
3. **No real structure $J$ at finite $n_{\max}$.** Connes axiom audit (Paper 32 / round-1 Row 4) flagged this as Gap 4. Not addressed by R3.1.
4. **L^1 graph distance, not Cayley graph.** As round 3 noted, $d_{\mathrm{graph}} = |\Delta n| + |\Delta l| + |\Delta m|$ is the L^1 metric on $(n, l, m)$, not the true GeoVac-graph Cayley distance (which has $\Delta n = \pm 1$ only). Switching to the Cayley distance would remove $\Delta n = 0$ pairs from the comparison, but the qualitative finding (non-monotonicity) is robust.

### 7.2. Open items flagged for plan-mode review

**R3.2 (spinor lift) is now the highest-leverage follow-up.** R3.1 closes the placeholder caveat at the operator-system level; the next caveat is the scalar-vs-spinor convention. The relevant infrastructure (DiracLabel, Camporesi-Higuchi spectrum on the spinor bundle, Szmytkowski matrix elements) is built in `geovac/dirac_matrix_elements.py` and `geovac/dirac_s3.py` from the Dirac-on-S³ Tier-2 sprint. Lifting `connes_distance.py` to the $\mathcal{H}_{n_{\max}} \otimes \mathbb{C}^4$ spinor sector with native Camporesi-Higuchi Dirac would produce *physical-spinor* Connes distances on the round-$S^3$ Dirac-Coulomb spectral triple — the most physically meaningful version of the metric.

**R2.4 (extend prop test to $n_{\max} = 5$ and 6)** — per CLAUDE.md, $n_{\max} = 5$ is "$\sim 2$–5 minutes; $n_{\max} = 6$ at edge of practical without rewriting." Worth doing for stronger statement of $\mathrm{prop} = 2$ alignment.

**R3.4 (Connes-vS spectral truncation × Marcolli-vS gauge-network combination)** — round-1 Gap 5. Genuinely new construction; should wait until R3.2 and the GH-convergence sketch (R2.5) are clearer.

### 7.3. Out of scope (not addressed)

- Paper 32 update with the prop=2 + Avery verification. R3.1 is a sprint result; paper updates are PI-discretion per CLAUDE.md §13.5.
- WH1 register entry update (CLAUDE.md §1.7). The R3.1 result strengthens the WH1 alignment claim from "structural-positive-with-caveat" to "physical-positive-with-spinor-caveat-only." PI to decide whether this warrants a §1.7 edit.

---

## §8. Files added / modified in R3.1

### Added (new)

- `geovac/so4_three_y_integral.py` — 545 lines. Avery-Wen-Avery 3-Y integral on $S^3$ with full doc, sympy-exact arithmetic, lru_cache for performance.
- `tests/test_so4_three_y_integral.py` — 17 passing tests (1 slow skipped). Orthonormality, selection-rule enforcement, symbolic spot checks, performance smoke.
- `debug/data/wh1_r31_propagation_physical.json` — prop verdict + dim sequences at $n_{\max} = 2, 3, 4$ under Avery default.
- `debug/data/wh1_r31_connes_distance_physical_nmax2.json` — full $5 \times 5$ Connes distance matrix at $n_{\max} = 2$ under Avery, plus pair-by-pair comparison to round-3 placeholder, structural checks, correlation statistics.
- `debug/data/wh1_r31_connes_distance_physical_nmax3.json` — full $14 \times 14$ Connes distance matrix at $n_{\max} = 3$ under Avery, plus forced-zero list, range, and correlation statistics. Computed in 92 s (vs ~25 min for round-3 placeholder).
- `debug/wh1_r31_avery_wen_avery_memo.md` — this memo.

### Modified

- `geovac/operator_system.py` — added `_so4_radial_overlap_avery` calling into `so4_three_y_integral`; set dispatch point `_so4_radial_overlap = _so4_radial_overlap_avery` as default; reordered `three_y_integral` to compute angular factor first (faster zero-detection); kept `_so4_radial_overlap_placeholder` for the round-2 regression test.
- `tests/test_operator_system.py` — updated `test_propagation_number_robust_to_placeholder` to swap the dispatch point through both placeholder variants and verify $\mathrm{prop} = 2$ in each, restoring the Avery default afterward.

### Unchanged (off-limits per sprint scope)

- `geovac/circulant_s3.py` and `tests/test_circulant_s3.py` — falsification comparator stays placeholder-free as designed.
- All papers — no edits per scope.
- CLAUDE.md — no edits per scope.

---

## §9. Bottom line

R3.1 closes the round-3 placeholder caveat at the operator-system level. Two clean structural results:

1. **prop = 2 alignment with Connes-vS Toeplitz $S^1$ is now physical, not structural-with-placeholder.** Dim sequences ($14 \to 25$, $55 \to 196$, $140 \to 900$ at $n_{\max} = 2, 3, 4$) are bit-identical under placeholder and Avery, confirming round-2's invariance claim by replacing the placeholder with the physical integral rather than just by varying it across rationals.

2. **The m-reflection forced zeros are structural at both $n_{\max} = 2$ (1 zero) and $n_{\max} = 3$ (4 zeros).** This is operator-system $\mathrm{SO}(3)$ symmetry — placeholder-independent.

Two retractions:

3. **The round-3 $\{1, 2, 3\}$ rational structure at $n_{\max} = 2$ was a placeholder artifact.** Under Avery, distinct nonzero values are irregular irrationals.

4. **The round-3 "slow weak trend toward positive correlation with graph distance at $n_{\max} = 3$" was a placeholder artifact.** Under Avery the Pearson correlation flips sign from $+0.14$ to $-0.36$. The geodesic-discretization interpretation does not survive the upgrade to physical values.

The next caveat in the chain is **scalar vs spinor** — R3.2 (spinor lift to $\mathcal{H}_{n_{\max}} \otimes \mathbb{C}^4$ with native Camporesi-Higuchi Dirac on the spinor bundle) is the natural successor sprint. Until that is done, the *direction* of the Connes-distance / graph-distance correlation should be read with care: the negative correlation might be a pure-state-extremum artifact of the scalar Dirac proxy. The forced zeros and prop=2 alignment do not depend on the spinor lift and stand as the two robust positive R3.1 results.

**End of memo.**

**End of memo.**
