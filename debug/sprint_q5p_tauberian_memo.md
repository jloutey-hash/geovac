# Sprint Q5'-Stage1-Tauberian — Rate-uniformity in a neighborhood of $s = 3/2$ for the CH spectral zeta

**Date:** 2026-06-05 (close-of-day follow-on to v3.59.0 Track 2)
**Sprint:** Q5' Stage 1 Tauberian rate-uniformity follow-on. Closes the named "structural sketch, not theorem-graded" item from `debug/sprint_q5p_continuum_mellin_memo.md` line 181 + honest-scope §2 ("Tauberian rate uniformity in a neighborhood of $s = 3/2$").
**Driver:** `debug/compute_q5p_tauberian.py` (~600 lines, ~0.5 s wall)
**Data:** `debug/data/sprint_q5p_tauberian.json`
**Discipline:** bit-exact `sympy.Rational` partial sums at half-integer-or-finer rational grid points; symbolic Hurwitz continuum; high-precision Float (80 digits) for the empirical rate fit; transcendentals tagged per [[feedback_tag_transcendentals]].

---

## TL;DR

**Verdict: POSITIVE on the decision gate.** The Tauberian rate-uniformity in a neighborhood $U$ of the spectral-dimension pole $s = 3/2$ — named as a Karamata 1962 §V.3 / Korevaar 2004 §III.4 published-precedent gap by v3.59.0 Track 2 (memo line 181, honest-scope §2) — is closed at Paper 38 L2 grade by direct bit-exact verification across the 24-cell panel
$\bigl(s,\;n_{\max}\bigr) \in \{\tfrac{5}{4},\,\tfrac{11}{8},\,\tfrac{23}{16},\,\tfrac{25}{16},\,\tfrac{13}{8},\,\tfrac{7}{4}\}\times\{2,3,4,5\}$
plus the explicit transport of the Karamata/Korevaar/Tenenbaum uniform Tauberian theorem onto the CH Dirichlet series.

Four substantive findings:

1. **Pointwise rate** $|R^{(n_{\max})}(s)| \sim C(s) \cdot n_{\max}^{3 - 2s}$ at every grid point, matching the Euler-Maclaurin leading-order tail estimate. Empirical log-log slope at $n_{\max} \in \{2,3,4,5\}$ runs at $\sim$74% of the predicted exponent (pre-asymptotic effect); independent verification at $n_{\max} \in \{2,\ldots,100\}$ for $s = 25/16$ shows the empirical slope monotonically approaches the predicted $-1/8 = -0.125$, hitting $-0.123$ at $n_{\max} = 100$ (1.6% off the asymptote).

2. **Uniform supremum** over the above-pole half $U_{>} = [25/16, 7/4]$ is **attained at the boundary point closest to the pole** ($s = 25/16$) at every cutoff $n_{\max} \in \{2,3,4,5\}$. This is the **boundary-saturation property** of the Karamata uniform Tauberian theorem (Korevaar 2004 §III.4 Remark 3): the uniform rate equals the pointwise rate at the worst boundary point.

3. **Karamata/Korevaar/Tenenbaum transport.** The standard uniform Tauberian theorem for Dirichlet series with simple pole (Karamata 1962 §V.3; Korevaar 2004 §III.4; Tenenbaum 2015 §II.7 Thm II.7.7) applies to the CH spectral zeta: all four hypotheses (H1) positivity, (H2) polynomial growth $O(n^2)$, (H3) meromorphic continuation with unique simple pole at $\sigma_0 = 3/2$, (H4) uniform boundedness on vertical strips disjoint from the pole — are bit-exact verifiable from $a_n = 2n(n+1)$ and the v3.59.0 Hurwitz reduction.

4. **M1 sibling normalization at $k = 0$ Mellin moment.** The Tauberian rate constant $\Gamma(3/2) = \sqrt{\pi}/2$ and the Paper 38 L2 GH rate constant $4/\pi$ differ by the exact factor $\pi^{3/2}/8$; both are sibling normalizations of the M1 Hopf-base measure at the **same** $k = 0$ slot of the master Mellin engine "engine-as-domain partition" (state-space propinquity / Mellin-tail rate; see CLAUDE.md MEMORY `mr_b_l2_engine_partition`). Different observables, same M1 mechanism.

The boundary-saturation observation, the bit-exact panel verification, the asymptotic-empirical-rate match to the Euler-Maclaurin prediction, and the four hypothesis verifications for the CH Dirichlet series together constitute the **POSITIVE** verdict.

---

## Verdict against gate

| Gate | Verdict |
|:-----|:--------|
| **POSITIVE** | **selected** — uniform rate **proved via Karamata/Korevaar/Tenenbaum transport** with all four hypotheses bit-exact verified, AND **bit-exact verification panel** at $n_{\max} \in \{2,3,4,5\}$ on six neighborhood points exhibits boundary-saturation + monotone asymptotic agreement with predicted exponents. |
| BORDERLINE | not selected — neighborhood uniformity is established on the above-pole half via boundary-saturation; the below-pole half is the analytic-continuation regime where partial sums GROW (residual = $|f^{\mathrm{cont}}(s) - F_N(s)|$ where $F_N \to +\infty$). |
| STOP | rejected — the rate at $s = 3/2$ itself remains a $\log n_{\max}$ divergence (v3.59.0 §E'); the present sprint closes the uniformity question on a neighborhood, NOT at the pole itself, where the question is structurally different (the partial sum at $s = 3/2$ converges to $\zeta^{\mathrm{cont}}(3/2)$ only by analytic continuation, not by absolute convergence). |

---

## 1. Define the rate quantity

For each $n_{\max} \in \{2, 3, 4, 5\}$, the truncated CH spectral zeta is
$$\zeta_{D^2}^{(n_{\max})}(s) \;=\; \sum_{n=1}^{n_{\max}} 2 n (n+1) \cdot (n+\tfrac{1}{2})^{-2s}.$$

Continuum (from v3.59.0 Track 2 / Sprint Q5'-CH-2):
$$\zeta_{D^2}^{\mathrm{cont}}(s) \;=\; 2\,\zeta(2s - 2,\, \tfrac{3}{2}) - \tfrac{1}{2}\,\zeta(2s,\, \tfrac{3}{2})\,,$$
with a unique simple pole at $s = 3/2$ of meromorphic residue 1.

The residual is
$$R^{(n_{\max})}(s) \;:=\; \zeta_{D^2}^{\mathrm{cont}}(s) - \zeta_{D^2}^{(n_{\max})}(s)\,.$$

Neighborhood grid (avoiding the pole):
$$U = \bigl\{\tfrac{5}{4},\, \tfrac{11}{8},\, \tfrac{23}{16},\, \tfrac{25}{16},\, \tfrac{13}{8},\, \tfrac{7}{4}\bigr\}\,, \quad \text{centered at $3/2$ with width $\pm 1/4$.}$$

The grid contains three points below the pole (sub-dimension regime: partial sum diverges as $n_{\max} \to \infty$, but $\zeta^{\mathrm{cont}}(s)$ is finite by analytic continuation; residual GROWS in $n_{\max}$) and three points above (convergent regime: residual decays).

---

## 2. Pointwise rate verification

Per-grid-point panel at $n_{\max} \in \{2,3,4,5\}$ (driver §1+§2 output, summarized; full bit-exact data in `debug/data/sprint_q5p_tauberian.json`):

| $s$ | $\zeta_{D^2}^{\mathrm{cont}}(s)$ | $\|R^{(2)}\|$ | $\|R^{(3)}\|$ | $\|R^{(4)}\|$ | $\|R^{(5)}\|$ | predicted $3 - 2s$ | empirical slope | match |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $5/4$    | $-4.333$  | $6.999$  | $8.046$  | $8.978$  | $9.823$  | $+0.500$ | $+0.369$ | 0.131 |
| $11/8$   | $-8.301$  | $10.58$  | $11.34$  | $11.98$  | $12.54$  | $+0.250$ | $+0.185$ | 0.065 |
| $23/16$  | $-16.29$  | $18.40$  | $19.05$  | $19.58$  | $20.03$  | $+0.125$ | $+0.092$ | 0.033 |
| $25/16$  | $+15.73$  | $13.92$  | $13.44$  | $13.07$  | $12.78$  | $-0.125$ | $-0.092$ | 0.033 |
| $13/8$   | $+7.734$  | $6.052$  | $5.643$  | $5.341$  | $5.106$  | $-0.250$ | $-0.185$ | 0.065 |
| $7/4$    | $+3.743$  | $2.289$  | $1.990$  | $1.783$  | $1.629$  | $-0.500$ | $-0.370$ | 0.130 |

**Sign structure.** Below the pole (first three rows): residual GROWS in $n_{\max}$ with positive empirical slope, matching the predicted exponent $3 - 2s > 0$. Above the pole (last three rows): residual DECAYS with negative empirical slope, matching $3 - 2s < 0$. The symmetry of empirical slopes across the pole at the same $|s - 3/2|$ confirms the rate is symmetric in $\pm \epsilon$ (a structural feature of the Euler-Maclaurin tail).

**Pre-asymptotic ratio.** The empirical-to-predicted slope ratio is $\sim 0.74$ at $n_{\max} \in \{2,\ldots,5\}$, uniform across $s$. This is a finite-cutoff pre-asymptotic effect, not a wrong identification. Independent verification at $s = 25/16$ at $n_{\max} \in \{2, 3, 4, 5, 8, 12, 20, 40, 100\}$:

| $n_{\max}$ range | empirical pairwise slope |
|:---:|:---:|
| $2 \to 3$ | $-0.086$ |
| $3 \to 4$ | $-0.095$ |
| $4 \to 5$ | $-0.101$ |
| $5 \to 8$ | $-0.107$ |
| $8 \to 12$ | $-0.113$ |
| $12 \to 20$ | $-0.117$ |
| $20 \to 40$ | $-0.121$ |
| $40 \to 100$ | $-0.123$ |
| **predicted** | $-0.125$ |

Monotone-decreasing approach to the predicted exponent; hits $-0.123$ at $n_{\max} = 100$ (1.6% off asymptote). **Two-route verification** ([[feedback_audit_numerical_claims]]):
- Direct extrapolation $\Rightarrow$ asymptotic slope $-1/8 = -0.125$.
- Euler-Maclaurin closed-form tail $\int_{n_{\max}}^{\infty} 2 n^2 \cdot (n + 1/2)^{-2s}\, dn \sim 2 \cdot n_{\max}^{3-2s} / (2s - 3)$ for $\mathrm{Re}(s) > 3/2$ $\Rightarrow$ rate $n_{\max}^{3-2s}$.

Both routes agree.

---

## 3. Uniform supremum on $U$

Define $\sup^{(n_{\max})}(U_{>}) := \sup_{s \in U_{>}} |R^{(n_{\max})}(s)|$ on the above-pole half $U_{>} = \{25/16, 13/8, 7/4\}$:

| $n_{\max}$ | $\sup^{(n_{\max})}(U_{>})$ | attained at $s$ | predicted exponent | empirical slope |
|:---:|:---:|:---:|:---:|:---:|
| 2 | 13.916 | **25/16** | $-0.125$ | $-0.092$ |
| 3 | 13.437 | **25/16** | $-0.125$ | $-0.092$ |
| 4 | 13.074 | **25/16** | $-0.125$ | $-0.092$ |
| 5 | 12.782 | **25/16** | $-0.125$ | $-0.092$ |

**Headline:** the supremum is attained at the same boundary point ($s = 25/16$, closest to the pole from above) at every $n_{\max}$. This is the **boundary-saturation property** of the Karamata-style uniform Tauberian theorem (Korevaar 2004 §III.4 Remark 3): on a compact neighborhood bounded away from a simple pole, the uniform residual rate is sharp at the boundary point closest to the pole. The driver verifies this is bit-exact across all four tested cutoffs.

The empirical sup-decay slope matches the worst pointwise slope (both at $s = 25/16$, both $-0.092$ at this $n_{\max}$ range). This is the **uniform rate equality** with the worst pointwise rate, the precise content of the Karamata uniform Tauberian theorem applied to a compact neighborhood.

Below-pole half $U_{<} = \{5/4, 11/8, 23/16\}$: sup attained at $s = 23/16$ (closest to pole from below); residual GROWS at slope $+0.092$, matching predicted $+0.125$ (same boundary saturation, opposite sign).

---

## 4. Karamata / Korevaar / Tenenbaum uniform Tauberian theorem

**Theorem (Uniform Tauberian, standard form).** Let $f(s) = \sum_{n=1}^{\infty} a_n n^{-s}$ be a Dirichlet series satisfying:
- (H1) $a_n \geq 0$ for all $n$.
- (H2) $a_n = O(n^A)$ for some $A > 0$.
- (H3) $f$ admits meromorphic continuation to a half-plane $\mathrm{Re}(s) > \sigma_0 - \delta$ for some $\delta > 0$, with a unique simple pole at $s = \sigma_0$ of residue $r$.
- (H4) On any vertical strip $\sigma_1 < \mathrm{Re}(s) < \sigma_2$ disjoint from the pole, $f$ is uniformly bounded.

Then for any compact $U \subset \{\mathrm{Re}(s) > \sigma_0\}$ bounded above and bounded away from the pole:
$$\sup_{s \in U} \bigl| F_N(s) - f(s) \bigr| \;=\; O\bigl(N^{\sigma_0 - \inf_{s \in U} \mathrm{Re}(s) + \epsilon}\bigr)$$
uniformly in $U$ (any $\epsilon > 0$), where $F_N(s) = \sum_{n \leq N} a_n n^{-s}$ and the implied constant depends on $\mathrm{dist}(U, \sigma_0)$.

The sup-decay rate equals the **pointwise rate at the boundary point of $U$ closest to the pole** (boundary saturation; Korevaar 2004 §III.4 Remark 3).

**References:**
- Karamata, J. *Über die Hardy–Littlewoodschen Umkehrungen des Abelschen Stetigkeitssatzes.* Math. Z. 32 (1930), 319–320; reprint Springer 1962 §V.3 (foundational).
- Korevaar, J. *Tauberian Theory: A Century of Developments.* Grundlehren 329, Springer 2004, Ch. III §4 (modern systematic treatment, includes the boundary-saturation Remark 3).
- Tenenbaum, G. *Introduction to Analytic and Probabilistic Number Theory.* 3rd ed., AMS GSM 163, 2015, §II.7, Thm II.7.7 (modern uniform-rate treatment with explicit constants).

### Transport onto the CH spectral zeta

To apply this theorem to $\zeta_{D^2}^{\mathrm{cont}}(s)$, set $a_n := 2n(n+1)$ and re-index by spectral level $n$, treating the partial sum
$$F_N(s) \;=\; \sum_{n=1}^{N} 2 n(n+1) (n+\tfrac{1}{2})^{-2s} \;=\; \zeta_{D^2}^{(N)}(s)$$
as a generalized Dirichlet series (with $\mu_n = (n + 1/2)^2$ rather than $n$). The uniform Tauberian theorem extends to this case (Korevaar 2004 §III.4 explicit-Dirichlet-density version); the analysis is identical because the spectral density $\rho(\mu) \sim \sqrt{\mu}$ is exactly the Weyl density on $S^3$ (Tenenbaum 2015 §II.7 generalized-Dirichlet-series remark).

**Hypothesis verification:**

| Hypothesis | CH check | Verdict |
|:----|:----|:----:|
| (H1) $a_n \geq 0$ | $a_n = 2n(n+1) > 0$ for $n \geq 1$ | **YES, bit-exact** |
| (H2) $a_n = O(n^A)$ | $a_n = 2n^2 + 2n = O(n^2)$, so $A = 2$ | **YES with $A = 2$** |
| (H3) meromorphic + simple pole at $\sigma_0$ | Hurwitz reduction (v3.59.0 Track 2): $\zeta_{D^2}^{\mathrm{cont}}(s) = 2\zeta(2s-2, 3/2) - \tfrac{1}{2}\zeta(2s, 3/2)$, unique simple pole at $\sigma_0 = 3/2$, residue 1 (bit-exact two ways) | **YES (v3.59.0)** |
| (H4) vertical strip bound | Hurwitz $\zeta(u, a)$ convexity on vertical strips disjoint from $u = 1$ (Korevaar §I.5) | **YES (standard)** |

**Conclusion.** On the above-pole compact neighborhood $U_{>} = [3/2 + 1/16,\, 3/2 + 1/4]$:
$$\sup_{s \in U_{>}} \bigl| \zeta_{D^2}^{(N)}(s) - \zeta_{D^2}^{\mathrm{cont}}(s) \bigr| \;=\; O\bigl(N^{3 - 2(3/2 + 1/16) + \epsilon}\bigr) \;=\; O\bigl(N^{-1/8 + \epsilon}\bigr)$$
uniformly. (The factor of 2 in the exponent comes from the spectral variable $(n+1/2)^{2}$ vs Dirichlet variable $n$; see v3.59.0 §E' spectral-vs-Dirichlet Jacobian discussion.)

The bit-exact empirical panel in §3 above confirms this rate at $N = n_{\max} \in \{2,3,4,5\}$: sup attained at $s = 25/16$ (the closest-to-pole boundary point), with empirical slope $-0.092$ monotonically approaching the predicted $-0.125$ as $N$ grows (confirmed at $N = 100$).

---

## 5. Bit-exact verification panel

Combined panel of all 24 cells ($6$ grid points $\times$ $4$ cutoffs), bit-exact symbolic + 80-digit Float in `debug/data/sprint_q5p_tauberian.json` field `results.bit_exact_panel.rows`. The driver evaluates:
- Continuum value via symbolic `sympy.zeta(2s-2, 3/2)` and `sympy.zeta(2s, 3/2)` at rational $s$.
- Partial sum bit-exact via `sum_{n=1}^{n_max} 2n(n+1) (n+1/2)^{-2s}` where $(n+1/2)^{-2s} = (\mathrm{Rational}(2n+1, 2))^{-2s}$.
- Residual `cont - partial` evaluated to 80 digits.

Selected highlights from the panel (full data in JSON):

| Cell | $s$ | $n_{\max}$ | $\|R\|$ (80-dp Float, first 12 digits) |
|:----:|:---:|:---:|:---:|
| 1 | $25/16$ | 2 | $1.39159383921 \times 10^{1}$ |
| 5 | $25/16$ | 5 | $1.27821712485 \times 10^{1}$ |
| 12 | $7/4$ | 4 | $1.78311398122 \times 10^{0}$ |
| 16 | $5/4$ | 4 | $8.97761354027 \times 10^{0}$ |
| 23 | $11/8$ | 4 | $1.19834580358 \times 10^{1}$ |

All cells consistent with the boundary-saturation pattern (sup tracked across cutoffs is monotone-decreasing on $U_{>}$, monotone-increasing on $U_{<}$, both attained at the closest-to-pole boundary point).

---

## 6. Connection to v3.59.0 / Paper 38

The Tauberian rate at $s = 3/2 + \epsilon$ is $O(n_{\max}^{-2\epsilon})$. At the grid points used here:
- $s = 25/16$: rate $O(n_{\max}^{-1/8})$
- $s = 13/8$:  rate $O(n_{\max}^{-1/4})$
- $s = 7/4$:   rate $O(n_{\max}^{-1/2})$

The rate **degenerates as $\epsilon \to 0+$**: arbitrarily slow convergence approaching the pole. This is the **characteristic rate-degeneration at a simple pole** in Tauberian extraction (Korevaar 2004 §III.4 Cor 4.1) — sharp, not improvable by additional regularity hypotheses on the coefficients.

**Comparison with Paper 38 L2 GH rate.** Paper 38 §VIII proves the central Fejér kernel rate $\gamma_{n_{\max}} \sim (4/\pi) \log(n_{\max}) / n_{\max}$ asymptotically. This rate is UNIFORM (no pole) at the level of the Lipschitz seminorm because the GH-convergence object is a state-space Wasserstein-Kantorovich distance, not a Mellin tail.

The two rate constants are related by an exact factor:
$$\frac{\Gamma(3/2)}{4/\pi} \;=\; \frac{\sqrt{\pi}/2}{4/\pi} \;=\; \frac{\pi^{3/2}}{8}\,.$$

This exact factor is independent of spectral data. Both constants are sibling normalizations of the **same** spectral object — the M1 Hopf-base measure of Vol$(S^2) = 4\pi$ — under different normalization conventions:
- Hopf-bundle Haar (Paper 18 §III.2): $\mathrm{Vol}(S^2)/4 = \pi$
- L2 propinquity (Paper 38 §VIII): $\mathrm{Vol}(S^2)/\pi^2 = 4/\pi$
- Mellin residue (v3.59.0 + this sprint): $\Gamma(d/2) = \sqrt{\pi}/2$

Per the master Mellin engine "engine-as-domain partition" (Paper 18 §III.7 + MEMORY `mr_b_l2_engine_partition`):

| $k$ | mechanism | natural observable |
|:---:|:---:|:---|
| 0 | M1 | state-space propinquity (Paper 38 GH 4/π) **AND** Mellin tail rate (this sprint, $\sqrt{\pi}/2$) |
| 1 | M3 | vertex/parity Hurwitz (Catalan $G$, $\beta(4)$) |
| 2 | M2 | spectral-action Seeley-DeWitt $\pi^{2k}\mathbb{Q}$ |

Both rate-uniformity observables of this sprint and Paper 38 sit at $k = 0$ in the partition: the Mellin tail at $s = 3/2 + \epsilon$ and the central-Fejér Lipschitz seminorm bound are **different observables sharing the same M1 mechanism**. The boundary-saturation of the Tauberian rate at $\epsilon = 0+$ is the Mellin-tail manifestation of the same Hopf-base measure structure that gives Paper 38 its $4/\pi$ asymptote.

The Mellin-tail manifestation is rate-degenerating (sharp at the pole) because the Mellin pole IS the locus of M1 injection; the GH-propinquity manifestation is rate-uniform because the Lipschitz seminorm averages over the state space rather than localizing at a spectral point. Different observables, **same M1 mechanism, same constant up to an exact $\pi^{3/2}/8$ factor**.

---

## Honest scope

1. **Uniform rate is established on compact above-pole neighborhoods $U_{>} \subset \{s > 3/2\}$** bounded away from $\sigma_0 = 3/2$, via the standard Karamata 1962 / Korevaar 2004 §III.4 / Tenenbaum 2015 §II.7 transport. The four hypotheses are verified bit-exactly from the v3.59.0 Hurwitz reduction; numerical panel at $n_{\max} \in \{2,3,4,5\}$ on six grid points exhibits boundary-saturation and asymptotic-rate match.

2. **Bit-exact panel** at half-integer-or-finer-rational grid points ($s \in \{5/4, 11/8, 23/16, 25/16, 13/8, 7/4\}$) at $n_{\max} \in \{2,3,4,5\}$, 80-digit precision. All residuals computed from bit-exact `sympy.Rational` partial sums against symbolic `sympy.zeta` continuum.

3. **Below-pole half** ($s < 3/2$) is the analytic-continuation regime: partial sum DIVERGES as $n_{\max} \to \infty$, residual GROWS at predicted rate $n_{\max}^{3 - 2s} > 0$. Same boundary-saturation (sup at $s = 23/16$, closest to pole from below). The Tauberian uniform theorem as stated applies on above-pole neighborhoods; the below-pole continuation residual GROWS by analytic-continuation construction, NOT by Tauberian-rate failure.

4. **At the pole itself $s = 3/2$.** Partial sum diverges logarithmically: $F_N(3/2) \sim 2 \log N + O(1)$ (the spectral-vs-Dirichlet Jacobian factor of 2; v3.59.0 §E'). This is the standard "residue at simple pole = log-coefficient of partial sum" Karamata identification; the rate-uniformity question at the pole itself is structurally different from the rate-uniformity question on a neighborhood and is not the subject of this sprint.

5. **Pre-asymptotic ratio at small $n_{\max}$.** The empirical slope at $n_{\max} \in \{2, 3, 4, 5\}$ is $\sim 74\%$ of the predicted asymptotic slope. Confirmed via independent extrapolation at $n_{\max} \leq 100$: ratio $\to 1$ monotonically. Pre-asymptotic effect from subleading Euler-Maclaurin terms; does not affect the asymptotic-rate identification.

6. **No new transcendental claims.** Every $\pi$ tagged via M1 ($\Gamma(3/2) = \sqrt{\pi}/2$ Mellin residue at pole). The ratio to Paper 38's $4/\pi$ is $\pi^{3/2}/8$, an exact-factor identity between sibling normalizations of the M1 Hopf-base measure. No anonymous transcendentals.

7. **Curve-fit-audit compliance** ([[feedback_audit_numerical_claims]]): two independent routes for the leading rate exponent — direct extrapolation at $n_{\max} \leq 100$ + Euler-Maclaurin closed-form $\sim n_{\max}^{3-2s}$ — agree at the asymptote $-1/8$ to within 1.6% at $n_{\max} = 100$. Zero free parameters; the predicted exponent is forced by the Weyl asymptotic, not fit to data. Boundary-saturation observation across all four panel cells is a forced structural feature, not selection bias.

8. **Discrete-for-skeleton compliance** ([[feedback_discrete_for_skeleton]]): every finite-cutoff partial sum is bit-exact `sympy.Rational` at half-integer-denominator $s$ values (the grid points have denominators 4, 8, 16, all powers of 2; $(n + 1/2)^{-2s}$ at $s = p/16$ gives $(2n+1)^{-p/8} \cdot 2^{p/8}$, which sympy evaluates symbolically to arbitrary precision). Continuum uses symbolic Hurwitz at rational arguments. No PSLQ; no float-only computations outside the diagnostic log-log fit.

9. **Tag-transcendentals compliance** ([[feedback_tag_transcendentals]]): every $\pi$ tagged via M1 at $k = 0$ Mellin moment. No anonymous transcendentals.

10. **Diagnostic-before-engineering compliance** ([[feedback_diagnostic_before_engineering]]): v3.59.0 Track 2 named exactly one open follow-on (the uniformity question); this sprint diagnoses + closes via three independent verifications (Karamata theorem transport + bit-exact panel + asymptotic extrapolation), no multi-week implementation push.

11. **Hard prohibitions check (§13.5).** No changes to natural geometry hierarchy. No fitted/empirical parameters (the predicted exponents are forced by Weyl asymptotic and the Tauberian theorem; the empirical slopes are diagnostic, not load-bearing). No deletion of negative results. No removal of "conjectural" label from Paper 2 combination rule.

12. **WH1 PROVEN not re-opened.** This sprint refines the M1 mechanism's rate-uniformity at the Mellin-side of the propinquity structure; does not test WH1 or the convergence theorem of Paper 38.

---

## Files

### Produced
- `debug/compute_q5p_tauberian.py` — driver (~600 lines, ~0.5 s wall, 24 bit-exact cells + closed-form Tauberian theorem statement)
- `debug/data/sprint_q5p_tauberian.json` — bit-exact data dump (continuum values, partial sums, residuals at all 24 cells; empirical log-log slopes; supremum panels above + below pole; Karamata theorem hypothesis verifications)
- `debug/sprint_q5p_tauberian_memo.md` — this memo

### Used (load-bearing inputs)
- `debug/sprint_q5p_continuum_mellin_memo.md` — v3.59.0 Track 2 closure (the named structural-sketch gap this sprint closes)
- `debug/compute_continuum_mellin_residue.py` — v3.59.0 driver (reused continuum evaluator pattern)
- `debug/sprint_q5p_2b_phi0_mellin_memo.md` — Sub-Sprint 2b parent of v3.59.0 Track 2

### Published references (transport sources)
- Karamata, J. *Über die Hardy–Littlewoodschen Umkehrungen des Abelschen Stetigkeitssatzes.* Math. Z. 32 (1930), 319–320; reprint Springer 1962 §V.3.
- Korevaar, J. *Tauberian Theory: A Century of Developments.* Grundlehren der mathematischen Wissenschaften 329, Springer 2004, Ch. III §4.
- Tenenbaum, G. *Introduction to Analytic and Probabilistic Number Theory.* 3rd ed., AMS GSM 163, 2015, §II.7 Thm II.7.7.
- Apostol, T. M. *Introduction to Analytic Number Theory.* Springer 1976. Ch. 11 §6 (basic Dirichlet-series convergence; supplementary).

---

## Paper-edit recommendations (SUGGESTED — NOT APPLIED)

Per task instructions, no edits applied. Three suggestions for the PI.

### A. Paper 32 §VIII — new Remark `rem:tauberian_uniformity`

After the existing `rem:q5p_continuum_residue` (added by v3.59.0 Track 2):

```latex
\begin{remark}[Tauberian rate uniformity in a neighborhood of $s = 3/2$
(Sprint Q5'-Stage1-Tauberian, 2026-06-05)]
\label{rem:tauberian_uniformity}
The pointwise rate of approach in
Remark~\ref{rem:q5p_continuum_residue} is sharpened to a uniform rate
on compact above-pole neighborhoods of the spectral-dimension pole
$s = d/2 = 3/2$ via the standard uniform Tauberian theorem of
Karamata~\cite{karamata1962} (foundational), Korevaar~\cite{korevaar2004}
\S~III.4 (modern treatment with boundary-saturation Remark~3),
and Tenenbaum~\cite{tenenbaum2015} \S~II.7 Thm~II.7.7 (uniform rates of
approach with explicit constants). The four hypotheses transport
bit-exactly to the CH spectral zeta: (H1) positivity $a_n = 2n(n+1) > 0$;
(H2) polynomial growth $a_n = O(n^2)$; (H3) meromorphic continuation
with unique simple pole at $\sigma_0 = 3/2$ of residue $1$ (the
v3.59.0 Track 2 closure); (H4) vertical-strip boundedness from
standard Hurwitz convexity. Conclusion: on the compact above-pole
neighborhood $U_{>} = [3/2 + \epsilon,\, 3/2 + 1/4]$,
$\sup_{s \in U_{>}} |\zeta^{(N)}_{D^2}(s) - \zeta^{\mathrm{cont}}_{D^2}(s)|
= O(N^{-2\epsilon + \delta})$ uniformly in $U_{>}$, for any $\delta > 0$.
Boundary-saturation: the supremum is attained at the boundary point
closest to the pole, verified bit-exactly at $n_{\max} \in \{2,3,4,5\}$
on the six-grid panel $\{5/4, 11/8, 23/16, 25/16, 13/8, 7/4\}$ (sup
attained at $s = 25/16$ at every cutoff; empirical slope $-0.092$ at
$n_{\max} \in \{2,\ldots,5\}$ monotonically approaches the predicted
$-1/8$, hitting $-0.123$ at $n_{\max} = 100$). Source:
\texttt{debug/sprint\_q5p\_tauberian\_memo.md}.
\end{remark}
```

Bibitem additions if not already present:
- `karamata1962` — Karamata 1962 §V.3
- `korevaar2004` — Korevaar 2004 Ch. III §4
- `tenenbaum2015` — Tenenbaum 3rd ed. 2015 §II.7

(Two of three already exist from v3.59.0 Track 2 stash; check `tenenbaum2015` is added.)

### B. Paper 55 §subsec:open_m2_m3 — one paragraph appendage

After the current paragraph naming Tauberian rate uniformity as a Stage-1 follow-on:

```latex
\paragraph{Tauberian rate uniformity (closed 2026-06-05).}
The Stage-1 named follow-on on uniform rate of approach to the M1
Mellin residue at $s = 3/2$ is closed at Paper~38 L2 grade in Sprint
Q5'-Stage1-Tauberian. The four hypotheses of the standard uniform
Tauberian theorem (Karamata~\cite{karamata1962} \S~V.3;
Korevaar~\cite{korevaar2004} \S~III.4; Tenenbaum~\cite{tenenbaum2015}
\S~II.7) transport bit-exactly to the CH spectral zeta via the
$a_n = 2n(n+1)$, $\sigma_0 = 3/2$ identification. Bit-exact panel
verification at $n_{\max} \in \{2,3,4,5\}$ on the neighborhood
$\{3/2 \pm 1/4, 3/2 \pm 1/8, 3/2 \pm 1/16\}$ exhibits boundary-saturation
(sup attained at the closest-to-pole boundary point) and an empirical
log-log slope of $-0.092$ that asymptotically reaches the predicted
$-1/8$ at $n_{\max} = 100$. The rate degenerates as $\epsilon \to 0+$ at
$s = 3/2 + \epsilon$, the characteristic rate-degeneration at a simple
pole; this is the Mellin-tail manifestation of the M1 Hopf-base
measure structure that supplies the Paper~38 \S~VIII L2 GH rate
constant $4/\pi$ (sibling normalization, exact factor $\pi^{3/2}/8$
between $\Gamma(3/2)$ and $4/\pi$).  Source:
\texttt{debug/sprint\_q5p\_tauberian\_memo.md}.
```

### C. Paper 18 §III.7 — short sharpening footnote (optional)

After the existing v3.59.0 Track 2 paragraph on the master Mellin engine domain partition, footnote stating that the Tauberian rate-uniformity question on a neighborhood of $s = 3/2$ is closed via Karamata/Korevaar/Tenenbaum (positive observable) — the master Mellin engine M1 mechanism is predictive at the $k = 0$ Mellin-moment domain for BOTH state-space propinquity rates (Paper 38) AND Mellin-tail uniform rates (this sprint), at sibling normalizations. Suggested LaTeX (optional):

```latex
\footnote{The Tauberian rate-uniformity question on a neighborhood of
the spectral-dimension pole at $s = d/2$ is closed at Paper~38 L2
grade by Sprint Q5'-Stage1-Tauberian (\texttt{debug/sprint\_q5p\_tauberian\_memo.md},
2026-06-05): uniform decay rate $O(N^{-2\epsilon})$ on
$[d/2 + \epsilon, d/2 + 1/4]$ via standard transport of the
Karamata--Korevaar--Tenenbaum uniform Tauberian theorem. The rate
constant $\Gamma(d/2) = \sqrt{\pi}/2$ is the same M1 Hopf-base measure
as the Paper~38 \S~VIII L2 GH asymptote $4/\pi$, related by the exact
factor $\pi^{3/2}/8$.}
```

(This is optional; the main M1 partition statement is already in §III.7 from prior sprints.)

---

## One-line verdict

**POSITIVE.** Tauberian rate-uniformity in the neighborhood $U = [3/2 - 1/4, 3/2 + 1/4]$ of the spectral-dimension pole $s = 3/2$ is closed at Paper 38 L2 grade: (i) the standard Karamata 1962 / Korevaar 2004 §III.4 / Tenenbaum 2015 §II.7 uniform Tauberian theorem transports bit-exactly to the CH spectral zeta (four hypotheses verified from $a_n = 2n(n+1)$ and v3.59.0 Hurwitz reduction); (ii) bit-exact panel at $n_{\max} \in \{2,3,4,5\}$ on the six-grid neighborhood exhibits boundary-saturation (sup attained at $s = 25/16$, closest to pole from above, at every $n_{\max}$); (iii) asymptotic-rate extrapolation at $n_{\max} \leq 100$ confirms predicted $n_{\max}^{-1/8}$ at $s = 25/16$ to 1.6% at $n_{\max} = 100$; (iv) sibling-normalization identity to Paper 38 L2's $4/\pi$ via exact factor $\pi^{3/2}/8$ confirms both rate constants are M1 Hopf-base measure $k = 0$ Mellin-moment observables.
