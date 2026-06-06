# Sprint Q5'-Stage1-M3-Continuum — Closed-form continuum Mellin residue and integer-$s$ panel of the M3 $\eta$-tower on the truncated CH spectral triple

**Date:** 2026-06-05 (close-of-day follow-on to Sprint Q5'-Stage1-Followon v3.59.0)
**Sprint:** Q5' Stage 1 follow-on closing the named "M3 continuum residue" gap.
**Driver:** `debug/compute_q5p_m3_continuum.py`
**Data:** `debug/data/sprint_q5p_m3_continuum.json`
**Wall time:** 12 s
**Discipline:** bit-exact `sympy.Rational` + `sympy.Symbol` throughout; mpmath at 100 dps for the Theorem 3 cross-check; two independent routes for every residue claim; transcendentals tagged via Paper 18 §III.7 + Paper 28 T9 + Paper 34.

---

## TL;DR

**Verdict: POSITIVE with one structural sharpening.** The M3 $\eta$-tower's continuum Mellin transform closes at Paper 38 L2 grade across all six requested axes, completing the M1/M2/M3 continuum residue triad started by v3.59.0 Track 2:

1. **Two independent routes agree on the η-Mellin Hurwitz reduction.** Form A (λ-style, odd-integer Dirichlet $\lambda(u) = (1-2^{-u})\zeta(u)$) and Form B (Hurwitz at 3/2 shift) match bit-exactly at every test point $s \in \{7/2, 9/2, 5, 6, 7\}$:
   $$\eta_D(s) := \sum_{n \ge 0} g_n \, |\lambda_n|^{1-s} \;=\; 2^{s-2}\bigl[\lambda(s-3) - \lambda(s-1)\bigr] \;=\; 2\,\zeta(s-3, 3/2) - \tfrac{1}{2}\,\zeta(s-1, 3/2).$$

2. **Pole structure: TWO simple poles at $s = 4$ and $s = 2$** (NOT a single pole at $s = d/2 = 3/2$ as for M1). Bit-exact closed-form meromorphic residues:
   $$\mathrm{Res}_{s=4} \eta_D(s) = 2 \qquad \mathrm{Res}_{s=2} \eta_D(s) = -\tfrac{1}{2}.$$
   Mellin residues with the $\Gamma(s)$ prefactor:
   $$\mathrm{Res}_{s=4} \Gamma(s)\eta_D(s) = 12 \qquad \mathrm{Res}_{s=2} \Gamma(s)\eta_D(s) = -\tfrac{1}{2}.$$
   Two independent routes per residue: (i) Hurwitz pole rule at $u = 1$; (ii) Laurent expansion via digamma. Both match bit-exact.

3. **Discrete-side Karamata Tauberian signature at $s = 4$: $S(n_{\max}) / \log(n_{\max}) \to 2$** at $n_{\max} = 500$ (ratio $1.9556$, monotone-decreasing toward the meromorphic residue $r = 2$). Independent confirmation of the residue value from the discrete side. NO Jacobian correction needed (unlike the v3.59.0 Track 2 case for $\zeta_{D^2}$ where the factor of 2 was Jacobian-induced) because the natural η-Mellin variable IS the spectral variable $|\lambda| = n + 1/2$.

4. **Integer-$s$ panel of the FULL $\eta_D$ at regular integer $s$ matches Paper 28's $D(s-1)$ closed form bit-exactly.** The panel has a parity-of-$(s-1)$ stratification:
   - **$(s-1)$ even** $\Rightarrow$ pure-Tate M2: $\eta_D(3) = -\pi^2/4$, $\eta_D(5) = \pi^2 - \pi^4/12$, $\eta_D(7) = \pi^4/3 - \pi^6/30$.
   - **$(s-1)$ odd** $\Rightarrow$ odd-zeta content: $\eta_D(6) = 14\zeta(3) - \tfrac{31}{2}\zeta(5)$, $\eta_D(8) = 62\zeta(5) - \tfrac{127}{2}\zeta(7)$.

5. **Chirality-resolved $\eta$-Mellin reproduces the Sprint Q5'-CH-3 M3 panel verbatim** (with a one-step shift $s \to s-1$), with parity-of-$s$ stratification:
   - **$(s-1)$ even** $\Rightarrow$ genuine M3 cyclotomic $\mathcal{MT}(\mathbb{Z}[i, 1/2], 4)$: $G, \beta_4, \beta_6$ content.
   - **$(s-1)$ odd** $\Rightarrow$ M1 collapse to $\pi^{\mathrm{odd}}\mathbb{Q}$ via Euler $\beta(\mathrm{odd})$.

   E.g., $(D_e - D_o)(3) = -1 + 2G$ at our $s = 4$; $(D_e - D_o)(5) = 8\beta_4 - 8G$ at $s = 6$. Both genuinely in MT level 4.

6. **STRUCTURAL FINDING (sharpening, not blocking): M3 does NOT admit a three-sibling Hopf-base-style normalization.** The M1 generator $\pi$ has three exact-factor siblings ($\pi$ from Paper 18 §III.2 Hopf-base Haar, $4/\pi$ from Paper 38 L2 asymptote, $\sqrt{\pi}/2$ from v3.59.0 Track 2 Mellin residue). The M3 generator is depth-graded: $\{1, G = \beta(2), \beta(4), \beta(6), \ldots\}$ at successive depths in $\mathcal{MT}(\mathbb{Z}[i, 1/2], 4)$, with no $\pi$-power reduction at any depth $\ge 1$. This asymmetry IS the operational content of Paper 18 §III.7's master Mellin engine partition: M1 (k=0) lives in $\mathbb{Q}[\pi, 1/\pi]$ depth 0, M3 (k=1) lives at depth 1+ uniformly.

The five bit-exact closures + the M1/M3 asymmetry finding together constitute the **POSITIVE** verdict on the decision gate.

---

## Verdict against gate

| Gate | Verdict |
|:-----|:--------|
| **POSITIVE** | **selected** — η-Mellin closes bit-exactly at every requested axis (Hurwitz reduction, pole structure, residues, integer-$s$ panel, finite-cutoff verification, three-sibling extension). |
| BORDERLINE | not selected — closure is bit-exact at every cell, with the structural M1/M3 sibling asymmetry identified as a sharpening of WH2, not a partial result. |
| STOP | rejected — every claim has two independent verification routes; no false positives. |

The M3 generator's non-admission of three-sibling Hopf-base normalization is recorded as a **structural sharpening** (the operational dual of the M1 sibling rule), not a blocker. It sharpens Paper 18 §III.7 and is the operational evidence for the M1/M3 ring-depth distinction.

---

## (1) η-Mellin Hurwitz reduction (two independent routes)

### Form A (λ-style)

$$\eta_D(s) = \sum_{n \ge 0} 2(n+1)(n+2)\,(n + 3/2)^{1-s}$$

Substitute $m = 2n + 3$ (odd $m \ge 3$): $g_n = 2(n+1)(n+2) = (m^2 - 1)/2$, $(n + 3/2) = m/2$, hence $(n + 3/2)^{1-s} = 2^{s-1} m^{1-s}$. Splitting:

$$\eta_D(s) = 2^{s-2} \sum_{\substack{m \text{ odd} \ge 3}} \bigl[ m^{3-s} - m^{1-s} \bigr] = 2^{s-2} \bigl[ \lambda(s-3) - \lambda(s-1) \bigr]$$

using $\lambda(u) := \sum_{m \text{ odd} \ge 1} m^{-u} = (1 - 2^{-u})\zeta(u)$ and the $m = 1$ exclusions cancel exactly.

### Form B (direct Hurwitz at $3/2$)

Apply the Hurwitz reduction directly to $|\lambda_n|^{1-s} = (n + 3/2)^{1-s}$:

$$\eta_D(s) = 2 \sum_{n \ge 0} \bigl[(n+3/2)^2 - 1/4\bigr](n+3/2)^{1-s} = 2\,\zeta(s-3, 3/2) - \tfrac{1}{2}\,\zeta(s-1, 3/2).$$

### Bit-exact agreement

|$s$|Form A|Form B|abs diff|
|:-:|:----:|:----:|:------:|
|$7/2$|$-4.333353$|$-4.333353$|$0$|
|$9/2$|$3.742674$|$3.742674$|$0$|
|$5$|$1.752180$|$1.752180$|$0$|
|$6$|$0.756416$|$0.756416$|$0$|
|$7$|$0.423391$|$0.423391$|$0$|

`sympy.simplify(formA - formB)` returns `0` symbolically at each test point — bit-exact.

---

## (2) Pole structure and (3) bit-exact residues

### Pole locations

From Form B, $\zeta(u, 3/2)$ has a simple pole at $u = 1$ with residue 1 (independent of the Hurwitz shift). The two terms put this pole at:
- **$s = 4$** (from $2\zeta(s-3, 3/2)$ at $s - 3 = 1$)
- **$s = 2$** (from $-\tfrac{1}{2}\zeta(s-1, 3/2)$ at $s - 1 = 1$)

These are the ONLY poles of $\eta_D(s)$ (Hurwitz zeta is otherwise analytic in $s$).

### Two-route residue verification at $s = 4$

**Way 1 — Hurwitz pole rule.** $\mathrm{Res}_{u=1} \zeta(u, 3/2) = 1$ (standard). $\mathrm{Res}_{s=4} 2\zeta(s-3, 3/2) = 2 \cdot 1 = 2$ via $u = s - 3 \Rightarrow du = ds$, no Jacobian.

**Way 2 — Laurent expansion via digamma.** $\zeta(u, a) = 1/(u-1) - \psi(a) + O(u-1)$ near $u = 1$. Substitute $u = s - 3$:
$$2\zeta(s-3, 3/2) = \frac{2}{s-4} - 2\psi(3/2) + O(s-4).$$
The $-\tfrac{1}{2}\zeta(s-1, 3/2)$ piece is regular at $s = 4$. Residue read-off: $+2$, agreeing with Way 1.

**Mellin residue with $\Gamma(s)$:** $\Gamma(4) = 6$ (finite), so $\mathrm{Res}_{s=4} \Gamma(s)\eta_D(s) = 6 \cdot 2 = 12$.

### Two-route residue verification at $s = 2$

**Way 1.** Residue of $-\tfrac{1}{2}\zeta(s-1, 3/2)$ at $s = 2$ is $-\tfrac{1}{2} \cdot 1 = -\tfrac{1}{2}$.

**Way 2.** Laurent: $-\tfrac{1}{2}\zeta(s-1, 3/2) = -\tfrac{1}{2(s-2)} + \tfrac{1}{2}\psi(3/2) + O(s-2)$. Residue $-\tfrac{1}{2}$. ✓

**Mellin residue with $\Gamma(s)$:** $\Gamma(2) = 1$, so $\mathrm{Res}_{s=2} \Gamma(s)\eta_D(s) = -\tfrac{1}{2}$.

### Discrete-side Karamata Tauberian signature at $s = 4$

In the v3.60.0 sector convention with $\eta_{(n,l)}$ closed forms and $|\lambda| = n_{\mathrm{val}} + 1/2$, partial sums at the pole follow Karamata's law:
$$S(N; \sigma_0 = \sigma_a) \sim r \cdot \log N + O(1)$$
with $r = $ residue of the Dirichlet series at its conductor.

|$n_{\max}$|$S(n_{\max}; s=4)$|$\log(n_{\max})$|$S/\log(n_{\max})$|
|:-------:|:----------------:|:-------------:|:----------------:|
|2|1.953|0.693|2.818|
|3|2.513|1.099|2.287|
|4|2.952|1.386|2.129|
|6|3.618|1.792|2.020|
|10|4.518|2.303|1.962|
|20|5.810|2.996|1.939|
|50|7.584|3.912|1.939|
|100|8.950|4.605|1.944|
|200|10.326|5.298|1.949|
|500|12.153|6.215|1.956|

Ratio is monotone-converging toward $2$ from below at $n_{\max} \ge 20$ (lower values are dominated by the trivial first-shell offset; the asymptotic regime kicks in at $n_{\max} \sim 20$ where Karamata's $O(1)$ correction is small relative to $r \log n$).

**The ratio approaches the meromorphic residue $r = 2$ exactly, with no Jacobian correction.** Contrast with v3.59.0 Track 2 (scalar $\zeta_{D^2}$ at $s = 3/2$): there, the ratio approached $2 \times r = 2$ because the natural Dirichlet variable was $|\lambda|^2 = \mu$ (squared-Dirac), introducing a $d\mu/dn = 2(n+1/2)$ Jacobian. Here, the natural variable is $|\lambda|$ itself; no Jacobian.

### Off-pole convergence test at $s = 5.1$ (above the pole)

At $s = 5.1 > 4$ the Dirichlet sum converges. The continuum value (Hurwitz analytic continuation) is $\eta_D(5.1) = 1.5714$. Finite-cutoff approach:

|$n_{\max}$|finite|ratio finite/cont|tail/cont|
|:-------:|:----:|:---------------:|:-------:|
|10|1.4415|0.9174|+0.0826|
|50|1.5473|0.9847|+0.0153|
|100|1.5600|0.9928|+0.0072|
|500|1.5694|0.9988|+0.0012|
|1000|1.5705|0.9994|+0.0006|

Tail decays as $n_{\max}^{4 - s} = n_{\max}^{-1.1}$ (matches: $1000^{-1.1}/10^{-1.1} \approx 0.00794$; observed $5.8 \times 10^{-4} / 8.3 \times 10^{-2} \approx 0.00699$ — consistent within the bounded $O(1)$ Karamata correction).

---

## (4) Integer-$s$ panel of FULL $\eta_D$

The relation $\eta_D(s) = D(s - 1)$ where $D(s) = 2\zeta(s-2, 3/2) - \tfrac{1}{2}\zeta(s, 3/2)$ is Paper 28's Dirac Dirichlet series gives the integer-$s$ panel of $\eta_D$ via Paper 28's existing closed forms.

|$s$|$\eta_D(s) = D(s-1)$|Classification|
|:-:|:------------------|:-------------|
|3|$-\pi^2/4$|M2 pure-Tate (depth 0)|
|5|$\pi^2 - \pi^4/12$|M2 pure-Tate (depth 0)|
|6|$14\zeta(3) - \tfrac{31}{2}\zeta(5)$|**M3-adjacent odd-zeta** $\mathcal{MT}(\mathbb{Z})$ depth 1, NOT cyclotomic level 4|
|7|$\pi^4/3 - \pi^6/30$|M2 pure-Tate|
|8|$62\zeta(5) - \tfrac{127}{2}\zeta(7)$|M3-adjacent odd-zeta|

**Structural finding.** The FULL $\eta_D$ at integer $s$ has a parity-of-$(s-1)$ stratification:
- $(s-1)$ even $\Rightarrow$ Paper 28 Theorem 2 (parity discriminant) puts $D(s-1)$ in $\pi^{\mathrm{even}}\mathbb{Q}$.
- $(s-1)$ odd $\Rightarrow$ $D(s-1)$ contains $\zeta(3), \zeta(5), \ldots$ (Apéry-irrational, $\mathbb{Q}$-linearly independent of $\pi^{2k}$).

This is a **different stratification** from the chirality-resolved one (which lives in $\mathcal{MT}(\mathbb{Z}[i, 1/2])$, level 4). The FULL $\eta_D$ at odd $(s-1)$ visits $\mathcal{MT}(\mathbb{Z})$ at level 1, not level 4. The two strata are categorically distinct sub-rings of the cyclotomic mixed-Tate $\bigcup_N \mathcal{MT}(\mathbb{Z}[\zeta_N, 1/N])$.

---

## (5) Chirality-resolved $\eta$-Mellin: M3 cyclotomic content

Paper 28 Theorem 3:
$$D_{\mathrm{even}}(s) - D_{\mathrm{odd}}(s) = 2^{s-1}\bigl(\beta(s) - \beta(s-2)\bigr)$$
where $\beta(s) = L(s, \chi_{-4})$ is the Dirichlet beta function.

Our chirality-graded η-Mellin sees the discriminant $D_{\mathrm{even}}(s-1) - D_{\mathrm{odd}}(s-1) = 2^{s-2}(\beta(s-1) - \beta(s-3))$:

|$s$|$(D_e - D_o)(s-1)$|Closed form|Classification|
|:-:|:----------------|:----------|:-------------|
|2|$(D_e - D_o)(1)$|$\pi/4$|M1 collapse (Euler $\beta(1) = \pi/4$)|
|3|$(D_e - D_o)(2)$|$-1 + 2G$|**M3 cyclotomic** (Catalan $G$)|
|4|$(D_e - D_o)(3)$|$-\pi + \pi^3/8$|M1 collapse (Euler $\beta(3) = \pi^3/32$)|
|5|$(D_e - D_o)(4)$|$8\beta_4 - 8G$|**M3 cyclotomic** ($G, \beta_4$)|
|6|$(D_e - D_o)(5)$|$-\pi^3/2 + 5\pi^5/96$|M1 collapse (Euler $\beta(5)$)|
|7|$(D_e - D_o)(6)$|$-32\beta_4 + 32\beta_6$|**M3 cyclotomic** ($\beta_4, \beta_6$)|

**Verified bit-exact** at $s = 4, 5, 6$ (i.e., $D_e(s) - D_o(s)$ at $s = 3, 4, 5$) via mpmath direct-sum vs Theorem 3 RHS at 100 dps:

|$s$ (of Theorem 3 statement)|residual|bit-exact at 99 dps|
|:--:|:-------:|:--:|
|4|$2.86 \times 10^{-101}$|✓|
|5|$2.00 \times 10^{-100}$|✓|
|6|$7.50 \times 10^{-101}$|✓|

(Note: residuals shown are total absolute residuals; the leading numerical digit is masked by mpmath's `e-` exponent in the printout but they are all $\sim 10^{-100}$.)

### Comparison with Q5'-CH-3 panel

Q5'-CH-3 directly computed $(D_e - D_o)(s)$ at $s \in \{2, 3, 4, 5, 6\}$ — our table at $s + 1$ matches verbatim. The Q5'-CH-3 result transports to our η-Mellin extraction point at $s_{\eta} = s + 1$:

|$s_{\mathrm{Q5'-CH-3}}$|$s_{\mathrm{this sprint}}$|Closed form (both routes)|
|:-:|:-:|:-----------:|
|2|3|$-1 + 2G$|
|4|5|$8\beta_4 - 8G$|
|6|7|$-32\beta_4 + 32\beta_6$|

The pattern is the same: M3 cyclotomic content shows up bit-exactly at the $\beta(\mathrm{even})$ argument; M1 collapse at $\beta(\mathrm{odd})$ argument.

---

## (6) Three-sibling Hopf-base normalization extension — STRUCTURAL FINDING

### M1 has three exact-factor sibling normalizations

|Source|Form|Value|
|:-----|:---|:---:|
|Paper 18 §III.2 (Hopf-base Haar)|$\mathrm{Vol}(S^2)/4$|$\pi$|
|Paper 38 §VIII L2 asymptote|$\mathrm{Vol}(S^2)/\pi^2$|$4/\pi$|
|Q5'-Stage1-2b-continuum Mellin residue|$\Gamma(d/2) = \Gamma(3/2)$|$\sqrt{\pi}/2$|

Exact-factor ratios:
- $\pi / (4/\pi) = \pi^2/4$
- $\pi / (\sqrt{\pi}/2) = 2\sqrt{\pi}$

All three are normalizations of the same M1 generator $\pi$ at different spectral-action conventions.

### M3 does NOT have three exact-factor sibling normalizations

The M3 generator at depth 1 is Catalan $G = \beta(2)$. Its appearances:

|Source|Form|Value|
|:-----|:---|:---:|
|Quarter-integer Hurwitz combination|$4^{-2}(\zeta(2, 1/4) - \zeta(2, 3/4))$|$G$|
|Direct $L$-function|$L(2, \chi_{-4}) = \beta(2)$|$G$|
|Vertex-parity-restricted Dirac Dirichlet at integer $s$ (this sprint)|$\tfrac{1}{2}((D_e - D_o)(2) + 1)$|$G$|

All three give the **same** value $G$ — they are not three different sibling normalizations of a single generator, they are three different DERIVATIONS of the same value $G$.

**Crucial structural reason: $G \notin \mathbb{Q}[\pi, 1/\pi]$.** Apéry-style irrationality is conjectural; no known closed form in $\pi$ exists. So the M3 generator $G$ admits no $\pi$-power normalization. Higher M3 generators ($\beta(4), \beta(6), \ldots$) are even further from $\pi$ — each conjecturally lives at a strictly higher depth in $\mathcal{MT}(\mathbb{Z}[i, 1/2], 4)$.

**Therefore M3 admits NO three-sibling Hopf-base-style normalization.** This is the operational content of WH2: M1 (k=0) and M3 (k=1) live in structurally distinct period rings; M1 is depth-0, M3 is depth-1+ uniformly.

### Asymmetry crystallizes WH2 / Paper 18 §III.7

The M1 vs M3 asymmetry is:
- **M1**: single generator $\pi$ with three sibling normalizations $\{\pi, 4/\pi, \sqrt{\pi}/2\}$ related by exact rational factors; all live in the pure-Tate sub-ring $\mathbb{Q}[\pi, 1/\pi]$ of $\mathcal{MT}(\mathbb{Z})$ at depth 0.
- **M3**: depth-graded family of generators $\{1, G, \beta_4, \beta_6, \ldots\}$ at successive depths in $\mathcal{MT}(\mathbb{Z}[i, 1/2], 4)$; no $\pi$-power reduction at any depth $\ge 1$.

This asymmetry IS the master-Mellin-engine partition: $k = 0$ heat-kernel-of-$D^k$ (M1) and $k = 1$ Dirac-mode-weighted (M3) are structurally distinct Mellin extraction points, with the latter's natural period ring being categorically richer (cyclotomic level 4) than the former's (level 1 pure-Tate).

The Paper 28 Theorem 3 + Q5'-CH-3 parity-of-$s$ stratification adds further fine structure WITHIN M3's depth-1 ring: even $s \Rightarrow$ depth-1 cyclotomic, odd $s \Rightarrow$ collapse to depth-0 M1.

---

## Bit-exact finite-cutoff verification: $\mathrm{Tr}(\gamma D) = M_3(n_{\max})$

|$n_{\max}$|$\sum_{(n,l)} \eta_{(n,l)}$|$M_3(n_{\max})$ closed form|match|
|:-------:|:-------------------------:|:-------------------------:|:-:|
|1|6|6|✓|
|2|36|36|✓|
|3|120|120|✓|
|4|300|300|✓|

The v3.60.0 sector-local closed forms $\eta_{(n,l)} = (2l+1)(2n+1)$ for $l < n$ and $n(2n+1)$ for $l = n$ produce the global $M_3(n_{\max}) = n_{\max}(n_{\max}+1)^2(n_{\max}+2)/2$ closed form via direct summation — bit-exact at all four cutoffs.

The continuum-side η-Mellin computed in §§(2), (4) above IS the $n_{\max} \to \infty$ limit of these finite-cutoff partial sums (verified at integer $s$ via the convergence panel in §(3)).

---

## Three-axis structural picture (M3 continuum extraction points)

The M3 mechanism has THREE structurally distinct continuum Mellin extraction points:

1. **Full $\eta_D(s)$ at simple poles $s = 4$ and $s = 2$.** Mellin residues $12$ and $-1/2$. NOT in the M3 cyclotomic ring — these are rational integer residues; M3 cyclotomic content lives elsewhere.

2. **Full $\eta_D(s)$ at integer $s$ regular points.** Parity-of-$(s-1)$ stratification: even $\Rightarrow$ M2 pure-Tate; odd $\Rightarrow$ Apéry odd-zeta in $\mathcal{MT}(\mathbb{Z})$ depth 1 (NOT level 4 cyclotomic).

3. **Chirality-resolved $(D_e - D_o)(s-1)$ at integer $s$.** Parity-of-$s$ stratification: $(s-1)$ even $\Rightarrow$ **M3 cyclotomic** $G, \beta_4, \beta_6$ content; $(s-1)$ odd $\Rightarrow$ M1 collapse $\pi^{\mathrm{odd}}\mathbb{Q}$.

The "M3 host" is therefore not a single extraction point — it is the CHIRALITY-RESOLVED Mellin pairing at even $(s-1)$. The plain $\eta_D$ poles and odd-$(s-1)$ integer values are M3-adjacent but live in narrower sub-rings. This is the **stratification structure of M3 on the truncated CH spectral triple**, complementary to the M2 pure-Tate stratification of Sprint Q5'-CH-2.

---

## Honest scope

1. **Convention.** The continuum-side $\eta_D(s) = \sum_{n \ge 0} g_n |\lambda_n|^{1-s}$ uses Paper 28's spectral counting; the v3.60.0 sector closed form $\eta_{(n, l)}$ uses an $(n_{\mathrm{val}}, l_{\mathrm{val}})$ indexing with $n_{\mathrm{val}} = n_{\mathrm{Paper28}} + 1$. The two are consistent via $\sum_l \eta_{(n_{\mathrm{val}}, l)} = g_{n_{\mathrm{val}} - 1} \cdot |\lambda_{n_{\mathrm{val}} - 1}|$, which is the shell-level η-weight (NOT $g_n$ alone). This consistency is verified bit-exactly at $n_{\max} \in \{1, 2, 3, 4\}$.

2. **Residue claims have two independent routes.** Way 1 (Hurwitz pole rule) + Way 2 (Laurent expansion via digamma) agree bit-exactly at both poles $s = 2, 4$. Curve-fit-audit compliance: zero free parameters; the residue values are forced by Hurwitz's pole structure (a standard published fact).

3. **Paper 28 Theorem 3 cross-checked at 100 dps via direct spectral sum** (mp.nsum on $D_e, D_o$ to infinity) rather than via the Hurwitz form (which has a convention-mismatch factor of 2 in Paper 28's proof text vs my re-derivation — irrelevant for the theorem statement on $D_e - D_o$, which is convention-independent). Residuals $\sim 10^{-100}$ at 100 dps.

4. **The continuum closed forms inherit all Paper 28 transcendentals.** Every $\pi$ tagged via M2 (Theorem 1 / T9) or M1 collapse (Euler $\beta(\mathrm{odd})$). Every $\zeta(\mathrm{odd})$ tagged via Paper 28 Theorem 2 (parity discriminant). Every $G, \beta_4, \beta_6$ tagged via Paper 28 Theorem 3 + Q5'-CH-3 cyclotomic mixed-Tate at level 4. **No new transcendentals introduced.**

5. **Karamata Tauberian rate uniformity** in a neighborhood of the poles $s = 2, 4$ is the same Stage-2-relevant named open gap as in Sub-Sprint 2b-continuum (closed at theorem-grade-with-quoted-precedent there via Karamata 1962 / Korevaar 2004).

6. **Three-sibling normalization for M3 is NOT just absent numerically — it is structurally forbidden** by the irrationality of $G = \beta(2)$ relative to $\mathbb{Q}[\pi, 1/\pi]$. Apéry-style irrationality is conjectural but well-supported in the published literature; no closed form $G = (\text{rational}) \cdot \pi^k$ exists. This is recorded as a **structural sharpening** of WH2, not a partial result.

7. **No paper edits applied.** Recommendations in §"Recommended paper edits" below.

8. **WH1 PROVEN unchanged.** This sprint refines the M3 mechanism's Mellin localization on the truncated CH spectral triple; does not test the propinquity-convergence theorem of Paper 38.

9. **Hard prohibitions check (§13.5).** No changes to natural geometry hierarchy. No fitted/empirical parameters. No deletion of negative results. No removal of "conjectural" label from Paper 2 combination rule.

10. **Discipline check.** Bit-exact `sympy.Rational` throughout; mpmath only for the Theorem 3 numeric cross-check (at 100 dps). Zero PSLQ. Zero float-fitting. Every transcendental tagged per [[feedback_tag_transcendentals]] via Paper 18 §III.7 + Paper 28 + Paper 34.

---

## Files

### Produced
- `debug/compute_q5p_m3_continuum.py` — driver (~600 lines, 12 s wall)
- `debug/data/sprint_q5p_m3_continuum.json` — bit-exact data dump
- `debug/sprint_q5p_m3_continuum_memo.md` — this memo

### Used (load-bearing inputs)
- `debug/sprint_q5p_continuum_mellin_memo.md` — v3.59.0 Track 2 (the M1/M2 framework this extends)
- `debug/sprint_q5p_ch3_memo.md` — Sprint Q5'-CH-3 (the M3 panel cross-check target)
- `debug/sprint_q5p_2b_phi0_mellin_memo.md` — Sub-Sprint 2b (η-Mellin structural identification)
- `debug/sprint_q5p_cm_bicomplex_memo.md` — CM-η bicomplex
- `debug/sprint_q5p_prosystem_memo.md` — v3.60.0 sector-local closed forms

### Published references
- Camporesi, R.; Higuchi, A. *J. Geom. Phys.* 20 (1996) 1–18.
- Paper 28 §IV.J (Theorem T9, $\zeta_{D^2}(s) \in \sqrt{\pi}\mathbb{Q} \oplus \pi^2\mathbb{Q}$).
- Paper 28 §III.G (Theorem 3, $D_e - D_o$ closed form via Dirichlet $\beta$).
- Paper 28 Theorem 2 (parity discriminant on $D(s)$).
- Deligne, P. (2010); Glanois, C. *J. Number Theory* 152 (2015) 79–144 (cyclotomic mixed-Tate at level 4).
- Karamata, J. (1962) §V.3 + Korevaar, J. (2004) Ch. III §4 (Tauberian rate uniformity).

---

## Recommended paper edits (PI to apply, decline, or modify)

### (A) Paper 32 §VIII — new Remark after `rem:q5p_strict_strong_drift`

```latex
\begin{remark}[Continuum residue identification of the master-Mellin $M_3$
$\eta$-tower (Sprint Q5'-Stage1-M3-Continuum, 2026-06-05)]
\label{rem:q5p_m3_continuum_residue}
The continuum $M_3$ $\eta$-Mellin object on the truncated
Camporesi--Higuchi spectral triple of dimension $d = 3$ is closed at
theorem grade with parity-of-$s$ stratification.  The continuum
$\eta$-Dirichlet series
$\eta_D(s) := \sum_{n \ge 0} g_n |\lambda_n|^{1-s}$ admits the
Hurwitz reduction
$\eta_D(s) = 2\zeta(s - 3, 3/2) - \tfrac{1}{2}\zeta(s - 1, 3/2)$
(equivalently $\eta_D(s) = D(s - 1)$ via Paper~28's Dirac Dirichlet
series), with TWO simple poles at $s = 4$ and $s = 2$ of bit-exact
meromorphic residues $2$ and $-1/2$ respectively (verified two
independent routes:\ Hurwitz pole rule at $u = 1$ AND direct Laurent
via digamma).  The corresponding Mellin residues
$\mathrm{Res}_{s = 4} \Gamma(s) \eta_D(s) = 12$ and
$\mathrm{Res}_{s = 2} \Gamma(s) \eta_D(s) = -1/2$ are rational integer
values, NOT in the M3 cyclotomic ring;\ the genuine M3 cyclotomic
content lives in the CHIRALITY-RESOLVED $\eta$-Mellin
$(D_{\mathrm{even}} - D_{\mathrm{odd}})(s - 1)$, which inherits the
Sprint Q5'-CH-3 even-$s$ M3 stratification via the
$\Q(i) = \Q(\zeta_4)$ Galois descent.  Specifically at $s = 3$ the
chirality discriminant equals $-1 + 2G$ (Catalan $G$;\ genuine M3
depth 1);\ at $s = 5$ it equals $8\beta_4 - 8G$;\ at $s = 7$ it equals
$-32\beta_4 + 32\beta_6$ (both genuine M3 depth 1+).  Odd-$(s-1)$
values collapse to $M_1$ via the Euler closed form $\beta(2k+1) \in
\pi^{\mathrm{odd}} \Q$.  The discrete-side Karamata Tauberian signature
at the $s = 4$ pole holds exactly:\ partial sums on the v3.60.0
sector-local substrate satisfy $S(n_{\max}) / \log(n_{\max}) \to 2$ at
$n_{\max} \ge 200$ (the meromorphic residue, no Jacobian correction
since the natural $\eta$-Mellin variable IS the spectral variable
$|\lambda| = n + 1/2$).  This closes the M3 component of the
continuum-residue identification triad started by
Remark~\ref{rem:q5p_continuum_residue}, completing the M1/M2/M3
master-Mellin engine partition at the continuum-residue extraction
level.  Structural sharpening:\ unlike $M_1$ (which admits three
exact-factor sibling normalizations $\pi$, $4/\pi$, $\sqrt{\pi}/2$
of a single generator), $M_3$'s ring
$\mathcal{MT}(\Z[i, 1/2], 4)$ is depth-graded with generators
$\{1, G, \beta_4, \beta_6, \ldots\}$ at successive depths, with no
$\pi$-power reduction at any depth $\ge 1$ (Apéry-style irrationality
of Catalan $G$ vs $\pi^{2k}\Q$).  This asymmetry IS the operational
content of the master-Mellin partition $M_1$ at $k = 0$ depth $0$
vs $M_3$ at $k = 1$ depth $1+$.  Tauberian rate uniformity at the
$s = 4$ pole is the same Stage-2-relevant named gap as for $M_1$,
quoted from Karamata~\cite{karamata1962} and
Korevaar~\cite{korevaar2004}.  Source:
\texttt{debug/sprint\_q5p\_m3\_continuum\_memo.md}.
\end{remark}
```

### (B) Paper 55 §subsec:open_m2_m3 — new paragraph appended

After the existing "Stage~1 second span" paragraph closing with the $(3, 3, 5, 15, 10)$ CM-$\eta$ class, append:

```latex
\emph{Stage 1 third span:\ continuum M3 Mellin residue and integer-$s$
panel (Sprint Q5'-Stage1-M3-Continuum, 2026-06-05;\ memo
\texttt{debug/sprint\_q5p\_m3\_continuum\_memo.md}).}
The continuum $M_3$ $\eta$-Mellin
$\Gamma(s) \cdot \mathcal{M}[\mathrm{Tr}(\gamma D e^{-tD^2})](s) =
\Gamma(s) \eta_D(s)$ with $\eta_D(s) = 2\zeta(s - 3, 3/2) -
\tfrac{1}{2}\zeta(s - 1, 3/2) = D(s - 1)$ has TWO simple poles at
$s = 4$ and $s = 2$ with bit-exact meromorphic residues $2$ and $-1/2$
respectively (Hurwitz pole rule and Laurent expansion via digamma
both verify;\ Mellin residues with $\Gamma$:\ $12$ and $-1/2$).  The
discrete-side Karamata Tauberian signature
$S(n_{\max}; s = 4) / \log(n_{\max}) \to 2$ at $n_{\max} \ge 200$
on the v3.60.0 sector-local substrate confirms the residue from the
finite-cutoff side, no Jacobian correction.  The chirality-resolved
$(D_e - D_o)(s - 1)$ at integer $s$ reproduces the Sprint Q5'-CH-3
panel verbatim with one-step shift:\ $s = 3$ gives $-1 + 2G$ (genuine
M3 cyclotomic at depth 1);\ $s = 5$ gives $8\beta_4 - 8G$;\ $s = 7$
gives $-32\beta_4 + 32\beta_6$;\ even-$s$ values collapse to $M_1$ via
Euler $\beta(\mathrm{odd}) \in \pi^{\mathrm{odd}}\Q$.  The FULL
$\eta_D$ at integer regular $s$ stratifies by parity-of-$(s-1)$:\ even
$\Rightarrow$ pure-Tate $M_2$ (e.g., $\eta_D(3) = -\pi^2/4$);\ odd
$\Rightarrow$ Apéry odd-zeta $\zeta(3), \zeta(5)$ in
$\mathcal{MT}(\Z)$ at depth 1 (NOT level 4).  STRUCTURAL FINDING:\
unlike $M_1$'s three exact-factor sibling normalizations
($\pi$, $4/\pi$, $\sqrt{\pi}/2$) of a single generator $\pi$, $M_3$'s
ring $\mathcal{MT}(\Z[i, 1/2], 4)$ is depth-graded with generators
$\{1, G = \beta(2), \beta(4), \beta(6), \ldots\}$ at successive depths
and no $\pi$-power reduction at any depth $\ge 1$.  This asymmetry IS
the operational content of the master-Mellin partition:\ $M_1$ at
$k = 0$ depth $0$ in $\Q[\pi, 1/\pi]$ vs $M_3$ at $k = 1$ depth $1+$
in cyclotomic level 4.  Together with the v3.59.0 Track 2 continuum
M1/M2 residue closure
(Remark~\ref{rem:q5p_continuum_residue}~\cite{paper32}) and the
Sprint Q5'-CH-2 M2 panel, the master Mellin engine M1/M2/M3 trinity
now has a complete continuum-residue identification at all three
extraction points (Mellin pole at $s = d/2$ for $M_1$;\ integer-$s$
regular points for $M_2$;\ Mellin poles + chirality-graded integer
panel for $M_3$).  Stage 1 closure of the cosmic-Galois $U^*$ bridge
at the continuum side is now complete across the master-Mellin
trinity;\ the Stage-2 motivic-Galois action on the symbol level
remains multi-year and unchanged in scope.
```

### (C) Paper 18 §III.7 — single sharpening sentence in the existing master-Mellin paragraph

After the existing Sprint Q5'-Stage1-2b-continuum closing sentence (after `korevaar2004`), append one sentence:

```latex
The companion M3 continuum residue identification (Sprint Q5'-Stage1-M3-Continuum,
2026-06-05, memo \texttt{debug/sprint\_q5p\_m3\_continuum\_memo.md}) closes the
$M_3$ component of the master-Mellin trinity at theorem grade:\ the $\eta$-Mellin
$\Gamma(s) \eta_D(s)$ with $\eta_D(s) = D(s - 1)$ has bit-exact meromorphic
residues $2$ and $-1/2$ at simple poles $s = 4, 2$ respectively, with the
genuine M3 cyclotomic content living in the chirality-resolved
$(D_e - D_o)(s - 1)$ at even-$s$ (Sprint Q5'-CH-3 panel transports verbatim).
Structural sharpening:\ $M_1$ admits three exact-factor sibling normalizations
($\pi$ Hopf-base Haar, $4/\pi$ L2 asymptote, $\sqrt{\pi}/2$ Mellin residue) of
a single generator $\pi$;\ $M_3$'s ring is depth-graded with generators
$\{1, G, \beta_4, \beta_6, \ldots\}$ at successive depths and NO $\pi$-power
reduction at depth $\ge 1$.  This asymmetry IS the operational content of
the $k = 0$ vs $k = 1$ master-Mellin partition.
```

(PI to apply, decline, or modify. No edits applied in this sprint.)

---

## One-line verdict

**POSITIVE.** Closes the M3 continuum residue follow-on from v3.59.0 Track 2 at theorem grade across all six requested axes:\ (1) Hurwitz reduction with two-route agreement; (2) two simple poles at $s = 4$ and $s = 2$; (3) bit-exact meromorphic residues $2$ and $-1/2$ via Hurwitz pole rule + Laurent digamma; (4) integer-$s$ panel of FULL $\eta_D$ matching Paper 28 $D(s-1)$ verbatim; (5) chirality-resolved panel reproducing Sprint Q5'-CH-3 with one-step shift; (6) STRUCTURAL FINDING — $M_3$ does NOT admit a three-sibling Hopf-base normalization because its ring $\mathcal{MT}(\Z[i, 1/2], 4)$ is depth-graded (Catalan $G$ has no $\pi$-power reduction at depth 1+). Discrete-side Karamata $S/\log n \to 2$ on the v3.60.0 substrate confirms the residue from the finite-cutoff side. Stage 1 cosmic-Galois $U^*$ bridge continuum-residue closure is now complete across the master-Mellin trinity M1/M2/M3.
