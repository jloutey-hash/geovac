# Sprint Q5'-Stage1-2a-JLO-nmax-sweep — bit-exact $\omega^{\mathrm{tri}}$ symbol across $n_{\max} \in \{1, 2, 3, 4\}$

**Date:** 2026-06-05
**Sprint:** Q5' Stage 1, Sub-Sprint 2a (cutoff sweep of the Sub-Sprint 1 JLO+CM-residue symbol)
**Driver:** `debug/compute_jlo_nmax_sweep.py`
**Data:** `debug/data/sprint_q5p_2a_jlo_nmax_sweep_data.json`
**Wall time:** 0.6 s (n_max=4 dominates at 0.44 s; n_max=1..3 negligible)
**Discipline:** bit-exact `sympy.Rational` throughout; no PSLQ; no floats; no transcendentals.

---

## TL;DR

**Verdict: CLEAN-POSITIVE.** The $\omega^{\mathrm{tri}}$ symbol triple $(M_1, M_2, M_3)$ of Sub-Sprint 1 extends bit-exactly across all four cutoffs. Closed-form polynomial-in-$n_{\max}$ expressions identified for each component (RE-DERIVED from explicit shell-degeneracy sums, NOT curve-fit). Pro-system rationality structure holds — every symbol value is a positive integer at every finite cutoff, no transcendentals appear at any $n_{\max}$.

The three closed forms:

$$
\boxed{
\;M_1(n_{\max}) \;=\; \dim\mathcal{H} \;=\; \frac{2\,n_{\max}(n_{\max}+1)(n_{\max}+2)}{3}\;
}
$$

$$
\boxed{
\;M_2(n_{\max}) \;=\; \mathrm{Tr}(\Lambda^2) \;=\; \frac{n_{\max}(n_{\max}+1)(n_{\max}+2)(2n_{\max}+1)(2n_{\max}+3)}{10}\;
}
$$

$$
\boxed{
\;M_3(n_{\max}) \;=\; \mathrm{Tr}(\gamma\Lambda) \;=\; \frac{n_{\max}(n_{\max}+1)^2(n_{\max}+2)}{2}\;
}
$$

Asymptotically: $M_1 \sim (2/3) n_{\max}^3$, $M_2 \sim (2/5) n_{\max}^5$, $M_3 \sim (1/2) n_{\max}^4$ — matching the JLO memo's named asymptotics (lines 190-194).

---

## Verdict against decision gate

| Gate | Verdict |
|:-----|:--------|
| CLEAN-POSITIVE | **selected** — symbol extends bit-exactly to $n_{\max} \in \{3, 4\}$; all three closed-form polynomial-in-$n_{\max}$ expressions are identified and verified at $n_{\max} \in \{1, 2, 3, 4\}$; pro-system rationality structure holds (all values integer). |
| POSITIVE-WITH-STRUCTURAL-FINDING | not selected — no new structural pattern emerged beyond the JLO memo's named asymptotics. The closed-form factorizations are richer than naive leading terms (each component has explicit polynomial factorization revealing common factors), but this is structural confirmation, not a new finding. |
| PARTIAL | not selected — even $n_{\max}=4$ completed in 0.44 s. |
| BLOCKED | rejected — bit-exact rational arithmetic closes without obstruction at every cutoff. |

---

## Computational summary

Bit-exact symbol triples at four cutoffs:

| $n_{\max}$ | $\dim\mathcal{H}$ | $\chi_+/\chi_-$ | $M_1$ | $M_2$ | $M_3$ | wall |
|:----------:|:-----------------:|:---------------:|:-----:|:-----:|:-----:|:----:|
| 1 | 4 | 2/2 | **4** | **9** | **6** | 0.00 s |
| 2 | 16 | 8/8 | **16** | **84** | **36** | 0.01 s |
| 3 | 40 | 20/20 | **40** | **378** | **120** | 0.09 s |
| 4 | 80 | 40/40 | **80** | **1188** | **300** | 0.42 s |

All values are positive integers. Cost scaling is roughly $O(N^3)$ in $\dim\mathcal{H}$ (one matrix multiplication for $D^2$ + scalar traces), well below sprint scale.

Sub-leading coefficients of $\phi_0^{\mathrm{odd}}(1; t)|_\Lambda$ expansion (rationals with denominators powers of 2):

| $n_{\max}$ | $c_0$ | $c_1$ | $c_2$ |
|:----------:|:-----:|:-----:|:-----:|
| 1 | 4 | -9 | 81/8 |
| 2 | 16 | -84 | 489/2 |
| 3 | 40 | -378 | 8181/4 |
| 4 | 80 | -1188 | 20493/2 |

$c_2 = +\mathrm{Tr}(\Lambda^4)/2!$ matches CH-1's $\mathrm{Tr}(\Lambda^4)$ at $n_{\max}=2,3,4$ (CH-1 reports $489$, $8181/2$, $20493$ — exactly $2 c_2$ in each case).

---

## Closed-form identification (curve-fit-audit-compliant)

**Methodology:** the closed forms were **RE-DERIVED** from the explicit shell-degeneracy sum given in the Sub-Sprint 1 JLO memo (lines 190-194), NOT fitted to four data points. Each derivation is reproduced below.

### $M_1$: $\dim\mathcal{H}$

**Starting point:** at each shell $n$, the chirality-balanced Camporesi-Higuchi degeneracy is $2n(n+1)$ (n(n+1) at each chirality). Summing:

$$
\dim\mathcal{H} \;=\; \sum_{n=1}^{n_{\max}} 2n(n+1) \;=\; 2\!\left[\frac{N(N+1)(2N+1)}{6} + \frac{N(N+1)}{2}\right] \;=\; \frac{N(N+1)(2N+4)}{3} \;=\; \frac{2N(N+1)(N+2)}{3},
$$

with $N \equiv n_{\max}$. **Verified bit-exact at all four cutoffs.** No free parameter; no fitting.

### $M_3$: $\mathrm{Tr}(\gamma\Lambda)$

**Starting point:** $\gamma$ flips the sign of $\Lambda$, so $\gamma\Lambda$ has eigenvalue $+(n+1/2) = +(2n+1)/2$ on every chirality block (the negative-chirality states sit at $-\lambda$ but $\gamma$ contributes $-1$, recovering $+\lambda$). Multiplicity at $|\lambda_n| = (n+1/2)$ is $2n(n+1)$. Summing:

$$
\mathrm{Tr}(\gamma\Lambda) \;=\; \sum_{n=1}^{n_{\max}} 2n(n+1)\!\left(n + \tfrac{1}{2}\right) \;=\; \sum_{n=1}^{n_{\max}} n(n+1)(2n+1).
$$

Using $n(n+1)(2n+1) = 2n^3 + 3n^2 + n$ and Faulhaber:

$$
\sum_{n=1}^{N} n(n+1)(2n+1) \;=\; 2\!\left[\tfrac{N(N+1)}{2}\right]^2 + \tfrac{3 N(N+1)(2N+1)}{6} + \tfrac{N(N+1)}{2} \;=\; \tfrac{N(N+1)}{2}\!\left[N(N+1) + (2N+1) + 1\right] \;=\; \tfrac{N(N+1)(N+1)(N+2)}{2}.
$$

So $M_3 = N(N+1)^2(N+2)/2$. **Verified bit-exact:** $N=1$: $1 \cdot 4 \cdot 3/2 = 6$; $N=2$: $2 \cdot 9 \cdot 4 /2 = 36$; $N=3$: $3 \cdot 16 \cdot 5/2 = 120$; $N=4$: $4 \cdot 25 \cdot 6/2 = 300$.

### $M_2$: $\mathrm{Tr}(\Lambda^2)$

**Starting point:** $\Lambda^2$ has eigenvalue $(n+1/2)^2 = n^2 + n + 1/4$ at multiplicity $2n(n+1)$. Summing:

$$
\mathrm{Tr}(\Lambda^2) \;=\; \sum_{n=1}^{n_{\max}} 2n(n+1)\!\left(n + \tfrac{1}{2}\right)^2.
$$

Expand $(n+1/2)^2 = (n^2 + n) + 1/4$, so $2n(n+1)(n+1/2)^2 = 2n(n+1)\cdot[n(n+1) + 1/4] = 2 n^2(n+1)^2 + n(n+1)/2$. Then:

$$
\mathrm{Tr}(\Lambda^2) \;=\; 2 \sum_{n=1}^{N} n^2 (n+1)^2 \;+\; \tfrac{1}{2}\sum_{n=1}^{N} n(n+1).
$$

Using $n^2(n+1)^2 = n^4 + 2n^3 + n^2$ and the standard Faulhaber identities $\sum n^4 = N(N+1)(2N+1)(3N^2+3N-1)/30$, $\sum n^3 = [N(N+1)/2]^2$, $\sum n^2 = N(N+1)(2N+1)/6$, $\sum n(n+1) = N(N+1)(N+2)/3$, the sum factors (sympy.factor) into:

$$
\mathrm{Tr}(\Lambda^2) \;=\; \frac{N(N+1)(N+2)(2N+1)(2N+3)}{10}.
$$

**Verified bit-exact:** $N=1$: $1 \cdot 2 \cdot 3 \cdot 3 \cdot 5/10 = 90/10 = 9$; $N=2$: $2 \cdot 3 \cdot 4 \cdot 5 \cdot 7/10 = 840/10 = 84$; $N=3$: $3 \cdot 4 \cdot 5 \cdot 7 \cdot 9/10 = 3780/10 = 378$; $N=4$: $4 \cdot 5 \cdot 6 \cdot 9 \cdot 11/10 = 11880/10 = 1188$.

### Curve-fit audit (per [[feedback_audit_numerical_claims]])

- **Free-parameter count:** zero. Each closed form is derived structurally from the CH-1 shell degeneracy $2n(n+1)$ and the CH eigenvalue $(n+1/2)$ — no fitted constants.
- **Selection-bias check:** the closed forms were proposed in the JLO memo (line 191) BEFORE the n_max=3,4 computations were performed in this sprint. The polynomial degrees (3 for $M_1$, 5 for $M_2$, 4 for $M_3$) match the leading asymptotic claims from that memo. The full polynomial structure (not just leading term) is a richer-than-named outcome.
- **Alternatives considered:** $M_1$ could in principle be $a N^3 + b N^2 + c N + d$ with four free parameters — at four data points, this is exactly determined regardless of mechanism. **However**, our derivation produces $M_1$ from the structural sum $\sum 2n(n+1)$, giving a single specific polynomial (factored as $\tfrac{2}{3} N(N+1)(N+2)$) without any fitted parameter. The 4-point match is a verification, not an inference.
- **Robustness:** the closed forms predict $M_1(5) = 140$, $M_3(5) = 270$, $M_2(5) = 5 \cdot 6 \cdot 7 \cdot 11 \cdot 13 / 10 = 30030/10 = 3003$ as falsifiable extrapolations. The driver can be extended to $n_{\max}=5$ in $\sim$5 s.
- **Independent test:** all three closed forms cross-check against the CH-1 panel (Sub-Sprint 1 memo lines 76-97) bit-exactly at $n_{\max} \in \{2, 3, 4\}$ (CH-1 reports the same integers via a different code path, `compute_ch_k_nmax_truncated.py`).

**Conclusion:** the closed forms come from re-derivation, not fitting. CLEAN-POSITIVE on curve-fit audit.

---

## Cross-check against CH-1 panel (memo lines 86, 96)

Sub-Sprint 1's CH-1 panel reports the relevant moments at $n_{\max} \in \{2, 3, 4\}$:

| $n_{\max}$ | CH-1 $\dim\mathcal{H}$ | CH-1 $\mathrm{Tr}(\Lambda^2)$ | CH-1 $\mathrm{Tr}(\gamma\Lambda)$ | this sprint $M_1$ | $M_2$ | $M_3$ |
|:----------:|:----------------------:|:-----------------------------:|:---------------------------------:|:-----------------:|:-----:|:-----:|
| 2 | 16 | 84 | 36 | 16 ✓ | 84 ✓ | 36 ✓ |
| 3 | 40 | 378 | 120 | 40 ✓ | 378 ✓ | 120 ✓ |
| 4 | 80 | 1188 | 300 | 80 ✓ | 1188 ✓ | 300 ✓ |

**All nine cells match bit-exactly across the two independent code paths.**

---

## Pro-system rationality structure

The pro-system $\{\mathcal{T}_{n_{\max}}\}_{n_{\max} \ge 1}$ carries the symbol triple $\omega^{\mathrm{tri}}(\mathcal{T}_{n_{\max}}) = (M_1, M_2, M_3)(n_{\max}) \in \mathbb{Z}^3_{>0}$. Three structural observations:

1. **Integrality.** At every finite cutoff, all three components are positive integers. No half-integers, no $\sqrt{2}$-style algebraic numbers, no transcendentals. This is consistent with the Layer-1 discrete-for-skeleton discipline ([[feedback_discrete_for_skeleton]]): the bare graph is integer-rational by construction; transcendentals appear only when Mellin-regulated periods are taken.

2. **Polynomial growth.** $M_k$ is a polynomial of degree $k+2$ in $n_{\max}$:
   - $M_1$: cubic, leading $(2/3) n_{\max}^3$ — heat-kernel scaling matches expected Weyl law for $S^3$ (volume scales as $\dim\mathcal{H}^{3/3} \sim n_{\max}^3$).
   - $M_2$: quintic, leading $(2/5) n_{\max}^5$ — matches $\dim\mathcal{H} \cdot \langle\lambda^2\rangle$ where $\langle\lambda^2\rangle \sim n_{\max}^2$.
   - $M_3$: quartic, leading $(1/2) n_{\max}^4$ — matches $\dim\mathcal{H} \cdot \langle\lambda\rangle_+$ where $\langle\lambda\rangle_+ \sim n_{\max}$.

3. **Common factor structure.** Both $M_1$ and $M_2$ contain the factor $N(N+1)(N+2)$ (the cubic "binomial-coefficient $\binom{N+2}{3}$" factor); $M_2$ multiplies this by $(2N+1)(2N+3)$; $M_3$ multiplies $N(N+1)(N+2)$ by $(N+1)$. The recurring $N(N+1)(N+2)/6 = \binom{N+2}{3}$ factor is the **GeoVac dimension generator** — appears in $\dim\mathcal{H}$ as $4\binom{N+2}{3}$, in $M_2$ as $\tfrac{6}{10}(2N+1)(2N+3)\binom{N+2}{3}$, in $M_3$ as $3(N+1)\binom{N+2}{3}$.

The pro-system data $\{M_k(n_{\max})\}$ diverges as $n_{\max} \to \infty$ (each component is a positive-coefficient polynomial). The continuum limit therefore lives in the **Mellin-regulated** quantities $\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})](s)$, not in the raw symbol triple. The transition from polynomial growth to convergent period is Sub-Sprint 2b territory (formal Mellin transform).

---

## Implications for Sub-Sprints 2b and 2c

### Sub-Sprint 2b — Mellin transform of $\phi_0^{\mathrm{odd}}$ (medium, ~2 weeks)

Now that $\phi_0^{\mathrm{odd}}(1; t)|_\Lambda$ is bit-exact at four cutoffs with full sub-leading data:

| $n_{\max}$ | $c_0$ | $c_1$ | $c_2$ |
|:----------:|:-----:|:-----:|:-----:|
| 1 | 4 | -9 | 81/8 |
| 2 | 16 | -84 | 489/2 |
| 3 | 40 | -378 | 8181/4 |
| 4 | 80 | -1188 | 20493/2 |

The formal Mellin transform $\int_0^\infty t^{s-1} \phi_0^{\mathrm{odd}}(1; t) dt$ requires the FULL series in $t$, which on the finite-dimensional Hilbert space is a polynomial in $t$ (eigenvalue decay times $e^{-tD^2}$, but on finite dimension the sum is finite). The natural Mellin object is therefore the spectral zeta $\zeta_{D^2}(s) = \mathrm{Tr}(D^{-2s})$ on the truncated spectrum — exactly the M2 master Mellin engine slot. **At $s = -1$ (analytic continuation), $\zeta_{D^2}(-1) = \mathrm{Tr}(D^2)$ matches $M_2$ bit-exactly.** This is the bridge to Paper 28 Theorem T9.

### Sub-Sprint 2c — bicomplex $(b, B)$ at $n_{\max} = 2$ (medium, ~1-2 weeks)

The sub-leading $c_1, c_2$ at $n_{\max} = 2$ give the heat-kernel expansion of $\phi_0^{\mathrm{odd}}$ to order $t^2$. The boundary maps $(b, B)$ act on this data via standard cyclic-cohomology formulae. The closed-form polynomial structure across cutoffs (this sprint) gives confidence that the bicomplex closure at $n_{\max} = 2$ extends predictably.

### Sub-Sprint 2d (NEW, named) — functoriality of restriction $\mathcal{T}_{n_{\max}+1} \to \mathcal{T}_{n_{\max}}$

The closed-form polynomial structure means the symbol shift $\omega^{\mathrm{tri}}(\mathcal{T}_{n_{\max}+1}) - \omega^{\mathrm{tri}}(\mathcal{T}_{n_{\max}})$ is itself a closed-form polynomial:
- $\Delta M_1 = M_1(n+1) - M_1(n) = 2(n+1)(n+2)$
- $\Delta M_3 = M_3(n+1) - M_3(n) = (n+1)(n+2)(2n+3)$
- $\Delta M_2 = M_2(n+1) - M_2(n) = 2(n+1)(n+2)(n+3/2)^2 = (n+1)(n+2)(2n+3)^2/2$

These are the **shell-by-shell increments**, each a positive integer (or half-integer for $\Delta M_2$ in this factorized form, but a positive integer when computed as the difference of the closed forms). The shell-wise functoriality (a candidate ingredient for the pro-system functor) is exactly the closed-form increment. **Named as Sub-Sprint 2d**, opt-in.

---

## Honest scope

1. **This sprint extends from one cutoff to three more cutoffs.** The bit-exact $(M_1, M_2, M_3)$ symbol is now known at $n_{\max} \in \{1, 2, 3, 4\}$ with three verified closed-form polynomial expressions in $n_{\max}$.

2. **No continuum limit computed.** The closed-form polynomials DIVERGE as $n_{\max} \to \infty$ (positive-coefficient polynomials in $n_{\max}$). The continuum object is the Mellin-regulated $\zeta_{D^2}(s)$ / heat-kernel coefficient ring, not the raw symbol triple. This is named Sub-Sprint 2b.

3. **No pro-limit functor constructed.** The pro-system $\{\mathcal{T}_{n_{\max}}\}$ has a clean cohomology-level pro-system (each $\omega^{\mathrm{tri}}(\mathcal{T}_{n_{\max}})$ is integer-valued and polynomial), but the *enriched fiber functor* on the Tannakian category requires a tensor structure across cutoffs that this sprint does NOT construct. Sub-Sprint 2c / Sub-Sprint 2d territory.

4. **CLEAN-POSITIVE on curve-fit audit.** The three closed forms are RE-DERIVED from the shell-degeneracy sum, NOT fitted to 4 data points. Selection-bias zero; free parameters zero; independent code path (CH-1 panel) confirms the integer values bit-exactly at $n_{\max} \in \{2, 3, 4\}$.

5. **No transcendentals introduced.** Zero $\pi$, zero $\zeta$, zero $G$, zero $\log$. The Layer-1 skeleton remains rational-integer at every $n_{\max}$. Tagging via Paper 18 §III.7 is therefore inapplicable at this stage; transcendentals enter only at Sub-Sprint 2b when the Mellin transform is taken against $\Gamma(s)$.

6. **Cost trivial.** Total wall time 0.6 s for the four-cutoff sweep; $n_{\max} = 5$ extension would add $\sim$5 s ($O(N^3)$ matrix multiplication for $D^2$ on $\dim\mathcal{H}_5 = 140$).

7. **WH1 PROVEN unchanged.** This sprint adds cohomological/symbolic data on top of WH1 PROVEN's propinquity-side foundation; it does not test or extend WH1 itself.

8. **No paper edits applied.** Recommendations below.

---

## Files

### Produced
- `debug/compute_jlo_nmax_sweep.py` — driver ($\sim$420 lines, $\sim$0.6 s wall).
- `debug/data/sprint_q5p_2a_jlo_nmax_sweep_data.json` — bit-exact rational dump ($\sim$5 KB) with symbol triples, closed forms, cross-checks at all four $n_{\max}$.
- `debug/sprint_q5p_2a_jlo_nmax_sweep_memo.md` — this memo.

### Used (load-bearing inputs)
- `debug/sprint_q5p_jlo_nmax2_memo.md` — Sub-Sprint 1 closure: explicit $\omega^{\mathrm{tri}}$ at $n_{\max} = 2$ and the named asymptotic closed-form expressions (lines 190-194).
- `debug/sprint_q5p_jlo_nmax2_data.json` — bit-exact reference for $n_{\max} = 2$ row.
- `debug/sprint_q5p_ch1_memo.md` — CH-1 panel (lines 76-97) providing independent verification at $n_{\max} \in \{2, 3, 4\}$.
- `debug/compute_ch_k_nmax_truncated.py` — CH-1 driver (independent code path for cross-check).
- `geovac/spectral_triple.py` — `FockSpectralTriple` class providing $\Lambda$, $\gamma$, $D$ in exact `sympy.Rational`.

---

## Recommended paper edits (PI to apply, decline, or modify)

### Paper 55 §subsec:open_m2_m3 — append the cutoff-sweep closure to the Stage-1 paragraph

The Sub-Sprint 1 paragraph (the Stage-1 explicit-symbol paragraph already in Paper 55 §subsec:open_m2_m3) ends with "Sub-Sprint 2 (continuum limit, functoriality across cutoffs, Tannakian-category construction) remains multi-week to multi-year." A one-paragraph follow-on closing Sub-Sprint 2a:

> *Stage 1 cutoff-sweep extension (Sprint Q5'-Stage1-2a-JLO-nmax-sweep, June 2026; memo \texttt{debug/sprint\_q5p\_2a\_jlo\_nmax\_sweep\_memo.md}; data \texttt{debug/data/sprint\_q5p\_2a\_jlo\_nmax\_sweep\_data.json}).* The $\omega^{\mathrm{tri}}$ symbol $(M_1, M_2, M_3)$ extends bit-exactly to $n_{\max} \in \{1, 2, 3, 4\}$ as
> \begin{align*}
> M_1(n_{\max}) &= \dim\mathcal{H} = \tfrac{2 n_{\max}(n_{\max}+1)(n_{\max}+2)}{3}, \\
> M_2(n_{\max}) &= \mathrm{Tr}(\Lambda^2) = \tfrac{n_{\max}(n_{\max}+1)(n_{\max}+2)(2n_{\max}+1)(2n_{\max}+3)}{10}, \\
> M_3(n_{\max}) &= \mathrm{Tr}(\gamma\Lambda) = \tfrac{n_{\max}(n_{\max}+1)^2(n_{\max}+2)}{2},
> \end{align*}
> giving the table $\{(4, 9, 6), (16, 84, 36), (40, 378, 120), (80, 1188, 300)\}$ at $n_{\max} = 1, 2, 3, 4$. Each closed form is RE-DERIVED from the chirality-balanced CH shell degeneracy $2n(n+1)$ and the CH eigenvalue $(n+1/2)$, not curve-fit; no free parameters. The pro-system $\{\mathcal{T}_{n_{\max}}\}$ carries integer-valued symbol triples at every cutoff with no transcendentals. The continuum limit and the pro-limit functor remain Sub-Sprints 2b/2c (Mellin transform, bicomplex closure, Tannakian-tensor structure).

### Paper 32 §VIII — no edit needed

The Sub-Sprint 1 Remark `rem:stage1_witness` recommended in the JLO memo (and not yet applied) would already cover Sub-Sprint 2a if applied. No new Paper 32 edit is needed for this sprint.

(Recommendations only; no edits applied.)

---

## One-line verdict

**CLEAN-POSITIVE.** The $\omega^{\mathrm{tri}}$ symbol of Sub-Sprint 1 extends bit-exactly to $n_{\max} \in \{1, 2, 3, 4\}$ with three structurally-derived closed-form polynomial-in-$n_{\max}$ expressions (degree 3, 5, 4 respectively); pro-system rationality structure holds; no transcendentals; 0.6 s wall; CH-1 cross-check bit-exact on 9/9 cells. Stage-1 first span of the cosmic-Galois $U^*$ bridge now has explicit polynomial pro-system data across four cutoffs.
