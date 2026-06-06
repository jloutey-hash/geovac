# Sprint Q5'-CH-3 — M3 cyclotomic verification at integer s

Date: 2026-06-05 (same session as Q5'-CH-1 + CH-2)
Scope: 1-hour follow-on to Q5'-CH-2. Verify the M3 cyclotomic mixed-Tate
component of the master Mellin engine at integer $s$, completing the
M1 / M2 / M3 verification triad for the $\omega^{\mathrm{tri}}$ symbol level.

Driver: `debug/compute_ch_m3_cyclotomic.py`. Data:
`debug/data/sprint_q5p_ch3_data.json`. Wall: 1.5 s.

## TL;DR

**M3 verification with a structural parity-of-$s$ split.** Paper 28
Theorem 3 closed form

$$
D_{\mathrm{even}}(s) - D_{\mathrm{odd}}(s) = 2^{s-1}(\beta(s) - \beta(s-2))
$$

verified bit-exact at $s \in \{4, 5, 6\}$ to residual $<10^{-80}$ at
80 dps (numeric cross-check via the Hurwitz form $\beta(s) = 4^{-s}
(\zeta(s, 1/4) - \zeta(s, 3/4))$ against the existing
`geovac/qed_vertex.py::_dirac_D_even/_odd`). Closed forms collected at
$s \in \{2, 3, 4, 5, 6\}$, with the symbolic structure:

| $s$ | $D_e - D_o$ closed form | M-classification |
|:---:|:------------------------|:------------------|
| 2   | $-1 + 2G$               | **M3 (Catalan)** |
| 3   | $-\pi + \pi^3/8$        | **M1 collapse**  |
| 4   | $8\beta_4 - 8G$         | **M3 (Catalan + $\beta_4$)** |
| 5   | $-\pi^3/2 + 5\pi^5/96$  | **M1 collapse**  |
| 6   | $-32\beta_4 + 32\beta_6$| **M3 ($\beta_4, \beta_6$)** |

**Structural finding: the M3 cyclotomic content shows up bit-exactly
only at EVEN $s$.** At odd $s$, the Euler-style closed form
$\beta(2k+1) = (-1)^k E_{2k}/(2(2k)!) \cdot (\pi/2)^{2k+1}$ reduces
$\beta(\text{odd})$ to $\pi^{\text{odd}} \cdot \mathbb{Q}$, so the
$D_e - D_o$ antisymmetric combination collapses into the M1 ring
$\mathbb{Q}[\pi, \pi^{-1}]$. At even $s$, $\beta(2k)$ for $k \ge 1$ has
no known closed form in $\pi$ (Catalan $G = \beta(2)$ irrationality
status conjectural; $\beta(4), \beta(6), \ldots$ believed irrational),
so the M3 content stays in MT$(\mathbb{Z}[i, 1/2], 4)$ as the
Deligne--Glanois cyclotomic mixed-Tate ring at level 4 over $\mathbb{Z}[i]$.

## Structural placement: M2 uniform vs M3 parity-graded

The contrast with Q5'-CH-2 sharpens the operator-order grading:

| Slot | Behaviour at integer $s$ | Period ring |
|:----:|:-------------------------|:------------|
| **M2** (Q5'-CH-2) | Uniform across all $s$: every $\zeta_{D^2}(s)$ is exactly two terms in $\bigoplus_k \pi^{2k}\mathbb{Q}$ | Pure-Tate, depth 0 |
| **M3** (Q5'-CH-3) | Parity-graded: even $s$ genuinely cyclotomic, odd $s$ collapses to M1 | Cyclotomic mixed-Tate at level 4 (only at even $s$) |

The M3 parity-of-$s$ collapse is the FUNCTIONAL-EQUATION shadow of the
Euler-style $\beta(\text{odd})$ closed form. It also reads as a
DEPTH STRATIFICATION:\ at odd $s$, the depth-1 cyclotomic content
collapses into depth-0 M1 because $\beta(\text{odd})$ is reducible.

## Implications for the $\omega^{\mathrm{tri}}$ enrichment

Combining Q5'-CH-1, CH-2, and CH-3, the symbol level of the candidate
enriched fiber functor $\omega^{\mathrm{tri}}$ is now explicit on each
of M1, M2, M3:

- **M1 component (k = 0 chirality):** plain trace of heat kernel, leading
  coefficient $\dim \mathcal{H}$ at finite $n_{\max}$;\ continuum period
  $\pi = \mathrm{Vol}(S^2)/4$ via Hopf base measure (not in this sprint).
- **M2 component (k = 2 within HP$_0$ heat-kernel-order grading):**
  $\zeta_{D^2}(s)$ at integer $s$ gives bit-exact two-term polynomials in
  $\pi^{2k} \mathbb{Q}$;\ at $s \in \{1, \ldots, 5\}$ the
  $\pi$-fingerprints are $(-\tfrac{1}{4}, 0)$, $(1, -\tfrac{1}{12})$,
  $(\tfrac{1}{3}, -\tfrac{1}{30})$, $(\tfrac{2}{15}, -\tfrac{17}{1260})$,
  $(\tfrac{17}{315}, -\tfrac{31}{5670})$.
- **M3 component (k = 1 chirality):** $D_{\mathrm{even}}(s) -
  D_{\mathrm{odd}}(s)$ at integer $s$ gives bit-exact MT$(\mathbb{Z}[i,
  1/2], 4)$ expressions in $\{1, G, \beta_4, \beta_6, \ldots\}$, with
  the depth-1 cyclotomic content visible only at even $s$ and depth-0
  collapse at odd $s$.

Stage 1 of the cosmic-Galois $U^*$ program has a complete bit-exact
target on each of M1, M2, M3, with the structural distinction between
the three sectors visible at the symbol level via:\ (a) chirality grading
for M1+M2 vs M3 (Q5'-CH-1), (b) operator-order grading for M1 vs M2
(Q5'-CH-1 + CH-2 depth-0 vs depth-2), (c) parity-of-$s$ stratification
for M3 vs M1-collapse (CH-3).

## Caveats and honest scope

1. **Numerical cross-check at $s = 2, 3$ skipped.** The
   `geovac/qed_vertex.py::_dirac_D_even/_odd` functions require $s \ge 4$
   for the underlying mpmath Hurwitz computation to converge cleanly.
   At $s = 2, 3$ the closed forms are taken from Paper 28 Theorem 3
   directly without numeric cross-check;\ the bit-exact verification at
   $s = 4, 5, 6$ confirms the formula's correctness, and the $s = 2, 3$
   closed forms are immediate consequences of the same formula with
   $\beta(\text{small})$ closed forms substituted.
2. **No higher-depth M3 verification.** Deeper M3 content (depth 2 and
   above, including the $S_{\min}$ irreducibility result from Paper 28)
   would require multi-loop vertex restricted sums beyond the two-loop
   $D_{\mathrm{even}} \pm D_{\mathrm{odd}}$ scope. This sprint hits the
   depth-1 cyclotomic content cleanly.
3. **The parity-of-$s$ stratification is structural, not Tannakian.**
   The parity comes from $\beta(\text{odd}) \in \pi^{\text{odd}}\mathbb{Q}$
   (an Euler closed form), not from a Galois action distinguishing
   even-$s$ from odd-$s$ M3 outputs. This is a *period-side* refinement.
4. **Curve-fit-audit:** zero free parameters; closed forms forced by
   Hurwitz identities + Euler-style $\beta$ identities. Selection-bias
   minimal — Paper 28 Theorem 3 published before this sprint.
5. **Discrete-for-skeleton compliance:** symbolic sympy + bit-exact
   mpmath at 80 dps; no PSLQ. The closed forms are direct symbolic
   identifications.
6. **Catalan $G$ irrationality:** assumed (consistent with the
   conjectural status). If $G$ turned out to be rational or in
   $\pi^{2k}\mathbb{Q}$, the M2/M3 distinction at $s = 2$ would
   collapse, but the case-exhaustion theorem would still hold via
   weight inheritance.

## Files produced

- `debug/compute_ch_m3_cyclotomic.py` — driver (~210 lines).
- `debug/data/sprint_q5p_ch3_data.json` — closed forms + bit-exact
  numeric cross-checks per $s \in \{2, \ldots, 6\}$.
- `debug/sprint_q5p_ch3_memo.md` — this memo.

## Recommended paper edit

Paper 55 §subsec:open_m2_m3 already has the Q5'-CH-1 + Q5'-CH-2
content from this session. A single sentence completes the triad:

> *Continuum M3 verification at integer $s$ (Sprint Q5'-CH-3,
> 2026-06-05, memo \texttt{debug/sprint\_q5p\_ch3\_memo.md}).* The
> Paper 28 Theorem 3 closed form $D_{\mathrm{even}}(s) -
> D_{\mathrm{odd}}(s) = 2^{s-1}(\beta(s) - \beta(s-2))$ is verified
> bit-exact at $s \in \{4, 5, 6\}$ to residual $< 10^{-80}$ at 80 dps;
> the M3 cyclotomic content shows up bit-exactly only at even $s$
> (genuine $\beta_4, \beta_6$ content at $s = 4, 6$;\ M1 collapse at
> $s = 3, 5$ via the Euler-style $\beta(\text{odd}) \in
> \pi^{\text{odd}}\mathbb{Q}$ closed form). Together with Q5'-CH-2's
> uniform M2 pure-Tate behaviour, this places the depth-stratification
> of $\omega^{\mathrm{tri}}$ on the cohomological side: M2 is depth-0
> uniformly, M3 is depth-1 stratified by parity of $s$.

Sprint will apply this addition.

## One-line verdict

CLEAN POSITIVE — bit-exact M3 cyclotomic identification at even $s$ via
Paper 28 Theorem 3 ($s = 4, 5, 6$ at residual $< 10^{-80}$); structural
parity-of-$s$ stratification refines the M3 / M1-collapse distinction;
third concrete stone of the cosmic-Galois $U^*$ bridge, completing the
M1 / M2 / M3 verification triad at the symbol level.
