# Sprint Q5'-CH-2 — Continuum spectral zeta and bit-exact M2 verification

Date: 2026-06-05 (same session as Q5'-CH-1)
Scope: 1-hour follow-on to Q5'-CH-1. Push the truncated CH master Mellin
data to the continuum and verify that the M2 component of the
$\omega^{\mathrm{tri}}$ fiber functor lands bit-exactly in
$\bigoplus_k \pi^{2k}\mathbb{Q}$ at integer $s$, as predicted by
Paper 32 §VIII Cor `cor:m2_mixed_tate`.

Driver: `debug/compute_ch_spectral_zeta_continuum.py`. Data:
`debug/data/sprint_q5p_ch2_data.json`. Wall: 0.5 s.

## TL;DR

**Five-for-five bit-exact verification of the M2 pure-Tate identification
at integer $s$.** Every continuum CH spectral zeta value at
$s \in \{1, 2, 3, 4, 5\}$ is exactly two terms in
$\bigoplus_k \pi^{2k}\mathbb{Q}$:

| $s$ | $\zeta_{D^2}(s)$ closed form | M2 fingerprint $(\pi^{2k}\text{ coeffs})$ |
|:---:|:----------------------------:|:------------------------------------------|
| 1   | $-\pi^2/4$                   | $\pi^2 : -\tfrac{1}{4}$                   |
| 2   | $\pi^2 - \pi^4/12$           | $\pi^2 : 1, \; \pi^4 : -\tfrac{1}{12}$    |
| 3   | $\pi^4/3 - \pi^6/30$         | $\pi^4 : \tfrac{1}{3}, \; \pi^6 : -\tfrac{1}{30}$ |
| 4   | $\tfrac{2}{15}\pi^6 - \tfrac{17}{1260}\pi^8$ | $\pi^6 : \tfrac{2}{15}, \; \pi^8 : -\tfrac{17}{1260}$ |
| 5   | $\tfrac{17}{315}\pi^8 - \tfrac{31}{5670}\pi^{10}$ | $\pi^8 : \tfrac{17}{315}, \; \pi^{10} : -\tfrac{31}{5670}$ |

Every closed form derived symbolically via the Hurwitz reduction

$$
\zeta_{D^2}(s) = 2\,\zeta(2s - 2,\, \tfrac{3}{2}) - \tfrac{1}{2}\,\zeta(2s,\, \tfrac{3}{2})
$$

combined with the standard identity $\zeta(s, 1/2) = (2^s - 1)\,\zeta(s)$
and the shift $\zeta(s, 3/2) = \zeta(s, 1/2) - 2^s$. Every nonzero
$\pi$-power has a rational coefficient with small denominator. The
verification uses zero free parameters and no PSLQ — closed forms
extracted by `sympy` polynomial collection.

## Structural finding: two-term exactness at integer $s$

**At every integer $s \geq 1$, the spectral zeta $\zeta_{D^2}(s)$ contains
exactly two non-zero $\pi^{2k}$ terms — no more, no fewer.** This matches
the Dirac SD two-term exactness on $S^3$ from Paper 51 Cor 2.1
($a_k^{D^2} = 0$ for $k \geq 2$ in the asymptotic heat-trace expansion).

The structural explanation: the CH spectral zeta reduces to a difference
of two Hurwitz zetas at shifted integer arguments, each of which is (via
the standard identity) a Riemann $\zeta(\text{even})$ times an integer
$2^j$ minus an integer. Riemann's $\zeta(2k) \in \pi^{2k}\mathbb{Q}$ then
forces each Hurwitz term into a single $\pi^{2k}\mathbb{Q}$ slot, and the
linear combination preserves the two-term shape.

This is a much sharper statement than "$\zeta_{D^2}(s)$ is in M2" — it's
"$\zeta_{D^2}(s)$ is in the M2 ring AT DEPTH TWO, with both depth slots
populated by rationals with denominators bounded by $\mathrm{poly}(s)$."

## Continuum vs truncated values

The truncated spectral zeta $\zeta_{n_{\max}}(s) = \sum_{n=1}^{n_{\max}}
2n(n+1)(n+1/2)^{-2s}$ converges to the continuum value for $s \geq 2$
(rate $\sim n_{\max}^{3 - 2s}$):

| $s$ | $n_{\max} = 2$ | $n_{\max} = 3$ | $n_{\max} = 4$ | continuum | tail at $n_{\max} = 4$ |
|:---:|:--------------:|:--------------:|:--------------:|:---------:|:---------------------:|
| 2   | 1.0973         | 1.2573         | 1.3548         | 1.7522    | $3.97 \times 10^{-1}$ |
| 3   | 0.4003         | 0.4134         | 0.4182         | 0.4234    | $5.20 \times 10^{-3}$ |
| 4   | 0.1639         | 0.1650         | 0.1652         | 0.1654    | $1.21 \times 10^{-4}$ |
| 5   | 0.0706         | 0.0707         | 0.0707         | 0.0707    | $3.32 \times 10^{-6}$ |

At $s = 1$ the truncated sums DIVERGE (each term $\sim 1$, the sum has a
pole of $\sum 1$ shape). The continuum value $-\pi^2/4$ is the analytic
continuation through Hurwitz zeta. This is the correct regularised value
(consistent with the Riemann zeta functional equation at $s = 0$:
$\zeta(0) = -1/2$, $\zeta(-2k) = 0$ for $k \geq 1$; the negative integer
values are the Bernoulli rational residues).

## Structural placement in the Q5' bridge

Sprint Q5'-CH-1 located the chirality-parity component of the master
Mellin k-slot at finite $n_{\max}$ in bit-exact rational form. Sprint
Q5'-CH-2 locates the heat-kernel-order M2 component at the continuum in
bit-exact $\pi^{2k}\mathbb{Q}$ form. Together they give the
$\omega^{\mathrm{tri}}$ fiber functor an explicit symbol-level realisation
at both ends:

- **At finite $n_{\max}$ (Q5'-CH-1):** $\omega^{\mathrm{tri}}$ sends the
  master Mellin source $\mathrm{Tr}(D^k\,e^{-tD^2})$ to a small-$t$
  rational Taylor series, with the chirality grading already factored.
  M1 and M2 share a chirality sector (HP$_0$); M3 occupies HP$_1$.
- **At the continuum Mellin transform (Q5'-CH-2):** the M2 component of
  $\omega^{\mathrm{tri}}$ lands in $\bigoplus_k \pi^{2k}\mathbb{Q}$ with
  two non-zero $\pi$-power slots per integer $s$. This is the bit-exact
  M2 period content predicted by the case-exhaustion theorem and
  Cor `cor:m2_mixed_tate`.

The natural Stage 1 of the cosmic-Galois program is to lift this
symbol-level data to a graded fiber functor on the dg-category. Q5'-CH-1
+ Q5'-CH-2 together supply the bit-exact targets:\ for each $k \in \{0,
2\}$, $\omega^{\mathrm{tri}}$ at integer $s$ produces a specific
two-term $\pi^{2k}\mathbb{Q}$ expression, and the operator-order
distinction $k = 0$ vs $k = 2$ is visible at the level of the leading
$\pi$-power (which depends linearly on $s$ via the Hurwitz shift).

## Caveats and honest scope

1. **No M3 verification this sprint.** M3 lives in cyclotomic mixed-Tate
   $\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$ via the vertex-restricted Dirichlet
   content (Paper 28 Thm 3: $D_{\mathrm{even}} - D_{\mathrm{odd}} = 2^{s-1}
   (\beta(s) - \beta(s-2))$). The CH spectral zeta evaluated at integer
   $s$ is the M2 sector; M3 lives in a different observable class.
2. **No M1 verification at the continuum level.** M1's $\pi =
   \mathrm{Vol}(S^2)/4$ enters through spectral action coefficient
   normalisations, not through $\zeta_{D^2}(s)$ directly. A separate
   M1 verification would target the spectral action or the L2 propinquity
   rate $4/\pi$.
3. **The pole at $s = 1$ is correctly handled** by Hurwitz analytic
   continuation, but the truncated sum does not converge there. The
   $-\pi^2/4$ value is the regularised analytic value, not a sum limit.
4. **Curve-fit-audit (per `feedback_audit_numerical_claims`):** zero free
   parameters; the closed forms are forced by Hurwitz identities and
   $\zeta(2k) \in \pi^{2k}\mathbb{Q}$. Selection-bias is minimal — the M2
   ring was named in Paper 32 §VIII before this sprint, and the prediction
   is bit-exact.
5. **Discrete-for-skeleton compliance:** symbolic sympy throughout; the
   M2 closed forms are exact closed-form identifications, not PSLQ +
   precision sweeps.
6. **Two-term-exactness has structural status** (Paper 51 Cor 2.1 Dirac
   SD), not just empirical pattern. The depth-2 structure is the
   continuum content of the operator-order grading at integer $s$.

## Files produced

- `debug/compute_ch_spectral_zeta_continuum.py` — driver (~170 lines).
- `debug/data/sprint_q5p_ch2_data.json` — exact rational + closed form
  data per $s \in \{1, \ldots, 5\}$.
- `debug/sprint_q5p_ch2_memo.md` — this memo.

## Recommended paper edit (already partly landed in Q5'-CH-1)

Paper 55 §subsec:open_m2_m3 received the Q5'-CH-1 sharpening. A small
extension would land the Q5'-CH-2 M2 verification as a single sentence
inside the same paragraph:

> The continuum M2 prediction is bit-exact at integer $s$: $\zeta_{D^2}(s)$
> for $s \in \{1, 2, 3, 4, 5\}$ each sits as a two-term polynomial in
> $\pi^{2k}$ with rational coefficients (memo
> `debug/sprint_q5p_ch2_memo.md`).

Sprint will apply this as a one-line addition.

## One-line verdict

CLEAN POSITIVE — bit-exact two-term M2 identification at integer
$s \in \{1, \ldots, 5\}$, supplying the continuum period-content side of
the $\omega^{\mathrm{tri}}$ symbol that Q5'-CH-1 located at the finite
cutoff; second concrete stone of the cosmic-Galois $U^*$ bridge.
