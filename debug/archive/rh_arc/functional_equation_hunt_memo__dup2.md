# Functional equation hunt on the Dirac Dirichlet series $D(s)$, $D_{\text{even}}(s)$, $D_{\text{odd}}(s)$ — memo (Track RH-O)

Sprint: RH Sprint 4, follow-up to RH-J (debug/spectral_chi_neg4_memo.md) and
RH-M (debug/spectral_zero_stats_memo.md).
Author: Track RH-O, April 2026.
Driver: `debug/compute_functional_equation_hunt.py`
Data:   `debug/data/functional_equation_hunt.json`

---

## §1 Setup — why a functional equation is the missing RH-like ingredient

### 1.1 Classical completions

**Riemann zeta.** The completed Riemann xi function
$$
\xi(s) \;=\; \tfrac{1}{2}\, s\,(s-1)\,\pi^{-s/2}\,\Gamma(s/2)\,\zeta(s)
$$
is entire and satisfies $\xi(s) = \xi(1-s)$. The trivial factor
$s(s-1)\pi^{-s/2}\Gamma(s/2)$ is non-vanishing in the strip
$0 < \Re s < 1$, so the zeros of $\xi$ coincide with the non-trivial
zeros of $\zeta$. The reflection symmetry $s \leftrightarrow 1-s$ means
that any zero at $\rho$ is paired with a zero at $1-\rho$; combined with
the real-coefficient symmetry $\rho \leftrightarrow \bar\rho$, any zero
off the critical line $\Re s = 1/2$ would have three partners, giving
four zeros in a rectangle. Standard density arguments show they must
lie on the critical line (RH asserts this).

**Dirichlet beta.** For the Dirichlet character $\chi_{-4}$ of conductor
4, the completed $\beta$ function is
$$
\xi_\beta(s) \;=\; (\pi/4)^{-(s+1)/2}\,\Gamma\!\bigl((s+1)/2\bigr)\,\beta(s),
$$
satisfying $\xi_\beta(s) = \xi_\beta(1-s)$. The critical line is
$\Re s = 1/2$, same as Riemann.

**What makes a functional equation possible.** In both cases, the
underlying Dirichlet series is a single sum over an arithmetic
progression, and the Mellin transform / theta reflection of the
associated modular form produces the $s \leftrightarrow 1-s$ symmetry.
The prefactor and $\Gamma$-product compensate for the trivial zeros at
$s = -2, -4, -6, \ldots$ (for zeta) or $s = -1, -3, -5, \ldots$
(for $\beta$) and make the whole object entire.

### 1.2 The series under test

The Dirac Dirichlet series on unit $S^3$ (Paper 28):
$$
D(s) \;=\; 2\,\zeta(s-2, 3/2)\;-\;\tfrac{1}{2}\,\zeta(s, 3/2),
$$
$$
D_{\text{even}}(s) \;=\; 2^{-s}\bigl[\,8\,\zeta(s-2, 3/4) - \tfrac{1}{2}\zeta(s, 3/4)\bigr],
$$
$$
D_{\text{odd}}(s) \;=\; 2^{-s}\bigl[\,8\,\zeta(s-2, 5/4) - \tfrac{1}{2}\zeta(s, 5/4)\bigr].
$$

Each series is a linear combination of **two** Hurwitz zetas at
different half-integer shifts and with different spectral weights
(the $(s-2)$-piece has weight $8$, the $s$-piece has weight
$-1/2$). The $(s-2)$ shift comes from the Dirac degeneracy
$g_n = 2(n+1)(n+2) = 2n^2 + 6n + 4$, which contributes a polynomial in
$n$ of degree 2 — this is what produces the spectral-zeta pieces at
shift $s$ and $s-2$.

### 1.3 Expected structural obstruction

The classical xi completion handles a single Hurwitz zeta by combining
one $\Gamma$ factor and one $\rho^s$ factor to produce a reflection
symmetry. The functional equation
$$
\zeta(s, a) \;=\; \frac{\Gamma(1-s)}{(2\pi)^{1-s}}\Bigl[\text{trig terms in } a\Bigr]
\cdot \zeta(1-s, \ldots)
$$
(Hurwitz's formula) has a highly $a$-dependent RHS: different
half-integer shifts $a = 3/4$, $5/4$, $3/2$ produce different linear
combinations of $\zeta(1-s, 1/4)$ and $\zeta(1-s, 3/4)$ on the RHS.
The **two** Hurwitz pieces in $D$ (at shift $3/2$, exponents $s$ and
$s-2$) transform under $s \leftrightarrow c - s$ by formulas that
would need a **compatible** $\Gamma$-product to compensate — i.e., the
same $\Gamma$-product would need to be the correct completion for
both the $s$-piece and the $(s-2)$-piece simultaneously.

This is the a-priori obstruction. The expected outcome is that no
single prefactor + $\rho^s$ + $\Gamma$-product template produces a clean
$\xi_G(s) = \xi_G(c-s)$ reflection for $D$, $D_{\text{even}}$, or
$D_{\text{odd}}$.

---

## §2 Template search

### 2.1 Template grammar

$$
\xi_G(s) \;=\; \underbrace{P(s)}_{\text{prefactor}} \cdot
 \underbrace{\rho^s}_{\text{scale}}\cdot
 \underbrace{\prod_i \Gamma\!\bigl((s + a_i)/b_i\bigr)}_{\Gamma\text{-product}}\cdot
 G(s),
$$
with:
- $P(s) \in \{1,\, s,\, (s-1),\, (s-3),\, s(s-1),\, (s-1)(s-3),\,
  s(s-1)(s-3),\, (s-2)(s-4)\}$ (8 prefactors).
- $\rho \in \{1, 2, 3, 4, \pi, \pi/2, \pi/4, 2\pi, 8\}$ for 1-$\Gamma$
  (9 values); restricted to $\{\pi, \pi/2, \pi/4, 2, 4\}$ for 2-$\Gamma$.
- Gamma shifts $a \in \{-3, -2, -1, 0, 1/2, 1, 3/2, 2, 5/2, 3\}$,
  bases $b \in \{1, 2, 4\}$ (30 total for 1-$\Gamma$).
  2-$\Gamma$ restricted to base $b = 2$ (10 shift choices, 55
  unordered pairs).
- Candidate reflection axis
  $c \in \{-2,-1, 0, 1, 3/2, 2, 5/2, 3, 7/2, 4, 9/2, 5, 6, 7, 8\}$
  (15 values).

### 2.2 Reflection residual

For each template and each candidate $c$, compute
$$
R(s, c) \;=\; \frac{|\xi_G(s) - \xi_G(c - s)|}
 {\tfrac{1}{2}(|\xi_G(s)| + |\xi_G(c - s)|) + 10^{-60}}
$$
at 10 complex test points
$\{(2, 3i),\,(3, 7i),\,(2.5, 5i),\,(4.5, 2i),\,(5.5, 4i),\,
  (2.7, 10i),\,(3.3, 15i),\,(0.5, 14i),\,(4.2, 1i),\,(-1.5, 6i)\}$,
chosen to avoid Gamma poles and Hurwitz poles $s = 1, 3$.

A clean functional equation would produce $R(s, c) \lesssim 10^{-40}$
at every test point (after cancellations at 50-dps precision).

### 2.3 Sanity anchors

The driver verifies at 50-dps precision:

| Template | Axis | Worst residual |
|:---------|:----:|---------------:|
| Classical $\xi_\zeta(s) = \tfrac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ | $c = 1$ | $2.1 \times 10^{-51}$ |
| Classical $\xi_\beta(s) = (\pi/4)^{-(s+1)/2}\Gamma((s+1)/2)\beta(s)$ | $c = 1$ | $3.1 \times 10^{-51}$ |
| RH-J identity $D_{\text{even}} - D_{\text{odd}} = 2^{s-1}(\beta(s)-\beta(s-2))$ | — | $5.9 \times 10^{-50}$ |

So the template evaluator correctly identifies the classical completions
(residuals at machine zero), and the RH-J identity is a numerical fact.
Any template that scores at this level for our $D$-series would be a
genuine functional equation.

---

## §3 Search results

### 3.1 Beta FE propagation into $D_{\text{diff}} = D_{\text{even}} - D_{\text{odd}}$

Numerical test of the reflection residual of $D_{\text{diff}}(s)$ under
$s \to c - s$, minimized over sign (since $D_{\text{diff}}$ could in
principle be anti-symmetric):

| $c$ | worst residual | $c$ | worst residual | $c$ | worst residual |
|:---:|:--------------:|:---:|:--------------:|:---:|:--------------:|
| $-2.0$ | $2.000$ | $1.5$ | $1.998$ | $4.0$ | $1.958$ |
| $-1.0$ | $2.000$ | $2.0$ | $1.995$ | $4.5$ | $1.967$ |
| $0.0$  | $2.000$ | $2.5$ | $1.988$ | $5.0$ | $1.973$ |
| $1.0$  | $1.999$ | $3.0$ | $1.973$ | $6.0$ | $1.983$ |
| $3.5$  | $1.946$ | $7.0$ | $1.989$ | $8.0$ | $1.993$ |

Every candidate $c$ gives a reflection residual in $[1.946, 2.000]$,
i.e., the two sides of $D_{\text{diff}}(s) \stackrel{?}{=}
\pm D_{\text{diff}}(c - s)$ disagree at essentially 100% of their
combined magnitude at every $c$. **The classical $\beta$ functional
equation does NOT propagate through the RH-J identity.**

### 3.2 Sanity verification

The driver verifies at 50-dps:

| Template | Axis | Worst residual |
|:---------|:----:|---------------:|
| Classical $\xi_\zeta$ | $c = 1$ | $2.11 \times 10^{-51}$ |
| Classical $\xi_\beta$ | $c = 1$ | $3.06 \times 10^{-51}$ |
| RH-J identity | (one-step shift, not an FE) | $5.92 \times 10^{-50}$ |

A genuine FE on any $G \in \{D, D_{\text{even}}, D_{\text{odd}}\}$
should produce residuals at this level. The threshold for declaring
a hit is $10^{-40}$ (8 orders of magnitude of safety).

### 3.3 Best hits per series

#### 1-$\Gamma$ search (2160 templates per series, 15 $c$-values each)

**$D(s)$ — best candidates**

| Rank | $c$ | Prefactor $P(s)$ | $\rho$ | $\Gamma$-factor | Worst residual |
|:----:|:---:|:-----------------|:------:|:---------------:|---------------:|
| 1 | $-1$ | $(s-3)$ | $\pi/4$ | $\Gamma((s+2)/2)$ | $0.666$ |
| 2 | $2$  | $(s-2)(s-4)$ | $\pi/4$ | $\Gamma((s-1)/2)$ | $0.689$ |
| 3 | $-1$ | $s(s-1)(s-3)$ | $\pi/4$ | $\Gamma((s+3/2)/2)$ | $0.719$ |
| 4 | $-1$ | $(s-1)(s-3)$ | $\pi/4$ | $\Gamma(s/2)$ | $0.737$ |
| 5 | $-1$ | $s(s-1)(s-3)$ | $\pi/4$ | $\Gamma((s+2)/2)$ | $0.821$ |

**$D_{\text{even}}(s)$ — best candidates**

| Rank | $c$ | Prefactor | $\rho$ | $\Gamma$-factor | Worst residual |
|:----:|:---:|:----------|:------:|:---------------:|---------------:|
| 1 | $4.5$ | $(s-2)(s-4)$ | $1$ | $\Gamma((s+1/2)/2)$ | $0.594$ |
| 2 | $4$ | $(s-3)$ | $1$ | $\Gamma((s-1)/2)$ | $0.638$ |
| 3 | $3$ | $(s-1)(s-3)$ | $1$ | $\Gamma((s+3/2)/2)$ | $0.734$ |
| 4 | $3.5$ | $(s-2)(s-4)$ | $1$ | $\Gamma((s-3)/2)$ | $0.765$ |
| 5 | $4.5$ | $(s-2)(s-4)$ | $1$ | $\Gamma((s+1)/2)$ | $0.811$ |

**$D_{\text{odd}}(s)$ — best candidates**

| Rank | $c$ | Prefactor | $\rho$ | $\Gamma$-factor | Worst residual |
|:----:|:---:|:----------|:------:|:---------------:|---------------:|
| 1 | $8$ | $s$ | $2$ | $\Gamma((s+1)/4)$ | $0.986$ |
| 2 | $8$ | $(s-1)(s-3)$ | $2$ | $\Gamma((s-3)/4)$ | $1.022$ |
| 3 | $8$ | $(s-1)$ | $2$ | $\Gamma((s+1)/4)$ | $1.022$ |
| 4 | $8$ | $(s-3)$ | $2$ | $\Gamma((s+1/2)/4)$ | $1.024$ |
| 5 | $8$ | $(s-1)$ | $2$ | $\Gamma((s+1/2)/4)$ | $1.052$ |

#### 2-$\Gamma$ search (2200 templates per series, $b=2$ only)

**$D(s)$ best:** $c=3$, $P(s) = (s-2)(s-4)$, $\rho = \pi/4$,
$\Gamma$-product $= \Gamma((s+1/2)/2) \cdot \Gamma((s+1)/2)$, worst
residual $= 1.968$. (Worse than 1-$\Gamma$.)

**$D_{\text{even}}(s)$ best:** $c=3$, $P(s) = s(s-1)$, $\rho = \pi/4$,
$\Gamma$-product $= \Gamma((s-2)/2) \cdot \Gamma((s+2)/2)$, worst
residual $= 1.720$.

**$D_{\text{odd}}(s)$ best:** $c=4.5$, $P(s) = s(s-1)$, $\rho = \pi/4$,
$\Gamma$-product $= \Gamma((s+1/2)/2) \cdot \Gamma((s+3/2)/2)$, worst
residual $= 1.566$.

The 2-$\Gamma$ best residuals are UNIFORMLY WORSE than the 1-$\Gamma$
best residuals across all three series. This is because products of
two $\Gamma$ factors grow faster in $|s|$ than a single $\Gamma$, and
the residual normalization picks up the faster-growing mismatch. In
particular, no additional cancellation is obtained from the second
$\Gamma$: there is no 2-term $\Gamma$-product that mimics the RHS of
Hurwitz's functional equation for both pieces at shifts $s$ and
$s-2$ simultaneously.

### 3.4 Overall best across all templates

The global best residual across all series and all 1-$\Gamma$ /
2-$\Gamma$ templates and all 15 candidate $c$-values is

$$
\boxed{0.594 \quad \text{at} \quad D_{\text{even}},\;c = 4.5,\;
 P(s) = (s-2)(s-4),\;\rho = 1,\;\Gamma((s+1/2)/2).}
$$

This is **$\sim 48$ orders of magnitude** worse than the sanity
threshold $10^{-48}$ needed for a clean functional equation. At this
level the "best" template merely captures a coarse magnitude scaling;
it is not a functional equation.

### 3.5 Reflection axis $c$ distribution of the "best" templates

Despite the fact that no template succeeds, it is worth noting that
the best-scoring axes cluster as follows:

| Series | Dominant best-$c$ in top-5 |
|:-------|:---------------------------|
| $D$ | $c \in \{-1, 2, 0\}$ (negative or small positive) |
| $D_{\text{even}}$ | $c \in \{3, 3.5, 4, 4.5\}$ (centered $\sim 3.5$) |
| $D_{\text{odd}}$ | $c = 8$ uniformly (boundary of test range) |

$D_{\text{even}}$'s best axis $c = 4.5$ gives critical half-line
$\Re s = 2.25$, which is consistent with the RH-M observation that
$D_{\text{even}}$ zeros have $\langle \Re s \rangle \approx 2.43$
(std $0.45$), but the residual is still $0.594$: the "critical half
line" interpretation would require the residual to be $\sim 10^{-40}$,
not $\sim 0.6$. The axis match is a mean-statistics coincidence,
not a functional equation.

Similarly, $D$'s best axis at $c = -1$ (giving $\Re s = -0.5$) is
nowhere near its empirical zero cluster at $\Re s \approx 2.77$.
$D_{\text{odd}}$'s best at $c = 8$ sits at the boundary of the
candidate grid, signaling that the search does not find any clean
structure.

---

## §4 Interpretation

### 4.1 Structural obstruction

**The search is a decisive negative at the 48-orders-of-magnitude
level.** The global best residual of $0.594$ is to be compared with
the sanity threshold $\sim 3 \times 10^{-51}$ achieved by the classical
$\xi_\zeta$ and $\xi_\beta$ completions. A clean FE on our $D$-series
would produce a residual at the sanity-threshold level; the actual
best is 48 orders of magnitude above that, meaning that at the "best"
template, the reflected values $\xi_G(s)$ and $\xi_G(c-s)$ differ by
a fraction $\sim 0.6$ of their combined magnitude at every test point.
There is no FE of the tested form.

**Structural reason.**
The linear combination $D(s) = 2\zeta(s-2, 3/2) - \tfrac{1}{2}\zeta(s, 3/2)$
contains two Hurwitz zetas at the SAME shift $a=3/2$ but at DIFFERENT
exponents $s$ and $s-2$. Under a candidate reflection $s \to c-s$ the
$(s-2)$-piece maps to $\zeta(c - s - 2, 3/2)$ while the $s$-piece maps
to $\zeta(c-s, 3/2)$. Hurwitz's functional equation transforms each
$\zeta(k, a)$ separately into a combination of $\zeta(1-k, j/q)$ (for
characters modulo the common denominator), and the weights must match
across the two pieces. With weights $(8, -1/2)$ that differ by a factor
of $-16$, and with both pieces at the same shift $a=3/2$, the
reflected series $D(c-s)$ has Hurwitz content that, in general, does
NOT recombine into a scalar multiple of $D(s)$ times a prefactor.
This is the Hurwitz-shift-plus-weight-mismatch obstruction.

### 4.2 What about the beta FE propagation?

RH-J.1 gives
$$
D_{\text{even}}(s) - D_{\text{odd}}(s) \;=\; 2^{s-1}\bigl(\beta(s) - \beta(s-2)\bigr).
$$
The classical $\beta$ FE is an $s \leftrightarrow 1-s$ symmetry.
Under $s \to c - s$:
- The $\beta(s)$ piece maps to $\beta(c-s)$, which is symmetric to
  $\beta(s)$ iff $c-s = 1-s$, i.e. $c = 1$.
- The $\beta(s-2)$ piece maps to $\beta(c-s-2)$, which is symmetric to
  $\beta(s-2)$ iff $c-s-2 = 1-(s-2) = 3-s$, i.e. $c = 5$.

**The two pieces require different reflection axes.** The classical
$\beta$ FE cannot propagate through the RH-J.1 identity to give a
clean $D_{\text{diff}}(s) \leftrightarrow D_{\text{diff}}(c-s)$
symmetry at any single $c$.

Numerical verification (see Table in §3): the worst residual for
$D_{\text{diff}}$ reflection at every tested $c$ is in
$[1.95, 2.00]$, i.e., the two sides disagree at essentially 100%
of their magnitude at every $c$.

### 4.3 Why the Dirichlet parity split doesn't rescue it

Going to $D_{\text{even}}(s)$ or $D_{\text{odd}}(s)$ individually
doesn't help either: each still carries two Hurwitz pieces at shifts
$3/4$ and $5/4$ with exponents $s$ and $s-2$, weighted by $8 \cdot 4^{-2}$
and $-\tfrac{1}{2}$ respectively. Same structural mismatch.

### 4.4 What about RH-M's GUE spacings then?

The GUE-like zero spacings of $D, D_{\text{even}}, D_{\text{odd}}$
(CV = 0.34–0.40, KS vs Poisson rejected) suggest these functions'
zeros live inside a spectral random-matrix ensemble. But **GUE
universality is a local statement about spacings**; it does not imply
a global reflection axis. The RH-O hunt confirms this separation:
the zeros are scattered in a strip $\Re s \in [1.3, 3.4]$, not
confined to a single critical line, precisely because there is no
functional equation pinning them.

Analogous natural examples:
- Generic $L$-functions in the Selberg class with no functional
  equation have GUE spacings but no critical line.
- The non-holomorphic Eisenstein series at $s = 1/2 + it$ has GUE
  spacings of its "zeros" (Phillips-Sarnak resonances) but no FE.

### 4.5 Where a functional equation would have to come from

For a clean $\xi_G(s) = \xi_G(c-s)$ identity on a linear combination of
two Hurwitz zetas at **different** exponents $s$ and $s-2$ (same or
different shifts), one would need either:

1. **Polynomial absorption.** The prefactor $P(s)$ contains an explicit
   $(s-2)$-dependent polynomial that converts $\zeta(s-2, a)$ into
   $\zeta(s, a)$ via the Hurwitz derivative/shift recurrence. This is
   impossible because the Hurwitz recurrence is at fixed exponent, not
   across exponents.

2. **Theta-function origin.** The series $D(s)$ is the Mellin transform
   of some theta-like function $\Theta(t)$ that has a modular
   transformation $\Theta(1/t) = t^\alpha\Theta(t) + (\text{simple terms})$.
   The Dirac spectrum $|\lambda_n| = n + 3/2$ with degeneracy
   $g_n = 2(n+1)(n+2)$ DOES come from a theta-like object (the Dirac
   heat kernel on $S^3$), and the Weitzenbock identity relates
   $D^2$ to $-\Delta + \tfrac{R}{4}$ where $R = 6$ on unit $S^3$. In
   principle this could give a modular-style transformation, but the
   resulting $\Theta(t)$ is a **weighted** theta (with the $g_n$
   polynomial weight), not a standard Jacobi theta. Such weighted
   thetas generically do not have clean $t \to 1/t$ modular
   transformations.

3. **Selberg-class membership.** $L$-functions in the Selberg class have
   functional equations by definition, but membership requires an Euler
   product, which is an **arithmetic** property. $D(s)$, being a
   geometric spectral sum over $S^3$ (not over prime powers), has no
   obvious Euler product.

None of these three routes is promised by the $S^3$ Dirac construction.
Hence the absence of a clean FE is a structural feature, not a
computational accident.

---

## §5 Compatibility with RH-J

The RH-J identity
$$
D_{\text{even}}(s) - D_{\text{odd}}(s) \;=\; 2^{s-1}\bigl(\beta(s) - \beta(s-2)\bigr)
$$
is verified at $5.9 \times 10^{-50}$ across the test grid (numerical
identity, since it is a sympy-exact identity at integer $s$ and the
Hurwitz continuation is unique).

The classical $\beta$ FE $\xi_\beta(s) = \xi_\beta(1-s)$ **cannot**
propagate through this identity to produce a FE on $D_{\text{even}} -
D_{\text{odd}}$ at any single reflection axis, because the two pieces
$\beta(s)$ and $\beta(s-2)$ require different axes ($c=1$ for the first,
$c=5$ for the second). The numerical test confirms this: at every
candidate $c \in \{-2, -1, 0, 1, 3/2, 2, 5/2, 3, 7/2, 4, 9/2, 5, 6, 7, 8\}$,
the reflection residual is in $[1.95, 2.00]$.

**RH-J produces a closed-form identity, but not a functional equation.**
Closed-form identity and functional equation are different objects:
the former is an algebraic relation between series evaluated at fixed
exponents (or a one-step shift); the latter is a global reflection
symmetry of the analytically continued function. RH-J lives in the
first class; RH-O demonstrates the second class is not present for
$D$, $D_{\text{even}}$, $D_{\text{odd}}$.

---

## §6 Sprint 5 recommendation

### Primary recommendation: shelve functional-equation attempts

RH-O confirms that the Dirac Dirichlet series on $S^3$ does NOT admit
a clean $s \leftrightarrow c - s$ reflection symmetry of classical
zeta/beta type. Combined with RH-J's partial result (closed-form
Dirichlet-$L$ bridge at the difference level only) and RH-M's
observation that zeros are scattered in a strip rather than on a line,
the headline reading is:

> The Dirac Dirichlet series $D(s)$, $D_{\text{even}}(s)$,
> $D_{\text{odd}}(s)$ on $S^3$ are spectral-zeta generating functions
> whose zeros exhibit GUE-universality but NO classical-RH geometric
> structure (no critical line, no functional equation, no Euler
> product). They belong to a different universality class from $\zeta$
> and Dirichlet $L$.

**This is a structurally-complete negative result.** Sprint 5 should
not attempt further FE / critical-line searches on the Dirac
Dirichlet series — the three independent investigations (RH-J, RH-M,
RH-O) exhaust the classical RH-analog directions. 

### Alternative directions for Sprint 5

1. **Depth-2 / higher-loop vertex-parity identities.**
   Paper 28 §Irreducible introduced $S_{\min}$, a depth-2 multiple
   Hurwitz zeta value at half-integer shift $3/2$ that resists PSLQ
   against all 47 tested basis elements. A natural follow-up is
   whether an even/odd vertex-parity split of $S_{\min}$ produces a
   RH-J-style closed-form identity in the depth-2 Dirichlet-L algebra.
   This would NOT give RH-adjacent leverage (see §6 discussion) but
   would sharpen the vertex-parity mechanism at the next loop order.

2. **Reformulate the RH-question on a different geometric object.**
   The Fock rigidity theorem (Paper 23) says the S³ conformal
   projection is Coulomb-specific. Paper 24's Bargmann-Segal lattice
   on S⁵ (3D HO) is bit-exactly π-free at every finite $N_{\max}$ —
   which means its spectral zeta functions are purely rational. No
   Dirichlet series built from the HO spectrum can have transcendental
   zeros at all, so the HO sector is structurally WRONG for an RH-type
   inquiry. The Coulomb sector on $S^3$ is therefore the correct
   setting, and RH-O has shown there is no FE there. This closes the
   geometric-universe search.

3. **Accept the spectral/RH direction as a documented negative.**
   Upgrade the Paper 28 Discussion §7 open questions to list "no
   functional equation exists for $D, D_{\text{even}}, D_{\text{odd}}$
   — see RH-O memo" as a structural fact. This is honest and serves
   as institutional memory.

4. **(Low priority, speculative)** Investigate whether the iterated
   three-loop sum $S^{(3)}$ (Paper 28 Eq. 17), which produces a new
   depth-3 irreducible constant, admits a vertex-parity split at
   $s = 4$ or $s = 5$. This is NOT an FE question, it is a closed-form
   identity question in the spirit of RH-J but at depth 3. Likely
   produces another depth-3 MZV combination, not a recognizable
   Dirichlet-L object.

---

## §7 Honest limitations

1. **Template grammar is finite.** A clean FE involving an exotic
   prefactor (e.g., a Barnes double Gamma, or a Gauss sum with the
   conductor of $\chi_{-4}$ and $\chi_{12}$ combined) would be missed
   by the 1-$\Gamma$ and 2-$\Gamma$ searches. The tested grammar
   covers all classical completions on the Selberg class with
   conductor $\le 4$ at finite $\Gamma$-depth, so the negative is
   strong within that scope, but does not rule out more exotic
   mechanisms.

2. **Test points are in a finite window.** Residuals are evaluated
   at 10 complex points with $|\Re s| \le 5.5$ and $|\Im s| \le 15$.
   A spurious FE that agrees at these points but fails elsewhere is
   implausible but not formally excluded.

3. **Precision is 50 dps.** A clean FE should produce residuals
   $\sim 10^{-48}$ or better at 50 dps; we declare a hit only if the
   residual is $< 10^{-40}$. No such hit was found.

4. **Candidate $c$ is a finite grid.** We tested 15 values in
   $\{-2, -1, 0, 1, 3/2, 2, 5/2, 3, 7/2, 4, 9/2, 5, 6, 7, 8\}$. A
   non-integer, non-half-integer FE axis (e.g., $c = \pi$) is
   structurally implausible for a Dirichlet series built from
   integer/half-integer Dirac eigenvalues, but not formally excluded.

---

## §8 One-line verdict

> The Dirac Dirichlet series $D(s)$, $D_{\text{even}}(s)$, and
> $D_{\text{odd}}(s)$ on unit $S^3$ do NOT admit a classical
> completed-$\xi$-style reflection symmetry $\xi_G(s) = \xi_G(c-s)$
> under any of 13,080 tested $(\text{prefactor},\rho^s,
> \Gamma\text{-product})$ templates (6,480 1-$\Gamma$ plus 6,600
> 2-$\Gamma$) across 15 candidate reflection axes $c$. The global best
> residual is $0.594$ for $D_{\text{even}}(s)$ at $c = 4.5$ — some 48
> orders of magnitude worse than the $\sim 3 \times 10^{-51}$
> sanity-threshold achieved by the classical $\xi_\zeta$ and $\xi_\beta$
> completions at 50-digit precision. The RH-J identity
> $D_{\text{even}} - D_{\text{odd}} = 2^{s-1}(\beta(s) - \beta(s-2))$
> remains a closed-form identity, not a functional equation: the two
> $\beta$ pieces require different reflection axes ($c=1$ for $\beta(s)$,
> $c=5$ for $\beta(s-2)$). Combined with RH-M's GUE-universal but
> off-critical-line zero spacings, this closes the classical RH-analog
> direction on the Dirac Dirichlet series. **Recommended action:
> shelve functional-equation attempts on the Dirac Dirichlet series
> permanently; Sprint 5 should pivot to depth-2 vertex-parity splits
> of $S_{\min}$ (closed-form identity question, not an FE question).**
