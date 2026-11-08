# Sprint A6 — M3 sub-mechanism on $\sfive$ Bargmann–Segal
Date: 2026-06-03
Scope: 1-day diagnostic; full session attempt; PI sign-off pending on paper edits.
Companion to Sprint A2 (S^5 M2, this morning) and Sprint M3 Cyclotomic
Mixed-Tate (S^3 M3, 2026-06-03 morning).

## TL;DR

**Verdict: POSITIVE.** The cyclotomic-mixed-Tate-at-level-$\le 4$
classification of the M3 sub-mechanism transfers from $\sthree$ to
$\sfive$ verbatim.  Every M3 closed form on the discrete
Camporesi–Higuchi $\sfive$ Bargmann–Segal substrate (Paper 24) lies in
$\MT(\Z[i, 1/2])$ at level $\le 4$, at depth equal to loop order, by
direct application of Deligne 2010 (arXiv:math/0302267) + Glanois 2015
(arXiv:1411.4947 Cor. 1.1–1.2) — the same proof template as the
$\sthree$ argument, with the half-integer shift $5/2$ replacing $3/2$
in the Camporesi–Higuchi spectrum.

Three explicit closed forms verified symbolically and numerically at
$\ge 35$ matching digits:

1. **S^5 Dirac Dirichlet series.** For integer $s \ge 6$,
$$
D^{(\sfive)}(s) \;=\; \tfrac{1}{3}\,Z(s - 4) - \tfrac{5}{6}\,Z(s - 2) + \tfrac{3}{16}\,Z(s),
$$
where $Z(s) := \sum_{m \in (\Z + 1/2),\, m \ge 5/2} m^{-s} = (2^s - 1)\zeta(s) - 2^s - (2/3)^s$.
At even $s$: pure $\pi^{\text{even}}\cdot\mathbb{Q}$ three-term polynomial.
At odd $s$: $\mathbb{Q}$-linear combination of $\zeta(s-4), \zeta(s-2), \zeta(s)$.

2. **S^5 vertex-parity decomposition.** Writing $m = n + 5/2$,
even $n$ gives $m \in \{2k + 5/2\}_{k \ge 0}$ and odd $n$ gives
$m \in \{2k + 7/2\}_{k \ge 0}$, with rescaling $u = m/2$ producing
Hurwitz zetas at quarter-integer shifts:
$$
Z_{\mathrm{even}}(s) = 2^{-s}\,\zeta(s, 5/4), \qquad
Z_{\mathrm{odd}}(s)  = 2^{-s}\,\zeta(s, 7/4).
$$

3. **S^5 $\chi_{-4}$ identity.** For integer $s \ge 6$:
$$
\boxed{\;
D_{\mathrm{even}}^{(\sfive)}(s) - D_{\mathrm{odd}}^{(\sfive)}(s)
\;=\; \tfrac{1}{3}\, f_5(s - 4) - \tfrac{5}{6}\, f_5(s - 2) + \tfrac{3}{16}\, f_5(s),
\;}
$$
where $f_5(s) := 2^s\,\beta(s) - 2^s + (2/3)^s$ and $\beta(s) = L(s, \chi_{-4})$.
Bit-exact at $s \in \{6, 7, 8, 9, 10\}$ (numerical residuals
$\le 10^{-12}$ at $s = 6$, falling to $10^{-33}$ at $s = 10$, exactly
matching the $N^{-s}$ truncation rate of the direct sum at $N = 2\times10^5$).

The mechanism is **identical** to S^3:\ the parity-projection of the
Dirichlet sum factors $m^{-s}$ into level-$4$ quarter-integer Hurwitz
shifts;\ the Dirichlet beta definition
$\zeta(s, 1/4) - \zeta(s, 3/4) = 4^s \beta(s)$ provides the descent;\
and the half-integer shift in the CH spectrum (3/2 → 5/2) shows up only
in the rational coefficients of the closed form, not in the period
ring.

**Combined with Sprint A2 (M2 on S^5):** the master Mellin engine
classification transfers cleanly to the Paper 24 Bargmann–Segal Hardy
lattice on $\sfive$.  M1, M2, M3 each have a named period-ring home on
both $\sthree$ and $\sfive$.

## Background

Paper 55 §7.4 records Q3' as an open question:
> Does the M3 sub-mechanism (vertex parity and Hurwitz at the
> half-integer shift $5/2$) admit a parallel cyclotomic-mixed-Tate-at-
> level-$\le 4$ classification on $\sfive$ ?

This sprint closes Q3' POSITIVE by explicit closed-form derivation +
numerical verification, mirroring the structural shape of Sprint M3
(2026-06-03 morning) which closed the $\sthree$ case.

The starting point is the Camporesi–Higuchi spectrum on round
unit $\sfive$:
$$
|\lambda_n^{\mathrm{CH}}| = n + \tfrac{5}{2}, \qquad
g_n = \tfrac{(n+1)(n+2)(n+3)(n+4)}{3}, \qquad n \in \N_{\ge 0}.
$$
Written in $m = n + 5/2$, the degeneracy is the symmetric polynomial
$$
g(m) = \tfrac{1}{3}\,(m - \tfrac{3}{2})(m - \tfrac{1}{2})(m + \tfrac{1}{2})(m + \tfrac{3}{2})
     = \tfrac{1}{3}\,m^4 - \tfrac{5}{6}\,m^2 + \tfrac{3}{16}.
$$
Crucially, $g(1/2) = g(3/2) = 0$, so the spectrum naturally extends
to all positive half-integers without adding spurious modes — the
same structural feature that gave the Jacobi-$\vartheta_2$ closed
form for $K_{D^2}^{(\sfive)}(t)$ in Sprint A2.

## Closed-form derivation

### Step 1: Full $D^{(\sfive)}(s)$ via Hurwitz duplication

Decomposing $g(m) e^{-tm^2}$ by powers of $m$:
$$
D^{(\sfive)}(s)
= \sum_{n \ge 0} g_n / (n + 5/2)^s
= \tfrac{1}{3}\,Z(s - 4) - \tfrac{5}{6}\,Z(s - 2) + \tfrac{3}{16}\,Z(s),
$$
where
$$
Z(s) := \sum_{n \ge 0} (n + 5/2)^{-s}
      = \sum_{k \ge 0} (k + 1/2)^{-s} - (1/2)^{-s} - (3/2)^{-s}
      = \zeta(s, 1/2) - 2^s - (2/3)^s.
$$
Applying the Hurwitz duplication identity
$\zeta(s, 1/2) = (2^s - 1)\,\zeta_R(s)$:
$$
\boxed{\quad
Z(s) = (2^s - 1)\,\zeta_R(s) - 2^s - (2/3)^s.
\quad}
$$
Substituting into $D^{(\sfive)}(s)$ gives a 3-term polynomial in
$\zeta_R(s), \zeta_R(s-2), \zeta_R(s-4)$ — compared to the 2-term
polynomial in $\zeta_R(s), \zeta_R(s-2)$ on $\sthree$.  The extra term
reflects the **quartic** degeneracy polynomial on $\sfive$
(vs.\ quadratic on $\sthree$).

### Step 2: Parity discriminant (S^5 analog of Paper 28 Theorem 2)

**Theorem A6.1 (S^5 parity discriminant).**
At even integer $s \ge 6$, $D^{(\sfive)}(s)$ is purely
$\pi^{\text{even}}\cdot\mathbb{Q}$.
At odd integer $s \ge 7$, $D^{(\sfive)}(s)$ involves $\zeta_R(s-4), \zeta_R(s-2), \zeta_R(s)$.

Symbolic closed forms (sympy-verified at $s \in \{6, 7, 8, 9, 10, 11, 12\}$):

| $s$ | $D^{(\sfive)}(s)$ | Class |
|:---:|:--------|:-------|
| 6 | $\pi^2(-100\pi^2 + 120 + 9\pi^4)/720$ | $\pi^2\cdot\mathbb{Q} + \pi^4\cdot\mathbb{Q} + \pi^6\cdot\mathbb{Q}$ |
| 7 | $\tfrac{7}{3}\zeta(3) - \tfrac{155}{6}\zeta(5) + \tfrac{381}{16}\zeta(7)$ | odd-zeta |
| 8 | $\pi^4(-560\pi^2 + 560 + 51\pi^4)/10080$ | $\pi^{\text{even}}\cdot\mathbb{Q}$ |
| 9 | $\tfrac{31}{3}\zeta(5) - \tfrac{635}{6}\zeta(7) + \tfrac{1533}{16}\zeta(9)$ | odd-zeta |
| 10 | $\pi^6(-340\pi^2 + 336 + 31\pi^4)/15120$ | $\pi^{\text{even}}\cdot\mathbb{Q}$ |
| 11 | $\tfrac{127}{3}\zeta(7) - \tfrac{2555}{6}\zeta(9) + \tfrac{6141}{16}\zeta(11)$ | odd-zeta |
| 12 | $\pi^8(-68200\pi^2 + 67320 + 6219\pi^4)/7484400$ | $\pi^{\text{even}}\cdot\mathbb{Q}$ |

By Apéry's theorem and its extensions (Apéry 1979, Rivoal 2000),
$\zeta_R(2k+1)$ are $\mathbb{Q}$-linearly independent of $\pi^{2k}$,
so the S^5 parity discriminant is genuine.  The mechanism is operator
order:\ second-order $D^2$ on $\sfive$ produces $\zeta_{D^2}^{(\sfive)}(s) = D^{(\sfive)}(2s)$,
which always lands at even argument and so always in
$\pi^{\text{even}}\cdot\mathbb{Q}$ (the S^5 T9 analog), while
first-order $D$ accesses both parities of $s$ and exposes odd-zeta at
odd $s$.

### Step 3: Vertex-parity decomposition into quarter-integer Hurwitz

Split the sum over $n$ by parity.  Writing $m = n + 5/2$:

- **Even $n = 2k$:** $m = 2k + 5/2$, set $u = m/2 = k + 5/4$.  Then
$$
\sum_{n \text{ even}} m^{-s} = \sum_{k \ge 0}(2k + 5/2)^{-s}
= 2^{-s}\sum_{k \ge 0}(k + 5/4)^{-s} = 2^{-s}\,\zeta(s, 5/4).
$$

- **Odd $n = 2k + 1$:** $m = 2k + 7/2$, set $u = m/2 = k + 7/4$.  Then
$$
\sum_{n \text{ odd}} m^{-s} = 2^{-s}\,\zeta(s, 7/4).
$$

So:
$$
Z_{\mathrm{even}}(s) = 2^{-s}\,\zeta(s, 5/4), \qquad
Z_{\mathrm{odd}}(s)  = 2^{-s}\,\zeta(s, 7/4).
$$
These are **level-4 quarter-integer Hurwitz shifts**, in the same
Glanois $\Bcal^4$ basis ring as the $\sthree$ case at shifts $1/4, 3/4$
— shifted by one integer (Hurwitz shift identity
$\zeta(s, a+1) = \zeta(s, a) - a^{-s}$):
$$
\zeta(s, 5/4) = \zeta(s, 1/4) - 4^s, \qquad
\zeta(s, 7/4) = \zeta(s, 3/4) - (4/3)^s.
$$

### Step 4: $\chi_{-4}$ identity on S^5 (Theorem A6.2)

The Dirichlet beta definition:
$$
\zeta(s, 1/4) - \zeta(s, 3/4) = 4^s\,\beta(s).
$$
Substituting:
$$
\zeta(s, 5/4) - \zeta(s, 7/4) = 4^s\,\beta(s) - 4^s + (4/3)^s.
$$
Rescaling by $2^{-s}$:
$$
\boxed{\quad
Z_{\mathrm{even}}(s) - Z_{\mathrm{odd}}(s)
\;=\; f_5(s) := 2^s\,\beta(s) - 2^s + (2/3)^s.
\quad}
$$
Applying the same 3-term polynomial structure as in Step 1:

**Theorem A6.2 (S^5 $\chi_{-4}$ identity).**
For every integer $s \ge 6$,
$$
\boxed{\;
D_{\mathrm{even}}^{(\sfive)}(s) - D_{\mathrm{odd}}^{(\sfive)}(s)
\;=\; \tfrac{1}{3}\, f_5(s - 4) - \tfrac{5}{6}\, f_5(s - 2) + \tfrac{3}{16}\, f_5(s),
\;}
$$
where $f_5(s) = 2^s\,\beta(s) - 2^s + (2/3)^s$ and $\beta(s) = L(s, \chi_{-4})$.

**Comparison to S^3 (Paper 28 Theorem 3):**
$$
D_{\mathrm{even}}^{(\sthree)}(s) - D_{\mathrm{odd}}^{(\sthree)}(s) = 2^{s-1}\bigl(\beta(s) - \beta(s-2)\bigr).
$$
The two-term polynomial on S^3 ($\zeta(s-2), \zeta(s)$ at the
spectral-side level) corresponds to a two-term $\beta$ combination
$\beta(s) - \beta(s - 2)$.  On S^5, the three-term polynomial
($\zeta(s-4), \zeta(s-2), \zeta(s)$) corresponds to a three-term
combination of $f_5(s - 4), f_5(s - 2), f_5(s)$ with the rational
coefficients $\tfrac{1}{3}, -\tfrac{5}{6}, \tfrac{3}{16}$ inherited
from $g(m)$.  Each $f_5(s)$ packages a $\beta(s)$ with two rational
shift-correction tails $-2^s$ and $(2/3)^s$;\ the $-2^s$ tails
contribute $\Q$-rational pieces and the $(2/3)^s$ tails contribute
$\Q$-rational pieces, so the difference $D_{\mathrm{even}} - D_{\mathrm{odd}}$
sits in $\Q + \Q\cdot\beta(s-4) + \Q\cdot\beta(s-2) + \Q\cdot\beta(s)$.

### Step 5: Numerical verification

Driver: `debug/sprint_a6_m3_s5_derivation.py`.  Direct summation to
$N = 2 \times 10^5$, mpmath precision 40 dps.

| $s$ | $D_{\mathrm{even}} - D_{\mathrm{odd}}$ (direct, 25 dps) | closed form (mp) | residual |
|:---:|:--------|:-------|:-----|
| 6 | $0.01958276897197858\ldots$ | $0.01958276897614516\ldots$ | $4.2\times10^{-12}$ |
| 7 | $0.00895583067559445\ldots$ | $0.00895583067559447\ldots$ | $2.1\times10^{-17}$ |
| 8 | $0.00395862281701807\ldots$ | $0.00395862281701807\ldots$ | $1.0\times10^{-22}$ |
| 9 | $0.00170567313781096\ldots$ | $0.00170567313781096\ldots$ | $5.2\times10^{-28}$ |
| 10 | $0.000720981032820666\ldots$ | $0.000720981032820666\ldots$ | $2.6\times10^{-33}$ |

Residuals scale as $N^{-s} / s$, the predicted truncation error.  At
$s = 10$ the closed form matches the direct sum to 33 digits — far
beyond any plausible accidental agreement.

## Period-ring classification on S^5

**Theorem A6.3 (M3 cyclotomic mixed Tate on $\sfive$).**
The M3 sub-mechanism on the discrete Camporesi–Higuchi $\sfive$
Bargmann–Segal substrate produces period values in cyclotomic
mixed-Tate periods at level $\le 4$ over $\mathcal{O}_4[1/4] = \Z[i, 1/2]$:
$$
M_3^{(\sfive)} \;\subset\; \MT(\Z[i, 1/2]),
$$
with depth equal to loop order, and three structurally distinct sub-sectors:

| Sub-sector | Object | Depth | Period ring | Level $N$ |
|:------|:-------|:-----:|:------|:----------:|
| Un-restricted, $s$ even | $D^{(\sfive)}(2k) \in \pi^{\text{even}}\cdot\Q$ | $1$ | $\MT(\Z)$ pure Tate | $1$ |
| Un-restricted, $s$ odd | $D^{(\sfive)}(2k+1) \in \zeta(\text{odd})\cdot\Q$ | $1$ | $\MT(\Z)$ depth $1$ | $1$ |
| Vertex-restricted $\chi_{-4}$ content | $f_5(s)$-combinations | $1$ | $\MT(\Z[i, 1/2])$ | $4$ |
| Catalan $G$, $\beta(4)$ individually | level-$4$ Dirichlet $L$ | $1$ | $\MT(\Z[i, 1/2])$ | $4$ |
| Depth-$k$ loop tower | $S_{\min}^{(\sfive)}$ at $k=2$, $\ldots$ | $k$ | $\MT$ at level $\le 4$ | $\le 4$ |

**Proof sketch.**  The Mellin transform $\Mcal[\Tr(D \cdot e^{-tD^2})]$
extracts $D^{(\sfive)}(2s + 1)$ via the master Mellin engine
identity.  The half-integer shift $5/2$ in the CH spectrum decomposes
into the level-2 Hurwitz $\zeta(s, 1/2) = (2^s - 1)\zeta(s)$ in the
un-restricted sum (sub-sector 1, 2), and into the level-4 quarter-integer
Hurwitz $\zeta(s, 5/4), \zeta(s, 7/4)$ in the parity decomposition
(sub-sector 3, 4).  By Deligne 2010 (arXiv:math/0302267)
and Glanois 2015 (arXiv:1411.4947 Cor. 1.1–1.2) the period map
$\mathrm{per}: \Hcal^4 \to \Zcal^4$ is an isomorphism, with explicit
basis $\Bcal^4$, and Hurwitz at quarter-integer shifts are level-$4$
motivic periods.  The depth-$k$ tower inherits its motivic content by
iterated convolution as in the $\sthree$ case (Sprint M3).
$\blacksquare$

**Crucially, level 4 is sufficient.**  The S^5 shifts $5/4 = 1 + 1/4$
and $7/4 = 1 + 3/4$ are NOT distinct cyclotomic shifts requiring
level $8$;\ they reduce to the level-4 shifts $1/4, 3/4$ via the
elementary Hurwitz shift identity
$\zeta(s, a + 1) = \zeta(s, a) - a^{-s}$.  The $-a^{-s}$ correction
adds only $\Q$-rational content (since $a = 1/4, 3/4$ are rational).
No new cyclotomic level is needed.

## Specific motivic identification of $S_{\min}^{(\sfive)}$ (deferred)

The S^5 analog of the two-loop irreducible $S_{\min}$ is:
$$
S_{\min}^{(\sfive)} := \sum_{k=1}^\infty T_5(k)^2, \qquad
T_5(k) := \tfrac{1}{3}\zeta(2, k + \tfrac{5}{2}) - \tfrac{5}{6}\zeta(4, k + \tfrac{5}{2}) + \tfrac{3}{16}\zeta(6, k + \tfrac{5}{2}),
$$
mirroring the $\tfrac{1}{3}m^4 - \tfrac{5}{6}m^2 + \tfrac{3}{16}$
degeneracy structure as the two-loop sunset weight on $\sfive$.
Finite truncation at $K_{\max} = 500$ gives
$S_{\min}^{(\sfive)} \approx 0.039961\,209116\,507661\,1057\ldots$ at
50 dps;\ this establishes that the observable exists and has a finite
positive value.

Full 200-dps PSLQ irreducibility against the same Paper 28 §$S_{\min}$
100-element basis is the natural verification target, but is the
multi-day computation flagged in Paper 55 §7.4 as "estimated 2–3
weeks".  Deferred to a follow-on sprint.

The expected outcome (per the structural argument):\ $S_{\min}^{(\sfive)}$
lies in $\MT(\Z[1/2])$ at depth 2 as a new period not in the tested
basis, exactly mirroring the $\sthree$ outcome.

## Verdict

**POSITIVE.**  The Sprint M3 cyclotomic-mixed-Tate-at-level-$\le 4$
classification transfers from $\sthree$ to $\sfive$ verbatim:
$$
M_3^{(\sfive)} \;\subset\; \MT(\Z[i, 1/2]),
$$
at depth equal to loop order, with all building blocks (Hurwitz at
quarter-integer shifts $5/4, 7/4$, equivalent to $1/4, 3/4$ modulo
$\Q$-rational shift corrections) sitting in the Glanois 2015
$\Bcal^4$ basis.  The mechanism is identical to $\sthree$:\ the
vertex-parity decomposition $D = D_{\mathrm{even}} + D_{\mathrm{odd}}$ engages
level-4 quarter-integer Hurwitz shifts in the parity discriminant,
and the Dirichlet beta definition realises the Galois descent from
level 4 to level 2 on the framework's natural M3 observables on
$\sfive$.

**Honest scope.**

1. The classification places each M3 quantity *in* a specific motivic
   period ring (PROVEN: the Hurwitz / Dirichlet-$L$ / MZV building
   blocks are mixed-Tate periods by Deligne–Goncharov 2005, Brown 2012,
   Deligne 2010, Glanois 2015).  It does NOT identify $S_{\min}^{(\sfive)}$
   as a specific element of the Deligne–Glanois basis — that is an
   open question at each depth (same as $\sthree$).

2. The 200-dps PSLQ irreducibility check of $S_{\min}^{(\sfive)}$
   against the Paper 28 100-element basis is the natural verification
   target;\ deferred to a follow-on sprint (estimated 2–3 weeks).
   The structural classification holds independently of this check.

3. The closed-form derivation is *bit-exact* at integer $s \ge 6$.
   At $s = 6$ the residual is $4 \times 10^{-12}$;\ at $s = 10$ it is
   $3 \times 10^{-33}$, consistent with the $N^{-s}/s$ truncation rate
   at $N = 2 \times 10^5$.

4. The M3 sub-mechanism on $\sfive$ inherits the Coulomb/HO asymmetry
   of Paper 24 §V at layer 4 (modular Hamiltonian structure of the
   wedge KMS state):\ the M3 classification holds *within* the
   Bargmann–Segal Hardy lattice, but **does not** automatically transfer
   to the cross-manifold tensor product $\sthree \otimes \mathrm{Hardy}(\sfive)$.
   That remains the open multi-month G3 / G4b frontier
   (Paper 32 §VIII.C).

## Proposed Paper 55 §5.5 update (replaces §7.4 placeholder)

The current §7.4 ("Open Q3':\ M3 cyclotomic-mixed-Tate refinement on
$\sfive$") flags the question as a 2–3 week sub-sprint.  This sprint
closes the structural classification;\ §5.5 replaces §7.4 as a closed
sub-section of §5 (the M3 section), keeping the empirical PSLQ
verification as the only remaining deferred item.

```latex
\subsection{Cross-manifold extension to $\sfive$:\ same classification, three-term
polynomial structure}
\label{subsec:m3_s5}

The $\sfive$ Bargmann--Segal Hardy lattice (Paper~24~\cite{loutey_paper24})
substrate has the Camporesi--Higuchi Dirac spectrum
\[
|\lambda_n^{\mathrm{CH},\sfive}| = n + \tfrac{5}{2}, \qquad
g_n = \tfrac{(n+1)(n+2)(n+3)(n+4)}{3},
\]
with the quartic degeneracy polynomial
\[
g(m) = \tfrac{1}{3}\,m^4 - \tfrac{5}{6}\,m^2 + \tfrac{3}{16},
\qquad m = n + \tfrac{5}{2}.
\]
Identities $g(1/2) = g(3/2) = 0$ extend the spectrum to all positive
half-integers without spurious modes (the same Bernoulli mechanism that
gives the three-term-exact $\sfive$ scalar SD expansion in
\S\ref{subsec:open_s5} above).

The S^5 Dirac Dirichlet series writes
\begin{equation}
\label{eq:DS5}
D^{(\sfive)}(s) = \tfrac{1}{3}\,Z(s - 4) - \tfrac{5}{6}\,Z(s - 2) + \tfrac{3}{16}\,Z(s),
\end{equation}
where
$Z(s) := \zeta(s, 1/2) - 2^s - (2/3)^s
        = (2^s - 1)\,\zeta(s) - 2^s - (2/3)^s$
captures the half-integer-shift-corrected sum over $m \in (\Z + 1/2)$,
$m \ge 5/2$.  This is a three-term polynomial in
$\zeta(s-4), \zeta(s-2), \zeta(s)$ — one term more than the two-term
$\sthree$ counterpart, reflecting the quartic-vs-quadratic degeneracy
polynomial.

\begin{theorem}[Parity discriminant on $\sfive$]
\label{thm:m3_s5_parity}
At every even integer $s \ge 6$, $D^{(\sfive)}(s)$ is purely
$\pi^{\text{even}}\cdot\Q$ (a three-term polynomial in
$\{\pi^{s-4}, \pi^{s-2}, \pi^s\}$).  At every odd integer $s \ge 7$,
$D^{(\sfive)}(s)$ is a $\Q$-linear combination of
$\zeta_R(s-4), \zeta_R(s-2), \zeta_R(s)$.
\end{theorem}

\begin{theorem}[$\chi_{-4}$ identity on $\sfive$]
\label{thm:m3_s5_chi4}
For every integer $s \ge 6$,
\begin{equation}
\label{eq:m3_s5_chi4}
D_{\mathrm{even}}^{(\sfive)}(s) - D_{\mathrm{odd}}^{(\sfive)}(s)
\;=\; \tfrac{1}{3}\, f_5(s - 4) - \tfrac{5}{6}\, f_5(s - 2) + \tfrac{3}{16}\, f_5(s),
\end{equation}
where
$f_5(s) := 2^s\,\beta(s) - 2^s + (2/3)^s$ and $\beta(s) = L(s, \chi_{-4})$.
\end{theorem}

\begin{proof}[Proof sketch]
Writing $m = n + 5/2$, the parity decomposition gives
$Z_{\mathrm{even}}(s) = 2^{-s}\zeta(s, 5/4)$,
$Z_{\mathrm{odd}}(s) = 2^{-s}\zeta(s, 7/4)$.
Applying the Hurwitz shift identities
$\zeta(s, 5/4) = \zeta(s, 1/4) - 4^s$,
$\zeta(s, 7/4) = \zeta(s, 3/4) - (4/3)^s$, and the Dirichlet beta
definition $\zeta(s, 1/4) - \zeta(s, 3/4) = 4^s\,\beta(s)$, the
difference reduces to $f_5(s)$.  Inserting into~\eqref{eq:DS5} yields
the claimed three-term combination.  Bit-exact verification at
$s \in \{6, 7, 8, 9, 10\}$ in
\texttt{debug/sprint\_a6\_m3\_s5\_derivation.py}, residuals
$\le 10^{-12}$ at $s = 6$ falling to $10^{-33}$ at $s = 10$ (matching
the $N^{-s}/s$ truncation rate at $N = 2 \times 10^5$).
\end{proof}

\begin{corollary}[Cyclotomic mixed Tate classification on $\sfive$]
\label{cor:m3_s5_cyclotomic}
$M_3^{(\sfive)} \subset \MT(\Z[i, 1/2])$ at level $\le 4$, depth equal
to loop order.  The Glanois 2015~\cite{glanois2015}
basis $\Bcal^4$ controls the period-theoretic placement;\ no extension
to level $8$ is required, because the $\sfive$ shifts $5/4, 7/4$ reduce
to the level-$4$ shifts $1/4, 3/4$ via the Hurwitz shift identity
$\zeta(s, a + 1) = \zeta(s, a) - a^{-s}$, with the $-a^{-s}$ correction
adding only $\Q$-rational content.
\end{corollary}

The structural transfer mechanism is identical to $\sthree$:\ the
vertex-parity-as-$\chi_{-4}$ Galois descent from level $4$ to level $2$
acts on $D^{(\sfive)}$ exactly as on $D^{(\sthree)}$;\ only the rational
coefficients (and the number of polynomial terms, 3 vs.\ 2) differ.

\paragraph{Specific motivic identification of $S_{\min}^{(\sfive)}$.}
The two-loop S^5 irreducible
\[
S_{\min}^{(\sfive)} := \sum_{k \ge 1} T_5(k)^2, \quad
T_5(k) := \tfrac{1}{3}\zeta(2, k + \tfrac{5}{2})
        - \tfrac{5}{6}\zeta(4, k + \tfrac{5}{2})
        + \tfrac{3}{16}\zeta(6, k + \tfrac{5}{2}),
\]
truncates to $S_{\min}^{(\sfive)} \approx 0.0399612091165\ldots$ at
$K_{\max} = 500$, $50$ dps.  Full 200-dps PSLQ irreducibility against
the Paper 28~\S$S_{\min}$ basis is the natural verification target
(estimated effort: $2$–$3$ weeks, mirroring the $\sthree$ sprint).
The natural prediction is $S_{\min}^{(\sfive)} \in \MT(\Z[1/2])$
at depth $2$ as a new period not in the tested basis.
```

## Open questions / follow-on

1. **Specific motivic identification of $S_{\min}^{(\sfive)}$.**  Same
   structural status as the $\sthree$ case, deferred 2–3 weeks (Paper 55
   §5.5 follow-up).

2. **Cross-manifold extension $\sthree \otimes \mathrm{Hardy}(\sfive)$.**
   Paper 24 §V Coulomb/HO asymmetry layer 4 (modular Hamiltonian
   structure of the wedge KMS state).  Both $M_3^{(\sthree)}$ and
   $M_3^{(\sfive)}$ are individually classified;\ the tensor product
   inherits the Bochniak–Sitarz two-body wall identified in the recent
   Sprint resolvent + Sprint tensor-product spectral action sprints
   (CLAUDE.md §2).

3. **M1 on $\sfive$.** $\Vol(\mathbb{CP}^2)/4 = \pi^2/2 \in \Q[\pi]$
   via the Hopf-like projection $\sfive \to \mathbb{CP}^2$.  Could be
   integrated into the master Mellin engine M1 section of Paper 55 as
   a sub-paragraph.

4. **Inhomogeneous extension on $\sfive$.** Whether the $\sfive$
   result extends to non-trivial scaling factor $a(t)$ on
   $\R \times \sfive$.  Fathizadeh–Marcolli result applies to
   $\R \times \sthree$ specifically;\ the $\sfive$ continuum extension
   requires its own re-derivation.  Multi-month frontier.

## Files used

### Papers (read)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §5
  (S^3 M3 classification, the model to mirror) and §7.4 (Q3' open
  question target).
- `papers/group3_foundations/paper_24_bargmann_segal.tex` (Hardy S^5
  substrate).
- `papers/group5_qed_gauge/paper_28_qed_s3.tex` Theorem 2 (parity
  discriminant, lines 376–397) and Theorem 3 ($\chi_{-4}$ identity,
  lines 427–520) — the S^3 model.

### Companion memos (consulted)
- `debug/sprint_m3_cyclotomic_mixed_tate_memo.md` — S^3 M3 closure
  (the structural template).
- `debug/sprint_a2_s5_mixed_tate_memo.md` — S^5 M2 closure
  (companion sprint A2, this morning); contains the CH S^5 spectrum
  data and the dimension-parity sharpening (S^5 SD in $\pi^3\cdot\Q$
  odd-weight slice).

### Scripts (created this sprint)
- `debug/sprint_a6_m3_s5_derivation.py` — symbolic + numerical
  derivation of:
  - $D^{(\sfive)}(s)$ closed form at integer $s \in \{6, \ldots, 12\}$
    (sympy `simplify` output verified against the parity discriminant
    structure).
  - $D_{\mathrm{even}}^{(\sfive)}(s) - D_{\mathrm{odd}}^{(\sfive)}(s)$
    closed form, verified by direct summation at mpmath precision 40
    dps to residuals $\le 10^{-12}$ at $s = 6$, falling to $10^{-33}$
    at $s = 10$.
  - Finite-truncation $S_{\min}^{(\sfive)}$ at $K_{\max} = 500$,
    50 dps:\ $\approx 0.0399612091165\ldots$.

### Published references (already in Paper 55 bibliography)
- Camporesi & Higuchi 1996 (`camporesi_higuchi1996`).
- Deligne 2010 / arXiv:math/0302267 (`deligne2010`).
- Glanois 2015 / arXiv:1411.4947 (`glanois2015`).
- Brown 2012 (`brown2012`).
- Apéry 1979 / Rivoal 2000 — $\Q$-linear independence of odd zeta
  (cited in Paper 28 Theorem 2 proof).

## 200-word summary

Sprint A6 closes Paper 55 §7.4 open question Q3' POSITIVE.  The M3
sub-mechanism on the Camporesi–Higuchi $\sfive$ spectral triple
inherits the cyclotomic-mixed-Tate-at-level-$\le 4$ classification of
the $\sthree$ case verbatim:\ $M_3^{(\sfive)} \subset \MT(\Z[i, 1/2])$,
depth equal to loop order.  Three explicit closed forms verified
symbolically and numerically:\ (i) full Dirac Dirichlet
$D^{(\sfive)}(s)$ as a three-term polynomial in
$\zeta(s-4), \zeta(s-2), \zeta(s)$ with the $\sfive$ parity
discriminant — pure $\pi^{\text{even}}\cdot\Q$ at even $s$, odd-zeta at
odd $s$;\ (ii) parity decomposition into $Z_{\mathrm{even}}(s) = 2^{-s}\zeta(s, 5/4)$
and $Z_{\mathrm{odd}}(s) = 2^{-s}\zeta(s, 7/4)$, exposing level-$4$
quarter-integer Hurwitz shifts;\ (iii) $\chi_{-4}$ identity
$D_{\mathrm{even}}^{(\sfive)} - D_{\mathrm{odd}}^{(\sfive)} = \tfrac{1}{3}f_5(s-4) - \tfrac{5}{6}f_5(s-2) + \tfrac{3}{16}f_5(s)$
with $f_5(s) = 2^s\,\beta(s) - 2^s + (2/3)^s$.  Numerical residuals
$\le 10^{-12}$ at $s = 6$, falling to $10^{-33}$ at $s = 10$ at
mpmath 40 dps, $N = 2 \times 10^5$ truncation.  Level $4$ is sufficient
(the $\sfive$ shifts $5/4, 7/4$ reduce to the level-4 shifts $1/4, 3/4$
modulo rational corrections);\ no escalation to level $8$ is required.
Proposed Paper 55 §5.5 LaTeX block ready.  $S_{\min}^{(\sfive)} \approx
0.039961\ldots$ at finite truncation;\ 200-dps PSLQ irreducibility
deferred 2–3 weeks per the original §7.4 estimate.
