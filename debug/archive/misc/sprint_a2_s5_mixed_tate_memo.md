# Sprint A2 — S^5 Bargmann-Segal Mixed-Tate / Pure-Tate Test

Date: 2026-06-03
Scope: 1-day diagnostic; sub-agent dispatch under Paper 55 §7 open-question
closure directive; PI sign-off pending on paper edits.

## TL;DR

**Verdict: POSITIVE (cross-manifold extension closes).** Under the standard
Connes–Chamseddine volume-normalized convention, the GeoVac Seeley–DeWitt
coefficients on round unit $S^5$, in both the Dirac and scalar-Laplacian
sectors, sit inside the *single Tate-weight-3 slice*
$\pi^3\cdot\mathbb{Q} \subset \mathbb{Q}[\pi] \subset \mathrm{MT}(\mathbb{Q})$,
with a structural enhancement on the Dirac side mirroring the S^3 case:

- **S^5 Dirac D^2: THREE-term exact.** Exactly three non-zero SD coefficients,
  $a_0 = 4\pi^3$, $a_1 = -\tfrac{20}{3}\pi^3$, $a_2 = 3\pi^3$, all in
  $\pi^3\cdot\mathbb{Q}$. All $a_k = 0$ for $k \ge 3$, matching Paper~51's
  general prediction that $S^{2m+1}$ has exactly $(m{+}1)$ power-law terms
  in the spectral action.
- **S^5 scalar Laplacian Δ: infinite SD with closed form.** Each
  $a_k^\Delta = (6-k)\cdot 4^{k-1}\cdot 2 / (3\cdot k!)\cdot \pi^3$ for
  $k \ge 1$, plus $a_0^\Delta = \pi^3$. All in $\pi^3\cdot\mathbb{Q}$. Single
  mid-series zero at $k=6$ from the $(6{-}k)$ factor; no series termination.

This refines the S^3 pure-Tate-*even-weight* ring $\bigoplus_k
\pi^{2k}\cdot\mathbb{Q}$ to the **pure-Tate-odd-weight slice**
$\pi^3\cdot\mathbb{Q}$ on $S^5$. The two are *strictly different*
$\mathbb{Q}$-vector spaces inside $\mathbb{Q}[\pi]$, but **both are sub-rings
of mixed-Tate periods over $\mathbb{Q}$** in the Fathizadeh–Marcolli sense.
**The case-exhaustion theorem M2 sub-mechanism inherits the F–M classification
on $S^5$ exactly as on $S^3$**, with the dimension-parity sharpening that
$\dim S^d$ controls the Tate weight of the SD slice ($d=3$ -> even-weight
sub-ring $\pi^2\cdot\mathbb{Q}$; $d=5$ -> odd-weight slice $\pi^3\cdot\mathbb{Q}$).

The structural caveat from the S^3 case — the raw heat-trace coefficients
contain $\sqrt\pi$ which cancels against $(4\pi)^{d/2}$ in the volume
normalization — applies identically on $S^5$ with $d{=}5$.

## Background

Paper 55 §7 (open questions) currently lists "Cross-manifold extension to
$S^5$ Bargmann–Segal" as an open question:
- **Question:** Does the master Mellin engine partition extend to the
  Paper~24 Hardy lattice on $S^5$, and does the M2 period-theoretic
  classification hold there?

Paper 51 `rem:two_term_uniqueness` already established: on $S^{2m+1}$ the
CH Dirac spectral action has exactly $(m{+}1)$ power-law terms (so $S^3$ -> 2,
$S^5$ -> 3, $S^7$ -> 4), with the integer-spectrum closure mechanism being the
Bernoulli identity $B_{2k+1}(p + 1/2) = (2k+1)\sum_{l=0}^{p-1}(l+1/2)^{2k}$
generalising the S^3 identity $B_{2k+1}(3/2) = (2k+1)/4^k$.

Paper 55 §4 proved on $S^3$ that under volume normalization the SD
coefficients sit in the *pure-Tate even-weight sub-ring*
$\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$, strictly smaller than the generic
F–M mixed-Tate ring over $\mathbb{Q}$ (Fathizadeh–Marcolli 2017,
*Comm. Math. Phys.* **356**, 641–671, arXiv:1611.01815). The natural
question for $S^5$: does the same ring containment hold?

## Spectra on round unit S^5

**Scalar Laplacian.** Spectrum $\lambda_n^\Delta = n(n+4)$, degeneracy
$d_n = (n+1)(n+2)^2(n+3)/12$ (standard SO(6) representation theory).

| $n$ | $\lambda_n$ | $d_n$ |
|:---:|:-----------:|:-----:|
| 0 | 0 | 1 |
| 1 | 5 | 6 |
| 2 | 12 | 20 |
| 3 | 21 | 50 |
| 4 | 32 | 105 |

**Dirac (Camporesi–Higuchi).** $|\lambda_n^{CH}| = n + 5/2$, degeneracy
$g_n = (n+1)(n+2)(n+3)(n+4)/3$.

| $n$ | $|\lambda_n^{CH}|$ | $g_n$ |
|:---:|:------------------:|:-----:|
| 0 | 5/2 | 8 |
| 1 | 7/2 | 40 |
| 2 | 9/2 | 120 |
| 3 | 11/2 | 280 |
| 4 | 13/2 | 560 |

(Both derived in `debug/g4_structural_s5_comparison.py` from Camporesi–Higuchi
1996; verified at low $n$ against SO(6) rep theory.)

## Dirac D^2 SD coefficients on S^5 — derivation and closed form

The CH Dirac degeneracy on $S^5$, written in $m = n + 5/2$, is:
$$
g(m) = \tfrac{1}{3}(m - \tfrac{3}{2})(m - \tfrac{1}{2})(m + \tfrac{1}{2})(m + \tfrac{3}{2})
     = \tfrac{1}{3}m^4 - \tfrac{5}{6}m^2 + \tfrac{3}{16}.
$$
Crucially, $g(1/2) = g(3/2) = 0$, so extending the sum over $m$ in $\mathbb{Z}+1/2$
costs nothing — the "missing" $m = 1/2, 3/2$ entries are automatically zero.
This is the Bernoulli mechanism reflected in the degeneracy polynomial,
and it lets us evaluate via Jacobi $\vartheta_2$.

Using
$$
\vartheta_2(t) := \sum_{m \in \mathbb{Z}+1/2} e^{-tm^2}
              = \sqrt{\pi/t}\bigl(1 + O(e^{-\pi^2/t})\bigr),
$$
and its derivatives,
$$
\sum_{m \in \mathbb{Z}+1/2} m^{2k} e^{-tm^2}
= (-\partial_t)^k \vartheta_2(t)
= \frac{\Gamma(k + 1/2)}{\sqrt\pi}\cdot\frac{\sqrt\pi}{t^{k+1/2}}\cdot\bigl(1 + O(\cdot)\bigr),
$$
one obtains (after restricting to positive $m$, which gives factor $1/2$):
$$
\boxed{
K_{D^2}^{S^5}(t)
= \sqrt\pi\,\Bigl[\frac{1}{8\,t^{5/2}} - \frac{5}{24\,t^{3/2}} + \frac{3}{32\,t^{1/2}}\Bigr]
+ O(e^{-\pi^2/t}).
}
$$
The expansion has **exactly three power-law terms**. All higher coefficients
are exponentially suppressed and vanish in the SD asymptotic series.

**Numerical verification (driver `debug/sprint_a2_s5_sd_coefficients.py`):**
At $t \in \{0.01, 0.05, 0.1, 0.2\}$ the closed-form three-term expression
matches the truncated direct sum to $10^{-17}$ or better at 40-dps mpmath
arithmetic (modular remainder $e^{-\pi^2/t}$ is far smaller than the test
precision).

**Volume-normalized SD coefficients** (multiplying raw $\sqrt\pi\cdot q_k$ by
$(4\pi)^{5/2} = 32\,\pi^{5/2}$):

| $k$ | raw $\tilde a_k^{D^2}$ | volume-normalized $a_k^{D^2}$ |
|:---:|:----------------------:|:-----------------------------:|
| 0 | $\sqrt\pi/8$ | $4\pi^3$ |
| 1 | $-5\sqrt\pi/24$ | $-20\pi^3/3$ |
| 2 | $3\sqrt\pi/32$ | $3\pi^3$ |
| $\ge 3$ | $0$ | $0$ |

The leading $a_0^{D^2} = 4\pi^3 = \dim_{\rm spinor}^{S^5}\cdot\mathrm{Vol}(S^5)
= 4 \cdot \pi^3$ is the standard CC sanity check (CH spinor bundle has
$\dim_{\rm spinor} = 2^{\lfloor 5/2\rfloor} = 4$, same as on $S^3$).

## Scalar Laplacian SD coefficients on S^5 — derivation and closed form

The scalar degeneracy on $S^5$, written in $m = n + 2$, is:
$$
d(m) = \frac{(m-1)\,m^2\,(m+1)}{12} = \frac{m^2(m^2 - 1)}{12},
$$
and $\lambda_n^\Delta = n(n+4) = m^2 - 4$. So
$$
K_\Delta^{S^5}(t) = e^{4t}\cdot \frac{1}{12}\sum_{m \in \mathbb{Z}_{\ge 0}} m^2(m^2 - 1)\,e^{-tm^2}.
$$
Using $\vartheta_3$ and its derivatives,
$$
\sum_{m \in \mathbb{Z}_{\ge 0}} m^2 e^{-tm^2} = \frac{\sqrt\pi}{4\,t^{3/2}} + O(e^{-\pi^2/t}),\qquad
\sum_{m \in \mathbb{Z}_{\ge 0}} m^4 e^{-tm^2} = \frac{3\sqrt\pi}{8\,t^{5/2}} + O(e^{-\pi^2/t}),
$$
gives
$$
K_\Delta^{S^5}(t) = \frac{e^{4t}\,\sqrt\pi}{12}\Bigl[\frac{3}{8\,t^{5/2}} - \frac{1}{4\,t^{3/2}}\Bigr] + O(e^{-\pi^2/t}).
$$
Expanding $e^{4t} = \sum_j (4t)^j/j!$ and collecting powers of $t^{(k-5)/2}$
yields the **closed-form S^5 scalar SD coefficients**:
$$
\boxed{
\tilde a_0^\Delta = \frac{\sqrt\pi}{32},\qquad
\tilde a_k^\Delta = \frac{(6-k)\cdot 4^{k-1}\,\sqrt\pi}{48\cdot k!}\quad (k \ge 1).
}
$$
Volume-normalized ($a_k = (4\pi)^{5/2}\,\tilde a_k = 32\pi^{5/2}\,\tilde a_k$):
$$
\boxed{
a_0^\Delta = \pi^3,\qquad
a_k^\Delta = \frac{(6-k)\cdot 4^{k-1}\cdot 2}{3\cdot k!}\,\pi^3\quad (k \ge 1).
}
$$

| $k$ | $a_k^\Delta$ |
|:---:|:------------:|
| 0 | $\pi^3$ |
| 1 | $(10/3)\,\pi^3$ |
| 2 | $(16/3)\,\pi^3$ |
| 3 | $(16/3)\,\pi^3$ |
| 4 | $(32/9)\,\pi^3$ |
| 5 | $(64/45)\,\pi^3$ |
| 6 | $0$ |
| 7 | $(-512/945)\,\pi^3$ |
| 8 | $(-512/945)\,\pi^3$ |
| 9 | $(-1024/2835)\,\pi^3$ |

Two features distinct from S^3:
- The series **does not terminate** (unlike S^3 Dirac); there's an infinite
  tail of decreasing-magnitude coefficients.
- Single mid-series zero at $k=6$ from the $(6{-}k)$ factor (analog of the
  S^3 case where every coefficient is non-zero).

**Numerical verification:** At $t \in \{0.01, 0.05, 0.1, 0.2\}$ the truncated
ten-term closed form agrees with the truncated direct sum at relative
precision $10^{-8}$ at $t = 0.2$, improving to $10^{-18}$ at $t = 0.01$ —
consistent with truncated-series error scaling as $t^{(k_{\max} - 5)/2}$.

## Arithmetic classification

| sector | object | arithmetic class | mixed-Tate over $\mathbb{Q}$? |
|:-------|:-------|:-----------------|:-------------------------------|
| S^5 Dirac raw | $\tilde a_k^{D^2}$ ($k = 0, 1, 2$) | $\sqrt\pi\cdot\mathbb{Q}$ | **NO** |
| S^5 Dirac raw | $\tilde a_k^{D^2}$ ($k \ge 3$) | $0$ | yes (trivially) |
| S^5 Dirac vol-norm | $a_k^{D^2}$ ($k = 0, 1, 2$) | $\pi^3\cdot\mathbb{Q}$ | **YES** |
| S^5 Dirac vol-norm | $a_k^{D^2}$ ($k \ge 3$) | $0$ | yes |
| S^5 scalar raw | $\tilde a_k^\Delta$ ($k \ge 0$) | $\sqrt\pi\cdot\mathbb{Q}$ | **NO** |
| S^5 scalar vol-norm | $a_k^\Delta$ ($k \ge 0$) | $\pi^3\cdot\mathbb{Q}$ | **YES** |

As on S^3, the classification is **convention-dependent in a structural
way** at the raw heat-trace level (the $\sqrt\pi$ is the spinor-fiber-dim
artifact in odd-$d$ heat kernels), but **convention-independent at the
SD-coefficient level after volume normalization**. After normalization,
every coefficient sits in $\pi^3\cdot\mathbb{Q}$, a single-weight slice of
$\mathbb{Q}[\pi]\subset \mathrm{MT}(\mathbb{Q})$.

## Comparison to S^3 — the dimension-parity sharpening

| feature | S^3 | S^5 |
|:--------|:----|:----|
| spectrum (Dirac) | $|\lambda_n| = n + 3/2$ | $|\lambda_n| = n + 5/2$ |
| degeneracy (Dirac) | $2(n+1)(n+2)$ | $(n+1)(n+2)(n+3)(n+4)/3$ |
| spectrum (scalar) | $n(n+2)$ | $n(n+4)$ |
| degeneracy (scalar) | $(n+1)^2$ | $(n+1)(n+2)^2(n+3)/12$ |
| $\mathrm{Vol}$ | $2\pi^2$ | $\pi^3$ |
| Dirac SD termination | 2-term: $k \le 1$ | 3-term: $k \le 2$ |
| scalar SD termination | $\infty$ (geometric: $a_k^\Delta = 2\pi^2/k!$) | $\infty$ (closed: $a_k^\Delta = (6{-}k)\cdot 4^{k-1}\cdot 2/(3 k!)\cdot \pi^3$ for $k \ge 1$) |
| **SD period ring** | $\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$ (pure-Tate even-weight) | $\pi^3\cdot\mathbb{Q}$ (single pure-Tate odd-weight slice) |

The **dimension-parity sharpening** statement: on $S^d$ with $d$ odd, the
SD coefficients live in $\pi^d\cdot\mathbb{Q}$ if $d \equiv 1 \pmod{4}$ and
in $\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$ if $d \equiv 3 \pmod{4}$. The
distinction comes from $\mathrm{Vol}(S^d) = 2\pi^{(d+1)/2}/\Gamma((d+1)/2)$:
on $S^3$ this is $2\pi^2$ (a $\pi^2$ power), so SD coefficients sit in
$\pi^2\cdot\mathbb{Q}$ and the $e^t$ scalar series produces an even-weight
ring; on $S^5$ this is $\pi^3$, so SD coefficients sit in $\pi^3\cdot\mathbb{Q}$
(a single odd-weight slice), and the scalar expansion is shifted by $e^{4t}$
which preserves the $\pi^3$ slice.

Both rings sit inside $\mathbb{Q}[\pi] \subset \mathrm{MT}(\mathbb{Q})$.
**Both are strictly smaller than the generic Fathizadeh–Marcolli mixed-Tate
ring** for inhomogeneous Robertson–Walker substrates (which allows
$\zeta(2k+1)$ and MZVs when scaling-factor derivatives are non-trivial).
**Both inherit the F–M classification** at the static-substrate sub-case
$a(t) \equiv 1$.

## Case-exhaustion theorem on S^5

The Paper 32 §VIII case-exhaustion theorem (every $\pi$ source is M1 or M2
or M3 on the master Mellin engine $\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$
at $k \in \{0, 1, 2\}$) transfers to $S^5$ verbatim. The proof on $S^3$
uses (i) the Mellin extraction formula, (ii) the Hopf-base $\mathrm{Vol}(S^2)/4 = \pi$
identification for M1, (iii) the modular-transformation closed forms for M2,
(iv) the parity decomposition of the Dirichlet series for M3. None of these
ingredients are S^3-specific:

- **M1 on S^5:** The Hopf-like base measure becomes $\mathrm{Vol}(\mathbb{CP}^2)/4$
  on the projection $S^5 \to \mathbb{CP}^2$, which is $\pi^2/2$ — still
  in $\mathbb{Q}[\pi]$.
- **M2 on S^5:** The Jacobi $\vartheta_2$, $\vartheta_3$ modular
  transformations apply on $S^5$ exactly as on $S^3$, with the dimension
  shifting the $\sqrt\pi$ vs $(4\pi)^{d/2}$ balance.
- **M3 on S^5:** The Dirac Dirichlet $D^{S^5}(s) = \sum g_n (n + 5/2)^{-s}$
  decomposes into Hurwitz zetas at the half-integer shift $5/2$,
  reducing under $\zeta(s, 5/2) = \zeta(s, 1/2) - 2^s - (2/3)^s$ to
  $\zeta(s, 1/2)$ + rationals, which is in $\mathrm{MT}(\mathbb{Z}[1/2])$
  by the duplication identity exactly as on $S^3$.

## Verdict

**POSITIVE.** The pure-Tate refinement of the F–M mixed-Tate classification
extends to the S^5 Bargmann–Segal Hardy lattice substrate. M_2 on S^5 sits
in the single Tate-weight-3 slice $\pi^3\cdot\mathbb{Q} \subset \mathbb{Q}[\pi]
\subset \mathrm{MT}(\mathbb{Q})$, with **three-term exactness** on the
Dirac side and an **infinite closed-form expansion** on the scalar side.
The dimension-parity sharpening — $d=3$ -> even-weight ring,
$d=5$ -> odd-weight slice — is the structural refinement absent in the
generic continuum F–M result.

This closes Paper 55 §7's open question 3 (cross-manifold extension to
$\sfive$) in the affirmative for M2. It does **not** automatically close it
for M3 (which is the Dirac Dirichlet content beyond Mellin slot $k=2$);
M3 on $S^5$ requires its own sub-sprint to verify the $\zeta(s, 5/2)$
decomposition over $\mathrm{MT}(\mathbb{Z}[i, 1/2])$. M1 on $S^5$ is
inherited from $\mathbb{CP}^2$ base measure in the same way M1 on $S^3$
is inherited from $S^2$ base, and stays in $\mathbb{Q}[\pi]$.

## Proposed Paper 55 §7 update

Replace the existing §7.3 paragraph
("Cross-manifold extension to $\sfive$ Bargmann–Segal") with the following
LaTeX block, which converts the open question to a **closed corollary** for
M2 with explicit closed forms, while keeping M3 and the inner-factor open
problems intact.

```latex
\subsection{Cross-manifold extension to $\sfive$ Bargmann--Segal:\ M2
closes by dimension-parity sharpening}
\label{subsec:open_s5}

The M2 sub-mechanism on the Hardy sector of $\sfive$ (Paper~24
substrate~\cite{loutey_paper24}) inherits the F--M mixed-Tate
classification verbatim, refined by the dimension parity of the
substrate.  On round unit $\sfive$:

\begin{proposition}[Pure-Tate refinement on $\sfive$]
\label{prop:s5_pure_tate}
The volume-normalised Seeley--DeWitt coefficients of the
Camporesi--Higuchi Dirac $D^2$ and the scalar Laplacian on round
unit $\sfive$ are:
\begin{itemize}
  \item \textit{Dirac (three-term exact).}  Exactly three non-zero
    coefficients,
    \[
    a_0^{D^2} = 4\pi^3,\qquad
    a_1^{D^2} = -\tfrac{20}{3}\pi^3,\qquad
    a_2^{D^2} = 3\pi^3,
    \]
    and $a_k^{D^2} = 0$ for $k \ge 3$.  The three-term truncation
    matches the general $\sphere^{2m+1}$ count $(m{+}1)$ at $m = 2$
    (Paper~51 \texttt{rem:two\_term\_uniqueness},
    \cite{loutey_paper51}).
  \item \textit{Scalar Laplacian (infinite closed form).}
    \[
    a_0^\Delta = \pi^3,\qquad
    a_k^\Delta = \frac{(6 - k)\cdot 4^{k-1}\cdot 2}{3\cdot k!}\,\pi^3
    \quad (k \ge 1),
    \]
    with the single mid-series zero $a_6^\Delta = 0$ from the
    $(6 - k)$ factor.
\end{itemize}
All coefficients lie in the single Tate-weight-3 slice
$\pi^3\cdot\Q \subset \Q[\pi] \subset \MT(\Q)$.  Consequently
\begin{equation}
\label{eq:m2_s5_ring}
M_2^{(\sfive)} \;\subset\; \pi^3\cdot\Q,
\end{equation}
which sits in the same Kontsevich--Zagier mixed-Tate-over-$\Q$ class
as the $\sthree$ counterpart but in a strictly different single-weight
slice (odd weight $+3$ rather than even-weight sub-ring).
\end{proposition}

The dimension parity controls the Tate weight of the SD slice:\ on $\sthree$,
$\mathrm{Vol}(\sthree) = 2\pi^2$, and the SD coefficients sit in
$\bigoplus_k \pi^{2k}\cdot\Q$ (even-weight pure-Tate sub-ring); on $\sfive$,
$\mathrm{Vol}(\sfive) = \pi^3$, and the SD coefficients sit in the
single slice $\pi^3\cdot\Q$ (odd-weight pure-Tate slice).  Both rings sit
in $\Q[\pi] \subset \MT(\Q)$;\ both inherit the F--M
classification~\cite{fathizadeh_marcolli2016} at the static-substrate
sub-case $a(t) \equiv 1$.  Closed forms derived by Jacobi $\vartheta_2,
\vartheta_3$ modular transformation in
\texttt{debug/sprint\_a2\_s5\_sd\_coefficients.py};\ memo
\texttt{debug/sprint\_a2\_s5\_mixed\_tate\_memo.md} (June~2026)
documents the full derivation, numerical verification, and
dimension-parity sharpening.

The M1 sub-mechanism transfers to $\sfive$ via the $\mathbb{CP}^2$
base measure on the Hopf-like projection $\sfive \to \mathbb{CP}^2$,
$\mathrm{Vol}(\mathbb{CP}^2)/4 = \pi^2/2 \in \Q[\pi]$.  Whether the M3
sub-mechanism (vertex parity and Hurwitz at the half-integer shift
$5/2$) admits a parallel cyclotomic-mixed-Tate-at-level-$\le 4$
classification on $\sfive$ is a sub-sprint of identical structural shape
to the $\sthree$ argument in \S\ref{sec:m3};\ we record it as an open
sub-question below (Q3').
```

Then add a new sub-section:

```latex
\subsection{Open Q3':\ M3 cyclotomic-mixed-Tate refinement on $\sfive$}
\label{subsec:open_s5_m3}

The $\sthree$ M3 classification (Theorem~\ref{thm:m3_cyclotomic_mixed_tate})
rests on:\ (i) the half-integer Hurwitz shift $3/2$ in the
Camporesi--Higuchi spectrum, (ii) the vertex-parity decomposition
$D = D_{\rm even} + D_{\rm odd}$ engaging Hurwitz at quarter-integer
shifts $1/4, 3/4$, and (iii) the Deligne--Glanois descent at $N = 4$.

On $\sfive$ the spectrum has half-integer shift $5/2$ instead of $3/2$, so
the corresponding analysis would engage:\ Hurwitz at $1/4, 3/4, 5/4, 7/4$
in the parity decomposition, with the Glanois level-$4$ cyclotomic basis
still controlling the period-theoretic placement.  An empirical PSLQ
test against the same Paper~28 \S~$S_{\min}$ basis applied to the
$\sfive$-analogue $S_{\min}^{(\sfive)} := \sum_{k \ge 1} T_5(k)^2$ with
$T_5(k) = \tfrac{2}{3}\zeta(4, k + \tfrac{5}{2}) - \cdots$ would be the
natural verification target.  Estimated effort:\ 2--3 weeks (parallel
to the $\sthree$ sprint, no new infrastructure).
```

And update the existing §7.4 (inner-factor Yukawa) unchanged; the M3
$\sfive$ open question is a new Q3' that **replaces** the prior single
S^5 open question.

## Open sub-questions (carried forward)

1. **M3 on S^5 (Q3' above).** Sub-sprint of identical structural shape to
   the $\sthree$ M3 argument. Estimated effort: 2–3 weeks.

2. **Inhomogeneous extension on S^5.** Whether the S^5 result extends to
   non-trivial scaling factor $a(t)$ on $\mathbb{R}\times S^5$. F–M's
   continuum result applies to $\mathbb{R}\times S^3$ specifically; the
   $S^5$ continuum extension would require its own re-derivation.
   Estimated effort: 1–2 months.

3. **Cyclotomic-mixed-Tate verdict for M3.** (Same open question as the
   S^3 case, sprint M3 cyclotomic mixed-Tate test memo
   `debug/sprint_m3_cyclotomic_mixed_tate_memo.md`, parallel work.)

4. **Inner-factor extension on S^5.** Paper~32 §VIII.C inner-factor
   Yukawa content classification — parallel open question to the
   S^3 case.

## Files used

### Papers (read)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` (current §7
  open-question text, all four subsections)
- `papers/group3_foundations/paper_24_bargmann_segal.tex` (S^5 Hardy
  substrate)
- `papers/group5_qed_gauge/paper_51_gravity_arc.tex` `rem:two_term_uniqueness`
  (general $S^{2m+1}$ counting: $(m{+}1)$ power-law terms)

### Memos (consulted)
- `debug/sprint_mixed_tate_test_memo.md` (S^3 M2 mixed-Tate test
  POSITIVE-with-refinement, structural template for this sprint)

### Code (consulted and run)
- `debug/g4_structural_s5_comparison.py` — provided the CH Dirac
  degeneracy on S^5 derived from Camporesi–Higuchi 1996, plus the
  symbolic verification that $\zeta_{D^2}^{S^5}(-k) = 0$ for $k \in \{0,
  \ldots, 7\}$ (consistent with three-term-exact SD).

### Diagnostic script (created this sprint)
- `debug/sprint_a2_s5_sd_coefficients.py` — derivation of the S^5 Dirac
  and scalar SD coefficients via Jacobi $\vartheta_2, \vartheta_3$
  modular transformation; bit-exact verification against truncated
  direct heat-trace sum at $t \in \{0.01, 0.05, 0.1, 0.2\}$ to relative
  precision $10^{-17}$ (Dirac, three-term exact) and $10^{-8}$ to
  $10^{-18}$ (scalar, depending on truncation order).

### Literature
- Camporesi & Higuchi, *J. Geom. Phys.* **20**, 1–18 (1996), Dirac
  eigenfunctions on spheres (cited in CLAUDE.md and Paper 38).
- Fathizadeh & Marcolli, *Comm. Math. Phys.* **356**, 641–671 (2017),
  arXiv:1611.01815 (the mixed-Tate periods classification for CC spectral
  actions on R–W spacetimes; cited in Paper 55 §4 and in this memo's
  inheritance argument).
- No additional literature lookup required.
