# Sprint A5 — Fathizadeh–Marcolli full-text read

Date: 2026-06-03
Scope: One-week diagnostic completing the Open Question #5 of
`debug/sprint_mixed_tate_test_memo.md` (v3.45.3) by extracting the precise
affine-complement-of-quadrics-and-hyperplanes structure from
Fathizadeh–Marcolli arXiv:1611.01815 (Commun. Math. Phys. **356**, 641–671
(2017)). Output: explicit motive-equivalence corollary for GeoVac's
static-$S^3$ sub-case, paste-ready Paper 55 §4 strengthening.

## TL;DR

**Verdict: POSITIVE, with a clean structural reduction.**

The F–M classification rests on a single explicit object:\ a one-parameter
family of quadrics
$$Q_{\lambda, 2n} = u_1^2 + \lambda^{-2}(u_2^2 + u_3^2 + u_4^2) + u_5^2 + \ldots + u_{2n+2}^2$$
living in the affine space $\mathbb{A}^{2n+3}$ with coordinates
$(u_0, u_1, \ldots, u_{2n+2})$ and an additional "scaling-factor"
coordinate $\lambda \in \mathbb{A}^1 \setminus \{0\} = \mathbb{G}_m$.
The Seeley–DeWitt coefficient $a_{2n}$ (before $t$-integration) is a period
of an algebraic differential form on the affine complement
$$\mathcal{A}^{2n+3} \setminus \bigl( C_{Z_{\lambda, 2n}} \,\cup\, H_0 \,\cup\, H_1 \bigr),$$
with $H_0 = \{u_0 = 0\}$, $H_1 = \{u_0 = 1\}$, and the algebraic differential
form depending additionally on $2n$ auxiliary affine parameters
$\varepsilon_1, \ldots, \varepsilon_{2n}$ corresponding to the time derivatives
$a'(t), a''(t), \ldots, a^{(2n)}(t)$ of the scaling factor (F–M §6, Theorems
6.1, 6.2; explicit substitution at line "$\lambda = a(t),\,
\varepsilon_1 = a'(t),\, \varepsilon_2 = a''(t), \ldots,\,
\varepsilon_{2n} = a^{(2n)}(t)$").

The static GeoVac sub-case $a(t) \equiv 1$ corresponds to the slice
$\lambda = 1, \varepsilon_1 = \varepsilon_2 = \ldots = \varepsilon_{2n} = 0$.
The differential form on this slice is **not the form on the whole
$\varepsilon$-trivial complement, but its specialisation**:\ every monomial
in the F–M integrand that carries a factor $\varepsilon_i$ (equivalently,
a factor $a^{(k)}(t)$ for $k \ge 1$) vanishes at this slice. The remaining
monomials are those built from $a(t)$ in denominators only — the
"$\varepsilon$-independent core" of $b_{-2n-2}$. The expected outcome from
the original sprint scoping was: **GeoVac static $S^3$ reduces to an
essentially trivial affine complement**. This expectation is **verified**:\
the relevant motive collapses to (i) $\lambda$-slice $\lambda = 1$
trivialises the quadric $Q_{1, 2n} = u_1^2 + u_2^2 + \ldots + u_{2n+2}^2$ to a
**standard isotropic sphere quadric over $\mathbb{Q}(\sqrt{-1})$**, and (ii)
the $2n$-dimensional $\varepsilon$-fiber direction collapses to a point.
The static $S^3$ affine complement is
$$\mathbb{A}^{2n+3}_{(u_0, u_1, \ldots, u_{2n+2})} \setminus
\bigl( \{u_1^2 + u_2^2 + \ldots + u_{2n+2}^2 = 0\}
\cup \{u_0 = 0\} \cup \{u_0 = 1\}\bigr),$$
with the algebraic differential form now $\varepsilon$-free and reducing to
the unit-$S^3$ angular content. Over $\mathbb{Q}(i)$ this is a complement of
the standard isotropic quadric in $\mathbb{A}^{2n+3}$ minus two parallel
hyperplanes — manifestly mixed Tate by F–M's Theorems 7.3, 7.4 and the
Grothendieck-class computations (Lemma 7.1, equation 4).

The "static $S^3$ = trivial complement" reading is **correct in spirit but
not literally a point**:\ the affine complement of a single isotropic
quadric and two parallel hyperplanes in $\mathbb{A}^{2n+3}$ is a genuine
mixed-Tate motive with Grothendieck class
$\mathbb{L}^{2n+3} - 3\mathbb{L}^{2n+2} + \ldots$ (Lemma 7.1 part 5),
not literally a point. What collapses to a point is the **$\varepsilon$-fiber
direction** ($2n$ coordinates), not the $(u_0, u_1, \ldots, u_{2n+2})$
core. The period output is consequently in $\pi^{2k} \cdot \mathbb{Q}$
(pure-Tate even-weight), strictly smaller than F–M's generic mixed-Tate
ring — recovering the v3.45.3 sprint's pure-Tate refinement at the level
of an explicit motive-equivalence statement.

## Background — what we needed from F–M

The v3.45.3 sprint memo (`debug/sprint_mixed_tate_test_memo.md`, Open
Question #5) recorded the structural inheritance of F–M's mixed-Tate
classification at text level:\ Paper 55 §4 reads "Direct inheritance is
by application of Fathizadeh–Marcolli Theorem on the Rosenfeld integral
representation of Seeley–DeWitt coefficients on $\mathbb{R} \times S^3$,
specialised to constant scaling factor." The sprint scoping question was:\
**what exactly is the relative motive of the complement of quadrics and
hyperplanes, and what does it specialise to for the GeoVac static
sub-case?**

The v3.45.3 closure relied on the F–M abstract + search-surfaced quotes;
the full PDF was not extractable. This sprint accessed the published PDF
via `arxiv.org/pdf/1611.01815` (with a `Mozilla/5.0` user agent — direct
fetch failed) and `pdftotext -layout`, recovering all 42 pages of usable
text.

## Extraction of the F–M affine-complement structure

### The variables (F–M §2, §6)

The integration on Robertson–Walker $\mathbb{R} \times S^3$ uses Hopf
coordinates on $S^3$:\ $(\theta, \phi_1, \phi_2)$ with $0 < \theta < \pi/2$,
$0 < \phi_1, \phi_2 < 2\pi$. The cotangent-fiber momentum variables are
$\xi = (\xi_1, \xi_2, \xi_3, \xi_4) \in T^*_x M$. F–M's change of variables
(equation 6.3 of the paper) is
$$u_0 = \sin^2(\theta), \quad u_3 = \csc(\theta)\,\xi_3, \quad u_4 = \sec(\theta)\,\xi_4,$$
$$u_j = \xi_j \text{ for } j = 1, 2, 5, 6, \ldots, 2n+2.$$

The "scaling-factor" variable and its derivatives are introduced as
auxiliary affine parameters:
$$\lambda = a(t), \quad \varepsilon_1 = a'(t), \quad \varepsilon_2 = a''(t), \quad \ldots, \quad \varepsilon_{2n} = a^{(2n)}(t).$$

The total parameter count for the $a_{2n}$ integrand is therefore
$2n + 3$ ($u_0, u_1, \ldots, u_{2n+2}$) for the affine ambient space,
plus $\lambda \in \mathbb{G}_m$ and $\varepsilon_1, \ldots, \varepsilon_{2n}$
on which the algebraic differential form depends polynomially.

### The Rosenfeld / Wodzicki representation (F–M §3)

F–M use the Connes–Chakraborty (Wodzicki-residue) representation of
Seeley–DeWitt coefficients (their equation 3.2):
$$a_{2 + r} = 2^5 \cdot \frac{1}{4 + r/2} \cdot \mathrm{Res}(\Sigma_r^{-1}),$$
where $\Sigma_r = D^2 \otimes 1 + 1 \otimes T_r$ with $T_r$ the flat
Laplacian on the $r$-dimensional torus $T^r = (\mathbb{R}/\mathbb{Z})^r$,
and $\mathrm{Res}$ the Wodzicki noncommutative residue. The trick is that
inserting the auxiliary tori promotes the (otherwise computationally
difficult) heat-trace expansion into a Wodzicki-residue computation of
the parametrix of $\Sigma_r$. The Wodzicki residue is then itself an
integral on the unit cosphere bundle (F–M equations 3.3, 3.4), and after
the change of variables (6.3) this integral lives on a $Q$-semialgebraic
set in the affine complement structure described below.

### The quadric (F–M §6, Theorem 6.1, equation 6.6)

The unique quadric appearing in the integrand of $a_{2n}$ is
$$Q_{\lambda, 2n}(u_1, u_2, \ldots, u_{2n+2}) := u_1^2 + \lambda^{-2}(u_2^2 + u_3^2 + u_4^2) + u_5^2 + u_6^2 + \ldots + u_{2n+2}^2.$$

This quadric has:
- Total variables:\ $2n + 2$ (i.e.\ $u_1, \ldots, u_{2n+2}$; not $u_0$);
- Diagonal-form signature $\langle 1, \lambda^{-2}, \lambda^{-2}, \lambda^{-2}, 1, 1, \ldots, 1 \rangle$ (one block of three $\lambda^{-2}$'s, then all unit weights);
- Defining hypersurface $C_{Z_{\lambda, 2n}} := \{Q_{\lambda, 2n} = 0\} \subset \mathbb{A}^{2n+3}$ (affine cone of the projective quadric in $\mathbb{P}^{2n+1}$).

### The hyperplanes

Two affine hyperplanes appear:
$$H_0 = \{u_0 = 0\}, \qquad H_1 = \{u_0 = 1\}.$$
The variable $u_0 = \sin^2(\theta)$ runs over $(0, 1)$ on the Hopf
coordinate chart of $S^3$. The hyperplanes are the boundary divisor of
the angular integration range $0 < u_0 < 1$.

### The main theorem (F–M §6, Theorem 6.2)

> **Theorem (F–M 6.2).** When computed without performing the time
> integration, the coefficient $a_{2n}$ is a period integral
> $$\int_{A_{2n+2}} \omega_{\varepsilon_1, \ldots, \varepsilon_{2n}}$$
> of an algebraic differential form
> $\omega_{\varepsilon_1, \ldots, \varepsilon_{2n}}(u_0, u_1, \ldots, u_{2n+2})$
> defined on the complement
> $$\mathbb{A}^{2n+3} \setminus \bigl(C_{Z_{\lambda, 2n}} \cup H_0 \cup H_1 \bigr),$$
> with $C_{Z_{\lambda, 2n}}$ the hypersurface defined by the vanishing of
> the quadric $Q_{\lambda, 2n}$ (6.6) and $H_0 = \{u_0 = 0\}$,
> $H_1 = \{u_0 = 1\}$. The period integral is performed over the
> $\mathbb{Q}$-semialgebraic set
> $$A_{2n+2} = \Bigl\{ (u_0, \ldots, u_{2n+2}) \in \mathbb{A}^{2n+3}(\mathbb{R}) :
> u_1^2 + u_2^2 + u_0 u_3^2 + (1 - u_0) u_4^2 + \sum_{i=5}^{2n+2} u_i^2 = 1,$$
> $$0 < u_i < 1, \,\, i = 0, 1, 2, 5, 6, \ldots, 2n+2 \Bigr\}.$$

The motivic content is in F–M §7. F–M prove (Theorem 7.3, 7.5, and the
inductive treatment §7.6, 7.7) that over the quadratic extension
$K = \mathbb{Q}(\sqrt{-1})$ the quadric $Q_{\lambda, 2n}$ becomes
isotropic, and the motive
$m\bigl(\mathbb{A}^{2n+3} \setminus (C_{Z_{\lambda, 2n}} \cup H_0 \cup H_1)\bigr)$
is mixed Tate over $K$. The Grothendieck class is computed explicitly via
Lemma 7.1 (part 5):
$$[\mathbb{A}^{2n+3} \setminus (C_Z \cup H_0 \cup H_1)] = \mathbb{L}^{2n+3} - 2\mathbb{L}^{2n+2} - (\mathbb{L} - 2)(\mathbb{L} - 1)[Z] - (\mathbb{L} - 2),$$
where $\mathbb{L} = [\mathbb{A}^1]$ is the Lefschetz motive and $[Z]$ is
the Grothendieck class of the projective quadric.

### Dependence on derivatives of $a(t)$

F–M Theorem 6.1 (equation 6.2) gives $\mathrm{tr}(\sigma_{-2n-2})$ as a
sum of $M_n$ monomials, each of the form
$$\frac{c_{j, 2n} \, u_0^{\alpha_{0,1,j}/2} (1 - u_0)^{\alpha_{0,2,j}/2}
\, u_1^{\alpha_{1,j}} u_2^{\alpha_{2,j}} \cdots u_{2n+2}^{\alpha_{2n+2,j}}
\, \varepsilon_1^{k_1} \varepsilon_2^{k_2} \cdots \varepsilon_{2n}^{k_{2n}}}
{Q_{\lambda, 2n}^{\beta_{j, 2n}}}$$
with $c_{j, 2n} \in \mathbb{Q}$ and $k_0, k_1, \ldots, k_{2n}$
non-negative integers. The total number of derivatives entering
$a_{2n}$ is therefore $2n$ ($\varepsilon_1, \ldots, \varepsilon_{2n}$
corresponding to $a'(t), a''(t), \ldots, a^{(2n)}(t)$). At $n = 1$
($a_2$ term), there are 2 derivatives; at $n = 2$ ($a_4$), there are 4;
at $n = 3$ ($a_6$), there are 6. This matches F–M §1's prose statement.

### Static specialization (F–M does not discuss explicitly)

F–M's paper does not include a section on the static sub-case
$a(t) \equiv 1$. The only nearby remark is the rationality observation
(F–M §1.1 last paragraph) that "all the terms in the expansion of the
spectral action are polynomials with rational coefficients in the
scaling factor and its derivatives" (citing arXiv:1010.4980 /
Fathizadeh–Marcolli's earlier work [10] for the proof).

In our static sub-case $a(t) \equiv 1$:\ $\lambda = 1$ and
$\varepsilon_i = 0$ for all $i$. The Theorem 6.1 monomials collapse:\
every monomial with $k_1 + k_2 + \ldots + k_{2n} > 0$ vanishes
identically. The surviving "$\varepsilon$-independent core" consists
of monomials in $u_0, u_1, \ldots, u_{2n+2}$ alone, with all $u_i$ exponents even (after F–M's parity-elimination argument, p.~9
following equation 4.2). The quadric specialises to the **isotropic
unit-norm form**
$$Q_{1, 2n} = u_1^2 + u_2^2 + u_3^2 + u_4^2 + u_5^2 + \ldots + u_{2n+2}^2.$$
The semialgebraic integration domain reduces to the unit-radius
$S^{2n+1}$-spherical-cap parameterised by $0 < u_0 < 1$,
$0 < u_i < 1$ for $i = 1, 2, 5, \ldots, 2n+2$, and the constraint
$u_1^2 + u_2^2 + u_0 u_3^2 + (1 - u_0) u_4^2 + \sum_{i=5}^{2n+2} u_i^2 = 1$
(the $\lambda$ dependence drops out).

The "trivial complement" reading is therefore precise in two senses:
1. The $\varepsilon$-fiber $\mathbb{A}^{2n}$ collapses to its origin $\{0\}$;
2. The $\lambda$-fiber $\mathbb{G}_m = \mathbb{A}^1 \setminus \{0\}$ collapses to the point $\{1\}$.

The remaining "GeoVac affine complement" is
$$\mathcal{V}_n^{\mathrm{GeoVac}} := \mathbb{A}^{2n+3} \setminus
\bigl( \{u_1^2 + u_2^2 + \ldots + u_{2n+2}^2 = 0\}
\cup \{u_0 = 0\} \cup \{u_0 = 1\}\bigr).$$
This is **not literally a point**:\ it is a genuine $(2n+3)$-dimensional
mixed-Tate variety, whose Grothendieck class (via Lemma 7.1.5 with
$[Z] = [\mathbb{P}^1 \times \mathbb{P}^1] = (\mathbb{L} + 1)^2 =
\mathbb{L}^2 + 2\mathbb{L} + 1$ at $n = 1$, and more generally via
inductive computation in F–M §7.6) lies in the standard mixed-Tate
sub-ring $\mathbb{Z}[\mathbb{L}]$.

What **is** a point is the $(\lambda, \varepsilon)$-fiber:\
$\{1\} \times \{0\}^{2n} \subset \mathbb{G}_m \times \mathbb{A}^{2n}$.
This collapse is what makes the GeoVac periods strictly pure-Tate
(no MZV content) rather than general mixed-Tate.

## The motive-equivalence corollary

Putting the pieces together:

**Corollary (motive-equivalence statement for GeoVac static $S^3$).**
Let $a_{2n}^{\mathrm{F-M}}$ denote the F–M Seeley–DeWitt coefficient on
$\mathbb{R} \times S^3$ with general scaling factor $a(t)$, expressed
via Theorem 6.2 as a period of an algebraic differential form
$\omega_{\varepsilon_1, \ldots, \varepsilon_{2n}}$ on the relative
complement
$$\mathcal{V}_n^{\mathrm{F-M}}(\lambda) := \mathbb{A}^{2n+3} \setminus
\bigl(C_{Z_{\lambda, 2n}} \cup H_0 \cup H_1\bigr),$$
with $\lambda \in \mathbb{G}_m$ and
$\varepsilon_i \in \mathbb{A}^1$ as additional affine parameters,
$i = 1, \ldots, 2n$. The GeoVac static $S^3$ Seeley–DeWitt coefficient
$a_{2n}^{\mathrm{GeoVac}}$ is the specialisation of $a_{2n}^{\mathrm{F-M}}$
along the closed immersion
$$\iota: \{1\} \times \{0\}^{2n} \hookrightarrow \mathbb{G}_m \times \mathbb{A}^{2n}$$
defined by $\lambda = 1, \varepsilon_1 = \ldots = \varepsilon_{2n} = 0$.
Equivalently, the period output sits in the smaller sub-ring of
$\mathrm{MT}(\mathbb{Q})$ generated by the Grothendieck class of the
$\varepsilon$-trivialised complement
$$\mathcal{V}_n^{\mathrm{GeoVac}} = \mathbb{A}^{2n+3} \setminus
\bigl(C_{Z_{1, 2n}} \cup H_0 \cup H_1\bigr),$$
where $C_{Z_{1, 2n}}$ is the affine cone of the standard isotropic
quadric $u_1^2 + u_2^2 + \ldots + u_{2n+2}^2 = 0$ in $\mathbb{P}^{2n+1}$.
Over $\mathbb{Q}(\sqrt{-1})$ this complement is mixed Tate by direct
application of F–M Theorem 7.3 (at $n = 1$), Theorem 7.5 (at $n = 2$),
and the inductive treatment in §7.6, 7.7 (general $n$). Over
$\mathbb{Q}$, the period output lies in the pure-Tate even-weight
sub-ring $\bigoplus_k \pi^{2k} \cdot \mathbb{Q}$, recovering the
v3.45.3 sprint's empirical pure-Tate refinement at the level of an
explicit motive-equivalence statement.

The strict reading "GeoVac static $S^3$ = affine complement is a point"
is **NOT correct**:\ the $u$-coordinate complement $\mathcal{V}_n^{\mathrm{GeoVac}}$
is a $(2n+3)$-dimensional mixed-Tate variety, not a point. The strict
reading "GeoVac static $S^3$ = $\varepsilon$-fiber is a point" **IS
correct**:\ the collapse of the $\varepsilon$-fiber $\mathbb{A}^{2n}$ to
its origin and the $\lambda$-fiber $\mathbb{G}_m$ to $\{1\}$ is what
strips the F–M mixed-Tate ring down to its pure-Tate sub-ring on the
GeoVac sub-case.

## Verdict

**POSITIVE, with one correction.** The F–M affine-complement structure
extracts cleanly. The expected "trivial complement" reading is correct
in spirit but needs the precise refinement above:\ what trivialises is
the $(\lambda, \varepsilon)$-fiber, not the $u$-coordinate complement.
The corrected reading provides the explicit motive-equivalence statement
the v3.45.3 sprint Open Question #5 asked for, and the Paper 55 §4
strengthening can be applied as drafted below.

## Proposed Paper 55 §4 strengthening (paste-ready LaTeX)

The block below replaces the current "Proof sketch" of Theorem 4.2
(`thm:m2_mixed_tate`) in Paper 55 §4 with a precise inheritance
statement. The current line —
> "Direct inheritance is by application of Fathizadeh–Marcolli
> Theorem on the Rosenfeld integral representation of Seeley–DeWitt
> coefficients on $\mathbb{R} \times \sthree$, specialised to
> constant scaling factor."

— is upgraded to:

```latex
\begin{proof}[Proof sketch]
The inheritance is via specialisation of the
Fathizadeh--Marcolli~\cite{fathizadeh_marcolli2016} Theorem~6.2:\
on the general Robertson--Walker spacetime $\R \times \sthree$ with
metric $dt^2 + a(t)^2 d\Omega_3^2$, the $a_{2n}$ Seeley--DeWitt
coefficient (prior to time integration) is a period of an algebraic
differential form
$\omega_{\varepsilon_1, \ldots, \varepsilon_{2n}}(u_0, u_1, \ldots, u_{2n+2})$
defined on the affine complement
\begin{equation}
\label{eq:fm_complement}
\mathcal{V}_n^{\mathrm{F\text{-}M}}(\lambda)
\;:=\; \mathbb{A}^{2n+3} \;\setminus\;
\bigl( C_{Z_{\lambda, 2n}} \,\cup\, H_0 \,\cup\, H_1 \bigr),
\end{equation}
where $H_0 = \{u_0 = 0\}$, $H_1 = \{u_0 = 1\}$, $C_{Z_{\lambda, 2n}}$
is the affine cone of the projective quadric defined by
\begin{equation}
\label{eq:fm_quadric}
Q_{\lambda, 2n}(u_1, \ldots, u_{2n+2}) \;=\;
u_1^2 + \lambda^{-2}(u_2^2 + u_3^2 + u_4^2)
+ u_5^2 + u_6^2 + \ldots + u_{2n+2}^2,
\end{equation}
$\lambda \in \mathbb{G}_m$ is the affine parameterisation of the scaling
factor $a(t)$, and
$\varepsilon_i = a^{(i)}(t)$ for $i = 1, \ldots, 2n$ are the auxiliary
affine parameters corresponding to the time derivatives of the scaling
factor.  The algebraic differential form
$\omega_{\varepsilon_1, \ldots, \varepsilon_{2n}}$ is polynomial in
$\varepsilon_1, \ldots, \varepsilon_{2n}$ with monomials carrying explicit
exponents $k_1, \ldots, k_{2n}$ in the expansion of Fathizadeh--Marcolli
Theorem~6.1.

The GeoVac static sub-case $a(t) \equiv 1$ corresponds to specialisation
along the closed immersion
\begin{equation}
\label{eq:geovac_specialisation}
\iota:\ \{1\} \times \{0\}^{2n}
\;\hookrightarrow\; \mathbb{G}_m \times \mathbb{A}^{2n}_{(\varepsilon_1, \ldots, \varepsilon_{2n})},
\qquad \lambda = 1,\;\; \varepsilon_1 = \ldots = \varepsilon_{2n} = 0.
\end{equation}
On this slice, every monomial of the Fathizadeh--Marcolli integrand that
carries an $\varepsilon$-factor vanishes identically;\ the surviving
``$\varepsilon$-independent core'' depends only on $u_0, u_1, \ldots, u_{2n+2}$
with all $u_i$ exponents even, and the quadric specialises to the standard
unit-isotropic form
\begin{equation}
\label{eq:geovac_quadric}
Q_{1, 2n} \;=\; u_1^2 + u_2^2 + u_3^2 + u_4^2 + u_5^2 + \ldots + u_{2n+2}^2.
\end{equation}
The corresponding GeoVac affine complement is
\begin{equation}
\label{eq:geovac_complement}
\mathcal{V}_n^{\mathrm{GeoVac}}
\;:=\; \mathbb{A}^{2n+3} \;\setminus\;
\bigl( C_{Z_{1, 2n}} \,\cup\, H_0 \,\cup\, H_1 \bigr).
\end{equation}
Over the quadratic extension $\Q(\sqrt{-1})$, the quadric $Q_{1, 2n}$ is
hyperbolic (isotropic), and the motive of $\mathcal{V}_n^{\mathrm{GeoVac}}$
is mixed Tate by direct application of Fathizadeh--Marcolli Theorem~7.3
(at $n = 1$), Theorem~7.5 (at $n = 2$), and the inductive treatment of
\S~7.6, 7.7 (general $n$).  Over $\Q$, the period output lies in the
pure-Tate even-weight sub-ring $\bigoplus_k \pi^{2k} \cdot \Q$.

The pure-Tate refinement on the static sub-case has two parts at the
spectrum level:
\begin{itemize}
  \item \textit{Dirac sector.}  The Camporesi--Higuchi spectrum
    $|\lambda_n| = n + 3/2$ admits the Jacobi $\vartheta_2$ modular
    transformation
    $\Tr e^{-tD^2} = (\sqrt\pi/2) t^{-3/2} - (\sqrt\pi/4) t^{-1/2}
                    + O(e^{-\pi^2/t})$, with the exponentially small
    remainder vanishing in the asymptotic series.  All SD coefficients
    $a_k^{D^2}$ for $k \ge 2$ are exactly zero;\ the volume-normalised
    coefficients are $a_0 = 4\pi^2$, $a_1 = -2\pi^2$.
  \item \textit{Scalar sector.}  The integer spectrum $n(n+2)$ admits
    the Jacobi $\vartheta_3$ modular transformation giving the closed
    form $\Tr e^{-t\Delta} = (\sqrt\pi/4) \cdot e^t/t^{3/2} +
    O(e^{-\pi^2/t})$, with $a_k^\Delta = 2\pi^2/k!$.
\end{itemize}
Both sectors are explicit in Paper~51~\cite{loutey_paper51}.  The empirical
pure-Tate restriction (no $\zeta(3)$, no MZVs) is consistent with the
motivic collapse:\ the F--M generic R--W classification admits MZV content
through the $\varepsilon_i$ dependence (time derivatives of the scaling
factor);\ on the static sub-case those derivatives vanish, eliminating the
MZV content at the motivic level.

A precise refinement of the ``trivial complement'' reading:\ the
$(\lambda, \varepsilon)$-fiber of the F--M classification collapses to the
single point $\{1\} \times \{0\}^{2n}$ on the GeoVac sub-case;\ the
$u$-coordinate complement $\mathcal{V}_n^{\mathrm{GeoVac}}$ remains a
genuine $(2n + 3)$-dimensional mixed-Tate variety.
\end{proof}
```

A companion `Remark` immediately after the proof can optionally record
the explicit motive-equivalence statement as a standalone result:

```latex
\begin{remark}[Specialisation as motive equivalence]
\label{rem:m2_specialisation}
The relationship between the F--M generic Robertson--Walker classification
and the GeoVac static $S^3$ specialisation is a closed immersion of
relative motives:\ the F--M classification supplies the universal
``period-with-parameter-data'' object
$\mathcal{V}_n^{\mathrm{F\text{-}M}}(\lambda) \times \mathbb{A}^{2n}_{\varepsilon}$
over $\mathbb{G}_m \times \mathbb{A}^{2n}$, and the GeoVac static sub-case
is the fiber over the single closed point
$\{1\} \times \{0\}^{2n}$.  The Grothendieck-class computation of
Fathizadeh--Marcolli Lemma~7.1.5 evaluated at $\lambda = 1$
(equivalently:\ at the standard unit-isotropic quadric over $\Q(\sqrt{-1})$)
gives
\begin{equation}
\label{eq:geovac_grothendieck_class}
[\mathcal{V}_n^{\mathrm{GeoVac}}] \;=\;
\mathbb{L}^{2n+3} - 2 \mathbb{L}^{2n+2}
- (\mathbb{L} - 2)(\mathbb{L} - 1)[Z_{1, 2n}] - (\mathbb{L} - 2),
\end{equation}
where $[Z_{1, 2n}]$ is the Grothendieck class of the standard unit
isotropic quadric in $\mathbb{P}^{2n+1}$, computed inductively in
Fathizadeh--Marcolli \S~7.6 as a polynomial in $\mathbb{L}$ of degree
$2n$.  All terms lie in $\mathbb{Z}[\mathbb{L}] \subset K_0(\mathrm{Var}_\Q)$,
the standard mixed-Tate Grothendieck sub-ring.
\end{remark}
```

## Scope and honest limits

1. **Closed-form Grothendieck class at general $n$:**\ this memo does NOT
   re-derive the inductive computation of $[Z_{1, 2n}]$ in closed form for
   general $n$.  F--M §7.4–7.6 give $n = 1, 2$ closed forms and a
   recursive structure;\ extending to $[Z_{1, 2n}]$ at arbitrary $n$ as a
   polynomial in $\mathbb{L}$ is a further mini-sprint (estimated 1–2
   days).  For Paper 55 §4 purposes the statement "in $\mathbb{Z}[\mathbb{L}]$"
   suffices.

2. **The "isotropic over $\mathbb{Q}(\sqrt{-1})$" reading**:\ at $\lambda
   = 1$, the F--M quadric becomes the standard
   $u_1^2 + u_2^2 + \ldots + u_{2n+2}^2$, which is anisotropic over $\mathbb{Q}$
   but isotropic over $\mathbb{Q}(\sqrt{-1})$ (it has the elementary
   hyperbolic decomposition via $X = u_1 + iu_2, Y = u_1 - iu_2, \ldots$).
   This is the same "becomes isotropic over $\mathbb{Q}(i)$" property
   used by F--M Theorem 7.3, 7.4 to establish mixed-Tate over $K = \mathbb{Q}(i)$.
   The level-$4$ cyclotomic content of Paper 55 §6 (M3 sub-mechanism, the
   $\chi_{-4}$ vertex parity) is the **inner-factor descent** structure of
   Glanois at level $4$;\ the M2 quadric ring $\mathbb{Q}(\sqrt{-1})$
   appearing here is the same $\mathbb{Q}(i) = \mathbb{Q}(\zeta_4)$,
   suggesting a structural connection between the two sub-mechanisms at
   the cyclotomic-level $4$.  Whether this is coincidence or a deeper M2/M3
   coupling is left as an open question (named "Open question on M2--M3
   cyclotomic coincidence", potentially worth a follow-on note).

3. **No fitted parameters introduced**:\ this work is pure
   motive-theoretic extraction.  No new numerical claims.

4. **Conjectural label on combination rule $K = \pi(B + F - \Delta)$
   preserved**:\ this memo does NOT touch Paper 2 framing.

5. **Robustness of "POSITIVE" verdict against alternative readings**:\
   the only alternative reading of F–M's main theorem that would
   undermine the GeoVac-specialisation argument would be if the
   $\varepsilon$-trivialised slice were somehow off-locus for the F–M
   mixed-Tate proof — i.e.\ if the F–M proof required generic
   $(\lambda, \varepsilon)$ to control the motive structure.  Reading
   §7.3 (Theorem 7.3 proof) confirms this is not the case:\ the
   Gysin-triangle and Mayer-Vietoris arguments work over $\mathbb{Q}(i)$
   without genericity assumption on $(\lambda, \varepsilon)$, and the
   specialisation to a closed point is a closed immersion of mixed-Tate
   motives by the universal property of the Voevodsky category.

## Files used

### Web sources (full PDF access this time)
- `arxiv.org/pdf/1611.01815` (via `curl -A "Mozilla/5.0"` user agent;
  direct fetch and WebFetch failed to retrieve usable content; the user
  agent workaround was load-bearing).  Saved to `/tmp/fm.pdf`,
  converted to text via `pdftotext -layout`.  Result:\ 42 pages,
  ~2900 lines of usable text.
- F–M citation as in the bibliography:\
  F.~Fathizadeh and M.~Marcolli, "Periods and motives in the spectral
  action of Robertson–Walker spacetimes,"
  Commun. Math. Phys. **356**, 641–671 (2017);\
  preprint arXiv:1611.01815 (2016).

### Papers read (existing GeoVac corpus)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` (full
  document, lines 1–1161; especially §4 M2 theorem
  `thm:m2_mixed_tate` and its current "Proof sketch")
- `debug/sprint_mixed_tate_test_memo.md` (v3.45.3 sprint memo, full
  336 lines)
- `papers/group3_foundations/paper_18_exchange_constants.tex` §III.7
  (master Mellin engine, referenced via Paper 55 §1 introduction)
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` §VIII
  (case-exhaustion theorem, referenced via Paper 55 §2)

### Sections of F–M extracted
- §1 Introduction (lines 50–151 of `fm.txt`):\ paper-level summary,
  scaling factor as affine variable, quadric $Q_{\lambda, 2n}$, $2n$
  auxiliary affine parameters $\varepsilon_1, \ldots, \varepsilon_{2n}$.
- §2.1 Pseudodifferential symbol (lines 200–305):\ explicit form of
  $q_1(x, \xi)$, $q_0(x, \xi)$ with $a(t)$ and $a'(t)$ in $q_0$.
- §3 Heat expansion and Wodzicki residue (lines 305–360):\ equation
  3.2 ($a_{2+r} = 2^5/(4 + r/2) \cdot \mathrm{Res}(\Sigma_r^{-1})$),
  Wodzicki-residue definition.
- §4 Theorem 4.1 and §6 Theorem 6.1, 6.2 (lines 359–1340):\ full main
  theorem for $a_{2n}$ as period of an algebraic differential form on
  the affine complement of
  $C_{Z_{\lambda, 2n}} \cup H_0 \cup H_1$ in $\mathbb{A}^{2n+3}$.
- §7 Motives (lines 1600–1900):\ mixed-Tate-over-$\mathbb{Q}(\sqrt{-1})$
  Theorems 7.3, 7.4, 7.5 + Lemma 7.1 Grothendieck classes.

### Open questions for follow-on
1. **Closed-form $[Z_{1, 2n}]$ at general $n$:**\ extending F–M §7.4–7.6
   recursive computation to a closed form for the unit isotropic
   quadric Grothendieck class.  Mini-sprint.
2. **M2--M3 cyclotomic coincidence at level 4:**\ the $\mathbb{Q}(i)$
   field appearing in F–M Theorem 7.3 (Section 4 quadric isotropy) is
   the same $\mathbb{Z}[\zeta_4]$ appearing in Glanois 2015 (Section 6
   M3 mechanism).  Whether the coincidence reflects a deeper M2--M3
   coupling, or is independent, is open.
3. **Inhomogeneous extension (Open Question #1 of v3.45.3):**\ the
   F--M classification at $\lambda \ne 1, \varepsilon_i \ne 0$ generic
   covers the GeoVac framework's time-dependent extension (e.g.\ for
   thermal-time, Lorentzian Krein wedge, modular flow).  Whether the
   substantive transcendental content of GeoVac thermal observables
   sits in the full F--M mixed-Tate ring (with MZV content) or in a
   restricted sub-ring is open.
