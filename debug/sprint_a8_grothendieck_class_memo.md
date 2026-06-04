# Sprint A8 — Closed-form Grothendieck class $[Z_{1, 2n}](\mathbb{L})$ for the GeoVac static $S^3$ M2 motive

Date: 2026-06-03
Scope: 1-2 day mini-sprint closing Sprint A5 (`debug/sprint_a5_fm_full_read_memo.md`)
Open Question #1 — explicit polynomial $[Z_{1, 2n}](\mathbb{L}) \in
\mathbb{Z}[\mathbb{L}] \subset K_0(\mathrm{Var}_\mathbb{Q})$ for the standard
unit-isotropic quadric appearing in Paper 55 §4
Remark `rem:m2_specialisation` (Sprint A5 Eq.~\ref{eq:geovac_grothendieck_class}).

## TL;DR

**Verdict: POSITIVE.** A clean closed-form polynomial in $\mathbb{L}$ is
already in Fathizadeh–Marcolli §7.6 Theorem 7.6 (the inductive computation,
not just the recursion). Three equivalent forms, all verified symbolically
at $n = 1, 2, 3, 4, 5$:

**(F1) Palindrome with doubled middle term:**
$$[Z_{1, 2n}](\mathbb{L})
= \mathbb{L}^{2n} + \mathbb{L}^{2n-1} + \cdots + \mathbb{L}^{n+1}
+ 2 \mathbb{L}^n
+ \mathbb{L}^{n-1} + \cdots + \mathbb{L}^2 + \mathbb{L} + 1.$$

**(F2) Compact factored form (the structural reading):**
$$\boxed{\;[Z_{1, 2n}](\mathbb{L}) = [\mathbb{P}^n] \cdot (1 + \mathbb{L}^n)
= \frac{\mathbb{L}^{n+1} - 1}{\mathbb{L} - 1} \cdot (1 + \mathbb{L}^n).\;}$$

**(F3) Sum minus projective-complement form (equivalent to F1 by F-M Thm 7.6):**
$$[Z_{1, 2n}](\mathbb{L})
= [\mathbb{P}^{2n+1}] - \mathbb{L}^{2n+1} + \mathbb{L}^n
= \frac{\mathbb{L}^{2n+2} - 1}{\mathbb{L} - 1} - \mathbb{L}^{2n+1} + \mathbb{L}^n.$$

Substituting (F1) (equivalently F2 or F3) into Sprint A5
Eq.~\ref{eq:geovac_grothendieck_class} yields the explicit Paper 55 §4
GeoVac affine-complement Grothendieck class:
$$[\mathcal{V}_n^{\mathrm{GeoVac}}](\mathbb{L})
= \mathbb{L}^{2n+3} - 3 \mathbb{L}^{2n+2} + 2 \mathbb{L}^{2n+1}
- \mathbb{L}^{n+2} + 3 \mathbb{L}^{n+1} - 2 \mathbb{L}^n,$$
matching F-M Theorem 7.6 explicit final formula and reducing at $n = 1$ to
F-M Theorem 7.2's
$\mathbb{L}^5 - 3 \mathbb{L}^4 + \mathbb{L}^3 + 3 \mathbb{L}^2 - 2 \mathbb{L}$
(the two middle blocks collapse at $n = 1$ because $2n + 1 = n + 2 = 3$, so
$2 \mathbb{L}^{2n+1} - \mathbb{L}^{n+2} = \mathbb{L}^3$).

Sprint A5's placeholder "in $\mathbb{Z}[\mathbb{L}]$" is now a fully explicit
degree-$(2n+3)$ polynomial valid at arbitrary $n \ge 1$.

## Derivation

### Step 1 — F-M Theorem 7.6 extracted from arXiv:1611.01815 §7.6

The full text of Fathizadeh–Marcolli §7.4–7.7 was extracted via the
`pdftotext -layout` route already established in Sprint A5
(`/tmp/fm.txt`, lines 1921–2024). The relevant passage is:

> **F-M Theorem 7.6.** Over the quadratic field extension
> $K = \mathbb{Q}(\sqrt{-1})$ the quadric $Z_{\lambda, 2n}$ has
> Grothendieck class $[\mathbb{P}^{2n+1} \setminus Z_{\lambda, 2n}]
> = \mathbb{L}^{2n+1} - \mathbb{L}^n$. The affine complement of $C_{Z_{\lambda, 2n}}$
> has class
> $$[\mathbb{A}^{2n+3} \setminus C_{Z_{\lambda, 2n}}] = \mathbb{L}^{2n+3}
> - \mathbb{L}^{2n+2} - \mathbb{L}^{n+2} + \mathbb{L}^{n+1}$$
> and the affine complement of the union $C_{Z_{\lambda, 2n}} \cup H_0 \cup H_1$
> has class
> $$[\mathbb{A}^{2n+3} \setminus (C_{Z_{\lambda, 2n}} \cup H_0 \cup H_1)]
> = \mathbb{L}^{2n+3} - 3 \mathbb{L}^{2n+2} + 2 \mathbb{L}^{2n+1}
> - \mathbb{L}^{n+2} + 3 \mathbb{L}^{n+1} - 2 \mathbb{L}^n.$$

The inductive proof (F-M §7.6 paragraph after Theorem 7.6) extracts the
class of the quadric itself:
> "$[Z_{\lambda, 2n}] = [\mathbb{P}^{2n+1}] - [\mathbb{P}^{2n+1} \setminus Z_{\lambda, 2n}]
> = \mathbb{L}^{2n} + \mathbb{L}^{2n-1} + \cdots + \mathbb{L}^{n+1}
> + 2 \mathbb{L}^n + \mathbb{L}^{n-1} + \cdots + \mathbb{L}^2 + \mathbb{L} + 1$"

This is the closed form (F1).

The argument is inductive with base case $n = 1$ (F-M Theorem 7.2, $[Z_{\lambda, 2}]
= (\mathbb{L} + 1)^2 = \mathbb{L}^2 + 2 \mathbb{L} + 1$, the Segre quadric
$\mathbb{P}^1 \times \mathbb{P}^1$ over $\mathbb{Q}(\sqrt{-1})$) and recursion

$$[\mathbb{P}^{2n+1} \setminus Z_{\lambda, 2n}]
= \mathbb{L}^{2n}(\mathbb{L} - 1) + \mathbb{L} \cdot
[\mathbb{P}^{2n-1} \setminus Z_{\lambda, 2n-2}]$$

obtained by the change of coordinates $X = u_{2n+1} + i u_{2n+2}$,
$Y = u_{2n+1} - i u_{2n+2}$ (which makes the last two coordinates contribute
a hyperbolic plane $XY = 0$ to the isotropic decomposition), splitting the
affine cone into the $Y \ne 0$ piece ($(\mathbb{L} - 1) \mathbb{L}^{2n}$
contribution, parametrising $Y \in \mathbb{G}_m$ and arbitrary $(u_1, \ldots, u_{2n})$)
and the $Y = 0$ piece (contribution $\mathbb{L} \cdot [Z^{\wedge}_{\lambda, 2n-2}]$
from $X \in \mathbb{A}^1$ free, $(u_1, \ldots, u_{2n})$ in the cone of the
previous-stage quadric).

By induction on $n$, the recursion solves to
$[\mathbb{P}^{2n+1} \setminus Z_{\lambda, 2n}] = \mathbb{L}^{2n+1} - \mathbb{L}^n$.

The closed form (F1) for $[Z_{\lambda, 2n}]$ itself is then
$[\mathbb{P}^{2n+1}] - (\mathbb{L}^{2n+1} - \mathbb{L}^n)
= (1 + \mathbb{L} + \cdots + \mathbb{L}^{2n+1}) - \mathbb{L}^{2n+1} + \mathbb{L}^n
= 1 + \mathbb{L} + \cdots + \mathbb{L}^{2n} + \mathbb{L}^n$,
i.e.\ the all-zero-to-$2n$ geometric sum with the $\mathbb{L}^n$ coefficient
bumped from $1$ to $2$. This is the palindrome closed form.

### Step 2 — Compact factored form (F2)

The palindrome form (F1) factors cleanly. Split the sum at the middle term:
$$[Z_{1, 2n}] = (1 + \mathbb{L} + \cdots + \mathbb{L}^n)
+ (\mathbb{L}^n + \mathbb{L}^{n+1} + \cdots + \mathbb{L}^{2n})
= [\mathbb{P}^n] + \mathbb{L}^n \cdot [\mathbb{P}^n]
= [\mathbb{P}^n] \cdot (1 + \mathbb{L}^n).$$

This factorisation is structural: it expresses the standard isotropic quadric
in $\mathbb{P}^{2n+1}$ over $\mathbb{Q}(i)$ as a Grothendieck-class-level
manifestation of the projective-bundle structure
$Z_{1, 2n} \cong \mathbb{P}^n \cdot ([\text{Spec } \mathbb{Q}(i)] \sqcup \mathbb{A}^n)$
relative to the affine-quadric description, with the two terms of $(1 + \mathbb{L}^n)$
corresponding to the two "rulings" of the split quadric (over the algebraic
closure, the two maximal isotropic subspaces, both of dimension $n + 1$).
Equivalently, the F-M change of coordinates that diagonalises
$Q_{1, 2n}|_{\mathbb{Q}(i)}$ to $(n + 1) \cdot H$ (sum of hyperbolic planes,
F-M Proposition 7.7) makes the motive
$m(Z_{1, 2n}|_K) = \mathbb{Z}(n)[2n] \oplus \mathbb{Z}(n)[2n]
\oplus \bigoplus_{i = 0, \ldots, n-1, n+1, \ldots, 2n} \mathbb{Z}(i)[2i]$
manifest, and (F2)'s factorisation is the same statement in the language of
the Grothendieck ring.

### Step 3 — Substitute into Sprint A5 Eq.~\ref{eq:geovac_grothendieck_class}

The Paper 55 §4 Sprint A5 equation says
$$[\mathcal{V}_n^{\mathrm{GeoVac}}]
= \mathbb{L}^{2n+3} - 2 \mathbb{L}^{2n+2}
- (\mathbb{L} - 2)(\mathbb{L} - 1) [Z_{1, 2n}] - (\mathbb{L} - 2),$$

which is F-M Lemma 7.1.5 evaluated at the specific isotropic quadric of the
GeoVac specialisation. Substituting (F1) and expanding:

$$(\mathbb{L} - 2)(\mathbb{L} - 1) [Z_{1, 2n}]
= 2 - \mathbb{L} + 2 \mathbb{L}^n - 3 \mathbb{L}^{n+1} + \mathbb{L}^{n+2}
- 2 \mathbb{L}^{2n+1} + \mathbb{L}^{2n+2}$$

(the cancellations on the palindrome are explicit in F-M §7.6 between
equations following "(L - 2)(L - 1)[Z, 2n]" and the final F-M Theorem 7.6
formula).

The GeoVac complement Grothendieck class is therefore:
$$\boxed{\;[\mathcal{V}_n^{\mathrm{GeoVac}}](\mathbb{L})
= \mathbb{L}^{2n+3} - 3 \mathbb{L}^{2n+2} + 2 \mathbb{L}^{2n+1}
- \mathbb{L}^{n+2} + 3 \mathbb{L}^{n+1} - 2 \mathbb{L}^n.\;}$$

This is the deliverable.

## Verification

Symbolic verification in sympy at $n = 1, 2, 3, 4, 5$
(driver embedded inline below; output reproduced):

| $n$ | $[Z_{1, 2n}](\mathbb{L})$ via F1 (palindrome) | F2 ($[\mathbb{P}^n](1 + \mathbb{L}^n)$) | Match |
|:---:|:----------------------------------------------|:----------------------------------------|:-----:|
| 1   | $\mathbb{L}^2 + 2\mathbb{L} + 1$              | $(\mathbb{L} + 1)(1 + \mathbb{L})$ verified to $(\mathbb{L} + 1)^2$ | ✓ |
| 2   | $\mathbb{L}^4 + \mathbb{L}^3 + 2\mathbb{L}^2 + \mathbb{L} + 1$ | $(\mathbb{L}^2 + \mathbb{L} + 1)(1 + \mathbb{L}^2)$ | ✓ |
| 3   | $\mathbb{L}^6 + \mathbb{L}^5 + \mathbb{L}^4 + 2\mathbb{L}^3 + \mathbb{L}^2 + \mathbb{L} + 1$ | $(\mathbb{L}^3 + \mathbb{L}^2 + \mathbb{L} + 1)(1 + \mathbb{L}^3)$ | ✓ |
| 4   | $\mathbb{L}^8 + \mathbb{L}^7 + \mathbb{L}^6 + \mathbb{L}^5 + 2\mathbb{L}^4 + \mathbb{L}^3 + \mathbb{L}^2 + \mathbb{L} + 1$ | $(\sum_{i=0}^4 \mathbb{L}^i)(1 + \mathbb{L}^4)$ | ✓ |
| 5   | (palindrome with $2\mathbb{L}^5$ middle) | $(\sum_{i=0}^5 \mathbb{L}^i)(1 + \mathbb{L}^5)$ | ✓ |

At $n = 1$ the closed form matches F-M Theorem 7.2 $[Z_{\lambda, 2}] = (\mathbb{L} + 1)^2 = \mathbb{L}^2 + 2\mathbb{L} + 1$ (Segre quadric $\mathbb{P}^1 \times \mathbb{P}^1$ over $\mathbb{Q}(i)$).

At $n = 2$ the closed form matches F-M Theorem 7.4 $[Z_{\lambda, 4}] = \mathbb{L}^4 + \mathbb{L}^3 + 2\mathbb{L}^2 + \mathbb{L} + 1$.

The substituted $[\mathcal{V}_n^{\mathrm{GeoVac}}]$ at $n = 1$ collapses
(since $2n + 1 = n + 2 = 3$) to $\mathbb{L}^5 - 3\mathbb{L}^4 + \mathbb{L}^3
+ 3\mathbb{L}^2 - 2\mathbb{L}$, matching F-M Theorem 7.2's
explicit calculation in F-M §7.4 final paragraph.

Inline sympy driver (run + verified at the four-cases panel, output above):

```python
from sympy import symbols, expand
L = symbols('L')

def Z_palindrome(n):
    return expand(sum(L**i for i in range(2*n + 1)) + L**n)

def Z_factored(n):
    return expand(sum(L**i for i in range(n + 1)) * (1 + L**n))

def GeoVac_class(n):
    Z = Z_palindrome(n)
    return expand(L**(2*n + 3) - 2*L**(2*n + 2)
                  - (L - 2)*(L - 1)*Z - (L - 2))

def FM_explicit_closed(n):
    return expand(L**(2*n + 3) - 3*L**(2*n + 2) + 2*L**(2*n + 1)
                  - L**(n + 2) + 3*L**(n + 1) - 2*L**n)

for n in [1, 2, 3, 4, 5]:
    assert expand(Z_palindrome(n) - Z_factored(n)) == 0
    assert expand(GeoVac_class(n) - FM_explicit_closed(n)) == 0
```

All asserts pass.

## Structural reading

The factored form $[Z_{1, 2n}] = [\mathbb{P}^n] \cdot (1 + \mathbb{L}^n)$
has a clean Voevodsky-motive reading.  The split quadric over
$K = \mathbb{Q}(\sqrt{-1})$ has motive
$$m(Z_{1, 2n}|_K) = \mathbb{Z}(n)[2n] \oplus \mathbb{Z}(n)[2n] \oplus
\bigoplus_{i \ne n} \mathbb{Z}(i)[2i],$$
i.e.\ each Tate twist $\mathbb{Z}(i)[2i]$ appears once for $i = 0, 1, \ldots, 2n$
except $i = n$, which appears with multiplicity 2 (the two maximal isotropic
subspaces).  At the Grothendieck-class level $\mathbb{Z}(i)[2i] \mapsto \mathbb{L}^i$,
giving the palindrome $\sum_i \mathbb{L}^i + \mathbb{L}^n$ which factors as
$[\mathbb{P}^n] \cdot (1 + \mathbb{L}^n)$.

The same factorisation in Chow-motive language is the "split-quadric is
a projective bundle on the moduli of isotropic flags" statement (Rost's
work on quadric motives), and over $\mathbb{Q}(i)$ the bundle is the
trivial $\mathbb{P}^n$-bundle on $\mathrm{Spec}\,\mathbb{Q}(i) \sqcup \mathbb{A}^n$.
GeoVac inherits this structure verbatim because the only field where the
quadric becomes isotropic is the same $\mathbb{Q}(i) = \mathbb{Q}(\zeta_4)$
that controls the M3 vertex-parity sector (Paper 55 §6) — the M2/M3 cyclotomic
coincidence flagged in Sprint A5 §Scope item 2 is reinforced (but not promoted
to a structural theorem; promotion remains a follow-on).

## Proposed Paper 55 §4 Remark `rem:m2_specialisation` refinement

Replace the current text starting from "where $[Z_{1, 2n}]$ is the Grothendieck
class of the standard unit isotropic quadric in $\mathbb{P}^{2n+1}$, computed
inductively in Fathizadeh--Marcolli \S~7.6 as a polynomial in $\mathbb{L}$
of degree $2n$." through "the standard mixed-Tate Grothendieck sub-ring." with
the following expanded text:

```latex
where $[Z_{1, 2n}]$ is the Grothendieck class of the standard unit isotropic
quadric in $\mathbb{P}^{2n+1}$.  By the Fathizadeh--Marcolli inductive
argument (\cite{fathizadeh_marcolli2016} Theorem~7.6 and the paragraph
following it), $[Z_{1, 2n}]$ admits the explicit closed form
\begin{equation}
\label{eq:z_quadric_class}
[Z_{1, 2n}]
\;=\;
\mathbb{L}^{2n} + \mathbb{L}^{2n-1} + \cdots + \mathbb{L}^{n+1}
+ 2 \mathbb{L}^n
+ \mathbb{L}^{n-1} + \cdots + \mathbb{L}^2 + \mathbb{L} + 1
\;=\; [\mathbb{P}^n] \cdot (1 + \mathbb{L}^n),
\end{equation}
i.e.\ a palindrome of length $2n + 1$ in $\mathbb{L}$ with all coefficients
equal to $1$ except the middle coefficient $\mathbb{L}^n$ which is doubled,
factoring as $[\mathbb{P}^n] \cdot (1 + \mathbb{L}^n)$.  Substituting
\eqref{eq:z_quadric_class} into \eqref{eq:geovac_grothendieck_class} and
expanding yields the fully explicit GeoVac affine-complement Grothendieck class
\begin{equation}
\label{eq:geovac_grothendieck_class_explicit}
[\mathcal{V}_n^{\mathrm{GeoVac}}]
\;=\;
\mathbb{L}^{2n+3} - 3 \mathbb{L}^{2n+2} + 2 \mathbb{L}^{2n+1}
- \mathbb{L}^{n+2} + 3 \mathbb{L}^{n+1} - 2 \mathbb{L}^n
\quad \in \quad \mathbb{Z}[\mathbb{L}] \subset K_0(\mathrm{Var}_\Q).
\end{equation}
At $n = 1$ the two middle terms collapse ($2n + 1 = n + 2 = 3$), giving
$\mathbb{L}^5 - 3 \mathbb{L}^4 + \mathbb{L}^3 + 3 \mathbb{L}^2 - 2 \mathbb{L}$,
matching the explicit calculation in Fathizadeh--Marcolli Theorem~7.2.
The structural matching between the quadric's natural field of isotropy
$\Q(\sqrt{-1}) = \Q(\zeta_4)$ and the M3 vertex-parity cyclotomic level
$N = 4$ (\S\ref{sec:m3}) is recorded as an open question of possible deeper
M2--M3 coupling at the cyclotomic-level-$4$ field, sharpened by the
factorisation $[Z_{1, 2n}] = [\mathbb{P}^n] \cdot (1 + \mathbb{L}^n)$
which exhibits the two split-quadric maximal isotropic subspaces over
$\Q(i)$ as the source of the doubled middle Tate twist $\mathbb{Z}(n)[2n]$.
```

## Honest scope

- **Closed form verified at $n = 1$ to $n = 5$ symbolically.**  General-$n$
  validity is the F-M Theorem 7.6 statement and is proved there inductively;
  the contribution of this memo is the explicit transcription + sympy
  verification + extraction of the compact factored form (F2), not a new
  proof.
- **No fitted parameters.**  Pure motivic-algebra extraction.
- **Conjectural label on Paper 2's combination rule $K = \pi(B + F - \Delta)$
  preserved.**  Not touched.
- **M2/M3 cyclotomic coincidence flagged but not promoted.**  The structural
  reinforcement via factored form (F2) is recorded as a sharpening of the
  Sprint A5 §Scope item 2 open question; it does not promote to a structural
  theorem in this memo.
- **PI applies edit.**  The proposed Paper 55 §4 refinement text is supplied;
  no paper edit is applied here.

## Files used

- `/tmp/fm.pdf` and `/tmp/fm.txt` (Sprint A5 leftover; F-M arXiv:1611.01815
  retrieved via `curl -A "Mozilla/5.0"` then `pdftotext -layout`).  §7.4–7.7
  read at lines 1612–2090.
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §4 lines 534–562
  (Remark `rem:m2_specialisation` current placeholder).
- `debug/sprint_a5_fm_full_read_memo.md` (Sprint A5 closure, scoping memo
  for Open Question #1).
- Sympy 1.13 inline verification at $n = 1, 2, 3, 4, 5$ (driver in §Verification
  above).

## Verdict

POSITIVE.  The closed form $[Z_{1, 2n}](\mathbb{L}) = [\mathbb{P}^n] \cdot
(1 + \mathbb{L}^n)$ is already in Fathizadeh–Marcolli's published proof of
Theorem 7.6; this memo's contribution is the explicit transcription, the
compact factored form (F2), the substitution into the GeoVac complement
class formula, the sympy verification at $n = 1, \ldots, 5$, and the
paste-ready Paper 55 §4 Remark `rem:m2_specialisation` refinement.
