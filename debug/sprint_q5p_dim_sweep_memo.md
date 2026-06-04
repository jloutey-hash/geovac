# Sprint Q5' — Dimension sweep on M2 and M3 over $S^3$, $S^5$, $S^7$

Date: 2026-06-03
Scope: 1-day scoping; structural reasoning on published periods facts
already in Paper 55 bibliography; no new code, no PSLQ, no numerical
checks. Sprint A7's S^5 robustness audit (memo
`debug/sprint_a7_m2_m3_cyclotomic_memo.md` lines 192–212) is the
starting point; Sprints A2 (S^5 M2) and A6 (S^5 M3) supply the
explicit S^5 closed forms; Sprint Mixed-Tate Test (S^3 M2 baseline)
supplies the comparison.

## TL;DR

**Verdict: NEGATIVE for the deeper-unification reading of Q5'.**

The dimensional sweep across $d \in \{3, 5, 7\}$ shows:

- **M2's field-of-output ring is controlled by $\mathrm{Vol}(S^d) =
  2\pi^{(d+1)/2}/\Gamma((d+1)/2)$ and dimension parity.** The
  arithmetic content shifts in a *predictable, generic-Euclidean* way:
  $d = 3 \to \bigoplus_k \pi^{2k}\Q$ (even-weight pure-Tate sub-ring),
  $d = 5 \to \pi^3\Q$ (single odd-weight slice), $d = 7 \to
  \bigoplus_k \pi^{4 + 2k}\Q$ (even-weight pure-Tate sub-ring shifted by
  $\pi^4$ leading volume). The Witt-isotropy field $\Q(i)$ depends
  ONLY on the Euclidean signature of $\R^{d+1}$ (which is $(+, \ldots,
  +)$ for every $d$), so $\Q(i)$ appears for every $d$ at the
  field-of-proof level, and the period output ring lives entirely
  in $\Q[\pi]$ for every $d$. **There is no dimension-parity move
  inside the M2 number field**: $\Q(i)$ is generic Witt-splitting at
  every $d$, with the dimension only controlling the Tate weight of
  the slice.

- **M3's field-of-output ring is controlled by the Camporesi–Higuchi
  half-integer shift $(d-1)/2$, but the cyclotomic level stays at $4$
  for every $d$ at which CH has half-integer spectrum.** $d = 3$:\
  shift $3/2$, quarter-integer Hurwitz $1/4, 3/4 \in \mu_4$,
  conductor of $\chi_{-4}$. $d = 5$:\ shift $5/2$, quarter-integer
  Hurwitz $5/4, 7/4$ reducing to $1/4, 3/4$ via Hurwitz shift
  $\zeta(s, a+1) = \zeta(s, a) - a^{-s}$ (Sprint A6 Cor.\
  `cor:m3_s5_cyclotomic`). $d = 7$:\ shift $7/2$, quarter-integer
  Hurwitz $7/4, 9/4$ reducing to $1/4, 3/4$ by the same shift
  identity. **No new cyclotomic level is engaged for any odd $d$ at
  which CH has a half-integer shift**; the cyclotomic field
  $\Q(\zeta_4)$ is **structurally tied to the parity of the
  half-integer shift modulo 2**, which is a fixed Z/2 input.

- **No shared dimension-dependent structure emerges.** The M2 ring
  shifts (in $\pi$ weight), but the underlying number field $\Q(i)$ is
  fixed as generic Euclidean Witt-splitting. The M3 ring shifts (in
  the number of polynomial terms $\zeta(s-2k)$, $k = 0, \ldots,
  (d-1)/2$ at integer $s$), but the underlying cyclotomic field
  $\Q(\zeta_4)$ is fixed by the parity of the CH shift's denominator.
  Both fields stay at $\Q(i)$ for **structurally independent reasons
  that are not strengthened by the sweep**. The "shared $\Q(i)$"
  observation persists across $d = 3, 5, 7$ for the **same structurally
  independent reasons** as at $d = 3$ — exactly Sprint A7's
  HALF-STRUCTURAL reading, with no dimension-parity escape hatch.

This **weakens Q5'**:\ if there were a deeper unification it would
plausibly manifest as a shared dimension-dependent mirror (e.g., the
M3 cyclotomic level rising to $8$ exactly where M2's slice shifts
parity, or M2's Witt field jumping to $\Q(\zeta_8)$ exactly where
M3 picks up an extra Hurwitz term). Nothing of the kind happens.
The two fields agree at $\Q(i)$ for one structural reason each, and
neither reason cares about the other.

The HALF-STRUCTURAL verdict from Sprint A7 is preserved AT A
TIGHTER level:\ the dimension sweep confirms the field-coincidence
is **field-of-input generic-Euclidean for M2** and **field-of-output
conductor-of-$\chi_{-4}$ for M3** in a way that does NOT change
with $d$.

## The four structural facts the sweep rests on

### Fact 1: M2's Witt-isotropy field is generic Euclidean for all $d$

Fathizadeh–Marcolli (arXiv:1611.01815) classify the SD coefficients of
the CC spectral action on a R–W spacetime via the quadric $Q_{\lambda,
2n} = \sum_{j=1}^{2n+2} u_j^2 - \lambda$ ($\lambda$ encodes the
scaling factor $a(t)$). On the GeoVac static specialisation $\lambda =
1$:
$$
Q_{1, 2n} = u_1^2 + u_2^2 + \ldots + u_{2n+2}^2.
$$
This form is **anisotropic over $\R$ and $\Q$, isotropic over $\Q(i)$**
via the elementary Witt-pair construction $X_j = u_{2j-1} + i u_{2j}$,
$Y_j = u_{2j-1} - i u_{2j}$, with $X_j Y_j = u_{2j-1}^2 + u_{2j}^2$.

This is dimension-INDEPENDENT in a strong sense:\ the rank $2n + 2$
needs to be even for the hyperbolic pairs to consume everything, and
the Witt-isotropy index over $\Q(i)$ reaches $n + 1$ regardless of $n$.
The argument transfers to $S^5, S^7, \ldots$ with no modification:\ on
$S^d$ the embedding $S^d \hookrightarrow \R^{d+1}$ has rank $d + 1$, and
positive-definite Euclidean signature is $(+, +, \ldots, +)$ for every
$d$. Witt isotropy index over $\Q(i)$ is then $\lfloor (d+1)/2 \rfloor$
(maximal for $d + 1$ even, one short of maximal for $d + 1$ odd, but
the residual anisotropic line lives over $\R$ and disappears under
volume normalisation anyway).

**The field $\Q(i)$ is therefore selected at every $d$ by signature
alone**, with no role for the parity of $d$ at the level of which
field appears. The dimension only controls the volume normalisation
factor $(4\pi)^{d/2}$, which fixes the Tate WEIGHT of the output
slice, not its number field.

### Fact 2: M2's output ring is dimension-parity controlled, but inside $\Q[\pi]$

Sprint A2 (`debug/sprint_a2_s5_mixed_tate_memo.md`) established the
dimension-parity sharpening at $d = 5$:

| $d$ | $\mathrm{Vol}(S^d)$ | M2 output ring | $S^d$ slice |
|:---:|:--------------------|:---------------|:------------|
| 3 | $2\pi^2$ | $\bigoplus_k \pi^{2k}\Q$ | even-weight pure-Tate sub-ring |
| 5 | $\pi^3$ | $\pi^3 \cdot \Q$ | single odd-weight Tate slice |
| 7 | $\pi^4 / 3$ | $\bigoplus_k \pi^{4 + 2k}\Q$ | even-weight pure-Tate sub-ring shifted by $\pi^4$ |

Explanation:\ on $S^d$ with $d$ odd, $\mathrm{Vol}(S^d) =
2\pi^{(d+1)/2}/\Gamma((d+1)/2)$. For $d = 3$ this is $2\pi^2$, for
$d = 5$ this is $\pi^3$, for $d = 7$ this is $\pi^4/3$.

The Mellin transform of $K_{D^2}^{(d)}(t) \sim \sqrt\pi \cdot t^{-d/2}$
(the leading singularity of the heat trace at odd $d$ on the CH spinor
bundle) produces a $\sqrt\pi$ prefactor in the raw SD coefficients
$\tilde a_k$; the volume normalisation $(4\pi)^{d/2}$ absorbs $\sqrt\pi$
exactly when $d$ is odd, leaving an even/odd-weight power of $\pi$ in
the volume-normalised $a_k$:
$$
a_k = (4\pi)^{d/2} \tilde a_k = (4\pi)^{d/2} \cdot \sqrt\pi \cdot q_k
    = (4\pi)^{(d-1)/2} \cdot 4\pi \cdot q_k \cdot \frac{1}{\sqrt{4\pi}} \cdot \sqrt\pi
    = \pi^{(d+1)/2} \cdot \text{(rational)}.
$$
Then the Bernoulli mechanism that gives the discrete CH Dirac
two-term-exact SD on $S^3$ (Sprint Mixed-Tate Test) gives, on $S^{2m+1}$,
$(m+1)$-term-exact SD (Paper 51 `rem:two_term_uniqueness`), with each
non-zero $a_k^{D^2}$ in $\pi^{(d+1)/2} \cdot \Q$ exactly. So:

- $d = 3$ Dirac: 2 non-zero terms in $\pi^2 \Q$.
- $d = 5$ Dirac: 3 non-zero terms in $\pi^3 \Q$ (Sprint A2 closed form).
- $d = 7$ Dirac: **4 non-zero terms expected in $\pi^4 \Q$**
  (by extrapolation of Paper 51 `rem:two_term_uniqueness`; not
  explicitly computed in this sprint but the structural argument is
  identical).

For the scalar Laplacian, the $e^{((d-1)/2)^2 t}$-type prefactor produces
an infinite series in the corresponding $\pi$-power slice:\ $S^3$ gives
the geometric collapse $a_k^\Delta = 2\pi^2/k!$, $S^5$ gives the
3-term-numerator closed form $a_k^\Delta = (6-k) 4^{k-1} \cdot 2/(3 k!)
\cdot \pi^3$ (Sprint A2), and $S^7$ would give a 5-term-numerator closed
form (4 = $(d-1)/2$ terms from the quartic-vs-quadratic structure of
the degeneracy polynomial, plus the $e^{((d-1)/2)^2 t}$ leading factor).

**The pattern is**: the underlying number field is fixed at $\Q$
(after volume normalisation, no $\sqrt\pi$ survives), with the Tate
weight of the slice given by $(d+1)/2$. The dimension parity controls
**which weight slice** the SD lives in, but **does NOT engage a new
number field** — $\Q(i)$ enters at the *proof* level via Witt isotropy
on $\R^{d+1}$ for every $d$, and the *output* descends to $\Q$ via
Galois invariance for every $d$.

### Fact 3: M3's cyclotomic level is fixed at $\le 4$ for all CH-shift parities

The Camporesi–Higuchi half-integer shift on $S^d$ for odd $d$ is
$\frac{d-1}{2}$. The dimension sweep gives:

| $d$ | CH shift | Hurwitz quarter-shifts in $D_{\rm even/odd}$ | Reduction to $\{1/4, 3/4\}$ | Cyclotomic level needed |
|:---:|:---------|:----------------------------------------------|:----------------------------|:-----------------------:|
| 3 | $3/2$ | $1/4, 3/4$ (directly) | trivial | $N = 4$ |
| 5 | $5/2$ | $5/4, 7/4$ | $\zeta(s, 5/4) = \zeta(s, 1/4) - 4^s$, $\zeta(s, 7/4) = \zeta(s, 3/4) - (4/3)^s$ | $N = 4$ |
| 7 | $7/2$ | $7/4, 9/4$ | $\zeta(s, 7/4) = \zeta(s, 3/4) - (4/3)^s$, $\zeta(s, 9/4) = \zeta(s, 1/4) - 4^s - (4/5)^s$ | $N = 4$ |

The argument (Sprint A6 Cor.\ `cor:m3_s5_cyclotomic`, transferred
verbatim):\ the parity decomposition $D = D_{\rm even} + D_{\rm odd}$
on the CH spectrum with shift $(d-1)/2$ produces sums of the form
$\sum_{k \ge 0} (2k + (d-1)/2 + \frac{1}{2})^{-s}$ and $\sum_{k \ge 0}
(2k + (d-1)/2 + \frac{3}{2})^{-s}$. Rescaling $u = m/2$:
$$
Z_{\rm even}^{(d)}(s) = 2^{-s} \zeta(s, a_d^{\rm even}), \qquad
Z_{\rm odd}^{(d)}(s) = 2^{-s} \zeta(s, a_d^{\rm odd}),
$$
with
$$
a_d^{\rm even} = \frac{1}{4} \cdot (d - 1) + \frac{1}{4} \cdot 1 = \frac{d}{4},
\qquad
a_d^{\rm odd} = \frac{1}{4} \cdot (d - 1) + \frac{1}{4} \cdot 3 = \frac{d + 2}{4}.
$$
At $d = 3$:\ $a^{\rm even} = 3/4, a^{\rm odd} = 5/4 = 1 + 1/4$.
At $d = 5$:\ $a^{\rm even} = 5/4 = 1 + 1/4, a^{\rm odd} = 7/4 = 1 + 3/4$.
At $d = 7$:\ $a^{\rm even} = 7/4 = 1 + 3/4, a^{\rm odd} = 9/4 = 2 + 1/4$.

In every case the fractional parts are in $\{1/4, 3/4\}$ — the
**level-4 quarter-integers** — with the integer parts adding only
$\Q$-rational shift corrections via the Hurwitz shift identity
$\zeta(s, a+1) = \zeta(s, a) - a^{-s}$.

**Consequently the cyclotomic level stays at $N = 4$ for every odd
$d \ge 3$.** The $\Q$-rational shift corrections add finitely many
$a^{-s}$ terms with $a$ a positive rational; these do not engage any
new cyclotomic field. The Dirichlet beta $\beta(s) = L(s, \chi_{-4})$
emerges from $\zeta(s, 1/4) - \zeta(s, 3/4) = 4^s \beta(s)$ at every
$d$.

**There is no dimension-parity move inside the M3 number field
either.** The Galois descent from level 4 to level 2 via $\chi_{-4}$
(Deligne–Glanois) is structurally the same at every $d$; the dimension
only controls the **number of polynomial terms** in the rational
combination, which is $(d-1)/2 + 1 = (d+1)/2$ terms reflecting the
degree of the CH degeneracy polynomial:

| $d$ | Degeneracy polynomial degree | $D^{(d)}(s)$ polynomial terms in $Z(s - 2k)$ |
|:---:|:----------------------------:|:--------------------------------------------:|
| 3 | 2 (quadratic) | 2 |
| 5 | 4 (quartic) | 3 |
| 7 | 6 (sextic) | 4 |

But every term lives in $\Q[\zeta(s), \zeta(s-2), \ldots]$ at the
un-restricted level, and in $\Q[f_d(s), f_d(s-2), \ldots]$ at the
vertex-restricted level with $f_d(s) = 2^s \beta(s) - 2^s + \text{(rat.\
corrections)}$. The number field is fixed at $\Q(\zeta_4)$ for every
odd $d \ge 3$.

### Fact 4: $S^7$ has no GeoVac-internal home

GeoVac has explicit constructions on $S^3$ (Fock-projected Coulomb,
Paper 7) and $S^5$ (Bargmann–Segal HO Hardy lattice, Paper 24).
$S^7$ has **no internal GeoVac home**:

- No published packing/projection construction would place a physical
  system on $S^7$;\ the closed-orbit (Bertrand) + Fock-rigidity arguments
  pick $S^3$ for $-Z/r$ and the holomorphic Hardy sector of $S^5$ for
  the 3D HO, with rigidity theorems (Paper 24 Thm 3) showing no other
  sphere is conformally compatible with either of those two operator
  types.

- Cross-manifold tensor products $\sthree \otimes \mathrm{Hardy}(\sfive)$
  are blocked at the cross-modular-Hamiltonian level by the four-layer
  Coulomb/HO asymmetry of Paper 24 §V (G3 frontier, Paper 32 §VIII.C);
  even if a $S^7$-internal construction existed it would not connect to
  the GeoVac-physical-substrate without crossing the asymmetry.

- The Hopf-tower analysis (Sprint Read 2 scoping, 2026-06-03)
  truncated the gauge-saturation at $S^5$ (complex Hopf $S^{2n-1} \to
  \SU(n)$ at $n \le 3$); the $S^7$ quaternionic Hopf $S^7 \to S^4$
  would engage $\SU(2)$ on the base and is not part of the SM
  truncation.

**Treatment in this sweep**:\ $S^7$ is included as a **generic-Euclidean
test** — what would the M2 / M3 patterns predict if the framework had
a $S^7$-internal home? — and the answer is that both rings would extend
in the manifestly generic way described above (M2 in $\pi^4 \Q$ even
weight, M3 still at cyclotomic level 4 with 4-term polynomial
structure). This **strengthens the NEGATIVE** verdict:\ the pattern is
not specific to the physically-realised $S^3, S^5$ dimensions; it's
structurally generic.

## The per-dimension structural table

Putting Facts 1–3 together, with $S^7$ included as the generic test:

| $d$ | Manifold | GeoVac home | M2 field-of-output | M2 Tate weight | M2 number-field signature | M3 CH shift | M3 cycl. level | M3 polynomial terms | M3 number-field signature |
|:---:|:--------:|:-----------:|:-------------------|:--------------:|:--------------------------|:------------|:---------------|:--------------------|:--------------------------|
| 3 | $S^3$ | Fock $-Z/r$ (Paper 7) | $\bigoplus_k \pi^{2k}\Q$ | even | $\Q(i)$ Witt, trivial Galois on output | $3/2$ | $N = 4$ | 2 in $\zeta(s-2), \zeta(s)$ | $\Q(\zeta_4)$ conductor, non-trivial Galois |
| 5 | $S^5$ | Bargmann–Segal HO Hardy (Paper 24) | $\pi^3 \Q$ | odd | $\Q(i)$ Witt, trivial Galois on output | $5/2$ | $N = 4$ (after Hurwitz shift) | 3 in $\zeta(s-4), \zeta(s-2), \zeta(s)$ | $\Q(\zeta_4)$ conductor, non-trivial Galois |
| 7 | $S^7$ | None (no GeoVac home) | $\bigoplus_k \pi^{4 + 2k}\Q$ | even | $\Q(i)$ Witt, trivial Galois on output | $7/2$ | $N = 4$ (after Hurwitz shift) | 4 in $\zeta(s-6), \zeta(s-4), \zeta(s-2), \zeta(s)$ | $\Q(\zeta_4)$ conductor, non-trivial Galois |

**The four columns "field-of-output", "Tate weight", "polynomial
terms", and "cyclotomic level" each follow a predictable pattern in
$d$:**

- M2 Tate weight scales linearly with $d$:\ $(d+1)/2$.
- M2 even-vs-odd-weight alternates with $d \pmod{4}$ (3 → even-weight
  ring, 5 → odd-weight slice, 7 → even-weight ring shifted by $\pi^4$).
- M3 polynomial-term count scales as $(d+1)/2$.
- M2 number field stays at $\Q(i)$ via signature.
- M3 number field stays at $\Q(\zeta_4)$ via conductor of $\chi_{-4}$.

The first two are M2-internal patterns; the third is an M3-internal
pattern; the last two are the field-coincidence we are testing. **The
last two do NOT mirror each other in any non-trivial way across the
sweep.** $\Q(i)$ appears in M2 for the same generic reason at every
$d$; $\Q(\zeta_4)$ appears in M3 for the same conductor reason at
every $d$.

## Camporesi–Higuchi shift $(d-1)/2$ — what this tells us

The half-integer shift $(d-1)/2$ feeds M3 via:

1. The denominator of $(d-1)/2 \pm 1/2$ being $2$, so after rescaling
   $u = m/2$ the quarter-integer Hurwitz shifts appear at fractional
   parts in $\{1/4, 3/4\}$ regardless of $d$.

2. The fractional part of $(d-1)/2$ being $0$ or $1/2$ depending on
   $d \pmod 4$ — specifically $(d-1)/2 \in \Z$ for $d \equiv 1
   \pmod 2$ (always, since CH applies on odd-dim spheres) and
   $(d-1)/2 + 1/2 = d/2$ is a half-integer for every odd $d$.

This is why $\Q(\zeta_4)$ stays the cyclotomic conductor for every
$d$:\ the shift modulo 2 is in $\{1/2\}$ for every odd $d$ (since
$(d-1)/2 + 1/2 = d/2$ and $d$ odd gives $d/2$ half-integer), and the
parity of the half-integer denominator is what selects $\chi_{-4}$.
A shift like $1/6$ or $1/3$ would engage level $12$ or $3$
respectively (different cyclotomic conductor); a shift like $1/2$ is
exactly the conductor-4 case.

**The structural read**:\ M3's $\Q(\zeta_4)$ is selected by the
**denominator of the CH shift modulo 2**, which is a *one-bit invariant*
of the substrate. It cannot encode any richer dimension-dependent
information, because there is only one bit to encode. M2's $\Q(i)$
is selected by the **Euclidean signature**, which is also a *one-bit
invariant* (signature parity for definite signature). These are two
DIFFERENT one-bit invariants — neither cares about the other.

The Sprint A7 sentence

> "Both occurrences trace to the dimension-3 parity of $S^3$, but via
> two structurally independent mechanisms"

is sharpened by the sweep to:

> Both occurrences are controlled by one-bit invariants of the
> substrate (Euclidean signature parity for M2, CH-shift denominator
> parity for M3); both one-bit invariants happen to evaluate to "even"
> on every odd-dimensional sphere, which is the structural reason the
> shared $\Q(i)$ persists across the dimension sweep. The two
> one-bit invariants are independent (no joint structure connects
> them); the agreement at $\Q(i)$ is a consequence of independent
> invariants happening to take the same value at all the tested
> substrates.

## Euclidean signature $(+, \ldots, +)$ on $\R^{d+1} \supset S^d$

For every odd $d$, the embedding $S^d \hookrightarrow \R^{d+1}$ has
positive-definite signature $(+, +, \ldots, +)$ with $d + 1$
positive entries. The Witt isotropy of the standard form
$\sum_{j=1}^{d+1} u_j^2$ over $\Q(i)$ is $\lfloor (d+1)/2 \rfloor$, with
the residual anisotropic line (when $d + 1$ is odd) living over $\R$
and disappearing under volume normalisation.

Lowest-dimensional rank-1 case ($d = 1, S^1$):\ the quadric $u_1^2 +
u_2^2 = 0$ has Witt index 1 over $\Q(i)$ (one hyperbolic pair, maximal).
Generic case ($d = 3$):\ Witt index 2. ($d = 5$):\ Witt index 3. ($d =
7$):\ Witt index 4. The Witt isotropy is **strictly monotone in $d$**;
the field $\Q(i)$ is FIXED.

This is the formal statement of the "$\Q(i)$ is generic Euclidean" for
the M2 substrate:\ the field of splitting is determined by signature
alone, with $\Q(i)$ being the canonical splitting field for any
positive-definite even-rank rational quadratic form.

## What a POSITIVE outcome would have looked like

A POSITIVE outcome of the dimension sweep would have required at least
one of the following:

1. **M2's number field shifts with $d \pmod 4$** — e.g., $\Q(i)$ at
   $d = 3$ and $d = 7$ but $\Q(\zeta_8)$ at $d = 5$ (the odd-weight
   slice). This would mirror M3's level-4-vs-level-8 distinction and
   suggest a shared dimensional Tannakian symmetry. **Does NOT happen**:
   M2's $\Q(i)$ is generic-signature, dimension-independent.

2. **M3's cyclotomic level shifts with $d$** — e.g., level $4$ at
   $d = 3$ but level $8$ at $d = 5$ or $d = 7$. **Does NOT happen**:
   Sprint A6 Cor.\ `cor:m3_s5_cyclotomic` proves level 4 is sufficient
   at $d = 5$; the same proof works at $d = 7$ (via the Hurwitz shift
   identity, no escalation).

3. **A dimension-dependent natural isomorphism between M2 and M3's
   $\Z/2$ Galois actions** — e.g., $\Q(i)$'s Galois acting trivially
   on M2 outputs but conjugating M3 outputs into themselves with the
   conjugation realising a known $S^d$ automorphism. **Does NOT
   happen**: M2's $\Z/2$ acts trivially on M2 outputs at every $d$
   (the output ring is $\Q[\pi^{\text{Tate-weight}}]$ which is
   $\Q$-rational); M3's $\Z/2$ acts non-trivially at every $d$ via the
   conductor-4 character. The two actions are independent of each
   other AND of $d$.

4. **A natural transformation between M2's pure-Tate sub-ring and M3's
   cyclotomic-mixed-Tate ring that respects the master Mellin slot
   index $k$ and depends on $d$** — would suggest a $d$-dependent
   functorial bridge. **Does NOT happen**: M2 sits at $k = 2$
   (second-order $D^2$) in $\Mcal[\Tr(D^k e^{-tD^2})]$, M3 sits at $k =
   1$ (first-order $D$ with vertex restriction), with no natural
   transformation in either direction at any $d$. Sprint MR-A
   (2026-05-06 evening) ruled out the $k = 0 \leftrightarrow k = 1
   \leftrightarrow k = 2$ engine-domain-mixing in the propinquity-rate
   direction (CLAUDE.md §1.7 WH1 MR-A entry); the same engine-domain
   partition prevents M2 ↔ M3 transformation at the period-ring level.

None of (1)–(4) is realised by the dimension sweep. The structural
NEGATIVE is robust.

## Honest scope

1. **No new computations**:\ this is pure motivic bookkeeping built on
   Sprint A2 (S^5 M2 closed forms), Sprint A6 (S^5 M3 closed forms),
   Sprint Mixed-Tate Test (S^3 M2 baseline), Sprint M3 cyclotomic mixed-Tate
   (S^3 M3 baseline), and Sprint A7's HALF-STRUCTURAL verdict.

2. **$S^7$ is included as a generic-Euclidean test**, with no claim that
   GeoVac has an internal $S^7$ home. The $S^7$ pattern follows from
   the F–M classification and Sprint A2/A6 structural arguments by
   direct extrapolation;\ no $S^7$-specific computation was performed
   in this sprint.

3. **Bertrand × Fock-rigidity caveat**:\ even if a $S^7$ extrapolation
   showed a deeper unification, the framework lacks a $S^7$-internal
   route to read it out as physics. The Paper 24 four-layer Coulomb/HO
   asymmetry is the structural obstruction beyond which a $S^7$
   classification would need a totally new construction.

4. **Q5' is not closed**:\ this sprint provides a NEGATIVE verdict on
   the **dimension-sweep entry point (i)** of the Q5' scoping. Other
   entry points remain possible — e.g., a depth-vs-loop-order sweep
   (sub-sprint entry point (ii)), or a higher-cyclotomic-level
   refinement under $\Q(\zeta_8)$ (the Bertrand $S^3$ + complex Hopf
   tower truncation question at the inner-factor level). The
   multi-year mathematical research project flagged in Sprint A7 lines
   232–242 — enriched motivic Galois group of the GeoVac discrete
   spectral triple — remains as the only non-sprint-scale path to Q5'.

5. **Combination rule $K = \pi(B + F - \Delta)$ unchanged**: this
   sprint affects only the period-ring side of Paper 55; Paper 2's
   numerical-observation status is preserved.

6. **No paper edits applied**:\ scoping only. The relevant Paper 55
   §subsec:open_m2_m3 already references Sprint A7's
   HALF-STRUCTURAL verdict (Paper 55 lines 1474–1499); the
   dimension-sweep NEGATIVE is a sharpening of that reading, not a
   change to it. If the PI wants, the §subsec:open_m2_m3 paragraph
   could be extended with one sentence:

   > Sprint Q5'-dim (2026-06-03; memo
   > \texttt{debug/sprint\_q5p\_dim\_sweep\_memo.md}) refines this to
   > NEGATIVE on the dimension-sweep entry point:\ the shared $\Q(i)$
   > persists across $d \in \{3, 5, 7\}$ for the same structurally
   > independent reasons (Witt-splitting at every $d$ for M2,
   > conductor-of-$\chi_{-4}$ at every $d$ for M3), with no
   > dimension-parity escape hatch. The HALF-STRUCTURAL verdict is
   > sharpened to "field coincidence is a consequence of two
   > independent one-bit invariants of the substrate that happen to
   > take the same value on every odd-dimensional sphere".

   This is a recommendation only;\ I have NOT applied it. PI may accept,
   modify, or decline per CLAUDE.md §13 sub-agent rule.

## Verdict

**NEGATIVE on the dimension-sweep entry point of Q5'.**

The dimensional sweep across $d \in \{3, 5, 7\}$ shows that M2's
$\Q(i)$ is **robustly generic-Euclidean** (Witt-splitting of the
standard positive-definite form over $\R^{d+1}$) at every $d$, and M3's
$\Q(\zeta_4)$ is **robustly conductor-of-$\chi_{-4}$** (cyclotomic
level 4) at every $d$. **No shared dimension-dependent structure
emerges**:\ the two number fields agree at $\Q(i)$ for two structurally
independent reasons at every $d$, neither of which depends on the
other. The HALF-STRUCTURAL Sprint A7 verdict is sharpened, not
upgraded:\ field coincidence is a consequence of two **independent
one-bit invariants** of the substrate (Euclidean signature for M2,
CH-shift denominator for M3) that happen to take the same value on
every odd-dimensional sphere. This **weakens Q5'**:\ if a deeper
Tannakian unification existed, the dimensional sweep would plausibly
have surfaced a mirror move; none is found.

The multi-year mathematical research project (enriched motivic Galois
group of the GeoVac discrete spectral triple) remains the only non-sprint-scale
path to closing Q5' in either direction.

## Files used

### Memos read
- `debug/sprint_a7_m2_m3_cyclotomic_memo.md` (HALF-STRUCTURAL verdict
  + S^5 robustness audit, lines 192–212)
- `debug/sprint_a2_s5_mixed_tate_memo.md` (S^5 M2 closed forms +
  dimension-parity sharpening, full memo)
- `debug/sprint_a6_m3_s5_memo.md` (S^5 M3 closed forms + Theorem A6.3
  cyclotomic mixed-Tate classification, full memo)
- `debug/sprint_mixed_tate_test_memo.md` (S^3 M2 baseline)

### Papers (consulted)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §4 (M2
  on S^3), §5 (M3 on S^3), §5.5 (M3 on S^5), §6 (joint structure),
  §7.5 (Q5' open question subsec)

### Code (none required)
This is structural reasoning on already-derived closed forms; no new
PSLQ, no new heat-trace expansion, no new periods identification.

### Published references (already in Paper 55 bibliography)
- Fathizadeh & Marcolli arXiv:1611.01815 Theorem 7.3 (Witt isotropy
  over $\Q(i)$, dimension-independent)
- Deligne arXiv:math/0302267 (period map at $N \in \{2, 3, 4, 6, 8\}$)
- Glanois arXiv:1411.4947 Cor. 1.1–1.2 (explicit level-4 basis)
- Camporesi & Higuchi *J. Geom. Phys.* 20 (1996) (CH spinor spectrum on
  odd-dim spheres, half-integer shifts $(d-1)/2$)

## One-line verdict

NEGATIVE on dimension-sweep entry point of Q5':\ the M2/M3 $\Q(i)$
coincidence persists across $d \in \{3, 5, 7\}$ via two independent
one-bit substrate invariants (Euclidean signature for M2;
CH-shift denominator for M3) that happen to evaluate identically on
every odd-dim sphere, with no dimension-parity escape hatch and no
mirror move between the two mechanisms.
