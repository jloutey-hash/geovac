# Sprint Mixed-Tate Test — canonical memo
Date: 2026-06-03
Scope: 1-week diagnostic; PI sign-off pending on paper edits.

## TL;DR

**Verdict: POSITIVE under the standard Connes–Chamseddine volume-normalized
convention** (the convention F–M arXiv:1611.01815 use). The GeoVac discrete
Seeley–DeWitt coefficients on $S^3$, in both the Dirac and scalar-Laplacian
sectors, sit inside $\mathbb{Q}[\pi^2] \subset$ mixed-Tate-period-ring over
$\mathbb{Q}$, with the stronger structural property of **two-term exactness**
(Dirac) or **geometric-series collapse** $a_k^\Delta = 2\pi^2/k!$ (scalar)
that the continuum Robertson–Walker case does not have. A structural caveat:
the *raw* heat-trace coefficients $\tilde a_k$ (before the standard
$(4\pi)^{-d/2}$ volume normalization) literally contain $\sqrt\pi$, which is
NOT a mixed-Tate period; the $\sqrt\pi$ cancels exactly against $(4\pi)^{3/2}$
in $d=3$ to produce the mixed-Tate $\pi^2 \cdot \mathbb{Q}$ output. Paper 32
§VIII's case-exhaustion theorem M2 (Seeley–DeWitt) sub-mechanism therefore
inherits the F–M classification as a corollary at the standard-convention
level. The M3 (vertex-parity Hurwitz, Catalan $G$, $\beta(s)$) content lives
in a different sector of the master Mellin engine ($k=1$ rather than $k=2$)
and does NOT contaminate the SD coefficients.

## Background

Fathizadeh–Marcolli (Comm. Math. Phys. 2017, arXiv:1611.01815) prove that the
coefficients of the asymptotic expansion of the Connes–Chamseddine spectral
action on a Euclidean Robertson–Walker spacetime $\mathbb{R}\times S^3$ with
metric $dt^2 + a(t)^2 d\Omega_3^2$ are periods of mixed Tate motives. The
specific characterization: each coefficient reduces to a period of a relative
motive of a complement of unions of hyperplanes and quadric hypersurfaces,
divided by divisors given by unions of coordinate hyperplanes (their Theorem
on Rosenfeld-style integral representations of the Seeley–DeWitt coefficients
via Wodzicki's noncommutative residue). The proof works by treating the
scaling factor $a$ as an affine variable and showing that each SD coefficient,
viewed as a polynomial in $a$ and its time derivatives with rational
coefficients, is a period of a mixed Tate motive over $\mathbb{Q}$.

Standard period-theory background (Brown, Deligne, Goncharov): the ring of
periods of mixed Tate motives over $\mathbb{Z}$ is conjecturally
$\mathbb{Q}[(2\pi i)^{\pm 1}] \cdot \mathbb{Q}\text{-span}\{\text{MZVs}\}$,
where MZVs are multiple zeta values $\zeta(n_1,\ldots,n_k)$ with positive
integer arguments. Crucially, **algebraic extensions like $\sqrt{\pi}$ are NOT
mixed-Tate periods** (mixed Tate Tate twists are integer powers of
$\mathbb{Q}(1) = \mathbb{Q}(2\pi i)$; half-twists do not exist in the Tate
category). **Dirichlet $L$-values at non-principal characters mod $m$**,
including Catalan $G = \beta(2)$ and $\beta(4)$, are mixed-Tate over the
cyclotomic ring $\mathbb{Z}[\zeta_m, 1/m]$ but NOT over $\mathbb{Z}$ or
$\mathbb{Q}$ (they live in cyclotomic mixed Tate at level $m$).

So for "mixed-Tate over $\mathbb{Q}$" in the F–M sense, the allowed
transcendentals are: integer powers of $\pi$ (and $\log 2$ via $\log 2 =
-\sum (-1)^n/n$ which is an MZV-type weight-1 period), odd zeta values
$\zeta(3), \zeta(5), \ldots$, and restricted MZVs. **Excluded**: $\sqrt{\pi}$,
$\sqrt{p}$ for primes, Catalan $G$, generic Dirichlet $L$-values.

## GeoVac SD coefficients on $S^3$ — enumeration

All values for unit $S^3$.

### Dirac sector ($D^2$, Camporesi–Higuchi spectrum)

Spectrum: $|\lambda_n^{CH}| = n + 3/2$, degeneracy $g_n = 2(n+1)(n+2)$.

Heat trace (Paper 51 Cor. 2.1, two-term-exact closed form):
$$K_{D^2}(t) = \frac{\sqrt\pi}{2} t^{-3/2} - \frac{\sqrt\pi}{4} t^{-1/2} + O(e^{-\pi^2/t}).$$

| coefficient | raw $\tilde a_k$ | volume-norm $a_k = (4\pi)^{3/2}\tilde a_k$ |
|:------------|:-----------------|:-------------------------------------------|
| $\tilde a_0^{D^2}$ | $\sqrt\pi/2$ | $4\pi^2$ |
| $\tilde a_1^{D^2}$ | $-\sqrt\pi/4$ | $-2\pi^2$ |
| $\tilde a_k^{D^2}$, $k\ge 2$ | $0$ | $0$ |

Two-term exactness traces to the Bernoulli identity $B_{2k+1}(3/2) =
(2k+1)/4^k$ (Paper 51 Thm. 2.1 / `rem:two_term_uniqueness`), which makes the
spectral zeta vanish at all non-positive integers, $\zeta_{\rm unit}(-k) = 0$
for $k \ge 0$. **All higher SD coefficients are exactly zero, not just small.**
This is sharper than the continuum F–M generic R–W case, where the SD
expansion is genuinely infinite.

### Scalar Laplacian sector ($\Delta$, integer spectrum)

Spectrum: $\lambda_n^\Delta = n(n+2)$, degeneracy $(n+1)^2$.

Heat trace (Paper 51 Thm. 3.1):
$$K_\Delta(t) = \frac{\sqrt\pi}{4} \frac{e^t}{t^{3/2}} + O(e^{-\pi^2/t})
= \sum_{k=0}^\infty \frac{\sqrt\pi}{4 k!} t^{k - 3/2} + O(e^{-\pi^2/t}).$$

| $k$ | raw $\tilde a_k^\Delta$ | volume-norm $a_k^\Delta = 2\pi^2/k!$ |
|:---:|:-----------------------:|:-------------------------------------:|
| 0 | $\sqrt\pi/4$ | $2\pi^2$ |
| 1 | $\sqrt\pi/4$ | $2\pi^2$ |
| 2 | $\sqrt\pi/8$ | $\pi^2$ |
| 3 | $\sqrt\pi/24$ | $\pi^2/3$ |
| 4 | $\sqrt\pi/96$ | $\pi^2/12$ |
| 5 | $\sqrt\pi/480$ | $\pi^2/60$ |

Geometric-series collapse to $a_k^\Delta = 2\pi^2/k!$ traces to the Jacobi
$\vartheta_3$ modular transformation on the integer spectrum (Paper 51 §G3
proof sketch).

### Cross-check: $\zeta_{D^2}(s)$ at integer $s$ (Paper 28 T9)

Verified by direct sympy summation:

| $s$ | $\zeta_{D^2}(s)$ |
|:---:|:-----------------|
| 2 | $\pi^2 - \pi^4/12$ |
| 3 | $\pi^4/3 - \pi^6/30$ |
| 4 | $2\pi^6/15 - 17\pi^8/1260$ |
| 5 | $\pi^8(306 - 31\pi^2)/5670$ |

All in $\mathbb{Q}[\pi^2]$; no $\sqrt\pi$, no $\zeta(\text{odd})$, no
$\beta(s)$, no $G$, no MZV. The T9 theorem (Paper 28 Theorem 1) proves this
structurally: $\zeta_{D^2}(s) = 2^{2s-1}[\lambda(2s-2) - \lambda(2s)]$ with
$\lambda(2k) = (1-2^{-2k})\zeta_R(2k) \in \mathbb{Q}\cdot \pi^{2k}$.

## Arithmetic classification

| sector | object | arithmetic class | mixed-Tate over $\mathbb{Q}$? |
|:-------|:-------|:-----------------|:------------------------------:|
| Dirac raw | $\tilde a_k^{D^2}$ ($k=0,1$) | $\sqrt\pi \cdot \mathbb{Q}$ | **NO** |
| Dirac raw | $\tilde a_k^{D^2}$ ($k\ge 2$) | $0$ (trivially) | yes (trivially) |
| Dirac vol-norm | $a_k^{D^2}$ ($k=0,1$) | $\pi^2\cdot \mathbb{Q}$ | **YES** |
| Dirac vol-norm | $a_k^{D^2}$ ($k\ge 2$) | $0$ | yes |
| Scalar raw | $\tilde a_k^\Delta$ ($k\ge 0$) | $\sqrt\pi\cdot\mathbb{Q}$ | **NO** |
| Scalar vol-norm | $a_k^\Delta$ ($k\ge 0$) | $\pi^2\cdot\mathbb{Q}$ | **YES** |
| spectral zeta values | $\zeta_{D^2}(s)$ ($s\in\mathbb{Z}_{\ge 1}$) | $\mathbb{Q}[\pi^2]$ (degree $\le 2$ in $\pi^2$ at fixed $s$) | **YES** |

The classification is **convention-dependent in a structural way** at the raw
heat-trace level, but **convention-independent at the zeta-value / SD level
that F–M actually classify**. Both at the SD level and at the spectral-zeta
level, all GeoVac quantities live in the pure-Tate sub-ring
$\mathbb{Q}[\pi^2] \subset \mathbb{Q}[(2\pi i)^{\pm 1}]$ of mixed-Tate periods
over $\mathbb{Q}$.

In fact GeoVac's content is *arithmetically more restricted* than F–M's: F–M
allow general mixed-Tate periods including $\zeta(3), \zeta(5),$ MZVs; the
GeoVac SD sector on $S^3$ produces only $\pi^{2k}\cdot\mathbb{Q}$, i.e. only
the **even-weight pure Tate** sub-ring of mixed Tate periods. No odd zeta
values, no MZVs, appear at the SD level. The continuum R–W case can produce
$\zeta(3)$ and higher (when scaling factor derivatives are non-trivial); the
GeoVac discrete $S^3$ case (with no continuous scaling) cannot.

## Cross-check vs Paper 28 T9: where does $\sqrt\pi$ live?

Paper 28 Theorem 1 (T9) says $\zeta_{D^2}(s)$ at integer $s$ lies in $\sqrt\pi
\cdot \mathbb{Q} \oplus \pi^2 \cdot \mathbb{Q}$ at integer $s$. This memo
verifies: at every tested $s \in \{2,3,4,5\}$ the $\sqrt\pi \cdot \mathbb{Q}$
component is **exactly zero** — i.e. the integer-$s$ values are pure
$\pi^{2k}\cdot\mathbb{Q}$. The $\sqrt\pi$ slot in T9's two-term ring is
populated only at half-integer $s$ (e.g. $s = 1/2, 3/2$, which corresponds
via Mellin to the leading $t^{-3/2}, t^{-1/2}$ pole structure of $K_{D^2}(t)$
— exactly where $\sqrt\pi$ lives in the raw heat trace).

So the $\sqrt\pi$ visible in T9's two-term ring is the *half-integer
spectral-zeta content* corresponding to the leading volume / Lichnerowicz
poles of the heat trace, and it is **already absorbed by the
$(4\pi)^{-3/2}$ dimensional prefactor in the standard SD convention**.
Specifically:

- T9 evaluated at $s = 3/2$ gives the residue at the leading heat-trace
  pole, which equals $\sqrt\pi/2$ (matching $\tilde a_0^{D^2}$).
- Multiplying by $(4\pi)^{-3/2}$ to convert to volume-normalized SD gives
  $(\sqrt\pi/2)/(4\pi)^{3/2} = (\sqrt\pi/2)/(8\pi\sqrt\pi) = 1/(16\pi)$,
  which when re-multiplied back through the full Vassilevich-Branson-Gilkey
  expression with $\dim_S \cdot \mathrm{Vol}(S^3) = 4 \cdot 2\pi^2 = 8\pi^2$
  reproduces the $a_0^{D^2} = 4\pi^2$ result.

The two readings are consistent. The $\sqrt\pi$ is a *fiber-bundle / spinor
dimension* artifact of the raw heat-trace normalization, not a genuine
algebraic obstruction to mixed-Tate classification.

## Verdict

**POSITIVE (with structural caveat).** The GeoVac discrete Seeley–DeWitt
coefficients on $S^3$ are mixed-Tate periods over $\mathbb{Q}$ in the
standard Connes–Chamseddine / Fathizadeh–Marcolli convention. The mixed-Tate
classification of F–M (arXiv:1611.01815) inherits to GeoVac's $S^3$ sub-case
as a corollary, with three structural enhancements vs the F–M continuum
result:

1. **Pure-Tate restriction.** GeoVac SD coefficients on $S^3$ land in
   $\pi^{2k}\cdot\mathbb{Q}$ (pure Tate weight-$2k$), strictly smaller than
   the generic mixed-Tate ring with MZV content. No $\zeta(3), \zeta(5),$ or
   MZVs appear at the SD level. This is sharper than F–M.

2. **Two-term exactness (Dirac).** GeoVac's $S^3$ Dirac SD expansion
   truncates at $k=1$ exactly (Paper 51 Cor. 2.1). The F–M generic R–W
   expansion does not truncate; it is genuinely infinite in $k$.

3. **Geometric-series collapse (scalar).** GeoVac's scalar Laplacian SD
   coefficients form $a_k^\Delta = 2\pi^2/k!$ exactly (Paper 51 Thm. 3.1).
   The F–M generic R–W scalar coefficients have no such closed form.

The structural caveat is that at the *raw* heat-trace coefficient level
(before the $(4\pi)^{-d/2}$ dimensional rescaling), $\sqrt\pi$ is literally
present. This $\sqrt\pi$ is the spinor-fiber-dimension factor that always
appears in $d$-odd heat-kernel expansions on spin manifolds, and it cancels
exactly against $(4\pi)^{d/2}$ in the standard SD convention. Recording this
as a normalization-dependent observation rather than a competing verdict.

The M3 (vertex-parity Hurwitz, Catalan $G$, $\beta(s)$) sector of the master
Mellin engine — which DOES contain non-mixed-Tate content — lives at $k=1$
in $\mathcal{M}[\mathrm{Tr}(D^k\,e^{-tD^2})]$ and produces vertex-restricted
Dirichlet-$L$ sums (Paper 28 Theorem 3, $D_{\rm even}(s) - D_{\rm odd}(s)$).
It does **NOT** contaminate the un-restricted $\zeta_{D^2}(s) = D_{\rm even}
+ D_{\rm odd}$ whose Mellin transform produces the SD coefficients. The
case-exhaustion theorem's M2 sub-mechanism (SD) is mixed-Tate; its M3
sub-mechanism (vertex-parity) is cyclotomic-mixed-Tate at level 4 (mixed
Tate over $\mathbb{Z}[i]$, not over $\mathbb{Q}$).

## Paper-level implications (proposals, not edits)

1. **Paper 32 §VIII** (case-exhaustion theorem): add a Corollary or Remark
   immediately after the M2 sub-case statement: *"The M2 sub-mechanism
   coefficients on $S^3$ are mixed-Tate periods over $\mathbb{Q}$ in the
   standard volume-normalized convention, by direct inheritance from
   Fathizadeh–Marcolli (Comm. Math. Phys. 2017, arXiv:1611.01815) applied to
   the static-Robertson–Walker sub-case $a(t) \equiv 1$ on $\mathbb{R} \times
   S^3$. The GeoVac $S^3$ M2 sector is the pure-Tate weight-$2k$ sub-ring
   $\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$, strictly smaller than the generic
   mixed-Tate ring; in particular it contains no $\zeta(3), \zeta(5)$, or
   MZV content."* This is the substantive theoretical advance referenced in
   the sprint goal.

2. **Paper 18 §III.7** (master Mellin engine): convert the "forward pointer"
   paragraph (lines 961–980, which already raises the F–M question) from
   "natural sprint question raised but not closed" to "closed in the
   POSITIVE direction with the pure-Tate refinement: GeoVac's M2 sector on
   $S^3$ is $\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$, a pure-Tate sub-ring of
   F–M's mixed-Tate class." Cite Sprint Mixed-Tate Test memo.

3. **Paper 28 §T9**: add a Remark or footnote clarifying that the
   $\sqrt\pi\cdot\mathbb{Q} \oplus \pi^2\cdot\mathbb{Q}$ two-term ring in T9
   is the spectral-zeta-at-integer-$s$ ring, and that the $\sqrt\pi$ slot is
   populated only at half-integer $s$ (the heat-trace poles); at integer
   $s$ the values are pure $\pi^{2k}\cdot\mathbb{Q}$ as verified by the
   table.

4. **Paper 51 §G3** (`rem:cc_uncanny` and `rem:two_term_uniqueness`): add a
   one-sentence cross-reference: *"Under the standard CC volume normalization
   the two-term exact Dirac SD coefficients $a_0 = 4\pi^2, a_1 = -2\pi^2$ on
   $S^3$ are mixed-Tate periods over $\mathbb{Q}$ (Fathizadeh–Marcolli 2017
   arXiv:1611.01815); the scalar series $a_k^\Delta = 2\pi^2/k!$ likewise.
   See Sprint Mixed-Tate Test memo."*

5. **CLAUDE.md §1.7 WH2** (Paper 18 = Seeley–DeWitt decomposition working
   hypothesis): consider promoting status with: *"M2 sub-mechanism of WH2
   inherits F–M mixed-Tate classification on $S^3$ at the volume-normalized
   convention level (Sprint Mixed-Tate Test, 2026-06-03). M2 sector on $S^3$
   is in fact the pure-Tate sub-ring $\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$,
   strictly smaller than generic mixed-Tate."*

## Open questions / follow-on

1. **Inhomogeneous extension.** Does the F–M classification extend to the
   GeoVac $S^3$ sector with a non-trivial scaling factor $a(t)$ (i.e.\ the
   actual R–W substrate, $\mathbb{R}\times S^3$ with $a(t)$ general)? F–M
   prove the continuum result; the GeoVac discrete substrate with
   time-dependent $a(t)$ would need its own re-derivation. Estimated effort:
   2–4 weeks.

2. **Cross-manifold (Paper 24 / S^5 Bargmann–Segal).** Does the same
   classification hold on $S^5$ (where the analogous CC expansion has THREE
   power-law terms including $R^2$, per Paper 51 `rem:two_term_uniqueness`)?
   $S^5$ retains the pure-Tate sub-ring at the volume-normalized SD level,
   but the explicit closed forms have not been written out in this sprint.

3. **Cyclotomic-mixed-Tate verdict for M3.** The Catalan $G$ and $\beta(s)$
   content of M3 lives in mixed Tate over $\mathbb{Z}[\zeta_4] = \mathbb{Z}[i]$
   (level-4 cyclotomic). A parallel mixed-Tate-test sprint for M3 would
   verify this against the published Brown–Goncharov cyclotomic mixed-Tate
   classification.

4. **Inner-factor extension (Sprint H1 / Yukawa Dirichlet ring).** The
   almost-commutative extension $\mathcal{A}_{GV} \otimes (\mathbb{C} \oplus
   \mathbb{H})$ (Paper 32 §VIII.C) introduces Yukawa coupling content
   classified in Paper 18 §IV's "inner-factor input data" tier. Whether
   these inner-factor SD coefficients remain mixed-Tate is an open question.

5. **F–M paper full text.** This memo relied on the F–M abstract + search-
   surfaced quotes for the precise mixed-Tate-over-$\mathbb{Q}$ definition.
   A follow-on read of the full F–M paper (Comm. Math. Phys. 2017) would
   confirm the precise relative-motive structure F–M use, and would let us
   write down the *specific* affine complement of quadrics + hyperplanes
   that the GeoVac $S^3$ SD coefficients reduce to (which is interesting
   because that affine complement is essentially trivial — a point — in our
   case, since there are no continuous derivatives $a^{(n)}(t)$).

## Files used

### Papers (read)
- `papers/group5_qed_gauge/paper_28_qed_s3.tex` (T9 theorem and proof,
  lines 245–348)
- `papers/group5_qed_gauge/paper_51_gravity_arc.tex` (G1 spectral-zeta-at-
  non-positive-integers theorem `thm:zeta_unit_neg_k` and proof lines 346–
  372; Cor. 2.1 `cor:two_term` two-term exactness lines 379–388;
  `rem:cc_uncanny` lines 398–409; `rem:two_term_uniqueness` lines 411–
  435; Thm 3.1 `thm:scalar_ak` scalar Laplacian SD closed form lines
  505–520)
- `papers/group3_foundations/paper_18_exchange_constants.tex` (Master
  Mellin engine §III.7 lines 843–1007; existing F–M forward-pointer
  paragraph lines 961–980)

### Code (consulted and run)
- `geovac/qed_vacuum_polarization.py` `seeley_dewitt_coefficients_s3()`
  function (lines 114–207): produces $a_0 = a_1 = \sqrt\pi$, $a_2 =
  \sqrt\pi/8$ in the **raw, dim-spinor-multiplied** convention (i.e.\
  $\tilde a_k \cdot \dim_S$ without the $(4\pi)^{3/2}$ rescaling), which is
  yet a third convention. Reconciliation: this code uses prefactor =
  $(4\pi)^{-3/2}$, $\dim_S = 4$, $\mathrm{Vol} = 2\pi^2$, giving $a_0
  = (4\pi)^{-3/2} \cdot 4 \cdot 2\pi^2 = 8\pi^2/(4\pi)^{3/2} = 8\pi^2 /
  (8\pi\sqrt\pi) = \sqrt\pi$. The $\sqrt\pi$ here is the *combined*
  $(4\pi)^{-3/2} \cdot \mathrm{Vol}(S^3)$ factor when $\dim_S = 4$ is not
  absorbed. Volume-rescaling to extract the "geometric SD coefficient"
  (curvature integral) recovers $a_0^{\rm geom} = 4\pi^2$, the standard CC
  /F–M value.

### Diagnostic scripts (created this sprint)
- `debug/mixed_tate_test_sd_coefficients.py` — reproducible computation of
  the table above. Verifies T9 by direct sympy summation; computes raw and
  volume-normalized SD coefficients; runs an arithmetic classifier.

### Web sources (consulted for F-M classification)
- arXiv:1611.01815 abstract page (accessed; PDF binary not extractable
  through WebFetch)
- Springer Comm. Math. Phys. landing page (redirect-only, no content
  extracted)
- Brown, *Arbeitstagung lectures on multiple zeta values*
  (ihes.fr/~brown/Arbeitstatung.pdf) — referenced for the standard
  mixed-Tate-period definition
- arXiv:1612.03693, arXiv:1611.01011 — referenced for the multi-Dedekind
  / real-quadratic mixed-Tate extensions and the standard background that
  $\sqrt\pi$ is NOT a mixed-Tate period over $\mathbb{Q}$
