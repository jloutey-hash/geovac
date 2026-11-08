# Galois / number-field structure of the Ihara zeros — memo (Track RH-F)

Sprint: RH Sprint 3 ("arithmetic / Galois structure")
Author: Track RH-F, April 2026
Driver: `debug/compute_galois_ihara.py`
Data: `debug/data/galois_ihara.json`
Tests: `tests/test_galois_ihara.py` (27/27 passing)
Sources: Paper 29 §5, `debug/ihara_zeta_memo.md` §3, `debug/ihara_zeta_dirac_memo.md` §4

## 0. Question

Every non-trivial factor of each $\zeta_G(s)^{-1}$ studied so far (Paper 29)
has integer coefficients, so every Ihara zero is an algebraic number. Track
RH-F asks: is there a single "framework-preferred" number field (one small
extension of $\mathbb{Q}$) that hosts all Ihara zeros of the GeoVac graphs,
and if so, which one, and why?

Candidate fields were motivated by GeoVac structure:

| Field | Motivation |
|-------|------------|
| $\mathbb{Q}(i)$ | Hopf $S^1$ action (Paper 25 Hopf-gauge) |
| $\mathbb{Q}(\omega) = \mathbb{Q}(\sqrt{-3})$ | 3-prism / triangle cycles in S³ max_n=3 |
| $\mathbb{Q}(\sqrt{2})$ | $\kappa = -1/16 = 2^{-4}$; 2-fold spin |
| $\mathbb{Q}(\sqrt{5})$ | Fibonacci / golden-ratio spectrum (Hopf base) |
| $\mathbb{Q}(\sqrt{3})$ | diagnostic completeness |
| $\mathbb{Q}(\sqrt{-7})$ | cubic-quadratic resolvent candidate |
| $\mathbb{Q}(\zeta_8) = \mathbb{Q}(i,\sqrt 2)$ | 8-th roots of unity |
| $\mathbb{Q}(\zeta_{12}) = \mathbb{Q}(i,\sqrt 3)$ | 12-th roots (cyclotomic cover of $\mathbb{Q}(i)$ and $\mathbb{Q}(\omega)$) |

## 1. Polynomials

The 25 non-trivial polynomial factors studied, taken from Paper 29 and the
spinor extension memo:

| Source | Factors |
|:-------|:--------|
| S³ Coulomb max_n=3 | $s^2 \pm s + 1$, $2s^3 \pm s^2 + s \pm 1$ |
| S⁵ BS N_max=2 | $2s^2 + 1$, $3s^4 + 3s^2 + 1$, $24s^6 + 21s^4 + s^2 - 1$ |
| S⁵ BS N_max=3 | $P_{12}(s)$ (degree 12 even), $P_{22}(s)$ (degree 22 even) |
| Dirac-S³ Rule A n_max=2 | $s^2 + 1 = \Phi_4$ |
| Dirac-S³ Rule B n_max=2 | $3s^2+1, 4s^2+1, 2s^2 \pm s+1, 12s^4+9s^2-1$ |
| Dirac-S³ Rule A n_max=3 | $s^2+1, s^2 \pm s+1$, four cubics, two quartics |
| Dirac-S³ Rule B n_max=3 | $9s^2+1, P_{22}(s), P_{24}(s)$ |

## 2. Per-polynomial arithmetic table

Summary of Galois-group and minimal-splitting-field structure
(sympy's `galois_group` works to degree 6; for higher degree even-in-$s$
polynomials we Galois-analyze the $u = s^2$ reduction):

| Polynomial | deg | disc (sq-free core in parens) | $\text{Gal}/\mathbb{Q}$ | order | solvable | abelian | minimal splitting field |
|:-----------|:---:|:------------------------------|:------------------------|:-----:|:---------|:--------|:------------------------|
| $s^2 + 1 = \Phi_4$ | 2 | $-4$ ($-1$) | $\mathbb{Z}/2$ | 2 | yes | yes | $\mathbb{Q}(i)$ |
| $s^2 - s + 1 = \Phi_6$ | 2 | $-3$ ($-3$) | $\mathbb{Z}/2$ | 2 | yes | yes | $\mathbb{Q}(\omega)$ |
| $s^2 + s + 1 = \Phi_3$ | 2 | $-3$ ($-3$) | $\mathbb{Z}/2$ | 2 | yes | yes | $\mathbb{Q}(\omega)$ |
| $2s^2 + 1$ | 2 | $-8$ ($-2$) | $\mathbb{Z}/2$ | 2 | yes | yes | $\mathbb{Q}(\sqrt{-2}) \subset \mathbb{Q}(\zeta_8)$ |
| $3s^2 + 1$ | 2 | $-12$ ($-3$) | $\mathbb{Z}/2$ | 2 | yes | yes | $\mathbb{Q}(\omega)$ |
| $4s^2 + 1$ | 2 | $-16$ ($-1$) | $\mathbb{Z}/2$ | 2 | yes | yes | $\mathbb{Q}(i)$ |
| $2s^2 \pm s + 1$ | 2 | $-7$ ($-7$) | $\mathbb{Z}/2$ | 2 | yes | yes | $\mathbb{Q}(\sqrt{-7})$ |
| $9s^2 + 1$ | 2 | $-36$ ($-1$) | $\mathbb{Z}/2$ | 2 | yes | yes | $\mathbb{Q}(i)$ |
| $2s^3 \pm s^2 + s \pm 1$ (S³) | 3 | $-83$ ($-83$) | $S_3$ | 6 | yes | no | degree-6, resolvent $\mathbb{Q}(\sqrt{-83})$ |
| $2s^3 \pm 2s^2 + 2s \pm 1$ (Dirac-A) | 3 | $-44$ ($-11$) | $S_3$ | 6 | yes | no | degree-6, resolvent $\mathbb{Q}(\sqrt{-11})$ |
| $3s^4 + 3s^2 + 1$ | 4 | $432$ ($3$) | order 8 ($D_4$) | 8 | yes | no | degree-8, resolvent $\mathbb{Q}(\omega)$ |
| $12s^4 + 9s^2 - 1$ | 4 | $-3195072$ ($-3$) | order 8 ($D_4$) | 8 | yes | no | degree-8, resolvent $\mathbb{Q}(\omega)$ |
| $2s^4 \pm 2s^3 + 2s^2 \pm s + 1$ | 4 | $788$ ($197$) | order 24 ($S_4$) | 24 | yes | no | degree-24 |
| $24s^6 + 21s^4 + s^2 - 1$ | 6 | $2.50 \times 10^{11}$ | order 48 | 48 | yes | no | degree-48 |
| $P_{12}(s) = P_6(u=s^2)$ | 12 | $1.83 \times 10^{35}$ | $\text{Gal}(P_6) = S_6$ | 720 | **NO** | no | degree $\le 1440$ (2·720) |
| $P_{22}(s) = P_{11}(u=s^2)$ (S⁵, Dirac-B) | 22, 22 | vast | sympy limit exceeded | $\ge 720$ | unknown | no | unknown |
| $P_{24}(s) = P_{12}(u=s^2)$ (Dirac-B) | 24 | vast | sympy limit exceeded | unknown | unknown | unknown | unknown |

All polynomials are **irreducible over $\mathbb{Q}$**.

## 3. Galois groups per polynomial

Distinct Galois groups appearing:

- **$\mathbb{Z}/2$** (order 2, abelian): all 12 quadratic factors. The
  minimal splitting field is the imaginary quadratic field
  $\mathbb{Q}(\sqrt{D})$ where $D$ is the square-free core of the
  discriminant.
- **$S_3$** (order 6, solvable, non-abelian): all 6 cubic factors, discs
  $-83$ and $-44$. Resolvent fields $\mathbb{Q}(\sqrt{-83})$ and
  $\mathbb{Q}(\sqrt{-11})$. Splitting fields have degree 6 over $\mathbb{Q}$.
- **$D_4$ or $\mathbb{Z}/2 \times \mathbb{Z}/2 \rtimes \mathbb{Z}/2$**
  (order 8, solvable, non-abelian): quartics $3s^4+3s^2+1$ and $12s^4+9s^2-1$.
  Both reduce to $\mathbb{Z}/2$ on $u = s^2$; splitting field is the
  compositum of $\mathbb{Q}(\sqrt 3)$ or $\mathbb{Q}(\sqrt{-3})$ with a
  further $\sqrt{u\text{-root}}$ layer.
- **$S_4$** (order 24, solvable, non-abelian): the pair of quartics
  $2s^4 \pm 2s^3 + 2s^2 \pm s + 1$ in Dirac-A n_max=3. This is the first
  genuinely non-Z/2-powered solvable Galois group appearing in the
  GeoVac data.
- **Order 48** (solvable): the sextic $24s^6 + 21s^4 + s^2 - 1$ in
  S⁵ N_max=2. Sympy identifies it as a specific solvable group.
  The $u = s^2$ reduction gives a cubic with Galois $S_3$ (order 6).
- **$S_6$** (order 720, **non-solvable**): the $u = s^2$ reduction of
  $P_{12}(s)$ in S⁵ BS N_max=3. The first non-solvable Galois group in
  the GeoVac Ihara data. Radical expressions for its roots do not exist.
- **Unknown but large** ($\ge 720$): degree-22 and degree-24 $P$-polynomials
  in S⁵ BS N_max=3 and Dirac-B n_max=3. Sympy's `galois_group` caps at
  degree 6; deeper analysis would require GAP or Magma.

## 4. Minimal splitting field per Ihara zero set

Per-graph analysis of what single number field contains ALL zeros of
$\zeta_G(s)^{-1}$:

### S³ Coulomb, max_n=3

Factors: $\Phi_6, \Phi_3$ (live in $\mathbb{Q}(\omega)$, degree 2),
plus two $S_3$-cubics with discs $-83$ (degree-6 splitting field $K_3$,
containing $\mathbb{Q}(\sqrt{-83})$).

**Minimal splitting field of the full S³ Ihara zero set:**
$\mathbb{Q}(\omega, K_3)$, where $K_3$ is the $S_3$ splitting field of
$2s^3 - s^2 + s - 1$. Since $K_3$ has degree 6 over $\mathbb{Q}$ and
contains $\mathbb{Q}(\sqrt{-83}) \ne \mathbb{Q}(\omega)$, the compositum
has degree $6 \cdot 2 = 12$ over $\mathbb{Q}$ (check: $\omega \notin K_3$
because the quadratic sub-field of $K_3$ is $\mathbb{Q}(\sqrt{-83})$,
not $\mathbb{Q}(\sqrt{-3})$).

### S⁵ Bargmann–Segal, N_max=2

Factors: $2s^2 + 1$ (splits in $\mathbb{Q}(\sqrt{-2})$),
$3s^4 + 3s^2 + 1$ (degree-8 splitting field, contains
$\mathbb{Q}(\omega)$), $24s^6 + \ldots$ (degree-48 splitting field).

**Minimal splitting field:** contains
$\mathbb{Q}(\sqrt{-2}, \omega) = \mathbb{Q}(\sqrt{-2}, \sqrt{-3}) = \mathbb{Q}(\sqrt 2, \sqrt 3, i)$
(degree 8) plus further extensions from the sextic, giving a total degree
$\le 48 \cdot 2 = 96$ over $\mathbb{Q}$ (exact structure not computed).

### S⁵ Bargmann–Segal, N_max=3

Contains $P_{12}$ with $S_6$ Galois on $u = s^2$. $S_6$ is simple-modulo-$A_6$,
**non-solvable**. The minimal splitting field of $P_{12}$ has degree at most
$2 \cdot 720 = 1440$ over $\mathbb{Q}$. Combined with the unknown but
large-degree $P_{22}$ factor, the full zero set of S⁵ BS N_max=3 lives in
a number field of degree at least 720 over $\mathbb{Q}$.

### Dirac-S³ Rule A, n_max=3

Factors include $\Phi_4$, $\Phi_6$, $\Phi_3$, four $S_3$-cubics with discs
$-83$ and $-44$, and two $S_4$-quartics with disc $788 = 4 \cdot 197$.
Minimal splitting field contains at least $\mathbb{Q}(i, \omega)$ (degree 4)
plus $\mathbb{Q}(\sqrt{-11}, \sqrt{-83})$ (disjoint resolvents, so extra
factor of 4), giving degree $\ge 16$, plus the $S_4$ extension of degree 24,
so total $\ge 96$ over $\mathbb{Q}$.

### Dirac-S³ Rule B, n_max=3

Contains $9s^2 + 1$ (in $\mathbb{Q}(i)$) and two huge-degree $P$-polynomials
with sympy-uncomputable Galois groups. Minimal splitting field is at least
of high degree, structure not determined.

## 5. Structural observations

### The zero sets do NOT live in a single small field

No single small number field ($\mathbb{Q}(i)$, $\mathbb{Q}(\omega)$,
$\mathbb{Q}(\sqrt 2)$, $\mathbb{Q}(\sqrt 5)$, or even a cyclotomic cover
$\mathbb{Q}(\zeta_M)$ for small $M$) hosts more than about $8{-}10$ of
the 25 polynomials studied. The "Paper 25 Hopf-$U(1)$" reading would
predict that all zeros should live in an abelian extension of $\mathbb{Q}$
(i.e. inside some $\mathbb{Q}(\zeta_M)$ by the Kronecker–Weber theorem).
This is **contradicted** by the $S_3, S_4, S_6$ Galois groups:

1. **$S_3$ is non-abelian**: every cubic factor has $S_3$ Galois. So the
   splitting fields of the cubics are genuinely non-abelian over
   $\mathbb{Q}$ and do not embed in any $\mathbb{Q}(\zeta_M)$.
2. **$S_4$ is non-abelian**: the Dirac-A quartic pair has $S_4$ Galois.
3. **$S_6$ is non-solvable** (in $P_{12}$, S⁵ N_max=3). This is a strong
   negative result: not merely non-abelian, the Galois group has no
   subnormal series with abelian quotients. The zeros of $P_{12}$ cannot
   be expressed in radicals over $\mathbb{Q}$.

### What IS true: each Ramanujan-critical quadratic lives in $\mathbb{Q}(i)$

When the Ramanujan bound is saturated (Dirac-B Rule B at $n_{\max}=2$ has
boundary zeros at $s = \pm i/2$), the critical-circle factor is $4s^2+1$
which splits over $\mathbb{Q}(i)$. Same for $9s^2+1$ at Dirac-B n_max=3
and $\Phi_4 = s^2+1$ in both Dirac-A cases. This is consistent with the
Paper-25 reading that the Hopf $U(1) \to \mathbb{Z}_2$ reduction (to real
integer adjacency) naturally produces $\pm i/n$ boundary zeros, whose
minimal field is exactly $\mathbb{Q}(i)$.

### The cyclotomic-polynomial factors are exactly the per-$\ell$ cycle
    counts

The Paper 29 §5.1 per-$\ell$-shell decomposition of the S³ Ihara zeta is
arithmetically mirrored: the S³ Coulomb max_n=3 zeta contains $\Phi_3$ and
$\Phi_6$ (the primitive 3rd and 6th roots of unity) because the $\ell=1$
component of the graph is the 3-prism $C_3 \times P_2$, which has
3-cycles (roots $\mapsto \omega, \bar\omega$) and 6-cycles
(roots $\mapsto \Phi_6$). **The cyclotomic factor pattern is a direct
fingerprint of the cycle structure of the underlying graph.**

### The discriminants $-83$ and $-44$ are not GeoVac-natural

The cubic discriminants $-83$ (S³ Coulomb and Dirac-A Rule A) and $-44$
(Dirac-A Rule A) are prime-ish (with 44 = 4·11) and carry no transparent
framework meaning. $83$ and $11$ do not appear in Paper 2's
$(B = 42, F = \pi^2/6, \Delta = 1/40)$ decomposition, nor in Paper 18's
exchange-constant taxonomy. They appear here as a numerical consequence
of the specific 3-prism and Dirac ladder geometry, not as injected
transcendentals. These are not "arithmetic constants of GeoVac"; they
are graph-specific invariants.

### The S_6 non-solvability closes the radical-expressions avenue

For the S⁵ Bargmann-Segal graph at N_max=3, the 12 Ihara zeros of $P_{12}$
have a non-solvable $S_6$ Galois group on the $u = s^2$ reduction.
**This means no closed-form radical expression for these zeros exists.**
Any attempt to express them via $\sqrt[n]{\ldots}$ chains is
mathematically obstructed, analogous to the Abel-Ruffini theorem for
degree-5 polynomials.

The analogous situation for $P_{22}$, $P_{24}$ in the larger graphs is
presumably at least as bad, and likely worse (degree 11 or 12 in $u$ with
generic Galois $S_{11}$ or $S_{12}$).

### Connection to Paper 25's Hopf-$U(1)$

Paper 25 §VII.1 reports that at N_max=5 the S⁵ Bargmann-Segal graph's
natural gauge group is $U(1)$, not $SU(3)$, because transitions between
$(N,0)$ and $(N+1,0)$ irreps are Clebsch-Gordan intertwiners, not group
elements. The $m_\ell \to -m_\ell$ reflection is the only sub-action of
the continuous $U(1)$ that commutes with a real integer adjacency — and
this is exactly what produces the $12+22$ dichotomy of Paper 29 §5.3.

**The $\mathbb{Z}_2$-block decomposition survives into the Galois
picture** in the limited sense that $P_{12}$ (the anti-symmetric block)
and $P_{22}$ (the symmetric block) are each separately irreducible over
$\mathbb{Q}$. But the actual arithmetic of the zeros within each block
is **much richer than a $\mathbb{Z}_2$ or any abelian structure would
allow**: $P_{12}$'s Galois is $S_6$ (order 720, non-solvable), not some
abelian quotient.

## 6. Open questions and next steps

1. **Compute the Galois group of $P_{22}$ and $P_{24}$ in GAP/Magma**
   (sympy is capped at degree 6). The even-in-$s$ reduction gives
   $u$-polynomials of degrees 11 and 12, still beyond sympy but within
   GAP's reach. If these are also $S_n$ (i.e. generic), then the
   non-solvability extends to the largest GeoVac Ihara zeros.

2. **Try the Frobenius element approach.** For a polynomial $P(s)$ with
   integer coefficients, Chebotarev's density theorem identifies
   Galois-conjugacy-class frequencies with prime-residue patterns. A
   statistical sampling (factoring $P_{12}(s) \bmod p$ for the first
   1000 primes and counting cycle types) could determine the Galois
   group of $P_{12}$ independently of sympy's built-in methods.

3. **Look for hidden cyclotomic structure at the Ramanujan boundary.**
   Every polynomial whose zeros are ON the critical circle
   $|s| = 1/\sqrt{q_{\max}}$ is a cyclotomic polynomial evaluated at
   $s = \zeta / \sqrt{q_{\max}}$: e.g. $4s^2 + 1 = \Phi_4(2s)$ up to
   normalization. A systematic catalog of which graphs have
   cyclotomic-boundary factors would refine Paper 29's Observation 1.

4. **Paper 25 $\mathbb{Z}_2$ action preserves irreducibility but does
   NOT reduce Galois complexity.** Is there a deeper (non-linear)
   framework symmetry that would block-diagonalize $P_{12}$ further into
   abelian factors? Conjecturally: no — the $S_6$ non-solvability is
   structurally intrinsic. But a Frobenius-based sampling test could
   confirm that the $S_6$ is not an artifact of a missed
   block-diagonalization.

5. **Cross-check with Paper 18's exchange-constant taxonomy.** Every
   Ihara zero is algebraic, not transcendental; yet the non-solvable
   $S_6$ says the algebraic structure is genuinely non-radical. Paper 18's
   classification (intrinsic / calibration / embedding / flow) would
   naturally append a fifth tier for "algebraic but non-radical": these
   are structurally intrinsic but cannot be expressed by a finite chain
   of root extractions. The Ihara zeros of $P_{12}$ are the first
   candidate for this tier.

## 7. Summary table — does a single framework-preferred ring exist?

**No.** No small number field hosts all Ihara zeros:

- $\mathbb{Q}(i)$ hosts 4 of 25 polynomial factors (pure-imaginary
  roots, disc = square × $(-1)$).
- $\mathbb{Q}(\omega)$ hosts 6 of 25 (cyclotomic $\Phi_3, \Phi_6$ and
  disc $= -3 \cdot (\text{square})$ cases).
- $\mathbb{Q}(\zeta_{12}) = \mathbb{Q}(i, \omega)$ hosts 10 of 25
  (union of the above).
- $\mathbb{Q}(\sqrt{-7})$ hosts 2 (Dirac-B n_max=2 cubic cuspids).
- No field hosts the $S_3$, $S_4$, or $S_6$-Galois factors, because
  these require non-abelian extensions not contained in any cyclotomic
  field by Kronecker–Weber.

**What IS universal:** every Ihara zero is an algebraic number; every
polynomial factor is irreducible over $\mathbb{Q}$ with integer
coefficients; **no transcendental ($\pi, \zeta(2), \zeta(3), G, \beta(4)$)
appears**. This is the Paper 29 Corollary 1 ($\pi$-freeness of zeros)
at the finest arithmetic level we can currently reach.

The headline negative: **the GeoVac Hopf Ihara zeros ARE NOT an arithmetic
object in the Kronecker–Weber sense.** The Paper 25 Hopf-$U(1)$ abelian
gauge structure controls the $\mathbb{Z}_2$ block decomposition of
$P_{12} + P_{22}$ but does not reduce the Galois group of each block to
an abelian quantity. The zeros sit in a genuinely non-abelian, and in
the S⁵ N_max=3 case non-solvable, region of the algebraic-number lattice.

## 8. Status

- Driver: `debug/compute_galois_ihara.py` — 27 polynomials processed.
- Data: `debug/data/galois_ihara.json` (33,892 bytes).
- Tests: `tests/test_galois_ihara.py` — 27/27 passing.
- Paper 29 addition proposed (not auto-applied): §5.4 or §6.x sub-section
  "Galois structure of the Ihara zeros." Flagged for plan-mode review.
- CLAUDE.md: no change required (Paper 29 is in `papers/observations/` as
  "synthesis observation"; the Galois analysis strengthens Corollary 1
  without modifying any core claim).
