# Sprint M1 Pure-Tate — canonical memo
Date: 2026-06-03
Scope: closure of the master Mellin engine trinity at the
Paper 32 §VIII corollary level; companion to Sprint Mixed-Tate
Test (M2) and Sprint M3 Cyclotomic Mixed-Tate (M3), both
2026-06-03.

## TL;DR

**Verdict: POSITIVE (trivially).** The M1 sub-mechanism of the
master Mellin engine (the Hopf-base measure factor
$\mathrm{Vol}(S^2)/4 = \pi$ on the $S^3 \to S^2$ Hopf bundle) produces
period values in the localised pure-Tate sub-ring
$\mathbb{Q}[\pi, \pi^{-1}]$ of mixed Tate periods over $\mathbb{Q}$,
at depth $0$ and Tate weight $\in \{-2, -1, 0, 1, 2, \ldots\}$.
Every M1 closed form is a $\mathbb{Q}$-linear combination of integer
powers of $\pi$.

This sprint is bookkeeping — the M1 closure is observation, not
investigation.  The motivic content of $\mathrm{Vol}(S^2)/4 = \pi$
has been clear since Hopf 1931; the Paper 32 §VIII proof of
Theorem~\ref{thm:pi_source_case_exhaustion} already writes
"M1's $\mathrm{Vol}(S^2)/4 = \pi$ identity" as one of the three
mechanism-side terms.  What this sprint adds is a formal Paper 32
§VIII corollary parallel to `cor:m2_mixed_tate` and
`cor:m3_cyclotomic_mixed_tate`, so the master Mellin engine
classification reads symmetrically across the three sub-mechanisms.

Combined with the morning's M2 and M3 sprints, the master Mellin
engine partition is now a **complete period-theoretic
classification** at the case-exhaustion-theorem-corollary level:

| Mellin index | Sub-mechanism | Period-ring home | Depth |
|:---:|:---|:---|:---:|
| $k = 0$ | M1 (Hopf-base measure) | $\mathbb{Q}[\pi, \pi^{-1}]$ | 0 |
| $k = 2$ | M2 (Seeley–DeWitt) | $\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$ pure Tate | 0 |
| $k = 1$ | M3 (vertex parity / Hurwitz / Dirichlet $L$) | $\mathcal{MT}(\mathbb{Z}[i, 1/2])$ at level $\le 4$ | $\ge 1$ |

Every transcendental that has appeared in GeoVac spectral data has
a named motivic home, indexed by the master Mellin engine slot
$k \in \{0, 1, 2\}$.

## Background

The M1 sub-mechanism appears in the master Mellin engine
$\mathcal{M}[\mathrm{Tr}(D^k \cdot e^{-tD^2})]$ at $k = 0$: the
trivial heat kernel collapses to the volume integral over the
substrate, with the Hopf-base measure
$\mathrm{Vol}(S^2)/4 = \mathrm{Vol}(S^3)/\mathrm{Vol}(S^1) = \pi$
appearing as the prefactor identified in Paper 32 §VIII proof
(line 1703 of the current `paper_32_spectral_triple.tex`:\ "M1's
$\mathrm{Vol}(S^2)/4 = \pi$ identity").

Motivic content of $\pi$ (standard Kontsevich–Zagier 2001,
Belkale–Brion 2003, Brown 2017 mixed Tate framework):\ $\pi$ is a
period of mixed Tate motives over $\mathbb{Z}$ at weight $1$, with
the canonical realisation $\pi = \frac{1}{2i}\oint_{|z|=1}\frac{dz}{z}$
(or equivalently $\pi = \int_{-1}^1 dx/\sqrt{1-x^2}$).  Its inverse
$1/\pi$ is NOT a period in the strict Kontsevich–Zagier sense
(periods are closed under multiplication but not division), but
sits in the localisation $\mathbb{Q}[\pi, \pi^{-1}]$ which is the
natural ambient algebra for ratios of periods.

The Tate weight grading of $\mathbb{Q}[\pi, \pi^{-1}]$ assigns
weight $n$ to $\pi^n$ for $n \in \mathbb{Z}$.  All M1 outputs sit
on a single Tate weight ray indexed by the integer power of $\pi$.
This is the simplest mixed-Tate sub-structure:\ pure Tate at depth
$0$, weight in $\mathbb{Z}$.

## GeoVac M1 enumeration

Every closed form witnessed across the corpus that engages the
Hopf-base measure mechanism, with its M1-content explicit:

| Witness | M1 factor | Source |
|:---|:---:|:---|
| Hopf base measure on $S^3 \to S^2$ Hopf bundle | $\mathrm{Vol}(S^2)/4 = \pi$ | Paper 25 §II.1 |
| Paper 38 L2 asymptotic propinquity rate constant on $SU(2)$ | $4/\pi = \mathrm{Vol}(S^2)/\pi^2$ | Paper 38 §L2 |
| Paper 40 universal rate $4/\pi$ on all compact Lie groups | $4/\pi$ | Paper 40 Thm 1 |
| Paper 43 §10.2 Pythagorean orthogonality prefactor | $1/\pi^2$ | Paper 43 `cor:pythagorean_orthogonality` |
| Paper 50 §3 F-theorem coefficients (M1×M3 cross-product) | $1/\pi^2$ in $-3\zeta(3)/(16\pi^2)$ | Paper 50 Thm 3.4 |
| K = $\pi(B + F - \Delta)$ Paper 2 combination rule | explicit $\pi$ factor | Paper 2 §IV; Paper 32 §VIII `rem:K_under_theorem` |

Every entry is $q \cdot \pi^n$ for $q \in \mathbb{Q}$ and $n \in \mathbb{Z}$
with $n \in \{-2, -1, 0, 1\}$ on the empirical witnesses logged so far.
No M1 closed form has produced $\pi^n$ for $|n| \ge 3$ at depth $0$;
higher powers of $\pi$ appear only through M2 (Seeley–DeWitt,
$\pi^{2k}$ at all $k$) or M3 (vertex-parity sums, $\pi^4$ etc.).

## Arithmetic classification

Trivial:\ every M1 closed form lies in $\mathbb{Q}[\pi, \pi^{-1}]$ at
Tate weight equal to the integer power of $\pi$.  The classification
is **convention-free** at the M1 level — unlike M2's
volume-normalization-dependent $\sqrt\pi$ vs $\pi^2$ structure (Sprint
Mixed-Tate Test §"raw vs volume-normalized"), M1 produces $\pi$
directly with no choice of dimensional rescaling.  This is because
the Hopf-base measure is itself a volume integral, not a heat-kernel
coefficient at $t \to 0^+$.

## Cross-check vs M2 and M3 sprints

The morning's three sprints partition the master Mellin engine:
- **M1** (this sprint) $\subset \mathbb{Q}[\pi, \pi^{-1}]$.  Pure Tate weight $\in \mathbb{Z}$, depth $0$.
- **M2** (Sprint Mixed-Tate Test) $\subset \bigoplus_k \pi^{2k}\cdot\mathbb{Q}$.  Pure Tate, even weight only, depth $0$.
- **M3** (Sprint M3 Cyclotomic Mixed-Tate) $\subset \mathcal{MT}(\mathbb{Z}[i, 1/2])$.  Mixed Tate at level $\le 4$, depth $\ge 1$.

Containment structure:\
$$M_2 \subset M_1 \cdot \mathbb{Q} \subset \mathcal{MT}(\mathbb{Z}),$$
$$\text{but } M_1 \not\subset M_2 \text{ (M2 has only even weights)},$$
$$M_3 \not\subset M_1 \text{ (M3 has depth} \ge 1 \text{)},$$
$$M_3 \cap M_1 = \{q \pi^n : n \in \mathbb{Z}\} \cap \mathcal{MT}(\mathbb{Z}[i, 1/2]) = \mathbb{Q}[\pi, \pi^{-1}].$$

The three sub-mechanisms are role-disjoint at the Mellin-engine
index $k \in \{0, 1, 2\}$ but produce overlapping period rings;
joint engagement (Paper 32 §VIII `rem:joint_engagement`) produces
products like $\mathrm{Vol}(S^2) \cdot a_0^{D^2}/\zeta_{D^2}(s)$
that pick up factors from multiple mechanisms.  The K formula
$K = \pi (B + F - \Delta)$ is exactly the canonical example:\ M1's
explicit $\pi$, M2's $\pi^2$ inside $F = \zeta(2)$ via the Fock
Dirichlet evaluation, and M3's $\Delta = 1/40$ Dirac mode count
(rational, weight-$0$ in M3).

## Verdict

**POSITIVE.**  M1 closure is bookkeeping; the substantive content
was already in the Paper 32 §VIII case-exhaustion theorem proof.
This sprint adds a formal Corollary so the master Mellin engine
classification reads symmetrically across M1, M2, M3.

**Honest scope.**

1. The M1 ring $\mathbb{Q}[\pi, \pi^{-1}]$ is the **localised** pure
   Tate algebra at $\pi$.  In strict Kontsevich–Zagier "periods are
   integrals of algebraic functions over algebraic domains," only
   $\pi^n$ for $n \ge 0$ are periods; $1/\pi$ is an inverse period.
   M1 witnesses include $1/\pi^2$ (Paper 43 Pythagorean
   orthogonality) and $4/\pi$ (Paper 38 L2 rate), which are inverse
   periods.  The localisation is the natural ambient algebra; the
   strict-periods statement is the sub-monoid generated by $\pi$.

2. M1 trivially partitions GeoVac's $\pi$-content.  The non-trivial
   structural content is the case-exhaustion theorem itself
   (Paper 32 §VIII):\ EVERY $\pi$ appearing in any GeoVac chain
   engages M1 or M2 or M3 (proven 208/208 in Paper 35 Prediction 1
   panel; corollary in `cor:m2_mixed_tate` + this sprint's
   `cor:m1_pure_tate` + `cor:m3_cyclotomic_mixed_tate`).

## Paper-level implications (proposals)

1. **Paper 32 §VIII**: add Corollary `cor:m1_pure_tate` immediately
   before `cor:m2_mixed_tate` (M1 should come first in $k$-index
   order).

2. **Paper 18 §III.7**: extend the existing M1 paragraph (around
   line 1026 of the current `paper_18_exchange_constants.tex`,
   "The master Mellin engine's $M_1$ Hopf-base measure mechanism
   produces the universal $4/\pi$ rate constant ...") with the
   period-ring statement.

3. **CLAUDE.md §2 entry**: brief one-liner.

## Open questions / follow-on

None at the M1 corollary level — M1 is structurally trivial.  The
genuine forward work is the broader Paper 18 §III.7 major rewrite
(separately tasked) that recasts the entire taxonomy as
periods-of-motives stratification with M1/M2/M3 as the three
period-ring spines.

## Files used

- `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
  Theorem~\ref{thm:pi_source_case_exhaustion} (lines 1620–1707;
  M1's "$\mathrm{Vol}(S^2)/4 = \pi$ identity" on line 1703).
- `papers/group3_foundations/paper_18_exchange_constants.tex` §III.7
  M1 paragraph (lines 1024–1040 in the current state, after this
  morning's M2/M3 edits).
- `debug/sprint_mixed_tate_test_memo.md` and
  `debug/sprint_m3_cyclotomic_mixed_tate_memo.md` (companion
  sprints from this morning).
