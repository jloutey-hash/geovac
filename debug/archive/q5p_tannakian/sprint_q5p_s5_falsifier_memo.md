# Sprint Q5'-S5-falsifier — Master Mellin engine on the truncated Bargmann-Segal Hardy-on-$S^5$ triple

Date: 2026-06-05
Scope: 1-day falsifier sprint of the Q5'-CH-1 chirality-parity finding
(2026-06-05, `debug/sprint_q5p_ch1_memo.md`).  Computes the master
Mellin engine moments $\mathrm{Tr}(L^k\,e^{-tL^2})$ and trial supertraces
$\mathrm{Tr}(\chi\,L^k\,e^{-tL^2})$ at $k \in \{0, 1, 2\}$ on the truncated
Bargmann-Segal Hardy-on-$S^5$ spectral triple at
$N_{\max} \in \{2, 3, 4\}$, using exact `sympy.Rational` arithmetic.

Driver: `debug/compute_s5_k_nmax_truncated.py`. Data:
`debug/data/sprint_q5p_s5_falsifier_data.json`. Wall time: ~0.3 s total
across the three cutoffs.  No floating point; no PSLQ.

## TL;DR

**Verdict: BREAKS-AS-EXPECTED.** The CH-1 chirality-parity selection
rule fails on the Bargmann-Segal Hardy-on-$S^5$ triple via *both*
mechanisms simultaneously, in the way Paper 24 §V.4 predicts:

| Mechanism | CH on $S^3$ | Hardy on $S^5$ |
|:----------|:------------|:----------------|
| Spectral $\pm$-symmetry | $\{\pm(n+1/2)\}$ symmetric | $\{N+3/2\}$ strictly positive |
| Chirality alignment $\gamma_{ii} = \chi_i$ | exact on Camporesi-Higuchi | no spinor lift on $(N,0)$ tower |

Bit-exact data: $\mathrm{Tr}(L^j)$ is NONZERO at every odd $j$
($j = 1, 3, 5, \ldots$, e.g. $\mathrm{Tr}(L^1) = 30$ at $N_{\max} = 2$,
$\mathrm{Tr}(L^1) = 75$ at $N_{\max} = 3$, $\mathrm{Tr}(L^1) = 315/2$ at
$N_{\max} = 4$) — the spectral-symmetry leg of the CH-1 argument is
absent on the Hardy sector exactly because the Bargmann transform
truncates to the holomorphic ($N \ge 0$) sector.

Both natural trial gradings $\chi_N(v) = (-1)^{N_v}$ and
$\chi_l(v) = (-1)^{l_v}$ produce *nonzero* even-$j$ supertraces
(e.g. $\mathrm{Tr}(\chi_N\,L^0) = 4 \ne 0$ at $N_{\max} = 2$).  The trial
gradings also collapse to one another: $\chi_N = \chi_l$ on every node
because on shell $N$ the orbital angular momentum has the same parity
as $N$ (the SU(3) symmetric branching rule, Paper 24 Eq. nearby
$\hat{N}$-shell decomposition).  This collapse is itself structural
evidence that the Hardy sector does not carry an independent
$\mathbb{Z}/2$ chirality grading orthogonal to shell parity.

The CH-1 finding is therefore **CH-specific**, and Stage 1 of Q5'
(operator-order grading enrichment $\omega^{\mathrm{tri}}$ on
$\mathrm{dg}(\mathcal{T}_{n_{\max}})$) correctly targets the CH triple.

## Verdict-against-gate

| Gate | Outcome |
|:-----|:--------|
| BREAKS-AS-EXPECTED | **THIS.** Spectral $\pm$-symmetry fails (spectrum positive) AND no natural trial chirality reproduces the CH-1 even-$j$ supertrace zeros.  Both structural legs of the CH-1 argument are absent. |
| BREAKS-DIFFERENTLY | Not observed. |
| HOLDS-UNEXPECTEDLY | Not observed.  No trial grading produces the rule. |
| BLOCKED | Not observed.  The Bargmann graph is structurally similar enough (rational arithmetic, diagonal spectrum-computing operator, dipole adjacency) that the computation runs at the same exact-rational grade as CH-1. |

## Computational summary

**$N_{\max} = 2$, $\dim\,\mathcal{H} = 10$ (= $1 + 3 + 6$).**

| $j$ | $\mathrm{Tr}(L^j)$ | $\mathrm{Tr}(\chi_N L^j)$ | $\mathrm{Tr}(\chi_l L^j)$ |
|:---:|:------------------:|:--------------------------:|:--------------------------:|
| 0   | 10                 | 4                          | 4                          |
| 1   | 30                 | 15                         | 15                         |
| 2   | 189/2              | 57                         | 57                         |
| 3   | 615/2              | 855/4                      | 855/4                      |
| 4   | 8181/8             | 3153/4                     | 3153/4                     |
| 5   | 27615/8            | 45855/16                   | 45855/16                   |
| 6   | 376749/32          | 164937/16                  | 164937/16                  |

**$N_{\max} = 3$, $\dim\,\mathcal{H} = 20$ (= $1+3+6+10$).**

| $j$ | $\mathrm{Tr}(L^j)$ | $\mathrm{Tr}(\chi_N L^j)$ |
|:---:|:------------------:|:--------------------------:|
| 0   | 20                 | $-6$                       |
| 1   | 75                 | $-30$                      |
| 2   | 297                | $-291/2$                   |
| 3   | 4875/4             | $-1395/2$                  |
| 4   | 20493/4            | $-26499/8$                 |
| 6   | 1516977/16         | $-2327331/32$              |

**$N_{\max} = 4$, $\dim\,\mathcal{H} = 35$.**

| $j$ | $\mathrm{Tr}(L^j)$ | $\mathrm{Tr}(\chi_N L^j)$ |
|:---:|:------------------:|:--------------------------:|
| 0   | 35                 | 9                          |
| 1   | 315/2              | 105/2                      |
| 6   | 32641323/64        | 21918753/64                |

Every entry is positive rational (or the sign matches the chirality
balance at the leading moment); none is identically zero at any $j$
across $N_{\max} \in \{2, 3, 4\}$.

## Structural argument for why the result is what it is

Both CH-1 facts that forced the rule on $S^3$ are *structurally*
absent on the Hardy sector:

1. **Spectral $\pm$-symmetry fails.** The Bargmann-Segal construction
   restricts to the *holomorphic* (Hardy) sector $H^2(S^5) =
   \bigoplus_{N \ge 0} \mathcal{H}_N$.  The Euler operator
   $\hat{N} + 3/2$ has eigenvalue $N + 3/2 > 0$ on every shell; there
   is no anti-holomorphic ($N < 0$) partner with eigenvalue
   $-(N + 3/2)$, because the Bargmann transform projects $L^2(\mathbb{R}^3)$
   to *only* the holomorphic sector (Paper 24 Eq. (1)).  This is not
   a truncation artifact:  the Hardy sector is the entire physical
   Hilbert space of the Bargmann construction, and it has been
   $\pm$-asymmetric by construction since Bargmann 1961.  The
   corresponding even powers of $L$ ($L^0, L^2, L^4, \ldots$) have
   trace $\sum_v (N_v + 3/2)^j > 0$; the odd powers similarly are
   positive sums of positive terms.  *No cancellation across $\pm$ is
   available.*

2. **No spinor lift, no canonical $\gamma$.** Paper 24 §V.4 states the
   block-level obstruction: *"the symmetric irreducible representations
   $(N, 0)$ of $\mathrm{SU}(3)$ are bosonic (integer $m_l$); the Hardy
   sector has no half-integer weight content, hence no spinor lift
   compatible with the Bargmann transform."*  Without a spinor sector
   there is no canonical $\mathbb{Z}/2$ grading $\gamma$ that would be
   the Hardy-sector analog of CH chirality.  Two natural trial
   gradings — shell-parity $\chi_N(v) = (-1)^{N_v}$ and angular-parity
   $\chi_l(v) = (-1)^{l_v}$ — are the only obvious candidates from the
   $(N, l, m_l)$ node labels, and *both collapse to each other*:
   on shell $N$, $l$ has parity $N$ (SU(3) symmetric branching), so
   $\chi_N = \chi_l$ everywhere.  Neither produces the even-$j$ zero
   pattern of the CH-1 rule.

The structural fingerprint of the CH-1 rule was the *coincidence* of
$\pm$-symmetric spectrum and chirality alignment.  The Hardy sector has
neither.  The CH-1 rule is therefore CH-specific:  it depends
load-bearingly on the Camporesi-Higuchi spinor bundle structure with
half-integer eigenvalues $\pm(n+1/2)$ and chirality alignment
$\gamma_{ii} = \chi_i$.

## Implications for Stage 1 and CH-1 scope-of-claim

**Stage 1 ($\omega^{\mathrm{tri}}$ on the truncated CH triple).** This
sprint confirms that the candidate enrichment $\omega^{\mathrm{tri}}$
proposed in CH-1's Stage 1 program lives naturally on the CH triple,
not on a generic discrete spectral triple in the framework's library.
The Bargmann-Segal triple, despite sharing the $\pi$-free rational
skeleton, *cannot* serve as a substrate for the same enrichment because
the operator-order/chirality factorization that CH-1 found at the
bit-exact level requires both spectral $\pm$-symmetry and chirality
alignment, and the Hardy sector provides neither.

**Scope sharpening for CH-1.**  The CH-1 memo's caveat 5 said:

> The robustness test would be to repeat on the Bargmann-Segal
> Hardy-on-$S^5$ triple (Paper 24) — where the spectrum is *not*
> symmetric (only positive eigenvalues) and the parity rule should
> *break*.  This is a clean follow-on test.

The current sprint executes that test bit-exactly and confirms the
prediction.  The CH-1 finding is the operator-system finite-cutoff
witness of a property of the Camporesi-Higuchi Dirac on $S^3$ — *not*
of generic spectral triples in the framework's library, and *not* of
all $\pi$-free discrete spectral triples.

**Cross-reference to Paper 24's four-layer asymmetry.** This sprint is
the master-Mellin-engine fifth slice of the Coulomb/HO asymmetry
already documented at four layers in Paper 24 §V:
(i) spectrum-computing $L_0$, (ii) calibration $\pi$, (iii) non-abelian
Wilson gauge with natural matter coupling, (iv) modular-Hamiltonian
structure of the wedge KMS state.  The fifth layer surfaced here is:

> **(v) Chirality-parity selection rule of the master Mellin engine.**
> The CH triple on $S^3$ exhibits a bit-exact factorization of the
> Mellin source into $\mathrm{Z}/2$ chirality $\times$ heat-kernel
> order; the Hardy triple on $S^5$ does not, because spectral
> $\pm$-symmetry fails AND the chirality grading is absent.

Whether to promote this to a formal layer-5 in Paper 24 §V is a
PI decision; the structural content is consistent with the four
existing layers and is already implicit in Paper 24 §V.4's spinor /
half-integer-wedge obstructions.  *Recommendation only*; no Paper 24
edits applied in this sprint.

## Curve-fit-audit (per `feedback_audit_numerical_claims`)

The claim "the chirality-parity rule breaks on the Bargmann-Segal Hardy
sector" has zero free parameters (no fitting; everything is exact
rational at finite $N_{\max}$).  Selection bias: the Hardy sector was
named in CH-1's caveats 5 as the structural-falsifier candidate, so
this sprint executes a pre-registered test, not a fishing expedition.
Alternative readings:

1. *Alternative 1: pick a non-natural grading $\chi$ that aligns.*
   One could attempt to construct an artificial $\chi$ that produces
   the even-$j$ zeros, e.g. by hand-tuning per-node signs.  This would
   require $\sum_v \chi_v\,(N_v + 3/2)^j = 0$ for every even $j$ —
   infinitely many polynomial constraints on a finite-dimensional sign
   assignment.  At $N_{\max} = 2$ (dim 10) this is an overdetermined
   linear system with no nontrivial rational solution; numerically
   running a least-squares solve and rounding would NOT yield exact
   rational zeros.  We did not formally verify this alternative is
   ruled out, but the polynomial-overdetermination argument is
   structurally tight.
2. *Alternative 2: Bargmann graph at higher $N_{\max}$ recovers the rule.*
   The mechanism is independent of $N_{\max}$: spectrum positivity is
   a structural property of the Hardy sector at every cutoff, and the
   $\chi_N = \chi_l$ collapse follows from $(N, 0)$ branching at every
   $N$.  Higher cutoffs would not change the verdict.
3. *Alternative 3: a different operator on the same node set might
   restore the rule.* One could ask whether a Dirac-like construction
   that is *not* the Euler operator (e.g., an off-diagonal `square-root`
   of the Hodge Laplacian) admits chirality alignment.  This is what
   Paper 24 §V.4 explicitly closes off — no spinor lift exists on the
   $(N, 0)$ tower compatible with the Bargmann transform.  It would
   require a categorically different spectral-triple construction on
   $S^5$, not a modification of the Bargmann-Segal one.

The robustness test (alternative reading 4) is to repeat on a third
discrete spectral triple: e.g., the Hopf-graph U(1) gauge structure
(Paper 25) which is bosonic, or the SU(2) Wilson lattice (Paper 30)
which has half-integer link-rep but no canonical chirality.  Both
would also produce BREAKS-AS-EXPECTED for the same load-bearing
reasons (Paper 25 has Hodge-1 $L_1 = B^T B$ as a graph Laplacian, not
a first-order Dirac; Paper 30 has SU(2) links but no spinor matter
sector at the level studied).  This is consistent with the CH-1
finding being CH-specific.

## Honest scope

1. **No periods at finite $N_{\max}$.** Every quantity computed here
   is bit-exact rational, just as in CH-1.  The M1/M2/M3 period
   content (which lives in the Mellin transform $\int_0^\infty t^{s-1}
   \mathrm{Tr}(L^k\,e^{-tL^2})\,dt$ via $\Gamma(s)$ and zeta values)
   enters only at the $N_{\max} \to \infty$ continuum limit, which
   this sprint does not compute.  The finite-cutoff data establishes
   the structural fingerprint, not the asymptotic content.
2. **Heat-kernel-order vs. chirality-grading distinction.** This
   sprint did not attempt to verify the M1/M2/M3 partition itself on
   the Bargmann graph — only the chirality-parity selection rule.  A
   separate sprint would be needed to compute the small-$t$ Mellin
   coefficient ring on the Bargmann graph and compare to Paper 32
   §VIII's case-exhaustion theorem.  Paper 24 §V.4 expects pure
   $\pi^{2k}\cdot\mathbb{Q}$ (M2 only, no M3 since no vertex-parity
   half-integer Hurwitz content), which would be the natural Q3'
   target (Paper 55 §7.4 — already opened by Sprint A6 closure for the
   S⁵ Casimir trace).
3. **Trial gradings tested are not exhaustive.** Only two natural Z/2
   gradings ($\chi_N$, $\chi_l$) were tested.  An exhaustive search
   over all Z/2 gradings on the Bargmann graph is computationally
   tractable at $N_{\max} = 2$ (2^10 = 1024 candidates) but irrelevant
   to the structural verdict: the polynomial-overdetermination
   argument in Alternative 1 above rules out *any* nontrivial $\chi$
   from making all even-$j$ supertraces vanish identically.
4. **Curve-fit-audit / diagnostic-only.** No structural claim about
   the Bargmann graph is promoted to a paper or to CLAUDE.md §2 by
   this sprint.  The verdict is *negative diagnostic*: CH-1 is
   CH-specific.  Recommendations for Paper 24 §V and Paper 32 §VIII
   are below as PI direction items, not applied.

## Files produced

- `debug/compute_s5_k_nmax_truncated.py` — driver (~330 lines).
- `debug/data/sprint_q5p_s5_falsifier_data.json` — exact rational
  data per $N_{\max} \in \{2, 3, 4\}$, including full-D moments.
- `debug/sprint_q5p_s5_falsifier_memo.md` — this memo.

## Recommended paper edits (PI to apply, decline, or modify)

No paper edits applied.  Recommendations:

### Paper 24 §V (Coulomb/HO asymmetry) — optional fifth layer

Consider adding a one-paragraph fifth-layer description after the
existing modular-Hamiltonian fourth layer:

> *Fifth layer: chirality-parity selection rule of the master Mellin
> engine (Sprint Q5'-S5-falsifier, June 2026; memo
> `debug/sprint_q5p_s5_falsifier_memo.md`).*  On the Camporesi-Higuchi
> $S^3$ spectral triple at every finite $n_{\max}$, the Mellin source
> $\mathrm{Tr}(\Lambda^k\,e^{-t\Lambda^2})$ exhibits a bit-exact
> $\mathrm{Z}/2$ chirality-parity factorization (Sprint Q5'-CH-1).  On
> the Bargmann-Segal Hardy sector at every finite $N_{\max}$ this rule
> fails, because the Hardy spectrum is strictly positive and no
> half-integer-weight grading exists on the $(N, 0)$ tower.  This
> propagates the Coulomb/HO asymmetry from the modular-Hamiltonian
> level to the master Mellin engine level.

### Paper 32 §VIII (master Mellin engine + case-exhaustion theorem) — optional remark

Consider an explanatory remark sharpening the CH-specificity of the
chirality-parity finding:

> *The bit-exact chirality-parity factorization of the master Mellin
> source on the truncated CH triple (Sprint Q5'-CH-1) is CH-specific:
> it requires spectral $\pm$-symmetry and chirality alignment of the
> Camporesi-Higuchi Dirac.  The Bargmann-Segal Hardy-on-$S^5$ triple
> at finite $N_{\max}$ has neither (Sprint Q5'-S5-falsifier).  The
> $\omega^{\mathrm{tri}}$ enrichment program (Q5'-CH-1 Stage 1)
> correctly targets the CH triple as substrate.*

(Both recommendations are optional; no paper edits applied this sprint.)

## One-line verdict

BREAKS-AS-EXPECTED — Bargmann-Segal Hardy-on-$S^5$ has neither spectral
$\pm$-symmetry nor a natural chirality grading, so the CH-1 master-Mellin
parity rule is CH-specific and Stage 1 of Q5' correctly targets the CH
triple.
