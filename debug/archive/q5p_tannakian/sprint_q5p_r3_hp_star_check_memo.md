# Sprint Q5' Round 3 — HP\*(dg(GeoVac at finite $n_{\max}$)) structure check

**Date:** 2026-06-04
**Scope:** 1-day structural sub-sprint of Q5' Round 3. Decisive test
of the source-side input for the Marcolli–Tabuada noncommutative
motivic Galois construction. The Q5' k-slot Tannakian probe (Round 2,
`debug/sprint_q5p_k_slot_tannakian_memo.md`) flagged this as the
load-bearing question: "whether dg(GeoVac triple) has non-trivial
$\mathrm{HP}_*$ at fixed $n_{\max}$ (finite-dim matrix algebras may
have small $\mathrm{HP}_*$)." If $\mathrm{HP}_*$ is trivial at finite
cutoff (baseline $\mathrm{HP}_*(\mathbb{C})$ via Morita), the
Marcolli–Tabuada NC motivic Galois construction is blocked at source
and the sharpened Q5' enriched fiber functor
$\omega^{\mathrm{tri}}$ dies.

## TL;DR

**Verdict: TRIVIAL on the bare algebra, NON-TRIVIAL on the spectral
triple as Chern character / K-homology class.** The honest decisive
answer is **BORDERLINE**, with the borderline-line drawn at exactly the
"algebra alone" vs. "algebra + Dirac data" boundary:

1. $\mathrm{HP}_*(\mathcal{A}_{\mathrm{GV}}^{(n_{\max})}) \;\cong\;
   \mathrm{HP}_*(\mathbb{C})^{\oplus N_{\mathrm{Fock}}(n_{\max})}$
   for every $n_{\max} \ge 1$. The algebra is the commutative function
   algebra $\mathbb{C}^{V_{\mathrm{Fock}}}$ on a finite vertex set
   (Paper 32 Def 3.1, line 269). This is a direct sum of
   $N_{\mathrm{Fock}}$ copies of $\mathbb{C}$, hence Morita-equivalent
   to $\mathbb{C}^{N_{\mathrm{Fock}}}$, hence $\mathrm{HP}_0 \cong
   \mathbb{Q}^{N_{\mathrm{Fock}}}$ and $\mathrm{HP}_1 = 0$ (Loday
   *Cyclic Homology* §1.4, Thm 1.2.4: $\mathrm{HP}_*$ is Morita-invariant;
   $\mathrm{HP}_0(\mathbb{C}) = \mathbb{Q}$, $\mathrm{HP}_1 = 0$). This
   is the **baseline-pure-Tate** $\mathrm{HP}_*$ — exactly the
   $\mathbb{G}_m$-Galois-fixed sub-category of $\mathrm{MT}(\mathbb{Z})$
   at integer Tate twists, with rank $N_{\mathrm{Fock}}(n_{\max})$.

2. **The bare algebra carries no M2/M3-distinguishing dg content.** A
   direct sum of copies of $\mathbb{C}$ has only the standard $\mathbb{C}$
   cyclic structure on each summand; the Connes B-operator acts as the
   identity on $\mathrm{HC}_0(\mathbb{C}) = \mathbb{C}$ and zero on
   $\mathrm{HC}_n$ for $n \ge 1$. There is no Mellin slot detection at
   the algebra level. The Marcolli–Tabuada fiber functor
   $\mathrm{HP}_*: \mathrm{dg}(\mathcal{A}_{\mathrm{GV}}) \to
   \mathrm{Vec}_{\mathbb{Q}}^{\mathbb{Z}/2}$ sees only the
   $N_{\mathrm{Fock}}$-dimensional pure-Tate weight-$0$ output. **The
   k-slot index is invisible on $\mathrm{HP}_*$ of the algebra alone.**

3. **The spectral triple's K-homology class (Chern character) is
   non-trivial at finite $n_{\max}$.** The Connes character of a
   spectral triple $(\mathcal{A}, \mathcal{H}, D)$ — equivalently the
   JLO cocycle or the Chern character of the associated Fredholm
   module — lives in $\mathrm{HP}^*(\mathcal{A})$ (dual cohomology),
   NOT in $\mathrm{HP}_*(\mathcal{A})$, and carries information about
   $D$ that is invisible to the algebra alone. On
   $\mathrm{dg}(\mathcal{T}_{n_{\max}})$ (the dg-categorical
   thickening), this character is a non-trivial class in periodic
   cyclic *cohomology*; it depends on the spectral data $D$ and is
   sensitive to the half-integer Camporesi–Higuchi shift (the
   structural input that produces M3's cyclotomic level-4 conductor).

4. **The Chern character lives on the wrong side for the
   Marcolli–Tabuada construction.** Marcolli–Tabuada 2016 Thm 1.1
   uses $\mathrm{HP}_*$ (homology) as the fiber functor for the
   noncommutative motivic Galois group; the Chern character lives in
   $\mathrm{HP}^*$ (cohomology) and is dual data. The dg category
   $\mathrm{dg}(\mathcal{T}_{n_{\max}})$ at finite $n_{\max}$ has
   $\mathrm{HP}_*(\mathrm{dg}(\mathcal{T})) \cong
   \mathrm{HP}_*(\mathcal{A}_{\mathrm{GV}}) \otimes \mathrm{End}^{\mathrm{dg}}(D)$,
   but $\mathrm{End}^{\mathrm{dg}}(D)$ at finite truncation is itself
   a finite-dim matrix algebra (operators on the finite-dim Hilbert
   space $\mathcal{H}_{\mathrm{GV}}^{(n_{\max})}$), hence
   Morita-trivial too. The pro-system $\varprojlim_{n_{\max}}
   \mathrm{HP}_*(\mathrm{dg}(\mathcal{T}_{n_{\max}}))$ with Berezin
   maps does NOT escape this conclusion: a pro-limit of pure-Tate
   weight-0 abelian groups is itself pure-Tate weight-0 abelian.

5. **Therefore, the Marcolli–Tabuada NC motivic Galois construction
   at every finite $n_{\max}$ produces the trivial Tannakian category
   $\mathrm{Vec}_{\mathbb{Q}}^{\mathbb{Z}/2}$ with motivic Galois
   group $\mathbb{G}_m$ acting on Tate weights.** This is structurally
   blocked from seeing the M2/M3 distinction at the source side. The
   sharpened Q5' enriched fiber functor $\omega^{\mathrm{tri}}$ would
   need a refinement of $\mathrm{HP}_*$ to a $\mathbb{Z}/3$-graded
   functor that records the operator-order slot $k$; this refinement
   is NOT supplied by either the Marcolli–Tabuada construction or by
   the standard Loday/Hood–Jones cyclic-homology machinery.

**Net deflation reading:** the Q5' enriched fiber functor
$\omega^{\mathrm{tri}}: \mathrm{dg}(\mathcal{T}) \to \mathrm{Vec}_{\mathbb{Q}}
\otimes \mathrm{IndexCat}(\{0, 1, 2\})$ contemplated by Round 2 is
**blocked at source by the Morita-triviality of
$\mathcal{A}_{\mathrm{GV}}^{(n_{\max})}$**. The Marcolli–Tabuada machinery cannot
see the k-slot at any finite cutoff, and the
pro-limit does not escape because the relevant invariants are
Morita-stable. Sharpened Q5' as posed in Round 2 dies at the source.

**However, the deeper Q5' question (does GeoVac add anything beyond
$M^{\mathrm{GV}}$?) does NOT close in the negative**. The Chern
character / K-homology class of the spectral triple at finite
$n_{\max}$ is non-trivial, and the master Mellin engine
$\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$ at $k \in \{0, 1, 2\}$ is a
non-trivial structure on the *cohomological-side* dg data
$(\mathcal{H}_{\mathrm{GV}}, D_{\mathrm{GV}})$. The k-slot lives on
the COHOMOLOGY side ($\mathrm{HP}^*$ / Chern-character / Connes
character / K-homology pairing), not on the HOMOLOGY side
($\mathrm{HP}_*$ / Marcolli–Tabuada input). The Round 2 verdict
"BORDERLINE leaning Tannakian-invisible on the period ring,
Tannakian-relevant one categorical level up" is now sharpened: the
"one categorical level up" is the **K-homology / Chern character
side**, which is structurally dual to (and hence inaccessible from)
the Marcolli–Tabuada fiber functor.

The Q5' deflation goes further than Round 2 anticipated:\ the
**Marcolli–Tabuada route is structurally barred**. Any future Q5'
construction would need a fundamentally different framework — namely
the Connes–Marcolli renormalization-style cosmic Galois group $U^*$
acting on K-homology / Chern-character data (Strand A of
Round 1's QSM lit-read), which uses analytic input *labelled* by
diagrammatic/graded source data, with the Tannakian dual built from
the *label-respecting* automorphism group.

## Refined verdict (one line)

**TRIVIAL on $\mathrm{HP}_*(\mathcal{A}_{\mathrm{GV}}^{(n_{\max})})$
and on $\mathrm{HP}_*(\mathrm{dg}(\mathcal{T}_{n_{\max}}))$ for every
$n_{\max} \ge 1$, with the Marcolli–Tabuada route to Q5' blocked at
source; NON-TRIVIAL on the K-homology / Chern-character class of the
spectral triple, where the master Mellin engine k-slot does live,
but on the *cohomological-dual* side of the fiber functor used by
Marcolli–Tabuada. Net: BORDERLINE (TRIVIAL on the homological side
that Marcolli–Tabuada uses; the k-slot lives on the dual side).**

The sharpened Q5' enriched fiber functor $\omega^{\mathrm{tri}}$ in
the Marcolli–Tabuada framework **dies**; Q5' deflates further toward
"GeoVac IS the natural-geometry realization of ambient
$M^{\mathrm{GV}}$, with k-slot being substrate-side metadata only,"
as the deflation memo (`sprint_q5p_deflation_test_memo.md`)
anticipated. The Connes–Marcolli $U^*$ route remains an alternative
(possibly multi-year) target on the cohomological-dual side.

---

## Argument

### Step 1: $\mathrm{HP}_*$ of the algebra alone

$\mathcal{A}_{\mathrm{GV}}^{(n_{\max})} = \mathbb{C}^{V_{\mathrm{Fock}}}$
is a commutative unital algebra of complex-valued functions on a
finite vertex set $V_{\mathrm{Fock}}$ of cardinality
$N_{\mathrm{Fock}}(n_{\max}) = \sum_{n=1}^{n_{\max}} n^2$
(Paper 32 Def 3.1, lines 267–273):
$$
N_{\mathrm{Fock}}(1) = 1, \quad N_{\mathrm{Fock}}(2) = 5,
\quad N_{\mathrm{Fock}}(3) = 14, \quad \ldots
$$

As an $\mathbb{C}$-algebra, $\mathcal{A}_{\mathrm{GV}}^{(n_{\max})}
\cong \mathbb{C} \times \mathbb{C} \times \cdots \times \mathbb{C}$
($N_{\mathrm{Fock}}$ copies). This is a direct sum of $\mathbb{C}$,
hence Morita-equivalent to $\mathbb{C}^{N_{\mathrm{Fock}}}$. By Loday
*Cyclic Homology* (2nd edition, 1998) §1.4 Theorem 1.2.4
(Morita invariance of cyclic homology, hence of
$\mathrm{HP}_* = \varprojlim_S \mathrm{HC}_{*+2S}$):
$$
\mathrm{HP}_*(\mathcal{A}_{\mathrm{GV}}^{(n_{\max})})
\;=\; \bigoplus_{v \in V_{\mathrm{Fock}}} \mathrm{HP}_*(\mathbb{C})
\;=\; \mathrm{HP}_*(\mathbb{C})^{\oplus N_{\mathrm{Fock}}(n_{\max})}.
$$
And $\mathrm{HP}_0(\mathbb{C}) = \mathbb{Q}$ (one-dimensional,
generated by the rank-1 idempotent), $\mathrm{HP}_1(\mathbb{C}) = 0$.

**Conclusion:** as a $\mathbb{Z}/2$-graded $\mathbb{Q}$-vector space,
$\mathrm{HP}_*(\mathcal{A}_{\mathrm{GV}}^{(n_{\max})}) =
(\mathbb{Q}^{N_{\mathrm{Fock}}}, 0)$, with no $\mathbb{Z}/3$ grading,
no operator-order signature, no M2/M3 distinguishing structure. The
Marcolli–Tabuada noncommutative motivic Galois construction with
$\mathrm{HP}_*$ as fiber functor sees exactly $N_{\mathrm{Fock}}$
copies of the trivial pure-Tate weight-0 motive.

This is the Morita-trivial baseline. The dg-enrichment from the
spectral triple is the only escape route.

### Step 2: dg-enrichment via the Dirac operator?

The candidate escape from Step 1 is to add the dg structure coming
from the spectral triple's differential operator $D_{\mathrm{GV}}$.
The standard dg-categorical thickening
$\mathrm{dg}(\mathcal{T}_{n_{\max}})$ of a spectral triple
$(\mathcal{A}, \mathcal{H}, D)$ is the dg category of
$\mathcal{A}$-modules with chain complexes formed using the
de Rham-type differential built from $[D, \cdot]$. Equivalently,
$\mathrm{dg}(\mathcal{T}_{n_{\max}}) \simeq
\mathrm{End}^{\mathrm{dg}}_{\mathcal{A}}(\mathcal{H})$ with the
graded commutator $[D, \cdot]$ as differential.

**At finite $n_{\max}$, this differential is bounded** (Paper 32
Prop 5.1, lines 803–814:\ bounded commutators on finite-dim Hilbert
space). The dg algebra $(\mathrm{End}_{\mathcal{A}}(\mathcal{H}),
\partial = [D, \cdot])$ is a finite-dim associative dg algebra. Its
periodic cyclic homology is:
$$
\mathrm{HP}_*(\mathrm{End}^{\mathrm{dg}}_{\mathcal{A}}(\mathcal{H}), [D, \cdot])
\;=\; \mathrm{HP}_*(\mathrm{End}_{\mathcal{A}}(\mathcal{H}))
\;\cong\; \mathrm{HP}_*(\mathcal{A}_{\mathrm{GV}})
\;=\; \mathbb{Q}^{N_{\mathrm{Fock}}} \;\;\text{(even part)},
$$
again by Morita invariance ($\mathrm{End}_{\mathcal{A}}(\mathcal{H})$
is Morita-equivalent to $\mathcal{A}$ as soon as $\mathcal{H}$ is a
finitely-generated projective $\mathcal{A}$-module, which it is at
finite cutoff). The dg differential $[D, \cdot]$ contributes to the
*cyclic* homology $\mathrm{HC}_*$ before stabilization but the
periodic version washes out the contribution of an exact bounded
derivation. Formally:\ the SBI exact sequence

$$
\cdots \to \mathrm{HC}_n \xrightarrow{S} \mathrm{HC}_{n-2}
\xrightarrow{B} \mathrm{HH}_{n-1} \xrightarrow{I} \mathrm{HC}_{n-1} \to \cdots
$$

stabilizes to $\mathrm{HP}_* = \varprojlim_S \mathrm{HC}_{*+2S}$,
and the inverse limit eats any finite-rank cancellation introduced by
bounded $[D, \cdot]$. The conclusion is that **$\mathrm{HP}_*$ does
not see the Dirac operator at finite cutoff** — it sees only the
underlying Morita class of the algebra.

This is the deeper Morita-stability statement: $\mathrm{HP}_*$ as a
fiber functor is robust against the *quantitative* spectral data
(which is what carries the M2/M3 distinction via the half-integer
Camporesi–Higuchi shift). It records only the *qualitative*
K-theoretic / Morita class.

**Conclusion for the dg side:** $\mathrm{HP}_*$ of the dg-enriched
GeoVac spectral triple at every finite $n_{\max}$ is the same as
$\mathrm{HP}_*$ of the bare algebra. The Dirac operator's
contribution lives in *cyclic cohomology* $\mathrm{HP}^*$ (the dual
side), not in the homology side that Marcolli–Tabuada uses.

### Step 3: The pro-limit $\varprojlim$ with Berezin maps

Round 1 of Q5' (`sprint_q5p_qsm_litread_memo.md` §4) flagged that
the right object might not be each $\mathrm{dg}(\mathcal{T}_{n_{\max}})$
individually, but the inverse system
$\{\mathrm{dg}(\mathcal{T}_{n_{\max}})\}_{n_{\max}}$ together with
Berezin reconstruction maps $B_{n_{\max}}: C(S^3) \to
\mathcal{O}_{n_{\max}}$ (Paper 38 §L4). The hope is that the pro-limit
sees the continuum spectral triple's non-trivial $\mathrm{HP}_*$ even
though each finite-truncation does not.

**This hope does not survive structural inspection:**

(a) $\mathrm{HP}_*$ commutes with filtered colimits (Loday §5.1.6)
    but not with inverse limits in general. The pro-system of
    Morita-trivial finite-dim algebras may have non-Morita-trivial
    limit, but the limit is the standard $\mathrm{HP}_*(C^\infty(S^3))$,
    which is $\mathrm{HP}_*(C^\infty(S^3)) \cong
    H_{\mathrm{dR}}^*(\smash{S^3}) = \mathbb{Q}$ in degree $0$ and
    $\mathbb{Q}$ in degree $3$ (Connes' Hochschild–Kostant–Rosenberg
    theorem, *Noncommutative Geometry* §III.2).

(b) **This is still a Tate-pure $\mathbb{Q}$-vector space**, of total
    dimension $2$ (one in $H^0$, one in $H^3$), with the standard
    Tate twist controlling the Galois action. It is NOT a refined
    $\mathbb{Z}/3$-graded structure tracking M1/M2/M3.

(c) The Berezin reconstruction maps make the pro-system into a
    coherent system on a single carrier ($C^\infty(S^3)$ via $B_{n_{\max}}$);
    the resulting pro-$\mathrm{HP}_*$ collapses to the continuum
    $\mathrm{HP}_*(C^\infty(S^3))$, which is the de Rham cohomology of
    $S^3$. This is well-known to be pure-Tate (no MZV content, no
    cyclotomic content), matching the M1 + Hopf-base-measure
    sub-mechanism but NOT M2 or M3.

**Conclusion for the pro-limit:** the continuum-limit
$\mathrm{HP}_*$ is the de Rham cohomology of $S^3$, which is
pure-Tate. The M2/M3 mixed-Tate / cyclotomic content does NOT appear
in $\mathrm{HP}_*$ of the pro-limit either. The Marcolli–Tabuada
fiber functor sees only the trivial pure-Tate sector of the master
Mellin engine — which is the M1 slot — and does not see M2's
Seeley–DeWitt or M3's vertex-parity Dirichlet content at all.

### Step 4: But the Chern character is non-trivial!

The escape route hinted by the task prompt — "the spectral triple
data ($D$, $\chi$ for $\mathbb{Z}/2$ grading) carries more than the
algebra alone" — points to the Connes character of the spectral
triple, equivalently the Chern character of the associated Fredholm
module. This is a class

$$
\mathrm{ch}(\mathcal{T}_{n_{\max}}) \in \mathrm{HP}^*(\mathcal{A}_{\mathrm{GV}}),
$$

where $\mathrm{HP}^*$ is periodic cyclic *cohomology*, dual to
$\mathrm{HP}_*$. The Chern character pairs with K-theory:

$$
\langle \mathrm{ch}(\mathcal{T}), - \rangle: K_*(\mathcal{A}_{\mathrm{GV}}) \to \mathbb{C}.
$$

On the finite-dim commutative $\mathcal{A}_{\mathrm{GV}}^{(n_{\max})}$,
$K_0 = \mathbb{Z}^{N_{\mathrm{Fock}}}$ (one rank-1 idempotent per
vertex) and the Chern character is the spectral-action evaluation:

$$
\mathrm{ch}(\mathcal{T}_{n_{\max}})(p_v) = \mathrm{Tr}(p_v \cdot e^{-tD^2}/\mathrm{Tr}(e^{-tD^2}))
$$

for $p_v$ the rank-1 idempotent at vertex $v$. The heat-trace
information of $D$ IS encoded in this character — and this is
precisely the input the master Mellin engine
$\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$ Mellin-transforms.

**Critically, the Chern character on
$\mathrm{HP}^*$ depends on the Camporesi–Higuchi half-integer shift
$3/2$**, because the heat-trace expansion is
$$
\mathrm{Tr}(e^{-tD^2}) = t^{-3/2}(a_0 + a_2 t + \cdots),
$$
with the half-integer shift encoded in the $t^{-3/2}$ leading
divergence (M2 input at $k = 2$) AND in the parity decomposition
$D = D_{\mathrm{even}} + D_{\mathrm{odd}}$ on the half-integer
spectrum (M3 input at $k = 1$). The k-slot index $k \in \{0, 1, 2\}$
*is* visible in the Chern character data once one decomposes the
class along the operator-order grading
$\mathrm{Tr}(D^k e^{-tD^2})$.

**But the Chern character lives in
$\mathrm{HP}^*$, dual to $\mathrm{HP}_*$**. The Marcolli–Tabuada
fiber functor uses $\mathrm{HP}_*$ as input, NOT $\mathrm{HP}^*$.
The pairing
$$
\mathrm{HP}^*(\mathcal{A}) \otimes \mathrm{HP}_*(\mathcal{A}) \to \mathbb{Q}
$$
is the natural duality:\ the character is a *functional* on
$\mathrm{HP}_*$, not an *element* of it. So the k-slot information
the Chern character carries is on the *wrong side* for the
Marcolli–Tabuada Tannakian construction.

This is the structurally decisive observation. The Marcolli–Tabuada
construction asks "what Tannakian symmetries act on
$\mathrm{HP}_*(\mathrm{dg}(\mathcal{T}))$?" — and the answer is "$\mathbb{G}_m$
on Tate weights of a finite-rank pure-Tate vector space."\ It does NOT
ask "what does the Chern character pair with?" — and that's where the
k-slot lives.

### Step 5: The Connes–Marcolli $U^*$ as the alternative

The cosmic Galois group $U^*$ of Connes–Marcolli (arXiv:math/0409306)
operates on the *characterization* side, not on the homology side. In
the renormalization context, $U^*$ acts on counterterm coefficients
that are *labelled* by the equisingular flat connections (= Feynman
diagrams) that produced them. The Tannakian dual is taken with
respect to a fiber functor that *respects the labels*.

This is exactly the shape of construction Q5' would need:\ a
"motivic Galois group of the GeoVac Chern character" with a
fiber functor that respects the Mellin slot $k \in \{0, 1, 2\}$. The
construction would be:

- **Input:** the Chern character $\mathrm{ch}(\mathcal{T}_{n_{\max}}) \in
  \mathrm{HP}^*(\mathcal{A}_{\mathrm{GV}})$, decomposed along the
  Mellin slot grading
  $\mathrm{ch}_k(\mathcal{T}_{n_{\max}}) := \mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]
  \in \mathrm{HP}^*(\mathcal{A}_{\mathrm{GV}})$ at $k = 0, 1, 2$.

- **Tannakian category:** the symmetric monoidal category generated
  by $(\mathrm{ch}_0, \mathrm{ch}_1, \mathrm{ch}_2)$ inside
  $\mathrm{HP}^*$, with the master-Mellin-engine partition theorem
  (Paper 32 §VIII case-exhaustion) as the structural exactness
  statement.

- **Motivic Galois group:** $\mathrm{Aut}^\otimes$ of the
  label-respecting fiber functor.

This is the Q5' enriched-construction frontier. It does NOT fit
into the Marcolli–Tabuada machinery (which uses
$\mathrm{HP}_*$ on dg categories of perfect modules) but DOES fit
into the Connes–Marcolli $U^*$ shape (which uses
characterization-side data with labels).

**The honest finding:** Q5' Round 3 closes the Marcolli–Tabuada route
in the negative (TRIVIAL via Step 1–3) but the Connes–Marcolli $U^*$
shape route remains open (NON-TRIVIAL on $\mathrm{HP}^*$ via Step 4–5,
matching Round 1's "Strand A is the gold-standard precedent for
motivic Galois from analytic data" finding).

### Step 6: Compatibility with prior Q5' rounds

The Round 3 verdict is the consistent refinement of all four prior
Round 1 + Round 2 closures:

- **Round 1 dim sweep (NEGATIVE):** the M2/M3 $\mathbb{Q}(i)$ field-coincidence
  is generic across $d \in \{3, 5, 7\}$. Round 3 explains: the dim
  sweep negative is because the ambient MT period ring is generic at
  every odd $d$, and the M2/M3 distinction lives only at the
  source-cohomological level. The dim sweep cannot rescue source-cohomological
  content from the period ring.

- **Round 1 QSM lit-read (PARTIALLY TRANSPORTABLE):** the
  Marcolli–Tabuada machinery is the "right scaffolding" only on the
  homological side. Round 3 sharpens to: it is the WRONG scaffolding
  for what Q5' actually needs, because the k-slot lives on the
  cohomological-dual side. Connes–Marcolli $U^*$ is the right shape.

- **Round 1 Tannakian obstruction (OPEN-POSSIBLE):** there is no
  cheap categorical obstruction. Round 3 sharpens to:\ there is a
  STRUCTURAL obstruction on the Marcolli–Tabuada *machinery side*
  (not on the abstract Tannakian framework side):\ the input fiber
  functor $\mathrm{HP}_*$ does not see the k-slot. The Round 1
  OPEN-POSSIBLE is preserved because the obstruction is to a
  particular construction route, not to the abstract enrichment.

- **Round 2 k-slot Tannakian probe (BORDERLINE):** the k-slot is
  Tannakian-invisible on the standard MT fiber functor and
  Tannakian-relevant one categorical level up. Round 3 sharpens
  "one categorical level up" to "on the $\mathrm{HP}^*$ /
  Chern-character / K-homology side" — i.e., on the dual side of
  the homology used by Marcolli–Tabuada. The Round 2 BORDERLINE
  is preserved with the dual-side identification made explicit.

- **Round 2 deflation test (PARTIAL with strong DEFLATES bias):**
  Q5' deflates at the period-ring level but does not deflate on the
  k-slot index. Round 3 sharpens:\ Q5' deflates at the Marcolli–Tabuada
  source level too; the k-slot index lives on the dual cohomological
  side, which is genuinely substrate-side metadata (not period-side
  content). The deflation memo's verdict ("Q5' deflates to
  natural-geometry realization, k-slot is substrate-side") is
  affirmed and sharpened.

The five Q5' Round 1–3 memos converge on the same picture:\ the M2/M3
output rings sit cleanly in the ambient mixed-Tate category, the
k-slot index is substrate-side metadata not visible to the standard
period-side Tannakian apparatus, and the enriched-construction frontier
is the Connes–Marcolli $U^*$ shape on the K-homology / Chern-character
side, not the Marcolli–Tabuada shape on the $\mathrm{HP}_*$ side.

## What this means for sharpened Q5'

### The Marcolli–Tabuada route dies

The sharpened Q5' construction Round 2 contemplated:

> Construct $\omega^{\mathrm{tri}}: \mathrm{dg}(\mathcal{T}) \to
> \mathrm{Vec}_{\mathbb{Q}} \otimes \mathrm{IndexCat}(\{0, 1, 2\})$ via
> Marcolli–Tabuada NC motivic Galois machinery

**is structurally blocked at source**:\ $\mathrm{HP}_*$ of the bare
algebra is pure-Tate weight-0 of rank $N_{\mathrm{Fock}}$;\ $\mathrm{HP}_*$
of the dg-enriched spectral triple is also pure-Tate weight-0 of rank
$N_{\mathrm{Fock}}$ (no Dirac contribution at periodic level); the
pro-limit is the de Rham cohomology of $S^3$, also pure-Tate. The
Marcolli–Tabuada fiber functor cannot distinguish M1, M2, M3 because
$\mathrm{HP}_*$ at every level returns the same trivial pure-Tate
data. The Tannakian dual is $\mathbb{G}_m$ acting on Tate twists, with
trivial action on M2's pure-Tate output ring and no machinery for
M3's cyclotomic content.

The Round 2 sharpened Q5' enriched fiber functor in the
Marcolli–Tabuada framework therefore **dies at the source**.

### The Connes–Marcolli $U^*$ route remains alive

The K-homology / Chern-character / Connes-character side of the
spectral triple **does** carry the k-slot index. The master Mellin
engine
$\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$ is naturally a decomposition
of the Chern character along the operator-order grading; the
case-exhaustion theorem is a structural statement about which
$k$-slot produces which period. The natural construction is:

- **Functorial input:** $\mathcal{T}_{n_{\max}} \mapsto
  (\mathrm{ch}_k(\mathcal{T}_{n_{\max}}))_{k \in \{0, 1, 2\}} \in
  \mathrm{HP}^*(\mathcal{A}_{\mathrm{GV}})^{\oplus 3}$.

- **Symmetric monoidal structure:** the $k$-grading is preserved by
  tensor products (the Mellin transform respects the tensor structure
  in the relevant sense; Connes–Moscovici JLO cocycle is graded by
  operator order).

- **Tannakian dual:** the motivic Galois group acting on this
  $\mathbb{Z}/3$-graded structure.

This is the genuine multi-year construction that **Round 1
Strand A** (Connes–Marcolli $U^*$) flagged as the right shape. It is
NOT the Marcolli–Tabuada NC motivic Galois construction.

### Net deflation reading at Round 3

Q5' deflates further than Round 2 anticipated:\ the
Marcolli–Tabuada route is blocked at source. The sharpened Q5' enriched
fiber functor in the Marcolli–Tabuada framework dies. The substantive
remaining content of Q5' is the question of whether the Connes–Marcolli
$U^*$ shape construction (on the K-homology / Chern-character side)
delivers a non-trivial motivic Galois group that distinguishes
M1/M2/M3.

This is a more honest framing of the open question:\ Q5' is NOT
"motivic Galois of dg(GeoVac) via Marcolli–Tabuada"; Q5' IS "motivic
Galois of the GeoVac Chern character via Connes–Marcolli $U^*$
shape." These are categorically distinct construction frameworks.

## Honest scope

1. **No fitted parameters, no PSLQ, no implementation, no code.** Pure
   structural reasoning on published cyclic-homology facts (Loday,
   Connes, Marcolli–Tabuada, Connes–Marcolli) plus the GeoVac structural
   data from Paper 32 §III–IV.

2. **No paper edits applied.** Whether to update Paper 55 §subsec:open_m2_m3
   Q5' wording (and possibly Paper 32 §VIII open-questions list) to
   reflect the Marcolli–Tabuada-route-dies finding is a PI call.

3. **Combination rule $K = \pi(B + F - \Delta)$ "conjectural" label
   preserved.** No edits to Paper 2 framing.

4. **The Morita-triviality of $\mathcal{A}_{\mathrm{GV}}^{(n_{\max})}$
   is load-bearing.** If Paper 32's algebra definition is ever revised
   (e.g., the operator-system replacement
   $\mathcal{O}_{n_{\max}}$ of Remark 3.2 lines 282–333), the
   conclusion would need to be re-checked. The current
   Definition 3.1 is the commutative function algebra, hence
   Morita-trivial; the operator system $\mathcal{O}_{n_{\max}}$ is
   $*$-closed but NOT multiplicatively closed (Paper 32 lines 301–311),
   so it is not a dg algebra at all in the standard sense and the
   Marcolli–Tabuada machinery does not directly apply. The
   conclusion is robust against this distinction:\ either choice of
   algebra gives no M2/M3-distinguishing $\mathrm{HP}_*$ input.

5. **The non-trivial Chern character claim depends on the spectral
   data being non-trivial.** At every finite $n_{\max}$, the
   Camporesi–Higuchi Dirac has spectrum $\{\pm(n + 3/2)\}$ with
   degeneracies; the heat trace is non-trivial and the Chern
   character is non-trivial. The $\mathrm{HP}^*$-side k-slot
   information is real, not vacuous.

6. **Sprint-pass depth caveat.** This probe was 1-day sprint-scale.
   A deeper examination might surface a non-standard $\mathrm{HP}_*$
   refinement (e.g., the Hood–Jones bivariant cyclic homology, or
   the operator-system cyclic homology constructed by van~Suijlekom
   and collaborators in the truncation context) that delivers
   $\mathbb{Z}/3$-graded source data. This was not checked here. The
   verdict is robust against the *standard* Marcolli–Tabuada
   construction; non-standard refinements remain a multi-year
   research target that this memo flags but does not pursue.

7. **The verdict TRIVIAL is on the Marcolli–Tabuada source side.**
   The framework is not claiming any global "GeoVac has no
   non-trivial motivic content" — only that the Marcolli–Tabuada
   route to the M2/M3 distinction is structurally blocked. The
   Connes–Marcolli $U^*$ route on the K-homology side may have
   non-trivial content; that is an open multi-year question, not
   closed by this probe.

8. **The dg algebra of operator-system multipliers is structurally
   different from the dg algebra of perfect modules.** The
   Connes–van~Suijlekom operator system $\mathcal{O}_{n_{\max}}$ is a
   $*$-closed subspace of $M_N(\mathbb{C})$ that is not
   multiplicatively closed; its dg structure (if any) would be
   different from the standard
   $\mathrm{End}^{\mathrm{dg}}_{\mathcal{A}}(\mathcal{H})$ used by
   Marcolli–Tabuada. The propagation number computation in Paper 32
   Proposition 3.3 (prop = 2 verbatim from Connes–vS Toeplitz $S^1$
   case) gives evidence that $\mathcal{O}_{n_{\max}}$ is a structured
   object with non-trivial operator-system content. Whether this
   operator-system content can be Tannakian-promoted via
   van~Suijlekom-style construction is a separate research target,
   not addressed here.

## One-line verdict

**TRIVIAL: $\mathrm{HP}_*(\mathcal{A}_{\mathrm{GV}}^{(n_{\max})}) =
\mathrm{HP}_*(\mathbb{C})^{\oplus N_{\mathrm{Fock}}}$ pure-Tate
weight-0 at every $n_{\max}$, dg-enrichment via bounded $[D, \cdot]$
does not escape Morita-stability, the pro-limit $\varprojlim
\mathrm{HP}_*(\mathrm{dg}(\mathcal{T}))$ equals $\mathrm{HP}_*(C^\infty(S^3))$
= pure-Tate de Rham cohomology of $S^3$, so the Marcolli–Tabuada NC
motivic Galois route to the sharpened Q5' enriched fiber functor
$\omega^{\mathrm{tri}}$ is blocked at source; the k-slot index lives
on the dual cohomological side $\mathrm{HP}^*$ via the Chern
character of the spectral triple, which is the natural input to the
Connes–Marcolli cosmic Galois group $U^*$ shape, not to Marcolli–Tabuada.**

Sharpened Q5' in the Marcolli–Tabuada framework **dies**;\ Q5'
deflates further toward "GeoVac IS the natural-geometry realization
of ambient $M^{\mathrm{GV}}$, with k-slot being substrate-side metadata
on the dual cohomological side."\ Connes–Marcolli $U^*$ remains a
multi-year open alternative.

## Files used

### Memos read
- `debug/sprint_q5p_deflation_test_memo.md` (Round 2 deflation test:\
  PARTIAL with DEFLATES bias)
- `debug/sprint_q5p_k_slot_tannakian_memo.md` (Round 2 k-slot
  Tannakian probe:\ BORDERLINE)
- `debug/sprint_q5p_qsm_litread_memo.md` (Round 1 QSM lit-read:\
  PARTIALLY TRANSPORTABLE)
- `debug/sprint_q5p_tannakian_obstruction_memo.md` (Round 1 cheap
  obstruction probe:\ OPEN-POSSIBLE)
- `debug/sprint_q5p_dim_sweep_memo.md` (Round 1 dim sweep:\ NEGATIVE
  on entry point (i))
- `debug/sprint_q5p_greenfield_marcolli_transport_memo.md` (Round 1
  Greenfield–Marcolli transport probe:\ NOT TRANSPORTABLE)

### Papers read
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
  §III (definition of $\mathcal{A}_{\mathrm{GV}}^{(n_{\max})}$, lines
  259–333), §IV (axiom audit, lines 851–1000), §V (representation,
  lines 717–725), §VI (bounded commutators, lines 803–814), real
  structure Prop 5.4 (lines 919–938), finite-$n_{\max}$ axiom audit
  (lines 941–1000).
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §6.1
  (joint Mellin engine, role-disjointness, eq:ring_nesting), §Q5'
  open subsection (lines 1474–1499).

### Published references (cited from Paper 32 + Paper 55 bibliographies)
- Loday, J.-L. *Cyclic Homology.* 2nd edition, Grundlehren der
  mathematischen Wissenschaften 301, Springer (1998). §1.4 Thm 1.2.4
  (Morita invariance), §5.1.6 (filtered colimits), §2.5
  ($\mathbb{Z}/2$ grading of $\mathrm{HP}_*$). Standard reference.
- Connes, A. *Noncommutative Geometry.* Academic Press (1994). §III.2
  (Hochschild–Kostant–Rosenberg theorem), §IV (Chern character of
  spectral triples). Standard reference.
- Marcolli, M.; Tabuada, G. "Noncommutative numerical motives,
  Tannakian structures, and motivic Galois groups." J. Eur. Math. Soc.
  18 (2016), 623–655. arXiv:1110.2438. Cor 4.2 ($\mathrm{HP}_*$ as
  $\mathbb{Z}/2$-graded symmetric monoidal functor); Thm 1.1
  (NC motivic Galois construction).
- Connes, A.; Marcolli, M. "Renormalization, the Riemann–Hilbert
  correspondence, and motivic Galois theory." arXiv:math/0409306
  (2004). The $U^*$ cosmic Galois group acting on counterterm
  coefficients labelled by Feynman diagrams — the structural shape
  Q5' Round 3 identifies as the alive alternative.
- Connes, A.; Moscovici, H. "Type III and spectral triples." in
  *Traces in number theory, geometry and quantum fields*, Aspects
  Math. E38, Vieweg (2008), 57–71. Background on twisted spectral
  triples / Chern character on type III triples (not directly used
  in this memo but referenced by Round 1's Strand B which Round 1.5
  ruled out).
- Marcolli, M.; van~Suijlekom, W. D. "Gauge networks in
  noncommutative geometry." J. Geom. Phys. 75 (2014), 71–91.
  arXiv:1301.3480. Gauge networks construction whose algebra
  $C(V)$ is the ambient category for $\mathcal{A}_{\mathrm{GV}}$.
- Connes, A.; van~Suijlekom, W. D. "Spectral truncations in
  noncommutative geometry and operator systems." Comm. Math. Phys.
  383 (2021), 2021–2067. Definition 2.39 (propagation number),
  Proposition 4.2 ($\mathrm{prop}(C(S^1)^{(n)}) = 2$ verbatim, matched
  by Paper 32 Prop 3.3 for $\mathcal{O}_{n_{\max}}$ at $n_{\max} \in
  \{2, 3, 4\}$).

### Scripts (none required)
Pure structural reasoning. No code, no PSLQ, no numerical checks.

### Round 3 takeaway in one sentence

The Marcolli–Tabuada noncommutative motivic Galois construction
cannot see GeoVac's k-slot index because $\mathrm{HP}_*$ of the bare
finite-truncation algebra is pure-Tate weight-0 of rank
$N_{\mathrm{Fock}}$, the dg-enrichment via bounded $[D, \cdot]$ does
not escape Morita-stability, and the pro-limit converges to the
pure-Tate de Rham cohomology of $S^3$;\ the k-slot lives on the dual
cohomological side ($\mathrm{HP}^*$ / Chern-character /
K-homology), which is the input for the Connes–Marcolli cosmic
Galois group $U^*$ shape, not for Marcolli–Tabuada.
