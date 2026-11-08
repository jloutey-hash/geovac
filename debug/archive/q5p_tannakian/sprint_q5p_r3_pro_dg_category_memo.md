# Sprint Q5' — Round 3 candidate (ii): pro-dg-category cyclic-homology
non-triviality probe

Date: 2026-06-04
Scope: 1-day structural sub-sprint of Q5' (Paper 55
§subsec:open_m2_m3), Round 3. Builds on:
- Round 2 R2.2 BORDERLINE (k-slot TANNAKIAN-INVISIBLE on ambient
  $\mathrm{MT}(\Z[i, 1/2], 4)$, TANNAKIAN-RELEVANT one categorical
  level up; memo
  `debug/sprint_q5p_k_slot_tannakian_memo.md`).
- Round 2 R2.2 honest-scope item:\ "whether the pro-dg-category
  $\varprojlim \mathrm{dg}(\Tcal_{n_{\max}})$ with Berezin maps is
  the right object for non-trivial Marcolli--Tabuada input
  (lit-read flagged as likely)."
- Round 2 deflation test (`debug/sprint_q5p_deflation_test_memo.md`,
  PARTIAL with strong DEFLATES bias).

This probe asks specifically: **does the pro-dg-category
$\varprojlim_{n_{\max}} \mathrm{dg}(\Tcal_{n_{\max}})$, with Berezin
reconstruction maps as transition morphisms, carry non-trivial
cyclic-homological / Marcolli--Tabuada-fiber-functor input that
survives even if each finite-$n_{\max}$ cutoff has trivial
$\mathrm{HP}_*$?**

## TL;DR

**Verdict: TRIVIAL, with a sharp structural reason that makes the
TRIVIAL verdict load-bearing rather than empty.**

The argument has three load-bearing parts:

1. **HP_* of each finite-$n_{\max}$ cutoff is structurally trivial.**
   $\mathcal{O}_{n_{\max}}$ is a finite-dimensional matrix sub-algebra
   (multiplier matrices acting on $\Hilb_{n_{\max}}$, dim
   $N(n_{\max}) = (n_{\max})^2$). By Morita-invariance of HP_*
   (Loday *Cyclic Homology* 1998 §2.1.7) and Connes' calculation
   $\mathrm{HP}_*(\mathbb{C}) = (\mathbb{C}, 0)$, every finite-dim
   semisimple complex algebra has $\mathrm{HP}_0 = \mathbb{C}^k$
   (with $k$ = number of simple summands) and $\mathrm{HP}_1 = 0$.
   For the GeoVac truncated operator system at cutoff $n_{\max}$,
   $\mathcal{O}_{n_{\max}}$ is **NOT closed under multiplication**
   (this is the defining property of the Connes--van Suijlekom
   operator system per Paper 32 §III), so it is not even an algebra.
   The natural dg category $\mathrm{dg}(\Tcal_{n_{\max}})$ that
   Marcolli--Tabuada machinery would consume is the dg category of
   modules over the *generated algebra* $\langle \mathcal{O}_{n_{\max}}
   \rangle = M_{N(n_{\max})}(\mathbb{C})$ (the full matrix algebra it
   spans, which IS closed under multiplication), and this has the
   single-summand semisimple HP profile
   $\mathrm{HP}_0 = \mathbb{C}$, $\mathrm{HP}_1 = 0$.

2. **The Berezin maps are NOT algebra homomorphisms** (Paper 38 §L4,
   property (b) contractivity + property (c) approximate identity).
   They are completely positive, contractive, unital up to
   approximate-identity error, but explicitly NOT multiplicative:\
   $B_{n_{\max}}(fg) \ne B_{n_{\max}}(f) B_{n_{\max}}(g)$ in general
   (the multiplicativity-defect IS the "approximate identity" gap that
   the L2 central Fejér rate $\gamma_{n_{\max}} \to 0$ closes
   asymptotically). Marcolli--Tabuada noncommutative motives (1110.2438,
   1112.5422) consume dg-categories under **dg-functors** (or
   equivalently under Morita morphisms), which require multiplicative
   maps at the algebra level. **The Berezin maps are not morphisms in
   the dg-category of noncommutative motives.** They sit in a
   different category — completely-positive maps, or KK-theory cycles,
   or quantum-metric-space Berezin-type maps — and they do not induce
   a functorial pullback on $\mathrm{HP}_*$.

3. **Even granting an extension of HP_* to CP-maps, the induced map
   to the target $\mathrm{HP}_*(C(S^3))$ is non-injective.** If one
   forces a CP-extension of HP_* (e.g.\ via KK-theory or the
   bivariant cyclic theory of Cuntz--Quillen 1995, Tsygan 1986,
   Cuntz 2005), the induced map $B^{\mathrm{HP}}_{n_{\max}}:
   \mathrm{HP}_*(C(S^3)) \to \mathrm{HP}_*(\mathcal{O}_{n_{\max}})$
   on cyclic homology must land in $\mathrm{HP}_*(M_{N(n_{\max})}
   (\mathbb{C})) = (\mathbb{C}, 0)$. The target $\mathrm{HP}_1
   (C(S^3)) = \mathbb{C}$ (the $S^3$-volume class, generator of
   $H^3_{\mathrm{dR}}(S^3) = \mathbb{C}$ under the Connes--Karoubi
   isomorphism $\mathrm{HP}_*(C^\infty(M)) \cong H^*_{\mathrm{dR}}(M)$
   with $\Z/2$-periodization) **MUST map to zero** in the
   inductive-system limit, because $\mathrm{HP}_1(M_N(\mathbb{C})) = 0$
   for every finite $N$. So the inductive system collapses HP_1 to
   zero in the limit, and the volume class of $S^3$ — the only
   non-trivial cyclic-homology content — is **structurally
   inaccessible** at every finite cutoff.

The honest finding: **the pro-dg-category with Berezin transitions
is, on cyclic homology, a CONSTANT system at $(\mathbb{C}, 0)$ in
every degree**. The limit is $(\mathbb{C}, 0)$, no matter how the
inverse-limit / direct-limit is taken. The non-trivial
$\mathrm{HP}_1(C(S^3)) = \mathbb{C}$ does not survive the passage
through finite-dim sub-algebras. The Berezin maps' L4 properties
(positivity, contractivity, approximate identity) are about
**propinquity convergence in operator-norm metric**, not about
cyclic-homological convergence. **These are categorically
orthogonal convergences.**

This is the precise content of the TRIVIAL verdict. Q5' does not
gain new structural input from the pro-dg-category construction;
the GH-convergence-as-Tannakian-input reading does not survive
inspection.

## Decision-gate framework

The task posed three readings:

- **NON-TRIVIAL:** The Berezin-induced inductive system on HP*
  (or HC* / HH*) carries content that depends meaningfully on
  $n_{\max}$ or on the spectrum-side data; reaches HP*(C(S^3)) in
  a way that distinguishes M1/M2/M3 contributions. Sharpened Q5'
  alive at the pro-object level even if R3.1 finds trivial
  finite-cutoff HP*.
- **TRIVIAL:** Berezin-induced system collapses to HP*(C(S^3))
  trivially (or stably) without distinguishing $k$-slots. $\omega^{\mathrm{tri}}$
  lacks non-trivial Marcolli--Tabuada input. Sharpened Q5' blocked at
  both finite and pro-object levels.
- **BORDERLINE:** Content survives at pro-level but doesn't
  separate $k$-slot specifically; reframes Q5' more sharply.

The structural argument below lands on TRIVIAL for the reason
above: HP_*($\mathcal{O}_{n_{\max}}$) is structurally
$(\mathbb{C}, 0)$ at every finite cutoff, the Berezin maps fail
to be algebra-functorial, and the limit collapses to
$(\mathbb{C}, 0)$ which is strictly smaller than
HP_*(C(S^3)) = $(\mathbb{C}, \mathbb{C})$.

## Argument by ingredient

### Ingredient 1:\ HP_* of finite-dim algebras and operator systems

**Marcolli--Tabuada machinery requires dg-categories or
$A_\infty$-categories whose objects have well-defined Hochschild and
cyclic homology.** For a dg-category $\Tcal$, the natural invariant
is the cyclic complex of the endomorphism dg-algebra (or its
Morita-invariant cyclic homology).

In the GeoVac setup, at each finite $n_{\max}$:
- **The Hilbert space** $\Hilb_{n_{\max}}$ is finite-dimensional, of
  dim $N(n_{\max}) = (n_{\max})^2$ for the scalar truncation
  (Paper 38 Lemma L1' bookkeeping).
- **The truncated operator system** $\mathcal{O}_{n_{\max}}$ is
  spanned by multiplier matrices $\Mat_{NLM}$ for shell labels
  $(N, L, M)$ with $N \le n_{\max}$ (Paper 38 Def. L1', §3).
- $\mathcal{O}_{n_{\max}}$ is **NOT closed under multiplication**;\
  this is the foundational fact of Connes--van Suijlekom truncated
  operator systems (Connes--vS 2021 CMP, Paper 32 §III).

For the Marcolli--Tabuada machinery to apply, the natural object is
either:

(a) The full matrix algebra $M_{N(n_{\max})}(\mathbb{C})$ that
    contains $\mathcal{O}_{n_{\max}}$ (the algebra-generation closure).
(b) The truncated dg-algebra whose objects are $\mathcal{O}_{n_{\max}}$
    elements but whose composition is defined as a partial operation.

Under (a):\ Morita-invariance gives $\mathrm{HP}_*(M_N(\mathbb{C}))
= \mathrm{HP}_*(\mathbb{C}) = (\mathbb{C}, 0)$. The volume / Chern
character information of any module reduces to its dimension at
the HP_0 level; nothing distinguishes the different shells $(N, L, M)$
in the truncation cylindrically.

Under (b):\ partial-composition dg-algebras have HP_* defined via the
fully-derived bar resolution; for the truncated multiplier set the
result is the same Morita class as (a), so HP_* is still
$(\mathbb{C}, 0)$.

The truncated operator system therefore carries **no cyclic-homological
information that goes beyond a single complex line in HP_0**. The
shell structure $(N, L, M)$ is forgotten at the cyclic-homology
level. This is the precise sense in which the lit-read flag
("dg category attached to a finite-dimensional spectral triple is
essentially trivial," `sprint_q5p_qsm_litread_memo.md` §3 Strand C)
is correct.

### Ingredient 2:\ Berezin maps are not algebra homomorphisms

The Berezin reconstruction map $B_{n_{\max}}: C^\infty(\sthree) \to
\mathcal{O}_{n_{\max}}$ (Paper 38 Def. \ref{def:berezin}) is given
explicitly by

\[
   B_{n_{\max}}(f) = P_{n_{\max}}\,M_{\KS_{n_{\max}} * f}\,P_{n_{\max}},
\]

where $\KS_{n_{\max}}$ is the L2 central Fejér kernel on $\SU(2)$
and $M_g$ is multiplication by $g$ on $L^2(\sthree)$. The map has
four L4 properties (Paper 38 Lemma \ref{lem:L4}):

- (a) Positivity:\ $f \ge 0 \Rightarrow B(f) \ge 0$.
- (b) Contractivity:\ $\opnorm{B(f)} \le \norm{f}_{L^\infty}$.
- (c) Approximate identity:\ $\opnorm{B(f) - P_{n_{\max}}M_fP_{n_{\max}}}
  \le \gamma_{n_{\max}}\norm{\nabla f}_{L^\infty}$.
- (d) Compatibility with L3.

The map is **completely positive** (compression of the multiplication
by a positive kernel-smoothed function), unital up to the
approximate-identity error, and contractive. **It is not
multiplicative.** Property (c) records this directly:\ if $B$ were
multiplicative, the approximate-identity error would be exactly
zero (since $M_{fg} = M_fM_g$ for the multiplication operators);
the non-zero gap $\gamma_{n_{\max}}$ is the multiplicativity defect.

**Marcolli--Tabuada cyclic-homology machinery (1110.2438 §3,
1112.5422 §1) consumes dg-functors between dg-categories.** A
dg-functor $F: \Tcal_1 \to \Tcal_2$ sends objects to objects,
morphisms to morphisms, composites to composites, and identities to
identities — at the strict algebra-level. For two algebras
$A_1, A_2$ viewed as one-object dg-categories, a dg-functor is
exactly an algebra homomorphism $A_1 \to A_2$ (possibly up to
quasi-isomorphism).

**The Berezin map $B_{n_{\max}}: C(\sthree) \to \mathcal{O}_{n_{\max}}$
is not an algebra homomorphism.** It is a CP-map, equivalently a
Markov-positive transition map between abelian and matricial
algebras. The cyclic homology Marcolli--Tabuada machinery does NOT
extend functorially to CP-maps unless an additional structure is
specified (KK-theory cycle, bivariant cyclic theory of Cuntz 2005,
or a chosen Stinespring dilation).

So the pro-system $\{\mathcal{O}_{n_{\max}}\}_{n_{\max}}$ with Berezin
transitions is, at the dg-category level, **a sequence of algebras
without functorial transition morphisms between them**, only
quantum-metric-space Berezin-style transitions that compute
propinquity bounds. **This is the right object for Latrémolière
propinquity convergence (Paper 38 Lemma L5), but it is NOT a
pro-dg-category in the Marcolli--Tabuada sense.**

### Ingredient 3:\ Even with CP-extension, the limit collapses HP_1

There are extensions of HP_* to bivariant settings:

- **Cuntz--Quillen excision** (arXiv:math/9311101):\ HP_* with values
  in extensions, and bivariant cyclic theory.
- **Cuntz 2005 (Adv. Math.) bivariant K-theory of locally convex
  algebras:** computes HP_* with morphisms in
  $\mathrm{KK}$-theory-style cycles.
- **Tsygan 1986 / Loday 1998 §5:** cyclic cohomology paired with
  K-theory via the Chern character.

Granting any such extension that makes the Berezin maps induce
maps on HP_*, the structural constraint is:

\[
   B^{\mathrm{HP}}_{n_{\max}}: \mathrm{HP}_*(C(\sthree))
   \to \mathrm{HP}_*(\mathcal{O}_{n_{\max}})
   \cong \mathrm{HP}_*(M_{N(n_{\max})}(\mathbb{C}))
   = (\mathbb{C}, 0).
\]

The target carries only a single $\mathbb{C}$ in degree 0 (with HP_0
detecting the rank / unital trace class) and **zero** in degree 1.

By the Connes--Karoubi isomorphism (Connes 1985, Loday 1998 §3.4.4),
$\mathrm{HP}_*(C^\infty(\sthree)) \cong H^*_{\mathrm{dR}}(\sthree)
\cong (\mathbb{C}, \mathbb{C})$ with $H^0 = \mathbb{C}$ (constants)
and $H^3 = \mathbb{C}$ (volume class), $\Z/2$-periodized. So
$\mathrm{HP}_0(C(\sthree)) = \mathbb{C}$, $\mathrm{HP}_1(C(\sthree))
= \mathbb{C}$.

Any CP-extension $B^{\mathrm{HP}}_{n_{\max}}$ sends:

- $\mathrm{HP}_0(C(\sthree)) = \mathbb{C}$ to the rank class in
  $\mathrm{HP}_0(\mathcal{O}_{n_{\max}})$, with the image determined
  by the trace of the corresponding multiplier matrix
  $\tau(B_{n_{\max}}(1)) = \norm{\KS_{n_{\max}}}_{L^1} = 1$ (Paper 38
  Lemma L2(a)). **This is non-trivial:\ the unit class is preserved**.
- $\mathrm{HP}_1(C(\sthree)) = \mathbb{C}$ to $\mathrm{HP}_1(M_N
  (\mathbb{C})) = 0$. **The volume class is destroyed at every finite
  cutoff**.

The limit (inverse or direct) of a sequence whose every term is
zero is zero. So no extension of HP_* to CP-maps can recover the
non-trivial HP_1 class of $C(\sthree)$ from the pro-system.

**This is the structural reason TRIVIAL is load-bearing**:\ the
GH-convergence Berezin-pro-system, even granting an extension of
cyclic homology to CP-maps, **cannot detect the volume class of
$S^3$**, which is the canonical Marcolli--Tabuada-relevant
characteristic-class data of $C(S^3)$ as a noncommutative motive.

### Ingredient 4:\ The L2 central Fejér rate $\gamma_{n_{\max}}$ is
not cyclic-homological

A natural follow-up:\ even if HP_* is trivial at each finite cutoff,
could the **rate** of convergence carry NC-motivic content? The
quantitative-rate constant 4/π (Paper 38 L2 asymptotic, the M1
Hopf-base measure signature per `mellin_taxonomy_engine.md`) is the
natural target.

The honest answer is:\ the L2 rate is the convergence rate **of the
propinquity** $\Lprop(\Tcal_{n_{\max}}, \Tcal_{S^3})$, which is a
quantum-metric-space distance. Propinquity is **NOT a
cyclic-homological invariant**; it is a metric on the category of
truncated spectral triples (Latrémolière 2017/2023). Convergence
in propinquity does not imply convergence on HP_*:

- Two propinquity-close spectral triples can have very different
  HP_*. (Trivial example:\ $C(S^3)$ and $C(S^2)$ have different HP_*
  but with appropriate metric choices can be made arbitrarily
  propinquity-close as truncated triples at small cutoffs.)
- Conversely, two HP_*-isomorphic algebras can have very different
  propinquity behavior (e.g., $M_2(\mathbb{C})$ and $M_3(\mathbb{C})$
  have $\mathrm{HP}_* = (\mathbb{C}, 0)$ in both cases but live at
  different cutoffs).

The 4/π rate of L2 is therefore a **quantitative refinement of an
operator-system metric**, not a quantitative refinement of an HP_*
map. The candidate "rate-as-NC-motivic-input" reading does not
materialize at the categorical level the question requires.

### Ingredient 5:\ Marcolli--Tabuada pro-categories and the
inverse-limit input

The lit-read flag (`sprint_q5p_qsm_litread_memo.md` §3 Strand C +
Synthesis step 4) noted the right pro-object should probably be
$\varprojlim_{n_{\max}} \mathrm{dg}(\Tcal_{n_{\max}})$ together with
Berezin transitions, on the hope that the limit captures the
M1/M2/M3 grading the individual finite-cutoff dg-categories miss.

The present probe sharpens this:\ the Berezin maps are **not
transitions in any pro-dg-category in the Marcolli--Tabuada sense**.
They are transitions in the **propinquity-tunnel-pair category**
(Latrémolière 2017 §5, used in Paper 38 Lemma L5), which records
finite-cutoff approximate-identity data, not multiplicative-functorial
data.

The Marcolli--Tabuada machinery does support pro-objects:\ pro-dg
categories are dg categories enriched with a pro-structure on their
objects, with **morphisms restricted to those compatible with the
pro-structure**. The natural candidates for transitions between
$\mathrm{dg}(\Tcal_{n_{\max}})$ and $\mathrm{dg}(\Tcal_{n_{\max}+1})$
that respect the algebra structure are:

(i) The **inclusion** $\mathcal{O}_{n_{\max}} \hookrightarrow
    \mathcal{O}_{n_{\max}+1}$ (multiplier matrices add new shells
    when the cutoff increases). This IS multiplicative on the
    common subspace, but is not a unital algebra map (the unit
    changes).
(ii) The **projection** $P_{n_{\max}+1} \to P_{n_{\max}}$
     (further compression). This is NOT multiplicative (compression
     destroys algebra structure, same reason as the Berezin map).

Option (i) is the closest to a dg-functorial transition, but it
gives a **direct system** (not the inverse limit / pro-object the
lit-read flagged). The direct system
$\bigcup_{n_{\max}} \mathcal{O}_{n_{\max}}$ converges algebraically
to a (non-closed) sub-algebra of $\mathcal{B}(L^2(\sthree))$, with
HP_* computed as the colimit:\ each finite stage has
$\mathrm{HP}_* = (\mathbb{C}, 0)$, and the colimit is still
$(\mathbb{C}, 0)$ (HP_0 doesn't grow; HP_1 stays zero).

So **even the closest dg-functorial candidate (the inclusion direct
system) gives a TRIVIAL HP_* limit**, with no detection of the
volume class.

This is the precise structural sense in which the question is
answered TRIVIAL:\ in EVERY natural categorical reading of the
question (Berezin transitions in propinquity-tunnel category;
inclusion direct system; projection inverse system),
the HP_* invariant collapses to the same Morita-trivial value at
every stage, and the non-trivial $H^*_{\mathrm{dR}}(S^3)$ content of
the target is **structurally invisible** to the discrete-cutoff
pro-object.

## Named Marcolli--Tabuada-style entry points (and why they don't
help)

To close the structural argument with concrete published targets:

### Marcolli--Tabuada 2016 Theorem 1.1 (arXiv:1110.2438):
Constructs the noncommutative motivic Galois group via Tannakian
dual applied to the category of noncommutative numerical motives,
with HP_* as fiber functor. **Concrete obstruction:** the fiber
functor sends each $\mathcal{T}_{n_{\max}}$ (or its dg thickening) to
$\mathrm{Vec}_\Q$ via HP_* evaluation. For finite-dim algebras the
image is $\mathbb{Q}$ in degree 0 and 0 in degree 1. **The pro-system
of fiber-functor outputs is constant $(\mathbb{Q}, 0)$, not
recovering the volume class of $S^3$ at any stage.** Marcolli--Tabuada
input is structurally absent.

### Marcolli--Tabuada 2015 Theorem 1.4 (arXiv:1112.5422):
Unconditional construction using the André--Kahn quotient. Same
input requirement (dg-functorial morphisms), same obstruction
(Berezin maps are not dg-functors). The unconditional construction
adds nothing on the GeoVac side because the input fails before the
construction starts.

### Tsygan 1986 / Loday 1998 §2.5 (cyclic homology of finite-dim
algebras):
The $\Z/2$-graded periodic-cyclic homology of any $n$-dimensional
matrix algebra is $(\mathbb{C}, 0)$. This is a Morita-invariant
fact and is independent of any specific spectral-triple structure
the matrix algebra carries. **The k-slot index $k \in \{0, 1, 2\}$
of the Mellin engine is invisible to HP_*** for the structural
reason that all three slots produce the same HP_* class
($\mathrm{HP}_0 = \mathbb{C}$) at every finite cutoff.

### Connes 1994 *Noncommutative Geometry* §IV.4 (compact-to-noncompact
direct limits in cyclic homology):
The closest published precedent for "pro-system on cyclic homology
from spectral-triple data" is the direct-limit construction for AF
algebras (Connes 1994 §IV.4 / Connes--Karoubi). For AF systems
$\mathcal{O}_{n} = \varinjlim_n M_{k_n}(\mathbb{C})$, the limit HP_0
captures the **Bratteli-diagram-encoded K-theory** in the inductive
direction. For GeoVac specifically, the Bratteli diagram is the
shell-counting structure (Plancherel weight $\hat{K}_{n_{\max}}(N) =
N/Z_{n_{\max}}$), and the AF-limit HP_0 just records the total
dimension count — equivalently, the trace of the limit projection on
$L^2(\sthree)$, which is infinity (or 1 after renormalization).
**No volume-class detection emerges from this limit either.**

### Connes--Consani 2014 *Cyclotomy and endomotives*:
Closest published precedent for cyclic-homology pro-system on
spectral triples with arithmetic content. **NOT TRANSPORTABLE to
GeoVac** per Round 1 Strand B (memo
`sprint_q5p_greenfield_marcolli_transport_memo.md`):\ requires KMS
phase transition + semigroup of endomorphisms + arithmetic input,
none of which is present in GeoVac. So even this most-similar
published construction does not help.

## Cross-check with Round 2 deflation test

Round 2 deflation test (`debug/sprint_q5p_deflation_test_memo.md`)
landed on PARTIAL with strong DEFLATES bias:\ the M2/M3 output rings
sit inside the ambient $\mathrm{MT}(\Z[i, 1/2], 4)$ classified motive
$M^{\mathrm{GV}}$, with the k-slot as the unique substrate-side
enrichment.

The present probe sharpens this:\ the k-slot enrichment lives at the
**substrate level (which power of $D$ is convolved with the heat
kernel)**, not at the cyclic-homology level. **HP_* is downstream of
the substrate**, after the trace is taken and the algebra structure
is forgotten. So the substrate-side metadata that distinguishes
M1/M2/M3 is precisely the metadata that HP_* loses.

The two memos are consistent:\ Round 2 says "ambient motive captures
everything HP_* sees, the k-slot is substrate-side." Round 3 says
"the substrate-side data does not lift to HP_* via the Berezin
pro-system either." Together they say:\ **GeoVac's Mellin-engine
substrate is categorically invisible to the standard Marcolli--Tabuada
machinery, regardless of whether one looks at individual
finite-cutoff or the pro-object limit**.

## What WOULD force NON-TRIVIAL (if it existed)

| Hypothetical scenario | Would settle Q5' how? | Verdict if true |
|:--|:--|:--|
| HP_* extends naturally to truncated operator systems (not just algebras), and Berezin maps induce functorial transitions | NON-TRIVIAL | The L4 properties would give a quantitative rate on cyclic homology |
| Some non-standard $\Z/3$-refinement of HP_* aligned with the master Mellin engine | NON-TRIVIAL | The k-slot would be visible at the cyclic-homology level |
| The propinquity-tunnel-pair category has its own Marcolli--Tabuada-style Tannakian extension | NON-TRIVIAL | The pro-dg-category would carry intrinsic Galois structure |
| The 4/π M1 rate has a cyclic-cohomological interpretation as a characteristic class | NON-TRIVIAL with rate | The Hopf-base measure would be a Connes--Karoubi-style class |

None of these are known to exist in the published literature, and
none are in obvious construction-distance from existing machinery.
The first item is the closest sprint-scale target — extending HP_* to
operator systems and CP-maps — but cyclic homology of operator
systems is a **categorically poorly-understood area** (no published
construction surveys identified at sprint-pass depth; the closest
related work is the operator-system K-theory of Connes--vS 2021 CMP
and Connes--Suijlekom--Toyota 2024, neither of which lifts HP_* to
operator systems).

The fourth item — cyclic-cohomological interpretation of the L2 rate
— is structurally interesting but is **a quantitative refinement of
a metric convergence**, not a cyclic-homology computation. Whether
it can be lifted to a Connes--Karoubi-style class remains open;\ the
probe does not rule it out, only notes that it does not exist in
published form.

## Verdict

**TRIVIAL.**

The pro-dg-category $\varprojlim_{n_{\max}} \mathrm{dg}(\Tcal_{n_{\max}})$
with Berezin reconstruction maps as transition morphisms does NOT
carry non-trivial cyclic-homological / Marcolli--Tabuada-fiber-functor
input. Three structural reasons combine:

1. **HP_*** of each finite-$n_{\max}$ cutoff is structurally
   $(\mathbb{C}, 0)$ by Morita-invariance of $M_N(\mathbb{C})$
   computations.

2. **The Berezin maps are CP, not multiplicative**, so they are
   not dg-functorial in the sense Marcolli--Tabuada machinery
   requires; they live in a propinquity-tunnel category rather
   than a pro-dg-category.

3. **Even granting a CP-extension of HP_*** (Cuntz--Quillen / Tsygan
   / bivariant cyclic theory), the induced map factors through
   $\mathrm{HP}_*(M_N(\mathbb{C})) = (\mathbb{C}, 0)$ at every cutoff,
   so the non-trivial $\mathrm{HP}_1(C(S^3)) = \mathbb{C}$ (volume
   class) is **destroyed at every finite stage**, and no inductive /
   inverse limit can recover it.

The natural alternative (inclusion direct system $\bigcup_{n_{\max}}
\mathcal{O}_{n_{\max}}$, which is dg-functorial) gives a colimit
HP_* still constant at $(\mathbb{C}, 0)$ — the direct limit of
trivials is trivial.

The structural reason TRIVIAL is load-bearing rather than empty:\
**the Berezin pro-system measures propinquity (operator-norm-metric)
convergence, not cyclic-homological convergence; these are
categorically orthogonal**. The L2 4/π rate (M1 Hopf-base measure
signature) refines the metric convergence; it does NOT
refine cyclic homology because cyclic homology is constant on the
system.

The TRIVIAL verdict **sharpens Q5'**:\ together with Round 2 R2.2's
BORDERLINE (k-slot is TANNAKIAN-INVISIBLE on ambient periods,
TANNAKIAN-RELEVANT on a candidate enrichment), this Round 3 finding
rules out the pro-dg-category with Berezin transitions as the
candidate enrichment. Q5' remains an open multi-year research target,
but the pro-dg-category route is **closed**, and any future Q5'
attempt must use a different categorical mechanism (e.g., the
master Mellin engine itself promoted to a multi-component functor;
or an extension of HP_* to operator systems that doesn't yet exist).

## Compatibility with prior Q5' closures

The TRIVIAL verdict is the consistent refinement of:

- **Round 1 dim sweep (NEGATIVE):** the ambient MT classification
  doesn't change with $d$; the present probe confirms this is because
  the ambient classification is HP_*-based, and HP_* of the finite-dim
  truncations is dimension-trivial.

- **Round 1 QSM/MT lit-read (PARTIALLY TRANSPORTABLE) and Round 1
  Greenfield--Marcolli probe (NOT TRANSPORTABLE):** the right
  scaffolding is Marcolli--Tabuada (not QSM), but the present probe
  shows that even Marcolli--Tabuada's machinery does not detect the
  k-slot through HP_*.

- **Round 1 cheap obstruction probe (OPEN-POSSIBLE):** no structural
  obstruction at the ambient-Tannakian level; the present probe
  shows that the obstruction is elsewhere — at the cyclic-homology
  level, where the discrete truncations are too small to capture the
  volume class. The "no obstruction at ambient level" verdict is
  preserved.

- **Round 2 R2.2 k-slot Tannakian-relevance probe (BORDERLINE):**
  the k-slot is TANNAKIAN-INVISIBLE on the ambient period ring,
  TANNAKIAN-RELEVANT on a candidate enriched functor. The present
  probe rules out the **pro-dg-category with Berezin transitions** as
  the candidate enrichment, sharpening the BORDERLINE finding to
  "the enrichment Q5' needs is NOT in the propinquity-pro-system
  direction; it must be sought elsewhere."

- **Round 2 deflation test (PARTIAL with DEFLATES bias):** the
  k-slot lives at the substrate level. The present probe confirms
  that the substrate-level information does NOT lift to cyclic
  homology via the Berezin pro-system, so no NC-motivic content
  reaches HP_*.

## Honest scope

1. **No fitted parameters introduced:** pure structural reasoning
   on published HP_* / Marcolli--Tabuada / Connes--Karoubi facts plus
   Paper 38 Berezin map properties.

2. **No PSLQ, no implementation, no code:** as scoped.

3. **No paper edits applied:** recommendation only, below.

4. **Sprint-pass depth caveat:** this probe was 1-day. A deeper
   examination might surface a non-standard cyclic-cohomology
   construction at the operator-system level that lifts the Berezin
   maps to functorial transitions. This was not exhaustively checked.
   The published-literature situation (no HP_* on operator systems
   yet constructed) suggests this is multi-year frontier; the TRIVIAL
   verdict at sprint-scale is the honest finding pending such future
   construction.

5. **What could NOT be checked at sprint-pass depth:**
   - Whether bivariant cyclic theory (Cuntz 2005) has a specific
     theorem ruling out non-trivial HP_* from CP-pro-systems with
     constant Morita type. The structural argument above suggests no
     non-trivial content can emerge, but a formal no-go theorem at
     this level was not located in single-session lit pass.
   - Whether the Camporesi--Higuchi Dirac operator equips
     $\mathcal{O}_{n_{\max}}$ with a non-trivial Chern character /
     index-class structure that bypasses the HP_* = $(\mathbb{C}, 0)$
     constraint. (The Connes--Karoubi character lives in HP_0 not
     HP_1 for the relevant index theorem on $S^3$, so this is
     unlikely to help, but was not fully checked.)
   - Whether the Connes--vS 2021 CMP truncated K-theory has a
     periodic-cyclic-style refinement that captures the M1/M2/M3
     partition. The published Connes--vS spectral truncations work
     does not engage HP_* directly at the cutoff level.

6. **The k-slot is, again, substrate-side metadata:** Round 2
   deflation memo identified this. The present probe confirms that
   the substrate-side metadata does not lift through HP_* of the
   Berezin pro-system. Any future Q5' attempt at the enriched
   functor would need to use a different categorical mechanism — not
   the GH-convergence-as-Tannakian-input reading.

7. **WH1 PROVEN is NOT re-opened by this finding.** Paper 38's
   GH-convergence theorem is a propinquity-convergence statement, not
   a cyclic-homology statement. The structural fact that propinquity
   convergence does not imply HP_* convergence is consistent with
   WH1 PROVEN being a fact about operator-norm metric convergence
   alone.

8. **Sprint A7's HALF-STRUCTURAL verdict is preserved.** The
   shared-field finding (M2 and M3 both engage $\mathbb{Q}(i)$) is
   purely a period-ring-level statement; the present probe says nothing
   about the period-ring shared field. It says only that the candidate
   enrichment cannot be the pro-dg-category with Berezin transitions.

9. **No combination rule $K = \pi(B + F - \Delta)$ framing edits.**
   Paper 2's "conjectural" label remains.

10. **Paper 18 §III.7 and Paper 32 §VIII are NOT modified.** The
    master Mellin engine partition and case-exhaustion theorem stand.

## Recommended paper edit (PI to apply, decline, or modify)

### Paper 55 §subsec:open_m2_m3 (Q5') sharpening

Round 2 R2.2 recommended adding one paragraph sharpening Q5' to the
Mellin-slot Tannakian-visibility question. Round 3 sharpens further
by ruling out one specific candidate enrichment (the
pro-dg-category with Berezin transitions). Suggested additional
paragraph:

> *Round 3 sharpening (Sprint Q5'-r3-pro-dg-category, June 2026; memo
> \texttt{debug/sprint\_q5p\_r3\_pro\_dg\_category\_memo.md}).*
> The candidate enrichment $\omega^{\mathrm{tri}}$ flagged in Round 2
> is NOT realised by the pro-dg-category
> $\varprojlim_{n_{\max}} \mathrm{dg}(\Tcal_{n_{\max}})$ with the
> Paper 38 Berezin reconstruction maps as transitions. Three
> structural reasons combine:\ (i) HP_* of each finite-$n_{\max}$
> cutoff is Morita-trivially $(\mathbb{C}, 0)$ (the truncated
> operator system spans the full matrix algebra
> $M_{N(n_{\max})}(\mathbb{C})$, whose periodic cyclic homology
> reduces to that of $\mathbb{C}$ by Morita invariance); (ii) the
> Berezin maps are completely positive but NOT multiplicative
> (property (c) of Paper 38 Lemma~\ref{lem:L4}, the
> approximate-identity error $\gamma_{n_{\max}}$ IS the
> multiplicativity defect), so they are not dg-functorial in the
> Marcolli--Tabuada sense; (iii) even granting a CP-extension of
> HP_* (bivariant cyclic theory, Cuntz--Quillen excision), the
> induced map factors through
> $\mathrm{HP}_*(M_N(\mathbb{C})) = (\mathbb{C}, 0)$ at every cutoff,
> destroying the non-trivial $\mathrm{HP}_1(C(\sthree)) = \mathbb{C}$
> (volume class) at every finite stage. The propinquity Berezin
> system measures operator-norm-metric convergence, not
> cyclic-homological convergence;\ the L2 quantitative-rate constant
> 4/π refines the metric convergence but not HP_*. Q5' therefore
> requires an enrichment of HP_* that lifts to operator systems and
> CP-maps — an extension not present in the published literature as
> of June 2026.

This is a recommendation only;\ no Paper 55 edits applied.

## Files used

### Memos read
- `debug/sprint_q5p_dim_sweep_memo.md` (Round 1 entry (i)).
- `debug/sprint_q5p_qsm_litread_memo.md` (Round 1 entry (ii)).
- `debug/sprint_q5p_tannakian_obstruction_memo.md` (Round 1
  entry (iii)).
- `debug/sprint_q5p_k_slot_tannakian_memo.md` (Round 2 R2.2).
- `debug/sprint_q5p_deflation_test_memo.md` (Round 2 deflation
  bonus).
- `debug/sprint_q5p_greenfield_marcolli_transport_memo.md` (Round 1
  Strand B depth-probe — relevant cross-reference).

### Papers read
- `papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex`
  (Lemma L4 Berezin reconstruction properties; central Fejér
  convolution form; L4(b) contractivity; L4(c) approximate identity;
  L2 rate $\gamma_{n_{\max}}$ asymptote 4/π).
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
  (§III truncated operator system not closed under multiplication,
  §VIII case-exhaustion theorem on Mellin engine).
- `papers/group3_foundations/paper_55_periods_of_geovac.tex`
  §subsec:open_m2_m3 (Q5' wording).

### Modules referenced
- `geovac/berezin_reconstruction.py` (BerezinReconstruction class,
  L4 properties at the operator level; approximate-identity residual,
  contractivity, positivity).

### Published references (no new bibitems required for Paper 55)
- Loday, *Cyclic Homology* (1998) §2.1.7 Morita invariance of HP_*;
  §2.5 $\Z/2$-graded structure of HP_*; §3.4.4 Connes--Karoubi
  isomorphism HP_*(C^∞(M)) ≅ H_dR^*(M).
- Marcolli--Tabuada arXiv:1110.2438 §3 dg-functorial requirement on
  morphisms.
- Cuntz, *Bivariant K-theory and the Weyl algebra*, K-theory 35
  (2005), the natural CP-extension of HP_*.
- Cuntz--Quillen arXiv:math/9311101 cyclic excision.
- Connes 1985 / Connes 1994 NCG book §IV.4 cyclic homology of AF
  algebras.

### Scripts (none required)
Pure structural reasoning;\ no code, no PSLQ, no numerical checks.

## One-line verdict

**TRIVIAL:** the pro-dg-category $\varprojlim \mathrm{dg}(\Tcal_{n_{\max}})$
with Berezin reconstruction maps as transitions does NOT carry
non-trivial Marcolli--Tabuada-fiber-functor input —
$\mathrm{HP}_*(\mathcal{O}_{n_{\max}}) = (\mathbb{C}, 0)$ by Morita
invariance at every cutoff, the Berezin maps are CP but not
multiplicative so they fail the dg-functorial requirement, and the
non-trivial $\mathrm{HP}_1(C(\sthree)) = \mathbb{C}$ volume class is
structurally inaccessible at every finite stage; the L2 quantitative
4/π rate refines propinquity convergence (operator-norm metric), not
cyclic-homological convergence, which is categorically orthogonal.
