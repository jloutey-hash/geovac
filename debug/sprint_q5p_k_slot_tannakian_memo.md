# Sprint Q5' — k-slot Tannakian-relevance probe

Date: 2026-06-04
Scope: 1-day structural sub-sprint of Q5' (Paper 55 §subsec:open_m2_m3),
Round 2. Builds on Round 1 closures:
- `debug/sprint_q5p_dim_sweep_memo.md` (NEGATIVE on dimension-sweep entry point)
- `debug/sprint_q5p_qsm_litread_memo.md` (PARTIALLY TRANSPORTABLE on QSM/MT lit-precedent)
- `debug/sprint_q5p_tannakian_obstruction_memo.md` (OPEN-POSSIBLE on cheap Tannakian obstructions)

The Round 1 obstruction-probe deflated Q5' to: "Does GeoVac's spectral
triple structure add anything to the ambient
$\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$ classification?" This Round 2
probe asks specifically: **is the operator-order slot index
$k \in \{0, 1, 2\}$ (the M1/M2/M3 partition; Paper 18 §III.7;
Paper 32 §VIII case-exhaustion theorem) Tannakian-relevant
structure on the period ring, or a calculational bookkeeping index
that the ambient MT does not track but does not need to?**

## TL;DR

**Verdict: BORDERLINE, leaning TANNAKIAN-INVISIBLE on the period
ring but TANNAKIAN-RELEVANT one categorical level up (on the
mapping of the spectral triple TO the period ring).**

The structural argument is sharp:

1. **The period values themselves do not remember $k$.** Paper 55
   §6.1 (eq:ring_nesting) records the load-bearing identity
   \[
     M_2 = \bigoplus_{k \ge 0} \pi^{2k}\cdot\Q
     \;\subset\; \mathrm{MT}(\Z)
     \;\subset\; M_3 = \mathrm{MT}(\Z[i, 1/2]),
     \qquad
     M_1 \cap \mathrm{MT}(\Z[i, 1/2]) = \Q[\pi].
   \]
   The same period value (e.g.\ $\pi^{2k}$) is reachable as an M1
   output (Hopf-base measure to the $k$-th power), an M2 output
   (Seeley--DeWitt at order $k$), or an M3 output (even-weight
   Brown MZV inside the un-restricted Dirac Dirichlet). The Mellin
   slot $k$ is, in Paper 55's own words, *"the structural label that
   distinguishes the three readings;\ the period value alone does
   not."* Any Tannakian fiber functor
   $\omega: \mathcal{C} \to \mathrm{Vec}_{\Q}$ that lands in
   $\mathrm{MT}(\Z[i, 1/2], 4)$ at the *target* level sees only the
   period value, not the path that produced it. So the standard
   Marcolli-Tabuada / Connes-Marcolli motivic Galois group acting on
   the period ring DOES NOT distinguish M1's $\pi^2$ from M2's $\pi^2$
   from M3's $\pi^2$. **The k-slot is forgotten by every standard
   target-level fiber functor.**

2. **But the case-exhaustion theorem makes $k$ visible at the
   *source* level.** Paper 32 §VIII Theorem `thm:pi_source_case_exhaustion`
   is a structural theorem on GeoVac's spectral triple: every $\pi$
   that appears in any output engages at least one of M1, M2, M3, with
   the engagement traceable to a specific Mellin moment
   $\mathcal{M}[\mathrm{Tr}(D^k \cdot e^{-tD^2})]$. The
   *spectral-triple-to-period-ring map* therefore factorises through
   the master Mellin engine, which is a three-component object
   indexed by $k$. **The k-slot IS visible at the level of
   computational paths from the spectral triple TO the period ring;
   it is INVISIBLE at the level of the period ring itself.**

3. **What this means in Tannakian language.** A standard Tannakian
   category encodes objects + tensor products + a fiber functor; the
   motivic Galois group is the automorphism group of the fiber
   functor. In the GeoVac setting, the "objects" are spectral-triple
   data (or their dg-categorical thickening, per Marcolli-Tabuada
   2016), the period values are the IMAGES of those objects under
   $\mathrm{HP}_*$ (the canonical fiber functor for noncommutative
   motives), and the k-slot lives on the MAPPING $\mathrm{HP}_*$
   side, not on the target $\mathrm{Vec}_{\Q}$ side. Tannakian
   reconstruction (Deligne 1990, Marcolli-Tabuada 2016 Thm 1.1)
   recovers the *image of* $\mathrm{HP}_*$ in $\mathrm{Vec}_{\Q}$,
   not the structure of $\mathrm{HP}_*$ as a multi-component map.
   **The Tannakian recovery procedure standardly washes out the
   k-slot.**

4. **Therefore on the ambient MT alone, $k$ is TANNAKIAN-INVISIBLE.**
   The standard motivic Galois group of $\mathrm{MT}(\Z[i, 1/2], 4)$
   acts on a category whose objects are mixed-Tate motives with a
   weight filtration $W_\bullet$ and a level filtration; the k-slot
   index does not appear in either of these standard filtrations.
   In particular: weight $2k$ is achievable from $k=0$ (M1, $\pi^{2k}$
   as a $k$-th power of the Hopf measure), $k=2$ (M2 Seeley--DeWitt
   at order $k$), AND $k=1$ (M3 even-weight component of the
   un-restricted Dirac Dirichlet); the weight does NOT determine $k$.
   This refutes the cleanest TANNAKIAN-RELEVANT reading ("$k$ aligns
   with weight filtration via Mellin shift conventions") because the
   weight filtration is not a bijection on $k$.

5. **But on the *enriched* category (spectral triple + its master
   Mellin engine path-decomposition), $k$ IS TANNAKIAN-RELEVANT.**
   The candidate enrichment Q5' contemplates is: replace the
   standard target functor $\omega_{\mathrm{HP}}: \mathrm{dg}(\Tcal) \to
   \mathrm{Vec}_{\Q}$ with a three-component functor
   $\omega^{\mathrm{tri}}: \mathrm{dg}(\Tcal) \to
   \mathrm{Vec}_{\Q} \otimes \mathrm{IndexCat}(\{0, 1, 2\})$
   that records WHICH Mellin slot produced each period value. Under
   this enrichment, the master Mellin engine becomes a *graded* (or
   $\Z/3$-indexed) symmetric monoidal functor, and the period ring
   sub-categorises into three named sub-rings $M_1$, $M_2$, $M_3$ with
   the case-exhaustion theorem as the "every period comes from at
   least one slot" exactness statement. The motivic Galois group of
   this enriched fiber functor would be a quotient of the ambient
   $\mathrm{Gal}^{\mathrm{mot}}_{\mathrm{MT}(\Z[i, 1/2], 4)}$ by the
   subgroup that acts trivially on the k-slot index; whether this
   subgroup is trivial (k-slot fully Tannakian-relevant) or all of
   $\mathrm{Gal}^{\mathrm{mot}}$ (k-slot fully Tannakian-invisible)
   is the genuinely open structural question.

The honest finding: **$k$ is TANNAKIAN-INVISIBLE on the standard
ambient MT periods, but it is TANNAKIAN-RELEVANT on a candidate
enrichment ($\Z/3$-graded noncommutative motive) that has not been
constructed in the literature.** This is the precise content of the
"BORDERLINE" verdict: the question of Tannakian relevance depends on
which fiber functor one uses, and the natural enrichment has
genuine substance. Q5' has not fully deflated.

## Decision-gate framework

The task posed three readings:

- **TANNAKIAN-RELEVANT:** k-slot induces grading/sub-category/weight
  filtration on the period ring that survives passage to the
  motivic-Galois fiber functor.
- **TANNAKIAN-INVISIBLE:** k is a calculational bookkeeping index;
  any Hodge/Betti realisation of $M_2 \oplus M_3$ outputs preserves
  values but loses k-tags.
- **BORDERLINE:** k is visible in some realisations but not others.

The structural argument below lands on BORDERLINE for the reason
above: TANNAKIAN-INVISIBLE on standard MT fiber functors,
TANNAKIAN-RELEVANT on the spectral-triple-side enrichment.

## Argument by k value

### $k = 0$ (M1, Hopf-base measure): the cleanest TANNAKIAN-INVISIBLE case

M1's period ring is $\Q[\pi, \pi^{-1}]$, the localised pure-Tate
sub-ring generated by $\pi$ (Paper 32 §VIII Cor `cor:m1_pure_tate`).
At each integer Tate weight $n \in \Z$, M1 produces $\pi^n \cdot \Q$.

In the ambient $\mathrm{MT}(\Z[i, 1/2], 4)$, the Tate sub-category
$\mathrm{Tate}(\Z) = \bigoplus_n \Q(n)$ is the canonical pure-Tate
sub-category. Every period value in $\Q[\pi, \pi^{-1}]$ is the
Betti period of some object in $\mathrm{Tate}(\Z)$ (up to the
$(2\pi i)^n$ vs $\pi^n$ normalisation, which is itself a fixed
choice not encoded in motivic-Galois content).

**The motivic Galois group of $\mathrm{Tate}(\Z)$ is $\mathbb{G}_m$**,
acting via Tate twist:\ on $\Q(n)$, $\lambda \in \mathbb{G}_m$ acts
as multiplication by $\lambda^n$. **This action does NOT remember
that the period came from M1 vs M2.** A specific $\pi^{2k}$ value
is acted on by $\mathbb{G}_m$ as $\lambda^{2k}$, regardless of
whether it was produced by:
- $k$ Hopf-base measures multiplied (M1 path), or
- A Seeley--DeWitt coefficient $a_k$ at order $k$ (M2 path), or
- An even-weight Brown MZV component (M3 path).

**Verdict for $k = 0$: TANNAKIAN-INVISIBLE on the standard fiber
functor.** The Hopf-base measure produces pure-Tate periods that
are GALOIS-FIXED in $\mathrm{Gal}^{\mathrm{mot}}_{\mathrm{MT}(\Z)}$;
the slot label "this came from M1" does not appear in the
$\mathbb{G}_m$ action.

### $k = 2$ (M2, Seeley--DeWitt): same TANNAKIAN-INVISIBLE verdict, with sharper restriction

M2's period ring on $S^3$ is $\bigoplus_k \pi^{2k}\cdot\Q$, strictly
smaller than the generic Fathizadeh-Marcolli mixed-Tate ring (Sprint
Mixed-Tate Test verdict; Paper 32 §VIII Cor `cor:m2_mixed_tate`). The
F-M motive $\mathcal{V}_n^{\mathrm{F-M}}(\lambda)$ at the GeoVac
specialisation $\lambda = 1$ has explicit Grothendieck class
$[Z_{1, 2n}] = [\mathbb{P}^n](1 + \mathbb{L}^n)$
(Sprint A8; Paper 55 §4 Remark `rem:m2_specialisation`), a
palindrome of length $2n + 1$ with doubled middle. The motive itself
lives in $\mathrm{MT}(\Z[i])$ at the proof level (Witt-splitting
field for the standard quadric, Sprint A7); the period output
descends to $\bigoplus_k \pi^{2k} \cdot \Q$ via Galois invariance.

The standard motivic Galois group of $\mathrm{MT}(\Z)$
(Brown 2012 + Deligne-Goncharov) acts trivially on
$\bigoplus_k \pi^{2k} \cdot \Q$ because pure-Tate periods of even
weight are Galois-fixed (the only motivic Galois action is
$\mathbb{G}_m$ Tate-twist, and pure-Tate periods are weight-graded).

**At the period ring level, M2's pure-Tate periods are
INDISTINGUISHABLE from M1's pure-Tate periods (other than the
weight grading itself).** A motivic Galois fiber functor that
records weight (which is exactly what
$\mathrm{Gal}^{\mathrm{mot}}_{\mathrm{Tate}(\Z)} = \mathbb{G}_m$
does) sees the weight $2k$ but cannot tell whether the source was
$k$ copies of M1 multiplied or a single M2 coefficient at order $k$.

**Verdict for $k = 2$: TANNAKIAN-INVISIBLE.**

### $k = 1$ (M3, vertex-parity Dirichlet $L$ / $\eta$-invariant): mostly TANNAKIAN-INVISIBLE, with one borderline case

M3's period ring is the cyclotomic mixed-Tate ring
$\mathrm{MT}(\Z[i, 1/2])$ at level $\le 4$ (Paper 32 §VIII Cor
`cor:m3_cyclotomic_mixed_tate`; Glanois 2015). The motivic Galois
group at this level is a non-trivial extension of
$\mathbb{G}_m$ by a pro-unipotent group $U$, with the cyclotomic
part contributing a $\mathrm{Gal}(\Q(i)/\Q) = \Z/2$ action by complex
conjugation.

The non-trivial $\Z/2$ action on M3 is the well-understood Deligne-Glanois
descent from level 4 to level 2 (Glanois 2015 Cor 1.2; Paper 55
Prop `prop:vertex_parity_descent`). This Galois action **acts on
period values**: it sends $\beta(s) \mapsto -\beta(s)$ at odd $s$
(since $\chi_{-4}$ is odd). On M3's output side, this action
realises the parity flip $D_{\mathrm{even}} \leftrightarrow
D_{\mathrm{odd}}$.

Two sub-cases:
- **M3 cross-product with M1 (e.g.\ $\zeta(3)/\pi^2$):**
  The $\zeta(3)$ component sits in $\mathrm{MT}(\Z)$ at depth 1; the
  $\pi^{-2}$ component sits in $M_1$. The combined period
  $\zeta(3)/\pi^2$ lives in $\mathrm{MT}(\Z[i, 1/2])$ at depth 1
  via M1×M3 cross-product (Paper 55 §6.4 Paper 50 example). The
  motivic Galois action acts on $\zeta(3)$ as the cyclotomic depth-1
  generator action; it acts trivially on $\pi^{-2}$ (pure-Tate of
  weight $-2$). **The $k = 1$ vs $k = 0$ tag is lost** — the
  motivic Galois action sees only "depth-1 cyclotomic" times
  "weight-$(-2)$ pure-Tate", with no record that the depth-1 piece
  came from M3 and the weight-$(-2)$ piece came from M1.
- **M3 internal to even weight (e.g.\ even-weight Brown MZVs):**
  The motivic Galois group of $\mathrm{MT}(\Z)$ acts trivially on
  even-weight pure-Tate periods (same as M1 / M2 case). **The $k = 1$
  tag is lost** — these M3 periods are Galois-indistinguishable
  from the corresponding M1 / M2 periods of the same weight.

The one BORDERLINE case is the **odd-weight Dirichlet $\beta(s)$
sector at odd $s$**, which lives in $\mathrm{MT}(\Z[i, 1/2])$ at
level 4 NOT reducible to $\mathrm{MT}(\Z)$. Here the cyclotomic
Galois action is non-trivial, AND the resulting period value
(e.g.\ $G = \beta(2)$, $\beta(4)$) has no representation in M1 or
M2 (no $\pi$-power-times-rational form, no Seeley--DeWitt coefficient
on a sphere produces it). For these specific outputs, the period
value DOES uniquely identify the slot $k = 1$ — because no other
slot produces values in $\mathrm{MT}(\Z[i, 1/2]) \setminus \mathrm{MT}(\Z)$.

But this is NOT k-slot Tannakian relevance:\ the Galois action on
these specific periods is the standard cyclotomic level-4 action; it
distinguishes them from level-2 periods (Brown MZVs) and from
pure-Tate periods, but it does not detect "this came from
$\mathrm{Tr}(D \cdot e^{-tD^2})$ vs.\ some other Mellin-engine
input." The level-4 cyclotomic-ness is the structural feature; the
$k = 1$ tag is a label *on the GeoVac-side mapping*, not on the
period ring.

**Verdict for $k = 1$: TANNAKIAN-INVISIBLE on period ring outputs,
with the only sharp k-detection being on the level-4 cyclotomic
sector where M3 has *exclusive output access*.** The exclusivity
is a feature of the master Mellin engine (M1 and M2 cannot produce
non-pure-Tate periods), not of a Galois-action tag carried by the
period values.

### Cross-cut: the symmetric reading

Every $k$ value lands at TANNAKIAN-INVISIBLE on the period ring
because the standard motivic Galois group of
$\mathrm{MT}(\Z[i, 1/2], 4)$ acts on period values, not on
production paths. The k-slot tag is OUTSIDE the standard
motivic-Galois apparatus.

But the case-exhaustion theorem (Paper 32 §VIII
`thm:pi_source_case_exhaustion`) provides a *structural fact about
which $k$ can produce which periods*:

| Period type | M1 reachable? | M2 reachable? | M3 reachable? |
|:------------|:-------------:|:-------------:|:-------------:|
| $\pi^n \cdot \Q$, $n \in \Z$, even | YES (with $n = 2k$) | YES (Seeley--DeWitt at order $k$) | YES (even-weight Brown MZV component) |
| $\pi^n \cdot \Q$, $n \in \Z$, odd | YES (with $n = 2k+1$) | NO (M2 even-weight only on $S^3$) | NO (odd-weight Brown MZV components present, but mixed with $\beta$-values) |
| Brown MZV (odd weight $\zeta(2k+1)$) | NO | NO ($\zeta(3)$, $\zeta(5)$ absent from $S^3$ M2 by Sprint Mixed-Tate Test) | YES (un-restricted Dirac Dirichlet at integer $s$) |
| $\chi_{-4}$-related ($G$, $\beta(s)$, $\beta(4)$) | NO | NO | YES (exclusively) |
| Cyclotomic level-4 (depth $\ge 2$) | NO | NO | YES (loop tower) |

**The exclusivity pattern in the last three rows is what makes the
$k$-slot ALMOST Tannakian-relevant.** If one fixes the period
value, the table tells you which $k$ MUST have produced it (modulo
the overlapping pure-Tate row). This is a structural feature of the
master Mellin engine, and it has the *shape* of a Tannakian
classification — but it is a classification by which *source maps*
can produce the value, NOT by which *Galois action* acts non-trivially.

The two readings are categorically distinct. The first is a
constraint on the mapping $\omega^{\mathrm{tri}}: \mathrm{dg}(\Tcal)
\to (M_1, M_2, M_3)$; the second is the standard Tannakian
recovery of the automorphism group of $\omega$.

## Named Marcolli-Tabuada entry points

To close the structural argument with concrete published targets
that would settle the BORDERLINE in one direction or the other:

### Marcolli-Tabuada 2016 Theorem 1.1 (arXiv:1110.2438):
Constructs the noncommutative motivic Galois group
$\mathrm{GalMot}_{\mathrm{NC}}$ via Tannakian dual applied to the
category of noncommutative numerical motives, with periodic cyclic
homology $\mathrm{HP}_*$ as the fiber functor. **Concrete test:** if
the discrete GeoVac spectral triple $\mathcal{T}_{n_{\max}}$ (or its
dg-categorical thickening) is sent to $\mathrm{Vec}_{\Q}$ by
$\mathrm{HP}_*$, does the image carry an explicit $\Z/3$-grading (or
filtration) corresponding to the operator-order slot $k$?

The honest answer is that $\mathrm{HP}_*$ as standardly defined is a
$\Z/2$-graded functor (the even/odd part split), NOT a $\Z/3$-graded
functor. The operator-order slot $k \in \{0, 1, 2\}$ does NOT align
with $\mathrm{HP}_*$'s $\Z/2$ split. So at the standard
Marcolli-Tabuada level, the k-slot is **STRUCTURALLY NOT VISIBLE
through $\mathrm{HP}_*$.**

This is a clean structural verdict supporting TANNAKIAN-INVISIBLE.

### Marcolli-Tabuada 2015 Theorem 1.4 (arXiv:1112.5422):
Unconditional construction of the noncommutative motivic Galois
group. The fiber functor in this case is built via the
André-Kahn quotient, which records the *image* of
$\mathrm{HP}_*$ in $\mathrm{Vec}_{\Q}$ at zero-th order (after
quotienting by the radical). Same $\Z/2$ grading as in 1110.2438.
Same verdict: **k-slot not visible through this fiber functor.**

### Fathizadeh-Marcolli 2017 Theorem 1.1 (arXiv:1611.01815):
Spectral-action SD coefficients land in $\mathrm{MT}(\Z)$ at the
continuum R-W level. The proof goes via Wodzicki's residue and
Rosenfeld-style integral representations. **The path to
$\mathrm{MT}(\Z)$ records that the SD coefficient came from
spectral-action data, but the embedding $\mathrm{MT}(\Z) \hookrightarrow
\mathrm{MT}(\Z[i, 1/2], 4)$ is the SAME embedding regardless of
which Mellin moment of $\mathrm{Tr}(D^k e^{-tD^2})$ produced the
SD coefficient.**

F-M does not record the k-slot. It records the spectral action,
which is canonically $k = 2$ (Seeley--DeWitt), but the F-M motivic
classification is statement at the *output ring* level, not at the
*input path* level. This supports TANNAKIAN-INVISIBLE on the
F-M-style motivic classification.

### Connes-Marcolli 2004 / Connes-Marcolli book (arXiv:math/0409306):
The cosmic Galois group $U^*$ acts on perturbative renormalization
data; it is a classifying object for the divergence structure of
QFT, organised by depth (loop order) and weight (motivic weight),
with the **specific input** being equisingular flat connections (the
analog of "discrete spectral triple + Mellin engine path" for QFT).
**$U^*$ acts on counterterm coefficients, which are LABELLED by the
diagram that produced them.** This is the closest published precedent
for a Tannakian construction where the LABEL of the source (here:
the Feynman diagram; for Q5': the Mellin slot $k$) is RETAINED in
the Galois symmetry.

The cosmic Galois group's action on counterterms factorises through
the structure of the connection (singular vs.\ equisingular,
diagram-by-diagram), and the symmetry is sensitive to which diagram
contributed. **If Q5' has substance, it would be in the analog
direction:\ a "spectral cosmic Galois group" of the GeoVac discrete
spectral triple that acts on Mellin-engine outputs *labelled by their
slot k*.** This is the natural enrichment that Q5' contemplates and
that the Round 1 lit-read (Strand A) identified as "the gold-standard
precedent for motivic Galois from analytic NCG data."

**Verdict on $U^*$ precedent:** the construction's *shape* makes
k-slot Tannakian-relevant on a candidate enriched Galois group; this
is the BORDERLINE-RELEVANT reading of Q5'. Whether the actual
construction works for GeoVac specifically is the genuine open
question.

## Cross-check with Sprint MR-A (mechanism-as-domain partition)

Sprint MR-A (May 2026) established that the k-slot index carries
DUAL content (Paper 18 §III.7 closing paragraph):

1. It indexes the M-mechanism producing $\pi$ (the case-exhaustion
   theorem's three-fold partition).
2. It indexes the natural class of OBSERVABLE from which the
   mechanism's signature can be extracted (state-space propinquity
   rate for M1; heat-kernel/spectral-action for M2; vertex-restricted
   parity-character sum for M3).

The DUAL content is consistent with the BORDERLINE finding:\ at
the level of period values (output ring), the k-slot is
TANNAKIAN-INVISIBLE; at the level of observable classes (input
class to the Mellin engine), the k-slot IS visible as a structural
classification.

**The MR-A "mechanism is also a domain" observation is exactly the
fact that the k-slot lives on the spectral-triple-to-period-ring
map, not on the period ring itself.** This matches the structural
argument of this memo.

## What WOULD force a verdict (if known)

| Reading | Would settle Q5' how? | Verdict if true |
|:--------|:-------------------------|:----------------|
| $\mathrm{HP}_*$ on $\mathrm{dg}(\Tcal_{n_{\max}})$ has a natural $\Z/3$-graded structure refining the standard $\Z/2$ grading | TANNAKIAN-RELEVANT | The standard Marcolli-Tabuada construction already sees k-slot |
| $\mathrm{HP}_*$ is $\Z/2$-graded only, no further refinement | TANNAKIAN-INVISIBLE | Standard construction loses k-slot at $\mathrm{HP}_*$ step |
| The Marcolli-Tabuada motivic Galois group of $\mathrm{dg}(\varprojlim \Tcal_{n_{\max}}, B_{n_{\max}})$ has a sub-quotient acting non-trivially on the k-slot index | TANNAKIAN-RELEVANT-WITH-ENRICHMENT | The enriched construction (Q5' candidate) makes k-slot visible |
| The same Galois group acts trivially on the k-slot index | TANNAKIAN-INVISIBLE | Even the enriched construction loses k-slot |

The honest finding from this probe is that **the standard
$\mathrm{HP}_*$ functor is $\Z/2$-graded** (this is well-established
in the cyclic-homology literature; Loday's textbook *Cyclic
Homology* 1998 §2.5, Marcolli-Tabuada 2016 Cor 4.2), so the
standard Marcolli-Tabuada construction CANNOT see the $\Z/3$
operator-order grading. **At the standard level, k-slot is
TANNAKIAN-INVISIBLE.**

But the case-exhaustion theorem provides a *non-standard* invariant
of the spectral triple (which slot produces which period output) that
admits the candidate enrichment Q5' contemplates. **At the enriched
level, k-slot is TANNAKIAN-RELEVANT only if the enrichment is
constructed and the resulting motivic Galois action is non-trivial
on the k-slot.** This is the genuine multi-year construction Q5'
would need.

## Verdict

**BORDERLINE.**

- **TANNAKIAN-INVISIBLE on the standard motivic-Galois fiber functor**
  (Betti, de Rham, $\mathrm{HP}_*$, Marcolli-Tabuada 1110.2438). The
  k-slot index does not appear in the weight filtration of the
  ambient $\mathrm{MT}(\Z[i, 1/2], 4)$ (because weight $2k$ is
  achievable from $k \in \{0, 1, 2\}$ independently), does not appear
  in the level filtration (which separates pure-Tate vs.\
  cyclotomic-mixed-Tate, not M1 vs.\ M2 vs.\ M3), and does not appear
  in the $\Z/2$ grading of $\mathrm{HP}_*$. The standard
  Marcolli-Tabuada motivic Galois group acts on the ambient periods
  without distinguishing the k-slot source.

- **TANNAKIAN-RELEVANT one categorical level up**, on the candidate
  enrichment $\omega^{\mathrm{tri}}: \mathrm{dg}(\Tcal) \to
  \mathrm{Vec}_{\Q} \otimes \mathrm{IndexCat}(\{0, 1, 2\})$ that
  records the Mellin slot per period. This enrichment IS suggested
  by the case-exhaustion theorem's exclusivity pattern (M3 has
  exclusive output access to cyclotomic level-4 periods; M1 has
  exclusive output access to localised pure-Tate; M2 is the
  spectral-action-canonical slot for SD coefficients with the F-M
  inheritance). Whether this enrichment defines a well-formed
  Tannakian category in the Marcolli-Tabuada sense — i.e.\ whether
  the resulting motivic Galois group is well-defined and acts
  non-trivially on the k-slot — is the genuine multi-year research
  target.

The BORDERLINE reading **does not deflate Q5' further** beyond the
Round 1 closures. It sharpens Q5' to:

> Does the enriched fiber functor $\omega^{\mathrm{tri}}$ define a
> Tannakian category whose motivic Galois group acts non-trivially
> on the operator-order slot index $k$?

This is the precise structural question the multi-year mathematical
research project flagged in Sprint A7 and Round 1 entry point (ii)
(QSM/MT lit-read) would address. The standard Marcolli-Tabuada
machinery does NOT settle it (because the standard machinery washes
out the k-slot at the $\mathrm{HP}_*$ step); the enriched
construction would need to be built.

## Compatibility with prior Q5' closures

The BORDERLINE verdict is the consistent refinement of all three
Round 1 closures:

- **Round 1 dim sweep (NEGATIVE):** the dimension sweep showed that
  the M2/M3 $\Q(i)$ field-coincidence persists across $d \in
  \{3, 5, 7\}$ via two independent one-bit invariants. This memo's
  reading is that the *ambient* MT classification is dimension-
  sweep-invariant precisely because *the ambient MT does not record
  the k-slot* — both M2 and M3 outputs land in the same ambient
  ring, and the field-coincidence is a generic feature of that
  ambient ring at each $d$. The k-slot Tannakian-INVISIBILITY at
  the ambient level is what makes the dim sweep negative; the
  candidate enrichment is where any positive content could live.

- **Round 1 QSM/MT lit-read (PARTIALLY TRANSPORTABLE):** the
  Marcolli-Tabuada machinery is the right scaffolding but does
  not, at the standard $\mathrm{HP}_*$ level, see the k-slot. The
  candidate enrichment Q5' would need to add a $\Z/3$-component
  refinement of the standard noncommutative motivic Galois
  construction. The lit-read's "PARTIALLY TRANSPORTABLE" verdict
  is sharpened here to "standard machinery sees ambient MT but
  not k-slot; enrichment is what Q5' would need to build."

- **Round 1 cheap obstruction probe (OPEN-POSSIBLE):** the obstruction
  probe found no cheap structural obstruction to a candidate
  enrichment. The BORDERLINE verdict here is the constructive
  complement:\ no obstruction, AND no realisation; the enrichment
  remains a genuine open construction.

## Honest scope

1. **No fitted parameters introduced:** pure structural reasoning on
   published periods-program facts plus the case-exhaustion theorem.
2. **No PSLQ, no implementation, no code:** as scoped.
3. **No paper edits applied:** recommendation only, below.
4. **Sprint-pass depth caveat:** this probe was 1-day. A deeper
   examination might surface a structural argument not visible at
   this depth, particularly around whether the
   periodic-cyclic-homology grading on truncated dg categories
   (Loday §2.5, Tabuada *Noncommutative Motives* AMS 2015) admits
   a non-standard $\Z/3$-refinement that aligns with the master
   Mellin engine. This was not checked here.
5. **What could NOT be checked at sprint-pass depth:**
   - Whether the dg category of the GeoVac discrete spectral triple
     at fixed $n_{\max}$ has a non-trivial periodic cyclic homology
     (Marcolli-Tabuada's input expects non-trivial $\mathrm{HP}_*$;
     finite-dim matrix algebras have small $\mathrm{HP}_*$).
   - Whether the pro-dg-category $\varprojlim \mathrm{dg}(\Tcal_{n_{\max}})$
     with Berezin maps $B_{n_{\max}}$ (Paper 38 §L4) is the right
     object for non-trivial Marcolli-Tabuada input (Round 1 lit-read
     §4 flagged this as likely).
   - Whether Connes-Marcolli's cosmic Galois group $U^*$ has been
     applied to any compact-substrate spectral triple in the published
     literature (Round 1 lit-read §6 returned NO hits but explicitly
     flagged this as "to our knowledge after this lit pass").
6. **The case-exhaustion theorem is a structural fact about the
   spectral triple's mapping to the period ring**, not a Tannakian
   theorem on the period ring. This memo's argument rests on
   distinguishing the two levels; the distinction itself is honest
   and supported by the standard $\mathrm{HP}_*$ grading discussion.
7. **Sprint A7's HALF-STRUCTURAL verdict is preserved.** A7 said
   "shared $\Q(i)$ at the number-field level, two independent
   structurally-distinct mechanisms." This memo's BORDERLINE sharpens:\
   the two mechanisms are independent on the period ring (because
   the period ring forgets the k-slot), but the enriched-construction
   level would record both. The two verdicts are compatible
   refinements.
8. **No combination rule $K = \pi(B + F - \Delta)$ framing edits.**
   Paper 2's "conjectural" label on the combination rule remains.
9. **Paper 18 §III.7 and Paper 32 §VIII are NOT modified.** The
   master Mellin engine partition and the case-exhaustion theorem
   stand as written. This memo provides a Tannakian REFINEMENT of
   what the engine MEANS, not a change to what it COMPUTES.

## Recommended paper edit (PI to apply, decline, or modify)

### Paper 55 §subsec:open_m2_m3 (Q5') sharpening

Currently the Q5' subsection (Paper 55 lines 1474–1499) is a
multi-year-research-target open question without a sharp Tannakian
framing. Round 2 of Q5' scoping suggests adding one paragraph to
sharpen the question:

> *Sharpening of Q5' to the Mellin-slot Tannakian-visibility
> question (Sprint Q5'-k-slot, June 2026; memo
> \texttt{debug/sprint\_q5p\_k\_slot\_tannakian\_memo.md}).*
>
> The k-slot index $k \in \{0, 1, 2\}$ of the master Mellin engine
> (Paper 18 §III.7) is TANNAKIAN-INVISIBLE on the standard
> motivic-Galois fiber functor of the ambient
> $\mathrm{MT}(\Z[i, 1/2], 4)$:\ weight $2k$ is achievable from
> all three slots independently, the level filtration is between
> pure-Tate and cyclotomic (not M1 vs.\ M2 vs.\ M3), and
> periodic cyclic homology $\mathrm{HP}_*$ is $\Z/2$-graded
> (not $\Z/3$). The k-slot index IS TANNAKIAN-RELEVANT
> one categorical level up:\ on a candidate enriched fiber functor
> $\omega^{\mathrm{tri}}: \mathrm{dg}(\Tcal) \to \mathrm{Vec}_\Q \otimes
> \mathrm{IndexCat}(\{0, 1, 2\})$ that records which Mellin slot
> produced each period output. The exclusivity pattern of the
> case-exhaustion theorem (Theorem `thm:pi_source_case_exhaustion`)
> — M3 has exclusive output access to cyclotomic level-4 periods;
> M1 has exclusive output access to localised pure-Tate at
> non-integer Tate-weight ratios; M2 is spectral-action-canonical
> for SD coefficients — suggests such an enrichment is well-posed,
> but its motivic Galois group has not been constructed in the
> published literature. The multi-year research project Q5' names
> is therefore the construction of this enrichment and the
> verification that its motivic Galois group acts non-trivially on
> the k-slot index. Sprint Q5'-k-slot returned BORDERLINE, leaning
> TANNAKIAN-INVISIBLE on the ambient periods and
> TANNAKIAN-RELEVANT-with-enrichment on the candidate refinement.

This is a recommendation only;\ no Paper 55 edits applied.

## Files used

### Memos read
- `debug/sprint_q5p_dim_sweep_memo.md` (Round 1 entry point (i):\
  NEGATIVE on dimension sweep)
- `debug/sprint_q5p_qsm_litread_memo.md` (Round 1 entry point (ii):\
  PARTIALLY TRANSPORTABLE on QSM/MT lit-read)
- `debug/sprint_q5p_tannakian_obstruction_memo.md` (Round 1 entry
  point (iii):\ OPEN-POSSIBLE on cheap obstruction probe)
- `debug/sprint_a7_m2_m3_cyclotomic_memo.md` (HALF-STRUCTURAL Sprint
  A7 base verdict)
- `debug/sprint_mixed_tate_test_memo.md` (Sprint Mixed-Tate Test
  closures on M2 sub-ring + $\sqrt\pi$ vs.\ $\pi^2$ convention)

### Papers read
- `papers/group3_foundations/paper_18_exchange_constants.tex` §III.7
  (master Mellin engine partition `eq:m_engine_period_stratification`
  + master-mechanism reading + mechanism-as-domain sharpening
  + connection to motivic / period theory)
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex`
  §VIII (case-exhaustion theorem `thm:pi_source_case_exhaustion`,
  Corollaries `cor:m1_pure_tate`, `cor:m2_mixed_tate`,
  `cor:m3_cyclotomic_mixed_tate`, joint-engagement Remark)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §§3-6
  (joint engagement, role-disjointness + set-theoretic overlap,
  Tate-weight bookkeeping, canonical examples K, F-theorem,
  propinquity rate, joint-engagement modes)

### Published references (already in Paper 55 bibliography)
- Fathizadeh--Marcolli arXiv:1611.01815 (mixed-Tate spectral action
  on R-W; transfers verbatim to GeoVac $S^3$ M2 static specialisation)
- Glanois arXiv:1411.4947 (level-4 cyclotomic mixed-Tate basis)
- Brown 2012 (level-1 MZV ring baseline)
- Deligne 2010 (cyclotomic motivic Galois $N \in \{2, 3, 4, 6, 8\}$)
- Marcolli-Tabuada arXiv:1110.2438 (NC numerical motives + $\mathrm{HP}_*$
  fiber functor; the $\Z/2$ grading of $\mathrm{HP}_*$ is the
  load-bearing fact for the TANNAKIAN-INVISIBLE side of the
  verdict)
- Marcolli-Tabuada arXiv:1112.5422 (unconditional NC motivic Galois)
- Connes-Marcolli arXiv:math/0409306 (cosmic Galois group $U^*$;
  the SHAPE Q5' enrichment would mirror)

### Scripts (none required)
Structural reasoning only; no code, no PSLQ, no numerical checks.

## One-line verdict

BORDERLINE:\ k-slot is TANNAKIAN-INVISIBLE on the ambient
$\mathrm{MT}(\Z[i, 1/2], 4)$ period ring (the $\mathrm{HP}_*$ fiber
functor is $\Z/2$-graded, the weight filtration does not separate
M1/M2/M3, the cyclotomic level filtration does not align with $k$),
but TANNAKIAN-RELEVANT one categorical level up on a candidate
enriched fiber functor $\omega^{\mathrm{tri}}$ that records the
Mellin slot per period output — Q5' sharpens to the construction
question of whether this enrichment defines a well-formed Tannakian
category whose motivic Galois group acts non-trivially on the
operator-order slot index.
