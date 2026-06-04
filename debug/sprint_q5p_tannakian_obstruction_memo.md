# Sprint Q5' — Tannakian obstruction probe for M2 / M3 ℤ/2 unification

Date: 2026-06-03
Scope: 1-day diagnostic sub-sprint of Q5' (Paper 55 §subsec:open_m2_m3).
Entry point (iii) of the three-way Q5' scoping recommendation: a
structural Tannakian obstruction probe, run *before* a Phase 1 lit
review (entry point (ii)) is justified. The question is whether any
Tannakian extension containing both M2's pure-Tate output ring
$\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$ and M3's level-4 cyclotomic
mixed-Tate ring over $\mathbb{Z}[i, 1/2]$ can realise the M2 ℤ/2
(trivial on outputs) and the M3 ℤ/2 (non-trivial, complex conjugation
on $\mathbb{Q}(\zeta_4)$) as the same ℤ/2 acting on a single enriched
object.

## TL;DR

**Verdict: OPEN-POSSIBLE.**

No clean structural Tannakian obstruction materialises at sprint-pass
depth. Three concrete candidate obstructions were considered:

1. **Weight-filtration obstruction.** M2 outputs live entirely in even
   Tate weights $\{0, +2, +4, \ldots\}$; M3 outputs live in mixed
   weights including odd weights $\{+1, +3\}$ (via $\log 2$ and
   $\zeta(3)$). One might hope the two ℤ/2 actions are forced to
   respect incompatible weight filtrations. **Does NOT force FALSIFIED.**
   A Tannakian category containing both rings already exists
   (the category $\mathrm{MT}(\mathbb{Z}[i, 1/2])$ of mixed Tate
   motives over $\mathbb{Z}[i, 1/2]$ at level 4 contains both as sub-rings),
   and the weight filtration on the mixed-Tate category is the
   *same* weight filtration that restricts trivially on the pure-Tate
   sub-category. There is no contradiction in having a single ℤ/2 act
   on the level-4 cyclotomic side and act through its quotient
   (trivially) on the pure-Tate sub-category — that is in fact the
   *standard* way pure-Tate sub-categories sit inside cyclotomic-Tate
   categories. The Galois-trivial-on-pure-Tate fact is a *consequence*
   of the categorical embedding, not an obstruction to it.

2. **Conductor-vs-signature obstruction.** M2's $\mathbb{Q}(i)$ is a
   Witt-splitting field selected by signature; M3's $\mathbb{Q}(\zeta_4)$
   is a cyclotomic conductor selected by $\chi_{-4}$. One might hope
   the two roles are forced to require structurally distinct ℤ/2
   actions (one geometric, one arithmetic). **Does NOT force
   FALSIFIED.** Both ℤ/2's are realised by *the same* underlying
   field automorphism $i \mapsto -i$. The fact that the M2 ℤ/2 acts
   trivially on the *output* ring while the M3 ℤ/2 acts non-trivially
   on the *output* ring is consistent with a single ℤ/2 acting on
   the *proof-time* level (the underlying field $\mathbb{Q}(i)$),
   with two different *fiber-functor images*. The trivial action on M2
   outputs is the action restricted to the Galois-fixed sub-ring
   $\mathbb{Q}(i)^{\mathbb{Z}/2} = \mathbb{Q}$; the non-trivial action
   on M3 outputs is the same ℤ/2 acting on a non-invariant
   sub-category. Different fiber-functor images of one Galois action
   is *exactly* the Tannakian story; it is not an obstruction.

3. **Fiber-functor compatibility obstruction.** One might hope that
   no single fiber functor $\omega: \mathcal{C}_{\mathrm{enr}} \to
   \mathrm{Vec}_{\mathbb{Q}}$ can simultaneously have M2 (pure-Tate)
   and M3 (level-4 mixed-Tate) as orthogonal sub-image sub-rings with
   the ℤ/2's identified. **Does NOT force FALSIFIED.** A single fiber
   functor on $\mathrm{MT}(\mathbb{Z}[i, 1/2])$ already has both as
   sub-categories; the restriction of the fiber functor to the
   pure-Tate sub-category and to the level-4 mixed-Tate sub-category
   is the *standard* Hodge / Betti realisation. The motivic Galois
   group acts compatibly via projection onto the relevant quotients.

The honest finding: **the published periods program already supplies
the categorical home that the Q5' question asks about.** The mixed
Tate motives over $\mathbb{Z}[i, 1/2]$ at level 4
($\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$ in Glanois 2015 notation) contain
the pure-Tate sub-category $\mathrm{MT}(\mathbb{Z})$ as a full
sub-category, and both M2 and M3 output rings embed into the larger
category compatibly. The motivic Galois group of the larger category
acts on both, and its action restricts trivially on the pure-Tate
piece because the pure-Tate piece is precisely the Galois-fixed
sub-category. This is not a unification of M2 and M3 in the deeper
sense Q5' asks for — it is rather the statement that *categorical
containment is trivial*; the substantive question (whether the
enriched motivic Galois group of the *GeoVac spectral triple
specifically* selects $\mathbb{Q}(\zeta_4)$ for both reasons
simultaneously) is *not* closed by this observation, because the
specifically-spectral-triple object has not been constructed.

The verdict is therefore **OPEN-POSSIBLE**: no obstruction at the
purely-categorical level; the obstruction (if any) lives at the
level of the spectral-triple-specific enrichment, which would
require constructing the enriched object first. This pushes back to
entry point (ii), the Phase 1 lit review, which is the legitimate
next step.

## Argument

### Setup: what the question actually asks

The question asks whether there is a Tannakian category
$\mathcal{C}_{\mathrm{enr}}$ with the following properties:

(a) $\mathcal{C}_{\mathrm{enr}}$ contains M2's output ring
    $\bigoplus_k \pi^{2k}\cdot\mathbb{Q}$ (pure-Tate over $\mathbb{Q}$)
    as the period ring of a sub-category $\mathcal{C}_2 \subset
    \mathcal{C}_{\mathrm{enr}}$.

(b) $\mathcal{C}_{\mathrm{enr}}$ contains M3's output ring (level-4
    cyclotomic mixed Tate over $\mathbb{Z}[i, 1/2]$) as the period
    ring of a sub-category $\mathcal{C}_3 \subset
    \mathcal{C}_{\mathrm{enr}}$.

(c) The motivic Galois group $G_{\mathrm{enr}} :=
    \mathrm{Aut}^\otimes(\omega_{\mathrm{enr}})$ for a fiber functor
    $\omega_{\mathrm{enr}}: \mathcal{C}_{\mathrm{enr}} \to
    \mathrm{Vec}_{\mathbb{Q}}$ contains a single ℤ/2 sub-group $H$
    such that:
    - $H$ acts via the M2 ℤ/2 (i.e., trivially on outputs) on
      $\omega(\mathcal{C}_2)$;
    - $H$ acts via the M3 ℤ/2 (i.e., as complex conjugation on
      $\mathbb{Q}(\zeta_4)$) on $\omega(\mathcal{C}_3)$.

Property (c) is the load-bearing part. (a) and (b) are weak and
trivially satisfied by taking $\mathcal{C}_{\mathrm{enr}}$ to be a
large enough ambient mixed-Tate category. The structural question is
whether the *same* ℤ/2 can do both jobs.

Sprint A7's finding can be restated in this language: A7 observed
that the two ℤ/2's act differently on outputs (trivial vs.
non-trivial). The Q5' question sharpens this to: is the "act
differently on outputs" a forced structural feature of any
Tannakian enrichment, or merely the standard feature that a single
ℤ/2 acts differently on different fiber-functor images of different
sub-categories?

### Candidate obstruction 1: weight-filtration incompatibility

**The hope:** M2's outputs live in even Tate weights; M3's outputs
include odd weights via $\log 2$ (weight $+1$) and $\zeta(3)$
(weight $+3$). The motivic weight filtration $W_\bullet$ on
$\mathcal{C}_{\mathrm{enr}}$ is preserved by any Tannakian
endomorphism. If we ask the same ℤ/2 to act trivially on the
even-weight sub-graded piece *and* non-trivially on the odd-weight
sub-graded piece, perhaps the requirement of preserving the weight
filtration produces an obstruction.

**Why it doesn't force FALSIFIED.** The weight filtration is a
*decreasing* filtration by sub-objects, and any Tannakian
endomorphism preserves it. But "acts trivially on the even-weight
piece" and "acts non-trivially on the odd-weight piece" is *exactly
the kind of thing* a Tannakian endomorphism can do — it permutes
sub-objects of the same weight, and the action on each weight slice
is independent. There is no constraint that the action on weight
$+2$ and the action on weight $+1$ be related: they live in
different graded pieces $\mathrm{gr}_2^W$ and $\mathrm{gr}_1^W$,
each of which is a separate representation of $G_{\mathrm{enr}}$.

A concrete check: in the standard $\mathrm{MT}(\mathbb{Z}[i, 1/2])$,
the motivic Galois group is an extension of the multiplicative
group $\mathbb{G}_m$ (acting on Tate twists, the pure-Tate part)
by a pro-unipotent group $U$ (acting on the mixed-Tate part). The
complex-conjugation ℤ/2 of $\mathrm{Gal}(\mathbb{Q}(i)/\mathbb{Q})$
acts as a quotient of the motivic Galois group; it acts trivially
on $\mathbb{G}_m$ (which has no non-trivial $\mathbb{Q}$-defined
ℤ/2 quotient), and it acts non-trivially on $U$ (specifically, on
the $\chi_{-4}$-isotypic component).

The decomposition "trivial on $\mathbb{G}_m$, non-trivial on
$\chi_{-4}$-isotypic" is *exactly the M2-vs-M3 split*. Far from
being an obstruction, this is the *standard* way the structure works.

**Verdict for candidate 1:** does NOT force FALSIFIED.

### Candidate obstruction 2: signature-vs-conductor mismatch

**The hope:** Sprint A7 emphasised that M2's $\mathbb{Q}(i)$ is a
Witt-splitting field (chosen by signature $(+,+,+,+)$) while M3's
$\mathbb{Q}(\zeta_4)$ is a cyclotomic conductor (chosen by
$\chi_{-4}$). The hope is that these two *roles* for $\mathbb{Q}(i)$
are structurally incompatible: a Tannakian enrichment that wants to
treat $\mathbb{Q}(i)$ as both a Witt-splitting field and a cyclotomic
conductor would have to specify *which role* the ambient field plays,
and the two roles produce different ℤ/2 actions on the underlying
field-extension category.

**Why it doesn't force FALSIFIED.** The two roles for $\mathbb{Q}(i)$
are descriptions of *where the field enters the computation*, not
descriptions of distinct mathematical objects. The Galois group
$\mathrm{Gal}(\mathbb{Q}(i)/\mathbb{Q}) = \mathbb{Z}/2$ is a single
ℤ/2 acting on a single field $\mathbb{Q}(i)$; the fact that this
ℤ/2 appears in two different proof contexts (Witt isotropy for M2's
SD integrand, cyclotomic-conductor for M3's $\eta$-invariant Mellin
slot) is a statement about *where in the GeoVac derivation* the
$\mathbb{Q}(i)$ shows up, not a statement that there are two
*different* ℤ/2's.

More precisely: in the F–M proof of M2's mixed-Tate classification
(Theorem 7.3 of FM 2016), the field $\mathbb{Q}(i)$ enters because
the standard quadric is isotropic over it; the Galois action
$i \mapsto -i$ swaps the two maximal isotropic subspaces. In
Glanois's classification of M3's level-4 cyclotomic mixed-Tate
ring (Glanois 2015, Cor. 1.1–1.2), the field $\mathbb{Q}(i) =
\mathbb{Q}(\zeta_4)$ enters because $\chi_{-4}$ is the unique
primitive character of conductor 4; the Galois action $i \mapsto -i$
swaps $\zeta_4 \leftrightarrow \zeta_4^{-1}$ and acts as the
non-trivial element of $(\mathbb{Z}/4)^\times / \{\pm 1\}$.

These are *the same field automorphism*. The Sprint A7 distinction
("two structurally distinct mechanisms") is a statement about *why*
the field $\mathbb{Q}(i)$ shows up in each case, not a statement
about *which automorphism* of the field acts. A Tannakian enrichment
that contains both contexts as sub-categories will have the same
$\mathrm{Gal}(\mathbb{Q}(i)/\mathbb{Q})$ acting on both, by
naturality. The "trivial action on M2 outputs" is the statement that
$\mathbb{Q}(i)$-Galois acts trivially on the period values of
$\mathcal{V}_n^{\mathrm{GeoVac}}$ because those period values
*happen to lie* in the Galois-fixed sub-ring $\mathbb{Q}$ (in fact in
$\bigoplus_k \pi^{2k} \cdot \mathbb{Q}$); the "non-trivial action on
M3 outputs" is the statement that $\mathbb{Q}(i)$-Galois acts
non-trivially on (e.g.) $\beta(s)$ for odd $s$ because those values
are *not* in the Galois-fixed sub-ring.

These are compatible — they are exactly what a single ℤ/2 *does* when
it acts on a category that contains both Galois-fixed and
Galois-non-fixed sub-categories.

**Verdict for candidate 2:** does NOT force FALSIFIED.

### Candidate obstruction 3: fiber-functor incompatibility

**The hope:** A Tannakian category is determined (up to equivalence)
by its fiber functor to $\mathrm{Vec}_{\mathbb{Q}}$. If no fiber
functor $\omega_{\mathrm{enr}}$ on $\mathcal{C}_{\mathrm{enr}}$ can
simultaneously identify M2's output ring and M3's output ring as
embedded sub-rings of period values with the ℤ/2 actions matching,
then the unification fails categorically.

**Why it doesn't force FALSIFIED.** The standard fiber functor on
$\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$ is the Betti realisation
(or equivalently the de Rham realisation), which produces period
values in $\mathbb{C}$. On the pure-Tate sub-category
$\mathrm{MT}(\mathbb{Z}) \subset \mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$,
this Betti realisation produces period values in
$\bigoplus_k (2\pi i)^k \cdot \mathbb{Q}$, which contains
$\bigoplus_k \pi^{2k} \cdot \mathbb{Q}$ at even Tate twists (and
becomes pure-Tate-real after fixing a real structure). On the full
level-4 cyclotomic sub-category, the same fiber functor produces
the M3 output ring.

A single fiber functor handles both. The motivic Galois group
$\mathrm{Aut}^\otimes(\omega_{\mathrm{Betti}})$ acts on both
sub-rings compatibly, restricting to the unipotent piece on
$\mathrm{MT}(\mathbb{Z})$ (which has no non-trivial ℤ/2 quotient
acting on the $\mathbb{Q}$-rational pure-Tate periods, as expected:
the pure-Tate periods are Galois-fixed) and acting non-trivially on
the level-4 cyclotomic piece.

This is in fact the *standard construction*; it does not require
any new enrichment of the spectral-triple-specific kind that Q5'
contemplates.

**Verdict for candidate 3:** does NOT force FALSIFIED.

### Why "OPEN-POSSIBLE" and not "FALSIFIED-CHEAPLY"

The three candidate obstructions above all fail to force FALSIFIED,
*and* they each provide a positive direction: the standard ambient
category $\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$ already supplies the
categorical home for both M2 and M3 outputs, with a single fiber
functor and a single motivic Galois group whose action restricts
correctly to each sub-ring. **The naive form of the Q5'
question—"is there a Tannakian category containing both?"—has answer
YES, with $\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$.**

But this does NOT close Q5' in the deep sense, because the deep Q5'
question is not "is there a Tannakian category containing both?"
(trivially yes, by ambient containment) but **"is there a Tannakian
enrichment of the discrete GeoVac spectral triple itself whose
motivic Galois group is the M2/M3 unification?"** This sharper
question requires constructing the spectral-triple-side motivic
Galois group as a *specific quotient or extension* of
$\mathrm{Aut}^\otimes(\omega_{\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)})$,
or as a *different group* mapping to it. Neither construction
exists in the published literature.

### Why "OPEN-POSSIBLE" and not "BORDERLINE"

A BORDERLINE verdict would require a partial obstruction — say, an
obstruction in one direction (M2 → M3 enrichment fails) but not the
other (M3 → M2 enrichment succeeds). None of the three candidates
produces such asymmetry: the categorical embedding $M_2 \hookrightarrow
M_3$ already works at the ambient mixed-Tate level, and there is no
known structural reason the spectral-triple-specific enrichment must
fail in one direction.

### What WOULD force FALSIFIED, if it existed

A clean FALSIFIED would require one of the following, none of which
materialised at sprint-pass depth:

(i) **Weight-filtration cohomological obstruction.** Show that any
    fiber functor of $\mathcal{C}_{\mathrm{enr}}$ restricting to
    Betti-Hodge on M2 and to motivic-iterated-integrals on M3 has a
    non-trivial Ext-class obstruction (perhaps an
    $\mathrm{Ext}^1_{G_{\mathrm{enr}}}(\mathbb{Q}(\chi_{-4}),
    \mathbb{Q}(0)) \ne 0$ statement). I see no such cohomological
    incompatibility in the published periods literature; the standard
    embeddings $\mathrm{MT}(\mathbb{Z}) \subset
    \mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$ are exact.

(ii) **GeoVac-specific weight obstruction.** Show that the
     GeoVac-specific weight filtration on the spectral triple
     (Camporesi-Higuchi half-integer shift, vertex-parity decomp)
     produces a different weight-grading than the standard motivic
     weight, and the two cannot be reconciled. The GeoVac empirical
     witness panel (Paper 55 §5.5 + Paper 28 Thm 3) shows
     vertex-parity-as-$\chi_{-4}$ matches the Deligne–Glanois
     descent *bit-exactly*, so the two weight gradings agree where
     they overlap; no obstruction here.

(iii) **No-go theorem in Glanois 2015 or follow-ups.** Glanois
      proves that level $N = 4$ cyclotomic mixed-Tate cannot be
      reduced to pure-Tate level $N = 1$ except by descent on a
      symmetry subgroup. This is a *structural* statement about the
      level filtration, not a Tannakian obstruction. The pure-Tate
      sub-category $\mathrm{MT}(\mathbb{Z}, 1)$ sits inside
      $\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$ as a *closed* full
      sub-category; the descent statement is precisely the
      compatibility of the ambient ℤ/2 with the embedding, not its
      failure.

None of (i)–(iii) is supplied by the existing periods literature.

### Cross-check: the deeper reading of Sprint A7

Sprint A7 wrote: "the two ℤ/2's are not the same ℤ/2 in the sense of
a single Tannakian symmetry of the underlying spectral triple". The
present probe sharpens what "in the sense of" means here:

- *At the categorical level (ambient mixed-Tate category):* the two
  ℤ/2's ARE the same ℤ/2 — namely the single
  $\mathrm{Gal}(\mathbb{Q}(i)/\mathbb{Q})$ acting on
  $\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$. This is what the present probe
  finds.

- *At the spectral-triple-specific enrichment level:* the two ℤ/2's
  are not (yet) known to be the same ℤ/2 — because the
  spectral-triple-specific enrichment has not been constructed. This
  is the open Sprint A7 question.

The two statements are compatible; Sprint A7 was making the second
statement, the present probe finds no obstruction at the first level
and re-confirms the open question is at the second level.

## Verdict

**OPEN-POSSIBLE.**

No clean structural Tannakian obstruction materialises at
sprint-pass depth. The three concrete candidate obstructions
considered (weight-filtration incompatibility; signature-vs-conductor
mismatch; fiber-functor incompatibility) all FAIL to force
FALSIFIED. The standard ambient category
$\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$ supplies a categorical home in
which both M2 and M3 outputs sit as embedded sub-rings with the
same $\mathrm{Gal}(\mathbb{Q}(i)/\mathbb{Q})$ acting compatibly.
The trivial-on-M2 / non-trivial-on-M3 action pattern is the
*expected* pattern of a single Galois acting on a category with
both Galois-fixed and Galois-non-fixed sub-categories — not an
obstruction.

The deeper Q5' question — whether the
**spectral-triple-specific** enrichment unifies the two — is *not*
closed by this probe, because the spectral-triple-specific enriched
motivic Galois group has not been constructed. The entry point (ii)
Phase 1 lit review remains justified as the legitimate next step.

## Concrete candidates for what would close Q5'

If a future sprint or external collaboration wants to push Q5' to
FALSIFIED or PROVEN, the following concrete targets are
sprint-scale:

(α) **Construct the GeoVac-specific motivic Galois group.** Define
    $G^{\mathrm{GeoVac}}_{\mathrm{mot}}$ as $\mathrm{Aut}^\otimes$ of a
    fiber functor on the Tannakian category generated by the
    Camporesi-Higuchi Dirac Dirichlet series and the SD heat-trace
    Mellin transforms. Compute its image in
    $G^{\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)}_{\mathrm{mot}}$.

(β) **Identify the natural ℤ/2 inside the GeoVac motivic Galois
    group.** Test if it factorises as the
    $\mathrm{Gal}(\mathbb{Q}(i)/\mathbb{Q})$ already discussed, or
    as a different ℤ/2 unique to the spectral-triple side. If the
    GeoVac-specific ℤ/2 is identified with $\mathrm{Gal}(\mathbb{Q}(i)/
    \mathbb{Q})$ acting on both M2 and M3, Q5' resolves PROVEN.

(γ) **Find an obstruction at the spectral-triple level.** Construct
    a specific period of the spectral triple whose
    motivic-Galois-orbit is forced to be larger than what the
    ambient $\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$ Galois orbit
    predicts. This would FALSIFY Q5'.

The present sprint does not pursue any of (α)–(γ); they are flagged
for downstream sprints if Q5' continues to interest.

## Honest scope

1. **No fitted parameters introduced**: this is pure structural
   reasoning.
2. **No paper edits applied**: this memo records a verdict; whether
   to update Paper 55 §subsec:open_m2_m3 Q5' wording is a PI call.
3. **Paper 2 framing preserved**: no edits to combination rule
   $K = \pi(B + F - \Delta)$ "conjectural" label.
4. **No PSLQ, no implementation, no code**: as scoped.
5. **Three concrete candidate obstructions considered** (above);
   each FAILS to force FALSIFIED. The verdict OPEN-POSSIBLE is the
   honest finding, not a punted verdict.
6. **The HALF-STRUCTURAL verdict of Sprint A7 is preserved.** Sprint
   A7 said "two structurally distinct mechanisms produce the same
   field at the number-field level, but no shared Tannakian symmetry
   is *known*." The present probe sharpens "known" to "the ambient
   mixed-Tate category gives a categorical home; the
   spectral-triple-specific enrichment is open." A7's verdict and
   the present verdict are compatible refinements of the same
   underlying state of knowledge.
7. **Phase 1 lit review (entry point (ii)) is now justified.** The
   present probe was the cheap-out attempt; it did not falsify, so
   the more expensive structural search is on the menu. This memo
   does not recommend launching it; that is a PI call after
   reviewing this verdict.
8. **Sprint-pass depth caveat.** This probe was 1-day sprint-scale.
   A deeper structural argument might exist that this probe did not
   surface. The OPEN-POSSIBLE verdict means "no obstruction *found
   at this depth*", not "no obstruction exists." A genuine no-go
   could surface in entry point (ii)'s lit review or in a multi-month
   construction attempt at (α)–(γ).

## Files used

### Memos read
- `debug/sprint_a7_m2_m3_cyclotomic_memo.md` (Sprint A7 base verdict
  HALF-STRUCTURAL; structural table at line 140; alternative readings
  (a)/(b)/(c) at lines 274–287).

### Papers read
- `papers/group3_foundations/paper_55_periods_of_geovac.tex`
  §subsec:open_m2_m3 (Q5' wording, lines 1474–1499); Remark 4.5
  `rem:m2_specialisation` (lines 534–605, including the
  distinction-from-M3 paragraph already added); Proposition 5.5
  `prop:vertex_parity_descent` (lines 905–924); §6.2 Tate-weight
  bookkeeping (lines 1067–1098); §6.4 Paper 50 F-theorem joint
  example with weight matching (lines 1140–1167).

### Published references (cited from Paper 55 bibliography)
- Deligne arXiv:math/0302267 (period map isomorphism at
  $N \in \{2, 3, 4, 6, 8\}$).
- Glanois arXiv:1411.4947 (level-4 cyclotomic mixed-Tate Hopf algebra
  structure; basis $\mathcal{B}^4$ with $\epsilon_p \in \{\pm 1, \pm i\}$).
- Brown 2012 (level-1 MZV ring baseline; the embedding
  $\mathrm{MT}(\mathbb{Z}, 1) \subset \mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$
  is the standard inclusion).
- Fathizadeh–Marcolli arXiv:1611.01815 Theorem 7.3 (Witt isotropy over
  $\mathbb{Q}(i)$ for M2 case).

### Scripts (none required)
The probe is purely structural and rests on the published categorical
embedding $\mathrm{MT}(\mathbb{Z}) \subset
\mathrm{MT}(\mathbb{Z}[i, 1/2], 4)$ and standard properties of
Tannakian fiber functors. No numerical computation required.
