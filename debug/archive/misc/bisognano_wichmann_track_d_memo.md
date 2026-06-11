# Track D — Bisognano–Wichmann reading of Sprint TD Track 4

**Date:** 2026-05-09
**Track:** Track D (Candidate A) of the post-multifocal-catalogue parallel sprint
**Builds on:** Track B Lorentz-boost scoping memo (`debug/lorentz_boost_scoping_memo.md`); Sprint TD Track 1 (`debug/sprint_td_track1_memo.md`); Sprint TD Track 4 (`debug/sprint_td_track4_memo.md`); Sprint MR-A/B/C master Mellin engine domain partition (`debug/mr_b_spectral_action_rate_memo.md`); Paper 38 GH-convergence theorem (WH1 PROVEN).
**Sprint window:** 1 week.
**Status:** **Closed positive**, structural-correspondence verdict (not literal identification).
**No production code modified.**

---

## Executive summary

The 2π factor entering Sprint TD Track 4's Hawking temperature
$T_H = 1/(8\pi M)$ on the Euclidean Schwarzschild cigar — identified there as
the M1 Hopf-base measure / $S^1_\tau$ circumference signature of the master
Mellin engine — is *the same numerical 2π* that appears as the
modular-automorphism period of the wedge / horizon vacuum in the
Bisognano–Wichmann theorem (1976) and its Schwarzschild generalization
(Sewell 1982). The standard Wick-rotation chain (Hartle–Hawking 1976
imaginary-time periodicity ↔ Bisognano–Wichmann modular flow ↔ Sewell
1982 KMS / Killing-time isometry) identifies the framework-side $2\pi$
with the Lorentzian-side boost-orbit period at the level of
correlation-function periodicity. **The framework's M1 mechanism IS the
Wick-rotated image of the Lorentz-boost / modular-flow mechanism on the
Schwarzschild horizon.**

This is a **structural correspondence**, not a literal identification.
Two distinct mathematical objects (a Hopf-bundle measure factor on the
Riemannian side; a modular-automorphism period on the Lorentzian side)
are equated by a published Wick-rotation prescription, not by an
operator-level proof inside the GeoVac spectral triple. The literal
identification — that the modular Hamiltonian on the truncated metric
spectral triple $\mathcal{T}_{n_{\max}}$ (Paper 38 framework) closes at
period $2\pi$ — is named below as the principal falsifier and
recommended as the natural 4–8 week follow-up.

**Top 3 deliverables of this 1-week sprint:**

1. **Memo (this document, ~7000 words)** documenting the
   Bisognano–Wichmann reading of Sprint TD Track 4, with verified
   literature citations and a five-section structural analysis.
2. **Three paper-update paragraphs** applied directly:
   Paper 32 §VIII master-Mellin-engine-domain remark gains a Bisognano–Wichmann
   pendant; Paper 35 §VIII gains a new Bisognano–Wichmann subsection;
   Paper 34 §VIII Lorentz-boost open question gains a paragraph noting
   that Track D closes the documentation loop.
3. **One off-precision catalogue row** added to Paper 34 §V.B for Hawking
   T (machinery-witness, error class C); one footnote added to Paper 34
   §III.15 (Observation/temporal-window projection) naming the
   Bisognano–Wichmann reading.

**Two corrections to Track B scoping memo data:** Sewell publication
year is 1982 (not 1980 as the scoping data and memo body had); Bisognano–Wichmann
1976 paper title is *On the duality condition for quantum fields* (not
*for a Hermitian scalar field*, which is the title of the *1975* paper).
Verified via WebSearch on AIP/ADS records; corrections applied to
Track D bibliography.

---

## §1. The structural identification — what is claimed, and what kind of claim it is

### §1.1 The two 2π's

The framework side. In Sprint TD Track 4 (`debug/sprint_td_track4_memo.md`,
`debug/data/sprint_td_track4.json`), the Euclidean Schwarzschild metric
near the horizon $r = 2M$ takes the local cigar form
$ds^2 \approx \rho^2\, d\phi^2 + d\rho^2 + (2M)^2 d\Omega^2$ with
$\phi = \tau / (4M)$. Smoothness at $\rho = 0$ — the cigar tip — requires
$\phi$ to have period $2\pi$ (any other period gives a conical defect, as
the standard Hawking–Gibbons argument). The $\tau$-circle therefore has
period
$$
\beta_{\rm cigar} \;=\; 4M \cdot 2\pi \;=\; 8\pi M,
\qquad T_H \;=\; \frac{1}{8\pi M},
\qquad T_H \;=\; \frac{\kappa_g}{2\pi},
$$
with surface gravity $\kappa_g = 1/(4M)$. The single $2\pi$ in $T_H$ sits
in the denominator and is *exactly* the circumference of the $\phi$-circle
in the cigar tip. In the master-Mellin-engine-domain partition
(Sprint MR-A/B/C, May 2026; Paper 32 §VIII Remark
`rem:master_mellin_domain`), this $2\pi$ is the **M1 Hopf-base measure
signature**: M1 corresponds to $\mathcal{M}[\mathrm{Tr}(D^0 \cdot e^{-tD^2})]$
which produces a $\mathrm{Vol} / \mathrm{Mellin\text{-}measure}$ factor;
on a temporal $S^1_\beta$ alone, this gives $\mathrm{Vol}(S^1_\beta) = \beta$
directly, with the bare $2\pi$ coming from the unit-circumference
normalization. The same M1 mechanism on the spatial $S^3$ Hopf base gives
$\mathrm{Vol}(S^2)/\pi^2 = 4/\pi$ (the L2 propinquity rate constant from
Sprint MR-B's reading; see Paper 38 §VIII).

The Lorentzian side. The Bisognano–Wichmann theorem (1976) — proven for
general Wightman-axiomatic QFT, including bosons and fermions — states
that the modular automorphism $\sigma_t$ of the algebra of observables
localized in a Rindler wedge, with respect to the vacuum state, *is* the
unitary representation of the Lorentz boost subgroup that maps the wedge
to itself. Concretely, $\sigma_t = e^{itK}$ where $K$ is the modular
Hamiltonian, and $K = -2\pi \cdot M$ where $M$ is the boost generator
(rapidity). The KMS condition (Tomita–Takesaki) gives
$\sigma_{i \cdot 2\pi} = $ identity: the modular-flow analytic period
along the imaginary axis is exactly $2\pi$.

For Schwarzschild, Sewell 1982 generalized this: the natural vacuum on
the Kruskal-extended Schwarzschild spacetime (the Hartle–Hawking vacuum
of Hartle–Hawking 1976), restricted to the static exterior region, is a
KMS state at inverse temperature $\beta = 8\pi M$ with respect to the
Killing time evolution, and the corresponding Tomita modular automorphism
is the surface-gravity-rescaled Killing flow. The same theorem applied to
a Rindler wedge gives the Unruh effect, with $T_U = a / (2\pi)$ for a
uniformly accelerated observer at proper acceleration $a$.

The numerical alignment. Both 2π's are the same number. The Hartle–Hawking
1976 path-integral derivation showed that the Wick-rotated Schwarzschild
geometry has *imaginary-time period* $8\pi M$, and that QFT correlators
on this Euclidean geometry are thermal at temperature $T_H = 1/(8\pi M)$
under the standard Wick rotation $\tau \to it$. Sewell 1982 then proved
that this thermal state IS a modular-flow image of a Killing-time
isometry. **The Wick-rotation chain
[Hartle–Hawking 1976] → [Sewell 1982] → [Bisognano–Wichmann 1976]
identifies the cigar's $\tau$-circle period $8\pi M$ with the Lorentzian
modular-flow / boost-orbit period.** The $2\pi$ inside $\beta = 8\pi M =
4M \cdot 2\pi$ is the $2\pi$ of the modular-flow analytic period.

### §1.2 The structural identification at the framework level

Sprint TD Track 4 reproduces $\beta = 8\pi M$ from the M1 mechanism on
the cigar's $S^1_\tau$ factor (sympy residual zero; `debug/data/sprint_td_track4.json`
`1_cigar_beta_symbolic` and `4_master_mellin_verdict.V1_M1_mechanism`).
The published Wick-rotation chain (Hartle–Hawking 1976, Sewell 1982,
Bisognano–Wichmann 1976) identifies this $\beta$ with the
modular-flow / Lorentz-boost orbit period.

**Therefore: the framework-internal M1 mechanism, instantiated on the
$S^1_\tau$ factor of the Schwarzschild cigar, IS the Wick-rotated image
of the Bisognano–Wichmann / Sewell modular-flow / Lorentz-boost mechanism.**

This is the headline claim of Track D. In the case-exhaustion-theorem
language of Paper 32 §VIII (Sprint TS-E1, May 2026), the M1 sub-mechanism
of the master Mellin engine has *two physical readings* on temporal
compactifications:

| Reading | Formal object | Physical content |
|:---|:---|:---|
| Riemannian | $\mathcal{M}[\mathrm{Tr}(D^0\, e^{-tD^2})]$ | Hopf-base measure / Vol(S¹) circumference factor |
| Lorentzian (via Wick rotation) | Bisognano–Wichmann modular flow | Lorentz-boost orbit period; Killing-time isometry |

These are not two different mechanisms; they are two different
*vocabularies* for the same M1 sub-mechanism.

### §1.3 Literal identification vs structural correspondence

Three readings of "the two $2\pi$'s are the same":

1. **Literal identification.** The $2\pi$ in M1 and the $2\pi$ in
   Bisognano–Wichmann are the *same analytic object* — i.e., one
   mathematical structure with two presentations, identified by an
   equality of operators on the spectral triple's Hilbert space. To
   *prove* this at the framework level would require constructing the
   modular Hamiltonian on the truncated metric spectral triple
   $\mathcal{T}_{n_{\max}}$ (Paper 38 framework) for a wedge-restricted
   state, and verifying directly that its imaginary-axis period is
   $2\pi$. This is not done in Track D and is the principal falsifier
   below (§3.1).
2. **Structural correspondence.** The two $2\pi$'s appear in two
   different mathematical objects (a Hopf-bundle measure factor on the
   framework side; a modular-automorphism period on the Lorentzian side)
   which standard Wick-rotation arguments (Sewell 1982; Hartle–Hawking
   1976; Bisognano–Wichmann 1976) identify at the level of
   correlation-function periodicity.
3. **Numerical coincidence.** The two $2\pi$'s happen to take the same
   numerical value but live in unrelated mechanisms. This reading is
   refuted by the Wick-rotation chain: Hartle–Hawking 1976 proved
   directly that the imaginary-time period of the Lorentzian thermal
   state IS $\beta = 8\pi M$, the same $\beta$ that the M1 mechanism
   reproduces on the cigar. The two $2\pi$'s are not coincident; they
   are the same $2\pi$ via a published derivation chain.

**Track D's verdict: structural correspondence, not literal
identification.** The Wick-rotation chain elevates the M1 reading to
a Lorentz-boost reading; what is missing for literal identification is
an operator-level proof of $\sigma_{i \cdot 2\pi} = $ identity inside
the framework's truncated spectral triple.

This is the honest scope of Track D as a 1-week documentation pass: name
the structural correspondence in the standard Bisognano–Wichmann
vocabulary, identify what would lift it to literal identification (§3.1),
and apply the appropriate paper-update paragraphs without
overclaiming.

### §1.4 What was already in the framework

It is important to be precise about what Track D *adds* versus what
Sprint TD Track 4 already contained. Sprint TD Track 4 already showed
that the M1 mechanism reproduces $T_H = 1/(8\pi M)$ on the Euclidean
cigar with sympy residual zero. The cigar is the Wick rotation of the
Lorentzian Schwarzschild horizon; this is textbook physics. What Track D
contributes is the *naming* of the Lorentz-boost / Bisognano–Wichmann
interpretation at the framework level — placing Sprint TD Track 4 inside
the published QFT machinery that interprets the Wick-rotated thermal
result as a modular-flow / boost-orbit-period statement. The new content
is the **framework-internal identification of M1 with the boost orbit
period**, embedded in the case-exhaustion theorem of Paper 32 §VIII.

The interpretation itself is published (Sewell 1982 says exactly this on
the Lorentzian side); what is new is the placement inside the master
Mellin engine partition that Paper 32 §VIII / Paper 18 §III.7 / Sprint
MR-A/B/C established. Saying it differently: Bisognano–Wichmann +
Sewell 1982 is the bridge from Lorentzian-side modular flow to the
Wick-rotated thermal partition function; Sprint TD Track 4 + Paper 32
§VIII case-exhaustion is the bridge from the Wick-rotated thermal
partition function to the framework's master Mellin engine. Track D is
the explicit composition of these two bridges.

---

## §2. Operator-level construction — what we have and what's missing

### §2.1 What we have (Sprint TD Track 1 spectrum-level construction)

Sprint TD Track 1 built the explicit tensor-product spectral triple
$\mathcal{T}_{S^3} \otimes \mathcal{T}_{S^1_\beta}$ at the spectrum level
(`geovac/thermal_tensor_triple.py`, ~580 lines, 36 tests). The Dirac
operator is
$$
D_{\rm total} \;=\; D_{S^3} \otimes 1 + \gamma_{S^3} \otimes D_\tau,
$$
with the anti-commutator
$\{D_{S^3} \otimes 1,\, \gamma_{S^3} \otimes D_\tau\} = 0$ making
$D_{\rm total}^2$ block-diagonal so spectra add. The KMS state on the
tensor-product algebra at inverse temperature $\beta = 8\pi M$ is the
standard Gibbs state $\rho = Z^{-1} e^{-\beta H}$. Track 1 verified this
for general $\beta$ (Stefan–Boltzmann, partition function, modular
residual). Track 4 plugged $\beta = 8\pi M$ into Track 1's
`matsubara_spectrum()` and obtained the standard Hawking spectrum
verbatim — bosonic lowest nonzero $\omega_1 = 2\pi/\beta = 1/(4M) = \kappa_g$
(surface gravity), fermionic lowest $\omega_0 = \pi/\beta = 1/(8M) =
\pi T_H$.

The framework's existing infrastructure includes:
- Truncated metric spectral triple at finite $n_{\max}$ (Paper 32 §III;
  `geovac/operator_system.py`).
- GH-convergence theorem $\mathcal{T}_{n_{\max}} \to \mathcal{T}_{S^3}$
  in Latrémolière propinquity (Paper 38, WH1 PROVEN).
- Connes distance SDP computation on $\mathcal{T}_{n_{\max}}$
  (`geovac/connes_distance.py`; Paper 32 §III rem:connes_distance).
- Real structure $J$ verified at finite $n_{\max}$ on truthful
  Camporesi–Higuchi (Paper 32 §IV).

### §2.2 What we don't have — the operator-level Bisognano–Wichmann statement

The literal identification of Track D's structural correspondence would
require:

1. **Construct the modular Hamiltonian $K$ on
   $\mathcal{T}_{n_{\max}}$** for a wedge-restricted state. The natural
   GeoVac-side analog of a Rindler wedge is a half-$S^3$ subalgebra (e.g.,
   the algebra of observables localized in one hemisphere of $S^3$, which
   is not a Lorentzian wedge but is the simplest non-trivial Tomita
   module). Tomita–Takesaki on truncated operator systems is studied in
   the Connes–van Suijlekom 2021 spectral-truncation framework but the
   modular flow has been explicitly constructed only in flat / abelian
   cases (Connes–vS 2021; Hekkelman–McDonald 2024); the SU(2) /
   Camporesi–Higuchi case is not done.

2. **Verify $\sigma_{i \cdot 2\pi} = $ identity** at finite $n_{\max}$.
   This is the Tomita–Takesaki analytic-continuation statement; if the
   modular flow is a one-parameter group on the operator system $O_{n_{\max}}$,
   its analytic continuation to imaginary parameter $t \to i\theta$ should
   close at $\theta = 2\pi$. Computing this requires either (a) the SDP
   framework of `geovac/connes_distance.py` extended to modular flow, or
   (b) a direct algebraic computation of the polar decomposition of $S$
   = closure of $a\Omega \to a^* \Omega$ on the GNS Hilbert space.

3. **Take the GH limit using Paper 38's L1'–L5 lemmas**, verify that
   the modular-flow period stays at $2\pi$ in the limit. Paper 38's
   propinquity machinery is for the metric spectral triple; extending it
   to the modular automorphism would require a "modular propinquity" —
   not published, but plausibly derivable from L4 (Berezin
   reconstruction) and the cb-norm bound that L2 produces.

4. **Identify the limiting period $2\pi$ with the M1 Hopf-base-measure
   $2\pi$** via the case-exhaustion theorem (Paper 32 §VIII).

Estimated cost: **4–8 weeks** of focused operator-system / Connes–vS
modular-flow work. This is *not* Track B's "multi-month full Lorentzian
extension" (which requires a Lorentzian propinquity, 6–12 months
original NCG-math) — it is the much narrower problem of verifying the
modular-flow period at finite $n_{\max}$ on the *Riemannian* truncated
spectral triple, using the published modular-truncation framework as
the starting point. **It is the single highest-leverage 4–8 week
follow-up to Track D**, and it would lift the structural correspondence
to literal identification.

### §2.3 What Track D therefore commits to

Track D's operator-level claim is restricted to:

- **The M1 mechanism produces $\beta = 8\pi M$ on the cigar.**
  Sympy-verified in Sprint TD Track 4.
- **This $\beta$ is the imaginary-time period of the Hartle–Hawking
  thermal state.** Standard Wick-rotation, published in Hartle–Hawking
  1976.
- **Sewell 1982 proves this thermal state IS a Bisognano–Wichmann
  modular flow.** Standard generalization of Bisognano–Wichmann 1976.
- **The framework's M1-side $2\pi$ and the Lorentzian-side modular-flow
  $2\pi$ are the same number, by the published Wick-rotation chain.**

We do not commit to operator-level closure of $\sigma_{i \cdot 2\pi} = $
identity inside the framework. The structural correspondence stands on
the published QFT machinery.

---

## §3. Spectral-triple-level falsifiers

Three falsifiers are named, each with operationalization, expected
outcome, and priority. The falsifiers test progressively stronger
versions of the Bisognano–Wichmann reading.

### §3.1 Falsifier #1 — Modular automorphism on the truncated triple fails to close at $2\pi$

**Claim being falsified:** The structural correspondence of §1 lifts to
literal identification at finite $n_{\max}$. That is, the modular
Hamiltonian $K$ on the truncated operator system $O_{n_{\max}}$ for a
wedge-restricted Camporesi–Higuchi state has imaginary-axis period
$2\pi$.

**Operationalization:** Extend the existing Connes-distance SDP framework
(`geovac/connes_distance.py`, Paper 32 §III rem:connes_distance) to
compute the modular automorphism on $O_{n_{\max}}$. Steps:

1. Choose a wedge-like sub-algebra of $\mathcal{A}_{\rm GV}$ at finite
   $n_{\max}$. Natural candidate: the algebra of multiplication operators
   by $f \in C^\infty(S^3)$ supported on one hemisphere (the
   Riemannian-side analog of a Rindler wedge; well-defined by SO(3) ⊂ SO(4)
   reduction to the rotation subgroup that fixes a single point on
   $S^2$ Hopf base).
2. Compute the GNS construction for the half-$S^3$ algebra acting on
   $\mathcal{H}_{n_{\max}}$ from the truncated CH spectral triple state.
3. Construct the polar decomposition of $S: a\Omega \to a^*\Omega$.
4. The modular operator $\Delta = S^*S$ generates the modular flow
   $\sigma_t = \Delta^{it}$.
5. Verify whether $\Delta^{2\pi i}$ is the identity at $n_{\max} = 2, 3, 4$.

**Expected outcome:** Pass at finite $n_{\max}$. The Bisognano–Wichmann
reading is published QFT for general Wightman-axiomatic theories; the
framework's M1 reproduces the right $\beta$; failure here would be
surprising and would indicate the Camporesi–Higuchi spectrum has a
different KMS structure than the standard Wightman-axiomatic vacuum.

**If it fails:** The Bisognano–Wichmann reading of Sprint TD Track 4 is
wrong; the M1 mechanism's $2\pi$ and Bisognano–Wichmann's $2\pi$ are
coincidentally equal at the numerical level but structurally unrelated.
This would be a major structural finding, indicating the framework's
KMS structure on the truncated triple is not compatible with the
Wightman-axiomatic vacuum's KMS structure. The follow-up sprint would
be a rederivation of the modular flow without assuming Wightman axioms.

**Priority:** **MEDIUM-HIGH** for a 4–8 week follow-up. Lifts the
structural correspondence to literal identification on the truncated
spectral triple — at finite $n_{\max}$, in the same framework where
WH1 PROVEN was established. Would make Sprint TD Track 4 a
*framework-internal* Bisognano–Wichmann statement at finite $n_{\max}$.

### §3.2 Falsifier #2 — The $2\pi$ in M1 has a different mechanism than the boost orbit period

**Claim being falsified:** The framework-side $2\pi$ in $\beta = 8\pi M$
comes from the master Mellin engine M1 sub-mechanism; the
Lorentzian-side $2\pi$ in Bisognano–Wichmann comes from the
modular-flow period; the two are identified by the Wick-rotation chain.

**Operationalization:** Trace the $2\pi$ inside Sprint TD Track 4's M1
mechanism step-by-step. The cigar metric near the tip is $ds^2 \approx
\rho^2 d\phi^2 + d\rho^2 + \ldots$; smoothness requires $\phi \in [0, 2\pi)$
to avoid conical defect. This $2\pi$ is the *circumference of the unit
$S^1$* in the Hopf-bundle measure, NOT the imaginary-time period of a
modular flow on the framework's spectral triple.

Compare with the L2 propinquity rate $4/\pi = \mathrm{Vol}(S^2)/\pi^2$
on $S^3$ (Paper 38 / Sprint MR-B). Here the Hopf-base measure is the
spatial $S^2$, not a temporal $S^1$. The $\pi^2$ in the denominator of
$4/\pi$ comes from $\mathrm{Vol}(S^2) = 4\pi$ + $\mathrm{Vol}(S^3)/\pi^2 =
2$.

The two M1 instantiations (cigar $S^1_\tau$ versus $S^3$ Hopf base) are
on different sub-manifolds, but both produce M1 transcendental ring
content (rational multiples of $\pi$). **They are different
*instantiations* of the same M1 sub-mechanism, not different
mechanisms.** The structural correspondence holds because the cigar's
M1 instance happens on the Wick-rotated horizon's $S^1_\tau$, which IS
the Lorentzian-side modular-flow circle.

**If a careful trace shows the $2\pi$ in M1 comes from a non-temporal
factor:** The Bisognano–Wichmann reading would not pass through the
$2\pi$ that Sprint TD Track 4 actually produces, and the structural
correspondence would fail.

**Verdict:** Sprint TD Track 4 explicitly identifies the cigar $\phi$-period
$2\pi$ as the source of $\beta = 8\pi M$. The Wick-rotation map from
$S^1_\tau$ Matsubara modes to Lorentz-boost orbit parameters is
structurally clean (Sewell 1982). Falsifier #2 is *not currently active*
— the $2\pi$ in M1 IS the temporal circle, by Sprint TD Track 4's
construction. We name this falsifier to make the structural
correspondence falsifiable in principle.

**Priority:** **CLOSED-LOW**. The mechanism trace is already in Sprint
TD Track 4 / Track 1; falsifier #2 is essentially a sanity check, not
an open question. Listed for completeness.

### §3.3 Falsifier #3 — A Krein-space / Lorentzian-NCG calculation produces a different $2\pi$

**Claim being falsified:** A direct Lorentzian-side calculation in the
Bizi–Brouder–Besnard $(m,n)$ framework at signature $(3,1)$ (cf. Track B
scoping memo §2.1; arXiv:1611.07062) produces a modular-flow period
that agrees with the framework's M1 $2\pi$.

**Operationalization:** Track B Candidate C — Krein-lift Connes-axiom
audit at $(m,n) = (3,1)$. New module
`debug/krein_lift_audit.py` (Track B scoping memo §5.3) extending the
Paper 32 §IV finite-$n_{\max}$ Connes-axiom audit to the Krein-space lift:

1. Lift the Camporesi–Higuchi $J$ on $S^3$ to the Krein-space charge
   conjugation on $S^3 \times \mathbb{R}_t$.
2. Verify the Bizi–Brouder–Besnard $(3,1)$ sign table.
3. Compute the modular Hamiltonian on the $(3,1)$ Krein lift for the
   Wightman-axiomatic vacuum state.
4. Verify the $2\pi$ period of the modular flow on the Krein lift, and
   compare to the framework's M1 $2\pi$.

**Expected outcome:** Pass. The Bizi–Brouder–Besnard sign table at
$(3,1)$ is well-developed; modular flow period $2\pi$ follows from the
standard KMS-Tomita-Takesaki algebra (Connes 1994 NCG book). The
framework-internal Camporesi–Higuchi $J$ audit at $n_{\max} = 1, 2, 3$
passes (Paper 32 §IV); the $(3,1)$ Krein lift should preserve this with
the modified sign structure.

**If it fails:** A discrepancy with Bizi–Brouder–Besnard would indicate
GeoVac's $J$ does not fit the standard Lorentzian NCG framework. This
would be a surprising negative — the literature is well-developed for
$(3,1)$ — and would refute the structural correspondence by exhibiting
a Lorentzian-side calculation that produces a *different* $2\pi$.

**Priority:** **LOW**. Track B already scoped this as Candidate C with
LOW priority (1.5 weeks of work, useful as a stress test for the
Bizi–Brouder–Besnard prescription against GeoVac at finite $n_{\max}$,
not directly load-bearing for the structural correspondence). Defer
until a multi-month Lorentzian extension is contemplated.

### §3.4 Summary table

| Falsifier | What it tests | Priority | Window |
|:---|:---|:---|:---|
| #1: $\sigma_{i \cdot 2\pi} = $ identity on $O_{n_{\max}}$ | Lifts correspondence to literal identification at finite $n_{\max}$ | MEDIUM-HIGH | 4-8 weeks |
| #2: Trace $2\pi$ source in M1 | Sanity check that M1 $2\pi$ IS on $S^1_\tau$ | CLOSED-LOW | done in Sprint TD Track 4 |
| #3: Bizi–Brouder–Besnard $(3,1)$ Krein audit | Sanity check from Lorentzian NCG side | LOW | 1.5 weeks |

---

## §4. Position vs Bisognano–Wichmann literature

### §4.1 What is already published

The Bisognano–Wichmann reading is published QFT, with a mature
literature spanning AQFT, NCG, gravitational thermodynamics, and (more
recently) holographic entanglement.

**Foundational layer (Wightman-axiomatic AQFT).**
Bisognano–Wichmann 1975 (J. Math. Phys. 16, 985) introduced the duality
condition for Hermitian scalar fields. Bisognano–Wichmann 1976 (J. Math.
Phys. 17, 303) extended to general bosonic and fermionic fields and
proved the modular-flow ↔ Lorentz-boost identification. (Note the title
correction: the *1976* paper is "On the duality condition for quantum
fields", general fields; the *1975* paper is "On the duality condition
for a Hermitian scalar field", scalar fields. The Track B scoping memo
data file had this swapped.)

**Gravitational generalization.** Sewell 1982 (Annals of Physics 141,
201–224 — note 1982, not 1980 as the scoping memo data and body had)
generalized to globally hyperbolic Lorentzian manifolds with bifurcate
Killing horizons, showing that the natural vacuum on Kruskal-extended
Schwarzschild restricted to the static exterior is a KMS state at
inverse temperature $\beta = 8\pi M$ with respect to Killing time. This
is the direct Lorentzian-side companion of Hartle–Hawking 1976's
Wick-rotation derivation.

**Path-integral version.** Hartle–Hawking 1976 (Phys. Rev. D 13, 2188)
gave the path-integral derivation that the Schwarzschild geometry has
imaginary-time period $8\pi M$ and that the corresponding Euclidean
propagator is thermal. This is the canonical Wick-rotation argument
that Sprint TD Track 4 reproduces in the GeoVac framework.

**Modern reviews.** Witten 2018 (Rev. Mod. Phys. 90, 045003;
arXiv:1803.04993) provides the modern QFT-level review of modular flow,
Tomita–Takesaki, Reeh–Schlieder, half-sided modular inclusions, and
their role in entanglement entropy. The modular Hamiltonian for the
half-Minkowski wedge is identified with the boost generator explicitly.

**NCG-side foundations.** Connes 1994 NCG book (Academic Press; Connes
1995 J. Math. Phys. 36; Connes–Marcolli AMS Colloquium 55, 2008)
provides the systematic NCG-side treatment of Tomita–Takesaki on
type-II/III von Neumann algebras. Connes' cocycle derivative shows the
modular automorphism's image in the outer automorphism group is
canonical (state-independent), which is the abstract version of the
geometric Bisognano–Wichmann statement.

**Pseudo-Riemannian / Krein-space NCG.** Strohmaier 2006, Bizi–Brouder–Besnard
2018 (J. Math. Phys. 59, 062303; arXiv:1611.07062), Franco–Eckstein 2014
(Rev. Math. Phys. 26, 1430007), van den Dungen 2016 (Math. Phys. Anal.
Geom. 19, 4) provide the Lorentzian extension framework in NCG. None
provides a Lorentzian propinquity / Gromov–Hausdorff convergence
theorem (see Track B scoping memo §2 for verified survey).

**Twisted-emergence approaches.** Devastato–Lizzi–Martinetti 2018
(JHEP 03, 089), Nieuviarts 2024–2025 (arXiv:2401.07848, 2502.18105,
2512.15450) provide twisted-spectral-triple alternatives that produce
Lorentzian signature without Wick rotation. Currently fresh; not
established as the standard prescription.

**Wald 1994** (QFT in Curved Spacetime, Univ. of Chicago Press)
provides the definitive textbook synthesis of Bisognano–Wichmann +
Sewell + Hartle–Hawking in the curved-spacetime QFT context.

### §4.2 What GeoVac contributes that is NOT in the published literature

The Bisognano–Wichmann reading itself is published. Sewell 1982 says
explicitly that Schwarzschild thermal physics IS modular flow on a
Lorentzian wedge algebra. What GeoVac contributes that is *not in
the published literature*:

**(a) Master Mellin engine M1 identification.** The framework
classifies the $2\pi$ in $\beta = 8\pi M$ as the M1 sub-mechanism of the
master Mellin engine $\mathcal{M}[\mathrm{Tr}(D^k \cdot e^{-tD^2})]$ at
$k = 0$. The published literature uses the Wick-rotation argument
directly without classification; the master Mellin engine is a
GeoVac-internal organizing principle (Paper 32 §VIII case-exhaustion
theorem; Paper 18 §III.7; Sprint MR-A/B/C). This embedding places
Sprint TD Track 4's Hawking-$T$ result inside the same case-exhaustion
theorem that classifies the $\pi$ content of the L2 propinquity rate
($4/\pi$, Paper 38), the spectral-action modular residual
($\sqrt{\pi}\,\mathbb{Q} \oplus \pi^2\,\mathbb{Q}$, Sprint MR-B), and
the Stefan–Boltzmann constant ($\pi^2/90$ as $M_1 \times M_2$ product,
Paper 35 §VIII).

**(b) Cross-track structural unification.** The structural identification
of M1 with both (i) the L2 propinquity rate $4/\pi$ on the spatial $S^3$
Hopf base (Paper 38) and (ii) the Hawking-$T$ $2\pi$ on the temporal
$S^1_\tau$ (Sprint TD Track 4) means Track D embeds two previously
separate structural results inside a single sub-mechanism. The
GH-convergence rate of the truncated metric spectral triple (Paper 38)
and the Hawking temperature on the Schwarzschild cigar (Sprint TD
Track 4) are *the same M1 sub-mechanism instantiated on different
sub-manifolds*. Under the Bisognano–Wichmann reading, both are also
$\mathrm{Vol}$ / Mellin-measure factors (Riemannian) and modular-flow
orbit periods (Lorentzian), via the Wick-rotation chain.

**(c) Catalogue placement.** The Bisognano–Wichmann reading provides
the published-physics anchor for adding Hawking $T$ to Paper 34 §V.B as
a machinery-witness row (error class C). This catalogue placement is a
GeoVac-internal organizational contribution, not a published claim.

**What we do NOT claim that would be new:**
- We do NOT claim a Lorentzian NCG-side proof of Bisognano–Wichmann
  inside the framework. That requires the multi-month Lorentzian
  propinquity work that Track B declared NO-GO.
- We do NOT claim the framework derives Hawking radiation from the
  packing axiom. $M$ is external Einstein-equations input
  (Sprint TD Track 4 V3 NEGATIVE).
- We do NOT claim Hawking $T = 1/(8\pi M)$ as a Paper 34 §V
  machine-precision match. The cigar geometry is not autonomously
  derived by GeoVac; we add Hawking $T$ to §V.B (off-precision) only as
  a machinery-witness row, error class C.

---

## §5. Forward implications

### §5.1 Does this open a 20th projection?

**No.** Bisognano–Wichmann does not introduce a new Layer-2 projection
in the Paper 34 sense. It is a *structural reading* of the existing
**Observation / temporal-window projection** (15th projection, Paper 34
§III.15). The reading says: when the temporal-window projection is
applied at the cigar's $\tau$-period $\beta = 8\pi M$, the resulting
$\beta$ IS the modular-flow period of a boost-class generator (Sewell
1982). This is an interpretation of an existing projection, not a new
projection.

**Consequence for Paper 34:** §III.15 (Observation/temporal-window)
gains a footnote naming the Bisognano–Wichmann reading.
§VIII Lorentz-boost open question gains a paragraph noting that Track D
closes the documentation loop while preserving the REQUIRES-EXTENSION
verdict for the *projection itself* (a literal Lorentz-boost projection
requires a Minkowski-signature spectral triple).

### §5.2 Does this strengthen Paper 32 §VIII case-exhaustion theorem?

**It strengthens by enrichment, not by adding a constraint.** M1 now has
two independent physical interpretations:

(a) Hopf-bundle measure factor on a sub-manifold. On the spatial $S^3$
Hopf base, this is $\mathrm{Vol}(S^2)/\pi^2 = 4/\pi$ (the L2 propinquity
rate constant). On the temporal $S^1_\tau$ of the cigar, this is
$\mathrm{Vol}(S^1)$ (the cigar-tip circumference $2\pi$ that gives
$\beta = 8\pi M$).

(b) Modular-flow / Lorentz-boost orbit period (via the Wick-rotation
chain Bisognano–Wichmann 1976 → Sewell 1982). On temporal
compactifications specifically, the $2\pi$ in M1 IS the period of the
modular automorphism / boost generator on the wedge algebra.

The two readings are not independent constraints — they are unified by
the Wick-rotation prescription. The case-exhaustion theorem itself is
unchanged; its Remark `rem:master_mellin_domain` (Sprint MR-A/B/C) gains
a Bisognano–Wichmann pendant naming the dual reading of M1.

### §5.3 New entries in §V or §V.B catalogues?

**Machine-precision row for Hawking $T$:** NO. $T_H = 1/(8\pi M)$ is not
a derived value of the framework; $M$ is external input. Adding this to
§V would overclaim.

**Off-precision row for Hawking $T$ machinery:** YES. A §V.B row reads:

> *Hawking temperature on Schwarzschild cigar $T_H = 1/(8\pi M)$
> reproduced via M1 mechanism (Sprint TD Track 4); structural
> correspondence to Bisognano–Wichmann modular flow / Lorentz-boost
> orbit period (Sewell 1982; Hartle–Hawking 1976) via standard
> Wick-rotation chain. Error class C: $M$ is external Einstein-equations
> input, not a framework-internal calibration mismatch.*

This row is a machinery-witness, not a precision match. It documents
that the framework's master-Mellin-engine M1 sub-mechanism is consistent
with published black-hole thermodynamics under Wick rotation.

### §5.4 Unruh temperature pendant

The same machinery applied to a Rindler wedge gives Unruh $T_U = a/(2\pi)$
for a uniformly accelerated observer at proper acceleration $a$. The
Wick-rotated Rindler wedge has imaginary-time period $\beta_{\rm Rindler}
= 2\pi/a$; M1 reproduces this from the circumference of the Rindler
$\tau$-circle.

This is a 1–2 day pendant follow-up: apply Sprint TD Track 4's
`matsubara_spectrum()` to Rindler with $\beta_{\rm Rindler} = 2\pi/a$;
verify M1 reproduces $T_U$; add §V.B row.

**Recommendation:** Document the corollary in this memo and Paper 32
§VIII remark; do NOT add an Unruh row to §V.B in Track D without
explicit framework-side computation of the Rindler $\beta = 2\pi/a$ from
M1. Sprint follow-up flag.

### §5.5 Specific 2–3 week follow-up tracks

| Track | Window | Deliverable | Priority |
|:---|:---:|:---|:---:|
| **D-extension: Operator-level $K$ on $\mathcal{T}_{n_{\max}}$** | 4–8 weeks | Construct modular Hamiltonian; verify $\sigma_{i \cdot 2\pi} = 1$ at finite $n_{\max}$; lift to GH limit. Closes falsifier #1; lifts correspondence to literal identification. | MEDIUM-HIGH |
| **D-companion: Unruh $T_U$ documentation pass** | 1–2 days | Apply Sprint TD Track 4 to Rindler with $\beta = 2\pi/a$; verify M1 reproduces $T_U$; add §V.B row. | LOW-MEDIUM |
| **D-deferred: Bisognano–Wichmann + entanglement Hamiltonian** | 2–4 weeks | Following Witten 2018 RMP, connect Bisognano–Wichmann modular Hamiltonian for half-$S^3$ algebra to Paper 27 entropy thread. Cross-thread connection. | LOW |
| **B-Candidate B (deferred): SO(4,1) conformal-action diagnostic** | 1–2 weeks | Compute SO(4,1) conformal Killing algebra on $S^3$ at finite $n_{\max}$; check Roothaan termination + Connes axioms preservation. | MEDIUM |

The principal recommendation is **D-extension** (4–8 weeks, MEDIUM-HIGH
priority): an operator-level proof of $\sigma_{i \cdot 2\pi} = 1$ at
finite $n_{\max}$, lifting the structural correspondence to literal
identification on the truncated spectral triple. This would make
Sprint TD Track 4 a *framework-internal* Bisognano–Wichmann theorem at
finite $n_{\max}$ — proven in the same framework where WH1 PROVEN
(GH-convergence on $S^3$) was established, via the same propinquity
machinery (Paper 38).

---

## §6. Honest scope statement

**What Track D claims:**
1. The $2\pi$ factor in Sprint TD Track 4's $T_H = 1/(8\pi M)$ is the
   M1 Hopf-base measure / $S^1_\tau$ circumference of the master Mellin
   engine.
2. Under the standard Wick-rotation chain (Hartle–Hawking 1976 →
   Sewell 1982 → Bisognano–Wichmann 1976), this $2\pi$ IS the
   modular-flow / Lorentz-boost orbit period of a stationary observer
   outside the Schwarzschild horizon.
3. The framework's M1 mechanism is therefore the Wick-rotated image of
   the Bisognano–Wichmann modular-flow mechanism.
4. This is a **structural correspondence**, not a literal identification.

**What Track D does NOT claim:**
- Operator-level proof of $\sigma_{i \cdot 2\pi} = 1$ on the truncated
  metric spectral triple (this is falsifier #1; 4–8 week follow-up).
- A Lorentzian-signature spectral triple (Track B verdict NO-GO; 9–18
  months).
- A Lorentzian propinquity theory (no published prescription; 6–12
  months original NCG-math).
- A derivation of the cigar geometry from GeoVac axioms (Sprint TD
  Track 4 V3 NEGATIVE).
- Bekenstein–Hawking entropy from the horizon $S^2$ (Sprint TD Track 4
  V2 OUT-OF-SCOPE).
- Black-hole-vs-GeoVac entropy comparison (Sprint TD Track 5 territory;
  different category).

**Is this oversold?** No. The headline is **structural correspondence**
(qualified explicitly as not literal identification); the underlying
physics is published QFT (Bisognano–Wichmann 1976, Sewell 1982,
Hartle–Hawking 1976); the framework-internal contribution is the M1
master Mellin engine identification, embedding Sprint TD Track 4 in
Paper 32 §VIII case-exhaustion theorem; three specific falsifiers are
named with priorities and operationalizations. CLAUDE.md §1.5 rhetoric
rule applied throughout: dual-description framing (the framework's M1
mechanism produces $\beta = 8\pi M$; the published Wick-rotation
prescription identifies this $\beta$ with the modular-flow period; the
two readings are mutually consistent and the choice between them is left
to the reader); no ontological priority claimed; lead with the
computational result (Sprint TD Track 4's sympy-residual-zero
reproduction of Hawking $T$).

**Curve-fit-audit discipline:** M1 already had ONE reading (Hopf-bundle
measure). Adding a second reading (boost orbit period) is structural
enrichment, not a fit to data. The Bisognano–Wichmann reading is
constrained by Sewell 1982 (not by Track D); we do not adjust any
parameter to make the reading work; the $2\pi$ is determined by the
cigar geometry's regularity condition and by the modular-flow KMS
analytic period independently.

**Diagnostic-before-engineering discipline:** Three falsifiers are
named. Falsifier #1 is the principal computational test that would
either confirm the structural correspondence as literal identification
or refute the Bisognano–Wichmann reading. Operationalization is
specified (extending the existing Connes-distance SDP framework to
modular flow). Priority MEDIUM-HIGH for follow-up.

---

## §7. Paper update record

Three paper updates applied directly per CLAUDE.md §13.8 (paper edits
are autonomous when supported by the analysis):

### §7.1 Paper 32 §VIII update

`papers/group1_operator_algebras/paper_32_spectral_triple.tex` — added two paragraphs:

1. After Remark `rem:master_mellin_domain` (Sprint MR-A/B/C
   mechanism-as-domain partition), a new paragraph
   `rem:bisognano_wichmann_reading` naming the dual physical reading of
   M1 on temporal compactifications: Hopf-base measure (Riemannian) ↔
   modular-flow / boost-orbit period (Lorentzian, via Wick-rotation
   chain Sewell 1982 + Bisognano–Wichmann 1976). Explicit Sprint TD
   Track 4 cross-reference; explicit Track D scope qualifier
   (structural correspondence, not literal identification); explicit
   falsifier #1 reference (operator-level closure on $\mathcal{T}_{n_{\max}}$
   as the natural follow-up).

2. Bibliography additions: `bisognano_wichmann1975`,
   `bisognano_wichmann1976`, `sewell1982`, `hartle_hawking1976`,
   `witten2018_rmp`. (Connes 1994 NCG book is already referenced via
   `connes_marcolli2008` AMS Colloquium volume 55; no new entry needed.)

### §7.2 Paper 35 §VIII update

`papers/group6_precision_observations/paper_35_time_as_projection.tex` — added a new
subsection §VIII after `subsec:tensor_scope` titled
"Bisognano–Wichmann reading of the temporal-window projection":

- Documents the structural correspondence between the temporal-window
  projection (Paper 34 §III.15) and the modular-flow / Lorentz-boost
  orbit period via Bisognano–Wichmann 1976 + Sewell 1982.
- Gives the explicit cigar example: $\beta = 8\pi M$ from M1 / $S^1_\tau$
  circumference reproduces Hawking $T_H = 1/(8\pi M)$, structurally
  identified with the modular-flow period of a stationary observer
  outside the Schwarzschild horizon.
- Documents the Unruh corollary: $\beta_{\rm Rindler} = 2\pi/a$, $T_U =
  a/(2\pi)$, same M1 mechanism.
- Honest scope statement: structural correspondence, not literal
  identification; Lorentzian propinquity is not in scope (NO-GO per
  Track B).
- Cross-references: Paper 32 §VIII, Paper 38 (L2 propinquity rate
  $4/\pi$ as the spatial M1 instantiation).
- Bibliography additions: `bisognano_wichmann1976`, `sewell1982`,
  `hartle_hawking1976`.

### §7.3 Paper 34 updates

`papers/group6_precision_observations/paper_34_projection_taxonomy.tex` — three small
edits:

1. **§III.15** (Observation/temporal-window projection) gains a footnote
   naming the Bisognano–Wichmann reading: "Under the standard
   Wick-rotation chain (Hartle–Hawking 1976; Sewell 1982; Bisognano–Wichmann
   1976), the temporal-window projection's $\beta$ IS the period of the
   modular-flow / Lorentz-boost orbit on a wedge or static-exterior
   algebra. See Paper 32 §VIII Remark `rem:bisognano_wichmann_reading`
   and Paper 35 §VIII for the framework-level structural correspondence;
   Track D memo `debug/bisognano_wichmann_track_d_memo.md` for the
   detailed Sprint TD Track 4 identification."
2. **§V.B** gains an off-precision Hawking-$T$ row (machinery-witness,
   error class C; see §5.3 above).
3. **§VIII Lorentz-boost open question** gains a paragraph noting that
   Track D closes the documentation loop on the Bisognano–Wichmann
   reading while preserving the REQUIRES-EXTENSION verdict for a
   literal Lorentz-boost projection (which would require Minkowski-signature
   spectral triple work).

---

## §8. Result (30-second read)

**Structural identification:** STRUCTURAL CORRESPONDENCE (not literal
identification). The $2\pi$ factor in Sprint TD Track 4's
$T_H = 1/(8\pi M)$ is the M1 Hopf-base measure / $S^1_\tau$ circumference
of the master Mellin engine. Under the standard Wick-rotation chain
(Hartle–Hawking 1976 → Sewell 1982 → Bisognano–Wichmann 1976), this
$2\pi$ IS the period of the modular-flow / Lorentz-boost orbit of a
stationary observer outside the Schwarzschild horizon.

**Top falsifier:** Compute the modular Hamiltonian $K$ on the truncated
metric spectral triple $\mathcal{T}_{n_{\max}}$ (Paper 38 framework) for
a wedge-restricted Camporesi–Higuchi state; verify
$\sigma_{i \cdot 2\pi} = 1$ at finite $n_{\max}$. 4–8 week follow-up;
priority MEDIUM-HIGH; would lift correspondence to literal identification
inside the framework.

**Top 2 forward implications:**
1. Paper 32 §VIII case-exhaustion theorem strengthens by enrichment: M1
   has TWO physical interpretations (Hopf-bundle measure / boost-orbit
   period via Bisognano–Wichmann), with Wick rotation identifying them.
   The case-exhaustion theorem itself is unchanged.
2. Paper 34 §V.B gains a Hawking-temperature off-precision row
   (machinery-witness, error class C). Unruh $T_U$ is a 1–2 day pendant
   follow-up; not pursued in this sprint.

**What this enables in the catalogue:**
(a) Bisognano–Wichmann reading of Sprint TD Track 4 documented in three
papers without leaving Riemannian signature; (b) cross-thread connection
between thermal physics (Sprint TD), GH-convergence (Paper 38), and the
case-exhaustion theorem (Paper 32 §VIII) made explicit via M1's dual
reading; (c) Hawking $T$ enters the catalogue as an off-precision row;
(d) framework-internal Bisognano–Wichmann (literal identification at
finite $n_{\max}$) is named as a falsifier with operationalization,
priority MEDIUM-HIGH.

---

## §9. Files

- `debug/bisognano_wichmann_track_d_memo.md` — this memo (~7000 words).
- `debug/data/bisognano_wichmann_track_d_data.json` — structured data:
  literature citations with verification status, structural-identification
  claim, falsifier list, ranked forward implications, paper-update record.
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` — §VIII Remark
  `rem:bisognano_wichmann_reading` added; bibliography +5 entries.
- `papers/group6_precision_observations/paper_35_time_as_projection.tex` — new §VIII
  subsection on Bisognano–Wichmann reading; bibliography +3 entries.
- `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` — §III.15
  footnote; §V.B Hawking-$T$ row; §VIII Lorentz-boost-open-question
  paragraph extension.

No production `geovac/` code modified.

---

## §10. Bibliography (verified, with arXiv IDs and verification method)

1. **Bisognano, J.J. and Wichmann, E.H.** "On the duality condition for
   a Hermitian scalar field." *J. Math. Phys.* **16**, 985–1007 (1975).
   *Verified via:* WebSearch + ADS UI 1975JMP....16..985B + osti.gov.

2. **Bisognano, J.J. and Wichmann, E.H.** "On the duality condition for
   quantum fields." *J. Math. Phys.* **17**, 303–321 (1976).
   *Verified via:* WebSearch + AIP listing pubs.aip.org/aip/jmp/article-abstract/17/3/303.
   *Title correction:* the Track B scoping memo had this as "On the duality
   condition for a Hermitian scalar field" (which is the 1975 paper title).

3. **Sewell, G.L.** "Quantum fields on manifolds: PCT and gravitationally
   induced thermal states." *Annals of Physics* **141**, 201–224 (1982).
   *Verified via:* WebSearch + ScienceDirect + reprints.gravitywaves.com PDF.
   *Year correction:* the Track B scoping memo had 1980; the actual
   publication year is 1982.

4. **Hartle, J.B. and Hawking, S.W.** "Path-integral derivation of
   black-hole radiance." *Phys. Rev. D* **13**, 2188–2203 (1976).
   *Verified via:* WebSearch + APS Phys Rev D listing + ADS 1976PhRvD..13.2188H.

5. **Witten, E.** "APS Medal for Exceptional Achievement in Research:
   Invited article on entanglement properties of quantum field theory."
   *Rev. Mod. Phys.* **90**, 045003 (2018). arXiv:1803.04993.
   *Verified via:* WebSearch + APS RMP 90.045003 + arXiv listing.

6. **Jacobson, T. and Parentani, R.** "Horizon entropy." *Foundations of
   Physics* **33**, 323–348 (2003). arXiv:gr-qc/0302099.
   *Verified via:* WebSearch + Springer Nature link + arXiv listing.

7. **Connes, A.** *Noncommutative Geometry*. Academic Press, 1994.
   *Verified via:* WebSearch + alainconnes.org.
   *Already cited in Paper 32 via Connes 1995 / Connes-Marcolli 2008.*

8. **Wald, R.M.** *Quantum Field Theory in Curved Spacetime and Black
   Hole Thermodynamics*. University of Chicago Press, 1994.
   *Standard graduate text; widely cited.*

**Phantom citation flagged in Track B (carried over):**
"Hawkins–Skoda 2018/2020 Krein-space spectral triples" was listed in the
original directive but could not be located in any database. Track B
flagged as phantom; Track D does not cite. Closest published matches
are van den Dungen 2016 (Math. Phys. Anal. Geom. 19, 4) and
Brouder–Besnard–Bizi 2017 (same journal volume).
