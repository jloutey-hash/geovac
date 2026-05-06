# Almost-Commutative Extension of the GeoVac Spectral Triple — Scoping Memo

**Sprint:** Step 3 of the May 3–4 2026 PI three-step spectral-triple commitment
(Step 1 = R2.5 GH convergence, Step 2 = real structure J, **Step 3 = this scoping**).
**Mode:** SCOPING ONLY. No construction, no theorem, no implementation.
**Date:** 2026-05-06.
**Status:** Design memo + interface stub. Tests deferred.
**Final venue (if construction succeeds):** A future Sprint H1 (Higgs-from-inner-fluctuation).
The natural paper home for the eventual result is Paper 32 §VIII.B (closes the G2 / G3
half of the four-gap analysis) or, if the construction grows into a self-contained
result on the S³ U(1) × SU(2) electroweak co-location, a separate Paper 37.
**Per CLAUDE.md §1.5 / §13.5:** this memo does NOT modify Paper 32 §VIII or the
combination-rule status of Paper 2; it scopes a future sprint.

---

## §0. One-paragraph verdict

Reachable, but conditional. The minimal almost-commutative extension
$\mathcal{A}_{\mathrm{GV}} \otimes \mathcal{A}_F$ with the electroweak finite algebra
$\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$ can be written down today on the
GeoVac infrastructure already in place (R3.5 full-Dirac chirality grading, Paper 32
real structure $J$, scalar function algebra $\mathcal{A}_{\mathrm{GV}}$). Inner
fluctuations $D \mapsto D + \omega + J\omega J^{-1}$ in this extension would
formally split into a gauge-1-form sector (recovering Papers 25 / 30) **plus** a
scalar Higgs sector living in the off-diagonal $(\mathbb{C}, \mathbb{H})$
component. The critical missing ingredient identified in the Marcolli–vS vs
Connes–Chamseddine literature is the **off-diagonal block of the finite Dirac
$D_F$ that connects the $\mathbb{C}$ and $\mathbb{H}$ summands**: this is what
turns "Yang–Mills without Higgs" into "Yang–Mills with Higgs". The $D_F$
off-diagonal component plays the role of the **Yukawa coupling matrix** in
Connes–Chamseddine; its zero-mode value is the Higgs vev. **GeoVac does not
yet supply a natural choice of this off-diagonal block.** Whether one can be
selected from GeoVac structure (rather than imposed) is the open question the
eventual Sprint H1 must answer; we identify two candidate ingredients and one
clean falsifier.

---

## §1. Scope and dependencies

### 1.1. What we have today (independent of Tracks 1 and 2)

The minimal data needed to write down an almost-commutative spectral triple
$\mathcal{T}_{\mathrm{GV}} \otimes \mathcal{T}_F$ exists in the codebase:

- $\mathcal{A}_{\mathrm{GV}}$: `geovac/operator_system.py` (commutative diagonal
  multipliers; truncated operator system with Avery overlap).
- $\mathcal{H}_{\mathrm{GV}}$: `geovac/spinor_operator_system.py` (Weyl sector)
  and `geovac/full_dirac_operator_system.py` (full Dirac, both chiralities).
- $D_{\mathrm{GV}}$: Camporesi–Higuchi spectrum, both `truthful` and `offdiag`
  variants, in `full_dirac_operator_system.py`.
- $J_{\mathrm{GV}}$: Paper 32 Prop. 4 (continuum CH charge conjugation) plus
  the finite implementation in `geovac/spectral_triple.py` (`j_type='kramers'`
  for $J^2 = -\mathbb{1}$ on half-integer $j$).
- Wilson SU(2) machinery: `geovac/su2_wilson_gauge.py` (already an inner-
  fluctuation-flavored construction at the level of Wilson links; per Perez-
  Sanchez 2024, it is YM **without** Higgs by default).

Everything in §1.1 is needed to **write down** the extension. No extension
work has been done; the building blocks have been kept ready.

### 1.2. What depends on Track 1 (R2.5 L4 Berezin reconstruction)

R2.5 is the GH-convergence proof of the truncated operator system to
$C^\infty(S^3)$. L1', L2, L3 are closed (offdiag operator system; central
Fejér kernel; Lipschitz spinor bound). L4 (Berezin reconstruction via Hawkins)
and L5 (propinquity assembly) are open.

The Higgs construction does **not formally require** R2.5 L4 to write down;
inner fluctuations are defined at finite $n_{\max}$. **However, the
*physical interpretation* of the Higgs scalar field — as a section of the
trivial $\mathbb{H}$-bundle over $S^3$ — requires that the algebraic content
of $\mathcal{A}_{\mathrm{GV}}$ converge in the GH sense to functions on
$S^3$.** Without R2.5 the Higgs is well-defined as a finite matrix; with R2.5
it converges to an honest scalar field. Tracks 1 and 3 are therefore
**parallel-independent at the construction level, sequentially-dependent at
the interpretation level**.

**Upshot:** The Sprint H1 mathematical computation can proceed before R2.5
L4 closes. The eventual paper-level claim ("the GeoVac Higgs is a continuum
scalar field on $S^3$") needs L4–L5 closed.

### 1.3. What depends on Track 2 (real structure J at finite $n_{\max}$)

This is the **load-bearing dependency**. The Connes–Chamseddine inner
fluctuation formula on a real spectral triple is

$$
D \;\longmapsto\; D \;+\; \omega \;+\; \epsilon' J \omega J^{-1},
\qquad \omega = \sum_i a_i [D, b_i],\; a_i, b_i \in \mathcal{A}.
$$

The $J \omega J^{-1}$ term encodes the **opposite-algebra action** that makes
the Higgs sector real. Without a verified $J$ at finite $n_{\max}$ with the
correct sign $\epsilon = -1$ (KO-dim 3 of $S^3$), the inner fluctuation
formula cannot be evaluated. Paper 32 §V Prop. 4 establishes $J$ in the
**continuum**; Paper 32 §VIII Open Question Q2 explicitly flags whether the
finite-$n_{\max}$ restriction satisfies the same $J^2 = -\mathbb{1}$ exactly.

**Upshot:** Sprint H1 cannot run before Track 2 (real structure verification
at $n_{\max} \le 5$) closes. If Track 2 finds a sign mismatch at finite
$n_{\max}$ that resolves only in the GH limit, the Higgs construction at
finite $n_{\max}$ acquires a controllable error term that needs to be tracked.

### 1.4. Critical-path summary

```
Track 2 (J at finite n_max, KO-dim 3 verification) ────► PREREQUISITE for H1
Track 1 (R2.5 L4 Berezin)                          ────► PARALLEL,
                                                          required only for
                                                          continuum-field
                                                          INTERPRETATION
This memo (scoping)                                ────► DONE today
Sprint H1 (compute the extension)                  ────► WAIT for Track 2
Paper-level Higgs claim                             ────► WAIT for Tracks 1 + 2
```

---

## §2. The minimal extension: $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$

### 2.1. Algebra and Hilbert space

The minimal electroweak almost-commutative extension uses the finite algebra

$$
\mathcal{A}_F \;=\; \mathbb{C} \;\oplus\; \mathbb{H},
$$

a $\mathbb{R}$-algebra of real dimension $1 + 4 = 5$ (or $\mathbb{C}$-dimension
$1 + 2 = 3$ when the quaternionic action is realized as $2 \times 2$ complex
matrices: $\mathbb{H} \otimes_\mathbb{R} \mathbb{C} \cong M_2(\mathbb{C})$).
This is the **electroweak slice** of the Connes–Chamseddine SM algebra
$\mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$; we explicitly **defer
the $M_3(\mathbb{C})$ color factor** as out-of-scope for Sprint H1 (it would
require committing to a cross-manifold extension to the $S^5$ Bargmann graph,
which Paper 32 §VIII.B G4 names as the dominant gap).

The combined algebra is

$$
\mathcal{A} \;=\; \mathcal{A}_{\mathrm{GV}} \otimes \mathcal{A}_F
\;=\; \mathbb{C}^{V_{\mathrm{Fock}}} \otimes (\mathbb{C} \oplus \mathbb{H})
\;\cong\; \mathbb{C}^{V_{\mathrm{Fock}}} \oplus M_2(\mathbb{C})^{V_{\mathrm{Fock}}}.
$$

The combined Hilbert space is

$$
\mathcal{H} \;=\; \mathcal{H}_{\mathrm{GV}} \otimes \mathcal{H}_F,
$$

where $\mathcal{H}_F$ is a finite Hilbert space carrying a faithful
$*$-representation of $\mathcal{A}_F$. The minimal choice for the electroweak
sector (matching Connes–Chamseddine's reduced model) is

$$
\mathcal{H}_F \;=\; \mathbb{C}^4_L \;\oplus\; \mathbb{C}^2_R
\quad \text{or} \quad
\mathbb{C}^2_L \;\oplus\; \mathbb{C}^2_R
$$

depending on whether one carries a left-handed weak-isospin doublet plus the
two right-handed singlets (the "leptonic minimum" $L = (\nu, e)$,
$R = e_R$, plus its right-handed neutrino partner $\nu_R$). For the
**electroweak-only-no-color, single-generation** scoping, $\mathcal{H}_F$ has
$\mathbb{C}$-dimension 4: two left-handed (one weak-isospin doublet) and two
right-handed.

The product

$$
\dim_\mathbb{C} \mathcal{H} \;=\; \dim \mathcal{H}_{\mathrm{GV}} \cdot 4
$$

stays small at finite $n_{\max}$: at $n_{\max} = 3$, $\dim
\mathcal{H}_{\mathrm{GV}} = 40$, so $\dim \mathcal{H} = 160$. Tractable.

### 2.2. The Dirac operator on $\mathcal{H}$

The standard CC ansatz is

$$
D \;=\; D_{\mathrm{GV}} \otimes \mathbb{1}_F \;+\; \gamma_{\mathrm{GV}} \otimes D_F,
$$

where $\gamma_{\mathrm{GV}}$ is a chirality grading on $\mathcal{H}_{\mathrm{GV}}$
(GeoVac supplies one via the R3.5 chirality $\chi = \pm 1$) and $D_F$ is a
Hermitian operator on $\mathcal{H}_F$.

**Critical observation:** $D_F$ is the source of the Higgs. Its non-zero
*off-diagonal* matrix elements between the $\mathbb{C}$ and $\mathbb{H}$
summands of $\mathcal{A}_F$ are the **Yukawa coupling matrix entries**;
inner fluctuations turn these off-diagonal entries into the Higgs scalar
field on $S^3$. If $D_F$ is block-diagonal w.r.t. the $\mathbb{C} \oplus
\mathbb{H}$ decomposition, no Higgs appears (this is the structural reason
Marcolli–vS gauge networks give YM-without-Higgs: they implicitly take
$D_F$ block-diagonal because the "fiber Dirac" comes only from the gauge
links).

### 2.3. Inner fluctuations at finite $n_{\max}$

In the Connes formalism, the connection 1-form is

$$
\omega \;=\; \sum_i a_i [D, b_i], \qquad a_i, b_i \in \mathcal{A}.
$$

For the product algebra $\mathcal{A} = \mathcal{A}_{\mathrm{GV}} \otimes
\mathcal{A}_F$ and product Dirac $D = D_{\mathrm{GV}} \otimes \mathbb{1}_F +
\gamma_{\mathrm{GV}} \otimes D_F$:

$$
[D, b_i] \;=\; [D_{\mathrm{GV}}, b_i^{\mathrm{GV}}] \otimes b_i^F
\;+\; \gamma_{\mathrm{GV}} b_i^{\mathrm{GV}} \otimes [D_F, b_i^F].
$$

The first term, summed against $a_i$, produces a **gauge 1-form** on
$\mathcal{A}_{\mathrm{GV}}$ tensored with an algebra element on
$\mathcal{A}_F$. Its restriction to traceless skew-Hermitian elements is
**exactly Papers 25 / 30**: $U(1)$ links from the abelian part, $SU(2)$
links from the quaternionic part. This recovers the Wilson-YM content
already in the code.

The second term, summed against $a_i$, produces a **scalar Higgs field**
$\Phi \in \mathcal{A}_{\mathrm{GV}} \otimes \mathrm{Hom}(\mathbb{C},
\mathbb{H})$. The Higgs scalar is the **off-diagonal block** of $\Phi$ as a
matrix in $\mathcal{H}_F$. Its zero-mode value (the constant function on
$\mathcal{A}_{\mathrm{GV}}$) is the Higgs vev; its fluctuations along
non-constant elements of $\mathcal{A}_{\mathrm{GV}}$ are the propagating
Higgs particle.

### 2.4. Real structure on the extension

For the combined real structure:

$$
J \;=\; J_{\mathrm{GV}} \otimes J_F,
$$

with $J_F$ the standard finite real structure on $\mathcal{A}_F = \mathbb{C}
\oplus \mathbb{H}$ (Connes–Marcolli 2008 Ch. 13). Sign multiplication: $J^2 =
J_{\mathrm{GV}}^2 \cdot J_F^2 = (-1)(-1) = +1$ in the electroweak slice. The
KO-dimension shifts: $\mathrm{KO}(\mathcal{T}) = \mathrm{KO}(\mathcal{T}_{
\mathrm{GV}}) + \mathrm{KO}(\mathcal{T}_F) \pmod 8 = 3 + 6 = 9 \equiv 1
\pmod 8$ for the leptonic-electroweak slice (Connes–Marcolli 2008 Table
13.1). KO-dim 1 is the **Lorentzian-flavored** dimension; this is consistent
with CC's observation that the SM almost-commutative geometry has
KO-dim 6 for the full algebra including color.

The order-one condition for the extension reads

$$
[[D, \pi(a)], J \pi(b) J^{-1}] \;=\; 0 \quad \forall a, b \in \mathcal{A}.
$$

For the gauge-1-form part this is automatic in the abelian-extended version.
For the off-diagonal Higgs part it is **nontrivial**: it pins the form of
$D_F$ and is the ingredient that, in CC, derives the Higgs Yukawa structure.

---

## §3. Critical reading: Marcolli–vS 2014 vs Connes–Chamseddine 1996/2010

### 3.1. The two frameworks side by side

| | Marcolli–vS gauge networks (2014) | Connes–Chamseddine SM (1996/2010) |
|---|---|---|
| Vertex algebra | $C(V)$ on a vertex set | $C^\infty(M) \otimes \mathcal{A}_F$, $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$ |
| Edge / gauge data | Connection $U_e \in G$ on edges | Inner fluctuation $\omega = \sum a_i [D, b_i]$ |
| Dirac structure | Hopping Dirac on the graph | $D = D_M \otimes \mathbb{1}_F + \gamma_M \otimes D_F$ |
| Continuum limit | Wilson YM **without Higgs** (Perez-Sanchez 2024) | YM + Higgs + gravity (Chamseddine–Connes 2010) |
| Higgs source | None by default | Off-diagonal block of $D_F$ between $\mathbb{C}$ and $\mathbb{H}$ |

### 3.2. The missing ingredient

The Perez-Sanchez 2024 / 2025 corrections (arXiv:2401.03705, arXiv:2508.17338)
clarified that Marcolli–vS gauge networks do **not** automatically reproduce
the CC Higgs because **the gauge-network Dirac operator is hopping-only**:
its action factorizes through edges that carry gauge data, with no internal
fiber Dirac coupling. The fiber Dirac $D_F$ in CC is what provides the
off-diagonal $\mathbb{C} \leftrightarrow \mathbb{H}$ matrix element that
fluctuations turn into the Higgs.

**The missing ingredient for GeoVac to acquire a Higgs is therefore:**

> A non-zero $D_F$ on the finite fiber $\mathcal{A}_F = \mathbb{C} \oplus
> \mathbb{H}$ with non-trivial off-diagonal blocks between the two
> summands, and a justification that the choice of $D_F$ is selected by
> GeoVac structure rather than imposed.

### 3.3. Two candidate sources of $D_F$ from GeoVac structure

**Candidate A (R3.5 chirality + scalar fiber bridging):** The R3.5 full Dirac
already provides a chirality grading $\chi = \pm 1$ on $\mathcal{H}_{\mathrm{GV}}$.
The "offdiag" CH variant of `geovac/full_dirac_operator_system.py` carries
**native cross-chirality coupling** via the SO(4) E1 selection rule. If we
identify the chirality grading with the left/right asymmetry of $\mathcal{H}_F
= \mathbb{C}^2_L \oplus \mathbb{C}^2_R$ via an embedding $\mathcal{H}_F
\hookrightarrow \mathcal{H}_{\mathrm{GV}}$ at fixed $(n_{\max}, l)$, the
offdiag CH structure already supplies a candidate $D_F$ off-diagonal block.
**Status:** Speculative. It is not yet clear whether the offdiag CH coupling
is a Higgs or a gauge fluctuation.

**Candidate B (Sturmian / DUCC bridging):** GeoVac's DUCC downfolding
(`geovac/downfolding.py`) and Sturmian basis (`geovac/sturmian_solver.py`)
both include effective two-body operators with off-diagonal core-valence
structure. These could potentially supply $D_F$ off-diagonal blocks on a
two-shell extension, with the "Higgs vev" being the energy gap between
shells. **Status:** More speculative. This route maps the Higgs onto an
effective inter-shell coupling and would need a non-trivial argument for
why it is the unique natural choice.

### 3.4. Why this is the load-bearing question

If $D_F$ off-diagonal can be selected from GeoVac structure, then the
electroweak Higgs appears as the inner-fluctuation natural completion of
Papers 25 / 30 (closing G2 within the electroweak sector). If $D_F$
off-diagonal must be imposed by hand, then **GeoVac is on the
Marcolli–vS-without-Higgs side of the 2014/2024 distinction**, and the
electroweak Higgs would be put in by hand in any GeoVac-SM construction —
which would be honest but would not constitute a derivation.

---

## §4. Expected scalar/vector structure on the GeoVac graph

Conditional on a non-trivial $D_F$ being selected (one of the candidate
sources of §3.3), the inner-fluctuation formula on the extended triple
predicts the following structure on the Fock-projected $S^3$ graph:

- **Gauge bosons.** From the diagonal-in-fiber part of $\omega$: a
  $\mathfrak{u}(1)$-valued 1-cochain $A_e^{Y}$ (hypercharge) and an
  $\mathfrak{su}(2)$-valued 1-cochain $W_e^a$ (weak isospin), both living
  on the same $S^3$ Fock graph as Papers 25 and 30. The gauge-coupling
  ratios $g$, $g'$ are inherited from the relative scales of the diagonal
  $\mathbb{C}$ and $\mathbb{H}$ blocks of $\mathcal{A}_F$.

- **Higgs scalar.** From the off-diagonal-in-fiber part: a complex doublet
  scalar $\Phi(v) = (\phi^+(v), \phi^0(v))$ on each Fock node $v$, valued
  in $\mathbb{C}^2 \subset \mathbb{H}$. At finite $n_{\max}$ this is a
  finite-dimensional complex vector space of dimension $2 \cdot
  N_{\mathrm{Fock}}$. In the GH limit (R2.5 L4–L5) it becomes a section of
  the trivial $\mathbb{C}^2$ bundle over $S^3$.

- **Higgs potential.** The leading-order spectral-action expansion gives
  $V(\Phi) = a |\Phi|^2 + b |\Phi|^4 + c |\nabla \Phi|^2$ where $a, b, c$
  are determined by spectral-action coefficients. The Mexican-hat shape
  (negative $a$, positive $b$) is generic in CC; the GeoVac question is
  whether the rationals/transcendentals in the coefficients reduce to
  GeoVac's exchange-constant taxonomy (Paper 18) — see §6 below.

- **Yukawa couplings.** Cross-terms $\bar\psi_L \Phi \psi_R$ from the
  $\gamma_{\mathrm{GV}} \otimes D_F$ part of the Dirac. Mass eigenvalues
  read off after Higgs symmetry breaking as $m_f = y_f \langle \Phi
  \rangle$; in GeoVac, $y_f$ would be entries of the $D_F$ off-diagonal
  block. This is exactly the place the candidate sources of §3.3 inject
  GeoVac structure.

- **Higgs vev.** A constant value $\langle \Phi \rangle = v_0 \in
  \mathbb{C}^2$ at the spectral-action minimum. In CC this comes out
  as $v_0 \sim \Lambda$ (the spectral cutoff); for GeoVac the natural
  cutoff is set by $n_{\max}$, and the Connes-vS spectral-action ceiling
  theorem on $S^3$ (Paper 28 §spectral_action) gives
  $\Lambda_\infty \approx 3.71$ at $n_{\max} \to \infty$. This is in
  natural ($S^3$ unit-radius) units; conversion to physical units would
  need a separate calibration.

---

## §5. The natural negative — what would prove Higgs construction fails on GeoVac

A clean falsifier for the Higgs-as-inner-fluctuation thread on GeoVac:

> **Statement.** Show that for *every* Hermitian $D_F$ on $\mathbb{C} \oplus
> \mathbb{H}$ derivable from GeoVac structure (§3.3 Candidate A or B), the
> resulting inner fluctuation produces *only* gauge 1-forms (Papers 25 / 30
> recovered) and the off-diagonal Higgs sector is identically zero.

This would close the Higgs direction by formally elevating Marcolli–vS-
without-Higgs from an external published correction to a GeoVac-internal
theorem. The mechanism would be:

- Either (a) the chirality grading $\chi$ commutes with every fiber-bridging
  candidate $D_F$, so $\gamma_{\mathrm{GV}} \otimes D_F$ has no off-diagonal
  matrix elements between $\mathbb{C}$ and $\mathbb{H}$;
- Or (b) the order-one condition forces all candidate $D_F$ to be block-
  diagonal in $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$, eliminating the
  Higgs sector before any spectral-action computation.

Either form of (a) or (b) would close the direction at the structural level
and generalize the Sprint 4H Track SM-B negative (no shell-to-generation
map) to a Higgs-sector negative.

A weaker (but still informative) negative would be:

> **Statement weak.** Show that the natural candidate $D_F$ from §3.3
> Candidate A (R3.5 offdiag CH) reproduces a Yukawa with all-equal couplings
> $y_e = y_\mu = y_\tau$, contradicting the SM mass hierarchy.

This would be a partial-positive-with-fatal-sub-feature, analogous to the
Sprint 4H Track SM-D positive on Δ⁻¹ = $g_3^{\mathrm{Dirac}}$ that did not
extend to derivation of $K$.

---

## §6. Cross-checks against existing GeoVac structural results

The Higgs direction interacts with several existing GeoVac structural results
in ways that constrain the possible sprint outcomes.

### 6.1. Paper 18 (exchange constants)

In CC, the Higgs vev is a calibration scale ($\Lambda$). In Paper 18's
taxonomy, this is a **calibration exchange constant** — the same tier as
$\pi$ in Paper 2. If the GeoVac Higgs construction succeeds, the Higgs vev
should appear as a Paper 18 calibration constant; if it appears as
intrinsic-rational, Paper 18 needs a new tier; if it appears as
intrinsic-transcendental, the case-exhaustion theorem of Paper 32 §VIII
needs an extension.

### 6.2. Paper 35 (time as projection)

Paper 35's load-bearing principle is "$\pi$ enters iff continuous integration
over a temporal/spectral parameter." The Higgs vev in CC enters via the
spectral-action heat-kernel computation, which **does** integrate over a
spectral parameter. So a Higgs vev with $\pi$ in it is **expected** under
Paper 35; a Higgs vev that is rational would be a small surprise that should
trigger a Paper 35 audit.

### 6.3. Paper 34 (projection taxonomy)

Sprint H1 introduces a candidate **16th projection** (almost-commutative
fiber promotion) on top of Paper 34's existing fifteen. This needs to be
tagged: it adds a variable (the fiber index $\sigma \in \{1, 2\}$ for
$\mathbb{C} \oplus \mathbb{H}$), it adds a dimension (1 internal dimension
in CC's bookkeeping), and its transcendental signature is determined by
the spectral-action expansion (likely $\pi$-bearing via the Mellin-of-heat-
kernel mechanism of Paper 18 §III.7).

### 6.4. WH4 (four-way $S^3$ coincidence)

WH4 says $S^3$ plays four roles: Fock projection, Hopf base, CH spin
carrier, SU(2) gauge manifold. The Higgs construction adds a **fifth** role:
$S^3$ as the base manifold of the Higgs $\mathbb{C}^2$ bundle. This is
formally consistent with the four-way coincidence; whether it also
*forces* the doublet structure of the Higgs is a Sprint H1 question.

### 6.5. R3.5 chirality grading

Paper 32 §VIII.B G3 explicitly notes: "The chirality grading introduced in
Sprint TS Track E2 (R3.5) is a $\mathbb{Z}_2$-grading of the spinor bundle,
not the SM weak-isospin chirality." This is the place Sprint H1 has the
most constrained move: either

- (i) an embedding $\chi_{\mathrm{R3.5}} \cong \chi_{\mathrm{weak isospin}}$
  is selected by GeoVac structure (closing G3 simultaneously with G2), or
- (ii) the two chiralities are independent and the weak-isospin chirality
  is added as an external $\mathbb{Z}_2$-grading on $\mathcal{H}_F$.

Outcome (i) would be a substantial structural finding; outcome (ii) is the
default and matches Connes–Chamseddine.

---

## §7. Sprint H1 scope statement (eventual sprint)

**Title.** Higgs as inner fluctuation on the GeoVac S³ electroweak triple.

**Inputs.** (a) Real structure $J_{\mathrm{GV}}$ at finite $n_{\max}$
verified (Track 2 closure). (b) R3.5 full-Dirac chirality grading
(`geovac/full_dirac_operator_system.py`, available). (c) The minimal
electroweak finite triple $\mathcal{T}_F = (\mathcal{A}_F, \mathcal{H}_F,
D_F)$ with $\mathcal{A}_F = \mathbb{C} \oplus \mathbb{H}$ and a candidate
$D_F$ from §3.3 (one of A or B).

**Tasks.**
1. Compute $\Omega_D^1(\mathcal{A})$ symbolically at $n_{\max} = 2$. Verify
   that the gauge-1-form sector reproduces Papers 25 / 30 (Cartan torus +
   non-abelian sector).
2. Identify the off-diagonal Higgs scalar field $\Phi$ in
   $\Omega_D^1(\mathcal{A})$. Verify that its zero-mode is a complex doublet
   on each Fock node.
3. Compute the leading-order spectral-action expansion
   $\mathrm{Tr}\, f(D_\omega^2 / \Lambda^2)$ for the fluctuated Dirac
   $D_\omega = D + \omega + J\omega J^{-1}$. Read off the coefficients of
   $|\Phi|^2$, $|\Phi|^4$, $|\nabla \Phi|^2$.
4. Check the Mexican-hat sign pattern (negative $|\Phi|^2$ coefficient,
   positive $|\Phi|^4$ coefficient) at the GeoVac spectral-action ceiling
   $\Lambda_\infty$.
5. Read off the Higgs vev $v_0$ in $S^3$-unit units. Compare to the natural
   GeoVac cutoff $n_{\max}$; note any unit-conversion factor.
6. Test the candidate $D_F$ from §3.3 against the natural negative of §5.

**Outputs.**
- Module `geovac/almost_commutative.py` (this scoping memo provides the
  interface skeleton).
- Computational sprint memo `debug/sprint_h1_higgs_memo.md`.
- Either positive (Higgs construction works on GeoVac, with $D_F$ selected
  from candidate A or B) or negative (one of §5 closes the direction).

**Estimated effort.** 1–2 weeks for steps 1–3 (computational); 2–4 weeks
for steps 4–6 (interpretive). Total ~6 weeks. **Not** the next sprint;
Track 2 must close first.

**Paper venue.** If positive and rich enough, Paper 37 (a self-contained
Higgs-on-GeoVac paper). If positive but thin, Paper 32 §VIII.C addendum
to the four-gap analysis. If negative, Paper 32 §VIII.C closing-statement
addendum that elevates Marcolli–vS-without-Higgs to a GeoVac-internal
result.

---

## §8. Specific recommendations to the PI

These are recommendations to apply post-sprint, **not auto-applied edits**:

- **Update Paper 32 §VIII.B G2 with this scoping result.** Specifically,
  replace "Whether multiplicative inner fluctuations can be defined on a
  non-algebra operator system is open in the NCG literature" (lines
  1636–1638) with a sharper statement: "The natural extension is to
  $\mathcal{A}_{\mathrm{GV}} \otimes (\mathbb{C} \oplus \mathbb{H})$, in
  which inner fluctuations split into a gauge-1-form sector (Papers 25 /
  30 recovered) and a candidate Higgs scalar sector. The off-diagonal
  $D_F$ ingredient that turns Marcolli–vS-without-Higgs into a Higgs
  construction is not yet selected from GeoVac structure; see Sprint H1
  scoping memo `debug/almost_commutative_scoping_memo.md`."
- **Add to Paper 32 §VIII Open Question Q1.** The question "whether
  GeoVac admits a natural rank-$\ge 2$ extension" should be sharpened with
  a forward reference to this memo and to the §3.3 candidate sources.
- **Do NOT update Paper 2 in any way.** This memo says nothing about the
  combination rule $K = \pi(B + F - \Delta)$ and does not change its
  conjectural status (CLAUDE.md §13.5).
- **Open question for the PI: Should Sprint H1 be opened immediately after
  Track 2 closes, or should it wait for Track 1 (R2.5 L4) too?** My
  recommendation: open H1 after Track 2 closure. The construction is
  defined at finite $n_{\max}$ without R2.5; only the continuum
  interpretation requires L4. Running H1 in parallel with L4 maximizes
  throughput.

---

## §9. Files produced

- `debug/almost_commutative_scoping_memo.md` (this file).
- `geovac/almost_commutative.py` (interface stub; ~100 lines, no
  implementation).

## §10. Files referenced

- `papers/synthesis/paper_32_spectral_triple.tex` (Paper 32, especially
  §VIII.B four-gap analysis).
- `papers/observations/paper_30_su2_wilson.tex` (Paper 30, SU(2) Wilson;
  YM-without-Higgs scope statement).
- `papers/synthesis/paper_25_hopf_gauge_structure.tex` (Paper 25, U(1)
  Wilson; Cartan-torus reading).
- `geovac/operator_system.py` (truncated commutative algebra).
- `geovac/spinor_operator_system.py` (Weyl sector).
- `geovac/full_dirac_operator_system.py` (full Dirac with chirality
  grading; pre-stages the off-diagonal Dirac that Sprint H1 would extend).
- `geovac/spectral_triple.py` (existing finite spectral triple data with
  $J$, $\gamma$ already in place).
- `geovac/su2_wilson_gauge.py` (Paper 30 implementation; Marcolli–vS-
  YM-without-Higgs by default).
- `debug/wh1_r35_full_dirac_memo.md` (R3.5 sprint memo; chirality-grading
  structure).
- `debug/sm_unified_gauge_synthesis_memo.md` (origin of Paper 32 §VIII.B
  four-gap analysis).

External references:
- A. Connes & A. Chamseddine, "The spectral action principle," Comm. Math.
  Phys. 186 (1997) 731–750.
- A. Chamseddine, A. Connes, M. Marcolli, "Gravity and the standard model
  with neutrino mixing," Adv. Theor. Math. Phys. 11 (2007) 991–1089.
- A. Connes & M. Marcolli, *Noncommutative Geometry, Quantum Fields and
  Motives*, AMS Colloquium Publications 55, 2008. Especially Ch. 13
  (Higgs as inner fluctuation; finite real structure on $\mathbb{C} \oplus
  \mathbb{H}$).
- W. D. van Suijlekom, *Noncommutative Geometry and Particle Physics*,
  Springer, 2015. Chapters 8 (inner fluctuations) and 9 (almost-commutative
  geometries).
- M. Marcolli & W. D. van Suijlekom, "Gauge networks in noncommutative
  geometry," J. Geom. Phys. 75 (2014), arXiv:1301.3480.
- C. I. Perez-Sanchez, arXiv:2401.03705 (2024), arXiv:2508.17338 (2025) —
  Yang–Mills-without-Higgs correction to the Marcolli–vS continuum limit.

---

**End of scoping memo.**
