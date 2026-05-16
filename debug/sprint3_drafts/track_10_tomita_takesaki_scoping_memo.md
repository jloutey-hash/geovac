# Track 10 Scoping Memo: Tomita–Takesaki Modular Flow as a §III Entry?

**Sprint:** Sprint 3 of Paper 34 dictionary-completion arc
**Date:** 2026-05-15
**Mode:** Diagnostic-only — NO §III subsection draft
**Author:** Track 10 worker (PM-dispatched scoping)
**Status when written:** Sprint 2 has just (today) promoted Wick rotation /
signature change to §III.27 and explicitly recast the residual content
as the "operator-system-level Lorentzian extension" entry in §VIII.

---

## 1. Question and structural test

Tomita–Takesaki (TT) theory associates to any von Neumann algebra
$\mathcal{M} \subset B(\mathcal{H})$ with a cyclic and separating
vector $\Omega$ a *modular operator* $\Delta$ (positive, self-adjoint,
generally unbounded) and a *modular conjugation* $J$ (antiunitary,
$J^2 = 1$). The pair satisfies $J \mathcal{M} J = \mathcal{M}'$
(commutant) and the modular flow $\sigma_t(A) = \Delta^{it} A
\Delta^{-it}$ is a one-parameter group of $*$-automorphisms of
$\mathcal{M}$ characterized by the KMS condition at inverse
temperature $\beta = 1$ (in suitable normalization). When the state
$\omega(A) = \langle \Omega, A \Omega\rangle$ is a KMS state at
inverse temperature $\beta$ for a physical Hamiltonian $H$, the
modular flow coincides with $\sigma_t(A) = e^{i\beta t H} A
e^{-i\beta t H}$. At the Bisognano–Wichmann wedge with the Minkowski
vacuum, $\sigma_{2\pi}(A) = \Lambda_B(2\pi) A \Lambda_B(-2\pi)$ where
$\Lambda_B(\eta)$ is the Lorentz-boost subgroup preserving the wedge,
and analyticity of the modular flow gives $\sigma_{i\cdot 2\pi} =$
identity.

The framework already invokes TT implicitly: §III.27 names "the same
KMS-Tomita-Takesaki algebra" as the bridge that ties Hawking, BW, and
Unruh together (line 2276), and the §III.15 footnote names the
modular flow as the Lorentzian-side reading of the M1 $2\pi$.

**Structural test for a separate §III entry.** A projection earns its
own §III slot if and only if it satisfies *both* of:

1. It introduces variables, dimensions, or transcendental content
   that no existing §III entry already covers (the V/D/P-tagging
   independence test, per §IV three-axis dictionary);
2. It produces empirical anchors — quantitative observables — that
   the existing chain does not already produce as Layer-2 outputs.

The default is *absorption*: if TT is the algebraic machinery that
*underlies* §III.15 + §III.27 without producing distinct V/D/P content
or distinct anchors, it stays inside the existing entries (with TT
named explicitly in their footnotes and the operator-system-level
extension named in §VIII).

The default is *requires-extension*: if TT produces structurally
distinct content but only at the operator-system level that §III.27
itself flags as the open multi-month NCG-mathematics target, then TT
sits naturally in the §VIII open-extension entry that already names
exactly this gap.

The default is *own-entry*: only if TT produces empirical or
structural content that neither §III.15 nor §III.27 captures at the
metric-functional level *and* that the framework already accesses
without the extension.

---

## 2. Comparison to §III.15 (observation / temporal-window)

§III.15 compactifies Euclidean time at finite period $\beta$ and
injects $2\pi/\beta$ per Matsubara mode. Its V/D/P signature:
$\beta$ (or $T = 1/\beta$); $[T] = [E]^{-1}$;
$2\pi \cdot \mathbb{Q}$ per Matsubara mode at the spectrum level and
$\pi^{2k} \cdot \mathbb{Q}$ from $\zeta_R(2k)$ in spectral integrals.
Its anchors: Casimir at $T \neq 0$, Stefan–Boltzmann $\pi^2/90$ on
$S^3 \times S^1_\beta$, all finite-temperature partition functions.
The footnote already names the Bisognano–Wichmann reading and
explicitly references the modular flow as the Lorentzian-side dual
interpretation, with $\sigma_{i\cdot 2\pi} = $ identity flagged as the
principal falsifier.

**What TT would distinctly add over §III.15.** At the metric-functional
level: nothing in V/D/P. The $2\pi$ in $\sigma_{i\cdot 2\pi}$ is the
*same* $2\pi$ that §III.15 already injects via Matsubara at the
generator level (M1 Hopf-base measure / $\mathrm{Vol}(S^1)$). TT
does not introduce a new variable beyond $\beta$ (or, for a wedge
algebra without a parameter, the wedge geometry itself supplies
$\beta = 2\pi/a$ as in Unruh), and the dimension and transcendental
ring are unchanged. TT *does* add structural objects that §III.15
does not — the modular operator $\Delta$, the conjugation $J$, the
KMS condition itself as a state-side characterization — but these
are operator-algebra objects on the Lorentzian-wedge algebra, not
*spectral content* of the truncated metric spectral triple.

The cleanest reading is that §III.15 currently spans two distinct
physical content classes via a single mechanism: thermal-state physics
(finite $T$, Matsubara) and wedge-algebra physics (KMS, modular
flow). The Bisognano–Wichmann footnote of §III.15 already
acknowledges this. The question for Sprint 3 is whether splitting
the entry helps the dictionary's clarity. Argument against splitting:
both readings are bound to the same $\beta = 2\pi/\kappa$ identity
and share a falsifier; the V/D/P tags are identical; the empirical
content reduces to the same M1 sub-mechanism. Argument for splitting:
§III.27 (Wick rotation) was just promoted out of the same footnote
for analogous reasons (distinct *mechanism*: analytic continuation
vs. compactification), so the precedent for fine-grained mechanism
distinctions inside the M1 family is fresh. The difference is that
§III.27 promotes by reading direction (Euclidean ↔ Lorentzian as a
signature toggle); a TT promotion would split *within* the
Lorentzian-side reading.

---

## 3. Comparison to §III.27 (Wick rotation / signature change)

§III.27 just promoted Sprint 2 today. Its V/D/P signature: signature
index (binary discrete choice between Euclidean and Lorentzian),
dimension-preserving, $2\pi \cdot \mathbb{Q}$ via the same M1
mechanism. Its anchors: the four-witness Wick-rotation theorem
(Hawking + Sewell + BW + Unruh), Sprint TD Track 4's reproduction of
$T_H = 1/(8\pi M)$, Sprint Unruh-pendant's reproduction of
$T_U = a/(2\pi)$. Its Honest-scope explicitly says the reading is
*structural correspondence*, not literal identification: the framework
computes a Euclidean spectral-action functional and the Wick-rotation
prescription reads it as a Lorentzian observable.

Critically, §III.27 is already the *Lorentzian-reading projection*.
Lines 2244–2255 say it cleanly: "at the metric functional level the
framework's Euclidean output IS the Lorentzian modular-flow /
boost-orbit period, with no extension required. At the
operator-system level (Tomita–Takesaki modular automorphism acting
on the truncated algebra at finite $n_{\max}$, …) the framework
currently sees the same $2\pi$ but cannot internally close the
period identity; that closure is the named open extension."

This is the cleanest possible diagnostic statement: §III.27 *already
identifies TT as the operator-system-level machinery whose framework-
internal realization is the open extension*. TT is not a parallel
projection sitting next to §III.27 — it is the algebraic structure
underneath the *Lorentzian side* of §III.27's correspondence, with
§III.15 supplying the corresponding Riemannian side. The V/D/P axes
of TT are identical to §III.27's; the falsifier is identical
($\sigma_{i \cdot 2\pi} = $ identity on $\mathcal{T}_{n_{\max}}$);
the four anchors are the same four-witness theorem.

A separate "TT projection" entry would essentially duplicate §III.27
at the operator-system level, with the same generator, same
falsifier, same anchors, and same V/D/P signature. This would
recursively split one entry into two whose only distinction is "the
algebraic vocabulary" — exactly the splinter pattern that the §13.8
splinter-file prohibition is the discipline-level analog of.

---

## 4. Empirical content test

The framework currently produces four witness anchors that
*structurally engage* TT: Hawking $T_H$, Sewell's KMS identification
of the Hartle–Hawking thermal state, BW's modular-flow identification
of the Lorentz boost, and Unruh $T_U$. The TT operator $\Delta$
shows up *named* in each via the Lorentzian-side bridge, but the
framework computes none of these via Tomita–Takesaki theory directly:
all four are computed via the master Mellin engine M1 sub-mechanism
($\mathcal{M}[\mathrm{Tr}(D^0 \cdot e^{-tD^2})]$ on a $\beta$-period
$S^1$ factor) and *interpreted* via the published Wick-rotation
chain.

**Would a TT entry produce new empirical anchors?** The natural
candidates are:

- **Connes' modular spectral triples.** Connes (1994) constructed
  spectral triples where the Dirac operator is the modular operator
  $\Delta^{1/2}$ of a type III von Neumann algebra. This would be a
  structurally distinct construction: a *new* spectral triple, not a
  new projection acting on the existing GeoVac triple. It sits at
  the Paper 32 spectral-triple-foundation level, not the Paper 34
  Layer-2 projection level. Promoting TT to a §III entry would
  conflate these levels.

- **KMS states for non-equilibrium QFTs.** TT applies to any
  cyclic-separating $(\mathcal{M}, \Omega)$ and the modular flow
  defines a thermal time. The Connes–Rovelli "thermal time
  hypothesis" (1994) is the program of taking *physical time itself*
  as the modular flow of a state. This is a striking foundational
  reading but produces no quantitative observable that the framework
  doesn't already compute via M1 — the $2\pi/\kappa$ identities are
  the same numerical content.

- **Type III von Neumann algebras / free-field wedge algebras.** The
  algebra of free QFT observables localized in a Rindler wedge is
  type III$_1$ (Bisognano–Wichmann + Driessler 1977; Buchholz et al.
  1990s). This is a deep structural feature with no immediate
  empirical anchor that the GeoVac framework computes. It is also
  the *target* of the operator-system-level extension named in
  §VIII: building the framework's truncated algebras at finite
  $n_{\max}$ in a way that genuinely captures type III character in
  the GH limit is part of what would close the BW correspondence
  literally.

- **Relative modular operators / Araki relative entropy.** The
  Araki entropy $S(\omega \| \phi) = -\langle \Omega, \log
  \Delta_{\phi|\omega} \Omega\rangle$ generalizes the von Neumann
  relative entropy to general $W^*$ algebras. This *does* produce a
  structurally distinct class of observables — but the catalogue
  already names this in §VIII's "candidate state-side dictionary
  enumeration" entry (relative entropy is item 3 in the explicit
  list, lines 6263–6264). It belongs to the state-side dictionary
  that Sprint TD Track 5's PSLQ negative opened, not to the
  spectral-side projection family that TT modular flow would join.

- **KMS-decoupling / split-property / Reeh–Schlieder.** These are
  structural features of QFT algebras, not framework observables.

In every case the empirical content either (i) reduces to M1
$2\pi/\kappa$ identities the framework already computes
(absorbed by §III.27), or (ii) requires the operator-system-level
construction that §III.27 itself names as the multi-month open
extension, or (iii) belongs to the state-side dictionary opened by
§III.28 and elaborated in §VIII rather than to the modular-flow
mechanism specifically.

**No distinct anchor is currently engaged by the framework that TT
alone would catalogue.** The two anchors that TT-promotion would
need to bring — a *quantitative* GeoVac result whose closed form
contains a TT object (modular operator eigenvalue at finite
$n_{\max}$, Araki relative entropy between two truncated CH states,
Connes–Rovelli thermal time on a $T_{S^3} \otimes T_{S^1_\beta}$
construction) — would all live downstream of the operator-system-
level extension that is currently open. Until that extension lands,
TT engages the framework only at the level of the bridge it supplies
between §III.15 and §III.27.

---

## 5. Verdict

**Absorbed into §III.15 + §III.27, with explicit deepening flagged
in §VIII (operator-system-level Lorentzian extension entry).** No
new §III entry. No requires-extension status as a new candidate.

**Justification.** TT is the algebraic structure underneath the
Bisognano–Wichmann reading that §III.15 and §III.27 *both* invoke:
§III.15's footnote names $\sigma_{i \cdot 2\pi}$ as the principal
falsifier, §III.27's Honest scope names TT explicitly ("Tomita–
Takesaki modular automorphism acting on the truncated algebra at
finite $n_{\max}$") as the operator-system-level content of the open
extension, and §III.27 line 2276 names "the same KMS-Tomita–Takesaki
algebra" as the bridge tying all four Wick-rotation witnesses
together. The V/D/P signature of TT modular flow is identical to
§III.27's. The falsifier is identical. The four anchors are the
same four-witness theorem. Splitting TT into a separate §III slot
would duplicate §III.27 with no distinct empirical, V/D/P, or
falsifier content — the splinter pattern §13.8 is the discipline-
level analog of.

**What to add to existing entries (not in scope for this scoping
memo, but flagged for any subsequent dictionary-completion
sprint).**

- §III.15 footnote already names BW and the modular-flow reading.
  A one-clause addition naming "Tomita–Takesaki modular flow" as
  the operator-algebra structure that underwrites the BW reading
  would make the relationship explicit; the substance is already
  present.
- §III.27 Honest scope already names "Tomita–Takesaki modular
  automorphism" explicitly. No edit needed.
- The §VIII operator-system-level Lorentzian-extension entry
  already names the modular Hamiltonian $K$ on $\mathcal{T}_{n_{\max}}$
  as the literal-identification verifier. A one-sentence note that
  "this is equivalently the framework-internal realization of
  Tomita–Takesaki modular flow on the truncated wedge algebra"
  would make the TT framing of the open extension explicit.
- If a future GeoVac sprint computes a quantitative observable
  whose closed form contains a TT object (modular operator
  spectrum at finite $n_{\max}$, Araki relative entropy between
  truncated states, Connes–Rovelli thermal time), *that* would
  trigger a revisit of this verdict and possibly motivate a
  TT-modular-flow §III slot. As of 2026-05-15, no such sprint
  is in flight, and the state-side complement opened by §III.28
  is the more natural near-term direction for non-spectral-side
  dictionary growth (per the §VIII state-side enumeration entry).

**Bottom line for the PM:** TT is correctly classified by the
current paper as "the algebraic vocabulary of the
Bisognano–Wichmann reading whose framework-internal closure is
the named open extension," not as an independent projection.
Promoting it to §III would duplicate §III.27 without adding
content; deferring it to §VIII (where it is already implicit in
the operator-system-level extension entry) keeps the dictionary
clean and the falsifier structure single-faced.
