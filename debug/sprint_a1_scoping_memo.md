# Sprint A1 — scoping memo: inhomogeneous Robertson–Walker extension

Date: 2026-06-03
Scope: 1-day diagnostic-only scoping sprint. NO implementation. NO new closed
forms. NO PSLQ. Output: verdict + plan/obstruction.

## TL;DR

**Verdict: COMPRESSIBLE — sprint-scale (2–3 weeks), with one substantive
caveat.** The naive "does GeoVac have an $a(t)$ parameter?" question has a
clean answer:\ the master Mellin engine sits on $S^3$ alone (no time
direction in the bare substrate), so a literal continuum $a(t)$ requires
extending GeoVac to a 4D substrate. **However, two structurally distinct
"discrete $a(t)$ surrogates" already exist in the corpus:**

1. **Constant warp $a(t) \equiv a_0 \ne 1$ (radius rescaling)** —
   directly accessible via uniform rescaling of the Camporesi-Higuchi
   spectrum $|\lambda_n| = (n + 3/2)/a_0$; this is Sprint G4-3a's
   "constant-warp $r(\rho) = r_h$" already done in Paper 51 §G4-3.
   Trivial extension of M2 (rescales SD coefficients by powers of $a_0$);
   inheritance argument is immediate.

2. **Matsubara-compactified $S^3 \times S^1_\beta$ (Paper 35)** — this
   IS the GeoVac framework's natural time direction, with $\beta$ playing
   a role structurally analogous to $a(t)$ in a thermal-time
   interpretation. The substrate is well-defined, the discrete spectrum is
   computed, and Paper 35 already classifies where $\pi$ enters
   ($2\pi/\beta$ Matsubara modes).

The substantive A1 question (with bounded sprint scope) is:\ **does the
F–M pure-Tate refinement survive on the Matsubara compactified substrate
$S^3 \times S^1_\beta$ at $\beta < \infty$, or does the Matsubara sum
introduce MZV content?** This is sprint-scale.

The full continuum $a(t)$ generic case (general time-derivative spectrum,
$\varepsilon_i \ne 0$ in F–M) is **structurally outside GeoVac's discrete
substrate** — it requires either (i) a continuum 4D substrate (which
GeoVac doesn't natively support) or (ii) a discrete time-dependent
Camporesi-Higuchi spectrum, which would re-engineer the spectral triple
foundation. That is the **NO-GO sub-case** within A1.

A1 therefore splits cleanly into a **sprint-scale POSITIVE half**
(constant warp + Matsubara) and a **NO-GO half** (generic continuum
$a(t)$). The COMPRESSIBLE verdict refers to the sprint-scale half; the
NO-GO half is recorded honestly as a structural boundary.

## Three-option analysis

### Option (a) NO-GO — discrete substrate is intrinsically static

**Argument FOR:** The GeoVac substrate is `Fock(S^3)` — the Fock
projection puts a single 3-manifold $S^3$ at the substrate level, with the
Camporesi-Higuchi Dirac $D_{\mathrm{CH}}$ acting on it. There is no time
direction in the bare substrate. Inserting an $a(t)$ requires either
(i) embedding the $S^3$ as a slice of a 4D Lorentzian/Euclidean spacetime
(continuum), or (ii) putting a discrete time index on the substrate and
making $D_{\mathrm{CH}}$ depend on it. (i) is outside GeoVac's discrete-
substrate principle; (ii) breaks Connes axiom verification at finite
$n_{\max}$ (the spectral triple structure of Paper 32 assumes
time-independent $D$).

**Argument AGAINST (this verdict alone is too restrictive):** The corpus
already has two structurally clean ways to inject a "time-like parameter"
into the discrete substrate without breaking the spectral triple
foundation:

- **Uniform radius rescaling** ($S^3_R$ instead of $S^3_1$, equivalently
  $\lambda = a_0$ in F–M). Paper 51 G1 (Path 1, spectral action on $S^3_R$)
  already verifies this. SD coefficients scale by integer powers of $R$,
  pure-Tate ring preserved.
- **Temporal compactification** ($S^3 \times S^1_\beta$, $\beta$ as
  compactification scale). Paper 35 is the entire paper on this construction.
  Substrate is well-defined; spectrum is `(n+1)^2 + (2\pi k/\beta)^2`;
  classification of where $\pi$ enters is known.

The NO-GO verdict, taken literally, would reject both of these as out of
scope, which contradicts existing corpus content. NO-GO is therefore the
correct verdict only for the **generic continuum $a(t)$ sub-case** with
non-trivial $a'(t), a''(t), \ldots$ (i.e.\ the part of F–M where the
$\varepsilon_i$ are non-zero).

**Verdict:** NO-GO applies to the generic continuum sub-case (the
$\varepsilon_i \ne 0$ slice of F–M), not to A1 as a whole.

### Option (b) COMPRESSIBLE — clean discrete surrogate exists

**Argument FOR:** Two surrogate substrates already exist and are
structurally tractable:

1. **Constant-warp $a(t) \equiv a_0$ (radius rescaling)** — Paper 51 G1
   already does this (`thm:zeta_unit_neg_k` on $S^3_R$). Specifically:
   SD coefficients on $S^3_R$ are $a_k(D^2) = a_k^{\mathrm{unit}}/R^{3-2k}$
   for the leading Dirac contributions, with the volume-normalised
   coefficients sitting in $\pi^{2k} R^{3-2k} \cdot \mathbb{Q}$.

   *Sprint A5 mapping:* this is the $\lambda$-slice $\lambda = a_0,
   \varepsilon_i = 0$ of F–M, which is a closed point in
   $\mathbb{G}_m \times \mathbb{A}^{2n}$. The $\lambda$-direction is the
   non-trivial part of the variation away from $\lambda = 1$; the
   $\varepsilon$-direction stays at the origin. Pure-Tate ring is
   preserved (multiplication by $R^{3-2k} \in \mathbb{Q}$ does not
   introduce MZV content).

2. **Matsubara compactification $S^3 \times S^1_\beta$ at finite $\beta$**
   — Paper 35 substrate. The spectrum is $n(n+2) + 4\pi^2 k^2/\beta^2$
   (scalar). Mellin transform of the heat trace produces a 2D Hurwitz
   structure where the $k$-sum is integer over $\beta$ and the $n$-sum is
   the standard $S^3$ Hurwitz at half-integer shift $3/2$. **The
   substantive A1 question:** does the cross-product structure
   $S^3 \times S^1$ keep the pure-Tate refinement, or does it inject
   $\zeta(\mathrm{odd})$ / MZV content at finite $\beta$?

   Bare structural prediction (not derived, recorded as the diagnostic
   target): the Matsubara sum is a Hurwitz at integer shift, which lives in
   $\MT(\mathbb{Z})$ at depth 1 (per Paper 55 §6 M3 sub-sector 1). At
   integer $s$, this resolves into $\pi^{\mathrm{even}}\cdot\mathbb{Q}$ at
   even $s$ and $\zeta(\mathrm{odd})\cdot\mathbb{Q}$ at odd $s$.
   **So at $\beta < \infty$, MZV content (specifically $\zeta(3)$) may
   enter the cross-product via the temporal Matsubara sum.** This would
   be a genuinely new finding — the pure-Tate refinement is special to
   the spatial-only ($\beta = \infty$) sub-case.

**Argument AGAINST:** Neither surrogate is literally a time-dependent
$a(t)$. Calling them "$a(t)$ extension" requires explanation. The constant-
warp case is a one-parameter rescaling, not a function-valued parameter;
the Matsubara case has a single scale $\beta$, not a function $a(t)$.
Naming convention matters for the literature audience (F–M's "$a(t)$" is
function-valued).

**Verdict:** COMPRESSIBLE applies cleanly to (i) the constant-warp
sub-case (already done in Paper 51 G1, modular needs only a paragraph in
Paper 55) and (ii) the Matsubara sub-case (substantive sprint, 2–3 weeks).
The COMPRESSIBLE verdict for A1-as-a-whole is conditional on naming the
sprint correctly: it's "A1 = compact time-direction surrogates of $a(t)$"
not "A1 = full F–M generic continuum".

### Option (c) MULTI-WEEK STANDS — original 2–4 week estimate confirmed

**Argument FOR:** The original v3.45.3 estimate ("2–4 weeks") was
implicitly assuming the F–M generic continuum sub-case with $a(t)$
function-valued. That sub-case genuinely requires either (i) building a
continuum 4D substrate analog of GeoVac with $a(t)$ baked into the metric,
or (ii) constructing a discrete time-dependent Camporesi-Higuchi spectrum
with proven spectral-triple axioms at every time slice. Both are
multi-month at minimum.

**Argument AGAINST:** The compression argument identifies two specific
sub-cases (constant warp, Matsubara) that the original estimate did not
distinguish from the generic case. With the sub-case split made explicit,
the generic part falls to NO-GO and the sub-cases compress to sprint-
scale. The original 2–4 week estimate was the un-stratified estimate.

**Verdict:** MULTI-WEEK STANDS for the generic continuum sub-case
(matches the original estimate). The stratified COMPRESSIBLE+NO-GO reading
is more honest about what's actually accessible.

## Structural argument for the chosen verdict (COMPRESSIBLE)

The case-exhaustion theorem (Paper 32 §VIII) and the master Mellin engine
on Camporesi-Higuchi $S^3$ are stated for the **static, single-radius
substrate**. They do not literally extend off this point. However, the
two natural one-parameter deformations of the substrate that GeoVac
**does** support (constant rescaling and temporal compactification) are
already in the corpus, and the F–M classification adapts naturally to
both:

| Deformation | F–M slice | Substrate | Status | Pure-Tate? |
|:--|:--|:--|:--|:--|
| **None** ($\lambda = 1$, $\varepsilon = 0$) | closed point | $S^3$ | Paper 55 §4 | YES |
| **Constant warp** ($\lambda = a_0$, $\varepsilon = 0$) | 1D slice in $\mathbb{G}_m$ | $S^3_{a_0}$ | Paper 51 G1 | YES (trivially) |
| **Matsubara** (no F–M analog, GeoVac-internal) | n/a — new substrate | $S^3 \times S^1_\beta$ | Paper 35 | OPEN (this sprint) |
| **Generic $a(t)$** ($\lambda$ varying with $t$, $\varepsilon \ne 0$) | full F–M family | continuum 4D | NO-GO | n/a |

The **constant-warp** row is structurally a corollary; it requires only a
paragraph in Paper 55 §7 cross-referencing Paper 51 G1.

The **Matsubara** row is the substantive A1 question:\ a thermal-time
generalisation of the static F–M classification. The Matsubara
contribution is structurally distinct from F–M's $\varepsilon$-derivatives:\
F–M's time-derivative content is purely **classical** (Lagrangian / metric
derivatives), while Matsubara is **thermal** (modular flow on KMS state).
They probe genuinely different sectors of the spectral structure. The A1
question on Matsubara therefore has a direct connection to:

- Paper 35 "time as projection" (single-paper coverage of the temporal
  compactification mechanism)
- Paper 32 §VIII case-exhaustion theorem (does the M2 sub-mechanism extend
  to thermal observables?)
- Paper 42 modular Hamiltonian (Tomita-Takesaki modular flow)
- Paper 43 Krein wedge (Lorentzian extension)
- Paper 50 F-theorem on $S^3$ + boundary (which already engages M1+M3
  cross-products via $\log 2$, $\zeta(3)/\pi^2$)

This places A1-Matsubara at a structural cross-roads in the corpus.

The **generic $a(t)$** row is a known structural boundary:\ GeoVac is a
3-manifold-substrate framework. Promoting to a 4-manifold substrate (with
$a(t)$ in the metric) is the kind of multi-month infrastructure build that
Paper 47 §6 already does for the Lorentzian non-compact case, and that the
G4-3 discrete-substrate program is currently underway for.

## Sprint plan (the COMPRESSIBLE half — Matsubara)

**Working title:** Sprint A1-Matsubara — temporal compactification of the
F–M mixed-Tate classification

**Sprint scale:** 2–3 weeks (diagnostic + computation + write-up)

**Phase 1: setup (3–4 days)**

- Read F–M paper through §7.7 with thermal product geometries in mind
- Identify what F–M says about $S^1$-fibered substrates if any (they do
  treat thermal $\R \times S^3$, with $a(t)$ — but is there a $\beta$
  variant?)
- Survey adjacent literature on thermal spectral actions:\ Chamseddine-
  Connes thermal spectral action, Dabrowski et al. on thermal modular
  flow, Hekkelman-McDonald (Riemannian spectral triple sequences) if any.
- Concurrent-work CLEAR check at math.OA / hep-th level via
  WebSearch+WebFetch.

**Phase 2: computation (1 week)**

- Compute the 2D heat trace
  $K_{S^3 \times S^1_\beta}(t) = K_{S^3}(t) \cdot K_{S^1_\beta}(t)$
  (factorisation already in Paper 51 G2, but at the Mellin-transform
  level not at the period level)
- Compute SD coefficients $a_k(D^2_{\mathrm{total}})$ at integer $k$ for
  $k = 0, 1, 2, 3, 4$. Check whether they stay in $\bigoplus_k
  \pi^{2k}\cdot\mathbb{Q}$ (pure-Tate) or introduce $\zeta(3), \zeta(5)$
  through the temporal Matsubara at odd integer arguments.
- Specifically: the temporal $\zeta_{S^1_\beta}(s) = 2\zeta(s)/\beta^s$
  introduces classical $\zeta$ values at integer $s$. At odd integer $s$,
  $\zeta(2k+1)$ enters. Whether this $\zeta$ content survives in the
  4D SD coefficient extraction or cancels structurally is the diagnostic
  target.
- Negative-result-friendly framing: the diagnostic should produce a clean
  "yes Matsubara injects $\zeta(3)$" or "no it cancels" answer either way.

**Phase 3: write-up and corpus integration (4–5 days)**

- Canonical sprint memo `debug/sprint_a1_matsubara_memo.md`
- Paper 55 §7 inhomogeneous extension paragraph upgrade
- Paper 35 cross-reference to the periods-classification
- Paper 32 §VIII case-exhaustion theorem either extends to thermal
  observables (positive) or has its scope clarified to spatial-only static
  $S^3$ (negative)
- Paper 51 G2 cross-reference: the Stefan-Boltzmann constant $\pi^2/90$
  classification (in Paper 35) gets re-framed in periods language

**Decision gate at end of Phase 2:**

- **POSITIVE-pure-Tate-survives:** Matsubara structurally preserves
  pure-Tate at integer $s$ (most likely via the high-temperature expansion
  being a sum of even-weight modes). Paper 55 §7 picks up a clean
  sub-sprint.
- **POSITIVE-MZV-enters:** Matsubara genuinely introduces $\zeta(3)$ at
  finite $\beta$. This is more interesting structurally — it would mean
  the temporal compactification step is the first GeoVac mechanism that
  pushes M2 out of pure-Tate into generic mixed-Tate. Genuinely new finding.
- **NEGATIVE-no-structural-engagement:** Matsubara sums trivialise in the
  4D SD coefficient extraction (the temporal $\zeta$ values cancel against
  the Stefan-Boltzmann normalisation). Sprint produces a clean negative
  with a clear corpus paragraph.

All three outcomes are publishable as a sprint memo + Paper 55 §7
paragraph.

## Sprint plan (the trivial half — constant warp)

**Scope:** Single afternoon, no sprint needed. Just add a paragraph to
Paper 55 §7 cross-referencing Paper 51 G1 (`thm:zeta_unit_neg_k`),
explaining that the constant-warp sub-case is the $\lambda$-slice of
F–M and trivially preserves the pure-Tate ring. Could be bundled into the
Matsubara sprint write-up.

## NO-GO half (full continuum generic $a(t)$)

**Recommendation:** Add a sentence to Paper 55 §7 inhomogeneous extension
section recording the structural obstruction:

> The full F–M generic continuum case with non-trivial time derivatives
> $a'(t), a''(t), \ldots$ (the $\varepsilon_i \ne 0$ slice in the notation
> of Theorem 4.2's specialisation) requires either a 4D continuum
> substrate or a time-dependent Camporesi–Higuchi spectral triple. Both
> are outside the scope of the discrete-substrate principle that defines
> GeoVac. The closest GeoVac analog is the discrete-warped-substrate
> program of Paper 51 §G4-3, which is currently a multi-month follow-on.

No active sprint on this row.

## Files used

- `papers/group3_foundations/paper_55_periods_of_geovac.tex` (full, §4
  thm:m2_mixed_tate proof sketch, §7 open questions inhomogeneous
  extension paragraph)
- `papers/group6_precision_observations/paper_35_time_as_projection.tex`
  (Matsubara spectrum, $\pi$-injection at exactly $2\pi k/\beta$, Casimir
  energy on $S^3$ classification)
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` §VIII
  (case-exhaustion theorem; static-$S^3$ specialisation paragraph at
  line 1739–1745)
- `papers/group5_qed_gauge/paper_51_gravity_arc.tex` §G1 (`thm:zeta_unit_neg_k`),
  §G2 (thermal product $S^3 \times S^1_\beta$ heat-trace factorisation),
  §G4-3 (constant-warp discrete substrate factorisation at lines
  1036–1068)
- `papers/group3_foundations/paper_24_bargmann_segal.tex` (HO rigidity
  theorem — different substrate, dimension-parity sharpening already
  treated in Paper 55 §7.3)
- `debug/sprint_mixed_tate_test_memo.md` (v3.45.3 sprint; open question
  #1 = inhomogeneous extension, this memo's parent)
- `debug/sprint_a5_fm_full_read_memo.md` (Sprint A5 explicit F–M
  structure; clarifies that the static $a(t) \equiv 1$ is the
  $\lambda = 1, \varepsilon_i = 0$ closed point, and that off-locus
  variation is parameterised by $\lambda \in \mathbb{G}_m$ and
  $\varepsilon \in \mathbb{A}^{2n}$)

## Honest scope

This memo is diagnostic-only. It does NOT:

- Run PSLQ on Matsubara SD coefficients
- Derive closed forms for $K_{S^3 \times S^1_\beta}(t)$ SD expansion
- Verify the Matsubara sub-mechanism actually does or does not preserve
  pure-Tate (that's the Phase-2 sprint deliverable)
- Edit any paper

It DOES:

- Identify the two-fold stratification (constant warp / Matsubara) +
  NO-GO generic
- Confirm that the COMPRESSIBLE verdict holds for the stratified sub-cases
- Confirm that the MULTI-WEEK estimate was correct for the un-stratified
  generic case
- Provide a sprint plan if PI authorises the Matsubara sub-sprint
- Provide a single-paragraph drop-in for Paper 55 §7 covering the
  constant-warp sub-case + NO-GO recording for the generic case (no sprint
  needed for those)

## Robustness of the COMPRESSIBLE verdict

The main risk to the COMPRESSIBLE verdict is that the Matsubara sub-case
might not be a genuine "$a(t)$" extension in the F–M sense:\ the F–M
parameter $a(t)$ is a function-valued metric coefficient
(R–W cosmology), while $\beta$ is a single compactification scale (thermal
field theory). The two are physically distinct. The COMPRESSIBLE verdict
holds **mathematically** (the period-theoretic question is sensible for
either substrate), but **physically** the connection between A1 (inhomo R-W
extension) and Matsubara is loose.

A cleaner framing for the sprint:\ "A1' — thermal-time compactification
of the F-M pure-Tate refinement" rather than "A1 — inhomo R-W extension".
The thermal-time framing maps directly onto the Paper 35 + Paper 42 +
Paper 43 + Paper 50 thermal-time arc, and the period-theoretic question
fits naturally as the missing classification piece of that arc.

If the PI wants the **physical** A1 (inhomo R-W with function-valued
$a(t)$), then the verdict is MULTI-WEEK STANDS at minimum, more likely
multi-month, and the structural prerequisite is the G4-3 discrete-warped-
substrate program of Paper 51.
