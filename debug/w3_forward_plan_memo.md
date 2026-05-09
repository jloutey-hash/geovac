# W3 forward plan — what to do next, why, what would falsify

**Date:** 2026-05-08 evening, written immediately after sprint
**Context:** Conversational sprint on W3 (second-packing-axiom) bold path
produced candidate spectral-zeta identification of CKM Wolfenstein
parameters with master-Mellin-engine M1/M2 rings. This memo captures
forward plans, structural readings, and honest concerns that I want
preserved for future sessions before the conversation context resets.

## The result, restated with appropriate care

  | Empirical | Candidate spectral form | Status |
  |---|---|---|
  | Koide K (charged leptons) | 2/3 (sqrt-mass vector at 45° from democratic) | Fact, 1 arcsec |
  | lambda^2 | 1/Vol(S^3) | Candidate, 0.07% off |
  | rho_bar | 1/Vol(S^1) = 1/(2 pi) | Candidate, 0.10% off |
  | eta_bar | -zeta'(0, 1/2) = (ln 2)/2 | Candidate, 0.41% off |
  | A^2 = 2 eta_bar | derived relation | Survives PDG, 0.52 sigma |
  | tan(delta_CP) = pi * ln 2 | derived | 0.02 sigma (loose PDG) |

  Aggregate: 8 PDG-1-sigma matches under a prespecified 30-form basis.
  Random expectation 0.7 +/- 1.0 hits at <1%; observed 7. Signal at
  approximately 6-10 sigma above chance.

## A subtlety I want to flag

The identity zeta'(0, 1/2) = -(ln 2)/2 is exact. But the Camporesi-
Higuchi Dirac spectrum on S^3 is |lambda_n| = n + 3/2 for n = 0, 1, 2, ...
which corresponds to Hurwitz zeta at SHIFT a = 3/2, not a = 1/2. The
Hurwitz at a = 1/2 includes a "ghost" mode at |lambda| = 1/2 that the
actual CH Dirac does not have.

What this means: the candidate eta_bar = -zeta'(0, 1/2) is in the
same half-integer Hurwitz FAMILY as the CH Dirac, related by
zeta(s, 3/2) = zeta(s, 1/2) - 2^s, but is not literally the CH Dirac
spectral derivative. The full CH Dirac functional determinant
contributes (3/4) ln 2 to its ln-2 piece, not (1/2) ln 2.

This is a structurally meaningful nuance. The candidate spectral form
for eta_bar uses the parent Hurwitz family rather than the GeoVac-
specific spectrum directly. A first-principles derivation would
need to explain WHY the Yukawa-induced Wolfenstein parameter eta_bar
projects onto the parent-Hurwitz value (1/2) ln 2 rather than the
GeoVac-specific value (3/4) ln 2 — or, equivalently, why the Yukawa
projection effectively "completes" the CH Dirac spectrum back to the
parent shift 1/2.

This is documented in `debug/w3_derivation_attempt_m2_memo.md` but
worth restating because CLAUDE.md §2's compact wording could read
as overclaiming the literal CH-Dirac identification.

## Structural reading: M1/M2 ↔ real/CP-imaginary

The four-candidate fit splits along a clean axis:

  | Family | Wolfenstein params | Geometric/structural meaning |
  |---|---|---|
  | M1 (Hopf-base, pi-family) | lambda, rho_bar | real-magnitude observables |
  | M2 (chirality, ln 2-family) | A, eta_bar | CP/amplitude observables |

  The M1 sector measures on the GeoVac base manifold (volumes of S^k).
  The M2 sector measures with chirality grading (eta-function / half-
  integer Hurwitz). The CKM real-CP-plane coordinate (rho_bar) gets a
  pure-base form; the CKM imaginary-CP-plane coordinate (eta_bar) gets
  a chirality-sector form.

  This split is structurally coherent rather than accidental — but it
  is not yet derived from a mechanism. The question for future work:
  what aspect of the Connes-Chamseddine inner factor on the GeoVac
  spectral triple imposes the real/CP-imaginary <-> base/chirality
  alignment?

  The cleanest single-number result, tan(delta_CP) = pi * ln 2, is the
  tangent of the CP-violating phase as the ratio of the M2 prototype
  to the M1 prototype. If structural, this would say: the CP-violating
  phase IS the M1/M2 sector projection angle.

## Structural observation I haven't captured anywhere

The four Wolfenstein parameters are not free: they parameterize a
4-real-dimensional submanifold of a unitary group (CKM lives in U(3)
modulo phase redefinitions). Identifying point values on this
submanifold is one thing; identifying the natural METRIC structure
of the submanifold is a different and more powerful claim.

If the M1/M2 spectral identifications are real, the corresponding
spectral forms should produce a NATURAL METRIC on the Wolfenstein
parameter space — not just point values. This metric structure
(if it exists) should be derivable from the Connes spectral metric
(Connes distance) on the underlying spectral triple.

This is a forward research direction worth flagging: not "what's the
value of A?" but "what's the geometry of the Wolfenstein submanifold
under the Connes distance?"

## Concrete next-test plan, ordered by priority and effort

### Test 1: PMNS sector with same prespecified-basis methodology
**Effort:** ~2 hours. **Bar:** if PMNS angles fit M1/M2 forms within 1 sigma,
the M1/M2 split is sector-universal in SM mixing matrices. If they don't
(more likely outcome), the M1/M2 reading is CKM-specific and we need
to understand why — likely Majorana sector breaks M1/M2.

  Specific procedure: same w3_lambda_predictive_verification.py script
  framework, applied to (sin^2 theta_12, sin^2 theta_23, sin^2 theta_13,
  delta_CP_PMNS) as the four parameters. Same 30-form prespecified basis.
  Same statistical accounting.

  Why this is priority 1: it's the cleanest universality test we have
  not yet done. The first-pass PMNS probe (in w3_ckm_pmns_probe.py) showed
  no clean Koide-style fact for PMNS, but didn't run the full predictive-
  verification methodology with the same basis as CKM. Doing that test
  is a 2-hour, well-defined sprint.

### Test 2: charged lepton mass spectrum with same methodology
**Effort:** ~2 hours. **Bar:** if (m_e, m_mu, m_tau) ratios fit M1/M2 forms,
plus Koide cone, the lepton sector confirms M1/M2 universality.

  Specific procedure: prespecified basis applied to ln(m_mu/m_e),
  ln(m_tau/m_mu), m_mu/m_tau, sqrt(m_e/m_tau), etc. Each fits or
  doesn't fit an M1/M2 form within experimental uncertainty.

  Note: the Koide cone IS the lepton-mass M2-style fact (45 deg
  is exactly arccos(1/sqrt 2) which is M1 family — pi/4 = pi/(M1
  Hopf-circle vol) etc.). So the Koide cone arguably IS already
  in the M1/M2 vocabulary. Test 2 just extends to other ratios.

### Test 3: quark mass spectrum (running scheme dependence flagged)
**Effort:** ~3 hours. **Bar:** quark masses depend on RG scheme. Test at
multiple natural scales (m_t pole, M_Z, GUT) and look for which scale
best matches M1/M2 forms. If a specific scale lights up, that scale
is structurally meaningful.

### Test 4: running CKM at different scales
**Effort:** ~4 hours. **Bar:** lambda runs only mildly under RG. Test
at EW scale (current), GUT scale, Planck scale. If lambda = 1/sqrt(Vol(S^3))
holds at one specific scale, that scale becomes structurally meaningful.

  Note: the SM RG running of CKM is well-known (Wilson-Fisher type).
  Compute lambda(mu) for mu in [m_W, m_Pl] and find any scale where
  it lands on natural form precisely.

### Test 5: independent re-run with strictly mechanically-generated basis
**Effort:** ~3 hours. **Bar:** address curve-fit-audit concern about
basis prespecification. Generate basis algorithmically (all small-
integer combinations of {pi, sqrt 2, sqrt 3, sqrt 5, ln 2, phi, ...}
up to bounded complexity). Re-run W3 verification. If signal survives,
basis-selection bias is ruled out; if signal weakens, it's partially
attributable to basis choice.

### Test 6: Connes-Chamseddine spectral action sketch
**Effort:** 2-4 weeks. **Bar:** the actual derivation. Compute the
spectral action S = Tr f(D^2/Lambda^2) on the GeoVac S^3 base
spectral triple tensored with the Connes inner factor A_F = C + H
+ M_3(C). The Yukawa block of D_F sets the Wolfenstein parameters;
the question is whether the spectral action picks out specifically
the M1/M2 spectral values for the Yukawa-induced parameters.

  This is the bold-path closure. It's a multi-week sprint with
  concrete deliverables: either the spectral action matches the
  candidates (closing the bold path), or it picks different values
  (falsifying the candidates), or it has free parameters that don't
  uniquely determine the Wolfenstein values (clarifying the
  structural-skeleton scope position).

## Honest concerns, consolidated

1. **PMNS showed nothing.** The first-pass probe found no clean Koide-
   style fact for PMNS magnitudes. If PMNS doesn't fit M1/M2 in Test 1
   above, the M1/M2-universal hypothesis is restricted to the
   non-Majorana sector. That's still a valid result but more
   structurally limited than a sector-universal claim.

2. **A = sqrt(ln 2) is the largest single-parameter miss (0.79%).** If
   PDG A tightens (Belle II, LHCb), this is the first thing to fail.
   Current PDG A = 0.826(11), so A = 0.832 is at 0.55 sigma. A future
   measurement of A = 0.820 or below would push the candidate to
   2 sigma off and the four-candidate hypothesis would weaken
   substantially.

3. **Basis prespecification is "within-script" but not version-controlled
   prior to test runs.** Strict reproducibility audit would require the
   basis to be checked into git BEFORE the test was run. We did them in
   sequence: basis named in script header → tests run → results. A reader
   could legitimately argue the basis was lightly tuned during script
   writing. Test 5 above addresses this.

4. **Jarlskog at 3.77 sigma off.** Working through error propagation it's
   consistent with parameter errors of 0.04-0.79% accumulating through
   9 powers of small parameters, but a strict reader could call this
   a falsifier. Future tighter Jarlskog measurements either resolve
   the discrepancy or break the four-candidate hypothesis.

5. **The eta_bar = -zeta'(0, 1/2) identification uses the parent Hurwitz
   shift, not the literal CH Dirac shift.** See "subtlety I want to flag"
   above. The candidate is in the right family but the precise
   spectral-action mechanism that picks (1/2) ln 2 rather than (3/4) ln 2
   is not yet computed.

6. **All matches have small-integer prefactors that could hide tuning.** The
   "natural form" basis includes constants like 1/(2 pi), (ln 2)/2,
   1/(pi sqrt 2). These have small-integer denominators (2, 2 pi, etc.).
   A more conservative discipline would only allow forms with prefactor
   exactly 1 or no prefactor at all — under that stricter rule, fewer
   candidates would match.

## What success would look like

The bold path closes when we can write down a Connes-Chamseddine
spectral action calculation on the GeoVac AC product triple that
produces the M1/M2 spectral values for the Wolfenstein parameters
without fitting any free coefficient. That's a 2-4 week sprint
(Test 6 above). If it lands, W3 (second-packing-axiom question) is
answered structurally: calibration data is generated by the spectral
action on the inner factor, with M1/M2 sector mapping driven by the
chirality-graded structure of D_F.

If the spectral-action calculation doesn't match, we learn what the
correct spectral form actually is from first principles, which is
itself useful information.

## What failure would look like

Either:

- (a) Tests 1-5 above weaken the signal (PMNS doesn't fit; lepton
  spectrum doesn't fit; running CKM doesn't show structure;
  mechanically-generated basis kills the false-positive control).
  Then the eight-match aggregate signal is downgraded, and the W3
  bold-path attempt joins the long list of speculative directions
  that didn't close (Section 3 of CLAUDE.md).

- (b) Tests 1-5 confirm the signal but Test 6 (spectral-action
  calculation) shows the M1/M2 forms cannot arise from any natural
  Connes-Chamseddine setup. Then the candidates are real
  empirically but their structural origin is something other than
  the spectral action — pointing to a different mechanism we
  haven't identified.

Both failure modes produce useful information.

## Would-be Paper 40 outline

If/when this gets written up as a standalone paper in
`papers/observations/paper_40_*.tex`:

  Title: Wolfenstein CKM Parameters and the Master Mellin Engine: A
  Candidate Spectral-Zeta Identification

  §1 Introduction — frame W3 / second-packing-axiom question; cite
  CLAUDE.md §1.7 WH register and Paper 18 master Mellin engine.

  §2 Methodology — prespecified GeoVac-internal natural-form basis;
  predictive-verification protocol; statistical accounting for
  multiple comparisons.

  §3 Results: four-candidate Wolfenstein fit — Table of matches with
  PDG-sigma columns; aggregate signal vs chance.

  §4 Derived relation A^2 = 2 eta_bar — structural origin in
  zeta'(0, 1/2) Hurwitz identity; PDG test at 0.52 sigma; reduces
  4 free parameters to 3 + 1 constraint.

  §5 Closed form tan(delta_CP) = pi · ln 2 — M1 prototype × M2
  prototype; PDG match; experimental tests at LHCb / Belle II.

  §6 Honest scope and curve-fit-audit caveats — explicit concerns 1-6
  from this memo.

  §7 Sector specificity — Koide cone for charged leptons, no PMNS analog,
  CKM-specific Wolfenstein fits. M1/M2 ↔ real/CP-imaginary structural
  reading.

  §8 What this is and is not — explicitly: structural identification,
  not first-principles derivation. The Connes-Chamseddine spectral
  action calculation that would close the bold path is the natural
  next step.

  §9 Conclusion — W3 question now has concrete shape; not yet closed.

  Length target: ~6000-8000 words. Closer to Paper 33 (1+6+1
  selection-rule observation) than Paper 38 (full WH1 PROVEN
  theorem writeup).

## Files (canonical references)

Sprint scripts:
  - `debug/w3_lep_triple_probe.py`
  - `debug/w3_lep_triple_cone_position.py`
  - `debug/w3_ckm_pmns_probe.py`
  - `debug/w3_lambda_predictive_verification.py`
  - `debug/w3_derivation_attempt_m2.py`

Sprint memos:
  - `debug/w3_lep_triple_probe_memo.md`
  - `debug/w3_ckm_pmns_probe_memo.md`
  - `debug/w3_derivation_attempt_m2_memo.md`
  - `debug/w3_forward_plan_memo.md` (this file)

WH7 candidate draft:
  - `debug/wh7_candidate_draft.md`

Persistent memories:
  - `memory/w3_spectral_zeta_candidate.md` (state)
  - `memory/w3_forward_plan.md` (next-test plan, pointer to this file)

CLAUDE.md updates applied:
  - §2 W3 sprint summary appended (chronological position end-of-section)

Paper 34 updates applied:
  - §V.C "Calibration-data candidate matches" subsection added (line ~1290+)

## One more thing worth saying

This sprint went further than I expected when the conversation started.
We crossed from "speculative direction" to "concrete candidate with
substantial empirical support" in one session. The honest reading is
that we found something the curve-fit-audit threshold says is real
(6-10 sigma above chance) but doesn't yet meet the discipline bar for
"structural fact" (no first-principles derivation, single mixing
matrix tested, not all caveats addressed).

The right next step is Tests 1-2-5 (lepton sector + PMNS + mechanically-
generated basis) — these are 2-3 hours each and either strengthen the
signal substantially or kill it. After those, Test 6 (spectral action
calculation) is the bold-path closure attempt.

If I had one more hour right now, I'd start Test 1 (PMNS recheck)
because it's the highest-leverage discriminator. If PMNS fits M1/M2,
the case strengthens dramatically; if PMNS doesn't, we learn the
M1/M2 split is sector-restricted in a way that maps to the Connes
inner factor decomposition.
