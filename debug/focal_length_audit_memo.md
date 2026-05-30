# Focal-length audit of Paper 34's variable axis

**Date:** 2026-05-30
**Question (Josh):** Is there a formula uniting the projection focal points? Audit the variable axis:
(1) which variables are genuine *periods* (compactification scales)? (2) which are dimensionless *couplings*?
(3) for each coupling, is it a *ratio* of two periods?
**Discipline:** audit-claim. Ratios are the numerology-prone part; flagged hard.

## The organizing claim being tested
A focal length = the size of a compact direction; the discrete quantum index is its dual lattice; level spacing
= 1/period. Open the aperture (period -> inf) and spacing -> 0, lattice melts to continuum. If true, the 28
projection variables should collapse to a few genuine periods, with couplings as a separate dimensionless class.

## Classification of all 28 (source: Paper 34 Table tab:projection_axes)

### CLASS P — genuine periods / compactification scales (dimensionful; a discrete index is conjugate)
| Projection | Variable | Compact direction it sizes | Dual index | Spacing ~ 1/period |
|:--|:--|:--|:--|:--|
| Fock conformal | E (via p0=sqrt(-2E)) | momentum 3-sphere radius | n | 1/n^3 yes |
| Bargmann-Segal | hbar*omega | HO phase-space scale | N | yes (rational ladder) |
| Connes-Chamseddine | Lambda (UV cutoff) | inverse lattice spacing | mode index | yes |
| Mol-frame hyperspherical | R | internuclear separation | channel | yes |
| Rest-mass | m | Compton wavelength 1/m = a compact temporal/rest scale | — | sets a scale |
| Observation/temporal | beta = 1/T | imaginary-time circle circumference | Matsubara k | 2pi/beta YES (cleanest) |
| Adiabatic / B-O | R (slow) | slow-coordinate scale | — | yes |
| Coupled-channel | R (or rho) | hyperradius | channel | yes |
| Stereographic | r | radial length (recoordinatization) | — | scale, not periodic |

  PERIOD COUNT: the *genuinely distinct* periods are fewer than the rows. They are all instances of ONE thing —
  "the size of a compact (or to-be-compactified) direction" — wearing 3 physical clothes:
   (i) a SPATIAL period:  Fock E/p0  ~  hyperspherical R  ~  coupled-channel rho  ~  adiabatic R  ~  stereographic r
   (ii) a TEMPORAL period: observation beta  (= the gravity/thermal aperture, Paper 47)
   (iii) a HO/oscillator period: Bargmann hbar*omega  ;  and the field-theory UV period: CC Lambda
  Rest-mass m sits at the spatial/temporal seam (Compton length = rest-energy period).
  => 28 variable-rows collapse to ~3 archetypal periods (spatial size, temporal circumference, oscillator scale).
  This is the compactness thesis (Paper 18 sec III) re-read as "the focal lengths are all compactification scales."

### CLASS C — dimensionless couplings (NOT periods; cannot set a size)
| Projection | Coupling | Role |
|:--|:--|:--|
| Hopf bundle | alpha | electromagnetic coupling across the S^2 base |
| Connes-Chamseddine | alpha (also carries Lambda period) | gauge coupling in spectral action |
| Camporesi-Higuchi | (alpha via (Z alpha)^2 in fine structure) | relativistic coupling |
| Wilson plaquette | beta_W = 1/g^2 | gauge coupling (INVERSE-square of g — note: dimensionless) |
| Two-body Breit | m_l/m_n | mass RATIO (already dimensionless — see Class R) |

  Couplings are a STRUCTURALLY SEPARATE class from periods. They do not size a compact direction; they set how
  strongly the two sides of a boundary talk. alpha is the recurring one (appears in Hopf, CC, CH, Breit-via-alpha^4).

### CLASS L — pure labels / recoordinatizations (no period, no coupling)
Sturmian (reuses Z,n), Wigner 3j, Wigner D, vector-photon, Drake-Swainson (transient K), nuclear charge-density
<r^2>, nuclear magnetization r_Z, nuclear tensor Q_N, Phillips-Kleinman, multipole order L, bipolar (k1,k2),
symmetry tableau lambda, gauge choice, Wick signature index, apparatus rho. These introduce structure-indices or
Layer-2 input data, not focal lengths and not couplings. (Nuclear r^2/r_Z/Q_N are dimensionful but they are
*input data about a fixed nucleus*, not tunable compactification scales — they don't open an aperture.)

## Answers to the three questions

### Q1 — which variables are periods?  ~9 rows, collapsing to ~3 archetypes.
The period archetypes: SPATIAL size (Fock/hyperspherical/adiabatic/stereographic), TEMPORAL circumference
(observation beta = the gravity aperture), OSCILLATOR/UV scale (Bargmann hbar omega, CC Lambda). Every "focal
length" in the framework is one of these three. The unification Josh sensed is real and is the compactness thesis.

### Q2 — which are couplings?  alpha (Hopf, CC, CH) and beta_W=1/g^2 (Wilson), plus the Breit ratio.
Dimensionless, separate class. Confirmed: the dictionary's variable axis splits cleanly into
periods (set WHERE the boundary is) and couplings (set HOW the two sides talk across it). Two roles, one boundary.

### Q3 — is each coupling a RATIO of two periods?  AUDIT RESULT BELOW.

  alpha:
   - TEXTBOOK FACT (not ours, not numerology): alpha = lambda_Compton / a_0 = (1/m) / (1/(m alpha)).
     i.e. alpha = (relativistic length) / (bound-state length) = ratio of the rest-mass period to the Fock period.
     This is a genuine ratio-of-two-periods, and BOTH periods are in Class P above (rest-mass m, Fock bound scale).
   - SHARP QUESTION (falsifiable, NOT yet done): do GeoVac's OWN two projections — the Camporesi-Higuchi/Dirac
     (relativistic) focal length and the Fock (bound) focal length — reproduce alpha = CH-scale / Fock-scale
     STRUCTURALLY, or is restating alpha=lambda_C/a_0 in our vocabulary just textbook redressed?
   - AUDIT VERDICT: PLAUSIBLE-BUT-UNPROVEN. alpha IS a ratio of two periods that both live in our Class P. But
     "alpha = ratio of our two focal lengths" would be CONTENT only if the ratio is FORCED by the projection
     geometry, not inserted by hand. Free-parameter check: alpha=lambda_C/a_0 has ZERO free parameters (it's a
     definition of a_0), so restating it is zero-content. A real result needs the ratio to fall out of the
     S^3-vs-Dirac-S^3 geometry WITHOUT using a_0's definition. Not established. Do NOT claim.
   - CROSS-LINK (real, and the better prize): Paper 2's 1/alpha = pi(B + F - Delta) is ALREADY a
     confinement-boundary decomposition: B = closed/discrete Casimir trace; F = zeta(2) = the SAME continuum-edge
     accumulation sum we computed in the aperture flag (threshold pile-up); Delta = Dirac boundary mode-count.
     So 1/alpha = pi*(discrete + continuum-edge - boundary). alpha's structure is discrete-continuous-boundary
     by construction. This is a structural RE-READING of an existing formula (combination rule stays a numerical
     observation per the hard prohibition); it is NOT a new derivation. But it says alpha-as-boundary-object is
     already in the corpus, independent of the ratio question.

  beta_W = 1/g^2:  NO — this is 1/coupling^2, not a ratio of two periods. Different structure.

  m_l/m_n (Breit):  YES trivially — it is by definition a ratio of two rest-mass periods (Class P / rest-mass).
     This one is clean and uncontroversial: a dimensionless coupling that IS a ratio of two periods.

## Net verdict
- Q1 YES: the focal lengths unify — 28 variable-rows -> ~3 period archetypes (spatial / temporal / oscillator).
  This is the compactness thesis stated as a focal-length-unification. Solid, inherited, not numerology.
- Q2 YES: variable axis splits cleanly into periods (where the boundary is) vs couplings (how sides talk). Clean.
- Q3 MIXED: m_l/m_n is a ratio of periods (trivially, clean). alpha is a ratio of periods in textbook physics
  (lambda_C/a_0) and BOTH periods are in our Class P — but "GeoVac derives alpha as its own focal-length ratio"
  is PLAUSIBLE-BUT-UNPROVEN and zero-content unless the ratio is geometry-forced. beta_W is NOT a ratio.
  => couplings are NOT uniformly ratios-of-periods. The clean general statement is the period unification (Q1),
     not a coupling=ratio law (Q3).

## The genuinely new structural reading (low-risk, capture-worthy)
"Discrete-continuous coupling as a geometric boundary problem" (Josh's phrasing) is correct and precise:
  - PERIODS locate the discrete<->continuum boundary (aperture position).
  - COUPLINGS quantify cross-boundary talk (alpha = EM cross-talk, g = gauge cross-talk).
  - and alpha's OWN defining decomposition 1/alpha = pi(B+F-Delta) is itself discrete + continuum-edge - boundary.
This unifies (a) today's aperture flag, (b) the compactness thesis, (c) the Paper 2 observation, under ONE picture:
the projection variables are coordinates on a confinement boundary — periods say where it is, couplings say how
permeable it is. This is a framing result (proposed-principle strength), not a theorem.

## Honest scope / what is NOT claimed
- No new alpha derivation. Combination rule K stays a numerical observation (hard prohibition respected).
- "alpha = ratio of GeoVac's own focal lengths" is flagged PLAUSIBLE-BUT-UNPROVEN; explicitly not asserted.
- The period unification is the compactness thesis (Paper 18 sec III) re-expressed, not a new theorem.
- The one clean Q3 hit (m_l/m_n) is trivial (a ratio by definition).

## Files
debug/focal_length_audit_memo.md ; reads from papers/group6_precision_observations/paper_34_projection_taxonomy.tex
Table tab:projection_axes ; cross-ref debug/open_confinement_aperture_flag_memo.md
