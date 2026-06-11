# Sprint W3-CKM-PMNS — second probe of the bold path

**Date:** 2026-05-08
**Author:** PM session, conversational sprint
**Context:** Universality follow-up to the lepton-Koide probe. We found a clean
sector-specific 45-deg cone in the lepton mass triple. Question for this probe:
do the quark mixing matrix (CKM) and the neutrino mixing matrix (PMNS) admit
their own packing-style structural constraints?

## TL;DR

1. **Quark sector: TWO candidate single-data-point structural facts surfaced.**
   Both are within experimental uncertainty and have GeoVac-internal natural form.
   Curve-fit-audit discipline applies: each is one data point against a
   parameterized form, so they are **candidates, not facts**.

   - **Wolfenstein lambda = 1/sqrt(Vol(S^3)) = 1/(pi sqrt 2) to 0.035 percent.**
     PDG: lambda = 0.22500(67). Candidate: 1/(pi sqrt 2) = 0.225079.
     Volume of the unit 3-sphere is the framework's natural measure on the
     Fock-projected manifold. If real, this would say: the leading CKM mixing
     scale is set by the natural amplitude on the Fock manifold.

   - **Ratio of log-spacings (|V_us| -> |V_cb| -> |V_ub|) = 1/ln(2) = log_2(e)
     to 0.004 percent.** Computed: 1.4427475. Candidate: 1.4426950. PDG
     uncertainty on the ratio is ~2 percent so this is well within 1 sigma.
     If real, this is an unusual relation: V_us / V_cb = (V_cb / V_ub)^{ln 2}.
     Information-theoretic flavor (ln 2 is the bit-vs-nat conversion).

   - Cabibbo angle theta_C is consistent with 13.000 deg exact to 0.023 percent
     (PDG gives 13.04 +/- 0.05 deg, so this is at the 1 sigma edge). Known
     speculation in the literature; not new.

2. **Neutrino sector (PMNS): no clean structural fact at this probe depth.**
   Tri-bimaximal benchmark is the closest "simple" pattern (sin^2 of three
   angles equal to 1/3, 1/2, 0); it is broken at the 2-5 percent level. No
   clean Koide-cone fact, no PSLQ hit, no GeoVac-internal near-miss.

3. **Sector asymmetry has more structure than I expected.** Each of the three
   sectors (lepton mass, quark mixing, neutrino mixing) shows DIFFERENT
   candidate structural behavior:
   - Lepton mass triple: ONE clean cone constraint (45 deg, 1 arcsec).
   - Quark mixing: TWO candidate scaling-style constraints (lambda and log
     spacing ratio), each at single-data-point level.
   - Neutrino mixing: NO constraint visible; closest pattern (TB) is broken.

## Pulling on these threads disciplined

The Koide-cone for leptons is on solid empirical ground: ONE constraint
predicts the third mass from the first two to 0.006 percent, well below
experimental uncertainty -- that is "fact" level.

The two quark candidates flagged here are NOT yet at that level. Each is a
single-data-point match against a parameterized form. By the curve-fit-audit
memo discipline (`docs/curve_fit_audit_memo.md`), the only escape from
"intriguing coincidence" status is **predictive verification**: the candidate
must imply another observable that we can test against independent data.

Predictive tests that would elevate each candidate:

  Candidate A: lambda = 1/sqrt(Vol(S^3)) = 1/(pi sqrt 2)
    Implies lambda^2 = 1/Vol(S^3) = 1/(2 pi^2).
    Numerically: lambda^2 = 0.05063 (PDG-derived), 1/(2 pi^2) = 0.05066.
    Diff = 0.06 percent. Consistent.
    Implies lambda^4 = 1/(2 pi^2)^2 = 1/(4 pi^4).
      Numerically: lambda^4 = 0.00256, 1/(4 pi^4) = 0.00256. Consistent.
    These don't add INDEPENDENT data points (lambda^k follows trivially from
    lambda). To get an independent test, the candidate would need to imply
    something OUTSIDE the Wolfenstein lambda hierarchy -- e.g., the A
    parameter, the rho/eta values, or a related quantity in another sector.

  Candidate B: ratio of log-spacings = 1/ln(2)
    Implies V_us = V_cb^a where a = ... [solving] ... = ln(R_2)/ln(R_1)+1.
    This is consistent with the data by construction but does not predict any
    independent observable.
    Independent test would need: same ratio relation in PMNS hierarchy?
      PMNS magnitudes |U_e2|, |U_mu3|, |U_e3| = 0.548, 0.731, 0.148.
      log(0.548/0.731) = -0.288; log(0.731/0.148) = 1.595.
      ratio = 1.595 / (-0.288) = -5.54 (negative because order is wrong).
      Take abs: 5.54. Far from 1/ln(2). NOT universal.
    Similar test on lepton mass log spacings: ratio = 0.529, not 1/ln(2).
    So this candidate is QUARK-SECTOR-SPECIFIC. Either a quark-side
    structural fact, or a coincidence.

## Sector-mapping reading

Putting both probes together produces a richer picture than either alone:

  | Sector         | Algebra factor       | Candidate pattern                         | Cleanness  |
  |----------------|----------------------|-------------------------------------------|------------|
  | Lepton mass    | (C + H)^3            | Koide 45-deg cone                         | 1 arcsec   |
  | Quark mixing   | M_3(C) on flavor     | lambda = 1/sqrt(Vol S^3); log-ratio 1/ln2 | within 1 sigma |
  | Neutrino mixing| (C + H)^3 + Maj seas | broken tri-bimaximal                       | 2-5 percent off |

This tabular reading is suggestive: the cleanest pattern is in the sector
where the inner algebra is simplest (C + H)^3 for charged leptons. The
neutrino sector is also (C + H)^3 plus the Majorana extension, and shows
structure (tri-bimaximal-broken) but not a clean fact. The quark sector lives
in M_3(C) and shows TWO candidate facts of different shape.

## What would a second packing axiom look like, after this probe?

Pre-probe expectation: a single packing rule on the inner factor that
generates calibration data.

After this probe: more nuanced. Each sector (lepton mass, quark mixing,
neutrino mixing) appears to have its OWN candidate constraint of distinct
form. The "second packing axiom" question splits into multiple sub-questions:

  1. Why does the lepton sqrt-mass vector lie at exactly 45 deg from
     democratic? (Koide derivation from (C + H)^3 algebra)
  2. Why is the Wolfenstein lambda set by 1/sqrt(Vol S^3)? (CKM derivation
     from M_3(C) and the natural manifold measure)
  3. Why does the log-spacing ratio of the CKM hierarchy land near 1/ln(2)?
     (CKM derivation, possibly information-theoretic)
  4. Why is PMNS near tri-bimaximal but broken? (PMNS derivation from
     (C + H)^3 + Majorana, possibly group-theoretic)

These are FOUR open questions at four different levels of cleanness. The
lepton case is the strongest claim and the cleanest target. The quark
candidates are structurally interesting because they cross-reference GeoVac
internal numbers (Vol S^3, ln 2) but are not yet predictively verified.

## Honest scope

Six sub-probes ran. Two single-data-point matches surfaced for CKM. PMNS
showed no candidate at this probe depth. The Koide cone for leptons remains
the single cleanest packing-style fact in the calibration data.

The pattern across all three probes (lepton + CKM + PMNS) is consistent
with the structural-skeleton-scope reading sharpened in May 2026: the
framework cleanly determines selection rules and scaling laws, while
calibration data shows scattered hints of structural origin without any
clean derivation. The bold path here -- chasing a "second packing axiom" --
is producing CANDIDATE structural facts in each sector but cannot at this
probe depth distinguish between (a) genuine structural facts awaiting
derivation, and (b) numerical coincidences in a parameter set that is in
fact externally input.

## Files

- `debug/w3_ckm_pmns_probe.py` -- main probe (8 sub-probes)
- `debug/data/w3_ckm_pmns_probe.json` -- raw data

## Recommended next moves

1. **Predictive verification of Candidate A (lambda = 1/sqrt(Vol S^3)).**
   Identify the natural cross-sector implications. If lambda is set by
   sqrt(Vol S^3) as a graph-amplitude scale, what does this imply for the
   A parameter, for the lepton-sector mixing, or for the leading CP-violation
   scale? Each independent prediction is a new data point.

2. **Decomposer-style scoping of the (C + H)^3 -> Koide 45 deg derivation.**
   Independent of the CKM candidates. Does the algebra of three copies of
   C + H admit a natural inner-product constraint forcing the sqrt-mass
   vector to lie at 45 deg from the diagonal? This is the lepton-sector
   second packing axiom hypothesis from the previous memo, unchanged.

3. **CKM at GUT scale.** The Wolfenstein parameters run with energy. Do the
   candidate facts hold at any specific running scale (EW, GUT, Planck)?
   If candidate A only holds at one specific scale, that scale becomes
   structurally significant.

4. **Apply same probe to charged-fermion mass anarchy (mass ratios across
   different fermion types).** Are there cross-sector coincidences between
   lepton masses and quark masses that we missed by treating them
   sector-by-sector?
