# WH7 candidate draft — FALSIFIED, NOT PROMOTED

**Date:** 2026-05-08 evening (drafted) → 2026-05-08 night (falsified)
**Status:** **FALSIFIED.** NOT promoted to CLAUDE.md §1.7. Retained as institutional record of what was almost added and why it was rejected.
**Source sprint:** W3 second-packing-axiom bold path (debug/w3_*.{py,md})

## Falsification summary (2026-05-08 night)

Three follow-up tracks dispatched in parallel the same evening as the
original sprint produced three independent negatives:

- **PMNS recheck:** zero matches at <2% in the same 30-form basis;
  signal at −0.68σ below random chance. PMNS sits closer to TBM
  pure rationals (1/3, 1/2, 0) than M1/M2 forms. W3 universality
  falsified. (`debug/w3_pmns_recheck_memo.md`)
- **Lepton mass spectrum:** 2 within-1σ matches across 9000 tests,
  both encoding the same Koide cone fact. Signal at +0.85σ,
  indistinguishable from random. W3-as-mass-spectrum falsified.
  (`debug/w3_lepton_mass_recheck_memo.md`)
- **Mechanical basis sweep:** mechanically-generated 21,448-form
  basis frozen before testing returns aggregate z = −0.52 across
  the four Wolfenstein parameters. Original M1↔real / M2↔CP-imaginary
  structural reading does not survive (η̄'s top mechanical match is
  M1+ALG, not M2). Selection-bias verdict: fully attributable.
  (`debug/w3_mechanical_basis_memo.md`)

**Three independent tests, three negatives. The W3 candidate is dead.**

The original "5–10σ above chance" calculation was statistically
correct for that basis size; the error was the implicit assumption
that the 30 hand-curated forms were drawn without knowledge of the
data. Mechanical generation removes that assumption, and the signal
vanishes.

This file is retained because it documents what was nearly added to
the WH register and the reasoning that was used to justify it. Future
sessions should treat it as a falsified candidate, not an open
proposal.

CLAUDE.md §2 W3 entry has been updated to mark the sprint falsified.
Paper 34 §V.C has been reframed as a documented negative. Memory
files updated.

---

## Original draft text (preserved unchanged below for the record)

---

## Proposed addition to §1.7

Suggested wording, following the format of existing WHs:

---

**WH7 — Calibration data lives in the master Mellin engine M1/M2 rings of the spectral triple.**

The Standard Model parameters that the framework does not autonomously
generate (W3 in the multi-focal-wall taxonomy) are observed to fit clean
spectral-zeta forms in the master Mellin engine vocabulary, with M1
(Hopf-base measure, $\pi \cdot \mathbb{Q}$) governing real-magnitude
observables and M2 (chirality-shifted half-integer Hurwitz zeta,
$\ln 2 \cdot \mathbb{Q}$) governing CP/amplitude observables. CKM
Wolfenstein four-candidate fit (Sprint W3, 2026-05-08): $\lambda^2 =
1/\mathrm{Vol}(S^3) = 1/(2\pi^2)$ and $\bar\rho = 1/\mathrm{Vol}(S^1) =
1/(2\pi)$ in the M1 family; $\bar\eta = -\zeta'(0, 1/2) = (\ln 2)/2$
and $A^2 = -2\zeta'(0, 1/2) = \ln 2$ in the M2 family. All four match
PDG within 1 sigma. The derived structural relation $A^2 = 2 \bar\eta$
follows and survives PDG at 0.52 sigma, reducing 4 free Wolfenstein
parameters to 3 + 1 constraint. Charged lepton mass triple satisfies
the Koide $45^\circ$ cone constraint (1 arcsec) with sector-specific
break to quark Koide angles. CP-violating phase satisfies the closed
form $\tan(\delta_{CP}) = \pi \cdot \ln 2$ at 0.02 PDG sigma.

*Falsifier:* tighter CKM measurements at LHCb / Belle II breaking the
$A^2 = 2\bar\eta$ relation; or a single Wolfenstein parameter (or
analogous lepton-sector parameter) shown to fit a clean form OUTSIDE
the M1 $\cup$ M2 ring at sub-PDG-sigma precision; or any first-principles
spectral-action calculation on the GeoVac AC factor proving the M1/M2
identifications cannot arise from the Connes--Chamseddine machinery.

*Status:* **CANDIDATE** as of single-session sprint 2026-05-08. The
empirical signal (8 PDG-1$\sigma$ matches, ~6--10$\sigma$ above chance
under prespecified basis, plus one new derived predictive relation that
also survives PDG, plus one clean closed form) is too strong to dismiss
but not yet at the level of "fact" because (a) only one mixing matrix
(CKM) tested with this discipline; (b) PMNS mixing showed no clean
analog at first-pass probe; (c) lepton mass spectrum not yet tested
with the same prespecified-basis methodology; (d) no first-principles
derivation from the Connes-Chamseddine spectral action attempted; (e)
basis-selection bias cannot be fully ruled out without independent
verification. The structural reading: the W3 question (what generates
calibration data) now has a concrete shape rather than being the
open hand-wave of prior framing.

---

## Why this WH and not just a Paper-34 catalogue addition

The catalogue addition (Paper 34 §V/V.B/V.C) records the empirical
matches. WH7 names the STRUCTURAL CLAIM behind those matches: that
calibration data has a structural origin in the spectral triple, located
specifically in the master Mellin engine ring vocabulary. Without WH7,
the eight matches sit as scattered coincidence rows in a catalogue.
With WH7, they are organized as evidence for one specific structural
hypothesis with named falsifiers.

This is exactly the role the WH register is designed for per §1.7's
governance description: bold internal claims that organize evidence,
distinct from the cautious paper-rhetoric.

## Why CANDIDATE status, not stronger

Existing WHs at MODERATE or STRONG status have multiple independent lines
of evidence (e.g., WH1 had Marcolli-vS literature lineage + R2 prop=2
matching Connes-vS Toeplitz S^1 + multiple PSLQ-grade computations).

WH7 has, after one session: 8 matches in one mixing matrix, one derived
relation, one closed form. That's a coherent body of evidence but it's
all from CKM. PMNS analog and lepton-mass analog are open. So CANDIDATE
is the appropriate starting status; promotion to MODERATE would require
the lepton sector to also fit the M1/M2 split.

## What I'd flag before you commit this

1. The basis was prespecified but not infinitely rigid. Independent
   re-run with a strictly mechanically-generated basis would tighten
   the false-positive control claim.

2. The Jarlskog discrepancy is at 3.77 sigma in the four-candidate
   verification (Step 4 of debug/w3_lambda_predictive_verification.py).
   Working through error propagation it's consistent with individual
   parameter errors of 0.04-0.79% accumulating through 9 powers of
   small parameters — but a strict reader could call this a falsifier.

3. The "natural form" matches for A and eta_bar both involve ln 2.
   Two parameters in one family is small-N. The M1/M2 split needs
   more data to be confirmed as a structural pattern rather than a
   2+2 coincidence.

4. The structural identification eta_bar = -zeta'(0, 1/2) is exact
   math, but our PDG-side eta_bar agrees with it only to 0.16 sigma
   (i.e., agrees within experimental uncertainty but is not 0.0 sigma
   exact). If future precision moves PDG eta_bar by even one sigma
   away, the identification is in trouble.

## Recommended PI action

Either:
(a) Add WH7 as drafted to CLAUDE.md §1.7 as CANDIDATE. The framework's
    standing position on W3 is updated; future sprints have a named
    structural claim to test against.
(b) Hold WH7 on the bench (this memo) until the lepton-mass sector
    has been tested with the same methodology. If the lepton sector
    confirms the M1/M2 pattern, promote at that point.
(c) Reject — treat the W3 results as Paper 34 catalogue additions only,
    without a §1.7 structural claim.

My recommendation: **(a)** with the explicit CANDIDATE status. The
register is designed for bold internal claims; the empirical signal
qualifies; the falsifiers are named. Sleeping on it doesn't change
that calculus.

But this is a §1.7 edit and per §13.5 access control PMs cannot make
it autonomously. Your call.
