# Yukawa patterns under renormalization-group running: literature scan

Date: 2026-10-29. Read-only scan, input to a later pre-registered GeoVac test
(Paper 32 sec VIII.C). No GeoVac quantity was computed. The only arithmetic is
ratios of published tabulated masses, marked DERIVED.

Verification key: **V** = checked against the INSPIRE record or arXiv abstract
(or, for PDG, the review text itself); **V-abs** = the abstract confirms the claim
but the number quoted is from the body and was not re-read; **SECONDARY** = taken
from a search-engine summary or a later paper, not the primary source.

Downloads: one file, `debug/lit_scan/dl_yukawa_rg/pdg_quark_masses/` (PDG 2024
Quark Masses review, PDF plus pdftotext output).

## Summary table

| Pattern | Set by / survives / destroyed by running | Loop order | Precision today (inputs, scale) | Primary citation | Verified? |
|---|---|---|---|---|---|
| 1a. Top Yukawa at the Pendleton-Ross ratio y_t^2/g_3^2 = 2/9 | SET (IR-attractive fixed ratio), but the flow is too slow and the SM top is not there | 1-loop | Fails: PR predicted m_t ~ 135 GeV; measured pole 172.57 GeV | Pendleton & Ross, PLB 98 (1981) 291 | V |
| 1b. Top at the Hill quasi-fixed point (y_t large at the GUT scale) | SET (quasi-fixed point: low-energy value insensitive to the UV value) | 1-loop (2-loop refinements) | SM: fails, Hill gave ~240 GeV. MSSM: m_t ~ (190-210 GeV) sin(beta), so tan(beta) ~ 1.5-2, a region disfavoured by LEP and in tension with m_h = 125 GeV | Hill, PRD 24 (1981) 691; Carena & Wagner hep-ph/9407208 | V (Hill); V-abs / SECONDARY (MSSM numbers) |
| 2. b-tau unification y_b = y_tau at M_GUT | DESTROYED in the SM (QCD running separates them); can SURVIVE in SUSY at small or large tan(beta), only with O(10-50%) threshold corrections | 1-loop QCD gives the main effect; 2-loop and thresholds decide | SM, 2e16 GeV: m_b/m_tau = 0.59 (DERIVED from XZZ). MSSM, no thresholds: 0.77 (tan beta = 10), 0.79 (tan beta = 50) | Buras, Ellis, Gaillard, Nanopoulos, NPB 135 (1978) 66; Xing, Zhang, Zhou, PRD 77 (2008) 113016 | V (BEGN abstract: m_b ~ 4.8-5.6 GeV); V (XZZ tables) |
| 2'. t-b-tau unification | Possible only in SUSY SO(10) at tan(beta) ~ 50; depends on the tan(beta)-enhanced threshold delta_b, up to ~50% | 2-loop plus 1-loop thresholds | Not a clean test: the outcome depends on the sign of mu and on the SUSY spectrum | Hall, Rattazzi, Sarid, PRD 50 (1994) 7048 | V-abs |
| 3. Koide Q = 2/3 | Holds for POLE masses. For MS-bar running masses it is DESTROYED by ~0.2% (QED pole/running conversion); the running-mass Q is then nearly scale-INVARIANT in the SM | Effect is 1-loop QED; Q's flow above M_Z is O(y_tau^2) | Pole: Q = 0.6666645(51), -0.43 sigma from 2/3 (DERIVED, PDG m_tau = 1776.93 +/- 0.09 MeV). Running: Q ~ 0.6679-0.6688 from M_Z up to 2e16 GeV | Koide, Lett. Nuovo Cim. 34 (1982) 201; PRD 28 (1983) 252; Xing & Zhang, PLB 635 (2006) 107 | V |
| 4a. Gatto-Sartori-Tonin, sin(theta_C) ~ sqrt(m_d/m_s) | SURVIVES: both sides are RG-invariant at leading order | 1-loop (QCD cancels exactly in mass ratios; Yukawa terms suppressed by light Yukawas or lambda^2) | sqrt(m_d/m_s) = 0.2249 (DERIVED, PDG 2024, Nf = 2+1+1, MS-bar at 2 GeV) vs \|V_us\| ~ 0.224, agreement ~0.3% before the m_s and m_d errors (~2%) | Gatto, Sartori, Tonin, PLB 28 (1968) 128 | V (title and abstract); relation form SECONDARY |
| 4b. Froggatt-Nielsen powers of lambda | Free initial data at the flavour scale; integer powers SURVIVE running, but O(1) coefficients absorb running | n/a (it is a texture, not a flow) | Not falsifiable at better than O(1) | Froggatt & Nielsen, NPB 147 (1979) 277 | V |
| 4c. Same-sector light mass ratios m_u/m_c, m_d/m_s, \|V_us\|, \|V_ud\|, \|V_cs\|, \|V_tb\| | SURVIVE (RG-invariant to leading order) | 1-loop | QCD: exact cancellation in quark mass ratios (PDG). Yukawa running of the Cabibbo angle is "tiny" | Olechowski & Pokorski, PLB 257 (1991) 388; Balzereit, Mannel, Pluemper, EPJC 9 (1999) 197; PDG 2024 Quark Masses review | V |
| 4d. \|V_cb\|, \|V_ub\|, \|V_td\|, \|V_ts\|, J, and ratios involving m_t or m_b | RUN (they increase sizably from M_Z to M_GUT), so they are not RG-protected | 1-loop, driven by y_t^2 | n/a | Balzereit et al. (above) | V-abs |
| 4e. Exact one-loop SM flavour invariants (traces built from Y_u, Y_d) | SURVIVE exactly by construction | 1-loop exact | Hold by construction; content lies in which combinations exist | Harrison, Krishnan, Scott, PRD 82 (2010) 096004 | V-abs |
| 5a. Reduction of couplings (y_t as a function of g_3) | SET by imposing RG invariance of the relation | 1-loop (2-loop later) | Failed: m_t = 81 GeV, m_h = 61 GeV (three generations, ~10-15%) | Kubo, Sibold, Zimmermann, NPB 259 (1985) 331 | V |
| 5b. Multiple-point criticality | SET by a boundary condition (degenerate vacuum at M_Pl) imposed on the flow | 2-loop | Predicted m_t = 173 +/- 5 GeV, m_h = 135 +/- 9 GeV (m_h is off by ~1.1 sigma of its own error) | Froggatt & Nielsen, PLB 368 (1996) 96 | V |
| 5c. Asymptotic safety (gravity) | SET (UV fixed point makes y_t, y_b, g_Y IR-predictive) | FRG truncation | m_t ~ 171 GeV (pole); a top-bottom mass splitting of ~170 GeV for some gravity couplings; m_h = 126 GeV (Shaposhnikov-Wetterich) | Eichhorn & Held, PLB 777 (2018) 217 [1707.01107]; PRL 121 (2018) 151302; Shaposhnikov & Wetterich, PLB 683 (2010) 196 | V-abs |
| 5d. SUSY-GUT and strong-unification fixed points | SET (IR fixed points of the unified theory) | 1-loop plus SUSY thresholds | e.g. m_t = 175 +/- 2 GeV (SU(5) fixed point, Ferrandis) | Lanzagorta & Ross, PLB 349 (1995) 319; Ghilencea, Lanzagorta, Ross, PLB 415 (1997) 253; Ferrandis, PRD 68 (2003) 015001 | V |

## 1. Top Yukawa near 1

**Statement.** The one-loop SM equation is
16 pi^2 dy_t/dt = y_t (9/2 y_t^2 - 8 g_3^2 - 9/4 g^2 - 17/12 g'^2).
Two different "fixed-point" readings exist, and they should not be conflated.

- **Pendleton-Ross (true IR fixed ratio).** Drop the electroweak terms and use
  b_3 = 7. Then R = y_t^2/g_3^2 has an IR-attractive fixed point at R* = 2/9.
  This follows from the beta functions, and its derivation is standard. The
  Pendleton & Ross abstract (V) predicts "a top quark mass of about 135 GeV"
  and a Higgs mass of about 72 GeV. Both are wrong. The measured pole mass
  (PDG) is ~172.6 GeV, i.e. y_t(m_t) ~ 0.94, well above the fixed ratio. The
  approach to R* is logarithmically slow, so between M_GUT and M_Z the real
  coupling does not reach it.
- **Hill (quasi-fixed point).** Take any large y_t(M_GUT). The low-energy value
  then collapses onto an upper envelope (the "quasi-fixed point") that is
  independent of the UV value. Hill's abstract (V) gives a fixed-point top
  mass of about 240 GeV in the SM. A later secondary summary quotes ~230 GeV
  pole mass as excluded (SECONDARY: Yeghiyan, hep-ph/9807411 body; its
  abstract does not say this). **SM verdict: the measured top is NOT at the
  quasi-fixed point.** Its y_t(M_GUT) ~ 0.4-0.5 is perturbative and is free
  initial data (SECONDARY: the magnitude follows from the XZZ table, m_t(2e16
  GeV) = 74 GeV).
- **MSSM.** m_t^QFP ~ (190-210 GeV) sin(beta) (Carena & Wagner; the number
  is SECONDARY). Matching m_t then forces tan(beta) ~ 1.5-2 (Abel & Allanach,
  PLB 415 (1997) 371, V-abs). That region gives m_h below ~105 GeV unless the
  stops are multi-TeV. Nevzorov et al. (Phys. Atom. Nucl. 65 (2002) 285, V):
  "moderate values of tan(beta) are almost excluded" by LEP II. After the
  125 GeV Higgs this region needs very heavy stops (SECONDARY: no LHC-era
  primary source was verified).

**Loop order.** The structure is one-loop. Two-loop terms and thresholds
shift the numbers by ~5%.

## 2. b-tau and t-b-tau unification

**Statement.** In SU(5) the b and the tau sit in the same 5-bar, so y_b = y_tau
at M_GUT. BEGN 1978 (V) obtained m_b ~ 4.8-5.6 GeV from m_tau. That was
historically successful, but it was a 6-flavour, one-loop, pre-precision result.

**Running.** QCD drives m_b/m_tau up by a factor ~2.3-3 from M_GUT down to m_b.
The relation is therefore not RG-invariant: it is a boundary condition whose
low-energy image depends on alpha_s, on the spectrum and on the thresholds.
- **SM (destroyed).** XZZ Table IV (V; m_H = 140 GeV input) gives
  m_b(2e16) = 1.00 +/- 0.04 GeV and m_tau(2e16) = 1.6856 GeV, so
  m_b/m_tau = 0.59 (DERIVED). Huang & Zhou (PRD 103 (2021) 016010, V, 3-loop
  Yukawa RGEs, 2-loop matching) give 0.689 at 1e12 GeV (DERIVED from their
  Table 2/3). Croon et al. (arXiv:1903.04977, V) say the SM Yukawas "are far
  from unification." A deficit of ~40% at M_GUT is not a threshold-sized gap.
- **MSSM (survives conditionally).** Without SUSY thresholds, XZZ give 0.77
  (tan beta = 10) and 0.79 (tan beta = 50) (DERIVED). Exact unification needs
  either tan(beta) ~ 1.5 (where y_t near its QFP drags y_b down) or
  tan(beta) ~ 40-60 with a tan(beta)-enhanced gluino/chargino correction
  delta_b of the right sign. Ross & Serna (PLB 664 (2008) 97, V) write: "not
  all the classic GUT mass relationships agree quantitatively at small to
  moderate tan beta", with the large-tan(beta) fix requiring gluino and wino
  masses of opposite sign. Croon et al. note tension between low-tan(beta)
  b-tau and m_h = 125 GeV.
- **t-b-tau** (Hall, Rattazzi, Sarid, V-abs): needs tan(beta) ~ 50, and the
  m_b prediction is controlled by delta_b, which can reach 50% (Croon et al.).
  It is not a sharp test of the flow. It tests the SUSY spectrum.

## 3. Koide relation

**Statement.** Q = (m_e + m_mu + m_tau)/(sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2 = 2/3.
Koide 1982 (Lett. Nuovo Cim. 34, 201, V: a subquark model) and Koide 1983
(PRD 28, 252, V) predicted "m_tau = 1.7866 GeV" from m_e and m_mu. That is the
1983 version; the modern 2/3 form gives 1776.969 MeV.

**Precision today.** Using PDG m_tau = 1776.93 +/- 0.09 MeV (PDG live average,
V) with CODATA m_e and m_mu gives Q = 0.6666645 +/- 0.0000051, -0.43 sigma
from 2/3 (DERIVED). The tau mass that Q = 2/3 implies is 1776.969 MeV.

**Running.** The relation **holds for pole masses only.** Xing & Zhang (V):
for running masses at M_Z, Q deviates by ~0.2%, and the running Q is "nearly
stable" from M_Z up to ~1e14 GeV. This agrees with the XZZ tables: Q = 0.6679
(SM, 2e16), 0.6683 (MSSM tan beta = 10), 0.6724 (MSSM tan beta = 50) (all
DERIVED). Huang-Zhou give 0.6679 at 1e12 GeV (DERIVED).

The mechanism has two parts. The MS-bar QED anomalous dimension is
flavour-universal, so mass ratios (and Q) do not run in pure QED above all
thresholds. What breaks Q = 2/3 is the pole-to-running conversion
m_pole = m(mu)[1 + (alpha/pi)(1 + (3/4) ln(mu^2/m^2))]. Its log depends on
the mass, so a relation exact for pole masses cannot be exact for running
masses at a common scale. Sumino (JHEP 0905:075, V) states that the QED
radiative correction violates the pole-mass relation unless it is cancelled.
He supplies a U(3) family-gauge model whose log cancels the QED log (the
size of the correction is in the body, not re-read). **No RG-invariant
reformulation that keeps Q = 2/3 for running masses was found.** The running
masses carry a scale-stable Q ~ 0.668 instead. That value is an observed
quantity without a theory.

**Read for GeoVac.** Koide is a pole-mass (on-shell, IR-endpoint) relation, so
it is not a statement about any MS-bar scale. If it is structural, it must be
attached to the on-shell projection, not to the flow.

## 4. Hierarchy and CKM pattern

- **Froggatt-Nielsen** (V): O(1) couplings times powers of a small ratio of a
  flavon vev to a mass scale, set by U(1) charges. The abstract states that
  the mixing angle between two quark types is "roughly the square root of
  their mass ratio." The integer exponents are discrete and pass through
  running intact. The O(1) prefactors are free, so FN is a texture with free
  data, not a prediction.
- **GST** (Gatto, Sartori, Tonin 1968, V: the Cabibbo angle from vanishing
  quadratic divergences of weak self-masses). Modern form
  |V_us| ~ sqrt(m_d/m_s) (form SECONDARY; the correction -e^{i phi} sqrt(m_u/m_c)
  is from Fritzsch-type textures and was not verified). The PDG 2024 review (V)
  states quark-mass ratios are "scheme and scale independent up to possible QED
  corrections". It gives m_s = 92.93 MeV and m_d = 4.70 MeV (Nf = 2+1+1, MS-bar,
  2 GeV), so sqrt(m_d/m_s) = 0.2249 (DERIVED) against |V_us| ~ 0.224
  (SECONDARY: the PDG CKM value was not re-fetched).
  Above M_Z the Yukawa running of m_d/m_s and of the Cabibbo angle is
  suppressed: Olechowski-Pokorski (V) derive the one-loop equations for |V_ij|
  and the Yukawa ratios, and Balzereit et al. (V) report that "the evolution of
  the Cabibbo angle is tiny". XZZ's m_d/m_s is 0.052 at 2e16 GeV (SM) against
  0.0506 at 2 GeV (DERIVED; the large errors make this consistent).
  **GST SURVIVES running at one loop**, because both sides are 1-2 generation
  quantities insensitive to y_t.
- **Which ratios run.** At one loop with y_t dominant, the mixing elements
  involving the third generation (|V_cb|, |V_ub|, |V_td|, |V_ts|) and J run with
  exp(+c integral of y_t^2). Same-generation-pair ratios (m_u/m_c, m_d/m_s) and
  the 2x2 Cabibbo block are invariant to leading order. |V_ub/V_cb| is invariant
  at leading order (SECONDARY: the standard result, not re-verified in the
  primary). Harrison-Krishnan-Scott (V-abs) build exact one-loop SM invariants
  from flavour-symmetric traces of the Yukawa matrices, but note they "do not
  appear to be relevant ... in the MSSM".

## 5. Other RG fixed-point or boundary-condition predictions

The numbers are in the table (rows 5a-5d). Notes:
- The multiple-point principle gets the top right but its Higgs is ~1 sigma
  high. The vacuum-stability status was not verified here.
- In Eichhorn & Held the gravity couplings are free over "a broad range", so
  their top mass depends on which trajectory is chosen.
- **"Lanyon-type"**: no author of that name was found. The likely referent is
  Lanzagorta & Ross (IR fixed points in unified theories) and
  Ghilencea-Lanzagorta-Ross "strong unification" (predicts m_t, m_b/m_tau and
  tan beta), both V.

## Reading for a forced/free test (not a verdict; for pre-registration design)

1. **Every fixed-point success needs one extra free boundary condition**
   (MSSM tan beta, the MPP vacuum degeneracy, the AS gravity trajectory). The
   flow alone fixes none of them. The SM top is not at its IR fixed point, so
   in the SM y_t is free initial data that the flow carries but does not set.
2. **Scale-invariant patterns** (GST; m_u/m_c, m_d/m_s; the Cabibbo block;
   the HKS exact invariants) are the only ones whose truth does not depend on
   the scale. Any claimed "forced" value for them would be falsifiable at every
   scale at once.
3. **Koide** is scale-independent for running masses (Q ~ 0.668) but equals 2/3
   only on shell. It cleanly separates "projection-forced" (pole) from
   "flow-forced" (MS-bar).
4. **b-tau** is a boundary condition destroyed by running in the SM. It has no
   clean SM test, and its SUSY test is controlled by threshold data.

Open items not done here: the LHC-era MSSM top-QFP exclusion source; the PDG
|V_us| re-fetch; the primary source for the |V_ub/V_cb| invariance; the size of
the QED correction to Koide in Sumino's body text.
