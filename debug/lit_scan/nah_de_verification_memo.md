# NaH X^1Sigma+ dissociation energy: primary-source verification

**Date:** 2026-09-06
**Trigger:** Paper 58 Table III footnote flags `D_e = 1.961 eV` for NaH as
"inherited from the framework's earlier compilation ... primary source not
yet verified"; not tabulated on the NIST WebBook page carrying
`r_e = 1.8874 Angstrom = 3.5667 a_0` (currently cited as `\cite{huberherzberg}`).

## Ranked sources checked

1. **Huang, Lu, Whang, Chang, Tsai, "Dissociation energy of the ground state
   of NaH," J. Chem. Phys. 133, 044301 (2010).**
   DOI: 10.1063/1.3458914 (confirmed via Crossref bibliographic search and
   cross-checked against PubMed PMID 20687644, which resolves to the same
   DOI). PRIMARY, experimental (stimulated-emission-pumping +
   fluorescence-depletion spectroscopy on 114 rovibrational levels,
   9 <= v'' <= 21, 1 <= J'' <= 14; highest level ~40 cm^-1 from the limit).
   Reported: **D_e = 15815 +/- 5 cm^-1**.
   Verified via WebSearch snippets of the AIP abstract page and ADS record
   (https://ui.adsabs.harvard.edu/abs/2010JChPh.133d4301H/abstract); full
   text is paywalled (AIP returned HTTP 403 to WebFetch), so the abstract
   wording is taken from search-engine-indexed snippets, not a directly
   fetched PDF/HTML body. Value could not be independently re-derived from
   raw data in this pass.

2. **Chu, Huang, Whang, Tsai, "Observation of double-well potential of NaH
   C^1Sigma+ state: Deriving the dissociation energy of its ground state,"
   J. Chem. Phys. 148, 114301 (2018).**
   DOI: 10.1063/1.5020827 (confirmed via Crossref; PubMed PMID 29566521
   located but full text/abstract not retrievable past a cookie wall in
   this pass). PRIMARY, experimental, same group, refines source 1 using
   the C-state double-well analysis (v = 6-42, inner+outer well +
   near-dissociation region, pulsed OODR fluorescence depletion).
   Reported: **D_e(X) = 15807.87 +/- 5 cm^-1** (and D_e(C) = 6595.10 +/- 5
   cm^-1 for the excited C state, not relevant to Paper 58).
   This is the **most recent primary determination found** (no post-2018
   revision located; searches through 2026 turned up nothing newer for the
   X-state D_e).

3. **Chu, He, Lin, Li, Whang, Tsai, "Spectroscopic determination of the
   ground-state dissociation energy and isotopic shift of NaD," J. Chem.
   Phys. 147, 024301 (2017).** DOI: 10.1063/1.4991036. Isotopologue (NaD),
   not NaH directly, but internally consistent: reports
   D_e(NaD) = 15822 +/- 5 cm^-1 and isotope shift
   delta D_e = D_e(NaH) - D_e(NaD) = -7 cm^-1, i.e. an implied
   D_e(NaH) = 15815 cm^-1 -- matching source 1 (this 2017 paper predates
   the 2018 refinement in source 2).

4. **Huber, K. P. and Herzberg, G., *Constants of Diatomic Molecules*
   (Van Nostrand Reinhold, New York, 1979)**, as mirrored on the NIST
   Chemistry WebBook diatomic-constants page for NaH
   (https://webbook.nist.gov/cgi/cbook.cgi?ID=C7646697&Mask=1000, data
   compiled through March 1977) and the NIST CCCBDB experimental-data page
   (https://cccbdb.nist.gov/exp2x.asp?casno=7646697). **Does NOT carry a
   confident experimental D_e/D_0 for NaH X^1Sigma+.** The WebBook page
   only offers a footnote: "extrapolation of the ground-state vibrational
   levels suggests 2.1 eV," alongside computed (not measured) estimates of
   1.92 eV and 1.88 eV from early theoretical work. This *confirms* the
   task's premise: 1.961 eV is genuinely absent from the Huber-Herzberg /
   NIST WebBook page that supplies `r_e`, so citing `\cite{huberherzberg}`
   for the D_e number (as Paper 58 currently does implicitly via the same
   footnote/table) is not defensible -- Huber-Herzberg predates the
   definitive 2010/2018 measurements by decades and only offers a rough
   extrapolated bound, not a value.

## Unit conversion (1 cm^-1 = 1.239841984e-4 eV, CODATA hc/e)

| Source | D_e (cm^-1) | D_e (eV) |
|---|---|---|
| Huang et al. 2010 | 15815 +/- 5 | 1.9608 -> rounds to **1.961** |
| Chu et al. 2018 (refined) | 15807.87 +/- 5 | 1.9599 -> rounds to **1.960** |

The two determinations differ by 7.13 cm^-1 (~0.0009 eV), i.e. ~1 sigma given
their combined +/-5 cm^-1 uncertainties -- a normal same-group refinement,
not a discrepancy.

**D_0 (zero-point-corrected) is not directly quoted in either paper's
abstract** (as located via search-engine snippets); it is derived here from
the tabulated Huber-Herzberg/NIST vibrational constants
(omega_e = 1172.2, omega_e x_e = 19.72, omega_e y_e = 0.16 cm^-1; CCCBDB
gives ZPE = 581.63 cm^-1 directly, consistent with
ZPE = 0.5*omega_e - 0.25*omega_e*x_e + 0.125*omega_e*y_e ~ 581.1-581.6
cm^-1):

D_0 = D_e - ZPE ~ 15807.87 - 581.6 = 15226.3 cm^-1 ~ **1.888 eV** (using the
2018 D_e; 1.889 eV using the 2010 D_e). Mark this D_0 figure as
**derived/UNVERIFIED against a directly-quoted D_0** -- it was not read off
either paper, only reconstructed from spectroscopic constants.

## Verdict

**1.961 eV in Paper 58 is D_e, not D_0, and it is correct** -- it reproduces
the Huang et al. (2010) spectroscopic value D_e = 15815 +/- 5 cm^-1 to the
last rounded digit. The number itself does not need to change. What needs
fixing is the **citation**: `\cite{huberherzberg}` does not carry this value
(confirmed above; the NIST/Huber-Herzberg page only has an extrapolated
2.1 eV footnote and 1.88-1.92 eV computed estimates), so the footnote's "not
yet verified" flag is justified as written and should be resolved by citing
Huang et al. 2010 (D_e = 15815(5) cm^-1 = 1.961 eV) directly, optionally
noting the refined Chu et al. 2018 value (D_e = 15807.87(5) cm^-1 =
1.960 eV) as the more current number if the PI wants Paper 58 to track the
latest determination instead of the originally-inherited one. Either choice
is a primary-source-backed number; the framework's currently-tabulated
1.961 eV happens to match the earlier (2010), not the later (2018),
determination.

## Caveats / what was not independently verified

- Full-text abstracts of both 2010 and 2018 papers were read only via
  search-engine-indexed snippets (AIP paywall returned HTTP 403 to WebFetch
  on both `pubs.aip.org/aip/jcp/...` and `pubs.aip.org/jcp/...` URL forms;
  PubMed pages returned a cookie-consent wall to WebFetch). The D_e numbers
  above are corroborated across >=3 independent search results plus the ADS
  abstract-page metadata and the internally-consistent NaD isotope-shift
  arithmetic (15822 - 7 = 15815), so confidence is high, but a direct PDF
  read was not achieved in this pass.
- No DOI, author, or value above was fabricated; all DOIs were confirmed via
  a live Crossref API bibliographic-search query
  (`api.crossref.org/works?query.bibliographic=...`), independently of the
  WebSearch summaries.
- D_0 = 1.888 eV is a derived quantity (D_e minus a Dunham-expansion ZPE),
  not a number read directly from either primary paper's text -- flagged
  UNVERIFIED-as-D_0 above.
