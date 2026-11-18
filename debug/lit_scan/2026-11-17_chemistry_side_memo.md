# Prior-art scan, chemistry / special-functions side, for the 2026-11-14..17 results (v6.7.26-v6.7.36)

Date: 2026-11-17. Scope: the R-operators of Paper 58 `sec:census`, exchange beyond 1s and at unequal exponents,
the (AA|BC) and four-centre rows of `docs/chemistry_period_inventory.md`, and the methods question.
Earlier memos are not repeated; this one only adds to them: `angle_ode_prior_art_memo.md`,
`matuzas2026_bond_angle_prior_art_memo.md`, `prior_art_bond_angle_series_memo.md`,
`water_contact_prior_art_memo.md`, `water_pinch_prior_art_memo.md`, `elliptic_eri_feynman_memo.md`.

Reading levels:
- **VERIFIED** = the arXiv abstract page or the full text (ar5iv) was opened and the quoted sentence read.
- **ABS** = a publisher or institutional landing page with the abstract was opened.
- **SEC** = known only from a search-engine summary or a citing source.

No files were downloaded.

## Q1. ODEs in R for two-centre Slater integrals

| Our result | Status | Reference | What it says | Credit to add | Shortcut / test |
|---|---|---|---|---|---|
| Two-centre two-electron Slater integrals obey a linear ODE in R | **PARTIALLY KNOWN** | K. Pachucki, "Two-center two-electron integrals with exponential functions", PRA 80, 032520 (2009), arXiv:0905.1635. **VERIFIED** (abstract) | "We show that the master integral with the single negative power of all distances can be obtained" from "the second order differential equation in r, the distance between nuclei" (abstract). | **Yes: Paper 58's R-operator paragraph.** Paper 59 already cites `pachucki2009` (l. 395, 835); Paper 58 does not. | See the next row. |
| (same) | **PARTIALLY KNOWN** | K. Pachucki, "Correlated exponential functions in high precision calculations for diatomic molecules", PRA 86, 052514 (2012), arXiv:1209.1258. **VERIFIED** (abstract + ar5iv text) | The general two-centre master integral f(r), with e^{-w1 r12}/r12 and four electron-nucleus exponentials each over its own distance, satisfies a fourth-order inhomogeneous ODE in r, Eq. (13). At w1 = 0 (the ordinary molecular integrals) it reduces to second order: "[16 u w x y d/dr r d/dr + 16 r (...)] f = F". The small-r expansion is f = sum r^k [f_k^(1) ln r + f_k^(2)] (Eq. 37), "absolutely convergent" (credited to Fromm-Hill). For the exchange case (w = u, x = y), "sigma00 vanishes", so his recursions degenerate there. | Same paragraph. Our exchange operator (order 4 at equal exponents, order 9 at unequal) is for the integral itself: homogeneous, with polynomial coefficients. His equation is for a generating master integral: inhomogeneous, low order. Exchange is a 4-fold derivative of his master in the u-parameters at w1 = 0. | **Shortcut:** start from his w1 = 0 second-order equation, differentiate in the parameters and eliminate the source. That gives a derivation route (an upper bound on the order) for our guessed order-9 operator. **Test:** check that the order-9 operator follows from the equation obtained this way. |
| ODEs for the Neumann building blocks; logs at small R | **PARTIALLY KNOWN** | M. Lesiuk, R. Moszynski, "...revisited. II. Neumann expansion of the exchange integrals", PRE 90, 063319 (2014), arXiv:1410.6113. **VERIFIED** (abstract + ar5iv) | "We derive two differential equations which are obeyed by the most difficult basic integrals" (abstract). Eq. (28): alpha^2 L_mu'' + 2 alpha L_mu' - [mu(mu+1) + alpha^2] L_mu = -e^{-alpha}, with alpha = (R/2)(zeta_a + zeta_b) (Paper I, 1410.6112, VERIFIED). The small-alpha particular solution carries log(2 alpha) times a power series, and "the presence of the logarithmic terms is necessary". There are no recurrences in R. | Mention it in the same paragraph as the per-term origin of the R = 0 logarithm. alpha is proportional to R, so each building block has a Bessel-type ODE in R. | A per-term order bound: the tau-series terms are inhomogeneous Bessel equations. |
| Order = number of distinct decay rates (exponential classes); an operator exists for any function built from these atoms | **KNOWN (generic mathematics)** | R. P. Stanley, "Differentiably finite power series", Eur. J. Combin. 1, 175-188 (1980). **SEC** (bibliographic record from two sources) | D-finite functions are closed under sum and product. | Paper 58 already says the operator "exists for any function built from these atoms" and claims only the measured order and singularities. Cite Stanley for the closure fact. | none |
| R = 0 regular singular with a log; integral non-analytic at the united atom | **PARTIALLY KNOWN (analogue)** | W. Byers Brown, E. Steiner, J. Chem. Phys. 44, 3934 (1966), doi:10.1063/1.1726555. **SEC** (abstract via search; ADS and Scilit refused fetch) | One-electron heteronuclear energy near the united atom contains R^5 log R: "not an analytic function of R at the united atom". | Optional context sentence. It concerns an energy, not an integral. | none |
| Exchange: R^9 ln R (1s), R^11 (2p), R^13 (3d, 2p-2p); coefficient (8/3375) a^5 b^5; ratios 32 and 128 | **NOT FOUND** | Searched: Sugiura small-R, united-atom exchange expansions, Lesiuk II, Pachucki 2012 (logs present but no coefficients for exchange), Rothstein arXiv:1207.5866 (VERIFIED abstract; nothing on small R). | — | Keep "not found", qualified: Pachucki's Eq. (37) guarantees log terms in the master integral, so log terms as such are not new. The onset powers and coefficients are new. | none |
| Order 4 / 9 exactly, apparent finite points, rank-1 irregular point at infinity with rates -(a+b), 2a, 2b, ...; no a - b rate | **NOT FOUND** | Harris 2002 (IJQC 88, 701; **SEC**, author abstract page unreachable) and Lesiuk 2016 (J. Math. Chem. 54, 572, arXiv:1507.05556, **VERIFIED** abstract) are about numerical stability and large-order Neumann terms. No D-finite or holonomic treatment of STO integrals was found in any search. | — | none | none |

**Large-R side.** Herring-Flicker and Kay-Silverstone VI (J. Chem. Phys. 53, 4269 (1970),
doi:10.1063/1.1673933, **ABS**) give asymptotic expansions at large R for four-centre r12^-1 STO integrals. These are
consistent with the rank-1 irregular point at infinity; neither states an ODE. Kay-Silverstone VI is optional context
for the four-centre row.

## Q2. Exchange at unequal exponents; exchange beyond 1s

| Our result | Status | Reference | What it says | Credit | Shortcut / test |
|---|---|---|---|---|---|
| The Neumann series terminates at equal exponents and not at unequal ones | **KNOWN** | Lesiuk-Moszynski II (above). **VERIFIED** | "The series in Eq. (3) is, in general, infinite"; "in the limit of beta1=0 or beta2=0 it becomes finite". For an (AB\|AB) density, beta = (R/2)(zeta_a - zeta_b), so beta = 0 is exactly equal exponents. | **Yes.** Paper 58's sentence "the Neumann series ... does not terminate" should cite lesiuk2014b for this. It is already in the bibliography. | none |
| Large-order structure of the non-terminating series | **PARTIALLY KNOWN** | M. Lesiuk, J. Math. Chem. 54, 572 (2016), arXiv:1507.05556. **VERIFIED** (abstract) | A "systematic approach to calculation of the higher-order terms" via large-order expansions; "all auxiliary quantities can be calculated analytically". | Optional. It is the nearest prior description of the tail that our order-9 operator encodes. | **Test:** the large-order asymptotics of the tau-terms should match the operator's behaviour. |
| Exchange beyond 1s stays in the {E1, ln, gamma} class | **KNOWN** | H. J. Silverstone, K. G. Kay, "...III. (2-2)-type three-center integrals", J. Chem. Phys. 48, 4116 (1968). **ABS** (Bar-Ilan CRIS page) | Fourier-convolution evaluation; "the two-center exchange integral is derived as a special case", valid for general n, l, m, zeta. The result is in exponential integrals and spherical-Bessel-type functions. Also Lesiuk II (general mu, sigma) and Maslen-Trefry (via Lesiuk). | Add Silverstone-Kay III next to lesiuk2014b for general (n, l, m) exchange. | none |
| Higher-l exchange has the same operator shape as 1s (same order, same rates, only the degree rises; log onset moves R^9 -> R^11 -> R^13) | **NOT FOUND** | No statement of shared analytic structure across l in Lesiuk II, Silverstone-Kay III (abstract) or Pachucki. Podolanski 1931 and Kotani's tables could not be reached online (Kotani is cited in Lesiuk II only as a recursive scheme). **UNVERIFIED** for those two. | — | Keep as new. | none |

## Q3. Three-centre integrals as functions of the bond angle

Most of this was settled by the earlier memos. Matuzas 2026 (Zenodo) has the 3/2 and 7/2 exponents, the location of
x*, and the rho^l l^-4 and l^-8 laws for (XY|XZ) and the one-electron integral. Silverstone PNAS 1967 has the Legendre
channel form of (1s_a1s_a|1s_b1s_c), which is our (AA|BC). Bouferguene-Jones 1998, Barnett-Coulson, Kutzelnigg and
Hill are already credited.

| Our result | Status | Reference | What it says | Credit | Shortcut / test |
|---|---|---|---|---|---|
| (AA\|BC): 21 fixed constants in every Legendre coefficient | **PARTIALLY KNOWN** | Silverstone, PNAS 58, 34 (1967). VERIFIED in `prior_art_bond_angle_series_memo.md` | Each channel is a finite expression in exp, E_n and Bessel-type functions for exactly this integral. | Already owed by that memo. Confirm that Paper 59's (AA\|BC) paragraph cites it. | none |
| Leading branch (x* - x)^{7/2} gives rho^L L^-4; 9/2 when it cancels; L^-5 in P_L^1; equal arms (1 - x)^{7/2} -> L^-8 | **PARTIALLY KNOWN** | Mechanism: Matuzas 2026 (the same 7/2 for (XY\|XZ)). Rate transfer: A. Sidi, "Asymptotic expansions of Legendre series coefficients for functions with interior and endpoint singularities", Math. Comp. 80, 1663-1684 (2011). **SEC** (AMS page 403; PDF certificate error) | Sidi: coefficient asymptotics "determined exclusively" by the expansions of f at its singular points, on [-1, 1]. For an exterior branch point (x* > 1) the rule is classical Darboux. No primary source for that case was opened. | Cite Matuzas for the mechanism. The singularity-to-decay rule is textbook (Darboux; Sidi for endpoints). The specific (AA\|BC) 7/2 at x*, the 3/2 cancellation, the zB D = 9 cancellation and the p/d cases were **not found**. | none |
| An ODE in the bond angle (order 8 one-electron; order 14 minimal for (AA\|BC)) | **NOT FOUND** | `angle_ode_prior_art_memo.md`, re-checked: no new source. Matuzas analyses local germs only. | — | Keep "not found". | none |

## Q4. Four-centre integrals

| Our result | Status | Reference | What it says | Credit | Shortcut / test |
|---|---|---|---|---|---|
| Channel/Fourier representation of four-centre STO integrals | **KNOWN** | K. G. Kay, H. J. Silverstone, "...V. Four-center integrals by Fourier-transform method", J. Chem. Phys. (1969). **SEC** (search summary of the abstract; volume/pages not confirmed; Part IV is J. Chem. Phys. 51, 956 (1969), doi:10.1063/1.1672162) | An infinite sum in which "the internuclear angles appear in spherical harmonics, and the internuclear distances in integrals over spherical Bessel functions and exponential-type integrals", evaluated as "convergent infinite expansions by contour integration". | Credit in Paper 63 `sec:four_centre` alongside Silverstone 1967 PNAS, once the volume and pages are confirmed. | **Test:** its radial integrals are where a depth-2 iterated E1 would surface. Reading the full text could confirm or refute FX's depth-2 key class directly. |
| Contact branch points (x* - x)^{7/2} at moving Y-W and Y-Z contacts | **NOT FOUND** | Hoggan, IJQC 110 (2010), doi:10.1002/qua.22213 (**SEC**; Wiley 403): reduces to two-centre overlap-like products. W. Chen, L. Silberman, U. Segal, C. Greif, arXiv:2610.12608 (Oct 2026, **VERIFIED** abstract): Fourier route, "only a one-dimensional quadrature problem" remains; no analytic or convergence-rate statement. Safouhi/Berlu: acceleration only (SEC). | — | Keep as new. | none |
| Depth-2 iterated exponential integrals and E1 x E1 product constants appear; the four-centre integral leaves the single-E1 class | **NOT FOUND** | No chemistry source classifies four-centre STO constants. Nearest is J. Matuzas, arXiv:2610.00036 (Sep 2026, **VERIFIED** abstract): the *cyclic four-electron atomic* Slater integral in closed form in log, Li2 and Li3 of rational functions of the exponents; at equal exponents "a weight-three period over" Q, ramified at 2 and 3. | Not a four-centre integral, but directly relevant to Paper 63's ring and chain results (a Coulomb ring of four?). Not cited anywhere in the corpus (grep). **Flag for Paper 63 / atlas.** | **Test:** compare with our ring-of-four STOP and the ring-of-three elliptic drop. Matuzas's cyclic integral being polylogarithmic (no elliptic) would be a direct check on the corpus's ring classification. |

## Q5. Molecular integrals as Feynman integrals / periods / D-modules; guessing mod primes

| Our result | Status | Reference | What it says | Credit | Shortcut / test |
|---|---|---|---|---|---|
| IBP / differential-equation methods applied to Slater integrals | **KNOWN (two-centre and atomic)** | Pachucki 2009/2012 (above); Puchalski-Kedziera-Pachucki (Li, arXiv:0909.5542); Pachucki-Puchalski-Remiddi (already in Paper 63) | IBP identities used for exponentially correlated integrals. | Paper 58 needs Pachucki. The rest is already credited in Papers 59 and 63. | none |
| Guessing ODEs from exact data modulo primes + rational reconstruction | **KNOWN as a method; NOT FOUND applied to molecular integrals** | Method: Kauers/BKKS (already in `angle_ode_prior_art_memo.md`); Belitsky-Smirnov-Yakovlev arXiv:2303.02511 (SEC) for modular reconstruction in IBP. Searches "guessing / modular / holonomic + Slater / Hylleraas / molecular" returned nothing in chemistry. | — | Cite the guessing method where the mod-p operators are stated. | none |
| Slater Coulomb integrals as periods | **PARTIALLY KNOWN** | Matuzas arXiv:2610.00036 (above) uses "period" language for an atomic four-electron integral. Brown and BBBG are already credited. Calderini/Aquilanti Sturmian work: the searches surfaced nothing on R-recurrences or ODEs (SEC). | — | Matuzas 2610.00036, as in Q4. | none |

## Top 5 things the PI should know

1. **Paper 58's R-operator paragraph is missing Pachucki.** PRA 80, 032520 (2009) and PRA 86, 052514 (2012) derive
   linear ODEs in the internuclear distance for the general two-centre two-electron exponential master integral:
   second order at w1 = 0, fourth order with e^{-w1 r12}, inhomogeneous. The 2012 paper also gives a small-r expansion
   with ln r terms. Paper 59 already cites the 2009 paper; Paper 58's "[SYMBOLIC] Each class also obeys a linear
   differential equation in R" paragraph does not. What stays new: homogeneous polynomial-coefficient operators for
   each integral class, their exact orders (4 and 9), the apparent singular points, the rates at infinity, and the
   log onset and coefficients.
2. **The non-terminating Neumann series at unequal exponents is stated by Lesiuk-Moszynski II.** The series is
   "in general, infinite" and finite when beta = 0, which for exchange means equal exponents. lesiuk2014b is already
   in Paper 58's bibliography; it should be cited at that sentence.
3. **Silverstone-Kay III (JCP 48, 4116, 1968) derives two-centre exchange for general n, l, m, zeta.** So "beyond 1s
   stays in the class" is known. The shared operator shape across l (same order and rates, the log onset rising by
   two per unit of l) was not found.
4. **Matuzas, arXiv:2610.00036 (Sep 2026), is uncited anywhere in the corpus.** It gives a closed form for the
   cyclic four-electron Slater integral in Li2 and Li3; at equal exponents it is a weight-three period. It bears on
   Paper 63's ring results and the atlas. This is the second Matuzas preprint overlapping the project's territory;
   a periodic check of his output is warranted.
5. **There is a concrete shortcut and test from Pachucki's equation.** Differentiating his w1 = 0 second-order
   equation in the exponent parameters and eliminating the source should give a *derived* upper bound for the
   unequal-exponent exchange operator, which is now guessed mod p. That could raise the order-9 claim from
   MEASURED toward SYMBOLIC. Separately, Kay-Silverstone IV/V (1969) is the place to look for whether depth-2
   iterated exponential integrals appear in the four-centre radial integrals; its full text is paywalled and was
   not read.

Unverified items to close before citing: Kay-Silverstone V volume and pages; Harris 2002 abstract; Byers
Brown-Steiner abstract wording; Sidi 2011 statement; Podolanski 1931 and Kotani (not reached).

## Corrections after primary-source verification (2026-11-17, citation-reviewer + PM)

- Hannesdottir-Mizera, arXiv:2204.02988, Eq. (7.12): the exponent is called rho, with d defined at Eq. (7.3). The
  source states the leading form only for rho <= 0 and for isolated non-degenerate saddles, and calls rho > 0
  non-singular. It gives no support for the positive half-integer exponents at x*; do not cite it for them.
- The memo's adapted count predicting a 5/2 term in (AA|BC) was wrong: it gave A's density cusp a power 1/2 in the
  same integral as C's cusp, but that cusp reaches the second electron only through 1/r12. Done correctly the count is
  the corpus's own segment-average mechanism, and the 5/2 coefficient vanishes identically (the potential of a bounded
  density has no s^-1 or s^1 term). See `debug/l52_design_memo.md`.
- Fresan-Sabbah-Yu (Duke Math. J. 171, 2022), Prop. 2.7: Sym^k Kl_2 is irreducible for all k, over C. Nothing found
  identifies that connection with Paper 63's F_n operator, so it does not lift the n <= 8 limit.
- Silverstone-Kay III is J. Chem. Phys. 48, 4108-4115 (1968), not 4116. Its abstract gives an infinite series and
  does not state a constant class.
- Pachucki 2009/2012 read in full: the exchange integral at any exponents lies inside the master integral, and at
  exchange exponents the second-order equation has p^2 = 0 (two quadratures of an explicit source). Followed up in
  `debug/pachucki_exchange_bound_memo.md`.
- Content not reached (exist as cited; nothing should rest on what they say): Podolanski 1931, Kotani tables,
  Harris 2002, Kay-Silverstone V, Sidi 2011 full text.
