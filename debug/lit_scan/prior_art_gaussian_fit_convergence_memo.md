# Prior art: convergence of Gaussian fits to exp(-gamma r), and the cost in a correlated integral

Date 2026-11-01. Read-only scan for Paper 59 sec:density_ring (GF and production-fit paragraphs), records
`docs/predictions/2026-10-25_gaussian_fit_geminal.md` and `2026-10-26_production_geminal_fit.md`. Downloads are in
`debug/lit_scan/dl_prior_gauss_conv/` (git-excluded). Paper 59's bibliography cites none of the items below except
Psi4 (`psi4_2020`). The GF record already names Kutzelnigg 1994 (recalled) and McKemmish-Gill 2012; the paper cites neither.

## 1. Is root-exponential convergence a theorem? Yes, as an upper bound, with constants in a weighted norm.

**[A] D. Braess, "Asymptotics for the approximation of wave functions by exponential sums," J. Approx. Theory 83,
93-103 (1995), doi:10.1006/jath.1995.1110.** VERIFIED from the author's PostScript (`braess1995/wave.ps`; only the
text strings were read, the file was not run). Substitute x = r^2. Then e^{-r} = e^{-sqrt x} is completely monotone,
and Gaussian sums become exponential sums in x. Rational approximation (Stahl's e^{-pi sqrt n} law for |x|) then
gives the rate.
- Section 5, Eq. (5.3), verbatim in substance: there is an n-term sum with ||e^{-sqrt x} - u_n|| <= c n^{3/2} e^{-pi sqrt(3n)}.
  The norm is the **weighted L1 norm** ||g|| = int_0^oo r^2 |g(r^2)| e^{-r} dr (Eq. 1.3). Braess calls this "the
  natural choice" because it tracks the hydrogen energy expectation.
- For 1/r (Theorem 4) the rate is gamma = 2 pi/sqrt 3, in the same norm. Eq. (4.4): in the corresponding
  e^{-r}-weighted L2 norm, the 1/r error is only about e^{-pi sqrt(n/6)}. "Our method provides only upper estimates."
  Braess's gamma is "twice the value" of Kutzelnigg's numerical exponent, because free exponents double the rate.
- **Bearing.** The e^{-b sqrt n} law in Paper 59 is a theorem (an upper bound) in the e^{-r}-weighted L1 norm. Braess
  does not state it for our plain norm int r^2 (.)^2 dr, for a pointwise error, or for an integrated functional such
  as Y. **Credit owed.**

**[B] W. Kutzelnigg, "Theory of the expansion of wave functions in a Gaussian basis," Int. J. Quantum Chem. 51,
447-463 (1994).** ABSTRACT ONLY, plus Braess's and Shaw's accounts. The method is the even-tempered trapezoid
discretisation of the Gaussian transforms of e^{-alpha r} and 1/r. The cutoff error goes as exp(-a h) and the
discretisation error as exp(-b/h), so the optimum is root-exponential. Braess quotes Kutzelnigg's numerical law
e^{-pi sqrt(n/3)} for 1/r. Also: Kutzelnigg, IJQC 113, 203 (2013), "...I. Local versus global approximation" (ABSTRACT
ONLY). It states an error of order exp(-a sqrt n), with the constant a depending on the property and on the
optimisation criterion. **Credit owed** as the origin.

**[C] L. K. McKemmish and P. M. W. Gill, "Gaussian expansions of orbitals," J. Chem. Theory Comput. 8, 4891 (2012),
doi:10.1021/ct300559t.** ABSTRACT ONLY (ANU portal). The energy errors of optimal even-tempered expansions are
root-exponential, of order exp(-pi sqrt(3n)) for s orbitals; sqrt(5n) and sqrt(7n) are the abstract's p and d
analogues (SECONDARY). They attribute this to "the orbital cusps at r = 0" and not to the Gaussian tail. A cuspless
Lorentzian also converges root-exponentially, which they put down to its "fat tail". **Credit owed:** the GF record
cites it but the paper does not.

**[D] R. A. Shaw, "The completeness properties of Gaussian-type orbitals in quantum chemistry," IJQC 120, e26264
(2020), doi:10.1002/qua.26264.** VERIFIED (accepted manuscript). It extends Kutzelnigg's trapezoid analysis to all
STOs. The total error goes as exp(-sqrt(K_max M)) in the L-infinity norm (Eq. 25) and the H1 norm (Eq. 34), and as
exp(-sqrt(n M)) in L1 (Eq. 28), with K_max = n - 1/2, so K_max = 1/2 for 1s. Its constants are loose (no pi). It also
says "the error at the origin always dominates". Use it as a secondary credit, not for a constant.

**[E] M. Bachmayr, H. Chen and R. Schneider, "Error estimates for Hermite and even-tempered Gaussian approximations in
quantum chemistry," Numer. Math. 128, 137-165 (2014).** ABSTRACT ONLY. It gives rigorous a priori H1 estimates,
root-exponential, for even-tempered Gaussian approximations of STOs, proved with sinc approximation theory.

**Our rate against the constants.** I ran a small check on `tests/data/paper59_gf_2026_11_01.json` (W1, n = 2..10).
- The fitted slope of ln|E_n| against sqrt(n) is b = -6.15, -6.74, -7.44, -7.82 at t = 1/2, 1, 2, 4.
- The local slopes are still rising at n = 10 (-6.5 to -8.0), so the fits are pre-asymptotic.
- The squared plain-norm residual of the fit itself has fitted b = -7.94, with local slopes from -6.3 up to -8.5.
- The fit error at r = 0 has b of about -1.7.
- So the error in Y follows the **squared** L2 residual, not its amplitude. This is consistent with Y depending
  bilinearly on a least-squares projection (my reading; no source states it).
- The nearest published constant is pi sqrt 3 = 5.44 (Braess's upper-bound rate and McKemmish-Gill's energy rate).
  Our |b| exceeds it, which is consistent with an upper bound on the error.
- No published constant exists for the plain norm or for an integrated functional. Do not claim agreement with pi sqrt 3.
- If Braess's lemma were run at alpha = 3/4 for e^{-sqrt x} with no cutoff penalty, the amplitude rate would be
  pi sqrt(3/2) and its square pi sqrt 6 = 7.70. That is close to our W1 b, but it is **my extrapolation, UNVERIFIED**.
- Under W2 the Y-slopes stay flat at about -1.6 to -2.1. That is the signature of tail-truncation error outside the
  weight's reach, not a fit-convergence law. It agrees with the paper's reading.

## 2. Gaussian-geminal fits in F12 practice

- **May and Manby, JCP 121, 4479 (2004); Tew and Klopper, JCP 123, 074101 (2005); Werner, Adler and Manby, JCP 126,
  164102 (2007)** are the origin of the weighted six-Gaussian fits. Full texts were blocked: the irsamc mirror has a
  bot wall, and the publisher pages returned 403. **The weight function stays UNVERIFIED**, as Paper 59 already caps.
  Molpro 2008 manual sec. 30.8 (SECONDARY, search snippet): GEM_M is "exponent for r in weighting function" (default
  0), and the weight exponent is derived from GEM_BETA.
- **Psi4 source** (`psi4/src/psi4/libmints/integralparameters.cc`, VERIFIED; downloaded to `psi4src/ip.cc`). There
  are six coefficients (-0.3144, -0.3037, -0.1681, -0.09811, -0.06024, -0.03726) and six exponents (0.2209 to 254.4),
  scaled by gamma^2. The code carries no literature reference. The same numbers appear in QCaml (SECONDARY). Their
  origin is still unverified; Tew-Klopper 2005 is the likely source, but that is a guess.
- **C. Wang, arXiv:2508.04635 (2025).** VERIFIED (text): "the equations for fitting the Slater geminal functions
  [Tew-Klopper; Werner-Adler-Manby] without the weighting functions will have the same functional forms as the
  one-electron orbitals". In other words, an unweighted geminal fit reuses the STO-nG exponents. **Credit owed** for
  Paper 59's "which at three and six terms gives the STO-3G and STO-6G exponents". Hehre, Stewart and Pople, JCP 51,
  2657 (1969) also deserves a citation there (standard; the record names it; the paper does not cite it).
- **Persson and Taylor, JCP 105, 5915 (1996)**, via Barca and Loos (arXiv:1704.08051, VERIFIED text): "a six- or
  ten-term fit introduces a 0.5 mEh or 20 microEh error, respectively". That fit is of linear r12, priced on an
  **energy**. Barca and Loos also report Ten-no's finding that "a 10-term GG fit of S12 yields very similar results".
- **Error propagation into three-electron integrals.** None found. The searches turned up three-electron
  integral machinery (Barca and Loos 2017; Persson and Taylor 1997), energy-level fit errors (Persson and Taylor), and
  geminal-exponent studies (Powell, Surjuse, Gaudel and Valeev, arXiv:2506.10780, 2025). None measures a geminal-fit
  error on an individual three-electron (triangle) integral against an exact reference. Paper 59's GF measurement
  appears new at the level of a brief search, not a novelty-grade one.

## 3. Build-on items

1. **Trapezoid (sinc) quadrature of the Gaussian transform instead of least squares.** VERIFIED in Powell et al.
   (arXiv:2506.10780, Eq. 3). Libint 2.11 represents the STG out of range as a Gaussian sum, using "trapezoidal
   quadrature of its integral representation exp(-beta r)/beta = (2/sqrt pi) int ds exp(s - beta^2 e^{2s} - r^2
   e^{-2s}/4)". The step is h = (0.2 - 0.5 log10 eps)^{-1}, with T = 26, giving relative precision 1e-12 for r in
   [1e-5, 1/beta] and absolute precision beyond. Sources: Harrison, Fann, Yanai and Beylkin, ICCS 2003, LNCS 2660,
   103; Beylkin and Monzon, ACHA 19, 17 (2005); used in Bischoff, JCP 141, 184106 (2014). This is a controlled,
   tail-uniform alternative to W1/W2 and can be put directly on Y_n. A test: it converts the tail problem into a
   known parameter (the s-interval).
2. **The rate constant to test.** The rigorous number is Braess's pi sqrt 3 in the e^{-r}-weighted L1 norm. Refit
   under that exact norm (n = 2..16), then compare ln of the fit error and of E_n against -pi sqrt(3n) + (3/2) ln n.
   A tail-uniform fit should also be checked against the sinc rate exp(-pi sqrt(...)).
3. **Best L-infinity fits (Remez/Chebyshev).** Braess and Hackbusch, IMA J. Numer. Anal. 25, 685 (2005),
   doi:10.1093/imanum/dri015 (VERIFIED, MPI preprint 75/2004). For 1/x on [1, infinity) they prove an
   exp(-c sqrt k) bound and conjecture c = pi sqrt 2 ("E_k ~ 6.7 log(2+k) e^{-pi sqrt(2k)} is a very good
   approximation"). Their table fits c2 from 4.01 to 4.36, still rising, which is the same pre-asymptotic drift seen
   in our b. The 2009 survey by Braess and Hackbusch (DeVore and Kunoth, eds., pp. 39-74) was NOT RETRIEVED.

## Paper 59 sentences touched

- "With the fit in the plain norm ..., which at three and six terms gives the STO-3G and STO-6G exponents": cite
  Hehre-Stewart-Pople 1969, and optionally Wang 2025 for the geminal reuse.
- "The error falls like e^{-b sqrt n} in the number of terms n, as the cusp at r = 0 predicts": cite Kutzelnigg 1994,
  Braess 1995 (theorem, weighted L1 upper bound pi sqrt 3), McKemmish-Gill 2012. Do not let b read as a known
  constant.
- "So on this integral the error is set by how far out the fit is trusted, not by the cusp": contrast with
  McKemmish-Gill, who credit the cusp. Both hold, because they price different quantities (orbital energy against a
  long-range integral).
