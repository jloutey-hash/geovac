# Pre-registered prediction: a production six-Gaussian geminal fit on the exact correlated triangle (PG1-PG3)

**Registered:** 2026-10-26, before any computation of Y in this file. PI direction: "I rather thought we were better
equipped to report on the cusp and tail as compared to gaussians", then, on the proposed test, "Let's do both."
**Prompted by:** record GF (`2026-10-25_gaussian_fit_geminal.md`), whose cap was that the production weight function
was not reproduced.
**Extends:** `2026-10-25_gaussian_fit_geminal.md`.
**Commit status:** committed on its own, before any computation of Y in this file.

## The question

Record GF fitted exp(-r) by six Gaussians under two weights of our choice and measured the error on the exact Y.
This record uses a fit that a production code actually ships, instead of reconstructing the weight.

## What I know before computing (disclosed)

- Source: Psi4, `psi4/src/psi4/libmints/integralparameters.cc`, class `FittedSlaterCorrelationFactor`, read today
  from the master branch on GitHub. It hard-codes six exponents 0.2209, 1.004, 3.622, 12.16, 45.87, 254.4 and
  coefficients -0.3144, -0.3037, -0.1681, -0.09811, -0.06024, -0.03726. The exponents are multiplied by the square of
  the Slater exponent. There is no comment on their source or weight. The sign is the code's convention for the
  correlation factor; this record fits exp(-r) with the coefficients' magnitudes.
- Not verified: the literature source of these numbers. Molpro's manual says its default is six Gaussians fitted
  with "a suitable weight function". TURBOMOLE's MP2-F12 paper says its six-term fit comes from a reference for
  gamma = 1, scaled. A search for all six numbers together found no source. Libint evaluates the Slater geminal
  directly, with a Gaussian fit only for exceptional parameters. So this is one code's shipped fit, used here as a
  stand-in for production practice. Which integral classes Psi4 uses it for is not checked.
- Pointwise, computed today from the three six-term fits, not from Y. Relative error of each fit against exp(-r):

  | fit | r = 0 | r = 1 | r = 2 | r = 3 | r = 4 | r = 6 | smallest exponent |
  |:--|:--|:--|:--|:--|:--|:--|:--|
  | W1 (plain norm, record GF) | -6.0e-2 | +1.6e-5 | +2.9e-4 | -1.2e-3 | +2.2e-3 | +1.1e-3 | 0.0651 |
  | W2 (weight exp(-r^2), record GF) | -2.4e-2 | -9.7e-6 | -4.4e-4 | -1.1e-2 | -1.9e-1 | -8.0e-1 | 0.1691 |
  | Psi4 | -1.8e-2 | -5.9e-5 | +5.8e-4 | -1.3e-1 | -5.0e-1 | -9.6e-1 | 0.2209 |

  So the Psi4 fit is the best of the three at the cusp and the worst in the tail.
- Record GF measured, at six terms and t = 1/2, 1, 2, 4: W1 -3.6e-5, -1.4e-5, -4.1e-6, -2.4e-6; W2 -5.4e-2, -4.1e-2,
  -2.8e-2, -1.7e-2. W1 has the worst cusp and the smallest error, which is why record GF read the error as set by
  the tail.
- The geminals in Y are exp(-(3/2) r) and exp(-3 r), so gamma r = 4 at r = 2.67 and at r = 1.33.
- Nothing about Y with the Psi4 fit has been computed.

## Model

Exactly record GF's. Y_6(t) = sum_ij c_i c_j Y_GG(alpha_i (3/2)^2, alpha_j 3^2; t) by `debug/gf_engine.py`, with
{c_i, alpha_i} the Psi4 numbers, against the exact Y(t) of the 60-digit grid, at t = 1/2, 1, 2, 4.
E(t) = Y_6/Y - 1.

## Predictions

**PG1 (size).** 1e-2 <= |E(t)| <= 3e-1 at each of the four t.
- Falsifier: outside that band at any t.
- Confidence: moderate.

**PG2 (the tail sets the order).** E(t) < 0 at each t, and |E(t)| > |E_W2(t)| > |E_W1(t)| at each t. The three fits
are ordered by their tail deficit, the reverse of their cusp deficit.
- Falsifier: a positive E at any t, or the order broken at any t.
- Confidence: moderate to high for the sign and for |E| > |E_W1|; moderate for |E| > |E_W2|.

**PG3 (diffuse apex, larger error).** |E(1/2)| > |E(1)| > |E(2)| > |E(4)|.
- Falsifier: the order broken anywhere.
- Confidence: moderate. It held for W1 and W2 in record GF.

## Novelty status

- That production geminal fits favour short range, since the correlation energy is dominated there: PROBABLY KNOWN.
  It is the stated purpose of a weight function; the primary sources were not opened.
- The size of a shipped fit's error on an exact correlated three-electron integral: NOT FOUND in record GF's search
  or today's, not to the standard of a novelty claim.

## Controls

- Engine: the same Y_GG controls as record GF, rerun on one pair of the Psi4 exponents at two node sets.
- Reproduction: record GF's W1 n = 6 row recomputed alongside, to its printed digits.

## Caps stated in advance

- One configuration, model densities, as in record GF.
- Y weights the geminal over its whole range. A correlation energy is dominated by short range, so a large error on
  Y does not imply a large error in an F12 energy. This record prices the fit on this integral, not on an energy.
- One code's fit. Its weight and its literature source are unverified, and how the code uses it is not checked.
