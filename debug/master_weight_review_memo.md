# Cold review: weight of the two-centre master integral f0 (debug/master_weight_memo.md)

Date 2026-11-19. Reviewer did not write the memo. No paper/test/doc/git state touched.
Script: `debug/master_weight_review_checks.py` (independent re-implementation; imports nothing from the
reviewed script; Phi and E2 by quadrature, not by the double series; own transcription of label 60).
Data: `debug/data/master_weight_review_checks.json` (+ `.log`), dps 50, 60 s.

## Verdicts

| Claim | Verdict |
|---|---|
| (CF) closed form | **SOUND** |
| (A) Theorem A | **SOUND** (holds for complex rates; even extends to shifted Ei arguments) |
| (B) Theorem B | **UNPROVED** (plausible sketch; must not carry any paper sentence) |
| (W) "weight exactly two" | **SOUND WITH CAVEATS** (function/iterated-integral sense only; not period sense) |
| (L) a -> b cancellation | **SOUND** (now re-run numerically, O(eps) approach) |
| (Lit) | Pachucki and Lesiuk-Moszynski claims **verified**; Ruedenberg/Harris/Maslen-Trefry unverified |

## (CF) Closed form -- tried, held

1. **Label 60 transcription.** Fetched arXiv:0905.1635 source (`rec_h2.tex`) myself. Label 60 at u = a, w = b:
   (e^{2bR} - e^{2aR}) Ei(-2(a+b)R) - e^{-2bR}[Ei(2(b-a)R) - 2Ei((b-a)R) + ln|(b-a)/(b+a)|]
   + e^{-2aR}[Ei(2(a-b)R) - 2Ei((a-b)R) + ln|..|], prefactor 1/(2(u^2 - w^2)). The memo's six (k, lam, mu) and the
   LC(e^{-2bs} - e^{-2as}) term match term by term, and match `Calc.FM1` in `tests/_pachucki_extension.py` (read
   by eye). Label 59: f = int_r^oo F^(-1)(r')/r' dr'. Held.
2. **Per-term T, re-derived.** Writing Ei(mu s) = gamma + ln|mu| + ln s + Ein(mu s) (real mu), splitting
   int_R^oo = FP int_0^oo - int_0^R, the subtraction int_0^R s^{eps-1}(g + ln s) ds produces exactly the
   -g ln R - ln^2 R/2 terms, and int_0^R [(e^{lam s}-1)(g + ln s) + e^{lam s} Ein(mu s)] ds/s
   = g Ein(lam R) + ln R Ein(lam R) - E2(lam R) + Phi(lam R, mu R). The Mellin finite part, from
   int_0^oo s^{eps-1} e^{-ps} E1(qs) ds = Gamma(eps) q^{-eps}[1/eps + eps Li2(-p/q) + O(eps^2)], is
   (gamma + ln q)^2/2 + pi^2/12 + Li2(-p/q). Matches the memo's formula and sign conventions.
3. **Numerics at fresh points**, quadrature vs closed form (my Phi/E2 by quadrature):
   (lam, mu) = (-5/3, -1/2), (3, -9/2) [mu < 0, one with lam > 0], (-7/2, 3/2), (-5/2, 5/4) [mu > 0, PV Ei],
   R = 0.25, 2.3: |diff| <= 3e-49 at dps 50. The constant K, obtained independently by inverting the closed
   form at r0 = 0.02 from quadrature, matches the Li2 formula to <= 2.4e-50 in all four cases, including the
   mu > 0 continuation (-5 pi^2/12, Re Li2(p/mu) with p/mu > 1).
4. **Assembled f0 at new pairs** (3/4, 2) and (1, 5), R = 0.2, 0.9, 2.2: closed form vs my own quadrature of
   label 60: rel <= 5e-46 except (1, 5), R = 2.2: 3.6e-40 (cancellation of e^{2bR}-size pieces; still 39
   digits). My quadrature vs the helper's `f0_quadrature`: <= 8e-49. Held at >= 39 digits (>= 45 except the
   one large-bR point).
5. **C(a, b) = 0** at five pairs including two new rational ones, (2, 3), and the irrational ratio (1, sqrt 2):
   |C| <= 1e-50. Sensitivity control: flipping the LC sign gives C = -0.225 at (3/4, 2), so the check bites.
   f0(1e-10) = 5.04e-11 by my quadrature (linear, f0(0) = 0). Note: C = 0 is not merely "observed": T is proved
   (derivative + Mellin limit), so C = lim_{R->0} f0, and f0(0) = 0 follows from label 58 (factor r, finite
   integrand at r = 0). C = 0 is therefore a consequence, i.e. a proved six-term Li2 identity; the memo
   undersells it. Not derived by hand as a dilogarithm identity by anyone.

## (A) Theorem A -- tried, held

Attacks:
* **Local form of Ei(nu R) for complex nu.** On R > 0, with any fixed continuous branch, Ei(nu R) = gamma +
  log(nu) + ln R + Ein(nu R) + (constant 2 pi i k or +-i pi from the convention); Ein is entire for every complex
  argument. All branch/convention choices change only a constant, absorbed into Q. Held.
* **Monodromy step.** Both sides continue to a slit punctured disc at 0; (P + Lam) ln R = M with M single-valued
  meromorphic near 0 forces P + Lam = 0 near 0, hence identically. Only local meromorphy at 0 of the non-log
  part is needed. Held. (The log coefficient of f0 does not even need the full closed form: f0' = -FM1/R and
  FM1's ln s coefficient is c sum k_i e^{lam_i s}, giving -Lam' = -L/R directly.)
* **Coinciding rates.** S is a set; any atoms with rate 2b are merged into rho_{2b}. Distinct complex
  exponentials are independent over C(R). Held.
* **rho' + 2b rho = -c/R.** A pole of order m >= 1 anywhere gives order m + 1 on the left; at R = 0 that is >= 2
  vs order 1 on the right; a polynomial rho gives a polynomial. No rational solution. Held (needs b != 0).
* **Extension found (strengthens A):** atoms e^{wR} Ei(nu(R - R0)) or e^{wR} ln(R - R0) with R0 != 0 (what a
  one-fold integral of a rational-times-exponential with poles off 0 produces) are analytic at 0, so they enter Q
  only; the argument is unchanged. So A covers the full "depth-one" class of int r(s) e^{ws} ds primitives.

No counterexample found.

## (B) Theorem B -- not checked to rigor; must not carry a claim

The shuffle obstruction in step (3) (the e^{2bR} Ei(-2(a+b)R)/R term needs its partner e^{-2(a+b)R}Ei(2bR)/R,
which FM1 lacks) is the right idea and is the standard integrability test for weight-two shuffle products. The
weak joints are those the memo names: the multi-primitive polynomial lemma, the handling of log-type primitives,
and the Li2 step (Li2(u_k) independence over F1). The literature citations for the tools are from memory and
unverified. Recommendation: no paper sentence may say "no dilogarithm form exists", "genuinely new", or "not
expressible through products of weight-one functions". Theorem A is enough for the sentence the paper needs.

## (W) The word "weight"

Paper 58's usage (abstract l.127; l.555-574; l.669-672) is informal and mixes functions and constants:
"weight one" = the seed set {exp, E1, ln, gamma}; "iterated integrals ... generically land at weight two
(dilogarithms, zeta(2))"; "weight-one filtration" = closure under C(R)/Q-linear combination of the seeds. The
census sentence the paper currently carries is a statement about **functions of R**.

* Lower bound: Theorem A is exactly "not in the C(R)-linear span of the weight-one seed functions, any rates".
  That is what "not weight one" must mean in the paper's own sentence ("not a combination of these seeds with
  rational coefficients"). Supported.
* Upper bound: f0 = c sum k_i int_R^oo e^{lam_i s} Ei(mu_i s) ds/s + c LC[E1 - E1], with Ei itself a one-fold
  integral of e^t dt/t: depth <= 2 over the kernels e^{w s} ds/s. Supported.
* So "weight two" is supported **in the sense of iterated-integral depth of functions of R**. It is not supported
  in a period/motivic sense: Theorem A says nothing about f0(R0) at a fixed algebraic R0 (whose numerical value
  could in principle reduce), and "weight" of values is not defined in the paper for E1 values anyway.
* "Exactly two" for the exchange integrals: K = c1 f0 + (C(R)-combination of g2..g9), the span of g2..g9 is
  closed under d/dR (table in `debug/exact9_independence_memo.md`), and order exactly 9 (proved in Paper 58)
  forces c1 != 0; so K is not in A's span. Sound for every case where order exactly 9 is proved **and** the
  integral lies on the nine functions: the 1s case, the four sigma = 0 cases, and the two m != 0 cases now
  listed in Paper 58 (the memo names only five; the paper now has seven). Weight two in the product sense
  (e.g. -Lam ln R, a product of weight-one functions) is not excluded and need not be.
* Wording to avoid: "no dilogarithm form exists", "genuinely new", "weight exactly two" without the sense stated,
  and "Li2 does not appear" unqualified -- correct form: "the dilogarithm constants of the six terms cancel".

## (L) a -> b -- re-run, held

With a = 1, b = 1 + eps, R = 0.7 and 1.9, W2 = -c sum k_i Phi(lam_i R, mu_i R) vs the memo's predicted limit
R/(2a) d_xPhi(2aR, -4aR), d_xPhi = e^x Ein(y)/x - (Ein(x + y) - Ein(x))/x (integration by parts checked against
quadrature to 5e-51):

| R | eps = 1e-2 | 1e-3 | 1e-4 | 1e-5 |
|---|---|---|---|---|
| 0.7 | (W2 - lim)/eps = -0.2881 | -0.2871 | -0.28704 | -0.28703 |
| 1.9 | -30.79 | -30.51 | -30.483 | -30.480 |

W2 tends to the weight-one limit at O(eps) (the 1/eps prefactor c is cancelled). Caveat: this shows the Phi part
of f0 degenerates to weight one; for K itself c1 may also grow as eps -> 0, so the statement about K rests on the
existing symbolic Sugiura-limit result (`debug/pachucki_extension_memo.md`), not on this check.

## (Lit)

* **Pachucki arXiv:0905.1635** (source opened): f at the exchange point is given only as label 59's integral of
  label 60; no closed form for u != w. "Li" appears only in label 10 (the three-body G, two Li terms). Verified.
* **Lesiuk & Moszynski arXiv:1410.6113** (source `neumann.tex` + supplemental opened): exchange integrals via the
  Neumann expansion; basic quantities in E1, gamma_E and finite/infinite power series (e.g. eq. near l.634,
  an e^{-alpha} sum_t M_t alpha^t series); grep for dilog / Li_2 / polylog / Spence across all source files: none.
  They evaluate the Neumann series numerically (convergence tables), not in closed form. Verified.
* Ruedenberg JCP 19, 1459 (1951) appears in their bibliography as cited; content not opened. Harris; Maslen &
  Trefry: not opened. UNVERIFIED -- the paper should not assert the closed form is new relative to these.

## Required caveats for the paper sentence

1. "Weight" = iterated-integral depth of functions of R (kernels e^{ws} ds/s); no claim about values at fixed R.
2. Lower bound is Theorem A only (linear, C(R) coefficients, any complex rates); no claim about dilogarithm or
   product representations (Theorem B is a sketch).
3. Novelty of the closed form relative to the classical exchange-integral literature (Ruedenberg, Harris,
   Maslen-Trefry) is unchecked.
4. The case list is the cases with order exactly 9 proved (seven, per Paper 58's caps), not "all".

**Proposed replacement** for "...which is not a combination of these seeds with rational coefficients (...);
its weight is not classified here.":

> ...which is a two-fold iterated integral, f0(R) = c sum_i k_i int_R^oo e^{lam_i s} Ei(mu_i s) ds/s
> + c ln((a+b)/(b-a)) [E1(2bR) - E1(2aR)], with an explicit closed form in which the dilogarithm constants of the
> six terms cancel. It is not a linear combination, with coefficients rational in R, of e^{wR}, e^{wR} ln R and
> e^{wR} Ei(nu R) for any finite sets of complex rates (its ln R coefficient at R = 0 has a derivative with an
> e^{2bR}/R component that no rational multiple of e^{2bR} differentiates to). As functions of R the
> unequal-exponent exchange integrals whose order 9 is proved therefore have iterated-integral weight two, not
> one; whether f0 admits a dilogarithm or product form, and what holds for values at fixed R, is not decided
> here. At a = b the weight-two part degenerates to weight one, consistent with Sugiura's closed form.

## Wrong answers a regression guard should reject

1. f0' = +FM1/R, or any LC sign flip (C(3/4, 2) becomes -0.225, not 0).
2. K with +pi^2/12 for mu > 0, or the complex Li2(p/mu) instead of its real part; ln^2 R with + sign.
3. Phi without the 1/(j + k) factor (fails T vs quadrature).
4. "C(a, b) = 0" checked only at rational pairs -- include an irrational ratio (1, sqrt 2).
5. "f0 has no ln R at R = 0" (Lam = 0) or Lam with the wrong sign: the log coefficient of f0 near 0 must be
   -Lam = R^2 + O(R^4).
6. a -> b: W2 diverging like 1/eps (Phi-sum not cancelling): (W2 - lim)/eps must stay bounded
   (about -0.287 at a = 1, R = 0.7).
7. Theorem-A guard: an exhibited rational rho with rho' + 2b rho = -c/R must be impossible (test the pole-order
   bookkeeping, e.g. a candidate rho = -c ln-free ansatz returns a residual).
