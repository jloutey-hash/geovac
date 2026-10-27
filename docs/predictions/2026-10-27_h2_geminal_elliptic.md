# Pre-registered prediction: the H2 Slater-geminal master integral and the elliptic curve of Paper 59 (HG0-HG5)

**Registered:** 2026-10-27, before any computation that scores a prediction in this file. PI-directed ("Yes
please. Let's do that H2 one then the water one after"; re-aimed after the prior-art check, "Yes please").
**Commit status:** committed on its own, before any computation in this file.
**Rests on:** Paper 59 `eq:pf` and its irreducibility statement (`sec:obstruction`).

## Prior art (read before registering)

- **Pachucki, arXiv:0905.1635 (2009).** The two-centre two-electron master integral with Coulomb
  r12 (no r12 exponent) obeys a second-order equation in the bond length r. Homogeneous solutions I0(pr), K0(pr);
  source built from Ei and exp. Already cited in Paper 59.
- **Lesiuk and Moszynski, PRA 86, 052513 (2012), arXiv:1209.0839** (read in full, `debug/lit_scan/dl_lesiuk2012/`,
  untrusted folder, text extracted with pdftotext). Master integral
  f(r) = r Int d3r1 d3r2 e^{-u3 r1A - u2 r1B - w2 r2A - w3 r2B - w1 r12} / (16 pi^2 r1A r1B r2A r2B r12).
  - Its Laplace transform g(u1) is the six-mass tetrahedron (three-electron one-centre generating integral, Fromm and
    Hill; Harris), Eqs. 11-13.
  - g obeys a FIRST-order equation in u1 (Eq. 35): Delta(u1) g' + (1/2) Delta'(u1) g + P = 0, with
    Delta = w1^2 u1^4 + s1 u1^2 + s2 (Eq. 31). s1 and s2 are polynomials in the five exponents (Eqs. 32-33).
  - Inverse Laplace gives the fourth-order equation (Eq. 54, homogeneous part):
    w1^2 r f'''' + 2 w1^2 f''' + s1 r f'' + s1 f' + s2 r f = 0.
  - Its solutions are named "hyper-Bessel functions", defined by series. The authors say they could not express them
    in known functions. At w1 = 0 they reduce to I0, K0 (Pachucki).
  - Special cases: s2 = 0 holds for the symmetric James-Coolidge pattern (u2 = u3 = w2 = w3) and the symmetric
    exchange pattern (u3 = w2, u2 = w3). There the solutions are 1F2 and 2F3 functions. s1 = 0 and
    s1^2 = 4 w1^2 s2 are also treated.
- **Pachucki, arXiv:1209.1258 (2012)** also derives a fourth-order equation for exponentially correlated diatomic
  integrals. Not read in full.
- **Not found:** a short search found no work that identifies the hyper-Bessel functions with periods of an elliptic
  curve, or relates Eq. 54 to Paper 59's three-centre operator. That identification is the question here.

## What I know before computing (disclosed)

- **By eye.** Paper 59 `eq:pf` is rho D L'''' + 2 rho L''' + (1 - 2 rho) D L'' + (1 - 2 rho) L' - (1 - rho) D L = 0.
  It has the same form as LM Eq. 54: a r f'''' + 2a f''' + b r f'' + b f' + c r f = 0. Under r -> lambda r only the
  ratio b^2/(a c) is invariant, so a one-parameter match is expected. Not checked symbolically.
- **By hand, from the pdftotext transcription of Eq. 33:** s2 = 0 on the two symmetric patterns above (it agrees with
  LM's statement). On the (AA|BB) pattern u3 = w3 = x, u2 = w2 = 0, s2 = x^4 (2x^2 - w1^2), which is nonzero.
  The transcription may garble superscripts.
- **Why the curve is elliptic (reasoning, not checked).** The homogeneous g is Delta^{-1/2}. So the r-space
  homogeneous solutions are Laplace integrals of du/sqrt(quartic). That is the form of Paper 59 `eq:laplace`.

## The questions

- **Question O (operator).** Is LM's geminal operator exactly Paper 59's three-centre Picard-Fuchs operator, up to
  a rescaling of r?
- **Question G (genus on physical H2).** Which parameter patterns that real H2 bases produce sit on the genus-zero
  locus, and which are genus one?

## Predictions

**HG0 (controls; run FIRST; HG1-HG4 void if any fails).**
- (a) Pachucki limit: at w1 = 0, Eq. 54 (as transcribed) reduces to r f'' + f' - p^2 r f = 0 with p^2 = -s2/s1.
- (b) Paper 59 `eq:pf` is annihilated (residual below 1e-20) by its own Laplace integral, `eq:laplace`, at
  rho = 1/3 and D = 1, 2 (re-check of the paper's statement on this machine).
- (c) The transcribed s1, s2 give s2 = 0 at u2 = u3 = w2 = w3 = 1, w1 = 3/2, and at u3 = w2 = 1, u2 = w3 = 2,
  w1 = 1 (LM's stated special cases).

**HG1 (operator identity).** For s2 != 0 and s1^2 != 4 w1^2 s2, let t1^2, t2^2 be the two roots of
w1^2 t^4 + s1 t^2 + s2 = 0. Then, with r = D / t1 and rho = t1^2 / (t1^2 - t2^2), LM's Eq. 54 operator equals a
nonzero constant times Paper 59's `eq:pf` operator. Exact symbolic identity in the five exponents, with t1 a symbol
constrained by its quartic.
- Falsifier: any nonzero remainder after normalizing the leading coefficient.
- Probability 0.9.

**HG2 (independent derivation of Delta).** Delta(u1) is derived without LM's Eqs. 32-33, as the polynomial in the
first-order u1-equation obeyed by the tetrahedron G(1,1,1,1,1,1). Route: integration by parts on the momentum-space
tetrahedron, or the leading coefficient read off a fitted first-order recurrence of g at exact parameter values. It
equals LM's w1^2 u1^4 + s1 u1^2 + s2 up to an overall constant.
- Falsifier: a different quartic, or no first-order equation.
- Probability 0.8 (most of the risk is transcription).
- If HG2 fails, HG1 and HG4 are scored with the independently derived Delta, and the mismatch is reported.

**HG3 (the master integral itself, independent route).** Evaluate f(r) by a route that does not use LM's
r-equation or their w1-series (Eq. 183). Examples: an expansion of e^{-w1 r12}/r12 in modified spherical Bessel
functions about centre A, or an inverse Laplace transform of a separately computed g(u1). It reproduces LM
Table I (u2 = w2 = 1, u3 = w3 = 2, w1 = 1, r = 2): 2.39121170158262983284e-4, and Table III (same, w1 = 2):
1.1754079...e-4 at the digits LM's converged entry gives (the extraction truncates it; use their last row).
- Falsifier: disagreement beyond 1e-15 relative (Table I).
- Probability 0.8.

**HG4 (genus classification of physical H2 patterns).** Genus zero means s2 = 0 (one root pair at zero, the analogue
of Paper 59's equal-scale collapse rho = 1) or s1^2 = 4 w1^2 s2 (double roots). Otherwise genus one.
- (a) Symmetric James-Coolidge u2 = u3 = w2 = w3: s2 = 0 identically. Probability 0.95.
- (b) Symmetric exchange / Kolos-Wolniewicz u3 = w2, u2 = w3: s2 = 0 identically. Probability 0.95.
- (c) Pattern u3 = w3 = x, u2 = w2 = y (gerade direct-like): genus one for generic x, y, w1. Probability 0.8.
- (d) Single-zeta 1s H2 classes, with the master's exponents read off each orbital product:
  (AA|AA) [u3 = w2 = 2z, u2 = w3 = 0] genus zero; (AB|AB) [all z] genus zero; (AA|BB) [u3 = w3 = 2z,
  u2 = w2 = 0] genus one; (AA|AB) [u3 = 2z, u2 = 0, w2 = w3 = z] genus one. Probability 0.6 for all four together.
- (e) Fully unequal exponents (heteronuclear, e.g. HeH+ or LiH valence): genus one generically. Probability 0.9.
- Falsifier per clause: s2 and s1^2 - 4 w1^2 s2 computed exactly do not match the stated verdict.

**HG5 (reported, not scored).** For the genus-one patterns, the dictionary value rho at a few physical points (H2 near
equilibrium, zeta = 1, w1 in {1/2, 1, 3/2}). Whether any lands at rho = 1/2 (the CM point of Paper 59) or near the
collapse rho = 1.
- Caveat stated in advance: physical integrals are parameter derivatives of the master, so a genus-one master does
  not by itself make a physical integral genus one. Paper 59 already shows one-centre-density (AA|BB) integrals
  with a geminal kernel are elementary. HG4(d) classifies the master's operator, not the physical integral.

## Decision gate

- **GO:** HG0 passes, HG1 holds. LM's hyper-Bessel functions are then Paper 59's elliptic deformation of K0, and HG4
  says where physical H2 sits.
- **BORDERLINE:** HG1 holds only after a change of variable other than r -> lambda r, or only with an extra gauge
  factor. Same curve, but the identification is weaker than stated.
- **STOP:** HG0 fails, or HG1 fails with no equivalence found.

## Caps

- Symbolic work (HG0, HG1, HG2, HG4): at most 1 h each. If unfinished, the item is not scored.
- HG3: at most 2 h wall, at most 8 processes. If over, HG3 is void and the gate still stands on HG1.
- Whole run: at most 4 h wall.

## Novelty status (before computing)

- Fourth-order equation, hyper-Bessel functions, s2 = 0 special cases: KNOWN (Lesiuk-Moszynski 2012).
- Tetrahedron form of the Laplace transform: KNOWN (LM; Fromm-Hill; Puchalski-Kedziera-Pachucki).
- Paper 59 `eq:pf`: this corpus.
- That the hyper-Bessel functions are Laplace periods of an elliptic curve, the same operator as the three-centre
  ERI's: NOT FOUND in a short search. Treat as OPEN, not as new, until a fuller search.

## Disposition

(Appended after computing. The predictions above are not edited.)
