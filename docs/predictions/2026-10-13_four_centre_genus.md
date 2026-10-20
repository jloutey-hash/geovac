# Pre-registered prediction: the genus ladder stops at one for electron repulsion

**Registered:** 2026-10-13, before any four-centre computation (PI-directed).
**Owner of the rule being extended:** Paper 59 (`sec:momentum`, `sec:genus`, `sec:scope`).

## The rule as Paper 59 states it

In momentum space each nucleus enters an electron-repulsion integral only as a
phase. Each two-centre charge density contributes one dispersion factor
Delta = sqrt(c k^2 + m^2) (c, m set by a Feynman parameter and the Slater
exponents). One factor gives genus 0 ({E1, ln, gamma}); two factors with
different scales give the genus-1 curve y^2 = (c1 k^2 + m1^2)(c2 k^2 + m2^2).
Paper 59 measured this for (XY|XZ) only (three centres, shared centre X).

## Predictions

**P1 (main).** The general four-centre integral (XY|ZW), four distinct
centres, reduces to the same form: a two-parameter Feynman integral over one
radial k-integral whose only square roots are Delta1(s) and Delta2(t). Its
fibre curve is the same genus-1 family. **No genus 2 at four centres.** More
generally: a two-electron integral with a rational kernel (Coulomb, Yukawa,
Slater geminal) has genus <= 1 whatever the number of centres, because it
couples exactly two densities.

**P2 (geometry is coefficients).** Positions, bond lengths and bond angles enter
only through the phase vector W(s,t) (algebraically, via |W|) and through the
pair distances D1, D2 in e^{-D Delta}. Bending water changes values, never the
curve.

**P3 (hunch, low confidence, NOT tested here).** The first new rung (genus 2)
appears where three distinct densities share one momentum integration,
i.e. in some three-electron integrals. The loop structure there differs (two
momenta), so the genus language may not transfer directly.

## Falsifiers (any one breaks P1)

1. The four-centre reduction requires a third independent square root.
2. The reduced (s, t, k) formula fails to reproduce an independent evaluation
   (Gaussian-fit McMurchie-Davidson) of the same integral.
3. The fibre at fixed (s, t) is not the Paper 59 curve.

## Prior art checked before registering

Avery & Avery, *Hyperspherical Harmonics and Their Physical Applications*
(WS 2017), Appendix B "Interelectron repulsion integrals with STO's": the
four-centre STO integral written as 4 pi int dp p^-2 rho1^t(-p) rho2^t(p), the
same momentum-space form, evaluated by Coulomb-Sturmian expansion. No genus
or elliptic reading of it appears in the book (searched: "elliptic", "genus";
the only "elliptic" hit is ch. 8's elliptic-coordinate Sturmians, a
two-focus/prolate false friend).

## Novelty classification (appended 2026-10-19)

Classes per `docs/predictions/README.md`, set after the fact by a literature audit (`debug/prediction_novelty_audit_memo.md`; this record predates the novelty field). FINDING = held, not found in the literature; CHECK = held, already known; MISS = failed.

- P1: FINDING, low weight (the momentum form is the Averys')
- P2: FINDING, routine
- P3: UNSCORED (superseded by P3a-P3e)
