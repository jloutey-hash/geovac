# Sprint memo — H2 prolate CI: 99.1% is a conditioning artifact, not a wall (2026-09-16)

**Status:** PRODUCTIONIZED 2026-09-17 (v5.12.9). Module `geovac/prolate_recondition.py`
+ backing test `tests/test_paper12_recondition.py` (fire-tested `debug/firetest_p12_recondition.py`);
the module reproduces the PoC headline bit-for-bit ((3,3,1)=99.216%, (4,4,2)=99.711%,
(5,5,2)=99.767%). Paper 12's "99.1% cap" corrected to a conditioning wall. Item-2
(Gegenbauer) measured with a correction to its premise — see "Owed / done" below.

## Question

Paper 12 reports H2 at 99.1% of D_e with the azimuthal channels restored, and states
the basis is "conditioning-capped" there (`cond(S)=2e16` at `(3,3)`, `mu<=1`; a direct
`eigh(H,S)` returns -79 Ha). Is that 99.1% a **structural** ceiling of the
prolate-spheroidal natural geometry (the electron-electron cusp), or a **fixable**
artifact of the raw-monomial basis `u = xi^j eta^l (xi^2-1)^{mu/2}(1-eta^2)^{mu/2} e^{-a xi}`?

## Diagnosis (`debug/h2_recondition_probe.py`)

The monomial radial factor `xi^j` is the classic ill-conditioned set (a Hankel moment
problem). Measured on the REAL `geovac.prolate_general_m` overlap, sigma-only:

| (j,l) | N | cond(S) |
|---|---|---|
| (1,1) | 6 | 2.9e4 |
| (3,3) | 72 | **3.0e14** (matches the module's recorded 2.6e14) |
| (4,4) | 175 | 4.6e18 |
| (5,5) | 342 | 1.1e23 |

Direct solve returns energies BELOW the true floor (-74, -2470, -110393 Ha at (2,2)/(3,3)/(4,4))
— provably broken. Even canonical orthogonalization fails at (4,4) (-27 Ha). So the basis
**cannot be pushed past ~(3,3)** in float64: that is what pins the paper at 99.1%.

Re-basing the SAME span to Laguerre(xi) x Legendre(eta) (mpmath, dps=60): normalized
condition number stays **~89 at n=6** where the monomial is 2.9e10 — ~8 orders better,
and grows linearly not exponentially. Same span => identical exact-arithmetic energy;
only the conditioning changes.

**Verdict: 99.1% is a conditioning artifact of the monomial basis, not a structural
ceiling. FIXABLE.**

## Float64 re-basing (`debug/h2_recondition_relift.py`)

Build S,H,V (monomial, float64 via `prolate_general_m` + `neumann_vee_general_m`),
change basis to Laguerre x Legendre, solve. Recovered the climb — and at the *identical*
(3,3), mu<=1 span (dim 144, verified) got **99.22%** where the monomial solve gets only
99.09% (Paper 12's own "last two decimals are contamination"). Climb: (3,3)pi 99.22 ->
(4,4)pi 99.33 -> (3,3)+delta 99.56. BROKE at (5,5)/(4,4)+delta (float64 change-of-basis
amplifies the 1e-16 entry error by ||C||^2 when the monomial matrices carry cond ~1e26).

## High-precision pipeline (`debug/h2_recondition_hp.py`)

Every matrix ENTRY built in mpmath (dps=40) so the change of basis stays clean.
Validations: mpf one_body/V_ee match the float64 modules to ~1e-13; the **factored**
(Kronecker) change of basis matches the dense one to **4e-36**; reproduces every
overlap point exactly.

Two optimizations were needed and both landed:
- **Solve (fixed).** The raw `cond(S_orth)` ~1e12 was mostly harmless NORM SPREAD.
  Rescaling each orth function to unit norm (a diagonal congruence that leaves the
  generalized eigenvalues unchanged) drops it to **~3e7**, and a plain float64
  canonical-orthogonalization is then correct and robust (flat across every discard
  threshold) in **~1 s**. This is what had produced the wrong (5,5)+delta = 99.465
  (a variational violation) under the earlier per-block mpf-eigsy + float64-reduction solve.
- **Assembly (vectorized).** V_ee depends on the basis only through the SUMS of quantum
  numbers, so a small mpf tensor `F_{mui,muj}[p1,q1,p2,q2]` is precomputed once and the
  whole matrix gathered by vectorized indexing — exact (matches the loop bit-for-bit),
  ~30-80 s vs ~40 min. **float64 V FAILS** (-46 Ha): V carries the same dynamic range as
  S, so its float64 error is amplified by the change of basis — the assembly must stay mpf.

Remaining bottleneck at N>~2000: the mpf change-of-basis contractions (~40 min at N=1944).

## Results — high-precision, variationally rigorous ladder (H2, R=1.4011, alpha=1.0)

| basis | channels | D_e % | error |
|---|---|---|---|
| monomial ceiling (Paper 12) | sigma+pi | 99.09 | 1.58 mHa |
| re-based (3,3) | sigma+pi | 99.22 | 1.37 |
| re-based (4,4) | sigma+pi | 99.33 | 1.17 |
| re-based (5,5) | sigma+pi | 99.37 | 1.10 |
| re-based (3,3) | +delta | 99.58 | 0.73 |
| re-based (4,4) | +delta | 99.71 | 0.505 |
| **re-based (5,5)** | +delta | **99.767** | **0.406** |
| re-based (4,4) | +delta+phi | 99.767 | 0.406 |
| field (grid, same coords, TMR 2010) | all | 99.97 | 0.05 |

**Error cut 1.58 -> 0.41 mHa (~4x), monotone, every point variationally consistent
(each contains its predecessor and lies below it).** 0.41 mHa is comfortably inside
chemical accuracy (1.6 mHa).

## Physics finding: the angular channels saturate

pi buys ~+7 pp, delta buys ~+0.38 pp (at (4,4)), phi buys only **+0.056 pp** ((4,4)+delta
99.71 -> (4,4)+delta+phi 99.767). (4,4)+phi equals (5,5)+delta to 0.4 uHa. So the residual
~0.4 mHa is NOT more angular channels — it is radial completeness + the true e-e cusp,
which converges as the slow L^-3 partial-wave crawl. Reaching a literal 99.9% needs bases
~(8,8)+ (N>~5000, impractical here) OR explicit correlation (geminals put r12 in directly;
the He R12-CI PoC hit 0.8 mHa at ~6 functions).

## Owed / done (2026-09-17, v5.12.9)

1. **DONE.** Productionized into `geovac/prolate_recondition.py` (`recondition_energy(...)`,
   `basis='laguerre_legendre'|'gegenbauer'`) + `tests/test_paper12_recondition.py`
   (3 fast algebra-guards + 5 slow physics tests) + fire test `debug/firetest_p12_recondition.py`
   (4/4 guards fire; one vacuous congruence-guard was caught by the fire test and replaced
   with a norm-spread discriminator). The module reproduces the PoC bit-for-bit. Gate
   cleared; Paper 12 corrected in place (abstract, new Sec. "The monomial cap is
   conditioning, not a ceiling", conclusion; registry keys p12_rebased_*).

2. **MEASURED, with a correction to the premise.** The associated-Laguerre L_n^{(mu)} x
   Gegenbauer C_n^{(mu+1/2)} family gives the **IDENTICAL energy** as Laguerre x Legendre
   at equal (j,l) — same span (7 digits at every truncation), so it does **NOT** reach a
   given accuracy with *fewer functions* (the memo's owed-item-2 phrasing was wrong: a
   change of polynomial basis at equal degree spans the same space). Its real payoff is
   **conditioning**: cond(norm) 326x better at (3,3,1) (1.07e3 vs 3.49e5), ~1025x at
   (5,5)+delta (9.14e4 vs 9.37e10) — which keeps the downcast solve trustworthy at large
   truncation (Laguerre's 9.4e10 at (5,5)+delta is where the earlier per-block solve
   produced a spurious non-variational 99.465). It does **NOT** enable a fast float64
   pipeline: a float64 change of basis breaks for BOTH families at (3,3)/(4,4)/(5,5)
   (float64_relift probe), because the build-precision ceiling is the *monomial* matrices'
   dynamic range (cond 1.8e16 -> 4.1e26), not the target basis. Reaching a literal 99.9%
   still needs a much larger basis or explicit correlation; a radial re-basing cannot
   reach the e-e cusp (ledger 2026-08-23, elliptic-basis row). The genuine fast route =
   a DIRECT recurrence build in the orthogonal basis (never forming the monomial matrices)
   — identified, not built; a separate diagnostic->implementation sprint.

## Scope / caveats

Single alpha=1.0, single geometry R=1.4011, homonuclear H2. The heteronuclear and
many-electron diatomic cases are untouched. The [[polyatomic_state_of_play]] ranked path #1
("rebuild the Paper 12 prolate 2e CI at literature resolution — extend neumann_vee to
general m") is what this is; the general-m V_ee already landed in v5.12.7, this adds the
re-conditioning that lets it be pushed.

Drivers: `debug/h2_recondition_{probe,relift,hp}.py`.
