# Pre-registered prediction: hydrogen p functions in the bond-angle series (PP)

**Registered:** 2026-10-29, before any p-function coefficient in this file is computed. PI-directed ("Yes please",
on drafting a pre-registration for hydrogen p functions as the next step).
**Commit status:** to be committed on its own, before any computation it scores.
**Rests on:** Paper 59 sec:scope (the cusp-cusp term and the contact reading); records
`docs/predictions/2026-10-27_water_zero_contact_fresh.md` (WA8: a cusp-free H function moves the equal-length term
from (1-x)^(7/2) to (1-x)^(9/2), and two of them to (1-x)^(11/2)) and `docs/predictions/2026-10-29_two_electron_large_l.md`
(B2X: 1s coefficients decay as rho^L L^-4). Prior art: Matuzas 2026 (doi:10.5281/zenodo.23448184) proves the 1s germ
exponents and decay laws and states that higher angular momentum is outside its scope.

## Why

Real water bases put p (polarization) functions on the hydrogens. The contact reading says the tail is set by the
short-distance form of the hydrogen functions. A p function has no 1s-type contact, so the reading predicts a
different tail. This is a prediction of the reading that the 1s work cannot test.

## The derivation (done before registering)

Integral (XY|XZ) as before: X the oxygen 1s (exponent zX) at the origin, Y at distance D1, Z at D2, x = cos(angle).
A p_sigma function on Y, pointing along the bond, is

    ((r - Y).yhat) e^{-z |r - Y|} = yhat . grad_Y f(|r - Y|),   f(s) = (s/z + 1/z^2) e^{-z s},

checked by differentiating: f'(s) = -s e^{-z s}. And yhat . grad_Y is d/dD1 at fixed angle. f has no linear term in s
(it is cusp-free; its first odd term is (z/3) s^3), and s e^{-zs} = -d/dz e^{-zs}. So every case is a derivative of
the 1s integral with separate hydrogen exponents z1 (on Y) and z2 (on Z), whose exact Legendre coefficients q_L are
computed by `debug/pp_exact.py`:

    O_j = 1/z_j^2 - (1/z_j) d/dz_j   (turns a 1s on that hydrogen into f)
    PA  (p_sigma on Y, 1s on Z):          F_L = d/dD1 O_1 q_L
    PB  (p_sigma on Y and on Z):          F_L = d/dD1 d/dD2 O_1 O_2 q_L
    PC  (p_pi out of plane on Y and Z):   I = (1/(D1 D2)) d/dx [O_1 O_2 I],  coefficients from O_1 O_2 q_L by the
                                          Legendre derivative rule, F_L = (2L+1)/(D1 D2) sum_{k = L+1, L+3, ...} c_k

(One out-of-plane p_pi alone gives zero by reflection symmetry, so it is not a case.)

Germ exponents. Near the coalescence point, with kappa = |Y - Z|^2 = D1^2 + D2^2 - 2 D1 D2 x, the 1s-f integral is
B kappa^(9/2) + analytic and f-f is B kappa^(11/2) + analytic (WA8 at equal lengths; Fourier counting in general).
Then, with d kappa/dD1 = 2(D1 - D2 x), d kappa/dD2 = 2(D2 - D1 x), d^2 kappa/dD1 dD2 = -2x, d kappa/dx = -2 D1 D2:

| Case | Unequal lengths: at x* the factors dkappa/dD are nonzero | Equal lengths: dkappa/dD1 = dkappa/dD2 = kappa/D |
|:--|:--|:--|
| PA | kappa^(7/2): F_L ~ rho^L L^-4 | kappa^(9/2): F_L ~ L^-10 |
| PB | kappa^(7/2): F_L ~ rho^L L^-4 | kappa^(9/2) (from the -2x term): F_L ~ L^-10 |
| PC | kappa^(9/2): F_L ~ rho^L L^-5 | not scored (see Caps) |

Rule used: a germ kappa^a gives L^(-a-1/2) rho^L at unequal lengths and L^(-2a-1) at equal lengths (Matuzas
Thm 4.1; B2X and WA8 confirm it for a = 7/2, 9/2, 11/2). For comparison, 1s on both hydrogens gives -4 and -8.

## What I know before comparing (disclosed)

- No p-function coefficient has been computed for any case in this file.
- Controls run on non-scored cases (`debug/pp_controls.py`, log `debug/data/pp_controls.log`):
  - C1: `pp_exact.coefficients2` with z1 = z2 reproduces `ao_exact.coefficients` term by term (L = 0..12, control
    geometry W1). PASS.
  - C2: with z1 = 1.1, z2 = 1.6, zX = 1.4, D1 = D2 = 1.809, it agrees with the independent float quadrature (the
    `F_combo` channel formula of `debug/wa7_score.py`) to 1e-16..3e-12 for L = 0..8. PASS.
  - C3: one d/dD1 at the W1 geometry (1s on both H) to L = 150, steps 1e-30 and 1e-31: the two derivatives agree
    to 7e-61 (L = 10) .. 1.1e-57 (L = 150). Cost about 1000 s per exact run at L = 150 with perturbed rationals.
- Extrapolation control, measured on the known 1s data (`debug/data/b2x_W1.json`, `b2x_N1.json`): fitting
  beta_eff(L) = beta + c/L + d/L^2 on L = 100..150 gives beta to 1.6e-4 (W1, true -4) and 3.6e-5 (N1, true -2). The raw
  beta_eff at L = 150 is off by 0.15 and 0.04, so the raw value is not scored.

## Cases (fresh; not computed before)

| Case | D1 | D2 | rho | zX | z (both H) |
|:--|:--|:--|:--|:--|:--|
| U | 9/5 | 12/5 | 3/4 | 13/10 | 6/5 |
| E | 19/10 | 19/10 | 1 | 13/10 | 6/5 |

Scored items: PA-U, PB-U, PC-U, PA-E, PB-E. Derivatives by central differences in exact rationals with step h = 1e-30,
coefficients evaluated at 6 Lmax + 400 digits; Lmax = 150.

## Predictions

**PP-0 (controls; run FIRST).**
- (a) PA-U recomputed with step 1e-31 agrees with step 1e-30 to 1e-20 relative for every L <= 150 (C3 measured
  1e-57 on a control; the bound is loose on purpose). If (a) fails, every item is void. Only PA-U is rechecked,
  because a full second set doubles a run of about 95 exact runs.
- (b) The 1s case (z1 = z2 = 6/5, no derivatives) at geometry U extrapolates to beta within 0.02 of -4, and at
  geometry E within 0.02 of -8. If (b) fails at a geometry, both items at that geometry are void.

**PP1 (decay powers).** For every non-void item, beta from the fit beta_eff(L) = beta + c/L + d/L^2 on L = 100..150
(G_L = F_L / rho^L) is within 0.02 of the table: PA-U -4, PB-U -4, PC-U -5, PA-E -10, PB-E -10.
- Falsifier: any non-void item outside 0.02 of its predicted beta.
- The bound 0.02 is about 100 times the measured 1s extrapolation error and 1/50 of the gap to the nearest integer
  alternative.
- Probability 0.75.

**PP2 (the shift is the p function's).** At geometry E, PA and PB both extrapolate at least 1.5 below the 1s value
from PP-0(b), i.e. beta < -9.5.
- Falsifier: either is above -9.5.
- Probability 0.8 (it is implied by PP1 at E but is stated so that a near miss of PP1 is still scored on the shift).

## Decision gate

- **GO:** PP-0(b) passes at both geometries, PP-0(a) passes, and PP1 holds.
- **BORDERLINE:** PP1 fails on exactly one item while PP2 holds.
- **STOP:** PP-0(a) fails, PP-0(b) fails at either geometry, PP2 fails, or PP1 fails on two or more items.

## Caps

- No change to `geovac/` or `tests/`. 1s on the oxygen; p functions only on the hydrogens; one exponent per hydrogen.
- PC at equal lengths is not scored: its derivative-rule sum converges only like k^-11 there, so truncation at Lmax
  would contaminate the coefficients near Lmax.
- The sizes of the leading terms are not predicted, only their powers.

## Disposition

(Appended after comparing. The predictions above are not edited.)

## Addendum 1 (appended 2026-10-29, after the registration commit and before any scored quantity is formed; the predictions above are not edited)

- **A defect in the PC plan.** PC's coefficients come from the Legendre derivative rule,
  F_L = (2L+1)/(D1 D2) sum_{k = L+1, L+3, ...} c_k. With c_k computed only to 150, F_L near L = 150 loses most of its
  sum (at L = 150 the sum is empty), so the fit window 100..150 cannot be used as registered. Found while writing the
  scorer; no PC quantity had been formed.
- **Fix.** PC's nine inputs O_1 O_2 q are computed to L = 230 (points named with suffix `_L230`, same step, digits
  6*230 + 400). The omitted part of each sum is about rho^(230 - L) of it, 1e-10 or less on L <= 150 at rho = 3/4.
  The PC fit window, bound and prediction are unchanged. The nine PC points already launched at L = 150 are not used
  for PC; the one with no offsets stays the PP-0(b) 1s control at U, as registered.

## Disposition, scored 2026-11-01 (appended after comparing; the predictions and Addendum 1 are not edited)

All 109 points on disk; `python debug/pp_score.py` (PC from the L = 230 points, per Addendum 1). Output
`debug/data/pp_score.log`.

| Item | beta (fit L = 100..150) | Predicted | Off by | Verdict |
|:--|--:|--:|--:|:--|
| PP-0(a) step 1e-31 vs 1e-30 | max rel dev 9.2e-58 | < 1e-20 | | PASS |
| PP-0(b) 1s at U | -3.99965 | -4 | 0.00035 | PASS |
| PP-0(b) 1s at E | -8.00004 | -8 | 0.00004 | PASS |
| PA-U | -3.99964 | -4 | 0.00036 | HELD |
| PB-U | -3.99797 | -4 | 0.00203 | HELD |
| PC-U | -4.99892 | -5 | 0.00108 | HELD |
| PA-E | -10.00006 | -10 | 0.00006 | HELD |
| PB-E | -9.99941 | -10 | 0.00059 | HELD |
| PP2 (PA-E, PB-E < -9.5) | -10.0001, -9.9994 | < -9.5 | | HELD |

No item changed sign in its raw beta_eff. Largest miss 0.002 (PB-U), a tenth of the bound.

**Gate: GO.** PP-0 passes at both geometries and PP1 holds on all five items. A p function on the hydrogens leaves the
decay power at unequal lengths unchanged (-4, PC -5) and lowers it by two at equal lengths (-8 to -10), as derived.
