# Pre-registered prediction: which large-L solutions the two-electron integral uses (B2X)

**Registered:** 2026-10-29, before any W1 coefficient is evaluated numerically. PI-directed ("Let's work that open
question").
**Commit status:** committed on its own, before the scored run.
**Rests on:** B2 (`2026-10-28_two_electron_bottom_up.md`, BORDERLINE). Its derived recurrence for case W1 has, at
lam = 17/25 (x = x* = 457/425), the exponents beta = -4, -5, -6, 3, 1, and decaying alternating roots
-1875/2057 (|lam| = 0.912, beta = -1/2) and -51/121. B2 left open whether the true coefficients use beta = 3 and 1.

## Why a direct test settles it

The true coefficients F_L are fixed sums of exact rationals times the 31 constants (`debug/ao_exact.py`). Evaluated at
high precision, their large-L form is lam^L L^beta for the slowest solution present. If beta = 3 or 1 is present with
any nonzero weight, it outgrows beta = -4 and dominates.

**Physical expectation.** Along the bond-angle line, |Y - Z|^2 = D1^2 + D2^2 - 2 D1 D2 x, which is linear in x and
vanishes at x*. So f(x) is a function of |Y - Z|^2, and the integral stays finite when the two hydrogen centres meet.
A finite value rules out negative alpha at x*. beta = 3 gives alpha = -7/2 and beta = 1 gives alpha = -3/2, so both
should be absent.

## Method

`debug/b2_large_l.py W1 301`. This evaluates F_L for L <= 301 exactly at two working precisions, dps = 6 LMAX + 200 and
1.5 times that, and reports the number of digits that agree. It also reports r_L = F_{L+1}/F_L and
beta_eff(L) = ln(G_{L+1}/G_L)/ln(1 + 1/L), where G_L = F_L (25/17)^L.

**Control, measured before registering** (the one-electron case N1 of NA; same script, lam = 9/13; NA's
exponent 3/2 means beta = -2):

| L | 100 | 200 | 300 |
|:--|:--|:--|:--|
| beta_eff | -2.062 | -2.032 | -2.022 |

At least 691 digits agree at every L. The drift is about 6/L.

## Predictions

- **X0 (precision, a validity check).** At least 20 digits agree at every L <= 300. If this fails, the run is void and
  is repeated at higher precision; no science changes.
- **X1 (no alternating tail).** r_L > 0 for 100 <= L <= 300, so the root -1875/2057 and the other negative roots are
  absent. Probability 0.8.
- **X2 (the open question).** |beta_eff(300) + 4| < 0.5. That means beta = 3 and 1 are absent, and the slowest term is
  the contact term (x* - x)^(7/2). The bound is about 20 times the control's drift at L = 300. Probability 0.85.
  This covers both the finiteness argument and the chance that the contact weight vanishes, in which case beta_eff
  would go to -5.

**Reported, not scored:** beta_eff and r_L at L = 100, 200 and 300.

**Gate.** GO if X1 and X2 hold. Then the exponents 3 and 1 belong to other solutions of B2's recurrence, as in NA, and
B2-3's physical content holds for the integral. STOP if beta_eff(300) lies within 0.5 of 3 or 1. BORDERLINE otherwise.

## Caps

This covers one geometry and 1s functions. It tests the true integral, not the generic sequence.
