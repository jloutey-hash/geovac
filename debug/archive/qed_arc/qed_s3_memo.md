# One-Loop QED Vacuum Polarization on S^3

## What was computed

Seeley-DeWitt heat kernel coefficients a_0, a_1, a_2 for the squared
Dirac operator D^2 on the unit three-sphere S^3 (radius R=1), using the
Camporesi-Higuchi spectrum from `geovac/dirac_s3.py`.  The vacuum
polarization coefficient and QED beta function were extracted and
verified against the standard flat-space QED result.

### Seeley-DeWitt coefficients on unit S^3

For the 4-component Dirac spinor (dim_S = 4) on S^3 with:
- Vol(S^3) = 2 pi^2
- R_scalar = 6,  |Ric|^2 = 12,  |Riem|^2 = 12
- Lichnerowicz endomorphism E = R_scalar/4 = 3/2

The coefficients are (exact sympy):

| Coefficient | Value | Simplified |
|-------------|-------|------------|
| a_0 | (4pi)^{-3/2} * 4 * 2pi^2 | sqrt(pi) |
| a_1 | (4pi)^{-3/2} * 4 * 1 * 2pi^2 | sqrt(pi) |
| a_2 | (4pi)^{-3/2} * 4 * (45/360) * 2pi^2 | sqrt(pi)/8 |

The ratio a_0 : a_1 : a_2 = 1 : 1 : 1/8 on unit S^3.

R-scaling: a_0 ~ R^3, a_1 ~ R, a_2 ~ 1/R (verified symbolically).

### Vacuum polarization and beta function

The coefficient of F_{mu nu} F^{mu nu} in the one-loop effective action
for a single Dirac fermion minimally coupled to U(1) is:

    Pi = 1/(48 pi^2)

This gives the one-loop QED beta function:

    beta(alpha) = 2 alpha^2 / (3 pi)

Both are verified symbolically to exact equality (sympy).

## Connection to flat-space QED

The Seeley-DeWitt coefficients are LOCAL invariants: they depend only on
the point-wise geometry and the operator structure, not on the global
topology.  Therefore a_2 (which encodes the vacuum polarization) is
identical on S^3 and on flat R^3.  The beta function is a UV quantity
that is insensitive to the global curvature in the R -> infinity limit.

The S^3 computation adds curvature-dependent terms to a_0 and a_1 (via
Vol(S^3) and R_scalar), but these do not affect charge renormalization.
The a_2 integrand on unit S^3 evaluates to the rational number 45; the
gauge-field piece (1/(48 pi^2)) sits on top of this curvature background.

## Transcendental classification (Paper 18 taxonomy)

| Quantity | Content | Tier |
|----------|---------|------|
| a_0 = sqrt(pi) | Vol(S^3) = 2pi^2 | Calibration pi (2nd-order Riemannian) |
| a_1 = sqrt(pi) | R_scalar = 6 (rational) x Vol | Calibration pi |
| a_2 curvature | Rational combo of R_sc, Ric, Riem x Vol | Calibration pi |
| Vacuum pol. | 1/(48 pi^2) | Calibration pi (2nd-order, gauge) |
| beta(alpha) | 2 alpha^2 / (3 pi) | Calibration pi |
| zeta_{D^2}(s) | T9 theorem: polynomial in pi^2 at each integer s | Calibration pi ONLY |
| Odd-zeta at 1-loop | ABSENT | Structural theorem (not numerical) |
| Odd-zeta at 2-loop | EXPECTED (zeta(3)) | Track D3 weight-m subchannel |

### Key structural result

The T9 theorem (Tier 3 sprint) guarantees:

    zeta_{D^2}(s) = 2^{2s-1} * [lambda(2s-2) - lambda(2s)]

where lambda(2k) = (1 - 2^{-2k}) * zeta_R(2k) = rational * pi^{2k}.

Since zeta_R(2k) is always rational * pi^{2k} (via Bernoulli numbers),
the squared Dirac spectral zeta at every integer s is a polynomial in
pi^2 with rational coefficients.  No zeta(3), zeta(5), or any odd-zeta
value appears.

Verified symbolically at s = 2, 3, 4, 5:
- zeta_{D^2}(2) = pi^2 - pi^4/12
- zeta_{D^2}(3) = pi^4/3 - pi^6/30
- zeta_{D^2}(4) = 2 pi^6/15 - 17 pi^8/1260
- zeta_{D^2}(5) = pi^even terms only (sympy evaluates all zeta(even) automatically)

Numerical cross-check: direct spectral sums at s = 2, 3 converge to the
T9 values (0.2% at n_max = 500 for s = 2; < 10^{-6} for s = 3).

## What comes next

### Two-loop: where zeta(3) should appear

Track D3 (Tier 1) showed that the Dirac degeneracy g_n^Dirac = 2(n+1)(n+2)
generates a Dirichlet series D_{g^Dirac}(4) = 2 zeta(2) + 2 zeta(3) at
the packing exponent s = d_max = 4.  The zeta(3) comes from the weight-m
subchannel (the linear piece of g_n).

At two loops, the QED effective action involves products of propagators
that probe the m-resolved structure of the Dirac degeneracy.  The Euler-
Maclaurin or heat-kernel expansion of such products will access D3's
zeta(3) subchannel.  This is consistent with the well-known appearance
of zeta(3) in the two-loop QED beta function (Rosner 1967, Laporta-
Remiddi 1996).

Paper 18's operator-order x bundle grid predicts:
- 1st-order operators (Dirac D): generate odd-zeta content (T9 corollary)
- 2nd-order operators (D^2, Laplace-Beltrami): generate even-zeta only

The one-loop effective action involves D^2 (second-order), confirming
the even-zeta-only result.  The two-loop effective action involves the
PRODUCT of two D-propagators, which effectively probes D at first order
in each propagator -- opening the odd-zeta channel.

### Future directions

1. Compute the two-loop vacuum energy numerically from the S^3 spectral
   data and verify the zeta(3) coefficient against Rosner/Laporta-Remiddi.
2. Extend to massive QED: the spectral_zeta_massive function provides
   the mass-dependent spectral zeta; the m -> 0 limit should reproduce
   the massless result smoothly.
3. Investigate whether the S^3 curvature modifies the coefficient of
   zeta(3) at two loops compared to flat space (likely no, by locality
   of the a_3 coefficient, but worth verifying).

## Files

- `geovac/qed_vacuum_polarization.py` -- implementation (6 public functions)
- `tests/test_qed_vacuum_polarization.py` -- 33 tests, all passing
- `debug/qed_s3_memo.md` -- this memo
