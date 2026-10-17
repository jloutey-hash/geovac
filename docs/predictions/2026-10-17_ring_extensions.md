# Pre-registered prediction: the correlated triangle, and a ring of four

**Registered:** 2026-10-17, before any computation in this file (PI-directed: "let's use some sub
agents to run as many of these next steps as is reasonable", after the ring result of
`2026-10-17_density_loop.md`).
**Extends:** `2026-10-17_density_loop.md` (PT1-PT4, Result, and "The chain control, settled");
Paper 59 `sec:density_ring`, `sec:four_density`, `sec:three_electron`; seam rows J9-J11.
**Commit status:** to be committed on its own, before any computation, if the PI authorizes it.
If computation starts uncommitted, the Result section must say so.

## Why these two

The Coulomb triangle sits on the elliptic rung, one below the chain of three, because a closed
loop of Coulomb kernels integrates out exactly in three dimensions (the star-triangle relation).
Two questions follow. (1) The LiH explicitly correlated engine meets the triangle with geminal
kernels, <f_ik f_jk / r_ij>, whose loop is massive and not "unique", so the exact step is gone.
Does that triangle still sit low? (2) Is "a ring sits one rung below" a rule, or a fact about
three? The smallest test is a ring of four Coulomb kernels.

## Conventions (as in PT1-PT3)

Fingerprint fibre: densities spherical at zero separation, rho~_v(k) = Delta_v(k)^(-3),
Delta_v = sqrt(c_v k^2 + 1), a_v = c_v^(-1/2); in position space rho_v is proportional to
K0(r / a_v). theta = t d/dt. Rungs: elementary or algebraic < elliptic (Legendre, F_2) <
K3 (F_3's operator 4 theta^3 - t^2 (theta+1)^3) < the four-density rung (AESZ #58).
"Sits on rung X" means: the minimal annihilator of the function, checked on the function
itself, has composition factors no higher than X, and at least one factor equivalent (rational
gauge, rational pullback allowed and stated) to X's operator. Lesson of 2026-10-17: an operator
from a Mellin recurrence is an upper bound; its minimality and any source term must be checked
on the function (the triangle's -3 t^2; the chain's spurious N2).

## QG: the correlated triangle

Object: Y(t) = <rho_i rho_j rho_k f(r_ik) f(r_jk) / r_ij> with the exponential geminal
f(r) = exp(-gamma r), at gamma = 1 and gamma = 3/2, scales (a_i, a_j, a_k) = (1, 1, t), t on the
apex k (the density carrying both geminals). Secondary: t on i. Route allowed: f(r) =
-d/dgamma [exp(-gamma r)/r], so Y is a mixed second gamma-derivative of the Yukawa-Yukawa-Coulomb
triangle with two masses set equal; the 3-D one-loop triangle with two massive and one massless
propagator is then the loop integral.

**QG1 (the massive loop is elementary, symbolic + measured).** The 3-D one-loop triangle with
masses (gamma, gamma, 0) is an elementary function of the three external momenta (inverse
tangents and square roots), in closed form; check it against direct 3-D integration at five
random external triples, relative < 1e-8. Prior art to verify, not assumed: Nickel, J. Math.
Phys. 19, 542 (1978).

**QG2 (main: the geminal does not raise the rung).** Y(t) sits on the elliptic rung or lower,
at both gamma: no factor equivalent to F_3's K3 operator and nothing higher. Reason: the loop
is still elementary, so after it the integrand lives on the same two momenta as the Coulomb
triangle, with an inverse tangent where there was a rational kernel; an inverse tangent adds an
integration, which can add weight, not geometry. Confidence: moderate. Also predicted: new
singular points that move with gamma (for example where c_k gamma^2 hits a fixed value), and
non-split attachment of the elliptic piece (an inverse tangent is an integral).

Outcomes: (G-low) elliptic or lower, as predicted. (G-chain) a K3 factor: the geminal triangle
sits with the chain, the mass undoes the drop. (G-climb) a factor above K3.

## QR: a ring of four

Object: S(t) = <rho_a rho_b rho_c rho_d / (r_ab r_bc r_cd r_da)>, four densities in a square of
Coulomb kernels, scales (1, 1, 1, t). Position-space route allowed: with K_v(x, y) =
int rho_v(z) / (|x - z| |z - y|) d^3z, S = int rho_a(x) rho_c(y) K_b(x, y) K_d(x, y), a
three-dimensional integral over |x|, |y| and the angle. The 3-D one-loop box with massless
propagators is not unique (indices sum to 4, not 3), so the triangle's exact step does not apply.

**QR1 (main, low confidence: a ring of n sits at the rung of a tree of n - 1).** S(t) has a
composition factor equivalent to F_3's K3 operator and nothing higher: the square sits on the
K3 rung, one below a tree of four. Reason: the only data point, the triangle, sits one below the
tree of its size; the guess is that closing a ring costs exactly one rung. This is a guess
from one point, registered so that it can fail.

Outcomes: (R-n-1) K3, as predicted. (R-ell) elliptic or lower: rings collapse to the bottom
whatever their size. (R-tree) the four-density rung: a ring of four sits with the tree of four.
(R-other) a factor not on the ladder.

**QR2 (one-centre limit, evidence only).** At one centre with 1s Slater densities the square is
not rational, by analogy with the triangle's dilogarithms. Report the smallest fitting basis.

## Tests and decision rules

1. QG1 as stated. For QR, the K_v kernel against direct 3-D quadrature at five (x, y) pairs,
   relative < 1e-8.
2. For each of Y and S: obtain an annihilating operator (Mellin route or series route), then
   check L f on the function at >= 3 values of t at >= 25 digits, including any source term;
   peel it exactly into factors (rational arithmetic); for each factor of order >= 2 test
   equivalence to Legendre and to F_3's operator (rational gauge, exact verification of any
   map found) and apply the quadric counter of `debug/linking_ceiling_C_quadrics.py`. Test
   whether the function needs every factor (the D11 check).
3. Bounded searches state their bounds. A "none" from a search needs a positive control that
   the same search finds a planted answer.

## Controls (the test is void if one fails)

- gamma -> 0 control for QG: with Yukawa kernels exp(-gamma r)/r on all three edges, the
  gamma -> 0 limit reproduces the Coulomb triangle G(t) of PT2 (numerically at two t).
- Square control for QR: with the fourth kernel removed (an open chain of four), the same
  pipeline must reproduce the tree-of-four rung of Paper 59 `sec:four_density`.
- The chain of three at the same fibre must still come out with its K3 piece (D9-D11).

## Falsifiers

1. QG2: a K3 factor or higher in Y at either gamma, with the controls passing.
2. QR1: the square sits on any rung other than K3, with the controls passing.

## Honest caps, stated before the result

- Fibre statements only; the integrated value is PI item 1, untouched.
- One geminal shape (exponential); the cusp-correct r exp(-gamma r) and the linear-exponential
  geminal of the LiH engine are not covered. Two values of gamma.
- One larger ring. Rings with trees attached are not covered.
- The ring-of-four computation may be too heavy for the precision a rung decision needs; a STOP
  for feasibility is an allowed outcome and must be reported as such, not as a result.
