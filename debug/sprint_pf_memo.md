# PF1-PF4: the price formula off the symmetric slice (memo)

**Record:** `docs/predictions/2026-10-23_price_formula_slots.md` (committed b84e1231 before computing). Scripts
`debug/pf_*.py`, data `debug/data/pf_*`, assembled in `debug/data/pf_pieces.json`. Nothing outside `debug/` was
touched. Nothing was committed.

**Outcome: STOP for feasibility (the record's decision gate).** Within the run cap, no annihilator was recovered
for the kite phi(z), and so none for Y(t). PF1-PF4 are **UNSCORED**. Every control the record asks to pass before
scoring did pass. The coordinator's discriminator (x = 25/4, 25; z = -25/4, -25) cannot be answered at the level of
local exponents, because there is no operator. It can be answered on the principal sheet, and there it separates the
true object from the direct-pairing object (below).

## Object

Y(t) is the exponential-geminal correlated triangle at the fingerprint fibre:

- vertex 1: density rho_i of scale a = 1;
- vertex 2: density rho_j of scale b = 2;
- vertex 3: apex density rho_k of scale t;
- geminal exp(-gamma_1 r_13) with gamma_1 = 3/2, and geminal exp(-gamma_2 r_23) with gamma_2 = 3;
- Coulomb 1/r_12.

Conventions are those of record QG. Y = (1/(4 pi^2)) t^3 int z^(1/2) (z+t^2)^(-3/2) phi(z) dz, with x = t^2 and
theta = t d/dt.

## Routes

**Kite (route B), re-derived off the slice** (`pf_common.py`). The Feynman parameter u weights the a line and v
weights the gamma_2 line:

- M^2 = u a^2 + (1-u) b^2 + u(1-u) z
- N^2 = (1-v) gamma_1^2 + v gamma_2^2 + v(1-v) z

After the shifts, the Coulomb line carries the momentum (1-u-v) L. So k2 = (1-u-v)^2 z, not the slice's (u-v)^2 z.
The two forms coincide only when sqrt(u(1-u)) and the masses are symmetric under u -> 1-u. The prefactor is
2 gamma_1 gamma_2 a^3 b^3, and the kernel K is the same as in QG.

**Route A (arbiter).** This is the momentum triangle with the closed-form loop I_221(gamma_1, gamma_2, 0) of
`ring_ext_QG_common` and the two density factors written out. It is built directly from the definition.

**Series.** The registered Cauchy route had to be abandoned (see Deviations). It was replaced by direct power-series
arithmetic in z at every Feynman node (flint `arb_series`, `pf_4b_series_direct.py`). The quadrature uses the
substitution mu = M0, sigma = N0, which makes the integrand rational in (mu, sigma) and moves the nearest singularity
away. A Horner form in 1/D cuts the number of series multiplications per node from 17 to 10.

## Controls (all PASS)

| control | result |
|:--|:--|
| loop (3/2, 3, 0) vs direct 3-D quadrature, five triples | worst 2.4e-15 (gate 1e-8) |
| route A vs route B (derived k2 = (1-u-v)^2 z), w = 0.5, 0.9, 2.1 | 0, 0, 1.7e-21 (gate 1e-14) |
| route A vs the slice form k2 = (u-v)^2 z, same w | 1.2e-2, 3.4e-2, 1.2e-1 (the slice form is wrong off the slice) |
| generalised kite at (1,1;2,2) and (1,1;3,3) vs `phi_gem_z`, 3 z each | <= 1.9e-47 at 50 digits |
| same members, Cauchy series vs stored 600-digit CF/MG series, 120 coefficients | 4.9e-77, 2.3e-77 (0.06, 0.03 of the expected bound) |
| two radii (4.5 and 6), 128 coefficients | agree within 0.20 of the expected aliasing + round-off bound |
| direct series vs Cauchy reference, N = 200 | 1.7e-76 (the reference's own node floor); 200x160 = 300x240 |
| Horner kernel vs original kernel | 8.0e-300 |
| N = 300 main (460x460 nodes) vs twin (400x400) | 4.5e-280 (n <= 150), 1.9e-254 (200), 1.6e-232 (250), 1.7e-216 (299) |
| Legendre-map search: stored gamma = 3 map (E = ((x-4)(x-16))^4 (x+8)) recovered from the general (B, gamma) = (1, 3) source | exact, 1/1000 breaks |
| independent map re-check (`pf_map_recheck.py`, rescaled `leg_op` + sympy residual): sources agree with `pf_map`'s for all 4 families and both twists; stored gamma = 3 map | exact, 1/1000 breaks |
| R_V map search (`pf_10_rvmap.py`): R_V into itself; R_V into GS's M14 in GS's box | identity exact; nullity 1, exact identity mod R_V (8 primes) |
| kite-first route (`pf_mellin.py` + `pf_5m_ylog_from_mellin.py`) on the stored gamma = 3 kite | M_Y order 27; Y's annihilator (9, 23) exact by CRT, divides M_Y, 2.2e-362 on 800 coefficients, **equal to the stored L9** |

The radius is the record's disclosed claim: radius 9 for the Cauchy route, to be checked by a second radius. It is
confirmed. The two radii agree, and c_n/c_(n+1) -> -9 (-8.99976 at n = 298). Off the slice, phi(0)/pi =
0.014217790637705248207... has no relation [phi0, pi] up to height 1e20, while the slice control at gamma = 2 recovers
phi0/pi = 173/15552. So the coefficients are not in pi Q, and the exact mod-p route of the slice does not apply.

## The guess (step 1 of the restart): no operator

From the N = 300 series at 280 digits, I searched for operators sum_{i<=d} z^i P_i(theta_z) with deg P_i <= r, at
scale 9 and with at least 40 check rows.

- **Bisection search** (r = 6..18, 20): every order "passed" at tolerance 1e-150, with residuals 1e-151 to 1e-174
  and U = (r+1)(d+1) close to the row limit.
- **Rerun scan** (`pf_5s_scan.py`, 12 workers, once): at orders 6, 9 and 12 the residual falls smoothly with U.
  - r = 6: 2.7e-29 (U = 63) down to 2.1e-191 (U = 259).
  - r = 9: 2.0e-43 (U = 90) down to 5.1e-194 (U = 250).
  - r = 12: 3.5e-62 (U = 117) down to 1.3e-196 (U = 247).
- **Reading:** there is no drop at any order to the data's precision (about 1e-216 or better). The bisection hits
  are the threshold solutions the QG memo describes (its defect 3). There is no genuine kite annihilator with
  (r+1)(d+1) <= 260.

The resolution is bounded by the precision of the series and by U <= N - 40. Under the stop rule, no larger run was
started. **Estimate for a successor run:** N about 720 coefficients at about 600 digits, with about 900x900
substituted nodes, would take about 5-9 hours on 15 cores. That allows U up to about 650. The kite's true size is
unknown; the slice's order 6 doubling to about 12 with x-degree about 50 would need U around 650.

## Discriminator (step 2): x = 25/4, 25 for Y; z = -25/4, -25 for the kite

**Operator level: UNDETERMINED for every point**, because no annihilator is in hand. This covers the two
discriminator points and the eight crossed points {1/4, 1, 9/4, 4, 9, 49/4, 16, 81/4}. No local exponents and no
apparent-versus-genuine verdicts can be given.

**Principal sheet: decided between the two objects.** The direct pairing is exactly the kite with its geminals
swapped. Mapping v -> 1-v turns the swapped kite into the slice-form k2 = (u-v)^2 z, which route A rejects by 1 to 12
percent. Its three-line cut {O1, 12, 23} would put a normal threshold at z = -(1 + 3/2)^2 = -25/4. Measured:

| object | Taylor radius | leading exponent at the radius |
|:--|:--|:--|
| the true kite (route A's) | 9 = (a+b)^2 | beta_n -> -1 (-1.051, -1.025, -1.012, -1.008 at n = 50, 100, 200, 298); possibly with a logarithm, not resolved |
| geminal-swapped kite (direct-pairing object) | 25/4 = (a + gamma_1)^2 | beta_n -> -1/2 (-0.5006, -0.5011, -0.5010, -0.5009 at n = 50, 75, 100, 125) |

So the true kite is **analytic at z = -25/4 on the principal sheet**. The direct-pairing object has a square-root
branch point exactly there. The same holds for the four pseudo-threshold points z = -1/4, -1, -9/4 and -4: they lie
inside the radius, so the true kite is analytic there on the principal sheet. That is consistent with their
predicted second-sheet character, but it does not test it.

Two things are not decided here:

- whether -25/4 (or any point) is a singularity of the kite on another sheet, which is what genuineness in the
  operator requires;
- anything about z = -25 or the cut points -49/4, -16, -81/4, which lie beyond the radius.

For Y, a singular point at x = s needs phi singular at z = -s on some sheet, so the same limits apply.

## Against the registration

| prediction | verdict | evidence |
|:--|:--|:--|
| PF1 (two Legendre factors, crossed slots) | **UNSCORED (STOP)** | no operator. The geometry of the kite is fixed by route A, and its principal sheet has no threshold at the direct pairing's -25/4. This is a property of the integral, not of a composition factor. |
| PF2 (eight genuine points) | **UNSCORED (STOP)** | no operator. Principal sheet: z = -9 singular (exponent about -1), -25/4 and the four pseudo-thresholds analytic, the rest beyond the radius. |
| PF3 (R_V factor) | **UNSCORED (STOP)** | no operator. R_V search tool built and controlled. |
| PF4 (no new type) | **UNSCORED (STOP)** | no operator. |

The decision gate's STOP clause names the kite's operator and its singular points as the deliverable. That operator
was not recovered, so even that clause is not met; what is delivered is the principal-sheet data above.

## Deviations from the record and the brief

1. **Series route.** The Cauchy route was started (3200 points, 700 digits, 1500x1200 nodes). I stopped it: the
   complex kernel costs about 0.9 ms per node pair at 700 digits, which put the run above 50 hours. It was replaced by
   direct power-series arithmetic with the mu/sigma substitution, validated against the Cauchy coefficients and by a
   node twin.
2. **The exact mod-p route does not exist off the slice** (phi_n is not in pi Q), so all coefficients are floats.
3. **The radius check was done at small size** (two radii, 128 coefficients), not on the main series. The main series
   is direct, so it has no radius.
4. **Not run, because nothing was accepted:** the quadrature of Y against an annihilator, the 1/200 controls on Y, the
   1/1000 controls on accepted maps, the planted-split controls, the minimality bounds, and the novelty search
   repetition. The Y grid itself was computed (32 digits, 105 nodes; Y(0.8) = 6.79440718201887e-4).
5. **Tooling built and controlled, ready for a successor run:**
   - `pf_map.py`: the two-scale Legendre sources z_1, z_2 and the direct d_1, d_2;
   - `pf_map_recheck.py`: a re-check for B != 1, which `pm_map_check.py` cannot do;
   - `pf_10_rvmap.py`: maps from R_V;
   - `pf_5m_ylog_from_mellin.py`: Y's exact annihilator from the kite's operator in under a minute, reproducing L9;
   - `pf_discriminator.py` and `pf_kite_local.py`: they report the 25/4 and 25 points explicitly.
   - The probe in `pf_lib.py` is parallelised, and the hyperexponential enumeration runs only for pieces of order
     <= 4.
6. **The first scan died** in a multiprocessing Pool started from a stdin script, which cannot re-import `__main__`
   on Windows. It was rerun once from a file.

## Closing

STOP for feasibility: 300 coefficients at 216-280 digits bound any kite annihilator to (r+1)(d+1) > 260; PF1-PF4 are
unscored. On the principal sheet the integral defined by the record is analytic at the direct pairing's -25/4, where
the geminal-swapped object has a square-root branch point.
