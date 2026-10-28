# Pre-registered prediction: Gaussian orbitals under the contact rule, the bond-angle series of symmetric water (WG)

**Registered:** 2026-10-28, before any F_l in this file is computed. PI-directed (a "Gaussian orbitals under the
contact rule" test). Drafted by a sub-agent; reviewed and committed by the PM.
**Commit status:** to be committed on its own, before any computation it scores.
**Rests on:** Paper 59's water passage (the contact form (pi^2/630) g_Y g_Z |Y - Z|^7 of the leading bond-angle term),
`docs/predictions/2026-10-27_water_cusp_next_term.md` (WA6, GO), `docs/predictions/2026-10-27_water_zero_contact.md`
(WA7, STOP on a misstated sign clause), `docs/predictions/2026-10-27_water_zero_contact_fresh.md` (WA8, GO).

## The question

Paper 59 says the slow large-l tail of the bond-angle Legendre series of (XY|XZ), x = cos(gamma), comes only from the
two hydrogen cusps, through their contacts g_Y, g_Z. A Gaussian function has no cusp anywhere. So the contact rule
says an all-Gaussian (XY|XZ) has no power-law tail at all. Two parts:

- (i) One s Gaussian on each centre. The integral is a Boys function of |P - Q|^2, which is linear in x. So it is an
  entire function of x and its Legendre coefficients fall faster than any geometric rate. The shape is derived below.
- (ii) STO-3G and STO-6G contractions of the 1s Slater functions used in WA5-WA6. A finite Gaussian sum, so its tail
  must also fall faster than any power. At moderate l it should follow the Slater l^-8 coefficients. Where it leaves
  them (a crossover l*) and how l* moves with n are predicted from a model, stated below.

Geometry: X at the origin, |XY| = |XZ| = D = 1.809 bohr, as in WA6.

## The derivation (done before registering)

Scripts `debug/wg_derive.py`, `debug/wg_derive2.py`, `debug/wg_derive3b.py`; logs `debug/data/wg_derive.log`,
`debug/data/wg_derive2.log`, `debug/data/wg_derive3b.log` (json beside them). None computes a scored F_l.

**Closed form.** For s Gaussians the product on X and Y is one Gaussian, exponent p = a_X + a_Y, centred at
P = (a_Y/p) Y, of charge M1 = c_X c_Y (pi/p)^(3/2) exp(-a_X a_Y D^2/p). Likewise Q, q, M2 for X and Z. Then
(rho1|rho2) = M1 M2 2 sqrt(mu/pi) F0(mu |P - Q|^2), mu = pq/(p + q), F0(T) = int_0^1 e^(-T u^2) du. With s1 = a_Y/p,
s2 = a_Z/q, mu |P - Q|^2 = T0 - T1 x, T0 = mu D^2 (s1^2 + s2^2), T1 = 2 mu D^2 s1 s2. Using
int_{-1}^{1} e^(b x) P_l(x) dx = 2 i_l(b) (i_l the modified spherical Bessel function), each quartet gives

  F_l = M1 M2 2 sqrt(mu/pi) (2l+1) int_0^1 e^(-T0 u^2) i_l(T1 u^2) du
      = M1 M2 2 sqrt(mu/pi) (2l+1) sum_k T1^(l+2k) / (2^k k! (2l+2k+1)!!) e^(-T0) M(1, l+2k+3/2, T0) / (2l+4k+1),

with M Kummer's function (from int_0^1 u^m e^(-T u^2) du = e^(-T) M(1, (m+3)/2, T)/(m+1)). For a contraction, sum over
all primitive quartets (a on X in rho1, b on X in rho2, j on Y, k on Z).

**Checks of the closed form, on toy values only (not a scored case).** At mu = 0.35, |P| = 1.1, |Q| = 0.8 the series,
the u-integral, a direct 40-digit Legendre projection of F0(T0 - T1 x), and an independent k-space route
(2/pi) int dk e^(-k^2/(4 mu)) j_l(k|P|) j_l(k|Q|) agree to 15 digits at l = 0, 3, 8. **Found while checking:** the
u-integral form, done by plain mpmath quadrature, loses accuracy at large l (u^(2l) piles up at u = 1; it gave a wrong
q_l by l = 40 on the toy). Phase 2 uses the series form, and the k-space route as the independent control.

**Signs.** Every term above is positive when all contraction coefficients are positive (i_l(z) > 0 for z > 0,
M > 0, T1 > 0 for every quartet since s1, s2 > 0). STO-3G and STO-6G 1s coefficients are all positive. So F_l > 0 at
every l, for every case here. No (1 - x)^a power appears, so the WA7 sign issue (the sign of Gamma(a + 1 - l)) does
not arise. The Slater reference S_l is also positive at l = 15..60 (guard `tests/test_paper59_water_bond_angle.py`).

**(i) Large-l shape.** Write F_l = Pref e^(-T0) T1^l / (2l+1)!! R_l. Expanding in 1/l (sympy):
R_l = 1 + (T0 + T1^2/4)/l + O(l^-2). The cleaner object is the step ratio q_l = F_(l+1) (2l+3) / (F_l T1) = R_(l+1)/R_l.
On the symmetric line T0 = T1 = T (both cases in (i)):

  q_l = 1 - T(T+4)/(4 l^2) + T(T+4)/l^3 + T(2T^3 + 12T^2 - 17T - 196)/(16 l^4) - T(9T^3 + 48T^2 + 76T - 272)/(8 l^5) + O(l^-6).

Checked on toy T = 1.2 against the series form: the residual after the l^-5 term is 6.5e-7, 1.2e-8, 2.0e-10 at
l = 20, 40, 80, falling like l^-6. (The general (T0, T1) form to l^-4 was checked the same way at (1.3, 1.1).)
So F_(l+1)/F_l = T1/(2l+3) to leading order: faster than any geometric rate, and no power law.

Cases for (i), exponents chosen, not fitted:

| Case | a_X | a_H | T0 = T1 = T | q2 | q3 | q4 | q5 |
|:--|:--|:--|:--|:--|:--|:--|:--|
| g1 | 1 | 1 | 1.6362405 | -2.305561 | 9.222245 | -18.70707 | -4.149909 |
| g2 | 0.3 | 2 | 5.6912713 | -13.78891 | 55.15565 | 165.2682 | -2400.561 |

**(ii) The contractions.** STO-3G and STO-6G for H, from the Psi4 basis library (files `psi4/share/psi4/basis/
sto-3g.gbs` and `sto-6g.gbs`, master branch, downloaded 2026-10-28; sha256 b472adfb... and 698150a1...). The STO-3G
file is the EMSL Basis Set Exchange export and cites Hehre, Stewart and Pople, J. Chem. Phys. 51, 2657 (1969) for H-Ne.
The STO-6G file is a Basis Set Exchange export marked "Data from Gaussian09"; the primary source is Hehre, Ditchfield,
Stewart and Pople, J. Chem. Phys. 52, 2769 (1970). **I have not checked the numbers against the 1969/1970 papers
themselves**; only against the Psi4 copies. Psi4's H rows are for zeta = 1.24; dividing exponents by 1.24^2 gives the
zeta = 1 expansions:
- STO-3G: alpha = 2.227661, 0.4057712, 0.1098175; d = 0.15432897, 0.53532814, 0.44463454.
- STO-6G: alpha = 23.10303, 4.235916, 1.185057, 0.4070989, 0.1580884, 0.06510954; d = 0.009163596281,
  0.04936149294, 0.1685383049, 0.3705627997, 0.4164915298, 0.1303340841.
- For exponent zeta: alpha_i -> zeta^2 alpha_i. Primitives normalized; the contraction renormalized exactly (raw norms
  0.9999999909 and 0.9999999999). Transcription check: overlap with the zeta = 1 Slater 1s is 0.9998347 (STO-3G) and
  0.99999938 (STO-6G).
- Sets as in WA6: E1 (zeta_X = zeta_H = 1) and E2 (zeta_X = 2.25, zeta_H = 1.24). All three centres contracted with
  the same n. Four cases: E1-3, E1-6, E2-3, E2-6.

**(ii) The crossover model.** Near Y the density chi_X chi_Y is about chi_X(Y) chi_Y(r), so the cusp-cusp part of
(XY|XZ) is about chi_X(D)^2 V(|Y - Z|), with V(w) the Coulomb integral of the H orbital with itself, each copy used as
a charge distribution, at separation w = D sqrt(2(1 - x)). The model ratio is

  rho_n(l) = [chi_X^(n)(D)/chi_X^S(D)]^2 Leg_l[V^(n)] / Leg_l[V^S].

V^S is the textbook 1s-1s Coulomb closed form; V^(n) is a sum of erf terms, so Leg_l[V^(n)] uses the series above.
What the model gets right and wrong, measured against the stored Slater S_l (WA5 data, known): its Slater version
chi_X^S(D)^2 Leg_l[V^S] over S_l is 1.42, 1.18, 1.06, 1.025, 1.016, 1.007 at l = 8, 12, 20, 32, 40, 60 (E1), and 1.44,
1.20, 1.07, 1.03, 1.02, 1.008 (E2). So it has the leading term exactly and misses the next order (the X-gradient
terms), which is about 20% at l = 12 and 2-3% at l = 32-40. In the ratio rho_n part of that error may cancel; how much
is not known. The model also ignores that the X factor shifts each Gaussian centre toward X (P = s1 Y, s1 < 1).

Model results (first l >= 6 with rho_n < 1/2, called l*):

| Case | alpha_max on H | D sqrt(alpha_max) | model l* | rho at l*-1, l* | max rho before l* (at l) | l*/(D sqrt(alpha_max)) | model local power at l*+10 |
|:--|:--|:--|:--|:--|:--|:--|:--|
| E1-3 | 2.22766 | 2.700 | 11 | 0.576, 0.341 | 1.133 (7) | 4.07 | 44.0 |
| E1-6 | 23.1030 | 8.695 | 32 | 0.573, 0.475 | 1.219 (23) | 3.68 | 25.5 |
| E2-3 | 3.42525 | 3.348 | 13 | 0.642, 0.417 | 1.221 (9) | 3.88 | 39.7 |
| E2-6 | 35.5232 | 10.782 | 40 | 0.528, 0.450 | 1.235 (28) | 3.71 | 23.7 |

Reading: a Gaussian fit copies the cusp only down to a length of about 1/sqrt(alpha_max), and the Legendre index l
probes |Y - Z| of about D/l. So l* is about 3.7-4.1 times D sqrt(alpha_max), and l*(6)/l*(3) should be near
sqrt(alpha_max(6)/alpha_max(3)) = 3.22 (model: 2.91 for E1, 3.08 for E2). Past l*, the model's local power
log(F_l/F_(l+4))/log((l+4)/l) is far above the Slater value 8 and keeps rising. Before l*, the STO-6G model ratio dips
to about 0.94 near l = 15-19 and then overshoots to about 1.22 near l = 23-28: the Gaussian fit's momentum tail
oscillates around the Slater k^-4 tail before it dies.

**What is not derived.** The model's accuracy at l* (so the tolerance on l* below is a judgment, not a bound). The
size of the overshoot. The low-l plateau (the model is poor below l of about 12). The asymptotic constant of the
contracted series (it is set by the quartet with the largest T1, and the second-largest T1 is within about 1-6% of
it, so that regime is far out). These are scored only as secondary items or left as reported-only.

## What I know before comparing (disclosed)

- No F_l for g1, g2 or any contracted case has been computed, at any l, by any route. The closed form was evaluated
  only on toy parameters (listed above).
- The stored Slater S_l (`debug/data/wa5_cusp.json`, "E1 (1.809, 1.809)" and "E2 (1.809, 1.809)", grid level 2,
  l = 0..60, relative noise floor below 4.4e-5 at l = 30..60) are known, and were used to calibrate the model (numbers
  above). Their local power has not been computed.
- (i) is close to certain once the closed form is right. It tests the code and a textbook identity, not new physics.
  It is kept because it is the cleanest form of "no cusp, no power".
- The model ratios for every l = 4..60 are in `debug/data/wg_derive2.log`.

## Model (what phase 2 computes)

- F_l^G (part i) and F_l^(n) (part ii): the series form above, mpmath at 50 digits, l = 0..60, summed over all
  primitive quartets (81 for STO-3G, 1296 for STO-6G).
- S_l: the stored WA5 values, and recomputed for control (a).
- ratio_l = F_l^(n) / S_l. Crossover l*_obs: the first l >= 6 with ratio_l < 1/2.
- Local power p_l = log(F_l/F_(l+4)) / log((l+4)/l).

## Predictions

**WG-0 (controls; run FIRST).**
- (a) The test file's `F_coeffs` (level 2) reproduces the stored S_l for E1 and E2 at l = 10, 20, 30, 40, 50, 60 to a
  relative 1e-10. If (a) fails, part (ii) is void.
- (b) Independent route: the k-space integral (2/pi) int dk e^(-k^2/(4 mu)) j_l(k|P|) j_l(k|Q|), per quartet, agrees
  with the series form to a relative 1e-10 for g1 and g2 at l = 0, 5, 10, 20, 30; for E1-3 at l = 0, 10, 20; for E1-6
  at l = 0, 20. A case failing (b) is void (for the contracted sets, a failure voids that set's two cases).
- (c) Resummation: sum_l F_l P_l(cos 104.5 deg), l = 0..60, equals the direct Boys evaluation at the bent geometry to a
  relative 1e-14, for g1 and E1-3.
- (d) The toy check of the q_l expansion (T = 1.2, l = 40: residual after l^-5 below 1e-7) re-run in the scorer.
- (e) Negative control: on the stored Slater S_l, the local power p_l at l = 21, 23, 42, 50 lies in 7..9. This shows
  the far-tail test WGb3 can fail on a cusped function.

**WGa (part i; no cusp, no power).** For g1 and g2:
- F_l > 0 for every l = 0..60;
- for every l = 40..60, |q_l - (1 + q2/l^2 + q3/l^3 + q4/l^4 + q5/l^5)| < 0.5 max(|q5|, 10) / l^5;
- p_30 > 50.
- Falsifier: any clause broken in a non-void case.
- Probability 0.95.
- Novelty: KNOWN. The Boys-function closed form is textbook; e^(z cos theta) = sum (2l+1) i_l(z) P_l(cos theta) is
  Abramowitz and Stegun 10.2.36 (location to be re-checked at scoring). A hit is a check of the code.

**WGb1 (part ii; sign).** F_l^(n) > 0 for every l = 0..60 in every non-void case.
- Falsifier: one negative F_l.
- Probability 0.97.
- Novelty: KNOWN (follows from positive coefficients). A check.

**WGb2 (part ii; the crossover).** In every non-void case l*_obs exists, it lies within 3 of the model l*
(E1-3: 8..14; E1-6: 29..35; E2-3: 10..16; E2-6: 37..43), and in each set with both cases live,
2.4 <= l*_obs(6) / l*_obs(3) <= 3.8.
- Falsifier: any l*_obs outside its window, no crossover by l = 60, or a ratio outside 2.4..3.8.
- Probability 0.55.
- Novelty: PROBABLY KNOWN as a mechanism (a Gaussian expansion copies a cusp only down to about 1/sqrt(alpha_max), so
  angular and momentum tails are cut there); this is general basis-set lore. Not searched for this integral or for
  the bond-angle series. Not searched means it cannot count as a finding.

**WGb3 (part ii; no power past the crossover).** In every non-void case, at l0 = model l* + 10 (E1-3: 21; E1-6: 42;
E2-3: 23; E2-6: 50): p_l0 > 16 (twice the Slater value) and p_(l0+5) > p_l0.
- Falsifier: either inequality fails in any non-void case.
- Probability 0.85.
- Novelty: KNOWN in principle (finite Gaussian sums are entire). A check of the contact reading.

**WGc (secondary; scored and reported, not in the gate).**
- (c1) STO-6G plateau: |ratio_l - 1| < 0.03 for every l = 0..12, in E1-6 and E2-6. Probability 0.7.
- (c2) STO-6G overshoot: max ratio_l over l in [model l* - 14, model l* - 4] (E1-6: 18..28; E2-6: 26..36) exceeds 1.08,
  in both. Probability 0.5.
- Novelty: NOT SEARCHED. Reported only, whatever the outcome.

## Decision gate

- **GO:** WG-0 (a), (c), (d), (e) pass; g1 and g2 both live under (b), and at least three of the four contracted cases
  live; WGa, WGb1, WGb2 and WGb3 hold. Then no cusp means no power-law tail, and a Gaussian contraction follows the
  Slater cusp series only up to a crossover set by its tightest exponent, at the predicted l.
- **BORDERLINE:** as GO, but WGb2 fails. The contact rule's "no cusp, no tail" holds; the crossover model does not
  place l*.
- **STOP:** WG-0 (a), (c), (d) or (e) fails; fewer than three contracted cases or either of g1, g2 live; or WGa, WGb1
  or WGb3 fails (a power-law tail, or a negative coefficient, in a cusp-free basis).

## Caps

- At most 45 min wall, at most 4 processes. No change to `geovac/`, `tests/` or `papers/`.
- Symmetric geometry, s functions only, 1s on X, the two published contractions. Not a claim about p or d Gaussians,
  even-tempered sets, or unequal bond lengths.

## Disposition

(Appended after comparing. The predictions above are not edited.)
