# FX design pass: exact L-recurrence for the four-centre 1s integral (XY|ZW), row 7 operator column

**Date:** 2026-11-16. **Status:** design + controls only; nothing scored, no recurrence guessed.
**Outcome: STOP at design.** The matched slice is real and rho2_L is exact on it (checked to 39 digits). F_L is
not in the constant class that `na_exact` / `aabc_exact` reduce to: it needs depth-2 iterated exponential
integrals. An AB-style split needs a new reducer, which this pass did not build. So `debug/fx_run.py` was **not**
written. The scored driver would have nothing exact to read.

Files: `debug/fx_engine.py` (exact rho2_L, residue polynomial), `debug/fx_controls.py` ->
`debug/data/fx_controls.json`, log `debug/data/fx_controls.log`.

## 1. Why the matched condition (re-derived; it is right, and it is the only one)

The slice is FC's: X at the origin, Y = D1 z, Z = a n, W = b n, x = n.z. Write s = |r - Z|, t = |r - W| and
mu = cos(r, n). Then

    s^2 = r^2 + a^2 - 2ra mu,   t^2 = r^2 + b^2 - 2rb mu   =>   b s^2 - a t^2 = (b - a)(r^2 - ab).

The right side does not depend on mu. At fixed r the angular integral therefore runs along the conic
S^2 - T^2 = c(r), with S = sqrt(b) s and T = sqrt(a) t. The exponent is

    zZ s + zW t = (zZ/sqrt b) S + (zW/sqrt a) T.

- When this is a light-cone coordinate S + T of the conic, the angular integral becomes a Laplace transform of a
  Laurent polynomial. That happens exactly when **zW/zZ = sqrt(a/b)**.
- The other light-cone direction, S - T, would need zW < 0.
- Off the matched line, the exponent is alpha u + beta c/u with u = S + T, which is a Bessel-K kernel. The angular
  integral then gives K_nu of sqrt(r-dependent) arguments and is not elementary.

So the design pass's condition holds, and nothing weaker works on this slice.

**Rationality.** kappa = sqrt(a/b) must be rational (a/b a ratio of squares), and then everything is rational. The
cases used are M1 (a, b = 1, 4; kappa = 1/2; zZ, zW = 6/5, 3/5; D1 = 6/5) and M2 (a, b = 4/5, 9/5; kappa = 2/3;
zZ, zW = 3/2, 1; D1 = 6/5). Both use zX, zY = 13/10, 1, as FC's G2.

**The b = a reduction sits outside the family.** When a = b the conic degenerates (s = t), so *any* exponents are
exact. That is AB's case. The matched line meets it only on its diagonal zW = zZ. FC's K4 control (AB case A1,
zZ != zW) is therefore a control of the a = b engine branch, not of the matched family.

## 2. The exact representation

With gamma(r) = (b - a)(r^2 - ab)/b and u = s + kappa t, we get gamma/u = s - kappa t,
s = (u + gamma/u)/2 and s ds = (u - gamma^2 u^-3)/4 du. Then

    rho2_L(r) = (2L+1)/(2ra) int_{u-}^{u+} e^{-zZ u} P_L(mu(u)) (u - gamma^2 u^-3)/4 du,
    u- = |r - a| + kappa|r - b|,   u+ = (r + a) + kappa(r + b).

The integrand is e^{-zZ u} times a Laurent polynomial whose powers are all odd, from u^-(2L+3) to u^(2L+1). So

    rho2_L(r) = q+(r) e^{-zZ u+} + q-(r) e^{-zZ u-} + R_L(r) [E1(zZ u-) - E1(zZ u+)],

where R_L(r) = (2L+1)/(2ra) Res_{u=0} e^{-zZ u} P_L(mu(u)) (u - gamma^2 u^-3)/4. The coefficients are rational in r.
u+- are piecewise linear in r, with breaks at a and b.

**Control K1.** The exact split was evaluated at 120 digits and compared with `fc_engine.rho2_point` (arb quadrature
at 400 bits). The comparison covers L in {0, 1, 2, 5, 12, 20, 30}, six radii (one per panel, near D1, and at
r = sqrt(ab), where gamma = 0), and both cases. Worst relative difference: **1.6e-39 (M1), 4.2e-40 (M2)**. At
rational r the split has 4 constants, {e^{-zZ u+-}, E1(zZ u+-)}, and 2 at r = sqrt(ab).

(Note for anyone reusing it: `fc_engine.rho2_point` does not set `ctx.prec` itself. Called cold it runs at 53 bits
and agrees only to about 1e-16, falling to 3e-4 at L = 30. FC's own controls do set it, so FC is unaffected.)

**Control K2: the E1 term does not vanish.** R_L(r) was computed symbolically for L = 0..4 in both cases, and it is
nonzero for every L. At M1, L = 0, it is -81 (r - 2)^2 (r + 2)^2 / (1600 r) (that is, -(1/(2ra)) zZ^2 gamma^2/8).
Every R_L carries the factor gamma^2, which is proportional to (b - a)^2 (r^2 - ab)^2. So **the E1 part of rho2_L
is there exactly because Z != W.** That is why AB (a = b) stayed in the {e, E1, ln, gamma} class.

**Control K3: the potential of rho1 carries a moving E1.** The potential multipole psi_L(r) of rho1 about X is,
up to normalisation, record NA's F_L with the point charge at D2 = r. The exact `na_exact` split, run at D2 = 2 D1
and D2 = D1/2 for L <= 30, has nonzero coefficients on E1((zX + zY) D2) for every L >= 1 (outside D1), and on
E1((zX +- zY) D2) inside. So psi_L(r) contains E1((zX + zY) r) as a function of r.

## 3. Consequence: F_L leaves the AB/NA class

    F_L = (4 pi/(2L+1))^2 int r^2 rho2_L(r) psi_L(r) dr.

Both orders of integration hit a depth-2 object:

- **rho2_L x psi_L:** the integrand contains R_L(r) E1(zZ u(r)) x c_L(r) E1((zX + zY) r) x rational, on
  [max(b, D1), inf) and elsewhere.
- **Potential of rho2, then rho1:** V2_L(r) needs int E1(c r + d) r^(1-L) dr with d != 0. That is a two-variable
  function (int_1^inf e^{-ds} E1(c r s) ds/s), which is then integrated against rho1_L.

The new base integrals have the form

    int_panel e^{lambda r} E1(beta r)/(r - p) dr,   int_panel e^{lambda r} E1(c r + d)/(r - p) dr,
    int_panel e^{lambda r} E1(beta r) E1(c r + d)/(r - p) dr.

The half-line versions include dilogarithms; for example int_0^inf e^{-ct} Ein(t)/t dt = -Li_2(-1/c). None of these is
in {e^{e0}, E1, ln, gamma}.

**Caveat, stated honestly.** I showed that the coefficients of the depth-2 integrands are nonzero. I did not prove
that the depth-2 constants survive in the summed F_L: a cancellation across panels is not excluded. That is
precisely what clause FX2 below would test.

The constant set is still finite and independent of L, so a representation F_L = sum_c q_{L,c} c with q_{L,c}
rational is expected. The rates are L-independent:
- from rho1_L: -zX +- zY;
- from rho2_L: -zZ u+-(r), per panel;
- the E1 functions listed above.

The poles are at r = 0 and at the zeros of the panel lines u+-(r) = 0, which lie outside their panels. The
breakpoints are 0, D1, a, b. The mod-p route (generic random combination of the q-sequences, `na_modguess`) then
applies unchanged.

## 4. What the missing reducer has to do (the build cost)

An integration-by-parts reducer over terms e^{lambda r} r^n (r - p)^-m E, with E in
{1, E1(beta r), E1(c r + d), E1(beta r) E1(c r + d)}, panel by panel:

- partial fractions over the poles {0, p_i};
- positive powers and higher poles reduced by IBP. Each derivative of E lowers its depth and introduces
  e^{-beta r}/r or e^{-(cr+d)}/(cr+d), that is, a new pole;
- bases named abstractly (`('I2', lambda, beta, p, panel)` and the like), evaluated numerically only in the
  controls;
- endpoint values at D1, a, b, which are products of exponentials and E1 values;
- the panel [0, min(D1, a)] needs a finite part at r = 0 with the log singularity. Writing E1 = Ein - ln - gamma
  on that panel keeps the factor entire, so the finite-part endpoint terms are rational (Ein's Taylor coefficients),
  plus E1(d) and e^{-d} from the regular rho2 factor. The ln r pieces cancel inside psi_L for r < D1 (the integrand
  is analytic at 0), but not separately per exponential; the reducer must carry them and check that they cancel.

Estimate: roughly 600-900 lines with flint fmpq/nmod fields as in `na_exact`. Per L, the work is a product of two
Laurent series of length O(L) followed by an IBP ladder of length O(L), so about O(L^2) per term, with hundreds of
terms per L. For L <= 1500 mod p that is an estimated several CPU-hours per prime (not measured). Building the
reducer and its controls is a sprint of its own.

## 5. Proposed test, once the reducer exists (not registrable now)

Cases M1 (dominant contact: Z, rho_Z = 5/6) and a replacement for M2. M2 has D1^2 = ab, so its two contacts
coincide (x_Z* = x_W* = 13/12, rho = 2/3). It served as a closed-form control, but it is unfit for scoring and must
be replaced before registering, for example by a = 4/9, b = 1, D1 = 6/5 (kappa = 2/3): x_Z* = 829/540, x_W* = 61/60, W dominant. Contact points as
FC's `landau()`: x* = (D1^2 + c^2)/(2 D1 c), c in {a, b}. M1: x_Z* = 61/60, x_W* = 109/60.

- **FX0 (gate; failure is VOID).** The engine's exact F_L matches the numeric `fc_engine` F_L at the M cases
  (level A and B) to <= 1e-25 relative, for L <= 30. At a = b with zZ != zW, the engine reproduces AB's
  `aabc_exact` A1 coefficients exactly. The mod-p splits from two primes agree with the fmpq split reduced mod p,
  for L <= 40. NA's (R, Q0) = (8, 31) is reproduced on both primes.
- **FX1 (constant set).** For L = 2..Lmax the set is L-independent and inside the candidate set C built by the
  rule in section 3. Probability 0.7.
- **FX2 (class).** At least one depth-2 base integral has a nonzero coefficient in F_L for every L >= 2.
  Probability 0.85. This is the clause that tests the caveat above.
- **FX3 (operator exists).** The generic sequence has an x-operator mod p within the row limit, at the minimal
  order found by NA's bisection. It must agree on both primes, and be checked at L beyond the guess window by 200.
  Probability 0.5. The order cannot be derived, so it is reported, not scored. A plausible range is 12-30, against
  NA's 8 and AB's operator (to be looked up). The rows needed are (R+1)(Q+2) + 30, so about 1500 coefficients at
  R = 24, Q = 40.
- **FX4 (contacts).** If FX3 holds, x_Z* and x_W* are singular points of the operator, with 7/2 among the local
  exponents (mod-p candidate roots, as in AB's `local()`), and the leading branch exponent 7/2 at the dominant
  contact (FC GO numerically). Probability 0.8 given FX3.
- **FX5 (no other finite real singular point in (1, x_dom*)).** This was FC's stated assumption. Here it is
  testable as the set of non-apparent finite singular points being contained in {+-1, x_Z*, x_W*, the complex
  contact partners}. Probability 0.5. On the matched slice the u-poles could add pseudo-thresholds.

**Decision gate (proposal).**
- **VOID:** FX0 fails, or an instrument fails.
- **GO:** FX1, FX2, FX3 and FX4 hold.
- **STOP:** FX3 fails up to the row limit on both primes, or 7/2 is absent at a contact.
- **BORDERLINE:** otherwise.

## 6. Recommendation

Two options for the PI:

- **(a)** Commission the depth-2 reducer as its own engine sprint (section 4), then register FX as in section 5.
- **(b)** Record row 7's operator column as blocked by a constant-class step. The exact slice exists, but F_L
  needs depth-2 exponential integrals, which no current tool reduces.

I recommend **(b) now, and (a) only if row 7's operator is wanted before the other owed rows**. The FC result
already fills the singular-point column, and the reducer is a larger build than any single prior record's engine.
