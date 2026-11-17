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

## 7. Reducer built (2026-11-17)

**Status:** the engine and its gate controls are done. Nothing is scored: no recurrence, operator, order or singular
point was guessed, and the constant set's composition was not read off.

**Files.** `debug/fx_reducer.py` is the IBP reducer, 473 lines. `debug/fx_exact.py` gives the exact
G_L = F_L/pi^2 = sum_c q_{L,c} c, 327 lines. `debug/fx_gate_controls.py` writes `debug/data/fx_gate_controls.json`
(log: `fx_gate_controls.log`). The frozen files `na_exact`, `aabc_exact`, `na_modguess` and `fc_engine` are imported,
not edited.

**What was built.**
- **Reducer.** It works panel by panel on e^{lam r} X(r) E(r), with E a product of zero to two factors
  E1(c r + d); the factor (beta, 0) is E1(beta r). The rational part X is kept as partial fractions over the poles
  {0, p_i}. The ladders are the antiderivatives S' + lam S = P and U' + lam U = T - C/t. Each derivative of E lowers
  the depth and adds the pole q = -d/c. Integrals are definite with exact endpoints. At r = 0 the finite part is the
  constant term of the expansion with every ln r dropped, using E1(beta r) = -gamma - ln beta - ln r + Ein and the
  Taylor series of E1(c r + d). Depth-0 integrals reduce to na_exact's {e, E1, ln, gamma} keys. The depth >= 1
  residues are named `('B', lam, E, p, A, B)`. Arithmetic is generic over fmpq or nmod, using flint polynomials.
- **rho2_L, matched case.** Built from the bivariate polynomial Q_l = (2ra)^l P_l(mu(u, r)), one Legendre step per L,
  with the u-ladders evaluated on the panel lines u(r).
- **rho2_L, a = b.** Taken as record AB's dY kernel, that is, NA's density with zX = 0, zH = mu.
- **psi_L.** Assembled from exact antiderivatives of s^{+-} rho1_L. The r-constant parts of psi are integrated once
  per panel, under a tag.

**Design changes.**
1. **Order of integration.** The order is the one in section 3, bullet 1: psi first, then int r^2 rho2_L psi_L. With
   it the a = b limit reproduces AB exactly. The base families are E1(beta r)/(r - p), E1(c r + d)/(r - p) and
   E1 E1/(r - p).
2. **The M2 replacement.** The memo fixes a = 4/9, b = 1, D1 = 6/5, kappa = 2/3 but does not fix the exponents.
   I used M2's zZ, zW = 3/2, 1 (matched: zW = (2/3) zZ). This is case `M2r`. **PI/PM to confirm before
   registering.**
3. **ln constants.** Logs of rationals are expanded over primes, `canon()`, because they are Q-dependent. This is
   applied to both sides of every exact comparison.
4. **Normalisation.** The output is G_L = F_L/pi^2. At a = b, G_L zA^3/2 = q_AB.

**FX0 control table** (every check passes):

| Check | Result |
|:--|:--|
| V1: exact vs fc_engine, L <= 30, <= 1e-25 | M1: max rel 8.7e-29 (A), 7.2e-32 (B); M2r: 3.0e-38 (A), 3.4e-42 (B) |
| V1: evaluation stability, dps 320 vs 400 | 3e-180 (M1), 8e-174 (M2r) |
| V2: a = b (K4A1 = AB A1) vs aabc_exact, L <= 40 | exactly equal at every L |
| V3: mod p (two primes) vs fmpq mod p, L <= 40 | equal at every L, for M1, M2r and K4A1 |
| V4: NA (R, Q0) = (8, 31), both primes | reproduced: order 7 none (Q up to 106), order 8 at Q0 = 31, nullity 1 |

Reading of V1: the agreement is limited by fc_engine's quadrature, since level B is about 1e3 better than level A.
The exact split cancels heavily (|q| reaches about 1e130 at L = 30), so values need about 300 digits.

**Cost** (mod p, case M1, one core):

| Measure | Prime 1 | Prime 2 |
|:--|:--|:--|
| Per L at L = 100 | 0.52 s | 0.50 s |
| Per L at L = 300 | 1.64 s | 1.59 s |
| Cumulative to L = 300 | 241 s | 238 s |

- The local exponent of the per-L cost is 1.05.
- **Extrapolated to L = 1500: about 1.8 CPU-hours per prime.** This is an extrapolation, not a measurement. It
  agrees with the section 4 estimate.
- M2r costs the same as M1 (0.51 s per L at L = 100).
- fmpq up to L = 40 takes 18 s per case.

**Can FX2 test the caveat as designed?** Only partly. FX2 asks whether a named depth-2 base has a nonzero coefficient.
That is a statement about this reducer's basis, and the basis is not proven free:
1. **The basis depends on the order.** In the other order (the potential of rho2 first), only
   E1(c r + d)/(r - p)-type bases arise.
2. **The bases can be related.** The product rule
   d[E1(beta r) E1(c r + d)] = -e^{-beta r} E1(c r + d)/r - c e^{-(cr+d)} E1(beta r)/(cr + d)
   links bases across families. The same holds for the split of a panel's interval and for ln-type relations.

So a nonzero coefficient in this basis is necessary for "depth-2 survives" but not sufficient. Before scoring, FX2
should add one of the following:
- **(a)** reduction modulo the known IBP and product-rule identities to a canonical basis, or
- **(b)** a numerical independence test at high precision: PSLQ over the depth <= 1 span plus the named bases
  appearing, with a positive control.

FX1's candidate set needs the same caveat. Everything FX3 needs is in place: the mod-p sequences and the
NA-validated `generic_sequence` together with `na_modguess` bisection.

## 8. Record-ready design (2026-11-17)

**Status:** instrument complete, controls measured, nothing scored. No operator was guessed on M1/M2r; on M1/M2r
only booleans were recorded (A3: value preserved, rewrite exercised). M2r confirmed matched: zW/zZ = 1/(3/2) = 2/3 =
sqrt(4/9) = sqrt(a/b); x_Z* = 829/540, x_W* = 61/60 (W dominant, rho 5/6). M1: x_Z* = 61/60 dominant, x_W* = 109/60.

### 8.1 FX1/FX2 instrument: route (a), canonical basis (`debug/fx_indep.py`)

**Choice.** Route (a) is exact (fmpq and mod p, every L) and costs ~2 s per 400 L. Route (b) cannot be applied to
F_L itself: its coefficients reach ~1e130, so a relation would need > n x 130 digits. (b) is therefore used only
at the basis level, as a control of (a)'s completeness.

**Identity used (proof in the module docstring).** For a multiset F of factors f = (c, d), x_f = c r + d,
q_f = -d/c: d/dr prod E1(x_f) = -sum_f e^{-d_f} e^{-c_f r}/(r - q_f) prod_{g != f} E1(x_g). Integrating over a
panel gives sum_f m_f e^{-d_f} B(-c_f, F - f, q_f) = -[prod E1]_A^B. This is valid when no q_f lies inside the
panel, every c_f > 0 at an infinite end, and A = 0 is a Hadamard finite part. Every B(lam, E, p) with lam != 0 lies in exactly
one set, F = E + {(-lam, lam p)}. A set is applied only when all its factors are natural for the panel (put there
by the engine). The member with the largest repr is eliminated.

**Controls (measured).**

| ID | Check | Result |
|:--|:--|:--|
| A1 | every valid natural pair set plus 6 triples, all panels (incl. FP at 0 and to infinity), C1 and C2, 60 digits | 38 and 42 sets, max rel 3.6e-60, pass |
| A2 | planted: whole identity sets plus depth-0 junk | depth-2 26 -> 0 (pos); one canonical member survives (neg); value and mod-p equal |
| A3 | real F_L, L = 2..5, raw vs canonical at 350 digits | C1/C2: 64 eliminations each, rel <= 3e-320. M1, M2r: value preserved, rewrite exercised (booleans) |
| A4 | canonical mod p == canonical fmpq mod p, L <= 10, two primes | pass (C1, C2) |
| A5 | a = b (K4A1): no named base; canonical() is the identity | pass |
| B1-B3 | LLL, C1 panel [0, 1], 300 digits, H = 1e6 | identity recovered exactly (pos); planted relation found; canonical set: none (also at 400 digits) |
| B4 | LLL on **all 77** canonical constants of C1 (L <= 6; 16 depth-2), 450 digits, H = 1e4 | none found (also at 550 digits); adding one eliminated member: relation found, involving it |

C1 = (a, b, D1) = (1, 9/4, 7/5), zZ, zW = 3/2, 1 (a < D1 < b, as M1). C2 = (1/4, 1, 7/5), zZ, zW = 2, 1
(D1 > b, as M2r). Both are matched, unscored, and use zX, zY = 13/10, 1.

**Alphabet (FX1's candidate set, `fx_run.alphabet`).** The section 3 rule is mechanised from the engine's input
structure (rates, factors, poles, breakpoints), closed under the reducer's derivative step, and the identity
partners are added. Calibrated on C1/C2 (`ctrl_alphabet`, L <= 60, two primes). The first pass missed 3 and 2 E1
atoms: the depth-0 base's E1(lam p) at the endpoint A = 0 was omitted, a rule bug. After that fix: 0 outside on
both. **Disclosed:** on C1/C2 the canonical set is L-independent for 2 <= L <= 60 and depth-2 keys are nonzero at
every such L. These are unscored geometries, but they bear on FX1 and FX2. The probabilities below take that into
account.

### 8.2 Clauses (M1, M2r; each clause on both cases)

- **FX0 (gate, VOID on failure):** section 7 V1-V4 as run, plus A1-A5, B1-B4, `ctrl_alphabet`, and the
  positive controls in 8.3, rerun with the frozen files. Also VOID if the natural factor sets are L-dependent
  (L = 2..6) or if any annihilator fails its held-out rows.
- **FX1 (constant set, P = 0.8):** in the canonical basis, the key set (present = nonzero mod either prime) is the
  same for every 2 <= L <= 1500, and every atom lies in `alphabet(case)`.
- **FX2 (class, P = 0.9):** for every 2 <= L <= 1500 some canonical depth-2 key is nonzero (either prime).
  Caveat as worded: this is modulo the product-rule family. B4 found no further relation at height <= 1e4 on C1.
- **FX3 (operator, P = 0.5):** the window is L <= 1300; 200 more rows (to L = 1500) are held out. For seeds 2026
  and 7 on both primes, NA's bisection finds an operator at order R0 <= 32, and every annihilator at R0..R0+2
  vanishes on all rows to L = 1500. The minimal operator (their gcrd, MO's construction) has the same
  (order, degree) in all 4 runs. The search shape and the minimal shape are reported, not scored.
- **FX4 (contacts, P = 0.8 given FX3):** in all 4 runs, x_Z* and x_W* are factors of the minimal operator's
  leading coefficient, and 7/2 is an indicial root at both (`na_modanalyze.local`). The memo's "leading branch
  exponent 7/2 at the dominant contact" cannot be decided mod p (AB's 3/2 is a root whose branch cancels). It is
  dropped from scoring and left to FC's numerical GO.
- **FX5 (FC's assumption, reported, outside the gate; tri-state):** genuine (class M/X/L/P) weight-fixed points of
  the minimal operator other than the predicted set {+-1, x_Z*, x_W*} are located by CRT over the two primes,
  height <= 1e6. The result is **False** if one lies in (1, x_dom*); **None** if some genuine factor is non-linear
  or not reconstructable; **True** otherwise. P(True) 0.4, P(None) 0.4, P(False) 0.2. Strict containment in the
  predicted set is reported only, and is expected to fail: NA's own operator has a genuine log point at
  -22633/21175 (reproduced below), and AB's has three extra log points. The "complex contact partners" of section
  5 were never derived, so they are not used.

**Gate.** VOID: FX0 fails. GO: FX1-FX4 hold on both cases. STOP: no operator on either prime up to R = 32, or a
contact present without 7/2. BORDERLINE: otherwise. Overall: GO only if both cases are GO; STOP if either is STOP.

**Caps.** LMAX = 1500, HOLD = 200, RMAX = 32 (at R = 32 the window allows degree <= 34; at R = 24, <= 45). Margin
30, two extra orders for the gcrd. Primes are AK's, 1073741441 and 1073741381 (< 2^31, as intclass_lib needs).

### 8.3 Scoring driver (`debug/fx_run.py`) and its positive controls (measured, AK primes, both)

| Control | Result |
|:--|:--|
| NA at S1's N3 geometry (`na_exact`, L <= 1000; window 800 + 200 held out) | search (8, 31), minimal (8, 25), weight-fixed points identical to frozen `s1_controls.json`, 200 rows verified, 3/2 a root at x*; the strict-FX5 analogue fires on the log point -22633/21175, which CRT recovers exactly. About 6 min per prime, 18 s of it the operator. |
| AB A1 frozen coefficients (`s1_coeffs_A1`) | search (17, 48), minimal (14, 111), weight-fixed points identical to frozen W0, 7/2 at x* = 125/117. About 62 s per prime. |
| K4A1 (a = b) through FX engine + canonical, L <= 1200 (window 1000 + 200) | search (17, 48), minimal (14, 111) = AB's; points identical to W0 per prime; verified 200. About 17 min per prime with 7 jobs running. |
| End-to-end `score` + `summary` on K4A1 (scratch dir) | FX1 holds, **FX2 False at every L** (correct: a = b is AB's class), FX3 holds (14, 111), FX4 holds, FX5 None (two AB log points not reconstructable at height <= 1e6; one at -179401/179400). Gate BORDERLINE, as it should be for this negative case. |

Fixed after these runs: a tuple-vs-list comparison, FX1's absent-key default, and FX5's single-prime ratrec
(bound sqrt(p/2) gives spurious positions; now CRT, height <= 1e6).

### 8.4 Runtime

- Coefficients: M1 runs 0.52 s per L at L = 100 and 1.64 s per L at L = 300, extrapolating to about 1.8 CPU-h to
  L = 1500 per prime (section 7). C1 matches: 224 s to L = 300.
- The canonical rewrite adds about 2 s per 400 L.
- The operator stage takes 18-44 s on the controls. For R <= 32 with 1300 rows, budget at most 20 min.
- **Per (case, prime) about 2.2 h; 4 jobs in parallel give about 2.5 h wall.** That is under the 6 h cap.

### 8.5 Commands (after registration)

```
python debug/fx_indep.py controls                       # FX0 part A/B rerun (~6 min)
python debug/fx_run.py ctrl_alphabet                    # alphabet calibration rerun (~30 s)
python debug/fx_run.py posctrl na 1073741441  (& 1073741381; ab; k4a1)   # FX0 positive controls
python debug/fx_run.py score M1 1073741441 --registered   (4 jobs: M1, M2r x 2 primes, in parallel)
python debug/fx_run.py summary                          # -> debug/data/fx_run.json, clauses + gate
```

Data: `debug/data/fx_indep_controls.json`,
`fx_runctl_*.json`, scored `fx_run_<case>_<p>.json` and `fx_seq_<case>_<p>.json`.

**For the PM.** Confirm: M2r's exponents; FX5 tri-state outside the gate; FX4 narrowed to "7/2 a root"; FX1/FX2
probabilities, which were raised after the C1/C2 controls. Side finding: the single-prime rational positions of
AB/MO's extra log points (20665/15246, 21025/14886) are ratrec artefacts. MO's verdict uses only the classes, so it
is unaffected.
