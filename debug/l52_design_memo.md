# L52 design memo: the contact-exponent count and the missing 5/2 in (AA|BC)

Date 2026-11-17. Design only; nothing committed, no paper or git touched, no SEALED file opened.
Draft record: `docs/predictions/2026-11-17_aa_bc_contact_amplitudes.md`. Code: `debug/l52_local.py` (formulas), `debug/l52_controls.py`,
`debug/l52_predict.py`, `debug/l52_run.py` (Phase 2 + gate fire test). Data: `debug/data/l52_*`.

## 1. The count, written out

Labelling note: the lit memo says x* is where |B - C|^2 vanishes. In AB's frame (B origin, C on the z axis at D1, A at
D2 and angle g) the quantity is |A - C|^2 = 2 D1 D2 (x* - x); the memo used NA's X, Y, Z letters. Same point.

Reduce (AA|BC) to one electron: f(x) = int V(r - A) rho_BC(r) d^3r, V the potential of rho_A. f is a convolution;
as a function of A it is non-analytic only where a singular point of V (at A) meets a singular point of rho_BC (at
C or B). Near those points:

- V(s) = E(s^2) + sum_{p odd} v_p s^p (E entire);
- rho_BC(r) = smooth + sum_{q odd} g_q(r) |r - C|^q, g_q = N (-zC)^q / q! e^{-zB |r - B|}, analytic near the segment.

The pinch of one term of each: I_pq(a) = int |s|^p |s - a|^q g(A + s) d^3s, a = C - A. Combine the two powers with
one Feynman parameter u, shift s = t + u a, integrate t over R^3:

    I_pq = analytic + C_pq |a|^{p+q+3} [ int_0^1 u^{(p+1)/2} (1-u)^{(q+1)/2} g(A + u a) du + O(|a|^2) ],
    C_pq = pi^{3/2} Gamma(-(p+q+3)/2) / (Gamma(-p/2) Gamma(-q/2)).

Each step is explicit: the Feynman formula (by analytic continuation in the powers), the t-integral
int (t^2 + M^2)^{-L} d^3t = pi^{3/2} Gamma(L - 3/2)/Gamma(L) M^{3-2L} with M^2 = u(1-u) a^2, and a Taylor
expansion of g in t whose even moments bring extra factors M^2, i.e. higher powers of |a|^2. The u-line is the
memo's Morse-Bott critical segment; no stationary-phase argument is needed. With |a|^2 = 2 D1 D2 (x* - x) the
pair (p, q) gives a branch (x* - x)^{(p+q+3)/2} whose leading amplitude is the segment average at a = a* (complex,
null); there the segment's distance to B is real, |A + u a* - B|^2 = (1-u) D2^2 + u D1^2. This is Paper 59's
Landau-pinch average for (XY|XZ), one dimension lower, and the same constant the corpus already uses:
C_pq Beta((p+3)/2, (q+3)/2) equals `aaeq2_local.K(q, p)` to 20 digits for all five pairs tested (control K2).

**Tier.** A derivation of the local germ, conditional on three things, which is where it could fail:
(a) g is analytic along the complex segment (true here: R_B(u)^2 > 0 on [0, 1]; it would fail if B's null cone met
the segment, i.e. at a second Landau point); (b) the O(|a|^2) remainder is a convergent series (not proved; the
corpus's local expansions assume the same); (c) the Legendre coefficients are governed by x* alone (Darboux/Lehmann;
f analytic inside the x* ellipse, AB P3b). It speaks about the function f only, never about which exponents an
annihilating operator has.

## 2. Why the lit memo got 5/2, and why AB P3a got 5/2

**The lit memo.** It entered (AA|BC) as "cusp x cusp, sum alpha = 1": alpha = 1/2 for A's density cusp |r1 - A| and
1/2 for C's. But A's cusp sits in r1 and reaches r2 only through 1/r12; integrating r1 out (Poisson) turns
|r1 - A| into |r2 - A|^3. In the one-electron form the A-side power is p = 3, not 1. Counting all three factors in
6-D gives the same: alpha = 1/2 - 1/2 + 1/2, plus 3/2 per 3-D integration, so gamma = 7/2. The memo's row "u^3 v,
Sigma alpha = 2 -> 7/2" is the true leading row; its "5/2" row does not exist.

**The derivation of the vanishing (function level).** The 5/2 coefficient of f is v_1 g_1 C_11 (seg) +
v_{-1} g_3 C_{-1,3} (seg). From (-lap + kappa^2) V = 4 pi rho_A with rho_A bounded, v_{-1} = 0 (no point charge)
and v_1 = 0 (lap(s) = 2/s would need rho ~ 1/s). Both vanish identically: every zA, every geometry, any bounded
spherical density on A, Coulomb or Yukawa line (sympy, control K1: v_-1 = v_1 = 0, v3 = 2 zA^4/3, v5 = 8 zA^6/45).
The same kills 3/2. So the 5/2 coefficient is zero at A1, A2 and everywhere, not by an accident of the geometry.

**This is already in the corpus.** AB's P3b (registered before computing, HELD) says "f has no (x*-x)^{3/2} or
^{5/2} branch: its leading branch is (x*-x)^{7/2}", by the Fourier form of the same argument (V_A^ ~ 64 pi zA^4/k^6,
the k^-2 point-charge term cancels); `aaeq2_local.py` builds the local expansion from v_n s^n, n odd, with v_1 = 0.
At function level the count adds nothing: it is the corpus mechanism, and the memo's 5/2 was an input error.

**AB P3a (the coordinator's pointer) is a different claim.** P3a predicted 5/2 as an *exponent of the minimal
operator*, by inheritance: the operator annihilates every constant's component, NA's components were taken to enter f
unchanged, and NA's operator has 3/2 and 5/2. AB's own record defect 1 shows the premise false (the four exponential
components also receive Y and W contributions; only NA's two E1 components are unchanged). Nothing in the count or
in Poisson kills an operator exponent: the operator sees every component separately, and an exponent can be present
with zero coefficient in f, which is exactly what happens to 3/2 (an exponent in MO's operator, absent from f). So:

- 5/2 absent from f: derived (above), and already scored (AB P3b).
- 5/2 absent from the operator while 3/2 is present: **not explained** by either argument. A per-component version
  of the count (each constant's region piece has its own singular terms, with endpoint contributions at D1, D2)
  might explain it; that is a different object and is not attempted here. Flag for the PI; it is the only open
  piece of the 5/2 question.

## 3. The 5/2 coefficient computed

Analytic (above): identically zero. Measured, indirectly and already on file: AB's beta_eff(220) = -4.106 at A1
and the 7/2 toy gives -4.119 at the same x*; a surviving 5/2 term would pull beta toward -3. No new computation
needed.

## 4. Fresh, nonzero predictions

The count yields no new exponent anywhere I looked: a different bounded density on A keeps p >= 3 (a cusp-free one
gives 9/2, already seen in AP's monopole piece and AS/AY's cancelling totals); a Yukawa line keeps v_-1 = v_1 = 0;
the four-centre contact gives 7/2 (FC) with 9/2 next and no integer exponent (no even power without a log, and no
log source); 1/r12^2 would give an even power with a log, (x*-x)^3 ln(x*-x), but needs a new kernel engine
(Legendre expansion of 1/r12^2) and is not under 3 h. What the count does give, at no new code cost, is
**amplitudes**:

- **P1. The 7/2 amplitude at unequal arms**, Phi7 = v3 N (-zC) (-2 pi/35) (2 D1 D2)^{7/2}
  int_0^1 u^2 (1-u) e^{-zB R_B(u)} du. Never computed for (AA|BC) away from D1 = D2 (AE2b had x* = 1, a real point,
  where the segment collapses). Values: A1 0.4840671140, A2 0.3115657420, A3 (fresh) 2.4426119168.
- **P2. A Yukawa line leaves the 7/2 amplitude unchanged.** v3 = -4 pi rho_1/12 does not involve kappa, so f_kappa
  and f share every v3-built term; the difference leads with 9/2 from Delta v5 = kappa^2 zA^4/45, amplitude
  DPhi9 = Delta v5 N (-zC) (-2 pi/63) (2 D1 D2)^{9/2} int u^3 (1-u) e^{-zB R_B} du (A1, kappa = 7/10: 0.0263619;
  A3, kappa = 1: 0.3493270). The Yukawa integral needs only the existing engine:
  f_kappa = mu^4 a (Y_kappa - Y_mu) - mu^3/(2(mu^2 - kappa^2)) W_mu (partial fractions of V_kappa^; kappa -> 0
  checked).

Neither is a new mechanism: P1 is Paper 59's segment average carried to (AA|BC); P2's kappa-independence is the new
derived statement. Cost: AB's engine ran 267 s per case at L = 240, 2000 digits; Phase 2 is three cases at two
precisions plus two Yukawa pieces and the per-L referee, about 30-50 min, at most two processes.

## 5. Controls (measured, `debug/data/l52_controls_*.json`)

- K1 symbolic V coefficients (above). K2 constant = corpus K, 20 digits.
- K3 equal arms (E3): segment formula 1.3356449576164494005e-6 = AE2b's closed form (rel 4e-20); amp x T_L
  = AE2b's pred_F at L = 220, 400, 580 to 1e-15.
- K4 NA's 3/2 amplitude at unequal arms (p = -1, weight 1 - u), exact F_L, 2000 digits: R_240 = 0.99616 (N1),
  0.99508 (N2); four-term extrapolation 0.99999987, 0.99999991. The unequal-arm segment formula works in the
  one-electron case to 1e-7 (this amplitude was never computed by record NA, which scored exponents).
- K5, K6 instrument toys (7/2 and 9/2 with nuisances up to 20): extrapolation error <= 1.2e-5 and <= 1.1e-4. A
  three-term fit was 10x worse (1.7e-3) and was replaced before any scored data existed.
- K7 the engine's Y at mass kappa runs (mod p, L <= 3) and its per-L referee agrees to 1e-51 at A1, L = 0, 1
  (Y_kappa alone; D not formed).
- Gate fire test 14/14; the first version let a missing case read BORDERLINE (fixed: missing = STOP).

AB's 7/2 at A1/A2 was not computed as a control on purpose: its amplitude is P1's content.

## 6. Recommendation

Do **not** register L52 as a test of the count: the count, done right, is the corpus's mechanism (AB P3b,
`aaeq2_local.py`), the memo's 5/2 was a miscount, and the function-level vanishing is derived and already scored.
Suggested one-line capture (PI call, no edit made): Paper 59's (AA|BC) paragraph could say why 5/2 is absent from f
(the potential of a bounded density has no s^{-1} or s^{1} term), and the lit memo's Q1 row should be corrected.

The slim amplitude record drafted here (P1 + P2, P(GO) about 0.72, under an hour) is worth registering only if the
PI wants row 6's singular-point cell to carry a measured amplitude at unequal arms, or wants the kappa-independence on
record. The genuinely open item is operator-level: why MO's operator has 3/2 but no 5/2 at x*. That needs a
per-component analysis, not this record.
