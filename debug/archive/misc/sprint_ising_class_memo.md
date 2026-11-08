# Sprint memo: the Ising-class moments on the period ladder (IS1-IS3)

Date: 2026-10-17. Contract: `docs/predictions/2026-10-17_ising_class.md` (commit 169a6ce8).
Object: c_{n,k} = int_0^oo x^k K0(x)^n dx; c_{n,0} = F_n(1) (Paper 59 sec:four_density).
Everything here lives in `debug/`. No commits.

## Verdicts

| Prediction | Verdict | Deciding evidence |
|:--|:--|:--|
| IS1 (parity rule) | CONFIRMED | Operator part was already done. Here: operators for n = 6, 7, 8 annihilate F_n by quadrature to <= 1.1e-42 (relative) at t = 0.5, 0.7, 1.3, 1.7. Perturbed controls: 6.8e-7 to 9.8e-2. No source term. |
| IS2 (odd n stay on rung) | CONFIRMED (within bounds) | c30 not Tate (bound 1e32). c50, c70: no relation with Tate, or with c30 x Tate. Bounds 1e17/1e14 (c50) and 1e9/1e7 (c70). All controls pass. |
| IS3 (even n drop) | CONFIRMED for c40 and c60, as registered | c40 is not Tate (bound 1e24). c40 = 4 pi^2 L(8.4.a.a, 2). c60 = 108 pi^2 L(24.6.a.a, 4). Both found by PSLQ and checked to 4.9e-261. |
| Not registered: n = 8 | The rule did not anticipate this | c80 shows no drop: no Tate relation (1e7), and no relation with any L-value tried (1e12 homogeneous). The literature puts c80 inside a 2x2 determinant of a rank-4 motive. |

The value-level rule is not a parity rule. In Lim-Tu-Yu's structure, c_{n,0} is a single L-value when the
period matrix D_{n,even} is 1x1. That holds for n = 1, 2, 3, 4 and for n = 6 (via D'_6). It fails for n = 5, 7, 8.
So n = 3 (odd) is a single L-value, and n = 8 (even) is not. The registration's heuristic matches the data at
n <= 7 only by coincidence of the cases it named. I state this plainly; I do not reinterpret IS3.

## 1. IS1 remainder (quadrature on F_n)

Script: `debug/ising_class_is1_quad.py`. Data: `debug/data/ising_class_is1_quad_6_7_8.json`. dps 45.
theta acts only on K0(tx), as in the n = 4 test. Control: the top-t-power constant coefficient gets +1.

| n | t = 0.5 | 0.7 | 1.3 | 1.7 | perturbed range |
|:-:|:-:|:-:|:-:|:-:|:-:|
| 6 | 5.7e-45 | 2.5e-44 | 3.4e-43 | 9.0e-43 | 8e-4 to 1e-1 |
| 7 | 4.2e-45 | 1.9e-44 | 3.2e-43 | 1.1e-42 | 7e-5 to 2e-2 |
| 8 | 6.0e-46 | 5.0e-45 | 2.3e-43 | 1.1e-42 | 7e-7 to 1e-3 |

That is >= 41 digits, well past the 25 the contract asks for. The operators are homogeneous: no source term.

## 2. Moments, two independent routes

- Route A (`ising_class_moments_mp.py`): mpmath tanh-sinh. K0 is our own code: the ascending series for x <= 2,
  and the trapezoid rule on int e^{-x cosh t} dt for x > 2. It was checked against mpmath.besselk at 7 points; the
  difference is 0 at 300 digits. Every quadrature error estimate is <= 1e-304.
- Route B (`ising_class_moments_arb.py`): on [0,1], exact termwise integration of the log-series. On [1,R],
  rigorous python-flint arb integration. Arb radii are <= 1e-350 for n >= 2.
- Agreement (`ising_class_compare.json`): >= 286 digits for every (n,k), n = 2..8, k = 0..3. That is the storage
  limit of Route A.
- n = 1: Route B gets only 93-101 digits (see tool bugs). Route A matches the closed forms to 290 digits.
- Closed forms reproduced: c10, c11, c12, c13, c20, c21, c30, c31, c32, c33, c41, c43. Route A matches them to
  >= 289 digits and Route B to >= 319.

c_{n,0}, n = 1..8: 1.5707963267948966192, 2.4674011002723396547, 6.9488227810796297894, 27.241338417805973407,
135.26830258086883759, 809.62084822486627594, 5662.3426897177563410, 45284.272554984611744.

## 3. Literature (primary sources downloaded to `debug/ising_class_lit_downloads/`, untrusted, read as text)

VERIFIED against the primary text and against our numbers:
- BBBG, arXiv:0801.0891, Table 1, proven:
  - c30 = 3 G(1/3)^6/(32 pi 2^{2/3}) and c31 = (3/4) L_{-3}(2); both match to >= 290 digits.
  - c32 and c33 = L_{-3}(2) - 2/3 match.
  - c40 = (pi^4/4) 4F3(1/2,1/2,1/2,1/2; 1,1,1; 1) matches to 60 digits (the series converges slowly).
  - c41 = 7 zeta(3)/8 and c43 = 7 zeta(3)/32 - 3/16 match.
  - For c50 and c60, BBBG give double integrals only; there is no closed form.
  - BBBG Conjecture 11 (c80 = pi^2 int I0^2 (6K0^2 - pi^2 I0^2) K0^4) holds to 6.6e-36 here.
- Broadhurst, arXiv:1604.03057, sec. 7. His L-series attach to odd-power moments with I0 factors, S_{n,s}, not to
  c_{n,0} for n >= 5. The forms are f_{3,15}, f_{4,6} = (eta1 eta2 eta3 eta6)^2, f_{6,6}, and L_7 (weight 3,
  level 525, quartic nebentypus).
  - Checked to 260 digits: S_{4,3} = 7 zeta(3) (that is c41), S_{3,2} = 3 L_3(2) (that is c31),
    L_5(2) = (2 pi/sqrt15) L_5(1), and L_6(3) = 2 zeta(2) L_6(1).
  - His conjecture L_8(5)/(zeta(2) L_8(3)) = 4/7 holds to 1.1e-260.
- Rogers-Wan-Zucker, arXiv:1303.2259. These were used to validate the L-value code. Three Gamma-value closed forms
  at s = 2 match to <= 3.3e-261: eta(4t)^6, eta(2t)^3 eta(6t)^3, and (eta3 eta5)^3 + (eta1 eta15)^3. The relation
  L(f2,3) = (pi/2) L(f1,2) = (pi^3/32) 4F3, with f2 = eta(2t)^4 eta(4t)^4, is consistent with c40.
- Lim-Tu-Yu, arXiv:2609.13576 (Sept 2026; quadratic-twisted Kloosterman moments). This is the source that owns
  c_{n,0}.
  - c30 = (pi^2/2) L(48.3.e.a, 1) and c40 = 4 pi^2 L(8.4.a.a, 2) are stated as equalities. We reproduce both to
    3.3e-260 and 3.2e-261.
  - c60 =? 108 pi^2 L(24.6.a.a, 4) is numerically checked by the authors. We reproduce it to 4.9e-261.
  - Theorem thm:EC parts (i) and (ii) hold here for every prime 3 <= p <= 619.
  - Part (iv): we rebuilt a_p of 24.6.a.a from tm_6(p). The paper prints a(5) = -34 and a(7) = -240; ours match.
    The functional equation then forces a_3 = -9 and eps = -1 (it closes to 3.9e-31; the other choices fail at
    1e-2). This agrees with their statement that L(24.6.a.a, 3) = 0.
  - Their Deligne table puts c50, c70 and c80 only inside 2x2 determinants (rank 3, 4 and 4 motives).

NOT VERIFIED:
- LMFDB data for 24.6.a.a and 60.3.b.a. Not downloaded; 24.6.a.a was rebuilt from Kloosterman sums.
- Evans's original papers; the AMNT and Zagier proofs cited by LTY.
- The registration's reading "attached to the conifold at x = 1/4". We find a weight-4 form, 8.4.a.a. We did not
  test whether it is the motive of the conifold fibre.
- Chudnovsky's algebraic independence of G(1/3) and pi. It is quoted from memory; with it, "c30 is not in Q(pi)"
  is a theorem.

## 4-5. PSLQ results (`ising_class_pslq.py`; data `debug/data/ising_class_pslq_{is2,is3,is2L,is3L,hom}.json`)

Settings:
- tol = 10^-(dps-40). Coefficient bound 10^min(35, floor(0.7(dps-40)/dim)).
- "none" means mpmath ended on its norm bound, which is at least the stated bound.
- A found relation is accepted only if its residual at full input precision is < 10^-(dps-5).
- Moments are stored to 330 digits; L-values to 260.

Relations found (all with residual <= 4.9e-261):
- c30 = (pi^2/2) L(48.3.e.a, 1)
- c40 = 4 pi^2 L(8.4.a.a, 2)
- c60 = 108 pi^2 L(24.6.a.a, 4)

No relation (target: basis, dim, bound):
- c30: Tate wt<=3 (6, 1e32).
- c50:
  - Tate wt<=5 (11, 1e17); + c30 x Tate wt<=2 (14, 1e14).
  - Extended with pi^a/c30 (22, 1e8).
  - Mixed, + pi^a L_5(2) and sqrt15 pi^a L_5(2) (19, 1e7); + Broadhurst C and 1/C (24, 1e6).
  - Homogeneous weight 5, Tate_5 + c30 Tate_2 + pi^{5-s} L(f,s) (7, 1e21), for f = 8.4.a.a, 6.4, 6.6, 24.6.a.a,
    15.3.
- c70:
  - Tate wt<=7 (20, 1e9); + c30 x Tate wt<=4 (27, 1e7); extended (37, 1e5).
  - L_5 and C mixed (32-37, 1e4).
  - Homogeneous weight 7 (10, 1e15), same five forms.
- c40:
  - Tate wt<=4 (8, 1e24).
  - Mixed with 6.4, 48.3.e.a and 15.3 L-values (16, 1e9). These runs show that 8.4.a.a is singled out.
- c60:
  - Tate wt<=6 (15, 1e13).
  - Mixed 6.4, 8.4.a.a and 6.6 (23, 1e6); homogeneous (9, 1e16).
- c80:
  - Tate wt<=8 (26, 1e7).
  - Mixed 6.6, 24.6.a.a, 6.4 and 8.4.a.a (34, 1e4); homogeneous weight 8 (12, 1e12).

Controls:
- Every mixed run recovers c20 = pi^2/4 and c41 = 7 zeta(3)/8 at its own precision and bound: 27/27 each.
- Every run recovers a planted relation with all coefficients nonzero, of height bound/100: 47/47 informative
  runs (27 mixed + 20 homogeneous).
- Homogeneous runs carry the planted control only. The c60 rediscovery is their positive control.
- The 4 homogeneous runs against 48.3.e.a are void: L(48.3.e.a, 1) is c30 x pi^-2 x 2, so they returned the
  basis-internal relation, and their planted control did too.

Not done, for feasibility (STOP):
- L_7 (newform on Gamma_0(525), quartic nebentypus, coefficients in Q(i, sqrt6, sqrt14)) needs modular symbols.
- L(chi Sym^2 g, s) for g = 60.3.b.a is degree 3, with unknown bad Euler factors and conductor.
- Neither is a q-expansion we can build. LTY attach both only to determinants, not to c50 or c70 alone.

## Transcendentals, by type

| Transcendental | Type | Where it occurs |
|:--|:--|:--|
| pi | Tate, weight 1 | c10, c20 |
| zeta(3) | Tate, odd zeta | c41, c43 |
| L_{-3}(2) | Dirichlet; mixed Tate over Q(sqrt-3); outside the registered Tate basis | c31, c33 |
| G(1/3)^6/pi | CM elliptic period (Chowla-Selberg), = a weight-3 CM L-value | c30, c32 |
| 4F3(1/2^4;1^3;1) = L(8.4.a.a, 2) x 16/pi^2 | weight-4 non-CM modular L-value | c40 |
| L(24.6.a.a, 4) | weight-6 modular L-value | c60 |
| G(1/4)^4/pi and the G(k/15) product | CM | validations only |

c50, c70 and c80 are unidentified.

## Caps

- Single-scale values only.
- PSLQ "none" is evidence at the stated dim and bound, not proof. The c70 and c80 mixed bounds are weak (1e4 to
  1e9); the homogeneous runs (1e12 to 1e21) carry the weight.
- The c60 relation rests on our rebuild of 24.6.a.a. Its checks are the printed a5 and a7, Ramanujan bounds and
  integrality at p <= 619, and functional-equation closure.
- IS3's "drops below its rung" has no operational form beyond the bases listed.

## Tool bugs found

1. mpmath.besselk costs about 0.45 s per call at 300 digits, which is impractical. Route A uses its own K0.
2. python-flint acb.integral on the K0^1 x^k integrand over [1, 730] stalls at radius about 1e-98, even subdivided.
   Every n >= 2 reaches 1e-350. Not diagnosed.
3. My PSLQ controls were wrong twice before the final runs:
   - A planted value rounded to run precision fails for a reason that has nothing to do with PSLQ.
   - tol = 10^-(dps-10) makes relations above height about 1e10 undetectable, so the norm bounds it printed above
     that were meaningless.
   - An earlier default-tol run with a too-generous bound returned a spurious c60 "Tate relation" (coefficients
     about 1e15). The full-precision residual check now rejects that class.
4. The first homogeneous bases contained FE-related L-values, which produced internal relations. Fixed by using
   one value per FE orbit.
5. The first IS1 run hit my own 1800 s timeout. The rerun memoises Bessel values.
