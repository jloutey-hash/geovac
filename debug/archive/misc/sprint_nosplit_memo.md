# Sprint nosplit: the chain's K3 piece does not split off (2026-10-17)

**Verdict: GO.** The bounded search of D12 is now a proof. The chain operator
W = R3 o R1 o (theta - 2) does not split: Sol(W) has no 3-dim sub mapped onto Sol(R3).
Two independent exact formulations, each searched exhaustively inside a proved bound, both
come out inconsistent; four planted-split controls and one natural split are all found by the
same code. New, and not expected: **the chain satisfies an order-4 right factor of W**. The
artanh piece of R1 o (theta-2) is absent from C (measured: the coefficient is 9e-45 at 40
digits and falls with precision; a perturbed control gives 2e-2). So C's minimal operator has
order 4, R3' o (theta - 2), and C is a K3 function integrated once:
C = c t^2 + t^2 int^t s^-3 (B' f)(s) ds, with f a solution of F_3's operator L3 and B' an
explicit rational map. That integral cannot be done inside the K3 module (the no-split).

Scripts: `debug/nosplit_{core,local_logs,search,order4,membership,remaining,maps}.py`.
Data: `debug/data/nosplit_*.json`, logs `debug/data/nosplit_{controls,order4,remaining}.log`.

## 1. Set-up

Work in x = t^2, delta = x d/dx (theta = 2 delta). All operators are even in t, so this loses
nothing (see step 3). E = R1 o (theta - 2) (order 2, solutions t^2 and
g2 = t^2 [t/(2(1-t^2)) + artanh(t)/2]). W = R3 o E.

**Split, as an equation.** The extension 0 -> Sol(E) -> Sol(W) -E-> Sol(R3) -> 0 splits as
differential modules over Q(t) iff there is S = s0 + s1 delta + s2 delta^2 with rational s_k and

    E o S = 1  (mod R3, right remainder).

(S is the section: E(S f) = f for f in Sol(R3).) The condition is affine-linear in S, so if a
solution exists over C(t), the average of S(t) and S(-t) is an even solution, i.e. one over
C(x); and the resulting linear system has rational coefficients, so a complex solution gives a
rational one. Searching over Q(x) is therefore complete.

## 2. Route 2 first: local monodromy does not obstruct

Exact generalised-Frobenius log counts (`nosplit_local_logs.py`, `nosplit_membership.py`):

| point | W | R3 | E | consistent with a split? |
|:--|:--|:--|:--|:--|
| t = 0 | group {1/2,3/2,3/2,5/2}: dim 4, max log 2 | MUM, log 2 | log-free | yes (Jordan 3+1) |
| t = inf | group {-1,0,1,2}: dim 4, max log 2 | log 2 | log-free | yes |
| t = +-1 | exps -1..3, max log 1 | ordinary | log 1 (artanh) | yes (rank-1, image t^2) |
| t = +-2, 4t^2+11 = 0 | E ordinary there, so lifts are local | | | yes |

The extension is locally split everywhere. Any obstruction is global, so Route 1.

## 3. Route 1: the bound (proof)

Let f_1..f_3 be a local basis of Sol(Top) at p (Top = R3 here), Phi the matrix of
delta^k f_i (or d^k f_i at finite p), y_i = S f_i in Sol(W). Then (s_k) = Phi^{-1} y.
With v the valuation (lowest real exponent; it is a valuation on Nilsson-class series):

- det Phi has valuation sum(rho) (delta-form) or sum(rho) - 3 (d-form) exactly (Fuchs);
- each cofactor has valuation >= (sum of two exponents) - (derivative orders), since
  v(d^j f) >= v(f) - j and v(delta^j f) >= v(f);
- v(y_i) >= mu_W(p), the least exponent of W at p.

Hence, with rho_max = largest exponent of Top at p:

    finite p != 0, Top singular:  ord_p s_k >= ceil(mu_W - rho_max + k)
    finite p, Top ordinary:       ord_p s_k >= ceil(mu_W)      (Phi^{-1} holomorphic)
    p = 0 and p = inf:            ord_p s_k >= ceil(mu_W - rho_max)

(the delta <-> d change is triangular with entries in Q[x, 1/x], harmless at p != 0). Poles
can sit only where W or Top is singular, apparent points included; the bound at infinity caps
the numerator degree. W and Top are Fuchsian everywhere (checked: indicial degree = order at
every point), so no irregular case arises. All exponents are exact (sympy indicial
polynomials; at irrational apparent points, exact arithmetic in Q[x]/(q), cross-checked
against the rational path at x = 1, 4, -11/4).

**The chain, full question (Top = R3, Bot = E):**

| p (x) | exps W | exps R3 | bound (s0,s1,s2) |
|:--|:--|:--|:--|
| 0 | 1/2, 1, 3/2, 3/2, 5/2 | 1/2, 3/2, 5/2 | -2, -2, -2 |
| 1 | -1, 0, 1, 2, 3 | ordinary | -1, -1, -1 |
| 4 | -3/2, 0, 1, 2, 3 | -7/2, 0, 1 | -2, -1, 0 |
| -11/4 | 0, 1, 2, 3, 5 | 0, 1, 3 (apparent) | -3, -2, -1 |
| inf | -1, 0, 1/2, 1, 2 | -1, 0, 1 | -2, -2, -2 |

So s_k = N_k/D_k with D_0 = x^2 (x-1) (x-4)^2 (4x+11)^3, D_1 = x^2 (x-1)(x-4)(4x+11)^2,
D_2 = x^2 (x-1)(4x+11), deg N_k <= 10, 8, 6: **27 unknowns**. In t this is denominators
t^4 (t^2-1) (t^2-4)^2 (4t^2+11)^3 and numerators up to degree 20.

**Search:** E o S mod R3 for each basis monomial, exact, then the polynomial identity
E o S = 1 coefficient by coefficient over QQ (DomainMatrix). **rank A = 27, rank [A|b] = 28:
inconsistent. No split.** Padding every pole order and degree by 2 (57 unknowns): rank
57 vs 58, still none. Run time 100 s.

**Independent cross-check (step 5).** Because the partial split exists (step 4), the full
question is equivalent to the split of the order-4 factor R3' o (theta - 2), with
R3' = ann(S1(Sol R3)) (apparent singular points x = 1 and 20x^2 - 7x + 5 = 0, no longer
-11/4). Same proof, different operators: bound gives 33 unknowns, rank 33 vs 34;
padded by 2, 69 vs 70. No split.

**Controls (same code path, bounds computed the same way):** planted splits
Top = ann(Bot(Sol L3)), so W_c = lclm(Bot, L3):

| Bot | singular points of Top | unknowns | found and verified exactly |
|:--|:--|:--|:--|
| delta | 0, 4 | 7 (25 padded) | yes, S = [-6, -12, -8(x-4)/x] |
| delta^2 | 0, 4 | 7 (25) | yes |
| R1 | -3, 1, 4, roots of x^2-10x-3 | 28 (70) | yes |
| E (the chain's own bottom) | 1, 4, roots of a quartic | 34 (82) | yes |

The natural split of step 4 is also found, and its pole order 2 at x = 1 meets the bound
exactly (the bound is tight there). D12's control (reversed composition) stands as before.

## 4. The order-4 question

**R3 o R1 splits** (exact, `nosplit_search.py chain`): S1 =
[(12x^3+229x^2-31x+12)/(3x(x-1)^2(4x+11)), -2(20x^3+619x^2-143x+80)/(15x(x-1)^2(4x+11)),
-4(x-4)(20x^2-7x+5)/(15x(x-1)^2(4x+11))] satisfies R1 o S1 = 1 mod R3, verified exactly.
Consequences:

- Sol(W) has exactly one 4-dim invariant subspace, U = {y : (theta-2) y in S1(Sol R3)};
  it contains t^2 and not g2. Its operator is P = (theta - 2) - S1 o E (order 4), a left
  multiple of R3' o (theta - 2), and **W = Q1 o P exactly**,
  Q1 = 15x(x-1)(x+3)(4x+11)/((x-4)(20x^2-7x+5)) + 30x(x-1)^2(4x+11)/((x-4)(20x^2-7x+5)) delta.
- For y in Sol(W), P y = alpha(y) h with h = t^3/(1-t^2)^2 and alpha = the g2-coefficient.

**Does C satisfy P?** alpha(C) = (P C)/h on the chain integral (step B2's evaluator):

| dps | t = 0.5 | 1.3 | 2.5 | 3.5 | control (S1 coeff x 1.01) |
|:--|:--|:--|:--|:--|:--|
| 20 | 2.7e-26 | 1.7e-24 | -4.0e-25 | 1.2e-23 | -0.020, -0.045, -0.062, -0.066 |
| 30 | 7.3e-33 | 3.7e-35 | -1.2e-34 | 7.5e-34 | same |
| 40 | | 8.9e-45 | | | |

alpha must be one constant at every t; it is zero to working precision and falls ten orders
per ten digits. **So C lies on the order-4 factor: P C = 0** [MEASURED, 40 digits; not
proved]. An exact proof via analyticity at t = 1 was tried and does not close: P itself has a
log at x = 1 (max log degree 1), so "single-valued at t = 1" does not single out U.

Lattice of invariant subspaces of Sol(W): 0, <t^2>, Sol(E), U, Sol(W). C is in U and not in
Sol(E), and U has no 3-dim K3 sub (the no-split). **C's minimal operator is P, order 4.**
D11's statement "W annihilates C" stands; W is not minimal.

What C is: (theta - 2) C lies in S1(A(Sol L3)) = B'(Sol L3), with B' = (S1 o A) mod L3 =
[-t^3(t^4-7), -t^3(t^2+5)(2t^2-5), -t(t^2+2)^2 (t^2-1)] / ((t^2-4)(t^2-1)^2) on
(f, theta f, theta^2 f) (common denominator). So

    C(t) = c t^2 + t^2 int^t s^-3 (B' f_C)(s) ds,   f_C in Sol(L3),

and the no-split says this primitive is not B'' f + c' t^2 for any rational B''.

## 5. What this means for the registered outcomes

- **(A) C = (R0 + R1 D + R2 D^2) F_3: refuted exactly.** It would put C in a 3-dim K3
  sub of Sol(W), which is a split.
- **(B) C = B F_3 + E, E elementary: refuted exactly**, in a stronger form: for any rational
  B (of any order) and any holonomic E whose D-module has no constituent equivalent to L3.
  Proof: D C = U sits in M_B (+) M_E, M_B = B(Sol L3) irreducible. If U contains M_B, U
  splits; if not, U meets M_B trivially and embeds in M_E, which has no K3 constituent. Both
  are impossible. (Uses: L3 = 3F2(1/2,1/2,1/2;1,1;x/4) is irreducible by Beukers-Heckman, and
  stays irreducible over Q(t) because its monodromy is Zariski dense in SO(3), which is
  connected, so the index-2 subgroup is dense too.)
- **The right shape is (C): a K3 function integrated once against an elementary kernel,**
  the extension of the K3 piece by t^2 being non-split. The artanh piece does not occur (40
  digits). This matches the memo's reading "K3 reached through integrals with elementary
  kernels" and makes it precise: one integral, one elementary constituent.

## 6. What a test should assert

1. Exact: S1 above gives R1 o S1 = 1 mod R3; W = Q1 o P with zero remainder.
2. Exact: the bound table of section 3 (exponents of W and R3 at 0, 1, 4, -11/4, inf), and the
   27-unknown system E o S = 1 mod R3 has rank 27 < 28. Fire-test: the same function on
   Top = ann(E(Sol L3)) must find S (rank 34 = 34). Rejects "a split exists" and a bound
   that shrinks when an exponent is miscomputed.
3. Exact cross-check: R3' o (theta - 2) system rank 33 < 34.
4. Slow numeric: |P C / h| < 1e-25 at dps 30 at two t; the control with S1's first
   coefficient times 1.01 gives > 1e-3. Rejects "C needs the artanh piece / order 5".

## 7. Caveats

- The order-4 membership is measured, not proved (section 4).
- The bound argument is the standard Wronskian-valuation one, written out above; its only
  inputs are exact exponents. A reader should check the cofactor inequality at the apparent
  point -11/4, where it is crude (allows pole order 3; the planted controls at apparent points
  find poles of order 1).
- Nothing outside debug/ was changed.
