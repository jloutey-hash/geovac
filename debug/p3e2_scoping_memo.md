# P3e-2 scoping: the generic three-electron chain fibre (2026-10-19)

**Status:** scoping only. No PSLQ of the generic fibre G was run against any new or widened
basis, and no new value of G was computed. One existing script was run: the positive control
`debug/three_electron_p3e_ctrl_exact.py` (the 150-digit reduction of Gctrl, from cached values).
It reproduced the recorded result: Gctrl = (11050792 V - 27272112 D1 - 52735824 D2 + 87003024 S12)/363825,
residual -1.57e-149, in 0.2 s. Its output file came out byte-identical. The per-call cost of K0, K1
and the Struve-based L was also timed, with no fibre integrand evaluated. Draft pre-registration:
`debug/p3e2_prereg_draft.md`.

## 1. State, checked against the files

- **Object.** The true D = 0 fibre of the chain <rho_a rho_b rho_c / (r12 r13)>. Here a is the
  central density, the one carrying both Coulomb kernels; b and c are the leaves. Up to a prefactor
  the fibre is G(a;b,c) = int_0^inf x K1(ax) J(bx) J(cx) dx, with J = 3L - 3XK0 - X^2 K1 and
  L(X) = int_0^X K0. G is symmetric under b <-> c (`docs/predictions/2026-10-13_three_electron_k3.md`, P3e).
- **What is known.**
  - At equal scales: 216 G(1,1,1) = 80 c30 - 6600 c32 + 81 pi^3.
  - At (1;2,4): the K3 block {V, D1, D2, S12} has rank 4. Gctrl, the part of G with no L in it,
    reduces exactly onto that block.
  - G itself has no relation that is stable from 260 to 320 digits. The bases tried were:
    - K3 block + pi^k: no relation below about 1e24.
    - adding pi times the residue-curve periods and quasi-periods: none below 1e12.
    - adding log 2, then log 2 and log 3: none below 1e9 to 1e13.

  (CHANGELOG v5.20.15; Paper 59 `sec:three_electron`, [OBSERVATION, evidence]; walls register,
  PERIODS-TRANSCENDENCE, member 7, OPEN-LEANING; ledger row P3e-2.)
- **What the registry entries do not yet reflect.** Paper 59 `sec:density_ring` (v5.20.21,
  2026-10-18) is newer than the walls entry of 2026-10-16. It settled the function-level class of a
  closely related object, the Coulomb chain C(t) with all three densities set to the fast-decay root
  Delta^{-3} and scales (1,1,t), t on the central density.
  - Its minimal annihilator is R3 R1 (theta - 2). R3 is equivalent to F_3's K3 operator; R1 and
    (theta - 2) are elementary.
  - The chain satisfies an order-4 factor P. On the slice this gives
    C = c t^2 + t^2 int^t s^{-3} (B'f)(s) ds, with f a K3 solution: "a K3 function integrated
    once against an elementary kernel".
  - The K3 piece provably does not split off: an exponent bound plus an exhaustive exact search.

  So on that slice the fibre is not a rational combination of K3 periods plus elementary terms, as
  a function. A null PSLQ against K3 periods + Tate numbers is exactly what this predicts. The
  P3e fibre uses the Slater weight Delta^{-5}, not Delta^{-3}. The two are related by mass
  derivatives (Delta^{-5} = Delta^{-3} + (2c/3) d/dc Delta^{-3}), but derivatives in the leaf masses
  leave the slice. So the transfer to G is an assumption, not a result.
- **What the old bases actually covered.** Gaps the record does not list:
  1. **(1,2,4) is a geometric progression.** P_uv = pi K((v-u)/(v+u))/(u+v) is homogeneous of
     degree -1, so P_bc(2,4) = P_ab(1,2)/2 and Q_bc(2,4) = Q_ab/4 exactly. The "(b,c) residue
     pair" that the record lists as omitted is therefore redundant at that triple, so that
     particular caveat had no content there. The K3 block's measured rank 4 (Picard rank 18) shows
     no other specialness.
  2. **No degree-3 Tate numbers.** The Tate part of every basis was {1, pi, pi^2, pi^3, log p,
     pi log p}. Count one degree for each factor of pi, K, E or log, and for each unit of polylog
     weight. Then every term of both exact relations has degree 3: the equal-mass collapse (c30
     and c32 are pi K K'-type, plus pi^3) and the Gctrl reduction. Yet pi^2 log p, zeta(3),
     log p x (elliptic periods), pi Li2 and Li3 were never in any basis. The degree count is a
     heuristic, not a motivic weight (Legendre's relation K E' + E K' - K K' = pi/2 already mixes
     degrees). Still, it says which numbers were missing.
  3. **Products of Bailey moduli carry a Galois twist.** At (1;2,4) the Bailey moduli are
     k+- = (sqrt21 +- sqrt5)/8. The individual products K(k-)K(k+') and K(k+)K(k-') are conjugate
     under sqrt5 -> -sqrt5. A Q-linear PSLQ can only see them through invariant combinations: the
     sum (which is V) and sqrt5 times the difference. A basis listing the raw products without
     that care would give false nulls.

## 2. Prior art

**Verified in the primary source this session:**
- Bonisch, Fischbach, Klemm, Nega, Safari, *Analytic structure of all loop banana integrals*,
  JHEP 05 (2021) 066, arXiv:2008.10574 (full text read).
  - The abstract treats "all banana amplitudes with arbitrary masses" on Calabi-Yau
    (l-1)-folds.
  - Introduction, p. 3-4: the Feynman integral is a relative period, because the integration
    domain has a boundary. Quoted: the "(full) Feynman integral will be a linear combination of
    the homogeneous solutions (a Frobenius basis for the closed periods) and a special
    inhomogeneous solution (regarded as a special relative period)."
  - Sec. 3.7: closed values are listed only at special points of the equal-mass integrals, at the
    on-shell point s = 1. Examples: F_sigma2(1) = (1 - 2^-2) zeta(2) = pi^2/4, and
    F_sigma3(1) = (12 pi/sqrt15) L(f3, 2) with f3 of weight 3 and level 15.
  - At generic kinematics the integral is given as a function (Frobenius data, Gamma-class and
    zeta-value coefficients). No closed value is given.
- Bloch, Kerr, Vanhove, *A Feynman integral via higher normal functions*, Compositio Math. 151
  (2015) 2329, arXiv:1406.2664 (abstract). The three-loop equal-mass banana in two dimensions is a
  sum of elliptic trilogarithms at sixth roots of unity. At a special kinematic point it equals a
  critical value of a K3 surface's L-function (Broadhurst's conjecture, proved there).
- Adams, Bogner, Weinzierl, arXiv:1405.5640 (abstract): the two-loop sunrise in D = 2 with
  arbitrary masses, in elliptic dilogarithms whose arguments change with the masses.
- BBBG 2008, arXiv:0801.0891 (abstract only): equal-argument moments c_{n,k}. Eq. (33) for three
  different arguments has been checked numerically in the corpus to 1e-34. I did not re-read it in
  the primary source.

**Corpus (read):**
- Paper 59 `sec:four_density` and `debug/sprint_banana_dictionary_memo.md` (exact arithmetic):
  - The chain's K3 is not the three-loop banana's K3. F_3 = (A) * bubble is the Clausen square of
    binom(2m,m)^3. The three-loop banana's K3 is W_4 / Verrill.
  - F_n is the n-propagator two-dimensional banana integrated over a line of external momentum,
    via int_0^inf J0(px) dp = 1/x. The same one-line identity gives
    V3(a,b,c) = int_0^inf dp S(p; a,b,c), with S the D = 2 unequal-mass sunrise. I derived this
    and did not check it numerically.
  - So the chain's holomorphic period is the unequal-mass sunrise integrated over momentum. The
    true fibre G adds the two Coulomb integrations, the L terms, on top.
- Paper 61 (cited in the walls register): one rung down, the analogous two-electron relative object
  T2 has no relation through weight three, Catalan-completed, at 64 digits.

**Not verified:** Davydychev-Tausk, Nucl. Phys. B 397 (1993) 123. A search did not confirm the
citation or its Clausen-function content. The expectation that the D = 2 vacuum sunrise
int x K0(ax)K0(bx)K0(cx) dx is a dilogarithm combination with Kallen-root arguments is therefore
mine, unverified.

**Your recollection, checked:**
- Bloch-Kerr-Vanhove: right. The three-loop banana is K3-related; its full integral is elliptic
  trilogarithms.
- Bonisch et al.: right. Unequal masses are treated as Calabi-Yau families.
- Correction: the chain is not the three-loop banana family. It is on the same rung, but it is the
  sunrise integrated over momentum.

**Does the literature settle P3e-2?** No.
1. No source treats this object, the sunrise integrated over momentum with two Coulomb
   integrations. The corpus searches found none, and neither did I in the sources above.
2. For bananas themselves, the literature gives the generic-mass full integral as a function: closed
   periods plus one special relative period. No source says whether its value at a generic rational
   point satisfies a Q-relation with closed periods and polylogarithms. That is a period-conjecture
   question, not a theorem.

What the literature does give is a strong expectation, by analogy, for both outcomes:
- At generic masses the fibre is a genuine relative period, and a null PSLQ is the expected result.
  It is a check of that picture, not a finding.
- Closed values occur at special points (on-shell or threshold) as zeta values or L-values. That
  points to threshold mass triples as the place where a collapse could be found. A collapse there
  is not in the literature for this object.

## 3. The widened basis, concretely

Conventions: a triple (a;b,c), with a the central density.

**Letters.** Take the nonzero values of a, b, c, a+-b, a+-c, b+-c, a+b+c, -a+b+c, a-b+c, a+b-c.
- Prime set Pr: the primes dividing them.
  - (1;2,4): Pr = {2, 3, 5, 7}. (3;4,6): Pr = {2, 3, 5, 7, 13}.
  - (1;3,7) and (3;1,7): Pr = {2, 3, 5, 7, 11}.
- Dilogarithm letters Lambda: the residue moduli k_uv = |v-u|/(v+u), and the threshold ratios
  (-a+b+c)/(a+b+c), (a-b+c)/(a+b+c), (a+b-c)/(a+b+c). Keep the distinct values with 0 < |l| < 1.
  At (1;2,4) this gives {1/3, 3/5, 5/7, 3/7, -1/7}.
  These are the PI's "mass ratios" in the form the problem's singular loci produce. Bare ratios
  like 1/2 give Li2(1/2), which is reducible.

**Blocks.** Counts are at (1;2,4).

| Block | Elements | n |
|:--|:--|:-:|
| K3 | V, D1, D2, S12. The two internal relations are removed; the rank is re-measured at each triple. | 4 |
| B, Bailey products | kappa = pi/(2 c_max). For X, Y in {K, E}: S_XY = kappa [X(k-)Y(k+') + X(k+)Y(k-')] and A_XY = kappa sqrt(d) [X(k-)Y(k+') - X(k+)Y(k-')], with sqrt(d) the root that k+ <-> k- negates (sqrt5 at (1;2,4); sqrt3 at (1;3,7)). For complex-conjugate moduli, (3;4,6), use 2 Re and 2 Im. S_KK = V. | 7 new |
| R, residue periods | pi P_ab, pi Q_ab, pi P_ac, pi Q_ac (as before) | 4 |
| RR, products of the two residue curves' periods | P_ab P_ac / pi, P_ab Q_ac / pi, Q_ab P_ac / pi, Q_ab Q_ac / pi | 4 |
| T3, degree-3 Tate | pi^3, zeta(3), pi^2 log p (p in Pr) | 6 |
| Li, polylogarithms | pi Li2(l), Li3(l) for l in Lambda | 10 |
| O, vacuum sunrise | U = int x K0(ax)K0(bx)K0(cx), U_a = int x^2 K1(ax)K0K0, U_b, U_c (odd parity; the dilogarithms with their true algebraic letters, so they need no guessing) | 4 |
| LR | log p x {P_ab, Q_ab, P_ac, Q_ac} | 16 |

The PI's two named candidates appear as follows:
- Dilogarithms of mass ratios: blocks Li and O.
- Products of elliptic periods: blocks B and RR.

**Sub-spans.** n counts the lead number (G, Gctrl or the planted control).

| Span | Contents | n |
|:--|:--|:-:|
| SA, elliptic squared | K3 + B + R + RR + pi^3 | 21 |
| SB, Tate degree 3 | K3 + T3 | 11 |
| SC, polylogarithms | K3 + T3 + Li | 21 |
| SD, vacuum sunrise | K3 + O + T3 + R | 19 |
| SE, log x elliptic | K3 + R + LR + pi^3 | 26 |
| SF, union (optional, separate tier) | about 60 distinct elements | about 61 |

The basis is guessed, not exhaustive. Sec. 5 explains why that matters.

**Precision.** A dense chance relation among n numbers at D digits has height about 10^(D/n). Set
the cap at 10^(floor(D/n) - 3), as the subspan script did. To admit genuine relations of the
control's scale (denominator about 4e5, numerators about 1e8) plus headroom up to 1e10-1e11, we need
D >= n x 13.
- Primary 400 digits, stability 340. For n = 26 the caps are 1e12 and 1e10; for n = 11 they are
  1e33 and 1e27.
- Any surviving relation must also reappear at a third precision, 460, and its residual must be
  confirmed by the independent adaptive route at 200 digits (the independent-route rule).
- SF needs about 800 digits.

## 4. Cost, from the existing scripts

- **Grid.** At 320 digits the shared double-exponential grid had 1709 nodes and took 2870 s, about
  1.7 s per node. The cost is dominated by 6 Bessel calls and 2 Struve-based L calls per node.
  Timed here at 320-450 digits: 0.1-0.3 s per call. Nodes grow as about D log D; the per-call cost
  as about D^1.5.
- **Per triple at 400 digits.** About 2200-2300 nodes at about 2.5 s each, so about 1.5-1.6 h. The
  extra integrands (O block, Gctrl, planted) reuse the same Bessel values, at no extra cost. The
  Bailey K/E, Li2, Li3 and zeta(3) take seconds. Two-mass triples are cheaper. The (1,1,1) control
  at 150 digits takes minutes.
- **Total.** Seven triples come to about 10 CPU hours (about 2 h wall in parallel). The PSLQ part
  is about 280 calls: 5 spans x {G, Gctrl, planted, basis-internal} x 2 precisions x 7 triples.
  At n of about 20-26 and 400 digits that is about 30-60 s each, so 2-4 CPU hours.
- **Optional tiers.** SF at about 800 digits is about 6-8 h per triple. The P3e-3 D-module test
  (below) reuses the `density_loop_D*` pipeline (Mellin-route operator, exact factorization,
  no-split search); by analogy with 2026-10-17 it takes several hours of sympy.
- **Code needed, in new drivers.** `three_electron_p3e_fast.py` hard-codes MASSES = (1,2,4) and a
  grid tail sized for a slowest mass of 1, so a parametrised copy is needed. The existing files are
  left untouched.

## 5. Risks

1. **The null is low-information.** The most plausible true class, a K3 period integrated once (or
   twice in the generic case) against algebraic kernels, has no finite closed basis. That is the
   lesson of the density_ring slice and of BFKNS. So no basis can turn the null into a refutation of
   "has a closed form", and widening yields diminishing evidence.
   The decisive test is function-level:
   - Take the P3e weights on the leaves-equal slice (t;1,1).
   - Find the minimal annihilator and its composition factors exactly, and run the no-split search.
   This is P3e-3 in the draft.
2. **Galois twist.** Omitting the sqrt(d)-twisted Bailey combinations gives false nulls.
3. **The degree heuristic may be wrong.** If any positive control (Gctrl at a new triple) needs an
   element of another degree, the hypothesis is dropped. SA-SE then rerun with the old
   inhomogeneous Tate set appended.
4. **Thresholds.** A collapse there may need L-values of weight-3/4 modular forms (as in BFKNS
   Sec. 3.7 and Paper 59's c_{4,0} and c_{6,0}) that are not in the threshold basis, so a null there
   does not refute a collapse. Of the planned triples, (2;1,1) has the CM modulus 1/sqrt2, so
   Chowla-Selberg puts its L-values inside the K3 block. (3;1,2) has the non-CM 1/sqrt3.
5. **Precision games.** Every returned relation must be identical at 340, 400 and 460 digits, with
   the planted control recovered in the same span. Otherwise it is noise (the 2026-10-14 lesson:
   dense relations at the cap).
6. **The transfer is assumed.** That the slice result for Delta^{-3} carries to the Slater weight
   Delta^{-5} is not shown. P3e-3 tests it.
7. **Outcome-blind scoring.** All bases, triples, caps and stability rules are fixed in the
   pre-registration before any value of G is computed. The basis-internal reductions (a pre-step)
   are run without G; they may change ranks but not spans.

## 6. Recommendation

Commit the pre-registration with four predictions:
- P3e-2W: the widened null at three generic triples, plus the swap at fixed K3.
- P3e-3: the slice D-module test, which carries the most information.
- P3e-4: collapse at threshold triples, the one outcome that would be a finding.
- P3e-5: the degree-homogeneity hypothesis, a cheap internal check carried by the controls.

Run P3e-3 and the threshold triples first. They decide more per CPU hour than another generic null.
