# Minimal order in R of the two-centre integrals of Paper 58: what is proved

Date 2026-11-18. Derivation memo; no paper, test, doc or git state touched. Unreviewed (one author, one pass).
Scripts: `debug/minorder_core.py` (machinery), `debug/minorder_run.py` (three panels), `debug/minorder_hybrid_scan.py`,
`debug/minorder_coincidence_check.py`, `debug/minorder_recon_check.py`, `debug/minorder_controls.py`.
Data: `debug/data/minorder_{p1s,r8,xeq}.json`, `minorder_hybrid_scan.json`, `minorder_coincidence_check.json`,
`minorder_recon_{p1s,r8,xeq}.json`, `minorder_recon_r8_dps70.json`, `minorder_controls.json` (+ `.log`).
Note: `minorder_r8.json` was written in two passes; the "minor is gamma- and log-free" field is recorded for the
ten hybrids and for the 1s and exchange panels, and for the other 38 R8 cases follows from their constants being
surds only.

## 0. Contradictions and scope findings

**No stated order is contradicted.** All 85 two-centre closed forms examined (30 of the 1s panel, 48 of record R8,
7 equal-exponent exchange integrals) have exactly the order Paper 58 states.

**One scope finding (hybrid class).** Paper 58 `sec:census`: "Only an $L\ge1$ multipole of $ab$ brings in $E_1$,
with $\mathrm{Ei}$ when $a_d>a_c$ and $\gamma_E+\ln R$ when $a_d=a_c$". There is a second rate coincidence the paper
(and inventory row 3b, and record R8's P3b) does not name: **a_d = a_p + a_c**, a_p the rate of the same-centre pair.
There the closed form carries Ei **and** gamma_E + ln R together, the number of rates falls from 4 to 3 and the order
from 5 to 4. Checked on (2p0 2p0|1s 1s_B) [Z_A = 3, Z_B = 6] and (2p1 2p1|2p1 2p1_B) [Z_A = 2, Z_B = 6]; the first
agrees with the module's quadrature to 1e-13 at R = 1.4 and both are continuous in Z_B across the point (neighbours at
+-1e-3, +-1e-5). The integral still has no ln R at R = 0 (exact), and "order = rates + 1" still holds. No panel case
sits on this coincidence, so no panel statement is false; the general sentence is incomplete.

**One guard weakness.** `no_lower_order` in `tests/test_paper58_r_operators.py` treats gamma and the logs of numbers
as independent symbols (or random rationals). Control K2: f = (ln 6 - ln 2 - ln 3) R e^{-R} + e^{-2R} has true order 1;
that routine reports 2. The panel constants are not of this kind, and Section 3 below removes the assumption.

## 1. Inventory of two-centre minimal-order claims in Paper 58

Paper wording quoted. "Today's backing" = what the committed tests do. The Caps paragraph applies to rows 1-9:
"the exchange operator is a common right divisor of found annihilators, not certified minimal ... $\gamma_E$, $\ln 2$
and $\ln 3$ are assumed independent over $\mathbb{Q}$; minimality assumes the exponential, logarithmic, $E_1$ and
$\mathrm{Ei}$ atoms are independent over $\mathbb{Q}(R)$". Elsewhere (Outlook) the paper asserts that independence as
"[SYMBOLIC (proof-by-argument)] a theorem (Liouville--Rosenlicht ...)" for $\{1, e^{-\lambda R}, E_1(\mu R), \ln R\}$,
without a proof for products $e^{cR}E_1(kR)$, several k, or Ei.

| # | class, parameters | claimed order | tier in paper | today's backing of the lower bound |
|---|---|---|---|---|
| 1 | one-electron 1s (S, T, <a\|1/r_A\|b>, <a\|1/r_B\|a>), exponent panel | "the operator's order is the number of distinct decay rates of the closed form, counting the rate $0$ of a $1/R$ tail" | [SYMBOLIC] paragraph | rank of f..f^(r-1) at R = 37/59, constants as symbols; atom independence assumed (3 pairs) |
| 2 | (AA\|BB) 1s | "order $2$ at equal exponents and $3$ otherwise" | [SYMBOLIC] | same, pairs (1,1), (1,2) |
| 3 | hybrid 1s; hybrid (2p0 2p0\|1s 1s_B) | elementary, order = rates; 4 at (3,3), 5 at (3,4) (claim matrix; (3,2) not guarded) | [SYMBOLIC] | same |
| 4 | exchange 1s, equal exponents | "at equal exponents its operator has order $4$" | [SYMBOLIC] | same, zeta = 1 |
| 5 | l <= 2 panel, one-electron (30) and (AA\|BB) (8) | "their operator order is the number of distinct decay rates ($38$ of $38$)" | [SYMBOLIC + MEASURED] | guard on 4 of 38 (E5VB, E4T, F5, C4), same rank test; other 34 in the debug driver only |
| 6 | l <= 2 hybrids (10) | "the order is then the number of rates plus one" (E1 present); s-type pair elementary | [SYMBOLIC + MEASURED] | guard: H2, H9 full; H5 lower bound (slow); other 7 debug only |
| 7 | exchange, one 2p on A, equal exponents: (2p sigma 1s\|2p sigma 1s), (2p+1 1s\|2p-1 1s) | "order $4$ at equal exponents, which is the most it can have there" | [SYMBOLIC + MEASURED]; upper bound shown | reconstructed operator (guess, gcrd, "not certified minimal") |
| 8 | exchange, 3d0 1s, 3d+-2 1s, 2p sigma 2p sigma, 2p+1 2p+1, equal exponents | "order $4$ at equal exponents, again the most it can have there" | [SYMBOLIC + MEASURED] | same |
| 9 | exchange 1s, unequal: (1,3/2), (1,2); bound for every a < b | "order $9$ and degree $16$"; "That the order is exactly $9$ rests on the guessed search, which excluded lower orders up to degree $60$ **[MEASURED]**" | bound [SYMBOLIC PROOF], exactness [MEASURED] | guess scan. **Proved in `debug/exact9_independence_memo.md` (reviewed; guard `tests/test_paper58_exact_order.py`); the paper sentence is not yet updated** |
| 10 | exchange sigma = 0, l > 0, unequal: 2p sigma 1s, 2p sigma 2p sigma, 2p+1 2p+1 at (1,2); 3d0 1s | "the order is $9$"; 3d0: "(that no lower order exists is not checked for $3d_0$)" | [SYMBOLIC + MEASURED] | guess scan; proved for every a < b in the same exact9 memo (mixed cases: the l > 0 orbital has the smaller exponent) |
| 11 | exchange m != 0, unequal: (2p+1 1s\|2p-1 1s), (3d+2 1s\|3d-2 1s) at (1,2) | "At $(a,b)=(1,2)$ the order is $9$" | [SYMBOLIC + MEASURED] | guess scan only; "The $2p_{\pm1}1s$ and $3d_{\pm2}1s$ integrals are not covered" by the nine-function closed form |

Already proved before this memo: rows 9, 10. Nothing else carried an independence argument: "order = number of
rates" was a rank computation under an assumed independence.

## 2. The two lemmas

Coefficient field C(R) throughout (an ODE with polynomial coefficients is a C(R)-relation among derivatives).

(F2) Exponentials e^{lam R}, distinct real lam, are linearly independent over C(R) (dominant rate at R -> +oo; a
nonzero rational function is not exponentially small). As in the exact9 memo.

**Theorem A (exponential-polynomial functions).** Let f = sum_{j=1..k} r_j(R) e^{lam_j R}, lam_j distinct real,
r_j in C(R) nonzero. Then the minimal order of a linear ODE for f over C(R) is exactly k.

Proof. Upper bound: each f_j = r_j e^{lam_j R} is killed by d/dR - f_j'/f_j, and the least common left multiple of k
first-order operators has order <= k. Lower bound: let M f = 0, M in C(R)[d/dR]. M f_j = m_j e^{lam_j R} with
m_j in C(R), because d/dR maps C(R) e^{lam R} to itself. So sum m_j e^{lam_j R} = 0 and by (F2) every m_j = 0: M
kills each f_j. The f_j are k functions linearly independent over C (F2 again) in the solution space of M, which has
dimension ord M. QED. No determinant, no parameter restriction, no exceptional ratio: when rates coincide they
are counted once and the order drops by exactly the number of merged rates.

**Lemma GI (general independence).** Let K be a finite set of distinct real numbers; y_0 = ln R and
y_kappa = Ei(kappa R) for kappa != 0 (so E1(kR) = -Ei(-kR); principal value for kappa > 0). The atoms
e^{cR}, and e^{cR} y_kappa, over distinct real c and distinct (c, kappa), are linearly independent over C(R).

Proof. Let E = C(R)[e^{lam R}, lam real] and FF = Frac(E), a differential field of meromorphic functions on C with
constant field C. A relation among atoms, grouped by y, reads a_0 + sum_kappa a_kappa y_kappa = 0 with a's in E. It
suffices to show 1 and the y_kappa are independent over FF: then each a is 0 in E and (F2) gives each rational
coefficient 0.

Step 1 (Ostrowski's reduction, the elementary linear case). y_kappa' = e^{kappa R}/R in FF. Take a relation with at
least one nonzero a_kappa and the fewest such; divide so one of them is 1; differentiate. The result is a relation
with fewer nonzero y-coefficients (the a_kappa'), so they all vanish: every a_kappa is a constant c_kappa, and
w := sum c_kappa y_kappa = -a_0 lies in FF, not all c_kappa zero.

Step 2 (elements of FF at R -> +oo). An exponential polynomial P = sum p_lam e^{lam R} equals
p e^{lam_0 R}(1 + O(e^{-delta R})), lam_0 its largest rate; so every nonzero w in FF is
rho(R) e^{mu R}(1 + O(e^{-delta R})) with rho in C(R) nonzero, mu real, delta > 0.

Step 3. Let kappa* be the largest kappa with c_kappa != 0. Ei(x) ~ (e^x/x) sum_k k!/x^k as x -> +-oo (DLMF 6.12;
number from memory), a divergent series.
* kappa* != 0: w e^{-kappa* R} = c Ei(kappa* R) e^{-kappa* R} + (exponentially smaller terms) has the asymptotic
  power series c sum k!/(kappa* R)^{k+1}. By Step 2, mu = kappa* and that series is the expansion of rho at infinity,
  which converges. A divergent series is not a convergent one (uniqueness of asymptotic power series). Contradiction.
* kappa* = 0: w = c_0 ln R + (exponentially small). By Step 2 a rational function would equal c_0 ln R + o(1).
  Impossible for c_0 != 0. QED.

GI does not use the local form at R = 0. It covers every atom family met below, at any real rates and arguments,
with no exceptional ratio (labels must be distinct; `minorder_core` asserts it and refuses a redundant pair, control
K3). It does not cover the master integral f0 of the unequal exchange; that needs the exact9 memo's argument.

**Consequence.** If f = sum q_a(R) atom_a with q_a in C(R), and the atom set is closed under d/dR, then the order
of f is the rank over C(R) of the coordinate matrix of f, f', f'', ... on the atoms.

## 3. The exact computation (leg D), with true constants

`minorder_core.analyse(f, r)`:
1. atoms and the derivative calculus from `tests/test_paper58_r_operators.py` (`decompose`, `to_field`, `d_atoms`;
   loaded, not copied). Recon check: the decomposition and its first two derivatives reproduce f, f', f'' at
   R = 7/5 to <= 1e-30 relative at dps 40 in 84 of 85 cases; H10 reads 1.0e-30, 2.0e-30, 7.0e-29 there
   (cancellation among large coefficients) and 3e-61, 1e-60, 9e-58 at dps 70, H6 likewise 1e-63 to 2e-62.
2. closure of the atom set under d/dR, size n; order <= n.
3. the TRUE coordinate matrix at R0 = 37/59: entries are polynomials over Q in the constants (gamma, logs,
   normalisation surds). An r x r minor is computed exactly as such a polynomial and then evaluated at the true
   constants in interval arithmetic (60 digits). Nonzero => rank r over C(R) => no ODE of order < r. One R0 is enough
   here because the exponents are fixed numbers (there is no parameter to be unlucky in).
4. where n > r: order <= r from the generic rank over Q(R) of the constant-split vectors of f..f^(r) (a relation
   among split vectors is a relation among the true functions).

**Result of step 3, stronger than needed:** in all 85 cases the minor is a single monomial: a nonzero rational
times a power of one normalisation surd, with **no gamma and no logarithm in it**. So the lower bounds do not depend
on any arithmetic property of gamma, ln 2, ln 3, ...: the paper's "assumed independent over Q" caveat is not needed
for minimality.

Controls (`minorder_controls.py`): claimed order one too high is refused (exchange 5, (AA|BB) 4); a function whose
constants are dependent in truth (K2 above) is refused at order 2 and accepted at order 1; a redundant atom pair is
refused.

## 4. Verdicts per class

| class | verdict | range | exceptional set |
|---|---|---|---|
| A. one-electron and (AA\|BB), 1s and l <= 2 (rows 1, 2, 5), and elementary hybrids (s-type pair; row 3 first part, H9) | **PROVED** (Theorem A; determinant agrees, 63/63: 24 of the 1s panel, 39 of R8) | every exponent choice for which the closed form is exponential-polynomial: order = number of distinct rates with a nonzero coefficient | none; coincidences merge rates (e.g. (AA\|BB) 1s: 3 -> 2 at a = b) |
| B. hybrids carrying E1/Ei/ln (rows 3, 6) | **PROVED at the computed exponents** (GI + exact minor + exact upper rank): 12 panel cases (HYB 33/34/32; H1-H8, H10) and 14 further scan points | fixed exponent values only; no formula in the ratio | order rates + 1 at every point computed; drops 5 -> 4 at a_d = a_c **and at a_d = a_p + a_c** (Section 0). Generic statement **NOT PROVED** |
| C. exchange, equal exponents (rows 4, 7, 8): 1s, 2p sigma 1s, 2p+-1 1s, 3d0 1s, 3d+-2 1s, 2p sigma 2p sigma, 2p+1 2p+1 | **PROVED**, order exactly 4 | a = b = 1 computed (1s also at zeta = 3/2, 2); every a = b by the scaling K(R; zeta) = zeta K(zeta R; 1) | none within a = b |
| D. exchange, unequal, sigma = 0 (rows 9, 10) | PROVED earlier (exact9 memo), every 0 < a < b | as stated there | a = b (order 4) |
| E. exchange, unequal, m != 0 (row 11) | **NOT PROVED** | - | gap: no closed form on a finite d/dR-closed list is established for these two integrals |

Details.

**A.** Function list: e^{lam R}, lam in the rate set (0 for a 1/R tail; -a, -b; -2a, -2b; sums for l > 0). Independence:
F2. Rank: not needed (Theorem A); the determinant of step 3 is the Wronskian-type check and is nonzero in all 63.
What Theorem A needs from the closed form is that each listed rate has a nonzero rational coefficient; step 3's
nonzero minor certifies that too (a vanishing coefficient would make the n x n determinant zero). For 1s (AA|BB)
"3 otherwise" is proved at the three unequal pairs computed; for all a != b it follows from Theorem A once the three
coefficients are nonzero, which Roothaan's formula shows (not re-derived here).

**B.** Function lists (closure, 6 or 7 atoms): e.g. Z = (3,2): e^{-2R}, e^{-3R}, e^{-6R}, e^{-2R}E1(R), e^{-2R}E1(4R),
e^{2R}E1(5R), e^{2R}E1(8R); in general rates -a_d, -a_c, -(a_p + a_c), and E-atoms e^{-a_d R}Ei((a_d - a_c)R),
e^{-a_d R}Ei((a_d - a_p - a_c)R), e^{a_d R}E1((a_p + a_d)R), e^{a_d R}E1((a_p + a_c + a_d)R), a zero argument
turning into gamma + ln R. Independence: GI. The integral spans only a 4- or 5-dimensional submodule of the 6- or
7-dimensional closure, so the upper bound is step 4 (generic rank r of r + 1 columns, exact; equals the frozen
operators where they exist), and the lower bound is an r x r minor. Scan (`minorder_hybrid_scan.py`), order / rates:

fam1 (2p0 2p0|1s 1s_B), Z_A = 3 (a_p = a_c = 3), Z_B = 1, 3/2, 2, 5/2: 5/4; 3: 4/3; 4, 5: 5/4; 6: 4/3; 7, 9, 12: 5/4.
fam2 (2p1 2p1|2p1 2p1_B), Z_A = 2 (a_p = 2, a_c = 1), a_d = 1/2: 5/4; 1: 4/3; 3/2, 2: 5/4; 3: 4/3; 4, 5: 5/4.

All three orderings of a_d against a_c and a_p + a_c are sampled, and both coincidences. The four rates are
-a_d, -a_c, -(a_p + a_c), +a_d, so these two are the only possible coincidences among them; that the order is 5 at
every other ratio is measured at these points, not proved.

**C.** Function list (same for all seven integrals, a = b = 1): e^{-2R}, e^{-2R} ln R, E1(2R), e^{2R}E1(4R), closed
under d/dR (n = 4 = r, so the upper bound is immediate and agrees with the paper's "the most it can have").
Independence: GI with labels (c, kappa) = (-2, none), (-2, 0), (0, -2), (2, -4). Rank: the 4 x 4 determinant at
R0 = 37/59 is an exact nonzero rational, free of gamma (gamma enters the e^{-2R} row only as gamma times the ln R
row): 1s -0.3437..., 2p sigma -1.190..., 2p pi -3.479e-3, 3d sigma 13.14..., 3d delta -1.665e-4, pp sigma
-3733.07..., pp pi0 0.4858... . The l > 0 closed forms are `geovac.qfd_core.exchange_closed_form` with the tau sum
checked to stop. This replaces "common right divisor ... not certified minimal" for the equal-exponent operators:
order 4 is now two-sided.

## 5. What a regression guard should reject (per proved class)

A. (i) a claimed order different from the number of distinct rates with nonzero coefficient, in either direction;
(ii) a closed form with a planted non-exponential atom passing as class A; (iii) a rate whose coefficient vanishes
identically being counted. Test the n x n determinant with true constants, and separately that the atom kinds are
"1" only. Parametric plant: two rates made equal must lower the order by one.

B. (i) order r - 1 accepted (plant: drop the minor to r - 1 columns must still pass, r columns of a function with one
E-atom coefficient zeroed must fail); (ii) a minor that depends on gamma or a log being evaluated with symbols as
independent (plant K2); (iii) the upper rank computed at a point instead of over Q(R); (iv) the second coincidence:
a guard for "ln R iff a_d = a_c" must fail at Z = (3, 6); (v) redundant atom labels (E1(kR) next to Ei(-kR)).

C. (i) order 5 accepted or order 3 accepted for any of the seven integrals; (ii) a closed form whose tau sum was cut
before it stops (plant tau_max = 2 for 3d: the function changes, the test must notice through the frozen value of
the determinant or an independent value check); (iii) the ln R atom dropped (would make E1 and the exponential a
three-function span); (iv) gamma treated as a free symbol in a way that changes the determinant (it must be
gamma-free: pin "no gamma in the minor").

GI and Theorem A themselves are arguments, not tests; a guard can pin their inputs only: the derivative rule of each
atom (already `test_atom_derivative_matches_sympy`), distinct labels, and the recon check of Section 3.

## 6. Where this could be wrong

1. **GI and Theorem A are hand proofs, unreviewed.** GI's Step 3 uses the asymptotic series of Ei at +-oo (textbook;
   DLMF number unverified) and uniqueness of asymptotic power series. Step 1 needs the coefficient field FF to be a
   differential field with constants C containing the y' (it is). If GI failed for some label set, every lower bound
   of classes B and C would be conditional again; class A does not use GI.
2. **Everything is conditional on the closed forms** (`qfd_core`, `two_center_eri`, `two_center_one_electron`), which
   I did not re-derive; their guards are the existing tests. The coincidence closed form at a_d = a_p + a_c is not
   covered by an existing test; I checked it against quadrature (one family, one R) and continuity only.
3. **`decompose` and `d_atoms` are trusted** (guarded in the r-operators test; recon check here to second
   derivative only, at one R).
4. **Interval evaluation**: constants are bracketed from sympy's 70-digit values, not from a certified enclosure.
   Since every minor turned out to be rational times a surd power, nothing depends on this.
5. **Class B is pointwise.** "Order = rates + 1" for all exponent ratios is not proved; a ratio at which the
   integral falls into a smaller submodule without a rate coincidence is not excluded.
6. **Class C scaling** to a = b != 1 for l > 0 is dimensional (hydrogenic orbitals, all charges scaled), not
   recomputed; 1s is recomputed at zeta = 3/2 and 2.
7. **Scope of "order"**: homogeneous linear ODE in R with rational coefficients at fixed exponents.
8. The R8 panel's 34 unguarded class-A cases and 7 unguarded hybrids were rebuilt here from the live engines with the
   record's panel definitions; the stored orders (`debug/data/row8_panel.json`, `two_centre_r_operators.json`) were
   the claims tested.
