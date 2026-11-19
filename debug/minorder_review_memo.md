# Cold review of `debug/minorder_two_centre_memo.md`

Date 2026-11-18. Reviewer did not write the memo. Nothing outside `debug/minorder_review_*` was written.
Scripts: `debug/minorder_review_minors.py` (own collector, sympy.diff derivatives, true constants as sympy numbers,
exact determinants at two points, exact generic rank for the upper bound), `debug/minorder_review_quadrature.py`
(from-scratch prolate quadrature, no geovac integral code), `debug/minorder_review_guard.py` (guard plants, Sugiura,
DLMF asymptotics), `debug/minorder_review_extra.py` (redundant-atom plant, scaling, (AA|BB) rates).
Data: `debug/data/minorder_review_{minors,quadrature,guard,extra}.json`, logs `minorder_review_minors*.log`.

## Verdicts

| claim | verdict |
|---|---|
| (A) Theorem A | SOUND |
| (GI) Lemma GI | SOUND WITH CAVEATS (two cosmetic gaps, no counterexample; Step 2 is true) |
| (C) equal-exponent exchange, order exactly 4 | SOUND |
| (B) hybrids, pointwise; second coincidence a_d = a_p + a_c | SOUND WITH CAVEATS (coincidence confirmed by independent quadrature and exact rank; paper sentence INCOMPLETE, not wrong) |
| (G) guard weakness in `no_lower_order` | SOUND, and UNDERSTATED: two further guard weaknesses found |
| (R) "85 of 85 minors rational x surd power" | SOUND for every case I recomputed (13); proof logic of lower/upper bounds correct |
| Inventory + "no class carried an independence proof" | SOUND WITH CAVEATS (quotes exact; the "nothing else" sentence must acknowledge the Outlook's Liouville-Rosenlicht claim, which the memo does cite in Section 1) |

## Attack table

| # | attack | outcome |
|---|---|---|
| A1 | "a minimal annihilator maps each line C(R)e^{lam R} to itself" | held. d/dR(r e^{lam R}) = (r' + lam r)e^{lam R}, so any M in C(R)[d/dR] sends r_j e^{lam_j R} to m_j e^{lam_j R}, m_j in C(R). Sum m_j e^{lam_j R} = 0 forces m_j = 0 by F2 (distinct real rates, dominant rate at +oo; a nonzero rational function is not exponentially small). The k functions are C-independent and live in ker M, dim ker M <= ord M on an interval free of singular points. Upper bound: LCLM of the k first-order operators D - (r_j'/r_j + lam_j). Correct; no minimality of M even needed. |
| A2 | are class-A closed forms exponential-polynomial with RATIONAL coefficients (no sqrt R, no ln R, no algebraic coefficient)? | held. My own collector checks every coefficient is a quotient of polynomials in R with constant coefficients (`coefficients_rational_in_R` True in all 13 cases recomputed); the memo's pipeline pushes every coefficient through `QQ.frac_field(R)`, which would raise on sqrt(R) or ln R in a coefficient, for all 85. Constants (sqrt 6 etc.) are constants, fine over C(R). |
| A3 | 1/R tails, rate 0 | held. Rate 0 is one more distinct real rate with rational coefficient; VB(aa) zA=3/2 recomputed: rows e^{-3R}, e^{0}, minor -15045/1369 at 37/59, -420/121 at 11/7. |
| A4 | Theorem A needs every listed rate's coefficient nonzero; "(AA\|BB) 3 for all a != b" | held at 5 extra unequal pairs incl. (1, 1001/1000): three rates 0, -2a, -2b each with nonzero coefficient. For ALL a != b the memo leans on Roothaan, not re-derived; that is what the memo says. |
| GI-1 | is FF = Frac(C(R)[e^{lam R}]) a differential field, constants C, y' in FF? | held. The ring is closed under d/dR (A1), its fraction field is a differential subfield of meromorphic functions on C (connected domain) so constants are C. Q-linear dependence among the rates (e.g. e^{2R} = (e^R)^2) is irrelevant: FF is defined as a subfield of meromorphic functions, not as a free polynomial ring, and F2 is the only independence used. y_kappa' = e^{kappa R}/R and (ln R)' = 1/R are in FF. Ei(kappa R) with kappa > 0 (principal value) and ln R are analytic on (0, oo); work in germs on (0, oo). |
| GI-2 | Step 1 (Ostrowski linear case) | held. Take a relation with the fewest nonzero y-coefficients among those with at least one; normalise one to 1; differentiate; the derivative is a relation with fewer nonzero y-coefficients, hence has none; so all y-coefficients are constants. Standard. |
| GI-3 | Step 2: "every nonzero element of FF is rho e^{mu R}(1 + O(e^{-delta R}))" | held, NO counterexample. For P = sum p_lam e^{lam R} with rational p_lam and real distinct rates, P = p_{lam0} e^{lam0 R}(1 + sum (p_lam/p_lam0) e^{(lam - lam0)R}); each ratio is rational, so the bracket is 1 + O(R^m e^{-delta R}). Hence P has no zeros for large R (so a denominator with infinitely many zeros has them all in a bounded set on the positive axis: e.g. Q = sin-type behaviour is impossible with real rates), and P/Q has the claimed form. 1/(e^R - e^{-R} + R) = e^{-R}(1 + O(R e^{-R})). The step would FAIL with complex rates (e^{iR} - 1 has infinitely many real zeros and no dominant term); all rates in Paper 58's atoms are real. Caveat 1: write O(R^m e^{-delta R}) or shrink delta; caveat 2: the statement is along the positive real axis only, which is all Step 3 uses. |
| GI-4 | Step 3, kappa* != 0 and = 0; Ei with positive vs negative argument | held. DLMF 6.12.2 (Ei(x) ~ e^x/x sum k!/x^k, x -> +oo) and 6.12.1 (E1); my check at x = 30, 60: remainder after optimal truncation below the smallest term (4.3e-13 < 1.3e-12; 5.7e-26 < 1.7e-25) for both signs. For kappa* < 0 the dominant term is c E1(|kappa*|R) e^{|kappa*|R}, alternating series, still divergent, argument identical. Gap: the case w = 0 (sum c_kappa y_kappa = 0 identically) is excluded by Step 2's "nonzero w" but not addressed; the same asymptotic dispatches it in one line. |
| GI-5 | E1(kR) = -Ei(-kR) redundancy | held in the memo's code (`minorder_core.atom_label` asserts distinct labels; control K3 refused). BROKE in the committed guard: see G3 below. |
| C1 | 1s equal-exponent closed form, independent route | held. geovac's closed form equals Sugiura's formula (my own transcription, mpmath) to 5e-41 at R = 0.5, 1.4, 3. |
| C2 | 4x4 determinant, own code path | held. 1s: 524742295998009665281408/1526641852514853284611875 = 0.34372... at 37/59 (memo -0.34372...; sign from row order and my E1 -> -Ei convention, magnitude identical), 398.69 at 11/7; 2p-pi: 0.0034787525784618... (memo -0.0034787...), 4.4048 at 11/7. Both exact rationals, no gamma, no log. Closure {e^{-2R}, e^{-2R}ln R, E1(2R), e^{2R}E1(4R)}, size 4 = r: d/dR-closed. |
| C3 | scaling to a = b | held. K(R; z) = z K(zR; 1) to 0 at 30 digits for 1s at z = 3/2, 2 and for 2p-sigma 1s at z = 2 (Z = (4, 2)); the z = 2 2p-sigma function passes order >= 4 independently. Minimal order is invariant under R -> zR. |
| B1 | a_d = a_p + a_c coincidence, independent quadrature | held. From-scratch prolate quadrature (potential of 2p0^2 by L = 0, 2 Laplace expansion with incomplete gamma; no geovac integral code) vs the closed form, Z = (3, 6): R = 1.4: 0.038323036074667335418 both (rel 0); R = 0.5: 0.359485944739692157 vs ...15695 (rel 1.5e-19, the quadrature's dps). Controls (3, 4), (3, 3) agree to <= 2.3e-19. |
| B2 | order drop at (3, 6), exact | held. Atoms e^{-3R}, e^{-6R}, e^{-6R} ln R, e^{-6R}Ei(3R), e^{6R}E1(9R), e^{6R}E1(12R); gamma present: Ei AND gamma + ln R together, three rates. Lower bound: 4x4 minor exact rational, -3752387138451207594488633568/61134704895495005 at 37/59 and nonzero at 11/7, gamma- and log-free. Upper bound: generic rank of f..f^(4) over Q(R)(constants as indeterminates) = 4 (proof: specialisation cannot raise rank); every 5x5 minor also zero at three points. Order 4 = rates + 1. |
| B3 | other coincidences possible? | held. Rates -a_d, -a_c, -(a_p + a_c), +a_d with a's > 0: only a_d = a_c and a_d = a_p + a_c can merge. Order 5 elsewhere is measured at sample points only (memo says so). |
| B4 | paper sentence | INCOMPLETE, not wrong: "with Ei when a_d > a_c and gamma_E + ln R when a_d = a_c" is literally true at (3, 6) (a_d > a_c, Ei present) but reads as an iff for the log, and the log also appears at a_d = a_p + a_c. Same for the Sec. "Two engine findings (i)" sentence ("at the coincidence a_d = a_c the E1 collapses ... gains the pair gamma_E + ln R"), which names one coincidence as if unique. "order = rates + 1" and "no integral carries ln R at R = 0" survive (lnR coefficient of f at 0 exact zero, memo's check). No guard currently asserts "log iff a_d = a_c" on a coincidence point (`test_function_content` H5 has a_p = 2, a_c = 1, a_d = 3/2, off it). |
| G1 | reproduce K2 against the committed `no_lower_order` | held (memo right). f = (ln6 - ln2 - ln3) R e^{-R} + e^{-2R}: `no_lower_order(f, 2)` True; true order 1. Also with an E1 atom: (ln6 - ln2 - ln3) e^{-R}E1(R) + e^{-2R}: True. sqrt(2)sqrt(3) - sqrt(6) and gamma - gamma are auto-simplified by sympy and caught. |
| G2 | could a current panel case be affected? | No committed guarded case is: their constant sets are Q-independent ({log3, log6, gamma}, {log10, log2, log7, sqrt3}, {gamma}, {log2, gamma}, H2/H5/H9 sets). BUT dependent log sets ARE reachable from the engine: the memo's own scan points fam1 Z_B = 9 ({log12, log15, log3, log6}) and Z_B = 12 ({log15, log18, log3, log6}) are Q-dependent over the primes {2, 3, 5}. Harmless there only because the true minors are free of logs. |
| G3 | redundant-atom plant through the committed guard (not in the memo) | BROKE the guard. f = e^{-R}E1(2R) + e^{-R}Ei(-2R) + e^{-3R} = e^{-3R} (checked numerically: difference 0 at 141 digits). `decompose` keeps (-1, E1, 2) and (-1, Ei, -2) as separate atoms; `no_lower_order(f, 2)` returns True: claims order >= 2 for an order-1 function. Current panel closed forms do not emit both spellings (the memo's label check across 85 cases), so no current case is affected. |
| G4 | pointwise upper bound in the committed guard (not in the memo) | BROKE the idiom. `test_one_electron_order_equals_number_of_rates` asserts "order-r exists" by `rank_at(deriv_vectors(f, r), 37/59) == r`, a rank AT A POINT. Plant: e^{-R} + e^{-2R} + (2R/3 + 103/177)e^{-3R} (true order 3) has rank 2 of its 3 columns at 37/59 and rank 3 at 11/7. A point rank only bounds the generic rank from below, so it cannot prove an upper bound. Harmless for current one-electron cases (exponentials only, closure size = rates), but the line proves nothing. The memo's own upper bounds use generic rank over Q(R) (`minorder_core` step 4, `generic_order`): correct. |
| R1 | recompute minors, 5+ cases across classes | held. 13 cases: VB(aa) 3/2, (AA\|BB) (1,2), C4, E5VB, exch 1s, exch 2p-pi, HYB (3,3), (3,4), (3,5), (3,6), H2 (+ sanity r = 5 at (3,6): minor exactly 0). Every minor at both R points is a rational, or rational x sqrt(3), sqrt(6), sqrt(15); none contains gamma or a log. Magnitudes equal the memo's to all printed digits (C4 -642319984.936..., E5VB -92244.68..., HYB(3,4) 1.66421967518217e14, H2 2.6736339e9). |
| R2 | "one R0 suffices because the exponents are fixed numbers" | held in substance, wrong reason. A nonzero minor at one point R0 (finite entries) means the minor is not identically zero in C(R): that is the whole argument, valid with or without parameters. The memo never uses a zero at a point for an upper bound. |
| I1 | quote-check 4 rows | held: row 2 (l. 388), row 5 "38 of 38" (l. 488), row 7 "order 4 ... the most it can have there" (l. 455-456), row 10 "(that no lower order exists is not checked for 3d_0)" (l. 434); Caps quote (l. 507-510) exact. |
| I2 | "Nothing else carried an independence argument" | caveat. The Outlook (l. 641-645) asserts as [SYMBOLIC (proof-by-argument)] that {1, e^{-lam R}, E1(mu R), ln R} are independent over Q(R), citing Liouville-Rosenlicht. That covers class A (exponentials alone) and plain E1/ln, but not the products e^{cR}E1(kR) / e^{cR}Ei(kR) the hybrid and exchange closed forms need; and Rosenlicht's "E1 is non-elementary" is a transcendence-of-one-function statement, not a linear-independence argument for several labels. So the memo's sentence is right about the minimality claims, but the paper is internally inconsistent: the Outlook calls a weaker independence a theorem while Caps lists the needed one as assumed. GI resolves both. |

## Independent numbers

* Quadrature vs closed form, (2p0 2p0|1s 1s_B): (3,6) R = 1.4: 0.038323036074667335418 / same; R = 0.5: rel 1.5e-19;
  (3,4) R = 1.4 rel 1.8e-19; (3,3) R = 1.4 rel 2.3e-19.
* Sugiura vs geovac, 1s exchange zeta = 1: rel 5.3e-41, 2.0e-41, 1.2e-41 at R = 1.4, 0.5, 3.
* Exact minors at 37/59 (reviewer path): exch 1s 524742295998009665281408/1526641852514853284611875;
  HYB (3,3) 1054598632128432968130304/284510835767601; HYB (3,6) -3752387138451207594488633568/61134704895495005;
  HYB (3,4) -400830899149150658265302252250623774171136 sqrt(3)/4171681737231293472783203125.
* Generic rank of r+1 columns with constants as indeterminates: HYB (3,3) 4, (3,4) 5, (3,5) 5, (3,6) 4, H2 4.

## Required caveats for any paper sentence

1. Class A (one-electron, (AA|BB), s-type-pair hybrids): "order = number of distinct rates with nonzero coefficient"
   is a theorem for every exponent choice; "for all a != b, (AA|BB) 1s has order 3" additionally needs the three
   coefficients nonzero (Roothaan's formula; not re-derived here).
2. Class C (equal-exponent exchange, seven integrals): order exactly 4 is proved, conditional on Lemma GI (a hand
   proof, now reviewed once) and on the closed forms (1s independently confirmed against Sugiura).
3. Class B (hybrids with an L >= 1 pair): proved at the computed exponent values only; "order = rates + 1 at every
   ratio" stays measured. The sentence on function content must name both coincidences: gamma_E + ln R appears at
   a_d = a_c and at a_d = a_p + a_c (where it comes with Ei), and the order drops by one at each.
4. Class E (m != 0 unequal exchange) remains unproved; do not fold it into "proved".
5. The Caps clause "gamma_E, ln 2 and ln 3 are assumed independent over Q" can be dropped for the minimality claims
   (every minor is free of them); the clause "minimality assumes the atoms are independent" can be replaced by a
   reference to the independence lemma. The Outlook's Liouville-Rosenlicht sentence should either point to the same
   lemma or be narrowed: as written it is cited for a set that lacks the e^{cR}E1(kR) products.
6. Lemma GI text must handle w = 0 and say the error term is O(R^m e^{-delta R}) along the positive axis; it relies
   on all rates being real.

## Wrong answers a regression guard should reject

* A claimed order one above or one below the truth (memo K1: exchange 5, (AA|BB) 4) - already in the memo's controls.
* Q-dependent constants read as independent: K2 above and the E1 variant must be refused (current `no_lower_order`
  accepts both). A fix: evaluate the minor at the true constants, or reduce logs to logs of primes before splitting.
* Redundant spelling E1(kR) next to Ei(-kR): f = e^{-R}E1(2R) + e^{-R}Ei(-2R) + e^{-3R} must be refused at order 2
  (current guard accepts).
* An upper bound read from a rank at a point: the plant e^{-R} + e^{-2R} + (2R/3 + 103/177)e^{-3R} must not be
  reported as order 2 (current idiom `rank_at(..., r) == r` would accept such a claim if the rate count said 2).
* "ln R iff a_d = a_c": must fail at (2p0 2p0|1s 1s_B), Z = (3, 6), and at (2p1 2p1|2p1 2p1_B), Z = (2, 6).
* Order 5 at (3, 6): the 5x5 minor is identically zero; a guard claiming order 5 there must fail.
* A closed form whose tau sum was cut early (memo's plant) and a dropped ln R atom for the equal-exponent exchange:
  the frozen exact determinant (e.g. |det| = 524742295998009665281408/1526641852514853284611875 at 37/59 for 1s)
  changes.

## Not done

Theorem A random-plant check (X4) did not run (600 s cap hit earlier in the script); Theorem A rests on the argument
in A1, which I checked by hand. I did not re-derive the l > 0 exchange closed forms beyond the scaling check and the
memo's tau-stop assertion, and did not recompute the 34 unguarded class-A R8 cases.
