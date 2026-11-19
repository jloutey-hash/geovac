# Cold review of `debug/tree_attain_memo.md`

Date 2026-11-18. Reviewer did not write the memo. Scripts `debug/tree_attain_review_main.py`,
`debug/tree_attain_review_modp.py`, `debug/tree_attain_review_p3.py`; data
`debug/data/tree_attain_review_{main,modp,p3}.json` (+ `.log`). No paper, test or tracked file edited.

## Verdicts

**Claim A: GAP REAL.** The proposition's statement is true (by the memo's Lemma 1 + the paper's own split); its
written proof is incomplete for n = 4, 6, 8 and wrong as worded for every n.

**Claim B: SOUND WITH CAVEATS.** I found no broken step. None of the caveats blocks a paper statement for s-type
trees, 3 <= n <= 8; they fix its wording and tier (list below).

## Claim A in detail

The failing sentences (proof of `prop:star_attain`): *"In the filtration of Proposition tree_rung the pieces with
c <= n-2 and the pieces of order one are nonsingular at t = -(n-1); quotients and extensions of such modules are
connections there, whose solutions are holomorphic on every branch. So a piece with c = n-1 must survive."*

The pieces of `prop:tree_rung` are modules over D_t = C(t)<theta_t>, and the proposition is about composition
factors over C(t). Over C(t) the inference fails twice:

1. Membership does not control poles. 1/(t-p) lies in C(t).1, a module whose only piece is the trivial one.
   A function can be non-holomorphic at p inside a module all of whose pieces are regular at p.
2. Extensions do not stay regular. Over the field of germs at p, Ext^1(trivial, trivial) = C dt/(t-p), not 0.
   The memo's toy is a correct instance: d/dt o (t-p) d/dt has the solutions 1 and log(t-p) [checked, P0], and
   its two factors d/dt and (t-p) d/dt are both equivalent to d/dt over C(t).

Reading the sentence over operators regular at p (where "extensions of connections are connections" is true)
does not rescue it: it then yields a top piece that is nonzero as a module over that local ring, which can be
torsion and vanish over C(t) (in the toy 1/(t-p): the quotient O_p[1/(t-p)]/O_p).

Is the toy a counterexample *as the paper uses the inference*? Yes. The paper's only input is "the fibre is not
holomorphic at p"; for n = 4 the non-analytic term is a pole (eps^-1), for n = 6 a logarithm, for n = 8
eps log eps. Poles come from (1), logarithms from (2). Nothing in the written proof excludes either.

What I looked for that the memo might have missed, and did not find: the paper does not use the singular
behaviour of a specific solution of a specific piece, does not bound apparent singularities, and does not state
a monodromy argument in the text. (The guard's docstring, `tests/test_paper63_star_attainment.py`, does say
"the monodromy argument" for the half-integer exponent; the paper's proof text does not.)

By n:
- n = 3, 5, 7: conclusion rescued in one line, not on the page. If every factor had only meromorphic local
  solutions at p, the fibre would lie in C((eps))[log eps]; its jump has leading term eps^(half-integer).
  So some factor has non-trivial local monodromy, and only the c = n-1 piece can.
- n = 4: not covered by the argument; covered by the paper's computed factorisation of the star of four, which is
  SYMBOLIC + MEASURED (the order-eight factor annihilates to 3e-45; minimality bounded mod a prime, not proved).
- n = 6, 8: not covered by anything written in the paper.

Statement vs proof: the statement holds for all 3 <= n <= 8 by Claim B applied to the star. This is a proof gap,
not a false proposition.

## Claim B: attack table

| # | attack | result | what I did |
|:-:|:--|:--|:--|
| 1 | weighted filtration, commutative graded, lift to elements | tried, held | [theta, z] = 2z has weight 2 < 3, so gr A = C[Z, Theta], Noetherian; F_w is finite-dimensional and zero for w < 0. theta F_w in F_(w+1) checked on e_i by hand. The lift is the standard statement Ch(A g) inside Ch(P) for the induced (good) filtration on a submodule, then the Nullstellensatz; the memo's "G_w in F_(w+w1) A g" is the usual comparison of two good filtrations and is correct for a filtration bounded below. The annihilator is of the ELEMENT g (a_g = X^m - b), not of the graded piece. Sums over S are handled by the direct sum. |
| 1b | conclusion on trees the memo could not finish | tried, held | own implementation of the module, mod p = 2147483629: fork of five (roots 0, 2, 4), path of five from the centre, path of six from an end. Symbol splits completely, roots integers |t| <= n-2, nonzero at -(n-1), a H pure and nonzero. Four memo rows reproduced as controls (same order, weight, root multiplicities). Not finished: fork rooted at the branch vertex (z-degree above my cap), broom of six (time cap). |
| 2 | (*) order = weight, leading coefficient t^-w sigma(1, -t) | tried, held | exact, 12 random b of weight 3..9; phi respects [theta, z] = 2z. |
| 3 | "a H = q K0^(n-1)", cyclic where; could phi(a) T(H) = 0 | tried, held | Cyclic in Sym^(n-1) B over C[z] (integer pivots m!/(m-a)!, `prop:tree_rung`). T is on [0, oo); boundary terms vanish (A.H is O(x^2 log^m) at 0; Sym^(n-1) B is O(log^(n-1)), times x). phi(a) T(H) = 0 is excluded only through (S); that is the right reason and it is the whole proof. EXACT: for the stars of three and four, a applied to chi^(n-1) x^(3-n) in the ring Q[x, 1/x][chi, K0, K1] contains no chi and equals the memo's g term by term (14 and 99 terms). |
| 4 | Lemma 2: lower pieces, [0,1], boundary at 1, "unbounded" | tried, held | Lower pieces obey phi(a_S) y = (Laurent polynomial) x theta^j K0(t), singular only at 0 and integers |t| <= n-2; the upper half plane avoids them, and the disc at p is free of them. [0,1] part: K0(tx) is analytic on the cover of C minus 0. Unboundedness is enough: a function holomorphic at p has bounded derivatives of every order along eps -> 0+, and the top term's J (n <= 6) or J' (n = 7, 8) is not, the rest being analytic at eps = 0. No statement about the type of singularity is needed, so the logarithmic cases are covered. |
| 5 | free-generator model vs true functions | tried, held | The model maps A-linearly onto the functions; a H lands in the pure span whatever relations hold, and its being nonzero is never assumed (it is the conclusion of the contradiction). Exact confirmation for two stars in row 3. |
| 6 | root dependence | tried, held | The proof does not use the root. Checked: top constant for the fork of five at four different roots; exact symbols for three rootings of the fork (two of them the same rooted tree) and both rootings of the path of five. |
| 7 | simplicity of D F_n, n <= 8 | inherited, guard passes | Paper: n = 2..5 argued in `sec:four_density` (Fuchs count); n = 6, 7, 8 only cited there and in `prop:tree_bound`, with the argument in the guard. I ran `test_paper63_four_density_rung.py -k "irreducible or thresholds"`: 4 passed. Tier SYMBOLIC by guard; I read the exponent argument (a right factor takes a sub-multiset of exponents) and accept it, not re-derived. |
| 8 | recompute constants and the n = 4 identity | tried, held | below |

## Independent numbers

**Top constant C** (own nested Gauss-Laguerre in exponentially scaled variables, x = 100 .. 1600, Richardson; no
code shared with the memo):

| tree, root | n | prod k_u^2 | predicted C | computed | rel. diff |
|:--|:-:|:-:|:--|:--|:-:|
| star 5 | 5 | 1 | 2.4674011 | 2.4674011 | 2e-12 |
| fork 5 (0-1, 1-2, 1-3, 3-4), root 0 | 5 | 64 | 0.03855314219 | 0.03855314219 | 7e-12 |
| fork 5, root 1 (branch vertex) | 5 | 4 | 0.6168502751 | 0.6168502751 | 3e-12 |
| fork 5, root 3 | 5 | 9 | 0.2741556778 | 0.2741556778 | 3e-12 |
| fork 5, root 4 | 5 | 144 | 0.01713472986 | 0.01713472986 | 8e-12 |
| path 5, centre | 5 | 16 | 0.1542125688 | 0.1542125688 | 6e-12 |
| spider 6 (0-1, 1-2, 1-3, 3-4, 3-5), root 0 | 6 | **225** | -0.01374412747 | -0.01374412747 | 6e-12 |
| branched 7 | 7 | 2304 | 0.001682198171 | 0.001682198171 | 1e-11 |
| branched 8 | 8 | 2304 | -0.002108322749 | -0.002108322749 | 8e-12 |

**Memo defect (typo, not a logic error):** its table lists prod k_u^2 = 64 for the spider of six. The sizes are
1, 1, 1, 3, 5, so the product is 225. The memo's predicted number 0.00245212 is the 225 value
(|C|/sqrt(10 pi)), so only the column entry is wrong.

**phi(a) T(H) = T(aH)**, phi under the integral sign, mpmath 40 digits, own fixed Gauss-Legendre panels, chi by
cumulative quadrature (checked against the Struve closed form to 7e-41):

| tree | t | phi(a) T(H) | T(aH) | rel. diff | largest node term |
|:--|:-:|:--|:--|:-:|:-:|
| star 3 | 3/2 | -93648.159753451791001834 | -93648.159753451791001797 | 4e-22 | 3.9e5 |
| star 3 | 5/2 | -1340.9060458453934938649 | -1340.9060458453934938285 | 3e-20 | 9.4e3 |
| star 4 | 3/2 | -6.285267392144423338e54 | -6.285267392144422217e54 | 1.8e-16 | 1.8e70 |
| star 4 | 5/2 | -2.4512823513869041942086e49 | -2.4512823513869041940986e49 | 4.5e-20 | 1.0e60 |

The star-4 values cancel 15 and 10 digits against the largest node term; the memo's float64 attempt could not
see them (it printed 2.56e55 against -2.45e49 at t = 5/2). T(aH) at t = 5/2 agrees with the memo's -2.4512823514e49.
The control I registered in the script (a with its lowest term dropped) is VOID: that term contributes below
1e-30 of the result, so the row in the JSON repeats the true value and tests nothing. The failing case on record
is my own slip:

my first version applied
t^(-2k) theta^(2k) for phi(z^k) instead of (t^-2 theta^2)^k = t^(-2k) prod_(i<k) (theta - 2i)^2, and the n = 3
identity failed at order one (21607.7 against -93648.2). The memo's kernel is the correct one.

## Caveats required in any paper sentence

1. **Scope.** s-type densities, Coulomb kernels, fingerprint weights, every non-root scale equal to one, t on one
   vertex. l > 0 (the rest of `conj:tree_attain`) is not touched. Unequal non-root scales are not covered.
2. **What is proved.** F_n lies in the module of the fibre, so the minimal operator has a composition factor
   equivalent to F_n's. Nothing about minimality, the other factors, or the integrated value.
3. **n <= 8 is inherited** from the irreducibility of F_n's operator; for n >= 9 only "meets D F_n".
4. **Tier.** The argument is a proof (SYMBOLIC PROOF), with Lemma 1 a standard characteristic-variety statement
   written out by the memo's author and read by one reviewer. Exact confirmation: integers for six rooted trees
   (memo), modulo one prime for five more (this review), including a branched tree and n = 6. No exact check at
   n = 7, 8 and none for the fork rooted at its branch vertex. Say so.
5. **The star proposition's proof must be replaced, not cited.** A paper sentence "the proof of
   prop:star_attain runs as written for trees" would import the gap of Claim A. The conjecture's text currently
   says exactly that ("would run as written") and has to change with it.
6. **Odd-n remark.** If the paper keeps the monodromy argument for n = 3, 5, 7 it must be stated; it is not a
   substitute for Lemma 1 at even n.
7. Fix the spider-of-six product (225) before any number from the memo's table is copied.

## Wrong answers a regression guard should reject

- "All composition factors regular at p implies holomorphic at p": the guard should carry the toy
  (d/dt o (t-p) d/dt, solution log(t-p)) and 1/(t-p) as cases the *old* inference would wrongly accept.
- An operator a whose symbol at (1, -t) vanishes at t = -(n-1) (e.g. a multiplied on the left by L_(n-1)): the
  guard must compute the symbol's roots and fail on a root of absolute value n-1.
- a H not pure: one non-pure coordinate left (drop one term of a).
- The wrong phi(z^k): t^(-2k) theta^(2k) in place of (t^-2 theta^2)^k (fails at order one, above).
- The top constant with k_u = 1 at every vertex (the star's value) on a non-star tree: off by prod k_u^2
  (64 for the fork from its end).
- The top constant with prod k_u in place of prod k_u^2, and with the wrong sign (-1)^n.
- A lower piece credited with rate n-1 (it would be singular at p): assert every S other than the full set has
  |S| <= n-2.
- For n = 7, 8: a check that tests J rather than its eps-derivative would see a bounded quantity and pass a
  vanishing top term; the guard must test the derivative there.

## Not done

- Lemma 1 exactly over Q for a branched tree (only modulo a prime). A wrong reduction would need the least order
  to drop modulo that prime; not excluded by a second prime.
- The n = 4 identity for a non-star tree (nested potentials in high precision).
- Any independent derivation of the irreducibility of F_6, F_7, F_8.
