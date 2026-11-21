# p63_pi_nonzero: does prop:tree_attain_hom extend to non-homogeneous P_v with Pi != 0?

Date 2026-11-21. Scripts `debug/p63_pi_nonzero_step1.py` (step (i), mod p, the guard model `tests/_tree_attain_l.py`
used unchanged, with four mixed-degree densities added at runtime) and `debug/p63_pi_nonzero_top.py` (step (ii)
top constant by its own exterior-tail quadrature, no model). Data `debug/data/p63_pi_nonzero_step1.json`,
`debug/data/p63_pi_nonzero_top.txt`. Replacement text `debug/p63_pi_nonzero_paragraph.tex`. No paper, test,
CLAUDE.md, CHANGELOG or git state touched.

## Verdict: GO

The proof of `prop:tree_attain_hom` never uses homogeneity except in one identity, Pi = (-(n-1))^{d_0} x
(principal-symbol average). That identity is what turns the proof's actual hypothesis, Pi != 0, into the
conjecture's. With Pi != 0 as the hypothesis, every step runs unchanged for P_v of mixed degree <= 2. The paper's
own "[ARGUMENT] ... not proved" (l.1177) is therefore stale. It is in tension with l.1169-1171 ("the proof needs
Pi != 0 instead"). The resolution is that the proof already covers Pi != 0.

## Step by step

**(0) Relocation.** This is term by term in P_0. Each term has degree <= 2, so the surface terms are the same as in
the homogeneous case: never two derivatives on K0, and the Laplacian term's cutoff difference is
4 pi t^3 r0 Psi(0) -> 0. Unchanged.

**(i) Module and purity.**
- *Symbol.* The filtration and good-filtration argument never look at degrees. Twists r^p, channel potentials
  and Gaunt constants all have weight 0. Unchanged.
- *Parity, the question asked.* The grading p + sum L_kept - L_channel (mod 2) is not about the degree d_v. A
  term h_l(d) lap^j has twist -l-2j in channel l, so its grade is -2j, which is even. This holds for each term
  separately, whatever other terms share the density. Mixing degree 0 and 1 (1 + d_z) gives two terms, (twist
  0, channel 0) and (twist -1, channel 1), and both have grade 0. Gaunt, Poisson, theta, z and the root's d_z
  (one step in twist and in channel) and lap preserve the grade. G therefore has grade 0, the pure top carries
  even powers of x, and z^J a G = q K0^{n-1} with q in A.
- *a G.* It is pure and even by the above. Nonzero as a function comes from step (iii), as in the s-type case.
- The existing paper sentence (l.1147-1150) is already stated as a grading and needs no change.

**(ii) Top term.**
- *Exponent.* Unchanged at (5-n)/2. Every term of P_v(d)K0 has the tail of K0 at leading order:
  d_z K0 = -mu K1, d_z^2 K0 = mu^2 K0 + O(K1/r), lap K0 = K0 - K1/r. So the leading density is
  P_v(-omega) K0 (1 + O(1/x)), and all degrees sit at the same order x^{-1/2} e^{-x}. This is exactly why the
  constant is the full polynomial Pi and not the principal part.
- *Constant.* Exactly C Pi, with C the s-type radial constant. The ingredients:
  - Channel by channel, Watson's lemma gives a leading tail -s_0/k^2 independent of L (below).
  - There are finitely many channels, so the leading angular factor of each potential is the product of the
    P_v(-omega) in its subtree times radial constants, pointwise and uniformly in omega. This holds even where a
    channel's leading coefficient vanishes, since the statement is only "+ O(x^{a-1} e^{-kx})".
  - At the root, P_0(-d) acting on Psi^top ~ f(omega) x^b e^{-(n-1)x} gives P_0((n-1)omega) f at leading order:
    each -d_z gives (n-1)mu, and angular derivatives cost 1/x.
  - G^top = x^2 <P_0((n-1)omega) prod P_v(-omega)> x (radial) = C Pi x^{(5-n)/2} e^{-(n-1)x}(1 + O(1/x)).
  - For homogeneous P_v, P_0((n-1)omega) = (-(n-1))^{d_0} sigma_0(-omega), which gives the current formula.
- *Rigour of the local rule.* The [ARGUMENT] paragraph derives it as the formal solution of
  (theta-L)(theta+L+1)Phi = -4 pi r^2 s. On its own that is a heuristic, because one still has to identify the
  formal series with the true tail. The rigorous route takes three lines:
  - The exterior tail is E_L = (1/(2L+1))[-r^{-L-1} int_r^oo s^{L+2} W + r^L int_r^oo s^{1-L} W].
  - Watson's lemma on both integrals, with W = s_0 x^a e^{-kx}(1 + s_1/x + O(x^{-2})), gives leading terms
    s_0 r^{a+1}/k that cancel. The s_1 terms also cancel.
  - The survivor is [(1-L+a) - (L+2+a)] s_0 r^a/k^2 / (2L+1) = -s_0 r^a/k^2, for every L.
  - The step needs W to have a two-term asymptotic expansion. This holds by induction from the leaves: every
    function is a finite sum of products of K_nu's and of tails defined by such integrals, all with Poincare
    expansions in powers of 1/x times x^a e^{-kx}, and no logs at infinity.
  - The root's differentiation of the asymptotics is legitimate because E' is given by the same two integrals,
    and E'' by the channel ODE.

  This is the **weakest step**, but only in wording. prop:star_attain already uses it as a proof step, and its
  "(1 + O(1/x))" there understates what is needed: with only that, the leading terms cancel and the O(1/x)
  remainder is as large as the answer. The replacement text says "Watson's lemma to second order". Termwise
  integration by parts, the ODE and the true potential to 40 digits are already guarded in
  `tests/test_paper63_local_rule.py`.

**(iii) Unchanged.** If phi(q)F_n = 0, T(G) would solve phi(z^J a)y = 0, which is regular at -(n-1), against
(ii). The simplicity of D F_n for n <= 8 does not involve P_v.

**n >= 9.** This is also unchanged. (ii) needs only C Pi != 0 and the leading non-analytic term of
int x^{(4-n)/2} e^{-eps x}: a half-integer power, or eps^m log eps, so some eps-derivative is unbounded. The
O(1/x) corrections give higher non-analytic powers. (iii) then gives a nonzero common submodule without
simplicity.

## Numerical checks run

**Step (i)**, least-order a with a G pure, computed modulo 2147483629 in the guard's model. Each row satisfies
all of the following:
- the symbol at (1,-t) splits into integer roots |r| <= n-2;
- the symbol is nonzero at t = -(n-1);
- a G has 0 non-pure coordinates and some nonzero pure ones;
- every pure r-power is even;
- the Z/2 grading is 0 on G, theta G, ..., theta^4 G.

| tree (root first) | n | Pi | principal avg | order | z-deg | weight | roots (mult) |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:--|
| N1: s root; leaves 1+d_z, 1+d_z (star) | 3 | 4/3 | 1/3 | 13 | 18 | 41 | 0^5, +-1^4 |
| N2: 1+d_z root; leaves 1+d_z, s (star) | 3 | 1/3 | 1/3 | 16 | 21 | 46 | 0^4, +-1^6 |
| N3: 1+d_z^2 root; s, s (star) | 3 | 7/3 | 1/3 | 7 | 8 | 19 | 0^3, +-1^2 |
| N4: s root - (1+d_z) - (1+d_z) (path) | 3 | 4/3 | 1/3 | 10 | 8 | 18 | 0^2, +-1^4 |
| N6: 1+d_z+d_z^2 root - (d_z+lap/2) - p_z (path) | 3 | 4/5 | 0 | 6 | 6 | 14 | 0^2, +-1^2 |
| N5: path 4, 1+d_z root - (1+d_z) - p_z - s | 4 | -2/3 | 0 | 10 | 9 | 22 | 0^4, +-1, +-2^2 |

N1, N2, N4 and N5 mix degrees 0 and 1 inside one density (the parity question). N3 and N6 mix even degrees, and
N6 mixes 0, 1 and 2 at the root.

**Step (ii)**, constant of G^top. The check uses its own nested scaled Gauss-Laguerre quadrature (90 nodes) and
the channel-L tails, with the root's derivatives done by hand through the d_z ladder and the channel ODE.
Richardson extrapolation is on x = 50..400 with the exponent 1 = (5-n)/2 fixed, and the local exponent is taken
from x = 400 -> 800.

| tree | predicted C Pi | computed | rel | "principal-symbol" formula | local exponent |
|:--|:--|:--|:-:|:--|:-:|
| N2 | pi/6 = 0.5235987756 | 0.5235987839 | 1.6e-8 | -pi/3 | 0.9969 |
| N3 | 7pi/6 = 3.6651914292 | 3.6651914432 | 3.8e-9 | 2pi/3 | 0.9968 |
| N6 | pi/10 = 0.3141592654 | 0.3141592651 | 8.0e-10 | 0 | 0.9966 |

The cases discriminate between the two formulas:
- N2 has the opposite sign and half the size of what the homogeneous formula (-(n-1))^{d_0} x principal would
  give with d_0 = 1.
- N3 differs from it by 7/4.
- N6 has principal integral 0 but Pi = 4/5. It is a tree outside the conjecture's hypothesis that the
  proposition covers.

The exponent drift of about 3e-3 is the O(1/x) correction at x = 400..800, the same size as in the guard's B1
and D rows. The n = 4 top constant (N5) was not computed numerically; there the check is the step-(i) row only.
Nothing beyond the constant (the branch and eps-fit) was re-run; that mechanism is the s-type one with C Pi as
its only input.

## What changes in the paper (loci for the PM)

Replacement text is in `debug/p63_pi_nonzero_paragraph.tex` (184 words). It contains:
1. the proposition, retitled "center derivatives reach the rung", with the label `prop:tree_attain_hom` kept;
2. the (ii) sentence on the constant;
3. the l.1169-1174 sentences, where "sufficient but not necessary" becomes "the hypothesis is Pi != 0", with the
   N6 example;
4. the conjecture's last sentence, which now says the conjecture is proved when Pi != 0 and leaves Pi = 0 open,
   with the d_z^2 - 1/3 example.

Other sentences that change and are not drafted:
- **[ARGUMENT] paragraph (l.1177-1192).** Delete "For non-homogeneous P_v this is not proved. The reason to
  expect it is local." and "With that coefficient nonzero, step (ii) ... would carry over." Also delete "What
  fails for non-homogeneous P_v is the leading coefficient, which is Pi ... rather than the principal-symbol
  average." The local-rule derivation can stay as the explanation of (ii), or be cut to its [MEASURED] records.
  Warning: `tests/test_paper63_local_rule.py::test_paper_sentences_present` pins phrases of this paragraph; check
  which ones before editing.
- **Pointers that still say "homogeneous center derivatives".** These are l.48 (abstract), l.76, l.296, l.343,
  l.932-934, l.2098-2099 ("for non-homogeneous ones it is not shown") and l.2197. Under the summary-surface rule
  they should read "center derivatives with Pi != 0".
- **The open list at l.2207.** It names "non-homogeneous center derivatives". It should become "center
  derivatives with Pi = 0".
- **Outside the paper.** `papers/synthesis/group3_foundations_synthesis.tex` l.1529,
  `papers/synthesis/period_ladder.tex` l.410, and `docs/claim_test_matrix.md` rows 907 (conjecture) and 910
  (proposition: "Pi != 0 sufficient, not necessary" has to go).
- **Guard.** `tests/test_paper63_tree_attainment_l.py` docstring says "Not asserted: non-homogeneous P_v". A
  guard for the extension would be a separate pass (§9): N3 and N6 top constants and one mixed-parity step-(i)
  row. Wrong answers to reject: the principal-symbol constant, and an odd pure r-power.

## What remains open

- Pi = 0 with a nonzero principal integral, e.g. a leaf d_z^2 - 1/3 next to s densities. The review memo argues
  that this case attains one order lower; that is not proved.
- Scale derivatives d_a.
- Exact step-(i) operators at n >= 4 for mixed degrees (only mod one prime here) and anything at n >= 5 with mixed
  degrees.
- Non-axial densities. They are covered by the proof but not computed.

## Where this could be wrong

1. Watson to second order needs two-term expansions of every source, argued by induction. A log at infinity would
   break it, and none can arise from K_nu products and their tails.
2. The angular factor at the root assumes P_0(-d) applied to an asymptotic series can be differentiated termwise.
   This is justified by the explicit integral formulas for E'; it was checked numerically only for n = 3 (the
   N2, N3, N6 constants to 1e-8 or better).
3. Step (i) computations are mod one prime; a wrong reduction would need the least order to drop mod that prime.
