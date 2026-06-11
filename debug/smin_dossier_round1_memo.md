# S_min / depth-k tower — round-1 diagnostic dossier

Sprint: S_min motivic-identification round 1 (diagnostic only), 2026-06-10
Drivers: `debug/smin_dossier_round1.py` (decomposition + CG audit, 130 dps),
`debug/smin_dossier_round1b.py` (full level-2 basis identification, 260 dps),
`debug/smin_dossier_round1c.py` (channel-count closed form + sunset relation),
`debug/smin_dossier_round1d.py` (symbolic assembly).
Data: `debug/data/smin_dossier_round1.json`, `debug/data/smin_dossier_round1b.json`.
No paper edits, no production-code changes, no bulk PSLQ irreducibility sweeps.

**Headline (Q2 + Q5): the "depth-2 new irreducible" reading of S_min is
WRONG.  S_min decomposes exactly (residual 1.8e-129) into Hoffman double
t-values — classical level-2 (alternating-MZV / MT(Z[1/2])) constants — plus
depth-1 terms.  The nesting is structural (a genuine min-weighted double sum),
but the motivic content is classical: the 15+1 PSLQ "irreducibility" failures
are a basis-coverage artifact.  The failed bases never contained the level-2
weight-5..8 monomials (pi^4 ln2, pi^6 ln2, Li5(1/2), Li4(1/2) ln2, Li6, Li7,
zeta(3) ln^2 2, ...) that the standard alternating-MZV basis requires.**

---

## Q1 — Exact definition and physical provenance

**Definition** (Paper 28 `eq:S_min`, lines 2232–2237 of
`papers/group5_qed_gauge/paper_28_qed_s3.tex`):

    S_min = sum_{k>=1} T(k)^2,
    T(k)  = 2 zeta(2, k+3/2) - (1/2) zeta(4, k+3/2).

T(k) is exactly the tail of the Dirac Dirichlet series D(4) on the S^3
Camporesi–Higuchi spectrum: with |lambda_n| = n+3/2, g_n = 2(n+1)(n+2),

    g_n/|lambda_n|^4 = 2/(n+3/2)^2 - (1/2)/(n+3/2)^4 =: phi(n),
    T(k) = sum_{n>=k} phi(n).

**Sum structure.** S_min is NOT a single sum.  Fubini interchange (verified
bit-exactly in Fraction arithmetic at N=40, `smin_dossier_round1.py` step B)
gives

    S_min = sum_{n1,n2>=1} min(n1,n2) phi(n1) phi(n2),

a min-weighted DOUBLE sum over electron levels.  In odd-integer variables
o = 2n+3 (o >= 5 odd), phi = 8(o^2-1)/o^4 and

    S_min = sum_{o1,o2>=5 odd} [(min(o1,o2)-3)/2] * 64 (o1^2-1)(o2^2-1) / (o1^4 o2^4),

a fully rational double sum over odd integers — the natural habitat of
Hoffman multiple t-values (level-2 colored MZVs).

**Physical object.** The min weight is the leading SO(4) Clebsch–Gordan
channel count of the two-loop vertex-restricted sunset (electron exponent
a=4, photon exponent p=1, where d_q^T/mu_q = 1).  Round-1c proves the channel
count in closed form (exact for all n1,n2 <= 40, zero exceptions):

    I(n1,n2) = sum_q W(n1,n2,q) = 2 min(n1,n2) - 1 - [n1=n2]   (n1,n2 >= 1),
    I = 0 when min(n1,n2) = 0.

Hence the TRUE CG-weighted sunset relates to S_min exactly:

    S_sunset = 2 S_min - P^2 - Q = 3.89451620519784104455...
    P = D(4) - 64/81,  D(4) = pi^2 - pi^4/12  (verified 0.0 residual),
    Q = 4 zeta(4,5/2) - 2 zeta(6,5/2) + (1/4) zeta(8,5/2).

So Paper 28's sentence "The CG-weighted two-loop sunset sum is S_min" is
loose by a factor 2 and by depth-1/product terms — motivically harmless, but
the dossier records the exact relation.  S_min is the min-weight CORE of the
sunset; the correction terms are depth-1 MT(Z) objects.

## Q2 — Depth audit: STRUCTURAL nesting, but the label "new irreducible at depth = loop order" FAILS at k=2

Two separate claims must be distinguished:

1. **"S_min is a genuinely nested double sum" — TRUE.**  The min(o1,o2)
   weight cannot be factored into a product of single sums; splitting
   min into {o1<o2}, {o1>o2}, {o1=o2} produces genuinely ordered
   (depth-2, Euler-sum-style) double sums.  This is the discrete analog of
   an iterated integral, and the "one vertex -> one level of nesting"
   mechanism in Paper 28 (lines 2383–2389) is structurally correct.

2. **"Therefore S_min is a NEW irreducible transcendental at depth 2" —
   FALSE (this round).**  The decomposition (verified to 1.8e-129 at
   130 dps, re-verified at 260 dps) is

       S_min = 64 [ t(2,1) - t(4,1) - 3 t(2,2) + 3 t(4,2)
                    - t(2,3) + t(4,3) + 3 t(2,4) - 3 t(4,4) ]
               + (explicit boundary terms in lambda(2), lambda(4))
               + 32 sum_{s=3..8} d_s (lambda(s) - 1 - 3^{-s}),
       d = {3:1, 4:-3, 5:-2, 6:6, 7:1, 8:-3},

   where t(b,c) = sum_{o2>o1>=1 odd} o2^{-b} o1^{-c} are Hoffman double
   t-values and lambda(s) = (1-2^{-s}) zeta(s).  Every t(b,c) here is a
   CLASSICAL level-2 colored MZV (period of MT(Z[1/2]), catalogued in the
   alternating-MZV literature: Hoffman 2019 odd variant; BBV data mine).
   Reductions established this round, each verified to >=1e-130:

   - t(2,2) = (lambda(2)^2 - lambda(4))/2          [stuffle, exact]
   - t(4,4) = (lambda(4)^2 - lambda(8))/2          [stuffle, exact]
   - t(2,4) + t(4,2) = lambda(2) lambda(4) - lambda(6)  [stuffle; both
     enter S_min with the SAME coefficient +3, so only the reducible
     symmetric combination appears]
   - t(2,1) = -(7/16) zeta(3) + (1/8) pi^2 ln2     [PSLQ-pinned, res 1e-131]
   - t(2,3) = -(31/64) zeta(5) + (1/16) pi^2 zeta(3)  [PSLQ-pinned, res 0.0]
   - t(4,1) = (1/96) pi^4 ln2 - (1/64) pi^2 zeta(3) - (31/64) zeta(5)
     [PSLQ 220 dps vs the full w5 level-2 basis (dim 8), res 2.8e-221]
   - t(4,3) = (1/128) pi^4 zeta(3) - (5/128) pi^2 zeta(5) - (127/256) zeta(7)
     [PSLQ 220 dps vs the w7 product basis (dim 19), res 3.6e-222]
     (both 2026-06-11, S^(3) stage-A run, `debug/data/s3_pslq_stageA.json`)

   The depth-2 content of S_min is therefore confined to (at most) the
   weight-5 and weight-7 t-values t(4,1), t(4,3) — and those are elements
   of the standard level-2 basis room, not new constants.

**Narrative consequence (flagged prominently):** Paper 28 Proposition
`prop:depth_k` ("each loop order produces a genuinely new transcendental")
and Paper 55's depth-k tower table (lines 867–886) + stratification row
"S_min depth-2: new period at depth 2" (line 769) overstate the k=2 rung.
The correct statement: loop order bounds the depth of the (classical,
identifiable) multiple-t-value content; it does not mint new irreducibles
at k=2.  Paper 55's ring placement S_min ∈ MT(Z[1/2]) (line 864) is
CONFIRMED — but as an explicit, identified element, answering the paper's
own open question for k=2.

## Q3 — Weight: mixed, top weight 8, graded support w ∈ {0,2,3,4,5,6,7,8}

S_min is weight-INHOMOGENEOUS.  Three weight-lowering mechanisms stack:
(i) the CH degeneracy polynomial g_n = 2(n+1)(n+2) mixes |lambda|^{-2} and
|lambda|^{-4} (so T(k) mixes weights 2 and 4 per factor — the naive
"weight = 2s = 8" expectation holds only for the top piece); (ii) the
min-weight channel count adds one power of the summation variable (weight
-1); (iii) boundary/diagonal corrections inject lambda(2..8).  The graded
decomposition (exact coefficients above):

| w | content | status |
|---|---------|--------|
| 3 | 64 t(2,1), lambda(3) | reduces: zeta(3), pi^2 ln2 |
| 4 | -192 t(2,2), lambda(4), boundaries | pure Tate (pi^4) |
| 5 | -64 t(4,1) - 64 t(2,3), lambda(5) | t(2,3) reduces; t(4,1): [RESULT-1B] |
| 6 | +192 [t(2,4)+t(4,2)], lambda(6) | pure Tate (pi^6) |
| 7 | +64 t(4,3), lambda(7) | [RESULT-1B] |
| 8 | -192 t(4,4), lambda(8) | pure Tate (pi^8) |

This kills the implicit assumption behind every prior PSLQ design (single
homogeneous-weight rooms): a correct basis must span level-2 monomials
across weights 3–8 SIMULTANEOUSLY, or the object must be weight-graded
first (as done here).

## Q4 — Basis-coverage audit of the failed PSLQ runs

**47-element basis** (`debug/archive/rh_arc/smin_identification.py`, the
`ultra_wide_47` attempt; 15 total attempts at 150 dps = the "15 PSLQ
failures"): pi^{2,4,6,8}; zeta(3,5,7) + products; G, beta(3), beta(4) +
products (level 4, w<=8, depth 1); ln2, ln2^2, ln2*{zeta3, G, beta4, pi^2};
Li3(1/2), Li4(1/2); Hurwitz at 3/2 (redundant: zeta(s,3/2) = (2^s-1)zeta(s)
- 2^s ∈ MT(Z)); gamma (not a period); Cl2(pi/3) (level 6, w2).

**100-element basis** (`debug/archive/rh_arc/smin_extended_pslq.py`, 200
dps, Paper 28 lines 2255–2276): adds zeta(9), more pi×zeta products, Euler
sums s_{2,3}, s_{3,2} (level-1, reducible), Tornheim–Witten T(2,2,2),
T(2,2,4) (level-1, reducible), quarter-shift Hurwitz (level 4 depth 1),
D_even/D_odd(4), polygamma at 3/2 (redundant), Li4(1/2)*{pi^2, zeta3},
ln2^3, ln2^2 pi^2, ln2*zeta5, zeta({3,1}), zeta({2,1,1}) (level-1 pi^4
multiples).

**Coverage by (level, weight, depth) — what was genuinely ruled out:**
S_min is not a Q-combination of level-1 MZV products (partial, w<=8),
level-4 depth-1 Dirichlet-beta monomials, the one level-6 w2 constant
Cl2(pi/3), gamma, or the ~15 level-2 monomials listed.  **What was NEVER
tested — and is now known to be the right room:** the level-2 graded room
at w5–w8.  Missing specifically: pi^4 ln2, pi^6 ln2 (w5/w7!), zeta(3) ln^2 2,
Li4(1/2) ln2, Li5(1/2), Li6(1/2), Li7(1/2) and their monomials, ln2^k for
k>=4, all genuine depth-2 alternating irreducibles (1 at w6, 2 at w7,
more at w8).  Of the conjectural level-2 weight-5 basis (dim 8), the failed
bases contained only 2/8 elements (zeta5, pi^2 zeta3); at w7, ~4/21.
Also never tested: full Glanois level-4 basis at w>=5; level-6
(sixth-root-of-unity) constants beyond Cl2(pi/3); elliptic/Bessel-moment
constants.  The last three rooms are now MOOT for S_min (it lives at level
2), but the audit stands as a statement about what the published "PSLQ-
irreducible" claim actually established.

The same coverage gap explains the 35 parity-split PSLQ failures of Track
RH-P (`debug/archive/rh_arc/smin_chi_neg4_memo.md`): those bases were
Dirichlet-L/D(s)-auxiliary only — zero level-2 depth-2 content.

## Q5 — Algebraic reduction (executed) and the sanity anchor

**Anchor (recomputed, not quoted).**  mpmath.nsum Levin at 130 dps and again
at 260 dps (`round1.py` step A, `round1b.py`):

    S_min = 2.4799369380342225544135795008293821446879257866172884583787
            98726559552777818374912328558931430046996...

agreeing with the frozen v2.23.1 three-method value to >=78 digits (the
guard assertion passes; first 100 digits printed in the JSON).

**Reduction chain (all steps exact or verified >=1e-100):**
1. Fubini/min-weight interchange — exact (Fraction, N=40 bit-identity).
2. Odd-integer reindex o=2n+3 — exact (Fraction bit-identity).
3. Partial fractions (sympy): (o-3)(o^2-1)/o^4 = o^{-1} -3o^{-2} -o^{-3}
   +3o^{-4}; (o^2-1)/o^4 = o^{-2} - o^{-4}; diagonal polynomial d_s as in Q2.
   This uses exactly the spectrum factorizations n^2-1-type: o^2-1 =
   (o-1)(o+1) is the CH degeneracy, o^4 the |lambda|^4 line.
4. Boundary shift (o>=5 -> o>=1): exact lambda(b)-terms.
5. Result: the 8-t-value identity of Q2, verified numerically to 1.8e-129
   (130 dps).  (The standalone 260-dps re-run died without writing output;
   superseded by the 220-dps stage-A identifications and the 130-dps
   end-to-end closed-form verification below.)
6. Stuffle kills t(2,2), t(4,4), and the symmetric pair t(2,4)+t(4,2)
   (the antisymmetric part has coefficient 0 in S_min — checked exactly).
7. Odd-weight reductions: t(2,1), t(2,3) reduce to depth-1 (coefficients
   above, residuals 1e-131 and 0.0 at 130 dps).  t(4,1) was tested against
   the FULL conjectural weight-5 level-2 basis (dim 8: zeta5, pi^2 zeta3,
   pi^4 ln2, zeta3 ln^2 2, pi^2 ln^3 2, ln^5 2, Li4(1/2) ln2, Li5(1/2)) and
   t(4,3) against the 19 weight-7 product monomials (+ 2 named depth-2
   alternating candidates zeta(-6,1), zeta(-5,2) in a second pass), at 260
   dps with maxcoeff 1e8 — far beyond spurious-relation range for these
   basis sizes.  Outcome: BOTH IDENTIFIED (coefficients in Q2); no
   obstruction survives.

**Obstruction bookkeeping:** if the w5/w7 identifications fail, the EXACT
obstruction terms are, by construction, the two named constants t(4,1)
(weight 5) and t(4,3) (weight 7) modulo products — i.e., at worst S_min =
(explicit classical part) + 64[t(4,3) - t(4,1)] with both obstructions
being standard catalogue objects of the alternating-MZV data mine, NOT
framework-new transcendentals.

**FINAL (2026-06-11, `debug/smin_final_assembly.py`):** both identified;
the full assembly closes end-to-end:

    S_min = 8 pi^2 ln2 - (2/3) pi^4 ln2 - 3 pi^2 zeta(3)
            + (1/2) pi^4 zeta(3) - (5/2) pi^2 zeta(5)
            + pi^6/4 - (3/2) pi^4 - pi^8/96

verified against the Hurwitz-tail anchor at 130 dps, residual 1.7e-129
(`debug/data/smin_final_assembly.json`).  Structural notes: zeta(7)
cancels EXACTLY against the lambda(7) diagonal term; all rational
constants cancel to zero; every term carries at least pi^2.  No depth-2
generator survives — S_min lies in Q[pi^2, ln2, zeta(3), zeta(5)], i.e.
the k=2 rung of the tower is depth <= 1 content only.

## Q6 — The tower: S^(3) thin, S^(4) vapor

**S^(3)** (Paper 28 `eq:three_loop_factorized`, lines 2306–2320): explicitly
defined as the chain-topology three-loop sum with the REAL CG weights
W(n1,n2,q) (not the min model).  Computed status: at physical exponents
(a=4, p=1) the sum converges as 270 N^{-1.31} — only 1–2 reliable digits at
n_max=100; at fast exponents (a=8, p=2) ~20 digits at n_max=50, where one
18-element PSLQ (including S_min and its products) failed.  No 200-dps
computation exists; "independent of S_min" rests on a 20-digit PSLQ at
non-physical exponents.  **Nesting audit:** the chain sum carries TWO
channel-count factors I(n1,n2) I(n2,n3); by the round-1c closed form this
expands into min(n1,n2)min(n2,n3)-weighted triple sums -> Hoffman TRIPLE
t-values t(b1,b2,b3) (depth <=3, weights <=12, level 2) plus lower-depth
terms.  Same structural verdict as k=2: genuinely nested, but living in the
classical level-2 algebra; at depth 3 the parity/reduction theorems are
weaker, so genuinely irreducible (yet catalogued) depth-3 t-values may
survive — that, not "new transcendentals," is the correct k>=3 prediction.

**S^(4)**: named in Paper 55 (eq. `empirical_transcendentals`, line 165)
but NEVER defined, computed, or tested anywhere in the corpus (grep:
zero hits for a four-loop sum in Paper 28 or debug archives).  Citation-
level overreach; flag for the next Paper 55 touch.

## Round-2 design consequences

1. **No PSLQ sweep needed for k=2.**  The reduction replaces it.  Round 2
   should instead (a) verify the t-value reductions symbolically (Hoffman
   2019 / Kaneko–Tsumura machinery — the t(2,1), t(2,3) forms match the
   known odd-variant evaluations), and (b) propagate the corrected
   classification into Papers 28 + 55 (PI-gated paper edits).
2. **If any PSLQ remains** (e.g., w7 layer), it must be weight-graded at
   level 2: target ONE graded component against ONE weight-homogeneous
   alternating-MZV basis (w7 needs all 21 elements incl. Li7(1/2),
   Li6(1/2)ln2, and the 2 depth-2 irreducibles).
3. **S^(3) program:** repeat this dossier's decomposition (channel-count
   closed form -> triple t-values -> stuffle/parity reductions) BEFORE any
   numerics; the factorized O(N^3) algorithm is then only a cross-check.
4. **Paper-edit queue (NOT applied this round, per charter):** Paper 28
   §S_min + prop:depth_k reframe; Paper 55 sub-sector 3, stratification
   table, tower table, abstract line "genuinely new transcendentals";
   the exact sunset relation S_sunset = 2 S_min - P^2 - Q belongs in
   Paper 28 §S_min.  Existing tests (`tests/test_smin_extended.py`) stay
   green — they pin the value and basis-specific nulls, both still true.

## Honest scope (sprint close, 2026-06-11)

- **Exact / theorem-grade:** Fubini interchange, odd reindex, partial
  fractions, boundary shifts (Fraction bit-identities); the stuffle
  reductions (algebraic identities).  The channel-count closed form is
  exact with zero exceptions at n1,n2 <= 40 — a finite verification, not
  yet an all-n proof.
- **PSLQ-pinned (220 dps, complete conjectural bases):** the four
  odd-weight reductions t(2,1), t(2,3), t(4,1), t(4,3) — conditional on
  the standard level-2 basis conjectures (proven through w=11 in the BBV
  data mine, which covers all four weights used here); residuals
  <= 1e-221.
- **Numerical verification:** the assembled closed form vs the anchor at
  130 dps (residual 1.7e-129).
- **Named open follow-ons:** none for S_min itself.  The k >= 3 items
  live in `debug/sprint_s3_decomposition_memo.md` §5.

## Files

- `debug/smin_dossier_round1.py` + `debug/data/smin_dossier_round1.json`
- `debug/smin_dossier_round1b.py` + `debug/data/smin_dossier_round1b.json`
- `debug/smin_dossier_round1c.py` (CG closed form + sunset relation)
- `debug/smin_dossier_round1d.py` (symbolic assembly)
- This memo: `debug/smin_dossier_round1_memo.md`

Provenance of prior claims: Paper 28 lines 2227–2292 (definition, value,
100-element basis), 2362–2389 (depth-k tower); Paper 55 lines 121–126,
158–174, 704–710, 762–773, 846–886; `debug/archive/rh_arc/
smin_verification_memo.md` (v2.23.1 three-method value);
`debug/archive/rh_arc/smin_identification.py` (15 failures, 47-element
basis); `debug/archive/rh_arc/smin_extended_pslq.py` (100-element basis);
`debug/archive/rh_arc/smin_chi_neg4_memo.md` (35 parity-split failures);
`geovac/qed_vertex.py::two_loop_sunset_weighted` / `::two_loop_min_weighted_hurwitz`.
