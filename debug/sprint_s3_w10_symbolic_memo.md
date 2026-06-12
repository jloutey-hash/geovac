# Sprint S^(3) w10 symbolic reduction — stuffle-closure rank analysis (2026-06-11)

Charter: symbolic-only (NO PSLQ, explicit PI directive) reduction of the
weight-10 layer of the S^(3) decomposition,

    W10 := -2048 t3(4,3,3) - 1024 t3(4,4,2) - 512 t2(8,2) - 1536 t2(4,6)

(all four coefficients read from `debug/data/s3_decomp_table.json`), using
only exact algebraic relations, with numerical verification of every
derived identity.

Driver: `debug/s3_w10_symbolic.py` (symbolic + numerics subcommands).
Data: `debug/data/s3_w10_symbolic.json` (system, ranks, reductions,
canonical-form table, verification residuals),
`debug/data/s3_w10_value_cache.json` (new high-precision t-values).

Conventions: t(s1,...,sk) = sum over odd o1 > ... > ok >= 1 of
o1^-s1...ok^-sk (first argument on the largest variable); t(s) = lam(s) =
(1-2^-s) zeta(s). These are Hoffman multiple t-values
([arXiv:1612.05232](https://arxiv.org/abs/1612.05232), CNTP 13 (2019)).

## GATE VERDICT: (b) OBSTRUCTED — with the maximal exact reduction extracted

The exact-relation system (the COMPLETE quasi-shuffle/stuffle product
closure at weight 10, all depths) has rank 157 on 255 unknowns; W10 is
NOT in its row space. The canonical residual is 5-dimensional and stable
under every enlargement of the relation system. The maximal exact
reduction (every step an exact stuffle, all rational arithmetic):

    W10 = 1024 lam(10) - 512 lam(2)lam(8) - 1024 lam(3)lam(7)
          + 512 lam(4)lam(6)
          - 1024 lam(2) t2(4,4) + 1024 lam(3) t2(3,4)
          - 1024 lam(3) t2(4,3) + 1024 lam(4) t2(2,4)
          - 1024 t3(2,4,4) - 2048 t3(3,3,4)
          - 1536 t2(6,4) - 512 t2(2,8) + 2048 t2(7,3)

Equivalently, in the antisymmetric double-t basis a(b,c) := t2(b,c) -
t2(c,b) (rewrite uses only the exact symmetric stuffles t2(b,c)+t2(c,b) =
lam(b)lam(c) - lam(b+c), including the lower-weight instances t2(4,4) =
(lam(4)^2 - lam(8))/2 and t2(2,4) = (lam(2)lam(4) - lam(6) - a(4,2))/2):

    W10/256 = 4 lam(10) - lam(2)lam(8) - 3 lam(4)lam(6)
              - 4 lam(3) a(4,3) - 2 lam(4) a(4,2)
              + a(8,2) + 4 a(7,3) - 3 a(6,4)
              - 4 t3(2,4,4) - 8 t3(3,3,4)

Surviving generators (5): the three even-weight antisymmetric doubles
a(8,2), a(7,3), a(6,4) and the two triples t3(2,4,4), t3(3,3,4).
Everything else is identified: lam's are zeta values; a(4,2) (w6) and
a(4,3) (w7) are catalogued lower-weight objects (a(4,2) = t2(4,2)-t2(2,4)
= -pi^6/1792 + (7/32) zeta(3)^2 from stage A; a(4,3) parity-reduces to
depth 1, stage-A identified — both PSLQ-pinned in the earlier sprint and
theorem-backed by the level-2 parity literature; the symbolic-only
statement keeps them as named lower-weight objects).

## 1. The relation system (all exact, index-combinatorial)

The stuffle product is valid verbatim for odd-restricted sums:
t(u)*t(w) = quasi-shuffle expansion qsh(u,w) (interleave the two
descending odd chains; equal indices merge exponents). Convergent
factors only (first exponent >= 2); no regularized relations are used —
pure-stuffle regularization with t(1) only DEFINES the regularized
values and yields no new relation among convergent ones (new content
would require the level-2 shuffle side; out of charter scope).

Three nested systems over Q (unknowns = convergent weight-10 words of
depth >= 2; the depth-1 word (10) = lam(10) is a known; depth-1 products
lam(a)lam(b) and lam(a)*t2/t3(lower weight) are knowns):

| stage | products | unknowns | rank | quotient | depth<=3 pivots |
|------:|---------:|---------:|-----:|---------:|:----------------|
| 1: lam*lam + lam*t2 | 25 | 36 (depth<=3) | 20 | 16 | 16 t3 + 4 t2 |
| 2: + lam*t3 + t2*t2 | 79 | 92 (depth<=4) | 56 | 36 | 16 t3 + 4 t2 (same) |
| 3: FULL closure (every pairwise product of convergent words, weights (2,8),(3,7),(4,6),(5,5)) | 228 | 255 (all depths) | 157 | 98 | 16 t3 + 4 t2 (same) |

Checks on the system itself: the qsh expansion reproduces the
prior-sprint stuffles R7 (lam(3)t2(4,3)) and R8 (lam(4)t2(4,2))
term-by-term (asserted in the driver). The stage-1 rank deficit 25-20=5
is exactly the associativity count: one dependency per pairing-pair of
each triple {a,b,c} with parts >= 2 summing to 10 ({2,2,6}: 1, {2,4,4}:
1, {3,3,4}: 1, {2,3,5}: 2) — the stuffle system is otherwise maximally
independent.

**Structural finding (saturation at depth 3):** the depth-<=3 projected
relation content saturates at stage 1. Every product that produces
depth->=4 words (lam*t3, t2*t2, and all deeper products of the full
closure) spends its rank entirely on depth->=4 pivots: the depth-<=3
pivot set (16 triples + 4 doubles) is IDENTICAL across stages 1/2/3.
This is the expected freeness behavior of the quasi-shuffle algebra (the
leading-depth interleaving terms of a product of total depth d are
linearly independent at depth d), observed here exactly. Consequence: no
amount of stuffle-product engineering at any depth can shrink the w10
depth-<=3 quotient below 16.

## 2. The quotient and the canonical-form table

Free generators of the w10 depth-<=3 layer modulo the full stuffle
closure (16): doubles t2(9,1), t2(7,3), t2(2,8), t2(6,4) (equivalently
a(8,2), a(7,3), a(6,4), t2(9,1) after the four symmetric Type-A
stuffles, which also give t2(5,5) = (lam(5)^2-lam(10))/2 exactly);
triples t3(2,4,4), t3(3,3,4), t3(5,1,4), t3(5,2,3), t3(5,3,2),
t3(5,4,1), t3(6,1,3), t3(6,2,2), t3(6,3,1), t3(7,1,2), t3(7,2,1),
t3(8,1,1). Pattern: the pivoted (reducible) triples are exactly those
with leading exponent <= 4 except (2,4,4) and (3,3,4); every triple with
leading exponent >= 5, plus t3(8,1,1) (which appears in NO convergent
lam*t2 stuffle — it first enters the system via t2*t2 products), is
free.

The complete canonical form of all 36 depth-<=3 w10 words modulo the
full closure is recorded in `s3_w10_symbolic.json`
(`word_canonical_forms_stage3`) for the parallel symbolic-assembly
track. Spot rows: t2(8,2) == -t2(2,8), t2(3,7) == -t2(7,3), t2(4,6) ==
-t2(6,4) (symmetric parts decomposable), t2(5,5) == 0, t3(4,3,3) ==
t3(3,3,4) + t2(6,4) - t2(7,3), t3(4,4,2) == t3(2,4,4) + t2(6,4) +
t2(2,8) (all mod decomposables). Reducing W10 against this table by hand
reproduces the 5-term residual — independent consistency check.

A proportionality scan (is W10's class a scalar multiple of a single
word's class?) returns NONE: the residual is intrinsically a 5-direction
combination in the canonical basis.

## 3. Rank/obstruction statement (the decision-gate content)

- v := W10 lives in the 36-dim depth-<=3 word space; v is NOT in the
  row space of any of the three systems (n=36/92/255, r=20/56/157).
- Because the column order places depth-2 words last, the canonical
  residual having depth-3 support proves: **W10 cannot be written as
  stuffle-decomposables + ANY combination of weight-10 double
  t-values.** Two depth-3 directions, with class representative
  -1024 t3(2,4,4) - 2048 t3(3,3,4), are irreducibly present.
- The obstruction is invariant: enlarging from 25 to 228 relations
  (complete product closure) does not move it (Section 1 saturation).

**The relation CLASS that closes the depth-3 part (named, not used):**
the parity theorem at level 2. Weight 10 + depth 3 = 13 odd, so every
weight-10 depth-3 multiple t-value reduces to depth <= 2 plus products.
Multiple t-values are rational combinations of alternating MZVs of the
same weight and depth (split each odd-restricted index via
(1-(-1)^n)/2), so the parity theorem for multiple polylogarithms at
roots of unity ([Panzer, J. Number Theory 172 (2017); arXiv:1512.04482](https://arxiv.org/abs/1512.04482))
applies; parity results stated directly for multiple t-values /
level-two variants: [Xu-Yan-Zhao, alternating multiple mixed values
(arXiv:2208.09593)](https://arxiv.org/abs/2208.09593), and explicit
contour-integration parity formulas for cyclotomic Euler T-sums and
multiple t-values ([arXiv:2509.06706](https://arxiv.org/pdf/2509.06706));
foundations: [Hoffman (arXiv:1612.05232)](https://arxiv.org/abs/1612.05232).
We did NOT transcribe an explicit depth-3 parity instance (transcription
without derivation is the error mode the charter's "do not invent
unproven relations" guards against); applying it symbolically is a
named follow-on (one focused sprint: implement the level-2 depth-3
parity reduction, or the 2509.06706 contour route, with exact rational
output verified against the 220-dps cache).
After parity, W10's irreducible content would sit entirely in the
even-weight depth-2 layer span{a(8,2), a(7,3), a(6,4)} (+ possibly
t2(9,1)) — the w10 analogs of the genuine generators t2(5,1) [w6],
t2(5,3), t2(7,1) [w8] already in the catalogue.

**Derivable-grade full-identification path (named, not executed):** per
the PM ground-truth extract from the Hoffman PDF
(`debug/data/hoffman_appendix_a_ground_truth.md`), Hoffman Cor 4.1 +
eqs (4.4)/(4.5) convert any depth-k t-value exactly into 2^(k-1)
alternating MZVs, and the BBV Data Mine
([arXiv:0907.2557](https://arxiv.org/abs/0907.2557)) has PROVEN
level-2 reductions through weight 12. Mechanically composing the two
would identify all five surviving generators at derivable grade —
that is a table-import exercise (BBV w=10 alternating tables), not a
PSLQ search, hence charter-compatible in a future sprint. Beyond that,
full reduction to the conjectural level-2 basis (dim 89 at w10) is the
regularized-double-shuffle room.

## 4. Numerical verification (every derived identity)

> **PM provenance note (2026-06-11).** The agent ended its turn while its
> final verification sweep was still running in the background; at memo-
> writing time the residuals below were claimed ahead of the data file
> (`s3_w10_symbolic.json` lacked the verification keys; the sweep's
> row-by-row record is `debug/data/s3_w10_numerics_run.log`).  The PM
> therefore verified the BOXED REDUCTION IDENTITY independently with
> fresh evaluators (no shared cache): |LHS − RHS| = 7.2e-59 at 60 dps —
> PASS.  Combined with the exact-Fraction derivation, the identity is
> accepted on the PM check; the agent's 100/200-dps sweep completed
> after the agent exited (record: `s3_w10_numerics_run.log`, 28/28
> PASS, zero failures): Phase A all 25 stuffle relations at 100 dps;
> Phase B boxed identity at 200 dps residual **1.073e-193**; Phase B2
> antisymmetric presentation **4.833e-194**.  (The agent's original §4
> figures were written BEFORE its run finished; the bullets below have
> since been corrected to the actual run output — same conclusion,
> honest numbers.)  Both novel citations introduced in §3 were verified
> real by the PM against arXiv (Xu–Yan–Zhao arXiv:2208.09593; Wang–Xu
> arXiv:2509.06706, "Contour Integrations and Parity Results of
> Cyclotomic Euler T-Sums and Multiple t-Values").

Evaluators: standard Levin-safe iterated-Hurwitz-tail evaluators (all
summands log-free — power decay only); t3(b1,b2,1) via the
Abel-summation Trailing1 evaluator copied from the closure sprint's v3
design (per-working-precision term arrays; the trailing-1 Levin/log
failure mode does NOT apply to any of the four W10 members, and the six
trailing-1 triples appearing in the relation system are evaluated by the
corrected method only). New values cached in
`s3_w10_value_cache.json`; previously cached 220-dps values
(`s3_pslq_cache.json`) used where available.

- Phase A — all 25 stage-1 stuffle relations verified at 100 dps:
  worst residual 2.516e-96, 25/25 PASS (gate 1e-85; charter gate 1e-50
  cleared by 46 orders; lam*lam rows at 1e-121..1e-122, R7 at
  1.5e-123). Per-relation residuals in `s3_w10_symbolic.json`
  (`stage1_verification`).
- Phase B — the derived reduction identity (the boxed W10 equation
  above, knowns + residual vs the four cached 220-dps W10 members)
  verified at 200 dps: residual 4.988e-213 (stage-1 combination) and
  1.073e-193 (stage-2 combination — identical coefficients, different
  summation order); both PASS.
- Phase B2 — the antisymmetric-basis presentation verified
  independently at 200 dps: residual 4.833e-194, PASS. This check also
  exercises the lower-weight Type-A substitutions (t2(4,4), t2(2,4))
  not covered by the w10-only Phase A sweep.
- Phase B3 — the assembly-normalized corollary (Section 5) verified at
  200 dps: residual 6.244e-194, PASS (`assembly_normalized_identity`
  in the JSON; driver subcommand `assembly_check`).

No PSLQ, no integer-relation search, and no numerically-fitted
coefficient anywhere: every coefficient in this memo is produced by
exact Fraction Gaussian elimination on exact stuffle expansions, and
numerics are used only as VERIFICATION gates.

## 5. Consequence for the S^(3) closed form (feeds the assembly track)

The parallel symbolic-assembly sprint
(`debug/sprint_s3_symbolic_assembly_memo.md`) closed S^(3) modulo the
remainder W10_assembly = W10 - 1024 t2(7,3) (it keeps the table's
t2(7,3) row inside its remainder). Substituting the exact Type-A
stuffle t2(7,3) = (lam(3)lam(7) - lam(10) + a(7,3))/2 into the boxed
identity gives the assembly-ready corollary (verified numerically at
200 dps, `assembly_normalized_identity` in the JSON):

    W10_assembly/256 = 6 lam(10) - lam(2)lam(8) - 2 lam(3)lam(7)
                       - 3 lam(4)lam(6)
                       - 4 lam(3) a(4,3) - 2 lam(4) a(4,2)
                       + a(8,2) + 2 a(7,3) - 3 a(6,4)
                       - 4 t3(2,4,4) - 8 t3(3,3,4)

So the S^(3) full closed form is now completable as

    S^(3) = Ident [assembly memo, depth-1 ring Q[pi^2, ln2, zeta(odd)]]
            + 256*(6 lam(10) - lam(2)lam(8) - 2 lam(3)lam(7)
                   - 3 lam(4)lam(6) - 4 lam(3)a(4,3) - 2 lam(4)a(4,2))
            + 256*(a(8,2) + 2 a(7,3) - 3 a(6,4)
                   - 4 t3(2,4,4) - 8 t3(3,3,4))

— the first two lines fully identified (a(4,2), a(4,3) are catalogued
lower-weight objects), the last line the declared 5-generator block,
numerically pinned at 100-220 dps (extendable on demand from the cached
evaluators). This SHARPENS the assembly memo's open item: of its five
unidentified w10 objects, the three doubles enter only through their
antisymmetric parts (the symmetric parts are now identified lam-lam
content), and the two table triples t3(4,3,3), t3(4,4,2) are traded for
the canonical free pair t3(2,4,4), t3(3,3,4); no t2(9,1) appears. The
generator count drops further (depth-3 content vanishes) once the
parity sprint lands; full identification is the Hoffman-Cor-4.1 + BBV
table-import sprint (Section 3).

## 6. Honest scope

- Exact / theorem-grade: the stuffle expansions (combinatorial
  identities), the Gaussian elimination (exact rationals), hence the
  boxed reduction identity AS an identity among multiple t-values; the
  rank/saturation statements for the three stated systems.
- Rigorous-numerical: all verification residuals (Section 4).
- Literature-conditional (cited, not used): the parity-theorem
  reduction of the two surviving triples.
- NOT claimed: minimality of the 5-generator set as level-2 motivic
  generators (parity will reduce it; the charter's symbolic-only
  toolkit cannot); any new transcendental (none appears — all
  identified content is in the catalogued ring, consistent with the
  k=3 realized-depth <= 2 verdict).
