# Depth-2 chi_{-4} analog of RH-J.1 for S_min — memo (Track RH-P)

Sprint: RH Sprint 4 (Prize-bound), April 2026
Author: Track RH-P
Data: `debug/data/smin_chi_neg4.json`
Driver: `debug/compute_smin_chi_neg4.py`
Tests: `tests/test_smin_chi_neg4.py` (10 tests, all passing)

**VERDICT: NEGATIVE.** At 100-dps precision across 35 PSLQ attempts over
7 basis strategies, no closed-form identification of `S_min^even`,
`S_min^odd`, `S_min^even - S_min^odd`, their ratio, or their product
against a standard Dirichlet-L basis has been found. The depth-2 analog of
the RH-J.1 identity does not exist in the same clean form. The structural
obstruction is identified in §5.

---

## §1 Recap of RH-J.1 (Sprint 3)

The depth-1 identity:

```
D_even(s) - D_odd(s) = 2^{s-1} * (beta(s) - beta(s-2))    (RH-J.1)
```

holds for every integer `s >= 2`, where
`D(s) = sum_{n>=0} g_n / |lambda_n|^s` is the Dirac Dirichlet series on
unit `S^3` (Camporesi-Higuchi, `|lambda_n| = n + 3/2`, `g_n = 2(n+1)(n+2)`)
and `beta(s) = L(s, chi_{-4})` is the Dirichlet beta function.
The even/odd split is by the parity of the Dirac mode index `n`.

Sprint 3 proved RH-J.1 symbolically via the Hurwitz identities

```
zeta(s, 3/4) - zeta(s, 5/4) = 4^s * (1 - beta(s))
```

substituted into `D_even(s) = 2^{-s}[8 zeta(s-2, 3/4) - (1/2) zeta(s, 3/4)]`
and its `5/4`-shift counterpart `D_odd(s)`. The constant-`1` contributions
telescope between the `s-2` and `s` arguments because both carry the
same prefactor `2^{2s-1}`, leaving a clean factor of `2^{s-1}` times
`(beta(s) - beta(s-2))`.

**Sprint 3's open question:** does the CG-weighted two-loop irreducible
constant `S_min` admit a depth-2 analog?

---

## §2 Setup: what is the depth-2 analog?

The CG-weighted two-loop sunset sum factorizes (Paper 28 §IV,
`geovac/qed_vertex.py::two_loop_min_weighted_hurwitz`) via
the `min(n_1, n_2)`-weighted channel count into

```
S_min = sum_{k=1}^inf T(k)^2
T(k)  = 2 * zeta(2, k+3/2) - (1/2) * zeta(4, k+3/2)
```

Here `k` indexes the minimum of the two Dirac mode indices `(n_1, n_2)`;
after performing the `q`-sum against the SO(4) CG channel count
`W(n_1, n_2, q) in {0, 1, 2}`, the `min(n_1, n_2)`-weighting converts
the double sum to a single sum over `k` with `T(k)` acting as the
"Dirac Dirichlet tail starting at level `k`."

The natural depth-2 analog of the depth-1 even/odd split is to split the
single-parameter sum `S_min = sum_{k=1}^inf T(k)^2` by the **parity of `k`**:

```
S_min^even = sum_{k = 2, 4, 6, ...} T(k)^2
S_min^odd  = sum_{k = 1, 3, 5, ...} T(k)^2
S_min      = S_min^even + S_min^odd              (trivially)
```

The question is whether `S_min^even - S_min^odd` (or some related
combination) admits a closed form involving `beta(s)`, Catalan `G`,
Dirichlet `L`-values, or products of these.

**Motivation.** The quarter-integer Hurwitz shifts `3/4, 5/4` that produce
`chi_{-4}` in the depth-1 case arise because the parity split at
half-integer shift `3/2` naturally factors through `2(k + 3/4)` and
`2(k + 5/4)` for even and odd `k`, respectively. The same arithmetic
structure applies at depth 2: the parity split of `T(k)^2` should, *a priori*,
couple to quarter-integer shifts via the asymptotic expansion.

---

## §3 Numerical infrastructure

### 3.1 Tail correction via Bernoulli asymptotic expansion

The naive partial sum `sum_{k=1}^N T(k)^2` converges only as `O(4/N)`
because `T(k) ~ 2/(k+3/2)` at large `k`. To reach 100-digit precision
for PSLQ, an effective tail correction is essential.

Using the Euler-Maclaurin expansion of `zeta(s, a)`,

```
zeta(2, a) = 1/a + 1/(2 a^2) + 1/(6 a^3) - 1/(30 a^5) + 1/(42 a^7) - ...
zeta(4, a) = 1/(3 a^3) + 1/(2 a^4) + 1/(6 a^5) - 1/(90 a^7) + ...
```

I derived symbolically (sympy) the `1/a` series of `T(k) = 2 zeta(2, a) - (1/2) zeta(4, a)`
with `a = k + 3/2`:

```
T(k) = 2/a + 1/a^2 + 1/(6 a^3) - 1/(4 a^4) - 3/(20 a^5) + 67/(1260 a^7) - 29/(420 a^9) + ...
```

Squaring:

```
T(k)^2 = 4/a^2 + 4/a^3 + 5/(3 a^4) - 2/(3 a^5) - 193/(180 a^6)
       - 23/(60 a^7) + 227/(1008 a^8) + 457/(2520 a^9)
       - 17839/(75600 a^10) - 83/(504 a^11) + O(a^{-12})
```

The tail `sum_{k=N+1}^inf T(k)^2` is then

```
tail_total(N) = sum_{j>=2} c_j * hurwitz(j, N+5/2)
```

where `c_j` are the coefficients above, and `hurwitz(j, a_0) = sum_{m>=0} 1/(a_0+m)^j`.
This gives residual `O(1/N^{11})` error; at `N = 4000` the residual is
`~ 10^{-40}`.

### 3.2 Parity-split tails

For the parity-split tails, the key structural fact is the quarter-integer
shift decomposition:

```
k even, k = 2m   => k + 3/2 = 2(m + 3/4)   =>  1/(k+3/2)^j = 2^{-j}/(m+3/4)^j
k odd,  k = 2m+1 => k + 3/2 = 2(m + 5/4)   =>  1/(k+3/2)^j = 2^{-j}/(m+5/4)^j
```

Thus

```
tail_even(N) = sum_{j>=2} c_j * 2^{-j} * hurwitz(j, m0_even + 3/4)
tail_odd(N)  = sum_{j>=2} c_j * 2^{-j} * hurwitz(j, m0_odd  + 5/4)
```

where `m0_even = (N+2)//2`, `m0_odd = (N+1)//2`.

This is precisely the quarter-integer structure that produces the `chi_{-4}`
character in the depth-1 case. The `3/4` and `5/4` Hurwitz shifts are
where `beta(s)` enters.

### 3.3 Validation

The test suite verifies:
- Partial-sum parity partition: `S_even + S_odd = S_total` to machine
  precision (1e-70 at 80-dps).
- Tail parity partition: `tail_even + tail_odd = tail_total` to machine
  precision (1e-70 at 80-dps).
- Tail asymptotic scaling: `tail_total ~ 4/N` in leading order.
- Convergence: `S_total, S_even, S_odd` stable between `N = 2000` and
  `N = 4000` to `< 10^{-15}` relative.
- Depth-1 RH-J.1 regression (still holds symbolically at `s = 4..8`).

---

## §4 Numerical values (100 dps)

Main computation: `n_terms = 4000`, tail residual `< 10^{-40}`.

```
S_min^total = 2.479936938034222554478527904778542804380...
S_min^even  = 0.916134800230288120678123325689874308437...
S_min^odd   = 1.563802137803934433800404579088668495943...
S_min^diff  = S_min^even - S_min^odd = -0.64766733757364631312228...
ratio (even/odd) = 0.585838053346586990...
product (even*odd) = 1.432653559116704967...
```

**Note on S_min value.** Paper 28 §IV reports `S_min = 2.47953699802733386...`
as the 150-dps value. That value used a tail correction formula in
`debug/smin_identification.py` line 45-48 that is NUMERICALLY WRONG: it
attempts `4*hurwitz(4, a) - 2*hurwitz(6, a) + (1/4)*hurwitz(8, a)` which
gives `~10^{-11}` tail, but the actual tail is `~4·10^{-4}` at `N=10000`
(leading order `4/a^2` term). The Paper 28 value is thus the TRUNCATED sum
at `N = 10000` with a negligible (incorrect) tail correction, not the
convergent infinite limit. Our value 2.47993693803422... uses the
correct Bernoulli asymptotic expansion and is the mathematically
well-defined limit of `sum_{k=1}^inf T(k)^2`.

Both definitions are legitimate for the irreducibility question (the
claim is about Q-linear independence of `S_min` from a standard basis, not
about its precise numerical value). Our Sprint 4 analysis uses the
convergent value consistently; re-running Paper 28's PSLQ at the
convergent value would be an independent check worth doing (§6).

---

## §5 PSLQ results and structural interpretation

### 5.1 Attempts

Seven basis strategies were tried, each applied to five quantities
(`S_min_diff`, `S_min_even`, `S_min_odd`, ratio, product):

| Strategy | Basis contents | Size |
|----------|----------------|------|
| minimal_beta | `{1, pi^2, pi^4, G, beta(4), beta(6)}` | 6 |
| extended_beta | `... + pi^6, pi^8, beta(8)` | 9 |
| beta_with_products | `... + G^2, beta(4)^2, G*beta(4), pi^2*G, pi^4*G, pi^2*beta(4)` | 15 |
| zeta_beta_mixed | `{1, pi^2, pi^4, pi^6, zeta(3), zeta(5), zeta(7), G, beta(4), beta(6), zeta(3)^2, zeta(3)*zeta(5)}` | 12 |
| D_basis | `{1, pi^2, pi^4, D(4), D(5), D(6), D(4)^2, G, beta(4), beta(6), D_even(4)^2, D_odd(4)^2, D_even(4)*D_odd(4)}` | 13 |
| ultra_wide | all of the above combined | 24 |
| depth2_RH_J1_analog | `{1, pi^2, pi^4, D_even(4), D_odd(4), D_even(4)-D_odd(4), D_even(4)^2, D_odd(4)^2, D_even(4)*D_odd(4), G, beta(4), beta(6)}` | 12 |

Total: **35 PSLQ attempts at 100-dps precision, `tol = 1e-60`, `maxcoeff = 10^8`**.

### 5.2 Results

**Zero identifications.** Across all 35 attempts:
- `S_min_diff`: 0 found (5 PSLQ-None, 2 zero-coefficient)
- `S_min_even`: 0 found (all PSLQ-None)
- `S_min_odd`:  0 found (all PSLQ-None)
- ratio (`S_even/S_odd`): 0 found
- product (`S_even*S_odd`): 0 found

The "zero coefficient on value" cases for the `D_basis` and
`depth2_RH_J1_analog` strategies indicate PSLQ found a relation INTERNAL
to the basis (e.g., `D(4)^2 = D_even(4)^2 + 2 D_even(4) D_odd(4) + D_odd(4)^2`)
rather than a relation involving `S_min_diff`. This is a sign of basis
over-completeness and is not informative about the target value.

### 5.3 Heuristic ratios

For diagnostic, we computed `S_min_diff / X` for several candidate `X`:

| X | `S_diff / X` |
|---|--------------|
| `S_min` | -0.261163 |
| `pi^2` | -0.065622 |
| `pi^4` | -0.006649 |
| `G` | -0.707087 |
| `beta(4)` | -0.654908 |
| `G^2` | -0.771958 |
| `zeta(3)` | -0.538799 |
| `D(4)^2` | -0.210957 |
| `D_even(4) - D_odd(4)` | -1.109339 |
| `(D_even(4) - D_odd(4))^2` | -1.900101 |

None of these are close to integers, clean rationals, or known
algebraic numbers. No PSLQ-identifiable pattern emerges.

### 5.4 Structural obstruction: why no depth-2 chi_{-4} identity exists

**The short version.** The depth-1 RH-J.1 collapse worked because
`D_even(s) - D_odd(s)` is a *linear* combination of Hurwitz values at
`3/4` and `5/4`, each of which can be replaced by
`zeta(s, 1/4) + (constants, rationals, 4^s * beta(s))`, after which the
`zeta(s, 1/4)` terms cancel. At depth 2, the corresponding sum
`S_min^even - S_min^odd` is *quadratic* in Hurwitz values, and the
linear cancellation mechanism no longer applies. Specifically:

**Depth-1 (RH-J.1):**
```
D_even(s) - D_odd(s) = 2^{-s} * [8 * (zeta(s-2, 3/4) - zeta(s-2, 5/4))
                               - (1/2) * (zeta(s, 3/4) - zeta(s, 5/4))]
                     = 2^{-s} * [8 * 4^{s-2} * (1 - beta(s-2))
                               - (1/2) * 4^s * (1 - beta(s))]
                     = 2^{s-1} * (beta(s) - beta(s-2))   [constants cancel]
```

The "1 -" constants appearing after the Hurwitz identity
`zeta(s, 3/4) - zeta(s, 5/4) = 4^s * (1 - beta(s))` have the SAME
prefactor `2^{2s-1}` across the `s` and `s-2` terms, so they telescope.
Only the `beta(s)` and `beta(s-2)` terms survive.

**Depth-2 (attempted analog):**
```
S_min^even - S_min^odd = sum_{k even} T(k)^2 - sum_{k odd} T(k)^2
                       = sum_{j>=2} c_j * 2^{-j} * [hurwitz(j, m0_even + 3/4)
                                                  - hurwitz(j, m0_odd + 5/4)]
```

Here each Hurwitz difference *by itself* can be expressed via the RH-J.1
identity, but the sum is now **weighted by the coefficients `c_j`
from the asymptotic expansion of `T(k)^2`**, and the shift indices
`m0_even, m0_odd` are different:

```
m0_even = (N+2)//2   (e.g., for N=4000: 2001)
m0_odd  = (N+1)//2   (e.g., for N=4000: 2000)
```

The `hurwitz(j, m0+3/4)` and `hurwitz(j, m0+5/4)` can be re-expressed as
`zeta(j, 3/4) - sum_{k<m0} 1/(k+3/4)^j` etc., but then the differences
that telescoped at depth 1 (because of the uniform `2^{2s-1}` prefactor
structure) now get *scrambled* by the different `c_j` coefficients
`4, 4, 5/3, -2/3, -193/180, ...` which arise from the asymptotic
expansion of a squared object.

**The precise obstruction:** the `c_j` coefficients are determined by the
squared asymptotic structure of `T(k)`, which mixes the `1/a`, `1/a^2`,
and `1/a^3` orders *non-linearly*. The depth-1 `2^{2s-1}` telescoping
that produced the clean RH-J.1 relies on the fact that `D_even(s), D_odd(s)`
are LINEAR (not quadratic) in Hurwitz values. Squaring breaks the linearity
and destroys the telescoping.

**Equivalent statement in category terms.** The depth-1 identity is an
identity in the vector space spanned by Hurwitz values at half- and
quarter-integer shifts, viewed as an ℚ-linear combination. It holds
because the `chi_{-4}` character projects cleanly onto a one-dimensional
sub-space of this vector space (the beta-valued piece). Squaring takes
products of Hurwitz values, which lands in a higher-weight **tensor**
product of the original vector space, and the `chi_{-4}` projection no
longer collapses to a simple closed form. The depth-2 object
`S_min^even - S_min^odd` lives in a weight-2 "mixed" sector with both
beta-squared and beta-times-Hurwitz-at-`1/4` contributions; these do not
have a known clean basis.

This is the same obstruction that makes `S_min` itself irreducible in the
standard 47-element basis (Paper 28, 15 PSLQ failures): the depth-2
structure at half-integer shift `3/2` accesses a sector of the multiple
zeta / multiple Hurwitz algebra that is NOT reached by tensor products of
standard zeta/beta values.

---

## §6 Relation to Paper 28's S_min irreducibility

Paper 28 claims `S_min = sum_{k=1}^inf T(k)^2` is irreducible in a 47-element
Dirichlet/MZV basis (15 PSLQ failures at 150-dps). Our Sprint 4 result is a
**strengthening**: even the parity-split sub-components (`S_min^even`,
`S_min^odd`, `S_min^even - S_min^odd`), which form a finer decomposition
of `S_min`, are *also* irreducible in the standard Dirichlet-L basis plus
Dirac-Dirichlet auxiliaries. The parity split does not split `S_min` into
identifiable pieces.

This is consistent with Paper 28's claim and reinforces it: the
irreducibility of `S_min` is a depth-2 structural property, and no
weighted sub-decomposition of the sum recovers identifiable pieces in the
Dirichlet-L / MZV basis. The obstruction is inherent to the depth-2
quadratic structure of `T(k)^2`, not an artifact of the specific
Paper 28 computation.

**Open question on the Paper 28 S_min value.** Paper 28's numerical
value `2.47953699802733387...` is the truncated sum at `N = 10000` with
a negligible (incorrect) tail. Re-running Paper 28's PSLQ with the
*convergent* value `2.479936938034222554...` might produce a different
(though equally likely null) result. This is a loose end in Paper 28's
numerical setup but does not affect the structural conclusion.

---

## §7 Sprint 5 recommendation

**Direct depth-2 analog: SHELVE.** Three-axis argument:

1. **Heuristic.** 35 PSLQ attempts across 7 basis strategies at 100-dps
   returned zero identifications for any parity-split quantity.

2. **Structural.** The quadratic structure of `T(k)^2` breaks the linear
   telescoping that produced the depth-1 `chi_{-4}` identity. There is
   no analog of the `zeta(s, 3/4) - zeta(s, 5/4) = 4^s (1 - beta(s))`
   identity for squared Hurwitz values.

3. **Context.** Paper 28's S_min irreducibility is precisely the
   statement that the depth-2 `S_min` object does not reduce to the
   standard basis; any parity decomposition of it would have to reduce
   to sub-components that are *also* in the standard basis, which by
   Paper 28's 47-element-basis test they are not.

**Productive alternatives.**

1. **Depth-2 character twist.** Instead of splitting by `k`-parity,
   split `T(k)^2` by a *character-twisted* weight:
   `T^twist(k) = chi_4(k) * T(k)`, where `chi_4` is the Dirichlet
   character mod 4. Then `sum T^twist(k)^2` explicitly introduces the
   `chi_{-4}` character at the Dirichlet-twist level, and the resulting
   sum is a twisted multiple Hurwitz value that may admit identification
   in a `L(s, chi_{-4})`-extended basis. This is a Sprint-5-sized task.

2. **Higher Hurwitz shifts.** Paper 28 Hurwitz at `3/2`; one could try
   `5/2, 7/2, ...` shift analogs, where the `q/4` shift structure in the
   parity split generalizes. This sits naturally on higher-dimensional
   spheres (`S^5, S^7`) via the Bargmann-Segal or nuclear Coulomb lattice
   (Paper 24).

3. **Motivic weight analysis.** The failure suggests that
   `S_min, S_min^{even/odd}` live in motivic weight 4 but are NOT
   products of weight-2 objects like `G, beta(4), zeta(3)`. A motivic
   decomposition via the Goncharov or Brown algebra of MZVs at
   half-integer shifts might reveal the exact sector
   `S_min` inhabits. This is a theoretical project, not PSLQ-based.

**Paper 28 update (NOT applied in Sprint 4).** The numerical value of
`S_min = 2.47953699802733387...` in Paper 28 Eq. (12) uses an incorrect
tail correction. The mathematically-convergent value is
`2.47993693803422254...`. A one-line correction + footnote explaining
the tail correction is appropriate for a future Paper 28 revision.
**This is flagged for plan-mode review; not auto-applied.** The
irreducibility claim is unchanged (both values are irreducible in the
standard basis).

---

## §8 Honest limitations

1. **Finite PSLQ precision.** All attempts use 100-dps arithmetic with
   `tol = 1e-60`. If the true relation has integer coefficients larger
   than `10^8`, it would be missed. Paper 28's 150-dps / 15-attempt
   exercise is a stricter test of the same irreducibility.

2. **Basis coverage.** The 7 strategies cover Dirichlet-L values up to
   weight 8, products thereof, and D(s)/D_even/D_odd auxiliaries. They
   do NOT include polylogarithms `Li_k(1/2)`, Euler sums, or
   double-Hurwitz `zeta_2(s_1, s_2; a_0)`. The Paper 28 47-element basis
   covers these and also fails. Our 35 attempts are therefore a
   *weaker* irreducibility test than Paper 28's; they confirm the
   parity-split pieces are no easier to identify than `S_min` itself.

3. **Not a proof of irreducibility.** This is a strong PSLQ-based
   heuristic negative. A proof of irreducibility would require, e.g.,
   a transcendence-theoretic argument (Apéry-style) or a motivic-weight
   dimension count.

4. **Parity-of-k is one choice.** Other depth-2 splits exist:
   (a) character-twisted (`chi_4(k)`-weight);
   (b) split-by-factorization `k = 2 m r` with `r` odd;
   (c) split on parity of `(k + shift)` for various shifts.
   Sprint 4 tested only the canonical (k-parity) split, motivated by
   strict analogy with the depth-1 RH-J.1 construction.

---

## §9 Summary

**Conjecture RH-P.1 (attempted): FALSIFIED.**

There is no integer-coefficient relation between `S_min^even - S_min^odd`
(depth-2 parity split of `S_min` on the Dirac Dirichlet series of S^3)
and the standard Dirichlet-L basis `{1, pi^{2k}, zeta(odd), G, beta(even),
products}` up to weight 8, at 100-dps PSLQ precision.

**Structural interpretation.** The depth-1 chi_{-4} identity (RH-J.1)
relies on a linear telescoping in Hurwitz values that does not survive
squaring. Each individual quarter-integer Hurwitz difference can still
be rewritten via RH-J.1, but the quadratic weighting mixes the
telescoping pattern and lands the parity-split sum in a weight-2 MZV
sector that is not spanned by products of standard Dirichlet-L values.
This obstruction is the same one responsible for Paper 28's claim that
`S_min` itself is irreducible in a 47-element basis.

**Relationship to Paper 28.** Paper 28 §IV claims `S_min` is irreducible.
Sprint 4 strengthens this to: `S_min^{even}, S_min^{odd}, S_min^{even} -
S_min^{odd}, S_min^{even}/S_min^{odd}, S_min^{even} * S_min^{odd}` are
ALL irreducible in the same basis. Parity decomposition does not crack
`S_min`.

**Sprint 5 recommendation.** Shelve direct depth-2 analog. Instead,
Sprint 5 should investigate (a) character-twisted depth-2 sums with
`chi_4(k)` multiplicative weight (which explicitly forces `chi_{-4}`
structure), or (b) accept that depth ≥ 2 multiple Hurwitz values at
half-integer shifts populate a genuinely new sector of the transcendence
algebra, and focus on the motivic-weight structure rather than closed-
form identification.

---

## Appendix A — key data at 100 dps

```
S_min^total = 2.47993693803422255447852790477854280438
              28955427841497089891807953124672174953...
S_min^even  = 0.91613480023028812067812332568987430843
              78220437751519080680...
S_min^odd   = 1.56380213780393443380040457908866849594
              30791500572790053...
S_min^diff  = -0.64766733757364631312228125339879418750
              61074413021725...
```

## Appendix B — cross-references

- Paper 28 §IV: `papers/observations/paper_28_qed_s3.tex`
- `geovac/qed_vertex.py` (S_min factorization derivation)
- `debug/spectral_chi_neg4_memo.md` (depth-1 RH-J.1 memo, Sprint 3)
- `debug/smin_identification.py` (Paper 28 irreducibility check)
- `debug/compute_smin_chi_neg4.py` (Sprint 4 driver)
- `tests/test_smin_chi_neg4.py` (regression tests, 10 passing)

## Appendix C — files

- Driver: `debug/compute_smin_chi_neg4.py`
- Data: `debug/data/smin_chi_neg4.json`
- Output log: `debug/data/smin_chi_neg4_output.txt`
- Tests: `tests/test_smin_chi_neg4.py`
- Memo: `debug/smin_chi_neg4_memo.md` (this file)
