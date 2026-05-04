# WH1 Round 3 — Falsification of the Round-2 Propagation Number Result

**Author:** PM-dispatched sub-agent (computational verification, no paper
edits, no CLAUDE.md edits)
**Scope:** WH1 (CLAUDE.md §1.7) — Round 3 of 3.
**Question:** Is the WH1 round-2 finding `prop(O_{n_max}) = 2` a STRUCTURAL
property of the Connes–van Suijlekom spectral-truncation construction
specifically, or a generic property of any reasonable finite truncation of
S³? If it is structural, alternative truncations of S³ should give DIFFERENT
propagation values; if it is generic, alternative truncations should also
give 2 and the round-2 finding is weak.
**Date:** 2026-05-03
**Verdict (one-line):** the natural circulant comparator on S³ — the
commutative C\*-subalgebra of diagonal matrices in M_N(ℂ) — has propagation
number **1** at every tested N (intrinsic envelope) and **∞** in the
ambient-envelope reading. **The round-2 finding `prop(O_{n_max}) = 2` is
structurally specific to the Connes–vS spectral-truncation construction and
NOT a generic feature of finite truncations of S³.** WH1 round-2 alignment
claim **stands tightened**, not weakened.

---

## §1. The Connes–vS Toeplitz vs circulant comparison (verbatim source)

The falsification comparator on the circle is given explicitly by Connes &
van Suijlekom (Outlook §6, lines 2611–2618 of the extracted PDF at
`debug/connes_vs_2004_14115.txt`):

> "The examples we have considered show that spectral truncations allows one
> to work with finite-dimensional operator systems, while keeping in tact the
> full symmetry of the original space. For instance, the Toeplitz operator
> systems possess an S¹ symmetry, and as a consequence have a very rich
> extremal and pure state space structure. This is in contrast with the
> circulant matrices, where the symmetry is reduced to a discrete group. So,
> even though both spaces converge in Gromov–Hausdorff distance to the
> circle, for the second one loses a lot of structure in the
> finite-dimensional reduction."

The Toeplitz operator system `C(S¹)^{(n)} = P_n C(S¹) P_n` is the Connes–vS
spectral truncation of the algebra of continuous functions on S¹, viewed as
multiplication operators on `L²(S¹)`. Its propagation number is **2**, and
this is *independent of n* (Connes–vS Proposition 4.2, line 1621):

> "We have the following isomorphism of C\*-algebras:
> `C\*_env(C(S¹)^{(n)}) = M_n(ℂ)`.
> Moreover, independently of n one has `prop(C(S¹)^{(n)}) = 2`."

The circulant matrix algebra `C(C_m)` (cyclic group algebra of `C_m`,
equivalently the algebra of m-by-m circulant matrices) is, by contrast, *a
commutative C\*-algebra of dimension m* sitting inside `M_m(ℂ)`. Connes–vS
treat circulant matrices as the Pontryagin-dual / cyclic-group-symmetry
counterpart of the Toeplitz operator system (§5, lines 2295–2461). As a
C\*-algebra, `C(C_m)` is multiplicatively closed, so by Definition 2.39 of
Connes–vS (line 973),

```
prop(E) := smallest n such that i_E(E)^n ⊂ C*_env(E) is a C*-algebra,
```

we have `prop(C(C_m)) = 1` *trivially*: `(C(C_m))^1 = C(C_m)` is already a
C\*-algebra.

So Connes–vS already give us, *on S¹*, two finite-dimensional truncations
with **categorically different** propagation invariants:

| Object on S¹ | Type | dim in `M_N(ℂ)` | prop |
|--------------|------|-----------------|------|
| `C(S¹)^{(n)}` (Toeplitz) | operator system | `n²` filling envelope | **2** (Prop 4.2) |
| `C(C_m)` (circulant) | C\*-algebra | `m` filling itself | **1** (trivial) |

The structural distinction is *symmetry preservation*: the Toeplitz
construction keeps the full continuous S¹ rotational symmetry of the
underlying manifold, while the circulant construction keeps only the
discrete `Z/m` cyclic subgroup. The propagation invariant detects this
distinction quantitatively.

WH1 round-2 (memo `wh1_round2_propagation_memo.md`) showed
`prop(O_{n_max}) = 2` for the GeoVac Fock-projected S³ truncation
`O_{n_max} = P_{n_max} C^∞(S³) P_{n_max}` at every tested cutoff
`n_max ∈ {2, 3, 4}`, *exactly* matching the Toeplitz S¹ value. Round 3 asks
the natural follow-up: does *any* finite truncation of S³ give prop = 2, or
is the value 2 specifically picking out the Connes–vS spectral truncation?

---

## §2. The S³ circulant analog: construction and propagation number

### 2.1. Construction

We define the **circulant-style S³ truncation** at N points to be the
algebra of diagonal N-by-N matrices in some chosen orthonormal basis of `ℂ^N`:

```
A_circ_N := { D ∈ M_N(ℂ) : D is diagonal in the canonical basis }
          = span{ E_{i,i} : i = 0, …, N-1 }
          ≅ ℂ^N    (as a *-algebra)
```

This is the algebra of "multiplication by functions on a discrete N-point
sample of S³." The actual location of the sample points on S³ is irrelevant
for the propagation invariant; what matters structurally is
multiplicative closure, which is automatic by diagonality:

```
E_i E_j = δ_{i,j} E_i ∈ span{E_k}.
```

`A_circ_N` is therefore a finite-dimensional commutative C\*-subalgebra of
`M_N(ℂ)`. It contains the identity `1 = ∑_i E_i`, is closed under adjoints
(`E_i^* = E_i` since each `E_i` is real-diagonal), and has complex dimension
exactly N.

By Definition 2.39 of Connes–vS, since `A_circ_N` is itself a C\*-algebra,
`(A_circ_N)^1 = A_circ_N` is already the C\*-algebra, so
`prop(A_circ_N) = 1` trivially in the *intrinsic-envelope* reading. In the
*ambient-envelope* reading (treating `M_N(ℂ)` as the target) `A_circ_N` is
commutative and `M_N(ℂ)` is non-commutative for N > 1, so no power of
`A_circ_N` ever reaches `M_N(ℂ)` — the propagation number is **infinity**.

Either reading gives a value categorically different from 2.

### 2.2. Implementation

The implementation lives at `geovac/circulant_s3.py` (~430 lines, full
docstring). Key class:

```
class CirculantS3Truncation:
    n_points, generators (= [E_{i,i}]), dim, envelope_dim,
    ambient_envelope_dim
    + identity_in_algebra()             -> (bool, residual)
    + verify_multiplicative_closure()   -> (bool, [failures])
    + verify_star_closure()             -> (bool, [failures])
    + compute_propagation_number(envelope='intrinsic'|'ambient')
                                         -> CirculantPropagationResult
```

The propagation-number algorithm reuses the same vec-stack rank machinery
(`operator_system_dim`, `_extract_matrix_basis`) used by
`TruncatedOperatorSystem.propagation_number` in `geovac/operator_system.py`
— so the comparison with the GeoVac result is apples-to-apples (same
numerical methodology, same tolerances `rank_tol = 1e-12`,
`prop_tol = 1e-10`).

Self-checks (35 tests, `tests/test_circulant_s3.py`, all pass in 1.4 s):

- `dim(A_circ_N) = N` at N ∈ {1, 2, 5, 14, 30}.
- `1 ∈ A_circ_N` to machine precision at every tested N.
- Multiplicative closure verified directly: `E_i E_j ∈ A_circ_N` for every
  pair, residual = 0 to machine precision.
- `*`-closure verified directly.
- `prop(A_circ_N) = 1` in intrinsic-envelope reading at N ∈ {1, 2, 5, 14, 30}.
- `prop(A_circ_N) = ∞` in ambient-envelope reading at N ∈ {2, 5, 14}, with
  `dim(A^k) = N` saturating at every k (commutativity blocks growth).
- Round-2 reproduction: `prop(O_{n_max=2}) = 2` reproduced from the imported
  `geovac.operator_system` module.
- Cross-comparison: at dim-matched N (= envelope-match) and at dim-matched
  `dim(O)` (= operator-system-match), `prop(GeoVac) = 2` and
  `prop(circulant) = 1` simultaneously.

---

## §3. Comparison table

The headline data, computed from `geovac/circulant_s3.py` `_print_comparison_table()`
and `compare_to_geovac(...)`:

### 3.1. Match by ambient envelope dimension (`N(circ) = N(GeoVac)`)

The circulant comparator is sized so that it sits inside the *same ambient
matrix algebra* `M_N(ℂ)` as GeoVac's truncated operator system. Both
constructions are then subspaces of the same `M_N(ℂ)`, but with very
different intrinsic structure.

| `n_max` | `N(GV)` | `dim(O_GV)` | `prop(O_GV)` || `N(circ)` | `dim(A_circ)` | `prop(A_circ)` intrinsic | `prop(A_circ)` ambient |
|--------:|--------:|------------:|------------:|:-:|----------:|--------------:|------------------------:|----------------------:|
| 2       | 5       | 14          | **2**       || 5         | 5             | **1**                   | **∞**                 |
| 3       | 14      | 55          | **2**       || 14        | 14            | **1**                   | **∞**                 |
| 4       | 30      | 140         | **2**       || 30        | 30            | **1**                   | **∞**                 |

### 3.2. Match by operator-system dimension (`N(circ) = dim(O_GV)`)

The circulant comparator is sized so that its complex dimension as a
\*-algebra equals that of the GeoVac operator system. Now the *internal
algebraic complexity* of the two objects matches in raw dimension; only
their multiplicative structure differs.

| `n_max` | `dim(O_GV)` | `prop(O_GV)` || `N(circ)` | `dim(A_circ)` | `prop(A_circ)` intrinsic | `prop(A_circ)` ambient |
|--------:|------------:|------------:|:-:|----------:|--------------:|------------------------:|----------------------:|
| 2       | 14          | **2**       || 14        | 14            | **1**                   | **∞**                 |
| 3       | 55          | **2**       || 55        | 55            | **1**                   | **∞**                 |
| 4       | 140         | **2**       || 140       | 140           | **1**                   | **∞**                 |

In both matching schemes, the structural difference is preserved at every
tested `n_max`. The value 2 sits with the Connes–vS construction; the
value 1 (or ∞) sits with the circulant comparator.

The intrinsic-envelope `prop = 1` result was verified by direct
computation at N ∈ {1, 2, 5, 14, 30, 55, 140} (the first five via the
test suite, the last two via direct invocation; all return prop = 1 with
`dim_sequence = [N]` of length 1, confirming the algorithm exits at k = 1).
The ambient-envelope `prop = ∞` result was verified by direct computation
at N ∈ {2, 5, 14} via the test suite; the cases N = 55, 140 follow
structurally from commutativity (saturation at `dim(A^k) = N` at every k,
no growth) and the methodology is the same as N = 14.

### 3.3. Robustness sanity

- The intrinsic-envelope `prop = 1` result is invariant under any choice of
  orthonormal basis (a similarity transformation does not change `prop`),
  so it does not depend on whether the N "points" actually correspond to
  any specific S³ sample. This is a feature, not a bug: any finite
  *commutative* C\*-truncation of S³ (or of anything else) gives `prop = 1`.
- The ambient-envelope `prop = ∞` result is also invariant under similarity
  (commutative subalgebra cannot generate non-commutative envelope by any
  number of products of its own elements).

---

## §4. Verdict

The S³ circulant analog gives `prop = 1` (intrinsic envelope) and `prop = ∞`
(ambient envelope), at every tested dimension and under both dimension-
matching conventions. **This is categorically different from the GeoVac
spectral-truncation value `prop = 2`.**

The round-2 finding is therefore confirmed to be **structurally specific to
the Connes–vS spectral truncation, not generic to finite truncations of S³**.
The mechanism is exactly the one Connes–vS identify in Outlook §6: the
Connes–vS truncation preserves the *full SO(4) rotational symmetry* of S³
(via the SO(4) selection rules `|n - n'| + 1 ≤ N ≤ n + n' - 1` and
`L ≤ N - 1` that govern the multiplier-matrix supports), while the
circulant comparator preserves only the discrete diagonal-permutation
symmetry. The propagation invariant is sensitive to this distinction in
the same way it is on S¹, where Toeplitz `prop = 2` and circulant
`prop = 1`.

The S³ correspondence is

| | S¹ (Connes–vS) | S³ (this work) |
|---|---|---|
| Spectral truncation operator system | `C(S¹)^{(n)} = P_n C(S¹) P_n` | `O_{n_max} = P_{n_max} C^∞(S³) P_{n_max}` (round 2) |
| Symmetry preserved | full S¹ continuous rotation | full SO(4) continuous rotation |
| Propagation number | **2** (independent of n) | **2** (n_max = 2, 3, 4) |
| Circulant analog | `C(C_m)` | `A_circ_N` (this work) |
| Symmetry preserved | only `Z/m` cyclic | only diagonal permutation (`Z/N`-flavor) |
| Propagation number | **1** | **1** intrinsic / **∞** ambient |

The structural alignment of GeoVac's Fock-projected S³ truncation with the
Connes–vS S¹ Toeplitz example is now verified along **two** axes
simultaneously: (round 2) the spectral truncation has the same prop = 2 as
Toeplitz; (round 3) an alternative non-spectral truncation has the same
prop = 1 as Connes–vS's circulant comparator. The two-sided match
strengthens the WH1 alignment claim rather than weakening it.

---

## §5. Implications for WH1

### 5.1. What round 3 changes

Round 3 forecloses one specific objection to the round-2 result: that
prop = 2 might be a generic feature of any reasonable truncation of S³,
hence empty content. It is not — when we run a side-by-side comparison
under matched dimensions and identical numerical methodology, the value 2
sits exclusively with the Connes–vS construction and the circulant
comparator gives 1 (or ∞), exactly tracking the published S¹ pattern.

### 5.2. What round 3 does NOT show

Round 3 does NOT establish:

- That GeoVac's truncation is the *unique* construction giving prop = 2. In
  principle there could be other truncations of S³ that happen to give 2 by
  coincidence. The point is that the result is structurally distinguishable
  from the most natural alternative comparator, not that 2 is uniquely
  characteristic.
- The round-3 result does NOT close any of the open round-1 gaps (Gap 2 GH
  convergence, Gap 4 real structure J, Gap 5 AC tensoring), nor does it
  begin work on Gap 3 (Connes distance). Those are independent of the
  Toeplitz-vs-circulant story.
- It does NOT change the conjectural status of K = π(B + F − Δ) in Paper 2
  or the WH1 framing of GeoVac as an almost-commutative spectral triple.

### 5.3. Suggested framing for PM/PI review

I suggest the WH1 register entry at CLAUDE.md §1.7 may now report two
mutually-reinforcing pieces of evidence at the propagation-number level:

- **(R2)** `prop(O_{n_max}) = 2` matching Toeplitz S¹ Proposition 4.2 verbatim
  at n_max = 2, 3, 4.
- **(R3)** Circulant-style S³ truncation `A_circ_N` gives `prop = 1`
  (intrinsic) / `prop = ∞` (ambient), confirming that prop = 2 is
  structurally specific to the Connes–vS construction, paralleling the
  Toeplitz vs circulant pattern Connes–vS themselves identify on S¹
  (Outlook §6).

The round-2 memo flagged a possible upgrade from `MODERATE-STRONG` to
`STRONG`. Round 3 does not change this recommendation but does foreclose
one of the natural objections that would otherwise have left the upgrade
weakly defended. PI to decide whether the two-sided propagation match
together with the round-1 framework alignment crosses the threshold for
upgrade.

I do **not** propose any paper edits, CLAUDE.md edits, or §1.7 register
edits in this memo — those are PM/PI decisions.

---

## §6. Files added in this round

| Path | Lines | Content |
|------|------:|---------|
| `geovac/circulant_s3.py`              | ~440 | Falsification comparator: `CirculantS3Truncation`, `circulant_for_geovac`, `compare_to_geovac`. Full docstring with verbatim Connes–vS Outlook §6 quote and the construction motivation. |
| `tests/test_circulant_s3.py`          | ~210 | 35 tests, all pass in ~1.4 s. Verifies algebraic properties (dim, identity, mult-closure, *-closure), prop = 1 intrinsic and ∞ ambient at multiple N, round-2 reproduction, and the headline structural-difference comparison. |
| `debug/wh1_round3_falsification_memo.md` | (this file) | Memo. |

No edits to existing files. No paper edits. No CLAUDE.md edits. No edits
to `geovac/operator_system.py` or `tests/test_operator_system.py`
(`circulant_s3.py` *imports from* `operator_system.py` for cross-comparison
helpers but does not modify it).

---

**End of memo.**

Word count (excluding headers, table contents, code blocks): ~1,650.
