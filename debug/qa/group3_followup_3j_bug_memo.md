# Group3 QA follow-up — Paper 24 / Bargmann–Segal rank-1 Wigner-3j factor-2 bug

**Date:** 2026-06-15 · **Scope:** analysis only (no code/test/paper edits) · **Trigger:** `/qa group3` surfaced a wrong value in `wigner3j_rank1_squared`.

## 1. The bug (confirmed)

`geovac/nuclear/bargmann_graph.py :: wigner3j_rank1_squared` returns **exactly half** the true squared rank-1 3j symbol in **every** admissible case (verified against `sympy.physics.wigner.wigner_3j` for all `l ≤ 5`, both `l'=l±1`, all `m,q`: 0 cases agree, all off by ratio true/code = 2).

Root cause: a substituted denominator factor.
- **`l'=l+1` branch** uses `(2l+1)(2l+2)(2l+3)`; the Edmonds form is `(2l+1)(l+1)(2l+3)`. `(2l+2)=2(l+1)` is twice `(l+1)` → factor ½.
- **`l'=l-1` branch** uses `(2l-1)(2l)(2l+1)`; the Edmonds form is `(2l-1)·l·(2l+1)`. `(2l)=2·l` is twice `l` → factor ½.

Example: `wigner3j_rank1_squared(1,2,0,0)` returns `1/15`; true `|(1 1 2; 0 0 0)|² = 2/15`.

The test `tests/test_bargmann_graph.py::test_wigner3j_rank1_specific_value` **hard-codes the wrong value** (`assert v == 1/15`) — it confirms the bug rather than catching it. The test docstring even shows the buggy derivation `2·2/(3·4·5)=1/15`, where the correct denominator is `(3)(2)(5)=30 → 2/15`.

**This is a genuine pre-existing bug, NOT a planted QA seed.** `debug/qa/group3_seed_key.json` lists seeds S1/S3/S5/S8/S9 (S3 is a *different* test, `test_ho_spectrum_matches_exact`). The 3j defect is independent — the QA process surfaced real, un-planted dirt.

## 2. Correct closed form (Edmonds §5.4 / Varshalovich §8.5)

For `T^1_q` between `|l m⟩` and `|l' m'⟩`, `q = m'−m`, the corrected `|3j|²`:

**`l'=l+1`** (denominator `(2l+1)(l+1)(2l+3)`):
- `q=0`:  `(l+1−m)(l+1+m) / [(2l+1)(l+1)(2l+3)]`
- `q=±1`: `(l±m+1)(l±m+2) / [2(2l+1)(l+1)(2l+3)]`

**`l'=l−1`** (denominator `(2l−1)·l·(2l+1)`):
- `q=0`:  `(l−m)(l+m) / [(2l−1)·l·(2l+1)]`
- `q=±1`: numerators unchanged, denominator `2(2l−1)·l·(2l+1)`

Verified: these match `wigner_3j` for all `l ≤ 5` (0 mismatches).

### Minimal code fix (6 token-level edits in `bargmann_graph.py`)
- `l'=l+1` branch: in the three `Fraction(...)` denominators **and** the shared `den` on line 197, replace `(2 * l + 2)` → `(l + 1)`.
  - lines **190, 195, 197** (and the docstring math on 185/188/193 if cleaned).
- `l'=l-1` branch: in the three denominators, replace `(2 * l)` → `l`.
  - lines **205, 210, 213**.

That is the entire correction. No structural/algorithmic change.

## 3. Blast radius — what changes vs. what is invariant

The bug is a **uniform global ×2** of the angular factor → every adjacency weight scales by exactly 2 (verified: `A_fixed == 2·A_buggy` bit-exact at `N_max=5`). Consequences split cleanly:

### Consumers
| Consumer | Uses | Affected? |
|:--|:--|:--|
| `bargmann_graph.build_bargmann_graph` (adjacency) | weighted edges | values ×2 |
| `geovac/su3_wilson_s5.py :: bargmann_adjacency_dense` | **0/1 connectivity only** (`if w != 0`) | **NO** |
| `tests/test_ihara_zeta.py` (S⁵ Ramanujan) | `(A != 0).astype(int)` — **support** | **NO** |
| `tests/test_su3_wilson_s5.py` (plaquettes, bipartite, β₁) | 0/1 adjacency | **NO** |
| `geovac/nuclear/ho_two_fermion.py` (→ Paper 27 zero-entropy) | imports `enumerate_nodes` **only**; CG via separate correct `_orbital_cg` in `nuclear_hamiltonian.py` | **NO** |
| `debug/compute_s5_k_nmax_truncated.py` | squared A as κ·A perturbation; its own docstring says verdict is weight-convention-independent (parity argument is on diagonal `Tr(L^j)`) | debug-only; verdict invariant |

### Paper 24 claims (verified by rebuild with fixed 3j at `N_max=5`)
| Claim | Status under fix |
|:--|:--|
| 56 nodes | **INVARIANT** (node enumeration, no 3j) |
| 165 edges | **INVARIANT** (support of A) |
| π-free certificate (`irrationals=[]`, all rational) | **INVARIANT** (2× a rational is rational) |
| HO spectrum `ℏω(N+3/2)` + degeneracies | **INVARIANT** (lives in diagonal, not A) |
| Selection rules ΔN=±1, Δl=±1 | **INVARIANT** (support) |
| Hermiticity of A | **INVARIANT** |
| Graph Laplacian PSD, 1 zero-mode (connected) | **INVARIANT** (L scales ×2; PSD & null-space preserved) |
| β₁ = E−V+c = 110, incidence B∈ℤ^{56×165}, Hodge decomp | **INVARIANT** (combinatorial, from 0/1 support) |
| Paper 27 HO zero-entropy corollary | **INVARIANT** (Moshinsky–Talmi quanta conservation; doesn't use these weights) |
| SU(3) Wilson β₁ / "no SU(3) Wilson" (Paper 30 / §V) | **INVARIANT** (uses connectivity) |

**The Laplacian *eigenvalues* and adjacency *entry values* change (×2)** — but the paper publishes **no specific edge-weight or Laplacian-eigenvalue number**, and explicitly states (§ "Reduced content", ~line 985–989) that the node Laplacian `L₀ = D−A` "plays no role in reconstructing the spectrum" (HO energies live in the diagonal). So the scaling is physically inert for every stated result.

**Verdict on the first-pass reviewer's "headline claims survive global rescaling":** CONFIRMED and made precise. Every published number is either support-derived (count/selection/β₁/Hodge), diagonal-derived (HO spectrum), or rationality-only (π-free) — all invariant under the ×2. No published number changes.

### Why the reviewer missed it
`review_paper24.md` (#6, #60) states it "verified independently" the rank-1 3j closed forms — but it verified **rationality** (true for both `1/15` and `2/15`), not the numeric values. The π-free framing made rationality the load-bearing property, and rationality is exactly what's invariant under the bug. Lesson for the matrix: a "rational/π-free" test must not be mistaken for a "correct-value" test.

## 4. Tests that break or need updating when corrected

Only **one** assertion breaks:
- `tests/test_bargmann_graph.py::test_wigner3j_rank1_specific_value` — `assert v == Fraction(1, 15)` must become `Fraction(2, 15)`; the docstring's worked denominator (`3·4·5`→`3·2·5`) should be corrected too.

All other tests in `test_bargmann_graph.py` (selection rules, hermiticity, rationality, π-free, HO spectrum, Laplacian PSD), and **all** of `test_su3_wilson_s5.py` and `test_ihara_zeta.py`, pass unchanged (confirmed: 17 + 54 passed / 1 skipped on current tree; the 73 non-target tests are weight-invariant by construction).

## 5. Second, separate discrepancy (documentation, not a number)

Paper 24 Eq. (`eq:edge-weight`) writes the angular factor as the **Clebsch–Gordan square** `|C^{l'm'}_{lm;1q}|²`, but the code computes the **3j square**. These differ by `(2l'+1)` (CG² = (2l'+1)·3j²). This is *independent* of the factor-2 bug and also touches **no published number** (no edge-weight values are printed; the π-free/structural claims are magnitude-blind). Worth a one-line reconciliation: either relabel the paper factor as the 3j (matching code) or note that the code uses the proportional 3j convention. The physical dipole reduced matrix element is convention-dependent up to the `(2l'+1)` normalization and an overall constant; for a discrete graph regulator the proportional choice is harmless, but prose and code should agree.

## 6. Recommendation

1. **Fix the code** (6 token edits, §2): `(2l+2)→(l+1)` in the `l'=l+1` denominators (lines 190, 195, 197) and `(2l)→l` in the `l'=l-1` denominators (lines 205, 210, 213). Clean the in-line docstring formulas to the Edmonds form.
2. **Fix the test assertion**: `test_wigner3j_rank1_specific_value` → `Fraction(2, 15)` (and docstring). Strongly consider adding a 3-4 case cross-check against `sympy.physics.wigner.wigner_3j` so a value bug can never again hide behind a rationality test.
3. **No Paper 24 numerical edit required** — 56 / 165 / π-free / spectrum / β₁ all stand verbatim. Optionally add a one-line note (or fix the CG-vs-3j label, §5) for prose–code agreement; this is cosmetic, not a correctness fix.
4. **Matrix/QA note**: log that the π-free certificate test family proves *rationality*, not *correctness of values* — the gap that let the bug survive ~since Sprint NK. A correct-value cross-check closes it.

**One-line risk note:** Low risk — the fix is a bit-exact global ×2 of inert edge weights; it changes only the one wrong unit-test value and the (non-published) adjacency/Laplacian magnitudes, leaving every published Paper 24/27/30 number invariant. The only way to *introduce* risk would be to leave the code wrong while a future result starts depending on the actual weight magnitude (e.g. a κ·A spectral claim), so correct it now while the blast radius is zero.
