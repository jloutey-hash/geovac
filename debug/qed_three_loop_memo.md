# Three-Loop QED on S³: Preliminary Results

**Date:** April 17, 2026
**Status:** INFRASTRUCTURE BUILT — convergence insufficient for PSLQ identification

## Diagram Topology

Iterated sunset ("chain"): three internal electron lines connected by two
internal photon lines:

    [n1] --q1-- [n2] --q2-- [n3]

Vertex 1: (n1, n2, q1) — SO(4) selection rule (triangle + parity)
Vertex 2: (n2, n3, q2) — SO(4) selection rule (triangle + parity)

Five-fold sum: O(n_max^5) with selection rules cutting ~90%.

## Spectrum Data

- Dirac: |λ_n| = n + 3/2, g_n = 2(n+1)(n+2)
- Hodge-1: μ_q = q(q+2), d_q^T = q(q+2), q ≥ 1

## Convergence Table (CG-weighted)

| n_max | S(n_max) | quintuples | time |
|-------|----------|------------|------|
| 3 | 0.2908 | 45 | 0.0s |
| 5 | 1.0616 | 566 | 0.0s |
| 7 | 2.0403 | 2,891 | 0.1s |
| 10 | 3.5721 | 16,257 | 0.2s |
| 12 | 4.5503 | 39,402 | 0.4s |
| 15 | 5.9083 | 116,823 | 1.0s |
| 20 | est ~8+ | ~480K | ~5s |

Unrestricted D(4)³ = 5.3794 (purely π^{even}).

## Convergence Table (vertex-restricted, no CG weights)

| n_max | S(n_max) | quintuples | time |
|-------|----------|------------|------|
| 10 | 1.7090 | 17,017 | 0.1s |
| 15 | 2.5453 | 119,288 | 1.2s |
| 20 | 3.2091 | 482,734 | 4.6s |
| 25 | 3.7451 | 1,437,605 | 20.8s |
| 30 | 4.1874 | 3,519,151 | 51.1s |

## Key Observations

1. **CG-weighted sum EXCEEDS unrestricted at n_max ≥ 15**: S_cg(15) = 5.91 > D(4)³ = 5.38.
   The CG channel count W ∈ {0,1,2} with W=2 for many triples amplifies the sum.
   This is structurally different from two loops where CG weights always reduce the sum.

2. **Vertex-restricted (no CG) converges slowly**: at n_max=30, still only 78% of D(4)³.
   The vertex parity constraint (n1+n2+q1 odd AND n2+n3+q2 odd) is more restrictive
   at three loops — it forbids more term combinations than at two loops.

3. **PSLQ cannot identify at current convergence**: even vertex-restricted at n_max=15
   fails PSLQ. The convergence deficit (~3 digits) is too large for 80-dps PSLQ.
   Need n_max ≈ 50-100 for clean identification (computationally expensive at O(n^5)).

## Structural Prediction

At k loops with vertex restrictions, expect depth-k MZVs:
- k=1: depth-1 → π^{even} ✓ (T9 theorem)
- k=2: depth-2 → S_min (irreducible, confirmed this session)
- k=3: depth-3 → PREDICTED but unverified (convergence insufficient)

The nested parity filtration at two vertices creates even-even / even-odd /
odd-even / odd-odd sub-sums, each involving quarter-integer Hurwitz shifts.
This is the mechanism expected to produce depth-3 structure.

## What Would Be Needed

To get clean PSLQ identification:
- Vertex-restricted: n_max ≈ 50+ (estimated ~10 min at O(n^5))
- CG-weighted: n_max ≈ 40+ (same)
- Alternative: Hurwitz closed-form decomposition (analytically subtract known terms)

## Key Files

- `geovac/qed_three_loop.py` — implementation (671 lines)
- `tests/test_qed_three_loop.py` — 18 tests (15 fast, 3 slow)
- This memo

## Verdict

INFRASTRUCTURE COMPLETE, RESULT PENDING. The three-loop module works correctly
and the convergence data is consistent. The depth-3 prediction is structurally
motivated but computationally unverified. Higher n_max runs or analytical
decomposition needed for definitive identification.
