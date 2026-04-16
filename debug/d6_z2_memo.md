# D6 Memo: Three Z_2 Symmetries on the S^3 Graph

**Track:** D6-A/B (Dirac-on-S^3 program)
**Date:** April 2026
**Status:** Complete

## Summary

Three candidate Z_2 symmetries were analyzed on the n_max=3 S^3 graph (14 scalar nodes, 13 edges). Two act only on the spinor-extended space; one acts on the scalar graph.

## Results

### sigma (m -> -m): Azimuthal reflection

- **Acts on:** scalar graph (14 nodes)
- **Is involution:** Yes (P^2 = I)
- **Is graph automorphism:** Yes (commutes with adjacency matrix A, degree matrix D, and Laplacian L)
- **Fock-compatible:** Yes (automatic from graph automorphism)
- **Fixed points:** 6 (all m=0 states: one per (n,l) pair)
- **2-cycles:** 4 (the m/-m pairs: (2,1,+/-1), (3,1,+/-1), (3,2,+/-1), (3,2,+/-2))
- **Eigenvalue decomposition:** 10 states in the +1 sector, 4 in the -1 sector

sigma is the only Z_2 that acts on the scalar graph. It is a genuine graph automorphism: the adjacency structure is symmetric under m-reflection because the angular transition rules (m <-> m+/-1) and radial rules (n <-> n+/-1, same l,m) are both m-sign-symmetric.

### P (kappa -> -kappa): Dirac kappa-parity

- **Acts on:** Dirac labels (n, kappa, two_m_j) — 28 labels at n_max=3
- **Maps label set to itself:** NO (12 labels unmapped)
- **Is valid permutation:** No
- **Is involution:** No (not a valid map on the finite label set)

**Critical finding:** kappa -> -kappa does NOT close on the Dirac label set at finite n_max. The problem: for kappa = -(l+1) with l = n_fock - 1 (the maximum l at each n_fock), the target -kappa = l+1 > 0 would require l' = l+1 = n_fock, which violates l < n_fock. Specifically, 12 states are unmapped:
- n=1, kappa=-1 (l=0): target kappa=+1 requires l=1 >= n_fock=1
- n=2, kappa=-2 (l=1): target kappa=+2 requires l=2 >= n_fock=2
- n=3, kappa=-3 (l=2): target kappa=+3 requires l=3 >= n_fock=3

**Physical interpretation:** kappa -> -kappa is the accidental degeneracy map, pairing states with the same (n, j) but different l (differing by 1). At the boundary of each n-shell, the highest-l states (l = n-1) have no partner because l+1 = n is not available. This is a finite-truncation artifact: in the n -> infinity limit, every state except the boundary gets a partner.

For the 16 states that DO have valid partners, kappa -> -kappa correctly identifies the Dirac accidental degeneracy pairs:
- 2s_1/2 (kappa=-1) <-> 2p_1/2 (kappa=+1)
- 3s_1/2 (kappa=-1) <-> 3p_1/2 (kappa=+1)
- 3p_3/2 (kappa=-2) <-> 3d_3/2 (kappa=+2)

All 6 pairs through n=4 confirmed: the accidental degeneracy map IS exactly kappa -> -kappa (simple negation), not kappa -> -kappa - sgn(kappa).

### C (chirality gamma^5 flip): Weyl sector swap

- **Acts on scalar graph:** No (chirality is spinor-only)
- **Acts on Dirac-on-S^3:** Yes, swapping +/- Weyl sectors
- **At n_max=3 (CH n=0,1,2):** 40 total Dirac states = 20 + 20 Weyl (exactly the Phase 4H SM-D split: Delta^{-1} = 40)

Chirality is structurally different from sigma and P: it doubles the state space and swaps the two copies.

## Key Question: Relation Between sigma and P

**sigma and P are genuinely independent Z_2 symmetries.** They cannot be identified through any natural embedding.

| Property | sigma | P |
|----------|-------|---|
| Space | scalar (n,l,m) | Dirac (n,kappa,m_j) |
| Action | m -> -m | kappa -> -kappa |
| What's preserved | n, l | n, j (= \|kappa\| - 1/2), m_j |
| What changes | sign of m | sign of kappa (hence l changes by +/-1) |
| Physical meaning | Time reversal (on magnetic QN) | Spin-orbit alignment flip |
| Closes on finite basis | Yes (always) | No (boundary states unmapped) |
| Graph automorphism | Yes | N/A (no scalar action) |

Under the natural scalar-to-Dirac lift (n,l,m) -> {(n, -(l+1), m_j), (n, +l, m_j)}:
- sigma lifts to m_j -> -m_j (time reversal on spinors)
- P has no scalar pre-image (it mixes different l values)

## Group Structure

On the full Dirac-on-S^3 space (in the n -> infinity limit where P closes):

    G = Z_2(sigma) x Z_2(P) x Z_2(C)

- sigma and P commute (act on different quantum numbers: m_j vs kappa)
- C acts on a separate sector index (chirality)
- On the scalar graph, only sigma acts; P and C are invisible

On the finite n_max=3 graph:
- sigma generates Z_2 on the 14-node scalar graph (the ONLY symmetry)
- P is broken by truncation (does not close)
- C is invisible (spinor-only)

## Data

- Script: `debug/d6_z2_enumeration.py`
- JSON: `debug/data/d6_z2_enumeration.json`
