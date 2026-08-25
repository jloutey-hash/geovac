# Sprint: four-center configuration operator -- a continuous four-body modulus (2026-08-25)

Canonical memo. Building the geometry-axis object flagged at the end of
`sprint_four_subspace_lambda_probe_memo.md`: does a FOUR-nuclei molecule's configuration
operator F = sum_i P_i carry a continuous four-body modulus (the D~4 four-subspace freedom)
that the three-center case (discrete Z2 conical intersections) lacked? YES.
Drivers `debug/four_center_config.py`, `debug/four_body_threshold.py`; test
`tests/test_paper32_four_center_modulus.py`.

## The object
Four planar centers, in-plane {s,px,pz} per center (rank 3) via the validated Slater-Koster
overlap block. G (12x12), whiten, four projectors, F = sum P_i. The genuine four-body invariant
is the **loop holonomy** W = Tr(S12 S23 S34 S41) = Tr(P1P2P3P4) -- the Wilson loop of the overlap
connection around the 4-cycle.

## Results
1. **Rank threshold (decisive, `four_body_threshold.py`):** local Jacobian-nullspace test of
   whether Tr(P1P2P3P4) is a function of pairwise+triple data:
   - rank 1 (four 1s): residual **2.6e-10 = PAIRWISE** (the 4x4 Gram IS the config; no four-body
     dof; Tr(P1P2P3P4) = product of pairwise overlaps, exact).
   - rank 2 / 3: residual **0.47 / 0.34 = GENUINE FOUR-BODY**. The four-body dof appears exactly
     at rank >= 2 (each center a sigma-pi pair) -- the rank-1->2 step of the tame four-subspace
     (D~4) problem.
2. **Gauge-invariant:** W is invariant under per-center O(r) gauge to machine precision (0, 3e-17,
   1e-16 at rank 1/2/3) -- a physical configuration invariant.
3. **Continuous modulus (molecular):** over a rhombus shape sweep (diagonal ratio t=0.5..1.8) W
   runs smoothly and monotonically **0.79 -> 0.042** (finite differences smooth) -- a genuine
   one-parameter continuous four-body modulus, the freedom the three-center Z2-CI case lacks.
4. **U(1) lift:** real orbitals -> W real (phase Z2); a small magnetic (Peierls) flux eps on the
   4-cycle -> arg W grows continuously from 0 (0.028/0.055/0.109/0.182 at eps=0.05/0.1/0.2/0.35,
   ~ linear; PSD preserved) -- the four-body Z2 -> U(1) lift.

## Interpretation / honest scope
The four-nuclei configuration DOES carry a genuine continuous four-body modulus -- the answer to
the lead's residual geometry-axis question is YES. But it is the **non-abelian loop holonomy**
(a matrix invariant), the *-representation analog of the tame D~4 cross-ratio, NOT a single scalar
cross-ratio: four ORTHOGONAL projections are *-WILD (>=3 projections wild, Halmos), so the tame
scalar cross-ratio (which is the GL/linear four-subspace modulus) does not literally apply; the
gauge-invariant holonomy is what survives. It is a GEOMETRY-axis invariant (function of nuclear
shape), **distinct from Paper 59's density-axis lambda** (confirming the earlier disambiguation:
different data, different axis; they are not the same number).

## Captured
- Paper 32 new `rem:four_center_modulus` (threshold, gauge-invariant holonomy, continuous modulus,
  U(1) lift, *-wild-vs-tame honesty, distinct-from-P59). Compiles clean (87 pp).
- Backing `tests/test_paper32_four_center_modulus.py` (3 fast: threshold rank1/rank2; gauge
  invariance; molecular modulus continuity + U(1) flux).
- CHANGELOG v5.1.0 + claim_test_matrix + memory.
