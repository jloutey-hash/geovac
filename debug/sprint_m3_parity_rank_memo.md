# Sprint: M3 period-map rank — vertex-parity resolution beats the rank-1 collapse
**Date:** 2026-07-08
**Context:** /ahha deep-dive follow-on. Grounding the cosmic-Galois injection
(Paper 56 `thm:injection_g4`) surfaced the **"sector-resolved rank-≥2 M3 period
map"** as the named open follow-on Paper 56 flags as *"of uncertain outcome /
not currently in the corpus."* This sprint supplies exactly that object — the
"per-(n,l) parity decomposition of D(s)" the paper names — and resolves it.

## Question
Paper 56 + `tests/test_paper56_injection_g4_periodmap.py` compute the M3-column
period-map Gram matrix and find **rank 1** (collinear) under the *only* per-sector
evaluation the corpus supplies: the parity-BLIND scalar η_{(n,l)}. Is rank-1
structurally forced, or an under-resolution artifact?

## Result (bit-exact; `debug/sprint_m3_parity_rank_probe.py`)

|                    | parity-BLIND | parity-RESOLVED |
|--------------------|:------------:|:---------------:|
| n_max=1 (N=2)      | rank 1       | rank 1          |
| n_max=2 (N=5)      | rank 1       | **rank 2**      |
| n_max=3 (N=9)      | rank 1       | **rank 2**      |
| n_max=4 (N=14)     | rank 1       | **rank 2**      |

- **Parity-blind** (single scalar η): rank 1 — reproduces Paper 56 / the periodmap test exactly.
- **Parity-resolved** (proven vertex parity, Paper 28 Thm 3): **rank 2** for every n_max ≥ 2, under both sector↔shell offset conventions (σ=(−1)^n and (−1)^{n−1}).
- **Rigor / anti-artifact check:** the two shift-classes are a *genuine* independent
  period direction — D_even(4)=π²/2−π⁴/24−4G+4β(4) and D_odd(4)=…+4G−4β(4) have
  **rank 2** in the motivic basis {π², π⁴, G=β(2), β(4)}. Not a manufactured sign.

## Mechanism (grounded, not invented)
The M3 content is the Dirac spectral zeta D(s)=2ζ(s−2,3/2)−½ζ(s,3/2). Paper 28
Thm 3 (χ₋₄ identity, PROVEN) splits it by CH-shell parity: **even shells → Hurwitz
shift 3/4, odd shells → 5/4**, and D_even − D_odd = 2^{s−1}(β(s)−β(s−2)) is pure
level-4 Dirichlet-L (Catalan G, β(4)). Each Coulomb sector (n,l) sits in its
shell's parity class, so it contributes at one shift-class. Period matrix rows =
η_{(n,l)} · (basis vector of its shift-class) ⇒ Gram rank = number of populated
shift-classes = **2** (for n_max ≥ 2; 1 at n_max=1 where only odd shells exist).

## Interpretation
- The rank-1 collapse in Paper 56 is a **parity-blind under-resolution artifact**.
- The true rank is **2 = the number of motivic LEVELS the vertex-parity χ₋₄
  mechanism accesses**: level-2 (un-restricted D(s)) ⊕ level-4 (χ₋₄ β-content).
- **Provably capped at 2:** there are only two parity classes, so this route
  *cannot* reach a faithful embedding (rank N) no matter how many sectors are
  added. GeoVac remains a **strict, now rank-2** subgroup — not rank-1, not faithful.
- This is the decisive answer to the named follow-on: the per-(n,l) parity
  decomposition exists, gives rank 2, and is capped — so Paper 56's implicit hope
  that "rank-≥2 ⇒ injectivity via Goncharov–Deligne" is unreachable by parity
  (rank 2 « N ⇒ not injective).

## Honest scope
- **Computed facts are bit-exact and solid:** parity-blind rank 1, parity-resolved
  rank 2, the two directions independent. Not in doubt.
- **One interpretive judgment (PI call):** whether the χ₋₄-graded trace
  Tr((−1)^N D e^{−tD²}) is admitted into *the* canonical Paper 56 period map. It
  is a legitimate M3-slot observable (same k=1 operator order; a *proven*
  GeoVac observable — it is what produces Catalan G in the QED vertex sums,
  Paper 28 Thm 3), but admitting it refines a load-bearing theorem, so it is
  the PI's adjudication, not the PM's.
- **Existing test not broken:** `test_paper56_injection_g4_periodmap.py` is about
  the parity-blind map (correctly rank 1); this is a different, richer map.
- **Big picture unchanged:** GeoVac is still a strict subgroup — "strict rank-2"
  instead of "strict rank-1." Depth-≥2 remains capped (JLO-Depth2, Reading A).
- **Open:** whether a finer-than-parity resolution beats rank 2 (suspect not —
  parity is the only grading χ₋₄ provides; 1/4-shift appears only via identities).

## Files
- `debug/sprint_m3_parity_rank_probe.py` — driver (exact sympy; prints table + rigor check).
- `tests/test_paper56_m3_parity_rank.py` — backing regression test.
- `papers/group3_foundations/paper_56_tannakian_substrate.tex` — added
  `rem:m3_parity_rank` (ADDITIVE; existing rank-1 statement of `thm:injection_g4`
  preserved). **Pending PI review of the χ₋₄-admission judgment before any
  headline rank-1→rank-2 change or propagation to the "uncertain outcome" pointers.**

## Decision
Banked for PI review. **HELD pending PI OK:** CHANGELOG entry, CLAUDE.md §2
one-liner, MEMORY index line, propagation to the other "uncertain outcome"
mentions (C4 proof; §sec:open_g4; §sec:open_na1), and any commit. The load-bearing
gate is the interpretive admission of the χ₋₄-graded observable into the period map.
