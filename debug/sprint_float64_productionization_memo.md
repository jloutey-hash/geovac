# #2 float64 productionization of the Route C ERI engine — progress + handoff

**Date:** 2026-09-22 (v5.15.18+). **PI-directed (3→2→1 plan, item 2).** Goal: turn the
Route C / core-enrichment engine from minutes/point (dps=60 mpf) into seconds, so it becomes
a *sweep* — and can break the dense M=16 wall for the last radial/valence mHa, and make the
eventual r₁₂ geminal build tractable.

## Finding 1 — dps-lowering is a DEAD END (measured)
`debug/data/dps_timing.log` (M=10 confirm baseline, bond(2,1) no π): the LiH energy is
**bit-identical from dps=60 down to dps=20** (dE = 0.00 µHa) while the time barely moves
(234s → 221s, ~6%). The engine downcasts mpf→float64 at the end, so dps was never
load-bearing for accuracy; and the cost is the **COUNT of mpf operations** (Python/mpmath
object overhead, ~flat in dps at small dps), not the arithmetic precision. So lowering dps
buys nothing. **The only real speedup is the float64 ASSEMBLY** (each mpf op → a hardware
float64 op, ~100× per op).

## Finding 2 — the float64 design already exists (single-exponent)
`geovac/neumann_vee_general_m.build_Xtab` productionizes the single-exponent general-m
Neumann X-table exactly the right way: mpf B-tables/moments at a small guarded dps
(`_DPS`=30 + the `_seed_guard(s)=8s+24` digits the unstable Q_l recurrence needs), then
**downcast to float64** and run the X assembly in float64, with a clean float64 IBP-tail
correction `ngm._corr`. C4's engine ignores this and does the WHOLE assembly in mpf.

## Deliverable — `build_Xtab_s_f` (float64 twin of C4's `build_Xtab_s`), VALIDATED
`debug/prolate_float_eri.py`. C4's `build_Xtab_s(m,s1,s2,l_hi,p_max,α1,α2)` is the mixed
generalization the ERI actually calls (independent weights s1,s2 AND rates c1,c2). The
float64 twin is a direct hybrid port: mpf seed tables at `SEED_DPS` (with the ngm guard) →
downcast → the two-ordering assembly in float64. **Key reuse:** `ngm._corr` works for the
two-rate case UNCHANGED — pass the inner rate c and the c1+c2 B-table, exactly as the mpf
`build_Xtab_s` passes them to `pr._corr_mp`.

**Gate (`python debug/prolate_float_eri.py`):** `build_Xtab_s_f` vs `float(mpf build_Xtab_s)`
across (m,s1,s2,c1,c2) incl. s1≠s2, big exponent ratio (core-valence 4.03 vs 1.0), m=0/1/2:
- **low-l (l ≤ m+4, energy-dominant): rel 1e-13 … 1e-10** — bit-exact, mechanism CORRECT.
- high-l: rel 1e-5 … 5e-2 — the DOCUMENTED Neumann-prefactor-suppressed float64 wall
  (ngm: "l=10 block degrades to ~4e-9, prefactor-suppressed in energy"); harmless in energy.
So raw X-table bit-matching is the wrong gate (dominated by the harmless high-l blocks); the
energy-level check is the right one — the low-l blocks that carry the energy are bit-exact.

## NEXT STEPS (the clear continuation)
1. **Float64 `eri_general` / `build_eri_tensor_m`.** Assemble the ERI from `build_Xtab_s_f` +
   the eta moments (`_Ymom_m`) + the Neumann prefactor + the (2π)² φ factor, all in float64
   (currently mpf in `prolate_allelectron_c4.eri_general`). An `engine='float'` switch on
   `eri_general`/`build_eri_tensor_m` is the cleanest (keep the mpf path as the reference).
2. **Energy-level validation:** float64 LiH energy vs the banked mpf **−8.02905** (M=16
   core-enriched) to ~1e-6, and vs H2 (the make-or-break control). This is the gate that
   matters (the high-l raw degradation washes out under the prefactor).
3. **Speedup measurement:** re-time an M=16 LiH point (mpf ~1079s → float64 target ~seconds).
4. **Then it enables:** fast core-enrichment SWEEPS (find the best M≤16 allocation), breaking
   the dense M=16 wall for the last radial/valence mHa (toward ~−8.04), and a tractable r₁₂
   geminal build (item 1 of the 3→2→1 plan).

## WIRING DONE + ENERGY-VALIDATED — but the speedup is only ~2× (the real walls are elsewhere)

`debug/prolate_float_eri.py` now has the full float64 path: `build_Xtab_s_f` +
`eri_general_f` + `build_eri_tensor_m_f`, monkeypatched into `L.assemble_rebased` and run
against the banked mpf energies (`validate` mode; `debug/data/float_validate.log`).

**Accuracy — PASS (energy-exact):** float64 ERI reproduces the banked mpf LiH energies to
**dE = +0.004 … +0.005 mHa (~5 µHa)** at M=6/10/16. The documented high-l float64 X-table
degradation IS prefactor-suppressed at the energy level, confirmed. The port is correct.

**Speedup — only ~2×** (M=6 61s vs 147s; M=10 101s vs 237s; M=16 636s vs 1079s). The
float-assembly floor was necessary but NOT the bottleneck. Two walls remain, and this is the
load-bearing finding:
1. **The mpf SEED tables are the ERI bottleneck (not the assembly).** A single X-table build
   is only 2.5× faster in float (0.27s→0.11s) because `ngm._B_table` (with the `_seed_guard`
   digits) dominates it, and those seeds are INTRINSICALLY mpf — the guard digits absorb the
   unstable forward Q_l recurrence's cancellation; drop them and the recurrence gives garbage
   at large l (documented: H2 (5,5)+δ at dps=40 flipped 99.767%→−220). So flooring the
   assembly caps at ~2.5× on the ERI.
2. **The dense FCI Python loop is the co-dominant M=16 wall.** `fci_energy` builds an nd×nd
   dense H with a pure-Python double loop over determinant pairs; at M=16 nd=14400 →
   ~10⁸ Slater-Condon evals ≈ **~340s of the 636s** (dps/float-independent — unchanged by any
   ERI work).

**So the transformative minutes→seconds needs the two real walls, not more assembly flooring:**
- **FCI (tractable, biggest single lever):** numba/Cython the Slater-Condon `_matel` loop, or
  a determinant-frugal solver (Davidson/selected-CI) — kills the ~340s Python wall.
- **Seeds (harder, research-y):** a float64-STABLE Neumann seed recurrence that avoids the
  guard-digit need (the Q_l instability is why mpf is used) — or global `_B_table` caching by
  (m,s,l,p,c) across the tensor (a modest win where exponent-pairs share a rate).

## FCI LEVER DONE — sparse solver, bit-exact, and it UNBLOCKS M>16

`debug/fci_fast.py` `fci_energy_fast` — the dense `fci_energy` builds all nd² determinant
pairs (measured: 3.3s/72s/280s at M=10/14/16, nd²-scaling), but the CI Hamiltonian is SPARSE
(dets connect only if they differ by ≤2 spin-orbitals). The fast solver: bitmask each
determinant → find connected pairs by a VECTORIZED SWAR popcount of the XOR → compute
Slater-Condon (reusing `_matel`) only for connected pairs → ground state via `eigsh` (Lanczos),
no dense nd×nd.
- **GATE (vs dense, random h1/eri):** bit-exact (dE 1e-13), speedup grows with size:
  2×(M=8) → 3× → 4× → **7×(M=14)**; ~14× at M=16 (280s→~10s).
- **Strategic win beyond speed:** the FCI scaling goes dense nd² → sparse (nd × connections),
  so **M>16 is now reachable** — the regime holding the last radial/valence mHa toward ~−8.04
  that the dense wall blocked.

## COMBINED (float64 ERI + sparse FCI) — energy-exact, ~3× at M=16
`debug/prolate_float_eri.py validate` monkeypatches BOTH into the engine.
`debug/data/float_fci_validate.log`:

| config | E (float+sparse) | mpf | dE | time |
|:--|:--:|:--:|:--:|:--:|
| M=6 (1 core) | −7.99468 | −7.99468 | +0.004 mHa | 61 vs 147s |
| M=10 (1 core) | −8.00329 | −8.00329 | +0.004 mHa | 99 vs 237s |
| M=16 (3 cores) | **−8.02905** | −8.02905 | +0.005 mHa | **366 vs 1079s (3×)** |

Both optimizations preserve the banked energy to ~5 µHa. FCI 280s→~10s; the remaining 366s at
M=16 is now almost entirely the ERI.

## REMAINING WALL + next lever
The ERI (~350s at M=16) is now the sole bottleneck. It splits:
1. **mpf SEED tables** (intrinsically mpf — the guard digits tame the Q_l recurrence).
2. **mpf ETA moments in the per-integral loop** — `_sum_Y_m` (via `_Ymom_m`) is still mpf and
   recomputed per unique integral (the eta_poly differs per orbital pair). These are polynomial
   integrals with NO bad cancellation, so they are SAFE to float64 (cache-then-downcast, or a
   float eta pipeline) — the next tractable ERI lever (~1.5–2× more).

**Verdict:** #2 delivered a validated, correct productionization — **~3× at M=16, energy-exact
to 5 µHa, the FCI wall removed and M>16 unblocked.** Not the literal minutes→seconds (the mpf
seeds cap the ERI), but the load-bearing outcome — a fast sparse FCI that lets core enrichment
push past M=16 — is in hand. Remaining levers, scoped: float64 eta moments (tractable),
float64-stable seed recurrence (research). Deliverables: `debug/prolate_float_eri.py`,
`debug/fci_fast.py`, `debug/data/{dps_timing,float_validate,float_fci_validate}.log`.
