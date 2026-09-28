# g_Vee defect audit: the general-m triangle weight (2026-09-26)

Follow-up to `debug/sprint_lih_exp_geminal_memo.md` Job 2 (the linexp c-offset). **Verdict: FOUND and
FIXED.** One weight in the general-m triangle of ⟨F²V_ee⟩ was wrong. After the fix the analytic g_Vee and c
agree with component VMC within 1σ at linexp γ=0.25 and 0.5. geovac/ is not affected and was not edited.

## Root cause

`debug/lih_pi_m16.py` `get_tri_matrix_m` (dense) and `debug/lih_gvee_sparse.py` `tri_entries` (sparse,
production) each carried the same weight for the middle-vertex mode combo (A, sa, B, sb):

```python
def wcomb(A, sa, B, sb):
    if A == 0 and B == 0:
        return 2 * np.pi
    return 0.5 * (2 * np.pi if A == 0 else np.pi)      # tests A only
```

The correct weight comes from f(|r1−r3|) = Σ_a F̂_a e^{ia(φ1−φ3)}, with F̂_0 = F_0 and F̂_{±A} = F_A/2, the same
expansion for f(|r2−r3|), and ∫dφ3 = 2π δ(a+b = μ_mid). That gives

w = 2π · (1 if A = 0 else ½) · (1 if B = 0 else ½).

For (A ≥ 1, B = 0) the code gave π/2; the correct value is π. Those combos exist only when the middle
(shared-vertex) pair density carries circulation μ_mid ≠ 0, i.e. only on a π reference. The σ marriage,
the 2-MO engine and every μ_mid = 0 triangle are therefore untouched. The error is L/R-asymmetric: A
belongs to the left f-edge. That asymmetry is why the dense==sparse guard stayed green (both used the same
weight) while |Ψ|²'s relabeling symmetry was broken. The (A≥1, B=0) contribution is signed and its sign
depends on γ, which explains the observed sign flip with γ.

**Fix:** one shared function, `lih_pi_m16.tri_wcomb` (new weight). Both call sites now use it. A
`_WCOMB_LEGACY` hook, used only by the fire test, restores the old value.

## How it was localized

**1. σ-only 2-MO reference, exact i.i.d. MC** (`gvee_audit_2mo_{analytic,mc,compare}.py`).
- Each spin block of |φ0²φ1²| is sampled exactly: marginal, then the exact conditional orbital,
  drawn by rejection. No grid, kernel or Neumann is involved.
- 3.6×10⁸ samples at linexp 0.25. The analytic g_Vee on the σ path is 0.877161; MC gives 0.877247 ± 0.000129
  (+0.67σ).
- Every one of the 9 topology groups is within 2.1σ, at the 1e-4 level.
- exp 0.5 gives 0.535896 vs 0.535889 ± 0.000075.

**Conclusion:** the m=0 machinery (PAIR / f²-coul / f-over-r / STAR / TREE4 / BRIDGE / TRI) is exact. The
defect must be m≠0.

**2. Symmetry-orbit test on the π trunc-1e-2 reference** (`gvee_audit_pi_triples.py`). All 216 triples were
evaluated per triple through the production sparse path. |Ψ|² is invariant under (12), (34) and (13)(24).
Of all 216 triples grouped into orbits, only one orbit was non-degenerate:

| orbit | members | old values | spread old → new |
|---|---|---|---|
| f₁₂f₁₃/r₂₃ (TRI, intra × inter f; μ_mid on the intra side) | 16 | α-middle 1.195730, β-middle 1.195322 | 4.1e-4 → 4.7e-9 |
| every other orbit | — | — | ≤ 1.4e-14 |

This pointed straight at the triangle's mode bookkeeping. The 4.7e-9 residual after the fix is
randomized-SVD noise.

**3. Independent π-reference Metropolis MC** (`gvee_audit_pi_{mc,compare}.py`).
- Lean sampler on |Ψ_CI|², with no guide and no geminal. 5 seeds × 4000 walkers × 38000 sweeps.
- Groups are compared in centered form Σ⟨(f_i−μ_i)(f_j−μ_j)v_k⟩.
- Block jackknife over 250 blocks; blocks are 760 sweeps, so the errors may be mildly optimistic.

linexp γ=0.25:

| group (centered) | MC | ± | analytic OLD | σ | analytic NEW | σ |
|---|---|---|---|---|---|---|
| BRIDGE | +0.035588 | 0.000072 | +0.035688 | −1.4 | +0.035688 | −1.4 |
| PAIR f²/r | +0.948069 | 0.000180 | +0.947707 | +2.0 | +0.947707 | +2.0 |
| STAR | −0.111573 | 0.000031 | −0.111549 | −0.8 | −0.111549 | −0.8 |
| TREE4 | −0.107471 | 0.000114 | −0.107580 | +1.0 | −0.107580 | +1.0 |
| **TRI** | **+0.059821** | **0.000046** | **+0.052693** | **+154** | **+0.059732** | **+1.9** |
| f/r × f disjoint | −0.221879 | 0.000134 | −0.221597 | −2.1 | same | −2.1 |
| f/r × f share | −0.914456 | 0.000236 | −0.914399 | −0.2 | same | −0.2 |
| f² × coul disjoint | +0.293385 | 0.000038 | +0.293383 | +0.1 | same | +0.1 |
| f² × coul share | +0.958549 | 0.000092 | +0.958503 | +0.5 | same | +0.5 |
| **all = g_Vee** | **+0.940033** | **0.000305** | **+0.932849** | **+23.5** | **+0.939888** | **+0.5** |

Only TRI disagreed, and the fix moves it by +0.00704, which is the whole offset. The residual +2σ on PAIR and
−2.1σ on "f/r × f disjoint" are opposite in sign and sit at the 2σ level with possibly optimistic errors. They
are below 4e-4 each, so no action was taken; flagged as open.

Traps checked and cleared:
- linexp f² vs exp kernels: the PAIR/Ef kernels are correct (d·e^{−2γd}), and PAIR agrees on both references.
- self-pair f²/r: handled as a pair term with the right kernel.
- F̄ cancellation: g is linear in the moments, and the centered comparison is independent of the F̄ ordering.
- Yukawa dressing: no bias on the exact path (exp 2-MO is within 0.1σ).
- the 2(l_bridge + l_leaf) leaf bound: BRIDGE agrees.
- 1s shortcut: none on the generalized path.

## Gates

**1. Analytic g_Vee vs component VMC (cvmc2, 12 seeds):**

| | VMC | old | new | new σ |
|---|---|---|---|---|
| linexp0.25 g_Vee | 0.93964 ± 0.00115 | 0.93286 (−5.9σ) | 0.93990 | +0.2 |
| linexp0.5 g_Vee | 0.64387 ± 0.0012 (v2) / 0.64443 (v1) | 0.64649 | 0.64398 | +0.1 / −0.4 |

**PASS.**

**2. c vs VMC:**

| | VMC c | old | new | new σ |
|---|---|---|---|---|
| linexp0.5 | 0.22287 ± 0.00037 | 0.221588 (−3.5σ) | 0.223191 | +0.9 |
| linexp0.25 | 0.15976 ± 0.00023 | 0.161465 (+7.4σ) | 0.159979 | +1.0 |

**PASS.**

**3. exp energies (trunc 1e-2, re-assembled with the delta tool):**

| γ | ΔF2Vee | E old | E new | ΔE (mHa) | VMC err (mHa) |
|---|---|---|---|---|---|
| 0.5 | −0.001302 | −8.039002 | −8.039066 | −0.064 | 0.39 |
| 1.0 | −0.000239 | −8.046711 | −8.046725 | −0.015 | 0.30 |
| 1.5 | −0.0000514 | −8.048138 | −8.048141 | −0.003 | — |
| 2.0 | −0.0000139 | −8.047573 | −8.047574 | −0.001 | 0.20–0.29 |

**PASS:** every shift is below the VMC error. The exp c values move by ≤ 1e-3 (γ=2: −0.228692 → −0.228699;
VMC −0.22834).

linexp energies:

| γ | E old | E new | ΔE (mHa) |
|---|---|---|---|
| 0.25 | −8.044831 | −8.044650 | +0.181 |
| 0.5 | −8.042181 | −8.042304 | −0.123 |
| 1.0 | −8.031278 | −8.031303 | −0.025 |
| 2.0 | −8.025205 | −8.025205 | 0.000 |

The trunc-1e-2 best is still exp γ=1.5, now −8.048141. The ±0.2 mHa systematic carried for the "open
g_Vee defect" can be dropped.

**4. σ reference −8.0372357:** unchanged by construction.
- A σ reference has mA ≡ 0, so μ_mid = 0 for every middle vertex, and the (A≥1, B=0) combos never occur.
- The delta tool skips μ_mid = 0 exactly, giving Δ = 0.
- Not re-run.

**5. Existing guards:**
- `tests/test_lih_gvee_sparse.py` (4, --slow): **PASS**. It compares sparse to dense, which share the weight,
  so it cannot see this class of bug.
- `tests/test_lih_r12ci*.py` (18, --slow) run on their own: **18 PASS**.
- **Pre-existing test-isolation defect (not caused by this fix):** run in the same session as
  `test_lih_gvee_sparse.py`, 4 r12ci tests FAIL. The `lih_gvee_sparse` fixture sets
  `KG.USE_EXACT_NEUMANN = True` module-globally and never restores it. The r12ci anchors then see the exact
  path, giving exp −7.943605 against the legacy −7.9420 anchor, and the default-switch test sees True.
  - Fix owed in that fixture: save and restore the flag, as the new guard does.
  - Owner call: I did not edit the existing test.

## geovac/lih_r12ci (the 2-MO engine)

It does not share the defect: `triangle.triangle_raw` is m-diagonal with a φ-independent middle density, so
there are no (A, B=0) combos. No patch was made. The banked 2-MO energies (linexp −7.917199 / exp −7.943605,
exact-Neumann path) are unchanged.

The old "+0.53834 vs VMC +0.53605" 2-MO discrepancy was the *legacy* Neumann path. On the exact path, exp
γ=0.5 gives g_Vee = 0.535896, and my exact MC gives 0.535889 ± 0.000075. Under the old VMC's noise that matches
0.53605.

## Regression guard (separate step, fire-tested)

`tests/test_lih_gvee_triangle_symmetry.py` (@slow, 4 tests, about 27 min, trunc-5e-2 π reference). **4 PASS.**
- `tri_wcomb(A,±,0,0) == tri_wcomb(0,0,A,±) == π`.
- The 8 α/β members of the f₁₂f₁₃/r₂₃ orbit are degenerate to 1e-6 relative.
- **Built-in fire test:** with `M._WCOMB_LEGACY = True` the orbit spread must exceed 1e-5 relative. It does,
  so the guard rejects the old weight.
- **Wrong answer rejected:** the (A≥1, B=0) half weight, which gave linexp0.25 g_Vee 0.93286 instead of
  0.93990.
- It restores `USE_EXACT_NEUMANN` on teardown.

## Re-assembly tool; the running trunc-1e-3 job

`debug/gvee_wcomb_delta.py --npz <run>.npz` replays only the 24 triangle triples. It evaluates only the
(A≥1, B=0) combos at the delta weight π/2, adds Δ to F2Vee/g_Vee and re-solves the 2×2. It writes
`<run>_wcombfix.npz` and leaves the original untouched.

Validated at linexp0.25: delta-tool g_Vee 0.9399013 vs full rerun 0.9399011 (2e-7, randomized-SVD noise). A run
takes 2–4 min at trunc 1e-2.

The exp γ=1.5 trunc-1e-3 job imported the old code, so its g_Vee carries the bug. When it finishes, run:

```
python debug/gvee_wcomb_delta.py --npz debug/data/lih_exp_geminal_exp_g1.5_72x44x28_p28-48_1e-03.npz
```

Expect |ΔE| of order 0.003 mHa, going by the trunc-1e-2 value. It was not run here: the job was still in
g_Vee/h_T, and the RAM budget ruled it out.

## Files

- **Fixed:**
  - `debug/lih_pi_m16.py`: new `tri_wcomb` + `_WCOMB_LEGACY`; `get_tri_matrix_m` uses it.
  - `debug/lih_gvee_sparse.py`: `tri_entries` uses `M.tri_wcomb`.
  - Both files are **untracked** in git; the legacy hook preserves the old behaviour.
- **New drivers:** `debug/gvee_audit_2mo_{analytic,mc,compare}.py`, `debug/gvee_audit_pi_{triples,mc,compare}.py`,
  `debug/gvee_wcomb_delta.py`.
- **Data:**
  - `debug/data/gvee_audit_*.npz`
  - `debug/data/lih_exp_geminal_*_1e-02*_wcombfix.npz`
  - logs in `debug/data/gvee_audit/`
- **Guard:** `tests/test_lih_gvee_triangle_symmetry.py`.
