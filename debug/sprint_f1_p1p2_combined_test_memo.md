# Sprint F1 Phase 1 + 2 — Combined W1c × multi-zeta unified path

**Date:** 2026-05-23 (post-α arc, same day as the Track β synthesis).
**Sprint position:** F1 Phase 1 (architectural unification) + Phase 2 (combined-architecture FCI diagnostic). Tests the **mutual-exclusivity hypothesis** raised in the Track β synthesis: did the α-PES Layer 3 finding (FCI invisibility of multi-zeta at NaH max_n=2) emerge from a real structural fact, or was it an artifact of running multi-zeta WITHOUT W1c screening?
**Cross-references:** `debug/sprint_modular_alpha_arc_synthesis_memo.md` (Track β synthesis with three named follow-on directions; this sprint executes P1+P3 minus the basis-size leg, which becomes the next sprint), `debug/sprint_alpha_3_pes_test_memo.md` (α-PES three-step gate with bit-zero finding), `debug/sprint_alpha_2_multizeta_memo.md` (multi-zeta architecture), `geovac/balanced_coupled.py`, `geovac/cross_center_screened_vne.py`, `geovac/shibuya_wulfman.py`, CLAUDE.md §3 W1c-residual entries.

---

## §0. Executive summary + verdict

**Verdict line: PARTIAL-CLOSURE-AT-MAX_N=2.**

The PI's mutual-exclusivity hypothesis is **CONFIRMED** at the substantive structural level, but the deeper expectation (binding emergence at max_n=2 once the W1c × multi-zeta combined architecture runs) is **NOT met**. Three new structural findings emerge from the diagnostic.

1. **The original Layer 3 framing was specific to the bare cross-V_ne case.** With BARE cross-V_ne, Na 3s sits at h1 eigenvalue index i=5 (eigenvalue ~ −0.79 Ha), well above the lowest 5 H-localized orbitals, so the 2-electron FCI ground state has zero Na 3s amplitude — bit-zero finding from α-PES Step 2 reproduces verbatim here. With **W1c screening on** (either with or without multi-zeta), the H-side eigenvalue floor lifts from −3.99 Ha to ~−0.80 Ha (the un-screened Z=11 attraction on H 1s is dramatically reduced) and the lowest h1 eigenvalues become **bonding-like superpositions of Na 3s and H 1s**: at W1c without multi-zeta, the lowest two h1 eigenvectors at R=3.5 are dominantly (|H 1s|² = 0.964, eigvalue −0.802) and (|Na 3s|² = 0.986, eigvalue −0.790). The 2-electron FCI ground state has natural occupations [1.0000, 1.0000] (two doubly-occupied orbitals, one Na-localized and one H-localized) with **Na 3s diagonal occupation = 0.986** — Na 3s IS occupied in the W1c FCI. Layer 3's "FCI variational state H-side localization" finding is therefore **bare-cross-V_ne-specific**, not a general property of the framework's NaH max_n=2 basis.

2. **Multi-zeta is NOT bit-zero on the FCI when combined with W1c.** At every tested R, the W1c+multi-zeta FCI energy differs from W1c-alone by O(0.05–0.2) Ha, with the differential scaling cleanly from +0.208 Ha at R=2.0 down to +8.5×10⁻⁵ Ha at R=10.0 (decays to zero at dissociation as expected). The "bit-zero" α-PES Step 2 finding was an artifact of multi-zeta being applied to an FCI state that doesn't occupy Na 3s; with W1c the FCI occupies Na 3s and the multi-zeta substitution is fully load-bearing.

3. **The combined architecture does NOT close the binding gap.** Two-point D_e = E(R=10)−E(R=3.5) = +0.336 Ha for W1c+mz looks positive, but the extended PES scan (R ∈ {2.0, 2.5, ..., 10.0}) shows the PES is still monotonically descending: R_min = 2.0 bohr (smallest tested), no internal equilibrium minimum. **Combined architecture reduces the spurious-overattraction descent depth from W1c-alone 0.898 Ha → W1c+mz 0.690 Ha (23% reduction)** but doesn't bring the framework into a binding regime. The remaining ~0.69 Ha descent is the W1c-residual orthogonality wall, beyond what W1c-screening + multi-zeta-basis can address together.

### Predicted D_e^combined: **+0.336 Ha (literal two-point) / monotonically-descending PES (true)**

The literal `D_e = E(R=10) - E(R=3.5)` is positive, but this is the same spurious-binding signature the W1c-only baseline produced (Track 3: descent depth 0.357 Ha; here: 0.690 Ha across a wider grid). The PES does NOT have an internal minimum; the framework still over-attracts at small R.

### Structural reading of the three findings together

Layer 3 (FCI variational state H-localization at NaH max_n=2 under bare cross-V_ne) is **real but bare-specific**. The mutual-exclusivity reframe is what the PI's question asked: with W1c × multi-zeta unified, the FCI engages Na 3s at ~0.98 occupation. The remaining residual lives in the **cross-V_ne integration kernel SHAPE** (Track 3's named follow-up P2) and/or in BASIS SIZE (P1's max_n=3 test still needed to falsify/confirm whether the descent depth shrinks further at larger basis).

The W1c residual is not a single wall; it is a hierarchy of effects, and the unified architecture (this sprint) now lets all of them be measured cleanly.

---

## §1. Phase 1 — NotImplementedError diagnosis

### 1.1 The original error location

The `NotImplementedError` is in `geovac/balanced_coupled.py` line 748, in the `build_balanced_hamiltonian` cross-V_ne dispatch loop:

```python
if use_screened:
    if sb_multi_zeta:
        raise NotImplementedError(
            "multi_zeta_basis=True with screened_cross_center=True "
            "is not yet supported.  Disable screened_cross_center "
            "or skip multi_zeta_basis."
        )
    vne_matrix = compute_screened_cross_center_vne(...)
```

This was **defensive code** added in the α-Multi-zeta sprint (2026-05-23) to handle the case where multi-zeta is wired on a sub-block whose off-center nucleus has a frozen core, since the screened-V_ne kernel did not yet accept a multi-zeta basis dict.

### 1.2 Diagnosis: architectural gap, not fundamental incompatibility

The two production kernels handle different physics:

- **Bare path** (`compute_cross_center_vne` in `shibuya_wulfman.py`): uses `_radial_split_integral` (hydrogenic orbital, default) OR `_radial_split_integral_multizeta` (multi-zeta orbital, when supplied) to evaluate the multipole radial integral against `-Z_nuc/rho`.
- **Screened path** (`compute_screened_cross_center_vne` in `cross_center_screened_vne.py`): splits the potential into bare Coulomb `-Z_nuc/rho` (handled by `_radial_split_integral`, hydrogenic-only as written) plus screening correction `f_screen(rho)` (handled by grid quadrature using `_radial_wf_grid`, hydrogenic-only as written).

The bare path's multi-zeta extension (α-PES Sprint, 2026-05-23) was implemented by replacing the hydrogenic radial wavefunction in `_radial_split_integral` with a sum over multi-zeta primitives in `_radial_split_integral_multizeta`. The screened path's `_screened_radial_integral` calls BOTH `_radial_split_integral` AND `_radial_wf_grid` — neither was multi-zeta-aware as written, hence the screened path could not accept a multi-zeta dispatch.

**The two paths compose cleanly at the integrand level.** The bare Coulomb part of the screened integral is structurally identical to what the bare path computes (just multiplied by `-Z_nuc` at the end), and the screening correction is a smooth `f_screen(rho)` against the same orbital radial wavefunctions. Multi-zeta orbital evaluation is a drop-in replacement for `_radial_wf_grid` (via `MultiZetaOrbital.evaluate(r_grid)`), and `_radial_split_integral_multizeta` is already a drop-in replacement for `_radial_split_integral`.

**Conclusion**: the `NotImplementedError` is a clean architectural gap, mechanical to remove. No physics conflict. No re-derivation needed.

### 1.3 Pre-emptive verification (the substantive new content of Phase 1)

Smoke test of the combined-flag run at NaH max_n=2 R=3.5 (with the original NotImplementedError still in place) revealed:

```
No error raised. Multi-zeta diagnostics:
[{'sub_block': 'NaH_bond_center', 'Z_nuc_center': 11.0, ...}]
```

**The NotImplementedError did NOT fire for NaH.** Reason: at NaH max_n=2, the only multi-zeta-flagged sub-block (Na center) sees an off-center nucleus (H) that has NO frozen core, so it takes the bare path (with multi-zeta dispatch). The screened path is only invoked on the H sub-block looking at Na (which has [Ne] core), and the H sub-block has no multi-zeta flag (H Z=1 < 11). So the iteration where BOTH `use_screened` AND `sb_multi_zeta` would be true never happens for NaH.

This means **for NaH max_n=2, the combined architecture had been functional all along** — the NotImplementedError was preventing only the hypothetical case of a frozen-core-on-frozen-core encounter (e.g., NaCl: Na center sees Cl which has [Ne] core; both sub-blocks would need multi-zeta basis simultaneously). The α-PES Step 2 BARE-multi-zeta test was therefore not the "without W1c" complement of a hypothetical "with W1c" run; rather, it was a fully isolated bare-cross-V_ne probe, and the bit-zero finding was a real structural fact of that isolated setup.

This is the genuine Phase 1 substantive content: **the α-PES Step 2 bit-zero finding was correct but specific to the isolated bare-V_ne run; the W1c+mz combined architecture was structurally always reachable for NaH, just not previously tested.**

---

## §2. Phase 1 — Unified path implementation

### 2.1 Code changes

Three files modified:

**`geovac/cross_center_screened_vne.py`** (~25 lines added, no existing code modified):

- New import of `_radial_split_integral_multizeta` from `geovac.shibuya_wulfman`.
- `_screened_radial_integral` extended with optional `orbital_bra` and `orbital_ket` kwargs (default `None`). When BOTH supplied as `MultiZetaOrbital` instances:
  - The bare Coulomb sub-integral routes through `_radial_split_integral_multizeta` (sum of K_bra × K_ket pair-integrals, each analytical via `_split_integral_analytical`).
  - The screening-correction grid quadrature uses `orbital.evaluate(r_grid)` instead of `_radial_wf_grid(Z_orb, n, l, r_grid)`.
  - Grid `r_max` is widened to `40/zeta_min` (multi-zeta orbitals can have outer primitives with very small `zeta` for the diffuse tail).
- `compute_screened_cross_center_vne` gains a `multi_zeta_basis: Optional[Dict] = None` kwarg; dispatches to the multi-zeta-aware `_screened_radial_integral` when both `(n1, l1)` and `(n2, l2)` are in the dict (same conservative no-mixed-dispatch rule as the bare path).
- `compute_screened_cross_center_vne_element` gains `orbital_bra` / `orbital_ket` kwargs for symmetry with the new internal kernel.

**`geovac/balanced_coupled.py`** (NotImplementedError replaced):

```python
if use_screened:
    # Unified path (Sprint F1 Phase 1, 2026-05-23): the screened
    # kernel now accepts a multi-zeta basis dict and applies it to
    # both the bare-Coulomb analytical sub-integral and the
    # screening-correction grid quadrature. When sb_multi_zeta is
    # None or missing the relevant (n, l) pair, the screened path
    # falls back bit-exactly to the hydrogenic Z_orb baseline.
    vne_matrix = compute_screened_cross_center_vne(
        Z_orb, states, Z_nuc, R_AB,
        L_max=L_max, n_grid=n_grid_vne,
        direction=direction,
        multi_zeta_basis=sb_multi_zeta,
    )
```

The `sb_multi_zeta` dict (built earlier in the function, identical to the bare path's wiring) now flows through to BOTH the bare and the screened cross-V_ne kernels uniformly.

### 2.2 Test verification

5 new tests added to `tests/test_balanced_coupled_multizeta.py` under `TestUnifiedW1cMultiZetaPath`:

1. `test_combined_path_runs_no_error_nah_max_n_2`: the combined flag runs (was previously running for NaH due to sub-block disjointness; now confirmed across the unified codepath).
2. `test_backward_compat_screened_alone_unchanged`: with `multi_zeta_basis=False`, the screened path is bit-exact to the pre-F1 behavior.
3. `test_screened_kernel_backward_compat_no_multi_zeta`: `multi_zeta_basis=None` to `compute_screened_cross_center_vne` gives bit-exact pre-F1 behavior at the kernel level.
4. `test_combined_path_fci_finite_and_reasonable`: FCI in all four (bare / bare+mz / W1c / combined) configurations is finite, with W1c much higher than bare (reproducing Track 3 baseline pattern), and combined within 1 Ha of W1c-only.
5. `test_combined_path_shifts_h1_cross_vne_h_side`: combined gives `h1_cross_vne` shift of >0.01 Ha relative to W1c-alone (multi-zeta is NOT bit-zero on the matrix elements).

**All 5 new tests pass.** Full regression:

```
$ pytest tests/test_balanced_coupled_screened_valence.py \
         tests/test_balanced_coupled_multizeta.py \
         tests/test_phillips_kleinman_cross_center.py \
         tests/test_multi_zeta_orbitals.py \
         tests/test_shibuya_wulfman.py \
         tests/test_cross_center_screened_vne.py \
         -q --no-header
150 passed, 1 skipped in 32.30s
```

Zero regression across 150 baseline tests (128 previously + 5 new + 17 cross_center_screened_vne baseline that wasn't in the original mandate's regression set).

Pre-existing failures in `tests/test_balanced_row2.py` (29 failures) are **unrelated to this sprint** (verified via `git stash` test): they come from the v2.3.0 multi-center reorganization that moved `hcl_spec` etc. from `composed_qubit` to `molecular_spec`.

---

## §3. Phase 2 — Combined-architecture FCI results

### 3.1 Two-point FCI at R = 3.5 and R = 10.0 bohr

System: NaH max_n=2 (Q=20, M=10 spatial orbitals, n_electrons=2, FCI dim = C(10,1)² = 100).

| Architecture | E(R=3.5) [Ha] | E(R=10.0) [Ha] | D_e [Ha] | Na 3s occ (at R=3.5) | H 1s occ (at R=3.5) |
|:---|---:|---:|---:|---:|---:|
| bare (PK+SV-class baseline) | −169.203970 | −164.672120 | +4.5319 | **0.000** | 1.086 |
| bare + multi-zeta | −169.203970 | −164.672120 | +4.5319 | **0.000** | 1.086 |
| W1c (screened, no mz) | −163.170854 | −162.778881 | +0.3920 | **0.986** | 0.964 |
| **COMBINED (W1c + multi-zeta)** | **−163.115069** | **−162.778796** | **+0.3363** | **0.981** | 0.964 |

Three structural observations from the two-point table alone:

(a) **Multi-zeta is bit-zero on bare-baseline** (E_bare = E_bare+mz to floating-point precision at both R values). This reproduces the α-PES Step 2 finding and confirms the structural reason: bare-baseline FCI has zero Na 3s amplitude.

(b) **W1c screening lifts the bare-baseline E by ~6 Ha** (from −169.2 → −163.2 at R=3.5). The un-screened Z=11 nuclear attraction on H 1s was the dominant spurious-binding contributor (Track 3 finding); W1c screening removes most of it.

(c) **Multi-zeta is fully visible when combined with W1c** (E_W1c=−163.171 vs E_W1c+mz=−163.115, differential **+0.056 Ha** at R=3.5; decays to +8.5×10⁻⁵ at R=10). The mutual-exclusivity reframe is confirmed at the substantive structural level: multi-zeta is load-bearing when the FCI engages Na 3s, and the FCI engages Na 3s when W1c removes the spurious-overattraction wall.

### 3.2 Extended PES scan — is there an internal minimum?

The two-point D_e looks positive in all four cases, but D_e from two points is misleading when the PES is monotonically descending and just sampled at two non-equilibrium points. Extended scan at R ∈ {2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0}:

| R [bohr] | bare | bare+mz | W1c | W1c+mz | mz-vs-W1c diff |
|:---:|---:|---:|---:|---:|---:|
| 2.00 | −172.941 | −172.941 | −163.677 | −163.468 | **+0.2084** |
| 2.50 | −171.216 | −171.216 | −163.436 | −163.299 | +0.1364 |
| 3.00 | −170.026 | −170.026 | −163.279 | −163.192 | +0.0878 |
| 3.50 | −169.204 | −169.204 | −163.171 | −163.115 | +0.0558 |
| 4.00 | −168.601 | −168.601 | −163.092 | −163.056 | +0.0351 |
| 5.00 | −167.671 | −167.671 | −162.984 | −162.970 | +0.0135 |
| 7.00 | −166.161 | −166.161 | −162.866 | −162.864 | +0.0019 |
| 10.00 | −164.672 | −164.672 | −162.779 | −162.779 | +8.5e-05 |

| Architecture | R_min | E_min | Descent depth (E(R=10) − E_min) | Internal min? |
|:---|---:|---:|---:|:---:|
| bare | 2.00 | −172.941 | 8.269 Ha | **NO** |
| bare+mz | 2.00 | −172.941 | 8.269 Ha | **NO** |
| W1c | 2.00 | −163.677 | 0.898 Ha | **NO** |
| **W1c+mz** | **2.00** | **−163.468** | **0.690 Ha** | **NO** |

All four architectures show monotonically descending PES (R_min at the smallest tested R = 2.0 bohr). The combined W1c+mz architecture has the SHALLOWEST descent depth (0.690 Ha = 23% reduction from W1c-alone's 0.898 Ha), confirming multi-zeta has structural content beyond W1c-alone — but the binding gap is **not closed**. Experimental NaH D_e ≈ 0.075 Ha; the combined architecture over-attracts by ~10× at R=2.0 vs the binding scale.

The **literal two-point D_e** (+0.336 Ha for W1c+mz) gave the misleading appearance of binding because R=3.5 and R=10.0 happen to be on opposite sides of the steep over-attraction region; the actual PES has no equilibrium.

### 3.3 Multi-zeta differential behavior

The `mz-vs-W1c` column in §3.2 is informative: the multi-zeta substitution gives a smooth R-dependent shift that's positive (multi-zeta is LESS over-attractive than hydrogenic Na 3s) and decays cleanly to zero at dissociation. This is consistent with the α-PES Step 1 finding (cross-V_ne kernel differential = −0.135 Ha at R_eq with the differential vanishing at R=10), but with a **sign flip** in interpretation: under bare baseline the differential was on the unoccupied Na 3s diagonal and bit-zero on the FCI; under W1c the differential is on the OCCUPIED Na 3s and propagates to a finite FCI shift.

The Step 1 finding said multi-zeta makes the H nucleus pull Na 3s MORE strongly than hydrogenic; under W1c (which screens the Na nucleus's pull on H), the Na 3s orbital is now the dominant Na-side bonding orbital, so multi-zeta IS load-bearing.

---

## §4. Na 3s occupation diagnostic — the headline structural finding

### 4.1 FCI 1-RDM analysis

Computed the spatial 1-RDM from the FCI ground-state eigenvector (driver: `debug/sprint_f1_p1p2_combined_test.py`'s `compute_1rdm` function; standard a†_p a_q construction summed over alpha and beta spin sectors). Diagonal element `rdm[0, 0]` corresponds to the Na 3s orbital (index 0 in the NaH spec at max_n=2).

| Architecture | Na 3s diag occ | H 1s diag occ | tr(RDM) = N | Natural occs (top 6) |
|:---|---:|---:|---:|:---|
| bare baseline | **0.000** | 1.086 | 2.000 | [1.993, 0.007, 0, 0, 0, 0] |
| bare + multi-zeta | 0.000 | 1.086 | 2.000 | [1.993, 0.007, 0, 0, 0, 0] |
| W1c (screened) | **0.986** | 0.964 | 2.000 | [1.000, 1.000, 0, 0, 0, 0] |
| COMBINED W1c+mz | **0.981** | 0.964 | 2.000 | [1.000, 1.000, 0, 0, 0, 0] |

The natural-occupation patterns reveal the structural picture clearly:

- **Bare**: dominant single-determinant character with one occupied orbital (occupation 1.993) — a single H-localized orbital is the doubly-occupied bonding pair. Na 3s is structurally absent from the wavefunction.
- **W1c (alone or with mz)**: two singly-occupied natural orbitals with occupations [1.000, 1.000] — a CI-like character with two near-degenerate orbitals, one Na-localized and one H-localized, each with ~1.0 occupation. This is exactly the open-shell-singlet-like situation expected for a separated-pair bonding regime.

### 4.2 h1 eigenspectrum at R=3.5 — where does Na 3s sit?

| Architecture | Lowest 5 h1 eigvals | Na 3s amplitude² in lowest 5 | H 1s amplitude² in lowest 5 |
|:---|:---|:---|:---|
| bare | [−3.986, −3.102, −2.211, −2.211, −1.636] | [0, 0, 0, 0, 0] | [0.647, 0.323, 0, 0, 0.030] |
| W1c | [−0.802, −0.790, −0.479, −0.419, −0.315] | [0, **0.986**, 0, 0.011, 0] | [**0.964**, 0, 0.032, 0, 0] |
| W1c+mz | [−0.802, −0.735, −0.479, −0.418, −0.315] | [0, **0.981**, 0, 0.015, 0] | [**0.964**, 0, 0.032, 0, 0] |

This is the genuinely new content: **Na 3s sits at h1 eigenvalue index i=1 (the SECOND-LOWEST eigenvector) with amplitude² ≈ 0.98** when W1c is on. Under bare cross-V_ne the lowest 5 eigenvalues are at −4.0 to −1.6 Ha (dominantly H-localized, dragged deep by un-screened Z=11), and Na 3s is at i=5 with eigenvalue ~−0.79 Ha (above the bonding manifold).

Under W1c screening, the H-side dragging is removed; the lowest h1 eigenvalues are at ~−0.80 Ha (dominantly H 1s) and ~−0.79 Ha (dominantly Na 3s). These two near-degenerate orbitals get singly occupied in the 2-electron singlet ground state.

The **multi-zeta substitution lifts the Na 3s eigenvalue from −0.790 to −0.735 Ha** (Na 3s becomes less bound when its compact physical-fit shape is substituted for the diffuse hydrogenic Z=1 placeholder). This is the algebraic kernel finding from α-PES Step 1 (the cross-V_ne differential, +0.056 Ha at R=3.5 on the Na 3s diagonal) propagating to the h1 eigenspectrum. Multi-zeta is structurally visible in this regime.

### 4.3 Reconciliation with α-PES Layer 3

The Track β synthesis memo §4 stated Layer 3 as a NaH max_n=2 basis-size finding: "the framework's NaH max_n=2 basis cannot represent the physical bonding configuration where electrons share between Na 3s and H 1s in a quantum superposition. With un-screened cross-V_ne the H side dominates; with W1c screening the H side is no longer artificially deep but the Na 3s diagonal at −0.78 Ha is still well below the H 1s + screened-cross-V_ne energy (~−0.5 Ha), so the FCI is still H-dominant."

**This is partially wrong.** The W1c FCI here shows natural occupations [1.000, 1.000] with Na 3s occupation 0.986, NOT H-dominant. Layer 3's "FCI is still H-dominant under W1c" claim from the synthesis memo was an inference, not a measurement; the actual FCI under W1c is an open-shell-singlet-like CI state with substantial Na 3s amplitude.

What IS true (and which the new diagnostic confirms): at max_n=2 the FCI lacks the structural flexibility to produce an internal-minimum bonding configuration. The natural-occupation pattern [1.0, 1.0] is a CI state, but it's not a bond — it's two near-degenerate orbitals occupied separately (one Na, one H). For a true Na-H bond, the framework would need to mix Na 3s + H 1s into a bonding orbital + antibonding orbital pair with one of each occupied, and the bonding orbital would need lower energy at R_eq than at R=∞. Currently the framework's bonding combinations have higher energy than the H-Na separated configuration, so the FCI prefers the separated-pair state at every R; the resulting PES is just the screened atomic-attraction over-attraction.

**Revised Layer 3 statement: at NaH max_n=2 under W1c + multi-zeta, the FCI engages Na 3s at full occupation (~0.98), but the structural-flexibility-limited basis cannot construct a bonding-vs-antibonding partition that places the bonding combination at a lower energy than the separated configuration. The wall is now at the level of orbital-orbital mixing, not at the level of FCI occupation.**

---

## §5. Verdict + corroborating evidence

### 5.1 Verdict line: PARTIAL-CLOSURE-AT-MAX_N=2

The verdict per the original gate rules:

> - If FCI Na 3s occupation > 0.1 AND $D_e^\text{combined} > 0$ (binding emerges): "W1c Layer 3" was a mutual-exclusivity artifact; the combined architecture closes (or substantially reduces) the wall at max_n=2.
> - If FCI Na 3s occupation > 0.1 BUT $D_e^\text{combined} \le 0$ (engagement without binding): partial closure of Layer 3 (FCI engages but cross-V_ne kernel still wrong); pivot to F2 (cross-V_ne kernel-shape substitution).
> - If FCI Na 3s occupation < 0.01 (still no engagement): "W1c Layer 3" corroborated at max_n=2; basis-size test (max_n=3) is the next experiment.

**Na 3s occupation under combined = 0.981 (well above 0.1) and literal D_e = +0.336 Ha (numerically > 0). The naive gate says W1C-LAYER3-EVAPORATES.**

However, the extended PES scan (§3.2) shows the W1c+mz PES is monotonically descending with no internal minimum; the "+0.336 Ha" looks like binding only because R=3.5 sits above the over-attraction floor at R=2.0. **The correct read is PARTIAL-CLOSURE-AT-MAX_N=2**: FCI engagement is achieved (the most substantial improvement multi-zeta delivers), but binding-with-equilibrium is not. The descent-depth reduction (0.898 → 0.690 Ha, 23%) is the measurable improvement; the residual ~0.69 Ha overattraction is the remaining wall.

### 5.2 What the mutual-exclusivity reframe DID close

(a) **The "bit-zero on FCI" mystery is solved.** Multi-zeta is bit-zero only when applied to an FCI state that has no Na 3s amplitude. With W1c the FCI has 0.98 Na 3s occupation, and multi-zeta is fully load-bearing (R-dependent shift +0.21 Ha at R=2.0 down to bit-zero at R=10).

(b) **The bimodule diagnostic (M-Y / α-1) is vindicated at the FCI level.** The d_R/d_L = 3.23 prediction from α-1 was that the right-action axis (Na-side wavefunction shape) is the dominant residual. Once the basis admits Na-3s engagement (via W1c), the bimodule's predicted axis shows up at the FCI eigenvalue level: multi-zeta on the Na side gives an O(0.05–0.2) Ha shift in the FCI ground state. The α-PES Step 2 bit-zero finding was real (under bare-V_ne the bimodule diagnostic predicted the right wavefunction-shape distance but the FCI didn't occupy that orbital); the unified architecture confirms the prediction.

(c) **The W1c residual is a hierarchy, not a single wall.** The α arc's three-layer taxonomy refines:

  - Layer 1 (H–Na orthogonality, PK barrier): produces a 14.6% reduction in W1c descent depth.
  - Layer 2 (Na-side wavefunction shape, multi-zeta basis): now CLEARLY load-bearing under W1c, producing an additional 23% descent-depth reduction.
  - Layer 3 (FCI variational state flexibility): the genuine remaining wall — the natural occupations [1, 1] indicate a CI of two separately-occupied orbitals, NOT a bonding-antibonding partition with one bonding pair. To get binding-with-equilibrium at this basis, the framework would need to produce orbital combinations where one mixed (Na 3s + H 1s) combination is at lower energy than each separate orbital; max_n=2 may not have enough basis flexibility.

### 5.3 What's NOT closed

The PES still over-attracts by ~10× at R=2.0 vs experimental D_e. The bonding-vs-antibonding orbital mixing is not achieved at max_n=2 with W1c + multi-zeta. To close: either (i) larger basis size at max_n=3 to let the FCI mix more orbital combinations (and possibly find a true bonding-vs-antibonding partition), (ii) cross-V_ne integration kernel SHAPE substitution (Track 3's named target P2, which addresses the bare cross-V_ne shape against multi-zeta orbital rather than just the on-site diagonal), or (iii) a structurally different solver (post-HF with explicit bonding ansatz) that doesn't rely on variational FCI to find a bond.

---

## §6. Recommended next sprint

Three named candidates ranked by expected new content / cost:

### Priority 1 — Increase basis to NaH max_n=3, with full W1c × multi-zeta unified architecture (NEW: now mechanical, 3-5 days)

**Why first:** with Phase 1 unification done, the W1c × multi-zeta combined run is mechanical at max_n=3. The structural question is whether the larger basis (Q=56, 5,349 balanced Pauli terms per CLAUDE.md §2 NaH n_max=3 v2.1.1) lets the FCI construct genuine bonding orbital combinations. Expected outcomes:

- If the W1c + multi-zeta combined PES at max_n=3 shows an **internal minimum** (R_eq near experimental 3.566 bohr, D_e closer to experimental 0.075 Ha): W1c-residual wall is closed in production; multi-zeta architecture is the right closure; the framework supports second-row hydride PES at appropriate basis size.
- If the PES is still monotonically descending at max_n=3: the wall is not basis-size limited at this scale; Priority 2 (kernel-shape substitution) is the next correct direction.
- If the PES is descending but with a substantially shallower descent depth (e.g., 0.2 Ha instead of 0.69): the wall is partially basis-size limited, larger basis (max_n=4) might close further but at exponentially-growing FCI cost; pragmatic answer is to declare W1c partially closable at modest basis sizes and document the n_max scaling.

This sprint tests the genuine Phase 2 hypothesis: does max_n=3 + combined architecture produce binding? Sub-week effort; one driver script + four PES scan runs at different R.

### Priority 2 — Cross-V_ne kernel-shape substitution (Track 3 diagnostic's named target, 1-2 weeks)

**Why second:** Track 3's diagnostic (2026-05-09) identified that the genuine load-bearing target is cross-V_ne integration shape, not diagonal h1 substitution. The α-PES Step 1 finding (-0.135 Ha differential for multi-zeta on Na 3s on-site) is a smaller-effect cousin of the Phase D D1 finding (-0.674 Ha when substituting physical-n hydrogenic shape in cross-V_ne). Mechanically this is a refactor of `shibuya_wulfman._radial_split_integral` to use physical-n orbital decomposition on the heavy-atom-side **partner** sub-block in cross-V_ne — distinct from the multi-zeta basis on the heavy-atom center sub-block.

**The wiring is parallel to the F1 Phase 1 work but on a different sub-block dimension.** The "physical-n hydrogenic" form is qualitatively different from the multi-zeta basis (it's a single-zeta with the correct radial node count, not a full multi-zeta expansion).

### Priority 3 — Pivot to a different system class (LiH+, BeH⁺, etc.) where the framework's basis structure is known to support binding

**Why third (or in parallel):** the W1c-residual wall is empirically Z-decreasing (NaH 5.4-6.0×, MgH₂ 2.99×, HCl 1.79× per CLAUDE.md §2). At MgH₂ and HCl the wall may already be smaller than the bond energy scale, so the framework may have effective binding even without a complete W1c closure. Test the combined W1c + multi-zeta architecture on HCl at max_n=2 to see if the residual is below the bond energy scale.

### Bundled recommendation

**Default**: open Priority 1 (max_n=3 NaH FCI with combined architecture) as the lead next sprint. Cost ~3-5 days. If POSITIVE → W1c closure achieved in production; if NEGATIVE → confirms the wall is not basis-size-limited at this scale and Priority 2 (kernel-shape) is the next correct direction.

**Parallel**: Priority 3 HCl/MgH₂ run with combined architecture (2-day scoping probe) to test whether the framework already supports second-row hydride binding for less-residual systems even without the W1c wall being fully closed.

---

## §7. Files

### Created (debug)

- `debug/sprint_f1_p1p2_combined_test.py` — Phase 2 main driver: two-point FCI + 1-RDM + h1 eigenspectrum + mini-PES + verdict logic. ~450 lines.
- `debug/sprint_f1_p1p2_extended_pes.py` — Phase 2 follow-on: 8-point PES scan across all four architectures to verify internal-minimum status. ~110 lines.
- `debug/data/sprint_f1_p1p2_results.json` — Main numerical results from the F1 P2 driver.
- `debug/data/sprint_f1_p1p2_extended_pes.json` — Extended PES scan data.
- `debug/sprint_f1_p1p2_combined_test_memo.md` — This memo (~3500 words).

### Modified (production)

- `geovac/cross_center_screened_vne.py` — Added `_radial_split_integral_multizeta` import; extended `_screened_radial_integral` with `orbital_bra`/`orbital_ket` kwargs; extended `compute_screened_cross_center_vne` and `compute_screened_cross_center_vne_element` with `multi_zeta_basis` kwarg. ~25 lines added; zero existing-function bodies modified.
- `geovac/balanced_coupled.py` — Replaced NotImplementedError-raising block in `build_balanced_hamiltonian` with unified dispatch to `compute_screened_cross_center_vne(multi_zeta_basis=sb_multi_zeta, ...)`. ~10 lines modified.

### Created (tests)

- `tests/test_balanced_coupled_multizeta.py` — Added `TestUnifiedW1cMultiZetaPath` class with 5 new tests covering: combined-path no-error run, backward-compat at balanced and kernel level, FCI sanity, h1_cross_vne shift verification. All 5 pass. 150/151 regression (1 pre-existing skip).

### Not modified

- Any paper `.tex` files (per sprint mandate: "this is a diagnostic sprint, paper edits come later").
- `geovac/multi_zeta_orbitals.py` (the Z=11 fits from α-2 are used as-is).
- `geovac/screened_valence_basis.py` (no changes needed; compose-with-screened-valence-basis worked).
- `geovac/phillips_kleinman_cross_center.py`, `geovac/cross_register_vne.py`, `geovac/neon_core.py`, `geovac/shibuya_wulfman.py` (none used new features).
- CLAUDE.md (per the deliverable rules: this is a diagnostic-arc memo, CLAUDE.md updates can be batched in a Track-β-style synthesis after the PI reviews).

---

**End of Sprint F1 Phase 1+2 memo.**
