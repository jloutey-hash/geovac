# Sprint MH Track D — Pachucki Dual Expansion (Verify-and-Document)

**Date:** 2026-05-08
**Status:** STRUCTURAL NEGATIVE on the directive's primary premise; documentation closure + diagnostic utility added.
**Verdict:** Regime dispatch in `geovac/cross_register_vne.py` is **structurally moot for the production path.** Track A's 24% SE gap is **not in this module**.
**Files:**
- `geovac/cross_register_vne.py` — added `roothaan_J0_taylor_expansion_dual` and `roothaan_recoil_shift_regime_aware` as parallel diagnostic utilities + module-level commentary block locating the SE gap correctly.
- `tests/test_cross_register_vne.py` — added `TestSprintMHTrackDDualExpansion` (7 tests, all pass).
- No paper edits (no observable improvement; the SE gap doesn't close from this work).

## Bottom line

**Did the dual expansion close Track A's 24% SE gap? No.** It cannot.

The 0.83 vs 0.668 meV gap (Track A's `self_energy_eides_lepton`, residual −0.16 meV vs Antognini 2013) is in the Eides Sec.3.2 leading-order SE bracket inside `debug/sprint_mh_track_a.py`. That formula does not call into `cross_register_vne.py` at all. Closing it requires either (a) extending the Eides bracket with the next-to-leading α(Zα)⁴(m_red/m_p) recoil-mixing terms (Bodwin–Yennie / Pachucki recoil-SE — same LS-8a wall as multi-loop QED, deferred per CLAUDE.md §2), or (b) implementing field-theoretic vertex renormalization at the spectral level (the LS-8a-renorm sprint, also deferred). Neither is a Roothaan-kernel adjustment.

**New muonic Lamb shift residual against CREMA: unchanged at −0.10%** (framework + literature) and −0.92% (framework-native). Track A's Antognini-1-ppm match on the Uehling kernel is preserved bit-identically. Track B's electronic recoil match (2.86% vs Bethe–Salpeter) and muonic recoil prediction are bit-identical.

## What the directive asked for vs. what's actually true

The directive was structured around the hypothesis: "Roothaan recoil kernel is regime-limited: λ_μ > λ_n breaks the large-nucleus expansion. Pachucki dual expansion exists but isn't specialized to muonic input." The directive itself anticipated the verify-and-document branch: *"If the kernel uses J_0 directly (not its expansion), the regime change is moot — verify this first and document the finding."*

The verification result is unambiguous, and the verify-and-document branch is the right answer:

### Finding 1: `_roothaan_J0` is symmetric and regime-agnostic

The closed form
```
J_0(λ_e, λ_n) = λ_e λ_n (λ_e² + 3 λ_e λ_n + λ_n²) / (λ_e + λ_n)³
```
is manifestly symmetric in (λ_e, λ_n). Numerical verification at the two regime endpoints:

| Inputs | `_roothaan_J0(λ_e, λ_n)` | Swapped | Symmetric? |
|---|---|---|---|
| Electronic: λ_e=1.0, λ_n=85.7 | 0.9997354669 | 0.9997354669 | ✓ to 1e-15 |
| Muonic: λ_lepton=185.85, λ_n=85.7 | 71.3223070476 | 71.3223070476 | ✓ to 1e-15 |

This is the closed form, exact and finite in **both** regimes. There is no "regime to dispatch on" at the level of the production output.

### Finding 2: Production code uses the closed form directly

The cross-register recoil pipeline:
- `cross_register_recoil_correction(spec, m_e_over_m_n)` — line 894 of `cross_register_vne.py`
- → `_roothaan_J0(spec.lam_e, spec.lam_n)` — line 946 (the closed form)

No asymptotic series is invoked on this path. `roothaan_J0_taylor_expansion` (line 1238) and `roothaan_recoil_shift_through_order` (line 1363) are diagnostic / verification utilities called only from `tests/test_cross_register_vne.py` (lines 546, 553, 566) and from `pachucki_higher_order_comparison`. No production observable depends on the asymptotic expansion.

Track B uses `_roothaan_J0` via `cross_register_recoil_correction`, with `lam_e` set to either 1.0 (electronic) or `M_RED_MUP` ≈ 185.85 (muonic) and `lam_n = LAM_NUCLEUS_DEFAULT` ≈ 85.7. Both regimes flow through the same closed-form code path. The "λ_μ > λ_n breaks the large-nucleus expansion" diagnosis was based on the asymptotic-series view, not on the production-path computation.

### Finding 3: Track A's SE gap is in a different module entirely

`debug/sprint_mh_track_a.py` line 168 defines `self_energy_eides_lepton`, which encodes the leading-order Eides Sec.3.2 SE bracket:

```python
common_dim = ALPHA**3 * Z**4 / (math.pi * n**3)
ha_lepton_meV = m_red_in_me * HA_TO_MEV
common_meV = common_dim * ha_lepton_meV
if l == 0:
    bracket = (4/3) ln(1/(Zα)²) - (4/3) ln_k0_overRy + 10/9
elif l == 1 and n == 2:
    bracket = -(4/3) ln_k0_overRy - 1/6
return common_meV * bracket
```

This formula:
- does not call `cross_register_vne.py`;
- scales linearly in m_red (via the `ha_lepton_meV` factor);
- contains no next-to-leading α(Zα)⁴(m_red/m_p) recoil-mixing terms.

Antognini 2013's canonical muon SE of −0.668 meV includes Salpeter recoil + Pachucki recoil-SE corrections (Eides Tab. 4.2, Pachucki 1996 Eq. 25–30). Those terms are **field-theoretic vertex corrections** (Bodwin–Yennie class) that Track A's leading-order Eides bracket structurally omits.

**The 0.16 meV gap is structurally outside the cross-register multi-focal kernel.** Adding a Pachucki dual expansion to `cross_register_vne.py` cannot close it — the gap doesn't live in the J_0 integral. It lives in the SE vertex topology, on the LS-8a wall.

## What was added in this sprint

Even though regime dispatch is moot for production, the directive's request for a dual symbolic expansion is a useful diagnostic addition: it parallels the existing `roothaan_J0_taylor_expansion` utility and makes the regime symmetry explicit in code. Two new functions and a module-level commentary block were added:

### `roothaan_J0_taylor_expansion_dual(n_terms, lam_n_value)`

Symbolic Taylor expansion of J_0 in eps' = 1/lam_e around the heavy-lepton limit lam_e → ∞ (the muonic regime). By the symmetry of J_0, this is the (lam_e ↔ lam_n) swap of the existing electronic-regime expansion:

```
J_0(1/eps', lam_n) = lam_n
                     − 2 lam_n^3 (eps')²
                     + 5 lam_n^4 (eps')³
                     − 9 lam_n^5 (eps')⁴
                     + 14 lam_n^6 (eps')⁵
                     − 20 lam_n^7 (eps')⁶ + …
```

The recoil shift `lam_n − J_0` has coefficients +2 lam_n³, −5 lam_n⁴, +9 lam_n⁵, … — same magnitudes as the electronic series, with lam_n replacing lam_e.

### `roothaan_recoil_shift_regime_aware(lam_lepton, lam_nucleus, Z, max_order)`

Auto-dispatch wrapper:
- **Electronic regime** (lam_lepton ≤ lam_nucleus): forwards to existing `roothaan_recoil_shift_through_order` (large-nucleus expansion in eps = 1/lam_n).
- **Muonic regime** (lam_lepton > lam_nucleus): uses the same numerical machinery with the labels swapped, and reports the natural muonic recoil estimator +Z(lam_lepton − J_0).

Both branches return the closed-form J_0 alongside the truncated series for cross-checking. The dispatcher labels the regime explicitly in its return value.

### Module-level commentary block

A new section headed `Sprint MH Track D: regime-aware dual expansion (May 2026)` documents:
- The structural finding that the closed-form `_roothaan_J0` is regime-agnostic and used directly in production;
- That the asymptotic-series utilities are diagnostic, not on the production path;
- Where Track A's SE gap actually lives (Eides bracket in `sprint_mh_track_a.py`), and why closing it requires the LS-8a-renorm sprint, not a Roothaan kernel adjustment.

This block is the durable institutional-memory artifact — future PMs reading `cross_register_vne.py` will not re-derive the regime-dispatch question.

## Validation

### Existing tests: bit-identical regression

```
tests/test_cross_register_vne.py: 56 baseline tests pass (was 56, still 56).
```

The 56 baseline tests cover the closed-form `_roothaan_J0`, the existing
electronic-regime Taylor expansion, the integer-only sub-sum, and the
Pachucki higher-order comparison. All 56 pass without modification.

### New tests: 7 added in `TestSprintMHTrackDDualExpansion`

```
test_J0_closed_form_is_symmetric_both_regimes              PASS
test_dual_expansion_symbolic_structure                      PASS
test_dual_expansion_evaluated_at_lam_n_one_matches_original PASS
test_regime_aware_electronic_matches_existing               PASS
test_regime_aware_muonic_dispatches_correctly               PASS
test_regime_aware_at_boundary_lam_lepton_equals_lam_nucleus PASS
test_dual_expansion_evaluated_at_specific_lam_n             PASS
```

Total: **63/63 pass** (56 baseline + 7 new).

### Numerical bit-identity on Track B

Re-running `debug/sprint_mh_track_b.py::predict_muonic_hyperfine()` after the changes:

| Component | Old value | New value | Bit-identical? |
|---|---|---|---|
| ep recoil estimate (Ha) | 2.645287836853250e-04 | 2.645287836853250e-04 | ✓ |
| mup recoil estimate (Ha) | 1.145189677907791e+02 | 1.145189677907791e+02 | ✓ |

Track A's Lamb shift was not re-run because `cross_register_vne.py` is not on Track A's call graph (the SE formula doesn't use it). Track A's full Uehling integration, Eides SE bracket, and Friar moment all live in `sprint_mh_track_a.py` and are untouched.

**Combined production observables unchanged:**
- Track A muonic Lamb shift: framework + lit total = +202.17 meV vs CREMA 202.37 meV → −0.10% (unchanged);
- Track A muonic Lamb shift framework-native = +200.50 meV → −0.92% (unchanged);
- Track A Uehling = +205.0074 meV vs Antognini 2013 +205.0074 meV → < 1 ppm (unchanged);
- Track B muonic Bohr–Fermi = 182.4433 meV vs Eides QED-only 182.443 meV → +2 ppm (unchanged);
- Track B Bethe–Salpeter electronic recoil match = 2.86% (unchanged).

## What would actually close Track A's 24% SE gap

Two routes, both outside the scope of `cross_register_vne.py`:

### Route 1: Extend `self_energy_eides_lepton` with literature inputs (Track A patch)

Add the next-to-leading α(Zα)⁴(m_red/m_p) recoil-mixing terms from Eides Tab. 4.2 / Pachucki 1996 Eq. 25–30 directly to the SE bracket. This is a **literature-input** patch (not framework-native): the m_red-dependent corrections are tabulated in published QED literature, and adding them brings Track A's SE from 0.83 → 0.668 meV (closes the 0.16 meV gap, matching Antognini at the printed precision).

This would close the framework-native total from 200.50 to ~200.66 meV, bringing the framework-native residual from −0.92% to ~−0.85% against CREMA. The combined framework + literature total would remain −0.10%, since the Track A literature inputs already include "recoil corrections" but use Antognini's central muon SE (so the gap closes only on the framework-native side, not on the combined total).

This is structurally a **paper edit** to Paper 36 §VIII.2 (Sprint MH Track A subsection): replace the Track A self-energy entry with the literature-supplemented version, and document the additional terms as inputs (Bodwin–Yennie + Pachucki recoil-SE, on the LS-8a wall).

This is straightforward and uncontroversial, but it is **not what this Track D was scoped for.**

### Route 2: LS-8a-renorm sprint (framework-native)

Implement field-theoretic vertex renormalization at the spectral level so the framework autonomously generates the next-to-leading recoil-SE terms. This is the LS-8a-renorm sprint flagged in CLAUDE.md §2 as deferred ("multi-sprint scope"). It would close the gap framework-natively, but it requires importing flat-space Z₂/δm conventions into the spectral-action machinery — explicitly deferred per the May 7 strategic decision ("would contaminate structural-skeleton purity").

Neither route is a Roothaan-kernel adjustment, so neither is reachable from `cross_register_vne.py`.

## Honest framing

The original directive's framing — "Pachucki dual expansion in cross_register_vne.py for the regime λ_lepton > λ_nucleus" — read as an extension of the May 7 D-Pachucki-higher-order sprint, which had specifically diagnosed sub-leading Roothaan terms as half-integer-aliased basis-truncation artifacts. The natural inference was that the muonic regime might have a structurally different sub-leading content recoverable via a dual expansion.

What's actually true: the half-integer aliasing was driven by the choice of focal-length calibration λ_n = 2√M_p (which makes 1/λ_n ~ √(m_e/m_p)). That choice is symmetric in regime — it's a property of the basis, not of the regime. Swapping to the dual expansion doesn't change the aliasing mechanism; it just relabels which axis is the asymptotic one. The sub-leading sign mismatch noted in the May 7 memo persists in either regime and would only close at multi-shell n_max ≥ 2.

Track A's 24% SE gap was attributed in its own memo to "leading-order m_red scaling that omits next-to-leading α(Zα)⁴ × (m_red/m_p) recoil-mixing terms." That attribution is correct; the omitted terms are field-theoretic vertex corrections (Bodwin–Yennie class), not asymptotic-series tail terms in J_0. The Roothaan kernel is the recoil **numerator**, not the SE vertex correction. Closing the SE gap requires reaching the vertex topology, which is the LS-8a wall.

This Track D therefore returns a clean structural negative on the directive's primary premise, plus a modest documentation-and-diagnostic addition that makes the regime symmetry explicit in code and prevents future re-derivation of the same diagnosis.

## Sprint provenance

- Parent context: Sprint MH Track A (debug/sprint_mh_track_a.py, debug/sprint_mh_track_a_memo.md) — Lamb shift at −0.10%, SE gap 0.16 meV;
- Parent context: Sprint MH Track B (debug/sprint_mh_track_b.py) — flagged the regime question;
- Antecedent: Phase C-Pachucki-higher-order sprint (cross_register_vne.py lines 1208–1336) — characterized half-integer aliasing in the electronic-regime asymptotic series;
- Antecedent: Phase C-W1a-physics (debug/multifocal_phase_c_w1a_physics_memo.md) — Roothaan 1951 cross-register kernel.

## Files produced

- `geovac/cross_register_vne.py` — appended ~250 lines: module commentary block + `roothaan_J0_taylor_expansion_dual` + `roothaan_recoil_shift_regime_aware`. No existing function modified. Backward compatibility preserved by construction.
- `tests/test_cross_register_vne.py` — appended `TestSprintMHTrackDDualExpansion` class with 7 tests covering the new utilities and the regime-symmetry finding.
- `debug/sprint_mh_track_d_memo.md` — this file.

No paper edits applied. No observable improvement. Track A's 24% SE gap is not closed by this sprint and cannot be closed by work in `cross_register_vne.py`.
