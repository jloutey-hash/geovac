# Gravity Campaign Phase 1 — Möbius α>1 sign-discrepancy diagnostic

**Date:** 2026-05-29
**Path:** Gravity Campaign Phase 1, diagnostic (recompute, not build). Curve-fit-audit discipline applied.
**Verdict:** **(A) ARTIFACT — of observable definition, not of driver/substrate.** The v3.19.0 −0.056 and the thread-9 +0.052 are two *different observables* of the *same* substrate `Δ_K(2) = +0.0843 > 0`. There is no real disagreement and no substrate artifact. The two signs are each independently consistent with the continuum Sommerfeld–Cheeger sign of their own observable. **Recommendation: scope-limited Paper 51 closure** (the substrate-level Möbius identification already stands; multi-axis A-extraction is not justified by this finding).

This diagnostic independently re-derives, from a clean single-convention recompute, the conclusion already reached in `debug/sprint_moebius_convention_audit_resolution_memo.md` (same day). It hardens that resolution with the explicit continuum sanity-anchor (Task 4) that the earlier audit did not run.

## 1. The two drivers and their EXACT convention difference (Tasks 1 & 3)

| | v3.19.0 driver | thread-9 driver |
|---|---|---|
| File | `debug/alpha_gt_1_ansatz_test.py` | `debug/sprint_moebius_reading_b_nrho_sweep.py` |
| docstring def | "slope = Delta_K / (1/alpha − alpha)" (line 5–6) | "slope = d(Delta_K)/d(alpha) … via central FD" (line 19, code line 79) |
| **observable** | **RATIO** `Δ_K(α) / (1/α − α)` | **DERIVATIVE** `dΔ_K/dα` |
| α=2 value | **−0.0562** | **+0.052** |

Everything else is identical between the two:
- **Same tip / bulk subtraction:** both use `Δ_K(α) = K_wedge(α) − α·K_disk`.
- **Same apex-angle convention:** both use `N_φ(α) = round(α·N_0)`, `m_eff = (k+½)/α` (anti-periodic spinor, half-integer azimuthal index).
- **Same baseline:** `a=0.05, N_0=120, N_ρ=200, t=1.0`.
- **Same substrate class:** FD vs spectral azimuthal gives the *same* number (verified below).

The pinned source of the sign flip is therefore **neither** the differentiation direction, **nor** the bulk-subtraction term, **nor** the apex-angle convention, **nor** FD-vs-spectral. It is **the choice of observable: ratio vs derivative.**

## 2. Clean recompute under ONE convention (spectral substrate, Task 2)

Driver: `debug/gravity_campaign_phase1_moebius_recompute.py`
Data: `debug/data/gravity_campaign_phase1_moebius_recompute.json`

`K_disk(t=1) = 41.6893`. At `a=0.05, N_0=120, N_ρ=200, t=1.0` on the **spectral** azimuthal substrate (G4-6d, `m_eff=(k+½)/α`):

| α | Δ_K | (1/α − α) | RATIO (v3.19 obs) | DERIV (thr9 obs) |
|---|---|---|---|---|
| 1.5 | +0.05388 | −0.8333 | **−0.06466** | **+0.0757** |
| 2.0 | +0.08433 | −1.5000 | **−0.05622** | **+0.0513** |
| 3.0 | +0.13020 | −2.6667 | **−0.04883** | **+0.0438** |

- The RATIO column reproduces v3.19.0's reported `{−0.06466, −0.05622, −0.04883}` **bit-exactly to 4 dp**.
- The DERIV column reproduces thread-9's `+0.0513` at α=2 **bit-exactly**.

Both come from the **single positive** quantity `Δ_K(2) = +0.0843`.

**Substrate-discretization invariance (FD vs spectral), α=2:**
- spectral RATIO = −0.056222
- FD RATIO = −0.056225
- difference 2.2×10⁻⁶ — confirms thread 9's "substrate-discretization-stable" finding (it was only the *observable definition*, not the substrate, that differed between the drivers).

## 3. The crisp α=2 statement (the resolution)

```
Δ_K(2)              = +0.0843   (POSITIVE — excess-angle tip, physically correct for a saddle cone)
(1/α − α) at α=2    = −1.5      (NEGATIVE for any α>1)
RATIO  = Δ_K/(1/α−α) = −0.0562  (v3.19.0)  — NEGATIVE because positive/negative
DERIV  = dΔ_K/dα     = +0.0513  (thread 9) — POSITIVE because Δ_K increases with α
```

Same substrate, same `Δ_K(2) > 0`. The minus sign in v3.19.0 is the kinematic sign of `(1/α − α) < 0` in the denominator of the *recovery ratio*; the plus sign in thread 9 is the genuine slope of `Δ_K` vs α. Neither is wrong.

## 4. Continuum Sommerfeld–Cheeger anchor (Task 4)

Pure scalar SC tip (deficit-angle form, analytically continued): `Δ_K^cont(α) = −(1/12)(1/α − α)`.

- **SC RATIO** = `Δ_K^cont/(1/α−α)` = **−1/12 = −0.08333** (α-independent).
- **SC DERIVATIVE** = `d/dα[−(1/12)(1/α−α)]` = `(1/12)(1/α² + 1)`; at α=2 = **5/48 = +0.10417**.

| observable | discrete sign (α=2) | continuum SC sign | consistent? |
|---|---|---|---|
| RATIO (v3.19.0) | −0.0562 | −1/12 = −0.0833 | **YES** (both negative) |
| DERIVATIVE (thread 9) | +0.0513 | 5/48 = +0.10417 | **YES** (both positive) |

**Both discrete signs match the continuum sign of their own observable.** This is the decisive sanity check: there is no sign that is "wrong relative to the continuum." The substrate sits *below* the pure SC value in both observables — RATIO at −0.0562 vs −0.0833 (the empirical Möbius suppression, matching `−1/12·α/(2α−1) = −0.0556` at 1.2%), and DERIVATIVE at +0.0513 vs +0.10417 (≈49% of the continuum SC derivative). The suppression is the *same physical Möbius deficit* viewed through two observables.

## 5. Curve-fit-audit pass (per memory/feedback_audit_numerical_claims.md)

The original "sign discrepancy" is itself a numerical claim about the substrate; the audit:

1. **Free-param count.** Both observables are parameter-free recomputes from the substrate. The match of the recompute to *both* prior reported values (RATIO bit-exact, DERIV bit-exact) leaves zero fit freedom. Pass.
2. **Selection bias.** No selection — both prior numbers reproduced under one driver. Pass.
3. **Alternatives.** Four candidate sources of a *real* sign flip (differentiation direction, bulk subtraction, apex-angle convention, FD-vs-spectral) were each tested and excluded; the only remaining difference is the ratio-vs-derivative observable definition. Pass.
4. **Robustness.** The ratio observable is N_ρ-stable (thread 9: CV 3.4%) and FD-vs-spectral-stable (this memo: 2×10⁻⁶). Pass.
5. **Independent test.** The continuum SC anchor (Task 4) — not used in either prior driver — independently confirms both signs are continuum-consistent. Pass.

The audit lands the same way as the Fursaev–Solodukhin grounding task #26 did for the *mechanism*: the **empirical observation is real and reproducible**; what needed resolving was an attribution / interpretation, here the observable definition behind the "sign discrepancy." (The *mechanism* of the Möbius factor — `−(1/12)·α/(2α−1)` — remains an OPEN structural question per task #26; this diagnostic does not reopen or close it.)

## 6. Verdict and downstream fork recommendation

**Verdict: (A) ARTIFACT of observable definition.** Resolved correct signs:
- **RATIO** `Δ_K/(1/α−α)` is **negative** at α>1 (v3.19.0 convention, the recovery quantity, continuum-anchored to −1/12).
- **DERIVATIVE** `dΔ_K/dα` is **positive** at α>1 (thread-9 convention, the slope, continuum-anchored to (1/12)(1/α²+1)).
- Pinned source: ratio-vs-derivative, both built on the single positive `Δ_K(2)=+0.0843`. No driver bug, no substrate artifact, no convention error in either driver.

**Downstream fork — recommend scope-limited Paper 51 closure (NOT multi-axis A-extraction).** The sign "discrepancy" was the load-bearing flag that gated G4-6a-refined (per `g4_6a_refined_multi_axis_scoping_memo.md` Option 2 → Option 1). It is now resolved as a false alarm: the substrate-level Möbius RATIO identification is real, reproducible, substrate-discretization-invariant, and continuum-sign-consistent. The 3–6 month multi-axis A-coefficient extraction (Option 1) was justified *only* if the Möbius foundation were undermined; it is not. The separate, genuine limitation — that the **A coefficient** (the leading `1/(24π)` Weyl term, a *different* quantity from the conical tip slope) does not converge under single-axis N_φ refinement (thread 8 clean negative) — stands independently and is not affected by this resolution. Paper 51 already carries the substrate-discretization-invariance verification (per the same-day convention-audit memo); the remaining honest scope note for Paper 51 is that the Möbius *mechanism* is open (task #26), not its sign. No further G4-6a-refined commitment is warranted by the sign question.

## 7. Honest scope

This diagnostic:
- **Recomputes** both observables under one convention on the spectral substrate; reproduces both prior numbers bit-exactly.
- **Pins** the sign flip to ratio-vs-derivative; excludes the four other candidate causes.
- **Anchors** both signs against the continuum SC ratio (−1/12) and derivative (5/48); both consistent.
- **Confirms** FD-vs-spectral substrate invariance at 2×10⁻⁶.

Does NOT:
- Edit papers or production `geovac/` code (diagnostic only).
- Reopen or close the Möbius *mechanism* question (open per task #26 — separate from the sign).
- Address the A-coefficient (Weyl-term) extraction limitation (thread-8 N_φ clean negative, separate and unaffected).

## 8. Files

- `debug/gravity_campaign_phase1_moebius_sign_memo.md` (this)
- `debug/gravity_campaign_phase1_moebius_recompute.py` (recompute driver, spectral substrate, both observables + continuum anchor)
- `debug/data/gravity_campaign_phase1_moebius_recompute.json` (data)

## 9. Cross-references

- `debug/sprint_moebius_convention_audit_resolution_memo.md` — same-day audit reaching the same conclusion (this memo hardens it with the explicit continuum anchor).
- `debug/sprint_moebius_reading_b_nrho_sweep.py` / `_memo.md` — thread 9 (the derivative observable, +0.052).
- `debug/alpha_gt_1_ansatz_test.py` — v3.19.0 (the ratio observable, −0.0562; source of truth on the ratio convention).
- `debug/fursaev_solodukhin_1995_grounding_memo.md` — task #26 (Möbius *mechanism* open; citation falsified). Distinct from the *sign* question resolved here.
- `debug/g4_6a_refined_multi_axis_scoping_memo.md` — the scoping memo whose Option-2 gate this diagnostic clears toward scope-limited closure.
