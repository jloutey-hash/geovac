# Sprint TD Track 4 — master Mellin engine M1 extensibility test on Euclidean Schwarzschild

**Date:** 2026-05-09
**Track:** Sprint TD Track 4
**Builds on:** Track 1 (commit 697666f, `geovac/thermal_tensor_triple.py`) — the explicit `T_{S^3} ⊗ T_{S^1_β}` tensor-product spectral triple that reproduces Stefan–Boltzmann as `M1 × M2` factorisation.
**Status:** **POSITIVE-LIMITED**, with two flagged scope boundaries (M2 deferred, Bertrand negative).

---

## §1. Framing — what this track tests, and what it does not

This track tests **one** narrow question: does the master Mellin engine's M1 mechanism (Hopf-base measure / Matsubara temperature from S¹ compactification) verified for Stefan–Boltzmann on `S^3 × S^1_β` (Track 1) extend cleanly to the Euclidean Schwarzschild cigar geometry, reproducing the Hawking temperature `T_H = 1/(8π M)`?

This is a **TEMPERATURE-MECHANISM test, not an entropy-shape test or a black-hole-derivation claim.** Per CLAUDE.md §1.5 rhetoric and the directive's honest-scope reminder, three things are explicitly **NOT** in scope:

1. **Deriving black hole thermodynamics from the GeoVac packing axiom.** Bertrand's theorem (closed-orbit selection of `−Z/r` among central potentials) is what selects `S^3` in GeoVac. There is no analog selecting the cigar; the cigar is forced by Einstein equations + Schwarzschild horizon at `r = 2M`. The mass `M` enters as external input.
2. **Deriving `S_BH = A/(4 ℓ_p²)` from the heat kernel.** That requires Connes–Chamseddine spectral action evaluated at the Schwarzschild saddle (`M2 / Seeley–DeWitt` content on the horizon `S²`), and is well-known multi-week NCG work. It is the natural follow-up but is explicitly out of scope for this track.
3. **Comparing GeoVac entropy to BH entropy structurally.** Track 5 (PSLQ probe of Paper 27 information entropy against the master Mellin ring) addresses whether GeoVac entropy is itself spectrally expressible, and is the cleaner setting for the entropy-shape comparison the user raised.

The user's earlier conversation (turn covering "GeoVac in S^3 vs BH in R^3") established that the entropy SHAPES of GeoVac (extensive-correlational on a compact ambient) and Schwarzschild (intensive-holographic on a compact submanifold of a non-compact ambient) are genuinely different. This track stays inside the temperature mechanism, where the Matsubara structure IS shared.

---

## §2. Setup — Euclidean Schwarzschild as a tensor product

Euclidean Schwarzschild has metric

$$
\mathrm{d}s^2 = \left(1 - \frac{2M}{r}\right) \mathrm{d}\tau^2 + \left(1 - \frac{2M}{r}\right)^{-1} \mathrm{d}r^2 + r^2 \, \mathrm{d}\Omega_2^2.
$$

Smoothness at the horizon `r = 2M` (no conical singularity at the cigar tip) is the standard Hawking–Gibbons argument. Set `ρ² = 4(2M)(r − 2M)`; then near the horizon `ds² ≈ ρ² dφ² + dρ² + (2M)² dΩ²` with `φ = τ/(4M)`. Smoothness at `ρ = 0` requires `φ` to have period `2π` (else conical defect). Hence the τ-circle has period

$$
\beta_{\rm cigar} = 4M \cdot 2\pi = 8\pi M, \qquad T_H = \frac{1}{\beta_{\rm cigar}} = \frac{1}{8\pi M}.
$$

Equivalently `T_H = κ_g/(2π)` with surface gravity `κ_g = 1/(4M)`.

The cigar factors **locally near the horizon** as `(ρ-φ disk) × S²_horizon`, and **globally** the spatial slice is non-compact in `r ∈ [2M, ∞)`. So strictly the geometry is not a compact-×-compact tensor product. This non-compactness is the principal scope difference from Track 1.

For this track, we keep attention on the **τ-circle factor** alone, which IS compact and identical in structure to Track 1's `S^1_β`.

---

## §3. The M1 mechanism applied to the τ-circle — does it reproduce 2π in T_H?

**Yes, bit-identically, with no new code path.**

`debug/sprint_td_track4.py` calls Track 1's `geovac/thermal_tensor_triple.matsubara_spectrum()` with `β = 8π M` and obtains the cigar Matsubara spectrum verbatim:

| Frequency | Value | Identification |
|:---|:---|:---|
| Bosonic lowest nonzero `ω_1 = 2π/β` | `1/(4M)` | **= surface gravity κ_g** (sympy residual = 0) |
| Fermionic lowest `ω_0 = π/β` | `1/(8M)` | **= π · T_H** (sympy residual = 0) |
| `T_H` from `β_cigar = 8π M` | `1/(8π M)` | matches `κ_g/(2π)` (sympy residual = 0) |

The single `2π` in `T_H = 1/(8π M)` is the M1 / Hopf-base measure signature of the τ-circle's circumference. This is precisely the same `2π` that appears in Track 1's `ω_k = 2π k/β` Matsubara modes on `S^1_β`, and was tagged in Track 1 as `PI_TAG_M1_MATSUBARA_CIRCLE = "M1: 2 pi from imaginary-time circle circumference (S^1_beta factor)"`.

The Track 1 K_temporal leading kernel `β/(2√(πs))` from Jacobi-θ inversion transports verbatim. K_spatial in Track 1 was on `S^3`; on the cigar K_spatial sits on the `(r, Ω)` non-compact-with-boundary geometry, and its full evaluation needs IR regularization at `r → ∞` — that part is out of scope. The τ-circle factor itself is unchanged.

**Net for §3:** the formal M1 mechanism extends cleanly. The 2π in T_H IS an M1 output.

---

## §4. Verdict

Three sub-verdicts (full text in `debug/data/sprint_td_track4.json`):

**V1. M1 mechanism (Matsubara temperature from S¹) — POSITIVE.**
The `2π` in `T_H = 1/(8π M)` is the same Hopf-base-measure / `M1` signature as in Track 1's Stefan–Boltzmann. Track 1's `matsubara_spectrum()` with `β = 8π M` reproduces the standard Hawking spectrum bit-identically. The lowest bosonic mode equals surface gravity; the lowest fermionic mode equals `π T_H`. All three sympy residuals are 0.

**V2. M2 mechanism (Seeley–DeWitt on horizon `S²`) — DEFERRED, OUT OF SCOPE.**
Recovering `S_BH = A/(4 ℓ_p²)` from the horizon `S²` heat kernel is the Connes–Chamseddine spectral action computation at the Schwarzschild saddle. Multi-week NCG work, explicitly out of scope per directive. Flagged as the natural follow-up.

**V3. Bertrand-class forcing of cigar geometry — NEGATIVE.**
The cigar is forced by Einstein equations + Schwarzschild horizon, NOT by GeoVac's packing axiom. There is no Bertrand-style closed-orbit theorem on the gravity side that singles out the cigar. The mass `M` is external input, set by the gravitational source. GeoVac does not "explain" black holes.

**Net verdict:** POSITIVE for the limited M1 extensibility claim; OUT-OF-SCOPE-DEFERRED for the M2 / `S_BH` claim; NEGATIVE for any claim that GeoVac autonomously produces the cigar geometry. Reading (a) of the directive's verdict structure, with the (a)(b)(c) obstruction list partially confirmed — **the cigar's non-compact r-direction (a) does sit outside the master Mellin engine's compact-spectrum scope**, but the τ-circle factor is unaffected, so M1 itself extends.

---

## §5. Honest scope statement

What this track shows:
- Track 1's `matsubara_spectrum()` function works verbatim with `β = 8π M`.
- The 2π in T_H is the same M1 mechanism as the 2π in Track 1's S¹_β Matsubara modes.
- Surface gravity, lowest fermionic mode, and `T_H = κ_g/(2π)` factorisation all reduce to identities verified symbolically (sympy residuals 0).
- Paper 35 Prediction 1 (π enters iff continuous integration over a temporal/spectral parameter) holds verbatim for the cigar — the single π in `T_H` comes from the τ-circle Matsubara sum.
- Paper 32 §VIII case-exhaustion theorem (M1/M2/M3 sub-cases of master Mellin engine) extends formally to the cigar at `k = 0` (M1).

What this track does NOT show:
- The cigar's non-compact r-direction was not handled — only the τ-circle factor.
- The horizon S² was not analysed at all (M2-class Seeley–DeWitt deferred).
- No statement is made about GeoVac entropy vs Bekenstein–Hawking entropy. Those have different shapes; Track 5 is the right venue for the comparison.
- M is external input; no GeoVac-internal mechanism selects it.

---

## §6. Implications for framework scope

This track is a small but real result. The structural-skeleton-scope language from CLAUDE.md §2 (Sprint of 2026-05-07) said GeoVac "determines selection rules, transcendental signatures, scaling laws, divergence structure, factorization theorems, and upper bounds." The narrow finding here is: **the M1 transcendental signature (2π from S¹ circumference) is universal across S¹ compactifications regardless of the host spacetime.** Coulomb-derived `S^1_β` (Track 1) and gravitational cigar (this track) both produce the same `2π` from the same Hopf-base mechanism.

This is consistent with Paper 32 §VIII's case-exhaustion theorem at `k = 0`: the M1 mechanism `Tr(D^0 e^{−tD^2}) = Vol/something` is the most kinematic-only of the three Mellin sub-mechanisms, and is the most likely to extend across geometries. M2 (Seeley–DeWitt on the horizon, would-be `S_BH`) is geometry-specific and requires per-system work; M3 (vertex parity / Hurwitz / Catalan G) is observable-restricted.

The right reading after this track:
- M1 extends across S¹ compactifications (this track).
- M2 needs case-by-case work for non-compact / non-Coulomb manifolds (out of scope).
- M3 is observable-restricted (Paper 32 §VIII Remark).

This is not a black-hole-physics claim. It is a **claim about the M1 sub-mechanism of the master Mellin engine**: namely, that it transports across S¹ compactifications without modification, exactly as the formal apparatus would predict.

---

## Files

- Driver: `debug/sprint_td_track4.py`
- Data: `debug/data/sprint_td_track4.json`
- Memo: this file (`debug/sprint_td_track4_memo.md`)

No production `geovac/` code modified. All sympy-symbolic checks pass (3 nontrivial residuals each = 0).
