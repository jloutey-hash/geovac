# Sprint WH7 Toeplitz probe — Step 1 (2026-06-10, v3.111.0)

**Goal:** first falsifier probe for WH7 (time-discreteness is observer-compactification,
registered same day): is time metrically visible to the translation seminorm when the
temporal algebra is built from genuine time-dependent multipliers (Toeplitz-compressed
e^{iωt}, the Connes–vS S¹ pattern) instead of the P45 momentum-diagonal functions g(D_t)?

**Verdict: POSITIVE-REBUILD, with the primary falsifier leaning "weakens-to-convention"
on the visibility leg.** Five checks, all PASS
(`debug/wh7_toeplitz_temporal_probe.py`, JSON in `debug/data/wh7_toeplitz_probe.json`,
frozen falsifier `tests/test_wh7_toeplitz_temporal.py`):

| Check | Result |
|:--|:--|
| A single-mode exactness | L(S_q) = 2πq/T = Lip(e_q) at max err 1.4×10⁻¹⁴ over T ∈ {1, 2π, 10} × K ∈ {4,8,16} × q ∈ {1..4} — the translation seminorm IS the continuum Lipschitz constant on temporal modes, exactly, at every finite cutoff |
| B P45 control | random non-constant g(D_t): L = 2.8×10⁻¹⁴ — the celebrated "Lipschitz invisibility of time" reproduced as a property of the momentum-diagonal ALGEBRA, not of time |
| C kernel condition | L = 0 ⟺ f constant on the band-limited algebra (q_max ≤ K); 20-random-panel min nonconst L = 1.52, const L = 5.7×10⁻¹⁴ |
| D Lipschitz domination | L(F) ≤ Lip(f) always; window ratio 0.752 → 0.995 monotone as K = 3 → 32 (Fejér-style convergence, same shape as the P38 spatial story) |
| E de-compactification | at fixed physical frequency ω = 2π: L = ω to 8.9×10⁻¹⁶ for T ∈ {1..16} — metric visibility does NOT degrade as the carrier de-compactifies at fixed bandwidth |

**Structural content.** The s→0 limit of the translation seminorm is the D_t-commutator
norm ‖[D_t, F]‖, and [D_t, S_q] = ω_q S_q ≠ 0, whereas [D_t, g(D_t)] = 0 identically.
So the P45 annihilation (and its "L3 structural identity") was an architecture artifact:
the v1 temporal algebra was built from functions OF the frequency operator, which no
translation can move. With the honest Toeplitz algebra, time is metrically visible on
the compact carrier with exact mode-level Lipschitz recovery — the compact-time wing of
the Lorentzian program is rebuilt at the operator level, and the joint S³×S¹ statement
is now a P38-machinery exercise (named follow-on B1).

**WH7 impact (status updated in CLAUDE.md §1.7 per rule 9).** Leg (ii) of WH7's inputs
("non-compact time is Lipschitz-invisible") is DOWNGRADED to "the momentum-diagonal
algebra is Lipschitz-invisible"; check E shows visibility survives de-compactification
at fixed bandwidth, so the primary falsifier currently leans toward the
"weakens-to-convention" branch for metric VISIBILITY. What it does NOT touch: leg (i) —
discreteness (the Matsubara ladder) and the π-injection appear exactly at
compactification (Paper 35), which is WH7's load-bearing content; and leg (iii)
spectral carrier-indistinguishability. WH7 sharpens to: *time is metrically visible
but only compactification makes it discrete — and the observer's window is what
compactifies.*

**Honest scope.** Step 1 is compact-carrier + fixed-bandwidth-scaling only. Step 2
(genuine non-compact construction: Paley–Wiener band-limit on ℝ_t inside the
Latrémolière pointed-proper framework, properness/total-boundedness bookkeeping) is the
remaining falsifier work — multi-week, named follow-on B2. The kernel condition holds
for the band-limited algebra within window reach (q_max ≤ K), mirroring Paper 38's
band N ≤ 2n_max − 1 condition.

**Follow-ons.**
- B1: joint S³ × S¹_T translation-seminorm GH statement (spatial P38 + temporal Toeplitz;
  rebuilds the descoped Paper 45 ambition on sound footing). Sprint-scale.
- B2: non-compact Step 2 (pointed-proper, PW_Λ(ℝ)). Multi-week; decides WH7's primary
  falsifier branch.
