# Sprint: quadrature-free diatomic + certified reference artifact + two Hamiltonian probes + u1/u2 — 2026-08-21

**Origin:** PI approved the "sell certainty" plan (QFD + certified table now, isoenergetic circuits later)
and the two probe killer-tests; all parallel. Canonical memo; per-track detail:
`debug/qfd_track1_findings.md`, `benchmarks/certified_reference/`, `debug/probeA_findings.md`,
`debug/probeB_findings.md`, `debug/beta2_u1_README.md`.

## Verdicts

| Track | Verdict |
|:--|:--|
| T1 QFD | **GO** — H₂ closed-form, truncation-free, **84 digits**; LiH zero-quadrature, **30 digits net of explicit τ-tail bound** |
| T2a Certified artifact | **GO** — 51 entries, all evidence-carrying; deterministic generator; 10 tests green |
| Probe A (spectrum-aware QSVT) | **CONSTANT-FACTOR** 8–10×, provably no more; the PSD-shift √κ variant is the one scoped follow-on (amplification-gated) |
| Probe B (analytic "THC") | **MIXED → mechanism dead** — the DF control matches it at rank ≤ n_orb² with ~10³× fewer leaves; not THC (full-rank leaves) |
| u1 (120-digit T2) | **DONE**: confirms the 66 certified digits exactly, refutes hi1's unclaimed 67+, extends single-method to ~139; **u2 cross-validation running** |

## T1 — the quadrature-free diatomic (promoted: geovac/qfd_core.py + qfd_assemble.py; tests/test_paper58_qfd.py)

- H₂ (1s/1s, ζ=1, R=1.4): E = −1.10655660609135850801940122475993772281832890689690719549979 Ha,
  84 digits (dps 60/90 agree 6e-85); R=1.6, 2.0, ζ=1.197 alongside. Homonuclear exchange τ-series
  TERMINATES symbolically (τ=1, τ>2 exact zeros) — no truncation anywhere. Exchange vs **Sugiura
  1927** 3.1e-61 (independent {E₁, ln, γ} corroboration); all integrals vs quadrature 1e-21..23.
- LiH (Li 1s,2s + H 1s, R=3.015): E_tot = −7.87261979241561217317558085835659691730258908 Ha,
  30 digits net of the τ-tail bound (measured monotone ratios → Σ|tail| ≤ 1.95e-32, amplified ×64
  through 2-RDM + ‖S^{-1/2}‖⁴ → 1.24e-30 Ha). Heteronuclear τ-series does NOT terminate.
- NEW: the one-electron two-center closed forms did not exist and were derived (Mulliken A/B
  auxiliaries; s-state radial Laplacian; two independent routes agree SYMBOLICALLY, diff exactly 0).
- Resolved: `step1_native_molecule.py`'s "DISAGREE — investigate" was the GAUSSIAN reference's
  one-electron cusp error; native block right to 60 digits. Investigation cancelled.
- Honest scope: minimal-basis numbers; certified = assembly + precision, NOT accuracy (exact≠accurate).
- Paper 58: new sec:qfd + 3 backing-table rows + sugiura1927 bibitem.

## T2a — certified reference artifact

- `benchmarks/certified_reference/` (generator + JSON, 51 entries) + `docs/certified_reference_values.md`
  + `tests/test_certified_reference_values.py` (10). Categories: T2 (66-digit + pending-120 placeholder
  at zero claimed digits), 26 two-center ERIs at 50 digits (worst quadrature agreement 9.5e-12),
  9 exact-rational Slater integrals, 9 resurgent-data entries, 5 anchors. 17 entries exact.
- **Engine finding:** shell-route hybrid (l>0) requires Z_B < Z_A strictly (NaN at equality =
  removable coincidence; Z_B>Z_A hits an unimplemented Ei branch). Documented in module + paper.
- No DOI action (PI-gated).

## Probe A — spectrum-aware QSVT (debug/probeA_*)

- Constrained interpolation on the known {1±σ_k} with the |p|≤1-on-[−1,1] QSP constraint:
  d_aware ≈ 0.19κ vs honest generic 2.3κ — ratio 8.5–10×, same exponent; **proof** it can only be
  constant (Bernstein Θ(κ) floors both ways; quadratic clustering spaces the bottom eigenvalues at
  exactly the x^{−1/2} variation scale). Water A₁: 10–11×. At P60's N=16: 1050 (model) → 209
  (honest generic) → 27 (aware).
- σ-law alone insufficient (11–21% off at λ_min ⇒ tailored polynomial fails); exact spectra are an
  O(N³) values-only classical SVD — no sparsity cost.
- **The scoped follow-on:** PSD shift → d ≈ 2.8√κ (endpoints from eq:sigma_law in closed form), but
  naive LCU 1-norm ≈3 restores linear-in-κ (worse than baseline); with ~3× amplification net ≈500
  vs ≈1900 at κ=864, growing with κ. Gating question = amplification pricing. Captured in P60
  sec:resource.
- Solver hygiene: LP vs exact Chebyshev minimax 1e-14; independent non-LP confirmation 3/3;
  convention trap (subnormalization head-room not a gauge: c=√λ_min grazes the QSP ceiling,
  convergence degrades to algebraic; all arms use c=½√λ_min).

## Probe B — analytic momentum factorization (debug/probeB_*)

- λ(spec convention): H₂ 1.105× (loses), LiH 0.937×, 4-orb 0.76–0.79× — but the **eigen-Cholesky DF
  control** on the same tensor matches within 1–13% at rank 3–15 vs 25k–43k momentum leaves ⇒ the
  gains are generic DF, not momentum. Proof: λ_spec ≤ λ_absCS ≤ λ_abs and λ_absCS ≥ λ_std for any
  positive-weight squared-one-body factorization ⇒ the stated |·|-convention could never win.
- Framing corrected: NOT THC (ρ̃_pq(k) is full-rank; rank-1-leaf version is the fitted grid one).
- Small-k audit: no singularity (Jacobian cancels 4π/k²); the k→0 leaf → identity ⇒ N̂², suggesting a
  normal-ordering variant (wins 0.54–0.82× but the control side wasn't normal-ordered — deferred, open).
- M(ε) basis-size-independent at fixed molecule (bandwidth/extent-set); classical data 228–9722× the
  unique-ERI list. Traps caught: mpmath.quad [0,∞] silently 38% wrong on the oscillatory integrand
  (fixed by zero-partitioning); the prolate direct-FT reference itself inaccurate above k≈60.
- Captured in P60 sec:resource; §3 row added.

## u1/u2 — the 120-digit T2 push

- u1 (dps 145, K=580, 12 chunks) assembled:
  T2_u1 = 0.3953557659017139643252292968048475642605639778670821089352342654695155086850672974258235932094071533877351035725537132194034178652479486089
  — **confirms all 66 certified digits; digits 67+ single-method until u2**; hi1's unclaimed 67+
  (…6836…) REFUTED (diverged exactly at the known cross-validation boundary — the claim discipline
  validated). u2 (dps 150, K=560, different panel/node/s-map config) launched + monitored; assemble:
  `python debug/beta2_assemble.py u2 150 560 240`, then the h≤1e4 PSLQ per beta2_u1_README.md.

## γ/ln tagging reconciliation (build plan §8.5.4)

The standing "ln and γ UNTAGGED" note is discharged at the resurgent level by v4.104.0: γ is
coordinate bookkeeping (never exists in the Borel plane — the origin-singularity cost of the R
variable), and the exchange/hybrid logs are SKELETON-FORCED boundary data (arguments = rational
monomials in the Borel positions: κ=Πλ_j^{q_j}, the cross-ratio Λ). Paper-34-chain placement:
not observation-projection transcendentals — connection-data boundary terms of Layer-2 objects.
Build-plan note updated.

## Follow-ons

- u2 assembly + 66→~120-digit adjudication + the h≤1e4 PSLQ (machinery ready; do NOT claim u1's
  extension before u2 agrees).
- Isoenergetic atomic algorithm at circuit level (the approved reach #3).
- PSD-shift amplification pricing (Probe A's scoped follow-on).
- Normal-ordered λ comparison (Probe B's deferred variant).
- Hybrid evaluator Z_B≥Z_A branch (coverage gap).
- Optional: DOI-stamp the certified artifact (PI).
