# Sprint: three-track parallel follow-on (β(2)-in-T2 / exchange-γ resurgence / Davidson decider) — 2026-08-21

**Origin:** PI direction "all parallel" on the v4.103.0 investigation queue. Three opus agents.
Canonical memo (one per sprint); per-track detail: `debug/beta2_track_a_findings.md`,
`debug/exch_gamma_track_b_findings.md`, `debug/davidson_track_c_findings.md` + drivers
`debug/beta2_t2_*.py`, `debug/exch_gamma_*.py`, `debug/davidson_*.py`.

## Verdicts

| Track | Question | Verdict |
|:--|:--|:--|
| A | T2 to ≥32 digits + the pre-registered β(2) (T-2) ring test | **GO**: 66 cross-validated digits via a NEW exact factorization; **PSLQ DECISIVE-NEGATIVE across every ring incl. the Catalan-completed wt≤3 ring (h≤10 @ 64 digits)** |
| B | Exchange-class {E₁, ln, γ} resurgence (the 5th pattern test) | **GO**: pattern **5/5**; Stokes charges INTEGER; **γ never exists in the Borel plane** (coordinate bookkeeping) |
| C | Matrix-free Davidson CI → banked n_max=4 A-vs-C decider | **GO**: decider = **A (irreducible wall)**, C falsified; + a production **sign bug found and fixed** |

## Track A — T2 = 0.395355765901713964325229296804847564260563977867082108935234265469…

- **The lever (new, paper-captured as eq:kw):** j₀(kb)=∫₀¹cos(kbw)dw decouples the b=s+t phase;
  P(s,k)=P(1−s,k) makes the inner integral real ⇒ **T2=(8/π)∫dk∫dw cos(kw)·R(k,w)²** — the (s,t)
  double integral is the SQUARE of a 1-D integral. The complex off-axis (s,t) singularity (the old
  13–14-digit wall) does not exist in this frame; no fibre is evaluated; the precision question is
  the closed-form Watson k-tail (truncation floor e^{−K/2}, K free). Guardrail: the §3 dead-end
  (brute Φ(ρ) quadrature toward the cusp) was NOT the route taken.
- **Certification (honest accounting, findings §9):** strictest two-complete-pipelines reading = 19
  digits; the 66 digits rest on a decomposed certification — identity vs an independent 2D
  evaluation to 69–96 digits; six parameter-disjoint runs (bit-identical at 50 digits); two parallel
  configurations to 1.7e-67; K-independence bounding the tail below 4.4e-33.
- **The frozen anchor was wrong in digit 19** (0.…9641 → …96432) — exactly the v4.97.0 ceiling
  warning. Paper 59 corrected.
- **PSLQ (64 digits, decoy-calibrated, two precisions):** corrected wt≤3 ring {π, K(½)^±1, G}
  DECISIVE-NEG at h≤10; disc-4 wt≤1/2/3 and disc-8 wt≤1/2 DECISIVE-NEG at h up to 1e12; eight
  targeted Catalan probes NEG at h≤1e12. **No conductor-4 (β(2)/Catalan) content in the physical
  integrated observable at clean heights** ⇒ seam test T-2 answered NEGATIVE; the finite
  wt-3-requiring-G element is confirmed shadow-only. Riders: (i) height budgets scale as 10^(D/n)
  in ring dimension n — the old "decidable at ≥32 digits" was never true for the dim-20 ring
  (pre-registrations must calibrate against ring dimension); (ii) the negative excludes CLEAN
  closed forms, not large-height ones (h>10 in that ring remains logically open).
- **Protocol catch:** the corpus PSLQ driver's false-positive height scale used 10^(D/(n−1));
  correct is the vector length 10^(D/n), with maxcoeff calibrated below it — with the old formula
  every leg returned decoy-matched garbage relations on the first pass; **the decoy guard caught
  it exactly as designed**.
- 120 digits (h≤1e4 on the dim-20 ring) is now a parallel compute question (K≈580, hours on 16 cores).
- **T4c reconciliation:** the v4.103.0 clause "≥32 digits is a resummation problem, not a
  quadrature problem" was true of every chart THEN tried; the (k,w) chart dissolves the digit
  question entirely (the invariant Gevrey structure persists but relocates to a controllable
  truncation floor). Digits and closed form are separate questions; the closed form stays [OPEN].

## Track B — exchange class: pattern 5/5, γ tagged

- **Reduced normal form (exact, 6/6 rate pairs):** F = [e^{−AR}(γ+ln κR) + e^{(a−b)R}E₁(2aR) +
  e^{(b−a)R}E₁(2bR) − e^{AR}E₁(2AR)]/(2abR²), A=a+b, κ=2ab/A — single sector (vs the hybrid's three).
- **Borel transform closed-form:** 𝔅(ξ)=(1/2ab)Σ_j q_j(ξ−ξ_j)ln(ξ−ξ_j), ξ_j∈{0,−2a,−2b,−2A},
  charge-neutral, **integer Stokes charges q=(−1,+1,+1,−1)** — no radical (hybrid needed ℚ(√d)).
  Verified: Laplace 90.4 digits / blind Prony-over-ℚ exact / branch cut 2.8e-60; blind extractor
  validated on synthetics (pole, log, u·ln u, √2 charges) BEFORE use.
- **γ:** never exists in the Borel plane. The R-dependent log moves one Borel singularity onto the
  sector's own origin (pole → u·ln u); γ = the bookkeeping cost of writing that origin singularity
  in the R variable. **Tags the γ left UNTAGGED since v4.77.0** (Increment 3c): coordinate
  bookkeeping, not a projection transcendental.
- **Boundary forcing, tighter than the hybrid:** κ = Π_j λ_j^{q_j} — the exponents ARE the charges.
- **Multivaluedness price:** 𝔅 analytic on (0,∞) ⇒ Borel-summable, no lateral ambiguity (gap
  4.6e-89), but F is genuinely multivalued in complex R (hybrid was single-valued) — costing no
  new transcendental.
- **π-power law refined:** integral local Borel exponent (pole OR log branch) → 2πi; half-integral
  → π⁰ / 1/√π. Decoys: 5/5 κ-decoys + 6 charge-vector decoys fail (2 information-free accidental
  hits at (a,b)=(1/2,7/2) documented — why the sweep uses six pairs).
- Open leg: the heteronuclear infinite τ-sum's resurgence.

## Track C — decider = A; sign bug found and FIXED

- **Bug:** `coupled_composition._double_excitation_phase` returned MINUS the correct same-spin
  double phase (computed for a†_s a†_r a_q a_p, contracted against a†_r a†_s a_q a_p). Confirmed
  vs brute-force Fock-space FCI (flip ⇒ max|ΔH| 1.4e-14). Effect: ~+4 mHa max_n-independent
  spurious over-binding on balanced LiH; geometry <0.2 pp. **FIXED in this close** (global −1 at
  the return; comment states the mechanism; pin test
  `test_library_same_spin_double_phase_is_correct` guards BOTH the corrected sign and the
  historical convention). All 79 consumer tests green post-fix. **Reconciliation:** Paper 19's
  published analytical values (−8.055 Ha, 0.20% at n_max=3) match the CORRECTED physics
  (−8.055224, 0.191%) — the fix brings the shipped code back into agreement with the published
  paper; no paper number changes. (Residual honesty: that the analytical run avoided the buggy
  path is inferred from the exact match, not from provenance.)
- **Engine:** `geovac/balanced_direct_ci.py` — string-based sigma + Davidson;
  validated ≤7e-15 vs brute force, 1.8e-15 Ha vs live library at n_max=2; protocol calibration
  reproduces the banked ABC/chem-error fits to 4–5 sig figs at n_max=2 AND 3. Throughput:
  n_max=3 8,300 s → 5.0 s (~1600×); n_max=4 (16,040,025 dets) 233–315 s/solve, 2.2 GB.
  The wall moved to the O(M⁴) integral build (~2,700–3,300 s/point, cached).
- **Decider (registered in the ABC memo): A.** curv(R_true)/k_true = 2.245×/2.219×/2.210× at
  n_max=2/3/4 → extrapolated 2.204× (99.5% of gap survives; robust across geometric/power-law/1/n
  models); tilt worsens to −0.0399; R_eq → +10.1% worsening; ω_e at own min 2040→1946→1901 →
  ~1861 (+32%) vs true 1406. Everything converges geometrically at a common ratio (~0.37) to
  WRONG geometry ⇒ irreducible orbital-basis-limit wall, NOT slow healing. Closing the curvature
  gap needs ~126 shells at undecaying step. Side-correction: "ω_e stays +45%" was an
  n_max=2-only fit (measured decline 45.1→38.4→35.3%).
- Caveats: 3-point extrapolation (robust, not a proof); single system (LiH); solver asserts
  N_α=N_β=2 (BeH₂ needs a general-N sigma). Cache `debug/data/ints/` (~3 GB) is deletable.

## Captures applied (papers compile clean; certified → compounds Phase-4 re-review OWED)

- **Paper 59:** eq:kw + 66-digit value + decisive PSLQ verdicts + both calibration riders; anchor
  marked correct-to-18; two frame-scoped corrections ("definitive evaluation awaited a different
  representation"; "brute quadrature not a route *in the (s,t)/co-area frames*"); the
  bessel_algebra PSLQ paragraph updated (decision made at 64 digits); resurgent-skeleton
  [OBSERVATION] upgraded to 5/5 with the exchange normal form, integer charges, γ-tagging,
  multivaluedness price, refined π-power law.
- **Paper 56 rem:paper59_cm:** T-2 run and negative; conductor-4 confined to family + shadow;
  "≥32 digits" rider corrected; seam = convergent on all tested legs.
- **Paper 19:** n_max=4 decider paragraph (irreducible wall) + reproducibility note (phase fix;
  no paper number changes); conjecture-status line sharpened. NIT (pre-existing, not fixed):
  3 dangling \ref's (sec:relationship, eq:Vpk, sec:conv_w1e_cross_domain_wall) pre-date this edit.
- **Tests (all green):** `tests/test_paper59_t2_value.py` (3: identity, tail law, digit-19),
  exchange legs added to `tests/test_paper59_resurgent_skeleton.py` (4 total),
  `tests/test_balanced_direct_ci.py` updated to the corrected convention (11).
- **Code:** `geovac/coupled_composition.py` sign fix; `geovac/balanced_direct_ci.py` (new tracked).

## Follow-ons (not started)

- 120-digit T2 (h≤1e4 on the pre-registered ring) — parallel compute, machinery built.
- Heteronuclear τ-sum resurgence (Track B's open leg).
- General-N sigma for `balanced_direct_ci` (BeH₂/H₂O deciders).
- Optional: fix Paper 19's three pre-existing dangling refs.
- The closed form of T2 (unchanged): second-cusp Borel–Lambert resummation, now decoupled from
  the digits question.
