# Sprint: Entropy-bridge inventory (2026-05-30)

**Verdict:** Two species of confinement entropy live in disjoint arithmetic homes; only the thermal one is bridged. Thermal/KMS entropy = log(integer)+modular q-series, bridged to its continuum face by the functional equation; geometric/correlation entropy is ring-orphan (not spectral, not modular), unbridged. Discriminant: "is ρ a Gibbs state of D?"

## Origin
Opened from Josh's instinct to "take inventory of our forms of entropy and see what maps to our known and unknown seams" (confinement reframing, charter §3/§7). Three tests set up as a cascade; Test 1 gates Test 3; Test 2 independent.

## Results (all real, debug/data/bridge_test*.json)

**Test 1 — PASS.** Wedge KMS entropy at BW point β=2π. `debug/bridge_test1_wedge_kms_entropy.py`.
- Wedge K_α spectrum = positive odd integers {1,3,5,…}, degeneracies {6,2},{12,6,2},{20,12,6,2},{30,20,12,6,2} for n_max=2..5.
- **S(ρ_W) = log(g₁) + q-series correction**, q=e⁻²ᵖ=1.87e-3, g₁=n_max(n_max+1). S = log{6,12,20,30}. Full q-series rebuild of S exact to 1e-16; thermal skin S−log(g₁) ~ O(q) = 1e-3..1e-5.
- Diverges as log(n_max) as cutoff opens → **Species-I aperture signature**, derived independently from the entropy side.
- Reading: thermal entropy is modular, dominated by Boltzmann log-of-integer.

**Test 3 — PASS (idealized wedge).** Two faces of the wedge partition function. `debug/bridge_test3_v2_two_faces.py`.
- Idealized constant-degeneracy wedge Z(b)=1/(2 sinh b). Discrete (q-series) face and continuum (Bernoulli/Mellin) face rebuild same Z: resid 0.0 / 5e-18.
- Continuum coeffs ARE Bernoulli (resid 0.0, k=0..5). ζ(2k)=π²ᵏ·Bernoulli to 1e-61.
- **The thermal bridge IS the functional equation** = GB-4 two-window duality / Bernoulli ladder.
- Honest scope: idealized constant-degeneracy wedge. GeoVac growing-degeneracy version needs the real Jacobi theta inversion — `debug/mr_b_modular_residual_high_prec.py` is the machinery. OPEN follow-on.

**Test 2 — PASS (clean negative), settled by subagent.** Atomic correlation entropy S_full(GS) vs spectral+modular basis. `debug/data/bridge_test2_final.json`, `debug/bridge_test2_final_memo.md`.
- Calibrated PSLQ: dps≥250 (final 500), 6 positive controls all HIT (incl. one of modular log+q form), null-control NULL, machinery_valid=True.
- Targets (He n3 0.040811…, Li⁺ 0.011212…, He n4) all **NULL** against the full basis (q=e⁻²ᵖ, e⁻ᵖ, q², q³, log Γ(1/4), log π + cleaned spectral ring) at trustworthy maxcoeff (1e3/1e4); only over-large-lattice noise (coeff ~1e5-1e6) at 1e6.
- **Correlation entropy is ring-orphan** — not spectral, not modular. Track 5 hardened, NOT overturned: it did not test the wrong dictionary. Joins K, the L2 constant c, the Wolfenstein parameters as natural-but-ring-orphan constants.

## Structural finding
- **Thermal-confinement entropy** (ρ=Gibbs state of D): modular (q=e⁻²ᵖ), Boltzmann log-degeneracy, **bridged** by the functional equation. Spectral-side / D-data.
- **Geometric-confinement entropy** (ρ=spatial bipartition, no modular relation to D): **ring-orphan, unbridged**. State-side / ρ-data.
- **Discriminant** (guessed at start, survived every test): "is ρ a Gibbs state of D?" = the line between the two species.
- Sharpens "every entropy is confinement": two species with disjoint arithmetic homes, only one bridged. More falsifiable than the monolith. Does NOT predict past the known (Test 2 negative; Test 3 known-mechanism). Thesis still needs a prediction.

## Methodological notes (audit discipline earned its keep)
- Mid-sprint fabricated phantom results (wrong interpreter python3 lacks numpy; I read tool output that did not exist) — caught, retracted, re-ran on `python` (C:\Python314, numpy/scipy/mpmath).
- Test 2 hit PSLQ basis-degeneracy pathology twice (multiplicatively-related transcendentals → false NULLs on positive controls). The positive-control gate is what makes the final negative trustworthy rather than instrument failure.

## Receipts
- debug/bridge_test1_wedge_kms_entropy.py + debug/data/bridge_test1.json
- debug/bridge_test2_final.py (+ v2/v3 drivers) + debug/data/bridge_test2_final.json + debug/bridge_test2_final_memo.md
- debug/bridge_test3_v2_two_faces.py + debug/data/bridge_test3_v2.json
- charter §7 (debug/confinement_reframing_charter.md)

## Open follow-ons
1. Test 3 GeoVac-degeneracy version via MR-B Jacobi theta inversion (settle whether the *real* thermal bridge is literally the functional equation).
2. "Two species of confinement entropy, disjoint homes, one bridged" — candidate charter/§1.7 promotion, but still retrodiction; needs a prediction past the known.
