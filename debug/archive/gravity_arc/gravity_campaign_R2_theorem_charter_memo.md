# Gravity Campaign — R2 theorem charter (α-differentiable spectral-functional convergence)

**Date:** 2026-05-29
**Type:** Theorem charter (main-session, no compute). Precise statement + proof skeleton for the single irreducible R2 deliverable. Follows the R2 scoping memo + phase-1 check.

## Framing correction (sharpened by writing the charter)

The floor was previously phrased "α-differentiable propinquity." Precise object: **the heat-trace spectral functional is what is differentiated in α, NOT the propinquity metric.** The entropy is the replica derivative of the heat trace. The propinquity is the *backbone* (certifies geometric faithfulness of the discrete triple at each α), undifferentiated. This is more correct and more tractable than differentiating a metric-on-triples.

## Objects

- 𝒯_{n,α}: discrete constant-warp cigar tip triple at truncation n, apex angle 2πα = (2-disk-with-cone D²_α) ⊗ (passive S²(r_h)). α-dependence enters only through spinor azimuthal momentum m_eff = (k+½)/α.
- 𝒯_α: continuum target (continuum disk-with-cone ⊗ S²).
- K_n(α,t) = Tr e^{-t D_{n,α}²}: discrete heat trace. Tip(α,t) = (dK/dα)|_{α=1} − K_disk (bulk cancels, Task 1).
- S_tip^{(n)} = M_{φ(0)}[ Tip_n(t) ] · (area normalization), where M_{φ(0)} is the zeroth Mellin moment (G4-5d).

## Theorem (target)

S_tip^{(n)} → A/4 as n→∞, with a rate, where A = 4πr_h² is the S² horizon area.

## Proof architecture — three layers (Paper-47-style two-rate hybrid + new α-layer)

**Layer 1 — Backbone: propinquity (TRANSPORTED, Paper 45 PURE_TENSOR).**
Λ(𝒯_{n,α}, 𝒯_α) ≤ C₃·γ_n → 0 at each fixed α. Certifies the discrete triple is a faithful NCG approximation. Assembled: azimuthal S¹ (Connes–vS Toeplitz / Paper 38) ⊗ radial interval (Paper 47, apex centrifugally screened) ⊗ S² (Hawkins Berezin, passive α-independent area). Warp cross-term zero at constant warp. C₃ k-independent (Paper 39/40 master theorem).

**Layer 2 — Spectral functional: norm-resolvent (ADAPTED, Paper 47).**
D_{n,α} → D_α in norm-resolvent sense ⇒ K_n(α,t) → K(α,t) at each fixed (α,t), with a rate. This is exactly Paper 47's two-rate-hybrid architecture: propinquity backbone (Layer 1) + norm-resolvent functional (Layer 2). The radial-interval norm-resolvent machinery is built; apex regular (screened) so resolvent is well-defined at ρ=0.

**Layer 3 — α-layer: NEW (L6, the prize).**
The Layer-2 convergence is uniform in α on [1−δ, 1+δ] and C¹ in α, so lim_n and d/dα|_{α=1} commute. The α-derivative is
dK_n/dα = Σ_modes (-t)(∂λ/∂m_eff)·(-(k+½)/α²)·e^{-tλ}.
The replica derivative pulls down a (k+½) weight — a heavier mode-moment than the heat trace. L6 = uniform-in-α dominated-convergence control of the (k+½)-weighted tail, with a rate. No transport source.

## Rate

O(max of the three layer rates); slowest gates. Layers 1,2 transport their rates; L6's rate is new.

## L6 difficulty assessment (honest)

Tractable, not a wall:
- For t>0, e^{-tλ} decays faster than any polynomial in k ⇒ dominates the (k+½) weight ⇒ absolute + uniform convergence for α bounded away from 0 (fine near α=1).
- The tip's φ(0) moment is IR-weighted (peaks t≈10, G4-5d table) — strongest kernel suppression.
- The small-t UV (where (k+½) could bite) is where the bulk A-coefficient lives, and that CANCELS in the tip (Task 1). L6's natural support is the benign regime.
- Genuine work: the **rate** and the **uniformity**, not bare convergence.

## Romantic hook (now concrete)

M1 signature 4/π = Vol(S²)/π² may enter through the **area itself**: A = 4πr_h² = Vol(S²)·r_h². The passive S² factor literally carries Vol(S²). So S_tip = A/4 may wear the M1 / Bernoulli-ladder signature by construction (through the area), not via a rate coincidence. Check once the rate is in hand → would tie cigar entropy to the master Mellin engine (gravity ↔ α).

## R2 execution sequence (when proof is authorized)

1. Assemble Layer 1 (fixed-α propinquity) explicitly for D²_α ⊗ S² via Paper 45 PURE_TENSOR + Paper 47 interval + Hawkins S². Confirm prop=2 on the assembled disk operator system (the phase-1-check follow-on compute).
2. Establish Layer 2 (fixed-α heat-trace convergence) via Paper 47 norm-resolvent on the radial factor.
3. Prove L6 (Layer 3): uniform-in-α C¹ convergence, the (k+½)-weighted dominated-convergence bound with rate. THE deliverable.
4. Rate-signature check (M1 / area).

## Honest scope (sprint-close)

- **Theorem grade closed this sprint: nothing.** This is the scoping/charter layer. The theorem is *stated*, not proven.
- **Structural sketch (solid analytical, not proof):** the three-layer architecture; the radial-apex centrifugal-screening argument; the Task 1 SEPARABLE verdict (bulk cancels in the tip). These rest on established facts (G4-6b data, G4-5d Mellin map, the G4-4 no-zero-mode fact) but are not formal proofs.
- **Numerical observation:** the Phase 1 Möbius recompute values (dispatched-agent compute).
- **Named open follow-ons:** (1) L6 proof — uniform-in-α C¹ convergence of the (k+½)-weighted replica derivative, with rate (THE prize, unproven, no transport source); (2) explicit Layer 1–2 assembly for D²_α ⊗ S² (transported pieces, not yet assembled); (3) prop=2 confirmation compute on the assembled disk operator system; (4) M1/area rate-signature check.
- **Not yet in any paper.** Findings remain in debug memos until L6 is proven; premature to seed Paper 51 §12.x with an R2 program section.

## Cross-references
- `debug/gravity_campaign_R2_scoping_memo.md` — lemma-transport map + phase-1 check
- `debug/gravity_campaign_phase1b_tip_bulk_independence_memo.md` — Task 1 SEPARABLE (bulk cancels in tip)
- Paper 45 (PURE_TENSOR propinquity), Paper 47 (two-rate hybrid: propinquity + norm-resolvent), Paper 38/40 (C₃, central-Fejér rate, 4/π = M1), G4-5d (φ(0) moment), G4-4 (anti-periodic/half-integer ⇒ no zero mode ⇒ apex screened)
