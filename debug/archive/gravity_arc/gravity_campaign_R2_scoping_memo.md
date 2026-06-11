# Gravity Campaign — R2 scoping: propinquity convergence of the cigar tip term

**Date:** 2026-05-29
**Type:** Scoping (main-session, no compute). Maps the Paper 38/39/40/45/47 five-lemma propinquity machinery onto the cigar tip-term geometry. Identifies transport-verbatim vs new-work, and the load-bearing difficulty.
**Deliverable question (R2):** Does the discrete cigar tip-term spectral triple converge to the continuum horizon triple in the Latrémolière propinquity, with a rate — so that S_tip → A/4?

## 1. The geometry (from G4-3 + G4-5d)

Cigar near-horizon = warped product, factorizing at constant warp r(ρ)=r_h:
$$\mathcal{G}_{\rm cigar} = \mathbb{Z}_+(a)|_{N_\rho} \times \mathbb{Z}/N_\phi \times \mathrm{Fock}(S^2, l_{\max}), \qquad K_{\rm cigar}(t) = K_{D^2}(t)\cdot K_{S^2_{r_h}}(t).$$
- **Disk D²** = (radial interval [0,R], measure ρdρ) × (azimuthal S¹). The **conical defect** α = N_φ/N_0 lives in the azimuthal sector (apex angle 2πα; α=1 smooth, α≠1 cone with singularity at ρ=0).
- **S²** = the spatial sphere; supplies the **area** A = 4π r_h².
- **Entropy** = φ(0) Mellin moment of the replica difference tip(t) = (dK/dα)|_{α=1} − K_disk(t) (Task 1 / G4-5d). The bulk A-coefficient cancels in the difference; the entropy is the cone contribution × area.

## 2. Factor-by-factor transport sources

| Factor | What it is | Propinquity transport source |
|:---|:---|:---|
| Azimuthal S¹ (smooth, α=1) | U(1) = compact Lie group | **Paper 38/39/40 verbatim** (U(1) is the abelian case; Leimbach–vS torus). |
| Azimuthal S¹ at α≠1 | S¹ of circumference 2πα | **Paper 47** (varying-period T-circle) — near transport, needs α-parameterization. |
| S² | SU(2)/U(1) homogeneous space (NOT a group) | **Berezin–Toeplitz / Hawkins** — the same machinery Paper 38 L4 already imports. NOT Paper 40 (S² isn't a group). |
| Radial interval [0,R], ρdρ | bounded interval, singular endpoint at ρ=0 (cone apex) | **Paper 47 norm-resolvent** on bounded interval — nearest source, but norm-resolvent ≠ full propinquity. **New-work-leaning.** |
| Tensor assembly | warped product of the above | **Paper 39 + Paper 45 PURE_TENSOR + k-fold master theorem** (joint rate = max of factor rates; C₃ k-independent). |

## 3. Lemma-by-lemma map

- **L1′ (operator-system substrate, propagation number).** Smooth factors transport (Paper 38 S¹/SU(2); Hawkins S²; Paper 39/44 tensor operator system). NEW: the radial-interval-with-apex operator system (singular endpoint at ρ=0). *Transport-mostly; new piece = radial apex.*
- **L2 (central-Fejér rate γ→0).** Joint rate = max of factor rates (Paper 39/45). U(1) and S² rates transport; radial-interval rate needed; conical α modifies the S¹ rate (Paper 47 varying-period). *Transport-mostly.* **Romantic hook:** does the joint rate carry the M1 4/π = Vol(S²)/π² signature? If yes → cigar entropy convergence sits on the Bernoulli ladder, tying gravity to the master Mellin engine.
- **L3 (Lipschitz bound C₃).** Paper 39/45: C₃ is k-independent, = 1 asymptotically. NEW: the **warp couples radial↔S²**, so the warped-product Lipschitz cross-term needs the Paper 45-style vanishing-cross-term check (it vanished there for time-chirality; here it's the warp Jacobian 2r′/r). *Transport-conditional-on-warp-cross-term.*
- **L4 (Berezin reconstruction).** Tensor-product Berezin (Paper 45 PURE_TENSOR). S² Berezin = Hawkins (already in Paper 38). NEW: radial-interval Berezin + the conical apex. *Transport-mostly.*
- **L5 (Latrémolière assembly).** Joint tunneling pair (Paper 45), conditional on L1′–L4. *Transports if the above close.*

## 4. The load-bearing new piece — NOT one of L1′–L5

The entropy is **not a fixed-triple convergence**. It is the **replica derivative** (dK/dα)|_{α=1} — on the discrete substrate a finite-difference in α (G4-5d: [K(α₊)−K(α₋)]/(α₊−α₋), ε=0.1). So R2's real theorem is one level up from Paper 38:

> The family of cone triples {𝒯_α} must converge in propinquity **uniformly in α on a neighborhood of α=1, and differentiably in α**, so that the α-derivative at α=1 (the replica = entropy) converges to A/4 with a rate.

The **fixed-α** convergence of each cone triple transports (factor-by-factor, §3). The **α-differentiability of the propinquity** — convergence of dΛ/dα, not just Λ — is the genuinely new mathematics. No published propinquity result has an α-family-differentiable statement (Mondino–Sämann has families but not replica derivatives). This is the "math ahead of the literature" piece, and it's the single load-bearing requirement.

## 5. R2 execution recommendation (for the convergence proof, when authorized)

1. **Smooth-factor convergence first** (transports): assemble the constant-warp tensor triple S¹ ⊗ S² ⊗ (radial), prove fixed-α propinquity convergence via Paper 39/45 PURE_TENSOR + Hawkins S² + Paper 47 interval. Surfaces the warp-cross-term (L3) and radial-apex (L1′/L4) as the only new sub-pieces.
2. **α-family uniform+differentiable convergence** (the new theorem): prove Λ(α) is C¹ in α near α=1 with a uniform rate, so dΛ/dα|_{α=1} converges. This is the load-bearing deliverable.
3. **Rate signature check** (romantic): is the joint rate's leading constant the M1 4/π? Connects gravity entropy convergence to the Bernoulli ladder / master Mellin engine.

## 6. Honest scope of this scoping

- The transport claims are first-pass structural reads, not proofs. The warp-cross-term (L3) and the radial-apex operator system (L1′) are the two transport-claims most likely to need genuine work, not just citation.
- S² is treated as Berezin–Toeplitz-Hawkins-transportable; G4-1 showed S² Dirac is NOT two-term exact (full SD series), but two-term exactness is not required for propinquity convergence — only the lemma structure is.
- The α-differentiable propinquity (§4) is named as new; it has no verbatim transport source. Estimating its difficulty is itself part of R2 phase-1.

## 8. Phase-1 setup + radial-apex check (2026-05-29 continuation)

Writing the constant-warp tip triple out explicitly sharpens R2 three ways:

**(a) Warp cross-term is zero at constant warp.** The entropy is the near-horizon (constant-warp r=r_h) object, where Δ_cigar = Δ_D² + (1/r_h²)Δ_S² is a clean tensor product (r′=0 ⇒ no Jacobian cross-term). The L3 warp-cross-term flag of §3 only applies to variable warp (G4-3b), which the entropy does not use. Paper 45 PURE_TENSOR applies directly. **L3 warp flag cleared.**

**(b) α lives entirely in the disk; S² is a passive area multiplier.** Disk Dirac in azimuthal modes has spinor momentum m_eff = (k+½)/α. The S² factor (spatial horizon sphere) is α-independent and supplies A = 4πr_h². So S_tip = (disk conical contribution, α-dependent) × A. The replica derivative d/dα acts only on the disk. **R2 reduces from the warped-product cigar to the 2-disk-with-cone D²_α ⊗ passive S².**

**(c) Radial apex is centrifugally screened — radial-apex L1' flag cleared.** Spinor radial modes ~ ρ^{|m_eff|} with m_eff=(k+½)/α ≠ 0, so every mode vanishes at ρ=0. The cone vertex is screened; no boundary condition at the apex; regular Sturm–Liouville on (0,R] with only the IR boundary at ρ=R (Paper 47). This is the SAME no-zero-mode fact that makes spinor SC extraction clean (G4-4; dead-ends table "anti-periodic + half-integer structurally essential") — the property that enables the entropy extraction is the property that regularizes the apex. *Regularity is solid; explicit prop=2 confirmation on the assembled disk operator system is a clean follow-on compute, not run here.*

**Net re-scope.** All scaffolding transports: azimuthal S¹ (Connes–vS Toeplitz / Paper 38), radial interval (Paper 47, apex screened), S² area (Hawkins, passive), tensor assembly (Paper 45 PURE_TENSOR), warp cross-term (zero), radial apex (screened). **The single irreducible new piece is α-differentiable propinquity** — convergence of dΛ/dα|_{α=1}, the replica derivative, for the 2-disk-with-cone. No transport source (Mondino–Sämann has families, not replica derivatives). This is R2's whole remaining content and will not compress.

## 7. Cross-references
- `debug/g4_3_warped_substrate_memo.md` — cigar substrate, warped-product factorization
- `debug/g4_5d_cutoff_dependence_memo.md` — entropy = φ(0) moment of tip = (dK/dα − K_disk)
- `debug/gravity_campaign_phase1b_tip_bulk_independence_memo.md` — Task 1 SEPARABLE verdict
- Paper 38 (single-factor SU(2)), Paper 39 (tensor product), Paper 40 (universal 4/π, compact Lie groups), Paper 45 (PURE_TENSOR + vanishing cross-term), Paper 47 (bounded interval / varying period, norm-resolvent)
