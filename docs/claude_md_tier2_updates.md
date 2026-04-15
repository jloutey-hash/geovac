# Proposed CLAUDE.md Updates — Dirac-on-S³ Tier 2 Sprint Closure

**Author:** Track T6 worker agent
**Date:** 2026-04-15
**Status:** PROPOSAL — NOT auto-applied. PM to apply after PI approval.
**Precedent:** Follows the format and edit-classification pattern of
`docs/claude_md_proposed_updates.md` (Tier 1 closure, 2026-04-15).

The edits below are ordered by CLAUDE.md section. Each is stated as a
specific "replace" or "append" instruction. Tier 1's edits are assumed
already applied (version v2.10.0 landed); this proposal bumps to v2.11.0
and appends Tier 2's Phase 4J / Section 10 / Section 11 / Section 12
rows alongside Tier 1's existing entries.

---

## §1. Project Identity — version bump

**Edit 1.1** (PM may apply per §13.5 — version number only).

In §1, replace:
```
**Version:** v2.10.0 (April 15, 2026)
```
with:
```
**Version:** v2.11.0 (April 15, 2026)
```

**Rationale:** Tier 2 is a completed multi-track sprint (T0–T6) that
delivers new capability (spinor composed qubit Hamiltonians), a new
Paper 22 result (d_spinor(l_max) table), and a Paper 18 taxonomic
extension (spinor-intrinsic subtier). Per CLAUDE.md §9, a minor
version is the correct granularity for completed sprints that add
features or paper additions.

---

## §2. Current Development Frontier — Tier 2 summary addition

**Edit 2.1** (PM may apply — Section 2 is PM-allowed).

Append a new bullet after the Phase 4I summary (placed at the end of
the α combination-rule backlog block that Tier 1's edit 2.1 added):

```
Dirac-on-S³ Tier 2 sprint complete (April 2026, seven-track sprint T0/T1/T2/T3/T4/T5/T6), ALL POSITIVE. Engineering upgrade: spin-ful composed qubit Hamiltonians for LiH, BeH, and CaH (first instance of a relativistic molecule in the composed pipeline), built end-to-end algebraically in the (κ, m_j) native Dirac-label basis. (T0) d_spinor(l_max) table in both pair-diagonal-m and full-Gaunt conventions, l_max=0..5: ratio d_spinor/d_scalar = 1/4 (l_max=0, pure spin-dilution) → 11/20 → 553/724 → 101/118 → 2533/2820 → 9611/10396 ≈ 0.92 (l_max=5). Paper 22 theorem extends verbatim to spinor basis with sparsity-exponent preserved and prefactor inflated. (T1) `geovac/dirac_matrix_elements.py`: DiracLabel dataclass, κ ↔ (l, σ) bridge, all Szmytkowski angular matrix elements in closed form, all diagonal radial ⟨r^k⟩ and ⟨1/r^k⟩ as Bethe-Salpeter rationals; off-diagonal radial via direct sympy integration fallback. 117 tests pass; 51/51 D1 regression tests preserved. (T2) `geovac/spin_orbit.py`: Breit-Pauli H_SO = −Z⁴α²(κ+1)/[4n³l(l+½)(l+1)], closed form in (n, κ) with exact Kramers cancellation at l=0 (κ=−1 forces numerator zero before the formally divergent ⟨1/r³⟩ is evaluated). Z⁴ scaling verified symbolically across Z ∈ {1, 3, 4, 38}; 2p doublet splitting α²/32 reproduces. 22 tests pass. (T3) `geovac/composed_qubit_relativistic.py`: dispatches from `build_composed_hamiltonian(spec)` when spec.relativistic=True. Three specs `lih_spec_relativistic`, `beh_spec_relativistic`, `cah_spec_relativistic` (BeH and CaH are one-bond reductions of BeH₂/CaH₂). Pauli ratio rel/scalar isostructural 1.00×/2.42×/5.89× at n_max=1/2/3 across all three molecules. λ_ni flat or decreasing vs scalar (QPE-favorable). 13 T3 regression tests + 164 pre-existing tests pass; scalar LiH/BeH₂/H₂O counts 334/556/778 bit-exact. (T4) Sunaga 2025 (PRA 111, 022817) head-to-head: only the RaH-18q cell (47,099 Pauli) is published in the main paper; per-molecule BeH/MgH/CaH/SrH/BaH at 18q are in SI Tables S1–S3, flagged DEFERRED. GeoVac native-Q ratios vs Sunaga RaH-18q: LiH 0.017×, BeH 0.017×, CaH 0.011× (at GeoVac Q ∈ {20, 30}). Projected matched-Q=18 advantage 150×–250× via Paper 14 §IV.B O(Q^2.5) scaling × rel/scalar 2.4×. Fine-structure sanity: He 2³P span −66%, Li 2²P doublet +211%, Be 2s2p ³P span −78% relative error — sign and OoM correct across all three atoms, 20–50% accuracy target NOT met (known limitation of leading-order Zα² + Slater Z_eff; Tier 3+ to lift via Darwin + mass-velocity + multi-electron SS/SOO). (T5) `geovac/spinor_certificate.py`: `verify_spinor_pi_free` walks the sympy expression tree of each T3 coefficient and enforces membership in the Paper 18 spinor-intrinsic ring R_sp := ℚ(α²)[γ]/(γ²+(Zα)²−1). T3's symbolic H_SO block passes cleanly for n_max ≤ 4 and across all κ branches; six negative-control contamination tests (bare π, π², ζ(3), log(Z), E₁, unregistered symbol) reject cleanly. 25 tests pass. `docs/paper18_spinor_subtier_proposal.tex` drafted (drop-in, awaiting PI approval). (T6) Drop-in proposals drafted (NOT auto-applied): `docs/tier2_verdict.md`, `docs/paper14_section5_proposal.tex`, `docs/paper22_spinor_section_proposal.tex`, `docs/paper20_tier2_table_proposal.tex`, `docs/claude_md_tier2_updates.md` (this file). **STRUCTURAL CONCLUSION:** The Paper 18 operator-order × bundle grid now has three of four cells populated (2nd-order/scalar = calibration π, 1st-order/scalar = Tier-1 odd-zeta, 1st-order/spinor = Tier-2 spinor-intrinsic α²+γ); 2nd-order/spinor is flagged as a conjectural slot for future QED-on-S³ work. Tier 3 roadmap: [Kr] frozen-core for SrH/RaH, Martínez-y-Romero γ radial corrections, Darwin + mass-velocity, multi-electron SS/SOO. Data: `debug/dirac_t0_memo.md`, `docs/dirac_matrix_elements_design_memo.md`, `docs/spin_orbit_design_memo.md`, `docs/spin_ful_composed_design_memo.md`, `docs/tier2_market_test.md`, `docs/tier2_t5_verdict.md`.
```

**Rationale:** Tier 2 is purely additive infrastructure + taxonomy. It
extends the framework into relativistic chemistry without closing or
re-opening any research program. The summary mirrors the style of the
Phase 4I entry: named tracks, concrete numbers, explicit honest
caveats (Sunaga SI deferred, fine-structure accuracy target not met),
file pointers.

---

## §3. Approaches That Failed — no new rows

Tier 2 produced no failed approaches. All seven tracks met their
success criteria. **No §3 update is required for Tier 2.**

**Flag for PI review:** should the existing "Hopf graph morphism"
entry (§3, from Phase 4B Track α-D) be updated to note Tier 1's
Hopf-Dirac D4 track is also negative for common-generator unification?
Recommendation: **no**. D4's negative is already captured in Tier 1's
§3 edit 3.1 ("Dirac-sector lift of Paper 2 α combination rule
ingredients B, F, Δ" row, count = 3). Adding a reference from the α-D
row would be redundant and muddle the distinction between the scalar
Hopf quotient (α-D) and the Hopf-equivariant Dirac decomposition (D4).

---

## §10. Validation Benchmarks — Tier 2 rows

**Edit 10.1** (PM may apply — Section 10 is PM-allowed for appending).

In §10's benchmark table, after the Tier 1 Dirac-on-S³ rows (added by
Tier 1 edit 10.1) and before the final "Speed regression" row, append:

| Test | Max Error | Purpose |
|:-----|:---------:|:--------|
| T0 d_spinor pair-diag at l_max=0..5 | exact rational | 1/4, 11/20, 553/724, 101/118, 2533/2820, 9611/10396 from sympy Wigner 3j |
| T0 d_spinor full-Gaunt at l_max=0..5 | exact rational | 1/4, ..., 0.923 (sec'dry table, physically correct Coulomb rule) |
| T0 d_spinor ≤ d_scalar | monotonic ∀ l_max | Spinor basis sparsity bounded by scalar sparsity |
| T0 scalar density reproduces Paper 22 | bit-exact | Pair-diag convention matches Paper 22 Table §III |
| `dirac_matrix_elements` module tests | 117 tests pass | Angular (Szmytkowski) + diagonal radial (Bethe-Salpeter) + κ↔(l,σ) bridge |
| σ·r̂ reduction identity | `(−κ, m_j, −1)` exact | Szmytkowski Eq. 2.7 verification |
| Diagonal ⟨1/r³⟩_{n,l} hydrogenic | Z³/[n³·l(l+½)(l+1)] | T1 closed form, verified to sympy machine precision |
| T1 ⟨1s\|r\|2s⟩ at Z=1 | −32√2/81 | Off-diagonal radial sympy integration sanity |
| `spin_orbit` module tests | 22 tests pass | H_SO = −Z⁴α²(κ+1)/[4n³l(l+½)(l+1)] closed form |
| SO Kramers cancellation at l=0 | H_SO = 0 exact | κ=−1 forces numerator zero before ⟨1/r³⟩ evaluated |
| SO Z⁴ scaling | (Z/Z_ref)⁴ symbolic | Verified at Z ∈ {1, 3, 4, 38} |
| 2p doublet splitting (Z=1) | α²/32 exact | Breit-Pauli fine-structure benchmark |
| `spin_ful_composed` module tests | 13 tests pass | Full composed-rel pipeline regression |
| LiH/BeH₂/H₂O scalar regression | 334/556/778 Pauli preserved | Bit-exact scalar path unchanged when relativistic=False |
| LiH relativistic Pauli at n_max=1 | exactly 9 | Matches scalar (no spin-orbit at l=0) |
| LiH rel/scalar Pauli ratio at n_max=2 | ∈ [1.9, 3.0] | Pinned at 2.42× in regression suite |
| Spinor composed block-diagonal ERI | zero cross-block entries | Factorization preserved in relativistic path |
| α → 0 zeroes H_SO diagonal | exact zero | Non-relativistic limit verification |
| `spinor_certificate` module tests | 25 tests pass | Ring R_sp = ℚ(α²)[γ]/(γ²+(Zα)²−1) enforcement |
| Contamination rejection (π, π², ζ(3), log, E₁, unregistered) | raises SpinorTaxonomyError | Six negative controls |
| T3 H_SO block R_sp membership | passes n_max ≤ 4 | Every κ-branch coefficient in ring |
| Sunaga RaH-18q baseline | 47,099 Pauli (published) | Single calibrated cell from PRA 111, 022817 |
| GeoVac rel/Sunaga RaH-18q ratio | 0.011×–0.017× native-Q | Resource advantage verified |
| Fine-structure Li 2²P splitting | sign + OoM correct | Breit-Pauli + Z_eff sanity |
| Fine-structure He 2³P span | sign + OoM correct | 66% relative error (accepted) |
| Fine-structure Be 2s2p ³P span | sign + OoM correct | 78% relative error (accepted) |

**Rationale:** Each row corresponds to a concrete test or a pinned
symbolic identity from Tier 2. Total new test count: 177 (117 + 22 + 13
+ 25). All tests pass in the T0–T5 memos.

---

## §11. Topic-to-Paper Lookup — Tier 2 rows

**Edit 11.1** (PM may apply — Section 11 is PM-allowed for appending).

In §11's topic-to-paper table, after the Tier 1 Dirac rows (added by
Tier 1 edit 11.1), append:

| Topic | Paper | Section | Tier |
|:------|:-----:|:--------|:----:|
| Spinor composed encoding (Tier 2) | 14 | §V (new) | Core |
| Dirac matrix elements in (κ, m_j) basis | 14 | §V | Core |
| Szmytkowski spinor matrix elements | 14 | §V | Core |
| Martínez-y-Romero radial recursion (reserved) | 14, 18 | §V, §IV | Core |
| jj-coupled Coulomb ERI (Dyall §9 separable form) | 14, 22 | §V, §III | Core |
| Breit-Pauli spin-orbit in (κ, m_j) | 14, 18 | §V, §IV | Core |
| Kramers cancellation at l=0 | 14 | §V | Core |
| d_spinor(l_max) sparsity density | 22 | spinor section (new) | Core |
| Spinor-block angular sparsity theorem | 22 | spinor section | Core |
| Pair-diagonal vs full-Gaunt convention | 22 | spinor section | Core |
| Universal partition (sharpened to spinor) | 22 | spinor section | Core |
| Sunaga 2025 head-to-head | 14, 20 | §V, Tier-2 table (new) | Core/Applications |
| OpenFermion-Dirac baseline (deferred SI cells) | 20 | Tier-2 table | Applications |
| Relativistic resource benchmarks (LiH/BeH/CaH) | 20 | Tier-2 table | Applications |
| Fine-structure sanity (He/Li/Be) | 14, 20 | §V, Tier-2 table | Core/Applications |
| α² and γ taxonomy (spinor-intrinsic) | 18 | §IV subtier (new) | Core |
| Ring R_sp = ℚ(α²)[γ]/(γ²+(Zα)²−1) | 18 | §IV | Core |
| Operator-order × bundle grid (4 cells) | 18 | §IV | Core |
| π-free spinor certificate | 18 | §IV | Core |
| First-order operator on spinor bundle | 18 | §IV | Core |

---

## §12. Algebraic Registry — Tier 2 rows

**Edit 12.1** (PM may apply — Section 12 is PM-allowed for appending).

Add a new subsection to §12 titled **"Level 5 — Spin-ful composed
(Tier 2)"** at the end of the existing Level 5 block, with rows:

| Matrix Element | Status | Notes |
|:---------------|:------:|:------|
| Szmytkowski angular σ·r̂, J², L·S, L², σ² | algebraic | Exact eigenvalues in (κ, m_j); integer / half-integer Kronecker deltas |
| Diagonal hydrogenic ⟨r^k⟩, ⟨1/r^k⟩ | algebraic | Bethe-Salpeter rationals; ⟨1/r⟩=Z/n², ⟨1/r³⟩=Z³/[n³l(l+½)(l+1)] (Z³ diverges l=0, suppressed by Kramers) |
| Off-diagonal ⟨n'l'\|r^k\|n l⟩ | algebraic | Direct sympy integration of assoc_laguerre; ~0.1–1s per call |
| Breit-Pauli H_SO in (κ, m_j) | algebraic | Closed form H_SO = −Z⁴α²(κ+1)/[4n³l(l+½)(l+1)], diagonal; exact Kramers at l=0 |
| jj-coupled angular coefficient X_k(κ_a, m_a, κ_c, m_c) | algebraic | sympy wigner_3j, cached, full-Gaunt selection |
| Spinor composed two-body ERI | algebraic | X_k·X_k·R^k factorization, radial via hypergeometric_slater (exact Fraction or machine-float) |
| Relativistic γ = √(1−(Zα)²) | algebraic | Degree-2 algebraic over ℚ(α²) via γ² + (Zα)² = 1; reserved symbol in T1/T2, not bound at Tier 2 |
| Dirac-Coulomb radial (Martínez-y-Romero) | algebraic-pending | Three-term recurrence identified; not wired in T2; Tier 3 deliverable |
| Darwin + mass-velocity α⁴ corrections | algebraic-pending | Closed-form in (n, κ); out of scope at Tier 2; Tier 3+ |

**Rationale:** Every Tier 2 deliverable is algebraic — no numerical-required
entries — making the Level 5 / spin-ful block cleanly separable in the
registry. Two pending rows (γ, Darwin + mass-velocity) mark Tier 3
entries that have identified algebraic routes but are not yet wired.

---

## Summary of edit classification

| Edit | CLAUDE.md section | PM access per §13.5 | Requires PI review? |
|:----:|:------------------|:-------------------|:-------------------|
| 1.1 | §1 version number | Allowed (version ONLY) | No |
| 2.1 | §2 append Tier 2 summary | Allowed (PM-allowed) | **YES** — PI should confirm the "Tier 2 Dirac-on-S³ sprint complete" framing and the Tier 3 roadmap language before commit |
| —   | §3 no new rows | — | No (PM flag noted above re: α-D / D4) |
| 10.1 | §10 append benchmark rows | Allowed | No |
| 11.1 | §11 append topic rows | Allowed | No |
| 12.1 | §12 append Level 5 spin-ful subsection | Allowed | No |

**Recommendation:** Edits 1.1, 10.1, 11.1, 12.1 can be applied by the
PM without further PI review (all mechanical, within §13.5 scope).
Edit 2.1 is PM-allowed in principle but contains a non-trivial Tier 3
roadmap recommendation; PI should review the summary language before
it lands, especially the "Tier 3 roadmap" bullet and the deferred-SI
framing.

Paper drop-ins (14 §V, 20 Tier-2 table, 22 spinor section, 18 spinor
subtier) are separate PI-review items in the T6 deliverables bundle.
If any of the paper drop-ins is rejected or substantially revised, the
relevant topic rows in Edit 11.1 should be adjusted to match.

---

## Notes for PM

- Do NOT apply any edit until the PI has reviewed
  `docs/tier2_verdict.md` and the four paper drop-ins.
- Apply in order 1.1 → 2.1 → 10.1 → 11.1 → 12.1.
- The Paper 18 §IV subtier drop-in
  (`docs/paper18_spinor_subtier_proposal.tex`, from T5) is the separate
  larger PI-review item. If it is rejected or revised, Edits 11.1's
  "α² and γ taxonomy" rows and 12.1's γ row should be revised to match.
- Tier 1 edits from `docs/claude_md_proposed_updates.md` are assumed
  already applied; if not, Tier 1's edits go first and this document's
  Edits 2.1/10.1/11.1 append after Tier 1's blocks rather than
  replacing them.
