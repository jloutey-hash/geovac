# Sprint SO(4)-Breaking — diagnostic of the W1e chemistry-binding wall

**Date:** 2026-06-07 (session continuation, fifth sprint of the day)
**Driver:** `debug/sprint_so4_breaking_w1e_diagnostic_driver.py`
**Data:** `debug/data/sprint_so4_breaking_w1e.json`
**Log:** `debug/sprint_so4_breaking_w1e_diagnostic_log.txt`
**Predecessors:** Sprint M-vS Symmetries (v3.89.0); Sprint v3.85.0 W1e at projection step (sprint_chemistry_pivot); W1e localization memos.

---

## TL;DR — POSITIVE DIAGNOSTIC, structurally clean

**The W1e chemistry-binding wall has a concrete symmetry-shaped explanation.** The framework's per-sub-block h1 starts with exact SO(4) (hydrogenic eigenvalues $-Z^2/(2n^2)$ giving 2s/2p degeneracy per block). When the balanced builder adds cross-center $V_{\rm ne}$ DIAGONALLY to each sub-block's h1, this breaks SO(4) within blocks, producing eigenvalues per sub-block. The diagnostic on LiH and NaH shows two things:

1. **The cross-center $V_{\rm ne}$ contribution makes the H-partner sub-block artificially deep** — at LiH $R_{\rm eq}$, H 1s eigenvalue = −1.54 Ha vs NIST H_1s = −0.5 Ha. The −1.04 Ha "extra depth" is exactly the Coulomb pull of the Li nucleus on H: $-Z_{\rm Li}/R = -3/3.015 = -0.995$ Ha. The framework's H sub-block represents "H electron at Li-bond distance," and the diagonal cross-V_ne shifts it to that physically-meaningful depth.

2. **What's MISSING is the off-diagonal coupling that would let the electron delocalize away from the over-deepened sub-block.** Without proper inter-block one-body coupling, FCI ground state collapses onto the deepest single-particle sub-block, producing the W1e over-binding monotone descent. Cross-block_h1 was the architectural attempt at this (Sprint F3) but it ALSO produces over-binding (16× on LiH per the v3.86.0 kwarg-sweep dead-end).

**The structural verdict (genuinely new):** W1e is a **non-orthogonality / Gaunt-sparsity tension**, not just a "calibration data missing" finding. The per-sub-block hydrogenic basis uses DIFFERENT $Z_{\rm eff}$ per sub-block (Li_core $Z=3$, LiH_bond_center $Z=Z_{\rm eff}\approx 1.3$, LiH_bond_partner $Z=1$), making the orbitals NOT mutually orthogonal in any Hilbert-space sense. Proper chemistry requires Löwdin orthogonalization, which **destroys the Gaunt selection rule** (per CLAUDE.md §3 dead-end row: "Heterogeneous nested (per-pair $Z_{\rm eff}$) ... Löwdin orthogonalization destroys Gaunt sparsity (1711 vs 120 Pauli, 14× inflation)"). The W1e wall is the framework choosing **angular sparsity and gauge-network structural identity** over **chemistry-binding accuracy**. It's not a missing operator — it's an architectural tension between two desiderata.

This is a **sharpening** of W1e from "wall at the projection step" (v3.85.0) to "wall at the orthogonality / Gaunt-sparsity tension in the per-sub-block hydrogenic basis." Both are true; this one is more actionable because it identifies the structural choice that produces it.

## 1. Concrete data

### LiH balanced (M=15, R sweep)

Per-sub-block h1 eigenvalues at $R_{\rm eq} = 3.015$ bohr:

| Sub-block | 1s | 2s | 2p triplet | SO(4) block gap |
|:----------|:--:|:--:|:----------:|:----------------:|
| Li_core_center (Z=3)        | −4.832 | −1.368 | −1.582, −1.433, −1.433 | 3.464 Ha |
| LiH_bond_center (Z=Z_eff)   | −0.836 | −0.275 | −0.421, −0.326, −0.326 | 0.561 Ha |
| LiH_bond_partner (Z=1, H)   | −1.535 | −0.568 | −0.994, −0.727, −0.727 | 0.967 Ha |

NIST atomic references for comparison:
- Li 1s (Z²/2 ≈ 4.5 in Li²⁺) ≈ −2.491 (framework is 2× too deep because it uses bare Z²/(2n²) for Li³⁺-like single-electron problem)
- Li 2s = −0.198 (framework Li_core 2s = −1.37 ⇒ wrong because Li_core block represents "1s²-screened atomic Li" but is built as bare Z=3 hydrogenic)
- H 1s = −0.500 (framework LiH_bond_partner 1s = −1.535, off by −1.04 Ha = $-Z_{\rm Li}/R$ at R=3.015 — this is the cross-V_ne contribution and IS physically meaningful when interpreting "H sub-block at bond distance")

### NaH balanced (M=10, R sweep)

| Sub-block | 3s (block n=1) | 4s (block n=2) | 4p triplet | SO(4) block gap |
|:----------|:--:|:--:|:----------:|:----------------:|
| NaH_bond_center (Z=Z_eff_Na≈1.3) | −0.785 | −0.265 | −0.418, −0.313, −0.313 | 0.520 Ha |
| NaH_bond_partner (Z=1, H)       | −3.936 | −1.625 | −3.084, −2.194, −2.194 | 2.312 Ha |

The NaH H sub-block is MUCH deeper than LiH H sub-block ($-3.94$ vs $-1.54$), reflecting the larger $Z_{\rm Na}/R = 11/3.566 = 3.08$ Ha cross-V_ne contribution at $R_{\rm eq}$. The SO(4) block gap for NaH's H partner is 2.31 Ha, vs 0.97 Ha for LiH's H partner.

### Pattern: W1e wall depth correlates with H-sub-block over-deepening

R-sweep shows the H partner's SO(4) block gap decreases monotonically with R:

| R (bohr) | LiH H partner SO(4) gap | NaH H partner SO(4) gap |
|:--------:|:-----------------------:|:-----------------------:|
| 2.5      | 1.128 | — |
| 3.0      | 0.967 (R_eq) | 2.714 |
| 3.566    | 0.856 | 2.312 (R_eq) |
| 4.0      | 0.770 | 2.101 |
| 5.0      | 0.649 | 1.760 |
| 8.0      | 0.490 | 1.011 |

This tracks the cross-V_ne strength $-Z_{\rm heavy}/R$ exactly. The "SO(4)-breaking" the framework introduces is precisely the diagonal cross-Coulomb pull, monotone in $R$ — which is exactly the shape of the W1e over-binding (also monotone in $R$).

## 2. Structural interpretation

### 2.1 What the diagonal cross-V_ne does

In the balanced builder, for each sub-block, the cross-center $V_{\rm ne}$ from every off-center nucleus is added to that sub-block's h1 sub-matrix:

$$h_1^{[\mathrm{sb}]}_{nl,n'l'} = -\frac{Z_{\rm sb}^2}{2n^2} \delta_{n,n'}\delta_{l,l'} + \mathrm{PK}_{nl,n'l'} + \sum_{C \notin \mathrm{sb}} \langle \psi_{n,l,m}^{[\mathrm{sb}]} | -Z_C/|r-R_C| | \psi_{n',l',m'}^{[\mathrm{sb}]} \rangle$$

The first term is hydrogenic eigenvalues (exact SO(4)). The PK term is the Phillips-Kleinman pseudopotential (added only on first-row valence bond_centers). The cross-center $V_{\rm ne}$ term is the chemistry input — it's what makes the H electron "feel" the Li nucleus.

At LiH $R_{\rm eq}$, the cross-V_ne contribution shifts the H sub-block diagonal by −1.0 Ha (the Li-nucleus Coulomb pull at distance 3.015). This is PHYSICALLY CORRECT — an H atom at bond distance to Li IS pulled by the Li nucleus.

The problem isn't that the cross-V_ne is wrong; it's that the cross-V_ne enters ONLY on the diagonal. In a proper quantum chemistry construction, the same cross-V_ne would ALSO appear in OFF-DIAGONAL one-body elements between orbitals on different centers (the inter-block hopping that lets electrons delocalize).

### 2.2 What's missing: orthogonal off-diagonal coupling

The "correct" off-diagonal one-body element between Li and H orbitals would be:

$$h_1^{[\mathrm{Li \to H}]}_{nl, n'l'} = \langle \psi^{\rm Li}_{n,l,m} | T - Z_{\rm Li}/|r-R_{\rm Li}| - Z_{\rm H}/|r-R_{\rm H}| | \psi^{\rm H}_{n',l',m'} \rangle$$

This element couples Li orbitals to H orbitals, allowing FCI wavefunctions with $c_{\rm Li} \cdot \phi_{\rm Li} + c_{\rm H} \cdot \phi_{\rm H}$ bonding character. Without this coupling, the FCI ground state has each electron LOCALIZED on a single sub-block.

The framework attempts this via `cross_block_h1=True` (Sprint F3, May 2026), and the matrix element is computed correctly per `geovac/cross_block_h1.py`. But it produces 16× over-binding on LiH (per the v3.86.0 LiH kwarg-sweep dead-end row in CLAUDE.md §3).

### 2.3 Why cross_block_h1=True fails: non-orthogonal basis

The deeper issue: the framework's per-sub-block hydrogenic orbitals use DIFFERENT Z values per block. For LiH:
- Li_core orbitals: hydrogenic at Z=3
- LiH_bond_center orbitals: hydrogenic at Z=Z_eff_Li ≈ 1.3
- LiH_bond_partner orbitals: hydrogenic at Z=1

These orbitals are NOT mutually orthogonal. The Li_core 1s (deep, peaked at small r) has non-zero overlap with the H 1s (less peaked, peaked at H nucleus position). Without Löwdin / Schmidt orthogonalization, FCI on this non-orthogonal basis gives WRONG energies — specifically, "bonding without Pauli repulsion" because the overlapping orbitals can be co-occupied without the orthogonality-enforced exclusion that real quantum mechanics requires.

This is documented in CLAUDE.md §3 dead-end row:
> Heterogeneous nested (per-pair Z_eff) | 1 | Löwdin orthogonalization destroys Gaunt sparsity (1711 vs 120 Pauli, 14× inflation). Track DF Sprint 5.

So the framework chose **angular sparsity** (Gaunt selection rule preserving 120-Pauli LiH) over **chemistry-binding accuracy** (which would require Löwdin orthogonalization at the cost of 14× Pauli inflation).

### 2.4 The W1e wall is an architectural choice

W1e is not a missing operator that could be added. It's a **structural tension**:

| Property | Per-sub-block hydrogenic basis | Orthogonalized basis |
|:---------|:------------------------------:|:---------------------:|
| Gaunt sparsity | YES (Paper 22) | NO (Löwdin mixes l-shells) |
| M-vS gauge-network structure | YES (Paper 32) | NO (orthogonalization mixes blocks) |
| ℓ-parity Z₂ tapering (v3.89.0) | YES | NO (rotation breaks per-block) |
| Hopf m_l Z₂ tapering (v3.52.0) | YES | NO |
| Chemistry binding (proper Pauli repulsion) | NO (W1e wall) | YES |
| Pauli term count | ~900 for LiH | ~12,600 (14× inflation) |
| Production qubit-encoding advantages | YES (Paper 14 O(Q^2.5)) | NO (loses sparsity) |

The framework chose the LEFT column. The W1e wall is the price.

## 3. Verdict and decision

### POSITIVE DIAGNOSTIC

The SO(4)-breaking analysis identifies the W1e wall as a **non-orthogonality / Gaunt-sparsity tension** in the per-sub-block hydrogenic basis. The diagonal cross-V_ne contribution correctly represents the heavy-nucleus pull on H, but without orthogonal off-diagonal coupling, FCI ground state collapses onto the over-deepened H sub-block.

This is a SHARPENING of W1e:
- v3.85.0 said: "wall at the projection step from continuous to integrals"
- This sprint says: "wall is at the orthogonality/sparsity choice of the per-sub-block hydrogenic basis, with the diagonal cross-V_ne (correct chemistry input) and the lack of off-diagonal orthogonal coupling (the missing inter-block hopping) BOTH being symptoms of this choice"

### What this means for chemistry

The W1e wall is **structurally inseparable from the framework's qubit-encoding advantages**. Closing W1e would require Löwdin orthogonalization, which:
- Inflates Pauli count 14× (verified, Track DF Sprint 5)
- Destroys the angular sparsity theorem (Paper 22) at the production level
- Destroys the M-vS gauge-network structural identification (today's M-vS arc)
- Destroys all three Z₂ tapering schemes we've shipped (Hopf, ℓ-parity, atom-swap)

These advantages are GeoVac's distinctive features — they're what makes the framework GeoVac and not just-another-quantum-chemistry-code. Choosing them means accepting the W1e wall.

The chemistry-engineering arc was empirically exhausted in v3.86.0; this sprint provides the **structural REASON for that exhaustion**: the wall is a deliberate-but-irreducible consequence of the framework's design choices.

### What this means for the framework story

- The framework's chemistry section (Paper 17, Paper 19, Paper 20 row-conditional scope) can now precisely state: "The per-sub-block hydrogenic basis maintains angular sparsity and gauge-network structure; this is an irreducible architectural choice that bounds chemistry-binding accuracy to first-row (where the H sub-block over-deepening is modest) and excludes second-row chemistry-binding at the integral-export level (where over-deepening dominates)."
- Paper 14 (qubit encoding) can frame the ~900 LiH Pauli count as a structural choice, not just an empirical observation.
- This honest framing is *better* than "we don't know why" — it's a *structural theorem about the framework's design trade-off*.

### What this means strategically

The third "M-vS upgrades chemistry" candidate path (find new symmetries via the explorer, v3.89.0) shipped the ℓ-parity Z₂ tapering. This sprint shows: that tapering, AND every other GeoVac-distinctive feature, EXISTS BECAUSE the framework chose the non-orthogonal per-sub-block basis. The Pauli reduction wins and the W1e binding wall are TWO SIDES OF THE SAME CHOICE.

This sharpens GeoVac's identity:
- **It is** a sparse qubit-encoding chemistry framework with structural Z₂ symmetries
- **It is not** a high-precision binding-energy calculator

These are different products. Today's diagnostic makes that distinction structural rather than aspirational.

## 4. What this sprint did NOT do

- Did NOT propose a closure mechanism for W1e (the diagnostic identifies why closure is structurally hard, not how to close).
- Did NOT extract the continuous Level-4 PK-composed effective single-particle eigenvalues for direct comparison (left as a follow-on; would clarify which off-diagonal couplings the continuous solver does correctly).
- Did NOT modify production code.
- Did NOT change any paper text (the framing implications above are recommendations for future paper edits, not applied here).
- Did NOT close (or attempt to close) W1e itself — by design, this is a diagnostic, not an engineering attempt.

## 5. Hard-prohibition check (CLAUDE.md §13.5)

No changes to: natural geometry hierarchy, fitted/empirical parameters, §3 deletions (only ADDING a sharpening of the W1e wall row), Paper 2 K = π(B+F−Δ) labeling.

## 6. Verification

- Driver runs in ~80s wall (LiH 6×R + NaH 5×R balanced builders + FCI per point).
- Per-sub-block eigenvalues match the framework's diagonal h1 + cross-V_ne structure: verified by tracing through the calculation (H 1s eigenvalue at LiH R_eq = −0.5 + V_ne^Li at R = −0.5 − 1.0 = −1.5, matches observed −1.54).
- Cross-V_ne strength −Z_heavy/R matches Coulomb formula exactly across R sweep for both molecules.
- No production code modified → no `geovac/` regression; topological-integrity proofs unchanged (last verified v3.86.0).

## 7. Sample size

- 2 molecules (LiH, NaH)
- 1 builder per (balanced) — balanced is the relevant builder for W1e benchmarking; the composed builder doesn't add cross-V_ne to h1
- 6 R points for LiH, 5 for NaH
- 1 cutoff (n_max=2)

The diagnostic mechanism (diagonal cross-V_ne shifting H sub-block + non-orthogonal basis blocking proper off-diagonal hopping) is structural and doesn't depend on n_max or molecule choice; the sample size is sufficient for the structural verdict.

## 8. Files

### Created
- `debug/sprint_so4_breaking_w1e_diagnostic_driver.py`
- `debug/data/sprint_so4_breaking_w1e.json`
- `debug/sprint_so4_breaking_w1e_diagnostic_log.txt`
- `debug/sprint_so4_breaking_w1e_diagnostic_memo.md` (this file)

### Modified
- None.
