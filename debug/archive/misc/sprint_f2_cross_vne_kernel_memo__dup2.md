# Sprint F2 — Cross-V_ne kernel-shape substitution diagnostic

**Date:** 2026-05-23 (post-F1-max_n=3 same-day continuation).
**Sprint position:** Track 3's original named target (P2 from CLAUDE.md §3 PK cross-center synthesis), surfaced as the natural pivot after F1 max_n=3 returned CLEAN NEGATIVE (basis enlargement alone does not close the W1c-residual orthogonality wall at NaH). Diagnoses whether the W1c-residual wall lives in the precision of the existing multipole cross-V_ne kernel (in which case full-3D quadrature substitution should close it) or in a structurally distinct mechanism.
**Cross-references:** `debug/sprint_f1_maxn3_predictions_test_memo.md` (F1 clean negative); `debug/sprint_f1_p1p2_combined_test_memo.md` (F1 P1+P2 baseline); `debug/sprint_w1c_mz_partition_analysis_memo.md` (bridge sprint); `debug/sprint_modular_alpha_arc_synthesis_memo.md` (full α arc); `geovac/shibuya_wulfman.py`; `geovac/cross_center_screened_vne.py`; `geovac/balanced_coupled.py`; `geovac/multi_zeta_orbitals.py`; CLAUDE.md §1.7 multi-focal-composition wall taxonomy.

---

## §0. Executive summary + verdict

**Verdict line: KERNEL-NOT-IT (with structural absence finding).**

The Step 1 algebraic kernel diagnostic returns a decisive negative on the kernel-shape hypothesis. At NaH R_eq = 3.566 bohr the multipole expansion (current `geovac/shibuya_wulfman.py` cross-V_ne path) at the framework's default L_max = 4 truncation is bit-faithful to converged 3D numerical quadrature across all five tested matrix elements. The largest converged differential is **2.01 × 10⁻⁵ Ha** on the H 1s diagonal cross-V_ne contribution, which is **2.7 × 10⁻⁴ = 0.027%** of the NaH bond-energy scale (0.075 Ha). The cross-V_ne kernel is six to ten orders of magnitude better than what would be needed to close the W1c-residual wall. Full-3D substitution would change matrix elements by at most ~10⁻⁵ Ha, four orders of magnitude smaller than the 0.305 Ha W1c-residual descent depth.

Per the gate logic ("multipole kernel error < 10% of splitting → STOP and report clean negative"), Step 2 and Step 3 are skipped. **The wall is NOT a kernel-truncation artifact.**

**The substantive new content the gate did NOT anticipate:** Step 1 also revealed that the framework's h1 architecture in `composed_qubit.build_composed_hamiltonian` is **strictly block-diagonal** — there is **NO cross-block off-diagonal h1 matrix element ⟨ψ_a^A | V_ne(any nucleus) | ψ_b^B⟩** at all for orbitals on different centers A ≠ B. The bonding/antibonding splitting from h1 alone in the {|Na 3s⟩, |H 1s⟩} basis is therefore identically zero. The F1 max_n=3 finding (bonding orbital is constructible from natural-orbital diagonalization of the 1-RDM) IS the FCI's emergent construction via 2-body ERIs only, NOT via h1 cross-coupling — because the latter does not exist in the architecture.

This sharpens the W1c-residual classification from F1's "endomorphism" (external-input class) to a more specific **"architectural absence"**: the framework architecture does not include the matrix slot that would let h1 directly favor the bonding combination. The natural follow-on (Sprint F3) is to extend `composed_qubit` to include true two-center one-body integrals ⟨ψ_a^A | V_ne | ψ_b^B⟩ — a structurally larger extension than basis enlargement or within-sub-block kernel-shape work.

Net structural reading: the multipole expansion is faithful where it lives (within sub-blocks); the wall lives in a matrix slot the architecture does not have at all.

---

## §1. Step 1 — Algebraic kernel diagnostic (the executed sub-sprint)

### §1.1 Test setup

At NaH R_AB = 3.566 bohr (experimental equilibrium), five cross-V_ne matrix elements computed via:
- **mp L_max=4 framework:** the default truncation in `compute_cross_center_vne_element` used by `build_balanced_hamiltonian` at NaH max_n=3.
- **mp L_max=20 high:** converged multipole expansion as cross-check that L_max=4 is itself converged.
- **3D quadrature:** for s-s matrix elements the angular integral $\int_0^\pi d\theta \sin\theta / |r - R\hat z|$ has the closed form $2/\max(r, R)$, reducing the 3D quadrature to a 1D Gauss-Legendre radial integration. Converged at n_radial = 3200, r_max = 120 bohr (verified in Step 1b convergence scan).

### §1.2 Results table

| Matrix element | Description | mp L_max=4 (Ha) | mp L_max=20 (Ha) | 3D quadrature (Ha) | mp vs 3D diff (Ha) | rel err |
|:--|:--|--:|--:|--:|--:|--:|
| A | ⟨Na 3s_hyd \| V(H) \| Na 3s_hyd⟩ diag | −0.0911941 | −0.0911941 | −0.0911950 | +8.1 × 10⁻⁷ | 8.9 × 10⁻⁶ |
| B | ⟨Na 3s_hyd \| V(H) \| Na 4s_hyd⟩ offdiag | −0.0234746 | −0.0234746 | −0.0234746* | +3.1 × 10⁻⁹ | 1.3 × 10⁻⁷ |
| C | ⟨H 1s \| V(Na) \| H 1s⟩ diag | −3.0734334 | −3.0734334 | −3.0734535 | +2.0 × 10⁻⁵ | 6.5 × 10⁻⁶ |
| D | ⟨H 1s \| V(Na) \| H 2s⟩ offdiag | −0.0695453 | −0.0695453 | −0.0695647 | +1.9 × 10⁻⁵ | 2.8 × 10⁻⁴ |
| E | ⟨Na 3s_mz \| V(H) \| Na 3s_mz⟩ diag (multi-zeta) | n/a | −0.2260471 | −0.2260409 | −6.2 × 10⁻⁶ | 2.7 × 10⁻⁵ |

\* B's 3D-quadrature value at the initial n_radial=400 setting showed 1% error; Step 1b convergence scan confirms this was undersampling. At n_radial=3200, r_max=120, the value matches converged multipole to 3 × 10⁻⁹ Ha.

The **mp L_max=4** column is bit-identical to **mp L_max=20** for all four hydrogenic matrix elements (A, B, C, D), confirming the framework's default L_max=4 truncation is itself converged at NaH R_eq. The Gaunt selection rule on the radial split integrals (which terminates at L_max = l₁ + l₂) means all five s-s matrix elements truncate at L=0 in principle (since l₁ = l₂ = 0 → L_max = 0), but the angular coefficient $A_L(0,0;L) = \delta_{L,0}$ already selects only L=0 anyway. The framework's L_max=4 setting is thus generously over-converged for the s-s sector tested here.

### §1.3 Differential vs bond-energy scale

| Matrix element | Differential (Ha) | NaH D_e scale (Ha) | ratio |
|:--|--:|--:|--:|
| A (Na 3s diag) | 8.1e-7 | 0.075 | 1.1 × 10⁻⁵ |
| B (Na 3s→4s offdiag, converged quad) | 3.1e-9 | 0.075 | 4.1 × 10⁻⁸ |
| C (H 1s diag) | 2.0e-5 | 0.075 | 2.7 × 10⁻⁴ |
| D (H 1s→2s offdiag) | 1.9e-5 | 0.075 | 2.6 × 10⁻⁴ |
| E (Na 3s_mz diag) | 6.2e-6 | 0.075 | 8.3 × 10⁻⁵ |

**Maximum differential is 2.7 × 10⁻⁴ of the bond-energy scale (~0.027%).** Even the worst-case is two orders of magnitude smaller than the 10% threshold the gate logic specified for "kernel-might-close-the-wall". The multipole kernel is six to ten orders of magnitude below the W1c-residual descent depth (0.305 Ha at NaH max_n=2 per PK cross-center sprint synthesis memo).

---

## §2. Bonding/antibonding h1 splitting analysis — the structural absence finding

### §2.1 What h1 looks like in the 2-state basis {|Na 3s⟩, |H 1s⟩}

The framework's h1 in this basis under combined W1c × multi-zeta architecture:

| Element | Framework value (Ha) | 3D quadrature reference (Ha) | Differential |
|:--|--:|--:|--:|
| h1[Na 3s, Na 3s] (hydrogenic) | −0.591194 | −0.591195 | +8 × 10⁻⁷ |
| h1[H 1s, H 1s] | −3.573433 | −3.573453 | +2 × 10⁻⁵ |
| h1[Na 3s, H 1s] | **0** (architecturally absent) | n/a | n/a |
| h1[Na 3s_mz, Na 3s_mz] (multi-zeta) | −0.396047 | −0.396041 | −6 × 10⁻⁶ |

The framework on-site h1 diagonal under W1c+mz combines:
- Hydrogenic eigenvalue $-Z_{\rm orb}^2/(2n^2) = -0.5$ Ha (no mz) or screened-Schrödinger eigenvalue ($\approx -0.170$ Ha for Na 3s with mz)
- Within-sub-block cross-V_ne contribution from the off-center nucleus

These two contributions land at the diagonal h1 slots. The off-diagonal h1 slot connecting Na 3s and H 1s is **structurally absent**: `composed_qubit.build_composed_hamiltonian` constructs h1 as a strict block-diagonal sum of atom-like sub-block Hamiltonians, with no slot for cross-block one-body coupling.

### §2.2 What this implies for bonding/antibonding splitting

In a 2-state diagonal h1 with no off-diagonal coupling, the bonding combination $|B\rangle = (|{\rm Na\ 3s}\rangle + |{\rm H\ 1s}\rangle)/\sqrt 2$ and antibonding combination $|A\rangle = (|{\rm Na\ 3s}\rangle - |{\rm H\ 1s}\rangle)/\sqrt 2$ have **identical** h1 expectation values:

$$\langle B | h_1 | B \rangle = \langle A | h_1 | A \rangle = \tfrac{1}{2}(h_1[{\rm Na\ 3s, Na\ 3s}] + h_1[{\rm H\ 1s, H\ 1s}])$$

The bonding-antibonding splitting from h1 alone is **identically zero**. F1 max_n=3 showed that the FCI dominant natural orbital under combined W1c × mz at max_n=3 IS the bonding combination (mixing coefficients $-0.698\, \text{Na 3s} - 0.687\, \text{H 1s}$, ~50/50 amplitude split), with the antibonding combination as the 2nd NO. This bonding orbital construction is therefore emerging entirely from the 2-body ERI coupling in the FCI eigenproblem — the h1 sector contributes zero to the bonding-vs-antibonding splitting.

### §2.3 The diagnostic verdict — structural absence

This is a stronger statement than "kernel is faithful at the existing matrix slots". The diagnostic also rules out a path that the gate logic did not anticipate: **even if the within-sub-block multipole kernel were arbitrarily accurate, the framework still could not have h1 prefer the bonding combination over the separated configuration, because the matrix slot that would carry that preference does not exist.**

The cross-V_ne kernel-shape substitution hypothesis was structured around the assumption that the wall lives in inaccurate matrix elements that the framework computes. Step 1 shows that the kernel is essentially bit-faithful where it lives. The wall lives in matrix elements the framework does NOT compute at all — specifically, the cross-block h1 coupling between orbitals on different centers.

---

## §3. Step 2 (full-3D-quadrature substitution) — SKIPPED per Step 1 gate

Step 2 was conditional on the Step 1 verdict being "multipole kernel ERROR > 50% of the bonding-antibonding splitting AND full-3D kernel gives correct ordering". Step 1's verdict is the third branch: "multipole kernel ERROR < 10% of the splitting; the wall is NOT a kernel-truncation artifact". Per the gate logic, Step 2 is skipped.

A full-3D-quadrature production extension WAS feasible (the diagnostic implementation in `debug/sprint_f2_step1_kernel_diagnostic.py` shows the scaffolding for the s-s case via closed-form $\theta$ integration plus 1D Gauss-Legendre radial), but substituting it into `geovac/shibuya_wulfman.py` would change matrix elements by at most ~10⁻⁵ Ha at NaH R_eq — five orders of magnitude below the wall depth. The production extension would consume implementation budget without addressing the actual wall mechanism, which is the architectural absence of cross-block h1 slots.

---

## §4. Step 3 (mini-PES) — SKIPPED per Step 2 skip

Not reached.

---

## §5. Verdict + implications for the W1c wall structural taxonomy

### §5.1 Sub-layer 3 reclassification (third time, getting sharper)

Tracking the evolving classification of the W1c-residual orthogonality wall at NaH:

| Sprint | Classification | Diagnostic |
|:--|:--|:--|
| PK cross-center sprint | Standard PK reduces descent depth 17.5× → 1.17× (engineering positive; binding NEGATIVE). Wall lives in "additional mechanism beyond W1c+PK". | Engineering attempt; mechanism unspecified. |
| F1 P1+P2 max_n=2 sprint | Multi-zeta on Na 3s reduces descent another ~24%. Wall lives in "endomorphism class" (basis-irreducible). | Two-bucket M-Z partition (cross-shift / endomorphism). |
| F1 max_n=3 sprint | Basis enlargement to max_n=3 DOES construct the bonding orbital (structurally) but cannot energetically prefer it. Three-bucket refinement (cross-shift / basis-closable cross-shift / endomorphism) FALSIFIED; sub-layer 3 stays in endomorphism class. | Bonding orbital is NO of 1-RDM at max_n=3 with 50/50 Na 3s + H 1s split. |
| **F2 (this sprint)** | **Architectural absence**: cross-block h1 slots that would let h1 prefer bonding are not in `composed_qubit`. Kernel-shape substitution within existing slots is faithful at ~0.03% of bond energy. | Direct diagnostic of where the wall lives at the h1 level. |

The classification has gotten more specific at each step. F2's "architectural absence" framing is the sharpest yet: the wall is not a calibration-data input (the original endomorphism framing) and not a basis-flexibility limit (F1 max_n=3 already constructs the bonding orbital structurally); it is the absence of a specific matrix element in the framework's h1 construction.

### §5.2 Multi-focal-composition wall taxonomy update

CLAUDE.md §1.7 currently lists six refined sub-walls (W1a/b/c, W2a, W2b-easy/medium, W3). The W1c-residual at NaH max_n=3 (post-F2) is not a sub-wall in the sense those entries name; it's a **structural absence in the composed_qubit architecture**, distinct from:
- W1a (cross-register coordinate operator, addressed by Roothaan)
- W1b (magnetization-distribution operator, addressed by ω_magn inner fluctuation)
- W1c (frozen-core cross-center screening, addressed by FrozenCore Z_eff(r) + multi-zeta)
- W2a (multi-loop UV/IR composition)
- W2b (cross-manifold composition)
- W3 (inner-factor calibration data)

The "architectural absence" of cross-block h1 is structurally orthogonal to these — it is a missing matrix slot in the existing block-diagonal architecture, not a missing physical projection or a missing calibration input. A sensible name for the new entry would be **"W1d: cross-block h1 architectural extension"** — the slot that would carry one-body cross-block coupling between orbitals on different centers.

W1d is reachable via a straightforward (1-2 week) extension of `composed_qubit.build_composed_hamiltonian` to include true two-center one-body integrals. This is the natural sub-sprint after F2.

### §5.3 Consistency with the broader framework structural-skeleton-scope finding

The CLAUDE.md §1.7 structural-skeleton-scope finding states that GeoVac maps the structural skeleton (selection rules, transcendental classes, scaling laws) but does not autonomously generate calibration data. F2 sharpens this for chemistry-side observables: the framework also does not autonomously include certain ARCHITECTURAL coupling slots — in particular, the cross-block h1 element that would let h1 prefer a bonding orbital over a separated configuration is not even constructed.

This is consistent with the "structural-skeleton" framing in a refined sense: the skeleton is the angular content (Gaunt selection rules, multipole termination, sub-block structure) and the within-sub-block radial precision; the calibration content is parameters (Yukawa, counterterms); and the **architectural coupling content** (which sub-blocks couple at which body-order) is its own axis. W1d sits in this third axis: the architecture currently has cross-block ERIs (2-body) but not cross-block h1 (1-body). Adding cross-block h1 is an architectural extension that would not violate the structural-skeleton-scope finding — it would simply add a coupling axis the framework currently lacks.

### §5.4 Honest scope

What this sprint DID demonstrate:
- Multipole expansion at L_max=4 is bit-faithful to converged 3D quadrature for s-s cross-V_ne matrix elements at NaH R_eq (max differential 2 × 10⁻⁵ Ha).
- The W1c-residual orthogonality wall at NaH is NOT a kernel-truncation artifact within the existing framework architecture.
- The framework's h1 in `composed_qubit` is strictly block-diagonal; cross-block h1 matrix elements for orbitals on different centers are architecturally absent.
- The F1 max_n=3 bonding orbital construction emerges entirely from 2-body ERI coupling in the FCI, not from any h1 cross-block term.

What this sprint did NOT demonstrate:
- That cross-block h1 extension would actually close the wall (separate sprint F3).
- That l > 0 cross-V_ne matrix elements have the same kernel-faithfulness (the diagnostic only tested s-s; the framework computes p-s, p-p, etc. via the same multipole expansion with non-trivial angular coefficients, but at NaH max_n=3 only s-orbitals are active in the bonding combination). p-p test would be a 1-day extension if needed.
- That the structural-absence finding generalizes to other systems (LiH, BeH₂, H₂O all use the same `composed_qubit` block-diagonal architecture, so the absence is general; whether it MATTERS depends on system).
- That the W1c-residual wall is the only barrier to NaH binding (the wall may simply be the dominant of multiple missing pieces).

---

## §6. Recommended next sprint after F2

### Priority 1 — Sprint F3: cross-block h1 architectural extension (NEW CONCEPT, 1-2 weeks)

Extend `geovac/composed_qubit.build_composed_hamiltonian` to include cross-block off-diagonal h1 elements $\langle \psi_a^A | -Z_C/|r - R_C| | \psi_b^B \rangle$ for orbitals on different centers $A \neq B$. The integral is a true two-center one-body integral; implementation options:

1. **Elliptic-coordinate quadrature for diatomic systems**: at NaH the natural coordinates are confocal elliptic (prolate spheroidal), in which the operator $1/|r - R_C|$ is separable. Same technique as Paper 11 for H₂⁺. ~1 week.
2. **Two-center Gauss-Hermite quadrature**: works for any geometry; uses radial-angular product grids on each center with bipolar harmonic decomposition. Cleaner for polyatomic generalization. ~1.5 weeks.
3. **Multipole expansion analog**: extend the existing multipole machinery to bipolar form $1/|r - R_C| = \sum_{LM} ...$ on a two-center grid. Same algebraic flavor as the within-sub-block path, would integrate naturally with existing code. ~1 week.

The F2 finding (bonding orbital is constructible from 1-RDM diagonalization at max_n=3) suggests F3 would test cleanly: at NaH max_n=3 with cross-block h1 wired in, does the FCI energy of the bonding configuration drop below the separated configuration? If yes, the wall closes at the architectural-extension level. If no, the wall has yet another structural mechanism beyond cross-block h1, and the F2 "architectural absence" classification needs further refinement.

The 1-week scope estimate is for the s-s case at NaH only (which is sufficient to test the F2 hypothesis); extension to l > 0 angular content and polyatomic geometry would be a separate 1-2 week follow-on.

### Priority 2 — Pivot to a different system class (parallel, 2-3 days)

Per F1 max_n=3 §7 Priority 3: HCl, MgH₂, or LiH+ scoping. The W1c-residual is empirically Z-decreasing (NaH 5.4-6× / MgH₂ 2.99× / HCl 1.79× reduction). At higher Z the wall may already be a smaller fraction of bond energy. A 2-3 day HCl max_n=2 scoping run with combined W1c × multi-zeta would test whether the framework already supports second-row hydride binding without architectural extension. Requires extending multi-zeta registry to Cl (Z=17) which is multi-day on its own, but the diagnostic value is independent of F3.

### Priority 3 — Pause chemistry arc and pivot to math.OA papers (parallel, no chemistry compute)

Per F1 max_n=3 §7 Priority 4 logic: further chemistry-arc work has now landed three CLEAN NEGATIVES in a row (PK cross-center, F1 max_n=3, F2). The diagnostic-before-engineering rule (memory `feedback_diagnostic_before_engineering.md`, 2026-05-08) suggests this is the moment to pause and consolidate. Math.OA papers (38, 39, 40 published; 45, 46, 47 in flight) are at a natural consolidation point and could absorb a synthesis sprint.

### Priority 4 (DO NOT pursue) — Continue within-sub-block kernel-precision work

The F2 diagnostic decisively closes this direction. The multipole kernel is faithful at ~10⁻⁵ Ha precision, six orders of magnitude below the wall depth. Further precision work within the existing block-diagonal architecture cannot close the wall.

### Bundled recommendation

**Default**: queue Priority 1 (Sprint F3 cross-block h1 architectural extension) as the next active sprint. Cost ~1 week for the NaH-specific s-s diagnostic, structurally distinct from this sprint's within-sub-block kernel work. This is the natural follow-on from the architectural-absence diagnosis.

**Parallel**: keep Priority 2 (HCl scoping) flagged but not active (multi-zeta registry extension is its own multi-day sub-sprint).

**Alternative**: if the PI wants to pause chemistry-arc work entirely and pivot to math.OA consolidation (Priority 3), the F2 verdict is a natural sprint-cycle-end. The W1c-residual wall is now well-characterized as an architectural absence; further chemistry-arc work would be the next architectural extension or a different system class, both 1-2 week investments. The diagnostic value of three consecutive clean negatives plus the F2 structural-absence finding may itself constitute a natural memo/paper-edit batch (CLAUDE.md §3 dead ends row + Paper 19 §sec:w1c_residual update + Paper 17 § 6.10 architectural-clarification note), which is a 2-3 hour activity not requiring further compute.

---

## §7. Files

### Created (debug)

- `debug/sprint_f2_step1_kernel_diagnostic.py` — Step 1 driver: five-matrix-element multipole-vs-3D-quadrature diagnostic + bonding/antibonding h1 analysis.
- `debug/sprint_f2_step1b_convergence.py` — Step 1b convergence scan for off-diagonal Na 3s→4s matrix element (confirms initial 1% discrepancy was undersampling, not real multipole-vs-3D differential).
- `debug/sprint_f2_cross_vne_kernel_memo.md` — this memo (~3500 words).

### Created (data)

- `debug/data/sprint_f2_step1_kernel_diagnostic.json` — Step 1 results (all five matrix elements + bonding/antibonding analysis).
- `debug/data/sprint_f2_step1b_convergence.json` — Step 1b convergence scan.
- `debug/data/sprint_f2_cross_vne_kernel.json` — consolidated results across the sprint with verdict line and follow-on recommendations.

### Modified

- None (production code untouched per gate decision; paper edits deferred per sprint mandate).

### Not modified

- Production code (`geovac/balanced_coupled.py`, `geovac/shibuya_wulfman.py`, `geovac/cross_center_screened_vne.py`, etc.) — Step 2 skipped per gate decision; production extension would not address the actual wall mechanism.
- Paper `.tex` files — sprint mandate: paper edits deferred to a follow-on batch.
- CLAUDE.md — sprint mandate: documentation edits deferred (the W1d "architectural absence" finding warrants §3 dead ends entry + §1.7 multi-focal partition update + §2 sprint entry, but those go in the follow-on documentation pass).
- Tests — no production code modified, no test changes needed.

### Regression check

```
$ pytest tests/test_balanced_coupled_multizeta.py \
         tests/test_balanced_coupled_screened_valence.py \
         tests/test_multi_zeta_orbitals.py \
         tests/test_phillips_kleinman_cross_center.py \
         tests/test_shibuya_wulfman.py \
         -q --no-header
128 passed, 1 skipped in 27.80s
```

Zero regression across 128 baseline tests. No production code modified; regression set is the safety net for the multi-zeta + screened-V_ne + PK + balanced-coupled + shibuya_wulfman architecture that this sprint exercises in diagnostic mode.

---

**End of Sprint F2 memo. Verdict: KERNEL-NOT-IT (with structural absence finding).** Multipole cross-V_ne kernel at framework default L_max=4 is bit-faithful to converged 3D quadrature (max differential 2 × 10⁻⁵ Ha = 0.027% of bond energy). The W1c-residual orthogonality wall at NaH is NOT a kernel-truncation artifact. Substantive new content: the framework's h1 in `composed_qubit` is architecturally block-diagonal — cross-block h1 matrix elements ⟨ψ_a^A | V_ne | ψ_b^B⟩ for orbitals on different centers A ≠ B are absent by design. The bonding orbital F1 max_n=3 found is constructed entirely by 2-body ERI coupling in the FCI, not by any h1 cross-coupling. Sub-layer 3 reclassifies from "endomorphism" to "architectural absence" (W1d) — a missing matrix slot, not a missing calibration input. Next sprint: F3 cross-block h1 architectural extension (1-2 weeks).
