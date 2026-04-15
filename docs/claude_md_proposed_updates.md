# Proposed CLAUDE.md Updates — Dirac-on-S³ Tier 1 Sprint Closure

**Author:** Track D5 worker agent
**Date:** 2026-04-15
**Status:** PROPOSAL — NOT auto-applied. PM to apply after PI approval per CLAUDE.md §13.5 (framing changes in Sections 1.5, 1.6, 13 require PI approval; Section 2 append, Section 3 append, Sections 10/11 append are PM-allowed but batched here for review).

The edits below are ordered by section of CLAUDE.md. Each edit is stated as a specific "replace X with Y" or "append Z" instruction. The full text of each appended block is given so the PM can copy-paste directly.

---

## §1. Project Identity — version bump

**Edit 1.1** (PM may apply per §13.5 — version number only).

In §1, replace:
```
**Version:** v2.9.2 (April 15, 2026)
```
with:
```
**Version:** v2.10.0 (April 15, 2026)
```

**Rationale:** A completed multi-track sprint (D1–D5) that closes the α structural-decomposition program is a minor version, not a patch. Per CLAUDE.md §9 changelog protocol: "minor versions (x.Y.0) for new features, completed diagnostic arcs, or paper additions."

---

## §2. Current Development Frontier — Phase 4I addition

**Edit 2.1** (PM may apply — Section 2 is PM-allowed).

In §2's long backlog entry on "**α combination rule derivation (Paper 2 → Paper 18)** — PAUSED (April 2026)", at the end of the Phase 4H summary (after "...The SM-running hypothesis is now a documented dead end — see Section 3. Data: `debug/data/track_alpha_sm/`."), append:

```
Phase 4I sprint complete (April 2026, five-track Dirac-on-S³ Tier 1 sprint D1/D2/D3/D4/D5), ALL NEGATIVE with three named closed-form structural obstructions. Tests whether the Dirac sector on S³ lifts B, F, or Δ to a common-generator spectral invariant. (D1) Infrastructure module `geovac/dirac_s3.py` with Camporesi-Higuchi spectrum (|λ_n| = n+3/2, g_n^Dirac = 2(n+1)(n+2)), spinor label dataclass, and bit-exact π-free certificate (51/51 tests in exact sympy rational arithmetic, analog of Paper 24 Bargmann-Segal certificate). (D2) Dirac analog of B = 42: NEGATIVE. The scalar Laplacian has a zero mode at (n=1, l=0) whose Casimir weight l(l+1) vanishes, giving B(m) = m(m-1)(m+1)(m+2)(2m+1)/20 with factor (m-1). The Dirac spectrum is gapped (|λ_n| ≥ 3/2), every level contributes to any truncated trace, and no cumulative Dirac/Weyl sum carries the (m-1) factor. Single-point coincidence at m=3: |λ_{m-1}|·g_{m-1}^Weyl = 42 via (m-1)(m+2)=10, uniquely at m=3 and at no other integer — elementary arithmetic coincidence, not spectral identity. (D3) Dirac analog of F = π²/6: NEGATIVE. At the packing exponent s=d_max=4, D_{g_m^Dirac}(4) = Σ 2m(m+1)/m^4 = 2ζ(2) + 2ζ(3) = π²/3 + 2ζ(3). Apéry's theorem gives ζ(2), ζ(3) Q-linearly independent; isolating F requires hand-subtracting 2ζ(3), equivalent to projecting onto the m² sub-weight (i.e. returning to the scalar case). Hurwitz-regularized form D_dirac^CH(4) = π² - π⁴/12 (closed form, inseparable). F lives on the scalar Fock-degeneracy Dirichlet at d_max, not on the Dirac. **POSITIVE BYPRODUCT:** first appearance of ζ(3) as a structural transcendental in the framework, via the weight-m subchannel of the Dirac degeneracy. (D4) Orduz S¹-equivariant Hopf decomposition: NEGATIVE. The 40 states at n_CH=3 decompose under the Hopf U(1) action as a clean 20⊕20 charge-parity split (4 half-integer charges × 5 mult + 5 integer charges × 4 mult), structurally cleaner than Paper 2's scalar factorization Δ⁻¹ = |λ_3|·N(2) = 8·5 — but does NOT produce B or F. Numerical near-miss Σ |λ_n|²·g_n^Dirac = 378 = 9·B at n_CH≤2 drops at n_CH≤3 (1188, not a multiple of 42), so the factor 9 is a coincidence. Categorically distinct from Phase 4B α-D (scalar Hopf quotient): the Dirac decomposition is a first-order operator decomposition with nontrivial charge-q twists, unavailable in the bosonic sector. Rules out the Hopf-equivariant Dirac as common generator. (D5) Three-tier coincidence formally documented. **STRUCTURAL CONCLUSION:** B, F, Δ have canonical spectral homes (B: scalar Casimir trace at m=3; F: scalar Fock Dirichlet at d_max; Δ: single-level Dirac degeneracy at n=3). The three homes span TWO sectors of S³ (scalar and spinor) and THREE object types (finite trace, infinite Dirichlet, single-level degeneracy). K = π(B + F - Δ) is confirmed as a genuine cross-sector sum with no common spectral origin. NINE mechanisms eliminated across Phases 4B-4I (six sphere-spectral, one arithmetic-Dirichlet for Δ, one SM-running, one Dirac-on-S³). FOUR positive structural identifications (κ↔B Fock link α-C; F = D_{n²}(d_max) α-J; Δ⁻¹ = g_3^Dirac SM-D; ζ(3) natively in Dirac sector D3 — new framework transcendental, NOT part of K). **Recommendation per sprint plan decision table:** α structural-decomposition program closed; Paper 2 §IV rewrite drafted (`docs/paper2_section4_rewrite.tex`, not auto-applied, flagged for PI review); NO Tier 1b proof sprint opened; move to Tier 2 (spin-ful composed qubit encoding, heavy-element relativistic chemistry) which builds on D1 infrastructure but is independent of the α program. Paper 18 taxonomy note flagged: distinction between even-zeta content (second-order operators, Jacobi-θ inversion) and odd-zeta content (first-order operators, half-integer Hurwitz) recommended but NOT applied. Data: `debug/data/dirac_d{2,3,4}_*.json`; `docs/dirac_s3_verdict.md`; `docs/paper2_section4_rewrite.tex`.
```

---

## §3. Approaches That Failed — new row

**Edit 3.1** (PM may apply — Section 3 is PM-allowed for appending).

In the §3 failed-approaches table, append a new row:

| Category | Count | Key Lesson |
|:---------|:-----:|:-----------|
| Dirac-sector lift of Paper 2 α combination rule ingredients B, F, Δ | 3 | Phase 4I Tier 1 sprint (D1-D5, 2026). Three obstructions proved closed-form: (1) B does not lift because the scalar Laplacian zero mode at (n=1,l=0) forces a (m-1) factor in B(m) that no Dirac/Weyl cumulative trace carries; single-point hit |λ_{m-1}|·g_{m-1}^Weyl = 42 at m=3 is the elementary coincidence (m-1)(m+2)=10. (2) F does not lift because g_m^Dirac = 2m²+2m mixes two homogeneity classes, giving D_dirac^Fock(4) = 2ζ(2) + 2ζ(3); Apéry Q-linear independence forbids isolating F without projecting back to the m² sub-weight (scalar case). Hurwitz form: D_dirac^CH(4) = π² - π⁴/12, inseparable. (3) Hopf-equivariant Orduz decomposition reproduces Δ⁻¹ = 40 as 20⊕20 charge-parity split (cleaner than Paper 2's 8·5) but does NOT produce B or F; does NOT factorize as 8·5 (the scalar "8" = \|λ_3\| of -Δ_LB is not a Dirac quantity). Positive byproduct: ζ(3) is the first non-π transcendental in the framework, appearing natively in the Dirac sector at s=d_max via the weight-m subchannel. K = π(B + F - Δ) is now a formally-documented cross-sector coincidence (scalar + spinor) with three spectral homes. NO Tier 1b opened per sprint plan D5 decision table; move to Tier 2 spin-ful composed qubit encoding. |

---

## §10. Validation Benchmarks — Dirac-on-S³ rows

**Edit 10.1** (PM may apply — Section 10 is PM-allowed for appending).

In §10's benchmark table, after the HO two-fermion entries and before the final Speed regression row, append:

| Test | Max Error | Purpose |
|:-----|:---------:|:--------|
| Dirac-on-S³ π-free certificate (Weyl sector) | 0 non-rationals, n_max=6 | Analog of Paper 24 §III Bargmann-Segal certification; every \|λ_n\| is exact sympy Rational (n + 3/2), every g_n^Weyl is positive int ((n+1)(n+2)) |
| Dirac-on-S³ π-free certificate (Dirac sector) | 0 non-rationals, n_max=6 | Every \|λ_n\| is Rational, every g_n^Dirac is positive int (2(n+1)(n+2)) |
| Dirac-on-S³ label generator | exactly g_n labels per level | `spinor_labels_at_n` generates exactly (n+1)(n+2) Weyl or 2(n+1)(n+2) Dirac labels; stronger invariant than bare π-free |
| Dirac Δ⁻¹ identity | g_3^Dirac = 40 exactly | Phase 4H SM-D identity (Δ = 1/(g_3^Dirac)) reproduced in D1 API as `delta_inverse_identity() == (40, Rational(1,40))` |
| Fock ↔ CH convention conversion | invertible at all n | `fock_to_ch(ch_to_fock(n)) == n` for n = 1..10; label compatibility with scalar graph (n_Fock = n_CH + 1) |
| D2 cumulative Dirac trace closed form | exact Rational | Σ g_n^Dirac = (N+1)(N+2)(2N+3)/3 symbolic identity, verified for N = 0..5 |
| D2 |λ_{m-1}|·g_{m-1}^Weyl = 42 at m=3 | (m-1)(m+2) = 10 uniquely | Single-point Dirac/Weyl coincidence with B = 42 occurs only at m=3, sympy exact |
| D3 Dirac Dirichlet at s=4 | 2ζ(2) + 2ζ(3) | `summation(2*m*(m+1)*m**(-4), (m,1,oo)) == pi**2/3 + 2*zeta(3)` symbolic |
| D3 Weyl Dirichlet at s=4 | ζ(2) + ζ(3) | Factor of 2 difference from Dirac, sympy exact |
| D4 Hopf charge partition sum | Σ mult = g_n^Dirac | Sympy-exact at n_CH = 0..5; 40 = 20 half-integer-charge + 20 integer-charge at n_CH = 3 |
| D4 Dirac Hurwitz spectral zeta at s=4 | π² − π⁴/12 | `summation(2*(n+1)*(n+2)/(n+Rational(3,2))**4, (n,0,oo))` symbolic closed form |

All D1-D4 tests together: 51 under `tests/test_dirac_s3.py` plus the symbolic checks inside `debug/dirac_d{2,3,4}_*.py`.

---

## §11. Topic-to-Paper Lookup — new Dirac rows

**Edit 11.1** (PM may apply — Section 11 is PM-allowed for appending).

In §11's topic-to-paper table, append rows for Dirac-on-S³ infrastructure and the Tier 1 sprint results. Group them at the end since they don't fit cleanly into the existing ordering:

| Topic | Paper | Section | Tier |
|:------|:-----:|:--------|:----:|
| Dirac-on-S³ infrastructure (Camporesi-Higuchi spectrum, spinor harmonics) | 2 (§IV rewrite), 23 | — | Conjecture/Applications |
| Camporesi-Higuchi spectrum on S³ | 2, 23 | §IV | Conjecture |
| Spinor spherical harmonics on S³ | 2 | §IV | Conjecture |
| π-free certificate for Dirac-on-S³ | 2 | §IV | Conjecture |
| Orduz Hopf-equivariant Dirac decomposition | 2 | §IV | Conjecture |
| 20⊕20 charge-parity split of Δ⁻¹ | 2 | §IV | Conjecture |
| ζ(3) in the Dirac sector at s = d_max | 2 | §IV.5 (new) | Conjecture |
| Odd-zeta vs even-zeta content (spinor vs scalar Laplacian) | 2 | §IV.5 | Conjecture |
| Three-homes theorem for B, F, Δ | 2 | §IV.1 (new) | Conjecture |
| Cross-sector structural coincidence (K combination rule) | 2 | §IV.6 | Conjecture |
| Dirac Fock rigidity (deferred) | 23 | — | Applications |
| First-order vs second-order spectral operators (π content) | 24 | §IV-V | Core |

---

## Summary of edit classification

| Edit | CLAUDE.md section | PM access per §13.5 | Requires PI review? |
|:----:|:------------------|:-------------------|:-------------------|
| 1.1 | §1 version number | Allowed (version ONLY) | No |
| 2.1 | §2 append Phase 4I summary | Allowed (PM-allowed) | **YES** — PI should confirm the "Phase 4I" framing and the closure language before it's committed |
| 3.1 | §3 append failed-approach row | Allowed (append only) | No |
| 10.1 | §10 append benchmark rows | Allowed | No |
| 11.1 | §11 append topic rows | Allowed | No |

**Recommendation:** Edits 1.1, 3.1, 10.1, 11.1 can be applied by the PM without further PI review (they are mechanical and fall within §13.5's allowed-edit scope). Edit 2.1 is PM-allowed in principle but the language closes an eight-phase research program; the PI should review the summary for tone and for the "move to Tier 2" recommendation before it lands.

The Paper 2 §IV rewrite (`docs/paper2_section4_rewrite.tex`) is the separate, larger PI-review item. If the §IV rewrite is rejected or substantially revised, Edit 2.1 should be revised to match.

---

## Notes for PM

- Do NOT apply any of these edits until the PI has reviewed `docs/dirac_s3_verdict.md` and `docs/paper2_section4_rewrite.tex`.
- Apply edits in order 1.1 → 2.1 → 3.1 → 10.1 → 11.1 (the version bump should land in the same commit as the Phase 4I summary, and the rest can follow).
- The Paper 18 taxonomy note (even-zeta vs odd-zeta content) is flagged in the verdict document and the §IV rewrite but is NOT part of these CLAUDE.md edits. A Paper 18 edit is a separate PI-review item if the framing is accepted.
