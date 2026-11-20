# Dirac-on-S³ Tier 3 Sprint Plan

**PI approved:** 2026-04-15. Direction 1 from Leader brief: γ corrections + Darwin/MV + D² empty-cell probe.
**Target:** hours under algebraic-first (Tier 1/2 precedent).
**Framing:** curiosity-driven — "does the squared Dirac operator on S³ produce the missing 2nd-order × spinor-bundle transcendental for Paper 18's taxonomic grid?" T7+T8 are mechanical accuracy upgrades shipping regardless.

---

## Governing philosophy

Algebraic-first. Paper 18 taxonomy. All worker dispatches carry the standard preamble.

## PI decisions (recorded)

1. **Scope:** Direction 1 (T7 + T8 + T9). Directions 2-4 deferred.
2. **If T9 closes the empty cell:** update both Paper 18 §IV and Paper 24 §V.
3. **Sunaga SI urgency:** can wait. No external pressure.
4. **Fine-structure target:** OoM sufficient for now. Direction 3 (SS/SOO) deferred.
5. **DIRAC build:** not available. Direction 2 planned around hand-extraction when needed.

---

## Tracks (all independent, launch in parallel)

### Track T7 — γ = √(1-(Zα)²) radial corrections

**Goal:** Bind Martínez-y-Romero three-term recursions for Dirac-Coulomb radial matrix elements ⟨n'κ'|r^s|nκ⟩ into `geovac/dirac_matrix_elements.py`. These carry γ as the relativistic correction factor. Currently γ is a reserved symbol; T7 makes it live.

**Success:** ⟨r^s⟩ matrix elements in exact Fraction(Z, α, γ) arithmetic for s ∈ {-3, -2, -1, 0, 1, 2}; match non-relativistic limit when α → 0; match published Dirac-Coulomb ⟨1/r⟩ = Z/(n²γ_nκ) for hydrogen.

**Deliverables:** `geovac/dirac_matrix_elements.py` extension; `tests/test_dirac_matrix_elements.py` new tests; memo.

### Track T8 — Darwin + mass-velocity α⁴ corrections

**Goal:** Add the two remaining α⁴ one-body diagonal corrections to the fine-structure ladder:
- Darwin: H_D = (Zα²/2) · |ψ_nκ(0)|² · δ_{l,0} — nonzero only for s-states.
- Mass-velocity: H_MV = -α²/(8) · ⟨p⁴⟩_{nκ} = -Z⁴α²/(8n⁴) · [8n/(2l+1) - 3] (hydrogenic closed form).

Combined with T2's SO, the full α⁴ fine-structure is: E_FS = E_SO + E_D + E_MV. Published result for hydrogen: E_FS = -(Zα)⁴/(2n⁴) · [n/(j+1/2) - 3/4]. Verify exact match.

**Success:** Full α⁴ fine-structure for hydrogen matches the Dirac formula exactly; fine-structure splittings for He/Li/Be improve from 66-211% (SO only) to estimated 10-30%.

**Deliverables:** `geovac/spin_orbit.py` extension (or new `geovac/fine_structure.py`); tests; memo; updated fine-structure comparison table.

### Track T9 — D² on S³: Paper 18 empty-cell probe

**Goal:** Compute the spectral zeta of the squared Dirac operator on unit S³ and check whether it produces a transcendental in Paper 18's empty (2nd-order × spinor-bundle) cell.

**Key facts already in hand:**
- D eigenvalues: |λ_n| = n + 3/2, with g_n = 2(n+1)(n+2) (full Dirac).
- D² eigenvalues: λ²_n = (n + 3/2)² = (2n+3)²/4.
- Spectral zeta: ζ_{D²}(s) = Σ_{n≥0} g_n · (λ²_n)^{-s} = Σ 2(n+1)(n+2) · [(2n+3)/2]^{-2s}.

**Evaluate at:**
- s = 1: the "first interesting" value (spectral dim / Casimir energy).
- s = 2: the Paper-24-calibration-π analog (where scalar Laplacian produces π²).
- s = 3, 4: for completeness.

**Also compute:** Hopf-equivariant decomposition of D² spectral zeta (Tier 1 D4 analog — decompose over U(1) fiber charges q).

**Success (MAJOR):** ζ_{D²}(s) at some natural s ∈ {1, 2} gives a closed form involving π², π⁴, or ζ(2) that is NOT already accounted for by the scalar (Paper 24 calibration) or Dirac-spectrum (Tier 1 ζ(3)) tiers. This fills the empty cell.

**Success (CLEAN NEGATIVE):** ζ_{D²}(s) reduces to the already-known scalar ∇*∇ spectral zeta (via Lichnerowicz D² = ∇*∇ + R/4 where R=6 is just a constant shift). The empty cell is then documented as "structurally degenerate with the scalar case on round S³" — extending Paper 24's asymmetry analysis.

**Failure:** ζ_{D²}(s) doesn't converge or gives a generic mixture with no clean closure. Unlikely (the sum is explicit).

**Deliverables:** `debug/tier3_t9_squared_dirac.py`; `debug/data/tier3_t9_squared_dirac.json`; `debug/dirac_t9_memo.md`.

---

## Guardrails

- T9 must NOT attempt to re-derive α combination rule K = π(B+F-Δ) from D² (Phase 4H closed those avenues).
- T9 CAN ask "does D² produce a new transcendental tier" — that's the taxonomic question, not the α question.
- T7 must NOT introduce numerical radial quadrature — Martínez-y-Romero gives recursions.
- T8 must NOT attempt multi-electron SS/SOO (Direction 3, deferred).
- Do NOT modify Paper 2 (Tier 1 closed that story).

---

## Papers affected

- **Paper 18** §IV — if T9 positive, new empty-cell closure paragraph.
- **Paper 24** §V — if T9 positive, spinor-Lichnerowicz corollary.
- **Paper 14** §V — updated fine-structure table from T8.
- **Paper 20** — updated accuracy table from T8.
