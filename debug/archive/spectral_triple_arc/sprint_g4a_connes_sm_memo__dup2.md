# Sprint G4a: Connes Standard Model on T_{S³}

**Date:** 2026-05-31
**Verdict:** POSITIVE-THIN
**Version:** v3.33.0

## §1. Goal

Extend Sprint H1's electroweak slice (A_F = ℂ ⊕ ℍ, dim H_F = 8) to the
full Connes-Chamseddine-Marcolli Standard Model finite algebra:

    A_F = ℂ ⊕ ℍ ⊕ M₃(ℂ)

with dim H_F = 32 (16 matter + 16 antimatter), and verify:
1. All load-bearing Connes axioms at finite n_max
2. Inner fluctuations decompose into U(1) × SU(2) × SU(3) + Higgs
3. The falsifier (does GeoVac force Yukawa = 0?) does NOT fire

## §2. Architecture

**Hilbert space:** H_F^matter = H_lepton (4D) ⊕ H_quark (12D = 4 flavor × 3 color).
Full H_F = matter ⊕ antimatter = ℂ^32.

**Ordering convention:** Leptons (0:4), quarks in flavor⊗color layout (4:16) with
color running faster. Within each sector: (particle_L, antiparticle_L, particle_R,
antiparticle_R) matching the H1 convention.

**Algebra action** π(λ, q, m) on matter:
- Leptons: block_diag(q, diag(λ, conj(λ)))  [same as H1]
- Quarks: block_diag(q, diag(λ, conj(λ))) ⊗ m  [M₃ acts on color]
- M₃(ℂ) acts trivially on leptons (color singlets)

**D_F:** block_diag(M_lepton, kron(M_quark_flavor, I₃)) on matter, conjugate on
antimatter. Lepton Yukawa Y_ℓ = diag(y_ν, y_e), quark Yukawa Y_q = diag(y_u, y_d).

**J_F:** matter↔antimatter swap with complex conjugation. J_F² = +I (KO-dim 6).

**γ_F:** L = +1, R = -1 on matter; sign-flipped on antimatter. {J_F, γ_F} = 0.

**Combined:** D = D_GV ⊗ 1_F + γ_GV ⊗ D_F. KO-dim = 3 + 6 = 9 ≡ 1 (mod 8).

## §3. Results

### Connes axioms (bit-exact at n_max ∈ {1, 2})

| Axiom | n_max=1 | n_max=2 |
|:------|:--------|:--------|
| J² = -I | 0 | 0 |
| JD = +DJ | 0 | 0 |
| D Hermitian | 0 | 0 |
| [a, JbJ⁻¹] = 0 (order-zero) | 0 | 0 |
| [[D,a], JbJ⁻¹] = 0 (order-one) | 0 | 0 |

{γ, D} ≠ 0 is expected (truthful CH is chirality-diagonal; KO-dim 1 is odd,
no chirality grading axiom). Same documented behavior as H1.

### Gauge group census (n_max=2)

| Sector | Present? |
|:-------|:---------|
| SU(2) lepton | ✓ |
| SU(2) quark | ✓ |
| SU(3) color | ✓ |
| U(1) hypercharge | ✓ |
| Lepton-quark mixing | ✗ (correct) |
| **Full gauge group** | **U(1) × SU(2) × SU(3)** |

At n_max=1, gauge is trivial (only one Fock shell → all multipliers ∝ I → [D_GV, M] = 0).

### Falsifier

| Config | Higgs zero? | Interpretation |
|:-------|:------------|:---------------|
| Y ≠ 0 (imposed Yukawa) | No (higgs/gauge = 0.003) | Construction works |
| Y = 0 (zero Yukawa) | Yes | Correctly trivial |

### H1 consistency

SM lepton-sector Yukawa and chirality match the H1 electroweak triple bit-exactly.
The SM extension adds quark content without modifying the lepton sector.

## §4. Inner fluctuation structure

omega = Σ a_i [D, b_i] decomposes as:

**Gauge piece:** Σ a_GV [D_GV, b_GV] ⊗ (a_F b_F)
- Fiber content a_F b_F is block-diagonal in L/R → gauge 1-forms
- Lepton gauge: SU(2)_L from ℍ + U(1)_Y from ℂ
- Quark gauge: SU(2)_L ⊗ SU(3)_c from ℍ ⊗ M₃(ℂ) + U(1)_Y ⊗ SU(3)_c

**Higgs piece:** Σ a_GV b_GV γ_GV ⊗ a_F [D_F, b_F]
- [D_F, b_F] has off-diagonal L↔R content from the Yukawa
- Non-zero iff Y ≠ 0 (the falsifier criterion)

**Matter-antimatter:** Decoupled (off-diagonal = 0). Automatic from π_F acting
on matter only + J_F conjugation structure.

**Lepton-quark:** Decoupled (off-diagonal = 0). M₃(ℂ) acts as identity on
leptons, preventing cross-sector mixing.

## §5. Structural interpretation

G4a is the G2-G3 corollary generalized to the full SM algebra:

- **G2 (Sprint H1):** No autonomous Yukawa selection
- **G3 (Sprint G3):** γ_GV and γ_F are independent ℤ₂'s
- **G4a (this sprint):** These are ONE structural fact at the SM level. The
  Yukawa lives in the off-diagonal D_F block that flips γ_F. Since γ_F is
  independent of γ_GV, no GeoVac-side constraint determines Y. Adding M₃(ℂ)
  introduces SU(3) color but does not change this structural fact.

The construction is **SM-consistent, not SM-selecting**. The framework admits
the SM gauge content given imposed Yukawa, hypercharge, and generation choices.
The GeoVac data is silent on each of those choices.

## §6. Honest scope

**Theorem-grade (bit-exact at finite n_max):**
- All five Connes axioms (J², JD, Hermiticity, order-zero, order-one)
- Finite-triple axioms (J_F², {J_F, γ_F}=0, {γ_F, D_F}=0)
- Matter-antimatter decoupling
- Lepton-quark decoupling
- H1 lepton-sector consistency

**Verified by census (numerical, n_samples=30):**
- Full U(1) × SU(2) × SU(3) gauge group identification
- SU(3) color content via Gell-Mann projection

**Structural observation:**
- POSITIVE-THIN verdict (G2-G3 corollary generalized)
- Higgs non-trivial iff Yukawa imposed

**Named open follow-ons:**
1. Multi-generation extension (3 generations with CKM/PMNS mixing matrices)
2. Spectral action computation on the full SM triple (Higgs potential, mass relations)
3. n_max=3 verification (dim_H = 1280; computationally heavier but no structural obstacle)
4. Paper 32 §VIII.C update with G4a results

## §7. Files

**Created:**
- `geovac/standard_model_triple.py` (~530 lines)
- `tests/test_standard_model_triple.py` (45 tests, all passing)

**Not modified:** No existing files changed. G4a is purely additive.
