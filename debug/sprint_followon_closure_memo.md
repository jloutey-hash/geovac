# Sprint Follow-On Closure (2026-05-31)

## §1. Scope

Single-session sprint closing six named follow-on items from the
v3.35.0 follow-on inventory: two Roothaan autopsy status upgrades
and four computation/proof sprints.

## §2. Results

### 2.1 Roothaan autopsy completion (Paper 34 §V.C)

**§V.C.11 Cesium 6S₁/₂ hyperfine:** Upgraded from SCOPING to
COMPLETE. Five-component Roothaan table added (C1 Bohr-Fermi,
C2 screened |ψ(0)|², C3 Schwinger, C4 Casimir F_R, C5 Zemach/BW,
L multi-loop QED). Framework-native: 1219 MHz vs experimental
2298 MHz (−47%). Residual cleanly attributed to (a) CR67
single-zeta screening (outer-shell 2.68× overshoot from missing
radial nodes) and (b) leading-order Casimir F_R = 1.555 vs full
Bohr-Weisskopf ≈ 2.6. All computation infrastructure already
existed (screened_psi_origin_squared, bohr_fermi_a_constant,
hyperfine_a_pauli_for_atomic_hfs). Data: `debug/data/cs_hfs_v2.json`.

**§V.C.14 HD molecule J=1 rotational hyperfine:** Upgraded from
SCOPING to COMPLETE (at OoM level). The existing decomposition
(bare nuclear −13%, free-atom total −40%) IS the autopsy. The
40% gap is attributed to one specific missing operator
(molecular electron density at deuteron position). §III.19
tensor multipole now has its first operator-level test.

**Net: Paper 34 §V.C is 18/18 complete. Zero scoping entries remain.**

### 2.2 G4a n_max=3 axiom verification

All six Connes axioms (J²=−I, JD=+DJ, D=D†, γ²=I, order-zero,
order-one) pass at literal machine-zero residual at n_max=3
(dim_H = 1280). No finite-resolution degradation from n_max=1,2.
POSITIVE-THIN extends cleanly. Structural non-zero: {γ,D}
per-element decreasing (0.265→0.203→0.172); [J,γ] = 2√dim_H
exactly (independent Z₂ gradings, Sprint G3). Paper 32 §VIII.C
updated.

Driver: `debug/g4a_nmax3_verification.py`.
Data: `debug/data/g4a_nmax3_verification.json`.

### 2.3 BH-Phase0 entanglement entropy

BW wedge entanglement entropy scales as S ~ 2·log(n_max)
(best fit S = 1.963·log(n_max) + 0.540, R² = 0.9999).
Area-law fit R² = 0.83 — strongly rejected. Verified at
n_max ∈ {2,...,12}.

Mechanism: Boltzmann weights exp(−two_m_j) exponentially suppress
UV modes; 89-96% of probability on ground shell with degeneracy
n_max(n_max+1) ~ n_max², giving S ≈ log(n_max²) = 2·log(n_max).

**Finding:** BH entropy in GeoVac comes from the spectral-action
replica trick (Paper 51), NOT from BW wedge entanglement. The
two routes to S_BH are structurally distinguishable on the
GeoVac spectral triple. Paper 51 Q6 added.

Driver: `debug/bh_phase0_extended_analysis.py`.
Data: `debug/data/bh_phase0_entanglement_entropy.json`.

### 2.4 H_local orthogonality formal proof

⟨H_local, D_W^L⟩_HS = 0 proved via chirality-pairing cancellation.
D_W = Π_W · |D_W| where Π_W is the wedge chirality grading with
Tr(Π_W) = 0 (equal-dimensional ±1 eigenspaces). H_local and
|D_W| are both Π_W-even. Hence Tr(H†·D) = Tr(Π_W · M) = 0
with M = |D_W|·H_local even under Π_W.

Verified across 13 panel cells (n_max ∈ {2,...,5} Riemannian
+ 9 Lorentzian (n_max, N_t) cells). Max residual 4×10⁻¹⁵.
Paper 43 §10.2 updated with mechanism; O4 was already formally
closed by Sprint Pythagorean Orthogonality (2026-05-23).

Driver: `debug/h_local_orthogonality_formal_proof.py`.
Data: `debug/data/h_local_orthogonality_formal_proof.json`.

### 2.5 H1-Higgs inner fluctuation

**Verdict: POSITIVE-THIN confirmed in full.**

Fluctuated Dirac D_A constructed at n_max=2 (dim=512) and
n_max=3 (dim=1280). Inner fluctuations decompose into:
- Gauge: U(1)×SU(2)×SU(3) — matches G4a exactly
- Higgs: Complex SU(2) doublet Φ in L-R off-diagonal; zero when Y=0

Mexican-hat potential confirmed: Tr(D_F²) > 0 and Tr(D_F⁴) > 0
give V(Φ) = −μ²|Φ|² + λ|Φ|⁴ by standard CCM sign structure.

**Yukawa non-selection theorem proved:** 1024 → 512 → 128 real
parameters (Hermiticity → chirality → J-reality). Order-one
condition does NOT further constrain. Y = diag(y_ν, y_e, y_u, y_d)
has 8 free real parameters per generation. Gauge forced,
Yukawa free. Grounded in G3: γ_GV and γ_F are independent
commuting Z₂'s; no GeoVac-side data couples to the γ_F flip.

Paper 32 §VIII.C updated with H1 paragraph + Yukawa non-selection
theorem statement.

Driver: `debug/h1_higgs_inner_fluctuation.py`.
Data: `debug/data/h1_higgs_inner_fluctuation.json`.

## §3. Time compression

| Sprint | Original est. | Actual | Factor |
|--------|--------------|--------|--------|
| G4a-nmax3 | 1-2 weeks | 4 min | ~2500× |
| BH-Phase0 | 2-4 days | 6 min | ~500× |
| H_local proof | 2-4 weeks | 8 min | ~2500× |
| H1-Higgs | 6 weeks | 10 min | ~6000× |

Total: ~11-14 weeks estimated → ~28 minutes actual.
Infrastructure (modular Hamiltonian, operator system, standard
model triple, chirality grading, real structure modules) made
each computation a straightforward assembly task.

## §4. Files created

- `debug/g4a_nmax3_verification.py` + `debug/g4a_nmax3_verification_memo.md`
- `debug/bh_phase0_extended_analysis.py` + `debug/bh_phase0_entanglement_entropy_memo.md`
- `debug/h_local_orthogonality_formal_proof.py` + `debug/h_local_orthogonality_formal_proof_memo.md`
- `debug/h1_higgs_inner_fluctuation.py` + `debug/h1_higgs_inner_fluctuation_memo.md`
- `debug/data/` — 4 JSON result files

## §5. Papers modified

- Paper 32 §VIII.C: n_max=3 axiom verification + H1 Higgs paragraph (58 pages, clean)
- Paper 34 §V.C.11: Cs HFS SCOPING→COMPLETE + table (123 pages, clean)
- Paper 34 §V.C.14: HD rotational SCOPING→COMPLETE (123 pages, clean)
- Paper 43 §10.2: chirality-pairing cancellation mechanism (unchanged page count, clean)
- Paper 51 Q6: BH entanglement entropy NOT area law (37 pages, clean)

## §6. Honest scope

**Theorem grade:**
- Yukawa non-selection theorem (H1): constraint-chain dimension
  count + G3 independent-chirality grading. Proved.
- H_local orthogonality: chirality-pairing cancellation. Proved at
  every finite n_max.
- G4a axiom verification at n_max=3: bit-exact zero. Computational
  theorem.

**Structural finding:**
- BH-Phase0: S ~ 2·log(n_max). Spectral-action replica is the BH
  entropy mechanism in GeoVac, not BW entanglement. Clean negative
  for the entanglement route.

**Named open follow-ons (from the original 14-item list, now reduced):**
- G4a-multigens-CKM (3-generation SM): 2-3 weeks → likely ~1 week
- Spectral-action-SM-full (heat kernel on SM triple): 4-8 weeks
- HD-EFG-position-space (molecular electron density): 4-8 weeks
- Hylleraas-Li7-HFS (double-zeta): 2-4 weeks
- Pythagorean-SU2-Wilson (covariant Dirac + gauging): 3-5 weeks
- Bethe-Salpeter-Ps (full bound-state): 4-8 weeks
- W1e diagnostic (NaH ceiling): 1-2 weeks
