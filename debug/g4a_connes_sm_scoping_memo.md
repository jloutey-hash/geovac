# Sprint G4a Scoping: Connes SM on T_{S³} alone, A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ)

**Sprint:** G4a (formal Connes Standard Model unification on the GeoVac S³ spectral triple, single-manifold flavor)
**Date:** 2026-05-06
**Author:** PM scoping fork
**Status:** Scoping only. No implementation. Recommendation at §7.
**Predecessors read:** `geovac/almost_commutative.py` (854 lines, H1 implementation), `debug/h1_ac_extension_memo.md` (H1 verdict POSITIVE-THIN), `debug/g3a_chirality_memo.md` (γ_GV = σ_x ⊗ I), `debug/g3b_chirality_F_audit_memo.md` (γ_F = diag(+1,+1,−1,−1) ⊕ −diag(...) at KO-dim 6), `debug/g4_cross_manifold_scoping_memo.md` (G4 split into G4a + G4b).

---

## §0. Executive verdict

**Reachable in 1.5–2 months. Predicted POSITIVE-THIN (same flavor as Sprint H1).**

The extension `A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ)` is the canonical Connes SM finite algebra. On top of T_{S³} alone (no S⁵ factor — see §3), the inner-fluctuation calculation produces U(1) × SU(2) × SU(3) gauge content from inner automorphisms of A_F, plus a richer Higgs sector with off-diagonal blocks across all three summands. The construction is reachable as a 1.5–2 month sprint extending `geovac/almost_commutative.py`. The expected verdict is positive-thin in the same sense H1 landed: the construction works, it produces the SM gauge group, **but no GeoVac mechanism autonomously selects the SM-distinguishing data** (Yukawas, hypercharges, generation count, neutrino mixing, the M_3(ℂ)-ℍ off-diagonal that gives quark masses).

The G3 negative (γ_GV / γ_F independence on H_GV ⊗ H_F, structurally one ℤ₂ short of weak-isospin co-location) **does not block G4a**, but it **does cap the verdict thinness**: without chirality co-location, the SU(2)_L vs SU(2)_R selection is a free choice on the AC side, not a GeoVac-derived asymmetry. That sharpens what G4a closes (the *gauge-group structure*) and what it cannot close (the *electroweak chirality assignment*).

This is a worthwhile sprint **if and only if** the right deliverable for the project is "GeoVac sits inside the Marcolli–vS / Connes-SM lineage as a single-manifold spectral triple producing the full SM gauge group via inner fluctuations, with all electroweak data identifiable as Paper-18 calibration constants." That deliverable is publishable, and it would be the most direct extension of WH1 PROVEN; but it is not a derivation of the SM. Recommendation at §7 is conditional: **open after Paper 38 lands and after PI confirms the deliverable framing.**

---

## §1. Architecture: J_AC, gradings, and Hilbert space doubling for M_3(ℂ)

### §1.1 Algebra and Hilbert space

**Algebra.**

  A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ)

with the standard Connes SM matter-sector representation acting on a 32-dim matter Hilbert space (Connes–Marcolli 2008 Ch. 13.4–13.7). One generation, in (ν_L, e_L, ν_R, e_R, u_L^a, d_L^a, u_R^a, d_R^a) basis with a = 1, 2, 3 the color index:

  H_F^matter = ℂ⁴ ⊕ ℂ²⊗ℂ³ ⊕ ℂ²⊗ℂ³  (lepton + quark-L + quark-R)
              = ℂ⁴ ⊕ ℂ⁶ ⊕ ℂ⁶
              = ℂ¹⁶

(The convention dimension count: 2 leptons L + 2 leptons R + 2 quarks L × 3 colors + 2 quarks R × 3 colors = 4 + 4 + 6 + 6 = 20. Connes–Marcolli normalize differently — they compute dim 16 for one-generation matter via L-R doublings inside M_2(ℂ) algebra — but the right number for our purposes is the matter sector for one generation, which is the dim_C(H_F^matter) = 16 entry on Connes–Marcolli Table 13.6.)

**Action of A_F on H_F^matter.** With (λ, q, m) ∈ ℂ ⊕ ℍ ⊕ M_3(ℂ):

  π_F(λ, q, m) = block_diag(
      π_lepton(λ, q),                             # 4×4 matter-leptons
      q ⊗ I_3,                                     # 6×6 quarks-L (ℍ acts on flavor, M_3(ℂ) acts on color via 1; 1 here is the *embedding* of color-trivial)
      diag(λ I_3, conj(λ) I_3 + something with m)  # 6×6 quarks-R (ℂ acts on (u_R, d_R) flavors as in lepton sector; M_3(ℂ) acts on color via *standard color triplet*)
  )

The exact form of the M_3(ℂ) embedding distinguishes Connes' SM ("particle physics" reference; Connes–Chamseddine 1996, Connes 1996; expanded in Connes–Marcolli 2008 Theorem 13.5) from naive "U(1) × SU(2) × SU(3) on H_F^matter" representations.

**Doubled Hilbert space.** As in Sprint H1:

  H_F = H_F^matter ⊕ H_F^antimatter = ℂ¹⁶ ⊕ ℂ¹⁶ = ℂ³²

A_F acts on the matter sector only. The antimatter sector picks up the opposite-algebra action via J_F.

### §1.2 D_F: matter-sector Dirac with full Yukawa structure

D_F = block_diag(M, M̄) with M the 16×16 matter-sector Dirac:

  M = [   0           Y_lepton^†      0           0       ]
      [   Y_lepton    0               0           0       ]   # leptons: L↔R
      [   0           0               0           Y_quark^†]
      [   0           0               Y_quark     0        ]   # quarks: L↔R

with Y_lepton = diag(y_ν, y_e) ∈ M_2(ℂ) and Y_quark a 6×6 Yukawa matrix mixing the (u, d) flavor doublet across (R, L), with structure

  Y_quark = U_CKM^† · diag(y_u I_3, y_d I_3) · U_CKM

where U_CKM is the CKM matrix on the L-quark side (color-trivial, acts on flavor). For a single generation U_CKM = I and Y_quark just gives diagonal up-type and down-type quark masses in the color-triplet sector.

The minimal Yukawa data for Connes' SM at one generation is **{y_ν, y_e, y_u, y_d}** — four real masses (or four complex Yukawa entries). With three generations, the data extends to **diagonal masses + CKM mixing matrix + PMNS mixing matrix + Majorana neutrino sector if applicable**. Sprint G4a should restrict to ONE GENERATION at the start; multi-generation is a downstream extension.

### §1.3 J_F (real structure on F factor) and KO-dim arithmetic

**J_F definition** (extends H1):

  J_F = U_F · K   with U_F = σ_x ⊗ I_16

where K is complex conjugation. As in H1, J_F² = +I (KO-dim 6 of T_F). And J_F D_F = +D_F J_F because D_F is block_diag(M, M̄) with the matter-block conjugated to the antimatter-block exactly under K.

**γ_F (chirality grading on F factor)** — adapting G3-B's KO-6 convention to the SM-extended A_F:

  γ_F^matter = diag(+I_2, −I_2, +I_2 ⊗ I_3, −I_2 ⊗ I_3)
             = diag(+1, +1, −1, −1, +1, +1, +1, +1, +1, +1, −1, −1, −1, −1, −1, −1)

  γ_F^antimatter = − γ_F^matter

  γ_F = block_diag(γ_F^matter, γ_F^antimatter)  (32×32)

This places L-doublets (lepton-L, quark-L) at γ_F = +1 and R-singlets (lepton-R, quark-R) at γ_F = −1, with antimatter sign-flipped to satisfy {J_F, γ_F} = 0 (KO-6 sign rule).

**γ_GV** stays as in G3-A: γ_GV = σ_x ⊗ I on the chi-doubled full Dirac sector.

**Combined real structure** J_AC = J_GV ⊗ J_F:
  - J_AC² = (J_GV²)(J_F²) = (−1)(+1) = −1 ⟹ J_AC² = −I.
  - J_AC D_AC = +D_AC J_AC (from H1 §1.2 generalized).
  - Combined KO-dim = 3 + 6 = 9 ≡ 1 (mod 8). ε = (−)(+) = −, ε' = (+)(+) = +. Same as H1.

**Combined chirality** γ_AC = γ_GV ⊗ γ_F:
  - {γ_AC, D_AC} structure: D_AC = D_GV ⊗ I_F + γ_GV ⊗ D_F. Compute
    {γ_AC, D_AC} = {γ_GV ⊗ γ_F, D_GV ⊗ I + γ_GV ⊗ D_F}
                  = {γ_GV, D_GV} ⊗ γ_F + γ_GV² ⊗ {γ_F, D_F}
                  = 0 ⊗ γ_F + I ⊗ {γ_F, D_F}.
    Need {γ_F, D_F} = 0. Since γ_F is L−R sign and D_F is L↔R off-diagonal, {γ_F, D_F} = 0 holds identically. ⟹ {γ_AC, D_AC} = 0. ✓
  - γ_AC² = γ_GV² ⊗ γ_F² = I ⊗ I = I. ✓
  - [γ_AC, π(a ⊗ b)] = 0 for any algebra element. ✓

So γ_AC is a well-defined chirality grading on the combined real spectral triple. The architecture is structurally clean.

### §1.4 What makes G4a different from H1

| Aspect | Sprint H1 | Sprint G4a |
|--------|-----------|------------|
| Algebra A_F | ℂ ⊕ ℍ | ℂ ⊕ ℍ ⊕ M_3(ℂ) |
| dim_C(H_F^matter) | 4 | 16 |
| dim_C(H_F) (doubled) | 8 | 32 |
| Yukawa parameters | 2 (y_ν, y_e) | 4 minimum (y_ν, y_e, y_u, y_d) per generation |
| Gauge groups recovered | U(1) × SU(2) | U(1) × SU(2) × SU(3) |
| Higgs sector | 1 block (ℂ ↔ ℍ) | 3 blocks (ℂ ↔ ℍ, ℂ ↔ M_3(ℂ), ℍ ↔ M_3(ℂ)) |
| KO-dim arithmetic | 3 + 6 = 9 ≡ 1 (mod 8) | unchanged |

The KO-dim arithmetic and J/γ axioms are identical to H1; the additions are *dimensional* (more Hilbert space), *algebraic* (a third summand of A_F), and *parametric* (more Yukawa entries).

---

## §2. Ingredients from Sprint H1: what transfers, what doesn't

### §2.1 Transfers cleanly (re-use Sprint H1 machinery directly)

| H1 Ingredient | G4a Re-use | Notes |
|---------------|-----------|-------|
| `geovac/almost_commutative.py` `AlmostCommutativeTriple` skeleton | Yes, extend with `M_3(ℂ)` summand | The `algebra_action`, `dirac_F`, `real_structure_F` methods need M_3(ℂ) extensions |
| GV factor: T_{S³} via `FullDiracTruncatedOperatorSystem` | Yes, unchanged | n_max ∈ {2, 3} sweep stays |
| Truthful CH Dirac D_GV | Yes, unchanged | Track 2 verdict locks this |
| γ_GV = σ_x ⊗ I (G3-A convention) | Yes, unchanged | G3-A delivers it ready-to-use |
| J_GV (real structure on H_GV) | Yes, unchanged | Already verified at KO-3 in `geovac/real_structure.py` |
| Inner fluctuation formula ω = Σ a_i [D, b_i] + ε' J ω J⁻¹ | Yes, unchanged | The decomposition into gauge + Higgs pieces is the same; the Higgs piece just has more blocks |
| Spectral action coherence test (Tr(D²) polynomial in Y) | Yes, generalize | Tr(D²) = Tr(D_0²) + 4·dim_GV·(|y_ν|² + |y_e|² + 3|y_u|² + 3|y_d|²) — the factor of 3 is the color count |
| Falsifier sweep at n_max ∈ {1, 2, 3} | Yes, extend Yukawa parameter space | Higgs ‖Φ‖_max should be > 0 for any non-zero Y; falsifier holds iff Y = 0 |
| Matter/antimatter localization (off-block = 0) | Yes, by construction | The 32×32 Hilbert space stays (matter, antimatter) doubled |
| Pauli sigma machinery (SIGMA_0..3, quaternion_to_matrix) | Yes, extends with Gell-Mann matrices | Need the 8 Gell-Mann generators of SU(3) ⊂ M_3(ℂ) for Cartan / SU(3)-traceless basis |

### §2.2 Requires new construction

| Component | Reason | Effort estimate |
|-----------|--------|-----------------|
| M_3(ℂ) action embedding (color-triplet rep on quark sectors) | Connes' SM uses a specific embedding that is NOT just M_3(ℂ) ⊗ I; it interleaves quark-L with the SU(2) doublet structure (Connes–Marcolli §13.5). Need careful rep theory | 1–2 weeks |
| Quark sector in H_F^matter (6+6=12 quark dim) | Currently H1 has only 4 lepton dim | 3 days |
| Three Higgs blocks in inner fluctuation decomposition | H1 has 1 Higgs block (ℂ↔ℍ); G4a has ℂ↔ℍ, ℂ↔M_3(ℂ), ℍ↔M_3(ℂ) | 1 week |
| Quark Yukawa Y_quark with CKM (one-gen: just diag(y_u, y_d) ⊗ I_3) | Yes, similar to lepton Y_lepton | 3 days |
| Order-zero / order-one verification at n_max ∈ {1, 2, 3} | The matter/antimatter doubling makes order-zero trivial; order-one needs to be checked at the new dim | 1 week |
| Test suite extension `tests/test_almost_commutative.py` | 53 → ~120 tests projected | 2 weeks |

### §2.3 What H1 ruled out that re-emerges (and is again ruled out)

The H1 candidates A and B that the scoping memo §3.3 named:
- **Candidate A** (offdiag CH bridging) — RULED OUT by Track 2's clean negative on offdiag CH J-symmetry. *Same status for G4a.* The offdiag CH SDP-bounding Dirac is still the wrong base for any spectral-action / inner-fluctuation work; G4a must use truthful CH.
- **Candidate B** (Sturmian / DUCC bridging) — not pursued in H1 and not pursued in G4a. The same ambiguity (whether inter-shell Sturmian couplings are a Yukawa) reappears, but solving it is not the G4a deliverable. Yukawa entries stay free inputs.

So the H1 free-Yukawa verdict carries directly: G4a will ALSO have free Yukawa entries (now 4+ instead of 2), and the same positive-thin verdict will hold.

---

## §3. The G3 negative (γ_GV / γ_F independence): does it constrain G4a?

### §3.1 What G3 closed

G3-C (per the closing memos, `debug/g3c_tensor_chirality_memo.md`) verified that γ_GV ⊗ I_F and I_GV ⊗ γ_F are **independent commuting ℤ₂ operators** on H_GV ⊗ H_F. At every n_max ∈ {1, 2, 3}:

  ‖γ_GV ⊗ I_F − I_GV ⊗ γ_F‖_op = 2 exactly.

Their common eigenspace decomposition is 1:2:1 multiplicities of {−2, 0, +2}, σ_x/σ_z convention-independent. In short: the two ℤ₂ gradings are independent, not co-located.

### §3.2 Why G3 does not block G4a

The G3 negative says GeoVac does not autonomously identify the GV-side chirality (which signs the Camporesi–Higuchi Dirac eigenvalue) with the AC-side weak-isospin chirality (which signs L vs R in the SU(2) ⊂ ℍ-matter sector). This is a structural finding about what GeoVac data does and doesn't determine — but it does NOT prevent the algebraic construction of the larger triple T_{S³} ⊗ T_F^{SM}.

G4a's deliverable is "the construction works, the gauge group is recovered, and the data not derived from GeoVac is on a finite list of Paper-18 calibration constants." G3 just tells us **which items go on that list**: γ_F as a free choice (which selects which sector of H_F is "L" and which is "R") is one of those free choices. It is not "more thin than expected" — H1 already had this issue with Y as a free input; G3 sharpens to "γ_F itself is a free choice."

### §3.3 What G3 *does* cap

G3's negative does cap how strong G4a's verdict can be:

| Without G3 (counterfactual) | With G3 (actual) |
|-----------------------------|-------------------|
| G4a closes GeoVac → SM gauge group up to chirality assignment, where GV chirality = weak-isospin | G4a closes GeoVac → SM gauge group, but weak-isospin chirality is still free input |
| Marcolli–vS-without-Higgs side, plus structural chirality co-location | Marcolli–vS-without-Higgs side, plus chirality independence |

This sharpens the §VIII.B G2/G3/G4 entries in Paper 32 simultaneously: G2 (no autonomous Higgs), G3 (no autonomous chirality co-location), G4a (no autonomous SU(3) — i.e., M_3(ℂ) is a free choice of A_F, not derived from GeoVac structure) are **three names for one structural fact**: GeoVac does not autonomously determine the SM-distinguishing finite data. The construction is consistent with the SM, but it is also consistent with any other choice of A_F, γ_F, and Y.

This is honest and useful: it tells us *exactly* where the SM lives in the GeoVac framework — as a Paper-18 calibration constant package, not as a derivation.

### §3.4 Could G4a close G3 retroactively?

Not by itself. The G3 negative is about whether γ_GV (a GeoVac-side object) identifies with γ_F (an AC-side object). G4a *adds* algebra-side structure (M_3(ℂ)) but doesn't change the GeoVac side. So G3 stays open / negative under G4a.

A *retroactive* G3 closure would require some kind of GeoVac-derived weak-isospin chirality on H_F — i.e., a derivation that the matter sector's L/R split *must* match the GV-side chirality. The G3-C ‖γ_GV ⊗ I_F − I_GV ⊗ γ_F‖_op = 2 result rules this out as an algebraic consequence of the existing axioms. So G4a leaves G3 in its sealed-negative state; G4a's verdict has G3 as a noted caveat.

---

## §4. Predicted positive-thin verdict: what specifically is "thin"

### §4.1 The thin items

G4a's positive-thin verdict places the following items on the **free choices / calibration constants** list:

1. **Choice of A_F** (ℂ vs ℂ ⊕ ℍ vs ℂ ⊕ ℍ ⊕ M_3(ℂ) vs Pati–Salam). GeoVac does not pick. — **major thin**
2. **Number of generations** (one vs three). The H_F dim scales linearly with generations. GeoVac does not pick. — **major thin**
3. **Yukawa entries y_ν, y_e, y_u, y_d** (per generation). GeoVac does not pick. — **major thin** (same as H1)
4. **CKM and PMNS mixing matrices** (multi-generation). GeoVac does not pick. — **major thin**
5. **Hypercharge assignments** within H_F^matter. The H_F internal structure (which states are R-singlets, which are L-doublets) is an external input. — **medium thin**
6. **γ_F convention** (sign of L vs R). G3-C established this is independent of γ_GV. — **medium thin**
7. **Dirac vs Majorana neutrino sector**. GeoVac does not pick. — **medium thin**
8. **Choice of M_3(ℂ) embedding into matter-action** (color-triplet rep details). Standard SM picks one; GeoVac does not autonomously discriminate. — **minor thin**

The list is 4 major + 3 medium + 1 minor = a substantial number of free choices. This is what "positive-thin" means in concrete: the construction works and produces the SM gauge group, but every SM-distinguishing data point is a free input.

### §4.2 The positive items

G4a closes:

1. **The construction itself** is well-defined at finite n_max on truthful CH GeoVac. (This is the H1 closure scaled up.) — **major positive**
2. **U(1) × SU(2) × SU(3) gauge content** is recovered as inner-automorphism gauge fluctuations from A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ). (This is the *structural* statement; gauge group is determined by A_F's center, and inner fluctuations of the GeoVac D produce the corresponding Lie algebras.) — **major positive**
3. **Three Higgs blocks** (ℂ↔ℍ, ℂ↔M_3(ℂ), ℍ↔M_3(ℂ)) appear as off-diagonal Yukawa data; non-zero for any non-zero Y; absent for Y = 0. (The H1 verdict generalizes to three blocks.) — **major positive**
4. **Spectral action coherence**: Tr(D²) is polynomial in Y of degree 2, with predictable color factor 3 multiplying quark Yukawa contributions. — **major positive**
5. **Order-zero and order-one** are automatic (matter/antimatter doubling). — **major positive**
6. **The combined real structure J_AC and chirality γ_AC** are well-defined at KO-dim 9 ≡ 1 (mod 8) with verified J²=−I, JD=+DJ, {γ, D}=0, [γ, a]=0 axioms. — **medium positive** (extends H1 §1.2)
7. **The Marcolli–vS-without-Higgs reading at SM scale**: GeoVac is consistently in this lineage at the full SM gauge content level, not just electroweak. — **major positive** (this is the G4a deliverable's structural punch)

### §4.3 What G4a does NOT close

In addition to the §4.1 thin items:

- **No native derivation of mass spectrum**: even if all Yukawas were fixed, the actual mass scales (m_e ≈ 0.5 MeV, m_t ≈ 173 GeV, etc.) are dimensional inputs that come from the spectral-action heat-kernel cutoff Λ — itself a calibration choice.
- **No native derivation of α(M_Z)**: the Paper 2 conjectural reading remains paused. G4a doesn't change anything for the α-derivation thread (CLAUDE.md WH5).
- **No structural reason for SM gauge content vs other consistent choices**: if a researcher argued for Pati–Salam (A_F = M_2(ℂ) ⊕ M_2(ℂ) ⊕ M_4(ℂ), say) instead, GeoVac would equally accept it. The SM is not selected on GeoVac-internal grounds.
- **G4b (cross-manifold)** is unchanged. Still a structural NCG framework gap.

### §4.4 The clean positive-thin headline

> *Sprint G4a closes the formal Connes SM unification on the GeoVac S³ spectral triple at the gauge-group level: U(1) × SU(2) × SU(3) all arise as inner fluctuations on the same underlying triple, with the same KO arithmetic and J/γ axioms as Sprint H1. The full Connes SM construction is consistent with GeoVac. However, GeoVac data does not autonomously select any of the SM-distinguishing finite data (algebra choice, generation count, Yukawas, mixings, hypercharges, neutrino sector, chirality assignment). The verdict is "GeoVac sits inside the Marcolli–vS / Connes-SM lineage; it is not a derivation of the SM."*

---

## §5. The explicit falsifier

Following Sprint H1's strong-falsifier protocol, the analogous statement for G4a:

> **Strong falsifier (G4a):** Show that for *every* Hermitian D_F on A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ) derivable from GeoVac structure alone — i.e., not specified by hand — the resulting inner fluctuation produces *only* gauge-1-forms with zero Higgs sector across all three blocks.

**Predicted result:** Falsifier holds iff D_F = 0 (zero Yukawas across all sectors). For any non-zero imposed Yukawa data, at least one Higgs block (ℂ↔ℍ for lepton, ℂ↔M_3(ℂ) for u-quark, ℍ↔M_3(ℂ) for d-quark) is non-trivial.

This *would* close G4a as a positive structural result (Higgs from inner fluctuation works at the full SM scale) without producing a derivation of the Yukawas. Same flavor as H1.

**Sweep design:** at n_max ∈ {2, 3} on the GeoVac side, scan Yukawa parameters (y_ν, y_e, y_u, y_d) across {0, 0.1, 0.3, 1.0} each. Verify ‖Φ_lepton‖, ‖Φ_u‖, ‖Φ_d‖ are zero iff the corresponding Y entry is zero. Verify Tr(D²) is polynomial in Y as predicted. Also stress-test order-zero and order-one numerically at the larger H_F dim (32×dim_H_GV).

A *secondary falsifier* worth running:

> **Secondary falsifier (auxiliary):** Show that for any choice of A_F that is *not* ℂ ⊕ ℍ ⊕ M_3(ℂ) (e.g., Pati–Salam, or a Higgs-less reduced Lie algebra), the GeoVac side produces equally well-defined inner-fluctuation gauge content. Confirms GeoVac autonomy at zero — i.e., GeoVac doesn't pick the SM over alternatives.

This secondary falsifier is the structural sharpening of the G4a verdict: the SM is one *choice* of A_F that GeoVac is consistent with, not the unique selection.

---

## §6. Sprint scope estimate

### §6.1 Modules to add / extend

| File | Action | Lines |
|------|--------|-------|
| `geovac/almost_commutative.py` | Extend `ElectroweakFiniteTriple` → `StandardModelFiniteTriple` (or new class); add color-triplet matter-sector representation, quark Yukawas, three-block Higgs decomposition | +400–600 |
| `geovac/standard_model_action.py` (new) | Optionally factor out the color-action machinery (Gell-Mann matrices, M_3(ℂ) embedding) for reuse | +200 |
| `geovac/inner_fluctuation_decomposition.py` (new, optional) | Extract the three-Higgs-block decomposition; cleaner factoring than packing into AlmostCommutativeTriple | +150 |
| `geovac/real_structure.py` | Add `extend_J_F_to_SM` helper if the M_3(ℂ) factor is non-trivial under K | +50 |
| `geovac/full_dirac_operator_system.py` | (No change expected; T_{S³} side is unchanged.) | 0 |

**Net new code: ~600–1000 lines.**

### §6.2 Tests to add

| Test file | Action | Tests |
|-----------|--------|-------|
| `tests/test_almost_commutative.py` | Extend with G4a section; add tests for M_3(ℂ) action, quark sector, three-Higgs decomposition, KO-9 axioms at extended dim | +60–80 |
| `tests/test_standard_model_action.py` (new) | Color-triplet rep tests, Gell-Mann basis closure, gauge group center identification | +20–30 |
| `tests/test_inner_fluctuation_decomposition.py` (new) | Three-Higgs block tests, gauge piece independence, falsifier sweep | +15–20 |

**Net new tests: ~95–130. Going from 53 to ~150–180 in `test_almost_commutative.py` total, a ~3× test count expansion.**

### §6.3 Driver / data files

| Driver | Action |
|--------|--------|
| `debug/g4a_sprint_analysis.py` | Falsifier sweep at n_max ∈ {2, 3} × Yukawa parameter grid |
| `debug/data/g4a_falsifier.json` | Output: ‖Φ‖ entries for each block × parameter combination |
| `debug/g4a_connes_sm_memo.md` | Closure memo (analog of `h1_ac_extension_memo.md`) |
| `debug/g4a_secondary_falsifier.py` | Run alternative-A_F sweep (Pati–Salam, Higgs-less, etc.) |

### §6.4 Time budget

| Phase | Duration |
|-------|----------|
| 1. M_3(ℂ) embedding rep theory + tests | 1.5 weeks |
| 2. Quark sector in H_F + Y_quark + tests | 1 week |
| 3. Three-Higgs decomposition implementation + tests | 1 week |
| 4. Combined J_AC, γ_AC at KO-9 verification at extended dim | 1 week |
| 5. Falsifier sweep + secondary falsifier + analysis | 2 weeks |
| 6. Memo writing + Paper 32 §VIII.D addendum draft | 1 week |
| **Total** | **6–7 weeks (≈1.5 months)** |

This is on the optimistic end of the G4 scoping memo §4.1's "1.5 months" estimate — that estimate had 3–6 weeks for algebra + 2–4 for falsifier, totaling 1.5 months. G4a is reachable in this window.

### §6.5 Risk factors

- **Color-triplet rep theory complexity** (Phase 1): Connes' embedding M_3(ℂ) → π_F is non-trivial; need to consult Connes–Marcolli §13.5 carefully. If the embedding is simpler than expected (just π_F(m) = m on the quark sector with trivial flavor), Phase 1 contracts to 0.5 weeks. If more subtle (twisting by hypercharge), could take 2.5 weeks. **Risk: ±1 week.**
- **Quark Yukawa multi-generation extension**: deferred to G4a-followup, but if PI directs one-generation only, this stays simple. **Risk: low if scoped to one-generation.**
- **Numerical conditioning at H_F dim 32**: at n_max=3, dim_H_GV ≈ 40 (full Dirac sector) × 32 = 1280, so 1280×1280 matrices. Should be tractable (H1 already runs at ≈320×320). Test sweep dimensions = (Yukawa grid size)^4 × 3 n_max points; tractable on a workstation. **Risk: low.**
- **Order-one verification at extended dim**: may surface finite-resolution failures analogous to H1's at 5–20% (which were noted as artifacts of the Connes–vS truncated operator system structure). If new failures appear at the SM scale that don't have an obvious analog, may need diagnostic work. **Risk: medium, would extend Phase 4 by 0.5–1 week.**
- **Paper 32 §VIII.D addendum scope**: PM should not pre-commit; PI decides whether the result merits §VIII.D, a Paper 39, or a §VIII.C extension. **Risk: low (just framing, after the implementation lands).**

---

## §7. Recommendation

### §7.1 The verdict

**Open the sprint, but conditionally.**

Specifically:

1. **Wait for Paper 38 to be drafted, reviewed, and submitted.** Paper 38 is the J. Geom. Phys. companion paper for WH1 PROVEN, the cleanest deliverable in the project's history. Opening G4a before Paper 38 lands divides attention; finishing Paper 38 first preserves the WH1 PROVEN claim as a focused submission.
2. **Confirm with PI that the deliverable framing is "GeoVac in the Marcolli–vS-without-Higgs lineage, full SM gauge content."** If PI is interested in a different framing (e.g., trying to derive Yukawas from GeoVac, or exploring G4b cross-manifold instead), G4a is not the right sprint. The framing matters because it determines whether the predicted positive-thin verdict is the *desired* outcome or a *disappointment*.
3. **Restrict to one generation initially.** Multi-generation extension is a follow-up sprint; opening G4a as a 3-generation construction expands scope by 30–50% and risks sprint-bloat.
4. **Do not write a standalone paper for G4a.** Per CLAUDE.md §1.5, the result is too thin to justify a self-contained paper. The natural venue is Paper 32 §VIII.D addendum extending §VIII.C (the H1 addendum), which keeps the spectral-triple synthesis paper as the home for all SM-related GeoVac structural claims.

### §7.2 Why open it eventually

G4a is the most direct extension of WH1 PROVEN to the SM scale. The deliverable — "GeoVac sits inside the Connes-SM lineage; the construction works at finite n_max; the full SM gauge group is recovered as inner fluctuations on a single S³ spectral triple; SM-distinguishing data is on a finite list of Paper-18 calibration constants" — is the strongest *honest* statement GeoVac can make about its connection to the Standard Model. It locks in the project's position vis-à-vis Connes' SM program: GeoVac is a single-manifold geometric realization of the same spectral-triple machinery, with the GeoVac packing axiom replacing the "M^4 × F" continuum-spacetime ansatz.

This is publishable in a Connes-program-aware venue (J. Geom. Phys., Adv. Theor. Math. Phys., or a CMP follow-up to Paper 38), and would close the §VIII.B SM gauge-content question of Paper 32 in the same way Sprint H1 closed the electroweak slice. Without G4a, the SM-gauge claim in Paper 32 §VIII.B remains forward-pointing rather than realized.

### §7.3 Why not open it now

Three reasons against opening G4a immediately:

1. **Paper 38 is the higher-priority deliverable.** WH1 PROVEN is the centerpiece of the May 2026 work; getting Paper 38 to arXiv and J. Geom. Phys. submission keeps the framing tight. Opening G4a in parallel risks the PI / PM split-attention failure mode.
2. **The deliverable is a *negative on autonomous SM emergence* dressed as a structural positive.** That's a valuable result, but it's a delicate framing exercise. Better to land Paper 38 first and let the WH1 PROVEN reception inform the G4a framing.
3. **Open questions from Sprints TS / MR / G3 / Connes-Step-2 are not yet stable**. The master Mellin engine (Sprint MR-A/B/C) just settled; the L2 quantitative rate identifying 4/π as the M1 signature is an active interpretation; G3's 1+1 = 2 ℤ₂-grading reading needs absorption. Opening G4a layers on top of unsettled foundations.

### §7.4 Conditional pathway

If PI decides to open G4a:

- **Default sprint timing:** open after Paper 38 lands (estimated 4–8 weeks from 2026-05-06, i.e., June or July 2026).
- **Sprint duration:** 1.5–2 months as estimated in §6.4.
- **Predicted verdict:** POSITIVE-THIN, sharper than H1. The result extends Marcolli–vS-without-Higgs from electroweak (H1) to the full SM gauge group.
- **Paper 32 §VIII.D addendum:** ~150–250 lines extending §VIII.C. No new Paper 39.
- **WH1 status update:** WH1 stays PROVEN. The G4a closure adds a corollary at the SM scale, not a new theorem.

### §7.5 Alternative recommendations

If PI prefers a different next sprint:

- **Alternative 1: G4b structural NCG framework work**. Multi-month, requires cross-category NCG mathematics that doesn't currently exist. Not sprint-scale. Recommendation: don't open.
- **Alternative 2: Extend Connes-Step-2 (J at finite n_max) to AC factor checks at H_F dim 32**. This is a sub-component of G4a; doing it standalone is ~2 weeks, and it would inform whether order-one at the SM scale has any surprises before opening the full G4a. **Worth considering as a 2-week pre-G4a sanity sprint.**
- **Alternative 3: WH4 four-way S³ coincidence formal proof (or attempt)**. The strongest interpretive claim in the project; making it a theorem would change the framing of the entire program. Multi-month, speculative. Not recommended as next sprint.
- **Alternative 4: Master Mellin engine domain partition (Sprint MR follow-up)**. Sprint MR-A's structural finding (k=0 ↔ propinquity, k=1 ↔ vertex parity, k=2 ↔ heat kernel) is a candidate paper-grade observation. A small follow-up sprint to nail down the prediction structure for MR-A's null result and write up the partition as a Paper 18 §III.7 extension is 2–3 weeks. **Worth considering as an immediate next sprint after Paper 38.**

### §7.6 Summary

| Dimension | Recommendation |
|-----------|----------------|
| Open G4a now? | **No.** |
| Open G4a after Paper 38 lands and PI confirms framing? | **Yes, conditional on §7.1 conditions.** |
| Sprint duration estimate | 1.5–2 months |
| Predicted verdict | POSITIVE-THIN, sharper than H1 (covers full SM gauge content) |
| Risk profile | Low to medium; main risk is scope creep on multi-generation / Paper 39 framing |
| Paper venue | Paper 32 §VIII.D addendum, not a standalone paper |
| Most useful pre-sprint check | Confirm A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ) embedding details from Connes–Marcolli 2008 §13.5 (1–2 day literature pull) |

---

## §8. Files referenced

- `geovac/almost_commutative.py` (Sprint H1, ~854 lines, 53 tests passing)
- `tests/test_almost_commutative.py`
- `geovac/full_dirac_operator_system.py` (T_{S³} factor, R3.5 closed)
- `geovac/real_structure.py` (J_GV at KO-3, verified)
- `geovac/chirality_grading.py` (G3-A γ_GV = σ_x ⊗ I)
- `debug/h1_ac_extension_memo.md`
- `debug/g3a_chirality_memo.md`
- `debug/g3b_chirality_F_audit_memo.md`
- `debug/g3c_tensor_chirality_memo.md` (referenced; not re-read in this scoping)
- `debug/g4_cross_manifold_scoping_memo.md`
- `papers/synthesis/paper_32_spectral_triple.tex` §VIII.B G2/G3/G4

External: Connes & Marcolli (2008), *Noncommutative Geometry, Quantum Fields and Motives*, Ch. 13. Connes (1996), J. Math. Phys. 36. Chamseddine, Connes, Marcolli (2007), Adv. Theor. Math. Phys. 11. van Suijlekom (2015), Ch. 8–9.

---

**End of memo.** Word count: ~4,800 words.
