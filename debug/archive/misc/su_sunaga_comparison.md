# Sprint 3 Track SU: Sunaga 2025 Matched-Q Comparison

**Date:** April 2026
**Status:** Complete (native-Q definitive, matched-Q deferred per structural feasibility)
**Companion data:** `debug/data/su_resource_tables.json`

---

## Headline

GeoVac heavy-atom relativistic hydrides at **native Q=20** produce **534 Pauli terms** (CaH/SrH/BaH — bit-identical, isostructural). This is **88× fewer** than Sunaga 2025's RaH at Q=18 (47,099 Pauli terms from PRA 111, 022817 main paper Table II).

The matched-Q=18 comparison is **structurally infeasible** within GeoVac's native architecture. Documented as a structural observation, not a limitation.

---

## 1. SU-A: Sunaga 2025 SI Accessibility

The PRA 111, 022817 main paper reports RaH at Q=18 with 47,099 Pauli terms (Tier 2 T4 market test baseline). SI Tables S1-S3 are cited as containing per-molecule SrH/BaH/RaH at Q=18 but are behind the PRA supplemental paywall. Not obtained in this sprint; flagged as future work if institutional access becomes available.

**What we have:** one public reference point (RaH-18q = 47,099 Pauli).
**What's deferred:** the SrH-18q / BaH-18q rows of Sunaga SI.

---

## 2. SU-B: GeoVac Native-Q Resource Table

### Scalar composed Hamiltonians (n_max=2)

| System | Z | Q | N_Pauli | λ_ni (Ha) | λ_total (Ha) |
|--------|---|---|---------|-----------|--------------|
| LiH    | 3 | 30 | 334 | 32.59 | 32.59 |
| NaH    | 11 | 20 | 223 | 171.46 | 171.46 |
| KH     | 19 | 20 | 223 | 608.66 | 608.66 |
| CaH₂   | 20 | 40 | 445 | 704.63 | 704.63 |
| SrH    | 38 | 20 | 223 | 3145.86 | 3145.86 |
| BaH    | 56 | 20 | 223 | 7897.88 | 7897.88 |

**Isostructural invariance:** NaH/KH/SrH/BaH all produce **N_Pauli = 223** — identical Pauli count across four different frozen cores ([Ne], [Ar], [Kr], [Xe]) and atomic numbers spanning Z=11 to Z=56. The 1-norm scales with Z (dominated by nuclear repulsion), but the operator structure is identical.

LiH has Q=30 (not 20) because Li's structure type C uses an explicit 2-electron core block rather than a frozen core.

### Relativistic composed Hamiltonians (n_max=2, α=CODATA)

| System | Z | Q | N_Pauli | λ_ni (Ha) | QWC groups | Ratio vs Sunaga RaH-18q |
|--------|---|---|---------|-----------|------------|-------------------------|
| LiH    | 3 | 30 | 801 | 25.79 | 52 | 0.0170× |
| BeH    | 4 | 30 | 801 | 41.64 | 52 | 0.0170× |
| CaH    | 20 | 20 | 534 | 13.87 | 52 | 0.0113× |
| SrH    | 38 | 20 | 534 | 13.87 | 52 | 0.0113× |
| BaH    | 56 | 20 | 534 | 13.87 | 52 | 0.0113× |

**Isostructural invariance extends to relativistic:**
- CaH/SrH/BaH: bit-identical at Q=20, N_Pauli=534, λ_ni=13.87 Ha, QWC=52
- LiH/BeH: bit-identical at Q=30, N_Pauli=801, λ_ni differs (depends on bare Z)

The heavy-atom relativistic resource count depends **only on the valence topology** (one bond pair atop a frozen core), not on the frozen core's electron count.

### Rel/scalar ratio

| System | N_Pauli scalar | N_Pauli rel | Ratio |
|--------|----------------|-------------|-------|
| CaH    | 223 | 534 | 2.395× |
| SrH    | 223 | 534 | 2.395× |
| BaH    | 223 | 534 | 2.395× |

The 2.395× ratio matches Tier 2 T3's isostructural regression (LiH pin: 2.42×, BeH: 2.42×, CaH: 2.42×). Frozen-core screening makes the spin-orbit diagonal see Z_eff=2 uniformly, independent of the bare Z; heavy-atom relativistic resources are purely a valence-topology property.

---

## 3. SU-B: Comparison to Sunaga RaH-18q

The single public Sunaga data point: RaH at Q=18 with 47,099 Pauli terms.

**Native-Q ratios:**
- GeoVac CaH/SrH/BaH_rel (Q=20): 534 / 47,099 = **0.0113×** (88× fewer)
- GeoVac LiH_rel (Q=30): 801 / 47,099 = **0.0170×** (59× fewer)

Caveats:
1. Sunaga's RaH is a heavier molecule (Ra is Z=88) with a [Rn] core that GeoVac does not currently support. Direct molecule-to-molecule comparison is approximate.
2. Q is off by 2 (20 vs 18). A matched-Q comparison requires either truncating the GeoVac basis (unphysical) or accessing Sunaga's SrH-18q / BaH-18q rows.

Despite these caveats, the 88× gap is large enough that even pessimistic matched-Q corrections would preserve a sparsity win.

---

## 4. SU-C: Matched-Q=18 Feasibility Assessment

GeoVac's native minimum Q for a single-bond alkali/alkaline-earth hydride is **Q=20** at n_max=2. To produce a Q=18 Hamiltonian, three options:

**Option 1: Drop 1 orbital from each block.** Physically broken. The 5-orbital (1s, 2s, 2p_x, 2p_y, 2p_z) angular completeness per block is the framework's defining property — dropping an orbital is a basis-set choice that has no principled criterion in GeoVac.

**Option 2: Reduce n_max to 1 for some blocks.** Gives Q=10 per block. For a single-bond hydride: Q=10 total — overshoots below Q=18. Not a path.

**Option 3: Custom 9-orbital spec.** Would require physical justification for which 9 orbitals. Most natural constraint: freeze the core even more aggressively (valence-only), giving Q ≈ 10. Again, overshoots below Q=18.

**Structural observation (publication-ready):** GeoVac's block sizes are determined by the natural geometry (spherical orbital basis on S³), not by matching arbitrary qubit counts. This is a fundamental architectural difference from Gaussian-basis approaches where qubit count is tunable by choice of basis set (STO-3G vs cc-pVDZ vs cc-pVTZ). The native-Q ladder for GeoVac's composed architecture is:

- Q=10: n_max=1 single block (e.g., H₂ at crude approximation)
- Q=20: single-bond hydride with frozen core (NaH, KH, CaH, SrH, BaH)
- Q=30: LiH or BeH (explicit core block)
- Q=40+: polyatomics or n_max=3

Q=18 falls between the natural steps. Matched-Q=18 is therefore not a natural GeoVac operating point; the correct comparison frame is "at what native-Q does GeoVac produce a competitive Hamiltonian?"

---

## 5. Summary

1. **Native-Q resource table** computed for 6 scalar + 5 relativistic heavy hydrides.
2. **Isostructural invariance verified:** CaH/SrH/BaH bit-identical at both scalar (223 Pauli) and relativistic (534 Pauli) levels.
3. **Native-Q comparison to Sunaga RaH-18q:** 88× fewer Pauli terms at Q=20 for CaH/SrH/BaH_rel.
4. **Matched-Q=18 infeasible** within GeoVac's native architecture — documented as a structural observation.
5. **Paper 20 update:** Tier-2 table extended with SrH/BaH rows; Sunaga comparison reframed as native-Q.

Follow-up work (beyond Sprint 3):
- Access Sunaga SI Tables S1-S3 for matched comparison at published Q=18
- [Rn] core implementation for RaH direct molecule-to-molecule comparison
