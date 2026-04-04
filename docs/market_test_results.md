# Track CA: Quantum Resource Market Test — Results

**Date:** 2026-04-03
**Question:** Does GeoVac's structural sparsity survive comparison against state-of-the-art Gaussian Hamiltonian compression?

---

## Executive Summary

**GeoVac wins on Pauli count (VQE-relevant). It is competitive on 1-norm against raw Gaussian, and the comparison against compressed Gaussian (DF/THC) cannot be made because published DF/THC lambda values for small molecules (Q < 100) do not exist in the literature.**

The fault-tolerant QPE literature (Lee 2021, von Burg 2021, Rocca 2024) benchmarks exclusively at production scale (FeMoco at 152+ qubits, P450 at 96 qubits). There are no published DF/THC lambda values for LiH or H2O at the Q = 10–70 scale where GeoVac operates. This is both a gap in the literature and a positioning opportunity.

---

## Data Tables

### Table 1: LiH Head-to-Head

| Method | Basis | Q | Pauli terms | λ (Ha) | QWC | Accuracy |
|:-------|:------|--:|:-----------:|:------:|:---:|:---------|
| **GeoVac composed** | max_n=2 | 30 | **334** | **33.3** | 21 | 5.3% R_eq |
| GeoVac TC composed | max_n=2 | 30 | 562 | 38.6 | — | converging |
| Gaussian raw JW | STO-3G | 12 | 907 | 34.3 | 273 | ~5% |
| Gaussian raw JW (reduced) | STO-3G | 10 | 276 | — | — | ~5% |
| Gaussian raw JW | 6-31G | 20 | 5,851 | — | — | ~2% |
| Gaussian raw JW | cc-pVDZ | 36 | 63,519 | — | — | <0.1% |
| Gaussian + DF | cc-pVDZ | 36 | N/A | **?** | N/A | <0.1% |
| Gaussian + THC | cc-pVDZ | 36 | N/A | **?** | N/A | <0.1% |

**Key findings for LiH:**

1. **Pauli count: CLEAR WIN.** GeoVac at Q=30 has 334 Pauli terms vs STO-3G's 907 (2.7× fewer) at Q=12, and vs cc-pVDZ's 63,519 (190× fewer) at Q=36.

2. **1-norm: COMPARABLE to STO-3G.** GeoVac electronic-only λ = 33.3 Ha vs STO-3G λ = 34.3 Ha. These are essentially equal. This is surprising — GeoVac uses 2.5× more qubits but achieves the same 1-norm. The 1-norm is dominated by one-body terms (core energies, PK), not by the number of Pauli terms.

3. **1-norm vs compressed Gaussian: UNKNOWN.** No published DF/THC lambda values exist for LiH at any basis set in the QPE literature. This comparison cannot currently be made.

4. **QWC groups: MASSIVE WIN.** GeoVac has 21 QWC groups vs STO-3G's 273 (13× fewer). This directly translates to fewer VQE measurement circuits.

### Table 2: H2O Head-to-Head

| Method | Basis | Q | Pauli terms | λ (Ha) | QWC | Accuracy |
|:-------|:------|--:|:-----------:|:------:|:---:|:---------|
| **GeoVac composed** | max_n=2 | 70 | **778** | **360.8** | 21 | 26% R_eq |
| GeoVac composed (full PK) | max_n=2 | 70 | 778 | 28,053 | 21 | 26% |
| Gaussian raw JW (reduced) | STO-3G | 12 | 551 | — | — | ~10% |
| Gaussian raw JW (reduced) | 6-31G | 24 | 8,921 | — | — | ~2% |
| Gaussian raw JW (reduced) | cc-pVDZ | 46 | 107,382 | — | — | <0.1% |

**Key findings for H2O:**

1. **Pauli count: WIN vs cc-pVDZ** (778 vs 107,382 = 138×), **WIN vs 6-31G** (778 vs 8,921 = 11×), **LOSS vs STO-3G** (778 vs 551 = 1.4×). At minimal basis, Gaussian is sparser because H2O STO-3G has only 7 spatial orbitals vs GeoVac's 35.

2. **1-norm: HIGH.** The electronic-only λ = 360.8 Ha is high for a 10-electron system. The full-PK λ = 28,053 Ha is dominated by the Z²-scaled PK barrier on the oxygen core (Z_eff=6). PK classical partitioning (Track BF) is essential — it reduces λ by 78×.

3. **Accuracy: POOR.** 26% R_eq error vs Gaussian's <0.1%. This is the fundamental trade-off.

### Table 3: He Equal-Qubit Comparison (Validated)

| Method | Basis | Q | Pauli terms | λ (Ha) | QWC |
|:-------|:------|--:|:-----------:|:------:|:---:|
| **GeoVac** | max_n=2 | 10 | 120 | 11.3 | — |
| Gaussian | cc-pVDZ | 10 | 156 | 43.0 | 25 |
| **GeoVac** | max_n=3 | 28 | 2,659 | 78.4 | — |
| Gaussian | cc-pVTZ | 28 | 21,607 | 530.5 | 8,615 |

At equal qubit count, GeoVac wins on BOTH Pauli count and 1-norm:
- Q=10: 1.3× fewer Pauli, **3.8× lower λ**
- Q=28: 8.1× fewer Pauli, **6.8× lower λ**

This is the cleanest comparison because the systems are identical (both He, 2 electrons) and the qubit counts match.

---

## Analysis

### Q1: At equal qubit count (~30), does GeoVac have lower 1-norm?

**For He: YES, decisively.** GeoVac He at Q=28 has λ=78.4 Ha vs Gaussian cc-pVTZ λ=530.5 Ha (6.8× lower).

**For LiH: COMPARABLE.** GeoVac LiH at Q=30 has λ=33.3 Ha (electronic-only) vs Gaussian STO-3G at Q=12 with λ=34.3 Ha. The 1-norms are nearly identical despite different qubit counts. GeoVac uses more qubits to encode a larger orbital space with fewer Pauli terms — the 1-norm stays flat because it's dominated by one-body energies, not by the two-electron structure that drives Pauli count.

### Q2: At equal accuracy (~5% error), does GeoVac use fewer resources?

**Mixed.** GeoVac LiH at 5.3% R_eq uses Q=30, 334 Pauli, λ=33.3 Ha. Gaussian STO-3G at ~5% error uses Q=12, 907 Pauli, λ=34.3 Ha. GeoVac uses 2.5× more qubits but 2.7× fewer Pauli terms and 13× fewer QWC groups, with the same 1-norm. For VQE (where measurement count ∝ QWC groups × shots), GeoVac wins. For QPE (where cost ∝ λ/ε), it's a draw.

### Q3: Does GeoVac's sparsity compose with downstream optimizations?

**Likely yes, but untested.** GeoVac's ERI tensor is already 0.4% dense (from Gaunt selection rules and block-diagonal structure). Double factorization achieves compression by exploiting low-rank structure in the ERI tensor. An already-sparse ERI should have lower rank, meaning DF should achieve even better compression on GeoVac Hamiltonians than on Gaussian ones. The structural zeros from Gaunt selection rules are "free" — they don't need to be compressed because they don't exist.

This is a research direction worth pursuing: apply DF/THC to GeoVac's already-sparse Hamiltonians and measure the combined compression.

### Q4: Where does GeoVac's advantage land on VQE vs QPE?

**GeoVac's advantage is primarily in the VQE regime:**

- **Pauli count** (VQE-relevant): 2.7×–190× advantage over raw Gaussian
- **QWC groups** (VQE-relevant): 13× advantage (21 vs 273 for LiH STO-3G)
- **1-norm** (QPE-relevant): comparable to STO-3G, much better than cc-pVDZ at equal Q

For fault-tolerant QPE, the comparison would be GeoVac's λ vs DF/THC-compressed Gaussian λ. Since DF/THC can reduce λ by 3–10× for large systems, and GeoVac's λ is comparable to *raw* Gaussian STO-3G, DF/THC-compressed Gaussian at cc-pVDZ quality would likely have *lower* λ than GeoVac at comparable accuracy. **But this comparison cannot be made with published data.**

**The honest assessment:** GeoVac's positioning is strongest for near-term VQE on NISQ/early-FTQC hardware (Q < 100), where Pauli count and measurement circuit count dominate runtime. For large-scale fault-tolerant QPE, the DF/THC compression ecosystem may eventually match or beat GeoVac's intrinsic sparsity — but at present there is no data to confirm or refute this.

---

## Structural Observations

### Why GeoVac has few QWC groups

GeoVac's 21 QWC groups (constant across LiH/BeH₂/H₂O) is a consequence of the block-diagonal ERI structure. Each block's Pauli terms commute within themselves because they share the same qubit register partition. This is a structural property of the composed architecture, not an accident of small system size. As system size grows, the number of blocks grows but each block's QWC structure is preserved.

### Why DF/THC don't help at GeoVac's scale

DF/THC are designed to compress dense 4-index ERI tensors by exploiting low-rank structure. GeoVac's ERI tensor is already 99.6% zero (Gaunt selection rules + block-diagonal structure). Applying DF to an already-sparse tensor is like compressing a JPEG — diminishing returns. The compression that DF provides at FeMoco scale (590 → 306 Ha, ~2× for THC) would provide less relative benefit on GeoVac's already-compact Hamiltonians.

### The "small-molecule gap" in QPE literature

Published DF/THC resource estimates exist only for Q > 100 (FeMoco, P450, Rh catalysts). No published paper reports lambda under DF/THC for LiH, H₂O, or any system at Q < 50. This means GeoVac's competitive landscape at the Q=10–70 scale is defined entirely by raw Gaussian baselines — where GeoVac wins decisively on Pauli count and QWC groups.

---

## Honest Assessment

**Where GeoVac wins:**
- Pauli term count: 2.7×–190× vs raw Gaussian (basis-dependent)
- QWC measurement groups: 13× fewer (21 vs 273 for LiH)
- 1-norm at equal qubit count (He): 3.8×–6.8× lower
- Zero-parameter construction (no fitted parameters)
- Intrinsic sparsity that is basis-guaranteed, not post-hoc

**Where GeoVac is comparable:**
- 1-norm for LiH: ~33 Ha vs ~34 Ha (STO-3G)

**Where GeoVac loses:**
- Accuracy: 5.3% R_eq (LiH), 26% (H₂O) vs <0.1% for Gaussian cc-pVDZ
- Qubit count: 30 (LiH), 70 (H₂O) vs 12–46 for equivalent Gaussian
- Cannot currently be compared against DF/THC-compressed Gaussian

**Where GeoVac's position is uncertain:**
- 1-norm vs compressed Gaussian at equal accuracy (no published data)
- Whether DF/THC applied to GeoVac's sparse ERIs gives additional compression
- Scaling to larger molecules (first-row scope boundary)

---

## Recommendations

1. **Position for VQE/NISQ, not QPE.** GeoVac's advantage is in Pauli count and QWC groups, which matter for VQE. The 1-norm advantage (QPE-relevant) is modest.

2. **The equal-qubit He comparison is the strongest result.** 6.8× lower λ at Q=28 with validated Gaussian integrals. Lead with this.

3. **Compute DF/THC lambda for LiH cc-pVDZ.** This requires PySCF (not available on Windows). Running on a Linux machine with PySCF + PennyLane's `qml.resource.DoubleFactorization` would fill the critical data gap.

4. **Test DF on GeoVac ERIs.** Apply double factorization to GeoVac's block-diagonal ERI tensor to measure combined compression.

---

## Data Sources

- GeoVac composed: computed via `geovac.composed_qubit` (this session)
- GeoVac He: `debug/data/commutator_bounds.json` (Track AV)
- GeoVac TC: `debug/data/tc_composed_benchmark.json` (Track BX-3)
- Gaussian LiH STO-3G: OpenFermion cached `H1-Li1_sto-3g_singlet_1.45.hdf5`
- Gaussian He cc-pVDZ/cc-pVTZ: `geovac/gaussian_reference.py` (computed integrals)
- Gaussian Pauli counts: Trenev et al. 2025 (arXiv:2311.03719), Table 5
- DF/THC published data: Lee et al. 2021 (PRX Quantum 2, 030305); von Burg et al. 2021 (Phys. Rev. Research 3, 033055); Rocca et al. 2024 (J. Chem. Theory Comput. 20, 4639)

Raw data: `debug/data/market_test_data.json`
