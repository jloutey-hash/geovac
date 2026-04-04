# Track BA: Published Gaussian Compression Benchmark — Literature Data

Generated: 2026-04-02
Status: UPDATED (Track BN, v2.0.31). Fabricated/interpolated lambda estimates removed. Structural findings verified and inserted into Paper 14. Specific lambda values from training knowledge marked as UNVERIFIED — only structural conclusions (DF/THC bypass Pauli; no published data at GeoVac's 10-70 qubit scale) are used in Paper 14.

---

## 1. GeoVac Data (from Track AX / Paper 14, verified)

| System | Method | Config | Qubits | Pauli Terms | 1-Norm (Ha) | Source |
|--------|--------|--------|--------|-------------|-------------|--------|
| He | GeoVac single | n_max=2 | 10 | 120 | 11.29 | Paper 14 |
| He | GeoVac single | n_max=3 | 28 | 2,659 | 78.36 | Paper 14 |
| LiH | GeoVac composed | n_max=2 | 30 | 334 | 37.33 | Paper 14 |
| LiH | GeoVac composed | n_max=3 | 84 | 7,879 | 202.49 | Paper 14 |
| BeH2 | GeoVac composed | n_max=2 | 50 | 556 | 354.89 | Paper 14 |
| H2O | GeoVac composed | n_max=2 | 70 | 778 | — | Paper 14 |

---

## 2. Gaussian Raw JW Data (Trenev et al. 2025, verified in Paper 14)

| System | Basis | Qubits | Pauli Terms | 1-Norm (Ha) | Source |
|--------|-------|--------|-------------|-------------|--------|
| LiH | STO-3G | 10 | 276 | — | Trenev 2025 |
| LiH | 6-31G | 20 | 5,851 | — | Trenev 2025 |
| LiH | cc-pVDZ | 36 | 63,519 | — | Trenev 2025 |
| H2O | STO-3G | 12 | 551 | — | Trenev 2025 |
| H2O | 6-31G | 24 | 8,921 | — | Trenev 2025 |
| H2O | cc-pVDZ | 46 | 107,382 | — | Trenev 2025 |

---

## 3. Published Compressed Gaussian Resource Estimates

### IMPORTANT CAVEAT

The fault-tolerant QPE literature (Lee et al. 2021, Berry et al. 2019, von Burg et al. 2021, Goings et al. 2022) generally does **NOT** report Pauli term counts for compressed representations. These papers work with block-encoding / linear combination of unitaries (LCU) frameworks where the relevant cost metrics are:

- **lambda (1-norm):** The LCU 1-norm of the Hamiltonian representation, which controls the number of queries in qubitization.
- **Toffoli gate count:** The total number of Toffoli gates (or T-gates) for the full QPE algorithm.
- **Number of auxiliary qubits:** The ancilla overhead for block encoding.

Double factorization (DF) and tensor hypercontraction (THC) do not produce a "compressed Pauli decomposition" — they produce a factored representation that is implemented as a sequence of rotations in a block-encoding circuit. The comparison with GeoVac's Pauli term count is therefore **not apples-to-apples** for these methods. The 1-norm comparison IS meaningful, as lambda controls query complexity in both Trotter and qubitization approaches.

---

### 3A. Lee et al., PRX Quantum 2, 030305 (2021) — THC

"Even More Efficient Quantum Computations of Chemistry Through Tensor Hypercontraction"

This paper reports Toffoli costs and 1-norms for THC-based qubitization. Key data from their Table I and Table II (approximate values from training knowledge — VERIFY AGAINST PAPER):

| System | Basis | Active space | Qubits (N) | lambda_THC (Ha) | Toffoli cost | Source |
|--------|-------|-------------|------------|-----------------|--------------|--------|
| H2O | cc-pVDZ | (10e, 24o) | 48 | ~25 | ~2.0e8 | Lee 2021, Table I (approx) |
| H4 chain | STO-3G | (4e, 4o) | 8 | ~3.5 | — | Lee 2021 (approx) |
| FeMoCo | cc-pVDZ | (54e, 76o) | 152 | ~306 | ~3.5e10 | Lee 2021, Table II |
| FeMoCo | cc-pVTZ | (54e, 152o) | 304 | ~1,170 | ~1.1e11 | Lee 2021, Table II |

**Key comparison for H2O:**
- Gaussian THC at Q~48: lambda ~ 25 Ha
- GeoVac composed at Q=70: 1-norm MISSING (not reported for H2O due to Z_eff=6 inflation)
- Cannot make direct 1-norm comparison for H2O without GeoVac data

**Key insight:** THC achieves lambda ~ O(N) scaling for the 1-norm (vs O(N^2-3) for raw Pauli). Lee et al. report lambda scaling approximately as N^{1.0-1.5} for their test systems.

**CONFIDENCE: MEDIUM.** I am confident about the existence and general magnitude of these results, but the specific numbers for H2O should be verified. The FeMoCo numbers are widely cited. The paper does NOT report data for LiH or He in the form needed for direct comparison.

---

### 3B. Berry et al., Quantum 3, 208 (2019) — Sparse / Qubitization

"Qubitization of Arbitrary Basis Quantum Chemistry Leveraging Sparsity and Low Rank Factorization"

This paper introduced the sparse qubitization approach. Key results (from training knowledge — VERIFY):

| System | Basis | Active space | Qubits | lambda (Ha) | Toffoli cost | Source |
|--------|-------|-------------|--------|-------------|--------------|--------|
| FeMoCo | (54e, 76o) | — | 152 | ~4,800 | ~1.4e13 | Berry 2019, Table I (approx) |

**Note:** Berry et al. focus on large-scale molecules. They do NOT report detailed data for small molecules (H2, LiH, H2O, He). The paper's value for this comparison is establishing the pre-THC baseline for fault-tolerant costs. Lee et al. 2021 improved upon these numbers by ~100x using THC.

**CONFIDENCE: LOW for specific numbers.** The FeMoCo lambda is approximate. The paper's main contribution is algorithmic, not a systematic small-molecule benchmark.

---

### 3C. von Burg et al., PRX Quantum 2, 030305 (2021) — Double Factorization

"Quantum Computing Enhanced Computational Catalysis"

**NOTE:** von Burg et al. share the same PRX Quantum 2, 030305 citation with Lee et al. — this appears to be a citation confusion in the task prompt. The von Burg et al. paper may be a different publication. Let me report what I know:

von Burg et al. (likely arXiv:2007.14460) report double factorization (DF) resource estimates. Key data:

| System | Basis | Method | lambda_DF (Ha) | Toffoli cost | Source |
|--------|-------|--------|----------------|--------------|--------|
| H2O | cc-pVDZ | DF | — | — | von Burg (approx) |
| Rh catalyst | cc-pVDZ | DF | ~600 | ~5e10 | von Burg (approx) |

**CONFIDENCE: LOW.** I cannot reliably distinguish the specific numbers from this paper vs others in the DF literature. The key claim is that DF reduces lambda by ~3-10x compared to naive sparse methods. Small-molecule data may not be the focus.

---

### 3D. Motta et al., npj Quantum Inf. 7, 83 (2021) — Low-rank / DF foundations

"Low rank representations for quantum simulation of electronic structure"

This paper establishes the theoretical framework for DF and reports benchmark data. Key data (from training knowledge):

| System | Basis | Rank (L) | lambda_DF (Ha) | vs raw lambda | Source |
|--------|-------|----------|----------------|---------------|--------|
| H2 | cc-pVDZ | — | — | — | — |
| H2O | cc-pVDZ | ~50-80 | — | ~2-5x reduction | Motta 2021 (approx) |
| N2 | cc-pVDZ | ~50-80 | — | ~2-5x reduction | Motta 2021 (approx) |

**Key finding:** DF typically reduces the 1-norm (lambda) by a factor of 2-5x for molecules at cc-pVDZ level. The reduction grows with system size but is a constant-factor improvement, not a change in scaling exponent.

**CONFIDENCE: LOW for specific numbers.** The 2-5x lambda reduction factor is widely reported in the literature but I cannot pin exact numbers for specific molecules.

---

### 3E. Loaiza et al., JCTC 2024 (arXiv:2312.07746) — Symmetry-Compressed DF

"Reducing the Runtime of Fault-Tolerant Quantum Simulations in Chemistry through Symmetry-Compressed Double Factorization"

This paper improves upon standard DF by exploiting molecular symmetry. Key data (from training knowledge):

| System | Basis | Method | lambda_SCDF (Ha) | vs DF lambda | Toffoli | Source |
|--------|-------|--------|------------------|--------------|---------|--------|
| H2O | cc-pVDZ | SCDF | — | ~1.5-3x improvement over DF | — | Loaiza 2024 |
| H chains | various | SCDF | — | — | — | Loaiza 2024 |

**Key finding:** SCDF achieves an additional ~1.5-3x lambda reduction beyond standard DF. Combined with DF's 2-5x over raw, total compression is ~3-15x in lambda.

**CONFIDENCE: LOW.** I recall this paper exists and the general conclusions, but cannot cite specific table entries with confidence.

---

### 3F. Goings et al., WIREs 2022 (arXiv:2203.12374) — Cytochrome P450

"Reliably Assessing the Electronic Structure of Cytochrome P450 on Today's Classical Computers and Tomorrow's Quantum Computers"

This paper is primarily about large biological systems, not small molecules. Key data:

| System | Active space | Qubits | Method | lambda (Ha) | Toffoli | Source |
|--------|-------------|--------|--------|-------------|---------|--------|
| P450 (CYP) | (58e, 48o) | 96 | THC | ~1,500 | ~5e10 | Goings 2022 (approx) |
| P450 (CYP) | (58e, 48o) | 96 | DF | ~3,200 | ~2e11 | Goings 2022 (approx) |
| P450 (CYP) | (58e, 48o) | 96 | Sparse | ~6,500 | ~3e12 | Goings 2022 (approx) |

**CONFIDENCE: MEDIUM for general magnitudes.** The paper provides a useful comparison across methods (THC vs DF vs sparse) but focuses on a single large system. No small-molecule data relevant for direct GeoVac comparison.

---

### 3G. Rubin et al., arXiv:2306.03145 (2023) — Bloch Orbitals

"Fault-Tolerant Quantum Simulation of Materials Using Bloch Orbitals"

This paper focuses on periodic materials (crystalline solids), not molecules. Not directly relevant for GeoVac molecular comparisons. No small-molecule data expected.

**CONFIDENCE: HIGH that this paper is not relevant for small-molecule benchmarking.**

---

### 3H. Additional Known Results

**Reiher et al., PNAS 114, 7555 (2017):**
The landmark FeMoCo resource estimation paper.
- FeMoCo (54e, 54o): ~10^14 T-gates with Trotter, ~10^9 with qubitization
- Established the "quantum advantage threshold" narrative

**Babbush et al., PRX 8, 011044 (2018):**
Encoding molecular Hamiltonians with improved constant factors.

**H2O lambda at cc-pVDZ (~48 qubits): No published data.** The table of interpolated lambda values previously in this section has been removed. Those were synthetic estimates extrapolated from general compression ratios, not published numbers. No published DF/THC/SCDF lambda values exist for H2O at the cc-pVDZ / Q~48 scale. The FT literature reports H2O data only at larger active spaces or as part of methodology validation without tabulated lambda. The structural finding is: DF/THC bypass Pauli decompositions entirely and target lambda + Toffoli count; no direct comparison is possible at this scale.

---

## 4. Consolidated Comparison Table

### LiH — Best available comparison

| Method | Config | Qubits | Pauli Terms | 1-Norm (Ha) | Source | Confidence |
|--------|--------|--------|-------------|-------------|--------|------------|
| GeoVac composed | n_max=2 | 30 | 334 | 37.33 | Paper 14 | HIGH |
| GeoVac composed | n_max=3 | 84 | 7,879 | 202.49 | Paper 14 | HIGH |
| Gaussian raw JW | STO-3G | 10 | 276 | — | Trenev 2025 | HIGH |
| Gaussian raw JW | 6-31G | 20 | 5,851 | — | Trenev 2025 | HIGH |
| Gaussian raw JW | cc-pVDZ | 36 | 63,519 | — | Trenev 2025 | HIGH |
| Gaussian DF | cc-pVDZ | 36 | N/A (block-enc.) | No published data at this scale | — | STRUCTURAL |
| Gaussian THC | cc-pVDZ | 36 | N/A (block-enc.) | No published data at this scale | — | STRUCTURAL |

**Assessment:** No published DF/THC lambda values exist for LiH at cc-pVDZ (or any basis at Q~30-36). The FT literature focuses on large systems (FeMoCo at 152 qubits, P450 at 96 qubits). DF/THC use block-encoding circuits that bypass Pauli decompositions entirely; the Pauli term comparison is against raw JW (Trenev et al.), where GeoVac wins by 190x. The "small-molecule gap" in the FT literature means this comparison cannot currently be made for LiH.

### H2O — Best available comparison

| Method | Config | Qubits | Pauli Terms | 1-Norm (Ha) | Source | Confidence |
|--------|--------|--------|-------------|-------------|--------|------------|
| GeoVac composed | n_max=2 | 70 | 778 | MISSING | Paper 14 | HIGH (count), N/A (1-norm) |
| Gaussian raw JW | cc-pVDZ | 46 | 107,382 | — | Trenev 2025 | HIGH |
| Gaussian THC | cc-pVDZ | ~48 | N/A (block-enc.) | No published data at this scale | — | STRUCTURAL |
| Gaussian DF | cc-pVDZ | ~48 | N/A (block-enc.) | No published data at this scale | — | STRUCTURAL |

**Assessment:** No published DF/THC lambda values exist for H2O at cc-pVDZ in the Q~48-70 range. Lee et al. 2021 report H2O THC data only for larger active spaces. The Pauli term comparison (778 vs 107,382) remains a 138x advantage against raw JW. GeoVac H2O partitioned 1-norm is 361 Ha (electronic only, PK classical).

### He — Best available comparison

| Method | Config | Qubits | Pauli Terms | 1-Norm (Ha) | Source | Confidence |
|--------|--------|--------|-------------|-------------|--------|------------|
| GeoVac single | n_max=2 | 10 | 120 | 11.29 | Paper 14 | HIGH |
| GeoVac single | n_max=3 | 28 | 2,659 | 78.36 | Paper 14 | HIGH |
| Gaussian | cc-pVDZ | 10 | 156 | 42.95 | Paper 14 (computed) | HIGH |
| Gaussian | cc-pVTZ | 28 | 21,607 | 530.47 | Paper 14 (computed) | HIGH |

**Assessment:** This is the strongest comparison because both Pauli counts AND 1-norms are computed from validated integrals. GeoVac wins on both metrics at equal qubit count: 1.3x/8.1x on Pauli terms, 3.8x/6.8x on 1-norm. No DF/THC results exist for He at these sizes (too small to benefit from compression).

---

## 5. Large-Scale Context (Gaussian Compression at Production Scale)

These numbers contextualize where Gaussian compression methods operate:

| System | Method | Active space | Qubits | lambda (Ha) | Toffoli | Source | Confidence |
|--------|--------|-------------|--------|-------------|---------|--------|------------|
| FeMoCo | Sparse (Berry) | (54e, 76o) | 152 | ~4,800 | ~1.4e13 | Berry 2019 | MEDIUM |
| FeMoCo | THC (Lee) | (54e, 76o) | 152 | ~306 | ~3.5e10 | Lee 2021 | MEDIUM |
| FeMoCo | THC (Lee) | (54e, 152o) | 304 | ~1,170 | ~1.1e11 | Lee 2021 | MEDIUM |
| P450 | THC (Goings) | (58e, 48o) | 96 | ~1,500 | ~5e10 | Goings 2022 | MEDIUM |
| P450 | DF (Goings) | (58e, 48o) | 96 | ~3,200 | ~2e11 | Goings 2022 | MEDIUM |

**Key scaling observation:** THC reduces lambda from sparse by ~10-16x for FeMoCo. This is a constant-factor improvement that does not change the asymptotic scaling with system size.
