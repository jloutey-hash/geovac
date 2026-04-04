# Track BA: Source Bibliography and Access Notes

Generated: 2026-04-02
Status: WebSearch and WebFetch were both DENIED. All data below is from training knowledge (pre-May 2025 cutoff). Numbers marked with confidence levels.

---

## Papers Surveyed

### 1. Loaiza et al., JCTC 2024 (arXiv:2312.07746)
**Title:** "Reducing the Runtime of Fault-Tolerant Quantum Simulations in Chemistry through Symmetry-Compressed Double Factorization"
**Access:** DENIED (WebFetch blocked)
**What I recall:** Introduces symmetry-compressed DF (SCDF). Reports lambda reductions of ~1.5-3x beyond standard DF for symmetric molecules. Benchmarks include H chains and possibly H2O. Reports Toffoli costs.
**What was extracted:** General compression factors only (LOW confidence). No specific table entries.
**What's needed:** Full text access to extract Tables with exact lambda values for H2O, H chains at various basis sets.

### 2. Lee et al., PRX Quantum 2, 030305 (2021) (arXiv:2011.03494)
**Title:** "Even More Efficient Quantum Computations of Chemistry Through Tensor Hypercontraction"
**Access:** DENIED (WebFetch blocked)
**What I recall:** The definitive THC resource estimation paper. Reports lambda_THC and Toffoli costs for H2O, H chains, and FeMoCo at cc-pVDZ and cc-pVTZ. Table I has small molecules; Table II has FeMoCo. THC achieves ~10-16x lambda reduction over sparse qubitization for FeMoCo.
**What was extracted:** Approximate FeMoCo lambda (~306 Ha at cc-pVDZ, ~1170 Ha at cc-pVTZ) — MEDIUM confidence. H2O lambda estimated ~20-30 Ha at cc-pVDZ — MEDIUM confidence. Toffoli costs approximate.
**What's needed:** Exact lambda values from Table I for H2O, H4, H10 at cc-pVDZ and cc-pVTZ. These would provide the strongest compressed-Gaussian comparison points.

### 3. Berry et al., Quantum 3, 208 (2019) (arXiv:1902.02134)
**Title:** "Qubitization of Arbitrary Basis Quantum Chemistry Leveraging Sparsity and Low Rank Factorization"
**Access:** DENIED (WebFetch blocked)
**What I recall:** Introduces sparse qubitization with DF. Reports lambda and T-gate costs primarily for FeMoCo. Establishes the baseline that THC (Lee 2021) later improves upon by ~100x in Toffoli count.
**What was extracted:** Approximate FeMoCo lambda (~4800 Ha) — LOW confidence. No small-molecule data recalled.
**What's needed:** Any data for small molecules (H2, H2O, LiH). May not exist in this paper.

### 4. von Burg et al. (arXiv:2007.14460, 2021)
**Title:** "Quantum Computing Enhanced Computational Catalysis"
**Access:** DENIED (WebFetch blocked)
**What I recall:** Reports DF resource estimates for Rh catalyst. May include H2O as a benchmark. The citation in the task prompt (PRX Quantum 2, 030305) conflicts with Lee et al. — these may be different papers with the same journal reference (unlikely; probably a prompt error).
**What was extracted:** Very little — LOW confidence for all numbers.
**What's needed:** Verify correct citation. Check for H2O or LiH lambda values.

### 5. Goings et al., WIREs 2022 (arXiv:2203.12374)
**Title:** "Reliably Assessing the Electronic Structure of Cytochrome P450 on Today's Classical Computers and Tomorrow's Quantum Computers"
**Access:** DENIED (WebFetch blocked)
**What I recall:** Comprehensive resource estimation comparing THC, DF, and sparse methods for cytochrome P450. Reports lambda and Toffoli costs for P450 at various active space sizes. Demonstrates THC ~2x better than DF for this system.
**What was extracted:** Approximate P450 numbers (lambda ~1500 Ha THC, ~3200 Ha DF) — MEDIUM confidence.
**What's needed:** Any small-molecule comparison data. This paper likely does NOT have it — focused on P450.

### 6. Rubin et al., arXiv:2306.03145 (2023)
**Title:** "Fault-Tolerant Quantum Simulation of Materials Using Bloch Orbitals"
**Access:** DENIED (WebFetch blocked)
**What I recall:** Focused on periodic/crystalline materials. Uses Bloch orbital encoding for solid-state Hamiltonians. Not directly relevant for molecular benchmarks.
**What was extracted:** Nothing relevant.
**What's needed:** Confirm no molecular benchmark data.

### 7. Motta et al., npj Quantum Inf. 7, 83 (2021) — already cited in Paper 14
**Title:** "Low rank representations for quantum simulation of electronic structure"
**Access:** DENIED (WebFetch blocked)
**What I recall:** Foundational DF paper. Reports rank (L) and lambda for small molecules including H2, H2O, N2 at cc-pVDZ. Key table compares lambda across methods.
**What was extracted:** General 2-5x lambda compression factor — LOW confidence for specific numbers.
**What's needed:** Exact lambda values from their benchmark tables. This is likely the BEST source for small-molecule DF lambda, as it explicitly benchmarks small systems.

### 8. Reiher et al., PNAS 114, 7555 (2017)
**Title:** "Elucidating Reaction Mechanisms on Quantum Computers"
**Access:** Not attempted (not in task list)
**What I recall:** Landmark FeMoCo QPE paper. ~10^14 T-gates baseline. Established resource estimation methodology.
**Relevance:** Historical context only. No small-molecule data.

---

## Papers NOT Searched (potential additional sources)

### Should search if internet access is granted:

1. **Rocca et al., arXiv:2404.xxxxx (2024)** — Recent compressed Gaussian benchmarks
2. **Oumarou et al., JCP 2024** — Linear T-complexity for electronic Hamiltonians
3. **Koridon et al., PRR 2024** — Orbital optimization for quantum chemistry
4. **Campbell et al., recent** — Any updated resource estimation surveys
5. **Google Quantum AI resource estimation papers (2023-2025)** — Often include small-molecule benchmarks

### The "small-molecule gap" problem:
The QPE resource estimation literature is biased toward large systems (50+ qubits) where fault-tolerant quantum advantage is plausible. At GeoVac's operating scale (10-70 qubits), the classical solution is trivially available, so the QPE community has little motivation to publish detailed resource estimates. This structural gap in the literature makes it genuinely difficult to find published compressed-Gaussian data for direct comparison.

---

## Data Confidence Summary

| Data Point | Confidence | Source | Verified? |
|-----------|------------|--------|-----------|
| GeoVac all numbers | HIGH | Paper 14, computed | Yes |
| Trenev et al. Pauli counts | HIGH | Published, in Paper 14 | Yes |
| GeoVac He vs Gaussian cc-pVDZ/pVTZ | HIGH | Computed integrals | Yes |
| Lee 2021 FeMoCo THC lambda ~306 Ha | MEDIUM | Training knowledge | No |
| Lee 2021 H2O THC lambda ~20-30 Ha | MEDIUM | Training knowledge | No |
| Goings 2022 P450 THC lambda ~1500 Ha | MEDIUM | Training knowledge | No |
| Motta 2021 DF 2-5x lambda reduction | MEDIUM | Training knowledge, widely cited | No |
| Loaiza 2024 SCDF 1.5-3x over DF | LOW | Training knowledge | No |
| Berry 2019 FeMoCo sparse lambda ~4800 | LOW | Training knowledge | No |
| All H2O DF/THC lambda estimates at cc-pVDZ | LOW | Extrapolation | No |
| All LiH DF/THC lambda estimates | LOW | Extrapolation, no published source | No |

---

## Action Items for PI

1. **Grant WebSearch/WebFetch permissions** and re-run this track to extract verified numbers from the papers listed above.
2. **Priority papers for manual lookup:**
   - Lee et al. 2021 Table I (H2O THC lambda) — most likely to have small-molecule data
   - Motta et al. 2021 benchmark tables (small-molecule DF lambda) — best source for direct comparison
   - Loaiza et al. 2024 tables (SCDF lambda for H chains, H2O)
3. **Compute GeoVac H2O 1-norm** — currently missing, prevents the most impactful comparison.
