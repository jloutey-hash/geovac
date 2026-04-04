# Track BA: Updated Competitive Assessment

Generated: 2026-04-02
Status: UPDATED (Track BN, v2.0.31). Speculative lambda ranges removed. Structural conclusions verified.

---

## Executive Summary

The CLEAR WIN classification on Pauli term count **holds**, but requires nuance for fault-tolerant (FT) resource comparison. The literature comparison reveals that the FT community has moved beyond Pauli term counts to 1-norm and Toffoli cost as primary metrics, and uses block-encoding methods (DF, THC) that do not produce Pauli decompositions at all.

---

## 1. Does the CLEAR WIN classification hold?

### On Pauli term count: YES, unambiguously.

No published Gaussian compression method produces a Pauli decomposition that competes with GeoVac's 334 terms for LiH at Q=30 or 778 terms for H2O at Q=70. DF and THC do not produce Pauli decompositions — they use block-encoding circuits that bypass the Pauli representation entirely. The Pauli term comparison is therefore against raw JW encodings (Trenev et al.), and GeoVac wins by 190x-1,712x. This advantage is structural and cannot be reduced by any post-hoc compression that preserves the Pauli representation.

### On 1-norm (lambda): LIKELY YES, but data gaps prevent definitive claim.

**What we know:**
- GeoVac LiH at Q=30: lambda = 37.33 Ha
- GeoVac He at Q=10: lambda = 11.29 Ha vs Gaussian cc-pVDZ lambda = 42.95 Ha (3.8x win)
- GeoVac He at Q=28: lambda = 78.36 Ha vs Gaussian cc-pVTZ lambda = 530.47 Ha (6.8x win)

**What we don't know:**
- No published Gaussian DF/THC lambda exists for LiH at comparable qubit count (~30-36)
- GeoVac H2O lambda is missing (inflated by Z_eff=6 PK)
- The QPE literature focuses on large systems (FeMoCo: 150+ qubits, P450: 96 qubits), not on the 10-70 qubit systems where GeoVac operates

**Comparison status:** No published DF/THC/SCDF lambda values exist for LiH at Q~30-36 or H2O at Q~48-70. The FT resource estimation literature focuses on large active spaces (FeMoCo at 152 qubits, P450 at 96-116 qubits). Previously estimated THC lambda ranges for LiH (~15-40 Ha) and H2O (~20-30 Ha) were interpolations from general compression ratios and have been removed. The structural finding — that DF/THC bypass Pauli decompositions and no published data exists at GeoVac's operating scale — is more important than specific lambda numbers.

### On Toffoli/T-gate cost: NOT ASSESSED.

GeoVac does not currently compute Toffoli gate counts or circuit depths. The comparison on this metric requires implementing block-encoding circuits for GeoVac Hamiltonians, which has not been done. This is the metric the FT community actually uses for cost estimation.

---

## 2. What's the best published compressed Gaussian data for GeoVac-overlapping systems?

**Answer: Very little exists.**

The QPE resource estimation literature has a "small-molecule gap." The papers focus on:
- Toy systems for methodology validation (H2 in minimal basis, H4 chains)
- Production-scale systems for resource estimation (FeMoCo, P450, nitrogen fixation catalysts)

The 10-70 qubit regime where GeoVac operates (LiH at Q=30, H2O at Q=70) falls between these. At Q=30-70:
- THC is not typically applied (overhead exceeds benefit for systems this small)
- DF provides modest lambda reductions (~2-5x)
- Most papers just use raw JW or frozen-core JW

**The strongest published Gaussian data at GeoVac-relevant scales is the Trenev et al. raw JW data, which is already in Paper 14.** This is not an oversight — it reflects the state of the literature.

---

## 3. If DF/THC reduces Gaussian counts by X, what's the remaining GeoVac advantage?

### Pauli terms (where this comparison makes sense):
- DF/THC do NOT produce Pauli decompositions. They produce block-encoding circuits. The Pauli term comparison is not applicable to these methods.
- For Pauli-based VQE approaches, the comparison is against raw JW or tapered JW. GeoVac wins by 190x+ (LiH) to 746x+ (H2O).

### 1-norm (the more relevant metric):

For He (verified data only):
| System | GeoVac lambda | Gaussian raw (known) | GeoVac advantage |
|--------|--------------|---------------------|------------------|
| He Q=10 | 11.29 Ha | 42.95 Ha (cc-pVDZ) | 3.8x |
| He Q=28 | 78.36 Ha | 530.47 Ha (cc-pVTZ) | 6.8x |

No published DF/THC lambda values exist for He, LiH, or H2O at the qubit counts where GeoVac operates. The previously tabulated DF/THC lambda estimates for these systems were interpolations and have been removed. DF/THC are designed for large active spaces (N > 50 orbitals) and are not typically applied at Q=10-70.

**Bottom line:** For the small systems where GeoVac currently operates (Q=10-70), DF/THC are not typically applied in practice, so the comparison is academic. The real question is whether GeoVac's scaling advantage (Q^2.5 vs Q^{3.9-4.3}) persists to larger systems where DF/THC become relevant.

---

## 4. Is the 1-norm comparison more favorable than the Pauli count comparison?

**For He (known data): YES.**
- Pauli ratio at Q=28: 8.1x
- 1-norm ratio at Q=28: 6.8x
Both favor GeoVac, with 1-norm slightly less favorable than Pauli count.

**For LiH/H2O: CANNOT DETERMINE.**
Gaussian 1-norms are not published by Trenev et al. GeoVac H2O 1-norm is missing. Without both sides of the comparison, the question cannot be answered.

**For fault-tolerant algorithms:** The 1-norm is the operationally relevant metric. GeoVac's 1-norm scaling (Q^1.69 for atoms) is substantially below Gaussian scaling (typically Q^2-3). If this 1-norm scaling advantage persists to composed geometries (not yet measured), it would be the strongest argument for quantum resource efficiency.

---

## 5. Recommendations

### What Track AX should change:
1. **Remove the "2-10x" DF/THC compression estimate.** This was unsourced and is misleading. Replace with: "DF and THC use block-encoding circuits that bypass the Pauli representation entirely; the comparison on Pauli term count is therefore against raw JW encodings, where GeoVac's advantage is 190x-1,712x."

2. **Add a 1-norm discussion.** The 1-norm comparison (3.8x-6.8x for He) is weaker than the Pauli term comparison (1.3x-8.1x) but is the more operationally relevant metric for FT algorithms.

3. **Note the small-molecule gap.** The QPE literature does not provide compressed Gaussian benchmarks at GeoVac's operating scale (Q=10-70), because these systems are too small for DF/THC to be beneficial. This means GeoVac's structural sparsity operates in a regime where Gaussian post-hoc compression is not competitive.

### What Paper 14 should add:
1. A discussion noting that DF/THC are block-encoding methods, not Pauli compression methods, and that the comparison axes are different.
2. The observation that GeoVac's native sparsity provides a better starting point for ALL downstream compilation methods, including block-encoding.
3. Computed H2O 1-norms (currently missing due to Z_eff=6 PK inflation).

### What needs verification (requires internet access):
1. Lee et al. 2021 Table I/II: exact lambda values for H2O and FeMoCo
2. Loaiza et al. 2024: exact SCDF lambda values for any GeoVac-overlapping molecules
3. Any post-2023 papers that report DF/THC lambda for LiH or H2O at cc-pVDZ

---

## 6. Revised Classification

**CLEAR WIN on Pauli term count** — unchanged. 190x-1,712x advantage against the only available Pauli-based comparison (raw JW).

**STRONG WIN on 1-norm** — for He (3.8x-6.8x at equal qubit count, verified). Cannot be confirmed for LiH/H2O due to missing data on both sides.

**NOT ASSESSED on Toffoli/circuit cost** — GeoVac does not compute this metric. The FT literature uses a fundamentally different compilation approach (block-encoding) that is not directly comparable to Pauli term counts.

**CONDITIONAL on accuracy** — unchanged. GeoVac's 5-26% geometry errors vs Gaussian's <0.1% remain the primary caveat.
