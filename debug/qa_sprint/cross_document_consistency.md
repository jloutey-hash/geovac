# Cross-Document Consistency Check

**Date:** 2026-03-15
**Status:** Complete

---

## Consistency Matrix

| Claim | README | Paper 13 | Paper 12 | Paper 1 | Code output | Status |
|:------|:------:|:--------:|:--------:|:-------:|:-----------:|:------:|
| He energy = -2.9052 Ha | -2.9052 | -2.9052 | -- | -- | -2.905173 | CONSISTENT |
| He error = 0.05% | 0.05% | 0.05% | -- | -- | 0.0499% | CONSISTENT |
| H2 Neumann D_e = 92.4% | 92.4% | -- | 92.4% | -- | 92.2% (27bf) / 92.4% (46bf) | MINOR* |
| omega_e = 4435 cm-1 | 4435 | 4435 | -- | -- | 4435 (Hylleraas) | CONSISTENT |
| nu_01 = 4157 cm-1 | -- | 4157 | -- | -- | 4157 (Hylleraas) | CONSISTENT |
| kappa = -1/16 | -1/16 | -- | -- | -- | -0.0625 | CONSISTENT |
| 18/18 symbolic proofs | 18/18 | -- | -- | -- | 18/18 | CONSISTENT |
| Berry phase k = 2.113 | -- | -- | -- | 2.113 | **0 (arg) / 1.0 (log)** | **INCONSISTENT** |

---

## Detailed Findings

### CONSISTENT: He Hyperspherical Energy

All sources agree on -2.9052 Ha (0.05% error). The solver produces -2.905173 Ha, which rounds to -2.9052 at 4 significant figures. The exact reference value -2.903724 Ha (Pekeris 1958) is consistently cited.

### CONSISTENT: kappa = -1/16

Hardcoded as `KINETIC_SCALE = -1/16` in `geovac/hamiltonian.py`. Validated by 18 symbolic proofs. No discrepancy across any document.

### CONSISTENT: Ab Initio Spectroscopy (Paper 13)

Paper 13's spectroscopic constants (omega_e = 4435, nu_01 = 4157) use Hylleraas+Neumann PES. The pipeline `run_pipeline.py` uses prolate CI PES (omega_e = 4918) -- a different method. The README matches Paper 13, not the pipeline. This is correct since Paper 13's method is more accurate.

**Note:** `benchmarks/ab_initio_nuclear/results.md` contains the pipeline numbers, not Paper 13's numbers. This file should be annotated.

### MINOR: H2 Neumann 92.4% vs 92.2%

The converged D_e percentage is 92.4% (at 46+ basis functions). At 27 bf it's 92.2%. The README header says "92.4% D_e" without specifying basis size. Paper 12 reports the convergence table correctly. No document explicitly claims "27 bf: 92.4%", so this is cosmetic.

### **INCONSISTENT: Berry Phase Exponent**

| Source | Claimed k | Actual computed k | Method |
|--------|:---------:|:-----------------:|--------|
| Paper 1 | 2.113 +/- 0.015 | 0 (arg of real product) | arg(product of CG coefficients) |
| berry_phase.py | 1.0 (exact) | 1.0 | log-holonomy of topological weights |
| SU(2)/SU(1,1) log-holonomy | -- | 0.54 | log-holonomy of CG coefficients |

Paper 1's claimed computation is impossible as described (arg of real positive product = 0). The figure is an explicit [Placeholder]. No reproducing code exists. See `debug/qa_sprint/berry_phase_reconciliation.md` for full analysis.

---

## Action Items

| Priority | Item | Status |
|:--------:|------|:------:|
| CRITICAL | Paper 1: Remove or correct Berry phase claims (k = 2.113 unvalidated) | PENDING |
| MINOR | `benchmarks/ab_initio_nuclear/results.md`: Annotate that Paper 13 uses Hylleraas PES | PENDING |
| COSMETIC | Verify "27 bf" is not associated with "92.4%" in any document | OK |
| -- | All other claims | VERIFIED |
