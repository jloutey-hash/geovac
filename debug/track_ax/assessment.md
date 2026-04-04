# Track AX: Competitive Assessment — GeoVac vs Compressed Gaussians

Generated: 2026-04-01

## Classification: CLEAR WIN on structural sparsity, with accuracy caveat

### Evidence

**He (equal-qubit, atoms):**
At Q=10, GeoVac has 120 Pauli terms vs Gaussian cc-pVDZ's 156 (1.3x). At Q=28, GeoVac
has 2,659 vs Gaussian cc-pVTZ's 21,607 (8.1x). The 1-norm ratios are even more favorable:
3.8x at Q=10, 6.8x at Q=28. This is a genuine comparison: same qubit count, same physical
system, independently computed integrals validated against published FCI energies.

**LiH (composed vs best available Gaussian):**
GeoVac at Q=30 produces 334 Pauli terms. The closest Gaussian (Trenev et al. cc-pVDZ at
Q=36, already including 2-qubit symmetry tapering) produces 63,519 terms. That is a 190x
advantage. At Q=84, the interpolated Gaussian count is ~2.4M terms vs GeoVac's 7,879 (307x).

**H2O (composed vs Gaussian interpolation):**
At equal qubit count Q=70, GeoVac produces 778 terms vs Gaussian-interpolated ~580,688 (746x).
At Q=196, the ratio reaches 1,712x.

**Scaling exponents:**
GeoVac composed: Q^2.50-2.52 (universal across LiH/BeH2/H2O, exponent spread 0.02).
Gaussian: Q^3.92-4.25 (system-dependent). The 1.5-2.0 exponent gap is structural, not
incidental. It arises from the block-diagonal ERI structure of composed geometries (cross-block
ERIs are exactly zero by construction). No amount of post-hoc compression changes the
scaling exponent.

### Accuracy Gap (must be stated honestly)

GeoVac composed geometries at l_max=2 produce:
- LiH: R_eq 5.3% error (equilibrium bond length)
- BeH2: R_eq 11.7% error
- H2O: R_eq 26% error

Gaussian baselines at cc-pVDZ or better produce <0.1% total energy error and sub-1%
geometry errors. This is not a close comparison on accuracy. GeoVac is not competitive
with Gaussians for classical potential energy surface computation.

The sparsity advantage is therefore most relevant for **quantum resource estimation**, where
the number of Pauli terms and the 1-norm determine circuit depth, Trotter step count, and
total gate count. A Hamiltonian with 190x fewer terms requires correspondingly fewer
measurements and shorter circuits, even if the Hamiltonian itself represents a less accurate
approximation to the exact electronic structure.

### Compression resistance

The Trenev et al. data already includes Z2 symmetry tapering (2-qubit reduction). Additional
Gaussian compressions include:

- **Frozen core:** Saves 2 qubits for LiH. Reduces Pauli terms by ~20-30%. Does not close
  a 190x gap.
- **Active space truncation:** Problem-dependent. Can dramatically reduce qubit count but
  introduces its own accuracy penalty.
- **Double factorization / tensor hypercontraction:** Can reduce Pauli terms by 2-10x for
  specific systems. Even a 10x reduction on cc-pVDZ LiH (63,519 -> 6,352) still leaves
  GeoVac 19x sparser.
- **Constant-factor compression cannot close a scaling exponent gap.** At Q=100+, the
  Q^2.5 vs Q^4.25 gap dominates any fixed-factor improvement.

### What this comparison does NOT show

1. **H2 Level 4:** No GeoVac qubit encoding exists for H2 in its natural coordinate system
   (mol-frame hyperspherical). The Level 4 solver achieves 96.0% D_e classically but has
   not been JW-encoded. This is the most favorable GeoVac system for classical accuracy,
   and its absence from the qubit comparison is a gap.

2. **GeoVac 1-norm for H2O:** Missing from Paper 14 (inflated by Z_eff=6 PK barrier).
   Without this, the quantum simulation cost comparison for H2O is incomplete.

3. **Gaussian 1-norms:** Trenev et al. report Pauli term counts but not 1-norms for their
   Gaussian encodings. The 1-norm comparison is therefore limited to He (computed integrals)
   and GeoVac composed systems. A full 1-norm comparison requires PySCF, which is not
   available on this platform.

### Verdict

**CLEAR WIN** on Pauli term count and scaling exponent. GeoVac composed geometries are
51x-1,712x sparser than the best available Gaussian baselines at equal or similar qubit
counts, and the scaling gap (Q^2.5 vs Q^3.9-4.3) means this advantage grows with system
size.

**CONDITIONAL** on accuracy: GeoVac's 5-26% geometry errors vs Gaussian's sub-1% mean
the sparser Hamiltonians represent a cruder approximation to the exact physics. The
value proposition is quantum resource efficiency (fewer gates, fewer measurements),
not classical accuracy parity. For production chemistry, the Gaussian basis is more
accurate; for near-term quantum hardware with limited circuit depth, GeoVac's structural
sparsity may be the enabling factor.
