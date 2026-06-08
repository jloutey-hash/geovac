# Sprint m_j-parity Z₂ rotated-basis — relativistic chemistry

**Date:** 2026-06-08
**Driver:** `debug/sprint_mj_parity_rotated_diagnostic.py`
**Data:** `debug/data/sprint_mj_parity_rotated.json`
**Log:** `debug/sprint_mj_parity_rotated_diagnostic_log.txt`
**Predecessors:** Sprint m_j-parity direct-basis (v3.94.0, NEGATIVE); Sprint CH κ-parity (v3.92.0, NEGATIVE).

---

## TL;DR

**Verdict: NEGATIVE-STRUCTURAL with sharper mechanism than v3.94.0.** The rotated-basis $P_{m_j}^{\rm rot}$ — the relativistic analog of the non-relativistic Hopf m_l → −m_l Z₂ that DOES work — also fails. Residuals 1.7×10⁻² to 6.8×10⁻² on LiH_rel / BeH_rel / CaH_rel.

The densification check passes BIT-EXACTLY (max residual 0.0 between the builder's qubit_op and the rebuild from dense h1 + symmetrized ERI), confirming the implementation is correct. The m_j Z₂ genuinely isn't a symmetry of the relativistic chemistry Hamiltonian, even in the symmetric/antisymmetric basis.

**Sharper structural mechanism (the load-bearing finding of this sprint):** under all-orbital m_j → −m_j flip, the four 3-j sign factors in the jj-coupled X_k angular coefficient product give:

$$(-1)^{j_a + k + j_c} \times (-1)^{j_b + k + j_d} = (-1)^{j_a + j_b + j_c + j_d}$$

For **half-integer j (relativistic)**, this is NOT forced to be +1. Example: quadruple $(j_a, j_b, j_c, j_d) = (3/2, 1/2, 1/2, 1/2)$ → sum = 3 (odd) → factor = −1. ERI elements with such j quadruples are NEGATED under the all-m_j flip, so the ERI tensor is NOT invariant. The rotation to (sym, antisym) basis therefore does NOT block-diagonalize H, and the Z-string on antisym qubits doesn't commute.

For **integer l (non-relativistic Hopf)**, the analogous calculation gives $(-1)^{l_a + l_b + l_c + l_d}$, which IS forced to be even (= +1) by the Gaunt parity selection rule: $(l_a + l_c + k)$ even AND $(l_b + l_d + k)$ even forces $l_a + l_b + l_c + l_d$ even. This is exactly the basis of the v3.89.0 ℓ-parity Z₂ tapering ("Gaunt-parity"). The non-rel Hopf m_l Z₂ inherits this protection.

The relativistic jj-coupling has NO analogous parity rule forcing half-integer-j sum to be even. **m_j-parity Z₂ is structurally absent in relativistic chemistry, in both original and rotated bases.**

## 1. What was tested

1. Refactored `geovac/composed_qubit_relativistic.build_composed_hamiltonian_relativistic` to expose dense h1 + sub_blocks in the result dict (backward-compatible additive change).
2. Densified eri_sparse to a Q⁴ tensor. Symmetrized for Hermiticity (mirrored the builder).
3. Built the m_j → −m_j (sym, antisym) rotation U from `orbital_table` enumerated by `enumerate_relativistic_orbital_table` from v3.92.0 κ-parity module.
4. Densification + rebuild check: built fermion_op from dense h1 + symmetrized ERI, applied JW, compared to the builder's qubit_op. **Bit-exact match (residual 0.0)** on all three systems — confirms implementation correctness.
5. Applied covariant rotation to h1, ERI.
6. Rebuilt fermion_op from rotated integrals + JW → rotated qubit op.
7. Built Z-string $P_{m_j}^{\rm rot} = \prod_{q: \text{parity}(q) = -1} Z_q$ over rotated antisym qubits.
8. Audited commutation.

## 2. Results

| System | Q | N_pauli_orig | N_pauli_rotated | global residual | per-block residuals | Verdict |
|:------:|:-:|:------------:|:---------------:|:---------------:|:-------------------:|:-------:|
| LiH_rel | 30 | 1413 | 1453 (+40) | **5.07 × 10⁻²** | 5.07e-2, 1.69e-2, 1.69e-2 | NEGATIVE |
| BeH_rel | 30 | 1413 | 1453 (+40) | **6.75 × 10⁻²** | 6.75e-2, 3.38e-2, 1.69e-2 | NEGATIVE |
| CaH_rel | 20 |  942 |  966 (+24) | **3.38 × 10⁻²** | 3.38e-2, 1.69e-2 | NEGATIVE |

The Pauli count INCREASES under rotation (~3% inflation), consistent with the ERI tensor not respecting the rotation cleanly — the cross-sym-antisym sector pickups generate new Pauli terms.

## 3. Comparison: why the non-rel Hopf Z₂ DOES work

The non-relativistic Hopf m_l → −m_l Z₂ works because of three properties that the relativistic case lacks:

| Property | Non-rel (Hopf m_l Z₂ works) | Relativistic (m_j Z₂ fails) |
|:--------|:----------------------------:|:----------------------------:|
| Angular quantum number | $l$ (INTEGER) | $j$ (HALF-INTEGER) |
| Per-orbital m → −m partner | (n, l, ±m_l) at same (n, l) | (n_fock, κ, ±two_m_j) at same (n_fock, κ) |
| ERI selection rule | Gaunt: $(l_a + l_c + k)$ even | jj-coupled: $(l_a + l_c + k)$ even (parity) AND m_j conservation |
| 4-orbital sum-parity invariance under all-m flip | $(-1)^{l_a + l_b + l_c + l_d}$ FORCED to be +1 by Gaunt parity | $(-1)^{j_a + j_b + j_c + j_d}$ NOT forced — sometimes +1, sometimes −1 |
| Rotated antisym sector closed under H | ✓ | ✗ |
| Z-string commutes | ✓ (Paper 29 Hopf Z₂, v3.52.0) | ✗ (this sprint) |

The deep structural reason: **integer-l Gaunt parity protects the non-rel Hopf Z₂; half-integer-j jj-coupling has no analogous protection.** This is a genuinely structural distinction between the non-relativistic and relativistic chemistry constructions.

## 4. Joint structural reading with κ-parity (v3.92.0) and m_j-direct (v3.94.0)

Three NEGATIVE-STRUCTURAL findings on relativistic-chemistry parity-style tapering this week:

| Sprint | Stabilizer | Mechanism | Verdict |
|:------:|:----------:|:---------:|:-------:|
| v3.92.0 κ-parity (direct) | $\prod_{q:\kappa < 0} Z_q$ | $\kappa$-branch is NOT conserved by Dirac-Coulomb | NEGATIVE |
| v3.94.0 m_j-parity (direct) | $\prod_{q:m_j < 0} Z_q$ | $M_J$ conservation is SUM not parity | NEGATIVE |
| v3.95.0 m_j-parity (rotated) | $\prod_{q:\text{parity}^{\rm rot} = -1} Z_q$ | Half-integer-j sum NOT forced to even by selection rules | NEGATIVE |

Combined verdict: **relativistic chemistry has NO direct Z₂ analog of the non-relativistic Hopf m_l Z₂ tapering.** The structural reason — half-integer-j angular coupling vs integer-l Gaunt parity — is fundamental, not a tunable parameter.

The relativistic Tier 2 chemistry-tapering frontier is now empirically exhausted at sprint scale via this systematic exploration: three candidate Z₂ stabilizers tested, three NEGATIVE-STRUCTURAL with crisp mechanisms identified. Further relativistic-symmetry exploitation would require either:

- Non-Z₂ stabilizer machinery (e.g., $j^2$ sector restriction analog of the v3.92.0 S² sector — "good quantum number" projection, gives variational-manifold reduction not qubit ΔQ).
- A different angular coupling scheme (e.g., LS-coupled instead of jj-coupled; would lose Dirac-spinor exactness).
- Genuinely new structural insight (TBD).

## 5. Production code change shipped

`geovac/composed_qubit_relativistic.py` line 605: result dict now includes `'h1'` (dense Q×Q one-body matrix), `'M'` (alias for Q), and `'sub_blocks'` (sub-block enumeration). Backward-compatible additive change.

This refactor enables future rotated-basis tapering probes on the relativistic builder (Kramers, time-reversal × parity, etc.) without re-writing the dense-extraction code.

## 6. What this sprint did NOT do

- Did NOT find a working relativistic Z₂ tapering — three structural NEGATIVEs in a row.
- Did NOT explore non-Z₂ relativistic symmetry exploitation (full $j$ sector restriction; Kramers-pair projection; etc.).
- Did NOT touch CLAUDE.md or CHANGELOG.md directly (this memo provides staging content).

## 7. Hard-prohibition check (CLAUDE.md §13.5)

No changes to natural geometry hierarchy, fitted/empirical parameters, §3 deletions (only adding the m_j-rotated NEGATIVE row), Paper 2 K = π(B+F−Δ) conjectural label.

## 8. Verification

- Densification + rebuild bit-exact match (residual 0.0) across all 3 systems — confirms code is correct.
- Rotation U is orthogonal (residual <1e-12).
- Rotated qubit op has predictable Pauli-count inflation (~3%) consistent with cross-sym-antisym pickups.
- Mechanism verified by hand on a (3/2, 1/2, 1/2, 1/2) j-quadruple example.

## 9. Sample size

- 3 relativistic systems × original + rotated bases × 1 cutoff (max_n=2)
- Sufficient for the structural verdict — the half-integer-j argument is independent of system or basis size.

## 10. Files

### Created
- `debug/sprint_mj_parity_rotated_diagnostic.py`
- `debug/data/sprint_mj_parity_rotated.json`
- `debug/sprint_mj_parity_rotated_diagnostic_log.txt`
- `debug/sprint_mj_parity_rotated_memo.md` (this file)

### Modified
- `geovac/composed_qubit_relativistic.py` — result dict additive extension (`h1`, `M`, `sub_blocks` keys).
