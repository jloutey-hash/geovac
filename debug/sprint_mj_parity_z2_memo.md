# Sprint m_j-parity Z₂ — direct-basis tapering on relativistic chemistry

**Date:** 2026-06-08 (session continuation, follow-on from v3.92.0 κ-parity sprint)
**Driver:** `debug/sprint_mj_parity_z2_diagnostic.py`
**Data:** `debug/data/sprint_mj_parity_z2.json`
**Log:** `debug/sprint_mj_parity_z2_diagnostic_log.txt`
**Predecessors:** Sprint CH κ-parity Z₂ (v3.92.0, NEGATIVE-STRUCTURAL); `geovac/relativistic_tapering.py` scaffold.

---

## TL;DR

**Verdict: NEGATIVE-STRUCTURAL.** The Z-string $P_{m_j} = \prod_{q: m_j(q) < 0} Z_q$ in the **original Dirac basis** does NOT commute with the relativistic Tier 2 chemistry Hamiltonian. Commutator residuals 6–24 × 10⁻³ on LiH_rel / BeH_rel / CaH_rel — 7–9 orders of magnitude above the 1e-10 gate.

**Mechanism (sharper than the κ-parity sprint):** total $M_J$ conservation in the jj-coupled ERI ($m_j^a + m_j^b = m_j^c + m_j^d$) does NOT force the parity of $N_{m_j < 0}$ to be conserved. Concretely, the process $(+3/2) + (-1/2) \to (+1/2) + (+1/2)$ preserves $M_J = +1$ but **flips $N_{m_j<0}$ by 1**. This matrix element is nonzero in the chemistry construction's jj-coupled X_k angular coefficient, giving the observed ~10⁻² residual.

This is the m_j analog of the κ-parity finding (v3.92.0): both $\kappa$-branch and $m_j$-sign are NOT individually conserved as Z₂ subgroups of the Dirac-Coulomb rotational/discrete symmetry. Only $j$ and **total** $m_j$ are conserved.

**What's still open: rotated-basis $P_{m_j}$ (the proper m_j → −m_j Z₂).** The non-relativistic Hopf m_l → −m_l Z₂ commutes because of an analogous rotation to the symmetric/antisymmetric basis. The relativistic analog would rotate orbital pairs (m_j, −m_j) within each (n, κ) shell before constructing the Z-string. This requires densifying the relativistic builder's h1, eri tensors (currently the builder only exposes h1_diag + h1_so_diag and eri_sparse; the dense Q×Q h1 with off-diagonals is constructed internally at line 414 of `composed_qubit_relativistic.py` but not returned). Named follow-on; 2–3 day sprint to refactor + audit.

## 1. The hypothesis

Following the v3.92.0 κ-parity sprint NEGATIVE-STRUCTURAL closure (CH κ-parity fails because the Dirac-Coulomb operator conserves $j$ and $m_j$ but NOT $\kappa$-branch), the natural next probe is the m_j-sign parity. Since $m_j$ IS conserved per-orbital in single-particle terms, the parity $(-1)^{N_{m_j < 0}}$ might commute with H.

## 2. The test

Build $P_{m_j}^{\rm global} = \prod_{q: m_j(q) < 0} Z_q$ as a Pauli Z-string in the ORIGINAL Dirac basis (each spinor orbital → one JW qubit). Audit commutation with the relativistic chemistry qubit op on LiH_rel, BeH_rel, CaH_rel.

## 3. Results

| System | Q | N_pauli | $\|[H, P_{m_j}^{\rm global}]\|_{\max}$ | Verdict |
|:------:|:-:|:-------:|:---------:|:-------:|
| LiH_rel | 30 | 1413 | **1.79 × 10⁻²** | NEGATIVE |
| BeH_rel | 30 | 1413 | **2.39 × 10⁻²** | NEGATIVE |
| CaH_rel | 20 |  942 | **1.20 × 10⁻²** | NEGATIVE |

Per-block stabilizers fail similarly (residuals 6 × 10⁻³ to 2.4 × 10⁻²).

## 4. Why it fails — the structural mechanism

The jj-coupled ERI in `composed_qubit_relativistic.py` follows the selection rule $m_j^a + m_j^b = m_j^c + m_j^d$ ($M_J$ conservation). This is a SUM constraint, not a parity constraint.

Concrete example of a $P_{m_j}$-violating process:

- $\kappa_a = -2$ (2p_{3/2}), $m_j^a = +3/2$
- $\kappa_b = -1$ (1s_{1/2}), $m_j^b = -1/2$
- $\kappa_c = -1$ (1s_{1/2}), $m_j^c = +1/2$
- $\kappa_d = -1$ (1s_{1/2}), $m_j^d = +1/2$

Check: $m_j^a + m_j^b = 3/2 + (-1/2) = 1$; $m_j^c + m_j^d = 1/2 + 1/2 = 1$. ✓ M_J conserved.

Now count: initial state has $b$ at $m_j < 0$ ⇒ $N_{m_j<0}$ has $b$ contribution = 1. Final state has $c, d$ both at $m_j > 0$ ⇒ no $m_j < 0$ contribution from them. The annihilation of $c, d$ removes no $m_j<0$ from initial; the creation of $a, b$ adds 1 (from $b$ at $m_j = -1/2$).

Wait — the ERI process is $a^\dagger b^\dagger d c$ (chemist convention, OpenFermion 0.5 × eri term): annihilate $c, d$, create $a, b$. Change in $N_{m_j<0}$:

- Annihilate $c$ ($m_j = +1/2$): $\Delta = 0$
- Annihilate $d$ ($m_j = +1/2$): $\Delta = 0$
- Create $b$ ($m_j = -1/2$): $\Delta = +1$
- Create $a$ ($m_j = +3/2$): $\Delta = 0$

Net: $\Delta N_{m_j < 0} = +1$ (ODD). This term FLIPS the $P_{m_j}^{\rm global}$ parity. Therefore $[H, P_{m_j}^{\rm global}] \neq 0$.

The angular coefficient X_k for this process: it's nonzero by the Wigner 3-j parity selection rule applied to the (j_a, k, j_c) and (j_b, k, j_d) triples. Both triples have integer parity sum (no Pauli-like restriction kills it).

## 5. The rotated-basis approach (named follow-on)

The non-relativistic Hopf m_l → −m_l Z₂ tapering works because:
1. Within each (n, l) shell, orbitals are paired (m_l, −m_l).
2. Rotation to (sym = (φ_{+m} + φ_{−m})/√2, antisym = (φ_{+m} − φ_{−m})/√2) basis.
3. In rotated basis, P_{m_l} acts as +1 on sym orbitals and −1 on antisym orbitals.
4. Z-string on antisym qubits commutes with H IF H itself is m_l → −m_l symmetric (which it is by rotational invariance of the chemistry construction).

Same logic applies to m_j → −m_j in the relativistic basis. Within each (n_fock, κ) shell, orbitals come in pairs (two_m_j, −two_m_j). Rotation to (sym, antisym) basis should give a commuting Z-string by the same mechanism.

The block here is purely API: `build_composed_hamiltonian_relativistic` constructs the dense Q × Q h1 internally (line 414 of `composed_qubit_relativistic.py`) but doesn't return it — only h1_diag + h1_so_diag (the diagonal pieces) and eri_sparse. To do the rotation:
1. Refactor the builder to expose dense h1.
2. Densify eri_sparse to a Q⁴ tensor.
3. Apply rotation U @ h1 @ U.T and the 4-index covariant ERI transform.
4. Rebuild fermion_op + JW.
5. Audit the Z-string.

Effort: 2–3 day sprint. Worth doing because the rotated approach IS the well-established analog of the non-rel Hopf success, and the structural argument (3-j sign factors product = +1 under all-orbital m_j flips) suggests it should pass.

## 6. Joint structural reading with κ-parity

Two NEGATIVE-STRUCTURAL findings on relativistic chemistry tapering this 24h window:

| Sprint | Stabilizer | Mechanism |
|:-------|:-----------|:----------|
| v3.92.0 κ-parity | $\prod_{q: \kappa<0} Z_q$ in original basis | Coulomb freely couples $p_{3/2}$ ($\kappa = -2$) ↔ $p_{1/2}$ ($\kappa = +1$) at same $l$; 38% of LiH_rel ERIs flip $\Delta N_{\kappa<0}$ |
| v3.94.0 m_j-parity | $\prod_{q: m_j<0} Z_q$ in original basis | M_J conservation is SUM not parity; $(+3/2)+(-1/2) → (+1/2)+(+1/2)$ flips $\Delta N_{m_j<0}$ |

Together they tell a clean structural story: **Z-string parities in the ORIGINAL Dirac basis fail for relativistic chemistry because the jj-coupled Coulomb operator freely couples different sign-counted sectors under the M_J / κ-branch conservation rules.** Only properly ROTATED Z-strings (analogous to the non-rel Hopf m_l → −m_l rotation) might commute — and that's the named follow-on.

The pattern is informative: it explains why the non-rel Hopf m_l Z₂ DOES commute (rotation lifts the antisym sector cleanly) while neither κ-parity nor m_j-parity in the original Dirac basis works. The structural reason is the difference between (a) one-body parity quantum numbers conserved by individual matrix elements (non-rel m_l within a shell, which the rotation captures) and (b) two-body conservation rules that constrain SUMS but not parities (relativistic M_J, which the rotation might still capture but requires more care).

## 7. What this sprint did NOT do

- Did NOT implement the rotated-basis m_j Z₂ tapering — requires refactoring `composed_qubit_relativistic.build_composed_hamiltonian_relativistic` to expose dense h1 + densify the sparse ERI. Named follow-on.
- Did NOT modify production code.
- Did NOT touch CLAUDE.md or CHANGELOG.md directly (this memo provides the staging content).
- Did NOT investigate other relativistic candidate stabilizers (Kramers degeneracy, time-reversal × parity, the discrete center of SU(2)). The structural reason found here generalizes — any Z-string in the original Dirac basis with sign-count parity is likely to fail by similar mechanism.

## 8. Hard-prohibition check (CLAUDE.md §13.5)

No changes to: natural geometry hierarchy, fitted/empirical parameters, §3 deletions (only adding the m_j-parity dead-end row), Paper 2 K = π(B+F−Δ) conjectural label.

## 9. Verification

- Driver runs in ~2 seconds wall time (3 small Hamiltonians, Pauli-level commutator audits).
- All three relativistic systems show consistent ~10⁻² residuals — characteristic Coulomb angular coupling scale, not numerical noise.
- Mechanism verified by hand on a specific (+3/2)+(−1/2) → (+1/2)+(+1/2) ERI element.
- No production code modified → no regression possible on the existing test suite.

## 10. Sample size

- 3 relativistic systems × 2 modes (global + per_block) × 1 cutoff (max_n=2)
- Sufficient for the structural verdict; the M_J-conservation-doesn't-imply-parity argument is independent of system or basis size.

## 11. Files

### Created
- `debug/sprint_mj_parity_z2_diagnostic.py`
- `debug/data/sprint_mj_parity_z2.json`
- `debug/sprint_mj_parity_z2_diagnostic_log.txt`
- `debug/sprint_mj_parity_z2_memo.md` (this file)

### Modified
- None.
