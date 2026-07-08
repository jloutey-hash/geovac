"""
Qubit-encoding wall diagnostic: does the non-commuting/basic-construction structure
permit a SPARSITY-PRESERVING fermionic encoding, or does JW force Loewdin either
way?  (PI-directed, 2026-07-07.)

JW/second quantization needs canonical anticommutation {a_i, a_j^dag} = delta_ij,
i.e. an ORTHONORMAL orbital basis. The genuine two-center basis is non-orthogonal
(overlap S). The standard route Loewdin-orthogonalizes (X = S^{-1/2}); Track DF
Sprint 5 / Paper 8 BU-2 measured that this destroys Gaunt sparsity (14x Pauli /
2.8-4.5x 1-norm). Question: is the destruction avoidable given the m-block /
finite-principal-angle structure?

Two selection rules make up the angular (Gaunt) sparsity:
  * m-conservation (aligned two centers: azimuthal integral -> delta_{m m'})
  * l-triangle (within an m-block, pure-l orbitals couple by Gaunt/3j rules)
This probe separates what Loewdin does to each.

Uses the VALIDATED m=0 two-center overlap block (hydrogenic Z=1 bond block, R_true;
from the commutator probe). Pure numpy, instant.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import sqrtm, block_diag

np.set_printoptions(precision=3, suppress=True)

# m=0 cross-center overlap block <A_i|B_j>, states (1s,2s,2p0), l=(0,0,1), R_true=3.015
S_AB = np.array([[+0.3455, -0.2518, -0.4950],
                 [-0.2518, +0.7993, +0.1268],
                 [+0.4950, -0.1268, +0.4786]])
l_A = np.array([0, 0, 1])                         # l of (1s,2s,2p0)

# full m=0 overlap: basis [1s_A,2s_A,2p0_A, 1s_B,2s_B,2p0_B]; within-center = I
n = 3
S = np.block([[np.eye(n), S_AB], [S_AB.T, np.eye(n)]])
l_full = np.concatenate([l_A, l_A])               # l labels of the 6 orbitals

Sm12 = np.real(sqrtm(np.linalg.inv(S)))           # Loewdin orthogonalizer S^{-1/2}
assert np.allclose(Sm12 @ S @ Sm12, np.eye(2 * n), atol=1e-10)   # X^T S X = I

# ---- TEST 1: is m-selection Loewdin-invariant? ----
# Structural: S block-diagonal in m => S^{-1/2} block-diagonal in m (sqrtm of a
# block-diagonal matrix is block-diagonal). Demo with a 2-block toy (m=0 + a m=1).
S_m1 = np.array([[1.0, 0.30], [0.30, 1.0]])       # a nonzero m=1 cross-overlap toy
S_two_m = block_diag(S, S_m1)                     # m=0 (6) + m=1 (2), no cross-m
X_two_m = np.real(sqrtm(np.linalg.inv(S_two_m)))
cross_m_mass = np.linalg.norm(X_two_m[:6, 6:])    # entries connecting m=0 <-> m=1
print("TEST 1  m-selection under Loewdin")
print(f"  cross-m mass in S^-1/2 = {cross_m_mass:.2e}  => m-selection PRESERVED "
      f"(exactly block-diagonal in m)")

# ---- TEST 2: is within-m l-selection destroyed? ----
# each Loewdin orbital = column of S^{-1/2}; its l-composition:
l0 = np.where(l_full == 0)[0]
l1 = np.where(l_full == 1)[0]
print("\nTEST 2  within-m l-mixing (m=0 block)")
print("  Loewdin orbital : weight on l=0 / weight on l=1  (pure => 1.0/0.0 or 0.0/1.0)")
labels = ["1s_A", "2s_A", "2p0_A", "1s_B", "2s_B", "2p0_B"]
max_contam = 0.0
for j in range(2 * n):
    col = Sm12[:, j]
    w0 = np.sum(col[l0] ** 2) / np.sum(col ** 2)
    w1 = np.sum(col[l1] ** 2) / np.sum(col ** 2)
    nominal_l = l_full[j]
    contam = w1 if nominal_l == 0 else w0        # weight on the "wrong" l
    max_contam = max(max_contam, contam)
    print(f"  {labels[j]:6s} (nominal l={nominal_l}): {w0:.3f} / {w1:.3f}"
          f"   contamination = {contam:.3f}")
print(f"  max off-l contamination = {max_contam:.3f}  => l-selection DESTROYED "
      f"within the m-block")

# ---- TEST 3: does any l-preserving orthogonalization exist? ----
# a block-diagonal-in-l orthogonalizer would leave the l0<->l1 overlap; show the
# residual off-l overlap cannot be removed without mixing l.
# l0-l1 overlap present in S (cross-center <1s|2p0> etc.):
S_l0l1 = S[np.ix_(l0, l1)]
print("\nTEST 3  is there an l-preserving (block-diagonal-in-l) orthogonalization?")
print(f"  ||S[l0,l1]|| (overlap between l=0 and l=1 orbitals) = {np.linalg.norm(S_l0l1):.3f}"
      f"  (nonzero => cross-center l-mixing overlap)")
print(f"  Any X with X^T S X = I that is block-diagonal in l would leave this")
print(f"  overlap intact (block-diagonal X preserves the off-l block of S), so it")
print(f"  CANNOT orthogonalize. => no l-preserving orthogonalization exists.")

# ---- basic-construction check: reaching conjugate-pair (2x2) form is a dense rotation ----
U, sig, Vt = np.linalg.svd(S_AB)
print("\nTEST 4  finite principal-angle (basic-construction) form")
print(f"  cross-block singular values (cos theta_k) = {np.round(sig,3)}")
print(f"  the conjugate-pair basis diagonalizes the coupling into {len(sig)} 2x2")
print(f"  blocks, BUT the rotation U (from SVD) mixes l:")
# does U mix l? U (3x3) maps the A-orbitals (rows l_A=[0,0,1]) to conjugate-pair
# vectors (columns). A column mixes l if it has weight on both l=0 and l=1 rows.
l0A = np.where(l_A == 0)[0]
l1A = np.where(l_A == 1)[0]
col_contam = [min(np.sum(U[l0A, c] ** 2), np.sum(U[l1A, c] ** 2)) for c in range(n)]
print(f"    per-conjugate-pair-vector l-mixing (min(l0,l1 weight)) = "
      f"{np.round(col_contam,3)}")
print(f"    max = {max(col_contam):.3f} > 0 => reaching the clean 2x2 form is itself")
print(f"    a dense l-mixing rotation. No free lunch.")

print("\n" + "=" * 84)
print("VERDICT")
print("=" * 84)
print("  * m-selection (the azimuthal/Hopf sparsity) is Loewdin-INVARIANT -> transfers.")
print("  * within-m l-selection (l-triangle Gaunt sparsity) is DESTROYED by any")
print("    orthogonalization, because the two-center overlap is l-dense within m")
print("    (<1s|2p0> etc. nonzero). This is intrinsic, not a Loewdin artifact.")
print("  * the finite-angle basic-construction form is clean (2x2 blocks) but")
print("    reaching it is itself a dense l-mixing rotation -> no sparsity rescue.")
print("  => JW forces the l-mixing. Prize #2 transfers to the qubit product ONLY")
print("     for the m-sparsity; the l-sparsity requires a NON-ORTHOGONAL fermionic")
print("     encoding (keep S in the anticommutators) -- a research direction, not")
print("     a quick fix. The classical binding route does not, by itself, rescue")
print("     the qubit sparsity.")
