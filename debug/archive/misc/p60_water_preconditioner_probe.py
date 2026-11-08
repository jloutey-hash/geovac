"""Does the band-Toeplitz breach reach a polyatomic block with symmetry-INEQUIVALENT
centers?  This is the case that matters: Paper 60's gerade lever already fails here
(water A_1, cond ~ N^1.97), so if the preconditioner also fails the breach is
confined to homonuclear diatomics and is much less interesting.

Setup (Paper 60 sec:molecular): s-only shared-scale basis on O + 2H, every s-s
block a function of inter-center distance only.  C_2v's only action is H1 <-> H2
(O sits on the axis), so the totally-symmetric A_1 block -- which holds the ground
state -- contains BOTH the O functions and the symmetric H combination, and with
them the O<->H coupling the group cannot remove:

    A_1 = [[ I ,  sqrt2 C_OH ],
           [ sqrt2 C_OH ,  I + C_HH ]]

Symbol analysis (the prediction this probe tests).  A_1 is a 2x2 MATRIX-valued
symbol A(chi) = [[1, sqrt2 a_OH],[sqrt2 a_OH, 1 + a_HH]].  At chi = pi every
a -> j0(0) = 1 regardless of separation, so A(pi) = [[1, sqrt2],[sqrt2, 2]], which
is singular with null vector v = (sqrt2, -1)/sqrt3 and trace 3.  So:

  * the degeneracy is RANK ONE and its direction is a FIXED vector, the same for
    every geometry and every n;
  * rotating by Q = kron([v, v_perp], I_n) puts the whole zero in one n x n block;
  * so the matching preconditioner is blockdiag(tri(1,2,1), I) in that frame.

If that works, the breach is not a homonuclear accident.
"""
import numpy as np
from geovac.sturmian_sigma_law import sw_cross_block

K = 1.0
R_OH, R_HH = 1.809, 2.862          # water, ~104.5 deg
M_QUAD = 100_001


def tridiag(n):
    return (np.diag(2.0 * np.ones(n)) + np.diag(np.ones(n - 1), 1)
            + np.diag(np.ones(n - 1), -1))


def inv_sqrt(M):
    ev, U = np.linalg.eigh(M)
    assert ev.min() > 0, f"not PD: lam_min = {ev.min():.3e}"
    return U @ np.diag(ev ** -0.5) @ U.T


def water_A1(n):
    C_OH = sw_cross_block(K * R_OH, n, M=M_QUAD)
    C_HH = sw_cross_block(K * R_HH, n, M=M_QUAD)
    I = np.eye(n)
    return np.block([[I, np.sqrt(2) * C_OH],
                     [np.sqrt(2) * C_OH.T, I + C_HH]])


# the fixed null direction of the symbol at chi = pi
A_pi = np.array([[1.0, np.sqrt(2)], [np.sqrt(2), 2.0]])
w, V2 = np.linalg.eigh(A_pi)
print(f"symbol at chi=pi: eigenvalues {w.round(12)}  (one must be 0, trace 3)")
print(f"null direction v = {V2[:, 0].round(6)}   (fixed: independent of geometry and n)\n")
R2 = V2[:, [0, 1]]                       # columns: null direction, then the rest


def preconditioned(n):
    A = water_A1(n)
    Q = np.kron(R2, np.eye(n))           # rotate the 2x2 block space only
    At = Q.T @ A @ Q
    P = np.zeros_like(At)
    P[:n, :n] = tridiag(n)               # quadratic zero along v ...
    P[n:, n:] = np.eye(n)                # ... identity on the healthy direction
    Pis = inv_sqrt(P)
    return A, Pis @ At @ Pis


print(f"  {'N=2n':>6} {'cond(A_1) raw':>15} {'growth':>8} {'cond(G) preconditioned':>24}")
prev = None
rows = []
for n in (6, 12, 24, 48, 96):
    A, G = preconditioned(n)
    cA, cG = np.linalg.cond(A), np.linalg.cond(G)
    g = "" if prev is None else f"{cA/prev:8.2f}"
    print(f"  {2*n:6d} {cA:15.1f} {g:>8} {cG:24.4f}")
    prev = cA
    rows.append((2 * n, cA, cG))

expo = np.polyfit(np.log([r[0] for r in rows]), np.log([r[1] for r in rows]), 1)[0]
print(f"\n  raw exponent  cond(A_1) ~ N^{expo:.2f}   (Paper 60 reports N^1.97)")
print(f"  preconditioned: {rows[0][2]:.4f} -> {rows[-1][2]:.4f} over N={rows[0][0]}..{rows[-1][0]}")

print("\n  control -- does the NAIVE (unrotated) preconditioner work?")
for n in (24, 96):
    A = water_A1(n)
    P = np.zeros_like(A)
    P[:n, :n] = tridiag(n)
    P[n:, n:] = tridiag(n)
    Pis = inv_sqrt(P)
    print(f"    N={2*n:4d}  blockdiag(T,T) without the rotation: "
          f"cond = {np.linalg.cond(Pis @ A @ Pis):12.1f}")
