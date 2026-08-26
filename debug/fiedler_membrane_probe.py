"""Does the Fiedler membrane MOVE with the nuclei? (the killer check)

PI direction 2026-08-22, from the field-line/magnet reading of the two-center
problem. Paper 8 sec:membrane claims:

    "The Fiedler vector zero-crossing is the discrete, combinatorial analog of
     Bader's bond critical point... The GeoVac graph gives the same object
     WITHOUT invoking the electron density: the membrane position is determined
     purely by the graph Laplacian spectrum."

with three sub-claims: homonuclear crossing at the midpoint by symmetry;
heteronuclear crossing shifted toward the lighter nucleus; membrane thinner for
strongly polar bonds.

WHY THIS IS WORTH CHECKING, AND WHY FIRST. This observation is SPECTRAL and
TOPOLOGICAL, not energetic -- it needs the shape of an eigenVECTOR, not
R-dependent eigenVALUES. So it might survive everything that killed the bond
sphere as an energy method (Papers 8-9: H and S share one SO(4) congruence, so
the generalized eigenvalues are R-INDEPENDENT and the construction cannot bind).
Topology does not need separability, which is why the analogous framework
(Bader/QTAIM) works for polyatomics where no separable coordinate system exists.

BUT the same R-blindness could infect the eigenvector. If the Fiedler vector
does not respond to R, its "bond critical point" cannot move with geometry and
the correspondence is decorative. That is the one question that decides whether
anything downstream is worth building, so it is asked first and alone.

CONSTRUCTION (Paper 8 sec:bond_sphere)
  p_R = 1/R;  cos gamma = (p0^2 - p_R^2)/(p0^2 + p_R^2)      [eq:cos_gamma]
  within-block edges : the single-atom lattice at each center
  cross-block edges  : w = D^(n)_{(l'm'),(lm)}(gamma(R))     [eq:cross_weight]
  D is block-diagonal in n by SO(4) representation theory (not by approximation).

PRE-REGISTERED GATES
  V1 symmetry. For a HOMONUCLEAR graph the Fiedler vector must be antisymmetric
     under swapping the two blocks, so the amplitude fraction on each center is
     0.5 to machine precision. If this fails the assembly is wrong.
  V2 connectivity. lambda_2 > 0 (the graph must be connected -- otherwise the
     "Fiedler vector" is just an indicator of a disconnected piece and the
     zero-crossing is meaningless).

  DECISION (the killer check)
    MOVES     the membrane position varies with R by > 1% of its range over
              R in [1, 8] a0, monotonically in a physically sensible direction
    R-BLIND   position constant to < 1e-8 across all R  => STOP, the
              correspondence is decorative and nothing downstream is worth it
    NOISY     varies but non-monotonically / erratically => inspect

Exploratory. No paper claim.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.lattice import GeometricLattice
from geovac.wigner_so4 import bond_angle, d_matrix_block


def block_states(max_n):
    """(n,l,m) in the SAME order d_matrix_block uses: l ascending, m ascending."""
    out = []
    for n in range(1, max_n + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                out.append((n, l, m))
    return out


def cross_block(max_n, gamma):
    """Cross-center weight matrix: D^(n)(gamma), block-diagonal in n."""
    states = block_states(max_n)
    idx = {s: i for i, s in enumerate(states)}
    W = np.zeros((len(states), len(states)))
    for n in range(1, max_n + 1):
        D = d_matrix_block(n, gamma)
        sub = [(n, l, m) for l in range(n) for m in range(-l, l + 1)]
        for a, sa in enumerate(sub):
            for b, sb in enumerate(sub):
                W[idx[sa], idx[sb]] = D[a, b]
    return W


def within_block(max_n, Z):
    """Single-atom lattice adjacency, reordered onto block_states."""
    lat = GeometricLattice(max_n, nuclear_charge=Z)
    A = lat.adjacency.toarray()
    order = [lat.states.index(s) for s in block_states(max_n)]
    return A[np.ix_(order, order)]


def molecular_laplacian(max_n, Z_A, Z_B, R, p0):
    gamma = bond_angle(R, p0)
    WA = within_block(max_n, Z_A)
    WB = within_block(max_n, Z_B)
    X = cross_block(max_n, gamma)
    # nuclear-charge modulation of the cross block (Paper 8: "modulated by
    # nuclear charges"); symmetric so the graph stays undirected
    X = X * np.sqrt(Z_A * Z_B)
    W = np.block([[WA, X], [X.T, WB]])
    W = np.abs(W)                      # graph weights are magnitudes
    np.fill_diagonal(W, 0.0)
    L = np.diag(W.sum(axis=1)) - W
    return L, W


def fiedler(L):
    w, v = np.linalg.eigh(L)
    order = np.argsort(w)
    return w[order[1]], v[:, order[1]], w[order[0]]


def membrane(max_n, Z_A, Z_B, R, p0=1.0):
    """Return (lambda2, frac_A, imbalance).

    frac_A  = share of Fiedler amplitude^2 sitting on center A. For a
              homonuclear graph symmetry forces 0.5; a shift away from 0.5 is
              the membrane moving off the midpoint.
    imbalance = signed measure of how the SIGN partition aligns with the
              center partition.
    """
    L, _W = molecular_laplacian(max_n, Z_A, Z_B, R, p0)
    lam2, v, lam1 = fiedler(L)
    nA = len(block_states(max_n))
    p = v ** 2
    frac_A = float(p[:nA].sum() / p.sum())
    imb = float(np.sign(v[:nA]).sum() - np.sign(v[nA:]).sum()) / (2 * nA)
    return lam2, frac_A, imb, lam1


def main():
    max_n = 3
    print("=== Does the Fiedler membrane move with R? ===")
    print(f"    max_n = {max_n}  ({len(block_states(max_n))} states per center)\n")

    # ---- V1 / V2 -------------------------------------------------------
    print("[V1/V2] homonuclear symmetry + connectivity (Z_A = Z_B = 1)")
    ok1 = ok2 = True
    for R in (1.0, 2.0, 4.0):
        lam2, fA, imb, lam1 = membrane(max_n, 1.0, 1.0, R)
        sym = abs(fA - 0.5)
        ok1 = ok1 and sym < 1e-8
        ok2 = ok2 and lam2 > 1e-10
        print(f"    R={R:<4} frac_A={fA:.10f}  |dev|={sym:.2e}  "
              f"lam2={lam2:.6e}  lam1={lam1:.2e}")
    print(f"    V1 (frac_A == 0.5): {'PASS' if ok1 else 'FAIL'}")
    print(f"    V2 (connected):     {'PASS' if ok2 else 'FAIL'}")
    if not (ok1 and ok2):
        print("\n*** assembly invalid -- stopping ***")
        return

    # ---- the killer check ----------------------------------------------
    print("\n[KILLER CHECK] heteronuclear LiH-like (Z_A=3, Z_B=1) vs R")
    print(f"    {'R':>5} {'lambda2':>14} {'frac_A':>14} {'imbalance':>11}")
    rows = []
    for R in (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        lam2, fA, imb, _ = membrane(max_n, 3.0, 1.0, R)
        rows.append((R, lam2, fA, imb))
        print(f"    {R:5.1f} {lam2:14.8f} {fA:14.10f} {imb:11.4f}")

    fr = np.array([r[2] for r in rows])
    l2 = np.array([r[1] for r in rows])
    spread_f = fr.max() - fr.min()
    spread_l = l2.max() - l2.min()
    print(f"\n    frac_A  range over R: {spread_f:.3e}")
    print(f"    lambda2 range over R: {spread_l:.3e}")

    print("\n=== verdict ===")
    if spread_f < 1e-8 and spread_l < 1e-8:
        print("  R-BLIND: neither the membrane position nor lambda_2 responds")
        print("  to the nuclei moving. The Fiedler/Bader correspondence cannot")
        print("  track geometry -> decorative. STOP.")
    elif spread_f < 1e-8:
        print("  PARTIAL: lambda_2 moves but the membrane POSITION does not.")
        print("  The spectrum sees R; the partition does not. A bond critical")
        print("  point that cannot move is still not usable.")
    else:
        mono = np.all(np.diff(fr) > 0) or np.all(np.diff(fr) < 0)
        print(f"  MOVES: membrane position varies by {spread_f:.3e} over R,")
        print(f"  {'monotonically' if mono else 'NON-monotonically'}.")
        print("  -> worth taking to the three-centre (water) test.")

    # heteronuclear direction check, at fixed R
    print("\n[direction] frac_A vs charge asymmetry at R=3.0")
    for ZA in (1.0, 2.0, 3.0, 6.0):
        _, fA, _, _ = membrane(max_n, ZA, 1.0, 3.0)
        print(f"    Z_A={ZA:<4} Z_B=1   frac_A = {fA:.10f}")


if __name__ == "__main__":
    main()
