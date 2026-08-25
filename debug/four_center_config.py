"""Four-center configuration operator F = sum_{i=1}^4 P_i: does a FOUR-nuclei molecule carry a
continuous four-body modulus (the D~4 cross-ratio) that the three-center case (discrete Z2 CIs)
lacked?  (2026-08-25, following the four-subspace<->elliptic-lambda lead.)

Honest scoping first:
 - REAL rank-1 (1s per center): the 4x4 Gram IS the configuration (up to O(4)), so EVERYTHING is
   pairwise -- no genuine four-body invariant, no cross-ratio. (Baseline, confirmed below.)
 - Genuine four-body content needs higher rank (sigma+pi {s,px,pz} per center = four rank-3
   subspaces in 12-dim) OR complexification (magnetic phase -> a continuous four-vertex phase).

Planar geometry (xz-plane), in-plane {s,px,pz} per center via the validated Slater-Koster block
(py decouples).  Charges all Z=1 (structural).  Invariants measured over a shape manifold:
  commutant dim (irreducible/wild?), pairwise principal angles, the genuinely-four-body nested
  commutator N4 = ||[[[P1,P2],P3],P4]||, the four-vertex Bargmann invariant Delta4 = Tr(P1P2P3P4)
  and its phase, and a reconstruction test (is Delta4 fixed by pairwise+triple data?).
"""
from __future__ import annotations
import sys, json
import numpy as np
sys.path.insert(0, "debug")
from fast_two_center_overlap import overlap_fast


def sk_block(R, u, Za=1, Zb=1):
    """lab-frame {s,px,pz}_A(0) vs {s,px,pz}_B(R u), u=(ux,uz) in-plane (xz)."""
    ss = overlap_fast(Za, 1, 0, Zb, 1, 0, 0, R)
    SP = overlap_fast(Za, 1, 0, Zb, 2, 1, 0, R)
    PS = overlap_fast(Za, 2, 1, Zb, 1, 0, 0, R)
    pps = overlap_fast(Za, 2, 1, Zb, 2, 1, 0, R)
    ppp = overlap_fast(Za, 2, 1, Zb, 2, 1, 1, R)
    Mloc = np.array([[ss, 0.0, SP], [0.0, ppp, 0.0], [PS, 0.0, pps]])
    ux, uz = u
    C = np.array([[1.0, 0.0, 0.0], [0.0, uz, ux], [0.0, -ux, uz]])
    return C @ Mloc @ C.T


def build_G(positions, rank=3):
    """positions: list of (x,z).  rank 3 = {s,px,pz}; rank 1 = {s} only (subblock)."""
    n = len(positions)
    B = np.eye(3)
    G = np.zeros((3 * n, 3 * n))
    for i in range(n):
        G[3*i:3*i+3, 3*i:3*i+3] = np.eye(3)
        for j in range(i + 1, n):
            d = np.array(positions[j]) - np.array(positions[i])
            R = np.hypot(d[0], d[1]); u = (d[0] / R, d[1] / R)
            Sij = sk_block(R, u)
            G[3*i:3*i+3, 3*j:3*j+3] = Sij
            G[3*j:3*j+3, 3*i:3*i+3] = Sij.T
    if rank == 1:
        idx = [3 * i for i in range(n)]                 # keep only the s row/col per center
        G = G[np.ix_(idx, idx)]
        bs = [1] * n
    else:
        bs = [3] * n
    return G, bs


def projectors(G, bs):
    if np.linalg.eigvalsh(G).min() < 1e-9:
        return None
    Xh = np.linalg.cholesky(G).T
    Ps, a = [], 0
    for s in bs:
        B = Xh[:, a:a + s]; a += s
        Ps.append(B @ np.linalg.pinv(B))
    return Ps


def commutant_dim(Ps, tol=1e-6):
    n = Ps[0].shape[0]
    rows = [np.kron(P, np.eye(n)) - np.kron(np.eye(n), P.T) for P in Ps]
    s = np.linalg.svd(np.vstack(rows), compute_uv=False)
    return int((s < tol * max(1.0, s[0])).sum()) + (n * n - len(s))


def cn(X):
    return float(np.linalg.norm(X, 2))


def nested4(Ps):
    a = Ps[0] @ Ps[1] - Ps[1] @ Ps[0]
    b = a @ Ps[2] - Ps[2] @ a
    c = b @ Ps[3] - Ps[3] @ b
    return cn(c)


def bargmann4(Ps):
    return complex(np.trace(Ps[0] @ Ps[1] @ Ps[2] @ Ps[3]))


def report(tag, positions, rank):
    G, bs = build_G(positions, rank)
    Ps = projectors(G, bs)
    if Ps is None:
        print("  %-40s: G not PSD" % tag); return None
    n = G.shape[0]
    dAp = commutant_dim(Ps)
    N4 = nested4(Ps)
    D4 = bargmann4(Ps)
    print("  %-40s: dim=%2d dim(A')=%d  N4=%.4f  Tr(P1P2P3P4)=%.4f (imag %.1e, arg=%.3f)"
          % (tag, n, dAp, N4, D4.real, D4.imag, np.angle(D4)))
    return dict(tag=tag, n=n, dimAp=dAp, N4=N4, D4_re=D4.real, D4_im=D4.imag)


if __name__ == "__main__":
    R = 2.5
    # a symmetric square (side ~R) as reference; then a shape sweep
    sq = [(-R/2, -R/2), (R/2, -R/2), (R/2, R/2), (-R/2, R/2)]     # 1-2-3-4 around the square
    print("=== BASELINE: four 1s (rank-1) -- expect fully pairwise, no four-body content ===")
    report("square, rank-1 {s}", sq, 1)
    # rank-1 four-body invariant should be a PRODUCT of pairwise overlaps (reconstruction check)
    G1, _ = build_G(sq, 1)
    Xh = np.linalg.cholesky(G1).T
    U = [Xh[:, i] / np.linalg.norm(Xh[:, i]) for i in range(4)]
    prod = (U[0] @ U[1]) * (U[1] @ U[2]) * (U[2] @ U[3]) * (U[3] @ U[0])
    Ps1 = projectors(G1, [1, 1, 1, 1])
    print("     rank-1 Tr(P1P2P3P4)=%.6f  vs product of pairwise <ui|uj>=%.6f  (equal => pairwise)"
          % (bargmann4(Ps1).real, prod))

    print("\n=== rank-3 {s,px,pz}: genuine four-body content? (square + rectangles) ===")
    for (w, h) in [(R, R), (R, 1.6*R), (R, 2.4*R), (1.4*R, 0.7*R)]:
        rect = [(-w/2, -h/2), (w/2, -h/2), (w/2, h/2), (-w/2, h/2)]
        report("rect w=%.1f h=%.1f" % (w, h), rect, 3)


# ---------------------------------------------------------------- continuous modulus + U(1)
def loop_holonomy(Ps):
    return complex(np.trace(Ps[0] @ Ps[1] @ Ps[2] @ Ps[3]))


def build_G_phase(positions, phase_bond=None, phi=0.0):
    """rank-3 G with an optional Peierls phase e^{i phi} on one off-diagonal bond block."""
    G, bs = build_G(positions, 3)
    G = G.astype(complex)
    if phase_bond is not None:
        i, j = phase_bond
        G[3*i:3*i+3, 3*j:3*j+3] *= np.exp(1j * phi)
        G[3*j:3*j+3, 3*i:3*i+3] *= np.exp(-1j * phi)
    return G, bs


def projectors_c(G, bs):
    if np.linalg.eigvalsh(G).min() < 1e-9:
        return None
    Xh = np.linalg.cholesky(G).conj().T
    Ps, a = [], 0
    for s in bs:
        B = Xh[:, a:a + s]; a += s
        Ps.append(B @ np.linalg.pinv(B))
    return Ps


def modulus_and_u1():
    R = 2.5
    print("\n=== continuous four-body MODULUS over shape (rhombus, half-diagonal ratio t) ===")
    print("  t     loop holonomy Tr(P1P2P3P4)   (a smooth continuous four-body modulus)")
    prev = None
    for t in np.linspace(0.5, 1.8, 14):
        rh = [(-R, 0.0), (0.0, -R*t), (R, 0.0), (0.0, R*t)]     # rhombus, diagonals 2R and 2Rt
        Ps = projectors(*build_G(rh, 3))
        h = loop_holonomy(Ps).real
        d = "" if prev is None else "  (d=%+.4f)" % (h - prev)
        print("  %.2f       %.5f%s" % (t, h, d)); prev = h

    print("\n=== U(1) lift: small magnetic flux on the 4-cycle -> holonomy phase continuous (Z2->U(1)) ===")
    sq = [(-R/2, -R/2), (R/2, -R/2), (R/2, R/2), (-R/2, R/2)]
    print("  eps     Tr(P1P2P3P4)          arg (Z2=0 at eps=0, continuous for eps!=0)   G_PSD")
    for eps in [0.0, 0.05, 0.10, 0.20, 0.35]:
        G, bs = build_G(sq, 3); G = G.astype(complex)
        for (i, j) in [(0, 1), (1, 2), (2, 3), (3, 0)]:            # flux circulation on the s-s cycle
            G[3*i, 3*j] *= np.exp(1j * eps); G[3*j, 3*i] *= np.exp(-1j * eps)
        evmin = np.linalg.eigvalsh(G).min()
        Ps = projectors_c(G, bs)
        if Ps is None:
            print("  %.2f    (G not PSD, evmin=%.1e)" % (eps, evmin)); continue
        h = loop_holonomy(Ps)
        print("  %.2f    %+.5f %+.5fj      arg=%+.5f   evmin=%.2e" % (eps, h.real, h.imag, np.angle(h), evmin))


if __name__ == "__main__":
    modulus_and_u1()
