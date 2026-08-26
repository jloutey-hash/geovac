"""
Paper 60 -- validated 2-electron H2 CI in the two-center shared-scale Sturmian basis, and the
STANDARD block-encoding 1-norm lambda = sum|h| + sum|(pq|rs)| (the qubitization cost driver),
vs basis size and vs a Gaussian basis for the SAME molecule.

This is the unambiguous, validatable molecular resource number.  (The isoenergetic-secular-matrix
1-norm -- the atomic sublinear object of eq:sublinear -- does not transfer to molecules for the
structural reason in sec:manyelectron; the STANDARD lambda below is what actually sets qubitization
cost and it IS measurable.)

Validation gates:
  * kinetic:   H atom (single 1s, scale 1) => E = -1/2.
  * CI:        R -> infinity => E_tot -> -1.0 Ha (two H atoms); interior minimum near R~1.4.
Block-encoding lambda in the Loewdin-orthonormalized basis, vs nmax, vs even-tempered s-Gaussians.
Diagnostic only.
"""
import warnings
warnings.filterwarnings("ignore")
import sys
import numpy as np
from itertools import combinations
from scipy.linalg import eigh, sqrtm

sys.path.insert(0, "debug")
from sturmian_goscinskian_integrals import GoscinskianIntegrals


def build_raw(gi, nmax, zeta):
    """One-electron h, overlap S, ERI tensor for shared-scale (decay zeta) s-Sturmians on A,B."""
    orbs = [(n, 'A', zeta) for n in range(1, nmax + 1)] + [(n, 'B', zeta) for n in range(1, nmax + 1)]
    N = len(orbs)
    S = np.zeros((N, N))
    h = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            S[i, j] = S[j, i] = gi.overlap(orbs[i], orbs[j])
            hij = gi.kinetic(orbs[i], orbs[j], zeta) + gi.nuclear(orbs[i], orbs[j])
            h[i, j] = h[j, i] = hij
    eri = np.zeros((N, N, N, N))
    for i in range(N):
        for j in range(i, N):
            for k in range(N):
                for l in range(k, N):
                    v = gi.eri(orbs[i], orbs[j], orbs[k], orbs[l])
                    for (a, b) in ((i, j), (j, i)):
                        for (c, d) in ((k, l), (l, k)):
                            eri[a, b, c, d] = v
    return h, S, eri


def loewdin(h, S, eri):
    w, U = np.linalg.eigh(S)
    keep = w > 1e-9
    X = U[:, keep] @ np.diag(1.0 / np.sqrt(w[keep]))       # S^{-1/2} (canonical, drops null)
    hm = X.T @ h @ X
    em = np.einsum('ip,jq,kr,ls,ijkl->pqrs', X, X, X, X, eri, optimize=True)   # C[raw,mo]=X
    return hm, em


def fci2e(h, eri):
    """Ground-state energy of 2 electrons, spin-orbital determinant FCI (orthonormal MOs)."""
    norb = h.shape[0]
    nso = 2 * norb
    dets = list(combinations(range(nso), 2))
    K = len(dets)

    def sp(i):
        return i // 2

    def spin(i):
        return i % 2

    def h1(i, j):
        return h[sp(i), sp(j)] if spin(i) == spin(j) else 0.0

    def g2(i, j, k, l):
        coul = eri[sp(i), sp(k), sp(j), sp(l)] if spin(i) == spin(k) and spin(j) == spin(l) else 0.0
        exch = eri[sp(i), sp(l), sp(j), sp(k)] if spin(i) == spin(l) and spin(j) == spin(k) else 0.0
        return coul - exch

    H = np.zeros((K, K))
    for a, (i, j) in enumerate(dets):
        for b, (k, l) in enumerate(dets):
            oa, ob = {i, j}, {k, l}
            diff = oa ^ ob
            if len(diff) == 0:
                val = h1(i, i) + h1(j, j) + g2(i, j, i, j)
            elif len(diff) == 2:
                m = (oa - ob).pop(); p = (ob - oa).pop(); c = (oa & ob).pop()
                sgn = (-1) ** ([i, j].index(m) + [k, l].index(p))
                val = sgn * (h1(m, p) + g2(m, c, p, c))
            elif len(diff) == 4:
                m1, m2 = sorted(oa - ob); p1, p2 = sorted(ob - oa)
                sgn = (-1) ** ([i, j].index(m1) + [i, j].index(m2) + [k, l].index(p1) + [k, l].index(p2))
                val = sgn * g2(m1, m2, p1, p2)
            else:
                val = 0.0
            H[a, b] = val
    return np.linalg.eigvalsh(0.5 * (H + H.T))[0]


# ---- Gaussian reference (same H2): even-tempered s-Gaussians, standard MO qubitization lambda
def gaussian_lambda(nmax, R, a0=0.10, ratio=2.0):
    exps = a0 * ratio ** np.arange(nmax)
    cen = [0.0] * nmax + [R] * nmax
    ex = list(exps) + list(exps)
    N = 2 * nmax

    def S1(a, b, d):
        return (4 * a * b / (a + b) ** 2) ** 0.75 * np.exp(-a * b / (a + b) * d ** 2)
    Sg = np.array([[S1(ex[i], ex[j], abs(cen[i] - cen[j])) for j in range(N)] for i in range(N)])
    # two-electron Gaussian ERI (s): standard closed form via the [0]^0 integral
    from scipy.special import erf

    def prim_eri(a, A, b, B, c, C, d, D):
        p = a + b; q = c + d
        P = (a * A + b * B) / p; Q = (c * C + d * D) / q
        Kab = (np.pi / p) ** 1.5 * np.exp(-a * b / p * (A - B) ** 2)
        Kcd = (np.pi / q) ** 1.5 * np.exp(-c * d / q * (C - D) ** 2)
        alpha = p * q / (p + q); T = alpha * (P - Q) ** 2
        F0 = 1.0 if T < 1e-12 else 0.5 * np.sqrt(np.pi / T) * erf(np.sqrt(T))
        return Kab * Kcd * 2.0 * np.sqrt(alpha / np.pi) * F0
    nrm = [(2 * ex[i] / np.pi) ** 0.75 for i in range(N)]
    er = np.zeros((N, N, N, N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    er[i, j, k, l] = (nrm[i] * nrm[j] * nrm[k] * nrm[l]
                                      * prim_eri(ex[i], cen[i], ex[j], cen[j], ex[k], cen[k], ex[l], cen[l]))
    # Loewdin orthonormalize (ERIs only; lambda_ee dominates) and sum|.|
    w, U = np.linalg.eigh(Sg); keep = w > 1e-9
    X = U[:, keep] @ np.diag(1.0 / np.sqrt(w[keep]))
    em = np.einsum('ip,jq,kr,ls,ijkl->pqrs', X, X, X, X, er, optimize=True)
    return np.abs(em).sum()


if __name__ == "__main__":
    np.set_printoptions(precision=5, suppress=True)

    print("=" * 74)
    print("VALIDATION")
    print("=" * 74)
    gfar = GoscinskianIntegrals(R=200.0, Lmax=8, nr=4000, nth=80, rmax=70.0)
    T = gfar.kinetic((1, 'A', 1.0), (1, 'A', 1.0), 1.0)
    V = -gfar.coulomb_center((1, 'A', 1.0), (1, 'A', 1.0), 'A')
    print(f"  H atom (1s, scale 1): <T>={T:.5f} (exact 0.5)  <V_A>={V:.5f} (exact -1.0)  E={T+V:.5f} (-0.5)")

    print("\n  H2 CI energy validation (shared-scale zeta swept; E_tot = E_elec + 1/R):")
    print(f"  {'R':>5} {'zeta*':>6} | {'E_elec':>9} {'E_tot':>9}   (ref: R->inf -> -1.0; min ~ -1.13 s-only)")
    for R in (1.4, 2.0, 6.0, 20.0):
        gi = GoscinskianIntegrals(R=R, Lmax=20, nr=2600, nth=160, rmax=max(60.0, 2 * R + 30))
        best = (1e9, None)
        for zeta in (0.9, 1.0, 1.1, 1.2, 1.3):
            h, S, eri = build_raw(gi, 2, zeta)
            hm, em = loewdin(h, S, eri)
            E = fci2e(hm, em)
            if E < best[0]:
                best = (E, zeta)
        print(f"  {R:>5.1f} {best[1]:>6.2f} | {best[0]:>9.4f} {best[0] + 1.0 / R:>9.4f}")

    print("\n" + "=" * 74)
    print("STANDARD block-encoding 1-norm lambda = sum|h| + sum|(pq|rs)|  (R=1.4, zeta=1.2)")
    print("=" * 74)
    R, zeta = 1.4, 1.2
    print(f"  {'nmax':>4} {'n_orb':>5} | {'lambda_SW (Sturmian)':>20} | {'lambda_Gauss':>13} | {'Gauss/SW':>9}")
    print("  " + "-" * 60)
    Ns, lamSW, lamG = [], [], []
    for nmax in (1, 2, 3):
        gi = GoscinskianIntegrals(R=R, Lmax=20, nr=2200, nth=150, rmax=55.0)
        h, S, eri = build_raw(gi, nmax, zeta)
        hm, em = loewdin(h, S, eri)
        lam = np.abs(hm).sum() + np.abs(em).sum()
        lg = gaussian_lambda(nmax, R)
        Ns.append(2 * nmax); lamSW.append(lam); lamG.append(lg)
        print(f"  {nmax:>4} {2 * nmax:>5} | {lam:>20.4f} | {lg:>13.4f} | {lg / lam:>9.2f}")
    Ns = np.array(Ns, float)
    pS = np.polyfit(np.log(Ns), np.log(lamSW), 1)[0]
    pG = np.polyfit(np.log(Ns), np.log(lamG), 1)[0]
    print(f"\n  lambda_SW ~ n_orb^{pS:.2f}   lambda_Gauss ~ n_orb^{pG:.2f}")
    print("\nDONE.")
