"""Elliptic-basis pilot -- do the bond's CM singular moduli sit at the variational
optimum of the H2 two-electron basis?  (2026-08-23)

The two-center two-electron ERI is a genus-1 period on the Legendre/Gamma(2) family
(Paper 59); modulus m = 1 - c_min/c_max (fock_f12_genus_probe.py), CM fibers at
singular moduli disc-4 m=1/2, disc-8 m=3-2sqrt2.  The elliptic geometry lives on the
INTERACTION.  Pilot question: is it also a useful STATE-SPACE basis principle -- if we
set the two 1s exponents so the dominant ERI's modulus is a CM value, is that exponent
ratio variationally special (near the CI minimum), or generic (inert)?

Key structural point (found while building): the SINGLE-zeta H2 basis has both densities
at one scale (c1=c2, m=0) => GENUS 0.  The elliptic regime only switches on with >=2
distinct scales -- which is exactly what the certified single-exponent-per-center closed
form cannot represent.  So the pilot runs on the mixed-exponent NUMERIC engine
(geovac.sturmian_integrals, ~1e-4), on n=1 1s Slaters (exact per-orbital kinetic).

System: H2, Z=1 each, fixed R.  Basis: K 1s-Slaters per center, exponents {e_i} shared
by symmetry.  E_tot = E_elec(FCI) + 1/R.
"""
from __future__ import annotations
import os, sys, time
import numpy as np
from itertools import combinations
from scipy.linalg import eigh

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from geovac.sturmian_integrals import GoscinskianIntegrals


def mixed_kinetic(gi, oi, oj):
    """<i|-1/2 grad^2|j> for n=1 1s Slaters, per-orbital decay (exact)."""
    ai, aj = oi[2], oj[2]
    S = gi.overlap(oi, oj)
    Ti = oi[0] * ai * gi.coulomb_center(oi, oj, oi[1]) - 0.5 * ai ** 2 * S
    Tj = oj[0] * aj * gi.coulomb_center(oi, oj, oj[1]) - 0.5 * aj ** 2 * S
    return 0.5 * (Ti + Tj)


def build_mixed(gi, exps):
    orbs = [(1, 'A', e) for e in exps] + [(1, 'B', e) for e in exps]
    N = len(orbs)
    S = np.zeros((N, N)); h = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            S[i, j] = S[j, i] = gi.overlap(orbs[i], orbs[j])
            hij = mixed_kinetic(gi, orbs[i], orbs[j]) + gi.nuclear(orbs[i], orbs[j])
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
    X = U[:, keep] @ np.diag(1.0 / np.sqrt(w[keep]))
    hm = X.T @ h @ X
    em = np.einsum('ip,jq,kr,ls,ijkl->pqrs', X, X, X, X, eri, optimize=True)
    return hm, em


def fci2e(h, eri):
    norb = h.shape[0]; nso = 2 * norb
    dets = list(combinations(range(nso), 2)); K = len(dets)
    sp = lambda i: i // 2; spn = lambda i: i % 2
    h1 = lambda i, j: h[sp(i), sp(j)] if spn(i) == spn(j) else 0.0

    def g2(i, j, k, l):
        coul = eri[sp(i), sp(k), sp(j), sp(l)] if spn(i) == spn(k) and spn(j) == spn(l) else 0.0
        exch = eri[sp(i), sp(l), sp(j), sp(k)] if spn(i) == spn(l) and spn(j) == spn(k) else 0.0
        return coul - exch
    H = np.zeros((K, K))
    for a, (i, j) in enumerate(dets):
        for b, (k, l) in enumerate(dets):
            oa, ob = {i, j}, {k, l}; diff = oa ^ ob
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


def h2_energy(exps, R, gi):
    h, S, eri = build_mixed(gi, exps)
    hm, em = loewdin(h, S, eri)
    return fci2e(hm, em) + 1.0 / R


if __name__ == "__main__":
    R = 1.4
    gi = GoscinskianIntegrals(R, Lmax=24, nr=3000, nth=200, rmax=60.0)

    print(f"=== validation (R={R}) ===")
    t0 = time.time()
    e1 = h2_energy([1.0], R, gi)
    print(f"  single-zeta e=1.0 : E_tot = {e1:.6f} Ha   (single-zeta H2 ~ -1.09)   [{time.time()-t0:.1f}s]")
    t0 = time.time()
    e2 = h2_energy([1.2, 1.2 * 2 ** 0.5], R, gi)
    print(f"  2-zeta {{1.2, 1.2*sqrt2}} : E_tot = {e2:.6f} Ha   [{time.time()-t0:.1f}s]")
    print(f"  correlation gain 2z vs 1z(1.2): {e2 - h2_energy([1.2],R,gi):+.6f} Ha")
