"""2-electron FCI over an arbitrary STO basis using TwoCenterLM. Validated by reproducing
the s-only pilot H2 energy. Ready for s+p once the cross-center ERI oracle blesses that class."""
import os, sys
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from two_center_grid_lm import TwoCenterLM
from elliptic_basis_pilot import loewdin, fci2e   # real FCI (all m=0 here)


def build_ci(orbs, eng):
    eng.clear_moments()
    N = len(orbs)
    S = np.zeros((N, N)); h = np.zeros((N, N)); eri = np.zeros((N, N, N, N))
    for i in range(N):
        for j in range(i, N):
            S[i, j] = S[j, i] = eng.overlap(orbs[i], orbs[j])
            h[i, j] = h[j, i] = eng.h_core(orbs[i], orbs[j])
    for i in range(N):
        for j in range(i, N):
            for k in range(N):
                for l in range(k, N):
                    v = eng.eri(orbs[i], orbs[j], orbs[k], orbs[l])
                    for (a, b) in ((i, j), (j, i)):
                        for (c, d) in ((k, l), (l, k)):
                            eri[a, b, c, d] = v
    return h, S, eri


def h2_energy_lm(orbs, eng, R):
    h, S, eri = build_ci(orbs, eng)
    hm, em = loewdin(h, S, eri)
    return fci2e(hm, em) + 1.0 / R


if __name__ == "__main__":
    R = 1.4
    eng = TwoCenterLM(R, nr=1400, nu=32, nphi=32, rmax=50.0, Lmax=14)
    # s-only 2-zeta basis, same as the pilot [1.2, 1.2*sqrt2]
    z = [1.2, 1.2 * 2 ** 0.5]
    orbs = [(z[0], 0, 0, "A"), (z[1], 0, 0, "A"), (z[0], 0, 0, "B"), (z[1], 0, 0, "B")]
    E = h2_energy_lm(orbs, eng, R)
    print(f"s-only 2-zeta H2 (R={R}): new engine E_tot = {E:.6f}")
    print(f"  pilot reference (Goscinskian): -1.152351 (high grid)")
    print(f"  diff = {E - (-1.152351):+.2e}  -> {'CI ASSEMBLY VALIDATED' if abs(E+1.152351)<2e-3 else 'CHECK'}")
    # single-zeta sanity
    orbs1 = [(1.0, 0, 0, "A"), (1.0, 0, 0, "B")]
    print(f"single-zeta (e=1.0): E_tot = {h2_energy_lm(orbs1, eng, R):.6f}  (pilot -1.106619)")
