"""Validation + basis-convergence for debug/p12_m_channel_probe.py.

Three controls with answers known independently of anything in this corpus:
  (V1) one H atom in the s set                 -> must approach -0.5 exactly
  (V2) H2 at R = 20 bohr                       -> must approach -1.0 (2 H atoms)
  (V3) H2 Hartree-Fock at R = 1.4011           -> must approach -1.1336 (HF limit)
and then the quantity of interest measured in a LARGER basis, to show the
sigma-only ceiling and the |m| >= 1 increment are converged as differences.

Run:  python debug/p12_m_channel_validate.py
"""

from __future__ import annotations

import time
import numpy as np

from geovac.noci_engine import (
    BasisFn, integral_set_md, lowdin_orbitals, transform_integrals, fci_ground,
)

R_BOHR = 1.4011
E_EXACT = -1.174475
DE_EXACT = 1.0 + E_EXACT * -1.0 - 1.0 + 0.174475   # = 0.174475, written explicitly below
DE_EXACT = -1.0 - E_EXACT


def shells(center, s_exp, p_exp, d_exp):
    orbs, mlab = [], []
    for a in s_exp:
        orbs.append(BasisFn(center, (0, 0, 0), np.array([a]), np.array([1.0])))
        mlab.append(0)
    for a in p_exp:
        for lmn, m in (((0, 0, 1), 0), ((1, 0, 0), 1), ((0, 1, 0), 1)):
            orbs.append(BasisFn(center, lmn, np.array([a]), np.array([1.0])))
            mlab.append(m)
    for a in d_exp:
        for lmn, m in (((0, 0, 2), 0), ((2, 0, 0), 9), ((0, 2, 0), 9),
                       ((1, 0, 1), 1), ((0, 1, 1), 1), ((1, 1, 0), 2)):
            orbs.append(BasisFn(center, lmn, np.array([a]), np.array([1.0])))
            mlab.append(m)
    return orbs, mlab


def m_transform(n, mlab):
    cols, mvals, i = [], [], 0
    while i < n:
        if mlab[i] == 9:
            v = np.zeros(n); v[i] = v[i + 1] = 1.0
            cols.append(v / np.sqrt(2)); mvals.append(0)
            v = np.zeros(n); v[i] = 1.0; v[i + 1] = -1.0
            cols.append(v / np.sqrt(2)); mvals.append(2)
            i += 2
        else:
            v = np.zeros(n); v[i] = 1.0
            cols.append(v); mvals.append(int(mlab[i]))
            i += 1
    return np.array(cols).T, np.array(mvals)


def fci_sub(C, s, h, g, n_elec=2):
    s2 = C.T @ s @ C
    h2 = C.T @ h @ C
    g2 = np.einsum("pi,qj,rk,sl,pqrs->ijkl", C, C, C, C, g, optimize=True)
    x = lowdin_orbitals(s2)
    ht, gt = transform_integrals(x, h2, g2)
    return fci_ground(ht, gt, n_elec)


def rhf_energy(s, h, g, n_elec=2, iters=60):
    """Closed-shell RHF in a non-orthogonal AO basis (n_elec = 2 -> 1 orbital)."""
    x = lowdin_orbitals(s)
    ht, gt = transform_integrals(x, h, g)
    n = ht.shape[0]
    d = np.zeros((n, n))
    e_old = 0.0
    for _ in range(iters):
        j = np.einsum("pqrs,rs->pq", gt, d, optimize=True)
        k = np.einsum("prqs,rs->pq", gt, d, optimize=True)
        f = ht + 2.0 * j - k
        w, c = np.linalg.eigh(f)
        occ = c[:, : n_elec // 2]
        d = occ @ occ.T
        e = np.einsum("pq,pq->", 2.0 * ht + 2.0 * j - k, d, optimize=True)
        if abs(e - e_old) < 1e-12:
            break
        e_old = e
    return e


S6 = np.array([0.0625, 0.1875, 0.5625, 1.6875, 5.0625, 15.1875])
S8 = np.array([0.0347, 0.0925, 0.2469, 0.6584, 1.7557, 4.6819, 12.485, 33.293])
P2 = np.array([0.35, 1.05])
P3 = np.array([0.25, 0.75, 2.25])
D1 = np.array([0.90])
D2 = np.array([0.55, 1.60])


def h2_system(R, s_exp, p_exp, d_exp):
    centers = [np.array([0.0, 0.0, -R / 2]), np.array([0.0, 0.0, +R / 2])]
    orbs, mlab = [], []
    for c in centers:
        o, m = shells(c, s_exp, p_exp, d_exp)
        orbs += o; mlab += m
    nuclei = [(c, 1.0) for c in centers]
    s, h, g = integral_set_md(orbs, nuclei)
    return orbs, np.array(mlab), s, h, g, 1.0 / R


print("=== V1: single H atom (must approach -0.500000) ===")
for name, sx in (("6s", S6), ("8s", S8)):
    orbs, _ = shells(np.zeros(3), sx, [], [])
    s, h, g = integral_set_md(orbs, [(np.zeros(3), 1.0)])
    x = lowdin_orbitals(s)
    ht, _gt = transform_integrals(x, h, g)
    e = np.linalg.eigvalsh(ht)[0]
    print(f"  {name}: E(H) = {e:.6f}   err = {1e6*(e+0.5):+8.2f} uHa")

print("\n=== V2: H2 at R = 20 bohr (must approach -1.000000) ===")
orbs, mlab, s, h, g, enuc = h2_system(20.0, S6, P2, [])
C, mv = m_transform(len(orbs), list(mlab))
e = fci_sub(C, s, h, g) + enuc
print(f"  6s2p FCI: E = {e:.6f}   err = {1e6*(e+1.0):+8.2f} uHa")

print("\n=== V3: H2 RHF at R = 1.4011 (HF basis limit = -1.13363) ===")
for name, sx, px, dx in (("6s2p", S6, P2, []), ("8s3p1d", S8, P3, D1)):
    orbs, mlab, s, h, g, enuc = h2_system(R_BOHR, sx, px, dx)
    e = rhf_energy(s, h, g) + enuc
    print(f"  {name}: E_RHF = {e:.6f}   vs -1.13363 -> {1000*(e+1.13363):+7.3f} mHa")

print("\n=== V4: m-channel decomposition, larger basis ===")
for name, sx, px, dx in (("6s2p", S6, P2, []),
                         ("8s3p", S8, P3, []),
                         ("8s3p2d", S8, P3, D2)):
    t0 = time.time()
    orbs, mlab, s, h, g, enuc = h2_system(R_BOHR, sx, px, dx)
    C, mv = m_transform(len(orbs), list(mlab))
    res = {}
    for mmax in sorted(set(mv)):
        sel = np.where(mv <= mmax)[0]
        res[mmax] = fci_sub(C[:, sel], s, h, g) + enuc
    ms = sorted(res)
    line = f"  {name:<8s} ({len(orbs):2d} fns, {time.time()-t0:5.1f}s): "
    line += f"sigma {res[0]:.6f} ({100*(-1.0-res[0])/DE_EXACT:5.2f}%)"
    for m in ms[1:]:
        line += f" | <={m} {res[m]:.6f} ({100*(-1.0-res[m])/DE_EXACT:5.2f}%)"
    line += f" | d(|m|>=1) = {1000*(res[ms[0]]-res[ms[-1]]):6.2f} mHa"
    print(line)

print(f"\n  Paper 12 sigma-only : -1.161304 (92.45%)   gap to exact = 13.17 mHa")
