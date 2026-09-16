"""Independent-route probe: how much H2 correlation energy lives in |m| >= 1?

Paper 12's prolate-spheroidal CI is sigma-only (every basis function is
phi-independent, m1 = m2 = 0) and plateaus at 92.4% of D_e, a 13.17 mHa gap
that the paper attributes to the electron-electron cusp requiring non-analytic
r12 terms.  `debug/sprint_tmr_method_memo.md` argues the gap is instead the
missing |m| >= 1 (pi, delta) channels.

This probe tests that in a completely different basis -- Cartesian Gaussians
through the corpus's own McMurchie-Davidson engine -- where the sigma/pi/delta
split is exact and the restriction is trivially imposed by orbital selection.
It touches none of the prolate-spheroidal machinery, so it is a genuinely
independent route (memory rule: feedback_independent_route_crosscheck).

Design
------
Same geometry as Paper 12: H2 at R = 1.4011 bohr, nuclei on z.
Uncontracted even-tempered s and p shells on each H, plus an optional d shell.

Cartesian components carry definite |m| about the bond axis z:
    s                       -> m = 0
    p: z                    -> m = 0        x, y            -> |m| = 1
    d: zz, xx+yy            -> m = 0        xz, yz          -> |m| = 1
       xy, xx-yy            -> |m| = 2
so a fixed linear transformation resolves the Cartesian basis by |m| exactly.

FCI is run in nested subspaces; the energy difference between the |m| = 0
subspace and the full space is the correlation energy that Paper 12's ansatz
structurally cannot reach.

Run:  python debug/p12_m_channel_probe.py
"""

from __future__ import annotations

import time
import numpy as np

from geovac.noci_engine import (
    BasisFn,
    integral_set_md,
    lowdin_orbitals,
    transform_integrals,
    fci_ground,
)

# --- Paper 12's geometry and reference values -------------------------------
R_BOHR = 1.4011
E_EXACT = -1.174475          # Kolos & Wolniewicz, as used by Paper 12
E_SEPARATED = -1.0           # two ground-state H atoms, non-relativistic BO
DE_EXACT = E_SEPARATED - E_EXACT          # 0.174475 Ha
E_PAPER12 = -1.161304        # Paper 12, N = 72 (j = 3, l = 3), Neumann V_ee

# --- basis ------------------------------------------------------------------
# Even-tempered uncontracted sets.  Kept modest so the dense FCI stays cheap;
# the quantity of interest is a DIFFERENCE between subspaces of the same basis,
# which converges far faster than either absolute energy.
S_EXPONENTS = np.array([0.0625, 0.1875, 0.5625, 1.6875, 5.0625, 15.1875])
P_EXPONENTS = np.array([0.35, 1.05])
D_EXPONENTS = np.array([0.90])


def build_basis(centers, use_d: bool):
    """Return (orbs, mlabels) with mlabels[i] = |m| of Cartesian component i."""
    orbs, mlab = [], []
    for c in centers:
        for a in S_EXPONENTS:
            orbs.append(BasisFn(c, (0, 0, 0), np.array([a]), np.array([1.0])))
            mlab.append(0)
        for a in P_EXPONENTS:
            for lmn, m in (((0, 0, 1), 0), ((1, 0, 0), 1), ((0, 1, 0), 1)):
                orbs.append(BasisFn(c, lmn, np.array([a]), np.array([1.0])))
                mlab.append(m)
        if use_d:
            for a in D_EXPONENTS:
                for lmn, m in (((0, 0, 2), 0), ((2, 0, 0), 9), ((0, 2, 0), 9),
                               ((1, 0, 1), 1), ((0, 1, 1), 1), ((1, 1, 0), 2)):
                    orbs.append(BasisFn(c, lmn, np.array([a]), np.array([1.0])))
                    mlab.append(m)   # 9 = "xx/yy, needs recombination"
    return orbs, np.array(mlab)


def m_resolved_transform(orbs, mlab, use_d: bool):
    """Fixed linear map from Cartesian components to |m|-pure combinations.

    Returns (C, mvals): C has one column per |m|-pure function.
    Only the xx/yy pair needs recombining: xx+yy is m = 0, xx-yy is |m| = 2.
    """
    n = len(orbs)
    cols, mvals = [], []
    i = 0
    while i < n:
        if mlab[i] == 9:
            # orbs[i] = xx, orbs[i+1] = yy on the same centre/exponent
            v = np.zeros(n); v[i] = 1.0; v[i + 1] = 1.0
            cols.append(v / np.sqrt(2.0)); mvals.append(0)
            v = np.zeros(n); v[i] = 1.0; v[i + 1] = -1.0
            cols.append(v / np.sqrt(2.0)); mvals.append(2)
            i += 2
        else:
            v = np.zeros(n); v[i] = 1.0
            cols.append(v); mvals.append(int(mlab[i]))
            i += 1
    return np.array(cols).T, np.array(mvals)


def fci_in_subspace(C, s, h, g, n_elec=2):
    """FCI energy (electronic) inside the column span of C."""
    s2 = C.T @ s @ C
    h2 = C.T @ h @ C
    g2 = np.einsum("pi,qj,rk,sl,pqrs->ijkl", C, C, C, C, g, optimize=True)
    x = lowdin_orbitals(s2)
    ht, gt = transform_integrals(x, h2, g2)
    return fci_ground(ht, gt, n_elec)


def report(label, e_elec, e_nuc, n_orb, n_det):
    e_tot = e_elec + e_nuc
    de = E_SEPARATED - e_tot
    pct = 100.0 * de / DE_EXACT
    print(f"  {label:<34s} {n_orb:3d} orb {n_det:5d} det   "
          f"E = {e_tot:12.6f}   D_e% = {pct:6.2f}")
    return e_tot


def main():
    centers = [np.array([0.0, 0.0, -R_BOHR / 2.0]),
               np.array([0.0, 0.0, +R_BOHR / 2.0])]
    nuclei = [(c, 1.0) for c in centers]
    e_nuc = 1.0 / R_BOHR

    for use_d in (False, True):
        tag = "s + p + d" if use_d else "s + p"
        print(f"\n=== basis: {tag} (uncontracted, even-tempered) ===")
        orbs, mlab = build_basis(centers, use_d)
        t0 = time.time()
        s, h, g = integral_set_md(orbs, nuclei)
        print(f"  integrals: {len(orbs)} Cartesian fns in {time.time()-t0:.1f}s")

        C_all, mvals = m_resolved_transform(orbs, mlab, use_d)

        results = {}
        # nested subspaces by |m|
        for mmax in sorted(set(mvals)):
            sel = np.where(mvals <= mmax)[0]
            C = C_all[:, sel]
            n_sp = C.shape[1]
            n_det = n_sp * (2 * n_sp - 1)      # C(2*n_sp, 2)
            t0 = time.time()
            e = fci_in_subspace(C, s, h, g)
            lab = f"|m| <= {mmax}" + ("   (sigma only)" if mmax == 0 else "")
            results[mmax] = report(lab, e, e_nuc, n_sp, n_det)
            print(f"      ({time.time()-t0:.1f}s)")

        ms = sorted(results)
        print(f"\n  Paper 12 (sigma only, N=72)        "
              f"          E = {E_PAPER12:12.6f}   "
              f"D_e% = {100*(E_SEPARATED-E_PAPER12)/DE_EXACT:6.2f}")
        print(f"  exact (Kolos-Wolniewicz)            "
              f"         E = {E_EXACT:12.6f}   D_e% = 100.00")
        print("\n  --- channel decomposition (mHa) ---")
        for a, b in zip(ms[:-1], ms[1:]):
            print(f"  |m| = {b} channels add            "
                  f"{1000*(results[a]-results[b]):8.2f} mHa")
        print(f"  total |m| >= 1 contribution     "
              f"{1000*(results[ms[0]]-results[ms[-1]]):8.2f} mHa")
        print(f"  Paper 12's unexplained gap      "
              f"{1000*(E_PAPER12-E_EXACT):8.2f} mHa")


if __name__ == "__main__":
    main()
