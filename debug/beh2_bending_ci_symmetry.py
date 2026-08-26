"""Parity-resolved CI census for BeH2 off the linear axis (Josh, 2026-08-25).

On the symmetric manifold (d1=d2, symmetric bend a) the H1<->H2 swap combined with the z->-z
reflection is a symmetry Sigma of the configuration operator F = P_Be+P_H1+P_H2:
  Sigma permutes P_H1<->P_H2 and fixes P_Be  =>  [F,Sigma]=0.
Sigma splits the 9-dim in-plane space into EVEN (5-dim) and ODD (4-dim) sectors.  A genuine
conical intersection (real-symmetric, PI Berry phase) is exactly an EVEN<->ODD crossing; a
same-sector closest approach is NOT a protected CI (generically avoided off the line, Berry 0).

This resolves the muddle in beh2_bending_ci.py where the GLOBAL min-gap hopped between the
sigma CI and an unprotected pi-sector touching.  Here we track even<->odd crossings only, and
census ALL of them vs bend angle -- cleanly answering "does the linear CI survive/connect, and
does bending open NEW ones?"

Uses the VALIDATED engine (metric/projectors) from beh2_bending_ci (validation: sigma-F eig
reproduces the mpmath linear driver to 1e-13; sigma<->pi decoupling exact at a=0).
"""
from __future__ import annotations
import json
import sys
import numpy as np
sys.path.insert(0, "debug")
import beh2_bending_ci as BC
from scipy.linalg import sqrtm
from scipy.interpolate import CubicSpline
from fast_two_center_overlap import funds_fast


def _install_fast_funds():
    """Precompute the 5 fundamentals for BeH & HH on a dense 1D R-grid and spline them, then
    OVERRIDE BC.funds so metric/projectors/Berry all run at interpolation speed (overlaps are
    the cost).  Cubic spline of these smooth STO overlaps is ~1e-9 accurate."""
    RG = np.linspace(0.8, 9.5, 600)
    tab = {"BeH": np.array([funds_fast(R, 2, 1) for R in RG]),
           "HH": np.array([funds_fast(R, 1, 1) for R in RG])}
    sp = {p: [CubicSpline(RG, tab[p][:, i]) for i in range(5)] for p in tab}

    def fast_funds(R, pair):
        s = sp[pair]
        return (float(s[0](R)), float(s[1](R)), float(s[2](R)), float(s[3](R)), float(s[4](R)))
    BC.funds = fast_funds
    return RG.min(), RG.max()


def projectors_lowdin(d1, d2, a):
    """LOWDIN (symmetric) whitening -> center-projectors.  Identical F-spectrum to BC's
    Cholesky version (whitening-independent), but SYMMETRIC whitening is Sigma-equivariant:
    Sigma G^{1/2} Sigma = G^{1/2} when Sigma G Sigma = G, so F commutes with the original-basis
    Sigma and parity projection is valid."""
    G = BC.metric(d1, d2, a)
    evmin = float(np.linalg.eigvalsh(G).min())
    if evmin < 1e-9:
        return None, evmin
    Gh = sqrtm(G).real                       # G = Gh^T Gh (Gh symmetric)
    Ps = [Gh[:, 3 * k:3 * k + 3] @ np.linalg.pinv(Gh[:, 3 * k:3 * k + 3]) for k in range(3)]
    return Ps, evmin


def swap_operator():
    """Sigma (9x9): Be-block reflection diag(1,1,-1) [s,px,pz]; H1<->H2 with pz-sign flip."""
    J = np.diag([1.0, 1.0, -1.0])              # s,px unchanged; pz odd under z->-z
    S = np.zeros((9, 9))
    S[0:3, 0:3] = J                            # Be -> Be (reflected)
    S[3:6, 6:9] = J                            # H1 -> H2
    S[6:9, 3:6] = J                            # H2 -> H1
    return S


SIGMA = swap_operator()
# even/odd orthonormal bases (Sigma is symmetric orthogonal, Sigma^2 = I)
_w, _V = np.linalg.eigh(SIGMA)
_BE = _V[:, _w > 0]                            # even subspace (5-dim)
_BO = _V[:, _w < 0]                            # odd  subspace (4-dim)


def F_of(d1, d2, a):
    Ps, evmin = projectors_lowdin(d1, d2, a)   # Lowdin: Sigma-equivariant
    if Ps is None:
        return None, evmin
    return Ps[0] + Ps[1] + Ps[2], evmin


def check_commute(d, a):
    F, _ = F_of(d, d, a)
    return float(np.linalg.norm(F @ SIGMA - SIGMA @ F, 2))


def eig_by_parity(d, a):
    """Return sorted [(eig, parity)] with parity in {+1 even, -1 odd}, block-diagonalized."""
    F, evmin = F_of(d, d, a)
    if F is None:
        return None
    Ee = np.linalg.eigvalsh(_BE.T @ F @ _BE)   # 5 even
    Eo = np.linalg.eigvalsh(_BO.T @ F @ _BO)   # 4 odd
    lev = [(float(e), +1) for e in Ee] + [(float(e), -1) for e in Eo]
    lev.sort()
    return lev


def protected_gap(d, a):
    """min gap between ADJACENT opposite-parity levels (=0 at a protected even<->odd CI)."""
    lev = eig_by_parity(d, a)
    if lev is None:
        return np.nan
    g = [abs(lev[i + 1][0] - lev[i][0]) for i in range(len(lev) - 1)
         if lev[i + 1][1] != lev[i][1]]
    return min(g) if g else np.nan


def all_even_odd_crossings(a, dspan=(1.6, 3.4), n=700, tol=8e-3):
    """Find every d where an even level crosses an odd level (protected CI) at fixed bend a.
    Detect sign changes of each even_i(d)-odd_j(d) tracked branch via the merged-spectrum
    adjacent opposite-parity gap dipping to a local min below tol, then refine."""
    ds = np.linspace(*dspan, n)
    pg = np.array([protected_gap(d, a) for d in ds])
    cross = []
    for i in range(1, len(ds) - 1):
        if np.isnan(pg[i]):
            continue
        if pg[i] < tol and pg[i] <= pg[i - 1] and pg[i] <= pg[i + 1]:
            # refine
            lo, hi = ds[i - 1], ds[i + 1]
            dd = np.linspace(lo, hi, 200)
            gg = np.array([protected_gap(d, a) for d in dd])
            j = int(np.nanargmin(gg))
            if gg[j] < tol:
                cross.append((float(dd[j]), float(gg[j])))
    # dedupe close roots
    out = []
    for d, g in sorted(cross):
        if not out or abs(d - out[-1][0]) > 0.02:
            out.append((d, g))
        elif g < out[-1][1]:
            out[-1] = (d, g)
    return out


def berry(dstar, a, rad, N=240):
    """Z2 Berry phase in the (symmetric,asymmetric)-stretch plane at fixed bend a."""
    return BC.berry_st_loop(dstar, a, rad, N)


def main():
    rlo, rhi = _install_fast_funds()
    print(f"installed fast interpolated funds (R in [{rlo:.1f},{rhi:.1f}])")
    print("=== [F,Sigma] on the symmetric line (must be ~0) ===")
    for d, adeg in [(2.445, 0), (2.6, 15), (2.4, 30), (2.0, 45)]:
        print(f"  d={d}, a={adeg}deg: ||[F,Sigma]|| = {check_commute(d, np.radians(adeg)):.2e}")

    print("\n=== sigma-CI at a=0: parity of the crossing pair (linear reference) ===")
    lev0 = eig_by_parity(2.445, 0.0)
    print("  eig(parity) at (d=2.445, linear):",
          [(round(e, 5), '+' if p > 0 else '-') for e, p in lev0])
    print(f"  protected (even<->odd) gap at 2.445: {protected_gap(2.445, 0.0):.2e}")

    print("\n=== full even<->odd CI census vs bend angle theta ===")
    census = []
    for adeg in range(0, 56, 5):
        a = np.radians(adeg)
        theta = 180 - 2 * adeg
        cr = all_even_odd_crossings(a)
        # Berry phase for each crossing (confirm PI)
        entries = []
        for dstar, g in cr:
            bs = {r: berry(dstar, a, r) for r in [0.06, 0.10, 0.15]}
            ph = {r: ("PI" if s is not None and s < 0 else "0" if s is not None else "s")
                  for r, s in bs.items()}
            entries.append(dict(dstar=round(dstar, 4), gap=g, berry=ph))
        census.append(dict(theta=theta, a_deg=adeg, n_ci=len(cr), cis=entries))
        tag = "  ".join(f"d*={e['dstar']:.3f}(g{e['gap']:.0e},{e['berry'][0.1]})"
                        for e in entries) or "(none in window)"
        print(f"  theta={theta:3d} (a={adeg:2d}): {len(cr)} protected CI:  {tag}")
    json.dump(dict(census=census), open("debug/data/beh2_bending_ci_symmetry.json", "w"),
              indent=2)
    print("\nwrote debug/data/beh2_bending_ci_symmetry.json")

    # trace the branch continuously connected to the LINEAR sigma-CI (d~2.445 @ theta=180)
    print("\n=== branch continuously connected to the linear sigma-CI (2.445 @ 180) ===")
    prev = 2.445
    for adeg in range(0, 56, 5):
        a = np.radians(adeg); theta = 180 - 2 * adeg
        cr = all_even_odd_crossings(a)
        if not cr:
            print(f"  theta={theta:3d}: branch LOST (no protected CI near {prev:.3f})"); continue
        d_near = min(cr, key=lambda x: abs(x[0] - prev))
        if abs(d_near[0] - prev) > 0.45:
            print(f"  theta={theta:3d}: nearest protected CI d*={d_near[0]:.3f} "
                  f"jumped >0.45 from {prev:.3f} (branch may end)")
        b = {r: berry(d_near[0], a, r) for r in [0.08, 0.12]}
        ph = "/".join("PI" if b[r] is not None and b[r] < 0 else "0" for r in [0.08, 0.12])
        print(f"  theta={theta:3d} (a={adeg:2d}): sigma-CI d*={d_near[0]:.4f}  "
              f"gap={d_near[1]:.1e}  Berry(0.08/0.12)={ph}")
        prev = d_near[0]


if __name__ == "__main__":
    main()
