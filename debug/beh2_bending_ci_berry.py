"""Robust two-band Z2 Berry phase for the BeH2 bent CIs + the sigma/pi-promotion mechanism
(Josh, 2026-08-25).  Isolates the crossing PAIR by energy-anchoring to the CI eigenvalue, so
transport can't hop to a neighbouring (unrelated) near-degeneracy -- the flakiness the global
min-gap tracker showed in the a=0 (sigma/pi-decoupled) limit.

A real-symmetric conical intersection has Z2 Berry phase = PI iff the loop encloses it; a
non-conical (glancing) touching gives 0.  So PI <=> genuine CI.
"""
from __future__ import annotations
import json
import sys
import numpy as np
sys.path.insert(0, "debug")
import beh2_bending_ci as BC
import beh2_bending_ci_symmetry as SY

SY._install_fast_funds()


def _F(d1, d2, a):
    F, evmin = SY.F_of(d1, d2, a)     # Lowdin
    return F, evmin


def berry_pair(dstar, a, rad, E_anchor=None, N=400):
    """Z2 phase of the LOWER of the two bands nearest E_anchor, transported around a loop in
    the (symmetric,asymmetric)-stretch plane (d1=dstar+s+t, d2=dstar+s-t) at fixed bend a."""
    if E_anchor is None:
        F0, _ = _F(dstar, dstar, a)
        w0 = np.linalg.eigvalsh(F0)
        k = int(np.argmin(np.diff(w0)))
        E_anchor = 0.5 * (w0[k] + w0[k + 1])
    v0 = vp = None
    for i in range(N + 1):
        th = 2 * np.pi * i / N
        s, t = rad * np.cos(th), rad * np.sin(th)
        F, evmin = _F(dstar + s + t, dstar + s - t, a)
        if F is None:
            return None
        w, V = np.linalg.eigh(F)
        # the crossing pair = the two levels nearest the anchor energy
        order = np.argsort(np.abs(w - E_anchor))
        pair = sorted(order[:2])
        lo = pair[0] if w[pair[0]] <= w[pair[1]] else pair[1]
        v = V[:, lo]
        if i == 0:
            v0 = v.copy(); vp = v
        else:
            v = v * np.sign(v @ vp); vp = v
    return float(np.sign(vp @ v0))


def ph(s):
    return "PI" if s is not None and s < 0 else ("0" if s is not None else "sing")


def main():
    out = {}
    print("=== robust Z2 Berry for the sigma-CI branch (connected to linear 2.445) ===")
    branch = [(180, 2.4451), (170, 2.4396), (160, 2.4219), (150, 2.3882),
              (140, 2.3313), (130, 2.2386), (120, 2.0889), (110, 2.0248)]
    rows = []
    for theta, dstar in branch:
        a = np.radians((180 - theta) / 2)
        F0, _ = _F(dstar, dstar, a); w0 = np.linalg.eigvalsh(F0)
        k = int(np.argmin(np.diff(w0))); Ea = 0.5 * (w0[k] + w0[k + 1])
        sgn = {r: berry_pair(dstar, a, r, Ea) for r in [0.04, 0.06, 0.09, 0.13]}
        rows.append(dict(theta=theta, dstar=dstar, E=float(Ea),
                         berry={r: ph(sgn[r]) for r in sgn}))
        print(f"  theta={theta:3d} d*={dstar:.4f} E={Ea:.4f}: " +
              "  ".join(f"r{r}:{ph(sgn[r])}" for r in [0.04, 0.06, 0.09, 0.13]))
    out["sigma_branch"] = rows

    print("\n=== the ~2.67 second CI branch ===")
    br2 = [(180, 2.671), (170, 2.671), (160, 2.672), (150, 2.669), (140, 2.652),
           (130, 2.596), (120, 2.438)]
    rows2 = []
    for theta, dstar in br2:
        a = np.radians((180 - theta) / 2)
        F0, _ = _F(dstar, dstar, a); w0 = np.linalg.eigvalsh(F0)
        k = int(np.argmin(np.diff(w0))); Ea = 0.5 * (w0[k] + w0[k + 1])
        sgn = {r: berry_pair(dstar, a, r, Ea) for r in [0.04, 0.06, 0.09]}
        rows2.append(dict(theta=theta, dstar=dstar, berry={r: ph(sgn[r]) for r in sgn}))
        print(f"  theta={theta:3d} d*={dstar:.4f}: " +
              "  ".join(f"r{r}:{ph(sgn[r])}" for r in [0.04, 0.06, 0.09]))
    out["branch_2p67"] = rows2

    print("\n=== the ~1.7 feature: NON-conical at theta=180 (Berry 0), PROMOTED to a real CI")
    print("    by the bend (sigma/pi mixing) -> Berry PI once bent ===")
    prom = [(180, 1.717), (175, 1.712), (170, 1.708), (165, 1.695), (160, 1.682),
            (150, 1.635)]
    rows3 = []
    for theta, dstar in prom:
        a = np.radians((180 - theta) / 2)
        F0, _ = _F(dstar, dstar, a); w0 = np.linalg.eigvalsh(F0)
        k = int(np.argmin(np.diff(w0))); Ea = 0.5 * (w0[k] + w0[k + 1])
        sgn = {r: berry_pair(dstar, a, r, Ea) for r in [0.03, 0.05, 0.07]}
        # sigma<->pi coupling strength at this geometry: max |G[sigma,pi]| off-block
        G = BC.metric(dstar, dstar, a)
        sig, pix = [0, 2, 3, 5, 6, 8], [1, 4, 7]
        cpl = float(np.max(np.abs(G[np.ix_(sig, pix)])))
        rows3.append(dict(theta=theta, dstar=dstar, sigpi_coupling=cpl,
                          berry={r: ph(sgn[r]) for r in sgn}))
        print(f"  theta={theta:3d} d*={dstar:.4f}  sig/pi-coupling={cpl:.3f}: " +
              "  ".join(f"r{r}:{ph(sgn[r])}" for r in [0.03, 0.05, 0.07]))
    out["promotion_1p7"] = rows3

    # sanity: far-away loop control (should be 0)
    print("\n=== far controls (should be 0) ===")
    for theta, dstar, cx in [(160, 2.42, 3.2), (120, 2.09, 3.3)]:
        a = np.radians((180 - theta) / 2)
        F0, _ = _F(cx, cx, a); w0 = np.linalg.eigvalsh(F0)
        k = int(np.argmin(np.diff(w0))); Ea = 0.5 * (w0[k] + w0[k + 1])
        s = berry_pair(cx, a, 0.1, Ea)
        print(f"  theta={theta} loop@d={cx} (away from CI): {ph(s)}")

    json.dump(out, open("debug/data/beh2_bending_ci_berry.json", "w"), indent=2)
    print("\nwrote debug/data/beh2_bending_ci_berry.json")


if __name__ == "__main__":
    main()
