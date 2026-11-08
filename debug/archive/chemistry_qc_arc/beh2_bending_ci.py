"""BeH2 three-center conical-intersection / Berry-phase analysis OFF the linear axis
(Josh, 2026-08-25).

The linear driver (beh2_three_way_invariant.py + beh2_ci_refine.py) found: three metric-
orthogonal center-projections P_Be, P_H1, P_H2 (sigma pair {1s,2p0} per center), configuration
operator F = sum P_i, a CONICAL INTERSECTION on the symmetric linear line at d* ~ 2.445 bohr
with a PI (Z2) Berry phase.

This driver extends that OFF the collinear axis -- the physical Renner-Teller bend.  When the
molecule bends the three centers are non-collinear, so p-orbital overlaps depend on the bond
ORIENTATION, not just the separation R.  Tractable, EXACT route = SLATER-KOSTER direction-
cosine decomposition: a rotation aligning the bond with z acts on the l=1 p-orbitals as the
l=1 rotation matrix (= direction cosines), so the two-center overlap of orbitals that factor
as (radial)x(real spherical harmonic) is reproduced exactly.

  Fundamentals (extracted from the axial engine at each separation R):
    ss  = <1s|1s>,  SP = <s_A|pz'_B>,  PS = <pz'_A|s_B>,  ppsig = <pz'|pz'>,  pppi = <px'|px'>.
  Bond-LOCAL 3x3 block (basis {s, px', pz'}, in-plane): M_local = [[ss,0,SP],[0,pppi,0],[PS,0,ppsig]].
  LAB block (basis {s,px,pz}): M_lab = C M_local C^T, C = [[1,0,0],[0,uz,ux],[0,-ux,uz]]
  (u = in-plane unit bond vector, uy=0).  p_y is out-of-plane and DECOUPLES for in-plane
  geometry, so we work in the 9-dim in-plane sector {s,px,pz} x 3 centers.

Overlaps: computed with the FAST float64 prolate-spheroidal engine (fast_two_center_overlap),
validated to <1e-9 vs the mpmath topos3 engine + the bargmann cache + analytic 1s-1s.  {1s,2p}
are single-exponential STOs so this is exact to ~1e-9 and runs in ms (the mpmath engine is too
slow for the diffuse Z=1 2p m=1 pi integral).

Geometry: C2 axis = x, molecule in xz-plane, so the LINEAR limit lies along z (matches the
linear driver's {1s,2p0=pz} basis exactly).  Symmetric half-angle a:
    Be at origin,  H1 = d1*(sin a,0, cos a),  H2 = d2*(sin a,0,-cos a),  theta(H-Be-H)=180-2a.
At a=0, d1=d2 the sigma subsector {s,pz}x3 reproduces the linear driver's 6-dim F EXACTLY
(the px orbitals decouple as a pi block); validated below against the mpmath linear driver.
"""
from __future__ import annotations
import json
import sys
import numpy as np
sys.path.insert(0, "debug")
from fast_two_center_overlap import funds_fast

_MEMO = {}


def funds(R, pair):
    Za, Zb = (2, 1) if pair == "BeH" else (1, 1)
    key = (round(float(R), 7), pair)
    if key not in _MEMO:
        _MEMO[key] = funds_fast(R, Za, Zb)
    return _MEMO[key]


# ---------------------------------------------------------------- SK oriented block

def C_of(u):
    """3x3 lab<-local map on {s,px,pz} for in-plane unit bond vector u=(ux,0,uz)."""
    ux, uz = u[0], u[2]
    return np.array([[1.0, 0.0, 0.0],
                     [0.0, uz, ux],
                     [0.0, -ux, uz]])


def block(R, u, pair):
    ss, SP, PS, pps, ppp = funds(R, pair)
    Mloc = np.array([[ss, 0.0, SP],
                     [0.0, ppp, 0.0],
                     [PS, 0.0, pps]])
    C = C_of(u)
    return C @ Mloc @ C.T


# ---------------------------------------------------------------- geometry -> projectors

def geometry(d1, d2, a):
    Be = np.zeros(3)
    H1 = d1 * np.array([np.sin(a), 0.0, np.cos(a)])
    H2 = d2 * np.array([np.sin(a), 0.0, -np.cos(a)])
    return Be, H1, H2


def metric(d1, d2, a):
    Be, H1, H2 = geometry(d1, d2, a)

    def bond(A, B):
        v = B - A
        R = float(np.linalg.norm(v))
        return R, v / R
    R1, u1 = bond(Be, H1)
    R2, u2 = bond(Be, H2)
    Rh, uh = bond(H1, H2)
    b1 = block(R1, u1, "BeH")
    b2 = block(R2, u2, "BeH")
    bh = block(Rh, uh, "HH")
    I = np.eye(3)
    return np.block([[I, b1, b2], [b1.T, I, bh], [b2.T, bh.T, I]])


def projectors(d1, d2, a):
    G = metric(d1, d2, a)
    evmin = float(np.linalg.eigvalsh(G).min())
    if evmin < 1e-9:
        return None, evmin
    Xh = np.linalg.cholesky(G).T
    Ps = [Xh[:, 3 * k:3 * k + 3] @ np.linalg.pinv(Xh[:, 3 * k:3 * k + 3]) for k in range(3)]
    return Ps, evmin


def F_spectrum(d1, d2, a):
    Ps, evmin = projectors(d1, d2, a)
    if Ps is None:
        return None, evmin
    return np.linalg.eigvalsh(Ps[0] + Ps[1] + Ps[2]), evmin


def min_gap(d1, d2, a):
    w, evmin = F_spectrum(d1, d2, a)
    return np.nan if w is None else float(np.min(np.diff(w)))


# ---------------------------------------------------------------- validation @ a=0

def validate_linear(mpmath_anchor=True):
    """(1) fast-engine sigma sector at a=0 == same construction with fast overlaps built as
    the linear driver does; (2) ANCHOR one d against the actual mpmath linear driver."""
    import beh2_three_way_invariant as B
    print("\n=== VALIDATION: bent SK engine at a=0 vs linear axial driver ===")
    ok = True
    # (2) mpmath anchor at a single d (slow ~1-2 min): reproduce the mpmath linear F.
    if mpmath_anchor:
        d = 2.445
        (PBe, PH1, PH2), ev0, _ = B.projectors(d, 2, 1)      # mpmath linear driver
        w_lin_mp = np.linalg.eigvalsh(PBe + PH1 + PH2)
        G9 = metric(d, d, 0.0)
        sig = [0, 2, 3, 5, 6, 8]                              # s,pz for Be,H1,H2
        Gs = G9[np.ix_(sig, sig)]
        Xh = np.linalg.cholesky(Gs).T
        Ps = [Xh[:, 2 * k:2 * k + 2] @ np.linalg.pinv(Xh[:, 2 * k:2 * k + 2]) for k in range(3)]
        w_bent = np.linalg.eigvalsh(Ps[0] + Ps[1] + Ps[2])
        derr = float(np.max(np.abs(np.sort(w_lin_mp) - np.sort(w_bent))))
        print(f"  [mpmath anchor d={d}] sigma-F eig: linear={np.round(w_lin_mp,6)}")
        print(f"                                     bent  ={np.round(w_bent,6)}")
        print(f"    |dEig| = {derr:.2e}  ({'OK' if derr < 1e-6 else 'FAIL'})")
        ok = ok and derr < 1e-6
    # sigma/pi decoupling at a=0: s/pz - px cross-overlaps must vanish
    G9 = metric(2.445, 2.445, 0.0)
    sig, px = [0, 2, 3, 5, 6, 8], [1, 4, 7]
    cpl = float(np.max(np.abs(G9[np.ix_(sig, px)])))
    print(f"  sigma<->pi decoupling at a=0: max|G[sig,px]| = {cpl:.2e} (must be ~0)")
    ok = ok and cpl < 1e-12
    print(f"  VALIDATION {'PASSED' if ok else 'FAILED'}")
    return ok


# ---------------------------------------------------------------- CI curve d*(theta)

def ci_on_symmetric_line(a, dspan=(1.7, 3.3), coarse=320, fine=400):
    ds = np.linspace(*dspan, coarse)
    g = np.array([min_gap(d, d, a) for d in ds])
    i = int(np.nanargmin(g))
    lo, hi = ds[max(i - 1, 0)], ds[min(i + 1, coarse - 1)]
    ds2 = np.linspace(lo, hi, fine)
    g2 = np.array([min_gap(d, d, a) for d in ds2])
    i2 = int(np.nanargmin(g2))
    return float(ds2[i2]), float(g2[i2])


# ---------------------------------------------------------------- Berry phase (Z2)

def berry_st_loop(dstar, a, rad, N=240):
    """Z2 Berry phase of the min-gap F-eigenvector around a loop in the (symmetric,
    asymmetric)-stretch plane (d1=dstar+s+t, d2=dstar+s-t) at fixed bend a."""
    v0 = vp = None
    for i in range(N + 1):
        th = 2 * np.pi * i / N
        s, t = rad * np.cos(th), rad * np.sin(th)
        Ps, evmin = projectors(dstar + s + t, dstar + s - t, a)
        if Ps is None:
            return None
        w, V = np.linalg.eigh(Ps[0] + Ps[1] + Ps[2])
        if i == 0:
            k = int(np.argmin(np.diff(w))); v = V[:, k]; v0 = v.copy(); vp = v
        else:
            ov = V.T @ vp; j = int(np.argmax(np.abs(ov)))
            v = V[:, j] * np.sign(ov[j]); vp = v
    return float(np.sign(vp @ v0))


# ---------------------------------------------------------------- main

def main():
    out = {}
    ok = validate_linear()
    out["validation_passed"] = bool(ok)
    if not ok:
        print("ABORT: validation failed; do not trust bent numbers.")
        json.dump(out, open("debug/data/beh2_bending_ci.json", "w"), indent=2)
        return

    print("\n=== CI curve on the symmetric line: d*(theta) as the molecule bends ===")
    adegs = list(range(0, 51, 5))
    curve = []
    for adeg in adegs:
        a = np.radians(adeg)
        dstar, gmin = ci_on_symmetric_line(a)
        theta = 180 - 2 * adeg
        g_lo = min_gap(dstar - 0.4, dstar - 0.4, a)
        g_hi = min_gap(dstar + 0.4, dstar + 0.4, a)
        curve.append(dict(theta=theta, a_deg=adeg, dstar=dstar, gmin=gmin,
                          gap_off_lo=g_lo, gap_off_hi=g_hi))
        print(f"  theta={theta:3d}deg (a={adeg:2d}): d*={dstar:.4f}  min-gap={gmin:.2e}   "
              f"(gap +-0.4 on line: {g_lo:.3f}/{g_hi:.3f})")
    out["ci_curve"] = curve

    print("\n=== Z2 Berry phase around each bent CI (loop in (d_sym,d_asym) at fixed bend) ===")
    berry = []
    for adeg in [0, 10, 20, 30, 40, 50]:
        a = np.radians(adeg)
        dstar, gmin = ci_on_symmetric_line(a)
        signs = {r: berry_st_loop(dstar, a, r) for r in [0.05, 0.08, 0.12, 0.18]}
        phase = {r: ("PI" if s is not None and s < 0 else "0" if s is not None else "sing")
                 for r, s in signs.items()}
        theta = 180 - 2 * adeg
        berry.append(dict(theta=theta, a_deg=adeg, dstar=dstar, phase=phase))
        print(f"  theta={theta:3d}deg: d*={dstar:.4f}  Berry(r) -> " +
              "  ".join(f"r{r}:{phase[r]}" for r in [0.05, 0.08, 0.12, 0.18]))
    out["berry_bent_ci"] = berry

    print("\n=== full 9-dim min-gap map over (d, theta) [symmetric]: new degeneracies? ===")
    dgrid = np.linspace(1.6, 3.4, 90)
    tgrid = np.linspace(80, 180, 90)
    M = np.full((len(tgrid), len(dgrid)), np.nan)
    Gmin = np.full_like(M, np.nan)
    for it, th in enumerate(tgrid):
        a = np.radians((180 - th) / 2)
        for idx, d in enumerate(dgrid):
            w, evmin = F_spectrum(d, d, a)
            Gmin[it, idx] = evmin
            M[it, idx] = np.nan if w is None else float(np.min(np.diff(w)))
    imin = np.unravel_index(np.nanargmin(M), M.shape)
    print(f"  global min-gap over map: {np.nanmin(M):.3e} at "
          f"theta={tgrid[imin[0]]:.1f}, d={dgrid[imin[1]]:.3f}")
    print(f"  metric stays PD over map: min eig(G) = {np.nanmin(Gmin):.2e}")
    print("  crossing valley (min-gap over d) at sampled theta:")
    for it in range(0, len(tgrid), 12):
        row = M[it]; j = int(np.nanargmin(row))
        print(f"    theta={tgrid[it]:6.1f}: min-gap={row[j]:.2e} at d={dgrid[j]:.3f}")
    out["gap_map"] = dict(dgrid=dgrid.tolist(), theta=tgrid.tolist(), gap=M.tolist(),
                          global_min=float(np.nanmin(M)),
                          global_min_at=[float(tgrid[imin[0]]), float(dgrid[imin[1]])])

    # a second, transverse Berry check right on the map's global minimum
    th0, d0 = tgrid[imin[0]], dgrid[imin[1]]
    a0 = np.radians((180 - th0) / 2)
    out["gapmap_min_berry"] = {r: ("PI" if (s := berry_st_loop(d0, a0, r)) is not None and s < 0
                                   else "0" if s is not None else "sing")
                               for r in [0.06, 0.1, 0.15]}
    print(f"  Berry around map-min (theta={th0:.1f}, d={d0:.3f}): {out['gapmap_min_berry']}")

    json.dump(out, open("debug/data/beh2_bending_ci.json", "w"), indent=2)
    print("\nwrote debug/data/beh2_bending_ci.json")


if __name__ == "__main__":
    main()
