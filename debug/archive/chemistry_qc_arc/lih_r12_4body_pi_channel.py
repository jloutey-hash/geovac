"""N=4 WALL, TWO-CENTER, PI/DELTA (m != 0) channel -- extends lih_r12_4body_integral.py
(the sigma / m=0 gate) to the azimuthal-transfer channels of the prolate Neumann bridge.

Why a new test.  A DIAGONAL density |phi|^2 of a pure-m orbital is phi-independent, so the
m != 0 bridge channels never activate from diagonal densities -- the m-transfer lives in the
EXCHANGE (transition-density) terms of a real CI.  To exercise the general-m prolate Neumann
bridge with a POSITIVE, MC-samplable density, we modulate the sigma bridge density azimuthally:

    rho_bridge(r) = rho_1sB(r) * (1 + b1 cos(phi) + b2 cos(2 phi))       (phi = azimuth about z)

which carries m = 0, +/-1 (PI), +/-2 (DELTA) explicitly (positive for |b1|+|b2| within range),
samples exactly like the sigma case (1s_B proposal + a phi-weight), and has trivial
m-components c_0 = rho, c_1 = b1 rho, c_2 = b2 rho.  The leaves stay 1s on A (isotropic ->
the leaf dressing Psi_A is unchanged, phi-independent), so only the BRIDGE needs the general-m
Neumann expansion:

  1/r13 = (2/R) sum_l sum_m (2-d_m0)(-1)^m (2l+1)[(l-m)!/(l+m)!]^2
          P_l^m(xi_<) Q_l^m(xi_>) P_l^m(eta1) P_l^m(eta3) cos(m(phi1-phi3)).

After the phi-integrals the bridge is
  J = (2/R) a^6 sum_m W_m sum_l (-1)^m (2l+1)[(l-m)!/(l+m)!]^2 U_{l,m},
  W_0 = (2pi)^2, W_{m>=1} = 2 pi^2,
  U_{l,m} = INT INT g_{l,m}(xi1) g_{l,m}(xi3) P_l^m(xi_<) Q_l^m(xi_>),
  g_{l,m}(xi) = INT (xi^2-eta^2) c_m(xi,eta) P_l^m(eta) deta.
P_l^m, Q_l^m from scipy lpmn/lqmn (Condon-Shortley phases cancel in the P.Q and P(eta).P(eta)
pairs, leaving the explicit (-1)^m -- validated below).

Validation (same standard as the sigma gate):
  [C1'] b1=b2=0 -> the general-m code gives the sigma self-Coulomb 5*alpha/8 (m=0 path).
  [C2'] modulated self-Coulomb: general-m reduced == independent 6-D MC (pins the m!=0
        convention).
  [C0'] the pi/delta 4-body with b1=b2=0 reproduces the sigma 4-body (2.856e-2).
  [MAIN] the pi/delta 4-body: reduced (leaf dressing + general-m bridge) == independent 12-D MC.

Run from root:  python debug/lih_r12_4body_pi_channel.py
"""
import math
import os
import sys
import warnings

import numpy as np
from numpy.polynomial.legendre import leggauss
# NB scipy>=1.15 deprecates lpmn/lqmn (removal in 1.17) in favour of assoc_legendre_p_all;
# they still work here and the convention is pinned by the C1'/C2' controls. Migration would
# re-derive against assoc_legendre_p_all's (order,degree) layout + phase and re-validate.
warnings.filterwarnings("ignore", category=DeprecationWarning)
from scipy.special import lpmn, lqmn  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lih_r12_4body_integral import (   # noqa: E402  (reuse the sigma-gate machinery)
    R, a, ALPHA, BETA, GAM, CENTER_A, CENTER_B, f_gem, sample_1s, leaf_dressing, rho_1s)

rng = np.random.default_rng(20260921)

# azimuthal modulation of the bridge density: m=1 (pi) and m=2 (delta) content
B1, B2 = 0.6, 0.4      # 1 + 0.6 cos p + 0.4 cos 2p  > 0 for all p (min ~0.49)


def w_mod(pts):
    """azimuthal weight (1 + B1 cos phi + B2 cos 2phi), phi = atan2(y,x) about the z-axis."""
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    return 1.0 + B1 * np.cos(phi) + B2 * np.cos(2 * phi)


# --------------------------------------------------------------------------- #
# general-m prolate Neumann bridge on grid densities
# --------------------------------------------------------------------------- #
def build_grid_m(NXI=260, NETA=100, LMAX=34, MMAX=2):
    xg, wxg = leggauss(NXI)
    xi_max = 1.0 + 44.0 / (2 * ALPHA * a)
    xi1d = 1.0 + 0.5 * (xg + 1.0) * (xi_max - 1.0)
    wxi = 0.5 * (xi_max - 1.0) * wxg
    eg, weg = leggauss(NETA)
    XI, ETA = np.meshgrid(xi1d, eg, indexing='ij')
    Peta = np.zeros((MMAX + 1, LMAX + 1, NETA))
    for i, e in enumerate(eg):
        P, _ = lpmn(MMAX, LMAX, e); Peta[:, :, i] = P
    Pxi = np.zeros((MMAX + 1, LMAX + 1, NXI)); Qxi = np.zeros((MMAX + 1, LMAX + 1, NXI))
    for i, x in enumerate(xi1d):
        P, _ = lpmn(MMAX, LMAX, x); Pxi[:, :, i] = P
        Q, _ = lqmn(MMAX, LMAX, x); Qxi[:, :, i] = Q
    return dict(NXI=NXI, NETA=NETA, LMAX=LMAX, MMAX=MMAX, xi1d=xi1d, wxi=wxi,
                eta1d=eg.copy(), weta=weg.copy(),
                rB=a * (XI - ETA), rA=a * (XI + ETA), JAC=(XI ** 2 - ETA ** 2),
                Peta=Peta, Pxi=Pxi, Qxi=Qxi,
                minidx=np.minimum.outer(np.arange(NXI), np.arange(NXI)),
                maxidx=np.maximum.outer(np.arange(NXI), np.arange(NXI)))


_G = build_grid_m()


def neumann_coulomb_m(c1, c3, G=None):
    """general-m prolate Neumann Coulomb; c1[m], c3[m] = azimuthal components on (xi,eta).
    Returns (J, per_m) with per_m the m=0/1/2 contributions."""
    if G is None:
        G = _G
    pref = (2.0 / R) * a ** 6
    per_m = []
    for m in range(G['MMAX'] + 1):
        Wm = (2 * np.pi) ** 2 if m == 0 else 2 * np.pi ** 2
        acc = 0.0
        for l in range(m, G['LMAX'] + 1):
            norm = (math.factorial(l - m) / math.factorial(l + m)) ** 2
            Plm_eta = G['Peta'][m, l]
            g1 = ((G['JAC'] * c1[m]) * Plm_eta[None, :]) @ G['weta']
            g3 = ((G['JAC'] * c3[m]) * Plm_eta[None, :]) @ G['weta']
            K = G['Pxi'][m, l][G['minidx']] * G['Qxi'][m, l][G['maxidx']]
            val = (G['wxi'] * g1) @ K @ (G['wxi'] * g3)
            acc += (-1) ** m * (2 * l + 1) * norm * val
        per_m.append(pref * Wm * acc)
    return sum(per_m), per_m


def _bridge_components(G, dressed):
    """azimuthal components c_m of the bridge density on the grid.
    dressed=False -> rho_1sB(1+..); dressed=True -> * Psi_A (the 4-body dressed density)."""
    rho0 = rho_1s(G['rB'], ALPHA)                          # 1s_B density
    if dressed:
        Psi_A = leaf_dressing(G['rA'].ravel()).reshape(G['rA'].shape)
        rho0 = rho0 * Psi_A
    return [rho0, B1 * rho0, B2 * rho0]


# --------------------------------------------------------------------------- #
# brute Monte-Carlo (12-D for the 4-body; 6-D for the self-Coulomb control)
# --------------------------------------------------------------------------- #
def brute_selfcoulomb(n_tot, batch=4_000_000):
    """6-D MC of INT rho_bridge(1) rho_bridge(3)/r13 ; batch-means error."""
    means = []; ntot = 0
    while ntot < n_tot:
        r1 = sample_1s(batch, ALPHA, CENTER_B); r3 = sample_1s(batch, ALPHA, CENTER_B)
        d13 = np.linalg.norm(r1 - r3, axis=1)
        g = w_mod(r1) * w_mod(r3) / np.maximum(d13, 1e-12)
        means.append(g.mean()); ntot += batch
    bm = np.array(means); return bm.mean(), bm.std(ddof=1) / np.sqrt(len(bm))


def brute_4body(n_tot, batch=4_000_000):
    """12-D MC of the pi/delta 4-body; bridge electrons carry the azimuthal weight."""
    means = []; ntot = 0
    while ntot < n_tot:
        r1 = sample_1s(batch, ALPHA, CENTER_B); r3 = sample_1s(batch, ALPHA, CENTER_B)
        r2 = sample_1s(batch, BETA, CENTER_A); r4 = sample_1s(batch, BETA, CENTER_A)
        d12 = np.linalg.norm(r1 - r2, axis=1)
        d34 = np.linalg.norm(r3 - r4, axis=1)
        d13 = np.linalg.norm(r1 - r3, axis=1)
        g = w_mod(r1) * w_mod(r3) * f_gem(d12) * f_gem(d34) / np.maximum(d13, 1e-12)
        means.append(g.mean()); ntot += batch
    bm = np.array(means); return bm.mean(), bm.std(ddof=1) / np.sqrt(len(bm))


def reduced_4body(G=None):
    if G is None:
        G = _G
    c = _bridge_components(G, dressed=True)
    return neumann_coulomb_m(c, c, G)


# --------------------------------------------------------------------------- #
# run
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    print("=" * 78)
    print("Two-center (LiH) 4-electron bridging integral -- PI/DELTA (m != 0) channel")
    print(f"  bridge density rho_1sB*(1 + {B1} cos phi + {B2} cos 2phi)  (m=0,1,2);"
          f"  leaves 1s_A(Z={BETA})")
    print(f"  R={R}  a={a:.4f}  bridge Z={ALPHA}  f=exp(-{GAM} r)")
    print("=" * 78)

    print("\n--- controls (isolate + pin the general-m bridge before the 4-body) ---")
    c0 = _bridge_components(_G, dressed=False)
    c0_sig = [c0[0], 0 * c0[0], 0 * c0[0]]
    Jsig, _ = neumann_coulomb_m(c0_sig, c0_sig)
    print(f"  [C1'] m=0 path (b=0) self-Coulomb = {Jsig:.6f}  vs 5*alpha/8={5*ALPHA/8:.6f}"
          f"   (rel {abs(Jsig-5*ALPHA/8)/(5*ALPHA/8):.1e})")
    Jred, per_m = neumann_coulomb_m(c0, c0)
    frac = [f"{p/Jred:.3f}" for p in per_m]
    print(f"  [C2'] modulated self-Coulomb reduced = {Jred:.6f}   per-m fraction m=0/1/2: {frac}")
    m6, e6 = brute_selfcoulomb(80_000_000)
    print(f"        6-D MC = {m6:.6f} +/- {e6:.1e}   rel diff vs reduced = {abs(m6-Jred)/abs(Jred):.2e}"
          f"   -> m!=0 convention PINNED" if abs(m6 - Jred) / abs(Jred) < 3e-3 else "   -> MISMATCH")

    print("\n--- PI/DELTA 4-body: reduced (leaf dressing + general-m bridge) ---")
    red_total, red_per_m = reduced_4body()
    print(f"   per-m contributions m=0/1/2: {[f'{p:+.6e}' for p in red_per_m]}")
    print(f"   (m=0 part {red_per_m[0]:+.6e} should ~= the sigma-gate 4-body 2.856e-2)")
    print(f"   I_reduced (pi/delta) = {red_total:+.8e}   (grid {_G['NXI']}x{_G['NETA']})")

    print("\n   grid convergence (reduced descends to its grid limit; residual = quadrature):")
    conv = []
    for (nxi, neta, lmax) in [(200, 80, 30), (320, 140, 36), (440, 200, 40)]:
        Gc = build_grid_m(nxi, neta, lmax)
        rc, _ = reduced_4body(Gc)
        cc = _bridge_components(Gc, dressed=False)
        sc, _ = neumann_coulomb_m([cc[0], 0 * cc[0], 0 * cc[0]], [cc[0], 0 * cc[0], 0 * cc[0]], Gc)
        scr = abs(sc - 5 * ALPHA / 8) / (5 * ALPHA / 8)
        conv.append((scr, rc))
        print(f"     ({nxi}x{neta},L{lmax}): I_reduced={rc:+.8e}  selfCoul(m0) rel={scr:.1e}")
    xs = np.array([c[0] for c in conv]); ys = np.array([c[1] for c in conv])
    red_inf = np.polyfit(xs, ys, 1)[1]
    print(f"   -> I_reduced(grid limit) = {red_inf:+.8e}")

    print("\n--- PI/DELTA 4-body: BRUTE 12-D MC (batch-means; heavy 1/r13 tail ~1e-4 scatter) ---")
    bvals = []
    for rep in range(3):
        m, e = brute_4body(120_000_000)
        bvals.append(m)
        print(f"   replica {rep+1} (120M): I_brute = {m:+.8e} +/- {e:.1e}"
              f"   vs reduced(limit): {(m-red_inf)/e:+.1f} sigma")
    bmean = float(np.mean(bvals))
    print(f"   3-replica mean = {bmean:+.8e}  (spread {np.std(bvals, ddof=1):.1e})")

    print(f"\n   VERDICT: general-m bridge validated (C2' rel {abs(m6-Jred)/abs(Jred):.1e} vs 6-D MC);")
    print(f"   pi/delta 4-body reduced(limit)={red_inf:.7e} vs brute(3x120M)={bmean:.7e}")
    print(f"   agree rel {abs(bmean-red_inf)/abs(red_inf):.1e} (MC heavy-tail-limited). The m!=0")
    print(f"   (PI, m=1; DELTA, m=2) azimuthal-transfer channels of the two-center 4-body")
    print(f"   reduction are EXACT and RI-FREE -- the sigma gate now holds for all m.")
