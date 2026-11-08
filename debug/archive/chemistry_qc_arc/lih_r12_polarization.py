"""Phase 1 (relaxable core): the CORE-POLARIZATION R_eq test for prolate LiH.

Phase 0 proved radial breathing is INERT (zc*(R) flat, R_eq moved 0.06 pp).  The
v5.15.2 mechanism predicts the live lever is SHAPE, not size: the core must lean
toward the bond (acquire a dipole) so its screening FOLLOWS the H atom.  This driver
adds exactly one new knob -- the core dipole d -- and minimizes E_tot over it at
each R, apples-to-apples with the frozen baseline (d=0).

The core dipole enters four consistent places:
  1. anisotropic screening felt by the valence:  + d * V_H_dip(r_A) cos(theta_A)
     (built by REUSING the validated vh_coupled with its pointwise potential swapped
      from the monopole V_H(r_A) to the dipole V_H_dip(r_A) cos(theta_A) -- same
      r12-coupled Phi_P machinery, no new integral code);
  2. distortion cost (resists polarization):     + d^2 / (2 alpha_c);
  3. H-nucleus drive (pulls the core toward H):   - Z_H * d * V_H_dip(R)   [on-axis];
  4. valence COUNTER-polarization: automatic, via d(e_val)/dd -- this is what encodes
     LiH's ionic (Li+ H-) screening of the driving field.

    E_tot(R,d) = e_val[H + VH_mono + d*VHDIP + proj]
                 + Z_Li Z_H / R - V_H_mono(R) + E_core
                 + d^2/(2 alpha_c) - Z_H d V_H_dip(R)

zc is FROZEN at 2.6875 (Phase 0: breathing inert).  alpha_c: run BOTH the parameter-
free model (9/zc^4 = 0.1725) and the true Li+ (0.1925) to bracket the magnitude.

Run from root:  python debug/lih_r12_polarization.py [jmax] [lmax] [lambda] [l_neumann]
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from geovac import prolate_recondition as pr                 # noqa: E402,F401
import prolate_r12_mpf as m                                  # noqa: E402
from heh_probe import build_basis_full                       # noqa: E402
from r12ci_first_energy import solve_canonical               # noqa: E402
from prolate_core_hartree import V_H_closed, ZC_LI           # noqa: E402
from prolate_core_dipole import (V_H_dip_closed, alpha_c_model,  # noqa: E402
                                 ALPHA_LIP_TRUE)
import lih_frozen_core_first as L                            # noqa: E402
import lih_r12_coupled as C                                  # noqa: E402

R_E_EXP = 3.015
JMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2
LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2
LAM = float(sys.argv[3]) if len(sys.argv) > 3 else 1000.0
LNEU = int(sys.argv[4]) if len(sys.argv) > 4 else 20
ALPHAS = [1.0, 1.4]
R_GRID = [2.70, 2.85, 3.015, 3.20, 3.45]
D_GRID = np.array([-0.02, 0.0, 0.005, 0.01, 0.02, 0.03, 0.045, 0.065, 0.09, 0.13])
ZC0 = ZC_LI  # frozen core exponent (2.6875)


def E_core_var(zc):
    return zc**2 - 5.375 * zc


E_CORE_ANCHOR = L.E_CORE - E_core_var(ZC0)  # anchor E_core(2.6875)=-7.2799 (validated)


def E_core(zc):
    return E_core_var(zc) + E_CORE_ANCHOR


def dip_pointwise(G, zc):
    """Per-unit-d dipole screening on the grid: (zc/4) V_H_dip(r_A) cos(theta_A)."""
    x, e, hR = G['x'], G['e'], G['hR']
    r_A = np.maximum(hR * (x + e), 1e-12)
    ct = (x * e + 1.0) / (x + e)             # cos(theta_A) at the Li focus
    return (zc / 4.0) * V_H_dip_closed(r_A, zc, 1.0) * ct


def VHdip_atR(R, zc):
    """On-axis (cos theta=1) dipole potential per unit d at the H nucleus (distance R)."""
    return (zc / 4.0) * V_H_dip_closed(R, zc, 1.0)


def assemble_mpf(R, alpha):
    basis = build_basis_full(JMAX, LMAX, alpha, p_set=(0, 1))
    S, H = m.assemble_hetero(basis, R, alpha, L.Z_LI, L.Z_H, l_neumann=LNEU, dps=30)
    return basis, S, H


def eval_over_d(basis, S, H, R, alpha):
    """e_val(d) over D_GRID for a fixed (R,alpha), frozen zc.  Returns (ev[d], drive)."""
    G = C.grid_arrays(R, alpha, N_xi=32, N_eta=24, xi_max=18.0, zc=ZC0)
    VH2 = C.vh_coupled(basis, G, alpha)
    Plam = LAM * C.projector_coupled(basis, G, alpha)
    Gd = dict(G); Gd['vH'] = dip_pointwise(G, ZC0)          # swap potential -> dipole
    VHDIP = C.vh_coupled(basis, Gd, alpha)                  # reuse validated builder
    base = H + VH2 + Plam
    ev = np.array([solve_canonical(S, base + d * VHDIP)[0] for d in D_GRID])
    return ev, VHdip_atR(R, ZC0)


def _min_over_d(ev, const_R, ac, drive):
    """min_d E_tot(d) by a deg-2 fit (E_tot ~ linear e_val + quadratic cost)."""
    Etot = ev + const_R + D_GRID**2 / (2.0 * ac) - L.Z_H * D_GRID * drive
    a, b, c = np.polyfit(D_GRID, Etot, 2)
    if a <= 0:
        i = int(np.argmin(Etot)); return Etot[i], D_GRID[i]
    dstar = float(np.clip(-b / (2 * a), D_GRID.min(), D_GRID.max()))
    return a * dstar**2 + b * dstar + c, dstar


def _req(R_grid, E):
    p = np.poly1d(np.polyfit(np.array(R_grid) - R_E_EXP, E, min(4, len(R_grid) - 1)))
    dp, ddp = p.deriv(1), p.deriv(2)
    roots = dp.r[np.isreal(dp.r)].real
    cand = [r + R_E_EXP for r in roots
            if ddp(r) > 0 and min(R_grid) < r + R_E_EXP < max(R_grid)]
    return float(min(cand, key=lambda rr: p(rr - R_E_EXP))) if cand else float('nan')


def main():
    t0 = time.time()
    ac_model = alpha_c_model(ZC0)
    print(f"CORE-POLARIZATION prolate LiH R_eq  (j,l)=({JMAX},{LMAX})  lambda={LAM:.0f} "
          f" l_neu={LNEU}   exp R_e={R_E_EXP}", flush=True)
    print(f"  frozen baseline: +3.4% (l=2)/+5.5% (l=3).  alpha_c: model={ac_model:.4f}, "
          f"true Li+={ALPHA_LIP_TRUE}.  GO if polarization pulls R_eq inward.", flush=True)
    rows = {'frozen': [], 'model': [], 'true': []}
    dstars = {'model': [], 'true': []}
    for R in R_GRID:
        const_R = L.Z_LI * L.Z_H / R - V_H_closed(R, ZC0) + E_core(ZC0)
        best = {'frozen': (None,), 'model': (None,), 'true': (None,)}
        for a in ALPHAS:
            basis, S, H = assemble_mpf(R, a)
            ev, drive = eval_over_d(basis, S, H, R, a)
            i0 = int(np.where(D_GRID == 0.0)[0][0])
            e_frozen = ev[i0] + const_R
            e_m, d_m = _min_over_d(ev, const_R, ac_model, drive)
            e_t, d_t = _min_over_d(ev, const_R, ALPHA_LIP_TRUE, drive)
            if best['frozen'][0] is None or e_frozen < best['frozen'][0]:
                best['frozen'] = (e_frozen,)
            if best['model'][0] is None or e_m < best['model'][0]:
                best['model'] = (e_m, d_m)
            if best['true'][0] is None or e_t < best['true'][0]:
                best['true'] = (e_t, d_t)
        rows['frozen'].append(best['frozen'][0])
        rows['model'].append(best['model'][0]); dstars['model'].append(best['model'][1])
        rows['true'].append(best['true'][0]); dstars['true'].append(best['true'][1])
        print(f"  R={R:.3f}  E_frozen={best['frozen'][0]:.5f}  "
              f"E_pol(model)={best['model'][0]:.5f} (d*={best['model'][1]:.4f})  "
              f"E_pol(true)={best['true'][0]:.5f} (d*={best['true'][1]:.4f})  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    Rf = _req(R_GRID, np.array(rows['frozen']))
    Rm = _req(R_GRID, np.array(rows['model']))
    Rt = _req(R_GRID, np.array(rows['true']))
    ef = (Rf - R_E_EXP) / R_E_EXP * 100
    em = (Rm - R_E_EXP) / R_E_EXP * 100
    et = (Rt - R_E_EXP) / R_E_EXP * 100
    print("\n" + "=" * 70, flush=True)
    print(f"  d*(R) model: " + " ".join(f"{r:.2f}:{d:.4f}" for r, d in zip(R_GRID, dstars['model'])),
          flush=True)
    print(f"  R_eq FROZEN        = {Rf:.4f} bohr   drift = {ef:+.2f}%", flush=True)
    print(f"  R_eq POL (model a) = {Rm:.4f} bohr   drift = {em:+.2f}%   "
          f"(moved {em-ef:+.2f} pp)", flush=True)
    print(f"  R_eq POL (true  a) = {Rt:.4f} bohr   drift = {et:+.2f}%   "
          f"(moved {et-ef:+.2f} pp)", flush=True)
    if np.isfinite(em) and np.isfinite(ef):
        v = ("GO: polarization pulls inward toward truth"
             if em < ef - 0.3 else
             "WEAK/NO: polarization inert or under-closes at this basis")
        print(f"  VERDICT: {v}", flush=True)
    print("=" * 70, flush=True)


if __name__ == '__main__':
    main()
