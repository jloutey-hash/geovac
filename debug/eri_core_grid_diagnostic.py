"""Route C, increment C2 SIZING diagnostic (measurement, NOT a build).

QUESTION: increment-7 showed the tight Li 1s core's ONE-BODY integrals swing
0.36 Ha across R on the two-center prolate GRID (r_A=(R/2)(xi+eta) resolves the
tight core R-dependently). C1 solved that ANALYTICALLY (1e-50, R-flat). C2 must
now decide the TWO-ELECTRON (ERI) core-involving integrals: are they ALSO
R-inaccurate on the grid (=> C2 must build a full analytic mixed-exponent Neumann
V_ee, hard), or does the ERI integration smooth the sharp core enough that the
grid is R-accurate for ERIs (=> reuse the existing grid ERIs for core classes,
easy hybrid)?

METHOD: put the tight Li core chi_c (STO, zc=2.6875) and a diffuse valence probe
chi_v (STO, zv=0.8), BOTH on the Li focus (r_A=(R/2)(xi+eta)), on the SAME grid
the real all-electron FCI uses (build_mo_integrals_multiexp: N_grid=44,
xi_max=15.0), and compute four ERI classes via the existing grid machinery
compute_vee_integral across R in {2.70,2.85,3.015,3.20,3.45}. Report R-spread.

Orbitals are ANALYTICALLY normalized (psi = sqrt(zeta^3/pi) e^{-zeta r_A}, the
exact 3D-normalized m=0 STO in this code's convention where psi_2d IS Psi_3D and
INT psi_2d^2 J dxi deta * 2pi = 1). This isolates the ERI *integration* error
(the C2 question) from the normalization error (already cured by C1). The grid
norm INT psi^2 J w * 2pi is reported separately as context.

Exact references (all four classes are spherical one-center distributions):
  (cc|cc) = 5 zc / 8   (R-independent);  (vv|vv) = 5 zv / 8.
  cross classes computed by high-precision mpmath 1D radial quadrature.

Run from root:  python debug/eri_core_grid_diagnostic.py
"""
import os
import sys
import numpy as np
from mpmath import mp, mpf, quad, e as mp_e, pi as mp_pi, sqrt as mp_sqrt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from geovac.prolate_scf import compute_vee_integral   # noqa: E402

ZC = 3.0 - 5.0 / 16.0      # tight Li^2+-like 1s core exponent (Slater), 2.6875
ZV = 0.8                   # diffuse valence probe exponent

R_LIST = (2.70, 2.85, 3.015, 3.20, 3.45)


# ---------------------------------------------------------------------------
# grid replicating get_orbital_on_grid EXACTLY, with an analytic STO on it
# ---------------------------------------------------------------------------
def sto_orbital(zeta, R, N_grid, xi_max):
    """orbital dict {psi,xi,eta,w_xi,w_eta,R} for chi = sqrt(z^3/pi) e^{-z r_A},
    r_A=(R/2)(xi+eta), evaluated on the code's standard quadrature grid.
    psi is ANALYTICALLY normalized (the exact 3D-normalized m=0 STO)."""
    eta, w_eta = np.polynomial.legendre.leggauss(N_grid)
    u_gl, w_u_gl = np.polynomial.legendre.leggauss(N_grid)
    t = (u_gl + 1) / 2
    xi = 1.0 + (xi_max - 1.0) * t ** 2
    w_xi = w_u_gl * (xi_max - 1.0) * t
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    r_A = (R / 2.0) * (XI + ETA)
    psi = np.sqrt(zeta ** 3 / np.pi) * np.exp(-zeta * r_A)
    return {'psi': psi, 'xi': xi, 'eta': eta, 'w_xi': w_xi, 'w_eta': w_eta, 'R': R}


def grid_norm(orb):
    """INT psi^2 J dxi deta * 2pi on the grid (should be 1 if grid resolves it)."""
    R = orb['R']
    XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
    J = (R / 2) ** 3 * (XI ** 2 - ETA ** 2)
    WX, WE = np.meshgrid(orb['w_xi'], orb['w_eta'], indexing='ij')
    return float(np.sum(orb['psi'] ** 2 * J * WX * WE) * 2 * np.pi)


# ---------------------------------------------------------------------------
# exact references via high-precision radial quadrature (mpmath)
# ---------------------------------------------------------------------------
def _V_exp(r, gamma, Q):
    """potential at r of a spherical density rho = A e^{-2 gamma r'} with total
    charge Q:  V(r) = Q (1/r)[1 - (1 + gamma r) e^{-2 gamma r}]."""
    return Q * (1.0 / r) * (1 - (1 + gamma * r) * mp_e ** (-2 * gamma * r))


def exact_cross(za, zb):
    """cross Coulomb (aa|bb): J = INT rho_b(r) V_a(r) d3r,
    rho_a=(za^3/pi)e^{-2za r} (Q=1), rho_b likewise. mpmath."""
    def integrand(r):
        rho_b = (mpf(zb) ** 3 / mp_pi) * mp_e ** (-2 * mpf(zb) * r)
        return 4 * mp_pi * r ** 2 * rho_b * _V_exp(r, mpf(za), mpf(1))
    return quad(integrand, [0, mp.inf])


def exact_overlap_self(za, zb):
    """(cv|cv): self-Coulomb of overlap density rho_cv = chi_a chi_b
    = sqrt(za^3 zb^3)/pi e^{-(za+zb) r}, a spherical A e^{-2 gamma r} with
    gamma=(za+zb)/2, charge S = 8 sqrt(za^3 zb^3)/(za+zb)^3.
    U = S^2 * 5 gamma / 8."""
    za, zb = mpf(za), mpf(zb)
    gamma = (za + zb) / 2
    S = 8 * mp_sqrt(za ** 3 * zb ** 3) / (za + zb) ** 3
    return S ** 2 * mpf(5) * gamma / 8


# ---------------------------------------------------------------------------
# the ERI classes via the existing grid machinery
#   compute_vee_integral(o1,o2,o3,o4) = INT (o1 o3)(1) 1/r12 (o2 o4)(2)
#   chemist (ab|cd)=INT a(1)b(1)/r12 c(2)d(2)  ->  call (a, c, b, d)
# ---------------------------------------------------------------------------
def eri_classes(c, v):
    cccc = compute_vee_integral(c, c, c, c)              # (cc|cc)
    ccvv = compute_vee_integral(c, v, c, v)              # (cc|vv): o1o3=cc, o2o4=vv
    cvcv = compute_vee_integral(c, c, v, v)              # (cv|cv): o1o3=cv, o2o4=cv
    vvvv = compute_vee_integral(v, v, v, v)              # (vv|vv)
    return cccc, ccvv, cvcv, vvvv


def run(N_grid, xi_max, log):
    def out(s=""):
        print(s, flush=True)
        log.write(s + "\n")

    mp.dps = 40
    ex_cccc = float(mpf(5) * mpf(ZC) / 8)
    ex_vvvv = float(mpf(5) * mpf(ZV) / 8)
    ex_ccvv = float(exact_cross(ZC, ZV))
    ex_cvcv = float(exact_overlap_self(ZC, ZV))

    out(f"=== ERI core-grid R-accuracy diagnostic  (N_grid={N_grid}, xi_max={xi_max}) ===")
    out(f"  zc={ZC:.4f} (tight Li core)   zv={ZV:.4f} (diffuse valence)")
    out(f"  EXACT refs (R-independent):  (cc|cc)={ex_cccc:.6f}  (cc|vv)={ex_ccvv:.6f}"
        f"  (cv|cv)={ex_cvcv:.6f}  (vv|vv)={ex_vvvv:.6f}  Ha")
    out("")
    hdr = (f"  {'R':>6} | {'(cc|cc)':>11} {'d_exact':>9} | {'(cc|vv)':>11} | "
           f"{'(cv|cv)':>11} | {'(vv|vv)':>11} | {'gridN_c':>9} {'gridN_v':>8}")
    out(hdr)
    out("  " + "-" * (len(hdr) - 2))

    vals = {k: [] for k in ('cccc', 'ccvv', 'cvcv', 'vvvv', 'Nc', 'Nv')}
    for R in R_LIST:
        c = sto_orbital(ZC, R, N_grid, xi_max)
        v = sto_orbital(ZV, R, N_grid, xi_max)
        cccc, ccvv, cvcv, vvvv = eri_classes(c, v)
        Nc, Nv = grid_norm(c), grid_norm(v)
        for k, x in zip(('cccc', 'ccvv', 'cvcv', 'vvvv', 'Nc', 'Nv'),
                        (cccc, ccvv, cvcv, vvvv, Nc, Nv)):
            vals[k].append(x)
        d_cccc = (cccc - ex_cccc) * 1e3   # mHa dev from exact
        out(f"  {R:6.3f} | {cccc:11.6f} {d_cccc:8.2f}m | {ccvv:11.6f} | "
            f"{cvcv:11.6f} | {vvvv:11.6f} | {Nc:9.5f} {Nv:8.5f}")

    out("")
    out("  R-SPREAD (max-min):")
    labels = {'cccc': '(cc|cc)', 'ccvv': '(cc|vv)', 'cvcv': '(cv|cv)', 'vvvv': '(vv|vv)'}
    exed = {'cccc': ex_cccc, 'ccvv': ex_ccvv, 'cvcv': ex_cvcv, 'vvvv': ex_vvvv}
    for k in ('cccc', 'ccvv', 'cvcv', 'vvvv'):
        arr = np.array(vals[k])
        spread_mHa = (arr.max() - arr.min()) * 1e3
        mean_dev_mHa = (arr.mean() - exed[k]) * 1e3
        flag = ">> NEEDS ANALYTIC" if spread_mHa > 1.0 else ("ok(<0.1)" if spread_mHa < 0.1 else "borderline")
        out(f"    {labels[k]:>9}:  spread = {spread_mHa:9.3f} mHa   "
            f"mean_dev_from_exact = {mean_dev_mHa:9.3f} mHa   [{flag}]")
    nc = np.array(vals['Nc']); nv = np.array(vals['Nv'])
    out(f"    grid-norm chi_c: {nc.min():.5f}..{nc.max():.5f}  spread={float(nc.max()-nc.min()):.2e}"
        f"   (analytic norm = 1; deviation = grid mis-resolution of the tight core)")
    out(f"    grid-norm chi_v: {nv.min():.5f}..{nv.max():.5f}  spread={float(nv.max()-nv.min()):.2e}")
    out("")
    return vals


if __name__ == '__main__':
    os.makedirs('debug/data', exist_ok=True)
    with open('debug/data/eri_core_grid_diagnostic.log', 'w') as log:
        # PRIMARY: the real LiH FCI resolution
        run(N_grid=44, xi_max=15.0, log=log)
        # secondary resolution (the H2 / build_mo_integrals resolution)
        run(N_grid=48, xi_max=14.0, log=log)
        # refinement check: does the (cc|cc) spread SHRINK with more points
        # (grid error, like the one-body core 0.36->0.05) or is it stable?
        run(N_grid=96, xi_max=15.0, log=log)
        run(N_grid=160, xi_max=15.0, log=log)
