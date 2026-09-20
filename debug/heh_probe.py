"""B-probe: does exact-algebraic r12 survive a heteronuclear (Z_A != Z_B) center?
HeH+ (He: Z=2 at focus (xi+eta), H: Z=1 at focus (xi-eta)), 2 electrons.
(1) validate the NEW heteronuclear V_ne block vs quadrature (the eta term);
(2) first HeH+ energy from the full engine (all r12 blocks reused unchanged).
Run from root:  python debug/heh_probe.py
"""
import os
import sys
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geovac import prolate_recondition as pr
import prolate_r12_mpf as m
from r12ci_first_energy import solve_canonical

ZA, ZB = 2, 1        # He at focus A=(xi+eta), H at focus B=(xi-eta)
R = 1.4632           # HeH+ equilibrium bond length (bohr)
E_REF = -2.97869     # HeH+ X^1Sigma+ BO total energy near R_e (LOAD-BEARING: verify)


def build_basis_full(j_max, l_max, alpha, p_set=(0, 1)):
    """FULL angular parity (l+m even AND odd) -- heteronuclear needs both."""
    rlist = [(j, k) for j in range(j_max + 1) for k in range(j_max + 1)]
    alist = [(l, mm) for l in range(l_max + 1) for mm in range(l_max + 1)]
    return [(pr.ProductFn(j, l, k, mm, 0, alpha), p)
            for p in p_set for (j, k) in rlist for (l, mm) in alist]


def vne_hetero_quad(basis_p, R, alpha, ZA, ZB):
    from geovac.hylleraas import build_quadrature_grids
    g = build_quadrature_grids(N_xi=26, N_eta=18, N_phi=20, xi_max=15.0)
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    dphi, wphi = g['dphi'], g['w_phi']
    hR = R / 2.0
    n = len(basis_p)
    V = np.zeros((n, n))
    for a in range(len(xi)):
        x1 = xi[a]
        for cc in range(len(xi)):
            x2 = xi[cc]
            ef = np.exp(-2.0 * alpha * (x1 + x2))
            for b in range(len(eta)):
                e1 = eta[b]; Jp1 = x1**2 - e1**2
                r1A = hR * (x1 + e1); r1B = hR * abs(x1 - e1)
                for d in range(len(eta)):
                    e2 = eta[d]; Jp2 = x2**2 - e2**2
                    r2A = hR * (x2 + e2); r2B = hR * abs(x2 - e2)
                    vne = -(ZA / r1A + ZB / r1B + ZA / r2A + ZB / r2B)
                    rho1 = np.sqrt(max((x1**2 - 1) * (1 - e1**2), 0.0))
                    rho2 = np.sqrt(max((x2**2 - 1) * (1 - e2**2), 0.0))
                    Aa = (x1*e1 - x2*e2)**2 + rho1**2 + rho2**2
                    Bc = 2 * rho1 * rho2
                    r12 = hR * np.sqrt(np.maximum(Aa - Bc * np.cos(dphi), 0.0))
                    wgt = wxi[a]*wxi[cc]*weta[b]*weta[d]*Jp1*Jp2*(hR**6)*ef*vne
                    for i in range(n):
                        bi, pi = basis_p[i]
                        gi = x1**bi.j * x2**bi.k * e1**bi.l * e2**bi.m
                        for jjj in range(n):
                            bj, pj = basis_p[jjj]
                            gj = x1**bj.j * x2**bj.k * e1**bj.l * e2**bj.m
                            Pp = pi + pj
                            rP = r12**Pp if Pp != 0 else np.ones_like(dphi)
                            V[i, jjj] += wgt * gi * gj * np.sum(rP * wphi) * 2 * np.pi
    return V


def validate_vne():
    a = 1.5
    mspec = [(0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1)]  # incl. odd parity
    mixed = ([(pr.ProductFn(j, l, k, mm, 0, a), 0) for (j, l, k, mm) in mspec]
             + [(pr.ProductFn(j, l, k, mm, 0, a), 1) for (j, l, k, mm) in mspec])
    Vm = m.vne_hetero_mpf(mixed, R, a, ZA, ZB)
    Vq = vne_hetero_quad(mixed, R, a, ZA, ZB)
    rel = np.max(np.abs(Vm - Vq) / np.maximum(np.abs(Vq), 1e-9))
    print("mpf hetero V_ne diag:", np.round(np.diag(Vm), 5))
    print("quad hetero V_ne diag:", np.round(np.diag(Vq), 5))
    print(f"max rel diff = {rel:.3e}   (quad grid-limited; expect ~1e-3)")
    print("HETERO V_ne", "VALIDATED" if rel < 5e-3 else "MISMATCH")
    return rel < 5e-3


def energy(j_max, l_max, alpha, p_set=(0, 1), use_mpf=True):
    basis = build_basis_full(j_max, l_max, alpha, p_set)
    t = time.time()
    if use_mpf:
        Smp, Hmp = m.assemble_hetero(basis, R, alpha, ZA, ZB, mpf_out=True)
        e = m.solve_canonical_mpf(Smp, Hmp)[0]
    else:
        S, H = m.assemble_hetero(basis, R, alpha, ZA, ZB)
        e = solve_canonical(S, H)[0]
    et = e + ZA * ZB / R      # nuclear repulsion Z_A Z_B / R
    print(f"  (j={j_max},l={l_max}) a={alpha:.2f} n={len(basis)} "
          f"E_tot={et:.6f}  err={(E_REF-et)*1000:+.3f} mHa  [{time.time()-t:.0f}s]",
          flush=True)
    return et


if __name__ == "__main__":
    ok = validate_vne()
    if ok:
        print("\nHeH+ energies (E_ref =", E_REF, "Ha, nuclear = ZA ZB/R =",
              round(ZA * ZB / R, 4), "):")
        for (jm, lm) in [(2, 1), (2, 2)]:
            for a in [1.3, 1.6, 2.0]:
                energy(jm, lm, a)
