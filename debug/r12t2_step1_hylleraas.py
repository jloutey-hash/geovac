r"""r12 t^2 STEP 1, TASK A -- how much of the He-like core cusp do the t-dependent terms carry?

Diagnostic only (2026-10-05).  Re-uses the validated He-like 2e Hylleraas machinery of
geovac/prolate_allelectron/helike_cusp.py (reproduced first, bit-for-bit, on its own grid), then
re-implements the same symmetric-gradient energy on a configurable grid so that
  (i) the basis can be decomposed term by term (t-free vs t-carrying),
 (ii) exponential geminal shapes g(u)=e^{-gam u} (the molecular marriage uses e^{-1.5 u}) can be
      used, and
(iii) the MARRIAGE-ANALOG ansatz can be posed:  a radially-correlated reference
      Psi_rad = P(r1,r2) Phi0  (the orbital-CI analog), multiplied by geminal multipliers
      {Psi_rad, g Psi_rad, t^2 g Psi_rad, u t^2 g Psi_rad, ...}.

Phi0 = 1s(zeta)^2 with the frozen single-zeta exponent (Li+: 2.6875, He: 27/16, H-: 11/16),
exactly as debug/lih_additive_f12.py.  u = r12, s = r1+r2, t = r1-r2 (from the nucleus).

Run:  python debug/r12t2_step1_hylleraas.py   (float64, ~1-2 min)
"""
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


# --------------------------------------------------------------------------- #
# configurable grid (same construction as helike_cusp)
# --------------------------------------------------------------------------- #
class Grid:
    def __init__(self, N=72, NU=56, RMAX=7.0):
        x1, w1 = leggauss(N)
        r = 0.5 * RMAX * (x1 + 1.0); wr = 0.5 * RMAX * w1
        xu, wu = leggauss(NU)
        self.R1 = r[:, None, None]; self.R2 = r[None, :, None]
        ulo = np.abs(self.R1 - self.R2); uhi = self.R1 + self.R2
        self.U = ulo + (uhi - ulo) * 0.5 * (xu[None, None, :] + 1.0)
        self.W = (8 * np.pi ** 2 * (wr[:, None, None] * self.R1) * (wr[None, :, None] * self.R2)
                  * (0.5 * (uhi - ulo) * wu[None, None, :]) * self.U)
        self.C1 = (self.R1 ** 2 - self.R2 ** 2 + self.U ** 2) / (2 * self.R1 * self.U)
        self.C2 = (self.R2 ** 2 - self.R1 ** 2 + self.U ** 2) / (2 * self.R2 * self.U)

    def integ(self, g):
        return np.sum(self.W * g)


# A basis "function" is a tuple (b, db/dr1, db/dr2, db/du) on the grid (prefactor of Phi0).
def prim(G, kind, gam=1.5):
    R1, R2, U = G.R1, G.R2, G.U
    O = np.ones_like(U); Z = np.zeros_like(U)
    t = (R1 - R2) * O; s = (R1 + R2) * O
    e = np.exp(-gam * U)
    tab = {
        '1': (O, Z, Z, Z),
        'u': (U, Z, Z, O),
        'u2': (U ** 2, Z, Z, 2 * U),
        't2': (t ** 2, 2 * t, -2 * t, Z),
        's': (s, O, O, Z),
        's2': (s ** 2, 2 * s, 2 * s, Z),
        'st2': (s * t ** 2, t ** 2 + 2 * s * t, t ** 2 - 2 * s * t, Z),
        's3': (s ** 3, 3 * s ** 2, 3 * s ** 2, Z),
        't4': (t ** 4, 4 * t ** 3, -4 * t ** 3, Z),
        'ut2': (U * t ** 2, 2 * t * U, -2 * t * U, t ** 2),
        'us': (U * s, U, U, s),
        # exponential geminal family  g = e^{-gam u}
        'g': (e, Z, Z, -gam * e),
        'ug': (U * e, Z, Z, (1 - gam * U) * e),
        't2g': (t ** 2 * e, 2 * t * e, -2 * t * e, -gam * t ** 2 * e),
        'ut2g': (U * t ** 2 * e, 2 * t * U * e, -2 * t * U * e, t ** 2 * (1 - gam * U) * e),
        'sg': (s * e, e, e, -gam * s * e),
    }
    return tab[kind]


def mul(f, h):
    """product rule for two prefactors."""
    return (f[0] * h[0], f[1] * h[0] + f[0] * h[1], f[2] * h[0] + f[0] * h[2],
            f[3] * h[0] + f[0] * h[3])


def lin(coefs, funcs):
    out = [np.zeros_like(funcs[0][k]) for k in range(4)]
    for c, f in zip(coefs, funcs):
        for k in range(4):
            out[k] = out[k] + c * f[k]
    return tuple(out)


def matrices(G, Z, zeta, B):
    Phi2 = (zeta ** 3 / np.pi * np.exp(-zeta * (G.R1 + G.R2))) ** 2
    V = -Z / G.R1 - Z / G.R2 + 1.0 / G.U
    n = len(B)
    S = np.zeros((n, n)); H = np.zeros((n, n))
    for i in range(n):
        bi, di1, di2, diu = B[i]
        for j in range(i, n):
            bj, dj1, dj2, dju = B[j]
            a_i1 = di1 - zeta * bi; a_j1 = dj1 - zeta * bj
            a_i2 = di2 - zeta * bi; a_j2 = dj2 - zeta * bj
            e1 = a_i1 * a_j1 + diu * dju + (a_i1 * dju + a_j1 * diu) * G.C1
            e2 = a_i2 * a_j2 + diu * dju + (a_i2 * dju + a_j2 * diu) * G.C2
            S[i, j] = S[j, i] = G.integ(Phi2 * bi * bj)
            H[i, j] = H[j, i] = G.integ(Phi2 * (bi * bj * V + 0.5 * (e1 + e2)))
    return S, H


def gsolve(S, H, tol=1e-12):
    sv, U = np.linalg.eigh(S)
    keep = sv > tol * sv.max()
    X = U[:, keep] / np.sqrt(sv[keep])
    w, v = np.linalg.eigh(X.T @ H @ X)
    return w[0], X @ v[:, 0], np.linalg.cond(S)


def energy(G, Z, zeta, kinds, gam=1.5):
    B = [prim(G, '1', gam)] + [prim(G, k, gam) for k in kinds]
    S, H = matrices(G, Z, zeta, B)
    E, c, cond = gsolve(S, H)
    return E, H[0, 0] / S[0, 0], cond


def radial_reference(G, Z, zeta, rad_kinds):
    """optimal Psi_rad = P(r1,r2) Phi0 over a purely radial (orbital-representable) basis."""
    B = [prim(G, '1')] + [prim(G, k) for k in rad_kinds]
    S, H = matrices(G, Z, zeta, B)
    E, c, _ = gsolve(S, H)
    return E, lin(c, B)


def marriage(G, Z, zeta, P, mults, gam):
    """{P Phi0, m_k P Phi0} linear variation (P fixed) -- the analog of (1 + c(F-Fbar)) Psi_CI."""
    B = [P] + [mul(prim(G, k, gam), P) for k in mults]
    S, H = matrices(G, Z, zeta, B)
    return gsolve(S, H)[0]


SYS = [  # name, Z, zeta (frozen, as lih_additive_f12), exact non-rel, literature note
    ("Li+", 3.0, 2.6875, -7.27991, "exact -7.279913 (Pekeris-type)"),
    ("He", 2.0, 27 / 16, -2.90372, "exact -2.903724; HF -2.86168; s-limit -2.87903"),
    ("H-", 1.0, 11 / 16, -0.52775, "exact -0.527751"),
]

if __name__ == "__main__":
    out = {}
    t0 = time.time()
    # ---------------- 1. reproduce the 30.1 mHa with the production module --------------- #
    from geovac.prolate_allelectron.helike_cusp import energy as prod_energy
    E_rad_p = prod_energy(3.0, 2.6875, ['t2', 's'])[1]
    E_rich_p = prod_energy(3.0, 2.6875, ['u', 't2', 's', 'u2', 'ut2'])[1]
    print(f"[repro] production helike_cusp: Li+ radial {E_rad_p:.6f} rich {E_rich_p:.6f} "
          f"cusp = {1e3*(E_rad_p-E_rich_p):.3f} mHa   ({time.time()-t0:.1f}s)")
    out['repro'] = dict(E_rad=E_rad_p, E_rich=E_rich_p, cusp_mHa=1e3 * (E_rad_p - E_rich_p),
                        seconds=time.time() - t0)

    grids = {'prod(72,56,7)': Grid(72, 56, 7.0), 'fine(96,72,10)': Grid(96, 72, 10.0),
             'wide(110,72,22)': Grid(110, 72, 22.0)}
    G0 = grids['prod(72,56,7)']
    E_chk = energy(G0, 3.0, 2.6875, ['t2', 's'])[0]
    print(f"[repro] re-implementation on the prod grid: radial {E_chk:.6f} "
          f"(diff {E_chk-E_rad_p:.1e})")

    LINEAR = [  # (label, kinds) -- polynomial Hylleraas, linear in the coefficients
        ("ref (HF-like 1s^2)", []),
        ("t-free {s}", ['s']),
        ("t-free {u}", ['u']),
        ("t-free {u,u2}", ['u', 'u2']),
        ("t-free {u,s}", ['u', 's']),
        ("t-free {u,u2,s}", ['u', 'u2', 's']),
        ("+t2   {u,u2,s,t2}", ['u', 'u2', 's', 't2']),
        ("+ut2  {u,u2,s,ut2}", ['u', 'u2', 's', 'ut2']),
        ("rich  {u,u2,s,t2,ut2}", ['u', 'u2', 's', 't2', 'ut2']),
        ("radial {t2,s} (old base)", ['t2', 's']),
        ("{u,t2,s}", ['u', 't2', 's']),
        ("{u,u2,t2,s} (no ut2)", ['u', 'u2', 't2', 's']),
    ]
    RAD = ['s', 't2', 's2', 'st2', 's3', 't4']   # richer purely radial (orbital-limit proxy)

    for name, Z, zeta, exact, note in SYS:
        Gs = {'prod(72,56,7)': G0} if name != 'H-' else {}
        Gs.update({k: grids[k] for k in ('fine(96,72,10)', 'wide(110,72,22)')} if name == 'H-'
                  else {'fine(96,72,10)': grids['fine(96,72,10)']})
        res = {}
        for gname, G in Gs.items():
            E0 = energy(G, Z, zeta, [])[0]
            rows = {}
            print("=" * 92)
            print(f"{name} (Z={Z:g}, zeta={zeta:.4f})  grid {gname}  E0={E0:.6f}  exact {exact}"
                  f"  corr={1e3*(E0-exact):.2f} mHa   [{note}]")
            for lab, ks in LINEAR:
                E, _, cond = energy(G, Z, zeta, ks)
                rows[lab] = E
                print(f"   {lab:28s} E={E:.6f}  rec={1e3*(E0-E):7.2f} mHa  "
                      f"left={1e3*(E-exact):7.2f}  cond={cond:.0e}")
            inc = {
                't2 on t-free {u,u2,s}': rows["t-free {u,u2,s}"] - rows["+t2   {u,u2,s,t2}"],
                'ut2 on t-free {u,u2,s}': rows["t-free {u,u2,s}"] - rows["+ut2  {u,u2,s,ut2}"],
                't2+ut2 on t-free': rows["t-free {u,u2,s}"] - rows["rich  {u,u2,s,t2,ut2}"],
                'ut2 on {u,u2,s,t2}': rows["+t2   {u,u2,s,t2}"] - rows["rich  {u,u2,s,t2,ut2}"],
                't2 on {u,u2,s,ut2}': rows["+ut2  {u,u2,s,ut2}"] - rows["rich  {u,u2,s,t2,ut2}"],
                'OLD cusp = {t2,s} - rich': rows["radial {t2,s} (old base)"]
                - rows["rich  {u,u2,s,t2,ut2}"],
            }
            print("   increments (mHa):  " + "  ".join(f"{k}: {1e3*v:.2f}" for k, v in inc.items()))

            # ---- marriage analog: richer radial ref P, then geminal multipliers ----
            Erad, P = radial_reference(G, Z, zeta, RAD)
            print(f"   [marriage analog] radial ref {RAD}: E_rad={Erad:.6f} "
                  f"(rad corr {1e3*(E0-Erad):.2f} mHa, left {1e3*(Erad-exact):.2f})")
            mar = {}
            for gam in (0.75, 1.5, 3.0, 2 * zeta):
                Eg = marriage(G, Z, zeta, P, ['g'], gam)
                Eug = marriage(G, Z, zeta, P, ['ug'], gam)
                Eg_t2 = marriage(G, Z, zeta, P, ['g', 't2g'], gam)
                Eg_ut2 = marriage(G, Z, zeta, P, ['g', 'ut2g'], gam)
                Eg_both = marriage(G, Z, zeta, P, ['g', 't2g', 'ut2g'], gam)
                Eg_u = marriage(G, Z, zeta, P, ['g', 'ug'], gam)
                Eg_u_ut2 = marriage(G, Z, zeta, P, ['g', 'ug', 'ut2g'], gam)
                Eg_u_both = marriage(G, Z, zeta, P, ['g', 'ug', 't2g', 'ut2g'], gam)
                mar[gam] = dict(g=Eg, ug=Eug, g_t2g=Eg_t2, g_ut2g=Eg_ut2, g_t2g_ut2g=Eg_both,
                                g_ug=Eg_u, g_ug_ut2g=Eg_u_ut2, g_ug_t2g_ut2g=Eg_u_both)
                print(f"   gam={gam:5.3f}: Psi_rad x {{1,g}} {1e3*(Erad-Eg):6.2f} | {{1,ug}} "
                      f"{1e3*(Erad-Eug):6.2f} | +t2g {1e3*(Eg-Eg_t2):5.2f} | +ut2g "
                      f"{1e3*(Eg-Eg_ut2):5.2f} | +both {1e3*(Eg-Eg_both):5.2f} || {{g,ug}} "
                      f"{1e3*(Erad-Eg_u):6.2f} +ut2g {1e3*(Eg_u-Eg_u_ut2):5.2f} +both "
                      f"{1e3*(Eg_u-Eg_u_both):5.2f}  [left {1e3*(Eg_u_both-exact):.2f}]")
            res[gname] = dict(E0=E0, rows=rows, inc_mHa={k: 1e3 * v for k, v in inc.items()},
                              E_rad=Erad, marriage={str(k): v for k, v in mar.items()})
        out[name] = res
    out['seconds'] = time.time() - t0
    os.makedirs(os.path.join(ROOT, 'debug', 'data'), exist_ok=True)
    with open(os.path.join(ROOT, 'debug', 'data', 'r12t2_step1_hylleraas.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f"\ntotal {time.time()-t0:.1f}s -> debug/data/r12t2_step1_hylleraas.json")
