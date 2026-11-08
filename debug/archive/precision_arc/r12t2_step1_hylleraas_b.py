r"""r12 t^2 STEP 1, TASK A (part b) -- the faithful MARRIAGE ANALOG on the He-like core.

Part (a) (r12t2_step1_hylleraas.py) used a purely radial reference.  But an orbital CI is not
radial-only: u^2 = r1^2 + r2^2 - 2 r1.r2 is ORBITAL-representable (p-wave products), as is any
EVEN power of u times radial functions.  Only ODD powers of u (u, u t^2, ...) lie outside every
finite orbital CI (the Kato cusp).  So the honest analog of the M=16 orbital-CI reference is
    Psi_orb = P_even(r1, r2, u^2) Phi0        (polynomial, even in u),
and the marriage is (1 + c (g - gbar)) Psi_orb with g = e^{-gam u}.  The r12 t^2 lever adds
t^2 g Psi_orb and/or u t^2 g Psi_orb as extra multipliers.  NOTE (load-bearing): because Psi_orb
already contains t^2 terms, g Psi_orb already contains t^2 g with a FIXED ratio; the new terms buy
only the freedom to change that ratio.

Also measured: RELAXING the reference coefficients in the presence of the geminal
({orbital primitives} U {g Psi_orb}), which is the (N_det + 1)-dimensional molecular analog.

Run:  python debug/r12t2_step1_hylleraas_b.py
"""
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r12t2_step1_hylleraas import Grid, prim, mul, lin, matrices, gsolve   # noqa: E402


def mono(G, i, j, k, gam=None):
    """s^i t^(2j) u^(2k)   (+ derivatives) -- orbital-representable for any i,j,k."""
    R1, R2, U = G.R1, G.R2, G.U
    s = (R1 + R2) * np.ones_like(U); t = (R1 - R2) * np.ones_like(U)
    S = s ** i; T = t ** (2 * j); Q = U ** (2 * k)
    dS = i * s ** (i - 1) if i else np.zeros_like(U)
    dT = 2 * j * t ** (2 * j - 1) if j else np.zeros_like(U)
    dQ = 2 * k * U ** (2 * k - 1) if k else np.zeros_like(U)
    b = S * T * Q
    d1 = dS * T * Q + S * dT * Q
    d2 = dS * T * Q - S * dT * Q
    du = S * T * dQ
    return (b, d1, d2, du)


def orbital_basis(G, deg):
    """all s^i t^2j u^2k with i + 2j + 2k <= deg, excluding the constant."""
    out = []
    for i in range(deg + 1):
        for j in range(deg // 2 + 1):
            for k in range(deg // 2 + 1):
                if 0 < i + 2 * j + 2 * k <= deg:
                    out.append(((i, j, k), mono(G, i, j, k)))
    return out


SYS = [("Li+", 3.0, 2.6875, -7.27991), ("He", 2.0, 27 / 16, -2.90372)]

if __name__ == "__main__":
    t0 = time.time()
    out = {}
    for gname, G in (('prod(72,56,7)', Grid(72, 56, 7.0)), ('fine(96,72,10)', Grid(96, 72, 10.0))):
        for name, Z, zeta, exact in SYS:
            one = prim(G, '1')
            E0 = gsolve(*matrices(G, Z, zeta, [one]))[0]
            print("=" * 96)
            print(f"{name} grid {gname}: E0={E0:.6f} exact {exact} corr {1e3*(E0-exact):.2f} mHa")
            res = {}
            for deg in (2, 3, 4):
                ob = orbital_basis(G, deg)
                B = [one] + [f for _, f in ob]
                Eorb, c, cond = gsolve(*matrices(G, Z, zeta, B))
                P = lin(c, B)
                row = dict(n_orb=len(ob), E_orb=Eorb)
                print(f" deg<={deg} ({len(ob)} orbital-even prims): E_orb={Eorb:.6f} "
                      f"(orb corr {1e3*(E0-Eorb):.2f}, left {1e3*(Eorb-exact):.2f}, cond {cond:.0e})")
                for gam in (1.5, 3.0, 2 * zeta):
                    def mar(mults, relax=False):
                        BB = [P] + [mul(prim(G, k, gam), P) for k in mults]
                        if relax:
                            BB = BB + [f for _, f in ob] + [one]
                        return gsolve(*matrices(G, Z, zeta, BB))[0]
                    Eg = mar(['g']); Et = mar(['g', 't2g']); Eut = mar(['g', 'ut2g'])
                    Eb = mar(['g', 't2g', 'ut2g'])
                    Er = mar(['g'], True); Erb = mar(['g', 't2g', 'ut2g'], True)
                    Eu = mar(['g', 'ug']); Eub = mar(['g', 'ug', 't2g', 'ut2g'])
                    row[str(gam)] = dict(g=Eg, g_t2g=Et, g_ut2g=Eut, g_both=Eb, relax_g=Er,
                                         relax_g_both=Erb, g_ug=Eu, g_ug_both=Eub)
                    print(f"   gam={gam:5.3f}: marriage {1e3*(Eorb-Eg):6.2f} [left {1e3*(Eg-exact):6.2f}]"
                          f" | +t2g {1e3*(Eg-Et):5.2f} +ut2g {1e3*(Eg-Eut):5.2f} +both {1e3*(Eg-Eb):5.2f}"
                          f" | relaxCI {1e3*(Eg-Er):5.2f} relax+both {1e3*(Er-Erb):5.2f}"
                          f" | {{g,ug}} {1e3*(Eg-Eu):5.2f} +both {1e3*(Eu-Eub):5.2f}"
                          f" [left {1e3*(Erb-exact):.2f}]")
                res[deg] = row
            out[f"{name}|{gname}"] = dict(E0=E0, exact=exact, res=res)
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'r12t2_step1_hylleraas_b.json')
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f"total {time.time()-t0:.1f}s -> {path}")
