r"""Independent quadrature check of the grid one-body e-n matrices  w_xy = <phi_x| w |phi_y>  and
(lap w)_xy, V_ne-free: kernels-grid (prolate, as en_factors uses) vs a fine Li-centred spherical
Gauss product quadrature (axial symmetry: the azimuth integrates to 2 pi delta_{m_x m_y}).
Usage: python debug/en_factors_wxy_check.py
"""
import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)

import geovac.lih_vmc as V  # noqa: E402
import en_factors_stage0_lm as E0  # noqa: E402
from geovac.lih_r12ci import pi_m16 as PM  # noqa: E402
from geovac.lih_r12ci.core_switch import CoreSwitch, grid_switch  # noqa: E402
import geovac.lih_r12ci.ref_tensor as P0  # noqa: E402
import geovac.lih_r12ci.enumerators as P2  # noqa: E402


def main():
    T_no, occ, m_no, R = E0.switch_reference()
    ref = PM.build_primitives(PM.CFG, PM.PRIM_CACHE_PI, verbose=False)
    prims = ref['prims']
    sw = CoreSwitch.from_reference(prims, T_no, occ, m_no, R)
    # ---- kernels grid ----
    Gg = P0.grid_setup(exact=True)
    pv = PM._prim_grad_complex  # noqa
    val, grad, lap = P0.prim_on_grid(prims, R, Gg['pts'])
    gradc = PM._prim_grad_complex(prims, R, Gg['pts']); lapc = PM._prim_lap_complex(prims, R, Gg['pts'])
    nval = T_no.T @ val; ngrad = np.einsum('ap,aGd->pGd', T_no, gradc); nlap = T_no.T @ lapc
    c = sw.no_index
    w, lw = grid_switch(nval[c], ngrad[c], nlap[c], sw.phi0)
    GI = P2.GI_W
    Wg = np.einsum('pg,qg,g->pq', nval, nval, GI * w)
    Lg = np.einsum('pg,qg,g->pq', nval, nval, GI * lw)
    # ---- fine Li-centred quadrature (r, theta), phi = 0 slice x 2 pi ----
    nr, nt = 600, 300
    s, ws = np.polynomial.legendre.leggauss(nr); s = 0.5 * (s + 1); ws = 0.5 * ws
    rr = 2.0 * s / (1 - s); wr = ws * 2.0 / (1 - s) ** 2
    ct, wt = np.polynomial.legendre.leggauss(nt)
    Rg, Cg = np.meshgrid(rr, ct, indexing='ij'); Wq = (wr[:, None] * wt[None, :]) * Rg ** 2 * 2 * np.pi
    st = np.sqrt(1 - Cg ** 2)
    pts = np.stack([Rg * st, np.zeros_like(Rg), -R / 2 + Rg * Cg], axis=-1).reshape(-1, 3)
    Phi = np.array([V.orbital_vgl(p, pts, R)[0].real for p in prims])   # (M, npts), phi=0 => real
    no = T_no.T @ Phi
    wq, _, lwq = sw.vgl(pts)
    Wq_ = Wq.reshape(-1)
    Wf = np.einsum('pg,qg,g->pq', no, no, Wq_ * wq)
    Lf = np.einsum('pg,qg,g->pq', no, no, Wq_ * lwq)
    S = np.einsum('pg,qg,g->pq', no, no, Wq_)
    same = (m_no[:, None] == m_no[None, :])
    Wf = Wf * same; Lf = Lf * same; S = S * same
    print(f"fine-quadrature overlap: max|S - I| = {np.max(np.abs(S - np.eye(len(S)))):.2e}")
    print(f"w_xy:     max|grid - fine| = {np.max(np.abs(Wg * same - Wf)):.2e}   (max |w_xy| {np.max(np.abs(Wf)):.3f})")
    print(f"lapw_xy:  max|grid - fine| = {np.max(np.abs(Lg * same - Lf)):.2e}   (max |lapw_xy| {np.max(np.abs(Lf)):.3f})")
    np.set_printoptions(precision=6, suppress=True, linewidth=150)
    print("w_xy grid - fine (first 6x6):"); print((Wg * same - Wf)[:6, :6])


if __name__ == "__main__":
    main()


def transition_W1(Wxy, ct, Mk):
    """<chi_I| sum_i w(r_i) |psi> by the direct 2x2-minor expansion (independent of the tensor code)."""
    from itertools import combinations
    from geovac.lih_r12ci.transition_tensor import determinants
    ia, ib, d, apairs = determinants(ct, Mk)
    E = np.eye(Mk)

    def pt(a, a2, va, vb):
        p, q = apairs[a]; p2, q2 = apairs[a2]
        return va[p, p2] * vb[q, q2] - va[p, q2] * vb[q, p2] - va[q, p2] * vb[p, q2] + va[q, q2] * vb[p, p2]
    out = np.zeros(len(d))
    for I in range(len(d)):
        for J in range(len(d)):
            for e in range(4):
                v = [Wxy if k == e else E for k in range(4)]
                out[I] += 0.25 * d[J] * pt(ia[I], ia[J], v[0], v[1]) * pt(ib[I], ib[J], v[2], v[3])
    return out, d
