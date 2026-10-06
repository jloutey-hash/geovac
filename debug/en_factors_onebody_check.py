r"""Deterministic, independent G2 route for every PRODUCT-OF-ONE-BODY e-n element (phase 1):
the kinetic and V_ne parts of <chi_I|H W_a|psi> and of <W_a psi|H|W_b psi>, plus all W overlaps.

Independent of the analytic engine on both axes:
  * one-electron matrices on a fine Li-centred spherical Gauss product quadrature (not the prolate
    kernels grid), azimuth -> 2 pi delta_{m_x m_y};  calibrated against the engine's T1_no / Vne_no;
  * the 4-electron contraction by the direct 2x2-minor expansion (not the transition-tensor code).
Operator content (en_algebra Ops): kind on electron i in {rho, grad, lapA, rhov} x weight tokens.
Usage: python debug/en_factors_onebody_check.py 1e-2
"""
import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)

import geovac.lih_vmc as V  # noqa: E402
import en_factors_stage0_lm as E0  # noqa: E402
from geovac.lih_r12ci import pi_m16 as PM  # noqa: E402
from geovac.lih_r12ci.core_switch import CoreSwitch  # noqa: E402
from geovac.lih_r12ci import en_algebra as A  # noqa: E402
from geovac.lih_r12ci.transition_tensor import determinants  # noqa: E402

REFS = {1e-2: "lih_pi_m16_energy_trunc1e-2.npz", 1e-3: "lih_pi_m16_energy_exp1.5.npz"}


def quadrature(R, nr=500, nt=260, npi_extra=True):
    s, ws = np.polynomial.legendre.leggauss(nr); s = 0.5 * (s + 1); ws = 0.5 * ws
    rr = 2.0 * s / (1 - s); wr = ws * 2.0 / (1 - s) ** 2
    ct, wt = np.polynomial.legendre.leggauss(nt)
    Rg, Cg = np.meshgrid(rr, ct, indexing='ij')
    W = ((wr[:, None] * wt[None, :]) * Rg ** 2 * 2 * np.pi).reshape(-1)
    pts = np.stack([Rg * np.sqrt(1 - Cg ** 2), np.zeros_like(Rg), -R / 2 + Rg * Cg], axis=-1).reshape(-1, 3)
    return pts, W


def main():
    t = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-2
    T_no, occ, m_no, R = E0.switch_reference()
    ref = PM.build_primitives(PM.CFG, PM.PRIM_CACHE_PI, verbose=False)
    prims = ref['prims']
    sw = CoreSwitch.from_reference(prims, T_no, occ, m_no, R)
    pts, Wq = quadrature(R)
    M = len(prims)
    val = np.empty((M, len(Wq)), complex); grd = np.empty((M, len(Wq), 3), complex)
    lap = np.empty((M, len(Wq)), complex)
    for a, p in enumerate(prims):
        val[a], grd[a], lap[a] = V.orbital_vgl(p, pts, R)
    nv = T_no.T @ val; ng = np.einsum('ap,aGd->pGd', T_no, grd); nl = T_no.T @ lap
    w, _gw, lw = sw.vgl(pts)
    rA = np.linalg.norm(pts - np.array([0, 0, -R / 2]), axis=1); rB = np.linalg.norm(pts - np.array([0, 0, R / 2]), axis=1)
    vne = -3.0 / rA - 1.0 / rB
    same = (m_no[:, None] == m_no[None, :])
    tokf = {'w': w, 'lapw': lw}

    def mat(kind, toks):
        f = Wq.copy()
        for x in toks:
            f = f * tokf[x]
        if kind == 'rho':
            m = np.einsum('pg,qg,g->pq', np.conj(nv), nv, f)
        elif kind == 'grad':
            m = np.einsum('pgd,qgd,g->pq', np.conj(ng), ng, f)
        elif kind == 'lapA':
            m = np.einsum('pg,qg,g->pq', np.conj(nv), nl, f) - np.einsum('pg,qg,g->pq', np.conj(nl), nv, f)
        elif kind == 'rhov':
            m = np.einsum('pg,qg,g->pq', np.conj(nv), nv, f * vne)
        return (m.real * same)
    # calibration vs the engine one-electron matrices
    print(f"quadrature calibration: |S-I| {np.max(np.abs(mat('rho', ()) - np.eye(len(m_no)))):.1e}  "
          f"|T - T1_no| {np.max(np.abs(0.5 * mat('grad', ()) - _engine(t)['T1_no'])):.1e}  "
          f"|Vne - Vne_no| {np.max(np.abs(mat('rhov', ()) - _engine(t)['Vne_no'])):.1e}")
    z = np.load(os.path.join(HERE, "data", REFS[t]), allow_pickle=True)
    ia, ib, d, apairs = determinants(z['ct'], int(z['Mk']))
    N = len(d)
    cache = {}

    def M_(kind, toks):
        k = (kind, toks)
        if k not in cache:
            cache[k] = mat(kind, toks)
        return cache[k]

    def pt(a, a2, va, vb):
        p, q = apairs[a]; p2, q2 = apairs[a2]
        return va[p, p2] * vb[q, q2] - va[p, q2] * vb[q, p2] - va[q, p2] * vb[p, q2] + va[q, q2] * vb[p, p2]

    def ev(op, kind=None):
        out = np.zeros(N)
        for toks, c in op.items():
            for ei in ((None,) if kind is None else (1, 2, 3, 4)):
                v = [M_(kind if (ei is not None and e == ei) else 'rho', toks[e - 1]) for e in (1, 2, 3, 4)]
                for I in range(N):
                    out[I] += c * 0.25 * sum(d[J] * pt(ia[I], ia[J], v[0], v[1]) * pt(ib[I], ib[J], v[2], v[3])
                                             for J in range(N))
        return out
    an = np.load(os.path.join(HERE, "data", f"en_factors_{t:.0e}.npz"), allow_pickle=True)
    pab = json.loads(str(an['parts_ab']))
    W = [A.W1_op(), A.W2_op()]; LW = [A.lap_op(x) for x in W]
    rows = []
    for a in range(2):
        S = ev(W[a]); T = 0.5 * ev(W[a], 'grad') - 0.25 * ev(LW[a]) + 0.25 * ev(W[a], 'lapA')
        Vn = ev(W[a], 'rhov')
        rows += [(f"S_I[W{a+1}]", S, an['S_I'][:, a]), (f"T_I[W{a+1}]", T, an['T_I'][:, a]),
                 (f"Vne_I[W{a+1}]", Vn, an['Vne_I'][:, a])]
    for (a, b) in ((0, 0), (0, 1), (1, 1)):
        AB = A.op_mul(W[a], W[b])
        s = d @ ev(AB); g = d @ ev(AB, 'grad')
        l1 = d @ ev(A.op_add(A.op_mul(W[a], LW[b]), A.op_mul(W[b], LW[a]))); v = d @ ev(AB, 'rhov')
        p = pab[str((a, b))]
        rows += [(f"S_ab[{a}{b}]", np.array([s]), np.array([p['S']])),
                 (f"T_ab[{a}{b}]", np.array([0.5 * g - 0.25 * l1]), np.array([0.5 * p['grad'] - 0.25 * p['lap']])),
                 (f"Vne_ab[{a}{b}]", np.array([v]), np.array([p['Vne']]))]
    worst = 0.0
    for name, ind, eng in rows:
        dif = np.max(np.abs(ind - eng)); rel = dif / np.max(np.abs(eng))
        worst = max(worst, rel)
        print(f"  {name:14s} max|indep - analytic| = {dif:.2e}   rel {rel:.1e}")
    print(f"worst relative deviation (one-body route) = {worst:.1e}")


def _engine(t):
    if not hasattr(_engine, "c"):
        z = np.load(os.path.join(HERE, "data", "en_switch_ref.npz"))
        ref = PM.build_primitives(PM.CFG, PM.PRIM_CACHE_PI, verbose=False)
        U = np.linalg.lstsq(ref['Tmap'], z['T_no'], rcond=None)[0]   # canonical -> NO rotation
        _engine.c = dict(T1_no=U.T @ ref['T1_f'] @ U, Vne_no=U.T @ ref['Vne_f'] @ U)
    return _engine.c


if __name__ == "__main__":
    main()
