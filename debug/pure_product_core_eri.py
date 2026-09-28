r"""Analysis #4 (debug/sprint_product_basis_tradeoff_scoping_memo.md): can a PURE prolate
product basis (xi-Laguerre x eta-Legendre, one shared alpha, as in Papers 11/12) carry the
Li 1s core on the TWO-BODY side, so that Route C's atom-centred STO core is unnecessary?

Basis (m = 0, full rectangular tensor, Nx x Ne):
    X_a(xi) = L_a(2 alpha (xi-1)) e^{-alpha xi},  Y_b(eta) = P_b(eta).
Core = lowest Galerkin eigenvector of the BARE one-centre h(Z_A, Z_B = 0) (H nucleus ABSENT),
whose exact eigenfunction is the hydrogenic 1s e^{-Z r_A}:
    Z = 3      -> (cc|cc) exact = 5Z/8 = 15/8  (step 2)
    Z = 2.6875 -> e^{-2.6875 r_A} = Route C's core STO exactly (step 3; 5Z/8 = 1.6796875)

Two-body side = Paper 12's Neumann route, general engine primitives REUSED from
geovac.neumann_vee_general_m (closed-form A moments, closed-form B seeds {E1, gamma, ln},
forward l-recurrence, the _W_poly IBP split) exactly as prolate_recondition._build_Xtab_mp /
debug/prolate_mixed_eri.build_Xtab_pair assemble them.  The only new thing is the ORDER of
the sums: the IBP tail correction is contracted as  sum_n Bsum[P+n] G_P[n]  with
G_P[n] = sum_{j>=n} W_j j!/(n! c^{j-n+1})  (an exact regrouping; gated bit-level against
_build_Xtab_mp / build_Xtab_pair below), and the X table is contracted straight against the
pair-density vectors instead of being stored.

    (ab|cd) = (2/R)(R/2)^6 (2 pi)^2 sum_l (2l+1)  w_l^{ab} . X_l(c_ab, c_cd) . w_l^{cd},
    w_l[P] = sum_Q rho[P-2,Q] Y_l[Q] - rho[P,Q] Y_l[Q+2],   Y_l[Q] = int eta^Q P_l,
    rho = monomial coeffs of phi_a phi_b / e^{-c xi}.
The l-sum TERMINATES at l = max Q (eta-polynomial densities): no l_neumann truncation.

Transcendentals (named, per the tagging rule): the basis-projected ERI is a finite
Q-combination of e^{-c}, e^{c} E_1(2c), Euler gamma, ln c, ln 2 at the density rates c
(the Neumann Q_0 log-moment seeds of Paper 12 / neumann_vee_general_m._L_moments); these
cancel to the rational 5Z/8 only in the complete-basis limit.

Usage: python debug/pure_product_core_eri.py [gate|cc|cross|all] ; log -> debug/data/pure_product_core_eri.log
"""
from __future__ import annotations

import json
import os
import sys
import time
from typing import Dict, List, Sequence, Tuple

import mpmath as mp
import numpy as np
from scipy.linalg import eigh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from geovac import neumann_vee_general_m as ngm  # noqa: E402
from geovac import prolate_recondition as pr  # noqa: E402

LOG: List[str] = []


def log(s: str) -> None:
    print(s, flush=True)
    LOG.append(s)


# ----------------------------------------------------------------------- eta moments
def _Me(q: int) -> mp.mpf:
    return mp.mpf(2) / (q + 1) if q % 2 == 0 else mp.mpf(0)


_YCACHE: Dict[Tuple[int, int, int], mp.mpf] = {}


def Ymom(l: int, Q: int) -> mp.mpf:
    """int_{-1}^1 eta^Q P_l(eta) d eta (exact rational)."""
    key = (l, Q, mp.mp.dps)
    if key not in _YCACHE:
        if l > Q or (Q + l) % 2:
            v = mp.mpf(0)
        else:
            v = pr._mom_eta(ngm._shift(list(ngm._leg_coeffs(l)), Q))
        _YCACHE[key] = v
    return _YCACHE[key]


# ----------------------------------------------------------------------- one-body (orthogonal basis)
def coeff_mats(Nx: int, Ne: int, alpha) -> Tuple[List[List[mp.mpf]], List[List[mp.mpf]]]:
    Cx = [pr.laguerre_coeffs(a, alpha, Nx) for a in range(Nx)]
    Ce = [pr.legendre_coeffs(b, Ne) for b in range(Ne)]
    return Cx, Ce


def core_orbital(Nx: int, Ne: int, alpha, R, Z) -> Tuple[List[List[mp.mpf]], float, float]:
    """Galerkin 1s of h(Z_A=Z, Z_B=0) in the Nx x Ne basis.
    Returns (K monomial coeff matrix [j][l] of the UNIT-normalized orbital / e^{-alpha xi},
             E_1s float, cond of the downcast overlap)."""
    alpha, R, Z = mp.mpf(alpha), mp.mpf(R), mp.mpf(Z)
    Cx, Ce = coeff_mats(Nx, Ne, alpha)
    A = ngm._mono_moments(2 * alpha, 2 * Nx + 4)
    # derivative coeffs of X_a e^{alpha xi}:  p' - alpha p
    Dx = []
    for a in range(Nx):
        p = Cx[a]
        d = [(k + 1) * p[k + 1] for k in range(Nx - 1)] + [mp.mpf(0)]
        Dx.append([d[k] - alpha * p[k] for k in range(Nx)])
    De = []
    for b in range(Ne):
        p = Ce[b]
        De.append([(k + 1) * p[k + 1] for k in range(Ne - 1)] + [mp.mpf(0)])

    def xmat(P, Qm, shift, minus1=False):
        M = np.zeros((Nx, Nx))
        for a in range(Nx):
            for b in range(a, Nx):
                s = mp.mpf(0)
                for i in range(Nx):
                    if P[a][i] == 0:
                        continue
                    for j in range(Nx):
                        if Qm[b][j] == 0:
                            continue
                        t = A[i + j + shift]
                        if minus1:
                            t -= A[i + j]
                        s += P[a][i] * Qm[b][j] * t
                M[a, b] = M[b, a] = float(s)
        return M

    def emat(P, Qm, shift, oneminus=False):
        M = np.zeros((Ne, Ne))
        for a in range(Ne):
            for b in range(a, Ne):
                s = mp.mpf(0)
                for i in range(Ne):
                    if P[a][i] == 0:
                        continue
                    for j in range(Ne):
                        if Qm[b][j] == 0:
                            continue
                        t = _Me(i + j + shift)
                        if oneminus:
                            t = _Me(i + j) - _Me(i + j + 2)
                        s += P[a][i] * Qm[b][j] * t
                M[a, b] = M[b, a] = float(s)
        return M

    Sx0, Sx1, Sx2 = (xmat(Cx, Cx, k) for k in (0, 1, 2))
    Tx = xmat(Dx, Dx, 2, minus1=True)                 # int (xi^2-1) X'X'
    Ey0, Ey1, Ey2 = (emat(Ce, Ce, k) for k in (0, 1, 2))
    Ty = emat(De, De, 0, oneminus=True)               # int (1-eta^2) Y'Y'
    h3 = float((R / 2) ** 3)
    tp = float(2 * mp.pi)
    Rf, Zf = float(R), float(Z)
    S = h3 * tp * (np.kron(Sx2, Ey0) - np.kron(Sx0, Ey2))
    T = 0.5 * (4 / Rf ** 2) * h3 * tp * (np.kron(Tx, Ey0) + np.kron(Sx0, Ty))
    V = -Zf * (Rf / 2) ** 2 * tp * (np.kron(Sx1, Ey0) - np.kron(Sx0, Ey1))
    H = T + V
    dS = np.sqrt(np.diag(S))
    Sn = S / np.outer(dS, dS)
    Hn = H / np.outer(dS, dS)
    ev = np.linalg.eigvalsh(Sn)
    w, U = eigh(Hn, Sn)
    v = U[:, 0] / dS
    v = v.reshape(Nx, Ne)
    vm = [[mp.mpf(float(v[a, b])) for b in range(Ne)] for a in range(Nx)]
    # K = Cx^T v Ce  (monomial coeffs), mpf
    tmp = [[sum((vm[a][b] * Ce[b][l] for b in range(Ne)), mp.mpf(0)) for l in range(Ne)]
           for a in range(Nx)]
    K = [[sum((Cx[a][j] * tmp[a][l] for a in range(Nx)), mp.mpf(0)) for l in range(Ne)]
         for j in range(Nx)]
    rho = conv2(K, K)
    nrm = overlap_from_rho(rho, 2 * alpha, R)
    sc = 1 / mp.sqrt(nrm)
    K = [[sc * x for x in row] for row in K]
    # sign: positive at the Li nucleus side (xi=1, eta=-1)
    val = sum((K[j][l] * (-1) ** l for j in range(Nx) for l in range(Ne)), mp.mpf(0))
    if val < 0:
        K = [[-x for x in row] for row in K]
    return K, float(w[0]), float(ev.max() / ev.min())


def conv2(K1: Sequence[Sequence], K2: Sequence[Sequence]) -> List[List[mp.mpf]]:
    """rho[P][Q] = sum K1[j][l] K2[P-j][Q-l]."""
    n1, m1, n2, m2 = len(K1), len(K1[0]), len(K2), len(K2[0])
    rho = [[mp.mpf(0)] * (m1 + m2 - 1) for _ in range(n1 + n2 - 1)]
    for j in range(n1):
        for l in range(m1):
            a = K1[j][l]
            if a == 0:
                continue
            for jj in range(n2):
                rj = rho[j + jj]
                for ll in range(m2):
                    b = K2[jj][ll]
                    if b != 0:
                        rj[l + ll] += a * b
    return rho


def overlap_from_rho(rho, c, R) -> mp.mpf:
    A = ngm._mono_moments(c, len(rho) + 3)
    s = mp.mpf(0)
    for P, row in enumerate(rho):
        for Q, x in enumerate(row):
            if x != 0:
                s += x * (A[P + 2] * _Me(Q) - A[P] * _Me(Q + 2))
    return 2 * mp.pi * (mp.mpf(R) / 2) ** 3 * s


def orbital_K(orb) -> List[List[mp.mpf]]:
    """Route C (prolate_mixed_eri.Orbital / OrbitalM, mu=0) -> monomial matrix K[j][l]."""
    K = [[mp.mpf(0)] * len(orb.eta_poly) for _ in range(orb.xi_power + 1)]
    for l, e in enumerate(orb.eta_poly):
        K[orb.xi_power][l] = orb.norm * mp.mpf(e)
    return K


# ----------------------------------------------------------------------- the Neumann contraction
def _Gvec(l: int, P: int, c) -> List[mp.mpf]:
    """G_P[n] = sum_{j>=n} W_j j!/(n! c^{j-n+1}),  W = xi^P P_l(xi) (= ngm._W_poly(l,0,0,P))."""
    W = ngm._W_poly(l, 0, 0, P)
    deg = len(W) - 1
    G = [mp.mpf(0)] * (deg + 1)
    ic = 1 / c
    for n in range(deg + 1):
        s = mp.mpf(0)
        f = ic                      # j!/(n! c^{j-n+1}) at j = n
        for j in range(n, deg + 1):
            if j > n:
                f = f * j * ic
            if W[j] != 0:
                s += W[j] * f
        G[n] = s
    return G


def eri_rho(rho1, c1, rho2, c2, R, verbose: bool = False) -> mp.mpf:
    """(1|2) for two m=0 pair densities rho_i(xi,eta) e^{-c_i xi} (monomial coeff matrices)."""
    c1, c2, R = mp.mpf(c1), mp.mpf(c2), mp.mpf(R)
    q_hi = max(len(rho1[0]), len(rho2[0])) + 1          # max eta power incl. eta^2 shift
    l_hi = q_hi
    p_max = max(len(rho1), len(rho2)) + 1                # max xi power incl. xi^2 shift
    csum = c1 + c2
    deg_extra = p_max + l_hi
    p_corr_max = p_max + deg_extra
    n_mono = p_max + l_hi + 2
    t0 = time.time()
    Am1 = ngm._mono_moments(c1, n_mono)
    Am2 = ngm._mono_moments(c2, n_mono)
    B1 = ngm._B_table(0, 0, l_hi, p_max, c1)
    B2 = B1 if c2 == c1 else ngm._B_table(0, 0, l_hi, p_max, c2)
    Bs = ngm._B_table(0, 0, l_hi, p_corr_max, csum)
    tB = time.time() - t0

    def wvec(rho, l):
        w = [mp.mpf(0)] * (p_max + 1)
        for P in range(p_max + 1):
            s = mp.mpf(0)
            if 0 <= P - 2 < len(rho):
                for Q, x in enumerate(rho[P - 2]):
                    if x != 0:
                        y = Ymom(l, Q)
                        if y != 0:
                            s += x * y
            if P < len(rho):
                for Q, x in enumerate(rho[P]):
                    if x != 0:
                        y = Ymom(l, Q + 2)
                        if y != 0:
                            s -= x * y
            w[P] = s
        return w

    total = mp.mpf(0)
    for l in range(l_hi + 1):
        w1 = wvec(rho1, l)
        w2 = wvec(rho2, l)
        if all(x == 0 for x in w1) or all(x == 0 for x in w2):
            continue
        A1 = [ngm._A_moment(l, 0, 0, P, Am1) for P in range(p_max + 1)]
        A2 = [ngm._A_moment(l, 0, 0, P, Am2) for P in range(p_max + 1)]
        wA1 = sum(w1[P] * A1[P] for P in range(p_max + 1))
        wA2 = sum(w2[P] * A2[P] for P in range(p_max + 1))
        wB1 = sum(w1[P] * B1[(l, P)] for P in range(p_max + 1))
        wB2 = sum(w2[P] * B2[(l, P)] for P in range(p_max + 1))
        # J1: electron 1 inner (P_l, rate c1), electron 2 outer (Q_l, rate c2)
        # corr1 = sum_P1 w1 sum_n G1_{P1}[n] sum_P2 w2 Bs[P2+n]
        nmax = p_max + l + 1
        bw2 = [sum(w2[P] * Bs[(l, P + n)] for P in range(p_max + 1)) for n in range(nmax)]
        bw1 = [sum(w1[P] * Bs[(l, P + n)] for P in range(p_max + 1)) for n in range(nmax)]
        corr1 = mp.mpf(0)
        corr2 = mp.mpf(0)
        for P in range(p_max + 1):
            if w1[P] != 0:
                G = _Gvec(l, P, c1)
                corr1 += w1[P] * sum(G[n] * bw2[n] for n in range(len(G)))
            if w2[P] != 0:
                G = _Gvec(l, P, c2)
                corr2 += w2[P] * sum(G[n] * bw1[n] for n in range(len(G)))
        Xl = wA1 * wB2 - corr1 + wA2 * wB1 - corr2
        total += (2 * l + 1) * Xl
    if verbose:
        log(f"      [eri_rho: p_max={p_max} l_hi={l_hi} B tables {tB:.1f}s total {time.time() - t0:.1f}s]")
    return (2 / R) * (R / 2) ** 6 * (2 * mp.pi) ** 2 * total


# ----------------------------------------------------------------------- gates
def gate() -> bool:
    """G1: my contraction == direct contraction of the corpus X tables (_build_Xtab_mp,
    build_Xtab_pair) on random small densities.  G2: Route C (cc|vv),(cv|cv),(cc|cc)
    reproduced from Route C's own orbital objects through eri_rho."""
    import prolate_mixed_eri as pmx
    ok = True
    mp.mp.dps = 50
    rng = np.random.default_rng(1)
    log("GATE G1: eri_rho vs corpus X tables (random rho, p<=5, q<=5)")
    for (a1, a2) in [(0.8, 0.8), (0.8, 1.3)]:
        r1 = [[mp.mpf(float(x)) for x in rng.standard_normal(5)] for _ in range(4)]
        r2 = [[mp.mpf(float(x)) for x in rng.standard_normal(5)] for _ in range(4)]
        R = mp.mpf(2.0)
        mine = eri_rho(r1, 2 * a1, r2, 2 * a2, R)
        l_hi, p_max = 6, 5
        caps = {(0, 0): l_hi}
        if a1 == a2:
            Xt = pr._build_Xtab_mp([(0, 0)], l_hi, p_max, a1, caps)
        else:
            Xt = pmx.build_Xtab_pair([(0, 0)], l_hi, p_max, a1, a2, caps)
        tot = mp.mpf(0)
        for l in range(l_hi + 1):
            X = Xt[(l, 0, 0)]
            for sg1, dP1, dQ1 in [(1, 2, 0), (-1, 0, 2)]:
                for sg2, dP2, dQ2 in [(1, 2, 0), (-1, 0, 2)]:
                    for p1 in range(4):
                        for q1 in range(5):
                            y1 = Ymom(l, q1 + dQ1)
                            if y1 == 0:
                                continue
                            for p2 in range(4):
                                for q2 in range(5):
                                    y2 = Ymom(l, q2 + dQ2)
                                    if y2 == 0:
                                        continue
                                    tot += (sg1 * sg2 * (2 * l + 1) * r1[p1][q1] * r2[p2][q2]
                                            * X[p1 + dP1][p2 + dP2] * y1 * y2)
        ref = (2 / R) * (R / 2) ** 6 * (2 * mp.pi) ** 2 * tot
        rel = float(abs(mine - ref) / abs(ref))
        log(f"  alpha=({a1},{a2}): mine={mp.nstr(mine, 20)} corpus={mp.nstr(ref, 20)} rel={rel:.1e}")
        ok &= rel < 1e-35
    log("GATE G2: Route C orbitals through eri_rho vs cached Route C ERIs (R=3.015)")
    d = np.load(os.path.join(HERE, "data", "prolate_separation_constant_stageC_ints.npz"))
    eri = d["eri"]
    R = mp.mpf("3.015")
    core = pmx.sto_orbital(pmx.ZC_LI, R, is_core=True)
    bond = pmx.valence_prolate_orbital(0, 0, mp.mpf("1.0"))
    Kc, Kb = orbital_K(core), orbital_K(bond)
    for name, (Ka, aa, Kb_, ab), (Kc_, ac, Kd, ad), ref in [
            ("(cc|vv)", (Kc, core.alpha, Kc, core.alpha), (Kb, bond.alpha, Kb, bond.alpha), eri[0, 0, 3, 3]),
            ("(cv|cv)", (Kc, core.alpha, Kb, bond.alpha), (Kc, core.alpha, Kb, bond.alpha), eri[0, 3, 0, 3]),
            ("(cc|cc)", (Kc, core.alpha, Kc, core.alpha), (Kc, core.alpha, Kc, core.alpha), eri[0, 0, 0, 0])]:
        val = eri_rho(conv2(Ka, Kb_), aa + ab, conv2(Kc_, Kd), ac + ad, R)
        log(f"  {name}: eri_rho={mp.nstr(val, 16)}  RouteC cache={ref:.16f}  diff={float(val) - ref:+.2e}")
        if name != "(cc|cc)":
            ok &= abs(float(val) - ref) < 1e-12
    log(f"GATE {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------- runs
def run_cc(sizes=(10, 15, 20), Rs=(2.5, 3.015, 3.5), Z=3.0, dps=None, out=None):
    res = []
    for N in sizes:
        for amode in ("R/2", "1.2"):
            for R in Rs:
                al = 0.5 * R if amode == "R/2" else 1.2
                use_dps = dps or (40 + 8 * N)
                mp.mp.dps = use_dps
                _YCACHE.clear()
                t0 = time.time()
                K, E, cS = core_orbital(N, N, mp.mpf(al), mp.mpf(R), mp.mpf(Z))
                t1 = time.time()
                rho = conv2(K, K)
                J = eri_rho(rho, 2 * mp.mpf(al), rho, 2 * mp.mpf(al), mp.mpf(R))
                t2 = time.time()
                exact = 5 * mp.mpf(Z) / 8
                err = float(J - exact)
                row = dict(N=N, alpha=amode, R=R, dps=use_dps, E1s_err=E + Z * Z / 2,
                           cc=float(J), err_mHa=err * 1e3, t_1body=t1 - t0, t_eri=t2 - t1,
                           condS=cS)
                res.append(row)
                log(f"  N={N:2d}x{N:<2d} alpha={amode:4s} R={R:5.3f} dps={use_dps}: "
                    f"E1s+Z^2/2={E + Z * Z / 2:+.2e}  (cc|cc)={mp.nstr(J, 14)}  "
                    f"err={err * 1e3:+.3e} mHa  [1b {t1 - t0:.1f}s, eri {t2 - t1:.1f}s]")
                if out:
                    with open(out, "w") as f:
                        json.dump(res, f, indent=1)
    return res


def run_cross(sizes=(10, 15, 20), dps=None, out=None):
    """Z = 2.6875 core (= Route C's STO exactly, in the complete-basis limit) x Route C's
    exact bond orbital (xi^0 eta^0 e^{-xi}) and H 1s; R = 3.015; vs cached Route C ERIs."""
    import prolate_mixed_eri as pmx
    import prolate_allelectron_analytic_fci as paf
    d = np.load(os.path.join(HERE, "data", "prolate_separation_constant_stageC_ints.npz"))
    eri = d["eri"]
    R = 3.015
    Z = float(pmx.ZC_LI)
    res = []
    for N in sizes:
        for amode in ("R/2", "1.2"):
            al = 0.5 * R if amode == "R/2" else 1.2
            use_dps = dps or (40 + 8 * N)
            mp.mp.dps = use_dps
            _YCACHE.clear()
            Rm = mp.mpf(R)
            t0 = time.time()
            K, E, cS = core_orbital(N, N, mp.mpf(al), Rm, mp.mpf(Z))
            rcc = conv2(K, K)
            ac = mp.mpf(al)
            Jcc = eri_rho(rcc, 2 * ac, rcc, 2 * ac, Rm)
            row = dict(N=N, alpha=amode, dps=use_dps, E1s_err=E + Z * Z / 2,
                       cc=float(Jcc), cc_err_mHa=float(Jcc - 5 * mp.mpf(Z) / 8) * 1e3)
            for vname, orb, idx in (("bond", pmx.valence_prolate_orbital(0, 0, mp.mpf("1.0")), 3),
                                    ("H1s", paf.sto_orbital_B(mp.mpf("1.0"), Rm), 2)):
                Kv = orbital_K(orb)
                av = orb.alpha
                rvv = conv2(Kv, Kv)
                rcv = conv2(K, Kv)
                ccvv = eri_rho(rcc, 2 * ac, rvv, 2 * av, Rm)
                cvcv = eri_rho(rcv, ac + av, rcv, ac + av, Rm)
                d1 = float(ccvv) - eri[0, 0, idx, idx]
                d2 = float(cvcv) - eri[0, idx, 0, idx]
                row[f"ccvv_{vname}"] = float(ccvv)
                row[f"cvcv_{vname}"] = float(cvcv)
                row[f"dccvv_{vname}_mHa"] = d1 * 1e3
                row[f"dcvcv_{vname}_mHa"] = d2 * 1e3
                log(f"  N={N:2d} alpha={amode:4s} v={vname:4s}: (cc|vv)={float(ccvv):.10f} "
                    f"(RouteC {eri[0, 0, idx, idx]:.10f}, d={d1 * 1e3:+.3e} mHa)  "
                    f"(cv|cv)={float(cvcv):.10f} (RouteC {eri[0, idx, 0, idx]:.10f}, d={d2 * 1e3:+.3e} mHa)")
            row["t_total"] = time.time() - t0
            log(f"  N={N:2d} alpha={amode:4s} Z=2.6875: E1s err {E + Z * Z / 2:+.2e}  "
                f"(cc|cc)-5Z/8 = {row['cc_err_mHa']:+.3e} mHa   [{row['t_total']:.1f}s]")
            res.append(row)
            if out:
                with open(out, "w") as f:
                    json.dump(res, f, indent=1)
    return res


def xtab_single_rate(c, p_max: int, l_hi: int):
    """Stored X_l (single density rate c, both electrons) via the G-regrouped IBP split --
    the object a production single-alpha engine would build ONCE for all ERIs."""
    c = mp.mpf(c)
    csum = 2 * c
    Am = ngm._mono_moments(c, p_max + l_hi + 2)
    B1 = ngm._B_table(0, 0, l_hi, p_max, c)
    Bs = ngm._B_table(0, 0, l_hi, 2 * p_max + l_hi, csum)
    X = {}
    for l in range(l_hi + 1):
        A1 = [ngm._A_moment(l, 0, 0, P, Am) for P in range(p_max + 1)]
        Gs = [_Gvec(l, P, c) for P in range(p_max + 1)]
        C = [[sum(G[n] * Bs[(l, P2 + n)] for n in range(len(G))) for P2 in range(p_max + 1)]
             for G in Gs]                                  # C[P1][P2] = corr(W[P1], P2)
        X[l] = [[A1[P1] * B1[(l, P2)] - C[P1][P2] + A1[P2] * B1[(l, P1)] - C[P2][P1]
                 for P2 in range(p_max + 1)] for P1 in range(p_max + 1)]
    return X


def contract_stored(X, rho1, rho2, R, p_max: int, l_hi: int) -> mp.mpf:
    tot = mp.mpf(0)
    for l in range(l_hi + 1):
        def wv(rho):
            w = [mp.mpf(0)] * (p_max + 1)
            for P in range(p_max + 1):
                s = mp.mpf(0)
                if 0 <= P - 2 < len(rho):
                    s += sum(x * Ymom(l, Q) for Q, x in enumerate(rho[P - 2]) if x != 0)
                if P < len(rho):
                    s -= sum(x * Ymom(l, Q + 2) for Q, x in enumerate(rho[P]) if x != 0)
                w[P] = s
            return w
        w1, w2 = wv(rho1), wv(rho2)
        Xl = X[l]
        tot += (2 * l + 1) * sum(w1[a] * sum(Xl[a][b] * w2[b] for b in range(p_max + 1))
                                 for a in range(p_max + 1) if w1[a] != 0)
    R = mp.mpf(R)
    return (2 / R) * (R / 2) ** 6 * (2 * mp.pi) ** 2 * tot


def run_timing(N=20, dps=120):
    """Measured cost split: one-time X_l build vs per-ERI contraction (single alpha)."""
    mp.mp.dps = dps
    R = mp.mpf("3.015")
    al = R / 2
    K, E, cS = core_orbital(N, N, al, R, mp.mpf(3))
    rho = conv2(K, K)
    p_max, l_hi = 2 * N, 2 * N
    t0 = time.time()
    X = xtab_single_rate(2 * al, p_max, l_hi)
    t1 = time.time()
    J = contract_stored(X, rho, rho, R, p_max, l_hi)
    t2 = time.time()
    J2 = eri_rho(rho, 2 * al, rho, 2 * al, R)
    log(f"TIMING N={N} dps={dps}: X_l build (once) {t1 - t0:.1f}s; per-ERI contraction "
        f"(incl. w vectors) {t2 - t1:.2f}s; stored-X (cc|cc)={mp.nstr(J, 15)} vs eri_rho "
        f"{mp.nstr(J2, 15)} (rel {float(abs(J - J2) / J2):.1e}); overlap cond {cS:.2e}")


def run_dps_check(N=20):
    """Precision convergence: same (cc|cc) at two dps."""
    for dps in (100, 140):
        run_cc(sizes=(N,), Rs=(3.015,), dps=dps)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    tag = sys.argv[2] if len(sys.argv) > 2 else mode
    D = os.path.join(HERE, "data")
    try:
        if mode in ("gate", "all"):
            gate()
        if mode in ("cc", "all"):
            log("STEP 2: (cc|cc), Z=3 bare Li 1s (H absent), exact 5Z/8 = 1.875")
            sizes = tuple(int(x) for x in sys.argv[3].split(",")) if len(sys.argv) > 3 else (10, 15, 20)
            run_cc(sizes=sizes, out=os.path.join(D, f"pure_product_core_{tag}.json"))
        if mode in ("cross", "all"):
            log("STEP 3: Z=2.6875 core x Route C valence, R=3.015")
            sizes = tuple(int(x) for x in sys.argv[3].split(",")) if len(sys.argv) > 3 else (10, 15, 20)
            run_cross(sizes=sizes, out=os.path.join(D, f"pure_product_core_{tag}.json"))
        if mode == "timing":
            run_timing(int(sys.argv[3]) if len(sys.argv) > 3 else 20)
        if mode == "dps":
            run_dps_check(int(sys.argv[3]) if len(sys.argv) > 3 else 20)
    finally:
        with open(os.path.join(D, f"pure_product_core_eri_{tag}.log"), "w") as f:
            f.write("\n".join(LOG))
