"""H2 Neumann-r12 build, increment 8 (START of the mpf physical-accuracy engine).

Extends Paper 12's mpf engine (geovac/prolate_recondition.py) with r12. Mapping:
  EVEN r12^(2s): (R/2)^(2s) <(A-BcosDphi)^s>_phi -> polynomial A^s, a SUM of
      separable products -> reuse _ov/_vne one-electron moments with power shifts.
  ODD  r12^(2s+1): (A K0 - B K1)-type -> reuse the vee_mp Neumann X-table (m=0,1).

This module starts with the EVEN OVERLAP block (r12^2, p1xp1) reusing _ov, and
validates it against the crude float64 overlap. Later increments add the odd
overlap, V_ne, V_ee, kinetic, then the re-basing + conditioned solve.

Run from debug/.
"""
import numpy as np
import mpmath as mp
from geovac import prolate_recondition as pr
from geovac import neumann_vee_general_m as ngm
from geovac.hylleraas import HylleraasBasisFunction, build_quadrature_grids, compute_overlap_matrix

mp.mp.dps = 40
R = 1.4011
ALPHA = 1.0

# A (the r12^2 geometry) as separable one-electron power shifts:
# A = (xi1 eta1 - xi2 eta2)^2 + (xi1^2-1)(1-eta1^2) + (xi2^2-1)(1-eta2^2)
# term = (coef, dxi1, deta1, dxi2, deta2)
A_TERMS = [
    (1, 2, 2, 0, 0), (-2, 1, 1, 1, 1), (1, 0, 0, 2, 2),
    (1, 2, 0, 0, 0), (-1, 2, 2, 0, 0), (-1, 0, 0, 0, 0), (1, 0, 2, 0, 0),
    (1, 0, 0, 2, 0), (-1, 0, 0, 2, 2), (-1, 0, 0, 0, 0), (1, 0, 0, 0, 2),
]


def overlap_even_mpf(basis_p, R, alpha, dps=40):
    """mpf overlap for a basis of r12-power-p functions, EVEN P=p_i+p_j only.
    Each entry: P=0 -> h6 cc o1 o2 ; P=2 -> h6 (R/2)^2 cc sum_A coef o1' o2'."""
    with mp.workdps(dps):
        c = mp.mpf(2.0 * alpha)
        # A_k(c) moments, generous upper index
        n_max = 4 * max(max(b.j, b.k) for b, _ in basis_p) + 8
        A = ngm._mono_moments(c, n_max)
        h6 = (mp.mpf(R) / 2) ** 6
        hR2 = (mp.mpf(R) / 2) ** 2
        n = len(basis_p)
        S = np.empty((n, n), object)
        for i in range(n):
            bi, pi = basis_p[i]
            for jj in range(i, n):
                bj, pj = basis_p[jj]
                mu = 0  # sigma
                cc = pr._phi_cc(mu)
                P = pi + pj
                j1, k1, l1, m1 = bi.j + bj.j, bi.k + bj.k, bi.l + bj.l, bi.m + bj.m
                if P == 0:
                    val = h6 * cc * pr._ov(j1, l1, mu, A) * pr._ov(k1, m1, mu, A)
                elif P == 2:
                    acc = mp.mpf(0)
                    for (co, d1, e1, d2, e2) in A_TERMS:
                        acc += co * pr._ov(j1 + d1, l1 + e1, mu, A) * pr._ov(k1 + d2, m1 + e2, mu, A)
                    val = h6 * hR2 * cc * acc
                else:
                    raise ValueError("odd P not in this increment")
                S[i, jj] = S[jj, i] = val
        return np.array([[float(S[i, j]) for j in range(n)] for i in range(n)])


def _cl_m0(l, Q):
    """mpf eta moment int eta^Q P_l(eta) deta (m=0)."""
    return pr._mom_eta(pr._pm(pr._shift([mp.mpf(1)], Q), list(ngm._leg_coeffs(l))))


def _cl_m1(l, Q):
    """mpf eta moment int eta^Q (1-eta^2) P_l'(eta) deta (m=1, s=1)."""
    poly = pr._pm(pr._meta2(1), list(ngm._RP_poly(l, 1)))
    return pr._mom_eta(pr._pm(pr._shift([mp.mpf(1)], Q), poly))


def vee_r12_odd_mpf(basis_p, R, alpha, l_neumann=30, dps=40):
    """mpf V_ee for a PURE p=1 basis: <g_i g_j r12^1> = A K0 - B K1, exact.
    Reuses pr._build_Xtab_mp for X_l^{0,0} (=build_X(0)) and X_l^{1,1} (=build_X(1))."""
    with mp.workdps(dps):
        p_max = 2 * max(max(b.j, b.k) for b, _ in basis_p) + 4
        q_max = 2 * max(max(b.l, b.m) for b, _ in basis_p) + 4
        l_neu = min(l_neumann, q_max)
        ms_pairs = [(0, 0), (1, 1)]
        l_caps = {(0, 0): l_neu, (1, 1): l_neu}
        Xtab = pr._build_Xtab_mp(ms_pairs, l_neu, p_max, alpha, l_caps)
        cl1 = [mp.mpf(0)] + [(2 * l + 1) * (mp.factorial(l - 1) / mp.factorial(l + 1)) ** 2
                             for l in range(1, l_neu + 1)]
        h6 = (mp.mpf(R) / 2) ** 6
        pref = h6 * (mp.mpf(R) / 2) ** 2 * (2 * mp.pi) ** 2 * (2 / mp.mpf(R))
        n = len(basis_p)
        V = np.empty((n, n), object)
        for i in range(n):
            bi, _ = basis_p[i]
            for jj in range(i, n):
                bj, _ = basis_p[jj]
                p1x, q1 = bi.j + bj.j, bi.l + bj.l
                p2x, q2 = bi.k + bj.k, bi.m + bj.m
                jac = [(1, p1x + 2, q1, p2x + 2, q2), (-1, p1x + 2, q1, p2x, q2 + 2),
                       (-1, p1x, q1 + 2, p2x + 2, q2), (1, p1x, q1 + 2, p2x, q2 + 2)]
                tot = mp.mpf(0)
                # A K0 term  (mono0 = jac x A_TERMS ; X_l^{0,0}, C_l m=0)
                for (jc, P1, Q1, P2, Q2) in jac:
                    for (ac, d1, e1, d2, e2) in A_TERMS:
                        co = jc * ac
                        PP1, QQ1, PP2, QQ2 = P1 + d1, Q1 + e1, P2 + d2, Q2 + e2
                        for l in range(0, l_neu + 1):
                            X = Xtab.get((l, 0, 0))
                            if X is None:
                                continue
                            c1 = _cl_m0(l, QQ1)
                            if c1 == 0:
                                continue
                            c2 = _cl_m0(l, QQ2)
                            if c2 == 0:
                                continue
                            tot += (2 * l + 1) * co * X[PP1][PP2] * c1 * c2
                # B K1 term  (mono1 = jac ; X_l^{1,1}, D_l m=1)
                for (jc, P1, Q1, P2, Q2) in jac:
                    for l in range(1, l_neu + 1):
                        X = Xtab.get((l, 1, 1))
                        if X is None:
                            continue
                        d1v = _cl_m1(l, Q1)
                        if d1v == 0:
                            continue
                        d2v = _cl_m1(l, Q2)
                        if d2v == 0:
                            continue
                        tot += 2 * cl1[l] * jc * X[P1][P2] * d1v * d2v
                V[i, jj] = V[jj, i] = pref * tot
        return np.array([[float(V[i, j]) for j in range(n)] for i in range(n)])


def _M_eta(q):
    return mp.mpf(2) / (q + 1) if q % 2 == 0 else mp.mpf(0)


# P_Vne = xi1(xi2^2-eta2^2) + xi2(xi1^2-eta1^2), applied to base powers (p1x,q1,p2x,q2)
def _pvne_terms(p1x, q1, p2x, q2):
    return [(1, p1x + 1, q1, p2x + 2, q2), (-1, p1x + 1, q1, p2x, q2 + 2),
            (1, p1x + 2, q1, p2x + 1, q2), (-1, p1x, q1 + 2, p2x + 1, q2)]


def vne_mpf(basis_p, R, alpha, l_neumann=30, dps=40, to_float=True, godd=None):
    """mpf V_ne for a mixed p={0,1} basis (homonuclear H2).  Even blocks
    (p0p0 r12^0, p1p1 r12^2) factorize via P_Vne; odd block (p0p1 r12^1) uses the
    shared A K0 - B K1 (G_odd) assembly with P_Vne in place of the Jacobian."""
    with mp.workdps(dps):
        A = ngm._mono_moments(mp.mpf(2.0 * alpha), 4 * max(max(b.j, b.k) for b, _ in basis_p) + 8)
        if godd is None:
            godd = _odd_context(basis_p, alpha, l_neumann)
        hR = mp.mpf(R) / 2
        twopi2 = (2 * mp.pi) ** 2
        pref_e0 = -(mp.mpf(R) ** 2 / 2) * hR ** 3 * twopi2          # V_ne r12^0
        pref_e2 = -(mp.mpf(R) ** 2 / 2) * hR ** 5 * twopi2          # V_ne r12^2 (x (R/2)^2 A)
        pref_odd = -mp.mpf(R) * hR ** 5 * twopi2                    # V_ne r12^1
        n = len(basis_p)
        V = np.empty((n, n), object)
        for i in range(n):
            bi, pi = basis_p[i]
            for jj in range(i, n):
                bj, pj = basis_p[jj]
                p1x, q1, p2x, q2 = bi.j + bj.j, bi.l + bj.l, bi.k + bj.k, bi.m + bj.m
                base = _pvne_terms(p1x, q1, p2x, q2)
                P = pi + pj
                if P == 0:      # even, r12^0
                    val = mp.mpf(0)
                    for (co, P1, Q1, P2, Q2) in base:
                        val += co * A[P1] * _M_eta(Q1) * A[P2] * _M_eta(Q2)
                    val *= pref_e0
                elif P == 2:    # even, r12^2 -> x A
                    val = mp.mpf(0)
                    for (co, P1, Q1, P2, Q2) in base:
                        for (ac, d1, e1, d2, e2) in A_TERMS:
                            val += co * ac * A[P1 + d1] * _M_eta(Q1 + e1) * A[P2 + d2] * _M_eta(Q2 + e2)
                    val *= pref_e2
                else:           # odd, r12^1 -> A K0 - B K1 (shared G_odd)
                    val = pref_odd * _odd_value(base, godd)
                V[i, jj] = V[jj, i] = val
        return V if not to_float else _tofloat(V)


# ==========================================================================
# HETERONUCLEAR V_ne (B probe: does exact-algebraic r12 survive Z_A != Z_B?).
# Per electron  Z_A/r_iA + Z_B/r_iB = (2/R)[(Z_A+Z_B)xi + (Z_B-Z_A)eta]/(xi^2-eta^2),
# still POLYNOMIAL after the Jacobian.  The (Z_B-Z_A) eta term is new -- it breaks
# gerade<->ungerade, so a heteronuclear basis needs BOTH angular parities and V_ne
# couples them.  Everything else in the engine is charge-independent, so this is
# the ONLY new block.  Focus A = (R/2)(xi+eta) carries Z_A; focus B = (R/2)(xi-eta)
# carries Z_B.  Reduces to vne_mpf at Z_A=Z_B=1.
# ==========================================================================
def _pvne_hetero_terms(p1x, q1, p2x, q2, cA, cB):
    """cA = (Z_A+Z_B)/2 (xi part), cB = (Z_B-Z_A)/2 (eta part), applied to base."""
    xi_terms = [(cA, p1x + 1, q1, p2x + 2, q2), (-cA, p1x + 1, q1, p2x, q2 + 2),
                (cA, p1x + 2, q1, p2x + 1, q2), (-cA, p1x, q1 + 2, p2x + 1, q2)]
    eta_terms = [(cB, p1x, q1 + 1, p2x + 2, q2), (-cB, p1x, q1 + 1, p2x, q2 + 2),
                 (cB, p1x + 2, q1, p2x, q2 + 1), (-cB, p1x, q1 + 2, p2x, q2 + 1)]
    return xi_terms + eta_terms


def vne_hetero_mpf(basis_p, R, alpha, Z_A, Z_B, l_neumann=30, dps=40,
                   to_float=True, godd=None):
    """mpf heteronuclear V_ne for a mixed p={0,1} basis (2 electrons, 2 centers).
    Same block structure as vne_mpf; base = _pvne_hetero_terms (xi + eta)."""
    with mp.workdps(dps):
        A = ngm._mono_moments(mp.mpf(2.0 * alpha), 4 * max(max(b.j, b.k) for b, _ in basis_p) + 8)
        if godd is None:
            godd = _odd_context(basis_p, alpha, l_neumann)
        cA = (mp.mpf(Z_A) + Z_B) / 2
        cB = (mp.mpf(Z_B) - Z_A) / 2
        hR = mp.mpf(R) / 2
        twopi2 = (2 * mp.pi) ** 2
        pref_e0 = -(mp.mpf(R) ** 2 / 2) * hR ** 3 * twopi2
        pref_e2 = -(mp.mpf(R) ** 2 / 2) * hR ** 5 * twopi2
        pref_odd = -mp.mpf(R) * hR ** 5 * twopi2
        n = len(basis_p)
        V = np.empty((n, n), object)
        for i in range(n):
            bi, pi = basis_p[i]
            for jj in range(i, n):
                bj, pj = basis_p[jj]
                p1x, q1, p2x, q2 = bi.j + bj.j, bi.l + bj.l, bi.k + bj.k, bi.m + bj.m
                base = _pvne_hetero_terms(p1x, q1, p2x, q2, cA, cB)
                P = pi + pj
                if P == 0:
                    val = mp.mpf(0)
                    for (co, P1, Q1, P2, Q2) in base:
                        val += co * A[P1] * _M_eta(Q1) * A[P2] * _M_eta(Q2)
                    val *= pref_e0
                elif P == 2:
                    val = mp.mpf(0)
                    for (co, P1, Q1, P2, Q2) in base:
                        for (ac, d1, e1, d2, e2) in A_TERMS:
                            val += co * ac * A[P1 + d1] * _M_eta(Q1 + e1) * A[P2 + d2] * _M_eta(Q2 + e2)
                    val *= pref_e2
                else:
                    val = pref_odd * _odd_value(base, godd)
                V[i, jj] = V[jj, i] = val
        return V if not to_float else _tofloat(V)


def assemble_hetero(basis_p, R, alpha, Z_A, Z_B, l_neumann=40, dps=40, mpf_out=False):
    """Full (S, H) for a 2e heteronuclear diatomic (e.g. HeH+).  Identical to
    assemble_mixed except V_ne is heteronuclear.  T, V_ee, S are charge-independent
    and reuse the homonuclear engine unchanged -- that reuse IS the B-probe result."""
    with mp.workdps(dps):
        gs = [b for b, _ in basis_p]
        ps = [p for _, p in basis_p]
        n = len(basis_p)
        godd = _odd_context(basis_p, alpha, l_neumann)
        Kev0 = _kern_even(gs, R, alpha, 0, dps)
        Kev2 = _kern_even(gs, R, alpha, 2, dps)
        Kod1 = _kern_odd1(gs, R, alpha, l_neumann, dps, godd=godd)
        Kinv = pr.vee_mp(gs, alpha, R, l_neumann)
        Tm = kinetic_mixed_mpf(basis_p, R, alpha, l_neumann, dps, to_float=False, godd=godd)
        Vne = vne_hetero_mpf(basis_p, R, alpha, Z_A, Z_B, l_neumann, dps, to_float=False, godd=godd)
        S = np.empty((n, n), object)
        Vee = np.empty((n, n), object)
        for i in range(n):
            for jj in range(n):
                P = ps[i] + ps[jj]
                S[i, jj] = Kev0[i, jj] if P == 0 else (Kod1[i, jj] if P == 1 else Kev2[i, jj])
                Vee[i, jj] = Kinv[i, jj] if P == 0 else (Kev0[i, jj] if P == 1 else Kod1[i, jj])
        H = Tm + Vne + Vee
        if mpf_out:
            return S, H
        return _tofloat(S), _tofloat(H)


# ==========================================================================
# TWO-EXPONENT (two-block) support -- closes the HeH+ single-alpha diagnosis.
# A cross-block <g_i(a_a)|Op|g_j(a_b)> has product exponent e^{-(a_a+a_b)(xi1+xi2)}
# on BOTH electrons (each g carries e^{-a(xi1+xi2)}), so S/V_ne/V_ee reuse the
# single-exponent machinery at c = a_a + a_b.  Only the KINETIC gradient differs:
# d/dxi hits each function's OWN exponent.  p=0 only (the diagnosis is radial,
# not r12).
# ==========================================================================
def _kin2_grad(ja, la, jb, lb, aa, ab, A):
    """Two-exponent sigma (mu=0) one-electron kinetic gradient moment; A = moments
    at c = aa + ab.  d/dxi of g_i uses exponent aa, of g_j uses ab."""
    aa, ab = mp.mpf(aa), mp.mpf(ab)

    def mx(jx, a):
        t = pr._shift([-a], jx)
        if jx > 0:
            t = pr._pa(t, pr._shift([mp.mpf(jx)], jx - 1))
        return t

    def my(lx):
        return pr._shift([mp.mpf(lx)], lx - 1) if lx > 0 else [mp.mpf(0)]

    xip = pr._pm(pr._pm(mx(ja, aa), mx(jb, ab)), pr._xi2m1(1))
    etap = pr._pm(pr._pm(my(la), my(lb)), pr._meta2(1))
    px = pr._shift([mp.mpf(1)], ja + jb)
    py = pr._shift([mp.mpf(1)], la + lb)
    return (pr._mom_xi(xip, A) * pr._mom_eta(py) + pr._mom_xi(px, A) * pr._mom_eta(etap))


def assemble_hetero_2block(A_funcs, B_funcs, R, ZA, ZB, l_neumann=40, dps=40,
                           mpf_out=False):
    """Full (S, H) for a 2e heteronuclear diatomic on a TWO-EXPONENT basis
    (A_funcs at alpha_a, B_funcs at alpha_b), p=0 only.  Per-pair exponent
    c = alpha_i + alpha_j.  V_ee stitched from three vee_mp calls."""
    with mp.workdps(dps):
        basis = list(A_funcs) + list(B_funcs)
        n = len(basis)
        na = len(A_funcs)
        aa = A_funcs[0].alpha
        ab = B_funcs[0].alpha
        aeff = (aa + ab) / 2.0
        jmax = max(max(b.j, b.k) for b in basis)
        # A-moment tables by exponent-sum c
        Atab = {}
        for c in {2 * aa, aa + ab, 2 * ab}:
            Atab[c] = ngm._mono_moments(mp.mpf(c), 4 * jmax + 8)
        h6 = (mp.mpf(R) / 2) ** 6
        cc = 4 * mp.pi ** 2
        pref_T00 = mp.mpf('0.5') * (4 / mp.mpf(R) ** 2) * h6
        hR = mp.mpf(R) / 2
        pref_e0 = -(mp.mpf(R) ** 2 / 2) * hR ** 3 * (2 * mp.pi) ** 2
        cAc = (mp.mpf(ZA) + ZB) / 2
        cBc = (mp.mpf(ZB) - ZA) / 2
        # V_ee: three blocks stitched (within-a, within-b, cross at aeff)
        Vaa = pr.vee_mp(list(A_funcs), aa, R, l_neumann)
        Vbb = pr.vee_mp(list(B_funcs), ab, R, l_neumann)
        Vcross = pr.vee_mp(basis, aeff, R, l_neumann)     # cross sub-block valid
        S = np.empty((n, n), object)
        H = np.empty((n, n), object)
        for i in range(n):
            bi = basis[i]
            ai = bi.alpha
            for jj in range(n):
                bj = basis[jj]
                aj = bj.alpha
                Ac = Atab[ai + aj]
                # overlap
                o1 = pr._ov(bi.j + bj.j, bi.l + bj.l, 0, Ac)
                o2 = pr._ov(bi.k + bj.k, bi.m + bj.m, 0, Ac)
                S[i, jj] = h6 * cc * o1 * o2
                # kinetic (two-exp)
                k1 = _kin2_grad(bi.j, bi.l, bj.j, bj.l, ai, aj, Ac)
                k2 = _kin2_grad(bi.k, bi.m, bj.k, bj.m, ai, aj, Ac)
                T = pref_T00 * cc * (k1 * o2 + o1 * k2)
                # V_ne (heteronuclear, P=0)
                base = _pvne_hetero_terms(bi.j + bj.j, bi.l + bj.l,
                                          bi.k + bj.k, bi.m + bj.m, cAc, cBc)
                vne = mp.mpf(0)
                for (co, P1, Q1, P2, Q2) in base:
                    vne += co * Ac[P1] * _M_eta(Q1) * Ac[P2] * _M_eta(Q2)
                vne *= pref_e0
                # V_ee (stitched)
                if i < na and jj < na:
                    vee = Vaa[i, jj]
                elif i >= na and jj >= na:
                    vee = Vbb[i - na, jj - na]
                else:
                    vee = Vcross[i, jj]
                H[i, jj] = T + vne + vee
        if mpf_out:
            return S, H
        return _tofloat(S), _tofloat(H)


def vne_quad(basis_p, R, alpha):
    """Direct quad <g_i g_j r12^(pi+pj) V_ne> (unsym single-term g), homonuclear."""
    g = build_quadrature_grids(N_xi=26, N_eta=18, N_phi=24, xi_max=15.0)
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
                    vne = -(1 / r1A + 1 / r1B + 1 / r2A + 1 / r2B)
                    rho1 = np.sqrt(max((x1**2 - 1) * (1 - e1**2), 0.0))
                    rho2 = np.sqrt(max((x2**2 - 1) * (1 - e2**2), 0.0))
                    Aa = (x1 * e1 - x2 * e2)**2 + rho1**2 + rho2**2
                    Bc = 2 * rho1 * rho2
                    r12 = hR * np.sqrt(np.maximum(Aa - Bc * np.cos(dphi), 0.0))
                    wgt = wxi[a] * wxi[cc] * weta[b] * weta[d] * Jp1 * Jp2 * (hR**6) * ef * vne
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


def vee_quad(basis_p, R, alpha):
    """Direct quad <g_i g_j r12^(pi+pj-1)> (unsym single-term g). V_ee = 1/r12."""
    g = build_quadrature_grids(N_xi=26, N_eta=18, N_phi=24, xi_max=15.0)
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
                for d in range(len(eta)):
                    e2 = eta[d]; Jp2 = x2**2 - e2**2
                    rho1 = np.sqrt(max((x1**2 - 1) * (1 - e1**2), 0.0))
                    rho2 = np.sqrt(max((x2**2 - 1) * (1 - e2**2), 0.0))
                    A = (x1 * e1 - x2 * e2)**2 + rho1**2 + rho2**2
                    Bc = 2 * rho1 * rho2
                    r12 = hR * np.sqrt(np.maximum(A - Bc * np.cos(dphi), 1e-300))
                    wgt = wxi[a] * wxi[cc] * weta[b] * weta[d] * Jp1 * Jp2 * (hR**6) * ef
                    for i in range(n):
                        bi, pi = basis_p[i]
                        gi = x1**bi.j * x2**bi.k * e1**bi.l * e2**bi.m
                        for jjj in range(n):
                            bj, pj = basis_p[jjj]
                            gj = x1**bj.j * x2**bj.k * e1**bj.l * e2**bj.m
                            Pm1 = pi + pj - 1
                            rP = r12**Pm1 if Pm1 != 0 else np.ones_like(dphi)
                            V[i, jjj] += wgt * gi * gj * np.sum(rP * wphi) * 2 * np.pi
    return V


def _kin_grad_mu0(ja, la, jb, lb, alpha, A, dxi=0, deta=0):
    """mpf sigma (mu=0) one-electron kinetic gradient moment
    int[(xi^2-1) dg_i dg_j] eta-overlap + xi-overlap [(1-eta^2) dg_i dg_j],
    with extra xi^dxi eta^deta multiplicative shifts (for the A-weighting)."""
    a = mp.mpf(alpha)

    def mx(jx):
        t = pr._shift([-a], jx)
        if jx > 0:
            t = pr._pa(t, pr._shift([mp.mpf(jx)], jx - 1))
        return t

    def my(lx):
        return pr._shift([mp.mpf(lx)], lx - 1) if lx > 0 else [mp.mpf(0)]

    xip = pr._pm(pr._pm(mx(ja), mx(jb)), pr._xi2m1(1))
    etap = pr._pm(pr._pm(my(la), my(lb)), pr._meta2(1))
    px = pr._shift([mp.mpf(1)], ja + jb)
    py = pr._shift([mp.mpf(1)], la + lb)
    return (pr._mom_xi(pr._shift(xip, dxi), A) * pr._mom_eta(pr._shift(py, deta))
            + pr._mom_xi(pr._shift(px, dxi), A) * pr._mom_eta(pr._shift(etap, deta)))


def kinetic_p1p1_mpf(basis_p, R, alpha, dps=40):
    """mpf kinetic for a PURE p=1 basis via the even-collapse:
    T = (R/2)^2 * KG_A - 2 <g_i|g_j>_p0,   KG_A = A-weighted g-kinetic."""
    with mp.workdps(dps):
        A = ngm._mono_moments(mp.mpf(2.0 * alpha), 4 * max(max(b.j, b.k) for b, _ in basis_p) + 8)
        h6 = (mp.mpf(R) / 2) ** 6
        cc = 4 * mp.pi ** 2
        pref_T = (2 / mp.mpf(R) ** 2) * h6
        hR2 = (mp.mpf(R) / 2) ** 2
        n = len(basis_p)
        T = np.empty((n, n), object)
        for i in range(n):
            bi, _ = basis_p[i]
            for jj in range(i, n):
                bj, _ = basis_p[jj]
                ja, la, ka, ma = bi.j, bi.l, bi.k, bi.m
                jb, lb, kb, mb = bj.j, bj.l, bj.k, bj.m
                KG = mp.mpf(0)
                for (co, d1, e1, d2, e2) in A_TERMS:
                    k1 = _kin_grad_mu0(ja, la, jb, lb, alpha, A, d1, e1)
                    o2 = pr._ov(ka + kb + d2, ma + mb + e2, 0, A)
                    o1 = pr._ov(ja + jb + d1, la + lb + e1, 0, A)
                    k2 = _kin_grad_mu0(ka, ma, kb, mb, alpha, A, d2, e2)
                    KG += co * (k1 * o2 + o1 * k2)
                KG *= pref_T * cc
                OV = h6 * cc * pr._ov(ja + jb, la + lb, 0, A) * pr._ov(ka + kb, ma + mb, 0, A)
                T[i, jj] = T[jj, i] = hR2 * KG - 2 * OV
        return np.array([[float(T[i, j]) for j in range(n)] for i in range(n)])


# ==========================================================================
# p0 x p1 kinetic (ODD).  IBP identity:  <grad phi_u . grad phi_v> collapses to
#   T = -1/2 int r12 g_v (Lap g_u) dV        (u = p0 fn, v = p1 fn)
# Piece A (grad.grad r12) cancels half of Piece B (g_v grad_gu.grad_r12) exactly,
# leaving a PURE odd-r12^1 element with NO grad_r12 term.  The Laplacian NUMERATOR
# L[g_u] (metric denominator cancels the electron's Jacobian) is a polynomial, so
# this routes through the same A K0 - B K1 machinery as the validated odd V_ee.
# ==========================================================================
def _pder(poly):
    """d/dx of a low->high coefficient list."""
    return [k * poly[k] for k in range(1, len(poly))] or [mp.mpf(0)]


def _L_polys(j, l, alpha):
    """Single-electron Laplacian numerator  d_xi((xi^2-1) d_xi g) + d_eta((1-eta^2) d_eta g)
    for g = xi^j eta^l e^{-a xi}, returned as two (xi_poly, eta_poly) pieces
    (the exp factor e^{-a xi} is implicit; it pairs with g_v -> e^{-2a xi})."""
    a = mp.mpf(alpha)
    # xi piece: Q = (xi^2-1) mx(j);  xi_poly = Q' - a Q ;  eta_poly = eta^l
    mxj = pr._pa(pr._shift([-a], j), pr._shift([mp.mpf(j)], j - 1)) if j > 0 else [-a]
    Q = pr._pm(pr._xi2m1(1), mxj)
    xi_piece = pr._ps(_pder(Q), [a * c for c in Q])
    eta_mono = pr._shift([mp.mpf(1)], l)
    # eta piece: Qe = (1-eta^2) my(l);  eta_poly = Qe' ;  xi_poly = xi^j
    myl = pr._shift([mp.mpf(l)], l - 1) if l > 0 else [mp.mpf(0)]
    Qe = pr._pm(pr._meta2(1), myl)
    eta_piece = _pder(Qe)
    xi_mono = pr._shift([mp.mpf(1)], j)
    return [(xi_piece, eta_mono), (xi_mono, eta_piece)]


def _expand(x1, y1, x2, y2, sign):
    """Cartesian product of four coeff lists -> [(coef, p1, q1, p2, q2), ...]."""
    out = []
    for p1 in range(len(x1)):
        if x1[p1] == 0:
            continue
        for q1 in range(len(y1)):
            if y1[q1] == 0:
                continue
            for p2 in range(len(x2)):
                if x2[p2] == 0:
                    continue
                for q2 in range(len(y2)):
                    if y2[q2] == 0:
                        continue
                    out.append((sign * x1[p1] * y1[q1] * x2[p2] * y2[q2],
                                p1, q1, p2, q2))
    return out


def _kin_odd_base(bu, bv, elec, alpha):
    """Monomials of  g_v * L_elec[g_u] * (Jacobian of the OTHER electron).
    elec=1: L on electron 1 (Jac cancelled), electron 2 carries (xi2^2-eta2^2)."""
    monos = []
    if elec == 1:
        for (xip, etap) in _L_polys(bu.j, bu.l, alpha):
            x1 = pr._pm(pr._shift([mp.mpf(1)], bv.j), xip)
            y1 = pr._pm(pr._shift([mp.mpf(1)], bv.l), etap)
            K2, M2 = bv.k + bu.k, bv.m + bu.m
            for (jc, dP, dQ) in ((1, 2, 0), (-1, 0, 2)):   # (xi2^2 - eta2^2)
                x2 = pr._shift([mp.mpf(1)], K2 + dP)
                y2 = pr._shift([mp.mpf(1)], M2 + dQ)
                monos += _expand(x1, y1, x2, y2, jc)
    else:
        for (xip, etap) in _L_polys(bu.k, bu.m, alpha):
            x2 = pr._pm(pr._shift([mp.mpf(1)], bv.k), xip)
            y2 = pr._pm(pr._shift([mp.mpf(1)], bv.m), etap)
            K1, M1 = bv.j + bu.j, bv.l + bu.l
            for (jc, dP, dQ) in ((1, 2, 0), (-1, 0, 2)):   # (xi1^2 - eta1^2)
                x1 = pr._shift([mp.mpf(1)], K1 + dP)
                y1 = pr._shift([mp.mpf(1)], M1 + dQ)
                monos += _expand(x1, y1, x2, y2, jc)
    return monos


def _odd_K0(monos, Xtab, l_neu):
    """sum over l, terms of  (2l+1) coef X_l^{0,0}[P1][P2] c_l(Q1) c_l(Q2)."""
    tot = mp.mpf(0)
    for (co, P1, Q1, P2, Q2) in monos:
        for l in range(0, l_neu + 1):
            X = Xtab.get((l, 0, 0))
            if X is None or P1 >= len(X) or P2 >= len(X):
                continue
            c1 = _cl_m0(l, Q1)
            if c1 == 0:
                continue
            c2 = _cl_m0(l, Q2)
            if c2 == 0:
                continue
            tot += (2 * l + 1) * co * X[P1][P2] * c1 * c2
    return tot


def _odd_K1(monos, Xtab, cl1, l_neu):
    """B-folded m=1 sum:  2 cl1[l] coef X_l^{1,1}[P1][P2] d_l(Q1) d_l(Q2)."""
    tot = mp.mpf(0)
    for (co, P1, Q1, P2, Q2) in monos:
        for l in range(1, l_neu + 1):
            X = Xtab.get((l, 1, 1))
            if X is None or P1 >= len(X) or P2 >= len(X):
                continue
            d1 = _cl_m1(l, Q1)
            if d1 == 0:
                continue
            d2 = _cl_m1(l, Q2)
            if d2 == 0:
                continue
            tot += 2 * cl1[l] * co * X[P1][P2] * d1 * d2
    return tot


def _mulA(monos):
    """base x A_TERMS  (for the A K0 part)."""
    out = []
    for (co, P1, Q1, P2, Q2) in monos:
        for (ac, d1, e1, d2, e2) in A_TERMS:
            out.append((co * ac, P1 + d1, Q1 + e1, P2 + d2, Q2 + e2))
    return out


# ==========================================================================
# OPTIMIZED odd routing.  The odd-r12^1 value of a unit monomial base
# (1, P1, Q1, P2, Q2) is a FIXED scalar  G_odd[P1,Q1,P2,Q2] = (A K0) + (B K1):
#   G0 = sum_A sum_l (2l+1) X_l^{0,0}[P1+d1][P2+d2] c_l(Q1+e1) c_l(Q2+e2)
#   G1 =        sum_l 2 cl1[l] X_l^{1,1}[P1][P2]   d_l(Q1)     d_l(Q2)
# Precomputed ONCE (memoized) so every matrix element is O(base) lookups instead
# of O(A_TERMS * l * moments).  Reproduces _odd_K0(_mulA(.)) + _odd_K1(.) exactly.
# ==========================================================================
def _precompute_cl(l_neu, q_max):
    """cl0_tab[l][Q] = int eta^Q P_l ;  cl1_tab[l][Q] = int eta^Q (1-eta^2) P_l'."""
    cl0 = [dict() for _ in range(l_neu + 1)]
    cl1 = [dict() for _ in range(l_neu + 1)]
    for l in range(l_neu + 1):
        leg = list(ngm._leg_coeffs(l))
        rp = pr._pm(pr._meta2(1), list(ngm._RP_poly(l, 1))) if l >= 1 else None
        for Q in range(q_max + 1):
            cl0[l][Q] = pr._mom_eta(pr._pm(pr._shift([mp.mpf(1)], Q), leg))
            cl1[l][Q] = pr._mom_eta(pr._pm(pr._shift([mp.mpf(1)], Q), rp)) if l >= 1 else mp.mpf(0)
    return cl0, cl1


def _make_godd(Xtab, cl1v, cl0_tab, cl1_tab, l_neu):
    cache = {}
    X0 = [Xtab.get((l, 0, 0)) for l in range(l_neu + 1)]
    X1 = [Xtab.get((l, 1, 1)) for l in range(l_neu + 1)]

    def godd(P1, Q1, P2, Q2):
        key = (P1, Q1, P2, Q2)
        v = cache.get(key)
        if v is not None:
            return v
        tot = mp.mpf(0)
        # G0: A K0
        for (ac, d1, e1, d2, e2) in A_TERMS:
            PP1, PP2, QQ1, QQ2 = P1 + d1, P2 + d2, Q1 + e1, Q2 + e2
            for l in range(l_neu + 1):
                X = X0[l]
                if X is None or PP1 >= len(X) or PP2 >= len(X):
                    continue
                c1 = cl0_tab[l].get(QQ1)
                if not c1:
                    continue
                c2 = cl0_tab[l].get(QQ2)
                if not c2:
                    continue
                tot += ac * (2 * l + 1) * X[PP1][PP2] * c1 * c2
        # G1: B K1
        for l in range(1, l_neu + 1):
            X = X1[l]
            if X is None or P1 >= len(X) or P2 >= len(X):
                continue
            d1v = cl1_tab[l].get(Q1)
            if not d1v:
                continue
            d2v = cl1_tab[l].get(Q2)
            if not d2v:
                continue
            tot += 2 * cl1v[l] * X[P1][P2] * d1v * d2v
        cache[key] = tot
        return tot

    return godd


def _odd_value(monos, godd):
    """Odd-r12^1 value of a monomial-list base = sum coef * G_odd[powers]."""
    return sum((co * godd(P1, Q1, P2, Q2) for (co, P1, Q1, P2, Q2) in monos), mp.mpf(0))


def _odd_context(basis_p_or_g, alpha, l_neumann, p_extra=8):
    """Shared odd-routing context (Xtab + G_odd), built once per assembly.
    Accepts a list of ProductFn or of (ProductFn, p) pairs."""
    gs = [b[0] if isinstance(b, tuple) else b for b in basis_p_or_g]
    jmax = max(max(b.j, b.k) for b in gs)
    lmax = max(max(b.l, b.m) for b in gs)
    p_max = 2 * jmax + p_extra
    l_neu = min(l_neumann, 2 * lmax + p_extra)
    Xtab = pr._build_Xtab_mp([(0, 0), (1, 1)], l_neu, p_max, alpha,
                             {(0, 0): l_neu, (1, 1): l_neu})
    cl1v = [mp.mpf(0)] + [(2 * l + 1) * (mp.factorial(l - 1) / mp.factorial(l + 1)) ** 2
                          for l in range(1, l_neu + 1)]
    cl0_tab, cl1_tab = _precompute_cl(l_neu, 2 * lmax + p_extra + 2)
    godd = _make_godd(Xtab, cl1v, cl0_tab, cl1_tab, l_neu)
    return godd


def kinetic_mixed_mpf(basis_p, R, alpha, l_neumann=40, dps=40, to_float=True, godd=None):
    """Unified mpf kinetic for a MIXED p={0,1} basis.
      p0 x p0  -> plain prolate kinetic (even, via _kin)
      p1 x p1  -> even-collapse  (R/2)^2 KG_A - 2 <g|g>
      p0 x p1  -> odd  T = -1/2 int r12 g_v (Lap g_u)  via A K0 - B K1 (G_odd)."""
    with mp.workdps(dps):
        jmax = max(max(b.j, b.k) for b, _ in basis_p)
        A = ngm._mono_moments(mp.mpf(2.0 * alpha), 4 * jmax + 12)
        h6 = (mp.mpf(R) / 2) ** 6
        cc = 4 * mp.pi ** 2
        hR2 = (mp.mpf(R) / 2) ** 2
        pref_T = (2 / mp.mpf(R) ** 2) * h6         # p1p1 even-collapse KG
        pref_T00 = mp.mpf('0.5') * (4 / mp.mpf(R) ** 2) * h6   # p0p0 plain kinetic
        pref_odd = h6 * (2 * mp.pi) ** 2 * (2 / mp.mpf(R))     # p0p1 odd
        if godd is None:
            godd = _odd_context(basis_p, alpha, l_neumann)
        n = len(basis_p)
        T = np.empty((n, n), object)
        for i in range(n):
            bi, ppi = basis_p[i]
            for jj in range(i, n):
                bj, ppj = basis_p[jj]
                P = ppi + ppj
                if P == 0:
                    ja, la, ka, ma = bi.j, bi.l, bi.k, bi.m
                    jb, lb, kb, mb = bj.j, bj.l, bj.k, bj.m
                    k1, _ = pr._kin(ja, la, jb, lb, 0, alpha, A)
                    k2, _ = pr._kin(ka, ma, kb, mb, 0, alpha, A)
                    o1 = pr._ov(ja + jb, la + lb, 0, A)
                    o2 = pr._ov(ka + kb, ma + mb, 0, A)
                    val = pref_T00 * cc * (k1 * o2 + o1 * k2)
                elif P == 2:
                    ja, la, ka, ma = bi.j, bi.l, bi.k, bi.m
                    jb, lb, kb, mb = bj.j, bj.l, bj.k, bj.m
                    KG = mp.mpf(0)
                    for (co, d1, e1, d2, e2) in A_TERMS:
                        k1 = _kin_grad_mu0(ja, la, jb, lb, alpha, A, d1, e1)
                        o2 = pr._ov(ka + kb + d2, ma + mb + e2, 0, A)
                        o1 = pr._ov(ja + jb + d1, la + lb + e1, 0, A)
                        k2 = _kin_grad_mu0(ka, ma, kb, mb, alpha, A, d2, e2)
                        KG += co * (k1 * o2 + o1 * k2)
                    KG *= pref_T * cc
                    OV = h6 * cc * pr._ov(ja + jb, la + lb, 0, A) * pr._ov(ka + kb, ma + mb, 0, A)
                    val = hR2 * KG - 2 * OV
                else:   # P == 1, odd: u = p0 fn, v = p1 fn
                    (bu, bv) = (bi, bj) if ppi == 0 else (bj, bi)
                    base1 = _kin_odd_base(bu, bv, 1, alpha)
                    base2 = _kin_odd_base(bu, bv, 2, alpha)
                    tot = _odd_value(base1, godd) + _odd_value(base2, godd)
                    val = mp.mpf('-0.5') * pref_odd * tot
                T[i, jj] = T[jj, i] = val
        return T if not to_float else _tofloat(T)


def quad_kinetic(basis_p, R, alpha):
    """Unsym single-term prolate-Green kinetic for phi = g r12^p (reference)."""
    from geovac.hylleraas import _eval_unsym_and_derivs, _r12_and_derivs_python
    g = build_quadrature_grids(N_xi=24, N_eta=18, N_phi=28, xi_max=15.0)
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    dphi, wphi = g['dphi'], g['w_phi']
    hR = R / 2.0
    Tpref = R / 4.0
    n = len(basis_p)

    def derivs(bf, p, x1, e1, x2, e2, dp):
        r12, drx1, dre1, drx2, dre2, drdp = _r12_and_derivs_python(x1, e1, x2, e2, dp, R)
        r12 = max(r12, 1e-15)
        gg, gx1, ge1, gx2, ge2 = _eval_unsym_and_derivs(bf.j, bf.k, bf.l, bf.m, alpha, x1, e1, x2, e2)
        if p == 0:
            return gg, gx1, ge1, gx2, ge2, 0.0
        rp = r12 ** p
        coup = 0.5 * p * r12 ** (p - 2)
        gc = gg * coup
        return (gg * rp, gx1 * rp + gc * drx1, ge1 * rp + gc * dre1,
                gx2 * rp + gc * drx2, ge2 * rp + gc * dre2, gc * drdp)

    T = np.zeros((n, n))
    for a in range(len(xi)):
        x1 = xi[a]; x1m1 = x1**2 - 1
        for c in range(len(xi)):
            x2 = xi[c]; x2m1 = x2**2 - 1; wac = wxi[a] * wxi[c]
            for b in range(len(eta)):
                e1 = eta[b]; om1 = 1 - e1**2; J1 = hR**3 * (x1**2 - e1**2)
                for d in range(len(eta)):
                    e2 = eta[d]; om2 = 1 - e2**2; J2 = hR**3 * (x2**2 - e2**2)
                    den1, den2 = x1m1 * om1, x2m1 * om2
                    pf1 = (x1**2 - e1**2) / den1 if den1 > 1e-30 else 0.0
                    pf2 = (x2**2 - e2**2) / den2 if den2 > 1e-30 else 0.0
                    acc = np.zeros((n, n))
                    for ip in range(len(dphi)):
                        D = [derivs(bf, p, x1, e1, x2, e2, dphi[ip]) for bf, p in basis_p]
                        for i in range(n):
                            _, ix1, ie1, ix2, ie2, idp = D[i]
                            for jj in range(n):
                                _, jx1, je1, jx2, je2, jdp = D[jj]
                                t1 = (x1m1 * ix1 * jx1 + om1 * ie1 * je1 + pf1 * idp * jdp) * J2
                                t2 = (x2m1 * ix2 * jx2 + om2 * ie2 * je2 + pf2 * idp * jdp) * J1
                                acc[i, jj] += (t1 + t2) * wphi[ip]
                    T += wac * weta[b] * weta[d] * Tpref * 2 * np.pi * acc
    return T


# ==========================================================================
# ASSEMBLY: full mixed-p (S, H) with H = T + V_ne + V_ee.
# Per entry the r12 power is dispatched by P = p_i + p_j:
#   S:     r12^P            (P in {0,1,2})
#   V_ee:  r12^(P-1)        (1/r12 for P=0; r12^0 for P=1; r12^1 for P=2)
#   V_ne:  vne_mpf (handles P internally)
#   T:     kinetic_mixed_mpf (handles P internally)
# ==========================================================================
def _tofloat(M):
    n = M.shape[0]
    return np.array([[float(M[i, j]) for j in range(n)] for i in range(n)])


def _kern_even(basis_g, R, alpha, power, dps=40):
    """<g_i g_j r12^power> (fullJac) for power in {0, 2}, ALL pairs -> mpf array."""
    with mp.workdps(dps):
        A = ngm._mono_moments(mp.mpf(2.0 * alpha),
                              4 * max(max(b.j, b.k) for b in basis_g) + 8)
        h6 = (mp.mpf(R) / 2) ** 6
        hR2 = (mp.mpf(R) / 2) ** 2
        cc = 4 * mp.pi ** 2
        n = len(basis_g)
        M = np.empty((n, n), object)
        for i in range(n):
            bi = basis_g[i]
            for jj in range(i, n):
                bj = basis_g[jj]
                j1, k1, l1, m1 = bi.j + bj.j, bi.k + bj.k, bi.l + bj.l, bi.m + bj.m
                if power == 0:
                    val = h6 * cc * pr._ov(j1, l1, 0, A) * pr._ov(k1, m1, 0, A)
                else:  # power == 2
                    acc = mp.mpf(0)
                    for (co, d1, e1, d2, e2) in A_TERMS:
                        acc += co * pr._ov(j1 + d1, l1 + e1, 0, A) * pr._ov(k1 + d2, m1 + e2, 0, A)
                    val = h6 * hR2 * cc * acc
                M[i, jj] = M[jj, i] = val
        return M


def _kern_odd1(basis_g, R, alpha, l_neumann=40, dps=40, godd=None):
    """<g_i g_j r12^1> (odd) for ALL pairs, as mpf objects (p ignored)."""
    with mp.workdps(dps):
        if godd is None:
            godd = _odd_context(basis_g, alpha, l_neumann)
        h6 = (mp.mpf(R) / 2) ** 6
        pref = h6 * (mp.mpf(R) / 2) ** 2 * (2 * mp.pi) ** 2 * (2 / mp.mpf(R))
        n = len(basis_g)
        M = np.empty((n, n), object)
        for i in range(n):
            bi = basis_g[i]
            for jj in range(i, n):
                bj = basis_g[jj]
                p1x, q1 = bi.j + bj.j, bi.l + bj.l
                p2x, q2 = bi.k + bj.k, bi.m + bj.m
                jac = [(1, p1x + 2, q1, p2x + 2, q2), (-1, p1x + 2, q1, p2x, q2 + 2),
                       (-1, p1x, q1 + 2, p2x + 2, q2), (1, p1x, q1 + 2, p2x, q2 + 2)]
                M[i, jj] = M[jj, i] = pref * _odd_value(jac, godd)
        return M


def assemble_mixed(basis_p, R, alpha, l_neumann=40, dps=40, mpf_out=False):
    """Full (S, H) for a mixed p={0,1} basis.  H = T + V_ne + V_ee.
    mpf_out=False -> float64 ndarrays;  True -> mpf object ndarrays (for the
    high-precision conditioned solve that breaks the float64 downcast wall)."""
    with mp.workdps(dps):
        gs = [b for b, _ in basis_p]
        ps = [p for _, p in basis_p]
        n = len(basis_p)
        godd = _odd_context(basis_p, alpha, l_neumann)       # built ONCE, shared
        Kev0 = _kern_even(gs, R, alpha, 0, dps)              # mpf
        Kev2 = _kern_even(gs, R, alpha, 2, dps)              # mpf
        Kod1 = _kern_odd1(gs, R, alpha, l_neumann, dps, godd=godd)          # mpf
        Kinv = pr.vee_mp(gs, alpha, R, l_neumann)            # mpf (<g g / r12>)
        Tm = kinetic_mixed_mpf(basis_p, R, alpha, l_neumann, dps, to_float=False, godd=godd)
        Vne = vne_mpf(basis_p, R, alpha, l_neumann, dps, to_float=False, godd=godd)
        S = np.empty((n, n), object)
        Vee = np.empty((n, n), object)
        for i in range(n):
            for jj in range(n):
                P = ps[i] + ps[jj]
                S[i, jj] = Kev0[i, jj] if P == 0 else (Kod1[i, jj] if P == 1 else Kev2[i, jj])
                if P == 0:
                    Vee[i, jj] = Kinv[i, jj]        # standard 1/r12
                elif P == 1:
                    Vee[i, jj] = Kev0[i, jj]         # r12^0
                else:
                    Vee[i, jj] = Kod1[i, jj]         # r12^1
        H = Tm + Vne + Vee
        if mpf_out:
            return S, H
        return _tofloat(S), _tofloat(H)


def solve_canonical_mpf(S, H, dps=50, thresholds=(1e-20, 1e-18, 1e-16, 1e-14, 1e-12)):
    """mpf canonical-orthogonalization eigensolve of the generalized problem
    (S, H are mpf object ndarrays).  Diagonalize S in mpf, drop directions below
    tol*max, form H_o = X^T H X in the retained subspace, mp.eigsy -> lowest E.
    Returns (E_min, tol, kept).  This is the wall-breaker: the float64 downcast
    of an ill-conditioned S loses the correlation energy; mpf does not."""
    with mp.workdps(dps):
        n = S.shape[0]
        Smp = mp.matrix([[S[i, j] for j in range(n)] for i in range(n)])
        Hmp = mp.matrix([[H[i, j] for j in range(n)] for i in range(n)])
        sval, svec = mp.eigsy(Smp)          # ascending eigenvalues, orthonormal vecs
        smax = sval[n - 1]
        best = None
        for tol in thresholds:
            cut = tol * smax
            idx = [k for k in range(n) if sval[k] > cut]
            if not idx:
                continue
            nk = len(idx)
            # X = svec[:, idx] * diag(1/sqrt(sval[idx]))
            X = mp.matrix(n, nk)
            for cc, k in enumerate(idx):
                inv = 1 / mp.sqrt(sval[k])
                for r in range(n):
                    X[r, cc] = svec[r, k] * inv
            Ho = X.T * Hmp * X
            # Ho is the orthonormalized Hamiltonian -> WELL conditioned, so the
            # final eigenvalue is safe in fast float64 (the mpf work that breaks
            # the wall is the S-orthogonalization above, not this solve).
            Hof = np.array([[float(Ho[i, j]) for j in range(nk)] for i in range(nk)])
            Hof = 0.5 * (Hof + Hof.T)
            e = float(np.linalg.eigvalsh(Hof)[0])
            if best is None or e < best[0]:
                best = (e, tol, nk)
        return best


def unsym_overlap_quad(basis_p, R, alpha):
    """Direct 5D quadrature of <g_i g_j r12^(pi+pj)> for UNsymmetrized single-term
    g (matching the ProductFn convention)."""
    g = build_quadrature_grids(N_xi=26, N_eta=18, N_phi=16, xi_max=15.0)
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    dphi, wphi = g['dphi'], g['w_phi']
    hR = R / 2.0
    n = len(basis_p)
    S = np.zeros((n, n))
    for a in range(len(xi)):
        x1 = xi[a]
        for cc in range(len(xi)):
            x2 = xi[cc]
            ef = np.exp(-2.0 * alpha * (x1 + x2))   # g_i g_j: each g has e^{-a(xi1+xi2)}
            for b in range(len(eta)):
                e1 = eta[b]; Jp1 = x1**2 - e1**2
                for d in range(len(eta)):
                    e2 = eta[d]; Jp2 = x2**2 - e2**2
                    rho1 = np.sqrt(max((x1**2-1)*(1-e1**2), 0.0))
                    rho2 = np.sqrt(max((x2**2-1)*(1-e2**2), 0.0))
                    A = (x1*e1 - x2*e2)**2 + rho1**2 + rho2**2
                    Bc = 2*rho1*rho2
                    r12 = hR*np.sqrt(np.maximum(A - Bc*np.cos(dphi), 0.0))
                    wgt = wxi[a]*wxi[cc]*weta[b]*weta[d]*Jp1*Jp2*(hR**6)*ef
                    for i in range(n):
                        bi, pi = basis_p[i]
                        gi = x1**bi.j * x2**bi.k * e1**bi.l * e2**bi.m
                        for jjj in range(n):
                            bj, pj = basis_p[jjj]
                            gj = x1**bj.j * x2**bj.k * e1**bj.l * e2**bj.m
                            rP = r12**(pi+pj) if (pi+pj) > 0 else np.ones_like(dphi)
                            S[i, jjj] += wgt * gi * gj * np.sum(rP * wphi) * 2*np.pi
    return S


def main():
    specs = [(0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 1, 1), (1, 1, 0, 0), (2, 0, 0, 0)]
    basis = [(pr.ProductFn(j, l, k, m, 0, ALPHA), 1) for (j, l, k, m) in specs]
    S_mpf = overlap_even_mpf(basis, R, ALPHA)
    S_q = unsym_overlap_quad(basis, R, ALPHA)   # unsymmetrized reference

    rel = np.max(np.abs(S_mpf - S_q) / np.maximum(np.abs(S_q), 1e-12))
    print("mpf even overlap (r12^2) diag:", np.round(np.diag(S_mpf), 6))
    print("unsym quad overlap     diag:", np.round(np.diag(S_q), 6))
    print(f"max rel diff = {rel:.3e}   (quad xi,eta-grid-limited; expect ~1e-3)")
    print("MPF EVEN OVERLAP", "VALIDATED" if rel < 5e-3 else "MISMATCH")

    # --- odd V_ee (r12^1) for the pure-p1 basis ---
    print()
    Vmpf = vee_r12_odd_mpf(basis, R, ALPHA)
    Vq = vee_quad(basis, R, ALPHA)
    relV = np.max(np.abs(Vmpf - Vq) / np.maximum(np.abs(Vq), 1e-12))
    print("mpf odd V_ee (r12^1) diag:", np.round(np.diag(Vmpf), 6))
    print("quad     V_ee (r12^1) diag:", np.round(np.diag(Vq), 6))
    print(f"max rel diff = {relV:.3e}   (quad xi,eta-grid-limited; expect ~1e-3)")
    print("MPF ODD V_ee", "VALIDATED" if relV < 5e-3 else "MISMATCH")

    # --- V_ne on a MIXED p={0,1} basis (tests p0p0, p1p1 even + p0p1 odd) ---
    print()
    mspec = [(0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 1, 1)]
    mixed = ([(pr.ProductFn(j, l, k, m, 0, ALPHA), 0) for (j, l, k, m) in mspec]
             + [(pr.ProductFn(j, l, k, m, 0, ALPHA), 1) for (j, l, k, m) in mspec])
    Vne_m = vne_mpf(mixed, R, ALPHA)
    Vne_q = vne_quad(mixed, R, ALPHA)
    relVne = np.max(np.abs(Vne_m - Vne_q) / np.maximum(np.abs(Vne_q), 1e-10))
    print("mpf V_ne (mixed) diag:", np.round(np.diag(Vne_m), 5))
    print("quad V_ne (mixed) diag:", np.round(np.diag(Vne_q), 5))
    print(f"max rel diff = {relVne:.3e}   (quad xi,eta-grid-limited; expect ~1e-3)")
    print("MPF V_ne (mixed p0p0/p1p1/p0p1)", "VALIDATED" if relVne < 5e-3 else "MISMATCH")

    # --- kinetic p1p1 (even-collapse) vs single-term prolate-Green quad ---
    print()
    Tm = kinetic_p1p1_mpf(basis, R, ALPHA)
    Tq = quad_kinetic(basis, R, ALPHA)
    relT = np.max(np.abs(Tm - Tq) / np.maximum(np.abs(Tq), 1e-10))
    print("mpf kinetic p1p1 diag:", np.round(np.diag(Tm), 5))
    print("quad kinetic p1p1 diag:", np.round(np.diag(Tq), 5))
    print(f"max rel diff = {relT:.3e}   (quad xi,eta-grid-limited; expect ~1e-3)")
    print("MPF KINETIC p1p1 (even-collapse)", "VALIDATED" if relT < 5e-3 else "MISMATCH")

    # --- unified MIXED kinetic (tests the p0xp1 ODD block) vs quad ---
    print()
    mspec2 = [(0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 1, 1)]
    mixed2 = ([(pr.ProductFn(j, l, k, m, 0, ALPHA), 0) for (j, l, k, m) in mspec2]
              + [(pr.ProductFn(j, l, k, m, 0, ALPHA), 1) for (j, l, k, m) in mspec2])
    Tmm = kinetic_mixed_mpf(mixed2, R, ALPHA)
    Tqm = quad_kinetic(mixed2, R, ALPHA)
    # focus on the p0xp1 off-diagonal block (rows 0-2 = p0, cols 3-5 = p1)
    off_m = Tmm[:3, 3:]
    off_q = Tqm[:3, 3:]
    relTo = np.max(np.abs(off_m - off_q) / np.maximum(np.abs(off_q), 1e-10))
    relTfull = np.max(np.abs(Tmm - Tqm) / np.maximum(np.abs(Tqm), 1e-10))
    print("mpf  p0xp1 kinetic block:\n", np.round(off_m, 5))
    print("quad p0xp1 kinetic block:\n", np.round(off_q, 5))
    print(f"p0xp1 block max rel diff = {relTo:.3e}   (quad grid-limited; expect ~1e-3)")
    print(f"full mixed matrix max rel diff = {relTfull:.3e}")
    print("MPF KINETIC p0xp1 (odd, IBP-collapsed)",
          "VALIDATED" if relTo < 5e-3 else "MISMATCH")


if __name__ == "__main__":
    main()
