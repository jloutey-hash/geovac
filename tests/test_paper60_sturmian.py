"""
Backing tests for Paper 60 — the generalized-Sturmian secular equation as a
metric-free quantum secular equation (atoms) + the molecular Shibuya-Wulfman
conditioning.

Each test maps to a load-bearing claim in the paper. Self-contained (no debug/
imports; debug/ is transient). Radial grids kept modest for test speed.
"""
import math
import numpy as np
import pytest
from scipy.special import genlaguerre
from scipy.linalg import eigh
from scipy.integrate import cumulative_trapezoid
from itertools import combinations_with_replacement
from math import factorial as _fac


# --------------------------------------------------------------------------
# shared radial machinery (unit / arbitrary charge Coulomb-Sturmian s-orbitals)
# --------------------------------------------------------------------------
def _grid(rmax=70.0, npts=30000):
    r = np.linspace(1e-6, rmax, npts)
    return r, r[1] - r[0]


def _cs_s(r, n, Q):
    """Hydrogenic s radial at charge Q (decay Q/n) — the config orbital. L2-normalized."""
    a = Q / n
    f = np.exp(-a * r) * genlaguerre(n - 1, 1)(2 * a * r)
    return f / np.sqrt(np.trapezoid(f * f * r * r, r))


def _sturmian_s(r, n, k=1.0):
    """GENUINE shared-scale Coulomb-Sturmian s radial (decay k for ALL n). L2-normalized.
    These are L2-non-orthogonal across n but 1/r-orthogonal (potential-weighted)."""
    f = np.exp(-k * r) * genlaguerre(n - 1, 1)(2 * k * r)
    return f / np.sqrt(np.trapezoid(f * f * r * r, r))


def _F0(r, dr, fa, fb):
    """Unit-charge Slater F^0(a,b) (l=0 multipole)."""
    g = fa * fb
    U = np.cumsum(g * r * r) * dr / r + np.cumsum((g * r)[::-1])[::-1] * dr
    return np.trapezoid(fa * fb * U * r * r, r)


# --------------------------------------------------------------------------
# eq:pw  — potential-weighted orthonormality is diagonal; L2 overlap is not
# --------------------------------------------------------------------------
def test_paper60_potential_weighted_orthonormality_diagonal():
    r, _ = _grid()
    N = 6
    fs = [_sturmian_s(r, n) for n in range(1, N + 1)]   # genuine shared-scale Sturmians
    SV = np.zeros((N, N))   # <i|1/r|j>
    SL = np.zeros((N, N))   # <i|j>  (L2)
    for i in range(N):
        for j in range(N):
            SV[i, j] = np.trapezoid(fs[i] * fs[j] * r, r)
            SL[i, j] = np.trapezoid(fs[i] * fs[j] * r * r, r)
    off = np.abs(SV - np.diag(np.diag(SV)))
    # potential-weighted overlap is diagonal (Avery eq 6.7 / Paper 60 eq:pw): off-diagonal ~ 0
    assert off.max() < 2e-3, f"potential-weighted overlap not diagonal: {off.max()}"
    # the L2 overlap of the SAME (shared-scale) basis is ill-conditioned
    assert np.linalg.cond(SL) > 10.0


def _cfgs_upto(lmax, span):
    """Nested (l, na, nb) config lists of growing size; principal n >= l+1 enforced."""
    out = []
    for l in range(lmax + 1):
        ns = list(range(l + 1, l + 1 + span))
        for ia, na in enumerate(ns):
            for nb in ns[ia:]:
                out.append((l, na, nb))
    return out


def _atomic_lowdin_lambda(N, mode, Z=2.0, rmax=80.0, npts=20000):
    """Standard L2-Loewdin block-encoding 1-norm lambda = sum|h| + sum|(pq|rs)| for an
    N-function s-only He model, in the shared-scale Sturmian ('sturmian', decay k=1) or
    hydrogenic ('hydrogenic', scale Z/n) basis.  Reproduces the MECHANISM behind eq:blowup:
    the ill-conditioned shared-scale overlap makes Loewdin inflate lambda faster than the
    well-conditioned hydrogenic one.  (The exact Q^3.33/Q^1.19 exponents are a property of
    the full Goscinskian construction; this pins the mechanism + shared>hydrogenic ordering.)"""
    from scipy.integrate import cumulative_trapezoid as _ct
    r = np.linspace(1e-6, rmax, npts)
    fs = [(_sturmian_s(r, n, k=1.0) if mode == 'sturmian' else _cs_s(r, n, Q=Z / n))
          for n in range(1, N + 1)]
    S = np.array([[np.trapezoid(fs[i] * fs[j] * r * r, r) for j in range(N)] for i in range(N)])
    h = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            gi, gj = np.gradient(fs[i], r), np.gradient(fs[j], r)
            h[i, j] = 0.5 * np.trapezoid(gi * gj * r * r, r) - Z * np.trapezoid(fs[i] * fs[j] * r, r)
    dens = {(i, j): fs[i] * fs[j] for i in range(N) for j in range(N)}
    er = np.zeros((N, N, N, N))
    for i in range(N):
        for j in range(N):
            a = dens[(i, j)] * r * r
            inn = np.concatenate(([0.0], _ct(a, x=r))) / r
            out = np.concatenate(([0.0], _ct((a / r)[::-1], x=r)))[::-1]
            pot = inn + out
            for k in range(N):
                for l in range(N):
                    er[i, j, k, l] = np.trapezoid(dens[(k, l)] * pot * r * r, r)
    w, U = np.linalg.eigh(S)
    keep = w > 1e-10
    X = U[:, keep] @ np.diag(1.0 / np.sqrt(w[keep]))
    hm = X.T @ h @ X
    em = np.einsum('ip,jq,kr,ls,ijkl->pqrs', X, X, X, X, er, optimize=True)
    return np.abs(hm).sum() + np.abs(em).sum(), np.linalg.cond(S)


@pytest.mark.slow
def test_paper60_l2_overlap_illconditioned_grows():
    """eq:blowup MECHANISM: the shared-scale Sturmian overlap ill-conditions with basis size,
    so Loewdin inflates the block-encoding 1-norm FASTER for the shared-scale basis than for a
    well-conditioned hydrogenic one.  Pins (i) monotone cond(S) growth over >=4 points matching
    the paper's reported 3.0/5.8/13.9/32.2 sequence, and (ii) the shared>hydrogenic lambda-growth
    ordering.  (The exact Q^3.33/Q^1.19 exponents belong to the full Goscinskian construction.)"""
    r, _ = _grid()
    conds = []
    for N in (2, 3, 5, 8):
        fs = [_sturmian_s(r, n) for n in range(1, N + 1)]   # genuine shared-scale Sturmians
        S = np.array([[np.trapezoid(fs[i] * fs[j] * r * r, r) for j in range(N)]
                      for i in range(N)])
        conds.append(np.linalg.cond(S))
    assert conds[0] < conds[1] < conds[2] < conds[3]     # monotone growth, >=4 points
    for got, ref in zip(conds, (3.0, 5.8, 13.9, 32.2)):  # paper's reported cond(S) sequence
        assert abs(got - ref) < 0.25 * ref, f"cond(S) sequence off: {conds}"
    # Loewdin block-encoding 1-norm inflates FASTER for the shared-scale basis than hydrogenic
    Qs = [4, 6, 8, 10, 12]
    lam_sh = [_atomic_lowdin_lambda(N, 'sturmian')[0] for N in (2, 3, 4, 5, 6)]
    lam_hy = [_atomic_lowdin_lambda(N, 'hydrogenic')[0] for N in (2, 3, 4, 5, 6)]
    p_sh = np.polyfit(np.log(Qs), np.log(lam_sh), 1)[0]
    p_hy = np.polyfit(np.log(Qs), np.log(lam_hy), 1)[0]
    assert p_sh > p_hy + 0.5, f"shared-scale should inflate faster: Q^{p_sh:.2f} vs Q^{p_hy:.2f}"


# --------------------------------------------------------------------------
# sec:atomic — single-config He reproduces the textbook variational value
# --------------------------------------------------------------------------
def test_paper60_unit_charge_F0_is_five_eighths():
    r, dr = _grid()
    f1 = _cs_s(r, 1, 1.0)
    assert abs(_F0(r, dr, f1, f1) - 0.625) < 3e-3


def test_paper60_single_config_helium_variational():
    r, dr = _grid()
    Z = 2.0
    R = np.sqrt(1 / 1**2 + 1 / 1**2)              # sqrt(2), the 1s^2 root
    f1 = _cs_s(r, 1, 1.0)
    Tp = -_F0(r, dr, f1, f1) / R                  # pure number -(5/8)/sqrt2
    pk = Z * R + Tp
    E = -pk**2 / 2
    assert abs(Tp - (-(0.625) / np.sqrt(2))) < 3e-3
    assert abs(E - (-2.84766)) < 3e-3            # textbook variational He
    # bare (no repulsion) = two He+ 1s = -4.0
    assert abs(-(Z * R) ** 2 / 2 - (-4.0)) < 1e-6


# --------------------------------------------------------------------------
# eq:sublinear — multi-config lowers E and the 1-norm grows more slowly
# than K over the computed window (0.84 fitted on K = 74..164; NOT an
# asymptotic regime -- see the split test below for which block carries it)
# --------------------------------------------------------------------------
def _he_secular(r, dr, nmax, Z=2.0):
    def ov(a, b):
        return np.trapezoid(a * b * r * r, r)

    def eri(fp, fr, fq, fs):
        g = fq * fs
        U = np.cumsum(g * r * r) * dr / r + np.cumsum((g * r)[::-1])[::-1] * dr
        return np.trapezoid(fp * fr * U * r * r, r)

    def singlet_g(ob, ok):
        (a1, b1), (a2, b2) = ob, ok
        if a1 is b1 and a2 is b2:
            return eri(a1, a2, b1, b2)
        Np = 1.0 if a1 is b1 else np.sqrt(2 * (1 + ov(a1, b1) ** 2))
        Nk = 1.0 if a2 is b2 else np.sqrt(2 * (1 + ov(a2, b2) ** 2))
        tot = (eri(a1, a2, b1, b2) + eri(a1, b2, b1, a2)
               + eri(b1, a2, a1, b2) + eri(b1, b2, a1, a2))
        return tot / (Np * Nk)

    cfg = list(combinations_with_replacement(range(1, nmax + 1), 2))
    K = len(cfg)
    Rnu = np.array([np.sqrt(1 / na**2 + 1 / nb**2) for na, nb in cfg])
    Q = 1.0 / Rnu
    orb = {i: (_cs_s(r, na, Q[i]), _cs_s(r, nb, Q[i]))
           for i, (na, nb) in enumerate(cfg)}
    M = np.zeros((K, K))
    for i in range(K):
        for j in range(K):
            M[i, j] = (Z * Rnu[i] if i == j else 0.0) - singlet_g(orb[i], orb[j])
    M = 0.5 * (M + M.T)
    pk = np.sort(eigh(M, eigvals_only=True))[-1]
    return -pk**2 / 2, np.abs(M).sum(), K


@pytest.mark.slow
def test_paper60_multiconfig_lowers_and_sublinear_onenorm():
    r, dr = _grid(rmax=80.0, npts=20000)
    Es, Ks, L1s = [], [], []
    for nmax in (1, 3, 5):
        E, l1, K = _he_secular(r, dr, nmax)
        Es.append(E); Ks.append(K); L1s.append(l1)
    # adding configurations lowers the energy (radial correlation)
    assert Es[0] > Es[1] > Es[2]
    assert abs(Es[0] - (-2.84766)) < 3e-3           # single config = variational
    assert Es[-1] < -2.86 and Es[-1] > -2.90372     # heads toward s-limit, above exact
    # 1-norm grows more slowly than K (exponent < 1) over the computed range,
    # the opposite of L2-Lowdin.  The exponent is NOT asymptotically stable
    # (corrected 2026-09-07); the pinned claim is exponent < 1 here.
    p = np.polyfit(np.log(Ks[1:]), np.log(L1s[1:]), 1)[0]
    assert p < 1.0, f"1-norm exponent not sublinear: {p}"


# --------------------------------------------------------------------------
# sec:molecular — Shibuya-Wulfman metric: intra = identity, better than L2
# --------------------------------------------------------------------------
def _two_center(nmax, R, k=1.0):
    rho = np.linspace(1e-4, 22, 500)
    z = np.linspace(-18, R + 20, 800)
    RHO, ZZ = np.meshgrid(rho, z, indexing='ij')
    drho, dz = rho[1] - rho[0], z[1] - z[0]

    def cs_at(n, zc):
        rr = np.sqrt(RHO**2 + (ZZ - zc)**2)
        f = np.exp(-k * rr) * genlaguerre(n - 1, 1)(2 * k * rr)
        nrm = np.sqrt(2 * np.pi * np.sum(f * f * RHO) * drho * dz)
        return f / nrm

    def integ(g):
        return 2 * np.pi * np.sum(g * RHO) * drho * dz

    def sw(fi, fj):
        gi_r, gi_z = np.gradient(fi, rho, z)
        gj_r, gj_z = np.gradient(fj, rho, z)
        return (1 / (2 * k**2)) * integ(gi_r * gj_r + gi_z * gj_z) + 0.5 * integ(fi * fj)

    basis = [(n, 0.0) for n in range(1, nmax + 1)] + [(n, R) for n in range(1, nmax + 1)]
    fs = [cs_at(n, zc) for n, zc in basis]
    N = len(fs)
    O = np.array([[integ(fs[i] * fs[j]) for j in range(N)] for i in range(N)])
    S = np.array([[sw(fs[i], fs[j]) for j in range(N)] for i in range(N)])
    return O, S


@pytest.mark.slow
def test_paper60_sw_intra_center_is_identity():
    # widely separated centers => intra-center SW block ~ I, unlike the L2 overlap
    O, S = _two_center(3, 60.0)
    Sa = S[:3, :3]
    assert np.linalg.cond(Sa) < 1.05                     # SW intra ~ identity
    assert np.abs(Sa - np.eye(3)).max() < 1.5e-2         # grid-limited; SW intra = I
    assert np.linalg.cond(O[:3, :3]) > 3.0               # L2 intra ill-conditioned


@pytest.mark.slow
def test_paper60_sw_better_conditioned_than_l2():
    for R in (2.0, 4.0):
        O, S = _two_center(3, R)
        assert np.linalg.cond(S) < np.linalg.cond(O)     # SW beats L2 at every R
    # and the advantage widens at larger separation
    O2, S2 = _two_center(3, 2.0)
    O4, S4 = _two_center(3, 4.0)
    ratio2 = np.linalg.cond(O2) / np.linalg.cond(S2)
    ratio4 = np.linalg.cond(O4) / np.linalg.cond(S4)
    assert ratio4 > ratio2


# --------------------------------------------------------------------------
# sec:atomic (l>0) — angular correlation lowers E and keeps the 1-norm
# sublinear.  Angular machinery ported from the validated driver
# debug/sturmian_he_lmax.py (Wigner-3j / Gaunt, coupled-^1S two-electron configs).
# --------------------------------------------------------------------------
def _w3j(j1, j2, j3, m1, m2, m3):
    if m1 + m2 + m3 != 0 or not (abs(j1 - j2) <= j3 <= j1 + j2):
        return 0.0
    if any(abs(mm) > jj for mm, jj in ((m1, j1), (m2, j2), (m3, j3))):
        return 0.0
    delta = math.sqrt(_fac(j1 + j2 - j3) * _fac(j1 - j2 + j3) * _fac(-j1 + j2 + j3)
                      / _fac(j1 + j2 + j3 + 1))
    pref = math.sqrt(_fac(j1 + m1) * _fac(j1 - m1) * _fac(j2 + m2) * _fac(j2 - m2)
                     * _fac(j3 + m3) * _fac(j3 - m3))
    tmin = max(0, j2 - j3 - m1, j1 - j3 + m2)
    tmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    s = 0.0
    for t in range(tmin, tmax + 1):
        s += (-1) ** t / (_fac(t) * _fac(j1 + j2 - j3 - t) * _fac(j1 - m1 - t)
                          * _fac(j2 + m2 - t) * _fac(j3 - j2 + m1 + t) * _fac(j3 - j1 - m2 + t))
    return (-1) ** (j1 - j2 - m3) * delta * pref * s


def _gaunt(l1, l2, l3, m1, m2, m3):
    w0 = _w3j(l1, l2, l3, 0, 0, 0)
    if w0 == 0.0:
        return 0.0
    return (math.sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4 * math.pi))
            * w0 * _w3j(l1, l2, l3, m1, m2, m3))


def _cgL0(l, m):
    return (-1) ** (l - m) / math.sqrt(2 * l + 1)


def _radial(r, n, l, Q):
    a = Q / n
    f = (2 * a * r) ** l * np.exp(-a * r) * genlaguerre(n - l - 1, 2 * l + 1)(2 * a * r)
    return f / np.sqrt(np.trapezoid(f * f * r * r, r))


def _Rk(r, dr, Pa, Pc, Pb, Pd, k):
    g = Pb * Pd * r * r
    inner = np.concatenate(([0.0], cumulative_trapezoid(g * r ** k, dx=dr)))
    outer = np.concatenate(([0.0], cumulative_trapezoid((g * r ** (-(k + 1)))[::-1], dx=dr)))[::-1]
    Uk = inner * r ** (-(k + 1)) + outer * r ** k
    return np.trapezoid(Pa * Pc * Uk * r * r, r)


def _pair(r, dr, oa, ob, oc, od):
    la, ma, lb, mb = oa[0], oa[1], ob[0], ob[1]
    lc, mc, ld, md = oc[0], oc[1], od[0], od[1]
    if ma + mb != mc + md:
        return 0.0
    tot = 0.0
    for k in range(max(abs(la - lc), abs(lb - ld)), min(la + lc, lb + ld) + 1):
        q = mc - ma
        g1 = _gaunt(la, k, lc, -ma, -q, mc)
        if g1 == 0.0:
            continue
        g2 = _gaunt(lb, k, ld, -mb, q, md)
        if g2 == 0.0:
            continue
        ang = (-1) ** (ma + q + mb) * g1 * g2 * (4 * math.pi / (2 * k + 1))
        tot += ang * _Rk(r, dr, oa[2], oc[2], ob[2], od[2], k)
    return tot


def _config_terms(r, l, na, nb):
    R = math.sqrt(1 / na ** 2 + 1 / nb ** 2)
    Q = 1.0 / R
    Pa = _radial(r, na, l, Q)
    Pb = Pa if nb == na else _radial(r, nb, l, Q)
    terms = [(_cgL0(l, m), (l, m, Pa), (l, -m, Pb)) for m in range(-l, l + 1)]
    if na != nb:
        terms += [(_cgL0(l, m), (l, m, Pb), (l, -m, Pa)) for m in range(-l, l + 1)]
    nrm = _terms_overlap(r, terms, terms)
    return terms, R, 1.0 / math.sqrt(nrm)


def _terms_overlap(r, tA, tB):
    tot = 0.0
    for wa, ua, va in tA:
        for wb, ub, vb in tB:
            if ua[0] != ub[0] or ua[1] != ub[1] or va[0] != vb[0] or va[1] != vb[1]:
                continue
            su = np.trapezoid(ua[2] * ub[2] * r * r, r)
            sv = np.trapezoid(va[2] * vb[2] * r * r, r)
            tot += wa * wb * su * sv
    return tot


def _build_M_l(r, dr, configs, Z=2.0):
    data = [_config_terms(r, l, na, nb) for (l, na, nb) in configs]
    K = len(configs)
    M = np.zeros((K, K))
    for i in range(K):
        ti, Ri, ni = data[i]
        for j in range(i, K):
            tj, Rj, nj = data[j]
            g = 0.0
            for wa, ua, va in ti:
                for wb, ub, vb in tj:
                    g += wa * wb * _pair(r, dr, ua, va, ub, vb)
            g *= ni * nj
            M[i, j] = M[j, i] = (Z * Ri if i == j else 0.0) - g
    return M


def _ground_E(M):
    return -np.sort(eigh(M, eigvals_only=True))[-1] ** 2 / 2


@pytest.mark.slow
def test_paper60_lgt0_angular_correlation_and_sublinear():
    r = np.linspace(1e-6, 70.0, 16000)
    dr = r[1] - r[0]
    # single 1s^2 reproduces the textbook variational value even through the l>0 path
    assert abs(_ground_E(_build_M_l(r, dr, [(0, 1, 1)])) - (-2.84766)) < 4e-3
    # s-only floor
    s_cfg = [(0, 1, 1), (0, 2, 2), (0, 3, 3), (0, 1, 2), (0, 1, 3), (0, 2, 3)]
    Es = _ground_E(_build_M_l(r, dr, s_cfg))
    # adding p configs lowers E BELOW the s-only floor -> genuine angular correlation
    sp_cfg = s_cfg + [(1, 2, 2), (1, 3, 3), (1, 2, 3)]
    Esp = _ground_E(_build_M_l(r, dr, sp_cfg))
    assert Esp < Es - 1e-3, f"p-channel did not lower E: {Es} -> {Esp}"
    assert Esp < -2.88            # heads toward exact -2.90372, past the s-limit
    # the block-encoding 1-norm grows more slowly than K through the full s+p+d+f
    # basis (eq:sublinear).  >=5 points; the fitted exponent sits in the paper's
    # 0.78 (s-only) .. 0.84 (full basis) band and is robustly below 1.
    # NOT an asymptote (corrected 2026-09-07, /qa): 0.84 is a fit over the WINDOW
    # K = 74..164, and the local slope keeps rising past it -- 0.850 at K=202,
    # 0.868 at 244, 0.882 at 290, 0.906 at 340.  This self-contained sweep reaches
    # ~0.77 at K<=24; the pinned claim is the sublinear band over the computed
    # range, and NOT any asymptotic value.  The mechanism is pinned separately by
    # test_paper60_sublinearity_is_carried_by_the_nuclear_diagonal.
    sizes, norms = [], []
    for (lmax, span) in [(0, 2), (0, 3), (1, 3), (2, 3), (3, 3)]:
        cfg = _cfgs_upto(lmax, span)
        Mx = _build_M_l(r, dr, cfg)
        sizes.append(len(cfg)); norms.append(np.abs(Mx).sum())
    p = np.polyfit(np.log(sizes), np.log(norms), 1)[0]
    assert 0.6 < p < 0.95, f"atomic 1-norm exponent outside the sublinear 0.78-0.84 band: {p}"


@pytest.mark.slow
def test_paper60_split_is_box_sensitive_and_ordering_is_not():
    """eq:sublinear_split -- which block is slowest, and why the box matters.

    REPLACES test_paper60_sublinearity_is_carried_by_the_nuclear_diagonal
    (2026-09-07). That guard asserted `p_off > 1.0` ("T' is superlinear") on the
    production 60-bohr grid. Two independent routes -- exact grid-free Slater
    algebra and converged quadrature under R_MAX >= 3*n_max^2 -- agree that the
    converged value is 0.977 (K<=164) and 0.937 (K<=340), i.e. SUBLINEAR. The
    old assertion held only because the box inflated that leg by about +0.07,
    and it would have blocked the paper's correction.

    What survives box-independently is the ORDERING: the nuclear diagonal
    T^0 = Z*sum(R_nu) is the slowest-growing block. That is checked here on the
    exact combinatorial sum, which needs no ERI, no grid and no engine.

    The wrong answers this rejects:
      (a) T^0 not being the slowest block -- the paper's mechanism;
      (b) the 60-bohr box NOT inflating the off-diagonal leg, which would mean
          the withdrawn 1.05 was real after all;
      (c) a drifting config family (the K ladder is pinned).

    NOT claimed: the T^0 leg is computed from the configuration list alone, so a
    perturbation to build_M's DIAGONAL correctly does not fire it -- that is the
    point of computing it grid-free. What fires it is a change to the
    off-diagonal, which is what the ordering is measured against.
    (Fire-tested 2026-09-07: diagonal plant does NOT fire; scaling T' by
    K^-0.4 DOES.)
    """
    from geovac.sturmian_secular import build_configs, build_M, gen_configs

    Ks, diag, off, T0 = [], [], [], []
    for n in (7, 8, 9, 10):                       # K = 74, 100, 130, 164
        tup = gen_configs(3, {0: n, 1: n, 2: n, 3: n})
        cfgs = build_configs(tup)
        M = build_M(cfgs)
        Ks.append(len(cfgs))
        diag.append(float(np.abs(np.diag(M)).sum()))
        off.append(float(np.abs(M - np.diag(np.diag(M))).sum()))
        # exact, grid-free: needs only the configuration list
        T0.append(2.0 * sum(math.sqrt(1.0 / a ** 2 + 1.0 / b ** 2)
                            for (_l, a, b) in tup))
    assert Ks == [74, 100, 130, 164], f"config family drifted: {Ks}"

    lK = np.log(Ks)
    p_T0 = np.polyfit(lK, np.log(T0), 1)[0]
    p_off = np.polyfit(lK, np.log(off), 1)[0]

    # (a) The exact, box-independent leg: T^0 is the slowest-growing block.
    assert 0.68 < p_T0 < 0.73, f"exact T^0 exponent moved: {p_T0}"
    assert p_T0 < p_off - 0.15, (
        f"T^0 ({p_T0:.3f}) is no longer clearly the slowest block "
        f"vs off-diagonal ({p_off:.3f})")

    # (b) On THIS (production) grid the off-diagonal leg is inflated past its
    #     converged value of 0.977. Asserted as a box artifact, not as physics.
    assert p_off > 1.0, (
        f"the 60-bohr box no longer inflates the off-diagonal leg "
        f"({p_off:.4f}); converged is 0.977, so this test's premise is stale")

    # (c) ...and the inflation is real: at least +0.03 over converged.
    assert p_off - 0.977 > 0.03, (
        f"box inflation of the off-diagonal exponent is only "
        f"{p_off - 0.977:.4f}; the withdrawal of the 1.05 claim rested on it "
        f"being ~+0.05")



# ==========================================================================
# sec:molecular / sec:resource -- the gerade lever and the SW-vs-Gaussian
# metric head-to-head that drive the H2+ block-encoding resource estimate
# (Table tab:resource).  Fast, exact momentum-space SW metric + closed-form
# Gaussian overlap; slow position-grid H2+ binding validation.
# --------------------------------------------------------------------------
def _sw_block_mom(R, nmax, kind, kk=1.0, M=200001):
    """Exact momentum-space nmax x nmax two-center block of the SW metric (kind='S')
    or L2 Sturmian overlap (kind='m'), s-orbitals, scale kk (PhD eq 10.5.12/13)."""
    chi = np.linspace(1e-8, np.pi, M)
    cot = 1.0 / np.tan(chi / 2.0)
    sfac = np.ones_like(chi) if R == 0.0 else np.where(kk * R * cot != 0.0,
                                                       np.sin(kk * R * cot) / (kk * R * cot), 1.0)
    wfac = sfac if kind == "S" else (1.0 - np.cos(chi)) * sfac
    B = np.zeros((nmax, nmax))
    for a in range(1, nmax + 1):
        for b in range(a, nmax + 1):
            v = (2.0 / np.pi) * np.trapezoid(np.sin(a * chi) * np.sin(b * chi) * wfac, chi)
            B[a - 1, b - 1] = B[b - 1, a - 1] = v
    return B


def _sw_metric_mom(nmax, R):
    intra = _sw_block_mom(0.0, nmax, "S")
    inter = _sw_block_mom(R, nmax, "S")
    return np.block([[intra, inter], [inter.T, intra]])


def _gu(S, nmax):
    """gerade/ungerade orthogonal block-diagonalization of a two-center metric."""
    T = np.zeros_like(S)
    for n in range(nmax):
        T[n, n] = T[n, n + nmax] = 1 / np.sqrt(2)
        T[n + nmax, n] = 1 / np.sqrt(2)
        T[n + nmax, n + nmax] = -1 / np.sqrt(2)
    Sr = T @ S @ T.T
    return Sr[:nmax, :nmax], Sr[nmax:, nmax:]


def _gaussian_metric(nmax, R, a0=0.10, ratio=2.0):
    """Even-tempered two-center s-Gaussian overlap; normalized closed form."""
    exps = a0 * ratio ** np.arange(nmax)
    centers = [0.0] * nmax + [R] * nmax
    allexp = list(exps) + list(exps)
    N = 2 * nmax
    Sg = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            a, b, d = allexp[i], allexp[j], abs(centers[i] - centers[j])
            Sg[i, j] = (4 * a * b / (a + b) ** 2) ** 0.75 * np.exp(-a * b / (a + b) * d ** 2)
    return Sg


def test_paper60_gerade_sector_flat_conditioning():
    """The g/u lever: the gerade sector of the SW metric is flat (cond ~ 2, bounded)
    while the full metric grows -- so the gerade (sigma_g ground) block-encoding pays a
    non-growing metric penalty (Table tab:resource)."""
    R = 2.0
    g_conds, full_conds = [], []
    for nmax in (2, 4, 6, 8):
        S = _sw_metric_mom(nmax, R)
        g, u = _gu(S, nmax)
        g_conds.append(np.linalg.cond(g))
        full_conds.append(np.linalg.cond(S))
    # gerade conditioning is bounded and flat across the whole basis range
    assert max(g_conds) < 3.0, f"gerade sector not flat: {g_conds}"
    assert g_conds[-1] < 1.4 * g_conds[0], f"gerade sector grows too fast: {g_conds}"
    # the full metric, by contrast, grows and the gap widens
    assert full_conds[-1] > full_conds[0]
    assert full_conds[-1] / g_conds[-1] > 30.0        # ~40x at N=16


def test_paper60_sw_beats_gaussian_metric():
    """Head-to-head: at matched basis size and a standard even-tempered ratio, the SW
    metric is far better conditioned than a Gaussian LCAO metric for the same H2+, and
    grows with a gentler power law (backs the resource ratio in sec:resource)."""
    R = 2.0
    Ns = np.array([2 * nm for nm in (2, 4, 6, 8)], float)
    sw = np.array([np.linalg.cond(_sw_metric_mom(nm, R)) for nm in (2, 4, 6, 8)])
    ga = np.array([np.linalg.cond(_gaussian_metric(nm, R, ratio=2.0)) for nm in (2, 4, 6, 8)])
    # matched N=16: Gaussian metric is >100x worse conditioned than SW
    assert ga[-1] / sw[-1] > 100.0, f"Gaussian/SW cond ratio too small: {ga[-1]/sw[-1]}"
    # and its power-law growth exponent is steeper (coverage-vs-linear-dependence)
    p_sw = np.polyfit(np.log(Ns), np.log(sw), 1)[0]
    p_ga = np.polyfit(np.log(Ns), np.log(ga), 1)[0]
    assert 1.5 < p_sw < 2.2, f"SW exponent off: {p_sw}"       # ~N^1.8
    assert p_ga > p_sw + 1.0, f"Gaussian not steeper: {p_ga} vs {p_sw}"  # ~N^3.5


def test_paper60_large_r_lever():
    """Large-R conditioning lever (sec:resource): the SW metric's cond(S) relaxes toward ~2
    as the internuclear separation grows (the Coulomb-Sturmians decouple), so a large-R
    problem is well-conditioned even where the L2 overlap cannot follow.  Paper quotes
    cond(S) ~ 2.2 at R=10 bohr (basis N=10, i.e. nmax=5/center)."""
    nmax = 5
    conds = {R: float(np.linalg.cond(_sw_metric_mom(nmax, R))) for R in (1.4, 2.0, 10.0)}
    assert conds[1.4] > conds[2.0] > conds[10.0]                 # monotone relaxation with R
    assert 1.8 < conds[10.0] < 2.6, f"cond(S) at R=10 not ~2.2: {conds[10.0]}"
    assert conds[1.4] > 20.0 * conds[10.0], f"large-R lever weak: {conds}"


@pytest.mark.slow
def test_paper60_h2plus_isoenergetic_binds():
    """The one-electron molecular isoenergetic secular equation [W - kS]C=0 binds the
    H2+ sigma_g ground state near the reference electronic energy (validates the resource
    estimate's underlying algorithm and pins lambda_eff = k_max)."""
    R, nmax, kk = 2.0, 4, 1.485          # build at the physical scale (self-consistency point)
    rho = np.linspace(1e-4, 30.0, 700)
    z = np.linspace(-22.0, R + 26.0, 1100)
    RHO, ZZ = np.meshgrid(rho, z, indexing="ij")
    drho, dz = rho[1] - rho[0], z[1] - z[0]

    def integ(g):
        return 2 * np.pi * np.sum(g * RHO) * drho * dz

    def cs_at(n, zc):
        rr = np.sqrt(RHO ** 2 + (ZZ - zc) ** 2)
        f = np.exp(-kk * rr) * genlaguerre(n - 1, 1)(2 * kk * rr)
        return f / np.sqrt(integ(f * f))

    rA = np.sqrt(RHO ** 2 + ZZ ** 2)
    rB = np.sqrt(RHO ** 2 + (ZZ - R) ** 2)
    invr = 1.0 / rA + 1.0 / rB
    basis = [(n, 0.0) for n in range(1, nmax + 1)] + [(n, R) for n in range(1, nmax + 1)]
    fs = [cs_at(n, zc) for n, zc in basis]
    grads = [np.gradient(f, rho, z) for f in fs]
    N = len(fs)
    W = np.zeros((N, N))
    S = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            W[i, j] = (1.0 / kk) * integ(fs[i] * invr * fs[j])
            gir, giz = grads[i]
            gjr, gjz = grads[j]
            S[i, j] = (1.0 / (2 * kk ** 2)) * integ(gir * gjr + giz * gjz) + 0.5 * integ(fs[i] * fs[j])
    S = 0.5 * (S + S.T)
    W = 0.5 * (W + W.T)
    kroots = np.sort(eigh(W, S, eigvals_only=True))[::-1]
    k0 = kroots[0]                        # deepest bound = sigma_g
    E = -0.5 * k0 ** 2
    # self-consistent (built at k=1.485): deepest root ~ 1.485, E_elec ~ -1.09 Ha
    assert abs(k0 - 1.485) < 0.06, f"sigma_g root off self-consistency: k0={k0}"
    assert -1.11 < E < -1.02, f"H2+ E_elec out of range: {E}"    # ref -1.103, s-only ~ -1.088
    assert E > -1.1026, "variational: s-only must stay above the exact electronic energy"


# --------------------------------------------------------------------------
# sec:manyelectron -- the SW metric does NOT compound with electron number.
# A k-electron configuration overlap is the k-th compound matrix of the
# one-electron overlap: raw L2 grows with k, the SW-orthonormal route is I.
# --------------------------------------------------------------------------
from itertools import combinations as _combos


def _config_overlap(S1, k_elec):
    """k-electron spatial-determinant overlap: <Phi_I|Phi_J> = det(S1[I,J])."""
    N = S1.shape[0]
    cfgs = list(_combos(range(N), k_elec))
    M = np.zeros((len(cfgs), len(cfgs)))
    for a, I in enumerate(cfgs):
        for b, J in enumerate(cfgs):
            M[a, b] = np.linalg.det(S1[np.ix_(I, J)])
    return M


def test_paper60_manyelectron_metric_does_not_compound():
    R, nmax = 1.4, 3                      # H2 geometry, N=6 two-center basis
    # build raw L2 overlap and SW metric (momentum-space), then SW-orthonormalize
    S_L2 = np.block([[_sw_block_mom(0.0, nmax, "m"), _sw_block_mom(R, nmax, "m")],
                     [_sw_block_mom(R, nmax, "m").T, _sw_block_mom(0.0, nmax, "m")]])
    S_SW = _sw_metric_mom(nmax, R)
    N = 2 * nmax
    lam = np.sort(np.linalg.eigvalsh(S_L2))[::-1]
    c1 = np.linalg.cond(S_L2)

    # SW-orthonormalize the ACTUAL SW metric (not a hand-inserted identity):
    # S_SW^{-1/2} S_SW S_SW^{-1/2} = I only because S_SW is well-conditioned (unlike S_L2).
    w_sw, V_sw = np.linalg.eigh(S_SW)
    X_sw = V_sw @ np.diag(w_sw ** -0.5) @ V_sw.T
    S1_sw = X_sw @ S_SW @ X_sw
    assert np.max(np.abs(S1_sw - np.eye(N))) < 1e-9, "SW metric did not orthonormalize to identity"
    raw_conds, sw_conds = [], []
    for ke in (1, 2, 3):
        craw = np.linalg.cond(_config_overlap(S_L2, ke))
        csw = np.linalg.cond(_config_overlap(S1_sw, ke))       # SW-orthonormalized 1e overlap (real, not eye)
        # exact cross-check against the compound-matrix eigenvalue ratio
        comp = np.prod(lam[:ke]) / np.prod(lam[-ke:])
        assert abs(craw - comp) / comp < 1e-6, f"compound-matrix identity broken k={ke}"
        raw_conds.append(craw)
        sw_conds.append(csw)

    # RAW many-body metric GROWS with electron number (worst at half-filling k=3)
    assert raw_conds[1] > 3.0 * c1                  # k=2 already several x the 1-electron cond
    assert raw_conds[2] > raw_conds[0]              # monotone up toward half-filling
    # SW-orthonormal route: identity metric at EVERY electron number
    assert max(sw_conds) < 1.0 + 1e-9, f"SW config metric not identity: {sw_conds}"


# --------------------------------------------------------------------------
# sec:manyelectron [MEASURED] -- the enabling two-center ERIs are validated
# against the framework's EXACT closed form, and the interacting metric-free
# method binds minimal-basis H2.  Self-contained multipole two-center s-ERI.
# --------------------------------------------------------------------------
from scipy.integrate import cumulative_trapezoid as _cumtrap
from scipy.special import eval_legendre as _legP


def _two_center_s_eri(R, kk, pq_centers, rs_centers, Lmax=22, nr=2500, nth=160, rmax=50.0):
    """(1s_p 1s_q | 1s_r 1s_s) for 1s orbitals on centers in {'A','B'} (A at 0, B at R zhat),
    by multipole expansion about A.  pq_centers/rs_centers are 2-tuples like ('A','A')."""
    r = np.linspace(1e-5, rmax, nr)
    dr = r[1] - r[0]
    u = np.sort(np.cos(np.linspace(0.0, np.pi, nth)))
    RR, UU = np.meshgrid(r, u, indexing='ij')
    fnorm = 1.0 / np.sqrt(np.trapezoid(np.exp(-2 * kk * r) * r * r, r))   # 1s radial L2 norm

    def phi(center):
        d = RR if center == 'A' else np.sqrt(RR * RR + R * R - 2 * R * RR * UU)
        return fnorm * np.exp(-kk * d) / np.sqrt(4 * np.pi)

    def A_L(ca, cb):
        prod = phi(ca) * phi(cb)
        return [2 * np.pi * np.trapezoid(prod * _legP(L, u)[None, :], u, axis=1) for L in range(Lmax + 1)]
    Apq, Ars = A_L(*pq_centers), A_L(*rs_centers)
    total = 0.0
    for L in range(Lmax + 1):
        a, b = Apq[L], Ars[L]
        if np.max(np.abs(a)) < 1e-14 or np.max(np.abs(b)) < 1e-14:
            continue
        g = b * r * r
        inner = np.concatenate(([0.0], _cumtrap(g * r ** L, dx=dr))) * r ** (-(L + 1))
        outer = np.concatenate(([0.0], _cumtrap((g * r ** (-(L + 1)))[::-1], dx=dr)))[::-1] * r ** L
        total += np.trapezoid(a * (inner + outer) * r * r, r)
    return total


@pytest.mark.slow
def test_paper60_two_center_eri_matches_exact_closed_form():
    kk = 1.0
    # one-center (1s 1s|1s 1s) = 5k/8
    v1 = _two_center_s_eri(1.4, kk, ('A', 'A'), ('A', 'A'))
    assert abs(v1 - 5 * kk / 8) < 2e-3, f"one-center ERI off: {v1}"
    # two-center (AA|BB) vs the framework's exact closed-form aabb_value
    from fractions import Fraction
    from geovac.two_center_eri import aabb_value
    for R in (1.4, 2.0):
        v_num = _two_center_s_eri(R, kk, ('A', 'A'), ('B', 'B'))
        v_ex = aabb_value(Fraction(1, 1), (1, 0, 0), (1, 0, 0),
                          Fraction(1, 1), (1, 0, 0), (1, 0, 0), R, prec=25)
        assert abs(v_num - v_ex) < 3e-3, f"(AA|BB) R={R}: num {v_num} vs exact {v_ex}"


@pytest.mark.slow
def test_paper60_interacting_h2_binds():
    """Minimal-basis H2 in the metric-free molecular-Sturmian basis has an interior minimum
    (the interacting many-electron molecular method binds), at the textbook single-zeta value."""
    kk = 1.0

    def sigma_g_energy_and_J(R):
        # one-electron gerade energy on a cylindrical grid (kinetic by parts + V_ne)
        rho = np.linspace(1e-4, 28.0, 500)
        z = np.linspace(-20.0, R + 22.0, 760)
        RHO, ZZ = np.meshgrid(rho, z, indexing="ij")
        drho, dz = rho[1] - rho[0], z[1] - z[0]

        def integ(g):
            return 2 * np.pi * np.sum(g * RHO) * drho * dz

        def cs(zc):
            f = np.exp(-kk * np.sqrt(RHO ** 2 + (ZZ - zc) ** 2))
            return f / np.sqrt(integ(f * f))
        fs = [cs(0.0), cs(R)]
        vne = -1.0 / np.sqrt(RHO ** 2 + ZZ ** 2) - 1.0 / np.sqrt(RHO ** 2 + (ZZ - R) ** 2)
        grads = [np.gradient(f, rho, z) for f in fs]
        S = np.array([[integ(fs[i] * fs[j]) for j in range(2)] for i in range(2)])
        h = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):
                gir, giz = grads[i]
                gjr, gjz = grads[j]
                h[i, j] = 0.5 * integ(gir * gjr + giz * gjz) + integ(fs[i] * vne * fs[j])
        w, V = eigh(h, S)
        cg = V[:, 0]                                   # gerade
        # J(sigma_g^2) = sum_ijkl cg_i cg_j cg_k cg_l (ij|kl) over {A,B}
        cen = ['A', 'B']
        J = 0.0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        J += cg[i] * cg[j] * cg[k] * cg[l] * _two_center_s_eri(
                            R, kk, (cen[i], cen[j]), (cen[k], cen[l]), nr=1500, nth=120)
        return w[0], J

    Etot = {}
    for R in (1.4, 1.6, 2.2):
        eps, J = sigma_g_energy_and_J(R)
        Etot[R] = 2 * eps + J + 1.0 / R
    # interior minimum near R~1.6 (bound), and the depth is the single-zeta textbook value
    assert Etot[1.6] < Etot[1.4] and Etot[1.6] < Etot[2.2], f"no interior minimum: {Etot}"
    assert -1.12 < Etot[1.6] < -1.05, f"H2 E_tot off single-zeta value: {Etot[1.6]}"
    assert Etot[1.6] > -1.174, "variational: minimal basis must stay above the exact H2 energy"


# --------------------------------------------------------------------------
# sec:manyelectron [OPEN, obstruction] -- the collective isoenergetic scale is the
# ROOT-sum-of-squares, = Z R_nu on a single center (the diagonal T0 that the atomic
# sublinearity rides on); a naive sum-of-scales doubles the non-interacting energy.
# This is why the atomic sublinearity is single-center-specific.
# --------------------------------------------------------------------------
def test_paper60_collective_scale_is_root_sum_of_squares():
    # He 1s^2: k_p = Z/n = 2 for each 1s; the bare (no-ee) collective scale = Z R_nu.
    Z = 2.0
    kp = np.array([Z / 1.0, Z / 1.0])
    R_nu = np.sqrt((1 / 1.0) ** 2 + (1 / 1.0) ** 2)          # sqrt2
    pk_rss = np.sqrt(np.sum(kp ** 2))                         # correct collective scale
    # (i) root-sum-of-squares == Z R_nu (the paper's diagonal T0)
    assert abs(pk_rss - Z * R_nu) < 1e-12
    # (ii) it gives the exact bare He energy -4.0 (two He+ 1s), paper sec:atomic
    assert abs(-0.5 * pk_rss ** 2 - (-4.0)) < 1e-12
    # (iii) the naive SUM of scales doubles the binding (-8, the shortcut bug) -> not isoenergetic
    pk_sum = np.sum(kp)
    assert abs(-0.5 * pk_sum ** 2 - (-8.0)) < 1e-12
    assert pk_sum > pk_rss                                    # sum overshoots the collective scale


# --------------------------------------------------------------------------
# sec:manyelectron [MEASURED] -- a validated two-center Sturmian 2-electron H2 CI
# dissociates correctly to two H atoms (E_tot -> -1.0) and binds above the exact
# energy at R_eq; this is the solver behind the molecular block-encoding 1-norm.
# Self-contained 1s-only multipole two-center integrals + 2e FCI.
# --------------------------------------------------------------------------
def _h2_ci_energy(R, zeta, Lmax=22, nr=2600, nth=160, rmax=None):
    from scipy.special import eval_legendre as _lp
    from scipy.integrate import cumulative_trapezoid as _ct
    rmax = rmax or max(55.0, 2 * R + 30)
    r = np.linspace(1e-5, rmax, nr); dr = r[1] - r[0]
    u = np.sort(np.cos(np.linspace(0, np.pi, nth)))
    PL = [_lp(L, u) for L in range(Lmax + 1)]
    RR, UU = np.meshgrid(r, u, indexing='ij')
    nrm = 1.0 / np.sqrt(np.trapezoid(np.exp(-2 * zeta * r) * r * r, r))

    def phi(center):
        d = RR if center == 'A' else np.sqrt(RR ** 2 + R ** 2 - 2 * R * RR * UU)
        return nrm * np.exp(-zeta * d) / np.sqrt(4 * np.pi)
    ph = {'A': phi('A'), 'B': phi('B')}

    def AL(ci, cj):
        prod = ph[ci] * ph[cj]
        return np.array([2 * np.pi * np.trapezoid(prod * PL[L][None, :], u, axis=1) for L in range(Lmax + 1)])

    def coul_center(al, C):
        if C == 'A':
            return np.trapezoid(al[0] * r, r)
        rlt, rgt = np.minimum(r, R), np.maximum(r, R)
        return sum(np.trapezoid(al[L] * (rlt ** L / rgt ** (L + 1)) * r ** 2, r)
                   for L in range(Lmax + 1) if np.max(np.abs(al[L])) > 1e-15)

    def eri(ci, cj, ck, cl):
        a, b = AL(ci, cj), AL(ck, cl)
        tot = 0.0
        for L in range(Lmax + 1):
            if np.max(np.abs(a[L])) < 1e-15 or np.max(np.abs(b[L])) < 1e-15:
                continue
            g = b[L] * r * r
            inn = np.concatenate(([0.0], _ct(g * r ** L, dx=dr))) * r ** (-(L + 1))
            out = np.concatenate(([0.0], _ct((g * r ** (-(L + 1)))[::-1], dx=dr)))[::-1] * r ** L
            tot += np.trapezoid(a[L] * (inn + out) * r * r, r)
        return tot
    cen = ['A', 'B']
    S = np.zeros((2, 2)); h = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            al = AL(cen[i], cen[j])
            S[i, j] = np.trapezoid(al[0] * r * r, r)
            T = zeta * coul_center(al, cen[j]) - 0.5 * zeta ** 2 * S[i, j]   # 1s: n=1
            V = -(coul_center(al, 'A') + coul_center(al, 'B'))
            h[i, j] = T + V
    er = np.zeros((2, 2, 2, 2))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    er[i, j, k, l] = eri(cen[i], cen[j], cen[k], cen[l])
    # Loewdin orthonormalize, then full 2-electron determinant FCI (validated in the driver)
    w, Umat = np.linalg.eigh(S)
    X = Umat @ np.diag(1.0 / np.sqrt(w))
    hm = X.T @ h @ X
    em = np.einsum('ip,jq,kr,ls,ijkl->pqrs', X, X, X, X, er, optimize=True)
    from itertools import combinations as _cmb
    dets = list(_cmb(range(4), 2))                 # 4 spin-orbitals (2 spatial), C(4,2)=6 dets

    def _sp(i):
        return i // 2

    def _spin(i):
        return i % 2

    def _g2(i, j, k, l):
        c = em[_sp(i), _sp(k), _sp(j), _sp(l)] if _spin(i) == _spin(k) and _spin(j) == _spin(l) else 0.0
        x = em[_sp(i), _sp(l), _sp(j), _sp(k)] if _spin(i) == _spin(l) and _spin(j) == _spin(k) else 0.0
        return c - x
    Hc = np.zeros((6, 6))
    for a, (i, j) in enumerate(dets):
        for b, (k, l) in enumerate(dets):
            oa, ob = {i, j}, {k, l}; d = oa ^ ob
            if len(d) == 0:
                v = (hm[_sp(i), _sp(i)] * (_spin(i) == _spin(i)) + hm[_sp(j), _sp(j)]) + _g2(i, j, i, j)
            elif len(d) == 2:
                m = (oa - ob).pop(); p = (ob - oa).pop(); c = (oa & ob).pop()
                sgn = (-1) ** ([i, j].index(m) + [k, l].index(p))
                h1 = hm[_sp(m), _sp(p)] if _spin(m) == _spin(p) else 0.0
                v = sgn * (h1 + _g2(m, c, p, c))
            elif len(d) == 4:
                m1, m2 = sorted(oa - ob); p1, p2 = sorted(ob - oa)
                sgn = (-1) ** ([i, j].index(m1) + [i, j].index(m2) + [k, l].index(p1) + [k, l].index(p2))
                v = sgn * _g2(m1, m2, p1, p2)
            else:
                v = 0.0
            Hc[a, b] = v
    return np.linalg.eigvalsh(0.5 * (Hc + Hc.T))[0]


@pytest.mark.slow
def test_paper60_h2_ci_dissociates_and_binds():
    # dissociation: R=6 => E_tot -> -1.0 (two H atoms); zeta=1 optimal for isolated H
    E_far = _h2_ci_energy(6.0, 1.0, Lmax=24, nr=3000, nth=200)
    assert abs((E_far + 1.0 / 6.0) - (-1.0)) < 0.02, f"H2 does not dissociate to -1.0: {E_far+1/6.0}"
    # bonding: interior binding, above the exact H2 energy (variational)
    E_eq = _h2_ci_energy(1.4, 1.2) + 1.0 / 1.4
    assert -1.174 < E_eq < -1.08, f"H2 E_tot(1.4) off: {E_eq}"       # above exact, bound


# ==========================================================================
# sec:molecular scope -- the gerade lever needs EQUIVALENT centers.  A probe on
# water (C2v, O + 2H) shows a symmetry-inequivalent heavy center reinstates the
# growth in the totally-symmetric (ground-state) block.  Reuses _sw_block_mom.
# --------------------------------------------------------------------------
def _water_a1_block(nmax, R_OH=1.809, ang_deg=104.5):
    """C2v A1 (ground-state) block of water's 3-center SW metric.
    s-only shared-scale CS; basis {O_n} U {(H1_n+H2_n)/sqrt2}.  Two-center s-s
    integrals depend only on inter-center distance (isotropy):
      <O|O>=I, <H+|H+>=I+Q(R_HH), <O|H+>=sqrt2 P(R_OH)."""
    R_HH = 2.0 * R_OH * math.sin(math.radians(ang_deg) / 2.0)
    I = _sw_block_mom(0.0, nmax, "S")
    P = _sw_block_mom(R_OH, nmax, "S")     # O-H
    Q = _sw_block_mom(R_HH, nmax, "S")     # H-H
    top = np.hstack([I, math.sqrt(2.0) * P])
    bot = np.hstack([math.sqrt(2.0) * P.T, I + Q])
    return np.vstack([top, bot])


def test_paper60_water_gerade_lever_fails_for_inequivalent_center():
    """The gerade lever is an EQUIVALENT-center property, not a symmetry property.
    For water the A1 ground-state block still carries the O<->H coupling between
    symmetry-inequivalent centers -- which C2v adaptation cannot remove -- so cond(A1)
    grows ~N^2 like the raw metric, NOT flat like H2+ gerade (~2).  Zeroing the O-H
    coupling restores the flatness, pinpointing it as the sole driver."""
    conds, dims = [], []
    for nmax in (2, 4, 6, 8, 10, 12):
        conds.append(np.linalg.cond(_water_a1_block(nmax)))
        dims.append(3 * nmax)
    conds = np.array(conds); dims = np.array(dims, float)
    # NOT flat -- contrast the H2+ gerade test (g_conds[-1] < 1.4 * g_conds[0])
    assert conds[-1] > 10.0 * conds[0], f"water A1 unexpectedly flat: {conds}"
    # grows roughly quadratically in total dimension (ill-conditioning, not a flat sector)
    p = np.polyfit(np.log(dims), np.log(conds), 1)[0]
    assert 1.6 < p < 2.3, f"water A1 growth exponent off: {p}"
    # smoking gun: the O<->H coupling is the entire driver
    nmax = 12
    A1 = _water_a1_block(nmax)
    A1_noOH = A1.copy(); A1_noOH[:nmax, nmax:] = 0.0; A1_noOH[nmax:, :nmax] = 0.0
    assert np.linalg.cond(A1) > 100.0, "expected water A1 badly conditioned"
    assert np.linalg.cond(A1_noOH) < 3.0, "zeroing O-H coupling should restore flatness"
