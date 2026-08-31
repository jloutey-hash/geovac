"""Backing for Paper 34's e-e partial-split remark (rem:ee_partial_split).

The exact split of the multipole Coulomb kernels via min/max:

    K_L = r_<^L / r_>^{L+1}
        = 1/2 [ phi(x)psi(y) + phi(y)psi(x) ]  -  1/2 |phi(x)psi(y) - phi(y)psi(x)|
          (separable head, rank 2)                (W_L >= 0, carries the kink)

with phi = r^L, psi = r^{-(L+1)}.  At L = 0 the head is pure labels x metric
(A = S, B = (k/n) diag); at L >= 1 the head is still rank-2 one-body x one-body
but B = <r^-(L+1)> is dense.

Legs:
  1. split identity + total scale-out g(k) = k g(1)          (L = 0)
  2. head rank 2 + labels form V = (k/n) delta               (L = 0)
  3. FCI with the separable part ONLY is exactly uncorrelated:
     equals the dressed one-body problem h1 + (N-1)/2 V bit-exactly
  4. payoff + teeth: rank-4 W reaches |dE| < 0.5 mHa of full FCI while
     rank-0 (separable only) is off by > 100 mHa
  5. general L: split exact, head = (A x B + B x A)/2, rank 2, for L = 1, 2

Provenance: debug/ee_split_probe.py, debug/ee_split_L_probe.py,
debug/ee_split_rank_sweep.py (2026-08-27) and
debug/sprint_minimal_presentation_memo.md section 8.
"""
from math import factorial

import numpy as np
import pytest
from scipy.linalg import eigh
from scipy.special import eval_genlaguerre

from geovac import transcorrelated_sturmian as TC


def _build(ns, k, Z=2.0, Ng=500):
    r, wr = TC.make_grid(k, Ng=Ng)
    S, h1, Rtab, W = TC.build_one_body(ns, r, wr, k, Z)
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    K_full = 1.0 / np.maximum(R1, R2)
    K_sep = 0.5 * (1.0 / R1 + 1.0 / R2)
    K_W = 0.5 * np.abs(1.0 / R1 - 1.0 / R2)
    D = {(i, kk): Rtab[i + 1][0] * Rtab[kk + 1][0] * W
         for i in range(ns) for kk in range(ns)}

    def eri(K):
        g = np.zeros((ns,) * 4)
        for i in range(ns):
            for j in range(ns):
                for kk in range(ns):
                    for ll in range(ns):
                        g[i, j, kk, ll] = D[(i, kk)] @ K @ D[(j, ll)]
        return g

    Vr = np.array([[np.sum(Rtab[i + 1][0] * Rtab[j + 1][0] * r * wr)
                    for j in range(ns)] for i in range(ns)])
    return S, h1, Vr, eri(K_full), eri(K_sep), eri(K_W)


def _pairmat(g, ns):
    return g.transpose(0, 2, 1, 3).reshape(ns * ns, ns * ns)


def _fci(ns, h1o, g_tensor, X, N_elec):
    nso = 2 * ns
    dets, didx = TC.make_dets(nso, N_elec)
    go = TC.transform_2(g_tensor, X)
    H = TC.build_H(dets, didx, TC.h_spin(h1o, nso),
                   TC.asym_from_phys(go, nso), nso)
    return float(eigh(H, eigvals_only=True)[0])


def test_split_identity_and_total_scale_out():
    ns = 4
    _, _, _, g2, gs2, gw2 = _build(ns, 2.0)
    _, _, _, g1, _, _ = _build(ns, 1.0)
    assert np.abs(g2 - (gs2 - gw2)).max() / np.abs(g2).max() < 1e-12
    assert np.abs(g2 / 2.0 - g1).max() / np.abs(g1).max() < 1e-12


def test_head_rank2_and_labels_form():
    ns, k = 4, 2.0
    S, _, Vr, _, gs, _ = _build(ns, k)
    n = np.arange(1, ns + 1)
    V = k * np.diag(1.0 / n)
    assert np.abs(Vr - V).max() < 1e-6                    # <1/r> = (k/n) delta
    head = 0.5 * (np.einsum("ik,jl->ijkl", V, S) + np.einsum("ik,jl->ijkl", S, V))
    assert np.abs(gs - head).max() / np.abs(gs).max() < 1e-6
    sv = np.linalg.svd(_pairmat(gs, ns), compute_uv=False)
    assert sv[2] / sv[0] < 1e-10                          # rank 2 exactly


@pytest.mark.parametrize("N_elec", [2, 3])
def test_separable_part_supports_no_correlation(N_elec):
    ns, k = 4, 2.0
    S, h1, Vr, _, gs, _ = _build(ns, k)
    X = TC.lowdin(S)
    h1o, Vo = TC.transform_1(h1, X), TC.transform_1(Vr, X)
    E_sep = _fci(ns, h1o, gs, X, N_elec)
    nso = 2 * ns
    dets, didx = TC.make_dets(nso, N_elec)
    h_eff = TC.h_spin(h1o + 0.5 * (N_elec - 1) * Vo, nso)
    H1b = TC.build_H(dets, didx, h_eff, np.zeros((nso,) * 4), nso)
    E_1b = float(eigh(H1b, eigvals_only=True)[0])
    assert abs(E_sep - E_1b) < 1e-10


def test_w_rank4_payoff_and_rank0_teeth():
    ns, k = 4, 2.0
    S, h1, _, g, gs, gw = _build(ns, k)
    X = TC.lowdin(S)
    h1o = TC.transform_1(h1, X)
    E_full = _fci(ns, h1o, g, X, 2)
    E_sep = _fci(ns, h1o, gs, X, 2)
    assert abs(E_sep - E_full) > 0.1                      # teeth: >100 mHa without W
    Wm = _pairmat(gw, ns)
    lam, U = np.linalg.eigh(Wm)
    o = np.argsort(-np.abs(lam))
    m = 4
    Wr = (U[:, o[:m]] * lam[o[:m]]) @ U[:, o[:m]].T
    gr = gs - Wr.reshape(ns, ns, ns, ns).transpose(0, 2, 1, 3)
    gr = 0.25 * (gr + gr.transpose(2, 1, 0, 3) + gr.transpose(0, 3, 2, 1)
                 + gr.transpose(2, 3, 0, 1))
    E_r = _fci(ns, h1o, gr, X, 2)
    assert abs(E_r - E_full) < 5e-4                       # rank 4 -> sub-0.5-mHa


def _R_nl(n, l, r, k):
    x = 2 * k * r
    norm = np.sqrt((2 * k) ** 3 * factorial(n - l - 1)
                   / (2 * n * factorial(n + l)))
    return norm * x ** l * np.exp(-x / 2) * eval_genlaguerre(n - l - 1, 2 * l + 1, x)


@pytest.mark.parametrize("L", [1, 2])
def test_general_L_split_and_rank2_head(L):
    k, Ng, nfun = 2.0, 500, 4
    r, wr = TC.make_grid(k, Ng=Ng)
    Wt = r * r * wr
    ns_list = list(range(L + 1, L + 1 + nfun))
    R = {n: _R_nl(n, L, r, k) for n in ns_list}
    A = np.array([[np.sum(R[m] * R[n] * r ** L * Wt) for n in ns_list]
                  for m in ns_list])
    B = np.array([[np.sum(R[m] * R[n] * r ** (-(L + 1)) * Wt) for n in ns_list]
                  for m in ns_list])
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    K_full = np.minimum(R1, R2) ** L / np.maximum(R1, R2) ** (L + 1)
    t1 = R1 ** L / R2 ** (L + 1)
    t2 = R2 ** L / R1 ** (L + 1)
    K_sep, K_W = 0.5 * (t1 + t2), 0.5 * np.abs(t1 - t2)
    # NB: the identity is checked at the INTEGRAL level below, not pointwise on the
    # kernels -- near the grid's r = 1e-9 floor the two split terms reach ~1e19 and
    # float cancellation cannot resolve their difference, while the r^2 measure makes
    # that region irrelevant to every integral.
    D = {(i, j): R[ns_list[i]] * R[ns_list[j]] * Wt
         for i in range(nfun) for j in range(nfun)}

    def eri(K):
        g = np.zeros((nfun,) * 4)
        for i in range(nfun):
            for j in range(nfun):
                for kk in range(nfun):
                    for ll in range(nfun):
                        g[i, j, kk, ll] = D[(i, kk)] @ K @ D[(j, ll)]
        return g

    g_full, gs, gw = eri(K_full), eri(K_sep), eri(K_W)
    assert np.abs(g_full - (gs - gw)).max() / np.abs(g_full).max() < 1e-10
    head = 0.5 * (np.einsum("ik,jl->ijkl", A, B) + np.einsum("ik,jl->ijkl", B, A))
    assert np.abs(gs - head).max() / np.abs(gs).max() < 1e-10
    sv = np.linalg.svd(_pairmat(gs, nfun), compute_uv=False)
    assert sv[2] / sv[0] < 1e-10


# ---------------------------------------------------------------------------
# the 2n-1 exact-rank law and the closed-form question (rem:ee_partial_split,
# second paragraph; provenance debug/ee_w_eigen_skeleton.py + _followup.py)
# ---------------------------------------------------------------------------
def test_pair_space_2n_minus_1_law():
    """s-wave pair densities span exactly 2n-1 dims => EVERY radial-kernel
    matricization has exact rank <= 2n-1 (g, W, and a random control), while
    the separable head stays rank 2."""
    rng = np.random.default_rng(1)
    for ns in (4, 6):
        r, wr = TC.make_grid(1.0, Ng=600)
        S, h1, Rtab, Wt = TC.build_one_body(ns, r, wr, 1.0, 2.0)
        prs = [(i, j) for i in range(1, ns + 1) for j in range(1, ns + 1)]
        Dm = np.array([Rtab[i][0] * Rtab[j][0] * Wt for (i, j) in prs])
        assert np.linalg.matrix_rank(Dm, tol=1e-10 * np.abs(Dm).max()) == 2 * ns - 1
        R1g, R2g = np.meshgrid(r, r, indexing="ij")
        c = rng.normal(size=4)
        # the law is a CEILING (rank <= 2n-1 for any kernel); kinked kernels
        # saturate it generically, very smooth kernels may sit below it.
        kernels = {
            "g": (1.0 / np.maximum(R1g, R2g), "eq"),
            "gsep": (0.5 * (1.0 / R1g + 1.0 / R2g), "two"),
            "W": (0.5 * np.abs(1.0 / R1g - 1.0 / R2g), "eq"),
            "rand_kinked": (c[0] + c[1] * np.exp(-0.3 * (R1g + R2g))
                            + c[2] / (1 + R1g + R2g)
                            + c[3] * np.exp(-0.1 * np.abs(R1g - R2g)), "eq"),
            "rand_smooth": (c[0] + c[1] * np.exp(-0.3 * (R1g + R2g))
                            + c[2] / (1 + R1g + R2g), "le"),
        }
        for nm, (K, mode) in kernels.items():
            M = np.array([[float(Dm[a] @ K @ Dm[b]) for b in range(len(prs))]
                          for a in range(len(prs))])
            got = int(np.linalg.matrix_rank(M, tol=1e-10 * np.abs(M).max()))
            if mode == "eq":
                assert got == 2 * ns - 1, f"ns={ns} {nm}: rank {got} != {2*ns-1}"
            elif mode == "two":
                assert got == 2, f"ns={ns} {nm}: rank {got} != 2"
            else:
                assert got <= 2 * ns - 1, f"ns={ns} {nm}: rank {got} > ceiling"


def test_exact_density_dependence_n3():
    """-rho_11 + 2 rho_12 - 3 rho_13 + 2 rho_22 == 0 identically (n = 3)."""
    r, wr = TC.make_grid(1.0, Ng=600)
    S, h1, Rtab, Wt = TC.build_one_body(3, r, wr, 1.0, 2.0)
    rho = {(i, j): Rtab[i][0] * Rtab[j][0] for i in (1, 2, 3) for j in (1, 2, 3)}
    combo = -rho[(1, 1)] + 2 * rho[(1, 2)] - 3 * rho[(1, 3)] + 2 * rho[(2, 2)]
    scale = max(np.abs(rho[(1, 1)]).max(), np.abs(rho[(2, 2)]).max())
    assert np.abs(combo).max() / scale < 1e-12


def test_w_rational_anchor_and_irreducible_charpoly():
    """Exact ns=2: W[1s^2,1s^2] = 3/8, and the active characteristic polynomial
    131072 x^3 - 50688 x^2 - 25200 x - 675 is IRREDUCIBLE over Q -- the finite
    eigenpairs have no low-degree closed form (the negative, pinned)."""
    import sympy as sp
    r1, r2 = sp.symbols("r1 r2", positive=True)

    def dens(m, n, var):
        p = (sp.Rational(2, m) * sp.assoc_laguerre(m - 1, 1, 2 * var)
             * sp.Rational(2, n) * sp.assoc_laguerre(n - 1, 1, 2 * var))
        return sp.expand(p * var ** 2)

    def w_entry(m1, n1, m2, n2):
        p1 = dens(m1, n1, r1) * sp.exp(-2 * r1)
        p2 = dens(m2, n2, r2) * sp.exp(-2 * r2)
        q1 = dens(m2, n2, r1) * sp.exp(-2 * r1)
        q2 = dens(m1, n1, r2) * sp.exp(-2 * r2)
        i1 = sp.integrate(p1 * (1 / r1 - 1 / r2) / 2, (r1, 0, r2))
        v1 = sp.integrate(sp.expand(p2 * i1), (r2, 0, sp.oo))
        i2 = sp.integrate(q1 * (1 / r1 - 1 / r2) / 2, (r1, 0, r2))
        v2 = sp.integrate(sp.expand(q2 * i2), (r2, 0, sp.oo))
        return sp.nsimplify(sp.simplify(v1 + v2), rational=True)

    pairs = [(1, 1), (1, 2), (2, 1), (2, 2)]
    Wx = sp.Matrix(4, 4, lambda a, b: 0)
    cache = {}
    for a, (m1, n1) in enumerate(pairs):
        for b, (m2, n2) in enumerate(pairs):
            if b < a:
                continue
            key = tuple(sorted([tuple(sorted((m1, n1))), tuple(sorted((m2, n2)))]))
            if key not in cache:
                cache[key] = w_entry(m1, n1, m2, n2)
            Wx[a, b] = Wx[b, a] = cache[key]
    assert Wx[0, 0] == sp.Rational(3, 8)
    lam = sp.symbols("lam")
    cp = sp.cancel(Wx.charpoly(lam).as_expr() / lam)      # one trivial zero mode
    cubic = sp.Poly(sp.expand(131072 * cp), lam)
    target = sp.Poly(131072 * lam ** 3 - 50688 * lam ** 2 - 25200 * lam - 675, lam)
    assert cubic == target
    factored = sp.factor_list(target.as_expr())
    degs = sorted(sp.Poly(f, lam).degree() for f, _ in factored[1])
    assert degs == [3], f"cubic unexpectedly factored: {degs}"


# ---------------------------------------------------------------------------
# l > 0: the Gaunt-coupled s+p payoff (rem:ee_partial_split, third paragraph;
# provenance debug/ee_split_sp_fci.py + debug/ee_split_sp_sweep.py)
# ---------------------------------------------------------------------------
def _sp_system(ns, npp, k=2.0, Z=2.0, Ng=500, Lmax=2):
    """Self-contained Gaunt-coupled s+p one-centre system in REAL harmonics."""
    from geovac.xtc_angular_sparsity import gA, gB
    r, wr = TC.make_grid(k, Ng=Ng)
    W2 = r * r * wr
    orbs = [(n, 0, 0) for n in range(1, ns + 1)]
    orbs += [(n, 1, m) for n in range(2, 2 + npp) for m in (-1, 0, 1)]
    nb = len(orbs)

    def Rf(n, l):
        x = 2 * k * r
        nrm = np.sqrt((2 * k) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
        return nrm * x ** l * np.exp(-x / 2) * eval_genlaguerre(n - l - 1, 2 * l + 1, x)

    def dRf(n, l):
        x = 2 * k * r
        nrm = np.sqrt((2 * k) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
        m = n - l - 1
        L = eval_genlaguerre(m, 2 * l + 1, x)
        dL = -eval_genlaguerre(m - 1, 2 * l + 2, x) if m >= 1 else np.zeros_like(x)
        e = np.exp(-x / 2)
        xlm1 = x ** (l - 1) if l >= 1 else np.zeros_like(x)
        return 2 * k * nrm * e * (l * xlm1 * L - 0.5 * x ** l * L + x ** l * dL)

    rad = {(n, l): Rf(n, l) for (n, l, _) in orbs}
    drad = {(n, l): dRf(n, l) for (n, l, _) in orbs}
    S = np.zeros((nb, nb))
    h1 = np.zeros((nb, nb))
    for a, (na, la, ma) in enumerate(orbs):
        for b, (nbb, lb, mb) in enumerate(orbs):
            if (la, ma) != (lb, mb):
                continue
            Ra, Rb = rad[(na, la)], rad[(nbb, lb)]
            S[a, b] = np.sum(Ra * Rb * W2)
            h1[a, b] = (0.5 * np.sum(drad[(na, la)] * drad[(nbb, lb)] * W2)
                        + 0.5 * la * (la + 1) * np.sum(Ra * Rb * wr)
                        - Z * np.sum(Ra * Rb * r * wr))
    R1g, R2g = np.meshgrid(r, r, indexing="ij")
    lo, hi = np.minimum(R1g, R2g), np.maximum(R1g, R2g)
    rps = sorted({(min((na, la), (nbb, lb)), max((na, la), (nbb, lb)))
                  for (na, la, _) in orbs for (nbb, lb, _) in orbs})
    rpi = {p: i for i, p in enumerate(rps)}
    P = np.array([rad[p0] * rad[p1] * W2 for (p0, p1) in rps])
    RL, RLs, RLw = {}, {}, {}
    for L in range(Lmax + 1):
        t1, t2 = R1g ** L / R2g ** (L + 1), R2g ** L / R1g ** (L + 1)
        RL[L] = P @ (lo ** L / hi ** (L + 1)) @ P.T
        RLs[L] = P @ (0.5 * (t1 + t2)) @ P.T
        RLw[L] = P @ (0.5 * np.abs(t1 - t2)) @ P.T
    idx = np.zeros((nb, nb), dtype=int)
    for a, (na, la, _) in enumerate(orbs):
        for c, (nc, lc, _) in enumerate(orbs):
            idx[a, c] = rpi[(min((na, la), (nc, lc)), max((na, la), (nc, lc)))]
    U = np.zeros((nb, nb), dtype=complex)
    for a, (na, la, ma) in enumerate(orbs):
        for b, (nbb, lb, mb) in enumerate(orbs):
            if (na, la) != (nbb, lb):
                continue
            if ma == 0 and mb == 0:
                U[a, b] = 1.0
            elif ma > 0:
                U[a, b] = (((-1) ** ma) / np.sqrt(2) if mb == ma
                           else (1 / np.sqrt(2) if mb == -ma else 0))
            elif ma < 0:
                mm = -ma
                U[a, b] = (-1j * ((-1) ** mm) / np.sqrt(2) if mb == mm
                           else (1j / np.sqrt(2) if mb == -mm else 0))

    def asm(RLd):
        g = np.zeros((nb,) * 4, dtype=complex)
        for L in range(Lmax + 1):
            RLf = RLd[L][idx[:, :, None, None], idx[None, None, :, :]]
            for M in range(-L, L + 1):
                A = np.array([[gA(la, ma, L, M, lc, mc) for (_, lc, mc) in orbs]
                              for (_, la, ma) in orbs])
                B = np.array([[gB(lb, mb, L, M, ld, md) for (_, ld, md) in orbs]
                              for (_, lb, mb) in orbs])
                if np.abs(A).max() < 1e-14 or np.abs(B).max() < 1e-14:
                    continue
                g += (4 * np.pi / (2 * L + 1)) * np.einsum(
                    "ac,bd,acbd->abcd", A, B, RLf, optimize=True)
        gr = np.einsum("ap,bq,cr,ds,pqrs->abcd", U.conj(), U.conj(), U, U, g,
                       optimize=True)
        return np.real(gr), float(np.abs(gr.imag).max())

    Sr = np.real(np.einsum("ap,bq,pq->ab", U.conj(), U, S, optimize=True))
    hr = np.real(np.einsum("ap,bq,pq->ab", U.conj(), U, h1, optimize=True))
    return orbs, Sr, hr, asm, RL, RLs, RLw, len(rps)


def _fci_sp(S, h1, g, n_elec=2):
    X = TC.lowdin(S)
    nso = 2 * S.shape[0]
    dets, didx = TC.make_dets(nso, n_elec)
    H = TC.build_H(dets, didx, TC.h_spin(TC.transform_1(h1, X), nso),
                   TC.asym_from_phys(TC.transform_2(g, X), nso), nso)
    return float(eigh(H, eigvals_only=True)[0])


def test_sp_engine_reproduces_s_only_and_is_real_symmetric():
    """s-only sector == the independent engine; s+p tensor real + 8-fold symmetric."""
    ns = 3
    orbs, S, h1, asm, RL, RLs, RLw, nrp = _sp_system(ns, 0, Lmax=0)
    g, im = asm(RL)
    r, wr = TC.make_grid(2.0, Ng=500)
    Sref, href, Rtab, W2 = TC.build_one_body(ns, r, wr, 2.0, 2.0)
    gref, _, _ = TC.two_body(ns, Rtab, W2, TC.build_kernels(r, 0.7, nx=64))
    assert np.abs(S - Sref).max() < 1e-12
    assert np.abs(h1 - href).max() < 1e-12          # analytic dR/dr required
    assert np.abs(g - gref).max() < 1e-12
    assert abs(_fci_sp(S, h1, g) - _fci_sp(Sref, href, gref)) < 1e-12
    # s+p: reality and permutational symmetry (catches a transposed real-Y transform)
    orbs, S, h1, asm, RL, RLs, RLw, nrp = _sp_system(3, 1, Lmax=2)
    g, im = asm(RL)
    assert im < 1e-12, f"tensor not real: {im:.1e}"
    for perm in ((2, 1, 0, 3), (0, 3, 2, 1), (1, 0, 3, 2)):
        assert np.abs(g - g.transpose(perm)).max() < 1e-10


def test_sp_split_exact_and_rank_payoff():
    """Split exact per channel; W-rank 4 reaches sub-mHa with L=0,1,2 all active;
    rank 0 (separable only) is off by >100 mHa."""
    orbs, S, h1, asm, RL, RLs, RLw, nrp = _sp_system(3, 1, Lmax=2)
    g, _ = asm(RL)
    gs, _ = asm(RLs)
    gw, _ = asm(RLw)
    assert np.abs(g - (gs - gw)).max() / np.abs(g).max() < 1e-12
    E_full = _fci_sp(S, h1, g)
    assert E_full < _fci_sp(*_sp_system(3, 0, Lmax=0)[1:3],
                            _sp_system(3, 0, Lmax=0)[3](
                                _sp_system(3, 0, Lmax=0)[4])[0])   # below s-only
    assert E_full > -2.9037243770                                   # above exact
    assert abs(_fci_sp(S, h1, gs) - E_full) > 0.1                   # teeth
    RLt = {}
    for L in range(3):
        lam, V = np.linalg.eigh(RLw[L])
        o = np.argsort(-np.abs(lam))[:4]
        RLt[L] = (V[:, o] * lam[o]) @ V[:, o].T
    gwt, _ = asm(RLt)
    assert abs(_fci_sp(S, h1, gs - gwt) - E_full) < 1e-3


def test_1norm_fold_and_degradation():
    """The exact one-body fold lowers the JW LCU 1-norm at bit-identical energy --
    AND the advantage degrades with basis size.  Both pinned: the degradation is
    part of the claim, so it cannot drift into an unqualified '15% lever'."""
    ratios = {}
    for ns in (4, 6, 8, 10, 12):
        S, h1, Vr, g, gs, gw = _build(ns, 2.0, Z=2.0, Ng=500)
        X = TC.lowdin(S)
        nso = 2 * ns
        dets, didx = TC.make_dets(nso, 2)

        def run(hx, gx):
            hso = TC.h_spin(TC.transform_1(hx, X), nso)
            asym = TC.asym_from_phys(TC.transform_2(gx, X), nso)
            E = float(eigh(TC.build_H(dets, didx, hso, asym, nso), eigvals_only=True)[0])
            return TC.lcu_lambda(hso, asym, nso)["lam"], E

        lamA, EA = run(h1, g)
        lamB, EB = run(h1 + 0.5 * Vr, -gw)          # N = 2  =>  (N-1)/2 = 1/2
        assert abs(EA - EB) < 1e-9, "the fold must be exact"
        ratios[ns] = lamB / lamA
        assert ratios[ns] < 1.0, f"ns={ns}: no 1-norm win ({ratios[ns]:.3f})"
    # TWO-SIDED per-ns pins (cert-2 M-1: the paper prints all five values
    # 0.83/0.85/0.89/0.91/0.94; the old one-sided ratios[4] < 0.90 +
    # ratios[8] < 1.0 admitted a 10% win at ns=4 and a 0.1% win at ns=8).
    expected_ratio = {4: 0.8278, 6: 0.8521, 8: 0.8897, 10: 0.9145, 12: 0.9367}
    for ns, want in expected_ratio.items():
        assert abs(ratios[ns] - want) < 5e-3, (
            f"ns={ns}: ratio {ratios[ns]:.4f} != {want}")
    # monotone degradation toward unity
    seq = [ratios[ns] for ns in (4, 6, 8, 10, 12)]
    assert all(b > a for a, b in zip(seq, seq[1:])), f"not monotone: {seq}"


# ---------------------------------------------------------------------------
# Coverage closure (2026-08-28 /qa delta, code dimension MATERIAL-3): these
# sub-claims of rem:ee_partial_split were cited to THIS file but asserted only
# in untracked debug/ drivers.  Paper claims must be regression-protected here.
# ---------------------------------------------------------------------------
def test_exact_rational_anchors_5k8_and_3_5():
    """<1s^2|g|1s^2> = 5k/8 and <1s^2|W>/<1s^2|g> = 3/5, both exact."""
    for k in (1.0, 2.0):
        _, _, _, g, gs, gw = _build(4, k, Z=2.0, Ng=700)
        # RELATIVE tolerance: the residual is grid quadrature and scales with k.
        # 1e-4 still discriminates hugely -- the nearest wrong rationals (4k/8,
        # 6k/8) are 20% away, i.e. 2000x outside this bound.
        assert abs(g[0, 0, 0, 0] / (5 * k / 8) - 1) < 1e-4, "5k/8 anchor"
        assert abs(gw[0, 0, 0, 0] / g[0, 0, 0, 0] - 0.6) < 1e-4, "3/5 ratio"
        assert abs(gw[0, 0, 0, 0] / (3 * k / 8) - 1) < 1e-4, "3k/8 anchor"


def test_W_frobenius_tail_improves_with_L():
    """Compressibility improves with L: the Frobenius tail after four terms
    falls monotonically across L = 0..3 (paper quotes 1.6/0.7/0.1/0.03%)."""
    from math import factorial
    k, Ng, nfun = 2.0, 700, 5
    r, wr = TC.make_grid(k, Ng=Ng)
    Wt = r * r * wr
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    tails = []
    for L in range(4):
        ns_list = list(range(L + 1, L + 1 + nfun))
        R = {n: _R_nl(n, L, r, k) for n in ns_list}
        D = np.array([R[a] * R[b] * Wt for a in ns_list for b in ns_list])
        t1, t2 = R1 ** L / R2 ** (L + 1), R2 ** L / R1 ** (L + 1)
        M = D @ (0.5 * np.abs(t1 - t2)) @ D.T
        lam = np.sort(np.abs(np.linalg.eigvalsh(M)))[::-1]
        tails.append(float(np.sqrt(np.sum(lam[4:] ** 2) / np.sum(lam ** 2))))
    # paper quotes 1.6%, 0.7%, 0.1%, 0.03% -- pin each (cert-2 N-12: the
    # old tails[0] < 0.05 admitted [4.9, 4.0, 3.0, 2.0]%)
    for t, want in zip(tails, (0.0163, 0.0065, 0.0014, 0.00031)):
        assert abs(t - want) / want < 0.30, f"tails {tails} vs paper"
    for a, b in zip(tails, tails[1:]):
        assert b < a, f"tail not improving with L: {tails}"


def test_ode_lambda_fpp_equals_w_f():
    """In u = 1/r the W-kernel is |u1-u2|/2, so eigenfunctions of any weighted
    compression solve lambda f'' = w f -- the surviving closed form."""
    U, Nu = 12.0, 1500
    u = np.linspace(0.0, U, Nu)
    du = u[1] - u[0]
    w = np.exp(-u)
    T = (0.5 * np.abs(u[:, None] - u[None, :])) * (w[None, :] * du)
    ev, V = np.linalg.eig(T)
    order = np.argsort(-np.abs(ev))
    for j in order[:2]:
        lam_j = float(np.real(ev[j]))
        f = np.real(V[:, j])
        fpp = np.gradient(np.gradient(f, du), du)
        lhs, rhs = lam_j * fpp[5:-5], (w * f)[5:-5]
        res = np.linalg.norm(lhs - rhs) / np.linalg.norm(rhs)
        assert res < 5e-5, f"ODE residual {res:.1e}"  # measured 3.0-4.6e-5


def test_sp_rank_flatness_two_basis_sizes():
    """'flat as the radial-pair space grows 10 -> 28': the pytest file tested
    only the low end.  Check that the 1 mHa rank does not grow with basis."""
    ranks = {}
    for ns, npp in ((3, 1), (4, 2), (5, 2)):  # nrp = 10, 21, 28
        orbs, S, h1, asm, RL, RLs, RLw, nrp = _sp_system(ns, npp, Lmax=2)
        g, _ = asm(RL)
        gs, _ = asm(RLs)
        E_full = _fci_sp(S, h1, g)
        got = None
        for m in range(1, nrp + 1):
            RLt = {}
            for L in range(3):
                lam, V = np.linalg.eigh(RLw[L])
                o = np.argsort(-np.abs(lam))[:m]
                RLt[L] = (V[:, o] * lam[o]) @ V[:, o].T
            gwt, _ = asm(RLt)
            if abs(_fci_sp(S, h1, gs - gwt) - E_full) * 1000 < 1.0:
                got = m
                break
        ranks[nrp] = got
        # paper claims rank 3-4; the old <= 6 admitted 5-6 (cert-2 M-5)
        assert got is not None and got <= 4, f"nrp={nrp}: rank {got}"
    ks = sorted(ranks)
    assert ranks[ks[-1]] <= ranks[ks[0]] + 1, f"rank not flat: {ranks}"


def test_split_beats_density_fitting_control():
    """The 13-17% 'beyond DF' claim: at matched rank the split's lambda must be
    below a plain eigen-truncation of g itself (the density-fitting control)."""
    ns, k = 4, 2.0
    S, h1, Vr, g, gs, gw = _build(ns, k, Z=2.0, Ng=600)
    X = TC.lowdin(S)
    nso = 2 * ns
    dets, didx = TC.make_dets(nso, 2)

    def lam_E(hx, gx):
        hso = TC.h_spin(TC.transform_1(hx, X), nso)
        asym = TC.asym_from_phys(TC.transform_2(gx, X), nso)
        E = float(eigh(TC.build_H(dets, didx, hso, asym, nso), eigvals_only=True)[0])
        return TC.lcu_lambda(hso, asym, nso)["lam"], E

    def trunc(t, m):
        M = t.transpose(0, 2, 1, 3).reshape(ns * ns, ns * ns)
        lam, U = np.linalg.eigh(M)
        o = np.argsort(-np.abs(lam))[:m]
        r = ((U[:, o] * lam[o]) @ U[:, o].T).reshape(ns, ns, ns, ns).transpose(0, 2, 1, 3)
        return 0.25 * (r + r.transpose(2, 1, 0, 3) + r.transpose(0, 3, 2, 1)
                       + r.transpose(2, 3, 0, 1))

    _, E_ex = lam_E(h1, g)
    lamC, EC = lam_E(h1 + 0.5 * Vr, -trunc(gw, 4))     # split, W at rank 4
    lamD, ED = lam_E(h1, trunc(g, 4))                  # DF control, g at rank 4
    assert abs(EC - E_ex) < 5e-4 and abs(ED - E_ex) < 5e-4, "not matched accuracy"
    assert lamC < lamD, f"split lambda {lamC:.2f} not below DF control {lamD:.2f}"
    # tightened 2026-08-28: paper claims 13-17%, measured 16.9%; 0.95 was 3x loose
    # TWO-SIDED 2026-08-28 (delta-2): the paper claims a 13-17% margin and the
    # measured ratio is 0.8311 (16.9%), stable to 1e-4 over Ng=400..1200.  A
    # one-sided `< 0.88` admitted a 12% margin -- BELOW the paper's own 13%
    # floor -- and let an implausibly good ratio (a control-side bug) pass.
    assert 0.82 < lamC / lamD < 0.87, (
        f"DF margin {1 - lamC / lamD:.1%} outside the paper 13-17% band "
        f"(ratio {lamC / lamD:.4f})"
    )
