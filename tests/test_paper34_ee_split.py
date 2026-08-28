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
