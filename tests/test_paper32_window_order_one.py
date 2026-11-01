"""Paper 32: the order-one defect on the true CH window is an exact edge-crossing sum.

Written as a separate step (CLAUDE.md Sec. 9) after a sub-agent derived and checked the identity
(debug/order_one_edge_crossing_memo.md; PM re-derived it).  With W the window, P its projection, Q = 1 - P and
D diagonal on Dirac eigenstates, [PDP, PaP] = P[D,a]P, so
    [[PDP, PaP], PbP] = P[[D,a],b]P - P[D,a] Q b P + P b Q [D,a] P.
For functions a, b the first term vanishes ([D,a] is Clifford multiplication by da and commutes with b: the
continuum order-one condition), so every surviving term passes through one state outside W:
    C_il = sum_{k not in W} [ -(lambda_i - lambda_k) a_ik b_kl + b_ik (lambda_k - lambda_l) a_kl ].
The gauge defect of tests/test_paper32_window_gauge_breaking.py is the same object at a = b = chi.

No debug/ import: the S^3 construction is the one ported into the gauge-breaking test.
"""
import importlib.util
import os

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('gb_for_o1', os.path.join(HERE, 'test_paper32_window_gauge_breaking.py'))
gb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gb)

PAIRS = [(-1, 1), (1, 1), (1, -1), (-1, -1)]      # the four D^{1/2}_{M M'} (doubled M, M')


def setup(n_max):
    S = gb.Sphere(2 * n_max + 2)
    funcs = [S.multiplier(M, Mp).toarray() for M, Mp in PAIRS]
    inW = np.abs(S.lam) <= n_max + 0.5
    return S, funcs, inW


def defect_direct(lam, a, b, inW):
    w = np.nonzero(inW)[0]
    Dw = np.diag(lam[w])
    aw, bw = a[np.ix_(w, w)], b[np.ix_(w, w)]
    Da = Dw @ aw - aw @ Dw
    return Da @ bw - bw @ Da


def defect_edge(lam, a, b, inW):
    w, out = np.nonzero(inW)[0], np.nonzero(~inW)[0]
    li, lk = lam[w], lam[out]
    a_wk, b_kw = a[np.ix_(w, out)], b[np.ix_(out, w)]
    b_wk, a_kw = b[np.ix_(w, out)], a[np.ix_(out, w)]
    t1 = -((li[:, None] - lk[None, :]) * a_wk) @ b_kw
    t2 = b_wk @ ((lk[:, None] - li[None, :]) * a_kw)
    return t1 + t2


def test_continuum_order_one_holds_on_window_rows():
    """P[[D,a],b]P = 0 for functions a, b (all intermediate states kept): the identity the edge form rests on.
    With |D| in place of D it fails at order one.

    Wrong answer rejected: |D| (the sign(D) confusion) in place of the Dirac operator.
    """
    for n_max in (1, 2, 3):
        S, funcs, inW = setup(n_max)
        w = np.nonzero(inW)[0]
        for lam, expect_zero in ((S.lam, True), (np.abs(S.lam), False)):
            worst = 0.0
            for a in funcs:
                for b in funcs:
                    Da = lam[:, None] * a - a * lam[None, :]
                    worst = max(worst, np.abs((Da @ b - b @ Da)[np.ix_(w, w)]).max())
            if expect_zero:
                assert worst < 1e-12, (n_max, worst)
            else:
                assert worst > 0.1, (n_max, worst)


def test_window_order_one_defect_is_the_edge_crossing_sum():
    """On the true window, [[D_W, a_W], b_W] equals the edge-crossing sum entry by entry, at n_max = 1, 2, 3 for
    all sixteen pairs of the four spin-1/2 coefficient functions, and the defect is not zero.

    Wrong answers rejected: an interior sum (k in W) in place of the crossing one, the defect vanishing on the
    window, and the gap weight evaluated with |lambda|.
    """
    for n_max in (1, 2, 3):
        S, funcs, inW = setup(n_max)
        biggest = 0.0
        for a in funcs:
            for b in funcs:
                C = defect_direct(S.lam, a, b, inW)
                E = defect_edge(S.lam, a, b, inW)
                scale = max(np.abs(C).max(), 1e-300)
                biggest = max(biggest, np.abs(C).max())
                assert np.abs(C - E).max() < 1e-12 * max(1.0, scale), (n_max, np.abs(C - E).max())
        assert biggest > 0.1, (n_max, biggest)


def test_gauge_defect_is_the_order_one_defect_at_a_equals_b():
    """With a = b = chi = i(b - b^dag), the gauge-breaking test's pure-gauge heat coefficient equals
    (1/2) Tr_W[ f'(D_W) [[D_W, chi_W], chi_W] ], f = exp(-t x^2): one object, two weights.

    Wrong answer rejected: a gauge defect unrelated to the order-one defect (any other weight or sign).
    """
    for tau, n_max in ((1.0, 2), (4.0, 3)):
        t = tau / (n_max + 0.5) ** 2
        _, kpg = gb.window_k(n_max, t)
        S, a, b = gb.sphere(2 * n_max + 2)
        chi = (1j * (b - b.getH())).toarray()
        inW = np.abs(S.lam) <= n_max + 0.5
        C = defect_direct(S.lam, chi, chi, inW)
        lw = S.lam[inW]
        fprime = -2 * t * lw * np.exp(-t * lw ** 2)
        val = 0.5 * float(np.real(np.sum(fprime * np.diag(C))))
        assert abs(val - kpg) < 1e-10 * max(1.0, abs(kpg)), (tau, n_max, val, kpg)
