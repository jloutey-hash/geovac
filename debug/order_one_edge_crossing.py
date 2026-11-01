"""Order-one defect of the true CH window as an edge-crossing sum.

Claim tested (walls register, SPECTRAL-WINDOW-COMPRESSION, owed step):
    C(a,b) := [[P D P, P a P], P b P]     (exactly the test's order-one quantity,
                                         tests/test_true_window_real_structure.py::_stats,
                                         with J b^* J^{-1} = P b P, asserted there)
equals, for i, l in W,
    C_il = sum_{k not in W} [ -(lam_i - lam_k) a_ik b_kl + b_ik (lam_k - lam_l) a_kl ],
because [D, P] = 0 and the continuum identity [[D, a], b] = 0 (Clifford multiplication by
da commutes with the function b).

Also tested: the gauge companion of tests/test_paper32_window_gauge_breaking.py is
    k_pg = (1/2) Tr_W [ f'(D_W) C(chi, chi) ],   chi = i(b - b^dag), f = exp(-t x^2).

Run with OMP/OPENBLAS/MKL_NUM_THREADS=1.  Ambient: j <= (3 n_max - 1)/2, which contains every
state a band <= 2 n_max multiplier reaches from the window (so the crossing sum is complete).
"""
from __future__ import annotations

import itertools
import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from math import factorial, sqrt

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(os.path.dirname(HERE), "tests")
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


# ------------------------------------------------------------------ Clebsch-Gordan (Racah)
@lru_cache(maxsize=None)
def cg2(j1, m1, j2, m2, J, M):
    """<j1 m1; j2 m2 | J M>, all arguments doubled integers."""
    if m1 + m2 != M or abs(m1) > j1 or abs(m2) > j2 or abs(M) > J:
        return 0.0
    if J < abs(j1 - j2) or J > j1 + j2 or (j1 + j2 + J) % 2:
        return 0.0
    f = lambda x: factorial(x // 2)
    pre = sqrt((J + 1) * f(J + j1 - j2) * f(J - j1 + j2) * f(j1 + j2 - J) / f(j1 + j2 + J + 2))
    pre *= sqrt(f(J + M) * f(J - M) * f(j1 - m1) * f(j1 + m1) * f(j2 - m2) * f(j2 + m2))
    s = 0.0
    for k in range(0, j1 + j2 + J + 1):
        args = [j1 + j2 - J - 2 * k, j1 - m1 - 2 * k, j2 + m2 - 2 * k,
                J - j2 + m1 + 2 * k, J - j1 - m2 + 2 * k]
        if min(args) < 0:
            continue
        den = factorial(k)
        for a in args:
            den *= factorial(a // 2)
        s += (-1) ** k / den
    return pre * s


def check_cg():
    from sympy.physics.quantum.cg import CG
    from sympy import Rational as R
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(60):
        j1, j2 = rng.integers(0, 8, 2)
        J = rng.choice(np.arange(abs(j1 - j2), j1 + j2 + 1, 2))
        m1 = rng.choice(np.arange(-j1, j1 + 1, 2)); m2 = rng.choice(np.arange(-j2, j2 + 1, 2))
        if abs(m1 + m2) > J:
            continue
        ref = float(CG(R(int(j1), 2), R(int(m1), 2), R(int(j2), 2), R(int(m2), 2),
                       R(int(J), 2), R(int(m1 + m2), 2)).doit())
        worst = max(worst, abs(ref - cg2(int(j1), int(m1), int(j2), int(m2), int(J), int(m1 + m2))))
    return worst


# ------------------------------------------------------------------ construction
class Ambient:
    """Labels (j, m, m', s) with j <= j2max/2 (doubled ints), Dirac eigenbasis per j block."""

    def __init__(self, j2max):
        self.labels = [(j2, m2, mp2, s) for j2 in range(j2max + 1)
                       for m2 in range(j2, -j2 - 1, -2) for mp2 in range(j2, -j2 - 1, -2) for s in (0, 1)]
        self.idx = {l: i for i, l in enumerate(self.labels)}
        n = len(self.labels)
        self.n = n
        V = np.zeros((n, n), complex)
        lam = np.zeros(n)
        for j2 in range(j2max + 1):
            d = j2 + 1
            Sm = spin_mats(j2)
            blk = 2 * (np.kron(Sm[0], SX) + np.kron(Sm[1], SY) + np.kron(Sm[2], SZ)) + 1.5 * np.eye(2 * d)
            w, U = np.linalg.eigh(blk)
            for a, m2 in enumerate(range(j2, -j2 - 1, -2)):
                ids = [self.idx[(j2, m2, mp2, s)] for mp2 in range(j2, -j2 - 1, -2) for s in (0, 1)]
                V[np.ix_(ids, ids)] = U
                lam[ids] = w
        self.V, self.lam = V, lam
        self.D = None

    def mult(self, J2, M2, Mp2):
        rows, cols, vals = [], [], []
        for (j1, m1, mp1, s) in self.labels:
            m2, mp2 = m1 + M2, mp1 + Mp2
            for j2 in range(abs(j1 - J2), j1 + J2 + 1, 2):
                if (j2, m2, mp2, s) not in self.idx:
                    continue
                c = cg2(j1, m1, J2, M2, j2, m2) * cg2(j1, mp1, J2, Mp2, j2, mp2)
                if c != 0.0:
                    rows.append(self.idx[(j2, m2, mp2, s)])
                    cols.append(self.idx[(j1, m1, mp1, s)])
                    vals.append(sqrt((j1 + 1) / (j2 + 1)) * c)
        return sp.csr_matrix((vals, (rows, cols)), shape=(self.n, self.n), dtype=complex)


def spin_mats(j2):
    ms = [j2 / 2 - k for k in range(j2 + 1)]
    j = j2 / 2
    jp = np.zeros((j2 + 1, j2 + 1), complex)
    for a, m in enumerate(ms):
        if a > 0:
            jp[a - 1, a] = np.sqrt(j * (j + 1) - m * (m + 1))
    return (jp + jp.conj().T) / 2, (jp - jp.conj().T) / (2j), np.diag(ms).astype(complex)


def band_labels(n_max):
    out = []
    for N in range(1, 2 * n_max + 1):
        J2 = N - 1
        Ms = [J2 - 2 * k for k in range(N)]
        for M2, Mp2 in itertools.product(Ms, Ms):
            out.append((N, J2, M2, Mp2))
    return out


def run(n_max, Dabs=False):
    A = Ambient(3 * n_max - 1)
    lam = A.lam
    W = np.nonzero(np.abs(lam) <= n_max + 0.5 + 1e-9)[0]
    Q = np.nonzero(np.abs(lam) > n_max + 0.5 + 1e-9)[0]
    VW, VQ = A.V[:, W], A.V[:, Q]
    lW, lQ = lam[W], lam[Q]
    labs = band_labels(n_max)
    TWW, TWQ, TQW = [], [], []
    for (_, J2, M2, Mp2) in labs:
        X = A.mult(J2, M2, Mp2)
        XVW = X @ VW                       # columns into W
        XhVW = X.getH() @ VW               # for rows W
        TWW.append(VW.conj().T @ XVW)
        TQW.append(VQ.conj().T @ XVW)
        TWQ.append((VQ.conj().T @ XhVW).conj().T)
    TWW, TWQ, TQW = np.array(TWW), np.array(TWQ), np.array(TQW)
    wW = np.abs(lW) if Dabs else lW
    wQ = np.abs(lQ) if Dabs else lQ
    # direct: C = [[D_W, T_a], T_b]
    XW = wW[None, :, None] * TWW - TWW * wW[None, None, :]
    # crossing pieces: [D,a]_{WQ} = (w_i - w_k) a_ik ; [D,a]_{QW} = (w_k - w_l) a_kl
    gWQ = wW[:, None] - wQ[None, :]
    gQW = wQ[:, None] - wW[None, :]
    n = len(labs)
    rel, absd, cmax = np.zeros((n, n)), np.zeros((n, n)), np.zeros((n, n))
    for i in range(n):
        Cd = XW[i] @ TWW - TWW @ XW[i]
        Cf = -((gWQ * TWQ[i]) @ TQW) + TWQ @ (gQW * TQW[i])
        d = np.abs(Cd - Cf).reshape(n, -1).max(1)
        m = np.abs(Cd).reshape(n, -1).max(1)
        absd[i], cmax[i] = d, m
        rel[i] = np.where(m > 1e-12, d / np.maximum(m, 1e-300), 0.0)
    return dict(A=A, W=W, Q=Q, lW=lW, labs=labs, TWW=TWW, TWQ=TWQ, TQW=TQW,
                rel=rel, absd=absd, cmax=cmax)


def continuum_identity(R, Dabs=False):
    """max over pairs of |([[D,a],b])_{WW}| with the sum over ALL ambient intermediates."""
    lW = R["lW"]; lQ = R["A"].lam[R["Q"]]
    wW = np.abs(lW) if Dabs else lW
    wQ = np.abs(lQ) if Dabs else lQ
    TWW, TWQ, TQW = R["TWW"], R["TWQ"], R["TQW"]
    n = len(TWW)
    worst, scale = 0.0, 0.0
    for i in range(n):
        DaWW = (wW[:, None] - wW[None, :]) * TWW[i]
        DaWQ = (wW[:, None] - wQ[None, :]) * TWQ[i]
        DaQW = (wQ[:, None] - wW[None, :]) * TQW[i]
        # ([D,a] b - b [D,a])_{WW} = DaWW bWW + DaWQ bQW - bWW DaWW - bWQ DaQW
        M = DaWW @ TWW + DaWQ @ TQW - TWW @ DaWW - TWQ @ DaQW
        worst = max(worst, float(np.abs(M).max()))
        scale = max(scale, float(np.abs(DaWW).max()))
    return worst, scale


def test_module_crosscheck(n_max, R):
    """Compare the crossing formula, mapped to the TEST's ambient label basis, to the test's
    own V @ C @ V^dag for every pair (elementwise), with C built from the test's Bs."""
    sys.path.insert(0, TESTS)
    import test_true_window_real_structure as T
    S = T._setup(n_max)
    Ts, Bs, V, Dc = S["Ts"], S["Bs"], S["V"], S["Dc"]
    tl = S["labels"]
    # map test labels (Fractions) -> my ambient index
    A = R["A"]
    rows = [A.idx[(int(2 * j), int(2 * m), int(2 * mp), s)] for (j, m, mp, s) in tl]
    VWm = A.V[rows][:, R["W"]]           # my window basis expressed on the test labels
    lW = R["lW"]; lQ = A.lam[R["Q"]]
    DaWQ_f = lambda i: (lW[:, None] - lQ[None, :]) * R["TWQ"][i]
    DaQW_f = lambda i: (lQ[:, None] - lW[None, :]) * R["TQW"][i]
    n = len(Ts)
    worst_rel = 0.0
    Xs = Dc @ Ts - Ts @ Dc
    for i in range(n):
        Ct = Xs[i] @ Bs - Bs @ Xs[i]
        At = V @ Ct @ V.conj().T                 # test's ambient-basis operator
        Cf = -np.einsum("ik,bkl->bil", DaWQ_f(i), R["TQW"]) + np.einsum("bik,kl->bil", R["TWQ"], DaQW_f(i))
        Af = VWm @ Cf @ VWm.conj().T
        m = np.abs(At).reshape(n, -1).max(1)
        d = np.abs(At - Af).reshape(n, -1).max(1)
        sel = m > 1e-12
        if sel.any():
            worst_rel = max(worst_rel, float((d[sel] / m[sel]).max()))
        assert np.all(d[~sel] < 1e-12)
    ma, op = T._stats(n_max, 1)
    return worst_rel, float(ma.max()), float(op.max()), int((ma > 1e-10).sum())


def gauge_identity():
    """k_pg (test's Duhamel eps^2 coefficient) vs (1/2) Tr_W f'(D_W) [[D_W, chi_W], chi_W]."""
    sys.path.insert(0, TESTS)
    import test_paper32_window_gauge_breaking as G
    out = []
    for tau in (1.0, 4.0):
        for n in (2, 3, 4, 5):
            t = tau / (n + 0.5) ** 2
            _, kpg = G.window_k(n, t)
            S, a, b = G.sphere(2 * n + 2)
            w = np.nonzero(np.abs(S.lam) <= n + 0.5)[0]
            chi = (1j * (b - b.getH()))[w][:, w].toarray()
            l = S.lam[w]
            Dchi = l[:, None] * chi - chi * l[None, :]
            C = Dchi @ chi - chi @ Dchi
            fp = -2 * t * l * np.exp(-t * l ** 2)
            tr = 0.5 * float(np.real(np.sum(fp * np.diag(C))))
            # control: |D| in place of D inside the double commutator (the sign(D) confusion)
            la = np.abs(l)
            Da = la[:, None] * chi - chi * la[None, :]
            Ca = Da @ chi - chi @ Da
            ctrl = 0.5 * float(np.real(np.sum(fp * np.diag(Ca))))
            out.append(dict(tau=tau, n=n, kpg=kpg, half_trace=tr, rel=abs(tr - kpg) / abs(kpg),
                            absD_half_trace=ctrl))
    return out


def main():
    res = {"cg_max_err_vs_sympy": check_cg(), "rows": []}
    nmax_list = [int(x) for x in sys.argv[1:]] or [1, 2, 3]
    for n in nmax_list:
        R = run(n)
        Rb = run(n, Dabs=True) if n <= 3 else None
        ci, sc = continuum_identity(R)
        cia = continuum_identity(R, Dabs=True)[0] if n <= 3 else None
        row = dict(n_max=n, dimW=len(R["W"]), dimQ=len(R["Q"]), pairs=int(R["rel"].size),
                   max_rel=float(R["rel"].max()), max_abs=float(R["absd"].max()),
                   max_defect=float(R["cmax"].max()),
                   absD_control_max_mismatch=None if Rb is None else float(Rb["absd"].max()),
                   absD_control_max_defect=None if Rb is None else float(Rb["cmax"].max()),
                   continuum_identity_D=ci, continuum_identity_absD=cia, gap_scale=sc)
        if n <= 3:
            wr, ma, op, nf = test_module_crosscheck(n, R)
            row.update(test_crosscheck_max_rel=wr, test_maxabs=ma, test_opnorm=op, test_nfail=nf)
        res["rows"].append(row)
        print(json.dumps(row), flush=True)
    if os.environ.get("SKIP_GAUGE"):
        res["gauge"] = []
    else:
        res["gauge"] = gauge_identity()
    for g in res["gauge"]:
        print(json.dumps(g), flush=True)
    out = os.path.join(HERE, "data", "order_one_edge_crossing_%s.json" % "_".join(map(str, nmax_list)))
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    print("cg err", res["cg_max_err_vs_sympy"])


if __name__ == "__main__":
    main()
