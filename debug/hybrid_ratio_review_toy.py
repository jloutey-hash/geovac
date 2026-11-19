"""Review, attack 1: the quotient count order(f) = d_E + k against a DIRECT minimal-order computation on toys.

Direct order = rank over Q(R) of the coordinate vectors of f, f', ..., f^(N) in the d/dR-closed atom basis
(N = closure size), using tests/_minimal_order.py's atom calculus (derivative, generic rank) -- import only.
Quotient count: d_E = sum over shifts of dim_Q span{q_j} (rank of the q's Wronskian matrix at that shift);
L = Wronskian operator of a basis {g_i e^{c_i R}} of the E-coefficient kernel, computed here independently of
debug/hybrid_ratio_analyse.py; k = number of rates in L f with nonzero coefficient.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import _minimal_order as M  # noqa: E402

R = M.R
FR, RG = M.FR, M.RG


class Toy:
    def __init__(self, terms):
        """terms: list of (coeff expr in R, (c, kind, kappa)) with kind '1' or 'y' (y = Ei(kappa R), or ln R if 0)."""
        self.coords = {}
        for cf, (c, kind, k) in terms:
            key = ((sp.Rational(c), kind, sp.Rational(k)), ())
            self.coords[key] = self.coords.get(key, FR.zero) + FR.from_sympy(sp.sympify(cf))
        self.coords = {k: v for k, v in self.coords.items() if v != 0}


def direct_order(f: Toy) -> int:
    labs = {lb for (lb, _m) in f.coords}
    N = len(M.closure(labs)) + 1
    ds = [f.coords]
    for _ in range(N):
        ds.append(M.derivative(ds[-1]))
    keys = sorted(set().union(*[set(D) for D in ds]), key=str)
    Mx = DomainMatrix([[D.get(k, FR.zero) for D in ds] for k in keys], (len(keys), len(ds)), FR)
    return int(Mx.rank())


def Dc(q, c, j):
    """(d/dR + c)^j q, q in FR."""
    for _ in range(j):
        q = q.diff(RG) + FR.from_sympy(c) * q
    return q


def quotient_count(f: Toy):
    E = {lb: v for (lb, _m), v in f.coords.items() if lb[1] == "y"}
    shifts = sorted({lb[0] for lb in E})
    basis = []      # (c, q) spanning the kernel
    dE = 0
    for c in shifts:
        qs = [v for lb, v in E.items() if lb[0] == c]
        n = len(qs)
        W = DomainMatrix([[q_.diff(RG) if False else Dc(q_, 0, j) for j in range(n)] for q_ in qs], (n, n), FR)
        # rank of the Wronskian matrix of the q's = dim_C span{q}; choose a basis greedily
        chosen = []
        for q_ in qs:
            trial = chosen + [q_]
            m = len(trial)
            Wt = DomainMatrix([[Dc(t, 0, j) for j in range(m)] for t in trial], (m, m), FR)
            if Wt.det() != 0:
                chosen = trial
        dE += len(chosen)
        basis += [(c, q_) for q_ in chosen]
    d = len(basis)
    # L y = det [ rows i: (D_{c_i}^j q_i)_{j=0..d} ; last row: y^(j) ], exponential factor dropped
    ds = [f.coords]
    for _ in range(d):
        ds.append(M.derivative(ds[-1]))
    A = []
    for j in range(d + 1):
        cols = [jj for jj in range(d + 1) if jj != j]
        minor = DomainMatrix([[Dc(q_, c, jj) for jj in cols] for (c, q_) in basis], (d, d), FR).det() if d else FR.one
        A.append(((-1) ** (d + j)) * minor)
    Lf = {}
    for j in range(d + 1):
        for key, v in ds[j].items():
            Lf[key] = Lf.get(key, FR.zero) + A[j] * v
    Lf = {k: v for k, v in Lf.items() if v != 0}
    assert all(lb[1] == "1" for (lb, _m) in Lf), f"L did not kill the E-part: {list(Lf)}"
    k = len({lb[0] for (lb, _m) in Lf})
    return dE, k, sorted(str(lb[0]) for (lb, _m) in Lf)


def run():
    q = R ** 2 + 1
    qm = q.subs(R, -R)
    cases = {
        # one shift, proportional E-coefficients, generic exponentials
        "T1 generic": [(q, (-1, "y", 1)), (-q, (-1, "y", 2)), (R, (-1, "1", 0)), (1, (-3, "1", 0))],
        # exponential coefficient IN L's kernel (r = 5q): rate -1 must drop
        "T2 r in kernel": [(q, (-1, "y", 1)), (-q, (-1, "y", 2)), (5 * q, (-1, "1", 0)), (1, (-3, "1", 0))],
        # an exponential at a CLOSURE rate (0 = -1 + 1) with a rational coefficient
        "T3 closure-rate exp": [(q, (-1, "y", 1)), (-q, (-1, "y", 2)), (R ** 3, (0, "1", 0))],
        # one shift, two INDEPENDENT coefficients: d_E = 2 from a single component
        "T4 independent at one shift": [(q, (-1, "y", 1)), (R, (-1, "y", 2)), (1, (-3, "1", 0))],
        # three atoms at one shift, coefficients q1, q2, q1 + q2: d_E = 2
        "T5 dependent triple": [(q, (-1, "y", 1)), (R, (-1, "y", 2)), (q + R, (-1, "y", 3)), (R, (-1, "1", 0))],
        # the real shape: two mirror shifts, log atom included, exp coefficient in the kernel at +1 shift? (rate +1
        # is not an exp rate of f; put q e^{-R} (kernel) and generic e^{-2R})
        "T6 mirror pair": [(q, (-1, "y", sp.Rational(-1, 2))), (-q, (-1, "y", sp.Rational(1, 3))),
                           (qm, (1, "y", -3)), (-qm, (1, "y", -sp.Rational(5, 2))),
                           (q, (-1, "1", 0)), (R + 2, (-2, "1", 0)), (1 / R, (-sp.Rational(3, 2), "1", 0))],
        # coincidence-like: log atom and Ei at the same shift, proportional
        "T7 log + Ei": [(q, (-1, "y", 0)), (-q, (-1, "y", 2)), (R, (-1, "1", 0)), (1, (1, "1", 0))],
        # adversarial: exp at rate c + kappa whose coefficient is chosen to make L f's coefficient there small?
        # (cannot vanish: would need w' + kappa w = -1/R rational) -- take w = -1/(kappa R) as near-miss
        "T8 near-miss": [(q, (-1, "y", 1)), (-q * sp.Rational(1, 1) / R * 0 + q * (-1 / R), (0, "1", 0)),
                         (-q, (-1, "y", 2))],
    }
    out = {}
    for name, terms in cases.items():
        f = Toy(terms)
        do = direct_order(f)
        dE, k, rates = quotient_count(f)
        out[name] = {"direct": do, "d_E": dE, "k": k, "Lf_rates": rates, "agree": do == dE + k}
        print(name, out[name], flush=True)
    (ROOT / "debug/data/hybrid_ratio_review_toy.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    run()
