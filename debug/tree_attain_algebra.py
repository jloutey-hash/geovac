"""tree_attain: exact module model of an s-type Coulomb tree (hat normalisation rho = K0, Phi W = int W/max).

For a rooted tree, H = z prod_{w child of root} P_w  (z = x^2), fibre = (2/pi)^n t^3 T(H), T(g) = int_0^oo K0(tx) g dx.
Function module: C[z]-span of monomials  K0^(k-i) D^i  x  prod_{y in frontier(S)} (P_y or theta P_y),
D = theta K0, S a set of non-root vertices with S + root connected, k = |S|.  Rules (theta = x d/dx):
  theta z = 2 z;  theta K0 = D;  theta D = z K0;  theta P_y = (theta P_y);
  theta (theta P_y) = -(theta P_y) - z K0 prod_{c child of y} P_c      (radial Poisson equation).
"Pure top" = S = all non-root vertices (n-1 Bessel factors, no potentials left).

Computed, exactly:
  (1) chi_c: the top-weight part of L_c (annihilator of K0^c), weights theta = 1, z = 2, is prod_j (theta - (c-2j) x);
  (2) for each tree: the least-order a = sum c_j(z) theta^j with a H in the pure top; its weighted symbol at
      (x = 1, theta = -t), which is the leading coefficient of phi(a) up to a power of t; the pure element g = a H;
  (3) numerically (float grid): phi(a) T(H) = T(g) and T(g) != 0, phi applied under the integral sign.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp
from scipy.integrate import cumulative_simpson
from scipy.special import k0, k1
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
import test_paper63_tree_rung as TR  # noqa: E402

z, th, tt, yy = sp.symbols("z theta t y")
T0 = time.time()


def log(*a) -> None:
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


# ------------------------------------------------------------------------------------------------ (1) symbols
def top_weight_symbol(c: int) -> sp.Expr:
    L = TR._sym_op(c)
    w = max(2 * k + j for (k, j) in L)
    X = sp.symbols("X")
    return w, sp.expand(sum(cf * X ** (2 * k) * th ** j for (k, j), cf in L.items() if 2 * k + j == w))


def check_symbols(cmax: int) -> Dict:
    X = sp.symbols("X")
    out = {}
    for c in range(1, cmax + 1):
        w, s = top_weight_symbol(c)
        target = sp.expand(sp.prod([th - (c - 2 * j) * X for j in range(c + 1)]))
        ok = (w == c + 1) and sp.expand(s - target) == 0
        out[c] = bool(ok)
        log(f"chi_{c}: weight {w}, top-weight part = prod (theta - (c-2j) x): {ok}")
    return out


# ------------------------------------------------------------------------------------------------ (2) module
class Tree:
    def __init__(self, n: int, edges: List[Tuple[int, int]], root: int = 0) -> None:
        self.n, self.root = n, root
        adj = {v: [] for v in range(n)}
        for a, b in edges:
            adj[a].append(b)
            adj[b].append(a)
        self.children: Dict[int, List[int]] = {}
        self.size: Dict[int, int] = {}

        def rec(v: int, p: int) -> int:
            self.children[v] = [w for w in adj[v] if w != p]
            self.size[v] = 1 + sum(rec(w, v) for w in self.children[v])
            return self.size[v]
        rec(root, -1)


def theta_vec(tree: Tree, vec: Dict) -> Dict:
    out: Dict = {}

    def add(key, p) -> None:
        q = out.get(key, 0) + p
        out[key] = sp.expand(q)

    for (S, kept, i), c in vec.items():
        k = len(S)
        dc = sp.expand(2 * z * sp.diff(c, z))
        if dc != 0:
            add((S, kept, i), dc)
        if k - i > 0:
            add((S, kept, i + 1), (k - i) * c)
        if i > 0:
            add((S, kept, i - 1), i * z * c)
        for idx, (yv, st) in enumerate(kept):
            rest = kept[:idx] + kept[idx + 1:]
            if st == 0:
                add((S, tuple(sorted(rest + ((yv, 1),))), i), c)
            else:
                add((S, kept, i), -c)
                newkept = tuple(sorted(rest + tuple((w, 0) for w in tree.children[yv])))
                add((S | frozenset([yv]), newkept, i), -z * c)
    return {k_: v for k_, v in out.items() if v != 0}


def analyse(tree: Tree, label: str, max_order: int = 60) -> Dict:
    n = tree.n
    v = {(frozenset(), tuple(sorted((w, 0) for w in tree.children[tree.root])), 0): z}
    vs = [v]
    is_pure = lambda key: len(key[0]) == n - 1  # noqa: E731
    z0 = sp.Rational(5, 7)
    N = None
    for j in range(1, max_order + 1):
        vs.append(theta_vec(tree, vs[-1]))
        cols = sorted({k_ for u in vs for k_ in u if not is_pure(k_)}, key=str)
        M = sp.Matrix([[sp.sympify(u.get(k_, 0)).subs(z, z0) for k_ in cols] for u in vs])
        if M.rank() < len(vs):
            N = j
            break
    assert N is not None, "no relation found"
    cols = sorted({k_ for u in vs for k_ in u if not is_pure(k_)}, key=str)
    log(f"{label}: least order of a = {N}; non-pure rank used = {len(cols)}")
    M = sp.Matrix([[u.get(k_, 0) for k_ in cols] for u in vs])            # (N+1) x cols over Q[z]
    dM = DomainMatrix.from_Matrix(M.T).to_field()
    ns = dM.nullspace().to_Matrix()
    assert ns.shape[0] == 1, ns.shape
    cvec = [sp.together(e) for e in ns.row(0)]
    den = sp.lcm([sp.fraction(e)[1] for e in cvec])
    cvec = [sp.expand(sp.cancel(e * den)) for e in cvec]
    cont = sp.gcd_list(cvec)
    cvec = [sp.expand(sp.cancel(e / cont)) for e in cvec]
    if sp.Poly(cvec[-1], z).LC() < 0:
        cvec = [-e for e in cvec]
    # residual check: a H is pure
    tot: Dict = {}
    for cj, u in zip(cvec, vs):
        for k_, p in u.items():
            tot[k_] = sp.expand(tot.get(k_, 0) + cj * p)
    assert all(p == 0 for k_, p in tot.items() if not is_pure(k_)), "a H not pure"
    g = {k_[2]: p for k_, p in tot.items() if is_pure(k_) and p != 0}
    # operator a as {(k, j): coeff}, weight and symbol
    a: Dict = {}
    for j, cj in enumerate(cvec):
        for (k,), cf in sp.Poly(cj, z).terms():
            a[(k, j)] = cf
    w = max(2 * k + j for (k, j) in a)
    lead = sp.expand(sum(cf * (-tt) ** j for (k, j), cf in a.items() if 2 * k + j == w))
    lead_f = sp.factor(lead)
    roots = sp.roots(sp.Poly(lead, tt))
    p = -(n - 1)
    res = {
        "label": label, "n": n, "order_theta": N, "weight": w, "nonpure_rank": len(cols),
        "leading_coefficient_of_phi_a_times_t^w": str(lead_f),
        "roots": {str(r): int(m) for r, m in roots.items()},
        "all_roots_rational_found": int(sum(roots.values())) == sp.Poly(lead, tt).degree(),
        "max_abs_root": str(max([abs(r) for r in roots], default=0)),
        "lead_at_p": str(lead.subs(tt, p)), "nonzero_at_p": bool(lead.subs(tt, p) != 0),
        "roots_within_n_minus_2": bool(all(abs(r) <= n - 2 for r in roots)),
        "g_nonzero": bool(len(g) > 0),
        "g_degrees": {str(i): int(sp.Poly(pp, z).degree()) for i, pp in g.items()},
    }
    log(f"{label}: weight {w}; t^w lead(phi a) = {lead_f}")
    log(f"{label}: nonzero at t={p}: {res['nonzero_at_p']}; roots within |t|<=n-2: {res['roots_within_n_minus_2']}; "
        f"g != 0: {res['g_nonzero']}")
    return res, a, g


# ------------------------------------------------------------------------------------------------ (3) numerics
U = np.linspace(np.log(1e-10), np.log(150.0), 600001)
R = np.exp(U)


def cum0(f):
    return cumulative_simpson(f * R, x=U, initial=0.0)


def coulomb_hat(W):
    c1 = cum0(W * R ** 2)
    c2 = cum0(W * R)
    return c1 / R + (c2[-1] - c2)


def H_hat(tree: Tree):
    def W(v: int):
        out = k0(R)
        for w_ in tree.children[v]:
            out = out * coulomb_hat(W(w_))
        return out
    out = R ** 2
    for w_ in tree.children[tree.root]:
        out = out * coulomb_hat(W(w_))
    return out


_KCACHE: Dict = {}


def kernel(k: int, j: int):
    """(y^-2 theta^2)^k (-1 - theta)^j K0(y) = al(y) K0 + be(y) K1."""
    if (k, j) in _KCACHE:
        return _KCACHE[(k, j)]
    al, be = sp.Integer(1), sp.Integer(0)

    def th_(al, be):
        return (sp.expand(yy * sp.diff(al, yy) - yy * be), sp.expand(yy * sp.diff(be, yy) - yy * al - be))
    for _ in range(j):
        a1, b1 = th_(al, be)
        al, be = sp.expand(-al - a1), sp.expand(-be - b1)
    for _ in range(k):
        a1, b1 = th_(*th_(al, be))
        al, be = sp.expand(a1 / yy ** 2), sp.expand(b1 / yy ** 2)
    f = (sp.lambdify(yy, al, "numpy"), sp.lambdify(yy, be, "numpy"))
    _KCACHE[(k, j)] = f
    return f


def numeric_identity(tree: Tree, a: Dict, g: Dict, tv: float) -> Dict:
    n = tree.n
    H = H_hat(tree)
    Y = tv * R
    K0y, K1y = k0(Y), k1(Y)
    terms = []
    for (k, j), cf in a.items():
        fa, fb = kernel(k, j)
        ker = fa(Y) * K0y + fb(Y) * K1y
        terms.append(float(cf) * np.trapezoid(ker * R ** (2 * k) * H * R, U))
    lhs = float(np.sum(terms))
    scale = float(np.max(np.abs(terms)))
    D = -R * k1(R)
    gfun = np.zeros_like(R)
    for i, p in g.items():
        gfun += sp.lambdify(z, p, "numpy")(R ** 2) * k0(R) ** (n - 1 - i) * D ** i
    rhs = float(np.trapezoid(K0y * gfun * R, U))
    return {"t": tv, "lhs_phi_a_T_H": lhs, "rhs_T_g": rhs, "abs_diff": abs(lhs - rhs),
            "largest_term": scale, "rel_to_rhs": abs(lhs - rhs) / abs(rhs)}


TREES = {
    "star3_root_centre": (3, [(0, 1), (0, 2)]),
    "path3_root_end": (3, [(0, 1), (1, 2)]),
    "star4": (4, [(0, 1), (0, 2), (0, 3)]),
    "path4_root_end": (4, [(0, 1), (1, 2), (2, 3)]),
    "path4_root_inner": (4, [(0, 1), (0, 2), (2, 3)]),
    "path5_root_end": (5, [(0, 1), (1, 2), (2, 3), (3, 4)]),
    "fork5_root_end": (5, [(0, 1), (1, 2), (1, 3), (3, 4)]),
    "path5_root_centre": (5, [(0, 1), (1, 2), (0, 3), (3, 4)]),
}

if __name__ == "__main__":
    which = sys.argv[1:] or list(TREES)
    outp = ROOT / "debug" / "data" / "tree_attain_algebra.json"
    data = json.loads(outp.read_text()) if outp.exists() else {}
    if "symbols" not in data:
        data["symbols"] = {str(k): v for k, v in check_symbols(7).items()}
    for name in which:
        n, edges = TREES[name]
        tree = Tree(n, edges)
        t1 = time.time()
        res, a, g = analyse(tree, name)
        res["seconds_algebra"] = round(time.time() - t1, 1)
        res["numeric"] = []
        for tv in ((1.5, 2.5) if res["weight"] <= 20 else ()):   # float64 cancels beyond this (star4: weight 56)
            r = numeric_identity(tree, a, g, tv)
            res["numeric"].append(r)
            log(f"{name}: t={tv}: phi(a)T(H) = {r['lhs_phi_a_T_H']:.10e}, T(g) = {r['rhs_T_g']:.10e}, "
                f"diff/|T(g)| = {r['rel_to_rhs']:.2e}, largest term {r['largest_term']:.2e}")
        data[name] = res
        outp.write_text(json.dumps(data, indent=1))
    log("done")
