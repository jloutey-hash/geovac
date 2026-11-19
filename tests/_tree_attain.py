"""Helpers for tests/test_paper63_tree_attainment.py (Paper 63, s-type Coulomb trees reach rung n, 3 <= n <= 8).

Hat normalisation: rho = K0, (Phi W)(x) = (1/x) int_0^x W s^2 ds + int_x^oo W s ds, theta(theta+1) Phi W = -x^2 W.
Root the tree at the vertex carrying t.  W_u = K0 prod_{w child of u} P_w, P_u = Phi W_u, k_u = |subtree(u)|,
H = x^2 prod_{w child of root} P_w, fibre = (2/pi)^n t^3 T(H), T(g) = int_0^oo K0(tx) g dx.

Module model (Lemma 1): A = C[z]<theta>, z = x^2.  Coordinates (S, kept, i): S a set of opened non-root vertices
(S + root connected), kept the frontier potentials each in state 0 (P_y) or 1 (theta P_y), i the power of
D = theta K0 in K0^(|S|-i) D^i.  Values: polynomials in z as coefficient lists (Python ints, or ints mod a prime).
  theta z = 2z; theta K0 = D; theta D = z K0; theta P_y = (theta P_y);
  theta (theta P_y) = -(theta P_y) - z K0 prod_{c child of y} P_c.
"Pure" = S = all non-root vertices (n - 1 Bessel factors, no potential left).

Nothing here imports from debug/.  The operators a in tests/data/paper63_tree_attain_2026_11_18.json were produced
once by exact (or mod-prime) linear algebra in this same model; the tests re-verify them from scratch.
"""
from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np
from scipy.special import k0e

Key = Tuple[frozenset, tuple, int]


class Tree:
    def __init__(self, n: int, edges: List[Tuple[int, int]], root: int = 0) -> None:
        self.n, self.root = n, root
        adj: Dict[int, List[int]] = {v: [] for v in range(n)}
        for a, b in edges:
            adj[a].append(b)
            adj[b].append(a)
        self.children: Dict[int, List[int]] = {}
        self.size: Dict[int, int] = {}

        def rec(v: int, p: int) -> int:
            self.children[v] = sorted(w for w in adj[v] if w != p)
            self.size[v] = 1 + sum(rec(w, v) for w in self.children[v])
            return self.size[v]
        assert rec(root, -1) == n, "edges do not form a tree on range(n)"

    def non_root(self) -> List[int]:
        return [u for u in range(self.n) if u != self.root]


# ------------------------------------------------------------------------------------------------ module model
def _padd(d: Dict, key: Key, poly: List[int], mod: int | None) -> None:
    cur = d.get(key)
    if cur is None:
        cur = d[key] = []
    if len(cur) < len(poly):
        cur.extend([0] * (len(poly) - len(cur)))
    for i, c in enumerate(poly):
        cur[i] = (cur[i] + c) % mod if mod else cur[i] + c


def _trim(d: Dict) -> Dict:
    out = {}
    for k, v in d.items():
        while v and v[-1] == 0:
            v.pop()
        if v:
            out[k] = v
    return out


def theta_vec(tree: Tree, vec: Dict, mod: int | None = None) -> Dict:
    out: Dict = {}
    for (S, kept, i), c in vec.items():
        k = len(S)
        _padd(out, (S, kept, i), [2 * m * cm for m, cm in enumerate(c)], mod)
        if k - i > 0:
            _padd(out, (S, kept, i + 1), [(k - i) * cm for cm in c], mod)
        if i > 0:
            _padd(out, (S, kept, i - 1), [0] + [i * cm for cm in c], mod)
        for idx, (yv, st) in enumerate(kept):
            rest = kept[:idx] + kept[idx + 1:]
            if st == 0:
                _padd(out, (S, tuple(sorted(rest + ((yv, 1),))), i), c, mod)
            else:
                _padd(out, (S, kept, i), [-cm for cm in c], mod)
                newkept = tuple(sorted(rest + tuple((w, 0) for w in tree.children[yv])))
                _padd(out, (S | frozenset([yv]), newkept, i), [0] + [-cm for cm in c], mod)
    return _trim(out)


def H_vector(tree: Tree) -> Dict:
    return {(frozenset(), tuple(sorted((w, 0) for w in tree.children[tree.root])), 0): [0, 1]}


def apply_operator(tree: Tree, a: Dict[Tuple[int, int], int], mod: int | None = None) -> Dict:
    """a H in the module model; a = {(k, j): c} = sum c z^k theta^j (z to the left)."""
    N = max(j for (_, j) in a)
    vs = [H_vector(tree)]
    for _ in range(N):
        vs.append(theta_vec(tree, vs[-1], mod))
    tot: Dict = {}
    for (k, j), c in a.items():
        for key, poly in vs[j].items():
            _padd(tot, key, [0] * k + [c * cm for cm in poly], mod)
    return _trim(tot)


def purity_report(tree: Tree, aH: Dict) -> Tuple[int, int]:
    """(number of non-pure coordinates left, number of pure coordinates nonzero)."""
    full = tree.n - 1
    nonpure = sum(1 for (S, _, _) in aH if len(S) != full)
    pure = sum(1 for (S, _, _) in aH if len(S) == full)
    return nonpure, pure


def weight(a: Dict[Tuple[int, int], int]) -> int:
    return max(2 * k + j for (k, j) in a)


def symbol_at_1_minus_t(a: Dict[Tuple[int, int], int], mod: int | None = None) -> List[int]:
    """Coefficients (ascending in t) of sigma(a)(1, -t) = sum_{2k+j=w} c (-t)^j: by (*) the leading coefficient of
    phi(a) is t^-w times this."""
    w = weight(a)
    N = max(j for (_, j) in a)
    lead = [0] * (N + 1)
    for (k, j), c in a.items():
        if 2 * k + j == w:
            lead[j] += c * (-1) ** j
    if mod:
        lead = [c % mod for c in lead]
    return lead


def a_mul(P: Dict, Q: Dict) -> Dict:
    """Product in A in normal form z^k theta^j, theta z = z (theta + 2); integer coefficients."""
    out: Dict = {}
    for (a, b), c1 in P.items():
        for (c, d), c2 in Q.items():
            # theta^b z^c = z^c (theta + 2c)^b
            for e in range(b + 1):
                cf = math.comb(b, e) * (2 * c) ** (b - e)
                if cf:
                    key = (a + c, e + d)
                    out[key] = out.get(key, 0) + c1 * c2 * cf
    return {k: v for k, v in out.items() if v != 0}


# ------------------------------------------------------------------------------------------------ phi
def phi_z_power(k: int, th):
    """phi(z^k) = (t^-2 theta_t^2)^k = t^(-2k) prod_{i<k} (theta_t - 2i)^2; this returns the theta-polynomial."""
    out = 1
    for i in range(k):
        out = out * (th - 2 * i) ** 2
    return out


# ------------------------------------------------------------------------------------------------ the split (Lemma 2)
def top_constant(tree: Tree):
    """C = (-1)^(n-1) (pi/2)^((n-1)/2) / prod_{u != root} k_u^2, as a sympy number."""
    import sympy as sp
    n = tree.n
    return (-1) ** (n - 1) * (sp.pi / 2) ** sp.Rational(n - 1, 2) / prod_k2(tree)


def prod_k2(tree: Tree) -> int:
    out = 1
    for u in tree.non_root():
        out *= tree.size[u] ** 2
    return out


def admissible_sets(tree: Tree) -> List[frozenset]:
    """All S of non-root vertices with S + root connected (the index set of the split H = sum q_S H^S)."""
    out = [frozenset()]
    frontier = [frozenset()]
    seen = {frozenset()}
    while frontier:
        nxt = []
        for S in frontier:
            cands = [c for v in (set(S) | {tree.root}) for c in tree.children[v] if c not in S]
            for c in cands:
                T = S | {c}
                if T not in seen:
                    seen.add(T)
                    out.append(T)
                    nxt.append(T)
        frontier = nxt
    return out


def claimed_rate(tree: Tree, S: frozenset) -> int:
    """Exponential rate of H^S claimed by Lemma 2: |S| (one e^-x per opened Bessel factor)."""
    return len(S)


GL_X, GL_W = np.polynomial.laguerre.laggauss(40)


def scaled_HS(tree: Tree, S: frozenset, x: float) -> float:
    """e^{|S| x} H^S(x) by nested Gauss-Laguerre on the exterior-tail formula
    E(x) = -(1/x) int_0^oo W(x + v) (x + v) v dv; children outside S contribute x^-1 (their Q is in q_S).
    Independent of the local rule and of the constant C."""
    def kS(u: int) -> int:
        return (1 if u in S else 0) + sum(kS(c) for c in tree.children[u])

    def wsc(u: int, s: np.ndarray) -> np.ndarray:          # e^{k_u^S s} W_u^(S)(s)
        out = k0e(s)
        for c in tree.children[u]:
            out = out * (esc(c, s) if c in S else 1.0 / s)
        return out

    def esc(u: int, xx: np.ndarray) -> np.ndarray:         # e^{k_u^S x} E_u^(S)(x)
        k = kS(u)
        v = GL_X / k
        s = xx[..., None] + v
        return -(1.0 / xx) * np.sum(GL_W / k * wsc(u, s) * s * v, axis=-1)

    xx = np.array(float(x))
    h = xx ** 2
    for w in tree.children[tree.root]:
        h = h * (esc(w, xx) if w in S else 1.0 / xx)
    return float(h)


def top_constant_numeric(tree: Tree, xs=(100.0, 200.0, 400.0, 800.0, 1600.0)) -> float:
    """C from e^{(n-1)x} H^top(x) x^{-(5-n)/2} -> C (1 + c1/x + ...), Richardson in 1/x with doubling."""
    full = frozenset(tree.non_root())
    a = (5 - tree.n) / 2
    tab = [scaled_HS(tree, full, x) * x ** (-a) for x in xs]
    lvl = 1
    while len(tab) > 1:
        tab = [(2 ** lvl * tab[i + 1] - tab[i]) / (2 ** lvl - 1) for i in range(len(tab) - 1)]
        lvl += 1
    return tab[0]


def nonholomorphy_derivative_order(n: int) -> int:
    """Least m with d^m J/d eps^m unbounded as eps -> 0+: the integrand of J is ~ x^((4-n)/2) e^(-eps x), and each
    eps-derivative multiplies by -x, so m is least with (4-n)/2 + m >= -1.  0 for n <= 6, 1 for n = 7, 8."""
    return max(0, math.ceil((n - 6) / 2))
