"""Paper 63, s-type part of conj:tree_attain (and the replacement for the last step of prop:star_attain's proof):
for every Coulomb tree of n s-type densities at the fingerprint weights, 3 <= n <= 8, slice (t; 1, ..., 1) with t
on any one vertex (the root), F_n lies in the D-module of the fibre, so the minimal operator has a composition
factor equivalent to F_n's operator.  Proof record: Lemma 1 (an operator a in A = C[z]<theta>, a product of pieces
with symbols chi_c^m, c <= n-2, with a H a pure (n-1)-Bessel polynomial) and Lemma 2 (H = sum_S q_S H^S, the lower
pieces holomorphic at t = -(n-1), one top term with constant
C = (-1)^(n-1) (pi/2)^((n-1)/2) / prod_{u != root} k_u^2, k_u = |subtree(u)|, and an unbounded jump).
Written as a separate pass from the proof; derived 2026-11-18 and cold-reviewed the same day.

Each test names the wrong answer it rejects.  Covered:
  1. the old inference "all composition factors regular at p => holomorphic at p" (toys log(t-p), 1/(t-p));
  2. a symbol of a with a root of absolute value n-1 (e.g. L_(n-1) a): roots of sigma(a)(1,-t) are integers with
     |t| <= n-2 and the symbol is nonzero at -(n-1) -- exact for star 3, path 3, star 4, path 4 (end and inner),
     path 5 (end); modulo the prime 2147483629 for the fork of five rooted at a leaf of its long arm and the path
     of six from an end (a branched tree, and n = 6);
  3. a H not pure (one non-pure coordinate left), and a H = 0;
  4. phi(z^k) = t^(-2k) theta^(2k) instead of (t^-2 theta^2)^k: phi(a) T(H) = T(a H) for the star of three to
     1e-9 in float64 (fast) and to 1e-18 at 30 digits (slow), and the star of four to 1e-15 at 45 digits (slow);
     plus (*) (order = weight, leading coefficient t^-w sigma(1,-t)) exactly for the star and path of three;
  5. the top constant with k_u = 1 everywhere, with prod k_u, or with sign (-1)^n: C against an independent nested
     Gauss-Laguerre evaluation of e^{(n-1)x} H^top (fork of five at two roots, spider of six with prod k_u^2 = 225,
     path of five from the centre, a branched seven and eight) and against the exact local-rule series lead;
  6. a lower piece credited with rate n-1: every S except the full set has |S| <= n-2, and the measured decay rate
     of each H^S equals |S| (several trees, every admissible S);
  7. for n = 7, 8 a check on J instead of its eps-derivative: J is bounded there, d J/d eps is unbounded.

Not asserted (the proof does not claim it): l > 0, minimality or the full factor list, unequal non-root scales,
anything for n >= 9 beyond "meets", irreducibility of F_6..F_8's operators (inherited, guarded by
tests/test_paper63_four_density_rung.py), the transfer's boundary terms (tests/test_paper63_tree_rung.py).

Fast tier does NOT cover: phi(a) T(H) = T(aH) at n = 4 (slow: star of four, ~10-15 digits of cancellation) and
the 30-digit star-of-three check (slow; the fast tier has its float64 twin at 1e-9);
the exact (non-modular) symbol of a branched tree (only mod one prime -- a wrong reduction would need the least
order to drop mod that prime; no second prime); no stored operator at n = 7, 8 and none for the fork rooted at its
branch vertex (Lemma 1 rests on its proof there); item 7 uses the star of seven and eight (true H^top through R),
not a branched tree -- the branched 7- and 8-trees enter only through C (item 5).

Data: tests/data/paper63_tree_attain_2026_11_18.json (the least-order operators a, produced once by exact / mod-p
linear algebra in the module model of tests/_tree_attain.py; every property is re-verified here).  Nothing is
imported from debug/.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, Tuple

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy.special import i0e, i1e, k1e, roots_laguerre

import _tree_attain as TA
import test_paper63_tree_rung as TR

DATA = Path(__file__).resolve().parent / "data" / "paper63_tree_attain_2026_11_18.json"
t, th = sp.symbols("t theta")


def _load() -> Dict:
    raw = json.loads(DATA.read_text())
    out = {}
    for name, rec in raw.items():
        if name.startswith("_"):
            continue
        a = {tuple(int(x) for x in key.split(",")): int(v) for key, v in rec["a"].items() if int(v) != 0}
        out[name] = (TA.Tree(rec["n"], [tuple(e) for e in rec["edges"]], rec["root"]), a, rec.get("prime"))
    return out


OPS = _load()


def _ops(name: str) -> Tuple[TA.Tree, Dict, int]:
    return OPS[name]


EXACT = ["star3_root_centre", "path3_root_end", "star4", "path4_root_end", "path4_root_inner", "path5_root_end"]
MODP = ["fork5_root4", "path6_end"]


# ================================================================================== 1. the old inference
def test_old_inference_has_counterexamples_the_new_criterion_rejects():
    """Rejects "all composition factors regular at p implies holomorphic at p" (the written last step of
    prop:star_attain).  Toy 1: d/dt o (t-p) d/dt, whose two factors d/dt and (t-p) d/dt both have only the
    constant solution (each equivalent to d/dt over C(t)), kills log(t-p), which is not holomorphic at p.
    Toy 2: 1/(t-p) lies in C(t).1, whose only piece is the trivial one, and has a pole at p.  The replacement
    argument asks instead that the LEADING COEFFICIENT of the operator not vanish at p; both toys fail that test
    (their leading coefficient is t - p), and phi(a) passes it (item 2)."""
    p = sp.symbols("p")
    y = sp.Function("y")
    # factors: solutions of d/dt y = 0 and (t-p) y' = 0 are the constants only (holomorphic at p)
    for fac in (sp.Derivative(y(t), t), (t - p) * sp.Derivative(y(t), t)):
        sol = sp.dsolve(sp.Eq(fac, 0), y(t))
        assert sol.rhs.free_symbols <= {sp.Symbol("C1")}
    # toy 1: the composite kills log(t - p), which is not holomorphic at p (its derivative has residue 1)
    f = sp.log(t - p)
    assert sp.simplify(sp.diff((t - p) * sp.diff(f, t), t)) == 0
    assert sp.residue(sp.diff(f, t), t, p) == 1
    # toy 2: g = (1/(t-p)) . 1 is in the module of 1 (annihilated by (t-p) d/dt + 1, gauge-equivalent to d/dt)
    g = 1 / (t - p)
    assert sp.simplify((t - p) * sp.diff(g, t) + g) == 0
    assert sp.residue(g, t, p) == 1
    # the criterion the new proof uses: leading coefficient nonzero at p.  The toy composite's is zero there.
    comp = sp.expand(sp.diff((t - p) * sp.diff(y(t), t), t))
    lead = comp.coeff(sp.Derivative(y(t), (t, 2)))
    assert sp.expand(lead - (t - p)) == 0 and lead.subs(t, p) == 0


# ================================================================================== 2/3. Lemma 1 on stored operators
def _int_roots_exact(coeffs_asc, n: int):
    P = sp.Poly(list(reversed(coeffs_asc)), t)
    _, facs = sp.factor_list(P)
    roots = {}
    for f, m in facs:
        assert f.degree() == 1, f"non-linear factor {f.as_expr()}"
        r = -f.coeff_monomial(1) / f.coeff_monomial(t)
        assert r.is_integer, r
        roots[int(r)] = roots.get(int(r), 0) + m
    return roots, P.eval(-(n - 1))


def _int_roots_modp(coeffs_asc, n: int, P: int):
    poly = [c % P for c in coeffs_asc]
    while poly and poly[-1] == 0:
        poly.pop()
    val = sum(c * pow(-(n - 1), i, P) for i, c in enumerate(poly)) % P
    roots = {}
    for r in range(-(n + 2), n + 3):
        while len(poly) > 1:
            # synthetic division by (t - r), ascending coefficients
            q = [0] * (len(poly) - 1)
            acc = 0
            for i in range(len(poly) - 1, 0, -1):
                acc = (acc * r + poly[i]) % P
                q[i - 1] = acc
            rem = (acc * r + poly[0]) % P
            if rem != 0:
                break
            poly = q
            roots[r] = roots.get(r, 0) + 1
    return roots, val, len(poly) - 1


@pytest.mark.parametrize("name", EXACT)
def test_symbol_roots_exact(name):
    """Rejects an operator whose symbol at (1, -t) has a root of absolute value n-1 (e.g. a left-multiplied by
    L_(n-1)): then phi(a)'s leading coefficient vanishes at p = -(n-1) and the contradiction in (S) is lost.
    Exact over Q: the symbol splits into integer roots, all |r| <= n-2, and is nonzero at -(n-1)."""
    tree, a, _ = _ops(name)
    n = tree.n
    roots, at_p = _int_roots_exact(TA.symbol_at_1_minus_t(a), n)
    assert roots and max(abs(r) for r in roots) <= n - 2, roots
    assert at_p != 0


@pytest.mark.parametrize("name", MODP)
def test_symbol_roots_modp(name):
    """Same wrong answer as test_symbol_roots_exact, modulo the prime: the symbol splits completely into roots
    among 0, +-1, ..., +-(n-2) and is nonzero at -(n-1) mod p."""
    tree, a, P = _ops(name)
    n = tree.n
    roots, val, deg_left = _int_roots_modp(TA.symbol_at_1_minus_t(a, P), n, P)
    assert deg_left == 0, (roots, deg_left)
    assert max(abs(r) for r in roots) <= n - 2, roots
    assert val != 0


@pytest.mark.parametrize("name", EXACT + MODP)
def test_aH_is_pure_and_nonzero(name):
    """Rejects a H not pure (one non-pure coordinate left, e.g. a with one term dropped) and a H = 0.  a H is
    recomputed in the module model from scratch (theta^j H, j <= order of a), exactly or mod p."""
    tree, a, P = _ops(name)
    aH = TA.apply_operator(tree, a, P)
    nonpure, pure = TA.purity_report(tree, aH)
    assert nonpure == 0
    assert pure >= 1


@pytest.mark.parametrize("c", [1, 2, 3, 4, 5])
def test_top_weight_part_of_L_c_is_chi_c(c):
    """The symbol chi_c = prod_j (Theta - (c-2j) X) of L_c (weights theta 1, z 2): the source of Lemma 1's root
    bound.  Rejects a weight above c+1 and a symbol with roots outside +-c, +-(c-2), ..."""
    L = TR._sym_op(c)
    w = max(2 * k + j for (k, j) in L)
    X = sp.symbols("X")
    sym = sp.expand(sum(cf * X ** (2 * k) * th ** j for (k, j), cf in L.items() if 2 * k + j == w))
    assert w == c + 1
    assert sp.expand(sym - sp.prod([th - (c - 2 * j) * X for j in range(c + 1)])) == 0


@pytest.mark.parametrize("name", ["star3_root_centre", "path3_root_end"])
def test_star_rule_order_equals_weight(name):
    """(*): phi(a) has order weight(a) in theta_t with leading coefficient t^-w sigma(a)(1,-t), exactly.  Rejects
    reading the leading coefficient off the theta-order of a instead of its weight."""
    _, a, _ = _ops(name)
    w = TA.weight(a)
    Pa = TR._phi({k: sp.Integer(v) for k, v in a.items()})        # {(k, j)}: u^k theta_t^j, u = t^-2
    order = max(j for (_, j) in Pa)
    lead = sum(c * t ** (-2 * k) for (k, j), c in Pa.items() if j == order)
    sig = sum(c * tt for c, tt in zip(TA.symbol_at_1_minus_t(a), [t ** i for i in range(w + 1)]))
    assert order == w
    assert sp.simplify(lead - t ** (-w) * sig) == 0


# ================================================================================== 4. the transfer identity
def _theta_K_tables(mmax: int):
    y = sp.symbols("y")
    A, B = [sp.Integer(1)], [sp.Integer(0)]
    for _ in range(mmax):
        al, be = A[-1], B[-1]
        A.append(sp.expand(y * sp.diff(al, y) - y * be))           # theta K0 = -y K1, theta K1 = -y K0 - K1
        B.append(sp.expand(y * sp.diff(be, y) - y * al - be))
    return y, A, B


def _chi(x):
    """chi(x) = int_0^x s K1(s) ds (closed form; chi = x P_leaf in the hat normalisation)."""
    return -x * mp.besselk(0, x) + (mp.pi * x / 2) * (mp.besselk(0, x) * mp.struvel(-1, x)
                                                      + mp.besselk(1, x) * mp.struvel(0, x))


def _transfer_sides(name: str, tq: sp.Rational, dps: int, xmax: int = 128, jmin: int = 16, ng: int = 24):
    """(relative difference, |T(aH)|) computed at dps digits; the global mpmath precision is restored."""
    with mp.workdps(dps):
        lhs, rhs = _transfer_sides_raw(name, tq, xmax, jmin, ng)
        return float(abs(lhs - rhs) / abs(rhs)), float(abs(rhs))


def _transfer_sides_raw(name: str, tq: sp.Rational, xmax: int, jmin: int, ng: int):
    """phi(a) T(H) with phi applied under the integral sign, and T(a H) with a H from the module model.  Star
    trees only (H = chi^(n-1) x^(3-n))."""
    tree, a, _ = _ops(name)
    n = tree.n
    pol = sp.Poly(sp.expand(sum(sp.Integer(c) * tq ** (-2 * k) * TA.phi_z_power(k, th) * (-1 - th) ** j
                                for (k, j), c in a.items())), th)
    coeffs = pol.all_coeffs()[::-1]
    y, A, B = _theta_K_tables(len(coeffs) - 1)

    def to_mp(c):
        c = sp.Rational(c)
        return mp.mpf(int(c.p)) / mp.mpf(int(c.q))
    alc = [to_mp(c) for c in sp.Poly(sp.expand(sum(c * A[m] for m, c in enumerate(coeffs))), y).all_coeffs()]
    bec = [to_mp(c) for c in sp.Poly(sp.expand(sum(c * B[m] for m, c in enumerate(coeffs))), y).all_coeffs()]
    aH = TA.apply_operator(tree, a)
    g = {key[2]: [mp.mpf(c) for c in poly] for key, poly in aH.items()}   # pure: key (S_full, (), i)
    assert all(len(key[0]) == n - 1 for key in aH)
    tv = mp.mpf(int(tq.p)) / int(tq.q)

    def lhs_f(x):
        yv = tv * x
        return (mp.polyval(alc, yv) * mp.besselk(0, yv) + mp.polyval(bec, yv) * mp.besselk(1, yv)) \
            * _chi(x) ** (n - 1) * x ** (3 - n)

    def rhs_f(x):
        k0v, D = mp.besselk(0, x), -x * mp.besselk(1, x)
        s = mp.fsum(mp.polyval(cs[::-1], x * x) * k0v ** (n - 1 - i) * D ** i for i, cs in g.items())
        return mp.besselk(0, tv * x) * s

    # fixed Gauss-Legendre panels: geometric below 1 (the integrands are O(x^2 log^m x) at 0), then up to 128
    edges = [mp.mpf(0)] + [mp.mpf(2) ** -j for j in range(jmin, 0, -1)] + [mp.mpf(v) for v in
                                                                           (1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24,
                                                                            28, 32, 40, 48, 56, 64, 80, 96, 112, 128) if v <= xmax]
    GX, GW = mp.gauss_quadrature(ng, "legendre")
    lhs = rhs = mp.mpf(0)
    for lo, hi in zip(edges[:-1], edges[1:]):
        h, c = (hi - lo) / 2, (hi + lo) / 2
        for xg, wg in zip(GX, GW):
            x = c + h * xg
            lhs += h * wg * lhs_f(x)
            rhs += h * wg * rhs_f(x)
    return lhs, rhs


def _transfer_sides_f64(name: str, tq: sp.Rational):
    """float64 twin of _transfer_sides for the fast tier (scipy Bessel/Struve, Gauss-Legendre panels)."""
    from scipy.special import k0 as K0, k1 as K1, modstruve
    tree, a, _ = _ops(name)
    n = tree.n
    pol = sp.Poly(sp.expand(sum(sp.Integer(c) * tq ** (-2 * k) * TA.phi_z_power(k, th) * (-1 - th) ** j
                                for (k, j), c in a.items())), th)
    coeffs = pol.all_coeffs()[::-1]
    y, A, B = _theta_K_tables(len(coeffs) - 1)
    alc = [float(c) for c in sp.Poly(sp.expand(sum(c * A[m] for m, c in enumerate(coeffs))), y).all_coeffs()]
    bec = [float(c) for c in sp.Poly(sp.expand(sum(c * B[m] for m, c in enumerate(coeffs))), y).all_coeffs()]
    aH = TA.apply_operator(tree, a)
    g = {key[2]: [float(c) for c in poly] for key, poly in aH.items()}
    edges = np.concatenate([[0.0], 2.0 ** -np.arange(12, 0, -1), np.arange(1.0, 65.0)])
    gx, gw = np.polynomial.legendre.leggauss(30)
    lo, hi = edges[:-1, None], edges[1:, None]
    x = ((hi + lo) / 2 + (hi - lo) / 2 * gx).ravel()
    w = ((hi - lo) / 2 * gw).ravel()
    tv = float(tq)
    chi = -x * K0(x) + (np.pi * x / 2) * (K0(x) * (modstruve(1, x) + 2 / np.pi) + K1(x) * modstruve(0, x))
    H = chi ** (n - 1) * x ** (3 - n)
    yv = tv * x
    lhs = np.sum(w * (np.polyval(alc, yv) * K0(yv) + np.polyval(bec, yv) * K1(yv)) * H)
    D = -x * K1(x)
    gx_ = sum(np.polyval(cs[::-1], x * x) * K0(x) ** (n - 1 - i) * D ** i for i, cs in g.items())
    rhs = np.sum(w * K0(yv) * gx_)
    return abs(lhs - rhs) / abs(rhs), abs(rhs)


def test_transfer_identity_star3_float64():
    """Fast twin of test_transfer_identity_star3 (float64, star of three, t = 3/2, agreement to 1e-9; the node
    terms reach 4e5 against 9.4e4).  Same wrong answers: phi(z^k) = t^(-2k) theta^(2k), and a sign slip in
    phi(theta) = -1 - theta (both are off by order one)."""
    rel, size = _transfer_sides_f64("star3_root_centre", sp.Rational(3, 2))
    assert size > 1e3
    assert rel < 1e-9


@pytest.mark.slow
def test_transfer_identity_star3():
    """phi(a) T(H) = T(a H) for the star of three at t = 3/2, 30-digit arithmetic, agreement to 1e-18
    (measured 2.7e-20), both sides nonzero.  Rejects
    phi(z^k) = t^(-2k) theta^(2k) in place of (t^-2 theta^2)^k = t^(-2k) prod_{i<k} (theta - 2i)^2 (the review's
    first attempt: 21607.7 against -93648.2), and any sign slip in phi(theta) = -1 - theta."""
    rel, size = _transfer_sides("star3_root_centre", sp.Rational(3, 2), 30, xmax=64, jmin=12, ng=20)
    assert size > 1e3
    assert rel < 1e-18


@pytest.mark.slow
def test_transfer_identity_star4():
    """Same wrong answer as test_transfer_identity_star3, at n = 4 (star of four, t = 5/2, 45 digits; the node
    terms reach 1e60 against a result 2.45e49, so float64 cannot see it)."""
    rel, size = _transfer_sides("star4", sp.Rational(5, 2), 45)
    assert size > 1e40
    assert rel < 1e-15


# ================================================================================== 5. the top constant
TREES_C = {
    "fork5_root_end": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 0, 64),
    "fork5_root_branch": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 1, 4),
    "spider6": (6, [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5)], 0, 225),
    "path5_centre": (5, [(0, 1), (1, 2), (0, 3), (3, 4)], 0, 16),
    "tree7_branched": (7, [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (5, 6)], 0, 2304),
    "tree8_branched": (8, [(0, 1), (0, 2), (2, 3), (2, 4), (4, 5), (4, 6), (6, 7)], 0, 2304),
}


@pytest.mark.parametrize("name", list(TREES_C))
def test_top_constant_against_nested_quadrature(name):
    """C = (-1)^(n-1) (pi/2)^((n-1)/2) / prod_{u != root} k_u^2 against e^{(n-1)x} H^top x^{-(5-n)/2} from nested
    Gauss-Laguerre (x = 100..1600, Richardson), to 1e-8.  Rejects k_u = 1 everywhere (off by 64 on the fork from
    its end, 225 on the spider), prod k_u instead of prod k_u^2, and the sign (-1)^n.  Also pins the spider's
    prod k_u^2 = 225 (the proof memo's table printed 64)."""
    n, edges, root, pk2 = TREES_C[name]
    tree = TA.Tree(n, edges, root)
    assert TA.prod_k2(tree) == pk2
    C = float(TA.top_constant(tree))
    num = TA.top_constant_numeric(tree)
    assert abs(num / C - 1) < 1e-8, (num, C)


@pytest.mark.parametrize("name", ["fork5_root_end", "fork5_root_branch", "spider6", "path5_centre"])
def test_top_constant_equals_local_rule_lead(name):
    """Second route, exact: the lead of the local-rule series for H^top (tests/test_paper63_local_rule.py Tail;
    each potential's exponential part from -(s_0/k^2)) equals C symbolically, and the series has rate n-1 and
    power (5-n)/2.  Same wrong answers as the nested-quadrature test."""
    import test_paper63_local_rule as LR
    n, edges, root, _ = TREES_C[name]
    tree = TA.Tree(n, edges, root)

    def full_series(v):
        out = LR.K_tail(0)
        for w in tree.children[v]:
            out = out * full_series(w).potential(0).times(1 / (4 * sp.pi))
        return out
    Hs = None
    for w in tree.children[tree.root]:
        e = full_series(w).potential(0).times(1 / (4 * sp.pi))
        Hs = e if Hs is None else Hs * e
    Hs = Hs.times(1, shift=2)
    assert Hs.k == n - 1 and Hs.a == sp.Rational(5 - n, 2)
    assert sp.simplify(Hs.lead() - TA.top_constant(tree)) == 0


# ================================================================================== 6. lower pieces
RATE_TREES = {
    "path4_inner": (4, [(0, 1), (0, 2), (2, 3)], 0),
    "fork5_root_end": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 0),
    "fork5_root_branch": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 1),
    "spider6": (6, [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5)], 0),
}


@pytest.mark.parametrize("name", list(RATE_TREES))
def test_lower_pieces_have_rate_at_most_n_minus_2(name):
    """Rejects a lower piece credited with rate n-1 (it would be singular at p and the top term would not be the
    only one): every admissible S other than the full set has |S| <= n-2, and the decay rate of H^S, measured
    from the nested-quadrature values at x = 200 and 400 (|S| minus the log-slope of e^{|S|x} H^S), equals the
    claimed rate |S| to 0.05, for every admissible S (the full set included: rate n-1)."""
    n, edges, root = RATE_TREES[name]
    tree = TA.Tree(n, edges, root)
    full = frozenset(tree.non_root())
    sets = TA.admissible_sets(tree)
    assert full in sets and frozenset() in sets
    for S in sets:
        claimed = TA.claimed_rate(tree, S)
        if S != full:
            assert claimed <= n - 2, (S, claimed)
        h1, h2 = TA.scaled_HS(tree, S, 200.0), TA.scaled_HS(tree, S, 400.0)
        assert h1 != 0 and h2 != 0
        measured = len(S) - math.log(abs(h2 / h1)) / 200.0
        assert abs(measured - claimed) < 0.05, (S, measured, claimed)
    assert TA.claimed_rate(tree, full) == n - 1


# ================================================================================== 7. the jump at n = 7, 8
UL, WL = roots_laguerre(80)


def _Rs(x):
    """e^x R(x), R(x) = int_x^oo s K1(s) ds."""
    X = x[:, None] + UL[None, :]
    return (X * k1e(X)) @ WL


def _witness(n: int, m: int, eps: float) -> float:
    """d^m J/d eps^m for the star of n (H^top = (-1)^(n-1) R^(n-1) x^(3-n)), J = int_1^oo I0((n-1-eps)x) H^top dx;
    m = 0 or 1 (d/d eps I0((n-1-eps)x) = -x I1((n-1-eps)x))."""
    y = np.linspace(0.0, np.log(1.0 + 80.0 / eps), 200001)
    x = np.exp(y)
    ker = i0e((n - 1 - eps) * x) if m == 0 else -x * i1e((n - 1 - eps) * x)
    f = ker * (-1.0) ** (n - 1) * _Rs(x) ** (n - 1) * np.exp(-eps * x) * x ** (3.0 - n)
    return float(np.trapezoid(f * x, y))


def test_derivative_order_of_the_witness():
    """The non-holomorphy witness is J for n <= 6 and dJ/d eps for n = 7, 8.  Rejects m = 0 at n = 7, 8."""
    assert [TA.nonholomorphy_derivative_order(n) for n in range(3, 9)] == [0, 0, 0, 0, 1, 1]


@pytest.mark.parametrize("n", [7, 8])
def test_jump_witness_unbounded_at_n7_n8(n):
    """Rejects a check on J instead of its eps-derivative at n = 7, 8: J is bounded there (it changes by < 10%
    between eps = 1e-3 and 1e-5 and its increments shrink, so a J-check would pass a vanishing top term), while the witness of order
    nonholomorphy_derivative_order(n) is unbounded with the predicted amplitude:
    c = C (2 pi (n-1))^(-1/2), C = (-1)^(n-1) (pi/2)^((n-1)/2) (star);
    n = 7: d^m J ~ -c Gamma(1/2) eps^(-1/2); n = 8: d^m J ~ c log eps.  Amplitude to 5%."""
    m = TA.nonholomorphy_derivative_order(n)
    C = (-1) ** (n - 1) * (np.pi / 2) ** ((n - 1) / 2)
    c = C / np.sqrt(2 * np.pi * (n - 1))
    e1, e2 = 1e-3, 1e-5
    w1, w2 = _witness(n, m, e1), _witness(n, m, e2)
    if n == 7:
        pred = -c * math.sqrt(math.pi) * (e2 ** -0.5 - e1 ** -0.5)
    else:
        pred = c * (math.log(e2) - math.log(e1))
    assert abs((w2 - w1) / pred - 1) < 0.05, (w1, w2, pred)
    assert abs(w2) > abs(w1) * (2.0 if n == 7 else 1.1)
    # J itself (order m - 1 = 0) is bounded and settles: a J-based check is blind here
    j1, j15, j2 = _witness(n, 0, e1), _witness(n, 0, 1e-4), _witness(n, 0, e2)
    assert abs(j2 - j1) < 0.1 * abs(j1)
    assert abs(j2 - j15) < 0.5 * abs(j15 - j1)


# ---------------------------------------------------------------------------------------------------------------
# Fire test, 2026-11-18 (plants in memory via a pytest plugin, or on a scratch copy of this file; nothing committed
# carried a plant).  "-k" is the selector run under the plant.
#
# | # | plant (wrong answer)                                         | -k selector                              | result        |
# |:--|:-------------------------------------------------------------|:-----------------------------------------|:--------------|
# | 1a| toy 1 replaced by d/dt o d/dt, solution t (inference "right") | old_inference                            | FIRED         |
# | 1b| toy 2 replaced by exp(t) (holomorphic member of its module)   | old_inference                            | FIRED         |
# | 2 | every stored a -> L_(n-1) a                                   | symbol_roots                             | FIRED 8/8     |
# | 2'| same plant                                                    | aH_is_pure                               | NOT FIRED (expected: L_(n-1) keeps the pure span; the root test is the guard) |
# | 3 | lowest nonzero term of every stored a dropped                 | aH_is_pure                               | FIRED 8/8     |
# | 3'| a H = 0                                                       | (not planted: no in-memory a with a H = 0 at hand) | --    |
# | 4 | phi(z^k) = t^(-2k) theta^(2k)                                 | transfer_identity_star3_float64          | FIRED (rel 1.23) |
# | 4 | same                                                          | transfer_identity_star3 --slow (30 dig.) | FIRED         |
# | 4 | same                                                          | transfer_identity_star4 --slow           | FIRED (rel 1.0) |
# |4b | phi(theta) = 1 + theta                                        | transfer_identity_star3_float64          | FIRED (rel 1.21) |
# |5a | C with k_u = 1 everywhere                                     | top_constant                             | FIRED 10/10   |
# |5b | C with prod k_u                                               | top_constant                             | FIRED 10/10   |
# |5c | C with sign (-1)^n                                            | top_constant                             | FIRED 10/10   |
# | 6 | claimed rate n-1 for every S with |S| = n-2                   | lower_pieces                             | FIRED 4/4 (|S| <= n-2 clause) |
# |6b | claimed rate |S|+1 for |S| = 1 (still <= n-2)                | lower_pieces                             | FIRED (measured-rate clause) |
# | 7 | derivative order 0 for every n (J checked, not dJ/d eps)      | derivative_order or jump_witness         | FIRED 3/3     |
# | - | control, no plant                                             | fast tier                                | passed        |
# ---------------------------------------------------------------------------------------------------------------
