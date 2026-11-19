"""Paper 63, homogeneous l > 0 part of conj:tree_attain: for a Coulomb tree of n densities at the fingerprint weights,
each P_v a combination of h_l(d) nabla^{2j} with l + 2j = d_v <= 2 (centre derivatives only), 3 <= n <= 8, t on
the root and every other scale 1, if the angular integral is nonzero then F_n lies in the D-module of the fibre.
Proof record (derived and cold-reviewed 2026-11-19; this guard is a separate pass):
  (0) relocation: move P_0(d) off the root density K0(t|r|) by 3-D integration by parts; surface terms vanish
      (never two derivatives on K0 in a surface term; for a Laplacian root the excised-ball difference is
      4 pi t^3 r0 Psi(0) -> 0), so the fibre is 4 pi t^3 T(G), G = x^2 <P_0(-d) Psi>;
  (i)-(iii) of prop:star_attain with G in place of H: an operator a in C[z]<theta> with a G in the pure top and
      symbol roots |t| <= n-2; the top term G^top = C x^((5-n)/2) e^{-(n-1)x}(1 + O(1/x)) with the exponent
      UNCHANGED by centre derivatives and C = (-1)^(n-1) (pi/2)^((n-1)/2) Pi / prod_{u != root} k_u^2,
      Pi = <P_0((n-1)omega) prod_{v != root} P_v(-omega)> = (-(n-1))^{d_0} x (the conjecture's angular average)
      for homogeneous P_v; a Z/2 grading (degree p + sum L_kept - L_channel mod 2) makes the pure top even in r.

Each test names the wrong answer it rejects.  Covered:
  1. relocation (own radial potentials, log grid, two cutoffs r0 = 1e-8, 1e-11): a Laplacian root (two s leaves,
     Psi(0) = 1) whose direct-minus-relocated difference is not 4 pi t^3 r0 Psi(0) (an O(1) or non-linear-in-r0
     boundary term); p_z and h_2 roots (Psi(0) = 0) where direct and relocated must agree with no r0 dependence;
  2. the model's relocated G against the hand relocation (h_2 root: (1/5)(theta+1)(theta+3) Phi_2; p_z root:
     -(1/3) x (theta+2) Phi_1) -- rejects P_0(+d) for the odd root and a wrong ladder;
  3. exact operators (over Q) for A (s root, two p_z leaves) and B1 (h_2 root - p_z - p_z): a G pure and nonzero,
     symbols exactly -4 t^3 (t-1)^2 (t+1)^2 and -4 t (t-1)^2 (t+1)^2; B1 theta-order exactly 5 (stored order 5
     and no relation of order 4); rejects a root at +-2, other multiplicities, other orders;
  4. mod 2147483629 for n = 4 (H, J, K) and n = 5 (P): symbol splits into integer roots |r| <= n-2 with the
     stored multiplicities, nonzero at -(n-1); a G pure and nonzero;
  5. parity: every pure coordinate of a G has even r-power (all stored trees); the grading holds on theta^j G;
  6. top constants by an independent exterior-tail quadrature: B1 = pi/15 (rejects pi/60, the missing
     (-(n-1))^{d_0}), D (s root, two h_2 leaves) = pi/10, A = pi/6; local exponent (5-n)/2 = 1 (rejects 1 + derivative
     degrees: 3 for B1, 5 for D); Pi = (-(n-1))^{d_0} x principal average on every stored homogeneous tree;
  7. RECORD only, no claim either way about the conjecture: leaf d_z^2 - 1/3 next to s densities has Pi = 0 and
     principal-symbol integral 4 pi/3; leaf 1 + d_z has Pi = 1 and principal integral 0.

Not asserted (outside the proved statement): non-homogeneous P_v beyond the record in 7; Pi = 0 (the cancelling
root B2 is deliberately absent); scale derivatives d_a; n >= 9; non-axial (m != 0) densities; n = 6..8 operators.

Gaps: the transfer phi(a) T(G) = T(a G) on true functions (debug N3) is not re-run here, so the model's agreement
with the real functions enters only through item 2 (B1 and one p_z root) and item 1; the n = 4, 5 operators are
checked modulo one prime only (a wrong reduction would need the least order to drop mod that prime); top constants
are checked for three n = 3 trees only (the n = 4..6 rows of the derivation record are not re-run); item 3's
minimality is checked for B1 only; the step (ii) branch / non-holomorphy witness is not re-run here (it is the s-type
mechanism of tests/test_paper63_tree_attainment.py, with C as the only l-dependent input).  Everything runs in the
fast tier (~7 s); there is no slow test.

Data: tests/data/paper63_tree_attain_l_2026_11_19.json (stored operators; every property re-verified here).
Nothing is imported from debug/.
"""
from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import cumulative_simpson, simpson
from scipy.special import k0, k1, kn, kve, roots_laguerre

import _tree_attain_l as TL

DATA = Path(TL.__file__).resolve().parent / "data" / "paper63_tree_attain_l_2026_11_19.json"
t = sp.symbols("t")


def _load() -> Dict:
    raw = json.loads(DATA.read_text())
    out = {}
    for name, rec in raw.items():
        if name.startswith("_"):
            continue
        a = {tuple(int(x) for x in key.split(",")): int(v) for key, v in rec["a"].items() if int(v) != 0}
        tree = TL.LTree(rec["n"], [tuple(e) for e in rec["edges"]], rec["root"],
                        {int(k): v for k, v in rec["densities"].items()})
        out[name] = (tree, a, rec["prime"])
    return out


OPS = _load()
EXACT = ["A_pz_s_pz_rootmid", "B1_h2root_pz_pz"]
MODP = ["H_path4_sroot_end_internal_pz", "J_path4_h2root_end", "K_path4_laproot_dz2_nested", "P_path5_h2root_end"]
# symbol root multiplicities on record (derivation table), as {root: multiplicity}
ROOTS = {
    "A_pz_s_pz_rootmid": {0: 3, 1: 2, -1: 2},
    "B1_h2root_pz_pz": {0: 1, 1: 2, -1: 2},
    "H_path4_sroot_end_internal_pz": {0: 4, 1: 2, -1: 2, 2: 2, -2: 2},
    "J_path4_h2root_end": {0: 3, 1: 2, -1: 2, 2: 2, -2: 2},
    "K_path4_laproot_dz2_nested": {0: 5, 1: 4, -1: 4, 2: 4, -2: 4},
    "P_path5_h2root_end": {0: 3, 1: 4, -1: 4, 2: 2, -2: 2, 3: 2, -3: 2},
}


# ================================================================================== 1. relocation, numerically
def _grid(r0: float, r1: float = 90.0, N: int = 400001):
    u = np.linspace(np.log(r0), np.log(r1), N)
    return u, np.exp(u)


def _cum0(f, u):
    return cumulative_simpson(f * np.exp(u), x=u, initial=0.0)


def _potential(W, L: int, u, r):
    """Hat channel-L potential Phi = (1/(2L+1)) [r^-L-1 int_0^r s^{L+2} W + r^L int_r^oo s^{1-L} W], and Phi'."""
    A = _cum0(r ** (L + 2) * W, u)
    c = _cum0(r ** (1 - L) * W, u)
    B = c[-1] - c
    Phi = (r ** (-L - 1) * A + r ** L * B) / (2 * L + 1)
    dPhi = (-(L + 1) * r ** (-L - 2) * A + L * r ** (L - 1) * B) / (2 * L + 1)
    return Phi, dPhi


def _integ(f, u):
    return simpson(f * np.exp(u), x=u)


@lru_cache(maxsize=None)
def relocation_table(r0: float) -> Dict:
    """direct (root derivative on K0(t r), no relocation) and relocated (4 pi t^3 int K0(t r) G dr) at t = 0.7, 1.5,
    2.6, for three n = 3 roots, on a log grid with lower cutoff r0 (own potentials; no model):
      h2:  root h_2 - p_z - p_z (B1): direct t^5 (4pi/5) int K2(tr) Phi_2 r^2, G = (1/5)(-r^2 W_2 + 3 r Phi_2' + 9 Phi_2)
      pz:  root p_z - s - p_z: direct -t^4 (4pi/3) int K1(tr) Phi_1 r^2, G = -(1/3) r (r Phi_1' + 2 Phi_1)
      lap: root lap, two s leaves: direct 4pi t^5 int (K0(tr) - K1(tr)/(tr)) Psi r^2, G = r^2 Psi'' + 2 r Psi'."""
    u, r = _grid(r0)
    phi1, _ = _potential(-k1(r), 1, u, r)                 # p_z leaf, channel 1 (d_z K0 = -mu K1)
    W2 = (2.0 / 3.0) * (-k1(r)) * phi1                    # middle p_z x leaf potential, channel 2
    Phi2, dPhi2 = _potential(W2, 2, u, r)
    G_h2 = (-r ** 2 * W2 + 3 * r * dPhi2 + 9 * Phi2) / 5.0
    Phi1m, dPhi1m = _potential(k0(r) * phi1, 1, u, r)     # middle s x leaf p_z potential, channel 1
    G_pz = -(1.0 / 3.0) * r * (r * dPhi1m + 2 * Phi1m)
    phi0, dphi0 = _potential(k0(r), 0, u, r)
    Psi = phi0 ** 2
    G_lap = 2 * r ** 2 * dphi0 ** 2 + 2 * phi0 * (r ** 2 * (-k0(r)))   # r^2 lap(phi^2), lap phi = -K0
    out = {"Psi0": float(Psi[0])}
    for tt in (0.7, 1.5, 2.6):
        out[tt] = {
            "h2": (tt ** 5 * 4 * np.pi / 5 * _integ(kn(2, tt * r) * Phi2 * r ** 2, u),
                   4 * np.pi * tt ** 3 * _integ(k0(tt * r) * G_h2, u)),
            "pz": (-tt ** 4 * 4 * np.pi / 3 * _integ(k1(tt * r) * Phi1m * r ** 2, u),
                   4 * np.pi * tt ** 3 * _integ(k0(tt * r) * G_pz, u)),
            "lap": (4 * np.pi * tt ** 5 * _integ((k0(tt * r) - k1(tt * r) / (tt * r)) * Psi * r ** 2, u),
                    4 * np.pi * tt ** 3 * _integ(k0(tt * r) * G_lap, u)),
        }
    return out


R0S = (1e-8, 1e-11)
TS = (0.7, 1.5, 2.6)


def test_relocation_laplacian_root_difference_is_the_vanishing_surface_term():
    """Rejects a relocation whose direct-minus-relocated difference for a Laplacian root (Psi(0) = 1 != 0) is not
    the excised-ball surface term 4 pi t^3 r0 Psi(0): an O(1) boundary term, or one that does not shrink linearly
    in r0.  At r0 = 1e-8 and 1e-11 and t = 0.7, 1.5, 2.6 the signed difference equals 4 pi t^3 r0 Psi(0) to 1e-3
    (measured <= 9e-5), and the two cutoffs' differences are in ratio 1000 to 1e-3."""
    tabs = {r0: relocation_table(r0) for r0 in R0S}
    for r0, tab in tabs.items():
        assert abs(tab["Psi0"] - 1.0) < 1e-6             # phi0(0) = int_0^oo s K0(s) ds = 1
        for tt in TS:
            d, rl = tab[tt]["lap"]
            pred = 4 * np.pi * tt ** 3 * r0 * tab["Psi0"]
            assert abs((d - rl) / pred - 1) < 1e-3, (r0, tt, d - rl, pred)
    for tt in TS:
        da = tabs[1e-8][tt]["lap"][0] - tabs[1e-8][tt]["lap"][1]
        db = tabs[1e-11][tt]["lap"][0] - tabs[1e-11][tt]["lap"][1]
        assert abs(da / db / 1000 - 1) < 1e-3, (tt, da, db)


@pytest.mark.parametrize("root", ["pz", "h2"])
def test_relocation_pz_and_d_roots_agree_without_cutoff_dependence(root):
    """For a p_z root (p_z - s - p_z) and an h_2 root (B1), where Psi(0) = 0, direct and relocated agree to 1e-9
    relative at BOTH cutoffs (measured 2e-15 for p_z, <= 3e-11 for h_2, grid-limited), so there is no r0-dependent
    term: a surviving term of size t^3 r0 would be >= 1e-8 relative at r0 = 1e-8.  Rejects a wrong relocated G
    (e.g. (theta+1)^2 in place of (theta+1)(theta+3), or (theta+1) in place of (theta+2))."""
    for r0 in R0S:
        tab = relocation_table(r0)
        for tt in TS:
            d, rl = tab[tt][root]
            assert abs(d) > 1e-2
            assert abs(d - rl) / abs(d) < 1e-9, (r0, tt, d, rl)


# ================================================================================== 2. model G = hand relocation
def _coords(vec: Dict) -> Dict:
    """{(k, i, kept, r-power): Fraction} for a model element (Fraction coefficients)."""
    out = {}
    for (k, i, kept, e), poly in vec.items():
        for m, c in enumerate(poly):
            if c:
                out[(k, i, kept, TL.Model.rpow(e, m))] = c
    return out


def test_model_G_matches_hand_relocation_h2_root():
    """B1 (h_2 root - p_z - p_z): the model's G = x^2 <h_2(-d) Psi> must equal the hand relocation
    (1/5)(theta+1)(theta+3) Phi_{1,2} = (9/5) Phi + (3/5) theta Phi - (1/5) r^2 W_{1,2}  (theta^2 Phi eliminated by
    the channel ODE); the r^2 W part is the next filtration piece.  Rejects a ladder slip (e.g. (theta+1)^2)."""
    tree, _, _ = OPS["B1_h2root_pz_pz"]
    mdl = TL.Model(tree)
    co = _coords(mdl.G())
    assert co.get((0, 0, ((1, 2, 0),), 0)) == sp.Rational(9, 5)
    assert co.get((0, 0, ((1, 2, 1),), 0)) == sp.Rational(3, 5)
    rest = {k: v for k, v in co.items() if k[0] == 0}
    assert len(rest) == 2                                    # only Phi_{1,2}, theta Phi_{1,2} at the k = 0 level


def test_model_G_matches_hand_relocation_pz_root():
    """p_z root - s - p_z (odd root): G = x^2 <(-d_z)(Phi_{1,1} P_1)> = -(1/3) x (theta + 2) Phi_{1,1}.
    Rejects P_0(+d) for an odd root (flips the sign) and a wrong ladder constant."""
    tree = TL.LTree(3, [(0, 1), (1, 2)], 0, {0: "pz", 1: "s", 2: "pz"})
    co = _coords(TL.Model(tree).G())
    assert co == {(0, 0, ((1, 1, 0),), 1): sp.Rational(-2, 3), (0, 0, ((1, 1, 1),), 1): sp.Rational(-1, 3)}, co


# ================================================================================== 3./4. operators: symbols, purity
def _int_roots_exact(coeffs_asc, n: int):
    P = sp.Poly([sp.Integer(c) for c in reversed(coeffs_asc)], t)
    _, facs = sp.factor_list(P)
    roots: Dict[int, int] = {}
    for f, m in facs:
        assert f.degree() == 1, f"non-linear factor {f.as_expr()}"
        r = -f.coeff_monomial(1) / f.coeff_monomial(t)
        assert r.is_integer, r
        roots[int(r)] = roots.get(int(r), 0) + m
    return roots, P.eval(-(n - 1)), P


def _int_roots_modp(coeffs_asc, n: int, P: int):
    poly = [c % P for c in coeffs_asc]
    while poly and poly[-1] == 0:
        poly.pop()
    val = sum(c * pow(-(n - 1), i, P) for i, c in enumerate(poly)) % P
    roots: Dict[int, int] = {}
    for r in range(-(n + 2), n + 3):
        while len(poly) > 1:
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


SYMBOLS = {"A_pz_s_pz_rootmid": -4 * t ** 3 * (t - 1) ** 2 * (t + 1) ** 2,
           "B1_h2root_pz_pz": -4 * t * (t - 1) ** 2 * (t + 1) ** 2}
ORDERS = {"A_pz_s_pz_rootmid": 7, "B1_h2root_pz_pz": 5}


@pytest.mark.parametrize("name", EXACT)
def test_exact_symbol(name):
    """Exact over Q: sigma(a)(1, -t) equals -4 t^3 (t-1)^2 (t+1)^2 (A) / -4 t (t-1)^2 (t+1)^2 (B1), splits into
    integer roots |r| <= n-2 = 1 with the stated multiplicities, is nonzero at -(n-1) = -2, and the theta-order is
    7 (A) / 5 (B1).  Rejects a root at +-2 (e.g. a left-multiplied by L_2), multiplicities other than 0^1, (+-1)^2
    for B1, and a different order."""
    tree, a, prime = OPS[name]
    assert prime is None
    n = tree.n
    roots, at_p, P = _int_roots_exact(TL.symbol_at_1_minus_t(a), n)
    assert sp.expand(P.as_expr() - SYMBOLS[name]) == 0, sp.factor(P.as_expr())
    assert roots == ROOTS[name]
    assert max(abs(r) for r in roots) <= n - 2
    assert at_p != 0
    assert max(j for (_, j) in a) == ORDERS[name]


def test_B1_has_no_relation_of_order_4():
    """Minimality of B1's theta-order 5: the non-pure coordinates of G, theta G, ..., theta^4 G are linearly
    independent over Q(z) (full row rank at z = 7/5; a specialisation can only lower the rank).  Rejects a claimed
    order below 5."""
    tree, _, _ = OPS["B1_h2root_pz_pz"]
    mdl = TL.Model(tree)
    vs = TL.theta_tower(mdl, mdl.G(), 5)
    zv = sp.Rational(7, 5)

    def ev(poly):
        return sum(sp.Rational(c.numerator, c.denominator) * zv ** m for m, c in enumerate(poly) if c)
    for N, full in ((4, True), (5, False)):
        cols = sorted({k for u in vs[:N + 1] for k in u if not mdl.pure(k)}, key=str)
        Mx = sp.Matrix([[ev(u.get(k, [])) for k in cols] for u in vs[:N + 1]])
        rank = Mx.rank()
        assert (rank == N + 1) == full, (N, rank)


@pytest.mark.parametrize("name", MODP)
def test_modp_symbol(name):
    """Modulo 2147483629 (n = 4: H, J, K; n = 5: P): the symbol splits completely into roots among 0, +-1, ...,
    +-(n-2) with the multiplicities on record, and is nonzero at -(n-1).  Rejects a root at +-(n-1) and an operator
    other than the one on record."""
    tree, a, P = OPS[name]
    n = tree.n
    roots, val, deg_left = _int_roots_modp(TL.symbol_at_1_minus_t(a, P), n, P)
    assert deg_left == 0, (roots, deg_left)
    assert max(abs(r) for r in roots) <= n - 2, roots
    assert roots == ROOTS[name], roots
    assert val != 0


@pytest.mark.parametrize("name", EXACT + MODP)
def test_aG_pure_nonzero_and_even(name):
    """a G recomputed in the model from scratch (exactly for n = 3, mod p for n = 4, 5): no non-pure coordinate
    left, at least one pure coordinate nonzero, and every pure coordinate has an EVEN power of r.  Rejects a G not
    pure (a term dropped), a G = 0, and a pure coordinate with an odd r-power (q would not lie in C[z]<theta>)."""
    tree, a, P = OPS[name]
    mdl = TL.Model(tree, P)
    aG = TL.apply_operator(mdl, a)
    nonpure, pure = TL.purity_report(mdl, aG)
    assert nonpure == 0
    assert pure >= 1
    assert TL.pure_r_parities(mdl, aG) == [0]


GRADING_TREES = {
    "A": (3, [(0, 1), (0, 2)], {0: "s", 1: "pz", 2: "pz"}),
    "B1": (3, [(0, 1), (1, 2)], {0: "h2", 1: "pz", 2: "pz"}),
    "C_pz_root": (3, [(0, 1), (0, 2)], {0: "pz", 1: "pz", 2: "s"}),
    "D": (3, [(0, 1), (0, 2)], {0: "s", 1: "h2", 2: "h2"}),
    "G_star4_pz_root": (4, [(0, 1), (0, 2), (0, 3)], {0: "pz", 1: "pz", 2: "s", 3: "s"}),
    "K_lap_root": (4, [(0, 1), (1, 2), (0, 3)], {0: "lap", 1: "s", 2: "dz2", 3: "s"}),
}


@pytest.mark.parametrize("name", list(GRADING_TREES))
def test_grading_even_on_theta_tower(name):
    """The Z/2 grading: every coordinate of theta^j G (j <= 6) has degree p + sum L_kept (channel 0) even.  Rejects
    a G with an odd twist (e.g. x^1 <...> in place of x^2 <...>), which would put odd r-powers in the pure top."""
    n, edges, dens = GRADING_TREES[name]
    mdl = TL.Model(TL.LTree(n, edges, 0, dens), TL.PRIME)
    for v in TL.theta_tower(mdl, mdl.G(), 6):
        assert TL.graded_degree_parities(v) == [0]


# ================================================================================== 6. top constants
VN, VW = roots_laguerre(80)


def _lag(f, x: float, k: float) -> float:
    """int_x^oo f(s) e^{-k(s-x)} ds for f already scaled."""
    v = VN / k
    return float(np.sum(VW / k * f(x + v)))


def _E1(s):
    """e^s x exterior tail of the p_z leaf potential (channel 1, W_1 = -K1)."""
    s = np.atleast_1d(s)
    return np.array([(si ** -2 * _lag(lambda w: w ** 3 * kve(1, w), si, 1.0) - si * kve(0, si)) / 3 for si in s])


def _E2(x: float) -> float:
    """e^x x exterior tail of an h_2 leaf potential (channel 2, W_2 = K2)."""
    A = _lag(lambda s: s ** 4 * kve(2, s), x, 1.0)
    B = _lag(lambda s: s ** -1 * kve(2, s), x, 1.0)
    return (-x ** -3 * A + x ** 2 * B) / 5


def _W2(s):
    return (2.0 / 3.0) * (-kve(1, s)) * _E1(s)


def gtop_scaled(name: str, x: float) -> float:
    """e^{2x} G^top(x), no model and no local rule:
      A  (s root, two p_z leaves):   x^2 <E_1(x) mu . E_1(x) mu> = x^2 E_1^2 / 3;
      D  (s root, two h_2 leaves):   x^2 E_2^2 <P_2^2> = x^2 E_2^2 / 5;
      B1 (h_2 root - p_z - p_z):     (1/5)(-x^2 W_2 + 3 x E_2' + 9 E_2), E_2 the middle's channel-2 tail."""
    if name == "A":
        return x ** 2 * _E1(x)[0] ** 2 / 3
    if name == "D":
        return x ** 2 * _E2(x) ** 2 / 5
    A = _lag(lambda s: s ** 4 * _W2(s), x, 2.0)
    B = _lag(lambda s: s ** -1 * _W2(s), x, 2.0)
    E = (-x ** -3 * A + x ** 2 * B) / 5
    dE = (3 * x ** -4 * A + 2 * x * B) / 5
    return (-x ** 2 * _W2(np.array([x]))[0] + 3 * x * dE + 9 * E) / 5


TOP_TREES = {
    "A": TL.LTree(3, [(0, 1), (0, 2)], 0, {0: "s", 1: "pz", 2: "pz"}),
    "B1": TL.LTree(3, [(0, 1), (1, 2)], 0, {0: "h2", 1: "pz", 2: "pz"}),
    "D": TL.LTree(3, [(0, 1), (0, 2)], 0, {0: "s", 1: "h2", 2: "h2"}),
}
TOP_EXPECT = {"A": sp.pi / 6, "B1": sp.pi / 15, "D": sp.pi / 10}


@lru_cache(maxsize=None)
def _top_numeric(name: str) -> Tuple[float, float]:
    """(C by Richardson on x = 50..400 with the exponent of TL.top_exponent, local exponent from x = 400 -> 800)."""
    b = float(TL.top_exponent(TOP_TREES[name]))
    xs = [50.0, 100.0, 200.0, 400.0]
    vals = [gtop_scaled(name, x) / x ** b for x in xs]
    r1 = [2 * vals[i + 1] - vals[i] for i in range(3)]
    r2 = [(4 * r1[i + 1] - r1[i]) / 3 for i in range(2)]
    C = (8 * r2[1] - r2[0]) / 7
    slope = math.log(gtop_scaled(name, 800.0) / gtop_scaled(name, 400.0)) / math.log(2.0)
    return C, slope


@pytest.mark.parametrize("name", list(TOP_TREES))
def test_top_constant_formula(name):
    """The formula C = (-1)^(n-1) (pi/2)^((n-1)/2) Pi / prod k_u^2 gives exactly pi/6 (A), pi/15 (B1), pi/10 (D).
    Rejects C built from the principal-symbol average without (-(n-1))^{d_0} (pi/60 for B1)."""
    C = TL.top_constant(TOP_TREES[name])
    assert sp.simplify(C - TOP_EXPECT[name]) == 0, C


@pytest.mark.parametrize("name", list(TOP_TREES))
def test_top_constant_and_exponent_numeric(name):
    """Independent exterior-tail quadrature of G^top (own code; no model, no local rule): C to 1e-6 (measured
    2e-9 A and B1, 7e-8 D) and the local exponent (5-n)/2 = 1 to 0.02 (measured 0.995 A, 0.9995 B1, 0.990 D).  Rejects C = pi/60
    for B1 (off by 4), any D constant other than pi/10, and the shifted exponent 1 + (derivative degrees): 3 for
    B1, 5 for D, 3 for A."""
    tree = TOP_TREES[name]
    C, slope = _top_numeric(name)
    pred = float(TOP_EXPECT[name])
    assert abs(C / pred - 1) < 1e-6, (C, pred)
    assert abs(slope - float(TL.top_exponent(tree))) < 0.02, slope
    assert TL.top_exponent(tree) == 1
    if name == "B1":
        assert abs(C / (math.pi / 60) - 1) > 1.0


@pytest.mark.parametrize("name", list(OPS) + ["D", "C_pz_root", "G_star4_pz_root"])
def test_Pi_equals_signed_principal_average_for_homogeneous(name):
    """For homogeneous P_v, Pi = (-(n-1))^{d_0} x <prod_v sigma_v(-omega)>, and Pi != 0 on every tree used here
    (B1: 8/15 = 4 x 2/15).  Rejects Pi read as the principal average alone."""
    if name in OPS:
        tree = OPS[name][0]
    else:
        n, edges, dens = GRADING_TREES[name]
        tree = TL.LTree(n, edges, 0, dens)
    assert all(TL.is_homogeneous(tree.dens[v]) for v in range(tree.n))
    d0 = TL.degree(tree.dens[tree.root])
    Pi = TL.Pi(tree)
    assert Pi != 0
    assert sp.simplify(Pi - (-(tree.n - 1)) ** d0 * TL.principal_average(tree)) == 0
    if name == "B1_h2root_pz_pz":
        assert Pi == sp.Rational(8, 15) and TL.principal_average(tree) == sp.Rational(2, 15)


# ================================================================================== 7. non-homogeneous record
def test_record_nonhomogeneous_leaves_pi_versus_principal_integral():
    """RECORD ONLY -- nothing here says the conjecture fails or holds for non-homogeneous P_v.  Next to s densities
    (s root, leaves X and s, n = 3): leaf X = d_z^2 - 1/3 has Pi = 0 but principal-symbol integral
    int omega_z^2 dOmega = 4 pi/3; leaf X = 1 + d_z has Pi = 1 but principal integral int (-omega_z) dOmega = 0.
    Rejects a guard (or a statement) that takes Pi and the principal-symbol integral to be the same quantity."""
    for leaf, pi_exp, princ_int in (("dz2_m_third", 0, 4 * sp.pi / 3), ("one_p_dz", 1, 0)):
        tree = TL.LTree(3, [(0, 1), (0, 2)], 0, {0: "s", 1: leaf, 2: "s"})
        assert not TL.is_homogeneous(leaf)
        assert TL.Pi(tree) == pi_exp
        assert sp.simplify(4 * sp.pi * TL.principal_average(tree) - princ_int) == 0
        assert sp.simplify(4 * sp.pi * TL.Pi(tree) - princ_int) != 0


# ---------------------------------------------------------------------------------------------------------------
# Fire test, 2026-11-19 (each plant = a scratch copy of this file with module-level overrides appended, or one text
# substitution, run with "-k"; nothing in the repo carried a plant).  Control (no plant): 40 passed.
#
# | #   | plant (wrong answer)                                                   | -k selector                  | result |
# |:----|:-----------------------------------------------------------------------|:-----------------------------|:-------|
# | P1  | lap root: O(1) boundary term (4 pi t^3 Psi0 x 1e-3, r0-independent)    | relocation_laplacian         | FIRED  |
# | P1b | lap root: boundary term ~ sqrt(r0) instead of r0                       | relocation_laplacian         | FIRED  |
# | P1c | lap root: surface term doubled (8 pi t^3 r0 Psi0)                      | relocation_laplacian         | FIRED  |
# | P2a | p_z root: relocated G with (theta+1) for (theta+2)                     | relocation_pz_and_d and pz   | FIRED  |
# | P2b | h_2 root: relocated G = (1/5)(theta+1)^2 Phi_2 (3 rPhi'+9 Phi -> 2, 3) | relocation_pz_and_d and h2   | FIRED  |
# | P2c | p_z, h_2 roots: an r0-linear term 4 pi t^3 r0 added                    | relocation_pz_and_d          | FIRED 2/2 |
# | P3a | model G with P_0(+d)                                                   | model_G_matches              | FIRED (p_z root); h_2 root NOT FIRED (even root, expected) |
# | P3b | model ladder slip: d_z down-coefficient (L+1)/(2L+1)                   | model_G_matches              | FIRED 2/2 |
# | P4  | B1: a -> L_2 a (symbol roots at +-2)                                   | exact_symbol or aG_pure      | FIRED exact_symbol[B1]; aG_pure NOT FIRED (expected: L_2 keeps the pure span; the symbol test is the guard) |
# | P5  | A, B1: a -> theta a (extra root 0, order + 1)                           | exact_symbol                 | FIRED 2/2 |
# | P5b | minimality: relation claimed at order 4                                | no_relation_of_order_4       | FIRED  |
# | P6a | mod p: a -> L_(n-1) a for J (n=4), P (n=5)                              | modp_symbol                  | FIRED 2/2 (H, K untouched, passed) |
# | P6b | every stored a: lowest nonzero term dropped                            | aG_pure                      | FIRED 6/6 |
# | P6c | every stored a: a G = 0                                                | aG_pure                      | FIRED 6/6 |
# | P7  | a G multiplied by r (odd r-power in the pure top)                      | aG_pure                      | FIRED 6/6 |
# | P7b | G = x^1 <...> (odd twist)                                              | grading                      | FIRED 6/6 |
# | P8  | Pi = principal average (B1: C = pi/60)                                 | top_constant_formula or Pi_equals | FIRED 7 (B1 formula; Pi test on every tree with d_0 > 0); d_0 = 0 trees NOT FIRED (expected, factor 1) |
# | P8b | C with sign (-1)^n (D -> -pi/10)                                       | top_constant_formula         | FIRED 3/3 |
# | P8c | numeric check expects pi/60 for B1                                     | top_constant_and_exponent_numeric | FIRED (B1) |
# | P8d | numeric check expects pi/9 for D                                       | top_constant_and_exponent_numeric | FIRED (D) |
# | P9  | top exponent 1 + derivative degrees (3 A, 3 B1, 5 D)                   | top_constant_and_exponent_numeric | FIRED 3/3 |
# | P10a| record: Pi taken = principal average                                   | record_nonhomogeneous        | FIRED  |
# | P10b| record: principal average taken = Pi                                   | record_nonhomogeneous        | FIRED  |
# ---------------------------------------------------------------------------------------------------------------
