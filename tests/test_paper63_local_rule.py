"""Guard for the [ARGUMENT]/[MEASURED] paragraph under Paper 63's conjecture conj:tree_attain ("trees reach their
rung"), papers/group3_foundations/paper_63_density_graph_ladder.tex.  Written as a separate pass after the paragraph.

The conjecture is NOT tested here (nothing in this file computes an operator or a composition factor).  What is
tested is the four statements the paragraph leans on:

 (1) the local rule: for a source e^{-kr} r^a sum_j s_j r^-j in angular channel L, the exponential part of its
     Coulomb potential is the formal solution of (theta-L)(theta+L+1) Phi = -4 pi r^2 s; leading term -4 pi s_0/k^2,
     independent of L (L first enters at relative order r^-2);
 (2) the rule reproduces the constants on record (P63L1 leaf tail; p_z chain 1/(3 pi^2), 47/8; s chain 1/pi^2, 31/8;
     star of five pi^2/(8 sqrt 2)) and gives 1/(12 pi), 15/8 for the chain with a d root;
 (3) star of pure harmonics: every density and potential is radial x harmonic (l <= 2), so the fibre is a radial
     integral times int prod h dOmega.  The factorisation and the angular integrals are tested; NO fibre is computed;
 (4) the cancelling root (lap - (5/2) h_2) has int cos^2 (1 - (5/2) P2) dOmega = 0 and its leading branch coefficient
     is exactly zero.

Independent routes (so the rule is not checked against itself):
  - TRUE potential: the exponential part of the channel-L potential is
        (4 pi/(2L+1)) [ r^L int_r^oo u^(1-L) s du - r^(-L-1) int_r^oo u^(L+2) s du ]
    (the full potential minus its multipole term r^(-L-1) x moment).  For s = e^{-ku} u^(a-j) this is a sum of upper
    incomplete gamma functions (mpmath, 40 digits): each coefficient phi_0..phi_5 of the rule is recovered from it;
  - termwise integration by parts of the same bracket (exact rationals) against the ODE recursion;
  - the rule's series substituted back into the differential operator by sympy differentiation;
  - P63L1's closed form E = (2/pi)[K1 + R/r^2], R = int_r^oo s K1(s) ds, against a direct quadrature of the bracket.
Transcribed (record constants, compared with what is derived here): 11/8, 97/128, -463/1024; 1/(3 pi^2), 47/8;
1/pi^2, 31/8; pi^2/(8 sqrt 2); 1/(12 pi), 15/8; 3/(8 pi).  Taken from tests/test_paper63_attainment_l2.py without
re-proof: the root-kernel jumps I0 + I1/y and I2 at t = -y + i0, and the root identities h_2 K0 = K2 P2.

Not duplicated: the measured branches (float/mp scans) are in test_paper63_attainment_l2.py and
test_paper63_star_attainment.py.  Nothing is imported from debug/ or read from any JSON.
"""
from __future__ import annotations

import itertools
import re
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import mpmath as mp
import pytest
import sympy as sp
from sympy.physics.wigner import gaunt

ROOT = Path(__file__).resolve().parent.parent
TEX = ROOT / "papers" / "group3_foundations" / "paper_63_density_graph_ladder.tex"

# ------------------------------------------------------------------------------------------- record constants
TAIL_E = [sp.Integer(1), sp.Rational(11, 8), sp.Rational(97, 128), sp.Rational(-463, 1024)]   # P63L1 leaf tail
A_PZ = 1 / (3 * sp.pi ** 2)            # p_z chain (record P63L1)
KAPPA_PZ = sp.Rational(47, 8)
A_S = 1 / sp.pi ** 2                   # s chain (record SA / P63L1 toy)
KAPPA_S = sp.Rational(31, 8)
C5 = sp.pi ** 2 / (8 * sp.sqrt(2))     # star of five, c_5 Gamma(1/2)
A_D = 1 / (12 * sp.pi)                 # chain with a d root (record P63L2, shape B1)
KAPPA_D = sp.Rational(15, 8)
A_CANCEL = 3 / (8 * sp.pi)             # cancelling root, eps^{-1/2} amplitude (record P63L2, shape B2)
ROOT_COEF = sp.Rational(5, 2)          # the cancelling root is lap - ROOT_COEF h_2(d)

NT = 7                                 # series terms kept
RHS_FACTOR = -4 * sp.pi                # right side of the radial Poisson equation: RHS_FACTOR r^2 s


# =========================================================================================== the local rule
def local_rule(k, a, L, s: Sequence) -> List:
    """Coefficients phi_j of the formal solution e^{-kr} r^a sum_j phi_j r^-j of
    (theta - L)(theta + L + 1) Phi = RHS_FACTOR r^2 e^{-kr} r^a sum_j s_j r^-j,  theta = r d/dr.
    (theta^2 + theta - L(L+1)) e^{-kr} r^b = [k^2 r^2 - 2k(b+1) r + b(b+1) - L(L+1)] e^{-kr} r^b."""
    phi: List = []
    for j in range(len(s)):
        acc = RHS_FACTOR * s[j]
        if j >= 1:
            b = a - (j - 1)
            acc += 2 * k * (b + 1) * phi[j - 1]
        if j >= 2:
            b = a - (j - 2)
            acc -= (b * (b + 1) - L * (L + 1)) * phi[j - 2]
        phi.append(sp.expand(acc / k ** 2))
    return phi


def upper_tail(p, k, n: int) -> List:
    """int_r^oo u^p e^{-ku} du = e^{-kr} r^p sum_j t_j r^-j, t_j = p(p-1)..(p-j+1)/k^(j+1) (integration by parts)."""
    out, ff = [], sp.Integer(1)
    for j in range(n):
        out.append(ff / k ** (j + 1))
        ff *= (p - j)
    return out


def bracket_series(k, a, L, s: Sequence) -> List:
    """Exponential part of the channel-L potential by termwise integration by parts of
    (4 pi/(2L+1)) [r^L int_r^oo u^(1-L) s - r^(-L-1) int_r^oo u^(L+2) s]: coefficients of e^{-kr} r^(a+1-m),
    m = 0..len(s) (the m = 0 coefficient must cancel).  Does not use the differential equation."""
    n = len(s) + 1
    out = [sp.Integer(0)] * n
    for j, sj in enumerate(s):
        t1 = upper_tail(1 - L + a - j, k, n)
        t2 = upper_tail(L + 2 + a - j, k, n)
        for i in range(n - j):
            out[j + i] += 4 * sp.pi / (2 * L + 1) * sj * (t1[i] - t2[i])
    return [sp.expand(x) for x in out]


def true_exp_part(k, a, L: int, s: Sequence, r) -> mp.mpf:
    """e^{kr} r^-a times the exponential part of the true channel-L potential of s = e^{-ku} u^a sum s_j u^-j at r,
    from upper incomplete gamma functions: int_r^oo u^p e^{-ku} du = Gamma(p+1, kr)/k^(p+1)."""
    k, a, r = mp.mpf(k), mp.mpf(a), mp.mpf(r)
    tot = mp.mpf(0)
    for j, sj in enumerate(s):
        p1, p2 = 1 - L + a - j, L + 2 + a - j
        g1 = mp.gammainc(p1 + 1, k * r) / k ** (p1 + 1)
        g2 = mp.gammainc(p2 + 1, k * r) / k ** (p2 + 1)
        tot += mp.mpf(sj) * (r ** L * g1 - r ** (-L - 1) * g2)
    return 4 * mp.pi / (2 * L + 1) * tot * mp.exp(k * r) * r ** (-a)


def _mpq(x) -> mp.mpf:
    with mp.workdps(45):
        return mp.mpf(str(sp.N(x, 45)))


# ------------------------------------------------------------------------------------------- (1) the rule itself
def test_leading_term_and_where_L_enters_symbolic():
    """Generic k, a, L, s_j: phi_0 = -4 pi s_0/k^2; phi_0 and phi_1 do not contain L; phi_2 does, through
    L(L+1) phi_0/k^2 only.  Restates the recursion's own algebra (it would accept any recursion with this shape);
    the tests against the true potential below are the ones that can tell a wrong rule."""
    k, a, L, s0, s1, s2 = sp.symbols("k a L s0 s1 s2")
    phi = local_rule(k, a, L, [s0, s1, s2])
    assert sp.simplify(phi[0] + 4 * sp.pi * s0 / k ** 2) == 0
    assert not phi[0].has(L) and not sp.simplify(phi[1]).has(L)
    assert sp.simplify(sp.diff(phi[2], L) - (2 * L + 1) * phi[0] / k ** 2) == 0
    assert sp.simplify(phi[1] - (-4 * sp.pi * s1 / k ** 2 + 2 * (a + 1) * phi[0] / k)) == 0


@pytest.mark.parametrize("k,a,L", [(1, sp.Rational(-1, 2), 1), (2, -1, 2), (sp.Rational(3, 2), sp.Rational(1, 3), 0),
                                   (3, 2, 3)])
def test_series_solves_the_radial_equation(k, a, L):
    """Substitute the truncated series into r^2 Phi'' + 2 r Phi' - L(L+1) Phi + 4 pi r^2 s by sympy differentiation:
    the residual has no power above r^(a+2-NT) (two orders of truncation remain).  Rejects a wrong coefficient in
    the recursion (2k(b+1), b(b+1), the sign of L(L+1)) and a wrong right side."""
    r = sp.symbols("r", positive=True)
    s = [sp.Integer(1), sp.Rational(3, 10), sp.Rational(-7, 10)] + [sp.Integer(0)] * (NT - 3)
    phi = local_rule(k, a, L, s)
    Phi = sp.exp(-k * r) * sum(c * r ** (a - j) for j, c in enumerate(phi))
    src = sp.exp(-k * r) * sum(c * r ** (a - j) for j, c in enumerate(s))
    res = sp.expand(sp.simplify((r ** 2 * sp.diff(Phi, r, 2) + 2 * r * sp.diff(Phi, r) - L * (L + 1) * Phi
                                 + 4 * sp.pi * r ** 2 * src) * sp.exp(k * r) * r ** (-a)))
    pw = [sp.Rational(t.as_coeff_exponent(r)[1]) for t in sp.Add.make_args(res) if t != 0]
    assert pw and max(pw) <= 2 - NT, pw
    assert len(pw) <= 2, pw            # only the truncation terms r^(2-NT), r^(1-NT) (one may vanish)


@pytest.mark.parametrize("L", [0, 1, 2, 3])
def test_rule_equals_integration_by_parts(L):
    """Exact, to NT orders, generic rational k and a: the recursion's phi_j equal the coefficients obtained by
    integrating the multipole bracket by parts term by term (no differential equation), and the bracket's
    r^(a+1) e^{-kr} term cancels.  Rejects -4 pi s_0/k (k = 3/2, 2), a sign flip, L entering at order one."""
    s = [sp.Integer(1), sp.Rational(3, 10), sp.Rational(-7, 10), sp.Rational(2, 3), 0, 0, 0]
    for k, a in ((sp.Rational(3, 2), sp.Rational(-1, 2)), (sp.Integer(2), sp.Integer(-1)),
                 (sp.Integer(1), sp.Rational(5, 3))):
        br = bracket_series(k, a, L, s)
        assert br[0] == 0
        assert [sp.simplify(x - y) for x, y in zip(local_rule(k, a, L, s), br[1:])] == [0] * NT, (k, a, L)


CASES = [(2, sp.Rational(-1, 2), [1]), (sp.Rational(3, 2), -1, [1, sp.Rational(3, 10), sp.Rational(-7, 10)]),
         (1, sp.Rational(1, 2), [sp.Rational(-2, 3), 1]), (3, 0, [1, 0, 2])]


@pytest.mark.parametrize("L", [0, 1, 2])
@pytest.mark.parametrize("case", range(len(CASES)))
def test_rule_against_true_potential(case, L):
    """The rule's coefficients, recovered one at a time from the TRUE exponential part (incomplete gamma, 40
    digits): with T(r) = e^{kr} r^-a x (true), r^N [T(r) - sum_{j<N} phi_j r^-j] -> phi_N, N = 0..5, at r = 200 and
    400, within 3 (|phi_{N+1}|/r + |phi_{N+2}|/r^2).  Rejects -4 pi s_0/k, +4 pi, an L(L+1) or 2k(b+1) slip, and a
    rule that drops the homogeneous multipole term incorrectly (the true object here is the full potential minus
    r^(-L-1) x moment)."""
    k, a, s = CASES[case]
    k, a = sp.Rational(k), sp.Rational(a)
    src = list(s) + [0] * (NT + 1 - len(s))
    phi = [_mpq(c) for c in local_rule(k, a, L, src)]
    with mp.workdps(40):
        for r in (200, 400):
            T = true_exp_part(_mpq(k), _mpq(a), L, [_mpq(x) for x in s], r)
            part = mp.mpf(0)
            for N in range(6):
                got = (T - part) * mp.mpf(r) ** N
                tol = 3 * (abs(phi[N + 1]) / r + abs(phi[N + 2]) / r ** 2) + mp.mpf(10) ** -25
                assert abs(got - phi[N]) <= tol, (case, L, r, N, got, phi[N])
                part += phi[N] * mp.mpf(r) ** (-N)
            assert abs(phi[0] + 4 * mp.pi * _mpq(s[0]) / _mpq(k) ** 2) < mp.mpf(10) ** -35


@pytest.mark.parametrize("case", range(len(CASES)))
def test_true_potentials_differ_between_channels_only_at_second_order(case):
    """From the true potentials alone (no recursion): T_L(r) - T_0(r) is O(r^-2), with
    r^2 (T_L - T_0) -> -4 pi s_0 L(L+1)/k^4 for L = 1, 2, 3 (so it is not zero either).  Rejects 'L enters at
    order r^-1' and 'the exponential part is L-independent to all orders'."""
    k, a, s = CASES[case]
    with mp.workdps(40):
        kk, aa, ss = _mpq(sp.Rational(k)), _mpq(sp.Rational(a)), [_mpq(x) for x in s]
        for L in (1, 2, 3):
            lim = -4 * mp.pi * ss[0] * L * (L + 1) / kk ** 4
            prev = None
            for r in (400, 1600):
                d = (true_exp_part(kk, aa, L, ss, r) - true_exp_part(kk, aa, 0, ss, r))
                assert abs(d * r) < 5 * abs(lim) / r * 1.5 + mp.mpf(10) ** -25        # no r^-1 term
                err = abs(d * r ** 2 - lim)
                assert err < 40 * abs(lim) / r, (case, L, r, d * r ** 2, lim)
                assert prev is None or err < prev / 2
                prev = err


# =========================================================================================== Bessel tails
def hankel(nu: int, n: int = NT) -> List[sp.Rational]:
    """a_m(nu): K_nu(x) ~ sqrt(pi/(2x)) e^-x sum a_m x^-m,  I_nu(x) ~ e^x (2 pi x)^-1/2 sum (-1)^m a_m x^-m."""
    out = [sp.Integer(1)]
    for m in range(1, n):
        out.append(out[-1] * sp.Rational(4 * nu ** 2 - (2 * m - 1) ** 2, 8 * m))
    return out


class Tail:
    """pref e^{-kx} x^a sum_j c_j x^-j, c_j rational."""

    def __init__(self, pref, k, a, c: Sequence) -> None:
        self.pref, self.k, self.a, self.c = pref, sp.Rational(k), sp.Rational(a), list(c)[:NT]

    def __mul__(self, o: "Tail") -> "Tail":
        c = [sum(self.c[i] * o.c[j - i] for i in range(j + 1)) for j in range(NT)]
        return Tail(self.pref * o.pref, self.k + o.k, self.a + o.a, c)

    def times(self, f, shift=0) -> "Tail":
        return Tail(self.pref * f, self.k, self.a + shift, self.c)

    def potential(self, L: int) -> "Tail":
        """Exponential part of the channel-L Coulomb potential (the local rule)."""
        phi = local_rule(self.k, self.a, L, [c / RHS_FACTOR for c in self.c])
        return Tail(self.pref * RHS_FACTOR, self.k, self.a, phi)

    def lead(self):
        return sp.simplify(self.pref * self.c[0])

    def monic(self) -> List:
        return [sp.nsimplify(c / self.c[0]) for c in self.c]


def K_tail(nu: int) -> Tail:
    return Tail(sp.sqrt(sp.pi / 2), 1, sp.Rational(-1, 2), hankel(nu))


def rho_s() -> Tail:
    """s density at the fingerprint weight, scale one: K0(r)/(2 pi^2)."""
    return K_tail(0).times(1 / (2 * sp.pi ** 2))


def f_p() -> Tail:
    """radial factor of the p_z density d_z [K0/(2 pi^2)] = f_p cos(theta): -K1/(2 pi^2)."""
    return K_tail(1).times(-1 / (2 * sp.pi ** 2))


def R_tail() -> Tail:
    """R(x) = int_x^oo s K1(s) ds, expanded by integrating K1's Hankel series by parts (no local rule)."""
    a1 = hankel(1)
    d = [sp.Integer(0)] * NT
    for m in range(NT):
        t = upper_tail(sp.Rational(1, 2) - m, sp.Integer(1), NT - m)
        for i in range(NT - m):
            d[m + i] += a1[m] * t[i]
    return Tail(sp.sqrt(sp.pi / 2), 1, sp.Rational(1, 2), d)


def same(x: Tail, y: Tail) -> bool:
    return (x.k == y.k and x.a == y.a
            and all(sp.simplify(x.pref * u - y.pref * v) == 0 for u, v in zip(x.c, y.c)))


def branch(nu_terms: Sequence[Tuple[object, List, Tail]], kthr: int, order: int = 2) -> Dict:
    """Singular part at eps -> 0 of sum_i w_i int_1^oo kern_i((kthr - eps) x) H_i(x) dx, kern ~ e^y (2 pi y)^-1/2
    sum c_m y^-m and H_i ~ e^{-kthr x}: uses int_1^oo x^q e^{-eps x} dx = Gamma(q+1) eps^(-q-1) + analytic
    (q not a negative integer).  Returns {power of eps: coefficient}, zero coefficients removed."""
    e = sp.symbols("epsilon", positive=True)
    tot = sp.Integer(0)
    for w, kern, H in nu_terms:
        assert H.k == kthr
        for m in range(order + 1):
            ser = sp.series((kthr - e) ** (-sp.Rational(1, 2) - m), e, 0, order + 1).removeO()
            for j in range(order + 1 - m):
                q = H.a - sp.Rational(1, 2) - m - j
                assert not (q.is_integer and q < 0)
                tot += w * H.pref * H.c[j] * kern[m] / sp.sqrt(2 * sp.pi) * ser * sp.gamma(q + 1) * e ** (-q - 1)
    out: Dict = {}
    for term in sp.Add.make_args(sp.expand(tot)):
        c, p = term.as_coeff_exponent(e)
        out[p] = out.get(p, 0) + c
    q0 = max(H.a for _, _, H in nu_terms) - sp.Rational(1, 2)
    keep = {p: sp.simplify(v) for p, v in out.items() if p <= -q0 - 1 + order}
    return {p: v for p, v in sorted(keep.items()) if v != 0}


def kern_I(nu: int, plus_I1_over_y: bool = False) -> List:
    c = [(-1) ** m * x for m, x in enumerate(hankel(nu))]
    if plus_I1_over_y:
        c1 = [(-1) ** m * x for m, x in enumerate(hankel(1))]
        c = [c[j] + (c1[j - 1] if j >= 1 else 0) for j in range(NT)]
    return c


def test_hankel_coefficients_against_mpmath():
    """The Bessel tails every constant below is built from: 7-term Hankel sums against mpmath's K_nu, I_nu at
    x = 60, nu = 0, 1, 2, relative error below 1e-9 (truncation ~ a_7/60^7).  Rejects a wrong a_m recurrence and a
    missing (-1)^m in I_nu."""
    with mp.workdps(30):
        x = mp.mpf(60)
        for nu in (0, 1, 2):
            a = [_mpq(c) for c in hankel(nu)]
            Ks = mp.sqrt(mp.pi / (2 * x)) * mp.exp(-x) * mp.fsum(c * x ** (-m) for m, c in enumerate(a))
            Is = mp.exp(x) / mp.sqrt(2 * mp.pi * x) * mp.fsum(_mpq(c) * x ** (-m) for m, c in enumerate(kern_I(nu)))
            assert abs(Ks / mp.besselk(nu, x) - 1) < 1e-9 and abs(Is / mp.besseli(nu, x) - 1) < 1e-9, nu


# ------------------------------------------------------------------------------------------- (2) T1, T1b
def test_T1_leaf_tail():
    """E = exponential part of the L = 1 potential of f_p: sqrt(2/pi) r^-1/2 e^-r [1 + 11/8r + 97/128r^2 -
    463/1024r^3], and its lead is -4 pi x (lead of f_p) (the 'r^(1/2) e^-r cancellation').  Independent route: the
    closed form E = (2/pi)[K1 + R/r^2] with R expanded by parts, to NT orders.  Rejects any recursion slip that
    moves 11/8 (e.g. L at order one gives 19/8)."""
    E = f_p().potential(1)
    assert sp.simplify(E.lead() - sp.sqrt(2 / sp.pi)) == 0
    assert sp.simplify(E.lead() + 4 * sp.pi * f_p().lead()) == 0
    assert E.monic()[:4] == TAIL_E
    R = R_tail()
    closed_c = [K_tail(1).c[j] + (R.c[j - 1] if j >= 1 else 0) for j in range(NT)]      # K1 + R/x^2, x^-1/2 e^-x
    assert same(E, Tail(2 / sp.pi * sp.sqrt(sp.pi / 2), 1, sp.Rational(-1, 2), closed_c))


def test_T1_closed_form_is_the_true_potential():
    """Numerically (22 digits), at r = 3 and 25: the multipole bracket of f_p = -K1/(2 pi^2) in channel 1, by direct
    quadrature, equals (2/pi)[K1(r) + R(r)/r^2] to 1e-15; and at r = 25 the four-term record tail is within
    2 |c_4| r^-4 of it.  Rejects a closed form or tail that is not the potential of this density."""
    with mp.workdps(22):
        fp = lambda u: -mp.besselk(1, u) / (2 * mp.pi ** 2)                                # noqa: E731
        c4 = abs(_mpq(f_p().potential(1).monic()[4]))
        for r in (mp.mpf(3), mp.mpf(25)):
            pts = [r, r + 5, r + 20, r + 60, mp.inf]
            br = 4 * mp.pi / 3 * (r * mp.quad(fp, pts) - r ** -2 * mp.quad(lambda u: u ** 3 * fp(u), pts))
            R = mp.quad(lambda u: u * mp.besselk(1, u), pts)
            closed = 2 / mp.pi * (mp.besselk(1, r) + R / r ** 2)
            assert abs(br / closed - 1) < mp.mpf(10) ** -15, (r, br, closed)
        r = mp.mpf(25)
        tail = mp.sqrt(2 / mp.pi) * r ** -0.5 * mp.exp(-r) * mp.fsum(_mpq(c) * r ** (-j) for j, c in enumerate(TAIL_E))
        assert abs(tail / closed - 1) < 2 * c4 * r ** -4


def test_T1b_s_leaf_tail():
    """Exponential part of the s-leaf potential (L = 0 potential of K0/(2 pi^2)) is -(2/pi) R/r, to NT orders, with R
    expanded by parts.  Rejects a wrong sign or normalisation of the rule on the s leaf."""
    assert same(rho_s().potential(0), R_tail().times(-2 / sp.pi, shift=-1))


# ------------------------------------------------------------------------------------------- angular integrals
MU = sp.symbols("mu")
P2 = (3 * MU ** 2 - 1) / 2


def axial(f) -> sp.Expr:
    """int_{S^2} f(cos theta) dOmega."""
    return 2 * sp.pi * sp.integrate(sp.expand(f), (MU, -1, 1))


def sphere_monomial(i: int, j: int, k: int) -> sp.Expr:
    """int_{S^2} x^i y^j z^k dOmega = 2 Gamma((i+1)/2) Gamma((j+1)/2) Gamma((k+1)/2)/Gamma((i+j+k+3)/2), else 0."""
    if i % 2 or j % 2 or k % 2:
        return sp.Integer(0)
    g = lambda n: sp.gamma(sp.Rational(n + 1, 2))                                           # noqa: E731
    return sp.simplify(2 * g(i) * g(j) * g(k) / sp.gamma(sp.Rational(i + j + k + 3, 2)))


X, Y, Z = sp.symbols("x y z")
HARM = {"1": sp.Integer(1), "x": X, "y": Y, "z": Z, "z2": (3 * Z ** 2 - (X ** 2 + Y ** 2 + Z ** 2)) / 2, "xz": X * Z,
        "yz": Y * Z, "xy": X * Y, "x2-y2": X ** 2 - Y ** 2}
DEG = {"1": 0, "x": 1, "y": 1, "z": 1, "z2": 2, "xz": 2, "yz": 2, "xy": 2, "x2-y2": 2}


def sphere_integral(poly) -> sp.Expr:
    return sp.simplify(sum(c * sphere_monomial(*m) for m, c in sp.Poly(sp.expand(poly), X, Y, Z).terms()))


SYM_S = sp.Integer(1)                       # symbol of an s density
SYM_PZ = MU                                 # d_z at -omega (overall sign (-1)^degree, even in every case below)
SYM_D = P2                                  # h_2(d) = (3 d_z^2 - lap)/2
SYM_CANCEL = 1 - ROOT_COEF * P2             # lap - (5/2) h_2(d)


def test_angular_integrals_named_in_the_paper():
    """s-type: 4 pi.  p_z chain (two p_z, s root): 4 pi/3, one third of 4 pi.  p_z - p_z - d chain: 8 pi/15, not
    zero.  Cancelling root: int cos^2 (1 - (5/2) P2) dOmega = 0 exactly, while each of its two channels is 4 pi/3.
    Rejects a root coefficient other than 5/2 and a zero that comes from both channels vanishing."""
    assert axial(SYM_S) == 4 * sp.pi
    assert axial(SYM_PZ * SYM_PZ * SYM_S) == 4 * sp.pi / 3
    assert axial(SYM_PZ * SYM_PZ * SYM_D) == 8 * sp.pi / 15
    assert axial(SYM_PZ * SYM_PZ * SYM_CANCEL) == 0
    assert axial(MU ** 2) == ROOT_COEF * axial(MU ** 2 * P2) == 4 * sp.pi / 3
    assert sp.expand(MU ** 2 - (sp.Rational(1, 3) + sp.Rational(2, 3) * P2)) == 0
    # the monomial formula used for the non-axial harmonics agrees with the axial integrals
    assert sphere_integral(Z ** 2 * HARM["z2"]) == 8 * sp.pi / 15 and sphere_integral(sp.Integer(1)) == 4 * sp.pi


# ------------------------------------------------------------------------------------------- (2) T2, T3, T4, d root
def chain_root_in_middle(leaf_pot: Tail, ang) -> Dict:
    """Chain leaf - root(t) - leaf, root t^3 K0(t r)/(2 pi^2): branch integral int I0(y) x^2 V V (ang/(2 pi^2)) dx
    at t = -(2 - eps), V the exponential part of the leaf potential's radial factor."""
    H = (leaf_pot * leaf_pot).times(ang / (2 * sp.pi ** 2), shift=2)
    return branch([(1, kern_I(0), H)], 2)


def test_T2_pz_chain_constants():
    """p_z chain: eps^{3/2} S -> A (1 + kappa eps), A = 1/(3 pi^2), kappa = 47/8, from E (the rule) and I0's Hankel
    series, with the angular weight int cos^2 dOmega computed.  Rejects 1/(2 pi^2), 45/8, an exponent other than
    -3/2."""
    P = chain_root_in_middle(f_p().potential(1), axial(MU ** 2))
    lead = min(P)
    assert lead == -sp.Rational(3, 2)
    assert sp.simplify(P[lead] - A_PZ) == 0, P[lead]
    assert sp.simplify(P[lead + 1] / P[lead]) == KAPPA_PZ


def test_T3_s_chain_constants():
    """s chain: A = 1/pi^2, kappa = 31/8, from the s-leaf potential (the rule) with angular weight 4 pi."""
    P = chain_root_in_middle(rho_s().potential(0), axial(SYM_S))
    lead = min(P)
    assert lead == -sp.Rational(3, 2)
    assert sp.simplify(P[lead] - A_S) == 0, P[lead]
    assert sp.simplify(P[lead + 1] / P[lead]) == KAPPA_S


def test_T4_star_of_five():
    """Star of five in the convention of prop:star_attain: eps^{1/2} int_1^oo I0((4 - eps) x) R^4 x^-2 dx ->
    pi^2/(8 sqrt 2), with R = -(pi/2) r x (exponential part of the s-leaf potential, from the rule).  Leading
    constant only.  Rejects pi^2/8, a lead at eps^-1 or eps^0, and R with the wrong normalisation."""
    Rr = rho_s().potential(0).times(-sp.pi / 2, shift=1)
    H = (Rr * Rr * Rr * Rr).times(1, shift=-2)
    P = branch([(1, kern_I(0), H)], 4, order=0)
    assert list(P) == [-sp.Rational(1, 2)]
    assert sp.simplify(P[-sp.Rational(1, 2)] - C5) == 0, P


def end_root_sources() -> Tuple[Tail, Tail]:
    """e^{-2x} part of the middle source of a chain leaf - middle - root(t): (radial) density x leaf potential, for
    the p_z - p_z chain (f_p E, to be multiplied by cos^2) and for the s - s chain."""
    return f_p() * f_p().potential(1), rho_s() * rho_s().potential(0)


def test_d_root_chain_constants():
    """Chain p_z leaf - p_z middle - d root, record P63L2's convention S1 = int I2(y) x^2 (2/3) Phi_2[f_p E] dx:
    eps^{3/2} S1 -> 1/(12 pi), kappa = 15/8 (recomputed here; also pinned in test_paper63_attainment_l2.py).
    Rejects 1/(6 pi), 1/(3 pi^2), the wrong Legendre weight."""
    w, _ = end_root_sources()
    leg2 = axial(MU ** 2 * P2) / axial(P2 * P2)                # coefficient of P2 in cos^2
    assert leg2 == sp.Rational(2, 3)
    P = branch([(leg2, kern_I(2), w.potential(2).times(1, shift=2))], 2)
    lead = min(P)
    assert lead == -sp.Rational(3, 2)
    assert sp.simplify(P[lead] - A_D) == 0, P[lead]
    assert sp.simplify(P[lead + 1] / P[lead]) == KAPPA_D


def test_leading_coefficient_is_radial_constant_times_angular_integral():
    """The paragraph's 'product of nonzero radial constants times the angular integral', on the shapes on record,
    all in one normalisation (root density t^(3+l) K_l/(2 pi^2), full angular integral):
      root in the middle:  A(p_z chain)/A(s chain) = int cos^2/(4 pi) = 1/3;
      root at the end:     A(p_z - p_z - d)/A(s - s - s) = int cos^2 P2/(4 pi) = 2/15, with A(s - s - s) = 1/(4 pi^2)
                           (derived here, not on record) and A(p_z - p_z - d) = (4 pi/5)/(2 pi^2) x 1/(12 pi).
    The first corrections differ (47/8 vs 31/8), so the rule is a statement about the leading coefficient only.
    Rejects an amplitude that is not proportional to the angular integral (e.g. 1/(6 pi) for the d root)."""
    Pp = chain_root_in_middle(f_p().potential(1), axial(MU ** 2))
    Ps = chain_root_in_middle(rho_s().potential(0), axial(SYM_S))
    lead = -sp.Rational(3, 2)
    assert sp.simplify(Pp[lead] / Ps[lead] - axial(SYM_PZ ** 2) / axial(SYM_S)) == 0
    assert sp.simplify(Pp[lead + 1] / Pp[lead] - Ps[lead + 1] / Ps[lead]) != 0
    wp, ws = end_root_sources()
    norm = 1 / (2 * sp.pi ** 2)
    Fs = branch([(norm * axial(SYM_S), kern_I(0), ws.potential(0).times(1, shift=2))], 2)
    Fd = branch([(norm * axial(MU ** 2 * P2), kern_I(2), wp.potential(2).times(1, shift=2))], 2)
    assert sp.simplify(Fs[lead] - 1 / (4 * sp.pi ** 2)) == 0
    assert sp.simplify(Fd[lead] - axial(P2 * P2) * norm * A_D) == 0
    assert sp.simplify(Fd[lead] / Fs[lead] - axial(SYM_PZ ** 2 * SYM_D) / axial(SYM_S)) == 0
    assert sp.simplify(Fd[lead] / Fs[lead]) == sp.Rational(2, 15)


def test_cancelling_root_lead_is_zero_and_next_is_not():
    """Root lap - (5/2) h_2(d): channel weights int cos^2 dOmega (kernel I0 + I1/y) and (5/2) int cos^2 P2 dOmega
    (kernel I2, opposite sign), each 4 pi/3.  The eps^{-3/2} coefficient is exactly zero; the eps^{-1/2} coefficient
    is (4 pi/3) x 3/(8 pi) = 1/2, not zero, so the branch is reached one order lower (the paragraph's 'second
    shape').  Rejects a non-cancelling root (5/2 -> 2: eps^{-3/2} survives)."""
    w, _ = end_root_sources()
    w0, w2 = axial(MU ** 2), ROOT_COEF * axial(MU ** 2 * P2)
    P = branch([(w0, kern_I(0, True), w.potential(0).times(1, shift=2)),
                (-w2, kern_I(2), w.potential(2).times(1, shift=2))], 2, order=3)
    assert -sp.Rational(3, 2) not in P, P
    assert min(P) == -sp.Rational(1, 2)
    assert sp.simplify(P[-sp.Rational(1, 2)] - w0 * A_CANCEL) == 0, P


# ------------------------------------------------------------------------------------------- (3) pure harmonics
def test_pure_harmonic_densities_and_potentials_factorise():
    """For each real harmonic h of degree l <= 2 and a radial f: h(d) f(r) = h(x) (r^-1 d/dr)^l f (Hobson), checked on
    f = e^{-c r} and f = r^p; and lap[phi(r) h(x)] = h(x) [phi'' + (2l+2) phi'/r] on phi = e^{-c r} r^p, so the
    Poisson equation keeps the form radial x h.  Hence in a star of pure harmonics about one centre every factor
    of the integrand is radial x harmonic and the fibre is (a radial integral) x int prod h dOmega.  The fibre
    itself is not computed here.  Rejects a harmonic that is not harmonic (e.g. z^2 in place of (3z^2 - r^2)/2)."""
    c, p = sp.symbols("c p", positive=True)
    r = sp.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    rr = sp.symbols("rr", positive=True)

    def apply(h, f):
        out = sp.Integer(0)
        for (i, j, k), cf in sp.Poly(h, X, Y, Z).terms():
            g = f
            for v, n in ((X, i), (Y, j), (Z, k)):
                if n:
                    g = sp.diff(g, v, n)
            out += cf * g
        return out

    lap = lambda g: sum(sp.diff(g, v, 2) for v in (X, Y, Z))                                # noqa: E731
    pt = {X: sp.Rational(3, 10), Y: sp.Rational(-7, 10), Z: sp.Rational(11, 10)}
    for name, h in HARM.items():
        l = DEG[name]
        assert sp.expand(lap(h)) == 0, name
        for F in (sp.exp(-c * rr), rr ** p):
            rad = F
            for _ in range(l):
                rad = sp.diff(rad, rr) / rr
            diff = apply(h, F.subs(rr, r)) - h * rad.subs(rr, r)
            assert sp.simplify(diff.subs(pt)) == 0, (name, F)
        phi = sp.exp(-c * rr) * rr ** p
        want = h * (sp.diff(phi, rr, 2) + (2 * l + 2) * sp.diff(phi, rr) / rr).subs(rr, r)
        assert sp.simplify((lap(phi.subs(rr, r) * h) - want).subs(pt)) == 0, name
    # the control: z^2 alone is not a pure harmonic and does not factorise
    assert sp.expand(lap(Z ** 2)) != 0


def test_star_angular_integrals_for_pure_harmonics():
    """int_{S^2} prod h dOmega for all triples and a sample of quadruples of the nine real harmonics of degree
    <= 2: zero whenever the total degree is odd or a coordinate appears to an odd power; the axial triples
    (P_l1 P_l2 P_l3) agree with sympy's Gaunt coefficients (nonzero exactly when l1 + l2 + l3 is even and the
    triangle holds).  Named zeros: (z, 1, 1), (z, z, z), (z2, 1, 1), (x, z, 1), (xz, z, z); named nonzeros:
    (z, z, 1) = 4 pi/3, (z, z, z2) = 8 pi/15, (z2, z2, z2) = 8 pi/35.  With the factorisation above these are the
    stars whose fibre vanishes identically, or does not vanish for that reason."""
    val = lambda *n: sphere_integral(sp.prod([HARM[k] for k in n]))                          # noqa: E731
    for names in itertools.combinations_with_replacement(HARM, 3):
        v = val(*names)
        if sum(DEG[k] for k in names) % 2:
            assert v == 0, names
    for names in (("z", "1", "1"), ("z", "z", "z"), ("z2", "1", "1"), ("x", "z", "1"), ("xz", "z", "z"),
                  ("x2-y2", "z", "z"), ("z", "z", "z2", "z"), ("xy", "x", "z", "1")):
        assert val(*names) == 0, names
    assert val("z", "z", "1") == 4 * sp.pi / 3
    assert val("z", "z", "z2") == 8 * sp.pi / 15
    assert val("z2", "z2", "z2") == 8 * sp.pi / 35
    assert val("xz", "x", "z") == 4 * sp.pi / 15
    assert val("x2-y2", "x", "x") == 8 * sp.pi / 15
    assert val("z", "z", "z2", "z2") != 0
    leg = {0: sp.Integer(1), 1: MU, 2: P2}
    for l1, l2, l3 in itertools.product(range(3), repeat=3):
        mine = axial(leg[l1] * leg[l2] * leg[l3])
        ref = gaunt(l1, l2, l3, 0, 0, 0) * sp.sqrt(sp.Rational(64, (2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1))
                                                     * sp.pi ** 3)
        assert sp.simplify(mine - ref) == 0, (l1, l2, l3)
        allowed = (l1 + l2 + l3) % 2 == 0 and abs(l1 - l2) <= l3 <= l1 + l2
        assert (mine != 0) == allowed, (l1, l2, l3)


# ------------------------------------------------------------------------------------------- paper sentences
PHRASES = [
    # the conjecture's hypothesis and its s-type case
    r"\int_{S^2}\prod_v\sigma_v(-\omega)\,d\omega\neq0",
    r"($\partial\to-\omega$, $\nabla^2\to1$)",
    r"In particular every $s$-type tree, where $\sigma_v\equiv1$, reaches rung $n$.",
    r"\textbf{[ARGUMENT]} This is not proved.",
    # (1)
    r"For a source $e^{-kr}r^{a}\sum_js_jr^{-j}$ in angular channel $L$, the exponential part of its Coulomb "
    r"potential is the formal solution of $(\theta-L)(\theta+L+1)\Phi=-4\pi r^2s$, whose leading term "
    r"$-4\pi s_0/k^2$ does not depend on $L$.",
    r"leading coefficient a product of nonzero radial constants times the angular integral above",
    # (2)
    r"\textbf{[MEASURED]} The local rule reproduces exactly the leading constants already on record (the $p_z$ chain's "
    r"$1/(3\pi^2)$, the star of five), and it predicted the amplitude $1/(12\pi)$ of the chain with a $d$ root "
    r"before that was computed",
    r"with amplitude $1/(3\pi^2)$, one third of the $s$ chain's",
    r"$c_5\Gamma(\tfrac12)=\pi^2/(8\sqrt2)",
    r"$d$ root gives $\varepsilon^{-3/2}$ with amplitude $1/(12\pi)$ and first correction "
    r"$(1+\tfrac{15}{8}\varepsilon)$",
    # (4)
    r"The second shape there lies outside the conjecture:\ its angular integral vanishes, the leading term "
    r"cancels, and the branch is reached one order lower.",
    r"$(\nabla^2-\tfrac52h_2(\partial))\,t^3K_0(tr)/(2\pi^2)$, whose angular weight against $\cos^2\theta$ vanishes",
    # (3)
    r"For a star of pure harmonics a vanishing angular integral makes the fibre vanish identically.",
    r"These records read a branch coefficient;\ none computes an operator",
]
RECORDS = ["2026-11-13_star5_attainment.md", "2026-11-16_paper63_fibre_l1.md", "2026-11-17_paper63_attainment_l2.md"]


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text)


def missing_phrases(text: str) -> List[str]:
    flat = _flat(text)
    return [p for p in PHRASES if _flat(p) not in flat]


def test_paper_sentences_present():
    """The sentences this file backs are in the paper as written (whitespace-insensitive), inside or after the
    conjecture for the paragraph's own sentences, and the three records the [MEASURED] sentence points to exist.
    Rejects a reworded or strengthened sentence (e.g. '[ARGUMENT] This is not proved' removed, 1/(12 pi) changed).
    The 'before that was computed' in the paper is a statement about the order of events; the record files exist
    but this test cannot check dates."""
    text = TEX.read_text(encoding="utf-8")
    assert missing_phrases(text) == []
    flat = _flat(text)
    start = flat.index(r"\label{conj:tree_attain}")
    end = flat.index(r"\textbf{[PRIOR ART]}", start)
    para = flat[start:end]
    for p in PHRASES[:7] + PHRASES[10:11] + PHRASES[12:]:
        assert _flat(p) in para, p
    for rec in RECORDS:
        assert (ROOT / "docs" / "predictions" / rec).exists(), rec


# TIMINGS (2026-11-18, this machine): 39 tests, 19 s in all; 12 s of it is the direct quadrature in
# test_T1_closed_form_is_the_true_potential, 3 s the Hobson check; every other test is below 1 s.
#
# FIRE TESTS (2026-11-18, a separate step after the guard was written).  22 source plants, each run on a scratch
# copy of this file so that every failing test is listed; 7 of them repeated with debug/qa/fire_test.py (all FIRED).
# All 22 fired.  Which tests caught each:
#   lead -4 pi s0/k^2 -> /k            true potential, by parts, radial equation, d root, cancelling root, symbolic
#                                      (T1-T4 do NOT see it: their leaves have k = 1)
#   -4 pi -> +4 pi                     true potential, by parts, radial equation, T1, T1b, symbolic
#                                      (T2, T3, d root do NOT: the sign enters squared)
#   L at order r^-1 (b+1 -> b+1+L)     true potential, by parts, radial equation, T1, T2, d root, cancelling root
#   L(L+1) sign; L(L+1) removed        true potential, by parts, radial equation, T1, symbolic
#   47/8 -> 45/8; 1/(3pi^2) -> 1/(2pi^2)          T2
#   1/(12pi) -> 1/(6pi)                d root, radial-constant-times-angular-integral
#   15/8 -> 13/8                       d root
#   1/pi^2 -> 2/pi^2; 31/8 -> 33/8     T3
#   pi^2/(8 sqrt 2) -> pi^2/8          T4
#   tail 11/8 -> 13/8; -463 -> +463    T1
#   root 5/2 -> 2                      angular integrals, cancelling root (eps^-3/2 survives)
#   3/(8pi) -> 1/(4pi)                 cancelling root
#   true-potential limit k^4 -> k^3    channels-differ test (3 of 4 cases; the k = 1 case cannot tell)
#   true potential: bracket sign       true potential, channels-differ
#   z2 -> z^2 (not harmonic)           Hobson, both angular tests
#   d symbol P2 -> cos                 angular integrals, radial-constant-times-angular-integral
#   I_nu without (-1)^m                Hankel vs mpmath, T2, T3, d root, cancelling root
#   by-parts 1/k^(j+1) -> 1/k^j        by parts
# Paper plants, on a temporary copy of the .tex (TEX re-pointed in a driver; the paper untouched), 9 of 9 fired:
#   each of the four backed sentences deleted or reworded; 1/(12 pi) -> 1/(6 pi) in the paragraph; "[ARGUMENT] This
#   is not proved." replaced by a proof tag; "none computes an operator" reversed; 5/2 -> 2 in the root; "neq 0"
#   dropped from the conjecture.
# Known blind spot, by design: test_leading_term_and_where_L_enters_symbolic restates the recursion and would accept
# any rule of that shape; the tests named 'true potential', 'by parts' and 'radial equation' carry claim (1).
