"""QFD validation references -- INDEPENDENT numerics for every closed form.

Nothing here is on the production path; this module exists only to falsify
`qfd_core` / `geovac.two_center_eri`.  Three kinds of reference are provided:

  (1) mpmath quadrature.  One-electron integrals are integrated raw in prolate
      spheroidal coordinates (no auxiliary-function expansion); the one-center
      and (AA|BB) and hybrid two-electron classes are integrated through the
      Newton spherical-potential route, which shares no code with the
      shell-kernel / multipole machinery under test.

  (2) LITERATURE closed forms for the H2 1s two-center integrals -- Sugiura's
      exchange integral and the classical Coulomb and hybrid formulas.  These are
      external to the corpus entirely and are the strongest available check.

  (3) a fully numeric Neumann evaluation of the exchange class (both the eta
      half and the ordered xi half by quadrature), which isolates
      `ordered_xi_general` + `integrate_poly_exp` from the rest of the assembly.
"""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp
from mpmath import mp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import qfd_core as Q  # noqa: E402


# ------------------------------------------------------------------- radial fns

def radial_fn(orb):
    """Numeric R_{n0}(r) for an orbital spec, plus its decay rate (mpf)."""
    coeffs, a = Q.s_radial(orb[1], orb[2])
    cs = [(k, mp.mpf(str(sp.N(v, mp.dps + 10)))) for k, v in coeffs]
    am = mp.mpf(str(sp.N(a, mp.dps + 10)))

    def f(r):
        return mp.fsum(c * r ** k for k, c in cs) * mp.e ** (-am * r)
    return f, am


# ------------------------------------------------ one-electron: raw 2D spheroidal

def one_electron_quad(oi, oj, kernel, R, ZA=1, ZB=1):
    """<i|K|j> by raw quadrature in prolate spheroidal coordinates.

    kernel in {"S", "inv_ra", "inv_rb", "T"}.  For "T" the Laplacian is applied
    analytically to the ket (a two-line derivative of r^j e^{-b r}), so the
    quadrature is of a smooth explicit function -- no numerical differentiation.
    """
    Rm = mp.mpf(str(sp.N(R, mp.dps + 10)))
    fi, _ai = radial_fn(oi)
    if kernel == "T":
        coeffs, b = Q.s_radial(oj[1], oj[2])
        lap = Q._laplacian_terms(coeffs, b)
        bm = mp.mpf(str(sp.N(b, mp.dps + 10)))
        lm = [(k, mp.mpf(str(sp.N(v, mp.dps + 10)))) for k, v in lap.items()]

        def fj(r):
            return mp.fsum(c * r ** k for k, c in lm) * mp.e ** (-bm * r)
    else:
        fj, _b = radial_fn(oj)

    def integrand(xi, eta):
        ra = Rm * (xi + eta) / 2
        rb = Rm * (xi - eta) / 2
        ri = ra if oi[0] == "A" else rb
        rj = ra if oj[0] == "A" else rb
        val = (xi ** 2 - eta ** 2) * fi(ri) * fj(rj)
        if kernel == "inv_ra":
            val /= ra
        elif kernel == "inv_rb":
            val /= rb
        return val

    pref = Rm ** 3 / 16                       # 2 pi (R^3/8) / (4 pi)
    if kernel == "T":
        pref *= mp.mpf(-1) / 2
    val = mp.quad(lambda xi: mp.quad(lambda eta: integrand(xi, eta),
                                     [-1, 0, 1]),
                  [1, 2, 4, 8, 16, mp.inf])
    return pref * val


def one_electron_quad_same_center(oi, oj, which, R):
    """<i|1/r_C|j> for i, j on the SAME center, C the other one -- Newton's
    theorem, a 1D quadrature that shares nothing with the spheroidal route."""
    Rm = mp.mpf(str(sp.N(R, mp.dps + 10)))
    fi, _ = radial_fn(oi)
    fj, _ = radial_fn(oj)
    assert oi[0] == oj[0] and oi[0] != which
    inner = mp.quad(lambda r: fi(r) * fj(r) * r ** 2, [0, Rm]) / Rm
    outer = mp.quad(lambda r: fi(r) * fj(r) * r, [Rm, 2 * Rm, 4 * Rm, mp.inf])
    return inner + outer


# ---------------------------------------------------------- spherical potentials

def spherical_potential(oa, ob):
    """V(r) of the one-center density R_a R_b Y00^2, by QUADRATURE:
    V(r) = (1/r) int_0^r P s^2 ds + int_r^oo P s ds  with P = R_a R_b."""
    fa, _ = radial_fn(oa)
    fb, _ = radial_fn(ob)

    def V(r):
        lo = mp.quad(lambda s: fa(s) * fb(s) * s ** 2, [0, r])
        hi = mp.quad(lambda s: fa(s) * fb(s) * s, [r, r + 2, r + 8, mp.inf])
        return lo / r + hi
    return V


def spherical_potential_closed(oa, ob):
    """Same V(r), Newton's theorem in closed form (elementary finite sums).

    Used to keep the (AA|BB) / hybrid references at 2D instead of 4D quadrature.
    It is Newton's shell theorem for a spherical density -- a different object
    from the shell-kernel / multipole machinery under test in
    `aabb_closed_form` / `hybrid_closed_form` -- and it is spot-checked against
    the pure-quadrature `spherical_potential` above.
    """
    ca, aa = Q.s_radial(oa[1], oa[2])
    cb, ab = Q.s_radial(ob[1], ob[2])
    P = Q._prod_coeffs(ca, cb)
    b = mp.mpf(str(sp.N(aa + ab, mp.dps + 10)))
    cs = [(k, mp.mpf(str(sp.N(v, mp.dps + 10)))) for k, v in P]

    def _lower(qq, x):        # int_0^x s^qq e^{-b s} ds
        return (mp.factorial(qq) / b ** (qq + 1)
                * (1 - mp.e ** (-b * x)
                   * mp.fsum((b * x) ** k / mp.factorial(k)
                             for k in range(qq + 1))))

    def _upper(qq, x):        # int_x^oo s^qq e^{-b s} ds
        return (mp.factorial(qq) / b ** (qq + 1) * mp.e ** (-b * x)
                * mp.fsum((b * x) ** k / mp.factorial(k)
                          for k in range(qq + 1)))

    def V(r):
        return mp.fsum(c * (_lower(k + 2, r) / r + _upper(k + 1, r))
                       for k, c in cs)
    return V


def one_center_eri_quad(oa, ob, oc, od):
    """One-center (ab|cd) by quadrature against the QUADRATURE potential --
    fully numeric, no closed form anywhere."""
    fa, _ = radial_fn(oa)
    fb, _ = radial_fn(ob)
    V = spherical_potential(oc, od)
    return mp.quad(lambda r: fa(r) * fb(r) * r ** 2 * V(r),
                   [0, 1, 2, 4, 8, 16, mp.inf])


def aabb_quad(oa, ob, oc, od, R):
    """(ab|cd), a,b on one center and c,d on the other: spherical density x
    spherical potential, bipolar 2D quadrature (Newton route)."""
    Rm = mp.mpf(str(sp.N(R, mp.dps + 10)))
    V = spherical_potential_closed(oa, ob)
    fc, _ = radial_fn(oc)
    fd, _ = radial_fn(od)

    def outer(r):
        # 2 pi int_-1^1 V(|r_vec - R zhat|) du  with u = cos angle
        return (2 * mp.pi * fc(r) * fd(r) * r ** 2
                * mp.quad(lambda u: V(mp.sqrt(r ** 2 + Rm ** 2 - 2 * r * Rm * u)),
                          [-1, 0, 1]))
    # (1/4pi) from the two Y00 factors of the (cd) pair
    return mp.quad(outer, [0, Rm, 2 * Rm, 4 * Rm, mp.inf]) / (4 * mp.pi)


def hybrid_quad(oa, ob, oc, od, R):
    """(ab|cd) with (a,b) a one-center pair and (c,d) a two-center pair.

    V_ab is the exact spherical potential of the one-center density; the outer
    integral is a raw 2D prolate-spheroidal quadrature.  Independent of the
    shell-kernel machinery in `hybrid_closed_form`.
    """
    Rm = mp.mpf(str(sp.N(R, mp.dps + 10)))
    assert oa[0] == ob[0]
    V = spherical_potential_closed(oa, ob)
    fc, _ = radial_fn(oc)
    fd, _ = radial_fn(od)

    def integrand(xi, eta):
        ra = Rm * (xi + eta) / 2
        rb = Rm * (xi - eta) / 2
        rv = ra if oa[0] == "A" else rb
        rc = ra if oc[0] == "A" else rb
        rd = ra if od[0] == "A" else rb
        return (xi ** 2 - eta ** 2) * V(rv) * fc(rc) * fd(rd)

    val = mp.quad(lambda xi: mp.quad(lambda eta: integrand(xi, eta), [-1, 0, 1]),
                  [1, 2, 4, 8, mp.inf])
    return Rm ** 3 / 16 * val


# --------------------------------------------------- exchange: numeric Neumann

def exchange_numeric_neumann(ZA, oa, ob, ZB, oc, od, R, tau_max=10):
    """Exchange with BOTH halves by quadrature -- isolates `ordered_xi_general`
    and `integrate_poly_exp` from the rest of the assembly."""
    from geovac.two_center_eri import (_deriv, _legendre_Q, _z, eta_s,
                                       two_center_spheroidal_product, xi_s)
    s = abs(oa[2] - ob[2])
    Rs = sp.nsimplify(R)
    P1, h1, p1, q1 = two_center_spheroidal_product(ZA, oa, ZB, ob, Rs)
    P2, h2, p2, q2 = two_center_spheroidal_product(ZA, oc, ZB, od, Rs)
    H1 = int(h1 + sp.Rational(s, 2))
    H2 = int(h2 + sp.Rational(s, 2))

    def cof(P):
        return dict(sp.Poly(sp.expand((xi_s ** 2 - eta_s ** 2) * P),
                            xi_s, eta_s).terms())

    c1, c2 = cof(P1), cof(P2)
    F = lambda e: mp.mpf(str(sp.N(e, mp.dps + 10)))  # noqa: E731
    p1m, p2m, q1m, q2m = F(p1), F(p2), F(q1), F(q2)
    Rm = F(Rs)
    total = mp.mpf(0)
    per_tau = []
    for tau in range(s, tau_max + 1):
        Dp = sp.lambdify(_z, sp.expand(_deriv(sp.legendre(tau, _z), s)), "mpmath")
        Dq = sp.lambdify(_z, sp.expand(_deriv(_legendre_Q(tau), s)), "mpmath")

        def B(k, H, q, _Dp=Dp):
            return mp.quad(lambda e: e ** k * (1 - e ** 2) ** H * _Dp(e)
                           * mp.e ** (-q * e), [-1, 0, 1])

        B1 = {k: B(k, H1, q1m) for k in {k for _j, k in c1}}
        B2 = {k: B(k, H2, q2m) for k in {k for _j, k in c2}}
        xicache: dict = {}

        def Xi(j1, j2, _Dp=Dp, _Dq=Dq):
            if (j1, j2) in xicache:
                return xicache[(j1, j2)]

            def outer(x1):
                lo = mp.quad(lambda x2: x2 ** j2 * (x2 ** 2 - 1) ** H2
                             * mp.e ** (-p2m * x2) * _Dp(x2), [1, x1])
                hi = mp.quad(lambda x2: x2 ** j2 * (x2 ** 2 - 1) ** H2
                             * mp.e ** (-p2m * x2) * _Dq(x2),
                             [x1, x1 + 2, x1 + 8, mp.inf])
                f1 = x1 ** j1 * (x1 ** 2 - 1) ** H1 * mp.e ** (-p1m * x1)
                return f1 * (_Dq(x1) * lo + _Dp(x1) * hi)
            v = mp.quad(outer, [1, 2, 4, 8, mp.inf])
            xicache[(j1, j2)] = v
            return v

        acc = mp.mpf(0)
        for (j1, k1), cc1 in c1.items():
            if B1[k1] == 0:
                continue
            for (j2, k2), cc2 in c2.items():
                if B2[k2] == 0:
                    continue
                acc += F(cc1) * F(cc2) * B1[k1] * B2[k2] * Xi(j1, j2)
        w = ((-1) ** (oa[2] - ob[2]) * (2 * tau + 1)
             * mp.mpf(str(sp.N(sp.factorial(tau - s) / sp.factorial(tau + s),
                               mp.dps + 10))) ** 2)
        per_tau.append(w * acc)
        total += w * acc
    C = (Rm ** 3 / 8) ** 2 * (2 * mp.pi) ** 2 * (2 / Rm)
    return C * total, [C * t for t in per_tau]


# ------------------------------------------------- literature closed forms (H2 1s)

def lit_coulomb(zeta, R):
    w = zeta * R
    return zeta * (1 / w - mp.e ** (-2 * w)
                   * (1 / w + mp.mpf(11) / 8 + 3 * w / 4 + w ** 2 / 6))


def lit_hybrid(zeta, R):
    w = zeta * R
    return zeta * (mp.e ** (-w) * (w + mp.mpf(1) / 8 + 5 / (16 * w))
                   - mp.e ** (-3 * w) * (mp.mpf(1) / 8 + 5 / (16 * w)))


def lit_exchange_sugiura(zeta, R):
    """Sugiura (1927): the two-center 1s exchange integral, in closed form.

    Content {exp, ln, gamma, Ei} -- the same {E_1, ln, gamma} seed set the corpus
    derives for the exchange class, arrived at by a completely different route.
    """
    w = zeta * R
    S = mp.e ** (-w) * (1 + w + w ** 2 / 3)
    Sm = mp.e ** (w) * (1 - w + w ** 2 / 3)
    t1 = -mp.e ** (-2 * w) * (mp.mpf(-25) / 8 + 23 * w / 4 + 3 * w ** 2
                              + w ** 3 / 3)
    t2 = (6 / w) * ((mp.euler + mp.log(w)) * S ** 2
                    + Sm ** 2 * mp.ei(-4 * w) - 2 * S * Sm * mp.ei(-2 * w))
    return zeta * (t1 + t2) / 5


def lit_overlap(zeta, R):
    w = zeta * R
    return mp.e ** (-w) * (1 + w + w ** 2 / 3)


def lit_kinetic_cross(zeta, R):
    """<1s_A|-1/2 nabla^2|1s_B> = (zeta^2/2)(1 + w - w^2/3) e^{-w}.

    From nabla^2 e^{-zeta r} = (zeta^2 - 2 zeta/r) e^{-zeta r}, so
    T_AB = -(zeta^2/2) S_AB + zeta <A|1/r_B|B>.
    """
    w = zeta * R
    return zeta ** 2 / 2 * (1 + w - w ** 2 / 3) * mp.e ** (-w)


def lit_nuc_same(zeta, R):
    """<1s_A|1/r_B|1s_A> = (1/R)[1 - (1 + w) e^{-2w}]."""
    w = zeta * R
    return (1 - (1 + w) * mp.e ** (-2 * w)) / R


def lit_nuc_cross(zeta, R):
    """<1s_A|1/r_A|1s_B> = zeta (1 + w) e^{-w}."""
    w = zeta * R
    return zeta * (1 + w) * mp.e ** (-w)
