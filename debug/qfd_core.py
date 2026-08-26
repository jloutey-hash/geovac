"""QFD -- Quadrature-Free Diatomic: the closed-form integral core.

Every production-path integral here is a CLOSED FORM: an exact sympy expression
built from rationals, sqrt, pi, exp, expint (E_1), log and EulerGamma. Nothing
in the production path calls a quadrature routine, so the value can be evaluated
to any requested precision by `sp.N(expr, dps)`.

Scope of this module: s-type (l = 0) hydrogenic orbitals, decay rate a = Z/n
(the `geovac.two_center_eri` convention -- HYDROGENIC, not Coulomb-Sturmian, not
Gaussian; see memory/polyatomic_state_of_play.md).

WHAT IS BUILT HERE (one-electron), and why it had to be:
  the two-center ERI engine in geovac/two_center_eri.py covers the two-ELECTRON
  classes only.  The one-electron closed forms (overlap / kinetic / both
  nuclear-attraction kernels) are derived here from the Mulliken auxiliary
  integrals

      A_m(p) = int_1^oo  xi^m  e^{-p xi} dxi        (finite, all-positive sum)
      B_n(q) = int_-1^1  eta^n e^{-q eta} deta      (finite, three-term recurrence)

  in prolate spheroidal coordinates (foci on the nuclei):

      r_A = R(xi+eta)/2 ,  r_B = R(xi-eta)/2 ,
      d3r = (R^3/8)(xi^2-eta^2) dxi deta dphi

  The classical cancellation (xi^2-eta^2) = (xi+eta)(xi-eta) is what makes the
  1/r_A and 1/r_B kernels polynomial, hence closed-form:

      I(i,j) = int r_A^i r_B^j e^{-alpha r_A - beta r_B} d3r
             = 2 pi (R/2)^{3+i+j} sum_{u,v} C(i+1,u) C(j+1,v) (-1)^v
                                            A_{i+j+2-u-v}(p) B_{u+v}(q)
      p = (alpha+beta)R/2 ,  q = (alpha-beta)R/2 ,  valid for i, j >= -1.

  Kinetic energy uses the radial Laplacian on s functions in closed form,
      nabla^2 [r^j e^{-b r}] = [ j(j+1) r^{j-2} - 2b(j+1) r^{j-1} + b^2 r^j ] e^{-b r}
  whose lowest surviving power is r^{-1} (the j = 0 term has coefficient 0), so
  it stays inside I(i, j >= -1).  No quadrature, no numerical differentiation.

WHAT IS REUSED (two-electron): geovac.two_center_eri
  (AA|BB)  aabb_closed_form   -- elementary, {exp}
  hybrid   hybrid_closed_form -- {exp, E_1, log}   (s-type pair: {exp} only)
  (AB|AB)  assembled HERE at arbitrary precision from the same symbolic pieces
           `exchange_value` uses (two_center_spheroidal_product,
           integrate_poly_exp, ordered_xi_general).  The library entry point
           casts to float, which caps certification at ~16 digits, so the loop is
           re-run keeping every factor symbolic.

Diagnostic / build module.  No paper, CLAUDE.md or tests edits.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path

import sympy as sp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.two_center_eri import (  # noqa: E402
    _deriv, _z, aabb_closed_form, eta_s, exchange_value, hybrid_closed_form,
    integrate_poly_exp, ordered_xi_general, radial_norm, radial_poly,
    two_center_spheroidal_product, xi_s,
)

# ---------------------------------------------------------------- radial data


@lru_cache(maxsize=None)
def s_radial(Z, n):
    """Normalized hydrogenic R_{n0}: ((power, exact coeff) tuple, rate a = Z/n)."""
    Zf = Z if isinstance(Z, Fraction) else Fraction(Z)
    c, a = radial_poly(Zf, n, 0)
    N = radial_norm(Zf, n, 0)
    return tuple(sorted((k, sp.nsimplify(N * v)) for k, v in c.items())), a


def _prod_coeffs(cA, cB):
    """Coefficients of the product of two radial polynomial lists."""
    out: dict = {}
    for k1, v1 in cA:
        for k2, v2 in cB:
            out[k1 + k2] = out.get(k1 + k2, 0) + v1 * v2
    return tuple(sorted(out.items()))


# ------------------------------------------------- Mulliken auxiliary integrals


@lru_cache(maxsize=None)
def A_aux(m: int, p):
    """A_m(p) = int_1^oo xi^m e^{-p xi} dxi, exact (all-positive finite sum)."""
    assert m >= 0
    acc = sp.Integer(0)
    term = sp.Integer(1)                       # m!/(m-i)!
    for i in range(m + 1):
        acc += term / p ** (i + 1)
        term *= (m - i)
    return sp.exp(-p) * acc


@lru_cache(maxsize=None)
def B_aux(n: int, q):
    """B_n(q) = int_-1^1 eta^n e^{-q eta} deta, exact."""
    assert n >= 0
    if q == 0:
        return sp.Integer(0) if n % 2 else sp.Rational(2, n + 1)
    if n == 0:
        return (sp.exp(q) - sp.exp(-q)) / q
    return (((-1) ** n * sp.exp(q) - sp.exp(-q)) / q
            + sp.Integer(n) * B_aux(n - 1, q) / q)


def I2c(i: int, j: int, alpha, beta, R):
    """int r_A^i r_B^j e^{-alpha r_A - beta r_B} d3r, exact.  Needs i, j >= -1."""
    assert i >= -1 and j >= -1, f"powers ({i},{j}) below the -1 floor"
    p = (alpha + beta) * R / 2
    q = (alpha - beta) * R / 2
    ii, jj = i + 1, j + 1
    tot = sp.Integer(0)
    for u in range(ii + 1):
        cu = comb(ii, u)
        for v in range(jj + 1):
            tot += (cu * comb(jj, v) * (-1) ** v
                    * A_aux(ii - u + jj - v, p) * B_aux(u + v, q))
    return 2 * sp.pi * (R / 2) ** (3 + i + j) * tot


def _gam(m: int, c):
    """int_0^oo r^m e^{-c r} dr = m!/c^{m+1}, exact."""
    return sp.factorial(m) / c ** (m + 1)


# ------------------------------------------------- one-electron matrix elements
#
# An "orbital" is (center, Z, n) with center in {"A", "B"}; A sits at the origin
# and B at R zhat.  All l = m = 0, so each function carries Y_00 = 1/sqrt(4 pi)
# and the pair carries 1/(4 pi), which is exactly the phi/eta measure I2c already
# integrates over.


def overlap(oi, oj, R):
    """<i|j>, closed form."""
    ci, ai = s_radial(oi[1], oi[2])
    cj, aj = s_radial(oj[1], oj[2])
    if oi[0] == oj[0]:                                   # same center
        return sp.expand(sum(v * _gam(k + 2, ai + aj)
                             for k, v in _prod_coeffs(ci, cj)))
    if oi[0] == "A":
        dA, aA, dB, aB = ci, ai, cj, aj
    else:
        dA, aA, dB, aB = cj, aj, ci, ai
    tot = sp.Integer(0)
    for k1, v1 in dA:
        for k2, v2 in dB:
            tot += v1 * v2 * I2c(k1, k2, aA, aB, R)
    return sp.expand(tot / (4 * sp.pi))


def _inv_r(oi, oj, which: str, R):
    """<i| 1/r_C |j> with C = `which` in {"A", "B"}, closed form."""
    ci, ai = s_radial(oi[1], oi[2])
    cj, aj = s_radial(oj[1], oj[2])
    if oi[0] == oj[0]:
        pc = _prod_coeffs(ci, cj)
        a = ai + aj
        if oi[0] == which:                               # genuine one-center
            return sp.expand(sum(v * _gam(k + 1, a) for k, v in pc))
        # density on one center, kernel on the other -> two-center, other rate 0
        tot = sp.Integer(0)
        for k, v in pc:
            tot += v * (I2c(k, -1, a, sp.Integer(0), R) if oi[0] == "A"
                        else I2c(-1, k, sp.Integer(0), a, R))
        return sp.expand(tot / (4 * sp.pi))
    if oi[0] == "A":
        dA, aA, dB, aB = ci, ai, cj, aj
    else:
        dA, aA, dB, aB = cj, aj, ci, ai
    sh = (-1, 0) if which == "A" else (0, -1)
    tot = sp.Integer(0)
    for k1, v1 in dA:
        for k2, v2 in dB:
            tot += v1 * v2 * I2c(k1 + sh[0], k2 + sh[1], aA, aB, R)
    return sp.expand(tot / (4 * sp.pi))


def _laplacian_terms(coeffs, b):
    """nabla^2 applied to sum_j c_j r^j e^{-b r}: {power: coeff}."""
    lap: dict = {}
    for j, v in coeffs:
        for dk, dc in ((j - 2, sp.Integer(j * (j + 1))),
                       (j - 1, -2 * b * (j + 1)),
                       (j, b ** 2)):
            if dc == 0:
                continue
            lap[dk] = lap.get(dk, 0) + v * dc
    return {k: v for k, v in lap.items() if v != 0}


def kinetic(oi, oj, R):
    """<i| -1/2 nabla^2 |j>, closed form via the s-state radial Laplacian."""
    ci, ai = s_radial(oi[1], oi[2])
    cj, b = s_radial(oj[1], oj[2])
    lap = _laplacian_terms(cj, b)
    if oi[0] == oj[0]:
        tot = sp.Integer(0)
        for k1, v1 in ci:
            for k2, v2 in lap.items():
                tot += v1 * v2 * _gam(k1 + k2 + 2, ai + b)
        return sp.expand(-tot / 2)
    tot = sp.Integer(0)
    for k1, v1 in ci:
        for k2, v2 in lap.items():
            if oi[0] == "A":
                tot += v1 * v2 * I2c(k1, k2, ai, b, R)
            else:
                tot += v1 * v2 * I2c(k2, k1, b, ai, R)
    return sp.expand(-tot / (8 * sp.pi))


def h_core(oi, oj, ZA, ZB, R):
    """<i| -1/2 nabla^2 - Z_A/r_A - Z_B/r_B |j>, closed form (Laplacian route)."""
    return sp.expand(kinetic(oi, oj, R)
                     - ZA * _inv_r(oi, oj, "A", R)
                     - ZB * _inv_r(oi, oj, "B", R))


def h_core_eigentrick(oi, oj, ZA, ZB, R):
    """Independent closed form: h_ij = E_j S_ij - Z_other <i|1/r_other|j>.

    Exact because every basis function here is a hydrogenic eigenfunction of its
    own center: (-1/2 nabla^2 - Z_j/r_j) chi_j = E_j chi_j, E_j = -Z_j^2/(2n_j^2).
    Shares no code with `h_core`'s explicit-Laplacian route, so agreement is a
    genuine cross-check of the kinetic-energy closed form.
    """
    Zj = sp.nsimplify(oj[1])
    Ej = -Zj ** 2 / (2 * oj[2] ** 2)
    other = "B" if oj[0] == "A" else "A"
    Zo = ZB if other == "B" else ZA
    return sp.expand(Ej * overlap(oi, oj, R) - Zo * _inv_r(oi, oj, other, R))


# ------------------------------------------------- one-center two-electron class


def one_center_eri(oa, ob, oc, od):
    """(ab|cd) with all four s functions on the SAME center.

    Only L = 0 survives for s x s densities, so this is the R^0 Slater integral

        int dr1 r1^2 P_ab(r1) [ (1/r1) int_0^r1 r2^2 P_cd + int_r1^oo r2 P_cd ]

    which is elementary (incomplete gammas with integer arguments).
    """
    ca, aa = s_radial(oa[1], oa[2])
    cb, ab = s_radial(ob[1], ob[2])
    cc, ac = s_radial(oc[1], oc[2])
    cd, ad = s_radial(od[1], od[2])
    P1, b1 = _prod_coeffs(ca, cb), aa + ab
    P2, b2 = _prod_coeffs(cc, cd), ac + ad
    r = sp.Symbol("r_qfd", positive=True)

    def lower(qq, bb, x):                      # int_0^x s^qq e^{-bb s} ds
        return (sp.factorial(qq) / bb ** (qq + 1)
                * (1 - sp.exp(-bb * x)
                   * sum((bb * x) ** k / sp.factorial(k) for k in range(qq + 1))))

    def upper(qq, bb, x):                      # int_x^oo s^qq e^{-bb s} ds
        return (sp.factorial(qq) / bb ** (qq + 1) * sp.exp(-bb * x)
                * sum((bb * x) ** k / sp.factorial(k) for k in range(qq + 1)))

    V = sum(v * (lower(k + 2, b2, r) / r + upper(k + 1, b2, r)) for k, v in P2)
    integrand = sp.expand(
        sum(v * r ** (k + 2) * sp.exp(-b1 * r) for k, v in P1) * V)

    tot = sp.Integer(0)
    for term in sp.Add.make_args(integrand):
        coeff, power, decay = sp.Integer(1), 0, sp.Integer(0)
        for f in sp.Mul.make_args(term):
            if f == r:
                power += 1
            elif f.is_Pow and f.base == r:
                power += int(f.exp)
            elif isinstance(f, sp.exp):
                decay += -sp.diff(f.args[0], r)
            else:
                coeff *= f
        assert not coeff.has(r), term
        assert power >= 0, f"negative power {power} reached the one-center ERI"
        tot += coeff * _gam(power, decay)
    return sp.expand(tot)


# ------------------------------------------- exchange class, arbitrary precision


def exchange_closed_form(ZA, oa, ob, ZB, oc, od, R, tau_max: int = 10,
                         return_terms: bool = False):
    """(ab|cd) exchange -- a & c on A, b & d on B -- EXACT symbolic, per tau.

    Same assembly as `geovac.two_center_eri.exchange_value(..., exact_xi=True)`
    with every float() cast removed.  Returns the summed sympy expression, or the
    per-tau list when `return_terms` (used for the truncation-tail audit).

    The tau sum TERMINATES exactly when the two centers of a density carry the
    same orbital exponent (q = (alpha-beta)R/2 = 0, Phase 0-e criterion) -- true
    for homonuclear H2 at zeta_A = zeta_B, false for LiH.
    """
    sig1, sig2 = oa[2] - ob[2], od[2] - oc[2]
    if sig1 != sig2:
        return [] if return_terms else sp.Integer(0)
    s = abs(sig1)
    R = sp.nsimplify(R)

    P1, h1, p1, q1 = two_center_spheroidal_product(ZA, oa, ZB, ob, R)
    P2, h2, p2, q2 = two_center_spheroidal_product(ZA, oc, ZB, od, R)
    H1, H2 = h1 + sp.Rational(s, 2), h2 + sp.Rational(s, 2)
    assert H1.is_integer and H2.is_integer, f"half powers {H1} {H2}"
    H1, H2 = int(H1), int(H2)

    def coeffs_of(P):
        return dict(sp.Poly(sp.expand((xi_s ** 2 - eta_s ** 2) * P),
                            xi_s, eta_s).terms())

    c1, c2 = coeffs_of(P1), coeffs_of(P2)
    C = (R ** 3 / 8) ** 2 * (2 * sp.pi) ** 2 * (2 / R)

    terms = []
    for tau in range(s, tau_max + 1):
        Dp = sp.expand(_deriv(sp.legendre(tau, _z), s))

        def beta(k, H, q, _Dp=Dp):
            integ = sp.expand(eta_s ** k * (1 - eta_s ** 2) ** H
                              * _Dp.subs(_z, eta_s) * sp.exp(-q * eta_s))
            return integrate_poly_exp(integ, eta_s, sp.Integer(-1), sp.Integer(1))

        B1 = {k: beta(k, H1, q1) for k in {k for _j, k in c1}}
        B2 = {k: beta(k, H2, q2) for k in {k for _j, k in c2}}
        if all(v == 0 for v in B1.values()) or all(v == 0 for v in B2.values()):
            terms.append(sp.Integer(0))
            continue

        xi_cache: dict = {}

        def Xi(j1, j2, _tau=tau):
            key = (j1, j2)
            if key not in xi_cache:
                xi_cache[key] = ordered_xi_general(_tau, s, H1, H2, j1, j2, p1, p2)
            return xi_cache[key]

        acc = sp.Integer(0)
        for (j1, k1), cc1 in c1.items():
            if B1[k1] == 0:
                continue
            for (j2, k2), cc2 in c2.items():
                if B2[k2] == 0:
                    continue
                acc += cc1 * cc2 * B1[k1] * B2[k2] * Xi(j1, j2)
        w = ((-1) ** sig1 * (2 * tau + 1)
             * (sp.factorial(tau - s) / sp.factorial(tau + s)) ** 2)
        terms.append(C * w * acc)
    return terms if return_terms else sp.Add(*terms)



def exchange_hp(ZA, oa, ob, ZB, oc, od, R, tau_max: int = 10, dps: int = 40):
    """Same assembly as `exchange_closed_form`, accumulated NUMERICALLY at `dps`.

    Every factor is still evaluated from its own closed form (`ordered_xi_general`
    for the ordered xi integral, `integrate_poly_exp` for the eta halves); only
    the products and the tau sum are done in mpmath.  This avoids building one
    enormous symbolic tree per tau, which is what makes the symbolic route slow
    at the tau values a heteronuclear pair needs.

    Returns (total, per_tau list) as mpmath mpf at the ambient precision.
    """
    from mpmath import mp

    sig1, sig2 = oa[2] - ob[2], od[2] - oc[2]
    if sig1 != sig2:
        return mp.mpf(0), []
    s = abs(sig1)
    R = sp.nsimplify(R)

    P1, h1, p1, q1 = two_center_spheroidal_product(ZA, oa, ZB, ob, R)
    P2, h2, p2, q2 = two_center_spheroidal_product(ZA, oc, ZB, od, R)
    H1 = int(h1 + sp.Rational(s, 2))
    H2 = int(h2 + sp.Rational(s, 2))

    def coeffs_of(P):
        return dict(sp.Poly(sp.expand((xi_s ** 2 - eta_s ** 2) * P),
                            xi_s, eta_s).terms())

    c1, c2 = coeffs_of(P1), coeffs_of(P2)

    def NM(e):
        return mp.mpf(str(sp.N(e, dps + 12)))

    with mp.workdps(dps + 20):
        C = NM((R ** 3 / 8) ** 2 * (2 * sp.pi) ** 2 * (2 / R))
        cc1n = {k: NM(v) for k, v in c1.items()}
        cc2n = {k: NM(v) for k, v in c2.items()}
        per_tau = []
        for tau in range(s, tau_max + 1):
            Dp = sp.expand(_deriv(sp.legendre(tau, _z), s))

            def beta(k, H, q, _Dp=Dp):
                integ = sp.expand(eta_s ** k * (1 - eta_s ** 2) ** H
                                  * _Dp.subs(_z, eta_s) * sp.exp(-q * eta_s))
                return integrate_poly_exp(integ, eta_s, sp.Integer(-1),
                                          sp.Integer(1))

            B1 = {k: NM(beta(k, H1, q1)) for k in {k for _j, k in c1}}
            B2 = {k: NM(beta(k, H2, q2)) for k in {k for _j, k in c2}}
            if all(v == 0 for v in B1.values()) or                all(v == 0 for v in B2.values()):
                per_tau.append(mp.mpf(0))
                continue
            xi_cache: dict = {}
            acc = mp.mpf(0)
            for (j1, k1), v1 in cc1n.items():
                if B1[k1] == 0:
                    continue
                for (j2, k2), v2 in cc2n.items():
                    if B2[k2] == 0:
                        continue
                    key = (j1, j2)
                    if key not in xi_cache:
                        xi_cache[key] = NM(ordered_xi_general(
                            tau, s, H1, H2, j1, j2, p1, p2))
                    acc += v1 * v2 * B1[k1] * B2[k2] * xi_cache[key]
            w = NM((-1) ** sig1 * (2 * tau + 1)
                   * (sp.factorial(tau - s) / sp.factorial(tau + s)) ** 2)
            per_tau.append(C * w * acc)
        total = mp.fsum(per_tau)
    return total, per_tau


__all__ = [
    "s_radial", "A_aux", "B_aux", "I2c", "overlap", "kinetic", "h_core",
    "h_core_eigentrick", "one_center_eri", "exchange_closed_form",
    "exchange_hp", "_inv_r",
    "aabb_closed_form", "hybrid_closed_form", "exchange_value",
]
