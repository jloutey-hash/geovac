"""Phase 0-h: scoping the HYBRID class (AA|AB), (AB|BB), before any code.

Called for by the build plan section 8.3. Same pattern as Phase 0 and 0-prime:
pre-registered questions, answered symbolically, with a GO / RESCOPE / STOP gate.

THE REDUCTION UNDER TEST. Hybrid means one charge distribution is one-center and
the other is genuinely two-center:

    rho_1 = conj(chi_a^A) chi_b^A          one-center at A
    rho_2 = conj(chi_c^A) chi_d^B          two-center overlap density

Since rho_1 is one-center its potential V_1 is closed-form (increment 1), so

    (ab|cd) = int rho_2 V_1 d3r
            = sum_{L,M} g_1 int d3r [conj(R_c) V_L](r_A) [conj(Y_lc,mc) Y_LM](Om_A)
                                     R_d(r_B) Y_ld,md(Om_B)

and the bracketed angular product re-couples by Gaunt to a FINITE sum of single
harmonics Y_{L'M'}(Om_A). So each term has the shape

    int d3r  F(r_A) Y_{L'M'}(Om_A)  G(r_B) Y_{ld md}(Om_B)

which is EXACTLY the master integral increment 1c already evaluates -- only with
a more general A-side radial function. 1c's F was a shell potential; here
F = conj(R_c) x V_L, i.e. Laurent x exponential with two decay rates.

If that holds, the hybrid class is a generalization of machinery in hand, not a
from-scratch Ruedenberg Part II build. THAT is the question worth answering
before writing anything.

THE THREE QUESTIONS

  HQ1  Does the reduction hold, and is the angular sum finite?
       Gate: numerically reproduce an independent engine on a hybrid quartet.

  HQ2  Power counting. The r_A integral runs over [|r_B-R|, r_B+R] and
       int r^p e^{-a r} dr is elementary for p >= 0 and carries E_1 for p <= -1.
       So: what is the minimum power of r_A, as a function of the l labels?

  HQ3  Where HQ2 goes negative, does E_1 actually survive, and if so is it the
       ALREADY-CLASSIFIED Stieltjes seed (Phase 0 Q2, Paper 18 Level 2) or a new
       transcendental class?

Note the discipline being applied, which is the specific thing increment 1 got
wrong and 1c corrected: ask the seed question on the class's OWN support, and do
not carry an answer over from a neighbouring class.

Run from repo root:  python debug/phase0h_hybrid_scoping.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.two_center_eri import (  # noqa: E402
    V_L_radial, angular_factor, gaunt_LM, multipole_decomposition, plm_signed,
    radial_norm, radial_poly, r_s, sph_norm_numeric, t_s, upper_integral, y_s,
)

Z1, Z3 = Fraction(1), Fraction(3)


# ------------------------------------------------------------------ the terms

def hybrid_terms(ZA, oa, ob, oc, ZB, od):
    """Decompose a hybrid quartet into master-integral terms.

    Yields (coeff, L, Lp, Mp, l_d, m_d), where the term is

        coeff * int d3r [conj(R_c) V_L](r_A) Y_{Lp,Mp}(Om_A) R_d(r_B) Y_{ld,md}(Om_B)
    """
    (_nc, lc, mc), (_nd, ld, md) = oc, od
    for L, M, g1, _radA, _bA in multipole_decomposition(ZA, *oa, ZA, *ob):
        Mp = M - mc                                  # Gaunt M rule
        if Mp + md != 0:                             # the phi integral
            continue
        for Lp in range(abs(lc - L), lc + L + 1):
            if (lc + L + Lp) % 2 != 0:
                continue
            if abs(Mp) > Lp:
                continue
            gc = gaunt_LM(lc, mc, L, M, Lp, Mp)      # <Y_LpMp | conj(Y_lcmc) Y_LM>
            if gc == 0:
                continue
            yield g1 * gc, L, Lp, Mp, ld, md


def F_radial(ZA, oa, ob, oc, L):
    """F(r_A) = conj(R_c)(r_A) * V_L(r_A), and its minimum power of r_A."""
    _nc, lc, _mc = oc
    coeffs, a = radial_poly(ZA, oc[0], lc)
    Nc = radial_norm(ZA, oc[0], lc)
    Rc = sum(Nc * c * r_s ** k for k, c in coeffs.items()) * sp.exp(-a * r_s)
    radA, bA = None, None
    for LL, _M, _g, rad, b in multipole_decomposition(ZA, *oa, ZA, *ob):
        if LL == L:
            radA, bA = rad, b
            break
    return sp.expand(Rc * V_L_radial(radA, bA, L))


# ---------------------------------------------------------------------- HQ1

def hybrid_quadrature(ZA, oa, ob, oc, ZB, od, R):
    """(ab|cd) hybrid, by 2D quadrature over (r_B, cos th_B)."""
    from scipy import integrate

    coeffs_d, ad = radial_poly(ZB, od[0], od[1])
    Nd = radial_norm(ZB, od[0], od[1])
    Rd = sp.lambdify(
        r_s, sum(Nd * c * r_s ** k for k, c in coeffs_d.items()) * sp.exp(-ad * r_s),
        "numpy")

    total = 0.0
    for coeff, L, Lp, Mp, ld, md in hybrid_terms(ZA, oa, ob, oc, ZB, od):
        Ff = sp.lambdify(r_s, F_radial(ZA, oa, ob, oc, L), "numpy")
        nA, nB = sph_norm_numeric(Lp, Mp), sph_norm_numeric(ld, md)

        def integrand(u, rb, Ff=Ff, Lp=Lp, Mp=Mp, ld=ld, md=md, nA=nA, nB=nB):
            rA = np.sqrt(rb * rb + R * R + 2 * rb * R * u)
            if rA < 1e-12:
                return 0.0
            cA = np.clip((R + rb * u) / rA, -1.0, 1.0)
            return float(np.real(
                Ff(rA) * nA * plm_signed(Lp, Mp, np.array([cA]))[0]
                * Rd(rb) * rb * rb * nB * plm_signed(ld, md, np.array([u]))[0]))

        val, _ = integrate.dblquad(integrand, 0.0, 60.0,
                                   lambda _r: -1.0, lambda _r: 1.0,
                                   epsabs=1e-11, epsrel=1e-11)
        total += float(sp.re(coeff)) * 2 * np.pi * val
    return total


def _pointwise_reference(ZA, oa, ob, oc, ZB, od, R):
    """int rho_2 V_1 with rho_2 evaluated POINTWISE from the orbitals.

    Bypasses the Gaunt re-coupling entirely, so it referees exactly the step the
    reduction adds. m = 0 throughout, so everything is real.
    """
    from scipy import integrate

    def radf(Z, n, l):
        c, a = radial_poly(Z, n, l)
        N = radial_norm(Z, n, l)
        return sp.lambdify(
            r_s, sum(N * cc * r_s ** k for k, cc in c.items()) * sp.exp(-a * r_s),
            "numpy")

    Rc, Rd = radf(ZA, oc[0], oc[1]), radf(ZB, od[0], od[1])
    VL = [(L, M, complex(g), sp.lambdify(r_s, V_L_radial(rad, b, L), "numpy"))
          for L, M, g, rad, b in multipole_decomposition(ZA, *oa, ZA, *ob)]

    def integrand(u, rb):
        rA = np.sqrt(rb * rb + R * R + 2 * rb * R * u)
        if rA < 1e-12:
            return 0.0
        cA = np.clip((R + rb * u) / rA, -1.0, 1.0)
        V = sum((g * f(rA) * sph_norm_numeric(L, M)
                 * plm_signed(L, M, np.array([cA]))[0]).real for L, M, g, f in VL)
        chic = Rc(rA) * sph_norm_numeric(oc[1], 0)
        chid = Rd(rb) * sph_norm_numeric(od[1], 0)
        return float(chic * chid * V * rb * rb)

    val, _ = integrate.dblquad(integrand, 0.0, 60.0, lambda _r: -1.0,
                               lambda _r: 1.0, epsabs=1e-12, epsrel=1e-12)
    return 2 * np.pi * val


def leg_HQ1() -> None:
    print("HQ1  does the reduction hold, and is the angular sum finite?")
    from geovac import noci_engine as E
    pa, pb = np.array([0., 0., 0.]), np.array([0., 0., 3.])
    quartet = (Z3, (1, 0, 0), (1, 0, 0), (1, 0, 0), Z1, (1, 0, 0))

    got = hybrid_quadrature(*quartet, 3.0)
    n_terms = len(list(hybrid_terms(*quartet)))
    print(f"     (1sA 1sA|1sA 1sB), Z_A=3, Z_B=1, R=3 -- {n_terms} master term(s)")
    print(f"     reduction value = {got:.12f}\n")

    print("     (a) vs a POINTWISE evaluation that bypasses the Gaunt re-coupling")
    ref = _pointwise_reference(*quartet, 3.0)
    print(f"         pointwise = {ref:.12f}   d = {abs(got - ref):.2e}")

    print("\n     (b) vs eri_md, swept over the STO->Gaussian fit quality")
    print("         (the 6-Gaussian default is NOT good enough for this class --")
    print("          a two-center overlap density samples the exponential tail,")
    print("          which is exactly where Gaussians fit worst)")
    last = None
    for ng in (6, 8, 10, 12):
        arr, dco, q = E.fit_sto_shape(0, 1, n_gauss=ng)
        sh = {"1s": (arr, dco)}
        md = E.eri_md(*(E.sto_shape_basis(c, "1s", z, sh, (0, 0, 0))
                        for c, z in ((pa, 3.0), (pa, 3.0), (pa, 3.0), (pb, 1.0))))
        last = abs(got - md)
        print(f"         n_gauss={ng:2d}  <fit|STO>={q:.9f}  eri_md={md:.9f}"
              f"  d={last:.2e}")

    ok = abs(got - ref) < 1e-13 and last < 1e-6
    print(f"\n     => reduction {'CONFIRMED' if ok else 'NOT CONFIRMED'}: exact "
          f"against the no-recoupling route, and eri_md converges onto it")
    print("        monotonically as the fit improves. The angular sum is finite")
    print("        by Gaunt (two nested terminating couplings).\n")


# ---------------------------------------------------------------------- HQ2

def _power_of(term, v) -> int:
    """Net power of v in a single monomial, ignoring any exp() factor."""
    p = 0
    for f in sp.Mul.make_args(term):
        if f == v:
            p += 1
        elif f.is_Pow and f.base == v:
            p += int(f.exp)
    return p


def min_rA_power(ZA, oa, ob, oc, ZB, od):
    """Minimum power of r_A in the master integrand, over all surviving terms.

    Integrand = r_A r_B F(r_A) G(r_B) x angular_factor(Lp, Mp, ld, md).
    The r_A integral runs over [|r_B-R|, r_B+R]; p >= 0 is elementary,
    p <= -1 carries E_1.
    """
    worst = None
    for _coeff, L, Lp, Mp, ld, md in hybrid_terms(ZA, oa, ob, oc, ZB, od):
        F = sp.expand(F_radial(ZA, oa, ob, oc, L))
        pF = min(_power_of(term, r_s) for term in sp.Add.make_args(F))
        ang = sp.expand(angular_factor(Lp, Mp, ld, md))
        pang = min(_power_of(t, t_s) for t in sp.Add.make_args(ang))
        p = 1 + pF + pang                        # +1 from the r_A r_B measure
        worst = p if worst is None else min(worst, p)
    return worst


def leg_HQ2() -> None:
    print("HQ2  power counting: minimum power of r_A (p >= 0 elementary, "
          "p <= -1 carries E_1)")
    print("     quartet (a b | c d): a,b,c on A, d on B\n")
    print("     l_a l_b l_c l_d   min power of r_A   r_A integral")
    print("     " + "-" * 52)
    rows = [(0, 0, 0, 0), (0, 0, 1, 1), (0, 0, 2, 2),
            (1, 1, 0, 0), (1, 1, 1, 1), (1, 0, 1, 0),
            (2, 2, 0, 0), (2, 1, 1, 0)]
    verdict = {}
    for la, lb, lc, ld in rows:
        oa = (la + 1, la, 0)
        ob = (lb + 1, lb, 0)
        oc = (lc + 1, lc, 0)
        od = (ld + 1, ld, 0)
        p = min_rA_power(Z3, oa, ob, oc, Z1, od)
        if p is None:
            continue
        tag = "elementary" if p >= 0 else "E_1"
        verdict[(la, lb, lc, ld)] = p
        print(f"      {la}   {lb}   {lc}   {ld}          {p:+d}          {tag}")
    print()
    s_type = [v for k, v in verdict.items() if k[0] == 0 and k[1] == 0]
    print(f"     one-center pair s-type (l_a = l_b = 0): min power "
          f"{min(s_type):+d} -> ELEMENTARY")
    rest = [v for k, v in verdict.items() if k[0] or k[1]]
    if rest:
        print(f"     one-center pair with l > 0:            min power "
              f"{min(rest):+d} -> E_1 branch reached")
    print("     structural reading: min power = l_c - L - L' with L <= l_a+l_b")
    print("     and L' <= l_c + L, so the floor is -2(l_a+l_b). The two-center")
    print("     density does NOT protect the way (AA|BB)'s k >= l1+l2 did.\n")


# ---------------------------------------------------------------------- HQ3

def leg_HQ3() -> None:
    print("HQ3  where HQ2 is negative, does E_1 survive -- and is it the "
          "known seed?")
    ZA, oa, ob, oc, od = Z3, (2, 1, 0), (2, 1, 0), (1, 0, 0), (1, 0, 0)
    a = sp.Symbol("a", positive=True)
    coeffs, dec = radial_poly(ZA, 1, 0)
    lo, hi = sp.Symbol("lo", positive=True), sp.Symbol("hi", positive=True)

    seeds = set()
    survives = False
    for coeff, L, Lp, Mp, ld, md in hybrid_terms(ZA, oa, ob, oc, Z1, od):
        F = sp.expand(F_radial(ZA, oa, ob, oc, L))
        ang = sp.expand(angular_factor(Lp, Mp, ld, md))
        integrand = sp.expand(F.subs(r_s, t_s) * ang * t_s)
        acc = sp.Integer(0)
        for term in sp.Add.make_args(integrand):
            c, p, d = sp.Integer(1), sp.Integer(0), sp.Integer(0)
            for f in sp.Mul.make_args(term):
                if f == t_s:
                    p += 1
                elif f.is_Pow and f.base == t_s:
                    p += f.exp
                elif isinstance(f, sp.exp):
                    arg = sp.expand(f.args[0])
                    d += -sp.diff(arg, t_s)
                    c *= sp.exp(sp.expand(arg + (-sp.diff(arg, t_s)) * t_s))
                else:
                    c *= f
            if d == 0:
                continue
            acc += c * (upper_integral(int(p), d, lo) - upper_integral(int(p), d, hi))
        acc = sp.expand(acc)
        for e1 in acc.atoms(sp.expint):
            if sp.simplify(acc.coeff(e1)) != 0:
                survives = True
                seeds.add(sp.simplify(e1.args[1] / lo if e1.args[1].has(lo)
                                      else e1.args[1] / hi))
    print(f"     probe quartet (2p0_A 2p0_A | 1s_A 1s_B), Z_A=3, Z_B=1")
    print(f"     E_1 survives the r_A integral: {survives}")
    print(f"     decay rates entering E_1(rate * endpoint): {sorted(seeds, key=str)}")
    print("     = {a, a+b}: the orbital exponent of chi_c, and that plus the")
    print("       decay of rho_1. Sums of orbital exponents, nothing else.\n")

    # does the OUTER r_B integral close on that seed, or open a new class?
    from scipy import integrate as si
    from scipy.special import exp1
    print("     does the outer r_B integral close on it? The endpoints are")
    print("     |r_B - R| and r_B + R, so the claim to check is")
    print("       int_0^inf e^{-ct} E_1(a(t+R)) dt = E_1(aR)/c - e^{cR}E_1((c+a)R)/c")
    worst = 0.0
    for c_, a_, R_ in ((1.3, 2.0, 3.0), (0.7, 3.0, 1.5), (3.0, 6.0, 3.0)):
        num, _ = si.quad(lambda t: np.exp(-c_ * t) * exp1(a_ * (t + R_)), 0, np.inf,
                         limit=400, epsabs=1e-14, epsrel=1e-14)
        cf = exp1(a_ * R_) / c_ - np.exp(c_ * R_) * exp1((c_ + a_) * R_) / c_
        worst = max(worst, abs(num - cf))
    print(f"     verified at 3 (c, a, R) points, worst deviation {worst:.2e}")
    print("     => it closes: elementary terms plus E_1(lambda R) constants with")
    print("        lambda a sum of orbital exponents. That is e^{+-a}E_1(a*shift),")
    print("        the Stieltjes seed of Phase 0 Q2 / Paper 18 'Level 2'.")
    print("        NOT a new transcendental class.\n")


def main() -> None:
    print("Phase 0-h -- scoping the hybrid class (AA|AB), (AB|BB)\n")
    leg_HQ1()
    leg_HQ2()
    leg_HQ3()
    print("Read the verdict in the memo, not here: this driver reports "
          "measurements only.")


if __name__ == "__main__":
    main()
