"""Increment 3b: exchange class assembled for general (l, m) -- eta half CLOSED.

Builds on 3a's verified (xi, eta) expansion. Phase 0-e validated the exchange
route end to end, but only for sigma = 0 and 1s orbitals; this generalizes the
assembly and lets the three scoping claims (termination criterion, seed set,
term count) be re-checked off that easy corner.

THE FACTORIZATION. With rho_1 = conj(chi_a^A) chi_b^B and rho_2 = conj(chi_c^A)
chi_d^B, 3a gives each as P_i(xi,eta) [(xi^2-1)(1-eta^2)]^{h_i} e^{-p_i xi - q_i eta}
times its azimuthal phase. The two phi integrals force

    sigma = m_a - m_b   AND   sigma = m_d - m_c

(the second from electron 2), so a single sigma is fixed by the labels and the
quartet vanishes unless the two agree -- that is M_L conservation.

The kernel contributes P_tau^{|s|}(xi_<) Q_tau^{|s|}(xi_>) P_tau^{|s|}(eta_1)
P_tau^{|s|}(eta_2). Writing P_tau^mu(x) = (x^2-1)^{mu/2} d^mu P_tau/dx^mu on
(1,oo) and (1-x^2)^{mu/2} d^mu P_tau/dx^mu on (-1,1), EACH electron picks up
exactly one (.)^{|s|/2} on each of its xi and eta -- whichever side of the
ordering it lands on. So the half-powers combine to

    H_i = h_i + |sigma|/2

which 3a's parity fact makes an INTEGER. After that nothing fractional survives,
the eta halves factor completely, and only the xi half stays coupled:

    (ab|cd) = C sum_tau w_tau sum_{j1 k1} sum_{j2 k2}
                 c1_{j1k1} c2_{j2k2} Beta_1(k1,tau) Beta_2(k2,tau) Xi(j1,j2,tau)

with c_i the coefficients of (xi^2 - eta^2) P_i (the volume factor folded in).

  Beta_i(k,tau) = int_{-1}^{1} eta^k (1-eta^2)^{H_i} D_tau(eta) e^{-q_i eta} d eta
                  -- polynomial x exponential, so CLOSED FORM here.
  Xi(j1,j2,tau)  = the ordered double integral over xi, still NUMERICAL. Closing
                  it is the remaining hard part of increment 3.

Run from repo root:  python debug/inc3b_exchange_assembly.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from scipy import integrate

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.two_center_eri import (  # noqa: E402
    eta_s, integrate_poly_exp, two_center_spheroidal_product, xi_s,
)

Z1, Z2, Z3 = Fraction(1), Fraction(2), Fraction(3)
_z = sp.Symbol("_z")


def _legendre_Q(tau: int):
    """Q_tau(z) on (1, oo): P_tau Q_0 - sum_{k=1}^{tau} P_{k-1} P_{tau-k}/k."""
    Q0 = sp.log((_z + 1) / (_z - 1)) / 2
    out = sp.legendre(tau, _z) * Q0
    for k in range(1, tau + 1):
        out -= sp.legendre(k - 1, _z) * sp.legendre(tau - k, _z) / k
    return out


def _deriv(expr, mu: int):
    return sp.diff(expr, _z, mu) if mu else expr


def exchange_value(ZA, oa, ob, ZB, oc, od, Rv: float, tau_max: int = 10):
    """Exchange (ab|cd): a,c on A; b,d on B. eta half closed, xi half numerical."""
    sig1 = oa[2] - ob[2]
    sig2 = od[2] - oc[2]
    if sig1 != sig2:
        return 0.0                      # M_L conservation kills it
    s = abs(sig1)

    R = sp.Rational(str(Rv))
    P1, h1, p1, q1 = two_center_spheroidal_product(ZA, oa, ZB, ob, R)
    P2, h2, p2, q2 = two_center_spheroidal_product(ZA, oc, ZB, od, R)
    H1, H2 = h1 + sp.Rational(s, 2), h2 + sp.Rational(s, 2)
    assert H1.is_integer and H2.is_integer, f"half-powers {H1}, {H2} not integral"
    H1, H2 = int(H1), int(H2)

    # fold the volume factor (xi^2 - eta^2) in, then split into monomials
    def coeffs_of(P):
        d = {}
        for (j, k), c in sp.Poly(sp.expand((xi_s ** 2 - eta_s ** 2) * P),
                                 xi_s, eta_s).terms():
            d[(j, k)] = c
        return d

    c1, c2 = coeffs_of(P1), coeffs_of(P2)
    p1f, p2f, q1f, q2f = float(p1), float(p2), float(q1), float(q2)

    total = 0.0
    for tau in range(s, tau_max + 1):
        Dp = sp.expand(_deriv(sp.legendre(tau, _z), s))
        Dq = sp.expand(_deriv(_legendre_Q(tau), s))

        # --- eta halves, CLOSED FORM (polynomial x exponential on [-1, 1])
        def beta(k, H, q):
            integrand = sp.expand(eta_s ** k * (1 - eta_s ** 2) ** H
                                  * Dp.subs(_z, eta_s) * sp.exp(-q * eta_s))
            return float(integrate_poly_exp(integrand, eta_s,
                                            sp.Integer(-1), sp.Integer(1)))

        B1 = {k: beta(k, H1, q1) for k in {k for _j, k in c1}}
        B2 = {k: beta(k, H2, q2) for k in {k for _j, k in c2}}
        if all(abs(v) < 1e-300 for v in B1.values()) or \
           all(abs(v) < 1e-300 for v in B2.values()):
            continue

        # --- xi half, NUMERICAL (ordered P/Q split)
        Dpf = sp.lambdify(_z, Dp, "numpy")
        Dqf = sp.lambdify(_z, Dq, "numpy")

        def xi_double(j1, j2):
            def outer(x1):
                def lo_int(x2):
                    return (x2 ** j2 * (x2 ** 2 - 1) ** H2 * np.exp(-p2f * x2)
                            * float(Dpf(x2)))

                def hi_int(x2):
                    return (x2 ** j2 * (x2 ** 2 - 1) ** H2 * np.exp(-p2f * x2)
                            * float(Dqf(x2)))
                a, _ = integrate.quad(lo_int, 1.0, x1, epsabs=1e-12, epsrel=1e-11,
                                      limit=120)
                b, _ = integrate.quad(hi_int, x1, np.inf, epsabs=1e-12,
                                      epsrel=1e-11, limit=120)
                f1 = x1 ** j1 * (x1 ** 2 - 1) ** H1 * np.exp(-p1f * x1)
                return f1 * (float(Dqf(x1)) * a + float(Dpf(x1)) * b)
            v, _ = integrate.quad(outer, 1.0, np.inf, epsabs=1e-11, epsrel=1e-10,
                                  limit=120)
            return v

        acc = 0.0
        for (j1, k1), cc1 in c1.items():
            if abs(B1[k1]) < 1e-300:
                continue
            for (j2, k2), cc2 in c2.items():
                if abs(B2[k2]) < 1e-300:
                    continue
                acc += (float(cc1) * float(cc2) * B1[k1] * B2[k2]
                        * xi_double(j1, j2))

        w = ((-1) ** sig1 * (2 * tau + 1)
             * float(sp.factorial(tau - s) / sp.factorial(tau + s)) ** 2)
        total += w * acc

    C = (Rv ** 3 / 8) ** 2 * (2 * np.pi) ** 2 * (2.0 / Rv)
    return C * total


# ------------------------------------------------------------------ validation

def main() -> None:
    print("Increment 3b -- exchange assembled for general (l, m)\n")

    print("V1  regression on the Phase 0-e sigma = 0 case")
    print("    (1s_A 1s_B | 1s_A 1s_B), Z_A=3, Z_B=1, R=3")
    got = exchange_value(Z3, (1, 0, 0), (1, 0, 0), Z1, (1, 0, 0), (1, 0, 0), 3.0)
    ref = 0.0063039386
    print(f"    general assembly = {got:.10f}")
    print(f"    Phase 0-e value  = {ref:.10f}   d = {abs(got - ref):.2e}"
          f"   {'OK' if abs(got - ref) < 1e-8 else 'FAIL'}\n")

    print("V2  M_L conservation: sigma_1 != sigma_2 must vanish")
    v = exchange_value(Z3, (2, 1, 1), (1, 0, 0), Z1, (1, 0, 0), (1, 0, 0), 3.0)
    print(f"    (2p+1_A 1s_B | 1s_A 1s_B) = {v:.3e}  "
          f"{'OK (killed)' if abs(v) < 1e-14 else 'NONZERO'}\n")

    print("V3  sigma != 0 against an INDEPENDENT engine")
    print("    V1 only exercises sigma = 0, so it cannot catch an error in any")
    print("    sigma-dependent factor -- the (-1)^sigma, the [(t-s)!/(t+s)!]^2,")
    print("    or the P^mu / Q^mu conventions. This leg is the one that does.\n")
    print("    Route: 2p_{+1} = -(px + i py)/sqrt(2), so")
    print("      (2p+1_A 1s_B | 1s_A 2p+1_B)")
    print("        = 1/2 [ (px_A 1s_B|1s_A px_B) + (py_A 1s_B|1s_A py_B) + i(...) ]")
    print("    and axial symmetry makes the px/py cross terms cancel and the two")
    print("    diagonal terms equal, leaving (px_A 1s_B | 1s_A px_B).\n")

    from geovac import noci_engine as E
    shapes = {}
    for kind, (l, nr) in (("1s", (0, 1)), ("2p", (1, 2))):
        arr, dco, q = E.fit_sto_shape(l, nr, n_gauss=10)   # >= 10, per Phase 0-h
        shapes[kind] = (arr, dco)
    pa, pb = np.array([0., 0., 0.]), np.array([0., 0., 3.])

    def B(c, kind, zeta, lmn):
        return E.sto_shape_basis(c, kind, zeta, shapes, lmn)

    got = exchange_value(Z3, (2, 1, 1), (1, 0, 0), Z1, (1, 0, 0), (2, 1, 1),
                         3.0, tau_max=10)
    md = E.eri_md(B(pa, "2p", 1.5, (1, 0, 0)), B(pb, "1s", 1.0, (0, 0, 0)),
                  B(pa, "1s", 3.0, (0, 0, 0)), B(pb, "2p", 0.5, (1, 0, 0)))
    print(f"    sigma=1  assembly = {got: .10f}")
    print(f"             eri_md   = {md: .10f}   d = {abs(got - md):.2e}")
    print(f"    {'OK (fit-limited)' if abs(got - md) < 1e-5 else 'FAIL'}\n")

    print("V4  a sigma = 0 control through the same Cartesian route")
    got0 = exchange_value(Z3, (2, 1, 0), (1, 0, 0), Z1, (1, 0, 0), (2, 1, 0),
                          3.0, tau_max=10)
    md0 = E.eri_md(B(pa, "2p", 1.5, (0, 0, 1)), B(pb, "1s", 1.0, (0, 0, 0)),
                   B(pa, "1s", 3.0, (0, 0, 0)), B(pb, "2p", 0.5, (0, 0, 1)))
    print(f"    sigma=0  assembly = {got0: .10f}")
    print(f"             eri_md   = {md0: .10f}   d = {abs(got0 - md0):.2e}")
    print(f"    {'OK (fit-limited)' if abs(got0 - md0) < 1e-5 else 'FAIL'}")


if __name__ == "__main__":
    main()
