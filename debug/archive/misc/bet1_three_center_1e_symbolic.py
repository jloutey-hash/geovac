"""Bet 1, the headline: the 3-centre 1e integral in CLOSED FORM -- what weight?

The pinned xi integral (validated numerically in bet1_..._assembly.py) is, per tau
and per density monomial xi^{j1}(xi^2-1)^{H1} e^{-p xi}:

    Q_tau(s0) * int_1^{s0} [poly * P_tau] e^{-p xi} dxi              (INNER, elementary)
  + P_tau(s0) * int_{s0}^inf [poly*(-W_tau)] e^{-p xi} dxi           (OUTER elem)
  + P_tau(s0) * int_{s0}^inf [poly * P_tau] Q0(xi) e^{-p xi} dxi     (OUTER log)

with s0 = xi_X, Q_tau = P_tau*Q0 - W_tau, Q0 = (1/2)ln((xi+1)/(xi-1)).

The transcendental content lives ONLY in the last line. Building it through the
engine's validated weight-1 moment machinery (finite_power_exp, upper_integral,
log_shift_moment) makes the seed set explicit:

  - INNER: finite-range poly*exp -> elementary (exp).
  - OUTER elem: poly*exp on [s0,inf) -> exp, E_1 via upper_integral.
  - OUTER log: substitute t = xi - s0; since s0 = xi_X > 1 the shifts are
    s0 +- 1 > 0, so log_shift_moment(m, p, s0+-1) applies -- {E_1, ln}, and
    crucially NO Euler gamma (the gamma of the two-centre exchange came from the
    xi = 1 endpoint, which a source pinned at xi_X > 1 never reaches).

So the prediction is sharper than registered: weight one, AND gamma-free for a
generic (off-foci) third centre.  This script (1) inspects the atoms with p, s0
SYMBOLIC -> the functional weight verdict, and (2) evaluates the full closed form
at the reference geometry -> correctness vs -0.341962.

Run:  python debug/bet1_three_center_1e_symbolic.py
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
    _legendre_Q, _z, eta_s, finite_power_exp, integrate_poly_exp,
    log_shift_moment, two_center_spheroidal_product, upper_integral, xi_s,
)

_xi = sp.Symbol("xi", positive=True)


def xi_pinned_symbolic(j1, H1, tau, p, s0):
    """Closed form of the pinned xi integral for one monomial, one tau.

    p = density xi-exponent, s0 = xi_X.  Returns a sympy expression whose only
    transcendental functions are meant to be {exp, expint(=E_1), log}.
    """
    poly = sp.expand(_xi ** j1 * (_xi ** 2 - 1) ** H1)
    Ptau = sp.legendre(tau, _xi)
    W = sum(sp.legendre(k - 1, _xi) * sp.legendre(tau - k, _xi) / sp.Integer(k)
            for k in range(1, tau + 1))                      # 0 when tau==0

    Q0_s0 = sp.log((s0 + 1) / (s0 - 1)) / 2
    Ptau_s0 = Ptau.subs(_xi, s0)
    Qtau_s0 = Ptau_s0 * Q0_s0 - (W.subs(_xi, s0) if tau else sp.Integer(0))

    # INNER: Q_tau(s0) * int_1^{s0} poly*P_tau e^{-p xi} dxi  (elementary)
    pP = sp.Poly(sp.expand(poly * Ptau), _xi)
    I_inner = sum(b * finite_power_exp(n, p, sp.Integer(1), s0)
                  for (n,), b in pP.terms())
    inner = Qtau_s0 * I_inner

    # OUTER elem: P_tau(s0) * int_{s0}^inf poly*(-W) e^{-p xi} dxi
    negW = sp.Poly(sp.expand(-poly * W), _xi) if tau else sp.Poly(0, _xi)
    I_outer_elem = sum(b * upper_integral(n, p, s0) for (n,), b in negW.terms())

    # OUTER log: P_tau(s0) * int_{s0}^inf poly*P_tau * (1/2)[ln(xi+1)-ln(xi-1)] e^{-p xi}
    I_outer_log = sp.Integer(0)
    for (n,), b in pP.terms():
        Lp = sp.exp(-p * s0) * sum(sp.binomial(n, m) * s0 ** (n - m)
                                   * log_shift_moment(m, p, s0 + 1)
                                   for m in range(n + 1))
        Lm = sp.exp(-p * s0) * sum(sp.binomial(n, m) * s0 ** (n - m)
                                   * log_shift_moment(m, p, s0 - 1)
                                   for m in range(n + 1))
        I_outer_log += b * sp.Rational(1, 2) * (Lp - Lm)

    return inner + Ptau_s0 * (I_outer_elem + I_outer_log)


def weight_inspection():
    print("=== (1) WEIGHT INSPECTION  (p, s0 symbolic) ===")
    p, s0 = sp.symbols("p s0", positive=True)
    for (j1, H1, tau) in [(0, 0, 0), (0, 0, 1), (1, 0, 2), (2, 1, 2), (0, 0, 3)]:
        e = xi_pinned_symbolic(j1, H1, tau, p, s0)
        funcs = {type(f).__name__ for f in e.atoms(sp.Function)}
        bad = funcs & {"polylog", "dilog", "lerchphi", "zeta"}
        gamma = e.has(sp.EulerGamma)
        print(f"  (j1={j1},H1={H1},tau={tau}): funcs={sorted(funcs)}"
              f"  gamma={gamma}  weight2={'YES' if bad else 'no'}")
    print("  -> expect {exp, expint, log}, gamma=False, weight2=no everywhere\n")


def correctness():
    print("=== (2) CORRECTNESS  (full closed form vs reference) ===")
    # reference geometry
    Yc, Zc, Xc, ZX = np.array([0., 0, 0]), np.array([0., 0, 2]), np.array([1.3, 0, .7]), 1.0
    R_YZ = 2.0
    rY, rZ = np.linalg.norm(Xc - Yc), np.linalg.norm(Xc - Zc)
    xiXf, etaXf = (rY + rZ) / R_YZ, (rY - rZ) / R_YZ
    s0 = sp.nsimplify(round(xiXf, 13), rational=True)     # rational xi_X for speed
    Z1 = Fraction(1)

    cases = {"A 1s_Y x 1s_Z": (((1, 0, 0), (1, 0, 0)), -0.341962099253),
             "B 2p0_Y x 1s_Z": (((2, 1, 0), (1, 0, 0)), -0.186874906193)}
    C_pinned = (R_YZ ** 3 / 8) * (2 * np.pi) * (2.0 / R_YZ)

    for name, ((oY, oZ), ref) in cases.items():
        P1, h1, p1, q1 = two_center_spheroidal_product(Z1, oY, Z1, oZ, sp.Rational(2))
        H1 = int(h1)
        # fold volume factor
        c1 = {}
        for (j, k), c in sp.Poly(sp.expand((xi_s ** 2 - eta_s ** 2) * P1),
                                 xi_s, eta_s).terms():
            c1[(j, k)] = c
        p1r = sp.nsimplify(p1)
        total = 0.0
        for tau in range(0, 16):
            Ptau = sp.legendre(tau, _z)
            PtauEtaX = float(Ptau.subs(_z, etaXf))
            # eta integrals (elementary, numeric ok -- not the weight question)
            B1 = {}
            for k in {kk for _j, kk in c1}:
                integrand = sp.expand(eta_s ** k * (1 - eta_s ** 2) ** H1
                                      * Ptau.subs(_z, eta_s) * sp.exp(-q1 * eta_s))
                B1[k] = float(integrate_poly_exp(integrand, eta_s,
                                                 sp.Integer(-1), sp.Integer(1)))
            acc = 0.0
            for (j1, k1), cc1 in c1.items():
                xip = xi_pinned_symbolic(j1, H1, tau, p1r, s0)
                acc += float(cc1) * B1[k1] * PtauEtaX * float(sp.re(sp.N(xip, 30)))
            total += (2 * tau + 1) * acc
        val = -ZX * C_pinned * total
        print(f"  {name}: closed form={val:+.9f}  ref={ref:+.9f}  |diff|={abs(val-ref):.2e}")


if __name__ == "__main__":
    weight_inspection()
    correctness()
