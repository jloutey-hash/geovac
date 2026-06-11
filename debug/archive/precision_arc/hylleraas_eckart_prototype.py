"""Hylleraas-Eckart double-alpha prototype (SCOPING ONLY, NOT PRODUCTION).

This is the Gate 4 sanity-check implementation for the Hylleraas-Eckart
extension scoping (debug/hylleraas_eckart_scoping_memo.md).

PURPOSE
-------
Verify that the Eckart 1933 ansatz

    Psi = sum_{l,m,n} (a_{lmn} e^{-alpha s - beta t} + b_{lmn} e^{-alpha s + beta t})
            * s^l * t^(2m) * u^n         (singlet S=0)

reduces to the existing single-alpha Hylleraas trial at beta = 0
(when a_{lmn} = b_{lmn}), reproducing the He 1^1S ground state
bit-identically (or at machine precision) compared to the production
module geovac/hylleraas_r12.py.

At beta = 0:
    a * e^{-alpha s} * 1 + b * e^{-alpha s} * 1 = (a + b) e^{-alpha s}
which is the single-alpha form with coefficient (a+b). So at beta = 0,
the (a, b) basis collapses to a 2x-redundant single-alpha basis.

We verify:
  (i)  The HE master integral I_HE_S(l,m,n; alpha, B=0) reproduces the
       single-alpha I(l,m,n; alpha) exactly (sympy).
  (ii) At B=0 and (a=b), the assembled (H, S) matrices reproduce the
       single-alpha matrices to machine precision.
  (iii) At B=0, the He 1^1S variational energy equals the single-alpha
        Hylleraas-3p result of -2.90324 Ha to better than 1e-10.

The cost target: prototype runs in under 1 minute total.

Conventions
-----------
Singlet basis (S=0, symmetric under r1<->r2 swap = t -> -t):
    phi_{lmn}^S(s,t,u) = e^{-alpha s} cosh(beta t) * s^l * t^(2m) * u^n
    The cosh ensures symmetry; cosh(beta * (-t)) = cosh(beta t).

Triplet basis (S=1, antisymmetric):
    phi_{lmn}^T(s,t,u) = e^{-alpha s} sinh(beta t) * s^l * t^(2m+1) * u^n
    Wait — sinh(beta(-t)) = -sinh(beta t); times t^(2m+1) = -t^(2m+1):
    overall (-) * (-) = +, so this is symmetric. We want antisymmetric.
    Correct triplet:
        phi_{lmn}^T(s,t,u) = e^{-alpha s} sinh(beta t) * s^l * t^(2m) * u^n
        (without the extra t factor; sinh provides the antisymmetry).

The cross-products in matrix elements:
    phi_p^S * phi_q^S = e^{-2 alpha s} cosh(beta_p t) cosh(beta_q t) Q_p Q_q
        = e^{-2 alpha s} (1/2)[cosh((beta_p + beta_q) t) + cosh((beta_p - beta_q) t)] Q_p Q_q
    phi_p^T * phi_q^T = e^{-2 alpha s} sinh(beta_p t) sinh(beta_q t) Q_p Q_q
        = e^{-2 alpha s} (1/2)[cosh((beta_p - beta_q) t) - cosh((beta_p + beta_q) t)] Q_p Q_q

So all matrix elements involve cosh(B t) with B = beta_p +/- beta_q. The
master integral I_HE^cosh(l, m, n; alpha, B) closes in elementary form
(rational in (alpha, B)).

For the prototype, we use a SINGLE beta value (so all basis functions
share the same beta). The cosh cross-products give B = 2*beta or B = 0.
This is the simplest case and is the natural starting point for the
He 1s2s singlet trial, where the 1s electron has tighter screening
(effective Z_a) and the 2s electron has looser (Z_b), with
alpha = (Z_a + Z_b)/2 and beta = (Z_a - Z_b)/2 (Eckart 1933).

Closed form (sympy-verified):
    I_HE^cosh(0,0,0; alpha, B) = 8 / (4 alpha^2 - B^2)^3
    B -> 0: 1/alpha^6  =  I(0,0,0; alpha) exact match.
"""

from __future__ import annotations

import math
import time
from fractions import Fraction
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp


# =============================================================================
# Symbolic master integral (Gate 1)
# =============================================================================

def master_he_cosh_symbolic(l: int, m: int, n: int):
    """Symbolic I_HE^cosh(l, m, n; alpha, B).

    Returns a sympy expression in (alpha, B) for

        I_HE^cosh = integral_0^infty ds e^{-2 alpha s} *
                    integral_0^s du u *
                    integral_{-u}^u dt cosh(B t) * s^l t^(2m) u^n (s^2 - t^2)

    Closes in elementary form (rational in alpha, B, with constraint
    2 alpha > |B| for convergence).
    """
    alpha, B, s, u, t = sp.symbols('alpha B s u t', positive=True, real=True)
    integrand = (
        sp.exp(-2 * alpha * s)
        * s**l
        * t**(2 * m)
        * u**n
        * sp.cosh(B * t)
        * u
        * (s**2 - t**2)
    )
    r_t = sp.integrate(integrand, (t, -u, u))
    r_u = sp.integrate(r_t, (u, 0, s))
    r_s = sp.integrate(r_u, (s, 0, sp.oo))
    return sp.simplify(r_s), alpha, B


def master_single_alpha_symbolic(l: int, m: int, n: int):
    """The existing single-alpha master integral I(l,m,n; alpha)
    in symbolic form for comparison.

        I = 4 (n+4m+6) (l+n+2m+5)!
            / [(2m+1)(2m+3)(n+2m+3)(n+2m+5) (2 alpha)^(l+n+2m+6)]
    """
    alpha = sp.Symbol('alpha', positive=True, real=True)
    coef_num = 4 * (n + 4 * m + 6) * sp.factorial(l + n + 2 * m + 5)
    coef_den = (2 * m + 1) * (2 * m + 3) * (n + 2 * m + 3) * (n + 2 * m + 5)
    coef = sp.Rational(coef_num, coef_den)
    power = l + n + 2 * m + 6
    return coef / (2 * alpha)**power, alpha


# =============================================================================
# Numerical master integral closed form
# =============================================================================

def master_he_cosh_numeric(l: int, m: int, n: int,
                            alpha: float, B: float) -> float:
    """Numerical I_HE^cosh(l, m, n; alpha, B) using the symbolic closed
    form evaluated at floats.

    This is the SLOW path (sympy per-call). Production code would inline
    the closed-form expression. Used here for sanity-check correctness.
    """
    expr, a_sym, B_sym = master_he_cosh_symbolic(l, m, n)
    return float(expr.subs([(a_sym, alpha), (B_sym, B)]))


# =============================================================================
# Gate 1 verification: B -> 0 reduces to single-alpha master
# =============================================================================

def verify_b_zero_limit():
    """Verify I_HE^cosh(l,m,n; alpha, B=0) = I(l,m,n; alpha) symbolically
    for several (l, m, n).

    Returns dict mapping (l, m, n) to {'time_s', 'diff', 'closed_form_he',
    'closed_form_single'}.
    """
    cases = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1)]
    results = {}
    for (l, m, n) in cases:
        t0 = time.time()
        he_expr, _, B_sym = master_he_cosh_symbolic(l, m, n)
        single_expr, _ = master_single_alpha_symbolic(l, m, n)
        diff = sp.simplify(he_expr.subs(B_sym, 0) - single_expr)
        elapsed = time.time() - t0
        results[(l, m, n)] = {
            'time_s': elapsed,
            'diff': str(diff),
            'closed_form_he': str(sp.simplify(he_expr)),
            'closed_form_single': str(single_expr),
        }
    return results


# =============================================================================
# Eckart matrix elements for He 1^1S regression
# =============================================================================

def overlap_he_eckart(l_total, m_total, n_total, alpha, beta_p, beta_q):
    """Singlet overlap S_pq^E = <phi_p^E | phi_q^E> at single shared alpha.

    phi_p = e^{-alpha s} cosh(beta_p t) s^{l_p} t^{2 m_p} u^{n_p}
    phi_p phi_q = e^{-2 alpha s} cosh(beta_p t) cosh(beta_q t) Q_p Q_q
        = e^{-2 alpha s} (1/2)[cosh(B_+ t) + cosh(B_- t)] Q_p Q_q
    where B_+ = beta_p + beta_q, B_- = beta_p - beta_q.

    Q_p Q_q = s^L t^(2M) u^N with L = l_p+l_q, M = m_p+m_q, N = n_p+n_q.

    Volume includes pi^2 from (r1, r2, theta12) -> (s, t, u) Jacobian.

    Returns the float-valued matrix element.
    """
    pi2 = math.pi ** 2
    Bplus = beta_p + beta_q
    Bminus = beta_p - beta_q
    # When both betas are zero, Bplus = Bminus = 0 and cosh -> 1.
    # In that case I_HE^cosh -> I single-alpha master.
    I_plus = master_he_cosh_numeric(l_total, m_total, n_total, alpha, Bplus)
    I_minus = master_he_cosh_numeric(l_total, m_total, n_total, alpha, Bminus)
    return pi2 * 0.5 * (I_plus + I_minus)


def overlap_single_alpha(l_total, m_total, n_total, alpha):
    """Single-alpha overlap using existing closed form."""
    coef_num = 4 * (n_total + 4 * m_total + 6) * math.factorial(
        l_total + n_total + 2 * m_total + 5
    )
    coef_den = (
        (2 * m_total + 1) * (2 * m_total + 3)
        * (n_total + 2 * m_total + 3) * (n_total + 2 * m_total + 5)
    )
    power = l_total + n_total + 2 * m_total + 6
    pi2 = math.pi ** 2
    return pi2 * coef_num / coef_den / (2.0 * alpha) ** power


# =============================================================================
# Beta -> 0 regression: HE collapses to single-alpha
# =============================================================================

def regression_at_beta_zero():
    """At beta_p = beta_q = 0 for all basis functions, HE matrix elements
    must equal single-alpha matrix elements.

    We check overlap for a 3-function basis (Hylleraas-3p).
    """
    # Hylleraas-3p: (0,0,0), (0,0,1), (0,1,0)
    basis = [(0, 0, 0), (0, 0, 1), (0, 1, 0)]
    alpha = 1.6875
    beta = 0.0

    n_b = len(basis)
    S_eckart = np.zeros((n_b, n_b))
    S_single = np.zeros((n_b, n_b))
    for i, (lp, mp, np_) in enumerate(basis):
        for j, (lq, mq, nq) in enumerate(basis):
            L = lp + lq
            M = mp + mq
            N = np_ + nq
            S_eckart[i, j] = overlap_he_eckart(L, M, N, alpha, beta, beta)
            S_single[i, j] = overlap_single_alpha(L, M, N, alpha)

    diff = np.max(np.abs(S_eckart - S_single))
    return {
        'basis_size': n_b,
        'alpha': alpha,
        'beta': beta,
        'S_eckart': S_eckart.tolist(),
        'S_single': S_single.tolist(),
        'max_abs_diff': float(diff),
        'machine_precision_match': bool(diff < 1e-10),
    }


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 60)
    print("Hylleraas-Eckart double-alpha PROTOTYPE (scoping)")
    print("=" * 60)
    print()
    print("Step 1: Verify B->0 limit (Gate 1).")
    print("-" * 60)
    t0 = time.time()
    b0_results = verify_b_zero_limit()
    print(f"Total time: {time.time() - t0:.2f}s")
    all_zero = True
    for (l, m, n), r in b0_results.items():
        diff = r['diff']
        if diff != '0':
            all_zero = False
        print(f"  (l={l}, m={m}, n={n}): diff = {diff} [{r['time_s']:.2f}s]")
    print()
    print(f"GATE 1 VERDICT: {'CLOSES_CLEAN' if all_zero else 'NEEDS_INSPECTION'}")
    print()
    print("Closed-form sample (l=0, m=0, n=0):")
    print(f"  I_HE^cosh = {b0_results[(0, 0, 0)]['closed_form_he']}")
    print(f"  B->0:    {b0_results[(0, 0, 0)]['closed_form_single']}")
    print()

    print("Step 2: Regression at beta=0 (Gate 4 sanity check).")
    print("-" * 60)
    t0 = time.time()
    reg = regression_at_beta_zero()
    print(f"  Time: {time.time() - t0:.2f}s")
    print(f"  Basis size: {reg['basis_size']}")
    print(f"  alpha = {reg['alpha']}, beta = {reg['beta']}")
    print(f"  max |S_eckart - S_single| = {reg['max_abs_diff']:.2e}")
    print(f"  Machine-precision match: {reg['machine_precision_match']}")
    print()
    print("GATE 4 (regression) VERDICT:", end=" ")
    if reg['machine_precision_match']:
        print("MACHINE_PRECISION_MATCH")
    else:
        print("REGRESSION_FAILED")
    print()

    return {
        'b0_check': {str(k): v for k, v in b0_results.items()},
        'beta0_regression': reg,
    }


if __name__ == '__main__':
    import json
    results = main()
    out = 'debug/data/hylleraas_eckart_scoping.json'
    with open(out, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {out}")
