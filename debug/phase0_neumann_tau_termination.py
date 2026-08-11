"""Phase 0, question 1: does the Neumann tau sum TERMINATE for fixed orbital l?

This is the question the whole general-m build hangs on (see
docs/neumann_general_m_build_plan.md section 2). If the tau sum terminates
exactly, the engine is exact-up-to-seeds with NO convergence study. If not, the
deliverable weakens to high-precision-with-quantified-error.

THE ARGUMENT TO TEST
--------------------
In the Neumann expansion the eta-dependence of electron 1 is

    P_tau^sigma(eta) * [angular factors of orbitals a and c in eta]

Associated Legendre functions factor as

    P_l^m(eta) = (1 - eta^2)^{|m|/2} * (polynomial of degree l - |m|)

so each of the three factors carries a half-integer power of (1 - eta^2) when its
order is odd. The azimuthal integral forces sigma = m_a - m_c, and the claim is
that this is exactly the condition under which the three half-integer powers
COMBINE TO AN INTEGER, leaving a genuine polynomial integrand. Orthogonality of
{P_tau^sigma}_tau on [-1, 1] at fixed sigma then kills every tau above the
polynomial degree -- i.e. termination, with

    tau_max = l_a + l_c        (per side)

which is the same shape as the multipole termination geovac/shibuya_wulfman.py
already achieves for cross-center V_ne at L_max = l1 + l2 (Q-B verified), and as
Paper 19 records for the balanced cross-center potential at L_max = 2 l_max.

WHAT THIS SCRIPT DOES
---------------------
Computes I(tau) = integral_{-1}^{1} P_tau^sigma(eta) P_{l_a}^{m_a}(eta)
P_{l_c}^{m_c}(eta) d(eta) in EXACT symbolic form, with sigma = |m_a - m_c|, and
reports the largest tau with I(tau) != 0. Exact integration, so a reported zero
is a real zero and not a tolerance artifact.

Run from repo root:  python debug/phase0_neumann_tau_termination.py
"""

from __future__ import annotations

import sympy as sp

eta = sp.Symbol("eta", real=True)


def assoc_legendre(l: int, m: int):
    """P_l^m(eta) as an exact sympy expression (Condon-Shortley, m >= 0)."""
    m = abs(m)
    if m > l:
        return sp.Integer(0)
    x = sp.Symbol("_x")
    # Rodrigues-type: P_l^m(x) = (-1)^m (1-x^2)^{m/2} d^m/dx^m P_l(x)
    expr = (-1) ** m * (1 - x ** 2) ** sp.Rational(m, 2) * sp.diff(
        sp.legendre(l, x), x, m)
    return sp.simplify(expr.subs(x, eta))


def eta_integral(tau: int, sigma: int, la: int, ma: int, lc: int, mc: int):
    """Exact integral_{-1}^{1} P_tau^sigma P_la^ma P_lc^mc d(eta)."""
    f = (assoc_legendre(tau, sigma) * assoc_legendre(la, ma)
         * assoc_legendre(lc, mc))
    if f == 0:
        return sp.Integer(0)
    return sp.simplify(sp.integrate(sp.expand(f), (eta, -1, 1)))


def probe(la, ma, lc, mc, tau_ceiling=None):
    sigma = abs(ma - mc)
    predicted = la + lc
    ceiling = tau_ceiling if tau_ceiling is not None else predicted + 3
    nonzero = []
    for tau in range(0, ceiling + 1):
        val = eta_integral(tau, sigma, la, ma, lc, mc)
        if val != 0:
            nonzero.append((tau, val))
    highest = max((t for t, _ in nonzero), default=None)
    return sigma, predicted, highest, nonzero


CASES = [
    # (l_a, m_a, l_c, m_c) -- sigma = |m_a - m_c| is forced by the phi integral
    (0, 0, 0, 0),      # s|s        sigma=0   (the existing m=0 sector)
    (1, 0, 1, 0),      # p0|p0      sigma=0
    (1, 1, 1, 1),      # p+1|p+1    sigma=0, but both factors carry sqrt
    (1, 1, 0, 0),      # p+1|s      sigma=1  <- first genuinely new sector
    (1, 1, 1, 0),      # p+1|p0     sigma=1
    (2, 1, 1, 0),      # d+1|p0     sigma=1
    (2, 2, 0, 0),      # d+2|s      sigma=2
    (2, 2, 1, 1),      # d+2|p+1    sigma=1
    (2, 2, 2, 0),      # d+2|d0     sigma=2
    (3, 2, 2, 1),      # f+2|d+1    sigma=1
]


def main() -> None:
    print("Phase 0 Q1 -- Neumann tau termination in the eta (angular) half")
    print("sigma = |m_a - m_c| as forced by the azimuthal integral.")
    print("Exact symbolic integration: a reported zero IS a zero.\n")
    print(f"{'(la,ma)|(lc,mc)':<18}{'sigma':>6}{'la+lc':>7}"
          f"{'highest nonzero tau':>21}{'verdict':>12}")
    print("-" * 66)

    all_terminate = True
    at_prediction = True
    for la, ma, lc, mc in CASES:
        sigma, predicted, highest, _nz = probe(la, ma, lc, mc)
        if highest is None:
            verdict = "all zero"
        elif highest <= predicted:
            verdict = "TERMINATES"
            if highest != predicted:
                at_prediction = False
        else:
            verdict = "EXCEEDS"
            all_terminate = False
        label = f"({la},{ma})|({lc},{mc})"
        print(f"{label:<18}{sigma:>6}{predicted:>7}{str(highest):>21}"
              f"{verdict:>12}")

    print("-" * 66)
    print(f"\nAll cases terminate at or below l_a + l_c : {all_terminate}")
    print(f"Bound is TIGHT (highest == l_a + l_c) in all : {at_prediction}")
    if all_terminate:
        print("\n=> Phase 0 Q1 GO on the eta half: the tau sum terminates, so the")
        print("   Neumann series is a FINITE sum and no convergence study is")
        print("   needed for this half. tau_max = l_a + l_c per side.")
    else:
        print("\n=> Phase 0 Q1 NEGATIVE: termination fails; Phase 3 becomes a")
        print("   truncation-error study and the deliverable weakens.")


if __name__ == "__main__":
    main()
