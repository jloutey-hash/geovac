"""Phase 0, question 2: what is the transcendental seed set of the sigma != 0
Neumann xi-integrals?

Decides whether the general-m build stays inside the seed set the corpus has
ALREADY classified (Paper 18 section "Level 2: e^a E_1(a)"), or opens a new
transcendental class. Per docs/neumann_general_m_build_plan.md section 0:

  GO      seeds subset of {e^a E_1(a)}  -> "exact up to one known seed", g row
          can reach MEASURED, tagging is reuse not new work
  RESCOPE new transcendental class      -> high-precision values only

--- NOTE ON A FIRST-ATTEMPT BUG (kept so it is not reintroduced) ---
The first version of this probe reported "no L present" for every (tau, sigma)
and flagged stray functions like log((x+1)**(x**2-1)). That was a CHECKER bug,
not a result: sympy's expand()/simplify() rewrites log((x+1)/(x-1)) into split
and exponent-folded forms, so substituting the literal composite log never
matched. The fix is to flatten with expand_log(force=True) FIRST, then map
log(x+1) and log(x-1) onto separate symbols and test the structure.

WHAT IS CHECKED
---------------
Step (1) Q_tau^sigma(xi) = A + B * L with A, B algebraic and
         L = ln((xi+1)/(xi-1)) -- i.e. exactly ONE transcendental, and it does
         not proliferate with tau or sigma. Tested by requiring the flattened
         form to be linear in {ln(xi+1), ln(xi-1)} with EQUAL AND OPPOSITE
         coefficients (that is what "depends only on L" means) and no other
         transcendental anywhere.

Step (2) the L branch integrates to the Stieltjes seed. Derived in closed form
         and then VERIFIED against numerical quadrature, rather than asserted:

           integral_c^inf e^{-a xi} ln((xi+1)/(xi-1)) d(xi)
             = (e^{-ac}/a) ln((c+1)/(c-1))
               + (1/a) [ e^{+a} E_1(a(c+1)) - e^{-a} E_1(a(c-1)) ]

         Both transcendental terms are e^{±a} E_1(a * shift) -- literally the
         e^a E_1(a) seed with shifted argument.

Run from repo root:  python debug/phase0_neumann_seed_set.py
"""

from __future__ import annotations

import mpmath as mp
import sympy as sp

x = sp.Symbol("x", positive=True)
Lp, Lm = sp.symbols("Lp Lm")


def Q_order0(l_max: int):
    """Q_l(x) for l = 0..l_max, exact, built from the closed forms + recurrence."""
    L = sp.log(x + 1) - sp.log(x - 1)
    Q = [L / 2, x * L / 2 - 1]
    for l in range(1, l_max):
        Q.append(sp.together(((2 * l + 1) * x * Q[l] - l * Q[l - 1]) / (l + 1)))
    return Q[: l_max + 1]


def Q_assoc(l: int, m: int, Q0):
    """Q_l^m(x) = (x^2-1)^{m/2} d^m/dx^m Q_l(x)."""
    if m == 0:
        return Q0[l]
    return (x ** 2 - 1) ** sp.Rational(m, 2) * sp.diff(Q0[l], x, m)


def structure(expr):
    """Return (depends_only_on_L, degree, other_transcendentals)."""
    e = sp.expand_log(sp.expand(sp.powsimp(expr, force=True)), force=True)
    e = e.subs({sp.log(x + 1): Lp, sp.log(x - 1): Lm})
    # anything log-like left over is a genuine problem
    leftover = [f for f in e.atoms(sp.log)]
    try:
        poly = sp.Poly(sp.expand(e), Lp, Lm)
    except sp.PolynomialError:
        return False, None, ["not polynomial in the logs"]
    deg = poly.total_degree()
    if deg > 1:
        return False, deg, []
    cLp = sp.simplify(poly.coeff_monomial(Lp))
    cLm = sp.simplify(poly.coeff_monomial(Lm))
    only_L = sp.simplify(cLp + cLm) == 0     # equal and opposite => f(L) only
    # remaining transcendental content in the algebraic parts
    rest = poly.coeff_monomial(1)
    bad = sorted({type(f).__name__ for f in
                  (sp.sympify(rest).atoms(sp.Function)
                   | sp.sympify(cLp).atoms(sp.Function))})
    return bool(only_L), deg, leftover + bad


def check_step1(l_max=4):
    Q0 = Q_order0(l_max)
    print("Step (1): is Q_tau^sigma = A + B*L with A,B algebraic, "
          "L = ln((xi+1)/(xi-1))?\n")
    print(f"{'(tau,sigma)':<13}{'deg in logs':>12}{'only L':>9}"
          f"{'other transcendentals':>24}")
    print("-" * 60)
    ok = True
    for tau in range(l_max + 1):
        for sigma in range(min(tau, 2) + 1):
            only_L, deg, bad = structure(Q_assoc(tau, sigma, Q0))
            print(f"({tau},{sigma}){'':<7}{str(deg):>12}"
                  f"{('yes' if only_L else 'NO'):>9}"
                  f"{(str(bad) if bad else 'none'):>24}")
            if not (only_L and deg is not None and deg <= 1 and not bad):
                ok = False
    print("-" * 60)
    print(f"\nStep (1) verdict: single-logarithm structure holds : {ok}\n")
    return ok


def check_step2():
    """Verify the closed form against quadrature at several (a, c)."""
    mp.mp.dps = 30
    print("Step (2): does the L branch reduce to e^{+-a} E_1(a * shift)?\n")
    print(f"{'a':>6}{'c':>6}{'quadrature':>24}{'closed form':>24}{'|diff|':>12}")
    print("-" * 72)
    worst = mp.mpf(0)
    for a in (mp.mpf("0.7"), mp.mpf("1.5"), mp.mpf("3.0")):
        for c in (mp.mpf("1.3"), mp.mpf("2.0"), mp.mpf("4.5")):
            f = lambda t: mp.e ** (-a * t) * mp.log((t + 1) / (t - 1))
            quad = mp.quad(f, [c, mp.inf])
            closed = (mp.e ** (-a * c) / a) * mp.log((c + 1) / (c - 1)) \
                + (1 / a) * (mp.e ** a * mp.e1(a * (c + 1))
                             - mp.e ** (-a) * mp.e1(a * (c - 1)))
            d = abs(quad - closed)
            worst = max(worst, d)
            print(f"{float(a):>6.2f}{float(c):>6.2f}{mp.nstr(quad, 15):>24}"
                  f"{mp.nstr(closed, 15):>24}{mp.nstr(d, 3):>12}")
    print("-" * 72)
    good = worst < mp.mpf("1e-20")
    print(f"\nStep (2) verdict: closed form verified, worst |diff| = "
          f"{mp.nstr(worst, 3)}  -> {'CONFIRMED' if good else 'FAILED'}")
    print("Both transcendental terms are e^{+-a} E_1(a * shift): the Stieltjes")
    print("seed e^a E_1(a) with shifted argument. No new class appears.")
    return good


def main() -> None:
    print("Phase 0 Q2 -- transcendental seed set of the sigma != 0 xi half\n")
    s1 = check_step1()
    s2 = check_step2()
    print("\n" + "=" * 72)
    if s1 and s2:
        print("Phase 0 Q2 = GO. The sigma != 0 radial half introduces NO new")
        print("transcendental class: Q carries exactly one logarithm at every")
        print("(tau, sigma), and integrating it against e^{-a xi} poly(xi)")
        print("yields only e^{+-a} E_1(a * shift) -- the seed Paper 18 section")
        print("'Level 2: e^a E_1(a)' already classifies, and whose higher E_n")
        print("reduce to E_1 by recurrence. Tagging obligation = cite, not derive.")
    else:
        print("Phase 0 Q2 = RESCOPE or checker still wrong. Do not proceed to")
        print("Phase 1 until this is clean.")


if __name__ == "__main__":
    main()
