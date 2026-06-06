"""
Sprint Q5'-CH-3 — M3 cyclotomic verification at integer s.

Verify the M3 component of the master Mellin engine sits bit-exactly in
the cyclotomic mixed-Tate ring MT(Z[i, 1/2], 4) at integer s, via the
vertex-parity Dirichlet content (Paper 28 Theorem 3):

    D_even(s) - D_odd(s) = 2^{s-1} . (beta(s) - beta(s-2))

This is the third stone of the cosmic-Galois U* bridge, completing the
M2 verification of Q5'-CH-2 with its M3 counterpart.

Two structural observations from the closed forms:

1. **Odd s:** beta(odd) has the Euler-style closed form
   beta(2k+1) = (-1)^k . E_{2k} / (2 (2k)!) . (pi/2)^{2k+1}
   which is pi^{odd} . Q (pure-Tate M1). The M3 expression at odd s
   therefore COLLAPSES to M1 — no genuinely cyclotomic content.

2. **Even s:** beta(even) has no closed form in pi (Catalan G = beta(2)
   irrationality status conjectural; beta(4), beta(6), ... believed
   irrational by Apery-style arguments but not proven). The M3 content
   at even s is GENUINELY cyclotomic mixed-Tate at level 4.

Methodology
-----------
For each s in {2, 3, 4, 5, 6}:
- Compute D_even(s) and D_odd(s) via the Hurwitz zeta formulas from
  geovac/qed_vertex.py at high precision (mpmath, 80 dps).
- Compute the antisymmetric M3 combination D_even - D_odd numerically.
- Cross-check against 2^{s-1} . (beta(s) - beta(s-2)).
- Express in closed form: at odd s, fully reduce to pi^{odd} . Q (M1);
  at even s, keep beta(even) symbolic and express as M3.

Sympy is used for the closed-form expressions; mpmath gives the high-
precision cross-checks.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import mpmath
import sympy as sp
from sympy import Rational, pi, Catalan, simplify, N as numeric, factorial

from geovac.qed_vertex import _dirac_D_even, _dirac_D_odd


# Closed forms for beta(s):
#   beta(0) = 1/2
#   beta(1) = pi/4
#   beta(2) = G (Catalan, kept symbolic)
#   beta(3) = pi^3 / 32
#   beta(4) = no closed form (kept symbolic as beta_4)
#   beta(5) = 5 pi^5 / 1536
#   beta(6) = no closed form (kept symbolic as beta_6)
#   beta(7) = 61 pi^7 / 184320

# Euler numbers E_{2k} for k = 0, 1, 2, 3: E_0=1, E_2=-1, E_4=5, E_6=-61
# beta(2k+1) = (-1)^k E_{2k} / (2 (2k)!) (pi/2)^{2k+1}

_EULER_NUMBERS = {0: 1, 2: -1, 4: 5, 6: -61, 8: 1385}


def beta_closed_form(s: int):
    """Return sympy closed form for Dirichlet beta(s).

    beta(0) = 1/2 (rational).
    beta(2k+1) for k >= 0: Euler-style pi^{2k+1} . Q.
    beta(2k) for k >= 1: kept symbolic as sympy Symbol("beta_{2k}").
    """
    if s == 0:
        return Rational(1, 2)
    if s == 1:
        return pi / 4
    if s == 2:
        return Catalan  # sympy.Catalan = G
    if s % 2 == 1:
        # odd s = 2k + 1, k = (s - 1) // 2
        k = (s - 1) // 2
        if 2 * k not in _EULER_NUMBERS:
            return sp.Symbol(f"beta_{s}")
        E_2k = _EULER_NUMBERS[2 * k]
        sign = Rational(-1) ** k
        return sign * Rational(E_2k) / (2 * factorial(2 * k)) * (pi / 2) ** (2 * k + 1)
    # even s >= 2 and s != 2: keep as opaque symbol (beta(s) for s = 4, 6, ...)
    # Plain Symbol (no positive=True) so the classifier matches it correctly.
    return sp.Symbol(f"beta_{s}")


def D_diff_closed_form(s: int):
    """D_even(s) - D_odd(s) = 2^{s-1} (beta(s) - beta(s-2)) closed form."""
    return 2 ** (s - 1) * (beta_closed_form(s) - beta_closed_form(s - 2))


def classify_M3_content(expr) -> Dict:
    """Classify expression as M1-only (pure pi^{odd|even}) or genuinely M3.

    M1: expression contains only pi powers + rationals (no beta_{even>=4}, no Catalan).
    M3: expression contains Catalan or beta_{2k} for k >= 1.
    """
    expr_e = sp.expand(expr)
    has_catalan = expr_e.has(Catalan)
    has_even_beta = any(
        expr_e.has(sp.Symbol(f"beta_{2 * k}"))
        for k in range(2, 6)  # check beta_4, beta_6, beta_8, beta_10
    )
    in_M3_genuine = has_catalan or has_even_beta

    return {
        "expr_simplified": str(sp.simplify(expr_e)),
        "in_M1_only": not in_M3_genuine,
        "in_M3_genuine": in_M3_genuine,
        "contains_Catalan": has_catalan,
        "contains_higher_beta_even": has_even_beta,
    }


def _beta_mpmath(s_val: int, dps: int = 80) -> mpmath.mpf:
    """Dirichlet beta(s) at high precision via the Hurwitz form
    beta(s) = 4^{-s} (zeta(s, 1/4) - zeta(s, 3/4))
    (NOT 2^{-s} — beta sums 1/(2n+1)^s starting from n=0, so the natural
    Hurwitz shift is at 1/4 and 3/4, with prefactor 4^{-s}).

    For s = 0, returns 1/2 directly (analytic continuation).
    """
    if s_val == 0:
        return mpmath.mpf("0.5")
    quarter = mpmath.mpf(1) / 4
    three_quarter = mpmath.mpf(3) / 4
    h1 = mpmath.zeta(s_val, quarter)
    h2 = mpmath.zeta(s_val, three_quarter)
    return mpmath.power(4, -s_val) * (h1 - h2)


def cross_check_numeric(s_val: int, dps: int = 80) -> Dict:
    """Cross-check via mpmath at high precision."""
    mpmath.mp.dps = dps

    D_even = _dirac_D_even(s_val)
    D_odd = _dirac_D_odd(s_val)
    D_diff = D_even - D_odd

    s_int = int(s_val)
    beta_s = _beta_mpmath(s_int)
    beta_s_minus_2 = _beta_mpmath(s_int - 2)

    closed_form_predict = mpmath.power(2, s_int - 1) * (beta_s - beta_s_minus_2)
    residual = abs(D_diff - closed_form_predict)

    return {
        "D_even_numeric": mpmath.nstr(D_even, 50),
        "D_odd_numeric": mpmath.nstr(D_odd, 50),
        "D_diff_numeric": mpmath.nstr(D_diff, 50),
        "closed_form_prediction_2^{s-1}(beta(s) - beta(s-2))": mpmath.nstr(
            closed_form_predict, 50
        ),
        "residual_abs": float(residual),
        "residual_below_1e-30": float(residual) < 1e-30,
    }


def main() -> None:
    output: Dict = {"sprint": "Q5'-CH-3", "panel": {}}

    s_values = [2, 3, 4, 5, 6]
    print("=== M3 verification at integer s via D_even(s) - D_odd(s) ===\n")

    odd_s_collapse_to_M1 = 0
    even_s_genuine_M3 = 0

    for s_val in s_values:
        print(f"--- s = {s_val} ---")

        # Closed form
        D_diff_cf = D_diff_closed_form(s_val)
        D_diff_simplified = sp.simplify(sp.expand(D_diff_cf))
        print(f"  D_even - D_odd = 2^{s_val-1} (beta({s_val}) - beta({s_val-2}))")
        print(f"                 = {D_diff_simplified}")

        # M3 classification
        classification = classify_M3_content(D_diff_simplified)
        if classification["in_M1_only"]:
            print(f"  Classification: M1 collapse (pi^... . Q only)")
            odd_s_collapse_to_M1 += 1
        else:
            print(f"  Classification: GENUINELY M3 (cyclotomic content)")
            if classification["contains_Catalan"]:
                print(f"    Contains Catalan G = beta(2)")
            if classification["contains_higher_beta_even"]:
                print(f"    Contains beta_{{2k}} for k >= 2 (no closed form)")
            even_s_genuine_M3 += 1

        # Numeric cross-check
        try:
            numeric_check = cross_check_numeric(s_val)
            print(f"  Numeric cross-check (80 dps):")
            print(f"    D_even - D_odd       = {numeric_check['D_diff_numeric']}")
            print(f"    2^(s-1)(beta(s)-beta(s-2)) = "
                  f"{numeric_check['closed_form_prediction_2^{s-1}(beta(s) - beta(s-2))']}")
            print(f"    residual = {numeric_check['residual_abs']:.3e} "
                  f"(below 1e-30: {numeric_check['residual_below_1e-30']})")
        except Exception as e:
            numeric_check = {"error": str(e)}
            print(f"  Numeric cross-check FAILED: {e}")

        output["panel"][str(s_val)] = {
            "closed_form": str(D_diff_simplified),
            "classification": classification,
            "numeric_check": numeric_check,
        }
        print()

    # Summary
    print("=" * 60)
    print(f"M1-collapse at odd s:   {odd_s_collapse_to_M1}/{len([s for s in s_values if s % 2 == 1])}")
    print(f"Genuine M3 at even s:   {even_s_genuine_M3}/{len([s for s in s_values if s % 2 == 0])}")
    print("=" * 60)

    output["summary"] = {
        "odd_s_M1_collapse": odd_s_collapse_to_M1,
        "even_s_genuine_M3": even_s_genuine_M3,
        "structural_finding": (
            "At odd s, beta(odd) reduces to pi^odd . Q via Euler-style "
            "closed form, so the M3 antisymmetric combination collapses to "
            "M1. At even s, beta(even) has no closed form in pi, so the M3 "
            "content is genuinely cyclotomic mixed-Tate at level 4."
        ),
    }

    out_path = Path("debug/data/sprint_q5p_ch3_data.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
    print(f"\nOutput written: {out_path}")


if __name__ == "__main__":
    main()
