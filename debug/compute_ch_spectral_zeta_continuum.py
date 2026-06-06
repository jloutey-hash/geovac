"""
Sprint Q5'-CH-2 — Continuum CH spectral zeta and M2 pure-Tate verification.

Goal
----
Pick up from Sprint Q5'-CH-1 (chirality-parity factorisation of the master
Mellin k-slot at finite n_max) and push to the continuum: evaluate the
Camporesi--Higuchi spectral zeta

    zeta_{D²}(s) = sum_n 2 n (n+1) (n + 1/2)^{-2s}

at integer s in {1, 2, 3, 4, 5}, verify bit-exact closed forms in the
M2 pure-Tate ring oplus_k pi^{2k} . Q (Paper 32 §VIII Cor cor:m2_mixed_tate),
and compare with truncated values at n_max in {2, 3, 4} from Sprint
Q5'-CH-1's moment data.

This is the period-data side of the Q5' bridge: M1, M2, M3 are
characterised by their period rings. Q5'-CH-1 located the chirality-parity
component of the k-slot at finite cutoff. Q5'-CH-2 locates the
heat-kernel-order component (M2) at the continuum, in bit-exact closed
form, in the M2 ring.

Methodology
-----------
The continuum CH spectral zeta reduces via the substitution m = n + 1/2:

    zeta_{D²}(s) = sum_{m = 3/2, 5/2, ...} 2 (m² - 1/4) m^{-2s}
                 = 2 . zeta(2s-2, 3/2) - (1/2) . zeta(2s, 3/2)

where zeta(s, a) is the Hurwitz zeta function. The Hurwitz identity
zeta(s, 1/2) = (2^s - 1) zeta(s) plus the shift
zeta(s, 3/2) = zeta(s, 1/2) - 2^s reduces every integer-s value to a
Q-linear combination of integer-s Riemann zeta values, hence (at even s)
to pi^{2k} . Q.

Sympy is used in exact rational arithmetic with sp.zeta() for the Riemann
values. Closed forms are extracted symbolically and verified at high
floating precision.

Output
------
debug/data/sprint_q5p_ch2_data.json with continuum closed forms,
truncated values, and the M2-ring fingerprint per integer s.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import sympy as sp
from sympy import Rational, Symbol, zeta, pi, simplify, nsimplify, N as numeric


def continuum_zeta_D2(s_val) -> sp.Expr:
    """Continuum CH spectral zeta at integer s.

    zeta_{D²}(s) = 2 zeta(2s-2, 3/2) - (1/2) zeta(2s, 3/2)
                 = 2 [(2^{2s-2} - 1) zeta(2s-2) - 2^{2s-2}]
                   - (1/2) [(2^{2s} - 1) zeta(2s) - 2^{2s}]

    Uses Hurwitz identity zeta(s, 1/2) = (2^s - 1) zeta(s), zeta(s, 3/2) =
    zeta(s, 1/2) - 2^s. Returns exact sympy expression in pi^{2k} . Q for
    integer s >= 1 with both arguments >= 0 (i.e., s >= 1).
    """
    s = Rational(s_val)
    two_s_minus_2 = 2 * s - 2
    two_s = 2 * s

    # zeta(2s-2, 3/2):
    # If 2s-2 = 0: this is the "constant" sum, divergent; the s=1 case has 2s-2=0
    # We handle the s=1 case by regularisation: zeta(0, 3/2) = 1/2 - 3/2 = -1
    # via the Hurwitz identity zeta(0, a) = 1/2 - a.
    if two_s_minus_2 == 0:
        hurwitz_low = Rational(1, 2) - Rational(3, 2)  # zeta(0, 3/2) = -1
    elif two_s_minus_2 < 0:
        # Hurwitz zeta at negative integer is a Bernoulli polynomial; handle via sympy
        hurwitz_low = sp.zeta(two_s_minus_2, Rational(3, 2))
    else:
        # Hurwitz identity:
        # zeta(s, 1/2) = (2^s - 1) zeta(s);  zeta(s, 3/2) = zeta(s, 1/2) - 2^s
        z_half = (2 ** two_s_minus_2 - 1) * sp.zeta(two_s_minus_2)
        hurwitz_low = z_half - 2 ** two_s_minus_2

    # zeta(2s, 3/2): always positive at s >= 1
    z_half_high = (2 ** two_s - 1) * sp.zeta(two_s)
    hurwitz_high = z_half_high - 2 ** two_s

    expr = 2 * hurwitz_low - Rational(1, 2) * hurwitz_high
    return sp.simplify(expr)


def truncated_zeta_D2(n_max: int, s_val) -> sp.Expr:
    """Truncated CH spectral zeta at finite n_max."""
    s = Rational(s_val)
    total = Rational(0)
    for n in range(1, n_max + 1):
        m = Rational(2 * n + 1, 2)  # n + 1/2
        total += Rational(2 * n * (n + 1)) * m ** (-2 * s)
    return sp.simplify(total)


def is_in_M2_ring(expr: sp.Expr, max_power: int = 10) -> Dict:
    """Check whether expr is a Q-linear combination of {1, pi^2, pi^4, ...}.

    Returns dict with:
        in_M2_ring : bool
        decomposition : list of (coefficient, power) tuples
        residual : sympy expression (zero if exactly in M2 ring)
    """
    expr_e = sp.expand(expr)
    expr_collected = sp.collect(expr_e, sp.pi)

    decomp = []
    residual = expr_collected
    for k in range(max_power + 1):
        coeff = residual.coeff(sp.pi, 2 * k)
        if coeff != 0:
            decomp.append((coeff, 2 * k))
            residual = residual - coeff * sp.pi ** (2 * k)

    residual = sp.simplify(residual)
    in_ring = (residual == 0)

    return {
        "in_M2_ring": in_ring,
        "decomposition": [(str(c), p) for c, p in decomp],
        "residual": str(residual),
    }


def compute_panel():
    """Compute the s-by-n_max panel of zeta_{D²}(s)."""
    output: Dict = {"sprint": "Q5'-CH-2", "panel": []}

    s_values = [1, 2, 3, 4, 5]
    n_max_values = [2, 3, 4]

    # Continuum closed forms
    print("=== Continuum CH spectral zeta zeta_{D²}(s) ===")
    continuum_table: Dict = {}
    for s_val in s_values:
        z_cont = continuum_zeta_D2(s_val)
        print(f"\n  s = {s_val}:")
        print(f"    zeta_{{D²}}({s_val}) = {z_cont}")
        m2_check = is_in_M2_ring(z_cont)
        print(f"    In M2 ring oplus pi^(2k) . Q ? {m2_check['in_M2_ring']}")
        for coeff, p in m2_check["decomposition"]:
            print(f"      coeff of pi^{p}: {coeff}")
        if not m2_check["in_M2_ring"]:
            print(f"      residual: {m2_check['residual']}")
        z_num = numeric(z_cont, 30)
        continuum_table[str(s_val)] = {
            "closed_form": str(z_cont),
            "numeric_30dps": str(z_num),
            "M2_ring_check": m2_check,
        }
    output["continuum"] = continuum_table

    # Truncated values + convergence
    print("\n=== Truncated values and convergence to continuum ===")
    truncated_table: Dict = {}
    for s_val in s_values:
        print(f"\n  s = {s_val}:")
        z_cont = continuum_zeta_D2(s_val)
        z_cont_num = float(numeric(z_cont, 30))
        per_n_max: Dict = {}
        for n_max in n_max_values:
            z_trunc = truncated_zeta_D2(n_max, s_val)
            z_trunc_num = float(numeric(z_trunc, 30))
            tail = float(numeric(z_cont - z_trunc, 30))
            print(
                f"    n_max={n_max}: zeta_trunc = {z_trunc} "
                f"(~ {z_trunc_num:.10f})"
            )
            print(f"      tail = zeta_cont - zeta_trunc = {tail:.3e}")
            per_n_max[str(n_max)] = {
                "truncated_rational": str(z_trunc),
                "numeric": z_trunc_num,
                "tail_to_continuum": tail,
            }
        truncated_table[str(s_val)] = {
            "continuum_numeric": z_cont_num,
            "per_n_max": per_n_max,
        }
    output["truncated"] = truncated_table

    return output


def main() -> None:
    output = compute_panel()

    out_path = Path("debug/data/sprint_q5p_ch2_data.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
    print(f"\nOutput written: {out_path}")

    # Final summary
    print("\n=== Summary ===")
    in_ring_count = 0
    total = 0
    for s_val in [1, 2, 3, 4, 5]:
        check = output["continuum"][str(s_val)]["M2_ring_check"]
        in_ring_count += int(check["in_M2_ring"])
        total += 1
        marker = "OK" if check["in_M2_ring"] else "FAIL"
        cf = output["continuum"][str(s_val)]["closed_form"]
        print(f"  s={s_val}: zeta_{{D²}}({s_val}) = {cf}  [M2 ring: {marker}]")
    print(f"\nIn M2 ring oplus pi^(2k) . Q: {in_ring_count}/{total}")


if __name__ == "__main__":
    main()
