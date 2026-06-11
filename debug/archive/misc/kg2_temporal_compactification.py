"""KG-2: Where 2*pi enters the Klein-Gordon spectrum.

Compares (a) open time t in R: spatial heat kernel only, no temporal modes,
no pi in the partial sum;  (b) periodic Euclidean time t ~ t+beta:
omega_{n,k}^2 = n(n+2) + (2 pi k / beta)^2 + m^2 introduces 2pi into the
spectrum at the compactification step.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp
import mpmath as mp

OUT = Path(__file__).parent / "data" / "kg2_temporal_compactification.json"

mp.mp.dps = 50


def open_time_partition(beta: float, m2_int: int, N_max: int) -> dict:
    """Z_open(beta) = sum_n (n+1)^2 exp(-beta * sqrt(n(n+2) + m^2)).

    No temporal modes; beta is just an external Boltzmann parameter.
    Spectrum = bare spatial spectrum, omega_n in Q[sqrt(d)].
    """
    beta_mp = mp.mpf(beta)
    leading_terms = []
    Z = mp.mpf(0)
    for n in range(0, N_max + 1):
        omega2 = n * (n + 2) + m2_int
        omega = mp.sqrt(omega2)
        deg = (n + 1) ** 2
        contribution = deg * mp.exp(-beta_mp * omega)
        Z += contribution
        if n < 10:
            # symbolic counterpart
            omega_sym = sp.sqrt(sp.Rational(n * (n + 2)) + sp.Rational(m2_int))
            leading_terms.append({
                "n": n,
                "degeneracy": deg,
                "omega_sym": str(omega_sym),
                "omega_numeric_50dps": mp.nstr(omega, 30),
                "contribution_50dps": mp.nstr(contribution, 30),
            })
    return {
        "beta": beta,
        "m2_int": m2_int,
        "N_max": N_max,
        "Z_partial_50dps": mp.nstr(Z, 30),
        "leading_terms": leading_terms,
        "transcendental_content": "none in spectrum; only exp(-beta*sqrt(d)) "
                                  "with d square-free integer.",
    }


def compactified_spectrum(beta: sp.Rational, m2_int: int, n_max: int, k_max: int):
    """Periodic Euclidean time: omega_{n,k}^2 = n(n+2) + (2 pi k / beta)^2 + m^2.

    Show the 2*pi appears here and not before.
    """
    beta_sym = beta
    rows = []
    pi_appears_at_step = None
    for n in range(0, n_max + 1):
        for k in range(0, k_max + 1):
            spatial = sp.Rational(n * (n + 2)) + sp.Rational(m2_int)
            temporal = (2 * sp.pi * k / beta_sym) ** 2
            omega2 = spatial + temporal
            contains_pi = sp.pi in omega2.atoms()
            if contains_pi and pi_appears_at_step is None:
                pi_appears_at_step = {"n": n, "k": k, "from": "temporal mode (2 pi k / beta)^2"}
            rows.append({
                "n": n,
                "k": k,
                "spatial_part": str(spatial),
                "temporal_part": str(temporal),
                "omega2_sym": str(omega2),
                "contains_pi": bool(contains_pi),
            })
    return {
        "beta": str(beta_sym),
        "m2_int": m2_int,
        "n_max": n_max,
        "k_max": k_max,
        "pi_appears_at": pi_appears_at_step,
        "rows": rows[:25],  # only need the leading block for the report
        "row_count_total": len(rows),
    }


def first_five_eigenvalues_each_case(beta: float, m2_int: int):
    """Report (a) first 5 omega_n in open-time case (k=0 only),
       (b) first 5 omega_{n,k} in compactified case (mixing n, k)."""
    # Open time: just spatial
    open_list = []
    for n in range(0, 5):
        omega2 = sp.Rational(n * (n + 2)) + sp.Rational(m2_int)
        if omega2 == 0:
            omega = sp.Integer(0)
        else:
            from kg1_algebraic_ring import squarefree_decompose
            c, d = squarefree_decompose(omega2)
            omega = c * sp.sqrt(d)
        open_list.append({
            "n": n,
            "omega_sym": str(omega),
            "omega_numeric": float(omega),
        })

    # Compactified: enumerate (n, k) by ascending omega^2
    beta_sym = sp.Rational(beta).limit_denominator(1000) if isinstance(beta, float) else sp.Rational(beta)
    pairs = []
    for n in range(0, 8):
        for k in range(0, 8):
            spatial = sp.Rational(n * (n + 2)) + sp.Rational(m2_int)
            temporal = (2 * sp.pi * k / beta_sym) ** 2
            omega2 = spatial + temporal
            omega2_num = float(omega2)
            pairs.append((omega2_num, n, k, omega2))
    pairs.sort()
    compact_list = []
    for omega2_num, n, k, omega2_sym in pairs[:5]:
        omega_sym = sp.sqrt(omega2_sym)
        compact_list.append({
            "n": n,
            "k": k,
            "omega2_sym": str(omega2_sym),
            "omega_sym": str(omega_sym),
            "omega_numeric": float(sp.N(omega_sym, 30)),
            "contains_pi": bool(sp.pi in omega2_sym.atoms()),
        })
    return {"open_time_first5": open_list, "compactified_first5": compact_list}


def main():
    open_result = open_time_partition(beta=1.0, m2_int=0, N_max=30)
    compact_result = compactified_spectrum(beta=sp.Integer(1), m2_int=0,
                                           n_max=4, k_max=4)
    first5 = first_five_eigenvalues_each_case(beta=1.0, m2_int=0)

    headline = {
        "open_time_pi_in_spectrum": "no",
        "compactified_pi_in_spectrum": "yes; first appearance: " +
                                       json.dumps(compact_result["pi_appears_at"]),
        "compactification_introduces_2pi_per_temporal_mode": True,
    }

    out = {
        "script": "kg2_temporal_compactification.py",
        "headline": headline,
        "open_time": open_result,
        "compactified": compact_result,
        "first5_each_case": first5,
    }
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps(headline, indent=2))


if __name__ == "__main__":
    main()
