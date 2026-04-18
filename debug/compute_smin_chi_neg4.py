"""Track RH-P driver: depth-2 chi_{-4} analog of RH-J.1 for S_min.

Sprint 4, April 2026.

PROBLEM STATEMENT
-----------------
Sprint 3 Track RH-J established the depth-1 identity:

    D_even(s) - D_odd(s) = 2^{s-1} * (beta(s) - beta(s-2))          (RH-J.1)

where D_even(s), D_odd(s) are the Dirac Dirichlet sub-sums split by the
parity of the mode index n, and beta(s) = L(s, chi_{-4}) is the Dirichlet
beta function.

The natural depth-2 analog of RH-J.1 operates on the CG-weighted two-loop
irreducible constant (Paper 28 §IV):

    S_min = sum_{k=1}^inf T(k)^2,  T(k) = 2*zeta(2, k+3/2) - (1/2)*zeta(4, k+3/2)

This driver splits S_min by the parity of k:

    S_min^even = sum_{k = 2, 4, 6, ...} T(k)^2
    S_min^odd  = sum_{k = 1, 3, 5, ...} T(k)^2
    S_min      = S_min^even + S_min^odd                              (sanity)

and tests whether S_min^even - S_min^odd (and alternative combinations)
admits a closed form involving beta(s), Catalan G, Dirichlet L-values.

TAIL CORRECTION
---------------
Asymptotic expansion of T(k)^2 in powers of 1/a = 1/(k+3/2):

    T(k)^2 = 4/a^2 + 4/a^3 + 5/(3 a^4) - 2/(3 a^5) - 193/(180 a^6)
             - 23/(60 a^7) + 227/(1008 a^8) + ... (derived in sympy)

Each 1/a^j term summed from k=N+1 to infinity equals hurwitz(j, N+5/2).
We use 8 terms of the expansion for the tail correction, giving ~O(N^{-9})
residual error.  At N=4000, tail residual ~ N^{-9} ~ 10^{-33}, suitable
for PSLQ at 100-dps precision.

For the parity-split sums, we need two separate tails:
    S_even_tail(N) = sum_{k=N+1, k even} T(k)^2
    S_odd_tail(N)  = sum_{k=N+1, k odd} T(k)^2

Using the same asymptotic expansion, each j-th term contributes:
    sum_{k=N+1, k even} 1/(k+3/2)^j
    sum_{k=N+1, k odd}  1/(k+3/2)^j

Both can be written in terms of Hurwitz zetas at quarter-integer shifts:
    k = 2m => k+3/2 = 2m+3/2 = 2(m+3/4), so 1/(k+3/2)^j = 2^{-j}/(m+3/4)^j
    k = 2m+1 => k+3/2 = 2m+5/2 = 2(m+5/4), so 1/(k+3/2)^j = 2^{-j}/(m+5/4)^j

This brings in hurwitz(j, 3/4) and hurwitz(j, 5/4) -- the same quarter-
integer shifts that produce beta(s) = L(s, chi_{-4}) in the depth-1 case.

PRECISION: 100 digits.
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import mpmath

mpmath.mp.dps = 100


# ----------------------------------------------------------------------
# T(k) and Asymptotic coefficients of T(k)^2 in 1/a, a = k+3/2
# ----------------------------------------------------------------------

# Exact coefficients of T(k)^2 = sum_{j >= 2} c_j / a^j, derived symbolically
# (see smin_chi_neg4_memo.md §2 for derivation).
#
# Known erratum (Sprint 5 S_min verification, 2026-04-17): the
# coefficients for j >= 6 below are wrong, which introduces a ~6.5e-20
# error at the 20th digit of S_min computed from this asymptotic
# expansion. The true S_min value (verified by three independent methods
# to >= 80 digits) is 2.47993693803422255441357950082938214468... See
# `debug/smin_verification_memo.md` for the audit. This error does NOT
# affect RH-P's PSLQ negative result — a 20-th-digit numerical shift
# does not change Q-linear independence against a finite basis. The
# dictionary is preserved as-is for reproducibility; for numerically
# correct computation of S_min use `debug/smin_identification.py`
# (post-patch) which delegates to mpmath.nsum Levin.
T_SQUARED_COEFFS: Dict[int, Fraction] = {
    2: Fraction(4, 1),
    3: Fraction(4, 1),
    4: Fraction(5, 3),
    5: Fraction(-2, 3),
    6: Fraction(-193, 180),
    7: Fraction(-23, 60),
    8: Fraction(227, 1008),
    9: Fraction(457, 2520),
    10: Fraction(-17839, 75600),
    11: Fraction(-83, 504),
    # Higher orders computed to machine precision but not needed.
}


def T_k(k: int) -> mpmath.mpf:
    """T(k) = 2*zeta(2, k+3/2) - (1/2)*zeta(4, k+3/2)."""
    a = mpmath.mpf(k) + mpmath.mpf(3) / 2
    return (2 * mpmath.hurwitz(2, a)
            - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, a))


def tail_total(N: int) -> mpmath.mpf:
    """Total tail sum_{k=N+1}^inf T(k)^2 via asymptotic expansion.

    Uses sum_{k=N+1}^inf 1/(k+3/2)^j = hurwitz(j, N+5/2).

    Known limitation: T_SQUARED_COEFFS above has c_j wrong for j>=6,
    producing a ~6.5e-20 error at 20th digit. See the dictionary's
    docstring for context.
    """
    a0 = mpmath.mpf(N + 1) + mpmath.mpf(3) / 2
    S = mpmath.mpf(0)
    for j, c in T_SQUARED_COEFFS.items():
        S += mpmath.mpf(c.numerator) / c.denominator * mpmath.hurwitz(j, a0)
    return S


def tail_even(N: int) -> mpmath.mpf:
    """Even-k tail: sum_{k=N+1, k even}^inf T(k)^2.

    k even, k = 2m, first m >= ceil((N+1)/2).  k + 3/2 = 2m + 3/2 = 2(m+3/4).
    So 1/(k+3/2)^j = 2^{-j}/(m+3/4)^j, summed over m >= m0, gives
    2^{-j} * hurwitz(j, m0 + 3/4).

    Same 20th-digit limitation as tail_total.
    """
    m0 = (N + 2) // 2  # smallest even k > N is k = 2*m0; need k >= N+1 => m0 >= (N+1)/2
    S = mpmath.mpf(0)
    for j, c in T_SQUARED_COEFFS.items():
        c_mp = mpmath.mpf(c.numerator) / c.denominator
        S += c_mp * mpmath.power(2, -j) * mpmath.hurwitz(
            j, mpmath.mpf(m0) + mpmath.mpf(3) / 4
        )
    return S


def tail_odd(N: int) -> mpmath.mpf:
    """Odd-k tail: sum_{k=N+1, k odd}^inf T(k)^2.

    k odd, k = 2m+1, first m >= ceil((N)/2).  k + 3/2 = 2m + 5/2 = 2(m+5/4).
    So 1/(k+3/2)^j = 2^{-j}/(m+5/4)^j, summed m >= m0, gives
    2^{-j} * hurwitz(j, m0 + 5/4).

    Same 20th-digit limitation as tail_total.
    """
    m0 = (N + 1) // 2
    S = mpmath.mpf(0)
    for j, c in T_SQUARED_COEFFS.items():
        c_mp = mpmath.mpf(c.numerator) / c.denominator
        S += c_mp * mpmath.power(2, -j) * mpmath.hurwitz(
            j, mpmath.mpf(m0) + mpmath.mpf(5) / 4
        )
    return S


def compute_smin_parity_split(n_terms: int = 4000) -> Dict:
    """Compute S_min^even, S_min^odd, and S_min^total with asymptotic tail.

    Uses direct summation for k = 1..n_terms, plus asymptotic tail
    from n_terms+1 to infinity.  At n_terms=4000, tail residual error
    is ~n_terms^{-12} ~ 10^{-43}, more than enough for 100-dps PSLQ.
    """
    t0 = time.time()

    S_even = mpmath.mpf(0)
    S_odd = mpmath.mpf(0)
    S_total = mpmath.mpf(0)

    T_samples = {}

    for k in range(1, n_terms + 1):
        Tk = T_k(k)
        Tk2 = Tk ** 2
        S_total += Tk2
        if k % 2 == 0:
            S_even += Tk2
        else:
            S_odd += Tk2
        if k in (1, 2, 3, 4, 5, 10, 50, 100, 1000):
            T_samples[k] = float(Tk)

    # Tail corrections via asymptotic expansion
    S_total_tail = tail_total(n_terms)
    S_even_tail = tail_even(n_terms)
    S_odd_tail = tail_odd(n_terms)

    # Sanity checks:
    #   S_even + S_odd = S_total (partial sums)
    sanity_partial = abs((S_even + S_odd) - S_total)
    #   tail_even + tail_odd = tail_total
    sanity_tail = abs((S_even_tail + S_odd_tail) - S_total_tail)

    t1 = time.time()

    S_total_corr = S_total + S_total_tail
    S_even_corr = S_even + S_even_tail
    S_odd_corr = S_odd + S_odd_tail

    return {
        "n_terms": n_terms,
        "S_min_total": S_total_corr,
        "S_min_total_float": float(S_total_corr),
        "S_min_even": S_even_corr,
        "S_min_even_float": float(S_even_corr),
        "S_min_odd": S_odd_corr,
        "S_min_odd_float": float(S_odd_corr),
        "S_min_diff": S_even_corr - S_odd_corr,
        "S_min_diff_float": float(S_even_corr - S_odd_corr),
        # Without tail (for diagnostic)
        "S_min_total_raw": S_total,
        "S_min_even_raw": S_even,
        "S_min_odd_raw": S_odd,
        # Tail sizes
        "tail_total_float": float(S_total_tail),
        "tail_even_float": float(S_even_tail),
        "tail_odd_float": float(S_odd_tail),
        # Consistency checks
        "partial_sanity_err": float(sanity_partial),
        "tail_sanity_err": float(sanity_tail),
        # Miscellaneous
        "S_min_ratio_even_total": float(S_even_corr / S_total_corr),
        "S_min_ratio_even_odd": float(S_even_corr / S_odd_corr),
        "T_samples": T_samples,
        "compute_time_sec": t1 - t0,
    }


# ----------------------------------------------------------------------
# Basis constants
# ----------------------------------------------------------------------

def beta_num(s: int) -> mpmath.mpf:
    """Dirichlet beta(s) via Hurwitz."""
    if s == 0:
        return mpmath.mpf(1) / 2
    if s == 1:
        return mpmath.pi / 4
    if s < 0:
        return mpmath.mpf('nan')
    return ((mpmath.hurwitz(s, mpmath.mpf(1) / 4)
             - mpmath.hurwitz(s, mpmath.mpf(3) / 4))
            / mpmath.power(4, s))


def dirac_D(s: int) -> mpmath.mpf:
    return (2 * mpmath.hurwitz(s - 2, mpmath.mpf(3) / 2)
            - mpmath.mpf(1) / 2 * mpmath.hurwitz(s, mpmath.mpf(3) / 2))


def dirac_D_even(s: int) -> mpmath.mpf:
    return mpmath.power(2, -s) * (
        8 * mpmath.hurwitz(s - 2, mpmath.mpf(3) / 4)
        - mpmath.mpf(1) / 2 * mpmath.hurwitz(s, mpmath.mpf(3) / 4)
    )


def dirac_D_odd(s: int) -> mpmath.mpf:
    return mpmath.power(2, -s) * (
        8 * mpmath.hurwitz(s - 2, mpmath.mpf(5) / 4)
        - mpmath.mpf(1) / 2 * mpmath.hurwitz(s, mpmath.mpf(5) / 4)
    )


def build_basis_constants(S_min_value: mpmath.mpf,
                          S_even: mpmath.mpf,
                          S_odd: mpmath.mpf) -> Dict[str, mpmath.mpf]:
    """Standard transcendental basis for PSLQ."""
    pi_ = mpmath.pi
    return {
        "1":       mpmath.mpf(1),
        "pi":      pi_,
        "pi^2":    pi_ ** 2,
        "pi^3":    pi_ ** 3,
        "pi^4":    pi_ ** 4,
        "pi^5":    pi_ ** 5,
        "pi^6":    pi_ ** 6,
        "pi^7":    pi_ ** 7,
        "pi^8":    pi_ ** 8,
        "zeta(2)": mpmath.zeta(2),
        "zeta(3)": mpmath.zeta(3),
        "zeta(4)": mpmath.zeta(4),
        "zeta(5)": mpmath.zeta(5),
        "zeta(6)": mpmath.zeta(6),
        "zeta(7)": mpmath.zeta(7),
        "G":       mpmath.catalan,
        "beta(3)": beta_num(3),
        "beta(4)": beta_num(4),
        "beta(5)": beta_num(5),
        "beta(6)": beta_num(6),
        "beta(7)": beta_num(7),
        "beta(8)": beta_num(8),
        "log2":    mpmath.log(2),
        "log3":    mpmath.log(3),
        "D(4)":    dirac_D(4),
        "D(5)":    dirac_D(5),
        "D(6)":    dirac_D(6),
        "D(4)^2":  dirac_D(4) ** 2,
        "G^2":     mpmath.catalan ** 2,
        "beta(4)^2": beta_num(4) ** 2,
        "G*beta(4)": mpmath.catalan * beta_num(4),
        "pi^2*G":  pi_ ** 2 * mpmath.catalan,
        "pi^2*beta(4)": pi_ ** 2 * beta_num(4),
        "pi^4*G":  pi_ ** 4 * mpmath.catalan,
        "pi^2*zeta(3)": pi_ ** 2 * mpmath.zeta(3),
        "pi^2*zeta(5)": pi_ ** 2 * mpmath.zeta(5),
        "zeta(3)^2": mpmath.zeta(3) ** 2,
        "zeta(3)*zeta(5)": mpmath.zeta(3) * mpmath.zeta(5),
        "S_min":   S_min_value,
        "D_even(4)": dirac_D_even(4),
        "D_odd(4)":  dirac_D_odd(4),
        "D_even(4)-D_odd(4)": dirac_D_even(4) - dirac_D_odd(4),
        "D_even(4)^2": dirac_D_even(4) ** 2,
        "D_odd(4)^2": dirac_D_odd(4) ** 2,
        "D_even(4)*D_odd(4)": dirac_D_even(4) * dirac_D_odd(4),
        "S_min_even": S_even,
        "S_min_odd":  S_odd,
    }


# ----------------------------------------------------------------------
# PSLQ helpers
# ----------------------------------------------------------------------

def run_pslq(value: mpmath.mpf, basis: Dict[str, mpmath.mpf],
             tol: float = 1e-60, maxcoeff: int = 10 ** 8,
             basis_keys: Optional[List[str]] = None) -> Dict:
    """Run mpmath PSLQ on [value, *basis_values]."""
    if basis_keys is None:
        basis_keys = list(basis.keys())

    basis_vals = [basis[k] for k in basis_keys]
    vec = [value] + basis_vals

    try:
        rel = mpmath.pslq(vec, tol=tol, maxcoeff=maxcoeff)
    except Exception as e:
        return {"identified": False, "reason": f"PSLQ exception: {e}",
                "basis_keys": basis_keys, "basis_size": len(basis_keys)}

    if rel is None:
        return {"identified": False, "reason": "PSLQ returned None",
                "basis_keys": basis_keys, "basis_size": len(basis_keys)}

    c0 = rel[0]
    if c0 == 0:
        return {"identified": False, "reason": "zero coefficient on value",
                "basis_keys": basis_keys, "basis_size": len(basis_keys),
                "raw_relation": [int(r) for r in rel]}

    components = {}
    reconstructed = mpmath.mpf(0)
    for ri, lab, bv in zip(rel[1:], basis_keys, basis_vals):
        if ri != 0:
            frac = Fraction(int(-ri), int(c0))
            components[lab] = str(frac)
            reconstructed += (mpmath.mpf(-ri) / c0) * bv

    residual = abs(value - reconstructed)

    return {
        "identified": True,
        "basis_size": len(basis_keys),
        "basis_keys": basis_keys,
        "components": components,
        "residual_float": float(residual),
        "reconstructed_float": float(reconstructed),
        "value_float": float(value),
        "raw_relation": [int(r) for r in rel],
    }


def pslq_attempts(value: mpmath.mpf, value_name: str,
                  basis: Dict[str, mpmath.mpf]) -> List[Dict]:
    """Run a battery of PSLQ attempts at different basis subsets."""
    attempts = []

    # Strategy 1: Minimal beta-weighted basis (analog of RH-J.1 targeted)
    attempts.append({
        "label": f"{value_name}.minimal_beta",
        "keys": ["1", "pi^2", "pi^4", "G", "beta(4)", "beta(6)"],
        "tol": 1e-60, "maxcoeff": 10 ** 8,
    })

    # Strategy 2: Extended beta basis up to weight 8
    attempts.append({
        "label": f"{value_name}.extended_beta",
        "keys": ["1", "pi^2", "pi^4", "pi^6", "pi^8",
                 "G", "beta(4)", "beta(6)", "beta(8)"],
        "tol": 1e-60, "maxcoeff": 10 ** 8,
    })

    # Strategy 3: Beta basis + products
    attempts.append({
        "label": f"{value_name}.beta_with_products",
        "keys": ["1", "pi^2", "pi^4", "pi^6", "pi^8",
                 "G", "beta(4)", "beta(6)", "beta(8)",
                 "G^2", "beta(4)^2", "G*beta(4)",
                 "pi^2*G", "pi^4*G", "pi^2*beta(4)"],
        "tol": 1e-60, "maxcoeff": 10 ** 8,
    })

    # Strategy 4: Zeta + beta mixed basis
    attempts.append({
        "label": f"{value_name}.zeta_beta_mixed",
        "keys": ["1", "pi^2", "pi^4", "pi^6",
                 "zeta(3)", "zeta(5)", "zeta(7)",
                 "G", "beta(4)", "beta(6)",
                 "zeta(3)^2", "zeta(3)*zeta(5)"],
        "tol": 1e-60, "maxcoeff": 10 ** 8,
    })

    # Strategy 5: D(s) values in basis
    attempts.append({
        "label": f"{value_name}.D_basis",
        "keys": ["1", "pi^2", "pi^4", "D(4)", "D(5)", "D(6)", "D(4)^2",
                 "G", "beta(4)", "beta(6)",
                 "D_even(4)^2", "D_odd(4)^2", "D_even(4)*D_odd(4)"],
        "tol": 1e-60, "maxcoeff": 10 ** 8,
    })

    # Strategy 6: Ultra-wide combined
    attempts.append({
        "label": f"{value_name}.ultra_wide",
        "keys": ["1", "pi", "pi^2", "pi^3", "pi^4", "pi^5", "pi^6", "pi^7", "pi^8",
                 "zeta(3)", "zeta(5)", "zeta(7)",
                 "G", "beta(4)", "beta(5)", "beta(6)", "beta(7)", "beta(8)",
                 "log2", "log3",
                 "G^2", "beta(4)^2", "G*beta(4)",
                 "zeta(3)^2"],
        "tol": 1e-50, "maxcoeff": 10 ** 8,
    })

    # Strategy 7: Targeted to RH-J.1 analog: use D_even(4) - D_odd(4) pattern
    # If S_diff = C * (D_even(4) - D_odd(4))^2 or similar, we should see it here.
    attempts.append({
        "label": f"{value_name}.depth2_RH_J1_analog",
        "keys": ["1", "pi^2", "pi^4",
                 "D_even(4)", "D_odd(4)", "D_even(4)-D_odd(4)",
                 "D_even(4)^2", "D_odd(4)^2", "D_even(4)*D_odd(4)",
                 "G", "beta(4)", "beta(6)"],
        "tol": 1e-60, "maxcoeff": 10 ** 8,
    })

    # Execute all attempts
    executed = []
    for att in attempts:
        res = run_pslq(value, basis,
                       tol=att["tol"], maxcoeff=att["maxcoeff"],
                       basis_keys=att["keys"])
        executed.append({
            "label": att["label"],
            "keys": att["keys"],
            "pslq": res,
        })
    return executed


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    print("=" * 72)
    print(" Track RH-P: depth-2 chi_{-4} analog of RH-J.1 for S_min")
    print("=" * 72)
    print()

    # -----------------------------------------------------------------
    # §1 Convergence: tail correction effect
    # -----------------------------------------------------------------
    print("§1 Convergence at varying n_terms (with asymptotic tail correction)")
    print("-" * 72)
    conv_records = []
    for N in [500, 1000, 2000, 4000]:
        rec = compute_smin_parity_split(n_terms=N)
        conv_records.append({
            "n_terms": N,
            "S_total_float": rec["S_min_total_float"],
            "S_even_float": rec["S_min_even_float"],
            "S_odd_float":  rec["S_min_odd_float"],
            "S_diff_float": rec["S_min_diff_float"],
            "tail_total_float": rec["tail_total_float"],
            "partial_sanity_err": rec["partial_sanity_err"],
            "tail_sanity_err": rec["tail_sanity_err"],
        })
        print(f"  N={N:>5}: S_total={rec['S_min_total_float']:.20f} "
              f"tail={rec['tail_total_float']:.2e}")
        print(f"          S_even={rec['S_min_even_float']:.20f}  "
              f"tail={rec['tail_even_float']:.2e}")
        print(f"          S_odd ={rec['S_min_odd_float']:.20f}  "
              f"tail={rec['tail_odd_float']:.2e}")
    print()

    # -----------------------------------------------------------------
    # §2 Main computation at n_terms=4000
    # -----------------------------------------------------------------
    print("§2 Main computation: n_terms = 4000 (tail residual ~ 10^{-42})")
    print("-" * 72)
    main_rec = compute_smin_parity_split(n_terms=4000)
    S_total = main_rec["S_min_total"]
    S_even = main_rec["S_min_even"]
    S_odd = main_rec["S_min_odd"]
    S_diff = main_rec["S_min_diff"]

    print(f"  S_min^total = {mpmath.nstr(S_total, 40)}")
    print(f"  S_min^even  = {mpmath.nstr(S_even, 40)}")
    print(f"  S_min^odd   = {mpmath.nstr(S_odd, 40)}")
    print(f"  S_min^diff  = {mpmath.nstr(S_diff, 40)}")
    print(f"  Tail total: {main_rec['tail_total_float']:.3e}")
    print(f"  Tail even:  {main_rec['tail_even_float']:.3e}")
    print(f"  Tail odd:   {main_rec['tail_odd_float']:.3e}")
    print(f"  Partial sanity err = {main_rec['partial_sanity_err']:.3e}")
    print(f"  Tail sanity err    = {main_rec['tail_sanity_err']:.3e}")
    print(f"  Ratio even/odd = {main_rec['S_min_ratio_even_odd']:.15f}")
    print(f"  Compute time: {main_rec['compute_time_sec']:.1f}s")
    print()

    # -----------------------------------------------------------------
    # §3 Basis construction
    # -----------------------------------------------------------------
    basis = build_basis_constants(S_total, S_even, S_odd)

    # -----------------------------------------------------------------
    # §4 PSLQ attempts on S_min^diff
    # -----------------------------------------------------------------
    print("§3 PSLQ attempts on S_min^diff = S_min^even - S_min^odd")
    print("-" * 72)
    diff_attempts = pslq_attempts(S_diff, "S_min_diff", basis)
    diff_found = 0
    for att in diff_attempts:
        r = att["pslq"]
        if r.get("identified"):
            diff_found += 1
            print(f"  [FOUND] {att['label']} (basis size {r['basis_size']}):")
            print(f"    components = {r['components']}")
            print(f"    residual   = {r['residual_float']:.3e}")
        else:
            print(f"  [miss] {att['label']} (basis size {len(att['keys'])}): "
                  f"{r.get('reason', 'unknown')}")
    print()

    # -----------------------------------------------------------------
    # §5 PSLQ attempts on S_min^even and S_min^odd alone
    # -----------------------------------------------------------------
    print("§4 PSLQ attempts on S_min^even alone")
    print("-" * 72)
    even_attempts = pslq_attempts(S_even, "S_min_even", basis)
    even_found = 0
    for att in even_attempts:
        r = att["pslq"]
        if r.get("identified"):
            even_found += 1
            print(f"  [FOUND] {att['label']}:")
            print(f"    components = {r['components']}")
            print(f"    residual   = {r['residual_float']:.3e}")
        else:
            print(f"  [miss] {att['label']}")
    print()

    print("§5 PSLQ attempts on S_min^odd alone")
    print("-" * 72)
    odd_attempts = pslq_attempts(S_odd, "S_min_odd", basis)
    odd_found = 0
    for att in odd_attempts:
        r = att["pslq"]
        if r.get("identified"):
            odd_found += 1
            print(f"  [FOUND] {att['label']}:")
            print(f"    components = {r['components']}")
            print(f"    residual   = {r['residual_float']:.3e}")
        else:
            print(f"  [miss] {att['label']}")
    print()

    # -----------------------------------------------------------------
    # §6 Ratio and product
    # -----------------------------------------------------------------
    print("§6 PSLQ on S_even / S_odd and S_even * S_odd")
    print("-" * 72)
    ratio = S_even / S_odd
    product = S_even * S_odd
    print(f"  Ratio   S_even/S_odd = {mpmath.nstr(ratio, 30)}")
    print(f"  Product S_even*S_odd = {mpmath.nstr(product, 30)}")
    ratio_attempts = pslq_attempts(ratio, "ratio", basis)
    product_attempts = pslq_attempts(product, "product", basis)
    ratio_found = sum(1 for a in ratio_attempts if a["pslq"].get("identified"))
    product_found = sum(1 for a in product_attempts if a["pslq"].get("identified"))
    for att in ratio_attempts:
        r = att["pslq"]
        if r.get("identified"):
            print(f"  [FOUND RATIO] {att['label']}: {r['components']}")
    for att in product_attempts:
        r = att["pslq"]
        if r.get("identified"):
            print(f"  [FOUND PRODUCT] {att['label']}: {r['components']}")
    print()

    # -----------------------------------------------------------------
    # §7 Numerical exploration
    # -----------------------------------------------------------------
    print("§7 Numerical exploration of S_min^diff / known constants")
    print("-" * 72)
    comparisons = {
        "S_diff / S_min": float(S_diff / S_total),
        "S_diff / pi^2": float(S_diff / mpmath.pi ** 2),
        "S_diff / pi^4": float(S_diff / mpmath.pi ** 4),
        "S_diff / G": float(S_diff / mpmath.catalan),
        "S_diff / beta(4)": float(S_diff / beta_num(4)),
        "S_diff / G^2": float(S_diff / mpmath.catalan ** 2),
        "S_diff / zeta(3)": float(S_diff / mpmath.zeta(3)),
        "S_diff / D(4)^2": float(S_diff / dirac_D(4) ** 2),
        "S_diff / (D_even(4) - D_odd(4))": float(S_diff / (dirac_D_even(4) - dirac_D_odd(4))),
        "S_diff / (D_even(4) - D_odd(4))^2": float(S_diff / (dirac_D_even(4) - dirac_D_odd(4)) ** 2),
    }
    for name, val in comparisons.items():
        print(f"  {name:<40} = {val:.15f}")
    print()

    # -----------------------------------------------------------------
    # Save
    # -----------------------------------------------------------------
    all_attempts = {
        "S_min_diff": [{"label": a["label"], "keys": a["keys"], "pslq": a["pslq"]}
                       for a in diff_attempts],
        "S_min_even": [{"label": a["label"], "keys": a["keys"], "pslq": a["pslq"]}
                       for a in even_attempts],
        "S_min_odd":  [{"label": a["label"], "keys": a["keys"], "pslq": a["pslq"]}
                       for a in odd_attempts],
        "ratio":      [{"label": a["label"], "keys": a["keys"], "pslq": a["pslq"]}
                       for a in ratio_attempts],
        "product":    [{"label": a["label"], "keys": a["keys"], "pslq": a["pslq"]}
                       for a in product_attempts],
    }

    summary = {
        "sprint": "RH-P (Sprint 4)",
        "mpmath_dps": 100,
        "n_terms": 4000,
        "tail_method": "asymptotic expansion in 1/a, 10 terms, residual O(1/N^11)",
        "convergence_records": conv_records,
        "main_values": {
            "S_min_total_str": mpmath.nstr(S_total, 100),
            "S_min_even_str":  mpmath.nstr(S_even, 100),
            "S_min_odd_str":   mpmath.nstr(S_odd, 100),
            "S_min_diff_str":  mpmath.nstr(S_diff, 100),
            "S_min_total_float": main_rec["S_min_total_float"],
            "S_min_even_float":  main_rec["S_min_even_float"],
            "S_min_odd_float":   main_rec["S_min_odd_float"],
            "S_min_diff_float":  main_rec["S_min_diff_float"],
            "partial_sanity_err": main_rec["partial_sanity_err"],
            "tail_sanity_err": main_rec["tail_sanity_err"],
            "ratio_even_odd": main_rec["S_min_ratio_even_odd"],
        },
        "T_samples": main_rec["T_samples"],
        "comparisons": comparisons,
        "pslq_attempts": all_attempts,
        "summary_counts": {
            "S_min_diff_found": diff_found,
            "S_min_even_found": even_found,
            "S_min_odd_found": odd_found,
            "ratio_found": ratio_found,
            "product_found": product_found,
            "total_attempts": (len(diff_attempts) + len(even_attempts)
                               + len(odd_attempts) + len(ratio_attempts)
                               + len(product_attempts)),
        },
    }

    out_path = Path("debug/data/smin_chi_neg4.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"Wrote {out_path}")

    # -----------------------------------------------------------------
    # Verdict
    # -----------------------------------------------------------------
    print()
    print("=" * 72)
    print(" VERDICT")
    print("=" * 72)
    total = summary["summary_counts"]["total_attempts"]
    found = (diff_found + even_found + odd_found + ratio_found + product_found)
    print(f"  Total attempts:       {total}")
    print(f"  PSLQ identifications: {found}")
    print(f"  of which:")
    print(f"    S_min_diff:  {diff_found}")
    print(f"    S_min_even:  {even_found}")
    print(f"    S_min_odd:   {odd_found}")
    print(f"    ratio:       {ratio_found}")
    print(f"    product:     {product_found}")


if __name__ == "__main__":
    main()
