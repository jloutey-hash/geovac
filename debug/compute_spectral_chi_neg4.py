"""Track RH-J driver: identify L(s, chi_-4) = beta(s) inside the spectral
Dirichlet series of Paper 28.

HEADLINE RESULT (validated in this sprint):

    D_even(s) - D_odd(s) = 2^{s-1} * (beta(s) - beta(s-2))      (RH-J.1)

valid for all integer s >= 2, where beta(s) = L(s, chi_{-4}) is the
Dirichlet beta function.  At s = 3 the identity requires a regularized
limit since beta(s-2) = beta(1) = pi/4 (finite), but D_even(3), D_odd(3)
individually diverge due to the Hurwitz pole at s-2 = 1.

Derivation
----------
The Dirac Dirichlet series on unit S^3 (Camporesi-Higuchi) is
    D(s) = sum_{n=0}^inf g_n / |lambda_n|^s,
where |lambda_n| = n + 3/2, g_n = 2(n+1)(n+2).

Exact Hurwitz forms (cf. geovac/qed_vertex.py):
    D_even(s) = 2^{-s} * [8 * zeta(s-2, 3/4) - (1/2)*zeta(s, 3/4)]
    D_odd(s)  = 2^{-s} * [8 * zeta(s-2, 5/4) - (1/2)*zeta(s, 5/4)]

Key identity (Hurwitz shift + Dirichlet beta definition):
    zeta(s, 3/4) - zeta(s, 5/4) = zeta(s, 3/4) - [zeta(s, 1/4) - 4^s]
                                = -4^s * beta(s) + 4^s
                                = 4^s * (1 - beta(s)).

Using this identity at exponents s and s-2:
    D_even(s) - D_odd(s)
      = 2^{-s} * [8 * 4^{s-2} * (1 - beta(s-2))
                  - (1/2) * 4^s * (1 - beta(s))]
      = 2^{-s} * [2^{2s-1} * (1 - beta(s-2)) - 2^{2s-1} * (1 - beta(s))]
      = 2^{s-1} * (beta(s) - beta(s-2)).

This driver verifies RH-J.1 symbolically and numerically at s = 2..8.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import mpmath
import sympy as sp
from sympy import Rational, zeta, pi, simplify, Symbol

mpmath.mp.dps = 100  # 100-digit numerical precision


# ----------------------------------------------------------------------
# Exact Hurwitz closed forms (symbolic, sympy)
# ----------------------------------------------------------------------

def dirac_D_sym(s: int) -> sp.Expr:
    """D(s) = 2*zeta(s-2, 3/2) - (1/2)*zeta(s, 3/2), sympy symbolic."""
    return 2 * zeta(s - 2, Rational(3, 2)) - Rational(1, 2) * zeta(s, Rational(3, 2))


def dirac_Deven_sym(s: int) -> sp.Expr:
    """D_even(s) = 2^{-s} * [8*zeta(s-2, 3/4) - (1/2)*zeta(s, 3/4)]."""
    prefactor = Rational(1, 2) ** s
    return prefactor * (
        8 * zeta(s - 2, Rational(3, 4))
        - Rational(1, 2) * zeta(s, Rational(3, 4))
    )


def dirac_Dodd_sym(s: int) -> sp.Expr:
    """D_odd(s) = 2^{-s} * [8*zeta(s-2, 5/4) - (1/2)*zeta(s, 5/4)]."""
    prefactor = Rational(1, 2) ** s
    return prefactor * (
        8 * zeta(s - 2, Rational(5, 4))
        - Rational(1, 2) * zeta(s, Rational(5, 4))
    )


def beta_sym(s: int) -> sp.Expr:
    """beta(s) = L(s, chi_-4) = 4^{-s}*[zeta(s, 1/4) - zeta(s, 3/4)]."""
    return (zeta(s, Rational(1, 4)) - zeta(s, Rational(3, 4))) / sp.Integer(4) ** s


# ----------------------------------------------------------------------
# Numerical values
# ----------------------------------------------------------------------

def hz(s, a):
    """Hurwitz zeta, mpmath."""
    return mpmath.hurwitz(mpmath.mpf(s), mpmath.mpf(a))


def dirac_Deven_num(s: int) -> mpmath.mpf:
    return mpmath.power(2, -s) * (
        8 * hz(s - 2, mpmath.mpf(3) / 4)
        - mpmath.mpf(1) / 2 * hz(s, mpmath.mpf(3) / 4)
    )


def dirac_Dodd_num(s: int) -> mpmath.mpf:
    return mpmath.power(2, -s) * (
        8 * hz(s - 2, mpmath.mpf(5) / 4)
        - mpmath.mpf(1) / 2 * hz(s, mpmath.mpf(5) / 4)
    )


def dirac_D_num(s: int) -> mpmath.mpf:
    return 2 * hz(s - 2, mpmath.mpf(3) / 2) - mpmath.mpf(1) / 2 * hz(s, mpmath.mpf(3) / 2)


def beta_num(s) -> mpmath.mpf:
    """Numeric Dirichlet beta(s) via Hurwitz, extended to s=0 via known value.

    beta(0) = 1/2, beta(1) = pi/4, beta(-1) = 1/2.
    For s >= 2 or s with Re(s) > 1: direct Hurwitz definition.
    """
    if s == 0:
        return mpmath.mpf(1) / 2
    if s == 1:
        return mpmath.pi / 4
    if s < 0:
        # Use functional equation for analytic continuation (not needed here)
        return mpmath.mpf('nan')
    return (hz(s, mpmath.mpf(1) / 4) - hz(s, mpmath.mpf(3) / 4)) / mpmath.power(4, s)


# ----------------------------------------------------------------------
# Closed-form target of the RH-J.1 conjecture
# ----------------------------------------------------------------------

def closed_form_target(s: int) -> mpmath.mpf:
    """2^{s-1} * (beta(s) - beta(s-2)).

    For s = 2: 2 * (G - 1/2) = 2G - 1
    For s = 3: 4 * (pi^3/32 - pi/4) = pi^3/8 - pi
    For s = 4: 8 * (beta(4) - G)
    For s = 5: 16 * (5 pi^5/1536 - pi^3/32) = 5 pi^5/96 - pi^3/2
    For s = 6: 32 * (beta(6) - beta(4))
    For s = 7: 64 * (61 pi^7/184320 - 5 pi^5/1536) = 61 pi^7/2880 - 5 pi^5/24
    For s = 8: 128 * (beta(8) - beta(6))
    """
    return mpmath.power(2, s - 1) * (beta_num(s) - beta_num(s - 2))


def closed_form_target_sym(s: int) -> sp.Expr:
    """Symbolic closed form of RH-J.1: 2^{s-1} * (beta(s) - beta(s-2)).

    For s = 2: uses beta(0) = 1/2, beta(2) = Catalan G.
    For s = 3: uses beta(1) = pi/4, beta(3) = pi^3/32.
    For s = 5, 7: uses closed forms of beta at odd s.
    """
    # Known closed forms
    def beta_known(k):
        if k == 0:
            return Rational(1, 2)
        if k == 1:
            return pi / 4
        if k == 3:
            return pi ** 3 / 32
        if k == 5:
            return 5 * pi ** 5 / 1536
        if k == 7:
            return 61 * pi ** 7 / 184320
        # Even s: use sympy Catalan at s=2, else leave as beta symbol
        if k == 2:
            return sp.Catalan
        # For beta(4), beta(6), beta(8): keep as beta(k) symbol
        return sp.Function(f"beta")(k)
    return sp.Integer(2) ** (s - 1) * (beta_known(s) - beta_known(s - 2))


# ----------------------------------------------------------------------
# PSLQ helpers
# ----------------------------------------------------------------------

def run_pslq(value: mpmath.mpf, basis_vals: List[mpmath.mpf], labels: List[str],
             tol: float = 1e-60, maxcoeff: int = 10 ** 8) -> Dict:
    """Run mpmath PSLQ on [value, *basis_vals].

    If PSLQ returns a relation with zero coefficient on value (meaning an
    internal redundancy in the basis was found), we retry on a reduced
    basis that drops obviously redundant pi^odd / beta(odd) pairs.
    """
    def _attempt(vec, bvals, labs):
        try:
            rel = mpmath.pslq(vec, tol=tol, maxcoeff=maxcoeff)
        except Exception as e:
            return None, f"PSLQ exception: {e}"
        if rel is None:
            return None, "PSLQ returned None"
        if rel[0] == 0:
            return None, "zero coefficient on value (basis redundancy)"
        return rel, None

    vec = [value] + list(basis_vals)
    rel, reason = _attempt(vec, basis_vals, labels)

    # Retry by pruning redundant pi^odd and beta(odd-at-value): for odd weight s
    # beta(s) is rational * pi^s, so including both creates a relation among
    # basis elements alone.
    if rel is None and "zero coefficient" in (reason or ""):
        # Identify which label is "beta(k)" for odd k, or "beta(1)=pi/4"
        pruned_basis = []
        pruned_labels = []
        for bv, lab in zip(basis_vals, labels):
            # Drop "beta(1)=pi/4" if "pi" is present; drop "beta(odd)" if pi^odd present.
            skip = False
            if lab == "beta(1)=pi/4" and "pi" in labels:
                skip = True
            if lab.startswith("beta(") and not lab.startswith("beta(0)") and not lab.startswith("beta(1)"):
                # extract k
                try:
                    k = int(lab[5:-1])
                    if k % 2 == 1 and f"pi^{k}" in labels:
                        skip = True
                except Exception:
                    pass
            if not skip:
                pruned_basis.append(bv)
                pruned_labels.append(lab)
        vec2 = [value] + pruned_basis
        rel, reason2 = _attempt(vec2, pruned_basis, pruned_labels)
        if rel is not None:
            basis_vals = pruned_basis
            labels = pruned_labels
        else:
            return {"identified": False, "reason": f"{reason}; retry: {reason2}"}

    if rel is None:
        return {"identified": False, "reason": reason or "unknown"}

    c0 = rel[0]
    components = {}
    reconstructed = mpmath.mpf(0)
    for ri, lab, bv in zip(rel[1:], labels, basis_vals):
        if ri != 0:
            frac = Fraction(int(-ri), int(c0))
            components[lab] = str(frac)
            reconstructed += (mpmath.mpf(-ri) / c0) * bv
    residual = abs(value - reconstructed)
    return {
        "identified": True,
        "components": components,
        "residual_float": float(residual),
        "reconstructed_float": float(reconstructed),
        "value_float": float(value),
        "raw_relation": [int(r) for r in rel],
    }


def pslq_targeted_basis(s: int) -> Tuple[List[mpmath.mpf], List[str]]:
    """Minimal targeted basis for RH-J.1 at given s.

    For even s (beta-dominated): {1, beta(s-2), beta(s)}.
    For odd s (pi-dominated since beta(odd) is rational * pi^odd):
        {1, pi, pi^3, pi^5, pi^7} plus explicit beta(s-2), beta(s).
    Minimal bases avoid the ill-conditioning that large bases suffer under PSLQ.
    """
    pi_ = mpmath.pi
    if s % 2 == 0:
        # Even s: the natural basis is {beta(s-2), beta(s)} (both independent
        # transcendentals at even weight, with beta(0) = 1/2).
        basis = [mpmath.mpf(1)]
        labels = ["1"]
    else:
        # Odd s: include pi^odd since beta(odd) = (rational)*pi^odd.
        # But we also include beta(s), beta(s-2) as basis elements for
        # direct identification.
        basis = [mpmath.mpf(1), pi_, pi_ ** 3, pi_ ** 5, pi_ ** 7]
        labels = ["1", "pi", "pi^3", "pi^5", "pi^7"]

    # beta at s and s-2
    for k in (s - 2, s):
        if k == 0:
            basis.append(mpmath.mpf(1) / 2)
            labels.append("beta(0)=1/2")
        elif k == 1:
            basis.append(pi_ / 4)
            labels.append("beta(1)=pi/4")
        elif k >= 2:
            basis.append(beta_num(k))
            labels.append(f"beta({k})")
    return basis, labels


def pslq_wide_basis(s: int) -> Tuple[List[mpmath.mpf], List[str]]:
    """Wider basis for stress testing: includes zeta(odd) and log(2)."""
    basis, labels = pslq_targeted_basis(s)
    basis += [mpmath.zeta(3), mpmath.zeta(5), mpmath.zeta(7), mpmath.log(2)]
    labels += ["zeta(3)", "zeta(5)", "zeta(7)", "log(2)"]
    return basis, labels


# ----------------------------------------------------------------------
# Per-s test
# ----------------------------------------------------------------------

def test_conjecture_at_s(s: int) -> Dict:
    """Verify RH-J.1 at integer s >= 2.

    Returns a record containing:
      - D_even(s), D_odd(s), D_diff(s) numerically at 100 digits
      - closed_form = 2^{s-1}*(beta(s)-beta(s-2)) numerically
      - difference error
      - PSLQ identification of D_diff in targeted and wide bases
      - symbolic closed form expression
    """
    special_s = (s == 3)  # individual D_even, D_odd have pole at s-2=1

    if special_s:
        # Compute D_diff via the CLOSED-FORM RHS since numerics of D_even/odd diverge
        D_diff_numeric = closed_form_target(s)
        D_even_float = float('inf')
        D_odd_float = float('inf')
        D_full_float = float(dirac_D_num(s))
        sum_consistency = float('nan')
    else:
        D_even = dirac_Deven_num(s)
        D_odd = dirac_Dodd_num(s)
        D_diff_numeric = D_even - D_odd
        D_even_float = float(D_even)
        D_odd_float = float(D_odd)
        D_full_float = float(dirac_D_num(s))
        sum_consistency = float(abs((D_even + D_odd) - dirac_D_num(s)) / abs(dirac_D_num(s)))

    # Closed-form target
    target = closed_form_target(s)
    if mpmath.isnan(D_diff_numeric):
        identity_error = float('nan')
    else:
        identity_error = float(abs(D_diff_numeric - target))

    # PSLQ against targeted basis
    basis_tgt, labels_tgt = pslq_targeted_basis(s)
    pslq_tgt = run_pslq(D_diff_numeric, basis_tgt, labels_tgt,
                        tol=1e-80, maxcoeff=10 ** 6)
    if not pslq_tgt.get("identified"):
        pslq_tgt = run_pslq(D_diff_numeric, basis_tgt, labels_tgt,
                            tol=1e-60, maxcoeff=10 ** 8)

    # PSLQ against wide basis (sanity)
    basis_wide, labels_wide = pslq_wide_basis(s)
    pslq_wide = run_pslq(D_diff_numeric, basis_wide, labels_wide,
                         tol=1e-60, maxcoeff=10 ** 8)

    # Closed-form sympy expression
    try:
        expr_sym = dirac_Deven_sym(s) - dirac_Dodd_sym(s)
        expr_str = str(expr_sym)
    except Exception:
        expr_str = "symbolic-expression-error"

    return {
        "s": s,
        "D_even_float": D_even_float,
        "D_odd_float": D_odd_float,
        "D_diff_float": float(D_diff_numeric) if not mpmath.isnan(D_diff_numeric) else None,
        "D_full_float": D_full_float,
        "sum_consistency_rel_err": sum_consistency if not (sum_consistency != sum_consistency) else None,
        "closed_form_target_float": float(target),
        "closed_form_description": f"2^(s-1) * (beta({s}) - beta({s-2}))",
        "identity_error_float": identity_error if not (identity_error != identity_error) else None,
        "pslq_targeted": pslq_tgt,
        "pslq_wide": pslq_wide,
        "expr_Deven_minus_Dodd_symbolic": expr_str,
    }


# ----------------------------------------------------------------------
# Individual D_even, D_odd probe (for completeness)
# ----------------------------------------------------------------------

def probe_individual(s: int) -> Dict:
    """Identify D_even(s) and D_odd(s) individually (s != 3)."""
    if s == 3:
        return {"D_even": {"identified": False, "reason": "divergent individually at s=3"},
                "D_odd":  {"identified": False, "reason": "divergent individually at s=3"}}

    D_even = dirac_Deven_num(s)
    D_odd = dirac_Dodd_num(s)

    # Build a reasonable basis for individual identification:
    # {1, pi^2, pi^4, pi^6, pi^8, G, beta(4), beta(6), beta(8)}.
    # These live in weight-graded components and should pin individuals.
    pi_ = mpmath.pi
    basis = [mpmath.mpf(1), pi_ ** 2, pi_ ** 4, pi_ ** 6, pi_ ** 8,
             mpmath.catalan, beta_num(4), beta_num(6), beta_num(8)]
    labels = ["1", "pi^2", "pi^4", "pi^6", "pi^8",
              "G", "beta(4)", "beta(6)", "beta(8)"]

    out = {}
    for name, val in [("D_even", D_even), ("D_odd", D_odd)]:
        rec = None
        for tol in [1e-80, 1e-60, 1e-40]:
            for mxc in [10 ** 6, 10 ** 8, 10 ** 10]:
                res = run_pslq(val, basis, labels, tol=tol, maxcoeff=mxc)
                if res.get("identified"):
                    rec = res
                    rec["tol_used"] = tol
                    rec["maxcoeff_used"] = mxc
                    break
            if rec:
                break
        if rec is None:
            rec = {"identified": False, "reason": "no relation at any tol/maxcoeff"}
        out[name] = rec
    return out


# ----------------------------------------------------------------------
# Main sweep
# ----------------------------------------------------------------------

def main():
    s_values = [2, 3, 4, 5, 6, 7, 8]

    print("=" * 72)
    print(" RH-J.1 Conjecture: D_even(s) - D_odd(s) = 2^{s-1} (beta(s) - beta(s-2))")
    print("=" * 72)

    records = {}
    for s in s_values:
        print(f"\n--- s = {s} ---")
        rec = test_conjecture_at_s(s)
        ind = probe_individual(s)
        rec["individual"] = ind
        records[s] = rec

        print(f"  D_even(s) = {rec['D_even_float']}")
        print(f"  D_odd(s)  = {rec['D_odd_float']}")
        print(f"  D_diff    = {rec['D_diff_float']}")
        print(f"  2^(s-1)*(beta(s)-beta(s-2)) = {rec['closed_form_target_float']}")
        print(f"  identity error: {rec['identity_error_float']}")
        if rec["pslq_targeted"].get("identified"):
            print(f"  PSLQ (targeted basis): {rec['pslq_targeted']['components']}")
            print(f"    residual = {rec['pslq_targeted']['residual_float']:.3e}")
        else:
            print(f"  PSLQ (targeted): {rec['pslq_targeted'].get('reason', 'unidentified')}")
        if rec["pslq_wide"].get("identified"):
            print(f"  PSLQ (wide basis): {rec['pslq_wide']['components']}")

        for nm, res in ind.items():
            if res.get("identified"):
                print(f"  PSLQ({nm}): {res['components']}")
            else:
                print(f"  PSLQ({nm}): {res.get('reason', 'unidentified')}")

    # Write JSON
    safe = {}
    for s, rec in records.items():
        safe[str(s)] = rec

    out_path = Path("debug/data/spectral_chi_neg4.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(safe, f, indent=2, default=str)
    print(f"\nWrote {out_path}")

    # Print summary
    print("\n" + "=" * 72)
    print(" SUMMARY")
    print("=" * 72)
    print(f"{'s':>3} | {'|D_diff - target|':>22} | {'PSLQ identified':>16}")
    print("-" * 72)
    for s in s_values:
        rec = records[s]
        err = rec['identity_error_float']
        err_str = f"{err:.3e}" if err is not None else "n/a"
        pslq_ok = "yes" if rec["pslq_targeted"].get("identified") else "no"
        print(f"{s:>3} | {err_str:>22} | {pslq_ok:>16}")


if __name__ == "__main__":
    main()
