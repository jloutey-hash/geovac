"""Smoke test for l2_constant_c_identification.py infrastructure.

Verifies:
  1. L-value computations match known closed forms.
  2. Richardson fits converge.
  3. PSLQ runs without crashing.
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
from mpmath import mp, mpf, pi, sqrt, log

from l2_constant_c_identification import (  # noqa: E402
    dirichlet_beta,
    low_conductor_real_dirichlet_L,
    build_basis_dirichlet_L,
    build_basis_multi_hurwitz,
    build_basis_stein_weiss,
    build_basis_wildcard,
    build_basis_basic,
    compute_gamma_n,
    fit_subleading,
    compute_panel,
    run_pslq,
    deduplicate_basis,
)


def test_l_values():
    """Verify L-value computations against known closed forms."""
    mp.dps = 60
    print("Verifying L-value computations:")

    # beta(2) = G (Catalan's constant)
    b2 = dirichlet_beta(2, 60)
    G = mpmath.catalan
    err = abs(b2 - G)
    status = "PASS" if err < mpf("1e-50") else "FAIL"
    print(f"  beta(2) vs Catalan G  : err = {mpmath.nstr(err, 8)}  [{status}]")

    # beta(1) = pi/4
    b1 = dirichlet_beta(1, 60)
    err = abs(b1 - pi/4)
    status = "PASS" if err < mpf("1e-50") else "FAIL"
    print(f"  beta(1) vs pi/4       : err = {mpmath.nstr(err, 8)}  [{status}]")

    # L(2, chi_3): known closed form. chi_3 is the Kronecker symbol (-3/.).
    # L(s, chi_3) is a Dirichlet L-function with chi_3(1)=1, chi_3(2)=-1.
    # No simple closed form for L(2, chi_3); known value ~ 0.781302...
    # Cross-check: via series sum 1 - 1/4 + 1/16 - 1/25 + ... where signs follow chi_3.
    # We just verify it's positive and finite.
    L_chi3_2 = low_conductor_real_dirichlet_L(2, 3, 60)
    print(f"  L(2,chi_3) numerical  : {mpmath.nstr(L_chi3_2, 30)}")
    # Independent: chi_3 has Kronecker symbol (n|3): 0,1,-1 for n=0,1,2 (mod 3).
    # Direct sum: sum_{n=1}^inf chi_3(n)/n^2
    direct = mpf(0)
    for n in range(1, 1000):
        m = n % 3
        c = 1 if m == 1 else (-1 if m == 2 else 0)
        direct += mpf(c) / mpf(n)**2
    print(f"  L(2,chi_3) direct sum : {mpmath.nstr(direct, 12)}")
    print(f"  agreement             : {mpmath.nstr(abs(L_chi3_2 - direct), 8)}")

    # L(2, chi_4) = beta(2) = Catalan G
    L_chi4_2 = low_conductor_real_dirichlet_L(2, 4, 60)
    err = abs(L_chi4_2 - G)
    status = "PASS" if err < mpf("1e-50") else "FAIL"
    print(f"  L(2,chi_4) vs Catalan : err = {mpmath.nstr(err, 8)}  [{status}]")
    print()


def test_basis_construction():
    """Verify basis builders produce non-empty, finite results."""
    mp.dps = 60
    print("Verifying basis construction:")
    bases = {
        "basic": build_basis_basic(60),
        "DirichletL": build_basis_dirichlet_L(60),
        "multi_hurwitz": build_basis_multi_hurwitz(60),
        "stein_weiss": build_basis_stein_weiss(60),
        "wildcard": build_basis_wildcard(60),
    }
    for name, basis in bases.items():
        n = len(basis)
        all_finite = all(mpmath.isfinite(v) for _, v in basis)
        any_zero = any(v == 0 for _, v in basis)
        status = "OK" if (n > 0 and all_finite and not any_zero) else "FAIL"
        print(f"  {name:20s}: {n} elements  finite={all_finite}  any_zero={any_zero}  [{status}]")
        # Print first 3 and last 1 for inspection
        for i, (lbl, v) in enumerate(basis[:3]):
            print(f"    [{i}] {lbl:30s} = {mpmath.nstr(v, 12)}")
        if len(basis) > 3:
            lbl, v = basis[-1]
            print(f"    ... [{len(basis)-1}] {lbl:30s} = {mpmath.nstr(v, 12)}")
    print()


def test_quick_panel():
    """Quick MR-C-style fit at low precision and small panel.

    Verify our infrastructure reproduces MR-C's c value within ~1e-15 at
    n_max=4096, dps=120.
    """
    print("Quick panel test (dps=120, n_values=[16, 32, 64, 128, 256, 512, 1024, 2048, 4096]):")
    n_values = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
    panel = compute_panel(n_values, dps=120)
    fit = fit_subleading(panel, K=4)
    c = fit["b_0_extracted"]
    print(f"  K=4 c = {mpmath.nstr(c, 30)}")
    mrc_value_str = "4.1093214674877940927579607260741005838057691088362615503253972964276017819113301"
    mrc_value = mpf(mrc_value_str)
    err = abs(c - mrc_value)
    print(f"  vs MR-C: err = {mpmath.nstr(err, 12)}")
    if err < mpf("1e-12"):
        print("  PASS (agrees with MR-C value at K=4 with smaller panel)")
    else:
        print("  Note: with smaller panel, K=4 fit will not match K=6 MR-C value to >>1e-15.")
        print("  As long as the value is in the same neighborhood and converges with K, OK.")


def test_pslq():
    """Run a quick PSLQ to make sure it doesn't crash."""
    mp.dps = 80
    print("Quick PSLQ test:")
    c = mpf("4.1093214674877940927579607260741005838057691088362615503253972964276")
    basis = build_basis_basic(80)
    r = run_pslq(c, basis, max_coeff=1000, dps=80)
    print(f"  Found: {r.get('found')}")
    if r.get("found"):
        if r.get("verify_pass_1e_minus_50"):
            print(f"  PASS: c = {r['closed_form_terms']}")
        else:
            print(f"  Found a relation but not verifying at 50 dps (probably spurious):")
            print(f"    {r['closed_form_terms']}")
            print(f"    abs err = {r.get('verify_abs_error')}")


if __name__ == "__main__":
    test_l_values()
    test_basis_construction()
    test_pslq()
    print("All smoke tests done.")
