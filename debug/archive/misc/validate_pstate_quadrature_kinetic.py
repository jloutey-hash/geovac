"""Validate the new universal P-state quadrature kinetic against the
algebraic symxsym implementation.

If the quadrature implementation is correct, symxsym kinetic via
quadrature should match the algebraic value at quadrature precision
(~1e-4 relative for n_r=32, n_theta=16).

Additionally validates:
- antisymxantisym at (l=m=n=0, beta=0) returns 0 (basis vanishes).
- cross-sector at beta=0 returns 0 (sinh→0).
- Hermiticity of T_pq at finite beta across all channel pairs.
"""

import math
import numpy as np

from geovac.hylleraas_eckart_pstate import (
    HylleraasPStateBasisFn,
    hylleraas_pstate_basis_total_degree,
    kinetic_element_pstate_eckart_sym_sym,
    _kinetic_via_quadrature_pstate,
)


def test_sym_sym_quadrature_vs_algebraic():
    """symxsym quadrature should agree with algebraic to ~1e-4."""
    basis = hylleraas_pstate_basis_total_degree(2, channel='sym')
    alpha = 1.35
    for beta in [0.0, 0.3]:
        worst = 0.0
        for bf_p in basis:
            for bf_q in basis:
                t_alg = kinetic_element_pstate_eckart_sym_sym(
                    bf_p, bf_q, alpha, beta
                )
                t_qua = _kinetic_via_quadrature_pstate(
                    bf_p, bf_q, alpha, beta, n_r=32, n_theta=16
                )
                denom = max(abs(t_alg), abs(t_qua), 1e-30)
                rel = abs(t_alg - t_qua) / denom
                if rel > worst:
                    worst = rel
                    worst_pair = (bf_p.label(), bf_q.label(), t_alg, t_qua)
        print(f"beta={beta}: worst symxsym rel diff = {worst:.2e}")
        if worst > 1e-3:
            print(f"  worst pair: {worst_pair}")
        assert worst < 1e-3, f"symxsym quadrature off: {worst:.2e}"
    print("  PASS: symxsym quadrature ~= algebraic")


def test_antisym_at_beta_zero():
    """antisym x antisym kinetic at beta=0 should be 0 (basis vanishes)."""
    basis = hylleraas_pstate_basis_total_degree(2, channel='antisym')
    alpha = 1.35
    for bf_p in basis:
        for bf_q in basis:
            t = _kinetic_via_quadrature_pstate(bf_p, bf_q, alpha, 0.0, n_r=24, n_theta=12)
            assert abs(t) < 1e-12, f"antisym at beta=0 not zero: {t}"
    print("  PASS: antisym x antisym kinetic = 0 at beta=0")


def test_cross_at_beta_zero():
    """sym x antisym cross kinetic at beta=0 should be 0 (sinh→0)."""
    basis_sym = hylleraas_pstate_basis_total_degree(2, channel='sym')
    basis_anti = hylleraas_pstate_basis_total_degree(2, channel='antisym')
    alpha = 1.35
    for bf_p in basis_sym:
        for bf_q in basis_anti:
            t = _kinetic_via_quadrature_pstate(bf_p, bf_q, alpha, 0.0, n_r=24, n_theta=12)
            assert abs(t) < 1e-12, f"cross at beta=0 not zero: {t}"
    print("  PASS: sym x antisym cross kinetic = 0 at beta=0")


def test_kinetic_hermitian():
    """T_pq = T_qp at finite beta across all channel combinations."""
    basis_sym = hylleraas_pstate_basis_total_degree(2, channel='sym')
    basis_anti = hylleraas_pstate_basis_total_degree(2, channel='antisym')
    alpha = 1.35
    beta = 0.3
    n_r, n_theta = 24, 12
    failures = []
    # Within symxsym
    for bf_p in basis_sym:
        for bf_q in basis_sym:
            t_pq = _kinetic_via_quadrature_pstate(bf_p, bf_q, alpha, beta, n_r, n_theta)
            t_qp = _kinetic_via_quadrature_pstate(bf_q, bf_p, alpha, beta, n_r, n_theta)
            if abs(t_pq - t_qp) > 1e-8:
                failures.append(('symxsym', bf_p.label(), bf_q.label(), t_pq, t_qp))
    # Within antisymxantisym
    for bf_p in basis_anti:
        for bf_q in basis_anti:
            t_pq = _kinetic_via_quadrature_pstate(bf_p, bf_q, alpha, beta, n_r, n_theta)
            t_qp = _kinetic_via_quadrature_pstate(bf_q, bf_p, alpha, beta, n_r, n_theta)
            if abs(t_pq - t_qp) > 1e-8:
                failures.append(('antixanti', bf_p.label(), bf_q.label(), t_pq, t_qp))
    # Cross: T(sym, antisym) vs T(antisym, sym) — these should also match by Hermiticity
    for bf_p in basis_sym[:3]:
        for bf_q in basis_anti[:3]:
            t_pq = _kinetic_via_quadrature_pstate(bf_p, bf_q, alpha, beta, n_r, n_theta)
            t_qp = _kinetic_via_quadrature_pstate(bf_q, bf_p, alpha, beta, n_r, n_theta)
            if abs(t_pq - t_qp) > 1e-8:
                failures.append(('cross', bf_p.label(), bf_q.label(), t_pq, t_qp))
    if failures:
        for f in failures[:5]:
            print(f"  HERMITICITY FAIL: {f}")
        assert False, f"{len(failures)} Hermiticity failures"
    print("  PASS: kinetic Hermitian across all channel pairs")


def test_antisym_at_lmn_zero_beta_finite():
    """antisymxantisym kinetic at (l=m=n=0, beta>0) is positive (kinetic energy)."""
    bf = HylleraasPStateBasisFn(0, 0, 0, channel='antisym')
    alpha = 1.5
    beta = 0.3
    t_pp = _kinetic_via_quadrature_pstate(bf, bf, alpha, beta, n_r=40, n_theta=20)
    print(f"  T(antisym,000)@beta=0.3: {t_pp:.6f}")
    # Should be positive and on order of magnitude of kinetic energy.
    assert t_pp > 0
    print("  PASS: antisym (000) kinetic positive at beta>0")


def main():
    print("=== Validating universal P-state quadrature kinetic ===")
    print()
    print("1) sym x sym: quadrature vs algebraic")
    test_sym_sym_quadrature_vs_algebraic()
    print()
    print("2) antisym x antisym at beta=0 (basis vanishes)")
    test_antisym_at_beta_zero()
    print()
    print("3) sym x antisym cross at beta=0 (sinh to 0)")
    test_cross_at_beta_zero()
    print()
    print("4) Hermiticity at finite beta across all channels")
    test_kinetic_hermitian()
    print()
    print("5) antisym (000) kinetic at beta>0 is positive")
    test_antisym_at_lmn_zero_beta_finite()
    print()
    print("=== ALL VALIDATION CHECKS PASSED ===")


if __name__ == '__main__':
    main()
