"""
Backing tests for Paper 60 eq:blowup -- the atomic L2-Loewdin Jordan-Wigner LCU
1-norm inflation of a naive encoding of the shared-scale Coulomb-Sturmian basis
vs a well-conditioned hydrogenic basis.

Promotes the diagnostic driver debug/io_ladder_sturmian_lambda.py (untracked,
not regression-protected) to a tracked module (geovac/sturmian_l2_encoding.py)
+ this regression test, so the claimed lambda ~ Q^3.33 (sturmian) / Q^1.19
(hydrogenic) exponents are reproduced under CI rather than only in a one-off
debug run.

Marked slow: builds and Jordan-Wigner-encodes several N-shell fermionic
Hamiltonians via openfermion.
"""
import numpy as np
import pytest

from geovac.sturmian_l2_encoding import (
    fit_lambda_exponent,
    lambda_sweep,
    validate_n1,
)

Z = 2.0  # He, matching the diagnostic driver


@pytest.mark.slow
def test_paper60_l2_encoding_n1_validation():
    """N=1 sanity: F0(1s,1s) = 5Z/8 and h1(1s) = -Z^2/2 (exact single-shell values)."""
    F0, h11 = validate_n1(Z=Z)
    F0_exact, h11_exact = 5 * Z / 8, -Z * Z / 2
    print(f"\nF0(1s,1s)={F0} (exact {F0_exact});  h1(1s)={h11} (exact {h11_exact})")
    # grid-resolution-limited (finite trapezoid quadrature, radial cutoff at r_max);
    # sub-0.3% relative agreement confirms the integral machinery, not exactness.
    assert abs(F0 - F0_exact) / F0_exact < 3e-3, f"F0(1s,1s)={F0} vs exact {F0_exact}"
    assert abs(h11 - h11_exact) / abs(h11_exact) < 3e-3, f"h1(1s)={h11} vs exact {h11_exact}"


@pytest.mark.slow
def test_paper60_l2_encoding_blowup_exponents():
    """eq:blowup: the faithful Jordan-Wigner LCU 1-norm lambda inflates as
    Q^1.19 (hydrogenic) vs Q^3.33 (shared-scale Sturmian) under naive Loewdin
    orthonormalization. Sweeps N=1..5 (Q=2N=2..10); the exponent fit excludes
    the degenerate N=1 point, matching the diagnostic driver's convention."""
    N_values = [1, 2, 3, 4, 5]

    hydro_sweep = lambda_sweep(N_values, "hydrogenic", Z=Z)
    sturm_sweep = lambda_sweep(N_values, "sturmian", Z=Z)

    p_hydro = fit_lambda_exponent(N_values, "hydrogenic", Z=Z)
    p_sturm = fit_lambda_exponent(N_values, "sturmian", Z=Z)

    print(f"\nhydrogenic lambda_excl: {[rec['lam_excl'] for rec in hydro_sweep]}")
    print(f"sturmian   lambda_excl: {[rec['lam_excl'] for rec in sturm_sweep]}")
    print(f"fitted exponents: hydrogenic p={p_hydro:.3f}  sturmian p={p_sturm:.3f}")

    # tight bands around the paper's reported exponents (1.19 / 3.33)
    assert 1.0 < p_hydro < 1.4, f"hydrogenic exponent out of band: {p_hydro}"
    assert 3.0 < p_sturm < 3.6, f"sturmian exponent out of band: {p_sturm}"

    # sanity: the shared-scale (sturmian) encoding must inflate faster than the
    # well-conditioned hydrogenic one, and lambda itself must be larger at
    # matched Q for every N > 1 (the mechanism eq:blowup documents).
    assert p_sturm > p_hydro + 1.0
    for h_rec, s_rec in zip(hydro_sweep, sturm_sweep):
        if h_rec["N"] > 1:
            assert s_rec["lam_excl"] > h_rec["lam_excl"], (
                f"sturmian lambda should exceed hydrogenic at N={h_rec['N']}: "
                f"{s_rec['lam_excl']} vs {h_rec['lam_excl']}"
            )

    # loose numerical anchor to the exact target values measured for this sprint
    # (Q=4,6,8,10): hydrogenic ~ 2.25/3.57/5.10/6.68; sturmian ~ 5.67/24.27/60.74/119.72
    hydro_by_Q = {rec["Q"]: rec["lam_excl"] for rec in hydro_sweep}
    sturm_by_Q = {rec["Q"]: rec["lam_excl"] for rec in sturm_sweep}
    hydro_targets = {4: 2.248, 6: 3.571, 8: 5.101, 10: 6.679}
    sturm_targets = {4: 5.669, 6: 24.267, 8: 60.735, 10: 119.724}
    for Q, target in hydro_targets.items():
        got = hydro_by_Q[Q]
        assert abs(got - target) / target < 0.15, f"hydrogenic Q={Q}: {got} vs target {target}"
    for Q, target in sturm_targets.items():
        got = sturm_by_Q[Q]
        assert abs(got - target) / target < 0.15, f"sturmian Q={Q}: {got} vs target {target}"
