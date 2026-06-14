"""
TRUNK QA — Claim 5: gearing ratio ||L+||/||T+|| -> 1.77 and [T+,L+] = 0 (Paper 1).

Paper 1 sec:conclusions states:
  - "The gearing ratio ||L+||/||T+|| converges to ~1.77, not 137 or 1/alpha."
  - "[T+, L+] = 0, indicating perfect integrability."

We build T+ and L+ from the EXACT Biedenharn-Louck matrix elements Paper 1
gives in its appendix (Eq. with the SU(2)/SU(1,1) coefficients):
  <l,m+1|L+|l,m>   = sqrt((l-m)(l+m+1))
  <n+1,l|T+|n,l>   = sqrt((n-l)(n+l+1)/4)
and compute the operator norms and the commutator directly. Nothing is
assumed about the answer.

FINDINGS (reported honestly):
  - [T+, L+] = 0 BIT-EXACT at every n_max -> COMMUTATOR CLAIM BACKED.
  - SPECTRAL-norm gearing ||L+||_2/||T+||_2 = 2.000 EXACTLY (not 1.77).
  - FROBENIUS-norm gearing ~1.74-1.84, n_max-dependent, NOT converging to a
    single 1.77 (1.84 at n_max=10, 1.74 at n_max=20). So the "1.77" is at best
    a Frobenius-norm value at a particular n_max, not a convergent ratio.
  => Commutator: BACKED. Gearing "1.77 convergence": NOT cleanly backed
     (spectral norm gives exactly 2; Frobenius drifts). DOWNGRADE the "1.77".
"""

from __future__ import annotations

import math

import numpy as np


def build_T_L(max_n: int):
    """T+ and L+ as dense matrices from Paper 1's documented matrix elements."""
    states = []
    for n in range(1, max_n + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
    idx = {s: i for i, s in enumerate(states)}
    N = len(states)
    Lp = np.zeros((N, N))
    Tp = np.zeros((N, N))
    for (n, l, m) in states:
        i = idx[(n, l, m)]
        # L+ : <l,m+1|L+|l,m> = sqrt((l-m)(l+m+1))
        if m < l:
            j = idx[(n, l, m + 1)]
            Lp[j, i] = math.sqrt((l - m) * (l + m + 1))
        # T+ : <n+1,l|T+|n,l> = sqrt((n-l)(n+l+1)/4)
        if (n + 1, l, m) in idx:
            j = idx[(n + 1, l, m)]
            Tp[j, i] = math.sqrt((n - l) * (n + l + 1) / 4.0)
    return states, idx, Lp, Tp


def test_commutator_is_exactly_zero():
    """[T+, L+] = 0 bit-exact (the integrability claim). BACKED."""
    for max_n in (5, 10, 15, 20):
        _, _, Lp, Tp = build_T_L(max_n)
        comm = Tp @ Lp - Lp @ Tp
        assert np.linalg.norm(comm, 2) < 1e-10, (
            f"n_max={max_n}: [T+,L+] norm {np.linalg.norm(comm,2):.2e} != 0"
        )


def test_spectral_gearing_is_exactly_two_not_177():
    """Spectral-norm gearing ||L+||_2/||T+||_2 = 2.000 exactly — NOT 1.77.
    This is the operator-norm reading; it contradicts the paper's 1.77."""
    for max_n in (5, 10, 15, 20):
        _, _, Lp, Tp = build_T_L(max_n)
        ratio = np.linalg.norm(Lp, 2) / np.linalg.norm(Tp, 2)
        assert abs(ratio - 2.0) < 1e-6, f"n_max={max_n}: spectral gearing {ratio}"


def test_frobenius_gearing_does_not_converge_to_177():
    """Frobenius-norm gearing drifts MONOTONICALLY DOWN with n_max
    (1.84 @10 -> 1.77 @15 -> 1.74 @20 -> 1.70 @30): it passes THROUGH 1.77 near
    n_max=15 but keeps decreasing toward ~1.70, so it does NOT converge to
    1.77. Document the drift (DOWNGRADE signal for the '~1.77 convergence').
    """
    ratios = {}
    for max_n in (10, 20, 30):
        _, _, Lp, Tp = build_T_L(max_n)
        ratios[max_n] = np.linalg.norm(Lp, "fro") / np.linalg.norm(Tp, "fro")
    # Monotone decreasing through the range — not converged to 1.77.
    assert ratios[10] > ratios[20] > ratios[30], (
        f"Frobenius gearing not monotone-decreasing: {ratios}"
    )
    assert ratios[10] > 1.77 > ratios[30], (
        f"1.77 is mid-drift, not a limit: {ratios}"
    )
    # The value is still moving by >0.05 across the range (not converged).
    assert abs(ratios[10] - ratios[30]) > 0.05, (
        f"Frobenius gearing barely moves; check claim: {ratios}"
    )


def test_177_is_not_reproduced_by_any_standard_operator_norm():
    """Neither spectral (=2.000) nor Frobenius (=1.74-1.84, drifting) gives a
    convergent 1.77. The '~1.77' claim is not cleanly backed by the documented
    operators under standard norms."""
    _, _, Lp, Tp = build_T_L(20)
    spec = np.linalg.norm(Lp, 2) / np.linalg.norm(Tp, 2)
    frob = np.linalg.norm(Lp, "fro") / np.linalg.norm(Tp, "fro")
    assert abs(spec - 1.77) > 0.2          # spectral is 2.0, far from 1.77
    # Frobenius happens to land ~1.735 at n_max=20 — close-ish, but it is a
    # moving target (see drift test), so not a *convergent* 1.77.
    assert frob < 2.0                       # at least it's bounded below 2
