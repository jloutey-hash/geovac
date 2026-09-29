"""Paper 16 §VI — Dirac instability is a METRIC, not topological, singularity.

The §VI headline (one of P16's five claims) had NO test across multiple QA certs.
It is an analytic claim, so the right backing is a numerical cross-check that the
paper's own §VI formulas reproduce its tabulated values (tab:metric) and that the
quantum-number side (nu, mu_free) is genuinely SMOOTH through Z = 1/alpha ~ 137
(no singularity / phase transition) while the metric side (delta_1s) diverges.

This is a formula/identity verification (CLAUDE.md §13.4a "symbolic identity" +
"numerical cross-check"), self-contained from the paper's equations.
"""
import math

ALPHA = 1.0 / 137.036  # fine-structure constant (paper §VI uses Z*alpha)


def _nu_min(N: int) -> int:
    """Paper 16 eq:nu_min (corrected 2026-09-28): lowest-spin Pauli filling,
    E(ceil(N/2)) + E(floor(N/2)), degree d holding (d+1)(d+2)/2 functions."""
    def e(m: int) -> int:
        tot, d = 0, 0
        while m > 0:
            k = min(m, (d + 1) * (d + 2) // 2)
            tot, m, d = tot + k * d, m - k, d + 1
        return tot
    return e((N + 1) // 2) + e(N // 2)


def _mu_free(N: int) -> int:
    """Paper 16 eq:mu_formula: SO(3N) Casimir nu(nu+3N-2)/2 at nu = nu_min."""
    nu = _nu_min(N)
    return nu * (nu + 3 * N - 2) // 2


def _delta_1s(Z: int) -> float:
    """Paper 16 eq:metric_correction: conformal metric warp of the 1s shell."""
    za2 = (Z * ALPHA) ** 2
    return za2 / (1.0 - za2)


def test_mu_free_smooth_through_Z137() -> None:
    """Quantum-number side is smooth at N=137 (no singularity): the integer formula
    gives nu=570, mu_free=279,015 (paper §VI, corrected 2026-09-28), finite and exact.
    Rejects the retracted nu = N-2 = 135 / mu_free = 36,720."""
    N = 137
    assert _nu_min(N) == 570                 # E(69) + E(68) = 288 + 282
    assert _mu_free(N) == 279015             # 570 * 979 / 2


def test_mu_free_monotone_finite_across_Z137() -> None:
    """No singularity / phase transition in the topology side through Z=1/alpha:
    mu_free is finite and strictly increasing across N = 130..140."""
    vals = [_mu_free(N) for N in range(130, 141)]
    assert all(math.isfinite(v) for v in vals)
    assert all(b > a for a, b in zip(vals, vals[1:]))  # strictly monotone, no pinch


def test_delta_1s_table_values() -> None:
    """Metric side reproduces tab:metric (paper §VI): perturbative -> singular.
    Every delta_1s row of the table is pinned (QA 2026-07-01 tightening)."""
    assert abs(_delta_1s(1) - 0.00005) < 0.00002  # non-relativistic
    assert abs(_delta_1s(29) - 0.047) < 0.002    # mildly relativistic
    assert abs(_delta_1s(47) - 0.133) < 0.003    # relativistic
    assert abs(_delta_1s(79) - 0.498) < 0.005    # strongly relativistic
    assert abs(_delta_1s(118) - 2.87) < 0.05     # near-singular
    assert _delta_1s(137) > 1.0e3                # conical singularity (~10^3), the metric pinch


def test_metric_diverges_while_topology_does_not() -> None:
    """The §VI thesis: as Z*alpha -> 1 the metric (delta_1s) diverges while the
    quantum numbers (mu_free) stay finite/integer -- metric singularity, not topological."""
    assert _delta_1s(137) > _delta_1s(79) > _delta_1s(29)   # metric blows up
    assert _delta_1s(137) > 1.0e3                           # ...and by MAGNITUDE (local self-sufficiency)
    assert _mu_free(137) == 279015 and isinstance(_mu_free(137), int)  # topology untouched
