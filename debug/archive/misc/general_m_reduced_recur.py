r"""Recurrence evaluator for the reduced associated-Legendre functions on the xi-grid,
the stable replacement for prolate_ci_general_m.py's differentiation route.

The V_ee assembly there uses
    pv = legendre_deriv_poly(l, m)   (= d^m P_l/dxi^m, a polynomial)
    qv = q_deriv(l, m, xi)           (= d^m Q_l/dxi^m)
whose polynomial coefficients reach ~1e10 at large l and overflow -> the mu=2
(d^4) instability. Both REDUCED functions R^P_l = d^m P_l/dxi^m and
R^Q_l = d^m Q_l/dxi^m satisfy the SAME associated-Legendre l-recurrence

    (l-m+1) R_{l+1} = (2l+1) xi R_l - (l+m) R_{l-1},

so we seed at the two lowest orders l=m, m+1 (where the differentiation is cheap
and stable) and RECUR upward -- no high-order differentiation, no 1e10
coefficients. This tests float64 stability (the pipeline's precision) of that
recurrence against the direct differentiation, for m up to 4 and l up to 24.
"""
import sys; sys.path.insert(0, 'debug')
import numpy as np
import numpy.polynomial.polynomial as P
from prolate_ci_general_m import legendre_deriv_poly, q_deriv


def reduced_recur(l_max, m, xi):
    """R^P_l(xi), R^Q_l(xi) for l=m..l_max via forward l-recurrence (float64)."""
    RP = {m: P.polyval(xi, legendre_deriv_poly(m, m)),
          m + 1: P.polyval(xi, legendre_deriv_poly(m + 1, m))}
    RQ = {m: q_deriv(m, m, xi), m + 1: q_deriv(m + 1, m, xi)}
    for l in range(m + 1, l_max):
        c = 1.0 / (l - m + 1)
        RP[l + 1] = c * ((2 * l + 1) * xi * RP[l] - (l + m) * RP[l - 1])
        RQ[l + 1] = c * ((2 * l + 1) * xi * RQ[l] - (l + m) * RQ[l - 1])
    return RP, RQ


def reduced_Q_backward(l_max, m, xi, l_top=None):
    """Reduced Q_l^m via BACKWARD (Miller) recurrence -- stable for the recessive
    (decaying) Q solution. Seed arbitrarily high, recur down, normalise each grid
    point to the direct q_deriv at l=m (low l, cheap and stable)."""
    if l_top is None:
        l_top = l_max + 30
    RQ = {l_top: np.zeros_like(xi), l_top - 1: np.ones_like(xi)}
    for l in range(l_top - 1, m, -1):
        RQ[l - 1] = ((2 * l + 1) * xi * RQ[l] - (l - m + 1) * RQ[l + 1]) / (l + m)
    scale = q_deriv(m, m, xi) / RQ[m]
    return {l: RQ[l] * scale for l in range(m, l_max + 1)}


def direct(l, m, xi):
    return P.polyval(xi, legendre_deriv_poly(l, m)), q_deriv(l, m, xi)


if __name__ == "__main__":
    # a representative xi grid on (1, ~large), like XiGrid produces
    xi = 1.0 + np.concatenate([np.linspace(1e-3, 2, 40), np.linspace(2.1, 40, 30)])
    def relerr(rec, drc):
        mask = np.isfinite(drc) & (np.abs(drc) > 1e-6 * np.max(np.abs(drc) + 1e-300))
        return float('nan') if not np.any(mask) else \
            float(np.max(np.abs(rec[mask] - drc[mask]) / np.abs(drc[mask])))

    print("=== stable evaluator: forward-P + backward(Miller)-Q vs direct diff (float64) ===")
    print("  m | l  | R^P fwd relerr | R^Q fwd relerr | R^Q BACKWARD relerr")
    worst_pfwd = worst_qbwd = 0.0
    for m in (0, 1, 2, 3, 4):
        RPr, RQfwd = reduced_recur(24, m, xi)
        RQbwd = reduced_Q_backward(24, m, xi)
        for l in (m, m + 4, m + 8, m + 14, m + 20):
            if l > 24:
                continue
            dP, dQ = direct(l, m, xi)
            eP, eQf, eQb = relerr(RPr[l], dP), relerr(RQfwd[l], dQ), relerr(RQbwd[l], dQ)
            worst_pfwd = max(worst_pfwd, eP)
            worst_qbwd = max(worst_qbwd, eQb)
            print(f"  {m} | {l:2d} | {eP:14.3e} | {eQf:14.3e} | {eQb:.3e}")
    print(f"\n  WORST: forward-P = {worst_pfwd:.2e}   backward-Q = {worst_qbwd:.2e}")
    print("  => stable evaluator = forward recurrence for P, backward(Miller) for Q.")
