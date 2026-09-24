"""Exact ordered-integral prolate-Neumann radial operator on a fixed Gauss-Legendre xi grid.

Phase 0b of the LiH "marriage" build (2026-09-22; plan debug/lih_marriage_build_plan.md,
memo debug/sprint_lih_marriage_memo.md).  Opt-in: switched on by kernels.USE_EXACT_NEUMANN
(default False keeps the legacy path bit-identical).

THE DEFECT IT REPLACES.  The prolate Neumann expansion of 1/r12 needs, for every target node
xi_i and every (l, m), the radial ORDERED integral of the eta-moment g of the density,

    radial_lm(xi_i) = Q_l^m(xi_i) INT_1^{xi_i} P_l^m(x) g(x) dx
                    + P_l^m(xi_i) INT_{xi_i}^{xi_max} Q_l^m(x) g(x) dx .

hVee.neumann_potential and triangle.coul_mode_potential (legacy) evaluate it as a cumulative
Gauss-Legendre sum over the SAME 72 nodes, K[i, j] = P(xi_<) Q(xi_>): the kernel has a KINK at
x = xi_i, so the GL rule is only O(NXI^-2) accurate there -- +5.0e-3 relative on the Li-1s
self-Coulomb (zeta = 2.69), +8.5e-3 at zeta = 4.5, and +8.09 mHa on the Phase-0 <V_ee>
(NXI 72 -> 144 -> 288 gives +8.09 -> +2.03 -> +0.51 mHa).

THE FIX (exact relative to the interpolant; no grid change).  All sigma densities are entire
in (xi, eta) (r_A = a(xi + eta), r_B = a(xi - eta) are linear), so g(x) is represented by the
degree-(NXI-1) Lagrange interpolant through the GL nodes (barycentric form in the mapped
variable t in [-1, 1]).  MEASURED (debug/lih_marriage_phase0b.py, "(diag)"; 400 off-node points,
l = 0..2): the interpolant error relative to max|g_l| is 1e-15 for 1s densities with zeta <= 2.69
and 1.7e-11..2.2e-11 at zeta = 4.5, where the n = 71 Legendre coefficient of e^{-2 zeta a (xi-1)}
is still 2e-9 of c_0 -- NOT the "~1e-13 by n = 72" the design memo assumed.  The integrals are
nevertheless exact to 1e-13 (5 zeta/8 at zeta = 4.5: 1.65e-13) because that residual is
oscillatory and the moments are weighted toward xi ~ 1.  The two partial integrals are then
evaluated ONCE per (l, m, i) as linear functionals of the node values by composite Gauss sub-rules:

  * P side, [1, xi_i]        : one n_gauss-point Gauss rule.  The integrand is a POLYNOMIAL of
                               degree <= l + m + NXI - 1 (<= 109 for LMAX 34, MMAX 4, NXI 72),
                               so n_gauss = 60 (exact to degree 119) integrates it exactly.
  * Q side, [xi_i, xi_max]   : panels graded geometrically in (x - 1) toward the lower endpoint
                               (ratio 4), n_gauss points each.  The only non-polynomial content
                               is the log((x+1)/(x-1)) of Q_l, singular at x = 1, which sits at
                               distance (xi_i - 1) below the first panel; geometric grading keeps
                               the Bernstein-ellipse parameter >= 2 on every panel, so the log
                               part converges like 2^(-2 n_gauss).

For m >= 1 the mode densities carry rho_cyl^m = a^m [(xi^2-1)(1-eta^2)]^{m/2}, so g(x) has a
(xi^2-1)^{m/2} factor (a square ROOT for odd m) that no polynomial represents.  It is divided
out of the node values and multiplied into the kernel functions:

    P_l^m(x) (x^2-1)^{m/2} = (x^2-1)^m d^m P_l / dx^m       (polynomial; numpy legder)
    Q_l^m(x) (x^2-1)^{m/2}                                  (smooth x log; scipy lqmn)

so the operator A_lm (NXI x NXI) acts on the RAW moments g(xi_j), exactly as the legacy
K.diag(WXI) does:  radial = A_lm @ g.  (A general mode-m density without the rho_cyl^m
factor is unphysical for an azimuthal Fourier coefficient and is NOT supported.)

Conventions match scipy for x > 1 (verified to 1e-15 against numpy legder and mpmath
legenq type 3, 2026-09-22): P_l^m(x) = (x^2-1)^{m/2} d^m P_l/dx^m, Q_l^m from lqmn; both real,
no Condon-Shortley phase -- i.e. the same tables the legacy code uses (lqmn is also more
accurate than the legacy lqn near x -> 1 at large l: 1e-15 vs 2e-10 relative at x = 1.02, l = 34).
Because lqmn is called with m_call = max(mmax, 1), the m = 0 block A[0] of a build with mmax = 0
differs from that of a build with mmax = 4 at the rounding level (3e-15 relative, from lqmn's
m-dependent recurrence); kernels.exact_neumann always builds mmax = kernels._EXACT_MMAX, so
production values are deterministic.

Falsifiers: tests/test_lih_r12ci_neumann_exact.py (fast, 7 tests: 5 zeta/8 at zeta = 4.5 and
2.6875 to 1e-7, the legacy error > 1e-3 pinned, the (aa|bb) Hartree anchor both dressing
directions, the m = 1 solid-harmonic closed form, the default switch off), fire-tested by
debug/firetest_lih_r12ci_neumann_exact.py (dead switches, P/Q interchanged, coarse sub-rule,
default flipped); driver debug/lih_marriage_phase0b.py (m = 0..4 solid-harmonic anchors on
both centres, integrated Richardson of the legacy formula on NXI = 144..1152 scratch grids and a
brute-force no-interpolant pointwise route for the real rho_(NO0,NO1) x rho_cyl^m, gate G0
re-run: <V_ee> grid vs engine +8.09 mHa -> 1e-6 mHa); legacy bit-identity vs `git HEAD` by
debug/lih_marriage_phase0b_headcheck.py.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
from numpy.polynomial import legendre as _leg
from numpy.polynomial.legendre import leggauss
from scipy.special import lqmn

__all__ = ["ExactNeumann", "barycentric_weights", "lagrange_matrix"]


def barycentric_weights(t: np.ndarray) -> np.ndarray:
    """Barycentric weights w_j = 1 / prod_{k != j} (t_j - t_k), in log form; common scale dropped."""
    t = np.asarray(t, dtype=float)
    diff = t[:, None] - t[None, :]
    np.fill_diagonal(diff, 1.0)
    logw = -np.sum(np.log(np.abs(diff)), axis=1)
    sign = np.prod(np.sign(diff), axis=1)
    w = sign * np.exp(logw - logw.max())
    return w


def lagrange_matrix(t_nodes: np.ndarray, w: np.ndarray, t_eval: np.ndarray) -> np.ndarray:
    """Lambda[s, j] = ell_j(t_eval[s]) for the Lagrange basis on t_nodes (barycentric 2nd form)."""
    d = t_eval[:, None] - t_nodes[None, :]
    hit = d == 0.0
    num = w[None, :] / np.where(hit, 1.0, d)
    L = num / np.sum(num, axis=1, keepdims=True)
    rows = np.any(hit, axis=1)
    if np.any(rows):
        L[rows] = hit[rows].astype(float)
    return L


def _legder_coeffs(lmax: int, mmax: int) -> Dict[Tuple[int, int], np.ndarray]:
    """Legendre-series coefficients of d^m P_l / dx^m for 0 <= m <= min(l, mmax)."""
    out: Dict[Tuple[int, int], np.ndarray] = {}
    for l in range(lmax + 1):
        e = np.zeros(l + 1)
        e[l] = 1.0
        for m in range(min(l, mmax) + 1):
            out[(l, m)] = e if m == 0 else _leg.legder(e, m)
    return out


class ExactNeumann:
    """Exact ordered-integral Neumann radial operators A[m, l] (NXI x NXI) on a fixed xi grid.

    Parameters
    ----------
    xi      : the Gauss-Legendre xi nodes of the grid (ascending), on [xi_lo, xi_hi]
    xi_lo   : lower end of the xi interval (1.0 for prolate coordinates)
    xi_hi   : upper end (xi_max of the grid)
    lmax    : highest Legendre degree l
    mmax    : highest azimuthal order m (0 for the sigma / m = 0 potential)
    n_gauss : Gauss points per sub-panel (60 makes the polynomial content exact for l+m+NXI-1 <= 119)
    ratio   : geometric panel ratio in (x - 1) on the Q side

    Attributes
    ----------
    A : ndarray (mmax+1, lmax+1, NXI, NXI);  radial_lm(xi_i) = (A[m, l] @ g)[i], g = raw eta-moments
    n_sub : total number of sub-quadrature points used in the build
    """

    def __init__(self, xi: np.ndarray, xi_lo: float, xi_hi: float, lmax: int, mmax: int = 0,
                 n_gauss: int = 60, ratio: float = 4.0) -> None:
        xi = np.asarray(xi, dtype=float)
        if np.any(np.diff(xi) <= 0.0) or xi[0] <= xi_lo or xi[-1] >= xi_hi:
            raise ValueError("xi must be strictly ascending and strictly inside (xi_lo, xi_hi)")
        if lmax < 0 or mmax < 0 or n_gauss < 2 or ratio <= 1.0:
            raise ValueError("bad build parameters")
        self.xi, self.xi_lo, self.xi_hi = xi, float(xi_lo), float(xi_hi)
        self.lmax, self.mmax, self.n_gauss, self.ratio = int(lmax), int(mmax), int(n_gauss), float(ratio)
        n = xi.size
        self.n = n

        # --- Lagrange basis on the GL nodes, in the mapped variable t in [-1, 1] ---
        scale = 2.0 / (xi_hi - xi_lo)
        t_nodes = scale * (xi - xi_lo) - 1.0
        w_bary = barycentric_weights(t_nodes)

        # --- sub-rules per target node: P side [xi_lo, xi_i], Q side [xi_i, xi_hi] ---
        gx, gw = leggauss(n_gauss)
        xs_all: List[np.ndarray] = []
        ws_all: List[np.ndarray] = []
        side: List[np.ndarray] = []            # 0 = P side, 1 = Q side
        owner: List[np.ndarray] = []           # target index i
        for i in range(n):
            lo, hi = xi_lo, xi[i]
            xP = 0.5 * (hi - lo) * (gx + 1.0) + lo
            wP = 0.5 * (hi - lo) * gw
            bounds = [xi[i]]
            while True:
                nxt = xi_lo + (bounds[-1] - xi_lo) * ratio
                if nxt >= xi_hi:
                    bounds.append(xi_hi)
                    break
                bounds.append(nxt)
            xQ: List[np.ndarray] = []
            wQ: List[np.ndarray] = []
            for b0, b1 in zip(bounds[:-1], bounds[1:]):
                xQ.append(0.5 * (b1 - b0) * (gx + 1.0) + b0)
                wQ.append(0.5 * (b1 - b0) * gw)
            xQa = np.concatenate(xQ)
            wQa = np.concatenate(wQ)
            xs_all += [xP, xQa]
            ws_all += [wP, wQa]
            side += [np.zeros(xP.size, dtype=int), np.ones(xQa.size, dtype=int)]
            owner += [np.full(xP.size, i), np.full(xQa.size, i)]
        xs = np.concatenate(xs_all)
        ws = np.concatenate(ws_all)
        sd = np.concatenate(side)
        ow = np.concatenate(owner)
        self.n_sub = int(xs.size)

        # --- kernel tables on the sub-points and on the target nodes ---
        coeffs = _legder_coeffs(lmax, mmax)
        x2m1 = xs * xs - 1.0
        x2m1_nodes = xi * xi - 1.0
        m_call = max(mmax, 1)                  # lqmn's m>=1 branch is the accurate one (also for m=0)
        Qsub = np.zeros((mmax + 1, lmax + 1, xs.size))
        for s, x in enumerate(xs):
            Qsub[:, :, s] = lqmn(m_call, lmax, x)[0][: mmax + 1]
        Qnode = np.zeros((mmax + 1, lmax + 1, n))
        for i, x in enumerate(xi):
            Qnode[:, :, i] = lqmn(m_call, lmax, x)[0][: mmax + 1]
        Ptil = np.zeros((mmax + 1, lmax + 1, xs.size))      # (x^2-1)^m d^m P_l/dx^m  (polynomial)
        Pnode = np.zeros((mmax + 1, lmax + 1, n))           # P_l^m(xi_i) = (xi^2-1)^{m/2} d^m P_l/dx^m
        Qtil = np.zeros_like(Qsub)                          # (x^2-1)^{m/2} Q_l^m(x)
        for (l, m), c in coeffs.items():
            der = _leg.legval(xs, c)
            Ptil[m, l] = (x2m1 ** m) * der
            Pnode[m, l] = (x2m1_nodes ** (0.5 * m)) * _leg.legval(xi, c)
        for m in range(mmax + 1):
            Qtil[m] = (x2m1 ** (0.5 * m))[None, :] * Qsub[m]
            Qtil[m, :m] = 0.0
            Qnode[m, :m] = 0.0

        # --- assemble A[m, l, i, :] = Q(xi_i) IP_i + P(xi_i) IQ_i, IP/IQ linear in the node values ---
        A = np.zeros((mmax + 1, lmax + 1, n, n))
        for i in range(n):
            selP = np.nonzero((ow == i) & (sd == 0))[0]
            selQ = np.nonzero((ow == i) & (sd == 1))[0]
            LamP = lagrange_matrix(t_nodes, w_bary, scale * (xs[selP] - xi_lo) - 1.0)
            LamQ = lagrange_matrix(t_nodes, w_bary, scale * (xs[selQ] - xi_lo) - 1.0)
            IP = (Ptil[:, :, selP] * ws[selP]) @ LamP        # (mmax+1, lmax+1, n)
            IQ = (Qtil[:, :, selQ] * ws[selQ]) @ LamQ
            A[:, :, i, :] = Qnode[:, :, i][:, :, None] * IP + Pnode[:, :, i][:, :, None] * IQ
        for m in range(1, mmax + 1):
            A[m] /= (x2m1_nodes ** (0.5 * m))[None, None, :]   # act on the raw moments g(xi_j)
            A[m, :m] = 0.0
        self.A = A

    def radial(self, g: np.ndarray, l: int, m: int = 0) -> np.ndarray:
        """radial_lm(xi_i) = sum_j A[m, l, i, j] g(xi_j)  -- the exact replacement of K @ (WXI * g)."""
        return self.A[m, l] @ g
