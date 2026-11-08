"""V_ee z-recurrence, step 1: the first-kind Neumann moment A_l[a] in the LAGUERRE
index, built by the associated-Legendre l-recurrence with the small-coeff Xi
operator -- float64-clean, never expanding P_l or L_a into large monomials.

A_l^{m,s}[a] = <L_a | (xi^2-1)^s d^m P_l / dxi^m>_{e^{-c xi}},  c = 2 alpha,
             L_a = L_a(c(xi-1)) (the re-basing radial functions).

Recurrence (assoc-Legendre, spectator weight (xi^2-1)^s):
   (l-m+1) A_{l+1} = (2l+1) Xi . A_l - (l+m) A_{l-1},   Xi = I + Z/c  (multiply by xi),
seeded at l = m, m+1.  Z is the symmetric tridiagonal multiply-by-z operator on
L^{(0)} (diag 2j+1, offdiag -(j+1)).

Validated 3 ways: (A) mpf ground truth via the monomial A-moment contraction
(exact); (B) the float64 RECURRENCE; (C) the float64 MONOMIAL contraction (the bad
route, to show it degrades where the recurrence stays exact).
"""
from __future__ import annotations
import numpy as np
import mpmath as mp
from geovac import prolate_recondition as pr
from geovac import neumann_vee_general_m as ngm

mp.mp.dps = 40


def _Z(n):
    Z = np.zeros((n, n))
    for j in range(n):
        Z[j, j] = 2 * j + 1
        if j + 1 < n:
            Z[j, j + 1] = Z[j + 1, j] = -(j + 1)
    return Z


def A_ground_truth(n_r, m, s, l_max, alpha, dtype):
    """A_l[a] = sum_p Lc_a[p] * int xi^p (xi^2-1)^s d^mP_l e^{-c xi} dxi.
    dtype=object -> mpf exact; dtype=float -> float64 (the bad monomial route)."""
    n = n_r + 1
    c = mp.mpf(2.0 * alpha)
    p_max = n_r + 2 * s + 2
    Amono = ngm._mono_moments(c, p_max + 2 * l_max + 4)
    Lc = [pr.laguerre_coeffs(a, alpha, n) for a in range(n)]   # L_a(2a(xi-1)) in xi
    out = {}
    for l in range(m, l_max + 1):
        Av = [ngm._A_moment(l, m, s, p, Amono) for p in range(p_max + 1)]
        vec = np.empty(n, object if dtype is object else float)
        for a in range(n):
            val = mp.mpf(0)
            for p in range(min(len(Lc[a]), p_max + 1)):
                if Lc[a][p] != 0:
                    val += Lc[a][p] * Av[p]
            vec[a] = val if dtype is object else float(val)
        out[l] = vec
    return out


def A_recurrence(n_r, m, s, l_max, alpha, pad=None):
    """A_l[a] via the float64 l-recurrence, mpf seeds at l=m, m+1.

    Xi couples a<->a+-1, so each l-step leaks the boundary; build PADDED (extra
    radial indices) and truncate to n_r+1 at the end -- pad >= l_max-m keeps the
    kept indices exact through l_max.
    """
    n = n_r + 1
    if pad is None:
        pad = l_max - m + 2
    npad = n_r + pad                                    # padded radial size (index 0..npad)
    c = 2.0 * alpha
    Xi = np.eye(npad + 1) + _Z(npad + 1) / c
    gt_seed = A_ground_truth(npad, m, s, m + 1, alpha, object)  # exact seeds, padded
    A = {m: np.array([float(x) for x in gt_seed[m]])}
    if m + 1 <= l_max:
        A[m + 1] = np.array([float(x) for x in gt_seed[m + 1]])
    for l in range(m + 1, l_max):
        A[l + 1] = ((2 * l + 1) * (Xi @ A[l]) - (l + m) * A[l - 1]) / (l - m + 1)
    return {l: A[l][:n] for l in A}                     # truncate to the basis size


def B_ground_truth(n_r, m, s, l_max, alpha):
    """B_l[a] = <L_a|(xi^2-1)^s d^mQ_l>_{e^{-c xi}}, mpf exact via ngm._B_table."""
    n = n_r + 1
    c = mp.mpf(2.0 * alpha)
    p_max = n_r + 2 * s + 2
    Bt = ngm._B_table(m, s, l_max, p_max, c)            # B[(l,p)] mpf
    Lc = [pr.laguerre_coeffs(a, alpha, n) for a in range(n)]
    out = {}
    for l in range(m, l_max + 1):
        vec = np.empty(n, object)
        for a in range(n):
            vec[a] = sum((Lc[a][p] * Bt[(l, p)] for p in range(min(len(Lc[a]), p_max + 1))
                          if Lc[a][p] != 0), mp.mpf(0))
        out[l] = vec
    return out


def B_recurrence(n_r, m, s, l_max, alpha, pad=None):
    """B_l[a] via the float64 forward l-recurrence (same shape as A), padded seeds."""
    n = n_r + 1
    if pad is None:
        pad = l_max - m + 2
    npad = n_r + pad
    c = 2.0 * alpha
    Xi = np.eye(npad + 1) + _Z(npad + 1) / c
    seed = B_ground_truth(npad, m, s, m + 1, alpha)     # exact seeds l=m, m+1
    B = {m: np.array([float(x) for x in seed[m]])}
    if m + 1 <= l_max:
        B[m + 1] = np.array([float(x) for x in seed[m + 1]])
    for l in range(m + 1, l_max):
        B[l + 1] = ((2 * l + 1) * (Xi @ B[l]) - (l + m) * B[l - 1]) / (l - m + 1)
    return {l: B[l][:n] for l in B}


def relerr(vec_f, vec_mpf):
    d = 0.0
    for a in range(len(vec_f)):
        b = float(vec_mpf[a])
        if abs(b) > 1e-14:
            d = max(d, abs(vec_f[a] - b) / abs(b))
    return d


if __name__ == "__main__":
    alpha = 1.0
    for (m, s) in [(0, 0), (1, 1), (2, 2)]:
        n_r, l_max = 11, 14
        gt = A_ground_truth(n_r, m, s, l_max, alpha, object)   # mpf exact
        bad = A_ground_truth(n_r, m, s, l_max, alpha, float)   # float64 monomial
        rec = A_recurrence(n_r, m, s, l_max, alpha)            # float64 recurrence
        print(f"m={m}, s={s}, n_r={n_r}:  A_l[a] relerr vs mpf, l = {m}..{l_max}")
        print(f"  {'l':>3} {'recurrence(f64)':>16} {'monomial(f64)':>16}")
        for l in range(m, l_max + 1):
            rr = relerr(rec[l], gt[l]) if l in rec else float('nan')
            rb = relerr(bad[l], gt[l])
            print(f"  {l:>3} {rr:>16.2e} {rb:>16.2e}")
        print()
