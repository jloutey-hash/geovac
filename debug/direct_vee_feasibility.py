"""V_ee follow-on, feasibility diagnostic (before building the engine).

Two make-or-break facts for a float64-fast V_ee in the orthogonal basis:

  (1) float64 re-basing of V FAILS -- V carries the monomial dynamic range (memo).
      => a DIRECT orthogonal build is required, not C V_mono C^T in float64.

  (2) the float64-clean route is the LINEARIZATION representation:
      L_a . L_c = sum_k lin_rad[a,c,k] L_k  (moderate coeffs, unlike the monomial
      product coeffs ~1e20), so V_orth assembles in float64 from an mpf-precomputed
      single-degree tensor FL.  This checks the linearization coeffs are moderate
      and that the monomial product coeffs are NOT (the thing to avoid).
"""
from __future__ import annotations
import numpy as np
import mpmath as mp
from geovac import prolate_recondition as pr
from geovac import neumann_vee_general_m as ngm

mp.mp.dps = 40
alpha = 1.0


# ---- (1) does a float64 re-basing of V reproduce the mpf one? ----
def fact1_float64_vee_rebasing_fails(j_max=3, l_max=3, mu_max=1):
    with mp.workdps(40):
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        l_neu = 2 * l_max + 4 * mu_max + 10
        V = pr.vee_mp(fns, alpha, pr.R_DEFAULT, l_neu)          # mpf monomial V
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1) if (l + m) % 2 == 0])
        Tr, Ta = pr._transforms_per_mu("gegenbauer", j_max, l_max, mu_max, alpha)
        V_o_mpf = pr._factored_cob(V, mu_max + 1, Nr, Na, Tr, Ta)   # exact re-basing
        # float64 re-basing: downcast V and transforms first
        Vf = np.array([[float(V[i, j]) for j in range(len(fns))] for i in range(len(fns))])
        Trf = [np.array([[float(x) for x in row] for row in T]) for T in Tr]
        Taf = [np.array([[float(x) for x in row] for row in T]) for T in Ta]
        # dense float64 C V C^T
        C = np.zeros((len(fns), len(fns)))
        for mu in range(mu_max + 1):
            sl = slice(mu * Nr * Na, (mu + 1) * Nr * Na)
            C[sl, sl] = np.kron(Trf[mu], Taf[mu])
        V_o_f64 = C @ Vf @ C.T
        gt = np.array([[float(V_o_mpf[i, j]) for j in range(len(fns))]
                       for i in range(len(fns))])
        rel = np.max(np.abs(V_o_f64 - gt) / np.maximum(np.abs(gt), 1e-12))
    print(f"(1) float64 V re-basing vs mpf: max rel err = {rel:.2e}  "
          f"({'FAILS -> direct build needed' if rel > 1e-3 else 'ok?!'})")


# ---- (2) linearization coeffs moderate vs monomial product coeffs huge ----
def _std_laguerre_z(n, mu, width):
    """Coeffs of L_n^{(mu)}(z) in z^k (standard argument): (-1)^k C(n+mu,n-k)/k!."""
    row = [mp.mpf(0)] * width
    for k in range(n + 1):
        row[k] = (-1) ** k * mp.binomial(n + mu, n - k) / mp.factorial(k)
    return row


def fact2_linearization_is_moderate(n_r=8, mu=0):
    with mp.workdps(40):
        n = n_r + 1
        W = 2 * n + 2
        Lz = [_std_laguerre_z(a, mu, W) for a in range(n)]
        # monomial product coeffs of L_a . L_c in z (the thing to AVOID in float64)
        max_mono = mp.mpf(0)
        for a in range(n):
            for c in range(n):
                for co in ngm._polymul(Lz[a], Lz[c]):
                    max_mono = max(max_mono, abs(co))
        # linearization L_a L_c = sum_k lin[k] L_k, lin[k] = <L_aL_c,L_k>_{z^mu e^-z}/h_k
        max_lin = mp.mpf(0)
        for a in range(n):
            for c in range(a, n):
                prod = ngm._polymul(Lz[a], Lz[c])
                for k in range(a + c + 1):
                    Lk = _std_laguerre_z(k, mu, W)
                    pk = ngm._polymul(prod, Lk)
                    ip = sum((co * mp.factorial(j + mu) for j, co in enumerate(pk)
                              if co != 0), mp.mpf(0))
                    hk = mp.factorial(k + mu) / mp.factorial(k)
                    max_lin = max(max_lin, abs(ip / hk))
    print(f"(2) n_r={n_r}, mu={mu}:  max |monomial product coeff| = {float(max_mono):.2e}  "
          f"(float64 kills this)")
    print(f"                    max |linearization coeff|         = {float(max_lin):.2e}  "
          f"({'MODERATE -> float64-clean path exists' if float(max_lin) < 1e6 else 'large'})")


if __name__ == "__main__":
    fact1_float64_vee_rebasing_fails()
    for nr in (5, 8, 11):
        fact2_linearization_is_moderate(nr, 0)
