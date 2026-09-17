"""Diagnostic for #3: can the prolate one-body RADIAL matrix be built DIRECTLY in
the associated-Laguerre basis via recurrences, float64-clean, never forming the
ill-conditioned monomial matrices?

Radial integral for the mu sector (the piece with the e^{-2 a xi} weight, where
conditioning bites):

    r0[i,j] = int_1^inf L_i^{(mu)}(z) L_j^{(mu)}(z) (xi^2-1)^mu e^{-2 a xi} dxi,
    r2[i,j] = same with an extra xi^2,     z = 2 a (xi - 1).

Substituting xi = 1 + t, z = 2 a t:
    (xi^2-1)^mu = (t(t+2))^mu = (2a)^{-2mu} z^mu (z+4a)^mu,   e^{-2a xi} = e^{-2a} e^{-z}.
So r0 = C0 * < L_i, (z+4a)^mu L_j >_{z^mu e^{-z}},  C0 = e^{-2a} (2a)^{-2mu-1}.
(z+4a)^mu is a degree-mu polynomial in z; applying it in the Laguerre basis via the
tridiagonal "multiply-by-z" recurrence and closing with orthogonality gives a BANDED,
float64-clean matrix -- no monomials, no large-coefficient cancellation.

Multiply-by-z:  z L_n^{(mu)} = (2n+mu+1) L_n - (n+1) L_{n+1} - (n+mu) L_{n-1}.
Orthogonality:  <L_i, L_j>_{z^mu e^{-z}} = delta_ij * h_i,  h_i = Gamma(i+mu+1)/i!.

Three builds compared vs mpf-exact ground truth (monomial expansion, exact):
  (A) float64 MONOMIAL  -- expand L to monomials, contract with monomial moments (the bad way)
  (B) float64 RECURRENCE -- the proposed direct build (banded)
Report accuracy and float64 conditioning at growing degree.
"""
from __future__ import annotations

import numpy as np
import mpmath as mp

from geovac import neumann_vee_general_m as ngm
from geovac import prolate_recondition as pr

mp.mp.dps = 50


# ---------- (ground truth, mpf) monomial-expansion radial integral -----------
def r_mpf(n_r: int, mu: int, alpha: float, extra_xi2: bool) -> np.ndarray:
    """r0 (extra_xi2=False) or r2 (True), exact in mpf, via monomial moments.
    L_i^{(mu)} expanded to xi-monomials (pr.assoc_laguerre_coeffs), integrand
    xi^p (xi^2-1)^mu e^{-2 a xi} integrated by the monomial moment recurrence A_n."""
    width = n_r + 1
    Lc = [pr.assoc_laguerre_coeffs(i, mu, alpha, width) for i in range(n_r + 1)]
    xi2m1 = list(ngm._xi2m1_poly(mu))                     # (xi^2-1)^mu, low->high
    A = ngm._mono_moments(2.0 * alpha, 4 * width + 4 * mu + 8)

    def mom(poly):
        return sum((poly[k] * A[k] for k in range(len(poly)) if poly[k] != 0), mp.mpf(0))

    R = np.empty((n_r + 1, n_r + 1), object)
    for i in range(n_r + 1):
        for j in range(i, n_r + 1):
            prod = ngm._polymul(Lc[i], Lc[j])
            prod = ngm._polymul(prod, xi2m1)
            if extra_xi2:
                prod = ngm._shift(prod, 2)                # * xi^2
            R[i, j] = R[j, i] = mom(prod)
    return R


# ---------- (A) float64 monomial build (the bad way) -------------------------
def r_f64_monomial(n_r: int, mu: int, alpha: float, extra_xi2: bool) -> np.ndarray:
    width = n_r + 1
    Lc = [np.array([float(c) for c in pr.assoc_laguerre_coeffs(i, mu, alpha, width)])
          for i in range(n_r + 1)]
    xi2m1 = np.array([float(c) for c in ngm._xi2m1_poly(mu)])
    A = np.array([float(a) for a in ngm._mono_moments(2.0 * alpha, 4 * width + 4 * mu + 8)])

    def mom(poly):
        return float(np.dot(poly, A[:len(poly)]))

    R = np.zeros((n_r + 1, n_r + 1))
    for i in range(n_r + 1):
        for j in range(i, n_r + 1):
            prod = np.polynomial.polynomial.polymul(Lc[i], Lc[j])
            prod = np.polynomial.polynomial.polymul(prod, xi2m1)
            if extra_xi2:
                prod = np.concatenate([[0.0, 0.0], prod])
            R[i, j] = R[j, i] = mom(prod)
    return R


# ---------- (B) float64 recurrence build (the proposed direct way) -----------
def _z_matrix(n: int, mu: int) -> np.ndarray:
    """Tridiagonal 'multiply-by-z' operator on the L^{(mu)} coefficient basis."""
    Z = np.zeros((n, n))
    for k in range(n):
        Z[k, k] = 2 * k + mu + 1
        if k + 1 < n:
            Z[k + 1, k] = -(k + 1)
        if k - 1 >= 0:
            Z[k - 1, k] = -(k + mu)
    return Z


def _apply_poly_in_z(coeffs_lowhigh, Z: np.ndarray) -> np.ndarray:
    """p(Z) for p given by monomial coeffs (low->high), via Horner in the L basis."""
    n = Z.shape[0]
    out = np.zeros((n, n))
    for c in reversed(coeffs_lowhigh):
        out = out @ Z + c * np.eye(n)
    return out


def r_f64_recurrence(n_r: int, mu: int, alpha: float, extra_xi2: bool) -> np.ndarray:
    n = n_r + 1
    a = float(alpha)
    Z = _z_matrix(n + 4, mu)                              # pad for the band to fit
    # polynomial factor in z:  (z + 4a)^mu  [* (1 + z/2a)^2 for r2]
    base = np.array([4.0 * a, 1.0])                       # 4a + z
    P = np.array([1.0])
    for _ in range(mu):
        P = np.polynomial.polynomial.polymul(P, base)
    if extra_xi2:
        q = np.array([1.0, 1.0 / (2.0 * a)])             # 1 + z/(2a)
        P = np.polynomial.polynomial.polymul(P, np.polynomial.polynomial.polymul(q, q))
    M = _apply_poly_in_z(P, Z)[:n, :n]
    h = np.array([float(mp.gamma(i + mu + 1) / mp.factorial(i)) for i in range(n)])
    C0 = float(mp.e ** (-2 * a) * (2 * a) ** (-2 * mu - 1))
    R = C0 * (h[:, None] * M)                             # r[i,j] = C0 h_i M[i,j]
    return 0.5 * (R + R.T)                                # enforce symmetry


def maxrel(F, Rmpf):
    n = F.shape[0]
    d = 0.0
    for i in range(n):
        for j in range(n):
            b = float(Rmpf[i, j])
            if abs(b) > 1e-14:
                d = max(d, abs(F[i, j] - b) / abs(b))
    return d


if __name__ == "__main__":
    alpha = 1.0
    print(f"Radial one-body build, alpha={alpha}.  (A) float64-monomial vs "
          f"(B) float64-recurrence, both vs mpf-exact.\n")
    for mu in (0, 1, 2):
        print(f"--- mu = {mu} ---")
        print(f"  {'n_r':>4} {'relA(mono)':>12} {'relB(recur)':>12} "
              f"{'cond(r0_B)':>11} {'cond(r0_mono)':>13}")
        for n_r in (3, 5, 7, 9, 11):
            for extra in (False, True):
                Rmpf = r_mpf(n_r, mu, alpha, extra)
                Ra = r_f64_monomial(n_r, mu, alpha, extra)
                Rb = r_f64_recurrence(n_r, mu, alpha, extra)
                relA = maxrel(Ra, Rmpf)
                relB = maxrel(Rb, Rmpf)
                if not extra:
                    condB = np.linalg.cond(Rb)
                    condA = np.linalg.cond(Ra)
                    ra_, rb_, cb_, ca_ = relA, relB, condB, condA
            print(f"  {n_r:>4} {ra_:>12.2e} {rb_:>12.2e} {cb_:>11.2e} {ca_:>13.2e}")
        print()
