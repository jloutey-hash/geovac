"""Numerical confirmation of the derived closed form for sum_nu R_nu.

Derivation (see debug/sturmian_T0_closed_form.py docstring for the setup):

  A_m(N) := sum_{m<=a<=b<=N} sqrt(a^-2 + b^-2)
          = N ln N + (c1 - H_{m-1}) N + O(log N),

  c1 = gamma + sqrt2 + ln 2 - 2 - ln(1+sqrt2) = -0.1967975...

obtained by writing sqrt(1+(a/b)^2) = 1 + [sqrt(1+(a/b)^2) - 1], summing the
first piece exactly (N H_{N-1} - N ~ N ln N + (gamma-1) N) and the second by
Euler-Maclaurin, which contributes phi*N from the b-integral minus
N * int_0^1 Psi(t) dt from the finite upper limit, with

  phi        = int_0^1 (sqrt(1+u^2)-1)/u^2 du = 1 - sqrt2 + ln(1+sqrt2),
  Psi(t)     = (1-sqrt(1+t^2))/t + asinh(t),
  int_0^1 Psi = phi - (sqrt2 - 1 - ln(1+sqrt2) + ln 2).

For lmax = 3 the configuration count is EXACTLY K = 2N^2 - 4N + 4, so

  ||T0||_1 = Z [ 4 N ln N + (4 c1 - H_0 - H_1 - H_2 - H_3) N + O(log N) ],
  N = 1 + sqrt(K/2 - 1),

i.e. ||T0||_1 ~ Z sqrt(2K) ln(K/2) -- K^{1/2} times a LOG, not a power law.
"""
from __future__ import annotations

import numpy as np

GAMMA = 0.5772156649015328606
C1 = GAMMA + np.sqrt(2) + np.log(2) - 2 - np.log(1 + np.sqrt(2))
H = [0.0]
for i in range(1, 200):
    H.append(H[-1] + 1.0 / i)


def A_exact(m: int, N: int) -> float:
    a = np.arange(m, N + 1, dtype=np.float64)
    inv2 = 1.0 / a ** 2
    tot = 0.0
    for i in range(len(a)):
        tot += np.sqrt(inv2[i] + inv2[i:]).sum()
    return float(tot)


def sum_Rnu(lmax: int, N: int) -> float:
    return sum(A_exact(l + 1, N) for l in range(lmax + 1))


def K_of(N: int, lmax: int = 3) -> int:
    return sum((N - l) * (N - l + 1) // 2 for l in range(lmax + 1))


def model2(N: int, lmax: int = 3) -> float:
    L = lmax + 1
    return L * N * np.log(N) + N * (L * C1 - sum(H[l] for l in range(L)))


if __name__ == "__main__":
    print("c1 (closed form) = %.12f" % C1)
    print()
    print("--- A_1(N): (A_1 - N lnN)/N   ->  c1 ? ---")
    for N in (1000, 3000, 10000, 30000, 60000):
        A = A_exact(1, N)
        print("  N=%6d  (A-NlnN)/N = %.8f   c1 = %.8f   diff*N/lnN = %.4f"
              % (N, (A - N * np.log(N)) / N, C1,
                 ((A - N * np.log(N)) / N - C1) * N / np.log(N)))
    print()
    print("--- K = 2N^2-4N+4 exact check (lmax=3) ---")
    print("  ", all(K_of(N) == 2 * N * N - 4 * N + 4 for N in range(4, 60)))
    print()
    print("--- ||T0||_1 local slope vs the derived asymptote 1/2 ---")
    print("      N        K        exact sum R_nu   2-term model   local slope"
          "   model slope")
    Ns = [10, 14, 20, 30, 50, 80, 130, 200, 320, 500, 800, 1300, 2000, 3200, 5000]
    prev = None
    for N in Ns:
        K = K_of(N)
        S = sum_Rnu(3, N)
        Mm = model2(N)
        if prev is not None:
            sl = np.log(S / prev[1]) / np.log(K / prev[0])
            slm = np.log(Mm / prev[2]) / np.log(K / prev[0])
            print("%7d %9d  %15.4f %14.4f      %.4f       %.4f"
                  % (N, K, S, Mm, sl, slm))
        else:
            print("%7d %9d  %15.4f %14.4f        --           --" % (N, K, S, Mm))
        prev = (K, S, Mm)
