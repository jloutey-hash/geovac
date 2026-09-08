"""Closed form for Paper 60's nuclear-diagonal 1-norm  ||T0||_1 = Z sum_nu R_nu.

R_nu = sqrt(1/na^2 + 1/nb^2), summed over the Goscinskian configuration set
{(l, na, nb) : 0 <= l <= lmax, l+1 <= na <= nb <= N}.  No ERI, no radial
integral, no quadrature: a pure combinatorial sum over quantum-number labels.

Claim (derived here, checked numerically below)
-----------------------------------------------
Per l-block, with m = l+1,

    A_m(N) := sum_{m<=a<=b<=N} sqrt(a^-2 + b^-2) = N ln N + c1(m) N + O(log N),

and the coefficient of N ln N is EXACTLY 1, independent of m.  Reason: the sum is
dominated by the region a << b where sqrt(a^-2+b^-2) -> 1/a, giving
sum_a (1/a)(N-a+1) = (N+1) H_N - N ~ N ln N.  The correction from a ~ b is
O(N): writing sqrt(1+(a/b)^2) = 1 + (sqrt(1+(a/b)^2) - 1) and summing the second
piece over b gives a * phi with

    phi = int_0^1 (sqrt(1+u^2) - 1)/u^2 du = 1 - sqrt(2) + ln(1+sqrt(2)),

so sum_a (1/a) * a * phi = phi * N, a pure N term.  Hence

    sum_nu R_nu = (lmax+1) N ln N + O(N),
    K           = sum_{l<=lmax} (N-l)(N-l+1)/2 = (lmax+1) N^2/2 + O(N).

Eliminating N (L := lmax+1):

    ||T0||_1 = Z sum_nu R_nu ~ Z * sqrt(L K / 2) * ln(2K / L),

which for lmax = 3 is  Z sqrt(2K) ln(K/2).  The log-log slope is therefore

    d ln||T0||_1 / d ln K  =  1/2  +  1/ln(2K/L)  +  O(1/ln^2 K)   ->   1/2.

So ||T0||_1 is NOT a power law: it is K^{1/2} times a logarithm.  The fitted
"K^0.70" is 1/2 + 1/ln K evaluated inside the paper's window; its slow decline
(the paper's "flat to falling 0.708 -> 0.697") is the 1/ln K term, not a
converged exponent.
"""
from __future__ import annotations

import numpy as np

PHI = 1.0 - np.sqrt(2.0) + np.log(1.0 + np.sqrt(2.0))
GAMMA = 0.5772156649015328606


def A_exact(m: int, N: int) -> float:
    """sum_{m<=a<=b<=N} sqrt(1/a^2 + 1/b^2), exact double sum (vectorised)."""
    a = np.arange(m, N + 1, dtype=np.float64)
    tot = 0.0
    inv2 = 1.0 / a ** 2
    for i, ai in enumerate(a):
        b = a[i:]
        tot += np.sqrt(inv2[i] + 1.0 / b ** 2).sum()
    return float(tot)


def sum_Rnu(lmax: int, N: int) -> float:
    return sum(A_exact(l + 1, N) for l in range(lmax + 1))


def K_of(lmax: int, N: int) -> int:
    return sum((N - l) * (N - l + 1) // 2 for l in range(lmax + 1))


def predicted_slope(K: int, lmax: int = 3) -> float:
    """1/2 + 1/ln(2K/L), L = lmax+1 -- the leading closed-form local slope."""
    L = lmax + 1
    return 0.5 + 1.0 / np.log(2.0 * K / L)


if __name__ == "__main__":
    print("phi = 1 - sqrt2 + ln(1+sqrt2) = %.12f" % PHI)
    print()
    print("--- leading coefficient of N ln N in A_1(N) (should -> 1) ---")
    print("     N        A_1(N)     (A_1 - N lnN)/N     A_1/(N lnN)")
    for N in (100, 300, 1000, 3000, 10000, 30000):
        A = A_exact(1, N)
        print("%7d %13.4f      %.8f       %.8f"
              % (N, A, (A - N * np.log(N)) / N, A / (N * np.log(N))))
    print()
    print("candidate constants:  gamma+phi-1 = %.8f   " % (GAMMA + PHI - 1.0))
    print()
    print("--- exact ||T0||_1 ladder (lmax=3, Z=2) and its local slope ---")
    print("   N     K     sum R_nu     ||T0||_1    local slope   1/2+1/ln(2K/L)")
    prev = None
    for N in range(7, 31):
        K = K_of(3, N)
        S = sum_Rnu(3, N)
        T0 = 2.0 * S
        if prev is not None:
            sl = np.log(T0 / prev[1]) / np.log(K / prev[0])
            Km = np.sqrt(K * prev[0])
            print("%4d %6d %11.4f %11.4f     %.4f        %.4f"
                  % (N, K, S, T0, sl, predicted_slope(Km, 3)))
        else:
            print("%4d %6d %11.4f %11.4f       --           --" % (N, K, S, T0))
        prev = (K, T0)
    print()
    print("--- global window fits ---")
    for (lo, hi, lab) in ((7, 10, "paper window K=74..164"),
                          (7, 14, "extended K=74..340"),
                          (7, 30, "K=74..1890"),
                          (20, 60, "K=880..7930 (asymptotic probe)")):
        Ks, Ts = [], []
        for N in range(lo, hi + 1):
            Ks.append(K_of(3, N))
            Ts.append(2.0 * sum_Rnu(3, N))
        p = np.polyfit(np.log(Ks), np.log(Ts), 1)[0]
        print("  %-34s  fitted exponent = %.4f" % (lab, p))
