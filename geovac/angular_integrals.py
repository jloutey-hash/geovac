"""
Angular integral utilities: Wigner 3j symbols and related coupling coefficients.

The Wigner 3j symbol (j1 j2 j3; m1 m2 m3) is computed via the Racah formula
for integer angular momenta. Results are cached for performance.

References:
  - Racah, Phys. Rev. 62, 438 (1942)
  - Varshalovich, Moskalev, Khersonskii, "Quantum Theory of Angular Momentum"
"""

from functools import lru_cache
from math import factorial, sqrt


@lru_cache(maxsize=4096)
def wigner3j(j1: int, j2: int, j3: int, m1: int, m2: int, m3: int) -> float:
    """
    Compute the Wigner 3j symbol (j1 j2 j3; m1 m2 m3) for integer j values.

    Uses the Racah formula:
        (j1 j2 j3; m1 m2 m3) = (-1)^{j1-j2-m3} * Delta(j1,j2,j3)
            * sqrt((j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!(j3+m3)!(j3-m3)!)
            * sum_t (-1)^t / [t! (j3-j2+t+m1)! (j3-j1+t-m2)!
                              (j1+j2-j3-t)! (j1-t-m1)! (j2-t+m2)!]

    where Delta(a,b,c) = sqrt((a+b-c)!(a-b+c)!(-a+b+c)! / (a+b+c+1)!)

    Parameters
    ----------
    j1, j2, j3 : int
        Angular momentum quantum numbers (non-negative integers).
    m1, m2, m3 : int
        Magnetic quantum numbers, must satisfy |mi| <= ji and m1+m2+m3 = 0.

    Returns
    -------
    float
        Value of the Wigner 3j symbol.
    """
    # Selection rules
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 > j1 + j2 or j3 < abs(j1 - j2):
        return 0.0

    # Triangle coefficient Delta(j1, j2, j3)
    a, b, c = j1, j2, j3
    triangle_num = factorial(a + b - c) * factorial(a - b + c) * factorial(-a + b + c)
    triangle_den = factorial(a + b + c + 1)

    # Prefactor
    prefactor_sq = (
        triangle_num
        * factorial(j1 + m1) * factorial(j1 - m1)
        * factorial(j2 + m2) * factorial(j2 - m2)
        * factorial(j3 + m3) * factorial(j3 - m3)
        / triangle_den
    )

    # Sum over t (Racah sum)
    t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
    t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)

    s = 0.0
    for t in range(t_min, t_max + 1):
        denom = (
            factorial(t)
            * factorial(j3 - j2 + t + m1)
            * factorial(j3 - j1 + t - m2)
            * factorial(j1 + j2 - j3 - t)
            * factorial(j1 - t - m1)
            * factorial(j2 - t + m2)
        )
        s += (-1) ** t / denom

    phase = (-1) ** (j1 - j2 - m3)
    return phase * sqrt(prefactor_sq) * s
