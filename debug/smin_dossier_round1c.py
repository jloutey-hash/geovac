"""Round-1c: Q1 exactness — CG channel count closed form and the exact
relation between S_min and the true CG-weighted two-loop sunset.

Claim 1 (combinatorial): for n1, n2 >= 1,
    I(n1,n2) = sum_q W(n1,n2,q) = 2*min(n1,n2) - 1 - [n1==n2],
and I = 0 when min(n1,n2) = 0.  Verified exactly for n1,n2 <= 40.

Claim 2 (exact consequence): the true CG-weighted sunset at (a=4, p=1)
    S_sunset = sum_{n1,n2} I(n1,n2) phi(n1) phi(n2),
    phi(n) = g_n/|lambda_n|^4 = 2/(n+3/2)^2 - (1/2)/(n+3/2)^4,
satisfies
    S_sunset = 2*S_min - P^2 - Q,
    P = sum_{n>=1} phi(n) = D(4) - 64/81,    D(4) = pi^2 - pi^4/12,
    Q = sum_{n>=1} phi(n)^2 = 4*zeta(4,5/2) - 2*zeta(6,5/2) + (1/4)*zeta(8,5/2).
So S_min and the sunset differ only by depth-1 / product terms (factor 2).
"""
from fractions import Fraction
import mpmath

mpmath.mp.dps = 80


def so4_count(n1, n2, q):
    j1_L = Fraction(n1 + 1, 2); j1_R = Fraction(n1, 2)
    j2_L = Fraction(n2, 2); j2_R = Fraction(n2 + 1, 2)
    cnt = 0
    for jg_L, jg_R in ((Fraction(q + 1, 2), Fraction(q - 1, 2)),
                       (Fraction(q - 1, 2), Fraction(q + 1, 2))):
        if jg_L < 0 or jg_R < 0:
            continue
        if (abs(j1_L - jg_L) <= j2_L <= j1_L + jg_L
                and abs(j1_R - jg_R) <= j2_R <= j1_R + jg_R):
            cnt += 1
    return cnt


def I_channel(n1, n2):
    return sum(so4_count(n1, n2, q)
               for q in range(max(1, abs(n1 - n2)), n1 + n2 + 1)
               if (n1 + n2 + q) % 2 == 1)


bad = []
for n1 in range(0, 41):
    for n2 in range(0, 41):
        pred = 0 if min(n1, n2) == 0 else (2 * min(n1, n2) - 1
                                           - (1 if n1 == n2 else 0))
        if I_channel(n1, n2) != pred:
            bad.append((n1, n2, I_channel(n1, n2), pred))
print("Claim 1 (I = 2min - 1 - delta), n<=40: %s"
      % ("EXACT, no exceptions" if not bad else "FAILS: %s" % bad[:10]))

# Claim 2 numbers
pi = mpmath.pi
D4 = pi ** 2 - pi ** 4 / 12
P = D4 - mpmath.mpf(64) / 81
Q = (4 * mpmath.hurwitz(4, mpmath.mpf(5) / 2)
     - 2 * mpmath.hurwitz(6, mpmath.mpf(5) / 2)
     + mpmath.mpf(1) / 4 * mpmath.hurwitz(8, mpmath.mpf(5) / 2))
S_min = mpmath.nsum(
    lambda k: (2 * mpmath.hurwitz(2, mpmath.mpf(int(k)) + mpmath.mpf(3) / 2)
               - mpmath.mpf(1) / 2
               * mpmath.hurwitz(4, mpmath.mpf(int(k)) + mpmath.mpf(3) / 2)) ** 2,
    [1, mpmath.inf], method='levin')
S_sunset = 2 * S_min - P ** 2 - Q
print("D(4) = pi^2 - pi^4/12 check: %s"
      % mpmath.nstr(abs(D4 - (2 * mpmath.hurwitz(2, mpmath.mpf(3) / 2)
                              - mpmath.mpf(1) / 2
                              * mpmath.hurwitz(4, mpmath.mpf(3) / 2))), 5))
print("P  = %s" % mpmath.nstr(P, 40))
print("Q  = %s" % mpmath.nstr(Q, 40))
print("S_sunset (exact relation) = %s" % mpmath.nstr(S_sunset, 60))

# independent slow-sum cross-check of S_sunset with tail estimate:
def phi(n):
    a = mpmath.mpf(n) + mpmath.mpf(3) / 2
    return 2 / a ** 2 - mpmath.mpf(1) / 2 / a ** 4


for N in (100, 200, 400):
    part = sum((2 * min(n1, n2) - 1 - (1 if n1 == n2 else 0))
               * phi(n1) * phi(n2)
               for n1 in range(1, N + 1) for n2 in range(1, N + 1))
    print("  direct partial N=%d: %s (diff %s)"
          % (N, mpmath.nstr(part, 20), mpmath.nstr(S_sunset - part, 8)))
