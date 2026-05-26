"""Sanity check: verify S^5 scalar matches v3.2.1 value, then confirm S^7 anomaly.

S^5 expected: zeta'(0) = log(2)/16 + zeta(3)/(16 pi^2) - 15 zeta(5)/(32 pi^4)
              ~~ 0.04594

S^7 my value: -0.001594... (no PSLQ relation in natural rings)

If S^5 reproduces the v3.2.1 bit-exact value, then the Hurwitz machinery is
right and the S^7 negative is structural. If S^5 doesn't match, there's a
bug in our Hurwitz expansion across the d-family.
"""

import mpmath as mp
import sympy as sp

mp.mp.dps = 100


def zeta_Sd_conf_prime_at_zero(d, k_max=100):
    """Conformal scalar zeta'(0) on round S^d for odd d.

    Spectrum: (n+(d-2)/2)(n+d/2) for n=0,1,2,...
    Re-indexed: eigenvalue v^2 - 1/4 with v = n + (d-1)/2, v=(d-1)/2, (d+1)/2, ...

    Multiplicity (polynomial in v): degree (d-2). For S^3 (d=3): v.  For S^5
    (d=5): v^2(v^2-1)/3 = (v^4-v^2)/3. For S^7 (d=7): v^2(v^2-1)(v^2-4)/360
                                                    = (v^6-5v^4+4v^2)/360.

    The zeta'(0) formula:
       zeta(s) = sum_k g_k(s) [appropriate combinations of zeta_R(2s+2k-...)]
    """
    if d == 3:
        # eigenvalue v^2 - 1/4, multiplicity v (re-indexed from (n+1)^2 with v=n+1)
        # Multiplicity polynomial: v (degree 1; eigenvalues start at v=1)
        # zeta(s) = sum g_k(s) [zeta_R(2s+2k-1) - ... hmm wait]
        # Actually S^3 multiplicity is (n+1)^2 in n indexing.
        # In v = n+1 indexing: v=1, 2, 3, ..., multiplicity = v^2
        # zeta(s) = sum_v v^2 / (v^2 - 1/4)^s
        #        = sum_k g_k(s) sum_v v^{2-2s-2k}
        #        = sum_k g_k(s) zeta_R(2s+2k-2)
        # Then zeta'(0) = 2 zeta_R'(-2) + sum_{k>=1} (1/(k 4^k)) zeta_R(2k-2)
        leading = 2 * mp.zeta(-2, derivative=1)
        series_sum = mp.mpf(0)
        for k in range(1, k_max + 1):
            term = mp.zeta(2*k - 2) / (k * mp.mpf(4)**k)
            series_sum += term
        return leading + series_sum
    elif d == 5:
        # eigenvalue v^2 - 1/4, multiplicity v^2(v^2-1)/3 = (v^4 - v^2)/3
        # zeta(s) = (1/3) sum_k g_k(s) [zeta_R(2s+2k-4) - zeta_R(2s+2k-2)]
        # zeta'(0) = (1/3) {2 [zeta_R'(-4) - zeta_R'(-2)] + sum_{k>=1} (1/(k 4^k)) [zeta_R(2k-4) - zeta_R(2k-2)]}
        leading = 2 * (mp.zeta(-4, derivative=1) - mp.zeta(-2, derivative=1))
        series_sum = mp.mpf(0)
        for k in range(1, k_max + 1):
            term = (mp.zeta(2*k - 4) - mp.zeta(2*k - 2)) / (k * mp.mpf(4)**k)
            series_sum += term
        return (leading + series_sum) / 3
    elif d == 7:
        # eigenvalue v^2 - 1/4, multiplicity (v^6 - 5v^4 + 4v^2)/360
        # zeta(s) = (1/360) sum_k g_k(s) [zeta_R(2s+2k-6) - 5 zeta_R(2s+2k-4) + 4 zeta_R(2s+2k-2)]
        leading = 2 * (mp.zeta(-6, derivative=1) - 5 * mp.zeta(-4, derivative=1) + 4 * mp.zeta(-2, derivative=1))
        series_sum = mp.mpf(0)
        for k in range(1, k_max + 1):
            term = (mp.zeta(2*k - 6) - 5 * mp.zeta(2*k - 4) + 4 * mp.zeta(2*k - 2)) / (k * mp.mpf(4)**k)
            series_sum += term
        return (leading + series_sum) / 360
    else:
        raise ValueError(f"d={d} not implemented")


# Sanity check: S^3 and S^5 against v3.2.1 published values
print("=" * 60)
print("Sanity check: S^3 and S^5 conformal scalar zeta'(0)")
print("=" * 60)
print()

# S^3 expected
s3 = zeta_Sd_conf_prime_at_zero(3, k_max=100)
s3_expected = -mp.log(2)/4 + 3 * mp.zeta(3) / (8 * mp.pi**2)
print(f"S^3 framework:  {mp.nstr(s3, 30)}")
print(f"S^3 expected:   {mp.nstr(s3_expected, 30)}")
print(f"Difference:     {mp.nstr(abs(s3 - s3_expected), 5)}")
print(f"  (S^3 expected = -log(2)/4 + 3 zeta(3)/(8 pi^2) from Paper 50 Sec 3)")
print()

s5 = zeta_Sd_conf_prime_at_zero(5, k_max=100)
s5_expected = mp.log(2)/16 + mp.zeta(3)/(16 * mp.pi**2) - 15 * mp.zeta(5) / (32 * mp.pi**4)
print(f"S^5 framework:  {mp.nstr(s5, 30)}")
print(f"S^5 expected:   {mp.nstr(s5_expected, 30)}")
print(f"Difference:     {mp.nstr(abs(s5 - s5_expected), 5)}")
print(f"  (S^5 expected = log(2)/16 + zeta(3)/(16 pi^2) - 15 zeta(5)/(32 pi^4)")
print(f"   from Paper 50 Theorem 7.1 PSLQ relation [32, -2, -2, 15])")
print()

s7 = zeta_Sd_conf_prime_at_zero(7, k_max=100)
print(f"S^7 framework: {mp.nstr(s7, 30)}")
print(f"  No PSLQ relation found in any natural ring tested in this session.")
print(f"  If S^3 and S^5 above match (bit-exact), the Hurwitz machinery is correct")
print(f"  and the S^7 non-match is structural, NOT computational.")
