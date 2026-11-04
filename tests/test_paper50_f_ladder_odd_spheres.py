"""Paper 50, Remark rem:general_d_conjecture: the conformally coupled scalar's
zeta'(0) on S^7 and S^11 has the printed closed forms in the ring
R_d = {log 2} u {zeta(2k+1)/pi^{2k}}.

Recomputed here from the spectrum, independently of the Door-1 drivers:
eigenvalues v^2 - 1/4 with v = n + (d-1)/2, multiplicity the degree-n harmonic
count (2n+d-1)(n+d-2)!/(n!(d-1)!), which is a polynomial p(v).  Then

    zeta(s) = sum_k binom(-s, k) (-1/4)^k sum_j p_j zeta_H(2s + 2k - j, (d-1)/2),

a binomial Hurwitz series converging like 4^-k, differentiated at s = 0.

Wrong answers rejected: any printed coefficient off (the closed form must match
to ~60 digits, far past the coefficient sizes); a wrong multiplicity or
eigenvalue convention (the S^3 positive control must reproduce the
Klebanov-Pufu-Safdi value F_s = log2/8 - 3 zeta(3)/(16 pi^2)); and a ring
missing an atom (each closed form uses every generator up to zeta(d)/pi^(d-1)).
"""
from __future__ import annotations

import mpmath as mp
import pytest
import sympy as sp

DPS = 80
TOL = mp.mpf(10) ** -60


def _mult_poly(d: int) -> list:
    """Coefficients p_j (ascending powers of v) of the S^d harmonic multiplicity."""
    v = sp.symbols("v")
    n = v - sp.Rational(d - 1, 2)
    m = (2 * n + d - 1) * sp.prod([n + i for i in range(1, d - 1)]) / sp.factorial(d - 1)
    poly = sp.Poly(sp.expand(m), v)
    return [sp.Rational(c) for c in reversed(poly.all_coeffs())]


def zeta_prime_scalar(d: int) -> mp.mpf:
    mp.mp.dps = DPS + 20
    a = mp.mpf(d - 1) / 2
    coeffs = [(j, mp.mpf(c.p) / c.q) for j, c in enumerate(_mult_poly(d)) if c != 0]
    k_max = int((DPS + 20) * mp.log(10) / mp.log(4)) + 2 * d

    def zeta(s):
        total = mp.mpf(0)
        for k in range(k_max):
            b = mp.binomial(-s, k) * mp.mpf(-0.25) ** k
            total += b * mp.fsum(c * mp.zeta(2 * s + 2 * k - j, a) for j, c in coeffs)
        return total

    val = mp.diff(zeta, 0)
    mp.mp.dps = DPS
    return +val


def _ring(d: int, log2_coeff, odd_coeffs) -> mp.mpf:
    """log2_coeff * log 2 + sum_k odd_coeffs[k-1] * zeta(2k+1)/pi^(2k)."""
    out = mp.mpf(log2_coeff.p) / log2_coeff.q * mp.log(2)
    for k, c in enumerate(odd_coeffs, start=1):
        out += mp.mpf(c.p) / c.q * mp.zeta(2 * k + 1) / mp.pi ** (2 * k)
    return out


R = sp.Rational

def test_s3_positive_control_matches_kps():
    """-1/2 zeta'(0) on S^3 = log2/8 - 3 zeta(3)/(16 pi^2)."""
    mp.mp.dps = DPS
    F = -zeta_prime_scalar(3) / 2
    kps = mp.log(2) / 8 - 3 * mp.zeta(3) / (16 * mp.pi ** 2)
    assert abs(F - kps) < TOL, F - kps


@pytest.mark.parametrize("d,log2_coeff,odd_coeffs", [
    (7, R(-1, 512), [R(-41, 15360), R(5, 1024), R(63, 2048)]),
    (11, R(-7, 131072), [R(-3897, 45875200), R(-485, 8257536), R(609, 1310720),
                         R(425, 262144), R(1023, 524288)]),
])
def test_printed_closed_form(d, log2_coeff, odd_coeffs):
    """The printed S^7 / S^11 relations hold to 60 digits, and use the full ring."""
    mp.mp.dps = DPS
    assert len(odd_coeffs) == (d - 1) // 2 and all(c != 0 for c in odd_coeffs)
    val = zeta_prime_scalar(d)
    assert abs(val - _ring(d, log2_coeff, odd_coeffs)) < TOL, val
