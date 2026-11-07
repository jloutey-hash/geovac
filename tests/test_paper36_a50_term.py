"""Paper 36: the one-loop alpha (Z alpha)^5 term the headline omits.

    A_50 = 4 pi (139/128 - (1/2) ln 2) ~ 9.29,
    Delta E(nS) = alpha (Z alpha)^5 m c^2 A_50 / (pi n^3)  ~ +6.9 MHz at 2S.

Added 2026-11-07 (group5 guard pass): /qa group5 FULL #1 added this term
to Paper 36's scope statements (abstract, residual section, conclusion)
as the reason the -5.65 MHz residual is not a measure of anything; it had
no backing test.  The A_50 value is the standard Karplus-Klein-Schwinger /
Baranger-Bethe-Feynman coefficient (Eides-Grotch-Shelyuto review); the
test checks the arithmetic of both printed numbers and rejects a missing
1/pi (which would print ~ +21.7 MHz) or a 1/n^2 scaling (~ +13.7 MHz).
"""
import math

ALPHA = 1.0 / 137.035999084
ME_C2_EV = 510998.95000
EV_TO_MHZ = 2.417989242e8
M_P_OVER_M_E = 1836.15267343


def test_a50_value():
    a50 = 4.0 * math.pi * (139.0 / 128.0 - 0.5 * math.log(2.0))
    assert abs(a50 - 9.29) < 0.005, a50


def test_a50_contribution_to_2s_mhz():
    a50 = 4.0 * math.pi * (139.0 / 128.0 - 0.5 * math.log(2.0))
    n, Z = 2, 1
    for reduced in (1.0, (1.0 + 1.0 / M_P_OVER_M_E) ** -3):
        dE = ALPHA * (Z * ALPHA) ** 5 * ME_C2_EV * a50 / (math.pi * n ** 3) * reduced
        mhz = dE * EV_TO_MHZ
        assert 6.8 < mhz < 7.0, mhz
