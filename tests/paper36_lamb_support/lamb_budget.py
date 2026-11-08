"""Order-by-order budget of the hydrogen 2S1/2-2P1/2 Lamb shift (Paper 36, Sec. "The residual").

Every 2S entry is a number printed in Eides, Grotch and Shelyuto, Phys. Rep. 342, 63 (2001),
arXiv:hep-ph/0002158 (table numbers as in the arXiv version).  The review prints no 2P1/2
column, so the 2P1/2 entries are evaluated here from its formulas.  All values in kHz.
"""
from math import log, pi
from typing import Dict

ALPHA = 1 / 137.03599958          # the review's alpha (Eq. (385) era)
M_E_OVER_M_P = 1 / 1836.1526675
ME_C2_KHZ = 1.2355899651e17       # m_e c^2 / h in kHz
LNK0_2S = 2.811769893120563       # Table I
LNK0_2P = -0.030016708630213      # Table I
EIDES_RP_FM = 0.862               # proton radius behind Table X and Table XX
EIDES_TABLE_XX_KHZ = 1057833.0    # theory at r_p = 0.862 fm, error (2)(4)

N = 2
MR = 1 / (1 + M_E_OVER_M_P)       # m_r / m
AP = ALPHA / pi
LOG_ZA = log(1 / (MR * ALPHA**2))  # ln[(m/m_r)(Z alpha)^-2]


def c4(reduced_mass: bool = True) -> float:
    """Table II prefactor 4 alpha (Z alpha)^4 / (pi n^3) (m_r/m)^3 m."""
    return 4 * ALPHA**5 / (pi * N**3) * (MR**3 if reduced_mass else 1.0) * ME_C2_KHZ


def c6() -> float:
    """Table V prefactor 4 alpha (Z alpha)^6 / (pi n^3) (m_r/m)^3 m."""
    return 4 * ALPHA**7 / (pi * N**3) * MR**3 * ME_C2_KHZ


def one_loop_4(reduced_mass: bool = True) -> float:
    """alpha (Z alpha)^4 one loop, SE + Uehling, 2S minus 2P1/2 (Table II)."""
    c = c4(reduced_mass)
    L = LOG_ZA if reduced_mass else log(1 / ALPHA**2)
    pauli_2p = -1 / 24 / (MR if reduced_mass else 1.0)
    level_2s = c * ((1 / 3) * L + 11 / 72 - (1 / 3) * LNK0_2S + 1 / 8 - 1 / 15)
    level_2p = c * (-(1 / 3) * LNK0_2P + pauli_2p)
    return level_2s - level_2p


def multi_loop_4() -> float:
    """alpha^n (Z alpha)^4, n = 2, 3, plus muonic and hadronic VP (Table II)."""
    two_2s = 443.48 - 77.50 - 238.83
    three_2s = -0.79 + 0.38 + 0.65
    mu_had_2s = -0.63 - 0.40
    # 2P1/2 anomalous-moment terms: Pauli l != 0 coefficient x (j(j+1)-l(l+1)-3/4)/(l(l+1)(2l+1)) = -1/3
    pauli_2p = c4() / MR * (-1 / 3) * (-0.0821197 * AP + 0.2953103 * AP**2)
    return two_2s + three_2s + mu_had_2s - pauli_2p


def one_loop_5() -> float:
    """alpha (Z alpha)^5, S states only (Table III)."""
    return 6886.29 + 242.55


def two_loop_5_6() -> float:
    """alpha^2 (Z alpha)^5 (Table III) and alpha^2 (Z alpha)^6 (Table V); the +-2 kHz log^2 term is unknown."""
    two5 = -0.33 + 2.75 + 3.31 - 0.39 - 0.66 - 41.78
    two6_2s = -3.55 + 0.22 - 0.06
    two6_2p = c6() * (N**2 - 1) / (27 * N**2) * AP * log(ALPHA**-2) ** 2
    return two5 + two6_2s - two6_2p


def one_loop_6() -> float:
    """alpha (Z alpha)^6 (Tables IV, V).  2P1/2: A61 = 103/180 log term, A60, l = 1 polarization."""
    level_2s = -208.16 - 172.21 - 7.10 - 4.47 + 0.31
    log_coef = (1 - 1 / N**2) * (1 / 30 + 1 / 12) + (6 - 2 * 2 / N**2) / (3 * 5 * 2 * 3)
    vp_2p = -(1 / 15) * (1 - 1 / N**2) * (1 / 14 + 1 / 4)
    level_2p = c6() * (log_coef * LOG_ZA + vp_2p) - 5.40
    return level_2s - level_2p


def one_loop_7() -> float:
    """alpha (Z alpha)^7: 2S from Table VII, 2P1/2 G_SE,7 from Table VI."""
    return (1.80 - 0.09 + 0.06 + 0.03 - 0.01) - 0.12


def recoil() -> float:
    """Table VIII recoil, plus the Barker-Glover (Z alpha)^4 m^2/M^2 term for 2P1/2."""
    level_2s = -73.75 - 305.46 - 0.06 + 720.56 - 0.92 - 0.05
    p = ALPHA**5 / (pi * N**3) * MR**3 * M_E_OVER_M_P * ME_C2_KHZ
    tt = -1 / 6                                   # transverse-transverse, l = 1
    ct = (8 / 3) * (-LNK0_2P - 1 / 12)            # Coulomb-transverse, l = 1
    gemk = 0.4 * (1 - 2 / 12) * pi * ALPHA        # (Z alpha)^6 m/M, l = 1
    level_2p = p * (tt + ct + gemk)
    barker_glover_2p = ALPHA**4 * MR**3 * M_E_OVER_M_P**2 * ME_C2_KHZ / (2 * N**3) * (1 / 1 - 1 / 1.5)
    return level_2s - level_2p - barker_glover_2p


def radiative_recoil() -> float:
    """Table IX, with Eq. (148) for alpha (Z alpha)^5 as Table XX's footnote 38 states."""
    return -2.52 - 0.05 + 0.58 + 0.19


def finite_size(rp_fm: float) -> float:
    """Table X (S states; the 2P1/2 terms are below 1 Hz).  r^2-proportional rows rescaled from 0.862 fm."""
    s = (rp_fm / EIDES_RP_FM) ** 2
    return 145.0 * s - 0.004 - 0.009 + 0.095 * s - 0.023 + 0.006


def budget(rp_fm: float = 0.8409) -> Dict[str, float]:
    """The rows of Paper 36's budget table, in kHz."""
    return {
        "one_loop_4_infinite_mass": one_loop_4(reduced_mass=False),
        "reduced_mass": one_loop_4() - one_loop_4(reduced_mass=False),
        "one_loop_5": one_loop_5(),
        "one_loop_6": one_loop_6(),
        "one_loop_7": one_loop_7(),
        "multi_loop_4": multi_loop_4(),
        "two_loop_5_6": two_loop_5_6(),
        "recoil": recoil(),
        "radiative_recoil": radiative_recoil(),
        "finite_size": finite_size(rp_fm),
    }


def total(rp_fm: float = 0.8409) -> float:
    return sum(budget(rp_fm).values())


if __name__ == "__main__":
    for k, v in budget().items():
        print(f"{k:28s} {v / 1e3:+12.3f} MHz")
    print(f"total r_p=0.8409 {total() / 1e3:.4f} MHz;  r_p=0.862 {total(EIDES_RP_FM) / 1e3:.4f} MHz")
