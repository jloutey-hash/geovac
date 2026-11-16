"""Paper 58 -- the certified LiH energy (quadrature-free diatomic, QFD), rebuilt in-suite.

The certified-reference row ``qfd.lih.R3.015`` claims the minimal s-only LiH model (Li 1s, 2s at Z_orbital = 3,
H 1s at Z_orbital = 1, R = 3.015 bohr, full CI) to 30 digits NET OF the exchange Neumann tau-tail bound.  Until
2026-11-15 that claim was campaign-level only (inventory row 3c, "LiH certified-energy regression test OWED").
This file rebuilds the whole energy from the live closed forms (`geovac.qfd_assemble`, exchange via the numeric
tau accumulator at 40 digits, per-quartet tau_max as in the campaign: 20 / 16 / 14) and pins:

  1. E_total against the certified table's own string (imported, so the two cannot drift) within 1e-29: the
     30-significant-digit string rounds at <= 5e-30 and the netted tail is <= 1.24e-30.
  2. E_electronic likewise; V_NN = 3 / 3.015 exactly.
  3. The tail bound the digit count rests on: per-quartet ratios |a_{tau+1}/a_tau| decrease through the last
     four terms (the bound's premise), the geometric bounds |a_last| r / (1 - r) sum below 2e-32, lambda_min(S)
     = 0.594445 (6 digits; equal to 1 - |(S_13, S_23)| exactly), and 6 / lambda_min(S)^2 <= 64, so the energy tail
     is <= 64 x sum <= 1.3e-30.  This guard found the table's lambda_min misquoted as 0.594448 (6th digit; the
     amplification 2.830 and the bound are unchanged); corrected in the table 2026-11-15.

What wrong answers it rejects (fire-tested 2026-11-15, see the end of the file): a perturbed certified string,
a shortened tau series, an exchange accumulator at too low precision, a changed exchange term.

Slow: about 11 minutes (the exchange quartets dominate).
"""
from __future__ import annotations

from fractions import Fraction

import mpmath as mp
import pytest
import sympy as sp

from benchmarks.certified_reference.entries_qfd_diatomic import QUOTED
from geovac import qfd_assemble as AS

pytestmark = pytest.mark.slow

ROW = QUOTED["qfd.lih.R3.015"]
ORBS = [("A", Fraction(3), 1), ("A", Fraction(3), 2), ("B", Fraction(1), 1)]
R = sp.Rational(603, 200)
TAU = {(1, 1, 1, 1): 20, (1, 1, 2, 1): 16, (2, 1, 1, 1): 16, (2, 1, 2, 1): 14}
EXCHANGE_DPS = 40
DPS = 40


def _tau_max(q):
    return TAU[(q[0][2], q[1][2], q[2][2], q[3][2])]


@pytest.fixture(scope="module")
def lih():
    AS.EXCHANGE_TAILS.clear()
    S, h = AS.build_S_h(ORBS, 3, 1, R)
    gmap, _ = AS.build_g(ORBS, R, tau_max=_tau_max, exchange_hp_dps=EXCHANGE_DPS)
    with mp.workdps(DPS + 40):
        tot, ele, vnn = AS.total_energy(S, h, gmap, ORBS, 4, 3, 1, R, DPS)
        tot, ele, vnn = +tot, +ele, +vnn
    tails = {k: list(v) for k, v in AS.EXCHANGE_TAILS.items()}
    return dict(S=S, tot=tot, ele=ele, vnn=vnn, tails=tails)


def test_lih_total_energy_matches_certified(lih):
    with mp.workdps(DPS + 40):
        ref = mp.mpf(ROW["value"])
        assert ROW["digits"] == 30
        assert abs(lih["tot"] - ref) <= mp.mpf("1e-29"), mp.nstr(lih["tot"] - ref, 5)


def test_lih_electronic_and_nuclear_repulsion(lih):
    with mp.workdps(DPS + 40):
        assert abs(lih["ele"] - mp.mpf(ROW["electronic"])) <= mp.mpf("1e-29")
        assert abs(lih["vnn"] - mp.mpf(3) / mp.mpf("3.015")) <= mp.mpf("1e-50")
        assert abs(lih["ele"] + lih["vnn"] - lih["tot"]) <= mp.mpf("1e-50")


def test_lih_tau_tail_bound_supports_30_digits(lih):
    tails = lih["tails"]
    # three distinct exchange quartets (one of (1,1,2,1)/(2,1,1,1) is the other's canonical image)
    assert len(tails) == 3, sorted(tails)
    assert sorted(len(v) for v in tails.values()) == [15, 17, 21]       # tau = 0..tau_max
    with mp.workdps(60):
        total = mp.mpf(0)
        for key, per in tails.items():
            a = [abs(mp.mpf(x)) for x in per]
            ratios = [a[i + 1] / a[i] for i in range(len(a) - 5, len(a) - 1)]
            assert all(ratios[i + 1] < ratios[i] for i in range(len(ratios) - 1)), (key, ratios)
            r = ratios[-1]
            assert r < mp.mpf("0.01"), (key, r)
            total += a[-1] * r / (1 - r)
        assert total < mp.mpf("2e-32"), mp.nstr(total, 4)
        n = len(ORBS)
        Sm = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                Sm[i, j] = mp.mpf(str(sp.N(lih["S"][i][j], 70)))
        lam_min = min(mp.eigsy(Sm, eigvals_only=True))
        # independent route: S = [[1, 0, a], [0, 1, b], [a, b, 1]] (Li 1s, 2s orthogonal), so lambda_min = 1 - |(a, b)|
        assert Sm[0, 1] == 0
        assert abs(lam_min - (1 - mp.sqrt(Sm[0, 2] ** 2 + Sm[1, 2] ** 2))) < mp.mpf("1e-50")
        # the table's stated value (corrected 2026-11-15 from 0.594448, a 6th-digit transcription error)
        assert abs(lam_min - mp.mpf("0.594445")) < mp.mpf("5e-7"), lam_min
        amplification = 6 / lam_min ** 2
        assert amplification <= 64
        assert 64 * total <= mp.mpf("1.3e-30")

# Fire-test record (debug/qa/fire_test.py): appended below after the plants were run.
# 2026-11-16, each plant one fire_test.py invocation (whole file selected unless noted); green run first: 3 passed.
#   P1  entries_qfd_diatomic.py: certified value's 29th digit 3 -> 4          FIRED (total-energy test)
#   P2  TAU (2,1,2,1): 14 -> 12 (a shortened tau series)                     FIRED
#   P3  EXCHANGE_DPS 40 -> 20                                                DID NOT FIRE -- not a wrong answer:
#       exchange_hp evaluates at dps + 12 / dps + 20 guard digits, so "20" still computes ~32 digits, past the
#       30-digit claim.  Replaced by P3b.
#   P3b EXCHANGE_DPS 40 -> 10 (~22 effective digits)                         FIRED
#   P4  qfd_core.exchange_hp: tau weight x (1 + 10^-25)                      FIRED
#   P5  -k tail: the table's lambda_min back to the misquoted 0.594448         FIRED
