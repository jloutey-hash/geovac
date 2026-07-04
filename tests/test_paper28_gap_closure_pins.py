"""Paper 28 gap-closure pins (2026-07-04, group5 certifying /qa run).

The completeness-critic found Paper 28's Sommerfeld D_k / three-loop
convergence-table region had no regression backing.  Two pin groups:

1. The Sommerfeld product-survival closed forms eq:D4_decomp,
   eq:D5_decomp, eq:D6_decomp evaluated against frozen 30-digit
   decimals.  D4/D5 were independently cross-verified during the review
   by direct partial summation of the paper's Euler-sum coefficient
   formulas (agreement ~1e-17); D6's frozen decimal pins the closed
   form as printed.  A derivation-level pin (the c_p(n) Euler-sum
   machinery, currently archived-only) is a named follow-up
   (docs/qa/group5.done.md).

2. The tab:three_loop_convergence rows, recomputed live from
   geovac/qed_three_loop.py (the printed rows were stale pre-refactor
   values until 2026-07-04).
"""

import pytest
from mpmath import mp, mpf, zeta

from geovac.qed_three_loop import three_loop_factorized


@pytest.fixture(autouse=True)
def _mp_dps_hermetic():
    saved = mp.dps
    yield
    mp.dps = saved


class TestSommerfeldDkClosedForms:
    """eq:D4_decomp / eq:D5_decomp / eq:D6_decomp of Paper 28."""

    def test_paper28_d4_closed_form_value(self) -> None:
        """D4 = -(205/64) z6 + (71/8) z7 - (9/2) z3 z4
        = -0.164140674050572146... (independently cross-verified by
        direct Euler-sum partial summation, 2026-07-04)."""
        mp.dps = 40
        d4 = (-mpf(205) / 64 * zeta(6) + mpf(71) / 8 * zeta(7)
              - mpf(9) / 2 * zeta(3) * zeta(4))
        assert abs(d4 - mpf("-0.1641406740505721461397122729203561439824")) \
            < mpf("1e-35")

    def test_paper28_d5_closed_form_value(self) -> None:
        """D5 = (497/128) z8 - (467/16) z9 + (385/32) z3 z6 + (75/8) z4 z5
        = -0.112930956644175611... (independently cross-verified)."""
        mp.dps = 40
        d5 = (mpf(497) / 128 * zeta(8) - mpf(467) / 16 * zeta(9)
              + mpf(385) / 32 * zeta(3) * zeta(6)
              + mpf(75) / 8 * zeta(4) * zeta(5))
        assert abs(d5 - mpf("-0.1129309566441756113970677383371376466718")) \
            < mpf("1e-35")

    def test_paper28_d6_closed_form_value(self) -> None:
        """D6 = -(2289/512) z10 + (1589/16) z11 - (1617/64) z3 z8
        - (1785/64) z4 z7 - (2065/64) z5 z6 = -0.084201027656549129...
        (pins the closed form as printed; derivation pin = follow-up)."""
        mp.dps = 40
        d6 = (-mpf(2289) / 512 * zeta(10) + mpf(1589) / 16 * zeta(11)
              - mpf(1617) / 64 * zeta(3) * zeta(8)
              - mpf(1785) / 64 * zeta(4) * zeta(7)
              - mpf(2065) / 64 * zeta(5) * zeta(6))
        assert abs(d6 - mpf("-0.0842010276565491292029368589888070306240")) \
            < mpf("1e-35")


class TestThreeLoopConvergenceTable:
    """tab:three_loop_convergence rows (a=4, p=1, CG-weighted),
    recomputed live.  Fast: n_max=50 takes ~2 s."""

    @pytest.mark.parametrize("n_max,expected", [
        (5, 1.061627220971082),
        (10, 3.5720844231573268),
        (15, 5.908257093542889),
        (20, 7.871873785801004),
    ])
    def test_paper28_convergence_row(self, n_max, expected) -> None:
        r = three_loop_factorized(n_max)
        assert abs(r["value_float"] - expected) < 1e-9

    @pytest.mark.slow
    def test_paper28_convergence_row_n50(self) -> None:
        r = three_loop_factorized(50)
        assert abs(r["value_float"] - 14.798403613123499) < 1e-9
