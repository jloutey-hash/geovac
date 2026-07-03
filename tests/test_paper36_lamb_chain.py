"""Paper 36 Lamb-chain reproducibility pins (LS-1 -> LS-6a, 1052.19 MHz).

Backs papers/group5_qed_gauge/paper_36_bound_state_qed.tex (group5 Tier-2
track E build, 2026-07-02).  The original chain lives in archived sprint
drivers under debug/archive/ (LS-1/LS-4/LS-6a in precision_arc, LS-2/LS-3
in qed_arc); this test recomputes the durable sub-chains from permanent
code so the paper's headline is continuously verified:

FAST (always run) -- the LS-6a Eides 3.2 convention assembly:
  closed-form arithmetic from the Drake--Swainson tabulated Bethe logs
  (ln k0(2S) = 2.8117698931, ln k0(2P) = -0.0300167089) through the full
  component table of Paper 36 sec. "Result and component breakdown":
  Lamb(LS-1) = 1025.06 MHz, Lamb(LS-6a) = 1052.19 MHz, shift +27.13 MHz,
  residual -5.65 MHz (-0.534%) vs experimental 1057.845 MHz, plus the
  exact-rational Uehling identity 10/9 - 38/45 = 4/15 (sympy).

SLOW (pytest --slow, ~60 s total) -- the GeoVac-native Bethe logarithms
  recomputed from the Sturmian machinery ported to
  tests/paper36_lamb_support/bethe_log_sturmian.py:
    * LS-4 Drake--Swainson 2P at N=40: ln k0 = -0.0307 (+2.40% vs Drake)
    * LS-4 Drake--Swainson 3D at N=40: ln k0 = -0.005236 (-0.24%)
    * LS-3 acceleration 2S at N=40:    ln k0 =  2.7858  (-0.92%)
    * LS-3 acceleration 1S:            ln k0 =  3.0019  (+0.60%) at N=20
      -- validated 2026-07-02: the 1S acceleration value is a CROSSING
      point at N=20, not a converged value; by N=40 the form has drifted
      to 3.268 (+9.5%).  Both facts are pinned (the paper's sec:drake
      text was corrected accordingly in this build).

Validate-before-pinning record (2026-07-02): every pinned number below was
first reproduced by running the archived drivers directly; the N=55 2P
value (-0.030389, |err| 1.24%, ~62 s) also reproduced but is left
archived-measured to keep the slow budget under ~1 minute.
"""

import math
import sys
from pathlib import Path

import pytest
import sympy as sp

# durable backing module (wh7_support import pattern; see the support README)
sys.path.insert(0, str(Path(__file__).resolve().parent / "paper36_lamb_support"))
from bethe_log_sturmian import (  # noqa: E402
    DRAKE_REF,
    bethe_log_acceleration,
    bethe_log_drake,
)

# ---------------------------------------------------------------------------
# Constants of the chain (identical across LS-1/LS-6a archived drivers)
# ---------------------------------------------------------------------------

ALPHA = 1.0 / 137.035999084          # CODATA 2018
HA_TO_MHZ = 6_579_683_920.502        # MHz / Ha (CODATA)
LAMB_EXP_MHZ = 1057.845              # experimental 2S-2P Lamb shift
BETHE_LOG_2S = 2.8117698931          # Drake & Swainson 1990, ln k0(2,0)
BETHE_LOG_2P = -0.0300167089         # Drake & Swainson 1990, ln k0(2,1)


def _lamb_components(se_2s_const: float) -> dict:
    """One-loop Lamb assembly at n=2, Z=1 with the given 2S SE constant
    (38/45 for LS-1, 10/9 for LS-6a).  Pure closed-form arithmetic."""
    Z, n = 1, 2
    common = ALPHA ** 3 * Z ** 4 / (math.pi * n ** 3)
    se_2s = common * ((4.0 / 3.0) * (math.log(1.0 / (Z * ALPHA) ** 2)
                                     - BETHE_LOG_2S) + se_2s_const)
    se_2p = common * (-(4.0 / 3.0) * BETHE_LOG_2P - 1.0 / 6.0)
    vp_2s = -4.0 * ALPHA ** 3 * Z ** 4 / (15.0 * math.pi * n ** 3)
    lamb = (se_2s + vp_2s) - se_2p
    return {
        "SE_2S_MHz": se_2s * HA_TO_MHZ,
        "SE_2P_MHz": se_2p * HA_TO_MHZ,
        "VP_2S_MHz": vp_2s * HA_TO_MHZ,
        "total_2S_MHz": (se_2s + vp_2s) * HA_TO_MHZ,
        "lamb_MHz": lamb * HA_TO_MHZ,
    }


# ---------------------------------------------------------------------------
# FAST: LS-6a assembly (the 1052.19 MHz headline)
# ---------------------------------------------------------------------------

class TestLS6aAssembly:
    """The full component table of Paper 36 'Result and component
    breakdown', recomputed closed-form.  Tolerance 0.005 MHz against the
    printed 2-decimal values (the assembly is deterministic float
    arithmetic; the archived LS-6a driver reproduces these to 1e-9)."""

    def test_paper36_uehling_constant_identity_exact(self) -> None:
        """Eq. (eq:uehling_kernel): 10/9 - 38/45 = 12/45 = 4/15, exact."""
        diff = sp.Rational(10, 9) - sp.Rational(38, 45)
        assert diff == sp.Rational(12, 45) == sp.Rational(4, 15)

    def test_paper36_ls1_baseline(self) -> None:
        """LS-1 column: SE 2S +1039.31, total 2S +1012.18, Lamb +1025.06."""
        r = _lamb_components(38.0 / 45.0)
        assert abs(r["SE_2S_MHz"] - 1039.31) < 0.005
        assert abs(r["total_2S_MHz"] - 1012.18) < 0.005
        assert abs(r["lamb_MHz"] - 1025.06) < 0.005

    def test_paper36_ls6a_headline(self) -> None:
        """LS-6a column: SE 2S +1066.44, VP 2S -27.13, total 2S +1039.31,
        SE 2P -12.88, Lamb +1052.19 MHz."""
        r = _lamb_components(10.0 / 9.0)
        assert abs(r["SE_2S_MHz"] - 1066.44) < 0.005
        assert abs(r["VP_2S_MHz"] - (-27.13)) < 0.005
        assert abs(r["total_2S_MHz"] - 1039.31) < 0.005
        assert abs(r["SE_2P_MHz"] - (-12.88)) < 0.005
        assert abs(r["lamb_MHz"] - 1052.19) < 0.005

    def test_paper36_convention_shift_is_pure_uehling_constant(self) -> None:
        """The LS-1 -> LS-6a shift is exactly (4/15)(alpha^3/(pi n^3)) in Ha
        = +27.13 MHz (Eq. eq:shift_value), i.e. purely the constant
        difference -- the Bethe-log and VP pieces cancel in the A/B."""
        shift = (_lamb_components(10.0 / 9.0)["lamb_MHz"]
                 - _lamb_components(38.0 / 45.0)["lamb_MHz"])
        predicted = (4.0 / 15.0) * ALPHA ** 3 / (math.pi * 8) * HA_TO_MHZ
        assert abs(shift - predicted) < 1e-9
        assert abs(shift - 27.13) < 0.005

    def test_paper36_residual(self) -> None:
        """Residual vs experiment: -5.65 MHz (-0.534%), the sub-percent
        headline of the abstract."""
        lamb = _lamb_components(10.0 / 9.0)["lamb_MHz"]
        residual = lamb - LAMB_EXP_MHZ
        assert abs(residual - (-5.65)) < 0.005
        assert abs(100.0 * residual / LAMB_EXP_MHZ - (-0.534)) < 0.001
        # sub-percent, and the LS-1 baseline was NOT (-3.10%)
        assert abs(100.0 * residual / LAMB_EXP_MHZ) < 1.0
        lamb_ls1 = _lamb_components(38.0 / 45.0)["lamb_MHz"]
        assert abs(100.0 * (lamb_ls1 - LAMB_EXP_MHZ) / LAMB_EXP_MHZ
                   - (-3.10)) < 0.005


# ---------------------------------------------------------------------------
# SLOW: GeoVac-native Sturmian Bethe logarithms (LS-3 / LS-4)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestBetheLogSturmian:
    """From-scratch Bethe logs on the Sturmian pseudostate basis at
    lam = Z/n, pinning the Paper 36 sec:drake values.  Percent residuals
    are quoted as |computed - Drake| / |Drake| (the paper prints the
    magnitudes)."""

    @staticmethod
    def _pct(value: float, ref: float) -> float:
        return 100.0 * abs(value - ref) / abs(ref)

    def test_paper36_ls4_2p_drake_n40(self) -> None:
        """LS-4 2P at N=40: ln k0 = -0.030737, 2.40% vs Drake--Swainson
        -0.0300167089 (closes the l>0 structural floor)."""
        res = bethe_log_drake(2, 1, N=40, lam=0.5)
        v = res["ln_k0_drake"]
        assert abs(v - (-0.030737)) < 5e-6
        assert abs(self._pct(v, DRAKE_REF[(2, 1)]) - 2.40) < 0.05
        # the velocity-form closure I_v vanishes for l>0 (the structural
        # floor LS-4 regularizes): |I_v| tiny vs the s-state scale 0.25
        assert abs(res["I_v"]) < 1e-3

    def test_paper36_ls4_3d_drake_n40(self) -> None:
        """LS-4 3D at N=40: ln k0 = -0.005236, 0.24% vs -0.0052491189
        (the cleanest from-scratch Bethe log in the framework)."""
        res = bethe_log_drake(3, 2, N=40, lam=1.0 / 3.0)
        v = res["ln_k0_drake"]
        assert abs(v - (-0.005236)) < 5e-6
        assert abs(self._pct(v, DRAKE_REF[(3, 2)]) - 0.24) < 0.05

    def test_paper36_ls3_2s_acceleration_n40(self) -> None:
        """LS-3 acceleration form 2S at N=40: ln k0 = 2.7858, 0.92% vs
        2.8117698931 (3.3x tighter than the velocity form)."""
        res = bethe_log_acceleration(2, 0, N=40, lam=0.5)
        v = res["ln_k0"]
        assert abs(v - 2.785799) < 1e-4
        assert abs(self._pct(v, DRAKE_REF[(2, 0)]) - 0.92) < 0.05

    def test_paper36_ls3_1s_acceleration_crossing(self) -> None:
        """LS-3 acceleration form 1S: the 3.002 (+0.60%) value is the N=20
        CROSSING of a non-converging drift, not a converged value -- by
        N=40 the form has drifted to 3.268 (+9.5%).  Both facts pinned
        (recompute 2026-07-02; Paper 36 sec:drake corrected to attribute
        1S to N=20 and disclose the drift)."""
        res20 = bethe_log_acceleration(1, 0, N=20, lam=1.0)
        v20 = res20["ln_k0"]
        assert abs(v20 - 3.001896) < 1e-4
        assert abs(self._pct(v20, DRAKE_REF[(1, 0)]) - 0.60) < 0.05
        res40 = bethe_log_acceleration(1, 0, N=40, lam=1.0)
        v40 = res40["ln_k0"]
        assert abs(v40 - 3.268050) < 1e-3
        assert self._pct(v40, DRAKE_REF[(1, 0)]) > 5.0  # the drift is real
