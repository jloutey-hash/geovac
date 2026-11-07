"""Paper 51 (sec. G3, TT-tensor sector): closed-form TT heat trace.

With spectrum m^2 - 3 (m >= 3) and multiplicity 2(m^2 - 4) on unit S^3,

    K_TT(t) = e^{3t} [ (sqrt(pi)/2) t^{-3/2} - 4 sqrt(pi) t^{-1/2} + 4 ]
              + 6 e^{2t} + O(e^{-pi^2/t}).

Added 2026-11-07 (group5 guard pass): the closed form was rewritten in
/qa group5 FULL #1 (theta_3 content, excluded m = 1, 2 modes) and had no
backing test.  The spectral sum is evaluated directly and compared with
the closed form; the spectrum and multiplicity are the paper's stated
input (not re-derived here).  Rejects a dropped +4 or 6 e^{2t} (the
excluded-mode terms) and a wrong t^{-1/2} coefficient.
"""
import math

import numpy as np


def _k_tt_sum(t: float) -> float:
    m = np.arange(3, 4000, dtype=float)
    return float(np.sum(2.0 * (m * m - 4.0) * np.exp(-(m * m - 3.0) * t)))


def _k_tt_closed(t: float) -> float:
    sp = math.sqrt(math.pi)
    return math.exp(3 * t) * (sp / 2 * t ** -1.5 - 4 * sp * t ** -0.5 + 4.0) + 6.0 * math.exp(2 * t)


def test_k_tt_closed_form_matches_spectral_sum():
    for t in (0.02, 0.05, 0.1, 0.2):
        s, c = _k_tt_sum(t), _k_tt_closed(t)
        # remainder O(e^{-pi^2/t}) is below 1e-20 relative at t <= 0.2
        assert abs(s - c) < 1e-10 * abs(s), (t, s, c)


def test_excluded_mode_terms_are_resolved():
    t = 0.2
    s = _k_tt_sum(t)
    assert abs(s - (_k_tt_closed(t) - 6.0 * math.exp(2 * t))) > 1.0
    assert abs(s - (_k_tt_closed(t) - 4.0 * math.exp(3 * t))) > 1.0
