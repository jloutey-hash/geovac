# -*- coding: utf-8 -*-
"""Frozen falsifier for the B3 Phase-3 band-exhaustion program (2026-06-10).

Multi-agent workflow (3 probes + rate-ID + adversarial; drivers
tests/wh7_support/wh7_band_exh_{legs,penalties,intervals,rates,adversarial}.py); this
falsifier RECOMPUTES the load-bearing facts directly from the shared substrate
library (tests/wh7_support/wh7_band_exhaustion_lib.py), independent of the agent drivers.

Pins:
  (1) the EXACT NORM CLOSED FORM: the folded null-class compression has
      ||H(j_max)||_2 = sqrt(6) * (2 j_max) / (2 j_max + 2)  at every window
      (limit sqrt(6) = sqrt(b(b+1)(2b+1)) at b = 1 -- algebraic, pi-free);
  (2) the b-PARITY STAIRCASE: integer-b cost data freezes bit-exactly across
      half-integer shell additions (c12 identical on the (1, 3/2) and
      (2, 5/2) window pairs, real jumps only at half-integer -> integer
      steps) -- the (-1)^{2b} spin-statistics grading appearing in the
      band-exhaustion dynamics;
  (3) the FLOW-TRANSLATION identity c12 = c23 at every window (exact);
  (4) RATE EXCLUSIONS: the mixed-reference (b = 1/2) deficit converges with
      strictly shrinking successive diffs whose ratios sit in (0.55, 0.85)
      -- excluding BOTH the metric-layer gamma rate (predicted ratios
      0.93-0.95) and the thermal rate (constant e^{-1} = 0.368); the same
      band holds at beta = 0.5 (the thermal-suppression attack fails:
      convergence is structural, not KMS tail suppression);
  (5) INTERVAL-LAYER STABILITY: operational interval recovery and additivity
      at machine precision at both the smallest and largest windows;
  (6) the ADMISSIBILITY TABLE per window (Sprint-3b window-edge content):
      mu'=0 classes commute at every window; (2,2) commutes only at
      j_max = 1; (2,1) is zero only at j_max = 1;
  (7) the QUALIFIED penalty statement (adversarial-corrected): the commuting
      (1,0) class has machine-zero excess at every window; non-commuting
      penalties do NOT stabilize monotonically across the ladder.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "wh7_support"))
import wh7_band_exhaustion_lib as lib  # noqa: E402

LADDER = lib.JMAX_LADDER                      # [1, 3/2, 2, 5/2, 3]


@pytest.fixture(scope="module")
def legs():
    """c12/c23/c13/deficit ladders for both references at beta in {1, 0.5}."""
    out = {}
    for ref in ("null", "mixed"):
        for beta in (1.0, 0.5):
            rows = []
            for jm in LADDER:
                cfg = lib.make_config(jm, lib.reference_H(jm, ref), beta=beta)
                c12 = lib.d_max(cfg["s1"], cfg["mid"])
                c23 = lib.d_max(cfg["mid"], cfg["s3"])
                c13 = lib.d_max(cfg["s1"], cfg["s3"])
                rows.append({"c12": c12, "c23": c23, "c13": c13,
                             "deficit": c12 + c23 - c13})
            out[(ref, beta)] = rows
    return out


def test_lib_selftest_vs_exact_arithmetic():
    st = lib.selftest()
    assert st["fold_ratio_max_dev_vs_3b_exact"] < 1e-12
    assert st["mu0_commute_1"] < 1e-12 and st["mu0_commute_32"] < 1e-12
    assert st["c22_commute_1"] < 1e-12 and st["c22_commute_32"] > 0.5
    assert st["c21_norm_1"] < 1e-12 and st["c21_norm_32"] > 0.5


def test_norm_closed_form_sqrt6():
    for jm in LADDER:
        n = np.linalg.norm(lib.reference_H(jm, "null"), 2)
        pred = np.sqrt(6.0) * float(2 * jm) / float(2 * jm + 2)
        assert abs(n - pred) < 1e-12, jm


def test_b_parity_staircase_pair_freeze(legs):
    c = [r["c12"] for r in legs[("null", 1.0)]]
    assert abs(c[0] - c[1]) < 1e-12          # (1, 3/2) frozen
    assert abs(c[2] - c[3]) < 1e-12          # (2, 5/2) frozen
    assert c[2] - c[1] > 0.2                 # real step at 3/2 -> 2
    assert c[4] - c[3] > 0.1                 # real step at 5/2 -> 3
    # steps shrink (the staircase is converging-type, ratio ~ 0.53)
    assert (c[4] - c[3]) < 0.7 * (c[2] - c[1])


def test_flow_translation_identity_every_window(legs):
    for key, rows in legs.items():
        for r in rows:
            assert abs(r["c12"] - r["c23"]) < 1e-11, key


def test_rate_exclusions_gamma_and_thermal(legs):
    """Mixed-ref deficit: strictly shrinking diffs with ratios in
    (0.55, 0.85) -- excludes gamma (0.93-0.95) and thermal (0.368);
    ratios INCREASING (the power-law signature, not geometric-constant)."""
    for beta in (1.0, 0.5):
        d = [r["deficit"] for r in legs[("mixed", beta)]]
        diffs = [d[i + 1] - d[i] for i in range(4)]
        assert all(x > 0 for x in diffs), beta
        ratios = [diffs[i + 1] / diffs[i] for i in range(3)]
        assert all(r1 < r2 for r1, r2 in zip(ratios, ratios[1:])), beta
        for r in ratios:
            assert 0.55 < r < 0.85, (beta, ratios)


def test_interval_layer_stable():
    """Operational recovery + additivity at machine precision at the
    smallest and largest windows (null ref: even-weight sector, period pi)."""
    for jm in (LADDER[0], LADDER[-1]):
        rho, _ = lib.kms_state(jm)
        from scipy.linalg import expm
        om0 = lib.conj(expm(0.3j * lib.reference_H(jm, "null")), rho)
        period = np.pi                       # null ref preserves weight parity
        assert lib.trace_dist(lib.conj(lib.flow_U(jm, period), om0), om0) < 1e-12

        def ell(om_a, om_b, n_grid=600):
            taus = np.linspace(0.0, period, n_grid, endpoint=False)
            dv = [lib.d_max(lib.conj(lib.flow_U(jm, t), om_a), om_b)
                  for t in taus]
            i0 = int(np.argmin(dv))
            return taus[i0]

        dt = 0.35 * np.pi                    # on-grid for n_grid = 600
        om_a = lib.conj(lib.flow_U(jm, 0.1 * np.pi), om0)
        om_b = lib.conj(lib.flow_U(jm, 0.1 * np.pi + dt), om0)
        assert abs(ell(om_a, om_b) - dt) < 1e-10, jm
        om_c = lib.conj(lib.flow_U(jm, 0.1 * np.pi + dt + 0.2 * np.pi), om0)
        assert abs(ell(om_a, om_b) + ell(om_b, om_c)
                   - ell(om_a, om_c)) < 1e-10, jm


def test_admissibility_table_per_window():
    for jm in (LADDER[0], LADDER[1], LADDER[-1]):
        _, _, w = lib.wedge(jm)
        K = np.diag(w)

        def comm(name):
            G = lib.class_gen_folded(name, jm)
            return float(np.linalg.norm(K @ G - G @ K, 2))

        assert comm("(1.0,0.0) spacelike") < 1e-12, jm
        assert comm("(2.0,0.0) spacelike") < 1e-12, jm
        c22 = comm("(2.0,2.0) timelike")
        n21 = float(np.linalg.norm(
            lib.class_gen_folded("(2.0,1.0) spacelike", jm), 2))
        if float(jm) == 1.0:
            assert c22 < 1e-12 and n21 < 1e-12
        else:
            assert c22 > 0.5 and n21 > 0.5


def test_penalties_qualified_statement():
    """(1,0) commuting class: machine-zero excess at every window; the
    (1.5,1.5) non-commuting class: signed E(0.2) diffs NOT monotonically
    shrinking across the ladder (the robust non-stabilization statement)."""
    E_10, E_15 = [], []
    for jm in LADDER:
        cfg = lib.make_config(jm, lib.reference_H(jm, "null"))
        base = lib.deficit(cfg, cfg["mid"])
        for name, store in (("(1.0,0.0) spacelike", E_10),
                            ("(1.5,1.5) timelike", E_15)):
            G = lib.class_gen_folded(name, jm)
            store.append(lib.deficit(cfg, lib.kicked(cfg, G, 0.2)) - base)
    assert max(abs(e) for e in E_10) < 1e-10
    d = [abs(E_15[i + 1] - E_15[i]) for i in range(4)]
    assert not all(d[i + 1] < d[i] for i in range(3))
