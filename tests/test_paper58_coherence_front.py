"""Backing test for Paper 58 Sec. "The continuous side" (the many-electron
front) -- the coherence-front paragraph, relocated here from Paper 60 on
2026-09-06 (it is bonding physics, not encoding cost).

Claim: correlation moves the SIGNED cross-center coherence front, and it does
so through the occupation numbers of the antibonding natural orbital, whose
per-center weight is amplified by (1+S)/(1-S) relative to the bonding one.
Two-orbital form:

    M2 = sum_k n_k w_k^2 cos(theta_k) / sum_k n_k w_k^2,
    w_g^2 = 1/(2(1+S)), w_u^2 = 1/(2(1-S)), cos_g = +S, cos_u = -S.

Limits (SYMBOLIC): molecular-orbital (n_u = 0) gives M2 = S; Heitler-London
occupations n_g = (1+S)^2/(1+S^2), n_u = (1-S)^2/(1+S^2) give M2 = S^2.

Rejected wrong answer: "the antibonding orbital enters with the SAME weight as
the bonding one" (w_u^2 = w_g^2).  Under that plant, the measured H2 point
(n_g = 1.969, n_u = 0.025, S = 0.674) would give 0.657, not the 0.591 the
formula gives and the 0.596 the FCI density measures.
"""

from __future__ import annotations

import pytest

sp = pytest.importorskip("sympy")


def _m2(n_g, n_u, S, w_g2, w_u2):
    return (n_g * w_g2 * S - n_u * w_u2 * S) / (n_g * w_g2 + n_u * w_u2)


def test_m2_two_orbital_limits_symbolic():
    S = sp.symbols("S", positive=True)
    w_g2 = 1 / (2 * (1 + S))
    w_u2 = 1 / (2 * (1 - S))
    # MO limit
    assert sp.simplify(_m2(2, 0, S, w_g2, w_u2) - S) == 0
    # Heitler-London limit
    n_g = (1 + S) ** 2 / (1 + S**2)
    n_u = (1 - S) ** 2 / (1 + S**2)
    assert sp.simplify(_m2(n_g, n_u, S, w_g2, w_u2) - S**2) == 0
    # the amplification factor itself
    assert sp.simplify(w_u2 / w_g2 - (1 + S) / (1 - S)) == 0


def test_m2_measured_point_needs_the_amplification():
    n_g, n_u, S = 1.969, 0.025, 0.674
    w_g2 = 1.0 / (2.0 * (1.0 + S))
    w_u2 = 1.0 / (2.0 * (1.0 - S))
    m2 = _m2(n_g, n_u, S, w_g2, w_u2)
    assert abs(m2 - 0.591) < 0.004, m2          # formula vs the paper's 0.591
    assert abs(m2 - 0.596) < 0.010, m2          # vs the FCI-density measurement
    # the plant: equal per-center weights -> no collapse toward S^2
    m2_plant = _m2(n_g, n_u, S, w_g2, w_g2)
    assert m2_plant > 0.65 and abs(m2_plant - m2) > 0.05, m2_plant
    # occupation-weighted UNSIGNED angle (M1) does not move at all
    m1 = (n_g * S + n_u * S) / (n_g + n_u)
    assert abs(m1 - S) < 1e-12


# ---------------------------------------------------------------------------
# Slow: drive the correlation ladder in its --quick mode (one R per molecule,
# no basis doubling) and check the pointwise residual claim of Paper 60:
# the occupied-space principal-angle cosine matches the 1s-1s law at the
# empirical zeta_eff to |Delta cos| <= 0.01 at both rungs, both molecules.
# The wrong answer this rejects: a front that moves under correlation by
# something other than the decay length (it would show up as a residual of
# order the M2 shift, ~0.05-0.1 in cos).
# ---------------------------------------------------------------------------

def _load_ladder():
    import importlib.util
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "debug" / "decompactification_correlation_ladder.py"
    if not path.exists():
        pytest.skip("exploratory ladder driver pruned; claim chronicled in CHANGELOG v5.10.2")
    spec = importlib.util.spec_from_file_location("decomp_ladder", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.slow
def test_ladder_quick_h2_pointwise_residual():
    mod = _load_ladder()
    out = mod.run(quick=True, do_double=False)
    checked = 0
    for mol in ("H2|base", "HeH+|base"):
        for rec in out["ladder"][mol]["records"]:
            if "error" in rec:
                continue
            for rung in ("HF", "FCI"):
                a = rec[rung]
                resid = a["M1"] - a["cos_pred_B"]
                assert abs(resid) <= 0.01, (mol, rung, rec["R"], resid)
                # the signed coherence never exceeds the unsigned angle measure
                assert a["M2"] <= a["M1"] + 1e-9, (mol, rung, a["M2"], a["M1"])
                checked += 1
    assert checked >= 4
