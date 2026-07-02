"""Paper 14 MPO bond rank (Thm 3.2.A.unified) — durable backing.

The theorem's empirical companions (LiH statement-A exactness 29/29, the
universal interior chi_H profile [4,16,16,9,9,9,6,3,3,2], the composed
chi=2 sub-block boundaries, and the balanced {9,16} = 2 + 7*N_cross
quantization across LiH/BeH2/H2O) were backed only by transient debug/
drivers until 2026-07-01 (6th group4 cert follow-on, PI-directed). The
drivers now live permanently in tests/mpo_rank_support/ (wh7_support
precedent); this file (a) RECOMPUTES the full LiH panel from the live
composed builder via the driver's own functions and (b) pins the shipped
panel JSONs so they cannot silently drift.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_SUPPORT = Path(__file__).resolve().parent / "mpo_rank_support"
sys.path.insert(0, str(_SUPPORT))

import sprint_s2_v2_unified_panel as panel  # noqa: E402

# Paper 14 §MPO bond rank: the universal interior CORE (local cuts 2..8 of a
# 10-qubit sub-block) is basis-determined; the edge cuts are position-
# dependent (first cut 4 for the leading sub-block / 5 with a completed block
# to the left; chain-end saturates to 2). Full per-window LiH profiles below.
# (The pre-2026-07-01 paper text claimed the FULL 10-cut profile was
# sub-block-role independent; the shipped data + live recompute show the
# first-cut 4-vs-5 split — the paper was corrected to the core+edges form.)
EXPECTED_WINDOWS = {
    1: [4, 16, 16, 9, 9, 9, 6, 3, 3, 2],   # leading sub-block, cuts 1..10
    11: [5, 16, 16, 9, 9, 9, 6, 3, 3, 2],  # middle sub-block, cuts 11..20
    21: [5, 16, 16, 9, 9, 9, 6, 3, 2],     # final sub-block, cuts 21..29
}
UNIVERSAL_CORE = [16, 16, 9, 9, 9, 6, 3]   # local cuts 2..8, every window
BOUNDARY_CHI = 2


def _load(name: str) -> dict:
    with open(_SUPPORT / "data" / name) as f:
        return json.load(f)


def test_unified_panel_json_pins_theorem_counts() -> None:
    """The shipped LiH panel carries the paper's 29/29 statement counts and
    the universal interior profile in every sub-block window."""
    d = _load("sprint_s2_v2_unified_panel.json")
    v = d["verification"]
    assert v["eq_h1_count"] == "29/29"      # (A) exact closed form
    assert v["ineq_sub_count"] == "29/29"   # (B) subadditivity
    assert v["ineq_perL_count"] == "29/29"  # (C) per-L subadditivity
    rows = {r["cut"]: r for r in d["rows"]}
    assert len(rows) == 29
    for lo, expected in EXPECTED_WINDOWS.items():
        prof = [rows[c]["chi_H"] for c in range(lo, lo + len(expected))]
        assert prof == expected, (lo, prof)
        assert prof[1:8] == UNIVERSAL_CORE, (lo, prof)
    for boundary in (10, 20):
        assert rows[boundary]["chi_H"] == BOUNDARY_CHI, boundary


def test_balanced_panel_json_pins_two_value_quantization() -> None:
    """Balanced-coupled negative control: 12/12 sub-block boundaries across
    LiH/BeH2/H2O break the composed chi=2 floor into chi in {9,16} with
    lift = 7 * N_cross, N_cross in {1,2}."""
    b = _load("sprint_s2_v2_balanced_library_panel.json")
    assert set(b) == {"LiH", "BeH2", "H2O"}
    n_boundaries = 0
    for mol, rec in b.items():
        for bd in rec["boundary_data"]:
            n_boundaries += 1
            assert bd["chi_composed"] == BOUNDARY_CHI, (mol, bd)
            assert bd["chi_balanced"] in (9, 16), (mol, bd)
            lift = bd["chi_balanced"] - bd["chi_composed"]
            assert lift % 7 == 0 and lift // 7 in (1, 2), (mol, bd)
    assert n_boundaries == 12, n_boundaries


@pytest.mark.slow
def test_unified_panel_recomputes_from_live_builder() -> None:
    """Statement (A) exactness + (B)/(C) subadditivity recomputed end-to-end
    from the live composed builder via the migrated driver's own functions,
    and cross-checked row-by-row against the shipped JSON (drift guard)."""
    from geovac.composed_qubit import build_composed_hamiltonian, _enumerate_states
    from geovac.molecular_spec import lih_spec

    spec = lih_spec()
    res = build_composed_hamiltonian(spec, verbose=False)
    h1, eri, H_qop = res["h1"], res["eri"], res["qubit_op"]
    Q = 2 * h1.shape[0]

    states_by_block, block_offsets, offset = [], [], 0
    for blk in spec.blocks:
        cs = _enumerate_states(blk.max_n, l_min=getattr(blk, "l_min", 0))
        states_by_block.append(cs)
        block_offsets.append(offset)
        offset += len(cs)
        if blk.has_h_partner:
            pn = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            ps = _enumerate_states(pn)
            states_by_block.append(ps)
            block_offsets.append(offset)
            offset += len(ps)

    h1_qop = panel.h1_to_qubit(h1)
    Vee_qop = panel.build_Vee_total_qubit(eri)
    Vee_L = panel.build_Vee_per_L_qubit(eri, states_by_block, block_offsets)
    shipped = {r["cut"]: r for r in _load("sprint_s2_v2_unified_panel.json")["rows"]}

    for cut in range(1, Q):
        chi_H = panel.operator_schmidt_rank(H_qop, cut, Q)
        chi_h1 = panel.operator_schmidt_rank(h1_qop, cut, Q)
        chi_Vee = panel.operator_schmidt_rank(Vee_qop, cut, Q)
        rk, LL, RR = panel.h1_cross_cut_rank(h1, cut)
        assert chi_h1 == 2 * rk + LL + RR, f"(A) fails at cut {cut}"
        assert chi_H <= chi_h1 + chi_Vee, f"(B) fails at cut {cut}"
        sum_L = sum(panel.operator_schmidt_rank(Vee_L[L], cut, Q) for L in Vee_L)
        assert chi_Vee <= sum_L, f"(C) fails at cut {cut}"
        s = shipped[cut]
        assert (chi_H, chi_h1, chi_Vee) == (s["chi_H"], s["chi_h1"], s["chi_Vee"]), (
            f"panel drift at cut {cut}: recomputed "
            f"{(chi_H, chi_h1, chi_Vee)} vs shipped "
            f"{(s['chi_H'], s['chi_h1'], s['chi_Vee'])}"
        )
