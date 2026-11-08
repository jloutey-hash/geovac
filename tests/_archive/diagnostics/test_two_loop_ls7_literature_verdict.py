"""ARCHIVED 2026-11-07 (diagnostics): LS-7 literature-target and verdict-tier tests.

Paper 36 no longer compares the two-loop model sum with a literature coefficient: the
C_2S = 3.63 target had no primary source and was removed (Eides-Grotch-Shelyuto Table II
gives +0.366 MHz two-loop SE and +0.127 MHz with two-loop VP on 2S). These tests check
module constants against their own literals and the verdict() tiering; they back no live
paper claim. Kept as institutional memory (CLAUDE.md sec. 14).
"""
from __future__ import annotations

import mpmath

from geovac.two_loop_self_energy import C_2S_LITERATURE, LS7_REFERENCE_2S_MHZ, verdict


def test_C_2S_literature_target():
    """Literature C_2S = +3.6266 from Eides 2001 Tab. 7.3."""
    assert abs(float(C_2S_LITERATURE) - 3.6266) < 1e-3


def test_LS7_reference_2S_MHz():
    """Eides 2001 Tab. 7.3 two-loop SE 2S = +0.857 MHz."""
    assert abs(float(LS7_REFERENCE_2S_MHZ) - 0.857) < 1e-3


def test_verdict_positive_tier():
    v = verdict(C_2S_LITERATURE * mpmath.mpf("1.05"))
    assert v["sign_correct"]
    assert "POSITIVE" in v["verdict_tier"]


def test_verdict_mixed_tier():
    v = verdict(C_2S_LITERATURE * mpmath.mpf("1.30"))
    assert v["sign_correct"]
    assert "WEAK: structural form correct" in v["verdict_tier"]  # renamed from MIXED, 1st-cert 2026-07-03


def test_verdict_weak_tier():
    """C_2S off by factor of ~5 lands in WEAK tier (within order of magnitude)."""
    v = verdict(C_2S_LITERATURE * mpmath.mpf("5.0"))
    assert "WEAK" in v["verdict_tier"]


def test_verdict_negative_tier():
    """C_2S off by 100x lands in NEGATIVE tier."""
    v = verdict(C_2S_LITERATURE * mpmath.mpf("100"))
    assert "NEGATIVE" in v["verdict_tier"]


def test_verdict_sign_check():
    """Verdict catches sign errors."""
    v = verdict(-C_2S_LITERATURE)
    assert not v["sign_correct"]
