"""Mirror test for the C17 headline-number registry gate (2026-07-02).

Two legs, matching the other deterministic-gate mirrors:
  1. the gate PASSES on the current corpus (exit 0), and
  2. every registry family demonstrably CATCHES a synthetic wrong value
     (and does not fire on the canonical/legitimate forms) -- the gate's
     own sensitivity/specificity calibration.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "debug" / "qa" / "check_headline_numbers.py"

sys.path.insert(0, str(SCRIPT.parent))
import check_headline_numbers as chn  # noqa: E402


def _entry(eid: str) -> dict:
    return next(e for e in chn.REGISTRY if e["id"] == eid)


def test_gate_passes_on_current_corpus() -> None:
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--gate", "group4"],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_z_span_family_catches_wrong_value_in_context() -> None:
    e = _entry("library-z-span")
    live, _ = chn.scan_entry(e, text_override=(
        "the composed molecular library (35 molecules spanning\n"
        "$Z = 1$--36, including the first transition series)"))
    assert len(live) == 1
    # canonical value: clean
    live, _ = chn.scan_entry(e, text_override=(
        "library spanning\n$Z = 1$--56 (H through Ba)"))
    assert not live
    # different quantity (periodic-row prose, no library context): clean
    live, _ = chn.scan_entry(e, text_override=(
        "First-row ($Z=1$--10) and second-row atoms are fully\nsupported."))
    assert not live


def test_floor_family_catches_non_51_floor() -> None:
    e = _entry("pauli-advantage-floor")
    live, _ = chn.scan_entry(
        e, text_override=r"advantage reported here (190$\times$--1{,}712$\times$)")
    assert len(live) == 1
    live, _ = chn.scan_entry(
        e, text_override=r"advantage reported here (51$\times$--1{,}712$\times$)")
    assert not live


def test_library_count_family_catches_retired_counts() -> None:
    e = _entry("library-count")
    live, _ = chn.scan_entry(e, text_override="a library of 38 systems across rows")
    assert len(live) == 1
    live, _ = chn.scan_entry(e, text_override="a library of 37 systems across rows")
    assert not live


def test_binds_family_catches_undisclosed_3015() -> None:
    e = _entry("balanced-lih-binds-at-3015")
    live, _ = chn.scan_entry(e, text_override="balanced coupled binds LiH at R_eq=3.015 bohr")
    assert len(live) == 1
    live, _ = chn.scan_entry(e, text_override=(
        "binds LiH at a computed R_eq=3.227 bohr (7.0\\% above the\n"
        "experimental 3.015 bohr)"))
    assert not live


def test_onenorm_family_catches_stale_pair_not_pauli_ratio() -> None:
    e = _entry("lih-onenorm-stale")
    live, _ = chn.scan_entry(e, text_override=r"essentially identical 1-norm ($0.97\times$).")
    assert len(live) == 1
    live, _ = chn.scan_entry(e, text_override=r"$\lambda = 33.3$~Ha, and 21 QWC groups")
    assert len(live) == 1
    # the l-parity Pauli-ratio 0.97 cells are a DIFFERENT quantity: clean
    live, _ = chn.scan_entry(
        e, text_override=r"The Pauli-count ratio is uniformly $\sim 0.97$ (a 3\% reduction)")
    assert not live
