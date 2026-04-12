"""Tests for Paper 14 revision: composed systems section."""
import os
import re
import pytest

PAPER_PATH = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    "paper_14_qubit_encoding_revised.tex",
)


def _read_paper() -> str:
    with open(PAPER_PATH, "r", encoding="utf-8") as f:
        return f.read()


def test_paper14_exists():
    """The revised .tex file exists."""
    assert os.path.isfile(PAPER_PATH), f"Revised paper not found at {PAPER_PATH}"


def test_paper14_contains_composed():
    """File contains 'Composed' or 'composed' (new section present)."""
    text = _read_paper()
    assert "omposed" in text, "Paper does not contain 'Composed' or 'composed'"


def test_paper14_contains_trenev():
    """File contains 'Trenev' (reference added)."""
    text = _read_paper()
    assert "Trenev" in text, "Paper does not contain 'Trenev' reference"


def test_paper14_contains_h2o():
    """File contains H2O data (as LaTeX H$_2$O)."""
    text = _read_paper()
    assert "H$_2$O" in text or "H₂O" in text or "H_2O" in text, (
        "Paper does not contain H2O data"
    )


def test_paper14_no_460():
    """File does NOT contain '4.60' without qualification.

    The old Gaussian exponent 4.60 is superseded by published data.
    If it appears, it should be in context like 'previously estimated'.
    """
    text = _read_paper()
    # Find all occurrences of '4.60'
    matches = [m.start() for m in re.finditer(r"4\.60", text)]
    for pos in matches:
        # Check surrounding context for qualification
        context = text[max(0, pos - 100) : pos + 100]
        qualified = any(
            word in context.lower()
            for word in ["earlier", "previously", "superseded", "estimated", "old"]
        )
        assert qualified, (
            f"Unqualified '4.60' found at position {pos}. "
            f"Context: ...{context}..."
        )
