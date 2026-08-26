"""Shared plumbing for the certified reference-value generators.

Certification discipline (see docs/certified_reference_values.md):

  * ``value``          -- the number, printed to AT MOST ``digits_claimed``
                          significant digits.
  * ``digits_claimed`` -- how many of those digits are certified.  Never larger
                          than what the evidence field actually establishes.
  * ``method``         -- how the number was produced.
  * ``evidence``       -- what certifies ``digits_claimed``, with the measured
                          agreements that back it.

Nothing in this package invents a digit count: every count is either exact
(rational / algebraic), or the measured agreement of two routes, or a value
quoted verbatim from a cross-validated corpus run whose accounting is cited.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any, Dict, List

import mpmath as mp

VALUE_KINDS = {"decimal", "rational", "algebraic", "complex_decimal",
               "integer_tuple", "expression"}


def entry(entry_id: str, category: str, label: str, value: str,
          digits_claimed: Any, method: str, evidence: str,
          value_kind: str = "decimal", **extra: Any) -> Dict[str, Any]:
    """Build one table row.  ``digits_claimed`` is an int or the string 'exact'."""
    assert value_kind in VALUE_KINDS, f"bad value_kind {value_kind}"
    row: Dict[str, Any] = {
        "id": entry_id,
        "category": category,
        "label": label,
        "value": value,
        "value_kind": value_kind,
        "digits_claimed": digits_claimed,
        "method": method,
        "evidence": evidence,
    }
    row.update(extra)
    return row


def nstr(x, digits: int) -> str:
    """Print an mpf to exactly ``digits`` significant digits, no exponent games."""
    return mp.nstr(mp.mpf(x), digits, strip_zeros=False)


def agree_digits(a, b) -> int:
    """Number of leading significant decimal digits on which a and b agree.

    Returns a conservative integer: floor(-log10(|a-b| / |a|)).  Returns a large
    sentinel (10**6) when the two are bit-identical.
    """
    a, b = mp.mpf(a), mp.mpf(b)
    if a == b:
        return 10 ** 6
    scale = abs(a) if a != 0 else mp.mpf(1)
    return int(mp.floor(-mp.log10(abs(a - b) / scale)))


def frac_str(f: Fraction) -> str:
    return f"{f.numerator}/{f.denominator}" if f.denominator != 1 else str(f.numerator)


def sort_entries(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    order = {"T2": 0, "two_center_eri": 1, "slater_rational": 2,
             "qfd_diatomic": 3, "helium_ci": 4, "resurgent": 5, "anchor": 6}
    return sorted(rows, key=lambda r: (order.get(r["category"], 99),
                                       r.get("sort_index", 0), r["id"]))
