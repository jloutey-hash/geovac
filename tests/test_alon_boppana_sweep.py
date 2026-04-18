"""Smoke test for the Alon-Boppana sweep JSON (Track RH-D, Sprint 2).

Verifies that `debug/data/alon_boppana_sweep.json` contains the expected
graph families, that the larger-size rows we added in Sprint 2 are
present, and that the recorded numerical invariants are internally
consistent (deviation = max|μ_non-trivial| − √q_max).

This is a smoke test only.  The real Ihara-zeta / Ramanujan logic is
exercised in `tests/test_ihara_zeta.py` and
`tests/test_ihara_zeta_dirac.py`.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


DATA = Path(__file__).parent.parent / "debug" / "data" / "alon_boppana_sweep.json"


def _load():
    if not DATA.exists():
        import pytest
        pytest.skip(
            f"{DATA} not present; run debug/compute_alon_boppana_sweep.py first"
        )
    with DATA.open() as f:
        return json.load(f)


def test_sweep_schema_and_expected_rows():
    """The JSON loads, has the expected graph families, and records
    all Sprint 2 larger-size computations that were requested.

    Expected Sprint 2 rows (at minimum):
      - S3_Coulomb   at n_label = 4 and 5.
      - S5_Bargmann  at n_label = 4 and 5.
      - Dirac_A      at n_label = 4.
      - Dirac_B      at n_label = 4.
    """
    data = _load()
    assert "sprint" in data, "missing 'sprint' key"
    assert data["sprint"] == "RH-D", f"expected RH-D, got {data['sprint']}"
    assert "rows" in data and len(data["rows"]) > 0, "missing rows"
    assert "fits" in data, "missing fits"

    # Index rows by (family, n_label)
    idx = {(r["label"], r["n_label"]): r for r in data["rows"]}

    expected = [
        ("S3_Coulomb", 4),
        ("S3_Coulomb", 5),
        ("S5_Bargmann", 4),
        ("S5_Bargmann", 5),
        ("Dirac_A", 4),
        ("Dirac_B", 4),
    ]
    for key in expected:
        assert key in idx, f"missing Sprint 2 row: {key}"

    # Internal consistency: deviation == max|mu_nontriv| - sqrt(q_max)
    # for every row where q_max > 0.
    for key, r in idx.items():
        if r.get("q_max", 0) <= 0:
            continue
        computed_dev = r["max_abs_nontrivial"] - r["sqrt_q_max"]
        recorded_dev = r["deviation"]
        assert math.isclose(computed_dev, recorded_dev, abs_tol=1e-6), (
            f"row {key}: deviation inconsistent: "
            f"max|mu|-sqrt(q_max) = {computed_dev:.6f} vs "
            f"recorded {recorded_dev:.6f}"
        )

    # Fits dictionary is populated for every family.
    for family in ("S3_Coulomb", "S5_Bargmann", "Dirac_A", "Dirac_B"):
        assert family in data["fits"], f"fits missing family {family}"
