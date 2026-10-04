"""Regression oracle for the promoted LiH / prolate stack (Paper 19).

`tests/data/lih_golden/golden.json` holds full-precision outputs captured on
2026-10-04 from the pre-promotion debug/ code (tag v5.20.2).  Each item here
re-runs one capture in a FRESH, single-threaded process through the promoted
`geovac.prolate_allelectron`, `geovac.lih_r12ci` and `geovac.lih_vmc` code and
requires a BIT-IDENTICAL match (`capture.py --compare-golden`, exit 0).

Fresh processes are required: several modules set mpmath precision or flip
`geovac.lih_r12ci.kernels.USE_EXACT_NEUMANN` at import, and module-level RNG
streams feed randomized eigh/SVD, so results depend on import/call order.  BLAS
threads are pinned to 1 because multi-threaded summation changes low bits.

Fast items run by default; slow items (minutes to ~3 h, up to ~8.8 GB, e.g. the
-8.049045089 Ha pi marriage) run under ``--slow``.

Rejected wrong answers: any change to a value, array byte or shape in the
promoted code (compare-golden fails on a 1-ulp change; checked 2026-10-04).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
GOLD_DIR = ROOT / "tests" / "data" / "lih_golden"
CAPTURE = GOLD_DIR / "capture.py"
_GOLDEN = json.loads((GOLD_DIR / "golden.json").read_text(encoding="utf-8"))


def _items() -> list:
    out = []
    for name, rec in _GOLDEN["items"].items() if "items" in _GOLDEN else _GOLDEN.items():
        if not isinstance(rec, dict) or "tier" not in rec:
            continue
        marks = [pytest.mark.slow] if rec["tier"] != "fast" else []
        out.append(pytest.param(name, marks=marks, id=name))
    return out


def _env() -> dict:
    env = dict(os.environ)
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[k] = "1"
    return env


@pytest.mark.parametrize("name", _items())
def test_promoted_code_reproduces_golden(name: str, tmp_path: Path) -> None:
    run = subprocess.run([sys.executable, str(CAPTURE), "--run-item", name, "--stage", str(tmp_path)],
                         cwd=ROOT, env=_env(), capture_output=True, text=True)
    assert run.returncode == 0, f"{name} capture failed:\n{run.stdout[-2000:]}\n{run.stderr[-2000:]}"
    cmp = subprocess.run([sys.executable, str(CAPTURE), "--compare-golden", name, "--stage", str(tmp_path)],
                         cwd=ROOT, env=_env(), capture_output=True, text=True)
    assert cmp.returncode == 0, f"{name} differs from golden.json:\n{cmp.stdout[-3000:]}"
