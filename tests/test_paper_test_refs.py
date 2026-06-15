"""Regression for the deterministic paper<->test reference integrity check
(the /qa C13 criterion, debug/qa/check_paper_test_refs.py).

Two guarantees:
  1. the trunk-gate exits 0 -> every inline `tests/test_*` reference in the
     trunk papers resolves to a live test (a renamed/typo'd/archived ref flips
     this to exit 1 and fails the test);
  2. the existence classifier discriminates -- a real test is LIVE, a glob ref
     resolves if any file matches, a fabricated name is MISSING.
"""
import importlib.util
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CHECK = ROOT / "debug" / "qa" / "check_paper_test_refs.py"


def _load():
    spec = importlib.util.spec_from_file_location("check_paper_test_refs", CHECK)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_trunk_gate_exit_zero():
    r = subprocess.run([sys.executable, str(CHECK)],
                       capture_output=True, text=True, encoding="utf-8")
    assert r.returncode == 0, (
        "a trunk paper cites a non-existent test:\n" + r.stdout)


def test_status_discriminates():
    mod = _load()
    # a real trunk-backing test
    assert mod.test_status("test_fock_projection") == "LIVE"
    # a glob family (>=1 file matches)
    assert mod.test_status("test_qed_*") == "LIVE"
    # a fabricated name
    assert mod.test_status("test_this_does_not_exist_zzz") == "MISSING"


def test_stem_normalizes_latex():
    mod = _load()
    assert mod.stem(r"tests/test\_fock\_projection.py") == "test_fock_projection"
    assert mod.stem(r"tests/test\_qed\_*.py") == "test_qed_*"
