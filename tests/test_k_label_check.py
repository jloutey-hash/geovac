"""Regression for the deterministic K-label hard-prohibition screen
(the /qa C5 backstop, debug/qa/check_k_label.py; run-#4 lesson, 2026-06-15).

Two guarantees:
  1. the live script exits 0 on the current corpus -> the gated TRUNK scope
     carries no SUSPECT K-tier assertion (a reintroduced trunk violation flips
     this to exit 1 and fails the test);
  2. the classifier flags a planted "K is now derived as a theorem" as SUSPECT
     and clears the compliant patterns (negation / hypothetical / future-goal /
     non-selection) -- so the screen genuinely discriminates.
"""
import importlib.util
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
CHECK = ROOT / "debug" / "qa" / "check_k_label.py"


def _load():
    spec = importlib.util.spec_from_file_location("check_k_label", CHECK)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _verdict(mod, text: str) -> str:
    with tempfile.TemporaryDirectory() as d:
        f = pathlib.Path(d) / "t.tex"
        f.write_text(text, encoding="utf-8")
        hits = mod.scan_file(f)
    got = [v for _, _, v, _ in hits]
    if "SUSPECT" in got:
        return "SUSPECT"
    return "COMPLIANT" if got else "NONE"


def test_trunk_gate_exit_zero():
    """The trunk scope is clean -> the screen exits 0 (audit hits are advisory)."""
    r = subprocess.run([sys.executable, str(CHECK)],
                       capture_output=True, text=True, encoding="utf-8")
    assert r.returncode == 0, (
        "a SUSPECT K-tier assertion entered the gated trunk scope:\n" + r.stdout)


def test_classifier_discriminates():
    mod = _load()
    # the run-#4 violation -> SUSPECT
    assert _verdict(mod,
        r"the Paper 2 combination rule is now \emph{derived} as a theorem "
        r"(no longer a mere observation).") == "SUSPECT"
    # a stale conjecture label -> SUSPECT
    assert _verdict(mod,
        r"Paper 2's conjectural combination rule has three components.") == "SUSPECT"
    # compliant patterns must NOT be SUSPECT
    for compliant in [
        r"It does not, however, derive the combination rule $K = \pi(B + F - \Delta)$.",
        r"Any future derivation of the combination rule must rely on a path integral.",
        r"No single-mechanism generation: no morphism derives the combination rule.",
        r"\item \textit{Derive the combination rule} from a spectral determinant.",
    ]:
        assert _verdict(mod, compliant) != "SUSPECT", compliant
