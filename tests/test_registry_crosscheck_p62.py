"""Registry cross-check for Paper 62's qubit row (kept out of the paper-backing test file).

`tests/test_paper62_dzp_headline.py` recomputes the Slater DZP and cc-pVDZ Jordan-Wigner counts and 1-norms and
pins them as literals.  This file checks that the numeric registry entries the paper's \\gvq macros cite carry the
same literals.  Wrong answer rejected: a registry edit that moves the printed count or 1-norm without the guard
moving with it.  It lives here, not in the paper test, because the registry sits under debug/qa/ and paper-backing
tests must not depend on debug/ (C22 check D).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

from test_paper62_dzp_headline import DZP_PAULI_FC, DZP_LAMBDA_FC, CCPVDZ_PAULI_FC  # noqa: E402

sys.path.insert(0, os.path.join(HERE, "..", "debug", "qa"))
import numeric_registry as NR  # noqa: E402


def test_registry_carries_the_dzp_pauli_count_paper62_pins():
    assert NR.MEASURED["p62_dzp_pauli"]["value"] == DZP_PAULI_FC
    assert NR.MEASURED["p62_dzp_lambda_fc"]["value"] == DZP_LAMBDA_FC
    assert NR.CITED["h2o_ccpvdz_pauli"]["value"] == CCPVDZ_PAULI_FC
