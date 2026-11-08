"""KNOWN-GAP 1/5, step 4 -- prove the six new C16 entries DISCRIMINATE.

REGISTRY DISCRIMINATION RULE (qa.md, hard): whenever you add or edit a C16
entry you must prove it fires on the RETIRED wording and stays SILENT on the
CORRECTED wording, and report both.  "An entry that fires on nothing is worse
than no entry at all, because the gate reports PASS and the class is now
believed guarded."  This has bitten this corpus twice, both times found only
because someone tested.

Runs each entry's compiled pattern against synthetic probes -- retired wording
that MUST match, corrected wording that MUST NOT -- rather than against the live
corpus, so it keeps working after the corpus is clean.
"""
from __future__ import annotations

import re
import sys

sys.path.insert(0, "debug/qa")
from check_retracted_terms import REGISTRY  # noqa: E402

PROBES = {
    "p60-removability-corollary": {
        "fire": [
            "a truncation-side price is a property of the matrix, which a "
            "preconditioner reaches",
            "a continuum-side price is a property of the symbol, which no "
            "congruence of the finite section can touch",
        ],
        "silent": [
            "a banded congruence multiplies the symbol by a trigonometric "
            "polynomial, which can cancel a zero of finite order but cannot "
            "alter a decay or smoothness class",
            "the preconditioner is built from the symbol: its matching "
            "polynomial is chosen precisely to share the symbol's zero",
        ],
    },
    "p60-frames-completeness": {
        "fire": [
            "overcompleteness is the price of completeness, not a defect of "
            "the basis",
            "completeness of the one-centre set alone forces lam_min to zero",
        ],
        "silent": [
            "the Bessel deficit plateaus at 0.380, 0.696 and 0.907, so the "
            "one-centre set is measurably far from complete in the molecular "
            "metric",
            "Ron-Shen is the surviving mechanism: the near-dependence is one "
            "direction",
        ],
    },
    "p60-sigma-law-derived": {
        "fire": [
            "The growth law is derived (a band-limited concentration rate)",
            "the conditioning grows polynomially, with a derived band-limited law",
            "Backs the derived conditioning law (CHANGELOG v4.103.0)",
        ],
        "silent": [
            "the growth law is not new: it is the Kac-Murdock-Szego asymptotic, "
            "and what is ours is the identification",
            "PRIOR ART: this asymptotic is NOT derived here",
        ],
    },
    "p60-prop-d-as-new": {
        "fire": [
            "Proposition D: a block-diagonal congruence cannot orthogonalize a "
            "metric that is not block diagonal",
            "our Proposition establishes that a block-diagonal congruence "
            "preserves the grading",
        ],
        "silent": [
            "this is Loewdin symmetry preservation specialised to the l "
            "grading, known since Slater-Koster; what the paper claims is the "
            "l-versus-m application",
        ],
    },
    "p60-translation-ours": {
        "fire": [
            "with symbol j0(kR cot(chi/2)); equivalently that the "
            "Shibuya--Wulfman\noperator is multiplication by the translation "
            "phase on the Fock sphere",
            "the translation identification is ours",
        ],
        "silent": [
            "that is prior art on three counts, and earlier versions of this "
            "paper claimed it; what is ours is the symbol",
        ],
    },
    "p60-tau-membership": {
        "fire": [
            "for matrices in the tau (DST-I)\nalgebra --- which is to say "
            "Toeplitz minus Hankel, the structure of this\nsection --- their "
            "result is an identity",
        ],
        "silent": [
            "our Toeplitz-minus-Hankel form does not put our matrices inside "
            "it: tau membership requires the coefficient sequence to terminate, "
            "which the chirp symbol's does not",
        ],
    },
}


def main() -> int:
    by_id = {e["id"]: e for e in REGISTRY}
    ok = True
    for eid, probe in PROBES.items():
        if eid not in by_id:
            print(f"  [MISSING] {eid} not in REGISTRY")
            ok = False
            continue
        rx = re.compile(by_id[eid]["pattern"])
        fired = [bool(rx.search(s)) for s in probe["fire"]]
        quiet = [not rx.search(s) for s in probe["silent"]]
        good = all(fired) and all(quiet)
        ok = ok and good
        print(f"  [{'ok' if good else 'FAIL'}] {eid:32s} "
              f"fires {sum(fired)}/{len(fired)}   silent {sum(quiet)}/{len(quiet)}")
        for s, f in zip(probe["fire"], fired):
            if not f:
                print(f"        DID NOT FIRE on: {s[:78]!r}")
        for s, q in zip(probe["silent"], quiet):
            if not q:
                print(f"        FIRED ON CORRECTED TEXT: {s[:78]!r}")
    print("\nRESULT:", "PASS -- all six discriminate both ways" if ok
          else "FAIL -- an entry does not discriminate")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
