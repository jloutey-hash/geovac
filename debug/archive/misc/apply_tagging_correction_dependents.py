"""Propagate the C23 run #3 correction to every document that restated it.

The claim was corrected in Paper 60 (its owner).  Per CLAUDE.md Sec. 9's
retraction->dependents rule, the documents whose ARGUMENT rests on it must move
too, or this becomes exactly the class that rule exists to stop: corrected in
the owner, left standing in the citers.

Dependents, enumerated from the argument rather than by grep:
  docs/claim_test_matrix.md  -- the row asserting the Bessel-free measurement is
                                an independent route and "the whole evidential
                                point"
  CHANGELOG.md v5.11.4       -- both the "no Bessel function of any kind" framing
                                and the false removability sentence
  debug/sprint_contraction_seam_memo.md Sec. 3 -- same two

Idempotent.
"""
from __future__ import annotations

import sys

EDITS = [
    # ---------------- claim matrix ----------------
    ("docs/claim_test_matrix.md",
     "BACKED-SOUND. The Bessel-free route is the "
     "whole evidential point and it is guarded against being vacuous: a "
     "companion guard raises the symbol's zero from quadratic to quartic and "
     "requires the constant to COLLAPSE, so the π^2 cannot be arriving by "
     "construction.",
     "BACKED-SOUND **as provenance** (corrected 2026-09-12, C23 run #3). The "
     "Bessel-free measurement is NOT an independent route: it is the `b == 1` "
     "case of the same Kac-Murdock-Szego theorem, so it is a change of "
     "REPRESENTATION. It is still guarded against being vacuous — a companion "
     "guard raises the symbol's zero from quadratic to quartic and requires the "
     "constant to COLLAPSE, so the π^2 cannot be arriving by construction — and "
     "it still shows the constant carries no Bessel content. Both the "
     "Dirichlet-eigenvalue reading and the extremal (Wirtinger-Sobolev) problem "
     "are Böttcher-Widom's, in the source the paper already cites; so is the "
     "independence of `c_alpha` from `b`."),

    ("docs/claim_test_matrix.md",
     "Second route for c_1: "
     "the exact tridiagonal spectrum, `tests/test_paper60_kms_attribution.py` "
     "(v5.10.18) — two representations, one constant.",
     "**The removability corollary this row originally carried is WITHDRAWN** "
     "(\"a truncation-side price is matrix-level and reachable; a "
     "continuum-side price is symbol-level and is not\"): the preconditioner is "
     "built FROM the symbol — its matching polynomial shares the symbol's zero "
     "— and preconditioning IS a congruence of the finite section. The "
     "surviving statement is about the symbol on both sides: a banded "
     "congruence cancels a finite-order zero but cannot alter a decay class."),

    # ---------------- CHANGELOG ----------------
    ("CHANGELOG.md",
     "- **`pi^2` of `eq:sigma_law` is TRUNCATION-side.** It is the "
     "Kac-Murdock-Szego `c_1`, the first Dirichlet eigenvalue of `-d^2/dx^2` on "
     "the unit interval, and the measurement above obtains it from a "
     "band-limited concentration problem **containing no Bessel function of any "
     "kind**. The price of the finite basis, not of the second centre. "
     "(`pi^2 . Q` half of M2. Second route: the exact tridiagonal spectrum, "
     "v5.10.18.)",
     "- **`pi^2` of `eq:sigma_law` comes from the symbol's zero order plus the "
     "finite section, and carries no Bessel content.** It is the "
     "Kac-Murdock-Szego `c_1`. **Corrected same day by C23 run #3:** the "
     "reading of `c_1` as a Dirichlet eigenvalue, the extremal "
     "(Wirtinger-Sobolev) problem behind it, and the independence of `c_alpha` "
     "from `b` are ALL Boettcher-Widom's — in the source this paper already "
     "cites, whose title names the inequality. And the Bessel-free measurement "
     "is the `b == 1` case of that same theorem, so it is a change of "
     "REPRESENTATION, not an independent route. (`pi^2 . Q` half of M2.)"),

    ("CHANGELOG.md",
     "Operational payoff, and it is the mechanism behind v5.11.0's split: "
     "**a truncation-side price is a property of the matrix, which a "
     "preconditioner reaches; a continuum-side price is a property of the "
     "symbol, which no congruence of the finite section can touch.**",
     "**The operational corollary first written here was WITHDRAWN the same day** "
     "(C23 run #3). It read \"a truncation-side price is a property of the "
     "matrix, which a preconditioner reaches; a continuum-side price is a "
     "property of the symbol, which no congruence of the finite section can "
     "touch\", and both halves are wrong: the v5.11.0 preconditioner is built "
     "FROM the symbol (its matching polynomial is chosen to share the symbol's "
     "zero, per Serra), and preconditioning IS a congruence, one that replaces "
     "the symbol by `f/g`. The surviving statement concerns the symbol on both "
     "sides and turns on the KIND of feature: a banded congruence multiplies the "
     "symbol by a trigonometric polynomial, which cancels a zero of finite order "
     "but cannot alter a decay class. The conditioning pole is a second-order "
     "zero; the locality pole is the chirp's `j^-5/4` envelope. Provenance of "
     "the constant predicts nothing here."),

    # ---------------- memo ----------------
    ("debug/sprint_contraction_seam_memo.md",
     "**The `pi^2` of `eq:sigma_law` is truncation-side.** It is the KMS constant\n"
     "`c_1`, the first Dirichlet eigenvalue of `-d^2/dx^2` on the unit interval, and\n"
     "F2 obtains it from a band-limited concentration problem **containing no Bessel\n"
     "function of any kind**. It is the price of the finite basis, not of the second\n"
     "centre. Paper 18 tier: calibration, M2 (`pi^2·Q` half).",
     "**The `pi^2` of `eq:sigma_law` carries no Bessel content.** It is the KMS\n"
     "constant `c_1`, fixed by the symbol's zero order together with the finite\n"
     "section. Paper 18 tier: calibration, M2 (`pi^2·Q` half). **Corrected by C23\n"
     "run #3 (same day):** the Dirichlet-eigenvalue reading, the extremal\n"
     "Wirtinger–Sobolev problem, and the independence of `c_alpha` from `b` are all\n"
     "Böttcher–Widom's, in the source already cited; and F2 is the `b == 1` case of\n"
     "that theorem, i.e. a change of representation, **not** an independent route.\n"
     "The independent-route claim below is withdrawn."),

    ("debug/sprint_contraction_seam_memo.md",
     "Consequence, which is the operational payoff: **a truncation-side price is a\n"
     "property of the matrix and a preconditioner reaches it; a continuum-side price\n"
     "is a property of the symbol and no congruence of the finite section can touch\n"
     "it.** That is the mechanism behind the v5.11.0 split — conditioning BREACHED,\n"
     "locality STANDING-but-capped.\n\n"
     "*Second route for `c_1`:* the exact tridiagonal spectrum, v5.10.18. Two\n"
     "representations, one constant (independent-route rule).",
     "**WITHDRAWN, same day (C23 run #3).** This section originally drew an\n"
     "operational corollary — that a truncation-side price is matrix-level and a\n"
     "preconditioner reaches it, while a continuum-side price is symbol-level and no\n"
     "congruence can touch it. Both halves are wrong. The v5.11.0 preconditioner is\n"
     "built FROM the symbol (matching polynomial chosen to share the symbol's zero,\n"
     "Serra), and preconditioning IS a congruence, replacing the symbol by `f/g`.\n"
     "The surviving statement is about the symbol on both sides and turns on the\n"
     "KIND of feature: a banded congruence cancels a finite-order zero but cannot\n"
     "alter a decay class — a second-order zero at one pole, a `j^-5/4` envelope at\n"
     "the other. **Open and load-bearing:** whether the two poles are independent\n"
     "features or two faces of one is unsettled; `debug/lit_scan/`"
     "`toeplitz_finite_section_memo.md` records both conclusions at different\n"
     "points, and the split adopted here is the weaker, measured one."),
]


def main() -> int:
    applied = skipped = 0
    for path, old, new in EDITS:
        with open(path, encoding="utf-8") as fh:
            t = fh.read()
        if new[:50] in t:
            skipped += 1
            continue
        if t.count(old) != 1:
            print(f"  MISS {path}: anchor count={t.count(old)}")
            print(f"       {old[:80]!r}")
            continue
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t.replace(old, new))
        applied += 1
        print(f"  ok   {path}")
    print(f"applied={applied} skipped={skipped} of {len(EDITS)}")
    return 0 if applied + skipped == len(EDITS) else 1


if __name__ == "__main__":
    sys.exit(main())
