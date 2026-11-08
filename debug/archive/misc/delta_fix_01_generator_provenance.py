"""DELTA F1 (LARGE) -- the prior-art credit lives only in a GENERATED file.

Claim-impact reviewer, verified by the PM before editing:

    docs/certified_reference_values.md   1 mention of Kac/Murdock  (the credit)
    benchmarks/.../entries_anchors.py    0
    benchmarks/.../certified_reference_values.json  0

The 2026-09-11 re-attribution was applied to the OUTPUT and never to the
generator.  Two live consequences, and the second is the worse one:

  (a) the distributed JSON -- the artifact this corpus explicitly offers
      outward as reference data -- already presents the Kac-Murdock-Szego
      asymptotic as an internal finding with no credit;
  (b) the next `python -m benchmarks.certified_reference.generate_table`
      SILENTLY DELETES the hand-applied credit from the markdown.

That is a regeneration bomb rather than a wording slip: the generated layer is
a standing eraser of corrections applied downstream of it.  `entry()` takes
`**extra`, so the fix carries the credit in the row itself, which reaches BOTH
artifacts and survives regeneration by construction.

Idempotent.
"""
from __future__ import annotations

import sys

G = "benchmarks/certified_reference/entries_anchors.py"
MARKER = "Kac-Murdock-Szego"

OLD = '''            "1 - sigma_max = (kR)^2 pi^2 / (24 n^2), so the rescaled quantity "
            "(1 - sigma_max)(n/kR)^2 tends to pi^2/24 as the per-centre basis "
            "size n grows.  This fixes the conditioning exponent at exactly 2, "
            "and reveals the previously fitted exponents 1.85 and 1.97 as "
            "pre-asymptotic windows of the same law.",'''

NEW = '''            "1 - sigma_max = (kR)^2 pi^2 / (24 n^2), so the rescaled quantity "
            "(1 - sigma_max)(n/kR)^2 tends to pi^2/24 as the per-centre basis "
            "size n grows.  This fixes the conditioning exponent at exactly 2, "
            "and reveals the previously fitted exponents 1.85 and 1.97 as "
            "pre-asymptotic windows of the same law.  PRIOR ART (recorded "
            "2026-09-11): the asymptotic is NOT derived in this corpus -- it is "
            "the Kac-Murdock-Szego extreme-eigenvalue law (c_1 = pi^2, J. "
            "Rational Mech. Anal. 2, 767 (1953); normal form in Boettcher-Widom "
            "arXiv:math/0412269), and pi^2/24 is c_1 times the symbol curvature "
            "b(1) = (kR)^2/24.  What is ours is the IDENTIFICATION of the "
            "two-centre Shibuya-Wulfman metric in the sine basis as such a "
            "finite section, Toeplitz minus Hankel with symbol "
            "j0(kR cot(chi/2)).  Do not describe the law as derived here.",'''

EXTRA_OLD = '''            transcendence_class="{pi^2} -- pure-Tate, weight 2",
            backing_test="tests/test_paper60_sigma_law.py",
        ))'''

EXTRA_NEW = '''            transcendence_class="{pi^2} -- pure-Tate, weight 2",
            backing_test="tests/test_paper60_sigma_law.py",
            provenance=("PRIOR ART: Kac, Murdock and Szego, J. Rational Mech. "
                        "Anal. 2, 767 (1953); c_1 = pi^2.  Ours is the "
                        "identification of the SW metric as such a finite "
                        "section, not the asymptotic."),
        ))'''


def main() -> int:
    with open(G, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    for nm, s in (("method", OLD), ("extra", EXTRA_OLD)):
        if t.count(s) != 1:
            print(f"  {nm} anchor count={t.count(s)}; ABORT")
            return 2
    t = t.replace(OLD, NEW).replace(EXTRA_OLD, EXTRA_NEW)
    with open(G, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("applied: prior-art credit moved into the GENERATOR row "
          "(method text + a provenance field)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
