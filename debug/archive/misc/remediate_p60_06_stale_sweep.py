"""REMEDIATION 6/6 -- the stale-text sweep across code, tests, and the synthesis.

/qa paper_60 FULL 2026-09-12.  All five loci carry a reading the owner document
corrected and the citer never received.

  A. SYNTHESIS molecular verdict is PRE-BREACH (synthesis M1, LARGE) + claims
     the KMS law as "derived" (synthesis M2) + names THREE levers and calls the
     gerade one "the sharpest", which sec:resource now contradicts in terms.
  B. SYNTHESIS ladder reductions bound to the wrong basis point (synthesis M3).
  C. geovac/sturmian_l2_encoding.py docstring carries the WITHDRAWN
     ill-conditioning mechanism (critic G1) -- C8#1 declares exactly this
     MATERIAL, and 0 of 63 C16 entries fire on the file.
  D. geovac/sturmian_sigma_law.py + tests/test_paper60_sigma_law.py present the
     KMS asymptotic as internal work (critic G2/G4).
  E. tests/test_paper60_sturmian.py carries the RETIRED rising-slope sequence
     ("0.850 at K=202 ... 0.906 at 340") and points at a test deleted on
     2026-09-07 (code-B M1); docs/claim_test_matrix.md still cites that deleted
     test as BACKED-SOUND (code-B M3).

Idempotent.
"""
from __future__ import annotations

import sys

SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
ENC = "geovac/sturmian_l2_encoding.py"
SIG = "geovac/sturmian_sigma_law.py"
TSIG = "tests/test_paper60_sigma_law.py"
TSTU = "tests/test_paper60_sturmian.py"
MTX = "docs/claim_test_matrix.md"

EDITS = [
 (SYN, "A-synthesis-verdict",
  """the conditioning grows
polynomially, with a derived band-limited law, and of its three levers
(the gerade sector, large separation, a metric confined to one electron)
the sharpest---the gerade sector's basis-independent conditioning---is an
\\emph{equivalent-center} property, and a symmetry-unique heavy atom
reinstates the \\emph{growth}, confining that lever to
homonuclear-diatomic-like systems.""",
  """the conditioning grows
polynomially, at a band-limited concentration rate that is the
Kac--Murdock--Szeg\\H{o} asymptotic (what is the paper's own is the
\\emph{identification} of the metric as such a finite section, not the law).
Of its four levers, the gerade sector's basis-independent conditioning is an
\\emph{equivalent-center} property, and a symmetry-unique heavy atom reinstates
the raw \\emph{growth} there;\\ but a band-Toeplitz preconditioner removes the
conditioning cost without requiring equivalent centers, reaching water's $A_1$
block, and a direct block-encoding of the preconditioned metric takes the
penalty from $n^3$ to $n$.  Established for $M=2$ and non-collinear $M=3$;\\ the
collinear case is open.  What remains standing is the \\emph{locality} cost,
capped by the symbol's other pole, and the $\\ell$-block sparsity loss, which is
structural rather than a conditioning effect."""),

 (ENC, "C-encoding-docstring",
  """    sturmian   : lambda ~ Q^3.33  (the shared-scale overlap ill-conditions,
                 and Loewdin whitening spreads that ill-conditioning into the
                 transformed integrals -- the motivation for the isoenergetic
                 reformulation, which removes the metric for atoms).""",
  """    sturmian   : lambda ~ Q^3.33  (the shared-scale overlap is strongly
                 NON-ORTHOGONAL, so Loewdin's S^-1/2 is DENSE and spreads the
                 coefficient distribution -- the motivation for the isoenergetic
                 reformulation, which removes the metric for atoms).

    MECHANISM, corrected 2026-09-07 [retracted 2026-09-12: p60-l2-metric-diverges]:
    this is NOT ill-conditioning.  The converged cond(S) is ordinary (~0.12*K,
    4.07 at K=10, 56.3 at K=452); the earlier "ill-conditions" reading was an
    artifact of a fixed 60-bohr radial domain.  The inflation is driven by the
    DENSITY of S^-1/2, not by numerical instability."""),

 (SIG, "D-sigma-module-docstring",
  """i.e. the conditioning exponent is exactly 2 (asymptotically; finite windows fit
lower slopes such as the N^1.85 / N^1.97 reported in Paper 60, which are
pre-asymptotic readings of this one law).""",
  """i.e. the conditioning exponent is exactly 2 (asymptotically; finite windows fit
lower slopes such as the N^1.85 / N^1.97 reported in Paper 60, which are
pre-asymptotic readings of this one law).

PRIOR ART (recorded 2026-09-11): this asymptotic is NOT derived here.  It is the
Kac-Murdock-Szego extreme-eigenvalue law (c_1 = pi^2, 1953; normal form in
Boettcher-Widom arXiv:math/0412269).  What is ours is the IDENTIFICATION -- that
the two-center Shibuya-Wulfman metric in the sine basis IS such a finite
section, Toeplitz minus Hankel with symbol j0(kR cot(chi/2)).  Do not describe
the law as derived here."""),

 (TSIG, "D-sigma-test-docstring",
  "Backs the derived conditioning law (CHANGELOG v4.103.0; Paper 60 sec:molecular):",
  """Backs the conditioning law (CHANGELOG v4.103.0; Paper 60 sec:molecular).
The ASYMPTOTIC is Kac-Murdock-Szego prior art, not a GeoVac derivation; what is
backed here is the identification and the measured spectrum:"""),

 (TSTU, "E-test-zombie",
  """  # NOT an asymptote (corrected 2026-09-07, /qa): 0.84 is a fit over the WINDOW
    # K = 74..164, and the local slope keeps rising past it -- 0.850 at K=202,
    # 0.868 at 244, 0.882 at 290, 0.906 at 340.  This self-contained sweep reaches
    # ~0.77 at K<=24; the pinned claim is the sublinear band over the computed
    # range, and NOT any asymptotic value.  The mechanism is pinned separately by
    # test_paper60_sublinearity_is_carried_by_the_nuclear_diagonal.""",
  """  # NOT an asymptote (corrected 2026-09-07; direction re-measured 2026-09-08
    # [retracted 2026-09-12: p60-tprime-superlinear]): the exponent is a WINDOW
    # fit, and on a CONVERGED box the local slope FALLS monotonically -- 0.827 at
    # K=74 down to 0.766 by K=514 -- so no window fit is stable.  The earlier
    # "keeps rising to 0.906" reading was measured on a truncated radial domain
    # and is retired.  This self-contained sweep reaches ~0.77 at K<=24; the
    # pinned claim is the sublinear band over the computed range, and NOT any
    # asymptotic value."""),
]


def main() -> int:
    applied, bad = 0, 0
    loaded = {}
    for path, name, old, new in EDITS:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        marker = {"C-encoding-docstring": "MECHANISM, corrected 2026-09-07",
                  "D-sigma-module-docstring": "PRIOR ART (recorded 2026-09-11)"}.get(name, new[:45])
        if marker in t:
            print(f"  skip {name} (already applied)")
            continue
        if t.count(old) != 1:
            print(f"  MISS {name}: count={t.count(old)}")
            bad += 1
            continue
        loaded[path] = t.replace(old, new)
        applied += 1
        print(f"  ok   {name}")
    if bad:
        print("ABORT -- no file written")
        return 2
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)

    # E (second half): the claim-matrix row citing a deleted test
    with open(MTX, encoding="utf-8") as fh:
        m = fh.read()
    dead = "::test_paper60_sublinearity_is_carried_by_the_nuclear_diagonal"
    if dead in m:
        m = m.replace(
            dead,
            "::test_paper60_lgt0_angular_correlation_and_sublinear`` "
            "(**corrected 2026-09-12: the row cited "
            "`test_paper60_sublinearity_is_carried_by_the_nuclear_diagonal`, "
            "retired 2026-09-07 with the T^0 mechanism it pinned; a labelled "
            "equation's claim row was pointing at nothing**)")
        with open(MTX, "w", encoding="utf-8") as fh:
            fh.write(m)
        applied += 1
        print("  ok   E-matrix-dangling-row")
    print(f"applied {applied}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
