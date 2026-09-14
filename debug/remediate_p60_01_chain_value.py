"""REMEDIATION 1/5 -- the wrong number in the abstract.

/qa paper_60 FULL 2026-09-12, claims-atomic M1 (LARGE), verified three ways by
the PM before editing:

  (a) the registry's own alias for K=164 is 7.289 mHa above the exact
      -2.903724377 Ha, i.e. E = -2.896435, which rounds to -2.896;
  (b) the paper's own K=244 pool value (7.06 mHa) plus its own eq:no_selection
      interlacing theorem forces E(164) >= -2.896667, because K=164 is a nested
      sub-family of K=244.  The printed -2.897 VIOLATES that bound; -2.8964
      satisfies it;
  (c) it sits in the abstract under a provenance paragraph certifying that every
      value in the section was recomputed on a converged radial domain.

Two loci, body and abstract.  The value is also REGISTERED here, because it was
an unregistered literal and therefore invisible to C21 -- the same class as the
4.43/4.40 mismatched pair caught on 2026-09-11.

NOT done here, and raised to the PI instead: the two earlier rungs (-2.873 s,
-2.894 +p) are unverified.  Only K=164 was confirmed.  Re-measuring them needs
the ladder driver and is its own pass.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
REG = "debug/qa/numeric_registry.py"
MARKER = "p60_he_chain_spdf_k164"

BODY_OLD = r"-2.847\ (1s^2)\;\to\;-2.873\ (s)\;\to\;-2.894\ (+p)\;\to\;-2.897\ (spdf,\,K{=}164),"
BODY_NEW = (r"-2.847\ (1s^2)\;\to\;-2.873\ (s)\;\to\;-2.894\ (+p)\;\to\;"
            r"\gvq{p60_he_chain_spdf_k164}{-2.896}\ (spdf,\,K{=}164),")

ABS_OLD = """the $s$, $+p$ and $spdf$ sectors---past the locked $s$-sector value $-2.873$~Ha to
$-2.897$~Ha, toward the exact non-relativistic $-2.90372$~Ha)"""
ABS_NEW = """the $s$, $+p$ and $spdf$ sectors---past the locked $s$-sector value $-2.873$~Ha to
$\\gvq{p60_he_chain_spdf_k164}{-2.896}$~Ha, toward the exact non-relativistic $-2.90372$~Ha)"""

REG_ANCHOR = '    "p60_window_richardson_pi2": dict(\n'
REG_NEW = '''    "p60_he_chain_spdf_k164": dict(
        value=-2.8964, convention="constant: Ha, the He ground-state energy "
                                  "reached by the LOCKED metric-free isoenergetic "
                                  "posing on the full s+p+d+f Goscinskian family "
                                  "at K=164 -- the endpoint of the convergence "
                                  "chain quoted in the abstract and in sec:atomic. "
                                  "Printed to three decimals as -2.896",
        q=None,
        provenance="DERIVED 2026-09-12 from the registered gap at the same K: "
                   "the extended-ladder alias gives 7.289 mHa above the exact "
                   "-2.903724377 Ha, so E = -2.896435. REGISTERED BECAUSE THE "
                   "PAPER CARRIED -2.897 -- a retired 60-bohr-domain value -- at "
                   "TWO loci (abstract and sec:atomic) as an unregistered literal "
                   "C21 could not see. It was also self-refuting: K=164 is a "
                   "nested sub-family of the K=244 pool whose own error is 7.06 "
                   "mHa, so eq:no_selection's interlacing forces E(164) >= "
                   "-2.896667, which -2.897 violates. Found by /qa paper_60 FULL "
                   "2026-09-12. OWED: the two earlier rungs (-2.873 s, -2.894 "
                   "+p) are NOT verified and are candidates for the same defect.",
        aliases={-2.895688: "K=74, from the 8.036 mHa alias"}),
'''


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    with open(REG, encoding="utf-8") as fh:
        r = fh.read()
    if MARKER in r and MARKER in t:
        print("ALREADY APPLIED")
        return 1
    for name, s, hay in (("body", BODY_OLD, t), ("abstract", ABS_OLD, t),
                         ("registry", REG_ANCHOR, r)):
        if hay.count(s) != 1:
            print(f"  {name} anchor count={hay.count(s)}; ABORT")
            return 2
    t = t.replace(BODY_OLD, BODY_NEW).replace(ABS_OLD, ABS_NEW)
    r = r.replace(REG_ANCHOR, REG_NEW + REG_ANCHOR)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    with open(REG, "w", encoding="utf-8") as fh:
        fh.write(r)
    print("applied: -2.897 -> -2.896 at 2 loci, registered as p60_he_chain_spdf_k164")
    return 0


if __name__ == "__main__":
    sys.exit(main())
