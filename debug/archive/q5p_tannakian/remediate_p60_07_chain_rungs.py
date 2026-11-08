"""OWED ITEM 1 -- recheck the two unverified He-chain rungs.

The /qa FULL run confirmed only the K=164 endpoint was stale.  The two middle
rungs were flagged as candidates for the same defect.  Recomputed here on a
converged grid (box 500, 40000 points -- the earlier apparent box drift was a
GRID artifact: at 40000 points the value is stable to 2e-5 across boxes
300..1200).  All four rungs taken at a COMMON n_max = 10, which is the basis the
chain's own last entry names (K=164):

    rung    l_max    K     E (Ha)        paper printed
    1s^2      -      1    -2.847651      -2.847   OK (truncation)
    s         0     55    -2.874468      -2.873   WRONG  <-- one more found
    +p        1    100    -2.894672      -2.894   OK (truncation)
    spdf      3    164    -2.896432      -2.896   OK (fixed earlier today)

The s value the paper printed, -2.873, is the n_max = 4 (K=10) value
(-2.873219).  So the chain silently MIXED basis sizes -- a violation of this
paper's own branch-defining criterion, which requires every quantity to name its
evaluation domain and basis-growth family.

Two fixes, not one:
  (i)  correct the s rung;
  (ii) print four decimals and label every rung's K, so the truncation
       convention cannot be misread (-2.894672 printed as "-2.894" looks wrong
       to a reader who rounds) and the common basis is stated.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
REG = "debug/qa/numeric_registry.py"
MARKER = "p60_he_chain_s_k55"

BODY_OLD = (r"-2.847\ (1s^2)\;\to\;-2.873\ (s)\;\to\;-2.894\ (+p)\;\to\;"
            r"\gvq{p60_he_chain_spdf_k164}{-2.896}\ (spdf,\,K{=}164),")
BODY_NEW = (r"-2.8477\ (1s^2)\;\to\;\gvq{p60_he_chain_s_k55}{-2.8745}\ (s,\,K{=}55)"
            r"\;\to\;-2.8947\ (+p,\,K{=}100)\;\to\;"
            r"\gvq{p60_he_chain_spdf_k164}{-2.8964}\ (spdf,\,K{=}164),")

CTX_OLD = "converging toward the exact non-relativistic $-2.90372$~Ha."
CTX_NEW = ("converging toward the exact non-relativistic $-2.90372$~Ha.  Every "
           "rung is taken at the same $n_{\\max}=10$ family, so the sequence "
           "is a pure $\\ell_{\\max}$ ladder;\\ values are truncated, not "
           "rounded.")

ABS_OLD = ("the $s$, $+p$ and $spdf$ sectors---past the locked $s$-sector value "
           "$-2.873$~Ha to\n$\\gvq{p60_he_chain_spdf_k164}{-2.896}$~Ha")
ABS_NEW = ("the $s$, $+p$ and $spdf$ sectors---past the locked $s$-sector value "
           "$\\gvq{p60_he_chain_s_k55}{-2.8745}$~Ha to\n"
           "$\\gvq{p60_he_chain_spdf_k164}{-2.8964}$~Ha")

REG_ANCHOR = '    "p60_he_chain_spdf_k164": dict(\n'
REG_NEW = '''    "p60_he_chain_s_k55": dict(
        value=-2.8745, convention="constant: Ha, the He ground-state energy from "
                                  "the LOCKED metric-free isoenergetic posing on "
                                  "the s-only (l_max=0) Goscinskian family at "
                                  "n_max=10, K=55 -- the second rung of the "
                                  "convergence chain. Printed truncated to four "
                                  "decimals",
        q=None,
        provenance="MEASURED 2026-09-12 on a converged grid (box 500, 40000 pts): "
                   "-2.874468, i.e. 29.256 mHa above the exact -2.903724377. "
                   "REGISTERED BECAUSE THE PAPER CARRIED -2.873, which is the "
                   "n_max=4 (K=10) value (-2.873219) -- so the chain silently "
                   "MIXED basis sizes, against this paper's own requirement that "
                   "every quantity name its basis-growth family. Found by the "
                   "owed-items recheck after /qa paper_60 FULL 2026-09-12, which "
                   "had confirmed only the K=164 endpoint. Grid note: an apparent "
                   "box drift (-2.8744 -> -2.8739 over boxes 300..1200) is a GRID "
                   "artifact; at 40000 points the value is stable to 2e-5 across "
                   "the same boxes.",
        aliases={-2.873219: "n_max=4, K=10 -- the value the paper had printed",
                 -2.894672: "the +p rung, n_max=10, K=100",
                 -2.847651: "the 1s^2 single-configuration rung"}),
'''


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    with open(REG, encoding="utf-8") as fh:
        r = fh.read()
    if MARKER in r and MARKER in t:
        print("ALREADY APPLIED")
        return 1
    for nm, s, hay in (("body", BODY_OLD, t), ("ctx", CTX_OLD, t),
                       ("abstract", ABS_OLD, t), ("registry", REG_ANCHOR, r)):
        if hay.count(s) != 1:
            print(f"  {nm} anchor count={hay.count(s)}; ABORT")
            return 2
    t = t.replace(BODY_OLD, BODY_NEW).replace(CTX_OLD, CTX_NEW).replace(ABS_OLD, ABS_NEW)
    r = r.replace(REG_ANCHOR, REG_NEW + REG_ANCHOR)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    with open(REG, "w", encoding="utf-8") as fh:
        fh.write(r)
    print("applied: s rung -2.873 -> -2.8745, all rungs labelled with K, "
          "common basis stated, s value registered")
    return 0


if __name__ == "__main__":
    sys.exit(main())
