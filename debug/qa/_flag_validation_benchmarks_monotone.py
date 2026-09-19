r"""Flag the graph-native He n_max ladder as internally inconsistent.

THE DEFECT.  docs/validation_benchmarks.md lines 139-141 read:

    | Graph-native CI n_max=7 error | < 0.20%  | Graph-native CI accuracy |
    | Graph-native CI n_max=8 error | 0.207% (2,262 configs) | ... |
    | Graph-native CI n_max=9 error | 0.201% (3,927 configs) | ... |

A variational ladder converging from above is monotone decreasing, so it
cannot be below 0.20 at n_max=7 and then 0.207 at n_max=8.  At least one of
the three rows is wrong.

WHY THIS IS FLAGGED AND NOT FIXED.  Deciding WHICH row is wrong needs a
graph-native build at n_max=6/7, which is out of reach here: n_max=5 alone
costs ~68 s and `hypergeometric_slater` dispatches to exact `Fraction`
arithmetic at n >= 5, and two attempts at n_max=6/7 during the /qa run died
without output.  The reproduced ladder is 5.294 / 0.4906 / 0.36555 / 0.28649 /
0.24963 % at n_max = 1..5, whose decrements shrink by ~0.47 each, extrapolating
to roughly 0.22-0.23 % at n_max=7 -- which is consistent with 0.201-0.207 at
n_max=8/9 being right and the "< 0.20%" at n_max=7 being the wrong row, but an
extrapolation is not a measurement and this registry does not record guesses
(numeric_registry.py: "a registry of guesses launders them into authority").

So: record the inconsistency, name what would settle it, and leave the values
alone.  The 0.19 % headline itself is ALREADY a declared coverage gap --
docs/claim_test_matrix.md:280 (NO-TEST, "tests stop at n_max=3"),
tests/test_casimir_ci.py:1008 ("NO-TEST (deliberate)"), and
benchmarks/certified_reference/entries_helium.py (N_MAX_LIST = (1,2,3,4,5),
certifying assembly rather than accuracy).  This adds the arithmetic
impossibility that those three records do not mention.

Run:  python debug/qa/_flag_validation_benchmarks_monotone.py
"""
from __future__ import annotations

import io
import sys

PATH = "docs/validation_benchmarks.md"

OLD = (
    "| Graph-native CI n_max=7 error | < 0.20% | Graph-native CI accuracy |\n"
    "| Graph-native CI n_max=8 error | 0.207% (2,262 configs) | Exact algebraic float integrals |\n"
    "| Graph-native CI n_max=9 error | 0.201% (3,927 configs) | Exact algebraic float integrals |\n"
)

NEW = (
    "| Graph-native CI n_max=7 error | < 0.20% ⚠ | Graph-native CI accuracy — **see inconsistency note below** |\n"
    "| Graph-native CI n_max=8 error | 0.207% (2,262 configs) ⚠ | Exact algebraic float integrals |\n"
    "| Graph-native CI n_max=9 error | 0.201% (3,927 configs) ⚠ | Exact algebraic float integrals |\n"
    "\n"
    "> **⚠ INTERNAL INCONSISTENCY, flagged 2026-09-19 (/qa group2 CODE run; NOT resolved).**\n"
    "> These three rows cannot all be right. A variational ladder converging from above is\n"
    "> monotone decreasing, so the error cannot be `< 0.20%` at n_max=7 and then `0.207%` at\n"
    "> n_max=8. At least one row is wrong, and **which one is unmeasured**: a graph-native\n"
    "> build at n_max=6/7 is expensive (n_max=5 alone ~68 s; `hypergeometric_slater`\n"
    "> dispatches to exact `Fraction` at n≥5) and two attempts during the run produced no\n"
    "> output. The reproduced ladder n_max=1..5 is 5.294 / 0.4906 / 0.36555 / 0.28649 /\n"
    "> 0.24963 %, whose shrinking decrements extrapolate to ~0.22–0.23 % at n_max=7 —\n"
    "> suggesting the n_max=7 row is the wrong one, but an extrapolation is not a\n"
    "> measurement and no value is being substituted on its strength.\n"
    ">\n"
    "> Related, and already declared elsewhere: the **0.19% @ n_max=7** headline has **no\n"
    "> backing test** — `docs/claim_test_matrix.md:280` (NO-TEST, tests stop at n_max=3),\n"
    "> `tests/test_casimir_ci.py:1008` (\"NO-TEST (deliberate)\"), and\n"
    "> `benchmarks/certified_reference/entries_helium.py` (`N_MAX_LIST = (1,2,3,4,5)`, which\n"
    "> certifies assembly, not accuracy). **To settle both:** one long-running\n"
    "> `build_graph_native_fci(Z=2, n_max=7)` job outside a QA pass.\n"
)


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    s = io.open(PATH, encoding="utf-8").read()
    if "INTERNAL INCONSISTENCY, flagged 2026-09-19" in s:
        print("note already present; nothing to do")
        return
    if OLD not in s:
        i = s.find("Graph-native CI n_max=7")
        print("ANCHOR NOT FOUND. Actual text:")
        print(repr(s[max(0, i - 120):i + 420]) if i >= 0 else "  (phrase absent)")
        raise SystemExit(1)
    io.open(PATH, "w", encoding="utf-8", newline="").write(s.replace(OLD, NEW, 1))
    print("validation_benchmarks.md: monotone inconsistency flagged (values untouched)")


if __name__ == "__main__":
    main()
