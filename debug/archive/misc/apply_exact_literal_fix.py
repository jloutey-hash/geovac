"""Fix the dropped-digit Kolos-Wolniewicz literal in tests/test_neumann_vee.py.

-1.17475 should be -1.174475.  The papers corrected this on 2026-09-13
("KW reference literal corrected -1.17475 -> -1.174475", claim_test_matrix
line for Paper 12), but the test kept the old value -- and it is the
DENOMINATOR of the D_e percentage the test asserts.

Effect: D_e_exact was 0.17475 instead of 0.174475, a 0.16% relative error in
the denominator, so every percentage the test computed read ~0.16% low.  The
assertion band [90, 94] is wide enough that nothing failed, which is why it
survived; the printed diagnostic was wrong throughout.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

T = "tests/test_neumann_vee.py"

with io.open(T, encoding="utf-8") as fh:
    t = fh.read()

before = t.count("-1.17475")
if before == 0:
    print("already fixed")
    sys.exit(0)

t = t.replace(
    "        E_exact = -1.17475   # Kolos-Wolniewicz, Paper 12 eq via Kolos1968\n"
    "        E_atoms = -1.0       # two H atoms, 2 x (-0.5) Ha\n"
    "        D_e_exact = E_atoms - E_exact   # = 0.17475 Ha",
    "        # Kolos-Wolniewicz.  NOTE the fourth decimal: -1.174475, not\n"
    "        # -1.17475.  The dropped digit made D_e_exact 0.17475 instead of\n"
    "        # 0.174475 -- a 0.16% error in the DENOMINATOR of every percentage\n"
    "        # below.  The [90, 94] band was wide enough to hide it.\n"
    "        E_exact = -1.174475  # Kolos & Wolniewicz, via Kolos1968\n"
    "        E_atoms = -1.0       # two H atoms, 2 x (-0.5) Ha\n"
    "        D_e_exact = E_atoms - E_exact   # = 0.174475 Ha",
    1)

t = t.replace("        E_exact = -1.17475\n",
              "        E_exact = -1.174475  # Kolos & Wolniewicz\n", 1)

after = t.count("-1.17475")
with io.open(T, "w", encoding="utf-8") as fh:
    fh.write(t)
print("replaced %d of %d occurrences of the dropped-digit literal"
      % (before - after, before))
if after:
    print("STILL PRESENT: %d" % after)
    sys.exit(1)
