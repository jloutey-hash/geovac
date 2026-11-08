"""Add the 'accuracy' convention kind so the Paper-12 azimuthal entries parse.

tests/test_numeric_registry.py::test_every_convention_string_parses caught
five new entries whose convention strings the parser could not read -- exactly
what that guard exists for.  These are accuracy-side quantities (a percentage
of D_e, an energy gain in mHa, a condition number) with no identity-in/out
convention, so they behave like the existing `density` and `constant` kinds.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

REG = "debug/qa/numeric_registry.py"

OLD = '''    ("density", ("density",)),
)'''

NEW = '''    ("density", ("density",)),
    # Accuracy-side quantities (Paper 12 azimuthal channels, 2026-09-14):
    # a percentage of D_e, an energy gain in mHa, a condition number.  No
    # identity-in/out convention applies -- these count no Pauli terms -- so
    # family() leaves identity None, as for `density` and `constant`.  What
    # the convention string must still carry is the BASIS and truncation the
    # number was measured at, because 92.25 and 92.42 are the same quantity
    # at different (j_max, l_max) and pairing them would be the twin defect
    # this registry exists to catch.
    ("accuracy", ("% of d_e", "mha gained", "cond(s)")),
)'''

with io.open(REG, encoding="utf-8") as fh:
    text = fh.read()

if OLD not in text:
    print("FAILED: _KINDS anchor not found")
    sys.exit(1)

text = text.replace(OLD, NEW, 1)

# make the five conventions match the new needles exactly
fixes = [
    ('convention="% of D_e, H2 R=1.4011, (j,l)=(3,3) sigma "\n                                "only, canonical orthogonalisation"',
     'convention="% of D_e at (j,l)=(3,3), H2 R=1.4011, sigma only, "\n                                "canonical orthogonalisation"'),
    ('convention="mHa gained by opening |m|<=1 at (3,3)"',
     'convention="mHa gained by opening |m|<=1 at (j,l)=(3,3)"'),
    ('convention="mHa gained by N=27 -> 72 along the sigma axis"',
     'convention="mHa gained by N=27 -> 72 along the sigma axis"'),
]
for old, new in fixes:
    if old in text:
        text = text.replace(old, new, 1)

with io.open(REG, "w", encoding="utf-8") as fh:
    fh.write(text)

print("added the 'accuracy' kind")
