"""Review the Paper-15 test file as a declared dependent of the withdrawn
Paper-12 diagnosis, and correct the docstrings that restate it.

The ASSERTIONS are untouched and remain valid: 92.4 < pct < 100 is still the
right numerical band, and still does double duty against the adiabatic
false-positive.  What changes is the interpretation the docstrings carry --
"exceeding Paper 12" was read as a coordinate-system advantage, and that
reading is withdrawn because the 94.1% includes pi channels while the 92.4%
does not.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

PATH = "tests/test_level4_multichannel.py"
EDITS = []


def edit(old, new, label):
    EDITS.append((old, new, label))


edit(
    '''        """Paper 15 Table III headline: 2D-variational sigma+pi (l_max=4, 29 ch)
        recovers 94.1% of D_e -- EXCEEDING Paper 12's 92.4%.

        Uses the 2D variational solver (n_coupled=-1), the paper's headline
        solver.  Previously this test ran the DEFAULT adiabatic solver
        (n_coupled=1), which returns 105.2% here (Table II) -- a
        variational-bound violation the paper disavows.  The band
        92.4 < pct < 100 therefore does double duty: it (a) confirms the
        paper's "exceeds Paper 12" claim and (b) rejects the adiabatic
        false-positive by enforcing the variational bound D_e <= D_e_exact.

        Recompute: D_e = 0.1642 Ha (94.1%), matching the paper exactly.
        """''',
    '''        """Paper 15 Table III headline: 2D-variational sigma+pi (l_max=4, 29 ch)
        recovers 94.1% of D_e, above the 92.4% of Paper 12's sigma-only CI.

        Uses the 2D variational solver (n_coupled=-1), the paper's headline
        solver.  Previously this test ran the DEFAULT adiabatic solver
        (n_coupled=1), which returns 105.2% here (Table II) -- a
        variational-bound violation the paper disavows.  The band
        92.4 < pct < 100 therefore does double duty: it (a) pins the headline
        value and (b) rejects the adiabatic false-positive by enforcing the
        variational bound D_e <= D_e_exact.

        SCOPE (2026-09-14).  The name of this test is historical.  Clearing
        92.4% is NOT evidence of a coordinate-system advantage over prolate
        spheroidal coordinates, and Paper 15 no longer claims one: the 94.1%
        here includes pi channels and Paper 12's 92.4% does not -- its basis
        is phi-independent.  At matched angular content the ordering reverses
        (Paper 12's own basis with |m| <= 1 reaches 99.09%).  The assertion
        below is unaffected; only the reading is.

        Recompute: D_e = 0.1642 Ha (94.1%), matching the paper exactly.
        """''',
    "test_mmax1_exceeds_paper12 docstring")

edit(
    '''        shift (~ -0.39 mHa at R_eq for l_max=4), so it raises D_e from 94.1%
        to 94.3%, above Paper 12's 92.4%.''',
    '''        shift (~ -0.39 mHa at R_eq for l_max=4), so it raises D_e from 94.1%
        to 94.3%, above the 92.4% of Paper 12's sigma-only CI (not a
        like-for-like comparison -- see test_mmax1_exceeds_paper12).''',
    "cusp-correction docstring")

edit(
    '''        87.0% of D_e -- BELOW Paper 12's 92.4%.''',
    '''        87.0% of D_e -- BELOW the 92.4% of Paper 12's sigma-only CI.''',
    "sigma-only docstring")

with io.open(PATH, encoding="utf-8") as fh:
    text = fh.read()

applied, failed = [], []
for old, new, label in EDITS:
    if old in text:
        text = text.replace(old, new, 1)
        applied.append(label)
    else:
        failed.append(label)

if failed:
    print("FAILED TO MATCH:")
    for f in failed:
        print("  -", f)
    sys.exit(1)

with io.open(PATH, "w", encoding="utf-8") as fh:
    fh.write(text)

print("applied %d docstring edits (no assertion changed)" % len(applied))
for a in applied:
    print("  +", a)
