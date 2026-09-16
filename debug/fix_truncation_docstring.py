"""Correct the truncation guard's docstring to match what the fire test showed.

Fire test 2026-09-14: removing the eta selection rule ALONE does not fire, and
removing the l cap ALONE does not fire; removing BOTH fires.  The two are
genuinely redundant, so the guard's subject is the conjunction, not either one.
The original docstring said "REJECTS: removal of the eta selection rule / the
l cap", which overstates what it discriminates.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

PATH = "tests/test_paper12_azimuthal_channels.py"

OLD = '''    """REJECTS: removal of the eta selection rule / the l cap.

    Without them the ~1e10 Legendre-derivative coefficients leave a
    floating-point residue that the radial integral amplifies; the observed
    failure was E = -3.2e8 Ha at l_neumann = 18.  Bit-level invariance across
    a wide truncation range is the property that distinguishes an imposed
    selection rule from a lucky cancellation.
    """'''

NEW = '''    """REJECTS: loss of BOTH exactness mechanisms at once.

    Exactness is protected twice over -- by the eta selection rule
    (l > Q + 2s - m gives an identically zero moment) and by the cap that
    stops the l sum at that cutoff so no overflow-prone block is built.
    Without either, the ~1e10 Legendre-derivative coefficients leave a
    floating-point residue that the radial integral amplifies; the observed
    failure was E = -3.2e8 Ha at l_neumann = 18.

    Scope, established by fire test on 2026-09-14 rather than asserted:
    removing the selection rule alone does NOT fire, and removing the cap
    alone does NOT fire, because each covers for the other; removing both
    fires.  So this guard discriminates the conjunction.  A reviewer wanting
    single-point coverage would need one of the two mechanisms removed on
    purpose, which no caller has reason to do -- the redundancy is
    deliberate, and the test says so rather than implying a sharper claim
    than it makes.
    """'''

with io.open(PATH, encoding="utf-8") as fh:
    text = fh.read()

if OLD not in text:
    print("FAILED TO MATCH docstring")
    sys.exit(1)

with io.open(PATH, "w", encoding="utf-8") as fh:
    fh.write(text.replace(OLD, NEW, 1))

print("docstring corrected to match the fire test")
