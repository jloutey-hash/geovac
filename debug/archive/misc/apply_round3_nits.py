r"""Round-3 code-review NITs: the two loose assertions in the headline guard,
and the omitted alpha = 1.00 value in Paper 12.

NIT-1  `len(pi_var) >= 2` of a hand-picked three-survivor alpha set tolerates
       one of the three being lost, and can never observe an envelope that
       WIDENS -- which is the direction the message claims to watch.  All
       three of these alpha are documented survivors, so the honest assertion
       is that all three survive.

NIT-2  `92.0 < _de_pct(e_sigma) < 93.0` is a one-point-wide band on a value
       the paper states to four significant figures (92.42).

PAPER  The envelope paragraph names alpha = 1.00 as failing but gives values
       only for 0.95, 1.15 and 1.20, omitting that the natural default returns
       -64.26 Ha -- the second-worst point on the grid, and the one a reader is
       most likely to try.  MEASURED here before writing it (N = 144, thresh
       1e-11, 102 of 144 vectors kept): -64.2611 Ha.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

TAZ = "tests/test_paper12_azimuthal_channels.py"
P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


edit(TAZ,
     '''    assert len(pi_var) >= 2, (
        f"only {len(pi_var)} of {len(pi_pts)} alpha points are variational at "
        f"|m|<=1, N=144; the conditioning envelope has widened beyond what "
        f"the paper documents"
    )''',
     '''    # All three of these alpha are DOCUMENTED survivors, so the assertion is
    # that all three survive.  An earlier version required only two of three,
    # which tolerated losing a survivor and -- the direction the message
    # actually claims to watch -- could never observe the envelope WIDENING.
    assert len(pi_var) == len(pi_pts), (
        f"only {len(pi_var)} of {len(pi_pts)} alpha points are variational at "
        f"|m|<=1, N=144; all three ({[a for _, a in pi_pts]}) are documented "
        f"survivors, so the conditioning envelope has widened beyond what the "
        f"paper documents"
    )
    assert len(sigma_var) == len(sigma_pts), (
        f"only {len(sigma_var)} of {len(sigma_pts)} alpha points are "
        f"variational at sigma-only, N=72, which the paper treats as the "
        f"well-conditioned case"
    )''',
     "NIT-1: the survivor count is now exact, in both sectors")

edit(TAZ,
     '''    assert 92.0 < _de_pct(e_sigma) < 93.0, (
        f"sigma-only at (3,3) is {_de_pct(e_sigma):.2f}%, paper says 92.42%"
    )''',
     '''    # Tightened from the original one-point-wide [92.0, 93.0], which was a
    # band a hundred times looser than the precision the paper states.
    assert 92.35 < _de_pct(e_sigma) < 92.50, (
        f"sigma-only at (3,3) is {_de_pct(e_sigma):.2f}%, paper says 92.42%"
    )''',
     "NIT-2: the sigma band matched to the stated precision")

edit(P12,
     r"""points fail---including $\alpha = 1.00$;\ $\alpha = 1.15$ and
$1.20$ return $-4.1$ and $-8.1$~Ha, and $\alpha = 0.95$ returns
$-277$~Ha.""",
     r"""points fail.  The failures are not marginal:\ $\alpha = 1.15$ and $1.20$
return $-4.1$ and $-8.1$~Ha, $\alpha = 0.95$ returns $-277$~Ha, and
$\alpha = 1.00$---the natural default, and the value a reader is most
likely to try first---returns $-64.3$~Ha.""",
     "PAPER: the alpha = 1.00 value stated, not just named as failing")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    for old, new, label in items:
        if old in t:
            t = t.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s   [%s]" % (label, path))
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)

print("applied %d of %d" % (len(applied), len(EDITS)))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
