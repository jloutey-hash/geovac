"""Close the F4 backing gap: no test pinned the 99.09% HEADLINE.

`test_azimuthal_channels_close_the_gap` asserts the (2,2)/N=54 row (98.96%);
the paper's headline is the (3,3)/N=144 value.  The claim matrix said
BACKED-SOUND for the headline, which was false.  This adds the missing
assertion as a @slow test and corrects the matrix row.

It matters more than usual here because N=144 is exactly the regime where the
paper records a direct eigensolve returning -79 Ha -- the headline's
trustworthiness rests on canonical orthogonalisation, and nothing exercised it.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

TEST = "tests/test_paper12_azimuthal_channels.py"
MATRIX = "docs/claim_test_matrix.md"

NEW_TEST = '''

@pytest.mark.slow
def test_headline_99_09_at_largest_basis():
    """REJECTS: a headline that only holds at the small truncation.

    The paper's abstract, conclusion, Paper 13, Paper 15, the group2 synthesis
    and docs/validation_benchmarks.md all carry 99.09%, which is the
    (j,l) = (3,3), |m| <= 1 value at N = 144 -- NOT the (2,2)/N = 54 value
    (98.96%) that the fast test above pins.  Until this test existed the
    headline was unpinned while the claim matrix read BACKED-SOUND.

    This is also the only regime the paper flags as numerically dangerous: at
    N = 144 cond(S) = 2.0e16 and a direct eigh(H, S) returns -79 Ha.  So the
    assertion below is doing two jobs -- pinning the published number, and
    standing guard over the canonical-orthogonalisation path that makes it
    meaningful.  A non-variational value here means that path regressed.
    """
    e_sigma, n_sigma = _energy(3, 3, 0, alpha=1.15)
    e_pi, n_pi = _energy(3, 3, 1, alpha=1.25)

    assert n_sigma == 72 and n_pi == 144, (
        f"basis sizes changed: {n_sigma}, {n_pi} (paper reports 72 and 144)"
    )
    assert e_pi > E_EXACT, (
        f"E = {e_pi:.6f} is BELOW the exact {E_EXACT:.6f}; the "
        f"canonical-orthogonalisation path has regressed (a direct solve "
        f"returns about -79 Ha at this basis)"
    )
    assert 92.0 < _de_pct(e_sigma) < 93.0, (
        f"sigma-only at (3,3) is {_de_pct(e_sigma):.2f}%, paper says 92.42%"
    )
    assert _de_pct(e_pi) > 98.9, (
        f"headline is {_de_pct(e_pi):.2f}% of D_e, paper claims 99.09%"
    )
    gain_mha = 1000.0 * (e_sigma - e_pi)
    assert gain_mha > 11.0, (
        f"azimuthal gain at the largest basis is {gain_mha:.2f} mHa, "
        f"paper claims 11.64"
    )
'''

with io.open(TEST, encoding="utf-8") as fh:
    text = fh.read()

marker = "# ======================================================================\n# 7. The conditioning caveat the paper states"
if marker not in text:
    print("FAILED: test anchor not found")
    sys.exit(1)
text = text.replace(marker, NEW_TEST.strip("\n") + "\n\n\n" + marker, 1)
with io.open(TEST, "w", encoding="utf-8") as fh:
    fh.write(text)
print("  + tests: test_headline_99_09_at_largest_basis (@slow)")

with io.open(MATRIX, encoding="utf-8") as fh:
    mat = fh.read()
old = ("`test_paper12_azimuthal_channels.py`: `test_azimuthal_channels_close_the_gap`, "
       "`test_gain_is_a_channel_effect_not_a_count_effect` (the discriminator), "
       "`test_sigma_only_reproduces_paper12` | **BACKED-SOUND** | new 2026-09-14.")
new = ("`test_paper12_azimuthal_channels.py`: `test_headline_99_09_at_largest_basis` "
       "(@slow -- pins the (3,3)/N=144 headline AND the canonical-orthogonalisation "
       "path), `test_azimuthal_channels_close_the_gap` (pins the (2,2) row), "
       "`test_gain_is_a_channel_effect_not_a_count_effect` (the discriminator), "
       "`test_sigma_only_reproduces_paper12` | **BACKED-SOUND** | new 2026-09-14; "
       "headline assertion added 2026-09-14 after the DELTA run found the matrix "
       "claimed BACKED-SOUND while only the (2,2) row was pinned.")
if old not in mat:
    print("  ! matrix row not matched -- check manually")
    sys.exit(1)
mat = mat.replace(old, new, 1)
with io.open(MATRIX, "w", encoding="utf-8") as fh:
    fh.write(mat)
print("  + claim_test_matrix: row corrected to name what is actually pinned")
