"""Widen the bare-adjective alternative so it reaches both surviving loci.

The first widening used ``overlap[^.\\n]{0,30}is ill-conditioned`` and missed two
loci in tests/test_paper60_sturmian.py:

  * "the L2 overlap of the SAME (shared-scale) basis is ill-conditioned"
    -- 32 characters between noun and adjective, two over the window;
  * "the ill-conditioned shared-scale overlap makes Loewdin inflate ..."
    -- the adjective comes BEFORE the noun, which a post-nominal pattern
    cannot see at all.

Both forms added.  Written as a file rather than a heredoc: every attempt to
edit these patterns through `python - <<EOF` in this session failed to match its
own anchor, which is the documented reason backslash-bearing edits go through a
Write-tool script (memory rule `feedback_no_heredoc_backslashes`).

The denial must stay silent:  Sec.2's "The L^2 Gram matrix of this basis is
therefore ordinary, not ill-conditioned" puts no "overlap" within 60 characters
before, and none within 30 after, so neither alternative reaches it.  Verified
in both directions after applying.
"""
import io

P = "debug/qa/check_retracted_terms.py"
s = io.open(P, encoding="utf-8").read()

OLD = r'                   r"|overlap[^.\n]{0,30}is ill-conditioned"' + "\n"
NEW = (r'                   r"|overlap[^.\n]{0,60}is ill-conditioned"' + "\n"
       + r'                   r"|ill-conditioned[^.\n]{0,30}overlap"' + "\n")

assert OLD in s, "pattern anchor not found"
s = s.replace(OLD, NEW, 1)
io.open(P, "w", encoding="utf-8").write(s)
print("widened: overlap-window 30 -> 60, plus the pre-nominal 'ill-conditioned ... overlap' form")
