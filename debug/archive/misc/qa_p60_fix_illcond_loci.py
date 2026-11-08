"""The three loci the widened pattern now reaches, in tests/test_paper60_sturmian.py.

Two are genuine zombies the earlier claim-wide sweep missed because C16 keyed on
the 4->3673 forms rather than the bare adjective;  the third is my own reframing
docstring, which names the retired value legitimately and takes the marker.

In every case the ASSERTIONS are untouched and were always true:  cond(S) does
grow 3.0 -> 32.2 over N=2..8, and the shared-scale Loewdin 1-norm does inflate
faster than the hydrogenic one.  What changes is the word for cond ~ 32, which
Paper 60 Sec.2 settled: ordinary, not ill-conditioned.
"""
import io

MARK = "[retracted 2026-09-07: p60-l2-metric-diverges]"
P = "tests/test_paper60_sturmian.py"
s = io.open(P, encoding="utf-8").read()

pairs = [
    # L64 -- a comment on a live assertion
    ("    # the L2 overlap of the SAME (shared-scale) basis is ill-conditioned\n",
     "    # the L2 overlap of the SAME (shared-scale) basis is NON-ORTHOGONAL;\n"
     "    # cond > 10 here is an ordinary Gram matrix, not ill-conditioning\n"),
    # L83 -- the helper docstring, adjective before the noun
    ("    the ill-conditioned shared-scale overlap makes Loewdin inflate lambda faster than the\n"
     "    well-conditioned hydrogenic one.",
     "    the progressively NON-ORTHOGONAL shared-scale overlap makes Loewdin's DENSE\n"
     "    S^-1/2 inflate lambda faster than for the orthonormal hydrogenic basis."),
    # L123 -- my own reframing text, naming the retired reading legitimately
    ("reading that the L2 metric is ill-conditioned (the 4 -> 3673 divergence was a",
     "reading that the L2 metric is ill-conditioned " + MARK + " (the 4 -> 3673 divergence was a"),
    # tidy the seam the earlier docstring merge left
    ("the growth of cond(S) and the lambda-growth ORDERING, both of which survive.  Pins (i) monotone cond(S) growth over >=4 points matching",
     "the growth of cond(S) and the lambda-growth ORDERING, both of which survive.\n\n"
     "    Pins (i) monotone cond(S) growth over >=4 points matching"),
]
for old, new in pairs:
    assert old in s, "not found: %.60s" % old.strip()
    s = s.replace(old, new, 1)
io.open(P, "w", encoding="utf-8").write(s)
print("fixed 2 zombie loci, marked 1 legitimate mention, tidied the docstring seam")
