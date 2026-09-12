"""A seventh locus of the ill-conditioned cluster, and a test NAMED for the zombie.

Found while wiring the new backing into docs/claim_test_matrix.md.  Row 525
describes eq:pw's claim as

    "... while the L2 overlap of the same shared-scale basis is ill-conditioned
     (the lambda-blowup source)"

which is both halves of what Paper 60 Sec.2 withdrew:  the characterization
(cond(S) reaching 32 is not ill-conditioning) AND the mechanism attribution (the
inflation is the DENSITY of S^-1/2, not conditioning).  Its backing test is
called ``test_paper60_l2_overlap_illconditioned_grows`` and its docstring says
the overlap "ill-conditions with basis size".

The test's CONTENT is sound and stays exactly as it is:  monotone cond(S) growth
matching 3.0/5.8/13.9/32.2, and shared-scale lambda inflating faster than
hydrogenic.  Both are true and both survive the withdrawal.  What is wrong is the
NAME and the framing -- this is the "test that ASSERTS the zombie keyword"
sub-flavor the corpus catalogued on 2026-06-22 (v4.43.5), where two tests
asserted `"Latremoliere" in s` for a withdrawn attribution.

C16 did not catch it because `p60-l2-metric-diverges` keys on the 4->3673 forms
and "ill-conditioned by construction", not on the bare adjective applied to the
overlap.  Pattern widened, and re-tested in both directions afterwards -- the
paper's own denial ("ordinary, not ill-conditioned") must stay silent.

No paper cites this test, so the rename is safe for C13.
"""
import io

# ---- 1. rename + reframe the test (content untouched).
T = "tests/test_paper60_sturmian.py"
s = io.open(T, encoding="utf-8").read()
s = s.replace("def test_paper60_l2_overlap_illconditioned_grows():",
              "def test_paper60_l2_overlap_condition_number_grows():", 1)
OLD_DOC = ('''    """eq:blowup MECHANISM: the shared-scale Sturmian overlap ill-conditions with basis size,
    so Loewdin inflates the block-encoding 1-norm FASTER for the shared-scale basis than for a
    well-conditioned hydrogenic one.''')
NEW_DOC = ('''    """eq:blowup MECHANISM: the shared-scale Sturmian overlap becomes progressively
    LESS ORTHOGONAL with basis size, so Loewdin's dense S^-1/2 inflates the block-encoding
    1-norm FASTER for the shared-scale basis than for an orthonormal hydrogenic one.

    NOT an ill-conditioning claim -- renamed and reframed 2026-09-11.  cond(S) here
    reaches 32 at N=8, which is an ordinary Gram matrix;  Paper 60 Sec.2 withdrew the
    reading that the L2 metric is ill-conditioned (the 4 -> 3673 divergence was a
    radial-box artifact) and the inflation is driven by the DENSITY of S^-1/2, not by
    numerical instability.  The assertions below are unchanged and were always about
    the growth of cond(S) and the lambda-growth ORDERING, both of which survive.''')
assert OLD_DOC in s, "docstring locus not found"
s = s.replace(OLD_DOC, NEW_DOC, 1)
io.open(T, "w", encoding="utf-8").write(s)
print("test renamed: ..._illconditioned_grows -> ..._condition_number_grows; docstring reframed")

# ---- 2. the claim-matrix row wording + both test references.
M = "docs/claim_test_matrix.md"
m = io.open(M, encoding="utf-8").read()
m = m.replace(
    "while the L² overlap of the same shared-scale basis is ill-conditioned (the λ-blowup source)",
    "while the L² overlap of the same shared-scale basis is progressively NON-ORTHOGONAL "
    "(cond(S) 3.0→32.2 over N=2..8 — ordinary, not ill-conditioned; the λ-blowup source is "
    "the DENSITY of Löwdin's S^-1/2, corrected 2026-09-11)", 1)
m = m.replace("test_paper60_l2_overlap_illconditioned_grows",
              "test_paper60_l2_overlap_condition_number_grows")
io.open(M, "w", encoding="utf-8").write(m)
print("claim matrix: row 525 reworded; both test references renamed")

# ---- 3. widen C16 so the bare adjective applied to the overlap is caught.
C = "debug/qa/check_retracted_terms.py"
c = io.open(C, encoding="utf-8").read()
OLD_PAT = ('                   r"|L[\\u00b22]-divergence"\n'
           '                   r"|only well-conditioned posing",\n')
NEW_PAT = ('                   r"|L[\\u00b22]-divergence"\n'
           '                   r"|only well-conditioned posing"\n'
           '                   # added 2026-09-11: the bare adjective applied to the\n'
           '                   # overlap, and the verb form.  Deliberately anchored on\n'
           '                   # "overlap"/"metric" so the paper\'s own DENIAL ("ordinary,\n'
           '                   # not ill-conditioned") stays silent -- the denial does not\n'
           '                   # put the noun within reach of the adjective.\n'
           '                   r"|overlap[^.\\n]{0,30}is ill-conditioned"\n'
           '                   r"|ill-conditions with basis"\n'
           '                   r"|illconditioned",\n')
assert OLD_PAT in c, "C16 pattern locus not found"
c = c.replace(OLD_PAT, NEW_PAT, 1)
io.open(C, "w", encoding="utf-8").write(c)
print("C16 p60-l2-metric-diverges: pattern widened to the bare-adjective form")
