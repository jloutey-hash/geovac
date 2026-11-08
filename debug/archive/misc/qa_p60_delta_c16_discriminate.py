"""STEP 2: the discrimination test, and the two pattern repairs it forced.

The REGISTRY DISCRIMINATION RULE says a new entry must be proven to FIRE on the
retired wording and stay SILENT on the corrected wording.  Run against live text,
the two new entries fired correctly on the zombies -- including the LARGE at
CLAUDE.md:119, which no gate could previously reach -- but ALSO on two loci that
say the opposite of what the entry retires:

  papers/.../paper_60...tex:438   "Neither $T'$ grouping is superlinear"
  papers/synthesis/group2...tex:708  "neither grouping of the pure-number block
                                      $T'$ is superlinear ($K^{0.88}$ full ...)"

Both are the CORRECTED text.  A guard that fires on the correction is worse than
no guard:  it makes the right answer unwritable, and the next author's cheapest
escape is to reword the correction until the gate goes quiet.  `is\\s+SUPERlinear`
is simply the wrong shape -- an assertion and its denial share the verb.  The
numeric anchors (K^1.05, K^0.70, "local slope 0.906") carry the discrimination on
their own, because a denial names the value it is denying only in a withdrawal
context, which the standardized marker then covers.

  docs/claim_test_matrix.md:281,292,295   "l_max-divergence"

matched `L\\S{0,4}-divergence`, which was meant for code_architecture.md's
"L2-divergence".  `\\S{0,4}` happily eats "_max".  Narrowed to the two spellings
that actually occur.

This is the second time in this run that an over-broad alternative has been the
defect rather than a missing one -- the first was the inherited
`exempt_if_nearby` vocabulary.  Both directions of a pattern are load-bearing.
"""
import io

P = "debug/qa/check_retracted_terms.py"
s = io.open(P, encoding="utf-8").read()

# -- repair 1: drop the assertion/denial-ambiguous alternative.
OLD1 = '        "pattern": r"is\\s+SUPERlinear"\n                   r"|is the sublinear block"\n'
NEW1 = ('        # `is\\s+SUPERlinear` REMOVED 2026-09-11 after fire-testing:  it fired on\n'
        '        # the CORRECTED text ("Neither $T\'$ grouping is superlinear") at two\n'
        '        # loci.  An assertion and its denial share the verb, so the verb cannot\n'
        '        # discriminate;  the retired VALUES can, and do.\n'
        '        "pattern": r"is the sublinear block"\n')
assert OLD1 in s, "repair 1 locus not found"
s = s.replace(OLD1, NEW1, 1)

# -- repair 2: narrow the divergence spelling so it cannot eat "l_max-divergence".
OLD2 = '                   r"|L\\S{0,4}-divergence"\n'
NEW2 = '                   r"|L[\\u00b22]-divergence"\n'
assert OLD2 in s, "repair 2 locus not found"
s = s.replace(OLD2, NEW2, 1)

io.open(P, "w", encoding="utf-8").write(s)
print("patterns repaired:")
print("  p60-tprime-superlinear : dropped `is SUPERlinear` (fired on the denial)")
print("  p60-l2-metric-diverges : L\\S{0,4}-divergence -> L[²2]-divergence")
