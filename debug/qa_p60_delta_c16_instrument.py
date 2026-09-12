"""/qa paper_60 DELTA remediation, STEP 1: the instrument, before any content fix.

The DELTA run returned 2 LARGE + 25 MATERIAL-SMALL.  Three of those findings are
not about the corpus at all -- they are about the gate:

  I1  No C16 entry lists CLAUDE.md, docs/claims_register.md,
      docs/code_architecture.md or geovac/sturmian_variational.py in `files`,
      so no gate could reach the LARGE (CLAUDE.md:119) or the two zombie
      docstrings in tracked code.

  I2  `p60-sublinear-as-regime` states its own CORRECTION in retired values:
      it says the fix is "T^0 at K^0.70 (stable), T' at K^1.05 (SUPERlinear)",
      both retired 2026-09-08 (now: T^0 has no power law, asymptotic 1/2;
      T' is 0.88, SUBlinear).  A gate whose correction text is wrong teaches
      the wrong answer to every reader who trusts it.

  I3  The same entry's `exempt_if_nearby` contains `nuclear diagonal|T\\^0` --
      the exact vocabulary of the retired mechanism.  Any locus asserting
      "the nuclear diagonal T^0 is the sublinear block" was exempted BY ITS
      OWN WRONG MECHANISM.  This is the authored-exemption failure declared
      fixed-as-a-class on 2026-09-08, recurring.

Design note on exemptions.  New entries here use `exempt_if_nearby = r"(?!)"`
(never matches) and rely solely on the standardized per-entry marker
`[retracted YYYY-MM-DD: <entry-id>]`.  The +-5-line legacy vocabulary CANNOT be
used for this class:  CLAUDE.md:119 sits two lines below a bullet containing
both "WITHDRAWN" and "2026-09-08", so any plausible retirement vocabulary would
have exempted the very LARGE this entry exists to catch.  An entry-specific
marker on the hit's own line is the only exemption that cannot be bought by an
accident of neighbouring text.

Per Sec.9 this is written and fire-tested as its OWN pass, before the content
remediation it protects -- not in the same pass.
"""
import io

P = "debug/qa/check_retracted_terms.py"
s = io.open(P, encoding="utf-8").read()

# ---------------------------------------------------------------- I2: the
# correction text of p60-sublinear-as-regime is itself retired.
OLD_RETIRED = (
    '                   "slope rises monotonically outside it (0.850, 0.868, 0.882, "\n'
    '                   "0.906 at K = 202, 244, 290, 340), and the split is the "\n'
    '                   "opposite of the printed mechanism -- nuclear diagonal T^0 = "\n'
    '                   "Z R_nu at K^0.70 (stable), pure-number block T\' at K^1.05 "\n'
    '                   "(SUPERlinear).  The paper\'s own Sec. 7 already said the "\n'
)
NEW_RETIRED = (
    '                   "slope was then thought to rise outside it.  CORRECTED AGAIN "\n'
    '                   "2026-09-08 on a converged domain -- see entry "\n'
    '                   "p60-tprime-superlinear:  the local slope FALLS (to 0.766 by "\n'
    '                   "K = 514), the window exponent is 0.82, the nuclear diagonal "\n'
    '                   "T^0 is not a power law at all (asymptotic 1/2 + O(1/log K)), "\n'
    '                   "and T\' is 0.88 -- SUBlinear, not superlinear.  The 2026-09-07 "\n'
    '                   "values K^0.70 / K^1.05 / rising-slope are THEMSELVES RETIRED "\n'
    '                   "and are guarded by that entry.  The paper\'s own Sec. 7 already said the "\n'
)
assert OLD_RETIRED in s, "I2: retired-prose locus not found"
s = s.replace(OLD_RETIRED, NEW_RETIRED, 1)

# ---------------------------------------------------------------- I3: the
# exemption vocabulary matches the retired mechanism.  Drop `nuclear
# diagonal` and `T^0` -- a locus must not be exempted by asserting the very
# thing the entry retires.
OLD_EXEMPT = (
    '        "exempt_if_nearby": r"K\\s*=\\s*74|window|over the computed|more slowly than "\n'
    '                            r"the (?:matrix dimension|configuration count)|nuclear "\n'
    '                            r"diagonal|T\\^0|NOT an asymptote|corrected 2026-09-07",\n'
)
NEW_EXEMPT = (
    '        # `nuclear diagonal` and `T^0` REMOVED 2026-09-11 (DELTA I3):  they are\n'
    '        # the vocabulary of the RETIRED mechanism, so a locus asserting "the\n'
    '        # nuclear diagonal T^0 is the sublinear block" exempted itself.  The\n'
    '        # surviving vocabulary names only the WINDOW scoping, which is the\n'
    '        # thing a corrected locus actually has to carry.\n'
    '        "exempt_if_nearby": r"K\\s*=\\s*74|window|over the computed|more slowly than "\n'
    '                            r"the (?:matrix dimension|configuration count)"\n'
    '                            r"|NOT an asymptote|corrected 2026-09-0[78]",\n'
)
assert OLD_EXEMPT in s, "I3: exemption locus not found"
s = s.replace(OLD_EXEMPT, NEW_EXEMPT, 1)

# ---------------------------------------------------------------- I1: widen
# `files` so the entry can reach the high-traffic live surfaces.
OLD_FILES = (
    '        "files": [\n'
    '            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",\n'
    '            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",\n'
    '            "papers/INDEX.md",\n'
    '            "docs/claim_test_matrix.md",\n'
    '            "tests/test_paper60_sturmian.py",\n'
    '            "geovac/sturmian_secular.py",\n'
    '            "geovac/sturmian_molecular_lambda.py",\n'
    '        ],\n'
)
NEW_FILES = (
    '        # Widened 2026-09-11 (DELTA I1).  The prior list was papers+docs only,\n'
    '        # so the two highest-traffic surfaces in the corpus -- CLAUDE.md, read\n'
    '        # by every session and every subagent dispatch, and the tracked geovac/\n'
    '        # docstrings -- were unreachable by this gate.\n'
    '        "files": [\n'
    '            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",\n'
    '            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",\n'
    '            "papers/INDEX.md",\n'
    '            "docs/claim_test_matrix.md",\n'
    '            "docs/claims_register.md",\n'
    '            "docs/code_architecture.md",\n'
    '            "CLAUDE.md",\n'
    '            "tests/test_paper60_sturmian.py",\n'
    '            "geovac/sturmian_*.py",\n'
    '        ],\n'
)
assert OLD_FILES in s, "I1: files locus not found"
s = s.replace(OLD_FILES, NEW_FILES, 1)

# ---------------------------------------------------------------- stale values
# inside p60-floor-is-angular's own note (4.49x / 1.12x are both retired;
# current DERIVED ratios are 6.8196/1.5936 = 4.28 and 1.7163/1.5936 = 1.08).
OLD_NOTE = (
    '                "GROUND-STATE specific: 4.49x chemical accuracy for the ground "\n'
    '                "state vs 1.12x for 2^1S at identical ||M||_1.  The genuine He "\n'
)
NEW_NOTE = (
    '                "GROUND-STATE specific: 4.28x chemical accuracy for the ground "\n'
    '                "state vs 1.08x for 2^1S at identical ||M||_1 (values corrected "\n'
    '                "2026-09-11 from the retired 4.49x/1.12x; both are DERIVED in "\n'
    '                "the numeric registry as p60_gnd_ratio_k452 / p60_exc_ratio_k452 "\n'
    '                "so they cannot drift from their energies again).  The genuine He "\n'
)
assert OLD_NOTE in s, "floor-is-angular note locus not found"
s = s.replace(OLD_NOTE, NEW_NOTE, 1)

# ---------------------------------------------------------------- NEW ENTRIES
ANCHOR = '    {\n        "id": "p60-sublinear-as-regime",'
assert ANCHOR in s, "anchor for new entries not found"

NEW_ENTRIES = '''    {
        "id": "p60-tprime-superlinear",
        "scope": "paper_60 group2 synthesis trunk",
        "severity": "fail",
        "retired": "2026-09-08 (v5.10.12, converged-domain re-measure).  The "
                   "2026-09-07 remediation of `p60-sublinear-as-regime` replaced one "
                   "wrong mechanism with another.  It asserted:  the nuclear diagonal "
                   "T^0 is the sublinear block at K^0.70; the pure-number block T' is "
                   "SUPERlinear at K^1.05; and the total local slope RISES with K "
                   "(0.850, 0.868, 0.882, 0.906 at K = 202..340).  All three are "
                   "retired.  Measured on a converged domain:  T^0 is NOT a power law "
                   "-- its exact closed form gives slope 1/2 + O(1/log K), running "
                   "0.709 at K=20 down to 0.603 at K=498004, still falling;  T' full "
                   "is K^0.88, SUBlinear;  and the TOTAL local slope FALLS, reaching "
                   "0.766 by K=514 against the 0.82 window fit.  The 0.906 figure was "
                   "the last point of a rising WINDOW artifact on the unconverged box. "
                   "WHAT SURVIVES: the 1-norm still grows more slowly than the matrix "
                   "dimension over every computable basis -- the encoding claim is "
                   "untouched;  only the block attribution and the direction of the "
                   "slope were wrong, twice.",
        "pattern": r"is\\s+SUPERlinear"
                   r"|is the sublinear block"
                   r"|K\\^\\{?0\\.70\\}?"
                   r"|K\\^\\{?1\\.05\\}?"
                   r"|local slope 0\\.906"
                   r"|slope rises monotonically"
                   r"|exponent trends toward 1",
        # NEVER-MATCH by design.  See the module docstring of the applying
        # script:  CLAUDE.md:119 sits two lines from a bullet containing both
        # "WITHDRAWN" and "2026-09-08", so ANY retirement vocabulary in a
        # +-5-line window would have exempted the LARGE this entry exists to
        # catch.  A locus that legitimately names these values must carry the
        # standardized per-entry marker on its own line.
        "exempt_if_nearby": r"(?!)",
        "cited_by": {
            "CLAUDE.md":
                "reviewed 2026-09-11 -- Sec.2 v5.10.10 bullet asserted all four "
                "retired values with no supersession marker; replaced per Sec.13.11 "
                "rule 9 (status updates replace, never append)",
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex":
                "reviewed 2026-09-11 -- Sec.4 carries the converged closed form "
                "(eq:T0_closed) and the falling-slope table",
            "docs/claims_register.md":
                "reviewed 2026-09-11 -- row 27 marks the 2026-09-07 values retired",
            "docs/claim_test_matrix.md":
                "reviewed 2026-09-11 -- rows carry the converged exponents",
        },
        "files": [
            "CLAUDE.md",
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "docs/code_architecture.md",
            "docs/topic_to_paper_lookup.md",
            "geovac/sturmian_*.py",
            "tests/test_paper60_*.py",
            "tests/test_sturmian_secular.py",
        ],
    },
    {
        "id": "p60-l2-metric-diverges",
        "scope": "paper_60 group2 synthesis",
        "severity": "fail",
        "retired": "2026-09-07/09-08.  Paper 60 Sec.2 claimed the L^2 Gram matrix of "
                   "the shared-scale Coulomb-Sturmian basis DIVERGES -- cond(S) 4 -> "
                   "3673 -- and that the metric-free posing is therefore the only "
                   "well-conditioned one.  Both are withdrawn.  The blow-up was a pure "
                   "radial-box artifact, switching on exactly where n_max^2 first "
                   "exceeds the domain size;  on a converged domain cond(S) = 16.0 at "
                   "K=100 and grows mildly and smoothly, approx 0.12*K (4.07 at K=10, "
                   "56.3 at K=452) -- an ORDINARY Gram matrix.  The paper's surviving "
                   "claim is a COST statement (the L^2 posing is expensive to encode), "
                   "which is a different assertion from numerical instability.  This "
                   "entry exists because the withdrawal reached the papers but not the "
                   "tracked code:  `solve_with_metric`'s docstring still said cond(S) "
                   "'grows into the thousands', and `_whiten` in the module promoted "
                   "the same week said the metric is 'ill-conditioned by construction, "
                   "so this truncation is required, not cosmetic' -- while the "
                   "truncation provably never fires (0 directions dropped at every "
                   "measured case).",
        "pattern": r"grows into the thousands"
                   r"|ill-conditioned by construction"
                   r"|4\\s*(?:->|-->|\\\\to)\\s*367[23]"
                   r"|cond\\(S\\)[^.\\n]{0,30}367[23]"
                   r"|L\\S{0,4}-divergence"
                   r"|only well-conditioned posing",
        "exempt_if_nearby": r"(?!)",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex":
                "reviewed 2026-09-11 -- Sec.2 carries the withdrawal and the "
                "converged cond(S) numbers",
            "docs/claim_test_matrix.md":
                "reviewed 2026-09-11 -- row 60/Sec.obstruction marked WITHDRAWN",
            "docs/claims_register.md":
                "reviewed 2026-09-11 -- row 26 carries the withdrawal",
            "tests/test_sturmian_secular.py":
                "reviewed 2026-09-11 -- test_L2_metric_conditioning_is_box_dependent "
                "pins the artifact AS an artifact",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
            "docs/claims_register.md",
            "docs/code_architecture.md",
            "CLAUDE.md",
            "geovac/sturmian_*.py",
            "tests/test_sturmian_secular.py",
            "tests/test_paper60_*.py",
        ],
    },
'''

s = s.replace(ANCHOR, NEW_ENTRIES + ANCHOR, 1)

io.open(P, "w", encoding="utf-8").write(s)
print("C16 instrument updated:")
print("  I2  p60-sublinear-as-regime: correction text re-priced to the converged values")
print("  I3  p60-sublinear-as-regime: T^0/nuclear-diagonal removed from exempt_if_nearby")
print("  I1  p60-sublinear-as-regime: files widened to CLAUDE.md + docs + geovac/sturmian_*.py")
print("      p60-floor-is-angular: note 4.49x/1.12x -> 4.28x/1.08x")
print("  NEW p60-tprime-superlinear      (the 2026-09-08 retirement class; reaches CLAUDE.md)")
print("  NEW p60-l2-metric-diverges      (the cond(S) zombie; reaches geovac/)")
