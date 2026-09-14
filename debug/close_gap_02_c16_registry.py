"""KNOWN-GAP 1/5, step 2 -- the six C16 entries for the 2026-09-12 withdrawals.

Each entry is proven to discriminate TWO WAYS by close_gap_03_prove_c16.py
before it is trusted (the REGISTRY DISCRIMINATION RULE): it must FIRE on the
retired wording and stay SILENT on the corrected wording.  An entry that fires
on nothing is worse than no entry, because the gate then reports PASS and the
class is believed guarded.

Only `p60-removability-corollary` needs an exemption, because the corrected
text quotes the retired wording in order to deny it; it uses the standardized
marker and nothing else.  The other five use (?!) and never exempt anything --
verified: every one of their retired phrases is already at 0 occurrences.

Idempotent.
"""
from __future__ import annotations

import sys

R = "debug/qa/check_retracted_terms.py"
MARKER = '"id": "p60-removability-corollary"'
ANCHOR = '    {\n        "id": "p60-l2-metric-diverges",\n'

REVIEWED = "reviewed 2026-09-12 -- /qa paper_60 FULL remediation"
PAPER = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
MTX = "docs/claim_test_matrix.md"

FILES = f'''        "files": [
            "{PAPER}",
            "{SYN}",
            "{MTX}",
            "CLAUDE.md",
            "geovac/sturmian_*.py",
            "tests/test_paper60_*.py",
        ],
'''

ENTRIES = f'''    {{
        "id": "p60-removability-corollary",
        "scope": "paper_60 group2 synthesis",
        "severity": "fail",
        "retired": "2026-09-12, /qa paper_60 FULL.  The transcendental tagging "
                   "paragraph drew an operational corollary -- that a "
                   "TRUNCATION-side price is 'a property of the matrix, which a "
                   "preconditioner reaches' while a CONTINUUM-side price is 'a "
                   "property of the symbol, which no congruence can touch'.  BOTH "
                   "HALVES ARE FALSE.  The preconditioner is built FROM the symbol "
                   "(its matching polynomial is chosen to share the symbol's zero, "
                   "per Serra), and preconditioning IS a congruence of the finite "
                   "section, one that replaces the symbol by f/g.  The surviving "
                   "statement concerns the symbol on BOTH sides and turns on the "
                   "KIND of feature: a banded congruence cancels a zero of finite "
                   "order but cannot alter a decay class.  Provenance of a constant "
                   "predicts nothing about removability.",
        "pattern": r"truncation-side price"
                   r"|continuum-side price"
                   r"|property of the matrix, which a preconditioner reaches"
                   r"|property of the symbol, which no congruence",
        "exempt_if_nearby": r"\\[retracted \\d{{4}}-\\d\\d-\\d\\d: p60-removability-corollary\\]",
{FILES}        "cited_by": {{
            "{PAPER}": "{REVIEWED} -- the paragraph carries the withdrawal inline",
            "{MTX}": "{REVIEWED} -- the tagging row records the corollary as WITHDRAWN",
        }},
    }},
    {{
        "id": "p60-frames-completeness",
        "scope": "paper_60 group2 synthesis",
        "severity": "fail",
        "retired": "2026-09-12 (adopted 2026-09-11, withdrawn 2026-09-12 v5.11.2).  "
                   "The frames reading -- that overcompleteness is the PRICE of "
                   "one-centre completeness, i.e. that completeness of the "
                   "one-centre set alone forces lam_min -> 0 -- is FALSE for this "
                   "basis.  Measured: the Bessel deficit of a displaced Sturmian "
                   "against the one-centre span PLATEAUS at 0.380/0.696/0.907 for "
                   "kR = 1/2/4, flat over N = 16..256, so the one-centre set is "
                   "measurably far from complete IN THE MOLECULAR METRIC and the "
                   "bound holds only vacuously.  Ron-Shen is the surviving "
                   "mechanism, and it says something narrower and more useful: the "
                   "near-dependence is ONE DIRECTION.",
        "pattern": r"overcompleteness is the price"
                   r"|price of one-cent(?:er|re) completeness"
                   r"|completeness of the one-cent(?:er|re) set[^.\\n]{{0,40}}forces",
        "exempt_if_nearby": r"(?!)",
{FILES}        "cited_by": {{
            "{PAPER}": "{REVIEWED} -- sec:molecular carries the plateau measurement",
            "{MTX}": "{REVIEWED} -- the one-direction row records the withdrawal",
        }},
    }},
    {{
        "id": "p60-sigma-law-derived",
        "scope": "paper_60 group2 synthesis",
        "severity": "fail",
        "retired": "2026-09-11, re-swept 2026-09-12.  eq:sigma_law is NOT derived "
                   "in this corpus: it is the Kac-Murdock-Szego extreme-eigenvalue "
                   "asymptotic (c_1 = pi^2, 1953).  What is ours is the "
                   "IDENTIFICATION of the two-centre Shibuya-Wulfman metric as such "
                   "a finite section.  The 2026-09-11 re-attribution reached the "
                   "paper body and MISSED the abstract, the group2 synthesis, the "
                   "sturmian_sigma_law module docstring and its backing test's "
                   "docstring -- all four fixed 2026-09-12.",
        "pattern": r"growth law is derived"
                   r"|derived band-limited law"
                   r"|[Bb]acks the derived conditioning law"
                   r"|conditioning law[^.\\n]{{0,24}}derived here",
        "exempt_if_nearby": r"(?!)",
{FILES}        "cited_by": {{
            "{PAPER}": "{REVIEWED} -- abstract now carries [PRIOR ART]",
            "{SYN}": "{REVIEWED} -- the molecular paragraph credits KMS",
            "{MTX}": "{REVIEWED} -- the KMS row states identification-only",
        }},
    }},
    {{
        "id": "p60-prop-d-as-new",
        "scope": "paper_60 group2 synthesis",
        "severity": "fail",
        "retired": "2026-09-12, C23 run #1.  'Proposition D' was demoted: the "
                   "block-diagonal-congruence result is Loewdin symmetry "
                   "preservation specialised to the l grading, known in this "
                   "paper's own field since Slater-Koster (1954), and in operator "
                   "terms the statement that block-diagonal matrices are a "
                   "commutant and therefore inverse-closed.  What the paper claims "
                   "is the l-vs-m APPLICATION, not a new proposition.",
        "pattern": r"Proposition~?D\\b"
                   r"|our Proposition[^.\\n]{{0,30}}block-diagonal congruence",
        "exempt_if_nearby": r"(?!)",
{FILES}        "cited_by": {{
            "{PAPER}": "{REVIEWED} -- sec:molecular credits Loewdin/Slater-Koster",
            "{MTX}": "{REVIEWED} -- the row records the re-attribution",
        }},
    }},
    {{
        "id": "p60-translation-ours",
        "scope": "paper_60 group2 synthesis",
        "severity": "fail",
        "retired": "2026-09-12.  The reading of the Shibuya-Wulfman operator as a "
                   "TRANSLATION is prior art on three counts -- Shibuya and "
                   "Wulfman's own 1965 abstract builds the molecular p0 operator "
                   "from 'a sum of unitary transformations, one for each nucleus in "
                   "the molecule'; Wulfman and Takahata gave the explicit "
                   "continuous-group formulation in 1967; Red and Weatherford "
                   "derived the general formula for that matrix in a "
                   "Coulomb-Sturmian basis in 2004.  What survives as ours is the "
                   "SYMBOL -- that in the sine basis the operator is the finite "
                   "section of multiplication by j0(kR cot(chi/2)).",
        "pattern": r"equivalently that the Shibuya--Wulfman\\s+operator is multiplication"
                   r"|what is ours[^.\\n]{{0,60}}translation phase"
                   r"|the translation identification is ours",
        "exempt_if_nearby": r"(?!)",
{FILES}        "cited_by": {{
            "{PAPER}": "{REVIEWED} -- the concession names all three counts",
            "{MTX}": "{REVIEWED} -- the symbol-not-translation split recorded",
        }},
    }},
    {{
        "id": "p60-tau-membership",
        "scope": "paper_60 group2 synthesis",
        "severity": "fail",
        "retired": "2026-09-12, /qa paper_60 FULL (code dimension).  The paper "
                   "glossed the tau (DST-I) algebra as 'which is to say Toeplitz "
                   "minus Hankel, the structure of this section', surrendering the "
                   "weight-independence result to Ahmad et al.'s identity.  FORM IS "
                   "NOT MEMBERSHIP: tau requires the coefficient sequence to "
                   "TERMINATE, which the matching polynomial's does and the chirp "
                   "symbol's does not.  Measured: the DST-I leaves 4.1/2.0/1.1 "
                   "percent off-diagonal weight on the cross block at n = 16/64/160 "
                   "against 2e-14 for a genuinely-tau matrix.  The identity covers "
                   "the tau IDEALISATION; the departure is exactly the residue the "
                   "paper reports.",
        "pattern": r"tau \\(DST-I\\)\\s*\\n?\\s*algebra --- which is to say Toeplitz minus Hankel"
                   r"|\\\\tau\\$? \\(DST-I\\)[^.\\n]{{0,40}}structure of this\\s*\\n?section"
                   r"|our matrices are in the \\$?\\\\?tau",
        "exempt_if_nearby": r"(?!)",
{FILES}        "cited_by": {{
            "{PAPER}": "{REVIEWED} -- the paragraph states the membership failure",
            "{MTX}": "{REVIEWED} -- the hypothesis-gap clause records it",
        }},
    }},
'''


def main() -> int:
    with open(R, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(ANCHOR) != 1:
        print(f"anchor count={t.count(ANCHOR)}; ABORT")
        return 2
    with open(R, "w", encoding="utf-8") as fh:
        fh.write(t.replace(ANCHOR, ENTRIES + ANCHOR))
    print("applied: 6 C16 entries added (discrimination NOT yet proven)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
