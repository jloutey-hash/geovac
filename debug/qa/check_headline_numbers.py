"""C17 -- headline-number registry gate (added 2026-07-02, PI direction).

The 7th group4 cert's meta-lesson: the MATERIAL classes that kept surviving
judgment review are MECHANICAL --
  (a) second-locus propagation: a decided headline value corrected at one
      locus and stale at another (Z=1-36 vs 56 across five loci; the
      memo-listed-but-never-applied KH fix), and
  (b) number-vs-source drift: a stated value contradicting its own cited
      source (the "190x" floor vs the cited table's 51x; the stale 33.3
      1-norm vs the live 32.6).
Both are registry-checkable. Each entry holds a headline FAMILY: either a
set of known-wrong variant patterns (C16 style) or a capture pattern plus
the CANONICAL value (any capture that disagrees is a live hit). Exempt
markers cover legitimately historical/withdrawn mentions.

MAINTENANCE RULE (mirrors C16): when a cert run corrects or demotes a
headline number, ADD/UPDATE its family here so the wrong value can never
silently re-surface at any locus.

Usage: python debug/qa/check_headline_numbers.py [--gate <branch>] [--all]
Exit 0 = PASS. Mirror test: tests/test_headline_numbers_check.py.
"""
from __future__ import annotations

import glob
import pathlib
import re
import sys
import os

# Shared --gate scope resolution (see debug/qa/qa_scopes.py): named
# scopes resolve to an explicit file list and every RESULT line carries
# the file count, so a gate can never report PASS on an empty scope.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qa_scopes  # noqa: E402


ROOT = pathlib.Path(__file__).resolve().parents[2]
WINDOW = 3  # +/- lines for the exemption window

GROUP4_FILES = [
    "papers/group4_quantum_computing/*.tex",
    "papers/synthesis/group4_quantum_computing_synthesis.tex",
]


# ---------------------------------------------------------------------------
# STANDARDIZED WITHDRAWAL MARKER (2026-09-03, /qa trunk FULL run #4 follow-up)
#
# `[retracted YYYY-MM-DD]` is accepted as a withdrawal flag for EVERY entry,
# globally, in addition to that entry's own `exempt_if_nearby`.
#
# Rationale, from measured defects rather than taste: three registry entries
# written this session exempted on vocabulary drawn from the surrounding
# CORRECTED text -- "misnomer", "corrected 2026-09-03",
# "bound|deficit|saturat" -- and correct text is exactly what surrounds a
# defect, so each reported clean on a live locus.  Authored exemption
# vocabulary is the failure mode; a fixed token removes the authoring step.
#
# New entries should set `exempt_if_nearby` to WITHDRAWAL_MARKER and nothing
# else.  Existing entries keep their lists (85+ loci depend on them) and the
# gate reports how many still do, so the debt shrinks instead of hiding.
# ---------------------------------------------------------------------------
# REDESIGNED 2026-09-04 (/qa DELTA #5).  The first version was a BARE
# token, OR'd into every entry's exemption -- so one
# `[retracted YYYY-MM-DD]` silenced ALL entries in its window, including
# entries written for unrelated claims.  A probe confirmed it: a comment
# withdrawing the Hopf-base label silenced a live volume-quotient defect
# in the next sentence.  A fixed token names nothing, which is exactly
# what the discrimination rule above forbids.
#
# The marker now carries the entry id it withdraws, so it exempts THAT
# entry and no other, while still requiring nothing to be invented --
# the id is looked up from the registry:
#
#     [retracted 2026-09-04: hopf-base-label-for-4-over-pi]
def withdrawal_marker(entry_id: str) -> str:
    """Regex matching the standardized marker FOR THIS ENTRY only."""
    return (r"\[retracted \d{4}-\d{2}-\d{2}:\s*"
            + re.escape(entry_id) + r"\]")


WITHDRAWAL_MARKER = "see withdrawal_marker()"  # sentinel, per-entry now

REGISTRY = [
    {
        "id": "trunk-saturation-rate-constant",
        "scope": "trunk group3",
        "severity": "fail",
        "canonical_note": "FULL run #4, 2026-09-03.  The l-block grid "
                          "structure makes the Laplacian spectrum closed "
                          "form, so the saturation deficit 2 d_max - "
                          "lambda_max has an EXACT asymptotic constant: "
                          "C = pi^2 (2 + 2^(1/3))^2 (1 + 2^(-2/3)) / 4 = "
                          "42.7397...  The printed 42.6 is the n_max = 320 "
                          "SAMPLE, not the limit (the sequence rises: 42.606 "
                          "at 320, 42.731 at 5000, 42.739 at 40000).  42.6 "
                          "stays legal as an explicitly-labelled finite-n "
                          "value; it is wrong as 'the rate' or with o(1).",
        "pattern": r"42\.6\s*\+\s*o\(1\)"
                   r"|\(\s*42\.6\s*\+\s*o\(1\)\s*\)"
                   r"|rate\s*\$?\(?42\.6"
                   r"|42\.6\s*(?:\+\s*o\(1\)\s*)?/\s*n",
        "require_nearby": r"n_\{?\\?max\}?|saturation|deficit|lambda_\{?\\?max|rate",
        "exempt_if_nearby": r"sample|finite-n|at 320|n_\{\\max\} = 320|"
                            r"corrected 2026-09|retired|superseded|42\.74",
        "files": [
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "tests/test_paper1_block_spectrum.py",
            "docs/claim_test_matrix.md",
        ],
    },
    {
        "id": "trunk-hydrogen-bound-not-accuracy",
        "scope": "trunk group3",
        "severity": "fail",
        "canonical_note": "The hydrogen graph number is a SATURATION DEFICIT "
                          "of a spectral bound, not an accuracy: E_0 = "
                          "kappa*lambda_max by construction, so the quantity "
                          "measured is lambda_max -> 2 d_max = 8.  CANONICAL: "
                          "0.574% at n_max = 30, 0.325% at 40, 0.107% at 70.  "
                          "RETIRED: '< 0.1%' as a hydrogen accuracy figure "
                          "(corrected 2026-09-03).",
        "pattern": r"H\s*\(hydrogen\)[^|\n]{0,40}\|\s*<\s*0\.1\\?%"
                   r"|hydrogen[^.\n]{0,60}accuracy[^.\n]{0,30}0\.1\\?%"
                   r"|\\?%\s*error[^.\n]{0,30}hydrogen graph",
        "require_nearby": r"hydrogen|graph|lambda|deficit|S\^?3",
        # NOTE (discrimination proof, 2026-09-03, third instance of
        # this error): this list held "bound|deficit|saturat|by
        # construction" -- the ordinary vocabulary of every
        # neighbouring benchmark row -- so a planted "< 0.1%" was
        # exempted by its own corrected siblings.  An exemption must
        # name the SPECIFIC retraction (a date, a retirement word, or
        # the replacement value), never the topic's vocabulary.
        "exempt_if_nearby": r"corrected 2026-09-03|retired|formerly|"
                            r"no longer|superseded",
        # benchmark TABLE: see the exempt_same_line note in scan_entry
        "exempt_same_line": True,
        "files": [
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "docs/validation_benchmarks.md",
            "CLAUDE.md",
        ],
    },
    {
        "id": "trunk-four-over-pi-normalisation",
        "scope": "trunk group1 group3 synthesis",
        "severity": "fail",
        "canonical_note": "DELTA #4 + FULL #4, 2026-09-03.  The state-space GH "
                          "rate constant is CONVENTION-DEPENDENT: 4/pi in the "
                          "corpus's declared dual-Coxeter rule Cas(ad) = h^v, "
                          "2/pi on the unit round S^3, 2 sqrt(2)/pi under the "
                          "field-standard Kac basic form (theta|theta) = 2.  "
                          "The unit-radius sphere-volume quotient is "
                          "Vol(S^2)/Vol(S^3) = 2/pi; 4/pi is TWICE that.  So "
                          "'4/pi is a quotient of unit-radius sphere volumes' "
                          "is FALSE as printed (Paper 32 cor:m1_pure_tate), "
                          "and any locus stating the constant must name its "
                          "normalisation.",
        # LaTeX markup may sit anywhere inside the phrase, so the words
        # are matched with markup-tolerant gaps rather than literally.
        "pattern": r"4\s*/\s*\\?pi[^.\n]{0,90}quotient\W{0,20}of\W{0,20}(?:\\emph\{)?unit-radius"
                   r"|quotient\W{0,20}of\W{0,20}(?:\\emph\{)?unit-radius[^.\n]{0,60}4\s*/\s*\\?pi"
                   r"|4\s*/\s*\\?pi\s*=\s*\\?(?:mathrm\{)?Vol\}?\(S\^\{?2\}?\)\s*/\s*\\?(?:mathrm\{)?Vol\}?\(\\?(?:s?three|S\^\{?3\}?)\)",
        "require_nearby": r"Vol|volume|sphere|quotient|normalis|rate",
        # NOTE (discrimination proof, 2026-09-03): "misnomer" and
        # "corrected 2026-09-03" were in this list and swallowed the
        # very locus the family exists for -- both appear in the
        # defective sentence, flagging the Hopf-label fix it carries
        # rather than the volume-quotient error it makes.  The only
        # honest exemption is evidence of THIS correction: the factor 2.
        "exempt_if_nearby": r"twice that|2\s*\\?(?:mathrm\{)?Vol|"
                            r"is\s+2\s*/\s*\\?pi|not a quotient",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "docs/qa/trunk.done.md",
            "docs/claims_register.md",
        ],
    },
    {
        "id": "trunk-forced-count-full-axiom",
        "scope": "trunk group1",
        "severity": "fail",
        "canonical_note": "The Forced-Count moduli chain endpoint is 32 "
                          "(2048 -> 1024 -> 512 -> 272 -> 32), stated at the "
                          "FULL-AXIOM count per the trunk C8 delta.  RETIRED: "
                          "260 (a representation bug, v5.4.1) and 128 (the "
                          "matter-sector subcount presented as the endpoint).",
        "pattern": r"full-axiom[^.\n]{0,60}\b(?:260|128)\b"
                   r"|\b(?:260|128)\b[^.\n]{0,40}full-axiom"
                   r"|forced[- ]count[^.\n]{0,40}\b(?:260|128)\b",
        "require_nearby": r"moduli|forced|axiom|chain|order-one",
        "exempt_if_nearby": r"retired|corrected|superseded|was |formerly|"
                            r"matter-sector subcount|no longer|2026-09",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group3_foundations/paper_57_forced_free_seam.tex",
            "docs/claim_test_matrix.md",
            "docs/qa/trunk.done.md",
        ],
    },
    {
        "id": "trunk-sp-splitting-residue-disclosure",
        "scope": "trunk group3 synthesis",
        "severity": "fail",
        "canonical_note": "FULL run #4, 2026-09-03.  The s/p splitting "
                          "endpoint 0.39% at n_max = 30 is SELECTION-BIASED: "
                          "lambda_2s = 3 EXACTLY whenever n_max = 0 (mod 3) "
                          "and only there, so the multiples of three lie on a "
                          "favourable branch (0.82/0.62/0.49/0.39% at "
                          "21/24/27/30) while their neighbours are 3-6x "
                          "larger (1.65% at 28, 2.58% at 29).  Any locus "
                          "quoting the endpoint must disclose the residue "
                          "class; quoting it as a convergence endpoint "
                          "without that disclosure is the defect.",
        "pattern": r"0\.39\d?\s*\\?%|0\.39\d?\\?%",
        "require_nearby": r"s/p|2s|2p|splitting|degenerac|n_\{?\\?max\}?\s*=\s*30",
        "exempt_if_nearby": r"mod 3|mod\\,3|residue|divisib|multiple of three|"
                            r"branch|selection-bias|favourable|2026-09-03",
        "files": [
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "docs/claim_test_matrix.md",
        ],
    },
    {
        "id": "p39-tensor-assembly-constant",
        "scope": "trunk group1",
        "severity": "fail",
        "canonical_note": "Trunk DELTA #6 remediation, 2026-09-04.  The "
                          "assembled tensor state-space GH constant is 2 "
                          "(reach-only: reach <= gamma_a + gamma_b <= 2 "
                          "max(gamma)), NOT 1 + 2 sqrt 2 ~ 3.828.  The "
                          "retired value added a Lipschitz-distortion height "
                          "leg (2 sqrt 2) that is refuted -- the joint height "
                          "is identically 1 by the finite-band witness "
                          "transported from Paper 38 -- and used max(gamma) "
                          "for the reach where L4(c) gives gamma_a + gamma_b. "
                          "The Lambda^full column printed 1.914x too large. "
                          "The cross-Stein-Weiss 2 sqrt 2 survives as a "
                          "unit-norm PANEL bound and is legal in that role.",
        "pattern": r"\(\s*1\s*\+\s*2\\sqrt\{2\}\s*\)\s*\\cdot"
                   r"|1\s*\+\s*2\s*\*?\s*sqrt\(2\)\s*\)?\s*(?:\*|max)"
                   r"|reach-plus-height"
                   r"|1\s*\+\s*2\\sqrt\{2\}\s*\\approx\s*3\.828",
        "require_nearby": r"max\(?\\?gamma|C_3|propinquity|Lambda|assembl|tensor",
        # DELTA #7 (CODE-B M2): "panel"/"formerly"/a cross-entry marker each
        # exempted the DEFECTIVE text -- the over-exemption class removed from
        # C16 the day before, reintroduced here in named form.  Own marker only.
        "exempt_if_nearby": r"(?!)",
        "files": [
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
            "geovac/gh_convergence_tensor.py",
            "tests/test_gh_convergence_tensor.py",
        ],
    },
    {
        "id": "composed-lih-market-test-retired",
        "scope": "all",
        "severity": "fail",
        "canonical_note": "Retired composed-LiH market-test figures and the "
                          "retired 'constant factor' framing.  CANONICAL "
                          "(exact rule): LiH composed 838 Pauli @ Q30 = near "
                          "parity with raw STO-3G 907, and 3.0x MORE than the "
                          "qubit-reduced 276; 1-norm 1.007x (34.5 vs 34.3 Ha). "
                          "RETIRED and now WRONG as live values: '334 Pauli' "
                          "(or 333) as the LiH composed count, the '13x fewer "
                          "QWC' market test, and the 'constant 2.51x "
                          "main-group / 3.25x d-block, scaling robust' "
                          "quantification -- the latter was a SINGLE-POINT "
                          "ARTIFACT: the re-pricing is not a constant factor "
                          "and every exponent moved.  See "
                          "debug/sprint_eri_evaluator_defects_memo.md + "
                          "debug/qa/delta4_run_notes.md.",
        "pattern": r"33[34]\s*~?\\?(Pauli|terms)|"
                   r"(constant|factor of)\s*\$?2\.51|"
                   r"2\.51\s*\\?times[^.]{0,40}(scaling robust|3\.25)|"
                   r"13\s*\\?times\s*fewer\s*QWC",
        "require_nearby": r"Pauli|LiH|composed|market|re-pric|scaling|QWC",
        "exempt_if_nearby": r"retired|corrected|Corrected 20|vintage|artifact|"
                            r"rule gave|previously|earlier|withdrawn|"
                            r"superseded|dissolve|DISSOLVED|former",
        "files": [
            "papers/group4_quantum_computing/paper_14_qubit_encoding.tex",
            "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex",
            "papers/synthesis/group4_quantum_computing_synthesis.tex",
            "docs/qa/group4.done.md",
        ],
    },
    {
        "id": "composed-retired-scaling-and-counts",
        "scope": "all",
        "severity": "fail",
        "canonical_note": "Second-locus forms of the retired pair-diagonal "
                          "rule that the composed-rule-a entry does not "
                          "match.  CANONICAL (exact rule, measured "
                          "2026-08-29/30): s+p block at n_max=2 = 279 "
                          "non-identity Pauli from 107 ERIs, 27.90 per qubit; "
                          "multi-center rows 1,953 (Q=70) / 2,790 (Q=100) / "
                          "1,395 (Q=50); term exponent Q^3.8; P20 "
                          "within-molecule two-point exponent 2.816 EXACTLY "
                          "for all six molecules (= 1 + log_5(18.6), forced by "
                          "exact linearity in Q), Q=100 Gaussian ratio ~370x.  "
                          "RETIRED and now WRONG as live values: 1,111 "
                          "multi-center rows, 111 Pauli per block, 65 ERIs "
                          "per block, 11.1 per qubit, Q^3.15, 51x-1712x "
                          "(any spacing), per-molecule 2.18-2.24, mean "
                          "2.21 +- 0.02, and the 6,000x extrapolation.  See "
                          "debug/qa/delta4_run_notes.md.",
        "pattern": r"1\{,\}111|1712\s*\\?times|Q\^\{?3\.15\}?|"
                   r"2\.21\s*\$?\\pm\$?\s*0\.02|6\{,\}000\s*\\?times|"
                   r"exponent\s+of\s*~?\$?\\sim\s*2\.2\$?",
        "require_nearby": r"Pauli|scaling|exponent|composed|GeoVac|term|ratio",
        "exempt_if_nearby": r"retired|corrected|Corrected 20|vintage|artifact|"
                            r"rule gave|previously|earlier|withdrawn|"
                            r"pair-diagonal",
        "files": [
            "papers/group4_quantum_computing/paper_14_qubit_encoding.tex",
            "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex",
            "papers/group6_precision_observations/paper_26_entanglement.tex",
            "papers/synthesis/group4_quantum_computing_synthesis.tex",
            "docs/qa/group4.done.md",
        ],
    },
    {
        "id": "composed-rule-a-retired-figures",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "Composed/atomic QC resource figures, CORRECTED "
                          "2026-08-29 (exact global-M_L rule; the pair-diagonal "
                          "'rule A' was a wrong-sign-q bug, CF-1 DISSOLVED).  "
                          "CANONICAL: N_Pauli = 27.90 x Q main-group / 30.03 x Q "
                          "d-block; within-molecule exponent 3.17 universal "
                          "(c(1)=3/2, c(2)=279/10, c(3)=7089/14; local slope "
                          "~3.8 at n_max=4); LiH composed 838 @ Q30, 1-norm "
                          "34.5 Ha; equal-qubit H2O 54x/297x/317x; cc-pVDZ "
                          "76x.  RETIRED and now WRONG as live values: "
                          "11.10 x Q, 9.23 d-block, O(Q^2.5) / Q^2.50-2.52 "
                          "exponents, 334/333 LiH, 51x-1712x, 190x cc-pVDZ, "
                          "32.6 Ha 1-norm, and the 'constant factor 2.51x/"
                          "3.25x scaling-unchanged' A->B claim.  See "
                          "debug/sprint_eri_evaluator_defects_memo.md; backed "
                          "by tests/test_paper14_eri_rule.py + "
                          "tests/test_paper20_balanced_lambda.py.",
        "pattern": r"11\.10\s*\\?times\s*Q|N_\{?\\?mathrm\{Pauli\}\}?\s*=\s*11\.10|"
                   r"O\(Q\^\{?2\.5\}?\)|Q\^\{2\.50\}|1\{,\}712\\?\$?\\times|"
                   r"190\$?\\times|coefficient\s+(of\s+)?11\.10|Pauli/\$?Q\$?\s*=\s*9\.23",
        "require_nearby": r"Pauli|scaling|coefficient|linear|composed|exponent",
        "exempt_if_nearby": r"retired|corrected|was measured|vintage|inverted|"
                            r"artifact|dissolve|earlier|Corrected 2026-08-29|"
                            r"pair-diagonal rule gave|rule gave|gave 11\.10|"
                            r"whose apparent",
        "files": [
            "papers/group4_quantum_computing/paper_14_qubit_encoding.tex",
            "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex",
            "papers/synthesis/group4_quantum_computing_synthesis.tex",
            "docs/qa/group4.done.md",
        ],
    },
    {
        "id": "p26-eri-sparsity-density",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Paper 26 Sec III basis-intrinsic ERI sparsity, "
                          "CORRECTED 2026-08-29.  The retired counts were taken "
                          "from an evaluator that omitted the Coulomb selection "
                          "rule m_a + m_b = m_c + m_d, so 59.6% of them were "
                          "physically zero.  CANONICAL: n_max=2 full-tensor "
                          "107/625 = 17.1%; canonical-unique 41/225 = 18.2%; "
                          "n_max=4 full-tensor 57,700/810,000 = 7.12%; "
                          "canonical-unique 15,293/216,225 = 7.07% "
                          "(Z-independent at Z = 2, 3, 10).  The density "
                          "IMPROVES with basis size (~M^-0.49).  RETIRED and now "
                          "WRONG: 265/625, 42.4%, 318,720, 39.3%, 79,465, 36.8%, "
                          "and the 'essentially flat / does not degrade with "
                          "basis size' sub-claim, which REVERSES.  See "
                          "debug/sprint_eri_evaluator_defects_memo.md; backed by "
                          "tests/test_paper26_entanglement.py.",
        "pattern": r"42\.4\\?%|265\s*/\s*625|265 nonzero|318\{,\}720|318720|"
                   r"79\{,\}465|79465|39\.3\\?%|36\.8\\?%",
        "require_nearby": r"ERI|sparsit|densit|nonzero|tensor",
        "exempt_if_nearby": r"retired|corrected|withdrawn|previously|earlier version|"
                            r"now wrong|Corrected 2026-08-29",
        "files": [
            "papers/group6_precision_observations/paper_26_entanglement.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
            "docs/qa/group6.done.md",
        ],
    },
    {
        "id": "p34-d-polarizability-sourced",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Deuteron 1S HFS nuclear structure, re-sourced 2026-08-29 "
                          "against Bonilla et al. arXiv:2508.18776 (electronic-D TPE). "
                          "CANONICAL: total TPE 44.5(1.1) kHz = +135.9 ppm of the "
                          "splitting; polarizability channel +110.16 kHz = +336.5 ppm; "
                          "all elastic-class pieces -200.6 ppm; a chain consuming the "
                          "total lands at +23.0 ppm (cf. H 21cm +18.4). RETIRED and now "
                          "WRONG: '+44 ppm deuteron polarizability' (a kHz-read-as-ppm "
                          "unit slip of the TOTAL) and the '~+200 ppm' polarizability "
                          "entry of the withdrawn PY-style budget. Backed by "
                          "tests/test_paper34_autopsy_baselines.py "
                          "(test_paper34_d_hfs_tpe_decomposition_and_sign_flip, "
                          "..._channel_target_and_h_parallel, ..._unit_slip_guard).",
        "pattern": r"\$\+44\$~?\s?ppm|\+44~ppm",
        "require_nearby": r"polariz|deuteron|QCD-internal",
        "exempt_if_nearby": r"retired|corrected|withdrawn|unit slip|previously|earlier version",
        "files": [
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "p27-ep2b-commutator-column",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Paper 27 tab:ep2b commutator column (and Paper 24 copy): "
                          "CANONICAL 0.74 / 0.63 / 0.67 at N_max = 2/3/4.  RETIRED/WRONG: "
                          "0.47 (a cert-2 digit transposition that the S-column-only family "
                          "missed).  Backed by the per-N_max pins in "
                          "tests/test_paper24_ho_entropy.py + tests/test_paper27_entropy.py.",
        "pattern": r"\b0\.47\b",
        "require_nearby": r"commutator|N_\{?\\?max|H_\{?HO|rel_norm|ep2b",
        "exempt_if_nearby": r"retired|corrected|transposition|previously|withdrawn",
        "files": [
            "papers/group6_precision_observations/paper_27_entropy_projection.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "p27-minnesota-contrast",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Minnesota relative-frame contrast (retraction record): CANONICAL "
                          "<00|V|10> = +17.3 MeV vs diagonal <00|V|00> = -0.55 MeV at S=0, b=1 "
                          "(ratio ~31x).  RETIRED/WRONG: +17.2 / -0.81 / ~20x (reproduced under "
                          "no convention of the production code; PM-adjudicated 2026-08-29).",
        "pattern": r"(-0\.81\$?~?MeV|\+17\.2\$?~?MeV|0,0 \\lvert V \\rvert 0,0 \\rangle = -0\.81)",
        "require_nearby": r"MeV|Minnesota|diagonal|n_\{?\\?mathrm\{rel",
        "exempt_if_nearby": r"corrected|earlier version|retired|withdrawn",
        "files": [
            "papers/group6_precision_observations/paper_27_entropy_projection.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "p34-li7-hfs-corrected",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Paper 34 Li-7 2S HFS autopsy, corrected 2026-08-28. The Track-5 Bohr-Fermi convention gives correct x 2*(m_p/m_N); at Li-7 that is 0.287204, so the baseline was low by 3.4818x. CANONICAL: baseline 288.913 MHz, final 289.13, residual -514.4 MHz = -64.0%, cliff factor 2.78x, Z_eff^effective 1.80, SCF enhancement ~2.8x. RETIRED and now WRONG: 82.977, 83.04, -720.5, -89.7%, 9.7x cliff, Z_eff 2.73, ~8x enhancement (8x is excluded by experiment, which caps it at 2.78x), and the Z_eff-scan values 82.98 / 97.58 / 1070.80 (now 288.91 / 339.75 / 3728.37). Backed by tests/test_paper34_autopsy_baselines.py.",
        "pattern": r"\b(82\.977|83\.04|720\.5|89\.7\\?%|9\.7\\times|1070\.80|97\.58)\b",
        "require_nearby": r"Li-?7|lithium|\\^7|cliff|Z_\\text\{eff\}|hyperfine|HFS",
        "exempt_if_nearby": r"retired|RETIRED|corrected|previously|withdrawn|superseded|was\\b|Track-5 convention",
        "files": [
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "p34-alkali-cliff-corrected",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Paper 34 alkali cliff sequence (tab:alkali_cliff_sequence), corrected 2026-08-28. The driver used m_e/m_N where the nuclear magneton requires m_e/m_p, suppressing every entry by m_N/m_p. CANONICAL A_fw: 144.6 / 219.8 / 10.5 / 23.5 / 36.5 MHz for Li/Na/K/Rb/Cs; cliffs -64.0 / -75.2 / -95.5 / -97.7 / -98.4%; density ratios 2.8 / 4.0 / 22.0 / 43.1 / 63.0.  The growth is NOT a power law in Z (2026-08-29 diagnosis): the closed form is (Z/Z_eff^3)(n/nu)^3 F_rel, reproducing all five to <=10%, while a fitted ~Z^1.2 misses Na by 126%; the fitted slope 1.16 survives only as a data pin discriminating the retired ~Z^2.5. RETIRED and now WRONG: 21.0 / 9.6 / 0.3 (x3) MHz, the -94.8 / -98.9 / -99.9 / -100.0% cliffs, the 19 / 92 / 851 / 3635 / 8306 ratios, the ~Z^2.5 growth, and the 'scales inversely with Z' reading (that compared bare-hydrogenic Li against FrozenCore Cs; like-for-like the cliff deepens monotonically).",
        "pattern": r"(8\{?,?\}?306|3\{?,?\}?635|851\\times|Z\^\{2\.5\}|41\.5|0\.665)",
        "require_nearby": r"alkali|cliff|Cs|Rb|contact|enhancement|Z_\\text\{eff\}",
        "exempt_if_nearby": r"retired|RETIRED|corrected|previously|withdrawn|superseded",
        "files": [
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
        ],
    },
    {
        "id": "p34-lamb-fns-2s",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Paper 34 H 2S-2P Lamb autopsy FNS row, corrected 2026-08-28. The printed +1.18 MHz is the 1S finite-size shift (r_p ~ 0.868 fm) in an n=2 table; the 2S value is +0.138 MHz at r_p = 0.8409 fm (n^3 = 8 smaller). Confirmed by the paper's own He+ autopsy, whose Lamb_FNS ratio 63.599 implies 8.7913/63.599 = 0.13823. CANONICAL: FNS +0.138, sum 1056.13, residual +1.72, Layer-2 net -1.06 MHz. RETIRED and now WRONG: +1.18, 1057.17, +0.68, and the -0.02 MHz 'empirical near-cancellation' reading, which is WITHDRAWN.",
        "pattern": r"(1057\.17|near-cancellation)",
        "require_nearby": r"Lamb|FNS|Layer-2 net|autopsy|sum",
        "exempt_if_nearby": r"retired|RETIRED|corrected|previously|withdrawn|superseded|was\\b|no longer",
        "files": [
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
        ],
    },
    {
        "id": "p27-ep2b-entropy-nmax2",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Paper 27 tab:ep2b corrected entropies, canonical after the v5.0.0 retraction (the moshinsky N_tot guard). CANONICAL S_full = 0.0671 / 0.0716 / 0.0833 nats at N_max = 2 / 3 / 4, with E_full = 21.6538 / 21.6279 / 21.5442 MeV and commutator 0.74 / 0.63 / 0.67. RETIRED and now WRONG: S = 0 (identically zero) and any other value in these rows. Backed by tests/test_paper27_entropy.py (expected dict) and tests/test_paper24_ho_entropy.py. Added 2026-08-28 after the group6 FULL run found the family MISSING despite the maintenance rule: the retraction corrected these numbers and C16 guards only the retraction phrases, so a numeric drift in the tab:ep2b row passed both gates.",
        "capture": r"21.6538\s*&(?:\s*[\d.]+\s*&)*?\s*(0\.0\d+)",
        "canonical": "0.0671",
        "exempt_if_nearby": r"retracted|RETRACTED|withdrawn|corrected|previously|retired|superseded|artifact|N_tot",
        "files": [
            "papers/group6_precision_observations/paper_27_entropy_projection.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "p27-ep2b-entropy-nmax3",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Paper 27 tab:ep2b corrected entropies, canonical after the v5.0.0 retraction (the moshinsky N_tot guard). CANONICAL S_full = 0.0671 / 0.0716 / 0.0833 nats at N_max = 2 / 3 / 4, with E_full = 21.6538 / 21.6279 / 21.5442 MeV and commutator 0.74 / 0.63 / 0.67. RETIRED and now WRONG: S = 0 (identically zero) and any other value in these rows. Backed by tests/test_paper27_entropy.py (expected dict) and tests/test_paper24_ho_entropy.py. Added 2026-08-28 after the group6 FULL run found the family MISSING despite the maintenance rule: the retraction corrected these numbers and C16 guards only the retraction phrases, so a numeric drift in the tab:ep2b row passed both gates.",
        "capture": r"21.6279\s*&(?:\s*[\d.]+\s*&)*?\s*(0\.0\d+)",
        "canonical": "0.0716",
        "exempt_if_nearby": r"retracted|RETRACTED|withdrawn|corrected|previously|retired|superseded|artifact|N_tot",
        "files": [
            "papers/group6_precision_observations/paper_27_entropy_projection.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "p27-ep2b-entropy-nmax4",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Paper 27 tab:ep2b corrected entropies, canonical after the v5.0.0 retraction (the moshinsky N_tot guard). CANONICAL S_full = 0.0671 / 0.0716 / 0.0833 nats at N_max = 2 / 3 / 4, with E_full = 21.6538 / 21.6279 / 21.5442 MeV and commutator 0.74 / 0.63 / 0.67. RETIRED and now WRONG: S = 0 (identically zero) and any other value in these rows. Backed by tests/test_paper27_entropy.py (expected dict) and tests/test_paper24_ho_entropy.py. Added 2026-08-28 after the group6 FULL run found the family MISSING despite the maintenance rule: the retraction corrected these numbers and C16 guards only the retraction phrases, so a numeric drift in the tab:ep2b row passed both gates.",
        "capture": r"21.5442\s*&(?:\s*[\d.]+\s*&)*?\s*(0\.0\d+)",
        "canonical": "0.0833",
        "exempt_if_nearby": r"retracted|RETRACTED|withdrawn|corrected|previously|retired|superseded|artifact|N_tot",
        "files": [
            "papers/group6_precision_observations/paper_27_entropy_projection.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "p23-nuclear-resource-counts",
        "scope": "group3 group4",
        "severity": "fail",
        "canonical_note": "Paper 23 nuclear qubit Hamiltonians, corrected "
                          "2026-08-22 after the N_tot truncation was removed "
                          "from geovac/nuclear/moshinsky.py (see the Paper 24 "
                          "retraction). CANONICAL: deuteron 688 non-I Pauli "
                          "(80 Z-only + 608 XY), 1-norm 383.7 MeV; He-4 828 "
                          "non-I Pauli, 1-norm 511.8 MeV (no Coulomb) / 507.2 "
                          "MeV (with). Qubit counts UNCHANGED at 16. The "
                          "structural claim survives exactly: 828/688 = "
                          "+20.3%, identical to the retired 712/592 = +20.3%. "
                          "RETIRED and now WRONG: 592, 712, 512 XY, 614 "
                          "(the composed nuclear-electronic total, now 710 = "
                          "688+10+12, measured + pinned), 342.2, "
                          "466.9, 462.4, and the HO ground-state energy "
                          "22.185 MeV (corrected to 21.6538 / 21.6279 / "
                          "21.5442 at N_max = 2 / 3 / 4). Backed by "
                          "tests/test_paper23_resource_counts.py and "
                          "tests/test_paper24_ho_entropy.py.",
        "pattern": r"\b(592|712|614|466\.9|462\.4|342\.2|22\.185)\b",
        "require_nearby": r"Pauli|1-norm|\$1\$-norm|deuteron|He-?4|"
                          r"helium|MeV|non-I|qubit|E_?0|ground[- ]state",
        "exempt_if_nearby": r"retracted|RETRACTED|withdrawn|corrected|"
                            r"previously published|artifact|retired|"
                            r"superseded|N_tot|truncation|old guard",
        "files": [
            "papers/group4_quantum_computing/paper_23_nuclear_shell.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/synthesis/*.tex",
            "docs/claim_test_matrix.md",
            "docs/validation_benchmarks.md",
        ],
    },
    {
        "id": "p58-census-builder",
        "scope": "group2",
        "severity": "fail",
        "canonical_note": "Paper 58 tab:census g row. Genuine 2,944 vs "
                          "builder 214 at n_max=2 (13.8x); 114,280 vs 7,600 "
                          "at n_max=3 (15.0x). BOTH columns use the same "
                          "axial rule m_p+m_r=m_q+m_s; the builder column is "
                          "the genuine column restricted to all-four-on-one-"
                          "center quartets. Backed by "
                          "tests/test_paper58_census.py. Earlier drafts of "
                          "the ratio as 13.7x/15.04x are rounding variants; "
                          "any OTHER builder count (e.g. 780 or 484, the "
                          "same-center readings refuted 2026-08-22) is wrong.",
        "pattern": r"\b(780|484)\b",
        "require_nearby": r"builder|census|inflation|permitted",
        "exempt_if_nearby": r"refuted|REFUTED|not the builder|wrong reading"
                            r"|hypothes|superseded",
        "files": [
            "papers/group2_quantum_chemistry/paper_58_abelian_residue.tex",
            "docs/claim_test_matrix.md",
        ],
    },
    {
        "id": "t2-collinear-anchor",
        "scope": "group2",
        "severity": "fail",
        "canonical_note": "T2 collinear = 0.395355765901713964325229296804847564... "
                          "(66 digits certified v4.104.0 via the (k,w) refactorization, "
                          "Paper 59 eq:kw; u1/u2 runs cross-validate 83). The pre-(k,w) "
                          "anchor 0.3953557659017139641 is WRONG from digit 19 "
                          "(...641 vs ...6432) and may appear only as an explicitly "
                          "superseded historical value.",
        "pattern": r"0\.3953557659017139641",
        "require_nearby": r"T2|collinear|anchor|ANCHOR",
        "exempt_if_nearby": r"superseded|18 digits|correct to 18|OLD|old anchor"
                            r"|frozen anchor|regression lock|Regression lock",
        "files": [
            "papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
            "tests/test_paper59_t2_value.py",
            "docs/claim_test_matrix.md",
        ],
    },
    {
        "id": "dirac-casimir-s3-sign",
        "scope": "group6",
        "severity": "fail",
        "canonical_note": "Dirac S^3 Casimir = +17/480 (POSITIVE). E = -1/2 "
                          "zeta_{|D|}(-1) = -1/2*(-17/240): the half-integer Dirac "
                          "shift makes zeta_{|D|}(-1) itself negative, so the fermion "
                          "-1/2 factor returns a POSITIVE Casimir -- same sign class "
                          "as the scalar +1/240 (Paper 35 KG-5 derivation, verified "
                          "numerically to 40 dps). The group6 DELTA run (2026-07-04) "
                          "caught a wrong-direction 1st-cert 'fix' that had flipped all "
                          "6 P35 loci + the code + the test to -17/480; reverted to "
                          "+17/480 across paper/code/tests. The NEGATIVE -17/480 is the "
                          "retired sign error.",
        # a NEGATIVE 17/480 (minus in front) is now the retired sign error
        "pattern": r"-\s*17/480|-\s*\\tfrac\{17\}\{480\}|-\s*\\frac\{17\}\{480\}",
        "require_nearby": r"Dirac|Casimir|zeta_\{\|D\|\}|KG-5|full.?[Dd]irac",
        "exempt_if_nearby": r"historical|stale|previously|corrected|heuristic"
                            r"|earlier|naive|retired|reverted|wrong-direction",
        "files": [
            "papers/group6_precision_observations/paper_35_time_as_projection.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "library-z-span",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "library span = Z=1--56 (H through Ba; SrH Z=38, BaH "
                          "Z=56 registry-probed). Decided v4.58.0 M-C; the 6th "
                          "cert found 5 stale Z=1--36 loci (second-locus class).",
        "capture": r"Z\s*=?\s*1\s*\$?\s*--\s*(\d{2})",
        "canonical": "56",
        # only the LIBRARY-span statements are in this family; bare periodic-row
        # prose ("First-row (Z=1--10) atoms...") is a different, legitimate quantity
        "require_nearby": r"librar|spanning|systems|H\s+through\s+Ba",
        "exempt_if_nearby": r"historical|previously|was\s+corrected|stale",
        "files": GROUP4_FILES,
    },
    {
        "id": "pauli-advantage-floor",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "raw-JW Pauli advantage floor = 51x (equal-qubit "
                          "table: 51/746/1712; P20: '51--1,712x'). The 7th cert "
                          "found a drifted '190x--1,712x' floor (M2).",
        "capture": r"(\d{2,4})\s*\$?\\times\$?\s*--\s*1\{?,\}?712",
        "canonical": "51",
        "exempt_if_nearby": r"historical|previously|stale",
        "files": GROUP4_FILES,
    },
    {
        "id": "library-count",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "library = 37 systems (35 composed + He + H2), "
                          "decided PI 2026-06-28; retired wrong counts 28/30/38/40.",
        "pattern": r"\b(?:28|30|38|40)\s+systems\b",
        "exempt_if_nearby": r"was|stale|historical|corrected|retired|previously",
        "files": GROUP4_FILES,
    },
    {
        "id": "balanced-lih-binds-at-3015",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "balanced LiH binds at the COMPUTED R_eq=3.227 bohr "
                          "(7.0% above the experimental 3.015); 'binds at 3.015' "
                          "was the v4.56.0 M2 finding (recurred v4.57.0 in the "
                          "synthesis).",
        "pattern": r"binds[^.\n]{0,60}3\.015",
        "exempt_if_nearby": r"experimental|7\.0\s*\\?%|above|computed",
        "files": GROUP4_FILES,
    },
    {
        "id": "lih-onenorm-stale",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "composed LiH 1-norm = 32.6 Ha live (0.95x vs STO-3G "
                          "34.3); the stale 33.3 / 0.97x pair retired v4.60.0 "
                          "(PI-directed corpus-wide 2026-07-01). 0.97 is scoped "
                          "to 1-norm proximity (the l-parity Pauli-ratio 0.97 "
                          "cells are a different, legitimate quantity).",
        "pattern": r"\b33\.3\b\s*~?Ha|\\lambda\s*=\s*33\.3"
                   r"|1-norm[^.\n]{0,40}\b0\.97\b|\b0\.97\b\$?\\times\$?[^.\n]{0,25}1-norm",
        "exempt_if_nearby": r"historical|stale|rested\s+on|retired",
        # widened beyond group4 2026-08-22: the full-run panel found a live 33.3 Ha
        # locus in papers/group2 (Paper 19), the second-locus class this family exists
        # to catch.  Any paper quoting the composed-LiH 1-norm is in scope.
        "files": GROUP4_FILES + [
            "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex",
            "papers/group2_quantum_chemistry/paper_58_abelian_residue.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
        ],
    },
    {
        "id": "beh2-h2o-qpe-onenorm-vintage",
        "scope": "group4",
        "severity": "fail",
        "canonical_note": "QPE-regime 1-norm cells (identity-included convention), "
                          "live-builder values pinned 2026-07-02 (8th cert): BeH2 "
                          "balanced 306.4 / composed-with-PK 373.4 (354.9 was the "
                          "deprecated legacy-builder vintage), H2O balanced 1,511. "
                          "Retired variants: 354.9, 304.7, 1{,}509 (as the balanced "
                          "H2O 1-norm).",
        "pattern": r"\b354\.9\b|\b304\.7\b|1\{,\}509~?Ha",
        "exempt_if_nearby": r"legacy|previously\s+printed|vintage|historical|stale",
        "files": GROUP4_FILES,
    },
    {
        "id": "paper60-atomic-sublinear-exponent",
        "scope": "group2",
        "severity": "fail",
        "canonical_note": "Atomic isoenergetic 1-norm sublinear exponent, HEADLINE form "
                          "\\|M\\|_1 ~ K^{0.84} (full s+p+d+f basis, eq:sublinear). 0.78 is the "
                          "legitimate s-only exponent (bare, in prose) and is NOT captured by "
                          "this family, which anchors on the \\|M\\|_1~K^{...} headline form. "
                          "W1 (2026-08-18 /qa paper 60) retired the abstract's headline K^{0.78}.",
        "capture": r"\\\|M\\\|_1\\sim\s*K\^\{(0\.\d+)\}",
        "canonical": "0.84",
        "require_nearby": r"sublinear|block-encoding|configuration|secular",
        "exempt_if_nearby": r"historical|stale|previously|s-only|retired|was|naive",
        "files": ["papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"],
    },
    {
        "id": "paper60-molecular-lambda-exponent",
        "scope": "group2",
        "severity": "fail",
        "canonical_note": "Molecular N-electron STANDARD block-encoding 1-norm exponent = "
                          "n_{\\rm orb}^{2.2} (polynomial, NOT sublinear; sec:manyelectron). "
                          "The point is polynomial-not-sublinear, so a wrong exponent here would "
                          "misstate the paper's honest molecular negative.",
        "capture": r"n_\{\\rm\s+orb\}\^\{(\d\.\d+)\}",
        "canonical": "2.2",
        "require_nearby": r"block-encoding|polynomial|1-norm|\\lambda|sublinear",
        "exempt_if_nearby": r"historical|stale|previously|retired|was",
        "files": ["papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"],
    },
]


def _gate_substr(argv: "list[str]") -> "str | None":
    for i, a in enumerate(argv):
        if a.startswith("--gate="):
            return a.split("=", 1)[1]
        if a == "--gate" and i + 1 < len(argv):
            return argv[i + 1]
    return None


def _resolve(globs: "list[str]") -> "list[pathlib.Path]":
    out: "list[pathlib.Path]" = []
    for g in globs:
        out.extend(sorted(ROOT.glob(g)))
    seen, uniq = set(), []
    for p in out:
        if p not in seen and p.is_file():
            seen.add(p)
            uniq.append(p)
    return uniq


def scan_entry(entry: dict, text_override: "str | None" = None):
    """Return (live_hits, exempt_hits); item = (relpath, line_no, snippet).

    text_override: scan the given text as a single pseudo-file (self-test hook).
    """
    exempt = re.compile(
        withdrawal_marker(entry["id"]) + "|" + entry["exempt_if_nearby"],
        re.IGNORECASE)
    # Table surfaces: every row lies inside every other row's +-WINDOW, so a
    # corrected sibling row exempts a defective one.  Opt in to same-line
    # exemption where the gated surface is a table (2026-09-03).
    exempt_same_line = entry.get("exempt_same_line", False)
    require = (re.compile(entry["require_nearby"], re.IGNORECASE)
               if "require_nearby" in entry else None)
    if "pattern" in entry:
        pat = re.compile(entry["pattern"], re.IGNORECASE)
        is_wrong = lambda m: True  # any match of a wrong-variant pattern
    else:
        pat = re.compile(entry["capture"], re.IGNORECASE)
        canonical = entry["canonical"]
        is_wrong = lambda m: m.group(1) != canonical

    def scan_lines(lines, rel):
        live, ok = [], []
        for i, line in enumerate(lines):
            m = pat.search(line)
            if not m or not is_wrong(m):
                continue
            lo, hi = max(0, i - WINDOW), min(len(lines), i + WINDOW + 1)
            window_txt = "\n".join(lines[lo:hi])
            if require is not None and not require.search(window_txt):
                continue  # outside this family's context (different quantity)
            snip = re.sub(r"\s+", " ", line.strip())[:160]
            exempt_txt = line if exempt_same_line else window_txt
            (ok if exempt.search(exempt_txt) else live).append((rel, i + 1, snip))
        return live, ok

    if text_override is not None:
        return scan_lines(text_override.splitlines(), "<override>")

    live_all, ok_all = [], []
    for path in _resolve(entry["files"]):
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        live, ok = scan_lines(lines, path.relative_to(ROOT))
        live_all.extend(live)
        ok_all.extend(ok)
    return live_all, ok_all


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    gate = _gate_substr(sys.argv)
    scope = f"scope '{gate}'" if gate else "ALL entries"

    # Entry selection is LOCUS-DERIVED as well as tag-based (2026-08-31
    # gate-scope audit; same fix as C16).  The hand-maintained `scope` tag
    # drifts from the loci it summarises, and the drift is invisible: a
    # --gate value matching no tag selects ZERO families and still prints
    # PASS -- the C19 failure shape.  An entry now runs whenever any of its
    # declared loci lies in the gated scope; the tag is kept as a widening
    # fallback.
    _in_scope, _scope_files, _scope_warnings = (
        qa_scopes.make_predicate(gate) if gate else (None, [], []))
    qa_scopes.emit_warnings(_scope_warnings)

    def _locus_gated(pattern: str) -> bool:
        if _in_scope is None:
            return True
        return any(_in_scope(h) for h in glob.glob(str(ROOT / pattern)))

    def selected(e: dict) -> bool:
        if gate is None or e["scope"] == "all" or gate in e["scope"]:
            return True
        return any(_locus_gated(f) for f in e.get("files", []))

    _selected = [e for e in REGISTRY if selected(e)]

    # GATE SELF-AUDIT (FULL run #4, 2026-09-03).  Family COUNT was not the
    # honest measure.  On `--gate trunk` this printed "3/25 families" and
    # PASS -- but all three were `scope: "all"` group4 families whose every
    # declared locus lies outside the trunk, so the gate examined nothing
    # while looking as though it had.  What matters is how many selected
    # families have a declared locus INSIDE the gated scope.
    def _has_gated_locus(e: dict) -> bool:
        if _in_scope is None:
            return True
        return any(_locus_gated(f) for f in e.get("files", []))

    _grounded = [e for e in _selected if _has_gated_locus(e)]
    if not _selected:
        print(f"   [scope] WARNING: --gate '{gate}' selected 0 of "
              f"{len(REGISTRY)} number families -- this run checks "
              f"NOTHING.")
    if gate is not None and not _grounded:
        print(f"   [scope] ERROR: --gate '{gate}' selected "
              f"{len(_selected)} family/families, but NONE of them declares "
              f"a locus inside this scope. The gate would examine nothing "
              f"and print PASS. Add a family for this branch (C17 "
              f"maintenance rule) rather than trusting this run.")
        return 1

    _authored = [e["id"] for e in REGISTRY
                 if not e.get("exempt_if_nearby", "").strip().startswith("[retracted")]
    n_live, n_exempt = 0, 0
    print(f"   [marker] {len(_authored)} of {len(REGISTRY)} families still "
          f"rely on hand-authored exemption vocabulary rather than "
          f"`[retracted YYYY-MM-DD]`.")
    print(f"headline-number registry gate (C17)   [{scope}: "
          f"{len(_selected)}/{len(REGISTRY)} families, "
          f"{len(_grounded)} with a locus in scope]\n")
    for e in REGISTRY:
        if not selected(e):
            continue
        live, ok = scan_entry(e)
        n_live += len(live)
        n_exempt += len(ok)
        status = "clean" if not live else f"{len(live)} LIVE"
        print(f"  [{'FAIL' if live else 'ok'}] {e['id']}: {status}"
              + (f"  (exempt: {len(ok)})" if ok else ""))
        for rel, ln, snip in live:
            print(f"      {rel}:{ln}  {snip}")

    if n_live:
        print(f"\nRESULT: FAIL ({n_live} live wrong-headline "
              f"occurrence(s) in {scope}; {len(_grounded)}/"
              f"{len(REGISTRY)} grounded families)")
        return 1
    print(f"\nRESULT: PASS (no live wrong headline value in {scope}"
          f"; {len(_grounded)} grounded of {len(_selected)} selected "
          f"of {len(REGISTRY)} families"
          + (f"; {n_exempt} exempt/historical mention(s))" if n_exempt else ")"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
