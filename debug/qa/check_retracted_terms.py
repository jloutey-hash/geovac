#!/usr/bin/env python3
r"""
Deterministic retracted-claims / zombie-drift screen  --  the /qa "withdrawn-as-live" backstop.

WHY THIS EXISTS (the 2026-06-23 lesson):
  Every recurring /qa miss across the Batch-2 and Batch-3 cert arcs was NOT a
  judgment failure -- it was a MECHANICAL consistency / stale-phrase drift sitting
  in a low-salience structured region (a status-table cell, a docstring, a
  footnote, one abstract clause). LLM reviewers reading a 1400-line document
  top-to-bottom systematically under-weight those regions, so a withdrawn claim
  re-surfaces and a judgment-mode panel walks past it. Concretely, the same class
  slipped repeatedly:
    * the withdrawn SU(2)xSU(2) "Pythagorean refinement" C_3<1 form -- left in
      docstrings THREE times (v4.43.5, v4.44.0, v4.46.0 each swept it incompletely);
    * the S^7 "structural negative" carried live in Paper 50's catalogue table +
      wall-list + the synthesis while the SAME paper's erratum retracts it;
    * "Latremoliere propinquity" used as the achieved metric where the result is
      van Suijlekom STATE-SPACE GH;
    * the Batch-2 "D_W/D_CH does not close the period closure" false-closure framing
      (found in body, then synthesis, then abstract, then Paper 32 -- one class,
      four locations).

  A human/LLM sweep keeps missing sites; a grep does not. This screen moves the
  RECURRING zombie classes OUT of expensive, high-variance judgment-mode review
  and INTO a cheap, exhaustive, zero-variance deterministic check that runs every
  /qa pass for ~0 tokens (the [[feedback_deferral_is_churn]] doctrine: duplicated
  /drifting fact -> single-source / deterministic check).

HOW IT WORKS:
  KNOWN LIMIT (measured 2026-08-31, /qa trunk). Patterns are matched
  LINE BY LINE and typically use `[^.\n]{0,45}` spans, so a targeted phrase
  straddling a LaTeX line wrap cannot match. Real: Paper 38 line 534 ends
  "stated the main theorem in the" and 535 begins "Latr\'emoli\`ere
  propinquity", so that alternative never fires there.
  BUT MEASURED IMPACT IS ZERO: joining adjacent line pairs across every
  entry and every gated file surfaces 0 additional live hits, 0 on
  fail-severity entries. The Paper 38 straddle sits inside a
  \begin{remark}[history] describing an EARLIER DRAFT -- disclosed history,
  which the exemption exists for. So this is latent fragility, not a live
  blind spot; do not re-derive it as a corpus alarm. (Scope of that
  measurement: single-wrap straddles. Spans are <=45 chars so that
  dominates; a two-wrap straddle was not tested.)

  THE REGISTRY (below) lists each retracted claim as {pattern, exempt_if_nearby,
  files, severity, scope}. For every pattern hit the screen checks whether a
  withdrawal marker (WITHDRAWN / retracted / "false" / Erratum / "state-space GH"
  / "named gap" / ...) appears within +-WINDOW lines. A hit WITHOUT a nearby marker
  is a live zombie:
    * severity "fail"     -> FAILS the gate (exit 1) when in --gate scope;
    * severity "advisory" -> printed for review, does NOT fail (for classes too
      noisy to gate on, where legitimate mentions abound).
      NOTE: "propinquity" used to be the example here.  It was PROMOTED to
      "fail" on 2026-08-31 (PI direction) -- the trunk criteria name that exact
      overclaim, so advisory severity left the gate soft on the one claim it
      most specifically guards.  The noise it was hedging against is now
      handled by exempt_if_nearby (0 live vs 2 correctly-exempted at promotion
      time), which is the general lesson: tighten the exemption, then gate --
      do not leave a named class permanently advisory.

  It BACKS the claims/code/synthesis reviewers (guarantees exhaustive enumeration
  of the known zombie phrases) -- it does not replace adjudication of NEW classes.
  When a /qa run retires a claim, ADD its phrase here so it can never silently
  re-surface.

Exit 0 = no live fail-severity zombie in gated scope. Exit 1 = >=1.

Usage:
  python debug/qa/check_retracted_terms.py --gate group1
  python debug/qa/check_retracted_terms.py            # all entries
  python debug/qa/check_retracted_terms.py --all      # also print exempt (compliant) hits
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

WINDOW = 5  # +- lines within which a legacy exempt_if_nearby vocabulary exempts a hit
MARKER_WINDOW = 2  # +- lines for the STANDARDIZED per-entry marker (FULL #6 fix,
# 2026-09-05): tight enough to cover a hard-wrapped sentence or an equation block
# whose marker lands 1-2 physical lines from the matched span, but NOT the +-5
# that let a marker shelter a DISTINCT occurrence of the same entry lines away
# (P7 L127 "recovering the exact Coulomb degeneracy" was sheltered by the L123
# marker, 4 lines up, for a different sub-claim).  Measured separation on the
# FULL #6 loci: every legitimate wrap/equation case sits at distance 1-2 from
# its marker; the one genuine zombie (P7 L127) at distance 4.

# ---------------------------------------------------------------------------
# THE REGISTRY -- append an entry whenever a /qa run retires/withdraws a claim.
# Each entry: a retracted phrase (pattern) that must NOT appear as LIVE; it is
# exempt only when a withdrawal marker (exempt_if_nearby) sits within +-WINDOW
# lines.  Patterns are case-insensitive raw regex.  `files` are ROOT-relative
# globs.  `severity`: "fail" (gates) | "advisory" (reports only).
# ---------------------------------------------------------------------------

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


# Entries predating the 2026-09-04 `cited_by` requirement.  Ratcheted so
# historical debt is visible without blocking;  every NEW entry must
# declare its dependents (or `cited_by: {}` if genuinely none).
# ---------------------------------------------------------------------------
# LaTeX-tolerant matching (2026-09-04, /qa DELTA #6).
#
# Four successive pattern rebuilds missed loci for one reason: the pattern is
# written in prose and the corpus is typeset.  `height_B <= gamma` matched
# while `\mathrm{height}_B \le \gamma` did not; `s/p-lift decay` matched while
# `$s/p$-lift decay` did not -- and THAT sweep reported clean with five loci
# live.  Stripping markup before matching removes the class, not the instance.
#
# Matching only: the snippet reported to the user is the original line, so
# quoted evidence and line numbers are unchanged.  Deliberately conservative --
# it removes inline math delimiters and the wrappers that split a word, and
# normalises the LaTeX comparison operators.  It does not render LaTeX.
_MARKUP_WRAPPER = re.compile(
    r"\\(?:emph|textbf|textit|texttt|mathrm|mathbf|mathbb|text|mathit)\s*\{([^{}]*)\}")
_MARKUP_OPS = ((r"\\leq", "<="), (r"\\le\b", "<="), (r"\\geq", ">="),
               (r"\\ge\b", ">="), (r"\\neq", "!="), (r"\\to\b", "->"))


# Unicode folding (2026-09-07, /qa paper_61 DELTA).  Same class as the LaTeX
# gap above, one alphabet over:  the patterns are ASCII and the corpus is
# typeset Unicode, so `Sp[_ ]?4 ... mathbb Z` could not see `Sp₄(ℤ)` and
# `Broadhurst--?Mellit` could not see `Broadhurst–Mellit`.  MEASURED: both
# Paper-61 entries written this session missed 4/4 of their live survivors,
# every miss traceable to this.
#
# En/em dash folds to `--`, NOT to `-`:  three entries (circle-fejer-constant-
# 4-over-pi, graph-spectrum-attribution, rename-pass-residue-p32) contain
# strict `--` tokens, and `X--?Y` matches `X--Y` and `X-Y` both, so folding UP
# is additive where folding DOWN would have broken them.
_UNICODE_FOLD = {
    "\u2013": "--", "\u2014": "--",           # en dash, em dash
    "\u2010": "-", "\u2011": "-", "\u2212": "-",   # hyphen, nb-hyphen, minus
    "\u2080": "0", "\u2081": "1", "\u2082": "2", "\u2083": "3", "\u2084": "4",
    "\u2085": "5", "\u2086": "6", "\u2087": "7", "\u2088": "8", "\u2089": "9",
    "\u2124": "Z", "\u2102": "C", "\u211a": "Q", "\u211d": "R", "\u2115": "N",
}


def _strip_markup(line: str) -> str:
    """Prose-ify a LaTeX/Markdown/Unicode line for pattern matching."""
    out = line
    for _ in range(3):                        # nested \emph{\textbf{...}}
        new = _MARKUP_WRAPPER.sub(r"\1", out)
        if new == out:
            break
        out = new
    for pat, rep in _MARKUP_OPS:
        out = re.sub(pat, rep, out)
    out = out.replace("$", "")                # inline math delimiters
    out = re.sub(r"\\[,;!]|\\ ", " ", out)     # thin spaces
    for _k, _v in _UNICODE_FOLD.items():      # typeset Unicode -> ASCII
        out = out.replace(_k, _v)
    out = out.replace("**", "").replace("`", "")   # markdown bold/code split
    return out


CITED_BY_BASELINE = {
    "all-compact-lie-groups-universality",
    "bare-graph-n2-minus-1-attribution",
    "batch2-false-closure",
    "brown-surjection-attribution",
    "cheeger-simons-cone-attribution",
    "circle-fejer-constant-4-over-pi",
    "combined-triple-ko1-additive",
    "cp2-50pct-vs-fit-floor",
    "dgv-graph-form-tautology",
    "drake-swainson-3d-mistranscription",
    "emn-catalan-negative-attribution",
    "fock-coupling-one-sixteenth-prefactor",
    "latremoliere-propinquity-named-for-gh-rate",
    "lorentzian-literal-identification-krein",
    "organics-in-library",
    "p24-entanglement-rigidity",
    "p27-ho-zero-entropy-rigidity",
    "p34-alkali-uniformity",
    "p34-ee-eigen-closed-form",
    "p34-lamb-near-cancellation",
    "p45-kplus-compression-theorem-live",
    "pair-diagonal-as-exact-sparsity",
    "pairdiag-composed-scaling-livesd",
    "pairdiag-density-as-the-angular-density",
    "pairdiag-density-pipeline-realizes",
    "propinquity-as-achieved-metric",
    "propinquity-as-achieved-metric-group3",
    "propinquity-as-achieved-metric-group5",
    "rename-pass-residue-p32",
    "s-p-splitting-retired-waypoints",
    "s7-structural-negative",
    "sbh-phi2-prefactor",
    "su2-kinetic-equals-L1",
    "tc-qubit-validated-success",
    "wald-factor2-cone-coefficient",
    "withdrawn-c3op-envelope-sqrt",
    "withdrawn-pythagorean-mechanism",
}

REGISTRY = [
    {
        "id": "sp-splitting-as-convergence-evidence",
        "note": "FULL run #4 + DELTA #5.  Distinct from the MECHANISM entry: "
                "this one guards the EVIDENTIAL reading -- the s/p decay "
                "listed as one of the two numerical indications supporting "
                "the graph -> S^3 continuum reading.  Paper 7 now states "
                "that the leg 'carries no evidential weight ... which now "
                "rests on the lambda_max saturation alone', and the synthesis "
                "that 'the saturation of the spectral bound is the only leg "
                "that carries weight'.  Any locus still pairing the two as "
                "joint support contradicts that.  The quantity itself is "
                "fine to report;\u00a0what is retired is citing it as evidence.",
        "pattern": r"s/p[- ]lift decay"
                   r"|degeneracy recovery"
                   r"|a degeneracy that is recovered"
                   r"|degeneracy lift[^.\n]{0,60}decays"
                   r"|saturation[^.\n]{0,40}and the (?:decay of the )?s/p",
        # Only this entry's own standardized marker exempts it:
        # `withdrawal_marker(id)` is OR'd in automatically, so the
        # per-entry pattern must never match.  A bare r"\[retracted"
        # here matches ANY marker, including one naming an unrelated
        # entry -- the exact over-exemption the entry-id redesign was
        # written to kill, reintroduced in that same commit and caught
        # by DELTA #6.
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group3 trunk synthesis",
                # Documents whose ARGUMENT rests on this claim (distinct
        # from `files`, which is where its wording may appear).
        "cited_by": {
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "reviewed 2026-09-04",
            "papers/synthesis/group3_foundations_synthesis.tex": "reviewed 2026-09-04",
            "docs/claim_test_matrix.md": "reviewed 2026-09-04",
        },
"files": [
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "docs/claim_test_matrix.md",
        ],
    },
    {
        "id": "graph-s3-convergence-established",
        "note": "Trunk FULL #3/#4 re-pricing, still live at DELTA #5.  The "
                "operator convergence L -> Delta_{S^3} and the identification "
                "of the limit manifold are NOT established:\u00a0Paper 7 says so "
                "at four loci ('not established here', 'no test in the "
                "repository establishes the limit manifold', 'neither proven "
                "nor established numerically').  What is measured is the "
                "lambda_max saturation.  Retired:\u00a0stating the convergence or "
                "the manifold identification as established/demonstrated/"
                "proven, and tagging it [SYMBOLIC PROOF] or [MEASURED].",
        "pattern": r"convergence[^.\n]{0,40}(?:is|was) demonstrated empirically"
                   r"|discrete graph converges[^.\n]{0,80}and that limit is"
                   r"|conformal equivalence between the discrete[^.\n]{0,60}established"
                   r"|discrete Fock graph's continuum limit is\s*\n?\s*conformally equivalent",
        # Only this entry's own standardized marker exempts it:
        # `withdrawal_marker(id)` is OR'd in automatically, so the
        # per-entry pattern must never match.  A bare r"\[retracted"
        # here matches ANY marker, including one naming an unrelated
        # entry -- the exact over-exemption the entry-id redesign was
        # written to kill, reintroduced in that same commit and caught
        # by DELTA #6.
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group3 trunk synthesis",
                # Documents whose ARGUMENT rests on this claim (distinct
        # from `files`, which is where its wording may appear).
        "cited_by": {
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex": "reviewed 2026-09-04",
            "papers/synthesis/group3_foundations_synthesis.tex": "reviewed 2026-09-04",
        },
"files": [
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
    },
    {
        "id": "lorentzian-strong-identification-live",
        "note": "FULL run #5, 2026-09-05.  Paper 32's Lorentzian section L2-E "
                "printed a live verdict STRONG_IDENTIFICATION_LORENTZIAN and a "
                "subsection title 'unified-strong four-witness theorem', reading "
                "as an achieved strong *Lorentzian* identification -- while the "
                "closure is signature-blind (the truncated BW boost is compact; "
                "WH7 de-compactification = convention), retracted two paragraphs "
                "below.  Relabelled 'compact_period_closure (signature-blind)'.",
        "pattern": r"strong[\\_ ]{0,2}identification[\\_ ]{0,2}lorentzian"
                   r"|unified-strong four-witness",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group1 trunk",
        "cited_by": {},
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
        ],
    },
    {
        "id": "p39-pythagorean-height-route",
        "note": "FULL run #6, 2026-09-05.  Paper 32 (L1453) described Paper 39's "
                "tensor result as proved by 'Paper 38's five-lemma machinery "
                "factor-by-factor ... joint Lipschitz-distortion height term ... "
                "via a Pythagorean refinement of the triangle inequality' -- the "
                "(B,P)-pair route Paper 39's discharge RETIRED (height leg "
                "refuted; superseded by the lifted-state route, constant 1).  A "
                "cited_by miss: P32's argument rested on Paper 39's proof route, "
                "the discharge changed it, and the dependent locus survived.  "
                "Rewritten to the lifted-state route; this entry guards the "
                "refuted-route phrasing.",
        "pattern": r"five-lemma machinery factor-by-factor"
                   r"|joint Lipschitz-distortion height",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group1 trunk",
        # Documents whose ARGUMENT rests on Paper 39's tensor proof route.
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex":
                "reviewed 2026-09-05 (FULL #6 remediation: rewritten to lifted-state)",
        },
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
        ],
    },
    {
        "id": "l5-assembly-listed-as-live-contribution",
        "note": "FULL run #5, 2026-09-05.  Paper 38's abstract listed Lemma L5 "
                "('assembly of the distance bound via an approximation pair') as "
                "a live contribution with no withdrawal marker, while L5 is "
                "withdrawn 2026-09-03 (the bound comes from the unconditional "
                "lifted-state theorem, not the approximation-pair assembly).  "
                "Marked withdrawn in the abstract.  "
                "WIDENED FULL run #7, 2026-09-06: Paper 32's thm:gh_convergence "
                "proof-sketch (L3471) opened 'The proof is a Latremoliere "
                "tunneling-pair assembly of five lemmas' -- the same L5-as-live "
                "class, in a trunk paper, restated in Paper 32's own words so "
                "the paper-38-only pattern missed it.  Reworded to the Paper 38 "
                "lifted-state route; 'tunneling-pair assembly of five' added as "
                "the discriminating pattern (fires on the retired wording, silent "
                "on the corrected 'historical five-lemma UCP-pair assembly ... "
                "recorded ... withdrawn').",
        # WIDENED 2026-09-07 (FULL cert of 58/59/60, group1 claim-impact pass).
        # The three alternatives above are Paper-38/32 wordings.  The citers
        # restate the same claim in their OWN words and were invisible:
        # Paper 52 -- "The proof is the five-lemma chain L1'-L2-L3-L4-L5 of
        # Paper 38" and "via Paper 38's five-lemma state-space Gromov--Hausdorff
        # convergence proof" -- attributes the WH1 keystone to a chain whose
        # fifth link is refuted; Paper 40's own main-theorem proof still reads
        # "Lemma L5 the assembly into the state-space GH bound", four sections
        # after its own L5 says the pair "does not supply an independent proof".
        # `files` was the other half of the failure: it listed only 38 and 32,
        # so the gate never opened 46/52/53/47/48/44 or the field guide.
        "pattern": r"assembly of the distance bound via an approximation pair"
                   r"|proves the five lemmas"
                   r"|tunneling-pair assembly of five"
                   r"|five-lemma chain"
                   r"|five-lemma state-space"
                   # NOT the bare phrase "the assembly into the state-space GH
                   # bound": Paper 38 uses it correctly, marked "in the
                   # *withdrawn* route only" and immediately followed by "Not
                   # the L5 assembly: its two height constituents are refuted".
                   # A first draft of this entry matched that mention and would
                   # have had me "fix" the one document that had it right.
                   # What is defective is the USE -- Paper 40's main-theorem
                   # proof combining L5's bound to obtain the theorem.  Lines
                   # 1831/1846 of the same file cite the same label correctly,
                   # so the label alone does not discriminate;  this anchors on
                   # the proof step.
                   r"|the bound~\\eqref\{eq:L5_bound_general\} with the asymptotic",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group1 trunk",
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex":
                "reviewed 2026-09-06 (FULL #7: proof-sketch reworded to the "
                "Paper 38 lifted-state unconditional route; five-lemma list "
                "reframed as historical/superseded, L5 withdrawn)",
        },
        "files": [
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            # Added 2026-09-07: the citer set the sweep never opened.
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex",
            "papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex",
            "papers/group1_operator_algebras/paper_47_two_rate_hybrid_convergence.tex",
            "papers/group1_operator_algebras/paper_48_krein_ms_bridge.tex",
            "papers/group1_operator_algebras/paper_52_category_iii_correspondence.tex",
            "papers/group1_operator_algebras/paper_53_disk_propinquity.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "papers/synthesis/geovac_field_guide.tex",
        ],
    },
    {
        "id": "l5-height-bound-achieved",
        "note": "FULL run #4, 2026-09-03.  Paper 38's Lemma L5 concluded "
                "height_B <= gamma_nmax and height_P = 0.  Both are FALSE by "
                "a one-line witness:\u00a0B_nmax is a finite-band reconstruction "
                "(envelope N <= 2 n_max - 1), so the unit-Lipschitz ball "
                "contains f with B(f) = 0, for which the height quantity is "
                "exactly 1 -- hence height_B == 1 at every cutoff, and the "
                "bound fails wherever gamma < 1 (every n_max >= 6).  The "
                "v5.4.3/v5.4.4 reading of the measured crossing as 'an open "
                "check on the panel-side quantity' is retired with it:\u00a0the "
                "panel entries normalise INTO the ball, so they lower-bound "
                "the supremum and an entry above gamma contradicts the bound. "
                "Two printed steps are separately invalid (the compressed-"
                "multiplier inequality runs the wrong way;\u00a0the good-kernel "
                "estimate bounds a sup-norm, not a Lipschitz-seminorm "
                "difference).  WHAT SURVIVES:\u00a0thm:main_unconditional uses "
                "only reach-type estimates and is untouched, so the "
                "convergence statement and WH1 stand -- do not read this as a "
                "retraction of the keystone.",
        # REBUILT 2026-09-04 (/qa DELTA #5): four of the five original
        # alternatives matched ZERO lines corpus-wide -- one required a
        # literal \gamma so no plain-text locus could match, two
        # required a thin-space \;\le spelling that occurs nowhere,
        # and one embedded \n under a line-by-line scanner.  The entry
        # guarding the run's biggest withdrawal caught nothing.
        # These are wording-tolerant and cover LaTeX, code and prose.
        # REBUILT AGAIN 2026-09-07 (group1 impact-set pass): the subscript
        # alternation `height\}?[_ ]?P` covered `height_P`, `height P` and
        # `height}P` but NOT the BRACED subscript `\mathrm{height}_{P}`,
        # which is how the corpus actually typesets it in prose.  The
        # group1 synthesis carried `contributes $\mathrm{height}_{P} = 0$`
        # -- the refuted claim, stated as live, in the document that
        # summarises the paper it was refuted in -- and this entry reported
        # `clean` on group1 for four days.  Exactly the spelling-defeats-the-
        # pattern class already recorded in CLAUDE.md
        # (`$s/p$-lift` vs `s/p-lift`).  `\{?` admits the braced form.
        "pattern": r"height\}?[_ ]?\{?\s*B[^\n]{0,40}?(?:<=|\\le|\\leq)[^\n]{0,25}?gamma"
                   r"|height\}?[_ ]?\{?\s*P[^\n]{0,30}?(?:=|==|is)\s*(?:0|0\.0|zero)\b"
                   r"|neither confirms nor contradicts L5"
                   r"|panel-side quantity is"
                   r"|exceeds gamma from n_max"
                   r"|height_B_theoretical\(\)[^\n]{0,40}upper bound"
                   # DELTA #7 (CLAIMS-A M1/M5): the anonymous four-term max
                   # whose last two arguments ARE the refuted heights, and the
                   # remark that still credits the height reduction with the
                   # main theorem.
                   # The display wraps over three lines, so anchor on its LAST line:
                   # a gamma-term followed by the trailing ", 0\bigr)".
                   r"|\\gamma[^\n]{0,24},\s*0\s*\\bigr\)"
                   r"|\\gamma[^\n]{0,24},\s*0\s*\)\s*(?:=|\\;=)"
                   r"|reduction of the height to[^\n]{0,60}is what gives"
                   r"|all other constituents are proved above",
        "exempt_if_nearby": r"withdrawn|REFUTED|refuted|false by|retired|"
                            r"until then|read \\`\\`this|printed here until",
        "severity": "fail",
        "scope": "group1 trunk",
                # Documents whose ARGUMENT rests on this claim (distinct
        # from `files`, which is where its wording may appear).
        "cited_by": {
            "geovac/gh_convergence.py": "reviewed 2026-09-04",
            "tests/test_gh_convergence.py": "reviewed 2026-09-04",
            # Stamped after the v5.8.3 edit; DELTA #7 found the four-term assembly
            # still live nine lines below (CLAIMS-A M1).  A stamp records a review
            # OUTCOME, so it reverts to None until a review confirms.
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex": "remediated 2026-09-04 (DELTA #7); review owed (DELTA #8)",
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex": "remediated 2026-09-04 (DELTA #7); review owed (DELTA #8)",
                    "geovac/gh_convergence_tensor.py": "remediated 2026-09-04 (DELTA #7); review owed (DELTA #8)",
            "geovac/central_fejer_su2.py": "remediated 2026-09-04 (DELTA #7); review owed (DELTA #8)",
            "geovac/lorentzian_propinquity_compact_temporal.py": "remediated 2026-09-04 (DELTA #7); review owed (DELTA #8)",
            "geovac/ecosystem_export.py": "remediated 2026-09-04 (DELTA #7); review owed (DELTA #8)",
        },
"files": [
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "geovac/gh_convergence.py",
            "tests/test_gh_convergence.py",
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
            "docs/qa/trunk.done.md",
            "docs/claim_test_matrix.md",
            "tests/test_gh_convergence_tensor.py",
            "geovac/lorentzian_propinquity_compact_temporal.py",
            # Added 2026-09-07 (FULL cert of 58/59/60).  Paper 46 states
            # `\mathrm{height}_{P} \;=\; 0.` as a PROPOSITION result and cites
            # "Paper 38 S L5" as its proof; Papers 47/48/52/53 and the field
            # guide carry the same claim in their own words.  None was in this
            # list, which is why the gate reported group1 clean while the
            # refuted heights stood in six documents.
            "papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex",
            "papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex",
            "papers/group1_operator_algebras/paper_47_two_rate_hybrid_convergence.tex",
            "papers/group1_operator_algebras/paper_48_krein_ms_bridge.tex",
            "papers/group1_operator_algebras/paper_52_category_iii_correspondence.tex",
            "papers/group1_operator_algebras/paper_53_disk_propinquity.tex",
            "papers/synthesis/geovac_field_guide.tex",
            "docs/claims_register.md",
        ],
    },
    {
        "id": "sp-splitting-aliasing-mechanism",
        "note": "FULL run #4, 2026-09-03.  The s/p node-amplitude proxy was "
                "attributed to 'spectral aliasing on a compact manifold' -- "
                "boundary reflections differentially shifting eigenvalues "
                "across angular-momentum sectors -- and its decay read as "
                "'recovering the exact Coulomb degeneracy', a hallmark of "
                "SO(4).  Both are refuted by the corpus's own result that the "
                "l-blocks are DISCONNECTED components:\u00a0there is no s/p "
                "degeneracy present to alias and no cross-sector reflection "
                "to shift it.  Measured: lambda_2s = 3 exactly iff n_max = 0 "
                "(mod 3), so the reported series tracks cutoff divisibility, "
                "and the 0.39% endpoint sits on a favourable branch (1.65% at "
                "28, 2.58% at 29).  The quantity IS closed form on that "
                "branch -- (2 - 2cos(pi/(n_max-1)))/3 -- which is an upgrade; "
                "what is retired is its use as evidence of convergence.",
        # WIDENED 2026-09-04 (/qa DELTA #5): matched one locus per paper
        # where the reviewers found eleven.  It caught the exact phrase
        # in front of me when I wrote it and none of the variants the
        # corpus actually uses.
        "pattern": r"spectral aliasing"
                   r"|finite-size aliasing"
                   r"|standing-wave reflection"
                   r"|absorbing wall that differentially"
                   r"|differential(?:ly)?[- ]weighted connectivity"
                   r"|differential connectivity between angular"
                   r"|recovered degeneracy is the hallmark"
                   r"|recovering the exact Coulomb degeneracy"
                   r"|continuing toward the exact Coulomb degeneracy"
                   # WIDENED FULL #8, 2026-09-06 (CLAIMS-A M1): Paper 1's intro
                   # roadmap paraphrased the withdrawn mechanism as
                   # "s/p degeneracy breaking by the graph Laplacian ...
                   # vanishes in the continuum limit" -- LIVE, no marker; the
                   # graph has NO s/p degeneracy (l-blocks disconnected).  A
                   # pure paraphrase the panel caught, invisible to the older
                   # patterns.  Discriminates: fires on the retired wording,
                   # silent on the corrected "node-amplitude proxy ... not a
                   # graph degeneracy".
                   r"|degeneracy breaking by the graph",
        # DELTA #7: exemption vocabulary REMOVED.  "residue" here let Paper 1's
        # Conclusion keep "confirms this is spectral aliasing" live, because
        # "set by that residue" sits in the same bullet (CLAIMS-B M1).  Only
        # the entry's own withdrawal marker exempts now.
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group3 trunk synthesis",
                # Documents whose ARGUMENT rests on this claim (distinct
        # from `files`, which is where its wording may appear).
        "cited_by": {
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "reviewed 2026-09-04",
            "papers/synthesis/group3_foundations_synthesis.tex": "reviewed 2026-09-04",
        },
"files": [
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "docs/claim_test_matrix.md",
        ],
    },
    {
        "id": "hopf-base-label-for-4-over-pi",
        "note": "Trunk DELTA #4, 2026-09-03.  Paper 38's asymptotic constant "
                "is a quotient of UNIT-RADIUS sphere volumes, "
                "2 Vol(S^2)/Vol(S^3) = 4/pi (numerically Vol(S^2)/pi^2).  It "
                "is NOT the Hopf fibration's base-to-total ratio:\u00a0the Hopf "
                "map is a Riemannian submersion S^3(r) -> S^2(r/2), so that "
                "ratio is 1/(2 pi) at unit radius and 1/(4 pi) in the "
                "dual-Coxeter metric -- neither is 4/pi or 2/pi.  The M1 "
                "MECHANISM (the k = 0 volume slot of the master Mellin "
                "engine) and its period ring Q[pi, 1/pi] are UNAFFECTED; "
                "only the name is a misnomer, and the name may still be used "
                "for the k = 0 slot or for Paper 25's Vol(S^2)/4 = pi.  What "
                "is retired is reading 4/pi ITSELF as a base-to-total ratio "
                "or as a Haar quotient SU(2)/U(1).  "
                "tests/test_p38_metric_convention.py.",
        # Widened 2026-09-03 (FULL #4): the pattern required the literal
        # "Hopf-base measure" adjacent to the constant, so Paper 32's
        # "a Vol/Mellin-measure factor on the spatial $S^3$ Hopf base
        # ($\Vol(S^2)/\pi^2 = 4/\pi$ ...)" -- the same retracted reading in
        # different words -- was invisible.  I had written the pattern around
        # the wording in front of me rather than around the class.
        "pattern": r"on the spatial \$?S\^?\{?3\}?\$? Hopf base"
                   r"|factor on the[^.\n]{0,30}Hopf base"
                   r"|Hopf[- ]base measure factor in the standard"
                   r"|Hopf[- ]base measure (?:factor )?of \$?\\?(?:s)?three"
                   r"|Hopf[- ]base measure of \$\\sthree"
                   r"|(?:as|is) the\s+Hopf[- ]base measure[^.]{0,60}"
                   r"(?:4\s*/\s*\\pi|\\Vol\(S\^\{?2\}?\)\s*/\s*\\?pi\^\{?2\}?)"
                   # WIDENED 2026-09-07 (group1 impact-set pass): the trailing
                   # form required "is the"/"as the" before the label, so the
                   # parenthetical apposition -- "asymptotic rate $4/\pi$ (the
                   # Hopf-base measure ...)", which is how Papers 20 and 42
                   # actually write it -- did not match.  `\(` admits it.
                   r"|(?:4\s*/\s*\\pi|\\Vol\(S\^\{?2\}?\)\s*/\s*\\?pi\^\{?2\}?)"
                   r"[^.]{0,80}(?:is|as|\()\s*the\s+Hopf[- ]base measure"
                   r"|Hopf[- ]base measure of \$?\\SU\(2\)"
                   r"|\\SU\(2\)\s*/\s*\\?U?one?\(1\)\s+Haar\s+normalisation",
        "exempt_if_nearby": r"misnomer|corrected 2026-09-03|withdrawn|not the Hopf"
                            r"|historically|k = 0 volume|volume-ratio note",
        "severity": "fail",
        "scope": "group1 group3 trunk",
                # Documents whose ARGUMENT rests on this claim (distinct
        # from `files`, which is where its wording may appear).
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "reviewed 2026-09-04",
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex": "reviewed 2026-09-04",
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "reviewed 2026-09-03",
            "papers/synthesis/group3_foundations_synthesis.tex": "reviewed 2026-09-03",
            "papers/synthesis/geovac_field_guide.tex": "reviewed 2026-09-03",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "reviewed 2026-09-03",
            # Added 2026-09-07 (group1 impact-set pass).  These three cite
            # the rate constant and their argument uses it, so they are
            # dependents, not merely places the wording appears.  All three
            # carried the retired reading and were fixed in that pass.
            "papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex": "reviewed 2026-09-07",
            "papers/group1_operator_algebras/paper_43_lorentzian_extension.tex": "reviewed 2026-09-07",
            "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex": "reviewed 2026-09-07",
        },
"files": [
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/synthesis/geovac_field_guide.tex",
            # ADDED 2026-09-07 (group1 impact-set pass).  These three carry
            # the retired reading and were never scanned: the entry's file
            # list was built from where the label was FIXED, not from where
            # the constant is CITED.  P42 states "4/pi = Vol(S^2)/pi^2
            # identifies as the same Hopf-base measure factor" in three
            # places; P20 writes "asymptotic rate 4/pi (the Hopf-base
            # measure)".  The gate reported this entry clean on group1 for
            # four days as a result.
            "papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex",
            "papers/group1_operator_algebras/paper_43_lorentzian_extension.tex",
            "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex",
        ],
    },
    {
        "id": "four-over-pi-canonical",
        "note": "Trunk DELTA #4, 2026-09-03.  The rule Cas(ad) = h^v that "
                "Papers 38 and 40 declare is the CORPUS's rule, not the "
                "field-standard one.  Kac's basic form (theta|theta) = 2 "
                "gives Cas(ad) = (theta|theta+2 rho) = 2 h^v, hence the "
                "SU(2) sphere of radius sqrt(2) and the constant "
                "2 sqrt(2)/pi;\u00a0the unit sphere gives 2/pi.  So 4/pi is the "
                "value in one consistently applied convention -- 'canonical' "
                "and the attribution of the rule to 'the standard "
                "normalisation used in representation theory and conformal "
                "field theory' are both retired.  The UNIVERSALITY statement "
                "is unaffected:\u00a0what makes it well-posed is that the rule is "
                "fixed independently of the answer, not which rule it is.  "
                "tests/test_p38_metric_convention.py.",
        "pattern": r"4\s*/\s*\\?pi is canonical"
                   r"|canonical,? not a convention artifact"
                   r"|standard \\emph\{dual-Coxeter normalisation\} used in\s*\n?\s*representation theory"
                   r"|the standard dual-Coxeter normalisation used in",
        # NOTE (gate self-audit, 2026-09-03): "convention-dependent" was in
        # this list and SWALLOWED the plant during the two-way
        # discrimination test -- the retired sentence sits next to the
        # phrase "metric-convention dependent".  Narrowed, then re-proved.
        "exempt_if_nearby": r"not\s+canonical|is retired|withdrawn|corrected 2026-09-03"
                            r"|Attribution note|Earlier drafts called",
        "severity": "fail",
        "scope": "group1 trunk",
                # Documents whose ARGUMENT rests on this claim (distinct
        # from `files`, which is where its wording may appear).
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "reviewed 2026-09-04",
            "papers/synthesis/group3_foundations_synthesis.tex": "reviewed 2026-09-04",
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "reviewed 2026-09-03",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "reviewed 2026-09-04",
        },
"files": [
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "geovac/gh_convergence.py",
            "geovac/central_fejer_su2.py",
        ],
    },
    {
        "id": "circle-fejer-constant-4-over-pi",
        "note": "Trunk FULL run #2, 2026-09-02 (P38 finding, propagated to "
                "P40 + group1 synthesis + outreach note N1).  The circle "
                "(Fejer on T^1, probability-normalised) first-moment "
                "constant is 2/pi, not 4/pi: m_n = pi/2 - (4/pi) "
                "sum_{k odd<n}(1-k/n)/k^2 ~ (2/pi) log n / n (exact closed "
                "form; tests/test_trunk_qa_fejer_4_over_pi.py).  The SU(2) "
                "constant 4/pi is TWICE the circle constant -- an "
                "Observation, not a derivation.  Retired: 'the same "
                "constant on both sides', '4/pi on each circle', "
                "'preserved across this transition', the 'Stein--Weiss "
                "SS I.1' citation, the 'classical Stein--Weiss estimate' / "
                "'real-line Stein--Weiss constant' naming for the circle "
                "value, the P40 torus-factor corollary extension built on "
                "it, and the arithmetic slip 2Vol(S^1)/Vol(SU(2)) (= 2/pi) "
                "for 4/pi (correct: 2Vol(S^2)/Vol(SU(2)) = Vol(S^2)/pi^2). "
                "'Stein--Weiss sharpening' as the name of App. A's "
                "sum-rule method stays.",
        "pattern": r"F_n\(\\theta\)[^\n]{0,40}\\sim[^\n]{0,12}(?:4\s*/\s*\\pi|\\frac\{4\}\{\\pi\})"
                   r"|(?:same|identical)\s+on\s+both\s+sides"
                   r"|on both sides,? which we read"
                   r"|4/\\pi\$? on each circle"
                   r"|(?:each|every)\s+circle\s+factor\s+carries\s+the\s+same"
                   r"|tori\s+included"
                   r"|preserved across this transition"
                   r"|stein_weiss1971\}\s*\\S\s*I\.1"
                   r"|classical Stein--Weiss estimate"
                   r"|same constant that appears"
                   r"|real-line Stein--Weiss"
                   r"|1\$?D Stein--Weiss"
                   r"|Stein--Weiss constant"
                   r"|Stein--Weiss circle"
                   r"|Fej\\'er--Stein--Weiss"
                   r"|2\s*\\?(?:mathrm\{)?Vol\}?\(S\^1\)\s*/\s*\\?(?:mathrm\{)?Vol\}?\(\\?(?:mathrm\{)?SU\}?\(2\)\)",
        "exempt_if_nearby": r"withdrawn|WITHDRAWN|corrected 2026-09|formerly|2/\\?pi|twice|half the",
        "severity": "fail",
        "scope": "group1 trunk",
        "files": [
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "docs/outreach/note_n1_su2_truncations.tex",
            "tests/test_trunk_qa_fejer_4_over_pi.py",
            "geovac/central_fejer_su2.py",
        ],
    },
    {
        "id": "dgv-graph-form-tautology",
        "note": "Trunk FULL run #2, 2026-09-02 (P32 finding).  Paper 32's "
                "Definition 'Graph form of D_GV' (def:D_GV_graph) and "
                "Proposition 'Equivalence of spectral and graph forms' "
                "(prop:D_equiv) were withdrawn: the proposition was "
                "tautological (weights chosen so the spectra agree, then "
                "unitary equivalence by the spectral theorem) and the "
                "module it cited (geovac/dirac_matrix_elements.py) holds "
                "closed-form matrix elements, not a hopping operator.  The "
                "paper carries ONE Dirac operator (def:D_GV_spectral); "
                "Remark rem:D_GV_no_graph_form records the withdrawal. "
                "Pinned by tests/test_paper32_dirac_module_contents.py.",
        "pattern": r"prop:D_equiv"
                   r"|def:D_GV_graph"
                   r"|Equivalence of spectral and graph forms"
                   r"|Graph form of \$D_\{\\mathrm\{GV\}\}\$"
                   r"|graph form is what one uses"
                   r"|D_\{\\mathrm\{GV\}\}\^\{\\mathrm\{graph\}\}",
        "exempt_if_nearby": r"No separate|Earlier drafts|the former|withdrawn|WITHDRAWN|retired|tautological",
        "severity": "fail",
        "scope": "group1 trunk",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "geovac/full_dirac_operator_system.py",
            "geovac/dirac_matrix_elements.py",
            "tests/test_paper32_dirac_module_contents.py",
        ],
    },
    {
        "id": "combined-triple-ko1-additive",
        "note": "Trunk FULL run #2, 2026-09-02 (P32 + code finding).  The "
                "H1 combined almost-commutative triple was labelled "
                "'KO-dimension 3 + 6 = 9 = 1 (mod 8)' while the MEASURED "
                "signs (J^2 = -I, JD = +DJ) are the (eps, eps') = (-, +) "
                "KO-3 pair; KO-1 is (+, -).  The additive rule is a "
                "theorem about the GRADED product D_GV (x) gamma_F + 1 (x) "
                "D_F (Dabrowski-Dossena 2011); the module builds the "
                "ungraded sum D_GV (x) 1 + gamma_GV (x) D_F, for which the "
                "GV factor's own signs survive.  fluctuated_dirac default "
                "epsilon_prime moved -1 -> +1 to match.  Relabelled KO-3, then "
                "(PI direction, same day) the finite-cutoff KO label was dropped "
                "in favour of the measured sign triple (-, +, +); "
                "the additive statement may not re-surface as the label.",
        "pattern": r"3\s*\+\s*6\s*(?:=\s*9\s*)?(?:\\equiv|=|\u2261)\s*1"
                   r"|combined KO-dim(?:ension)?\s*=?\s*1\b"
                   r"|KO-dim(?:ension)?~?\s*1\s*[:(]\s*\(?\\?(?:varepsilon|epsilon)"
                   r"|at KO-dim(?:ension)?~?\s*1\b"
                   r"|epsilon_prime:\s*int\s*=\s*-1",
        "exempt_if_nearby": r"does not apply|NOT the additive|not KO-1|is \\emph\{not\} the KO-dimension|rather than|withdrawn|WITHDRAWN|relabel|graded product|GRADED product|KO-3",
        "severity": "fail",
        "scope": "group1 trunk",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "geovac/almost_commutative.py",
            "geovac/standard_model_triple.py",
            "tests/test_almost_commutative.py",
            "tests/test_standard_model_triple.py",
        ],
    },
    {
        "id": "pairdiag-composed-scaling-livesd",
        "note": "Trunk FULL 2026-09-01 (C9 finding).  The retired pair-diagonal composed-scaling claims -- O(Q^{2.5}) Pauli scaling and the 51x-1712x advantage range -- were live at ~20 loci across 7 documents, including a Paper 22 Corollary TITLED with the retired exponent, and survived the 2026-08-29 group3 certification.  Canonical: N_Pauli = 27.90 x Q exactly linear across molecules at fixed basis; equal-qubit advantage 54x-317x.",
        "pattern": r"Q\^\{2\.5\}\)?\$?[^.\n]{0,45}(?:Pauli|scal)"
                   r"|(?:Pauli|scal\w*)[^.\n]{0,45}Q\^\{2\.5\}"
                   r"|51\\times\$?\s*(?:to|--|-)\s*\$?1\{?,?\}?712"
                   r"|51\$?\s*(?:to|--|-)\s*\$?1\{?,?\}?712\s*\$?\\times",
        "exempt_if_nearby": r"retired|RETIRED|pair-diagonal (?:rule |count[s]? )?gave|superseded|withdrawn|corrected 2026-08|wrong-sign|under the retired",
        "severity": "fail",
        "scope": "group2 group3 group4 synthesis trunk",
        "files": [
            "papers/group3_foundations/paper_22_angular_sparsity.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex",
            "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex",
            "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex",
            "papers/group2_quantum_chemistry/paper_58_abelian_residue.tex",
            "papers/group4_quantum_computing/paper_14_qubit_encoding.tex",
            "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/synthesis/group4_quantum_computing_synthesis.tex",
            "papers/synthesis/geovac_field_guide.tex",
        ],
    },
    {
        "id": "pairdiag-density-pipeline-realizes",
        "note": "CF-1 zombie, retired 2026-08-30.  The pair-diagonal density "
                "D_pd was described as 'the density the composed pipeline "
                "realizes' in P22 tab:sparsity, the P22 Theorem-3 note, and "
                "the group3 synthesis.  It was never an intent: the "
                "production enumerator returned D_pd only because the "
                "wrong-sign-q ck_coefficient zeroed every m-changing "
                "multipole (fixed 2026-08-29).  The pipeline realizes the "
                "global-M_L D.",
        "pattern": r"density the composed pipeline realizes"
                   r"|D_\{\\mathrm\{pd\}\} realized by",
        "exempt_if_nearby": r"RETIRED|retired|withdrawn|WITHDRAWN"
                            r"|corrected 2026-08-30|Erratum",
        "severity": "fail",
        "scope": "group3",
        "files": [
            "papers/group3_foundations/paper_22_angular_sparsity.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
    },
    {
        "id": "pairdiag-density-as-the-angular-density",
        "note": "Retired 2026-08-30 (M9).  2.76% / 97.24% at l_max=2 is the "
                "stricter pair-diagonal D_pd, valid only for "
                "axially-symmetric or m-decoupled bases; the Coulomb "
                "(global-M_L) figures are 8.52% / 91.48%.  P26 quoted the "
                "former as 'the angular ERI density', overstating the "
                "symmetry-enforced sparsity 3.1x.",
        "pattern": r"angular ERI density is \$2\.76",
        "exempt_if_nearby": r"pair-diagonal|D_\{\\mathrm\{pd\}\}|stricter",
        "severity": "fail",
        "scope": "group6",
        "files": [
            "papers/group6_precision_observations/paper_26_entanglement.tex",
            "papers/group4_quantum_computing/paper_14_qubit_encoding.tex",
            "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex",
        ],
    },
    {
        "id": "bare-graph-n2-minus-1-attribution",
        "note": "Standing corpus tripwire (kappa_observation_not_derived; "
                "cert-3 item 10): -(n^2-1) is the CONTINUUM S^3 "
                "Laplace-Beltrami spectrum, which the graph encodes as "
                "labels; the bare graph Laplacian is a different, "
                "positive-semidefinite operator.  Text attributing the "
                "n^2-1 eigenvalues to the graph Laplacian itself may not "
                "re-surface without the continuum qualifier.",
        "pattern": r"eigenvalues are integers \(\$n\^2 ?- ?1\$"
                   r"|graph Laplacian are integers",
        "exempt_if_nearby": r"continuum|Laplace--?Beltrami|encodes as labels"
                            r"|positive-semidefinite|converge",
        "severity": "fail",
        "scope": "group6",
        "files": [
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "p27-ho-zero-entropy-rigidity",
        "note": "EP-2b RETRACTED v5.0.0: 'HO zero-entropy rigidity' may not appear "
                "as a live headline (the INDEX.md one-liner carried it un-flagged "
                "through three cert runs; cert-2 G11).",
        "pattern": r"HO zero-entropy rigidity|zero-entropy rigidity",
        "exempt_if_nearby": r"RETRACTED|retracted|withdrawn|WITHDRAWN|Erratum",
        "severity": "fail",
        "scope": "group6",
        "files": [
            "papers/INDEX.md",
            "papers/group6_precision_observations/paper_27_entropy_projection.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "p34-lamb-near-cancellation",
        "note": "Withdrawn 2026-08-28/29: the H Lamb Layer-2 inputs net to -1.06 MHz "
                "and do NOT cancel; 'netting to $-0.02$' and 'near-cancellation' "
                "may not re-surface un-flagged (cert-2 M3).",
        "pattern": r"netting to \$-0\.02|empirical near-cancellation",
        "exempt_if_nearby": r"withdrawn|WITHDRAWN|retired|corrected|Until 2026-08-28",
        "severity": "fail",
        "scope": "group6",
        "files": [
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
    },
    {
        "id": "p34-alkali-uniformity",
        "note": "Refuted 2026-08-28/29: the alkali undershoot is NOT a common factor "
                "('scales similarly for Na, K, Rb, Cs'); it grows monotonically "
                "-64.0% -> -98.4% (cert-2 M4).",
        "pattern": r"scales similarly for Na",
        "exempt_if_nearby": r"refuted|REFUTED|withdrawn|corrected|retired",
        "severity": "fail",
        "scope": "group6",
        "files": [
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
        ],
    },
    {
        "id": "p34-ee-eigen-closed-form",
        # A NEGATIVE result, added 2026-08-28 while extending the group6 C8
        # enumeration for the post-cert P34 remarks.  Negatives are the class most
        # prone to zombie back into positives (cf. the P24 rigidity corollary and
        # the Lorentzian "literal identification"), and this one is attractive:
        # "the coalescence object has a closed skeleton form" is a nicer sentence
        # than the truth.  Guarded pre-emptively rather than after the drift.
        "scope": "group6",
        "severity": "fail",
        "retired": "2026-08-27 (v5.1.5, rem:ee_partial_split sharpening): the "
                   "question 'do W's leading eigenvectors have closed skeleton "
                   "forms?' was answered NO, decisively and three ways. W(1) is "
                   "exactly rational but the active characteristic polynomials are "
                   "IRREDUCIBLE over Q (n_s=2 cubic 131072 l^3 - 50688 l^2 - 25200 l "
                   "- 675; quintic at n_s=3), so the finite eigenpairs are algebraic "
                   "numbers of full degree 2n-1 with no low-degree closed form; and "
                   "the Loewdin-compressed spectrum GROWS with n_s (1.03 -> 7.81 over "
                   "n_s = 2..8), so the finite eigenvectors do not converge to any "
                   "fixed continuum eigenfunctions. The ONLY surviving closed form is "
                   "DIFFERENTIAL: in u = 1/r the Coulomb kernel is min(u1,u2) (Green's "
                   "function of -d^2/du^2) and eigenfunctions of any weighted "
                   "compression solve lambda f'' = w f. A sentence asserting closed "
                   "FORMS for the eigenvectors/eigenvalues themselves is the zombie. "
                   "Pinned by tests/test_paper34_ee_split.py::"
                   "test_w_rational_anchor_and_irreducible_charpoly.",
        # 2026-09-07 (/qa paper_61 DELTA): both lookaheads were anchored
        # on `.*`, which scan_entry re-scans from every position of the
        # newline-joined file -- C16 did not terminate on Paper 34.
        # Bounded to the clause; discrimination unchanged, fire-tested.
        "pattern": "(?i)(?![^.\\n]{0,140}\\b(?:no|not|without|irreducible|lacks|lacking)\\b[^.\\n]{0,50}closed)(?![^.\\n]{0,140}closed[- ]form is differential)(?:eigen(?:vector|pair|value)s?[^.\\n]{0,140}?(?:closed[- ]form|closed skeleton form)|(?:closed[- ]form|closed skeleton form)[^.\\n]{0,140}?eigen(?:vector|pair|value)s?)",
        "exempt_if_nearby": "irreducible|no low-degree|no closed|differential|grows with|RETRACTED|withdrawn",
        "files": ["papers/group6_precision_observations/"
                  "paper_34_projection_taxonomy.tex"],
    },
    {
        "id": "p24-entanglement-rigidity",
        # spans group3 (P24), group4 (P23 resource counts + nuclear code) AND
        # group6 (P27 EP-2b retraction + the group6 synthesis) -- the 2026-08-24
        # retraction re-review found a LIVE zombie of this claim in the group6
        # synthesis that the original group3+group4 scope never checked.
        "scope": "group3 group4 group6",
        "severity": "fail",
        "retired": "2026-08-22 (FULL cert #2, code dimension): Paper 24's "
                   "two-fermion entanglement-rigidity corollary claimed the "
                   "closed-shell HO ground state has IDENTICALLY ZERO spatial "
                   "1-RDM entropy for any central V(r_12), on the mechanism "
                   "that a central potential preserves the total HO quantum "
                   "number N_tot. Both claim and mechanism are false. The "
                   "Moshinsky-Talmi BRACKET conserves N; the potential MATRIX "
                   "ELEMENT does not -- a central V(r_rel) is diagonal in the "
                   "CM quantum numbers and in l_rel but COUPLES different "
                   "relative-n, and n - n' = (N_bra - N_ket)/2 makes the "
                   "N-changing elements exactly the n != n' ones. They are "
                   "large: +17.3 MeV against a -0.55 MeV diagonal for the "
                   "Minnesota singlet. The published S = 0 came from an "
                   "undocumented `if N_bra != N_ket: return 0.0` guard in "
                   "geovac/nuclear/moshinsky.py. Corrected values: S = "
                   "0.0671 / 0.0716 / 0.0833 (spatial trace-1 1-RDM, "
                   "nats) and commutator ratio 0.74 / 0.63 / 0.67 "
                   "at N_max = 2 / 3 / 4. The paper's own text should have "
                   "exposed it -- it argued the Coulomb case is nontrivial "
                   "because 1/r_12 preserves no HO-like total-quanta number, "
                   "but 1/r_12 IS central. Retraction: Paper 24 "
                   "sec:entanglement-rigidity.",
        # The last two alternatives (added 2026-08-24) match the SYNTHESIS
        # phrasing the original pattern missed. They are written to fire on the
        # retracted TWO-BODY claim WITHOUT firing on the legitimate ONE-BODY
        # statement ("a non-degenerate ground state of a purely one-body
        # Hamiltonian has identically zero entropy"): "identically zero for"
        # (two-body: "... for any central two-body potential") never appears in
        # the one-body sentence ("... is identically zero (the ground state ...)"),
        # and "total-quanta conservation" is the withdrawn mechanism specifically.
        # `= 0(?![.\d])` = a BARE zero, so "S_full = 0.902" (a real nonzero value)
        # is not a hit; "entropy is identically zero for" is the two-body-entropy
        # claim specifically (does not match "identically zero for the diagonal
        # kinetic term" in P23, nor the one-body "is identically zero (the ...").
        "pattern": (r"zero\s+(von\s+Neumann\s+)?entanglement\s+entropy"
                    r"|entanglement\s+rigidity\s+(corollary|theorem)"
                    r"|S_?\{?\\?mathrm\{?full\}?\}?\s*=\s*0(?![.\d])"
                    r"|preserves\s+the\s+total\s+HO\s+quantum\s+number"
                    r"|entropy\s+is\s+identically\s+zero\s+for"
                    r"|total-quanta\s+conservation"
                    # the "zero-entropy rigidity"/"entropy-side rigidity/dual"
                    # phrasing the 2026-08-24 group6 FULL cert found live in the
                    # synthesis's "What is robust" list (does NOT match the VALID
                    # spectral "HO rigidity theorem" -- that phrase has no
                    # "zero-entropy"/"entropy-side" qualifier)
                    r"|zero-entropy\s+rigidity"
                    r"|entropy[-\s]side\s+(?:rigidity|dual)"),
        # NB: bare "artifact" removed 2026-08-24 -- it is too broad ("projection
        # artifact" is common non-retraction prose and false-exempted a live
        # zombie 4 lines away in the group6 synthesis); "guard" covers the real
        # P27 "undocumented-guard artifact" retraction discussion.
        "exempt_if_nearby": r"RETRACTED|retracted|withdrawn|guard|"
                            r"is false|do NOT commute|does not conserve",
        "files": [
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/group4_quantum_computing/paper_23_nuclear_shell.tex",
            "papers/group6_precision_observations/paper_27_entropy_projection.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
            "geovac/nuclear/*.py",
            "tests/test_paper24_ho_entropy.py",
        ],
    },
    {
        "id": "brown-surjection-attribution",
        "scope": "group3",
        "severity": "fail",
        "retired": "2026-08-22 (DELTA run, citation dimension): Paper 56 "
                   "asserted at two loci that Brown ESTABLISHES the "
                   "surjection U*_CM ->> G_MT(Z), and used it to flag a "
                   "'common framing slip' in other people's work. The "
                   "claim is in NEITHER candidate Brown paper: the cited "
                   "ICM survey (arXiv:1407.5165) contains no mention of "
                   "Connes, Marcolli, cosmic Galois, surjection or "
                   "renormalisation; and Brown's actual cosmic-Galois "
                   "paper (arXiv:1512.06409, section 0.5 Relation to "
                   "other work) explicitly declines the connection: 'It "
                   "is not clear if it is at all related to the groups "
                   "defined here.' The surjection is an elementary "
                   "comparison of two free pro-unipotent groups; the "
                   "paper now states it as its own observation, with the "
                   "two structural inputs sourced separately "
                   "(Connes-Marcolli 2004; Brown 2012 Ann. Math. 175).",
        "pattern": (r"Brown[^.\n]{0,60}establishes\s+a\s*"
                    r"(\\emph\{)?surjection"
                    r"|establishes\s+a\s*\n?\s*\\emph\{surjection\}"),
        "exempt_if_nearby": r"declines|not a theorem|elementary "
                            r"comparison|own observation|does \\emph\{not\}",
        "files": [
            "papers/group3_foundations/paper_56_tannakian_substrate.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
    },
    {
        "id": "propinquity-as-achieved-metric-group3",
        "scope": "group3",
        "severity": "fail",
        "retired": "2026-08-22 (FULL certifying run, synthesis "
                   "dimension): the group3 synthesis carried 'the GeoVac "
                   "propinquity convergence rate' and 'M1 governs "
                   "propinquity convergence rates', the label retired for "
                   "this result. What is established is van Suijlekom "
                   "STATE-SPACE Gromov-Hausdorff convergence (Paper 38, "
                   "unconditional via the translation-seminorm "
                   "metrization); Latremoliere propinquity is a different, "
                   "strictly stronger metric that is NOT achieved. The "
                   "same document already used the correct label "
                   "elsewhere. The existing propinquity entries were "
                   "scoped group1/group5 only, so no file list covered "
                   "group3 -- this entry closes that gap.",
        "pattern": r"(GeoVac|governs|the)\s+propinquity\s+convergence",
        "exempt_if_nearby": r"named\s+gap|state-space|strictly\s+stronger"
                            r"|NOT\s+the|historical|retract|not\s+achieved",
        "files": [
            "papers/group3_foundations/*.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
    },
    {
        "id": "latremoliere-propinquity-named-for-gh-rate",
        "scope": "trunk",
        "severity": "fail",
        "retired": "2026-09-01 (trunk FULL cert #2, registry-coverage "
                   "note): a planted S9 wrong-metric seed -- 'an "
                   "unconditional Latremoliere quantum-propinquity "
                   "convergence at rate (4/pi) log n/n' in the group3 "
                   "synthesis -- was caught by the LLM reviewer but NOT "
                   "by C16, for two reasons: the two existing propinquity "
                   "entries anchor on words BEFORE 'propinquity' "
                   "('the/GeoVac/governs', 'converge/established/proves' "
                   "+ Latr) and neither matches 'Latr...propinquity "
                   "convergence' with the rate AFTER it; and the group3 "
                   "entry exempts on 'state-space' within +-5 lines, which "
                   "the seeded sentence carried in its own next clause "
                   "(the scalar half correctly labelled). This entry "
                   "anchors on the Latremoliere name / the 'at rate' "
                   "clause and exempts only on an explicit denial. The "
                   "achieved result is van Suijlekom STATE-SPACE GH "
                   "convergence (Paper 38); Latremoliere propinquity is "
                   "strictly stronger and NOT achieved. Widened 2026-09-02 "
                   "(trunk Part F, F1.3): also fires on 'propinquity "
                   "rate(s)' / 'propinquity-rate' used as a name for the "
                   "paper's own state-space GH rate (six Paper 32 loci "
                   "renamed to 'state-space GH rate'). Widened again the "
                   "same day after the DELTA found 'propinquity asymptote' "
                   "at Paper 18:1120: rate(s) / asymptote(s) / constant(s).",
        "pattern": (r"Latr[^.\n]{0,30}(?:quantum-)?propinquity\s+"
                    r"convergence"
                    r"|(?:quantum-)?propinquity\s+convergence\s+at\s+rate"
                    r"|\bpropinquity[\s-]+(?:rates?|asymptotes?|constants?)\b"),
        "exempt_if_nearby": r"no\s+published|not\s+achieved|strictly\s+"
                            r"stronger|NOT\s+the|named\s+gap|retract"
                            r"|historical|different\s+metric|is\s+not\s+"
                            r"Latr|not\s+a\s+Latr",
        "files": [
            "geovac/gh_convergence.py",
            "tests/test_gh_convergence.py",
            "papers/group3_foundations/*.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group1_operator_algebras/"
            "paper_38_su2_propinquity_convergence.tex",
            "papers/group1_operator_algebras/"
            "paper_40_unified_propinquity_convergence.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
        ],
    },
    {
        "id": "p45-kplus-compression-theorem-live",
        "scope": "trunk",
        "severity": "fail",
        "retired": "2026-09-01 (trunk FULL cert #2, registry-coverage "
                   "note): Paper 45's K^+-compression theorem was "
                   "WITHDRAWN 2026-06-09 (falsifier "
                   "tests/test_p45_kplus_degeneracy.py): the K^+ "
                   "compression annihilates the Lipschitz seminorm, a "
                   "degeneracy theorem, not a Lorentzian convergence; "
                   "Lorentzian quantum-metric convergence is DESCOPED "
                   "(Papers 45-49 Status notes). Two planted S8 zombie "
                   "seeds citing the theorem as established (group3 "
                   "synthesis; Paper 32 'Continuum scope' paragraph) "
                   "were caught by the LLM reviewers but by NO registry "
                   "entry -- the krein entry guards the literal Krein "
                   "identification and exempts on 'K^+', the propinquity "
                   "entries never mention the theorem. The legitimate "
                   "wordings are 'K^+-compression DEGENERACY theorem' "
                   "and 'descoped'; a bare 'K^+-compression theorem' or "
                   "a live 'Lorentzian quantum-metric convergence' is "
                   "the zombie.",
        "pattern": (r"(?:K\^\{?\+\}?\$?|\\Kplus\$?)-?\s*compression\s+"
                    r"theorem"
                    r"|Lorentzian\s+quantum[- ]metric\s+convergence"),
        "exempt_if_nearby": r"degenerac|descope|withdrawn|retract"
                            r"|annihilat|intended\s+to\s+assert"
                            r"|open\s+question|convergence\s+question"
                            r"|open\s+named|not\s+a\s+Lorentzian",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex",
            "papers/group1_operator_algebras/"
            "paper_38_su2_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_4[2-9]_*.tex",
            "papers/group1_operator_algebras/paper_5[0-3]_*.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/synthesis/geovac_field_guide.tex",
        ],
    },
    {
        "id": "emn-catalan-negative-attribution",
        "scope": "group3",
        "severity": "fail",
        "retired": "2026-08-22 (FULL certifying run, citation dimension): "
                   "Paper 56 asserted in four places that "
                   "Eskandari-Murty-Nemoto 2025 (arXiv:2510.20648) PROVE "
                   "Catalan G is NOT a period of mixed Tate motives over Q "
                   "(or MT(Z)), and used that negative as the FORCING "
                   "argument for adopting G_4 over Brown's G_MT(Z). The "
                   "source establishes only the POSITIVE half. Its abstract "
                   "reads in full: 'We first give a geometric construction "
                   "of a 2-dimensional mixed motive over Q with the Catalan "
                   "constant G as a period. We then use this motive to "
                   "obtain a supply of linear forms in 1 and G. We also "
                   "explicitly compute the coefficients of 1 and G in these "
                   "linear forms.' No negative result is claimed; "
                   "non-membership in MT(Q) is expected but open. The paper "
                   "now says the level-4 choice is MOTIVATED by where G is "
                   "known to live, not FORCED by a proven exclusion.",
        "pattern": r"(prove[sd]?|provabl[ey])[^.]{0,80}not a period of "
                   r"mixed Tate|prove[sd]? Catalan \$?G\$?[^.]{0,40}is not "
                   r"a\s*\n?\s*period",
        "exempt_if_nearby": r"expected but|open|conjectur|not proven"
                            r"|to our knowledge",
        "files": [
            "papers/group3_foundations/paper_56_tannakian_substrate.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
    },
    {
        "id": "lorentzian-literal-identification-krein",
        "scope": "group3 group6 trunk group1",
        "severity": "fail",
        "retired": "2026-07-04 (group6 first-cert run): P34 III.29's pre-descope "
                   "Lorentzian claims -- 'literal identification at the Krein "
                   "operator-system level (finite cutoff)' and 'genuine Lorentzian "
                   "extension of Paper 42' (Sprint L2-E, 2026-05-17) -- are WITHDRAWN. "
                   "The 2026-06-09 P45 K+ annihilation theorem + the 2026-06-19 "
                   "compact-boost closure show the truncated BW boost is compact "
                   "(integer spectrum, e^{2 pi i K}=I), so the period closure is the "
                   "compact KMS beta=2pi circle and the Lorentzian signature is "
                   "metrically invisible at finite cutoff (Euclidean/convention). "
                   "NOTE: the Riemannian operator-system-level closure is NOT retired. "
                   "Completeness-critic catch. SCOPE WIDENED 2026-08-24 (group3 "
                   "re-cert): the SAME withdrawn reading survived in Paper 31 §9 "
                   "(sec:sig_l2_verification) -- 'as literal identification ... not "
                   "just structural correspondence' and 'the Lorentzian closure is "
                   "complete' -- because the entry was group6-scoped and P31 (group3) "
                   "was never gated. Added group3 scope + P31 file + two single-line "
                   "tells of the strong reading.",
        "pattern": r"genuine\s+Lorentzian\s+(?:\\emph\{)?extension"
                   r"|literal\s+identification\s+at\s+the\s+Krein"
                   r"|Krein-level\s+four-witness\s+Wick-rotation\s+theorem\s+closes"
                   r"|not\s+just\s+structural\s+correspondence"
                   r"|Lorentzian\s+closure\s+is\s+complete"
                   r"|literal\s+identification\s+at\s+the\s+operator-system\s+level\s*\(Lorentzian"
                   r"|the\s+Lorentzian\s+\\emph\{extension\}\s+of\s+Paper",
        "exempt_if_nearby": r"withdrawn|WITHDRAWN|signature-blind|compact[- ]boost"
                            r"|compact\s+KMS|K\^?\+|descope|convention|period[- ]closure"
                            r"|Euclidean|not\s+constitute",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex",
        ],
    },
    {
        "id": "fock-coupling-one-sixteenth-prefactor",
        "scope": "trunk group3 group5 synthesis",
        "severity": "fail",
        "retired": "2026-09-03 (trunk FULL run #3, I.1.1): the inter-shell coupling "
                   "c^2(n,l) = |<n+1,l|cos chi|n,l>|^2 is (1/4)[1 - l(l+1)/(n(n+1))] "
                   "(Chebyshev amplitude 1/2), not (1/16)[...]; the 1/16 is the inverse "
                   "Fock Jacobian, a different quantity.  c^2(4,3) = 1/10; Delta = 1/40 "
                   "is the composite (2/5)/Omega^4(0).",
        "pattern": r"c\^2\(n,\s*l\)\s*(?:\\;\s*)?(?::)?=\s*(?:\\;\s*)?\\t?frac\{1\}\{16\}"
                   r"|\(1/16\)\[1\s*-\s*l\(l\+1\)"
                   r"|c\^2\(n,\s*0\)\s*=\s*1/16"
                   r"|c\^2\(4,\s*3\)\s*=\s*(?:\\t?frac\{1\}\{16\}|1/40)"
                   r"|squared\s+Chebyshev\s+transition\s+amplitude\s+\$?\(1/4\)\^2",
        "exempt_if_nearby": r"corrected 2026-09-03|retired|withdrawn|different quantity|printed before|earlier version|before 2026-09-03",
        "files": [
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/group5_qed_gauge/paper_2_alpha.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "tests/test_trunk_qa_kappa.py",
            "tests/test_trunk_qa_c2_delta.py",
        ],
    },
    {
        "id": "p39-chirality-grading-on-s3",
        "scope": "trunk group1",
        "severity": "fail",
        "retired": "2026-09-04 (DELTA #7, sub-agent finding F1, verified symbolically "
                   "in the parent session).  Paper 39 described gamma_a as 'the "
                   "chirality grading on L^2(S^3, Sigma)' and asserted "
                   "{D_CH, gamma_a} = 0.  No such operator exists: a grading "
                   "anticommuting with every generator of Cl(3) also anticommutes "
                   "with the volume element omega = sigma1 sigma2 sigma3 = i I, which "
                   "in odd dimension is CENTRAL, forcing gamma = 0.  Paper 38 ('KO-3 "
                   "carries no chirality', 2026-09-03) and Paper 32 (2026-09-02) "
                   "already said so;  the correction never reached Paper 39.  The "
                   "displayed D_{a,b} is the even-x-odd product formula while the "
                   "paper's KO arithmetic 3+3=6 is the odd-x-odd case.  Repair "
                   "(Clifford doubling) is NAMED, NOT ADOPTED: whether L3-T and L4-T "
                   "hold for it is open.",
        "pattern": r"chirality grading on \$?L\^2\(\\?sthree"
                   r"|gamma_a\$? is the chirality grading"
                   r"|\\\{\\DCH\^\{\(a\)\}, \\gamma_a\\\} = 0[^\n]{0,40}$"
                   r"|satisfies the KO-dim-6 chirality and real-structure axioms",
        "exempt_if_nearby": r"(?!)",
        "cited_by": {
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex":
                "repair adopted 2026-09-05 (lem:real_structure-T, Clifford doubling)",
            "geovac/gh_convergence_tensor.py":
                "disclosed 2026-09-04 (DELTA #7) at 4 loci; repair open",
        },
        "files": [
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
            "geovac/gh_convergence_tensor.py",
        ],
    },
    {
        "id": "p39-lambda-placement",
        "scope": "trunk group1",
        "severity": "fail",
        "retired": "2026-09-04 (DELTA #7, sub-agent finding F2, verified two ways in "
                   "the parent session).  The focal-length placement was inverted: "
                   "gamma/lambda where it should be lambda*gamma.  D -> lambda^-1 D "
                   "divides the Lipschitz seminorm by lambda, so the MK unit ball -- "
                   "and every distance -- GROWS by lambda;  gamma is a first moment of "
                   "the geodesic distance.  Internal witness: Paper 38 S2.1's "
                   "dual-Coxeter Dirac is half the CH Dirac and its moment is TWICE "
                   "the unit-S^3 one.  gh_convergence_tensor stated the error "
                   "explicitly, conflating the seminorm's scaling with the moment's.",
        "pattern": r"\\frac\{\\gamma_\{\\nmax[ab]\}\}\{\\lambda_[ab]\}"
                   r"|gamma_\{\\nmax[ab]\}/\\lambda_[ab]"
                   r"|gamma_\{n_[ab]\}\s*/\s*lambda_[ab]"
                   r"|gamma_[ab]\s*/\s*lambda_[ab]"
                   r"|lambda_[ab]\^\{-1\}\s*gamma",
        "exempt_if_nearby": r"(?!)",
        "cited_by": {
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex":
                "corrected 2026-09-04 (DELTA #7)",
            "geovac/gh_convergence_tensor.py": "corrected 2026-09-04 (DELTA #7)",
        },
        "files": [
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
            "geovac/gh_convergence_tensor.py",
        ],
    },
    {
        "id": "sp-ratio-false-precision",
        "scope": "trunk group3",
        "severity": "fail",
        "retired": "2026-09-04 (trunk DELTA #7, CLAIMS-B headline).  Two successive "
                   "statements of the off-branch/on-branch s/p ratio, both stated as "
                   "theorems from three sampled cutoffs.  'Three to six times larger' "
                   "(retired 2026-09-04 morning; ratio is unbounded, ~n/(sqrt3 pi)) and "
                   "its replacement 'tracking n/(sqrt3 pi) to within 1%' -- false at two "
                   "of its own three printed anchors (19.7% high at 30, 4.6% at 120), "
                   "silently using the n-1 neighbour (the n-2 neighbour gives 4.2x at 30, "
                   "not 6.6x), and 'off-branch is larger' inverts below n ~ 12 (0.65% at "
                   "7 vs 12.7% at 6).  The guard-asymptotics failure, twice, in the same "
                   "sentence.",
        "pattern": r"three to six times"
                   r"|to within 1\\%[^\n]{0,80}n_\{\\max\}/\(\\sqrt"
                   r"|tracking n_\{\\max\}/\(\\sqrt\{3\}[^\n]{0,20}to within"
                   r"|within \$?1\\%\$?:\\? *\$?6\.6",
        "exempt_if_nearby": r"(?!)",
        "cited_by": {
            "papers/group3_foundations/paper_1_spectrum.tex": "remediated 2026-09-04 (DELTA #7); review owed (DELTA #8)",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "remediated 2026-09-04 (DELTA #7); review owed (DELTA #8)",
        },
        "files": [
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
    },
    {
        "id": "p53-gradient-non-expansive-height",
        "scope": "group1",
        "severity": "fail",
        "retired": "2026-09-07 (/qa group1, PI-flagged height leg). Paper 53 "
                   "step (iii) claimed the plane Bochner-Riesz Berezin is "
                   "GRADIENT-NON-EXPANSIVE, ||grad B f|| <= ||grad f||, "
                   "'verified numerically, ratio <= 1 at every Lambda, rising "
                   "to 1 as Lambda -> infinity'. FALSE over the unit-Lipschitz "
                   "ball: B is a radial convolution, so the sharp constant is "
                   "the kernel L1 norm, which is > 1 at every finite Cesaro "
                   "order (5.72 at s=0.6, 2.01 at s=1, 1.23 at s=2) and is "
                   "LAMBDA-INDEPENDENT by scaling -- so the reported "
                   "Lambda-dependence was a property of the test functions. "
                   "NOT Paper 38 L5 transported (different quantity; the "
                   "finite-band witness satisfies non-expansiveness rather "
                   "than breaking it). Consequence for the assembly is OPEN "
                   "(PI adjudication); the reach leg is untouched.",
        "pattern": r"gradient-non-expansive"
                   r"|ratio \$?\\le\s*1\$? at every \$?\\Lambda"
                   r"|rising to \$?1\$? as \$?\\Lambda\\to\\infty",
        "exempt_if_nearby": r"NOT|not\b|false|corrected 2026-09-07|Lebesgue"
                            r"|rem:height_constant|>\s*1|exceeds",
        "cited_by": {
            "papers/group1_operator_algebras/paper_53_disk_propinquity.tex":
                "reviewed 2026-09-07 -- step (iii) restated, rem:height_constant added",
            "tests/test_paper53_height_constant.py":
                "reviewed 2026-09-07 -- new, fire-tested both directions",
        },
        "files": [
            "papers/group1_operator_algebras/paper_53_disk_propinquity.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "tests/test_paper53_height_constant.py",
        ],
    },
    {
        "id": "p40-pythagorean-cross-manifold-future-work",
        "scope": "group1 synthesis",
        "severity": "fail",
        "retired": "2026-09-07 (/qa group1 DELTA, carryforward U.2). Paper 40 "
                   "sec:cross_manifold proposed extending to G != H 'with the "
                   "Pythagorean Leibniz constant generalising to a Cartan-"
                   "product-type bound'. Paper 39 WITHDREW that refinement "
                   "(graded anticommutation does not give operator-norm "
                   "orthogonality; Remark rem:no_pythagorean) and proves its "
                   "theorem by the lifted-state assembly instead. Future work "
                   "aimed at an abandoned route. The C_3^(2) <= sqrt2 constant "
                   "survives with its value unchanged but is NOT a rate "
                   "constant of Paper 39's theorem.",
        "pattern": r"Pythagorean Leibniz\s*\n?constant generalising"
                   r"|extend mechanically from Paper~?39 with the Pythagorean",
        "exempt_if_nearby": r"withdrawn|WITHDRAWN|abandoned|corrected 2026-09-07"
                            r"|lifted-state",
        "cited_by": {
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex":
                "reviewed 2026-09-07 -- paragraph redirected to the lifted-state route",
        },
        "files": [
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
        ],
    },
    {
        "id": "k2-proved-case-carries-sqrt2-constant",
        "scope": "group1 synthesis",
        "severity": "fail",
        "retired": "2026-09-07 (/qa group1 DELTA, carryforward U.2 'Upgrade'). "
                   "The group1 synthesis presented the PROVED k=2 case as "
                   "carrying the sqrt(k) triangle-bound constant. Paper 39 "
                   "states that C_3^(2) <= sqrt2 'belongs to the ABANDONED "
                   "(B,P)-pair route and is not a rate constant of this "
                   "theorem; its value is unchanged, only its role'. NOTE THE "
                   "DIRECTION: this is the MIRROR of ordinary staleness -- the "
                   "owner STRENGTHENED its claim (lifted-state, no Berezin map, "
                   "no partial inverse) and the summary kept the weaker form. "
                   "C16 and cited_by cannot catch that class in general, "
                   "because nothing is retracted and nothing is misspelt; this "
                   "entry guards only THIS occurrence, now that its wording is "
                   "known. The class stays open (CLAUDE.md S9, mirror "
                   "direction).",
        # `\$?` around k=2 and `\s*` across the wrap: the synthesis writes
        # "Only the $k=2$ case is\nproved (Paper~39)", and the first draft of
        # this pattern (no math delimiters) matched none of it.
        # `[^.\n]` NOT `[^\n]`: C16 also scans text with newlines replaced by
        # spaces, where `[^\n]{0,60}` stops bounding anything and can cross a
        # paragraph break.  It did -- matching "$k=2$ case is proved; see
        # below)." plus the next heading's "$\sqrt{k}$" 52 chars later.
        "pattern": r"\$?k\s*=\s*2\$?\s+case\s+is\s+proved[^.\n]{0,60}sqrt"
                   r"|proved \$?k\s*=\s*2\$? case[^.\n]{0,60}sqrt"
                   r"|\$?k\s*=\s*2\$? case carries the[^.\n]{0,30}sqrt",
        "exempt_if_nearby": r"abandoned|not a rate constant|Upgrade|only\s+its\s+role"
                            r"|lifted-state|corrected 2026-09-07",
        "cited_by": {
            "papers/synthesis/group1_operator_algebras_synthesis.tex":
                "reviewed 2026-09-07 -- k-fold paragraph now records the constant's ROLE",
        },
        "files": [
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
        ],
    },
    {
        "id": "p61-w0-transcendence-cancels",
        "scope": "paper_61 group3",
        "severity": "fail",
        "retired": "2026-09-07 (/qa paper_61 FULL). Paper 61 asserted that "
                   "'the transcendence of the individual masters cancels in "
                   "their determinant'. FALSE. Abel's identity fixes only the "
                   "D-dependence; the constant W_0 was never computed in the "
                   "paper or its backing test, which normalises it to 1 by "
                   "construction -- so no artifact in the corpus could have "
                   "seen the cancellation fail. It is now computed from the "
                   "branch-point data: W_0 = pi^2/rho^2 (Paper 61 eq:W0, "
                   "[SYMBOLIC]), verified symbolically in sympy step-by-step "
                   "and numerically at rho = 0.3/0.5/0.71 to 40 digits. So "
                   "pi^2 does not cancel -- it is SQUARED -- and W_0 is "
                   "rho-dependent. The withdrawal is a net gain: an exact "
                   "value replaced a false cancellation.",
        "pattern": r"transcendence[^.\n]{0,70}(?:individual\s+)?masters?"
                   r"[^.\n]{0,40}cancels"
                   r"|cancels\s+in\s+their\s+determinant"
                   r"|transcendence\s+cancels",
        "exempt_if_nearby": r"(?!)",   # standardized marker only
        "cited_by": {
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex":
                "reviewed 2026-09-07 -- withdrawn in-paper and REPLACED by the "
                "closed form eq:W0 = pi^2/rho^2, tier [SYMBOLIC]; the quoted "
                "retired sentence in the correction note carries the marker",
            "docs/qa/paper_61.done.md":
                "reviewed 2026-09-07 -- W3 row records the withdrawal; carries "
                "the marker where it quotes the retired sentence",
        },
        "files": [
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex",
            "papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
            "docs/qa/paper_61.done.md",
            "docs/claim_test_matrix.md",
            "tests/test_paper59_bessel_moment_algebra.py",
        ],
    },
    {
        "id": "p61-every-cm-fibre-universal",
        "scope": "paper_61 group3 group2",
        "severity": "fail",
        "retired": "2026-09-07 (/qa paper_61 FULL, claims dimension). The "
                   "corpus repeatedly claimed that sweeping the physical base "
                   "rho over (0, inf) makes the family pass through EVERY CM "
                   "fibre. FALSE: lambda = 1 - rho, so the sweep covers only "
                   "the REAL locus lambda < 1 of X(2). Fibres with non-real "
                   "lambda -- discriminant -3 among them -- are off the "
                   "contour, as is lambda = 2. Correct form: INFINITELY MANY "
                   "CM fibres, not every one. This mattered because the "
                   "universal was used deflationarily (to argue tau = i is not "
                   "special), so it read as harmless while being false. "
                   "SURVIVED THE FIRST REMEDIATION at Paper 56:1848, an ACTIVE "
                   "paper citing Paper 61 fourteen lines above -- the measured "
                   "reason this entry exists.",
        "pattern": r"(?:every|all)\s+(?:of\s+)?(?:the\s+)?CM[\s_-]*fib"
                   r"|CM[\s_-]*fib\w*\s+is\s+hit"
                   r"|hits?\s+EVERY\s+CM",
        "exempt_if_nearby": r"(?!)",   # standardized marker only
        "cited_by": {
            "papers/group3_foundations/paper_56_tannakian_substrate.tex":
                "reviewed 2026-09-07 -- :1848 corrected to 'infinitely many "
                "(not every one)'; its paragraph is deflationary so the "
                "conclusion is unchanged",
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex":
                "reviewed 2026-09-07 -- owner; states the exclusion explicitly "
                "at L173-177 (non-real lambda, disc -3, lambda=2)",
            "memory/cosmic_galois_elliptic_rung1.md":
                "reviewed 2026-09-07 -- corrected; this file is loaded into "
                "EVERY session, so a false universal here reseeds itself",
            "debug/routeC_cosmic_galois_rung2.py":
                "reviewed 2026-09-07 -- docstring and stdout both corrected "
                "(live driver named by claim-matrix row 496)",
        },
        "files": [
            "papers/group3_foundations/paper_56_tannakian_substrate.tex",
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex",
            "papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "memory/cosmic_galois_elliptic_rung1.md",
            "debug/routeC_cosmic_galois_rung2.py",
            "docs/qa/paper_61.done.md",
            "docs/claim_test_matrix.md",
        ],
    },
    {
        "id": "p61-galois-in-sp4-z",
        "scope": "paper_61 group3",
        "severity": "fail",
        "retired": "2026-09-07 (/qa paper_61 FULL, TWO reviewers converging "
                   "independently). Paper 61 said the DIFFERENTIAL GALOIS "
                   "group of the rank-4 Picard-Fuchs operator lies in "
                   "Sp_4(Z), in its abstract, introduction and body, plus the "
                   "test docstring and claim-matrix row 499. FALSE: Sp_4(Z) is "
                   "DISCRETE, so a Zariski-closed subgroup of GL_4(C) inside "
                   "it is finite, forcing every solution algebraic -- "
                   "contradicting the connection's irregularity and the "
                   "exponential torus (C*)^2 that Paper 59 establishes in the "
                   "local Galois group (Ramis density). The backing variable "
                   "is literally `_L4_MONODROMY` and the test asserts "
                   "M0^T Omega M0 = Omega, a MONODROMY statement. Paper 59:585 "
                   "states the Galois containment correctly as Sp_4(C). "
                   "CORRECT FORM: monodromy in Sp_4(Z), Galois in Sp_4(C). "
                   "Nothing numeric changes -- B = pi*Omega and the integral "
                   "structure stand.",
        "pattern": r"differential Galois group[^.\n]{0,80}Sp[_ ]?\{?4\}?"
                   r"[^.\n]{0,20}\\?mathbb\{?Z"
                   r"|Galois group into \$?\\mathrm\{Sp\}_\{4\}\(\\mathbb\{Z\}\)"
                   r"|forcing the differential Galois group into"
                   # bare ASCII spelling used across docs/ and debug/, reachable
                   # only after the 2026-09-07 Unicode folding (Sp4(Z) <- Sp_4(Z)).
                   # The `\(\s*Z` is load-bearing: WITHOUT it this alternative
                   # fired on test_routeC_momentum.py:492/:504, "self-adjoint =>
                   # differential Galois group in Sp4", which is TRUE over C.
                   # Only Sp4(Z) is the retired claim.
                   r"|Galois[^.\n]{0,40}Sp[_ ]?\{?4\}?\s*\(\s*(?:\\?mathbb\{?)?Z",
        "exempt_if_nearby": r"monodromy|FALSE|false|Correction|corrected 2026-09-07"
                            r"|Sp_\{4\}\(\\mathbb\{C\}\)|discrete|withdrawn",
        "cited_by": {
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex":
                "reviewed 2026-09-07 -- all three loci corrected, with an "
                "in-paper correction note giving the discreteness argument",
            "docs/claim_test_matrix.md":
                "reviewed 2026-09-07 -- row 499 scanned; Galois/monodromy "
                "wording corrected where present",
            "tests/test_routeC_momentum.py":
                "reviewed 2026-09-07 -- :900 docstring and :907 inline comment "
                "corrected to MONODROMY; :504 left as-is, self-adjointness "
                "genuinely does put the Galois group in Sp4 over C",
        },
        "files": [
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex",
            "papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/INDEX.md",
            "docs/claim_test_matrix.md",
            "docs/qa/paper_59.done.md",
            "docs/qa/paper_61.done.md",
            "tests/test_routeC_momentum.py",
            "debug/routeC_intersection_form.py",
        ],
    },
    {
        "id": "p61-broadhurst-mellit-quadratic",
        "scope": "paper_61 group3 group2",
        "severity": "fail",
        "retired": "2026-09-07 (/qa paper_61 FULL, citation dimension). The "
                   "quadratic relations between Bessel moments are the "
                   "Broadhurst-ROBERTS relations (proved by Fresan-Sabbah-Yu "
                   "2023 and independently by Zhou, CNTP 15(4) 651-741 (2021), "
                   "arXiv:2012.03523). Broadhurst-MELLIT names the DETERMINANT "
                   "formulae -- Zhou's cited paper is literally titled "
                   "'Wronskian factorizations and Broadhurst-Mellit determinant "
                   "formulae', while Paper 61's body cited that same work for "
                   "'Broadhurst-Mellit quadratic period relations' at four loci "
                   "INCLUDING THE ABSTRACT. The bibliography contradicted the "
                   "text. No literature usage of 'Broadhurst-Mellit quadratic "
                   "relations' exists. "
                   "REBUILT 2026-09-07 (same-day DELTA): v1 matched ASCII "
                   "hyphens only and missed every en-dashed locus in docs/; it "
                   "also exempted on the word 'determinant' within +-5 LINES, "
                   "so claim_test_matrix row 498's legitimate determinant "
                   "mention sheltered row 499's live defect one line below. "
                   "The determinant sense is CORRECT usage, so it is now "
                   "discriminated in the PATTERN (tempered: the match dies at "
                   "'determinant'/'Wronskian') rather than by a line window.",
        # Tempered on both sides: Broadhurst-Mellit is legitimate for the
        # DETERMINANT formulae, so a span reaching `determinant`/`Wronskian`
        # before it reaches `quadratic`/`relation` is correct usage and must
        # not fire.  Corrective mentions carry the standardized marker.
        # Linear-cost tempering (2026-09-07): the exclusion is asserted ONCE
        # by a bounded lookahead instead of per-character, because the gate
        # also scans the file newline-joined and the per-character form
        # (?:(?!X)[^.\n]){0,200} backtracked without terminating there.
        "pattern": r"(?<![-\w])Broadhurst--?Mellit"
                   r"(?![^.\n]{0,90}(?:determinant|Wronskian|Wro\\?'?nskian))"
                   r"[^.\n]{0,90}(?:quadratic|period relation|relations between)"
                   r"|(?:quadratic(?:\s+period)?\s+(?:relation|pairing)"
                   r"|period\s+relation)"
                   r"(?![^.\n]{0,200}(?:determinant|Wronskian|Wro\\?'?nskian))"
                   r"[^.\n]{0,200}(?<![-\w])Broadhurst--?Mellit",
        "exempt_if_nearby": r"(?!)",   # standardized marker only, per the
        # 2026-09-03 rationale: authored exemption vocabulary drawn from the
        # surrounding correct text is the documented failure mode, and v1 of
        # THIS entry is a measured instance of it.
        "cited_by": {
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex":
                "reviewed 2026-09-07 -- 4 loci renamed, zhou_quadratic2021 "
                "bibitem added; the 3 remaining Mellit mentions are the "
                "determinant sense and are correct",
            "papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex":
                "reviewed 2026-09-07 -- only locus is bibitem :1017, Zhou's "
                "actual title (determinant sense); correct, left as-is",
            "docs/claim_test_matrix.md":
                "reviewed 2026-09-07 -- row 499 renamed to Broadhurst-Roberts; "
                "row 498 is the determinant sense and is correct",
            "docs/qa/paper_59.done.md":
                "reviewed 2026-09-07 -- :192/:252 renamed; :419 NIT marked "
                "APPLIED (it had named this defect on 2026-08-21 and sat "
                "unapplied at the 2026-09-07 run)",
            "tests/test_paper59_bessel_moment_algebra.py":
                "reviewed 2026-09-07 -- :3/:55/:120 renamed; :25 is the "
                "determinant sense and is correct",
            "debug/routeC_bessel_moment_algebra.py":
                "reviewed 2026-09-07 -- :3/:111 renamed (live driver); "
                ":16/:50 determinant sense, correct",
        },
        "files": [
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex",
            "papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
            "papers/INDEX.md",
            "docs/claim_test_matrix.md",
            "docs/qa/paper_59.done.md",
            "docs/qa/paper_61.done.md",
            "tests/test_paper59_bessel_moment_algebra.py",
            "debug/routeC_bessel_moment_algebra.py",
        ],
    },
    {
        "id": "p60-sublinear-as-regime",
        "scope": "paper_60 group2 synthesis",
        "severity": "fail",
        "retired": "2026-09-07 (/qa 58/59/60 FULL, code L1).  Paper 60 stated "
                   "||M||_1 ~ K^0.84 as a REGIME ('grows sublinearly') and "
                   "attributed the mechanism to T', 'the pure-number matrix "
                   "elements decrease with quantum number'.  Re-measured with the "
                   "paper's own gen_configs/solve past its largest fitted point: "
                   "0.84 is a fit over the WINDOW K = 74..164 (0.8399), the local "
                   "slope rises monotonically outside it (0.850, 0.868, 0.882, "
                   "0.906 at K = 202, 244, 290, 340), and the split is the "
                   "opposite of the printed mechanism -- nuclear diagonal T^0 = "
                   "Z R_nu at K^0.70 (stable), pure-number block T' at K^1.05 "
                   "(SUPERlinear).  The paper's own Sec. 7 already said the "
                   "advantage 'rides on the clean diagonal T^0', so Sec. 4 "
                   "contradicted Sec. 7 inside one document.  WHAT SURVIVES: the "
                   "1-norm does grow more slowly than the matrix dimension over "
                   "every computable basis, and the contrast with the L2 "
                   "superlinear inflation is real -- the encoding claim stands, "
                   "only its asymptotic reading and its mechanism were wrong.",
        # (a) the regime reading; (b) the backwards mechanism.  Both are
        # wording-tolerant: the synthesis restated each in its own words.
        "pattern": r"grows \\emph\{sublinearly\}"
                   r"|with a sublinear\s+\n?\$1\$-norm"
                   r"|\$1\$-norm grows sublinearly"
                   r"|1-norm grows sublinearly"
                   r"|pure-number matrix elements decrease\s*\n?with quantum number"
                   r"|large-K~?164 asymptote",
        # The corrected text carries the window or the diagonal attribution.
        "exempt_if_nearby": r"K\s*=\s*74|window|over the computed|more slowly than "
                            r"the (?:matrix dimension|configuration count)|nuclear "
                            r"diagonal|T\^0|NOT an asymptote|corrected 2026-09-07",
        # Documents whose ARGUMENT rests on the sublinearity claim.
        "cited_by": {
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex":
                "reviewed 2026-09-07 -- block rewritten with window + T^0/T' split",
            "tests/test_paper60_sturmian.py":
                "reviewed 2026-09-07 -- 'asymptote' comment corrected; new "
                "test_paper60_sublinearity_is_carried_by_the_nuclear_diagonal "
                "fire-tested in both directions",
            "docs/claim_test_matrix.md":
                "reviewed 2026-09-07 -- row re-tiered to the windowed claim",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/INDEX.md",
            "docs/claim_test_matrix.md",
            "tests/test_paper60_sturmian.py",
            "geovac/sturmian_secular.py",
            "geovac/sturmian_molecular_lambda.py",
        ],
    },
    {
        "id": "graph-spectrum-attribution",
        "scope": "trunk group3 synthesis",
        "severity": "fail",
        "retired": "2026-09-04 (trunk DELTA #7, SYNTH M1/M2).  The graph Laplacian's "
                   "spectrum is closed form and bounded in [0, 8] (Paper 0 S VI); it "
                   "does NOT approach n^2-1, and 'integer eigenvalues on the unit S^3' "
                   "is a statement about the CONTINUUM operator, not the discrete graph "
                   "(Paper 7:15, :77).  The group3 synthesis said both -- 'the spectrum "
                   "of (D-A) approaches ... n^2-1' and listed integer eigenvalues under "
                   "'statements about the discrete graph' in its abstract and "
                   "conclusion -- contradicting its own S1.5-compliant sentence at :434.",
        "pattern": r"spectrum of \$?\(D\s*-\s*A\)\$? approaches"
                   r"|approaches the continuum \$?S\^3\$? Laplace--Beltrami spectrum"
                   r"|spectrum-generating\s+operator"
                   r"|integer\s+eigenvalues on the unit \$?S\^3\$?",
        "exempt_if_nearby": r"(?!)",
        "cited_by": {
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated 2026-09-04 (DELTA #7); review owed (DELTA #8)",
        },
        "files": [
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex",
        ],
    },
    {
        "id": "saturation-approach-monotone",
        "scope": "trunk group3",
        "severity": "fail",
        "retired": "2026-09-04 (trunk DELTA #6 remediation): the approach of the "
                   "l-block saturation constant to C = 42.7397... is NOT monotone -- "
                   "195 decreasing steps over n = 10..600, e.g. C_20 = 40.7285 > "
                   "C_21 = 40.6470.  Every finite sample understates C because the "
                   "estimate is a ONE-SIDED BOUND, not because the sequence rises. "
                   "Paper 0 corrected this 2026-09-04;  Paper 7 was still arguing "
                   "from monotonicity a day later, which is what put it here.",
        "pattern": r"approach\s+is\s+monotone"
                   r"|monotone\s+from\s+below"
                   r"|monotonically\s+approach(?:es|ing)?\s+\$?C\$?\b"
                   r"|rises\s+monotonically\s+to(?:ward)?s?\s+\$?C\$?\b"
                   r"|rises\s+TOWARD\s+C\b",
        "exempt_if_nearby": r"\[retracted \d{4}-\d{2}-\d{2}: saturation-approach-monotone\]"
                            r"|not monotone|NOT monotone|is \\emph\{not\} monotone",
        # Documents whose ARGUMENT rests on this claim.
        "cited_by": {
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex": "reviewed 2026-09-04",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "reviewed 2026-09-04",
            # CODE-B M6: the test still ASSERTS monotonicity on the doubling
            # grid Paper 0 names as the trap.  Insert 21 into ns and it fails.
            "tests/test_paper1_block_spectrum.py": "remediated 2026-09-04 (DELTA #7); review owed (DELTA #8)",
        },
        "files": [
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "tests/test_paper1_block_spectrum.py",
        ],
    },
    {
        "id": "forced-count-260-endpoint",
        "scope": "trunk group1 group3",
        "severity": "fail",
        "retired": "2026-09-03 (trunk FULL run #3, I.0.3): the Forced-count endpoint 260 "
                   "and the matter-sector 128 came from a degenerate algebra sample "
                   "(18/24 zero elements; non-linear quark action kron(ew, m)).  With a "
                   "linear *-representation the chain ends at 32 (matter projection rank 16).",
        "pattern": r"272\s*(?:\\to|\\rightarrow|->|→|\\xrightarrow\{(?:[^{}]|\{[^{}]*\})*\})\s*260"
                   r"|full-axiom\s+(?:\$?D_F\$?\s+)?moduli\s+dimension\s+(?:at\s+\$n_\{\\max\}\s*=\s*2\$\s+)?is\s+\$?260"
                   r"|512\s*(?:\\to|→|->)\s*256\s*(?:\\to|→|->)\s*128\s*(?:\\to|→|->)\s*128"
                   r"|\\mathrm\{matter\}\}\s*(?:\\;)?=\s*(?:\\;)?128"
                   r"|128\s*per\s+generation\s+is\s+forced",
        "exempt_if_nearby": r"artefact|artifact|degenerate|corrected 2026-09-03|retired|printed before|before 2026-09-03|reproduces the retired",
                # Documents whose ARGUMENT rests on this claim (distinct
        # from `files`, which is where its wording may appear).
        "cited_by": {
            "docs/claim_test_matrix.md": "reviewed 2026-09-04",
            "papers/group3_foundations/paper_57_forced_free_seam.tex": "reviewed 2026-09-04",
        },
"files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group3_foundations/paper_57_forced_free_seam.tex",
            "papers/synthesis/*.tex",
            "docs/claim_test_matrix.md",
            "tests/test_trunk_qa_forced_count_moduli.py",
        ],
    },
    {
        "id": "all-compact-lie-groups-universality",
        "scope": "trunk group1 group3 synthesis",
        "severity": "fail",
        "retired": "2026-09-02/03 (v5.3.0 + trunk FULL run #3, I.3.1): Paper 40's class is "
                   "compact connected SEMISIMPLE; a circle/torus factor carries the circle "
                   "constant, so '4/pi universal across all compact (connected) Lie groups' "
                   "is retired.  The constant is also metric-normalisation dependent.",
        "pattern": r"(?:to|across|over|for)\s+all\s+compact\s+(?:connected\s+)?Lie\s+groups"
                   r"|rank-invariant\s+across\s+all\s+compact\s+Lie\s+groups"
                   r"|universal\s+across\s+the\s+class\s+via\s+a\s+Plancherel",
        "exempt_if_nearby": r"semisimple|torus|tori|circle factor|2/\\?pi|withdrawn|excluded|retired|corrected",
        "files": [
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/group3_foundations/paper_57_forced_free_seam.tex",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/synthesis/*.tex",
        ],
    },
    {
        "id": "rename-pass-residue-p32",
        "scope": "trunk group1",
        "severity": "fail",
        "retired": "2026-09-03 (trunk FULL run #3, I.2.1-3): variants that evaded the "
                   "patterns written for their retired twins -- 'graph Dirac operator' "
                   "(rem:D_GV_no_graph_form withdrew the graph form), 'provably disjoint' "
                   "(F1.5 retired 'provably non-overlapping'), 'Peter--Weyl propinquity' "
                   "(the paper's own object is the state-space GH truncation).",
        "pattern": r"\bgraph\s+Dirac\s+operator\b|\bDirac\s+graph\s+operator\b"
                   r"|provably\s+(?:disjoint|non-overlapping)"
                   r"|Peter--Weyl\s+propinquity",
        "exempt_if_nearby": r"withdrawn|retired|different operator|no graph form|is not a graph|Observation-level|sharing no generator|corrected",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
        ],
    },
    {
        "id": "s-p-splitting-retired-waypoints",
        "scope": "trunk group3 synthesis",
        "severity": "fail",
        "retired": "2026-09-03 (trunk FULL run #3, I.1.3/I.1.4): the binary production "
                   "lattice gives 37% at n_max=5 (the maximum), 4.1% at 20 and 0.39% at 30; "
                   "'peaking near n_max = 8', '13% at 5', '0.3% at 20', '0.005% at 30' are retired.",
        "pattern": r"peaking\s+near\s+\$?n_\{\\max\}\s*=\s*8"
                   r"|n_\{\\max\}\s*=\s*30\}?\$?:?\s*\$?\\Delta\s*E_\{\\mathrm\{rel\}\}\s*\\approx\s*0\.005"
                   r"|n_\{\\max\}\s*=\s*20\}?\$?:?\s*\$?\\Delta\s*E_\{\\mathrm\{rel\}\}\s*\\approx\s*0\.3\\%"
                   r"|n_\{\\max\}\s*=\s*5\}?\$?:?\s*\$?\\Delta\s*E_\{\\mathrm\{rel\}\}\s*\\approx\s*13",
        "exempt_if_nearby": r"did not reproduce|retired|withdrawn|corrected|measured 2026-09-03|printed before",
        "files": [
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
    },
    {
        "id": "su2-kinetic-equals-L1",
        "scope": "group5",
        "severity": "fail",
        "retired": "2026-07-04 (group5 cert run): the P30 weak-coupling kinetic form "
                   "is the co-exact completion B2 B2^T with support complementary to "
                   "L1 (L1 vanishes on the cycle space); the 'kinetic term = L1 / "
                   "proportional to L1' reading is withdrawn (P30 Prop 2 correction, "
                   "v4.65.0; the synthesis echo was caught by the cert panel).",
        "pattern": r"L_1\s*=\s*B\^?\{?\\top\}?\s*B\$?\s+as\s+the\s+kinetic\s+term"
                   r"|returns\s+Paper~?25's\s+\$?L_1.{0,30}kinetic\s+term"
                   r"|kinetic\s+form\s+reduces\s+to\s+\$?L_1"
                   r"|kinetic\s+term.{0,40}proportional\s+to\s+\$?L_1"
                   r"|\$?L_1\$?.{0,20}up\s+to\s+a\s+positive\s+scalar\s+multiple",
        "exempt_if_nearby": r"withdrawn|WITHDRAWN|wrong|incorrect|complementary"
                            r"|completes?|completion|co-exact|corrected|not\s+a\s+multiple",
        "files": [
            "papers/group5_qed_gauge/*.tex",
            "papers/synthesis/group5_qed_gauge_synthesis.tex",
            "tests/test_su2_wilson_gauge.py",
        ],
    },
    {
        "id": "wald-factor2-cone-coefficient",
        "scope": "group5",
        "severity": "fail",
        "retired": "2026-07-04 (group5 cert run): attributing the G7/G4-2 factor of 2 "
                   "to a scalar-vs-Dirac conical-defect coefficient calibration is the "
                   "2026-05-30-rejected reading (the cones share the coefficient "
                   "magnitude; the factor is Wald-forced bookkeeping). The Q2 open-"
                   "question echo was caught by the cert panel.",
        "pattern": r"scalar\s+vs\s+Dirac\s+conical-defect\s+coefficient\s+calibration"
                   r"|factor[- ]of[- ]2.{0,60}scalar[- ]vs[- ]Dirac\s+co(?:ne|efficient)",
        "exempt_if_nearby": r"incorrect|rejected|WRONG|Wald|bookkeeping|does\s+NOT|not\s+come",
        "files": [
            "papers/group5_qed_gauge/paper_51_gravity_arc.tex",
            "papers/group5_qed_gauge/paper_28_qed_s3.tex",
            "papers/synthesis/group5_qed_gauge_synthesis.tex",
        ],
    },
    {
        "id": "sbh-phi2-prefactor",
        "scope": "group5",
        "severity": "fail",
        "retired": "2026-07-04 (group5 cert run): the naive 'S_BH prefactor ~ "
                   "phi(2)/phi(1) for arbitrary cutoff' prediction was REJECTED by "
                   "G4-5d (65% deviation; the tip term lives at the log-regulated "
                   "phi(0) moment). A live echo in the Connection-to-G8 subsection "
                   "was caught by the cert panel.",
        "pattern": r"prefactor\s*\$?\\propto\s*\\phi\(2\)/\\phi\(1\)\$?\s+for\s+arbitrary\s+cutoff",
        "exempt_if_nearby": r"REJECTED|rejected|naive|\\phi\(0\)|not\s+\\phi\(2\)",
        "files": [
            "papers/group5_qed_gauge/paper_51_gravity_arc.tex",
            "papers/synthesis/group5_qed_gauge_synthesis.tex",
        ],
    },
    {
        "id": "cp2-50pct-vs-fit-floor",
        "scope": "group5",
        "severity": "fail",
        "retired": "2026-07-03 (group5 delta-2): the P25 CP^2 'no rescaling can leave "
                   "less than ~50% maximum residual against the fit' floor is INVALID as "
                   "a one-sided bound (sqrt(rmax/rmin)-1 is not a vs-fit floor; the LS "
                   "fit itself achieves 40.8%). Sharp one-sided floor = "
                   "(rmax-rmin)/(rmax+rmin) ~ 38%; the >=33% vs-data floor stands.",
        "pattern": r"50\\?\%\$?\s+maximum\s+(?:relative\s+)?residual\s+against\s+the\s+fit"
                   r"|minimax-optimal\s+rescaling\s+still\s+leaves\s+50",
        "exempt_if_nearby": r"overstat|supersed|WITHDRAWN|retract|invalid|not\s+a\s+valid"
                            r"|corrected",
        "files": [
            "papers/group5_qed_gauge/paper_25_hopf_gauge_structure.tex",
            "papers/synthesis/group5_qed_gauge_synthesis.tex",
            "tests/test_su3_wilson_s5.py",
        ],
    },
    {
        "id": "drake-swainson-3d-mistranscription",
        "scope": "group5",
        "severity": "fail",
        "retired": "2026-07-03 (group5 delta-2): the P36 3D Bethe-log reference "
                   "-0.005249 was a mistranscription of Drake--Swainson 1990 Table I "
                   "(actual -0.0052321481; residual +0.07%, not -0.24%).",
        "pattern": r"0\.005249",
        "exempt_if_nearby": r"mistranscrib|corrected|supersed|earlier\s+printing|WITHDRAWN",
        "files": [
            "papers/group5_qed_gauge/paper_36_bound_state_qed.tex",
            "tests/test_paper36_lamb_chain.py",
            "tests/paper36_lamb_support/*.py",
        ],
    },
    {
        "id": "cheeger-simons-cone-attribution",
        "scope": "group5",
        "severity": "fail",
        "retired": "2026-07-03 (group5 delta-2): 'Cheeger--Simons' is a phantom "
                   "co-author pair for the spinor conical-defect heat-kernel "
                   "coefficient (cites are Cheeger 1983 solo + Solodukhin 1995 solo; "
                   "Cheeger--Simons is the unrelated differential-characters work). "
                   "Correct label: Cheeger--Solodukhin.",
        "pattern": r"Cheeger--?Simons",
        "exempt_if_nearby": r"differential\s+character|corrected|phantom|WITHDRAWN",
        "files": [
            "papers/group5_qed_gauge/*.tex",
            "papers/synthesis/group5_qed_gauge_synthesis.tex",
        ],
    },
    {
        "id": "withdrawn-pythagorean-mechanism",
        "scope": "group1",
        "severity": "fail",
        "retired": "2026-06-18 (Paper 39): the SU(2)xSU(2) Pythagorean operator-norm "
                   "identity C_3<1->1 is FALSE on the real CH harmonics; live bound is "
                   "the triangle C_3>=1->sqrt(2).",
        # the zombie SIGNATURES (the legit product-metric 'Pythagorean d^2=d_a^2+d_b^2'
        # / 'Pythagorean triangle inequality' / 'sup-norm Pythagorean' are NOT matched)
        "pattern": r"Pythagorean\s+refinement"
                   r"|graded\s+Pythagorean\s+(?:operator-norm|Leibniz)"
                   r"|Pythagorean\s+operator-norm\s+(?:formula|identity)",
        # exempt: a withdrawal flag nearby, OR the DIFFERENT (live, legit) Paper-43
        # *modular* Pythagorean HS-orthogonality (||H-D||^2 = ||H||^2 + ||D||^2), which
        # is a genuine result, not the withdrawn tensor-C_3 refinement.
        "exempt_if_nearby": r"WITHDRAWN|withdrawn|retract|\bfalse\b|operator-norm-false"
                            r"|historical|do\s+NOT\s+use|not\s+the\s+live|triangle"
                            r"|modular|Hilbert--Schmidt|HS-orthogonal|\bHS\b|orthogonal",
        "files": [
            "geovac/gh_convergence_tensor.py",
            "geovac/gh_convergence.py",
            "papers/group1_operator_algebras/paper_39_*.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
        ],
    },
    {
        "id": "s7-structural-negative",
        "scope": "group1",
        "severity": "fail",
        "retired": "2026-06-23 (Paper 50 Erratum, S8): the S^7 scalar 'structural "
                   "non-match' was a FALSE NEGATIVE (30-dps under-resolved search); the "
                   "ladder GENERATES in-ring closed forms at every odd rung S^3..S^11.",
        "pattern": r"S\^?\{?7\}?[^.\n]{0,70}(?:structural\s+non-match|PSLQ\s+fails"
                   r"|scalar\s+negative)"
                   r"|S\^?\{?7\}?[^&\n]{0,40}UNKNOWN",
        "exempt_if_nearby": r"Erratum|false\s+negative|generates|in-ring|\\mathcal\{R\}"
                            r"|DONE|ladder|earlier\s+draft|earlier\s+version",
        "files": [
            "papers/group1_operator_algebras/paper_50_*.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
        ],
    },
    {
        "id": "batch2-false-closure",
        "scope": "group1",
        "severity": "fail",
        "retired": "2026-06-22 (Papers 42/32): 'D_W (or D_CH) does not close the period "
                   "closure' is FALSE -- conjugation by the scalar -I closes; the real "
                   "distinction is operator-level e^{i2pi D_W}=-I (double cover) vs "
                   "e^{i2pi K_alpha}=+I.",
        "pattern": r"(?:D_?\{?CH\}?|D_W)[^.\n]{0,55}(?:would\s+not|does\s+not|cannot)"
                   r"[^.\n]{0,35}(?:close|closure|produce\s+the\s+bit-exact)",
        "exempt_if_nearby": r"corrected|operator-level|double\s+cover|-I\b|\+I\b|scalar\s+-?I",
        "files": [
            "papers/group1_operator_algebras/paper_42_*.tex",
            "papers/group1_operator_algebras/paper_32_*.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
        ],
    },
    {
        "id": "propinquity-as-achieved-metric-group5",
        "scope": "group5",
        "severity": "fail",
        "retired": "group5 1st cert (2026-07-03): P51 carried two "
                   "'Latremoliere propinquity' keystone-metric restatements "
                   "(:159 abstract-adjacent, :2233 inside a Lemma) that the "
                   "line-based group1 pattern missed (converge/Latr split "
                   "across a line break). Post-fix the gated group5 scope has "
                   "ZERO legitimate Latremoliere mentions, so the bare "
                   "pattern is safe at fail severity here.",
        "pattern": r"Latr[^\n]{0,30}propinquity",
        "exempt_if_nearby": r"named\s+gap|state-space|strictly\s+stronger"
                            r"|NOT\s+the|historical|retract",
        "files": [
            "papers/group5_qed_gauge/*.tex",
            "papers/synthesis/group5_qed_gauge_synthesis.tex",
        ],
    },
    {
        "id": "propinquity-as-achieved-metric",
        "scope": "group1",
        # PROMOTED advisory -> fail 2026-08-31 (PI direction).  The trunk
        # criteria name this exact overclaim, but at advisory severity the
        # gate could only print a note about it -- soft on the one claim it
        # most specifically guards.  The old "noisy" rationale is stale:
        # exempt_if_nearby now covers the legitimate framework / named-gap /
        # descope / state-space mentions, and the corpus shows 0 live hits
        # against 2 correctly-exempted ones -- a no-op today, a guard later.
        "severity": "fail",
        "retired": "Papers 38/39/40 establish van Suijlekom STATE-SPACE GH, NOT the "
                   "strictly-stronger Latremoliere quantum-GH propinquity (a named gap). "
                   "Flag 'propinquity' asserted as the ACHIEVED convergence metric.",
        "pattern": r"(?:converge\w*|established?|proves?|in\s+the)[^.\n]{0,45}"
                   r"Latr[^.\n]{0,25}propinquity"
                   r"|propinquity\s+(?:sense|convergence)\s+at\s+(?:quantitative|explicit)",
        "exempt_if_nearby": r"not\s+claimed|named\s+gap|descoped|WITHDRAWN|degenerac"
                            r"|state-space|strictly\s+stronger|open|target|annihilat"
                            r"|historical|retract|weak-form|NOT\s+a",
        "files": [
            # 38 and 32 added 2026-08-31: the /qa trunk run found Paper 38
            # carried NO C16 entry at all, while its trunk C7 criterion names
            # this exact overclaim. 32 is the other group1 trunk root.
            "papers/group1_operator_algebras/paper_38_*.tex",
            "papers/group1_operator_algebras/paper_32_*.tex",
            "papers/group1_operator_algebras/paper_39_*.tex",
            "papers/group1_operator_algebras/paper_40_*.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            # code-docstring scope (the 2026-06-23/24 docstring-gate: C16 had scanned
            # papers only; the recurring code-docstring stale-prose class -- state-space
            # GH mislabeled "Latremoliere propinquity", retracted-convergence-as-live --
            # lived in these backing modules' docstrings. Advisory severity = the
            # fix-on-sight NIT bar set at group1 certification, v4.49.0):
            "geovac/lorentzian_propinquity_compact_temporal.py",
            "geovac/gh_convergence.py",
            "geovac/gh_convergence_tensor.py",
        ],
    },
    {
        "id": "withdrawn-c3op-envelope-sqrt",
        "scope": "group1",
        "severity": "fail",
        "retired": "2026-06-23 (Papers 45/46; surfaced live in 47+synthesis on the "
                   "first whole-group /qa): the operator-norm 'C3^op / Cthreejoint' "
                   "envelope constant sqrt(1 - 1/n_max) (= sup_{N<=2n_max-1} "
                   "sqrt((N-1)/(N+1))) is operator-norm-FALSE; the correct Paper-38 "
                   "Lemma-L3 (gradient/translation seminorm) value is C_3 = 1. NOTE: "
                   "Paper 38's own per-harmonic sqrt((N-1)/(N+1)) gradient ratio is the "
                   "LEGIT form and is a DIFFERENT expression -- not matched here.",
        # zombie signature = the ENVELOPE form sqrt(1 - 1/n_max) specifically
        # (matches \sqrt{1 - 1/\nmax}, \sqrt{1-1/n_{\max}}, sqrt(1 - 1/n_max)).
        "pattern": r"\\?sqrt\s*[\{(]\s*1\s*-\s*1\s*/\s*\\?n_?\{?\\?max",
        "exempt_if_nearby": r"withdrawn|WITHDRAWN|retract|operator-norm-false"
                            r"|\bfalse\b|earlier|historical|App\.?~?A\.3|do\s+NOT",
        "files": [
            # Paper 38 added 2026-09-01: it is the home of the LEGITIMATE
            # gradient-normalised cousin, so it was never in scope -- and
            # that scope gap is exactly the class the GATE SELF-AUDIT RULE
            # warns about (a gate silent by construction is indistinguishable
            # from a gate that passes).  P38 now prints the envelope
            # expression explicitly, to warn readers off substituting it for
            # the sub-envelope values; two-way discrimination proven at the
            # time of widening (fires on a bare occurrence, exempt on the
            # disclosed one via the "false" trigger).
            "papers/group1_operator_algebras/paper_38_*.tex",
            "papers/group1_operator_algebras/paper_44_*.tex",
            "papers/group1_operator_algebras/paper_45_*.tex",
            "papers/group1_operator_algebras/paper_46_*.tex",
            "papers/group1_operator_algebras/paper_47_*.tex",
            "papers/group1_operator_algebras/paper_48_*.tex",
            "papers/group1_operator_algebras/paper_49_*.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "geovac/lorentzian_propinquity_compact_temporal.py",
        ],
    },
    {
        "id": "tc-qubit-validated-success",
        "scope": "group2",
        "severity": "fail",
        "retired": "2026-06-27 (Papers 15/17; surfaced live on the group2 re-cert): the "
                   "transcorrelated (TC) qubit pipeline 'has been validated / succeeds / "
                   "eliminates the basis divergence (5.3->8.2 pct)' is FALSE -- the Track "
                   "BX-3 benchmark was a qubit-(Fock)-space-diagonalization false positive "
                   "(wrong-particle-number sectors below the variational bound). Under "
                   "particle-number-projected FCI (Track TC-V) the standard pipeline "
                   "converges to ~2.0 pct and TC plateaus at ~3.4 pct (WORSE). The cusp is "
                   "an energy-evaluation, not a wavefunction, problem.",
        "pattern": r"transcorrelated[^.\n]{0,90}(?:has\s+been\s+validated|is\s+validated|succeeds)"
                   r"|\bTC\b[^.\n]{0,40}(?:has\s+been\s+validated|\bvalidated\b|succeeds)"
                   r"|eliminat\w+\s+the\s+basis\s+divergence"
                   r"|from\s+divergent[^.\n]{0,45}to\s+convergent",
        "exempt_if_nearby": r"false\s+positive|dead\s+end|non-Hermitian|wrong-particle-number"
                            r"|Fock\)?\s+space|not\s+pursued|worse\s+than|particle-number-projected"
                            r"|TC-V|WITHDRAWN|retract|\bfalse\b",
        "files": [
            "papers/group2_quantum_chemistry/paper_15_*.tex",
            "papers/group2_quantum_chemistry/paper_17_*.tex",
            "papers/group2_quantum_chemistry/paper_fci_*.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
        ],
    },
    {
        "id": "pair-diagonal-as-exact-sparsity",
        "scope": "group4",
        "severity": "fail",
        "retired": "2026-06-28 (Papers 14/20; CF-1, the A/B dual-rule framing -- "
                   "criteria.md 'Dual-rule ERI framing'): the composed/atomic Pauli-sparsity "
                   "advantage is realized under the PAIR-DIAGONAL ERI approximation (rule A, "
                   "q=mc-ma, m_a=m_c & m_b=m_d), NOT the exact global-M_L Coulomb selection "
                   "rule (rule B). A sparsity claim that presents a pair-diagonal number as the "
                   "EXACT/full Gaunt-selection-rule value -- without disclosing the "
                   "pair-diagonal approximation -- is a framing zombie. Under rule B the LiH "
                   "market test is PARITY (838 vs 907) and the d-block is DENSER (30.0 vs 27.9).",
        # zombie signatures: the 2.7x-vs-STO-3G market test; the d-block-sparser / 9.23-as-
        # genuine-selection-rule claim. (The LEGIT disclosed forms carry 'pair-diagonal' /
        # 'approximation' nearby and are exempted.)
        # NOTE: broadened 2026-06-28 after the group4 first-cert FAIL surfaced
        # ~10 C16-dodging phrasings ("cheaper to encode", "more economical",
        # "structurally sparser than s/p", bare "2.7x Pauli") -- the reviewers
        # caught them; these patterns now backstop the recurrence.
        "pattern": r"2\.7\s*(?:x|×|\\times|\$\\times\$)?\s*(?:fewer|less|Pauli)"
                   r"|9\.23[^.\n]{0,80}(?:Gaunt|restrictive|sparser|selection\s+rule|economical|cheaper)"
                   r"|(?:lower|sparser)[^.\n]{0,40}9\.23"
                   r"|d-?orbital[^.\n]{0,60}(?:sparser|more\s+restrictive\s+Gaunt|cheaper|economical)"
                   r"|d-electron[^.\n]{0,40}cheaper"
                   r"|structurally\s+sparser\s+than\s+\$?s\$?/?\$?p\$?"
                   r"|more\s+economical\s+angular",
        "exempt_if_nearby": r"pair-diagonal|pair\s+diagonal|approximation|global-M_L"
                            r"|global\s+rule|global-?ML|exact\s+rule|\bparity\b|CF-1|disclos"
                            r"|left\s+on\s+the\s+table|artifact\s+of",
        "files": [
            "papers/group4_quantum_computing/paper_14_*.tex",
            "papers/group4_quantum_computing/paper_20_*.tex",
            "papers/synthesis/group4_quantum_computing_synthesis.tex",
        ],
    },
    {
        "id": "organics-in-library",
        "scope": "group4",
        "severity": "fail",
        "retired": "2026-06-28 (v4.52.0 library decision): CH2O/C2H2/C2H6 are "
                   "non-buildable and were REMOVED from the shipping library "
                   "(37 systems = 35 composed + He + H2). The 6th cert (v4.60.0) "
                   "found + removed surviving P14 tab:multi_center rows and prose "
                   "counts; this entry backstops any re-surfacing of the organics "
                   "as live library members / with live Pauli counts.",
        # zombie signature: an organic presented with a live count or as a
        # library row (a mention with a removed/dropped qualifier is exempt)
        "pattern": r"C\$?_\{?2\}?\$?H\$?_\{?[26]\}?\$?\s*(?:&|both\s+yield|at\s+\$?Q|yields?)"
                   r"|CH\$?_\{?2\}?\$?O\s*(?:&|both\s+yield|at\s+\$?Q|yields?)",
        "exempt_if_nearby": r"removed|dropped|non-buildable|not\s+(?:in|part\s+of)\s+the"
                            r"|historical|retired|de-shipped",
        "files": [
            "papers/group4_quantum_computing/*.tex",
            "papers/synthesis/group4_quantum_computing_synthesis.tex",
        ],
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
    # de-dup, preserve order
    seen, uniq = set(), []
    for p in out:
        if p not in seen and p.is_file():
            seen.add(p)
            uniq.append(p)
    return uniq


def scan_entry(entry: dict) -> "tuple[list, list]":
    """Return (live_hits, exempt_hits); each item = (relpath, line_no, snippet)."""
    pat = re.compile(entry["pattern"], re.IGNORECASE)
    # Two-tier exemption (FULL #6 fix, 2026-09-05).  The standardized per-entry
    # marker exempts only the hit's OWN line -- a marker attaches to the
    # specific withdrawn instance.  The legacy `exempt_if_nearby` vocabulary
    # keeps the +-WINDOW window (85+ loci depend on it).  BEFORE this, a marker
    # anywhere in +-WINDOW sheltered EVERY occurrence of the same entry, so a
    # live verbatim zombie a few lines from a legitimately-withdrawn instance
    # escaped -- P7 L127 "recovering the exact Coulomb degeneracy" was exempted
    # by the DIFFERENT-sub-claim marker on L123, within the window.
    marker_re = re.compile(withdrawal_marker(entry["id"]), re.IGNORECASE)
    legacy_re = re.compile(entry["exempt_if_nearby"], re.IGNORECASE)
    live, ok = [], []
    for path in _resolve(entry["files"]):
        text = path.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        hit_lines = set()
        # Match on a markup-stripped copy; report the original (2026-09-04).
        stripped = [_strip_markup(l) for l in lines]
        # Match on the markup-stripped copy OR the raw line (DELTA #7,
        # CODE-B M4): stripping alone silently disabled three declared
        # alternatives written in typeset form ($D_{\mathrm{GV}}^{\mathrm{graph}}$,
        # $N_{\mathrm{matter}} = 128$).  Either spelling now fires.
        for i, line in enumerate(stripped):
            if pat.search(line) or pat.search(lines[i]):
                hit_lines.add(i)
        # 2026-09-03 (trunk FULL run #3, I.4.1): a phrase wrapped across a
        # line break was invisible to the per-line scan -- the entry written
        # for `2 Vol(S^1) / Vol(SU(2))` could not fire on its own declared
        # locus.  Scan the newline-joined text too (offsets preserved, so a
        # match maps back to the line where it starts).
        joined = text.replace("\n", " ")
        for m in pat.finditer(joined):
            hit_lines.add(text.count("\n", 0, m.start()))
        for i in sorted(hit_lines):
            line = lines[i]
            rel = path.relative_to(ROOT)
            snip = re.sub(r"\s+", " ", line.strip())[:160]
            # (a) standardized marker: within +-MARKER_WINDOW (covers a
            #     hard-wrapped sentence whose marker lands on the adjacent
            #     physical line; does NOT shelter a distinct occurrence lines
            #     away -- the +-5 cross-shelter bug).
            mlo, mhi = max(0, i - MARKER_WINDOW), min(len(lines), i + MARKER_WINDOW + 1)
            marker_txt = "\n".join(stripped[mlo:mhi]) + "\n" + "\n".join(lines[mlo:mhi])
            same_line = bool(marker_re.search(marker_txt))
            # (b) legacy exempt_if_nearby vocabulary: +-WINDOW window.
            lo, hi = max(0, i - WINDOW), min(len(lines), i + WINDOW + 1)
            window_txt = "\n".join(stripped[lo:hi]) + "\n" + "\n".join(lines[lo:hi])
            legacy = bool(legacy_re.search(window_txt))
            if same_line or legacy:
                ok.append((rel, i + 1, snip))
            else:
                live.append((rel, i + 1, snip))
    return live, ok


def _dependency_report(gate: str | None) -> int:
    """Which documents DEPEND on each retracted claim, and were they revisited?

    See the `cited_by` note at the top of the REGISTRY.  A dependent is a
    document whose argument rests on the claim -- not merely one where the
    claim's wording might appear (that is `files`).  Dependents are listed
    once, from the argument, so this report cannot be defeated by spelling,
    which is what defeated four successive pattern rebuilds.
    """
    rows, unreviewed = [], 0
    for e in REGISTRY:
        dep = e.get("cited_by")
        if not dep:
            continue
        if gate and gate not in str(e.get("scope", "")):
            continue
        pending = [d for d, stamp in dep.items() if not stamp]
        rows.append((e["id"], dep, pending))
        unreviewed += len(pending)

    scope = f"scope '{gate}'" if gate else "ALL entries"
    print(f"claim-dependency report (C16 cited_by)   [{scope}: "
          f"{len(rows)} entr(y/ies) with declared dependents]\n")
    for eid, dep, pending in rows:
        mark = "ok " if not pending else "GAP"
        print(f"  [{mark}] {eid}")
        for d, stamp in dep.items():
            print(f"        {'REVISITED  ' if stamp else 'NOT REVISITED'}  "
                  f"{d}" + (f"   ({stamp})" if stamp else ""))
    if unreviewed:
        print(f"\nRESULT: FAIL ({unreviewed} dependent document(s) rest on a "
              f"retracted claim and have not been revisited since it was "
              f"retracted). Review each, then stamp it in `cited_by`.")
        return 1
    print(f"\nRESULT: PASS (every declared dependent of a retracted claim has "
          f"been revisited)")
    return 0


def _coverage_report(gate: str) -> int:
    """Which gated documents does each in-scope entry NOT declare?

    A `files` list that omits a gated document is a silent scope hole: the
    entry runs, reports clean, and never looked there.  The FULL #4 critic
    found four of these by hand; this makes the question a command.
    """
    _in_scope, scope_files, _w = qa_scopes.make_predicate(gate)
    gated = sorted({p.replace("\\", "/") for p in scope_files})
    print(f"C16 scope coverage   [--gate {gate}: {len(gated)} document(s)]\n")
    holes = 0
    for e in REGISTRY:
        if gate not in str(e.get("scope", "")):
            continue
        declared = set(e.get("files", []))
        missing = [g for g in gated
                   if not any(g.endswith(d) or d.endswith(g.split("/")[-1])
                              for d in declared)]
        mark = "ok " if not missing else "GAP"
        print(f"  [{mark}] {e['id']}"
              + (f"\n        not declared: "
                 + ", ".join(m.split("/")[-1] for m in missing) if missing else ""))
        holes += bool(missing)
    print(f"\nRESULT: {holes} entry/entries do not declare every gated "
          f"document. That is advisory -- a claim that cannot surface in a "
          f"document need not name it -- but each GAP is a place the entry "
          f"reports clean without looking.")
    return 0


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    show_all = "--all" in sys.argv
    gate = _gate_substr(sys.argv)
    if "--coverage" in sys.argv:
        return _coverage_report(gate or "trunk")
    if "--dependencies" in sys.argv:
        return _dependency_report(gate)
    scope = f"scope '{gate}'" if gate else "ALL entries"

    # Entry selection is LOCUS-DERIVED as well as tag-based (2026-08-31
    # gate-scope audit).  The hand-maintained `scope` tag had drifted from
    # the `files` loci it is supposed to summarise, and the drift was
    # invisible: NO entry carried the tag "trunk", "synthesis", "paper_58",
    # "paper_59" or "paper_60", so `--gate trunk` selected 0 of 27 entries
    # and printed PASS having checked nothing -- the C19 failure shape, in
    # the gate the protocol leans on hardest.  The sharpest instance: the
    # `propinquity-as-achieved-metric` entry lists paper_38 and paper_32 --
    # both TRUNK papers, added by the 2026-08-31 trunk run itself -- under
    # the tag "group1", so the trunk gate could never see its own fix.
    #
    # An entry now runs whenever ANY of its declared loci lies in the gated
    # scope.  The tag is kept as a widening fallback so a locus pattern that
    # matches nothing on disk cannot silently narrow an entry out.
    _in_scope, _scope_files, _scope_warnings = (
        qa_scopes.make_predicate(gate) if gate else (None, [], []))
    qa_scopes.emit_warnings(_scope_warnings)

    def _locus_gated(pattern: str) -> bool:
        if _in_scope is None:
            return True
        hits = glob.glob(str(ROOT / pattern))
        return any(_in_scope(h) for h in hits)

    def selected(e: dict) -> bool:
        if gate is None or e["scope"] == "all" or gate in e["scope"]:
            return True
        return any(_locus_gated(f) for f in e.get("files", []))

    # Coverage accounting (2026-08-31 gate-scope audit).  This gate selects
    # REGISTRY ENTRIES by their own `scope` field, not files -- so a --gate
    # value matching no entry would run zero checks and still print PASS,
    # the same shape as the C19 bug.  Count what actually ran and say so.
    _selected = [e for e in REGISTRY if selected(e)]
    _files_seen: set = set()
    if not _selected:
        print(f"   [scope] WARNING: --gate '{gate}' selected 0 of "
              f"{len(REGISTRY)} registry entries -- this run checks "
              f"NOTHING. Add a '{gate}' entry scope, or use a scope "
              f"name that exists.")

    # GATE SELF-AUDIT (ported from C17, 2026-09-07).  C17 grew this at FULL
    # run #4; C16 -- the gate the protocol leans on hardest -- never got it.
    # Measured on the 58/59/60 FULL run:  `--gate paper_59` and
    # `--gate paper_60` selected ZERO entries, printed the WARNING above, and
    # still returned 0, so the run reported "C16 PASS" for two papers the
    # gate had not looked at.  Entry COUNT is not the honest measure; what
    # matters is how many selected entries declare a locus INSIDE the scope.
    def _has_gated_locus(e: dict) -> bool:
        if _in_scope is None:
            return True
        return any(_locus_gated(f) for f in e.get("files", []))

    _grounded = [e for e in _selected if _has_gated_locus(e)]
    if gate is not None and not _grounded:
        print(f"   [scope] ERROR: --gate '{gate}' selected "
              f"{len(_selected)} entry/entries, but NONE of them declares a "
              f"locus inside this scope. The gate would examine nothing and "
              f"print PASS. Add an entry for this target (C16 maintenance "
              f"rule) rather than trusting this run.")
        return 1

    # An entry is marker-only if its exemption is the never-match "(?!)" or
    # begins with the (escaped) standardized marker.  The previous test
    # compared against a bare "[retracted", which no escaped string ever
    # starts with, so this counter could never decrease (DELTA #7, N4).
    _authored = [e["id"] for e in REGISTRY
                 if not (e.get("exempt_if_nearby", "").strip() in ("(?!)", r"(?!)")
                         or e.get("exempt_if_nearby", "").strip()
                            .startswith(("\\[retracted", "[retracted")))]
    fail_hits, advisory_hits, exempt_total = [], [], 0
    print(f"retracted-claims / zombie-drift screen   [{scope}: "
          f"{len(_selected)}/{len(REGISTRY)} entries]\n")
    for e in REGISTRY:
        if not selected(e):
            continue
        for _pat in e.get("files", []):
            _files_seen.add(_pat)
        live, ok = scan_entry(e)
        exempt_total += len(ok)
        tag = "FAIL" if e["severity"] == "fail" else "ADVISORY"
        status = "clean" if not live else f"{len(live)} LIVE"
        print(f"  [{tag}] {e['id']}: {status}  (exempt/withdrawn-flagged: {len(ok)})")
        for rel, ln, snip in live:
            (fail_hits if e["severity"] == "fail" else advisory_hits).append(
                (e["id"], rel, ln, snip))
        if show_all:
            for rel, ln, snip in ok:
                print(f"        [exempt] {rel}:{ln}  {snip}")

    if fail_hits:
        print(f"\n*** LIVE RETRACTED CLAIM(S) ({len(fail_hits)}) -- a withdrawn claim "
              f"re-surfaced WITHOUT a withdrawal flag in {scope}: ***")
        for cid, rel, ln, snip in fail_hits:
            print(f"  [{cid}] {rel}:{ln}\n      {snip}")

    if advisory_hits:
        print(f"\n--- ADVISORY ({len(advisory_hits)}) -- review (does NOT fail the gate): ---")
        for cid, rel, ln, snip in advisory_hits:
            print(f"  [{cid}] {rel}:{ln}\n      {snip}")

    if fail_hits:
        print(f"\nRESULT: FAIL ({len(fail_hits)} live retracted "
              f"claim(s) in {scope}; {len(_selected)}/{len(REGISTRY)} "
              f"entries over {len(_files_seen)} declared locus "
              f"pattern(s))")
        return 1
    # cited_by enforcement (2026-09-04).  Two failure modes, both of which
    # previously printed PASS:  dependents declared but never revisited, and
    # entries declaring nothing at all -- the latter indistinguishable from
    # "no dependents", which is how C17 came to examine nothing.
    # The ratchet must use the SAME selector as the scan.  It previously read
    # the `scope` TAG only, while entry selection above is locus-derived -- so
    # on any target whose name is not literally a scope tag (every single-paper
    # target: paper_58/59/60) `_pending` and `_undeclared` were both empty and
    # the two FAILs added 2026-09-04 could not fire.  Same class as the C11
    # path-convention defect: a check proven to fire in one selection path,
    # silently inert in another.  Found by the 58/59/60 FULL run.
    _in = selected
    _pending = [(e["id"], d) for e in REGISTRY if _in(e)
                for d, v in e.get("cited_by", {}).items() if not v]
    _undeclared = [e["id"] for e in REGISTRY if _in(e) and "cited_by" not in e
                   and e["id"] not in CITED_BY_BASELINE]
    if _pending:
        print(f"\n*** UNREVISITED DEPENDENTS ({len(_pending)}) -- documents whose "
              f"argument rests on a retracted claim, not yet revisited: ***")
        for eid, d in _pending:
            print(f"  [{eid}] {d}")
        print("  A pattern cannot catch this class:\u00a0a citing document restates "
              "the claim in its own words.  Review each, then stamp it in "
              "`cited_by`.")
    if _undeclared:
        print(f"\n*** ENTRIES NOT DECLARING DEPENDENTS ({len(_undeclared)}) -- "
              f"silence here is indistinguishable from 'nothing depends on "
              f"this'. Declare them, or `cited_by: {{}}` if genuinely none: ***")
        for eid in _undeclared:
            print(f"  [{eid}]")
    print(f"   [marker] {len(_authored)} of {len(REGISTRY)} entries still rely "
          f"on hand-authored exemption vocabulary rather than the standardized "
          f"`[retracted YYYY-MM-DD]` token -- the class that produced three "
          f"false-clean entries on 2026-09-03.  New entries use the token only.")
    if _pending or _undeclared:
        print(f"\nRESULT: FAIL ({len(_pending)} unrevisited dependent(s), "
              f"{len(_undeclared)} entry/entries not declaring dependents, in "
              f"{scope}). No live retracted wording -- but the dependency "
              f"class is what patterns cannot see.")
        return 1
    print(f"\nRESULT: PASS (no live fail-severity retracted claim in "
          f"{scope}; {len(_selected)}/{len(REGISTRY)} entries over "
          f"{len(_files_seen)} declared locus pattern(s); {exempt_total} "
          f"occurrence(s) correctly carry a withdrawal flag)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
