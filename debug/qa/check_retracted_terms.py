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
        "id": "p11-h2plus-0002pct-prose",
        # The CLAIM in words, for when the numeral is gone but the magnitude
        # survives as prose.  Deliberately narrow: only phrasings that assert a
        # ~1e-6-relative accuracy for H2+, which is the retired magnitude.
        "pattern": r"(two parts in a million|2 parts in 10\^6|"
                   r"two[- ]in[- ]a[- ]million)",
        "require_nearby": r"H\$?_2\^?\{?\+|H2\+|prolate|spectral|Laguerre",
        "exempt_if_nearby": withdrawal_marker("p11-h2plus-0002pct-prose"),
        "severity": "fail",
        "scope": "paper_11 paper_12 paper_13 paper_15 paper_17 group2 synthesis trunk",
        "note": "Registered 2026-09-19 (/qa group2 CODE run).  The NUMERIC half "
                "is C17 family p11-h2plus-0002pct-retired (32 live loci).  "
                "MEASURED: the spectral solver reproduces E_ref to 3.6e-14 Ha "
                "at n_basis=20, R=2.0 -- i.e. to the precision at which E_ref "
                "is conventionally quoted -- so 0.0002% (1.21e-6 Ha) understates "
                "the method by ~7.6 orders.  Even n_basis=5 is 57x better than "
                "the published claim.  Canonical form (PI direction): state it "
                "QUALITATIVELY as machine precision; do NOT substitute another "
                "percentage, because the mantissa is reference-limited.  "
                "Registry keys: p11_h2plus_err_ha, p11_h2plus_req_bohr.",
        "files": [
            "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/synthesis/geovac_field_guide.tex",
            "README.md",
            "CLAUDE.md",
        ],
        # Documents whose ARGUMENT rests on H2+'s accuracy figure -- each cites
        # it to motivate its own geometry/level choice.  Unstamped dependents
        # fail the gate; stamping is a REVIEW outcome, not an edit.
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "reviewed 2026-09-19 -- abstract L34 + intro L120 + hierarchy L1349 -> machine precision",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-19 -- L101 prose + tab:hierarchy Level-2 cell + caption (marker added)",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-19 -- L102 prose + L1418 table; FD figure 0.70% -> 1.01% per P11's own table",
            "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex": "reviewed 2026-09-19 -- L72 prose + L1542 table; FD figure 0.70% -> 1.01%",
            "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex": "reviewed 2026-09-19 -- L1349 -> machine precision; the R_eq 2.001/0.21% FD claim left intact",
            "papers/group2_quantum_chemistry/paper_fci_molecules.tex": "reviewed 2026-09-19 -- L766 -> machine precision; guardrail negative untouched",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-19 -- L87, L169, L279 (the -0.6026 vs -0.6026 form), L1068",
            "papers/synthesis/geovac_field_guide.tex": "reviewed 2026-09-19 -- L224 hierarchy table cell",
        },
    },
    {
        "id": "graph-eigenvalues-integer-or-rational",
        # The CLAIM: the S^3 graph's OWN eigenvalues are integers/rationals, or the
        # graph spectrum converges to the integer ladder n^2-1.  Both false: the
        # graph (block diagonal in l, each block a product of two paths) has the
        # closed-form spectrum {2-2cos(j pi/(n_max-l)) + 2-2cos(k pi/(2l+1))}
        # (Paper 0 Sec. VI) -- algebraic, generally irrational (e.g. 2-sqrt 2),
        # bounded in [0,8].  The integer ladder -(n^2-1) is the CONTINUUM S^3
        # Laplace-Beltrami spectrum.  pi-free survives: integer matrix entries
        # give algebraic eigenvalues, never transcendental.
        "pattern": r"(eigenvalues,\s+degeneracies,\s+and\s+coupling\s+coefficients\s+are\s+integers\s+or\s+rationals"
                   r"|graph\s+and\s+its\s+integer\s+eigenvalue\s+spectrum"
                   r"|graph:\s+rational\s+eigenvalues"
                   r"|graph\s+Laplacian\s+spectrum\s+to\s+the\s+pure\s+integers"
                   r"|that\s+the\s+graph\s+Laplacian\s+converges\s+to:\s+integer)",
        "exempt_if_nearby": withdrawal_marker("graph-eigenvalues-integer-or-rational"),
        "severity": "fail",
        "scope": "trunk group3 synthesis",
        "note": "Registered 2026-09-25 (PI-directed correction, marriage/MC session). "
                "Paper 0 Sec. VI [INTERNAL THEOREM] gives the graph spectrum in closed "
                "form; group3 synthesis L1419 already stated it correctly.  Corrected "
                "at: CLAUDE.md S1.5 (PI-authorized edit of a PI-owned section), Paper 18 "
                "obs:pi-free + Sec. origin-P24 sentence, Paper 22 L95 (was: graph "
                "spectrum CONVERGES to the pure integers -- it cannot, it is bounded) "
                "and L571 bullet, field guide L367, group3 synthesis L887.  Canonical "
                "form: matrix entries integer/rational => eigenvalues algebraic.",
        "files": [
            "CLAUDE.md",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/group3_foundations/paper_22_angular_sparsity.tex",
            "papers/synthesis/geovac_field_guide.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
        # Documents reviewed 2026-09-25 whose nearby wording was checked; each
        # speaks of the CONTINUUM integer spectrum or of rational invariants
        # and was found correct as written.
        "cited_by": {
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "reviewed 2026-09-25 -- L674 '(ii) Integer spectrum' is the continuum LB ladder; correct as written",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex": "reviewed 2026-09-25 -- L181 bare continuum LB n(n+2); correct as written",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "reviewed 2026-09-25 -- L3029 graph-derived integer invariants (Z, n); correct as written",
            "papers/synthesis/group5_qed_gauge_synthesis.tex": "reviewed 2026-09-25 -- L298 'pi-free rationals or algebraic numbers'; correct as written",
            "papers/synthesis/group6_precision_observations_synthesis.tex": "reviewed 2026-09-25 -- L97 Bargmann-Segal S^5 rational entries; L332 V_ee irrational; correct as written",
        },
    },
    {
        "id": "graph-laplacian-converges-to-s3",
        # The CLAIM: the GeoVac edge Laplacian D-A converges to / is equivalent to the S^3
        # Laplace-Beltrami operator, or has spectrum n^2-1.  REFUTED 2026-09-25: no edge changes l,
        # so D-A splits into n_max disconnected l-blocks (n_max zero modes, zero cross-block heat
        # kernel).  What survives exactly: the (n,l,m) label set + the diagonal SO(4) Casimir.
        "pattern": r"graph Laplacian converges to (the )?(continuous|continuum)?\s*\$?S\^3|converges to the (unit )?\$?S\^3\$? Laplace|whose continuum limit is (the unit three-sphere|\$?-?\\Delta_\{S\^3\})|continuum limit of the (positive-semidefinite )?graph Laplacian|discrete graph Laplacian (for hydrogen|on \$?S\^3\$?) is mathematically equivalent|graph (that )?converges to the unit three-sphere|graph discretizes \$?S\^3|has spectrum \$?n\^2 ?- ?1\$? exactly|identified with the continuum \$?S\^3\$?\s*Laplace|carries the Coulomb spectrum",
        "exempt_if_nearby": withdrawal_marker("graph-laplacian-converges-to-s3"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 synthesis",
        "note": "Registered 2026-09-25 (PI-directed re-pricing; independently re-derived twice + adversarial pass).  Evidence: tests/test_graph_s3_nonconvergence.py (fire-tested 4/4).  The S^3 equivalence is Fock's CONTINUUM projection; Paper 38/WH1 (spectral truncations) untouched; Paper 24 Coulomb/HO asymmetry Layer 1 withdrawn as a consequence (7 catalogued, 6 standing).",
        "files": [
            "CLAUDE.md",
            "docs/claims_register.md",
            "papers/group1_operator_algebras/paper_29_ramanujan_hopf.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex",
            "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex",
            "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex",
            "papers/group2_quantum_chemistry/paper_fci_atoms.tex",
            "papers/group2_quantum_chemistry/paper_fci_molecules.tex",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/paper_22_angular_sparsity.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex",
            "papers/group3_foundations/paper_57_forced_free_seam.tex",
            "papers/group4_quantum_computing/paper_23_nuclear_shell.tex",
            "papers/group5_qed_gauge/paper_25_hopf_gauge_structure.tex",
            "papers/group5_qed_gauge/paper_28_qed_s3.tex",
            "papers/group5_qed_gauge/paper_2_alpha.tex",
            "papers/group5_qed_gauge/paper_30_su2_wilson.tex",
            "papers/group5_qed_gauge/paper_41_rule_b_wilson_u1.tex",
            "papers/group5_qed_gauge/paper_51_gravity_arc.tex",
            "papers/group6_precision_observations/paper_26_entanglement.tex",
            "papers/group6_precision_observations/paper_27_entropy_projection.tex",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "papers/synthesis/geovac_field_guide.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/synthesis/group4_quantum_computing_synthesis.tex",
            "papers/synthesis/group5_qed_gauge_synthesis.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_29_ramanujan_hopf.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-25 (S3 sweep) -- L96 already scoped; no convergence claim",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group2_quantum_chemistry/paper_fci_atoms.tex": "remediated + reviewed 2026-09-25 (residual sweep) -- L90-93 matched-by-construction wording; L182 'matched scale'",
            "papers/group2_quantum_chemistry/paper_fci_molecules.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group3_foundations/paper_1_spectrum.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group3_foundations/paper_22_angular_sparsity.tex": "remediated + reviewed 2026-09-25 (PM) -- L102, L675, L709: graph Laplacian no longer called Coulomb-specific",
            "papers/group3_foundations/paper_24_bargmann_segal.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group3_foundations/paper_57_forced_free_seam.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group4_quantum_computing/paper_23_nuclear_shell.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group5_qed_gauge/paper_25_hopf_gauge_structure.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group5_qed_gauge/paper_28_qed_s3.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group5_qed_gauge/paper_2_alpha.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group5_qed_gauge/paper_30_su2_wilson.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group5_qed_gauge/paper_41_rule_b_wilson_u1.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/group5_qed_gauge/paper_51_gravity_arc.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/synthesis/geovac_field_guide.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/synthesis/group4_quantum_computing_synthesis.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
            "papers/synthesis/group5_qed_gauge_synthesis.tex": "remediated + reviewed 2026-09-25 (S3 re-pricing sweep)",
        },
    },
    {
        "id": "he-graph-native-hop-headline",
        # The CLAIM: He graph-native CI 0.216% (and 0.19%), H^- over-binding by 21%, a Z_c~1.84
        # graph-validity boundary, and graph off-diagonals carrying 10-39% of the correlation.  All
        # are energies/features of a HOP-AUGMENTED model: casimir_ci._build_graph_h1 adds a
        # non-physical (1/16)A one-body hop.  Physical fixed-k CI: 1.881% at n_max=7.
        "pattern": r"graph-native[^.]{0,80}(0\.216|0\.19)\s?\\?%|Z_c\s*\\?(approx|\u2248|~)\s*1\.84|over-binds?\s+by\s+21|graph[- ]validity[- ]boundary|10\s*(--|\u2013|-)\s*39\s*\\?%",
        "require_nearby": r"\bHe\b|helium|H\$?\^-|H\u207b|\bCI\b|correlation",
        "exempt_if_nearby": withdrawal_marker("he-graph-native-hop-headline"),
        "severity": "fail",
        "scope": "trunk group2 group3 group6 synthesis",
        "note": "Registered 2026-09-25 (PI-approved retraction).  In the k=Z hydrogenic basis T+V_ne is exactly diagonal, so the (1/16)A hop is not part of the atom's Hamiltonian (-48 mHa, more than He's correlation energy).  Registry keys he_gn_phys_n{2..7}_pct, he_gn_hop_shift_mha, hminus_phys_pct.  Guards tests/test_he_physical_hamiltonian.py + tests/test_he_hop_guard.py.  CLAUDE.md is deliberately NOT scanned: its S3 CUSP-2 row is an append-only historical record (corrected by the appended RE-INTERPRETATION row) and must not be modified, even by a token; the PM is the single writer of the corrected CLAUDE.md He rows.",
        "files": [
            "docs/claims_register.md",
            "docs/validation_benchmarks.md",
            "papers/INDEX.md",
            "papers/group1_operator_algebras/paper_29_ramanujan_hopf.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex",
            "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex",
            "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex",
            "papers/group2_quantum_chemistry/paper_fci_atoms.tex",
            "papers/group2_quantum_chemistry/paper_fci_molecules.tex",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/paper_22_angular_sparsity.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex",
            "papers/group3_foundations/paper_57_forced_free_seam.tex",
            "papers/group4_quantum_computing/paper_23_nuclear_shell.tex",
            "papers/group5_qed_gauge/paper_25_hopf_gauge_structure.tex",
            "papers/group5_qed_gauge/paper_28_qed_s3.tex",
            "papers/group5_qed_gauge/paper_2_alpha.tex",
            "papers/group5_qed_gauge/paper_30_su2_wilson.tex",
            "papers/group5_qed_gauge/paper_41_rule_b_wilson_u1.tex",
            "papers/group5_qed_gauge/paper_51_gravity_arc.tex",
            "papers/group6_precision_observations/paper_26_entanglement.tex",
            "papers/group6_precision_observations/paper_27_entropy_projection.tex",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "papers/synthesis/geovac_field_guide.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/synthesis/group4_quantum_computing_synthesis.tex",
            "papers/synthesis/group5_qed_gauge_synthesis.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "remediated + reviewed 2026-09-25 (He regression)",
            "papers/group2_quantum_chemistry/paper_fci_atoms.tex": "remediated + reviewed 2026-09-25 (He regression)",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "remediated + reviewed 2026-09-25 (residual sweep)",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated + reviewed 2026-09-25 (residual sweep)",
            "papers/group6_precision_observations/paper_26_entanglement.tex": "remediated + reviewed 2026-09-25 (He regression)",
            "papers/group6_precision_observations/paper_27_entropy_projection.tex": "remediated + reviewed 2026-09-25 (He regression)",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + reviewed 2026-09-25 (residual sweep) -- incl. He spectrum re-measured on the physical H",
            "papers/synthesis/geovac_field_guide.tex": "remediated + reviewed 2026-09-25 (residual sweep)",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "remediated + reviewed 2026-09-25 (residual sweep)",
            "papers/synthesis/group6_precision_observations_synthesis.tex": "remediated + reviewed 2026-09-25 (residual sweep)",
        },
    },
    {
        "id": "coulomb-ho-layer1-and-three-discriminators",
        # The CLAIM: the Coulomb/HO asymmetry has seven STANDING layers incl. a "spectrum-computing
        # role of L_0" on S^3, and Paper 31 has three independent discriminators.  Withdrawn
        # 2026-09-25/26 with the graph->S^3 re-pricing: L_0 is spectrally inert on both spheres.
        # The Wilson layer survives, re-scoped to representation theory (claims review, verdict B).
        "pattern": r"modular-?\s*Hamiltonian\s+structure\s+of\s+the\s+wedge\s+KMS\s+state|four-layer\s+Coulomb/HO\s+asymmetry|three independent discriminators|matter-propagator role for the node Laplacian|spectrum-computing role of \$?L_0\$?(?![^.]{0,80}withdrawn)|seven asymmetry layers|current count is seven layers|seven layers in the current count|Yes \(eigenvalues = spectrum\)|nontrivial on the gauge side|along the same\s+operator-order divide|U\(1\) only \(Sprint~?5 negative\)",
        "exempt_if_nearby": withdrawal_marker("coulomb-ho-layer1-and-three-discriminators"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 synthesis",
        "note": "Registered 2026-09-26 (PI-approved).  Seven catalogued, six standing; Paper 31: two standing discriminators.  Consequence of graph-laplacian-converges-to-s3.",
        "files": [
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/group5_qed_gauge/paper_25_hopf_gauge_structure.tex",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "CLAUDE.md",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_24_bargmann_segal.tex": "remediated + reviewed 2026-09-26 (layer-count agent) -- L454 count; L1199-1207 Wilson layer re-scoped to representation theory",
            "papers/group5_qed_gauge/paper_25_hopf_gauge_structure.tex": "remediated + reviewed 2026-09-26 (layer-count agent) -- ~L457, ~L921-926, ~L1009",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex": "remediated + reviewed 2026-09-26 (layer-count agent) -- abstract, intro, five-layers, table, checklist, conclusion; MATERIAL zombie L728-730 fixed",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + reviewed 2026-09-26 (layer-count agent) -- L1803, L1924, L3719, L3943, L4074, L4946, tab:coulomb_ho (live zombie row fixed)",
            "papers/synthesis/group3_foundations_synthesis.tex": "reviewed 2026-09-26 (PM) -- L182 already reads 'seven layers catalogued, six standing'",
            "CLAUDE.md": "reviewed 2026-09-26 (PM) -- S2 L273 already reads '7 layers catalogued, 6 standing'",
        },
    },
    {
        "id": "p11-pi-orbital-selection-rule",
        # The CLAIM (Paper 11 Limitations, H2 paragraph): <sigma_g sigma_g|pi_+ pi_-> = 0 by azimuthal
        # symmetry, so pi orbitals cannot improve the 1Sigma_g+ ground state.  FALSE: the state
        # constrains only total M; restoring azimuthal channels is worth +11.6 mHa (Paper 12, v5.11.18).
        "pattern": r"orbital selection rule|\\pi\$? orbitals cannot improve|pi orbitals cannot improve",
        "exempt_if_nearby": withdrawal_marker("p11-pi-orbital-selection-rule"),
        "severity": "fail",
        "scope": "group2 synthesis",
        "note": "Withdrawn 2026-09-26 (PI-approved Paper 11 cleanup).  Found live in Paper 11 L1188 two weeks after Paper 12 disproved it (sigma-only restriction, v5.11.18).",
        "files": [
            "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex",
        ],
        "cited_by": {},
    },
    {
        "id": "kappa-maps-graph-spectrum",
        # The CLAIM: kappa = -1/16 maps / converts the graph's eigenvalues onto the Rydberg
        # (physical, hydrogen) spectrum.  FALSE (trunk DELTA 2026-09-27): the graph spectrum lies in
        # [0,8] (n_max disconnected l-blocks); kappa matches only lambda_max -> 8 to -1/2 (an
        # Observation).  E_n comes from n^2-1 on the labels through Fock's energy shell p0^2 = -2E.
        "pattern": r"operational\s+role\s+is\s+to\s+map|maps\s+(the\s+)?(dimensionless\s+)?(GeoVac\s+)?graph\s+eigenvalues|Rydberg-to-graph-eigenvalue\s+conversion|converts\s+a\s+discrete\s+spectrum|mediating\s+between\s+graph\s+eigenvalues\s+and\s+physical|kinetic\s+scale[^.]{0,60}fixes\s+the\s+energy\s+unit|hydrogen\s+atom,\s+to\s+within\s+an\s+overall\s+scale",
        "exempt_if_nearby": withdrawal_marker("kappa-maps-graph-spectrum"),
        "severity": "fail",
        "scope": "trunk group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-28 (PI-approved remediation of the 2026-09-27 /qa trunk DELTA). Sibling of graph-laplacian-converges-to-s3; that entry's 2026-09-25 cited_by stamps were premature for these paraphrases.",
        "files": [
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/paper_22_angular_sparsity.tex",
            "papers/synthesis/geovac_field_guide.tex",
        ],
        "cited_by": {
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L294/L279/L176 kappa role",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L440 kappa signature",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L147, L3063, L3326",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L669 band-edge only",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L113, L117, L741",
            "papers/group3_foundations/paper_1_spectrum.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L54-59 kappa plays no role",
            "papers/group3_foundations/paper_22_angular_sparsity.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L660",
            "papers/synthesis/geovac_field_guide.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L158 energy-shell map",
        },
    },
    {
        "id": "graph-hamiltonian-carries-spectrum",
        # The CLAIM: the graph (edge Laplacian, graph Hamiltonian kappa*L + W, graph topology) carries,
        # encodes, discretizes or converges to the hydrogen / S^3 / HO spectrum, or has integer
        # eigenvalues.  FALSE: the v4.76.0 audit (2026-09-28) confirmed the model is unchanged since v4 and
        # no construction of it yields -1/(2n^2) (AtomicSolver n_max=20: -0.4936, -0.4936, -0.4932, ...).
        # TRUE and kept: the label set + diagonal SO(4) Casimir; Paper 1's SU(1,1) algebra on interior
        # labels; lambda_max -> 8; algebraic (pi-free) graph eigenvalues; Paper 38 (CH Dirac, not D-A).
        "pattern": r"emerge\s+from\s+the\s+full\s+graph\s+Hamiltonian|graph\s+topology\s+correctly\s+encodes|(topology|connections)[^.]{0,40}generates\s+physical\s+corrections|lattice\s+that\s+discretizes\s+the\s+Laplace|Laplacian\s+converges\s+to\s+the\s+HO\s+Hamiltonian|whose\s+spectrum\s+converges,\s+under\s+the\s+Fock\s+projection|is\s+the\s+matter\s+propagator\s+\(it\s+converges|admits\s+a\s+\S{0,8}-free\s+graph\s+Laplacian\s+interpretation|integer\s+eigenvalues,\s+rational",
        "exempt_if_nearby": withdrawal_marker("graph-hamiltonian-carries-spectrum"),
        "severity": "fail",
        "scope": "trunk group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-28 (PI-approved remediation of the 2026-09-27 /qa trunk DELTA; v4.76.0 audit: refutation stands).",
        "files": [
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/group4_quantum_computing/paper_23_nuclear_shell.tex",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex",
            "papers/synthesis/group5_qed_gauge_synthesis.tex",
            "papers/group3_foundations/paper_22_angular_sparsity.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L174 + Class-S definition",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L743 deleted; L695 T+- vs Casimir",
            "papers/group3_foundations/paper_1_spectrum.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L31, L379, L388, L464",
            "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L48, L113, L160, L1128, L1247",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L121, L1737",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L135, L208",
            "papers/group4_quantum_computing/paper_23_nuclear_shell.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L303, L1216 -> Paper 24 (L0 inert)",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L146",
            "papers/synthesis/group5_qed_gauge_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L97, L137",
            "papers/group3_foundations/paper_22_angular_sparsity.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L670 Fock rigidity",
            "papers/synthesis/group6_precision_observations_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L39",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L38",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex": "remediated + summary surfaces reread 2026-09-28 (trunk DELTA remediation agent; PM-verified gate 0 live) -- L108",
        },
    },
    {
        "id": "graph-carries-physics-residue",
        # Round-2 residue (trunk re-check 2026-09-28) of the graph/kappa-carries-the-spectrum
        # family: paraphrases and table cells no round-1 pattern anticipated -- the retired
        # "<0.1%" Level-1 accuracy cell, "graph-Laplacian limit", "Laplacian encodes the kinematic
        # manifold", "D-A is SO(4)-covariant" (false: D-A preserves l), "bare combinatorial graph
        # is right", the pre-remediation natural-lattice / discretization principle.
        "pattern": r"graph-Laplacian\s+limit|encodes\s+the\s+\S{0,12}kinematic\S{0,2}\s+(structure|manifold)|D\s*[-−]\s*A\$?\s+is\s+covariant|bare\s+combinatorial\s+graph\s+(being|is)\s+right|bare\s+combinatorial\s+graph\s+plus\s+three|natural\s+lattice-{1,3}\s*a\s+discrete\s+graph|discretizing\s+it\s+as\s+a\s+sparse\s+graph|by\s+discretizing\s+the\s+Laplace|graph\s+Laplacian\s+FCI\s+approach|that\s+defines\s+GeoVac's\s+graph\s+Laplacian|quantum\s+mechanics\s+emerges\s+from\s+discrete\s+graph|\(Fock(\s+projection)?\)\s*&\s*\$?<\s?0\.1|Accuracy:\s*\$?<\s?0\.1",
        "exempt_if_nearby": withdrawal_marker("graph-carries-physics-residue"),
        "severity": "fail",
        "scope": "trunk group2 group3 group5 synthesis",
        "note": "Registered 2026-09-28 (PI-approved remediation round 2 after the /qa trunk re-check). Sibling of kappa-maps-graph-spectrum / graph-hamiltonian-carries-spectrum.",
        "files": [
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex",
            "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex",
            "papers/group5_qed_gauge/paper_36_bound_state_qed.tex",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/group2_quantum_chemistry/paper_fci_atoms.tex",
            "papers/group2_quantum_chemistry/paper_fci_molecules.tex",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_1_spectrum.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live) -- guardrail theorem rests on the Sturmian identity + D^(n) blocks, not on the D-A sentence (rewritten: L,A generate SO(4) per shell; D-A not SO(4)-covariant)",
            "papers/group5_qed_gauge/paper_36_bound_state_qed.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/group2_quantum_chemistry/paper_fci_atoms.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/group2_quantum_chemistry/paper_fci_molecules.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (round-2 agent; PM-verified gate 0 live)",
        },
    },
    {
        "id": "graph-carries-physics-residue-r3",
        # Round-3 residue (narrow re-check #2, 2026-09-28): P13 "LB discretizes cleanly on a sparse
        # graph"; P18 "kappa absorbs the accumulated conformal rescaling"; group2 synthesis "discrete
        # graph and the continuous LB operator are dual ... proven conformal maps"; and the Level-2
        # knock-on in P0 -- "(l,m) skeleton with 2l+1 at every level" is FALSE at Level 2 (only m).
        "pattern": r"discretizes\s+cleanly\s+on\s+a\s+sparse\s+graph|absorbs\s+the\s+accumulated\s+conformal\s+rescaling|discrete\s+graph\s+and\s+the\s+continuous\s+Laplace-{1,2}Beltrami\s+operator\s+are\s+dual|same\s+discrete\s+quantum\s+number\s+skeleton|recurs\s+as\s+a\s+universal\s+structural\s+factor\s+inside\s+every\s+natural\s+geometry|lattice\s+is\s+the\s+invariant;\s+the\s+continuum\s+geometries\s+are\s+projections",
        "exempt_if_nearby": withdrawal_marker("graph-carries-physics-residue-r3"),
        "severity": "fail",
        "scope": "trunk group2 group3 group5 synthesis",
        "note": "Registered 2026-09-28 (PI-approved round 3: fix + whole-paper fix-mode reads).",
        "files": [
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group5_qed_gauge/paper_2_alpha.tex",
            "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "remediated 2026-09-28, round 3: WHOLE-PAPER fix-mode read (not locus-only); PM-verified gate 0 live",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated 2026-09-28, round 3: WHOLE-PAPER fix-mode read (not locus-only); PM-verified gate 0 live",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "remediated 2026-09-28, round 3: WHOLE-PAPER fix-mode read (not locus-only); PM-verified gate 0 live",
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex": "remediated 2026-09-28, round 3: WHOLE-PAPER fix-mode read (not locus-only); PM-verified gate 0 live",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "remediated 2026-09-28, round 3: WHOLE-PAPER fix-mode read (not locus-only); PM-verified gate 0 live",
            "papers/group3_foundations/paper_1_spectrum.tex": "remediated 2026-09-28, round 3: WHOLE-PAPER fix-mode read (not locus-only); PM-verified gate 0 live",
            "papers/group5_qed_gauge/paper_2_alpha.tex": "remediated 2026-09-28, round 3: WHOLE-PAPER fix-mode read (not locus-only); PM-verified gate 0 live",
            "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex": "remediated 2026-09-28, round 3: WHOLE-PAPER fix-mode read (not locus-only); PM-verified gate 0 live",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated 2026-09-28, round 3: WHOLE-PAPER fix-mode read (not locus-only); PM-verified gate 0 live",
        },
    },
    {
        "id": "energy-not-function-of-casimir",
        # Paper 0 abstract said E_n 'is not a function of ... the diagonal Casimir n^2-1 on the labels alone',
        # contradicting its corrected body: on the labels E_n = -Z^2/(2(C+1)), C = n^2-1, i.e. H = -(Z^2/2)(1-Delta)^{-1},
        # the function fixed by Fock's kernel eigenvalue 2pi^2/n (PI-approved 2026-09-29). Registered DELTA #17.
        "pattern": r"not\s+a\s+function\s+of\s+the\s+graph\s+Laplacian,\s+nor\s+of\s+the\s+diagonal\s+Casimir",
        "exempt_if_nearby": withdrawal_marker("energy-not-function-of-casimir"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #17; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
        },
    },
    {
        "id": "framework-reproduces-continuum-cft",
        # Crediting 'the framework' with reproducing continuum CFT-on-sphere observables: Paper 50's F values are
        # continuum zeta determinants the truncation never enters (v5.17.9). Registered DELTA #17.
        "pattern": r"framework's\s+reproduction\s+of\s+continuum\s+CFT|closures:[^.]{0,120}bit-exact\s+CFT-on-sphere",
        "exempt_if_nearby": withdrawal_marker("framework-reproduces-continuum-cft"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #17; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
        },
    },
    {
        "id": "alpha-observation-at-cubic-precision",
        # 'the alpha observation at 8.8e-8' / 'alpha reproduced to 8.8e-8': K matches alpha^-1 to 4.8e-7; 8.8e-8 is the
        # precision of the cubic's root. Registered DELTA #17 (sibling of k-equals-alpha-inverse-to-cubic-precision).
        "pattern": r"\\alpha\$?\s+observation\s+at\s+\$?8\.8|\\alpha\$?\s*---?\s*\$?\\alpha\$?\s+reproduced|\\alpha\$\s+reproduced\s+to\s+\$?8\.8|number---\$\\alpha\$\s+reproduced",
        "exempt_if_nearby": withdrawal_marker("alpha-observation-at-cubic-precision"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #17; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
            "papers/synthesis/group5_qed_gauge_synthesis.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
        },
    },
    {
        "id": "posited-k-boost-class-bw-vacuum",
        # K_alpha^W = 2J_z is POSITED (O0 open); calling it the 'wedge-restricted boost-class generator', 'the chosen
        # K_boost', or the wedge KMS state 'the operator-system analog of the Wightman BW vacuum', or rho_W = e^{-K}/Z
        # 'the BW canonical state', is retired: rho_W is the unit-beta Gibbs state; BW-normalised would be e^{-pi K}/Z
        # (Paper 50 rem:bw_normalisation). Registered DELTA #17.
        "pattern": r"boost-class\s+generator|operator-system\s+analog\s+of\s+the\s+Wightman\s+BW\s+vacuum|chosen\s+K\\?_\\?\{?boost|BW\s+canonical\s+(state|\$?s\s*=|\$?\\beta)",
        "exempt_if_nearby": withdrawal_marker("posited-k-boost-class-bw-vacuum"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #17; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
            "papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
        },
    },
    {
        "id": "p16-spin-resolved-formulation-as-new",
        # Paper 16's minimal-quanta nu_min is prior art (Rittenhouse et al. 2006 Thm 1, proof Cavagnero 1986; Fabre de la
        # Ripelle et al.); the fixed-S S_N-irrep statement is the same count as fixed (n_up, n_down). What is new is the
        # S_N/Schur-Weyl proof route and the atomic tabulation. Registered DELTA #17.
        "pattern": r"contribution\s+is\s+(\S+\s+){0,3}the\s+spin-resolved\s+\$?S_N\$?\s+formulation|narrower:\\?\s*the\s+spin-resolved\s+statement",
        "exempt_if_nearby": withdrawal_marker("p16-spin-resolved-formulation-as-new"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #17; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group4_quantum_computing/paper_16_periodicity.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
            "papers/synthesis/group4_quantum_computing_synthesis.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
        },
    },
    {
        "id": "k-match-inside-theorem",
        # Paper 2 thm:three_homes ('Three homes, structural mystery') placed the empirical K ~ alpha^-1 agreement inside
        # a Theorem environment (touches the sec. 13.5 Observation-only rule). The theorem keeps the three spectral homes
        # and the no-common-generator statement; the approximate match moves to an Observation. Registered DELTA #17.
        "pattern": r"\\begin\{theorem\}\[Three\s+homes,\s+structural\s+mystery\]",
        "exempt_if_nearby": withdrawal_marker("k-match-inside-theorem"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #17; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group5_qed_gauge/paper_2_alpha.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
        },
    },
    {
        "id": "t9-half-integer-claims",
        # Paper 28 T9 is an INTEGER-s theorem: zeta_{D^2}(s) = 2^{2s-1}[lambda(2s-2) - lambda(2s)] in (+)_k pi^{2k} Q at integer
        # s. At half-integer s it carries ODD zeta values (zeta_{D^2}(5/2) = 14 zeta(3) - (31/2) zeta(5), PM-verified to 30
        # digits) and never sqrt(pi); sqrt(pi) lives only in the heat-trace coefficients a_k before division by Gamma(s).
        # Retires 'no odd-zeta content at any s' and 'the sqrt(pi) slot is populated only at half-integer s'. DELTA #20.
        "pattern": r"No\s+odd-zeta\s+content\s+appears\s+at\s+any\s*~?\s*\$s\$|populated\s+only\s+at\s+(?:\\emph\{)?half-integer|slot\s+is\s+populated\s+only\s+at",
        "exempt_if_nearby": withdrawal_marker("t9-half-integer-claims"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #20, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group5_qed_gauge/paper_28_qed_s3.tex": "remediated 2026-09-30 (DELTA #20 idea-level sweeps X/Y; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated 2026-09-30 (DELTA #20 idea-level sweeps X/Y; PM-verified gate 0 live)",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated 2026-09-30 (DELTA #20 idea-level sweeps X/Y; PM-verified gate 0 live)",
        },
    },
    {
        "id": "m1-alone-derives-hawking-beta",
        # IDEA-LEVEL: that the framework's M1 mechanism (alone) produces / reproduces beta = 8 pi M, T_H, beta_U or the shared
        # 2pi. M comes from the Einstein equations and at finite cutoff the generator and beta are POSITED (P34's own honest
        # scope: 'a consistency reproduction of the continuum chain, not a derivation of beta from M1'). Correct: the M1 volume
        # slot re-expresses / is consistent with the continuum value. Registered DELTA #20 (PI-approved).
        "pattern": r"M1[^.]{0,120}?\balone\b|M1\s+mechanism\s+produces\s+\$\\beta|reproduc\w*\s+(?:this\s+)?\$?\\beta(?:_\{?\w+\}?)?\s*=\s*[^$]{0,20}\$?\s+from\s+(?:the\s+)?(?:framework's\s+)?M1|reproduc\w*[^.]{0,40}?\$T_[HU]\b[^.]{0,80}?\bfrom\s+(?:the\s+)?(?:framework's\s+)?M1|reproduc\w*\s+a\s+shared\s+\$2\\pi|from\s+the\s+same\s+M1|M1\s+mechanism[^.]{0,60}?\breproduces",
        "exempt_if_nearby": withdrawal_marker("m1-alone-derives-hawking-beta"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #20, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (DELTA #20 idea-level sweeps X/Y; PM-verified gate 0 live)",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated 2026-09-30 (DELTA #20 idea-level sweeps X/Y; PM-verified gate 0 live)",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex": "remediated 2026-09-30 (DELTA #20 idea-level sweeps X/Y; PM-verified gate 0 live)",
        },
    },
    {
        "id": "torus-rate-cesaro2-p38",
        # Paper 38's kernel is Cesaro-1 (Fejer); its log n/n moment is a kernel property already present on the circle, and Leimbach-van Suijlekom state no torus rate. Retires 'Cesaro-2-style log n/n rate rather than the unweighted 1/n rate of the torus'. FULL #10.
        "pattern": r"rate\s+of\s+(?:the\s+)?torus|unweighted\s+\$1\s*/\s*n\$\s+rate|Ces\\`?aro-2|Cesaro-2",
        "exempt_if_nearby": withdrawal_marker("torus-rate-cesaro2-p38"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex": "remediated 2026-09-30 (FULL #10 mirror check; PM-verified gate 0 live)",
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "p32-connes-placeholder-findings",
        # Paper 32 R3.1 retracts the n_max=2 1:2:3 rational grouping (placeholder artifact) and the 'trend toward correlation' reading (sign flips under physical values); they survived as live 'robust' findings. FULL #10.
        "pattern": r"trend\s+is\s+weakly\s+toward\s+correlation|rational\s+ratios\s+at\s+\$n_\{\\max\}\s*=\s*2\$\)\s+are\s+robust|group\s+into\s+ratios\s+\$1\s*:\s*2\s*:\s*3\$",
        "exempt_if_nearby": withdrawal_marker("p32-connes-placeholder-findings"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "categorically-disjoint-rings",
        # The inner-factor ring 'categorically disjoint' from M1/M2/M3 is the withdrawn strong form ('provably disjoint', 2026-09-03) and false as stated (both contain Q); Observation-level: shares no generator. FULL #10.
        "pattern": r"categorically\s+disjoint\s+from\s+(?:where|the\s+outer|M1|\$?M_?1)|\\emph\{categorically\s+disjoint\}\s+from|forced\s+to\s+a\s+sub-ring\s+categorically",
        "exempt_if_nearby": withdrawal_marker("categorically-disjoint-rings"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
            "papers/synthesis/geovac_field_guide.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "catalan-irrationality-asserted",
        # The irrationality of Catalan's G is open; 'Apery-style irrationality of Catalan G' asserted it and made an M1/M3 separation look proved. FULL #10.
        "pattern": r"irrationality\s+of\s+Catalan|Catalan(?:'s)?\s+(?:\$G\$\s+)?is\s+irrational",
        "exempt_if_nearby": withdrawal_marker("catalan-irrationality-asserted"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "f-as-zeta-of-laplacian",
        # F = pi^2/6 = zeta_{1-Delta}(2) = sum n^2 n^-4; zeta_{-Delta}(2) = sum_{n>=2} n^2/(n^2-1)^2 = 0.885 is not F. Retires 'F equals zeta_{Delta_{S^3}}(2)'. FULL #10.
        "pattern": r"equals\s+\$\\zeta_\{\\Delta_\{S\^3\}\}\(2\)|F\s*=\s*\\zeta_\{\\Delta_\{S\^3\}\}\(2\)",
        "exempt_if_nearby": withdrawal_marker("f-as-zeta-of-laplacian"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "m3-output-carries-pi",
        # M3's closed form beta(s)-beta(s-2) at s=2 gives G - 1/2, with no pi; so 'each mechanism's own output carries pi' is false and the case-exhaustion (iii)=>(i) leg holds only for M1/M2. FULL #10.
        "pattern": r"[Ee]ach\s+mechanism's\s+own\s+output\s+carries|every\s+mechanism[^.]{0,30}carries\s+\$\\pi",
        "exempt_if_nearby": withdrawal_marker("m3-output-carries-pi"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "theorem-grade-non-selection",
        # The Paper 32 VIII non-selection results are all [PROOF-BY-ARGUMENT] (one with an empirical leg), and the Forced-Count Theorem is a dimension count, not a non-selection result. Retires 'theorem-grade non-selection results' / 'fully theorem-bound'. FULL #10.
        "pattern": r"theorem-grade\s+non-selection|at\s+theorem-grade\s+level|fully\s+theorem-bound|\btheorem-bound\s+at",
        "exempt_if_nearby": withdrawal_marker("theorem-grade-non-selection"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_57_forced_free_seam.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
            "papers/synthesis/geovac_field_guide.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "m2-spectral-zeta-unscoped",
        # Paper 28 T9 / Paper 55 M2: spectral-zeta values lie in (+)_k pi^{2k} Q only at INTEGER s (half-integer s >= 5/2 carries odd zeta). The synthesis stated it unscoped. FULL #10.
        "pattern": r"spectral-zeta\s+values\s+it\s+produces\s+lie|zeta\s+values\s+(?:it\s+)?produces?\s+lie\s+in\s+the\s+pure-Tate",
        "exempt_if_nearby": r"integer|" + withdrawal_marker("m2-spectral-zeta-unscoped"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "p38-falsifier-doubled-model",
        # test_p38_action_seminorm.py's T1 checks per-band rank on the DOUBLED model and cannot see the top band; the true-window injectivity premise is backed by test_p38_true_ch_window.py. FULL #10.
        "pattern": r"action\\_seminorm\.py\}\s+backs\s+the\s+injectivity",
        "exempt_if_nearby": withdrawal_marker("p38-falsifier-doubled-model"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "pearson-081-stale",
        # After the 2026-08-29 ERI-evaluator fix the best resolvent Pearson values are 0.846 (n_max=2) and 0.796 (n_max=3), asserted in tests/test_paper31_two_body.py; '<= 0.81' and '0.81 / 0.75' are retired. FULL #10.
        "pattern": r"Pearson[^.]{0,40}?\\le\s*0\.81|Pearson\s+correlation\s+\(\$0\.81\$|\$0\.81\$\s+at\s*\n?\s*\$n_\{\\max\}\s*=\s*2\$,\s*\$0\.75\$",
        "exempt_if_nearby": r"retired|" + withdrawal_marker("pearson-081-stale"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_57_forced_free_seam.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "berezin-mass-distribution-weight",
        # Paper 38: the Berezin map with the mass distribution N/Z over N <= n_max is not unital (B(1) = Z^-1 I); the correct weight is the multiplier symbol sigma(N) over N <= 2n_max-1. Paper 32 still printed the N/Z form as current. FULL #10.
        "pattern": r"Plancherel\s+symbol\s+\$\\hat\{K\}_\{n_\{\\max\}\}\(N\)\s*=\s*N/Z|\\sum_\{N\s*\\le\s*n_\{\\max\}\}\\sum_\{L,M\}\s*\\hat\{K\}",
        "exempt_if_nearby": withdrawal_marker("berezin-mass-distribution-weight"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #10, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (FULL #10 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "case-exhaustion-equivalence-or-complete",
        # Paper 32's pi-source case-exhaustion is [CONDITIONAL]: (ii)=>(iii) argued over the first fifteen Paper 34 projections, (i)<=>(ii) empirical, and the converse (iii)=>(i) WITHDRAWN (M2 at half-integer s and M3 at s=4 give pi-free values). Retires 'the following are equivalent', 'shown equivalent', 'stated and proven', 'the classification is complete', 'every appearance of pi in any ... observable', 'subsumes the prediction'. DELTA #21.
        "pattern": r"finite\s+composition\s+of\s+GeoVac\s+projections[\s\S]{0,250}?following\s+are\s+equivalent|are\s+shown\s+equivalent|stated\s+and\s+proven\s+in\s+Paper~32|classification\s+is\s+complete|Every\s+appearance\s+of\s+\$\\pi\$\s+in\s+any|subsumes\s+the\s+prediction",
        "exempt_if_nearby": withdrawal_marker("case-exhaustion-equivalence-or-complete"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #21, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/synthesis/geovac_field_guide.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "cb-norm-mass-distribution-conflation",
        # The central smoothing maps are UCP, cb-norm 1; (2j+1)/Z is the mass distribution, not the multiplier symbol; 2/(n+1) and 4/[(n_a+1)(n_b+1)] are mass-distribution maxima, not cb-norms. Widens central-multiplier-cb-norm-2-over-n (which could not see the \cbnorm macro). DELTA #21.
        "pattern": r"\\cbnorm\{[^$]{0,40}?\s*(?:\\;)?\s*=\s*(?:\\;)?\s*(?:\\frac\{[24]\}|[24]\s*/)|Plancherel\s+symbol\s+on\s+the\s+Avery|joint\s+Plancherel\s+symbol|joint\s+cb-norms?\s+(?:are|decreases|is\s+one\s+order)",
        "exempt_if_nearby": withdrawal_marker("cb-norm-mass-distribution-conflation"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #21, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "discreteness-is-compactness-biconditional",
        # Retired 2026-09-06: compactness FORCES discreteness; the converse is false (discrete spectrum = compact resolvent; Paper 24's oscillator on non-compact R^3). DELTA #21 found the biconditional live in a synthesis abstract and in Paper 57.
        "pattern": r"discreteness\s+is\s+compactness|\\emph\{Discreteness\}\s*\$\\leftrightarrow\$\s*\\emph\{compactness\}",
        "exempt_if_nearby": withdrawal_marker("discreteness-is-compactness-biconditional"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #21, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/synthesis/group6_precision_observations_synthesis.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_57_forced_free_seam.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "paper7-18-proofs-independent",
        # Paper 7's 18 sympy items are CHECKS of Fock's continuum chain, ten of them limiting-case, definitional or non-discriminating; never 'independent proofs'. DELTA #21.
        "pattern": r"18\s+independent\s+proofs|18\s+symbolic\s+proofs|18/18\s+symbolic\s+proofs|The\s+18\s+symbolic\s+proofs",
        "exempt_if_nearby": withdrawal_marker("paper7-18-proofs-independent"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #21, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
            "papers/group5_qed_gauge/paper_2_alpha.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "f-evaluation-on-scalar-laplacian",
        # F = pi^2/6 = zeta_{1-Delta}(2); 'the Fock-Dirichlet evaluation at s = 2 on the scalar Laplacian' is zeta_{-Delta}(2) = 0.885. Sibling of f-as-zeta-of-laplacian. DELTA #21.
        "pattern": r"evaluation\s+at\s+\$s\s*=\s*2\$\s+on\s+the\s+scalar\s+Laplacian",
        "exempt_if_nearby": withdrawal_marker("f-evaluation-on-scalar-laplacian"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #21, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated 2026-09-30 (DELTA #21 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "code-layer-mass-max-labelled-cb-norm",
        # DELTA #22: the cb-norm correction (central smoothing maps are UCP, cb-norm 1; 2/(n_max+1) and
        # 4/((n_a+1)(n_b+1)) are Plancherel MASS maxima) had not reached tests/ and geovac/. Identifiers
        # containing cb_norm are kept for API stability; prose/docstrings may not call the mass max a cb-norm.
        "pattern": r"cb-norm\s*(?:is\s+)?(?::|=)?\s*(?:2\s*/\s*\(n_max|4\s*/\s*\(\(n_a|2/\(n_max|decreas)|Cesaro-2\s+cb-norm|central-multiplier\s+cb-norm\s*:?\s*(?:2|4)\s*/",
        "exempt_if_nearby": withdrawal_marker("code-layer-mass-max-labelled-cb-norm"),
        "severity": "fail",
        "scope": "trunk group1 group3 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #22, code-layer gate; PI-approved).",
        "files": [
            "tests/*.py",
            "geovac/*.py",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex": "remediated 2026-09-30 (DELTA #22 code-layer sweep; PM-verified gate 0 live)",
            "papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex": "remediated 2026-09-30 (DELTA #22 code-layer sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "p1-ladder-bit-exact",
        # Paper 1: the ladder operators reproduce n (interior shells) and l(l+1) to floating-point tolerance 1e-9 [MEASURED] (the backing test asserts < 1e-9); 'bit-exact' / '[PANEL-VERIFIED] ... exactly' retired. DELTA #24.
        # NARROWED 2026-10-02 (DELTA #25, PI-approved decision B): the identity N = n on interior
        # shells is a two-line [SYMBOLIC PROOF], so 'exactly' is now correct; only the wrong
        # [PANEL-VERIFIED] tier on it stays retired.
        "pattern": r"PANEL-VERIFIED\]\}?\s*(They do, on interior shells|Paper~1~\\cite\{loutey_paper1\}\s+demonstrates)",
        "exempt_if_nearby": withdrawal_marker("p1-ladder-bit-exact"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "higher-l-higher-connectivity",
        # Node degree is not monotone in l on either construction (binary: plateau then drop; CG: 3d m=+-2 3.23 < 3p m=0 5.41). Retires 'higher-l states have higher connectivity'. DELTA #24.
        "pattern": r"higher-\$?l\$?\s+states\s+have\s+higher\s+connectivity|higher-l\s*=\s*higher\s+connectivity",
        "exempt_if_nearby": withdrawal_marker("higher-l-higher-connectivity"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "tests/*.py",
            "geovac/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "door4-forced-paraphrases",
        # Door 4: inner algebra reached only under the Hopf-tower organizing reading plus the adopted Upgrade B axiom; H admitted not forced; seam Observation-level. Retires 'closes the n=2 fork', 'elimination forcings', 'relative necessity', 'seam theorem', 'determines the algebra structure', 'largely/essentially forced', 'forcing of structure', 'packing-reach theorem'. DELTA #24.
        "pattern": r"closes\}?\s+the\s+\$?n\s*=\s*2\$?\s+fork|elimination\s+forcings|relative\s+necessity|is\s+a\s+\\emph\{necessity\}|seam\s+theorem|determines\s+the\s+algebra\s+structure|(largely|essentially)\s+forced|forcing\s+of\s+\\emph\{structure\}|packing-reach\s+theorem",
        "exempt_if_nearby": withdrawal_marker("door4-forced-paraphrases"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "tests/*.py",
            "geovac/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "framework-forces-moduli",
        # Forced-Count: given the admitted A_F (organizing reading + adopted Upgrade B) and the imported CCM representation, the standard axioms fix dim M(D_F); 'the framework forces/determines the moduli space' retired. DELTA #24.
        "pattern": r"framework\s+forces\s+the\s+\\emph\{dimension\}|framework\s+determines\s+what\s+the\s+moduli\s+space|content\s+already\s+proved\s+across\s+the\s+sub-sprints",
        "exempt_if_nearby": withdrawal_marker("framework-forces-moduli"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "three-homes-no-single-generator",
        # The 'no spectral construction generates all three of B, F, Delta' negative is empirical (twelve mechanisms eliminated; none exhibited), withdrawn as a theorem (P32 rem:e6_scope); it may not sit inside a theorem. DELTA #24 (PI-approved).
        "pattern": r"simultaneously\s+generates\s+all\s+three|generates\s+all\s+three\s+at\s+once",
        "exempt_if_nearby": withdrawal_marker("three-homes-no-single-generator"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "k-prediction-label",
        # K = pi(B+F-Delta) is only ever an Observation (CLAUDE.md 13.5): 'This is a prediction' about a K-based consequence is a prohibited tier word. DELTA #24.
        "pattern": r"(\\pi\s*\(B|K\s*=\s*\\pi|combination\s+rule)[\s\S]{0,500}?This\s+is\s+a\s+prediction",
        "exempt_if_nearby": withdrawal_marker("k-prediction-label"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "pi-source-complete-stratification",
        # The pi-source case-exhaustion is [CONDITIONAL] over the first fifteen projections; the period stratification holds over the cases it treats and slot-indexing is a reading. Retires 'complete stratification', 'every transcendental has a motivic home indexed by', 'classifies the transcendentals', 'became a classified program'. DELTA #24.
        "pattern": r"complete\s+(period-theoretic\s+|period-ring\s+)?stratification|[Ee]very\s+transcendental\s+(that\s+has\s+appeared|observation)[^.]{0,120}(motivic\s+home|indexed\s+by)|classif(ies|y)\s+the\s+transcendentals|which\s+transcendentals\s+can\s+appear|became\s+a\s+classified\s+program",
        "exempt_if_nearby": withdrawal_marker("pi-source-complete-stratification"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "m1-governs-gh-rate",
        # Paper 38: the appearance of 4/pi in the GH rate is a reading, not a derivation; M1 does not govern or produce the rate. DELTA #24.
        "pattern": r"M1\s+governs\s+state-space|identifies\s+as\s+the\s+\$?M_?\{?1|appears\s+as\s+the\s+leading\s+constant\s+in\s+the\s+GeoVac\s+state-space",
        "exempt_if_nearby": withdrawal_marker("m1-governs-gh-rate"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "fk-on-s3-density-overlap",
        # Paper 7: the S^3 density-overlap formula is the Fourier-Parseval form of F^0 for s-orbital pairs only [PRIOR ART]; general single-centre F^k are exact rationals by a separate closed form, not 'on S^3'. DELTA #24.
        "pattern": r"Slater\s+\$?F\^k\$?\s+integrals\s+on\s+\$S\^3\$|\$S\^3\$\s+density\s+overlap\s+formula|exact\s+rational\s+\$S\^3\$\s+density\s+overlap|algebraic\s+on\s+\$S\^3\$|Slater\s+\$F\^0\$\s+integrals\s+on\s+the\s+\$S\^3\$",
        "exempt_if_nearby": r"\$s\$-orbital pairs|s-orbital pairs|Parseval|" + withdrawal_marker("fk-on-s3-density-overlap"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "he-0004-called-variational",
        # Paper 13: 0.004% is a Schwartz cusp EXTRAPOLATION lying marginally below exact (non-variational); 0.022% raw is the variational 2D figure. Retires 0.004% presented as a variational / '2D variational' result. DELTA #24.
        "pattern": r"0\.004\\?%--0\.022\\?%|2D\s+variational\s+solver\s+achieves\s+0\.004|0\.004\\?%[^\n]{0,90}via\s+the\s+2D\s+variational\s+solver|He\s+0\.004\\?%,",
        "exempt_if_nearby": r"non-variational|extrapolat|below the exact|below exact|" + withdrawal_marker("he-0004-called-variational"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "p38-almost-commutative-proved",
        # Paper 38 proves the COMMUTATIVE factor only; the almost-commutative / gauge-network convergence is OPEN; the retired display Lambda <= C_3 gamma is superseded by d <= gamma. DELTA #24.
        "pattern": r"is\s+proved\s+at\s+the\s+state-space\s+Gromov|\\le\s*C_\{?3\}?\s*(\\cdot\s*)?\\gamma|<=\s*C_3\s*\*\s*gamma|Theorem\s+L5\s+\(GH\s+convergence",
        "exempt_if_nearby": r"withdrawn|refuted|retired|REFUTED|" + withdrawal_marker("p38-almost-commutative-proved"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "tests/*.py",
            "geovac/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "uniform-bound-for-all-n",
        # gamma_n <= 6 log n / n is a finite-range [MEASURED] check, not proved for all n_max. Retires 'holds for all n >= 2' and 'pinned rigorously'. DELTA #24.
        "pattern": r"6\s*\\?,?\s*\(?\\?log[^.]{0,30}?/\s*n[^.]{0,40}?(holds\s+)?for\s+all\s+\$?n|6\s*log\s*\(?n\)?\s*/\s*n\s+for\s+all\s+n|pinned\s+the\s+asymptote\s+rigorously",
        "exempt_if_nearby": r"not proved|finite[- ]range|checked range|MEASURED|" + withdrawal_marker("uniform-bound-for-all-n"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "tests/*.py",
            "geovac/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "p40-general-g-settled",
        # Paper 40 proves general-G convergence only CONDITIONALLY on the per-group spin-window decomposition (unconditional at rank 1); scalar compact-group convergence is Gaudillot-Estrada--van Suijlekom's. Retires 'resolved by Paper 40', 'settled/closed by the present author', 'closed it on the universal family', 'to all simple compact Lie groups' unscoped. DELTA #24.
        "pattern": r"has\s+been\s+resolved\s+by\s+Paper~?40|first\s+non-abelian[^.]{0,80}(closed|settled)\s+by\s+the\s+present\s+author|closed\s+it\s+on\s+the\s+universal|machinery\s+to\s+all\s+simple\s+compact\s+Lie\s+groups",
        "exempt_if_nearby": withdrawal_marker("p40-general-g-settled"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "four-pi-is-unit-quotient",
        # 4/pi is TWICE the unit-radius quotient Vol(S^2)/Vol(S^3) = 2/pi. Retires '4/pi ... a quotient of unit-radius sphere volumes'. DELTA #24.
        "pattern": r"constant\s+is\s+a\s+quotient\s+of\s+(\\emph\{)?unit-radius|\\pi\^2\$?,\s*a\s+quotient\s+of\s+unit-radius",
        "exempt_if_nearby": withdrawal_marker("four-pi-is-unit-quotient"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "k-within-1e-7",
        # K itself matches 1/alpha to 4.8e-7; only the cubic's root is within 1e-7. DELTA #24.
        "pattern": r"remain\s+within\s+\$10\^\{-7\}\$",
        "exempt_if_nearby": withdrawal_marker("k-within-1e-7"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "cvs-outlook-misattribution",
        # Connes-van Suijlekom 2021's Outlook does not discuss Peter-Weyl on SU(2) or S^3 GH convergence; they defer GH questions generically ('will be studied elsewhere', Sec. 3). DELTA #24 citation review.
        "pattern": r"left\s+open\s+in\s+the\s+outlook\s+of\s+Connes",
        "exempt_if_nearby": withdrawal_marker("cvs-outlook-misattribution"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {},
    },
    {
        "id": "code-retired-p38-labels",
        # Code/docs layer: GH bound is not a basis-truncation error estimate; rate tier is INTERNAL THEOREM + MEASURED, not 'proven unconditionally'; 'Hopf-base signature' label retired. DELTA #24.
        "pattern": r"basis-truncation\s+error\s+estimate|rate[^.\n]{0,40}is\s+proven\s+unconditionally|Hopf-base\s+(measure\s+)?signature|Hopf[- ]base\s+measure",
        "exempt_if_nearby": r"not\s+a\s+basis-truncation|retired|misnomer|withdrawn|" + withdrawal_marker("code-retired-p38-labels"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #24, idea-level gate; PI-approved).",
        "files": [
            "tests/*.py",
            "geovac/*.py",
            "docs/validation_benchmarks.md",
        ],
        "cited_by": {},
    },
    {
        "id": "p38-uniform-bound-panel-verified",
        # Paper 38 Lemma L2(d.iii): gamma_n <= 6 log n / n is a FINITE-RANGE monotonicity check, [MEASURED], not proved for all n_max. Retires restatements tagging it [PANEL-VERIFIED] or stating it for all n_max >= 2. FULL #11 self-check.
        "pattern": r"6\s*(?:\\,)?\s*\\log\s*\\?n(?:max|_\{\\max\})[^.;]{0,80}?PANEL-VERIFIED|PANEL-VERIFIED\]\}\s*A\s+uniform\s+bound",
        "exempt_if_nearby": withdrawal_marker("p38-uniform-bound-panel-verified"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #11 pre-checkpoint self-check; idea-level gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (FULL #11 self-check sweep; PM-verified gate 0 live)",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated 2026-09-30 (FULL #11 self-check sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "p1-higher-l-more-connections",
        # Paper 1: on the binary label lattice the node degree depends on m and on which shell neighbours exist, not monotonically on l (n=3 interior: l=0,1,2 -> 2, 3-4, 2-3); only D_{2p,m=0}=3 > D_{2s}=2 holds. Retires 'higher l => more connections' / 'high-l states occupy higher-degree nodes'. FULL #11.
        "pattern": r"higher\s+angular\s+momentum\s+\$?l?\$?\s*have\s+\\textit\{more\s+connections|high\s+angular\s+momentum\s+states\s+occupy\s+higher-degree|higher-l\s+(states|nodes)\s+have\s+more\s+connections",
        "exempt_if_nearby": withdrawal_marker("p1-higher-l-more-connections"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #11, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
        ],
        "cited_by": {
        },
    },
    {
        "id": "p7-he-bracket-denied",
        # Paper 7: the fixed-k CI (1.88% above exact) and the cusp-corrected 0.004% (a Schwartz extrapolation lying marginally BELOW exact, Paper 13) DO bracket the exact He energy. Retires 'both methods sit above the exact value, so they do not bracket it'. FULL #11.
        "pattern": r"do\s+not\s+bracket\s+it|both\s+methods\s+sit\s+above\s+the\s+exact",
        "exempt_if_nearby": withdrawal_marker("p7-he-bracket-denied"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #11, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
        ],
        "cited_by": {
        },
    },
    {
        "id": "p7-parseval-credited-to-s3",
        # Paper 7: the single-integral F^0 master formula is the standard Fourier-Parseval form of the Coulomb integral, rewritten in the stereographic variable t = q/(2Z) [PRIOR ART]; the single-integral structure is Parseval's, not a consequence of S^3 geometry. Retires restatements crediting the Fock projection / S^3 with it. FULL #11.
        "pattern": r"key\s+structural\s+finding\s+is\s+that\s+Eq\.?~?\\eqref\{eq:f0_s3\}|Fock\s+projection\s+converts\s+position-space\s+double|derives\s+the\s+two-electron\s+repulsion\s+integral\s+as\s+a\s+density-overlap\s+integral\s+on",
        "exempt_if_nearby": withdrawal_marker("p7-parseval-credited-to-s3"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #11, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
        ],
        "cited_by": {
        },
    },
    {
        "id": "inner-algebra-mostly-forced",
        # Paper 32 Door 4: the inner algebra C+H+M3(C) is an organizing reading (Bertrand x complex-Hopf tower, n<=3) plus the adopted Upgrade B axiom: factor count and the C, M3(C) summands are FIXED UNDER THAT READING, H is admitted but not forced, the rings share no generator (no 'disjointness'), and the seam is Observation-level. Retires 'mostly forced', 'relative theorem of necessity', 'disjointness forbids selection', 'the one inner datum the construction does reach'. FULL #11.
        "pattern": r"mostly\s+forced|relative\s+theorem\s+of\s+necessity|(?:disjointness|seam)(?:\s+already)?\s+forbids\s+selection|inner\s+datum\s+the\s+construction\s+does\s+reach|inner\s+algebra\s+structure\s+is\s+(?:mostly\s+)?forced",
        "exempt_if_nearby": withdrawal_marker("inner-algebra-mostly-forced"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk FULL #11, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_57_forced_free_seam.tex": "remediated 2026-09-30 (FULL #11 sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "pi-source-converse-stated-false",
        # The pi-source converse (iii)=>(i) is WITHDRAWN as NOT ESTABLISHED: no mechanism output is PROVED pi-free (B = 42's M1 tag is a Paper 55 convention). Retires 'does not force pi', 'not every output of a mechanism carries pi', 'B = 42 ... proven witness'. DELTA #23.
        "pattern": r"does\s+not\s+force\s+\$?\\?pi|not\s+every\s+output\s+of\s+a\s+mechanism\s+carries|proven\s+witness|B\s*=\s*42\s*\n?\s*\(proven\)",
        "exempt_if_nearby": withdrawal_marker("pi-source-converse-stated-false"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #23, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "tests/*.py",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-30 (DELTA #23 sweep; PM-verified gate 0 live)",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex": "remediated 2026-09-30 (DELTA #23 sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated 2026-09-30 (DELTA #23 sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "lorentzian-hybrid-labelled-mass-max",
        # Paper 45: the SU(2) x U(1) value 2/(n_max+1)*1 pairs the SU(2) MASS maximum with the U(1) SYMBOL maximum -- a hybrid, neither a mass maximum (2/(n_max+1)/N_t) nor a cb-norm (1). DELTA #22 had mislabelled it 'joint mass maximum'. DELTA #23.
        "pattern": r"mass_max\(U\(1\)\)|[Jj]oint\s+mass\s+maximum\s*(?:=|\(not)",
        "exempt_if_nearby": withdrawal_marker("lorentzian-hybrid-labelled-mass-max"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #23, idea-level gate; PI-approved).",
        "files": [
            "tests/*.py",
            "geovac/*.py",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex": "remediated 2026-09-30 (DELTA #23 sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "bozejko-fendler-as-support",
        # Bozejko-Fendler (Arch. Math. 57, 1991) concerns DISCRETE groups; it supports neither the compact-group central kernel, nor a cb-norm equality, nor the C_3 <= 2 Leibniz bound, nor a UCP partial inverse. Withdrawn as support corpus-wide (P32/P40/P39/group1 synthesis/code). DELTA #23.
        "pattern": r"Bo\S{0,6}?ejko--?Fendler\s+(?:central\s+kernel|estimate\b|bound\b|symbol-side|product\s+cb-norm\s+formula|central-multiplier\s+equality)|generic\s+Bozejko-Fendler|conservative\s+Bozejko-Fendler",
        "exempt_if_nearby": r"withdrawn|mislabel|historical|discrete groups|retired|" + withdrawal_marker("bozejko-fendler-as-support"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-30 (/qa trunk DELTA #23, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "tests/*.py",
            "geovac/*.py",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex": "remediated 2026-09-30 (DELTA #23 sweep; PM-verified gate 0 live)",
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex": "remediated 2026-09-30 (DELTA #23 sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "hopf-base-measure-label",
        # The M1 mechanism's name 'Hopf-base measure' is a misnomer: M1 is the k=0 VOLUME SLOT of the master
        # Mellin engine (Paper 18 sec. III.7), not a Hopf base-to-total ratio, and the 4/pi GH-rate constant is the
        # unit-sphere volume quotient (convention-dependent), not a Hopf-base quantity. Swept corpus-wide 2026-09-29
        # (77 loci, 16 papers); one explained historical mention in Paper 18 carries the marker. Genuine geometry
        # ('a quarter of the area of the Hopf base S^2') stays silent. Registered DELTA #19 (PI-approved).
        "pattern": r"Hopf[- ](?:base|bundle)\s+(?:measure|volume)|Hopf[- ]measure|quarter\s+of\s+the\s+area\s+of\s+the\s+Hopf\s+base|Hopf[- ]base[- ](?:Weyl|coefficient|calibration|family|constant|measure)|M1\s+Hopf[- ]base|Hopf/\$?S\^\{?3\}?\$?-measure|Hopf\s+bundle[^.]{0,60}?produces\s+\$?\\(?:Vol|mathrm\{Vol\})|Hopf[- ]bundle\s+topological\s+measure|Hopf-base\s+solid-angle",
        "exempt_if_nearby": withdrawal_marker("hopf-base-measure-label"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #19 label sweep; PI-approved). Docs lookup/vocabulary tables added 2026-10-15 after a survivor was found in topic_to_paper_lookup.md (Mellin single-sourcing sweep). Code layer (tests/, geovac/) added 2026-10-16 after the claim-search DELTA found nine comment/docstring survivors.",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "docs/topic_to_paper_lookup.md",
            "docs/vocabulary_translation.md",
            "docs/forced_free_seam.md",
            "docs/atlas.md",
            "tests/**/*.py",
            "geovac/**/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "framework-2pi-is-boost-period",
        # IDEA-LEVEL: any sentence whose subject is the framework / its projection / M1 / 'this 2pi' and which says it IS the
        # modular-flow / boost-orbit period or the Wick-rotated image of the BW mechanism. The framework's M1 2pi MATCHES that
        # period (an analogy); identifying the framework's generator with a boost is OPEN (Paper 42 O0). Continuum statements
        # (Sewell: the cigar KMS state IS the boost modular flow; Unruh) are correct and stay silent. Registered DELTA #19.
        "pattern": r"(?:framework's|projection's|\bM1\b|this\s+\$2\\pi\$)[^.;]{0,140}?\bIS\s+(?:the\s+)?(?:Lorentzian\s+)?(?:period\s+of\s+the\s+)?(?:modular|Lorentz|boost|Wick)|is\s+the\s+Wick-rotated\s+image\s+of\s+the\s+(?:Bisognano|modular)|\\textbf\{The\s+structural\s+identification\.\}",
        "exempt_if_nearby": withdrawal_marker("framework-2pi-is-boost-period"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #19, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-29 (DELTA #19 idea-level sweep; PM-verified gate 0 live)",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated 2026-09-29 (DELTA #19 idea-level sweep; PM-verified gate 0 live)",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex": "remediated 2026-09-29 (DELTA #19 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "sqrtpi-ring-anywhere",
        # IDEA-LEVEL: the two-term ring sqrt(pi) Q (+) pi^2 Q anywhere. At integer s the values are pure-Tate (+)_k pi^{2k} Q
        # (Paper 28 T9); heat-kernel coefficients are sqrt(pi) x (+) pi^{2k} Q (Paper 18 sec. III.7). The one correctly scoped
        # mention (Paper 28's half-integer-s remark) carries a withdrawal marker. Registered DELTA #19.
        "pattern": r"\\sqrt\{?\\pi\}?\s*(?:\\cdot|\\,)?\s*(?:\\mathbb\{Q\}|\\Q)\s*\\oplus\s*\\pi\^\{?2\}?",
        "exempt_if_nearby": withdrawal_marker("sqrtpi-ring-anywhere"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #19, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-09-29 (DELTA #19 idea-level sweep; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated 2026-09-29 (DELTA #19 idea-level sweep; PM-verified gate 0 live)",
        },
    },
    {
        "id": "posited-k-is-bw-local-hamiltonian",
        # 'D_W is NOT the right local Hamiltonian for the BW vacuum; H_local := K_alpha/beta is' asserts the posited
        # generator IS the BW vacuum's Hamiltonian; its beta = 2pi Gibbs state is rho_W, and whether that is a BW vacuum
        # is OPEN (Paper 42 O0). Registered DELTA #18.
        "pattern": r"right\s+local\s+Hamiltonian\s+for\s+the\s+BW\s+vacuum",
        "exempt_if_nearby": withdrawal_marker("posited-k-is-bw-local-hamiltonian"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #18; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #18; PM-verified gate 0 live)",
        },
    },
    {
        "id": "energies-not-output-of-casimir",
        # 'The energies are not an output of these operators' next to the SO(4) Casimir (and P1's 'Casimir built from
        # them ... not an output') retires the Casimir from the energies; on the labels E_n = -Z^2/(2(C+1)) through the
        # function Fock's kernel fixes. Registered DELTA #18 (sibling of energy-not-function-of-casimir).
        "pattern": r"not\s+an\s+output\s+of\s+these\s+operators|imposed\s+energy-shell\s+map\s+of\s+\$n\$,\s+not\s+an\s+output\.",
        "exempt_if_nearby": withdrawal_marker("energies-not-output-of-casimir"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #18; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #18; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_1_spectrum.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #18; PM-verified gate 0 live)",
        },
    },
    {
        "id": "k-defined-inside-theorem",
        # A Theorem environment that opens by defining K = pi(B+F-Delta) (P18 thm:k_decomposition) puts the combination
        # rule under a Theorem label; state the three ingredients and keep K (an Observation) outside. Registered DELTA #18
        # (twin of k-match-inside-theorem).
        "pattern": r"\\begin\{theorem\}[^\n]*\n(?:\\label\{[^}]*\}\s*\n)?\s*Let\s+\$K\s*=",
        "exempt_if_nearby": withdrawal_marker("k-defined-inside-theorem"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #18; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #18; PM-verified gate 0 live)",
        },
    },
    {
        "id": "framework-is-boost-physics",
        # 'the temperature mechanism IS already Lorentz-boost physics' / 'already does (more) boost-class physics': the
        # framework's posited generator shares the 2pi period of the BW boost orbit; identifying it with a boost is OPEN
        # (Paper 42 O0). Registered DELTA #18.
        "pattern": r"IS\s+already\s+Lorentz-boost|already\s+does\s+(more\s+)?(Lorentz-)?boost-class\s+physics|already\s+captures\s+the\s+Lorentz-boost",
        "exempt_if_nearby": withdrawal_marker("framework-is-boost-physics"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #18; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #18; PM-verified gate 0 live)",
        },
    },
    {
        "id": "m2-ring-sqrtpi-at-integer-s",
        # M2 at integer s is the pure-Tate ring (+)_k pi^{2k} Q; the sqrt(pi) slot is populated only at half-integer s
        # (Paper 28 remark; Paper 18 sec. III.7: overall sqrt(pi) x (+) pi^{2k} Q for heat-kernel coefficients).
        # 'zeta_{D^2}(s) in/constrained to sqrt(pi) Q (+) pi^2 Q' at integer s is false (zeta_{D^2}(2) = pi^2 - pi^4/12).
        # Registered DELTA #18.
        "pattern": r"(constrained\s+to|\\in|coefficients|residual\s+is)\s+\$?\\sqrt\{?\\pi\}?\s*\\cdot\s*(\\mathbb\{Q\}|\\Q)\s*\\oplus\s*\\pi\^\{?2\}?\s*\\cdot\s*(\\mathbb\{Q\}|\\Q)",
        "exempt_if_nearby": withdrawal_marker("m2-ring-sqrtpi-at-integer-s"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #18; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group5_qed_gauge/paper_28_qed_s3.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #18; PM-verified gate 0 live)",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #18; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_57_forced_free_seam.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #18; PM-verified gate 0 live)",
        },
    },
    {
        "id": "cubic-root-transcendence-understated",
        # The cubic's root alpha_K is PROVABLY transcendental: alpha^3 - K alpha + 1 = 0 gives K = alpha^2 + 1/alpha,
        # so algebraic alpha_K would make K algebraic; but K = (1679/40) pi + pi^3/6 is a nonconstant Q-polynomial in pi
        # (Lindemann). The earlier 'likely (though not proven)' understated it, and its counterexample
        # x^2 - pi x + pi(pi-1) = 0 is false at x = 1 ((pi-1)^2 != 0). Nothing follows for the physical alpha.
        # Registered 2026-09-29 (DELTA #17 self-check; PI-delegated).
        "pattern": r"\\pi\s*\(\s*\\pi\s*-\s*1\s*\)\s*=\s*0|\(though\s+not\s+proven\)|root\s+is\s+likely\s+transcendental",
        "exempt_if_nearby": withdrawal_marker("cubic-root-transcendence-understated"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (DELTA #17 self-check).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated 2026-09-29 (DELTA #17 self-check; one-line Lindemann proof, PM-verified; gate 0 live)",
            "papers/group5_qed_gauge/paper_2_alpha.tex": "remediated 2026-09-29 (DELTA #17 self-check; one-line Lindemann proof, PM-verified; gate 0 live)",
        },
    },
    {
        "id": "p50-log2-as-m2",
        # Paper 50 labelled log 2 an M2 (Seeley-DeWitt) transcendental and its F_D + 2F_s = log(2)/2 combination
        # 'the M2 sector'. Paper 18's period-ring stratification puts M2 inside (+)_k pi^{2k} Q (pure Tate, no log 2,
        # no zeta(3)); log 2 is a level-2 cyclotomic period (M3's ring) and zeta(3)/pi^2 is M1 x M3, as Paper 55 had it.
        # In Paper 50's own Dirac computation log 2 enters through the 2^s factor of the half-integer Hurwitz
        # identity (its own M3 mechanism). Registered 2026-09-29 (PI-approved: adopt Paper 18's stratification).
        "pattern": r"isolates\s+(the\s+)?M2\b|\\log\s*2\$?\s*\(an\s+M2\s+transcendental|M2--M3\s+orthogonal\s+decomposition|orthogonal\s+(\$?M_?\{?2\}?\$?)\s*/\s*\$?M_?\{?3\}?\$?\s+(dual-basis\s+)?decomposition|M_\{2\}\$/\$M_\{3\}\$\s+orthogonal\s+decomposition|log\(2\)/2\$?\s*\(M2",
        "exempt_if_nearby": withdrawal_marker("p50-log2-as-m2"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (M-label reconciliation, PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_50_cft3_partition_function.tex": "remediated + summary surfaces reread 2026-09-29 (M-label reconciliation per Paper 18 period-ring stratification; PM-verified S5 value; gate 0 live)",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated + summary surfaces reread 2026-09-29 (M-label reconciliation per Paper 18 period-ring stratification; PM-verified S5 value; gate 0 live)",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated + summary surfaces reread 2026-09-29 (M-label reconciliation per Paper 18 period-ring stratification; PM-verified S5 value; gate 0 live)",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-29 (M-label reconciliation per Paper 18 period-ring stratification; PM-verified S5 value; gate 0 live)",
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "remediated + summary surfaces reread 2026-09-29 (M-label reconciliation per Paper 18 period-ring stratification; PM-verified S5 value; gate 0 live)",
        },
    },
    {
        "id": "p53-bochner-riesz-positivity-and-apex-backbone",
        # Paper 53: (1) the 2D Bochner-Riesz kernel J_{s+1}(r)/r^{s+1} oscillates, so it is negative somewhere
        # at EVERY finite order (min/max -0.07, -0.06, -0.03, -0.009 at s = 0.75, 1, 2, 4) and the BR map is
        # NOT positivity-preserving -- the paper's own rem:height_constant said so while thm:plane and the
        # abstract claimed classical positivity. (2) The 'plane' is the flat radial sector (order-0 Hankel);
        # the cone R^2_alpha has no translation group, so no apex backbone for Paper 51 is supplied (and
        # Paper 51 never cited it). Registered 2026-09-29 (PI-approved rebuild).
        "pattern": r"Bochner--?Riesz\S*\s+(\S+\s+){0,6}positivity-preserving|Bochner--?Riesz\s+means\s+are\s+positivity|kernel\s+is\s+positive\s+for\s+s\s*>=\s*1/2|propinquity\s+backbone\s+for\s+the\s+cigar-tip|upgrading\s+its\s+Lemma\s+on\s+geometric\s+faithfulness|supplies\s+the\s+propinquity\s+backbone|upgraded\s+from\s+cited\s+to\s+proved\s+for\s+this\s+purpose|obstructed\s+by\s+(its|the)\s+Dirichlet\s+boundary|positivity\s+and\s+contractivity\s+at\s+Ces|which\s+the\s+interior\s+L4\s+supplies",
        "exempt_if_nearby": withdrawal_marker("p53-bochner-riesz-positivity-and-apex-backbone"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (Paper 53 rebuild; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "tests/test_paper53_plane_bochner_riesz.py",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_53_disk_propinquity.tex": "remediated + summary surfaces reread 2026-09-29 (P53 rebuild; P51 lem:L6_backbone -> OPEN; PM-reviewed; gate 0 live)",
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "remediated + summary surfaces reread 2026-09-29 (P53 rebuild; P51 lem:L6_backbone -> OPEN; PM-reviewed; gate 0 live)",
            "papers/group5_qed_gauge/paper_51_gravity_arc.tex": "remediated + summary surfaces reread 2026-09-29 (P53 rebuild; P51 lem:L6_backbone -> OPEN; PM-reviewed; gate 0 live)",
        },
    },
    {
        "id": "p50-kps-from-truncated-triple",
        # Paper 50's 'bit-exact KPS' theorems sum the FULL continuum conformal-scalar and CH Dirac spectra
        # (n = 0..inf) with standard zeta regularisation -- exactly Klebanov-Pufu-Safdi's computation. The
        # finite truncation never enters, and state-space GH convergence plays no role in the F values, so
        # 'reproduced from the discrete truncated triple' and the Category III 'correspondence mechanism'
        # thesis (Paper 52, archived) are retired. Registered 2026-09-29 (PI-approved correct-in-place).
        "pattern": r"reproduction\s+from\s+the\s+discrete\s+truncated|discrete\s+spectral\s+triple\s+reproduces|from\s+the\s+spectra\s+of\s+the\s+truncated\s+scalar|spectral\s+asymptotics\s+of\s+the\s+discrete\s+spectral\s+triple|F-theorem\s+coefficients\s+on\s+the\s+truncated|Category\s+III''?\s+discretization\s+correspondence|State-space\s+GH\s+convergence\s+at\s+rate\s+\$4/\\pi\$\}",
        "exempt_if_nearby": withdrawal_marker("p50-kps-from-truncated-triple"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (Paper 50 correct-in-place + Paper 52 archive; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_50_cft3_partition_function.tex": "remediated + summary surfaces reread 2026-09-29 (P50 correct-in-place; PM-reviewed abstract + not-holographic para; gate 0 live)",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated + summary surfaces reread 2026-09-29 (P50 correct-in-place; PM-reviewed abstract + not-holographic para; gate 0 live)",
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "remediated + summary surfaces reread 2026-09-29 (P50 correct-in-place; PM-reviewed abstract + not-holographic para; gate 0 live)",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + summary surfaces reread 2026-09-29 (P50 correct-in-place; PM-reviewed abstract + not-holographic para; gate 0 live)",
        },
    },
    {
        "id": "fock-delta-overcorrection",
        # v5.17.6/7 OVERCORRECTED DELTA #15's genuine Green's-function defect into 'not a free-particle
        # problem' / 'Delta_{S^3} is not the link'. Exact statement (PI-approved 2026-09-29): hydrogen's bound
        # states are the S^3 harmonics (Delta_{S^3} eigenfunctions); H_bound = -(Z^2/2)(1-Delta_{S^3})^{-1}, i.e.
        # E_n = -Z^2/(2(C+1)), C = n^2-1; the kernel ~ (1-Delta)^{-1/2} (eigenvalue 2pi^2/n) fixes p0 = Z/n; the
        # kernel is NOT Delta's Green's function. Registered 2026-09-29, DELTA #16.
        "pattern": r"not\s+a\s+free-particle|is\s+not\s+what\s+links|not,\s+by\s+itself,\s+the\s+link",
        "exempt_if_nearby": withdrawal_marker("fock-delta-overcorrection"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #16 family sweeps; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 family sweeps; PM-verified gate 0 live corpus-wide)",
            "papers/group3_foundations/paper_24_bargmann_segal.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 family sweeps; PM-verified gate 0 live corpus-wide)",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 family sweeps; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "k-equals-alpha-inverse-to-cubic-precision",
        # K = pi(B+F-Delta) = 137.036064 vs alpha^-1 = 137.035999: relative gap 4.8e-7. 8.8e-8 is the accuracy
        # of the CUBIC root (1/alpha = K - alpha^2), not of K. Registered 2026-09-29, DELTA #16.
        "pattern": r"\(B\s*\+\s*F\s*-\s*\\Delta\)[^.]{0,200}equals\s+\$?\\alpha\^\{-1\}\$?\s+to\s+\$?8\.8|combination\s+rule(?:(?!cubic|root|4\.8|4\.77)[^.]){0,160}?8\.8|B\s*\+\s*F\s*-\s*\\Delta\s*(?:\\right)?\)(?:(?!cubic|root|4\.8|4\.77)[^.]){0,200}?(?:reproduces|matches|match)\s(?:(?!cubic|root|4\.8|4\.77)[^.]){0,80}?8\.8|Delta\s*\)\s*=\s*\\alpha\^\{-1\}",
        "exempt_if_nearby": withdrawal_marker("k-equals-alpha-inverse-to-cubic-precision"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #16 family sweeps; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 family sweeps; PM-verified gate 0 live corpus-wide)",
            "papers/group5_qed_gauge/paper_2_alpha.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 self-check remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group5_qed_gauge/paper_25_hopf_gauge_structure.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 self-check remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 self-check remediation; PM-verified gate 0 live corpus-wide)",
            "papers/synthesis/geovac_field_guide.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 self-check remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 self-check remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "bbb-mn-signature-mislabel",
        # Bizi-Brouder-Besnard (arXiv:1611.07062 p.3): m = t+s, n = t-s (mod 8). Riemannian S^3 (t,s)=(0,3) is
        # (m,n) = (3,5); Lorentzian (t,s)=(1,3) is (4,6). Writing the signature as '(m,n) = (3,0)' / '(3,1)'
        # conflates (s,t) with BBB's (m,n). Registered 2026-09-29, DELTA #16.
        "pattern": r"\(m,\s*n\)\s*=\s*\(3,\s*[01]\)",
        "exempt_if_nearby": withdrawal_marker("bbb-mn-signature-mislabel"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-09-29 (/qa trunk DELTA #16 family sweeps; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 family sweeps; PM-verified gate 0 live corpus-wide)",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 family sweeps; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "coulomb-kernel-as-laplacian-green-function",
        # Paper 7 called 1/|n-n'|^2 on S^3 'the Green's function of the Laplace-Beltrami operator'. It is
        # not: its proved eigenvalue 2pi^2/n (tests/test_paper7_fock_kernel.py) makes it proportional to
        # (1 - Delta_{S^3})^{-1/2}; a Green's function of Delta_{S^3} goes like 1/d on a 3-manifold and has
        # eigenvalues ~1/(n^2-1). (It is the restriction of the R^4 Laplacian's Green's function.)
        # Registered 2026-09-28, DELTA #15.
        "pattern": r"Green's\s+function\s+for\s+the\s+Laplace|free-particle\s+Green's\s+function\s+on",
        "exempt_if_nearby": withdrawal_marker("coulomb-kernel-as-laplacian-green-function"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #15; PI-approved remediation).",
        "files": [
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "finite-truncation-even-unscoped",
        # Only Paper 32's chirality-doubled MODEL truncation is even; the TRUE Camporesi-Higuchi window
        # that converges is ODD like the continuum (commutant of O = scalars). Unscoped 'the finite
        # truncation is even' is false for the converging object. Registered 2026-09-28, DELTA #15.
        "pattern": r"(the\s+)?finite\s+truncation\s+is\s+even|cutoff\s+the\s+truncation\s+is\s+even|J\$?-relations;\s+an\s+even",
        "exempt_if_nearby": withdrawal_marker("finite-truncation-even-unscoped"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #15; PI-approved remediation).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "dw-off-diagonal-wedge-basis",
        # The wedge Dirac D_W is chirality-DIAGONAL and sign-paired, not off-diagonal (Paper 42 §11.7,
        # proof of Thm 11.14). The old 'diagonal/off-diagonal subspace decomposition' mechanism is
        # withdrawn; the orthogonality is a chirality-parity identity. Registered 2026-09-28, DELTA #15.
        "pattern": r"D_W\^\{?L\}?\$?\s+(is\s+)?off-diagonal",
        "exempt_if_nearby": withdrawal_marker("dw-off-diagonal-wedge-basis"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #15; PI-approved remediation).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "docs/topic_to_paper_lookup.md",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_24_bargmann_segal.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
            "docs/topic_to_paper_lookup.md": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "central-multiplier-cb-norm-2-over-n",
        # The central-multiplier cb-norm is ||T_K||_cb = 1; 2/(n_max+1) is the mass-distribution maximum,
        # not a cb-norm (Paper 32 rem:gh_status_caveat (i); Paper 38). Registered 2026-09-28, DELTA #15.
        "pattern": r"_\{\\mathrm\{cb\}\}\s*=\s*2/",
        "exempt_if_nearby": withdrawal_marker("central-multiplier-cb-norm-2-over-n"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #15; PI-approved remediation).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "lvs-torus-rate-attribution",
        # Leimbach-van Suijlekom (Adv. Math. 439, 2024) state NO explicit rate (Thm 3.10: gamma -> 0 only);
        # the log is already in the circle's Fejer kernel. Registered 2026-09-28, DELTA #15.
        "pattern": r"Leimbach--van~?Suijlekom\s+torus\s+rate",
        "exempt_if_nearby": withdrawal_marker("lvs-torus-rate-attribution"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #15; PI-approved remediation).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "posited-k-as-truncated-bw-boost",
        # The finite K = diag(2 m_j) is POSITED; calling it 'the truncated Bisognano-Wichmann boost
        # generator' asserts the OPEN identification (Paper 42 O0). Registered 2026-09-28, DELTA #15.
        "pattern": r"truncated\s+Bisognano--Wichmann\s+boost|Bisognano--Wichmann\s+boost\s+generator\s+\$?K\s*=|truncated\s+modular\s+flow\s+is|the\s+truncated\s+boost\s+is",
        "exempt_if_nearby": withdrawal_marker("posited-k-as-truncated-bw-boost"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #15; PI-approved remediation).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 family sweeps; PM-verified gate 0 live corpus-wide)",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #16 family sweeps; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "p50-cft3-entangling-s2",
        # For a CFT_3 (Lorentzian R x S^2, Euclidean S^3) a hemisphere's entangling surface is a CIRCLE
        # (S^1); the Euclidean rotation fixes a great circle. 'S^2_equator' is the CFT_4 picture, and the
        # 'dim(dW) = 2' reading of Paper 50 headline (B)'s entropy slope falls with it (the 2 is the
        # ~n_max^2 degeneracy of the lowest K shell). Registered 2026-09-28, DELTA #15.
        "pattern": r"S\^\{?2\}?_\{\\mathrm\{equator\}\}",
        "exempt_if_nearby": withdrawal_marker("p50-cft3-entangling-s2"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #15; PI-approved remediation).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_50_cft3_partition_function.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #15 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "energy-follows-from-casimir-via-shell",
        # MIRROR-direction staleness: Paper 7 Stage 3 now PROVES p0 = Z/n from the S^3 Coulomb-kernel
        # eigenvalue 2pi^2/n (tests/test_paper7_fock_kernel.py). Summaries still say E_n 'follows from
        # n^2-1 through the energy shell' (crediting the label Casimir; P7's own Stage 3 denies it) or that
        # the energies are 'imposed' (understating P7). Paper 1's 'imposed map' wording about ITS OWN
        # operators is correct and is not in scope. Registered 2026-09-28, DELTA #14.
        "pattern": r"follows?\s+(instead\s+)?from\s+(the\s+diagonal\s+Casimir\s+)?\$?n\^2\s*-\s*1\$?[^.]{0,60}energy\s+shell|imposed\s+energy-shell\s+map|are\s+an\s+imposed\s+map\s+from\s+\$?n\$?\s+through\s+Fock",
        "exempt_if_nearby": withdrawal_marker("energy-follows-from-casimir-via-shell"),
        "severity": "fail",
        "scope": "trunk group1 group3 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #14; PI-approved remediation).",
        "files": [
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "papers/synthesis/geovac_field_guide.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/Paper_0_Geometric_Packing.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "p56-injection-theorem-live",
        # Paper 56's level-4 map (thm label injection_g4 kept) is a comparison HOMOMORPHISM, not
        # injective (rank-2 M3 column). Citers still say 'the injection theorem ... \hookrightarrow'.
        # Registered 2026-09-28, DELTA #14.
        "pattern": r"injection\s+theorem|cosmic-Galois\s+injection\b",
        "exempt_if_nearby": withdrawal_marker("p56-injection-theorem-live"),
        "severity": "fail",
        "scope": "trunk group1 group3 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #14; PI-approved remediation).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group3_foundations/paper_56_tannakian_substrate.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "dirac-seminorm-degeneracy-motivation",
        # Paper 38 (owner) assigns the 10/14, 26/55 Dirac-commutator kernel counts to the DOUBLED Fock
        # model; on the true CH window the kernel is the scalars (n_max <= 3). Citers still present the
        # counts as a failure of 'the truthful chirality-diagonal CH Dirac' that the translation seminorm
        # resolves. Registered 2026-09-28, DELTA #14.
        "pattern": r"fails\}?\s+for\s+the\s+truthful\s+chirality-diagonal",
        "exempt_if_nearby": withdrawal_marker("dirac-seminorm-degeneracy-motivation"),
        "severity": "fail",
        "scope": "trunk group1 group3 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #14; PI-approved remediation).",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "p50-chm-integer-spectrum",
        # Paper 50 described the continuum Casini-Huerta-Myers hemisphere modular Hamiltonian as having
        # integer spectrum 2pi J and '2pi-periodic' flow, and tagged 'both spectra integer-valued'
        # [SYMBOLIC]. The continuum modular Hamiltonian generates a boost-like flow with continuous spectrum
        # (P45 Q1 leg (a) says the same); only the Euclidean rotation generator has integer spectrum.
        # Introduced by the v5.17.4 remediation. Registered 2026-09-28, DELTA #14.
        "pattern": r"Both\s+spectra\s+integer-valued|in\s+the\s+single-particle\s+basis\s+are\s+\$?2\\pi\s*J",
        "exempt_if_nearby": withdrawal_marker("p50-chm-integer-spectrum"),
        "severity": "fail",
        "scope": "trunk group1 group3 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #14; PI-approved remediation).",
        "files": [
            "papers/group1_operator_algebras/paper_50_cft3_partition_function.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_50_cft3_partition_function.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "l3-per-harmonic-ratio-bound",
        # The per-harmonic L3 bound (N-1)/sqrt(N^2-1) -> 1 is FALSE (Paper 38 owner; L3 now proved by
        # compression, C3 = 1, [INTERNAL THEOREM]). Paper 32 still states the old bound.
        # Registered 2026-09-28, DELTA #14.
        "pattern": r"\(N\s*-\s*1\)\s*/\s*\\sqrt\{N\^2\s*-\s*1\}",
        "exempt_if_nearby": withdrawal_marker("l3-per-harmonic-ratio-bound"),
        "severity": "fail",
        "scope": "trunk group1 group3 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk DELTA #14; PI-approved remediation).",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "p16-nu-universal-n-minus-2",
        # Paper 16's 'universal ground-state angular quantum number nu = N-2' (nu_min = N - lambda_1).
        # Verified false 2026-09-28 by two routes (S_N character generating function; Young-symmetrizer
        # brute force) + O(3) check: the true rule is the Pauli filling nu_min = E(n_up) + E(n_down),
        # degree d holding (d+1)(d+2)/2 functions. N-2 holds only for the lowest spin with N <= 8; Ne
        # singlet is 10; S = N/2 is spatially ANTISYMMETRIC with large nu (not 0); Li 2S (L=0 even) needs 2;
        # large N grows ~ (3N)^{4/3}/4. mu_free = 2(N-2)(N-1), '4 units per pair', nu/(3N-2) -> 1/3 fall with it.
        "pattern": r"universal\s+(result|angular\s+quantum\s+number)\s+\$?\\?nu\s*=\s*N\s*-\s*2|\\boxed\{\\nu\s*=\s*N\s*-\s*2\}|\\nu_\{\\min\}\s*=\s*N\s*-\s*\\lambda_1|Laughlin-like\s+state\s+with\s+no\s+Pauli|exactly~?\s*4\s+units|\(N-2\)/\(3N-2\)\s*\\to\s*1/3|\\nu\s*=\s*135|\\nu\s*=\s*N\s*-\s*2\$?\s+is\s+independent\s+of\s+the\s+spin|supplies\s+the\s+universal\s+angular",
        "exempt_if_nearby": withdrawal_marker("p16-nu-universal-n-minus-2"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-28 (FULL #9 open-items pass; PI-directed).",
        "files": [
            "papers/group4_quantum_computing/paper_16_periodicity.tex",
            "papers/synthesis/group4_quantum_computing_synthesis.tex",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "geovac/atomic_classifier.py",
        ],
        "cited_by": {
            "papers/group4_quantum_computing/paper_16_periodicity.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- OWNER: eq:nu_min = Pauli filling (two routes + O(3)); Lemma 1; mu_free, tables, N=137 (570/279,015), S=N/2, Cr fixed; tests/test_paper16_nu_min.py",
            "papers/synthesis/group4_quantum_computing_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- eq rewritten to nu_min = E(n_up)+E(n_down); classifier nu described as metadata",
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- L647 'S_N-only Pauli-filling count of Paper 16'",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- L1598 spin diagram vs conjugate spatial irrep corrected",
        },
    },
    {
        "id": "pythagorean-inv-pi2-as-m1-output",
        # The Pythagorean HS-orthogonality residual's 1/pi^2 was read as an M1 (Mellin k=0, Hopf-base
        # measure) output 'predicted' by the case-exhaustion theorem and as 'the same M1 signature' as P38's
        # 4/pi. It enters only through the posited beta = 2pi/kappa_g (H_local = K/beta): the continuum KMS
        # period imported as a normalisation and squared. The orthogonality is a chirality-parity identity
        # independent of beta. Retired 2026-09-28 (four-witness demotion follow-on).
        "pattern": r"exactly\s+the\s+M1\s+output|same\s+\$?1/\\pi\^2\$?\s+factor\s+reappears|predicts\s+the\s+exact\s+M1\s+ring|two\s+appearances\s+of\s+the\s+same\s+\$?k\s*=\s*0\$?\s+volume|Pythagorean[^.]{0,120}M1\s+signature|\$?1/\\pi\^2\$?\s+master-Mellin-engine|lives\s+entirely\s+in\s+the\s+master\s+Mellin\s+engine\s+M1|class\s+of\s+the\s+\$?k\s*=\s*0\$?\s+volume\s+slot",
        "exempt_if_nearby": withdrawal_marker("pythagorean-inv-pi2-as-m1-output"),
        "severity": "fail",
        "scope": "trunk group1 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-09-28 (FULL #9 open-items pass; PI-directed).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- rem:pythagorean_m1_closure: 1/pi^2 = posited beta squared, not M1; L2215 witness list",
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- L314-324 link withdrawn; 'therefore not isolated' -> 'a reading'",
            "papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- tagging sentence",
            "papers/group1_operator_algebras/paper_43_lorentzian_extension.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- 10.2 intro, 'Transcendental tag' para, N_t>1 para, O4",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- two Born-rule loci",
            "papers/group5_qed_gauge/paper_51_gravity_arc.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- false '4/pi sets the 1/pi^2 prefactor' removed",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- witness list",
            "papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- witness-table row daggered class-only",
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- abstract, 8.x, Mellin section, T5, G3(iv), concluding remark",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 open-items pass; PM-verified gate 0 live) -- witness mention",
            "papers/group3_foundations/paper_24_bargmann_segal.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/group3_foundations/paper_57_forced_free_seam.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "finite-four-witness-literal-identification",
        # The finite-cutoff 'four-witness literal identification' (Paper 42 owner; P32/P34/P31/P38/
        # P43/P50/P52 citers). geovac/modular_hamiltonian.py sets K = diag(2 m_j): an integer-spectrum
        # generator, so sigma_{2pi}(O) = O 'exactly for any O' (its own docstring), and the four
        # witnesses are the same K with kappa_g relabelled. A trivial SYMBOLIC fact; it identifies no
        # Bisognano-Wichmann modular flow. The real falsifier (K from a wedge-restricted CH state) was
        # never run. Retired 2026-09-28, /qa trunk FULL (PI-approved).
        "pattern": r"strongest\s+possible\s+level|literal\s+identification\s+(at|of)\s+(the\s+)?(operator[- ]system|finite\s+cutoff|four-witness)|Wick-rotation\s+literal\s+identification|to\s+``literal\s+identification|lift\S*\s+(\S+\s+){0,8}to\s+``?literal\s+identification|four-witness\s+modular-Hamiltonian\s+identification|close\s+the\s+four-witness|canonical\s+strong\s+identification|also\s+gives\s+strong\s+identification|PRIMARY\s+FALSIFIER",
        "exempt_if_nearby": withdrawal_marker("finite-four-witness-literal-identification"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk FULL; PI-approved remediation).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "tests/*.py",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- L1 / L1-tighten paras, sec:modular_hamiltonian_l1, L2-E, BH remark, honest scope, Q/conclusion; bibitem retitled",
            "papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- OWNER: retitled ('the four-witness 2pi period at finite cutoff'); thms tagged [SYMBOLIC]; O0 falsifier OPEN; Tomita circularity stated",
            "papers/group1_operator_algebras/paper_43_lorentzian_extension.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- abstract, intro, six-witness corollary by construction, new sec:scope_automatic, bibitem",
            "papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- L1747-1754 restated as trivial integer-spectrum fact",
            "papers/group1_operator_algebras/paper_50_cft3_partition_function.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- abstract (C), contribution (iii), item list, sec:chm -> Observation (CHM analogy; identification OPEN; GH-rate-as-operator-rate withdrawn)",
            "papers/group1_operator_algebras/paper_52_category_iii_correspondence.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- sec 6.1-6.2 rewritten; boundary-modular content scoped to continuum",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- wick-map-free bullet; Riemannian closure 'not established'",
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- III.27 L1/L1-tighten, table rows, V.E, L2-E, bibitem",
            "papers/group6_precision_observations/paper_35_time_as_projection.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- status note: falsifier still open",
            "papers/synthesis/group1_operator_algebras_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- sec 8.1 retitled + rewritten, scope paragraph, bibitem",
            "papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- header, L188-197/328-331/1249/1487, bibitem",
            "papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- L211/450, Remark 1522, 8.1, Q1 generator labelled posited, bibitem",
            "papers/synthesis/geovac_field_guide.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "order-zero-one-decaying-artifact",
        # Paper 32 called the order-zero/one failures at finite n_max a decaying 'finite-resolution
        # artifact' ('residuals should decay'). Measured: max residual 0.055 -> 0.079 -> 0.096 and failing
        # fraction 63% -> 86% -> 92% (n_max 2,3,4); fixed degree-1 generators approach a nonzero limit.
        # Retired 2026-09-28, /qa trunk FULL.
        "pattern": r"residuals\s+should\s+decay|finite-resolution\s+artifacts?\s+of\s+the\s+same\s+class|Fails\s+at\s+\$?n_\{?\\max\}?\s*\\ge\s*2\$?\s*\(artifact\)|structural-artifact\s+reading\s+of\s+Connes|finite-resolution\s+artifacts?\s+of\s+the\s+truncated|expected\s+to\s+vanish\s+in\s+the\s+GH\s+limit",
        "exempt_if_nearby": withdrawal_marker("order-zero-one-decaying-artifact"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk FULL; PI-approved remediation).",
        "files": [
            "papers/group1_operator_algebras/paper_39_tensor_propinquity_convergence.tex",
            "papers/group1_operator_algebras/paper_43_lorentzian_extension.tex",
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/group1_operator_algebras_synthesis.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- audit item 5 + Reading + Net status + tables: residuals grow (0.055/0.079/0.096), [MEASURED]",
            "papers/group1_operator_algebras/paper_43_lorentzian_extension.tex": "ARCHIVED 2026-09-28 (PI decision): moved to papers/archive/; surviving content is Paper 42 sec. 11 (rem:order_zero_one corrected there)",
        },
    },
    {
        "id": "e6-no-single-mechanism-theorem",
        # Paper 32 thm:no_single_mechanism_K. Fact 2 files the rationals B=42 (M2) and Delta=1/40 (M1) as
        # if the sub-rings were exclusive; Q lies in both, so B+F-Delta lies wholly in the M2 ring and
        # 'mixes M1 and M2' is false. Fact 1 invokes the [CONDITIONAL] case-exhaustion theorem unscoped.
        # Only the empirical Fact 3 survives -> Observation. Adjacent to the S13.5 K prohibition.
        # Retired 2026-09-28, /qa trunk FULL (PI-approved).
        "pattern": r"\\begin\{theorem\}\[No\s+single-mechanism|mixes\s+\$?M_1\$?\s+and\s+\$?M_2|mixes\s+outputs\s+of\s+\$?M_2|eight\}?\s+theorem-grade|thm:no\\?_single\\?_mechanism\\?_K",
        "exempt_if_nearby": withdrawal_marker("e6-no-single-mechanism-theorem"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk FULL; PI-approved remediation).",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- thm -> obs:no_single_mechanism_K + rem:e6_scope",
            "papers/group3_foundations/paper_57_forced_free_seam.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
            "papers/synthesis/geovac_field_guide.tex": "remediated + summary surfaces reread 2026-09-28 (DELTA #14 remediation; PM-verified gate 0 live corpus-wide)",
        },
    },
    {
        "id": "diagonal-commutator-lipschitz",
        # Paper 32 says 'the commutator norm equals the Lipschitz constant on the graph metric'. D_GV
        # and pi(f) are both diagonal in |n_D,kappa,m_j>, so [D_GV, pi(f)] = 0 identically: the diagonal
        # triple carries no metric. Retired 2026-09-28, /qa trunk FULL.
        "pattern": r"commutator\s+norm\s+equals\s+the\s+Lipschitz\s+constant",
        "exempt_if_nearby": withdrawal_marker("diagonal-commutator-lipschitz"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk FULL; PI-approved remediation).",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- bounded-commutators para: [D_GV,pi(f)]=0, spectral not metric",
        },
    },
    {
        "id": "truncation-inherits-odd",
        # Paper 32 prop:odd: 'The truncation inherits the same statement' (no grading at finite n_max).
        # The same paper measures one: gamma_GV = sigma_x (x) 1 anticommutes with D and commutes with
        # [2026-09-28 scope note: this holds for the chirality-doubled MODEL; the TRUE CH window that
        # converges is odd like the continuum (commutant of O = scalars, n_max <= 3, two routes). The
        # retired wording is retired because P32 applied it to its model triple, which is even.]
        # pi(A) (chirality doubling), so the finite triple is EVEN. Odd holds for the continuum factor only.
        # Retired 2026-09-28, /qa trunk FULL.
        "pattern": r"truncation\s+inherits\s+the\s+same\s+statement|Z\}?_2\$?-grading\s*&\s*Absent\s+\(odd\)\s*&\s*Absent",
        "exempt_if_nearby": withdrawal_marker("truncation-inherits-odd"),
        "severity": "fail",
        "scope": "trunk group1 group3 group6 synthesis",
        "note": "Registered 2026-09-28 (/qa trunk FULL; PI-approved remediation).",
        "files": [
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
            "papers/synthesis/*.tex",
        ],
        "cited_by": {
            "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- prop:odd split (a) continuum odd / (b) finite even; tables + captions",
        },
    },
    {
        "id": "pair-diagonal-density-realized",
        # Paper 31 says the composed pipeline REALIZES the pair-diagonal density D_pd = 1.44% and
        # 'at least 98.56% sparse'. Retired 2026-08-30: the realized figure is the global-M_L density
        # D = 6.06% at l_max=3 (93.94% sparse); D_pd survives only as the stricter restriction.
        "pattern": r"realizes\s+the\s+\S{0,12}pair-diagonal|is\s+the\s+figure\s+the\s+composed\s+qubit\s+pipeline\s+realizes|pair-diagonal\s+realized|at\s+least\s+\$?98\.56",
        "exempt_if_nearby": withdrawal_marker("pair-diagonal-density-realized"),
        "severity": "fail",
        "scope": "trunk group2 group3 group5 synthesis",
        "note": "Registered 2026-09-28 (final trunk re-check; PI-approved list-closing pass).",
        "files": [
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
        "cited_by": {
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex": "remediated + summary surfaces reread 2026-09-28 (list-closing agent; PM-verified gate 0 live) -- L37, eq:eri_density rebuilt on global-M_L D, L827 93.94%",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (list-closing agent; PM-verified gate 0 live) -- reviewed: tab:eri_density already global D; no change",
        },
    },
    {
        "id": "alpha-cubic-stated-as-law",
        # Paper 2 states the K combination rule as a LAW ('the physical alpha satisfies',
        # 'the fine structure constant satisfies the cubic'), the p-value as 'establishing ... not
        # accidental', and the Delta sign as 'the (-1)^F grading, not an ad hoc sign choice'.
        # PI 2026-09-28: approximately satisfies (Observation); p-value restated; sign 'consistent with'.
        # Residual 8.8e-8 vs CODATA ~1.5e-10: excluded as an exact law at ~600 sigma. S13.5-adjacent.
        "pattern": r"physical\s+\$?\\?alpha\$?\s+satisfies|fine\s+structure\s+constant\s+satisfies\s+the\s+cubic|unlikely\s+to\s+be\s+accidental|distinguishes\s+the\s+formula\s+from\s+numerology|not\s+an\s+ad\s+hoc\s+sign\s+choice|sign\s+on\s+\$?\\?Delta\$?[^.]{0,60}is\s+the\s+standard|\\?alpha\$?\s+satisfies\s+the\s+cubic|\\alpha\$?\s+is\s+a\s+root\s+of\s+the\s+cubic|combination\s+rule\s+\$?K\s*=\s*\\alpha\^\{-1\}|\\alpha\$\s+is\s+(algebraic\s+over|a\s+root\s+of\s+a\s+polynomial)",
        "exempt_if_nearby": withdrawal_marker("alpha-cubic-stated-as-law"),
        "severity": "fail",
        "scope": "trunk group2 group3 group5 synthesis",
        "note": "Registered 2026-09-28 (final trunk re-check; PI-approved list-closing pass).",
        "files": [
            "papers/group3_foundations/paper_1_spectrum.tex",
            "papers/group3_foundations/paper_24_bargmann_segal.tex",
            "papers/group5_qed_gauge/paper_2_alpha.tex",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/synthesis/group5_qed_gauge_synthesis.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
        ],
        "cited_by": {
            "papers/group5_qed_gauge/paper_2_alpha.tex": "remediated + summary surfaces reread 2026-09-29 (DELTA #17 sweeps A/B/C; PM-verified gate 0 live)",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "remediated + summary surfaces reread 2026-09-28 (list-closing agent; PM-verified gate 0 live) -- L2690 approximately; L3701 sign reading + conditional",
            "papers/synthesis/group5_qed_gauge_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (list-closing agent; PM-verified gate 0 live) -- L57, L320, L338",
            "papers/synthesis/group3_foundations_synthesis.tex": "remediated + summary surfaces reread 2026-09-28 (list-closing agent; PM-verified gate 0 live) -- L1411 plain sum 4.8e-7 vs root 8.8e-8",
            "papers/group3_foundations/paper_24_bargmann_segal.tex": "remediated + summary surfaces reread 2026-09-28 (list-closing agent; PM-verified gate 0 live) -- L1258 approximately / Observation",
            "papers/group3_foundations/paper_1_spectrum.tex": "remediated + summary surfaces reread 2026-09-28 (FULL #9 remediation; PM-verified gate 0 live) -- L505 approximately satisfies (8.8e-8; Observation)",
        },
    },
    {
        "id": "p13-graph-lattice-he-row",
        # Paper 13 Table II carried 'S^3 lattice (n_max=5, Z_eff) -2.851, 1.82%, dim 50' -- the v0.4.x
        # H = k(D-A) + W graph-kinetic model (the retired (1/16)-hop family) -- beside a Full-CI row
        # '-2.854, 1.71%, 215,820' whose dimension is LITHIUM's n_max=5 3e count. Physical fixed-k
        # CI at n_max=5 is 1.93% (dim 55 orbitals).
        "pattern": r"-2\.851\s*&\s*1\.82|-2\.854\s*&\s*1\.71\s*&\s*\$?215",
        "exempt_if_nearby": withdrawal_marker("p13-graph-lattice-he-row"),
        "severity": "fail",
        "scope": "trunk group2 group3 group5 synthesis",
        "note": "Registered 2026-09-28 (final trunk re-check; PI-approved list-closing pass).",
        "files": [
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "remediated + summary surfaces reread 2026-09-28 (list-closing agent; PM-verified gate 0 live) -- Table II -> physical fixed-k CI 1.93% (dim 266); Li-dimension row deleted; L146-148, L566-567",
        },
    },
    {
        "id": "benchmark-graph-equivalence-verdict",
        # benchmark_vs_pyscf.py (and its saved outputs) printed '18/18 symbolic proofs that graph
        # Laplacian is conformally equivalent' and 'Single-electron systems converge to <0.1%' as an
        # HONEST VERDICT -- the retracted graph->S^3 claim and a false win (own table: 0.574%).
        "pattern": r"graph\s+Laplacian\s+is\s+conformally\s+equivalent|Single-electron\s+systems\s+converge\s+to\s+<\s?0\.1|Converges\s+to\s+exact\s+H,\s+He\+",
        "exempt_if_nearby": withdrawal_marker("benchmark-graph-equivalence-verdict"),
        "severity": "fail",
        "scope": "trunk group2 group3 group5 synthesis",
        "note": "Registered 2026-09-28 (final trunk re-check; PI-approved list-closing pass).",
        "files": [
            "benchmarks/scripts/benchmark_vs_pyscf.py",
            "benchmarks/figures/pyscf_comparison.txt",
            "benchmarks/figures/pyscf_comparison_ci.txt",
        ],
        "cited_by": {
            "benchmarks/scripts/benchmark_vs_pyscf.py": "remediated + summary surfaces reread 2026-09-28 (list-closing agent; PM-verified gate 0 live) -- SUMMARY block rewritten; Sat.deficit% column",
            "benchmarks/figures/pyscf_comparison.txt": "remediated + summary surfaces reread 2026-09-28 (list-closing agent; PM-verified gate 0 live) -- text corrected in place, numeric rows unchanged",
            "benchmarks/figures/pyscf_comparison_ci.txt": "remediated + summary surfaces reread 2026-09-28 (list-closing agent; PM-verified gate 0 live) -- text corrected in place, numeric rows unchanged",
        },
    },
    {
        "id": "p60-floor-is-angular",
        "note": "Paper 60's accuracy floor was attributed to l_max truncation "
                "and declared unreachable at any K.  Both halves are refuted "
                "(2026-09-08).  The floor is the SCALE LOCK: metric-free holds "
                "iff E = -lambda^2/2, hence lambda = p_kappa, and that is not "
                "the variational optimum -- freeing lambda over the IDENTICAL "
                "span reaches 1.28 mHa at K=130 against 7.46 locked.  And it is "
                "GROUND-STATE specific: at K=452, 4.28x chemical accuracy for the "
                "ground state vs 1.08x for 2^1S at identical ||M||_1 (at K=202 the "
                "same pair reads 4.49x / 1.12x).  This note named NO K until "
                "2026-09-11, which is a violation of this paper's own branch "
                "criterion -- every quantity names its evaluation domain -- and it "
                "made two correct values look like a contradiction of two others. "
                "Both K=452 figures are DERIVED in the numeric registry "
                "(p60_gnd_ratio_k452 / p60_exc_ratio_k452) so they cannot drift "
                "from their energies.  The genuine He "
                "l>=4 partial-wave tail is 0.187 mHa (ground) and 0.008 mHa "
                "(2^1S) -- 3% and 0.5% of the respective floors -- measured from "
                "the l-increments themselves, so angular truncation cannot be "
                "the mechanism.  "
                "The floor VALUE (6.44 mHa) is unretired and still cited; what "
                "is retired is its attribution and the 'at any K' universality.",
        "pattern": r"cannot reach chemical accuracy at any"
                   r"|angular detail that would close the gap"
                   r"|known and unimprovable"
                   r"|genuine basis incompleteness"
                   r"|(?:cheap growth|Cost growth)[^.\n]{0,80}same fact seen twice"
                   r"|saturat\w*[^.\n]{0,60}at any basis size",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        # Documents whose ARGUMENT rests on this claim (distinct from `files`,
        # which is only where its wording might appear).
        "cited_by": {
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-08",
            "docs/qa/paper_60.done.md": "reviewed 2026-09-08",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/qa/paper_60.done.md",
        ],
    },
    {
        "id": "p12-grid-floor-955",
        "note": "2026-09-14 (round-3 DELTA, code review).  Paper 12 said the "
                "weakest variational value 'anywhere on the full (alpha, "
                "threshold) grid is 95.5%' and that the conclusion is 'robust "
                "across the grid'.  REFUTED by re-measurement: 95.5% is the "
                "cell at alpha = 1.10, threshold 1e-8, and is the minimum only "
                "if alpha is restricted to >= 1.10 -- excluding the low-alpha "
                "half the same paragraph calls pathological.  Over the declared "
                "range alpha in [0.90, 1.30] the weakest variational value is "
                "76.11% (alpha = 0.90, 1e-8), where the sigma-only solve gives "
                "91.06% -- so the azimuthal channels LOSE fifteen points at "
                "that cell.  A restricted-evaluation artifact: a clean floor "
                "produced by deleting the part of the object that breaks it.  "
                "Recorded also as process: this was the fourth consecutive "
                "round whose largest defect sat in the previous round's "
                "remediation, and the sentence was written the same day by the "
                "PM to fix the round-2 defect next to it.",
        "pattern": r"weakest variational value\s+anywhere on the full"
                   r"|ranges over\s+\$?95\.5\$?[^0-9]{0,12}99\.1"
                   r"|robust across the grid",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "reviewed 2026-09-14 -- owner, corrected in place",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14 -- quotes the 99.0-99.1 envelope, never the grid floor",
            "docs/claim_test_matrix.md": "reviewed 2026-09-14 -- no row asserts a grid floor",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
        ],
    },
    {
        "id": "p15-sigma-pi-decoupled",
        "note": "2026-09-14 (round-3 DELTA).  Paper 15 stated that the sigma "
                "(m=0) and pi (|m|=1) sectors are 'completely decoupled in the "
                "angular eigenvalue problem', on the grounds that the nuclear "
                "coupling is diagonal in m and the e-e coupling conserves "
                "M = m1 + m2.  WITHDRAWN: the reason is a non-sequitur -- "
                "(0,0) and (+1,-1) both carry M = 0, so conserving M is "
                "precisely what puts them in the SAME block, and the e-e "
                "multipole expansion couples them (Paper 12 says so explicitly: "
                "its m != 0 terms 'are what couple different mu').  The "
                "measurement agrees: decoupled blocks would make the ground "
                "state a minimum over blocks, so adding pi channels could not "
                "lower it by a near-constant amount, which is what is observed. "
                "The empirical offset survives as an observation; the "
                "decoupling explanation does not.",
        "pattern": r"completely decoupled in the angular"
                   r"|decoupling in the angular eigenvalue problem",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-14 -- owner, withdrawn at both loci",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "reviewed 2026-09-14 -- its m != 0 coupling result is the contradicting evidence",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14 -- does not restate the decoupling claim",
            "docs/claim_test_matrix.md": "reviewed 2026-09-14",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/synthesis/geovac_field_guide.tex",
            "docs/claim_test_matrix.md",
        ],
    },
    {
        "id": "p13-hyperspherical-cusp-advantage",
        "note": "2026-09-14 (round-3 DELTA).  Paper 13 claimed at three loci "
                "that placing the e-e cusp on a coalescence manifold is a "
                "'crucial structural advantage' of hyperspherical coordinates "
                "and 'the key advantage over single-electron coordinate "
                "systems'.  WITHDRAWN as an ADVANTAGE claim: it is a structural "
                "DIFFERENCE.  Paper 15 records that it does not translate into "
                "an accuracy advantage, and Paper 12 reaches 99.1% of H2's D_e "
                "in prolate spheroidal coordinates, where the same locus is not "
                "a coordinate surface at all.  Papers 12 and 15 were corrected "
                "in round 2; Paper 13 was in the same review set and was not "
                "reached until round 3 -- the owner-corrected / citer-stale "
                "class again.",
        "pattern": r"crucial structural advantage of these coordinates"
                   r"|key advantage over single[- ]electron\s+coordinate systems"
                   r"|cusp at a boundary\s+condition rather than a coordinate singularity",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14 -- owner, all three loci corrected",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-14 -- carries the corrected framing already",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "reviewed 2026-09-14",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/synthesis/geovac_field_guide.tex",
        ],
    },
    {
        "id": "p13-sparsity-causes-005",
        "note": "2026-09-14 (round-3 DELTA).  Paper 13 said the SO(6) angular "
                "sparsity 'is ultimately responsible for the 0.05% accuracy' of "
                "the single-channel adiabatic approximation.  WITHDRAWN: the "
                "same paper's closing section attributes that number to "
                "fortuitous cancellation between the adiabatic approximation "
                "error and the finite-difference discretization error.  The two "
                "readings are incompatible -- one makes the number a structural "
                "property of the framework, the other an artifact -- and the "
                "cancellation reading is the one the evidence supports.",
        "pattern": r"ultimately responsible for the 0\.05",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14 -- owner",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14 -- states the non-variational caveat correctly",
            "docs/validation_benchmarks.md": "reviewed 2026-09-14",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/validation_benchmarks.md",
        ],
    },
    {
        "id": "p13-graph-reproduces-hydrogenic",
        "note": "2026-09-14 (round-3 DELTA).  Paper 13 said 'The discrete graph "
                "Laplacian reproduces hydrogenic energies to < 0.1%'.  "
                "WITHDRAWN: this states as a measured accuracy what is a "
                "property of the construction.  E_0 = kappa * lambda_max holds "
                "BY CONSTRUCTION (CLAUDE.md Sec. 5, Level 1), so any such "
                "figure is a bound on the spectral truncation deficit, not an "
                "accuracy against experiment.  Same class as the 'reproduces "
                "the Coulomb degeneracy' paraphrase retired from Paper 1 in "
                "DELTA #10 (v5.10.1).",
        "pattern": r"reproduces\s+hydrogenic energies to",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis trunk",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14 -- owner",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
        ],
    },
    {
        "id": "p18-mu-rho-r-algebraic",
        "note": "2026-09-14 (round-3 DELTA).  Paper 18's sec:algebraic_curve "
                "bullet said the adiabatic eigenvalues mu(R) AND mu(rho,R) are "
                "algebraic functions whose R-dependence is algebraic at every "
                "truncation.  SCOPED: the construction given is a LINEAR MATRIX "
                "PENCIL H_0 + R*V_C, which is what supplies the single global "
                "P(R,mu) = 0.  The Level-4 sweep is not of that form, and the "
                "same section's own later text says so.  What the Level-4 "
                "eigenvalues are instead is left open rather than reclassified "
                "-- the round-1 and round-2 defects in this subsection were "
                "both replacement readings written in place of a withdrawal.",
        "pattern": r"mu\(\\?rho\s*,\s*R\)[^.]{0,70}are\s+algebraic",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group3",
        "cited_by": {
            "papers/group3_foundations/paper_18_exchange_constants.tex": "reviewed 2026-09-14 -- owner",
            "docs/algebraic_registry.md": "reviewed 2026-09-14",
        },
        "files": [
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "docs/algebraic_registry.md",
        ],
    },
    {
        "id": "p18-mu-needed-for-sub-01",
        "note": "2026-09-14 (round-3 DELTA).  Paper 18 said 'the mu(R) "
                "parameterization is needed to achieve sub-0.1% accuracy'.  "
                "WITHDRAWN -- false in both directions: the adiabatic route "
                "that CARRIES mu(R) floors at 0.19-0.20% and never reaches "
                "sub-0.1%, while sub-0.1% IS reached without it by the 2D "
                "variational solver (0.022% raw at l_max = 7), which treats R "
                "and alpha simultaneously.  The Class-S to Class-C reading in "
                "the following sentence does not depend on the necessity claim "
                "and stands.",
        "pattern": r"parameterization is needed\s+to achieve sub-0\.1",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group3",
        "cited_by": {
            "papers/group3_foundations/paper_18_exchange_constants.tex": "reviewed 2026-09-14 -- owner",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14 -- source of both floor figures",
        },
        "files": [
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
        ],
    },
    {
        "id": "p12-envelope-insensitive",
        "note": "2026-09-14 (round-3 DELTA).  Paper 12's stability-envelope "
                "paragraph -- itself written as the round-2 remediation -- said "
                "the conclusion that the azimuthal channels close essentially "
                "all of the 7.6% gap 'is insensitive to the choice', and then "
                "quoted the number that refutes it: the weakest variational "
                "value on the full (alpha, threshold) grid is 95.5%, which "
                "closes about 41% of the gap, not essentially all.  WITHDRAWN. "
                "The qualitative effect is robust; the quantitative value ranges "
                "95.5-99.1% across the grid.  Recorded because this is the third "
                "consecutive round in which the largest defect was in the "
                "previous round's own remediation.",
        "pattern": r"is insensitive to the choice",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": "reviewed 2026-09-14 -- owner",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14 -- carries the envelope, not the insensitivity claim",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
        ],
    },
    {
        "id": "p15-delta-cbs-above-97",
        "note": "2026-09-14 (round-3 DELTA).  Paper 15's remaining-gap budget "
                "said 'CBS extrapolation including delta would push the limit "
                "above 97%'.  WITHDRAWN: the +0.65 pp delta gain it rests on "
                "was measured against SIGMA-ONLY at l_max = 4, not against the "
                "sigma+pi baseline the budget uses, so the transfer is invalid. "
                "The row it comes from is itself unreconciled -- its N_ch = 37 "
                "implies sigma+pi+delta, which cannot return 87.6% when its own "
                "sigma+pi subset returns 93.6%.",
        "pattern": r"push the limit above 97",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-14 -- owner",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14 -- does not carry the >97% figure",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
        ],
    },
    {
        "id": "p12-cusp-as-the-h2-gap",
        "note": "2026-09-14.  Paper 12 attributed its 7.6% H2 residual to the "
                "electron-electron cusp requiring non-analytic r12^(1/2) and "
                "r12 ln r12 terms, and read that as identifying hyperspherical "
                "coordinates as the next natural geometry.  Both are "
                "WITHDRAWN.  The gap is the sigma-only restriction: the basis "
                "is phi-independent and spans only m1 = m2 = 0, while a "
                "1Sigma_g+ state constrains only the TOTAL M.  Restoring the "
                "channels in the same basis reaches 99.09% of D_e (+11.6 mHa, "
                "against 0.34 mHa from tripling the sigma basis); confirmed "
                "independently in a Gaussian basis (99.10%) and externally by "
                "Tao-McCurdy-Rescigno, PRA 82, 023423 (2010), who reach 0.05 "
                "mHa in these same coordinates with a polynomial angular "
                "basis and no r12 factors.  The cusp is still real physics "
                "and still Paper 18's embedding-tier transcendental; it is "
                "not what capped this calculation.  Memo: "
                "debug/sprint_tmr_method_memo.md.",
        "pattern": r"cusp requires (?:the )?non-analytic"
                   r"|requires non-analytic terms"
                   r"|no polynomial prolate\s+spheroidal basis can represent"
                   r"|diagnosis identifies the next natural geometry"
                   r"|confirming the cusp[- ]resolution advantage"
                   r"|cusp[- ]resolution advantage of the"
                   r"|gap (?:was )?attributed to the\s+electron[-–]electron cusp",
        "exempt_if_nearby": r"(?!)",
        "severity": "fail",
        "scope": "group2 synthesis",
        # Documents whose ARGUMENT rests on this claim (distinct from `files`,
        # which is only where its wording might appear).
        # AMENDED after the 2026-09-14 DELTA run.  Paper 18 was MISSING from
        # the original list and is the largest dependent of all: its
        # `sec:mu_level4` builds a TAXONOMY classification on this claim.  Two
        # independent reviewers found it; no pattern could, because Paper 18
        # restates the claim in its own words and its own numbers (92.5%,
        # 7.5%).  Re-priced 2026-09-14 under PI direction.
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-14",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14",
            "papers/group3_foundations/paper_18_exchange_constants.tex": "reviewed 2026-09-14 -- Level-4 subsection re-priced",
            "docs/validation_benchmarks.md": "reviewed 2026-09-14",
            "docs/claim_test_matrix.md": "reviewed 2026-09-14",
            "docs/qa/group2.done.md": "reviewed 2026-09-14",
            "README.md": "reviewed 2026-09-14",
            "docs/topic_to_paper_lookup.md": "reviewed 2026-09-14",
            "docs/paper_notes_archive.md": "reviewed 2026-09-14",
            "tests/test_level4_multichannel.py": "reviewed 2026-09-14",
            "viz/public/papers/paper_12_algebraic_vee.html": "reviewed 2026-09-14 -- regenerated",
            "viz/public/papers/paper_13_hyperspherical.html": "reviewed 2026-09-14 -- regenerated",
            "viz/public/papers/paper_15_level4_geometry.html": "reviewed 2026-09-14 -- regenerated",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/validation_benchmarks.md",
            "docs/claim_test_matrix.md",
        ],
    },
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
        "id": "p61-k12-is-lemniscate",
        "scope": "paper_61 group3 group2",
        "severity": "fail",
        "retired": "2026-09-07 (/qa paper_61, claims dimension); corpus-swept "
                   "2026-09-13 (DELTA). K(1/2) = varpi = Gamma(1/4)^2/(4 sqrt pi) "
                   "was labelled 'the lemniscate constant'. FALSE: the classical "
                   "lemniscate constant is sqrt(2)*varpi = Gamma(1/4)^2/(2 sqrt(2 "
                   "pi)), a DIFFERENT number (off by exactly sqrt 2). The SYMBOL "
                   "K(1/2) is the corpus convention; the NAME was wrong. Retired "
                   "claim #6 got no C16 entry (paper-prose-only fix), so it "
                   "survived in the tracked driver debug/routeC_pslq_fit.py:6 "
                   "through two passes -- the locus-by-locus-remediation class. "
                   "This entry makes the conflation gateable.",
        "pattern": r"K\(\s*1\s*/?\s*2\s*\)[,;]?\s*(?:the\s+)?lemniscate\s+const",
        "exempt_if_nearby": r"differ|NOT\s+the\s+lemniscate|corpus\s+convention|name.{0,8}was\s+wrong",
        "cited_by": {
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex":
                "reviewed 2026-09-13 -- owner; L113-115 states varpi = K(1/2) is "
                "the corpus symbol and the classical lemniscate constant is sqrt2 "
                "varpi (carries the sqrt2 distinction)",
            "debug/routeC_pslq_fit.py":
                "reviewed 2026-09-13 -- the survivor; comment corrected to '= "
                "K(1/2); this is NOT the lemniscate constant, which is "
                "sqrt(2)*varpi'",
        },
        "files": [
            "papers/group3_foundations/paper_61_bessel_moment_periods.tex",
            "papers/group2_quantum_chemistry/paper_59_elliptic_bessel_moment.tex",
            "papers/synthesis/group3_foundations_synthesis.tex",
            "docs/qa/paper_61.done.md",
            "debug/routeC_pslq_fit.py",
            "debug/routeC_cosmic_galois_rung3c.py",
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
        # `is\s+SUPERlinear` REMOVED 2026-09-11 after fire-testing:  it fired on
        # the CORRECTED text ("Neither $T'$ grouping is superlinear") at two
        # loci.  An assertion and its denial share the verb, so the verb cannot
        # discriminate;  the retired VALUES can, and do.
        "pattern": r"is the sublinear block"
                   r"|K\^\{?0\.70\}?"
                   r"|K\^\{?1\.05\}?"
                   r"|local slope 0\.906"
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
        "exempt_if_nearby": r"\[retracted \d{4}-\d\d-\d\d: p60-removability-corollary\]",
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
            "CLAUDE.md",
            "geovac/sturmian_*.py",
            "tests/test_paper60_*.py",
            # widened 2026-09-13 (DELTA #2): these six entries shipped with a
            # narrower list than the OLDER p60 entries next to them, which is why
            # the "Proposition D" loci in the walls register and the single-value
            # floor in claims_register/INDEX were invisible to this gate.
            "docs/claims_register.md",
            "docs/code_architecture.md",
            "docs/topic_to_paper_lookup.md",
            "docs/walls/register.md",
            "docs/qa/paper_60.done.md",
            "papers/INDEX.md",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- the paragraph carries the withdrawal inline",
            "docs/claim_test_matrix.md": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- the tagging row records the corollary as WITHDRAWN",
        },
    },
    {
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
                   "near-dependence is ONE DIRECTION.  SCOPE CORRECTION 2026-10-09: "
                   "every number above is the s-SECTOR sine basis; on a full sigma "
                   "basis (l < n) the deficit falls without plateau and the "
                   "near-null directions multiply (see p60-incomplete-in-molecular-"
                   "metric), so the withdrawal holds for the s-sector only.  The "
                   "phrase stays blocked: a general restatement would still need its "
                   "own measured statement, not this one.",
        "pattern": r"overcompleteness is the price"
                   r"|price of one-cent(?:er|re) completeness"
                   r"|completeness of the one-cent(?:er|re) set[^.\n]{0,40}forces",
        "exempt_if_nearby": r"\[retracted \d{4}-\d\d-\d\d: p60-frames-completeness\]",
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
            "CLAUDE.md",
            "geovac/sturmian_*.py",
            "tests/test_paper60_*.py",
            # widened 2026-09-13 (DELTA #2): these six entries shipped with a
            # narrower list than the OLDER p60 entries next to them, which is why
            # the "Proposition D" loci in the walls register and the single-value
            # floor in claims_register/INDEX were invisible to this gate.
            "docs/claims_register.md",
            "docs/code_architecture.md",
            "docs/topic_to_paper_lookup.md",
            "docs/walls/register.md",
            "docs/qa/paper_60.done.md",
            "papers/INDEX.md",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- sec:molecular carries the plateau measurement",
            "docs/claim_test_matrix.md": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- the one-direction row records the withdrawal",
        },
    },
    {
        "id": "p60-incomplete-in-molecular-metric",
        "scope": "paper_60 group2 synthesis",
        "severity": "fail",
        "retired": "2026-10-09.  'The one-centre shared-scale set is far from complete "
                   "in the molecular metric' (Paper 60 sec:molecular, 2026-09-12) was "
                   "measured on the s-SECTOR only: its 0.380/0.696/0.907 Bessel-deficit "
                   "plateau is the l=0 projection of an off-centre function, which an "
                   "s-only span cannot hold.  On a full sigma basis (n<=N, l<n, m=0; "
                   "debug/sturmian_sigma_basis_sw.py, position-space prolate quadrature "
                   "reproducing the s-sector plateau to 3 digits) the deficit falls "
                   "without plateau (0.158 -> 0.037 over N=2..6 at kR=1), cond(S) is "
                   "exponential in the number of shells, gerade sector included, and the "
                   "near-null directions multiply.  Prompted by Avery & Avery 2017 ch. 9 "
                   "(two-centre expansion ill-conditioned by overcompleteness).",
        "pattern": r"far from complete[^.\n]{0,40}molecular metric"
                   r"|not complete in the molecular metric"
                   r"|one-cent(?:er|re) set is not complete",
        "exempt_if_nearby": r"\[retracted \d{4}-\d\d-\d\d: p60-incomplete-in-molecular-metric\]",
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
            "CLAUDE.md",
            "geovac/sturmian_*.py",
            "tests/test_paper60_*.py",
            "docs/claims_register.md",
            "docs/code_architecture.md",
            "docs/topic_to_paper_lookup.md",
            "docs/walls/register.md",
            "docs/qa/paper_60.done.md",
            "papers/INDEX.md",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex": "reviewed 2026-10-09 -- Avery-book pass -- sec:molecular rewritten: plateau attributed to the l=0 projection, sigma-basis paragraph added, abstract/conclusion/resource scoped to the s-sector",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-10-09 -- Avery-book pass -- molecular paragraph scoped to the s-sector + sigma-basis sentence",
            "docs/walls/register.md": "reviewed 2026-10-09 -- Avery-book pass -- dated scope-correction paragraph; row relabelled BREACHED (s-sector) / OPEN (sigma basis), PI-approved",
            "docs/claim_test_matrix.md": "reviewed 2026-10-09 -- Avery-book pass -- one-direction row scoped; sigma-basis row added",
        },
    },
    {
        "id": "p60-sw-matrix-is-v0-ours",
        "scope": "paper_60 group2 synthesis",
        "severity": "fail",
        "retired": "2026-10-09.  Two Paper 60 provenance/accuracy claims, read against "
                   "Avery & Avery 2017 directly: (i) 'the identification of the "
                   "Shibuya-Wulfman matrix AS the V_0-weighted overlap is ours' -- the "
                   "general-V_0 secular equation (eq:general_v0) is the book's eq. "
                   "(7-14)/(12.48) (Goscinski), the molecular identity-form secular "
                   "equation is its 'Molecular spectra' section, and the book writes each "
                   "SW integral as (n/k)<chi'|1/r_a|chi>; (ii) 'the molecular metric IS "
                   "V_0' for V_0 = all nuclei -- for an atom-centred Coulomb-Sturmian "
                   "basis the eq:sw metric is weighted by the KET'S OWN nucleus; the "
                   "full-V_0-weighted overlap's intra block departs from the identity by "
                   "1.16 at kR=2 (geovac/sturmian_sigma_basis.py), against 2e-12 for S.",
        "pattern": r"metric IS V_?0"
                   r"|metric (?:\\emph\{)?is\}? \$(?:\\mathbf\{)?V\}?_0\$"
                   r"|\\emph\{is\} \$(?:\\mathbf\{)?V\}?_0\$"
                   r"|identification of (?:that|the (?:Shibuya.{1,3}Wulfman|SW)) matrix"
                   r"[^.\n]{0,40}weighted overlap",
        "exempt_if_nearby": r"\[retracted \d{4}-\d\d-\d\d: p60-sw-matrix-is-v0-ours\]",
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
            "CLAUDE.md",
            "geovac/sturmian_*.py",
            "tests/test_paper60_*.py",
            "docs/claims_register.md",
            "docs/code_architecture.md",
            "docs/topic_to_paper_lookup.md",
            "docs/walls/register.md",
            "docs/qa/paper_60.done.md",
            "papers/INDEX.md",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex": "reviewed 2026-10-09 -- Avery-book pass -- general-V0 paragraph re-tiered to [ESTABLISHED, from Avery]; atom-centred metric stated as own-centre-weighted; acknowledgments' 'ours' list corrected",
            "docs/claim_test_matrix.md": "reviewed 2026-10-09 -- Avery-book pass -- eq:general_v0 row re-tiered",
            "docs/qa/paper_60.done.md": "reviewed 2026-10-09 -- Avery-book pass -- DoD item 13 provenance reversed",
        },
    },
    {
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
                   r"|conditioning law[^.\n]{0,24}derived here",
        "exempt_if_nearby": r"(?!)",
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
            "CLAUDE.md",
            "geovac/sturmian_*.py",
            "tests/test_paper60_*.py",
            # widened 2026-09-13 (DELTA #2): these six entries shipped with a
            # narrower list than the OLDER p60 entries next to them, which is why
            # the "Proposition D" loci in the walls register and the single-value
            # floor in claims_register/INDEX were invisible to this gate.
            "docs/claims_register.md",
            "docs/code_architecture.md",
            "docs/topic_to_paper_lookup.md",
            "docs/walls/register.md",
            "docs/qa/paper_60.done.md",
            "papers/INDEX.md",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- abstract now carries [PRIOR ART]",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- the molecular paragraph credits KMS",
            "docs/claim_test_matrix.md": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- the KMS row states identification-only",
        },
    },
    {
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
        "pattern": r"Proposition[~\s]?D\b"
                   r"|our Proposition[^.\n]{0,30}block-diagonal congruence",
        "exempt_if_nearby": r"(?!)",
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
            "CLAUDE.md",
            "geovac/sturmian_*.py",
            "tests/test_paper60_*.py",
            # widened 2026-09-13 (DELTA #2): these six entries shipped with a
            # narrower list than the OLDER p60 entries next to them, which is why
            # the "Proposition D" loci in the walls register and the single-value
            # floor in claims_register/INDEX were invisible to this gate.
            "docs/claims_register.md",
            "docs/code_architecture.md",
            "docs/topic_to_paper_lookup.md",
            "docs/walls/register.md",
            "docs/qa/paper_60.done.md",
            "papers/INDEX.md",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- sec:molecular credits Loewdin/Slater-Koster",
            "docs/claim_test_matrix.md": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- the row records the re-attribution",
        },
    },
    {
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
        "pattern": r"equivalently that the Shibuya--Wulfman\s+operator is multiplication"
                   r"|what is ours[^.\n]{0,60}translation phase"
                   r"|the translation identification is ours",
        "exempt_if_nearby": r"(?!)",
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
            "CLAUDE.md",
            "geovac/sturmian_*.py",
            "tests/test_paper60_*.py",
            # widened 2026-09-13 (DELTA #2): these six entries shipped with a
            # narrower list than the OLDER p60 entries next to them, which is why
            # the "Proposition D" loci in the walls register and the single-value
            # floor in claims_register/INDEX were invisible to this gate.
            "docs/claims_register.md",
            "docs/code_architecture.md",
            "docs/topic_to_paper_lookup.md",
            "docs/walls/register.md",
            "docs/qa/paper_60.done.md",
            "papers/INDEX.md",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- the concession names all three counts",
            "docs/claim_test_matrix.md": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- the symbol-not-translation split recorded",
        },
    },
    {
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
        "pattern": r"tau \(DST-I\)\s*\n?\s*algebra --- which is to say Toeplitz minus Hankel"
                   r"|\\tau\$? \(DST-I\)[^.\n]{0,40}structure of this\s*\n?section"
                   r"|our matrices are in the \$?\\?tau",
        "exempt_if_nearby": r"(?!)",
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "docs/claim_test_matrix.md",
            "CLAUDE.md",
            "geovac/sturmian_*.py",
            "tests/test_paper60_*.py",
            # widened 2026-09-13 (DELTA #2): these six entries shipped with a
            # narrower list than the OLDER p60 entries next to them, which is why
            # the "Proposition D" loci in the walls register and the single-value
            # floor in claims_register/INDEX were invisible to this gate.
            "docs/claims_register.md",
            "docs/code_architecture.md",
            "docs/topic_to_paper_lookup.md",
            "docs/walls/register.md",
            "docs/qa/paper_60.done.md",
            "papers/INDEX.md",
        ],
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- the paragraph states the membership failure",
            "docs/claim_test_matrix.md": "reviewed 2026-09-12 -- /qa paper_60 FULL remediation -- the hypothesis-gap clause records it",
        },
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
                   r"|4\s*(?:->|-->|\\to)\s*367[23]"
                   r"|cond\(S\)[^.\n]{0,30}367[23]"
                   r"|L[\u00b22]-divergence"
                   r"|only well-conditioned posing"
                   # added 2026-09-11: the bare adjective applied to the
                   # overlap, and the verb form.  Deliberately anchored on
                   # "overlap"/"metric" so the paper's own DENIAL ("ordinary,
                   # not ill-conditioned") stays silent -- the denial does not
                   # put the noun within reach of the adjective.
                   r"|overlap[^.\n]{0,60}is ill-conditioned"
                   r"|ill-conditioned[^.\n]{0,30}overlap"
                   r"|ill-conditions with basis"
                   r"|illconditioned",
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
                   "slope was then thought to rise outside it.  CORRECTED AGAIN "
                   "2026-09-08 on a converged domain -- see entry "
                   "p60-tprime-superlinear:  the local slope FALLS (to 0.766 by "
                   "K = 514), the window exponent is 0.82, the nuclear diagonal "
                   "T^0 is not a power law at all (asymptotic 1/2 + O(1/log K)), "
                   "and T' is 0.88 -- SUBlinear, not superlinear.  The 2026-09-07 "
                   "values K^0.70 / K^1.05 / rising-slope are THEMSELVES RETIRED "
                   "and are guarded by that entry.  The paper's own Sec. 7 already said the "
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
        # `nuclear diagonal` and `T^0` REMOVED 2026-09-11 (DELTA I3):  they are
        # the vocabulary of the RETIRED mechanism, so a locus asserting "the
        # nuclear diagonal T^0 is the sublinear block" exempted itself.  The
        # surviving vocabulary names only the WINDOW scoping, which is the
        # thing a corrected locus actually has to carry.
        "exempt_if_nearby": r"K\s*=\s*74|window|over the computed|more slowly than "
                            r"the (?:matrix dimension|configuration count)"
                            r"|NOT an asymptote|corrected 2026-09-0[78]",
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
        # Widened 2026-09-11 (DELTA I1).  The prior list was papers+docs only,
        # so the two highest-traffic surfaces in the corpus -- CLAUDE.md, read
        # by every session and every subagent dispatch, and the tracked geovac/
        # docstrings -- were unreachable by this gate.
        "files": [
            "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/INDEX.md",
            "docs/claim_test_matrix.md",
            "docs/claims_register.md",
            "docs/code_architecture.md",
            "CLAUDE.md",
            "tests/test_paper60_sturmian.py",
            "geovac/sturmian_*.py",
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
    {
        "id": "pairdiag-p19-balanced-prose-values",
        "note": "group2 baseline Batch 3 (2026-09-13).  Paper 19's PROSE carried "
                "the retired pair-diagonal balanced-coupled figures while its own "
                "re-measured tables (2026-09-01, exact global-M_L) carried the "
                "current ones -- a direct prose-vs-table contradiction the "
                "trunk-level pairdiag entry did not reach.  Retired -> current: "
                "n_max=3 balanced Pauli 19,959 -> 127,855 / 1-norm 448.9 -> 436.9 "
                "/ QWC 2,298 dropped; Pauli exponent 3.03 -> 3.74; LiH "
                "balanced/composed ratio 2.53x~2.63x -> 3.25x; polyatomic ratio "
                "sequence 2.63/4.77/7.45x -> 3.25/6.35/10.10x; ERI census 130/195 "
                "-> 214/321 (test_cross_block_eri_count).  Also (2026-09-20, "
                "/qa DELTA): Paper 20 sec:basis-completeness prose carried the "
                "same 19,959 (balanced) PLUS the retired composed 7,879 at "
                "n_max=3 (Q84) -> exact-rule 127,855 / 42,535 (Paper 14 "
                "tab:composed_pauli LiH n3 = 42,535; 7,879 = 42,535/5.40, the "
                "paper's own n3 re-pricing factor).  Scope widened group2 -> "
                "group2 group4 so the Paper 20 locus is gated.",
        "pattern": r"19\{,\}959|19,959|7\{,\}879|7,879"
                   r"|(?:4\.77|7\.45)\$?\\times\$?\s*\(\d-block"
                   r"|2\.53\$?\\times\$?\)? is consistent",
        "exempt_if_nearby": r"retired|RETIRED|pair-diagonal|withdrawn"
                            r"|corrected|under the retired|superseded",
        "severity": "fail",
        "scope": "group2 group4",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex":
                "reviewed 2026-09-20 -- owner; its tables already exact-rule",
            "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex":
                "reviewed 2026-09-20 (/qa DELTA) -- L1158 prose 19,959/7,879 -> 127,855/42,535",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex",
            "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex",
        ],
    },
    {
        "id": "tmr-h2-groundstate-figures",
        "note": "2026-09-20 (/qa DELTA, claims + citations).  Paper 12's abstract "
                "(softened v5.14.11) and CLAUDE.md dropped the Tao-McCurdy-Rescigno "
                "'ground-state H2' figures -- -1.17442 Ha at l_max=6, 0.05 mHa, "
                "99.97% of D_e, and an l_max=2 point -1.1732 Ha -- as UNVERIFIABLE: "
                "the cited paper (PhysRevA.82.023423) computes DOUBLE-PHOTOIONIZATION "
                "cross-sections, NOT a bound H2 ground-state energy / D_e (citation "
                "reviewer fetched the APS abstract this pass; three prior "
                "verification attempts also failed).  The body retained them at "
                "three loci (sec:gap, sec:implications, sec:conclusion) and Paper 15 "
                "L725, and the surviving '99.97% of TMR' collided with Paper 12's "
                "OWN new explicit-r12 99.97% / 0.053 mHa (whose energy -1.1744215 is "
                "nearly identical to the TMR '-1.17442').  Corrected: bodies brought "
                "in line with the already-softened abstract -- TMR cited only as "
                "'an accurate two-electron H2 calculation in these coordinates', no "
                "specific number.  The paper's OWN 99.97% / 0.053 mHa (explicit r12, "
                "registry p12_r12_de_pct / p12_r12_err_mha) is legitimate and "
                "untouched; the pattern keys the retired numerals to the TMR cite so "
                "it never fires on the paper's own result.",
        # Discrimination is IN-PATTERN via joined-text proximity (scan_entry
        # matches the newline-joined document, so `.{0,220}` spans lines).
        # Bare 99.97 / 0.05 mHa ALSO describe the paper's OWN explicit-r12
        # result (L1372/L1390), so they fire ONLY within ~220 chars of a TMR
        # marker.  The TMR-unique energies -1.17442 / -1.1732 fire on their own
        # (the paper's own r12 energy is -1.1744215, guarded by the negative
        # lookahead).  NB: `require_nearby` is NOT read by scan_entry -- proximity
        # must live in the pattern itself.
        "pattern": r"-?1\.17442(?![0-9])"
                   r"|-?1\.1732(?![0-9])"
                   r"|(?:tao_mccurdy_rescigno|McCurdy|Rescigno).{0,220}99\.97"
                   r"|99\.97.{0,220}(?:tao_mccurdy_rescigno|McCurdy|Rescigno)"
                   r"|(?:tao_mccurdy_rescigno|McCurdy|Rescigno).{0,220}0\.05(?![0-9])~?mHa"
                   r"|0\.05(?![0-9])~?mHa.{0,220}(?:tao_mccurdy_rescigno|McCurdy|Rescigno)",
        "exempt_if_nearby": withdrawal_marker("tmr-h2-groundstate-figures"),
        "severity": "fail",
        "scope": "paper_12 paper_15 group2",
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex":
                "reviewed 2026-09-20 -- sec:gap/sec:implications/sec:conclusion TMR "
                "numbers removed; own 99.97% kept; sec:implications understatement fixed",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex":
                "reviewed 2026-09-20 -- L725 TMR '99.97%' -> qualitative",
        },
        "files": [
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex",
        ],
    },
    {
        "id": "p45-panel-column-as-joint-rate",
        # Paper 45's panel column tabulates the spatial rate gamma_SU (C_3 = 1), not C3joint * gammajoint (at N_t = 1, gamma_U = pi/2 exceeds gamma_SU at n_max = 4); the Riemannian-limit identity holds on the spatial leg only, and the panel ratio match is a formula identity, not SU(2) dominance. Self-check #10, DELTA #25.
        "pattern": r"closed-form\s+expression\s+\$\\Cthreejoint\\cdot\\gammajoint\$|confirming\s+that\s+the\s+SU\(2\)\s+factor\s+dominates|they\s+quote\s+\$\\gammajoint\$\s+directly|\\Cthreejoint\(\\nmax,\s*1\)\\,\\gammajoint_",
        "exempt_if_nearby": withdrawal_marker("p45-panel-column-as-joint-rate"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "fejer-circle-rate-T-over-N",
        # The circle Fejer first moment is Theta(T log N_t / N_t), not O(T/N_t) (Paper 45, corrected 2026-10-02). Self-check #10, DELTA #25.
        "pattern": r"T\(\\Nt\)/\\Nt|O\(T/\\Nt\)|\\tfrac\{T\}\{\\Nt\}|\+\s*T\s*/\s*\\Nt\s*\\?b?r?\)|rate\s+\$T/\\Nt|\$T/\\Nt\$\s+Fej",
        "exempt_if_nearby": withdrawal_marker("fejer-circle-rate-T-over-N"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "p45-p46-propinquity-credited",
        # Paper 45 proves a degeneracy theorem, not a propinquity convergence; Paper 46 (descoped) preserves no propinquity bound -- its panel coincidence is the same spatial rate formula. Self-check #12, DELTA #25.
        "pattern": r"natural\s+setting\s+for\s+the\s+Paper~?45\s+propinquity\s+proof|leaving\s+the\s+Latr\S*\s+propinquity\s*\n?\s*bound\s+bit-exactly\s+preserved|\$\\Lambda\$\s+is\s+the\s+standard\s+Latr",
        "exempt_if_nearby": withdrawal_marker("p45-p46-propinquity-credited"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "market-test-1norm-097",
        # Composed LiH vs STO-3G (exact rule): near parity on raw Pauli counts (838 vs 907) and 1-norm (1.007x), 3.0x MORE than the qubit-reduced 276; the '1-norm parity (0.97x), 2.7x fewer Pauli' framing is retired. Self-check #12, DELTA #25.
        "pattern": r"1-norm\s+parity\s+with\s+STO-3G\s+\(0\.97|2\.7\s*[x×]\s+fewer\s+Pauli",
        "exempt_if_nearby": withdrawal_marker("market-test-1norm-097"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "gh-rate-credited-to-qubit-library",
        # Paper 38's state-space GH truncation rate is a property of the S^3 spectral triple, not of the molecular qubit Hamiltonians and not an energy-error bound. Self-check #13, DELTA #25.
        "pattern": r"library\s+of\s+qubit\s*\n?\s*Hamiltonians\s+with\s+a\s+proven\s+state-space",
        "exempt_if_nearby": withdrawal_marker("gh-rate-credited-to-qubit-library"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "lih-1norm-34p5",
        # Composed LiH 1-norm is 34.0 Ha including identity (registry lih_composed_lambda alias 34.036; Paper 20), vs STO-3G 34.3; the 34.5 Ha / 1.007x figure is superseded. Self-check #13, DELTA #25.
        "pattern": r"34\.5\s*(?:Ha\s*)?vs\.?\s*34\.3|1\.007\s*[x×]|34\.5\s*Ha\s*\|\s*34\.3",
        "exempt_if_nearby": withdrawal_marker("lih-1norm-34p5"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "p39-closes-propinquity-gap",
        # Paper 39 proves van Suijlekom state-space GH convergence for the S^3 x S^3 tensor triple; the propinquity question stays open. 'The present paper closes that gap' (after the propinquity-literature paragraph), 'Open prior to this work', and crediting Paper 39 with a tensor-product propinquity are retired. Self-check #8, DELTA #25.
        "pattern": r"present\s+paper\s+closes\s+that\s+gap|Open\s+prior\s+to\s+this\s+work|tensor-product\s+propinquity\s+of\s+Paper~?39|Tensor-product\s+propinquity\s+convergence\s+for|Paper\s+39\s+closes\s+the\s+same-manifold|is\s+closed\s+in\s+a\s*\n?\s*companion\s+standalone\s+paper",
        "exempt_if_nearby": withdrawal_marker("p39-closes-propinquity-gap"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "code-layer-maxgamma-as-bound",
        # Self-check #8, DELTA #25: max(gamma_a, gamma_b) (the propinquity_bound field) is a reference value, not a proved bound and not propinquity; the module proves nothing numerically.
        "pattern": r"The\s+joint\s+propinquity\s+bound\s+\(factorized|return\s+the\s+scalar\s+joint\s+propinquity\s+bound|confirms\s+the\s+joint\s*\n?\s*propinquity\s+bound|Verify\s+the\s+joint\s+bound\s+at\s+the\s+largest|absorb\s+into\s+the\s+C_3\^\{\(2\)\}\s*=\s*2|joint\s+propinquity\s+bound\s+qualitative-rate|Joint\s+propinquity\s+rate\s+bound|L5-T\s+assembly:\s+joint\s+propinquity\s+bound",
        "exempt_if_nearby": withdrawal_marker("code-layer-maxgamma-as-bound"),
        "severity": "fail",
        "scope": "trunk group1 group3 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25 self-check #8, code-layer gate; PI-approved).",
        "files": [
            "tests/*.py",
            "geovac/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "p45-as-master-instance",
        # Paper 45's SU(2) x U(1) product carrier is outside Paper 40's semisimple class (a circle factor carries 2/pi) and is not an instance of Paper 39's multi-factor form. Self-check #7, DELTA #25.
        "pattern": r"instantiates\s+the\s+\$k\s*=\s*2\$|U\(1\)\$\s+instance\s+of\s+the\s+multi-factor",
        "exempt_if_nearby": withdrawal_marker("p45-as-master-instance"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "fact-column-as-bound",
        # Paper 39's Lambda^fact = max(gamma_a, gamma_b) is a reference value, not a proved bound (the proved gamma^(ab) >= max(gamma)). Self-check #7, DELTA #25.
        "pattern": r"\\Lambda\^\{\\mathrm\{fact\}\}\s*\\le\s*\\max|joint\s+\$\\max\$-rate\s+is\s+dominated|factorised-panel\}?\s+bound\s+\$\\Lambda",
        "exempt_if_nearby": r"withdrawn|no\s+lemma|" + withdrawal_marker("fact-column-as-bound"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "p40-interior-closure-established",
        # Paper 40's interior-summand (INT) is panel-verified on weight-bounded panels and an analytical sketch at general rank; 'analytic closure ... is established', 'What follows proves (DT)/(INT)', 'carries through verbatim' overclaim it. Self-check #7, DELTA #25.
        "pattern": r"analytic\s+closure\s+of\s+\\eqref\{eq:INT\}\s+on\s+this\s+class\s+is\s+established|follows\s+proves\s+the\s+Dirac-triangle|structure\s+of\s+Steps\s+1--5\s+carries\s+through\s+verbatim",
        "exempt_if_nearby": withdrawal_marker("p40-interior-closure-established"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "code-layer-c3-live-legacy",
        # Self-check #7, DELTA #25: the code layer must not present the triangle / per-irrep / '<= 2' estimates as the live joint constant (C_3^(2) = 1 by compression, Paper 39 Lemma L3-T), nor say the C_3 / Berezin estimates are open for the doubled operator, nor give the (B,P) route's max(gamma) as the theorem's rate.
        "pattern": r"survive\s+it\s+is\s+OPEN|live\s+(?:value|bound)\s+is\s+the\s+triangle|live\s+=\s+triangle|R1-tightened|PROOF-SKETCHED\s+WITH\s+NUMERICAL|assembly\s+with\s+rate\s*\n?\s*max\(gamma",
        "exempt_if_nearby": withdrawal_marker("code-layer-c3-live-legacy"),
        "severity": "fail",
        "scope": "trunk group1 group3 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25 self-check #7, code-layer gate; PI-approved).",
        "files": [
            "tests/*.py",
            "geovac/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "p39-lemmas-open-for-doubled",
        # Paper 39's L3-T sharp bound and L4-T(d) are proved for the Clifford-doubled operator (2026-10-02); 'whether L3-T/L4-T hold for it is open' / 'what those lemmas establish for the repaired operator is open' contradicts that. Self-check #6, DELTA #25.
        "pattern": r"hold\s+for\s+it\s+is\s+open|establish\s+for\s+the\s+repaired\s+operator\s+is\s+open",
        "exempt_if_nearby": withdrawal_marker("p39-lemmas-open-for-doubled"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "kfold-rate-max-form",
        # The k-fold / master bound's rate is the joint moment (<= sum of factor rates), not max_j gamma, which is below the proved k=2 bound gamma^(ab) >= max(gamma_a, gamma_b). Self-check #6, DELTA #25.
        "pattern": r"C_\{?3\}?\^\{\(k\)\}\s*\\cdot\s*\\max_\{?j|joint\s+rate\s+\$?\\gamma\^\{\(k\)\}\s*=\s*\\max|joint\s+rate\s+is\s+the\s+maximum\s+of\s+the\s+factor\s+rates",
        "exempt_if_nearby": r"[Ee]arlier\s+version|withdrawn|" + withdrawal_marker("kfold-rate-max-form"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "interior-int-all-ranks",
        # Paper 40's interior-summand (INT) lemma is panel-verified at A2/A3/C2/G2 and a sketch at general rank (Step 4's cancellation asserted); 'holds on every interior summand' unconditionally, 'signed-sum mechanism succeeds at general rank', 'PRV case closes via the dominance-order fact' (it falls short by |rho|^2) and 'established analytically' overclaim it. Self-check #6, DELTA #25.
        "pattern": r"Hence\s+\\eqref\{eq:INT_interior\}\s+holds\s+on\s+every|mechanism\s+succeeds\s+at\s+general\s+rank|closes\s+via\s+the\s+standard\s+dominance-order|is\s+established\s+analytically\s+in\s+Lemma~\\ref\{lem:L3_interior\}|promotes\s+the\s+per-PRV-orbit",
        "exempt_if_nearby": withdrawal_marker("interior-int-all-ranks"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "triangle-controls-joint-constant",
        # In Paper 39 the triangle bound gives only sqrt(2); compression gives the joint constant 1 (Lemma L3-T), and propinquity is not claimed. 'the triangle bound ... is what controls the joint Lipschitz comparison constant' / '... suffices to make the tensor-product propinquity convergent' are retired. Self-check #6, DELTA #25.
        "pattern": r"triangle\s+bound\s+on\s+that\s+sum\s+is\s+what\s+controls|triangle\s+bound,\s+which\s+suffices\s+to\s+make|controlled\s+by\s+the\s+\(tight\)\s+operator-norm",
        "exempt_if_nearby": withdrawal_marker("triangle-controls-joint-constant"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "c3-2-sqrt2-as-constant",
        # Paper 39's joint Lipschitz comparison constant is C_3^(2) = 1 by compression on the Clifford-doubled module (Lemma L3-T, 2026-10-02); sqrt(2) is only the elementary triangle + Cauchy-Schwarz bound, and was used downstream against the classical Lipschitz norm it does not bound. Self-check #5, DELTA #25.
        "pattern": r"C_3\^\{\(2\)\}\s*\\le\s*\\sqrt\s*\{?2\}?",
        "exempt_if_nearby": r"triangle|[Ee]lementary|superseded|Cauchy|only|withdrawn|record|" + withdrawal_marker("c3-2-sqrt2-as-constant"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "c3-k-sqrtk-growth",
        # The k-fold constant is C_3^(k) = 1 by compression (Paper 39 Sec. 6(iii)); 'C_3^(k) <= sqrt(k)' / 'a genuine sqrt(k) factor' / 'grows with the number of factors' was the triangle route, and the growth claim rested on Paper 38's false per-harmonic bound. Self-check #5, DELTA #25.
        "pattern": r"genuine\s+\$?\\sqrt\s*\{?k|C_\{?3\}?\^\{\(k\)\}\s*(?:\\;)?\s*\\le\s*(?:\\;)?\s*\\sqrt\s*\{?k|C_3\^\{\(k\)\}\s*\\nearrow|constant\s+grows\s+with\s+the\s+number\s+of\s+factors",
        "exempt_if_nearby": r"Cauchy|triangle|withdrawn|only|" + withdrawal_marker("c3-k-sqrtk-growth"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "int-from-prv-theorem-b",
        # Paper 40 Theorem prv_summand_bound(b) as printed falls short of (INT) by |rho|^2, so (INT) on the PRV class rests on the exact-arithmetic panels, and the corollary is a (DT)/(INT) status, not an 'L3 closure at all ranks'. Self-check #5, DELTA #25.
        "pattern": r"holds\s+on\s+this\s+entire\s+family\s+by|L3\}?\s+closure\s+at\s+all\s+ranks\]|Lemma~\\ref\{lem:L3\}\s+closure\s+at\s+all\s+ranks|establishes\s+\(INT\)\s+on\s+the\s+asymptotic\s+family|rigorously\s+covered\s+at\s+all\s+ranks|PRV-summand\s+bound\s+at\s+all\s+ranks",
        "exempt_if_nearby": withdrawal_marker("int-from-prv-theorem-b"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "c3-per-summand-bridge",
        # Paper 40's Dirac-triangle / PRV / Brauer-Klimyk analysis is representation theory; its bridge to the commutator norm (entrywise -> operator norm, eigenvalue ratio read as C_3(G,sigma)) is the rank-1 step Paper 38 refuted. It is not a proof or a check of L3 (proved by compression). Self-check #4, DELTA #25.
        "pattern": r"independent\s+per-summand\s+check|per-summand\s+route\s+also\s+reaches|independently\s+via\s+PRV|independently\s+is\s+asymptotically|analysis\s+is\s+an\s+independent\s+check|C_3\(G,\s*\\sigma\)\s*\\le\s*1\$?\s+asymptotic|asymptotic-tight\s+\$?C_3\s*=\s*1\$?\s+(?:at\s+all\s+ranks\s+)?via",
        "exempt_if_nearby": withdrawal_marker("c3-per-summand-bridge"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "c3-2-closed-form",
        # Paper 39's diagonal closed form C_3^(2)(n,n) = sqrt(2n/(n+2)), its values and the per-cutoff k-fold factor sqrt((N*-1)/(N*+1)) rest on Paper 38's false per-harmonic bound and are withdrawn; C_3^(2) = 1 (compression, Lemma L3-T) stands. Self-check #4, DELTA #25.
        "pattern": r"\\sqrt\{2\s*n\s*/\s*\(n\s*\+\s*2\)\}|\\sqrt\{6/5\}\s*\\approx|\\sqrt\{\\frac\{N_\{\*\}\s*-\s*1\}|\\sqrt\{\(N_\{\*\}\s*-\s*1\)",
        "exempt_if_nearby": r"[Ee]arlier\s+versions?|withdrawn|[Cc]orrected 2026-10-02|draft claimed|" + withdrawal_marker("c3-2-closed-form"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "berezin-mass-weight-form",
        # Paper 38's Berezin map uses the multiplier symbol sigma_nmax(N); the mass-distribution weight N/Z_nmax is not unital (B(1) = Z^{-1} I) and fails the approximate identity. Self-check #3, DELTA #25.
        "pattern": r"\\widehat\{K\}\(N\)\s*=\s*N\s*/\s*Z",
        "exempt_if_nearby": withdrawal_marker("berezin-mass-weight-form"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "heh-931-as-result",
        # HeH+ hyperspherical 93.1% D_e is an adiabatic, non-variational over-binding coincidence (variational ~5%); it may not be quoted as a validation result. Self-check #3, DELTA #25.
        "pattern": r"HeH[^.\n]{0,80}93\.1\s*\\?%",
        "exempt_if_nearby": r"adiabatic|non-variational|over-binding|coincidence|" + withdrawal_marker("heh-931-as-result"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "he-hybrid-h1-example",
        # h1_method='hybrid' adds the retired non-physical (1/16)A graph hop; examples and docs must use the physical default. Self-check #3, DELTA #25.
        "pattern": r"h1_method\s*=\s*['\"]hybrid['\"]",
        "exempt_if_nearby": r"[Ww]ithdrawn|HOP MODEL|non-physical|retired|" + withdrawal_marker("he-hybrid-h1-example"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "p5-accuracy-label",
        # P5 packing-reachability's 98.3% is an AGREEMENT rate (an internal-consistency check, on the 2026-06-08 labelling), not an accuracy of a principle. Self-check #3, DELTA #25.
        "pattern": r"98\.3\s*\\?%[- ]accura|98\.3\s*\\?%\s+accuracy",
        "exempt_if_nearby": withdrawal_marker("p5-accuracy-label"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "ho-closures-as-magic-numbers",
        # 2, 8, 20, 40, 70, 112 are the spin-orbit-free HO shell closures, not the physical nuclear magic numbers (2, 8, 20, 28, 50, 82, 126). Self-check #3, DELTA #25.
        "pattern": r"magic\s+numbers\s*\$?\\?\{?\s*2,\s*8,\s*20,\s*40",
        "exempt_if_nearby": r"not\s+the\s+physical|spin-orbit-free|NOT the physical|" + withdrawal_marker("ho-closures-as-magic-numbers"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "transcendental-classification-label",
        # The master Mellin engine / Paper 18 / Paper 55 TAG the transcendentals of a census (a census and a reading, not a classification). Retires 'transcendental classification', 'the master Mellin engine's classification', 'Paper 55 classifies'. Self-check #3, DELTA #25.
        "pattern": r"transcendental(-tier)?\s+classification|Mellin\s+engine'?s\s+classification|Paper~?\s?55(~\\cite\{[^}]*\})?\s+classifies|Paper~?\s?18(~\\cite\{[^}]*\})?\s+classifies\s+transcendental",
        "exempt_if_nearby": withdrawal_marker("transcendental-classification-label"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (DELTA #25 self-check #3, idea-level gate; PI-approved remediation).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "README.md",
            "papers/INDEX.md",
            "docs/claims_register.md",
        ],
        "cited_by": {},
    },
    {
        "id": "composed-req-range-26pct",
        # The composed R_eq error range is 5.3-19.4% (LiH 5.3, BeH2 11.7, H2O 19.4; Paper 17). The '26%' upper end was H2O's earlier 1.34-bohr value, superseded by solver evolution (Paper 17, H2O PES paragraph). Distinct from the 5-26% Z^2-scaled PK-PARAMETER errors (SCOPE_BOUNDARY, atomic_classifier), which are a different quantity and stay. 2026-10-02, PI-directed.
        # widened 2026-10-02 (self-check #3): per-molecule forms ('26% for H2O', 'H2O carries a 26%')
        "pattern": r"5(\.3)?\s*(--|\u2013|-)\s*26\s*\\?\s*%|26\s*\\?%\s*for\s+H|H\$?_?\{?2\}?\$?O[^.;]{0,30}26\s*\\?%",
        "require_nearby": r"R_?\{?\\?(rm\s*)?eq|R_eq|equilibrium|geometry|bond length|composed",
        "exempt_if_nearby": r"Z\^?2-scaled|Z\u00b2-scaled|PK parameters|earlier|retired|historical|" + withdrawal_marker("composed-req-range-26pct"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (PI-directed after /qa trunk DELTA #25: 'if CLAUDE.md is wrong, change CLAUDE.md').",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "README.md",
            "CLAUDE.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/claims_register.md",
            "papers/INDEX.md",
        ],
        "cited_by": {},
    },
    {
        "id": "depth-equals-loop-order",
        # Loop order BOUNDS motivic depth (realized depths drop, Paper 55); 'depth equal to (the) loop order' and 'loop order = motivic depth' are retired. DELTA #25 self-check.
        "pattern": r"depth\s+equal\s+to\s+(the\s+)?loop\s+order|loop\s+order\s+\$?=\$?\s+motivic\s+depth|depth\s+\$?k\s*=\s*\$?\s*loop\s+order",
        "exempt_if_nearby": withdrawal_marker("depth-equals-loop-order"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "engine-predictive-within-domain",
        # The engine 'predictive within each domain' re-asserts the downgraded k=0 GH-rate domain; the domain partition is a reading. DELTA #25 self-check.
        "pattern": r"engine\s+is\s+predictive|predictive\s+\\emph\{within\}|Mechanism-as-domain\s+sharpening",
        "exempt_if_nearby": withdrawal_marker("engine-predictive-within-domain"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "circulant-explains-cubic",
        # Any depressed cubic with unit constant term is realizable by a traceless 3x3 circulant, so the circulant realizes but does not select or explain Paper 2's cubic. DELTA #25 self-check.
        "pattern": r"circulant\s+explains|explain\s+the\s+cubic\s+form|cubic\s+form\s+is\s+natural|cubic\s+form\s+arises\s+naturally|suggesting\s+democratic\s+coupling",
        "exempt_if_nearby": withdrawal_marker("circulant-explains-cubic"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "y0-empirically-falsified",
        # The strict Y = 0 slice is disfavoured (the empirical fact excluding it is the nonzero measured fermion masses); not 'empirically falsified by the sphaleron transition at T_C ~ 160 GeV' (a lattice-SM output). DELTA #25.
        "pattern": r"empirically\s+falsified\s+by\s+the\s+sphaleron|Empirical\s+falsification\s+of\s+the\s+strict",
        "exempt_if_nearby": withdrawal_marker("y0-empirically-falsified"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "inner-factor-period-ring-forced",
        # eta-trivialization removes the M3 SLOT on the inner factor; it does not force an inner-factor period ring M1 u M2 (the Yukawa values carry their own generators y_i). DELTA #25 self-check.
        "pattern": r"forces\s+(the\s+)?inner-factor\s+period\s+ring|only\s+periods\s+allowed\s+on\s+Connes|inner-factor\s+period\s+ring\s+\$?M_1\s*\\cup\s*M_2",
        "exempt_if_nearby": withdrawal_marker("inner-factor-period-ring-forced"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "m1-governs-gh-rate-d25",
        # The GH-rate constant 4/pi shares its value with the M1 volume slot -- a reading, not a derivation (Paper 38). Retires 'M1 produces / owns / is the natural domain of' the GH rate, 'pure M1 expression', 'lives in M1 alone', 'exactly what the M1 mechanism predicts', '4/pi identifies (the GH rate) with the M1 sub-mechanism'. DELTA #25.
        "pattern": r"natural\s+domain\s+of\s+(state-space\s+)?G|M1\s+owns|mechanism\s+produces\s+the\s+\$?4/|pure\s+M1\s+expression|lives\s+in\s+\$?M_?1\$?\s+alone|exactly\s+what\s+the\s+M1\s+mechanism\s+predicts|identifies\s+the\s+state-space\s+GH\s+rate\s+with|4/\\pi\$?\s+identifies\s+with|M1\s+\$?\\leftrightarrow\$?\s+state-space\s+GH\s+rates",
        "exempt_if_nearby": withdrawal_marker("m1-governs-gh-rate-d25"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "p40-closes-positively",
        # Paper 40's rank>=2 universality is numerically supported and fit-sensitive (rigorous at rank 1); not 'a class-wide theorem' that 'closes it positively', and not 'closed in Paper 40'. DELTA #25.
        "pattern": r"closes\s+it\s+positively|to\s+a\s+class-wide\s+theorem|was\s+closed\s+in\s+Paper~?\s?40",
        "exempt_if_nearby": withdrawal_marker("p40-closes-positively"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "c3-gamma-display-live",
        # The Paper 38/40 bound is d_GH <= gamma with constant 1 (lifted-state route); the C_3 * gamma display belongs to the withdrawn L5 approximation-pair route. DELTA #25.
        "pattern": r"\\le\\?;?\s*C_\{?3\}?(\(G\))?\s*\\cdot\s*\\gamma|\\le\\?;?\s*\\Cthree(?:SU|joint)(?:\([^)]*\))?\s*\\cdot\s*\\gamma",
        "exempt_if_nearby": r"withdrawn|L5|earlier display|abandoned|record|retracted|" + withdrawal_marker("c3-gamma-display-live"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "per-harmonic-ratio-live",
        # Paper 38's per-harmonic Lipschitz ratio sqrt((N-1)/(N+1)) / (N-1)/sqrt(N^2-1) is FALSE (withdrawn 2026-09-28); C_3 = 1 follows by compression. DELTA #25.
        "pattern": r"\(N-1\)/\\sqrt\{N\^\{?2\}?\s*-\s*1\}|\\sqrt\{\(N-1\)/\(N\+1\)\}|\\sqrt\{\(2j-1\)/\(2j\+1\)\}",
        "exempt_if_nearby": r"false|withdrawn|earlier version|" + withdrawal_marker("per-harmonic-ratio-live"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "vinberg-1990-attribution",
        # The PRV-summand attribution to Vinberg 1990 carried no locator and is withdrawn (Paper 38); cite Kumar 1988. DELTA #25.
        "pattern": r"Vinberg",
        "exempt_if_nearby": r"withdrawn|no\s+locator|" + withdrawal_marker("vinberg-1990-attribution"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "p38-sole-cvs-credit",
        # The Connes--van Suijlekom deferred convergence question is answered on the Riemannian side for scalar compact groups by Gaudillot-Estrada--van Suijlekom and for the spinor commutative factor by Paper 38; not 'closed by Paper 38'. DELTA #25.
        "pattern": r"which\s+Paper~38\s+closed\s+on\s+the|closed\s+on\s+the\s+Riemannian\s+side\s+by\s+Paper~38",
        "exempt_if_nearby": withdrawal_marker("p38-sole-cvs-credit"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "yukawa-structural-reason",
        # The Yukawa / inner-factor residue is an Observation-level reading: eta-trivialization restricts the inner factor to the M1/M2 SLOTS; it does not force the period ring, supply 'the structural reason', or yield a 'forced' inner-factor basis. DELTA #25.
        "pattern": r"forces\s+the\s+inner-factor\s+period\s+ring|This\s+is\s+the\s+structural\s+reason\s+F1|supplies\s+the\s+structural\s+reason|forced\s+inner-factor\s+pure-Tate",
        "exempt_if_nearby": withdrawal_marker("yukawa-structural-reason"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "three-homes-in-theorem",
        # That no common generator produces B, F, Delta is an empirical record (twelve eliminated), not theorem content and not an exclusion. DELTA #25 (Paper 18 Three-homes theorem).
        "pattern": r"do\s+not\s+arise\s+from\s+a\s+single\s+common\s+generator|exclude\s+any\s+common|No\s+common\s+generator\s+produces",
        "exempt_if_nearby": withdrawal_marker("three-homes-in-theorem"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "p55-all-transcendentals",
        # Paper 55 is a census of S^3 spectral-zeta/heat-trace data; chemistry seeds, spectral determinants and Paper 61's elliptic period lie outside it. Retires 'all GeoVac transcendentals' framings and 'classification'. DELTA #25.
        "pattern": r"We\s+show\s+that\s+the\s+answer\s+is\s+affirmative|across\s+the\s+entire\s+GeoVac\s+corpus,?\s+only|full\s+classification\s+fits|master\s+Mellin\s+engine\s+classification|classification\s+of\s+the\s+transcendentals\s+appearing|do\s+all\s+of\s+GeoVac'?s\s+transcendentals",
        "exempt_if_nearby": withdrawal_marker("p55-all-transcendentals"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "delta-equals-c2-32",
        # Delta = 1/40 = 1/g_3^Dirac; c^2(3,2) = 1/8 under Paper 2's formula, so 'Delta = c^2(3,2)' is a false equation (Paper 57, DELTA #25).
        "pattern": r"\\Delta\s*=\s*1/40\s*=\s*c\^2\(3,\s*2\)",
        "exempt_if_nearby": withdrawal_marker("delta-equals-c2-32"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "k-conjectural-front",
        # K = pi(B+F-Delta) is an Observation (section 13.5); 'remains conjectural' and 'eight theorem-grade non-selection results' (E6 is an Observation) are retired on the front door. DELTA #25.
        "pattern": r"remains\s+conjectural|retains\s+its\s+\\emph\{conjectural\}|Conjectural\s+status\s+of\s+K|eight\s+theorem-grade\s+non-selection",
        "exempt_if_nearby": r"Historical|historical|\[Historical|" + withdrawal_marker("k-conjectural-front"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "door4-series-theorems",
        # The Door 4 record is partial-door and its seam an Observation-level reading; not 'the Door 4 series theorems'. DELTA #25.
        "pattern": r"Door~?\s?4\s+series\s+theorems",
        "exempt_if_nearby": withdrawal_marker("door4-series-theorems"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "hm-as-convergence",
        # Hekkelman--McDonald construct a noncommutative integral on spectrally truncated triples (tori are examples); it is not a GH/propinquity convergence result and does not 'live on T^d'. DELTA #25 citation review.
        "pattern": r"Their\s+construction\s+lives\s+on\s+\$T\^d\$|Hekkelman--McDonald~\\cite\{hekkelman_mcdonald2024b\}\s+all\s+extend|Hekkelman--McDonald--Toyota",
        "exempt_if_nearby": withdrawal_marker("hm-as-convergence"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #25, idea-level gate; PI-approved).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/geovac_positioning.md",
            "docs/geovac_onepager.md",
            "docs/validation_benchmarks.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "p40-semisimple-universality",
        # Paper 40's 4/pi universality holds (at its tiers) for SIMPLE G only; the semisimple corollary (factor-wise 'c/2 + c/2') is withdrawn -- the joint moment is not additive (Paper 39 measures ~0.74 of the sum). DELTA #26.
        "pattern": r"universal\s+(?:across|over)\s+(?:the\s+)?compact\s+(?:connected\s+)?(?:\*\*)?semisimple|c\(G_1\)/2\s*\+\s*c\(G_2\)/2|[Uu]niversality\s+is\s+a\s+forcing|rank-invariant\s+across\s+(?:the\s+)?(?:compact\s+connected\s+)?semisimple|universality\s+across\s+the\s+semisimple\s+class|universal\s+across\s+the\s+class\s+of\s+compact\s+connected\s+semisimple",
        "exempt_if_nearby": withdrawal_marker("p40-semisimple-universality"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #26, family gate; remediation PI-directed).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
            "docs/topic_to_paper_lookup.md",
        ],
        "cited_by": {},
    },
    {
        "id": "rigorous-at-rank-1",
        # Paper 38's 4/pi rate is [INTERNAL THEOREM + MEASURED] (sub-coefficients pinned by doubling estimators), not 'rigorous'; ~40 loci said 'rigorous at rank 1'. DELTA #26.
        "pattern": r"rigorous\s+(?:at\s+)?rank[~\s]?\$?1|rigorous\s+rank-1|rigorous\s+\$4/\\pi\$",
        "exempt_if_nearby": withdrawal_marker("rigorous-at-rank-1"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #26, family gate; remediation PI-directed).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
            "docs/topic_to_paper_lookup.md",
            "tests/test_paper40_universal_rate.py",
        ],
        "cited_by": {},
    },
    {
        "id": "sd-ratio-four",
        # The 'Seeley-DeWitt cancellation theorem' (ratio dim(S) at every order; perturbative supertrace vanishes identically) is false: exact S^3 heat traces give ratios 2, -1, 0, ... The 'K/pi lives entirely in the non-perturbative sector' inference fell with it. DELTA #26.
        "pattern": r"=\s*4\s*=\s*\\dim\(S\)|Seeley--DeWitt\s+cancellation\s+theorem|SD\s+cancellation\s+theorem|vanishes\s+identically\s+in\s+the\s+asymptotic\s+expansion|lives\s+\\emph\{entirely\}\s+in\s+the\s+non-perturbative",
        "exempt_if_nearby": r"An earlier|withdrawn|" + withdrawal_marker("sd-ratio-four"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #26, family gate; remediation PI-directed).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
            "docs/topic_to_paper_lookup.md",
        ],
        "cited_by": {},
    },
    {
        "id": "odd-zeta-apery-independence",
        # Apery proved zeta(3) irrational; its Q-linear independence from pi^2 (and zeta(2k+1) outside Q(pi)) is CONJECTURAL. 'independent ... by Apery' and an unqualified 'zeta(2k+1) perp Q(pi)' are retired. DELTA #26.
        "pattern": r"independent\s+of[^.;]{0,40}?by\s+Ap|\\perp\s+\\mathbb\{Q\}\(\\pi\)|Ap(?:é|\\'e|\\'\{e\})ry-style\s+irrationality",
        "exempt_if_nearby": r"conjectur|" + withdrawal_marker("odd-zeta-apery-independence"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #26, family gate; remediation PI-directed).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
            "docs/topic_to_paper_lookup.md",
        ],
        "cited_by": {},
    },
    {
        "id": "sd-coeff-sign-error",
        # The production Seeley-DeWitt function had E = +R/4 (should be -R/4), dropped E from a_1 and omitted the spin-curvature term: a_0 = a_1 = sqrt(pi), a_2 = sqrt(pi)/8 was wrong; correct (4-spinor) a_0 = sqrt(pi), a_1 = -sqrt(pi)/2, a_k = 0 (k >= 2). DELTA #26.
        "pattern": r"a_0\s*=\s*a_1\s*=\s*\\sqrt|a_2\s*=\s*\\sqrt\{?\\pi\}?/8|a_2\$\s+on\s+\$S\^3\s*=\s*\\sqrt|a_0\s*=\s*a_1\s*=\s*sqrt\(pi\)|a_2\s*=\s*sqrt\(pi\)/8|a_0\s*=\s*a_1\s*=\s*\\sqrt\\pi",
        "exempt_if_nearby": r"corrected|sign error|pre-2026-10-02|withdrawn|" + withdrawal_marker("sd-coeff-sign-error"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #26, family gate; remediation PI-directed).",
        "files": [
            "docs/validation_benchmarks.md",
            "docs/algebraic_registry.md",
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
            "docs/topic_to_paper_lookup.md",
            "geovac/*.py",
            "tests/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "m1-census-overreach",
        # M1 placement of the GH rate is a reading; 'sibling normalizations of a single generator', 'predictive engine', 'every GeoVac period', 'GeoVac sits in this corner', 'empirically rejects' (a null PSLQ), and 'gauge structure forces Hermiticity' over-state. DELTA #26.
        "pattern": r"sibling\s+normali[sz]ations?|predictive\s+engine|Every\s+GeoVac\s+period\s+(?:=|is|sits)|places\s+GeoVac'?s\s+periods|GeoVac\s+sits\s+in\s+this\s+corner|loop\s+tower,\s+all\s+of\s+which|empirically\s+rejects|gauge\s+structure\s+forces\s+Hermiticity|Hermiticity\s+arises\s+from\s+the\s+compatibility",
        "exempt_if_nearby": withdrawal_marker("m1-census-overreach"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #26, family gate; remediation PI-directed).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
            "docs/topic_to_paper_lookup.md",
        ],
        "cited_by": {},
    },
    {
        "id": "k-family-surface-d26",
        # K = pi(B+F-Delta) is an Observation: no 'cubic identity' (the cubic holds to 8.8e-8, not exactly), no 'new derivation of alpha' (none exists), no 'explains why pi appears in K', no Paper 2 'Conjecture' tier row. DELTA #26.
        "pattern": r"cubic\s+identity|new\s+derivation\s+of\s+(?:\$\\alpha\$|the\s+fine)|in\s+the\s+Paper~?2\s+sense|explains\s+why\s+\$?\\pi\$?\s+appears\s+in\s+\$?K|^\|[^|]*\|\s*2(?:\s|,|\s*\()[^|]*\|[^|]*\|\s*Conjectur|exact-rational\s+[—-]+\s+OBSERVATION",
        "exempt_if_nearby": withdrawal_marker("k-family-surface-d26"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #26, family gate; remediation PI-directed).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
            "docs/topic_to_paper_lookup.md",
        ],
        "cited_by": {},
    },
    {
        "id": "seam-counts-c2-d26",
        # Paper 57 catalogue is 34 F / 1 Obs / 1 A / 24 C (not 35 F); the Fock coupling prefactor is 1/4 (c^2(4,3) = 1/10, not 1/40 = Delta); the skeleton gives HO shell closures, not the physical magic numbers. DELTA #26.
        "pattern": r"\b35\s*F\s*/|c(?:\^2|²)\(4,\s*3\)\s*=\s*1/40|\(1/16\)\[1\s*[−-]\s*l\(l\+1\)|Fock\s+coupling\s+c²\(n,0\)|[Nn]uclear\s+magic\s+numbers\s+from\s+the\s+(?:packing\s+)?skeleton",
        "exempt_if_nearby": r"superseded|withdrawn|corrected|before 20\d\d|earlier version|" + withdrawal_marker("seam-counts-c2-d26"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #26, family gate; remediation PI-directed).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
            "docs/topic_to_paper_lookup.md",
        ],
        "cited_by": {},
    },
    {
        "id": "propinquity-rate-d26",
        # Paper 38 proves state-space GH (not propinquity); Paper 45's thm:main is a degeneracy theorem; the circle Fejer rate is T log N_t / N_t (not T/N_t or 1/N_t). DELTA #26.
        "pattern": r"propinquity\s+(?:limit|machinery)|SU\(2\)\s+propinquity\s+(?:convergence|bound)|O\(T/\\?N_?t\)|\+\s*T/N_t\)|\+\s*1\s*/\s*N_t\)|C\s*T\s*/\s*N_t",
        "exempt_if_nearby": r"not claimed|not a propinquity|withdrawn|abandoned|degenera|closed|" + withdrawal_marker("propinquity-rate-d26"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-02 (/qa trunk DELTA #26, family gate; remediation PI-directed).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "README.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
            "docs/topic_to_paper_lookup.md",
            "geovac/*.py",
            "tests/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "one-loop-odd-zeta-unscoped",
        # T9 makes the INTEGER-s values zeta_{D^2}(s) pi^even; the one-loop determinant zeta'_{D^2}(0) = 3 zeta(3)/(4 pi^2) + (ln 2)/2 carries zeta(3), and vertex-restricted one-loop sums carry G, beta(4). Pi = 1/(48 pi^2) and beta(alpha) are ENTERED flat-space values, not computed from S^3 data. DELTA #26 self-check.
        "pattern": r"no\s+odd-zeta\s+(?:content\s+)?at\s+one\s+loop|one-loop\s+QED\s+on\s+S(?:\^|³)\{?3?\}?\s+(?:has\s+no\s+odd|is\s+entirely\s+within|produces\s+only)|one-loop\s+quantities\s+on~?\$?S\^3\$?\s+are\s+polynomials|T9\s+theorem\s+constrains\s+\$?c_3|Reproduces\s+standard\s+QED\s+beta|VP\s+from\s+S³\s+spectral\s+data|reproduce\s+the\s+standard\s+flat-space\s+one-loop\s+QED\s+results\s+from",
        "exempt_if_nearby": r"earlier version|withdrawn|integer-\$?s|" + withdrawal_marker("one-loop-odd-zeta-unscoped"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk DELTA #26 self-check, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "docs/validation_benchmarks.md",
            "docs/topic_to_paper_lookup.md",
            "docs/forcing_catalogue.md",
            "geovac/*.py",
            "tests/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "single-cutoff-as-proof-by-argument",
        # FULL #12: Paper 32's single-cutoff spectral action is true by stipulation (an
        # OBSERVATION after inspection); the multi-loop corollary rests on an unargued
        # premise (on-shell scheme a candidate counterexample) and is CONDITIONAL. The
        # C-arc tally is six proof-by-argument results, not seven.
        "pattern": r"seven\s+(?:\\emph\{)?proof-by-argument|single-cutoff\s+spectral\s+action\s+theorem|"
                   r"closes\s+(?:the\s+)?multi-loop\s+QED\s+renormali[sz]ation\s+sub-sector|"
                   r"single-cutoff\s+machinery\s+cannot\s+produce",
        "exempt_if_nearby": withdrawal_marker("single-cutoff-as-proof-by-argument"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #12, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "README.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
        ],
        "cited_by": {"papers/group3_foundations/paper_57_forced_free_seam.tex": "reviewed 2026-10-03 (FULL #12 remediation; self-check rounds 1-2 read the tally sections in full)",
                     "papers/synthesis/geovac_field_guide.tex": "reviewed 2026-10-03 (FULL #12 remediation; self-check rounds 1-2 read the tally sections in full)"},
    },
    {
        "id": "spatial-wall-angular-reproduced",
        # FULL #12: Paper 54 backs Gaunt + m-conservation structurally and a PURE k = 0
        # monopole at n_max in {2,3} only -- the k >= 1 multipoles of Coulomb are absent,
        # so 'angular structure (incl. monopole hierarchy) autonomously reproduced' is false.
        "pattern": r"monopole\s+(?:structure|hierarchy)[^.;]{0,80}(?:reproduc|forced)|"
                   r"angular\s+structure\s+(?:is\s+)?forced|reproduces\s+angular\s+structure|"
                   r"angular\s+content\s+is\s+forced",
        "exempt_if_nearby": withdrawal_marker("spatial-wall-angular-reproduced"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #12, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "README.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
        ],
        "cited_by": {"papers/group3_foundations/paper_57_forced_free_seam.tex": "reviewed 2026-10-03 (FULL #12 remediation; self-check rounds 1-2 read the tally sections in full)"},
    },
    {
        "id": "ko-dim-alternatives-unexhibited",
        # FULL #12: thm:ko_dim_non_selection proves non-selection only; alternative inner
        # data with another KO-dimension satisfying every axiom are NOT exhibited (OPEN).
        "pattern": r"each\s+producing\s+a\s+valid\s+GeoVac|alternative\s+(?:KO-dim\w*|specifications)[^.;]{0,100}each\s+consistent|"
                   r"different\s+inner\s+KO-dim\w*\s+signatures,\s+each",
        "exempt_if_nearby": withdrawal_marker("ko-dim-alternatives-unexhibited"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #12, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "README.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "docs/forced_free_seam.md",
            "docs/forcing_catalogue.md",
        ],
        "cited_by": {},
    },
    {
        "id": "polar-generator-preserves-wedge",
        # FULL #12: K = 2J_z anticommutes with R_polar, so P_W (2J_z) P_W = 0; it does NOT
        # preserve the 'hemispheric' wedge. K_alpha^W (Paper 42 Def 5.1) is a diagonal
        # RELABELLING (two_m_j on the m_j > 0 states), not a restriction of 2J_z.
        "pattern": r"e\^\{-2\\pi\s*J_\{?z\}?\}/Z\s*=\s*e\^\{-\\pi|(?:rotation|boost)\s+generator\s+(?:preserving|that\s+preserves)\s+the\s+(?:hemispheric\s+)?wedge|"
                   r"wedge-preserving\s+\$?SO\(2\)|rotation\s+around\s+the\s+Hopf-base\s+axis\s+preserves\s+the\s+wedge|"
                   r"generator\s+around\s+the\s+wedge\s+axis|its\s+wedge\s+restriction\s+\$?\\KL|"
                   r"inheriting\s+integer\s+eigenvalues|which\s+preserves\s+the\s+wedge\s+by\s+construction|"
                   r"K_alpha\s+restricted\s+to\s+the\s+wedge|"
                   r"K_\\?\{?\\?alpha\}?\^\{?(?:\\alpha,\s*)?W\}?\s*=\s*(?:\\mathrm\{diag\}\(2\s*m_j\)\s*=\s*)?2\s*J_z|"
                   r"K_alpha\^W\s*=\s*2J_z|wedge-preserving\s+boost\s+is\s+played|"
                   r"wedge-block\s+restriction\s+of\s+K_alpha|K_alpha\s+restricted\s+to\s+wedge|"
                   r"natural\s+local\s+Hamiltonian\s+for\s+the\s+wedge",
        "exempt_if_nearby": withdrawal_marker("polar-generator-preserves-wedge"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #12, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "geovac/*.py",
            "tests/*.py",
            "docs/claim_test_matrix.md",
        ],
        "cited_by": {},
    },
    {
        "id": "w1c-tooling-addressable",
        # FULL #12: W1c was BUILT and the NaH PES still fails (W1d/W1e); calling it mechanical/tooling-addressable suppressed a negative.
        "pattern": r"three\s+walls\s+are\s+tooling-addressable|PES\s+extension\s+is\s+mechanical|closing\s+W1a/b/c\s+is\s+internal-effort|W1c[^.;]{0,80}not\s+implemented",
        "exempt_if_nearby": withdrawal_marker("w1c-tooling-addressable"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #12, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "README.md",
            "docs/claims_register.md",
            "docs/forced_free_seam.md",
            "docs/topic_to_paper_lookup.md",
            "docs/forcing_catalogue.md",
        ],
        "cited_by": {},
    },
    {
        "id": "g4a-live-path",
        # FULL #12: G4a is BUILT and CLOSED positive-thin; it is not the live/remaining electroweak path.
        "pattern": r"G4a\s+is\s+the\s+only\s+sprint-scale|only\s+sprint-scale\s+path[^.;]{0,60}G4a|remaining\s+electroweak-unification\s+path\s+is\s+G4a|G4a\s+remains\s+worth\s+pursuing|G4a\s+the\s+remaining",
        "exempt_if_nearby": withdrawal_marker("g4a-live-path"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #12, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "README.md",
            "docs/claims_register.md",
            "docs/forced_free_seam.md",
            "docs/topic_to_paper_lookup.md",
            "docs/forcing_catalogue.md",
        ],
        "cited_by": {},
    },
    {
        "id": "three-homes-structurally-independent",
        # FULL #12: Paper 2's three spectral homes are of DISTINCT TYPE, not 'structurally independent'
        # (B and Delta are rationals lying in every Mellin sub-ring; the exclusive-sub-ring argument fails).
        "pattern": r"three\s+structurally\s+independent\s+(?:objects|projections|spectral)|each\s+have\s+an\s+independent\s+spectral\s+home|three\s+independent\s+spectral\s+homes",
        "exempt_if_nearby": withdrawal_marker("three-homes-structurally-independent"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #12, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "README.md",
            "docs/claims_register.md",
            "docs/forced_free_seam.md",
        ],
        "cited_by": {},
    },
    {
        "id": "kappa-universal-topological-constant",
        # FULL #12 (PI decision 3): kappa = -1/16 is the matched energy scale, an OBSERVATION
        # (Papers 0/7), not a 'universal topological constant'.
        "pattern": r"universal\s+(?:topological|kinetic)\s+(?:constant|scale)|universal\s+constant\s+(?:K\s*=\s*)?-1/16|Origin\s+of\s+the\s+Universal\s+Constant",
        "exempt_if_nearby": withdrawal_marker("kappa-universal-topological-constant"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #12, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "geovac/*.py",
            "geovac-hamiltonians/src/geovac_hamiltonians/*.py",
            "docs/claims_register.md",
            "docs/forced_free_seam.md",
            "README.md",
            "tests/*.py",
            "demo/*.py",
            "benchmarks/scripts/*.py",
            "docs/topic_to_paper_lookup.md",
        ],
        "cited_by": {},
    },
    {
        "id": "fejer-constant-as-gh-rate",
        # FULL #13 (2026-10-03): 4/pi is the Fejer kernel's moment constant; what is proved is
        # the bound d_GH <= gamma_n.  The GH convergence rate is not 4/pi (a taper gives O(1/n),
        # MEASURED).  Retires 'asymptotic/universal state-space GH rate' and 'explicit asymptotic
        # rate' for the convergence (P32 conclusion, P55 table/subsection).
        "pattern": r"(?:asymptotic|universal)\s+state-space\s+GH\s+rate|explicit\s+asymptotic\s+rate|"
                   r"rate\s+constant\s+\$?4/\\pi\$?\s+in\s+the\s+van|with\s+asymptotic\s+rate\s+\$?4/\\pi|"
                   r"GH\s+metric\s+rate|sense,\s+with\s+rate\s+constant\s+\$?4/\\pi|GH-rate\s+entry|Paper~38\s+rate\s+is\s+placed|"
                   r"sets\s+the\s+L2\s+GH-convergence\s+rate|universal\s+rate\s+4/π|the\s+4/π\s+rate\s+is|"
                   r"rate\s+now\s+established|independent\s+rate\s+coincidence",
        "exempt_if_nearby": withdrawal_marker("fejer-constant-as-gh-rate"),
        "severity": "fail",
        "scope": "trunk group1 group2 group3 group4 group5 group6 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #13, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "README.md",
            "docs/*.md",
            "geovac/*.py",
            "tests/*.py",
        ],
        "cited_by": {"papers/group3_foundations/paper_55_periods_of_geovac.tex": "remediated 2026-10-03 (FULL #13; table row + canonical example 3 reread in full)",
                     "papers/group1_operator_algebras/paper_32_spectral_triple.tex": "remediated 2026-10-03 (FULL #13; Conclusion reread in full)",
                     "papers/group4_quantum_computing/paper_20_resource_benchmarks.tex": "remediated 2026-10-03 (FULL #13 self-check; metadata paragraph reread)",
                     "papers/group5_qed_gauge/paper_41_rule_b_wilson_u1.tex": "remediated 2026-10-03 (FULL #13 self-check; Relation-to-38/40 paragraph reread)",
                     "papers/group5_qed_gauge/paper_51_gravity_arc.tex": "remediated 2026-10-03 (FULL #13 self-check round 2; cigar-entropy paragraph reread)"},
    },
    {
        "id": "multi-loop-corollary-confirmed-closure",
        # FULL #13: Paper 32's multi-loop renormalization corollary is CONDITIONAL on an
        # unargued premise and the single-cutoff statement is an OBSERVATION; LS-8a is
        # consistent with it, not a confirmation, and Paper 57 records a conditional
        # classification, not a closure.
        "pattern": r"empirically\s+confirmed\s+this\s+structural\s+verdict|"
                   r"renormali[sz]ation\s+sub-sector\s+closure|"
                   r"record\s+this\s+theorem\s+as\s+the\s+renormali[sz]ation|"
                   r"wall\s+is\s+now\s+bound\s+by\s+argument|now\s+bound\s+by\s+argument\s+at\s+the\s+canonical-rep\s+level\s+under\s+two",
        "exempt_if_nearby": withdrawal_marker("multi-loop-corollary-confirmed-closure"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #13, family gate).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "README.md",
            "docs/*.md",
        ],
        "cited_by": {"papers/group3_foundations/paper_57_forced_free_seam.tex": "reviewed 2026-10-03 (FULL #13; L613/L621/L631 already conditional / Observation)"},
    },
    {
        "id": "lcao-bridge-graph-validated",
        # FULL #13: the LCAO bridged-graph (graph concatenation) H2 model is superseded and
        # guard-railed (FCI-M: R-independent kinetic, no equilibrium).  Its bridge counts and
        # optimized-R '0.43%' are not validated physics.
        "pattern": r"experimentally\s+validated\s+for\s+H|H.?\s+optimal\s+range|QUANTITATIVE\s+ACCURACY\s+ACHIEVED|"
                   r"Bridges\s*=\s*molecular\s+bonds|Binding\s+energy\s*=\s*Eigenvalue\s+lowering|validated\s+with\s+H.?.?\+|confirms\s+topology\s+is\s+correct",
        "exempt_if_nearby": withdrawal_marker("lcao-bridge-graph-validated"),
        "severity": "fail",
        "scope": "trunk group2",
        "note": "Registered 2026-10-03 (/qa trunk FULL #13, family gate).",
        "files": [
            "geovac/*.py",
            "demo/*.py",
            "benchmarks/scripts/*.py",
            "tests/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "graph-model-credited-with-hydrogen-dynamics",
        # FULL #13: the Rabi / dynamics tests run hamiltonian="graph" (non-physical excited
        # states); they check the integrator, not hydrogen.
        "pattern": r"lattice\s+produces\s+quantitatively\s+correct|Build\s+the\s+hydrogen\s+system|"
                   r"converge\s+to\s+1:4:9:16|property\s+of\s+the\s+graph'?s\s+continuum\s+limit|"
                   r"unique\s+constant\s+that\s+maps\s+these\s+dimensionless|encode\s+the\s+nuclear-electron\s+attraction",
        "exempt_if_nearby": withdrawal_marker("graph-model-credited-with-hydrogen-dynamics"),
        "severity": "fail",
        "scope": "trunk group3 synthesis",
        "note": "Registered 2026-10-03 (/qa trunk FULL #13, family gate).",
        "files": [
            "tests/*.py",
            "geovac/*.py",
            "demo/*.py",
            "docs/*.md",
            "benchmarks/scripts/*.py",
        ],
        "cited_by": {},
    },
    {
        "id": "trunk-small-sweep-2026-10-03",
        # Retired wordings from the trunk SMALL sweep after FULL #13 (one backstop entry; each
        # alternative fired on HEAD and has zero live hits after the sweep).  1/Omega^4 is the
        # fourth power of the inverse conformal factor (the Jacobian is Omega^3); a negative real
        # plaquette product gives a pi Berry phase; P0's choices are motivated, not forced; M1 is
        # the volume factor (4/pi only shares the value); W1 walls are GeoVac-internal; CvS P_n
        # is a spectral projection; G2/G3 are related, not one fact; Paper 28 is a reading; S^3
        # (odd-dim) has no almost-complex structure; J flips m_j, not kappa; dynamics.py radial
        # factors are heuristic.
        "pattern": r"inverse\s+Fock\s+Jacobian|Fock\s+Jacobian\s+factor|nonzero\s+Berry\s+phase\s+requires\s+complex|"
                   r"Two\s+axioms\s+force\s+shell\s+capacity|collapses\s+to\s+a\s+sphere-volume\s+quotient|"
                   r"each\s+of\s+these\s+is\s+a\s+frontier-of-field|Bloch--?Siegert-corrected\s+two-level\s+period|"
                   r"basis-index\s+cutoff|one\s+structural\s+fact,\s+not\s+two|"
                   r"it\s+is\s+the\s+spectral\s+action\s+of\s+the\s+present\s+triple|GeoVac\s+is\s+rank~1|"
                   r"rate-controlled\s+restatement|no\s+(?:integrable|compatible)\s+almost-complex|J\s+flips\s+kappa|"
                   r"and\s+hydrogen\s+radial\s+integrals",
        "exempt_if_nearby": withdrawal_marker("trunk-small-sweep-2026-10-03"),
        "severity": "fail",
        "scope": "trunk group1 group3 group5 synthesis",
        "note": "Registered 2026-10-03 (trunk SMALL sweep after /qa FULL #13).",
        "files": [
            "papers/group*/*.tex",
            "papers/synthesis/*.tex",
            "geovac/*.py",
            "tests/*.py",
            "docs/claim_test_matrix.md",
            "README.md",
        ],
        "cited_by": {},
    },
    {
        "id": 'p59-chain-one-centre-densities',
        "pattern": 'three\\s+one-cent(?:er|re)\\s+\\$?1s\\$?\\s+Slater\\s+densit',
        "exempt_if_nearby": withdrawal_marker('p59-chain-one-centre-densities'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). Paper 59 sec:three_electron defined the chain with one-centre densities (the rational case).',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'genus-jump-bars-finite-closed-form',
        "pattern": 'genus\\s+jump\\s+bars\\s+(?:a\\s+finite\\s+closed\\s+form|it\\b)|finite\\s+closed\\s+form\\s+is\\s+barred|closed\\s+form\\s+is\\s+barred\\s+by|closed\\s+forms\\s+stop\\s+at\\s+genus\\s+zero|a\\s+closed\\s+form\\s+past\\s+genus\\s+zero|rules\\s+out\\s+a\\s+closed\\s+form\\s+of\\s+the\\s+two-cent|closed\\s+forms?\\s+of\\s+the\\s+two-cent(?:re|er)\\s+kind',
        "exempt_if_nearby": withdrawal_marker('genus-jump-bars-finite-closed-form'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). J3/J6 walls and Paper 59 sec:scope overstated a bar on any finite closed form.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'bridge-transcendence-class',
        "pattern": 'by\\s+J9\\s+(?:its|the)\\s+class|two-density\\s+one\\s+\\(genus|two-cent(?:re|er)\\s+class\\s+of\\s+door\\s+11|weight-1\\s+\\(one\\s+kernel\\)|bridge\\s+has\\s+one\\s+kernel|has\\s+one\\s+kernel\\s+and\\s+sits\\s+on|one\\s+kernel\\s*&\\s*three-body,\\s+then\\s+two\\s*&\\s*the\\s+two-cent|borne\\s+out\\s+constructively\\s+at\\s+two\\s+centers,\\s+where\\s+the\\s+genus\\s+is\\s+zero',
        "exempt_if_nearby": withdrawal_marker('bridge-transcendence-class'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). An unbacked transcendence class for the four-body r12 bridge (seam J8, Paper 57, atlas, field guide).',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p18-three-centre-elliptic-polylog',
        "pattern": 'full\\s+integral\\s+is\\s+an\\s+elliptic\\s+polylogarithm|closed\\s+form\\s+in\\s+elliptic\\s+polylogarithms\\s+is\\s+not\\s+derived|an\\s+elliptic\\s+polylogarithm,\\s+a\\s+new\\s+ring|elliptic\\s+polylog\\$\\\\\\}\\s+at\\s+genus',
        "exempt_if_nearby": withdrawal_marker('p18-three-centre-elliptic-polylog'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). Paper 18 sec:level2_genus restated the three-centre integral as an elliptic polylog.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'paper35-equivalence-attribution',
        "pattern": 'Paper~?\\s?35\\s+states[^.]{0,140}?\\bas\\s+equivalent|Paper~?\\s?35\\s+states\\s+as\\s+equivalent|restatement\\s+is\\s+a\\s+PI\\s+call|Restating\\s+the\\s+claim\\s+is\\s+a\\s+PI\\s+call|\\(equivalently:\\s*finite-time\\s+observation\\s+window|is\\s+structurally\\s+what\\s+an\\s+observer\\s+does\\s+when\\s+integrating|and\\s+compactified\\s+time\\s+is\\s+discrete',
        "exempt_if_nearby": withdrawal_marker('paper35-equivalence-attribution'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). Mirror-direction staleness: citers kept the pre-restatement Paper 35 status.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'rung-rule-paraphrase',
        "pattern": 'kernels\\s+sharing\\s+a\\s+vertex,\\s+equivalently|rung\\s+(?:variable\\s+)?is\\s+(?:set\\s+by\\s+)?the\\s+number\\s+of\\s+Coulomb\\s+kernels|(?:class|genus)\\s+counts\\s+the\\s+densities\\s+a\\s+kernel|two-cent(?:er|re)\\s+transition\\s+densities---zero\\s+or|indexed\\s+by\\s+kernels\\s+sharing\\s+a\\s+vertex|genus\\s+\\$=\\$\\s+densities\\s+coupled|J9,\\s+the\\s+rung\\s+rule\\s+\\(F\\)|\\|\\s*\\*\\*F\\*\\*\\s*\\[SYMBOLIC\\s+at\\s+one\\s+centre|FORCED\\s+\\(symbolic\\s+at\\s+one\\s+centre|stay\\s+genus\\s+zero\\s+\\(Sec\\.~\\\\ref\\{sec:genus\\}\\)|adds\\s+no\\s+genus,\\s+because\\s+the\\s+class\\s+counts|adds\\s+nothing,\\s+because\\s+genus\\s+counts|no\\s+genus\\s+\\[MEASURED\\],\\s+because|transcendence\\s+class\\s+is\\s+thus\\s+set\\s+by|centers:\\\\\\s+the\\s+transcendence\\s+class\\s+is\\s+set|The\\s+genus\\s+is\\s+set\\s+by\\s+the\\s+number|genus\\s+is\\s+set\\s+by\\s+the\\s+number\\s+of\\s+distinct|densities\\s+set\\s+the\\s+genus|genus\\s+set\\s+by\\s+distinct|the\\s+two\\s+that\\s+sit\\s+at\\s+Observation\\s+status\\s+\\(J2,\\s+J8\\)|consistent\\s+with\\s+genus\\s+counting',
        "exempt_if_nearby": withdrawal_marker('rung-rule-paraphrase'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). Rung rule stated as kernel-count-equivalent and Forced; re-tiered J9 F -> Obs.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p55-ring-nesting-equality',
        "pattern": 'M_1\\s*\\\\cap\\s*\\\\MT\\(\\\\Z\\[i,\\s*1/2\\]\\)\\s*=\\s*\\\\Q\\[\\\\pi\\]|M_1\\s*\\\\cdot\\s*\\(M_0\\s*\\\\oplus\\s*M_2\\s*\\\\oplus\\s*M_3\\)|M_1\\s*\\\\cdot\\s*M_3\\s*\\\\subset\\s*\\\\MT\\(\\\\Z\\[i,\\s*1/2\\]\\)(?!\\[)',
        "exempt_if_nearby": withdrawal_marker('p55-ring-nesting-equality'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). Paper 55 ring nesting kept an equality and two table inclusions the theorems do not give.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'catalan-via-half-integer-hurwitz',
        "pattern": 'Catalan~?\\$?G\\$?\\s+via\\s+the\\s+half-integer|half-int\\.\\\\?\\s+Hurwitz\\s+\\$\\\\to\\$\\s+odd~\\$\\\\zeta\\$;\\s+\\$G\\$',
        "exempt_if_nearby": withdrawal_marker('catalan-via-half-integer-hurwitz'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). Catalan G is quarter-integer vertex-parity (M3) content, not half-integer Hurwitz.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 's5-m3-half-integer-shift',
        "pattern": 'half-integer\\s+shift\\s+\\$?5/2|base\\s+measure\\s+on\\s+the\\s+Hopf',
        "exempt_if_nearby": withdrawal_marker('s5-m3-half-integer-shift'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": "Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). Paper 55 S^5 paragraph: the S^5 M3 shifts are 5/4, 7/4; 'base measure on the Hopf' is the retired label.",
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'engine-k1-plain-trace',
        "pattern": '\\\\mathrm\\{Tr\\}\\(D\\s*\\\\cdot\\s*e\\^\\{-tD\\^2\\}\\)\\$\\s+produces|arises\\s+by\\s+Mellin\\s+transform\\s+of\\s+\\$\\\\Tr\\(D\\s+e\\^\\{-tD\\^2\\}\\)\\$,\\s+with|\\\\Mcal\\[\\\\Tr\\(D\\s*\\\\cdot\\s*e\\^\\{-tD\\^2\\}\\)\\]\\(s\\)\\$,\\s+extracts|commits\\s+to\\s+the\\s+\\$\\\\eta\\$-invariant\\s+Mellin\\s+family|\\\\Tr\\\\bigl\\(D\\^k\\\\,\\s*e\\^\\{-tD\\^2\\}\\\\bigr\\)\\\\right\\]\\(s\\)|\\\\Mcal\\[\\\\Tr\\(D\\^k\\s*e\\^\\{-tD\\^2\\}\\)\\]|\\\\Mcal\\\\!\\\\left\\[\\\\,\\\\Tr\\\\bigl\\(D\\^k\\\\,\\s*e\\^\\{-tD\\^2\\}\\\\bigr\\)|Mellin\\s+slot\\s+of\\s+\\$\\\\Tr\\(D\\^k\\s+e\\^\\{-tD\\^2\\}\\)\\$|\\\\mathcal\\{M\\}\\[\\\\mathrm\\{Tr\\}\\(D\\^k\\s+e\\^\\{-tD\\^2\\}\\)\\]\\$\\s+at\\s+\\$k|\\\\mathcal\\{M\\}\\[\\\\Tr\\(D\\^\\{k\\}\\s+\\\\cdot\\s+e\\^\\{-t\\s+D\\^\\{2\\}\\}\\)\\]',
        "exempt_if_nearby": withdrawal_marker('engine-k1-plain-trace'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). The headline k=1 slot was stated as the plain trace, which vanishes on the symmetric spectrum.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'k3-period-product-not-sum',
        "pattern": 'product\\s+of\\s+two\\s+elliptic\\s+ones|product\\s+of\\s+two\\s+elliptic\\s+periods',
        "exempt_if_nearby": withdrawal_marker('k3-period-product-not-sum'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). Bailey/BBBG eq. 33 is a SUM of two products of complete elliptic integrals.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p55-period-map-isomorphism',
        "pattern": 'proves\\s+the\\s+period\\s+map\\s+is\\s+an\\s+isomorphism|no\\s+new\\s+motivic\\s+conjecture\\s+required|the\\s+period\\s+map\\s+\\$\\\\mathrm\\{per\\}:\\s*\\\\Hcal\\^N\\s*\\\\to\\s*\\\\Zcal\\^N\\$\\s+is\\s+an\\s+isomorphism|GeoVac\\s+IS\\s+the\\s+abelianis|is\\s+precisely\\s+the\\s+Deligne--Glanois|literal\\s+Deligne--Glanois\\s+Galois|\\[Vertex\\s+parity\\s+is\\s+Deligne--Glanois\\s+descent\\]',
        "exempt_if_nearby": withdrawal_marker('p55-period-map-isomorphism'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": "Registered 2026-10-16 (/qa claim-search DELTA on the v5.20.16-18 claim-impact set). Deligne/Glanois give motivic bases whose periods SPAN the real values; injectivity is the period conjecture (Glanois arXiv:1411.4947 abstract: 'a generating family').",
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p59-chi4-only-topological-string',
        "pattern": 'nearest\\s+\\$\\\\chi_\\{-4\\}\\$\\s+appearance\\s+being|levels/conductors\\s+\\$3,\\s*6,\\s*8,\\s*15\\$',
        "exempt_if_nearby": withdrawal_marker('p59-chi4-only-topological-string'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (/qa claim-search DELTA on the v5.20.21 claim-impact set). chi_{-4} enters the twisted Bessel moments as the twist 12.3.c.a -> 48.3.e.a for c_{3,0} (Lim-Tu-Yu; LMFDB: CM field stays Q(sqrt-3)); what stays absent from the cited sources is the Q(i) / Catalan content, not chi_{-4} itself.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p59-loop-mechanism-exact-integration',
        "pattern": 'loop\\s+lowers\\s+the\\s+rung\\s+when\\s+it\\s+integrates\\s+out\\s+exactly|whose\\s+loop\\s+does\\s+not(\\s+integrate\\s+out\\s+exactly)?\\s*,|drop\\s+comes\\s+from\\s+the\\s+exact\\s+integration\\s+of\\s+the\\s+Coulomb\\s+loop|because\\s+its\\s+loop\\s+integrates\\s+out\\s+exactly',
        "exempt_if_nearby": withdrawal_marker('p59-loop-mechanism-exact-integration'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (/qa claim-search DELTA on the v5.20.21 claim-impact set). The geminal loop also integrates out exactly (a solid angle, Paper 59 sec:density_ring) and does not lower the rung; the discriminator is uniqueness, which factors the Coulomb loop into pairwise inverse distances.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p57-every-tree-one-centre-rational',
        "pattern": 'one-cent(re|er)\\s+limit\\s+of\\s+every\\s+tree\\s+integral',
        "exempt_if_nearby": withdrawal_marker('p57-every-tree-one-centre-rational'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (/qa claim-search DELTA on the v5.20.21 claim-impact set). Paper 59 shows a rational one-centre limit for the three-density chain only; for four densities it cites closed forms (Fromm-Hill, Padhy), not rationality.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'seam-catalogue-stale-total',
        "pattern": '\\b4[01]\\s+forced\\s+entries|60-entry\\s+forced/free\\s+catalogue|seventy-(three|four|six)\\s+named\\s+entries|\\b7[46]-entry\\b|\\b7[46]\\s+entries',
        "exempt_if_nearby": withdrawal_marker('seam-catalogue-stale-total'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (/qa claim-search DELTA on the v5.20.21 claim-impact set). The seam catalogue has 74 entries (41 F / 4 Obs / 1 A / 25 C / 3 W) since J11 (2026-10-17). Widened 2026-10-18: 76 entries (43 F / 4 Obs / 1 A / 25 C / 3 W) since J12-J13; then 77 (26 C) since J14.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p59-whole-tree-ladder',
        "pattern": 'whole\\s+tree\\s+ladder|relation\\s+to\\s+four\\s+is\\s+open',
        "exempt_if_nearby": withdrawal_marker('p59-whole-tree-ladder'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (/qa claim-search DELTA on the v5.20.21 claim-impact set). The Hadamard/banana relation is checked in exact arithmetic for n = 2 to 5 only; n = 5 gives no new type (its operator is printed with AESZ #39).',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'ring-unscoped-coulomb',
        "pattern": '(?<!Coulomb\\s)\\bring\\s+of\\s+three\\s+densities\\s+(sits|does\\s+not\\s+climb)|\\ba\\s+ring\\s+does\\s+not\\s+climb|\\ba\\s+ring\\s+behaves\\s+differently|\\ba\\s+ring\\s+of\\s+three\\s+does\\s+not\\s+sit|rung\\s+is\\s+set\\s+by\\s+how\\s+many\\s+densities\\s+one\\s+integral\\s+links',
        "exempt_if_nearby": withdrawal_marker('ring-unscoped-coulomb'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (/qa claim-search DELTA on the v5.20.21 claim-impact set). Only the Coulomb ring of three drops to the elliptic rung; the correlated triangle with exponential geminals sits with the chain on the K3 rung at gamma = 1, and the count of linked densities orders trees, not rings.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'correlated-triangle-gamma-scope',
        "pattern": 'correlated\\s+(\\(geminal\\)\\s+)?triangle\\s+(sits\\s+)?on\\s+the\\s+K3\\s+rung(?!\\s+at)|sits\\s+with\\s+the\\s+chain,\\s+on\\s+the\\s+K3\\s+rung|So\\s+the\\s+correlated\\s+triangle\\s+sits|by\\s+local\\s+data\\s+at\\s+\\$\\\\gamma=3/2\\$,\\s+which',
        "exempt_if_nearby": withdrawal_marker('correlated-triangle-gamma-scope'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (/qa claim-search DELTA on the v5.20.21 claim-impact set). The K3 placement is exact at gamma = 1 only; at gamma = 3/2 it rests on local data mod p plus rigidity, with a factor above K3 not excluded.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'rung-rule-unscoped-trees',
        "pattern": '(?<!tree of densities, )the\\s+class\\s+is\\s+set\\s+by\\s+the\\s+number\\s+of\\s+distinct\\s+two-centre\\s+densities\\s+a\\s+kernel\\s+couples,\\s+and\\s+a\\s+second\\s+kernel',
        "exempt_if_nearby": withdrawal_marker('rung-rule-unscoped-trees'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (/qa claim-search DELTA on the v5.20.21 claim-impact set). The rung rule holds for trees; the Coulomb ring of three (J11) shares vertices and drops to the elliptic rung.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'one-centre-rationals-unscoped',
        "pattern": 'At\\s+one\\s+centre\\s+the\\s+classes\\s+are\\s+exact\\s+rationals|one\\s+centre,\\s+rationals\\s+and\\s+\\$\\\\pi\\$\\s+\\(|One\\s+centre:\\\\\\s+rationals\\s+and\\s+\\$\\\\pi\\$\\.',
        "exempt_if_nearby": withdrawal_marker('one-centre-rationals-unscoped'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 1 of the v5.20.21 claim-search DELTA remediation). One-centre TWO-electron classes are exact rationals (Paper 58); the one-centre Coulomb ring of three densities carries dilogarithms (Remiddi; Paper 59 sec:density_ring).',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'rung-rule-table-cell-unscoped',
        "pattern": 'class\\s+set\\s+by\\s+distinct\\s+two-centre\\s+densities\\s+&',
        "exempt_if_nearby": withdrawal_marker('rung-rule-table-cell-unscoped'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 1 of the v5.20.21 claim-search DELTA remediation). The rung rule holds for trees; the Coulomb ring of three (J11) drops to the elliptic rung.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'banana-operators-known-unscoped',
        "pattern": 'the\\s+operators\\s+are\\s+known,\\s+tied\\s+to\\s+the\\s+banana\\s+ladder',
        "exempt_if_nearby": withdrawal_marker('banana-operators-known-unscoped'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 1 of the v5.20.21 claim-search DELTA remediation). The Hadamard/banana relation is checked for n = 2 to 5 only.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'uniqueness-mechanism-untagged',
        "pattern": 'drop\\s+is\\s+(a\\s+property\\s+of\\s+)?uniqueness,\\s+not\\s+(of\\s+)?rings(?!\\s*\\(an\\s+Observation)',
        "exempt_if_nearby": withdrawal_marker('uniqueness-mechanism-untagged'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": "Registered 2026-10-18 (self-check round 1 of the v5.20.21 claim-search DELTA remediation). Uniqueness as the mechanism of the Coulomb ring's drop is an [OBSERVATION] in Paper 59 (one geminal shape, two gamma values); restatements inside F rows must carry the tier.",
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p59-conductor4-content-absent',
        "pattern": 'the\\s+conductor-\\$4\\$\\s+content\\s+this\\s+family\\s+carries',
        "exempt_if_nearby": withdrawal_marker('p59-conductor4-content-absent'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 2 of the v5.20.21 claim-search DELTA remediation). chi_{-4} itself enters the twisted Bessel moments (c_{3,0} via the twist 12.3.c.a -> 48.3.e.a); what does not appear is the Z[i] / CM-by-Q(i) content (theta_3^2, Catalan).',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'uniqueness-mechanism-as-cause',
        "pattern": 'because\\s+its\\s+loop\\s+factors\\s+by\\s+uniqueness',
        "exempt_if_nearby": withdrawal_marker('uniqueness-mechanism-as-cause'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 2 of the v5.20.21 claim-search DELTA remediation). Uniqueness as the mechanism is an [OBSERVATION] (Paper 59 sec:density_ring); a summary may not state it as the cause.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'bridge-reduction-called-exact',
        "pattern": 'reduction\\s+to\\s+geminal-dressed\\s+densities\\s+is\\s+exact|exact\\s+reduction\\s+(does\\s+not\\s+lower|changes\\s+the\\s+cost)',
        "exempt_if_nearby": withdrawal_marker('bridge-reduction-called-exact'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": "Registered 2026-10-18 (self-check round 2 of the v5.20.21 claim-search DELTA remediation). The r12 bridge reduction is validated against Monte Carlo, not derived (Paper 12); and its four-density rung goes through J10's assumed three-dimensional reduction.",
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'one-centre-rationals-and-pi',
        "pattern": 'rationals\\s+and\\s+\\$\\\\pi\\$\\s+for\\s+two-electron',
        "exempt_if_nearby": withdrawal_marker('one-centre-rationals-and-pi'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 2 of the v5.20.21 claim-search DELTA remediation). One-centre two-electron classes are exact rationals (Paper 58 L417).',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'lg-composition-factors-three',
        "pattern": 'composition\\s+factors\\s+of\\s+\\$L_G\\$\\s+are\\s+therefore\\s+elementary',
        "exempt_if_nearby": withdrawal_marker('lg-composition-factors-three'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 2 of the v5.20.21 claim-search DELTA remediation). L_G = M o (theta - r) has two composition factors (algebraic, elliptic); the elementary piece of G is the source term.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'correlated-gamma32-on-k3-rung',
        "pattern": 'by\\s+local\\s+data\\s+at\\s+\\$\\\\gamma=3/2\\$,\\s+where\\s+a\\s+factor\\s+above\\s+K3\\s+is\\s+not\\s+excluded|and\\s+at\\s+\\$\\\\gamma=3/2\\$\\s+by\\s+local\\s+data\\s+mod\\s+\\$p\\$\\s+plus\\s+rigidity,\\s+a\\s+factor\\s+above',
        "exempt_if_nearby": withdrawal_marker('correlated-gamma32-on-k3-rung'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": "Registered 2026-10-18 (self-check round 3 of the v5.20.21 claim-search DELTA remediation). 'Sits on rung X' means no composition factor above X (record 2026-10-17_ring_extensions.md); at gamma = 3/2 the correlated triangle reaches the K3 rung by local data mod p plus rigidity, and whether it climbs further is open.",
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'seam-ladder-header-stale',
        "pattern": 'no\\s+single\\s+paper\\s+states\\s+the\\s+whole\\s+chain|rationals\\s+and\\s+(π|\\$\\\\pi\\$)\\s+at\\s+one\\s+centre',
        "exempt_if_nearby": withdrawal_marker('seam-ladder-header-stale'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": "Registered 2026-10-18 (self-check round 3 of the v5.20.21 claim-search DELTA remediation). One-centre two-electron classes are exact rationals (Paper 58); Paper 57's sub-domain J now states the whole ladder, including J10.",
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'bridge-four-to-three-uncapped',
        "pattern": "four\\s+reduces\\s+to\\s+three''\\),\\s+while",
        "exempt_if_nearby": withdrawal_marker('bridge-four-to-three-uncapped'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": "Registered 2026-10-18 (self-check round 3 of the v5.20.21 claim-search DELTA remediation). The r12 bridge's validated reduction changes the cost, not the value; its unreduced form sits on the four-density rung through J10's assumed 3-D reduction.",
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'gamma-singular-points-meet',
        "pattern": '\\$,\\s+which\\s+meet\\s+at\\s+\\$t=2\\$\\s+for',
        "exempt_if_nearby": withdrawal_marker('gamma-singular-points-meet'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 3 of the v5.20.21 claim-search DELTA remediation). At gamma = 1, |1 - gamma| = 0; three of the four points t = 2, 1 + gamma, 2 gamma meet at t = 2.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'qg-reuses-chain-r3r1',
        "pattern": 'divides\\s+exactly\\s+as\\s+\\$R_3R_1\\$',
        "exempt_if_nearby": withdrawal_marker('qg-reuses-chain-r3r1'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": "Registered 2026-10-18 (self-check round 4 of the v5.20.21 claim-search DELTA remediation). The correlated triangle's order-four annihilator factors as an order-three times an order-one operator that are not the chain control's R_3, R_1 (Y's does not split; the chain's R_3 R_1 does); they need their own names.",
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'fg-five-densities-open',
        "pattern": 'whether\\s+five\\s+linked\\s+densities\\s+go\\s+beyond\\s+four',
        "exempt_if_nearby": withdrawal_marker('fg-five-densities-open'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 4 of the v5.20.21 claim-search DELTA remediation). Paper 59 sec:four_density settles it: five linked densities give no new type (the operator printed with AESZ #39).',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'open-questions-closed-set',
        "pattern": 'are\\s+the\\s+two\\s+open\\s+questions|(?<!Among\\s)the\\s+open\\s+questions\\s+the\\s+ladder\\s+leaves\\s+are\\s+whether|that\\s+is\\s+the\\s+open\\s+question\\s+the\\s+ladder\\s+leaves',
        "exempt_if_nearby": withdrawal_marker('open-questions-closed-set'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": "Registered 2026-10-18 (self-check round 5 of the v5.20.21 claim-search DELTA remediation). The ladder's open questions include whether the correlated triangle climbs at gamma = 3/2 and whether the real 3-D four-density integral reduces to the 1-D family (Paper 57 J10-J11); a closed list omitting them reads them as settled.",
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p57-p1-misclassification-count-unscoped',
        "pattern": 'Misclassifications\\.\\}\\s+Eight\\s+in\\s+all',
        "exempt_if_nearby": withdrawal_marker('p57-p1-misclassification-count-unscoped'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check of the v5.20.21 SMALL sweep). The eight P1 misclassifications are counted on the 60-entry catalogue; the released-axis rows J7, J10 and J11 are three more forced rows at MF > 1.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'c80-parity-definite',
        "pattern": 'so\\s+parity\\s+is\\s+not\\s+the\\s+rule\\s+for\\s+values',
        "exempt_if_nearby": withdrawal_marker('c80-parity-definite'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": "Registered 2026-10-18 (self-check of the v5.20.21 SMALL sweep). c_{8,0}'s negative is a bounded integer-relation search (evidence, not proof); the searches give parity no support as a rule for values.",
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'ising-parity-overclaim',
        "pattern": 'the\\s+value\\s+rule\\s+is\\s+not\\s+parity|is\\s+false\\s+as\\s+a\\s+general\\s+rule|so\\s+the\\s+searches\\s+give\\s+parity\\s+no\\s+support',
        "exempt_if_nearby": withdrawal_marker('ising-parity-overclaim'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 2 of the v5.20.21 SMALL sweep). The registered value cases hold; only the bounded n = 8 search counts against parity, so the n = 8 search gives parity no support as a rule for values (evidence, not a disproof).',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
            'docs/predictions/2026-10-17_ising_class.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p57-four-candidate-principles',
        "pattern": 'four\\s+candidate\\s+(discriminator\\s+)?principles',
        "exempt_if_nearby": withdrawal_marker('p57-four-candidate-principles'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 2 of the v5.20.21 SMALL sweep). Paper 57 tests five candidate discriminator principles (P1-P5); the fifth is introduced after the first four.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
            'docs/predictions/2026-10-17_ising_class.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p57-released-axis-not-multifocal',
        "pattern": 'they\\s+are\\s+not\\s+consumed\\s+and\\s+not\\s+multi-focal',
        "exempt_if_nearby": withdrawal_marker('p57-released-axis-not-multifocal'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 3 of the v5.20.21 SMALL sweep). J7, J10 and J11 sit at MF > 1 (more than one momentum integration), which is not multi-focal composition.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
            'docs/predictions/2026-10-17_ising_class.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'p57-c2-closes-6-1',
        "pattern": 'principle\\s+hunt\\s+that\\s+closes\\s+the\\s+question',
        "exempt_if_nearby": withdrawal_marker('p57-c2-closes-6-1'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 3 of the v5.20.21 SMALL sweep). Paper 57 sec:open_family2 frames the P5 question as resolved-as-tautology plus an open meta-theorem question, not closed.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
            'docs/predictions/2026-10-17_ising_class.md',
        ],
        "cited_by": {},
    },
    {
        "id": 'ising-record-all-three-confirmed',
        "pattern": 'all\\s+three\\s+confirmed\\s+as\\s+registered|IS3:\\s+CONFIRMED\\s+for\\s+the\\s+registered\\s+cases',
        "exempt_if_nearby": withdrawal_marker('ising-record-all-three-confirmed'),
        "severity": 'fail',
        "scope": 'trunk group1 group2 group3 group4 group5 group6 synthesis',
        "note": 'Registered 2026-10-18 (self-check round 3 of the v5.20.21 SMALL sweep). IS3 as registered includes the conifold attachment of c_{4,0}, which was not tested; the registered VALUE cases hold.',
        "files": [
            'papers/group*/*.tex',
            'papers/synthesis/*.tex',
            'docs/atlas.md',
            'docs/forced_free_seam.md',
            'docs/walls/register.md',
            'README.md',
            'papers/INDEX.md',
            'docs/claim_test_matrix.md',
            'docs/predictions/2026-10-17_ising_class.md',
        ],
        "cited_by": {},
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
