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

WINDOW = 5  # +- lines within which a withdrawal marker exempts a hit

# ---------------------------------------------------------------------------
# THE REGISTRY -- append an entry whenever a /qa run retires/withdraws a claim.
# Each entry: a retracted phrase (pattern) that must NOT appear as LIVE; it is
# exempt only when a withdrawal marker (exempt_if_nearby) sits within +-WINDOW
# lines.  Patterns are case-insensitive raw regex.  `files` are ROOT-relative
# globs.  `severity`: "fail" (gates) | "advisory" (reports only).
# ---------------------------------------------------------------------------
REGISTRY = [
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
        "pattern": "(?i)(?!.*\\b(?:no|not|without|irreducible|lacks|lacking)\\b[^.\\n]{0,50}closed)(?!.*closed[- ]form is differential)(?:eigen(?:vector|pair|value)s?[^.\\n]{0,140}?(?:closed[- ]form|closed skeleton form)|(?:closed[- ]form|closed skeleton form)[^.\\n]{0,140}?eigen(?:vector|pair|value)s?)",
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
        "scope": "group3 group6",
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
                   r"|Lorentzian\s+closure\s+is\s+complete",
        "exempt_if_nearby": r"withdrawn|WITHDRAWN|signature-blind|compact[- ]boost"
                            r"|compact\s+KMS|K\^?\+|descope|convention|period[- ]closure"
                            r"|Euclidean|not\s+constitute",
        "files": [
            "papers/group6_precision_observations/paper_34_projection_taxonomy.tex",
            "papers/synthesis/group6_precision_observations_synthesis.tex",
            "papers/group3_foundations/paper_31_universal_coulomb_partition.tex",
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
    exempt = re.compile(entry["exempt_if_nearby"], re.IGNORECASE)
    live, ok = [], []
    for path in _resolve(entry["files"]):
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        for i, line in enumerate(lines):
            if not pat.search(line):
                continue
            lo, hi = max(0, i - WINDOW), min(len(lines), i + WINDOW + 1)
            window_txt = "\n".join(lines[lo:hi])
            rel = path.relative_to(ROOT)
            snip = re.sub(r"\s+", " ", line.strip())[:160]
            if exempt.search(window_txt):
                ok.append((rel, i + 1, snip))
            else:
                live.append((rel, i + 1, snip))
    return live, ok


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    show_all = "--all" in sys.argv
    gate = _gate_substr(sys.argv)
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
    print(f"\nRESULT: PASS (no live fail-severity retracted claim in "
          f"{scope}; {len(_selected)}/{len(REGISTRY)} entries over "
          f"{len(_files_seen)} declared locus pattern(s); {exempt_total} "
          f"occurrence(s) correctly carry a withdrawal flag)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
