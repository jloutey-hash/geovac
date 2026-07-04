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
  THE REGISTRY (below) lists each retracted claim as {pattern, exempt_if_nearby,
  files, severity, scope}. For every pattern hit the screen checks whether a
  withdrawal marker (WITHDRAWN / retracted / "false" / Erratum / "state-space GH"
  / "named gap" / ...) appears within +-WINDOW lines. A hit WITHOUT a nearby marker
  is a live zombie:
    * severity "fail"     -> FAILS the gate (exit 1) when in --gate scope;
    * severity "advisory" -> printed for review, does NOT fail (used for the noisier
      classes, e.g. bare "propinquity", where legit mentions abound).

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

import pathlib
import re
import sys

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
        "severity": "advisory",   # noisy: legit framework/named-gap/descope mentions abound
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

    def selected(e: dict) -> bool:
        return gate is None or e["scope"] == "all" or gate in e["scope"]

    fail_hits, advisory_hits, exempt_total = [], [], 0
    print(f"retracted-claims / zombie-drift screen   [{scope}]\n")
    for e in REGISTRY:
        if not selected(e):
            continue
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
        print(f"\nRESULT: FAIL ({len(fail_hits)} live retracted claim(s) in {scope})")
        return 1
    print(f"\nRESULT: PASS (no live fail-severity retracted claim in {scope}; "
          f"{exempt_total} occurrence(s) correctly carry a withdrawal flag)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
