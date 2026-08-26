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

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
WINDOW = 3  # +/- lines for the exemption window

GROUP4_FILES = [
    "papers/group4_quantum_computing/*.tex",
    "papers/synthesis/group4_quantum_computing_synthesis.tex",
]

REGISTRY = [
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
    exempt = re.compile(entry["exempt_if_nearby"], re.IGNORECASE)
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
            (ok if exempt.search(window_txt) else live).append((rel, i + 1, snip))
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

    def selected(e: dict) -> bool:
        return gate is None or e["scope"] == "all" or gate in e["scope"]

    n_live, n_exempt = 0, 0
    print(f"headline-number registry gate (C17)   [{scope}]\n")
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
        print(f"\nRESULT: FAIL ({n_live} live wrong-headline occurrence(s) in {scope})")
        return 1
    print(f"\nRESULT: PASS (no live wrong headline value in {scope}"
          + (f"; {n_exempt} exempt/historical mention(s))" if n_exempt else ")"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
