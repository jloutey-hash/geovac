r"""Register the retired H2+ "0.0002%" wording in C16, with its dependents.

GATE-FIRST ORDER (qa.md hard rule): register the CLASS, prove it discriminates
two ways, then let the gate enumerate the loci and fix every one.  C17 already
carries the numeric family (`p11-h2plus-0002pct-retired`, 32 live).  C16 is the
prose/zombie half: it catches the CLAIM re-surfacing in words after the number
is gone, e.g. "achieves near-exact accuracy of two parts in a million", which
no numeric pattern sees.

WHY `cited_by` IS THE LOAD-BEARING FIELD HERE.  `files` is only where the
wording might appear.  `cited_by` is the set of documents whose ARGUMENT rests
on the claim, enumerated from the argument rather than from grep -- and the
retraction->dependents rule exists because eight of the ten recurring defect
classes in the v5.4.4..v5.7.3 arc had exactly one shape: a claim corrected in
its owner and left standing in its citers, restated in the citers' own words so
no pattern could reach it.  "Monotone from below" was corrected in Paper 0 and
left false in Paper 7 IN THE SAME COMMIT.

For this claim the dependents are not decorative.  Each of these papers cites
H2+'s accuracy to justify its OWN choice of geometry -- "the two-centre
one-electron problem is essentially solved, so the residual must be
correlation" is the inference Papers 12/13/15/17 build on.  That inference
survives the correction (it gets STRONGER: the one-electron problem is solved
to machine precision, not to 0.0002%), but each locus states it with a number
that is now retired, so each must be revisited and stamped.

Written as a file, not a heredoc: two C17 patterns shipped dead today from
heredoc escape mangling, and the standing rule routes backslash edits through
Write.

Run:  python debug/qa/_register_c16_p11_retraction.py
"""
from __future__ import annotations

import io
import re
import sys

PATH = "debug/qa/check_retracted_terms.py"

ENTRY = '''    {
        "id": "p11-h2plus-0002pct-prose",
        # The CLAIM in words, for when the numeral is gone but the magnitude
        # survives as prose.  Deliberately narrow: only phrasings that assert a
        # ~1e-6-relative accuracy for H2+, which is the retired magnitude.
        "pattern": r"(two parts in a million|2 parts in 10\\^6|"
                   r"two[- ]in[- ]a[- ]million)",
        "require_nearby": r"H\\$?_2\\^?\\{?\\+|H2\\+|prolate|spectral|Laguerre",
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
            "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex": None,
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": None,
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": None,
            "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex": None,
            "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex": None,
            "papers/group2_quantum_chemistry/paper_fci_molecules.tex": None,
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": None,
            "papers/synthesis/geovac_field_guide.tex": None,
        },
    },
'''


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    s = io.open(PATH, encoding="utf-8").read()
    if "p11-h2plus-0002pct-prose" in s:
        print("entry already present; nothing to do")
    else:
        m = re.search(r"^REGISTRY\s*=\s*\[\s*\n", s, re.M)
        if not m:
            raise SystemExit("REGISTRY opening not found")
        s = s[:m.end()] + ENTRY + s[m.end():]
        io.open(PATH, "w", encoding="utf-8", newline="").write(s)
        print("C16: added p11-h2plus-0002pct-prose")

    import importlib.util
    spec = importlib.util.spec_from_file_location("crt", PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    e = [x for x in mod.REGISTRY if x.get("id") == "p11-h2plus-0002pct-prose"][0]
    pat = re.compile(e["pattern"])
    exm = re.compile(e["exempt_if_nearby"])
    assert "\\\\." not in e["pattern"], "DOUBLED ESCAPE in pattern"

    # DISCRIMINATION, both directions -- an entry that fires on nothing is
    # worse than no entry, because the gate then reports PASS.
    retired = "the spectral solver reaches two parts in a million for H2+"
    corrected = "the spectral solver reproduces H2+ to machine precision"
    assert pat.search(retired), "does NOT fire on the retired prose"
    assert not pat.search(corrected), "fires on the CORRECTED prose"
    assert exm.search("note [retracted 2026-09-19: p11-h2plus-0002pct-prose] here"), \
        "per-entry withdrawal marker does not exempt"
    print("  fires on retired prose: yes;  silent on corrected: yes;  marker works")
    print(f"  cited_by declares {len(e['cited_by'])} dependents, all unstamped (correct for a new entry)")


if __name__ == "__main__":
    main()
