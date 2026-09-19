r"""Register the C17 family for the retired graph-native He "0.19% @ n_max=7"
headline, gate-first, and enumerate its loci.

DIAGNOSIS (2026-09-19), which established the sweep is SAFE and not a
mis-correction:
  * build_graph_native_fci (the function the tests use; hybrid h1 = exact
    diagonal + graph-Laplacian off-diagonal + exact rational Slater at k=Z)
    gives, at n_max = 5/6/7: 0.2496% / 0.22864% / 0.21559%.
  * NON-CIRCULAR anchor: paper_fci_atoms states the graph-native error is
    "0.23% at n_max=6" -- which MATCHES the measured 0.22864%.  So the pipeline
    is the right construction, confirmed against the paper's OWN stated value.
  * At n_max=7 the papers say 0.19% (E = -2.8983 Ha); production gives 0.216%
    (E = -2.89746, same dim=1218).  0.19% is a NON-MONOTONE drop from the
    paper's own 0.23% at n_max=6 -- an outlier no current path reproduces, and
    NO-TEST.  So 0.19% @ n_max=7 is stale (pre-correction); the correct value
    is 0.216%.

WHAT MUST NOT BE SWEPT (kept distinct):
  * the adiabatic "0.19-0.20%" FLOOR (a different solver; the pattern excludes
    it because the % is not adjacent to 0.19 in "0.19--0.20\%");
  * the n_max=6 = 0.23% value (correct, matches measurement);
  * dozens of unrelated 0.19x numbers (polarizabilities, D_e, Wilson loops) --
    excluded by require_nearby = graph-native context.

Run:  python debug/qa/_register_c17_graphnative019.py
"""
from __future__ import annotations

import io
import re

PATH = "debug/qa/check_headline_numbers.py"

ENTRY = r'''    {
        "id": "p13-he-graphnative-019-nmax7",
        "scope": "paper_13 paper_18 paper_7 group2 group3 synthesis trunk",
        "severity": "fail",
        "canonical_note": "Registered 2026-09-19. The graph-native He CI "
                          "'0.19% @ n_max=7' (E = -2.8983 Ha) is STALE. Measured "
                          "on the production build_graph_native_fci (hybrid h1, "
                          "exact rational Slater at k=Z): n_max=5/6/7 = "
                          "0.2496/0.22864/0.21559 %. Non-circular anchor: "
                          "paper_fci_atoms's own 'n_max=6 = 0.23%' matches the "
                          "measured 0.22864%, so the pipeline is the right "
                          "construction; the n_max=7 value 0.19% is a "
                          "non-monotone outlier no current path reproduces "
                          "(production E = -2.89746, dim 1218; NO-TEST). Correct "
                          "value: 0.216% (E = -2.89746 Ha). The adiabatic "
                          "'0.19-0.20%' FLOOR is a DIFFERENT solver and is NOT "
                          "this claim (the pattern excludes the range form). "
                          "Driver debug/qa/_graph_native_he_nmax67.py.",
        "pattern": r"0\.19\s*\\?%",
        "require_nearby": r"graph.native|graph.consistent|Graph-native",
        "exempt_if_nearby": r"\[retracted \d{4}-\d{2}-\d{2}:\s*p13-he-graphnative-019-nmax7\]|"
                            r"0\.216|falsifies|MEASURED 2026-09-19|retired|superseded|"
                            r"stale|0\.19\s*-+\s*\$?-*\$?\s*0\.20",
        "files": [
            "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex",
            "papers/group2_quantum_chemistry/paper_fci_atoms.tex",
            "papers/group3_foundations/paper_18_exchange_constants.tex",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex",
            "papers/INDEX.md",
            "docs/claims_register.md",
            "docs/claim_test_matrix.md",
            "docs/qa/group2.done.md",
            "docs/qa/synthesis.done.md",
            "CLAUDE.md",
        ],
    },
'''


def main() -> None:
    s = io.open(PATH, encoding="utf-8").read()
    if "p13-he-graphnative-019-nmax7" in s:
        print("family already present")
    else:
        m = re.search(r"^REGISTRY\s*=\s*\[\s*\n", s, re.M)
        assert m, "REGISTRY opening not found"
        s = s[:m.end()] + ENTRY + s[m.end():]
        io.open(PATH, "w", encoding="utf-8", newline="").write(s)
        print("C17: added p13-he-graphnative-019-nmax7")

    import importlib.util
    spec = importlib.util.spec_from_file_location("chn", PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    e = [x for x in mod.REGISTRY if x["id"] == "p13-he-graphnative-019-nmax7"][0]
    pat, req, exm = (re.compile(e["pattern"]), re.compile(e["require_nearby"]),
                     re.compile(e["exempt_if_nearby"]))
    assert "\\\\." not in e["pattern"], "DOUBLED ESCAPE in pattern"

    # discrimination, both directions
    fires = "graph-native FCI ... achieves 0.19\\% at $n_{\\max}=7$"
    silent_corr = "graph-native FCI ... achieves 0.216\\% at $n_{\\max}=7$"
    silent_floor = "the adiabatic structural floor of 0.19--0.20\\%"
    assert pat.search(fires), "does not fire on the retired form"
    assert not pat.search(silent_corr), "fires on the corrected 0.216"
    assert not pat.search(silent_floor), "fires on the adiabatic FLOOR range"
    print(f"  discrimination: fires on 0.19%%={bool(pat.search(fires))}  "
          f"silent on 0.216={not pat.search(silent_corr)}  "
          f"silent on floor-range={not pat.search(silent_floor)}")

    # enumerate live loci (require_nearby present, exempt absent, in window)
    win = getattr(mod, "WINDOW", 3)
    print(f"  WINDOW={win}; enumerating live loci:")
    n = 0
    for f in e["files"]:
        try:
            lines = io.open(f, encoding="utf-8", errors="replace").read().splitlines()
        except FileNotFoundError:
            print(f"    (missing: {f})"); continue
        for k, ln in enumerate(lines):
            if not pat.search(ln):
                continue
            ctx = "\n".join(lines[max(0, k - win):k + win + 1])
            if req.search(ctx) and not exm.search(ctx):
                print(f"    LIVE  {f}:{k+1}  {ln.strip()[:95]}")
                n += 1
    print(f"  total live: {n}")


if __name__ == "__main__":
    main()
