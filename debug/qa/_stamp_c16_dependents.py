r"""Stamp the eight `cited_by` dependents of `p11-h2plus-0002pct-prose` as
revisited, recording WHAT was done at each.

WHY STAMPING IS NOT A FORMALITY.  The retraction->dependents rule exists
because eight of the ten recurring defect classes in the v5.4.4..v5.7.3 arc had
one shape: a claim corrected in its owner and left standing in its citers,
restated in the citers' own words so no pattern could reach it.  "Monotone from
below" was corrected in Paper 0 and left false in Paper 7 IN THE SAME COMMIT.
`cited_by` is enumerated from the ARGUMENT, and the gate fails while any
dependent is unstamped -- so stamping is the review outcome, not the intention.

The rule is explicit that an edit is not a review ("stamp outcomes, not
intentions").  Each stamp below therefore names the specific locus changed in
that document during the 2026-09-19 sweep, which is the evidence that the
dependent was actually read.  Every one of these eight was edited in
`_sweep_p11_h2plus_headline.py` (23 anchored replacements, all applied).

WHAT THE DEPENDENTS' ARGUMENT RESTS ON, and why it SURVIVES: each cites H2+'s
accuracy to justify its own geometry choice -- "the two-centre one-electron
problem is essentially solved, so the residual must be correlation."  That
inference is unchanged and in fact strengthened: the one-electron problem is
solved to the reference's own precision, not to 0.0002%.  Only the numeral
moved, so no dependent's conclusion is withdrawn.

Run:  python debug/qa/_stamp_c16_dependents.py
"""
from __future__ import annotations

import io
import re
import sys

PATH = "debug/qa/check_retracted_terms.py"
STAMP = "reviewed 2026-09-19"

# dependent -> what was actually changed there (evidence of the read)
OUTCOMES = {
    "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex":
        "abstract L34 + intro L120 + hierarchy L1349 -> machine precision",
    "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex":
        "L101 prose + tab:hierarchy Level-2 cell + caption (marker added)",
    "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex":
        "L102 prose + L1418 table; FD figure 0.70% -> 1.01% per P11's own table",
    "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex":
        "L72 prose + L1542 table; FD figure 0.70% -> 1.01%",
    "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex":
        "L1349 -> machine precision; the R_eq 2.001/0.21% FD claim left intact",
    "papers/group2_quantum_chemistry/paper_fci_molecules.tex":
        "L766 -> machine precision; guardrail negative untouched",
    "papers/synthesis/group2_quantum_chemistry_synthesis.tex":
        "L87, L169, L279 (the -0.6026 vs -0.6026 form), L1068",
    "papers/synthesis/geovac_field_guide.tex":
        "L224 hierarchy table cell",
}


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    s = io.open(PATH, encoding="utf-8").read()
    i = s.find('"id": "p11-h2plus-0002pct-prose"')
    if i < 0:
        raise SystemExit("entry not found")
    j = s.find("\n    },\n", i)
    block = s[i:j]

    n = 0
    for dep, outcome in OUTCOMES.items():
        old = f'"{dep}": None,'
        if old not in block:
            print(f"  MISS (already stamped or path differs): {dep}")
            continue
        block = block.replace(old, f'"{dep}": "{STAMP} -- {outcome}",', 1)
        n += 1
    s = s[:i] + block + s[j:]
    io.open(PATH, "w", encoding="utf-8", newline="").write(s)
    print(f"stamped {n}/{len(OUTCOMES)} dependents")

    import importlib.util
    spec = importlib.util.spec_from_file_location("crt", PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    e = [x for x in mod.REGISTRY
         if x.get("id") == "p11-h2plus-0002pct-prose"][0]
    unstamped = [k for k, v in e["cited_by"].items() if not v]
    print(f"  remaining unstamped: {unstamped or 'none'}")


if __name__ == "__main__":
    main()
