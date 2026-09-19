r"""Widen `require_nearby` on the p11 H2+ family, and prove it now sees the
loci it was missing.

WHY.  The first version of the family matched `0.0002\%` but gated on a
context alternation that omitted the very tokens the worst loci use.  Measured
misses (WINDOW = 3):

  * papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex:492
        `$E_{\mathrm{total}}$ error & 1.01\% & 0.0002\% \\`
    The table cell that prints the wrong energy error AND the wrong R_eq side
    by side -- the single most damning locus in the corpus.  Its nearest
    context token is the column header `Spectral ($N_b = 20$)` two lines up;
    my alternation had `spectral Laguerre` but not bare `Spectral`, and no
    `N_b`, and no `tab:spectral_convergence`.

  * docs/paper_notes_archive.md:108, docs/project_closeout_plan.md:26,
    docs/development_frontier_archive.md:190, docs/qa/synthesis.done.md:185
    Plain-prose "H2+ 0.0002%" with no LaTeX `H$_2^+$` and often no context
    token inside +-3 lines.

So the family under-reported by five loci: 27 named, >=32 real.  A pattern
that silently drops the worst locus is the same defect as one that fires on
nothing -- just quieter, because the gate prints a number and the number looks
like diligence.

FIX.  Add the table/archive vocabulary to `require_nearby`, and add `H2\+`
without LaTeX plus `H_2^+`/`H₂⁺` spellings.  Keep it a CONTEXT gate rather
than dropping it: the four innocent loci (a `within 0.0002 Ha` He comparison
in paper_14, an ERI value beginning 0.00025 in the certified-reference tables,
a cc-pVTZ tolerance in validation_benchmarks) must stay silent, and they will,
because none of them sits near H2+/prolate/spectral vocabulary.

Written as a file, not a heredoc: two patterns already shipped dead today from
heredoc escape mangling.

Run:  python debug/qa/_widen_c17_p11_context.py
"""
from __future__ import annotations

import io
import re

PATH = "debug/qa/check_headline_numbers.py"

WIDER_REQUIRE = (
    r'r"H\$?_2\^?\{?\+|H2\+|H₂⁺|H_2\^\+|prolate|Prolate|spectral|Spectral|'
    r'Laguerre|N_b|n_\{?\\?rm basis|paper11|paper_11|loutey_paper11|Bates|'
    r'0\.6026|spectral_convergence|Molecular Fock"'
)

# Loci that MUST now be seen (regression anchors for this widening).
MUST_SEE = [
    ("papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex", 492),
    ("docs/paper_notes_archive.md", 108),
    ("docs/project_closeout_plan.md", 26),
    ("docs/development_frontier_archive.md", 190),
    ("docs/qa/synthesis.done.md", 185),
]

# Loci that MUST stay silent (innocent 0.0002 strings).
MUST_NOT_SEE = [
    ("papers/group4_quantum_computing/paper_14_qubit_encoding.tex", 916),
    ("docs/validation_benchmarks.md", 30),
]


def main() -> None:
    s = io.open(PATH, encoding="utf-8").read()
    i = s.find('"id": "p11-h2plus-0002pct-retired"')
    if i < 0:
        raise SystemExit("p11 family not found")
    j = s.find("\n    },\n", i)
    block = s[i:j]
    m = re.search(
        r'"require_nearby":\s*(r?"(?:[^"\\]|\\.)*"(?:\s*\n\s*r?"(?:[^"\\]|\\.)*")*)',
        block)
    if not m:
        raise SystemExit("require_nearby not found")
    s = s[:i] + block[:m.start(1)] + WIDER_REQUIRE + block[m.end(1):] + s[j:]
    io.open(PATH, "w", encoding="utf-8", newline="").write(s)

    import importlib.util
    spec = importlib.util.spec_from_file_location("chn", PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    e = [x for x in mod.REGISTRY if x["id"] == "p11-h2plus-0002pct-retired"][0]
    pat, req = re.compile(e["pattern"]), re.compile(e["require_nearby"])
    assert "\\\\." not in e["require_nearby"], "DOUBLED ESCAPE in require_nearby"
    win = getattr(mod, "WINDOW", 3)

    def seen(path: str, lineno: int) -> bool:
        lines = io.open(path, encoding="utf-8", errors="replace").read().splitlines()
        k = lineno - 1
        if not pat.search(lines[k]):
            return False
        ctx = "\n".join(lines[max(0, k - win):k + win + 1])
        return bool(req.search(ctx))

    print("MUST SEE (were missed before):")
    ok = True
    for path, ln in MUST_SEE:
        try:
            v = seen(path, ln)
        except Exception as exc:
            v = f"ERR {type(exc).__name__}"
        print(f"  {path}:{ln}  -> {v}")
        ok = ok and (v is True)
    print("MUST STAY SILENT (innocent 0.0002):")
    for path, ln in MUST_NOT_SEE:
        try:
            v = seen(path, ln)
        except Exception as exc:
            v = f"ERR {type(exc).__name__}"
        print(f"  {path}:{ln}  -> {v}")
        ok = ok and (v is False)
    print("WIDENING", "VERIFIED both directions" if ok else "FAILED -- do not sweep")


if __name__ == "__main__":
    main()
