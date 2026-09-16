"""Register the withdrawn Paper-12 cusp diagnosis in check_retracted_terms.py,
and sweep the remaining documentation surfaces.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


GATE = "debug/qa/check_retracted_terms.py"
BENCH = "docs/validation_benchmarks.md"
MATRIX = "docs/claim_test_matrix.md"
FIGREADME = "papers/group2_quantum_chemistry/paper_15_figures/README.md"

# ------------------------------------------------- new registry entry
edit(
    GATE,
    r"""    {
        "id": "sp-splitting-as-convergence-evidence",""",
    r'''    {
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
        "cited_by": {
            "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex": "reviewed 2026-09-14",
            "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex": "reviewed 2026-09-14",
            "papers/synthesis/group2_quantum_chemistry_synthesis.tex": "reviewed 2026-09-14",
            "docs/validation_benchmarks.md": "reviewed 2026-09-14",
            "docs/claim_test_matrix.md": "reviewed 2026-09-14",
            "tests/test_level4_multichannel.py": "reviewed 2026-09-14",
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
        "id": "sp-splitting-as-convergence-evidence",''',
    "registry: p12-cusp-as-the-h2-gap entry")

# ------------------------------------------------- validation benchmarks
edit(
    BENCH,
    r"""| H2 Neumann V_ee | 92.4% D_e | Algebraic integral accuracy |""",
    r"""| H2 Neumann V_ee (sigma only) | 92.4% D_e | Algebraic integral accuracy |
| H2 Neumann V_ee, \|m\| <= 1 | 99.09% D_e | Azimuthal channels restored (Paper 12, 2026-09-14) |""",
    "benchmarks: scope the 92.4% + add the extended value")

# ------------------------------------------------- figure README
edit(
    FIGREADME,
    r"""  - Horizontal dashed line at 92.4% (Paper 12 Neumann V_ee)""",
    r"""  - Horizontal dashed line at 92.4% (Paper 12 Neumann V_ee, SIGMA-ONLY --
    not a like-for-like reference line for a sigma+pi curve; Paper 12's own
    basis with |m| <= 1 reaches 99.09%)""",
    "figure README: scope the reference line")

with io.open(MATRIX, encoding="utf-8") as fh:
    matrix = fh.read()
old_row = "| 12 | full V_ee \"exact\" + **92.4% D_e headline** |"
if old_row in matrix:
    matrix = matrix.replace(
        old_row,
        "| 12 | full V_ee \"exact\" + **92.4% D_e headline (SIGMA-ONLY; "
        "the |m| <= 1 value is 99.09%, see Paper 12 Sec. 'Restoring the "
        "Azimuthal Channels')** |", 1)
    with io.open(MATRIX, "w", encoding="utf-8") as fh:
        fh.write(matrix)
    print("  + claim_test_matrix: P12 row scoped")
else:
    print("  ! claim_test_matrix P12 row not matched")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

failed, applied = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        text = fh.read()
    ok = True
    for old, new, label in items:
        if old in text:
            text = text.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append(label)
            ok = False
    if ok:
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write(text)

if failed:
    print("FAILED TO MATCH:")
    for f in failed:
        print("  -", f)
    sys.exit(1)

for a in applied:
    print("  +", a)
print("done")
