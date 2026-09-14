"""group2 Batch-1 remediation: the MATERIAL reference value + the prose NITs.

P12-KW (MATERIAL, citations) -- Kolos-Wolniewicz exact H2 is -1.174475 Ha /
   D_e 0.174475 Ha; the paper printed -1.17475 / 0.17475 (dropped digit) at
   three loci. The D_e% column already used the correct denominator (the
   numerical column 79.6/80.1 matches 0.174475, not 0.17475), so ONLY the
   printed literal is fixed; the column and headline are the paper's own
   computed values and are left untouched.
P12-92.4 (NIT, 3 reviewers) -- abstract says "92.4% with 27 basis functions";
   Table I gives 92.2% at N=27 (92.4% is the N=72 plateau). Corrected to 92.2%;
   "plateaus at 92.4%" already follows and is correct.
P8-fiedler (NIT, claims) -- "The GeoVac graph gives the same object" overclaims
   identity; L541 already calls it a "combinatorial analog". Softened to "a
   discrete analog of this object".
P11-quad (NIT, code) -- abstract "are computed algebraically ... eliminating
   quadrature entirely" implies the algebraic path is the default; it is opt-in
   (default = Gauss-Laguerre quadrature, agreeing to 1e-13). Softened to "can be
   computed", matching the body's "can alternatively".

Write-first; LaTeX raw strings; idempotent.
"""
from __future__ import annotations
import sys

P8 = "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex"
P11 = "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex"
P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"

EDITS = [
    (P12, "kw-Eexact", r"E_{\mathrm{exact}} = -1.174475",
     r"$E_{\mathrm{exact}} = -1.17475$~Ha~\cite{Kolos1968}",
     r"$E_{\mathrm{exact}} = -1.174475$~Ha~\cite{Kolos1968}"),
    (P12, "kw-De-L571", r"$D_e = 0.174475$~Ha relative",
     r"$D_e = 0.17475$~Ha relative",
     r"$D_e = 0.174475$~Ha relative"),
    (P12, "kw-De-L584", r"$D_e = 0.174475$~Ha recovered",
     r"$D_e = 0.17475$~Ha recovered",
     r"$D_e = 0.174475$~Ha recovered"),

    (P12, "abstract-92.2", r"yields 92.2\% of",
     r"""demonstrates that the Neumann algebraic $V_{ee}$ yields 92.4\% of
the exact $D_e$ with 27 basis functions, compared to 79.6\% for""",
     r"""demonstrates that the Neumann algebraic $V_{ee}$ yields 92.2\% of
the exact $D_e$ with 27 basis functions, compared to 79.6\% for"""),

    (P8, "fiedler-analog", r"gives a discrete analog of this object",
     r"The GeoVac graph gives the same object",
     r"The GeoVac graph gives a discrete analog of this object"),

    (P11, "quad-hedge", r"can be computed algebraically from three-term",
     r"are computed algebraically from three-term recurrence relations,\neliminating numerical quadrature entirely.",
     r"can be computed algebraically from three-term recurrence relations,\neliminating numerical quadrature entirely (the default path uses Gauss--Laguerre\nquadrature, which agrees to machine precision)."),
]


def main() -> int:
    loaded: dict[str, str] = {}
    applied, skipped, missed = [], [], []
    for path, name, marker, old, new in EDITS:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            skipped.append(name); continue
        if t.count(old) != 1:
            missed.append((name, t.count(old))); continue
        loaded[path] = t.replace(old, new); applied.append(name)
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied: print(f"  ok    {n}")
    for n in skipped: print(f"  skip  {n} (already applied)")
    for n, c in missed: print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
