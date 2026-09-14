"""group2 Batch-2 remediation, papers: the P17 contradiction, prose NITs, and
the two citation defects that must not ship.

P17-VIII.B (contradiction) -- L1607 "molecular equilibrium exists without PK ...
   a structural result" + "5x overbinding" is the ADIABATIC reading; the
   authoritative 2D-variational result (conclusion table + CLAUDE §5) is UNBOUND
   (D_e<0). Reconcile: the minimum is a feature of the adiabatic PES, not a bound
   molecule.
P15-abstract-cusp (F2) -- the 96.0% includes a ~1pp Schwartz cusp correction the
   abstract omits (body says "2D+cusp"); state it, and that the pure-variational
   value ~95% already exceeds 92.4%.
P15-kolos (truncation) -- exact H2 D_e printed 0.17447 (truncation of 0.174475);
   correct to 0.174475 (same class as Batch 1's KW value).
P17-category (F5) -- "first realization of a category ... morphisms" is an
   unverified category-theory novelty claim; the paper already has the correct
   fiber-bundle framing. Drop the category claim.
Mitnik (WRONG-ID) -- cited to CPC 269/108145 (a QCD code, REvolver); real paper
   verified at source: D. M. Mitnik, F. A. Lopez, L. U. Ancarani, "Generalized
   Sturmian Functions in prolate spheroidal coordinates," Mol. Phys. 119(8),
   e1881179 (2021) [arXiv:2006.06616]. Third author Gasaneo -> Lopez.
Madden (MISATTRIBUTION) -- the 3.29 doubly-excited ^1S width ratio is attributed
   to Madden-Codling, who measured the dipole-allowed ^1P^o series; reframe.

Write-first; LaTeX raw strings; idempotent.
"""
from __future__ import annotations
import sys

P13 = "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex"
P15 = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"
P17 = "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex"

EDITS = [
    (P17, "p17-viiib", "the minimum is a feature of the adiabatic potential curve",
     r"""This is a \emph{structural} result: molecular equilibrium exists
without PK.  The 64\% $R_{\rm eq}$ error and $5{\times}$ overbinding
indicate that $l_{\max} = 2$ is far from convergence for the""",
     r"""This is a \emph{structural} result: an adiabatic PES minimum exists
without PK.  (The $l_{\max}=2$ 2D-variational treatment is itself \emph{unbound},
$D_e<0$, as the conclusion table records;\ the minimum is a feature of the adiabatic potential curve, not a bound molecule.)  The 64\% $R_{\rm eq}$ error
and $5{\times}$ \emph{adiabatic} overbinding
indicate that $l_{\max} = 2$ is far from convergence for the"""),

    (P15, "p15-abstract-cusp", "and a Schwartz cusp correction",
     r"""recovers 96.0\% of the exact dissociation energy $D_e$ at
$l_{\max}=6$ with $\sigma$ and $\pi$ channels (61~channels,
CBS extrapolation ${\sim}97\%$),""",
     r"""recovers 96.0\% of the exact dissociation energy $D_e$ at
$l_{\max}=6$ with $\sigma$ and $\pi$ channels and a Schwartz cusp correction
(61~channels, CBS extrapolation ${\sim}97\%$;\ the pure-variational value is
${\sim}95\%$, itself above 92.4\%),"""),

    (P15, "p15-kolos", r"D_e^{\rm exact} = 0.174475",
     r"$D_e^{\rm exact} = 0.17447$~Ha~\cite{kolos1968}.}",
     r"$D_e^{\rm exact} = 0.174475$~Ha~\cite{kolos1968}.}"),

    (P17, "p17-cat-L111", "a fiber-bundle composition of natural",
     r"""This is the first realization of a \emph{category} of natural
geometries with morphisms preserving lattice structure.""",
     r"""This is a fiber-bundle composition of natural
geometries with typed inter-layer maps ($R\to\infty$ dissociation,
$Z\to Z'$ isoelectronic substitution)."""),

    (P17, "p17-cat-L618", "The layers are natural",
     r"""This is a mathematical structure (a category) where objects are natural
geometries and morphisms are the continuous deformations""",
     r"""The layers are natural
geometries and the typed inter-layer maps are the continuous deformations"""),

    (P13, "mitnik-wrongid", "e1881179",
     r"""D.~M.~Mitnik, L.~U.~Ancarani, and G.~Gasaneo,""",
     r"""D.~M.~Mitnik, F.~A.~L\'opez, and L.~U.~Ancarani,"""),

    (P13, "mitnik-venue", r"Mol.\ Phys.\ \textbf{119}",
     r"Comput.\ Phys.\ Commun.\ \textbf{269}, 108145 (2021).",
     r"``Generalized Sturmian Functions in prolate spheroidal coordinates,'' \textit{Mol.\ Phys.}\ \textbf{119}(8), e1881179 (2021);\ arXiv:2006.06616."),

    (P13, "madden-misattr", "the dipole-allowed ${}^1\\!P^{\\rm o}$ series",
     r"(experiment~\cite{Madden1963})---probe physics beyond the adiabatic",
     r"(the autoionization phenomenon was observed by Madden--Codling~\cite{Madden1963}, there in the dipole-allowed ${}^1\!P^{\rm o}$ series;\ the doubly-excited ${}^1\!S$ ratio quoted here is illustrative)---probe physics beyond the adiabatic"),
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
