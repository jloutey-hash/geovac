r"""Sweep the retired H2+ "0.0002%" and He "0.019%" headlines, 40 loci.

PI direction 2026-09-19: state the H2+ accuracy QUALITATIVELY ("machine
precision"); do NOT substitute another percentage, because the residual is
reference-limited and a mantissa reinvites the same defect in smaller form.

MEASURED basis for every replacement below (production, this session):
  ProlateSpheroidalLattice(R=2.0, radial_method='spectral', n_basis=20)
    -> E_total = -0.602634214494936 Ha
       |err| vs E_ref = -0.6026342144949 Ha  ->  3.6e-14 Ha  (6.0e-12 %)
       n_basis=5 -> 2.13e-8 Ha (3.5e-6 %);  n_basis=10 -> 2.53e-12 Ha
  fit_spectroscopic_constants on a 0.005 grid -> R_eq = 1.99726 bohr
       (+0.013 % vs 1.997), E_min = -0.60263464, D_e = 0.10263464
  Retired: 0.0002 % (= 1.21e-6 Ha) and R_eq 2.005 / 0.38 %.

FOUR DERIVED FIGURES DIE WITH THE HEADLINE, and this script retires them too:
  * "5000x accuracy improvement" -- computed as 1.01%/0.0002%.  With the
    numerator a discretization artifact and the denominator retired, no
    defensible ratio exists, so the claim is dropped rather than recomputed.
  * R_eq 2.005 bohr and its 0.38 % -- coarse-grid FIT artifacts (MATERIAL-2).
  * "0.70 % with the original finite-difference lattice" in Papers 15 and 17,
    which contradicts Paper 11's own table (1.01 %).  Pre-existing, swept here
    because leaving it would replace one inconsistency with another.

Written as a FILE, not a heredoc: every replacement carries LaTeX backslashes,
and two C17 regexes shipped dead today from heredoc escape mangling.

Run:  python debug/qa/_sweep_p11_h2plus_headline.py
"""
from __future__ import annotations

import io
import sys

MACH = "machine precision"

EDITS: list[tuple[str, str, str]] = []


def add(path: str, old: str, new: str) -> None:
    EDITS.append((path, old, new))


P11 = "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex"
P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
P13 = "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex"
P15 = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"
P17 = "papers/group2_quantum_chemistry/paper_17_composed_geometries.tex"
P8 = "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex"
PFM = "papers/group2_quantum_chemistry/paper_fci_molecules.tex"
SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
FG = "papers/synthesis/geovac_field_guide.tex"
P34 = "papers/group6_precision_observations/paper_34_projection_taxonomy.tex"

# ---------------------------------------------------------------- Paper 11
# L66-70: the abstract-adjacent result sentence.  Both numbers wrong, and the
# "-0.6026 vs -0.6026" form made the claim unfalsifiable as printed.
add(P11,
    r"""$E_{\min} = -0.6026$~Ha (0.0002\% error vs.\ exact $-0.6026$~Ha),
equilibrium bond length $R_{\mathrm{eq}} = 2.005$~bohr
(0.38\% error), and correct dissociation to H~+~p at $R \to \infty$.
The spectral solver achieves a $250\times$ dimension reduction,
$270\times$ speedup, and $5000\times$ accuracy improvement over""",
    r"""$E_{\min} = -0.602634214494936$~Ha, reproducing the two-center
reference $-0.6026342144949$~Ha~\cite{bates1953} to
$3.6\times10^{-14}$~Ha---at or below the precision to which that
reference is conventionally quoted, so the residual is
reference-limited rather than a property of the method
(\textbf{[MEASURED 2026-09-19]}).  Equilibrium bond length
$R_{\mathrm{eq}} = 1.9973$~bohr ($0.013\%$ error), and correct
dissociation to H~+~p at $R \to \infty$.
The spectral solver achieves a $250\times$ dimension reduction and a
$270\times$ speedup over""")

add(P11,
    r"radial solver that achieves $0.0002\%$ accuracy.",
    r"radial solver that reaches " + MACH + r".")

add(P11,
    r"bottleneck, achieving $0.0002\%$ accuracy with $N_b = 20$.",
    r"bottleneck, reaching " + MACH + r" with $N_b = 20$.")

add(P11,
    r"""achieving $0.0002\%$ error with $N_b = 20$ basis functions---a
$5000\times$ accuracy improvement over the FD solver at $250\times$
smaller matrix dimension.""",
    r"""reaching """ + MACH + r""" with $N_b = 20$ basis
functions, at $250\times$ smaller matrix dimension than the FD solver.""")

add(P11,
    r"""energy to $0.0002\%$ accuracy and the equilibrium bond length to
$0.38\%$ accuracy with zero free parameters.""",
    r"""energy to """ + MACH + r""" and the equilibrium bond length to
$0.013\%$ accuracy with zero free parameters.""")

add(P11,
    r"""wavefunction's asymptotic decay, achieving $0.0002\%$ accuracy
with $N_b = 20$.""",
    r"""wavefunction's asymptotic decay, reaching """ + MACH + r"""
with $N_b = 20$.""")

# ---------------------------------------------------------------- Paper 12
add(P12,
    r"""geometry for diatomic molecules, solving H$_2^+$ to $0.0002\%$ energy
error (spectral Laguerre solver) with zero free parameters.""",
    r"""geometry for diatomic molecules, solving H$_2^+$ to """ + MACH + r"""
(spectral Laguerre solver) with zero free parameters.""")

add(P12,
    r"""molecule H$_2^+$, achieving $0.0002\%$ energy error (spectral Laguerre
solver) with zero free parameters.""",
    r"""molecule H$_2^+$, reaching """ + MACH + r""" (spectral Laguerre
solver) with zero free parameters.""")

add(P12,
    r"    Accuracy: 0.0002\% (spectral Laguerre)~\cite{loutey_paper11}.",
    r"    Accuracy: " + MACH + r" (spectral Laguerre)~\cite{loutey_paper11}.")

# ---------------------------------------------------------------- Paper 13
add(P13,
    r"""    prolate spheroidal lattice achieves 0.0002\% error for H$_2^+$
    with a spectral Laguerre radial solver and zero free parameters.""",
    r"""    prolate spheroidal lattice reaches """ + MACH + r""" for H$_2^+$
    with a spectral Laguerre radial solver and zero free parameters.""")

# ---------------------------------------------------------------- Paper 15
add(P15,
    r"""  Paper~11~\cite{paper11} achieves $0.0002\%$ energy error for H$_2^+$
  with a spectral Laguerre radial solver (0.70\% with the original
  finite-difference lattice).""",
    r"""  Paper~11~\cite{paper11} reaches """ + MACH + r""" for H$_2^+$
  with a spectral Laguerre radial solver ($1.01\%$ with the original
  finite-difference lattice).""")

add(P15,
    r"2 & H$_2^+$ & Prolate spheroid    & 0.0002\% error (spectral) \\",
    r"2 & H$_2^+$ & Prolate spheroid    & " + MACH + r" (spectral) \\")

# ---------------------------------------------------------------- Paper 17
add(P17,
    r"""  exactly.  Paper~11~\cite{paper11} achieves $0.0002\%$ energy error
  (spectral Laguerre radial solver; 0.70\% with the original
  finite-difference lattice) and zero free parameters.""",
    r"""  exactly.  Paper~11~\cite{paper11} reaches """ + MACH + r"""
  (spectral Laguerre radial solver; $1.01\%$ with the original
  finite-difference lattice) and zero free parameters.""")

add(P17,
    r"2 & H$_2^+$ & Prolate spheroid       & 0.0002\% error (spectral) \\",
    r"2 & H$_2^+$ & Prolate spheroid       & " + MACH + r" (spectral) \\")

# ------------------------------------------------- Paper 8 / FCI-molecules
for _p in (P8, PFM):
    add(_p,
        r"""reaches $0.0002\%$ energy error, both with zero fitted
parameters~\cite{loutey_paper11}.""",
        r"""reaches """ + MACH + r""", both with zero fitted
parameters~\cite{loutey_paper11}.""")

# ------------------------------------------------------------ synthesis
add(SYN,
    r"""$0.0002\%$ for $\mathrm{H}_2^+$, $0.022\%$ (properly variational, raw""",
    MACH + r""" for $\mathrm{H}_2^+$, $0.022\%$ (properly variational, raw""")

add(SYN,
    r"2 & $\mathrm{H}_2^+$ (2, 1) & prolate spheroid & $0.0002\%$ energy (spectral Laguerre) & \cite{loutey_paper11} \\",
    r"2 & $\mathrm{H}_2^+$ (2, 1) & prolate spheroid & " + MACH + r" (spectral Laguerre) & \cite{loutey_paper11} \\")

add(SYN,
    r"""$E_{\min} = -0.6026$~Ha against the exact $-0.6026$~Ha~\cite{bates1953}---a $0.0002\%$
error, the most accurate single result in the chemistry
arc---with the correct dissociation to $\mathrm{H}+p$ as""",
    r"""$E_{\min} = -0.602634214494936$~Ha against the reference
$-0.6026342144949$~Ha~\cite{bates1953}---agreement to
$3.6\times10^{-14}$~Ha, i.e.\ to the precision at which that reference
is quoted, making it the most accurate single result in the chemistry
arc---with the correct dissociation to $\mathrm{H}+p$ as""")

add(SYN,
    r"for $\mathrm{H}_2^+$ ($0.0002\%$)~\cite{loutey_paper11}, hyperspherical",
    r"for $\mathrm{H}_2^+$ (" + MACH + r")~\cite{loutey_paper11}, hyperspherical")

# ------------------------------------------------------------ field guide
add(FG,
    r"2 & H$_2^+$ (2-center, 1e) & Prolate spheroid & $0.0002\%$ & 11 \\",
    r"2 & H$_2^+$ (2-center, 1e) & Prolate spheroid & " + MACH + r" & 11 \\")

# ------------------------------------------------------------ Paper 34 (He)
add(P34,
    r"""$0.019\%$ error (Track~DI v2.6.0, Sprint~1; 2D variational solver on""",
    r"""$0.022\%$ raw error (properly variational; $0.004\%$ after a
non-variational cusp extrapolation at $l_\text{max}=4$) (Track~DI
v2.6.0, Sprint~1; 2D variational solver on""")

add(P34,
    r"""$0.019\%$ at $\ell_\text{max}=7$ & T & Track DI Sprint 1; cusp-corrected
$0.004\%$ at $\ell_\text{max}=4$. \\""",
    r"""$0.022\%$ at $\ell_\text{max}=7$ & T & Track DI Sprint 1; cusp-corrected
$0.004\%$ at $\ell_\text{max}=4$. \\""")


def main() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    by_file: dict[str, list[tuple[str, str]]] = {}
    for path, old, new in EDITS:
        by_file.setdefault(path, []).append((old, new))

    failures: list[str] = []
    for path, pairs in by_file.items():
        s = io.open(path, encoding="utf-8").read()
        n = 0
        for old, new in pairs:
            if old not in s:
                failures.append(f"{path}: anchor NOT FOUND -> {old[:70]!r}")
                continue
            if s.count(old) > 1:
                failures.append(f"{path}: anchor AMBIGUOUS ({s.count(old)}x) -> {old[:60]!r}")
                continue
            s = s.replace(old, new, 1)
            n += 1
        if n:
            io.open(path, "w", encoding="utf-8", newline="").write(s)
        print(f"  {path}: {n}/{len(pairs)} applied")

    print()
    if failures:
        print("FAILURES (nothing was written for these):")
        for f in failures:
            print("   " + f)
        raise SystemExit(1)
    print(f"all {len(EDITS)} anchored replacements applied cleanly")


if __name__ == "__main__":
    main()
