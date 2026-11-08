"""group2 Batch-1 remediation: propagate the CORRECTED Sturmian Structural
Theorem scope (v4.72.x, PI-approved) to the stale summary surfaces of Papers 8
and 11. PI-confirmed 2026-09-13.

The no-go is for the SINGLE-n / shared-p0 / block-diagonal APPROXIMATION only;
the genuine cross-n Coulomb-Sturmian overlap couples n, escapes the theorem, and
binds (supplying the R-dependence algebraically via Shibuya-Wulfman). The body
Remark rem:overlap_imposed and Corollary cor:binding already say this and it is
already in paper_fci_molecules L513-523 (the approved template mirrored here);
only the summary surfaces were never updated -- and P8's conclusion cited
cor:binding while asserting exactly what cor:binding supersedes.

Two independent reviewers (claims chunk + code Paper 8) flagged this LARGE.

Write-first; LaTeX raw strings; idempotent.
"""
from __future__ import annotations
import sys

P8 = "papers/group2_quantum_chemistry/Paper_8_Bond_Sphere_Sturmian.tex"
P11 = "papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex"

EDITS = [
    # --- Paper 8 abstract -------------------------------------------------
    (P8, "p8-abstract", "not Sturmian molecular binding",
     r"""Binding then requires the $R$-dependent prolate
spheroidal spectrum $\beta_k(R)$, reintroducing continuous geometry.
The resolution---discretizing the Laplace--Beltrami operator in prolate
spheroidal coordinates---is developed in the companion
Paper~11~\cite{loutey_paper11}.""",
     r"""In this single-$n$, shared-$p_0$ form, binding would require the
$R$-dependent prolate spheroidal spectrum $\beta_k(R)$;\ but the genuine
cross-$n$ Coulomb--Sturmian overlap couples $n$, escapes the theorem, and
supplies that $R$-dependence algebraically via Shibuya--Wulfman
(Remark~\ref{rem:overlap_imposed}), so what is ruled out is this approximation,
not Sturmian molecular binding.  The continuous-geometry route---discretizing
the Laplace--Beltrami operator in prolate spheroidal coordinates---is developed
in the companion Paper~11~\cite{loutey_paper11}."""),

    # --- Paper 8 section-intro -------------------------------------------
    (P8, "p8-structural", "in this single-$n$ approximation without solving",
     r"""cannot produce
geometry-dependent binding in the shared-$p_0$ framework without solving a
continuous PDE at every~$R$.""",
     r"""cannot produce
geometry-dependent binding in this single-$n$ approximation without solving a
continuous PDE at every~$R$---while the genuine cross-$n$ Coulomb--Sturmian
overlap escapes the theorem and binds (Remark~\ref{rem:overlap_imposed}), so the
obstruction is to the approximation, not to Sturmian molecular binding."""),

    # --- Paper 8 conclusion (the self-contradiction) --------------------
    (P8, "p8-conclusion", "the genuine\ncross-$n$ Coulomb--Sturmian overlap escapes the theorem and supplies it",
     r"""Binding requires the full
$R$-dependent prolate spheroidal eigenvalue spectrum $\beta_k(R)$, which
reintroduces continuous spatial geometry (Corollary~\ref{cor:binding}).""",
     r"""In this single-$n$, shared-$p_0$ form, binding would require the full
$R$-dependent prolate spheroidal eigenvalue spectrum $\beta_k(R)$;\ the genuine
cross-$n$ Coulomb--Sturmian overlap escapes the theorem and supplies it
algebraically (Corollary~\ref{cor:binding}), so what is ruled out is this
approximation, not Sturmian molecular binding."""),

    # --- Paper 11 abstract ----------------------------------------------
    (P11, "p11-abstract", "in that single-$n$ shared-$p_0$ approximation",
     r"""the single-$p_0$ momentum scale cannot encode the
$R$-dependent bonding physics (Paper~8--9, negative theorem).""",
     r"""the single-$p_0$ momentum scale cannot encode the
$R$-dependent bonding physics in that single-$n$ shared-$p_0$ approximation
(Paper~8--9, negative theorem;\ the genuine cross-$n$ Coulomb--Sturmian basis
escapes it and binds)."""),

    # --- Paper 11 intro -------------------------------------------------
    (P11, "p11-intro", "one route to it, while the genuine cross-$n$",
     r"""Binding in this
framework requires the full $R$-dependent eigenvalue spectrum
$\beta_k(R)$ from the continuous prolate spheroidal PDE---reintroducing
the very continuous spatial geometry the graph was meant to replace.""",
     r"""Binding in this
single-$n$ shared-$p_0$ framework would require the full $R$-dependent
eigenvalue spectrum $\beta_k(R)$;\ the continuous prolate spheroidal PDE below is
one route to it, while the genuine cross-$n$ Coulomb--Sturmian overlap supplies
it algebraically (Paper~8, Remark), so the obstruction is to this approximation,
not to Sturmian molecular binding."""),

    # --- Paper 11 sec 8.2 -----------------------------------------------
    (P11, "p11-sec82", "cannot produce\n$R$-dependent binding \\emph{in the single-$n$, block-diagonal approximation}:",
     r"""Paper~8--9 proved that the molecular Sturmian basis on a single
$S^3$ with shared momentum scale $p_0$ cannot produce
$R$-dependent binding:""",
     r"""Paper~8--9 proved that the molecular Sturmian basis on a single
$S^3$ with shared momentum scale $p_0$ cannot produce
$R$-dependent binding \emph{in the single-$n$, block-diagonal approximation}:"""),

    # --- Paper 11 conclusion --------------------------------------------
    (P11, "p11-conclusion", "cannot\nproduce $R$-dependent binding \\emph{in its single-$n$ form}",
     r"""the single-$p_0$ Sturmian basis on $S^3$ cannot
produce $R$-dependent binding, but the prolate spheroidal coordinate""",
     r"""the single-$p_0$ Sturmian basis on $S^3$ cannot
produce $R$-dependent binding \emph{in its single-$n$ form} (the genuine
cross-$n$ basis escapes this;\ Paper~8, Remark), while the prolate spheroidal coordinate"""),
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
