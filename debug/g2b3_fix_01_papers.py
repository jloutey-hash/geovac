"""group2 Batch-3 remediation, papers: the Paper 19 pair-diagonal zombie
(prose contradicting the paper's own re-measured exact-rule tables) + the
FCI-molecules M3/M4 softenings.

Paper 19 -- the re-measured tables (2026-09-01, exact global-M_L rule) carry
the CURRENT values; the prose was never updated and cites RETIRED pair-diagonal
figures that contradict the tables directly:
  * Step-3 (L594-600): "19,959 Pauli / 448.9 Ha / 2,298 QWC" and exponents
    3.03/1.75 and ratio 2.53x~2.63x  ->  tab:convergence gives 127,855 Pauli,
    436.9 Ha at n_max=3; 2726 Pauli / 3.25x at n_max=2 (tab:resources,
    tab:balanced_census).  Exponents recomputed from the table values:
      Pauli:  ln(127855/2726)/ln(84/30) = 3.74
      1-norm: ln(436.9/75.2)/ln(84/30)  = 1.71
  * Step-4 (L850): same exponents 3.03/1.75 -> 3.74/1.71.
  * Polyatomic ratio (L917): 2.63/4.77/7.45x -> 3.25/6.35/10.10x
    (= the census table's own ratio column). Growth exponent recomputed
    B^1.14 -> B^1.24 (ln(10.10/3.25)/ln(5/2)); per-block increments
    2.14/1.34 -> 3.10/1.88.
  * Census (L144, L399): 130/195 nonzero ERIs -> 214/321 (exact-rule,
    already the GREEN pin in test_cross_block_eri_count; ratio stays 2:3).

FCI-molecules:
  * M3 (L556-559): "the only sound discrete framework for heteronuclear FCI"
    -> scope the "only" to the alternatives surveyed here.
  * M4 (L773-778): "the correct weighting" -> "a qualitatively correct weighting".

Write-first; LaTeX raw strings; idempotent.
"""
from __future__ import annotations
import sys

P19 = "papers/group2_quantum_chemistry/paper_19_coupled_composition.tex"
FCM = "papers/group2_quantum_chemistry/paper_fci_molecules.tex"

EDITS = [
    # --- Paper 19 Step-3 prose (the core zombie) ---
    (P19, "p19-step3", "127{,}855 Pauli terms and a 436.9",
     r"""has 19{,}959 Pauli terms, 448.9~Ha 1-norm, and 2{,}298 QWC
groups.  The Pauli scaling exponent (2-point fit from $Q = 30$ to
$Q = 84$) is 3.03; the 1-norm exponent is 1.75.  The
Pauli-to-composed ratio (2.53$\times$) is consistent with
$n_{\max} = 2$ (2.63$\times$), indicating that the cross-block
overhead is a constant factor, not a growing penalty.""",
     r"""has 127{,}855 Pauli terms and a 436.9~Ha 1-norm
(Table~\ref{tab:convergence}).  The two-point Pauli scaling exponent
(from $Q = 30$ to $Q = 84$) is 3.74 and the 1-norm exponent is 1.71;
at $n_{\max} = 2$ the balanced Hamiltonian carries 2{,}726 Pauli terms,
$3.25\times$ the composed count of 838."""),

    # --- Paper 19 Step-4 exponents ---
    (P19, "p19-step4", "Pauli scaling exponent 3.74, 1-norm exponent 1.71.",
     r"Pauli scaling exponent 3.03, 1-norm exponent 1.75.",
     r"Pauli scaling exponent 3.74, 1-norm exponent 1.71."),

    # --- Paper 19 polyatomic ratio sequence ---
    (P19, "p19-ratio", r"$3.25\times$ (2-block LiH)",
     r"""$2.63\times$ (2-block LiH) $\to$ $4.77\times$ (3-block BeH$_2$)
$\to$ $7.45\times$ (5-block H$_2$O), because each block pair adds
cross-block ERIs and cross-center $V_{ne}$ terms.  The growth is
sub-\emph{quadratic} in the block count over these three systems
($\sim B^{1.14}$---slightly super-linear;\ the per-block increment
\emph{decelerates}, $2.14$ then $1.34$ per block), so an $O(B^2)$ envelope""",
     r"""$3.25\times$ (2-block LiH) $\to$ $6.35\times$ (3-block BeH$_2$)
$\to$ $10.10\times$ (5-block H$_2$O), because each block pair adds
cross-block ERIs and cross-center $V_{ne}$ terms.  The growth is
sub-\emph{quadratic} in the block count over these three systems
($\sim B^{1.24}$---slightly super-linear;\ the per-block increment
\emph{decelerates}, $3.10$ then $1.88$ per block), so an $O(B^2)$ envelope"""),

    # --- Paper 19 census table (nonzero ERI counts) ---
    (P19, "p19-census-tab", r"Nonzero ERIs & 321 & 214 \\",
     r"Nonzero ERIs & 195 & 130 \\",
     r"Nonzero ERIs & 321 & 214 \\"),

    # --- Paper 19 census in-prose reference ---
    (P19, "p19-census-prose", r"ERI count (214 vs.~321 for LiH",
     r"ERI count (130 vs.~195 for LiH at $n_{\max} = 2$).",
     r"ERI count (214 vs.~321 for LiH at $n_{\max} = 2$)."),

    # --- FCI-molecules M3: scope the "only" ---
    (FCM, "fcm-m3", r"remains, among the alternatives surveyed here, the",
     r"""Topological LCAO remains the only sound discrete framework for
heteronuclear FCI.""",
     r"""Topological LCAO remains, among the alternatives surveyed here, the
only sound discrete framework for heteronuclear FCI."""),

    # --- FCI-molecules M4: qualitatively correct ---
    (FCM, "fcm-m4", r"provides a qualitatively correct weighting",
     "provides\nthe correct weighting for kinetic repulsion",
     "provides\na qualitatively correct weighting for kinetic repulsion"),
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
