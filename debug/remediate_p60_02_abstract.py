"""REMEDIATION 2/5 -- the abstract sweep.

/qa paper_60 FULL 2026-09-12.  Four findings converging on one cause: every
re-tiering and scoping of v5.10.18..v5.11.4 reached the body paragraph that owns
the claim and never reached the abstract.

  A. growth law claimed as DERIVED     (claims-atomic M2/F1, synthesis M2, LARGE)
     -> body L852 carries [PRIOR ART] + "we claim only the identification".
  B. molecular verdict is PRE-BREACH   (claims-molecular F3, synthesis M1, LARGE)
     -> "essentially removed ... (a water probe shows a symmetry-inequivalent
        heavy center reinstates it)" is the reading sec:resource reverses; the
        preconditioned water A_1 column is BOUNDED (38.45 -> 44.06).
  C. "essentially removed" vs tab:resource            (claims-atomic M7)
     -> the gerade row is kappa 2.3 / d_inv 18, not kappa 1 / d_inv 0. FLAT,
        not removed.
  D. lever inventory names THREE       (claims-molecular F3, code-A MATERIAL-3)
     -> sec:resource names a fourth and calls it stronger than either of the
        first two on the conditioning axis.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "Kac--Murdock--Szeg\\H{o} asymptotic, and what is ours"

EDITS = [
    # ---- A ----
    ("""\\textbf{[SYMBOLIC + MEASURED]} The growth law is
derived (a band-limited concentration rate,
$1-\\sigma_{\\max}\\propto(kR)^2/n^2$, giving asymptotic exponent $2$---the
leading order is derived, the remainder measured rather than bounded;\\ the
fitted window exponents $1.85$/$1.97$ are pre-asymptotic readings).""",
     """\\textbf{[PRIOR ART]} The growth law is not new:\\ the band-limited
concentration rate $1-\\sigma_{\\max}\\propto(kR)^2/n^2$, asymptotic exponent
$2$, is the Kac--Murdock--Szeg\\H{o} asymptotic, and what is ours is the
\\emph{identification} of the Shibuya--Wulfman metric as such a finite section
(the fitted window exponents $1.85$/$1.97$ are pre-asymptotic readings)."""),

    # ---- B, C, D: the molecular verdict ----
    ("""\\textbf{[OPEN]} The atomic case is a clean, apparently novel,
metric-free quantum secular equation; the molecular case is a genuine
frontier in which the metric cost is mitigated for excited/ungerade states, essentially removed
for the ground state of \\emph{equivalent-center} systems (a water probe shows a symmetry-inequivalent heavy center reinstates it), and confined to one electron.""",
     """\\textbf{[MEASURED]} A third lever removes the conditioning cost outright
where it applies:\\ the ill-conditioning is a symbol zero of known order and
location, so a band-Toeplitz preconditioner bounds $\\mathrm{cond}$ flat in basis
size, and---unlike the gerade lever---it does not require equivalent centers,
reaching water's $A_1$ block (bounded at $44$ out to $N=192$, against a raw
$N^{1.97}$).  A direct block-encoding of the preconditioned metric then takes the
penalty from $n^3$ to $n$.  \\textbf{[OPEN]} The atomic case is a clean,
apparently novel, metric-free quantum secular equation;\\ the molecular case is a
frontier on which the \\emph{conditioning} axis is removable for $M=2$ and
non-collinear $M=3$, the \\emph{locality} cost is capped by the symbol's other
pole, and the $\\ell$-block sparsity loss is structural and not a conditioning
effect at all."""),

    # ---- D: the closing inventory ----
    ("""The real molecular savings are the metric levers (gerade sector, large
separation, one-electron-only metric), not a sublinear matrix.""",
     """The real molecular savings are the metric levers (the band-Toeplitz
preconditioner, the gerade sector, large separation, one-electron-only metric),
not a sublinear matrix."""),
]


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    for i, (old, _) in enumerate(EDITS, 1):
        if t.count(old) != 1:
            print(f"  edit {i} anchor count={t.count(old)}; ABORT")
            return 2
    for old, new in EDITS:
        t = t.replace(old, new)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print(f"applied: {len(EDITS)} abstract edits (prior-art tier, "
          f"post-breach verdict, four-lever inventory)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
