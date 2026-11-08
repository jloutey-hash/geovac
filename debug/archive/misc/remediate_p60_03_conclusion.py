"""REMEDIATION 3/5 -- the conclusion + Acknowledgments sweep.

/qa paper_60 FULL 2026-09-12, same root cause as remediation 2.

  A. pre-breach molecular verdict            (synthesis M1, claims-molecular F3,
     code-A MATERIAL-3; LARGE) -- "confining the lever to homonuclear-diatomic-
     like systems" and "a frontier with a well-conditioned ground-state corner".
  B. "perfectly conditioned"                 (claims-atomic M7) -- cond -> 2.555,
     which is well-conditioned, not perfect; tab:resource gives kappa 2.3.
  C. ladder reductions bound to the wrong K  (synthesis M3, claims-atomic M4) --
     the 4.3x/3.0x/2.6x rungs are measured at K=105 s-only; adjacency binds them
     to the K=452 spdf pair.  The abstract states this correctly; the conclusion
     does not.
  D. "four results we claim as our own" is an incomplete enumeration offered as
     complete (claims-molecular F10, claims-atomic, code-A MATERIAL-3) -- the
     DoD's own lesson: "an enumeration offered as complete is a stronger claim
     than the literals it lists".

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "removable on the conditioning axis"

EDITS = [
    # ---- C: name the sector and the basis point ----
    ("""$1.72$~mHa for $2\\,^{1}S$ at identical encoding cost, and falling by $4.3\\times$,
$3.0\\times$ and $2.6\\times$ across the first three rungs of the $^{1}S$ ladder.""",
     """$1.72$~mHa for $2\\,^{1}S$ at identical encoding cost, and---on the $s$-only
ladder at $K=105$---falling by $4.3\\times$, $3.0\\times$ and $2.6\\times$ across
the first three rungs of the $^{1}S$ ladder."""),

    # ---- A + B: the molecular verdict ----
    ("""The one qualification that lifts this from a caveat to a payoff
is symmetry: because the $\\mathrm{H}_2^+$ ground state is gerade and the gerade sector
is perfectly conditioned and flat in basis size, the ground-state block-encoding pays a
metric penalty that does not grow (Table~\\ref{tab:resource})---the growth is confined to
the excited/ungerade states---though this is an \\emph{equivalent}-center property: a probe on water ($C_{2v}$) shows a symmetry-unique heavy center reinstates the growth in the ground-state ($A_1$) block (Sec.~\\ref{sec:molecular}), confining the lever to homonuclear-diatomic-like systems.""",
     """Two qualifications lift this from a caveat to a payoff.  The first is
symmetry:\\ because the $\\mathrm{H}_2^+$ ground state is gerade and the gerade
sector's conditioning is \\emph{flat} in basis size (at the constant
$2.555041\\ldots$, so flat rather than absent---Table~\\ref{tab:resource} still
charges it), the ground-state block-encoding pays a metric penalty that does not
grow, and the growth is confined to the excited/ungerade states.  That first
lever is an \\emph{equivalent}-center property, and a probe on water ($C_{2v}$)
shows a symmetry-unique heavy center reinstates the raw growth in the
ground-state ($A_1$) block (Sec.~\\ref{sec:molecular}).  The second qualification
is what removes that limit:\\ because the ill-conditioning is a symbol zero of
known order and location, a band-Toeplitz preconditioner bounds it flat in basis
size without requiring equivalent centers, and it reaches water's $A_1$ block
(Sec.~\\ref{sec:resource}).  Established for $M=2$ and non-collinear $M=3$;\\ the
collinear case is open and is not claimed."""),

    # ---- A: the closing sentence ----
    ("""The atomic
case is a result; the molecular case is a frontier with a well-conditioned ground-state corner.""",
     """The atomic case is a result;\\ the molecular case is a frontier that is
\\emph{removable on the conditioning axis} within that scope, capped on locality
by the symbol's other pole, and structurally closed on $\\ell$-sparsity."""),

    # ---- D: the enumeration ----
    ("""and four results we claim as our own:""",
     """and the results we claim as our own, which are (the list is enumerated,
not exhaustive of the paper's smaller observations):"""),
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
            print(f"    {old[:90]!r}")
            return 2
    for old, new in EDITS:
        t = t.replace(old, new)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print(f"applied: {len(EDITS)} conclusion/acknowledgments edits")
    return 0


if __name__ == "__main__":
    sys.exit(main())
