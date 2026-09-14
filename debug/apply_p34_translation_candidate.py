"""Log the inter-centre translation as a Paper 34 §VIII projection candidate.

NOT a §III entry.  The tag-transcendentals STOP rule (memory/feedback_tag_
transcendentals.md, CLAUDE.md §4) says a candidate projection that does not fit
the existing inventory is flagged for §VIII open-question review and is NOT
silently written into §III without PI direction.  This does the flagging half.

Verified before writing: the string "translat" appears NOWHERE in Paper 34, and
the only Shibuya-Wulfman mentions (§III.22 multipole expansion) name the SW
BASIS expansion as a contrast -- something that does NOT terminate -- never as a
projection.  §III.11 covers Wigner-D ROTATION between centres and explicitly
preserves rationality; nothing covers TRANSLATION.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group6_precision_observations/paper_34_projection_taxonomy.tex"
MARKER = "inter-center translation on the Fock sphere"

ANCHOR = """Decision pending on (a), (b); structural progress on (c), (d) named
in their respective entries below.
"""

NEW = r"""Plus one candidate surfaced by the Paper~60 molecular arc (2026-09-12):
(e)~\textbf{inter-center translation on the Fock sphere.}  The two-center
Shibuya--Wulfman operator is multiplication by the translation phase
$e^{i\vec{p}\cdot\vec{R}}$ on the Fock sphere, and its zonal average
$j_0(pR)$ is the \emph{Euclidean} $\mathrm{E}(3)$ zonal spherical function
rather than the $\mathrm{SO}(4)$ Gegenbauer one.  \S~\ref{sec:proj_wignerD}
covers the Wigner-$D$ \emph{rotation} between molecular centers and
explicitly preserves rationality up to
$\mathbb{Q}[\sqrt2,\sqrt3,\sqrt6]$;\ no entry covers the
\emph{translation}, which is where the Bessel content enters.  Its
transcendental signature is calibration-tier and \emph{two-sided}, which is
what makes it interesting as a slot rather than a relabel:\ the
finite-section truncation prices $\pi^2$ (M2, the $\pi^2\cdot\mathbb{Q}$
half --- the Kac--Murdock--Szeg\H{o} constant, obtainable with no Bessel
function present at all), while the symbol's asymptotics at the opposite
pole price $(2\pi)^{-1/2}$ and a $\pi/4$ branch phase (M2, the
$\sqrt\pi\cdot\mathbb{Q}$ half).  Paper~60 \S~molecular;\ CHANGELOG
v5.11.4.  \emph{Caution if promoted:}\ the reading of $j_0$ as the
contraction limit of the $S^3$ zonal function is \textbf{prior art} ---
Clerc, \textit{Studia Math.}\ \textbf{57}, 27 (1976);\ D\'iaz Mart\'in and
Pacharoni, arXiv:1807.03904;\ lineage from In\"on\"u and Wigner,
\textit{PNAS} \textbf{39}, 510 (1953) --- so the entry would cite that
reading, never claim it, and would have to state the antipodal parity
caveat ($n^{-1}U_{n-1}(\cos(\pi-z/n))\to(-1)^{n+1}j_0(z)$, convergence
along parities only).
Decision pending on (a), (b), (e);\ structural progress on (c), (d) named
in their respective entries below.
"""


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(ANCHOR) != 1:
        print(f"anchor count={t.count(ANCHOR)}; aborting")
        return 2
    t = t.replace(ANCHOR, NEW)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("applied: candidate (e) logged in Paper 34 §VIII")
    return 0


if __name__ == "__main__":
    sys.exit(main())
