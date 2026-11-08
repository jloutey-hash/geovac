"""DELTA #2 remediation: the defects both reviewers found in round-1's own fixes.

All verified against primary text by the PM first.  The recurring shape is the
same one the corpus keeps recording: a claim corrected at the locus where it
was reported and left standing where it was restated.

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

P12 = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"
P13 = "papers/group2_quantum_chemistry/paper_13_hyperspherical.tex"
P15 = "papers/group2_quantum_chemistry/paper_15_level4_geometry.tex"
SYN = "papers/synthesis/group2_quantum_chemistry_synthesis.tex"
FG = "papers/synthesis/geovac_field_guide.tex"
EDITS = []


def edit(path, old, new, label):
    EDITS.append((path, old, new, label))


# ============ LARGE: the reverse-ordering claim, recreated at three loci ====
edit(P15,
     r"""At matched angular content the ordering reverses---Paper~12's own basis
with $|m| \le 1$ reaches 99.1\%---so no coordinate-system advantage is
claimed here.""",
     r"""At matched $m_{\max}$ no comparison between the two is available:\ the
two published figures are not a matched pair, so this paper claims no
coordinate-system advantage in either direction.""",
     "D1a: P15 abstract -- reverse-ordering claim removed")

edit(P15,
     r"""Run at matched angular content the ordering reverses.  Paper~12's own
basis with $|m| \le 1$ reaches 99.1\%, against 94.1\% here at
$l_{\max}=4$ and 96.0\% at $l_{\max}=6$ with a cusp correction;\ and a
grid-based prolate spheroidal calculation reaches
99.97\%~\cite{tao_mccurdy_rescigno2010}.""",
     r"""Nor is a corrected comparison available in the other direction.  The
figures differ in more than $m_{\max}$:\ Paper~12's $99.1\%$ is
$(j_{\max},l_{\max}) = (3,3)$ at $|m| \le 1$, while the $94.1\%$ here is
$l_{\max}=4$ and the $96.0\%$ is $l_{\max}=6$ with a Schwartz cusp
correction and a different solver class.  What the external grid-based
prolate spheroidal result~\cite{tao_mccurdy_rescigno2010} does settle is
the narrower point:\ those coordinates are not cusp-limited.""",
     "D1b: P15 SCOPE -- reverse-ordering claim removed")

edit(SYN,
     r"""At matched angular content the
ordering reverses---Paper~12's own basis with $|m| \le 1$ reaches
${\sim}99.1\%$.""",
     r"""At matched $m_{\max}$ no comparison is
available:\ Paper~12's own basis with $|m| \le 1$ reaches ${\sim}99.1\%$,
but that figure and the $96.0\%$ differ in truncation, cusp treatment and
solver class alike.""",
     "D1c: synthesis -- reverse-ordering claim removed")

edit(SYN,
     r"""for He, ${\sim}99.1\%$ of $D_e$ for $\mathrm{H}_2$ in prolate
spheroidal coordinates with the azimuthal channels open (against
$96.0\%$ in the molecule-frame hyperspherical treatment, cusp-corrected),
and---in the""",
     r"""for He, ${\sim}99.1\%$ of $D_e$ for $\mathrm{H}_2$ in prolate
spheroidal coordinates at $|m| \le 1$ and $96.0\%$ (cusp-corrected) in
the molecule-frame hyperspherical treatment---two routes to the same
system, not a matched pair---and---in the""",
     "D1d: synthesis abstract -- 'against' removed, |m|<=1 restored")

# ============ LARGE: the cusp-advantage zombie, synthesis + its P15 source ==
edit(SYN,
     r"""a structural advantage over
prolate spheroidal coordinates, where the cusp is a coordinate
singularity""",
     r"""a structural difference from
prolate spheroidal coordinates, where the same locus is not a coordinate
surface---a difference in where the cusp sits, which does not translate
into an accuracy advantage""",
     "D0: synthesis cusp-advantage zombie removed")

edit(P15,
     r"""This is the structural advantage of Level~4 coordinates: the cusp
is always a boundary condition, never a coordinate singularity""",
     r"""This is the structural difference of Level~4 coordinates: the cusp
is always a boundary condition rather than a coordinate surface""",
     "B7: P15 body -- the source sentence the synthesis copied")

# ============ LARGE: the coined level index ================================
edit(SYN,
     r"""2$'$ & $\mathrm{H}_2$ (2, 2) & prolate spheroidal, $|m|\le1$ & ${\sim}99.1\%$ of $D_e$ & \cite{loutey_paper12} \\""",
     r"""4$^{\ast}$ & $\mathrm{H}_2$ (2, 2) & prolate spheroidal, $|m|\le1$ & ${\sim}99.1\%$ of $D_e$ (envelope $99.0$--$99.1$) & \cite{loutey_paper12} \\""",
     "D2: coined level '2-prime' replaced by a same-cell marker")

# ============ factual: not all three discrepancies vanish at m = 0 =========
edit(P12,
     r"""Eq.~\eqref{eq:neumann_sigma}, in which all three discrepancies
vanish identically;\ Sec.~\ref{sec:azimuthal} is the first use of the
general form.""",
     r"""Eq.~\eqref{eq:neumann_sigma}, which was derived independently---from
the Legendre completeness relation, as stated below---rather than
specialised from the printed general form.  Two of the three
discrepancies do vanish at $m = 0$;\ the $(2l+1)$ does not, and
Eq.~\eqref{eq:neumann_sigma} carries it because that independent
derivation supplies it.  Sec.~\ref{sec:azimuthal} is the first use of the
general form, which is why the error surfaced only now.""",
     "D6: the provenance note corrected -- (2l+1) does not vanish at m=0")

# ============ D7: conclusion quotes a precision the envelope forbids =======
edit(P12,
     r"""The $\sigma$-only value reproduces independently in a Gaussian
    basis to $0.2$~mHa, and the extended value to $0.01$
    percentage points.""",
     r"""The $\sigma$-only value reproduces independently in a Gaussian
    basis to $0.2$~mHa, and the extended value to within the stated
    stability envelope.""",
     "D7: conclusion no longer claims 0.01 pp cross-method agreement")

# ============ D3: the fifth quadrature surface ============================
edit(P12,
     r"""The
six-dimensional repulsion integral is thereby reduced to a finite
sum of algebraic operations on quantum number labels plus this one
1D quadrature.""",
     r"""The
six-dimensional repulsion integral is thereby reduced, \emph{in the
$\sigma$ sector}, to a finite sum of algebraic operations on quantum
number labels plus this one 1D quadrature;\
Sec.~\ref{sec:azimuthal} states what the $|m| \ge 1$ extension costs.""",
     "D3: intro quadrature claim scoped -- the fifth surface")

# ============ synthesis D4: B_l is adaptive quadrature, not a recurrence ===
edit(SYN,
     r"""Legendre moments all satisfy exact three-term recurrences, so the
$V_{ee}$ matrix is exact within the Neumann truncation order with no
quadrature grids and no fitted parameters~\cite{loutey_paper12}.""",
     r"""Legendre moments satisfy exact three-term recurrences, so the
$V_{ee}$ matrix is exact within the Neumann truncation order with no
six-dimensional quadrature and no fitted parameters~\cite{loutey_paper12}.
One seed resists:\ the log-singular second-kind moment $B_l$ is evaluated
by a single one-dimensional adaptive quadrature, its recurrence being
numerically unstable at the base case.""",
     "D4: synthesis no-quadrature claim corrected to match Paper 12")

edit(SYN,
     r"""Prolate spheroidal $V_{ee}$ ($\sigma$) & algebraic & Neumann expansion, recurrence moments \cite{loutey_paper12} \\""",
     r"""Prolate spheroidal $V_{ee}$ ($\sigma$) & algebraic (one seed) & Neumann recurrences $+$ 1-D quadrature for $B_l$ \cite{loutey_paper12} \\""",
     "D4b: synthesis registry row uses the corpus's own 'one seed' vocabulary")

# ============ D3 (synthesis): the conclusion never got the tier fix ========
edit(SYN,
     r"""coordinates for $\mathrm{H}_2$ (${\sim}99.1\%$ in prolate spheroidal
coordinates with the azimuthal channels open~\cite{loutey_paper12};\
$96.0\%$ of""",
     r"""coordinates for $\mathrm{H}_2$ ($96.0\%$, cusp-corrected;\ prolate
spheroidal at $|m| \le 1$ separately reaches
${\sim}99.1\%$~\cite{loutey_paper12})  ($96.0\%$ of""",
     "D3s: synthesis conclusion -- tier restored, attribution unnested")

# ============ B1: the field guide =========================================
edit(FG,
     r"""4 & H$_2$ (2-center, 2e) & Mol-frame hyperspherical & $96.0\%$ $D_e$ & 15""",
     r"""4 & H$_2$ (2-center, 2e) & Mol-frame hyperspherical & $96.0\%$ $D_e$ (cusp-corrected;\ prolate spheroidal at $|m|\le1$ reaches ${\sim}99.1\%$, Paper 12) & 15""",
     "B1: field guide H2 row -- tier + the better route")

by_path = {}
for path, old, new, label in EDITS:
    by_path.setdefault(path, []).append((old, new, label))

applied, failed = [], []
for path, items in by_path.items():
    with io.open(path, encoding="utf-8") as fh:
        t = fh.read()
    for old, new, label in items:
        if old in t:
            t = t.replace(old, new, 1)
            applied.append(label)
        else:
            failed.append("%s   [%s]" % (label, path))
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(t)

print("applied %d edits" % len(applied))
for a in applied:
    print("  +", a)
if failed:
    print("UNMATCHED (%d):" % len(failed))
    for f in failed:
        print("  -", f)
    sys.exit(1)
