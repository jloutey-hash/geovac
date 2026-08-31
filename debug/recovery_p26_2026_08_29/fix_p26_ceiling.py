r"""Paper 26: correct the N/O/F hub values and RETIRE the ceiling-attainment
claim (both rested on the diagonal-only entropy routine fixed 2026-08-29).

Measured with the corrected full-RDM entropies:
  member-representative I(2p-1,2p+1):  0.566 (N)  0.014 (O)  0.982 (F)
     [retired: 1.837 / 1.785 / 1.644]
  attainable max over sampled members: 1.431 (N)  1.386 (O)  1.386 (F)
     -- O and F land on ln 4 = 1.38629 EXACTLY; ln 8 / ln 16 are NOT attained.
  I_cv: 4.3e-4 (N, approximate), 0 (O, F, exact)
"""
import io

P = "papers/group6_precision_observations/paper_26_entanglement.tex"
s = io.open(P, encoding="utf-8").read()
n = 0


def rep(old, new, label):
    global s, n
    assert s.count(old) == 1, f"{label}: {s.count(old)} matches"
    s = s.replace(old, new)
    n += 1
    print(f"  ok  {label}")


rep(r"""The hub transitions to the $2p(\pm 1)$ pair, which carries the
mutual information: $I(2p_{-1}, 2p_{+1}) = 1.837$ (N),
$1.785$ (O), $1.644$ (F) \emph{for the multiplet member the eigensolver
returns} --- these ground states are 4-, 6- and 4-fold degenerate, and
the value ranges from $0$ (single-determinant members) up to the
multiplet's information-theoretic ceiling, which is attained ---
$\ln 8 \approx 2.08$ for N/F, $\ln 16 \approx 2.77$ for O
(\S\ref{sec:corevalence}, degeneracy caveat).""",
    r"""The hub transitions to the $2p(\pm 1)$ pair, which carries the
mutual information: $I(2p_{-1}, 2p_{+1}) = 0.566$ (N),
$0.014$ (O), $0.982$ (F) \emph{for the multiplet member the eigensolver
returns} --- these ground states are 4-, 3- and 4-fold degenerate, and
the value ranges from $0$ (single-determinant members) up to a
member-dependent maximum measured at $1.43$ (N), $1.386$ (O), $1.386$ (F)
(\S\ref{sec:corevalence}, degeneracy caveat).  (Values corrected
2026-08-29:\ the retired $1.837/1.785/1.644$, and a claimed
$\ln 8$/$\ln 16$ ceiling, came from a single-orbital entropy routine that
built $\rho_i$ diagonal-only and so discarded the spin coherence
$\langle{\uparrow}\rvert\rho_i\lvert{\downarrow}\rangle$ --- nonzero for
exactly these sector-mixed multiplet members.)""",
    "N-F hub values")

rep(r"""determinant basis states are themselves exact ground eigenstates with
every orbital mutual information $0$, verified directly --- (for O,
$|1s^2 2s^2 2p_{-1}^2 2p_0^2\rangle$ is the unique determinant in its
$(M_L, M_S)$ sector and hence an exact ground eigenstate) for which
\emph{every} orbital mutual information is $0$, alongside strongly
correlated members with $I(2p_{-1}, 2p_{+1})$ approaching the
multiplet's information-theoretic ceiling ($\ln 8 \approx 2.08$ for
N/F, $\ln 16 \approx 2.77$ for O; both ceilings are \emph{attained
exactly} by explicit members --- the uniform superposition of the
support determinants for N/F, and amplitudes $(1/2, 1/2,
1/(2\sqrt{2})\times 4)$ on the paired/mixed determinants for O ---
pinned to $10^{-9}$ by the backing test); the
values previously quoted here ($1.837/1.785/1.644$) are
member-representative, not basis constants, and even the qualitative
``the valence network is alive'' reading is member-dependent.""",
    r"""determinant basis states are themselves exact ground eigenstates with
every orbital mutual information $0$ (for O,
$|1s^2 2s^2 2p_{-1}^2 2p_0^2\rangle$ is the unique determinant in its
$(M_L, M_S)$ sector and hence an exact ground eigenstate), alongside
strongly correlated members reaching a measured maximum of $1.43$ (N),
$1.386$ (O), $1.386$ (F) --- the last two being $\ln 4$ to all quoted
digits.  The member-representative values ($0.566/0.014/0.982$) are not
basis constants, and even the qualitative ``the valence network is
alive'' reading is member-dependent.

\emph{Retracted 2026-08-29:\ ceiling attainment.}  This paragraph
previously asserted an information-theoretic ceiling of
$\ln 8 \approx 2.08$ (N/F) and $\ln 16 \approx 2.77$ (O), ``attained
exactly by explicit members \dots\ pinned to $10^{-9}$ by the backing
test''.  Both the quoted ceilings and their attainment were artifacts of
a single-orbital entropy routine that built $\rho_i$ diagonal-only,
discarding the spin coherence
$\langle{\uparrow}\rvert\rho_i\lvert{\downarrow}\rangle$; because a
diagonal Shannon entropy majorizes the true von Neumann entropy, the
routine over-reported $s_i$ (e.g.\ $0.798$ against a true $0.313$) and
the resulting ``mutual information'' violated subadditivity on half the
orbital pairs.  With the corrected full reduced density matrix the
explicit constructions reach $1.398$/$0$/$1.386$, not the claimed
ceilings.  No analytic ceiling is asserted in its place:\ the measured
maxima above are MEASURED, and the qualitative claim they support ---
that the value is member-dependent and ranges from $0$ up to an
$O(1)$ maximum --- is unaffected.""",
    "ceiling retraction")

io.open(P, "w", encoding="utf-8").write(s)
print(f"\n{n} edits applied")
