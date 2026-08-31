r"""RECOVERY step 1: re-author Paper 26 SS II.C from the pre-truncation PDF.

This restores the PI-delegated LARGE correction that the over-reversion lost:
HEAD still asserts the naive 1/(8Z^2) scaling "matches" the measured
f_E^h1 column.  It does not -- neither in magnitude nor exponent.
Source of truth: debug/recovery_p26_2026_08_29/paper_26_PRE_TRUNCATION_text.txt
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


rep(r"""The fractional energy contribution of $H_{\text{h1,off}}$ scales
as $\sim 1/(8Z^2)$, consistent with the graph validity scaling
identified in Paper~7.  The graph off-diagonal coupling
$\kappa = -1/16$ is $Z$-independent, while the diagonal energies
scale as $Z^2$.  The relative importance of the graph topology
is therefore
\begin{equation}
  \frac{|\kappa|}{Z^2/n^2} = \frac{1}{16 Z^2/n^2},
\end{equation}
which is $12.5\%$ at $Z = 1$, $3.1\%$ at $Z = 2$, and $0.13\%$
at $Z = 10$.  This scaling matches the $f_E^{\text{h1}}$ column
in Table~\ref{tab:decoupling}.""",
    r"""The graph off-diagonal coupling $\kappa = -1/16$ is
$Z$-independent, while the diagonal energies scale as $Z^2$, so a
naive estimate of the relative importance of the graph topology is
\begin{equation}
  \frac{|\kappa|}{Z^2/n^2} = \frac{1}{16 Z^2/n^2},
\end{equation}
which predicts a $Z^{-2}$ decay ($3.1\%$ at $Z = 2$ for $n = 1$).
The measured $f_E^{\text{h1}}$ column of
Table~\ref{tab:decoupling} decays only as $\sim Z^{-0.85}$
(log--log fit over $Z = 2$--$10$) and sits $12\times$--$79\times$
\emph{above} the naive estimate across the table.

\emph{Corrected 2026-08-29:}\ an earlier version of this subsection
asserted that the naive scaling ``matches'' the $f_E^{\text{h1}}$
column; it does not, in either magnitude or exponent.  The
discrepancy is structural:\ $f_E^{\text{h1}}$ aggregates
off-diagonal coupling over the full spectrum, not a single gap, and
no quantitative bridge from $\kappa$ to this column is currently
established.""",
    "SS II.C false-match claim replaced")

rep(r"""the graph-native CI over-binds (violating the variational bound).
The energy-entanglement decoupling sharpens this picture: at
$Z < Z_c$, the graph topology not only overcounts the energy but
may also distort the entanglement structure by mixing shells
beyond what the exact $V_{ee}$ interaction warrants.""",
    r"""the graph-native CI over-binds (violating the variational bound).
Whether the slow $Z^{-0.85}$ decay measured here connects
quantitatively to that boundary is open:\
Table~\ref{tab:decoupling} starts at $Z = 2$ and does not probe the
sub-$Z_c$ regime.""",
    "SS II.C Z_c paragraph")

io.open(P, "w", encoding="utf-8").write(s)
print(f"\n{n} edits applied -- live false claim removed")
