# Paper 51 §12.7 update draft — G4-5 closure (placeholder)

**Date:** 2026-05-29
**Purpose:** Proposed §12.7 (and adjacent §12.8) subsection text for `papers/group5_qed_gauge/paper_51_gravity_arc.tex`, to be applied once the four parallel G4-5 sub-sprints (G4-5a-refined, G4-5b, G4-5c, G4-5d) close. Mirrors the existing §12 G4-4 subsection structure (12.1 architectural inputs, 12.2 G4-4a constant warp, 12.3 G4-4b variable warp, 12.4 G4-4c conical defect, 12.5 G4-4d/e/f Seeley–DeWitt + replica derivative, 12.6 implications).

**This is a planning document.** Do NOT edit `paper_51_gravity_arc.tex` directly until parallel sprints close and `{{TBD from G4-5x}}` placeholders are filled. Once ready, the §12.7 subsection slots in immediately after §12.6 (`subsec:implications` at line ~1167) and before §13 (`sec:discussion` at line ~1204).

---

## Proposed insertion: §12.7 G4-5 — Discrete replica method for $S_{\rm BH}$

```latex
% =====================================================================
\subsection{G4-5: discrete replica method for $S_{\BH}$}
\label{subsec:g4_5}
% =====================================================================

The G4-4c tip coefficient $\Delta_K^{\rm Dirac, tip}(\alpha) = -\tfrac{1}{12}(1/\alpha - \alpha)$ (Eq.~\ref{eq:spinor_sc_tip}) and G4-4f replica derivative $+1/6$ (Eq.~\ref{eq:replica_derivative}) are the load-bearing structural inputs for the replica-method derivation of $S_{\rm BH}$ on the discrete substrate. G4-5 integrates these over $t$ with a Connes--Chamseddine cutoff function $f(t\Lambda^2)$ to produce the BH entropy directly from the warped Dirac heat trace.

\subsubsection{Inputs from G4-3 and G4-4}
\label{subsubsec:g4_5_inputs}

G4-5 reuses the G4-3a-cleanup Hermitian polar Laplacian, the G4-3b variable warp $r(\rho) = r_h\sqrt{1 + (\rho/r_h)^2}$, the G4-4a constant-warp factorization $K_{\rm cigar}^{\rm Dirac}(t) = K_{D^2}^{\rm Dirac}(t) \cdot K_{S^2}^{\rm Dirac}(t)$, the G4-4c \texttt{DiscreteWedgeDirac} class for the conical-defect spinor, and the G4-4f central finite-difference replica derivative protocol. Substrate convention: $R = 10$, $a = 0.05$, $N_\rho = 200$, $N_0 = 120$, $r_h = 2$, $\varepsilon = 12/120 = 0.1$ — the G4-4c week 3 sweet spot.

The replica formula
\begin{equation}
  S_{\BH}
  = -\frac{dI_E}{d\alpha}\bigg|_{\alpha = 1},
  \quad
  I_E(\alpha) = -\frac{1}{2}\int_0^\infty \frac{dt}{t}\, \tilde{f}(t\Lambda^2)\, K^{\rm Dirac}_{\rm cigar}(\alpha, t)
  \label{eq:replica_full}
\end{equation}
splits into a \emph{tip contribution} (the $\alpha$-derivative piece, identified at the G4-4 tip coefficient) and a \emph{bulk contribution} (the disk part, $\alpha$-independent). G4-5 closes both.

\subsubsection{G4-5a: tip-only replica integration}
\label{subsubsec:g4_5a}

The tip-only $S_{\rm tip}$ piece:
\begin{equation}
  S_{\rm tip}(\Lambda) = +\frac{1}{2}\int_{t_{\rm UV}}^{t_{\rm IR}} \frac{dt}{t}\, \tilde{f}(t\Lambda^2)\, \Delta'(t), \quad \Delta'(t) = \frac{dK_{\rm wedge}^{\rm Dirac}}{d\alpha}\bigg|_{\alpha=1} - K_{\rm disk}^{\rm Dirac}(t).
  \label{eq:s_tip_integral}
\end{equation}
G4-5a first move (sprint scale) verified the integration is operational: $S_{\rm tip}$ extracted as finite, positive, $\Lambda$-decreasing across $\Lambda \in \{0.5, 1.0, 1.5, 2.0\}$. The discrete/continuum ratio (0.13--0.37) reflected the sprint-scale $t$-grid starting at $t = 0.1$, missing the UV regime $t \in [a^2, 0.1]$ where the Gaussian cutoff peaks at large $\Lambda$.

G4-5a-refined extends the $t$-grid into the substrate UV with $4/\pi^2$ azimuthal-truncation correction (G4-3d-UV finding) and higher-order quadrature. {{TBD from G4-5a-refined}}: ratio closes to {{TBD}}\% across $\Lambda \in \{0.5, 1.0, 1.5, 2.0\}$, with the $\Lambda$-trend becoming flat after proper integration.

\subsubsection{G4-5b: bulk Weyl extraction ($\Lambda^4 + \Lambda^2$)}
\label{subsubsec:g4_5b}

The disk part $K_{\rm disk}^{\rm Dirac}(t)$ has Seeley--DeWitt expansion
\begin{equation}
  K_{\rm disk}^{\rm Dirac}(t) \sim \frac{A_{D^2}}{2\pi t} + a_1^{D^2}\, t^{-1/2} + a_2^{D^2} + O(t)
\end{equation}
with $A_{D^2}$ the disk area and $a_1^{D^2}$ a boundary contribution from the Robin/Dirichlet BC at the cigar tip. Mellin-transforming against the Gaussian cutoff:
\begin{equation}
  S_{\rm Weyl}^{(\Lambda^4)} = \frac{A_{D^2}\Lambda^4}{2\pi}\phi(2), \quad S_{\rm EH}^{(\Lambda^2)} = a_2^{D^2}\Lambda^2 \phi(1) + \text{(boundary)}.
\end{equation}

{{TBD from G4-5b}}: bulk $\Lambda^4$ coefficient extracted to {{TBD}}\% recovery against continuum $A_{D^2}/(2\pi)$; bulk $\Lambda^2$ Einstein--Hilbert coefficient at {{TBD}}\% recovery against the continuum prediction $r_h^2/3$ at $A_{\rm horizon} = 4\pi r_h^2$.

\subsubsection{G4-5c: joint warp + conical defect}
\label{subsubsec:g4_5c}

The physical cigar geometry combines smooth-warp at the asymptotic end with a conical defect at the horizon. The G4-4b factorization-loss $\sim (r'/r)^2$ and the G4-4c tip coefficient $-\tfrac{1}{12}(1/\alpha - \alpha)$ must be combined on a single discrete substrate.

A new class \texttt{VariableWarpConicalDirac} extending \texttt{VariableWarpDirac} (G4-4b) with the conical-defect azimuthal BC of \texttt{DiscreteWedgeDirac} (G4-4c) implements this. The load-bearing F11 falsifier:
\begin{itemize}
  \item F6 extension: at $\alpha = 1$ and constant warp, reduce bit-exact to G4-4a's \texttt{WarpedDiracConstant}. {{TBD from G4-5c}}: residual $\sim 10^{-16}$.
  \item F4 tip regularity at variable warp + conical defect. {{TBD from G4-5c}}: {{TBD}}\% recovery at the sweet spot.
  \item Joint factorization-loss + tip coefficient: $\Delta_K^{\rm Dirac}(\alpha, r_h, t) \approx -\tfrac{1}{12}(1/\alpha - \alpha) \cdot K^{\rm Dirac}_{S^2}(t; r_h) + O((r'/r)^2)$.
\end{itemize}

{{TBD from G4-5c}}: F11 closure verdict.

\subsubsection{G4-5d: cutoff-function dependence}
\label{subsubsec:g4_5d}

G8 (\S\ref{sec:g8}) predicts $S_{\BH}(f) \propto \phi(2)$ at the continuum level. G4-5d verifies this at the discrete replica-method level by computing $S_{\BH}^{\rm discrete}$ for three named cutoff classes:
\begin{center}
  \begin{tabular}{lccc}
    \toprule
    Cutoff & $f(x)$ & $\phi(2)$ & discrete $S_{\BH}^{(f)} / S_{\BH}^{(\rm Gaussian)}$ \\
    \midrule
    Gaussian & $e^{-x}$ & $1$ & $1$ \\
    Sharp & $\Theta(1 - x)$ & $1/2$ & {{TBD from G4-5d}} \\
    Polynomial & $e^{-x^2}$ & $1/2$ & {{TBD from G4-5d}} \\
    \bottomrule
  \end{tabular}
\end{center}
{{TBD from G4-5d}}: F12 closure verdict. Expected: $\phi(2)$ ratios within {{TBD}}\% at the sweet spot, lifting G8's continuum prediction to a discrete-substrate theorem.

\subsubsection{Headline closure}
\label{subsubsec:g4_5_headline}

The discrete-substrate BH entropy:
\begin{equation}
  \boxed{\;
  S_{\BH}^{\rm discrete}(r_h, \Lambda; f) = \frac{r_h^2 \Lambda^2}{3} \cdot \phi(2)
  + O\!\left(\frac{a^2}{r_h^2}\right)
  + O\!\left(\frac{r_h^2}{R^2}\right)\;
  }
  \label{eq:s_bh_discrete}
\end{equation}
matches the continuum G4-2 derivation (Eq.~\ref{eq:s_bh_continuum_g4_2}) at {{TBD}}\% recovery at the sweet spot. The cutoff-function $\phi(2)$ dependence is verified at {{TBD}}\% across the three named cutoff classes.

\subsubsection{Implications for the full $S_{\BH}$ closure (G4-6)}
\label{subsubsec:g4_5_implications}

G4-5 closes the LEADING-order $r_h^2 \Lambda^2 / 3 \cdot \phi(2)$ piece of $S_{\BH}$. The remaining work for G4-6 (full closure):
\begin{itemize}
  \item Subleading $O(a^2/r_h^2)$ UV corrections via multi-substrate extrapolation;
  \item Subleading $O(r_h^2/R^2)$ IR corrections via boundary regularization;
  \item Higher-curvature contributions to Wald entropy ($\Lambda^0$ Euler density, $\Lambda^{-2}$ Riemann squared);
  \item $\alpha > 1$ branch closure (excess-angle regime structural asymmetry, G4-4c open).
\end{itemize}
Total G4-6 estimate post G4-5 closure: 7--11 months.

\paragraph{Structural-skeleton-scope reading.} G4-5 confirms the predicted division: the framework predicts the structural form $r_h^2 \Lambda^2 / 3 \cdot \phi(2)$ of $S_{\BH}$ (framework prediction in the strong sense) but DOES NOT autonomously fix the cutoff function $f$ (Class 1 calibration, G8 / \S\ref{sec:g8}). The discrete replica method functions as designed: structural quantities are bit-exact, calibration data is external. Consistent with the structural-skeleton-scope reading (\S\ref{sec:discussion}).
```

---

## Cross-reference updates required (when applying)

When this section lands, the following cross-references in adjacent Paper 51 sections should be checked:

1. **§13 Discussion** (`sec:discussion`, line ~1204): add G4-5 to the list of cutoff-dependent calibrated predictions.
2. **§14 Open questions** (`sec:open`, line ~1310): remove "G4-5 multi-month discrete replica method" if it appears, replace with "G4-6 sub-leading corrections to $S_{\BH}$."
3. **§12.6 Implications for multi-month $S_{\BH}$ program** (`subsec:implications`, line ~1167): update the closing paragraph to point forward to §12.7 instead of "remaining multi-month work."
4. **CHANGELOG.md**: add G4-5 version bump entry with sub-sprint roll-call.
5. **CLAUDE.md §2**: add G4-5 one-liner per [`feedback_changelog_for_chronicle`](memory/feedback_changelog_for_chronicle.md) discipline.

---

## Placeholder count

This draft contains the following placeholders, to be filled at G4-5 closure:

| Marker | Source | Count |
|---|---|---|
| `{{TBD from G4-5a-refined}}` | G4-5a-refined sub-sprint | 5 |
| `{{TBD from G4-5b}}` | G4-5b sub-sprint | 4 |
| `{{TBD from G4-5c}}` | G4-5c sub-sprint | 5 |
| `{{TBD from G4-5d}}` | G4-5d sub-sprint | 3 |
| `{{TBD}}` (recovery percentages) | aggregate | ~12 |

Total: ~29 placeholders. When all four parallel sub-sprints close, the synthesis is a mechanical fill-in operation taking ~30 minutes.

---

## Verification checklist before applying

Before the actual §12.7 update is applied to Paper 51, verify:

1. ☐ All four parallel sub-sprint memos exist in `debug/`.
2. ☐ G4-5a-refined ratio closure ≥ 90% at sweet spot.
3. ☐ G4-5b $\Lambda^4$ and $\Lambda^2$ recoveries ≥ {{decision-gate %}}.
4. ☐ G4-5c F11 closure with bit-exact F6 reduction.
5. ☐ G4-5d $\phi(2)$ ratios within {{decision-gate %}}.
6. ☐ Aggregate $S_{\BH}^{\rm discrete}$ recovery ≥ {{decision-gate %}} against continuum $r_h^2 \Lambda^2 / 3$.
7. ☐ Synthesis memo `g4_5_synthesis_memo.md` populated.
8. ☐ Paper 51 §12.6 cross-reference updated to point forward to §12.7.
9. ☐ All `{{TBD}}` placeholders replaced with numerical values.
10. ☐ LaTeX three-pass compile clean (zero substantive warnings).

If any item fails, the failed sub-sprint(s) get a follow-on diagnostic per [`feedback_diagnostic_before_engineering`](memory/feedback_diagnostic_before_engineering.md) before the §12.7 update applies.

---

## Files referenced

- `papers/group5_qed_gauge/paper_51_gravity_arc.tex` (target file; do NOT edit until closure)
- `debug/g4_5_synthesis_memo.md` (companion synthesis memo, placeholder)
- `debug/g4_5_scoping_memo.md` (scope)
- `debug/g4_5a_first_move_tip_replica_memo.md` (G4-5a baseline)
- {{TBD from G4-5a-refined / G4-5b / G4-5c / G4-5d memos}}
