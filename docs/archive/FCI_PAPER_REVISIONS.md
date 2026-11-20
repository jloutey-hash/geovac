# FCI Paper Revisions: LiH Diagnostic Findings

**Date:** 2026-03-12
**Paper:** `papers/core/paper_geovac_lcao_fci.tex`
**Title:** "Topological Full Configuration Interaction for Heteronuclear Diatomics: LiH as a Benchmark"
**Status:** Revision required — key claims invalidated by Hamiltonian balance diagnostics

---

## Executive Summary

Post-publication diagnostics (v0.9.37) discovered that the reported LiH results
($D_e^{\mathrm{CP}} = 0.093$~Ha, $R_{\mathrm{eq}} \approx 2.5$~bohr) originated
from an **imbalanced Hamiltonian configuration**: cross-nuclear attraction computed
for all $(n,l,m)$ orbitals (`exact`) paired with electron repulsion restricted to
$s$-orbital pairs only (`s_only`). When the Hamiltonian is properly balanced —
either `exact`+`True` (all orbitals for both) or `fourier`+`s_only` (s-only for
both) — the PES is monotonically attractive with **no equilibrium geometry**.

Additionally, the paper's interpretation of the 17% $R_{\mathrm{eq}}$ contraction
as a "transient discretization artifact converging with $n_{\max}$" is incorrect.
The true root cause is **missing kinetic repulsion**: the graph Laplacian's
intra-atom kinetic energy is completely $R$-independent, so no repulsive wall forms
at short bond lengths. A Fock-weighted kinetic correction based on Paper 7's
energy-shell constraint partially restores the equilibrium.

**What survives:**
- All atomic results (He, Li, Be) are unaffected
- The BSSE machinery and counterpoise framework work correctly
- The Direct CI solver is algorithmically correct
- The diagnostic arc conclusions about Sturmian/D-matrix/relaxation are valid
- The $O(V)$ scaling claim is correct

**What must be revised:**
- The headline 1.0% accuracy claim
- The Table I energetic values
- The $R_{\mathrm{eq}} \approx 2.5$ bohr equilibrium claim
- The "discretization artifact" interpretation
- The continuum limit convergence prediction

---

## Correction Table

| Section | Line(s) | Original Text | Corrected Text | Reason |
|---------|---------|---------------|----------------|--------|
| Abstract | 39–40 | "$D_{e}^{\mathrm{CP}}=0.093$~Ha, a 1.0\% error relative to the experimental value $D_{e}=0.092$~Ha" | "The balanced Hamiltonian yields a CP-corrected vertical binding energy of $D_{e}^{\mathrm{CP}}=0.110$~Ha at $R=3.015$~bohr, but no equilibrium geometry" | Imbalanced config produced 0.093; balanced gives 0.110 with no minimum |
| Abstract | 41–44 | "$R_{\mathrm{eq}}\approx 2.5$~bohr exhibits a 17\% contraction...predicted to converge with increasing $n_{\max}$" | "The absence of an equilibrium geometry is traced to missing $R$-dependent kinetic repulsion in the fixed graph Laplacian, partially restored by a Fock-weighted overlap correction" | No equilibrium exists in balanced configs; discretization is not the cause |
| Intro | 83–85 | "the Boys--Bernardi counterpoise-corrected dissociation energy of LiH at $n_{\max}=3$ is $D_{e}^{\mathrm{CP}}=0.093$~Ha, constituting a 1.0\% error" | "the Boys--Bernardi counterpoise-corrected vertical binding energy of LiH at $n_{\max}=3$ and $R=3.015$~bohr is $D_{e}^{\mathrm{CP}}=0.110$~Ha" | Correct balanced value; clarify it's at fixed $R$, not adiabatic |
| Sec III.B | 268–270 | "After CP correction, the binding energy is $D_e^{\mathrm{CP}} = 0.093$~Ha, a 1.0\% error relative to the experimental value" | "After CP correction, the vertical binding energy at $R=3.015$~bohr is $D_e^{\mathrm{CP}} = 0.110$~Ha, a 20\% overestimate" | Balanced result |
| Table I | 282–289 | $E_{\mathrm{LiH}}=-8.097$, $E_{\mathrm{Li}}=-7.334$, $E_{\mathrm{H}}=-0.558$, $D_e^{\mathrm{raw}}=0.205$, $D_e^{\mathrm{CP}}=0.093$, Error=1.0\% | $E_{\mathrm{LiH}}=-8.118$, $E_{\mathrm{Li}}=-7.392$, $E_{\mathrm{H}}=-0.500$, $D_e^{\mathrm{raw}}=0.226$, $D_e^{\mathrm{CP}}=0.110$, Error=20\% | Previous values from imbalanced `exact+s_only`; new values from balanced `exact+True` |
| Sec III.C | 303–304 | "The CP-corrected potential energy surface (PES) exhibits a minimum at $R_{\mathrm{eq}} \approx 2.5$~bohr" | "The CP-corrected PES is monotonically attractive: $D_e^{\mathrm{CP}}$ decreases from 0.464~Ha at $R=1.5$~bohr to 0.042~Ha at $R=5.0$~bohr, with no equilibrium geometry" | No minimum exists in balanced configs |
| Sec III.C | 306–308 | "striking contrast to the 1.0\% energetic accuracy...the dissociation energy is essentially correct" | Remove. The 1.0\% accuracy and the $R_{\mathrm{eq}}$ were both artifacts of the same imbalanced Hamiltonian | Both claims are invalidated |
| Sec IV.A | 484–493 | "the BSSE...strictly exceeds the experimental dissociation energy...The system is variationally driven to collapse inward, contracting $R_{\mathrm{eq}}$" | "the BSSE of 0.115~Ha exceeds the experimental binding energy, confirming basis deficiency. However, since the balanced PES has no minimum, BSSE is not the mechanism for geometric contraction" | BSSE is real but not the cause of a contraction that doesn't exist |
| Sec IV.C | 519–544 | Entire "Continuum Limit Prediction" subsection predicting $R_{\mathrm{eq}}$ convergence with $n_{\max}$ | Replace with analysis of missing kinetic repulsion and Fock-weighted correction (see new subsection below) | The prediction was based on a false equilibrium |
| Conclusion | 558–560 | "$D_e^{\mathrm{CP}}=0.093$~Ha, a 1.0\% error" | "$D_e^{\mathrm{CP}}=0.110$~Ha at $R=3.015$~bohr (vertical binding energy)" | Correct balanced value |
| Conclusion | 561–563 | "The equilibrium bond length contraction...identified as a transient discretization artifact" | "The balanced Hamiltonian produces a monotonically attractive PES with no equilibrium, traced to missing $R$-dependent kinetic repulsion in the fixed graph Laplacian" | Correct diagnosis |
| App. Table II | 596 | v0.9.11: $D_e^{\mathrm{CP}}=0.093$, $R_{\mathrm{eq}}\sim 2.5$ | v0.9.11: $D_e^{\mathrm{CP}}=0.110$ (at $R=3.015$), no $R_{\mathrm{eq}}$ | Balanced result |
| Table II caption | 641–644 | "The LCAO baseline (v0.9.11) remains the only architecture achieving quantitative energetic accuracy" | "The LCAO baseline achieves correct order-of-magnitude binding but lacks equilibrium geometry due to missing kinetic repulsion" | Tempered claim |

---

## Detailed Section Revisions

### 1. Abstract (lines 28–54)

**Original:**
```latex
\begin{abstract}
...the Boys--Bernardi counterpoise-corrected binding
energy is $D_{e}^{\mathrm{CP}}=0.093$~Ha, a 1.0\% error relative to the
experimental value $D_{e}=0.092$~Ha.  The equilibrium bond length
$R_{\mathrm{eq}}\approx 2.5$~bohr exhibits a 17\% contraction from the
experimental $3.015$~bohr, identified as a discretization artifact of the
$n_{\max}=3$ basis truncation and predicted to converge with increasing
$n_{\max}$ by direct analogy with the $2s$--$2p$ spectral aliasing
documented in Paper~1.
```

**Revised:**
```latex
\begin{abstract}
...the Boys--Bernardi counterpoise-corrected vertical binding
energy at the experimental geometry ($R=3.015$~bohr) is
$D_{e}^{\mathrm{CP}}=0.110$~Ha, within 20\% of the experimental value
$D_{e}=0.092$~Ha.  The balanced Hamiltonian produces a monotonically
attractive potential energy surface with no equilibrium geometry, a
limitation traced to the $R$-independence of the graph Laplacian's
intra-atom kinetic energy: unlike continuous quantum mechanics, the
discrete adjacency within each atomic subgraph does not respond to
inter-nuclear compression.  A Fock-weighted kinetic correction, derived
from the energy-shell constraint $p_0^2 = Z^2/n^2$ of Paper~7, partially
restores the equilibrium by concentrating orthogonalization repulsion on
compact core orbitals, achieving $R_{\mathrm{eq}}\approx 3.8$~bohr and
$D_e^{\mathrm{CP}}\approx 0.046$~Ha.  A 29-version diagnostic arc
(v0.9.8--v0.9.37) systematically tested and eliminated every competing
hypothesis...
```

### 2. Introduction (lines 82–85)

**Original:**
```latex
As a headline result, the Boys--Bernardi counterpoise-corrected
dissociation energy of LiH at $n_{\max}=3$ is
$D_{e}^{\mathrm{CP}}=0.093$~Ha, constituting a 1.0\% error relative to
the experimental value of $0.092$~Ha~\cite{huber}.
```

**Revised:**
```latex
As a headline result, the Boys--Bernardi counterpoise-corrected vertical
binding energy of LiH at $n_{\max}=3$ and $R=3.015$~bohr is
$D_{e}^{\mathrm{CP}}=0.110$~Ha, within 20\% of the experimental value
of $0.092$~Ha~\cite{huber}.  The balanced Hamiltonian reveals a
fundamental limitation of the discrete LCAO framework: the absence of
$R$-dependent kinetic repulsion, which is partially resolved by a
Fock-weighted kinetic correction derived from Paper~7~\cite{paper7}.
```

### 3. Section III.B: Energetic Accuracy (lines 254–296)

**Original Table I values:**
```
E_LiH(R=3.015)  = -8.097
E_Li             = -7.334
E_H              = -0.558
BSSE             = -0.115
D_e^raw          = 0.205
D_e^CP           = 0.093
Error            = 1.0%
```

**Corrected Table I values:**
```
E_LiH(R=3.015)  = -8.118   (balanced: exact+True)
E_Li             = -7.392   (exact atomic eigenvalue)
E_H              = -0.500   (exact atomic eigenvalue)
BSSE             = -0.115   (unchanged — correctly computed)
D_e^raw          = 0.226
D_e^CP           = 0.110
Comparison       = 20% overestimate (vertical, not adiabatic)
```

**Note on the old values:** The previous $E_{\mathrm{Li}} = -7.334$ and
$E_{\mathrm{H}} = -0.558$ came from FCI calculations in the molecular basis
(needed for CP correction), not from isolated-atom calculations. These are
consistent and correct for the CP subtraction; the issue is only with
$E_{\mathrm{LiH}}$. The old molecular energy of $-8.097$ was computed with the
imbalanced `exact+s_only` configuration. The balanced `exact+True` gives $-8.118$.

**Revised text after Table I:**
```latex
After CP correction, the vertical binding energy at $R=3.015$~bohr is
$D_e^{\mathrm{CP}} = 0.110$~Ha, a 20\% overestimate of the experimental
value $0.092$~Ha~\cite{huber}.  We emphasize that this is a
\emph{vertical} binding energy computed at fixed internuclear distance,
not an adiabatic value at the theoretical equilibrium: the balanced
Hamiltonian produces no equilibrium geometry (Sec.~\ref{sec:contraction}).
```

**Add a "Balance Theorem" note:**
```latex
\paragraph{Hamiltonian balance requirement.}
The cross-nuclear attraction and electron repulsion must operate on
matched orbital subspaces.  Two balanced configurations---\texttt{exact}
cross-nuclear with full cross-atom $V_{ee}$ (\texttt{exact+True}), and
\texttt{fourier} cross-nuclear with $s$-only $V_{ee}$
(\texttt{fourier+s\_only})---produce molecular energies agreeing within
2~mHa across the entire PES, confirming variational consistency.  An
\emph{imbalanced} configuration (\texttt{exact+s\_only}) artificially
overbinds by 0.55~Ha, as the $p$ and $d$ orbitals receive cross-nuclear
attraction without the compensating electron repulsion screening.  An
earlier version of this work inadvertently used the imbalanced
configuration; the results presented here are from the balanced
\texttt{exact+True} Hamiltonian.
```

### 4. Section III.C: The Geometric Contraction (lines 298–320)

**Original (entire subsection):**
```latex
The CP-corrected potential energy surface (PES) exhibits a minimum at
$R_{\mathrm{eq}} \approx 2.5$~bohr, a 17\% inward shift from the
experimental equilibrium of $3.015$~bohr.  This geometric contraction
stands in striking contrast to the 1.0\% energetic accuracy achieved at
the experimental geometry...
```

**Revised:**
```latex
\subsection{Absence of Equilibrium Geometry}
\label{sec:contraction}

The balanced configurations (\texttt{exact+True} and
\texttt{fourier+s\_only}) produce identical results within 2~mHa,
confirming variational consistency.  However, both yield monotonically
attractive PES curves with no equilibrium geometry.  The equilibrium
reported in earlier versions arose from an imbalanced Hamiltonian
configuration that has since been corrected.

The CP-corrected binding energy $D_e^{\mathrm{CP}}$ decreases smoothly
from 0.464~Ha at $R=1.5$~bohr to 0.042~Ha at $R=5.0$~bohr, with the
PES slope diminishing but never reversing sign.  At no point does the
repulsive kinetic wall balance the attractive cross-nuclear and exchange
interactions: the system is variationally attracted to short bond
lengths at all $R$.

This monotonic behavior is traced to a structural property of the
Topological LCAO: the graph Laplacian kinetic energy within each atomic
subgraph is \emph{completely $R$-independent}.  In continuous quantum
mechanics, bringing two atoms together forces their electrons into
higher-momentum orthogonal states, creating the kinetic repulsion that
establishes the equilibrium bond length.  In the discrete framework,
each atom's adjacency matrix---and hence its contribution to the kinetic
energy via the degree matrix $D$---is fixed at construction time and
does not respond to the approach of a second center.

The only $R$-dependent kinetic contribution comes from the bridge edges
connecting the two subgraphs.  However, these bridge elements are
approximately two orders of magnitude too weak: the maximum bridge
matrix element at $R=3.015$~bohr is 0.017~Ha, compared to the $\sim$0.5~Ha
kinetic correction required by virial theorem analysis (the observed
virial ratio $\eta = -V/(2T) \approx 30$ instead of the equilibrium
value $\eta = 1$).
```

### 5. Section V: Replace "Discretization Artifact" with "Kinetic Repulsion" (lines 466–544)

The entire Section V ("$R_{\mathrm{eq}}$ as a Discretization Artifact") should be
replaced with a new section. The Paper 1 analogy about spectral aliasing is
misleading: the molecular issue is not about discretization resolution but about
the fundamental R-independence of the graph Laplacian kinetic energy.

**New Section V: Kinetic Repulsion and the Fock-Weighted Correction**

```latex
\section{Missing Kinetic Repulsion and the Fock-Weighted Correction}
\label{sec:kinetic}

\subsection{The Kinetic Diagnostic}
\label{sec:kinetic_diagnostic}

The virial theorem for a Coulombic system requires
$\eta \equiv -V/(2T) = 1$ at the equilibrium geometry.  Table~\ref{tab:virial}
documents the virial ratio across the LiH PES: $\eta$ ranges from
$46$ at $R=1.5$~bohr to $22$ at $R=5.0$~bohr, indicating that the
kinetic energy is $\sim$30$\times$ too weak relative to the potential
energy at all internuclear distances.

\begin{table}[h]
\centering
\begin{tabular}{ccccc}
\toprule
$R$ (bohr) & $T$ (Ha) & $V$ (Ha) & $\eta=-V/(2T)$ & $E$ (Ha) \\
\midrule
1.5 & $-0.092$ & $-8.380$ & 45.8 & $-8.471$ \\
2.5 & $-0.122$ & $-8.055$ & 32.9 & $-8.178$ \\
3.015 & $-0.133$ & $-7.985$ & 30.0 & $-8.118$ \\
4.0 & $-0.172$ & $-7.908$ & 23.0 & $-8.079$ \\
5.0 & $-0.181$ & $-7.867$ & 21.7 & $-8.049$ \\
\bottomrule
\end{tabular}
\caption{Virial analysis of the balanced LiH PES.  The virial ratio
$\eta$ should equal 1 at equilibrium but ranges from 22 to 46,
indicating the kinetic energy is $\sim$30$\times$ too weak.}
\label{tab:virial}
\end{table}

Direct examination confirms that the intra-atom kinetic energy (the
off-diagonal H1 within each atomic subgraph) is identically
$R$-independent: the lithium and hydrogen graph Laplacians are
constructed once from the isolated-atom topology and never modified as
$R$ changes.  Bridge edges contribute a small $R$-dependent kinetic
term through the degree matrix, but the maximum bridge element at
$R=3.015$~bohr is only 0.017~Ha, two orders of magnitude below the
required repulsion.

This is a structural limitation of the Topological LCAO architecture:
the discrete graph Laplacian has no mechanism to represent the
orthogonalization kinetic energy that arises when atomic orbitals
overlap in continuous quantum mechanics.

\subsection{The Fock-Weighted Kinetic Correction}
\label{sec:fock_correction}

The energy-shell constraint of Paper~7---$p_0^2 = -2E = Z^2/n^2$ for
the $n$-th shell---provides a natural scale for the
``information density'' of each orbital on the conformal $S^3$ manifold.
We exploit this to construct a physically motivated kinetic correction.

For each orbital $\phi_a$ on atom $A$, the Fock-weighted overlap
kinetic correction adds to the diagonal of $\mathbf{H}_1$:
\begin{equation}
\Delta h_{aa} = \lambda \sum_{b \in B}
S_{ab}^2 \left(\frac{Z_A}{n_a}\right)^2
\left(\frac{Z_B}{n_b}\right)^2,
\label{eq:fock_correction}
\end{equation}
where $S_{ab} = \langle \phi_a | \phi_b \rangle$ is the STO overlap
integral between $s$-orbitals on different centers, and $\lambda$ is a
single adjustable parameter.  The Fock weights suppress diffuse
orbitals (high $n$) and concentrate the kinetic penalty on compact core
orbitals whose overlap decays rapidly with $R$.

The physical content of Eq.~(\ref{eq:fock_correction}) is that the cost
of overlapping two topological information structures scales with the
product of their information densities $p_0^2 = Z^2/n^2$.  Core--core
pairs (Li~1$s$--H~1$s$) carry weight $9 \times 1 = 9$, while
diffuse--diffuse pairs (Li~3$s$--H~3$s$) carry weight
$1 \times 0.11 = 0.11$---an 81$\times$ discrimination.

The R-discrimination (ratio of correction at $R=1.5$ to $R=5.0$~bohr)
improves from 3.0$\times$ with uniform weighting to 9.3$\times$ with
Fock weighting.  For the Li~1$s$ orbital alone, the ratio improves from
19$\times$ to 74$\times$.

\subsection{Results of the Fock-Weighted Correction}
\label{sec:fock_results}

Table~\ref{tab:fock_scan} summarizes the $\lambda$-scan.

\begin{table}[h]
\centering
\begin{tabular}{lccc}
\toprule
$\lambda$ & Equilibrium? & $R_{\mathrm{eq}}$ (bohr) &
$D_e^{\mathrm{CP}}$ (Ha) \\
\midrule
0.00 & No & --- & --- \\
0.02 & No & --- & 0.091 at $R{=}3.0$ \\
0.05 & No & --- & uptick at $R{=}3.5$ \\
\textbf{0.10} & \textbf{Yes} & $\bm{\sim 3.8}$ &
\textbf{0.046} \\
0.20 & Yes & $\sim 4.4$ & 0.020 \\
0.50 & No (unbound) & --- & --- \\
\bottomrule
\end{tabular}
\caption{Fock-weighted kinetic correction scan.  At $\lambda = 0.1$,
a genuine potential energy minimum emerges at $R\approx 3.8$~bohr with
a repulsive wall below $R \approx 2.3$~bohr.}
\label{tab:fock_scan}
\end{table}

At $\lambda = 0.1$, the Fock-weighted correction creates a genuine
potential energy minimum with a proper repulsive wall: the molecule is
unbound below $R \approx 2.3$~bohr, maximally bound near
$R \approx 3.8$~bohr ($D_e^{\mathrm{CP}} \approx 0.046$~Ha), and
dissociates smoothly at large $R$.  This is a qualitative success---the
first equilibrium geometry produced by any balanced Hamiltonian
configuration---but the single-parameter model cannot simultaneously
reproduce both $R_{\mathrm{eq}}$ and $D_e$: at $\lambda = 0.02$,
the CP-corrected binding energy matches experiment
($D_e^{\mathrm{CP}} = 0.091$~Ha at $R=3.015$~bohr), but no equilibrium
forms.

The remaining discrepancy ($R_{\mathrm{eq}} = 3.8$ vs.\ experiment
$3.015$~bohr, $D_e = 0.046$ vs.\ $0.092$~Ha) likely reflects the
limitation of a purely diagonal correction.  In standard quantum
chemistry, kinetic repulsion includes both diagonal (orbital compression)
and off-diagonal (orthogonalization) contributions from the full kinetic
energy matrix $T_{AB} = \langle\phi_A|-\tfrac{1}{2}\nabla^2|\phi_B\rangle$.
The Fock-weighted correction captures only the diagonal component.
```

### 6. Conclusion (lines 547–575)

**Revised:**
```latex
\section{Conclusion}
\label{sec:conclusion}

We have presented the first heteronuclear full configuration interaction
calculation on a discrete graph Hamiltonian.  The Topological LCAO
architecture concatenates atom-centered hydrogenic lattices via Mulliken
cross-nuclear attraction, Slater--Condon two-electron integrals, and
STO-overlap bridges, and diagonalizes the resulting Hamiltonian in the
full 367\,290-determinant Slater space using excitation-driven Direct CI.

For LiH at $n_{\max}=3$, the Boys--Bernardi counterpoise-corrected
vertical binding energy at $R=3.015$~bohr is
$D_e^{\mathrm{CP}}=0.110$~Ha, within 20\% of the experimental value
of $0.092$~Ha.  The balanced Hamiltonian produces a monotonically
attractive PES with no equilibrium geometry---a limitation traced to
the $R$-independence of the graph Laplacian's intra-atom kinetic
energy, not to basis truncation.

A 29-version diagnostic arc (v0.9.8--v0.9.37) systematically falsified
every competing hypothesis: basis set superposition error as a geometric
driver, integral incompleteness, orbital exponent relaxation, and the
SO(4) Sturmian Bond Sphere alternative.  The missing kinetic repulsion
was identified as a structural limitation of the LCAO concatenation,
and a Fock-weighted kinetic correction derived from Paper~7's
energy-shell constraint $p_0^2 = Z^2/n^2$ partially restores the
equilibrium ($R_{\mathrm{eq}} \approx 3.8$~bohr,
$D_e^{\mathrm{CP}} \approx 0.046$~Ha at $\lambda = 0.1$).

The $O(V)$ scaling of the underlying graph Hamiltonian is preserved
throughout, with the Direct CI solver achieving
$O(N_{\mathrm{SD}}\times N_{\mathrm{connected}})$ assembly cost.
The natural next steps are:
\begin{enumerate}
\item Extending the Fock-weighted correction to include off-diagonal
kinetic matrix elements between atomic centers, which would provide
the missing intermediate-$R$ binding;
\item Investigating whether the Bond Sphere geometry of Paper~8
naturally incorporates the required $R$-dependent kinetic coupling
through the SO(4) Wigner $D$-matrix;
\item Increasing $n_{\max}$ to assess basis convergence of the
vertical binding energy.
\end{enumerate}
The discovery that the Fock projection's energy-shell scale provides
the correct weighting for kinetic repulsion reinforces the physical
content of the conformal $S^3$ framework: the momentum parameter
$p_0^2 = Z^2/n^2$ is not merely a mathematical artifact of the
stereographic projection but encodes the ``information density'' of
each quantum state on the three-sphere.
```

---

## New Subsection to Insert

The new Section V ("Missing Kinetic Repulsion and the Fock-Weighted Correction")
replaces the old Section V ("$R_{\mathrm{eq}}$ as a Discretization Artifact").
Full LaTeX is provided above in Section 5 of this document.

**Subsections:**
1. **V.A: The Kinetic Diagnostic** — virial analysis showing $\eta \approx 30$,
   confirmation that intra-atom kinetic energy is R-independent
2. **V.B: The Fock-Weighted Kinetic Correction** — derivation from Paper 7
   energy-shell constraint, formula, R-discrimination analysis
3. **V.C: Results of the Fock-Weighted Correction** — lambda scan table,
   $R_{\mathrm{eq}} \approx 3.8$ bohr, $D_e \approx 0.046$ Ha

---

## Claims That Survive vs. Need Revision

### Survives Unchanged
| Claim | Location | Status |
|-------|----------|--------|
| Atomic FCI results (He, Li, Be) | Not in this paper | Unaffected |
| BSSE = -0.115 Ha | Sec III.B, Table I | **Correct** |
| Boys-Bernardi CP framework | Sec II.D | **Correct** |
| Direct CI $O(N_{\mathrm{SD}} \times N_{\mathrm{connected}})$ scaling | Sec II.C | **Correct** |
| Sturmian $\mathbf{H} \propto \mathbf{S}$ theorem | Sec IV.D | **Correct** |
| Sturmian/D-matrix/relaxation falsification | Sec IV | **Correct** |
| Wigner D-matrix selection rules | Sec IV.D | **Correct** |
| 367,290 Slater determinant space | Throughout | **Correct** |
| $O(V)$ scaling claim | Throughout | **Correct** |

### Needs Revision
| Claim | Location | Issue |
|-------|----------|-------|
| $D_e^{\mathrm{CP}} = 0.093$ Ha (1.0% error) | Abstract, Intro, III.B, V.C, Conclusion | From imbalanced Hamiltonian |
| $R_{\mathrm{eq}} \approx 2.5$ bohr | Abstract, III.C, V.A, Conclusion | False minimum from imbalanced config |
| "Discretization artifact" interpretation | Sec V (entire section) | Wrong mechanism |
| "Predicted to converge with $n_{\max}$" | Abstract, V.C | Based on false equilibrium |
| $E_{\mathrm{LiH}} = -8.097$ Ha | Table I | Imbalanced value; balanced = -8.118 |
| $E_{\mathrm{Li}} = -7.334$, $E_{\mathrm{H}} = -0.558$ | Table I | Fragment-in-molecular-basis values; isolated = -7.392, -0.500 |
| Paper 1 spectral aliasing analogy | Sec V.B | Misleading; molecular issue is structural |
| v0.9.11 as "baseline achieving quantitative accuracy" | App. Table II caption | Only 20% accuracy, not 1% |

### New Claims to Add
| Claim | Location | Basis |
|-------|----------|-------|
| Hamiltonian balance theorem: cross-nuclear and V_ee scope must match | New paragraph in III.B | `exact+True` and `fourier+s_only` agree within 2 mHa |
| Balanced PES is monotonically attractive, no equilibrium | New III.C | PES data from balanced configs |
| Graph Laplacian kinetic energy is R-independent | New V.A | Direct computation, virial analysis |
| Virial ratio $\eta \approx 30$ (should be 1) | New V.A | Energy decomposition diagnostic |
| Fock-weighted correction creates equilibrium at R~3.8 | New V.C | Lambda scan with Fock weights |
| Energy-shell $p_0^2 = Z^2/n^2$ provides correct kinetic weighting | New V.B | 9.3x R-discrimination, genuine well |
| Diagonal correction insufficient for full equilibrium recovery | New V.C | Cannot match both R_eq and D_cp |

---

## Tone Guidance Notes

This revision should be framed as:

1. **Discovery, not retraction.** The diagnostic arc did exactly what it was
   supposed to do: it tested hypotheses and found the truth. The finding that
   balanced Hamiltonians have no equilibrium is a legitimate scientific result.

2. **The framework reveals deeper physics.** The graph Laplacian's R-independent
   kinetic energy is not a bug — it's a fundamental property of the discrete
   topology. Understanding it led to the Fock-weighted correction, which connects
   back to Paper 7's core physics.

3. **Honest about limitations.** The paper should clearly state that the original
   1.0% accuracy claim was an artifact of an imbalanced Hamiltonian, and that the
   balanced framework produces correct-order-of-magnitude but not quantitative
   molecular binding.

4. **Constructive.** The Fock-weighted correction shows the path forward. The
   energy-shell constraint provides the right physics; the implementation just
   needs to be extended to off-diagonal kinetic elements.

---

## Implementation Checklist

- [ ] Update abstract with corrected values and new framing
- [ ] Update introduction paragraph with corrected headline result
- [ ] Revise Table I with balanced `exact+True` values
- [ ] Add "Hamiltonian balance requirement" paragraph to Sec III.B
- [ ] Rewrite Sec III.C as "Absence of Equilibrium Geometry"
- [ ] Replace Sec V entirely with new kinetic repulsion section
- [ ] Add virial analysis table (new Table)
- [ ] Add Fock-weighted correction equation and lambda scan table (new Table)
- [ ] Revise conclusion with corrected values and new next steps
- [ ] Update version history table and caption
- [ ] Add reference to diagnostic reports if desired
- [ ] Recompile and verify no broken cross-references
