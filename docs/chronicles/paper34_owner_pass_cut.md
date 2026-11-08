# Paper 34 owner pass (2026-11-07): material cut, verbatim

Owner pass on Paper 34 (`papers/group6_precision_observations/paper_34_projection_taxonomy.tex`)
for three PI-approved items (2026-11-07), with its citers (Paper 31, Paper 35, the group6
synthesis, `docs/topic_to_paper_lookup.md`).  Every removed line is preserved verbatim below
as the `-` lines of a zero-context unified diff against commit d0d110b8 (the `+` lines are the
replacements, kept for orientation).

## What was cut or replaced

| Item | Material | Disposition |
|:-----|:---------|:------------|
| 1 | Projection 13 name "Drake--Swainson asymptotic subtraction"; its Source/Target text; the "combined hydrogen Lamb shift via this projection ... $-0.39\%$ at the sweet-spot $N=40$" sentence; the base-unit row "$K$ (cancels), $[E]$, transient scale"; the $[E]$ count "four" | Renamed "Drake--Swainson normalization of Bethe-log sums" (division by $2Z^4/n^3$ for every $\ell$, $(2\ell+1)$ bookkeeping, Paper 36 Sec. 5); the K-split is cut from Paper 36, so the transient $[E]$ scale and the sweet-spot sentence are cut; $[E]$ now three |
| 2 | Gauge-choice "second independence witness"; "the transcendental-only corner is the gauge sector"; "$1/(4\pi)$ per loop" signatures (dictionary table, Lorentzian audit row 26, duality item, vd-correlation remark, Obs. axes_correlated, open questions, conclusion); "Feynman-parameter denominators of the Lorenz-gauge computation" | Vector-photon $1/(4\pi)$ recorded as a vertex-model normalization (Paper 33); gauge choice is not a witness (no transcendental measured); corner witnessed by the $k=0$ volume-slot factors of the Hopf bundle base sphere and the Wick-rotation time circle |
| 3 | Hydrogen/He$^+$ Lamb rows tagged Fock $\circ$ Sturmian $\circ$ spectral action $\circ$ Drake--Swainson and "framework-native"; the LS-4 ($1054$ MHz, $-0.39\%$) and LS-3 ($1029$ MHz, $-2.77\%$) Lamb rows of the matches catalogue; the LS-4 and LS-3 Lamb rows of the off-precision table; the Lamb shift in Theorem 1's tested rows ("eight" -> "seven"), in the depth-4 falsification, in the L2-count = 0 list and the open-question datapoint, and in Remark signature != value; muonium Lamb comparison "vs Mariam 1982 $1054(22)$ MHz" (Mariam et al. PRL 49, 993 (1982) is a ground-state hyperfine measurement); "Eides Tab. 7.4 / Karshenboim 2005: $\Delta_{\rm el-VP}\approx+1.5$ meV" (Eides Phys. Rep. 342 numbers tables in roman and treats the muonic-hydrogen 2S, not 1S, hyperfine structure); Paper 35 "Eides--Grotch--Shelyuto 2001 Table 7.4 cumulative central value $+7.10$ MHz"; the 3D Bethe-log residual "$-0.24\%$" (Paper 36: $+0.07\%$) | Retagged to the textbook assembly of Paper 36; comparison moved to Mu-MASS (Ohayon et al., PRL 128, 011802 (2022), $1047.2(2.3)(1.1)$ MHz); el-VP size unsourced; Paper 35 now cites Table III, $+7.13$ MHz, via Paper 36 |

## Verbatim diff

```diff
diff --git a/docs/topic_to_paper_lookup.md b/docs/topic_to_paper_lookup.md
index 1fa84bd9..eff6390f 100644
--- a/docs/topic_to_paper_lookup.md
+++ b/docs/topic_to_paper_lookup.md
@@ -376 +376 @@
-| Drake-Swainson asymptotic subtraction (13th projection) | 34 | Sec III.13 | Observation |
+| Drake-Swainson normalization of Bethe-log sums (13th projection) | 34 | Sec III.13 | Observation |
diff --git a/papers/group3_foundations/paper_31_universal_coulomb_partition.tex b/papers/group3_foundations/paper_31_universal_coulomb_partition.tex
index 68d1056d..f46825d3 100644
--- a/papers/group3_foundations/paper_31_universal_coulomb_partition.tex
+++ b/papers/group3_foundations/paper_31_universal_coulomb_partition.tex
@@ -908,6 +908,6 @@ The two charter members of this central tier are:
-\item \textbf{Drake--Swainson asymptotic subtraction}
-  (Paper~34 projection~13).  The intermediate
-  cutoff $K$ in the Bethe-log split-and-cancel scheme is a
-  transient additive shift of the spectral integration measure
-  along $\mathbb{R} \cdot \mathbb{1}$ in the spectral-density
-  domain; the shift cancels in the converged result.
+\item \textbf{Drake--Swainson normalization of Bethe-log sums}
+  (Paper~34 projection~13).  The Bethe-log spectral sum is
+  divided by the fixed rational $2Z^4/n^3$ for every $l$ (with a
+  factor $2l+1$ when the initial magnetic substates are summed;
+  Paper~36~\S~5): a scalar rescaling, central by construction,
+  that introduces no variable.
@@ -924 +924 @@ therefore three-tier rather than two-tier:
-  shifts, flow-tier subtractions.  Universal across central
+  shifts, Bethe-log normalizations.  Universal across central
diff --git a/papers/group6_precision_observations/paper_34_projection_taxonomy.tex b/papers/group6_precision_observations/paper_34_projection_taxonomy.tex
index 374f2736..8a99b5eb 100644
--- a/papers/group6_precision_observations/paper_34_projection_taxonomy.tex
+++ b/papers/group6_precision_observations/paper_34_projection_taxonomy.tex
@@ -80,6 +80,2 @@ form on the catalogue itself (depth-3 chains exist at residual~$0$
-and at $-211$~ppm; depth-4 chains exist at $+2$~ppm and at
-$-0.534\%$).  The Lamb shift example (depth-4 chain, $-0.534\%$
-residual) is the leading one-loop assembly with tabulated Bethe
-logarithms, truncated before the one-loop $\alpha(Z\alpha)^5$ term
-(Paper~36), so its residual measures an order truncation, not a basis
-limit; the H~1S polarizability example
+and at $-211$~ppm; a depth-4 chain exists at $+2$~ppm).  The
+H~1S polarizability example
@@ -88 +84,3 @@ zero-Layer-2-input class for rational-analytic-answer
-observables.
+observables.  The hydrogen Lamb shift is outside the catalogue's
+projection chains: it is a textbook one-loop assembly with tabulated
+Bethe logarithms (Paper~36).
@@ -638,4 +636,6 @@ Paper~33; the physical vertex is Clifford multiplication by a coclosed
-\textbf{Transcendental signature:} $1/(4\pi)$, the normalization of
-spherical harmonics on the unit $S^2$ that the model's angular vertex
-carries (Paper~33: $\Sigma_{l=1} = 1/(4\pi)$ at $n_{\max} = 2$); the
-physical vertex, normalized to volume one on $S^3$, does not carry it.\\
+\textbf{Transcendental signature:} $1/(4\pi)$, a vertex-model
+normalization: the normalization of spherical harmonics on the unit
+$S^2$ that the model's angular vertex carries (Paper~33:
+$\Sigma_{l=1} = 1/(4\pi)$ at $n_{\max} = 2$).  The physical vertex,
+normalized to volume one on $S^3$, does not carry it, so no
+transcendental of the physical promotion is identified.\\
@@ -660 +660 @@ all H$_2$/HeH$^+$/H$_2$O calculations.
-\subsection{Drake--Swainson asymptotic subtraction}
+\subsection{Drake--Swainson normalization of Bethe-log sums}
@@ -662,4 +662,3 @@ all H$_2$/HeH$^+$/H$_2$O calculations.
-\textbf{Source:} finite-basis Coulomb-Sturmian pseudostate sum
-$J_v(nl)$ for the Bethe logarithm of an $\ell > 0$ state, where the
-velocity-form sum rule $I_v(nl) = (2Z^4/n^3)\,\delta_{\ell,0}$ vanishes
-for $\ell > 0$.\\
+\textbf{Source:} the Bethe-logarithm spectral sum (the numerator of
+the definition in Paper~36~\S~5.1), evaluated as a finite-basis Coulomb-Sturmian
+pseudostate sum $J_v(n,\ell)$.\\
@@ -668 +667,6 @@ $\ln k_0(n,\ell) = J_v(n,\ell) / D_{\text{drake}}(n,\ell)$ in the uniform
-normalization of the Drake--Swainson table, with
+normalization of the Drake--Swainson table: the sum is divided by
+$2Z^4/n^3$ for every $\ell$.  For $s$-states this denominator equals the
+velocity-form sum rule; for $\ell > 0$ that sum rule vanishes while the
+numerator stays finite, and the same denominator is kept.  The code also
+sums the numerator over initial magnetic substates, which multiplies it
+by $2\ell+1$, so it divides by
@@ -672,3 +676,2 @@ normalization of the Drake--Swainson table, with
-the factor $2\ell+1$ is the code's sum over initial magnetic substates
-(Paper~36~\S~5.1).  For $\ell > 0$ the velocity-form numerator converges
-directly in the basis, so there is no $0/0$.\\
+the factor $2\ell+1$ is that bookkeeping and nothing more
+(Paper~36~\S~5.1).\\
@@ -683,6 +686,3 @@ logarithms $\ln |\Delta E_m / \mathrm{Ry}|$.\\
-Drake--Swainson 1990 tabulated values).  The combined hydrogen Lamb shift via this projection
-chained with Sturmian and spectral-action gives $-0.39\%$ residual at
-the sweet-spot $N=40$ (transient cancellation; converges to LS-1
-baseline $-3.10\%$ at $N \to \infty$ --- the pre-LS-6a convention,
-superseded by the LS-6a fix at $-0.534\%$, whose residual is not
-attributed; Paper~36).
+Drake--Swainson 1990 tabulated values).  The hydrogen Lamb shift of
+Paper~36 does not use these computed values: it is assembled with the
+tabulated Bethe logarithms.
@@ -798,3 +798,3 @@ follow-up beyond sprint scale.
-\textbf{Source:} framework-native single-particle Eides~\S~3.2
-self-energy / Lamb-shift bracket on the Sturmian-reparameterized
-Fock graph, valid in the $\lambda_l \ll \lambda_n$ hierarchy where the
+\textbf{Source:} the textbook single-particle one-loop
+self-energy / Lamb-shift bracket of Eides \emph{et al.}\ (Paper~36),
+valid in the $\lambda_l \ll \lambda_n$ hierarchy where the
@@ -1130,2 +1130,2 @@ sense Phillips--Kleinman is in the same dimensionless class as
-Drake--Swainson asymptotic subtraction (\S\ref{sec:proj_drake_swainson})
---- both are bookkeeping projections that re-organize an already-
+the Drake--Swainson normalization of Bethe-log sums
+(\S\ref{sec:proj_drake_swainson}) --- both are bookkeeping projections that re-organize an already-
@@ -2324,20 +2324,15 @@ data and does not move the chain across the
-This is the \emph{second independence witness} in the
-(no-variable, no-dimension, transcendental-only) corner of the
-variable/dimension/transcendental axis grid.  Until Sprint~2 the only
-witness in that corner was \S~\ref{sec:proj_vector_photon} (vector-photon
-promotion, which introduces $1/(4\pi)$ without any new variable or
-dimension); the gauge-choice projection is its structural sibling.
-Both witnesses operate on the same object (the framework's $U(1)$
-connection) at successive structural steps: vector-photon promotion
-transforms the scalar $1$-cochain into a vector harmonic with
-$(q, m_q)$ labels, and the present projection fixes a specific gauge
-representative of the resulting connection before propagators are
-computed.  The fact that the two witnesses are \emph{both
-gauge-related} sharpens the V/D-correlation observation of
-\S~\ref{sec:dictionary} and \S~\ref{sec:axiomatization}: the
-transcendental-only corner of the axis grid is the gauge sector of the
-framework.  It is not an accident that the only two projections
-preserving both V and D while moving the transcendental axis are
-gauge-construction steps; the structural content of that corner is
-``add transcendental content to the gauge sector without adding any
-physical variable or dimension.''
+The projection is \emph{not} a witness of the (no-variable,
+no-dimension, transcendental-only) corner of the
+variable/dimension/transcendental axis grid.  It preserves V and D, but
+the transcendental axis has not been seen to move: the one factor
+identified in Coulomb gauge, $1/(4\pi)$, is the vertex model's
+unit-$S^2$ normalization, which the physical vertex does not carry
+(Paper~33), and no Lorenz- or Feynman-gauge content has been computed.
+Vector-photon promotion (\S~\ref{sec:proj_vector_photon}) carries the
+same $1/(4\pi)$ and the same caveat, so the corner is not the gauge
+sector of the framework.  The corner's occupants in the dictionary are
+volume factors of the $k = 0$ volume slot: the base-sphere factor
+$\pi = \mathrm{Vol}(S^2)/4$ of the Hopf bundle (\S~\ref{sec:proj_hopf}) and the
+$2\pi = \mathrm{Vol}(S^1)$ of the Euclidean time circle under Wick
+rotation (\S~\ref{sec:proj_wick_rotation}, counting its binary signature
+index, like a gauge selection, as no physical variable).
@@ -2346 +2341 @@ Categorically distinct from \S~\ref{sec:proj_vector_photon} despite
-the shared axis pattern: vector-photon promotion is the structural
+the shared V/D pattern: vector-photon promotion is the structural
@@ -2367 +2362 @@ input.
-Sibling in the calibration tier of Paper~18: the $1/(4\pi)$
+Relation to the calibration tier of Paper~18: the $1/(4\pi)$
@@ -2369,13 +2364,10 @@ of the vertex model is the normalization of spherical harmonics on the
-unit $S^2$ (Paper~33), which the physical vertex does not carry.
-Under the master-Mellin reading of Paper~32~\S~VIII and
-Paper~18~\S~III.7, this calibration $\pi$ is the M1 sub-mechanism
-signature ($\mathrm{Vol}(S^2)/4$, the $k = 0$ volume slot).  A
-Lorenz-gauge calculation, if performed, would produce different
-intermediate $\pi$-content (retardation kernels rather than static
-Hopf-base normalisations) but would necessarily land on the same
-physical observable.  The gauge-choice projection therefore sits at
-the interface of the framework's gauge construction
-(\S~\ref{sec:proj_vector_photon} +
-\S~\ref{sec:proj_wilson}) and the framework's transcendental taxonomy
-(Paper~18), and it is the structural step at which a Layer-2 chain
-commits to a specific way of accumulating $\pi$ at each loop.
+unit $S^2$ (Paper~33), which the physical vertex does not carry, so it
+is not a per-loop factor of Coulomb gauge.  A Lorenz-gauge calculation,
+if performed, would produce different intermediate $\pi$-content
+(retardation kernels) but would necessarily land on the same physical
+observable.  Whether the choice of gauge changes the intermediate
+transcendental content of a framework computation is therefore open;
+the projection sits at the interface of the framework's gauge
+construction (\S~\ref{sec:proj_vector_photon} +
+\S~\ref{sec:proj_wilson}) and the transcendental taxonomy (Paper~18)
+without, so far, any measured transcendental of its own.
@@ -2392,4 +2384,4 @@ content that the Coulomb-gauge construction is silently encoding;
-if it succeeds, the gauge-choice projection is verified as a
-genuine independence witness in the sense of
-Observation~\ref{obs:axes_correlated} rather than a tautological
-relabelling.
+if it succeeds and the intermediate transcendental content differs
+between the gauges, the gauge-choice projection becomes a witness of
+the transcendental-only corner in the sense of
+Observation~\ref{obs:axes_correlated} rather than a relabelling.
@@ -3165 +3157,2 @@ Wilson plaquette & $\beta = 1/g^2$ & dimensionless & SU(2) Haar \\
-Vector-photon promotion & none & dimensionless & $1/(4\pi)$ per loop \\
+Vector-photon promotion & none & dimensionless & $1/(4\pi)$, a
+vertex-model normalization (absent from the physical vertex) \\
@@ -3167,2 +3160,3 @@ Molecule-frame hypersph. & $R$ & length & piecewise-smooth in $R$ \\
-Drake--Swainson asymptotic subtraction & none & dimensionless &
-flow tier; $D_{\text{drake}} = 2(2\ell+1)Z^4/n^3$ rational \\
+Drake--Swainson normalization of Bethe-log sums & none & dimensionless &
+rational normalization $2Z^4/n^3$ (code: $2(2\ell+1)Z^4/n^3$); flow-tier
+logarithms sit in the sum it normalizes \\
@@ -3206,3 +3200,3 @@ none (selection of representative) & preserves &
-gauge-dependent at intermediate-propagator level
-($1/(4\pi)$ per loop in Coulomb gauge = M1 volume slot);
-gauge-invariant on physical observables (Ward identity) \\
+gauge-dependent at intermediate-propagator level in QED; none measured
+in the framework (the Coulomb-gauge $1/(4\pi)$ is the vertex model's
+normalization); gauge-invariant on physical observables (Ward identity) \\
@@ -3226,2 +3220,2 @@ axes proposed by this paper.  The thirteenth row (Drake--Swainson
-asymptotic subtraction) was added in sprint LS-4 for the $\ell > 0$
-Bethe logarithms; the fourteenth and
+normalization of Bethe-log sums) was added in sprint LS-4 for the
+$\ell > 0$ Bethe logarithms; the fourteenth and
@@ -3250,3 +3244,2 @@ identity / state-side reduction) as axis-completeness additions:
-gauge choice is the second independence witness in the
-transcendental-only corner alongside vector-photon promotion, Wick
-rotation promotes the Bisognano--Wichmann reading from a \S~VIII
+gauge choice preserves V and D with no transcendental measured in the
+framework, Wick rotation promotes the Bisognano--Wichmann reading from a \S~VIII
@@ -3268 +3261,2 @@ tier:
-Hopf bundle / vector-photon promotion projections.
+Hopf bundle projections (vector-photon promotion carries only the
+vertex model's $1/(4\pi)$ normalization).
@@ -3279,3 +3273,4 @@ Wigner $3j$ / $D$ alone (no projection beyond Layer~1).
-\item Flow tier (Paper~18 \S~IV) $\leftrightarrow$ Drake--Swainson
-asymptotic subtraction.  The flow content sits in the logarithms
-$\ln|\Delta E_m/\mathrm{Ry}|$ of the Bethe-logarithm spectral sum.
+\item Flow tier (Paper~18 \S~IV) $\leftrightarrow$ the Bethe-logarithm
+spectral sum that the Drake--Swainson normalization divides.  The flow
+content sits in the logarithms $\ln|\Delta E_m/\mathrm{Ry}|$ of that
+sum; the normalization itself is rational.
@@ -3294,2 +3289,2 @@ spectral and temporal observables catalogued here traces to one of (a)~an M1 ($k
-$\mathrm{Vol}(S^3) = 2\pi^2$, or the vector-photon $1/(4\pi)$
-unit-$S^2$ Weyl normalisation), (b)~spectral-action even-zeta content
+$\mathrm{Vol}(S^3) = 2\pi^2$, or the vertex model's $1/(4\pi)$
+unit-$S^2$ normalization), (b)~spectral-action even-zeta content
@@ -3302,2 +3297,3 @@ enters through embedding-tier routes outside this list (e.g.\ the
-Feynman-parameter denominators of the Lorenz-gauge computation, and the
-genus-one elliptic periods of the embedding tier).  The
+genus-one elliptic periods of the embedding tier; a covariant-gauge
+computation would add Feynman-parameter denominators, but none has
+been performed in the framework).  The
@@ -3344,9 +3340,12 @@ further strengthen the V/D-correlation observation: gauge choice and
-Wick rotation both preserve V and D while moving the transcendental
-axis (joining vector-photon promotion in the (no-variable,
-no-dimension, transcendental-only) corner --- giving \emph{three}
-independent witnesses there, all gauge-related or signature-related;
-the structural content of that corner is the gauge-and-signature
-sector of the framework), and apparatus identity introduces $\rho$ as
-the first state-side variable while preserving dimension.  The
-transcendental axis remains the only genuinely independent axis; V/D
-are still perfectly correlated.  The three intra-nuclear projections of
+Wick rotation both preserve V and D, and apparatus identity introduces
+$\rho$ as the first state-side variable while preserving dimension.
+Of these, only Wick rotation moves the transcendental axis
+($2\pi = \mathrm{Vol}(S^1)$); gauge choice has no transcendental
+measured in the framework (\S~\ref{sec:proj_gauge_choice}).  The
+(no-variable, no-dimension, transcendental-only) corner is witnessed by
+volume factors of the $k = 0$ volume slot --- the base-sphere factor
+$\mathrm{Vol}(S^2)/4$ of the Hopf bundle and the Wick-rotation time
+circle --- not by the gauge sector: the
+$1/(4\pi)$ of vector-photon promotion is a vertex-model normalization
+(Paper~33).  The transcendental axis remains the only genuinely
+independent axis; V/D are still perfectly correlated.  The three intra-nuclear projections of
@@ -3357,3 +3356,3 @@ its associated dimension ($[L]$, $[L]$, $[L]^2$ respectively).  The transcendent
-independence witness --- vector-photon promotion (\S~\ref{sec:layer2},
-projection~11) introduces $1/(4\pi)$ without introducing any new
-variable or dimension.  This is the structural content of
+independence witness --- the Hopf bundle (\S~\ref{sec:proj_hopf},
+projection~2) introduces $\pi = \mathrm{Vol}(S^2)/4$ without
+introducing any new variable or dimension.  This is the structural content of
@@ -3421,3 +3420,2 @@ scale-extension \\
-Drake--Swainson & $K$ (cancels) & $[E]$ & transient $[E]$ scale
-(cancels in $\beta_\text{low}(K) + \beta_\text{high}(K)$) &
-scale-extension (transient) \\
+Drake--Swainson normalization & none & dim'less & none (rational
+normalization $2Z^4/n^3$ of a dimensionless sum) & pure \\
@@ -3496 +3494 @@ and requires $l \geq 1$ on the electronic side.
-\item $[E]$: \emph{four} projections.  Fock or Bargmann--Segal
+\item $[E]$: \emph{three} projections.  Fock or Bargmann--Segal
@@ -3498,3 +3496 @@ anchors the foundational energy (Ry or $\hbar\omega$); spectral
-action introduces UV cutoff $\Lambda$; Drake--Swainson introduces a
-transient subtraction scale $K$ (which cancels in the regularized
-result); two-body Breit retardation introduces $\alpha^4 \cdot
+action introduces UV cutoff $\Lambda$; two-body Breit retardation introduces $\alpha^4 \cdot
@@ -3524 +3520 @@ admits \emph{six} scalar instances plus one rank-2 tensor instance,
-while $[E]$ admits four.  The framework records more length scales
+while $[E]$ admits three.  The framework records more length scales
@@ -3527 +3523 @@ distribution scales) than energy scales (Rydberg, UV cutoff,
-multi-loop subtraction, two-body Breit correction).  The structural
+two-body Breit correction).  The structural
@@ -3532 +3528 @@ selection rules; spectral content has fewer distinct roles
-(foundational, regulator, subtraction, correction).  The conjugate
+(foundational, regulator, correction).  The conjugate
@@ -3614,2 +3610,2 @@ Remark~\ref{rem:vd_correlation}: V--D correlation holds at the
-compound dimension level, but at the base-unit level $[L]$ and
-$[E]$ each admit four instances while $[M]$, $[T]$ admit one each
+compound dimension level, but at the base-unit level $[L]$ admits
+six scalar instances and $[E]$ three, while $[M]$, $[T]$ admit one each
@@ -3786,3 +3782,4 @@ $-0.25$ ppm vs Mu-MASS \\
-Mu 2$S_{1/2}$-2$P_{1/2}$ Lamb shift $1047.63$ MHz (framework one-loop
-QED with full Uehling, sprint precision-catalogue 2026-05-08) &
-Fock $\circ$ spectral action $\circ$ Sturmian $\circ$ rest-mass projection
+Mu 2$S_{1/2}$-2$P_{1/2}$ Lamb shift $1047.63$ MHz (the Paper~36
+leading one-loop assembly with full Uehling, at the muonium reduced
+mass; sprint precision-catalogue 2026-05-08) &
+textbook assembly (tabulated Bethe logarithms) $\circ$ rest-mass projection
@@ -3790 +3787 @@ Fock $\circ$ spectral action $\circ$ Sturmian $\circ$ rest-mass projection
-$-0.6\%$ vs Mariam 1982 experiment $1054(22)$ MHz (a theory comparison is withdrawn) \\
+$+0.04\%$ vs Mu-MASS (Ohayon~2022) $1047.2(2.3)(1.1)$ MHz; not attributed \\
@@ -3847,6 +3844 @@ Bethe log $\ln k_0(3D) = -0.005236$ (LS-4) & Drake--Swainson $\circ$ Sturmian
-$\circ$ Fock & $Z$ & dimensionless & flow & $-0.24\%$ at $N=40$ \\
-Lamb shift $1054$ MHz (LS-4 N=40 sweet spot) & Fock $\circ$ spectral action
-$\circ$ Sturmian $\circ$ Drake--Swainson & $Z, \alpha$ & energy & calibration
-$+$ flow & $-0.39\%$ (T$+$A cancellation) \\
-Lamb shift $1029$ MHz (LS-3) & Fock $\circ$ spectral action $\circ$
-Sturmian (acc.) & $Z, \alpha$ & energy & calibration & $-2.77\%$ (T$+$A)\\
+$\circ$ Fock & $Z$ & dimensionless & flow & $+0.07\%$ at $N=40$ \\
@@ -3920,4 +3912,3 @@ hydrogen Bohr level $E_n = -Z^2/(2n^2)$ (Fock alone, signature
-$\{$energy, no~$\pi\}$) and the hydrogen Lamb shift (Fock $+$
-spectral action $+$ Sturmian $+$ Drake--Swainson, signature
-$\{$energy, $\pi\}$) both have dimension $\{$energy$\}$, yet differ
-by 13 orders of magnitude.  The Layer-2-presence bound of
+$\{$energy, no~$\pi\}$) and the $2S$ Uehling shift (entered $\Pi$
+$\circ$ Fock, signature $\{$energy, $\pi\}$) both have dimension
+$\{$energy$\}$, yet differ by seven to eight orders of magnitude.  The Layer-2-presence bound of
@@ -4045,34 +4036,15 @@ $\ell>0$ case; sprint LS-4. \\
-Lamb shift $1057.85$ MHz & Fock $\circ$ spectral action $\circ$
-Sturmian (acc.) $\circ$ (Drake--Swainson, l=1) & $-0.39\%$ at LS-4
-sweet spot $N=40$ & T cancels A & At $N=40$, $T_{2S} = +28.8$~MHz
-truncation cancels $A = -32.78$~MHz one-loop ceiling; converges to
-LS-1 baseline $-3.10\%$ at $N\to\infty$ (pre-LS-6a convention;
-superseded by LS-6a at $-0.534\%$).  Sprint LS-4. \\
-Lamb shift $1057.85$ MHz (LS-3 baseline) & Fock $\circ$ spectral action
-$\circ$ Sturmian (acc.) & $-2.77\%$ (LS-3 prediction $1029$ MHz) &
-\textbf{C} (resolved by LS-6a) $+$ A & Sprint LS-6a (May 2026)
-re-derived the same one-loop physics in the Eides \S~3.2 canonical
-convention and found the LS-1 implementation had inadvertently
-subtracted the Uehling kernel constant from the canonical $10/9$ SE
-coefficient (Uehling already enters separately as a VP contribution).
-The convention fix supersedes this row: see the LS-6a row below for
-the clean result. The $\sim +32.78$~MHz residual that had been
-classified A is structurally a $\sim +27.13$~MHz one-loop convention
-correction (\textbf{C} resolved), leaving a $\sim +5.65$~MHz residual
-that is not attributed (Paper~36). \\
-Lamb shift $1057.85$ MHz (LS-6a, Eides \S~3.2 convention) & Fock
-$\circ$ spectral action $\circ$ Sturmian (acc.) $\circ$ (Drake--Swainson,
-$\ell=1$) & $-0.534\%$ (LS-6a prediction $1052.19$~MHz) & residual not
-attributed (one-loop $\alpha(Z\alpha)^5$ omitted) & Eides \S~3.2 canonical one-loop
-SE convention.  Convention shift from LS-1: $+27.13$~MHz, exactly
-$(4/15) \cdot \alpha^3 Z^4 / (\pi n^3) \cdot \text{Ha-to-MHz}$ at
-$n=2$, $Z=1$ (the Uehling kernel constant $10/9 - 38/45 = 4/15$).
-Components: SE 2$S_{1/2}$ $+1066.44$~MHz (Bethe $+953.40$~MHz +
-canonical $10/9$ piece $+113.04$~MHz); VP 2$S_{1/2}$ $-27.13$~MHz;
-SE 2$P_{1/2}$ $-12.88$~MHz; VP 2$P_{1/2}$ $0$.  $5.8\times$
-reduction in absolute error vs the LS-1 one-loop ceiling
-($-32.78$~MHz); vs LS-3 ($-28.85$~MHz) the reduction is $5.1\times$.
-(Baseline label corrected 2026-08-28: the factor $5.8$ is the ratio
-to LS-1, not to LS-3.)  The residual is not attributed.  An LS-7 decomposition into about $+1.20$~MHz of multi-loop QED and about $+4.4$~MHz of non-loop physics is withdrawn (/qa group5, 2026-11-07): its itemized entries reached only about $70\%$ of the residual, it omitted the one-loop $\alpha(Z\alpha)^5$ term ($A_{50} \approx 9.29$, about $+6.9$~MHz on 2S, by itself larger than the residual), and its ``hyperfine averaging'' row had no source, since the measured interval is defined between hyperfine centroids (Paper~36).  Sprint provenance:
-\texttt{debug/ls6a\_eides\_convention\_memo.md},
-\texttt{debug/ls7\_two\_loop\_se\_memo.md}. \\
+Lamb shift $1057.85$ MHz (leading one-loop assembly, Paper~36) &
+textbook assembly: one-loop expressions of Eides \emph{et al.}\ with
+tabulated Drake--Swainson Bethe logarithms; no GeoVac projection enters
+& $-0.534\%$ (assembly $1052.19$~MHz) & residual not
+attributed (one-loop $\alpha(Z\alpha)^5$ omitted) & Components:
+SE 2$S_{1/2}$ $+1066.44$~MHz (Bethe $+953.40$~MHz + canonical $10/9$
+piece $+113.04$~MHz); VP 2$S_{1/2}$ $-27.13$~MHz; SE 2$P_{1/2}$
+$-12.88$~MHz; VP 2$P_{1/2}$ $0$.  An earlier assembly used the
+constant $38/45$ in place of $10/9$ and double-counted the Uehling
+term; the correction is exactly $(4/15)\,\alpha^3 Z^4/(\pi n^3)$~Ha
+$= +27.13$~MHz (Paper~36).  The residual is not attributed: the
+one-loop $\alpha(Z\alpha)^5$ term ($A_{50} \approx 9.29$, about
+$+6.9$~MHz on 2S) is by itself larger than it, and an earlier
+decomposition into multi-loop and non-loop pieces is withdrawn
+(Paper~36). \\
@@ -4157,3 +4129,3 @@ Layer-2 input $r_Z = 1.045$~fm & $-7689$ ppm vs Krauth full theory
-& Dominant residual is electron-VP (Uehling) in the muonic potential,
-the leading muonic QED correction (Eides Tab.~7.4 / Karshenboim 2005:
-$\Delta_\text{el-VP} \approx +1.5$ meV).  Same LS-8a wall as HF-5; in
+& Dominant residual is expected to be electron-VP (Uehling) in the
+muonic potential, the leading muonic QED correction (the size of that
+correction is not sourced here).  Same LS-8a wall as HF-5; in
@@ -4288,5 +4260,8 @@ Compton wavelength.  Sprint MH Track A, 2026-05-08. \\
-Mu 2$S_{1/2}$-2$P_{1/2}$ Lamb $1047.63$ MHz vs Mariam~1982 experiment
-$1054(22)$ MHz & Fock $\circ$ spectral action $\circ$ Sturmian $\circ$
-rest-mass projection & $-0.6\%$ (within $\pm 22$ MHz exp uncertainty)
-& A & $-0.6\%$ residual against Mariam, within the 1982-era
-experimental uncertainty.  A comparison with a muonium theory value
+Mu 2$S_{1/2}$-2$P_{1/2}$ Lamb $1047.63$ MHz vs Mu-MASS experiment
+$1047.2(2.3)(1.1)$ MHz (Ohayon \emph{et al.}, PRL \textbf{128}, 011802
+(2022)) & textbook assembly (Paper~36, tabulated Bethe logarithms)
+$\circ$ rest-mass projection & $+0.04\%$ (within the experimental
+uncertainty) & A & The leading one-loop assembly at the muonium reduced
+mass; it omits the one-loop $\alpha(Z\alpha)^5$ term and recoil, each
+larger than the $+0.4$~MHz difference, so the agreement is not
+attributed.  A comparison with a muonium theory value
@@ -4298,4 +4273,4 @@ He$^+$ 2$S_{1/2}$-2$P_{1/2}$ Lamb $14041.13(17)$ MHz (van
-Wijngaarden--Holuj--Drake 2000) & Fock $\circ$ spectral action $\circ$
-Sturmian $\circ$ (Drake--Swainson, $\ell=1$) $\circ$
-\S\ref{sec:proj_charge_density} & $-1.29\%$ (framework $+$ Drake-Yan
-Layer-2 lit total $13\,859.58$ MHz; framework-native $13\,841.11$ MHz at
+Wijngaarden--Holuj--Drake 2000) & textbook assembly (Paper~36,
+tabulated Bethe logarithms) $\circ$
+\S\ref{sec:proj_charge_density} & $-1.29\%$ (assembly $+$ Drake-Yan
+Layer-2 lit total $13\,859.58$ MHz; assembly alone $13\,841.11$ MHz at
@@ -4303,2 +4278,2 @@ $-1.42\%$) & A & \textbf{Z-scaling test of the Paper~36
-one-loop assembly (Z=1 $\to$ Z=2)}.  Framework's projection-chain machinery $Z$-portable:
-three Z-scaling structural tests pass (Lamb$_\text{SE}$ ratio $13.21$
+one-loop assembly (Z=1 $\to$ Z=2)}.  The implementation is $Z$-portable:
+three Z-scaling tests pass (Lamb$_\text{SE}$ ratio $13.21$
@@ -4308 +4283 @@ match to $0.12\%$ against $Z^4 \times (r_E({}^4\text{He})/r_p)^2 = 63.52$).  The
-$-1.42\%$ framework-native residual is structurally the Z-amplified analog
+$-1.42\%$ assembly residual is structurally the Z-amplified analog
@@ -4832,3 +4807 @@ SE 2$S_{1/2}$ (Bethe + Eides 10/9) & $+1066.44$ &
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_sturmian} $\circ$
-\S\ref{sec:proj_spectral_action} $\circ$ \S\ref{sec:proj_drake_swainson}
-& FN \\
+textbook (Eides one-loop; tabulated Bethe log) & TA \\
@@ -4836 +4809 @@ VP 2$S_{1/2}$ (Uehling) & $-27.13$ &
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_spectral_action} & FN \\
+textbook (Uehling $-4/15$) & TA \\
@@ -4838,3 +4811 @@ SE 2$P_{1/2}$ & $+12.88$ &
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_sturmian} $\circ$
-\S\ref{sec:proj_spectral_action} $\circ$ \S\ref{sec:proj_spinor}
-$\circ$ \S\ref{sec:proj_drake_swainson} & FN \\
+textbook (Eides one-loop; tabulated Bethe log) & TA \\
@@ -4842 +4813 @@ VP 2$P_{1/2}$ & $0$ &
-\S\ref{sec:proj_fock} (parity zero) & FN \\
+textbook (zero at this order) & TA \\
@@ -4863,2 +4834,6 @@ $63.599$, which implies an H reference of $8.7913/63.599 = 0.13823$~MHz.
-Each component is
-tagged to its \S\ref{sec:layer2} projection chain.  All entries are
+The four leading one-loop rows are the textbook assembly of Paper~36:
+the one-loop expressions of Eides \emph{et al.}\ with the tabulated
+Drake--Swainson Bethe logarithms, and no GeoVac projection enters them
+($S^3$ never enters the number; the Coulomb-Sturmian basis of Paper~36
+is a numerical device for computing Bethe logarithms, not used here).
+All entries are
@@ -4867,2 +4842 @@ self-energy is a $-12.88$~MHz level shift, entering the splitting as
-$+12.88$~MHz.  Status: FN
-$=$ framework-native; L2 $=$ Layer-2 input.  Source memo:
+$+12.88$~MHz.  Status: TA $=$ textbook assembly; L2 $=$ Layer-2 input.  Source memo:
@@ -4873 +4847 @@ $=$ framework-native; L2 $=$ Layer-2 input.  Source memo:
-\textbf{Framework-native subtotal:} $+1052.19$ MHz, the leading one-loop
+\textbf{One-loop subtotal:} $+1052.19$ MHz, the leading one-loop
@@ -4876 +4850 @@ assembly ($-0.534\%$).  An earlier version of this table also listed an
-averaging'' row of $\sim +5.00$~MHz, giving a framework-native subtotal of
+averaging'' row of $\sim +5.00$~MHz, giving a subtotal of
@@ -4885,4 +4859,5 @@ attributed (Paper~36).
-\textbf{Structural reading.}  The autopsy tags each listed component to its
-projection chain; FNS is stated in projection-chain language rather than as
-a black-box external constant.  A complete order-by-order budget is owed
-before the components can be compared with the measurement.
+\textbf{Structural reading.}  The leading one-loop components carry no
+GeoVac projection chain; only the FNS and recoil rows are stated in
+projection-chain language, as Layer-2 inputs.  A complete order-by-order
+budget is owed before the components can be compared with the
+measurement.
@@ -5261 +5236 @@ one-loop assembly (\S\ref{sec:autopsy_lamb}): identical lepton
-(electron), identical projection chains, only the nuclear charge $Z=2$
+(electron), identical assembly, only the nuclear charge $Z=2$
@@ -5274,4 +5249,4 @@ to $\sim$ppm and the full Uehling kernel matches contact form at 0.5\%.
-$E(2S_{1/2}) - E(2P_{1/2})$).  The framework-native subtotal
-$+13\,841.11$ MHz comes from the same projection-chain machinery that
-closed the H Lamb shift (\S\ref{sec:autopsy_lamb}) with $Z = 2$
-substituted throughout: Bohr-Fock Dirac (zero by 2$S_{1/2}=2P_{1/2}$
+$E(2S_{1/2}) - E(2P_{1/2})$).  The subtotal $+13\,841.11$ MHz
+(components 1--4) is the textbook one-loop assembly of the H Lamb shift
+(\S\ref{sec:autopsy_lamb}; Paper~36), in which no GeoVac projection
+enters the one-loop terms, with $Z = 2$ substituted throughout: Bohr-Fock Dirac (zero by 2$S_{1/2}=2P_{1/2}$
@@ -5295,3 +5270,2 @@ Component & Value (MHz) & Projection chain & Status \\
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_sturmian} $\circ$
-\S\ref{sec:proj_spectral_action} $\circ$ \S\ref{sec:proj_drake_swainson}
-& FN ($\alpha^3 Z^4$ truncation) \\
+textbook (Eides one-loop; tabulated Bethe logs)
+& TA ($\alpha^3 Z^4$ truncation) \\
@@ -5299,2 +5273,2 @@ Component & Value (MHz) & Projection chain & Status \\
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_spectral_action}
-& FN ($\beta = 137$, contact form $0.5\%$ check) \\
+textbook (Uehling kernel)
+& TA ($\beta = 137$, contact form $0.5\%$ check) \\
@@ -5315,2 +5289,2 @@ Layer-2 input & W3 inner-factor (smallest in catalogue) \\
-\textbf{Framework-native subtotal} & \textbf{+13\,841.11} & (1+2+3+4) & \\
-Framework $+$ literature total & $+13\,859.58$ & (1+2+3+4+5a+5b+5c+6+7) & \\
+\textbf{One-loop assembly subtotal} & \textbf{+13\,841.11} & (1+2+3+4) & \\
+Assembly $+$ literature total & $+13\,859.58$ & (1+2+3+4+5a+5b+5c+6+7) & \\
@@ -5326,2 +5300,2 @@ Framework $+$ literature total & $+13\,859.58$ & (1+2+3+4+5a+5b+5c+6+7) & \\
-Three component ratios He$^+$/H test whether the projection-chain
-machinery transfers cleanly Z=1 $\to$ Z=2:
+Three component ratios He$^+$/H test whether the implementation
+transfers cleanly Z=1 $\to$ Z=2:
@@ -5351,5 +5325,5 @@ architectural modification.
-All three structural tests pass.  The projection chains are
-\textbf{$Z$-portable}: substituting $Z=2$ into the existing pipeline
-gives He$^+$ with no architectural changes.  This is what Paper~34's
-projection taxonomy was designed to guarantee --- projections describe
-physics, not specific $Z$ values.
+All three tests pass: substituting $Z=2$ into the existing pipeline
+gives He$^+$ with no architectural changes.  For the SE and VP rows this
+checks the implementation of textbook formulas, whose $Z$-dependence is
+known in closed form; only the FNS row tests a Paper~34 projection
+(\S\ref{sec:proj_charge_density}).
@@ -5357 +5331 @@ physics, not specific $Z$ values.
-\textbf{Residual attribution.}  The framework-native subtotal is
+\textbf{Residual attribution.}  The assembly subtotal is
@@ -5397,2 +5371,2 @@ total, different per-component itemization.  See
-First, the framework's projection-chain machinery is \emph{Z-portable}
-at the structural level: all three Z-scaling tests pass cleanly, with
+First, the implementation is \emph{Z-portable}: all three Z-scaling
+tests pass cleanly, with
@@ -5402 +5376 @@ kernel-mixed, see above).
-Second, the $-1.42\%$ native residual is the Z-amplified analog of
+Second, the $-1.42\%$ assembly residual is the Z-amplified analog of
@@ -7871,2 +7845,4 @@ on it are withdrawn until a verifiable reference is supplied.
-The framework total, $1047.63$~MHz, stands against the Mariam~1982
-experiment $1054(22)$~MHz (\S\ref{sec:offprecision}).
+The assembly total, $1047.63$~MHz (recomputed by
+\texttt{tests/paper36\_precision\_support/precision\_catalogue\_muonium\_lamb.py}),
+lies $+0.04\%$ from the Mu-MASS measurement $1047.2(2.3)(1.1)$~MHz
+(\S\ref{sec:offprecision}).
@@ -9095 +9071 @@ rationale; for full discussion see \texttt{debug/lorentzian\_l0\_audit\_memo.md}
-13 & Drake--Swainson asymptotic subtraction & T.F. & H & Central tier (Paper~31 TS-C); signature-blind. \\
+13 & Drake--Swainson normalization of Bethe-log sums & T.F. & H & Central tier (Paper~31 TS-C); signature-blind. \\
@@ -9108 +9084 @@ rationale; for full discussion see \texttt{debug/lorentzian\_l0\_audit\_memo.md}
-26 & Gauge choice (Coulomb/Lorenz/Feynman) & W.M.F. & MH & Coulomb-gauge $1/(4\pi) = $ M1; Ward identity holds. \\
+26 & Gauge choice (Coulomb/Lorenz/Feynman) & W.M.F. & MH & Coulomb-gauge $1/(4\pi)$ is the vertex model's normalization (Paper~33); Ward identity is standard QED, not established on the graph. \\
@@ -9273,3 +9249,5 @@ The catalogue rows of \S~\ref{sec:catalog} satisfy this consistency
-check across all eight fully-tagged entries with non-trivial chains
-(Lamb shift, the Paper~2 $\alpha$ Observation row, Stefan-Boltzmann, He CI,
-hydrogen $E_n$, scalar/Dirac spatial Casimir, Uehling VP).
+check across all seven fully-tagged entries with non-trivial chains
+(the Paper~2 $\alpha$ Observation row, Stefan-Boltzmann, He CI,
+hydrogen $E_n$, scalar/Dirac spatial Casimir, Uehling VP).  The
+hydrogen Lamb shift is not among them: it is a textbook assembly with
+no GeoVac projection chain (Paper~36).
@@ -9283 +9261 @@ expected target signature matches the chain signature; this passes
-for all eight tested rows.
+for all seven tested rows.
@@ -9358,3 +9336,3 @@ retardation.  The structurally distinguishing fact is that they
-operate on disjoint Layer-2 content (Drake--Swainson regularizes
-divergent Bethe-log spectral sums; Breit retardation adds two-body
-kinematic content to the framework's single-particle Lamb-shift
+operate on disjoint Layer-2 content (Drake--Swainson normalizes
+Bethe-log spectral sums; Breit retardation adds two-body
+kinematic content to the single-particle Lamb-shift
@@ -9439,3 +9417,6 @@ $(\text{D-only})$ bucket.  The transcendental axis admits an
-independence witness: vector-photon promotion (\S~\ref{sec:layer2},
-projection~11) introduces $1/(4\pi)$ without adding any variable
-or dimension.  Two projections (\texttt{spectral\_action},
+independence witness: the Hopf bundle (\S~\ref{sec:proj_hopf},
+projection~2) introduces $\pi = \mathrm{Vol}(S^2)/4$ without adding any
+variable or dimension, and Wick rotation does the same with
+$2\pi = \mathrm{Vol}(S^1)$.  Vector-photon promotion is not a witness:
+it carries only the vertex model's $1/(4\pi)$ normalization
+(Paper~33).  Two projections (\texttt{spectral\_action},
@@ -9488,2 +9469,3 @@ ceilings ($\approx 0.024\%$ for He NR variational; the H Lamb shift's
-$-0.534\%$ is an order truncation with tabulated Bethe logarithms, not a
-basis ceiling, Paper~36).  The residual is bounded
+$-0.534\%$ is not in this class: it is a textbook one-loop assembly with
+tabulated Bethe logarithms and no GeoVac projection chain, Paper~36).
+The residual is bounded
@@ -9597,4 +9579,3 @@ of the retired baseline), and the He oscillator strength ($+3.4\%$).  At depth~4
-the catalogue contains $\mu$H Bohr--Fermi ($+2$~ppm) and the H
-Lamb shift at leading one-loop order ($-0.534\%$).  Residuals at every
-depth $\geq 2$ are scattered across more than four orders of
-magnitude.  Depth correlates with \emph{opportunity} for
+the catalogue contains $\mu$H Bohr--Fermi ($+2$~ppm).  Residuals at
+depth~3 alone are scattered across more than four orders of
+magnitude, from exact to percent level.  Depth correlates with \emph{opportunity} for
@@ -9703,2 +9684,2 @@ matches do not.
-are now documented.  The thirteenth (Drake--Swainson asymptotic
-subtraction) was added in sprint LS-4 for the $\ell > 0$ Bethe
+are now documented.  The thirteenth (Drake--Swainson normalization of
+Bethe-log sums) was added in sprint LS-4 for the $\ell > 0$ Bethe
@@ -9749,2 +9730,4 @@ the V/D/transcendental axis-grid completeness identified by the
-TX-A audit: gauge choice provides a second independence witness next
-to vector-photon promotion in the transcendental-only corner; Wick
+TX-A audit: gauge choice fills the no-variable, no-dimension row
+beside vector-photon promotion (neither carries a transcendental of
+physical content: the $1/(4\pi)$ they share is the vertex model's
+normalization, Paper~33); Wick
@@ -9842,2 +9825 @@ $\epsilon_\text{basis}(\text{framework})$ for each major chain
-class (one-loop QED Lamb-shift class, NR variational class,
-graph-native CI class) so that the L2-count $= 0$ side of the
+class (NR variational class, graph-native CI class) so that the L2-count $= 0$ side of the
@@ -9854,10 +9836,6 @@ basis-quality limit.
-The LS-4 Lamb shift datapoint adds historical context.  The LS-4
-chained projection (Fock + spectral action + Sturmian +
-Drake--Swainson) achieves $-0.39\%$ residual at the sweet-spot
-$N=40$.  At converged $N$, the LS-4 result reverts to the LS-1
-ceiling $-3.10\%$ (LS-4, pre-LS-6a convention); the LS-6a fix supersedes this with $-0.534\%$.  Under
-Prediction~\ref{pred:l2_presence} the $-0.534\%$ is \emph{not} a
-basis-quality value: the assembly uses tabulated Bethe logarithms and
-stops before the one-loop $\alpha(Z\alpha)^5$ term (Paper~36); the
-$-0.39\%$ at $N=40$ is a transient finite-basis cancellation.  Prediction~\ref{pred:l2_presence} applies to
-converged-basis residuals.
+The hydrogen Lamb shift is not a datapoint for
+Prediction~\ref{pred:l2_presence}: its $-0.534\%$ is a textbook
+one-loop assembly with tabulated Bethe logarithms that stops before the
+one-loop $\alpha(Z\alpha)^5$ term, and no GeoVac projection chain
+enters it (Paper~36).  Prediction~\ref{pred:l2_presence} applies to
+converged-basis residuals of projection chains.
@@ -9922,2 +9900,3 @@ is open.  Two candidate additions surfaced by the TX-A axiomatization
-sprint: (a)~a gauge-choice projection (Coulomb vs Lorenz) that would
-sit in the pure-transcendental bucket alongside vector-photon, (b)~a
+sprint: (a)~a gauge-choice projection (Coulomb vs Lorenz), since added
+as \S~\ref{sec:proj_gauge_choice}, which preserves V and D but has no
+transcendental measured in the framework, (b)~a
@@ -10375 +10354 @@ admits six scalar instances plus one rank-2 tensor instance, while
-$[E]$ admits four.  Sprint~1 (2026-05-15) added six more
+$[E]$ admits three.  Sprint~1 (2026-05-15) added six more
@@ -10386,3 +10365,6 @@ change; apparatus identity / state-side reduction), closing the
-axis-grid completeness gap: gauge choice provides a second
-independence witness in the transcendental-only corner of the
-V/D/transcendental grid, Wick rotation promotes the
+axis-grid completeness gap: gauge choice preserves V and D with no
+transcendental measured in the framework (the transcendental-only corner
+is witnessed by volume-slot factors, the base-sphere $\mathrm{Vol}(S^2)/4$
+of the Hopf bundle and the Wick-rotation time circle, not by the gauge
+sector), Wick rotation
+promotes the
diff --git a/papers/group6_precision_observations/paper_35_time_as_projection.tex b/papers/group6_precision_observations/paper_35_time_as_projection.tex
index 4e589392..5f5b164a 100644
--- a/papers/group6_precision_observations/paper_35_time_as_projection.tex
+++ b/papers/group6_precision_observations/paper_35_time_as_projection.tex
@@ -1350,5 +1350,3 @@ to close the multi-loop ceiling.  \emph{Second}, the original test
-target was misidentified.  The $\alpha^5$-order contributions of all
-kinds to the 2$S$ Lamb shift total approximately $+7.10$~MHz
-(Eides--Grotch--Shelyuto 2001 Table 7.4 cumulative central value), most
-of it the one-loop $\alpha(Z\alpha)^5$ term (about $+6.9$~MHz), not
-the $\sim +29$~MHz LS-3 residual.  Sprint LS-6a (May 2026)
+target was misidentified.  The complete one-loop $\alpha(Z\alpha)^5$
+correction to the 2$S$ level is $+7.13$~MHz (Eides--Grotch--Shelyuto
+2001, Table~III; Paper~36), not the $\sim +29$~MHz LS-3 residual.  Sprint LS-6a (May 2026)
diff --git a/papers/synthesis/group6_precision_observations_synthesis.tex b/papers/synthesis/group6_precision_observations_synthesis.tex
index fc8ae9de..f0ab0ad7 100644
--- a/papers/synthesis/group6_precision_observations_synthesis.tex
+++ b/papers/synthesis/group6_precision_observations_synthesis.tex
@@ -173 +173 @@ Wigner $3j$ & dimensionless & $\mathbb{Q}[\sqrt{2k{+}1}]_k$ \\
-Vector-photon & dimensionless & $1/(4\pi)$ per loop \\
+Vector-photon & dimensionless & $1/(4\pi)$, a vertex-model normalization \\
@@ -197,4 +197,4 @@ anonymous transcendental:\ each is pinned to a projection row, which is
-the discipline this group exists to model.  The per-loop
-$S^2$ Weyl constant $1/(4\pi)$ (M1 volume slot) is tagged to the vector-photon
-promotion; the Paper~2~\cite{loutey_paper2} \gvc{k-combination-rule}{combination rule
-$K = \pi(B+F-\Delta)$ is a separate object (it does not consume the per-loop $1/(4\pi)$) and
+the discipline this group exists to model.  The $1/(4\pi)$ tagged to
+vector-photon promotion is the vertex model's unit-$S^2$ harmonic
+normalization, which the physical vertex does not carry (Paper~33); the Paper~2~\cite{loutey_paper2} \gvc{k-combination-rule}{combination rule
+$K = \pi(B+F-\Delta)$ is a separate object (it does not consume that $1/(4\pi)$) and
@@ -425,5 +425,5 @@ chain and a framework-native (FN) versus Layer-2-input (L2) status.  The
-hydrogen $2S_{1/2}$--$2P_{1/2}$ Lamb-shift autopsy is representative:\ the framework-native
-components (self-energy via Fock\,$\circ$\,Sturmian\,$\circ$\,spectral
-action\,$\circ$\,Drake--Swainson, Uehling vacuum polarization via
-Fock\,$\circ$\,spectral action, $\ldots$) sum to $+1052.19$~MHz, the
-leading one-loop assembly with tabulated Bethe logarithms; an earlier
+hydrogen $2S_{1/2}$--$2P_{1/2}$ Lamb-shift autopsy is the limiting
+case:\ its leading components (one-loop self-energy and Uehling
+vacuum polarization) are a textbook assembly with tabulated Bethe
+logarithms, in which no GeoVac projection and no $S^3$ enters (Paper~36),
+and they sum to $+1052.19$~MHz; an earlier
@@ -433,2 +433,2 @@ term (about $+6.9$~MHz) is not itemized, so the residual is not
-attributed (Papers~34, 36).  The autopsy makes the projection-chain
-dictionary's reach visible at the observable level.
+attributed (Papers~34, 36).  The autopsies make the projection-chain
+dictionary's reach, and its limits, visible at the observable level.
@@ -464,2 +464,2 @@ synthesis pass falsified that form on the catalogue itself (depth-3
-chains exist at residual~$0$ and at $-211$~ppm; depth-4 chains at
-$+2$~ppm and at $-0.534\%$).  The replacement is the \emph{Layer-2-presence
+chains exist at residual~$0$ and at $-211$~ppm; a depth-4 chain at
+$+2$~ppm).  The replacement is the \emph{Layer-2-presence
@@ -475,2 +475,3 @@ hydrogen static polarizability $\tfrac92 a_0^3$, $S^3$ Casimir $1/240$) to
-basis-quality ceilings (the one-loop Lamb shift's $-0.534\%$ is instead
-an order truncation with tabulated Bethe logarithms);
+basis-quality ceilings (the one-loop Lamb shift's $-0.534\%$ is not
+such a chain: it is a textbook assembly with tabulated Bethe
+logarithms, Paper~36);
```

## 2026-11-08: muonic Lamb rows (textbook assembly)

Owner pass on Paper 34 applying the hydrogen Lamb ruling above (item 3) to the three muonic
Lamb-shift entries (muonic hydrogen, muonic deuterium, muonic helium-4 ion): the one-loop Uehling
term (flat-space vacuum-polarization function entered as an input) and the Eides one-loop
self-energy bracket are a textbook assembly (TA), not framework-native (FN); closed-form
Friar/Foldy finite-size rows keep their charge-density chain with the radius as a Layer-2 input
(L2), as the hydrogen FNS row does.  No number changed.  Loci: tab:catalog (muonic Uehling, mu-d
Lamb, mu-H Lamb 202.17 rows), tab:offprecision (mu-H Lamb 200.50 and mu-4He+ Lamb rows), the three
autopsy subsections (captions, Status/chain columns, headline and structural-reading paragraphs),
and the L2-count >= 2 item of the Layer-2-presence section (also corrected "deuteron
polarizability" to "proton polarizability" for mu-H, per its own autopsy table).  Every removed
line is preserved verbatim below as the `-` lines of a zero-context unified diff against the
working-tree state immediately before this pass (the `+` lines are the replacements).

```diff
--- a/papers/group6_precision_observations/paper_34_projection_taxonomy.tex
+++ b/papers/group6_precision_observations/paper_34_projection_taxonomy.tex
@@ -3736 +3736 @@ $\mu$H 2$S$--2$P$ Uehling $+205.0074$ meV full kernel
-(Sprint MH Track A) & spectral action $\circ$ Fock $\circ$ rest-mass
+(Sprint MH Track A) & entered $\Pi$ (flat-space) $\circ$ Fock $\circ$ rest-mass
@@ -3739 +3739 @@ $\pi^2{\cdot}\mathbb{Q}$ & $<1$ ppm vs Antognini~2013 \\
-$\mu$d 2$S$--2$P$ Lamb $+202.63$ meV framework $+$ Layer-2
+$\mu$d 2$S$--2$P$ Lamb $+202.63$ meV textbook assembly $+$ Layer-2
@@ -3741,2 +3741,2 @@ $\mu$d 2$S$--2$P$ Lamb $+202.63$ meV framework $+$ Layer-2
-Roothaan autopsy in \S\ref{sec:autopsy_mud_lamb}) & spectral action
-$\circ$ Fock $\circ$ rest-mass projection ($m_e \to m_\mu$, nucleus
+Roothaan autopsy in \S\ref{sec:autopsy_mud_lamb}) & entered $\Pi$
+(flat-space) $\circ$ Fock $\circ$ rest-mass projection ($m_e \to m_\mu$, nucleus
@@ -3756 +3756 @@ $\mu^4$He$^+$ 2$S$--2$P$ Uehling $+1665.77$ meV full kernel
-autopsy in \S\ref{sec:autopsy_muhe4_lamb}) & spectral action $\circ$
+autopsy in \S\ref{sec:autopsy_muhe4_lamb}) & entered $\Pi$ (flat-space) $\circ$
@@ -3845,3 +3845,3 @@ $\circ$ Fock & $Z$ & dimensionless & flow & $+0.07\%$ at $N=40$ \\
-$\mu$H 2$S$--2$P$ Lamb shift $202.17$ meV (framework $+$ literature,
-Sprint MH Track A) & Fock $\circ$ spectral action $\circ$ Sturmian
-$\circ$ Drake--Swainson $\circ$ rest-mass projection ($m_e \to m_\mu$)
+$\mu$H 2$S$--2$P$ Lamb shift $202.17$ meV (textbook assembly $+$ literature,
+Sprint MH Track A) & textbook assembly (entered flat-space $\Pi$, Eides
+self-energy bracket, closed-form Friar) $\circ$ rest-mass projection ($m_e \to m_\mu$)
@@ -4245,3 +4245,3 @@ Sprint precision-catalogue, 2026-05-08.  \\
-$\mu$H 2$S$--2$P$ Lamb framework-native $200.50$ meV (no literature
-inputs, Sprint MH Track A) & Fock $\circ$ spectral action $\circ$
-Sturmian $\circ$ Drake--Swainson $\circ$ rest-mass projection
+$\mu$H 2$S$--2$P$ Lamb one-loop assembly $200.50$ meV (no literature
+inputs, Sprint MH Track A) & textbook assembly (entered flat-space
+$\Pi$, Eides self-energy bracket, closed-form Friar) $\circ$ rest-mass projection
@@ -4250,2 +4250,2 @@ Sturmian $\circ$ Drake--Swainson $\circ$ rest-mass projection
-& Framework reproduces dominant Uehling at $<1$ ppm vs Antognini and
-Friar moment leading-order at $4.3\%$.  $-0.92\%$ residual decomposes to within $0.20$~meV into missing K\"all\'en--Sabry two-loop VP ($+1.51$ meV),
+& The full Uehling kernel integration reproduces Antognini at $<1$ ppm and the
+closed-form Friar moment leading order at $4.3\%$.  $-0.92\%$ residual decomposes to within $0.20$~meV into missing K\"all\'en--Sabry two-loop VP ($+1.51$ meV),
@@ -4715 +4715 @@ HFS (closes Li atomic accuracy floor $1.10\%$ $\to$ sub-$0.1\%$
-$\mu^4$He$^+$ 2$S$--2$P$ Lamb $+1358.66$ meV framework $+$ Layer-2
+$\mu^4$He$^+$ 2$S$--2$P$ Lamb $+1358.66$ meV textbook assembly $+$ Layer-2
@@ -4718,2 +4718,2 @@ $\mu^4$He$^+$ 2$S$--2$P$ Lamb $+1358.66$ meV framework $+$ Layer-2
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_spectral_action}
-$\circ$ \S\ref{sec:proj_sturmian} $\circ$ \S\ref{sec:proj_restmass}
+textbook assembly (entered flat-space $\Pi$, Eides self-energy
+bracket) $\circ$ \S\ref{sec:proj_restmass}
@@ -4731 +4731 @@ mass-hierarchy $\times$ Z combination tested: $m_\text{red} \sim
-at $\mu$d).  Framework-native VP at $-0.032\%$ vs Diepold 2018
+at $\mu$d).  Full-kernel Uehling VP at $-0.032\%$ vs Diepold 2018
@@ -4969,3 +4969,2 @@ Pohl~et~al.; muonic convention $E(2P_{1/2}) - E(2S_{1/2})$).
-The CREMA 2010 muonic hydrogen Lamb shift decomposes under five
-named projection-chain components plus three companion Layer-2
-inputs.
+The CREMA 2010 muonic hydrogen Lamb shift decomposes into five
+named components plus three companion Layer-2 inputs.
@@ -4977,9 +4976,9 @@ Lamb shift (muonic convention $E(2P_{1/2}) - E(2S_{1/2})$).  The
-framework-native subtotal $+200.50$ meV comes from full Uehling
-integration via \S\ref{sec:proj_spectral_action}, the self-energy
-bracket under rest-mass projection
-(\S\ref{sec:proj_restmass}~+~\S\ref{sec:proj_sturmian}), and the
-operator-level Friar moment via the magnetization-density module
-(\S\ref{sec:proj_magnetization_density}; cf.\
-\S\ref{sec:proj_charge_density} for the parent charge-density
-projection).  The $+1.67$ meV gap to the $+202.17$ meV
-framework-plus-literature total is closed by Antognini~2013
+subtotal $+200.50$ meV is a textbook assembly in which no GeoVac
+projection enters the one-loop terms: the full Uehling integration with
+the flat-space vacuum-polarization function $\Pi$ entered as an input,
+the Eides one-loop self-energy bracket at the muonic reduced mass, and
+the closed-form Friar moment (Eides Eq.~2.35) with $r_p$ a Layer-2
+input (the operator-level Friar route via
+\S\ref{sec:proj_magnetization_density} is compared below).
+The $+1.67$ meV gap to the $+202.17$ meV
+assembly-plus-literature total is closed by Antognini~2013
@@ -4988 +4987,2 @@ mixed VP-VP/VP-SE, $\alpha^7$) and QCD-internal proton
-polarizability (W3 inner-factor).  Source memo:
+polarizability (W3 inner-factor).  Status: TA $=$ textbook assembly;
+L2 $=$ Layer-2 input.  Source memo:
@@ -4996,2 +4996,2 @@ Component & Value (meV) & Projection chain & Status \\
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_spectral_action}
-& FN ($-0.10$ ppm vs Antognini) \\
+textbook (full Uehling kernel; entered flat-space $\Pi$)
+& TA ($-0.10$ ppm vs Antognini) \\
@@ -4999,3 +4999,2 @@ Component & Value (meV) & Projection chain & Status \\
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_sturmian}
-$\circ$ \S\ref{sec:proj_restmass}
-& FN ($+24\%$ vs Antognini, recoil-mixing gap) \\
+textbook (Eides one-loop bracket at $m_{\text{red},\mu}$)
+& TA ($+24\%$ vs Antognini, recoil-mixing gap) \\
@@ -5004,2 +5003,2 @@ $\circ$ \S\ref{sec:proj_restmass}
-$/$ \S\ref{sec:proj_magnetization_density}
-& FN ($-4.3\%$ vs Antognini) \\
+$+$ Layer-2 $r_p$
+& L2 ($-4.3\%$ vs Antognini) \\
@@ -5019 +5018 @@ vs CREMA $+202.3706(23)$ meV} \\
-\textbf{Headline framework computation.}  The full Uehling integration
+\textbf{Headline computation.}  The full Uehling integration
@@ -5046 +5045 @@ as ``VP-VP mixed'' ($+0.150$~meV), with one-loop reported at
-$+205.0074$~meV.  The framework's $-0.10$~ppm match against Antognini
+$+205.0074$~meV.  The full-kernel $-0.10$~ppm match against Antognini
@@ -5074,12 +5073,11 @@ calibration for the HFS Zemach channel.  See
-\textbf{Structural reading.}  The autopsy realizes the
-structural-skeleton-scope statement at the operator level for a
-cross-register multi-focal observable in the muonic regime.  Three
-framework-native components reproduce $+200.50$ meV (the
-five-projection-chain skeleton); three Layer-2 inputs close the
-remaining $+1.67$ meV gap, attributable cleanly to LS-8a multi-loop
-renormalization (Component 4) and W3 inner-factor (Component 5).
-The same architecture that closes the H 1S Lamb shift at
-\S\ref{sec:autopsy_lamb} closes $\mu$H 2$S$--2$P$ here, with
-Uehling-kernel kinematics (contact-form $\to$ full kernel) and
-Friar-profile calibration (first moment $\to$ RMS) being the only
-two ingredients that change between regimes.
+\textbf{Structural reading.}  Components 1--3 are a textbook
+assembly, $+200.50$ meV: the flat-space $\Pi$ is entered rather than
+produced by the spectral-action projection, the self-energy is the
+Eides one-loop bracket, and the Friar term is closed-form with $r_p$ a
+Layer-2 input; $S^3$ never enters these numbers.  Three Layer-2 inputs
+close the remaining $+1.67$ meV gap, attributable cleanly to LS-8a
+multi-loop renormalization (Component 4) and W3 inner-factor
+(Component 5).  Relative to the electronic one-loop assembly of
+\S\ref{sec:autopsy_lamb}, Uehling-kernel kinematics (contact-form
+$\to$ full kernel) and Friar-profile calibration (first moment $\to$
+RMS) are the only two ingredients that change between regimes.
@@ -5109,7 +5107,7 @@ deuterium Lamb shift (muonic convention $E(2P_{1/2})-E(2S_{1/2})$).
-The framework-native subtotal $+198.91$ meV comes from full Uehling
-integration via \S\ref{sec:proj_spectral_action}, the self-energy
-bracket under rest-mass projection
-(\S\ref{sec:proj_restmass}~+~\S\ref{sec:proj_sturmian}), and the
-closed-form Friar moment at the deuteron's measured RMS charge
-radius via \S\ref{sec:proj_charge_density}.  The $+3.71$ meV gap
-to the $+202.63$ meV framework-plus-literature total is closed by
+The subtotal $+198.91$ meV is a textbook assembly in which no GeoVac
+projection enters the one-loop terms: the full Uehling integration with
+the flat-space $\Pi$ entered as an input, the Eides one-loop
+self-energy bracket at the muonic reduced mass, and the closed-form
+Friar moment at the deuteron's measured RMS charge radius (a Layer-2
+input; \S\ref{sec:proj_charge_density}).  The $+3.71$ meV gap
+to the $+202.63$ meV assembly-plus-literature total is closed by
@@ -5120 +5118,2 @@ contribution---deuteron polarizability (Carlson--Gorchtein--
-Vanderhaeghen 2014).  Source memo:
+Vanderhaeghen 2014).  Status: TA $=$ textbook assembly; L2 $=$
+Layer-2 input.  Source memo:
@@ -5128,2 +5127,2 @@ Component & Value (meV) & Projection chain & Status \\
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_spectral_action}
-& FN ($-0.50\%$ vs Krauth) \\
+textbook (full Uehling kernel; entered flat-space $\Pi$)
+& TA ($-0.50\%$ vs Krauth) \\
@@ -5131,3 +5130,2 @@ Component & Value (meV) & Projection chain & Status \\
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_sturmian}
-$\circ$ \S\ref{sec:proj_restmass}
-& FN ($+41\%$ vs Krauth, same recoil-mixing gap as $\mu$H widened
+textbook (Eides one-loop bracket at $m_{\text{red},\mu d}$)
+& TA ($+41\%$ vs Krauth, same recoil-mixing gap as $\mu$H widened
@@ -5137 +5135,2 @@ at heavier nucleus) \\
-& FN ($+0.024\%$ vs Krauth at $r_d$, bit-tight) \\
+$+$ Layer-2 $r_d$
+& L2 ($+0.024\%$ vs Krauth at $r_d$, bit-tight) \\
@@ -5153 +5152 @@ $-0.12\%$ vs CREMA $+202.8785(34)$ meV} \\
-\textbf{Headline framework computation.}  The Friar moment leading
+\textbf{Headline computation.}  The Friar moment leading
@@ -5156 +5155 @@ tabulated leading line at the deuteron's measured radius at
-$+0.024\%$ --- the cleanest framework-vs-literature match for the
+$+0.024\%$ --- the cleanest match for the
@@ -5165 +5164 @@ The scaling
-is reproduced by the framework's bare Eides calculation to 4 digits.
+is reproduced by the bare Eides calculation to 4 digits.
@@ -5208,4 +5207,5 @@ Friar--Payne~2005 and Pachucki--Yerokhin~2010+).
-\textbf{Structural reading.}  The autopsy validates the multi-focal
-architecture under nuclear-isotope substitution at the sub-percent
-precision level on a precision-frontier observable.  Three
-framework-native components reproduce $+198.91$ meV; four Layer-2
+\textbf{Structural reading.}  Components 1--3 are the same textbook
+assembly as at $\mu$H (entered flat-space $\Pi$, Eides one-loop
+bracket, closed-form Friar with a Layer-2 radius), $+198.91$ meV; no
+GeoVac projection chain produces them and $S^3$ never enters these
+numbers.  Four Layer-2
@@ -5216,2 +5216,2 @@ multi-loop QED that dominated at $\mu$H.  The rest-mass projection
-(\S\ref{sec:proj_restmass}) transports the entire skeleton (Uehling
-integrator, SE bracket, Friar leading) verbatim; only the nuclear
+(\S\ref{sec:proj_restmass}) carries the assembly (Uehling
+integrator, SE bracket, Friar leading) over unchanged; only the nuclear
@@ -5219 +5219 @@ charge-radius input ($r_p \to r_d$) and the polarizability
-Layer-2 input change between isotopes.  Same architecture, same
+Layer-2 input change between isotopes.  Same assembly, same
@@ -6983,6 +6983,6 @@ ion Lamb shift (muonic convention $E(2P_{1/2})-E(2S_{1/2})$).  The
-framework-native subtotal $+1358.24$ meV comes from full Uehling
-integration via \S\ref{sec:proj_spectral_action}, the self-energy
-bracket under rest-mass projection
-(\S\ref{sec:proj_restmass}~+~\S\ref{sec:proj_sturmian}), and the
-closed-form Foldy charge-density at the $\alpha$-particle's measured
-RMS charge radius via \S\ref{sec:proj_charge_density}.  No
+subtotal $+1358.24$ meV is a textbook assembly in which no GeoVac
+projection enters the one-loop terms: the full Uehling integration with
+the flat-space $\Pi$ entered as an input, the Eides one-loop
+self-energy bracket at the muonic reduced mass, and the
+closed-form Foldy charge-density term at the $\alpha$-particle's measured
+RMS charge radius (a Layer-2 input; \S\ref{sec:proj_charge_density}).  No
@@ -6994 +6994,2 @@ CREMA 2021 measurement modulo two convention surfaces (SE bracket
-recoil-mixing, Friar HO absorption) documented below.  Source memo:
+recoil-mixing, Friar HO absorption) documented below.  Status: TA $=$
+textbook assembly; L2 $=$ Layer-2 input.  Source memo:
@@ -7002,2 +7003,2 @@ Component & Value (meV) & Projection chain & Status \\
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_spectral_action}
-& FN ($-0.032\%$ vs Diepold 2018) \\
+textbook (full Uehling kernel; entered flat-space $\Pi$)
+& TA ($-0.032\%$ vs Diepold 2018) \\
@@ -7005,3 +7006,2 @@ Component & Value (meV) & Projection chain & Status \\
-\S\ref{sec:proj_fock} $\circ$ \S\ref{sec:proj_sturmian}
-$\circ$ \S\ref{sec:proj_restmass}
-& FN ($+8.4\%$ vs Diepold, recoil-mixing gap shrunken at higher Z) \\
+textbook (Eides one-loop bracket at $m_{\text{red},\mu^4\text{He}}$)
+& TA ($+8.4\%$ vs Diepold, recoil-mixing gap shrunken at higher Z) \\
@@ -7010,2 +7010,2 @@ $\circ$ \S\ref{sec:proj_restmass}
-$\circ$ \S\ref{sec:proj_restmass}
-& FN ($+2.9\%$ vs Diepold leading, $+4.0\%$ vs Diepold full FNS) \\
+$\circ$ \S\ref{sec:proj_restmass} $+$ Layer-2 $r_\alpha$
+& L2 ($+2.9\%$ vs Diepold leading, $+4.0\%$ vs Diepold full FNS) \\
@@ -7027 +7027 @@ $-1.44\%$ vs CREMA $+1378.521(48)$ meV (Krauth 2021)} \\
-\textbf{Headline framework computation.}  The full Uehling integration
+\textbf{Headline computation.}  The full Uehling integration
@@ -7032,3 +7032,2 @@ the catalogue at heavy-Z, transporting bit-identical kernel
-infrastructure from the $\mu$H regime to the $\mu^4$He$^+$ regime.
-This is the architectural transport theorem of the muonic Lamb
-catalogue: the same Uehling-kernel integral closes verbatim under
+infrastructure from the $\mu$H regime to the $\mu^4$He$^+$ regime:
+the same Uehling-kernel integral closes verbatim under
@@ -7129,5 +7128,6 @@ leading directly.
-\textbf{Structural reading.}  The autopsy validates the multi-focal
-architecture under the most extreme mass-hierarchy / Z combination in
-the precision catalogue.  Three framework-native components reproduce
-$+1358.24$ meV; four Layer-2 inputs (one less than $\mu$H since no
-Zemach is needed) close the remaining gap, with component-level
+\textbf{Structural reading.}  Components 1--3 are the same textbook
+assembly as at $\mu$H and $\mu$d (entered flat-space $\Pi$, Eides
+one-loop bracket, closed-form Foldy term with a Layer-2 radius),
+$+1358.24$ meV; no GeoVac projection chain produces them and $S^3$
+never enters these numbers.  Four Layer-2 inputs close the remaining
+gap, with component-level
@@ -7136,2 +7136,2 @@ mixing $+0.92$ meV; r$_\alpha$ extraction signal $\sim 7$ meV) summing
-to the $-19.86$ meV residual within 1 meV.  None of these is a framework
-architectural failure; all are documented convention surfaces with sibling
+to the $-19.86$ meV residual within 1 meV.  All are documented
+convention surfaces with sibling
@@ -7139,2 +7139,2 @@ entries in $\mu$H and $\mu$d.  The rest-mass projection
-(\S\ref{sec:proj_restmass}) transports the entire five-projection-chain
-skeleton bit-identically; only the nuclear charge ($Z = 1 \to 2$), nuclear
+(\S\ref{sec:proj_restmass}) carries the assembly over unchanged; only
+the nuclear charge ($Z = 1 \to 2$), nuclear
@@ -7143,7 +7143,3 @@ r_\alpha$) inputs change between the muonic catalogue entries.  The
-spinless-nucleus simplification eliminates one projection
-(\S\ref{sec:proj_magnetization_density}) entirely; the remaining four
-projections give the same architecture at all three muonic systems.
-The same architecture that closes the H 1S Lamb shift at
-\S\ref{sec:autopsy_lamb} closes the muonic Lamb shifts at
-\S\ref{sec:autopsy_muh_lamb}, \S\ref{sec:autopsy_mud_lamb}, and
-$\mu^4$He$^+$ here, with full Uehling kernel kinematics (contact-form
+spinless nucleus removes the magnetization-density channel
+(\S\ref{sec:proj_magnetization_density}) entirely.  Across the three
+muonic systems, full Uehling kernel kinematics (contact-form
@@ -7151 +7147 @@ $\to$ full kernel) and Foldy r$^2$ scaling ($r_p \to r_d \to r_\alpha$)
-being the only two ingredients that change between the muonic regimes.
+are the only two ingredients that change.
@@ -9473,2 +9469,3 @@ closed by the $\alpha^4$ Breit Layer-2 input as the
-$\mu$H Lamb at $-0.10\%$ closes by Antognini~2013 multi-loop QED
-plus deuteron polarizability plus recoil~NLO with the K\"all\'en--
+$\mu$H Lamb at $-0.10\%$ (a textbook one-loop assembly, not a
+projection chain) closes by Antognini~2013 multi-loop QED
+plus proton polarizability plus recoil~NLO with the K\"all\'en--
```

## 2026-11-08: rZG fit withdrawn, D HFS re-sourced

Owner pass (PI-approved follow-on). The Calc-rZG global Zemach-radius fit is
withdrawn (circular hydrogen Layer-2; unsourced deuterium Layer-2 budgets of
-150 and -286 ppm; unsourced +-25 ppm Eides vs Pachucki--Yerokhin drift).
D HFS re-sourced to Friar--Sick, Phys. Lett. B 579, 285 (2004) (r_Z(D) =
2.593(16) fm, Zemach -98.0(6) ppm, nuclear effects about 138 ppm) and Eides
et al., Phys. Rep. 342, 63 (2001), Eq. (414) (45.2 kHz). Verbatim before
(-) / after (+) of every cut or rewritten passage, as unified diffs with
zero context; line numbers refer to the pre-edit / post-edit files.

### papers/group6_precision_observations/paper_34_projection_taxonomy.tex

```diff
@@ -989,7 +989 @@ overshoot at LO-only).  $r_Z$ values themselves are Layer-2 inputs
-lattice or atomic spectroscopy).  The W1b recoil-mixing extension
-sprint (May 2026, see \texttt{debug/w1b\_recoil\_mixing\_extension\_memo.md})
-showed that this kernel improvement does NOT close the $\sim 0.22$~fm
-muH-alone vs Eides offset diagnosed by Sprint Calc-rZG-extended;
-that offset is now attributable to Layer-2 itemization-convention
-mismatch between Krauth~2017, Eides--Grotch--Shelyuto~2007~\cite{eides_grotch_shelyuto2007}, and Antognini-Krauth-Pohl~2015,
-not to kernel approximation.  Operator-level verified at I=1 in
+lattice or atomic spectroscopy).  Operator-level verified at I=1 in
@@ -4167,7 +4161,8 @@ Cumulative residual $-211$ ppm (corrected 2026-08-28; see
-\S\ref{sec:autopsy_d_hfs}).  The framework chain now sits $\sim 72$~ppm
-\emph{below} the Pachucki--Yerokhin~2010 / Karshenboim~2005 total theory
-($327.339$~MHz), the expected direction for a chain that omits deuteron
-polarizability (total TPE $44.5$~kHz $= +136$~ppm, QCD-internal;
-corrected 2026-08-29 from a $+44$~ppm unit slip), multi-loop QED ($\sim$few ppm,
-LS-8a wall), recoil NLO and finite-size charge --- but the itemization
-does \emph{not} close at the $\pm 25$~ppm level and is an open item.
+\S\ref{sec:autopsy_d_hfs}).  The negative sign is the expected one for a
+chain that omits deuteron polarizability: the total nuclear-structure
+effect is positive, about $+138$~ppm (Friar--Sick~2004~\cite{friar_sick2004};
+total TPE $44.5$~kHz $= +136$~ppm, QCD-internal; experiment exceeds
+nonrecoil theory by $45.2$~kHz, Eides et al.\ Phys.\ Rep.\ Eq.~(414)).
+The chain also omits multi-loop QED ($\sim$few ppm, LS-8a wall), recoil
+NLO and finite-size charge; the itemization does \emph{not} close and is
+an open item.
@@ -4187,4 +4182,2 @@ the PY~2010 budget magnitude, are \textbf{withdrawn} --- they rested on
-a baseline that was high by $496.5$~ppm).  Polarizability
-sub-component convention exposure (Friar--Payne~2005 vs Pachucki--Yerokhin
-2010+) catalogued in V.D.5 (\S\ref{sec:conv_d_polarizability}; V.D-prediction
-sprint, 2026-05-10), refining \S\ref{sec:conv_rZG_DHFS}. \\
+a baseline that was high by $496.5$~ppm).  The omitted polarizability
+channel is recorded in V.D.5 (\S\ref{sec:conv_d_polarizability}). \\
@@ -4382,2 +4375 @@ Source: \texttt{debug/hylleraas\_r12\_track1\_memo.md} \S3. \\
-Self-consistent $r_Z$ global extraction (Sprint Calc-rZG, 2026-05-09;
-Layer-2 corrected 2026-05-09; extended-fit diagnostic 2026-05-09)
+Self-consistent $r_Z$ global extraction (Sprint Calc-rZG, 2026-05-09)
@@ -4385,22 +4377,8 @@ Layer-2 corrected 2026-05-09; extended-fit diagnostic 2026-05-09)
-$\mu$H 1S HFS $+$ D 1S HFS & \emph{Original 3-obs fit} (Layer-2 corrected,
-$\sigma_{\mu H}$=1500 ppm, ``LS-8a wall scale''): $r_Z(p) = 1.180(166)$ fm
-[$+0.80\sigma$ above Eides $1.045(20)$, $+1.00\sigma$ above lattice
-$1.013(10)(12)$], $r_Z(D) = 2.583(1588)$ fm [$-0.01\sigma$ from
-Friar--Payne $2.593(16)$], $\chi^2/\text{dof} = 0.43/1$.
-\emph{Extended-fit diagnostic} (sprint Calc-rZG-extended, $\sigma_{\mu H}$
-hardened to literature kernel scale ${\sim}367$ ppm): $r_Z(p) =
-1.257(51)$~fm, $\sigma$ closes to $51$~mfm but extracted central
-value sits $+3.85\sigma$ above Eides and $+4.55\sigma$ above lattice.
-& C & \textbf{Framework $r_Z(p)$ extraction is currently kernel-limited,
-not Layer-2-itemization-limited}.  Original 3-obs fit at $\sigma=166$~mfm
-(driven by $\sigma_{\mu H}=1500$~ppm conservative) lands within $1\sigma$
-of all three literature determinations after Layer-2 budget correction
-(2026-05-09 bug fix; original rZG memo specified Layer-2$_D$ = $-150$ ppm
-based on a mis-itemization, corrected to $-286$ ppm per Pachucki--Yerokhin
-2010 + Friar--Payne 2005 closure analysis; Layer-2$_H$ tightened from
-$0$ to $-18.5$ ppm); $\chi^2/\text{dof} = 0.43/1$
-confirms internal fit consistency.  \textbf{The original $r_Z(D) = 6.18$
-fm artifact was traced to a Layer-2 budget misspecification, NOT to a
-recoil double-counting in the production cross-register $V_{eN}$
-kernel.}  Diagnostic verification
-(\texttt{debug/rzg\_bug\_diagnosis\_compute.py}): the production
+$\mu$H 1S HFS $+$ D 1S HFS & \textbf{Withdrawn}: no extracted radius
+stands & C & The global fit is withdrawn.  Its hydrogen Layer-2 term
+was the framework residual evaluated at the literature radius
+$r_Z(p) = 1.045$~fm, so the hydrogen point returned that radius by
+construction, and its deuterium Layer-2 budget had no source.  A radius
+extraction with sourced Layer-2 budgets is owed (the sourced deuterium
+inputs are collected in \S\ref{sec:autopsy_d_hfs}).  Independent of
+the fit, a code check stands: the production
@@ -4410,34 +4388,4 @@ hydrogen, $2.03\%$ for deuterium, and $8.18\%$ for muonium without
-any species-specific code path; it is species-agnostic and
-correct.  Six new D and muonium tests added to
-\texttt{tests/test\_cross\_register\_vne.py} document and protect
-this property going forward.  \textbf{Extended-fit diagnostic
-finding (2026-05-09 follow-on)}: tightening $\sigma_{\mu H}$ from
-1500~ppm to either 150~ppm (Krauth 2017 multi-loop QED itemization) or
-to 367~ppm (5\% kernel approximation in muonic systems) drives
-$\sigma(r_Z(p))$ below 50 mfm but extracts a central value
-$r_Z(p)\sim 1.18$--$1.26$~fm, INCONSISTENT with both Eides
-$1.045(20)$ and lattice $1.013(16)$ at $>3.85\sigma$.  The +0.22~fm
-offset between the muH-alone extraction ($r_Z=1.265$~fm) and Eides
-($1.045$~fm) is the structural fingerprint of recoil-mixing in the
-literature Zemach formula that the framework's leading-order Eides
-kernel $-2Zm_\text{red}r_Z$ does not yet capture
-(Sprint MH Track B \S2.4--2.5 W1b operator-level extension point).
-\textbf{Verdict: framework currently \emph{cannot} resolve the
-Eides-vs-lattice 30~mfm tension at the operator-level precision of
-the leading-order Eides Zemach kernel.}  Path to sub-30~mfm: extend
-\texttt{magnetization\_density.py} to include recoil-mixing
-(Friar--Payne 2005 \S III treatment), expected to close the +0.22~fm
-muH-alone offset and bring $\sigma(r_Z(p))$ to ${\sim}30$~mfm.
-\textbf{Genuine prospective contribution}: §III.18-based global fit
-operator-level-localizes the ``Layer-2 itemization is the load-bearing
-systematic for atomic $r_Z$ extraction'' finding --- a unified
-analytical lens not currently produced by standard atomic-physics
-compilations.  Note that the diagnostic exposes a SECOND load-bearing
-systematic: the kernel approximation order of the framework's W1b
-operator dominates over Layer-2 itemization in muonic regime
-extraction.  §III.18 itself reproduces Eides analytical
-$-2 Z\alpha m_e r_Z$ to floating-point precision.  Source memos:
-\texttt{debug/calc\_track\_rZG\_global\_zemach\_memo.md},
-\texttt{debug/rzg\_bug\_diagnosis\_memo.md} (bug fix),
-\texttt{debug/calc\_track\_rZG\_extended\_memo.md} (extended-fit
-diagnostic).  \textbf{Sprint W1b recoil-mixing extension
+any species-specific code path (\texttt{tests/test\_cross\_register\_vne.py}).
+\S III.18 itself reproduces Eides analytical
+$-2 Z\alpha m_e r_Z$ to floating-point precision.
+\textbf{Sprint W1b recoil-mixing extension
@@ -4452,14 +4400 @@ fm reproduces Krauth~2017's full-theory line to $\sim 1\%$ (vs $\sim
-3\%$ overshoot at LO).  \textbf{rZG-extended-v2 fit result}: with
-NLO kernel and Layer-2 re-anchored at the Eides reference,
-$\sigma(r_Z(p))$ tightens from 51~mfm to \textbf{16~mfm} (below the
-30~mfm Eides-vs-lattice gap), but the central value \emph{remains}
-at $r_Z(p) \approx 1.286$~fm --- only $\sim 10\%$ of the v1
-$+0.22$~fm offset closes.  \textbf{Diagnostic verdict revised}: the
-$+0.22$~fm offset is NOT kernel-approximation; it is a Layer-2
-itemization-convention mismatch between Krauth~2017 (used for muH
-Layer-2) and Eides--Grotch--Shelyuto~2007~\cite{eides_grotch_shelyuto2007} Ch.~11 (used for H~21cm Layer-2), with
-the convention difference at $\sim 0.01\%$ of $\nu_F$ propagating
-to $\sim 25$~mfm in extracted $r_Z(p)$.  Next-step recommendation:
-Layer-2 convention reconciliation sprint, NOT further kernel
-extensions.  See \S~\ref{sec:convention_exposures} for catalogue
-context (entry \S~\ref{sec:conv_rZG_DHFS}).  Source memo:
+3\%$ overshoot at LO).  Source memo:
@@ -5195,3 +5130,2 @@ V.D-prediction candidate}, class~(i) inter-compilation itemization,
-sibling to \S\ref{sec:conv_d_polarizability} (which addresses
-sub-component itemization of deuteron HFS polarizability across
-Friar--Payne~2005 and Pachucki--Yerokhin~2010+).
+sibling to \S\ref{sec:conv_d_polarizability} (which records the
+deuteron polarizability channel of the D HFS).
@@ -5597,4 +5531,3 @@ Cross-anchoring the three baselines against each other now reproduces T
-to $0.02$~ppm.  \textbf{The prior claim that the $+285.6$~ppm residual
-matched the $-286$~ppm PY 2010 Layer-2 budget magnitude is withdrawn};
-the itemization does not close at the $\pm 25$~ppm level.  What replaces
-it is a sharper structural statement.
+to $0.02$~ppm.  The corrected residual does not close against the
+sourced nuclear-structure budget; what it yields instead is a sharper
+structural statement.
@@ -5692,6 +5625,8 @@ cited Layer-2 evaluation; the $+23.0$~ppm attribution to higher-order
-QED is an inference from the H parallel, not an itemization.  The Layer-2 inputs are the scalar
-$r_Z(D)$ value; the PY 2010 budget itemization
-\textbf{[2026-08-28: listed here historically, but the chain does not
-consume it --- the only Layer-2 scalar entering is $r_Z(D)$.  This
-matters for Prediction~\ref{pred:l2_presence}: with $r_Z(D)$ alone as
-the input, the corrected $-210.9$~ppm residual is not bounded by
+QED is an inference from the H parallel, not an itemization.  The only
+Layer-2 scalar the chain consumes is $r_Z(D) = 2.593(16)$~fm, whose
+Zemach term is $-98.0(6)$~ppm (Friar--Sick~2004~\cite{friar_sick2004});
+the same source puts the total nuclear effects at about $+138$~ppm,
+consistent with the $+135.9$~ppm TPE total above and with the
+$45.2$~kHz by which experiment exceeds nonrecoil theory (Eides et al.\
+Phys.\ Rep.\ Eq.~(414)).  With $r_Z(D)$ alone as the input, the
+$-210.9$~ppm residual is not bounded by
@@ -5699,4 +5634,2 @@ $\max(\varepsilon_\text{basis}, |L_2|)$ on the ppm-class reading, so this
-row is an open case for that Prediction rather than a supporting one.]}
-($-286$~ppm convention; alternate Eides \S~12.2.2 gives $-150$~ppm,
-the source of the well-known $\sim 5$~fm extraction artifact in the
-global rZG fit, see \S\ref{sec:conv_rZG_DHFS}).  The operator
+row is an open case for Prediction~\ref{pred:l2_presence} rather than a
+supporting one.  The operator
@@ -5707,16 +5640,7 @@ autopsy.
-\textbf{Why D Layer-2 is deeper than H Layer-2.}  The cumulative
-residual was itemized per Pachucki--Yerokhin~2010 as follows.
-\textbf{[Withdrawn 2026-08-28: these entries sum to $\approx +285$~ppm,
-the retired residual.  The corrected residual is $-210.9$~ppm and the
-itemization does not close; the list is retained only as the historical
-budget it was matched against.  The $+200$~ppm polarizability entry
-below has since been re-sourced (2026-08-29): its near-coincidence with
-the $+210.9$~ppm requirement was accidental, and the primary value is
-$+336.5$~ppm with the \emph{total} TPE at $+135.9$~ppm.]}
-\begin{itemize}\setlength\itemsep{0pt}
-\item Deuteron polarizability ($\sim +200$~ppm) is the dominant
-component, $\sim 30\times$ the proton's contribution.
-\textbf{[Re-sourced 2026-08-29: the primary value is
-$+110.16(38)$~kHz $= +336.5$~ppm (Bonilla \emph{et al.},
-arXiv:2508.18776); the $\sim +200$ figure is neither that nor the
-$+135.9$~ppm total.  Dominance and sign stand; the magnitude does not.]}  This is
+\textbf{Why D Layer-2 is deeper than H Layer-2.}  The sourced budget is
+the one above: Zemach $-98.0(6)$~ppm at $r_Z(D) = 2.593(16)$~fm and
+total nuclear effects of about $+138$~ppm
+(Friar--Sick~2004~\cite{friar_sick2004}; Eides et al.\ Phys.\ Rep.\
+Eq.~(414)), with the deuteron polarizability channel ($+110.16(38)$~kHz
+$= +336.5$~ppm, Bonilla \emph{et al.}, arXiv:2508.18776) dominant and
+opposite in sign to the elastic terms.  Polarizability is
@@ -5725,15 +5649,5 @@ fundamentally because the deuteron is a weakly bound n+p system
-with much larger spatial extent ($r_Z(D) = 2.5\,r_Z(p)$).
-\item Multi-loop QED ($\alpha^2(Z\alpha)$, $\sim +40$~ppm) is the
-LS-8a renormalization gap, $\sim 7\times$ H scaling because of the
-larger nuclear extent's $\langle r\rangle$ couplings.
-\item Recoil NLO beyond reduced-mass ($\sim +30$~ppm) is the W1a-D
-Roothaan kernel-level recoil-mixing wall.
-\item Finite-size charge ($\sim +10$~ppm) and higher Friar moments
-($\sim +5$~ppm) are Layer-2 inputs via \S\ref{sec:proj_charge_density}
-and \S\ref{sec:proj_magnetization_density} NLO opt-in.
-\item Zemach NLO recoil-mixing $m_e/(m_e+m_d) \approx 2.7 \times 10^{-4}$
-is structural noise in the electronic regime ($+0.028$~ppm).
-\item Convention drift ($\pm 25$~ppm; Eides \S~12.2.2 vs Pachucki--Yerokhin
-2010) is the largest pre-existing literature-itemization
-sensitivity, documented in \S\ref{sec:conv_rZG_DHFS}.
-\end{itemize}
+with much larger spatial extent ($r_Z(D) = 2.5\,r_Z(p)$).  The Zemach
+NLO recoil-mixing factor $m_e/(m_e+m_d) \approx 2.7 \times 10^{-4}$
+is structural noise in the electronic regime ($+0.028$~ppm).  The
+remaining omitted components (multi-loop QED, recoil beyond reduced
+mass, finite-size charge) are not itemized here from a sourced budget.
@@ -5747,2 +5661,3 @@ adding one qubit and 7 Pauli terms to the I·S encoding), (c) the
-Layer-2 budget magnitude ($-286$ ppm vs $-18$ ppm, scaling with the
-nuclear extent and complexity).  The operator-level autopsy verifies
+nuclear-structure budget (about $+138$~ppm in total for D, sign-reversed
+by polarizability, against a Zemach-dominated negative correction for
+H, scaling with the nuclear extent and complexity).  The operator-level autopsy verifies
@@ -5758,2 +5673,2 @@ sentence below refers only to the itemization axis surveyed here.]}
-No new convention mismatch
-surfaced beyond the pre-existing \S\ref{sec:conv_rZG_DHFS} exposure.
+No convention mismatch
+surfaced on that axis.
@@ -7601 +7516 @@ compilations (Eides--Grotch--Shelyuto~2007~\cite{eides_grotch_shelyuto2007}, Kar
-Pachucki--Yerokhin~2010, Antognini~2013, Friar--Payne~2005)
+Pachucki--Yerokhin~2009, Antognini~2013)
@@ -7605,2 +7520,2 @@ inputs in their native conventions cannot expose these mismatches
-structurally.  This subsection consolidates the nine exposures
-catalogued to date (V.D.1--V.D.4 from the 2026-05-09 multi-track
+structurally.  This subsection consolidates the eight exposures
+catalogued to date (V.D.2--V.D.4 from the 2026-05-09 multi-track
@@ -7616 +7531 @@ Convention difference~$|$~Magnitude~$|$~Sprint reference $+$ Closure verdict.}
-\caption{Literature convention exposures catalogued to date (V.D.1--V.D.6
+\caption{Literature convention exposures catalogued to date (V.D.2--V.D.6
@@ -7629,5 +7543,0 @@ Magnitude & Sprint ref. \\
-\S\ref{sec:conv_rZG_DHFS} & D 1S HFS Layer-2 budget &
-Eides--Grotch--Shelyuto~2007~\cite{eides_grotch_shelyuto2007} \S~12.2.2 vs Pachucki--Yerokhin~2010 &
-itemization of recoil NLO + multi-loop QED + finite-size correction &
-$\sim 0.01\%$ of $\nu_F$ $\to$ $\sim 25$ mfm in $r_Z(p)$ &
-Calc-rZG (\S\ref{sec:offprecision}) \\
@@ -7653,4 +7563,4 @@ Sprint MH / Calc-muH Lamb autopsy v1 \\
-\S\ref{sec:conv_d_polarizability} & D 1S HFS polarizability sub-component &
-Friar--Payne~2005 vs Pachucki--Yerokhin~2010+ &
-low-energy aggregate (FP) vs explicit polarizability $+$ Zemach split (PY) &
-$\sim 240$~ppm of $E_F$; $\sim 80$~kHz absolute &
+\S\ref{sec:conv_d_polarizability} & D 1S HFS polarizability channel &
+Bonilla \emph{et al.}~2025~\cite{bonilla2025} (single evaluation) &
+none established; the channel is not consumed by the chain &
+$+336.5$~ppm polarizability vs $+135.9$~ppm TPE total &
@@ -7702,80 +7611,0 @@ Sub-sprint M-Z, 2026-05-23 \\
-\subsubsection{D 1S HFS Layer-2 itemization (Eides vs Pachucki--Yerokhin)}
-\label{sec:conv_rZG_DHFS}
-
-\textbf{Observable.}  Deuterium 1S hyperfine splitting
-$\nu_\text{HFS}(D) = 327\,384\,352.522$~Hz (Wineland--Ramsey 1972),
-used in the Sprint Calc-rZG global fit (\S\ref{sec:offprecision}, row
-``Self-consistent $r_Z$ global extraction'') to extract
-$r_Z(D)$ via the framework's \S\ref{sec:proj_magnetization_density}
-operator.
-
-\textbf{Source compilations.}  Eides--Grotch--Shelyuto~2007~\cite{eides_grotch_shelyuto2007} Ch.~11 (used as the
-canonical Layer-2 budget for H 21cm in this paper) and
-Pachucki--Yerokhin~2010 (the standard atomic-physics itemization for
-deuterium HFS).
-
-\textbf{Convention difference.}  The two compilations use different
-groupings for the recoil-NLO + multi-loop QED + finite-size
-contributions to $\nu_F^D$.  Eides \S~12.2.2 itemizes the leading
-Zemach line at $-150$~ppm; Pachucki--Yerokhin~2010, when fully
-expanded with the Friar--Payne~2005 closure analysis, places the
-itemization at $-286$~ppm.  Both compilations close to the same
-total $\nu_\text{HFS}(D)$ at the precision the deuterium HFS
-measurement reports; the difference is in how sub-leading recoil and
-multi-loop contributions are aggregated.
-
-\textbf{Magnitude.}  $\sim 0.01\%$ of $\nu_F$, propagating to
-$\sim 25$~mfm in extracted $r_Z(p)$ via the global fit.
-
-\textbf{Sprint reference and closure.}  The original Sprint Calc-rZG
-fit (2026-05-09 morning) used the Eides $-150$~ppm itemization for
-deuterium; this gave $r_Z(D) = 6.18(132)$~fm, the well-known
-$\sim 5$~fm artifact.  The diagnostic memo
-\texttt{debug/rzg\_bug\_diagnosis\_memo.md} verified that the
-production cross-register $V_{eN}$ kernel
-(\texttt{geovac.cross\_register\_vne}) reproduces leading-order
-Bethe--Salpeter recoil at $2.86\%$ for hydrogen, $2.03\%$ for
-deuterium, and $8.18\%$ for muonium, species-agnostically and
-correctly.  The artifact came from the Layer-2 budget specification
-for D HFS, not from the framework kernel.  Re-anchoring to the
-Pachucki--Yerokhin~2010 budget at $-286$~ppm closed the extraction
-to $r_Z(D) = 2.583(1588)$~fm (within $0.01\sigma$ of Friar--Payne
-$2.593(16)$~fm).  \textbf{Closure verdict:} resolved by Layer-2
-budget re-anchoring.  The original ``recoil double-counting''
-diagnosis was a misdiagnosis; the issue was Layer-2 itemization
-specification, not a kernel bug.
-
-\textbf{Operator-level autopsy refinement (May 9 2026).}  The Track 5
-operator-level autopsy of D 1S HFS (\S\ref{sec:autopsy_d_hfs})
-operates exclusively at the Pachucki--Yerokhin~2010 / Friar--Payne~2005
-$-286$~ppm Layer-2 budget convention.  Component-by-component the
-framework reproduces (a) the $\hat{I}\cdot\hat{S}$ Hamiltonian at I=1
-on the 6-dim joint nuclear-electronic Hilbert space at machine
-precision (multiplicity ratio D/H $= 3/2$ to 12 digits); (b) the
-\S\ref{sec:proj_magnetization_density} operator-level Zemach at
-$r_Z(D) = 2.593$~fm to $\sim 1.45 \times 10^{-14}\%$ of the LO shift;
-and (c) the cumulative chain residual, corrected on 2026-08-28 to
-$-210.9$~ppm.  \textbf{The previously reported match between a
-$+285.6$~ppm residual and the $-286$~ppm PY 2010 Layer-2 budget
-magnitude is withdrawn}: it rested on a Bohr--Fermi baseline that was
-high by $496.5$~ppm (\S\ref{sec:autopsy_d_hfs}).  With the corrected
-baseline the residual changes sign and the budget itemization does not
-close; whether it closes under either compilation's convention is now
-an open item.  No new convention mismatch surfaced at the operator level.  This sharpens the $\sim 0.01\%$ /
-$\sim 25$~mfm magnitude estimate to: framework-side architecture is
-identical for H and D (same operators, same chains, same scope
-boundaries); the Layer-2 budget magnitudes differ by $\sim 16\times$
-because the deuteron's nuclear-structure complexity (NN binding
-dynamics inside a weakly bound n+p system) makes its polarizability,
-multi-loop QED, and recoil NLO contributions much larger than the
-proton's at the same precision.  The convention exposure is
-intrinsically a Layer-2 itemization split, NOT a framework
-architecture issue.
-
-\textbf{Cross-references.}  Originating documentation:
-\S\ref{sec:offprecision} Sprint Calc-rZG row.  Operator-level
-elaboration: \S\ref{sec:autopsy_d_hfs}.  Source memos:
-\texttt{debug/calc\_track\_rZG\_global\_zemach\_memo.md},
-\texttt{debug/rzg\_bug\_diagnosis\_memo.md},
-\texttt{debug/d\_hfs\_autopsy\_track5\_memo.md}.
-
@@ -7889 +7719 @@ Source memo: \texttt{debug/muH\_Lamb\_autopsy\_v1\_memo.md}
-\subsubsection{D 1S HFS polarizability sub-component (Friar--Payne vs Pachucki--Yerokhin)}
+\subsubsection{D 1S HFS polarizability channel}
@@ -7894,41 +7724,18 @@ $\nu_\text{HFS}(D) = 327\,384\,352.522$~Hz, with focus on the
-\emph{deuteron polarizability} sub-component of the Layer-2 budget
-entering the framework's $-210.9$~ppm chain residual
-(\S\ref{sec:autopsy_d_hfs} Track~5 autopsy).  This refines the broader
-\S\ref{sec:conv_rZG_DHFS} entry to the dominant polarizability piece.
-
-\textbf{Source compilations.}  Friar--Payne~2005 (PRC \textbf{72},
-014002; building on Friar~1979 PRC \textbf{20}, 325) versus
-Pachucki--Yerokhin~2010 (the modern itemization framework adopted in
-Kalinowski--Pachucki--Yerokhin~2018 PRA \textbf{98}, 062513 for
-muonic D and applied analogously for electronic D), and the recent
-Bonilla et al.~2025 chiral-EFT deuteron nuclear-structure update
-(arXiv:2508.18776).
-
-\textbf{Convention difference.}  Friar--Payne~2005 aggregates the
-``low-energy nuclear-structure $+$ polarizability $+$ dispersion $+$
-Friar-moment recoil'' contributions into a single low-energy
-aggregate term (e.g., $\delta E_\text{Low}(eD) = 87.3$~kHz in
-Friar--Payne's reporting style).  Pachucki--Yerokhin~2010+ separates
-polarizability ($\sim +240 \times 10^{-6}\,E_F \approx +78.6$~kHz)
-and Zemach ($\sim -100 \times 10^{-6}\,E_F \approx -32.7$~kHz)
-explicitly as distinct itemized lines, using chiral-EFT-derivable
-nuclear-structure operators rather than Friar--Payne effective
-Friar moments.  Both compilations close to similar bottom-line
-totals at the experimental precision (Wineland--Ramsey~1972 $\sim$few
-ppb), but per-component itemization differs structurally.  The
-disagreement reaches operational consequences: Bonilla et al.
-2024/2025 documents a $5\sigma$ disagreement with previous
-calculations on muonic D HFS, while the underlying convention
-difference applies to electronic D HFS as well.
-
-\textbf{Magnitude.}  $\sim 240$~ppm of $E_F$ for the polarizability
-term in the PY itemization \textbf{[note 2026-08-29: the chiral-EFT
-evaluation of Bonilla \emph{et al.}\ gives the polarizability channel
-as $+110.16(38)$~kHz $= +336.5$~ppm and the \emph{total} TPE as
-$44.5(1.1)$~kHz $= +135.9$~ppm; the spread between these itemizations
-is exactly the convention exposure this entry documents, and the
-\emph{total} is the convention-independent quantity]}; the convention split (bundled-vs-separated) is
-$\sim$tens-to-hundreds of ppm of the observable; $\sim 80$~kHz
-absolute on the 327~MHz observable.  Comparable to V.D.1's
-$\sim 25$~mfm in $r_Z(p)$ precedent; both at the $10^{-4}$ to
-$10^{-3}$ relative scale.
+\emph{deuteron polarizability} channel of the nuclear-structure budget
+behind the framework's $-210.9$~ppm chain residual
+(\S\ref{sec:autopsy_d_hfs} Track~5 autopsy).
+
+\textbf{Sources.}  Bonilla \emph{et al.}~2025~\cite{bonilla2025}
+chiral-EFT deuteron nuclear-structure evaluation (arXiv:2508.18776);
+Friar--Sick~2004~\cite{friar_sick2004} for the Zemach radius and the
+total nuclear effects; Eides et al.\ Phys.\ Rep.\ Eq.~(414) for the
+difference between experiment and nonrecoil theory.
+
+\textbf{Convention difference.}  None is established: no second
+sourced itemization of this channel for electronic D HFS is available
+here to compare against.  What is sourced is a single decomposition:
+polarizability $+110.16(38)$~kHz $= +336.5$~ppm against a TPE total
+of $44.5(1.1)$~kHz $= +135.9$~ppm (Bonilla \emph{et al.}), the total
+consistent with the $\approx +138$~ppm of
+Friar--Sick~2004~\cite{friar_sick2004} and with the $45.2$~kHz of
+Eides et al.\ Phys.\ Rep.\ Eq.~(414).
@@ -7937,3 +7744,2 @@ $10^{-3}$ relative scale.
-(\texttt{debug/v\_d\_prediction\_sprint\_memo.md} Candidate~A); refines
-\S\ref{sec:conv_rZG_DHFS} to the polarizability-specific sub-component,
-addressing the dominant Layer-2 budget item.  \textbf{Closure verdict \textnormal{(corrected 2026-08-29)}:}
+(\texttt{debug/v\_d\_prediction\_sprint\_memo.md} Candidate~A).
+\textbf{Closure verdict \textnormal{(corrected 2026-08-29)}:}
@@ -7944,2 +7750,2 @@ $-210.9$~ppm and why this row is an open case for
-Prediction~\ref{pred:l2_presence}.  The convention choice documented
-here therefore does not yet propagate into the residual; it becomes
+Prediction~\ref{pred:l2_presence}.  The channel's itemization therefore
+does not yet propagate into the residual; it becomes
@@ -7948,4 +7754 @@ the itemization spread ($+336.5$~ppm inelastic vs the $+135.9$~ppm
-total) is the leading systematic on the input.  Currently anchored to Pachucki--Yerokhin~2010+
-itemization (matching the framework's Track~5 autopsy assumption); the
-Friar--Payne~2005 convention would shift the budget bucketing without
-changing bottom-line totals.
+total) is the leading systematic on the input.
@@ -7953,3 +7756 @@ changing bottom-line totals.
-\textbf{Cross-references.}  Originating documentation: this sprint memo.
-Refines \S\ref{sec:conv_rZG_DHFS} (broader D HFS itemization) to the
-polarizability-specific axis.  Cross-references
+\textbf{Cross-references.}  Cross-references
@@ -7959,9 +7759,0 @@ paragraph (deuteron polarizability is the dominant component there).
-\textbf{Honest caveat: overlap with V.D.1.}  V.D.1 documents the
-broader Eides \S~12.2.2 vs Pachucki--Yerokhin~2010 D HFS itemization
-mismatch.  V.D.5 refines V.D.1 to the polarizability-specific
-sub-component on a distinct convention axis (Friar--Payne~2005 vs
-Pachucki--Yerokhin~2010+).  The two entries differ in which
-compilation pair they compare; kept separate per V.D protocol since
-the convention axes differ.  Future audit may merge V.D.5 into V.D.1
-as a sub-paragraph if the two axes prove operationally equivalent.
-
@@ -8078 +7870 @@ distinct from \S\ref{sec:conv_d_polarizability} (deuteron HFS
-polarizability sub-component, Friar--Payne vs Pachucki--Yerokhin):
+polarizability channel, no convention claim):
@@ -8080,2 +7872 @@ this entry is between two QCD nuclear-structure evaluations of
-the deuteron polarizability rather than between atomic-physics
-itemization frameworks.  Cross-references the
+the deuteron polarizability.  Cross-references the
@@ -8155 +7946 @@ architecture limitation}, distinct from class~(i) inter-compilation
-itemization (entries V.D.1, V.D.2, V.D.3, V.D.5, V.D.6, V.D.7)
+itemization (entries V.D.2, V.D.3, V.D.6, V.D.7)
@@ -8688,2 +8479,2 @@ The exposures partition into \emph{four} structural classes.
-\S\ref{sec:conv_rZG_DHFS}, \S\ref{sec:conv_muH_VP},
-\S\ref{sec:conv_karshenboim_recoil}, \S\ref{sec:conv_d_polarizability},
+\S\ref{sec:conv_muH_VP},
+\S\ref{sec:conv_karshenboim_recoil},
@@ -8692,3 +8483,3 @@ The exposures partition into \emph{four} structural classes.
-six of the nine exposures (one of them, the muonium recoil
-aggregation, now withdrawn) are between Friar--Payne~2005, Eides--Grotch--Shelyuto~2007~\cite{eides_grotch_shelyuto2007}, Krauth~2017, Karshenboim~2005,
-Drake~1990 / Drake--Yan~1995, Pachucki--Yerokhin~2009/2010,
+four of the eight exposures (one of them, the muonium recoil
+aggregation, now withdrawn) are between Krauth~2017, Karshenboim~2005,
+Drake~1990 / Drake--Yan~1995, Pachucki--Yerokhin~2009,
@@ -8743,2 +8534 @@ itemization choices that single-observable analyses cannot see.
-(D HFS polarizability sub-component) is HFS-class, consistent with the
-original four entries.  V.D.6 (He $2{}^3P$ relativistic Bethe-log) is
+(D HFS polarizability channel) is HFS-class.  V.D.6 (He $2{}^3P$ relativistic Bethe-log) is
@@ -8748 +8538 @@ observables with multi-component Layer-2 decomposition'' rather than
-``HFS-only.''  All six entries remain in the same general category:
+``HFS-only.''  All five entries (V.D.2--V.D.6) remain in the same general category:
@@ -8774 +8564 @@ exposures.} % Source: synthesis pass §3.4 Pattern C predictions.
-The four entries above all involve hyperfine-class or muonic-Lamb
+The original entries all involve hyperfine-class or muonic-Lamb
@@ -8779 +8569 @@ polarizability, hadronic VP), each compiled differently across
-Eides--Grotch--Shelyuto~2007~\cite{eides_grotch_shelyuto2007}, Karshenboim~2005, Krauth~2017, Pachucki--Yerokhin~2010,
+Eides--Grotch--Shelyuto~2007~\cite{eides_grotch_shelyuto2007}, Karshenboim~2005, Krauth~2017
@@ -8788,2 +8578,4 @@ polarizability contribution to D HFS.  Expected magnitude: tens of
-ppm.  \textbf{Verdict (V.D-prediction sprint, 2026-05-10):
-CONFIRMED $\to$ V.D.5 (\S\ref{sec:conv_d_polarizability}).}
+ppm.  \textbf{Verdict: NOT CONFIRMED} --- no Pachucki--Yerokhin
+itemization of D HFS could be sourced (their 2009--2010 papers treat
+helium fine structure); the sourced polarizability channel is recorded
+in V.D.5 (\S\ref{sec:conv_d_polarizability}) without a convention claim.
@@ -8820 +8612 @@ The protocol classifies entries by predicted observable class:
-class-(i) literature-itemization (V.D.1, V.D.2, V.D.3, V.D.5, V.D.6
+class-(i) literature-itemization (V.D.2, V.D.3, V.D.6
@@ -8826,2 +8618,2 @@ should tag which class the new entry falls under.
-6/6 of the \emph{then-catalogued} \S~V.D entries (V.D.1--V.D.6; the
-2026-08 catalogue has nine, and the later architecture/wall-class
+5/5 of the \emph{then-catalogued} \S~V.D entries (V.D.2--V.D.6; the
+2026-08 catalogue has eight, and the later architecture/wall-class
@@ -8829 +8621 @@ additions are not spectroscopy observables, so Pattern~C is a statement
-about the original six) are precision-spectroscopy class observables
+about these five) are precision-spectroscopy class observables
@@ -8831 +8623 @@ about the original six) are precision-spectroscopy class observables
-\textbf{CONFIRMED at n=6 and STRENGTHENED}; the pattern broadens
+\textbf{CONFIRMED at n=5 and STRENGTHENED}; the pattern broadens
@@ -10417,0 +10210 @@ Friar--Payne~2005~\cite{friar_payne2005},
+Friar--Sick~2004~\cite{friar_sick2004},
@@ -10431,9 +10224,3 @@ that the autopsies of \S\ref{sec:autopsies} actually consume.
-still lack bibliography entries.  Two deserve naming because they are
-load-bearing and ambiguous rather than merely unlisted.  (i)~The
-``Pachucki--Yerokhin~2010'' itemization invoked for the deuterium
-hyperfine budget: that author pair's 2009/2010 papers are helium
-fine-structure calculations, and the modern deuterium-HFS itemization is
-Kalinowski--Pachucki--Yerokhin~2018 (and now
-Bonilla~\emph{et al.}~2025~\cite{bonilla2025}), so the attribution needs
-re-sourcing rather than a bibitem transcribed from it.  (ii)~The
-$-(11/48)$ positronium Breit coefficient: the section anchor is now
+still lack bibliography entries.  One deserves naming because it is
+load-bearing and ambiguous rather than merely unlisted: the
+$-(11/48)$ positronium Breit coefficient.  Its section anchor is now
@@ -10700,0 +10488,5 @@ Phys.\ Rev.\ C \textbf{72}, 014002 (2005).
+\bibitem{friar_sick2004}
+J.~L. Friar and I. Sick,
+``Zemach moments for hydrogen and deuterium,''
+Phys.\ Lett.\ B \textbf{579}, 285 (2004), nucl-th/0310043.
+
```

### papers/synthesis/group6_precision_observations_synthesis.tex

```diff
@@ -440,2 +440,2 @@ differences between independent precision compilations
-(Eides, Karshenboim, Krauth, Pachucki--Yerokhin, Antognini,
-Friar--Payne) at the sub-percent level the framework operates at.  Nine
+(Eides, Karshenboim, Krauth, Pachucki--Yerokhin, Antognini)
+at the sub-percent level the framework operates at.  Eight
@@ -443,3 +443 @@ exposures are catalogued (the later entries include architecture-class
-exposures beyond the convention class), including:\ the $\sim0.01\%$-of-$\nu_F$ recoil/multi-loop
-itemization difference that propagates to $\sim25$~mfm in an extracted
-Zemach radius; the Friar-moment convention difference between an RMS
+exposures beyond the convention class), including:\ the Friar-moment convention difference between an RMS
```

### tests/test_cross_register_vne.py

```diff
@@ -928,2 +928,3 @@ class TestSprintMHTrackDDualExpansion:
-# extraction was traced to the rZG global fit's Layer-2 budget for D HFS
-# (-150 ppm specified vs ~-286 ppm correct), NOT to a bug in this module.
+# extraction came from the rZG global fit's Layer-2 budgets, NOT from a bug
+# in this module. Those budgets (-150 and -286 ppm for D HFS) were
+# unsourced, and the fit is withdrawn (Paper 34).
```
