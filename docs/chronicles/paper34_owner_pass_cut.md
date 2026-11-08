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
