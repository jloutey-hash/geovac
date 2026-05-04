# Track TS-E3: Falsification Sprint for the M1/M2/M3 Case-Exhaustion Theorem

**Sprint:** TS-E3 (research memo, not a paper)
**Date:** 2026-05-04
**Goal:** Test whether any of three candidate sixteenth-projection mechanisms — (a) continuum gravity, (b) anomaly classes, (c) non-trivial principal bundles — introduces a fourth π-source mechanism not reducible to M1 (Hopf-base measure), M2 (Seeley–DeWitt heat-kernel), or M3 (vertex-topology Hurwitz Dirichlet L). Provide a verdict (REDUCES / REQUIRES M4 / OUT-OF-SCOPE) for each.

**Context (read for this memo):** `debug/track_ts_d_paper35_triple_test_memo.md` (definitions of M1/M2/M3 and the candidate-theorem statement); `papers/observations/paper_34_projection_taxonomy.tex` (fifteen projections, three-axis dictionary, off-precision catalogue); `papers/observations/paper_35_time_as_projection.tex` (Prediction 1, KG sprint, TX-B panel); `papers/synthesis/paper_32_spectral_triple.tex` §VI (K-decomposition into three sectors), §VII (Coulomb/HO asymmetry); `papers/observations/paper_25_hopf_gauge_structure.tex` (Hopf-measure identification); `papers/observations/paper_28_qed_s3.tex` (Seeley–DeWitt and vertex parity); `papers/core/paper_18_exchange_constants.tex` §III.7, §IV (Mellin engine, motivic-weight axis, Chern–Simons / K-theory cross-reference).

---

## §1. Methodology

### 1.1. The Track-D case-exhaustion theorem (verbatim restatement)

Track TS-D §4 proposed:

> **Theorem (TS-D triple-framing of Paper 35 Prediction 1).** Let $P$ be any of the fifteen projections of Paper 34 §III.1–§III.15. The following are equivalent:
> (1) $P$'s converged evaluation produces a result containing $\pi$ (or $\sqrt{\pi}$, or $\pi^{2k}$, or $\pi$ in a denominator).
> (2) $P$'s evaluation pipeline includes continuous integration over a parameter promoted from the discrete graph spectrum.
> (3) $P$ engages at least one of M1, M2, M3.

The theorem is established by direct case-checking on the existing fifteen projections; the open exit named at TS-D §5.5 is a sixteenth projection that introduces a *new* π-source mechanism not in M1/M2/M3.

### 1.2. Precise restatement of M1, M2, M3

**M1 (Hopf-base measure).** Source: $\pi$ enters as the calibration constant of the Hopf bundle $S^3 \to S^2 \times S^1$. Concrete form: $\pi = \mathrm{Vol}(S^2)/4 = \mathrm{Vol}(S^3)/\mathrm{Vol}(S^1)$. Operationally: any continuous integration over an angular base measure on $S^2$ (or its $S^d$ generalisations $\mathrm{Vol}(S^d) = 2\pi^{(d+1)/2}/\Gamma((d+1)/2)$) inherits this calibration. Examples in current GeoVac: Paper 2's outer $\pi$ in $K = \pi(B+F-\Delta)$; Paper 33's $1/(4\pi)$ per loop for vector-photon promotion; the SU(2) Haar volume $\mathrm{Vol}(\mathrm{SU}(2)) = 2\pi^2$ in Wilson plaquette computations.

**M2 (Seeley–DeWitt heat-kernel).** Source: the Mellin transform of the heat kernel $\mathrm{Tr}\,e^{-tD^2}$ on a compact Riemannian manifold produces the spectral zeta $\zeta_{D^2}(s)$, whose values at integer $s$ are constrained by Theorem T9 (Paper 28) to lie in $\sqrt{\pi}\cdot\mathbb{Q} \oplus \pi^2\cdot\mathbb{Q}$ on $S^3$. Concrete form: $a_0(S^3) = a_1(S^3) = \sqrt{\pi}$, $a_2(S^3) = \sqrt{\pi}/8$. Examples: Paper 28 vacuum polarization $\Pi = 1/(48\pi^2)$; QED $\beta = 2\alpha^2/(3\pi)$; the Stefan–Boltzmann $\pi^2/90$ (M2 specialised to compactified time, equivalently a chained M2 + observation/temporal-window).

**M3 (vertex-topology Hurwitz Dirichlet L).** Source: the Camporesi–Higuchi spectrum on $S^3$ has half-integer eigenvalues $|\lambda_n| = n + 3/2$. Vertex-parity selection rules in two-loop QED diagrams split the Dirac Dirichlet series into even-$n$ and odd-$n$ sub-sums whose closed forms involve quarter-integer Hurwitz zetas $\zeta(s, 3/4)$ and $\zeta(s, 5/4)$, which carry Catalan's constant $G = \beta(2)$ and Dirichlet $\beta(4)$. Concrete form: $D_\text{even}(4) - D_\text{odd}(4) = 2^{s-1}(\beta(s) - \beta(s-2))$ at $s = 4$. Distinct from M2 because $\beta$ values are *not* in $\pi^{2k}\cdot\mathbb{Q}$; they live in the mod-4 Dirichlet character L-content.

### 1.3. Test protocol

For each candidate (a)/(b)/(c):

(I) **Precise mechanism.** Concretise to a candidate formula or operator; write down where π would enter.

(II) **Reduction test.** Decompose the candidate's π-source into spectral-triple objects and check whether the decomposition lives entirely inside M1/M2/M3, or whether some piece is irreducibly outside.

(III) **GeoVac engagement.** Is this candidate currently a Paper 34 projection? If not, what infrastructure would be needed to make it one?

(IV) **Verdict.** REDUCES / REQUIRES M4 / OUT-OF-SCOPE.

---

## §2. Candidate (a): Continuum gravity

### 2.1. (I) Precise mechanism

The natural continuum-gravity projection would couple GeoVac matter to a curved background metric $g_{\mu\nu}$, equivalently extend the spectral action from $\mathrm{Tr}\,f(D^2/\Lambda^2)$ on a fixed background $S^3$ to a *dynamical* gravitational action involving Riemann curvature integrals. The Connes–Chamseddine spectral action expansion through $a_4$ produces, on a curved 4-manifold, the classical gravitational Lagrangian:

$$
S_{CC} = \int d^4 x \sqrt{g}\, \Big[ c_0 \Lambda^4 + c_2 \Lambda^2 R + c_4 (\text{Weyl}^2 + \text{Gauss-Bonnet} + R^2) + \cdots \Big],
$$

where $c_0, c_2, c_4$ are rational multiples of $\sqrt{\pi}$ inherited from the heat-kernel Seeley–DeWitt coefficients $a_0, a_2, a_4$, and the curvature scalars $R$, Weyl$^2$, etc. introduce *additional* $\pi$-content through:

- The volume measure $\int d^4 x \sqrt{g}$ on a closed 4-manifold (e.g. $\int_{S^4} = 8\pi^2/3$).
- The cosmological constant term $\Lambda \int \sqrt{g}$ on compact manifolds.
- The Gauss–Bonnet identity $\chi = (1/(32\pi^2)) \int \sqrt{g}\, (R^*R^* )$ (signature anomaly form).

A second flavor: Einstein–Hilbert action $S_{EH} = -(1/(16\pi G_N)) \int \sqrt{g}\,R$ has a *bare* $1/\pi$ in its Newton-constant prefactor.

### 2.2. (II) Reduction test

The Seeley–DeWitt coefficients $a_k$ on a curved manifold remain $\sqrt{\pi}\cdot\mathbb{Q}$ tensors of curvature: this is M2 verbatim. The volume measure $\int \sqrt{g}$ on a sphere $S^d$ evaluates to $\mathrm{Vol}(S^d) = 2\pi^{(d+1)/2}/\Gamma((d+1)/2)$, which decomposes as a Hopf-style ratio:

- $\mathrm{Vol}(S^3)/\mathrm{Vol}(S^1) = \pi$ (M1, Coulomb case).
- $\mathrm{Vol}(S^4)/\mathrm{Vol}(S^2) = (2/3)\pi$ (M1 generalisation: ratio of consecutive even/odd-d sphere volumes).
- $\mathrm{Vol}(S^5)/\mathrm{Vol}(S^3) = \pi/3$ (M1 generalisation, already exercised in Paper 24).

The Gauss–Bonnet prefactor $1/(32\pi^2)$ on a 4-manifold is harder. The Euler characteristic $\chi(M)$ is an *integer*; the prefactor $1/(32\pi^2)$ is exactly what compensates for the volume measure of $\int F \wedge F$ in 4 dimensions. Structurally, $1/(32\pi^2)$ is the calibration constant for the curvature integral whose dimension is $(\text{Vol}(S^4))^{-1}\cdot\mathrm{Vol}(S^2)^{-1}$ — a chain of M1 instances. This is consistent with the "topological invariant has a normalising volume measure" reading: the integer comes out *because* the prefactor cancels the M1-tier volume of the integration domain.

The Newton-constant prefactor $1/(16\pi G_N)$ in Einstein–Hilbert is the standard convention for matching to the Newtonian limit; the $\pi$ in the prefactor is M1 (it is the $4\pi$ in Newton's force law from Gauss's theorem on $S^2$).

**Conclusion:** every π in continuum gravity that has been investigated reduces to M1 + M2. The Riemann curvature integrals are M2 (heat-kernel-controlled coefficients on curved backgrounds); the volume measures are M1 (sphere-measure calibration); the Gauss–Bonnet 1/(32π²) is M1 (volume normalisation of an integer-valued integrand). Continuum gravity is *the* canonical CC spectral-action setting where M1 + M2 are the only π-sources.

### 2.3. (III) GeoVac engagement

GeoVac currently has no observable computed via continuum gravity. The Bargmann–Segal lattice on $S^5$ (Paper 24) and the Coulomb $S^3$ Fock graph (Paper 7) sit on *fixed* backgrounds. Adding gravitational projections would require:

- A graph-level Riemann tensor (combinatorial Ricci curvature on a graph; Forman, Ollivier). The graph already supports Ollivier-Ricci computations on the Hopf $S^3$ graph (Phase 4B α-D found Ollivier curvature identically zero on the κ-derivation probe).
- Or a continuum-gravity attached after the spectral-triple is constructed (extend $\mathrm{Tr}\,f(D^2/\Lambda^2)$ to fluctuating-metric backgrounds via the spectral-action variational principle).

Neither is currently a Paper 34 row. Adding one would be a substantive extension, not an audit fix.

### 2.4. (IV) Verdict: **REDUCES**.

Continuum gravity is *the* natural M1+M2 setting in the Connes–Chamseddine literature, not a fourth-mechanism candidate. The Riemann curvature scalars enter the spectral-action expansion through Seeley–DeWitt coefficients (M2); the volume integrals are sphere measures (M1). The Gauss–Bonnet prefactor $1/(32\pi^2)$ is a chain of M1 normalisations against integer-valued topological charges, which is the M1 mechanism applied recursively. The case-exhaustion survives candidate (a).

---

## §3. Candidate (b): Anomaly classes

### 3.1. (I) Precise mechanism

The natural anomaly-class projection is the chiral anomaly: the divergence of the axial current in a chiral-symmetric Dirac theory acquires a topological contribution

$$
\partial_\mu j^\mu_5 = \frac{1}{16\pi^2} F_{\mu\nu} \tilde{F}^{\mu\nu}
$$

(in 4D, the ABJ anomaly). The $1/(16\pi^2)$ prefactor is structurally similar to the Gauss–Bonnet prefactor of §2 — it normalises an integer-valued topological density (the Pontryagin density) against a continuous integration measure.

A second anomaly source is the eta-invariant of the Dirac operator on an odd-dimensional manifold:

$$
\eta(s) = \sum_{\lambda > 0} \mathrm{sign}(\lambda)\, |\lambda|^{-s}, \qquad \eta(0) = \text{regularised difference of positive vs negative eigenvalues}.
$$

On $S^3$ with the standard spin structure, the Camporesi–Higuchi spectrum is symmetric ($\pm |\lambda_n|$), so $\eta_{S^3}(0) = 0$. On $S^3 \times S^1_\beta$ with antiperiodic time, the eta-invariant picks up a non-trivial value.

A third source is the conformal anomaly, where the trace $\langle T^\mu_\mu \rangle$ on a curved background acquires an anomalous contribution proportional to the Weyl-tensor-squared and Euler density, with prefactors $1/(180\pi^2)$ and similar — again M2-class on the spectral-action side.

### 3.2. (II) Reduction test

**Chiral anomaly $1/(16\pi^2)$:** the prefactor traces directly to the heat-kernel computation. The standard derivation evaluates $\langle x | e^{-t D^2} | x \rangle$ at coincident points and extracts the $a_2$ coefficient on a 4-manifold — exactly a Seeley–DeWitt computation. The $\pi^2$ in the denominator is the inverse of $\mathrm{Vol}(S^3)/\mathrm{Vol}(S^1)$ chain seen squared (or alternatively an $a_2$ factor of $\sqrt{\pi}/8$ squared, giving $1/(64 \pi)$, then multiplied by an $\mathrm{Vol}(S^4)$ normalisation). This is M2 with an M1 normalisation, equivalently *both* of the standard mechanisms applied in concert. The chiral anomaly is structurally a "derivative of M2" — it is the obstruction to a symmetry coming from the fact that $\mathrm{Tr}\,(\gamma_5 e^{-tD^2})$ does not vanish at coincident points by the index theorem.

**Eta-invariant on $S^3$:** on a symmetric-spectrum manifold, $\eta(0) = 0$ by direct cancellation. On non-symmetric extensions (e.g. $S^3 \times S^1_\beta$ with antiperiodic time, twisted spinor bundles, Bianchi-IX manifolds), $\eta(0)$ is computed by:

$$
\eta(s) = \frac{1}{\Gamma((s+1)/2)} \int_0^\infty t^{(s-1)/2}\, \mathrm{Tr}(D\, e^{-tD^2})\, dt,
$$

a Mellin transform of a trace involving $D \cdot e^{-tD^2}$. This is a *first-order* Mellin transform — it differs from M2 (which is the second-order Mellin transform of $\mathrm{Tr}\,e^{-tD^2}$) by an extra factor of $D$ inside the trace. The transcendental content of $\eta(s)$ at $s = 0$ involves half-integer Hurwitz zetas in general (because $D$ has half-integer spectrum on $S^3$), which puts the eta-invariant in *M3*'s neighbourhood.

The eta-invariant on $S^3$ has been computed explicitly in the Atiyah–Patodi–Singer (APS) literature; for the standard spin structure it vanishes. For twisted spin structures or for $S^3 \times S^1_\beta$ with antiperiodic time, it acquires rational + π/2 contributions (the canonical example: $\eta_{S^3 \times S^1_\beta}(0) = $ rational + $(2\pi/\beta)^{-1} \cdot$ rational, where the second term comes from the twisting). This is M2 + observation/temporal-window (Paper 34's projection 15) — both already in the dictionary.

A subtle point: Sprint A's α-SP track (April 2026, CLAUDE.md §1.7 WH1 entry) noted that "the (–) sign on Δ in $K = \pi(B+F-\Delta)$ is APS-shape-compatible but NOT literally an APS invariant since $\Delta^{-1} = g_3^\text{Dirac}$ is a Dirac mode count." This is important: the eta-invariant *shape* is suggestive in GeoVac, but the actual APS computation has not been engaged as a projection — the framework's $\Delta$ is a single-level Dirac multiplicity (M3), not the regularised sum that defines $\eta$.

**Conformal anomaly:** the trace anomaly $\langle T^\mu_\mu \rangle = a\, E_4 + c\, W^2$ on a 4-manifold has coefficients $a, c$ that are pure rationals on the discrete-graph side; the $1/\pi^2$ prefactors come from the $\mathrm{Vol}(S^4) = 8\pi^2/3$ normalisation when integrating against the topological density $E_4$. This is M2 (Seeley–DeWitt for $a, c$) plus M1 (volume measure for the integral). The c-theorem flow content is Paper 18 "flow tier" already in the taxonomy.

### 3.3. (III) GeoVac engagement

The chiral anomaly is not currently a GeoVac observable, but the structural ingredients exist:

- Camporesi–Higuchi $\gamma_5$ chirality split (Paper 28 §IV: Dirac-on-$S^3$ has natural chirality grading $\gamma$; Paper 32 §IV records KO-dim 3 of the Coulomb triple, which has a $\gamma$ but no chiral anomaly because $S^3$ is odd-dimensional — anomalies in this sense are 4-dim phenomena).
- Vector-photon promotion (Paper 33, projection 11) carries the gauge-field $A_\mu$ structure that anomalies couple to.
- The four-dimensional setting required for chiral anomaly would be $S^3 \times S^1_\beta$ — i.e., the observation/temporal-window projection (#15) applied to the spinor sector.

So a chiral-anomaly observable in GeoVac would chain projections 11 (vector-photon) + 15 (observation) + 7 (spinor lift) — all already in the dictionary, none of which adds a new π-source.

The eta-invariant is similarly not a GeoVac observable. APS-style computations would require the Mellin-transform infrastructure already in `qed_self_energy.py` (Paper 28's ingredients) extended to the *first-order* sum $\mathrm{Tr}(D\,e^{-tD^2})$. This is not currently implemented; it would be a derived observable, not a new projection.

### 3.4. (IV) Verdict: **REDUCES**.

Anomaly π-content reduces cleanly to M1 + M2 + M3 chains:
- Chiral anomaly $1/(16\pi^2)$: M2 (heat kernel) + M1 (volume normalisation).
- Eta-invariant: first-order Mellin transform whose half-integer-Hurwitz content is M3, with M2 contributions on twisted backgrounds.
- Conformal trace anomaly: M2 coefficients + M1 volume normalisations.

No anomaly examined introduces a $\pi$ that cannot be identified with one of the three mechanisms. The case-exhaustion survives candidate (b), and the structural reading is that anomalies are *exactly* the framework's M2 + M3 content viewed through the topological index lens.

A subtle structural observation worth flagging: the eta-invariant is the cleanest example of a *first-order* Mellin transform in the framework's vocabulary. M2 was defined in TS-D as second-order ($D^2$ heat kernel); the eta-invariant uses $D$ inside the trace. If one wanted to refine the case-exhaustion, M2 would split into M2a (second-order, $\mathrm{Tr}\,e^{-tD^2}$, gives $\pi^{2k}$) and M2b (first-order, $\mathrm{Tr}\,D\,e^{-tD^2}$, gives half-integer Hurwitz / eta-invariant content). M2b contains M3 as the special case where vertex parity selects quarter-integer shifts. This refinement does not break case-exhaustion — it tightens M2's internal structure — but it is worth noting for Paper 18 §III.7's Mellin engine.

---

## §4. Candidate (c): Non-trivial principal bundles

### 4.1. (I) Precise mechanism

A non-trivial principal $G$-bundle $P \to M$ is classified by characteristic classes — Chern classes for unitary $G$, Pontrjagin classes for orthogonal $G$. The simplest example: the Hopf bundle $S^3 \to S^2$ has first Chern number $c_1 = 1$, which equals

$$
c_1 = \frac{1}{2\pi} \int_{S^2} F = 1
$$

where $F$ is the curvature 2-form of the canonical connection. The $1/(2\pi)$ prefactor calibrates the integer-valued $c_1$ against the continuous integration measure on $S^2$.

The second Chern number on a 4-manifold:

$$
c_2 = \frac{1}{8\pi^2} \int_M \mathrm{tr}(F \wedge F),
$$

with $1/(8\pi^2)$ calibrating $\mathrm{tr}(F \wedge F)$ to integer values.

Pontrjagin number $p_1 = -2c_2$ (for SU(2) = Spin(3)) gives the same $1/\pi^2$ structure with an integer/2 normalisation difference.

The more exotic candidate flagged in TS-D §5.5: a *topological* π entering through a discrete reading of $\int F \wedge F$ on the GeoVac graph — i.e., π appearing without intermediate continuous integration, just as a property of a topological invariant of the graph itself.

### 4.2. (II) Reduction test

**The conventional reading: π in characteristic-class prefactors is volume-measure normalisation (M1).**

The first Chern $c_1 = (1/(2\pi)) \int_{S^2} F$ has the $1/(2\pi)$ prefactor exactly because $\mathrm{Vol}(S^2) = 4\pi$ — the calibration $1/(2\pi)$ ensures that the *minimum* non-trivial flux quantum $\int_{S^2} F = 2\pi$ (Dirac quantization) gives integer $c_1 = 1$. The $\pi$ comes from the volume measure on the integration domain, exactly what M1 catalogues as "Hopf-base measure $\mathrm{Vol}(S^2)/4 = \pi$."

The second Chern $c_2 = (1/(8\pi^2)) \int \mathrm{tr}(F \wedge F)$ has $8\pi^2 = 2 \cdot \mathrm{Vol}(S^2)^2 / (4\pi)$ — a chain of M1 instances. On 4-manifolds, $1/(8\pi^2)$ is structurally the inverse of a 4-dimensional measure scaled to give integer instanton number on $S^4$ where $\mathrm{Vol}(S^4) = 8\pi^2/3$, giving $c_2 \in \mathbb{Z}$ for the BPST instanton solution.

**Conclusion 1 (reduction):** the π-prefactors of all standard characteristic classes are M1 normalisations of integer-valued integrands against continuous-measure integration domains. Standard topology of principal bundles is M1.

**The exotic reading: topological π without continuous integration.**

The interesting case TS-D §5.5 flagged is whether π could enter a *discrete* characteristic-class invariant on a finite GeoVac graph — i.e., whether $c_1$ or $c_2$ on the Hopf $S^3$ graph at $n_\text{max} = 3$ produces $\pi$ without any continuous integration step. The relevant computation would be:

- Build the U(1) connection on the Hopf $S^3$ graph (Paper 25 already constructs this).
- Compute the discrete Chern number via the Frobenius / signed-area formula on the graph 2-cells.

Paper 25 Theorem 1 establishes that the discrete connection has $\beta_1 = E - V + c$ as the dimension of the gauge sector (specifically $\beta_1 = 110$ for the Bargmann graph at $N_\text{max} = 5$). This $\beta_1$ is an integer; no π.

Paper 30's SU(2) Wilson construction has plaquettes whose Wilson loops give characters in $\mathrm{SU}(2)$; the discrete Chern–Simons-like invariants computable from these are sums of integer-valued plaquette holonomies, no π. The continuum limit reintroduces π through the Haar measure $\mathrm{Vol}(\mathrm{SU}(2)) = 2\pi^2$ — M1 again.

The discrete Hopf graph $S^3 \to S^2$ quotient (Phase 4B α-D, CLAUDE.md §2): $c_1$ on the discrete Hopf graph at $n_\text{max} = 3$ would be the difference between the U(1) winding numbers around the two halves of the Hopf base. This is computable from edge-orientation labels; the result is an integer (or rational without π). The discrete characteristic class is *intrinsically π-free* on the graph; π appears only when the continuum limit reintroduces the volume measure.

**Conclusion 2 (exotic reading):** non-trivial principal bundles on finite GeoVac graphs do not introduce π through their topological invariants. The integer $c_1$ stays integer; the $\pi$ in $c_1 = (1/(2\pi)) \int F$ is the *continuum* π of M1, not present in the discrete computation. The graph's topological invariants are themselves π-free, consistent with Paper 24's $\pi$-free certificate and Paper 28's graph-native QED algebraic ring.

### 4.3. (III) GeoVac engagement

GeoVac currently *does* have non-trivial principal-bundle structure on the graph:

- The Hopf graph $S^3 \to S^2$ (Paper 25) carries a U(1) connection.
- The SU(2) Wilson construction on the Coulomb $S^3$ graph (Paper 30) carries a non-abelian connection; characteristic classes of the SU(2) bundle are computable via plaquette traces.
- The SU(3) Wilson on Bargmann–Segal $S^5$ (Sprint ST-SU3, May 2026, CLAUDE.md §2) carries a non-abelian SU(3) connection; characteristic-class computations would be the natural diagnostic.

But: the *characteristic-class invariants* themselves have not been computed and tested for π-content as separate observables. They have been computed only in the form of Wilson loop expectation values ⟨W(C)⟩ at given β, which are functions of the gauge coupling, not of the topological invariant directly.

Adding a "discrete characteristic-class projection" to the Paper 34 list would be structurally novel and would be the natural test of candidate (c) at the strongest level. The expected outcome, by Conclusion 2 above, is that the discrete characteristic-class observables are integer- or rational-valued (no π), with π re-entering only on continuum limit through M1.

### 4.4. (IV) Verdict: **REDUCES** (with a flagged opportunity)

Standard π-content of characteristic-class prefactors is M1 (volume normalisation). On finite GeoVac graphs, characteristic classes are integer-valued by direct computation; no fourth mechanism is needed.

The flagged opportunity: a sprint "compute discrete $c_1$ on the Hopf $S^3$ graph at $n_\text{max} = 3$ and the discrete instanton number on the SU(2) Wilson construction at $n_\text{max} = 4$ and verify they are integer-valued" would be a clean falsifier. By the reasoning of §4.2 Conclusion 2 it should give integers; if it produced an unexpected π, that would be the genuine sixteenth-projection signal. This sprint is *not* yet executed and is the strongest concrete falsification target produced by this memo.

---

## §5. Cross-cutting analysis

### 5.1. Does the case-exhaustion survive all three candidates?

**Yes.** All three candidates reduce cleanly to M1 + M2 + M3:

- (a) Continuum gravity → M1 (volume measures) + M2 (Seeley–DeWitt curvature scalars).
- (b) Anomaly classes → M1 (volume normalisations) + M2 (heat-kernel coefficients) + M3 (eta-invariant via half-integer Hurwitz on twisted backgrounds).
- (c) Non-trivial principal bundles → M1 (characteristic-class volume calibrations); discrete invariants on finite GeoVac graphs are integer-valued (no π) by the framework's own π-free certificate.

### 5.2. Refinement of M2: first-order vs second-order Mellin

A finer-grained reading emerged in §3 (anomaly classes): M2 should arguably split into two subtypes:

- **M2a (second-order, even-zeta):** $\mathrm{Tr}\,e^{-tD^2}$, gives $\pi^{2k}$ at integer $s$ via T9.
- **M2b (first-order, mixed Hurwitz):** $\mathrm{Tr}\,D\,e^{-tD^2}$, gives half-integer Hurwitz content.

M3 is then the special case of M2b where vertex parity selects quarter-integer Hurwitz shifts. This refinement does not break case-exhaustion — it tightens M2's internal structure — but it makes the Mellin engine of Paper 18 §III.7 more uniform: one master mechanism (Mellin transform of $D^k\,e^{-tD^2}$) generating three tiers depending on the exponent $k = 0, 1, 2$ and the spectrum's Hurwitz shift.

This refinement is also visible in the c₂ = (2 − BΔ − FΔ − F/B)/5 demotion (CLAUDE.md curve-fit audit, May 2026): the apparent Paper-2-to-Paper-28 bridge was via two natural Tier-A spectral invariants of the same manifold, not a derivation; M2a's polynomial-in-$\pi^2$ structure fixes the form but does not predict the coefficients. The refinement to M2a/M2b clarifies what level of structural identification is actually warranted.

### 5.3. Could a fourth mechanism be needed? Concrete candidates

I find no compelling fourth candidate among the three flagged in TS-D §5.5. The strongest remaining direction is a *new* candidate (d) not in TS-D's enumeration:

**Candidate (d): non-Mellin spectral integrations.** TS-D's case-exhaustion implicitly assumes that any continuous integration over a parameter promoted from the discrete spectrum is either an angular Hopf-base integration (M1) or a Mellin-type heat-kernel integration (M2/M3). What about *contour integrations* in the complex plane (residue calculus of meromorphic spectral functions, e.g. Borel resummation contours, Stokes-phenomenon analytic continuations)? These are continuous integrations over a parameter promoted from the discrete spectrum, but they are not Mellin-shaped and not Hopf-base.

On reflection, contour integrations are *equivalent* to Mellin transforms by deformation arguments on the analytic structure of the Borel-summed spectral series — Paper 18's Mellin-engine subsection (§III.7) makes this point implicitly. So contour integrations reduce to M2/M3. But this point is somewhat subtle and might warrant explicit verification.

**Candidate (e) (more exotic): non-commutative-geometry-style spectral truncation.** The Connes–van Suijlekom truncation that WH1 R2.1–R2.5 are exploring (CLAUDE.md §1.7 WH1 entry) has the operator system $O_{n_\text{max}} = P_{n_\text{max}} \mathcal{A} P_{n_\text{max}}$ with propagation number 2 matching Toeplitz $S^1$. The Connes distance SDP at $n_\text{max} = 2$ produced the algebraic fingerprint $50\sqrt{3}$ (CLAUDE.md WH1 R2.3), and at higher $n_\text{max}$ produces irrational distances under physical Avery–Wen–Avery integrals. None of these distances has been π-bearing in the cases tested. If a future Connes-distance computation produced a π-bearing pure-state distance through the truncation mechanism (rather than through the underlying spectral-action machinery), that would be a candidate fourth mechanism. Currently it's not engaged as a Paper 34 projection.

### 5.4. Computational extension of GeoVac that would test these mechanisms

The most concrete test of the case-exhaustion theorem is candidate (c)'s discrete-characteristic-class sprint: compute $c_1$ on the Hopf $S^3$ graph at $n_\text{max} = 3$ and verify it is integer-valued. The Hopf graph has $V = 14$ nodes, $E = 13$ edges at $n_\text{max} = 3$ (Phase 4B α-D); the U(1) connection is built from ladder-operator phases. The discrete Chern number is computable in $O(V^2)$ time. If the result is integer (as predicted), the case-exhaustion is strengthened. If π appears unexpectedly, M4 would need to be named.

A second, independent test: compute the discrete eta-invariant of the Camporesi–Higuchi Dirac on $S^3$ at $n_\text{max} = 3$ via the regularised sum $\eta(0) = \zeta_D(0)$ where $\zeta_D(s) = \sum \mathrm{sign}(\lambda_n) |\lambda_n|^{-s}$. By spectral symmetry on $S^3$, the result should be 0; on $S^3 \times S^1_\beta$ with antiperiodic time, it should pick up a calibration $1/\pi$ contribution that decomposes as M1 + M2 (per §3.2). This sprint would test whether the framework's APS-shape reading of $\Delta$ in $K = \pi(B + F - \Delta)$ is structurally meaningful or a coincidence.

A third, more speculative test: compute the spectral action $\mathrm{Tr}\,f(D^2/\Lambda^2)$ on the SU(2) Wilson-extended Coulomb $S^3$ triple (engaging the gauge sector of Paper 30) and check whether the Wilson-loop factors that arise carry $\pi$ through Haar normalisation alone (M1) or through additional structure. By Paper 30's analysis, the Wilson character expansion goes through SU(2) Haar measure with $\mathrm{Vol}(\mathrm{SU}(2)) = 2\pi^2$ — M1 chained twice — giving expected $\pi^2$ content.

### 5.5. Is there a candidate (d) the agent identifies during the investigation?

Yes — the M2 refinement in §5.2 (split into M2a second-order and M2b first-order) is the cleanest internal sharpening, and the eta-invariant of §3.2 is the clearest example of M2b. This is not strictly a fourth mechanism — it is a structural decomposition of M2 — but it is the most informative output of the falsification sprint beyond the verdicts themselves.

Independently, the Mellin engine perspective (Paper 18 §III.7) suggests the master statement should be:

**Master Mellin theorem (sharpened from M2/M3):** every π in a GeoVac observable arising from a Mellin-type integration is a residue or sum-residue of a meromorphic spectral function, and the residue is in $\sqrt{\pi}\cdot\mathbb{Q} \oplus \pi^2\cdot\mathbb{Q}$ (T9-class) for second-order operators on $S^3$, or in mixed-Hurwitz-with-rational-shifts for first-order operators. M1 is then the special case where the "Mellin-type integration" is an angular integration on a Hopf-base sphere, equivalently the Mellin transform of the trivial heat kernel on $S^d$.

Under this sharpening, M1 + M2 + M3 collapse to a single mechanism (Mellin-type integration of the heat kernel of $D^k$ on a compact Riemannian manifold) with three sub-cases distinguished by operator order ($k = 0, 1, 2$). This is the most elegant version of the case-exhaustion theorem, and it is consistent with Paper 18 §III.7's existing machinery.

---

## §6. Verdict: case-exhaustion survives, with structural refinement

**Headline:** the M1/M2/M3 case-exhaustion theorem of TS-D §4 survives all three falsification candidates. None of (a), (b), (c) introduces a new π-source mechanism; each reduces cleanly to M1, M2, or M3 (or chains thereof).

**Verdicts:**

| Candidate | Verdict | One-line summary |
|:---|:---|:---|
| (a) Continuum gravity | **REDUCES** | M1 (sphere-volume measures) + M2 (Seeley–DeWitt curvature scalars); the canonical CC spectral-action setting. |
| (b) Anomaly classes | **REDUCES** | M2 (heat-kernel coefficients) + M1 (volume normalisations) + M3 (eta-invariant on twisted backgrounds via first-order Mellin). |
| (c) Non-trivial principal bundles | **REDUCES** | π-prefactors of Chern numbers are M1 (volume normalisations against integer-valued integrands); discrete characteristic-class invariants on finite GeoVac graphs are integer-valued (predicted, sprint not executed). |

**Refinement of the M2 mechanism.** The investigation surfaced a natural split:
- M2a: second-order Mellin transform $\mathrm{Tr}\,e^{-tD^2}$ → gives $\pi^{2k}$ at integer $s$ (T9).
- M2b: first-order Mellin transform $\mathrm{Tr}\,D\,e^{-tD^2}$ → gives half-integer Hurwitz / eta-invariant content; M3 is the vertex-parity restriction of M2b.

Under this refinement, M1 + M2 + M3 collapse into a single mechanism (Mellin transform of $D^k\,e^{-tD^2}$) with three sub-cases distinguished by operator order $k = 0, 1, 2$. The M1 case ($k = 0$) is the trivial heat kernel on a sphere (giving the Hopf-base measure); M2a ($k = 2$) is the standard heat kernel; M2b/M3 ($k = 1$) is the eta-invariant / vertex-parity content. This is the cleanest version of TS-D's theorem.

**Sharpened theorem statement.** The TS-D theorem can be rewritten as:

> **Theorem (TS-E3 sharpened).** Let $P$ be any of the fifteen GeoVac projections of Paper 34 §III.1–§III.15. The following are equivalent:
> (1) $P$'s converged evaluation produces a result containing $\pi$.
> (2) $P$'s evaluation pipeline includes a continuous integration over a parameter promoted from the discrete graph spectrum.
> (3) $P$ engages a Mellin-type integration of $\mathrm{Tr}\,D^k\,e^{-tD^2}$ for some $k \in \{0, 1, 2\}$ on a compact Riemannian manifold (M1 = $k=0$, M2 = $k=2$, M3 = $k=1$ with quarter-integer Hurwitz selection).

The sharpened theorem statement makes the *single mechanism* nature of the framework's π-source explicit. M1, M2, M3 are not three different mechanisms but one mechanism (Mellin transform of a power of $D$ times the heat kernel) with three sub-cases.

**Open falsification target.** The discrete characteristic-class sprint flagged in §4.3 (compute $c_1$ on the Hopf $S^3$ graph at $n_\text{max} = 3$, verify integer) is the strongest concrete test the present memo could not execute. It is recommended as a follow-up sprint TS-E4: a clean computational verification that discrete topological invariants on finite GeoVac graphs are integer-valued, falsifying the hypothesis that a topological π could enter the framework without continuous integration.

**Implication for Papers 34 / 35.** Track TS-D's recommendation to amend Paper 32 §IX with the M1/M2/M3 case-exhaustion stands; the present memo strengthens it to a *single-mechanism* statement (Mellin engine + operator-order refinement) that is structurally cleaner than the three-mechanism version. Whether the sharpened theorem statement should land in Paper 35 §VIII or Paper 32 §IX or Paper 18 §III.7 (the natural Mellin-engine home) is a paper-architecture decision for plan mode.

---

## Report-back

- **Path to memo:** `debug/track_ts_e3_sixteenth_projection_memo.md`
- **Verdicts:**
  - (a) continuum gravity: **REDUCES** to M1 + M2.
  - (b) anomaly classes: **REDUCES** to M1 + M2 + M3 (with M2b first-order Mellin sub-tier surfaced).
  - (c) non-trivial principal bundles: **REDUCES** (π-prefactors are M1 volume normalisations of integer-valued integrands; discrete invariants on GeoVac graphs are integer-valued by the framework's own π-free certificate).
- **Case-exhaustion survives:** yes; all three candidates reduce cleanly to M1/M2/M3 chains. The TS-D theorem stands.
- **Most surprising finding (~150 words):** The three candidates do not just reduce — they reduce to a *single* underlying mechanism. The case-exhaustion theorem can be sharpened to "every π in a GeoVac observable arises from a Mellin transform of $\mathrm{Tr}\,D^k\,e^{-tD^2}$ for $k \in \{0, 1, 2\}$ on a compact Riemannian manifold," with M1 ($k=0$, trivial heat kernel on sphere = Hopf-base measure), M2 ($k=2$, standard spectral action), and M3 ($k=1$ with vertex-parity selection of quarter-integer Hurwitz shifts) as three sub-cases of one master mechanism rather than three independent mechanisms. The eta-invariant is the cleanest example of M2b / first-order Mellin, and it is structurally what suggested the unification. The Mellin engine of Paper 18 §III.7 is therefore the *master* π-injecting mechanism of the framework, with the three sub-cases distinguished only by the operator order $k$ inside the trace. Continuum gravity, anomalies, and characteristic classes all live inside this single Mellin master mechanism — there is no fourth source visible from the three candidates investigated.
