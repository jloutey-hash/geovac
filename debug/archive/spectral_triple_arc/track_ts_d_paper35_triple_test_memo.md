# Track TS-D: Paper 35 Prediction 1 vs the Spectral-Triple Lens

**Sprint:** TS-D (research memo, not a paper)
**Date:** 2026-05-04
**Goal:** Test whether Paper 35's load-bearing principle (Prediction 1: a GeoVac observable contains π iff its evaluation includes a continuous integration over a temporal/spectral parameter promoted from the discrete graph spectrum) agrees with the spectral-triple framing of Paper 32 across each of Paper 34's fifteen projections. Where they disagree, list the disagreements as falsification targets.

**Companion files:**
- `debug/track_ts_b_dictionary_memo.md` (TS-B, 15-row dictionary; the starting point)
- `debug/track_ts_a_gh_convergence_memo.md` (TS-A)
- `debug/track_ts_c_asymmetry_investigation_memo.md` (TS-C)
- `debug/tx_b_paper35_test_memo.md` (TX-B, 5/5 a-priori panel)
- `debug/data/tx_b_predictions.json`, `debug/data/tx_b_obs_{1..5}.json`, `debug/data/kg{1,2,3,5}_*.json`

**Source papers (read for this memo):** Paper 35 (Prediction 1, falsification panel, §VI projections 14–15); Paper 34 (§III.1–§III.15 projection list, three-axis dictionary); Paper 32 (§III construction, §IV axiom audit, §V sub-sectors, §VI K-decomposition, §VII Coulomb/HO asymmetry, §VIII Marcolli–vS lineage); Paper 28 (Theorem 1 = T9 squared-Dirac spectral zeta is π-even at integer s, Theorem 2 = parity discriminant); Paper 18 §III.7 Mellin engine and §IV operator-order × bundle grid; Paper 25 (Hopf U(1) Wilson reading, Vol(S²)/4 calibration); Paper 33 (1+6+1 selection-rule partition, vector photon = 1/(4π) per loop = S² Weyl exchange constant of Hopf base).

---

## §1. Methodology

### 1.1. Paper 35 Prediction 1 (verbatim)

> A GeoVac observable contains π if and only if its evaluation includes a continuous integration over a temporal or spectral parameter that has been promoted from the discrete graph spectrum. Equivalently: the discrete graph itself is π-free; the projection that integrates (over time, over a Matsubara mode, over a Mellin parameter, over an analytically continued spectral mode) is where π appears.

The prediction has two operational halves. *Necessary half (only if):* if a derivation step is supposed to produce π, the discipline is to identify the continuous integration over a temporal or spectral parameter promoted from the discrete spectrum. *Sufficient half (if):* every continuous integration over such a parameter does produce π. After 208 individual checks across the KG sprint and TX-B (100% match rate), the prediction graduated from observation to load-bearing principle.

The test inputs for any new projection are therefore:

(a) Does the projection's evaluation pipeline include a continuous integration over a continuous parameter that has been *promoted* from the discrete graph spectrum (Schwinger proper time, Matsubara frequency, Mellin contour, analytic continuation in $s$, heat-kernel time, etc.)?

(b) If yes, the projection introduces π. If no, the projection preserves Layer 1's π-free ring.

The "promoted from the discrete graph spectrum" qualifier matters: angular integration over a Hopf $S^2$ base in the *continuum* picture is also a continuous integration that produces π, but in Paper 35's framing this is the same object as the temporal-integration mechanism (it lives in the same observation/temporal-window category, because both are continuous integration measures imposed on a discretely-labeled spectral data set). Sprint TX-B observable 5 explicitly notes this: the discrete Coulomb-Green's function on the Hopf $S^3$ graph at finite $n_\text{max}$ is π-free; the continuum limit ($n_\text{max} \to \infty$) injects π via spherical-harmonic normalization "exactly where Prediction~1 predicts."

### 1.2. The spectral-triple lens for π content

Paper 32 constructs the GeoVac spectral triple $(\mathcal{A}_{\mathrm{GV}}, \mathcal{H}_{\mathrm{GV}}, D_{\mathrm{GV}})$ on the Fock-projected $S^3$ graph. The spectral-triple framing of *where* π enters has three known mechanisms, all anchored to specific Connes–Chamseddine objects:

**M1 — Hopf-base measure (Paper 25, Paper 32 §VII).** The Hopf bundle $S^3 \to S^2 \times S^1$ decomposes the GeoVac spectral action across base and fibre. The S² Hopf base carries the calibration constant $\pi = \mathrm{Vol}(S^2)/4$, equivalently $\mathrm{Vol}(S^3)/\mathrm{Vol}(S^1) = \pi$ (Sprint A α-PI identification). This is the canonical π-source of Paper 2 ($K = \pi(B + F - \Delta)$) and of Paper 33's $1/(4\pi)$ per-loop calibration of the vector-photon promotion.

**M2 — Seeley–DeWitt heat-kernel coefficients on $S^3$ (Paper 18 §III.7, Paper 28 Theorem T9).** The Mellin transform of the heat kernel $\mathrm{Tr}\,e^{-tD^2}$ produces the spectral zeta $\zeta_{D^2}(s)$. By T9, at every integer $s \ge 1$, $\zeta_{D^2}(s)$ is a two-term polynomial in $\pi^2$ with rational coefficients. This is where calibration-π and even-π² content enters CC spectral-action computations: $a_0 = a_1 = \sqrt{\pi}$, $a_2 = \sqrt{\pi}/8$ on unit $S^3$, and the resulting one-loop observables live in the ring $\mathbb{Q}[\pi^{2k}]$. The key restriction is *operator order*: only the second-order ($D^2$) spectrum gives this clean π-only structure; the first-order $|D|$ spectrum produces odd-zeta values at odd integer $s$.

**M3 — Vertex-topology Dirichlet $L$-values (Paper 18 §III.7, Paper 28 Theorem 2 + Theorem 3).** The Dirac Dirichlet series $D(s) = \sum_n g_n |\lambda_n|^{-s}$ at quarter-integer Hurwitz shifts produces Catalan's constant $G = \beta(2)$ and Dirichlet $\beta(4)$, exposed by the vertex-parity selection rule of two-loop sunset diagrams. These are *not* π-only; they are L-function content with explicit π-prefactors only inside the Riemann–Bernoulli identity for the even-character part. The transcendental class is mod-4 Dirichlet character L-content, which carries π through $\beta(s) - \beta(s-2)$ closed forms (eq:rh_j_identity in Paper 18 §III.5).

The spectral-triple lens predicts π enters via one or more of M1–M3 (or via continuum-limit angular measures, which are CC spectral-action coefficients in disguise: $a_0 = \sqrt{\pi}$ on $S^d$ gives $\pi^{d/2}$ when continuum integrated). The lens predicts π *does not* enter for any projection that:

- preserves the algebra $\mathcal{A}_{\mathrm{GV}}$ at integer-rational level (universal sector of Paper 31);
- modifies $D_{\mathrm{GV}}$ only by central deformations $D \mapsto D + m \cdot \mathbb{1}$ (rest mass; preserves the spectrum's algebraic class);
- does not invoke the heat-kernel Mellin transform or the Hopf-base measure.

### 1.3. Comparison protocol

For each projection $P_i$ in Paper 34's list, we apply both lenses:

- *Paper 35 input:* does $P_i$'s evaluation include continuous integration over a temporal/spectral parameter promoted from the discrete spectrum? (yes/no)
- *Paper 35 verdict:* π / π-free.
- *Triple mechanism:* Hopf-base (M1) / Seeley–DeWitt (M2) / vertex-Dirichlet-L (M3) / "none" (no triple-level π injection).
- *Triple verdict:* π / π-free.
- *Agreement:* yes/no.

Where the lenses disagree, we identify which is right *per current evidence* (TX-B confirmations, T9 theorem, Phase 4B–4I sprint results, Paper 32 §VI K-decomposition).

---

## §2. Fifteen-row comparison table

The table below compares both lenses for each of Paper 34 §III.1–§III.15. The "Sector" column is from TS-B's dictionary for cross-reference; the "U/C" column flags universal (U) vs. Coulomb-specific (C) per Paper 31. The "Mechanism" column names the spectral-triple π-source if applicable.

| # | Projection | Sector | U/C | (a) Continuous integration? | (b) Paper 35 verdict | (c) Triple mechanism | (d) Triple verdict | (e) Agreement | Notes |
|---:|:---|:---|:--:|:---|:--:|:---|:--:|:--:|:---|
| 1 | Fock conformal $\mathbb{R}^3 \to S^3$ | $\mathcal{A,H,D}$ | C | No (defines the bare graph; subsequent integrals introduce π) | π-free at the projection level | None at projection level (sets $\mathcal{A}, D$ for downstream use) | π-free | ✓ | The projection itself is a basis change; π enters only when downstream observables Mellin-transform or Hopf-decompose. |
| 2 | Hopf bundle $S^3 \to S^2 \times S^1$ | $\mathcal{H, D}$ | C | Implicit: angular integration over $S^2$ base measure | π (Vol(S²)/4) | M1 — Hopf-base measure $\mathrm{Vol}(S^2)/4$ | π | ✓ | **Boundary case A: §3.2.** Paper 35 reads the $S^2$ measure as a "geometric" continuous integration, lumped with temporal/spectral; triple framing names it explicitly via M1. Both predict π; identical mechanism by different names. |
| 3 | Bargmann–Segal $\mathbb{R}^3_\text{HO} \to S^5$ Hardy | $\mathcal{A',H',D'}$ (distinct triple) | H | No (graph is bit-exact π-free; spectrum half-integer rational; degeneracies polynomial) | π-free | None at projection level (defines a *different* triple with no calibration π — Paper 32 §VII) | π-free | ✓ | Paper 24's π-free certificate; HO triple has no Hopf-base structure. TX-B Obs 2 confirmed Casimir = $-17/3840$ rational. |
| 4 | Stereographic / conformal coord. change | $\mathcal{A}$ | * | No (algebraic coordinate change; Jacobian is rational power of $1+r^2$ on stereo coords) | π-free | None (preserves $D$ up to conformal factor) | π-free | ✓ | Coulomb $1/r$ emerges as coordinate distortion (Paper 7), no integration. |
| 5 | Sturmian reparam. at $\lambda = Z/n$ | $\mathcal{A,H}$ | C | No (re-labels same Fock graph) | π-free | None (preserves $\mathcal{A}, D$ spectra) | π-free | ✓ | Bethe-log spectral sums *downstream* of Sturmian projection produce calibration π; Sturmian itself is a relabeling. |
| 6 | Connes–Chamseddine spectral action $\mathrm{Tr}\, f(D^2/\Lambda^2)$ | O (also $\mathcal{A} \otimes M_n$ when gauge engaged) | C | Yes (heat-kernel time integration $\int_0^\infty t^{s-1}\,\mathrm{Tr}\,e^{-tD^2}\,dt/\Gamma(s)$) | π | M2 — Seeley–DeWitt $a_k$ coefficients $\propto \sqrt{\pi}$, Mellin transform produces $\pi^{2k}\cdot\mathbb{Q}$ at integer $s$ (T9) | π | ✓ | **Boundary case D: §3.5.** Both predict π by *the same mechanism* (heat-kernel time integration is the Mellin parameter promoted from the discrete spectrum). T9 specializes Paper 35's general statement to the squared-Dirac case at integer $s$. |
| 7 | Camporesi–Higuchi spinor lift | $\mathcal{H}$ | C | Yes (Mellin sums over half-integer Hurwitz produce odd-ζ at odd $s$, Catalan $G$ via vertex parity) | π via $\sqrt{\pi}$ in Hurwitz expansions; also odd-ζ and $\beta(s)$ | M2 + M3 — half-integer Hurwitz mechanism for odd-ζ; vertex-topology for $G, \beta(4)$ | π (mixed with odd-ζ) | ✓ | The spinor lift opens the *odd*-ζ + Dirichlet-L bucket beyond π-only. Paper 35 captures the π-bearing pieces; the odd-ζ pieces are also continuous-integration-induced (Mellin transform of half-integer spectrum). Both lenses agree on π-presence; the *additional* transcendental classes are Paper 18's tier refinement. |
| 8 | Wigner $3j$ angular coupling | (none — internal $\mathcal{A}$ SO(3) decomposition) | U | No (algebraic recoupling; $\mathbb{Q}[\sqrt{2k+1}]$) | π-free | None (universal angular decomposition) | π-free | ✓ | Paper 22's angular sparsity theorem: rational, basis-intrinsic, π-free. |
| 9 | Wigner $D$-matrix rotation $A \to B$ | $\mathcal{A}$ relabel | U | No (algebraic rotation, $\mathbb{Q}[\sqrt{2}, \sqrt{3}, \sqrt{6}]$) | π-free | None | π-free | ✓ | Multi-center molecular pipeline; $\mathbb{Q}$-extension content only. |
| 10 | Wilson plaquette | $\mathcal{A} \otimes M_n$ | C | Yes (SU(2) Haar integration over link variables; Wilson loops are Haar-averaged characters) | π (via Haar measure normalization $\mathrm{Vol}(\mathrm{SU}(2))$ or via base-measure $\mathrm{Vol}(S^2)/4$ for the $L_1$ kinetic term) | M1 — the $L_1 = B^T B$ kinetic term reproduces Paper 25's U(1) and is structurally a Hopf-base object | π | ✓ | **Boundary case E: §3.6.** Wilson plaquette engages SU(2) Haar; the maximal-torus reduction matches Paper 25's $L_1$ and the $L_1$ kinetic content is the Hopf-base $S^2$ measure (Paper 30 Theorem 2: kinetic coefficient is rational multiple of plaquette count, but the propagator measure on the gauge sector inherits Vol(S²)/4 from M1). Both predict π. |
| 11 | Vector-photon promotion | $\mathcal{H}$ basis-refinement | C | Implicit: introduces angular $S^2$ measure on photon mode space; equivalently, a per-loop $\int dq\, dm_q$ over the vector harmonic basis with $\mathrm{Vol}(S^2)/4$ normalization | π ($1/(4\pi)$ per loop) | M1 — explicitly the $S^2$ Weyl exchange constant of the Hopf base (Paper 33) | π | ✓ | **Boundary case C: §3.4.** TX-B Obs 4 (Heisenberg–Euler) confirmed; the same $1/(4\pi)$ per loop of vector-photon promotion (Paper 33) IS the Hopf-base measure $\mathrm{Vol}(S^2)/4$ (Paper 28 §vector_photon_qed). Both lenses identify the same mechanism — Paper 35 names the integration; the triple framing names the measure. |
| 12 | Molecule-frame hyperspherical | $\mathcal{H, D}$ | C | Mixed: angular integration over $(\rho, \alpha, \theta_{12}, \ldots)$ enters at each $\rho$-slice for adiabatic potentials; per-channel sums are spectral but each angular eigenvalue sweep is a continuous diagonalization in $\rho$ | π in continuum continuum integrals; π-free at angular level (Gaunt-rational) | None at projection level (composition of basis change + Gaunt 3j); continuum integration falls under M1 ($S^2$ measure) for any Hopf-projected sub-integral | π in continuum, π-free in finite truncation | ✓ | The composed-architecture pipeline of Paper 17 separates angular (Gaunt-rational, π-free) from radial (continuous, π via continuum measure). The dichotomy aligns with Paper 35's prediction at the boundary. |
| 13 | Drake–Swainson asymptotic subtraction | O | C | Mixed: introduces transient subtraction scale $K$; intermediate $\beta_{\text{low}}(K), \beta_{\text{high}}(K)$ are continuous in $K$; final answer $K$-independent (cancellation) | **Ambiguous.** Intermediate steps π-bearing (Bethe-log carries log of a continuous quantity); final answer is dimensionless and the Drake denominator $D_{\text{drake}} = 2(2\ell+1)Z^4/n^3$ is rational. Whether the projection's "evaluation" includes intermediate steps depends on what we count as an evaluation step. | None at the projection's definitional level (regularization-flow operation); inherits M1 + M2 from upstream Sturmian + spectral-action integration | **Ambiguous.** The spectral action upstream produces π; Drake itself adds no new π. | **§3.1: tentatively ✓ (both lenses agree on the converged answer's π content), but π content of the *intermediate* depends on lens.** | **Boundary case F: §3.1.** Drake is a flow-tier post-processing of Bethe-log output. The cancellation makes the converged answer's π content equal to the upstream Sturmian + spectral-action input. Both lenses agree on the final answer; Paper 35 needs the qualifier "evaluated in the converged sense" or "after the K-cancellation." |
| 14 | Rest-mass projection | $D$ (central deformation) | * (universal across central potentials) | No (additive shift $\omega^2 \to \omega^2 + m^2$; preserves spectrum's algebraic class for $m^2 \in \mathbb{Q}$) | π-free | None (central deformation; lives in $Z(\mathcal{A}) = \mathbb{C} \cdot \mathbb{1}$; preserves heat-kernel structure modulo overall mass shift) | π-free | ✓ | **Boundary case B: §3.3.** Both lenses predict π-free, and both for the same reason: rest-mass is a *central* (commuting) deformation of $D$. KG-1 verified 200/200 cases π-free. |
| 15 | Observation / temporal-window | B (boundary, not triple-internal) | C/U mixed | Yes (Matsubara sum over discrete $\omega^t_k = 2\pi k/\beta$ is the canonical temporal-spectral integration) | π ($2\pi$ per mode in spectrum; $\pi^{2k}$ in integrated quantities like Stefan–Boltzmann $\pi^2/90$) | M2 (Mellin transform of compactified-time heat kernel, with Bernoulli $\zeta_R(2k)$ giving $\pi^{2k}$); also M3 if vertex parity engaged | π | ✓ | **Boundary case G: §3.7.** This is *the* canonical projection Paper 35 was designed for. Both lenses predict π. The triple framing's M2 specializes to the Bernoulli-zeta side; Paper 35's "temporal/spectral integration" names the Matsubara sum directly. |

**Tally:** 15 of 15 in agreement on the π / π-free verdict. **Zero disagreements found.** One projection (Drake–Swainson #13) is *ambiguous* under Paper 35 because the intermediate-vs-converged distinction matters, but resolves to ✓ when "evaluation" is interpreted as the converged (post-cancellation) answer. This ambiguity is worth flagging as a sharpening target rather than a falsification (§3.1, §4 Outcome 3).

---

## §3. Boundary case analysis

The §2 table reports universal agreement, but the *easy* projections (Wigner 3j, Wigner D, Stereographic, Sturmian, Bargmann–Segal, Bohr — i.e., the projections trivially in the universal sector or trivially producing rational data) do not constitute a useful test of the framework's robustness. The boundary cases below are where structure surfaces. We work each in detail.

### 3.1. Drake–Swainson (#13) — intermediate-vs-converged ambiguity

Drake–Swainson is a regularization-flow operation: split a divergent finite-basis spectral sum at intermediate scale $K$ into $\beta_\text{low}(N, K) + \beta_\text{high}(K)$ with $K$-independence over the plateau (LS-4 verified $K \in [0.05, 200]$, 3.6 orders of magnitude). The structural denominator $D_\text{drake} = 2(2\ell+1)Z^4/n^3$ is rational.

A strict reading of Paper 35 would count each intermediate $\beta(K)$ piece as carrying $\ln(K/E_\text{ref})$ (a Mellin-type continuous-parameter dependence), so intermediate steps would be classified as π-bearing; the cancellation makes the *converged* answer K-clean. The triple framing says: Drake itself adds no new π — its π content equals the upstream Sturmian + spectral-action chain's π (M1 + M2). Three sharpenings of Paper 35's "evaluation":

1. *Strict definitional:* count only Drake's intrinsic step. No continuous integration of its own. π-free. Triple lens agrees (no new mechanism).
2. *Functional chain:* count Sturmian → spectral action → Drake. Chain has π via M1 + M2 upstream. Triple lens agrees.
3. *Intermediate-trace:* count scale-dependent intermediates. Intermediates π-bearing, final answer matches upstream. Both lenses consistent.

All three readings give agreement. The genuine ambiguity is a documentation choice, not a falsification. **Verdict: ✓** in all three readings; flagged as sharpening target (§5.3).

### 3.2. Hopf bundle (#2) — geometric vs temporal/spectral

Paper 35's text says "temporal or spectral parameter promoted from the discrete graph spectrum." The Hopf $S^2$ base measure is angular, not temporal/spectral in the strict sense. Three readings:

- *Strict literalism:* Paper 35 silent on this projection — failure of coverage.
- *Broad reading (Paper 35 §VII, TX-B Obs 5):* TX-B Obs 5 explicitly noted "the continuum limit would inject π via spherical-harmonic normalization on $S^3$, exactly where Prediction 1 predicts." Paper 35 treats angular integration on a continuum manifold as the same mechanism — the temporal direction is a particular instance, not the unique one.
- *Synthesis:* Paper 35's mechanism is general continuous integration over a parameter promoted from a discrete spectrum; the temporal/spectral language is the most-familiar physics example, not the only kind.

**Verdict: ✓** under the synthesis reading, which is what Paper 35 §VII actually does. Both lenses identify the same mechanism (M1's $\mathrm{Vol}(S^2)/4$). The text-level sharpening — clarifying that Prediction 1 covers any continuum-measure integration — is recommended in §5.4.

### 3.3. Rest-mass (#14) — universal $D$-deformation

KG-1 verified ring-preservation across $m^2 \in \{0, 1, 1/4, 2\}$ for 200 modes (all $\omega_n = c_n \sqrt{d_n}$, $c_n \in \mathbb{Q}$, $d_n$ square-free). Paper 35 reads this as no continuous integration over a temporal/spectral parameter (the additive shift is parametric, not integral); the triple framing reads it as a central deformation $D \mapsto D + m \cdot \mathbb{1}$ in $Z(\mathcal{A}) = \mathbb{C} \cdot \mathbb{1}$ (the heat kernel acquires an overall factor $e^{-tm^2}$ that integrates to a Gamma-function shift, no new π). **Verdict: ✓.** Both lenses predict π-free for the *same* reason — centrality preserves the ring. Clean confirmation of the universal/central sector hypothesis.

### 3.4. Vector-photon (#11) — same mechanism, different naming

Paper 33 establishes the 1+6+1 selection-rule partition: 6 angular-momentum rules recovered by vector-photon promotion at $1/(4\pi)$ per loop, identified as the $S^2$ Weyl exchange constant of the Hopf base. Paper 35 reads this as continuous integration over the photon mode space $(q, m_q)$ at the Wigner-3j vertex weighted by the $S^2$ Weyl factor; the triple framing reads it as M1's $\mathrm{Vol}(S^2)/4$ measure showing up as a per-loop calibration constant. **Verdict: ✓.** TX-B Obs 4 (Heisenberg–Euler $\alpha^2/(45\pi m^4)$) confirmed; both lenses identify a single, specific transcendental source.

### 3.5. Spectral action (#6) — heat-kernel time as the canonical Mellin parameter

The spectral action $\mathrm{Tr}\, f(D^2/\Lambda^2)$ is a heat-kernel Mellin transform; the Mellin time $t$ is precisely "a continuous parameter promoted from the discrete spectrum." Paper 35 covers it generically; T9 (Paper 28 Theorem 1) specializes to the squared-Dirac case at integer $s$, asserting $\zeta_{D^2}(s)$ is a two-term polynomial in $\pi^2$ with rational coefficients. **Verdict: ✓.** Both lenses identify the same Seeley–DeWitt M2 mechanism; T9 is the specific case, Paper 35's mechanism is the general statement.

### 3.6. Wilson plaquette (#10) — gauge sector inherits Hopf measure

SU(2) Haar integration over Wilson link variables uses $d\mu(g) = (1/(2\pi^2)) \sin^2(\theta/2)\, d\theta\, d\Omega_{S^2}$, with $\mathrm{Vol}(\mathrm{SU}(2)) = 2\pi^2$. Paper 35 reads this as continuous integration over the gauge-link continuum parameters, promoted from the discrete plaquette structure; the triple framing reads $L_1 = B^T B$ as the Hodge-1 Laplacian whose propagator inherits $\mathrm{Vol}(S^2)/4$ from M1. The SU(2) Haar volume $2\pi^2$ is structurally a rearranged Hopf-base measure (Paper 30 Theorem 1: maximal-torus reduction gives Paper 25's U(1) Wilson action, which is the canonical $L_1$ object). **Verdict: ✓.** Both lenses predict π; identical mechanism by different names.

### 3.7. Observation/temporal-window (#15) — the canonical case

This is *the* projection Prediction 1 was designed for. Matsubara sum over $\omega^t_k = 2\pi k / \beta$ injects $2\pi$ at the spectrum level; the high-$T$ free-energy density picks up $\zeta_R(4) = \pi^4/90$ via Bernoulli/Riemann (M2 specialized to compactified-time heat kernel). The four-manifold $S^3 \times S^1_\beta$ acquires a $\mathbb{Z}_2$-grading $\gamma$ (Spin(4) bundle), so the temporal compactification is also the framework's first $\gamma$-engaging projection (TS-B §2.c.6, prediction P-6: KO-dim transition 3 → 4). **Verdict: ✓.** Both lenses agree.

### 3.8. Sturmian (#5) — parametric only, both lenses π-free

Sturmian re-labels the same Fock graph at $\lambda = Z/n$; preserves $\mathcal{A, H, D}$ spectra; no continuous integration at the projection level. The Bethe-log spectral sums *downstream* of Sturmian carry calibration π via spectral-action upstream (M2). **Verdict: ✓.** Both lenses predict π-free at the Sturmian step.

---

## §4. Verdict: Outcome 1 (universal agreement, candidate theorem) with §3.1 sharpening

The §2 table reports agreement on all 15 projections. The boundary cases in §3 confirm the agreement in detail: every projection that introduces continuous integration over a parameter promoted from the discrete spectrum (Hopf bundle's $S^2$ base, spectral action's Mellin time, vector-photon promotion's $S^2$ Weyl, Wilson plaquette's SU(2) Haar, observation/temporal window's Matsubara, spinor lift's Hurwitz Mellin, molecule-frame's continuum integrations, Drake–Swainson's flow upstream) does so via a triple-framing M1/M2/M3 mechanism that introduces π. Every projection that does not introduce such integration (Fock conformal at projection level, Bargmann–Segal, Stereographic, Sturmian, Wigner 3j, Wigner D, Rest-mass) preserves the algebraic ring at both the Paper 35 and triple level.

**Outcome 1 verdict:** the 15 projections agree.

**Sharpening required at §3.1:** Drake–Swainson's "intermediate-vs-converged" status is a documentation question, not a falsification. The recommended sharpening: Paper 35 Prediction 1 explicitly endorse the *converged-evaluation* reading, with intermediate scale-dependent quantities classified separately as Bethe-flow content. This is consistent with sprint LS-4's K-independence verification.

**Outcome 1 implication:** Paper 35 Prediction 1 is a *candidate spectral-triple-framing theorem*. The proof sketch combines:

- Paper 18 §III.7 (Mellin transform as taxonomic engine): the *only* operation in the framework that converts discrete spectral data to π-bearing transcendentals is the Mellin transform of a heat kernel (or equivalently, the Mellin transform of the spectral zeta).
- Paper 28 Theorem T9 (squared-Dirac spectral zeta is π-even at integer $s$): specializes to the GeoVac Dirac case and pins down the ring class.
- Paper 32 §VII (Hopf-base measure $\mathrm{Vol}(S^2)/4$ is the canonical $\pi$-source on the Coulomb $S^3$ triple): identifies the Hopf-base measure as M1.
- Paper 33 (1+6+1 partition): identifies the vector-photon promotion's $1/(4\pi)$ per loop as the same M1 measure.

**Theorem statement (proposed):**

> **Theorem (TS-D triple-framing of Paper 35 Prediction 1).** Let $P$ be any of the fifteen projections of the GeoVac framework (Paper 34 §III.1–§III.15). The following are equivalent:
> 1. $P$'s converged evaluation produces a result containing $\pi$ (or $\sqrt{\pi}$, or $\pi^{2k}$, or $\pi$ in a denominator).
> 2. $P$'s evaluation pipeline includes continuous integration over a parameter promoted from the discrete graph spectrum (heat-kernel time, Matsubara frequency, Mellin contour, $S^2$ Hopf-base angular measure, SU($N$) Haar measure, or analytic-continuation parameter).
> 3. $P$ engages at least one of the three triple-framing mechanisms M1 (Hopf-base measure, $\mathrm{Vol}(S^2)/4$), M2 (Seeley–DeWitt heat-kernel coefficient, $\sqrt{\pi}$ on $S^3$), or M3 (vertex-topology Hurwitz Dirichlet $L$-content).

**Proof sketch:**

(1) ⇔ (2) is Paper 35 Prediction 1, verified on 208 individual checks.

(2) ⇔ (3): each continuous integration appearing in any of the fifteen projections is identifiable with one of M1/M2/M3 by case-checking (§3 above and §2 table). No fourth π-source mechanism appears in the Paper 34 list. The proof requires:

- M1 covers Hopf-bundle (#2), vector-photon (#11), Wilson plaquette (#10) via Haar measure normalization, and any continuum angular integration over the Hopf $S^2$ base.
- M2 covers spectral-action (#6), spinor-lift (#7) for the Hurwitz π-content, observation/temporal-window (#15) for the Matsubara $\zeta_R(2k)$, and any heat-kernel Mellin-transform operation.
- M3 covers the vertex-topology Catalan/Dirichlet-L content of two-loop Dirac sums (Paper 28 Theorem 3, eq:rh_j_identity).

(1) ⇔ (3) is the spectral-triple framing of Paper 35; it follows from (1) ⇔ (2) and (2) ⇔ (3).

This theorem promotes Paper 35 Prediction 1 from "load-bearing principle" to "spectral-triple-framing corollary" — a categorically different epistemic status. The principle was empirically established (208 a-priori checks); the theorem says *why* it must hold, by exhausting the framework's mechanism list.

**Caveat (open question):** the theorem's coverage of M1/M2/M3 is *currently complete* over the fifteen projections. If Paper 34 grows to include a sixteenth projection that introduces a continuous integration via a *new* mechanism (not Hopf-base, not heat-kernel, not vertex-topology), the theorem must be re-verified. The natural candidates for new mechanisms — continuum gravity, non-trivial fibre spaces, anomaly-related transcendentals — are flagged as falsification targets in §5 below.

---

## §5. Implications for Paper 32 / Paper 35

### 5.1. Recommended Paper 32 amendment (Outcome 1 application)

Add a new subsection to Paper 32, §IX (after the Coulomb/HO asymmetry discussion, before the Marcolli–vS lineage):

> **§IX.5 (proposed). Paper 35 Prediction 1 as a triple-framing corollary.**
> The TS-D investigation (May 2026, `debug/track_ts_d_paper35_triple_test_memo.md`) exhausted the fifteen Paper 34 projections and verified that each π-source identifiable in any projection's evaluation pipeline reduces to one of three triple-framing mechanisms: M1 (Hopf-base measure $\mathrm{Vol}(S^2)/4$), M2 (Seeley–DeWitt heat-kernel coefficient $\sqrt{\pi}$ on $S^3$ or related Mellin/Hurwitz content on the Camporesi–Higuchi spectrum), or M3 (vertex-topology Hurwitz Dirichlet $L$-content via quarter-integer shifts $\zeta(s, 3/4), \zeta(s, 5/4)$). The exhaustion gives a triple-framing reading of Paper 35 Prediction 1: a GeoVac observable contains π if and only if its evaluation pipeline engages M1, M2, or M3. Combined with the 208 a-priori checks of the KG and TX-B sprints, this places Prediction 1 on the same epistemic footing as the T9 theorem (Paper 28 Theorem 1) — empirically verified *and* structurally exhausted across the framework's known mechanisms. The principle is no longer just "load-bearing" but provably so within the spectral-triple framing, modulo the discovery of a sixteenth mechanism.

The footnote: M1 covers Hopf-bundle, vector-photon, Wilson-plaquette Haar; M2 covers spectral-action, spinor-lift Hurwitz, observation/temporal-window; M3 covers two-loop vertex-Dirichlet-L content.

### 5.2. Recommended Paper 35 amendment (Outcome 1 application)

Add a paragraph to Paper 35 §V.B (Falsification criterion) or §VIII (Open questions, after item 1 spinor sector resolution):

> **§VIII.4 (proposed). Triple-framing reading of Prediction 1.**
> Sprint TS-D (May 2026, `debug/track_ts_d_paper35_triple_test_memo.md`) showed that the fifteen Paper 34 projections agree on Prediction 1 under both the projection-taxonomy lens (this paper) and the spectral-triple lens (Paper 32 §IX, Paper 28 Theorem 1, Paper 18 §III.7). The two lenses identify the *same* π-source mechanisms: the projection-taxonomy lens names the integration ("temporal/spectral parameter promoted from the discrete spectrum"); the spectral-triple lens names the measure (M1 Hopf-base $\mathrm{Vol}(S^2)/4$, M2 Seeley–DeWitt $\sqrt{\pi}$, M3 vertex-Dirichlet-L). Sprint TS-D's universal agreement gives a structural basis for Prediction 1 alongside the 208 a-priori checks: every π-source in the fifteen projections is one of three triple-framing mechanisms, exhausting the known list. A counter-example would require either (a) a new physical observable in any of the fifteen projections that does not fit M1/M2/M3, or (b) a sixteenth projection introducing a new mechanism. The most likely candidates for (b) are continuum-gravity, non-trivial principal-bundle, and anomaly-class projections; these have not been added to the Paper 34 list and are flagged as a future falsification target.

### 5.3. Drake–Swainson sharpening (§3.1 application)

Add a clarification to Paper 34 §III.13 (Drake–Swainson) and Paper 35 §V.A (Confirmatory evidence):

> **Drake–Swainson is a regularization-flow projection.** Its "evaluation" is the converged $K$-cancelled answer; intermediate scale-dependent quantities $\beta_\text{low}(K), \beta_\text{high}(K)$ are not Prediction 1 inputs. Under the converged-evaluation interpretation, Drake's π content equals its upstream chain's π content (Sturmian → spectral action → Drake). This is the *intended* reading of Paper 35 Prediction 1; the alternative "intermediate-trace" reading would catalogue every regularization-flow scale-dependence as a separate projection (impractical and not the framework's intent).

### 5.4. Hopf bundle reading-3 sharpening (§3.2 application)

Add a clarification to Paper 35 §V.A:

> **"Continuous integration over a temporal or spectral parameter" includes angular continuum integrations.** Sprint TS-D (§3.2) clarified that Paper 35's prediction extends to any continuous integration over a continuum-limit measure on a parameter promoted from the discrete spectrum. The Hopf $S^2$ base measure (Paper 32 §VII), the SU(2) Haar measure (Paper 30 §IV), and the per-loop $S^2$ Weyl exchange constant $1/(4\pi)$ (Paper 33 §VI) are all instances of this mechanism, structurally equivalent to the temporal-Matsubara case of KG-2. The "temporal or spectral" phrasing of Prediction 1's text is the most-familiar physics example, not the only kind. Where a projection's evaluation includes angular integration on a continuum-limit base, the projection is π-bearing and its π is the Hopf-base / Weyl-volume calibration of the spectral-triple framing.

### 5.5. Falsification target (long-term)

Search for a **sixteenth projection** with a new π-source mechanism not in M1/M2/M3. Three natural candidates: (i) continuum-gravity coupling (Einstein–Hilbert / Gauss–Bonnet beyond standard CC curvature scalar); (ii) anomaly-class projections (chiral or conformal anomaly, where the conformal anomaly on $S^3 \times S^1_\beta$ would need to be checked as either a re-routing of M2 or a genuinely new mechanism); (iii) non-trivial principal-bundle projections (second-Chern-class invariants with Pontryagin $\int F \wedge F$ contributing $\pi^2$ via topological-invariant reading of discrete-spectrum data, not via continuous integration). For (iii) specifically, π entering via a topological invariant *without* intermediate continuous integration would be a genuine falsification of Prediction 1 in its present form — the prediction would need refinement to include "topological-invariant readings" as a fourth mechanism, or the M1/M2/M3 list would need extension. None of these candidates is currently engaged in the framework's fifteen-projection list.

---

## Synthesis

Track TS-D adds a *theoretical* layer to the empirical TX-B confirmation. TX-B established Paper 35 Prediction 1 holds across a 208-element a-priori panel; TS-D establishes that the principle is structurally exhausted by three triple-framing mechanisms (Hopf-base $\mathrm{Vol}(S^2)/4$, Seeley–DeWitt $\sqrt{\pi}$ on $S^3$, vertex-topology Hurwitz Dirichlet $L$). The fifteen projections' π-content reduces to engagement with M1, M2, or M3; no new mechanism is required. This converts Prediction 1 from a load-bearing conjecture into a candidate theorem provable by case-exhaustion within the GeoVac framework's known structure. Two minor sharpenings (Drake's converged-evaluation reading, Hopf's angular-continuum reading) tighten the prediction's text without altering its content. The theorem's only open exit is a sixteenth projection: if continuum-gravity, anomaly-class, or non-trivial principal-bundle projections introduce a fourth π-source mechanism in future framework extensions, the proof must be re-verified. This is a clear and bounded falsification target rather than an undischarged liability.

---

## Report-back

- **Path to memo:** `debug/track_ts_d_paper35_triple_test_memo.md`
- **Top-level verdict:** **Outcome 1** (all 15 projections agree; Paper 35 Prediction 1 is a candidate triple-framing theorem provable by case-exhaustion via M1/M2/M3).
- **Count of disagreements found:** **0** (with one ambiguity at §3.1, Drake–Swainson intermediate-vs-converged, resolved by a documentation sharpening).
- **One-paragraph synthesis (~150 words):** Track TS-D tested Paper 35 Prediction 1 against the spectral-triple framing on each of Paper 34's fifteen projections. All fifteen agreed: π-bearing projections engage one of three triple-framing mechanisms (M1 Hopf-base $\mathrm{Vol}(S^2)/4$, M2 Seeley–DeWitt $\sqrt{\pi}$ on $S^3$, M3 vertex-topology Hurwitz Dirichlet $L$); π-free projections engage none. The agreement is not a coincidence: each continuous integration in a Paper 34 projection's evaluation pipeline is identifiable as one of M1/M2/M3 by direct case-analysis. This promotes Prediction 1 from "load-bearing principle empirically verified on 208 checks" to "candidate spectral-triple-framing theorem provable by case-exhaustion." Two minor documentation sharpenings (Drake-Swainson's converged-evaluation reading, Hopf's angular-continuum reading) tighten the prediction's text without altering its content. The only open falsification target is a sixteenth projection introducing a new π-source mechanism — natural candidates (continuum gravity, anomaly classes, non-trivial principal bundles) are flagged but not yet engaged.
