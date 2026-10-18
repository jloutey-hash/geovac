# The Forced/Free Seam

*A catalogue of what the GeoVac framework fixes (relative to its stated axioms and organizing readings) and what it consumes as external input.*

> **Re-synced 2026-10-02 to Paper 57 and Paper 32 §VIII** (/qa trunk DELTA #25). **Widened 2026-10-15 (PI-directed):** status **W (Wall)**, sub-domain **J (released axes)**, rows A6/A7/I4, and the period column's `weight-1` / `elliptic` / `K3` classes; the door-indexed view is `docs/atlas.md`, gated by `tests/test_atlas_register.py` (which also checks the summary counts below against the rows). Paper 57 is authoritative; where a row and the paper differ, the paper wins. "F" for the gauge rows means *forced relative to the Bertrand × Hopf-tower organizing reading plus the adopted Upgrade B axiom*, not derived.

Started 2026-06-08. Supersedes `docs/archive/forcing_catalogue.md` (archived 2026-10-15) for the forced-vs-free question; the prior catalogue's forward-run / probe structure remains the canonical record for *what does the same rigidity force that we haven't looked at*.

---

## Purpose

GeoVac is bit-exact on a wide class of structural facts (Paper 0 packing, the S³ label set and SO(4) Casimir spectrum, dim M(D_F) given the admitted A_F and the imported CCM representation, the master Mellin engine census, …; the gauge group relative to the adopted Upgrade B axiom) and consumes specific numerical content as external input (Yukawa values, generation count, calibration cutoff f, atomic-physics Z_eff fits, …). The boundary between these two regimes — the **forced/free seam** — is itself a structural property of the framework, repeatedly named (`memory/multi_focal_wall_pattern.md`, `memory/external_input_three_class_partition.md`, `memory/geovac_structural_skeleton_scope_pattern.md`, Paper 32 §VIII Forced-Count Theorem, Paper 18 §IV.6 inner-factor input-data tier) but never consolidated.

This catalogue consolidates it. Companion paper: `papers/group3_foundations/paper_57_forced_free_seam.tex`.

## Status conventions

- **F (Forced)** — fixed by the Paper 0 packing axiom + standard spectral-triple machinery, relative to any adopted axiom or organizing reading named in the row. Witness: a theorem or bit-exact verification.
- **A (Admitted)** — the framework structurally allows this but does not autonomously select it among multiple structurally-equivalent options. Witness: existence proof + non-selection theorem or negative search.
- **Obs (Observation)** — a matched or coincident value with no derivation bridge (e.g. κ = −1/16); neither forced nor consumed as calibration input.
- **C (Calibration)** — the framework consumes this as external input. Witness: a non-selection theorem (preferred) or a documented negative search (acceptable).
- **W (Wall)** — a proven obstruction with its mechanism pinned (the third verdict of the CLAUDE.md §1 mission statement): neither forced nor consumed. Witness: a theorem or a bit-exact/high-precision demonstration of the obstruction, plus the named mechanism. Added 2026-10-15 with sub-domain J; the wall's working record stays in `docs/walls/register.md`, this row only indexes it.

## Discriminator columns

For each entry we record three properties to be tested in the paper §5 candidate-principle probe:

- **MF** — multi-focal depth. `1` if the observable is the output of a single Paper-34 projection; `>1` if its definition requires composing multiple projections.
- **P** — period class. `M1` / `M2` / `M3` / `outside` / `none` per Paper 18 §III.7 master Mellin engine; `meta` for the engine itself; `N/A` for non-period observables. Widened 2026-10-15 for the released-axis rows (J): `weight-1` (the two-centre seed set {e^a E₁(a), ln, γ_E}, π-free, Paper 18 `sec:level2_seed_set`), `elliptic` (the three-centre genus-one fibre, Paper 59 `sec:genus`), `K3` (the three-electron chain's two-momentum fibre, `docs/predictions/2026-10-13_three_electron_k3.md`). These sit outside the Mellin engine's M1/M2/M3 and outside the polylogarithm/MZV weight tower (Paper 18 `sec:level2_genus`); `outside` is kept for the calibration rows whose class is not identified.
- **Dim** — dimensional character. `dimensionless` (ratio) / `dimensionful` (scale).

---

## A. Foundations / skeleton

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| A1 | Label set $(n,l,m)$ with the diagonal SO(4) Casimir $n^2-1$ ($-\Delta_{S^3}$ on the labels) | **F** | Packing axiom; Paper 0 + Paper 7 (18 symbolic checks of the continuum chain). The graph edge Laplacian $D-A$ does NOT have this spectrum and does not converge to $\Delta_{S^3}$ (WALL, 2026-09-25) | Skeleton | 1 | none | dimensionless |
| A2 | $S^3$ as forced outer manifold | **F** | Bertrand + SO(4) closure; Paper 7 | Skeleton | 1 | none | dimensionless |
| A3 | $\kappa = -1/16$ matched kinetic scale | **Obs** | Matched to the graph's band edge; its coincidence with the Fock value $1/\Omega^4(0)$ is an Observation, not a derivation | Skeleton | 1 | none | dimensionless |
| A4 | Shell degeneracy $g_n = n^2$ | **F** | Representation theory; exact integer | Skeleton | 1 | none | dimensionless |
| A5 | Compactness ⇒ discreteness (Peter–Weyl; the converse fails, e.g. the HO of Paper 24 — discrete spectrum ⇔ compact resolvent) | **F** | Peter–Weyl; `memory/geovac_structural_skeleton_scope_pattern.md` | Skeleton | 1 | N/A | N/A |
| A6 | Edge Laplacian $D-A$ of the label lattice as a route to $\Delta_{S^3}$ | **W** | Does not converge: $l$-preserving, $n_{\max}$ disconnected blocks, $\dim\ker L = n_{\max}$ at every cutoff, zero cross-block heat kernel; spectrum algebraic ($2-2\cos(j\pi/N)$ sums). Paper 0 `sec:graph`, Paper 7 `sec:attribution`; `tests/test_graph_s3_nonconvergence.py` (2026-09-25) | Skeleton | 1 | none | dimensionless |
| A7 | Bargmann–Segal lattice of the 3D harmonic oscillator: π-free, bit-exact rational at every $N_{\max}$ (56 nodes, 165 edges, zero irrationals at $N_{\max}=5$) | **F** | Paper 24 `thm:pi-free`; the only π is π⁻³ in the continuous Bargmann measure, never invoked by the graph; the Maz'ya–Shubin counterexample to the converse of A5 (stated in Paper 18 `sec:pattern`); `tests/test_bargmann_graph.py` | Skeleton | 1 | none | dimensionless |

## B. Gauge structure

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| B1 | Gauge group $U(1) \times SU(2) \times SU(3)$ | **F** (relative to adopted Upgrade B) | Bertrand × Hopf-tower organizing reading truncated at $n \le 3$, plus the adopted Upgrade B axiom; Paper 32 §VIII.B | Gauge | 1 | none | dimensionless |
| B2 | Inner algebra factor count (= 3) | **F** (relative to the organizing reading) | Hurwitz + Hopf-tower; Door 4b/4d | Inner | 1 | none | dimensionless |
| B3 | $\mathbb{C}$ factor at $n=1$ | **F** | Door 4d (unique division algebra of $S^1$) | Inner | 1 | none | dimensionless |
| B4 | $M_3(\mathbb{C})$ factor at $n=3$ | **F** | Door 4d (Hurwitz fallback at $S^5$) | Inner | 1 | none | dimensionless |
| B5 | $\mathbb{H}$ vs $M_2(\mathbb{C})$ at $n=2$ | **A** | Door 4c NEGATIVE on $J$-sign-table; ℍ admitted, not forced; Door 4e Upgrade B closes if adopted | Inner | 1 | none | dimensionless |
| B6 | $\dim \mathcal{M}(D_F) = 32$ at one generation (matter-sector projection rank 16) | **F** (given the admitted $\mathcal{A}_F$ and imported CCM representation) | Forced-Count Theorem, Paper 32 §VIII (`thm:forced_count`); the retired 128 / 260 were a degenerate-sample artefact (corrected 2026-09-03) | Inner | 1 | none | dimensionless |
| B7 | Higgs admission (Mexican-hat structurally allowed) | **F** | Sprint H1 POSITIVE-THIN; Paper 32 §VIII.C | Inner | 1 | none | dimensionless |
| B8 | Lower-bound forcing on gauge content (saturated, not subset) | **C** | Bertrand × Hopf gives upper bound only; structural-skeleton scope pattern | Gauge | 1 | none | dimensionless |

## C. Spectral / period content

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| C1 | Master Mellin engine $\mathcal{M}[\mathrm{Tr}(D^k \cdot e^{-tD^2})]$, $k \in \{0,1,2\}$ (at $k=1$ read through the absolute sum and its vertex-parity grading; the plain trace vanishes on the symmetric spectrum) | **F** (CONDITIONAL) | Case-exhaustion theorem, [CONDITIONAL], first fifteen Paper 34 projections; converse withdrawn as not established; Paper 32 §VIII | Periods | 1 | meta | dimensionless |
| C2 | M1 prefactor = $\mathrm{Vol}(S^2)/4 = \pi$ | **F** | Unit-sphere volume slot (the earlier "Hopf-base" label is retired) | Periods | 1 | M1 | dimensionless |
| C3 | Seeley–DeWitt (Dirac, 4-spinor) $a_0 = \sqrt\pi$, $a_1 = -\sqrt\pi/2$, $a_{k\ge2} = 0$ on $S^3$ (corrected 2026-10-02 from $a_0 = a_1 = \sqrt\pi$, $a_2 = \sqrt\pi/8$, a production sign error) | **F** | Closed form; Paper 28 (two-term exactness) | Periods | 1 | M2 | dimensionless |
| C4 | $\zeta_{D^2}(s) = 2^{2s-1}[\lambda(2s-2) - \lambda(2s)]$ | **F** | T9 theorem; Paper 28 | Periods | 1 | M2 | dimensionless |
| C5 | $\chi_{-4}$ identity $D_{\rm even} - D_{\rm odd} = 2^{s-1}(\beta(s) - \beta(s-2))$ | **F** | Paper 28 §IV; vertex-parity Hurwitz | Periods | 1 | M3 | dimensionless |
| C6 | $M_3 \subset \mathcal{MT}(\mathbb{Z}[i, 1/2])$ at level $\le 4$ | **F** | Sprint Mixed-Tate; Paper 55 | Periods | 1 | M3 | dimensionless |
| C7 | Fejér-moment constant $4/\pi$ of the state-space GH bound (dual-Coxeter; $2/\pi$ unit $S^3$; a kernel constant, not the convergence rate) | **F** | Paper 38, commutative factor on the true window, [INTERNAL THEOREM + MEASURED] (not propinquity); Paper 40 reports it universal over compact simple groups (derived at rank 1; general-G convergence conditional on the spin-window decomposition); its M1 placement is a reading | OA | 1 | M1 (reading) | dimensionless |
| C8 | Pythagorean orthogonality $\langle H_{\rm local}, D_W^L\rangle_{\rm HS} = 0$ | **F** | Bit-exact; Paper 42 §11.7 Thm 11.14 (former Paper 43 §10.2, archived 2026-09-28); β-independent chirality-parity identity — its residual's 1/π² is the posited β squared, not an M1 output | OA | 1 | M1 (class only; re-tiered 2026-09-28) | dimensionless |
| C9 | $S_{\min} = 8\pi^2\ln 2 - \tfrac{2}{3}\pi^4\ln 2 - 3\pi^2\zeta(3) + \tfrac{1}{2}\pi^4\zeta(3) - \tfrac{5}{2}\pi^2\zeta(5) + \tfrac{\pi^6}{4} - \tfrac{3}{2}\pi^4 - \tfrac{\pi^8}{96}$ | **F** | Paper 28 (identified 2026-06-11; prior "irreducible" = basis-coverage artifact) | Periods | >1 | M3 | dimensionless |

## D. Gravity

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| D1 | Spectral action two-term-exact on $S^3$ | **F** | Bernoulli identity $\zeta_{\rm unit}(-k) = 0$; Paper 51 | Gravity | 1 | M2 | dimensionless |
| D2 | $S_{\rm BH} = A\Lambda^2/(12\pi)$; $G_N = 3\pi/\Lambda^2$ | **F** | Closed form; Paper 51 §G4 | Gravity | 1 | M1 | dimensionful |
| D3 | Cone coefficient $-1/12$ (Sommerfeld–Cheeger) | **F** | Bit-exact discrete extraction; Paper 51 G4-4c | Gravity | 1 | M2 | dimensionless |
| D4 | Replica derivative $d\Delta_K/d\alpha\rvert_{\alpha=1} = +1/6$ | **F** | Closed form; Paper 51 G4-4f | Gravity | 1 | M2 | dimensionless |
| D5 | Cutoff function $f$ moments $\varphi(0), \varphi(1), \varphi(2)$ | **C** | No autonomous selection; CC fine-tuning $\varphi(2)/\varphi(1)^2 \approx 10^{-124}$ required | Gravity | >1 | outside | dimensionful |
| D6 | Cosmological constant value $\Lambda_{cc}$ | **C** | Inherits D5 | Gravity | >1 | outside | dimensionful |
| D7 | Relations among $G_N, \Lambda_{cc}, S_{\rm BH}$ (action-$G \equiv$ entropy-$G$) | **F** | Wald forces; two-term-exactness $\Rightarrow$ pure Einstein | Gravity | 1 | M1 | dimensionless |

## E. QED and the $\alpha$ coincidence

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| E1 | Self-energy structural zero $\Sigma(n_{\rm ext}=0) = 0$ | **F** | Proven; Paper 28 | QED | 1 | none | dimensionless |
| E2 | Selection rules 8/8 (1+6+1 partition) | **F** | Paper 33; vector-photon promotion at $1/(4\pi)$ per loop | QED | 1 | M1 | dimensionless |
| E3 | $F_2 = 5\sqrt{2}/3$ graph-native ($\pi$-free) | **F** | Paper 33 graph-native QED | QED | 1 | none | dimensionless |
| E4 | $F_2/(\alpha/2\pi) = 1.084$ Parker–Toms ($R/12$ curvature) | **F** | Bit-exact; 0.5% match | QED | 1 | M2 | dimensionless |
| E5 | $B = 42$ (Casimir), $F = \pi^2/6$ (Fock Dirichlet), $\Delta = 1/40$ ($g_3^{\rm Dirac}$) | **F** | Three distinct spectral homes (distinct type proven); no common generator found [OBSERVATION]; 12 mechanisms eliminated | $\alpha$ | 1 | M1/M2/M2 | dimensionless |
| E6 | Combination rule $K = \pi(B + F - \Delta)$ for $\alpha^{-1}$ | **C** (an Observation, consumed as calibration) | Numerical coincidence: $K$ matches $\alpha^{-1}$ to $4.8 \times 10^{-7}$ (its cubic's root to $8.8 \times 10^{-8}$); no derivation; stays an Observation per CLAUDE.md §13.5 | $\alpha$ | >1 | M1+M2 | dimensionless |
| E7 | Two-loop SE counterterms $Z_2, \delta m$ | **C** | LS-8a wall: framework cannot generate autonomous renormalization | QED | >1 | outside | dimensionful |
| E8 | Multi-loop QED corrections beyond LS-7 | **C** | HF-5 wall (same mechanism as E7) | QED | >1 | outside | dimensionful |

## F. Standard Model inner factor (the load-bearing free row)

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| F1 | Yukawa values $y_e, y_\mu, y_\tau, y_u, \ldots$ (9 charged-lepton+quark) | **C** | Sprint H1 non-selection theorem + 162-cell PSLQ clean negative against M1∪M2 at $M \le 10^3$ | Inner | >1 | outside | dimensionless |
| F2 | Generation count $N_{\rm gen} = 3$ | **C** | Direction 2 NO-GO + Read 2 NO-GO (3 algebra factors $\neq$ 3 generations in standard rep) | Inner | 1 | none | dimensionless |
| F3 | Inner KO-dim signature | **C** | Direction 2 NO-GO; packing-unreachable | Inner | 1 | none | dimensionless |
| F4 | Higgs VEV $v$ | **C** | Mass scale; calibration input | Inner | >1 | outside | dimensionful |
| F5 | CKM mixing angles | **C** | Cross-generation; Boyle–Farnsworth gives shape not values | Inner | >1 | outside | dimensionless |
| F6 | PMNS mixing angles | **C** | Same as F5 | Inner | >1 | outside | dimensionless |
| F7 | Neutrino mass values | **C** | Inherits F1 + F4 | Inner | >1 | outside | dimensionful |

## G. Multi-focal-composition walls (the cross-domain unity row)

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| G1 | HF-3 recoil cross-register $V_{eN}(\hat r_e, \hat R_n)$ | **C** | Multi-focal wall pattern instance 1; spatial $\times$ spatial composition | QED | >1 | outside | dimensionful |
| G2 | HF-4 Zemach magnetization density | **C** | Multi-focal wall pattern instance 2 (W1b extended) | QED | >1 | outside | dimensionful |
| G3 | HF-5 multi-loop QED on hyperfine | **C** | Multi-focal wall pattern instance 3 | QED | >1 | outside | dimensionful |
| G4 | LS-8a two-loop SE renormalization (= E7) | **C** | Multi-focal wall pattern instance 4 | QED | >1 | outside | dimensionful |
| G5 | W1e chemistry inner-region overattraction | **C** | Multi-focal wall instance 5; chemistry-side analog of H1 (W1e period-class memo) | Chemistry | >1 | outside | dimensionful |
| G6 | H1 Yukawa non-selection (= F1) | **C** | Multi-focal wall pattern instance 6 | Inner | >1 | outside | dimensionless |

## H. Chemistry / QC scaling (forced architectural)

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| H1 | Angular sparsity (6.06% ERI density at $l_{\max}=3$, global-M_L rule) depends only on $l_{\max}$ | **F** | Paper 22 potential-independent theorem | Chemistry | 1 | none | dimensionless |
| H2 | $N_{\rm Pauli} = 27.90 \times Q$ exactly linear across the main-group composed molecules ($30.03 \times Q$ d-block; fixed basis) | **F** | Paper 14; isostructural invariance | QC | 1 | none | dimensionless |
| H3 | Pauli ratio balanced/composed grows as $O(B^2)$ | **F** | Structural-architectural level (Paper 57; the earlier per-cell ratios predate the 2026-08-29 exact rule and are not quoted) | QC | 1 | none | dimensionless |
| H4 | HO shell closures 2, 8, 20, 40, 70, 112 (spin-orbit-free; NOT the physical magic numbers 2, 8, 20, 28, 50, 82, 126) | **F** | Shell counting on the Bargmann–Segal lattice (shell degeneracy (N+1)(N+2) with spin; cumulative (N+1)(N+2)(N+3)/3); Paper 24 | Nuclear | 1 | none | dimensionless |
| H5 | Fock projection rigidity ($S^3$ unique to $-Z/r$) | **F** | Paper 23 NJ theorem | Nuclear | 1 | none | dimensionless |
| H6 | Clementi–Raimondi $Z_{\rm eff}$ exponents | **C** | External fits; W1e period-class memo classifies as inner-factor input-data tier | Chemistry | 1 | outside | dimensionless |
| H7 | Multi-zeta coefficients in second-row valence basis | **C** | Inherits H6 | Chemistry | 1 | outside | dimensionless |

## I. Foundational calibration (Class 1)

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| I1 | Fine-structure constant $\alpha$ (value) | **C** | E6 combination rule is an Observation; $\alpha$ value external | $\alpha$ | >1 | M1+M2 | dimensionless |
| I2 | Born rule probability rule $p = \lvert\langle a\lvert\psi\rangle\rvert^2$ | **C** | Gleason via Hilbert inheritance; framework does not improve on Gleason; `external_input_three_class_partition.md` Class 1 | Foundations | 1 | N/A | N/A |
| I3 | Higgs direction $\hat n \in S^2$ | **C** | Boyle–Farnsworth input | Inner | 1 | M1 | dimensionless |
| I4 | Temporal observation window β (KMS / Matsubara circle); modes $2\pi k/\beta$ | **C** | The window is observer-supplied; once time is compactified to a circle the $2\pi$ per mode is a volume-slot (M1) reading, Paper 35 `sec:where_pi` (π enters at temporal compactification; WH7 registered, not proven). Open at register level: a compact *interval* is not a compact *circle* (CLAUDE.md §1.7 WH7, 2026-09-06); Paper 35 reads the finite real-time window as the Euclidean circle, an identification it adopts rather than derives (restated 2026-10-16), with both caveats (interval vs circle; Euclidean vs real time) and Kempf 2010 as the published cap; whether a real-time window forces the circle stays open. Falsifiers `tests/test_paper35_kg_panel.py`, `tests/test_wh7_toeplitz_temporal.py` | Foundations | 1 | M1 (reading) | dimensionful |

## J. Released axes: the decompactification doors (added 2026-10-15)

Each row is one axis released from the compact skeleton toward the continuum, with the number class that enters as its price. The ladder these rows trace — exact rationals at one centre for two-electron integrals; weight-one logarithms at two centres; an elliptic period at three; a K3 period for a three-electron chain; a period not built from elliptic ones for four linked densities; an elliptic period again, at the fibre, for a Coulomb ring of three — is stated piecewise in Paper 18 `sec:level2_seed_set` / `sec:level2_genus`, Paper 59 `sec:genus` and `sec:beyond`, and the K3 prediction record; this section and Paper 57's sub-domain J index it together. The door-indexed view with falsifiers is `docs/atlas.md`.

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| J1 | Second nucleus (prolate, diatomic): $l$ released across centres, $m$ survives as additive $U(1)$ conservation; two-centre ERI classes (AA\|BB) elementary and π-free, (AA\|AB)/(AB\|BB) {E₁, ln}, (AB\|AB) {E₁, ln, γ_E}: weight one, π-free, no dilogarithm | **F** (the $m$ rule is a theorem; the seed classes are counted, the elementary block decided) | Paper 58 `thm:abelian_residue`, `sec:mechanism`, `sec:census`; Paper 12 `sec:Bl_closed`; Paper 18 `sec:level2_seed_set`. Prior art: Curie 1894, Bethe 1929 (the $m$ rule); Neumann 1878, Roothaan 1951, Ruedenberg 1951 (the seeds). `tests/test_paper58_abelian_residue.py`, `tests/test_paper58_eri_census.py`, `tests/test_two_center_eri_aabb.py`, `tests/test_paper12_general_m_neumann.py` | Chemistry | 1 | weight-1 | dimensionless |
| J2 | Decompactification front: the coupling that $l$ forbade turns on continuously in $R$, set by the decay length $n/Z$; the label stays Boolean | **Obs** | Paper 58 `sec:continuous` (front $R^*/\ell_n$ and exponent MEASURED; rate bound $0 < dA/dc^2 < 1$ via DLMF spheroidal, SYMBOLIC + MEASURED). Prior art: Herring 1962, Herring–Flicker 1964. `tests/test_paper58_decompactification_front.py`, `tests/test_paper58_dlmf_spheroidal.py`, `tests/test_paper58_coherence_front.py` | Chemistry | 1 | none (seed $e^aE_1(a)$, $a = R$) | dimensionful |
| J3 | Third nucleus, shared centre (XY\|XZ): the fibre is the genus-one curve $y^2 = (c_1k^2+1)(c_2k^2+1)$, period $K(m)$; CM fibre at $\rho = 1/2$ ($\tau = i$): $K(1/2) = \Gamma(1/4)^2/(4\sqrt\pi)$; the radial core $N(D)$ is a period of an irreducible connection, so it has no closed form through any factorization, and the elliptic-dilogarithm route is obstructed; whether the integrated three-centre integral has a finite closed form is open, and accurate evaluation is unaffected | **W** (closed forms of the radial core only) | Paper 59 `sec:genus`, `sec:obstruction`, `sec:scope`; Paper 61 `sec:modular`; Paper 18 `sec:level2_genus`; walls register GEOMETRY-2FOCI. Prior art: Shibuya–Wulfman 1965, Barnett–Coulson 1951, Chowla–Selberg 1967; the sunrise literature. `tests/test_routeC_momentum.py`, `tests/test_paper59_euclidean_dictionary.py`, `tests/test_pairwise_three_center.py` | Chemistry | 1 | elliptic | dimensionless |
| J4 | Fourth distinct centre (XY\|ZW): still one radial $k$-integral with the same two dispersion roots; no genus two at four centres (pre-registered P1, 2026-10-13; confirmed against McMurchie–Davidson; nearly a corollary of Paper 59) | **F** [PRE-REGISTERED, CONFIRMED; promoted from Obs 2026-10-15 when the paper home and guard landed] | Paper 59 `sec:four_centre`; record `docs/predictions/2026-10-13_four_centre_genus.md`; `tests/test_paper59_four_centre_genus.py` (values recomputed; radicand set equal to the two quadrics, symbolically). Prior art: Avery & Avery 2017 App. B (same momentum form, no genus reading) | Chemistry | 1 | elliptic | dimensionless |
| J5 | Bending (water, AH₂): $m$ released; $\lvert m\rvert = \Lambda$ opens at order $s^{2\Lambda}$; the opening order and labels are exact rationals (6/5, 5/6), the slopes exponentials in $R$ (no π) | **F** (order and labels; the slopes are measured, not forced) | Paper 32 `rem:multicenter_composition`; Paper 59 `sec:scope`. Prior art: Renner 1934, Walsh 1953 (the CHANGELOG records "not novel as a method"). `tests/test_renner_teller_m_splitting.py`, `tests/test_renner_teller_rung2.py`, `tests/test_renner_teller_rung4.py`, `tests/test_walsh_ah2.py`, `tests/test_paper32_composition_pins.py` | Chemistry | 1 | none | dimensionless |
| J6 | Water's three-centre ERIs: shared-centre (Barnett–Coulson) channels converge to $\sim 10^{-9}$ by $L \approx 12$; a closed form meets the J3 wall | **W** (closed forms of the radial core only; the evaluation is accurate) | Paper 59 `sec:scope`. Prior art: Barnett–Coulson 1951 / Löwdin α-function. `tests/test_pairwise_general.py` | Chemistry | 1 | elliptic | dimensionless |
| J7 | Three-electron chain $\langle\rho_a\rho_b\rho_c/(r_{12}r_{13})\rangle$: a second momentum integration (three two-centre densities, three roots); the fibre is a K3 (P² on six lines); the holomorphic period is the Bessel moment $\bar V_3 = \int K_0K_0K_0$ (BBBG 2008 eq. 33, from an identity of Bailey), a sum of products $K\,K'$ ("elliptic squared"); at equal scales the fibre collapses to K3 periods plus π³; at generic scales a PSLQ search at one mass triple finds no relation in a K3, elliptic and Tate basis that omits dilogarithms and products of elliptic periods (evidence, not proof); the one-centre limit is rational | **F** for the K3 fibre and its elliptic-squared period [PRE-REGISTERED; P3a–d HOLD, guarded; promoted from Obs 2026-10-15]; the equal-scale collapse (P3e-1) and the generic content beyond K3 (P3e-2) stay evidence, unguarded | Paper 59 `sec:three_electron`; record `docs/predictions/2026-10-13_three_electron_k3.md`; `tests/test_paper59_three_electron_k3.py` (fibre agreement, elementary angle integral, rational one-centre limits, the closed form vs quadrature, K3 period vs $\bar V_3$, six-line branch locus, non-isotrivial slice); CHANGELOG v5.20.14–15. Prior art: Bailey 1936; Bailey–Borwein–Broadhurst–Glasser 2008 eq. 33; the three-loop banana literature. GeoVac's contribution is the reduction of the Slater chain to $\bar V_3$, not the closed form | Chemistry | >1 (two momenta) | K3 | dimensionless |
| J8 | Explicit-$r_{12}$ four-body bridge $f_{12}f_{34}/r_{13}$: four electrons linked by one bridging Coulomb kernel; reduces RI-free to a finite $(L,M)$ contraction of three-body vertex kernels (terminating at $L \le 2(l_{\rm bridge}+l_{\rm leaf})$) and, for $s$ leaves, to a two-electron Slater–Condon integral of dressed densities; $N \ge 5$ conjectural | **Obs** (validated reduction, no derivation; audit 2026-09-22) | Paper 12 `sec:r12`; `debug/sprint_n4_wall_diagnostic_memo.md`. Over two-centre densities the unreduced bridge is a four-density tree (kernels 1-2, 1-3, 3-4), so its holomorphic period is the four-density rung of J10, through J10's assumed three-dimensional reduction; the reduction to geminal-dressed densities, validated rather than derived, changes the cost, not the value. What stays Obs is the reduction itself. Prior art: one-centre Hylleraas-CI (King 2014, Sims–Hagstrom 2015, Ruiz–Sims–Padhy 2021); only the two-centre extension is the corpus's. `tests/test_lih_r12ci_bridge.py`, `tests/test_paper12_r12.py` | Chemistry | 1 | CY-type over two-centre densities (via J10) | dimensionless |
| J9 | Rule: the transcendence class of a two-electron integral with a rational kernel is set by the number of distinct two-centre transition densities the kernel couples (zero or one: genus zero; two: genus one; two densities on the same centre pair, the two-centre exchange class, close in the two-centre seed set of Paper 18 although they build the same fibre family, so the fibre genus alone does not decide the class), not by the number of nuclei; three densities over two momenta give a surface (K3), not a higher-genus curve, and four give the first rung not built from elliptic periods (J10); this ordering holds for trees, and a Coulomb ring of three densities sits on the elliptic rung while the correlated triangle with exponential geminals sits on the K3 rung at $\gamma=1$ (J11) | **Obs** [a pattern across the tested instances, not a derivation: symbolic at one centre, test-backed multicentre; re-tiered from F 2026-10-16 by the /qa claim-search DELTA] | Paper 59 `sec:f12`, `sec:momentum`, `sec:genus`; CHANGELOG v5.20.14. `tests/test_routeC_momentum.py` (`test_f12_multicenter_shares_eri_elliptic_genus`, `test_f12_geminal_pair_is_elementary`) | Chemistry | 1 | meta | N/A |
| J10 | Four linked densities: the holomorphic period of any tree of $n$ two-centre densities joined by $n-1$ rational kernels is $2^n/(\pi\sqrt{\prod c_v})\int_0^\infty\prod_v K_0(a_vx)\,dx$, independent of the tree's shape; the one-scale family $F_4(t)=\int K_0(x)^3K_0(tx)\,dx$ satisfies an irreducible fourth-order Calabi–Yau-type operator (AESZ #58, a Hadamard product) with no quadratic relation among its solutions and a non-constant Yukawa coupling, so its monodromy is Sp4 or SL4 and it is not built from elliptic periods; five densities give an orthogonal rank-five structure and no new type | **F** for the tree identity and the non-elliptic four-density operator [PRE-REGISTERED 2026-10-16, written before computing but not committed beforehand; guarded]; the reduction of the real three-dimensional four-density integral into the one-dimensional moment family is assumed, and $n = 5$ gives no new type (its operator is printed with AESZ #39, prior art); for two to five densities the tree ladder is the banana ladder under a Hadamard product with the Legendre family (Almkvist–Zudilin 2007; Vanhove 2014), so the operators are prior art and the route from a molecular integral is the contribution | Paper 59 `sec:four_density`; record `docs/predictions/2026-10-16_linking_ceiling.md`; `tests/test_paper59_four_density_rung.py` (tree identity, star and chain; Mellin-route operator vs eq:four_density_op and vs quadrature; quadric counts 0, 1, 0, 1 with controls; irreducibility; Yukawa); CHANGELOG v5.20.20. Prior art: AESZ tables (#58); Broadhurst 2009 (random walks, Bessel moments, CY equations); Fresán–Sabbah–Yu 2023 (symplectic/orthogonal alternation); Bourjaily–McLeod–von Hippel–Wilhelm 2019 (a ceiling set by loop order). GeoVac's contribution is the route from a molecular integral to the operator, not the operator | Chemistry | >1 (three momenta) | CY-type (Sp4) | dimensionless |
| J11 | Ring of three densities: the Coulomb triangle $\langle\rho_a\rho_b\rho_c/(r_{12}r_{13}r_{23})\rangle$, every pair of three densities joined by a kernel; the loop momentum integrates out exactly by the three-dimensional uniqueness (star-triangle) relation, leaving the chain's two momenta and three densities; at the fingerprint fibre the triangle's operator factors into elementary, algebraic ($(1-t^2/4)^{-3/2}$) and elliptic pieces, the last equivalent to the Legendre equation of the two-density period, so the ring sits on the elliptic rung, below the three-density chain's K3; at one centre it carries dilogarithms where the chain is rational | **F** for the loop reduction, the operator and its exact factorization, and the one-centre closed form [PRE-REGISTERED 2026-10-17 and committed before computing; guarded]; the registered three-way decision rule did not anticipate a landing below the chain (the record states both readings); the Coulomb kernels' drop is read as a property of uniqueness, not of rings (an Observation): the correlated triangle with exponential geminals on two edges, whose massive loop is a solid angle but not unique, sits with the chain on the K3 rung at $\gamma=1$ (exact peeling and map to $F_3$'s operator; its annihilator measured on the function; pre-registered as elliptic and falsified, record `docs/predictions/2026-10-17_ring_extensions.md`; guarded), and reaches the K3 rung at $\gamma=3/2$ by local data mod $p$ plus rigidity, whether it climbs further being open; the chain control is settled: its only non-elementary piece is the K3 family of $F_3$ (exact map), the Mellin route's top order-two factor is spurious, and the K3 piece provably does not split off the elementary ones (an exponent bound plus an exhaustive exact search) | Paper 59 `sec:density_ring`; record `docs/predictions/2026-10-17_density_loop.md`; `tests/test_paper59_density_ring.py` (star-triangle constant by Feynman parameters; one-centre closed form against the reduced formula and Remiddi's master; the (1,2,3) dilogarithm relation; $L_GG=-3t^2$ on the momentum form; exact right division and the rational map from Legendre's equation; $F_2$ solves that equation; the chain control's exact peeling, its K3 map, the spurious top factor and the bounded no-split search). Prior art: Remiddi 1991 (the one-centre dilogarithms); Fromm–Hill 1987; Pachucki–Puchalski–Remiddi 2004; the uniqueness relation (classical). GeoVac's contribution is the fibre reading, not the one-centre values | Chemistry | >1 (two momenta after the loop) | elliptic (fibre); weight-2 (one centre) | dimensionless |

---

## Summary counts

- **Forced (F):** 41 entries (several relative to a named organizing reading / adopted axiom, or CONDITIONAL, as marked in the rows)
- **Observation (Obs):** 4 entries (A3 κ, re-classified from F on 2026-10-02; J2, a measured front; J8, a validated reduction without a derivation; J9, the rung rule, a pattern across the tested instances)
- **Admitted (A):** 1 entry
- **Calibration (C):** 25 entries
- **Wall (W):** 3 entries (A6, J3, J6)

Total: 74 entries across ten domains (60 on the 2026-06-08 catalogue that Paper 57 §5–6 audited; 12 rows added 2026-10-15: A6, A7, I4, J1–J9; J10 added 2026-10-16; J11 added 2026-10-17).

## Preliminary discriminator read

Looking at the table cold, before §5's formal principle test:

- **Multi-focal depth (MF)** is the strongest single predictor. Every multi-focal-composition entry (G1–G6, plus E6, E7, E8, F1, F4, F5, F6, F7, D5, D6, I1; 17 entries) sits at MF > 1. Every gauge / spectral / forced entry sits at MF = 1, with four exceptions tagged MF > 1: C9 (the two-loop S_min period) and the released-axis rows J7, J10 and J11 (more than one momentum integration).
- **Period class (P)** (from the [CONDITIONAL] case-exhaustion, first fifteen projections) tracks the C-side multi-focal walls cleanly (all at "outside") but is the wrong axis for forced gauge content — the M1 tag is about the *prefactor* in the case-exhaustion theorem, not a *forcing mechanism* for gauge content.
- **Dimensional character (Dim)** correlates partially: most dimensionful entries sit on the C side, but D2 ($S_{\rm BH}$ / $G_N$) is F and dimensionful, while F2, F3, F5, F6 (N_gen, KO-dim, CKM, PMNS) are dimensionless and C.

**Two distinct families on the C side, by first read:**
- **Multi-focal composition** (G1–G6, E7, E8, F1, F4, F7, D5, D6): "framework can BUILD it but cannot autonomously pick a number." MF > 1, period outside, mostly dimensionful.
- **Inner-factor input data** (F2, F3, B8, H6, H7, I2, I3): "categorical to the framework's outer-factor mechanism; lives elsewhere by structure." MF = 1, period "none" or N/A, dimensionless.

This is the split already named in `memory/external_input_three_class_partition.md` (Classes 1/2/3). The catalogue confirms it as a quantitative pattern across the table: the calibration side decomposes into two structurally distinct families with different discriminator-column signatures.

The §5 candidate-principle probe in Paper 57 will therefore test principles *per family* in addition to globally.

---

*Source: this catalogue consolidates `docs/archive/forcing_catalogue.md` (2026-06-01; archived 2026-10-15), CLAUDE.md §1.7 working hypotheses, CLAUDE.md §3 dead-end rows, `debug/sprint_h1_*_memo.md`, `debug/sprint_yukawa_pslq_memo.md`, `debug/sprint_read2_n_gen_scoping_memo.md`, `debug/sprint_w1e_period_class_memo.md`, `debug/sprint_forced_count_synthesis_memo.md`, `debug/seam_packing_scoping_memo.md`, and the memory files indexed in MEMORY.md under "Working hypotheses and core structural findings" and "Sprint outcomes (durable findings)". Companion paper: `papers/group3_foundations/paper_57_forced_free_seam.tex`.*
