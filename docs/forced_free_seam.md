# The Forced/Free Seam

*A catalogue of what the GeoVac framework autonomously generates and what it consumes as external input.*

Started 2026-06-08. Supersedes `docs/forcing_catalogue.md` for the forced-vs-free question; the prior catalogue's forward-run / probe structure remains the canonical record for *what does the same rigidity force that we haven't looked at*.

---

## Purpose

GeoVac is bit-exact on a wide class of structural facts (Paper 0 packing, S³ spectrum, gauge group, dim M(D_F), master Mellin engine classification, …) and consumes specific numerical content as external input (Yukawa values, generation count, calibration cutoff f, atomic-physics Z_eff fits, …). The boundary between these two regimes — the **forced/free seam** — is itself a structural property of the framework, repeatedly named (`memory/multi_focal_wall_pattern.md`, `memory/external_input_three_class_partition.md`, `memory/geovac_structural_skeleton_scope_pattern.md`, Paper 32 §VIII Forced-Count Theorem, Paper 18 §IV.6 inner-factor input-data tier) but never consolidated.

This catalogue consolidates it. Companion paper: `papers/group3_foundations/paper_57_forced_free_seam.tex`.

## Status conventions

- **F (Forced)** — the framework autonomously determines this structural feature from the Paper 0 packing axiom + standard spectral-triple machinery. Witness: a theorem or bit-exact verification.
- **A (Admitted)** — the framework structurally allows this but does not autonomously select it among multiple structurally-equivalent options. Witness: existence proof + non-selection theorem or negative search.
- **C (Calibration)** — the framework consumes this as external input. Witness: a non-selection theorem (preferred) or a documented negative search (acceptable).

## Discriminator columns

For each entry we record three properties to be tested in the paper §5 candidate-principle probe:

- **MF** — multi-focal depth. `1` if the observable is the output of a single Paper-34 projection; `>1` if its definition requires composing multiple projections.
- **P** — period class. `M1` / `M2` / `M3` / `outside` / `none` per Paper 18 §III.7 master Mellin engine; `meta` for the engine itself; `N/A` for non-period observables.
- **Dim** — dimensional character. `dimensionless` (ratio) / `dimensionful` (scale).

---

## A. Foundations / skeleton

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| A1 | Graph Laplacian spectrum $\lambda_n = n^2-1$ | **F** | Packing axiom; Paper 0 + Paper 7 (18 symbolic proofs) | Skeleton | 1 | none | dimensionless |
| A2 | $S^3$ as forced outer manifold | **F** | Bertrand + SO(4) closure; Paper 7 | Skeleton | 1 | none | dimensionless |
| A3 | $\kappa = -1/16$ universal kinetic scale | **F** | Fock projection; exact rational | Skeleton | 1 | none | dimensionless |
| A4 | Shell degeneracy $g_n = n^2$ | **F** | Representation theory; exact integer | Skeleton | 1 | none | dimensionless |
| A5 | Discreteness ↔ compactness | **F** | Peter–Weyl; `memory/geovac_structural_skeleton_scope_pattern.md` | Skeleton | 1 | N/A | N/A |

## B. Gauge structure

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| B1 | Gauge group $U(1) \times SU(2) \times SU(3)$ | **F** | Bertrand × Hopf-tower truncated at $n \le 3$; Paper 32 §VIII.B | Gauge | 1 | M1 | dimensionless |
| B2 | Inner algebra factor count (= 3) | **F** | Hurwitz + Hopf-tower; Door 4b/4d | Inner | 1 | none | dimensionless |
| B3 | $\mathbb{C}$ factor at $n=1$ | **F** | Door 4d (unique division algebra of $S^1$) | Inner | 1 | none | dimensionless |
| B4 | $M_3(\mathbb{C})$ factor at $n=3$ | **F** | Door 4d (Hurwitz fallback at $S^5$) | Inner | 1 | none | dimensionless |
| B5 | $\mathbb{H}$ vs $M_2(\mathbb{C})$ at $n=2$ | **A** | Door 4c NEGATIVE on $J$-sign-table; ℍ admitted, not forced; Door 4e Upgrade B closes if adopted | Inner | 1 | none | dimensionless |
| B6 | $\dim \mathcal{M}(D_F) = 128$ per generation | **F** | Forced-Count Theorem, Paper 32 §VIII (`thm:forced_count`); $1024 \to 512 \to 128$ chain | Inner | 1 | none | dimensionless |
| B7 | Higgs admission (Mexican-hat structurally allowed) | **F** | Sprint H1 POSITIVE-THIN; Paper 32 §VIII.C | Inner | 1 | none | dimensionless |
| B8 | Lower-bound forcing on gauge content (saturated, not subset) | **C** | Bertrand × Hopf gives upper bound only; structural-skeleton scope pattern | Gauge | 1 | none | dimensionless |

## C. Spectral / period content

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| C1 | Master Mellin engine $\mathcal{M}[\mathrm{Tr}(D^k \cdot e^{-tD^2})]$, $k \in \{0,1,2\}$ | **F** | Case-exhaustion theorem; Paper 32 §VIII | Periods | 1 | meta | dimensionless |
| C2 | M1 prefactor = $\mathrm{Vol}(S^2)/4 = \pi$ | **F** | Hopf-base measure; M1 signature | Periods | 1 | M1 | dimensionless |
| C3 | Seeley–DeWitt $a_0 = a_1 = \sqrt\pi$, $a_2 = \sqrt\pi/8$ on $S^3$ | **F** | Closed form; Paper 28 T9 | Periods | 1 | M2 | dimensionless |
| C4 | $\zeta_{D^2}(s) = 2^{2s-1}[\lambda(2s-2) - \lambda(2s)]$ | **F** | T9 theorem; Paper 28 | Periods | 1 | M2 | dimensionless |
| C5 | $\chi_{-4}$ identity $D_{\rm even} - D_{\rm odd} = 2^{s-1}(\beta(s) - \beta(s-2))$ | **F** | Paper 28 §IV; vertex-parity Hurwitz | Periods | 1 | M3 | dimensionless |
| C6 | $M_3 \subset \mathcal{MT}(\mathbb{Z}[i, 1/2])$ at level $\le 4$ | **F** | Sprint Mixed-Tate; Paper 55 | Periods | 1 | M3 | dimensionless |
| C7 | Universal propinquity rate $4/\pi$ | **F** | WH1 PROVEN; Paper 38 + Paper 40 (rank-invariant across compact Lie groups) | OA | 1 | M1 | dimensionless |
| C8 | Pythagorean orthogonality $\langle H_{\rm local}, D_W^L\rangle_{\rm HS} = 0$ | **F** | Bit-exact; Paper 43 §10.2 | OA | 1 | M1 | dimensionless |
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
| E5 | $B = 42$ (Casimir), $F = \pi^2/6$ (Fock Dirichlet), $\Delta = 1/40$ ($g_3^{\rm Dirac}$) | **F** | Three independent spectral homes proven; 12 mechanisms eliminated | $\alpha$ | 1 | M1/M2/M2 | dimensionless |
| E6 | Combination rule $K = \pi(B + F - \Delta)$ for $\alpha^{-1}$ | **C** | Numerical coincidence $8.8 \times 10^{-8}$; no derivation; stays conjectural per CLAUDE.md §13.5 | $\alpha$ | >1 | M1+M2 | dimensionless |
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
| H1 | Angular sparsity 1.44% depends only on $l_{\max}$ | **F** | Paper 22 potential-independent theorem | Chemistry | 1 | none | dimensionless |
| H2 | $N_{\rm Pauli} = 11.10 \times Q$ across 38 molecules | **F** | Paper 14; isostructural invariance | QC | 1 | none | dimensionless |
| H3 | Pauli ratio balanced/composed $O(B^2)$: 2.63 / 4.77 / 7.45 | **F** | Proven structural | QC | 1 | none | dimensionless |
| H4 | Nuclear magic numbers 2, 8, 20, 40, 70, 112 | **F** | HO shell closures from graph counting; Paper 23 | Nuclear | 1 | none | dimensionless |
| H5 | Fock projection rigidity ($S^3$ unique to $-Z/r$) | **F** | Paper 23 NJ theorem | Nuclear | 1 | none | dimensionless |
| H6 | Clementi–Raimondi $Z_{\rm eff}$ exponents | **C** | External fits; W1e period-class memo classifies as inner-factor input-data tier | Chemistry | 1 | outside | dimensionless |
| H7 | Multi-zeta coefficients in second-row valence basis | **C** | Inherits H6 | Chemistry | 1 | outside | dimensionless |

## I. Foundational calibration (Class 1)

| # | Name | Status | Mechanism / Witness | Domain | MF | P | Dim |
|---|------|:------:|---------------------|--------|:--:|:--:|:--:|
| I1 | Fine-structure constant $\alpha$ (value) | **C** | E6 combination rule conjectural; $\alpha$ value external | $\alpha$ | >1 | M1+M2 | dimensionless |
| I2 | Born rule probability rule $p = \lvert\langle a\lvert\psi\rangle\rvert^2$ | **C** | Gleason via Hilbert inheritance; framework does not improve on Gleason; `external_input_three_class_partition.md` Class 1 | Foundations | 1 | N/A | N/A |
| I3 | Higgs direction $\hat n \in S^2$ | **C** | Boyle–Farnsworth input; possible Hopf-base identification flagged | Inner | 1 | M1 | dimensionless |

---

## Summary counts

- **Forced (F):** 38 entries
- **Admitted (A):** 1 entry
- **Calibration (C):** 23 entries

Total: 62 entries across nine domains.

## Preliminary discriminator read

Looking at the table cold, before §5's formal principle test:

- **Multi-focal depth (MF)** is the strongest single predictor. Every multi-focal-composition wall entry (G1–G6, plus E6, E7, E8, F1, F4, F5, F6, F7, D5, D6) sits at MF > 1. Every gauge / spectral / forced entry sits at MF = 1.
- **Period class (P)** tracks the C-side multi-focal walls cleanly (all at "outside") but is the wrong axis for forced gauge content — gauge group at B1 carries an M1 signature in its derivation, but the M1 classification is about the *prefactor* in the case-exhaustion theorem, not the *forcing mechanism*.
- **Dimensional character (Dim)** correlates partially: most dimensionful entries sit on the C side, but D2 ($S_{\rm BH}$ / $G_N$) is F and dimensionful, while F2, F3, F5, F6 (N_gen, KO-dim, CKM, PMNS) are dimensionless and C.

**Two distinct families on the C side, by first read:**
- **Multi-focal composition** (G1–G6, E7, E8, F1, F4, F7, D5, D6): "framework can BUILD it but cannot autonomously pick a number." MF > 1, period outside, mostly dimensionful.
- **Inner-factor input data** (F2, F3, B8, H6, H7, I2, I3): "categorical to the framework's outer-factor mechanism; lives elsewhere by structure." MF = 1, period "none" or N/A, dimensionless.

This is the split already named in `memory/external_input_three_class_partition.md` (Classes 1/2/3). The catalogue confirms it as a quantitative pattern across the table: the calibration side decomposes into two structurally distinct families with different discriminator-column signatures.

The §5 candidate-principle probe in Paper 57 will therefore test principles *per family* in addition to globally.

---

*Source: this catalogue consolidates `docs/forcing_catalogue.md` (2026-06-01), CLAUDE.md §1.7 working hypotheses, CLAUDE.md §3 dead-end rows, `debug/sprint_h1_*_memo.md`, `debug/sprint_yukawa_pslq_memo.md`, `debug/sprint_read2_n_gen_scoping_memo.md`, `debug/sprint_w1e_period_class_memo.md`, `debug/sprint_forced_count_synthesis_memo.md`, `debug/seam_packing_scoping_memo.md`, and the memory files indexed in MEMORY.md under "Working hypotheses and core structural findings" and "Sprint outcomes (durable findings)". Companion paper: `papers/group3_foundations/paper_57_forced_free_seam.tex`.*
