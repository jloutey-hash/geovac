# Track TS-B: Spectral-Triple ↔ Projection-Taxonomy Dictionary

**Sprint:** TS-B (research memo, not a paper)
**Date:** 2026-05-04
**Goal:** Bridge the spectral-triple framing (Paper 32 / Paper 31, Connes–Marcolli–vS lineage) with the projection-taxonomy framing (Paper 34 / Paper 35, GeoVac-native two-layer architecture). Produce a 15-row dictionary table, then a synthesis identifying patterns, gaps, asymmetries, and predictions. The asymmetries are the deliverable that feeds Track C investigation seeds.

**Companion files:**
- `debug/tx_a_dimensional_axiomatization_memo.md` (TX-A, 15-row commutation analysis)
- `debug/tx_b_paper35_test_memo.md` (TX-B, Paper 35 falsification panel)
- `debug/data/tx_a_projection_table.json` (per-projection axis flags)
- `debug/data/tx_a_composition_table.json` (15×15 relations)

**Source papers (read for this memo):** Paper 34 (projection list, §III.1–§III.15), Paper 35 (rest-mass + observation projections), Paper 31 (universal/Coulomb partition), Paper 32 (spectral-triple synthesis, §3 construction, §4 axiom audit, §5 sub-sectors, §6 K-decomposition, §7 Coulomb/HO asymmetry), Paper 18 (transcendental-tier names: intrinsic / calibration / embedding / flow / composition).

---

## 0. Setup: what the two framings provide

**Spectral-triple side (Papers 31, 32).**
GeoVac is a finite spectral triple $(\mathcal{A}_{\mathrm{GV}}, \mathcal{H}_{\mathrm{GV}}, D_{\mathrm{GV}})$ on the Fock-projected $S^3$ graph (Paper 32 §3):
- $\mathcal{A}_{\mathrm{GV}} = \mathbb{C}^{V_{\mathrm{Fock}}}$ — commutative function algebra on graph nodes (Paper 32 Definition 1).
- $\mathcal{H}_{\mathrm{GV}} = \mathcal{H}^{\mathrm{scalar}} \oplus \mathcal{H}^{\mathrm{spinor}}$ — scalar Fock states + Camporesi–Higuchi spinor sector (Paper 32 §3.2).
- $D_{\mathrm{GV}}$ — Camporesi–Higuchi Dirac operator with eigenvalues $|\lambda_n^D| = n + 3/2$ and degeneracies $g_n^D = 2(n+1)(n+2)$ (Paper 32 §3.3).

Connes axioms audited at finite $n_{\max}$ (Paper 32 Table 1):
- **Order zero** (faithful $*$-rep): holds rigorously.
- **Bounded commutators** $[D, \pi(a)]$: holds rigorously.
- **Compact resolvent / 3⁺-summability**: automatic at finite truncation; CH summability in continuum.
- **Order one** $[[D, \pi(a)], \pi(b)^\circ] = 0$: holds *trivially* because $\mathcal{A}_{\mathrm{GV}}$ is abelian.
- **Real structure** $J$: holds via Camporesi–Higuchi charge conjugation; signs $(\epsilon, \epsilon') = (-1, +1)$ for KO-dim 3.
- **$\mathbb{Z}_2$-grading $\gamma$**: *absent* (correct: $S^3$ is odd, the triple is odd).

Universal/Coulomb partition (Paper 31): features depending only on $\mathcal{A}$ are *universal* (transfer to any spherical-fermion system); features depending on the specific $D$ are *potential-specific*.

**Projection-taxonomy side (Paper 34, expanded by Paper 35).**
Layer 1 = bare combinatorial graph (π-free, integer/rational matrix elements). Layer 2 = a finite family of named projections, each adding (variable, dimension, transcendental). The fifteen documented projections are the rows of the dictionary below. TX-A established that of the three axes, only the transcendental axis is genuinely independent (V/D perfectly correlated). TX-B graduated Paper 35 Prediction 1 — π enters iff continuous integration over a temporal/spectral parameter promoted from the discrete spectrum — to load-bearing principle (208/208 a priori).

---

## 1. The 15-row dictionary

Column legend.
- **Sector modified** (column 2): which piece of $(\mathcal{A}, \mathcal{H}, D)$ — or none of them — does the projection touch?
  - `A` = extends or restricts $\mathcal{A}_{\mathrm{GV}}$
  - `H` = tensors / restricts / re-bases $\mathcal{H}_{\mathrm{GV}}$
  - `D` = deforms $D_{\mathrm{GV}}$ (changes spectrum, changes operator)
  - `B` = adds boundary structure (truncation, periodic time, temporal compactification)
  - `O` = observable-extraction map outside the triple (e.g. spectral-action evaluation, regularization flow)
- **Connes axioms** (column 3): which axioms are preserved (`+`) vs. broken (`−`) vs. trivially preserved (`triv`) vs. not applicable (`n/a`). Five-letter code for {Reality $J$, $\gamma$-grading, order-one O1, orientation, KO-dim}.
- **Universal vs. Coulomb-specific** (column 4): from Paper 31's A/D split. `U` = universal (depends only on $\mathcal{A}$ + angular structure); `C` = Coulomb-specific (depends on the specific $D$ or on the Coulomb projection); `H` = HO-specific; `*` = either-side or ambiguous.
- **Transcendental class** (column 5): from Paper 18's tier names + Paper 35's new tiers. Categories used: rational / algebraic-extension / intrinsic / calibration-π / odd-ζ-Dirac / Catalan G + Dirichlet L / vertex-topology / observation-2π / flow / composition.
- **TX-A axis tag** (column 6): (V, D, T) each + or –. Match to `tx_a_projection_table.json`.
- **Commutation flags** (column 7): non-commuting partners (from TX-A 22 one-way pairs + 6 variable conflicts + 3 idempotents + 3 inverse pairs).

Note on column-3 encoding: the *order-one* axiom is trivially preserved (`triv`) for any projection that leaves $\mathcal{A}_{\mathrm{GV}}$ commutative, because trivial commutativity makes the order-one bracket vanish identically (Paper 32 Proposition 5). Order-one becomes substantive only under non-abelian almost-commutative extension — which is where Wilson plaquette and the spectral action sit.

| # | Projection (Paper 34 §III) | Sector | Axioms preserved | Univ./Coul. | Transcendental class | (V, D, T) | Non-commuting partners / inverse / idempotent |
|---:|:---|:---|:---|:---:|:---|:---:|:---|
| 1 | Fock conformal $(\mathbb{R}^3 \to S^3)$ | $\mathcal{A}, \mathcal{H}, D$ (defines the triple) | $J^+$, $\gamma^{\mathrm{n/a}}$, O1$^{\mathrm{triv}}$, KO-dim $3$ | C | calibration-π via $\mathrm{Vol}(S^3)$, rational $\kappa = -1/16$ | (+, +, –) | one-way w/ drake_swainson; commutes w/ all 13 others |
| 2 | Hopf bundle $S^3 \to S^2 \times S^1$ | $\mathcal{H}, D$ (decomposition into base+fibre; injects U(1) inner-derivation algebra) | $J^+$, O1$^{\mathrm{triv}}$ at scalar level; **breaks O1 substantively in $\mathcal{A}\otimes M_n$ extension** | C | calibration-π $= \mathrm{Vol}(S^2)/4$ | (+, –, +) | **variable_conflict w/ spectral_action AND spinor_lift** (all three carry α); one-way w/ drake_swainson |
| 3 | Bargmann–Segal $\mathbb{R}^3_{\mathrm{HO}} \to S^5$ Hardy | $\mathcal{A}', \mathcal{H}', D'$ — defines a *distinct* triple (Paper 31 Table 2) | $J$ status open (Paper 32 §11), $\gamma^{\mathrm{n/a}}$ at scalar, O1$^{\mathrm{triv}}$ | H | rational (π-free certificate, Paper 24) | (+, +, –) | one-way w/ drake_swainson |
| 4 | Stereographic / conformal coord. change $S^n \to \mathbb{R}^n$ | $\mathcal{A}$ (re-bases function space, preserves $D$ up to conformal factor) | preserves all (conformal, not new $D$) | * | calibration: conformal factor producing Coulomb $1/r$ | (+, +, –) | inverse pair (self-inverse on open dense subset); one-way w/ drake_swainson |
| 5 | Sturmian reparam at $\lambda = Z/n$ | $\mathcal{A}, \mathcal{H}$ (re-labels same Fock graph) | preserves all | C | rational (preserves Layer 1 ring); flow-tier transcendentals enter via subsequent Bethe-log sums | (–, –, –) | **idempotent**; one-way w/ drake_swainson |
| 6 | Connes–Chamseddine spectral action $\mathrm{Tr}\, f(D^2/\Lambda^2)$ | O (observable extraction over triple); also $\mathcal{A}\otimes M_n$ when gauge sector included | preserves $J$, $\gamma$, O1; introduces UV cutoff $\Lambda$ as auxiliary | C | calibration: $\sqrt{\pi}\cdot\mathbb{Q}$ SD coeffs → $\pi^{2k}\cdot\mathbb{Q}$ observables | (+, +, +) | **variable_conflict w/ hopf AND spinor_lift** (all carry α); one-way w/ drake_swainson |
| 7 | Camporesi–Higuchi spinor lift | $\mathcal{H}$ (extends scalar to spinor sector); makes $D_{\mathrm{GV}}$ act on spinor fibre | activates $J$ (CH conjugation, $J^2 = -1$ for KO-dim 3); preserves O1$^{\mathrm{triv}}$ at abelian core; **engages $\gamma$ choice (odd, not present at $S^3$)** | C | spinor-intrinsic odd-ζ + Catalan G + Dirichlet $\beta(4)$; ring $R_{\mathrm{sp}} = \mathbb{Q}(\alpha^2)[\gamma]/(\gamma^2 + (Z\alpha)^2 - 1)$ | (+, –, –) | **variable_conflict w/ hopf AND spectral_action**; one-way w/ drake_swainson |
| 8 | Wigner $3j$ angular coupling | (none — internal to $\mathcal{A}$'s SO(3) decomposition) | preserves all | U | rational $\mathbb{Q}[\sqrt{2k+1}]$ | (–, –, –) | **idempotent**, **inverse pair** (CG ↔ decoupling); one-way w/ drake_swainson |
| 9 | Wigner $D$-matrix rotation $A \to B$ | $\mathcal{A}$ (relabels nodes between molecular centres) | preserves all | U | rational $\mathbb{Q}[\sqrt{2}, \sqrt{3}, \sqrt{6}]$ | (+, +, –) | **inverse pair** (rotation by $R^{-1}$); one-way w/ drake_swainson |
| 10 | Wilson plaquette (edge → 2-cell) | $\mathcal{A}\otimes M_n$ (almost-commutative extension; $n=2$ for SU(2)); also adds 2-cells over $\mathcal{H}$ | **breaks O1 substantively** (the canonical place where order-one constrains gauge fluctuations); preserves $J$ choices | C (full SU(2) only on Coulomb $S^3$; only U(1) survives on Bargmann $S^5$) | calibration via SU(2) Haar; transverse photon propagator | (+, –, +) | one-way w/ drake_swainson; commutes w/ all others |
| 11 | Vector-photon promotion (scalar 1-cochain → vector $(q, m_q)$) | $\mathcal{H}$ (lifts photon 1-cochain to vector-harmonic) | preserves all | C (recovers 7/8 selection rules — Paper 33 §V) | calibration $1/(4\pi)$ per loop = $S^2$ Weyl exchange constant of Hopf base | (–, –, +) | one-way w/ drake_swainson; **unique witness of pure-transcendental axis (no V, no D, only T)** |
| 12 | Molecule-frame hyperspherical | $\mathcal{H}, D$ (re-bases multi-electron Hilbert space onto $(\rho, \alpha, \theta_{12}, \ldots)$) | preserves all at angular level; $D$ becomes piecewise-smooth in $R$ | C (composed architecture of Paper 17) | embedding tier; Gaunt-rational angular content | (+, +, –) | one-way w/ drake_swainson |
| 13 | Drake–Swainson asymptotic subtraction | O (observable extraction; introduces transient $K$ that cancels) | preserves all | C (closes the 2P Bethe-log floor; Coulomb-Sturmian content) | flow tier; rational structural denominator $D_{\mathrm{drake}} = 2(2\ell+1)Z^4/n^3$ | (+, –, –) | **idempotent**; **one-way reverse w/ every other projection (11 cells: drake → all 11 non-idempotent others, never the other direction)** |
| 14 | Rest-mass projection (Paper 35 §VI; new in TX-B) | $D$ (additive shift $\omega^2 \to \omega^2 + m^2$); does not touch $\mathcal{A}$ | preserves all | * (universal across central potentials; relativistic mass is generic) | **trivial / ring-preserving** (Paper 35 KG-1: 200/200 cases, $m^2 \in \mathbb{Q} \Rightarrow \omega \in \mathbb{Q}[\sqrt d]$) | (+, +, –) | one-way w/ drake_swainson; commutes w/ all 13 others |
| 15 | Observation / temporal-window projection (Paper 35 §VI; new in TX-B) | B (adds periodic Euclidean time $S^1_\beta$; equivalently a finite-time observation window) | preserves $J$, KO-dim of spatial sector; introduces fermionic vs bosonic Matsubara distinction (boson $2\pi k/\beta$, fermion $(2k+1)\pi/\beta$) | C/U mixed (heated triple is generic; partition function structure is universal) | observation-2π: $2\pi\cdot\mathbb{Q}$ per Matsubara mode in spectrum; $\pi^{2k}\cdot\mathbb{Q}$ in integrated quantities (Stefan–Boltzmann $\pi^2/90$) | (+, +, +) | one-way w/ drake_swainson |

---

## 2. Synthesis

### 2.a Patterns: do `A`-projections cluster into the universal sector?

The hypothesis to test was: *all "modify $\mathcal{A}$" projections sit in Paper 31's universal sector; all "modify $D$" projections are Coulomb-specific.* The dictionary's column-2/column-4 cross-tabulation gives a partial confirmation with three explicit exceptions.

Cross-table (15 rows):

| Sector touched | Universal (`U`) | Coulomb-specific (`C`) | Mixed (`*`) | HO-specific (`H`) |
|:--|:-:|:-:|:-:|:-:|
| `A` only or `A`-relabel | **2** (Wigner 3j #8, Wigner D #9) | 0 | 1 (stereographic #4) | 0 |
| `H` only or `A, H` | 0 | 3 (spinor lift #7, vector photon #11, molframe #12) | 0 | 0 |
| `D` (or `D`-deforming) | 0 | 1 (Hopf #2, decomposes $D$ into base+fibre) | 1 (rest mass #14) | 0 |
| `A, H, D` (defines triple) | 0 | 1 (Fock conformal #1) | 0 | 1 (Bargmann–Segal #3) |
| `O` (observable extraction) | 0 | 3 (spectral action #6, Wilson #10, drake #13) | 0 | 0 |
| `B` (boundary) | 1 (observation #15: partition function structure is generic) | 0 | 0 | 0 |
| `C` (Coulomb) total | — | **8** | — | — |

The hypothesis *holds with three caveats*:

1. **`A`-only projections are exactly the two universal cases (Wigner 3j #8, Wigner D #9).** Stereographic #4 is `*` because the SO(3) sub-action on angular labels is universal but the radial $1/r$ specialization is Coulomb-specific. The hypothesis predicts these two as universal; confirmed.

2. **`D`-touching projections are not all Coulomb-specific.** Rest-mass projection #14 is `D`-deforming (additive $\omega^2 \to \omega^2 + m^2$) but is *universal across central potentials* — it works on Coulomb $S^3$ exactly as well as on HO $S^5$ or any other central system, and Paper 35 verified the ring-preserving property generically. This is the dictionary's first counter-example to the hypothesis: a `D`-projection that lives in Paper 31's universal sector. Mechanism: an *additive scalar shift* of $D$ that preserves the spectrum's algebraic class is an $\mathcal{A}$-level property in disguise, because the shift is by an element of the centre $\mathbb{C}\cdot\mathbb{1}$, which is contained in $\mathcal{A}_{\mathrm{GV}}$.

3. **Observation/temporal-window #15 is `B`-modifying (boundary, not $\mathcal{A}/\mathcal{H}/D$).** This is structurally outside the spectral-triple data $(\mathcal{A}, \mathcal{H}, D)$ — the triple does not natively encode the temporal direction; periodicity in Euclidean time is a *trace boundary condition* on $\mathcal{H}$ rather than a modification of any of the three pieces. This forces a modest extension of the spectral-triple framing: Paper 32 audits the triple at fixed-spatial-slice; finite-temperature observables live in a different category (e.g., a Banach-algebra crossed product or a spectral triple with a separate temporal action). The dictionary correctly flags this as non-trivially related to the universal/Coulomb partition: the *partition-function structure* (boson/fermion Matsubara) is universal (any $D$ on any spherical sector accepts a thermal trace); the *specific transcendental content* (Stefan–Boltzmann $\pi^2/90$) depends on the spatial spectrum, hence is Coulomb-specific in the sense of Paper 31.

The pattern is therefore: $\mathcal{A}$-only projections are universal (consistent with the hypothesis); `D`-substantive projections are Coulomb-specific (consistent), with rest-mass as the witness that *trivial* `D`-modifications can stay universal; `O` (observable-extraction) projections inherit the Coulomb-specificity of the underlying $D$ (consistent); `B` (boundary) is a new sector the spectral-triple language does not natively cover.

### 2.b Gaps

**Connes-axiom gap.** Of the five Connes axioms (R, $\gamma$, O1, orientation, KO-dim), the only one *substantively broken* by any projection in the dictionary is order-one — and only in the almost-commutative extensions (Hopf #2 with $\mathcal{A}\otimes\mathbb{C}$ at U(1), Wilson plaquette #10 with $\mathcal{A}\otimes M_2$ at SU(2), spectral action #6 when the gauge sector is engaged). All other Connes axioms are preserved by every projection, either substantively or trivially. This is a strong structural statement: the projection family is engineered to live entirely inside the Connes axiomatic framework, with order-one as the single substantive axiomatic discriminator.

**Sector-coverage gap.** No projection touches the *grading* $\gamma$ — because $S^3$ is odd and the GeoVac triple is correctly classified as odd (Paper 32 Proposition 6). But the framework has not yet exhibited a projection that *evens* the triple by tensoring with a $\mathbb{Z}_2$-graded fibre (e.g., $\mathcal{A}\otimes\mathrm{Cl}_1$). This is the closest thing to a "sector no projection touches" in the dictionary. The companion question — whether the framework should have a projection that introduces an even fibre — is a candidate Track C investigation seed (see §2.c.6 below).

**Universal-sector transcendental gap.** Among the universal-sector projections in the dictionary (Wigner 3j #8, Wigner D #9), neither introduces a transcendental beyond $\mathbb{Q}[\sqrt{2k+1}]$. This is consistent with Paper 31's prediction: the universal sector is determined by $\mathcal{A}$ and the SO(3) decomposition, neither of which carries calibration $\pi$. *No projection in the dictionary surprises by injecting transcendental content into the universal sector.* This is a clean confirmation of Paper 31, not a gap.

**Coulomb-specific without π gap.** Of the eight Coulomb-specific projections, two do *not* introduce π: spinor lift #7 (introduces odd-ζ + Catalan G + Dirichlet β instead) and Drake–Swainson #13 (flow tier, $D_{\mathrm{drake}}$ rational). This is consistent with Paper 18's tier taxonomy (spinor-intrinsic and flow are distinct from calibration-π) but it does say that Coulomb-specificity is *necessary but not sufficient* for π-injection. Vector-photon #11 is the converse case: a Coulomb-specific projection that is pure-transcendental (no V, no D), which by Paper 35 Prediction 1 must come from a continuous integration somewhere — and indeed, the $1/(4\pi)$ comes from the $S^2$-Weyl integration of the Hopf base.

### 2.c Asymmetries (Track C investigation seeds)

This section is the deliverable. Each entry below is a candidate place where the dictionary makes a structural prediction that has not been verified, or where the two framings disagree, or where a TX-A "innocent" projection has hidden spectral-triple content.

#### C-1. *The TX-A memo–data inconsistency*: spectral_action ↔ observation_window

**Asymmetry kind:** Internal inconsistency between TX-A memo (narrative) and TX-A data (JSON).

The TX-A memo (`debug/tx_a_dimensional_axiomatization_memo.md` §4.4) names `spectral_action ; observation_window` as the *most surprising entry* in the composition table, claiming both directions are flagged one-way and that this reproduces the field-theoretic UV/IR ordering anomaly ($\beta \to \infty$ does not commute with $\Lambda \to \infty$).

But the actual JSON data (`debug/data/tx_a_composition_table.json`) returns `commute` for both `spectral_action ; observation_window` and `observation_window ; spectral_action`. There are *zero* one-way pairs other than those involving `drake_swainson` — the 22 one-way cells in the JSON are *all* `drake_swainson ↔ X` for $X$ ranging over the other 11 non-idempotent projections (11 forward + 11 reverse).

The memo's claim is therefore *narrative*, not metadata-supported. The JSON's actual structure is much simpler than the memo describes: drake_swainson is the unique non-trivial Layer-2-output projection, and its asymmetric position in the layer ordering produces all 22 one-way cells; the other 14 projections all mutually commute at the metadata level. The headline "UV/IR ordering anomaly captured by signature arithmetic" is not in the JSON.

**Track C seed:** Re-derive the composition table with the layer ordering refined to distinguish *which* Layer 2 output a projection produces (UV-cutoff Layer 2 vs thermal-trace Layer 2 vs flow-regulated Layer 2). The expectation is that with a finer grading, spectral_action and observation_window genuinely *do* fail to commute at the metadata level — but the current TX-A grading is too coarse to see it. This Track C exercise would be: enrich the layer indices and re-run `tx_a_signature_arithmetic.py` to verify whether the memo's claim holds under any reasonable refinement.

#### C-2. The K-chain non-commutation: variable_conflict, not operator-system

**Asymmetry kind:** Metadata mismatch with structural reality.

The TX-A memo §3 cites `(hopf_bundle, spinor_lift)` as the *unique non-commuting pair* in the worked examples and identifies it with the WH1 R3.2 finding that the Connes–vS truncated operator system is sensitive to the order of Hopf-then-spinor vs spinor-then-Hopf.

The JSON encodes this differently. `hopf_bundle ; spinor_lift = variable_conflict`, not `one_way_AB`. The conflict is over the variable α — both projections introduce α, so the metadata flags the composition as ill-defined under the dictionary's "no double-introduction" rule. The same conflict exists for `(hopf_bundle, spectral_action)` and `(spinor_lift, spectral_action)`: all three projections that touch α are mutually variable_conflict.

This is a *real asymmetry between framings*, not a memo error. From the spectral-triple side, the operator-system order matters because Hopf and spinor act on different Hilbert-space tensor factors; the order of tensor factor introduction is not commutative for $\mathcal{H}$-modifications when the second factor acts non-trivially on the first. From the projection-taxonomy side, the order is invisible at the signature level (both orderings produce the same union of {α, …}) and the metadata correctly flags "two projections share the α slot, composition is ill-defined" — but not "operator-system ordering matters."

**Track C seed:** The variable-conflict mechanism captures a *necessary* condition for operator-system non-commutation (both projections must touch the same Hilbert-space factor) but not a *sufficient* one. Build the explicit operator-system non-commutator for each variable-conflict pair $(P_i, P_j)$ in the JSON and tabulate which of the three pairs produce genuine operator-system inequivalences (the WH1 R3.2 result for hopf+spinor) vs. which merely double-introduce a variable that could be re-named. This would refine the "variable_conflict" cell into "double-naming" vs "operator-system genuine non-commutation."

#### C-3. Drake–Swainson is the unique Layer-2 output

**Asymmetry kind:** Taxonomically singular, structurally unremarkable.

In the JSON, the *only* projection whose composition with every other projection is asymmetric (one-way) is drake_swainson. All 22 one-way cells involve drake_swainson. This makes drake_swainson the unique projection with input_layer = output_layer = 2 (according to the TX-A `tx_a_signature_arithmetic.py` layer assignment).

But spectral_action and observation_window also have output_layer = 2 in the JSON. So the asymmetric layer behavior of drake_swainson is *not* generic to "Layer 2 output" projections; it's specific to drake_swainson's role as a *flow* tier projection that consumes a regulated Layer 2 output and produces the same Layer 2 output (idempotent). The TX-A memo's framing of drake as "Layer 2 → Layer 2" idempotent that absorbs the K-cancellation is consistent with the JSON, but the *uniqueness* of drake's position should be recognized as structurally informative.

In spectral-triple language: drake_swainson is the only projection in the dictionary that operates *outside* the triple proper — it is a regularization-flow operation on observables already extracted by spectral_action / spectral sums. It is the framework's only post-observable projection.

**Track C seed:** Investigate whether other regularization-flow operations (Borel resummation, Mellin-Barnes contour rotation, Padé approximation, asymptotic-series re-summation) deserve their own projection rows in Paper 34. If so, the "flow tier" is a sector of the projection family and not a singleton. Currently drake_swainson is the only flow-tier projection in the dictionary. This is a candidate gap.

#### C-4. Rest-mass commutes with everything but is `D`-deforming

**Asymmetry kind:** Innocent at metadata level, surgical at spectral-triple level.

Rest-mass #14 is one of the simplest projections in the dictionary: it commutes with all 13 non-drake projections, introduces no transcendental (Paper 35 KG-1 verified ring-preserving for $m^2 \in \mathbb{Q}$), and adds a single variable + dimension (m, mass). At the TX-A signature-arithmetic level it is "innocent."

But at the spectral-triple level it modifies $D$: $D \mapsto D + m\cdot\mathbb{1}$ for the Klein-Gordon case, and a more substantive deformation for the Dirac case. *Any* modification of $D$ is, by Paper 31's partition, structurally Coulomb-specific (it depends on the choice of Dirac operator). Yet rest-mass modification is universal across central potentials — every central potential admits an additive mass shift. The asymmetry is: *the dictionary correctly tags rest-mass as universal at the V/D/T level, but the spectral-triple framing classifies it as `D`-touching, and Paper 31's partition would naively classify `D`-touching as Coulomb-specific.*

The resolution is that *additive scalar deformations of $D$* are an exception to Paper 31's partition: they are $D$-modifications that commute with everything in $\mathcal{A}$ and therefore live in the centre $Z(\mathcal{A}) = \mathbb{C}\cdot\mathbb{1}$. The centre is the place where universal-sector content can hide as a nominal $D$-modification.

**Track C seed:** Identify and characterize the centre of the partition. Is rest-mass the *only* projection that lives in $Z(\mathcal{A})$? Or are there others — perhaps a global-phase projection on $\mathcal{H}$, a Wick-rotation projection that affects only the time-direction sign, a coupling-constant rescaling — that share rest-mass's "trivially universal $D$-touch" status? The dictionary may be missing an entire sector of "central deformations."

#### C-5. Vector-photon is the unique pure-transcendental witness

**Asymmetry kind:** Axiomatic singularity.

Vector-photon #11 is the unique projection in the dictionary with axis pattern (V = –, D = –, T = +). It introduces no variable, no dimension, only a transcendental — $1/(4\pi)$ per loop. This is the projection that proves the transcendental axis is independent of V and D in the TX-A theorem.

But on the spectral-triple side, vector-photon is `H`-modifying (it lifts the photon 1-cochain to a vector harmonic with $(q, m_q)$ labels), which means it does change $\mathcal{H}_{\mathrm{GV}}$. So at the spectral-triple level it is structurally non-trivial; at the TX-A axis level it is "pure π."

This asymmetry is *informative*: it says that "pure transcendental" projections are exactly the ones that modify the *Hilbert space's basis* (in this case, the photon mode basis) without changing the *physical content* (no new variable, no new dimension). The transcendental enters through the *measure on the new basis* (the $S^2$ Weyl factor for vector harmonics). This is a structural prediction: any future "pure transcendental" projection should be a basis-refinement of $\mathcal{H}$ that introduces a measure-induced transcendental factor.

**Track C seed:** Search for other basis refinements of $\mathcal{H}_{\mathrm{GV}}$ that fit the (V, D, T) = (–, –, +) pattern. Candidates: (a) vector-photon analogs on the gauge sector (lift Wilson 1-cochain to vector basis; would inject $1/(4\pi)$ per gauge loop separately from QED); (b) tensor-photon promotion (rank-2 photon, would inject $1/(4\pi)^2$ or similar from $S^2$ tensor harmonics); (c) graviton promotion (would inject Gauss–Bonnet-like topological invariants of $S^2$). Each is a candidate new dictionary row in the (–, –, +) bucket.

#### C-6. The grading $\gamma$ has no projection that engages it

**Asymmetry kind:** Connes-axiom slot with no resident projection.

The GeoVac triple is correctly *odd* (no $\gamma$) because $S^3$ is odd-dimensional. But this means the dictionary's projections never interact with $\gamma$: there is no projection whose effect is to add a $\mathbb{Z}_2$-grading to the triple, and no projection's analysis gets blocked by the absence of $\gamma$.

The asymmetry: in the Connes–Chamseddine Standard Model spectral triple, $\gamma$ is *essential* — it implements chirality and distinguishes left- vs right-handed particles. If GeoVac were to acquire a chirality sector (e.g., by considering $S^3 \times S^1$ with the $S^1$ as Euclidean time, where the temporal compactification produces a Spin(4) bundle on the four-manifold $S^3 \times S^1$, which is even-dimensional and *does* carry $\gamma$), then a new projection — call it "Wick rotation" or "temporal-chirality" — would engage the grading slot.

This is precisely the structural location of observation/temporal-window #15. Adding the periodic time direction *does* even the triple (the ambient four-manifold is Spin(4)). But Paper 35 does not analyze the resulting four-manifold spectral triple's $\gamma$ structure — it stays on the spatial $S^3$ and treats time as an additive boundary.

**Track C seed:** Build the four-manifold spectral triple $(\mathcal{A}_{S^3 \times S^1}, \mathcal{H}_{S^3 \times S^1}, D_{S^3 \times S^1})$ with Spin(4) structure and verify that observation/temporal-window #15 acts as the projection that *promotes* the odd $S^3$ triple to the even $S^3 \times S^1$ triple. If so, observation/temporal-window is the dictionary's first $\gamma$-engaging projection, and Paper 35's "calibration-2π" tier acquires a structural reading: *the 2π in Stefan–Boltzmann is the $\gamma$-signature of the temporal compactification*. This is a falsifiable prediction at the spectral-action level.

#### C-7. K = π(B + F − Δ) crosses three sectors of the same triple

**Asymmetry kind:** Cross-sector composition with no common generator.

Paper 32 §6 reads K = π(B + F − Δ) as a sum across three structurally distinct sectors of the *same* GeoVac triple:
- B = finite Casimir trace at $m = 3$ (scalar $\mathcal{A}$ sector, finite trace).
- F = $\zeta_R(2) = \zeta_{\Delta_{S^3}}(2)$ (scalar $\mathcal{A}$ sector, infinite spectral zeta).
- $\Delta^{-1} = g_3^D = 40$ (spinor $\mathcal{H}$ sector, single-level Dirac degeneracy).

In Paper 34's projection-language reading, these are three different chains:
- B: Fock #1 alone.
- F: Fock #1 + spectral action #6 (or its Dirichlet sub-projection).
- $\Delta$: Fock #1 + spinor lift #7.

The K combination requires three different projection chains to evaluate the same $\alpha^{-1}$ — *and the chains share Fock #1 but differ in their second projection*. This is the dictionary-level reading of WH5 ("α is a projection constant, not a derivable number"): the three pieces live in three different sectors of the same triple under three different chains, and their combination has no common-generator chain.

Paper 34's depth-prediction (§VI) says a three-projection match should carry ~3% error; K matches at $8.8 \times 10^{-8}$. The dictionary therefore predicts that K is an under-counted chain — that the three chains share *more structure* than the dictionary currently records. Possible mechanisms: (i) Hopf #2 enters all three chains implicitly through the Vol(S²)/4 factor; (ii) the three chains share an $\mathcal{A}$-level invariant that the dictionary misses.

**Track C seed:** Re-derive the K-decomposition with explicit projection-chain accounting and identify whether the three chains share a fourth implicit projection (most likely Hopf #2). If so, K is a *single-projection*-depth coincidence in disguise (the shared Hopf base), and the depth-prediction is consistent. If not, K is a genuine anomaly under the depth-prediction — which would be the dictionary's first internal inconsistency, candidate for falsification.

#### C-8. Sturmian #5 is the only projection that adds zero

**Asymmetry kind:** Idempotent with zero V, zero D, zero T.

Sturmian #5 has axis pattern (–, –, –). It introduces no variable, no dimension, no transcendental. It is also idempotent. Its projection-taxonomy reading is "this is a relabelling of the same Fock graph at the Sturmian exponent $\lambda = Z/n$."

But on the spectral-triple side, the Sturmian relabelling is *the construction* by which the bound-state Fock graph becomes a complete Sturmian basis spanning bound states + discretized continuum. Without Sturmian, the Bethe-log spectral sums of Paper 36 are formally divergent. Sturmian is essential to the bound-state QED arc.

The asymmetry: a dictionary projection that adds *nothing* (V = D = T = –) is doing essential work at the spectral-triple level. The reason is that Sturmian is a *change of $\mathcal{H}$-basis* that preserves the algebra, the spectrum, and the ring — but rearranges what "summing over modes" means. Continuous integrations on the Sturmian-relabelled basis can converge where they diverge on the bare Fock basis.

**Track C seed:** Audit which other dictionary projections do "essential work" without registering on the V/D/T axes. Wigner 3j #8 is another candidate — it preserves all axes and yet implements the entire Gaunt selection rule structure. The "(–, –, –) but essential" cells are the dictionary's *invisible-but-load-bearing* projections, and they form a special class deserving its own tag.

#### C-9. Two energy-bearing projections (Fock + Bargmann) define different triples

**Asymmetry kind:** Same column-2 metadata, structurally distinct outputs.

Fock conformal #1 and Bargmann–Segal #3 both modify $(\mathcal{A}, \mathcal{H}, D)$ — both literally *define the triple* on which everything else operates. Both add an energy variable + dimension. At the dictionary level they look symmetric.

But they produce *distinct, non-isomorphic spectral triples* (Paper 31 Table 2; Paper 24): Coulomb $S^3$ vs HO $S^5$, second-order Riemannian vs first-order complex Euler, calibration-π vs π-free. This is *the* central asymmetry of the framework — and the dictionary's three-axis tagging is too coarse to see it. The TX-A composition table flags `(fock_conformal, bargmann_segal) = commute`, which is technically correct (both add an energy variable + dimension, so signature union is fine) but structurally misleading: composing the two would require working on two different sphere-quotient spectral triples simultaneously, which is structurally not what the framework does.

**Track C seed:** This is the dictionary's coarsest blind spot. Refine the column-2 sector code to distinguish "defines a new triple" from "modifies the existing triple" and re-run the composition table. Expectation: Fock and Bargmann–Segal will be flagged as *not composable* because they construct distinct triples; the variable_conflict will be redirected to the right place.

### 2.d Predictions

**P-1. Two new dictionary entries.** The "central deformations" sector (§2.c.4) and the "essential-but-axis-invisible" tag (§2.c.8) are missing structure. Specifically:
- A *global phase / Wick rotation* projection ((V, D, T) = (–, –, –) but $D$-deforming on the centre). This would join rest-mass #14 in the universal $D$-touching cell.
- A *temporal-chirality* projection that explicitly engages $\gamma$ on the four-manifold $S^3 \times S^1$ (§2.c.6).

**P-2. Operator-system inequivalence implies variable-conflict, but not converse.** The variable_conflict cells in the composition table (hopf+spinor, hopf+spectral, spinor+spectral) are a *necessary* condition for operator-system non-commutation. WH1 R3.2 verified this for hopf+spinor. The other two pairs (hopf+spectral, spinor+spectral) should likewise produce operator-system inequivalences when their respective truncated triples are computed. Specifically: at the operator-system level, $P_{n_{\max}} (\mathrm{Hopf} \circ \mathrm{spectral}) P_{n_{\max}}$ should differ from $P_{n_{\max}} (\mathrm{spectral} \circ \mathrm{Hopf}) P_{n_{\max}}$. This is testable via the same SDP framework as WH1 R3.

**P-3. Paper 35's temporal-window principle hits a triple-framing wall at $S^3 \times S^1$.** When time is compactified, the four-manifold becomes even-dimensional and the spectral triple acquires a natural $\gamma$. The four-manifold triple is structurally distinct from the spatial triple — it has different Connes axiom signs (KO-dim 4 instead of 3) and different real-structure conventions. Paper 35 does not analyze the four-manifold triple; the framing wall is at the moment temporal compactification is performed. The prediction: any spectral-action computation on $S^3 \times S^1_\beta$ should *not* match spectral-action on $S^3$ alone with a thermal trace tacked on; the two should differ at the spectral-action coefficient level. This is testable.

**P-4. The K-chain shares Hopf #2 implicitly.** Per §2.c.7, the three chains for B, F, Δ should be re-tagged to include Hopf #2 in all three. After re-tagging, K would be a *unified-three-projection-chain* match (Fock + Hopf + something) at $8.8 \times 10^{-8}$, which still violates the depth-prediction by 8 orders of magnitude — so this prediction does *not* save the depth-prediction, but it does identify where the under-counting is. The genuine question is whether the chains share *more than* Hopf #2.

**P-5. The (–, –, +) bucket should grow.** Vector-photon #11 is currently the unique pure-transcendental witness. The dictionary predicts that any *basis refinement of $\mathcal{H}$ that introduces a measure-induced transcendental factor* fits this bucket. Two specific candidates: tensor-photon (rank-2 photon) and graviton (rank-2 traceless photon on the Hopf base), each potentially injecting $1/(4\pi)^2$ from $S^2$ tensor harmonics. If these are added, the (–, –, +) bucket has multiple witnesses and the transcendental-axis independence becomes a tier rather than a singleton.

**P-6. KO-dim transition under temporal compactification.** Paper 32 audits the GeoVac triple at KO-dim 3. Adding observation/temporal-window #15 should transition the four-manifold triple to KO-dim 4 (mod 8). The prediction: the *signs* in the real-structure axiom flip from $(\epsilon, \epsilon') = (-1, +1)$ at KO-dim 3 to whatever the KO-dim 4 sign-table prescribes (Connes 1995 sign table: $(\epsilon, \epsilon', \epsilon'') = (-1, +1, +1)$ at KO-dim 4). Paper 35 does not check this. If verified, observation/temporal-window is the framework's first projection that *changes the KO-dimension* of the ambient triple.

---

## 3. Summary

**Dictionary completeness:** 15 rows tagged across 7 columns (sector, axioms, U/C, transcendental class, V/D/T axes, commutation flags). Every Paper 34 projection has a place in the spectral-triple framing; no projection is left untyped.

**Hypothesis test (§2.a):** The "all `A` is universal, all `D` is Coulomb-specific" hypothesis holds with three explicit caveats: (i) `*` projections (stereographic) are mixed; (ii) rest-mass is `D`-touching but universal because it lies in the centre of $\mathcal{A}$; (iii) observation/temporal-window is `B` (boundary), structurally outside the spectral-triple data and outside Paper 31's partition.

**Gaps (§2.b):** Order-one is the unique substantively-engaged Connes axiom; $\gamma$ is absent and untouched (correctly, given $S^3$ is odd); the universal sector is exactly the no-π sector (consistent with Paper 31).

**Asymmetries / Track C seeds (§2.c):** **9 distinct asymmetries flagged.** Roughly: one internal TX-A inconsistency (memo claim not in JSON, C-1); two metadata vs spectral-triple mismatches (variable_conflict ≠ operator-system non-commutation, C-2; "(–, –, –) but essential" projections, C-8); two singularities (drake_swainson is the unique flow-tier projection, C-3; vector-photon is the unique pure-T witness, C-5); one centre-of-algebra exception (rest-mass, C-4); one un-engaged Connes-axiom slot ($\gamma$, C-6); one cross-sector compositional puzzle (K-chain, C-7); one structurally-distinct-but-metadata-equivalent pair (Fock vs Bargmann–Segal, C-9).

**Predictions (§2.d):** Six concrete predictions, the most testable being P-3 (the four-manifold triple at $S^3 \times S^1_\beta$ has different spectral-action coefficients than $S^3$ + thermal trace), P-6 (KO-dim transitions under temporal compactification), and P-2 (the other two variable_conflict pairs should produce operator-system non-commutation at the SDP level, paralleling WH1 R3.2).

---

## 4. Report-back

- **Path to memo:** `debug/track_ts_b_dictionary_memo.md`
- **Asymmetries flagged in §2.c:** 9 (numbered C-1 through C-9, each with explicit Track C investigation seed)
- **Most surprising finding (one sentence):** The TX-A composition table's actual JSON data shows that the *only* non-commuting partner for any projection is drake_swainson — the spectacular UV/IR ordering anomaly the TX-A memo names as the "headline entry" (spectral_action ↔ observation_window mutual one-way) is *not in the data*, which flags both directions as commute, making C-1 (the memo–data discrepancy) the most likely single Track C investigation: either the layer grading is too coarse or the memo overstates what signature arithmetic can capture.
