# Follow-on register

**Last updated:** 2026-06-04 (post verification-cleanup workflow `wf_39fc232e-906`; A9–A14 added for cross-arc drift)

This file is the durable list of carved-out follow-ons from recent sprints that
the PI has not yet authorized for execution. New items are added at the top of
each section; items get crossed out (or removed) when completed. Each entry
should be self-contained enough that a fresh PM session can pick it up cold,
including the originating sprint reference for context.

---

## A. Substantive research sprints (PI decision required)

### ~~A9.~~ ~~Paper 34 — 22 of 28 named projections still without numerical tests~~ — **CLOSED 2026-06-04 (A9 batches 1+2+3)**
All 28 named Paper 34 §III projections now have numerical spot-checks
across four files (`tests/test_paper34_projection_spot_checks.py` for 6,
`..._batch1.py` for 8 load-bearing, `..._batch2.py` for 7
gauge/symmetry/separation, `..._batch3.py` for 7 remaining). 108 tests
total, all passing. §13.4a equation-verification coverage on Paper 34 is
complete. **Side discovery:** Paper 34 §III.16 had a transcription error
in the $R^0_\text{BP}(1s,2s;1s,2s)$ closed form (paper said
$-4\log 2 - 19/9 + 9\log 3/2 \approx 0.060$, correct value is the pure
rational $4/81 \approx 0.0494$ matching the production module's existing
regression test). Paper 34 corrected; autopsy memo
`debug/ps_1s2s_autopsy_track4_memo.md` line 81 corrected with provenance
note. Also surfaced CLAUDE.md §7 drift: claimed `geovac/constants.py` does
not exist (KAPPA_SCALAR lives at `geovac/graph_qed_propagator.py:103`);
filed as C4 below.

### A10. Paper 40 §sec:open trio — log-power ansatz, rank-3 SU(4), compact symmetric spaces G/K
**Source:** Verification cleanup math.OA drift scan (2026-06-04); Paper 40
§sec:open.
**(i) Log-power ansatz test (sprint-scale, 1–2 weeks):** is γ_n ~
(4/π)·log(n)/n single-log or double-log at rank ≥ 2? Clean discriminator
needs Λ² ≳ 10⁵. Analytic Plancherel × Vandermonde already shows
rank-invariance (Paper 40 Thm L2); a sharper Stein–Weiss derivation at
rank 2 would settle the question.
**(ii) Rank-3 SU(4) numerical extraction (sprint-scale, 1–2 weeks):**
current 29.3 % error; pre-asymptotic structure dominates at rank 3.
Larger n_max sweep would verify 4/π constant directly.
**(iii) Compact symmetric spaces G/K (medium, 3–6 months):** extend from
compact Lie groups to S^n, CP^n, Grassmannians, flag varieties. Zonal
spherical-function machinery expected to lift; 4/π presumed to persist.

### A11. Paper 42 §10 O3 + O4 — Pythagorean formal proof, higher-rank modular Hamiltonian
**Source:** Verification cleanup math.OA drift scan (2026-06-04); Paper 42 §10.
**(O3) Pythagorean diagonal / off-diagonal formal proof (sprint-scale, 2–3 weeks):**
the H_local diagonal + D_W^L off-diagonal decomposition under Δn = ±1 sectoring
is sketched (rem:pythagorean_underlies_collapse, Sprint L2-F.1) but not formally
proved at general n_max.
**(O4) Higher-rank modular Hamiltonian (sprint-scale, 2–3 weeks):** verify
integer-spectrum property of BW-α generator analog on SU(n) for n ≥ 3.

### A12. Paper 50 §8 — Maxwell, squashed-S³, quantitative propinquity rate
**Source:** Verification cleanup math.OA drift scan (2026-06-04); Paper 50 §8
Catalogue and open targets.
**(i) Maxwell partition function p = 1 (sprint-scale, 2–4 weeks):** apply
F-theorem coefficient + Pythagorean M2/M3 decomposition machinery to Maxwell
(vector p = 1) and graviton (p = 2) on S³/S⁵. Closed forms known in literature.
**(ii) Squashed-S³ generalization (sprint-scale, 1–3 weeks):** does the
deformed half-integer Hurwitz machinery still reduce to the MR-B
Jacobi-θ_2 closed form, or does the deformation introduce a non-M2/M3
transcendental?
**(iii) Quantitative propinquity rate (medium, 3–6 weeks):** prove
|F_{n_max} − F^{KPS}| ≤ C·γ_{n_max} via Mellin-domain Berezin reconstruction.

### A13. Paper 53 — Q1 / Q2 / Q3 / Q4 cluster (disk-with-cone open questions)
**Source:** Verification cleanup gravity-arc drift scan (2026-06-04);
Paper 53 §Open.
**(Q1) Sharp self-contained propinquity constant C (medium, 3–6 weeks):**
closed-form constant in the disk-with-cone full-plane propinquity bound.
**(Q2) Sharp interior Berezin reconstruction rate (medium, 3–6 weeks):**
plus closed-form Cesàro-order dependence.
**(Q3) Boundary-adapted operator system (medium, 3–6 weeks):** for
finite-R reconstruction of boundary observables.
**(Q4) Higher-dimensional + warped carriers (multi-month, 3–6 months):**
extend to cigar D²_α × S² as a single carrier.

### A14. Gravity arc — G4-6a Q1–Q3 dichotomy, Möbius v3.19.0 audit, L6 write-out, L6 sharp rate
**Source:** Verification cleanup gravity-arc drift scan (2026-06-04); items
flagged in `papers/group5_qed_gauge/paper_51_gravity_arc.tex` §sec:open and in
L6 closure memo lines 135–138.
**(i) G4-6a Q1–Q3 dichotomy decision (sprint-scale, 1–2 weeks):** resolve
whether A-coefficient extraction is (Q1) a genuine multi-month substrate
refinement, (Q2) a compute-budget limit, or (Q3) an observable mismatch.
Currently unresolved after the N_φ-axis clean negative.
**(ii) Möbius v3.19.0 sign-convention audit (sprint-scale, 1–2 weeks):**
substrate slope +0.052 (Track α'' thread 9) vs v3.19.0 reported −0.0562.
Note: the v3.19.0 Fursaev–Solodukhin mechanism attribution was already
retired in §3 dead-ends (B4, v3.24.0); this is the remaining sign-convention
audit only.
**(iii) L6 B3-standalone L1'–L5 write-out (medium, 3–6 weeks):** transport
the framework-grade proof sketch from the L6 closure memo into a
self-contained document (Paper 53 §appendix-B or standalone memo). No new
analytical content; bookkeeping only.
**(iv) L6 sharp lattice rate + B1 constant (medium, 3–6 weeks):** current
empirical scaling is O(a^{1.1}); the expected exact O(a) from Lichnerowicz
scalar curvature is not yet recovered. Resolve and extract the sharp B1
constant.

### ~~A1.~~ ~~Inhomogeneous Robertson–Walker extension~~ — **STRATIFIED 2026-06-03 (Sprint A1 scoping)**
A1 splits cleanly into three sub-cases. (i) Constant warp $a(t) \equiv a_0$:
trivial corollary of Paper 51 G1 `thm:zeta_unit_neg_k`; closed by drop-in
paragraph in Paper 55 §7.2. (ii) Generic continuum $a(t)$ with non-trivial
$\varepsilon_i = a^{(i)}(t) \ne 0$: NO-GO at discrete-substrate level; tied
to Paper 51 G4-3 discrete-warped-substrate multi-month program. (iii)
Thermal-time compactification $S^3 \times S^1_\beta$ (Paper 35 Matsubara
substrate, the natural GeoVac time-direction surrogate, NOT literally an
$a(t)$): sprint-scale 2-3 weeks, named A1-Matsubara (new item below).
Memo: `debug/sprint_a1_scoping_memo.md`. Paper 55 §7.2 updated with the
three-sub-case stratification.

### ~~A1-Matsubara.~~ ~~Thermal-time compactification of M2 pure-Tate refinement~~ — **CLOSED-NEGATIVE 2026-06-03 (Sprint A1-Matsubara)**
Sub-agent returned NEGATIVE-no-structural-engagement. Temporal compactification
multiplies SD coefficients by Tate-weight-0 rational factor $\beta$ and changes
no other structural feature. Explicit closed forms (sympy-verified bit-exact):
scalar $a_k^{\mathrm{4D,scalar}}(\beta) = 2\pi^2\beta/k!$ for all $k \ge 0$;
Dirac $a_0 = 4\pi^2\beta$, $a_1 = -2\pi^2\beta$, $a_k = 0$ for $k \ge 2$
(two-term exactness inherited verbatim). Mechanism: $m \ne 0$ windings of
$K_{S^1_\beta}(t)$ are exponentially suppressed in small-$t$ asymptotic; the
Matsubara mode sum lives in the orthogonal $t \to \infty$ regime reached via
Jacobi $\vartheta_3$ inversion. Two regimes don't exchange transcendental
content at SD level. **No conflict with Paper 35 Stefan-Boltzmann** $\pi^2/90$
— that's the M1×M2 Mellin-at-integer-$s$ extraction, structurally distinct
observable. Paper 55 §7.2 thermal-time bullet upgraded to closure-in-the-negative
+ new Proposition `prop:thermal_pure_tate`. Memo:
`debug/sprint_a1_matsubara_memo.md`.

### ~~A2.~~ ~~S^5 Bargmann–Segal mixed-Tate analog~~ — **CLOSED 2026-06-03 (Sprint A2)**
Sub-agent dispatch returned POSITIVE with structural refinement: Dirac D² on S⁵ is
three-term exact (a₀ = 4π³, a₁ = -20π³/3, a₂ = 3π³); scalar Laplacian gives infinite
closed form a_k^Δ = (6-k)·4^(k-1)·2/(3·k!)·π³ with single mid-series zero at k=6.
Both sit in single Tate-weight-3 slice π³·Q. **Dimension-parity sharpening** (new
structural finding): S³ (even-dim Vol=2π²) → even-weight ⊕π^(2k)·Q; S⁵ (odd-dim
Vol=π³) → single odd-weight slice π³·Q. Closed forms in
`debug/sprint_a2_s5_sd_coefficients.py`; memo at
`debug/sprint_a2_s5_mixed_tate_memo.md`. Paper 55 §7.3 updated with Proposition
`prop:s5_pure_tate` + dimension-parity remark; new Q3' subsection records
M3-on-S⁵ as new sprint-scale open question (2-3 weeks).

### ~~A3.~~ ~~M3 cyclotomic-mixed-Tate verdict at level 4 over Z[i]~~ — **CLOSED 2026-06-03 (Sprint M3 + Paper 55 §5)**
Sprint M3 Cyclotomic Mixed-Tate verified M3 ⊂ MT(Z[i, 1/2]) at level ≤ 4 via
direct application of Deligne 2010 + Glanois 2015. Vertex-parity-as-χ_{-4} identified
as the literal Deligne-Glanois Galois descent from level 4 to level 2 on the
framework's natural QED observables. Closed by Paper 55 §5 fill (Theorem 5.1
`thm:m3_cyclotomic_mixed_tate` + Proposition 5.3 `prop:vertex_parity_descent`).
Memo: `debug/sprint_m3_cyclotomic_mixed_tate_memo.md`.

### ~~A6.~~ ~~M3-on-S⁵ cyclotomic-mixed-Tate refinement~~ — **CLOSED 2026-06-03 (Sprint A6)**
Sub-agent returned POSITIVE. Cyclotomic-mixed-Tate-at-level-≤4 transfers from
S³ to S⁵ verbatim. Explicit S⁵ Dirac Dirichlet:
$D^{(\sfive)}(s) = (1/3)Z(s-4) - (5/6)Z(s-2) + (3/16)Z(s)$ (three-term, vs S³'s
two-term). χ_{-4} identity:
$D_{\mathrm{even}}^{(\sfive)} - D_{\mathrm{odd}}^{(\sfive)} = (1/3)f_5(s-4) -
(5/6)f_5(s-2) + (3/16)f_5(s)$ with $f_5(s) := 2^s\beta(s) - 2^s + (2/3)^s$,
verified bit-exact at $s \in \{6,7,8,9,10\}$. **Level 4 sufficient** — S⁵ shifts
$5/4, 7/4$ reduce to level-4 shifts $1/4, 3/4$ via Hurwitz $\mathbb{Q}$-rational
corrections. Paper 55 §5.5 new subsection `subsec:m3_s5` with two Theorems +
Corollary; §7.4 updated to mark structural closure. Memo:
`debug/sprint_a6_m3_s5_memo.md`. Driver: `debug/sprint_a6_m3_s5_derivation.py`.

### A6-followon. S_min^(S⁵) PSLQ irreducibility (new from A6)
**Source:** Sprint A6 closure (Paper 55 §5.5 paragraph "Specific motivic
identification").
**Estimated effort:** 2-3 weeks (mirrors S³ sprint of same effort).
**Question:** 200-dps PSLQ irreducibility of $S_{\min}^{(\sfive)} := \sum_k
T_5(k)^2 \approx 0.0399612091165\ldots$ against the Paper 28 §$S_{\min}$
basis. Natural prediction: $S_{\min}^{(\sfive)} \in \MT(\Z[1/2])$ at depth 2.

### ~~A7.~~ ~~M2-M3 cyclotomic coincidence at level 4~~ — **HALF-STRUCTURAL CLOSED 2026-06-03 (Sprint A7)**
Sub-agent returned HALF-STRUCTURAL. M2's $\Q(i)$ is the Witt-splitting field of
the standard Euclidean quadratic form $Q_{1,2n}$ (Gal acts trivially on M2
outputs, pure-Tate descends to $\Q$). M3's $\Q(i) = \Q(\zeta_4)$ is the
cyclotomic conductor of $\chi_{-4}$ (Gal acts non-trivially via Glanois descent
$N=4 \to N=2$). Both trace to dimension-3 parity of S³ via INDEPENDENT
mechanisms — no shared Tannakian symmetry. Sharpened (not promoted) by A8
factored form $[Z_{1,2n}] = [\mathbb{P}^n] \cdot (1 + \mathbb{L}^n)$. Paper 55
§4 Remark `rem:m2_specialisation` extended with Distinction paragraph; new §7.6
Q5' subsection opened for deeper unification (recorded as multi-year research
project under Section A.M2-M3-unification below). Memo:
`debug/sprint_a7_m2_m3_cyclotomic_memo.md`.

### A.M2-M3-unification. M2/M3 motivic Galois unification (multi-year Marcolli–Tabuada construction of $\omega^{\mathrm{tri}}$)
**Source:** Sprint A7 (2026-06-03; memo
`debug/sprint_a7_m2_m3_cyclotomic_memo.md`). Six-sprint scoping arc closed
2026-06-03/04 (memos `debug/sprint_q5p_dim_sweep_memo.md`,
`debug/sprint_q5p_qsm_litread_memo.md`,
`debug/sprint_q5p_tannakian_obstruction_memo.md`,
`debug/sprint_q5p_deflation_test_memo.md`,
`debug/sprint_q5p_k_slot_tannakian_memo.md`,
`debug/sprint_q5p_greenfield_marcolli_transport_memo.md`).
**Estimated effort:** multi-year Marcolli–Tabuada-lineage construction;
no sprint-scale handle. Viable lineage strands A/C/D/E only (Connes–Marcolli
motivic Galois, arXiv:math/0409306 / Marcolli–Tabuada noncommutative motives,
arXiv:1110.2438 + arXiv:1112.5422 / Fathizadeh–Marcolli periods, arXiv:1611.01815
/ Deligne–Glanois cyclotomic descent, arXiv:math/0302267 + arXiv:1411.4947).
QSM / Bost–Connes / Greenfield–Marcolli–Teh lineage (strand B,
arXiv:1305.5492) RULED OUT on four structural grounds.
**Sharpened target:** construct the enriched fiber functor
$\omega^{\mathrm{tri}}: \mathrm{dg}(\Tcal) \to \mathrm{Vec}_\Q \otimes
\mathrm{IndexCat}(\{0, 1, 2\})$ recording the master Mellin engine
slot index $k$ on the pro-dg category $\varprojlim_{n_{\max}}
\mathrm{dg}(\Tcal_{n_{\max}})$ with Berezin reconstruction maps
(Paper 38 §L4), and check whether its associated motivic Galois group
reads $k$ as a Tannakian invariant. Period-ring level deflates
($M^{\mathrm{GV}} \subset \mathrm{MT}(\Z[i,1/2], 4)$ already
suffices); dg-category level is the precise multi-year target.
$k$-slot Tannakian status BORDERLINE (TANNAKIAN-INVISIBLE on ambient
$M^{\mathrm{GV}}$, TANNAKIAN-RELEVANT on the candidate enriched
$\omega^{\mathrm{tri}}$).

### ~~A8.~~ ~~F-M general-n Grothendieck-class closed form~~ — **CLOSED 2026-06-03 (Sprint A8)**
Sub-agent retrieved closed form from F-M §7.6 Theorem 7.6: three equivalent
forms. **Factored (structural reading):
$[Z_{1,2n}] = [\mathbb{P}^n] \cdot (1 + \mathbb{L}^n)$** — two maximal isotropic
subspaces over $\Q(i)$ as the source of the doubled middle Tate twist.
Palindrome + complement forms also given. Substituted into Sprint A5 GeoVac
Grothendieck-class equation:
$[\mathcal{V}_n^{\mathrm{GeoVac}}] = \mathbb{L}^{2n+3} - 3\mathbb{L}^{2n+2} +
2\mathbb{L}^{2n+1} - \mathbb{L}^{n+2} + 3\mathbb{L}^{n+1} - 2\mathbb{L}^n$.
Verified symbolically at $n=1,\ldots,5$. Matches F-M Thm 7.2 at $n=1$. Paper 55
§4 Remark `rem:m2_specialisation` expanded with explicit polynomial
(`eq:z_quadric_class` + `eq:geovac_grothendieck_class_explicit`). Memo:
`debug/sprint_a8_grothendieck_class_memo.md`.

### ~~A4.~~ ~~Inner-factor mixed-Tate extension~~ — **CLOSED-BY-PRIOR 2026-06-03 (Sprint A4 scoping)**
Three theorem-grade prior ingredients jointly settle A4 without new
computation: (i) Paper 18 §IV.6 `thm:eta_trivialization` (inner M3
vanishes); (ii) Paper 18 §IV.6 `thm:ac_factorization` (heat trace
factorises cleanly into outer × inner); (iii) Sprint Yukawa-PSLQ
(2026-06-03; 162-cell sweep, zero hits, Yukawas not in low-coefficient
pure-Tate). Combined: AC-extension SD coefficients factor as (outer
pure-Tate $\bigoplus_k \pi^{2k}\mathbb{Q}$) × (inner Dirichlet
$\mathbb{Q}[y_i^{-2s}]$); inner factor generically NOT mixed-Tate over
$\mathbb{Q}$ because Yukawas are Class 1 calibration. Inner-factor
content lives in Paper 18 §IV.6 sixth tier (parameter-tied Dirichlet),
categorically disjoint from M1/M2/M3. **No new sub-mechanism required.**
Memo: `debug/sprint_a4_scoping_memo.md`. Paper 55 §7.4 updated with the
factorisation + sixth-tier escape route.

### ~~A5.~~ ~~Fathizadeh–Marcolli full paper read~~ — **CLOSED 2026-06-03 (Sprint A5)**
Sub-agent dispatch returned POSITIVE with one correction. F-M Thm 6.2 gives the
explicit affine complement $\mathcal{V}_n^{F\text{-}M}(\lambda) = \mathbb{A}^{2n+3}
\setminus (C_{Z_{\lambda,2n}} \cup H_0 \cup H_1)$ with quadric
$Q_{\lambda,2n} = u_1^2 + \lambda^{-2}(u_2^2+u_3^2+u_4^2) + u_5^2 + \ldots$;
GeoVac static specialisation is the closed immersion $\{1\} \times \{0\}^{2n}$.
**Correction to v3.45.3 expected reading:** "trivial complement (a point)" is half-right
— the $(\lambda,\varepsilon)$-fiber collapses, but the $u$-coordinate complement
$\mathcal{V}_n^{GeoVac}$ remains a genuine $(2n+3)$-dim mixed-Tate variety.
Grothendieck-class formula computed (eq:geovac_grothendieck_class). **Side observation
flagged as open question:** the M2 quadric ring $\mathbb{Q}(\sqrt{-1}) = \mathbb{Q}(\zeta_4)$
matches the M3 vertex-parity Glanois level-4 ring — possible deeper M2-M3 cyclotomic
coupling at level 4. Memo at `debug/sprint_a5_fm_full_read_memo.md`. Paper 55 §4
proof sketch + new Remark `rem:m2_specialisation` applied.

---

## B. Collaboration outreach (PI decision required; user has flagged
   as a separate later comprehensive push)

### B1. Garofalo / Urbach group (Bonn / DESY)
**Source:** Phase 4 long-tail lit survey (v3.45.2, 2026-06-03; memo
`debug/lit_survey_phase4_longtail_memo.md`); now cited in Paper 30
(Garofalo et al. arXiv:2311.15926 + Jakobs et al. arXiv:2503.03397).
**Why:** they have built SU(2) Hamiltonian lattice gauge theory using
$S_3$ symmetric-group partitionings — same Wilson Hamiltonian on the
same group manifold as Paper 30, from a quantum-computing-digitization
motivation. Two communities (math-phys QFT-on-graphs vs lattice-QCD-for-QC)
not aware of each other's work.

### B2. Brian Hall (Notre Dame) on coherent states for $S^3$ / Bargmann–Segal
**Source:** Phase 4 long-tail lit survey (v3.45.2); now cited in Paper 24
(Higgs–Pickrell arXiv:2503.23549; Hall's coherent-states-on-spheres
papers also cited).
**Why:** Hall's group's open question for $d \ge 3$ on SO(d)-invariant HO
on d-spheres is exactly Paper 24's $d=5$ case. Natural conversation partner.

### B3. Stetcu / Sarma-Stevenson nuclear VQE
**Source:** Phase 4 long-tail lit survey (v3.45.2); now cited in Paper 23
(Sarma–Stevenson arXiv:2510.02124, Pillai et al. arXiv:2509.08948).
**Why:** direct head-to-head competitors for Paper 23's deuteron 16-qubit /
592-Pauli headline.

### B4. Hekkelman / van Suijlekom / Latrémolière / Farsi cluster
**Source:** Phase 2 connection-finding lit survey (v3.45.2; memo
`debug/lit_survey_phase2_connections_memo.md`).
**Why:** 4 directly-adjacent papers in 2024 alone (Hekkelman-McDonald
arXiv:2412.00628, van Suijlekom arXiv:2409.02773, Farsi-Latrémolière
arXiv:2404.00240, Farsi-Latrémolière-Packer Adv. Math. 437). GeoVac math.OA
arc (Papers 38–49) is in citation distance.

---

## C. Mechanical / technical follow-ons (PM-actionable, low-substantive)

### ~~C1.~~ ~~Paper 14 figure regeneration~~ — **CLOSED 2026-06-04**
Created `benchmarks/qubit_encoding/generate_paper14_figures.py` (figure
generator dedicated to Paper 14, using canonical Table I/II data so the
figures match the paper exactly). Generates `pauli_scaling.png` and
`eri_density.png` directly into `papers/group4_quantum_computing/paper_14_figures/`.
Extracted GeoVac scaling exponent 3.147 matches Paper 14's claimed
$\alpha = 3.15$. Paper 14 `\fbox` placeholders replaced with
`\includegraphics`. Three-pass clean compile, 21 pages.

### C3. Paper 14 missing bibitems + section labels (NEW from C1 compile audit)
**Source:** Closing C1 surfaced 8 cite keys used in Paper 14 but missing
from the inline `\thebibliography` block: `reiher2017`, `BJL`,
`ChildsBerry`, `Szmytkowski2007`, `MartinezYRomero2004`, `Dyall`,
`BreitPauli`, `Sunaga2025`. Plus 3 missing section labels (cross-refs to
sections that don't exist in the current Paper 14): `sec:pk_partitioning`,
`sec:frozen_core`, `sec:paper22_spinor`. The missing `sec:*` labels
suggest text was migrated out without updating refs; redirect to the
correct papers (likely Paper 22 for sec:paper22_spinor; Paper 17 / 19 for
PK and frozen core) or remove the refs.
**Estimated effort:** 2-4 hours. Requires looking up each bibitem in
existing project corpus (most are likely already cited in other papers).

### C4. CLAUDE.md §7 — geovac.constants module doesn't exist (NEW from A9 batch 1)
**Source:** A9 batch 1 test ModuleNotFoundError. CLAUDE.md §7 Code
Architecture table lists `Physical constants | geovac/constants.py |
HBAR, C, ALPHA, etc.` but the file does not exist. `KAPPA_SCALAR` lives
at `geovac/graph_qed_propagator.py:103`.
**Decision needed:** (a) create `geovac/constants.py` as a true central
constants module and migrate `KAPPA_SCALAR` and others; (b) update
CLAUDE.md §7 to reflect reality (drop the row or point to the correct
module). Per CLAUDE.md §13.5 access control, PM may edit §7 directly
to update module paths.
**Estimated effort:** 15 minutes for (b); 1-2 hours for (a) plus
downstream import-path updates.

### ~~C2.~~ ~~Paper 14 §sec:ft_gaussian BLISS-THC comparison table~~ — **CLOSED 2026-06-04 (as paragraph not table)**
Added a contextualizing paragraph ("BLISS-THC at the production-scale
ceiling") to Paper 14 §sec:ft\_gaussian after the "Computed 1-norm
comparison" paragraph. Pulled specific numbers from Caesura et al. 2025
(arXiv:2501.06165) via PDF extraction: P450 (63e/58o, Q=116) BLISS-THC at
999 memory qubits, λ ≈ 133 Ha, 1.71×10⁹ Toffolis; prior-art THC at
λ = 388.9 Ha, 7.79×10⁹ Toffolis; ~233× total speedup with breakdown
(25× AV + 8× BLISS + 1.1× circuit). Framed honestly per CLAUDE.md §1.5:
the GeoVac C₂H₆ entry at Q=160 is the framework's library ceiling but
the metrics aren't directly comparable (Pauli/QWC for VQE vs Toffoli/λ
for FT qubitization-LCU). Paper 14 compiles 22 pages, zero undefined
references. A literal head-to-head table is structurally impossible
because the molecule sets and circuit-family targets don't overlap.

---

## D. Closed in v3.45.3 (no longer follow-ons)

- ~~Mixed-Tate period test sprint~~ → POSITIVE-with-pure-Tate-refinement;
  Paper 32 Corollary applied; Paper 18 forward pointer closed; Paper 28 T9
  footnote added; Paper 51 rem:cc_uncanny cross-ref added; CLAUDE.md WH2
  status updated.
- ~~Paper 14 dcolumn shorthand bug~~ → fixed (preamble macro added).
- ~~Paper 14 missing figure compile error~~ → fixed (placeholders;
  see C1 above for actual figure regeneration).
- ~~Paper 30 remark environment undefined~~ → fixed (preamble macro added).
- ~~Paper 41 marcolli_vs2014 orphan bibitem~~ → closed (added "Relation to
  Marcolli–vS gauge networks" paragraph at end of related-work section
  citing all three: Marcolli-vS 2014, Perez-Sanchez 2024, Perez-Sanchez 2025;
  fixes the orphan AND closes the Phase 2 recommendation to add
  Perez-Sanchez follow-ups across Papers 25/30/41 — 25 and 30 already had
  them, 41 didn't).

## E. Closed in v3.45.2 (no longer follow-ons)

- ~~Comprehensive 4-phase literature audit~~ → done; memos in
  `debug/lit_survey_phase{1,2,3,4}_*_memo.md`.
- ~~Citation verification gate~~ → done; 18% hallucination rate caught,
  6 mismatches corrected; table in `debug/citation_verification_table.md`.
- ~~21 verified citations added across 14 papers~~ → done.
- ~~5 framing edits (Papers 35, 36, 48, 50, 51)~~ → done.
- ~~Pre-existing Paper 51 `\nmax` macro bug~~ → fixed.
- ~~Pre-existing Paper 36 `\to` in `\text{}` bug~~ → fixed.

---

## Maintenance notes

- Add a row to the appropriate section when a new follow-on is carved out
  by a sprint that does not authorize execution.
- Move completed rows to section D / E / etc. with a strikethrough and
  brief note rather than deleting (preserves the audit trail).
- Reference the originating sprint memo for each entry so the next PM
  session can find context cold.
- This file is referenced from CLAUDE.md §2 most-recent-sprint entry; keep
  it in sync.
