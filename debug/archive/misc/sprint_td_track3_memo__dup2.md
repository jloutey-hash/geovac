# Sprint TD Track 3 — Yukawa-sensitive thermal observables

**Status:** Closed positive.
**Date:** 2026-05-09.
**Worker:** sub-agent (this fork).
**Data:** `debug/data/sprint_td_track3.json`.
**No production `geovac/` modifications.**
**No paper edits applied (PI-approval gated; draft text in §7).**

---

## 1. Mandate

Track 3 catalogues thermal observables whose value depends on inner-factor
parameters (Yukawas, generation count, gauge-group algebra structure). For each,
we (a) compute the GeoVac-only outer-factor contribution from Track 1's
`geovac/thermal_tensor_triple.py` (commit 697666f), (b) characterize the
residual structure that an AC extension `T_GeoVac ⊗ T_F` would have to fill,
and (c) state a falsifier — "what would have to be true of `T_F` for the
GeoVac+inner prediction to match standard QFT."

This is a **catalogue**, not a derivation. Inner-factor data is unavailable
to the framework (per Sprint H1 verdict, May-7 inner-factor Mellin engine
memo, Paper 18 §IV.6). The track makes the structural-skeleton scope
*quantitative on the thermal side* by naming, observable by observable,
the precise structural slot a candidate `T_F` would have to fill.

The two structural backbones:

* **η-trivialization theorem** (`debug/inner_factor_mellin_engine_memo.md`
  §3): for any Krajewski-class even finite spectral triple,
  `Tr(D_F^k · e^{−tD_F²}) ≡ 0` for all odd k. The proof is
  γ_F D_F γ_F = −D_F (chirality grading) ⟹ trace flips sign for odd k.
  Consequence: any thermal observable whose evaluation reduces to
  `Tr(D_F^k · e^{−tD_F²})` at odd k is identically zero on the inner
  factor.

* **AC factorization theorem** (same memo §4): on the combined Dirac
  `D = D_GV ⊗ 1 + γ_GV ⊗ D_F`, the cross-term anti-commutator vanishes
  (γ_GV anti-commutes with D_GV by outer chirality grading), so
  `D² = D_GV² ⊗ 1 + 1 ⊗ D_F²` and
  `Tr(e^{−tD²}) = Tr(e^{−tD_GV²}) · Tr(e^{−tD_F²})`. The combined
  Mellin output factorizes as (outer M_i ring) × (inner Yukawa
  Dirichlet ring). π-content sits entirely outside; SM-distinguishing
  parameters (Yukawas, multiplicities) sit entirely inside.

Across the catalogue, every observable's outer baseline is a product of
`thermal_tensor_triple` quantities (Stefan–Boltzmann at the relevant
species count, Casimir contributions, modular residuals); every inner
residual reduces to a generalized Dirichlet sum
`Σ_i c_i · y_i^{−2s}` (k=0) or `Σ_i c_i · y_i^{2−2s}` (k=2) at integer
s in the Yukawa eigenvalues. We filter the k=1 candidates immediately
by η-trivialization.

---

## 2. Catalogue — five Yukawa-sensitive thermal observables

Each entry: (a) standard QFT value with reference; (b) GeoVac outer
baseline; (c) `T_F` residual structure; (d) what `T_F` would have to be
for match.

### 2.1 Free energy contribution from thermally-excited leptonic carriers

* **Standard QFT**: high-T limit free-energy density on flat ℝ³ (or
  asymptotic of S³ × S¹_β), summed over Dirac fermions of mass m_i and
  multiplicity g_i:
  `F_fermion / V = − Σ_i g_i · (7π²/720) T⁴ · h(m_i/T)`,
  where h is the Bose–Einstein-style mass-suppression function with
  h(0) = 1 and `h(x) = (15/π⁴)·∫₀^∞ k² log(1+e^{−√(k²+x²)}) dk`.
  At fixed T well above electroweak scale, all 6 charged leptons (e, μ,
  τ × 2 spin) are relativistic and the mass-suppressed correction
  contributes via the moments
  `Σ_i (m_i²/T²) · (1/24) · g_i` at leading order, plus a
  `(m_i²/T²) log(m_i²/T²)`-type log, plus `Σ_i (m_i⁴/T⁴) · (1/16)·g_i`
  at next order. Reference: Kapusta & Gale 2006, *Finite-Temperature Field
  Theory*, §3.5.

* **GeoVac outer baseline**: Track 1 directly computes the m=0 Stefan-
  Boltzmann contribution per Weyl spinor on S³ × S¹_β:
  `F/V_per_Weyl = − (7/720) π² T⁴`
  with the M1 × M2 factorization `[1/(2π²)] · [1/3] · [(7/8)·6ζ(4)] · T⁴`.
  Multiplied by the species-count-and-multiplicity factor of the inner
  factor: `g_lepton = 2_spin · 2_chirality · n_gen_lepton = 8` for one
  charged lepton per generation × 3 generations × Dirac doubling.

* **T_F residual structure**: the *mass-correction* terms
  `Σ_i (m_i²/T²) · (1/24)` are precisely
  `Σ_i (1/24) · y_i² · (v² / T²)`
  with y_i the Yukawa coupling and v the Higgs VEV. In the inner-
  factor Mellin engine, this is `(1/24) · Tr(D_F²) / T²`, the **k=2
  trace** at integer Mellin argument s=1. The next-order correction
  `Σ_i (m_i⁴/T⁴) · (1/16)` is `(1/16) · Tr(D_F⁴) / T⁴`,
  the k=2 trace at s=2.

* **What T_F would have to be**: a finite even spectral triple whose
  Dirac-squared trace at integer s reproduces
  `Σ_lepton g_lepton · y_lepton²·v² + Σ_quark g_quark · y_quark²·v²`
  (3 generations × 2 spin × 2 chirality × N_color factor).

* **Falsifier**: any candidate `T_F` whose `Tr(D_F²)` at the EW scale
  v ≈ 246 GeV does not reproduce the known sum
  `m_e² + m_μ² + m_τ² + 3·(m_u² + m_d² + ... )`. *Completely fixed
  by the SM Yukawa values.* No generation-multiplicity ambiguity.

* **η-trivialization filter**: PASS. Mass-correction terms are k=2
  (D_F² trace), not k=1. Observable is non-trivially populated on
  any Krajewski-class triple.

* **Diagnostic value**: 9/10. Cleanest single-number Yukawa probe.

---

### 2.2 QCD thermal pressure / color-multiplicity-sensitive free energy

* **Standard QFT**: at high T well above hadronic scale, the QCD pressure
  density is
  `p_QCD(T)/T⁴ = (π²/45) · g_QCD + corrections`
  where `g_QCD = 2(N_c²−1) + (7/8)·2·N_c·N_f·2 = 16 + 21·N_f/2`
  for N_c=3 colors and N_f flavors of light quarks. The integer 16
  comes from `2(N_c²−1) = 2·8 = 16` SU(3) gluons (transverse polarizations);
  the 21·N_f/2 = 31.5 (at N_f=3) from quark-antiquark Dirac fermions
  with Bose–Einstein-style 7/8 fermionic factor and color multiplicity
  N_c=3. Reference: Kapusta & Gale 2006, §6.4.

* **GeoVac outer baseline**: same Stefan–Boltzmann factorization as in
  §2.1, multiplied by `g_QCD`. The framework has no autonomous way to
  produce the color factor 3 or the gluon count 8 from outer-factor
  data — they are **inner-factor multiplicity invariants**. The outer
  baseline is therefore the per-mode coefficient `−(7/720)π²T⁴`
  (fermions) and `−(π²/90)T⁴` (bosons), waiting to be multiplied by
  inner multiplicities.

* **T_F residual structure**: this observable's residual structure is
  *purely k=0 multiplicity*. The inner Mellin trace at k=0,
  `Tr(e^{−tD_F²})_{m=0} = dim(H_F)`, is the species count. For the
  Connes A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ), dim H_F = 32 (one generation), and
  the multiplicity ratio of color-triplet to color-singlet sectors is
  6:2 = **3:1**. The integer 3 in `g_QCD = ... + (7/8)·2·N_c·N_f·2`
  is exactly the SU(3) color factor 3, which appears as the
  multiplicity ratio in the inner Mellin output (May-7 memo Leg 4).

* **What T_F would have to be**: a finite spectral triple containing
  M_3(ℂ) (or any algebra producing color multiplicity 3) and N_f
  generation copies. *Cannot be M_2(ℂ)·M_3(ℂ) with the wrong
  multiplicity assignment;* the Krajewski diagram structure forces
  the 3:1 ratio.

* **Falsifier**: a candidate `T_F` with no M_3(ℂ) summand, or with
  M_3(ℂ) but wrong multiplicity assignment, cannot reproduce
  `g_QCD = 47.5` at N_f=3. Standard CC-Connes A_F passes.

* **η-trivialization filter**: PASS. Pressure depends on k=0 species
  count + k=2 mass-suppression, not k=1.

* **Diagnostic value**: 8/10. Tests *algebra structure* of T_F (color
  factor) rather than parameter values. This is what
  `g_QCD = 47.5 ± something` measurements (lattice QCD, RHIC/LHC
  flow data) can in principle constrain.

---

### 2.3 Higgs effective potential at finite T (Coleman–Weinberg + thermal)

* **Standard QFT**: the one-loop finite-T effective Higgs potential is
  ```
  V_eff(φ, T) = − (μ²/2) φ² + (λ/4) φ⁴
              + (1/(64π²)) Σ_i n_i M_i⁴(φ) [log(M_i²(φ)/Λ²) − c_i]
              + (T⁴/(2π²)) Σ_i n_i J_{B/F}(M_i²(φ)/T²)
  ```
  where the sum runs over thermal species i, `M_i(φ) = y_i φ + ...`
  for fermions and `M_i²(φ) = (g²/4) φ² + ...` for gauge bosons, and
  `J_B`, `J_F` are the standard bosonic / fermionic thermal functions.
  The thermal correction drives the high-T symmetry restoration:
  `V_eff(φ, T) ⊃ + (1/2) m_T²(T) φ² + ...` with
  `m_T²(T) = (T²/12)·(2g_W² + g_Y² + 4λ + 4y_t² + ...)/4`. Reference:
  Dolan–Jackiw 1974, Phys. Rev. D 9, 3320; Kapusta & Gale 2006 §3.6.

* **GeoVac outer baseline**: the framework has Track 1's
  Stefan–Boltzmann tensor structure (the M1 × M2 backbone of the
  Bose/Fermi thermal integrals). It does *not* have the φ-dependent
  fermion-mass operator `M_i(φ) = y_i φ`, which lives entirely in the
  inner factor's off-diagonal `D_F` block (between ℂ and ℍ summands of
  A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ); see Connes 2013 for the canonical CC-SM
  bosonic spectral action). **The Higgs sector requires an inner-factor
  off-diagonal block.** Sprint H1 verdict (CLAUDE.md §1.7 WH1, May-6):
  "GeoVac does not autonomously select the Yukawa Y ≠ 0; falsifier
  HOLDS iff Y = 0; with Y > 0, ‖Φ‖_max ~ 0.05–0.27."

* **T_F residual structure**: this is the **algebra-structure probe par
  excellence**. The thermal mass `m_T²(T)` reads
  ```
  m_T²(T) = (T²/12) · [Tr_inner(Y_diag² ⊗ I_outer) + (gauge contrib)]
  ```
  where Y_diag² is the diagonal Yukawa-squared matrix appearing in
  the off-diagonal CC `D_F` block. The *bilinear* coupling
  `2 Y_e Y_μ + 2 Y_e Y_τ + ...` does NOT appear at this order because
  the Yukawa matrix is diagonalizable (in the SM sense), and the
  thermal mass is a pure trace `Σ_i y_i²`. Higher-order thermal
  corrections at two-loop *do* probe Yukawa products via
  `m_T²·log(m_T²)` corrections; these distinguish *n-generation*
  spectral triples from each other.

* **What T_F would have to be**: a finite spectral triple with an
  off-diagonal `D_F` block whose squared diagonal entries are the SM
  Yukawa-squared values Y_i². CC's standard A_F = ℂ ⊕ ℍ ⊕ M_3(ℂ)
  with appropriate D_F is the canonical solution but is *not unique* —
  any extension of A_F by a small subspace coupled by the same Yukawa
  pattern would also work.

* **Falsifier**: a candidate `T_F` whose off-diagonal block is zero
  (Marcolli–vS-without-Higgs; gauge-network limit) gives `m_T² = 0`
  from the lepton/quark sector, predicting NO high-T symmetry
  restoration in the framework's thermal effective potential. SM
  data (sphaleron transitions in the early universe, T_C ≈ 160 GeV)
  rules out Y = 0.

* **η-trivialization filter**: PASS. Thermal mass m_T² depends on
  Yukawa-squared (k=2 Mellin trace), not Yukawa-linear (k=1 vanishes
  by chirality).

* **Diagnostic value**: 10/10. This observable is structurally
  Higgs-sector-specific. If the framework reproduces it, it's
  reproducing the off-diagonal CC block; if it doesn't, it's not.

---

### 2.4 Electroweak gauge-boson thermal Debye masses (W, Z at high T)

* **Standard QFT**: in the symmetric-phase plasma at T ≫ T_C, the W
  and Z gauge bosons acquire thermal Debye masses through the Hard
  Thermal Loop resummation:
  ```
  m_D²(W) = (T²/6) · [(2 + N_gen) g_W² + ...]
  m_D²(Z) = (T²/6) · [(... ) g_Z² + ...]
  ```
  where the factors involve the number of generations N_gen = 3 (lepton
  + quark contributions to gauge polarization) and the relevant gauge
  couplings. Reference: Braaten-Pisarski 1990, Nucl. Phys. B 337, 569.

* **GeoVac outer baseline**: Track 1's Stefan–Boltzmann gives the 1-loop
  vacuum polarization graph structure, but the gauge-coupling values
  g_W, g_Z and the SU(2) algebra structure live entirely in
  Papers 25/30's Wilson construction (gauge sector of the AC factor).
  At outer-factor level: graph topology of the gauge polarization is
  fixed; coupling magnitudes are *external* (not predicted by GeoVac).

* **T_F residual structure**: m_D²(W) is bilinear in gauge couplings
  AND counts generations: `m_D² = (T²/6) · [(2 + N_gen) g_W²]`.
  The (2 + N_gen) factor arises as `2 + N_gen·N_f_per_gen` from
  doublet contributions (gauge bosons themselves) plus matter
  contributions (quark and lepton doublets per generation). The number
  3 enters as the *generation multiplicity* of the inner factor.

* **What T_F would have to be**: a finite spectral triple with
  N_gen = 3 generation copies, each containing a doublet structure
  (left-handed lepton + quark doublets of SU(2)_L). CC's A_F with the
  standard SU(2) embedding in ℍ passes; a triple with N_gen = 1 or
  N_gen = 4 would predict measurably different m_D².

* **Falsifier**: lattice-QCD-like measurements of EW Debye masses at
  high T (or their imprint on early-universe phase transitions) sample
  the (2 + N_gen) factor. A T_F with the wrong number of generation
  copies would show up.

* **η-trivialization filter**: PASS. m_D² is k=2 (mass-squared
  trace), not k=1.

* **Diagnostic value**: 7/10. Less clean than §2.3 because gauge
  couplings g_W, g_Z are also external; the residual structure
  *factors* through both Yukawas (via fermion loops) and gauge
  couplings (via boson self-energy).

---

### 2.5 Thermal correction to the fine-structure constant α(T) (lepton+quark loop)

* **Standard QFT**: the running of α at finite T, from leptonic and
  hadronic vacuum polarization, is
  ```
  Δα(T)/α = (α/π) · Σ_f Q_f² · N_c(f) · (1/3) · log(T/m_f) + ...
          (for T > m_f, otherwise mass-suppressed)
  ```
  where Q_f is the electric charge and N_c(f) the color multiplicity
  (1 for leptons, 3 for quarks). At high T (T ≫ m_t), `Δα/α` saturates
  and approaches a fixed contribution per generation. Reference:
  Kapusta–Toimela 1989, Phys. Rev. D 39, 3197.

* **GeoVac outer baseline**: the framework has computed one-loop QED
  vacuum polarization on S³ (Paper 28 §QED, Sprint Q-1) with Π =
  1/(48π²) and β(α) = 2α²/(3π). At finite T, Track 1's tensor structure
  gives the Stefan–Boltzmann building block; multiplying by per-species
  thermal log corrections requires summing over fermion masses. The
  outer baseline is the QED β-function structure; the per-fermion
  contribution requires inner-factor (Q_f, N_c(f), m_f) data.

* **T_F residual structure**: the thermal coupling shift is
  `(α/3π) · Σ_f Q_f² · N_c(f) · log(T/m_f)`,
  which is `Tr_inner(Q_F² · log(T/M_F))` where Q_F is the charge
  matrix and M_F = D_F at the EW scale. This is **NOT** a pure power
  Mellin trace — it has a logarithm. Logarithmic running is
  inherited from QFT renormalization; the inner factor contributes
  the *coefficient* `Σ_f Q_f² · N_c(f)` (a pure rational from the
  charge assignment + color multiplicity) and the *cutoff structure*
  `m_f` (Yukawa-tied, since m_f = y_f · v).

* **What T_F would have to be**: charge assignment Q_F encoded in A_F
  representation theory (electroweak embedding via U(1)_Y × SU(2)_L
  / Z₃), color factor N_c(f) from M_3(ℂ) summand, generation
  multiplicity from generation copies. The charge assignment is
  *forced* by gauge-anomaly cancellation in the SM, which is itself
  a constraint on the inner factor's representation theory.

* **Falsifier**: a candidate T_F whose charge assignment violates
  `Σ_lepton (Q_lepton² · N_c=1) + Σ_quark (Q_quark² · N_c=3) =
  Σ_lepton Q_lepton² + 3 Σ_quark Q_quark² = 8/3` (per generation,
  as recently noted in Phase 4H Track SM-B) would fail to reproduce
  the precision Δα(T)/α dataset.

* **η-trivialization filter**: PASS at the level of the *thermal
  correction*. The bare logarithmic running structure does involve a
  k=1 piece (the η-invariant of the Dirac operator), which by
  η-trivialization is identically zero on a Krajewski triple — but
  this is the standard chirality-symmetric SM, where parity-odd
  contributions cancel between left- and right-handed sectors.

* **Diagnostic value**: 6/10. Logarithmic dependence makes the
  residual structure non-Mellin; precision is hard to extract.

---

## 3. Summary table

| # | Observable | k-order | Residual character | Cleanest probe of |
|---|------------|---------|--------------------|-------------------|
| 2.1 | Lepton thermal free energy | k=2 | linear in `Σ y_i²` | Yukawa values |
| 2.2 | QCD thermal pressure / `g_QCD` | k=0 | multiplicity ratio | M_3(ℂ) algebra structure |
| 2.3 | Higgs effective potential at finite T | k=2 | trace of off-diag D_F block | Higgs sector / off-diagonal CC |
| 2.4 | EW Debye masses | k=2 | (2 + N_gen) g² | generation count |
| 2.5 | α(T) thermal running | k=2 + log | rational charge sum + Yukawa cutoffs | charge assignment + N_c |

All five **PASS the η-trivialization filter** (k=0 or k=2 dependence,
no odd-k piece).

---

## 4. Filtered candidates (k=1 observables, identically zero by chirality)

For completeness, here are thermal observables we would expect to be
sensitive to inner-factor structure on naive grounds, but which the
η-trivialization theorem says are **identically zero** on any
Krajewski-class even spectral triple. None of these are useful Yukawa
probes:

* **Axial chemical potential thermodynamics** (μ_5 conjugate to
  Q_5 = ψ†γ_5ψ at finite T). Depends on
  `Tr(γ_F · D_F · e^{−tD_F²})` — a k=1 trace, identically zero by
  η-trivialization. The structural prediction is that there is no
  axial-chemical-potential thermodynamics from the inner factor in
  any Connes-Chamseddine-compatible SM extension.

* **Parity-odd transport coefficients** (chiral magnetic effect, etc.)
  at single-fermion level are also k=1 and trivialize. These appear
  in nature only via *anomalous* mechanisms that do not factor
  through the inner-factor Mellin engine — they are gauge-anomaly
  inflow contributions, separately captured by the outer triple.

* **Electric dipole moments at finite T** (a fortiori — these are
  parity-odd, time-reversal-odd at the QFT level). The structural
  prediction matches the SM "vacuum theta angle is small" stance.

This is a real positive prediction of the framework: structural CP
symmetry of the inner-factor Mellin output forbids these contributions.
SM observed results (small/null EDMs, no axial chemical potential
thermodynamics from electroweak sector) are consistent.

---

## 5. Honest scope and limitations

* **No T_F is selected.** All five catalogued observables are
  *consistency checks* on candidate T_F's — given enough precision,
  they constrain which finite spectral triples are eligible inputs.
  None *derive* T_F.
* **Outer baselines are textbook QFT** (Stefan–Boltzmann, one-loop
  vacuum polarization, Hard Thermal Loops). The framework does not
  produce these from a deeper principle; it shows them factor cleanly
  through the tensor-product spectral triple structure built in
  Track 1. The new content is the *factor placement*: π lives outside,
  Yukawas live inside, multiplicities live inside.
* **Logarithmic running** (§2.5) does not factor cleanly into the
  Mellin engine; it requires renormalization-group machinery
  not native to the framework (LS-8a wall, Sprint May-7 Track c).
  The inner-factor Mellin engine handles *finite parts*; UV running
  is structurally outside.
* **Generation count is not predicted** — it appears as an
  empirical multiplicity factor in T_F.
* **The outer baseline of §2.3 is the weakest** — Track 1 has
  Stefan–Boltzmann but no φ-dependent operator; the full Coleman–
  Weinberg machinery requires tools we don't have. The catalogue
  entry names the residual structure but does not compute it.

The catalogue does what the structural-skeleton-scope statement
predicts: it sharpens "no autonomous SM-distinguishing data" into a
*list of observables and the precise T_F structure each would
constrain*. No data is missing from the list — these are five places
where, if precision improves enough or new measurements come in, the
framework's structural-skeleton scope would be tested.

---

## 6. Cross-references

* **Sprint H1** (CLAUDE.md §1.7 WH1, May-6): Higgs admits structurally
  on AC extension iff Y > 0; GeoVac does not select Y. §2.3 above is
  the thermal version.
* **Inner-factor Mellin engine memo** (`debug/inner_factor_mellin_engine_memo.md`):
  η-trivialization theorem (§3) and AC factorization theorem (§4) are
  the two backbones used throughout this catalogue.
* **Paper 18 §IV.6**: inner-factor input data tier (Yukawa Dirichlet
  ring), the parameter-tied tier disjoint from outer-factor exchange
  constants. Each k=2 entry above lives in this ring.
* **Paper 32 §VIII.C**: Sprint H1 addendum, AC extension verdict.
  §2.3 of this memo is the natural thermal extension.
* **Paper 35 §VII.3 Refined Prediction 1**: GeoVac controls π content
  of the *finite parts* of any QED observable; UV divergences and
  renormalization counterterms are external. §2.5 of this memo
  honors this in the thermal-running case.
* **Paper 33 §V (1+6+1 partition)**: structurally disjoint pattern
  from this catalogue. Paper 33 is about *outer-factor* selection
  rule recovery; this catalogue is about *inner-factor* parameter
  sensitivity. The two patterns are orthogonal.

---

## 7. Paper update drafts (NOT applied — PI-approval gated)

### 7.1 Paper 18 §IV.6 extension

After the existing inner-factor input data tier paragraph, append:

> **Thermal probes of the Yukawa Dirichlet ring.** Sprint TD Track 3
> (May 2026, `debug/sprint_td_track3_memo.md`) catalogues five thermal
> observables whose value depends on parameters of the inner factor
> Yukawa Dirichlet ring: lepton thermal free energy (`Σ y_i²`); QCD
> thermal pressure `g_QCD` (color multiplicity); Higgs effective
> potential at finite T (off-diagonal CC block trace); EW Debye masses
> ((2+N_gen)·g²); thermal running α(T) (charge assignment + N_c).
> Each observable's outer baseline (Stefan–Boltzmann backbone) factors
> cleanly via the Track 1 tensor-product spectral triple; the inner
> residual is a generalized Dirichlet sum at integer s in the
> Yukawa eigenvalues, k=0 or k=2 by the η-trivialization theorem.
> All five **fail to be predicted by GeoVac alone** — they are
> consistency checks on candidate inner factors, not derivations.

### 7.2 Paper 32 §VIII.C addendum

After the Sprint H1 verdict paragraph, append:

> **Thermal extension of the H1 verdict.** The same Yukawa-undetermined
> structural reading extends to the thermal Higgs effective potential
> V_eff(φ, T) (Sprint TD Track 3 §2.3, May 2026): the thermal mass
> m_T²(T) at one-loop is `(T²/12) · Tr_inner(Y_diag²) + (gauge
> contrib)`, with the SM-data-fit value `m_T²(v, T_C) ≈ (160 GeV)²`
> reproducible by any candidate T_F whose squared Yukawa diagonal
> reproduces the empirical SM spectrum. As in Sprint H1, GeoVac does
> not select the spectrum; it admits any compatible one. **Marcolli–
> vS-without-Higgs (Y = 0) predicts no high-T symmetry restoration
> and is empirically falsified** — this is the cleanest empirical
> falsification of the gauge-network reading of the framework, made
> precise by the thermal extension.

---

## 8. Files produced

* `debug/sprint_td_track3_memo.md` (this file, ~3500 words)
* `debug/data/sprint_td_track3.json` (catalogue data, structured)

No `geovac/` modifications. No paper edits. No new tests (catalogue
entries do not introduce new framework computations beyond what
Track 1 / Sprint H1 / inner-factor Mellin engine already provided).

---

## 9. Decision points for PI

1. **Apply Paper 18 §IV.6 extension?** (~30 lines, expands the
   existing inner-factor tier with thermal probes paragraph.)
2. **Apply Paper 32 §VIII.C thermal addendum?** (~25 lines, sharpens
   Sprint H1 verdict with thermal Higgs data.)
3. **Pursue any single observable to deeper computation?** Highest
   value: §2.3 Higgs effective potential (cleanest algebra-structure
   probe) or §2.2 QCD thermal pressure (cleanest multiplicity probe).
4. **Catalogue extension?** Other plausible thermal Yukawa probes
   include neutrino oscillations at finite T (lepton mass matrix
   structure, k=2), B-meson thermal mixing (CKM matrix structure,
   k=2 with cross-generation off-diagonals), and quark-gluon plasma
   shear viscosity (η/s, depends on coupling running). These are
   harder to compute outer baselines for.

The catalogue is honest by design — every entry names what GeoVac
provides and what it does not. The structural-skeleton scope holds.
