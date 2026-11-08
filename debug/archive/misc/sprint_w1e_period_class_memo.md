# Sprint W1e period-class diagnostic — canonical memo

**Date:** 2026-06-04 (same-day execution).
**Sprint position:** Diagnostic test of the hypothesis that the W1c→W1e chemistry-correlation wall (NaH inner-region overattraction) lives in a higher-weight cyclotomic mixed-Tate period class (M3) while the diagonal sits in pure-Tate (M1), under the master Mellin engine partition of Paper 18 §III.7 / Paper 32 §VIII case-exhaustion theorem.
**Verdict line: STOP because 0/11 W1e correction terms identify with low-coefficient outer-factor periods (M1, M2, or M3) after audit filtering; period-class framing was the wrong axis — W1e corrections are FCI eigenvalues and Hartree integrals built from external chemistry-side calibration data (Clementi-Raimondi exponents), structurally in the Class 1 / inner-factor input-data tier (Paper 18 §IV.6 chemistry-side analog), categorically disjoint from the outer-factor M1/M2/M3 engine.**
**Cross-references:** `debug/sprint_f2_cross_vne_kernel_memo.md` (W1d closure, W1e named), `debug/sprint_f4_bonding_pk_memo.md` (W1e PK insufficient), `debug/sprint_f5_explicit_core_memo.md` (W1e Hartree partial), `debug/sprint_f6_maxn4_nah_memo.md` (W1e basis enlargement partial), `papers/group3_foundations/paper_18_exchange_constants.tex` (§III.7 master Mellin engine, §IV.6 inner-factor input data), `papers/group2_quantum_chemistry/paper_19_coupled_composition.tex` §sec:w1c_residual, CLAUDE.md §1.7 multi-focal-composition wall pattern, MEMORY.md `external_input_three_class_partition`.

---

## §0. Executive summary

The hypothesis under test: the W1c→W1e chemistry wall lives in a structurally different period class than the diagonal of the composed Hamiltonian. If true, this would *structurally explain* why six modification attempts (CLAUDE.md §3 entries: PK cross-center, screened-Schrödinger valence, multi-zeta substitution, three-bucket M-Z partition, kernel-shape substitution, bonding-PK rank-1) all failed at sprint scale: they would be trying to fit a higher-weight cyclotomic-mixed-Tate-level-4 transcendental with pure-Tate machinery.

The diagnostic gathered all eleven W1e correction terms from the F4/F5/F6 sprint outputs at NaH ($R = R_{\rm eq} = 3.566$ bohr) and ran PSLQ at 100 dps against four bases:
- **M1** (Hopf-base measure, $\mathbb{Q}[\pi, 1/\pi]$, 4 transcendental generators)
- **M2** (Seeley-DeWitt, $\bigoplus_k \pi^{2k} \cdot \mathbb{Q}$ on unit $S^3$, 3 generators)
- **M3** (vertex-parity Hurwitz, $\mathcal{MT}(\mathbb{Z}[i, 1/2])$ at level $\le 4$, 5 generators: $\log 2$, Catalan $G$, $\beta(4)$, $\zeta(3)$, $\zeta(5)$)
- **INNER** (chemistry-side analog of Paper 18 §IV.6 inner-factor input data: Clementi-Raimondi $Z_{\rm eff}$ exponents for the Na [Ne] core + hydrogenic energies built from them, 6 generators)

Audit filtering rejects three classes of spurious "hit":
1. **Rational-only fits** (transcendental coefficient sum = 0): float64 input precision admits low-denominator rational identifications that are not period assignments (e.g.\ 0.194 = 97/500, target=$-7$ constant=$-16$ ...).
2. **Basis-internal fits** (target coefficient = 0): PSLQ finding a Q-linear dependency within the basis (the basis itself is over-complete) rather than a target identification.
3. **Curve-fit-audit failure mode** (basis size × input precision): for an $n$-element basis at $p$-digit input precision, PSLQ at ceiling $C$ will find spurious relations when $p/(n-1) < \log_{10} C$. For 6-element INNER basis + 3-4 digit float inputs (the F4 PK barrier "0.194", the F6 baseline "4.374", the experimental NaH $D_e$ "0.0746") at audit ceiling 100, this gate is violated.

Results at audit ceiling $C = 100$ (precision-safe):

| Basis | hits (raw) | spurious_rational | basis_internal | genuine target-coefficient hits |
|:------|---:|---:|---:|---:|
| **M1** | 0 | 0 | 0 | **0** |
| **M2** | 3 | 3 | 0 | **0** |
| **M3** | 0 | 0 | 0 | **0** |
| **INNER** | 8 | 0 | 5 | 3 (all at curve-fit-audit failure regime) |

**Zero genuine outer-factor (M1, M2, M3) identifications across 11 correction terms.** The three INNER "real hits" are at audit ceiling 100 with 6-7 basis elements vs 3-4 digit input precision — the curve-fit-audit failure mode (any value fits at this basis-to-precision ratio).

**Net structural finding:** W1e correction terms are **not periods of the outer-factor Mellin engine**. They are FCI eigenvalues and Hartree integrals computed from external chemistry-side calibration data (Clementi-Raimondi 1963 atomic Hartree-Fock fits). They live in the **Class 1 / inner-factor input-data** tier (Paper 18 §IV.6 chemistry-side analog), categorically disjoint from M1/M2/M3.

The original hypothesis — "W1e lives in M3, diagonal in M1" — is structurally rejected. The reason no rational/algebraic modification of the framework has closed W1e is NOT that the framework is trying to fit a higher-weight period; it is that the W1e content lies entirely outside the outer-factor period ring. This is the precise chemistry-side analog of the H1 Yukawa non-selection theorem (Sprint H1, Paper 32 §VIII.C): the framework's outer-factor machinery cannot autonomously generate inner-factor data.

This is a sharpening, not a deflation, of the structural-skeleton-scope finding (CLAUDE.md §1.7, MEMORY.md `external_input_three_class_partition`). The wall's structural location is now explicit: it is in the **inner-factor input-data tier**, and W1e's six failed closure attempts are accounted for as attempts to derive Class-1 calibration content from Class-2 outer-factor machinery.

---

## §1. Method

### §1.1 Correction terms extracted

Eleven NaH correction terms loaded from prior-sprint JSONs (precision ~15 digits = float64 from JSON round-trip):

| Label | Source sprint | Value (Ha) | Mechanism |
|:------|:-----|---:|:------|
| `F4_baseline_De_Ha` | F4 §1.4 | 4.37346 | F3 architecture FCI well depth (no closure) |
| `F4_PK_barrier_Ha` | F4 §1.2 | 0.194 | Predicted Phillips-Kleinman rank-1 barrier |
| `F4_De_at_predicted_PK_Ha` | F4 §1.4 | 4.37117 | FCI well depth at predicted PK shift |
| `F4_correction_at_PK_Ha` | derived | 0.00229 | FCI sensitivity to PK rank-1 (0.05% wall closure) |
| `F5_J_total_Ha` | F5 §1.3 | 1.13139 | Hartree mean-field of [Ne] core on bonding orbital |
| `F5_K_total_Ha` | F5 §1.3 | 0.00815 | Slater-Condon exchange estimate (K/J = 0.72%) |
| `F5_JmK_correction_Ha` | F5 §1.3 | 1.12324 | Predicted explicit-core J − K correction (25.7% wall closure) |
| `F6_well_depth_at_max_n_4_Ha` | F6 §3 | 3.92883 | NaH well depth at max_n=4 enlarged basis |
| `F6_F3_baseline_Ha` | F6 | 4.374 | F3 baseline (round number from sprint metadata) |
| `F6_closure_Ha` | derived | 0.44517 | Wall closure from basis enlargement (10.2%) |
| `NaH_De_exp_Ha` | external | 0.0746 | Huber-Herzberg experimental binding energy |

### §1.2 Period bases (independent generators only)

Constants `1` are excluded from M1/M2/M3 bases because basis-with-constant + float64 input ALWAYS finds spurious rational fits at low coefficient ceiling. INNER basis has 6 generators chosen as the three CR exponents plus three hydrogenic energies built from them.

**M1** (Hopf-base measure, $k=0$ Mellin):
$\{\pi, 1/\pi, \pi^2, 1/\pi^2\}$ — 4 generators, pure-Tate $\mathbb{Q}[\pi, 1/\pi]$.

**M2** (Seeley-DeWitt, $k=2$ Mellin):
$\{\pi^2, \pi^4, \pi^6\}$ — 3 generators. (Raw heat-trace has $\sqrt\pi$ which cancels against $(4\pi)^{3/2}$ on unit $S^3$ to give pure-Tate output; Sprint Mixed-Tate Test 2026-06-03 verified GeoVac discrete S³ M2 sits in this sub-ring with no $\zeta(3)$.)

**M3** (vertex-parity Hurwitz, $k=1$ Mellin):
$\{\log 2, G, \beta(4), \zeta(3), \zeta(5)\}$ — 5 generators in $\mathcal{MT}(\mathbb{Z}[i, 1/2])$ at level $\le 4$ (Sprint M3 Cyclotomic Mixed-Tate 2026-06-03; Deligne 2010, Glanois 2015 Cor.~1.1–1.2).

**INNER** (chemistry-side inner-factor analog, Paper 18 §IV.6 reading):
$\{\zeta_{1s}^{CR}, \zeta_{2s}^{CR}, \zeta_{2p}^{CR}, E_{1s}^{CR}, E_{2s}^{CR}, E_{2p}^{CR}\}$ — 6 generators built from Clementi-Raimondi 1963 Table II exponents for Na [Ne] core. The hydrogenic energies are $E_{1s} = \zeta_{1s}^2/2$, $E_{2s} = \zeta_{2s}^2/8$, $E_{2p} = \zeta_{2p}^2/8$, which are the closed-form rationals-in-$Z_{\rm eff}$ that the Hartree J integral produces algebraically.

### §1.3 Audit filtering

Three filter classes (added incrementally during the sprint as initial runs returned obvious artifacts):

1. **Rational-only fit filter**: a PSLQ hit with $\sum_i |c_i| = 0$ across all NON-constant basis elements means PSLQ matched the float64 target value to a low-denominator rational. This is an input-precision artifact (e.g.\ 0.194 = 97/500), not a period identification.

2. **Basis-internal-only filter**: a PSLQ hit with target coefficient = 0 means PSLQ found a Q-linear dependency WITHIN the basis (the basis is over-complete) rather than identifying the target. The INNER basis has six generators that are not Q-linearly independent at audit precision — only three are truly independent (the three CR exponents); the three energies are $\zeta^2 / (2n^2)$ which is quadratic and admits low-coefficient relations in the basis at audit ceiling.

3. **Curve-fit-audit failure mode** (advisory, not auto-filtered): for $n$-element basis at $p$-digit precision and ceiling $C$, spurious fits dominate when $p/(n-1) < \log_{10} C$. For 6-element INNER basis at 3-4 digit input precision (the F4 PK = "0.194", the F6 baseline = "4.374", the experimental $D_e$ = "0.0746") and audit ceiling 100, this gate is violated. These 3 INNER hits are advisory-flagged in the JSON output (`audit.curve_fit_audit_failure_mode_flag` named follow-on).

---

## §2. Results

### §2.1 PSLQ run at audit ceiling 100

| Term | M1 hit | M2 hit | M3 hit | INNER hit |
|:-----|:------:|:------:|:------:|:---------:|
| F4_baseline_De_Ha | — | — | — | **basis-internal** (target=0) |
| F4_PK_barrier_Ha | — | rational-only | — | **CFA-failure** (target=11, max\|c\|=41) |
| F4_De_at_predicted_PK_Ha | — | — | — | — |
| F4_correction_at_PK_Ha | — | — | — | basis-internal |
| F5_J_total_Ha | — | — | — | — |
| F5_K_total_Ha | — | — | — | basis-internal |
| F5_JmK_correction_Ha | — | — | — | — |
| F6_well_depth_at_max_n_4_Ha | — | — | — | basis-internal |
| F6_F3_baseline_Ha | — | rational-only | — | **CFA-failure** (target=10, max\|c\|=26) |
| F6_closure_Ha | — | — | — | basis-internal |
| NaH_De_exp_Ha | — | rational-only | — | **CFA-failure** (target=39, max\|c\|=70) |

Summary: **0 genuine outer-factor hits**, 3 rational-only M2 artifacts (filtered), 5 basis-internal INNER artifacts (filtered), 3 INNER hits at curve-fit-audit failure regime (advisory).

The three rational-only M2 "hits" are diagnostic: 0.194 = 97/500, 4.374 = 2187/500, 0.0746 = 373/5000. These are round-number float artifacts (the JSONs stored short decimal strings), not framework-derived periods.

### §2.2 PSLQ at permissive ceiling $10^6$

Same pattern: zero genuine M1/M2/M3 hits, INNER hits all in CFA failure regime. Permissive ceiling adds more spurious INNER fits but no new M1/M2/M3 content.

### §2.3 Cross-check against M3 specifically

The original hypothesis was that W1e sits in M3 (cyclotomic mixed-Tate at level 4) and the diagonal in M1. The M3 column shows **zero hits at any ceiling, with or without filtering**. Catalan $G$, $\beta(4)$, $\zeta(3)$, $\zeta(5)$ have nothing to do with the W1e correction terms.

This is consistent with M3's structural reading (Paper 18 §III.7, §IV.6): M3 lives at $k=1$ in the master Mellin engine, ties to **vertex-parity** spectral sums on the GeoVac graph, and on Krajewski-class inner factors $\eta$-trivializes (Theorem `thm:eta_trivialization`, Paper 18). The W1e correction terms have NO vertex-parity content (they are FCI eigenvalues and Hartree integrals over hydrogenic + multi-zeta orbitals); they cannot pick up M3 content even structurally.

---

## §3. Structural classification — where does W1e actually live?

### §3.1 The terms-by-type table

| Term | Underlying object | Free parameters | Closed-form status |
|:---|:---|:---:|:---|
| F4 PK barrier | $\sum_c (E_v - E_c) \|S_{bc}\|^2$ rank-1 projector eigenshift | 5 CR exponents + 4 bonding-orbital coefficients | Empirical float (input data not derived) |
| F5 Hartree J | $\int \|\phi_b\|^2 V_{\rm core}$ Coulomb integral against [Ne] mean-field | 5 CR exponents + 4 bonding coefficients | Closed-form rational-in-$Z_{\rm eff}$, but $Z_{\rm eff}$ are external fits |
| F5 Exchange K | $\sum_c \|S_{bc}\|^2 F^0(c,c)$ Slater-Condon | 5 CR exponents + 4 bonding coefficients | Same — closed-form-in-fits |
| FCI eigenvalues (F4/F6 baselines, F6 max_n=4) | Ground state of 100-3600 determinant FCI matrix | 5-20 multi-zeta exponents + framework Z, R inputs | **Algebraic-implicit** (Paper 18 algebraic-implicit tier): root of a characteristic polynomial with transcendental-input matrix elements, no closed form even with rational inputs |
| Experimental $D_e$ | Measured binding energy | — | External observation (Class 1) |

### §3.2 The structural classification

The W1e correction terms fall into three structurally distinct sub-categories, none of which are outer-factor M1/M2/M3 periods:

1. **External input data (Class 1 / inner-factor input-data tier, Paper 18 §IV.6 chemistry-side analog)**: the Clementi-Raimondi 1963 atomic Hartree-Fock fits for the [Ne] core exponents. The Hartree J integral is rational in these exponents, but the exponent values themselves are external calibration data, structurally categorical to GeoVac's outer-factor mechanism. This is the chemistry-side analog of Sprint H1's SM Yukawa values: rational structure in external generators, but the generators are not framework-derived.

2. **Algebraic-implicit content (Paper 18 algebraic-implicit tier)**: the FCI eigenvalues themselves. These are roots of a characteristic polynomial whose entries are rational + multi-zeta-fitted. The eigenvalues do not have closed form in general (Paper 18 §III: algebraic-implicit means defined by polynomial $P = 0$ with known coefficient ring, but pointwise diagonalization is computational convenience, not mathematical necessity for closed form to exist; here NO closed form exists even in principle for the 100-dim FCI characteristic polynomial roots).

3. **Class-1 external observations**: the experimental $D_e$.

None of (1), (2), (3) are outer-factor periods in the master Mellin engine sense. The original hypothesis confused "deep transcendental wall" with "wrong period class". The actual structural reason W1e cannot be closed by rational/algebraic modification is **not** that the framework needs a higher-weight cyclotomic-mixed-Tate level-4 period; it is that the wall content **is not a period at all** — it is computer-extracted FCI output built from externally-fitted atomic-physics data.

### §3.3 Consistency with prior structural findings

The diagnostic confirms three prior framework-internal findings:

1. **Multi-focal-composition wall pattern** (MEMORY.md, CLAUDE.md §1.7, 6 instances including H1 Yukawa, LS-8a, HF-3/4/5, W1e chemistry): W1e is the chemistry-side analog of the LS-8a renormalization-counterterm wall (LS-8a gap in QED two-loop self-energy needs UV-completion counterterms; W1e needs calibration-data atomic-physics input). Both walls share the structural reason: the framework can compose pre-computed external input cleanly via its bimodule machinery, but cannot autonomously generate the external input.

2. **External-input three-class partition** (MEMORY.md `external_input_three_class_partition`): W1e cleanly classifies as **Class 1 (calibration data)** + algebraic-implicit FCI eigenvalue tier. The Class 2 (multi-focal composition) and Class 3 (multi-determinant correlation) labels apply to the *mechanism* of how W1e composes; the *data* tier is Class 1.

3. **H1 Yukawa non-selection theorem** (Sprint H1, Paper 32 §VIII.C): the SM Yukawa values are gauge-forced but value-free; the framework's $A_F = \mathbb{C} \oplus \mathbb{H} \oplus M_3(\mathbb{C})$ admits any Yukawa $Y$ for non-zero Higgs but cannot autonomously select $Y$. The chemistry-side analog: the framework's composed_qubit + balanced_coupled architecture admits any Clementi-Raimondi exponent set for the [Ne] core but cannot autonomously select the exponents. The "non-selection" extends from inner-factor SM data to chemistry-side atomic-physics data, with the same structural cause (outer-factor machinery is silent on inner-factor data).

### §3.4 Why the M3 hypothesis was tempting (and wrong)

The "W1e sits in M3" reading was tempting because:
- M3 is the third pillar of the master Mellin engine (the "harder" one).
- F5 found that the Hartree J correction has the right sign but wrong magnitude (25.7% closure); F4/F6 saturated at 43% / 10.2%.
- A natural reading: the framework supplies pure-Tate (M1) machinery; the wall needs M3 content; the residual is the M1↔M3 mismatch.

The reading was wrong because M3 lives at $k=1$ in $\mathcal{M}[{\rm Tr}(D^k e^{-tD^2})]$, which is the **vertex-parity / parity-character spectral sum** sector. The W1e mechanism has zero vertex-parity content — it is a Hartree integral over hydrogenic orbital products + an FCI eigenvalue of a Coulomb-kernel matrix. There is no spectral parity sum anywhere in the W1e computation pipeline. M3 was the wrong target by structure, not by precision.

The correct structural location for W1e is the **inner-factor input-data tier** added to Paper 18 in v3.X (§IV.6) precisely to name where framework-external calibration data lives. The chemistry-side analog of the SM inner-factor data is the atomic-physics calibration set (Clementi-Raimondi exponents, multi-zeta coefficients, FrozenCore Z_{\rm eff}(r) profiles).

---

## §4. PROPOSED Paper 18 §IV.6 entry (do not apply to .tex in this sprint)

Suggested text block for Paper 18 §IV.6 ("Inner-factor input data" section, after the SM-side discussion):

> **Chemistry-side analog: external atomic-physics calibration data.** The same structural classification applies to GeoVac's composed-natural-geometry chemistry side. The composed_qubit architecture builds heavy-atom molecules (Z ≥ 11) using a frozen-core treatment: the [Ne], [Ar], [Kr], [Xe] cores are reduced to a screened $Z_{\rm eff}(r)$ profile via the FrozenCore solver, with hydrogenic + multi-zeta valence orbitals built on top. The five Clementi-Raimondi 1963 exponents per row, the multi-zeta coefficients, and the screened-core orbital shapes are external atomic-physics inputs to the framework — they are not derived from any outer-factor Mellin mechanism (M1/M2/M3).
>
> Sprint W1e period-class diagnostic (2026-06-04, `debug/sprint_w1e_period_class_memo.md`) tested whether the W1c→W1e chemistry-correlation wall sits in a higher-weight cyclotomic-mixed-Tate (M3) period class than the diagonal of the composed Hamiltonian. Eleven NaH correction terms extracted from sprints F4/F5/F6 were PSLQ-classified at 100 dps against M1, M2, M3, and an INNER analog basis built from the Clementi-Raimondi exponents. Result: zero genuine outer-factor (M1/M2/M3) identifications across all 11 terms after audit filtering for rational-only fits, basis-internal Q-linear dependencies, and the curve-fit-audit failure mode (basis size × input precision). The corrections classify cleanly into the chemistry-side analog of the inner-factor input-data tier: FCI eigenvalues are algebraic-implicit roots of polynomials with externally-fitted matrix elements; Hartree integrals are rational in external $Z_{\rm eff}$ values; PK barriers depend on the same external orbital fits.
>
> The chemistry-side data tier shares two structural features with the SM-side inner-factor data: (i) it factorizes cleanly through the framework's outer-factor architecture (composed_qubit + balanced_coupled compose external input via well-defined matrix-element machinery), and (ii) the framework provides no autonomous derivation principle for the external values. The chemistry analog of Sprint H1's Yukawa non-selection theorem is the W1e wall: the framework can use any Clementi-Raimondi exponent set but cannot autonomously select them. Six W1c/W1e closure attempts (CLAUDE.md §3 entries for PK cross-center, screened-Schrödinger valence, multi-zeta substitution, kernel-shape substitution, bonding-PK, basis enlargement, explicit-core Hartree) are explained at the structural level by this non-selection: they attempt outer-factor mechanisms (rational/algebraic modifications) on a wall whose content lies entirely in the external input-data tier.
>
> The five chemistry-side observables in Sprint TD Track 3's analog role (alkali-hydride binding energies, ionization energies, fine-structure splittings, hyperfine splittings, polarizabilities) are not predicted by GeoVac alone; they are consistency checks that constrain candidate atomic-physics calibration sets, not framework-derived predictions. This is the chemistry-side analog of the SM-side empirical finding that the Marcolli–van Suijlekom gauge-network reading without Higgs is falsified by sphaleron physics: the framework's outer-factor sector is consistent with — but does not select among — multiple inner-factor / calibration-data realizations.

---

## §5. Honest scope

What this sprint demonstrated:
- Eleven W1e correction terms at NaH ($R = R_{\rm eq}$) yield zero genuine outer-factor period identifications at audit ceiling 100 with 3-5 transcendental-content basis elements.
- INNER-basis hits are all in the curve-fit-audit failure regime (6-element basis vs 3-4 digit input precision); they support the structural reading but do not constitute identifications.
- M3 specifically (the original hypothesis) returns zero hits at any ceiling, with or without audit filtering — consistent with the structural reason that M3 lives in the vertex-parity sector that W1e cannot touch.

What this sprint did NOT demonstrate:
- High-precision (100 dps) recomputation of the FCI eigenvalues: the inputs were loaded from sprint JSONs at float64 precision. Recomputing at 100 dps would not change the structural reading (FCI eigenvalues are algebraic-implicit, not periods, regardless of computational precision), but would be needed for a final-cut PSLQ verdict at 100-dps precision.
- Cross-system robustness: the diagnostic ran only on NaH. LiH (W1c small), MgH₂ (W1e visible), H₂O (W1c+W1e mixed) would each need a separate F4/F5/F6 sprint cycle to extract analogous correction terms.
- Construction of the chemistry-side analog of the SM η-trivialization theorem (Paper 18 Thm `thm:eta_trivialization`): the analog claim would be that the chemistry inner-factor data ring has structural properties (e.g., closed under certain operations on the atomic-physics side) that mirror the chirality-grading axiom on the SM side. This is a multi-month NCG-research target, not a sprint-scale move.
- Falsification of the alternative reading: "the correction terms ARE outer-factor periods but at higher coefficient ceiling / higher precision than sprint-scale can probe." This is a non-falsifiable null hypothesis at sprint scale; the structural argument (M3 lives at $k=1$ vertex-parity, W1e has no vertex-parity content) is the actual close.

Decision-gate outcome: **STOP** (clean negative on the original M3 hypothesis; positive identification of W1e in the inner-factor input-data tier).

---

## §6. Follow-on register

| Priority | Item | Cost | Status |
|:--:|:-----|:----|:-------|
| 1 | Apply Paper 18 §IV.6 chemistry-side text block (this memo §4) | ~30 min mechanical | PROPOSED (do not apply autonomously per sprint mandate) |
| 1 | CLAUDE.md §3 row append: "M3 cyclotomic-mixed-Tate level-4 hypothesis for W1e correction terms — clean negative" | ~10 min | PROPOSED |
| 2 | Cross-system test: extract F4/F5/F6 analog correction terms at LiH max_n=3 (sub-percent W1c-residual at LiH); rerun W1e period-class diagnostic | 1-2 days | NAMED follow-on |
| 2 | Cross-system test on MgH₂ max_n=2 with the F5 explicit-core Hartree (already implemented infrastructure-wise; just needs F4-class sensitivity scan) | 2-3 days | NAMED follow-on |
| 3 | Construct chemistry-side analog of η-trivialization: identify a structural property of the FrozenCore $Z_{\rm eff}(r)$ pipeline that mirrors $\{\gamma_F, D_F\} = 0$ for the SM side | 2-4 weeks NCG-research | NAMED open question for Paper 18 §IV.6 |
| 4 | (DO NOT pursue) M3-engineering of chemistry corrections: any attempt to manually construct an M3-content correction (Catalan G times rationals) for W1e is wrong by structure | — | PROHIBITED |

---

## §7. Files

### Created (drivers)
- `debug/sprint_w1e_period_class.py` (~350 lines): PSLQ classification of 11 W1e correction terms against four bases (M1, M2, M3, INNER) at audit ceiling 100 + permissive ceiling 10^6 with three-class audit filtering.

### Created (data)
- `debug/data/sprint_w1e_period_class.json`: full numerical output, per-term PSLQ verdicts in each basis, audit-filter breakdown, verdict line.

### Created (memo)
- `debug/sprint_w1e_period_class_memo.md` (this memo).

### NOT modified
- Production `geovac/` modules — diagnostic-only sprint per sprint mandate.
- Tests — no production code modified, regression preserved.
- Papers (`papers/group3_foundations/paper_18_exchange_constants.tex` §IV.6 chemistry-side analog text) — drafted in §4 above, not applied per sprint mandate ("PROPOSED Paper 18 §IV.6 entry as text block in the memo. Do NOT edit .tex directly").
- CLAUDE.md — sprint mandate: documentation edits deferred.

---

**End of Sprint W1e period-class diagnostic memo. Verdict: STOP because 0/11 W1e correction terms identify with low-coefficient outer-factor periods (M1, M2, or M3) after audit filtering; period-class framing was the wrong axis — the W1e correction terms live in the inner-factor input-data tier (Paper 18 §IV.6 chemistry-side analog), categorically disjoint from outer-factor M1/M2/M3.**
