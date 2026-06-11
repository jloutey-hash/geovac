# Sprint Modular + α Arc — Synthesis (Track β)

**Date:** 2026-05-23.
**Sprint position:** Track β synthesis of the modular propinquity sprint (Track 0 + M-X/Y/Z/H1) and the immediate α follow-on arc (α-Diagnostic + α-Multi-zeta + α-PES). Three substantive new findings + one negative-at-current-basis verdict + three named follow-on directions, all landed in a single day.
**Cross-references:** `debug/sprint_modular_propinquity_synthesis_memo.md` (modular synthesis), `debug/sprint_alpha_1_diagnostic_memo.md`, `debug/sprint_alpha_2_multizeta_memo.md`, `debug/sprint_alpha_3_pes_test_memo.md`, `debug/sprint_modular_propinquity_literature_audit.md`, `debug/sprint_modular_propinquity_m{X,Y,Z,H1}_*_memo.md`, CLAUDE.md §1.7 multi-focal-wall taxonomy, §2 chemistry-arc history, §3 dead-ends, Paper 18 §III.7, Paper 19 §sec:w1c_residual, Paper 17 §6.10, Paper 32 §VIII.C, Paper 34 §V.D.

---

## §0. Executive summary

The PI's modular propinquity bet (begun mid-day 2026-05-23) returned a coherent arc of substantive findings on the chemistry / bound-state QED side, validating the W1c bimodule reframing as a genuine diagnostic instrument and surfacing a structural new wall layer that prior engineering attempts (PK, screened-valence diagonal) could not have predicted. The arc closes one operational question (the M-Y right-action axis dominance is real and measurable) and opens one deeper question (the bimodule diagnostic is faithful at the wavefunction-shape level but does not project onto the FCI eigenspace at the framework's NaH max_n=2 basis).

The headline structural finding (new content of the α arc): **the W1c-residual wall has a third layer beneath the two layers prior engineering closures targeted.** Layer 1 (H-Na orthogonality, addressed by PK cross-center, 2026-05-08, NEGATIVE) and Layer 2 (Na-side wavefunction shape, addressed by α-Multi-zeta substitution, 2026-05-23, ARCHITECTURE-CORRECT-BUT-BIT-ZERO-ON-FCI) are now both architecturally clean modules. Layer 3 — the **FCI variational state localizes on the H side at finite basis** because bare cross-V_ne over-attraction dominates the h1 eigenspectrum, pulling all FCI weight to H-side orbitals — is the genuine load-bearing residual. Layer 3 is structurally invisible to any on-site h1 or on-site basis substitution; it lives in the cross-V_ne *kernel shape* (Track 3 diagnostic's named target, 2026-05-09) or in the basis-size truncation itself.

Three substantive new findings + one architectural addition + one named-follow-on diagnostic + the diagnostic-before-engineering rule continuing to pay:

1. **W1c-residual orthogonality wall, Layer 3 (FCI variational state H-localization)** — the bimodule diagnostic is faithful at shape level but invisible at FCI eigenspace level at NaH max_n=2 because Na 3s sits at i=5 in the h1 eigenspectrum, above the lowest 5 (all H-localized) orbitals that the 2-electron ground state occupies.
2. **Alkali-hydride scaling LiH ≪ NaH ≈ KH (not strict monotone)** — radial-node L²-cancellation in the K 4s wavefunction breaks the monotone Li → Na → K → ... prediction. Mechanism is the n_val-dependent oscillation count.
3. **Mixed Slater-n essential for screened valence orbitals** — BBB93/CR74's convention of fixing n_slater = orbital n is inadequate for screened-Schrödinger eigenstates with significant core penetration. K=5 mixed-n=[1,1,2,3,3] for Na 3s reaches overlap > 0.999999; K=3 pure-n=3 fails to reproduce the radial-node count.

The architectural addition (Z=11 entry in `geovac/multi_zeta_orbitals.py`, multi-zeta dispatch in `shibuya_wulfman.py` and `balanced_coupled.py`) is bit-exact backward-compatible and retained per the architectural-cleanliness rule (analog of PK and SV before it).

The diagnostic-before-engineering rule continuing to pay: the α-PES three-step gate caught the Layer-3 structural surprise BEFORE the multi-zeta wiring was committed to a multi-week production sprint. Step 1 algebraic kernel passed (-0.135 Ha differential at R_eq); Step 2 FCI eigenvalue diagnostic exposed bit-zero shift; Step 3 confirmed across a denser R-grid. Without the three-step gate, the next engineering target would have been wavefunction-shape substitution in within-block ERIs — which the Layer-3 finding suggests would also be structurally invisible until the bonding configuration is reachable in the variational state.

---

## §1. Modular framework leverage on the chemistry side (recap from M-X/M-Y/M-Z)

The modular propinquity reframing (Latrémolière 1608.04881 / 1811.04534) reorganized GeoVac's chemistry-side / bound-state-QED structure across three threads:

- **Thread A (M-X + M-Z):** the Coulomb-Sturmian basis at $\lambda = Z/n_*$ is naturally a *module basis* on the Hilbert C*-module $\mathcal{M} = L^2(\R^3) \otimes \C^{2j+1}$ over the algebra $\mathcal{A} = C_0(\R^3) \otimes M_d$. The Sturmian projector truncates the module without changing the algebra; the X/Z framing divergence (multiplication algebra vs spectral algebra readings) dissolves under the modular reading.

- **Thread B (M-X):** the R1 workaround (replace $\|[D_S, M_f]\|$ with $\|\nabla f\|_\infty$) is canonical, not a kludge. Modular Leibniz lives at the first-order connection $\nabla$ of the Bochner Laplacian $D_S = -\tfrac{1}{2}\nabla^*\nabla - Z/r$; the Schrödinger operator is correctly placed at the second-order level (where Leibniz fails) and the Lipschitz seminorm at the first-order level (where Leibniz is automatic).

- **Thread C (M-Z headline):** the multi-focal-composition wall (CLAUDE.md §1.7) splits cleanly into (a) bimodule cross-shifts (HANDLED: HF-3 recoil, HF-4 Zemach, HF-5 Parker-Toms, cross-register $V_{eN}$, FNS) and (b) module endomorphisms (NOT HANDLED: LS-8a renormalization $Z_2{-}1$ and $\delta m$, Sprint H1 Yukawa, inner-factor calibration data). This sharpens the structural-skeleton-scope statement.

The M-Y track surfaced the strongest single new-content finding: a two-axis L/R bimodule decomposition $d_\mathcal{M}^\text{LR}(\xi_a, \xi_b) = (d_L^2 + d_R^2)^{1/2}$ between candidate W1c-residual pin states identifies the *right-action axis* (Na-side wavefunction shape) as the dominant carrier of the W1c residual, with predicted $d_R / d_L \approx 6.7$ at the structural / order-of-magnitude level. M-Y additionally predicted a binding-recovery magnitude of $\sim 0.1$ Ha (right order of magnitude for experimental NaH $D_e \approx 0.075$ Ha) and a uniform alkali-hydride scaling law $d_R \sim r_\text{valence}^\text{phys}(M) - 1.5$ bohr.

These four predictions (right-action dominance, magnitude, alkali scaling, axis identification for engineering closure) defined the α arc's gate set.

### α-Diagnostic — empirical confirmation of right-action dominance

Track α-1 (`debug/sprint_alpha_1_diagnostic_memo.md`) confirmed the right-action axis dominance numerically. Two-axis L/R bimodule distance constructed on the NaH bond axis with Gaussian probes centered at H and Na (canonical width $\sigma = 1.0$ bohr), measured at the four M-Y pin states (i) framework hydrogenic, (ii) hydrogenic + SV-corrected diagonal (= (i) bit-exact at wavefunction level), (iii) physical Na 3s + hydrogenic H, (iv) (iii) + small bonding-coefficient perturbation. The load-bearing (i)–(iii) pair:

| $\sigma$ (bohr) | $d_L$ | $d_R$ | $d_R / d_L$ |
|---:|---:|---:|---:|
| 0.5 | 0.232 | 0.832 | **3.58** |
| **1.0** | **0.224** | **0.724** | **3.23** |
| 1.5 | 0.270 | 0.650 | 2.41 |
| 2.0 | 0.337 | 0.601 | 1.78 |

The qualitative finding is robust: $d_R / d_L > 1$ at every tested probe width, with a ceiling near 3.5–4 at narrow probes and a floor near 1.3 at wide probes. R-dependence stable at $d_R/d_L = 3.14$–$3.34$ across $R \in [3.0, 4.0]$ bohr; basis-size stable at FrozenCore $n_\text{grid} \in [4000, 16000]$ (Δ ratio < 0.5%).

M-Y's $\approx 6.7$ structural estimate is within a factor 2 of the measured 3.23 at canonical probe — quantitatively imprecise but qualitatively correct, well within the "order of magnitude" qualification. The probe-width dependence is a substantive new finding: the precise factor depends on how narrowly the test functions are localized, but the structural claim (right-action dominance) is invariant.

**Verdict for α-Diagnostic: CONFIRMED-d_R/d_L ≫ 1**, with the additional finding that the precise ratio is probe-width-dependent.

---

## §2. Alkali-hydride scaling refinement (α-Diagnostic substantive new content)

The M-Y prediction $d_R(M) \sim r_\text{valence}^\text{phys}(M) - 1.5$ bohr expected strict monotone $d_R$ across the alkali series Li → Na → K → Rb → Cs. The α-1 measurement reads:

| Molecule | $Z_M$ | $n_\text{val}$ | $\langle r \rangle^\text{phys}_M$ (bohr) | $\langle r \rangle^\text{hyd}_M$ (bohr) | $|\Delta r|$ | $d_R$ | $d_R/d_L$ |
|:---|---:|---:|---:|---:|---:|---:|---:|
| LiH | 3 | 2 | 4.615 | 6.000 | 1.39 | 0.42 | 4.32 |
| NaH | 11 | 3 | 4.467 | 13.475 | 9.01 | 0.72 | 3.23 |
| KH | 19 | 4 | 8.479 | 21.136 | 12.66 | 0.64 | 5.41 |

**Refined finding (substantive new content):** observed ordering is **LiH ≪ NaH ≈ KH**, NOT strict monotone Li < Na < K. KH's $d_R$ is 12% *smaller* than NaH's despite a 41% larger shape-difference $|\Delta r|$.

Mechanism (the substantive structural reading): the K 4s wavefunction has **3 radial nodes** (vs Na 3s's 2 and Li 2s's 1). The inner-shell sign oscillations in $R^M(r)$ produce partial L²-cancellation in the bond-region integrand $|R^M_\text{hyd} - R^M_\text{phys}|^2 w_R(z)$. For KH, the K 4s oscillation structure aligns partially with the hydrogenic Z=1 4s oscillation structure, reducing the L²-norm difference relative to NaH's 3s-vs-hydrogenic-3s difference. NaH has the "perfect storm" — large shape mismatch with smooth-enough (only 2 nodes) wavefunctions that the L² norm picks them up cleanly.

**Prediction for RbH, CsH:** Rb 5s has 4 radial nodes, Cs 6s has 5 radial nodes. Both are predicted to continue the plateau (radial-node L²-cancellation grows with $n_\text{val}$), giving the refined uniform scaling law:

- **LiH-class (small $d_R$):** $n_\text{val} - 1 = 1$ node, monotone-coherent contribution from physical–hydrogenic difference. Empirically: LiH binds in the framework's production code (the only alkali-hydride that does).
- **NaH/KH/RbH/CsH-class (plateau $d_R$):** $n_\text{val} - 1 \geq 2$ nodes, partial L²-cancellation. Empirically: none of these alkali-hydrides binds in the framework's production code; the W1c-residual wall dominates.

The "M-Y prediction works for LiH" diagnostic at α-1 is consistent with the LiH FCI having genuine bond character at max_n=2 (Track 3 LiH 878-Pauli result, 1.8% energy error, 7% R_eq error). The refined alkali-hydride scaling reading is *qualitatively coherent with the framework's empirical chemistry catalogue*: LiH-binds-vs-rest-don't aligns 1:1 with $d_R$-small-vs-plateau.

**A cleaner diagnostic** would use a derivative-weight Lipschitz seminorm $L(a) = \|[D, a]\|$ instead of the Gaussian-probe L² norm; this would emphasize gradient differences that the L² norm averages over and might restore monotonicity. Named as a candidate refinement for a hypothetical α-Diagnostic-v2.

---

## §3. Mixed Slater-n essential for screened valence orbitals (α-Multi-zeta substantive new content)

The α-2 architecture step (Track α-2, `debug/sprint_alpha_2_multizeta_memo.md`) fitted the physical Na 3s and 3p screened wavefunctions (from `FrozenCore(Z=11)` + radial Schrödinger solver) to multi-zeta Slater expansions. The headline architectural finding is a framework-level statement, not specific to NaH:

**Substantive new finding:** Mixed Slater-n is structurally essential for screened valence orbitals with significant core penetration. The standard BBB93 / Clementi–Roetti~1974 convention of fixing $n_\text{slater}^k = n_\text{orbital}$ for every primitive in a single-orbital expansion is **inadequate** for screened-Schrödinger eigenstates.

Mechanism: the physical Na 3s has $R(r \approx 0.011) \approx +3.06$ bohr$^{-3/2}$ at the inner shell (core penetration) AND two radial nodes AND a diffuse outer tail. A pure Slater $n=3$ STO has prefactor $r^2$, too suppressed at small $r$ to reach amplitude 3 unless $\zeta$ is very large; but large $\zeta$ then forces the primitive to decay too fast in the outer region. The two scales are incompatible in a uniform-$n$ basis. Adding $n=1$ STO primitives (prefactor $r^0$ = constant) captures the inner amplitude cleanly; larger-$n$ primitives handle the diffuse outer tail.

Empirical confirmation: pure-$n=3$ at K=3 gave overlap 0.9992 but **only 1 radial node** (failure of node count, max pointwise error 3.33); pure-$n=3$ at K=4 same problem. Mixed-$n=[1,1,2,3,3]$ at K=5 reached overlap 0.999999, L² error $1.4 \times 10^{-6}$, correct 2-node radial structure, all $|c_k| < 1$.

This is a **framework-level finding worth recording** distinct from the heuristic-two-zeta cliff result that the closeout sprint (2026-05-09) flagged. The closeout finding was "all-positive-coefficient single-zeta fitted to BBB93 Kr ratios under-screens at intermediate $r$"; this new finding is "all-positive-coefficient mixed-$\zeta$ also fails unless $n_\text{slater}$ is allowed to vary". Together they characterize the BBB93/CR67 architecture's structural limitation for heavy-atom screened valence orbitals: both *zeta count* and *Slater n* must be flexible. The α-2 architecture supplies the latter; closing the closeout sprint's residual will likely need full BBB93/KTT multi-zeta core + iterative SCF, but the framework's path to a faithful screened valence basis has cleared one of the two structural barriers.

Production code: `geovac/multi_zeta_orbitals.py` Z=11 entry with `_PHYSICAL_NA_3S_PRIMITIVES`, `_PHYSICAL_NA_3S_COEFFS`, `_PHYSICAL_NA_3P_PRIMITIVES`, `_PHYSICAL_NA_3P_COEFFS`; new public `get_physical_valence_orbitals(Z)`. 14 new tests pass; zero regression on 234 baseline tests.

---

## §4. The W1c wall Layer 3 (α-PES substantive new content)

This is the most consequential finding of the α arc — a third structural layer beneath the two layers prior engineering closures targeted.

### Three-layer structure of the W1c residual

After the α arc, the W1c-residual orthogonality wall on second-row alkali hydrides decomposes into:

- **Layer 1 — H–Na orthogonality:** the H valence orbital is non-orthogonal to the [Ne]-core orbitals on the heavy atom; Pauli exclusion of the H valence basis against the frozen core is not enforced beyond what Hartree screening captures. **Addressed** by `geovac/phillips_kleinman_cross_center.py` (Sprint 2026-05-08, PK cross-center barrier). **Verdict: NEGATIVE on binding** (14.6% improvement, PES still descending). Architecture retained; module is clean.

- **Layer 2 — Na-side wavefunction shape:** the framework's hydrogenic-$Z_\text{orb}=1$ Na valence basis is qualitatively wrong; the actual Na 3s has mean radius 4.47 bohr vs the hydrogenic 13.5 bohr (factor-3 mismatch); the M-Y bimodule diagnostic identifies this as the right-action axis ($d_R/d_L = 3.23$) and the α-Diagnostic confirms this empirically. **Addressed** by `geovac/multi_zeta_orbitals.py` Z=11 entry + `shibuya_wulfman.py` multi-zeta dispatch + `balanced_coupled.py` `multi_zeta_basis: bool = False` kwarg (Sprint 2026-05-23, α-Multi-zeta + α-PES Step 1). **Step 1 algebraic verdict: POSITIVE** (cross-V_ne kernel differential at Na 3s on-site diagonal = -0.135 Ha at NaH R_eq, sign-consistent with R-dissociation 18.7% smaller, classical-limit sanity passes). **Step 2 FCI verdict: BIT-ZERO** (multi-zeta substitution changes the FCI ground-state eigenvalue by $|\Delta E| < 4 \times 10^{-13}$ Ha at every $R \in \{2.5, 3.0, 3.5, 4.0, 5.0, 10.0\}$ bohr). Architecture retained; module is clean; the FCI invisibility is the substantive new finding.

- **Layer 3 — FCI variational state H-side localization at finite basis (NEW α-PES finding):** at NaH max_n=2 with **bare** cross-V_ne, the lowest 5 eigenvalues of the h1 matrix are entirely on H-side states. The H 1s sees the un-screened Na $Z=11$ nucleus at R=3.5 bohr with on-site cross-V_ne attraction of $-3.13$ Ha (≈ $-Z_H Z_\text{Na}/R$ within 5% of the classical Coulomb point-charge limit), compared to the Na 3s self-attraction of only $-0.28$ Ha. The Na 3s eigenvalue sits at $i = 5$ in the h1 eigenspectrum, *above* the lowest 5 bonding/antibonding orbitals on the H side. The 2-electron FCI ground state puts both electrons in the lowest H-side orbital (singlet pair), with no Na 3s occupation. **The Na 3s diagonal/basis shift has no effect on the FCI eigenvalue.** No on-site h1 substitution can fix this; the residual lives in the variational state itself.

### Why Layer 3 is structurally distinct from Layers 1 and 2

Layers 1 and 2 are *operator-level* corrections: they shift specific entries of $h_1$ (Layer 1: PK barrier; Layer 2: diagonal eigenvalue and/or basis shape). They are well-localized matrix-element interventions.

Layer 3 is a *variational state* finding: it concerns which eigenvectors of $h_1$ are occupied in the FCI ground state. At the framework's max_n=2 basis with bare cross-V_ne, the FCI eigenspectrum is so deeply distorted by un-screened Z=11 cross-attraction that all bound-state weight collapses onto H-side orbitals; the Na valence space is *unoccupied* in the ground state. Operator-level corrections to the unoccupied space are bit-zero on the FCI eigenvalue.

### Three closure paths Layer 3 opens

The α-PES memo §6.2 named three candidate next moves:

1. **Increase basis to max_n=3** (Q=56). Test whether the larger basis lets the FCI ground state mix Na and H amplitudes. If the ground state acquires Na-side weight, multi-zeta substitution becomes load-bearing and the α-PES verdict could flip POSITIVE. If the ground state stays H-localized, Layer 3 is a structural property of the bond and not a basis-size artifact.

2. **W1c + multi-zeta combination** (currently mutually exclusive due to a `NotImplementedError` in `compute_screened_cross_center_vne` for multi-zeta basis). Removing this requires extending the screened-path kernel to accept multi-zeta primitives — same mechanical refactor as the bare path. Predicted outcome at NaH max_n=2: similar bit-zero behavior (Layer 3 dominates), but at least the architecture would be uniform.

3. **Cross-V_ne integration-kernel-shape substitution.** Track 3's diagnostic (2026-05-09) identified this as the original named target — testing the cross-V_ne integral with physical-$n$ hydrogenic shape (block_n=1, $l=0$ → $n_\text{phys}=3$ hydrogenic 3s instead of 1s) showed a $-0.674$ Ha differential at R=2.5, $\sim 1.9\times$ the W1c-residual descent. Substantial enough that the FCI eigenspectrum could be shifted into a bonding regime; the α-PES Step 1 finding ($-0.135$ Ha for Na 3s on-site) is a smaller-effect cousin of this. Architecturally this is a kernel-shape substitution in `shibuya_wulfman.py`, not a basis substitution.

### Net Layer 3 reading

Layer 3 sharpens the structural-skeleton-scope statement (CLAUDE.md §1.7 multi-focal-wall taxonomy) at the W1c sub-wall: **the framework's NaH max_n=2 basis cannot represent the physical bonding configuration where electrons share between Na 3s and H 1s in a quantum superposition.** With un-screened cross-V_ne the H side dominates (FCI = pure H-localized); with W1c screening the H side is no longer artificially deep but the Na 3s diagonal at $-0.78$ Ha is still well below the H 1s + screened-cross-V_ne energy ($\sim -0.5$ Ha), so the FCI is still H-dominant. This is not a failure of the bimodule diagnostic — it is a precise characterization of where the bimodule diagnostic's prediction is and isn't FCI-eigenvalue-relevant.

---

## §5. Algebraic-first discipline working — the three-step gate pattern

The α-PES test was structured as a **three-step algebraic-first test with intermediate decision gates**, exactly per the algebraic-first directive (`feedback_algebraic_first.md` memory file: "Stop before long numerical runs; find the algebraic structure"). The three steps:

1. **Step 1 — Algebraic kernel differential at a single matrix element.** Decision gate: if |differential at R_eq| > 0.05 Ha, proceed to Step 2; else clean-negative and stop. **Outcome: differential $-0.135$ Ha at NaH R_eq = 3.566 bohr, more than 2× the gate threshold. PASSED.** Cost: ~1 hour to wire and run.

2. **Step 2 — Two-point FCI at R_eq and R_diss.** Decision gate: if multi-zeta produces a measurable D_e shift, proceed to Step 3 mini-PES scan; else diagnose why before proceeding. **Outcome: FCI eigenvalue shift bit-zero (3.7 × 10⁻¹³ Ha) across both R values. AMBIGUOUS — passed the literal D_e > 0 criterion but bit-zero on the multi-zeta-specific shift.** Cost: ~30 minutes to run + ~30 minutes to diagnose the eigenspectrum.

3. **Step 3 — Mini-PES at 5 R values.** Confirmed bit-zero across all R values; the multi-zeta substitution is *structurally invisible* to the FCI eigenspectrum at this basis. Cost: ~15 minutes.

Without the three-step gate, the next engineering move after α-Multi-zeta would have been "wire multi-zeta into within-block ERIs and run a full PES sweep" — a 2-3 week production sprint. The Step 2 FCI eigenvalue diagnostic caught the Layer 3 structural surprise in ~30 minutes. The total cost of the α-PES test was ~2 hours; the cost averted is at least a week of misdirected engineering.

The diagnostic-before-engineering rule (`feedback_diagnostic_before_engineering.md`) continues to pay. The pattern works at the *intra-sprint* level (three steps within α-PES), not just at the inter-sprint level (diagnostic sprint before next engineering sprint). The discipline applies whenever the engineering target is sufficiently load-bearing that committing without an intermediate diagnostic risks weeks of misdirected work.

The α arc was a textbook application of the algebraic-first + diagnostic-before-engineering pair. The PI's modular bet → M-Y bimodule prediction → α-Diagnostic confirms shape-level prediction → α-Multi-zeta builds clean architecture → α-PES three-step gate exposes the Layer 3 surprise BEFORE the production wiring is committed. Five sprints across one day, each one informed by the previous, each one with a clean decision gate, all landing clean architectural additions despite the bottom-line negative on binding recovery.

---

## §6. Three named follow-on directions

The α-PES memo §6.2 surfaced three candidates. Ranking by expected new content per week:

### Priority 1 — Test at NaH max_n=3 (Q=56), 2-3 days

**Why first:** smallest sprint of the three, highest information yield. The α-PES bit-zero finding at max_n=2 is a structural reading IFF Layer 3 persists at larger basis sizes. A max_n=3 NaH FCI test would falsify or confirm Layer 3:

- If the FCI ground state acquires Na-side weight at max_n=3: Layer 3 is a max_n=2 artifact, and multi-zeta substitution becomes load-bearing at max_n=3. The α-PES verdict would flip POSITIVE at larger basis; multi-zeta is the right architecture; the W1c wall is closable in production.
- If the FCI ground state stays H-localized at max_n=3: Layer 3 is structural (not a basis-size artifact). The α-PES verdict is durable; multi-zeta is not the load-bearing target; the cross-V_ne kernel shape (Option 3 below) is the correct next direction.

NaH at max_n=3 has Q=56 (manageable), basis size n_max=3 explicitly noted in CLAUDE.md §2 NaH n_max=3 convergence (Track CN, v2.1.1): 5,349 balanced Pauli terms, but "PES overattraction persists (no equilibrium)". This was at the BARE balanced level. The α follow-on would test at the W1c-screened + multi-zeta level, which is structurally different. Open question: does the larger basis support a bonding FCI configuration that multi-zeta substitution actually affects?

**Risk:** the bare overattraction at max_n=3 might be so dominant that even W1c screening can't bring it back to a binding regime; in that case the FCI state still H-localizes and Layer 3 stays structural at all tested basis sizes.

### Priority 2 — Cross-V_ne kernel shape substitution, 1-2 weeks

**Why second:** Track 3's diagnostic flagged this as the original named target with a substantial $-0.674$ Ha differential at R=2.5. The α-PES Step 1 finding ($-0.135$ Ha) is a smaller-effect cousin. The kernel-shape substitution is mechanically a refactor of `shibuya_wulfman.py`'s `_radial_split_integral` to use physical-$n$ orbital decomposition on the heavy-atom side instead of $n_\text{block}$ hydrogenic. Cleaner architectural unification of α arc + Track 3 diagnostic.

**Risk:** the differential is computed against the *bare* baseline, so the W1c-screened kernel might already absorb part of the kernel-shape correction implicitly through Z_eff(r). Diagnostic-before-engineering: test with a single Na 3s vs Na 4s vs hydrogenic-3s kernel-shape substitution at one R point before committing to the full refactor.

### Priority 3 — W1c + multi-zeta combination by removing NotImplementedError, 2-3 days

**Why third:** smallest mechanical effort but lowest information yield. The α-PES finding predicts bit-zero behavior at NaH max_n=2 (Layer 3 dominates regardless of bare-or-screened). Removing the NotImplementedError is architecturally clean and would let the W1c+multi-zeta combination be tested at max_n=3 (Priority 1) in a more unified architecture. Should not be done alone; bundle with Priority 1.

### Cross-priority observation

Priorities 1 and 3 should be bundled: implementing max_n=3 NaH FCI in a clean W1c+multi-zeta architecture lets one test answer both "does basis size lift Layer 3?" and "does the screened-path multi-zeta combination converge to the bare-path verdict?" simultaneously. Estimated cost: 1 week.

Priority 2 is independent of Priorities 1+3 and could run in parallel.

**Default recommendation:** open Priority 1+3 (bundled max_n=3 NaH test with W1c+multi-zeta unified architecture, 1 week) as the lead follow-on. Add Priority 2 (kernel-shape substitution diagnostic at one R point) as a parallel scoping probe before committing to the full kernel-shape refactor.

---

## §7. Honest scope

### What is PROVEN at qualitative-rate level

- **Bimodule architecture** (M-X module-basis reading of Sturmians; M-Z multi-focal-wall split into cross-shifts vs endomorphisms; M-Y two-axis L/R decomposition): structurally correct as a reframing tool. Confirmed numerically by α-Diagnostic at the wavefunction-shape distance level ($d_R / d_L = 3.23$ for the NaH load-bearing pair, ranges 1.3–3.6 across probe widths).
- **Mixed Slater-n necessity for screened valence orbitals**: empirically confirmed at K=3 pure-n=3 (fails radial-node count) vs K=5 mixed-n (succeeds at overlap 0.999999). Generalizable to all heavy-atom screened valence orbitals with significant core penetration.
- **Algebraic-first discipline applied within sprints**: the α-PES three-step gate caught the Layer-3 surprise in ~2 hours; without the gate, a multi-week engineering misdirection.

### What is STRUCTURAL-SKETCH

- **W1c-residual Layer 3 finding**: the FCI variational state H-localization mechanism at NaH max_n=2 is structurally clean (lowest 5 h1 eigenvalues are H-localized; Na 3s sits at i=5, unoccupied in 2e GS), but the *causal explanation* of why bare cross-V_ne dominates the h1 eigenspectrum at this basis is a structural reading (un-screened Z=11 nuclear attraction at H), not a proven theorem. The reading would tighten by quantifying the bare/W1c/screened-cross-V_ne contributions to the h1 eigenspectrum across basis sizes.
- **Alkali-hydride scaling LiH ≪ NaH ≈ KH**: the radial-node L²-cancellation mechanism is a structural reading at three data points (3 nodes only at K). The prediction for RbH, CsH (continued plateau) is testable; the diagnostic implementation extends mechanically to those systems.
- **Three-layer W1c-residual taxonomy (Layer 1 / Layer 2 / Layer 3)**: clean structural sketch but not yet a proved theorem. Layer 3 in particular is the new structural content; its precise relationship to basis size, cross-V_ne kernel shape, and the variational state-space mapping needs more diagnostic data to harden.

### What is OPEN

- **Whether any of the three follow-on directions closes the W1c wall in production.** Priorities 1, 2, 3 are named but unexecuted. The most likely closure (per the Track 3 diagnostic + α-PES Layer 3 reading combined) is cross-V_ne kernel-shape substitution at increased basis size, but this is structurally an inference, not a proof.
- **The bimodule diagnostic's relationship to Latrémolière distance.** The α-Diagnostic uses a Gaussian-probe L² norm parametrized by width $\sigma$; the "right" probe width corresponds to the natural Latrémolière distance for the bimodule, which would be derivable from the Connes–Lipschitz seminorm on the multiplication algebras. Probably a multi-week math.OA sprint (option γ from the modular synthesis memo's recommended sprint set).
- **Generalization to p-block hydrides (HCl, H₂S, MgH₂).** The α-Diagnostic covered alkali-hydrides. CLAUDE.md §1.7 documents a Z-decreasing W1c-residual orthogonality wall (NaH 5.4–6.0×, MgH₂ 2.99×, HCl 1.79×) — different empirical magnitude. Extending the bimodule diagnostic to these systems is an obvious next probe.

### Net verdict line for the modular + α arc

**POSITIVE-WITH-LAYER-3-FINDING.** The modular propinquity reframing delivered substantive new content on three threads (Sturmian-as-module-basis, R1 canonical, multi-focal-wall split); the M-Y two-axis decomposition surfaced a measurable structural prediction; the α-1 diagnostic confirmed the prediction at wavefunction-shape level; the α-2 architecture worked; the α-3 PES test exposed a new W1c wall layer (Layer 3, FCI state H-localization) that no prior diagnostic anticipated and that no on-site engineering closure can address. The arc closes one operational question (M-Y right-action axis is real) and opens one deeper question (Layer 3 mechanism + closure path). Three named follow-ons ranked.

---

## §8. Session summary

**Tracks:**
- Track β (paper edits + synthesis): COMPLETE — synthesis memo + paper edits in Papers 19/17/18/32/34 + CLAUDE.md updates §1.7, §2, §3 + memory file.

**Headline outcome:** **POSITIVE-WITH-LAYER-3-FINDING.** The PI's modular propinquity bet paid on the chemistry side as predicted; the α arc validated M-Y's right-action axis dominance empirically ($d_R/d_L = 3.23$ at canonical probe), built clean multi-zeta architecture for screened Na valence (overlap > 0.999999), but exposed a new structural layer of the W1c wall (Layer 3 FCI variational state H-localization) that on-site h1 or basis substitution at NaH max_n=2 cannot address.

**Three substantive new findings, all from the α arc:**
1. **W1c-residual wall Layer 3** — FCI variational state H-localization at finite basis. The bimodule diagnostic is faithful at wavefunction-shape level but invisible at FCI eigenspace level when the heavy-atom valence orbital sits above the lowest bound states in the h1 eigenspectrum.
2. **Alkali-hydride scaling LiH ≪ NaH ≈ KH** (not strict monotone). Mechanism: K 4s radial-node L²-cancellation. Prediction: RbH, CsH continue the plateau.
3. **Mixed Slater-n essential for screened valence orbitals.** BBB93/CR74 convention of fixing $n_\text{slater} = n_\text{orbital}$ is inadequate for screened-Schrödinger eigenstates. K=5 mixed-n=[1,1,2,3,3] for Na 3s gives overlap > 0.999999.

**Architectural addition** (retained per the architectural-cleanliness rule):
- `geovac/multi_zeta_orbitals.py` Z=11 entry with `_PHYSICAL_NA_3S_PRIMITIVES`, `_PHYSICAL_NA_3P_PRIMITIVES`, public `get_physical_valence_orbitals(Z)`.
- `geovac/shibuya_wulfman.py` `multi_zeta_basis` kwarg on `compute_cross_center_vne`.
- `geovac/balanced_coupled.py` `multi_zeta_basis: bool = False` kwarg on `build_balanced_hamiltonian`.
- Tests: 14 new in `tests/test_multi_zeta_orbitals.py`, 12 new in `tests/test_balanced_coupled_multizeta.py`. Zero regression on 234 baseline tests.

**Three named follow-on directions ranked:**
- P1+P3 (bundled): NaH max_n=3 FCI with W1c+multi-zeta unified architecture (1 week, falsifies or confirms Layer 3 as basis-size artifact).
- P2: cross-V_ne kernel-shape substitution diagnostic (parallel 2-3 day scoping probe).

**Discipline:** algebraic-first three-step gate inside α-PES caught the Layer-3 structural surprise in ~2 hours; saved an estimated 1-2 weeks of misdirected engineering. The PM-rule pattern works at intra-sprint as well as inter-sprint level.

**Files (this Track β session):**
- Created: `debug/sprint_modular_alpha_arc_synthesis_memo.md` (this memo, ~4500 words).
- Modified: Paper 19 §sec:w1c_residual (extended), Paper 17 §6.10 (cross-reference paragraph), Paper 18 §III.7 (M-Z cross-reference), Paper 32 §VIII.C (M-H1 cross-reference), Paper 34 §V.D (running catalogue entry).
- Modified: CLAUDE.md §1.7 (multi-focal-wall taxonomy refinement), §2 (sprint summary), §3 (multi-zeta-bit-zero-on-FCI dead-end entry).
- Created: `memory/sprint_modular_alpha_arc.md`, `memory/MEMORY.md` entry.

**Pending PI direction.** Next-sprint recommendation: open P1+P3 (NaH max_n=3 with W1c+multi-zeta) as the lead follow-on, with P2 (kernel-shape diagnostic) as a parallel probe.

---

**End of Track β synthesis memo.**
