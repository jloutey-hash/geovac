# Sprint L3e-P3 synthesis — three sub-sprints landed, net verdict

**Date:** 2026-05-23.
**Sprint position:** Synthesis of Sub-sprints X (Sturmian-as-Latrémolière verification), Y (W1c chemistry pin-state diagnostic), Z (Bethe log Latrémolière interpretation). All three from the Phase A.2' deep-read of Latrémolière arXiv:2512.03573.

**Outcome:** **MIXED across all three. The deep-read memo's structural identification has a key choice of $\Acal$ that determines whether the identification is rich or trivial. Two framings emerged — both legitimate, with different physics implications.** Net Latrémolière framework gives the GeoVac framework rigorous math.OA vocabulary for its Sturmian-FCI machinery, but does NOT autonomously close any open physics walls. Net consistent with the structural-skeleton-scope statement (CLAUDE.md §1.7) — the framework gives skeleton structure, calibration data (and now: convergence-rate constants, pin-state choices, renormalization counterterms) remains external input.

---

## §1. The X/Z divergence — choice of C*-algebra $\Acal$

The deep-read memo §2 claimed: "GeoVac's Sturmian truncation at $n_{\max}$ is an L-Lipschitz μ-pinned exhaustive sequence in Latrémolière's Def 1.29." Sub-sprints X and Z both verified this claim, but **at different choices of $\Acal$** that give structurally different framings.

### X reading (Coulomb multiplication algebra)

$\Acal = C_0(\R^3) \otimes M_d(\C)$ (multiplication algebra on $\R^3$, possibly with spinor structure). $L(f) = \|[D, M_f]\|$ for $D$ = Dirac (Schrödinger fails Leibniz; needs R1 gradient workaround).

Under this reading:
- **The Sturmian projector $P_{n_{\max}}^{\mathrm{Sturmian}}$ is NOT in $\Acal$** (it's a projector on Hilbert space, not a multiplication operator).
- **The exhaustive sequence is smooth radial cutoffs $\chi_n \in C_c^\infty(\R^3)$**, and the three Def 1.29 axioms are non-trivially verified (X memo §4).
- **The Sturmian truncation lives at the operator-system-truncation level** (tunnel quotient morphism in Latrémolière §2.1 sense).

### Z reading (spectral algebra)

$\Acal$ includes projectors on Hilbert space (e.g., a $C^*$-subalgebra of $\mathcal{B}(\mathcal{H})$ containing Sturmian projectors). $L(T) = \|[D, T]\|$ for $T$ in this algebra.

Under this reading:
- **The Sturmian projector $P_{n_{\max}}^{\mathrm{Sturmian}}$ IS in $\Acal$.**
- **The three Def 1.29 axioms hold but two are TRIVIAL:**
  - (ii) $\mu(h_n) = 1$ identically (ground state is in Sturmian span for any $n_{\max} \ge 1$, by the Sturmian saturation property at $\lambda = Z/n_{\mathrm{state}}$)
  - (iii) $\|h_n\| = 1$ identically (projector operator norm)
  - Only (i) $L(h_n) \to 0$ is non-trivial.
- **Leibniz holds for first-order Dirac** (Z's bound-state QED scope is exactly this).
- **The Sturmian truncation IS the exhaustive sequence** (in the spectral framing).

### Are both framings valid?

YES — they are different but compatible structural framings of the same physical content. The choice between them depends on what the framework is trying to extract:

- **X framing (Coulomb-mult algebra) is structurally richer** — it uses the natural algebra, gives a non-trivial exhaustive sequence, places Sturmian truncations at the operator-system-truncation level (where propinquity bookkeeping happens), and exposes the Leibniz axiom failure for Schrödinger.
- **Z framing (spectral algebra) is structurally trivial-at-two-axioms but applies cleanly to Dirac** — the Sturmian saturation property makes two axioms automatic, leaving only the L-decay rate as the substantive question.

**Both framings give Sprint Z's load-bearing finding:** Paper 38's $4/\pi$ rate does NOT transport to the non-compact Coulomb case. The Sturmian convergence rate is $\sim 1/N^4$ empirically (Paper 36 LS-3 acceleration form, 1S contributions: 8.62%/3.13%/0.60% at N=12/16/20), which is a categorically different rate-class than Paper 38's compact SU(2) Plancherel-weighted $\log n/n$ asymptote.

**Synthesis verdict:** the framework gains a math.OA-grounded error model for Sturmian truncations under EITHER reading, but **the explicit convergence rate is not determined by either Latrémolière 2512.03573 or by analogy to Paper 38**. The rate is observable empirically (Paper 36 LS-3 sub-percent precision at finite N) but Latrémolière's framework only proves existence-of-decay, not a specific rate.

## §2. The Y verdict — pin-state choice as Latrémolière diagnostic

Sub-sprint Y reframed the W1c-residual NaH binding wall via Latrémolière pin-state-distance diagnostic. Four candidate pin states:
- (i) Hydrogenic Z=1 on Na valence (default — gives the wall)
- (ii) SV-corrected diagonal eigenvalue (Track 3 — bit-exact at wavefunction level to (i), local metametric distance = 0)
- (iii) Cross-center bonding-orbital pin state (THE RIGHT CANDIDATE)
- (iv) MP2-perturbed pin state (smaller effect than (iii))

Structural ordering at bond-region radii: $\delta_r((i), (iii)) > \delta_r((i), (iv)) > \delta_r((i), (ii)) = 0$.

**Headline:** the Latrémolière reframing **CONFIRMS Track 3's named target at a math.OA-grounded level** (the genuine residual lives in cross-V_ne integration shape — specifically the bonding-orbital pin state — not in diagonal eigenvalue substitution). The reframing **does NOT autonomously close the wall**; the implementation work (Option 1 physical-n hydrogenic relabeling, or Option 2 strict-Schmidt) remains as named in Track 3.

**Net Y contribution:** structural vocabulary upgrade from "diagnostic intuition" to "quantitative pin-state distance ordering". Real new content but bounded by the structural-skeleton-scope statement.

**Recommendation:** open Y.1 (compute explicit $\delta_r$ values for the four candidates at bond-region radii, ~1-2 weeks) AFTER Sub-sprint X lands (X's Sturmian-as-Latrémolière framing is a prerequisite for Y's pin-state-distance framing).

## §3. The Z verdict — Bethe log / Lamb shift Latrémolière interpretation

Sub-sprint Z interpreted Paper 36's bound-state QED Lamb shift result through the Latrémolière framework. The Sturmian truncation in Paper 36's LS-3 acceleration form maps to the spectral-framing Def 1.29 exhaustive sequence (Z's reading from §1).

**What the framework GIVES Paper 36:**
- Rigorous math.OA error bar replacing the empirical convergence statement (Sturmian convergence at $N = 12, 16, 20$)
- Uniform convergence vocabulary across the atomic precision catalogue (H, He, Li, Be, ...)
- Cleaner scope statement for the LS-7/LS-8a residual decomposition (~+1.20 MHz multi-loop QED + ~+4.4 MHz non-loop physics — Latrémolière now provides the "non-loop physics" envelope as the operator-system-content beyond the pin-state local metametric)

**What the framework does NOT give Paper 36:**
- LS-3 acceleration-vs-velocity 3.3× speedup (this is a spectral-density-reweighting choice, not in the framework's vocabulary)
- LS-4 Drake-Swainson regularization (external regularization scheme, not a Latrémolière construction)
- LS-7 two-loop bracket $C_{2S} = +3.63$ (operator-system content beyond the local metametric)
- LS-8a Z_2/δm renormalization counterterms (the renormalization wall confirmed by Sprint H1 / LS-8a renormalization gap, see CLAUDE.md §2)

**Net Z contribution:** the framework gives Paper 36 a math.OA-grounded error model, but the LS-8a renormalization gap is reframed (operator-system-content language) rather than closed. Respects the structural-skeleton-scope statement.

**Recommendation:** open Z.1 (Paper 36 §IX or §X vocabulary update, ~10 working days = 2 weeks) AFTER Sub-sprint X lands.

## §4. Cross-cutting net findings

Three findings consistent across all three sub-sprints:

### Finding 1: the Latrémolière 2512.03573 framework applies to GeoVac's Sturmian-FCI machinery

Under either the X or Z reading of $\Acal$, the framework's pinned-QLCMS structure (Def 1.26) and pinned propinquity (Def 4.X) provide a rigorous math.OA model for:
- Sturmian truncation convergence (atomic FCI, Hylleraas r₁₂, Hylleraas-Eckart)
- Cross-system comparison (different atoms, different geometries)
- Error bounds on bound-state QED observables (Paper 36 Lamb shift)
- Pin-state-shape diagnostic (Paper 17 W1c residual NaH binding)

This is **substantively new content** for the framework — none of the prior math.OA infrastructure (Papers 38, 39, 40, 42, 43, 44, 45, 46) provides a tool for these non-compact / open-system / pin-state-dependent contexts.

### Finding 2: the framework does NOT autonomously close any open physics wall

Across the three sub-sprints, none of the framework results closes an open physics wall:
- X corrects the structural identification, doesn't close any wall
- Y confirms Track 3's named target structurally, doesn't autonomously close the NaH binding wall
- Z reframes the LS-7/LS-8a residual decomposition, doesn't close the renormalization gap

**This is consistent with the structural-skeleton-scope statement** (CLAUDE.md §1.7): the framework provides math.OA vocabulary and structural understanding, but calibration data (rate constants, pin-state choices, renormalization counterterms) remains external input.

### Finding 3: Schrödinger Leibniz failure requires R1 gradient workaround

Sub-sprint X surfaced a load-bearing structural point: $L(f) = \|[D_\mathrm{S}, M_f]\|$ for second-order Schrödinger $D_\mathrm{S}$ FAILS the Leibniz axiom (Latrémolière Def 1.18). The cross-term $2 \nabla f \cdot \nabla g$ in $[\nabla^2, M_{fg}]$ doesn't cancel.

**R1 workaround:** use $L(f) = \|\nabla M_f\| = \|\nabla f\|_\infty$ (first-order gradient operator) as the Lipschitz seminorm. This is standard NCG (Connes-Marcolli spin-Dirac approach) and satisfies Leibniz cleanly.

**Implication:** the framework's Schrödinger-based atomic FCI calculations (`casimir_ci`, `hylleraas_r12`, `hylleraas_eckart_pstate`) inherit Latrémolière hypertopology ONLY via the R1 gradient-Dirac substitution, not via the direct Schrödinger commutator. The framework's Dirac-based calculations (Tier 2 spinor composed, Paper 36 Lamb shift) inherit it directly.

This is a structural finding that **applies uniformly to all GeoVac Schrödinger-based work** — not just to the three sub-sprints. It should be added to Paper 18 §III.7 (master Mellin engine) or to a new paper that treats the gradient-Dirac structure explicitly.

## §5. Comparison to deep-read memo's claims

Deep-read memo §2-§4 made three physics-application claims. Status after Sub-sprints X+Y+Z:

| Deep-read claim | Status | Note |
|:----------------|:-------|:-----|
| **#1 Sturmian-FCI inherits Latrémolière convergence model** | PARTIAL — structural framing has two valid readings; explicit rate is not determined by Latrémolière | X+Z findings |
| **#2 W1c chemistry pin-state has Latrémolière diagnostic** | CONFIRMED-STRUCTURAL — pin-state-distance ordering confirms Track 3 named target, doesn't autonomously close | Y finding |
| **#3 Bethe log / Lamb shift have Latrémolière error model** | CONFIRMED — rigorous error bar replaces empirical convergence statement, LS-8a renormalization gap reframed but not closed | Z finding |

**Net deep-read assessment after verification:** the deep-read memo's identification was structurally on the right track but mapped Sturmian truncations to the WRONG Latrémolière slot under one of the two reasonable choices of $\Acal$ (the X reading). Both X (structurally rich, Coulomb-mult algebra) and Z (structurally trivial-at-two-axioms, spectral algebra) readings give legitimate identifications. The framework gives the GeoVac framework new math.OA vocabulary for Sturmian-FCI machinery under either reading.

## §6. Recommendations for next sprints

Three sprint options surfaced from the synthesis:

### Option A: Open Sub-sprints X.1, Y.1, Z.1 in parallel — ~3-6 weeks total

- **X.1**: Compute explicit Latrémolière propinquity bound for He 1¹S Hylleraas-Eckart sequence (ω=2,3,4 at empirical 3.6%, 0.04%, 0.0006%). Identify the explicit Lipschitz constant and the asymptotic rate. ~1-2 weeks.
- **Y.1**: Compute explicit local metametric $\delta_r$ values for the four NaH pin-state candidates. Verify the structural ordering, confirm Track 3 named target. ~1-2 weeks.
- **Z.1**: Paper 36 §IX or §X vocabulary update with Latrémolière error model, LS-8a scope refinement. ~10 working days = 2 weeks.

These three can run in parallel. Net deliverable: three additions to the math.OA arc + one paper update.

### Option B: Open Sub-sprint X.1 only, then sequential Y.1 then Z.1

If sequential pacing is preferred, X.1 should land first (prerequisite for both Y.1 and Z.1), then Y.1 and Z.1 can run in parallel afterward.

### Option C: Defer all three and return to math.OA Krein-lift Phase A.2'

The Phase A.2'-A.4' Krein-lift sprint (re-scoped, `debug/sprint_l3e_p3_rescope_memo.md`, ~10 weeks, single merged Paper 48) is still on the program. If the physics-side applications don't justify the next sprint cycle, return to the Krein-lift work.

**Sub-agent recommendation:** Option A or B. The framework gains real math.OA-vocabulary upgrades for Sturmian-FCI machinery, pin-state-shape diagnostic, and bound-state QED error model. Each of the three sub-sprints landed substantively (MIXED verdicts, not negatives). Net program progression: from "Latrémolière 2512.03573 = scoop of Phase B" (Phase A.1 audit finding) to "Latrémolière 2512.03573 supplies vocabulary for three concrete physics-side applications" (Phase A.2' deep-read).

## §7. Open structural questions surfaced

Three structural questions remain open after the three sub-sprints:

### Q1: What is the right $\Acal$ for the Sturmian framework?

X uses $\Acal = C_0(\R^3) \otimes M_d$ (Coulomb-mult algebra); Z uses $\Acal$ containing projectors (spectral algebra). Both give valid framings but with different structural implications:
- X is structurally richer (non-trivial exhaustive sequence)
- Z is structurally cleaner (Sturmian = exhaustive sequence directly, but two of three axioms trivial)

The framework's natural choice depends on what's being extracted. If the focus is on propinquity bookkeeping (operator-system truncation level), X is right. If the focus is on direct exhaustion-rate calculation, Z is right.

**Recommendation:** keep both framings explicitly in the framework's vocabulary; use X for cross-system propinquity work, use Z for explicit Sturmian-FCI error model.

### Q2: Why doesn't Paper 38's $4/\pi$ rate transport to non-compact Coulomb?

Sub-sprint Z surfaced this clean negative: SU(2) Plancherel-weighted compact rate ($4/\pi \cdot \log n/n$ asymptote, Paper 38) is categorically different from Sturmian's empirical $1/N^4$ rate (Paper 36 LS-3).

**Structural reading:** the framework's compact-side master theorem (Paper 40, $4/\pi$ universal across all simple compact Lie groups) is genuinely compact-specific. Non-compact extensions need fresh rate analysis (which Latrémolière 2512.03573 does NOT provide — it proves existence-of-decay, not rate).

This is a NEW open question for the math.OA arc — what is the non-compact analog of the $4/\pi$ universality? Could it be $1/N^4$ (Sturmian-specific) or something more general? Sub-sprint X.1 would test this empirically for He 1¹S Hylleraas-Eckart.

### Q3: Does the R1 gradient-Dirac workaround give a uniform Latrémolière framework for both Schrödinger and Dirac GeoVac work?

R1 says: use $L(f) = \|\nabla f\|_\infty$ instead of $\|[D_\mathrm{S}, M_f]\|$ for Schrödinger-based work. This gives Leibniz cleanly and matches the Connes-Marcolli spin-Dirac standard.

**Open question:** does the R1 framework give the same physical predictions as the direct Schrödinger commutator framework? Probably yes for atomic FCI calculations (where the natural seminorm is the gradient anyway), but should be verified explicitly.

**Recommendation:** add R1 verification to Sub-sprint X.1 (~1 day of work).

## §8. Honest scope

This synthesis:
- IS a synthesis of three independent sub-sprint verdicts
- IS NOT a new structural result (the structural identifications and findings are from X, Y, Z individually)
- DOES surface a load-bearing X/Z divergence (choice of $\Acal$) that the deep-read memo missed
- DOES propose three sprint options for next steps

**Confidence:**
- HIGH on the X/Z framing divergence (genuine structural distinction)
- HIGH on the R1 gradient-Dirac workaround (standard NCG)
- HIGH on each sub-sprint's individual verdict (MIXED with named candidates)
- MEDIUM on whether the deep-read memo §2 needs full re-write or just a §1 footnote correction

**Files:**
- `debug/sprint_l3e_p3_synthesis_memo.md` (this memo, ~3500 words)
- Cross-references:
  - `debug/subsprint_x_sturmian_latremoliere_verification.md` (X memo)
  - `debug/subsprint_y_w1c_pinstate_diagnostic_memo.md` (Y memo)
  - `debug/subsprint_z_bethe_log_latremoliere_memo.md` (Z memo)
  - `debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md` (deep-read memo with §2 correction needed)
  - `debug/sprint_l3e_p3_rescope_memo.md` (Phase A.2'-A.4' Krein-lift plan)
  - `debug/l3e_p3_phase_a1_literature_audit.md` (Phase A.1 literature audit — scoop of Phase B)
  - `debug/l3e_p3_phase_a2_operator_algebraic_eps_net.md` (Phase A.2 F1+F2+F3 findings)

## §9. Session summary

**Three sub-sprints X+Y+Z landed in parallel. Net verdict: MIXED across all three. Framework gains math.OA-grounded vocabulary for Sturmian-FCI / pin-state-shape / bound-state QED, does NOT autonomously close any open physics wall. Three next-sprint options surfaced (Option A: parallel X.1+Y.1+Z.1, Option B: sequential X.1 then Y.1+Z.1, Option C: return to Krein-lift work).**

Pending PI direction.
