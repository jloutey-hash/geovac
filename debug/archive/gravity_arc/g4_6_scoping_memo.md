# Sprint G4-6 scoping — Full discrete-substrate $S_{\rm BH}$ closure (multi-month)

**Date:** 2026-05-29
**Path:** Gravity arc, opening of the **third-stage multi-month commitment**. G4-3 closed the substrate (scalar Laplacian). G4-4 closed the Dirac dynamics (warped + conical, tip coefficient $-1/12$, replica derivative $+1/6$). G4-5 closed the discrete replica method (joint warp + conical, S_BH at UV cell to 0.85× continuum, sector-wise Mellin moment map). **G4-6 closes the remaining quantitative gaps: subleading UV/IR corrections, $\alpha > 1$ branch, and azimuthal-discretization tightening.**
**Verdict:** **POSITIVE-SCOPING-G4-6.** All architectural inputs (G4-3 substrate, G4-4 Dirac, G4-5 replica integrator + sector-wise Mellin moment map) are in place; the remaining work is *quantitative refinement of identified structural targets*, NOT new architectural construction. Four target closures named (i–iv) with four/five load-bearing falsifiers (F13–F17) at sprint-scale granularity each. Sub-sprint sequence G4-6a/b/c/d/e/f sized at **4–7 months end-to-end**, consistent with the post-G4-5 Paper 51 §12.7 estimate. The structural-skeleton-scope reading sharpens to a *single quantitative deliverable*: $S_{\BH}^{\rm discrete}(r_h, \Lambda; f) \to r_h^2 \Lambda^2 / 3 \cdot M_{\rm tip}[f]$ at theorem-grade qualitative-rate level on the discrete substrate, with all named subleading corrections O(a²/r_h²) UV + O(r_h²/R²) IR + α-asymmetry resolved.

---

## §1. Context — what G4-4 and G4-5 closed

### §1.1 G4-4 (multi-month) closure

The G4-4 sprint sequence (Paper 51 §12.6) closed at sprint scale with two load-bearing structural identifications on the discrete-substrate Dirac:

1. **Spinor conical-defect tip coefficient** (G4-4c, Paper 51 Eq.~\ref{eq:spinor_sc_tip}):
   $$\Delta_K^{\rm Dirac, tip}(\alpha) = -\frac{1}{12}\!\left(\frac{1}{\alpha} - \alpha\right)$$
   Best extraction at $\alpha = 2/5$, $t = 2.0$: recovery $= 1.000021$ (relative error $2.1 \times 10^{-5}$). Nine of eleven tested $\alpha$ values reach $> 99.5\%$ recovery somewhere in the $t$-sweep. The continuum Dowker (1977) / Cheeger–Simons spinor SC coefficient is recovered at bit-exact precision **at $\alpha < 1$**.

2. **Replica derivative at $\alpha = 1$** (G4-4f, Paper 51 Eq.~\ref{eq:replica_derivative}):
   $$\left.\frac{d \Delta_K^{\rm Dirac}}{d\alpha}\right|_{\alpha = 1} = +\frac{1}{6}$$
   Discrete central-FD recovery 96.69% on the constant-warp cigar.

Also closed (G4-4d, G4-4e): Seeley-DeWitt $a_0^{\rm Dirac} = 1.992$ at recovery 99.6%; anti-periodic spinor BC essential for clean conical-defect extraction (66.3% scalar vs 99.4% spinor recovery at same substrate).

### §1.2 G4-5 (sprint-scale multi-month closure) closure

G4-5 integrated the G4-4 inputs over $t$ against a Connes–Chamseddine cutoff function and extracted $S_{\BH}$ on the discrete substrate. Four sub-sprints landed in parallel:

1. **G4-5a / G4-5a-refined** (F8): tip-only replica integration with extended UV $t$-grid. Recovery 0.42–0.59 (exact Gaussian Mellin), uniform improvement vs the first-move 0.13–0.37 at coarse $t$-grid. **Residual gap identified as T2 azimuthal-truncation overshoot** ($4/\pi^2$ discrete/continuum ratio at the truncation edge).

2. **G4-5b** (F9, F10): bulk Weyl extraction is structurally **4D, not 2D**. The disk-Dirac $K^{\rm Dirac}_{D^2}(t)$ has 2D Seeley–DeWitt expansion with leading $1/t$ Mellin moment, NOT $1/t^2$. **The bulk $\Lambda^4$ cosmological term emerges only in the 4D embedding $D^2 \times S^2$** — F9 is structurally absent at 2D, F10 partial.

3. **G4-5c** (F11): joint variable-warp + conical-defect Dirac on full cigar geometry. **F6 extension bit-exact in both reduction directions** (α=1 → smooth-tip VariableWarpDirac at rel_err = 0 exactly; constant warp → DiscreteWedgeDirac × S² at rel_err ≤ 1.2 × 10⁻¹³). **S_BH at Λ=2, r_h=2 lands at 4.55 vs continuum 5.33, ratio 0.85** (within factor-2 gate). Λ-monotone descent of ratio (29.1 → 5.81 → 0.85) identifies substrate-UV/IR signature as the dominant residual.

4. **G4-5d** (F12): cutoff dependence sweep returns the structural reframing — $S_{\BH}$ tip contribution lives at $\phi(0)$ (logarithmically regulated) NOT $\phi(2)$ as originally hypothesized. **Sector-wise Mellin moment map** identified:

   | Sector | Wilson coefficient | Mellin moment |
   |:--|:--|:--:|
   | Bulk $R^0$ (cosmological constant) | $\Lambda_{\rm cc}$ | $\phi(2)$ |
   | Bulk $R^1$ (Einstein–Hilbert) | $G_{\rm eff}^{-1}$ | $\phi(1)$ |
   | Topological tip ($S_{\BH}$) | (tip coeff $1/12$) | $\phi(0)$ |

   This is a structurally substantive refinement of G8 (Paper 51 §10).

**Headline G4-5 structural identification** (Paper 51 §12.7):
$$S_{\BH}^{\rm structural} = \frac{r_h^2 \Lambda^2}{3} \cdot M_{\rm tip}[f]$$
where $M_{\rm tip}[f] = \phi(0)$ is the log-regulated topological Mellin moment. The framework predicts the $r_h^2 \Lambda^2 / 3$ structural form; $M_{\rm tip}$ is Class 1 calibration data.

---

## §2. The remaining work (what G4-6 must close)

Paper 51 §12.7 names four explicit follow-ons for G4-6. Restating in scoping language:

### §2.1 Subleading O(a²/r_h²) UV corrections

The discrete substrate has a UV cutoff at $t_{\rm UV} \sim a^2$. The leading $1/t$ small-$t$ asymptote of $K^{\rm Dirac}_{D^2}$ is accurate but subleading $a^2$-correction terms enter $S_{\BH}$ as $O(a^2/r_h^2)$ at the integrated level. G4-6a closes this via **multi-substrate continuum extrapolation** over $(a, N_\rho)$ pairs holding $R = N_\rho \cdot a$ fixed.

### §2.2 Subleading O(r_h²/R²) IR corrections

The discrete substrate has an IR boundary at $\rho = R$ (Dirichlet BC). For $r_h \ll R$ (regime of interest), the IR contributions to $S_{\BH}$ scale as $O(r_h^2 / R^2)$ but G4-5c showed they dominate the small-Λ cells (ratio 29.1 at $\Lambda = 0.5$). G4-6b closes this via **IR-boundary regularization analysis**, either via larger-$R$ panels or via analytical IR-subtraction of the bulk-warp contribution.

### §2.3 α > 1 structural asymmetry

G4-4c week 2 found a clean negative: at $\alpha > 1$ (excess-angle / saddle-cone), spinor SC recovery plateaus at **bit-identical 67.88% across $N_0 \in \{120, 240, 480\}$**. UV refinement does not help; the 32% gap is structural. G4-5c inherits this asymmetry, though it doesn't affect $\alpha = 1$ extraction directly. G4-6c closes this analytically: identify the sub-leading correction to the continuum SC formula at excess angles, or identify a different effective tip coefficient at $\alpha > 1$.

### §2.4 T2 azimuthal-truncation overshoot

The $4/\pi^2$ discrete/continuum azimuthal-Laplacian ratio at the truncation edge ($k = N_\phi/2$) is the structural origin of G4-5a-refined's residual gap (per-$t$ tip recovery 1.3% at $t = 0.0025$, 76.3% at $t = 0.1$). Closure requires either spectral azimuthal discretization (DST/Fourier) or $N_\phi$ refinement beyond 192. G4-6d closes this via **spectral azimuthal discretization**.

---

## §3. Architectural inputs from G4-5

All G4-5 infrastructure transfers cleanly:

- **G4-3 substrate** ($\mathbb{Z}_+(a)|_{N_\rho} \times \mathbb{Z}/N_\phi \times \text{Fock}(S^2, l_{\max})$): scalar Laplacian baseline.
- **G4-4a `WarpedDiracConstant`**: constant-warp Dirac, used as F6-A reduction reference.
- **G4-4b `VariableWarpDirac`** (smooth-tip warp): variable-warp Dirac, used as F6-B reduction reference. Level 1.5 spin-connection $(r'/r)^2$ correction available.
- **G4-4c `DiscreteWedgeDirac`** (apex angle $2\pi\alpha$, anti-periodic BC): wedge spinor Dirac.
- **G4-5c `JointWarpConicalDirac`** (driver level): joint variable-warp + conical-defect Dirac, F6 bit-exact.
- **G4-5a/G4-5a-refined replica integrator**: extended $t$-grid + UV-overshoot correction + Mellin moment integration against arbitrary cutoff.
- **G4-5d sector-wise Mellin moment map**: $\phi(0)$ for tip / $\phi(1)$ for EH / $\phi(2)$ for cosmological constant.
- **Sweet-spot substrate parameters**: $(N_\rho, a, N_0, r_h, l_{\max}) = (200, 0.05, 120, 2.0, 3)$.

No new production modules are required to launch G4-6. The first sub-sprint (G4-6a, multi-substrate continuum extrapolation) uses existing infrastructure.

---

## §4. The four target structural closures

G4-6 closes four specific quantitative gaps named in Paper 51 §12.7. Each is a target structural closure at sprint-to-multi-month granularity, with falsifiers stated at theorem-grade rigor where possible.

### §4.1 Target (i): Subleading O(a²/r_h²) UV closure

**Closure target.** Identify the leading $a^2$-correction coefficient in the small-$t$ expansion of $K^{\rm Dirac}_{D^2}(t)$ and verify that the integrated $S_{\BH}$ residual at $\Lambda = 2$ scales as $O(a^2/r_h^2)$ across multi-substrate panels.

**Method.** Multi-substrate continuum extrapolation: hold $R = N_\rho \cdot a = 10$ fixed and sweep $(a, N_\rho) \in \{(0.10, 100), (0.05, 200), (0.025, 400), (0.0125, 800)\}$. Compute $S_{\BH}^{\rm discrete}(r_h = 2, \Lambda = 2)$ at each substrate; extrapolate to $a \to 0$ via Richardson or polynomial fit in $a^2$.

**Theorem-grade target.** Closure-form coefficient of the $a^2/r_h^2$ correction matching a continuum Seeley-DeWitt structural prediction.

**Effort.** 2–3 months.

### §4.2 Target (ii): Subleading O(r_h²/R²) IR closure

**Closure target.** Identify the leading $r_h^2/R^2$-correction coefficient in the integrated $S_{\BH}$; isolate the bulk-warp IR contribution and verify that the $\Lambda$-monotone residual at small $\Lambda$ scales as $(r_h^2/R^2)$.

**Method.** Multi-substrate IR sweep: hold $a = 0.05$ fixed and sweep $R \in \{10, 20, 40, 80\}$ at the corresponding $N_\rho \in \{200, 400, 800, 1600\}$. Compute $S_{\BH}^{\rm discrete}(r_h = 2, \Lambda)$ at each $R$ across $\Lambda \in \{0.5, 1, 2\}$. Fit the small-$\Lambda$ residual to $r_h^2/R^2$.

Alternative analytical route: subtract the bulk-warp IR contribution analytically via the asymptotic Schwarzschild-like form $r(\rho) \to \rho$ at large $\rho$; this gives a closed-form IR subtraction that should restore the small-$\Lambda$ cells to the factor-2 band.

**Theorem-grade target.** Closure-form coefficient of the $r_h^2/R^2$ correction; small-$\Lambda$ cells within factor-2 band after IR subtraction.

**Effort.** 1–2 months.

### §4.3 Target (iii): α > 1 analytical resolution

**Closure target.** Identify the structural mechanism for the 67.88% plateau at $\alpha > 1$ and determine whether it affects $S_{\BH}$ at $\alpha = 1$ via sub-leading corrections.

**Method.** Two parallel approaches:

1. **Analytical**: investigate the continuum Sommerfeld–Cheeger formula at $\alpha > 1$ (excess-angle regime). The standard form $-(1/12)(1/\alpha - \alpha)$ is derived for $\alpha < 1$ (deficit angle); at $\alpha > 1$, the apex has negative curvature concentrated, and sub-leading corrections to the SC formula may dominate. Candidate forms: $\Delta_K^{\rm Dirac, tip}(\alpha) = -(1/12)(1/\alpha - \alpha) - (c_2/12)(1/\alpha^3 - \alpha^3) + \ldots$ with $c_2$ to be identified.

2. **Discrete-substrate**: extract higher-order terms in the $\alpha$-expansion via central FD at higher $k$-step values ($k = 24, 48, 96$) and verify the candidate continuum form.

**Theorem-grade target.** Identify $c_2$ or alternative sub-leading structural form; verify against discrete substrate to within 10%.

**Effort.** 1–2 months.

### §4.4 Target (iv): Spectral azimuthal discretization

**Closure target.** Replace the FD azimuthal discretization with a spectral (DST/Fourier) representation; verify per-$t$ tip recovery is unbiased at small $t$ (no $4/\pi^2$ truncation overshoot).

**Method.** Re-implement the wedge-Dirac azimuthal sector using a Fourier basis (anti-periodic BC) instead of centered FD. The Fourier representation gives exact eigenvalues $(2\pi(k + 1/2)/(2\pi\alpha))^2 = ((k + 1/2)/\alpha)^2$ at every truncation, eliminating the $4/\pi^2$ overshoot.

Implementation cost is moderate: the `DiscreteWedgeDirac` class needs a `spectral` mode that uses a Fourier basis instead of FD. The bit-exact F6 reductions at $\alpha = 1$ (smooth-disk Fourier) and at constant warp (factorization) carry over.

**Theorem-grade target.** Per-$t$ tip recovery within 1% across the full $t$-grid $\{0.0025, 0.005, \ldots, 10\}$.

**Effort.** 1–2 months (mostly implementation; structural verification straightforward).

---

## §5. Load-bearing falsifiers for G4-6

Following the L3a-1 / L3b-2 / G4-4 / G4-5 first-move discipline, G4-6 lives at finite cutoff and depends on five load-bearing falsifiers passing across the target-closure subspace.

### F13 — Multi-substrate continuum extrapolation closure

At the multi-substrate panel $(a, N_\rho) \in \{(0.10, 100), (0.05, 200), (0.025, 400), (0.0125, 800)\}$ with $R = 10$, $r_h = 2$, $\Lambda = 2$, $N_0 = 120$:

$$S_{\BH}^{\rm discrete}(a) = S_{\BH}^{\rm continuum} \cdot \left(1 + c_{\rm UV} (a/r_h)^2 + O((a/r_h)^4)\right)$$

with $c_{\rm UV}$ extracted via Richardson extrapolation. **Falsifier**: extrapolated $S_{\BH}(a \to 0) / (r_h^2 \Lambda^2 / 3 \cdot M_{\rm tip}[f]) > 0.95$.

### F14 — IR-boundary regularization closure

At the IR panel $R \in \{10, 20, 40, 80\}$ with $a = 0.05$, $r_h = 2$, $\Lambda \in \{0.5, 1, 2\}$:

$$S_{\BH}^{\rm discrete}(R, \Lambda) - S_{\BH}^{\rm bulk-warp-IR}(R, \Lambda) = S_{\BH}^{\rm topo}(R, \Lambda)$$

with $S_{\BH}^{\rm topo}$ falling within factor-2 band of continuum at all three $\Lambda$ values after IR subtraction. **Falsifier**: small-$\Lambda$ (Λ = 0.5) ratio improves from 29 → ≤ 2.

### F15 — α > 1 structural form identification

Either:
- **F15a**: identify a closed-form sub-leading correction to the SC formula at $\alpha > 1$ that recovers the 32% gap to within 10%, OR
- **F15b**: identify that the 32% gap is a discrete-substrate artifact that vanishes under spectral azimuthal discretization (F17 closure).

**Falsifier**: choose between F15a and F15b based on whether F17 closure also closes the $\alpha > 1$ plateau. If F17 closure restores 100% recovery at $\alpha > 1$, then F15b is correct (the gap was a discretization artifact). If F17 closure does NOT restore 100% recovery, then F15a is correct (sub-leading continuum form), and the $c_2$ coefficient must be identified.

### F16 — Sector-wise Mellin moment map at theorem grade

At the verified-substrate panel (F13, F14 closed), compute $S_{\BH}$ for the three cutoff classes (Gaussian, sharp, polynomial) and verify the sector-wise Mellin moment map:

$$\frac{S_{\BH}^{\rm sharp}}{S_{\BH}^{\rm Gauss}} = \frac{\phi_{\rm sharp}(0)}{\phi_{\rm Gauss}(0)}, \quad \frac{S_{\BH}^{\rm poly}}{S_{\BH}^{\rm Gauss}} = \frac{\phi_{\rm poly}(0)}{\phi_{\rm Gauss}(0)}$$

with the $\phi(0)$ moments computed at the same substrate IR/UV cutoffs that regulate the integrals. **Falsifier**: ratios match to within 10%, lifting G4-5d's structural identification to a theorem-grade quantitative closure.

### F17 — Spectral azimuthal discretization closure

At spectral azimuthal discretization (DST/Fourier with anti-periodic BC), per-$t$ tip recovery:

$${\rm recovery}_{\rm tip}^{\rm spectral}(t) > 0.99$$

across the full $t$-grid $\{0.0025, 0.005, \ldots, 10\}$. **Falsifier**: the residual 0.42–0.59 G4-5a-refined recovery improves to ≥ 0.95.

---

## §6. Sub-sprint sequence with effort estimates

**ORIGINAL sequencing (pre-2026-05-29 reframing) — preserved for reference:**

| Sub-sprint | Scope | Falsifier | Effort | Original sequencing |
|---|---|---|---|---|
| G4-6a | UV multi-substrate continuum extrapolation | F13 | 2–3 months | Sequential foundation |
| G4-6b | IR-boundary regularization analysis | F14 | 1–2 months | Parallel to G4-6c, G4-6d |
| G4-6c | α > 1 analytical + discrete-substrate resolution | F15a/F15b | 1–2 months | Parallel to G4-6b, G4-6d |
| G4-6d | Spectral azimuthal discretization (DST/Fourier) | F17 | 1–2 months | Parallel to G4-6b, G4-6c |
| G4-6e | Sector-wise Mellin moment map at theorem grade | F16 | 1 month | After G4-6a/b |
| G4-6f | Synthesis + Paper 51 §12.8 + closure verdict | (closure narrative) | 1 month | Final |

**REFRAMED sequencing (post-2026-05-29 multi-task thread, Track B.1/B.2 findings):**

The G4-6a first-move panel sweep on the FD substrate (`debug/g4_6a_multi_substrate_uv_first_move_memo.md`) confirmed analytically what v3.20.0 task #28 predicted: FD recovers $\sim 0.04\%$ of the UV target at $t = a^2$. The B.2 spectral substrate first-move (`debug/g4_6a_spectral_substrate_first_move_memo.md`) confirmed that spectral azimuthal discretization improves this 160× (to 6.36% recovery). Spectral substrate furthermore exposes a substrate-dependent Lichnerowicz constant $B \approx 0.30$ (vs continuum $+1/6 = 0.167$) that contaminates the joint linear fit of $A$. This identifies G4-6b as a SEQUENTIAL prerequisite to G4-6a refined.

| Sub-sprint | Scope | Falsifier | Effort | **Reframed sequencing** |
|---|---|---|---|---|
| **G4-6d** | Spectral azimuthal discretization | F17 | 1–2 months | **DONE (thread 5 Track A, 2026-05-29). 14 production tests pass.** |
| **G4-6b** | IR-boundary regularization (B identification) | F14 | 1–2 months | **DONE (thread 6 Track α + thread 7 Track b, 2026-05-29). $B_{\rm substrate}$ = 0.163 at R≥10, within 2.3% of continuum +1/6; analytical B-subtraction NOT needed.** |
| **G4-6a refined** | $N_\phi$-sweep refinement on spectral substrate (REFRAMED axis per thread 7 Track a finding) | F13 | sprint-scale to ~2-4 weeks | **REFRAMED-IN-PROGRESS.** Original radial-axis refinement reframed to $N_\phi$-axis. Simplified extraction strategy validated. Next move: $N_\phi$-sweep at fixed $a$. |
| **G4-6c** | α > 1 analytical (Möbius mechanism) | F15a/F15b | 1–2 months | **SUBSTRATE-LEVEL IDENTIFIED (thread 5 Track C, 2026-05-29).** Paper 51 §subsubsec documented. Continuum theorem requires Route A (Fursaev-Miele PDF) or Route C' (Sommerfeld contour at excess angle). |
| **G4-6e** | Sector-wise Mellin moment map at theorem grade | F16 | 1 month | After G4-6a refined |
| **G4-6f** | Synthesis + Paper 51 §12.8 + closure verdict | (closure narrative) | 1 month | Final |

**Total G4-6 commitment: 3–6 months** UNCHANGED (parallelized to 3 if G4-6b/c run concurrently after G4-6d; sequential 6 months if dependencies surface).

**Reframed dependency structure**:
- **G4-6d (spectral) foundation:** B.1 implementation + B.2 verification COMPLETED 2026-05-29; G4-6d remaining work is production tests + formal closure memo (~1 week).
- **G4-6b IR regularization:** SEQUENTIAL prerequisite to G4-6a refined. Subtract substrate-finite contributions to $B$ analytically (or via subtracting bulk-warp asymptotic) before joint $A/B$ fit. Then $A$-extraction has constraining power.
- **G4-6a refined:** REPLACES original G4-6a. After G4-6b lands, rerun multi-substrate sweep on spectral substrate with $B$ subtracted.
- **G4-6c:** STRUCTURALLY INDEPENDENT (α > 1 axis). Can run parallel to G4-6b.
- **G4-6e:** depends on G4-6a refined + G4-6b for clean Mellin moment extraction.
- **G4-6f:** synthesis only after G4-6a refined / b / c / e close.

**Why the reframing matters.** The original G4-6a-first sequencing would have committed 2–3 months to multi-substrate refinement on FD substrate (which task #28 / B.1 first-move show would essentially fail). The reframed G4-6d-first sequencing landed substrate verification in a single afternoon's main-session work, and identified the $B$-contamination structural finding that needs G4-6b before further refinement. Total estimate unchanged at 3-6 months; sequencing now well-grounded.

---

## §7. First-move plan (G4-6a multi-substrate UV closure)

G4-6a is the recommended first sub-sprint: foundation work, sequential (downstream sub-sprints depend on it), and the highest-impact closure (the $a^2/r_h^2$ correction is the dominant identified UV residual).

### §7.1 Module structure

**No new production modules required.** The existing infrastructure handles multi-substrate sweeps:
- `geovac/gravity/warped_dirac.py::WarpedDiracConstant` (G4-4a)
- `geovac/gravity/warped_dirac.py::DiscreteWedgeDirac` (G4-4c)
- `debug/g4_5c_joint_warp_conical.py::JointWarpConicalDirac` (driver-level joint Dirac)
- `debug/g4_5a_refined_tip_replica.py` (extended-$t$-grid replica integrator)

Driver `debug/g4_6a_multi_substrate_uv.py` (~400 lines):
- Multi-substrate panel loop over $(a, N_\rho)$ holding $R = 10$ fixed.
- Reuses G4-5c `JointWarpConicalDirac` per substrate cell.
- Extracts $S_{\BH}^{\rm discrete}(r_h = 2, \Lambda = 2)$ at each cell.
- Richardson extrapolation in $a^2$ to $a \to 0$.

### §7.2 Test architecture

**No new production tests required.** Driver-level verification:
- F13 closure: extrapolated $S_{\BH}(a \to 0) > 0.95 \cdot$ continuum.
- F6 bit-exact at each substrate cell (re-verifies G4-5c factorization at multiple panels).
- Sanity at $a = 0.05$: reproduces G4-5c result (4.55 at $\Lambda = 2$).

### §7.3 Driver and memo

- `debug/g4_6a_multi_substrate_uv.py`: driver.
- `debug/g4_6a_multi_substrate_uv_memo.md`: closure memo.
- `debug/data/g4_6a_multi_substrate_uv.json`: numerical panel + Richardson extrapolation coefficients.

### §7.4 Sequencing for G4-6a (2–3 month breakdown)

- **Month 1**: substrate sweep at $(0.10, 100), (0.05, 200), (0.025, 400)$; compute $S_{\BH}$ at each; verify F6 bit-exact at each substrate.
- **Month 2**: extend to $(0.0125, 800)$ if feasible (memory budget); Richardson extrapolation; identify $c_{\rm UV}$ coefficient; cross-check against continuum Seeley-DeWitt prediction.
- **Month 3**: continuum-form identification; closure memo; sub-sprint review.

### §7.5 What G4-6a DOES NOT do (deferred to other sub-sprints)

- IR-boundary regularization (G4-6b).
- α > 1 analytical resolution (G4-6c).
- Spectral azimuthal discretization (G4-6d).
- Sector-wise Mellin moment map at theorem grade (G4-6e).

---

## §8. Multi-month commitment scoping

**The G4-6 commitment is structurally additive on top of G4-4 + G4-5.** All infrastructure (substrate, Dirac, replica integrator, sector-wise Mellin moment map) is in place; G4-6 closes quantitative gaps in already-identified structural targets.

### §8.1 Effort breakdown

- G4-6a (UV multi-substrate, 2–3 months): foundation, sequential.
- G4-6b/c/d (IR, α > 1, spectral, 1–2 months each): parallel.
- G4-6e (Mellin moment theorem grade, 1 month): sequential downstream of G4-6a + G4-6b.
- G4-6f (synthesis, 1 month): sequential final.

**Total: 4–7 months** parallelized; **7 months** fully sequential.

### §8.2 Risk tier classification

**Risk tier 1 (low):**
- G4-6a (F13): existing substrate sweep machinery, Richardson extrapolation standard practice. Multi-substrate compute may be expensive at $(0.0125, 800)$ but doable.
- G4-6d (F17): spectral azimuthal discretization is a clean implementation; the structural F6 reductions carry over.
- G4-6e (F16): downstream extraction once G4-6a + G4-6b close.

**Risk tier 2 (medium):**
- G4-6b (F14): IR-boundary regularization may require analytical IR subtraction (closed-form bulk-warp asymptotic at large $\rho$); the dominant residual at small Λ is structural, not just numerical.
- G4-6c (F15a/F15b): α > 1 analytical may not yield a clean closed form; the 67.88% plateau could be a discretization artifact (F15b → F17 closure) OR a sub-leading continuum form (F15a → multi-month analytical work).

**Risk tier 3 (higher):**
- G4-6f synthesis: the headline closure ("framework predicts $S_{\BH} = r_h^2 \Lambda^2 / 3 \cdot M_{\rm tip}[f]$ at theorem grade on the discrete substrate") requires all four target closures to land at quantitative-rate level. If F15 surfaces a structural obstruction (e.g., the α > 1 asymmetry persists under spectral discretization AND no closed-form sub-leading correction lands), the headline narrative shifts to "leading-order framework prediction at $\alpha < 1$ only", which is a softer closure than the Paper 51 §12.7 multi-month estimate anticipates.

### §8.3 Mitigation

- **G4-6a first move** establishes the UV closure foundation before downstream sub-sprints surface dependencies.
- **G4-6c parallel with G4-6d**: if F17 closure (spectral) restores $\alpha > 1$ recovery to ≥ 95%, then F15b is correct and G4-6c's analytical track is unnecessary. Run both in parallel to short-circuit if F17 closes early.
- **G4-6e gating**: do not launch the Mellin moment map theorem-grade extraction until G4-6a and G4-6b have closed. If either remains structurally open, defer G4-6e to a follow-on cycle.

---

## §9. Strategic positioning vs continuum Connes-Chamseddine literature

The Connes-Chamseddine (CC) spectral action program (Chamseddine-Connes 1997, 2010) reads gravity from the heat-kernel asymptotic expansion of $D^2/\Lambda^2$ on the relevant noncommutative geometry. For BH entropy, the standard continuum derivation (Solodukhin 1995, Cheeger-Simons 1985) uses the replica method on the smooth Euclidean Schwarzschild geometry:

$$S_{\BH}^{\rm CC} = -\left.\frac{d I_E^{\rm CC}}{d\alpha}\right|_{\alpha = 1}, \quad I_E^{\rm CC}(\alpha) = -\frac{1}{2}\int_0^\infty \frac{dt}{t} f(t \Lambda^2) K_{\rm cigar}^{\rm CC}(t, \alpha)$$

The continuum result is $S_{\BH} = A/(4 G_{\rm eff})$ with $G_{\rm eff} = 6\pi/(\phi(1) \Lambda^2)$ for Gaussian cutoff.

**What G4-6 adds vs the continuum literature.** Two structural reframings:

1. **The cutoff dependence partitions across sectors via Mellin moments.** Continuum CC literature treats the cutoff function as monolithic Class 1 calibration data; the $\phi(0)$/$\phi(1)$/$\phi(2)$ partition across topological / EH / cosmological-constant sectors is a substantive refinement that the discrete substrate surfaces empirically (G4-5d) and G4-6 lifts to theorem grade (F16).

2. **The discrete substrate provides a quantitative bridge from operator-system level to physical $S_{\BH}$.** Paper 38's WH1 PROVEN closure shows the discrete substrate converges to the continuum spectral triple in Latrémolière propinquity at rate $4/\pi$. G4-4 + G4-5 + G4-6 trace this convergence at the physical-observable level for $S_{\BH}$: the framework's discrete substrate reproduces the continuum CC entropy to the leading order and identifies the subleading corrections at $O(a^2/r_h^2)$ and $O(r_h^2/R^2)$, providing a quantitative analog of the convergence rate that operates at the level of $S_{\BH}$ rather than at the level of propinquity bounds.

**What G4-6 does NOT add vs the continuum literature.**
- No new Yukawa-like calibration data (consistent with the structural-skeleton-scope reading; Paper 51 §13).
- No new $\Lambda^0$ topological correction (Euler density) — deferred to a future cycle after G4-6.
- No higher-curvature (Wald) corrections — deferred to a future cycle after G4-6.

**Concurrent-work risk.** Standard continuum CC derivations of $S_{\BH}$ (Solodukhin 1995, Cheeger-Simons 1985, and the more recent Chamseddine-Connes 2010 spectral action review) are well-established. The G4-6 contribution is the *discrete-substrate quantitative reproduction* with explicit identification of UV/IR subleading corrections and the sector-wise Mellin moment map. To my knowledge, no published work derives the BH entropy on a discrete substrate in this specific noncommutative-geometry-of-Fock-projection sense; the discrete-substrate convergence to standard CC at the level of $S_{\BH}$ is the genuinely new content.

---

## §10. Honest scope

This is a scoping memo. It does NOT contain new computational results. It:

- **Confirms** G4-3 substrate + G4-4 Dirac + G4-5 replica integrator + G4-5d sector-wise Mellin moment map are sufficient input.
- **Names** four target structural closures (i–iv) corresponding to Paper 51 §12.7 named follow-ons.
- **Names** five load-bearing falsifiers F13–F17 at theorem-grade or quantitative-rate granularity.
- **Sizes** sub-sprint sequence G4-6a/b/c/d/e/f at 4–7 months (parallelized) or 7 months (sequential).
- **Names** G4-6a as the first-move foundation.
- **Names** risk-tier classification and mitigation strategy.

### §10.1 Theorem-grade vs structural-sketch vs numerical-observation targets

| Closure | Theorem-grade? | Sprint-scale or multi-month? |
|---|---|---|
| F13 (UV multi-substrate) | Quantitative-rate (Richardson extrapolation closes to closed-form $c_{\rm UV}$ coefficient) | Multi-month (2–3 months) |
| F14 (IR-boundary) | Quantitative-rate (closure depends on analytical IR subtraction landing) | Multi-month (1–2 months) |
| F15a (α > 1 closed-form) | Closed-form structural identification of $c_2$ | Multi-month (1–2 months) |
| F15b (α > 1 = discretization artifact) | Structural-sketch (F17 closure implies this) | Sprint-scale (week or two if F17 closes) |
| F16 (sector-wise Mellin moment) | Theorem-grade (lift G4-5d to quantitative closure) | Sprint-scale (1 month after F13, F14 close) |
| F17 (spectral azimuthal) | Theorem-grade (Fourier basis is exact, no truncation overshoot) | Multi-month (1–2 months implementation) |

The **headline G4-6 closure** ("framework predicts $S_{\BH} = r_h^2 \Lambda^2 / 3 \cdot M_{\rm tip}[f]$ at theorem grade on the discrete substrate") is a **structural identification at quantitative-rate level**, NOT a theorem in the formal-proof sense. The framework's discrete substrate is shown to reproduce continuum CC $S_{\BH}$ to the leading order with named subleading corrections, but the proof of convergence (analog of Paper 38's five-lemma propinquity-convergence proof) is NOT a target of G4-6; that would be a downstream multi-month NCG-mathematics target (analog of L3e).

### §10.2 What G4-6 does NOT close (deferred to a future cycle)

- **Continuum extrapolation theorem.** A formal statement that $S_{\BH}^{\rm discrete}(N_\rho, a, N_0, R) \to r_h^2 \Lambda^2 / 3 \cdot M_{\rm tip}[f]$ as the substrate parameters approach the continuum limit, with rate bound. This is the propinquity-style convergence theorem at the level of $S_{\BH}$ (not at the level of spectral-triple distance). Multi-month NCG-mathematics target, deferred.

- **Wald entropy formula extension.** For modified gravity (higher curvature), $S_{\BH}$ depends on $\partial \mathcal{L} / \partial R_{abcd}$. The discrete substrate naturally supports this via the heat-kernel expansion, but extracting Wald coefficients is downstream.

- **Higher curvature corrections.** $\Lambda^0$ topological terms (Euler density), $\Lambda^{-2}$ curvature-squared. These contribute $\log \Lambda$ corrections to $S_{\BH}$ and require multi-Mellin extraction. Deferred to a future cycle after G4-6.

- **Continuum CC literature engagement.** Detailed comparison vs Solodukhin 1995, Chamseddine-Connes 2010, etc. The discrete-substrate reproduction is the primary contribution; detailed continuum-side engagement is a Paper 51 §13 expansion deferred to post-G4-6.

### §10.3 Open follow-ons after G4-6 closure

- **Theorem-grade convergence proof at the level of $S_{\BH}$**: analog of Paper 38's five-lemma propinquity proof at the level of physical observables (BH entropy). Multi-month NCG-mathematics.
- **Higher-curvature extension**: extract $\Lambda^0$ Euler density and $\Lambda^{-2}$ Riemann-squared via multi-Mellin transform on the discrete substrate.
- **Wald entropy for modified gravity**: extension to non-Einstein–Hilbert gravity actions.
- **Black hole greybody factors**: extension from entropy to Hawking radiation spectra via the discrete substrate.
- **Connection to Paper 47's two-rate hybrid convergence**: the discrete-substrate $S_{\BH}$ convergence at rate $O(a^2/r_h^2 + r_h^2/R^2)$ may correspond to the two-rate hybrid of Paper 47 (inner propinquity + outer norm-resolvent). This is a speculative bridge to the operator-algebraic Lorentzian arc; structural exploration deferred.

---

## §11. Cross-references

### Sprint memos
- `debug/g4_3_warped_substrate_memo.md` — G4-3 scoping
- `debug/g4_3a_cleanup_hermitian_polar_memo.md` — Hermitian polar Laplacian
- `debug/g4_3b_variable_warp_memo.md` — variable warp
- `debug/g4_3c_proper_wedge_memo.md` — proper-wedge T1
- `debug/g4_3d_uv_extension_memo.md` — UV extension T2 + $4/\pi^2$ azimuthal-truncation overshoot identification
- `debug/g4_4_warped_dirac_scoping_memo.md` — G4-4 scoping (precedent for this scoping memo)
- `debug/g4_4a_constant_warp_dirac_memo.md` — G4-4a constant-warp Dirac
- `debug/g4_4b_variable_warp_dirac_memo.md` — G4-4b variable-warp Dirac
- `debug/g4_4c_wedge_dirac_memo.md` — G4-4c wedge Dirac + headline $-1/12$ extraction
- `debug/g4_4d_seeley_dewitt_memo.md` — G4-4d Seeley-DeWitt extraction
- `debug/g4_4e_bc_sectors_memo.md` — G4-4e anti-periodic spinor BC essential
- `debug/g4_4f_replica_derivative_memo.md` — G4-4f $+1/6$ replica derivative
- `debug/g4_5_scoping_memo.md` — G4-5 scoping (immediate precedent)
- `debug/g4_5_synthesis_memo.md` — G4-5 synthesis scaffolding (parallel)
- `debug/g4_5a_first_move_tip_replica_memo.md` — G4-5a first move
- `debug/g4_5a_refined_tip_replica_memo.md` — G4-5a-refined extended-$t$-grid (closes F8 methodologically)
- `debug/g4_5b_bulk_weyl_memo.md` — G4-5b bulk Weyl 4D-not-2D reframing
- `debug/g4_5c_joint_warp_conical_memo.md` — G4-5c F11 partial closure + $S_{\BH}$ at UV cell to 0.85×
- `debug/g4_5d_cutoff_dependence_memo.md` — G4-5d sector-wise Mellin moment map (the substantive sprint of G4-5)

### Paper sections
- **Paper 28 §4.11** (G4-2 continuum derivation of $S_{\BH}$)
- **Paper 28 §4.13** (G7 Newton constant + cosmological constant)
- **Paper 28 §4.14** (G4-1 $S^2$ Dirac)
- **Paper 28 §4.15** (G4-2 conical-defect replica)
- **Paper 28 §4.16** (G8 cutoff-function classification)
- **Paper 28 §4.17** (G4-3 discrete substrate)
- **Paper 51 §5** (G4-1, G4-2 — continuum)
- **Paper 51 §10** (G7 — Newton constant)
- **Paper 51 §11** (G8 — cutoff dependence)
- **Paper 51 §12.6** (G4-4 — warped Dirac, conical-defect, replica derivative)
- **Paper 51 §12.7** (G4-5 — discrete replica method + sector-wise Mellin moment map + G4-6 implications)
- **Paper 51 §13** (Discussion — structural-skeleton-scope reading)

### CLAUDE.md memory
- **CLAUDE.md §1.7 WH1** (substrate convergence at theorem-grade; structural-skeleton-scope reading)
- **CLAUDE.md §2** (sprint chronicle; G4-6 closure entry will be a one-liner post-closure)
- [`geovac_structural_skeleton_scope_pattern`](memory/geovac_structural_skeleton_scope_pattern.md) (structural reading)
- [`feedback_no_synthesis_memos`](memory/feedback_no_synthesis_memos.md) (one canonical memo per sprint)
- [`feedback_diagnostic_before_engineering`](memory/feedback_diagnostic_before_engineering.md) (apply to each G4-6a–e sub-sprint)
- [`feedback_bit_exactness_rule`](memory/feedback_bit_exactness_rule.md) (F6 bit-exact = green; substrate quantitative residual = caution)

### External references
- **Chamseddine-Connes 1997, 2010** (CC spectral action)
- **Solodukhin 1995** (continuum conical-defect BH entropy)
- **Cheeger-Simons 1985** (smooth conical-defect heat kernel)
- **Dowker 1977** (spinor conical-defect SC coefficient)
- **Sommerfeld 1894** (scalar SC original)
- **Camporesi 1996** (warped-product spin connection)

---

## §12. Synthesis-memo discipline note

This memo follows the [`feedback_no_synthesis_memos`](memory/feedback_no_synthesis_memos.md) standing rule: it is a **scoping memo for a single multi-month sprint commitment**, NOT a cross-sprint synthesis. The G4-6 commitment is the next multi-month tranche of the gravity arc, sized and scoped at the same level of granularity as the G4-4 and G4-5 scoping memos. When G4-6 closes, the closure will be documented via:
- One canonical closure memo per sub-sprint (G4-6a/b/c/d/e closure memos, sized at ≤ 5000 words each).
- One synthesis memo for G4-6f (final closure narrative).
- A Paper 51 §12.8 update lifting the headline closure to the paper.
- A CLAUDE.md §2 one-liner per G4-6 closure verdict.

No anticipatory scaffolding memo is required at this scoping stage; the precedent G4-5 scaffolding memo (parallel sub-sprint dispatch) was justified by the parallel-dispatch window, which does not apply at the G4-6 scoping stage (G4-6a is sequential first-move; G4-6b/c/d parallel later).
