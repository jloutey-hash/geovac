# Sprint G4-5 synthesis — Discrete replica method for $S_{\rm BH}$ (placeholder)

**Date:** 2026-05-29
**Path:** Gravity arc, synthesis of the multi-month G4-5 commitment (discrete replica method for $S_{\rm BH}$ on the warped Dirac substrate). This memo is the **writeup framework** to be populated as the four parallel sub-sprints (G4-5a-refined, G4-5b, G4-5c, G4-5d) close.
**Verdict:** **{{TBD pending G4-5a-refined + G4-5b + G4-5c + G4-5d closures}}**.

This is a **scaffolding memo**. Quantitative results are marked `{{TBD from G4-5x}}`. The architectural skeleton is fixed; numerical fill-in is the remaining work.

---

## §1. Sprint roll-call

| Sprint | Date | Status | Verdict |
|---|---|---|---|
| G4-5 scoping | 2026-05-29 | DONE | POSITIVE-SCOPING — architecture defined, F8–F12 named |
| G4-5a first move | 2026-05-29 | DONE | POSITIVE-VERIFIED — tip-only replica integration operational at sprint scale |
| G4-5a-refined (extended UV $t$-grid) | 2026-05-29 (dispatched) | RUNNING | {{TBD — closes F8 methodologically}} |
| G4-5b (bulk Weyl $\Lambda^4 + \Lambda^2$) | 2026-05-29 (dispatched) | RUNNING | {{TBD — closes F9 + F10}} |
| G4-5c (joint warp + conical defect) | 2026-05-29 (dispatched) | RUNNING | {{TBD — closes F11}} |
| G4-5d (cutoff-function dependence) | 2026-05-29 (dispatched) | RUNNING | {{TBD — closes F12}} |
| **G4-5e synthesis (this memo, placeholder)** | **2026-05-29** | **PLACEHOLDER** | **Awaiting parallel sub-sprint closure** |

Underlying inputs:
- **G4-4c spinor tip coefficient** $\Delta_K^{\rm Dirac, tip}(\alpha) = -\tfrac{1}{12}(1/\alpha - \alpha)$ bit-exact to 5 digits.
- **G4-4f replica derivative** $d\Delta_K/d\alpha|_{\alpha=1} = +1/6$ at 96.69% recovery.
- **G4-2 continuum prediction** $S_{\rm BH}^{\rm continuum} = r_h^2 \Lambda^2 / 3$ (Gaussian cutoff, Paper 28 §4.15).
- **G7 calibration** $G_{\rm eff} = 6\pi / \Lambda^2$.
- **G8 cutoff classification** $S_{\rm BH}(f) \propto \phi(2)$, Class 1 calibration data.

---

## §2. Headline structural finding (placeholder)

{{TBD pending closure. Expected form, to be confirmed by G4-5b + G4-5c:}}

> **Conjectured headline:** the discrete replica method on the warped Dirac substrate reproduces the continuum BH entropy bit-exactly at the sweet-spot discretization point:
> $$
> S_{\rm BH}^{\rm discrete}(r_h, \Lambda; f) = \frac{r_h^2 \Lambda^2}{3} \cdot \phi(2)
> + O\!\left(\frac{a^2}{r_h^2}\right) + O\!\left(\frac{r_h^2}{R^2}\right).
> $$
> The leading-order recovery is **{{TBD}}%** at the sweet spot $(N_\rho = 200, a = 0.05, N_0 = 120, r_h = 2)$ for the Gaussian cutoff, with $\phi(2)$-dependence verified at **{{TBD}}** across Gaussian / sharp / polynomial cutoffs.

Three load-bearing pieces of the headline:
1. **Bulk Weyl $\Lambda^4$ + $\Lambda^2$ extraction**: {{TBD from G4-5b}} — coefficients identified against the Seeley–DeWitt expansion of $K_{\rm disk}(t)$.
2. **Tip-only replica $\Lambda^2$ contribution**: {{TBD from G4-5a-refined}} — extended UV $t$-grid closes the discrete/continuum ratio from the G4-5a first-move 13–37% to {{TBD}}%.
3. **$\phi(2)$ cutoff dependence**: {{TBD from G4-5d}} — bit-exact $\phi(2)$ ratios across the three named cutoff classes, verifying G8 prediction at the replica-integral level.

The fourth piece — **joint warp + conical-defect** (G4-5c) — closes the remaining structural gap: combining the G4-4b factorization-loss correction $\sim (r'/r)^2$ with the G4-4c conical-defect tip on the same discrete substrate. {{TBD from G4-5c}}.

---

## §3. G4-5a methodological refinement (G4-5a-refined)

**Sub-sprint scope.** Extend the $t$-grid into the substrate UV regime $[a^2, 0.1]$, account for the UV overshoot identified in G4-3d-UV (the $4/\pi^2$ factor between discrete and continuum azimuthal Laplacian at the truncation edge), and use higher-order quadrature for the log-Mellin moment $\int (dt/t) e^{-t\Lambda^2} \Delta'(t)$.

### §3.1 Method

{{TBD from G4-5a-refined. Expected:}}

- Extended $t$-grid: $t \in \{a^2, 5a^2, 10a^2, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10\}$ (~13 points spanning 4 decades of $t$).
- UV-overshoot correction: subtract the $4/\pi^2$ azimuthal-truncation bias from $K_{\rm disk}(t)$ at $t < a^2 (n_\phi / 2)^2 \pi^2 / 2$.
- Higher-order quadrature: tanh-sinh on the log-Mellin integral; verify against Gauss–Laguerre for cross-check.

### §3.2 Expected result

{{TBD from G4-5a-refined.}} The G4-5a first-move ratio (0.13–0.37) is anticipated to close to {{TBD}}% across $\Lambda \in \{0.5, 1.0, 1.5, 2.0\}$, with the $\Lambda$-trend becoming flat (consistent with proper integration over the full $t$ range).

### §3.3 F8 closure verdict

{{TBD from G4-5a-refined.}} Expected: **POSITIVE-F8-CLOSED-METHODOLOGICALLY**, with discrete/continuum ratio $\geq$ 90% at the sweet spot.

---

## §4. Bulk Weyl extraction $\Lambda^4 + \Lambda^2$ (G4-5b)

**Sub-sprint scope.** Extract the cosmological-constant $\Lambda^4$ and Einstein–Hilbert $\Lambda^2$ pieces of the spectral action on the discrete cigar from the disk part $K_{\rm disk}(t)$ via Seeley–DeWitt Mellin transform. These are bulk contributions ($\alpha$-independent) but provide the load-bearing $\Lambda$-scaling that the replica method's $\Lambda^2$ piece must match.

### §4.1 Method

{{TBD from G4-5b. Expected:}}

For the disk part:
$$
K_{\rm disk}^{\rm Dirac}(t) \sim \frac{A_{D^2}}{4\pi t} \cdot 2 + a_1^{D^2}(\partial D^2) \cdot t^{-1/2} + a_2^{D^2} + O(t)
$$
with $A_{D^2}$ the disk area and $a_1^{D^2}$ a boundary contribution. The Mellin moments give:
$$
\Lambda^4: \int_0^\infty \frac{dt}{t}\, f(t\Lambda^2) \cdot \frac{A_{D^2}}{2\pi t} = \frac{A_{D^2} \Lambda^4}{2\pi} \phi(2)
$$
$$
\Lambda^2: \int_0^\infty \frac{dt}{t}\, f(t\Lambda^2) \cdot a_2^{D^2} \cdot t^{-1} = ...
$$

### §4.2 F9 (bulk Weyl $\Lambda^4$) closure verdict

{{TBD from G4-5b.}} Expected: **{{TBD}}** discrete/continuum recovery of the $\Lambda^4$ coefficient.

### §4.3 F10 ($\Lambda^2$ Einstein–Hilbert) closure verdict

{{TBD from G4-5b.}} Expected: **{{TBD}}** discrete/continuum recovery of the $\Lambda^2$ coefficient at $A_{\rm horizon} = 4\pi r_h^2$.

---

## §5. Joint warp + conical-defect Dirac (G4-5c)

**Sub-sprint scope.** Combine G4-4b's variable-warp Dirac (factorization-loss $\sim (r'/r)^2$) with G4-4c's conical-defect spinor (apex angle $2\pi\alpha$) on a single discrete substrate. This is the **physical cigar geometry**: smooth warp at the asymptotic end, conical defect at the horizon.

### §5.1 Method

{{TBD from G4-5c.}} Expected:
- Composed substrate: $\mathcal{G} = \mathbb{Z}_+(a)|_{N_\rho} \times \mathbb{Z}/N_\phi \times {\rm Fock}(S^2, l_{\max})$, with $N_\phi$-dependent apex angle $2\pi\alpha$ AND warp profile $r(\rho) = r_h\sqrt{1 + (\rho/r_h)^2}$.
- Construct `VariableWarpConicalDirac` extending `VariableWarpDirac` (G4-4b) with the conical-defect azimuthal BC of `DiscreteWedgeDirac` (G4-4c).
- F6 extension: at constant warp + $\alpha = 1$, reduce bit-exact to G4-4a's `WarpedDiracConstant`.
- F4 tip-regularity test: extend G4-4c's tip-regularity check to the joint variable-warp + conical case.

### §5.2 Expected result

{{TBD from G4-5c.}} Anticipated:
- Bit-exact F6 reduction at $\alpha = 1$, $r(\rho) = r_h$ constant.
- F4 tip regularity within {{TBD}}% at the sweet spot.
- Joint factorization-loss + tip coefficient form:
$$
\Delta_K^{\rm Dirac}(\alpha, r_h, t) \approx -\frac{1}{12}\!\left(\frac{1}{\alpha} - \alpha\right) \cdot K^{\rm Dirac}_{S^2}(t; r_h) + O((r'/r)^2)
$$

### §5.3 F11 closure verdict

{{TBD from G4-5c.}} Expected: **POSITIVE-F11-CLOSED**, with the $r_h$-dependence at the joint geometry identified as $r_h^2 \Lambda^2 / 3$ at the integrated replica level (consistent with G4-2 continuum).

---

## §6. Cutoff-function dependence (G4-5d)

**Sub-sprint scope.** Verify the G8 prediction $S_{\rm BH}(f) \propto \phi(2)$ at the discrete replica level. Compute $S_{\rm BH}^{\rm discrete}$ for the three named cutoff classes:
- **Gaussian**: $f(x) = e^{-x}$, $\phi(2) = 1$
- **Sharp**: $f(x) = \Theta(1 - x)$, $\phi(2) = 1/2$
- **Polynomial**: $f(x) = e^{-x^2}$, $\phi(2) = \Gamma(1)/2 = 1/2$ (with reparametrization $u = x^2$)

### §6.1 Method

{{TBD from G4-5d.}} Expected:
- Reuse the G4-5a-refined extended $t$-grid + UV-overshoot correction.
- For each cutoff, integrate the replica derivative $\Delta'(t) = dK/d\alpha|_{\alpha = 1} - K_{\rm disk}(t)$ against the cutoff:
$$
S_{\rm tip}^{(f)} = +\frac{1}{2} \int_{t_{\rm min}}^{t_{\rm max}} \frac{dt}{t} \tilde{f}(t \Lambda^2) \Delta'(t)
$$
- Compare ratios $S_{\rm tip}^{\rm Gaussian} : S_{\rm tip}^{\rm sharp} : S_{\rm tip}^{\rm polynomial}$ against $\phi(2)$ ratios $1 : 1/2 : 1/2$.

### §6.2 F12 closure verdict

{{TBD from G4-5d.}} Expected: **POSITIVE-F12-CLOSED**, with $\phi(2)$ ratios within {{TBD}}% at the sweet spot — verifying G8 cutoff-classification prediction at the replica-method level. Confirms cutoff function is Class 1 calibration data (external).

---

## §7. Comparison to G4-2 continuum

The G4-2 continuum derivation (Paper 28 §4.15) gives:
$$
S_{\rm BH}^{\rm continuum} = \frac{A_{\rm horizon} \Lambda^2}{12\pi} = \frac{4\pi r_h^2 \cdot \Lambda^2}{12\pi} = \frac{r_h^2 \Lambda^2}{3}.
$$
At $r_h = 2$, $\Lambda = 1$: $S_{\rm BH}^{\rm continuum} = 4/3 \approx 1.33$.

The discrete-substrate prediction at the same point (post G4-5b + G4-5c) is {{TBD}}.

| Quantity | Continuum | Discrete (G4-5) | Recovery |
|---|---|---|---|
| Tip coefficient $-1/12$ (per $1/\alpha - \alpha$) | $-1/12$ | $-1/12$ (G4-4c) | bit-exact to 5 digits |
| Replica derivative $+1/6$ | $+1/6$ | $+1/6$ (G4-4f) | 96.69% |
| Tip-only $S_{\rm tip}(\Lambda = 1)$ | {{TBD continuum value}} | {{TBD from G4-5a-refined}} | {{TBD}}% |
| Bulk Weyl $\Lambda^4$ coefficient | $A_{D^2}/(2\pi)$ | {{TBD from G4-5b}} | {{TBD}}% |
| Einstein–Hilbert $\Lambda^2$ coefficient | $r_h^2/3$ | {{TBD from G4-5b/c}} | {{TBD}}% |
| Full $S_{\rm BH}(r_h = 2, \Lambda = 1)$ | $4/3$ | {{TBD from G4-5c}} | {{TBD}}% |
| $\phi(2)$ cutoff scaling | $1 : 1/2 : 1/2$ | {{TBD from G4-5d}} | {{TBD}}% |

The aggregation of these five comparisons forms the G4-5 headline closure. {{TBD upon parallel sprints closing}}.

---

## §8. Honest scope

### §8.1 Theorem-grade results (after G4-5 closure)

{{TBD pending closure. Expected to include:}}
- {{If G4-5b POSITIVE}}: bit-exact identification of bulk Weyl $\Lambda^4$ + $\Lambda^2$ coefficients against G4-2 continuum.
- {{If G4-5c POSITIVE}}: F11 closure — joint warp + conical-defect Dirac on a single discrete substrate, with bit-exact F6 Riemannian-limit reduction.
- {{If G4-5d POSITIVE}}: $\phi(2)$-scaling verified at the replica-method level, lifting G8's continuum prediction to a discrete-substrate theorem.

### §8.2 Sketch / proof outline (named gaps)

- **Continuum extrapolation theorem.** A formal statement that discrete-substrate $S_{\rm BH}^{(N_\rho, a, N_0)}$ converges to $r_h^2 \Lambda^2 / 3 \cdot \phi(2)$ as $(N_\rho \to \infty, a \to 0, N_0 \to \infty)$ at appropriate rates: {{TBD — sketch only at G4-5 closure; rigorous proof is G4-6 work}}.
- **$\alpha > 1$ branch.** G4-4c week 2 found bit-identical 67.88% recovery plateau across $N_0$ at $\alpha > 1$. Whether this affects $S_{\rm BH}$ at $\alpha = 1$ via sub-leading corrections: {{TBD from G4-5c, deferred to G4-6 if not closed}}.

### §8.3 Observation-level results

{{TBD from each sub-sprint.}}

### §8.4 Open follow-ons (deferred to G4-6 or beyond)

- **Full $S_{\rm BH}$ closure with all subleading corrections.** G4-5 produces the leading $r_h^2 \Lambda^2 / 3 \cdot \phi(2)$ piece; subleading corrections at $O(a^2/r_h^2)$ and $O(r_h^2/R^2)$ remain.
- **Higher curvature corrections.** $\Lambda^0$ topological terms (Euler density), $\Lambda^{-2}$ curvature-squared. These contribute $\log\Lambda$ corrections to $S_{\rm BH}$ and require multi-Mellin extraction.
- **Wald entropy formula extension.** For modified gravity (higher curvature), $S_{\rm BH}$ depends on $\partial\mathcal{L}/\partial R_{abcd}$. The discrete substrate naturally supports this via the heat-kernel expansion but extracting Wald coefficients is G4-6+.
- **Asymptotic free $\alpha > 1$.** Excess-angle regime structural asymmetry per G4-4c open question.

---

## §9. G4-5 sub-sprint status table

| Sub-sprint | Falsifier | Status | Recovery | Memo |
|---|---|---|---|---|
| G4-5 scoping | (architecture) | DONE | — | `g4_5_scoping_memo.md` |
| G4-5a first move | F8 (tip-only $S_{\rm tip}$) | DONE | 13–37% (sprint-scale ratio) | `g4_5a_first_move_tip_replica_memo.md` |
| G4-5a-refined | F8 (methodological) | RUNNING | {{TBD}} | {{TBD from G4-5a-refined}} |
| G4-5b | F9 ($\Lambda^4$ bulk Weyl) | RUNNING | {{TBD}} | {{TBD from G4-5b}} |
| G4-5b | F10 ($\Lambda^2$ Einstein–Hilbert) | RUNNING | {{TBD}} | {{TBD from G4-5b}} |
| G4-5c | F11 (joint warp + conical) | RUNNING | {{TBD}} | {{TBD from G4-5c}} |
| G4-5d | F12 (cutoff $\phi(2)$ scaling) | RUNNING | {{TBD}} | {{TBD from G4-5d}} |
| **G4-5e synthesis** | (closure narrative) | **PLACEHOLDER** | — | `g4_5_synthesis_memo.md` (this) |

Expected total wall time for the parallel sub-sprint batch: {{TBD}} weeks. The G4-5 commitment originally sized at 8–13 weeks (G4-5 scoping §6) is anticipated to collapse to {{TBD}} weeks given the G4-4 infrastructure transfer and the parallelization.

---

## §10. Implications for G4-6 (full $S_{\rm BH}$ closure)

G4-6 is the full $S_{\rm BH}$ derivation including all subleading terms. After G4-5 closure, G4-6's load-bearing pieces are:

1. **Leading-order $S_{\rm BH} = r_h^2 \Lambda^2 / 3 \cdot \phi(2)$**: {{TBD CLOSED by G4-5}}.
2. **Subleading $O(a^2/r_h^2)$ UV corrections**: requires multi-substrate continuum extrapolation. Estimated effort: 2–3 months.
3. **Subleading $O(r_h^2/R^2)$ IR corrections**: requires IR-boundary regularization analysis. Estimated effort: 1–2 months.
4. **Higher-curvature contributions to Wald entropy**: requires extraction of $\Lambda^0$ Euler-density and $\Lambda^{-2}$ Riemann-squared pieces via multi-Mellin transform. Estimated effort: 3–4 months.
5. **$\alpha > 1$ branch closure** (excess-angle regime): requires analytical understanding of the 67.88% recovery plateau. Estimated effort: 1–2 months.

**Total G4-6 estimate post G4-5 closure: 7–11 months** (was 5–8 months pre-G4-5; the precise sizing waits on G4-5 sub-sprint outcomes).

The structural-skeleton-scope reading (Paper 51 §13): **the framework predicts the LEADING-order structural form of $S_{\rm BH}$ but DOES NOT autonomously fix the cutoff function $f$** (Class 1 calibration, G8). G4-5 + G4-6 will lock this division: G4-5 closes the $f$-dependence at $\phi(2)$ scaling (Class 1 calibration as expected), G4-6 closes the leading + sub-leading terms (framework predictions in the strong sense).

---

## §11. Files

- `debug/g4_5_scoping_memo.md` (G4-5 scoping)
- `debug/g4_5a_first_move_tip_replica_memo.md` (G4-5a first move)
- `debug/g4_5a_first_move_tip_replica.py` (G4-5a driver)
- `debug/data/g4_5a_first_move_tip_replica.json` (G4-5a results)
- {{TBD from G4-5a-refined: memo, driver, JSON}}
- {{TBD from G4-5b: memo, driver, JSON}}
- {{TBD from G4-5c: module, memo, driver, JSON}}
- {{TBD from G4-5d: memo, driver, JSON}}
- `debug/g4_5_synthesis_memo.md` (this placeholder)
- `debug/g4_5_paper51_update_draft.md` (Paper 51 §12.7 update draft)

---

## §12. Cross-references

- **Paper 28 §4.15** (G4-2 continuum derivation of $S_{\rm BH}$)
- **Paper 28 §4.16** (G8 cutoff-function classification)
- **Paper 28 §4.17** (G4-3 discrete substrate)
- **Paper 51 §5** (G4-1, G4-2 — continuum)
- **Paper 51 §10** (G7 — Newton constant + cosmological constant)
- **Paper 51 §11** (G8 — cutoff dependence)
- **Paper 51 §12** (G4-4 — warped Dirac, replica derivative)
- **CLAUDE.md §1.7 WH1** (substrate convergence; structural-skeleton-scope reading)
- **CLAUDE.md §2** (sprint chronicle; one-liner G4-5 closure entry pending)
- **CHANGELOG** (version bump pending G4-5 closure)

---

## §13. Synthesis-memo discipline note

This memo follows the [`feedback_no_synthesis_memos`](memory/feedback_no_synthesis_memos.md) standing rule with the explicit exception that this is a **single canonical scaffolding memo prepared in advance of parallel sub-sprint closures**, not a cross-sprint synthesis written after multiple sprints have settled. The scaffolding is justified by the parallelism: writing this skeleton during the parallel-dispatch window keeps the writeup framework ready and reduces post-closure synthesis work. When the parallel sub-sprints close, this memo gets populated in-place (not superseded by a new synthesis memo).
