# Sprint G4-5 scoping — Discrete replica method for $S_{\rm BH}$

**Date:** 2026-05-29
**Path:** Gravity arc, opening of the multi-month G4-5 commitment ($\sim 2$–4 months) built on top of the G4-4 sprint-scale closure (this session). G4-5 implements the **replica-method derivation of $S_{\rm BH}$** on the discrete substrate, integrating the wedge-Dirac heat trace $K(\alpha, t)$ over $t$ with a Connes–Chamseddine cutoff function.
**Verdict:** **POSITIVE-SCOPING-G4-5.** Architecture well-defined; load-bearing structural quantities already validated at sprint scale by G4-4c (tip coefficient $-1/12$ bit-exact to 5 digits) and G4-4f (replica derivative $+1/6$ at 96.69%). Sub-sprint sequence G4-5a/b/c/d/e named with five falsifiers and effort estimates. **Multi-month $S_{\rm BH}$ derivation has a quantitatively validated foundation** — remaining work is integration over $t$ with cutoff regularization, not foundational uncertainty.

## §1. Context

The Connes–Chamseddine spectral action for gravity is
$$
S_{\rm CC}[D, \Lambda] = \mathrm{Tr}\, f(D^2/\Lambda^2)
$$
with cutoff function $f$. In the heat-kernel representation, this is
$$
S_{\rm CC} = \int_0^\infty \frac{dt}{t} \tilde{f}(t \Lambda^2)\, K(t)
$$
where $\tilde{f}$ is related to $f$ by inverse Mellin transform. The Seeley–DeWitt expansion of $K(t)$ at small $t$ gives a $\Lambda^4$ cosmological term + $\Lambda^2$ Einstein–Hilbert term + $\Lambda^0$ higher curvature terms (G7, this paper §10).

For the conical-defect $\alpha$-deformed cigar, the **Euclidean action**
$$
I_E(\alpha) = -\frac{1}{2} S_{\rm CC}(\alpha)
= -\frac{1}{2}\int_0^\infty \frac{dt}{t} \tilde{f}(t\Lambda^2)\, K^{\rm Dirac}_{\rm cigar}(\alpha, t)
$$
gives the **BH entropy** via the replica formula
$$
S_{\rm BH} = -\frac{dI_E}{d\alpha}\bigg|_{\alpha = 1}.
$$

The G4-4c + G4-4f sprint-scale identifications give:
1. Tip coefficient: $\Delta_K^{\rm Dirac, tip}(\alpha) = -\tfrac{1}{12}(1/\alpha - \alpha)$ (bit-exact to 5 digits)
2. Replica derivative: $d \Delta_K / d\alpha|_{\alpha = 1} = +1/6$ (96.69% recovery)

What G4-5 must do: integrate these structural pieces over $t$ with the cutoff function to produce $S_{\rm BH}$.

## §2. Continuum prediction

Per Paper 28 §4.15 / G4-2, the continuum derivation gives
$$
S_{\rm BH}^{\rm continuum} = \frac{A_{\rm horizon} \Lambda^2}{12 \pi}
$$
where $A_{\rm horizon} = 4\pi r_h^2$ is the area of the spatial $S^2$ slice at the cigar tip. Substituting:
$$
S_{\rm BH}^{\rm continuum} = \frac{4\pi r_h^2 \cdot \Lambda^2}{12\pi} = \frac{r_h^2 \Lambda^2}{3}.
$$

For a Gaussian cutoff $f(x) = e^{-x}$ (Mellin moments $\phi(s) = \Gamma(s)$): $\phi(2) = 1$, $\phi(1) = 1$. Identifying with Eq.~(eq:replica\_S\_BH) of the Paper 51 G4-4 section requires integrating the heat-kernel expansion against the cutoff and extracting the $\Lambda^2$ piece — exactly what G4-5 does.

For the discrete substrate at $r_h = 2$, $\Lambda = 1$: $S_{\rm BH} = 4/3 \approx 1.33$.

## §3. Architectural inputs from G4-4

All G4-4 infrastructure transfers cleanly:
- `DiscreteWedgeDirac` (G4-4c): wedge spinor heat trace at $\alpha \ne 1$
- `WarpedDiracConstant` (G4-4a): constant-warp cigar heat trace
- `VariableWarpDirac` (G4-4b): variable-warp cigar heat trace
- `verify_F4_tip_regular`, `verify_F6_riemannian_limit`: load-bearing falsifiers
- Sweet-spot $t$ window identification: $t \in [0.5, 5]$ for clean extraction (per G4-4f cross-$t$ analysis)
- Substrate convention: $R = 10$, $a = 0.05$, $N_0 = 120$ (per G4-4c week 3 sweet spot)

## §4. Three computational ingredients

### §4.1 The heat-trace data
$$
K^{\rm Dirac}_{\rm wedge \times S^2}(\alpha, t; r_h)
= \text{(joint wedge-Dirac on $D^2_\alpha \times S^2_{r_h}$)}.
$$

For the **constant-warp** case (G4-4a + G4-4c product), this factorizes:
$$
K^{\rm Dirac}_{\rm wedge \times S^2}(\alpha, t; r_h)
= K^{\rm Dirac}_{\rm wedge}(\alpha, t) \cdot K^{\rm Dirac}_{S^2}(t; r_h)
$$
where the second factor is the Camporesi–Higuchi $S^2$ spinor heat trace (already implemented).

For the **variable-warp** case (G4-4b), the factorization breaks at $O((r'/r)^2)$ and we need the full `VariableWarpDirac.heat_trace`. Defer this to G4-5c.

### §4.2 The cutoff function and Mellin moments
Standard choices:
- **Gaussian**: $f(x) = e^{-x}$, $\phi(s) = \Gamma(s)$, $\phi(1) = 1$, $\phi(2) = 1$
- **Sharp**: $f(x) = \Theta(1 - x)$, $\phi(s) = 1/s$, $\phi(1) = 1$, $\phi(2) = 1/2$
- **Polynomial**: $f(x) = e^{-x^2}$, $\phi(s) = \Gamma(s/2)/2$

Per G8 (this paper §10), cutoff choice is Class 1 calibration data. G4-5 should reproduce the cutoff-dependent prediction
$$
S_{\rm BH}(f) = \frac{A_{\rm horizon} \Lambda^2}{12\pi} \cdot \phi(2).
$$
(with appropriate normalization).

### §4.3 The integration
$$
I_E(\alpha) = -\frac{1}{2} \int_{t_{\rm IR}}^{t_{\rm UV}} \frac{dt}{t}\, \tilde{f}(t\Lambda^2)\, K^{\rm Dirac}_{\rm cigar}(\alpha, t).
$$

Cutoffs $t_{\rm IR}$ and $t_{\rm UV}$ are dictated by the discrete substrate ($t_{\rm UV} \sim a^2$ from the lattice spacing, $t_{\rm IR} \sim R^2$ from the IR boundary). The cutoff function $\tilde{f}$ regularizes both ends.

## §5. Load-bearing falsifiers (G4-5 first move)

### F8 — Replica-method tip-only integration

The **tip contribution** is the load-bearing entropy piece:
$$
S_{\rm tip} = -\frac{1}{2} \int_0^\infty \frac{dt}{t}\, \tilde{f}(t \Lambda^2) \cdot \frac{d \Delta_K^{\rm Dirac, tip}}{d\alpha}\bigg|_{\alpha = 1}
$$
$$
= -\frac{1}{2} \cdot \frac{1}{6} \cdot \int_0^\infty \frac{dt}{t} \tilde{f}(t\Lambda^2)
$$
$$
= -\frac{1}{12} \cdot M_0(\tilde{f}; \Lambda^2)
$$
where $M_0$ is the log-Mellin moment (a $1/\epsilon$ divergence regulated by $\Lambda$). For a Gaussian Mellin moment + the cutoff $\tilde{f}$ matching $\phi(s)$ convention, the integral collapses to a $\log(\Lambda)$-divergent expression that, after proper renormalization, contributes a constant per topological tip.

**Test**: discrete computation of the LHS via integration of `DiscreteWedgeDirac` data; verify against the continuum tip-entropy contribution per Cheeger–Simons.

### F9 — Bulk Weyl Λ⁴ extraction

The $K_{\rm disk}$ (linear-in-$\alpha$ bulk) part gives the cosmological-constant divergence at $\Lambda^4$:
$$
S_{\rm Weyl} \sim \Lambda^4 \cdot \frac{A_{D^2}}{16\pi^2} \cdot \phi(2)
$$
**Test**: extract this from discrete-substrate $K_{\rm disk}(t)$ via Mellin transform; reproduce continuum.

### F10 — Λ² Einstein–Hilbert extraction

The boundary $a_1$ term gives the $\Lambda^2$ Einstein–Hilbert contribution. **Test**: extract via $K_{\rm disk}(t) - K_{\rm Weyl}(t)$ at small $t$; reproduce continuum.

### F11 — Joint warp + conical defect

For the **physical cigar** (variable warp + conical defect at the horizon), $S_{\rm BH}$ from the joint geometry. **Test**: bit-exact factorization at constant warp (F6 extension); identify the $r_h$-dependence as $r_h^2 \Lambda^2/3$.

### F12 — Cutoff-function dependence

Verify the G8 prediction: $S_{\rm BH}(f) \propto \phi(2)$, with $\phi(2)$ varying by cutoff choice. **Test**: compute $S_{\rm BH}$ for Gaussian, sharp, polynomial cutoffs; verify $\phi(2)$ ratios.

## §6. Sub-sprint sequence

| Sub-sprint | Scope | Effort |
|---|---|---|
| **G4-5a first move** | F8 tip-only replica integration on constant-warp cigar | 1-2 weeks |
| G4-5b | F9 + F10 bulk Weyl extraction (Λ⁴ + Λ²) | 2-3 weeks |
| G4-5c | F11 joint warp + conical defect (variable warp) | 3-4 weeks |
| G4-5d | F12 cutoff-function dependence | 1-2 weeks |
| G4-5e | Synthesis: compare to continuum $S_{\rm BH}$; closure narrative | 1-2 weeks |

**Total G4-5 commitment: 8-13 weeks** (sprint-scale within multi-month). Less than the original 2-4 month estimate because G4-4 has already validated the load-bearing structural quantities.

## §7. First-move plan (G4-5a)

### §7.1 Module structure

New driver `debug/g4_5a_first_move_tip_replica.py` (~250 lines). No new production module needed at first move — use the integration of G4-4f data with cutoff function.

### §7.2 Numerical recipe

For sweet-spot panel ($N_\rho = 200$, $a = 0.05$, $N_0 = 120$, $r_h = 2$):

1. Compute $K_{\rm wedge}^{\rm Dirac}(\alpha, t)$ at $\alpha = 1 \pm \varepsilon$ via `DiscreteWedgeDirac`, $\varepsilon = 1/120$
2. Compute the derivative $dK/d\alpha|_{\alpha = 1}$ at each $t$ on a grid $t \in \{0.1, 0.2, 0.5, 1, 2, 5, 10\}$
3. Subtract the bulk: $\Delta'(t) = dK/d\alpha - K_{\rm disk}$
4. Integrate: $J = \int_{t_{\rm min}}^{t_{\rm max}} (dt/t) \tilde{f}(t \Lambda^2) \Delta'(t)$
5. Compare $S_{\rm BH} = -J/2$ against the continuum prediction

### §7.3 Falsifier exit gate

- $S_{\rm BH}^{\rm discrete}$ within 10% of the continuum prediction
- Cutoff-dependence consistent with $\phi(2)$ scaling

### §7.4 What this DOES NOT do (deferred to subsequent sub-sprints)

- Variable-warp + conical-defect joint Dirac (G4-5c)
- Sub-leading $\alpha^0, \alpha^{-2}$ contributions to $S_{\rm BH}$
- Full Seeley–DeWitt + Mellin integration analytical derivation
- α > 1 branch open question (inherited from G4-4c week 2)

## §8. Open questions inherited from G4-4

1. **α > 1 branch structural asymmetry**: at $\alpha > 1$, the spinor SC recovery plateaus at 67.88% across $N_0$. Does this affect the replica derivative at $\alpha = 1$? G4-4f showed 97% recovery at $\alpha = 1 \pm 0.1$ which is in the cleaner regime, but extending to wider $\alpha$ may surface the asymmetry.

2. **Joint variable-warp + conical-defect Dirac**: G4-4 separated these axes; G4-5c brings them together.

3. **Cutoff-function calibration**: per G8, cutoff is external Class 1 data. G4-5 should make the $f$-dependence explicit but not solve it.

## §9. Honest scope

This is a scoping memo. It does NOT contain new computational results. It:
- Confirms G4-4c + G4-4f outputs are sufficient input
- Names the architectural extension G4-5 needs (integration over $t$ with cutoff)
- Names five load-bearing falsifiers F8 through F12
- Sizes sub-sprint sequence at 8-13 weeks total
- Identifies G4-5a first move plan

The G4-5 multi-month commitment is **structurally additive on top of G4-4's closure** rather than starting from scratch. The replica derivative is already extracted at sprint scale; what remains is its integration against a cutoff function.

## §10. Files
- `debug/g4_5_scoping_memo.md` (this)
- (Future, when G4-5a launches): driver, JSON, memo
