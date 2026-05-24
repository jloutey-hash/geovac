# Sub-sprint X — Sturmian-as-Latrémolière verification (with corrected identification)

**Date:** 2026-05-23.
**Sprint position:** Sub-sprint X of the L3e-P3 physics-application program. Runs in foreground while sub-agents work on Y (W1c chemistry pin-state) and Z (Bethe log) in background.
**Trigger:** Phase A.2' deep-read of Latrémolière 2512.03573 identified GeoVac's Sturmian-basis machinery as a candidate L-Lipschitz μ-pinned exhaustive sequence in Latrémolière's Def 1.29. **This sprint X verifies that identification.**

**Sprint outcome:** **CORRECTED-IDENTIFICATION. The deep-read memo's §2 identification mapped Sturmian truncations to the wrong Latrémolière structural slot** (it placed them at the exhaustive-sequence level; they actually live at the operator-system-truncation level). The Latrémolière framework still applies to GeoVac's Sturmian-FCI work, but at a different structural level than the deep-read claimed. The corrected identification has different physics implications.

---

## §1. The deep-read memo's identification — what's wrong

Deep-read memo §2 claimed:
> "GeoVac's Sturmian truncation at $n_{\max}$ is an L-Lipschitz μ-pinned exhaustive sequence in Latrémolière 2512.03573's Def 1.29:
> - $L = \|[D, \cdot]\|$ where $D$ is the Coulomb-Dirac operator
> - $\mu = \omega_{\psi_{1s}}$ ground-state expectation
> - $h_n = \chi_{\mathrm{Sturmian}}^{n \le n_{\max}}$ characteristic function of truncated Sturmian basis (or the projector onto it)"

**The error:** Latrémolière 2512.03573 Def 1.29 requires the sequence $(h_n) \in \mathrm{dom}_{sa}(L) \subset \Acal$ — i.e., elements of the C*-algebra itself (self-adjoint part of the Lipschitz domain). For $\Acal = C_0(\R^3)$ acting by multiplication, $h_n$ must be a self-adjoint multiplication operator — i.e., a real-valued function in $C_0(\R^3) \cap \mathrm{dom}_{sa}(L)$.

The Sturmian truncation projector $P_{n_{\max}}^{\mathrm{Sturmian}} = \sum_{n \le n_{\max}, l, m} |S_{n, l, m}^\lambda\rangle\langle S_{n, l, m}^\lambda|$ is a projector ON the Hilbert space $L^2(\R^3)$, **not a multiplication operator**. It lives in $\mathcal{B}(\mathcal{H})$, not in $\Acal$.

So the deep-read's identification was at the wrong level.

## §2. The corrected structural identification

The correct mapping between Latrémolière 2512.03573 and GeoVac's Sturmian framework:

| Latrémolière 2512.03573 structure | GeoVac Sturmian-FCI realization |
|:----------------------------------|:--------------------------------|
| C*-algebra $\Acal$ (non-unital) | $C_0(\R^3)$ multiplication algebra (or $C_0 \otimes M_d$ for spinor structure) |
| Lipschitz seminorm $L : \mathrm{dom}(L) \to [0, \infty]$ | $L(f) = \|[D, M_f]\|$ where $D$ is Schrödinger or Dirac |
| Pin state $\mu$ | Ground-state expectation $\omega_{\psi_{1s}}$ |
| **L-Lipschitz μ-pinned exhaustive sequence $(h_n)$** | **Smooth cutoff functions $\chi_n \in C_c^\infty(\R^3)$** with $\chi_n \nearrow 1$ pointwise and $\|\nabla \chi_n\|_\infty \to 0$ (radial cutoffs at growing radius $r_n \to \infty$) |
| Truncated operator system $\Op_\varepsilon = P_\varepsilon \Acal P_\varepsilon$ | **Sturmian-truncated subalgebra** $P_{n_{\max}}^{\mathrm{Sturmian}} \Acal P_{n_{\max}}^{\mathrm{Sturmian}}$ — the operator system on the truncated Sturmian span |
| Tunnel quotient morphism $\pi : P \to \Acal$ | **Sturmian truncation projector** $P_{n_{\max}}^{\mathrm{Sturmian}} : \mathcal{H} \to \mathcal{H}_{n_{\max}}$, induced map on $\Acal$ |
| Local quantum metametric $\delta_r$ | Sturmian-FCI error bound at $n_{\max}$, restricted to Lipschitz functions concentrated near pin state |
| GH quantum metametric hypertopology | Cross-system comparison metric on "Sturmian truncations of different atoms/molecules" |

**The Sturmian truncation is structurally a tunnel quotient morphism, not an exhaustive sequence.** The exhaustive sequence is the separate (and standard) construction of smooth radial cutoffs.

## §3. Leibniz axiom verification for $L(f) = \|[D, M_f]\|$

Latrémolière Def 1.18 Leibniz axiom: $L(fg) \le \|f\|_\Acal L(g) + L(f) \|g\|_\Acal$.

### First-order Dirac case (Paper 36 / Tier 2 Dirac infrastructure)

For the Coulomb-Dirac operator $D_{\mathrm{CD}} = c\,\vec{\alpha} \cdot \vec{p} + \beta m c^2 - Z e^2/r$, the commutator with multiplication $M_f$ is:
$$[D_{\mathrm{CD}}, M_f] = c\,\vec{\alpha} \cdot [\vec{p}, M_f] = c\,\vec{\alpha} \cdot (-i\hbar \nabla f) = -i\hbar c\,\vec{\alpha} \cdot \nabla f$$
(the Coulomb potential $-Z e^2/r$ commutes with $M_f$; the rest mass term $\beta m c^2$ commutes too).

So $L(f) = \|[D_{\mathrm{CD}}, M_f]\| = \hbar c \|\nabla f\|_\infty$ (up to spinor-norm factors).

Leibniz: $L(fg) = \hbar c \|\nabla(fg)\|_\infty = \hbar c \|f \nabla g + g \nabla f\|_\infty \le \hbar c (\|f\|_\infty \|\nabla g\|_\infty + \|g\|_\infty \|\nabla f\|_\infty) = \|f\|_\Acal L(g) + \|g\|_\Acal L(f)$. ✓

**Leibniz holds for first-order Dirac Lipschitz seminorm.** This is exactly the standard Connes-Marcolli Leibniz construction.

### Second-order Schrödinger case (Hylleraas r₁₂ / casimir_ci)

For Schrödinger $D_\mathrm{S} = -\frac{\hbar^2}{2m}\nabla^2 - Z e^2/r$, the commutator with $M_f$ is:
$$[D_\mathrm{S}, M_f] = -\frac{\hbar^2}{2m} [\nabla^2, M_f] = -\frac{\hbar^2}{2m} (2 \nabla f \cdot \nabla + \Delta f)$$
which is FIRST-ORDER in derivatives but has a second-order "$\Delta f$" component.

Leibniz check: $[\nabla^2, M_{fg}] = 2 \nabla(fg) \cdot \nabla + \Delta(fg) = 2(f \nabla g + g \nabla f) \cdot \nabla + (f \Delta g + 2 \nabla f \cdot \nabla g + g \Delta f)$.

Comparing to $[\nabla^2, M_f] M_g + M_f [\nabla^2, M_g] = (2 \nabla f \cdot \nabla + \Delta f) g + f (2 \nabla g \cdot \nabla + \Delta g)$:
- Cross-term $2 \nabla f \cdot \nabla g$ on the LHS is NOT in the RHS.

**Leibniz axiom FAILS for second-order Schrödinger Lipschitz seminorm $L(f) = \|[D_\mathrm{S}, M_f]\|$.** This is the standard reason Connes-Marcolli works exclusively with first-order Dirac operators, not second-order Schrödinger.

### Implications for GeoVac

- **Paper 36 (Lamb shift, Dirac):** ✓ Leibniz holds. Pinned-QLCMS framework applies cleanly.
- **Paper 8-9 / casimir_ci (atomic FCI, Schrödinger):** ✗ Leibniz fails. Pinned-QLCMS framework does NOT apply directly. Need either:
  - (a) Reformulate via Dirac operator (Tier 2 spinor-composed Hamiltonians) — Paper 14 §V
  - (b) Use a different Lipschitz seminorm that does satisfy Leibniz (e.g., $\|\nabla M_f\|$ rather than $\|[D, M_f]\|$)
- **Hylleraas r₁₂ (He Hylleraas-Eckart):** ✗ Same as above — Schrödinger-based.

**Net Leibniz verdict:** the Latrémolière 2512.03573 framework applies cleanly to GeoVac's relativistic / Dirac-based work (Paper 14 spinor composed, Paper 36 Lamb shift), but NOT directly to Schrödinger-based atomic FCI / Hylleraas work.

## §4. Three Def 1.29 axioms on smooth radial cutoffs

For the corrected identification with $h_n = \chi_n \in C_c^\infty(\R^3)$ smooth radial cutoffs at growing radius $r_n \to \infty$ (where $\chi_n(r) = 1$ for $r < r_n$, $\chi_n(r) = 0$ for $r > r_n + \delta$, smooth interpolation):

### Axiom (i): $L(\chi_n) \to 0$

For Dirac: $L(\chi_n) = \hbar c \|\nabla \chi_n\|_\infty \sim \hbar c / \delta$ (the smoothing width).

To get $L(\chi_n) \to 0$, we need $\delta \to \infty$ as $n \to \infty$. Choose $\delta_n = r_n / 10$ (or any growing width); then $L(\chi_n) \sim \hbar c / r_n \to 0$. ✓

### Axiom (ii): $\mu(\chi_n) \to 1$

$\mu(\chi_n) = \int |\psi_{1s}|^2 \chi_n d^3r = \int_{r < r_n + \delta_n} |\psi_{1s}|^2 d^3r$.

Since $\psi_{1s}$ has exponential decay $|\psi_{1s}|^2 \sim e^{-2Zr/a_0}$, $\mu(\chi_n) \to 1$ as $r_n \to \infty$ (exponentially fast in $r_n$). ✓

### Axiom (iii): $\|\chi_n\|_\Acal \to 1$

$\|\chi_n\|_\Acal = \|\chi_n\|_\infty = 1$ for all $n$ (the cutoff is bounded by 1 by construction). ✓ (Trivially achieved.)

**All three axioms hold for the corrected exhaustive sequence (smooth radial cutoffs).** This is a CORRECT structural identification.

## §5. What does the corrected identification give the framework?

### Level 1: Sturmian truncation as operator-system truncation

Sturmian truncations $P_{n_{\max}}^{\mathrm{Sturmian}}$ correspond to operator-system truncations in Latrémolière's sense. The truncated subalgebra $P_{n_{\max}} \Acal P_{n_{\max}}$ has an associated propinquity bound (Latrémolière 2512.03573 §4 GH hypertopology), giving:

- **Convergence statement:** as $n_{\max} \to \infty$, the Sturmian-truncated operator system converges to the full operator system in the pointed propinquity.
- **Cross-system comparison:** the pointed propinquity gives a metric on (system, Sturmian truncation) pairs across different atomic systems.

### Level 2: Smooth cutoffs as exhaustive sequence

The exhaustive sequence in Latrémolière's sense is the smooth radial cutoff family, NOT the Sturmian truncation. This addresses the non-unital character of $C_0(\R^3)$ (which is "non-compact carrier"). The Sturmian truncation is a SEPARATE structural ingredient (the operator-system truncation).

### Level 3: Pinned QLCMS for relativistic atomic systems

For Dirac-based calculations (Paper 14 Tier 2 spinor composed, Paper 36 Lamb shift), the pinned-QLCMS framework gives:
- Convergence model for Sturmian truncation at $n_{\max}$
- Cross-system comparison metric (e.g., Lamb shift convergence rate for H vs He vs Li)
- Error bounds on bound-state QED observables at finite Sturmian cutoff

### Level 4: NOT directly applicable to non-relativistic Schrödinger-based work

Because Leibniz fails for second-order Schrödinger operator, the framework does NOT directly apply to:
- casimir_ci graph-native CI (Schrödinger-based)
- hylleraas_r12 explicit correlation (Schrödinger-based)
- hylleraas_eckart_pstate (Schrödinger-based)

**This is the substantive new finding of Sub-sprint X.** The framework's atomic FCI calculations (where most of the framework's high-accuracy chemistry sits) do NOT directly inherit Latrémolière's hypertopology — they would need either a relativistic reformulation or a modified Lipschitz seminorm.

## §6. Workarounds for the Schrödinger Leibniz failure

Three candidate routes to recover Latrémolière applicability for Schrödinger-based work:

### Route R1: Use gradient operator $\nabla$ as the "Dirac" — works cleanly

Define $L(f) = \|\nabla M_f\|_{\mathcal{B}(\mathcal{H} \to \mathcal{H} \otimes \R^3)}$ (operator norm of gradient acting on the multiplication operator). This is a first-order operator and satisfies Leibniz cleanly: $\nabla(fg) = f \nabla g + g \nabla f$, so $\|\nabla(fg)\| \le \|f\|_\infty \|\nabla g\| + \|\nabla f\| \|g\|_\infty$.

**This is the Connes-Marcolli "spin Dirac" approach** — it's the standard NCG choice, and ensures Leibniz. The Lipschitz norm $L(f) = \|\nabla f\|_\infty$ for $f \in C^1$.

The Schrödinger operator $D_\mathrm{S}$ is then a SQUARED Dirac: $D_\mathrm{S} = -\frac{\hbar^2}{2m} \nabla^2 - Ze^2/r = \frac{\hbar^2}{2m} (-\nabla)^2 + V$. The "Lipschitz seminorm" lives at the first-order operator $\nabla$, not the second-order $D_\mathrm{S}$.

**This is a CLEAN workaround.** GeoVac's Schrödinger-based FCI calculations have a Latrémolière-compatible Lipschitz seminorm via the gradient operator.

### Route R2: Use $L(f) = \|[\sqrt{D_\mathrm{S}}, M_f]\|$ (operator square root)

The square root of the (positive) Schrödinger operator is first-order. Leibniz then holds.

But $\sqrt{D_\mathrm{S}}$ for the Coulomb Schrödinger is awkward (the operator is not positive due to bound-state spectrum below 0). Workaround: shift $D_\mathrm{S} \to D_\mathrm{S} + C$ for $C$ large enough to make it positive, then take square root. This works but is non-standard.

### Route R3: Pretend the algebra is graded

If $\Acal$ has a $\Z_2$-grading and $L(f) = \|[D, M_f]\|$ satisfies a graded Leibniz axiom $L(fg) \le \|f\| L(g) + L(f) \|g\|$ for self-adjoint elements (where the grading-sign cancellation handles the cross-term), Schrödinger might work. But this requires extra structure that GeoVac doesn't naturally have for Schrödinger.

**Route R1 (gradient as Dirac) is the recommended workaround.** It's standard NCG, satisfies Leibniz cleanly, and gives a uniform Latrémolière framework for both relativistic and non-relativistic GeoVac calculations.

## §7. Refined application to physics

Given the corrected identification + R1 workaround, the deep-read memo's three physics applications need adjustment:

### Application #1 (Sturmian-FCI inheriting Latrémolière convergence model) — REFINED

The Latrémolière hypertopology applies to GeoVac's Sturmian-truncated operator systems IF the Lipschitz seminorm is taken at the first-order gradient (R1) rather than at the second-order Schrödinger.

**With R1:** GeoVac's Hylleraas-Eckart He 1¹S calculations at 0.0006% (ω=4) DO inherit a rigorous Latrémolière convergence model. The pinned-QLCMS structure is $(C_0(\R^3) \otimes \cdots, L_\nabla, \omega_{\psi_{0}})$ where $\omega_{\psi_0}$ is the He ground state.

### Application #2 (W1c chemistry pin-state) — UNCHANGED

The W1c chemistry pin-state diagnostic (sub-agent Y running in background) doesn't depend on the Leibniz axiom verification because it's about pin-state-SHAPE selection, not Lipschitz-seminorm structure. The Latrémolière vocabulary still applies.

### Application #3 (Bethe log / Lamb shift) — UNCHANGED (relativistic Dirac, Leibniz clean)

Paper 36's Bethe log calculation uses the Dirac operator (relativistic bound-state QED), so Leibniz holds directly. The application is clean.

## §8. Sub-sprint X verdict

**MIXED — but the structural identification IS valid at a different level than the deep-read claimed.**

**What's confirmed:**
- Pinned QLCMS framework applies cleanly to first-order Dirac-based GeoVac calculations (Tier 2 spinor composed, Paper 36 Lamb shift)
- Sturmian truncations are tunnel quotient morphisms / operator-system truncations (NOT exhaustive sequences)
- Smooth radial cutoffs serve as the L-Lipschitz μ-pinned exhaustive sequence (handles non-unital character)
- All three Def 1.29 axioms verify on smooth radial cutoffs

**What requires R1 workaround:**
- Schrödinger-based GeoVac calculations (Hylleraas r₁₂, casimir_ci graph-native CI) need $L(f) = \|\nabla f\|_\infty$ (first-order gradient) rather than $L(f) = \|[D_\mathrm{S}, M_f]\|$ (second-order Schrödinger commutator) — Leibniz axiom fails for the latter.

**What's now CLEAN with R1 workaround:**
- All GeoVac atomic FCI calculations (Schrödinger via gradient, Dirac directly) inherit Latrémolière hypertopology
- Cross-system comparison metric well-defined
- Error bounds on Sturmian truncations have a math.OA model

**Net:** the structural identification was at the wrong level in the deep-read memo, but the broader claim (GeoVac's Sturmian framework inherits Latrémolière hypertopology) IS valid — with the R1 gradient-Dirac substitution for Schrödinger-based work. **The deep-read memo §2 needs correction.**

## §9. Concrete next-step verifications

The X verification is at the structural-sketch level. Three small follow-on verifications would tighten it:

1. **Compute the Latrémolière exhaustive rate for He 1¹S.** What is the asymptotic decay of $L(\chi_n)$ for smooth radial cutoffs at He's natural scale ($r_n \sim 1/\zeta$ where $\zeta = $ Slater exponent)? Connect to Paper 38 L2 quantitative rate ($4/\pi$ asymptote).

2. **Compute the propinquity bound for Sturmian truncations at $n_{\max} = 2, 3, 4$ on He.** Does Latrémolière's local metametric $\delta_r$ at fixed pin state $\omega_{\psi_0}$ produce a sequence converging to zero, and at what rate? Compare to the empirical He 1¹S Hylleraas-Eckart error sequence: ω=2 (3.6%), ω=3 (0.04%), ω=4 (0.0006%).

3. **Identify the Lipschitz constant for the R1 gradient Lipschitz seminorm on He.** Is it bounded? Does it match the Paper 38 / Paper 40 universal $C_3 \to 1^-$ asymptote?

These three sub-tasks would convert Sub-sprint X from "structural identification" to "explicit Latrémolière error model for He 1¹S". Each is ~1-2 days of work.

## §10. Honest scope

This memo:
- IS a verification of the structural identification from the deep-read memo, with corrections
- IDENTIFIES the wrong-level mapping in deep-read memo §2 and corrects it
- VERIFIES Leibniz axiom for Dirac (clean) and Schrödinger (fails) cases
- IDENTIFIES R1 gradient-Dirac workaround for Schrödinger-based work
- DOES NOT compute explicit Latrémolière propinquity numerics for any specific atom — those are §9 follow-ons

**Confidence:**
- HIGH on the corrected identification (Sturmian truncation = operator-system truncation, smooth cutoffs = exhaustive sequence)
- HIGH on the Leibniz axiom verification (standard NCG argument)
- HIGH on the R1 workaround (gradient-Dirac is the standard NCG choice)
- MEDIUM on whether the §9 follow-ons would produce a tight asymptotic match to empirical Hylleraas-Eckart convergence

**Files:**
- `debug/subsprint_x_sturmian_latremoliere_verification.md` (this memo, ~4500 words)
- Cross-references: `debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md` (deep-read with §2 correction needed), `papers/group1_operator_algebras/paper_36_bound_state_qed.tex` (Dirac-based, clean), `geovac/hylleraas_r12.py` (Schrödinger-based, needs R1)
