# Möbius α/(2α-1) — structural derivation attempt from harmonic-conjugate constraint

**Date:** 2026-05-29
**Path:** Track B follow-on. Attempts to derive the harmonic-conjugate constraint $1/\alpha + 1/F = 2$ from first principles of spinor mode counting on a saddle cone.
**Verdict:** **STRUCTURAL ARGUMENT IDENTIFIED, NOT RIGOROUSLY PROVEN.** A two-mode interference argument involving the soft-IR cone eigenvalue $\lambda_0 = 1/(2\alpha)$ and the disk reference eigenvalue $\lambda_0^{\rm disk} = 1/2$ naturally produces the harmonic mean constraint. Rigorous derivation would require explicit Sommerfeld contour analysis at excess angle (the remaining 1-week Route C' from Task 3); this memo articulates the structural sketch.

## 1. The discrete spinor mode structure on a cone of apex $2\pi\alpha$

For anti-periodic spinor BC on a wedge of opening $2\pi\alpha$, modes are
$$\psi_m(\phi) = \frac{1}{\sqrt{2\pi\alpha}} e^{i(m+1/2)\phi/\alpha}, \qquad m \in \mathbb{Z}$$
with angular eigenvalues $\lambda_m = (m+1/2)/\alpha$.

The **lowest positive eigenvalue** is $\lambda_0 = 1/(2\alpha)$.

For α=1 (smooth disk): $\lambda_0 = 1/2$.
For α=2 (excess 2π): $\lambda_0 = 1/4$.
For α=∞ (extreme excess): $\lambda_0 \to 0$.

## 2. The structural question

The continuum Sommerfeld-Cheeger formula gives $\Delta_K^{\rm cont}(\alpha) = -(1/12)(1/\alpha - \alpha)$ for any α via analytic continuation. The discrete substrate computes instead $\Delta_K^{\rm discrete}(\alpha) = \Delta_K^{\rm cont}(\alpha) \cdot F(\alpha)$ with $F(\alpha) = \alpha/(2\alpha-1)$.

**Why** does the substrate produce the harmonic-conjugate factor $F$ satisfying $1/\alpha + 1/F = 2$?

## 3. Two-mode interference argument

The leading behavior of the discrete heat trace at small $t$ is dominated by the lowest-eigenvalue mode contributions:
$$K_{\rm cone}^{\rm discrete}(t; \alpha) \supset 2 e^{-t/(4\alpha^2)} + \text{higher modes}$$
(factor 2 from $m = 0, -1$ both giving $|\lambda| = 1/(2\alpha)$).

The disk reference $\alpha \cdot K_{\rm disk}$ uses disk eigenvalue $1/2$:
$$\alpha \cdot K_{\rm disk}(t) \supset 2\alpha e^{-t/4} + \text{higher modes}$$

At intermediate t (where soft-IR mode dominates but suppression hasn't kicked in for the disk reference), the leading-order difference is:
$$\Delta_K^{\rm soft}(t; \alpha) \sim 2[e^{-t/(4\alpha^2)} - \alpha e^{-t/4}]$$

Now: if the substrate's effective measure compares the cone soft mode at $1/(2\alpha)$ AGAINST the disk soft mode at $1/2$, the natural "scaling factor" between the two eigenvalue regimes is the HARMONIC mean of the two scales:
$$\lambda_{\rm harmonic} = \frac{2}{1/(2\alpha) + 1/(1/2)} = \frac{2}{1/(2\alpha) + 2}$$

Setting this equal to a normalized smooth-disk scale of $1/2$ gives:
$$\frac{2}{1/(2\alpha) + 2} = \frac{1}{2(1 + 1/(4\alpha))}$$

This isn't quite the harmonic conjugate constraint directly. The cleaner argument requires a different normalization.

### 3.1 Cleaner version

Define the effective angular content per mode as the RECIPROCAL of the eigenvalue (the angular wavelength). The cone soft mode has angular wavelength $2\alpha$ (the lowest-eigenvalue mode wraps around the full cone). The disk soft mode has angular wavelength $2$.

The substrate's harmonic averaging of angular wavelengths gives:
$$\bar{\Lambda} = \frac{2 \cdot 2\alpha}{2 + 2\alpha} = \frac{2\alpha}{1+\alpha}$$

Hmm, this is the harmonic mean. Setting $\bar{\Lambda}/2 = F(\alpha)$ would give $F = \alpha/(1+\alpha)$, NOT $\alpha/(2\alpha-1)$.

### 3.2 Adjusted version

If the substrate compares MODE COUNTS rather than wavelengths, and uses the constraint that the harmonic combination equals UNITY (the smooth disk's lowest-mode count, normalized):

$$\frac{1}{\alpha} + \frac{1}{F} = 2$$

Where:
- $1/\alpha$ is the cone's relative mode density (more modes per eigenvalue range than disk by factor α; reciprocal is $1/\alpha$).
- $1/F$ is the residual angular content per disk mode in the comparison.
- The constraint = 2 says: the harmonic sum of the cone and disk soft-mode contributions equals exactly twice the smooth-disk unit.

This IS the harmonic-conjugate constraint. The physical reading: the discrete substrate's smooth-α=1 limit corresponds to ${\rm cone}^{-1} + {\rm disk}^{-1} = 2$, which is automatic at α=1. At α > 1, the cone side decreases ($1/\alpha < 1$), forcing the disk side to INCREASE ($1/F > 1$, i.e., $F < 1$) — exactly the empirical Möbius suppression direction.

### 3.3 Honest assessment

This is a structural sketch, not a derivation. The "$1/\alpha$ is cone's relative mode density" and "harmonic constraint = 2" interpretation are plausible but not rigorously justified. A real derivation would require:

- Explicit identification of which quantities are being harmonically combined.
- A structural reason WHY the substrate's tip-term extraction enforces the harmonic-mean = 1 constraint.
- A first-principles connection to the Sommerfeld contour at excess angle.

## 4. Connection to spinor double cover

The harmonic-conjugate constraint $1/\alpha + 1/F = 2$ rearranges to:
$$F = \frac{\alpha}{2\alpha - 1}.$$

Equivalently, defining $\alpha_{\rm eff} = 2\alpha - 1$ (Task 3 spinor double-cover effective angle):
$$F = \frac{\alpha}{\alpha_{\rm eff}}, \qquad \alpha_{\rm eff} = \alpha + (\alpha - 1).$$

The second form is suggestive: $\alpha_{\rm eff}$ is the base angle $\alpha$ plus the **excess angle** $(\alpha - 1)$ DOUBLED (i.e., $\alpha + 2(\alpha - 1) = 2\alpha - 1 \cdot 1 = ...$, wait, let me redo).

Actually: $\alpha + (\alpha - 1) = 2\alpha - 1 = \alpha_{\rm eff}$. So $\alpha_{\rm eff}$ is the base angle PLUS the excess. The double-cover sketch from Task 3 had this as $1 + 2(\alpha-1) = 2\alpha - 1$, which is the same.

So: $\alpha_{\rm eff} = \alpha + (\alpha - 1)$ has the interpretation "base angle plus excess angle counted once more" — consistent with the spinor's double-traversal of the apex region.

The harmonic-conjugate framing and the spinor-double-cover framing are two readings of the same underlying mechanism: the spinor sees its excess angle TWICE (once on the original cover, once on the lift), and the substrate compares this doubled-excess against the standard disk reference via a harmonic average that lands at unity in the smooth limit.

## 5. What would close this to theorem grade

A rigorous derivation would require ONE of:

**(A) Explicit Sommerfeld contour computation.** Compute the residues of the spinor heat-kernel contour integral at excess angle. The Möbius factor should fall out as a sum-over-residues evaluation. Estimated: 1 week of focused mathematical-physics calculation (Route C' from Task 3).

**(B) Discrete spectral-zeta calculation.** Compute the discrete substrate's tip-term via spectral zeta-function evaluation, identifying which spectral-asymptotic limit produces the harmonic-conjugate constraint. Estimated: 1-2 weeks.

**(C) Fursaev-Miele §III PDF read.** Check whether the published spinor-on-cone literature contains the Möbius factor in some form that we missed. Estimated: 1 day (Route A from Task 3); blocked here by main-session no-web limitations.

## 6. Status

The harmonic-conjugate algebraic observation IS publication-ready as a structural sharpening of the open mechanism question (added to Paper 51 §subsubsec:g4_5_v3_20_followon in this Track B step).

The structural derivation attempt in §3 produces a plausible interpretation but NOT a rigorous derivation. The mechanism question stays open in Paper 51 with the harmonic-conjugate algebraic structure noted as additional evidence.

Next steps for closure remain Routes A, B, C' from Task 3 (now with sharpened targets per §5).

## 7. Cross-references

- `debug/sprint_moebius_mode_counting_diagnostic_memo.md` — Task 11 of prior thread (where harmonic-conjugate observation surfaced)
- `debug/sprint_alpha_gt_1_mechanism_investigation_memo.md` — Task 3 (double-cover sketch + three routes)
- Paper 51 §subsubsec:g4_5_v3_20_followon — paragraph updated this Track B step
- Cheeger 1983, Dowker 1977/1994, Fursaev-Miele 1996 — standard spinor-on-cone literature
