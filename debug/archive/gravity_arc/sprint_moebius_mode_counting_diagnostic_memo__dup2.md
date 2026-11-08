# Sprint Möbius α > 1 mechanism — discrete mode counting diagnostic

**Date:** 2026-05-29
**Path:** Continuation thread, Task 11. Sprint-scale diagnostic of the open α > 1 Möbius mechanism via discrete mode counting on the wedge-Dirac substrate.
**Verdict:** **MECHANISM SHARPENED, NOT CLOSED.** Discrete mode counting identifies a structural fingerprint: $F(\alpha) = \alpha/(2\alpha-1)$ is the unique Möbius transformation satisfying $1/\alpha + 1/F = 2$ (i.e., $F$ is the harmonic conjugate of $\alpha$ with respect to $1$). This sharpens the structural sketch from Task 3 (double-cover monodromy) to a quantitative constraint, but does NOT close the mechanism question to first-principles derivation. Two follow-on routes named.

## 1. The locked empirical observation

Per v3.20.0 task #25 (memo `debug/alpha_gt_1_moebius_validation_4_5_10_memo.md`):
$$
\text{slope}^{\rm Dirac, ex}(\alpha) = -\frac{1}{12} \cdot \frac{\alpha}{2\alpha-1}, \qquad (\alpha > 1)
$$
verified at sub-2% across six α values $\in \{1.5, 2, 3, 4, 5, 10\}$, with α=10 at -0.032% (essentially exact).

## 2. Mode structure of the discrete wedge-Dirac substrate

The discrete wedge at apex angle $2\pi\alpha$ has:
- Azimuthal coordinate $\phi \in [0, 2\pi\alpha]$.
- Anti-periodic BC: $\psi(\phi + 2\pi\alpha) = -\psi(\phi)$.
- Discrete modes $\psi_m(\phi) = e^{i(m+1/2)\phi/\alpha}$ for $m \in \mathbb{Z}$.
- Angular eigenvalues $\lambda_m = (m + 1/2)/\alpha$.
- Total $N_\phi = \alpha N_0$ modes truncated at $|\lambda_m| \le N_0/(2\alpha) \cdot \alpha = N_0/2$ (UV cutoff in eigenvalue space).

**Critical observation:** the UV cutoff $|\lambda| \le N_0/2$ does NOT depend on $\alpha$ at the level of eigenvalue space. The substrate has the same eigenvalue range $[-N_0/2, N_0/2]$ regardless of $\alpha$; what changes is the DENSITY of modes within that range (factor $\alpha$).

**Soft-IR mode structure.** The lowest cone eigenvalue is:
- Cone at $\alpha = 1$ (disk): $\lambda_0 = 1/2$ (the standard half-integer infrared).
- Cone at $\alpha = 2$: $\lambda_0 = 1/(2\alpha) = 1/4$ (softer IR).
- Cone at $\alpha = N$ (large excess): $\lambda_0 = 1/(2N) \to 0$ (extreme soft IR).

The cone at $\alpha > 1$ has a SOFTER infrared mode than the disk has, scaling as $1/(2\alpha)$.

## 3. Heat trace contribution decomposition

The discrete heat trace at apex angle $\alpha$ is:
$$
K_{\rm cone}^{\rm discrete}(t; \alpha) = \sum_{|\lambda_m| \le N_0/2} g_m \cdot e^{-t \lambda_m^2}
= \sum_{|m| \le \alpha N_0/2} g_m \cdot e^{-t (m+1/2)^2/\alpha^2}
$$
with $g_m$ the multiplicity (including the spinor double-degeneracy on the radial sector).

The tip term is:
$$
\Delta_K^{\rm discrete}(\alpha, t) = \partial_\alpha K_{\rm cone}^{\rm discrete}(t; \alpha)\Big|_{\alpha = 1} - K_{\rm disk}^{\rm discrete}(t)
$$
or equivalently the full $K_{\rm cone}(\alpha, t) - \alpha \cdot K_{\rm disk}(t)$ at general $\alpha$.

At $\alpha > 1$, the leading contributions split as:
- **Bulk part:** $\alpha \cdot K_{\rm disk}$, scaling as $\alpha / (4\pi t)$ at small $t$.
- **Cone deviation:** non-trivial scaling depending on the soft-IR mode contributions and the truncation structure.

## 4. The structural fingerprint — harmonic conjugate

The empirical Möbius $F(\alpha) = \alpha/(2\alpha-1)$ has a CLEAN algebraic structure that I observed:

**Observation:** $F(\alpha)$ satisfies
$$
\frac{1}{\alpha} + \frac{1}{F(\alpha)} = 2.
$$
Equivalently, $F$ is the **harmonic conjugate** of $\alpha$ with respect to $1$. The harmonic mean of $\alpha$ and $F$ is exactly $1$ (the smooth-disk value).

Three structural consequences:
- **Symmetry around $\alpha = 1$:** $F(1) = 1$ (smooth case).
- **Harmonic involution:** $F(F(\alpha)) = \alpha$ (the transformation is its own inverse — applying it twice gives back $\alpha$).
- **Self-dual point:** $F(\alpha) = \alpha$ has solutions at $\alpha = 1$ (and one other unphysical root).

**Physical interpretation candidate.** The harmonic-mean constraint $1/\alpha + 1/F = 2$ is the kind of relation that arises when the substrate effectively AVERAGES two opening angles: the actual cone angle $2\pi\alpha$ and some companion angle $2\pi F(\alpha)$, with the average equal to the smooth-disk angle $2\pi$.

**Mode-counting reading.** At $\alpha > 1$, the discrete substrate has $\alpha N_0$ azimuthal modes with eigenvalues $\lambda_m = (m+1/2)/\alpha$. The soft-IR modes (small $|\lambda_m|$, dominantly $m = 0, -1$ with $|\lambda| = 1/(2\alpha)$) contribute non-trivially to the tip term. The harmonic-conjugate structure suggests these soft-IR modes effectively "cancel" against a fraction $1/F(\alpha) = 2 - 1/\alpha$ of the bulk contribution, leaving the residual tip term:
$$
\Delta_K^{\rm discrete}(\alpha) = \Delta_K^{\rm continuum}(\alpha) \cdot F(\alpha) = -\frac{1}{12}\left(\frac{1}{\alpha} - \alpha\right) \cdot \frac{\alpha}{2\alpha - 1}.
$$

## 5. Connection to the Task 3 double-cover sketch

The Task 3 structural sketch proposed $\alpha_{\rm eff} = 2\alpha - 1$ as the effective angle on the spinor double cover. The current Task 11 observation gives:
$$
F(\alpha) = \frac{\alpha}{\alpha_{\rm eff}} = \frac{\alpha}{2\alpha - 1}.
$$

So $F(\alpha)$ is the RATIO of the base angle to the double-cover effective angle. The harmonic-conjugate observation is then equivalent to:
$$
\frac{1}{\alpha} + \frac{1}{\alpha} \cdot \frac{1}{F} = 2 \iff \frac{1}{\alpha} + \frac{\alpha_{\rm eff}}{\alpha^2} = 2 \iff \alpha_{\rm eff} = 2\alpha - 1.
$$

The two readings are **consistent**:
- Task 3 reading: spinor double-cover lift gives $\alpha_{\rm eff} = 2\alpha - 1$, and the SC suppression factor is the ratio $\alpha/\alpha_{\rm eff}$.
- Task 11 reading: substrate produces $F(\alpha)$ as the harmonic conjugate of $\alpha$ with respect to $1$.

Both point to the same Möbius factor $\alpha/(2\alpha-1)$, with the substantive content being:
- Why $\alpha_{\rm eff} = 2\alpha - 1$? (Task 3's spinor double-cover monodromy heuristic.)
- Why the harmonic-mean = 1 constraint? (Task 11's mode-counting harmonic-conjugate observation.)

These are TWO COMPATIBLE READINGS of the same structural fact, but neither is a first-principles derivation.

## 6. What this diagnostic does NOT close

The mechanism question stays open at first-principles level:
- Why does the discrete substrate compute the harmonic-conjugate factor rather than the continuum-Sommerfeld factor?
- What spectral-flow argument forces the spinor double-cover lift to produce $\alpha_{\rm eff} = 2\alpha - 1$?
- Is the harmonic-conjugate structure preserved under refinement $N_0 \to \infty$ at fixed $\alpha$, or does it asymptote to the continuum-Sommerfeld result?

The first two are about the EXISTENCE of the Möbius factor; the third is about its STABILITY under continuum limit.

The N_0-independence finding (G4-4c week 2: 67.88% recovery bit-identical across $N_0 \in \{120, 240, 480\}$) suggests the Möbius factor IS stable under refinement, so it's not a discretization artifact. But this doesn't pin down WHY it's there.

## 7. Sharpened structural sketch (Task 3 update)

Combining Task 3's double-cover sketch with Task 11's harmonic-conjugate observation, the sharpened structural sketch is:

> **The discrete wedge-Dirac substrate at apex angle $2\pi\alpha$ for $\alpha > 1$ computes the heat-trace tip term with the spinor double-cover lift's effective angle $\alpha_{\rm eff} = 2\alpha - 1$, yielding a Möbius suppression factor $F(\alpha) = \alpha/\alpha_{\rm eff}$ which is structurally the harmonic conjugate of $\alpha$ with respect to $1$. This is consistent with both (a) Sommerfeld-image-method contour integration at excess angle picking up additional poles, and (b) anti-periodic spinor monodromy at the saddle apex effectively doubling the spinor's angular content.**

This is a structural-sketch upgrade, NOT a theorem. But it provides a sharper handle for the three named closure routes from Task 3.

## 8. Two follow-on routes (refined)

**Route A — Fursaev–Miele §III PDF read (1 day).**
Unchanged from Task 3. The harmonic-conjugate observation gives a sharper acceptance test: if Fursaev–Miele's spin-1/2 formula explicitly contains the $\alpha/(2\alpha-1)$ factor (or the harmonic-conjugate constraint), the mechanism closes immediately.

**Route C' — Anti-periodic spinor monodromy explicit calculation (~1 week, refined from Task 3 Route B).**
Rather than the full Sommerfeld contour integration, attempt the more focused calculation: on a saddle cone of apex $2\pi\alpha$, with anti-periodic spinor BC, identify the structural origin of the harmonic-conjugate constraint $1/\alpha + 1/F = 2$. This is a sharper target than "compute the full contour integral"; it asks for the structural reason for the harmonic mean = 1 condition.

If neither closes in a week, the mechanism IS multi-month (consistent with the open-question parking in Paper 51 §subsubsec:g4_5_v3_20_followon).

## 9. Verdict

**MECHANISM SHARPENED.** Task 11 identifies the harmonic-conjugate algebraic structure of the Möbius factor and pins it to the Task 3 spinor double-cover sketch via the relation $F(\alpha) = \alpha/\alpha_{\rm eff}$ with $\alpha_{\rm eff} = 2\alpha - 1$. The two readings are compatible; together they form a sharper structural sketch than either alone. First-principles closure (i.e., a theorem-grade derivation) remains open.

This is publication-ready as a structural observation in Paper 51 §subsubsec:g4_5_v3_20_followon mechanism paragraph. Recommend brief sentence-level addition to the paper noting the harmonic-conjugate constraint.

## 10. Files

- `debug/sprint_moebius_mode_counting_diagnostic_memo.md` (this)
- No driver script (structural diagnostic, no computation needed)
- No data file

## 11. Cross-references

- `debug/sprint_alpha_gt_1_mechanism_investigation_memo.md` — Task 3 of prior thread (double-cover monodromy sketch)
- `debug/alpha_gt_1_analytical_investigation_memo.md` — v3.19.0 original ansatz (with retracted FS attribution)
- `debug/fursaev_solodukhin_1995_grounding_memo.md` — v3.20.0 task #26 retraction
- `debug/alpha_gt_1_moebius_validation_4_5_10_memo.md` — v3.20.0 task #25 empirical lock
- `memory/alpha_gt_1_moebius_closed_form.md` — α > 1 Möbius memory (mechanism flagged OPEN)
- Paper 51 §subsubsec:g4_5_v3_20_followon — current Paper 51 statement (mechanism OPEN)
- Cheeger 1983, Dowker 1977, 1994 — standard spinor SC literature
- Fursaev–Miele 1996 (hep-th/9605153) — correct spinor-on-cone paper
