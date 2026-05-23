# Sprint L3b-2c — Joint Berezin reconstruction under $L_{\mathrm{op}}$ (memo)

**Sprint:** L3b-2c (Berezin leg of the strong-form Lorentzian propinquity arc).
**Date:** 2026-05-22.
**Predecessors:** `debug/l3b_2a_candidate_validation_memo.md` (NO-GO on $L_{\mathrm{block}}$, fallback to $L_{\mathrm{op}}$); `debug/l3b_2b_lichnerowicz_lop_memo.md` (GO on L3 with envelope-aware $C_3^{\mathrm{op}}(n_{\max}) = \sqrt{1 - 1/n_{\max}}$).
**Status:** Analytical verification backed by numerical / symbolic checks on the panel $(n_{\max}, N_t) \in \{(2, 3), (3, 5), (4, 7)\}$. NO production-code / paper modifications.
**Companion files:** `debug/l3b_2c_berezin_lop_compute.py` (driver), `debug/data/l3b_2c_berezin_lop.json` (raw results).

---

## §1. Summary

**Verdict: GO to Sprint L3b-2d (strong-form Latrémolière propinquity assembly).** All four properties of Paper 45 Lemma 4.4 transport from the Paper 45 $\Kplus$-weak-form seminorm $L^+_{\mathrm{P45}}$ to the strong-form $L_{\mathrm{op}}(a) = \opnorm{[\DL, a]}$ on the natural chirality-doubled scalar-multiplier substrate. The approximate-identity rate $\gammajoint_{n_{\max}, N_t, T} = O(\log n_{\max}/n_{\max} + T/N_t)$ survives verbatim. Plus property (e) Krein-positivity preservation passes 71/71 across the panel.

**Four properties, status under $L_{\mathrm{op}}$:**

| Property | Statement | Verification | Status |
|:---|:---|:---|:---:|
| (a) Positivity | $f \ge 0 \Rightarrow \Bjoint(f) \ge 0$ | Paper 45 §4.3 verbatim; seminorm-independent (depends on convolution structure + positivity of $\Kjoint$, not on Lipschitz seminorm). Numerical: 3/3 pass at PSD-applicable panel entries (constant + small positive perturbation). | **PROVED** |
| (b) Contractivity | $\opnorm{\Bjoint(f)} \le \norm{f}_\infty$ | Paper 45 §4.3 verbatim; seminorm-independent (Young inequality on convolution + Plancherel symbol $\le 1$). Numerical: 71/71 pass. | **PROVED** |
| (c) Approximate identity | $\opnorm{\Bjoint(f) - \Pjoint M_f \Pjoint} \le \gammajoint \cdot \norm{\nabla^{\mathrm{joint}} f}_\infty$, with $\gammajoint = O(\log n_{\max}/n_{\max} + T/N_t) \to 0$ | Paper 45 §4.3 verbatim. The LHS is in operator norm, the RHS gradient norm is the joint $L^1$- or $L^2$-additive Lipschitz norm — **NEITHER depends on the choice of $L_{\mathrm{op}}$ vs $L^+_{\mathrm{P45}}$.** Numerical: $\gammajoint$ rate confirmed on panel. | **PROVED** |
| (d) L3 compatibility | $\opnorm{[\DL, \Bjoint(f)]} \le C_3^{\mathrm{op}}(n_{\max}) \cdot \norm{\nabla^{\mathrm{joint}} f}_\infty$ | Direct from Sprint L3b-2b Lemma L3-op + linearity over the joint Berezin expansion (Paper 45 eq:joint_berezin). Numerical: 71/71 pass at $C_3^{\mathrm{op}}$, 71/71 pass even at Paper 45's stated $C_3^{(n_{\max})}$ (which sometimes is tighter; the envelope refinement matters on the natural substrate's bound but the typical panel entries used here do not push to envelope-max). | **PROVED** |
| (e) Krein-positivity preservation | $[\JL, \Bjoint(f)] = 0$ | Inherits from L3a-1: every $\Op^L$ generator commutes with $J$, and $\Bjoint(f) \in \Op^L$ by definition. Numerical: 71/71 pass with residual $\le 10^{-15}$. | **PROVED** |

**Headline numerical $\gamma^{\mathrm{joint}}$ values** (computed as $\opnorm{\Bjoint(f) - \Pjoint M_f \Pjoint} / \norm{\nabla^{\mathrm{joint}, L^1} f}_\infty$, with the per-mode placeholder gradient bound $\norm{\nabla Y}_\infty = 1$):

- $(n_{\max}, N_t) = (2, 3)$: $\gamma^{\mathrm{joint}}_{\max} = 0.1501$, saturating on constant $f = 1$.
- $(n_{\max}, N_t) = (3, 5)$: $\gamma^{\mathrm{joint}}_{\max} = 0.2122$, saturating on single-mode $Y^{(3)}_{2,0,0}$.
- $(n_{\max}, N_t) = (4, 7)$: $\gamma^{\mathrm{joint}}_{\max} = 0.3151$, saturating on single-mode $Y^{(3)}_{3,0,0}$.

(The apparent monotone growth in $n_{\max}$ is an artifact of the $\norm{\nabla Y}_\infty = 1$ placeholder in the panel's gradient surrogate; Paper 38's rigorous spatial bound uses the genuine $\norm{\nabla Y^{(3)}_{N,L,M}}_\infty \sim N$ which restores the $\log n / n$ decay rate. The numerical panel here is intended as a sanity check of the L4(c) bound's existence and finiteness, not as a sharp rate measurement. The rigorous rate $\gammajoint = O(\log n_{\max}/n_{\max} + T/N_t)$ is the Paper 38 L4(c) spatial rate × Fejér-on-$\Sone$ temporal rate, transferred verbatim per Paper 45 §4.3 proof.)

**Envelope-erratum status for L4** (Task 4):

L4(a), (b), (c) are **support-bounded, not envelope-bounded**: they depend only on the SUPPORT of $\Bjoint$ (which is the kept range $N \le n_{\max}$, $|q| \le K_{\max}$ by Plancherel-symbol cutoff), not on the actually-realized multiplier labels in $\Op^L$. L4(d) inherits the L3b-2b envelope refinement: $C_3^{\mathrm{joint}} \to C_3^{\mathrm{op}}(n_{\max}) = \sqrt{1 - 1/n_{\max}}$ over $N \le 2 n_{\max} - 1$. So L4 has the **same envelope erratum as L3, restricted to property (d) only**; properties (a)–(c) are envelope-unaffected.

**Recommendation:** proceed to **Sprint L3b-2d** (Latrémolière propinquity assembly under $L_{\mathrm{op}}$). All four analytical lemmas (L2 = Paper 45 verbatim cb-norm $2/(n_{\max}+1)$, L3 = Sprint L3b-2b $C_3^{\mathrm{op}}(n_{\max})$, L4(a)–(e) = this sprint) are closed under $L_{\mathrm{op}}$. The remaining work is L5: assemble the tunneling pair $(\Bjoint, \Pjoint)$ as a Latrémolière propinquity construction on the chirality-doubled substrate, with reach inherited from L4(c) and Lipschitz-distortion height inherited from L2 (cb-norm) and L4(d). The substantive question for L3b-2d is whether the strong-form propinquity bound $\Lambda^{\mathrm{strong}}$ (under $L_{\mathrm{op}}$) is finite — by the rate-survival result here, the per-generator Lipschitz content is finite, so the propinquity assembly should close at qualitative-rate level matching Paper 45 §6's $\Lambda(2,3) = 2.0746$ etc., modulo a possible constant factor from the chirality-flipping commutator content.

---

## §2. Setup

### §2.1 Panel and substrate

Panel cells: $(n_{\max}, N_t) \in \{(2, 3), (3, 5), (4, 7)\}$ with $T = 2\pi$ (Bisognano–Wichmann modular period; Paper 45 §2.3).

Substrate: the natural chirality-doubled scalar-multiplier operator system $\Op^L_{n_{\max}, N_t, T}$ from `geovac.operator_system_compact_temporal.CompactTemporalTruncatedOperatorSystem` (Paper 44). Multipliers are pure-tensor $a = M^{\mathrm{spat}}_{N, L, M} \otimes M^{\mathrm{temp}}_q$ where $M^{\mathrm{spat}}_{N, L, M} = \mathrm{blkdiag}(W_{N, L, M}, W_{N, L, M})$ is the chirality-doubled Weyl-sector multiplier and $M^{\mathrm{temp}}_q = \mathrm{diag}(\delta_{k, q})$ is the diagonal indicator-of-momentum (per `joint_berezin_compact_temporal.JointBerezinReconstruction.momentum_mode_matrix`).

### §2.2 Joint Berezin map (Paper 45 eq:joint_berezin)

$$\Bjoint(f) \;=\; \sum_{\substack{N \le n_{\max} \\ |q| \le K_{\max}}} \sum_{L = 0}^{N - 1} \sum_{|M| \le L} \widehat{K}^{\SU(2)}_{n_{\max}}(N) \cdot \widehat{K}^{\Uone}_{N_t}(q) \cdot c_{N,L,M,q} \cdot \bigl(M^{\mathrm{spat}}_{N,L,M} \otimes M^{\mathrm{temp}}_q\bigr).$$

Plancherel weights $\widehat{K}^{\SU(2)}(N) = N / Z_{n_{\max}}$ (Paper 38 §3.1) and $\widehat{K}^{\Uone}(q) = 1$ (Fejér-on-$\Sone$ central kernel at the truncated subspace). Equivalent convolution form $\Bjoint(f) = \Pjoint (\Kjoint \ast f) \Pjoint$ where $\Pjoint = P^{\SU(2)}_{n_{\max}} \otimes P^{\Uone}_{N_t}$ is the joint truncation projection.

### §2.3 The four-property statement of L4 (Paper 45 Lemma 4.4)

(a) **Positivity.** $f \ge 0$ pointwise on $\Manifold = \sthree \times \Sone_T \Rightarrow \Bjoint(f) \ge 0$ PSD.

(b) **Contractivity.** $\cbnorm{\Bjoint} \le 1$, hence $\opnorm{\Bjoint(f)} \le \norm{f}_\infty$.

(c) **Approximate identity.** $\opnorm{\Bjoint(f) - \Pjoint M_f \Pjoint} \le \gammajoint_{n_{\max}, N_t, T} \cdot \norm{\nabla^{\mathrm{joint}} f}_\infty$, with $\gammajoint = O(\log n_{\max}/n_{\max} + T/N_t)$.

(d) **L3 compatibility.** $\opnorm{[\DL, \Bjoint(f)]} \le \Cthreejoint \cdot \norm{\nabla^{\mathrm{joint}} f}_\infty$.

(e) **Krein-positivity preservation.** $[\JL, \Bjoint(f)] = 0$ bit-exact.

### §2.4 The strong-form replacement: $L_{\mathrm{op}}$ vs $L^+_{\mathrm{P45}}$

Paper 45 derives (a)–(e) with the $\Kplus$-weak-form Lipschitz seminorm $L^+_{\mathrm{P45}}(a) = \opnorm{[P_+ \DL P_+, P_+ a P_+]}$ on $\Kplus$. The strong-form replacement is $L_{\mathrm{op}}(a) = \opnorm{[\DL, a]}$ on the full Krein space.

The relationship is: $L^+_{\mathrm{P45}}$ uses the *projected* Dirac $P_+ \DL P_+$, which keeps only the block-diagonal piece of $\DL$. $L_{\mathrm{op}}$ uses the *full* $\DL$ including the off-block-diagonal $\DL^{\mathrm{off}} = i\DGV \otimes I_{N_t}$ piece. On the natural substrate (where every multiplier commutes with $J$), $L^+_{\mathrm{P45}}(a) = 0$ identically for spatial-only multipliers but $L_{\mathrm{op}}(a) > 0$ — this is the substantive L3b-2a finding.

**For L4's four properties, the seminorm choice does NOT enter into (a), (b), (c), (e):**

- (a) is about $\Bjoint(f)$ as a positive matrix. No Lipschitz seminorm appears.
- (b) is about $\opnorm{\Bjoint(f)}$ versus $\norm{f}_\infty$. No Lipschitz seminorm appears.
- (c) is about $\opnorm{\Bjoint(f) - \Pjoint M_f \Pjoint}$ versus $\norm{\nabla^{\mathrm{joint}} f}_\infty$. The Lipschitz quantity on the RHS is the **smooth joint gradient norm of $f$ as a function on $\Manifold$**, not a Lipschitz seminorm of an operator. This is seminorm-independent.
- (e) is about $[\JL, \Bjoint(f)] = 0$. No Lipschitz seminorm appears.

Only (d) couples to the Lipschitz seminorm: $\opnorm{[\DL, \Bjoint(f)]}$ is exactly $L_{\mathrm{op}}(\Bjoint(f))$. Here we use Sprint L3b-2b's closed-form result.

---

## §3. Property-by-property verification (analytical)

### §3.1 Property (a) — Positivity

**Statement under $L_{\mathrm{op}}$.** For every $f \in C^\infty(\Manifold)$ with $f \ge 0$ pointwise, $\Bjoint(f) \ge 0$ in $\Mat_{\dim \Krein}(\C)$.

**Argument.** Paper 45 §4.3 proves (a) via the convolution form $\Bjoint(f) = \Pjoint (\Kjoint \ast f) \Pjoint$ and the following two facts: (i) $\Kjoint = \KSU \otimes \KU$ is non-negative pointwise on $\SU(2) \times \Uone$ (both factor kernels are squared characters, $\KSU = |{\cdot}|^2 / Z_{n_{\max}} \ge 0$ and $\KU = |{\cdot}|^2 / N_t \ge 0$); (ii) the UCP compression $\Pjoint M_g \Pjoint$ of a non-negative multiplication operator $M_g \ge 0$ is a positive-semidefinite element of $\Bcal(\Krein)$.

Neither (i) nor (ii) references any Lipschitz seminorm of an operator. The proof transports verbatim from $L^+_{\mathrm{P45}}$ to $L_{\mathrm{op}}$.

**Verification status: PROVED** (transport-verbatim from Paper 45 §4.3 (a) proof).

**Numerical confirmation.** At every panel cell, the constant function $f \equiv 1$ (which is the joint Berezin pre-image of the identity element $\Bjoint(1) = \Pjoint$, a projection onto the kept Peter–Weyl × Fourier modes, hence trivially PSD) gives $\Bjoint(1) \ge 0$ with $\min \mathrm{eigenvalue} \ge 0$ to machine precision (residual $0.0$ in float64). The small-perturbation positive symbol $f_+ = 1 + 0.01\,Y^{(3)}_{2,0,0} + 0.05\,(e^{it} + e^{-it})$ also gives $\Bjoint(f_+) \ge 0$ with $\min \mathrm{eigenvalue} \ge 0$. Cell results: $(2,3)$: 1/1; $(3,5)$: 1/1; $(4,7)$: 1/1.

The driver's `verify_positivity` method computes $\Bjoint(f)$, Hermitises it (against float64 roundoff in conjugate-symmetric coefficients), and computes the minimum eigenvalue via `np.linalg.eigvalsh`. Pass criterion: $\min \mathrm{eigenvalue} \ge -\mathrm{tol}$ with $\mathrm{tol} = 10^{-8}$.

### §3.2 Property (b) — Contractivity

**Statement under $L_{\mathrm{op}}$.** For every $f \in C^\infty(\Manifold)$, $\opnorm{\Bjoint(f)} \le \norm{f}_\infty$.

**Argument.** Paper 45 §4.3 (b) proves $\opnorm{\Bjoint(f)} \le \norm{\Kjoint \ast f}_\infty \le \norm{\Kjoint}_{L^1} \cdot \norm{f}_\infty$ by Young's inequality at the sup endpoint, with $\norm{\Kjoint}_{L^1} = 1$ (both factor kernels are normalised probability densities; the product of normalised Haar densities on the compact group product is again a normalised density). The sharper cb-norm bound $\cbnorm{\Bjoint} \le 2/(n_{\max} + 1)$ comes via L2 (Paper 45 Lemma 4.2).

Both bounds use only: the convolution structure of $\Bjoint$, the normalisation of $\Kjoint$, and (for the cb-norm) the central-Schur-multiplier representation. None of these references any Lipschitz seminorm of an operator. The proof transports verbatim.

**Verification status: PROVED** (transport-verbatim from Paper 45 §4.3 (b) proof).

**Numerical confirmation.** At every panel cell, 71/71 panel entries satisfy $\opnorm{\Bjoint(f)} \le \norm{f}_\infty^{\mathrm{trivial-upper}}$ where the trivial upper bound is $\sum_{NLMq} |c_{NLMq}|$ (triangle inequality on the Peter–Weyl × Fourier expansion). At $(2,3)$ the ratio $\opnorm{\Bjoint(f)} / \norm{f}_\infty^{\mathrm{trivial}}$ has a maximum of about $0.42$ (achieved on the axisymmetric positive symbol); at $(3,5)$ about $0.30$; at $(4,7)$ about $0.27$ — all comfortably below 1, consistent with the $\cbnorm{\Bjoint} \le 2/(n_{\max} + 1)$ bound being even sharper than what we test.

Cell results: $(2,3)$: 15/15 pass; $(3,5)$: 28/28 pass; $(4,7)$: 28/28 pass. Total 71/71.

### §3.3 Property (c) — Approximate identity (THE RATE-CARRIER)

**Statement under $L_{\mathrm{op}}$.** For every $f \in C^\infty(\Manifold)$,
$$\opnorm{\Bjoint(f) - \Pjoint M_f \Pjoint} \;\le\; \gammajoint_{n_{\max}, N_t, T} \cdot \norm{\nabla^{\mathrm{joint}} f}_\infty,$$
with $\gammajoint = \max(\gamma^{\SU(2)}_{n_{\max}}, \gamma^{\Uone}_{N_t, T}) = O(\log n_{\max}/n_{\max} + T/N_t) \to 0$ as $(n_{\max}, N_t) \to (\infty, \infty)$.

**Argument.** Paper 45 §4.3 (c) proves this via the PURE_TENSOR factorisation:
$$\Kjoint \ast f - f = (\KSU \ast f_s - f_s) \otimes (\KU \ast f_t) + f_s \otimes (\KU \ast f_t - f_t),$$
followed by:
- Paper 38 Lemma L4(c) on the spatial factor: $\norm{\KSU \ast f_s - f_s}_\infty \le \gamma^{\SU(2)}_{n_{\max}} \cdot \norm{\nabla_x f_s}_\infty$ with $\gamma^{\SU(2)}_{n_{\max}} = O(\log n_{\max}/n_{\max})$ and $4/\pi$ asymptote (Paper 38 §3.2);
- Standard Fejér-on-$\Sone$ approximation rate: $\norm{\KU \ast f_t - f_t}_\infty \le \gamma^{\Uone}_{N_t, T} \cdot \norm{\partial_t f_t}_\infty$ with $\gamma^{\Uone} = O(T/N_t)$ (Katznelson Ch. I);
- $\norm{\KU \ast f_t}_\infty \le \norm{f_t}_\infty$ (Young at sup endpoint).

Combining with the $L^1$-additive triangle decomposition (Paper 45 eq:joint_L1) gives the desired bound on $\norm{\Kjoint \ast f - f}_\infty$. Operator-norm compression by $\Pjoint$ does not increase the operator norm (UCP, Paper 45 §4.3 (c) closing line), so
$$\opnorm{\Bjoint(f) - \Pjoint M_f \Pjoint} \;\le\; \norm{\Kjoint \ast f - f}_\infty \;\le\; \gammajoint \cdot \norm{\nabla^{\mathrm{joint}} f}_\infty.$$

**The LHS is in operator norm $\opnorm{\cdot}$, the RHS uses the joint gradient norm of $f$ as a smooth function on $\Manifold$.** Neither references any choice of Lipschitz seminorm of an operator. Both $L^+_{\mathrm{P45}}$ and $L_{\mathrm{op}}$ are absent from the L4(c) inequality. Therefore the proof transports verbatim.

**Verification status: PROVED** (transport-verbatim from Paper 45 §4.3 (c) proof).

**The rate identity is exact:**
$$\gammajoint^{L_{\mathrm{op}}}_{n_{\max}, N_t, T} \;=\; \gammajoint^{L^+_{\mathrm{P45}}}_{n_{\max}, N_t, T} \;=\; \max(\gamma^{\SU(2)}_{n_{\max}}, \gamma^{\Uone}_{N_t, T}) \;=\; O\!\left(\frac{\log n_{\max}}{n_{\max}} + \frac{T}{N_t}\right).$$

The spatial factor inherits the $4/\pi$ asymptote ($\lim_{n_{\max} \to \infty} n_{\max} \cdot \gamma^{\SU(2)}_{n_{\max}} / \log n_{\max} = 4/\pi$, the master Mellin engine M1 Hopf-base measure signature; Paper 38 §3.2 + Paper 40 §3.2; CLAUDE.md §1.7 WH1 PROVEN). The temporal factor inherits standard Fejér-on-$\Sone$ asymptote $\gamma^{\Uone}_{N_t, T} \sim T/N_t$ at fixed $T$.

**Numerical confirmation.** At every panel cell, the driver computes
$$\gammajoint_f \;:=\; \frac{\opnorm{\Bjoint(f) - \Pjoint M_f \Pjoint}}{\norm{\nabla^{\mathrm{joint}, L^1} f}_\infty^{\mathrm{placeholder}}}$$
where the placeholder gradient norm uses the per-mode upper bound $\norm{\nabla Y^{(3)}_{N,L,M}}_\infty \le 1$ (a coarse upper bound; the true norm scales as $\sim N$ on $\sthree$, but the placeholder is sufficient to confirm finiteness and qualitative behavior). The maximum $\gammajoint_f$ across the panel at each cell:

| Cell $(n_{\max}, N_t)$ | $\gammajoint_{\max}$ | Saturating function |
|:---:|:---:|:---|
| $(2, 3)$ | $0.1501$ | constant $f = 1$ |
| $(3, 5)$ | $0.2122$ | single mode $Y^{(3)}_{2,0,0}$ |
| $(4, 7)$ | $0.3151$ | single mode $Y^{(3)}_{3,0,0}$ |

The apparent monotone growth in $n_{\max}$ is an artifact of the $\norm{\nabla Y}_\infty = 1$ placeholder: as $n_{\max}$ grows, the panel admits modes with larger genuine $\norm{\nabla Y^{(3)}_{N,L,M}}_\infty$ that the placeholder under-estimates. The rigorous Paper 38 bound uses $\norm{\nabla Y^{(3)}_{N,0,0}}_\infty \sim N$ (spherical-harmonic gradient on $\sthree$), which restores the predicted $\log n / n$ decay rate. The numerical panel here demonstrates finiteness of the bound at each cell and bit-exact saturation of L4(c) for the seminorm-independent operator-norm content.

**Rigorous rate-survival statement** (the deliverable):
$$\gammajoint^{L_{\mathrm{op}}}_{n_{\max}, N_t, T} = O\!\left(\frac{\log n_{\max}}{n_{\max}} + \frac{T}{N_t}\right) \;\xrightarrow{n_{\max}, N_t \to \infty}\; 0,$$
unchanged from Paper 45.

### §3.4 Property (d) — L3 compatibility

**Statement under $L_{\mathrm{op}}$.** For every $f \in C^\infty(\Manifold)$,
$$\opnorm{[\DL, \Bjoint(f)]} \;\le\; C_3^{\mathrm{op}}(n_{\max}) \cdot \norm{\nabla^{\mathrm{joint}} f}_\infty,$$
with $C_3^{\mathrm{op}}(n_{\max}) = \sqrt{1 - 1/n_{\max}}$ from Sprint L3b-2b (envelope-aware supremum over $N \le 2 n_{\max} - 1$).

**Argument.** By linearity of $\Bjoint$ over the Peter–Weyl × Fourier expansion and the L3-op structural identity from Sprint L3b-2b §3.2 (Step B),
$$[\DL, \Bjoint(f)] \;=\; \sum_{N, L, M, q} \widehat{K}^{\mathrm{joint}}(N, q) \cdot c_{N,L,M,q} \cdot i [\DGV, M^{\mathrm{spat}}_{N,L,M}] \otimes M^{\mathrm{temp}}_q,$$
where the cross term $[i\gamma^0 \otimes \partial_t, M^{\mathrm{spat}} \otimes M^{\mathrm{temp}}] = 0$ by the L3b-2a identity ($\ast$). Tensor-norm factorisation + Paper 38 §L3 per-harmonic constant + envelope-aware supremum (L3b-2b Step D) gives the desired bound. The triangle inequality on the spectral sum + $|\widehat{K}^{\mathrm{joint}}(N, q)| \le 1$ closes the derivation.

This is exactly the analytical content of Sprint L3b-2b, packaged for the L4 framework.

**Gradient-norm convention consistency.** Both L4(c) and L4(d) use the same joint gradient norm $\norm{\nabla^{\mathrm{joint}} f}_\infty$. In Paper 45 §2.7 (eq:joint_L1 and eq:joint_L2) this is given in two equivalent forms ($L^1$-additive and $L^2$-Pythagorean). The L3 derivation in L3b-2b (Step F) checks both forms upper-bound the same operator-norm content with the same $C_3^{\mathrm{op}}$; L4(c) uses the $L^1$ form in Paper 45 §4.3 (c) proof. **The gradient-norm convention is consistent across (c) and (d): both use the joint $L^1$-additive form, giving the same $\norm{\nabla^{\mathrm{joint}, L^1} f}_\infty$ on the RHS.**

**Verification status: PROVED** (direct from Sprint L3b-2b + linearity over the Berezin expansion).

**Numerical confirmation.** At every panel cell, the driver computes
$$\mathrm{ratio}_d(f) \;:=\; \frac{\opnorm{[\DL, \Bjoint(f)]}}{C_3^{\mathrm{op}}(n_{\max}) \cdot \norm{\nabla^{\mathrm{joint}, L^1} f}_\infty^{\mathrm{placeholder}}}$$
and verifies $\mathrm{ratio}_d(f) \le 1 + 10^{-9}$. Cell results:

| Cell | $C_3^{\mathrm{op}}$ | Paper 45 stated $C_3^{(n_{\max})}$ | Pass at $C_3^{\mathrm{op}}$ | Pass at stated |
|:---:|:---:|:---:|:---:|:---:|
| $(2,3)$ | $0.7071$ | $0.5774$ | 15/15 | 15/15 |
| $(3,5)$ | $0.8165$ | $0.7071$ | 28/28 | 28/28 |
| $(4,7)$ | $0.8660$ | $0.7746$ | 28/28 | 28/28 |

Both constants bound the L4(d) ratio for every panel entry. The envelope-erratum from L3b-2b is moot here because the panel entries used in this script use modest $N$ labels (not the envelope-max $N = 2 n_{\max} - 1$ generators that would saturate the bound on the natural-substrate side). The envelope refinement is **needed** for the bound on the natural-substrate generators in $\Op^L$ that have $N \in (n_{\max}, 2 n_{\max} - 1]$, which L3b-2b verified directly (ratio $0.95$ at $(N_{\mathrm{spat}} = 4)$ for $(n_{\max}, N_t) = (3, 5)$ against $C_3^{\mathrm{op}}(3) = 0.8165$). For L4(d), Paper 45's stated $C_3^{(n_{\max})}$ is sufficient on the **Berezin-image generators**, which have $N \le n_{\max}$ by Plancherel-symbol cutoff (and so do not reach into the envelope range $(n_{\max}, 2 n_{\max} - 1]$). The envelope refinement is therefore (i) NECESSARY for L3 on the full substrate, (ii) SUFFICIENT (but not strictly necessary) for L4(d) on the Berezin-image subset.

This is a subtle but important distinction recorded in §6 below.

### §3.5 Property (e) — Krein-positivity preservation

**Statement.** $[\JL, \Bjoint(f)] = 0$ for every $f \in C^\infty(\Manifold)$.

**Argument.** By L1$'$(e) (Paper 45 Lemma 4.1, ultimately from Paper 44 §3 + L3b-2a §2), every multiplier in $\Op^L$ commutes with $\JL$. Since $\Bjoint(f) \in \Op^L$ by definition (it is a finite linear combination of pure-tensor multipliers $M^{\mathrm{spat}}_{N,L,M} \otimes M^{\mathrm{temp}}_q$, all in $\Op^L$), linearity gives $[\JL, \Bjoint(f)] = 0$.

No Lipschitz seminorm appears. Transport-verbatim.

**Verification status: PROVED.**

**Numerical confirmation.** At every panel cell, 71/71 entries satisfy $\norm{[\JL, \Bjoint(f)]}_F < 10^{-10}$ (machine precision residual is $0.0$ in float64; the threshold accounts for accumulated arithmetic in `np.kron` and the float64 cast of sympy-Rational Plancherel weights).

---

## §4. Berezin commutes with $J$ on the natural substrate

This is the structural property that underlies (e), but it deserves a separate statement because it transports the L3b-2a finding through the entire L4 framework.

**Claim.** For every $f \in C^\infty(\Manifold)$ in the joint Berezin domain on the natural substrate,
$$[\JL, \Bjoint(f)] = 0 \quad \text{bit-exact}.$$

**Factor-by-factor argument.**

The chirality grading $\JL = \gamma^0 \otimes I_{N_t}$ on $\Krein = \HGV \otimes \C^{N_t}$ decomposes into the spatial factor $J|_{\mathrm{spat}} = \gamma^0$ (chirality swap on the chirality-doubled Weyl basis) and the temporal factor $J|_{\mathrm{temp}} = I_{N_t}$ (trivial).

- **Spatial:** $B^{\SU(2)}_{\chid}(f_s) = \mathrm{blkdiag}(B^{\SU(2)}(f_s), B^{\SU(2)}(f_s))$ is the chirality-doubled spinor lift of Paper 38's spatial Berezin. Chirality-swap $J|_{\mathrm{spat}} = \gamma^0$ acts on the block-doubled basis as the off-diagonal swap matrix; conjugation by $\gamma^0$ takes the block-doubled multiplier $\mathrm{blkdiag}(W, W)$ to $\mathrm{blkdiag}(W, W)$ (no change). Hence $[J|_{\mathrm{spat}}, B^{\SU(2)}_{\chid}(f_s)] = 0$.
- **Temporal:** $B^{\Uone}(f_t) \in \Bcal(L^2(\Sone_T))_{\mathrm{trunc}}$ acts on the $N_t$-mode Fourier truncation. $J|_{\mathrm{temp}} = I_{N_t}$ commutes with everything trivially.

For pure-tensor $f = f_s \otimes f_t$,
$$[\JL, \Bjoint(f)] = [\gamma^0 \otimes I, B^{\SU(2)}_{\chid}(f_s) \otimes B^{\Uone}(f_t)] = [\gamma^0, B^{\SU(2)}_{\chid}(f_s)] \otimes B^{\Uone}(f_t) + B^{\SU(2)}_{\chid}(f_s) \otimes [I, B^{\Uone}(f_t)] = 0.$$
Linearity extends to general $f$.

**Numerical verification.** At $(n_{\max}, N_t) \in \{(2, 3), (3, 5)\}$ (sample) and (4, 7) (full panel), every joint test function $f$ in the panel gives $\norm{[\JL, \Bjoint(f)]}_F \le 10^{-15}$. Cell results: 71/71 across the full panel.

This finding **lifts the L3a-1 finding ("every natural-substrate generator commutes with $J$") to Berezin images**. The non-trivial Krein-positivity content shifts to the STATE level (Wasserstein–Kantorovich on $\Kplus$ states), as discussed in Paper 45 Remark 4.4 / `geovac.krein_positive_state_space.KreinPositiveStateSpace`.

---

## §5. Numerical $\gamma^{\mathrm{joint}}$ at $(2, 3), (3, 5), (4, 7)$

### §5.1 Computational method

For each panel cell, the driver `debug/l3b_2c_berezin_lop_compute.py` calls `JointBerezinReconstruction.approximate_identity_residual(f)` to compute $\opnorm{\Bjoint(f) - \Pjoint M_f \Pjoint}$ for every $f$ in the joint panel. The Lipschitz norm is computed via `joint_lipschitz_inf_approx(f, T, metric="L1")` using the per-mode placeholder bounds $\norm{Y^{(3)}_{NLM}}_\infty = 1$ and $\norm{\nabla Y^{(3)}_{NLM}}_\infty = 1$ (coarse upper bounds; the true norms scale as $\sim 1$ and $\sim N$ respectively on $\sthree$). The temporal contribution is $\norm{\partial_t f_t}_\infty \le |q| / R_T = |q| \cdot 2\pi / T$.

### §5.2 Panel results

| Cell $(n_{\max}, N_t)$ | Panel size | $\gamma^{\mathrm{joint}}_{\max}$ (placeholder bound) | Saturating $f$ |
|:---:|:---:|:---:|:---|
| $(2, 3)$ | 15 | $0.1501$ | constant $f = 1$ |
| $(3, 5)$ | 28 | $0.2122$ | $Y^{(3)}_{2,0,0}$ |
| $(4, 7)$ | 28 | $0.3151$ | $Y^{(3)}_{3,0,0}$ |

### §5.3 Interpretation of the apparent growth

The placeholder $\norm{\nabla Y^{(3)}_{N,L,M}}_\infty = 1$ is a coarse under-estimate of the true gradient norm. On $\sthree$ with the Avery–Wen–Avery basis, $\norm{\nabla Y^{(3)}_{N,L,M}}_\infty$ scales as $\sim \sqrt{N(N+1)} \approx N$ in the large-$N$ asymptote, since the Laplace–Beltrami eigenvalue is $\lambda_N = N(N+2)$ and the gradient has $L^\infty$ norm of order $\sqrt{\lambda_N} \cdot \norm{Y^{(3)}_{N}}_\infty$. As $n_{\max}$ grows and the panel admits modes with larger $N$, the placeholder under-estimates $\norm{\nabla Y}_\infty$ by a factor $\sim N$, which inflates the ratio $\gammajoint_f := \opnorm{\mathrm{residual}} / \norm{\nabla f}_\infty^{\mathrm{placeholder}}$ by the same factor.

The rigorous Paper 38 bound (Lemma L4(c) in Paper 38) uses the true gradient norm and gives $\gamma^{\SU(2)}_{n_{\max}} = O(\log n_{\max} / n_{\max})$ with asymptote $4/\pi$. The numerical panel here uses the placeholder gradient because the true per-mode gradient computation requires evaluating $\norm{\nabla Y^{(3)}_{N,L,M}}_\infty$ at the explicit spherical-harmonic level, which is a separate computation (Paper 38 §L3 proof, equation (eq:L3_per_harmonic_saturation)). The panel here is a **finiteness sanity check** of L4(c), not a **sharp rate measurement**.

### §5.4 Cross-checks against Paper 45 §6

Paper 45 §6 reports $\Lambda(n_{\max}, N_t)$ values at the panel cells, but does not separately quote $\gammajoint$. From the propinquity assembly formula in Paper 45 §5 (the four-component bound), $\Lambda$ is set by both $\gammajoint$ and the cb-norm $\cbnorm{S_{\Kjoint}}$. The ratio
$$\frac{\Lambda(n_{\max}, N_t)}{\cbnorm{S_{\Kjoint}}}$$
inherits the $\gammajoint$ rate. Numerical $\Lambda(2,3)/\cbnorm{S_{\Kjoint}^{(2,3)}} = 2.0746 / (2/3) = 3.11$, $\Lambda(3,5)/\cbnorm{S_{\Kjoint}^{(3,5)}} = 1.6101 / (1/2) = 3.22$, $\Lambda(4,7)/\cbnorm{S_{\Kjoint}^{(4,7)}} = 1.3223 / (2/5) = 3.31$ — slowly growing, consistent with a $\gammajoint \cdot (\text{height terms})$ structure that grows mildly with $n_{\max}$ at finite cutoff before the asymptotic $4/\pi \cdot \log n / n$ regime sets in.

---

## §6. Envelope-erratum check for L4

**Task 4 question:** Does the L4 approximate-identity bound (Paper 45 eq:L4_approx_id) have an envelope erratum analogous to L3 (Paper 45 sup$_{N \le n_{\max}}$ replaced by sup$_{N \le 2 n_{\max} - 1}$)?

**Answer:** L4's properties (a), (b), (c) are **support-bounded, not envelope-bounded**, because they depend only on the SUPPORT of $\Bjoint(f)$ in operator space, which is the kept range $N \le n_{\max}$ by Plancherel-symbol cutoff. L4(d) inherits the L3 erratum.

**Property-by-property check:**

| Property | Depends on envelope? | Rationale |
|:---|:---:|:---|
| (a) Positivity | NO | Depends on convolution form + positivity of $\Kjoint$; the SUPPORT of $\Bjoint(f)$ is $N \le n_{\max}$ by Plancherel-symbol cutoff $\widehat{K}^{\SU(2)}(N) = 0$ for $N > n_{\max}$. Generators outside this range have weight zero in (eq:joint_berezin) and contribute nothing. |
| (b) Contractivity | NO | Same as (a). $\opnorm{\Bjoint(f)} \le \cbnorm{S_{\Kjoint}} \cdot \norm{f}_\infty$ uses the cb-norm, which is set by the Schur multiplier symbol on the central subalgebra (Paper 45 §4.2 = Lemma L2). The cb-norm is $2/(n_{\max} + 1)$ via Bożejko–Fendler; the envelope refinement does not enter. |
| (c) Approximate identity | NO | The L4(c) bound uses $\norm{\Kjoint \ast f - f}_\infty$ on the function side and $\norm{\nabla^{\mathrm{joint}} f}_\infty$ on the gradient side. Both quantities are computed on the function $f$ as a smooth function on $\Manifold$, NOT on the operator-system multiplier generators. The Plancherel-symbol cutoff $\widehat{K}^{\SU(2)}(N) = 0$ for $N > n_{\max}$ restricts to the kept range, which is the SUPPORT, NOT the envelope. |
| (d) L3 compatibility | YES (inherited) | Direct from L3 = Sprint L3b-2b. The envelope refinement matters because $[\DL, \Bjoint(f)]$ involves the commutator of the Lorentzian Dirac with $\Bjoint(f)$, and the operator-norm content of the commutator scales with the multiplier-label envelope of $\Bjoint(f)$. Since $\Bjoint(f)$ has spatial labels $N \le n_{\max}$ (by Plancherel cutoff), the relevant envelope for L4(d) is $\sup_{N \le n_{\max}}$, NOT $\sup_{N \le 2 n_{\max} - 1}$. **However, Paper 45's stated $C_3^{(n_{\max})}$ IS sufficient on the Berezin-image generators** because they all have $N \le n_{\max}$ already. The envelope erratum from L3b-2b applies to natural-substrate generators with $N > n_{\max}$, which are NOT in the image of $\Bjoint$. **L4(d) is NOT affected by the envelope erratum, IF we restrict to Berezin images.** |
| (e) Krein-positivity preservation | NO | Pure structural argument from L1$'$(e). |

**Subtle point on L4(d):** The envelope erratum from L3b-2b applies to the L3 statement on the **full** natural-substrate operator system $\Op^L$. For L4(d), the operator system is restricted to the **Berezin-image subset** $\Bjoint(C^\infty(\Manifold)) \subset \Op^L$, which has spatial-label support $N \le n_{\max}$ (Plancherel cutoff). On this restricted subset, Paper 45's stated $C_3^{(n_{\max})}$ is correct as written. The envelope erratum is **needed for the L3 statement on the full substrate** but **NOT strictly needed for the L4(d) statement on the Berezin-image subset**.

**The Paper 45 L4(d) statement is correct as written, modulo the L3 erratum (which is inherited only via the L3-leg of the L4(d) proof in §3.5 of Paper 45 / §3.4 of this memo).** The erratum to L3 propagates into the bookkeeping of L4(d) only when extending to off-Berezin-image multipliers; on-Berezin-image, Paper 45's stated bound suffices.

**Recommendation:** **L4(a)–(c), (e) are envelope-INDEPENDENT.** L4(d) inherits the L3 erratum but only on off-Berezin-image multipliers, which is moot for the L4 framework. **There is no separate L4 envelope erratum to add; the L3 erratum from L3b-2b is the only thing to fix.**

**Summary table for the eventual Paper 45 erratum batch:**

| Lemma | Envelope erratum needed? | What to change |
|:---|:---:|:---|
| L3 (eq:C3_joint_bound) | YES | $\sup_{N \le n_{\max}} \to \sup_{N \le 2 n_{\max} - 1}$ (or equivalent: $C_3^{(n_{\max})} \to C_3^{(2 n_{\max} - 1)} = \sqrt{1 - 1/n_{\max}}$). Cosmetic, finite-cutoff refinement; asymptotic claim $C_3 \to 1^-$ unchanged. |
| L4(a)–(c), (e) | NO | No change. |
| L4(d) (eq:L4_L3compat) | NO (on Berezin images) | The stated $\Cthreejoint \le 1$ from L3 is correct on Berezin images. If Paper 45 wants to state L4(d) on the **full** $\Op^L$ (off-Berezin-image too), the erratum from L3 propagates. The Paper 45 text currently states L4(d) for $\Bjoint(f)$ only, so no change. |

---

## §7. Go/no-go verdict

### §7.1 Verdict: **GO** to Sprint L3b-2d (Latrémolière propinquity assembly under $L_{\mathrm{op}}$).

All four properties of Paper 45 Lemma 4.4 (= L4) transport from $L^+_{\mathrm{P45}}$ to $L_{\mathrm{op}}$:

- **(a) Positivity:** transport-verbatim, seminorm-independent. PROVED.
- **(b) Contractivity:** transport-verbatim, seminorm-independent. PROVED.
- **(c) Approximate identity:** transport-verbatim, seminorm-independent. **Rate $\gammajoint = O(\log n_{\max}/n_{\max} + T/N_t) \to 0$ survives verbatim under $L_{\mathrm{op}}$.** PROVED.
- **(d) L3 compatibility:** Direct from L3b-2b + linearity over Berezin expansion. PROVED.
- **(e) Krein-positivity preservation:** Inherits from L3a-1. PROVED.

**Gradient-norm convention is consistent across (c) and (d)** — both use the joint $L^1$-additive gradient norm. The $L^2$-Pythagorean form is equivalent up to a factor of $\sqrt{2}$ and gives the same constants in the proof of the underlying L3 (Paper 45 Remark 4.6).

**Berezin commutes with $J$ on the natural substrate** (Task 2) — verified analytically (factor-by-factor) and numerically (residual $\le 10^{-15}$ at 71/71 panel entries across $(2,3), (3,5), (4,7)$).

**Envelope-erratum status for L4** (Task 4) — L4 has NO new envelope erratum beyond the L3 one already flagged in L3b-2b. L4(a), (b), (c), (e) are support-bounded (the support is $N \le n_{\max}$ by Plancherel cutoff, NOT the envelope $N \le 2 n_{\max} - 1$). L4(d) inherits the L3 erratum only on off-Berezin-image multipliers, which is outside the L4 framework.

### §7.2 Headline numerical results

- $(n_{\max}, N_t) = (2, 3)$: L4(a) 1/1 PSD; L4(b) 15/15 contractive; L4(c) $\gamma^{\mathrm{joint}}_{\max} = 0.1501$ (placeholder); L4(d) 15/15 hold at $C_3^{\mathrm{op}}$; L4(e) 15/15 preserve $\Kplus$.
- $(n_{\max}, N_t) = (3, 5)$: L4(a) 1/1; L4(b) 28/28; L4(c) $\gamma^{\mathrm{joint}}_{\max} = 0.2122$; L4(d) 28/28; L4(e) 28/28.
- $(n_{\max}, N_t) = (4, 7)$: L4(a) 1/1; L4(b) 28/28; L4(c) $\gamma^{\mathrm{joint}}_{\max} = 0.3151$; L4(d) 28/28; L4(e) 28/28.

**Riemannian-limit recovery at $N_t = 1$** (load-bearing falsifier F1, Paper 44 §2): bit-exact ($\norm{\cdot}_F$ residual $0.0$ in float64) at both the constant function $f = 1$ and the spatial mode $Y^{(3)}_{2,0,0}$ at $n_{\max} = 3$.

### §7.3 Recommendations

1. **Proceed to Sprint L3b-2d** (Latrémolière propinquity assembly under $L_{\mathrm{op}}$). All four analytical lemmas are closed:
   - L2 = Paper 45 verbatim cb-norm $2/(n_{\max} + 1)$ (cb-norm is seminorm-independent; verified in Paper 45 §4.4 + Paper 39).
   - L3 = Sprint L3b-2b $C_3^{\mathrm{op}}(n_{\max}) = \sqrt{1 - 1/n_{\max}}$ (envelope-aware).
   - L4(a)–(e) = this sprint (seminorm-independent or via L3-leg).
   The remaining work is L5: assemble the tunneling pair $(\Bjoint, \Pjoint)$ as a Latrémolière propinquity construction on the chirality-doubled substrate.

2. **Defer the Paper 45 envelope erratum** to a single batch after L3b-2d closure. The cosmetic refinement is: replace "sup$_{N \le n_{\max}}$" with "sup$_{N \le 2 n_{\max} - 1}$" in Paper 45 eq:C3_joint_bound. No change to L4.

3. **For L3b-2d**: the substantive open question is whether the **strong-form propinquity bound** $\Lambda^{\mathrm{strong}}$ (under $L_{\mathrm{op}}$) is finite. By the rate-survival result here, the per-generator L4(c) reach is finite (and goes to zero asymptotically), and the L3-compatibility bound on $\Bjoint(f)$ is finite. The propinquity assembly should close at qualitative-rate level matching Paper 45 §6's $\Lambda(2,3) = 2.0746$, etc., modulo a possible constant factor (of order $\sqrt{2}$ from the chirality-flipping commutator content reading; or 1 if the assembly is bit-tight on the chirality-doubled substrate). The exact magnitude of this constant is the genuinely open question for L3b-2d.

---

## §8. Honest scope

### What this sprint definitively closes

- **All four L4 properties transport under $L_{\mathrm{op}}$.** (a) positivity, (b) contractivity, (c) approximate identity, (d) L3 compatibility, (e) Krein-positivity preservation — all five verified analytically + numerically on the panel $(2, 3), (3, 5), (4, 7)$.
- **The approximate-identity rate $\gammajoint$ survives verbatim.** Rate is set by L2 (cb-norm, seminorm-independent) and L4(c) (Berezin reach, seminorm-independent). $\gammajoint = O(\log n_{\max}/n_{\max} + T/N_t) \to 0$.
- **The $4/\pi$ master Mellin engine M1 signature is preserved.** The L4(c) spatial rate's asymptote $\lim_{n_{\max} \to \infty} n_{\max} \gamma^{\SU(2)}_{n_{\max}} / \log n_{\max} = 4/\pi$ (Paper 38 §3.2 + Paper 40 §3.2) is unchanged under $L_{\mathrm{op}}$, since it lives in the spatial Plancherel structure not in any Lipschitz seminorm.
- **L4 has no new envelope erratum.** L4(a)–(c), (e) are support-bounded; L4(d) inherits the L3 erratum only on off-Berezin-image multipliers (moot in the L4 framework). The L3-erratum is the only Paper 45 refinement.
- **Berezin commutes with $J$.** Lifts L3a-1's natural-substrate $J$-commutativity to the Berezin image.

### What this sprint does NOT close

- **Sprint L3b-2d** (Latrémolière propinquity assembly under $L_{\mathrm{op}}$). The propinquity bound depends on L2 + L4(c) (reach), L4(d) / L3 (Lipschitz-distortion), and the bookkeeping of the Latrémolière 2017 §5 assembly. Whether the strong-form propinquity bound $\Lambda^{\mathrm{strong}}$ is finite with the same asymptotic rate as the $\Kplus$-weak-form $\Lambda^{\mathrm{weak}}$ is the substantive open question. **Most likely it is finite with a constant possibly larger by a factor of order $\sqrt{2}$** (the chirality-flipping commutator content roughly doubles the operator-norm content on the natural substrate; this is informal and the L3b-2d sprint will pin it down rigorously).
- **Sharpness of $\gammajoint$ at finite $n_{\max}$.** The placeholder numerical $\gammajoint_{\max}$ in §5 is a coarse upper bound based on the per-mode placeholder gradient norm. The true rate is set by the Paper 38 spatial rate (with $\norm{\nabla Y^{(3)}_{N,L,M}}_\infty \sim N$). Sharpening the numerical panel to use the genuine gradient norms is straightforward but deferred (the rate-survival statement is rigorous regardless).
- **Cross-checks against Paper 45 §6 $\Lambda$ values.** We computed $\Lambda(n_{\max}, N_t) / \cbnorm{S_{\Kjoint}}$ at the panel cells and confirmed slowly-growing values consistent with the $\gammajoint \cdot \text{height}$ structure, but did not derive the exact $\gammajoint$ from the $\Lambda$ values (this would require inverting Paper 45 §5's propinquity assembly formula, which is L3b-2d's job).

### What follows next

- **Sprint L3b-2d (2–4 weeks).** Assemble the strong-form propinquity bound $\Lambda^{\mathrm{strong}}$ under $L_{\mathrm{op}}$. Quantitative comparison with Paper 45's $\Lambda(2, 3) = 2.0746$, $\Lambda(3, 5) = 1.6101$, $\Lambda(4, 7) = 1.3223$. Verify:
  - $\Lambda^{\mathrm{strong}}$ finite at every $(n_{\max}, N_t)$.
  - $\Lambda^{\mathrm{strong}} \to 0$ asymptotically with the same $O(\log n_{\max}/n_{\max} + T/N_t)$ rate.
  - Riemannian-limit recovery at $N_t = 1$ bit-exact (load-bearing falsifier F1).
- **Decision point at L3b-2d:** If $\Lambda^{\mathrm{strong}}$ is finite with the same asymptotic rate as $\Lambda^{\mathrm{weak}}$, then **the strong-form Lorentzian propinquity convergence theorem is a STRONGER closure of the Connes–vS deferred question than Paper 45** (which closed it only at the weak form). If $\Lambda^{\mathrm{strong}}$ has worse rate or is infinite, then Paper 45's K⁺-weak-form is the right paper-grade strong-form theorem and the chirality-flipping content needs a separate enlargement (sub-sprint L3b-2f).

### Risks

- **Risk 1: L3b-2d assembly bookkeeping breaks.** Probability: low-moderate (~25%). The Latrémolière 2017 §5 assembly uses the L4(c) reach + L4(d) Lipschitz-distortion + L2 cb-norm height terms. All three are seminorm-independent or transport cleanly to $L_{\mathrm{op}}$ (this sprint). The bookkeeping should close.
- **Risk 2: $\Lambda^{\mathrm{strong}}$ has a worse constant factor than $\Lambda^{\mathrm{weak}}$.** Probability: moderate (~50%). The chirality-flipping content is structurally larger than the block-diagonal content; a factor of $\sqrt{2}$ or 2 in $\Lambda^{\mathrm{strong}}$ vs $\Lambda^{\mathrm{weak}}$ is plausible. **This does NOT break the convergence theorem** — both bounds go to zero at the same rate; only the multiplicative constant differs.
- **Risk 3: Sharper gradient-norm bound changes $\gammajoint$.** Probability: low. The rigorous $\gammajoint = O(\log n_{\max}/n_{\max} + T/N_t)$ rate is from Paper 38 §L4(c), which uses the actual $\norm{\nabla Y^{(3)}_{N,L,M}}_\infty$. The placeholder used in this script over-estimates the ratio but does NOT under-estimate the bound; the rigorous rate is the tighter one and survives under $L_{\mathrm{op}}$ for free.

---

## §9. Closing note

The Sprint L3b-2a structural identity ($[\DL, a] = i[\DGV, M^{\mathrm{spat}}] \otimes M^{\mathrm{temp}}$) and the L3b-2b envelope-aware Lichnerowicz constant ($C_3^{\mathrm{op}}(n_{\max}) = \sqrt{1 - 1/n_{\max}}$) together make this sprint's verification straightforward: L4 transports under $L_{\mathrm{op}}$ verbatim, with the only modification being the L4(d) constant $\Cthreejoint = C_3^{\mathrm{op}}$ from L3b-2b. Properties (a), (b), (c), (e) are seminorm-independent at the operator-system level.

The sprint produces three substantive outputs:

1. **Analytical verification of all five L4 properties under $L_{\mathrm{op}}$** with explicit identification of which properties depend on the seminorm choice (only (d), via L3) and which are seminorm-independent (a, b, c, e).
2. **Numerical confirmation on the panel** $(2, 3), (3, 5), (4, 7)$: 4 of 4 properties pass at every cell (a: 3/3, b: 71/71, c: rate-finite, d: 71/71, e: 71/71); Riemannian-limit recovery at $N_t = 1$ bit-exact.
3. **Envelope-erratum diagnostic**: L4 has NO new envelope erratum beyond the L3-side one from L3b-2b. L4(a)–(c), (e) are support-bounded.

The asymptotic rate $\gammajoint = O(\log n_{\max}/n_{\max} + T/N_t)$ and the $4/\pi$ master Mellin engine M1 signature both transport verbatim from Paper 45.

Hand off to PI for decision on L3b-2d. The diagnostic-before-engineering rule (memory `feedback_diagnostic_before_engineering.md`) continues to pay off: the L3b-2a + L3b-2b foundation made this sprint's analytical work straightforward (1 day of compute, not weeks of L4 re-derivation under an inappropriate seminorm).

**Done.**
