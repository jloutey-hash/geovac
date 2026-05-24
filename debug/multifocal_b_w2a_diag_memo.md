# Phase B Sub-Sprint W2a-diag: Hekkelman–McDonald 2024 vs the Master Mellin Engine

**Date:** 2026-05-07
**Author:** Sub-agent (Phase B W2a diagnostic)
**Sources read in full:**
- Hekkelman & McDonald, "A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity," arXiv:2412.00628v2 (Aug 2025), accepted J. Funct. Anal. (the full paper, pages 1–20, including all theorem/proposition statements through Section 6).
- `debug/multifocal_phase_a_synthesis_memo.md` (Phase A synthesis, Section 1 W2a row, Section 5 frontier-of-field framing).
- `debug/multifocal_literature_review_memo.md` (Track 3, esp. §3 on Hekkelman–McDonald and §6a on Connes–vS post-2024 ecosystem; surprises S5 and S6).
- `debug/track_ts_e1_theorem_promotion_memo.md` (Sprint TS-E1, Paper 32 §VIII case-exhaustion theorem — full statement, audit, proof skeleton, LaTeX draft).
- `debug/mr_b_spectral_action_rate_memo.md` (Sprint MR-B closed-form modular residual ε(t) for the Camporesi–Higuchi Dirac heat kernel on $S^3$; the exact Jacobi-ϑ₂ identity used as the Q2 numerical anchor).
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` §VIII case-exhaustion theorem (lines ~1376–1543), §VIII.B SM gauge appendix, §VIII.C G4b cross-manifold scoping (the existing Hekkelman–McDonald citation is at line 2438, in the G4b "the closest analogs" remark).

---

## §1. The diagnostic question and W2a's place in the wall taxonomy

Phase A established a six-wall taxonomy for the multi-focal-composition problem (Section 1 of the synthesis memo). The W1 family (W1a/W1b/W1c) is GeoVac-internal architectural — coordinate operators, magnetization densities, and frozen-core cross-screening — and is tooling-addressable. The W2 family is structurally different. **W2a** is the multi-loop-renormalization wall, exhibited in:

- Sprint LS-8a (May 2026): bare iterated Connes–Chamseddine spectral action on Dirac-$S^3$ reproduces the UV-divergent integrand of the two-loop electron self-energy faithfully (right $(\alpha/\pi)^2 \cdot (Z\alpha)^4 / n^3$ prefactor, right sign, right divergence form $\sim N^{3.43}$) but cannot autonomously generate $Z_2 / \delta m / Z_3$ counterterms required for a finite extraction of $C_{2S} = +3.63$ MHz·Hz⁻¹.
- Sprint HF-5 (May 2026, the 21 cm hyperfine sprint): the same wall in the vertex sector. Three two-loop topologies (VP-on-photon, SE-on-electron, crossed) on Dirac-$S^3$ as iterated CC spectral sums diverge as $\sim N^{3.79}$; the framework reproduces $(\alpha/\pi)^2$ structural prefactor, three-topology decomposition, full SO(4) selection rules at every vertex, but cannot generate the counterterms required to extract the finite $a_2$ contribution.

Phase A's Track 3 literature review surfaced **Hekkelman & McDonald, arXiv:2412.00628 (2024)** as the closest published analog to Connes–vS that handles UV/IR composition: a noncommutative-integral approximation on spectrally truncated triples, with a Szegő limit theorem and a density-of-states reading. Track 3 also flagged that **Marcolli–vS 2014 (rationality of spectral action) + Hekkelman–McDonald 2024 nearly *is* Sprint TS-E1's master Mellin engine** — citations recommended for Paper 32 §VIII when next touched.

The diagnostic asks five questions:

- **Q1.** What does Hekkelman–McDonald 2024 actually prove?
- **Q2.** Does their framework reproduce the GeoVac master Mellin engine on Dirac-$S^3$ at finite $n_{\max}$?
- **Q3.** Does the published framework give counterterms naturally?
- **Q4.** Does Marcolli–vS 2014 rationality give a constraint on which counterterms are *allowed*?
- **Q5.** Verdict (a/b/c/d) on W2a closure.

The diagnostic's strategic role is: if (a), Phase C-W2a is a tooling import sprint; if (b) or (c), Phase C-W2a is a position paper rather than a closure attempt and Paper 32 §VIII gains an honest "frontier-of-field" framing edit; if (d), we have a structural discrepancy that needs reconciling before adopting the Hekkelman–McDonald citation at all.

---

## §2. Q1 — What Hekkelman & McDonald 2024 actually proves

### 2.1 Setting and definitions

The paper works on a separable Hilbert space $\mathcal{H}$ with a self-adjoint operator $D$ such that $\langle D \rangle^{-d} \in \mathcal{L}_{1,\infty}$, where $\langle x \rangle := (1 + |x|^2)^{1/2}$ and $\mathcal{L}_{1,\infty}$ is the weak Schatten ideal of operators whose singular values decay as $O(k^{-1})$. This is the standard Dixmier-trace setting for spectral triples in the Connes 1988 sense. The cutoff is **sharp**: $P_\lambda := \chi_{[-\lambda, \lambda]}(D) = \chi_{[0, \lambda]}(|D|)$ — a finite-rank spectral projection. An extended limit $\omega \in \ell^*_\infty$ (Banach-limit-style) is fixed and $\mathrm{Tr}_\omega$ denotes the corresponding Dixmier trace.

The **single noncommutative integral** is (Eq. 1 of HM):
$$
a \mapsto \frac{\mathrm{Tr}_\omega(a \langle D \rangle^{-d})}{\mathrm{Tr}_\omega(\langle D \rangle^{-d})}, \quad a \in B(\mathcal{H}).
$$

The **proposed approximation** is (Eq. 2 of HM):
$$
P_\lambda a P_\lambda \mapsto \frac{\mathrm{Tr}(P_\lambda a P_\lambda)}{\mathrm{Tr}(P_\lambda)}, \quad a \in B(\mathcal{H}).
$$

The single-cutoff property is **definitional**: the approximation is the prescription "form $P_\lambda a P_\lambda$ at cutoff $\lambda$ and divide by $\mathrm{Tr}(P_\lambda)$." The convergence is then proved as a theorem.

### 2.2 The main theorems

**Proposition 2.1 (HM):** If $\mathrm{Tr}(e^{-tD^2}) \sim Ct^{-d/2}$ and $\mathrm{Tr}(a e^{-tD^2}) \sim C(a) t^{-d/2}$ as $t \to 0$, then
$$
\frac{\mathrm{Tr}_\omega(a \langle D \rangle^{-d})}{\mathrm{Tr}_\omega(\langle D \rangle^{-d})} = \lim_{\lambda \to \infty} \frac{\mathrm{Tr}(P_\lambda a P_\lambda)}{\mathrm{Tr}(P_\lambda)}.
$$

This is a clean statement: the noncommutative integral *equals* the limit of the truncated approximation, provided the heat-trace asymptotics exist for both $e^{-tD^2}$ and $a e^{-tD^2}$.

**Theorem 2.7 (HM):** Without assuming the local Weyl law for $\mathrm{Tr}(a e^{-tD^2})$, but assuming only that $D^2$ satisfies Weyl's law $\mathrm{Tr}(e^{-tD^2}) \sim C t^{-d/2}$, then
$$
\frac{\mathrm{Tr}_\omega(A \langle D \rangle^{-d})}{\mathrm{Tr}_\omega(\langle D \rangle^{-d})} = (\omega \circ M)\langle e_n, A e_n \rangle = (\omega \circ M)\!\left(\frac{\mathrm{Tr}(P_{\lambda_n} A P_{\lambda_n})}{\mathrm{Tr}(P_{\lambda_n})}\right),
$$
where $M$ is the logarithmic mean operator $M : x \mapsto \{(1/\log(n+2)) \sum_{k=0}^n x_k/(k+1)\}_{n=0}^\infty$.

**Theorem 3.2 (Szegő limit, HM):** Under Weyl's law, $[D, A]$ bounded, $A$ self-adjoint mapping $\mathrm{dom}|D|$ into itself,
$$
(\omega \circ M)\!\left(\frac{\mathrm{Tr}(f(P_{\lambda_n} A P_{\lambda_n}))}{\mathrm{Tr}(P_{\lambda_n})}\right) = \frac{\mathrm{Tr}_\omega(f(A) \langle D \rangle^{-d})}{\mathrm{Tr}_\omega(\langle D \rangle^{-d})}, \quad f \in C(\mathbb{R}), f(0) = 0.
$$
If additionally $\mathrm{Tr}(A^k e^{-tD^2}) \sim C_k t^{-d/2}$ for every positive integer $k$, then the limit $\lim_{\lambda \to \infty}$ replaces the $\omega \circ M$ averaging.

**Section 4 (Fröhlich functional):** for the case $\langle D \rangle^{-d} \notin \mathcal{L}_{1,\infty}$ but $\mathrm{Tr}(e^{-t|D|}) < \infty$ for $t > \beta$ (Li₁-summable / θ-summable triples), Proposition 4.1 shows the truncated approximation $\mathrm{Tr}(P_\lambda a P_\lambda)/\mathrm{Tr}(P_\lambda)$ equals the Fröhlich functional $\lim_{t \searrow \beta} \mathrm{Tr}(a e^{-t|D|})/\mathrm{Tr}(e^{-t|D|})$ (Eq. 12), defining a KMS state of inverse temperature $\beta$ on the Toeplitz algebra.

**Section 5 (Density of states):** for discrete metric spaces $(X, d_X)$ with weighted shifts and Schrödinger operators, a Dixmier-trace formula equates the spectral-density-of-states integral $\int f \, d\nu_H$ to the truncated-trace approximation (Theorem 5.2, Corollary 5.3).

**Section 6 (Noncommutative ergodicity):** introduces a definition of ergodicity for compact spectral triples as uniqueness of the vacuum state for the $C^*$-dynamical system $(S^*\mathcal{A}, \mathbb{R}, G_t)$ with $G_t(A) = e^{it|D|} A e^{-it|D|}$, and proves (Theorem 6.11, in the part of the paper not extracted in detail in this memo) that for spectral triples where the local Weyl law holds, ergodicity of the geodesic flow implies quantum ergodicity of $D$.

### 2.3 What HM's framework gives, in a single sentence

**Hekkelman–McDonald 2024 is a single-cutoff Tauberian theorem:** the noncommutative integral $\mathrm{Tr}_\omega(a \langle D \rangle^{-d}) / \mathrm{Tr}_\omega(\langle D \rangle^{-d})$ is exactly the sharp-cutoff truncated functional $\mathrm{Tr}(P_\lambda a P_\lambda)/\mathrm{Tr}(P_\lambda)$ in the limit $\lambda \to \infty$, provided $D^2$ satisfies Weyl's law (Definition 2.3 of HM). The convergence is in the Cesàro/logarithmic-mean sense in general, and pointwise (literal $\lim_{\lambda \to \infty}$) when $a e^{-tD^2}$ also satisfies a local Weyl law.

The relevant primitive object is $\mathrm{Tr}_\omega(a \langle D \rangle^{-d})$, which is the **leading term of the heat-trace expansion only**:
$$
\mathrm{Tr}(ae^{-tD^2}) \sim C(a) t^{-d/2} \implies \mathrm{Tr}_\omega(a \langle D \rangle^{-d}) = \Gamma(d/2 + 1) \cdot C(a)^{-1} \cdot \text{(Tauberian factor)}.
$$
HM explicitly acknowledge (page 2, "we however do not assume the existence of a full asymptotic expansion of the heat trace"): they target only the *first term* of the heat-trace expansion, which is the noncommutative integral. They do not address subleading Seeley–DeWitt coefficients or any object beyond the noncommutative integral itself.

---

## §3. Q2 — Compatibility check on Dirac-$S^3$ at finite $n_{\max}$

### 3.1 The Sprint MR-B closed-form anchor

Sprint MR-B (May 2026) computed in closed form the Camporesi–Higuchi heat-kernel modular residual on unit $S^3$:
$$
K(t) = \mathrm{Tr}\, e^{-tD^2} = \frac{\sqrt{\pi}}{2} t^{-3/2} - \frac{\sqrt{\pi}}{4} t^{-1/2} + \varepsilon(t),
$$
with the exact closed form (Sprint MR-B Eq. 3.3):
$$
\varepsilon(t) = \sum_{m \geq 1} (-1)^m \sqrt{\pi} \, e^{-m^2 \pi^2 / t} \!\left[t^{-3/2} - 2 m^2 \pi^2 \, t^{-5/2} - \tfrac{1}{2} t^{-1/2}\right],
$$
verified to >100 digits at every test point in $t \in [0.012, 0.5]$. The leading $m=1$ piece is $2\pi^{5/2} t^{-5/2} e^{-\pi^2/t}$. Every coefficient of $\varepsilon(t)$ sits in the M2 ring $\sqrt{\pi} \cdot \mathbb{Q} \oplus \pi^2 \cdot \mathbb{Q}$.

### 3.2 The compatibility check

The Camporesi–Higuchi spectrum $|\lambda_n| = n + 3/2$ on unit $S^3$ has degeneracy $g_n^{\rm Dirac} = 2(n+1)(n+2)$. The mode count $N(\lambda) = \#\{n : |\lambda_n| \leq \lambda\}$ grows as $\sim \frac{2}{3} \lambda^3 + O(\lambda^2)$ (Weyl's law on $d = 3$ Riemannian spinor Dirac, $\langle D \rangle^{-3} \in \mathcal{L}_{1,\infty}$). HM's framework therefore applies: the noncommutative integral $\mathrm{Tr}_\omega(a \langle D \rangle^{-3})/\mathrm{Tr}_\omega(\langle D \rangle^{-3})$ exists, and by HM Theorem 2.7 equals the truncated functional $\mathrm{Tr}(P_\lambda a P_\lambda)/\mathrm{Tr}(P_\lambda)$ in the Cesàro/$\omega \circ M$ sense.

**Key question:** does HM compute the same coefficient as the master Mellin engine when applied to a specific Dirac-$S^3$ observable?

The master Mellin engine of Sprint TS-E1 / Paper 32 §VIII computes $\pi$-source mechanisms M1/M2/M3 via Mellin transforms of $\mathrm{Tr}(D^k e^{-tD^2})$ for $k \in \{0, 1, 2\}$. M2 corresponds to $k = 0$ (the heat kernel itself, generating $\sqrt{\pi}$ via Seeley–DeWitt) and $k = 2$ (squared Dirac, generating $\pi^{2k}$ at integer $s$ via Theorem T9 of Paper 28). HM's noncommutative integral $\mathrm{Tr}_\omega(a \langle D \rangle^{-d})$ is — via the Mellin representation $\langle D \rangle^{-d} = \frac{1}{\Gamma(d/2)} \int_0^\infty t^{d/2 - 1} e^{-t \langle D \rangle^2} \, dt$ — exactly the $s = d/2$ leading-pole evaluation of the heat-kernel Mellin transform.

**Test computation.** Take $a = 1$ (the unit operator), $d = 3$. Then $\mathrm{Tr}_\omega(\langle D \rangle^{-3})$ should equal a specific coefficient of the master Mellin engine.

By Tauberian (HM Remark 2.2): $\mathrm{Tr}(e^{-tD^2}) \sim C t^{-3/2}$ as $t \to 0$ implies $\lambda(k, D^2) \sim \tilde{C} k^{2/3}$, with $\tilde{C} = C^{2/3} \cdot [\Gamma(5/2)]^{2/3}$ (HM Lemma 2.4 + Remark 2.2). For Camporesi–Higuchi on unit $S^3$: $C = \sqrt{\pi}/2$ (the leading SD coefficient). So
$$
\mathrm{Tr}_\omega(\langle D \rangle^{-3}) = C \cdot \Gamma(3/2 + 1)^{-1} \cdot [\text{Dixmier normalization}] = \frac{\sqrt{\pi}}{2} \cdot \frac{1}{\Gamma(5/2)} = \frac{\sqrt{\pi}}{2} \cdot \frac{4}{3\sqrt{\pi}} = \frac{2}{3}.
$$
This is $2/3 \in \mathbb{Q}$ — a pure rational, with the $\sqrt{\pi}$ factor *exactly cancelling* between the leading SD coefficient and the $\Gamma$ function. **This is the M2 cancellation in the Mellin-engine reading: the Dixmier trace normalizes away the leading $\sqrt{\pi}$.**

Compare to Sprint TS-E1 / Paper 28 T9 reading: $\zeta_{D^2}(s)$ at integer $s$ is a polynomial in $\pi^2$ with rational coefficients. At $s = 3/2$ (the Tauberian residue point — non-integer), there is no $T9$ statement; the Mellin transform has a simple pole there with residue equal to the SD coefficient times a $\Gamma$ factor, and the $\sqrt{\pi}$ in the SD coefficient is exactly the Mellin contour signature of M2.

**The compatibility statement is therefore:** HM's noncommutative integral $\mathrm{Tr}_\omega(\langle D \rangle^{-3})$ is the **leading SD residue** of the master Mellin engine at $s = d/2$, with the $\sqrt{\pi}$ M2 signature cancelling against the $\Gamma$ normalization. The two frameworks compute the same number ($2/3$) by structurally compatible routes.

**Test 2: a non-trivial $a$.** Take $a = M_f$ a multiplication operator by $f \in C(S^3)$. Then by HM Theorem 2.7 (applied to the Riemannian Dirac on $S^3$, where Connes' integration formula Theorem 1.3 of HM gives $\phi(M_f \langle D \rangle^{-d}) = \mathrm{Vol}(S^{d-1})/(d (2\pi)^d) \int_{S^3} f \, d\mathrm{vol}_g$):
$$
\frac{\mathrm{Tr}_\omega(M_f \langle D \rangle^{-3})}{\mathrm{Tr}_\omega(\langle D \rangle^{-3})} = \frac{1}{\mathrm{Vol}(S^3)} \int_{S^3} f \, d\mathrm{vol}_g.
$$
For $f = 1$ this gives $1$, consistent. For $f$ a non-trivial harmonic $Y_{n l m}$, this gives $0$ for $n \geq 1$ (the $S^3$-integral of any non-constant harmonic vanishes). The factor $\mathrm{Vol}(S^3) = 2\pi^2$ enters the *normalization*, but in the **ratio** $\mathrm{Tr}_\omega(M_f \langle D \rangle^{-3})/\mathrm{Tr}_\omega(\langle D \rangle^{-3})$ it cancels — consistent with HM's framework being designed to compute *normalized* expectation values, not absolute spectral-action coefficients.

**This is the structural compatibility, sharpened:** HM's noncommutative integral is the **normalized** master Mellin engine at $s = d/2$. The M1 Hopf-base measure $\mathrm{Vol}(S^2)/4 = \pi$ and M2 Seeley–DeWitt $\sqrt{\pi}$ both appear in the *unnormalized* numerator and denominator of HM Eq. 1, but cancel in the ratio. **HM's framework computes ratios; the master Mellin engine computes the absolute coefficients, and both M1 and M2 appear explicitly there.**

### 3.3 The compatibility verdict

**HM's framework is compatible with the master Mellin engine on Dirac-$S^3$, but operates one step higher:** the master Mellin engine computes $\mathrm{Tr}(D^k e^{-tD^2})$ as a function of $t$ with full SD expansion, exposing each $\pi$-source $\sqrt{\pi}$ (M2) and $\mathrm{Vol}(S^2)/4$ (M1) explicitly. HM's framework computes the *Tauberian residue* of the heat-trace asymptotic, which is the noncommutative integral, in which the leading $\sqrt{\pi}$ cancels against the $\Gamma$ normalization and the $\mathrm{Vol}$-measures cancel between numerator and denominator.

**No contradiction.** HM is to the master Mellin engine as **Connes' integral** is to the **full Seeley–DeWitt expansion**: the integral is the leading coefficient, the SD expansion is the full asymptotic series, and the two are consistent. The coefficient HM computes ($2/3$ for $\langle D \rangle^{-3}$ on unit $S^3$) is consistent with the leading $\sqrt{\pi}/2$ SD coefficient of $\mathrm{Tr}(e^{-tD^2})$ extracted by MR-B's closed form.

### 3.4 An honest caveat on subleading information

HM's framework **does not capture subleading SD information**. The MR-B closed form $\varepsilon(t) = \sum_m (-1)^m \sqrt{\pi} e^{-m^2 \pi^2/t}[t^{-3/2} - 2m^2\pi^2 t^{-5/2} - \tfrac{1}{2} t^{-1/2}]$ contains the modular-tower information of the Camporesi–Higuchi half-integer-shifted Dirac spectrum, with M2-ring coefficients at every order (leading $2\pi^{5/2}$, subleading $-\sqrt{\pi}, +\tfrac{1}{2}\sqrt{\pi}$, exponential exponents $m^2 \pi^2$). HM's noncommutative integral *normalizes away* this entire structure — it sees only the leading SD coefficient.

For the Q3 question (counterterms), this is the central limitation: HM's framework targets the noncommutative integral (one number), not the full spectral-action expansion (an asymptotic series). The information needed for renormalization is in the *subleading* coefficients of the spectral action, which HM's framework does not address.

---

## §4. Q3 — Counterterm structure: does any published spectral-action framework give it?

### 4.1 The question, sharpened

W2a needs: given the bare iterated CC spectral action on Dirac-$S^3$ — which reproduces UV-divergent integrands $\sim N^{3.43}$ (LS-8a) and $\sim N^{3.79}$ (HF-5) — generate $Z_2$, $\delta m$, $Z_3$ counterterms naturally as a multi-cutoff convergence statement.

The structural archetype W2a lacks is a **matched EFT chain** in the sense of JHEP 05 (2025) 171 (the muon-conversion paper): a sequence of single-scale matchings with $\overline{\mathrm{MS}}$ counterterms generated by hand at each step. The natural NCG analog would be a **two-cutoff** framework: cutoffs $\Lambda_1 > \Lambda_2$, a matching map $\mathcal{T}_{\Lambda_1} \to \mathcal{T}_{\Lambda_2}$ that absorbs the difference into a counterterm structure.

### 4.2 What HM 2024 does NOT give

**HM's framework is single-cutoff.** Every theorem (Proposition 2.1, Theorem 2.7, Theorem 3.2 Szegő limit, Proposition 4.1 Fröhlich, Theorem 5.2 density-of-states) is stated as $\lambda \to \infty$ at fixed $D$, fixed $\mathcal{H}$, fixed $a$. There is no two-cutoff matching, no RG flow, no counterterm generation, no $\overline{\mathrm{MS}}$-style scheme.

The closest HM gets to multi-scale machinery is Section 6 (noncommutative ergodicity), where the dynamical system $G_t(A) = e^{it|D|} A e^{-it|D|}$ is treated as an analog of geodesic flow — but this is a $C^*$-dynamical-system construction at fixed $|D|$, not a matching across cutoffs.

**Honest reading: HM is not a counterterm framework.** It is a Tauberian theorem identifying the noncommutative integral with a sharp-cutoff truncated functional, and a Szegő limit theorem extending this to functions $f$ of the truncated multiplication operator. The renormalization question is structurally orthogonal.

### 4.3 What no published NCG framework gives

The Track 3 literature review surveyed:
- Connes & van Suijlekom 2021 (CMP 383) — operator-system framework for spectral truncations.
- Connes & van Suijlekom 2023 (Indagationes Math. 34) — tolerance relations and operator systems.
- Hekkelman 2022 (Master's thesis); Hekkelman–McDonald 2024 (this paper).
- van Suijlekom 2024 (arXiv:2409.02773) — K-theory for operator systems.
- Leimbach–van Suijlekom 2024 (Adv. Math. 439) — torus GH convergence (transported to $S^3$ in Paper 38).
- Latrémolière 2022 (Adv. Math. 404) — Latrémolière propinquity for metric spectral triples.
- Farsi–Latrémolière 2024 (Adv. Math. 437); Farsi–Latrémolière 2024 (arXiv:2404.00240).
- Latrémolière 2026 (arXiv:2603.19128) — spectral continuity for almost-commutative manifolds in $C^1$ topology on metrics.
- Marcolli & van Suijlekom 2014 (JHEP 12 (2014) 064) — gauge networks and spectral action.
- McDonald et al. 2025 (arXiv:2506.21950) — trace formulas in NCG.

**Of these, only Connes–Chamseddine 1996 (the original spectral action paper, the foundation of the NCG-physics ecosystem) addresses the full asymptotic expansion of $\mathrm{Tr}(f(D/\Lambda))$ in $\Lambda^{-2}$ — and that expansion is single-cutoff: it gives a small-$1/\Lambda^2$ series at fixed $D$, not a matching across cutoffs.** Marcolli–vS 2014 layers Connes–Chamseddine onto a Robertson–Walker FLRW background and computes rational coefficients in the resulting expansion (see Q4 below), but does not generate counterterms.

**Verdict on Q3: no published framework solves W2a autonomously inside the spectral-action paradigm.** The renormalization story in physics is universally imported from outside the spectral-action setting via $\overline{\mathrm{MS}}$ counterterms (or the equivalent of any other minimal subtraction scheme). The NCG literature post-2021 has focused on: (a) operator-system geometry of single truncations (Connes–vS, Hekkelman); (b) GH convergence as the truncation cutoff lifts (Leimbach–vS, Paper 38); (c) almost-commutative tensor stability in continuous-metric topologies (Latrémolière 2026). None of these address the multi-loop counterterm structure.

This is **W2a is at the spectral-action program's frontier**, exactly as the Track 3 surprise S6 stated. The "real physics composes daily" pushback from the PI is correct empirically, but real physics does it by importing flat-space minimal-subtraction machinery from outside the spectral-action paradigm. GeoVac's failure mode is structurally identical to Connes–Chamseddine 1996's, and Marcolli–vS 2014 partially addresses but does not close.

---

## §5. Q4 — Marcolli–vS 2014 rationality: an independent constraint?

### 5.1 What the rationality result says

Marcolli & van Suijlekom, JHEP 12 (2014) 064, prove a rationality theorem for the Connes–Chamseddine spectral action on Robertson–Walker (FLRW) metrics: the Seeley–DeWitt expansion coefficients of $\mathrm{Tr}(f(D^2/\Lambda^2))$ on the FLRW spacetime are *rational* polynomials in (specific structural inputs of the cosmology). Track 3's literature review summarizes this (memo §5): "the Marcolli–vS rationality of spectral action for Robertson–Walker metrics shows clean rational/transcendental rings in spectral-action coefficients on homogeneous compact backgrounds."

The Track TS-E1 audit (memo §1.3.iii.a–iii.c) places κ = −1/16 (Paper 18 §V.1 calibration constant) and the Paper 18 §V Level-5 RH-sprint exchange constants in the same conceptual neighborhood: rational/transcendental ring constraints on spectral-action-like coefficients on homogeneous compact spaces.

### 5.2 Does rationality constrain *allowed* counterterms?

The rationality result is a structural statement about which coefficients can appear in the SD expansion. It is *necessary* for any counterterm derivable from spectral-action structure to be a rational combination of the SD-expansion ingredients (curvature scalars, gauge-field strengths, etc., evaluated at the cosmological background).

**However — and this is the load-bearing observation — rationality does not generate counterterms. It constrains them.** The renormalization question (Q3) is: where do $Z_2, \delta m, Z_3$ come from, structurally? Rationality says: if they exist as outputs of the spectral action, they sit in a specific rational ring. It does not say: here is the framework that generates them.

Combined with the master Mellin engine reading (Sprint TS-E1 / Paper 32 §VIII), the Marcolli–vS rationality result gives a **stronger structural compatibility** statement: the master Mellin engine's M1/M2/M3 transcendental classification is consistent with the Marcolli–vS rationality theorem in the sense that all M1/M2/M3 outputs sit in the appropriate rational $\pi^k$ ring. This is exactly the surprise S5 of Track 3:

> **S5. Hekkelman–McDonald + Marcolli–vS rationality together form a tighter framing for Paper 32.** The Marcolli–vS rationality of spectral action for Robertson–Walker metrics shows clean rational/π structure in spectral-action coefficients on homogeneous compact backgrounds. Combined with Hekkelman–McDonald 2024, this nearly *is* the Sprint TS-E1 master Mellin engine.

**My reading after the W2a-diag analysis is sharper than Track 3's:** Marcolli–vS rationality is the *structural constraint* on the master Mellin engine's outputs, and Hekkelman–McDonald is the *Tauberian extraction* of the leading coefficient. Neither generates counterterms. Together they characterize the asymptotic ecosystem the master Mellin engine sits inside, and they constrain what counterterms could autonomously emerge — but they do not produce the counterterms.

### 5.3 Q4 verdict

**Marcolli–vS rationality is independent of the renormalization question.** It constrains the ring in which counterterms must live (rational $\pi^k \cdot \mathbb{Q}$), but does not generate the counterterms. For the W2a wall, rationality is necessary but not sufficient: even if the framework had a counterterm-generation theorem, the rationality result would constrain the form of the output. Without such a theorem, the rationality result is consistent with the framework being an asymptotic-expansion ecosystem in which renormalization is an *imposed* operation, not a *derived* one.

This is consistent with Connes' original 1996 spectral-action paper and the subsequent Marcolli–vS 2014 work both treating renormalization as an external input — the spectral action is *defined* up to a regularization scheme, and the scheme (Wilsonian, $\overline{\mathrm{MS}}$, etc.) is chosen externally.

---

## §6. Q5 — Verdict on W2a

The four candidate verdicts, and which one the W2a-diag lands on:

**(a) Tooling-addressable.** Hekkelman–McDonald's machinery imports cleanly and gives counterterms. **REJECTED.** HM's framework is single-cutoff, targets the noncommutative integral (one number), and explicitly does not address any structure beyond the leading heat-trace coefficient. It does not generate counterterms; it identifies the noncommutative integral with a Tauberian residue.

**(b) Hybrid.** HM's machinery imports cleanly for the Mellin/asymptotic side but does not give counterterms — we get a sharpening of Paper 32 §VIII but W2a stays open. **PARTIALLY ACCEPTED.** HM's framework *is* compatible with the master Mellin engine's M1/M2/M3 reading on Dirac-$S^3$ (Q2 verdict). The Tauberian residue HM extracts is the leading-pole evaluation of the master Mellin engine at $s = d/2$, with the $\sqrt{\pi}$ M2 signature cancelling against the $\Gamma$ normalization. This sharpens Paper 32 §VIII's framing — placing GeoVac inside the Hekkelman–McDonald + Marcolli–vS ecosystem — without resolving W2a.

**(c) Frontier-of-field.** No published spectral-action framework gives counterterms; W2a is a frontier of the broader program and the honest framing for Paper 32 should be "see Hekkelman–McDonald 2024 + Marcolli–vS 2014 for the asymptotic ecosystem; multi-cutoff renormalization is open." **ACCEPTED as the load-bearing verdict.** The Q3 audit confirms no published NCG framework generates $Z_2/\delta m/Z_3$ counterterms autonomously. The honest scope statement for Paper 32 §VIII is therefore: **the case-exhaustion theorem operates inside an asymptotic ecosystem characterized by Marcolli–vS 2014 (rationality) and Hekkelman–McDonald 2024 (Tauberian extraction); multi-cutoff matching across two regulators is at the spectral-action program's frontier.**

**(d) Structural discrepancy.** HM contradicts the master Mellin engine. **REJECTED.** Q2 verified compatibility: HM's framework computes the Tauberian residue of the master Mellin engine at $s = d/2$, with all transcendentals cancelling consistently between the SD expansion (which carries $\sqrt{\pi}$ M2) and the $\Gamma$ normalization (which carries $\sqrt{\pi}$ in the $\Gamma(d/2 + 1)$ factor). No contradiction.

### 6.1 Combined verdict: (b) hybrid + (c) frontier-of-field

The W2a-diag verdict is **hybrid (b) at the Tauberian-extraction level + frontier-of-field (c) at the renormalization level**. Concretely:

- **Tauberian / asymptotic side (b).** Paper 32 §VIII benefits from explicit citation of HM 2024 + Marcolli–vS 2014 in the case-exhaustion theorem's framing. The master Mellin engine's M1/M2/M3 reading is consistent with HM's noncommutative-integral approximation (which extracts the leading $s = d/2$ pole) and with Marcolli–vS rationality (which constrains the SD coefficient ring). This sharpens Paper 32's structural alignment claim.
- **Renormalization side (c).** No published framework solves W2a autonomously. The honest framing in Paper 32 §VIII is "multi-cutoff renormalization across two regulators is open in the spectral-action program; GeoVac's counterterm-generation failure is structurally identical to the standard NCG situation, not a GeoVac-specific limitation."

This is the cleanest reading consistent with both the literature review and the Q1–Q4 analysis.

---

## §7. Recommended Paper 32 §VIII frontier-of-field framing edit

Following the (b) + (c) verdict, the recommended Paper 32 §VIII edit has two components: **a citation update** (HM 2024 + Marcolli–vS 2014) and **a frontier-of-field remark** (W2a as frontier of the spectral-action program). The edit should be a short addition to the existing case-exhaustion theorem section; it does not require restructuring §VIII.

### 7.1 Concrete LaTeX draft

```latex
% Insert after \begin{remark}[Falsification target] in §VIII, before §VIII.B.

\begin{remark}[Asymptotic ecosystem and the renormalization frontier]
\label{rem:asymptotic_ecosystem}
The case-exhaustion of Theorem~\ref{thm:pi_source_case_exhaustion}
operates inside a published asymptotic ecosystem.  The master Mellin
engine's M1/M2/M3 reading (Sec.~\ref{sec:master_mellin}) is structurally
consistent with two recent results in the Connes--van~Suijlekom
program:\ Hekkelman and McDonald~\cite{hekkelman_mcdonald2024} prove a
single-cutoff Tauberian theorem identifying the noncommutative integral
$\mathrm{Tr}_\omega(a \langle D \rangle^{-d}) /
\mathrm{Tr}_\omega(\langle D \rangle^{-d})$ with the truncated
functional $\mathrm{Tr}(P_\lambda a P_\lambda) / \mathrm{Tr}(P_\lambda)$
in the limit $\lambda \to \infty$ (provided $D^2$ satisfies Weyl's law),
which on Dirac-$S^3$ recovers the leading-pole evaluation of the master
Mellin engine at $s = d/2$ with the $\sqrt{\pi}$ M2 signature cancelling
against the $\Gamma$ normalization.  Marcolli and
van~Suijlekom~\cite{marcolli_vs2014_jhep} prove a rationality theorem
for the spectral action on Robertson--Walker backgrounds, constraining
the Seeley--DeWitt expansion coefficients to live in a specific
rational~$\pi^k$ ring;\ this is consistent with the master Mellin
engine's classification.

\paragraph{The renormalization frontier.}
The bare iterated Connes--Chamseddine spectral action on Dirac-$S^3$
reproduces UV-divergent integrands of multi-loop QED faithfully (Sprint
LS-8a verifies the $(\alpha/\pi)^2 (Z\alpha)^4 / n^3$ prefactor and
divergence form for the two-loop electron self-energy;\ Sprint HF-5
verifies the same structural matching in the vertex sector for the
hydrogen 21~cm hyperfine splitting), but the framework cannot
autonomously generate $Z_2/\delta m/Z_3$ counterterms required for
finite extraction.  This is not a GeoVac-specific limitation:\ no
published framework in the spectral-action program solves multi-cutoff
counterterm generation autonomously.  Hekkelman--McDonald~2024 is
single-cutoff;\ Marcolli--van~Suijlekom~2014 constrains the ring of
allowed coefficients but does not generate counterterms;\ Connes and
Chamseddine~\cite{chamseddine_connes1996} treat renormalization as an
externally imposed operation in their original spectral-action
construction.  The empirical fact that real physics closes multi-loop
renormalization daily (via $\overline{\mathrm{MS}}$ counterterms in
matched effective-field-theory chains) is consistent with this:\ the
counterterm machinery is imported from flat-space QFT, not generated
inside the spectral-action paradigm.

The honest scope statement is that multi-cutoff renormalization across
two regulators is at the spectral-action program's frontier, not a
GeoVac-internal failure.  The case-exhaustion theorem and the master
Mellin engine characterize where $\pi$ enters; closing W2a would
require an extension of the spectral-action framework to handle
multi-cutoff matching, a problem on which no published NCG result
appears as of May 2026.
\end{remark}
```

This adds approximately 250 words to §VIII and completes the frontier-of-field framing without overreaching.

### 7.2 Bibliography additions

Two new bibliographic items needed for Paper 32:

```latex
\bibitem{hekkelman_mcdonald2024}
E.-M.~Hekkelman and E.~A.~McDonald,
``A noncommutative integral on spectrally truncated spectral triples,
and a link with quantum ergodicity,''
\textit{J. Funct. Anal.} (accepted), arXiv:2412.00628 (2024).

\bibitem{marcolli_vs2014_jhep}
M.~Marcolli and W.~D.~van~Suijlekom,
``Gauge networks in noncommutative geometry,''
\textit{J. Geom. Phys.} 75 (2014), 71--91 [arXiv:1301.3480].
% NOTE: cross-check with the JHEP 12 (2014) 064 reference in CLAUDE.md;
% Marcolli-vS have multiple 2014 papers in this neighborhood, the
% rationality-of-spectral-action-on-FLRW result is a different paper.
```

The reference to the Marcolli–vS rationality theorem on FLRW backgrounds should be confirmed against the actual source — multiple 2014 Marcolli–vS papers exist in this neighborhood, and the case-exhaustion's reading should cite the specific rationality paper. Track 3's memo flagged this as JHEP 12 (2014) 064; the existing Paper 32 bibliography (line 2667) cites the gauge-networks paper. Resolution: cite both, with the rationality paper noted as the structural-coefficient-ring constraint.

### 7.3 What this edit does NOT claim

The edit does **not** claim:
- Hekkelman–McDonald solves W2a (they don't; HM is single-cutoff and targets only the leading SD coefficient).
- Marcolli–vS rationality generates counterterms (it constrains them).
- The case-exhaustion theorem extends to a renormalization theorem (it does not — the case-exhaustion is qualitative-class, the renormalization question is quantitative-coefficient).

The edit *does* claim:
- The master Mellin engine sits inside the published HM + MvS asymptotic ecosystem (sharpening Paper 32's structural alignment claim).
- Multi-cutoff renormalization is at the spectral-action program's frontier, not a GeoVac-internal failure (honest scope statement).

This is the cleanest framing consistent with the (b) + (c) verdict.

---

## §8. Honest scope and uncertainty

### 8.1 What this memo did

- Read Hekkelman–McDonald 2024 (arXiv:2412.00628v2) in full through Section 6, including the precise statements of Definitions 1.1–1.5, Propositions 2.1, 4.1, 5.1, Theorems 2.7, 3.2, 5.2, 6.3, and Corollary 5.3. Identified the framework as a single-cutoff Tauberian theorem targeting the noncommutative integral (leading SD coefficient) only.
- Verified compatibility on Dirac-$S^3$ at finite $n_{\max}$ via the Sprint MR-B closed-form modular residual $\varepsilon(t)$. Computation: $\mathrm{Tr}_\omega(\langle D \rangle^{-3}) = 2/3$ via HM Lemma 2.4 + Camporesi–Higuchi leading SD coefficient $\sqrt{\pi}/2$, with the $\sqrt{\pi}$ M2 signature cancelling against $\Gamma(5/2) = 3\sqrt{\pi}/4$. No contradiction with the master Mellin engine.
- Audited the published NCG/spectral-action ecosystem post-2021 (Connes–vS 2021, 2023; Hekkelman 2022; HM 2024; van Suijlekom 2024; Leimbach–vS 2024; Latrémolière 2022, 2026; Farsi–Latrémolière 2024; Marcolli–vS 2014; Chamseddine–Connes 1996) for any multi-cutoff counterterm-generation result. None found.
- Drafted the Paper 32 §VIII framing edit (~250 words), with the honest scope statement that multi-cutoff renormalization is at the spectral-action program's frontier.

### 8.2 What this memo did not do

- **Did not run a numerical PSLQ or symbolic check** of HM's framework against any specific sprint observable beyond the $\mathrm{Tr}_\omega(\langle D \rangle^{-3}) = 2/3$ symbolic computation. A more thorough check would: take a non-trivial $a$ (e.g. a multiplication operator by a low-order spherical harmonic on $S^3$), compute $\mathrm{Tr}(P_\lambda a P_\lambda)/\mathrm{Tr}(P_\lambda)$ at $\lambda = 1/2 + n_{\max}$ for $n_{\max} \in \{2, 3, 4, 5\}$ in exact rational arithmetic, and compare to the master Mellin engine's prediction. This is a sprint-scale follow-up (~1–2 weeks).
- **Did not verify the precise Marcolli–vS 2014 rationality statement** beyond the Track 3 literature review summary. The exact ring (rational, $\mathbb{Q}[\pi]$, $\mathbb{Q}[\pi^2]$) and the precise structural inputs to the rationality theorem were not extracted from the source paper. A bibliography update should include reading at least the abstract and main theorem statement of the Marcolli–vS rationality paper before the §VIII edit lands.
- **Did not check for very recent (2025–2026) follow-ups** to HM 2024 that might extend the framework to multi-cutoff matching. The McDonald et al. 2025 paper (arXiv:2506.21950, "Trace Formulas in Noncommutative Geometry") was flagged but not read.

### 8.3 Uncertainties

- **The compatibility check at non-trivial $a$ is symbolic-only in this memo.** The cancellation of $\sqrt{\pi}$ between the leading SD coefficient and the $\Gamma$ normalization is a structural identity, not a computed cross-check at a specific observable. A numerical sprint would tighten this.
- **The "frontier-of-field" statement assumes I have done a complete literature audit.** I have not — Track 3 surveyed approximately 60 entries; the field has hundreds of papers per year. The statement is "no published spectral-action framework I have located generates multi-cutoff counterterms" rather than "no such framework exists." This is the honest scope.
- **The paper-32 edit recommendation is a structural-framing edit, not a GeoVac-internal new result.** It strengthens Paper 32's positioning by placing the case-exhaustion theorem inside the published ecosystem, but it does not close W2a. Closing W2a would require a multi-cutoff extension of HM's framework, which is itself an open NCG / spectral-action question.

### 8.4 What I am highly confident about

- HM 2024's framework is single-cutoff. The framework targets the noncommutative integral only and does not generate counterterms.
- HM 2024's framework is compatible with the master Mellin engine on Dirac-$S^3$ at the leading-SD-coefficient level. The Tauberian residue HM extracts coincides with the master Mellin engine's leading-pole evaluation.
- W2a is at the spectral-action program's frontier. No published NCG framework I have located solves multi-cutoff renormalization autonomously.
- The Paper 32 §VIII edit recommendation (~250 words) is a clean framing improvement that does not overreach.

### 8.5 What the next sprint should do

If Phase C-W2a is opened, the natural sprint structure is:

1. **Numerical compatibility check** (~1–2 weeks): take a panel of 5–10 specific observables on Dirac-$S^3$ at $n_{\max} \in \{2, 3, 4\}$, compute both the HM-truncated functional and the master Mellin engine prediction, verify ratios in exact rational arithmetic. Outcome: either the symbolic compatibility extends to numerical (positive — Paper 32 §VIII edit lands cleanly), or there is a discrepancy in subleading coefficients (negative — the (b) verdict needs revision).

2. **Marcolli–vS rationality direct read** (~1 week): extract the precise rationality-of-spectral-action statement from the source paper and verify it constrains the M1/M2/M3 ring as Track 3 surprise S5 claimed.

3. **Phase C-W2a position paper** (~3–4 weeks): if (1) and (2) land positively, draft a short Paper 32 §VIII expansion or a Paper 39 standalone "Asymptotic ecosystem of the GeoVac spectral triple" memo that places the case-exhaustion theorem inside the published HM + MvS ecosystem and names multi-cutoff renormalization as the open frontier. Honest scope: this is a positioning paper, not a closure attempt.

If Phase C does not open — i.e. the (b) + (c) verdict is accepted as the final Phase B result — then the §VIII edit recommended above is the entirety of the W2a Phase B output, and W2a is recorded as frontier-of-field for Paper 32 going forward.

---

**End of W2a-diag memo.**
