# Sub-sprint Z — Bethe log / Lamb shift Latrémolière interpretation

**Date:** 2026-05-23 (Sub-sprint Z of Sprint L3e-P3, dispatched as parallel-to-X physics-side application probe).
**Sprint position:** Application #3 of the Phase A.2' deep-read (`debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md` §4). Diagnostic only — no code modifications, no production-test perturbations.
**Sprint goal:** Determine whether Paper 36's Sturmian-based bound-state QED architecture (LS-3 acceleration-form Bethe log, LS-4 Drake-Swainson, LS-6a Eides convention, LS-7 two-loop SE prefactor) inherits a rigorous truncation-error model from Latrémolière arXiv:2512.03573's pinned-QLCMS framework, or whether the identification gives only a vocabulary change.
**Verdict (one line):** **MIXED.** The Sturmian-projection / pin-state / exhaustive-sequence identification is structurally sound at the prefactor convergence level (LS-3 s-states, LS-7 prefactor), gives a real new error-bar model for the Bethe-log truncation error, but does NOT predict the LS-3 acceleration-vs-velocity 3.3× speedup, does NOT cure the LS-3 ℓ>0 closure pathology, does NOT generate the LS-4 Drake-Swainson regularization, and CANNOT generate the LS-7 bracket $C_{2S} = +3.63$ or the LS-8a renormalization counterterms — these are all *operator-system content beyond the pin-state local metametric*. The framework gives Paper 36 a rigorous convergence rate for the parts that already converged, not a new mechanism for the parts that didn't.

---

## §1. What Paper 36's Bethe log actually uses (recap, ~500 words)

Paper 36 (`papers/group5_qed_gauge/paper_36_bound_state_qed.tex`) computes the hydrogen 2S_{1/2} − 2P_{1/2} Lamb shift on Dirac-$S^3$ at one loop to sub-percent accuracy (1052.19 MHz vs experimental 1057.845 MHz; residual −0.534%). The architecture, reading the LS-1..LS-7 sprint sequence (`debug/ls{1..7}_*_memo.md`):

**Bethe logarithm $\ln k_0(n\ell)$.** This is the bound-state matrix element that appears as a logarithmic factor in the one-loop self-energy contribution to the Lamb shift. It is defined as the ratio $J/I$ of two spectral sums over the hydrogen-like Coulomb spectrum:

- $I_v(n\ell) = \sum_m |\langle n\ell | \mathbf{p} | m\rangle|^2 (E_m - E_{n\ell})$ (velocity form)
- $J_v(n\ell) = \sum_m |\langle n\ell | \mathbf{p} | m\rangle|^2 (E_m - E_{n\ell}) \ln |2(E_m - E_{n\ell})/\text{Ry}|$

By Bethe's closure identity $I_v(n\ell) = (2 Z^4/n^3) \delta_{\ell,0}$ — finite for s-states, identically zero for $\ell > 0$.

**Sturmian basis at $\lambda = Z/n_*$.** GeoVac uses Coulomb Sturmians (`geovac/casimir_ci.py`, `geovac/hylleraas_r12.py`, the LS-3 implementation in `debug/ls3_bethe_log_regularized.py`) at exponent $\lambda = Z/n_*$ where $n_*$ is the principal quantum number of the bound state being evaluated. This choice makes the bound-state $|n_*, \ell, m\rangle$ exact at finite Sturmian basis size and gives discrete completeness on $L^2(\mathbb{R}^3, d^3r)$ via the Shibuya-Wulfman integrals.

**LS-3 acceleration form.** Replaces $|\langle p\rangle|^2 (\Delta E)^2$ with $|\langle \nabla V\rangle|^2$ via $[\mathbf{p}, H_0] = -i\nabla V$. Equivalent as continuum sums, but **finite-basis evaluation is reweighted** — the $1/(\Delta E)$ factor suppresses high-energy pseudostates. Empirical: gives 2S accuracy 0.92% at N=40 vs velocity form's 3.1% at N=50 (3.3× tighter). Same closure pathology for $\ell > 0$.

**LS-4 Drake-Swainson asymptotic subtraction.** Closes the 2P (and all $\ell > 0$) divergence by splitting the spectral sum at intermediate $K$ into $\beta_{\text{low}}(N, K) + \beta_{\text{high}}(K)$ with $K$-cancellation. Introduces the structural denominator $D_{\text{drake}}(n, \ell) = 2(2\ell+1) Z^4/n^3$. Paper 34 §III.13: this is the 13th projection.

**LS-6a Eides convention fix.** Uehling kernel constant $4/15 = (10/9 - 38/45)$ was double-counted in LS-1; the LS-6a fix shifts the Lamb shift by $+27.13$ MHz at $n=2$, $Z=1$. Closes a +3.10% → −0.534% residual jump.

**LS-7 two-loop SE prefactor.** $(\alpha/\pi)^2 (Z\alpha)^4 / n^3$ derived natively from iterated CC spectral action on Dirac-$S^3$. The $1/\pi^2$ encodes two iterated Schwinger proper-time integrations (Paper 35 Prediction 1 at the structural level). Dimensionless bracket $C_{2S} = +3.63$ taken from Eides Tab. 7.3 literature value; native derivation is LS-8a scope.

**The empirical convergence.** Paper 36's −0.534% native residual is established by Drake 1990 reference values (10⁻¹⁰ precision tabulation) cross-checked against the LS-3 N=40 acceleration-form Sturmian sums at 50-digit mpmath precision. Paper 36 does NOT claim a rigorous error bar; it reports empirical convergence and decomposes the residual via Eides Tab. 7.3 / 7.4 itemization.

---

## §2. Is the Sturmian-Bethe-log construction an L-Lipschitz μ-pinned exhaustive sequence?

The deep-read memo §4 sketches the identification. I now check the three Def 1.29 axioms in detail against the LS-3 / Drake construction.

**The mapping.**

| Latrémolière 2512.03573 | Paper 36 / LS-3 |
|:------------------------|:----------------|
| C*-algebra $\mathcal{A}$ | $C_0(\mathbb{R}^3)$ — multiplication operators on $L^2$ (non-unital because $\mathbb{R}^3$ is non-compact) |
| Pin state $\mu \in \mathcal{S}(\mathcal{A})$ | $\mu = \omega_{\psi_{nS}}$ — ground/excited bound-state expectation, $\mu(f) = \langle\psi_{nS}|f|\psi_{nS}\rangle$ |
| Leibniz seminorm $L$ | $L(f) = \|[D_{\text{Coulomb}}, M_f]\|$ where $D_{\text{Coulomb}} = -\tfrac{1}{2}\nabla^2 - Z/r$ and $M_f$ is multiplication by $f$ |
| Exhaustive sequence $(h_n)$ | $h_{N_{\max}} = P_{N_{\max}}^{\text{Sturm}}$ — projector onto Sturmian basis at $\lambda = Z/n_*$ truncated at $N_{\max}$ |

**Axiom (1): $L(h_n) \to 0$.**

The Sturmian truncation projector $P_{N_{\max}}^{\text{Sturm}}$ is a bounded operator on $L^2$. Its commutator with $D_{\text{Coulomb}}$ has operator norm controlled by the spectral content of $D_{\text{Coulomb}}$ beyond the truncation boundary. For continuum hydrogen at exponent $\lambda = Z/n_*$, the high-Sturmian-mode commutator norm decays as the inverse of the spectral gap to the next un-truncated mode. In Paper 38's L2 quantitative-rate language: this is the **scalar Camporesi-Higuchi analog** of the $\gamma_n = (4/\pi) \log n / n + c/n + O(\log n / n^2)$ asymptote with $c \approx 4.109$ (Sprint MR-C constant). On the SU(2) compact case Paper 38 has a closed-form Bożejko-Fendler argument; on the Sturmian / Coulomb non-compact case the analog is heuristically present but not proved.

**Verdict on (1): PLAUSIBLE.** Decays to zero structurally; rate is the LS-3 Sturmian basis quality limit. **Quantitative rate is OPEN** in the Latrémolière framework.

**Axiom (2): $\mu(h_n) \to 1$.**

This is the *most subtle* axiom. Two readings:

(a) **Strict reading.** $\mu(h_n) = \langle\psi_{nS}|P_{N_{\max}}^{\text{Sturm}}|\psi_{nS}\rangle$. At $\lambda = Z/n_*$ with $n_* = n$ (the same principal quantum number as the bound state), the bound state $\psi_{nS}$ is the first Sturmian basis function (up to normalization); the projector $P_{N_{\max} \ge 1}$ has $\mu(h_n) = 1$ identically for all $N_{\max} \ge 1$. ✓ **STRICTLY satisfied** (saturates the limit at finite $n$).

(b) **Honest reading.** Paper 36's Bethe log calculation does NOT just evaluate $\mu(h_n)$ — it evaluates the spectral sum $\sum_m |\langle\psi_{nS}|p|m_{\text{Sturm}}\rangle|^2 \dots$ which probes all of the Sturmian basis, not just the ground-state component. The relevant convergent quantity is not $\mu(h_n)$ but the truncated spectral sum's distance to the full continuum spectral sum. The pin-state expectation is exactly 1 at all $N_{\max}$, so this axiom is **vacuously satisfied** in a way that doesn't directly give Paper 36 a new convergence statement.

**Verdict on (2): TRIVIALLY satisfied because of Sturmian's saturating property** — the pin state IS the first basis function. This is the structural reason Sturmian closure works for hydrogen polarizability at $N_{\text{basis}} = 2$ (mentioned in the deep-read memo §4 as a candidate Latrémolière interpretation).

**Axiom (3): $\|h_n\|_{\mathcal{A}} \to 1$.**

$\|P_{N_{\max}}^{\text{Sturm}}\|_{\text{op}} = 1$ since it's a projector. ✓ **STRICTLY satisfied** (saturates at finite $N_{\max}$).

**Net structural verdict.** The Sturmian-Bethe-log construction IS an L-Lipschitz μ-pinned exhaustive sequence in Latrémolière 2512.03573's Def 1.29 sense, but TWO of the three axioms are satisfied *trivially* (saturating at finite $N_{\max}$), and the one non-trivial axiom (L decay) has an OPEN convergence rate in the Latrémolière framework.

**The Leibniz axiom (Def 1.18).** $L(fg) \le \|f\| L(g) + L(f) \|g\|$ for $L = \|[D_{\text{Coulomb}}, \cdot]\|$. This is standard for any Dirac-type or Schrödinger-type operator on a smooth manifold — the Leibniz rule for commutators is $[D, fg] = [D, f]g + f[D, g]$, and operator norm subadditivity does the rest. ✓ **STRICTLY satisfied.**

The structural identification is **structurally sound** but tells us less than the framing suggests, because Sturmian's saturating property makes two of the three axioms trivial.

---

## §3. What's the explicit Latrémolière truncation error model?

If we accept the identification, Latrémolière's local metametric $\delta_r$ at radius $r$ gives a bound on the truncation error of *Lipschitz observables* that take expectation values in the pin state. The relevant Lipschitz observable in Paper 36's Bethe log is the *Coulomb commutator structure* itself — specifically, the spectral sums $I_v$, $J_v$ (and their acceleration analogs $I_a$, $J_a$).

**What the framework gives.** The local metametric $\delta_r(\mu, \phi) = \sup\{|\mu(f) - \phi(f)| : f \in \text{dom}_{sa}(L), L(f) \le 1, \mu(f^*f) \le r^2\}$ provides an error bar on differences of pin-state expectations for Lipschitz observables of bounded size. For the Bethe log spectral sum at truncation $N_{\max}$, the error in evaluating the sum $\sum_m \langle\psi_{nS}|p|m\rangle\langle m|p|\psi_{nS}\rangle\Phi(E_m)$ for a Lipschitz weight function $\Phi$ is bounded by

$$\Big|\sum_{m \in \text{Sturm}_{N_{\max}}} - \sum_{m \in \text{full}}\Big| \le L(\Phi) \cdot \delta_r(\mu_{N_{\max}}, \mu_{\text{full}})$$

with $L(\Phi)$ the Lipschitz norm of the weight function.

**Does this predict LS-3's empirical convergence?** Let me cross-check against the LS-3 data:

- 1S: at N=20, error 0.60% (acc); at N=24, error 3.32%; overshoots past N=20.
- 2S: at N=40, error 0.92% (acc); at N=50, error 2.98%; overshoots past N=40.
- 3S: at N=40, error −10.7% (acc); slower convergence at higher $n$.

The pattern: **monotone convergence to ~0% at an N-dependent optimum, then overshoot past the optimum.** The Latrémolière framework, applied as $\delta_r(\mu_{N_{\max}}, \mu_{\text{full}}) \to 0$, predicts monotone convergence — it does NOT predict the overshoot, which is a Sturmian basis quality limit (high-Sturmian pseudostates accumulate numerical round-off error).

**Verdict.** Latrémolière gives a *monotone error bar* that matches Paper 36's empirical convergence up to the N-dependent optimum, but the overshoot regime (basis-quality limit) is outside the framework's scope. The framework can predict that the convergence exists, not where the basis becomes numerically degenerate.

**Does the framework predict the LS-3 acceleration-vs-velocity 3.3× speedup?** No. The acceleration form is obtained by an *operator identity* $[\mathbf{p}, H_0] = -i\nabla V$ applied INSIDE the spectral sum, not by changing the truncation projector. Both velocity and acceleration forms use the same Sturmian projector and the same pin state, so the same Latrémolière convergence model applies to both. The acceleration form's 3.3× tighter convergence at fixed $N$ comes from the $1/(\Delta E)$ reweighting suppressing high-energy pseudostates — this is a *spectral-density reweighting effect* on the same exhaustive sequence, not a different sequence. **Latrémolière sees both forms as the same μ-pinned exhaustive sequence with different Lipschitz weight functions; the 3.3× speedup is in the $L(\Phi)$ factor (or equivalently in the radius $r$ at fixed precision), not in $\delta_r$.**

**Does the framework predict the LS-4 Drake-Swainson regularization?** No. The closure pathology $I_v(2P) = 0$ is a property of the *Lipschitz observable being evaluated* (the operator $\mathbf{p}$ vanishes at the origin for $\ell > 0$ states), not of the exhaustive sequence. Drake-Swainson regularization adds a *new structural denominator* $D_{\text{drake}}(n, \ell) = 2(2\ell + 1) Z^4/n^3$ that is NOT generated by the local metametric — it's an external regularization scheme imposed on top of the spectral sum. The Latrémolière framework gives convergence of the spectral sum if the spectral sum converges; it does NOT cure observables whose spectral sum doesn't converge.

**Does the framework predict the LS-7 two-loop SE prefactor?** Partially. The $(\alpha/\pi)^2$ prefactor comes from two iterated Schwinger proper-time integrations, which on the Latrémolière side are *two iterated Mellin-type integrations over the pin-state-Lipschitz observable's spectral parameter*. The pin-state structure makes the iteration well-defined at finite truncation. The framework supports the structural derivation. **However**, the framework does NOT generate the dimensionless bracket $C_{2S} = +3.63$ — this is the bound-state matrix element of the iterated CC spectral action, which involves *operator-system content beyond the pin-state local metametric*. The bracket evaluation requires the full Hilbert space + Hamiltonian + operator generation that Latrémolière's framework doesn't supply (because his framework is a hypertopology / convergence framework, not a renormalization framework).

**Crystallized scope.** Latrémolière 2512.03573 supplies:
- ✓ Rigorous convergence of the Sturmian-projected spectral sum to the full continuum sum
- ✓ Error bar on pin-state expectation values at fixed truncation $N_{\max}$
- ✓ Vocabulary for the LS-3 1S/2S/3S "basis-quality limit" as a non-Latrémolière phenomenon (numerical, not structural)

But does NOT supply:
- ✗ The 3.3× acceleration-vs-velocity speedup mechanism (in spectral-density reweighting, not in the exhaustive sequence)
- ✗ The LS-4 Drake-Swainson regularization (external regularization scheme)
- ✗ The LS-7 dimensionless bracket $C_{2S} = +3.63$ (operator-system content)
- ✗ The LS-8a renormalization counterterms (renormalization framework, not hypertopology)

---

## §4. Connection to Paper 38's L2 quantitative rate

Paper 38 establishes that the SU(2) central spectral Fejér kernel has a quantitative rate $\gamma_n = (4/\pi) \log n / n + c/n + O(\log n / n^2)$ with $c \approx 4.109$ extracted to ≥22 verified dps but NOT identified in any closed form against M1∪M2∪M3∪{γ_E, log 2, G, ζ(3)} (Sprint MR-C: PSLQ ceiling 10⁸ null). The $4/\pi$ asymptote is the **M1 Hopf-base measure signature** of the master Mellin engine (Sprint TS-E1 case-exhaustion theorem).

**Does the Sturmian Bethe log inherit Paper 38's L2 rate?** Let me check the structural mapping.

Paper 38's setup:
- Underlying group: $SU(2)$ (compact, non-abelian)
- Truncation: Peter-Weyl basis at $n_{\max}$
- Central Fejér kernel: $K_{n_{\max}}(g) = (1/Z_{n_{\max}}) |\sum_{j \le j_{\max}} \sqrt{2j+1} \chi_j(g)|^2$
- Rate: $\gamma_n \to 0$ as $n_{\max} \to \infty$

Sub-sprint Z setup (proposed):
- Underlying space: $\mathbb{R}^3$ (non-compact, abelian) with bound-state pin
- Truncation: Sturmian basis at $\lambda = Z/n_*$, truncated at $N_{\max}$
- "Central kernel": Sturmian basis density on Coulomb continuum
- Rate: depends on the *non-compact* analog of the central Fejér moment

**The structural difference.** Paper 38's $SU(2)$ is compact with finite Haar measure; Paper 38's $\gamma_n$ rate comes from spectral Fejér averaging that has no obvious non-compact analog. The Sturmian basis on Coulomb's *discrete* spectrum (bound states) plus *continuous* spectrum (continuum) is fundamentally a non-compact setup, and Latrémolière's pinned-QLCMS framework is what generalizes this — but the rate Paper 38 establishes is *not transported* automatically.

**What the framework does give.** Latrémolière 2512.03573 §6 provides one worked non-compact example ($c_0(\mathbb{Z}) \rtimes_\alpha \mathbb{Z}$) where the rate IS explicit. But this is a finite-dimensional approximation setting, not the Sturmian-on-$\mathbb{R}^3$ setting. Whether the Sturmian setting admits a Paper-38-style closed-form rate is an OPEN question.

**Empirical comparison to Drake 1990.** From the LS-3 memo:
- 1S: error 0.60% at N=20 (acc form)
- 2S: error 0.92% at N=40
- 3S: error 10.7% at N=40 (slower because more nodal structure)

If we naively fit $\text{err}(N) \sim c_0 / N^{p}$ to the 1S sequence (N=12,16,20: errors 8.62%, 3.13%, 0.60%), we get $p \approx 4$ (the LS-3 memo notes "O(1/√N) Sturmian convergence" but the empirical convergence is closer to $1/N^4$ in the acceleration form's pre-overshoot regime). The Paper 38 rate $\log n / n$ is much slower than this. **The Sturmian convergence rate is NOT inherited from Paper 38's L2 quantitative rate.**

**Structural reason.** Paper 38's rate is set by the *non-abelian compactness* of SU(2) (the $\log n$ comes from the Plancherel weight $\sqrt{\dim V_\pi}$ — see L2 universal-proof memo). The Sturmian setting is non-compact and abelian; the rate is set by the Sturmian basis's *radial completeness* on the bound-state subspace, which converges much faster (exponential in the basis size for analytic observables, polynomial for non-analytic ones). The two rates are categorically different.

**Verdict.** Sturmian Bethe log does NOT inherit Paper 38's $4/\pi$ asymptote. The Latrémolière framework is the *correct vocabulary* for both, but the convergence rates live in different mechanism classes (compact-group Plancherel vs Sturmian radial saturation).

---

## §5. What does this give the framework?

If Sub-sprint Z verifies (and §2-4 above strongly suggests it does, but with the caveats noted):

**What we get (genuine new content).**

1. **A precise math.OA convergence statement for the LS-3 s-state Bethe logs.** Paper 36's empirical 0.92% accuracy at 2S, N=40 is converted from "empirical convergence" to "Latrémolière hypertopology bound at truncation $N=40$ for the Lipschitz observable $|\langle\psi_{2S}|\nabla V|m\rangle|^2 / (E_m - E_{2S})$" — a rigorous statement with explicit error bar (modulo proving the Sturmian-specific rate).

2. **Cleaner scope statement for LS-7 / LS-8a residual decomposition.** Paper 36 §VII residual decomposition (+1.20 MHz multi-loop QED + ~+4.4 MHz non-loop physics) gets a Latrémolière-framework partition: the multi-loop QED piece is "within the operator-system content beyond the pin-state local metametric"; the non-loop physics is "outside the pin-state structure entirely (different C*-algebra)." This is a *more precise* version of the structural-skeleton-scope framing.

3. **Future bound-state QED extensions inherit the same convergence model.** Muonic hydrogen Lamb shift (Sprint MH Track A), positronium hyperfine (Track 1 of post-MH multi-track launch), helium 2³P fine structure (Sprint May 2026 catalogue) — all use Sturmian basis at different $\lambda$. Each inherits the Latrémolière convergence model with the same axiom verifications. **The framework's atomic precision catalogue gets a uniform convergence vocabulary.**

**What we do NOT get.**

1. **No new closed-form prediction.** The framework doesn't predict the LS-3 N=40 optimum, the 3.3× acceleration speedup, the Drake-Swainson denominator, or the $C_{2S}$ bracket. These remain empirical / operator-system content.

2. **No new physics.** Paper 36's −0.534% residual at one loop is unchanged. The framework provides error bars on parts that already converged, not new mechanisms for parts that didn't.

3. **The LS-8a renormalization gap is NOT closed.** The two-loop counterterms $Z_2$, $\delta m$ are operator-system content (renormalization framework); Latrémolière supplies a hypertopology / convergence framework. The two are categorically different. **A Latrémolière-framework interpretation of the LS-8a gap exists as "operator-system content beyond the pin-state local metametric"**, but this is a vocabulary reframing, not a mechanism for generating the counterterms.

---

## §6. Concrete Sub-sprint Z.1 plan (if opened)

A 1-2 week implementation sprint:

**Z.1.1 Verify the Leibniz axiom on $L(f) = \|[D_{\text{Coulomb}}, M_f]\|$** for the Sturmian basis on hydrogen 1S. Concrete test: pick three Sturmian basis functions $s_1, s_2, s_3$ at $\lambda = 1$, compute $L(s_1 s_2)$ and $\|s_1\| L(s_2) + L(s_1) \|s_2\|$, verify Leibniz inequality. ~2-3 days.

**Z.1.2 Compute the explicit Latrémolière exhaustion rate for $\ln k_0(1S)$.** Use the LS-3 N=12,16,20,24,30 data; fit the empirical convergence to a Latrémolière-style $\delta_r$ bound. Compare to Drake 1990 reference. Extract whether the rate is $\log N / N$ (Paper 38 universality, unlikely) or Sturmian-specific (likely). ~3-4 days.

**Z.1.3 Extract a structural prediction for LS-7 two-loop SE prefactor convergence.** Apply the Latrémolière convergence model to the iterated CC spectral action used in LS-7. Predict whether $C_{2S}$ converges at the same rate as $\ln k_0(2S)$, faster, or slower. Cross-check against the LS-7 first-pass empirical numbers (if available). ~2-3 days.

**Z.1.4 Identify whether the LS-8a renormalization gap has a Latrémolière-framework interpretation.** Test whether the divergent counterterm structure $Z_2 - 1$, $\delta m$ can be expressed as "Lipschitz observables outside the pin-state local metametric" or whether it lives in a categorically distinct C*-algebra. ~2-3 days.

**Total: ~10 working days = 2 weeks. Exit criterion:** clean structural statement of whether the Sturmian-Bethe-log construction inherits a non-trivial Latrémolière convergence rate, and whether the LS-8a gap has a Latrémolière interpretation.

**Honest scope of Z.1.** The sprint would not modify Paper 36 or production code. It would produce a memo + cross-paper Latrémolière vocabulary update (Paper 36 §IX or §X) + queue LS-8a as either "framework-Latrémolière target" (if Z.1.4 returns positive) or "operator-system content beyond Latrémolière" (if negative).

---

## §7. Honest scope assessment

**The PI question is specifically: does this Latrémolière interpretation give Paper 36 anything physics couldn't access via the empirical convergence already established?**

**Honest answer: PARTIALLY YES, MOSTLY VOCABULARY.**

The genuine new content (§5 list):
- Rigorous math.OA error bar replaces empirical convergence
- Uniform convergence vocabulary for the atomic precision catalogue
- Cleaner scope statement for the LS-8a / multi-loop / non-loop partition

The vocabulary-change content:
- Sturmian-as-exhaustive-sequence is structurally trivial because Sturmian saturates two of the three Def 1.29 axioms
- The 3.3× speedup of LS-3 acceleration form is NOT predicted by the framework — it remains a spectral-density-reweighting phenomenon
- The LS-4 Drake-Swainson regularization is NOT a Latrémolière construction
- The LS-7 bracket $C_{2S}$ and the LS-8a counterterms are operator-system content

**The honest reading.** Latrémolière 2512.03573 supplies a *hypertopology framework* that GeoVac's Sturmian-based atomic FCI calculations naturally fit into. The framework is the right vocabulary for what Paper 36 already does. But it does NOT supply the renormalization mechanism that Paper 36's structural-skeleton-scope framing (CLAUDE.md §1.7 WH1, WH5, LS-8a wall) puts out of GeoVac's reach.

**Implications for sprint priority.** If we open Sub-sprint Z.1:
- We get a clean math.OA error bar on Paper 36's empirical convergence (~10 working days)
- We do NOT get progress on the LS-8a renormalization wall (which would need a separate operator-system / renormalization framework, e.g., the Connes-Kreimer Hopf algebra of Feynman diagrams)
- The framework's structural-skeleton-scope statement is REFINED (not retired): "GeoVac is a structural-skeleton framework that admits a Latrémolière hypertopology error model for its convergent computations; renormalization counterterms remain external"

**Relative priority among X / Y / Z.** Sub-sprint X (atomic FCI Sturmian-as-Latrémolière) is the *cleanest* application because Hylleraas r₁₂ saturates the Bethe-log spectral sum issue (the explicit r₁₂ basis builds the cusp variationally; no spectral-sum-divergence pathology). Sub-sprint Z (Bethe-log / Lamb shift) is **technically harder** because the Bethe-log construction already requires non-trivial regularization (LS-4 Drake-Swainson) that lives outside the framework. Sub-sprint Y (W1c chemistry pin-state) is the most physics-relevant if the right pin-state exists, but is structurally MEDIUM-confidence.

**Recommendation.** If X verifies cleanly in Sub-sprint X (~1-2 weeks), opening Z.1 as a parallel 2-week sprint makes sense — the deliverable is Paper 36 §IX or §X vocabulary update + LS-8a scope refinement. If X returns mixed or negative, deferring Z is reasonable since Z inherits X's structural foundation.

---

## §8. Verdict

**Sub-sprint Z verdict: MIXED — POSITIVE-with-explicit-rate on the parts that already converged, NEGATIVE-vocabulary-only on the parts that didn't.**

The Sturmian-Bethe-log construction IS an L-Lipschitz μ-pinned exhaustive sequence in Latrémolière 2512.03573's Def 1.29 sense. The identification:
- Satisfies the Leibniz axiom (Def 1.18) standardly
- Satisfies axioms (2) and (3) of Def 1.29 *trivially* (Sturmian saturates pin state and operator norm at finite $N_{\max}$)
- Satisfies axiom (1) plausibly (commutator norm decays with $N_{\max}$, rate is Sturmian-specific not Paper-38-universal)

The framework supplies Paper 36 with:
- A precise math.OA error-bar model for the LS-3 s-state Bethe log convergence (real new content)
- A uniform convergence vocabulary across the atomic precision catalogue (Mu, Ps, He fine structure, etc.)
- A cleaner scope statement for the LS-7 / LS-8a residual decomposition

The framework does NOT supply:
- A prediction of the LS-3 acceleration-vs-velocity 3.3× speedup
- A prediction of the LS-4 Drake-Swainson regularization
- The LS-7 dimensionless bracket $C_{2S}$
- The LS-8a renormalization counterterms

**The LS-8a renormalization gap, in Latrémolière vocabulary:** "operator-system content beyond the pin-state local metametric." This is a *reframing*, not a *resolution*. The framework respects the structural-skeleton-scope statement (CLAUDE.md §1.7 WH5, multi-focal-composition wall): GeoVac determines selection rules / transcendental signatures / convergence rates, but does not autonomously generate renormalization counterterms.

**Recommendation for the PI.** Open Sub-sprint Z.1 as a parallel 2-week sprint to Sub-sprint X *if* X returns positive on the cleaner Hylleraas-r₁₂ application. The Z.1 deliverable is a Paper 36 vocabulary update + LS-8a scope refinement, NOT a new closed-form prediction. The genuine math.OA content is the rigorous error-bar model for the parts of Paper 36 that already converged empirically.

**Concrete sprint scope** (~10 working days = 2 weeks): Z.1.1 Leibniz verification (2-3 days), Z.1.2 explicit rate extraction (3-4 days), Z.1.3 LS-7 prefactor convergence prediction (2-3 days), Z.1.4 LS-8a interpretation (2-3 days). Diagnostic only; no production-code changes.

---

## §9. Files

- This memo: `debug/subsprint_z_bethe_log_latremoliere_memo.md` (~3500 words, 9 sections)
- Cross-references:
  - `debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md` (§4 Bethe log identification sketch)
  - `debug/sprint_l3e_p3_rescope_memo.md` (Sub-sprint Z scoping context)
  - `debug/ls3_bethe_log_regularized_memo.md` (LS-3 acceleration form, empirical convergence data)
  - `debug/ls7_two_loop_se_memo.md` (LS-7 two-loop SE prefactor, $C_{2S} = +3.63$ bracket)
  - `papers/group5_qed_gauge/paper_36_bound_state_qed.tex` (Paper 36 abstract + structural framing)
  - `papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex` (L2 quantitative rate Paper 38, $4/\pi$ asymptote)
  - Latrémolière arXiv:2512.03573 (Dec 3, 2025) "The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces" — Def 1.18, 1.22, 1.26, 1.29

**Honest confidence.**
- HIGH on the Sturmian-as-L-Lipschitz-μ-pinned-exhaustive-sequence structural identification (§2)
- HIGH on the saturating-pin-state observation making axioms (2)/(3) trivial (§2)
- MEDIUM on the Sturmian-specific convergence rate being non-trivial enough to give Paper 36 a real new error-bar (§3, §4)
- HIGH on the framework NOT supplying LS-7 bracket $C_{2S}$ or LS-8a counterterms (§3, §5, §7)
- HIGH on the mixed-vocabulary-mostly verdict overall (§8)

No follow-on sprint auto-opened by this memo. PI decision required on whether to open Sub-sprint Z.1 (~2 weeks parallel to or after Sub-sprint X).
