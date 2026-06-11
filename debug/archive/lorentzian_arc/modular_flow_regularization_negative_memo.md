# Modular flow on the Krein wedge as a regularization scheme — NEGATIVE

**Date:** 2026-05-18
**Context:** Post-Paper 45 conversational sprint asking whether the new Lorentzian machinery (Krein space, hemispheric wedge KMS state, Tomita-Takesaki modular flow) moves the LS-8a two-loop QED renormalization gap.

## Question

LS-8a verdict (2026-05-07): iterated CC spectral action on Dirac-S³ reproduces the UV-divergent integrand of two-loop QED faithfully — right prefactor $(\alpha/\pi)^2(Z\alpha)^4/n^3$, right sign, divergence $\sim N^{3.43}$ — but cannot autonomously generate wavefunction $Z_2$ or mass $\delta m$ counterterms to extract a finite finite-part coefficient.

Sprint L1 + L2-E (Papers 42, 43) built modular Hamiltonian machinery on the Krein hemispheric wedge with bit-exact $\sigma_{2\pi}(O) = O$ closure and integer-spectrum $K_\alpha^W$. Speculative hypothesis: could the modular flow provide a renormalization scheme that generates $Z_2 / \delta m$ counterterms structurally?

Two parallel agents dispatched 2026-05-18: literature survey + GeoVac-internal code feasibility.

## Verdict: NEGATIVE on both sides, converging

### Literature side (no published precedent + structural obstruction)

No paper in the published literature uses modular flow $\sigma_t = \mathrm{Ad}(\Delta^{it})$ to replace bare cutoff regularization or to generate $Z_2 / \delta m$ counterterms in perturbative QFT. The thermal time hypothesis (Connes-Rovelli 1994, arXiv:gr-qc/9406019) is interpretive — it identifies modular flow with time at the algebraic level but does not deploy $\Delta^{it}$ as a UV regulator.

The closest adjacent program is Fredenhagen-Lindner pAQFT with KMS states (CMP 332, 895; arXiv:1306.6519). This program does combine modular structure and perturbative renormalization, but the arrow runs the *opposite* direction: standard Epstein-Glaser renormalization first (which produces $Z_2$, $\delta m$), then exhibit KMS / modular structure of the resulting interacting state.

**Published structural obstruction** (load-bearing): Bostelmann-Cadamuro-Minz 2025, "Non-local modular flows across deformed half-spaces" (arXiv:2501.02998). Modular flow is generically *non-local* for interacting theories. $\Delta^{it}$ preserves the local algebra of the wedge but does not implement scale transformations in the Wilsonian sense — there is no published wedge-to-momentum-shell correspondence. Modular flow is conjugate to KMS-temperature, not to RG-scale.

Other adjacent programs (Casini-Huerta entropic c-theorems, Faulkner et al. modular toolkit in holography, Hollands-Wald QFT-in-CGS) all use modular Hamiltonians as *probes* of RG fixed points or as bulk-reconstruction tools, never as regularization machinery generating counterterms.

### Code side (three independent structural obstructions)

1. **Wrong variable.** $K_\alpha^W$'s spectrum is $m_j$ (angular projection on the wedge); the LS-8a divergence is in internal-line level $n$. Damping by $e^{-K_\alpha^W/\beta}$ suppresses high-$m_j$ components of the external state, not the internal sum over $(n_1, n_2, n_3, q_a, q_b)$.

2. **Wrong action space.** Modular flow $\sigma_t$ acts on operators in $O_{n_{\max}}$ (the truncated operator system at the algebra-action level). $\Sigma_{2L}$ in LS-8a is a c-number bound-state matrix element, not an operator the modular flow can act on. There is no operator in $O_{n_{\max}}$ whose modular evolution generates the renormalization counterterms.

3. **Wedge restriction is orthogonal to the UV divergence.** The hemispheric wedge $P_W = (1/2)(I + R_{\text{polar}})$ cuts the spinor basis by $m_j$-parity. LS-8a's divergence is in the bulk — every internal mode contributes regardless of $m_j$ sign. Restricting to the wedge halves $\dim \mathcal{H}$ but leaves the $n^{3.43}$ growth rate untouched.

The $\sigma_{2\pi}(O) = O$ closure that powers Sprint L1 is a *kinematic identity* of the rotation generator (integer spectrum gives $e^{i \cdot 2\pi \cdot n} = 1$), not a regularization mechanism. No adjustable scale parameter.

## Synthesis

Two independent angles converge on the same negative. The literature has a published structural obstruction (modular flow non-locality for interacting theories, Bostelmann-Cadamuro-Minz 2025). The code has three independent GeoVac-internal obstructions (wrong variable, wrong action space, orthogonal to UV divergence).

**The LS-8a renormalization wall is established as genuinely signature-independent.** The Lorentzian machinery (Papers 42, 43, 44, 45) does not move it. The structural-skeleton-scope framing of GeoVac is intact and tightened at the renormalization-counterterm boundary: counterterm generation is not autonomous, and is not unlocked by signature extension.

## Recommendation

1. **Do not pursue** modular-flow-as-regulator as an active sprint direction. The negative is clean and converged from two independent angles.
2. **Citation.** Bostelmann-Cadamuro-Minz 2025 (arXiv:2501.02998) is the published obstruction. Cite when the question recurs.
3. **CLAUDE.md §3 dead ends:** add a row noting "modular flow as QFT regularization scheme — NEGATIVE, both literature (non-locality obstruction) and code (wrong variable, wrong action space)." PI to apply if desired.
4. **LS-8a-renorm extension** (multi-loop QED counterterm generation in a way the framework natively supports) remains an open named follow-on, but the modular-flow path is now ruled out as a route to it.

## Files

- This memo: `debug/modular_flow_regularization_negative_memo.md`
- Original LS-8a memo: `debug/ls8a_two_loop_self_energy_memo.md`
- Two-loop self-energy infrastructure: `geovac/two_loop_self_energy.py`
- Modular Hamiltonian infrastructure: `geovac/modular_hamiltonian.py`, `geovac/modular_hamiltonian_lorentzian.py`
