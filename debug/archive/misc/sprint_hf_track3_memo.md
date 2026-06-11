# Sprint HF Track 3 — Recoil structural-derivability diagnostic

**Date:** 2026-05-07
**Goal:** Diagnose whether GeoVac's cross-register architecture (Track NI, `geovac/nuclear/nuclear_electronic.py`) natively produces reduced-mass / recoil corrections to the hyperfine Bohr-Fermi prediction, or whether the (1 + m_e/m_p)^{-3} factor that HF-1 applied by hand is an external prescription.
**Verdict:** **NEGATIVE.** Track NI's V_ne is a classical scalar parameter, not a two-body coordinate-coupling operator on the joint register. Recoil at every order remains an external focal-length input.
**Status:** Closed. Recommendation for HF-2 baseline at the bottom.

---

## 1. Why this diagnosis matters

HF-1 produced two numbers:

- **Strict Bohr-Fermi**, $g_e = 2$, no recoil: $A_{\rm hf}^{\rm BF, strict} = 1421.16$ MHz, residual $+0.754$ MHz / $+531$ ppm.
- **Reduced-mass-corrected**, $|\psi(0)|^2 \mapsto |\psi(0)|^2 \cdot (1 + m_e/m_p)^{-3}$: $A_{\rm hf}^{\rm BF, recoil} = 1418.84$ MHz, residual $-1.566$ MHz / $-1102$ ppm.

The recoil step was applied by hand: HF-1's agent multiplied $|\psi_{1s}(0)|^2$ by the textbook factor $(1 + m_e/m_p)^{-3}$ using the standard reduced-mass replacement $m_e \to \mu_{\rm red} = m_e m_p / (m_e + m_p)$ in the Bohr radius. That replacement is universally correct in NRQM/Bethe-Salpeter, but it imports the recoil structure as an *external* prescription rather than deriving it from the framework.

The PI question for HF-3 is sharp: Track NI's architecture *looks like* exactly the multi-focal seam where recoil should be derivable natively. The proton lives on its own register (HO basis, $(n_r, l, m_l, m_s)$), the electron lives on its own register (hydrogenic, $(n, l, m, m_s)$), and the two are tensored together with explicit cross-register couplings. If V_ne were implemented as a two-body Coulomb operator $-Z/|r_e - R_n|$ that quantum-couples electron position to nuclear position, then variational diagonalization on the joint register would automatically pick up the reduced-mass shift in the binding energy and the $(1 + m_e/m_p)^{-3}$ shift in $|\psi(0)|^2$ — without any hand-coding.

That would be a clean structural-skeleton win for multi-focal composition. The diagnosis below shows it isn't there.

## 2. The four checks

I read `geovac/nuclear/nuclear_electronic.py` and `geovac/nuclear/form_factor.py` end-to-end and confirmed each of the following with a runnable check (`debug/sprint_hf_track3.py`).

### Check 1 — V_fs (the only V_ne-like term in Track NI) is a one-body electronic operator

`finite_size_correction(Z, R_nuc, n, l)` (form_factor.py L121) returns a *real number*: $\Delta E_{ns} = (2/5) Z^4 R_{\rm nuc}^2 / n^3$. It is a closed-form perturbation evaluated at the classical scalar `R_nuc = R_PROTON_BOHR ≈ 1.59 × 10^{-5}` bohr. Not a Pauli dict, not an operator on either register.

`finite_size_coupling_pauli` (nuclear_electronic.py L219) consumes that scalar and emits a Pauli operator $V_{fs} = \Delta E \cdot (n_{1s,\uparrow} + n_{1s,\downarrow})$ — three Pauli strings (an identity-on-everything term plus two single-Z's on the electron 1s up/down qubits). Confirmed numerically: 0 of these Pauli strings have any nontrivial Pauli on the nuclear register.

The docstring is explicit and honest:

> "This simplest implementation treats R_nuc as a classical parameter (not a nuclear operator). A fully quantum coupling would replace R_nuc with an operator on the nuclear register (the nuclear charge distribution), but that is deferred as a future extension."

### Check 2 — The hyperfine coupling crosses registers, but only on spin indices

`hyperfine_coupling_pauli` (nuclear_electronic.py L303) builds $A_{\rm hf} \mathbf{I} \cdot \mathbf{S}$ as a four-qubit operator on (proton 0s up, proton 0s down, electron 1s up, electron 1s down). All 12 Pauli strings have weight 2 or 4 and act on exactly those four qubits. `R_nuc` does not appear in the function signature or body. The proton's *spatial* quantum numbers ($n_r, l, m_l$) are pinned to $(0,0,0)$ throughout — only $m_s$ is dynamical, and it couples only to the electron's $m_s$.

This is a spin-spin coupling, not a position-position coupling. It correctly recovers the 21 cm singlet-triplet gap at the level of `Track NI`'s validation test (CLAUDE.md §2 records the hyperfine validation at $3 A_{\rm hf}/4 = 1.62 \times 10^{-7}$ Ha), but it carries zero information about the proton's spatial wavefunction.

### Check 3 — The full Track NI Hamiltonian has zero spatial-spatial cross-register terms

I built the full composed Hamiltonian via `build_deuterium_composed_hamiltonian(N_shells=2, hw=10.0, n_max_elec=2)` and partitioned each Pauli string by which register it touches and whether it touches *spatial* qubits (anything other than the proton 0s spin pair and electron 1s spin pair). Result:

| Bucket | Count |
|---|---|
| Identity | 1 |
| Nuclear-only | 49 |
| Electronic-only | 9 |
| Cross-register, spin-spin only | 12 |
| Cross-register, spatial-spatial | **0** |

The framework has no $V_{ne}(r_e, R_n)$ two-body coordinate operator. There is nothing in the joint-register Hamiltonian that would couple the electron's hydrogenic $(n, l, m)$ to the proton's HO $(n_r, l, m_l)$ — at the level of energy expectation, the two registers are spectator-additive with the sole exception of the spin-spin hyperfine and the (one-body, electronic-only) finite-size correction.

### Check 4 — The reduced-mass factor agrees, but is not derived

The standard reduced-mass replacement gives $|\psi_{1s}(0)|^2 \to (Z \mu_{\rm red}/m_e)^3 / \pi$, so the multiplicative factor on $|\psi(0)|^2$ is $(\mu_{\rm red}/m_e)^3 = (1 + m_e/m_p)^{-3}$. Numerically at CODATA $m_e/m_p = 1/1836.15$: factor = $0.998368$, deviation = $-1632$ ppm. Multiplying HF-1's strict-BF $1421.16$ MHz by this factor reproduces $1418.84$ MHz to displayed precision.

This is consistent with the framework but external to it. The Fock projection $p_0 = \sqrt{-2 E_n}$ in CLAUDE.md §4 is a stereographic scaling of a *single-particle* Hamiltonian; the energy is the electron's binding energy and the implicit mass is $m_e$. Making the Fock projection two-body would require a center-of-mass / relative coordinate split (write the joint Hamiltonian as $\mathbf{P}_{\rm CM}^2/(2M) + \mathbf{p}_{\rm rel}^2/(2\mu) - Z/r$, identify $r$ as the relative coordinate, and apply Fock to the relative-coordinate sector with mass $\mu$ in the conformal scale). Track NI does not host this split. The proton is on its own register but its register holds *occupation numbers* of HO orbitals, not a position operator that enters V_ne.

I want to be careful here, because this is the "alternative positive route" the briefing flagged. The framework *permits* a recoil-aware reading: declare that the conformal scale uses $\mu_{\rm red}$, propagate that through the wavefunction normalization, and the $(1 + m_e/m_p)^{-3}$ factor follows. But that declaration is no different in content from the textbook reduced-mass replacement — it does not buy any additional structural derivation. The factor is a Layer-2 calibration in Paper 34's vocabulary, and the framework permits it but does not produce it from Layer-1 data.

## 3. Verdict

**NEGATIVE.** Track NI's cross-register architecture is set up such that:

- V_ne is a *classical* parameter $R_{\rm PROTON\_BOHR}$ fed into a closed-form perturbation, applied as a one-body electronic operator.
- The proton register holds nuclear-spin and nuclear-occupation states (HO orbitals) for the spin-spin hyperfine coupling; the proton's *spatial coordinate* is not a dynamical variable that couples to the electron's position.
- There is no $V_{ne}(r_e, R_n)$ two-body coordinate operator on the joint register.

The framework cannot natively produce reduced-mass effects, and recoil at every order remains an external focal-length input. This confirms the PI's prior.

What would be required for "positive" (the briefing's machinery (a)–(c)):

- (a) A V_ne term implemented as a two-body coordinate operator: in second quantization, $V_{ne} = -\sum_{ij} \langle e_i, p_j | 1/|r_e - R_n| | e_k, p_l \rangle a_{e_i}^\dagger a_{p_j}^\dagger a_{p_l} a_{e_k}$ on the joint register.
- (b) Variational diagonalization of $H_{\rm joint} = T_e + T_p + V_{ne}$ in a basis spanning both registers, with $T_p$ scaling as $1/m_p$.
- (c) Verification that $m_p \to \infty$ recovers the standard hydrogen 1s and finite $m_p$ shifts $|\langle r_e = 0\rangle|^2$ by $(1 + m_e/m_p)^{-3}$.

Per the briefing's instructions, I do not propose extending Track NI to add (a) — that's a multi-week separate sprint. Track NI's architectural choice (proton occupation numbers + classical R_nuc) is the right level of fidelity for what Track NI was built to demonstrate (a 26-qubit composed nuclear-electronic POC at the scale of MeV-Ha-Hz). It's just not the level needed to derive recoil natively.

## 4. Updated A_hf prediction status

The HF-1 number set stands. For HF-2's baseline, **the right target is $A_{\rm hf}^{\rm BF, recoil} = 1418.84$ MHz with the recoil factor applied externally**, residual $-1.566$ MHz / $-1102$ ppm to attack. The recoil correction stays as a hand-applied factor — Track NI cannot produce it. HF-2's job is to derive $a_e$ from graph-native vertex correction and reduce the residual to $\sim 33$ ppm; that path uses the GeoVac vertex machinery (Paper 28 §curvature_coefficients, the $F_2(\kappa)$ pipeline), which is genuinely framework-internal and a fair test of the structural-skeleton claim.

A defensible alternative is to flag both numbers in HF-2 and report residuals against both, with a note that the framework hosts the recoil prescription but does not produce it. That keeps the framework's scope honest without burying the strict-BF baseline.

## 5. PI-facing framing

This result confirms that the Track NI cross-register architecture is "single-focal calibration plus parametric coupling," not a multi-focal composition that resolves recoil natively. The architecture is exactly what its name says — a composed *qubit Hamiltonian* whose nuclear and electronic blocks are each Hamiltonians at their natural focal length, glued by spin-spin (hyperfine) and parametric (V_fs at fixed R_nuc) cross-terms. The MeV-Ha 10^13 dynamic range that makes Track NI a useful PoC for quantum simulation is the *same* dynamic range that prevents single-pass variational diagonalization from coupling the registers' spatial degrees of freedom. The framework is therefore consistent with what we already understood from the structural-skeleton scope memo (`geovac_structural_skeleton_scope_pattern.md`): GeoVac determines the structural skeleton (selection rules, transcendental class, scaling laws) but recoil — like Yukawa selection in H1, like renormalization counterterms in LS-8a, like the proton g-factor — is a calibration input.

The honest one-line: GeoVac's cross-register architecture provides a *vehicle* for two-body recoil (you can put the right operators on it), but the operators have not been built. Building them is a separate sprint with its own scope; using the textbook reduced-mass replacement is the appropriate calibration in the meantime.

A small concrete suggestion for HF-2: when the agent reports the residual after deriving $a_e$, also report the recoil correction explicitly as an external Layer-2 input in the Paper 34 vocabulary, not as a piece of the framework's prediction. That keeps the structural claim crisp.

## 6. Files

- `debug/sprint_hf_track3.py` — diagnostic script (no production code modified).
- `debug/data/sprint_hf_track3.json` — structured outputs of the four checks.
- `debug/sprint_hf_track3_memo.md` — this memo.
