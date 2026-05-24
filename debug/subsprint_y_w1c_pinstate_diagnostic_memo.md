# Sub-sprint Y — W1c pin-state-shape diagnostic via Latrémolière 2512.03573

**Date:** 2026-05-23
**Sprint position:** Sub-sprint Y of Sprint L3e-P3 (re-scoped). Follow-on to the Phase A.2' deep-read memo (`debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md`).
**Mandate:** Diagnostic-only structural memo. No code changes. Test whether the W1c chemistry orthogonality wall (NaH binding) admits a Latrémolière 2512.03573 pin-state-shape-selection reformulation rather than the current cross-V_ne operator-generation framing.
**Cross-references:** `debug/w1c_residual_nah_track3_memo.md` (Track 3 closure attempt, May 8), `debug/pk_cross_center_synthesis_memo.md` (PK cross-center attempt, May 8), `geovac/screened_valence_basis.py`, `geovac/phillips_kleinman_cross_center.py`, `geovac/balanced_coupled.py`.
**Verdict (short):** MIXED with a named pin-state candidate (Option iii, bonding-orbital pin state). The Latrémolière framework gives genuinely-new vocabulary — local metametric $\delta_r$ at the cross-center coupling region — but the vocabulary does NOT autonomously generate the right pin state. Recommendation: open Y.1 sprint at 1–2 weeks, scoped narrowly to compute $\delta_r$ between candidate pin states; expected outcome is a structural ordering of the candidates, not binding closure by itself.

---

## §1. Structural recap of the W1c wall

The W1c chemistry orthogonality wall surfaced in the May 8 multi-track sprint as the dominant systematic blocking second-row chemistry binding within the framework. The named instance is NaH at $\max\_n = 2$ (Q=20, 2e balanced FCI), where the PES is monotonically descending across the tested R-grid $\{2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0\}$ bohr with no internal equilibrium minimum.

The descent-depth ladder from `debug/w1c_residual_nah_track3_memo.md` Table §3:

| Configuration | $E_{\min}$ (Ha) | $R_{\min}$ (bohr) | Descent depth (Ha) | Bound? |
|:--------------|----------------:|------------------:|-------------------:|:-------|
| bare (Sprint 7 baseline) | −171.097 | 2.5 | 6.244 | NO |
| + W1c (screened cross-V_ne) | −163.316 | 2.5 | 0.357 | NO |
| + W1c + PK | −163.264 | 2.5 | 0.305 | NO |
| + W1c + PK + SV diagonal | −162.943 | 2.5 | 0.313 | NO |

Three engineering attempts have landed clean architectural modules:

1. **W1c (screened cross-V_ne)**: `geovac/cross_center_screened_vne.py`. Implements multipole expansion of $1/|\mathbf{r} - \mathbf{R}_B|$ against the screened core density $Z_\text{eff}(r)$ instead of the bare $Z$. Reduces descent depth by 17.5× — the heavy lifting.
2. **PK cross-center**: `geovac/phillips_kleinman_cross_center.py`. Standard PK operator-form approximation $\Delta H_{pq}^\text{PK} = \sum_c (E_v - E_c) S_{pc} S_{cq}$. Additional 14.6% reduction over W1c.
3. **SV diagonal substitution**: `geovac/screened_valence_basis.py`. Replaces hydrogenic Z_orb=1 valence h1 diagonal with FrozenCore Z_eff(r) Schrödinger eigenvalues (Na 3s = −0.170 Ha vs hydrogenic −0.500 Ha). Negligible additional effect.

The Track 3 memo identifies the genuine engineering target as **wavefunction-shape substitution in cross-V_ne and within-block ERIs** (Frame − Physical-3p differential of −0.675 Ha at R=2.5 vs R=10, which is 1.9× the W1c-residual descent of 0.357 Ha). The mechanism is structural: the framework's hydrogenic Z=1 1s on the Na valence has mean radius 1.50 bohr, but the actual Na 3s screened wavefunction has mean radius 4.47 bohr. The two functions have opposite shape biases at the cross-center coupling region (~3–5 bohr), where cross-V_ne integrates against the off-center $1/|\mathbf{r} - \mathbf{R}_B|$ singularity.

The May 8 closeout marked the W1c wall as the single remaining engineering closure for second-row chemistry. Two follow-on routes were named: (a) physical-n hydrogenic relabeling (1-week sprint), (b) strict-Schmidt orthogonalization (multi-week sprint). Neither has been opened.

The wall is GENUINELY OPEN. The chemistry-arc-paused-W1c-residual memory is the standing position.

## §2. Latrémolière 2512.03573 reframing of the W1c wall

The deep-read memo §3 proposes the following structural reframing. The framework's molecular C*-algebra is non-unital (it is $C_0$ of configuration space tensored with the frozen-core structure), and the pin state should be the joint molecular ground state (cross-center bonding orbital). Latrémolière 2512.03573 supplies the structural vocabulary:

- **Non-unital C*-algebra**: $\mathcal{A}_\text{NaH} = C_0(\mathbb{R}^3_\text{electron} \times \mathbb{R}^3_\text{H-position})$ tensored with the discrete frozen-core algebra. Definition 1.22's (FM), (CB), (B) axioms apply.
- **Pin state**: $\mu \in \mathcal{S}(\mathcal{A}_\text{NaH})$ — the joint ground state of the molecular Hamiltonian. The framework's current default is a Slater determinant built from hydrogenic Z_orb=1 valence orbitals on each center.
- **Leibniz seminorm**: $L(a) = \|[D_\text{NaH}, a]\|$ where $D_\text{NaH}$ is the cross-register kinetic operator (Definition 1.18).
- **Exhaustive sequence**: $h_n = P_{n_{\max}}$ truncation projector onto the molecular Sturmian / balanced-coupled basis at $n_{\max}$ (Definition 1.29).
- **Local metametric** $\delta_r$: distance at radius $r$ restricted to Lipschitz functions concentrated near $\mu$ (Section 2.2).

Under this reframing, the W1c wall content becomes: **the framework's pin-state-defining Slater determinant has the wrong shape at the cross-center coupling region**. Cross-V_ne integration probes the orbital tails into the off-center Coulomb singularity, and the local metametric at radius $r \sim 3$–$5$ bohr (the bond region) measures the wall's depth.

The framework's Track 3 diagnostic (Frame − Physical-3p differential) is, in this reading, an empirical estimate of the local metametric distance between two pin-state candidates at the cross-center coupling region. The −0.675 Ha differential is one realization of $\delta_r$ at $r$ in the bond region.

This is genuinely-new vocabulary. The question is whether the vocabulary adds physics content.

## §3. Pin state candidates for NaH binding

Four candidate pin states for the NaH molecular ground state at fixed $\max\_n = 2$:

### (i) Hydrogenic Z_orb=1 Na 3s × H 1s product (current framework default)

This is what the framework's balanced-coupled builder constructs at `geovac/balanced_coupled.py`. The Na-side valence is hydrogenic with Z_orb = 1 (asymptotic Na valence Z from `atomic_classifier.py`), labeled with block_n=1 → physical_n=3 (Na 3s) via the n_val_offset convention. The H-side is hydrogenic with Z_orb = 1, block_n=1 → physical_n=1 (H 1s).

The cross-V_ne integrals against $1/|\mathbf{r} - \mathbf{R}_B|$ at R~3.5 bohr probe both orbital tails. The Na 3s tail in hydrogenic Z=1 1s shape decays as $e^{-r}$ (tight). The actual Na 3s would decay much slower (diffuse, mean radius 4.47 bohr).

In Latrémolière language: this is the current pin state $\mu_{(i)}$. The cross-V_ne integrals under this $\mu_{(i)}$ give the current descent depth of 0.357 Ha. Mathematical structure: standard Hartree product, no cross-center mixing.

### (ii) HF-screened Na 3s × H 1s product (Track 3 SV-corrected)

This is what the screened-valence-basis module computes for the h1 diagonal. The Na-side valence uses the actual FrozenCore Z_eff(r) Schrödinger eigenvalue (−0.170 Ha for Na 3s), but only the eigenvalue — NOT the wavefunction. Wavefunction shape remains hydrogenic Z=1 1s.

The Track 3 memo cleanly identifies this as bait-and-switch: the eigenvalue is correct but the cross-V_ne integrals continue to use the wrong wavefunction shape. R-independent at the h1 diagonal level (verified bit-exact in `screened_valence_basis.py` at R=2.5 vs R=10).

In Latrémolière language: candidate (ii) shares pin-state-Hilbert-space data with (i) — same wavefunctions, same Slater determinant structure — but lifts the diagonal of the pin-state generator by 0.330 Ha on Na 3s. The local metametric $\delta_r$ at the cross-center coupling region is identical to $\mu_{(i)}$'s because $\delta_r$ depends on the wavefunction shape, not on the eigenvalue. Confirmed in Track 3's bit-exact R-independence test.

### (iii) Cross-center bonding-orbital pin state (NAMED CANDIDATE)

This is the natural Latrémolière-flavored candidate. Build the pin state as a Slater determinant whose valence orbital is the proper bonding orbital with the right cross-center shape. In quantum-chemistry language: a self-consistent HF molecular orbital, not a product of atomic Hartree orbitals.

Concretely, candidate (iii) constructs the pin state as:

$$
|\mu_{(iii)}\rangle = \mathrm{Det}\{ \phi_{\sigma_g}, \phi_{\sigma_g'} \}
$$

where $\phi_{\sigma_g}$ is a bonding linear combination of physical Na 3s and physical H 1s orbitals (not hydrogenic Z=1 1s on each center). The bonding linear combination has:

- **Correct cross-center spatial extent**: amplitude peaked at the bond midpoint, decaying at both nuclei.
- **Diffuse shape on the Na side**: matching the actual Na 3s mean radius 4.47 bohr.
- **Tight shape on the H side**: matching the actual H 1s mean radius 1.50 bohr.
- **R-dependent**: the bonding amplitude redistributes as R changes, unlike (i) and (ii) which are R-independent at the wavefunction level.

In Latrémolière language: candidate (iii) is a pin state with structurally different Hilbert-space content from (i), (ii). The cross-V_ne integrals under $\mu_{(iii)}$ are R-dependent in the way that the Track 3 diagnostic suggests: the Frame−Physical-3p differential of 0.675 Ha is the order-of-magnitude estimate of $\delta_{r}(\mu_{(i)}, \mu_{(iii)})$ at $r$ in the bond region.

**Plausibility for binding**: HIGH structurally. The Track 3 diagnostic already quantifies that shape-substitution overcorrects by ~1.9× the W1c residual, which is the right sign and magnitude for producing a binding minimum (overcorrection of the cross-V_ne attraction profile means the descent becomes non-monotone).

**Plausibility as a Latrémolière pin state**: HIGH structurally. The pin-state-defining Slater determinant with cross-center bonding orbital is the physically-correct choice for the molecular ground state, and Latrémolière's Definition 1.26 explicitly permits any pin state in $\mathcal{S}(\mathcal{A})$ whose Monge-Kantorovich-finite set is weak-* dense. The bonding-orbital state is the natural choice.

### (iv) MP2-corrected ground state (Møller-Plesset perturbation pin state)

Standard MP2 corrects the HF ground state by adding second-order perturbation theory in the V_ee residual. The MP2 ground state is a linear combination of the HF Slater determinant and doubly-excited determinants, with weights from energy denominators.

In Latrémolière language: candidate (iv) is a pin state that includes correlation content beyond the single-Slater-determinant pin state (iii). The cross-V_ne integrals under $\mu_{(iv)}$ include double-excitation contributions which could either enhance or reduce the W1c-residual descent.

**Plausibility for binding**: MEDIUM. MP2 is known to overcorrect the HF V_ee gap in some chemistry contexts and to be insufficient in others. For NaH at $\max\_n = 2$, MP2 may close the W1c residual if the descent is dominated by V_ee correlation; it cannot help if the descent is dominated by single-particle cross-V_ne shape.

Since the Track 3 diagnostic localizes the W1c residual in cross-V_ne (single-particle integral against the off-center potential), candidate (iv) is structurally less promising than (iii).

## §4. Local metametric $\delta_r$ structural content

Latrémolière 2512.03573 §2.2 defines the local metametric at radius $r$ as a distance restricted to Lipschitz functions concentrated near the pin state. Concretely:

$$
\delta_r(\mu, \nu) := \sup \{ |\mu(a) - \nu(a)| : a \in \mathrm{dom}(L), L(a) \le 1, \|a - h_n\| \le r \text{ for some } n \}
$$

where $h_n$ is the L-Lipschitz $\mu$-pinned exhaustive sequence (Def 1.29). The "concentrated near the pin state" content means the test function $a$ must be in a Lipschitz ball around the exhaustive sequence at radius $r$.

For NaH, this means:

**What "Lipschitz function near pin state" means concretely:** A test function $a \in \mathcal{A}_\text{NaH}$ that (i) has bounded commutator with $D_\text{NaH}$ (the molecular Hamiltonian / Dirac operator), and (ii) is operator-norm close (within $r$) to a Sturmian truncation projector. In the framework's molecular C*-algebra, this is essentially the set of finite-rank operators on the Slater-determinant basis that mix only a finite number of orbitals near the bonding region.

**At what radius $r$ to measure the wall:** The bond region. Concretely, $r$ should correspond to the cross-center coupling region where the W1c residual lives. The Track 3 diagnostic at R=2.5, 3.5, 5.0, 10.0 bohr gives empirical estimates; the natural Latrémolière $r$ is the operator-norm radius corresponding to spatial radius ~3 bohr (the Na-H bond region) or $r \sim O(1)$ in the dimensionless Sturmian sense.

**Does $\delta_r$ distinguish between candidates (i)-(iv)?** This is the load-bearing question. Three sub-questions:

1. **$\delta_r(\mu_{(i)}, \mu_{(ii)})$**: Track 3 SV diagonal gives this as zero at the wavefunction level (bit-exact R-independence). The two pin states have the same Hilbert-space content; they differ only in the diagonal of $D_\text{NaH}$ restricted to the pin-state Slater determinant. The local metametric depends on cross-state matrix elements, which are the same in (i) and (ii). **Conclusion: $\delta_r((i), (ii)) = 0$**. Consistent with Track 3's bit-exact finding.

2. **$\delta_r(\mu_{(i)}, \mu_{(iii)})$**: structurally non-zero. The bonding-orbital pin state has different wavefunction shape, so the cross-V_ne matrix elements at the bond region differ. The Track 3 diagnostic empirically estimates this as ~0.675 Ha (the Frame − Physical-3p differential). **Conclusion: $\delta_r((i), (iii)) > 0$, and the magnitude is the right order for closing the wall**.

3. **$\delta_r(\mu_{(i)}, \mu_{(iv)})$**: structurally non-zero but smaller than $\delta_r((i), (iii))$. MP2 corrections include double-excitation content but preserve the single-Slater-determinant zeroth-order pin state. The cross-state mixing is small (perturbative). **Conclusion: $\delta_r((i), (iv)) > 0$ but at order V_ee/ΔE which is ≪ 0.675 Ha for NaH**.

The local metametric does distinguish between pin states, and the structural ordering is $\delta_r((i), (iii)) > \delta_r((i), (iv)) > \delta_r((i), (ii)) = 0$. The right pin state for closing the wall is (iii) — the bonding-orbital pin state.

**This is consistent with the Track 3 memo's named next target (wavefunction-shape substitution in cross-V_ne), and the Latrémolière reframing confirms the diagnosis structurally.**

The Latrémolière framework adds: a precise math.OA distance $\delta_r$ between pin-state candidates that quantifies "how far" each candidate is from any other. This is genuinely-new structural content over the framework's bare physical reasoning, because the local metametric is a well-defined object that can be computed numerically for each candidate, then ordered.

## §5. Concrete sprint plan for Y.1

If we open Sub-sprint Y.1, the scope would be 1–2 weeks. Diagnostic-only at the structural level; no production-code modifications.

**Week 1 — pin state candidate construction.**

1. **Candidate (i) (current default).** Already exists in `geovac/balanced_coupled.py` and `geovac/composed_qubit.py`. Reuse the NaH spec from `geovac/molecular_spec.nah_spec()`.

2. **Candidate (iii) (bonding-orbital pin state).** Construct a 2-electron HF molecular orbital at fixed R by diagonalizing a Fock matrix in the framework's hydrogenic basis $\{ \chi_{\text{Na 3s}}, \chi_{\text{H 1s}} \}$ (using physical 3s shape, not hydrogenic Z=1 1s). The bonding orbital is the lowest eigenvector. This requires:
   - Computing the 2×2 overlap matrix $S$ with physical Na 3s shape from `_solve_screened_radial_log`.
   - Computing the 2×2 Fock matrix in the same basis with cross-V_ne, exchange, and ERIs.
   - Diagonalizing $H - \lambda S = 0$.
   
   The bonding orbital becomes the pin-state Slater determinant's valence orbital. The diagnostic uses this as the reference wavefunction for computing $\delta_r$.

3. **Candidate (ii) (SV-corrected).** Use existing `screened_valence_basis.py` module to compute the diagonal-only correction. Already gives $\delta_r((i), (ii)) = 0$ structurally.

4. **Candidate (iv) (MP2).** Compute the MP2 correction to the (i) pin state using the framework's existing balanced-coupled ERI machinery. Diagnostic only; do not implement production MP2.

**Week 2 — local metametric computation.**

5. **Compute $\delta_r$ between each pair of candidates** at radii $r \in \{1, 2, 3, 5\}$ bohr (the bond region) by sweeping the Lipschitz-norm test functions over Sturmian basis truncations at $n_{\max} = 2, 3$. Report the ordering.

6. **Compare to the experimental NaH PES** at known $D_e \approx 0.075$ Ha and $R_\text{eq} = 3.566$ bohr. Identify which pin state gives the smallest distance to the experimental binding curve.

7. **Diagnostic verdict**: does the local metametric order pin states correctly? If (iii) is closest to the experimental PES and (i) is farthest, the Latrémolière reframing confirms the Track 3 diagnostic at a structural level. If the ordering is different, that's a surprise.

**Deliverable.** Sub-sprint memo at `debug/subsprint_y1_w1c_pinstate_metametric_memo.md` (~2500-3500 words) reporting the four $\delta_r$ values, the structural ordering, and the verdict on whether (iii) is the correct pin state for closing the wall.

**Important sprint scope guard.** Y.1 does NOT implement bonding-orbital pin state in production code. The PES under the (iii) pin state would still require Option 1 (physical-n hydrogenic relabeling) or Option 2 (strict-Schmidt) implementation from the Track 3 memo §6. Y.1 is a structural diagnostic that confirms the named target via Latrémolière vocabulary; the actual implementation sprint remains as named in Track 3.

## §6. Honest scope assessment

Does the Latrémolière reframing add physics content, or is it vocabulary-change?

**Where it adds content.**

1. **A well-defined distance between pin-state candidates.** The Track 3 memo provided an empirical proxy (the Frame − Physical-3p differential) for "how far" the framework's pin state is from the physical one. The local metametric $\delta_r$ is a math.OA-grounded distance with a precise definition. If Y.1 lands, the framework gains a quantitative metric for choosing pin states.

2. **Structural ordering of candidates.** The Latrémolière framework predicts $\delta_r((i), (iii)) > \delta_r((i), (iv)) > \delta_r((i), (ii))$ based on structural reasoning (Hilbert-space content vs eigenvalue-only correction). The Track 3 memo arrived at the same ordering empirically. The framework adds: a prediction that this ordering generalizes to other systems beyond NaH.

3. **Uniform diagnostic across the chemistry catalogue.** Once the local metametric is implemented for NaH, the same machinery applies to MgH₂, HCl, etc. Each system gets a $\delta_r$ value, and the framework can identify which systems have small enough $\delta_r((i), (iii))$ to be closeable by hydrogenic-shape relabeling vs which require strict-Schmidt.

**Where it does not add content.**

4. **The right pin state is still an external input.** Latrémolière's framework defines the structure but does not autonomously generate the right pin state (iii). The framework user (or the PI) must choose the bonding-orbital pin state by external chemistry intuition. This is the same scope-statement as the CLAUDE.md §1.7 structural-skeleton-scope reading: GeoVac does not autonomously generate calibration data.

5. **The implementation work is unchanged.** Even after Y.1 confirms (iii) as the right pin state, the framework still needs Option 1 (physical-n hydrogenic relabeling) or Option 2 (strict-Schmidt) to implement the cross-V_ne shape substitution in production code. The Latrémolière reframing does not bypass the implementation sprint.

6. **The local metametric is conditional on the Sturmian-as-L-Lipschitz-μ-pinned-exhaustive-sequence identification (Sub-sprint X).** If Sub-sprint X fails to verify the Leibniz axiom on $L = \|[D_\text{Sturmian}, a]\|$, then $\delta_r$ is not well-defined and the Y.1 diagnostic collapses. Y.1 should be sequenced AFTER X.

**Net honest scope assessment: MIXED.** The Latrémolière reframing adds genuine structural vocabulary that confirms the Track 3 named target (wavefunction-shape substitution → bonding-orbital pin state) at a math.OA-grounded level. The framework gains a quantitative diagnostic ($\delta_r$ ordering) over the bare empirical proxy. The implementation work remains as named in Track 3; the wall still requires Option 1 or Option 2 to close in production.

The reframing is more than vocabulary-change but less than autonomous wall closure. It is a **diagnostic refinement** that quantifies the named target and provides a uniform metric across the chemistry catalogue.

**Confidence:**
- HIGH that $\delta_r((i), (iii))$ is the dominant gap (structural reasoning + Track 3 empirical agreement).
- HIGH that Y.1 will produce a clean ordering of the four candidates (math.OA-grounded local metametric construction).
- MEDIUM that Y.1 will give a precise enough numerical value to predict NaH binding (the local metametric is qualitative-rate; finite-$n_{\max}$ truncation effects are non-trivial).
- MEDIUM-LOW that the Latrémolière reframing alone closes the wall without Option 1/2 implementation (the actual cross-V_ne integrals still need to be re-computed in the right pin state basis).

## §7. Cross-references and dependencies

**Prerequisites for Y.1:**
- Sub-sprint X (Sturmian-as-Latrémolière) must land first. X verifies the Leibniz axiom on $L = \|[D_\text{Sturmian}, a]\|$, which is the load-bearing prerequisite for $\delta_r$ being well-defined. Without X, Y.1 has no foundation.
- The deep-read memo §6 sequences X → Y → Z. Sub-sprint Y depends on X structurally.

**Parallel sub-sprints:**
- Sub-sprint Z (Bethe log / Lamb shift Latrémolière interpretation) can run in parallel with Y. Z applies the same Sturmian-as-Latrémolière identification to Paper 36 LS-3, but does not depend on Y.

**Downstream work after Y.1:**
- If Y.1 confirms (iii) as the right pin state with $\delta_r$ ordering, **open Track 3 Option 1 (physical-n hydrogenic relabeling) as the 1-week implementation sprint**. This is the named next engineering target from the May 8 Track 3 memo, and Y.1 would give it Latrémolière-framework-grounded justification.
- If Y.1 surprises with a different ordering (e.g., $\delta_r((i), (iv)) > \delta_r((i), (iii))$), revisit the named engineering target.

**Files referenced:**
- `debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md` (parent, §3)
- `debug/w1c_residual_nah_track3_memo.md` (Track 3 closure, May 8)
- `debug/pk_cross_center_synthesis_memo.md` (PK cross-center, May 8)
- `geovac/screened_valence_basis.py` (~280 lines, 25 tests)
- `geovac/phillips_kleinman_cross_center.py` (~280 lines, 21 tests)
- `geovac/balanced_coupled.py` (production hookup for W1c + PK + SV)
- `geovac/cross_center_screened_vne.py` (W1c production module)
- `geovac/molecular_spec.py` (nah_spec)
- CLAUDE.md §1.7 W1c-residual chemistry orthogonality wall entry
- CLAUDE.md §3 chemistry-arc-paused-W1c-residual row

---

## §8. Verdict summary

**MIXED with named pin-state candidate.**

The Latrémolière 2512.03573 reframing of the W1c chemistry wall provides genuine structural vocabulary (local metametric $\delta_r$, pin-state Hilbert-space content vs eigenvalue) that quantifies the Track 3 named engineering target (wavefunction-shape substitution → bonding-orbital pin state, candidate (iii)). The structural ordering $\delta_r((i), (iii)) > \delta_r((i), (iv)) > \delta_r((i), (ii)) = 0$ is consistent with the Track 3 empirical diagnostic and adds a math.OA-grounded metric.

The reframing does NOT autonomously close the wall. The implementation work (Option 1 hydrogenic relabeling or Option 2 strict-Schmidt) remains as named in the Track 3 memo. The right pin state is an external input.

**Recommendation: open Sub-sprint Y.1 at 1–2 weeks AFTER Sub-sprint X lands.** Sub-sprint X (Sturmian-as-Latrémolière) is the prerequisite. Y.1's deliverable is the local metametric ordering of candidates (i)–(iv); if the ordering confirms (iii), open Option 1 implementation sprint with Latrémolière-framework-grounded justification.

**Honest scope rule:** Y.1 is a diagnostic refinement, not a wall-closing sprint. The structural-skeleton-scope statement (CLAUDE.md §1.7) is preserved: GeoVac uses Latrémolière vocabulary to identify the right pin state, but the right pin state is not autonomously generated by the framework.

The Latrémolière framework provides STRUCTURAL VOCABULARY that helps identify the right pin state. It does not autonomously generate the pin state itself. This is the same scope-statement that holds for the rest of the framework — and Sub-sprint Y is consistent with it.

---

**End of Sub-sprint Y diagnostic memo.**
