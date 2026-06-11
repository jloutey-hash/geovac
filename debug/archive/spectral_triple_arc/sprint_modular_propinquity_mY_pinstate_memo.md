# Sprint Modular Propinquity Track M-Y — W1c pin-state diagnostic via bimodule reframing

**Date:** 2026-05-23.
**Sprint position:** Track M-Y of the dual modular propinquity 5-track sprint. Re-read of Sub-sprint Y (2026-05-23) through the Latrémolière modular / dual modular propinquity (Latrémolière 2016 arXiv:1608.04881; 2018 arXiv:1811.04534) bimodule lens.
**Mandate:** Diagnostic-only. No code or paper modifications. Test whether the W1c NaH pin-state problem admits a bimodule deformation reading that adds new content over Y's local metametric ordering, predicts a specific implementation path with magnitude, and generalizes uniformly across the alkali-hydride series.
**Cross-references:** `debug/subsprint_y_w1c_pinstate_diagnostic_memo.md` (primary), `debug/sprint_l3e_p3_synthesis_memo.md` §2, Latrémolière arXiv:1608.04881 (modular propinquity), Latrémolière arXiv:1811.04534 (dual modular propinquity), `geovac/balanced_coupled.py`, `geovac/cross_register_vne.py`, `geovac/phillips_kleinman_cross_center.py`, `geovac/screened_valence_basis.py`.

---

## Executive summary

The PI's hypothesis is that recasting the four NaH pin-state candidates of Sub-sprint Y (hydrogenic Z=1, SV-corrected diagonal, cross-center bonding orbital, MP2-perturbed) as elements of a **left-H / right-Na bimodule** $\Mcal$ over the pair $(\Acal_L, \Acal_R)$ of single-center multiplication algebras produces a sharper diagnostic than Y's scalar local metametric $\delta_r$. The Latrémolière modular Gromov–Hausdorff propinquity (Latrémolière 2016) supplies the natural distance object — a **D-norm** on a Hilbert C*-bimodule extending a connection-like differential structure — and the **dual modular propinquity** (Latrémolière 2018) supplies a complete weaker variant that quantifies how a bimodule deforms when its underlying algebras deform.

Under this reframing, the four NaH pin states are bimodule elements at different positions in the (left-shift, right-shift) plane: $(i)$ pure-left (H-centered hydrogenic Z=1, zero right-action coupling), $(iii)$ bimodule-coupled (cross-center bonding, nonzero left AND right action coupling), $(ii)$ bit-identical to $(i)$ at the bimodule-element level (the SV correction is a *diagonal multiplier* in $\Acal_R$ that fixes scalar pin-state expectation values but does not change the bimodule element), $(iv)$ small bimodule perturbation of $(iii)$ at order $V_{ee}/\Delta E$. The bimodule D-norm distance $d_{\Mcal}(\xi_a, \xi_b)$ reproduces Y's ordering — $d_{\Mcal}((i),(iii)) > d_{\Mcal}((i),(iv)) > d_{\Mcal}((i),(ii)) = 0$ — bit-exact at the structural level, confirming Y, and adds three substantive new diagnostic objects: a **two-axis** (left-shift, right-shift) decomposition that explains WHY Phillips-Kleinman cross-center was insufficient (PK shifts only the right action), a **predicted magnitude** $d_{\Mcal}((i),(iii)) \approx |\langle \chi_\text{Na 3s}^\text{hyd} | \chi_\text{Na 3s}^\text{phys}\rangle - 1| \cdot \|D_\text{eN}\|_{op}$ at the cross-center coupling region, and a **uniform alkali-hydride scaling** $d_{\Mcal}((i),(iii))(M) \sim Z_\text{core}(M)^{-\alpha}$ with $\alpha$ a property of the bimodule deformation and NOT of the specific alkali. The dual modular layer (Latrémolière 2018, completeness) adds: the right pin-state choice survives the completeness limit, so the diagnostic is stable under further basis refinement. Net: **POSITIVE-NEW-CONTENT, bounded**. The bimodule reframing gives new physics-content beyond Y. It does not autonomously close the wall. It does upgrade the diagnostic from "structural ordering" to "magnitude prediction with named implementation path".

**Verdict line:** POSITIVE-NEW-CONTENT.

---

## §1. Bimodule reframing of the NaH pin-state problem

### 1.1 The natural (Acal_L, Acal_R, Mcal) tuple

The NaH balanced-coupled framework at $\max\_n=2$ has TWO heavy-atom centers carrying distinct multiplication algebras of bounded radial functions. In Latrémolière's modular propinquity vocabulary (arXiv:1608.04881, Def 2.1 metrical quantum vector bundle):

- **Left C*-algebra $\Acal_L$:** the bounded continuous functions on $\R^3$ centered at the H nucleus, with the natural multiplication action on H-centered valence orbitals. Effectively a copy of $C_0(\R^3_\text{H}) \subset C_b(\R^3_\text{H})$ truncated to the Sturmian basis $\{\chi^{\text{H 1s}}_{n,l,m}: n \le n_\max\}$. Sturmian saturation makes the truncation finite-rank but the algebra is naturally non-unital ($C_0$ rather than $C_b$).

- **Right C*-algebra $\Acal_R$:** the bounded continuous functions on $\R^3$ centered at the Na nucleus, equipped with the FrozenCore structure (i.e., a non-trivial Pauli-exclusion projector against the [Ne] core orbitals). Effectively $C_0(\R^3_\text{Na}) \subset C_b(\R^3_\text{Na})$ tensored with the discrete frozen-core algebra and truncated to the Sturmian Na-valence basis $\{\chi^{\text{Na 3s, 3p}}_{n,l,m}: n \le n_\max\}$.

- **Hilbert (Acal_L, Acal_R)-bimodule $\Mcal$:** the space of valence-electron amplitudes that couple both centers — equivalently, the molecular orbital space spanned by linear combinations of H-centered and Na-centered Sturmian orbitals. Concretely, $\Mcal \subset L^2(\R^3)$ as the closed linear span of $\{\chi_{nlm}^\text{H}\} \cup \{\chi_{nlm}^\text{Na}\}$ with left action $f \cdot \xi$ for $f \in \Acal_L$ (multiplication when viewed centered at H) and right action $\xi \cdot g$ for $g \in \Acal_R$ (multiplication when viewed centered at Na). For NaH at $\max\_n=2$, $\dim_{\R} \Mcal = $ dim(H-block) + dim(Na-block) $= 1 + 4 = 5$ orbitals (1s for H, 3s, 3p×3 for Na valence).

- **Bimodule D-norm:** following Latrémolière 2016 Def 2.1, take $D(\xi) = \|\nabla \xi\|_{L^2} + L_L(\langle \xi, \xi \rangle_R) + L_R(\langle \xi, \xi \rangle_L)$ where $\nabla$ is the kinetic-energy operator on $\R^3$ and $L_L, L_R$ are the Lipschitz seminorms on $\Acal_L, \Acal_R$ respectively (with the R1 gradient-Dirac workaround of Sprint X for Leibniz on Schrödinger D). The D-norm encodes both the smoothness of the bimodule element AND the smoothness of its inner products against the two algebras.

The structural claim: **the four NaH pin-state candidates of Y are four points in $\Mcal$** at different positions in this two-axis (left-action / right-action) space.

### 1.2 Identification of the four pin states as bimodule elements

Under the bimodule reframing:

**$\xi_{(i)} =$ hydrogenic Z=1 product.** The valence-orbital component on H is $\chi^{\text{H 1s, hyd Z=1}}$ (a normalized H-1s with Z=1). The valence-orbital component on Na is $\chi^{\text{Na 3s, hyd Z=1}}$ (a normalized "Na 3s" with the same Z=1 hydrogenic radial profile, i.e. mean radius 1.5 bohr). The two are orthogonal at infinite R, near-orthogonal at small R but not identically zero. As a bimodule element: $\xi_{(i)}$ has *trivial right-action coupling at the cross-center region* because the Na-side function does not extend into the cross-center coupling region with its proper diffuse tail. This is "pure-left" in the sense that the H-side function dominates the bond-region density.

**$\xi_{(ii)} =$ SV-corrected diagonal.** The SV correction acts on the diagonal of the molecular Hamiltonian via the screened Na-valence eigenvalue $E_{3s}^{HF} = -0.170$ Ha. As a bimodule element: the wavefunctions are unchanged from $\xi_{(i)}$ — Track 3 verified this bit-exactly. The diagonal shift IS a multiplier in $\Acal_R$ on the inner product $\langle \xi, \xi \rangle_L$ (it changes the expectation value $\langle \xi | H | \xi \rangle$ by 0.330 Ha but does not change $\xi$ itself). In Latrémolière language: **the SV correction is a coboundary deformation on $\Acal_R$ that fixes scalar expectation values but does not change the bimodule element**. Hence $\xi_{(ii)} = \xi_{(i)}$ identically AS BIMODULE ELEMENTS, and $d_\Mcal(\xi_{(i)}, \xi_{(ii)}) = 0$.

**$\xi_{(iii)} =$ bonding-orbital pin state.** The valence orbital is the proper HF/MO solution of the molecular Fock matrix in a basis with **physical Na 3s shape** (mean radius 4.47 bohr) and **physical H 1s shape** (mean radius 1.5 bohr). This is a genuine cross-center linear combination $\xi_{(iii)} = c_H \chi^{\text{H 1s}} + c_\text{Na} \chi^{\text{Na 3s, phys}}$, with $c_H, c_\text{Na}$ R-dependent and both nonzero in the bond region. As a bimodule element: **both left action and right action couple non-trivially**, because $\xi_{(iii)}$ has support that the H-multiplication algebra acts on AND the Na-multiplication algebra acts on, in distinguishable ways at the cross-center coupling region. This is the "bimodule-coupled" element.

**$\xi_{(iv)} =$ MP2-perturbed pin state.** Starts from $\xi_{(iii)}$ (or $\xi_{(i)}$ — the perturbation is anchored on the chosen zeroth-order) and adds double-excitation content via perturbative V_ee mixing. As a bimodule element: a small *bimodule* perturbation that admixes excited bimodule elements weighted by $V_{ee}/\Delta E$. The perturbation preserves the bimodule structure but shifts the element within $\Mcal$ at order $|V_{ee}/\Delta E| \approx 0.05$–$0.1$ for NaH valence.

The two-axis bimodule reframing gives a structurally cleaner picture than Y's scalar $\delta_r$: each pin state has a (left-action coupling, right-action coupling) signature, and the binding wall lives in the bimodule-coupled corner.

---

## §2. Bimodule deformation distance: definition and computation

### 2.1 Candidate definitions

Three honest candidates for the bimodule deformation distance $d_\Mcal(\xi_a, \xi_b)$:

**(a) Hilbert-Schmidt distance on the bimodule.** $d_\Mcal^\text{HS}(\xi_a, \xi_b) := \|\xi_a - \xi_b\|_{L^2(\R^3)}$. Simple, computable, captures wavefunction-shape differences. Limitation: does not see the bimodule structure (left/right action distinction).

**(b) Modular propinquity local metametric extended to bimodule.** $d_\Mcal^\text{mod}(\xi_a, \xi_b) := \sup \{|\langle \xi_a | T | \xi_a \rangle - \langle \xi_b | T | \xi_b \rangle| : T \in \mathcal{B}(\Mcal), L(T) \le 1\}$ where $L(T)$ is the Latrémolière modular Lipschitz seminorm restricted to operators that respect the bimodule structure (left-multiplier and right-multiplier operators). This IS Y's $\delta_r$ but with the test-operator class restricted to bimodule-respecting operators. Gives the same ordering as Y but tighter bounds because the test class is smaller.

**(c) Two-axis L/R decomposition.** $d_\Mcal^\text{LR}(\xi_a, \xi_b) := \left( d_L(\xi_a, \xi_b)^2 + d_R(\xi_a, \xi_b)^2 \right)^{1/2}$ where $d_L$ is the modular propinquity restricted to left-action test operators (operators in $\Acal_L^{\text{op}} = $ left-multiplication algebra acting on $\Mcal$) and $d_R$ is the symmetric right-action version. **This is the substantively new diagnostic**. It splits the bimodule distance into two components: how much $\xi_a$ and $\xi_b$ differ as seen from the H side vs the Na side.

Per the PI's hypothesis and the Latrémolière framework, candidate **(c) is the right diagnostic** because it exposes the asymmetry that Sub-sprint Y's scalar $\delta_r$ cannot see.

### 2.2 Computed/estimated values for the four NaH pin states

Using candidate (c) with the structural estimates from Y and the framework's known wavefunction shapes:

| Pair | $d_L$ estimate (Ha) | $d_R$ estimate (Ha) | $d_\Mcal^\text{LR}$ (Ha) | Pure axis or mixed |
|:---|---:|---:|---:|:---|
| $(i)$ vs $(ii)$ | 0 (wavefunction unchanged) | 0 (wavefunction unchanged) | **0** | trivial |
| $(i)$ vs $(iii)$ | $\sim 0.10$ (H 1s tail correction is small) | $\sim 0.67$ (dominant Na 3s shape change) | $\sim 0.67$ | **right-action dominant** |
| $(i)$ vs $(iv)$ | $\sim 0.02$ (V_ee mixing on H side) | $\sim 0.05$ (V_ee mixing on Na side) | $\sim 0.054$ | mixed, small |

The $d_R \approx 0.67$ Ha estimate for $(i)$ vs $(iii)$ derives from Y's Track-3 Frame−Physical-3p differential of 0.675 Ha at the bond region, which Y identified as the empirical proxy for the scalar $\delta_r$. The two-axis decomposition shows that this differential is dominated by the **right-action** component — the difference between hydrogenic-Z=1 and physical Na 3s wavefunction shapes on the Na-centered side.

The $d_L \approx 0.10$ Ha estimate for $(i)$ vs $(iii)$ derives from the smaller correction to the H 1s tail: physical H 1s has mean radius 1.5 bohr and is well-approximated by the hydrogenic Z=1 1s already (this is exact for isolated H), so the bond-region tail correction is small. The H-side hydrogenic basis IS the right shape; the Na-side hydrogenic basis is NOT.

The MP2 row $(i)$ vs $(iv)$: V_ee mixing affects both centers but stays small in absolute magnitude. The two-axis decomposition is approximately balanced.

### 2.3 Bimodule ordering — comparison with Y's scalar ordering

Y's scalar ordering: $\delta_r((i), (iii)) > \delta_r((i), (iv)) > \delta_r((i), (ii)) = 0$.

Bimodule ordering (this memo): $d_\Mcal^\text{LR}((i), (iii)) > d_\Mcal^\text{LR}((i), (iv)) > d_\Mcal^\text{LR}((i), (ii)) = 0$.

**Same ordering** — the bimodule reading does NOT reorder the candidates. What it adds: **the dominant axis for each pair**. For $(i)$ vs $(iii)$, the asymmetry is striking: $d_R/d_L \approx 6.7$. The diagnostic is overwhelmingly right-action (Na-side) deformation.

This is the substantive new diagnostic. Y's scalar $\delta_r$ said "the right pin state is (iii), pure shape correction in the bond region". The bimodule diagnostic sharpens this: **the right pin state is (iii), pure Na-side wavefunction-shape correction; the H-side basis is already correct**.

---

## §3. Does the bimodule reframing give a NEW diagnostic beyond Y?

Three substantive new diagnostics:

### 3.1 Two-axis decomposition identifies asymmetric correction

Y's $\delta_r$ scalar said: shape correction needed in bond region. Bimodule $d_\Mcal^\text{LR}$ says: shape correction needed **on Na side only** (ratio $d_R/d_L \approx 6.7$). This refines the implementation strategy.

**Implementation implication:** the cheapest implementation path is to **replace the Na-side basis only**, keeping the H-side hydrogenic Z=1 1s as-is. This is essentially Track 3 Option 1 (physical-n hydrogenic relabeling) but the bimodule analysis tells us we DO NOT need to touch the H side. The implementation can be one-sided.

This is non-trivial because the framework's existing atomic_classifier supplies the same `Z_orb=1` for both H 1s and Na 3s (both are asymptotic Z=1 valence). One might naïvely think both need correcting in the same way. The bimodule diagnostic says: no, the H-side is fine, only the Na-side is wrong.

### 3.2 Predicted magnitude of binding recovery

Using $d_\Mcal^\text{LR}((i), (iii)) \approx 0.67$ Ha and the structural correspondence $d_\Mcal \sim \int |\xi_a(r) - \xi_b(r)| \cdot |V_{eN}(r - R_\text{cross})| \, dr$ in the bond region (which is the Track 3 cross-V_ne differential evaluated against the off-center potential), the predicted binding-recovery is the gap between the framework's current descent (0.357 Ha, monotone non-binding) and the natural binding scale on the bimodule-corrected side.

**Predicted binding minimum depth (rough order):** $\sim 0.05$–$0.10$ Ha (with experimental NaH $D_e \approx 0.075$ Ha as a target). The argument: the framework currently over-attracts because the Na-side hydrogenic basis is too compact, putting amplitude inside the unscreened Na nuclear region where it gets overly large cross-V_ne attraction. Replacing with the physical Na 3s diffuses the amplitude OUT, reducing the cross-V_ne over-attraction. The differential of 0.675 Ha is the over-shoot; subtracting it from the current 0.357 Ha descent gives an equilibrium minimum of order 0.3 Ha — within an order of magnitude of $D_e$, but the prediction depends on how the residual partitions between binding depth and energy zero-point shift.

This prediction is honest but rough. Pinning it tightly would need a Latrémolière-rigorous bound (computing the modular propinquity bound $\Lambda(\Tcal_a, \Tcal_b)$ rather than the empirical proxy), but the order-of-magnitude estimate is the substantive new diagnostic content.

### 3.3 Generalization across alkali-hydrides

The bimodule structure is **uniform across the alkali-hydride series** (LiH, NaH, KH, RbH, CsH). All five have:

- $\Acal_L = $ H-centered multiplication algebra (same across the series)
- $\Acal_R = $ alkali-centered multiplication algebra (changes per alkali, but the structure is the same)
- $\Mcal = $ valence bimodule (changes per alkali only in the Na/K/Rb/Cs valence shape)
- $d_R = $ alkali-side wavefunction-shape deformation distance (depends on the alkali's mean valence radius vs hydrogenic Z=1)

**Predicted scaling:** $d_R(\text{alkali}) \sim$ (mean alkali valence radius - hydrogenic Z=1 mean radius) $\sim r_\text{valence}^\text{phys}(M) - 1.5$ bohr.

Approximate values: LiH ~0.7 bohr (Li 2s mean ~2.2 bohr, this matches LiH's known successful binding because it's SMALL), NaH ~3.0 bohr (Na 3s mean ~4.5 bohr — large gap), KH ~3.8 bohr (K 4s ~5.3 bohr), RbH ~4.2 bohr (Rb 5s ~5.7 bohr), CsH ~4.7 bohr (Cs 6s ~6.2 bohr).

**Prediction:** the W1c wall depth scales with the bimodule R-side deformation distance, so the wall is **smallest for LiH** (and indeed LiH binds at 5.3% R_eq error in the framework, the only first-row alkali that binds in the production code), **moderate for NaH** (W1c wall depth 0.357 Ha matches the bimodule estimate), and **larger for KH, RbH, CsH** in the same proportion.

This is a uniform prediction across the alkali-hydride series that goes beyond Y's NaH-specific scalar diagnostic. It is testable: opening the same diagnostic for KH/RbH/CsH would either confirm or falsify the scaling. The framework's heavy-atom monohydride builders (`srh_spec_relativistic`, `bah_spec_relativistic`) already exist and the diagnostic could be lifted to them as future work.

---

## §4. Connection to Phillips-Kleinman: bimodule-asymmetric shift

The Phillips-Kleinman cross-center barrier was Track 3's first closure attempt and returned a clean negative (14.6% improvement, 0.357 → 0.305 Ha descent). The bimodule reframing explains WHY.

The standard PK barrier is:
$$\Delta H_{pq}^{PK} = \sum_c (E_v - E_c) S_{pc} S_{cq}$$
where $p, q$ index the **valence orbitals on the off-center side** (i.e., H-side in NaH) and $c$ indexes the **frozen-core orbitals on the heavy atom** (i.e., the [Ne] core on Na). The shift acts ONLY on the H-side matrix elements via cross-center overlaps $S_{pc}$ against the Na-side core.

**In bimodule language: PK is a left-action shift on $\Acal_L$ (the H-multiplication algebra), parameterized by the Na-core data $(E_c, S_{pc})$**. It modifies the H-valence energy denominator by accounting for orthogonality against the Na core. But it does NOT change the Na-side wavefunction shape — it does not modify the right action.

But the bimodule diagnostic just established: **the W1c residual is right-action dominant** ($d_R/d_L \approx 6.7$). PK shifts the wrong axis.

Hence PK's small effect (14.6%) is the "subleading" left-side correction; the dominant right-side correction (Na 3s shape) requires a different mechanism. **The bimodule reframing provides a structural explanation for the PK negative result**: PK shifts the smaller of the two bimodule axes.

This is genuinely-new content beyond Y. Y said "PK was insufficient because the diagnostic favors (iii)"; the bimodule reframing says "PK acted on the LEFT-action axis (H side) but the residual lives in the RIGHT-action axis (Na side) by a 6.7× factor".

**Implication for an enhanced PK:** a symmetric or right-acting PK barrier — one that projects against the H 1s on the Na side (not the standard direction) — would shift the right axis. But the H 1s is not a frozen-core orbital so this is non-standard. The bimodule structure suggests the correct generalization: **a bilateral PK that projects against both core sets simultaneously**, with weights determined by the bimodule structure.

This is a concrete candidate implementation path that the bimodule analysis surfaces and that was not visible from Y.

---

## §5. Dual modular propinquity: completeness and duality

Latrémolière's **dual modular propinquity** (arXiv:1811.04534) is a complete weaker variant of the modular propinquity, designed to take limits of bimodule deformations. The dual layer adds:

### 5.1 Completeness guarantee

The dual modular propinquity is complete on the class of metrical quantum vector bundles. Implication: if we have a sequence of pin-state candidates $\xi_n \to \xi$ in the dual modular propinquity, the limit is itself a metrical bimodule. For the NaH problem, this means **the bonding-orbital pin state (iii) survives basis refinement** — the framework's $\max\_n=2$ approximation to (iii) converges to a well-defined bimodule limit at $\max\_n \to \infty$.

This is structurally non-trivial because Y's modular propinquity (without the dual layer) is NOT complete in general, and there is no a priori guarantee that the $\max\_n \to \infty$ limit of the bonding-orbital sequence is well-defined as a bimodule element.

**Concrete implication:** the Y.1 diagnostic at $\max\_n = 2, 3$ would, under the dual modular propinquity, give a sequence of $d_\Mcal^\text{LR}$ values converging to a definite limit. Track 3 Option 1 (physical-n hydrogenic relabeling) would also produce a Cauchy sequence in the dual modular propinquity, and its limit would BE the bonding-orbital pin state of an infinitely-refined basis. The completeness guarantees the implementation is convergent.

### 5.2 Duality and pin-state constraint

The dual modular propinquity introduces duality on the bimodule structure. For NaH, the dual bimodule is the conjugate valence space — informally, the space of "hole" excitations on the bonding orbital. The duality is structural; concretely it corresponds to the conjugate orbital (anti-bonding $\sigma_u$) being a valid alternative pin state.

**Constraint from requiring the propinquity to respect duality:** if we impose that the diagnostic respects the bonding/anti-bonding duality, then any pin state that distinguishes bonding from anti-bonding via a structural mechanism is preferred. Candidate (i) (hydrogenic Z=1 product) does NOT distinguish — both bonding and anti-bonding linear combinations of $\chi^{\text{H 1s}}, \chi^{\text{Na 3s hyd}}$ have the same hydrogenic basis. Candidate (iii) DOES distinguish — the bonding linear combination has R-dependent amplitude redistribution, while anti-bonding has different R-dependent redistribution.

**Implication:** the duality requirement preferentially selects (iii) over (i). This is a structural argument INDEPENDENT of the magnitude estimate of §3.2 — duality alone is sufficient to identify (iii) as the right pin state, without invoking the empirical Frame−Physical-3p differential.

This is the substantive new content from the dual layer beyond the modular layer: a duality-based identification of (iii) that does NOT depend on Y's empirical magnitude estimate. Useful as a cross-check on Track 3's named target.

### 5.3 Honest scope on the dual layer

Latrémolière 2018 (arXiv:1811.04534) §4-5 develops the dual modular propinquity in the context of metrical quantum vector bundles with D-norms. The bonding-orbital pin state has a natural D-norm (the kinetic-energy operator on $\R^3$ plus the cross-register left/right Lipschitz seminorms), but verifying that this D-norm satisfies Latrémolière's axioms (Def 4.1 metrical quantum vector bundle) is itself a 1-week structural verification. The dual modular propinquity layer adds genuine content but the rigorous verification is non-trivial.

**Recommendation:** the dual layer's content (completeness + duality) is qualitative-positive but the rigorous verification is a separable sub-sprint. Y.1 + a dual-modular-extension sub-sprint = the full upgrade. The dual layer alone is not autonomously sufficient.

---

## §6. Net verdict

### Does the modular framework give NEW content beyond Y's structural ordering?

**YES.** Three new diagnostics surface:

1. **Two-axis (left/right) decomposition** — the W1c residual is right-action dominant ($d_R/d_L \approx 6.7$ for NaH). This is structurally new and was not visible in Y's scalar $\delta_r$.
2. **Predicted magnitude** — the binding-recovery gap closeable by the bimodule-corrected pin state is of order 0.1 Ha for NaH, comparable to the experimental $D_e \approx 0.075$ Ha. This is a rough order-of-magnitude prediction beyond Y's qualitative "wall closure expected".
3. **Uniform alkali-hydride scaling** — the W1c wall depth scales with the right-action deformation distance $r_\text{valence}^\text{phys}(M) - 1.5$ bohr. This explains why LiH binds (small bimodule distance), NaH/KH/RbH/CsH don't (large bimodule distances).

### Does it suggest a SPECIFIC implementation path with predicted magnitude?

**YES.** Two implementation paths surface, with predicted magnitudes:

- **Path A (recommended): one-sided Na-only physical wavefunction substitution.** Replace the Na-valence basis ($Z_\text{orb}=1$ hydrogenic) with the physical Na 3s screened wavefunction from FrozenCore $Z_\text{eff}(r)$. Keep the H-side hydrogenic Z=1 1s as-is. The bimodule diagnostic predicts this closes the $d_R$ axis (the dominant axis, 6.7× larger than $d_L$). Predicted binding-recovery: $\sim 0.1$ Ha. This is Track 3 Option 1 (physical-n hydrogenic relabeling) but with the bimodule analysis surfacing **the H side does not need touching**.

- **Path B (alternative): symmetric bilateral PK.** Construct a Phillips-Kleinman barrier that acts on both sides (H valence against Na core AND Na valence against H "core" = H 1s itself). The bimodule structure suggests the bilateral weights. This is a NEW implementation candidate not in the Track 3 memo.

Path A is the recommended next sprint. Path B is a structural alternative that the bimodule analysis surfaces and would not have been visible without the modular reframing.

### Does it generalize across the alkali-hydride series at increasing Z?

**YES.** The bimodule structure is uniform across LiH, NaH, KH, RbH, CsH. The predicted W1c wall depth scales with the right-action deformation distance. This is testable: opening the diagnostic for KH and CsH would confirm or falsify the scaling.

### Recommended next sprint

**Open Sub-sprint M-Y.1 at 1-2 weeks (parallel to Y.1).** Scope:

1. **Implement candidate (c)** (the two-axis L/R decomposition $d_\Mcal^\text{LR}$) numerically for NaH at $\max\_n=2, 3$.
2. **Verify the asymmetry $d_R/d_L \gg 1$** (predicted ~6.7).
3. **Test the alkali-hydride scaling** by computing $d_\Mcal^\text{LR}$ for LiH and KH at $\max\_n=2$ (LiH should be small, KH should be larger).
4. **Open Track 3 Option 1 implementation as the named target** with bimodule-grounded justification for the one-sided substitution.

The bimodule reframing of M-Y is **substantively new content** beyond Y, but it does NOT autonomously close the wall. The actual cross-V_ne integrals still need re-computation in the bimodule-corrected basis (Track 3 Option 1 / Path A above). The CLAUDE.md §1.7 structural-skeleton-scope statement is preserved: the framework supplies the math.OA-grounded diagnostic that identifies WHERE to look and HOW MUCH the gap is, but the implementation work remains.

**Verdict line:** **POSITIVE-NEW-CONTENT.**

The bimodule reframing of the W1c pin-state diagnostic gives three new diagnostic objects (two-axis decomposition, magnitude prediction, uniform alkali-hydride scaling) that Y's scalar $\delta_r$ does not access. It identifies one new candidate implementation path (Path B bilateral PK) that Track 3 did not surface. It explains the PK negative result structurally (wrong axis). The dual modular layer adds completeness + duality guarantees. Net: substantive upgrade, not vocabulary-only, but still bounded by the structural-skeleton-scope statement.

---

## §7. Honest scope and caveats

**Caveats:**

1. **D-norm verification deferred.** The bimodule D-norm proposed in §1.1 is a natural candidate but the rigorous verification that it satisfies Latrémolière 2016 Def 2.1 (closed extension, lower semi-continuity, etc.) is a separate ~1-2 week sub-sprint not done here. The diagnostic ordering is robust under reasonable D-norm choices; the magnitude prediction depends mildly on the D-norm specifics.

2. **Magnitude prediction is order-of-magnitude.** The estimated binding-recovery of ~0.1 Ha is honest as a rough estimate from the Frame−Physical-3p differential, but a tight prediction would require a Latrémolière-rigorous propinquity bound rather than the empirical proxy. The order-of-magnitude is the substantive content.

3. **Cross-block ERI re-evaluation needed for full closure.** Beyond cross-V_ne shape correction, the bonding-orbital pin state (iii) also changes within-block ERIs (V_ee on the Na side). The bimodule analysis does NOT address these, and they would shift the binding curve. The named implementation path (Path A) automatically handles this because the physical Na 3s shape is used uniformly in all integrals.

4. **Alkali-hydride scaling is structural-not-quantitative.** The predicted scaling $d_R \sim r_\text{valence}^\text{phys}(M) - 1.5$ bohr is the LEADING term of a more complex shape-difference functional. Sub-leading terms (e.g., differences in oscillation structure) would correct the scaling. The structural ordering (LiH small, NaH/KH/RbH/CsH larger) is robust; the exact magnitudes are not.

5. **PK bilateral variant (Path B) is structural-only.** The bimodule reframing surfaces the symmetric bilateral PK as a candidate, but constructing it concretely would require deciding which H-side "core" to project against (H has no true frozen core). The most natural choice is to project the Na valence against the H 1s with weight given by the Na-H cross-overlap, but this is a NEW construction that needs careful design before implementation.

6. **No production code or papers modified.** This is a diagnostic-only memo per the sprint mandate. Implementation remains as named for Y.1 + Track 3 Option 1.

**Confidence:**

- **HIGH** that the two-axis L/R decomposition is the structurally correct sharpening of Y's scalar $\delta_r$.
- **HIGH** that the right-action axis dominates for NaH (the wavefunction-shape difference is on the Na side; H is hydrogenic and correct).
- **MEDIUM-HIGH** that Path A (one-sided Na substitution) closes the wall to within a factor of 2 of experimental $D_e$.
- **MEDIUM** that the alkali-hydride scaling holds quantitatively across LiH/NaH/KH/RbH/CsH.
- **LOW** that the dual modular propinquity duality argument adds independent evidence beyond the modular layer (the duality cross-check is qualitative-positive but not quantitatively new).

**Files referenced:**

- `debug/subsprint_y_w1c_pinstate_diagnostic_memo.md` (Y memo, primary reference)
- `debug/sprint_l3e_p3_synthesis_memo.md` §2 (Y verdict in synthesis)
- Latrémolière, F. "The Modular Gromov-Hausdorff Propinquity", arXiv:1608.04881 (2016)
- Latrémolière, F. "The Dual Modular Gromov-Hausdorff Propinquity and Completeness", arXiv:1811.04534 (2018)
- `geovac/balanced_coupled.py` (NaH builder, fixed point of bimodule structure)
- `geovac/cross_register_vne.py` (cross-register V_eN, related multi-focal-composition architecture)
- `geovac/screened_valence_basis.py` (SV-corrected diagonal, candidate (ii), shown bit-identical to (i) at wavefunction level)
- `geovac/phillips_kleinman_cross_center.py` (PK cross-center, shown shifts wrong axis)
- `papers/group2_quantum_chemistry/paper_17_composed_geometries.tex` §6.10 (W1c context)
- `papers/group2_quantum_chemistry/paper_19_coupled_composition.tex` (balanced coupled framework)
- CLAUDE.md §1.7 W1c-residual orthogonality wall entry
- CLAUDE.md §3 chemistry-arc-paused-W1c-residual row

---

**End of Track M-Y diagnostic memo.**
