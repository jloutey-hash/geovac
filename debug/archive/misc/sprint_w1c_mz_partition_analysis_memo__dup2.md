# Sprint W1c × M-Z partition bridge — structural analysis

**Date:** 2026-05-23 (same-day bridge sprint, post-F1 P1+P2 and post-modular synthesis).
**Sprint position:** One-shot bridge sprint applying the modular propinquity M-Z partition (cross-shift vs endomorphism, Track M-Z 2026-05-23) to the W1c-residual three-layer hierarchy (Sprint F1 P1+P2 2026-05-23). Both findings landed today; both are structural; the bridge is the right scope for a structural-only analysis sprint.
**Cross-references:** `debug/sprint_modular_propinquity_mZ_bethe_log_memo.md` (M-Z primary), `debug/sprint_modular_propinquity_synthesis_memo.md` (modular synthesis §3 cross-shift vs endomorphism partition), `debug/sprint_f1_p1p2_combined_test_memo.md` (F1 primary, three-layer W1c hierarchy), `debug/sprint_modular_alpha_arc_synthesis_memo.md` (α arc synthesis with the original Layer 3 sketch), `debug/sprint_alpha_1_diagnostic_memo.md` (bimodule distance d_R/d_L = 3.23), CLAUDE.md §1.7 multi-focal-composition wall taxonomy.
**Mandate:** Diagnostic/structural only. No production code modifications, no paper edits, no CLAUDE.md edits. Classify each W1c sub-layer under M-Z; produce falsifiable predictions for F1 max_n=3.

---

## §0. Executive summary + verdict

### Verdict line: **PARTITION-EXTENDS-CLEANLY** (with one sub-classification refinement)

All three F1 W1c sub-layers classify cleanly under the M-Z partition:

| Sub-layer | Mechanism | M-Z classification |
|:---|:---|:---|
| 1 — H ↔ Na core orthogonality | Phillips-Kleinman cross-center barrier | **Bimodule cross-shift** (handled by framework, empirically 14.6% W1c-descent reduction) |
| 2 — Na-side wavefunction shape | Multi-zeta substitution of Na 3s/3p | **Module endomorphism** (external input from FrozenCore + Clementi-Roetti fits; bit-zero without cross-shift activator, load-bearing once cross-shift activates Na occupancy — 23% additional W1c-descent reduction) |
| 3 — Orbital-pair flexibility (no bonding/antibonding mixing at max_n=2) | FCI cannot construct [H 1s ± Na 3s]_± combinations at limited basis | **Basis-closable bimodule cross-shift** — class is cross-shift (the linear combination is structurally bimodule), obstruction is basis dimensionality, predicted to close at max_n=3 |

The substantive new content of this sprint is the **third bucket**: the M-Z partition as originally stated has two buckets (cross-shift / endomorphism); sub-layer 3 introduces a sub-classification *within* the cross-shift bucket. The basis-closable refinement was not visible from the bound-state QED side (M-Z's original context) because Sturmian basis truncations on $L^2(\R^3)$ saturate the bound-state space at low $N_\max$ for the specific pin state being studied. On the chemistry side with multi-determinant FCI on a multi-center bimodule, the basis dimensionality directly limits which cross-shift combinations the variational state can realize.

### Three falsifiable predictions for F1 max_n=3

**P1 (Internal minimum emerges):** NaH F1 max_n=3 with combined W1c + multi-zeta produces an internal PES minimum at $R_\text{eq} \in [3.0, 4.5]$ bohr. **Falsifier:** R_min at the smallest tested R value (e.g., 2.0 bohr).

**P2 (Binding within 2×):** If P1 holds, the binding energy $D_e \in [0.0375, 0.150]$ Ha (factor-2 around experimental NaH $D_e \approx 0.075$ Ha). **Falsifier:** $D_e$ outside this range.

**P3 (Multi-zeta differential persists):** $|E(\text{W1c+mz}) - E(\text{W1c})|$ at $R_\text{eq}$ is $\in [0.02, 0.30]$ Ha. **Falsifier:** differential <0.005 Ha (mz absorbed by basis expansion) or >0.30 Ha (anomalous basis-set scaling).

Confidence: P1 MEDIUM-HIGH, P2 MEDIUM, P3 HIGH.

### Recommendation

Open F1 max_n=3 sprint (3-5 days) testing P1/P2/P3 simultaneously. Promote M-Z partition to CLAUDE.md §1.7 and Paper 34 §V.D.7 only if P1 confirmed; otherwise keep as local-vocabulary refinement.

---

## §1. Bimodule structure of the NaH balanced architecture

### 1.1 The (A_L, A_R, M) tuple

The NaH balanced-coupled Hamiltonian (`geovac/balanced_coupled.py`) decomposes the molecular valence space into two sub-blocks (per the `_get_block_geometry` partition):

- **`NaH_bond_partner`** — H-side sub-block, dimension 5 at max_n=2, basis is hydrogenic Z=1 at H center.
- **`NaH_bond_center`** — Na-side sub-block, dimension 5 at max_n=2, basis is hydrogenic Z=1 at Na center (placeholder) OR multi-zeta-fitted physical Na 3s/3p (when `multi_zeta_basis=True`).

The natural bimodule tuple:

- **$\mathcal{A}_L = C_0(\R^3_H)$** — multiplication algebra of coordinate functions localized near H. Acts on H sub-block states via diagonal multiplication; cross-V_ne mixing terms involving H-side bra/ket states are A_L-active.
- **$\mathcal{A}_R = C_0(\R^3_\text{Na}) \otimes \text{FrozenCore}(Z=11)$** — multiplication algebra of coordinate functions localized near Na, compressed by the [Ne] frozen-core orthogonalization. Acts on Na sub-block states; cross-V_ne mixing terms involving Na-side bra/ket states are A_R-active.
- **$\mathcal{M}$ = valence bimodule** — 10-dimensional Hilbert space spanned by 5 H-side + 5 Na-side hydrogenic orbitals at max_n=2; 2-electron FCI dimension = $C(10,1)^2 = 100$.

### 1.2 Are the actions truly two-sided?

**YES, in the M-Y / M-Z sense.** The decomposition is two-sided because:

- $\mathcal{A}_L$ acts non-trivially on H sub-block matrix elements via on-site potential (H nucleus attraction) and cross-V_ne (Na nucleus attraction reaches H states via multipole expansion centered at H).
- $\mathcal{A}_R$ acts non-trivially on Na sub-block matrix elements via on-site potential (Na nucleus attraction with FrozenCore screening) and cross-V_ne (H nucleus attraction reaches Na states via multipole expansion centered at Na).
- $\mathcal{A}_L$ does NOT act on Na sub-block matrix elements (and vice versa) — coordinate functions centered at H do not have direct matrix elements between Na-centered orbitals.
- Cross-block ERIs (two-body Coulomb integrals between H-side and Na-side orbitals) constitute the bimodule cross-shift kernel — they are the canonical signature of a true bimodule (each integrand is bilinear in left-multiplication × right-multiplication).

This matches the M-Y description verbatim (`debug/sprint_modular_propinquity_synthesis_memo.md` §3): "$\mathcal{A}_L$ = H-centered multiplication algebra... $\mathcal{A}_R$ = Na-centered with FrozenCore structure... $\mathcal{M}$ = valence bimodule (dim 5 at max_n = 2)" — dim 5 per sub-block, dim 10 total for the bimodule basis.

### 1.3 Tensor product vs bimodule

One could alternatively read the structure as a tensor product $\mathcal{M} = \mathcal{M}_H \otimes \mathcal{M}_\text{Na}$ where each factor is a single-center module over a single multiplication algebra. **The bimodule reading is structurally more accurate** because:

- The cross-V_ne kernel is NOT a tensor product of single-center operators — it is a bilinear coupling $\int \rho_H(r_H) \cdot (-Z_\text{Na}/|\mathbf{r}_H - \mathbf{R}_\text{Na}|) \cdot \rho'_H(r_H) \, dr_H$ which involves the H-sub-block density evaluated against the Na-centered Coulomb kernel.
- In bimodule language, this is precisely the $\mathcal{A}_L$-action of "H-localized density" composed with the $\mathcal{A}_R$-input of "Na-localized point-charge potential evaluated in H coordinates."
- The Wigner-D rotation machinery (`geovac/shibuya_wulfman.py`) that handles non-collinear cross-V_ne is the algebraic concretization of how $\mathcal{A}_L$ and $\mathcal{A}_R$ communicate through the cross-shift kernel.

The tensor-product reading would be appropriate only if the cross-V_ne were absent (which would correspond to NO bimodule structure, just two non-interacting modules). The presence of cross-V_ne with non-trivial matrix elements between H and Na sub-blocks is what makes this a genuine bimodule.

### 1.4 Connection to the bound-state QED bimodule (M-Z)

M-Z's bimodule reading of hydrogen Lamb shift had:
- $\mathcal{A}_e = C_0(\R^3)$ (electron coordinate algebra)
- $\mathcal{M}_e = L^2(\R^3) \otimes \C^4$ (Dirac spinor bundle)
- Cross-register $V_{eN}$ as bimodule cross-shift to nuclear module $\mathcal{M}_n$

The NaH bimodule is the *molecular* analog: two coordinate algebras ($C_0(\R^3_H)$ and the Na-centered FrozenCore-compressed analog), one molecular valence bimodule, cross-V_ne as the cross-shift kernel. The structural parallel is faithful — both bimodules have the same shape (two algebras × one bimodule × cross-shift kernel), differing only in the physical scale (single-particle bound-state QED vs many-body molecular FCI).

---

## §2. Sub-layer 1 (PK orthogonality) classification

### 2.1 The mechanism

Phillips-Kleinman cross-center barrier (`geovac/phillips_kleinman_cross_center.py`, May 2026):

$$\Delta H^\text{PK}_{pq} = \sum_c (E_v - E_c) S_{pc} S_{cq}$$

where $p, q$ are H-valence orbital indices (in $\mathcal{M}_H$), $c$ ranges over Na [Ne]-core orbitals, $S_{pc} = \langle p_H | c_\text{Na} \rangle$ is the cross-center overlap, and $(E_v - E_c)$ is the energy gap between H valence energy and Na core eigenvalue (taken from Clementi-Roetti HF tabulation).

### 2.2 Classification: BIMODULE CROSS-SHIFT

The PK kernel is bilinear in:
- $S_{pc}$: cross-overlap between H-valence orbital and Na-core orbital. This integrand has an H-localized component (from $p$) and a Na-localized component (from $c$), so the overlap is structurally bimodule.
- $(E_v - E_c)$: a scalar (eigenvalue difference), which weights the contribution but does not change the structural type.

The full $\Delta H^\text{PK}_{pq}$ has matrix elements between H-valence indices $p, q$ (i.e., it lives in $\mathcal{M}_H^* \otimes \mathcal{M}_H$, an endomorphism of the H module). **But** its construction is bilinear in the cross-overlap $S_{pc}$, which is the canonical bimodule cross-shift bilinear: it mixes the left-action algebra ($\mathcal{A}_L$, H coordinate functions, acting on $p, q$ via their location) with the right-action algebra ($\mathcal{A}_R$, Na coordinate functions, acting on $c$ via its location).

The framework's bimodule machinery composes this: $\Delta H^\text{PK}_{pq}$ is built from a sum over $c \in $ Na core, each term contributing a bilinear $\mathcal{A}_L \otimes \mathcal{A}_R$ matrix element. This is the structural signature of a bimodule cross-shift in the Latrémolière framework, exactly as M-Z identified for cross-register $V_{eN}$ in bound-state QED.

### 2.3 Consistency check: framework handles it

PK is implemented in production code (`geovac/phillips_kleinman_cross_center.py`, ~280 lines, 21 tests). The F1 P1+P2 result confirms it is effective: 14.6% W1c-descent reduction at NaH max_n=2 (F1-P2 §3.2 ladder: bare 6.244 → +W1c 0.357 → +W1c+PK 0.305 Ha). This empirical match is what M-Z predicts for handled cross-shifts.

### 2.4 Mapping to M-Y's prediction

M-Y predicted that PK acts on the "wrong axis" (left-action only) because the W1c residual is right-action dominant ($d_R / d_L = 3.23$ per α-1 diagnostic). The bimodule classification of PK as a cross-shift is consistent with this — PK *is* a bimodule cross-shift, but its specific construction (cross-overlap with Na core, energy weighting with H valence) sets up matrix elements on the H sub-block only. The right-action axis (Na-side wavefunction shape) needs a separate bimodule cross-shift mechanism (multi-zeta + cross-V_ne integration), which is sub-layer 2's domain.

**Net for sub-layer 1: cross-shift classification is unambiguous. Framework handles it. The 14.6% improvement is the legitimate contribution of one specific bimodule cross-shift mechanism (PK against Na core) to the W1c residual.**

---

## §3. Sub-layer 2 (multi-zeta Na shape) classification

### 3.1 The mechanism

Multi-zeta substitution (`geovac/multi_zeta_orbitals.py` Z=11 entry, α-2 sprint 2026-05-23): replace the hydrogenic-Z=1 Na 3s/3p radial wavefunction with a K=5 mixed-Slater-n expansion that fits the FrozenCore-screened physical Na valence to overlap > 0.999999.

Concretely: the Na sub-block basis functions $\phi_\text{Na}(r) = R^\text{hyd}_{Z=1,n=3,l}(r) \cdot Y_{l,m}(\hat{r})$ are replaced by $\phi^\text{phys}_\text{Na}(r) = R^\text{phys}_\text{Na}(r) \cdot Y_{l,m}(\hat{r})$ where $R^\text{phys}_\text{Na}$ is the mixed-Slater-n fit (K=5 primitives with $n_\text{slater} = [1, 1, 2, 3, 3]$ and coefficients fit to the FrozenCore eigenstate).

### 3.2 Classification: MODULE ENDOMORPHISM

Multi-zeta changes the radial shape of basis functions in the Na sub-block:
- It is internal to $\mathcal{M}_R$ (the Na sub-module) — the H sub-module is untouched.
- It does NOT directly couple H and Na coordinates — the substitution is on a single sub-module's basis function shape.
- The bimodule structure $(\mathcal{A}_L, \mathcal{A}_R, \mathcal{M})$ is preserved; what changes is the choice of basis for $\mathcal{M}_R$ within that structure.

Compare to PK (sub-layer 1): PK adds a *new bilinear matrix element* mixing H and Na coordinates. Multi-zeta does NOT add new matrix elements; it changes the radial integrand within existing matrix elements that were already cross-shift (cross-V_ne) or already endomorphism (on-Na ERIs).

This is the canonical signature of a **module endomorphism** in M-Z's vocabulary: a self-deformation of a single sub-module that does not generate new bimodule cross-coupling.

### 3.3 The "bit-zero without W1c, load-bearing with W1c" pattern

The most distinctive empirical signature of the endomorphism classification is the mutual-exclusivity behavior found in F1:

- **Without W1c (bare cross-V_ne):** Multi-zeta is bit-zero on FCI eigenvalues (α-PES Step 2 finding, reproduced in F1-P2 §3.1: $E_\text{bare} = E_\text{bare+mz}$ to floating-point precision at every tested R).
- **With W1c:** Multi-zeta is load-bearing (W1c+mz vs W1c-alone gives +0.056 Ha differential at R=3.5, scaling to +0.208 Ha at R=2.0, F1-P2 §3.2).

**Structural explanation under the endomorphism reading:** the Na 3s endomorphism content (shape difference between hydrogenic Z=1 and physical Na 3s) only becomes visible at the FCI level when the Na 3s orbital is occupied. Under bare cross-V_ne, the un-screened Z=11 attraction pulls all FCI weight to H-side orbitals (Na 3s sits at h1 eigenvalue index i=5, unoccupied at 2-electron singlet ground state; F1-P2 §4.2). The cross-shift (cross-V_ne) does not activate the Na sub-module's occupation, so the endomorphism content on that sub-module is invisible to the eigenvalue.

W1c screening removes the un-screened Z=11 attraction. The Na 3s eigenvalue shifts to i=1 (second-lowest h1 eigenvector at R=3.5, amplitude² = 0.986; F1-P2 §4.2), and the FCI naturally engages it with ~0.98 occupation. Now the endomorphism content (multi-zeta substitution) propagates to the FCI eigenvalue, producing the +0.056 to +0.208 Ha differential.

**This is precisely what M-Z predicts for endomorphism content:** it is structurally external input to the bimodule machinery, but the bimodule cross-shifts (here W1c + cross-V_ne) compose it into observable matrix elements once the relevant sub-module is in play.

### 3.4 Framework handles partially

The endomorphism classification means: the framework does NOT generate the multi-zeta exponents and coefficients from its own data. They come from:
- FrozenCore radial Schrödinger eigenstates (a separate atomic solver)
- Clementi-Raimondi 1963/1967 or BBB93/Clementi-Roetti 1974 atomic basis fits (literature input)

This is the LS-8a-class wall on the chemistry side: external input is required for the endomorphism content. **HOWEVER**, the chemistry-side endomorphism is *computable* (via FrozenCore + atomic-physics fits), whereas the QED-side endomorphism (LS-8a counterterms $Z_2 - 1$, $\delta m$) requires UV-completion theory. See §6 for the structural-vs-contingent argument.

### 3.5 Consistency with the α arc's bit-zero finding

The α-PES Step 2 BARE-multi-zeta test (`debug/sprint_alpha_3_pes_test_memo.md`) found multi-zeta bit-zero on FCI eigenvalues. The α arc synthesis (`debug/sprint_modular_alpha_arc_synthesis_memo.md` §4) interpreted this as "Layer 3 — FCI variational state H-side localization at finite basis." The endomorphism classification refines this interpretation:

- The α arc's reading: "the bimodule diagnostic predicts wavefunction-shape distance correctly but the FCI doesn't occupy that orbital, so the prediction is invisible to FCI."
- The M-Z partition reading: "multi-zeta is endomorphism content; endomorphism is only visible when a cross-shift activates the module's occupation. Bare cross-V_ne does not activate Na occupation; W1c does."

The two readings agree on the mechanism but differ in framing. The M-Z partition reading is structurally cleaner because it places the bit-zero finding inside a general partition (any endomorphism is bit-zero on FCI without an activator cross-shift, regardless of basis size or system).

**Net for sub-layer 2: endomorphism classification is unambiguous. Framework does not autonomously generate the content but composes it into observable matrix elements via bimodule cross-shifts. Empirical mutual-exclusivity behavior (bit-zero without cross-shift, load-bearing with cross-shift) is the structural signature.**

---

## §4. Sub-layer 3 (orbital-pair flexibility) classification

### 4.1 The mechanism (F1's revised statement)

F1-P2 §4.3 revised the original α arc Layer 3 framing. The empirical finding:

- At NaH max_n=2 with combined W1c + multi-zeta, FCI engages Na 3s at full occupation (~0.98).
- Natural occupations are [1.000, 1.000] for the two near-degenerate orbitals (one Na-localized, one H-localized).
- The framework's max_n=2 basis does NOT construct bonding/antibonding linear combinations $[H 1s \pm Na 3s]_\pm$ with one bonding pair occupied; instead it has two separately-occupied orbitals.
- The PES is monotonically descending (no internal minimum) because the bonding-vs-antibonding partition would require placing the symmetric combination at lower energy than the separated configuration, and at max_n=2 this combination is not lower.

**The wall is at the level of orbital-orbital mixing, not at the level of FCI occupation.**

### 4.2 Classification: BIMODULE CROSS-SHIFT (basis-closable)

The bonding/antibonding combinations are explicitly:

$$\phi_\pm = \frac{1}{\sqrt{2}} (\phi_{H,1s} \pm \phi_{\text{Na},3s})$$

These are *bimodule linear combinations*: they couple a left-action ($\mathcal{A}_L$) basis element ($\phi_{H,1s}$, H-localized) with a right-action ($\mathcal{A}_R$) basis element ($\phi_{\text{Na},3s}$, Na-localized). The $\pm$ sign coefficients are the cross-shift weights.

**The class is cross-shift.** A bonding orbital is structurally a left+right combination in the bimodule basis. The framework's bimodule machinery DOES support such combinations — they are eigenvectors of the h1 matrix when off-block elements (cross-V_ne, in particular) are non-trivial.

**The obstruction is basis dimensionality.** At max_n=2, the basis has 10 single-particle orbitals (5 H + 5 Na). The FCI variational state can in principle construct any linear combination of these via the CI coefficient vector, including bonding/antibonding pairs. But the *lowest-energy* configuration at max_n=2 happens to be two separate orbitals each at ~−0.79 Ha, not a bonding combination at lower energy. The latter would require additional orbital flexibility — H 2s/2p polarization, Na 3p mixing — to construct a bonding combination that is actually lower in energy than the separated configuration.

At max_n=3, the basis adds H 2s, 2p, Na 3p, 3d (and additional radial flexibility), giving 5 → 14 orbitals per center, 28 total. The cross-V_ne + ERI matrix elements between these additional orbitals provide the dimensional richness for the FCI to construct a polarized bonding orbital combination with energy lower than the separated configuration.

### 4.3 Why "basis-closable" is a sub-classification of cross-shift, not a new class

In the M-Z partition's original framing, cross-shifts are "HANDLED by the framework" and endomorphisms are "NOT HANDLED." Sub-layer 3 is structurally cross-shift but is conditionally handled — handled only at sufficient basis size. This is a sub-classification:

- **Cross-shift, basis-handled:** the cross-shift mechanism is realized in the FCI variational state at the basis size in use. PK (sub-layer 1) is in this class.
- **Cross-shift, basis-closable:** the cross-shift mechanism exists in the framework but the basis size is insufficient for the FCI to realize it. Sub-layer 3 is in this class.
- **Endomorphism:** structurally external input required, not generated by the framework. Sub-layer 2 (multi-zeta) is in this class.

The "basis-closable" qualifier is genuinely new content beyond M-Z's original two-bucket partition. It emerges naturally from FCI variational state analysis on a multi-determinant bimodule, which is a feature M-Z's bound-state QED context did not have (Lamb shift uses fixed-state diagrams, not variational FCI).

### 4.4 Predicted behavior at max_n=3

If the basis-closable classification is correct, max_n=3 should provide the dimensional richness for the FCI to construct bonding orbital combinations as eigenvectors of the h1 + W1c + multi-zeta + cross-V_ne effective single-particle Hamiltonian. The FCI ground state at max_n=3 should:

- Have natural occupations closer to [2.0, 0.0] for a true bonding/antibonding partition (the bonding pair is doubly occupied; the antibonding is empty), as opposed to max_n=2's [1.0, 1.0] separate-orbital occupancy.
- Have one occupied orbital that is a polarized bonding combination with substantial weight on both H 1s and Na 3s (mixed Slater-n by multi-zeta + augmented by H 2s/2p, Na 3p polarization functions).
- Produce an internal PES minimum at $R_\text{eq}$ closer to experimental 3.566 bohr.

If the classification is incorrect (sub-layer 3 is actually a basis-NON-closable endomorphism), max_n=3 will NOT produce a binding minimum and the framework's chemistry-side analog of the LS-8a wall will be deeper than just "external atomic structure input."

This is the falsifiable test for the bridge prediction.

---

## §5. Synthesis: testable predictions for F1 max_n=3

### 5.1 Predictions matrix

| Sub-layer | Expected behavior at max_n=3 vs max_n=2 | Mechanism |
|:---|:---|:---|
| 1 (PK cross-shift) | Same content, comparable contribution | Cross-shift mechanism is basis-size-independent at structural level; new max_n=3 orbitals add small additional PK contributions but the dominant H-valence × Na-core overlap is dominated by the n=1, 2 core states which are present at both max_n levels |
| 2 (multi-zeta endomorphism) | Endomorphism content fixed (Na shape is what it is); cross-shift bilinear (cross-V_ne) integrates against more orbital pairs, so total visible content slightly larger | Endomorphism content per matrix element is bounded by the Na shape data; max_n=3 adds more matrix elements where this content shows up, but the per-cell contribution is bounded |
| 3 (basis-closable cross-shift) | PES qualitatively changes — internal minimum should emerge if classification is correct | The basis adds H 2s, 2p, Na 3p, 3d giving sufficient dimensional richness for the FCI to construct polarized bonding orbital combinations |

### 5.2 The three falsifiable predictions

#### P1: Internal PES minimum emerges
- **Statement:** NaH F1 max_n=3 with combined W1c + multi-zeta architecture produces an internal PES minimum at $R_\text{eq} \in [3.0, 4.5]$ bohr.
- **Falsifier:** $R_\min$ at the smallest tested R value (e.g., 2.0 bohr), i.e., monotonically descending PES.
- **Underlying classification:** Sub-layer 3 = basis-closable cross-shift; max_n=3 provides the dimensional richness for bonding orbital construction.
- **Confidence: MEDIUM-HIGH.** Depends on the basis-closable hypothesis. Alternative interpretation: max_n=3 also provides more diffuse orbitals on H-side that further dominate the bonding manifold, leaving the wall unchanged. CLAUDE.md §2 NaH n_max=3 v2.1.1 was tested at BARE balanced and overattracted; this prediction tests at W1c-screened + multi-zeta which is structurally different — should be more favorable.

#### P2: Binding energy within 2× of experimental
- **Statement:** If P1 holds, the binding energy $D_e \in [0.0375, 0.150]$ Ha (factor-2 around experimental NaH $D_e \approx 0.075$ Ha).
- **Falsifier:** $D_e$ (internal min to dissociation) outside this range. $D_e < 0.0375$ Ha indicates incomplete bonding mechanism; $D_e > 0.150$ Ha indicates residual over-attraction.
- **Underlying classification:** Sub-layer 2 endomorphism content is mitigable via FrozenCore input; under combined architecture should give physical-scale energetics.
- **Confidence: MEDIUM.** Even if P1 holds, the absolute $D_e$ value depends on the cross-V_ne kernel-shape accuracy (Track 3's named P2 target, addressed by neither this sprint nor F1 max_n=3 sprint). Factor-2 tolerance is conservative.

#### P3: Multi-zeta differential persists at max_n=3
- **Statement:** $|E(\text{W1c+mz}) - E(\text{W1c})|$ at $R_\text{eq}$ is in $[0.02, 0.30]$ Ha (similar to max_n=2's +0.056 Ha at R=3.5).
- **Falsifier:** $|E(\text{W1c+mz}) - E(\text{W1c})|$ at $R_\text{eq}$ < 0.005 Ha (mz absorbed by basis expansion, suggesting basis expansion provides Na-side flexibility that multi-zeta was capturing) or > 0.30 Ha (anomalous basis-set scaling).
- **Underlying classification:** Sub-layer 2 endomorphism content is fixed; cross-shift content (cross-V_ne) provides visibility; larger basis adds matrix elements but per-cell endomorphism content is bounded.
- **Confidence: HIGH.** The endomorphism-vs-cross-shift partition predicts this composition behavior cleanly. Violation would require revising the classification.

### 5.3 What the predictions test together

If **P1 + P2 + P3 all hold**: the M-Z partition extends cleanly to the W1c sub-layer hierarchy with the basis-closable refinement; the framework's chemistry-side endomorphism mitigability is empirically validated; the W1c wall is closable in production at max_n=3.

If **P1 fails (no internal min)**: the basis-closable classification of sub-layer 3 is wrong; the wall is structurally deeper. Pivot to Track 3 named P2 (cross-V_ne kernel-shape substitution) which addresses a distinct cross-shift mechanism.

If **P1 holds but P2 fails (D_e wildly off)**: the bonding mechanism is realized but the cross-V_ne kernel is not accurate enough. Same pivot to Track 3 P2.

If **P3 fails (mz absorbed)**: the endomorphism classification of multi-zeta needs revision — perhaps at larger basis the framework can construct adequate radial-shape representations of Na 3s from hydrogenic-Z=1 primitives without external multi-zeta input. Would imply the chemistry endomorphism wall is *less* structural than thought, which would weaken the LS-8a analog reading.

---

## §6. LS-8a wall analog comparison

### 6.1 The structural parallel

Both the QED-side LS-8a wall and the chemistry-side multi-zeta wall are module endomorphisms in M-Z vocabulary:

| Aspect | QED-side (LS-8a) | Chemistry-side (multi-zeta sub-layer 2) |
|:---|:---|:---|
| Module being deformed | Electron Dirac module $\mathcal{M}_e = L^2(\R^3) \otimes \C^4$ | Na valence sub-module $\mathcal{M}_\text{Na}$ |
| Endomorphism content | $Z_2 - 1$ wavefunction renormalization; $\delta m$ mass shift; multi-loop QED counterterms | Multi-zeta exponents and coefficients for Na 3s/3p radial functions |
| Source of external input | UV-completion theory (perturbative QED expansion + cutoff regulators + counterterm extraction) | FrozenCore radial Schrödinger solver + atomic-physics basis fits (Clementi-Raimondi / BBB93 / Clementi-Roetti) |
| Framework autonomy | NOT autonomously generable from spectral action data | Computable from autonomous atomic-physics solvers (HF on frozen-core potential) |

### 6.2 Is the mitigability difference structural or contingent?

**STRUCTURAL.** The argument:

- **QED-side mitigability problem.** Renormalization counterterms parameterize divergent loop integrals. The framework's spectral action (Connes-Chamseddine) does not regulate these divergences autonomously; counterterms come from a deeper UV-completion structure (e.g., asymptotic safety, string theory UV-completion, lattice QFT cutoff). As long as QED is treated as an effective field theory with a UV cutoff, the counterterm values are not autonomously fixed by the spectral action alone. This is structural to the framework: the spectral action is a *low-energy* construction, not a UV-completion.

- **Chemistry-side mitigability solution.** Atomic structure of the frozen core is a closed-form computable eigenvalue problem on a well-posed potential. For Z up to ~80 (Pb), the FrozenCore + Schrödinger or Dirac-Coulomb radial solvers give converged eigenstates to arbitrary precision. The multi-zeta fit is a finite-dimensional approximation that can be systematically improved by adding more primitives (K=5 mixed-Slater-n already gives overlap > 0.999999). The framework's chemistry-side endomorphism content is in principle computable from autonomous atomic-physics methods.

The distinction: QED counterterms parameterize an *open-ended* UV physics problem (the SM is not UV-complete in the strict sense; the Landau pole is a real feature of QED at very high energies). Atomic structure parameterizes a *closed-form* eigenvalue problem on a frozen-core potential. The two problems have categorically different mathematical structure, regardless of current technology.

### 6.3 Implication for the structural-skeleton-scope statement

The framework's structural-skeleton-scope statement (CLAUDE.md §1.7) says: "Framework determines selection rules, transcendental signatures, scaling laws, divergence structure, factorization theorems, and upper bounds; calibration data (parameter values, renormalization counterterms, gauge lower bounds, inner-factor selection) is empirical input."

The chemistry-side endomorphism wall is *consistent with* this statement but with a softer interpretation: the calibration data is computable from atomic-physics methods (which themselves use the same framework's atomic FCI machinery on isolated atoms). The chemistry-side wall is therefore *internally mitigable* via the framework's own infrastructure — it does not need external compositions outside the GeoVac architecture.

The QED-side wall is structurally external: counterterms require UV-completion physics that is genuinely outside any extension of the spectral action.

**Net:** chemistry-side endomorphism content (sub-layer 2 multi-zeta) is *intrinsically external* in the sense that it comes from a frozen-core solver, but *framework-internally computable* in the sense that the frozen-core solver is itself a GeoVac construction. QED-side endomorphism content (LS-8a counterterms) is *intrinsically external* in a strictly stronger sense: no GeoVac construction generates it.

### 6.4 Where this leaves the multi-focal-composition wall

The multi-focal-composition wall (CLAUDE.md §1.7) bundled five Layer-2 walls (Sprint H1, LS-8a, HF-3, HF-4, HF-5) as "framework couples discrete labels cleanly, doesn't compose multiple Fock-style projections." Under the M-Z partition + W1c sub-layer hierarchy + LS-8a analog, the picture refines:

- **HF-3, HF-4, HF-5 (and W1c sub-layer 1 = PK):** bimodule cross-shifts, handled by the framework.
- **LS-8a, Sprint H1 Yukawa (QED-side endomorphisms):** module endomorphisms requiring external UV-completion input. Structurally not mitigable.
- **W1c sub-layer 2 (multi-zeta, chemistry-side endomorphism):** module endomorphism requiring external atomic-physics input. Structurally mitigable via FrozenCore.
- **W1c sub-layer 3 (orbital-pair flexibility, chemistry-side basis-closable cross-shift):** bimodule cross-shift, handled at sufficient basis size. Sub-classification of cross-shift.

The five-observable multi-focal-composition wall splits into three distinct mechanism types, not two. The chemistry-side adds the basis-closable cross-shift class (sub-layer 3) and the mitigable-endomorphism class (sub-layer 2) — neither of which the original M-Z partition (derived from QED-side) named explicitly.

This is genuine framework-level new content from the bridge analysis.

---

## §7. Verdict and recommendations

### 7.1 Verdict

**The M-Z partition extends cleanly to the W1c sub-layer hierarchy: PARTITION-EXTENDS-CLEANLY.**

All three F1 W1c sub-layers classify under M-Z's cross-shift / endomorphism partition:

1. PK (sub-layer 1) = bimodule cross-shift ✓
2. Multi-zeta (sub-layer 2) = module endomorphism ✓ (chemistry-side, mitigable)
3. Orbital-pair flexibility (sub-layer 3) = bimodule cross-shift with **basis-closable** sub-qualifier ✓

The sub-classification "basis-closable cross-shift" is the substantive new content beyond M-Z's original two-bucket framing. It emerges from FCI variational state analysis on a multi-determinant bimodule — a feature M-Z's bound-state QED context did not have.

The structural-skeleton-scope statement is preserved and slightly refined: the framework's chemistry-side endomorphism content (sub-layer 2) is intrinsically external (FrozenCore input) but framework-internally computable, in contrast to QED-side endomorphism (LS-8a counterterms) which is strictly externally specified.

### 7.2 Are the predictions well-grounded?

| Prediction | Grounding | Confidence |
|:---|:---|:---|
| P1 (internal min emerges at max_n=3) | Basis-closable cross-shift classification; max_n=3 adds H 2s/2p, Na 3p/3d providing dimensional richness | MEDIUM-HIGH |
| P2 (D_e within 2× of experimental) | Chemistry endomorphism mitigability; combined architecture should give physical-scale energetics | MEDIUM |
| P3 (mz differential persists in [0.02, 0.30] Ha) | Endomorphism-vs-cross-shift partition predicts bounded per-cell endomorphism content | HIGH |

All three predictions are structurally falsifiable, with named numerical thresholds. The most discriminating prediction is P1 (binary: internal min or not); P2 and P3 are quantitative refinements that test the strength of the classification.

### 7.3 Should M-Z partition be promoted in CLAUDE.md and Paper 34?

**CONDITIONAL on F1 max_n=3 outcome.**

- **If P1 confirmed**: the M-Z partition has predicted a specific physics outcome on a system class (alkali hydrides) it was not originally designed for. This is the structural test passed; the partition deserves promotion as a general framework tool. Recommended edits:
  - CLAUDE.md §1.7: refine the multi-focal-composition wall taxonomy with the three-bucket sub-classification (cross-shift handled / cross-shift basis-closable / endomorphism mitigable-or-not).
  - Paper 34 §V.D.7: add the running-catalogue entry "Multi-focal wall splits into bimodule cross-shifts (some basis-closable) vs module endomorphisms (some mitigable). Per Sub-sprint M-Z + F1 + W1c bridge, 2026-05-23."
  - Paper 19 §sec:w1c_residual: extend with the three-sub-layer hierarchy + M-Z classifications + F1 max_n=3 binding closure.

- **If P1 refuted**: keep the partition as local-vocabulary refinement (bound-state QED context where it was originally derived). Do NOT promote as general framework tool. The chemistry-side basis-closable sub-classification was incorrect; sub-layer 3 is a genuine endomorphism (basis truncation alone cannot generate the bonding combination from within), and the wall is structurally deeper.

### 7.4 Recommended sprint sequence

**Step 1 (HIGH PRIORITY, 3-5 days):** F1 max_n=3 NaH FCI with W1c + multi-zeta architecture. Tests P1/P2/P3 simultaneously. This is the operational test the bridge analysis points to.

**Step 2 (CONDITIONAL on P1 confirmed):** Write up M-Z partition + three-bucket refinement (cross-shift handled / cross-shift basis-closable / endomorphism mitigable-or-not) in Paper 19 + Paper 34 + CLAUDE.md. The two-bucket → three-bucket refinement is the genuinely new structural content; locking it in writing prevents memory decay.

**Step 3 (CONDITIONAL on P1 refuted):** Pivot to Track 3's named target P2 (cross-V_ne kernel-shape substitution), which addresses a distinct cross-shift mechanism on the partner-side basis. The Track 3 May-9 diagnostic predicted -0.674 Ha differential at R=2.5 — substantial enough that the FCI eigenspectrum could be shifted into a bonding regime even at limited basis.

**Step 4 (PARALLEL, independent):** Extend the bridge analysis to other chemistry-side walls (W1c on MgH₂, HCl, H₂S where Z-decreasing residual gives different empirical magnitudes; lone pair coupling at H₂O with Z_eff > 4 where Slater integrals are unphysical). Tests whether the three-bucket M-Z partition refinement generalizes across the chemistry catalogue.

### 7.5 Honest scope

**This memo:**
- IS a structural analysis applying the M-Z partition (Track M-Z, 2026-05-23) to the W1c sub-layer hierarchy (F1, 2026-05-23).
- IS NOT a numerical verification of the predictions P1/P2/P3 — those require F1 max_n=3 execution.
- DOES introduce one substantive new structural element: the "basis-closable cross-shift" sub-classification of sub-layer 3.
- DOES NOT modify production code, papers, or CLAUDE.md. Recommendations gated on F1 max_n=3 outcome.

**Confidence:**
- HIGH on the bimodule structure (A_L, A_R, M) for NaH balanced.
- HIGH on sub-layer 1 PK classification as cross-shift (mechanism explicit).
- HIGH on sub-layer 2 multi-zeta classification as endomorphism (mutual-exclusivity behavior matches signature).
- MEDIUM-HIGH on sub-layer 3 classification as basis-closable cross-shift (depends on bonding-vs-antibonding interpretation; alternative interpretation as genuine endomorphism is testable).
- MEDIUM-HIGH on prediction P1 (internal min emerges).
- MEDIUM on prediction P2 (binding within 2× of experimental).
- HIGH on prediction P3 (mz differential persists).
- HIGH on the LS-8a analog mitigability argument (structural difference between open-ended UV physics and closed-form atomic eigenvalues).

The bridge analysis is well-grounded as a structural framework application; the numerical test (F1 max_n=3) is the natural next step to validate or refute the bridge predictions.

---

## §8. Files

### Created (this sprint)
- `debug/sprint_w1c_mz_partition_analysis_memo.md` (this memo, ~3800 words)
- `debug/data/sprint_w1c_mz_partition_predictions.json` (structured predictions for programmatic checking)

### Cross-references (read-only)
- `debug/sprint_modular_propinquity_mZ_bethe_log_memo.md` (M-Z primary, partition source)
- `debug/sprint_modular_propinquity_synthesis_memo.md` (modular synthesis §3 with partition statement)
- `debug/sprint_f1_p1p2_combined_test_memo.md` (F1 primary, W1c sub-layer hierarchy)
- `debug/sprint_modular_alpha_arc_synthesis_memo.md` (α arc synthesis with original Layer 3 sketch)
- `debug/sprint_alpha_1_diagnostic_memo.md` (d_R/d_L = 3.23 bimodule diagnostic)
- `debug/sprint_alpha_2_multizeta_memo.md` (multi-zeta architecture, mixed Slater-n)
- `geovac/balanced_coupled.py` (NaH builder, bimodule arena)
- `geovac/multi_zeta_orbitals.py` (Z=11 entry)
- `geovac/cross_center_screened_vne.py` (W1c × mz unified path)
- `geovac/phillips_kleinman_cross_center.py` (PK implementation)
- `geovac/neon_core.py` (FrozenCore Z_eff(r))
- CLAUDE.md §1.7 multi-focal-composition wall taxonomy
- CLAUDE.md §2 v2.59+ modular arc + F1 entries
- Paper 32 §VIII.C (Sprint H1 verdict — LS-8a wall context)
- Paper 19 §sec:w1c_residual (W1c chemistry context)

### Not modified
- Any paper `.tex` files (recommendations only; PI to decide on conditional promotion based on F1 max_n=3).
- CLAUDE.md (recommendations only; PI to decide).
- Any production code (this is structural analysis).
- Any tests.

---

**End of W1c × M-Z partition bridge analysis memo.**
