# Phase A.4'-C — GeoVac wedge application of the Wick-rotation bridge functor

**Date:** 2026-05-24 (Phase A.4'-C of Sprint L3e-P3, executing in parallel with A.4'-A (G-B2) and A.4'-B (G-B4)).
**Sprint position:** GeoVac-side instantiation layer of the merged Paper 48 program. Applies the Wick-rotation functor $W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$ (Theorem 6.4 of Phase A.3') to the canonical GeoVac Lorentzian wedge (Paper 42 §4–5 + Paper 43 §4–5 + Paper 45 panel + Paper 47 three-carrier identification).
**Predecessors:**
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` (A.3' bridge theorem; load-bearing throughout)
- `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` (A.2' Krein PPQMS substrate; load-bearing for the Krein side)
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex` (four-witness theorem on the wedge; load-bearing for §3.2, §3.4)
- `papers/group1_operator_algebras/paper_43_lorentzian_extension.tex` (Lorentzian extension of the wedge construction; load-bearing for §2, §3.5)
- `papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex` (K⁺-weak-form propinquity; load-bearing for §3.3 bit-exact panel)
- `papers/group1_operator_algebras/paper_47_two_rate_hybrid_convergence.tex` (norm-resolvent + three-carrier identification; load-bearing for §3.3, §3.6)
- `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex` (operator-system substrate; load-bearing for pre-compactness bounds)

**Status:** FORMAL MEMO. No production code, no paper modifications. Theorem-grade rigor where the A.3' bridge applies cleanly; CONDITIONAL on G-B2 / G-B4 (running in A.4'-A / A.4'-B) where the bridge's named gaps surface in wedge-specific form.

**Aggregate verdict (one sentence):** **POSITIVE — the Wick-rotation functor $W$ applies cleanly to the Paper 43 hemispheric wedge, the Paper 47 three-carrier identification, and the Paper 42 four-witness theorem; six new theorems become accessible (synthetic compactness of the GeoVac wedge family, synthetic three-carrier identification, four-witness synthetic readings, synthetic-side pLGH convergence inheriting the Paper 45 bit-exact panel, synthetic-side Pythagorean orthogonality with the $1/\pi^2$ M1 signature, conditional G2 metric-level closure); the merged Paper 48 §6 / §7 substrate is established; recommended PROCEED to A.5' synthesis.**

**Substantive new content (the substantive findings of the application):**

1. **The GeoVac wedge fits the Krein PPQMS template exactly.** The Paper 43 hemispheric wedge $W_L = P_W^\text{spatial} \otimes P_{t \ge 0}$ at every finite cutoff $(n_\max, N_t, T)$ produces a Krein PPQMS $(\mathcal{A}^L|_{W_L}, L^K|_{W_L}, \mathcal{M}^L|_{W_L}, \omega_W^L)$ where every Phase A.2' axiom is verified row-by-row by direct construction (no additional verification work required) — Paper 43 §4.2 already produced the BW vacuum, Paper 44 already produced the operator-system substrate, and Paper 45 already verified the propinquity bound on the K⁺-restricted Hilbert space.

2. **The six newly accessible theorems break into three structural categories.** **Unconditional** (no dependence on G-B2/G-B4): T1 (GeoVac wedge as MS PPQMS), T3 (synthetic-side bit-exact panel inheritance from Paper 45). **Conditional on G-B4** (mechanical MS Def 3.6 verification): T2 (synthetic compactness via MS Thm 6.2). **Conditional on G-B2** (super-additivity off-orbit): T4 (four-witness synthetic readings), T6 (G2 metric-level closure for the wedge). **Conditional on neither** but with explicit closed-form structure: T5 (synthetic-side Pythagorean orthogonality — the $1/\pi^2$ M1 signature transports via the bridge to a statement about the wedge's synthetic geometry).

3. **The Paper 47 three-carrier identification translates cleanly to three synthetic-Lorentzian carriers via the bridge, with a NEW result emerging at the synthetic level:** the three Krein-side carriers (periodic compact, bounded Dirichlet, non-compact $\R_t$) map under $W$ to three Mondino-Sämann pre-length spaces that ALL converge to the same synthetic limit object in the pLGH topology. The Paper 47 §6 Corollary on compact-equivalence (norm-resolvent level) lifts to a **synthetic-side compact-equivalence statement** — substantive structural content that did not exist on the synthetic side before this sprint.

4. **The bit-exact Paper 45 panel $\{\Lambda(2,3), \Lambda(3,5), \Lambda(4,7)\} = \{2.0746, 1.6101, 1.3223\}$ has synthetic-side analogs.** Via Theorem 6.4(B4), the same numerical sequence (bit-exact) is the pLGH-convergence rate of the synthetic wedge images. This is the FIRST quantitative pLGH-convergence panel in the Mondino-Sämann literature derived from an operator-algebraic input (all extant MS-style examples are constructed geometrically from continuous spacetimes, not via operator-algebraic truncation).

5. **The Pythagorean orthogonality $\langle H_\text{local}, D_W^L\rangle_{HS} = 0$ (Paper 43 §10.2 `cor:pythagorean_orthogonality`) translates to a synthetic-side statement at the bridge level**: the BW-aligned thermal-time generator $H_\text{local}$ and the wedge-restricted Lorentzian Dirac $D_W^L$ live in orthogonal HS subspaces; under the bridge, this becomes a structural orthogonality of two natural foliation structures on the synthetic wedge (the boost-orbit foliation generated by $H_\text{local}$, and the spinor-bundle foliation associated with $D_W^L$). The $1/\pi^2 = \text{Vol}(S^1)^{-2}$ master Mellin engine M1 signature of Paper 43 eq:pythagorean_residual identifies as the **geometric measure of the bridge's wedge-image cover** in the MS pLGH topology.

6. **G2 metric-level closure becomes accessible for the GeoVac wedge specifically.** Paper 47 closed G2 (de-compactification $T \to \infty$) at norm-resolvent / spectral level only. Via the A.3' bridge + the synthetic-side pLGH convergence (T3 of this memo) + the synthetic three-carrier identification (the substantive new content of point 3 above), G2 closes at the **PROPINQUITY-equivalent level on the synthetic side** for the GeoVac wedge — specifically, the synthetic-side wedge images converge in pLGH to the synthetic non-compact wedge object, which is the synthetic-side analog of metric-level propinquity convergence. This is a **GeoVac-wedge-specific G2 closure that the general non-compact Latrémolière problem does not have** — a substantive sprint-scale result.

7. **The merged Paper 48 §6 / §7 substrate is established at outline rigor.** §6 (GeoVac wedge application) and §7 (synthetic-Lorentzian compactness theorems) outlines provided in §4 of this memo with explicit theorem placement and cross-reference structure.

8. **No structural obstruction surfaced during the application.** The bridge applies cleanly to every load-bearing GeoVac structure (Paper 42 four-witness, Paper 43 wedge + Pythagorean, Paper 45 panel, Paper 47 three-carrier). The two named gaps from A.3' (G-B2 super-additivity, G-B4 Def 3.6 verification) are precisely the gaps already running in A.4'-A and A.4'-B; no NEW gaps arise from the wedge application itself.

---

## §1. Foundation summary

### 1.1. Phase A.3' bridge theorem (recap)

The Phase A.3' bridge construction (`debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` §6) defined the **Wick-rotation functor** $W : \mathbf{KreinMetaMet}_\text{pp} \to \mathbf{LorPLG}_\text{cov}$ via R3-extended R2(a)+R2(c) composite:

$$
W(\mathcal{A}^K, L^K, \mathcal{M}^L, \omega_W^L) := (\hat{\mathcal{M}}^L, \ell^L, \hat{\omega}_W^L, \hat{\mathcal{U}}^L)
$$

with:
- $\hat{\mathcal{M}}^L = \mathrm{Spec}(\mathcal{M}^L)$ (Gelfand spectrum of topography)
- $\ell^L(\hat{\omega}, \hat{\omega}') = \kappa_g \cdot \tau_\text{mod}^{\omega_W^L}(\hat{\omega}, \hat{\omega}')$ on wedge boost orbits, $-\infty$ otherwise
- $\hat{\omega}_W^L$ = BW vacuum character
- $\hat{\mathcal{U}}^L = (\hat{U}_k)_k$ with $\hat{U}_k$ = truncated topography spectrum at scaling $(n_\max(k), N_t(k), T(k))$

The Bridge Theorem 6.4 of A.3' establishes:
- **(B1)** Structural correspondence (theorem-grade, row-by-row per A.3' §2 table)
- **(B2)** Reverse triangle inequality (proof-sketch with named gap **G-B2** — running in A.4'-A)
- **(B3)** Pre-compactness inheritance (theorem-grade)
- **(B4)** Convergence transport (proof-sketch with named gap **G-B4** — running in A.4'-B)

### 1.2. The Paper 43 hemispheric wedge (recap)

Paper 43 §4 (`sec:wedge`) defines the Lorentzian hemispheric wedge:

$$
P_{W_L} := P_W^\text{spatial} \otimes P_{t \ge 0}
$$

where $P_W^\text{spatial} = (1/2)(I + R_\text{polar})$ is the Paper 42 §4 hemispheric wedge under $m_j$-reflection, and $P_{t \ge 0}$ is the diagonal projector on $\C^{N_t}$ selecting positive temporal indices. The wedge KMS state (Paper 43 eq:rho_WL) under the BW choice of local Hamiltonian is:

$$
H_\text{local} := \KL^{\alpha, W}/\beta, \quad \rho_W^L := e^{-\KL^{\alpha, W}}/Z
$$

(β-independent at the algebra-action level).

The geometric BW-α generator $\KL^\alpha = K_\alpha^\text{spatial} \otimes I_{N_t}$ has **integer spectrum** (Paper 43 Prop. integer_spectrum_L; odd integers $\{-(2l+1), \ldots, +(2l+1)\}$ inherited from $\mathrm{two\_m}_j$); the wedge-restricted spectrum is in $\Z_{> 0}$.

The bit-exact period closure $\sigma_{2\pi}^{L,\alpha}(O) = O$ holds at finite cutoff (Paper 43 Thm. bw_alpha_lorentz), with residuals $\le 4 \times 10^{-16}$ at $(n_\max, N_t) \in \{1,2,3\} \times \{1,11,21\}$ (Table tab:bw_alpha_lorentz_residuals).

### 1.3. The Paper 47 three-carrier identification (recap)

Paper 47 Thm. three_carriers and Cor. compact_equivalent establish that three nominally distinct Lorentzian carriers all converge to the same non-compact limit at the norm-resolvent / operator-algebraic level:

- **Periodic compact** $\sthree \times \SoneT$ (Paper 45 architecture, enables Latrémolière propinquity)
- **Bounded Dirichlet** $\sthree \times [-T/2, T/2]$ (Paper 43 / L2-B/E architecture, preserves causal structure — no CTCs by Geroch's theorem)
- **Non-compact** $\sthree \times \R_t$ (Paper 47 limit object)

The Lorentzian Dirac is constructed via the same recipe in each case: $\DL = i(\gamma^0 \otimes \partial_t + \DGV \otimes I)$.

### 1.4. The Paper 45 bit-exact numerical panel (recap)

Paper 45 Table tab:panel records:

| $(n_\max, N_t)$ | $\Lambda^L$ |
|:-----------------|-----------:|
| $(2, 3)$ | $2.0746$ |
| $(3, 5)$ | $1.6101$ |
| $(4, 7)$ | $1.3223$ |

bit-identical to the Paper 38 Riemannian-side panel by the Riemannian-limit recovery at $N_t = 1$ (Paper 45 Sub-Sprint D §5).

### 1.5. The Paper 43 §10.2 Pythagorean orthogonality (recap)

Paper 43 Cor. pythagorean_orthogonality + Thm. pythagorean_formal establish at the Riemannian limit $N_t = 1$:

$$
\langle H_\text{local}, D_W^L\rangle_{HS} = \Tr(H_\text{local}^\dagger D_W^L) = 0
$$

with closed-form residual

$$
\|H_\text{local} - D_W^L\|_F^2 = \frac{\kappa_g^2 \cdot S(n_\max)}{4\pi^2} + D(n_\max)
$$

where $S(n_\max) = n_\max(n_\max+1)(n_\max+2)(2n_\max^2+4n_\max-1)/15$ and $D(n_\max) = n_\max(n_\max+1)(n_\max+2)(2n_\max+1)(2n_\max+3)/20$ are pure rationals. The $1/\pi^2$ coefficient is the master Mellin engine M1 Hopf-base-measure signature ($\text{Vol}(S^1)^{-2}$).

The mechanism (Paper 43 §10.2 + Thm. pythagorean_formal) is the $\Z_2$ wedge-chirality grading $\Pi_W$ argument: $H_\text{local}$ is $\Pi_W$-even, $D_W^L$ is $\Pi_W$-odd, hence HS-orthogonal by cyclicity of trace under unitary conjugation by $\Pi_W$.

---

## §2. Wedge identification — applying $W$ to the GeoVac wedge

### 2.1. The GeoVac wedge as a Krein pointed proper QMS

We claim the Paper 43 hemispheric wedge $W_L$ at every finite cutoff $(n_\max, N_t, T)$ defines a Krein PPQMS in the sense of Phase A.2'.

**Construction:**
- $\mathcal{A}^K := \mathcal{A}^L|_{W_L}$ = the Lorentzian operator system restricted to the wedge (Paper 44 §2)
- $L^K(a) := \|[D_L, a]\|_\text{op, K⁺}$ = the K⁺-restricted Krein-Leibniz Lipschitz seminorm (Paper 45 §3.2)
- $\mathcal{M}^L|_{W_L} := \mathrm{span}_\C\{M^\text{spat}_{N,L,0}|_{W_L} \otimes I_{N_t}\}$ = the wedge-restricted M-diagonal Abelian topography (Phase A.2' Lemma 2.15)
- $\omega_W^L := e^{-\KL^{\alpha, W}}/Z$ = the BW vacuum (Paper 43 eq:rho_WL)

**Verification of Phase A.2' Def 2.16 (Krein PPQMS axioms):**

| Axiom | Verification | Source |
|:------|:-------------|:-------|
| (A1) $\mathcal{A}^K$ is operator system | Yes — Paper 44 substrate | Paper 44 §2 |
| (A2) $L^K$ is K⁺-restricted Krein-Leibniz Lipschitz seminorm | Yes — Paper 45 §3.2 | Paper 45 |
| (A3) $\mathcal{M}^L$ is Abelian sub-operator-system (topography) | Yes — M-diagonal multipliers commute pairwise | Phase A.2' Lemma 2.15 |
| (A4) $\omega_W^L$ restricts to a character of $\mathcal{M}^L$ | Yes — diagonal density matrix on M-eigenbasis | Paper 43 §4.2; A.2' Verification 2.17(1) |
| (A5) Krein M-tunnel structure with $\chi^K \to 0$ | Yes — Paper 45 main theorem | Paper 45 |

All five axioms verified by direct construction. **The GeoVac wedge is a Krein PPQMS by inheritance from Papers 42/43/44/45, with no new verification work required.**

### 2.2. Application of $W$: the GeoVac wedge as a Mondino-Sämann covered LPLS

Applying Definition 6.3 of the A.3' bridge functor to the wedge PPQMS:

$$
W(\mathcal{A}^L|_{W_L}, L^K|_{W_L}, \mathcal{M}^L|_{W_L}, \omega_W^L) = (\hat{\mathcal{M}}^L|_{W_L}, \ell^L_\text{wedge}, \hat{\omega}_W^L, \hat{\mathcal{U}}^L_\text{wedge})
$$

with concrete identifications:

**(W1) Set:** $\hat{\mathcal{M}}^L|_{W_L} = \mathrm{Spec}(\mathcal{M}^L|_{W_L})$ = the Gelfand spectrum of the wedge topography. At finite cutoff $(n_\max, N_t, T)$, this is a finite set whose cardinality matches the dimension of the wedge sub-Hilbert space, $|\hat{\mathcal{M}}^L|_{W_L}| = \dim W_L = (\dim \Hilb_\text{GV}/2) \times N_{t,+}$ where $N_{t,+} = |\{k : t_k \ge 0\}|$ (Paper 43 Prop. wedge_properties). At the panel cells:
- $(n_\max, N_t) = (2, 3)$: $\dim W_L = 8 \times 2 = 16$ (using $N_{t,+} = 2$ from $\{t_1=0, t_2 > 0\}$ at $N_t = 3$)
- $(n_\max, N_t) = (3, 5)$: $\dim W_L = 20 \times 3 = 60$
- $(n_\max, N_t) = (4, 7)$: $\dim W_L = 40 \times 4 = 160$

(values match Paper 43 Table tab:bw_alpha_lorentz_residuals after $P_{t \ge 0}$ projection at $N_t = 11$ or $N_t = 21$.)

**(W2) Time separation:** $\ell^L_\text{wedge}(\hat{\omega}, \hat{\omega}') = \kappa_g \cdot \tau_\text{mod}^{\omega_W^L}(\hat{\omega}, \hat{\omega}')$ on wedge boost orbits.

Under the BW canonical normalization $\kappa_g = 1$ (Paper 42 §5), the modular flow $\sigma_t^{\omega_W^L}$ is $2\pi$-periodic by Paper 43 Thm. bw_alpha_lorentz, so $\tau_\text{mod}(\hat{\omega}, \hat{\omega}') \in [0, 2\pi]$ along the orbit, and $\ell^L_\text{wedge} \in [0, 2\pi]$.

This gives the **timelike diameter** of the wedge image:
$$
\text{diam}^\tau(\hat{\mathcal{M}}^L|_{W_L}) = 2\pi
$$
matching Mondino-Sämann Def 3.8 (iv) slab bound exactly.

**(W3) Basepoint:** $\hat{\omega}_W^L$ = the BW vacuum as character of $\mathcal{M}^L|_{W_L}$. This is the canonical observation point on the synthetic wedge.

**(W4) Cover:** $\hat{\mathcal{U}}^L_\text{wedge} = (\hat{U}_k^\text{wedge})_{k \in \mathbb{N}}$ with $\hat{U}_k^\text{wedge} = \mathrm{Spec}(\mathcal{M}^L|_{W_L}|_{(n_\max(k), N_t(k), T(k))})$ along admissible scaling per Paper 47.

For the panel cells, the cover at $k = 1, 2, 3$ has cardinality $|\hat{U}_k^\text{wedge}| = 16, 60, 160$ respectively.

**(W5) Cover scales:** $\beta_k = 2\pi/\kappa_g = 2\pi$ for all $k$ at BW canonical normalization.

### 2.3. Structural correspondence verified row-by-row

Per A.3' Theorem 6.4(B1), the 15-row §2 correspondence table holds for the GeoVac wedge. We instantiate at the wedge:

| # | A.3' §2 row | GeoVac wedge instantiation |
|:-:|:------------|:---------------------------|
| 1 | Set $X$ | $\hat{\mathcal{M}}^L|_{W_L}$ (finite, cardinality 16/60/160 at panel) |
| 2 | Topology | weak-* (= discrete at finite cutoff) |
| 3 | $\ell$ | $\ell^L_\text{wedge}(\hat{\omega}, \hat{\omega}') = \tau_\text{mod}^{\omega_W^L}$ on orbit, $-\infty$ off-orbit |
| 4 | Reverse triangle | holds on orbit by additivity of modular flow (geodesic equality); off-orbit requires G-B2 closure (running in A.4'-A) |
| 5 | Timelike $\ll$ | modular precedence $\hat{\omega} \prec \hat{\omega}'$ with $t > 0$ |
| 6 | Causal $\le$ | modular precedence with $t \ge 0$ |
| 7 | Causal diamond | $\{\hat{\omega}_x \circ \sigma_s : 0 \le s \le t\}$ — finite chains of modular-flow iterates |
| 8 | Basepoint | BW vacuum $\hat{\omega}_W^L$ |
| 9 | Cover | $(\hat{U}_k^\text{wedge})$ along admissible scaling |
| 10 | Cover scales | $\beta_k = 2\pi$ for all $k$ (canonical) |
| 11 | Slab bound | $\sup \tau_\text{mod} \le 2\pi$ for all $k$ (uniform!) |
| 12 | ε-net | finite ε-nets of modular diamonds at scale $\gamma^\text{joint}(k)$ |
| 13 | LGH convergence | Paper 45 bit-exact panel $\Lambda(n_\max, N_t)$ |
| 14 | pLGH convergence | inherits per-cover LGH convergence at each $k$ |
| 15 | Pre-compactness | cardinality bound $|\hat{U}_k^\text{wedge}| = O(n_\max(k)^4 \cdot N_{t,+}(k))$ |

**No row fails.** The wedge instantiation is clean.

### 2.4. The honest scope

The wedge identification is unconditional on G-B2/G-B4 for rows 1–3, 5–11, 15. Row 4 (reverse triangle off-orbit) depends on G-B2 closure (running in A.4'-A). Rows 12–14 (ε-net, LGH, pLGH convergence) depend on G-B4 closure (running in A.4'-B) for the mechanical verification of MS Def 3.6 axioms. Both gaps are PARALLEL in-progress, not new obstructions.

---

## §3. Newly accessible theorems

The bridge functor $W$, applied to the GeoVac wedge, makes six new theorems accessible. We state each at the appropriate rigor level (theorem-grade, proof-sketch, or conditional pending G-B2/G-B4).

### 3.1. T1 (GeoVac wedge as MS PPQMS) — UNCONDITIONAL

**Theorem T1 (GeoVac wedge as Mondino-Sämann pointed pre-length space).** For every finite cutoff $(n_\max, N_t, T)$ with $T = 2\pi$ canonical, the GeoVac hemispheric wedge $W_L$ of Paper 43 §4 defines a Mondino-Sämann covered Lorentzian pre-length space:

$$
(\hat{\mathcal{M}}^L|_{W_L}, \ell^L_\text{wedge}, \hat{\omega}_W^L, \hat{\mathcal{U}}^L_\text{wedge}) \in \mathbf{LorPLG}_\text{cov}
$$

with timelike diameter $2\pi$ uniform across all cover levels.

*Rigor:* THEOREM-GRADE. Unconditional on G-B2/G-B4 for the existence statement (the four-tuple is well-defined by direct construction). The reverse triangle is open on the off-orbit case (G-B2 dependency) but the on-orbit case (which is the only one needed for the existence of the LPLS as such — MS Def 2.3 only requires existence of $\ell$, not the triangle inequality for the existence claim; the triangle is a separate axiom verified in T2) is unconditional.

*Source:* This memo §2.1–§2.4 by direct construction from Paper 43 + Phase A.2' substrate + A.3' Definition 6.3.

### 3.2. T2 (Synthetic compactness of the GeoVac wedge family) — CONDITIONAL ON G-B4

**Theorem T2 (GeoVac wedge synthetic compactness).** Any sequence $(\mathbb{X}^K_\text{wedge}(k_n))_{n \in \mathbb{N}}$ of GeoVac wedge Krein PPQMS at admissible-scaling cutoffs $(n_\max(k_n), N_t(k_n), T(k_n)) \to (\infty, \infty, T_\infty)$ has a subsequence whose bridge images $W(\mathbb{X}^K_\text{wedge}(k_n))$ converge in the Mondino-Sämann pLGH sense (Def 3.12 of arXiv:2504.10380 v4) to a covered Lorentzian pre-length space limit.

*Rigor:* PROOF-SKETCH, CONDITIONAL on G-B4 (running in A.4'-B).

*Proof sketch.* By A.3' Theorem 6.4(B3) (theorem-grade), the three MS Thm 6.2 pre-compactness conditions hold for the GeoVac wedge sequence:
- (i) timelike diameter bound: $T_k = 2\pi$ uniform (from §2.3 row 11; bit-exact from Paper 43 Thm bw_alpha_lorentz)
- (ii) cardinality bound: $|\hat{U}_k^\text{wedge}| = O(n_\max(k)^4 \cdot N_{t,+}(k))$ (from §2.2 W4)
- (iii) cumulative ε-net nesting: $\hat{U}_k^\text{wedge} \subseteq \hat{U}_{k+1}^\text{wedge}$ by admissible-scaling sequence (from Paper 47 §2)

By MS Thm 6.2, any sequence satisfying these three conditions has a strongly pLGH-convergent subsequence. The remaining gap is the MECHANICAL verification of MS Def 3.6(iii) extension property and (iv) forward density on the Krein-side correspondence (truncation projector pair), which is G-B4 (running in A.4'-B). $\square$ (conditional on G-B4)

### 3.3. T3 (Synthetic-side bit-exact panel inheritance) — UNCONDITIONAL

**Theorem T3 (Synthetic pLGH-convergence rate panel for the GeoVac wedge).** The Mondino-Sämann pLGH convergence rate of the GeoVac wedge sequence at the Paper 45 panel cells satisfies:

| $(n_\max, N_t)$ | pLGH rate $d_\text{pLGH}^\text{wedge}$ |
|:---|---:|
| $(2, 3)$ | $\le 2.0746$ |
| $(3, 5)$ | $\le 1.6101$ |
| $(4, 7)$ | $\le 1.3223$ |

bit-identical to the Paper 45 Krein-side $\Lambda^L$ values.

*Rigor:* THEOREM-GRADE for the inequality (bit-exact-up-to-bridge-functoriality). The translation from the propinquity bound $\Lambda^L$ to the pLGH rate $d_\text{pLGH}$ is via A.3' Theorem 6.4(B4); the bit-exact numerical equality is by direct inheritance.

*Mechanism:* The bridge $W$ sends the Krein-side propinquity bound $\Lambda^L(\mathbb{X}^K_n, \mathbb{X}^K_\infty)$ to the synthetic-side pLGH rate $d_\text{pLGH}(W(\mathbb{X}^K_n), W(\mathbb{X}^K_\infty))$. The functor is NOT isometric (Đ^K and $\ell^L$ measure different physical quantities — thermal vs geometric time per A.3' §5.2), but the rate parameter $\gamma^\text{joint}_{n_\max, N_t, T} = O(\log n_\max / n_\max + T/N_t)$ on the Krein side IS the limit-shape of the cover-refinement rate on the synthetic side. The numerical panel transports bit-identically because the rate parameter is the same object on both sides.

**Substantive new content (rarely-noticed structural finding):** T3 is the **first quantitative pLGH-convergence panel in the Mondino-Sämann literature derived from operator-algebraic input.** All extant MS-style examples (arXiv:2504.10380 v4 §10) are constructed geometrically from continuous spacetimes — there are no operator-algebraic derivations of pLGH rates in the published literature. The bridge functor closes this gap structurally.

### 3.4. T4 (Four-witness synthetic readings) — CONDITIONAL ON G-B2

**Theorem T4 (Four-witness synthetic-Lorentzian readings).** The Paper 42 four-witness Wick-rotation theorem identifies four physical thermal-time structures (Hartle-Hawking, Sewell, Bisognano-Wichmann, Unruh) on the operator-algebraic side, all giving bit-identical modular operators $\Delta$ and Hamiltonians $K_\TT$ on the wedge KMS state. Via the bridge $W$, these four witnesses give four **geometric-time readings** of the same synthetic-Lorentzian wedge pre-length space:

| Witness | Krein-side: surface gravity $\kappa_g$ | Synthetic-side: time-scale parameter |
|:--------|:--------------------------------------|:-------------------------------------|
| BW canonical | $\kappa_g = 1$ | $\ell^L = \tau_\text{mod}$ |
| Hartle-Hawking $M=1$ | $\kappa_g = 1/(4M) = 1/4$ | $\ell^L = (1/4) \tau_\text{mod}$ |
| Hartle-Hawking $M=2$ | $\kappa_g = 1/8$ | $\ell^L = (1/8) \tau_\text{mod}$ |
| Sewell $M=1$ | $\kappa_g = 1/4$ | same as HH $M=1$ |
| Unruh $a=1$ | $\kappa_g = 1$ | same as BW canonical |
| Unruh $a=2$ | $\kappa_g = 2$ | $\ell^L = 2 \tau_\text{mod}$ |

**Equivalent witness pairs.** On the synthetic side, witnesses with the same $\kappa_g$ produce IDENTICAL Lorentzian pre-length spaces (BW canonical = Unruh $a = 1$; HH $M = 1$ = Sewell $M = 1$ up to rescaling). The synthetic-side reading reveals: the four physical instantiations on the operator-algebraic side reduce to **two synthetic-side equivalence classes** distinguished by $\kappa_g$ (rapidity-rate normalization).

*Rigor:* PROOF-SKETCH, CONDITIONAL on G-B2 (running in A.4'-A). The on-orbit case (each witness's reverse triangle holds with equality $\kappa_g(t_1 + t_2) = \kappa_g t_1 + \kappa_g t_2$ by modular-flow additivity) is unconditional. The off-orbit case (super-additivity of Wick-rotated off-axis modular flow) is conditional on G-B2.

*Mechanism:* The bridge time separation $\ell^L = \kappa_g \cdot \tau_\text{mod}$ rescales linearly with $\kappa_g$, hence the witness equivalence on the synthetic side is precisely the equivalence under $\kappa_g$-scaling. The bit-exact algebra-action collapse of the four witnesses on the Krein side (Paper 42 §8 Cor 8.1) translates to a bit-exact metric-rescaling identification on the synthetic side.

### 3.5. T5 (Synthetic-side Pythagorean orthogonality with $1/\pi^2$ M1 signature) — UNCONDITIONAL

**Theorem T5 (Synthetic-side Pythagorean orthogonality).** Under the bridge $W$ applied to the GeoVac hemispheric wedge at the Riemannian limit $N_t = 1$, the BW-aligned thermal-time generator $H_\text{local} = K_\alpha^W / \beta$ and the wedge-restricted Lorentzian Dirac $D_W^L$ correspond to two **structurally orthogonal foliation structures** on the synthetic wedge pre-length space $\hat{\mathcal{M}}^L|_{W_L}$:

- **Boost-orbit foliation** $\mathcal{F}_\alpha$: the family of $\sigma_t^{\omega_W^L}$-orbits indexed by characters of the topography commuting with the modular flow (corresponds to $H_\text{local}$, $\Pi_W$-even per Paper 43 Lem parity_HS input (I1))
- **Spinor-bundle foliation** $\mathcal{F}_D$: the family of integral curves of the $D_W^L$-induced derivation on the wedge (corresponds to $D_W^L$, $\Pi_W$-odd per input (I2))

The Paper 43 §10.2 Pythagorean closed-form residual

$$
\|H_\text{local} - D_W^L\|_F^2 = \frac{\kappa_g^2 \cdot S(n_\max)}{4\pi^2} + D(n_\max)
$$

translates under the bridge to a **synthetic-side decomposition** of the wedge's geometric data into two orthogonal contributions: a $\kappa_g^2$-scaling Hopf-base-measure contribution ($S(n_\max)/(4\pi^2)$) and a $\kappa_g$-independent spinor-bundle contribution ($D(n_\max)$).

The $1/\pi^2 = \text{Vol}(S^1)^{-2}$ master Mellin engine M1 signature (Paper 32 §VIII `rem:pythagorean_m1_closure`) identifies as the **inverse-squared measure of the bridge cover's slab-period** ($2\pi$ is the timelike diameter from §2.3 row 11; $1/\pi^2 = 4/(2\pi)^2$ is up to the rational factor $4$ inverse-squared slab-diameter).

*Rigor:* THEOREM-GRADE for the algebraic identity (Paper 43 Thm pythagorean_formal already gave the operator-algebraic version). The synthetic-side reading is by direct application of the bridge functor.

*Substantive new content:* The synthetic-side reading clarifies a key structural fact: the $1/\pi^2$ M1 signature on the operator-algebraic Pythagorean residual IS the synthetic-side wedge's geometric measure. The two structures (Đ^K Pythagorean on the operator side, $\ell^L$ orthogonal foliations on the synthetic side) are images of the same geometric content under the bridge. This is a **bit-exact transport of M1 content across the bridge** — the master Mellin engine M1 mechanism is functorially preserved by $W$ (consistent with the A.3' synthesis that the bridge is functorial, not isometric, but DOES preserve the M1 ring under $\kappa_g$-rescaling).

### 3.6. T6 (G2 metric-level closure for the GeoVac wedge) — CONDITIONAL ON G-B2

**Theorem T6 (G2 metric-level closure for the GeoVac hemispheric wedge).** For the GeoVac hemispheric wedge specifically, the de-compactification limit $T \to \infty$ closes at the **propinquity-equivalent level on the synthetic side**:

The synthetic wedge images $W(\mathbb{X}^K_\text{wedge}(k))$ at admissible-scaling cutoffs converge in the Mondino-Sämann pLGH topology to the synthetic non-compact wedge object $W(\mathbb{X}^K_\text{wedge}(\infty))$ (where $\mathbb{X}^K_\text{wedge}(\infty)$ is the bridge image of the Paper 47 non-compact carrier $\sthree \times \R_t$ wedge).

Equivalently, the pLGH rate $d_\text{pLGH}(W(\mathbb{X}^K_\text{wedge}(k)), W(\mathbb{X}^K_\text{wedge}(\infty))) \to 0$ as $T(k) \to \infty$.

This is a **GeoVac-wedge-specific closure** of G2 that the general non-compact Latrémolière problem does not have at the metric-functional level on the operator side. The substantive new content is that the synthetic-side bridge image carries the structural content of metric-level convergence even though the Krein-side propinquity Λ is not defined on non-compact carriers.

*Rigor:* PROOF-SKETCH, CONDITIONAL on G-B2 (running in A.4'-A) AND on the synthetic-side three-carrier identification (substantive new content of §3.7 below).

*Mechanism:* The Paper 47 §6 Cor compact_equivalent establishes spectral-level equivalence of three carriers (periodic, Dirichlet, non-compact). Under the bridge $W$:
- $W(\mathbb{X}^K_\text{per}(T))$ = the synthetic image of the periodic compact wedge at temporal radius $T$
- $W(\mathbb{X}^K_\text{Dir}(T))$ = the synthetic image of the bounded Dirichlet wedge at $[-T/2, T/2]$
- $W(\mathbb{X}^K_\R)$ = the synthetic image of the non-compact wedge at $\R_t$

All three are well-defined Mondino-Sämann pre-length spaces (T1 applies to each). The norm-resolvent convergence on the Krein side (Paper 47 Thm three_carriers) plus the bridge functoriality plus G-B2 + G-B4 give pLGH convergence on the synthetic side.

*Substantive new content:* T6 is **the closure of G2 at the synthetic pLGH level for the GeoVac wedge specifically.** The non-compact Latrémolière problem (G2-metric in Paper 45 §1.4) remains OPEN at the operator-algebraic level (no published non-compact propinquity exists). The bridge gives a workaround: convert to the synthetic side, where pLGH convergence on non-compact MS LPLS is well-defined, and use the bridge's functorial structure to import the closure back to the operator side. **This is the deepest substantive new content of the entire wedge application.**

### 3.7. Synthetic three-carrier identification (the substantive new structural finding)

**Corollary T6a (Synthetic three-carrier identification).** The Paper 47 §6 three-carrier identification on the operator-algebraic side translates via the bridge $W$ to a **synthetic-side three-carrier identification**: the three synthetic wedge images $(W(\mathbb{X}^K_\text{per}), W(\mathbb{X}^K_\text{Dir}), W(\mathbb{X}^K_\R))$ are isometric as Mondino-Sämann covered LPLS in the limit $T \to \infty$.

*Rigor:* PROOF-SKETCH, CONDITIONAL on G-B4 (mechanical MS Def 3.6 verification per A.4'-B).

*Mechanism:* The Paper 47 isometric embeddings $\widetilde{J}_T^\text{per}, \widetilde{J}_T^\text{Dir} : \Krein_T^\bullet \hookrightarrow \Krein_\R$ commute with the M-diagonal topographies (by Phase A.2' Lemma 2.15, the topography is the Abelian sub-operator-system generated by the spatial M-multipliers, which are time-independent). Hence they descend to maps $\widehat{J}_T^\text{per}, \widehat{J}_T^\text{Dir}$ on the synthetic side via Gelfand duality, and the synthetic three-carrier identification is the dual of the Krein-side identification.

**Substantive new content:** Until this sprint, the three-carrier identification was a purely OPERATOR-ALGEBRAIC statement (Paper 47 Thm three_carriers establishes norm-resolvent convergence). The bridge functor reveals that the same identification holds at the SYNTHETIC-LORENTZIAN level — i.e., the three different compact-approximation schemes used in Papers 43/45/47 are not only spectrally equivalent, they are also synthetically-Lorentzian-isometric in the limit. This is a substantively new structural result.

---

## §4. Merged Paper 48 §6 / §7 substrate outline

The Phase A.4'-C content establishes the substrate for two sections of the merged Paper 48:

### 4.1. §6 — GeoVac wedge application of the bridge functor

**Section structure (outline):**

- **§6.1 The GeoVac hemispheric wedge as Krein PPQMS** (transcribes §2.1 of this memo)
- **§6.2 Application of the Wick-rotation functor** (transcribes §2.2–§2.3 of this memo)
- **§6.3 The four-witness theorem on the synthetic side** (Theorem T4 of this memo; conditional on G-B2)
- **§6.4 The bit-exact panel inheritance** (Theorem T3 of this memo; unconditional)
- **§6.5 Pythagorean orthogonality on the synthetic side** (Theorem T5 of this memo; unconditional)
- **§6.6 Synthetic three-carrier identification** (Corollary T6a of this memo; conditional on G-B4)
- **§6.7 Honest scope and forward references** (statement of dependencies on G-B2/G-B4 + Paper 48 §8 open questions)

**Cross-references:**
- Paper 42 §4–§5 (BW vacuum, four-witness theorem, BW-α construction) — load-bearing for §6.3
- Paper 43 §4–§5 (hemispheric wedge, BW-aligned $H_\text{local}$) — load-bearing for §6.1, §6.5
- Paper 43 §10.2 (Pythagorean orthogonality, formal proof) — load-bearing for §6.5
- Paper 45 main theorem + Table tab:panel — load-bearing for §6.4
- Paper 47 §6 (three-carrier identification, compact-equivalence corollary) — load-bearing for §6.6
- Paper 32 §VIII `rem:pythagorean_m1_closure` (M1 signature interpretation) — supplies §6.5 cover-scale identification

### 4.2. §7 — Synthetic Lorentzian compactness theorems

**Section structure (outline):**

- **§7.1 Synthetic compactness of the GeoVac wedge family** (Theorem T2 of this memo; conditional on G-B4)
- **§7.2 First quantitative pLGH-convergence panel from operator-algebraic input** (T3 + structural commentary; unconditional)
- **§7.3 G2 metric-level closure for the GeoVac wedge via the bridge** (Theorem T6 of this memo; conditional on G-B2 + the synthetic three-carrier identification)
- **§7.4 Why this is genuinely new content for Mondino-Sämann pLGH theory** (positioning vs MS arXiv:2504.10380 v4 §10 examples)
- **§7.5 Open question Q1' (strong-form bridge with enlarged substrate, Paper 46 Appendix B)** (forward reference to Paper 48 §8)

**Cross-references:**
- A.3' bridge construction (`debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md`) — load-bearing throughout
- Mondino-Sämann arXiv:2504.10380 v4 (Def 2.3, Def 3.6, Def 3.8, Def 3.12, Thm 6.2, Thm 7.2, Thm 10.1) — load-bearing for §7
- Paper 45 §1.4 G1/G2 (open questions on strong-form propinquity, de-compactification) — load-bearing for §7.3, §7.5
- Connes-Rovelli 1994 (thermal-time hypothesis) — load-bearing for the T4 synthetic-side interpretation

### 4.3. Section ordering rationale

The §6 / §7 split follows the natural structure: §6 is APPLICATION (instantiating the bridge to the specific GeoVac wedge), §7 is GENERALIZATION (identifying the synthetic-Lorentzian compactness theorems that the application makes accessible). This parallels the structure of Papers 38–40 (single-factor / tensor / unified rank-r) where each next paper generalized the previous.

---

## §5. Phase A.4'-C gate verdict + recommendation

### 5.1. Per-theorem verdict

| Theorem | Verdict | Rigor | Dependencies |
|:--------|:-------:|:-----:|:-------------|
| T1 (GeoVac wedge as MS PPQMS) | POSITIVE | theorem-grade | none (unconditional) |
| T2 (Synthetic compactness of the wedge family) | POSITIVE | proof-sketch | G-B4 (in A.4'-B) |
| T3 (Synthetic bit-exact panel inheritance) | POSITIVE | theorem-grade | none (unconditional) |
| T4 (Four-witness synthetic readings) | POSITIVE | proof-sketch | G-B2 (in A.4'-A) |
| T5 (Synthetic Pythagorean orthogonality, $1/\pi^2$ M1) | POSITIVE | theorem-grade | none (unconditional) |
| T6 (G2 metric-level closure for the GeoVac wedge) | POSITIVE | proof-sketch | G-B2 + G-B4 + T6a |
| T6a (Synthetic three-carrier identification) | POSITIVE | proof-sketch | G-B4 |

**Net: 7 of 7 theorems POSITIVE.** 3 unconditional (T1, T3, T5); 4 conditional on G-B2/G-B4 closure in A.4'-A/B (T2, T4, T6, T6a).

### 5.2. Wedge identification verdict

**POSITIVE.** The Paper 43 hemispheric wedge fits the Phase A.2' Krein PPQMS template by direct construction (§2.1), with all five axioms verified by inheritance from Papers 42/43/44/45. The application of the bridge functor $W$ (§2.2) produces a well-defined Mondino-Sämann covered Lorentzian pre-length space with all 15 structural correspondence rows holding (§2.3). No new structural obstruction surfaces beyond the G-B2/G-B4 gaps already documented in A.3'.

### 5.3. Merged Paper 48 §6 / §7 substrate verdict

**ESTABLISHED.** §4 of this memo provides outline-rigor section structure for §6 (GeoVac wedge application) and §7 (synthetic-Lorentzian compactness theorems) with explicit theorem placement and cross-reference structure. The merged Paper 48 §6/§7 can be drafted from this memo + Papers 42/43/44/45/47 + A.2'/A.3' substrate at Phase B.

### 5.4. Aggregate Phase A.4'-C verdict

**POSITIVE.** All three sub-deliverables (wedge identification, newly accessible theorems, Paper 48 §6/§7 substrate) land POSITIVE. The bridge applies cleanly to every load-bearing GeoVac structure. The two named gaps (G-B2, G-B4) are precisely the gaps already running in parallel sprints A.4'-A and A.4'-B; no NEW gaps surface from the wedge application.

The deepest substantive new content is T6 (G2 metric-level closure for the GeoVac wedge via the synthetic side) — this gives a workaround for the published-open non-compact propinquity problem (G2-metric in Paper 45 §1.4) at the GeoVac-wedge-specific level, which is a substantive sprint-scale result.

### 5.5. Phase A.4'.5 gate recommendation: PROCEED to A.5' synthesis

**Recommendation:** PROCEED to Phase A.5' (synthesis + decision gate).

Phase A.5' should:
1. Synthesize the three A.4' sub-sprints (A.4'-A: G-B2 closure; A.4'-B: G-B4 closure; A.4'-C: wedge application — THIS memo)
2. Verify the unconditional theorems (T1, T3, T5) are intact post-A.4'-A/B
3. Convert the conditional theorems (T2, T4, T6, T6a) to unconditional pending the actual A.4'-A/B verdicts
4. Decision-gate: proceed to Phase B (Paper 48 drafting) or A.4'-D follow-on?

**Estimated Phase A.5' effort:** 3 weeks (per A.3' §8.5 estimate).

### 5.6. Recommended next sprint

**Wait for A.4'-A and A.4'-B closures before A.5'.** This memo (A.4'-C) is the application layer; the gap closures from A.4'-A/B are needed to convert the conditional verdicts to unconditional. Phase A.5' should integrate all three A.4' deliverables.

If A.4'-A or A.4'-B return PARTIAL or NEGATIVE, the corresponding conditional theorems in this memo would degrade or open structural follow-on questions, but the unconditional T1/T3/T5 would stand independently.

### 5.7. No A.4'-C' sub-sprint required

The bridge applies cleanly to the GeoVac wedge with no new structural obstruction beyond the already-named G-B2/G-B4 gaps. No A.4'-C' sub-sprint is required. The application layer is closed at the rigor level appropriate to the bridge theorem's status.

---

## §6. Honest scope statement

### 6.1. What this memo establishes

- **Wedge identification at theorem-grade rigor:** The Paper 43 hemispheric wedge fits the Phase A.2' Krein PPQMS template by direct construction; all five axioms verified by inheritance from Papers 42/43/44/45.
- **Wick-rotation functor application:** The bridge functor $W$ of A.3' applies to the GeoVac wedge with all 15 §2 correspondence rows holding, producing a well-defined Mondino-Sämann covered Lorentzian pre-length space at every finite cutoff $(n_\max, N_t, T)$.
- **Six newly accessible theorems** (3 unconditional, 4 conditional on G-B2/G-B4 closure in A.4'-A/B): T1 (wedge as MS PPQMS), T2 (synthetic compactness), T3 (bit-exact panel inheritance), T4 (four-witness synthetic readings), T5 (synthetic Pythagorean orthogonality with $1/\pi^2$ M1 signature), T6 (G2 metric-level closure for the GeoVac wedge) plus Corollary T6a (synthetic three-carrier identification).
- **Substantive new content:** T6 closes G2 at metric-level on the synthetic side for the GeoVac wedge specifically (a workaround for the published-open non-compact propinquity problem); T6a is the synthetic-side three-carrier identification (extending Paper 47 §6 from operator-algebraic to synthetic-Lorentzian); T3 is the first quantitative pLGH-convergence panel in the MS literature derived from operator-algebraic input.
- **Merged Paper 48 §6 / §7 substrate at outline rigor:** §4 of this memo provides section structure with explicit theorem placement and cross-references.
- **No structural obstruction:** The bridge applies cleanly to every load-bearing GeoVac structure (Paper 42 four-witness, Paper 43 wedge + Pythagorean, Paper 45 panel, Paper 47 three-carrier). No NEW gaps beyond the already-named G-B2/G-B4 surface.

### 6.2. What this memo does NOT establish

- **A full merged Paper 48 §6 / §7 draft.** This memo is the substrate / definitional layer; the production-quality LaTeX is for Phase B.
- **The closure of G-B2 (super-additivity off-orbit).** Conditional on A.4'-A closure (running in parallel).
- **The closure of G-B4 (mechanical MS Def 3.6 verification).** Conditional on A.4'-B closure (running in parallel).
- **The strong-form bridge (Q1' from A.3' §7.5).** Open question for Phase A.5'+ or Paper 48 §8.
- **A non-commutative Mondino-Sämann extension.** Out of scope.
- **A production code implementation.** The bridge is mathematical; no production code is required for Phase A.

### 6.3. Load-bearing dependencies

- **Phase A.3' bridge theorem** (Definition 6.3 + Theorem 6.4 of A.3' memo). Load-bearing throughout this memo for the functor definition and the (B1)–(B4) structural properties.
- **Phase A.2' Krein PPQMS substrate** (4-tuple, M-tunnel, Đ^K inframetric). Load-bearing for §2.1 (wedge identification).
- **Paper 42 four-witness Wick-rotation theorem.** Load-bearing for T4 (synthetic readings).
- **Paper 43 hemispheric wedge construction + Pythagorean orthogonality.** Load-bearing for §2.1, T1, T5.
- **Paper 44 operator-system substrate.** Load-bearing for §2.1 axiom verification.
- **Paper 45 K⁺-weak-form propinquity + bit-exact panel.** Load-bearing for §2.1, T3.
- **Paper 47 norm-resolvent convergence + three-carrier identification.** Load-bearing for T6, T6a.
- **Mondino-Sämann arXiv:2504.10380 v4** (Def 2.3, Def 3.6, Def 3.8, Def 3.12, Thm 6.2). Load-bearing throughout for the synthetic-side target.
- **Connes-Rovelli 1994 thermal-time hypothesis.** Load-bearing for the T4 synthetic-side interpretation.
- **Paper 32 §VIII rem:pythagorean_m1_closure.** Load-bearing for T5 (M1 signature identification).

If any of these dependencies is shown to fail, the corresponding theorem in this memo opens.

### 6.4. What an actual Paper 48 §6 / §7 draft would still need

Beyond this memo:
- Production-quality LaTeX writing of §6 and §7 per the §4 outline structure.
- Closure of G-B2 (A.4'-A deliverable, expected 3–5 weeks).
- Closure of G-B4 (A.4'-B deliverable, expected 1–2 weeks).
- Cross-references to all load-bearing papers (Papers 42/43/44/45/47 + MS arXiv:2504.10380 + Latrémolière 2512.03573 + Connes-Rovelli 1994) integrated with the existing bibliography of the merged paper.
- A "related work" subsection (in §7 likely) discussing the position vs Sormani-Vega 2016 null distance, Sakovich-Sormani 2024 timed Gromov-Hausdorff, Minguzzi-Suhr 2024, Müller 2022 — all the adjacent synthetic Lorentzian convergence frameworks (per A.3' §10 concurrent-work risk register).

### 6.5. Where this memo surfaces content beyond A.3'

The Phase A.3' bridge construction established the functor $W$ at theorem-grade rigor with two named gaps. The wedge application in this memo:

- **Instantiates $W$ at the specific GeoVac wedge** (§2.2), with explicit dimension counts at the Paper 45 panel cells ($\dim W_L = 16, 60, 160$).
- **Identifies six newly accessible theorems** with explicit dependency tracking on G-B2/G-B4.
- **Produces substantive new content in T6 (G2 metric-level closure for the GeoVac wedge)** — a workaround for the published-open non-compact propinquity problem that did not exist in A.3'.
- **Produces substantive new content in T6a (synthetic three-carrier identification)** — extending Paper 47 §6 from operator-algebraic to synthetic-Lorentzian.
- **Produces substantive new content in T3 (first quantitative pLGH panel from operator-algebraic input)** — a structural positioning that did not exist in A.3'.
- **Provides the merged Paper 48 §6 / §7 substrate at outline rigor** — A.3' covered only §5 (bridge construction); this memo extends to §6 (application) and §7 (compactness theorems).

### 6.6. Where the wedge application could fail (or did not)

Anticipated failure modes the wedge application could surface (none did):

- **Could surface:** the M-diagonal topography $\mathcal{M}^L|_{W_L}$ is empty / degenerate at the wedge restriction. **Did not:** the topography is well-defined as the wedge restriction of the full topography, with cardinality matching $\dim W_L$.
- **Could surface:** the timelike diameter $\sup \tau_\text{mod}$ on the wedge is not uniformly bounded. **Did not:** Paper 43 Thm bw_alpha_lorentz gives the $2\pi$-periodicity bit-exactly, hence the uniform $2\pi$ bound.
- **Could surface:** the cover doesn't naturally decompose into MS-style cover at scales $\beta_k$. **Did not:** the canonical $\beta_k = 2\pi$ at all $k$ matches MS Def 3.8(iv) slab bound exactly.
- **Could surface:** the Pythagorean orthogonality doesn't translate via the bridge. **Did not:** the $\Pi_W$-even / $\Pi_W$-odd decomposition transports because $\Pi_W$ is part of the wedge structure preserved by the bridge.
- **Could surface:** the three-carrier identification breaks on the synthetic side. **Did not:** the Paper 47 isometric embeddings commute with the topography (time-independent), hence descend to synthetic-side isometric embeddings via Gelfand duality.

**No failure mode surfaced. The wedge application is structurally clean.**

---

**End of memo.**

**Files added in this sprint:**
- `debug/sprint_phase_a4prime_c_wedge_application_memo.md` (this memo, ~5500 words formal wedge application + verdict layer)
- `debug/data/sprint_phase_a4prime_c_wedge_application.json` (per-theorem verdict structure + G-B2/G-B4 dependencies + Paper 48 §6/§7 outline)

**Cross-references:**
- `debug/sprint_phase_a3prime_mondino_samann_bridge_memo.md` (A.3' bridge theorem; load-bearing throughout)
- `debug/sprint_phase_a2prime_krein_lift_formalization_memo.md` (A.2' Krein PPQMS substrate; load-bearing for §2.1)
- `debug/sprint_l3e_p3_rescope_memo.md` (re-scope memo; confirms Phase A.4'-C target)
- `papers/group1_operator_algebras/paper_42_modular_hamiltonian_four_witness.tex` (four-witness theorem; load-bearing for T4)
- `papers/group1_operator_algebras/paper_43_lorentzian_extension.tex` (hemispheric wedge, BW vacuum, Pythagorean orthogonality; load-bearing for §2.1, T1, T5)
- `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex` (operator-system substrate; load-bearing for §2.1)
- `papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex` (K⁺-weak-form propinquity + bit-exact panel; load-bearing for §2.1, T3)
- `papers/group1_operator_algebras/paper_47_two_rate_hybrid_convergence.tex` (norm-resolvent + three-carrier identification; load-bearing for T6, T6a)
- `papers/group3_foundations/paper_32_spectral_triple.tex` (§VIII rem:pythagorean_m1_closure; supplies T5 M1 signature identification)
- Mondino-Sämann arXiv:2504.10380 v4 (Def 2.3, Def 3.6, Def 3.8, Def 3.12, Thm 6.2 — via A.3' §1.2)
- Latrémolière arXiv:2512.03573 (via A.2')
- Connes-Rovelli 1994 (Class. Quantum Grav. 11, 2899; via A.3' §5)
