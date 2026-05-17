# Paper 43 (outline) — Lorentzian extension of the four-witness Wick-rotation theorem at finite cutoff

**Status:** OUTLINE ONLY (no .tex draft yet).
**Date:** 2026-05-17.
**Sprint origin:** Sprint L2 closure (sub-sprints L2-A through L2-G). See `debug/sprint_l2_synthesis_memo.md`.
**Companion paper(s):** Paper 42 (Riemannian closure of the four-witness theorem on SU(2) truncated triple). When drafted, **Paper 43 does NOT supersede Paper 42**. Paper 42 is the Riemannian closure; Paper 43 is the Lorentzian extension; together they constitute the four-witness Wick-rotation arc at the operator-system level.

**Target venue:** math.OA / J. Geom. Phys. — sibling to Papers 38, 39, 40, 42 in `papers/standalone/`. Natural Zenodo deposit when drafted (per CLAUDE.md §6 papers_zenodo_not_journals rule).

**Estimated draft length:** ~25-30 pages, ~12,000-15,000 words.

**Estimated draft scope:** 2 weeks of writing + 1 week of review.

---

## §1. Background and statement of the main theorem

### §1.1 Background

- Paper 38 (WH1 PROVEN): qualitative-rate Latrémolière-propinquity GH-convergence theorem on the Riemannian truncated triple $\mathcal{T}_{n_{\max}} \to \mathcal{T}_{S^3}$. The continuum spectral triple is the round-S³ Camporesi–Higuchi triple.
- Paper 39: tensor-product propinquity convergence ($\mathcal{T}_{S^3}^{\lambda_a} \otimes \mathcal{T}_{S^3}^{\lambda_b}$, same manifold, distinct focal lengths).
- Paper 40: unified rate constant $4/\pi$ across compact connected Lie groups with bi-invariant metric.
- Paper 42: Tomita-Takesaki modular structure on $\mathcal{T}_{n_{\max}}$; four-witness Wick-rotation theorem closes at the operator-system level (Riemannian, finite cutoff).
- Sprint TD Track 4 (2026-05-09) + Sprint Unruh-pendant (2026-05-10): four-witness Wick-rotation theorem at the metric-functional level ($\beta = 2\pi/\kappa_g$ for HH/Sewell/BW/Unruh).

### §1.2 The Lorentzian-extension question

Paper 42 §10 open question O1: extend the operator-system-level closure of the four-witness theorem to signature $(s, t) = (3, 1)$ Lorentzian. The natural construction recipe is van den Dungen 2016 Prop 4.1 (arXiv:1505.01939):

> Let $(M, g)$ be a pseudo-Riemannian spin manifold of signature $(s, t)$ with a spacelike reflection $r$. Then $(C^\infty_c(M), L^2(\mathbb{S}), i^t \, \not{\!\!D}_{g_r}, \mathcal{J}_M)$ is an even Krein spectral triple with $\mathcal{J}_M = \gamma(e_0)$.

For $(M, g) = (S^3 \times \mathbb{R}, ds^2_{S^3} - dt^2)$ at $(s, t) = (3, 1)$: $g_r = ds^2_{S^3} + dt^2$ (positive-definite Wick-rotated 4-manifold), $i^t = i$, $\mathcal{J} = \gamma^0$.

The signature classification of Bizi–Brouder–Besnard 2018 (arXiv:1611.07062) places Lorentzian $(s, t) = (3, 1)$ West-coast at $(m, n) = (4, 6)$ via BBB Table 3 ($m = t+s = 4$, $n = t-s = -2 \equiv 6 \pmod 8$).

### §1.3 Main theorem (Krein-level four-witness Wick-rotation theorem at finite cutoff)

**Theorem 1.1 (Krein-level four-witness Wick-rotation theorem).** On the Lorentzian Krein space $\mathcal{K}_{n_{\max}, N_t} := \mathcal{H}_{\mathrm{GV}}^{n_{\max}} \otimes \mathbb{C}^{N_t}$ at signature $(s, t) = (3, 1)$ West-coast in the Peskin-Schroeder chiral basis, with fundamental symmetry $J = \gamma^0$, Lorentzian Dirac $D_L = i \cdot [\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t}]$, hemispheric wedge $W_L = P_W^{\mathrm{spatial}} \otimes P_{t \geq 0}$, and wedge KMS state $\rho_W^L = e^{-K_L^{\alpha, W}}/Z$ under the BW choice $H_{\mathrm{local}} = K_L^{\alpha, W}/\beta$, the following hold bit-exactly at every $(n_{\max}, N_t) \in \{1, 2, 3\} \times \{1, 11, 21\}$:

1. **BW-α period closure (geometric):** $\sigma_{2\pi}^{L, \alpha}(O) = O$ for all $O \in B(\mathcal{K}_{W_L})$, residual $\leq 4 \times 10^{-16}$.
2. **BW-γ Tomita period closure:** $\sigma_{2\pi}^{L, \mathrm{TT}}(a) = a$ on Krein-GNS Hilbert-Schmidt space, same residual.
3. **Flow conjugacy at general $t$:** $\sigma_t^{L, \mathrm{TT}}(a) = \sigma_{-t}^{L, \alpha}(a)$, residual $\leq 4 \times 10^{-16}$.
4. **Six-witness collapse:** all six instantiations (BW, $\mathrm{HH}_{M=1}$, $\mathrm{HH}_{M=2}$, $\mathrm{Sew}_{M=1}$, $\mathrm{Unruh}_{a=1}$, $\mathrm{Unruh}_{a=2}$) give bit-identical $\Delta_L$ and $K_L^{\mathrm{TT}}$.

**This lifts the four-witness Wick-rotation theorem from structural correspondence at the metric-functional level to literal identification at the operator-system level (Lorentzian, finite cutoff).** Framework's first operator-system-level literal identification of a Lorentzian QFT theorem at finite cutoff.

---

## §2. Setup: Krein space + Lorentzian Dirac + Connes axioms at (m, n) = (4, 6)

### §2.1 The Camporesi-Higuchi spatial spinor bundle (recap from Paper 32 §III + Paper 42 §2)

- $\mathcal{H}_{\mathrm{GV}}^{n_{\max}}$ at truncation $n_{\max}$: $\dim = N_{\mathrm{Dirac}}(n_{\max}) = (2/3) n_{\max} (n_{\max} + 1)(n_{\max} + 2)$.
- Basis labels $(n_{\mathrm{fock}}, l, m_j, \chi)$ via `FullDiracLabel`, with $\chi \in \{+1, -1\}$ the chirality grading.
- Camporesi-Higuchi spectrum $|\lambda_n| = n + 1/2$ (Fock convention) or $n + 3/2$ (CH convention), with multiplicities $g_n^{\mathrm{Dirac}} = 2 n (n+1)$.

### §2.2 The Cl(3, 1) chiral-basis gamma matrices

- Standard Peskin-Schroeder chiral basis (West-coast metric $\eta = \mathrm{diag}(+, -, -, -)$):

$$\gamma^0 = \begin{pmatrix} 0 & I_2 \\ I_2 & 0 \end{pmatrix}, \quad \gamma^i = \begin{pmatrix} 0 & \sigma^i \\ -\sigma^i & 0 \end{pmatrix}, \quad \gamma^5 = i \gamma^0 \gamma^1 \gamma^2 \gamma^3 = \mathrm{diag}(-I_2, +I_2).$$

- Clifford algebra $\{\gamma^\mu, \gamma^\nu\} = 2 \eta^{\mu\nu} I_4$ verified bit-exact.
- $(\gamma^5)^2 = +I_4$; $\{\gamma^5, \gamma^\mu\} = 0$ for all $\mu$.

### §2.3 The Krein space (Sprint L2-B)

- $\mathcal{K}_{n_{\max}, N_t} = \mathcal{H}_{\mathrm{GV}}^{n_{\max}} \otimes \mathbb{C}^{N_t}$ with bounded uniform temporal grid $t_k = -T_{\max} + k(2T_{\max})/(N_t - 1)$.
- Fundamental symmetry $J = J_{\mathrm{spatial}} \otimes I_{N_t}$ with $J_{\mathrm{spatial}}$ the lift of $\gamma^0$ (the chirality-swap in the chiral basis) to chirality-doubled $\mathcal{H}_{\mathrm{GV}}$.
- Krein inner product $\langle \psi, \phi \rangle_{\mathcal{K}} = \langle \psi, J \phi \rangle$.
- Bit-exact axioms: $J^2 = +I$, $J^* J = I$, $J = J^*$, $\mathcal{K} = \mathcal{K}^+ \oplus \mathcal{K}^-$ with $\dim \mathcal{K}^\pm = \dim \mathcal{K} / 2$.

### §2.4 The Lorentzian Dirac (Sprint L2-C)

- $D_L = i \cdot [\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t}]$ per vdD Prop 4.1.
- $\partial_t$ centered FD + Dirichlet zero BC: anti-Hermitian on the finite grid. NOT periodic (would create CTCs on $S^3 \times S^1$ per Geroch's theorem).
- Sign $i^t = +i$ derived from Krein-self-adjointness requirement $D_L^\times = D_L$.
- Krein-self-adjointness bit-exact at every $(n_{\max}, N_t)$.
- Riemannian-limit recovery (load-bearing) bit-identical at $N_t = 1$: $D_L|_{N_t = 1} = i \cdot D_{\mathrm{GV}}$.

### §2.5 Connes axiom audit at (m, n) = (4, 6) (Sprint L2-D)

- BBB Table 1 signs at $(m, n) = (4, 6)$ verified from arXiv:1611.07062 v2 directly: $\varepsilon = +1$, $\varepsilon'' = -1$, $\kappa = -1$, $\kappa'' = +1$.
- $J_L$ on Krein space via lift of 4-spinor charge conjugation $U_4 = i \gamma^2 = (i \sigma_y)_{\mathrm{chir}} \otimes (i \sigma^2)_{\mathrm{spin}}$.
- Four BBB-predicted-sign axioms hold bit-exactly: $J_L^2 = +I$, $\{J_L, \gamma^5\} = 0$, $\{J_L, \eta\} = 0$, $J_L D_L = +D_L J_L$.
- Structural finding: BBB universal axiom $\chi D = -D\chi$ fails on truthful Camporesi-Higuchi $D_{\mathrm{GV}}$. Construction uses R1 resolution (truthful $D_{\mathrm{GV}}$, preserves Riemannian-limit bit-exactness, accept structural finding).

---

## §3. BW-α construction (geometric)

### §3.1 The Lorentzian hemispheric wedge

$$W_L := P_W^{\mathrm{spatial}} \otimes P_{t \geq 0}$$

- $P_W^{\mathrm{spatial}} = (1/2)(I + R_{\mathrm{polar}})$ — Paper 42 Definition 4.1 hemispheric wedge ($m_j$-reflection involution on the spinor basis).
- $P_{t \geq 0}$ — diagonal projector on $\mathbb{C}^{N_t}$ selecting $t_k \geq 0$.
- Properties: $P_{W_L}^2 = P_{W_L}$, $P_{W_L}^* = P_{W_L}$.
- $\dim W_L = (\dim \mathcal{H}_{\mathrm{GV}}/2) \times N_{t, +}$.
- Riemannian-limit: $P_{W_L}|_{N_t = 1} = P_W^{\mathrm{spatial}}$ bit-identically.

### §3.2 The BW choice of local Hamiltonian

- $H_{\mathrm{local}} := K_L^{\alpha, W} / \beta$ at $\beta = 2\pi$.
- Wedge KMS state $\rho_W^L = e^{-\beta H_{\mathrm{local}}}/Z = e^{-K_L^{\alpha, W}}/Z$.
- $\beta$-independent at the algebra-action level — same property as Paper 42 §4.2.

### §3.3 The BW-α geometric generator

- $K_L^\alpha := K_\alpha^{\mathrm{spatial}} \otimes I_{N_t}$ with $K_\alpha^{\mathrm{spatial}} = \mathrm{diag}(\mathrm{two}\_m_j)$ inheriting integer eigenvalues from Paper 42 Definition 5.1.
- Wedge-restricted $K_L^{\alpha, W}$ via the "unfolded" basis (Paper 42 §5.2 lifted).
- Spectrum $\subset \mathbb{Z}$ (odd integers) — preserved at $(3, 1)$.

### §3.4 The modular flow

- $\sigma_t^{L, \alpha}(O) := e^{i t K_L^{\alpha, W}} \, O \, e^{-i t K_L^{\alpha, W}}$.
- **Theorem 3.1 (BW-α period closure, bit-exact at finite cutoff):** $\sigma_{2\pi}^{L, \alpha}(O) = O$ identically (machine precision residual due to round-off accumulation only).
- Proof: $K_L^{\alpha, W}$ has integer eigenvalues, so $e^{i \cdot 2\pi \cdot K_L^{\alpha, W}} = I$ identically.

---

## §4. BW-γ construction (Tomita-Takesaki on Krein-GNS)

### §4.1 The Krein-GNS Hilbert-Schmidt space

- $\mathcal{K}_{\mathrm{GNS}} := M_{\dim W_L}(\mathbb{C}) \cong \mathbb{C}^{\dim W_L^2}$, inner product $\langle X, Y \rangle_{\mathrm{GNS}} = \mathrm{Tr}(X^* Y)$.
- Cyclic vector $\Omega := (\rho_W^L)^{1/2}$.
- Convention (Krein-positive completion, per van den Dungen 2016 §2): tracial-Gibbs structure on a Type I factor.

### §4.2 The modular $S$ operator and polar decomposition

- Tomita $S$ operator $S: a\Omega \mapsto a^* \Omega$.
- Polar decomposition $S = J_L^{\mathrm{TT}} \cdot \Delta_L^{1/2}$ with $\Delta_L = ((\rho_W^L)^{-1})^T \otimes \rho_W^L$ (column-stacked vec representation).
- $J_L^{\mathrm{TT}}(X) = (\rho_W^L)^{1/2} \, X^* \, (\rho_W^L)^{-1/2}$ acting antilinearly.

### §4.3 The modular Hamiltonian

- $K_L^{\mathrm{TT}} := -\log \Delta_L$.
- **Proposition 4.2 (Integer spectrum of $K_L^{\mathrm{TT}}$):** $\mathrm{Spec}(K_L^{\mathrm{TT}}) \subset \mathbb{Z}$ (integer differences of $K_L^{\alpha, W}$ eigenvalues).
- **Proposition 4.3 ($J_L^{\mathrm{TT}}{}^2 = +I$, categorical distinction from Connes $J_{\mathrm{GV}}$ with $J_{\mathrm{GV}}^2 = -I$).**

### §4.4 The modular flow on the algebra

- $\sigma_t^{L, \mathrm{TT}}(a) := \Delta_L^{it}(a) = e^{-i t K_L^{\alpha, W}} \, a \, e^{+i t K_L^{\alpha, W}}$ (using $\rho_W^L = e^{-K_L^{\alpha, W}}/Z$).

### §4.5 Period closure (bit-exact)

- **Theorem 4.4 (BW-γ period closure at finite Krein cutoff):** $\sigma_{2\pi}^{L, \mathrm{TT}}(a) = a$ bit-exact for all $a \in B(\mathcal{K}_{W_L})$.

---

## §5. Flow conjugacy

### §5.1 The conjugacy theorem

- **Theorem 5.1 (Flow conjugacy at general $t$):** $\sigma_t^{L, \mathrm{TT}}(a) = \sigma_{-t}^{L, \alpha}(a)$ for every $t \in \mathbb{R}$ and every $a \in B(\mathcal{K}_{W_L})$.
- Proof: identical operator actions $e^{-i t K_L^{\alpha, W}} \cdot e^{+i t K_L^{\alpha, W}}$ in both flows.

### §5.2 Bit-exact numerical verification

- Cross-flow difference at $t \in \{1, \pi, 2\pi\}$ on five test multipliers: residual $\leq 4 \times 10^{-16}$ at every tested cell.

---

## §6. Six-witness collapse

### §6.1 The Krein-level collapse

- **Corollary 6.1 (Single operator-system construction for the six witnesses at $(3, 1)$):** Under the unit normalisation $\kappa_g = 1$ and the wedge KMS state choice $\rho_W^L = e^{-K_L^{\alpha, W}}/Z$, the six witness instantiations (BW, $\mathrm{HH}_{M=1}$, $\mathrm{HH}_{M=2}$, $\mathrm{Sew}_{M=1}$, $\mathrm{Unruh}_{a=1}$, $\mathrm{Unruh}_{a=2}$) collapse to a single construction. $\Delta_L$ and $K_L^{\mathrm{TT}}$ are bit-identical across the six instantiations.
- Proof: $\rho_W^L$ is $\beta$-independent under the BW choice; the modular construction depends only on $\rho_W^L$.

### §6.2 Bit-exact cross-witness consistency

- Cross-witness residual difference is exactly $0.0$ at every tested $(n_{\max}, N_t)$.

### §6.3 Structural reading

- The witness-specific physical content ($M$ for HH/Sew, $a$ for Unruh, $\kappa_g = 1$ for BW canonical) parameterises the correspondence to continuum physical observables but does NOT modify the underlying modular-Hamiltonian construction.
- The $2\pi$ in the modular period is the M1 sub-mechanism of the master Mellin engine (Hopf-base measure $\mathrm{Vol}(S^1)$, Paper 32 §VIII case-exhaustion theorem).
- **Lifts the four-witness Wick-rotation theorem from "structural correspondence at the metric-functional level" (Sprint TD Track 4, Unruh-pendant) to "literal identification at the operator-system level (Lorentzian, finite cutoff)"** at every finite Krein cutoff.

---

## §7. Structural results

### §7.1 The $H_{\mathrm{local}}$ signature-independence finding (the headline)

- **Theorem 7.1 ($H_{\mathrm{local}}$ signature-independence at Riemannian limit):** $\| H_{\mathrm{local}}|_{(3, 1), N_t = 1} - D_L^W|_{N_t = 1} \|_F = \| H_{\mathrm{local}}|_{(3, 0)} - D_{\mathrm{GV}}^W \|_F$ bit-exact at every tested $n_{\max} \in \{1, 2, 3\}$. Values: $2.1332$, $6.5275$, $13.854$.
- At $N_t > 1$ the Lorentzian residual is refined upward by the temporal-derivative content of $D_L$.

### §7.2 Structural reading

- Paper 42 §7.2 / §10 O3 distinction (spectral-action Dirac vs. modular-Hamiltonian generator on BW vacuum) is **signature-INDEPENDENT at the Riemannian limit**.
- NOT a peculiarity of the round-$S^3$ Riemannian truncation.
- A **deeper structural feature** of the framework's modular content.
- Refined upward at $N_t > 1$ by temporal-derivative content of $D_L = i(\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t})$.

### §7.3 The BBB universal-axiom structural finding (from Sprint L2-D)

- BBB universal $\chi D = -D \chi$ (Sec 5(v)) fails on truthful Camporesi-Higuchi $D_{\mathrm{GV}}$.
- Structural reason: GeoVac's chirality-as-$\gamma^5$ identification + chirality-diagonal $D_{\mathrm{GV}}$ + BBB universal axiom is a mutually-inconsistent triple.
- Not a basis-convention bug; pick any two of the three, the third must go.
- Sprint L2-E uses R1 resolution (truthful $D_{\mathrm{GV}}$, preserves Riemannian-limit bit-exactness); the wedge construction is $K_L^{\alpha, W}$-driven, not $D_L$-driven, so the period closure is $D_L$-independent at the operator-action level.

---

## §8. Honest scope

### §8.1 Finite cutoff only

- Theorem 1.1 closure is at finite $n_{\max}$ and finite $N_t$.
- NOT in the continuum / GH-limit / Lorentzian-propinquity limit.
- No published Lorentzian propinquity construction exists as of May 2026.

### §8.2 Lorentzian-propinquity extension (Sprint L3) is open

- Latrémolière 2017/2026, Hekkelman-McDonald 2024 a/b, Toyota 2023, Farsi-Latrémolière 2024/2025 all strictly Riemannian.
- Constructing Lorentzian propinquity is original NCG-mathematics work at 6-12 month scale.
- Nieuviarts 2025 twist-morphism (arXiv:2502.18105) is a candidate shortcut, pending L2-Nieuviarts-scoping verification of applicability to $S^3 = \mathrm{SU}(2)$ (odd-dim caveat).

### §8.3 Cross-manifold W2b not addressed

- $\mathcal{T}_{S^3} \otimes \mathcal{T}_{\mathrm{Hardy}(S^5)}$ remains structurally open.
- Distinct frontier from the (3, 1) extension within $\mathcal{T}_{S^3}$.
- Blocked at NCG-framework level by Riemannian-vs-Hardy-sector asymmetry (Paper 24 §V Coulomb/HO category mismatch).

### §8.4 Calibration data (W3, second packing axiom) not addressed

- Sprint L2 is structurally additive at the spectral-triple-machinery level.
- Inner-factor-input-data question (Yukawa Dirichlet ring, calibration constants) untouched.

---

## §9. Open questions

1. **Lorentzian propinquity construction (Sprint L3).** Multi-month NCG-mathematics work. May be shortened via Nieuviarts 2025 morphism if applicability to $S^3 = \mathrm{SU}(2)$ verifies.
2. **Cross-manifold W2b.** Structural extension to $\mathcal{T}_{S^3} \otimes \mathcal{T}_{\mathrm{Hardy}(S^5)}$.
3. **Deeper $H_{\mathrm{local}}$ question.** Why is the spectral-action Dirac not the modular generator on the BW vacuum, structurally? (Residual O3 from Paper 42 §10, now refined to "signature-independent" but not closed at the structural-explanation level.)
4. **Non-tracial wedge states.** Paper 42 O2 lifted to (3, 1).
5. **Higher-rank compact non-abelian extensions.** Paper 42 O4 lifted to (3, 1).
6. **BBB universal-axiom-compatible construction.** Sprint L2-D R3 resolution (redefine $\gamma^5$ as off-diagonal grading) remains a possible direction for a "full BBB Krein spectral triple in the strict Sec 5(v) sense" with truthful $D_{\mathrm{GV}}$.

---

## §10. References

Core citations expected in Paper 43:

- **van den Dungen 2016** (Math. Phys. Anal. Geom. 19:4, arXiv:1505.01939) — Prop 4.1 construction recipe.
- **Bizi–Brouder–Besnard 2018** (J. Math. Phys. 59:062303, arXiv:1611.07062) — $(m, n)$ classification, Table 1 signs at (4, 6).
- **Hartle–Hawking 1976** — KMS state on Euclidean Schwarzschild cigar at $\beta = 8\pi M$.
- **Bisognano–Wichmann 1976** — modular automorphism / Lorentz-boost identification on Rindler wedge.
- **Sewell 1982** — KMS state on bifurcate Killing horizon.
- **Unruh 1976** — accelerated-observer thermal bath at $T = a/(2\pi)$.
- **Paper 42** — Riemannian closure of the four-witness theorem on $\mathcal{T}_{n_{\max}}$.
- **Paper 38** — WH1 PROVEN, qualitative-rate propinquity on the Riemannian truncated triple.
- **Paper 39** — tensor-product propinquity ($\mathcal{T}_{S^3}^{\lambda_a} \otimes \mathcal{T}_{S^3}^{\lambda_b}$).
- **Paper 40** — unified compact-Lie-group propinquity convergence with $4/\pi$ rate.
- **Paper 32** — GeoVac spectral triple synthesis, §VIII Mellin engine + §VIII.E Lorentzian transfer audit.
- **Paper 34** — Projection taxonomy, §III.27 Wick rotation entry.
- **Paper 31** — Universal/Coulomb partition, §8 Signature Partition.
- **Connes-Rovelli 1994** — Thermal time hypothesis (background for the H_local discussion).
- **Latrémolière 2017/2026** — Riemannian propinquity (referenced for the L3 distinction).
- **Nieuviarts 2025** (arXiv:2502.18105) — Twist-morphism approach to pseudo-Riemannian spectral triples (referenced for L3 shortcut).

---

## §11. Draft sequencing

Estimated 2-3 week draft schedule:

- **Week 1:** §§1-2 (background, statement, setup) + §6 (six-witness collapse, the structural summary).
- **Week 2:** §§3-5 (BW-α + BW-γ + flow conjugacy, the main construction).
- **Week 3:** §§7-9 (structural results, honest scope, open questions) + §10 (references) + first internal review pass + LaTeX cleanup.

The numerical-table content is already in `debug/data/l2_e_modular_hamiltonian_lorentzian.json` and can be transcribed directly. The proof sketches in Paper 42 §§5-8 transport verbatim with the temporal slot added; Paper 43 should mirror Paper 42's structure section by section.

---

## §12. Relationship to Paper 42 (CRITICAL — NOT a supersession)

**When Paper 43 is drafted, it does NOT supersede Paper 42.** The relationship is:

- **Paper 42:** Riemannian closure of the four-witness Wick-rotation theorem on the SU(2) truncated triple $\mathcal{T}_{n_{\max}}$. Header section: "Tomita–Takesaki modular structure on truncated SU(2) spectral triples: a four-witness Wick-rotation theorem at finite cutoff."

- **Paper 43:** Lorentzian extension of the four-witness theorem to signature $(3, 1)$ on the Krein space $\mathcal{K}_{n_{\max}, N_t}$. Header section: "Krein-level Tomita–Takesaki modular structure at signature (3, 1) and the operator-system-level four-witness Wick-rotation theorem at finite cutoff."

Together, Papers 42 + 43 constitute the four-witness Wick-rotation arc at the operator-system level — Paper 42 the Riemannian half, Paper 43 the Lorentzian half. The Riemannian closure stands on its own (it is the WH1-PROVEN result for the modular-Hamiltonian on $\mathcal{T}_{n_{\max}}$); the Lorentzian extension extends it without altering the Riemannian foundation.

Paper 43's §3-§6 should reference Paper 42 §§5-8 explicitly for proof transport ("the construction mirrors Paper 42 §X verbatim with the temporal slot added"), keeping Paper 43 readable as a Lorentzian-extension paper without requiring full re-derivation of the Riemannian-side machinery.

The math.OA quartet (Papers 38, 39, 40, 42, 43 — actually a quintet after Paper 43 is drafted) collectively documents the framework's NCG-mathematical content: Paper 38 propinquity convergence on $\mathcal{T}_{S^3}$ (the WH1 PROVEN result); Paper 39 tensor-product propinquity ($\mathcal{T}_{S^3}^{\lambda_a} \otimes \mathcal{T}_{S^3}^{\lambda_b}$); Paper 40 unified compact-Lie-group propinquity with $4/\pi$ rate; Paper 42 Riemannian Tomita-Takesaki modular structure + four-witness Wick-rotation theorem; Paper 43 Lorentzian extension of the four-witness theorem at finite cutoff.

---

End of outline.
