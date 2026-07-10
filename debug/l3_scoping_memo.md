# Sprint L3 Scoping Memo — Continuum Lorentzian Propinquity on the GeoVac Camporesi–Higuchi Spectral Triple

**Date:** 2026-05-17
**Author:** L3-scoping Track 1 PM (Claude)
**Status:** Forward-looking architecture document. NOT a proof attempt; NOT a theorem; NOT a commitment to open Sprint L3.
**Inputs:** CLAUDE.md §1.7 (WH1 PROVEN + Sprint L2 closure paragraphs); Paper 38 §3–§7 (five-lemma proof); Paper 42 §7.2, §11; Paper 43 outline; `debug/sprint_l2_synthesis_memo.md`; `debug/lorentzian_l0_audit_memo.md`; `debug/lorentzian_partition_transfer_memo.md`.
**Sister deliverable:** Track 2 literature parallel (not produced by this PM).

---

## §1. L3 deliverable definition — what does "continuum Lorentzian propinquity convergence" actually mean?

### §1.1 The Riemannian template (what we are trying to extend)

Paper 38 proves: the truncated Camporesi–Higuchi spectral triples $\mathcal{T}_{n_{\max}}$ converge to the round-$S^3$ Camporesi–Higuchi spectral triple $\mathcal{T}_{S^3}$ in the **Latrémolière quantum Gromov–Hausdorff propinquity for metric spectral triples** of Latrémolière 2017 (arXiv:1703.07073):

$$\Lambda(\mathcal{T}_{n_{\max}}, \mathcal{T}_{S^3}) \;\le\; C_3 \cdot \gamma_{n_{\max}} \;\to\; 0 \quad \text{as } n_{\max}\to\infty,$$

with $C_3 = 1$ from L3 and $\gamma_{n_{\max}} = (4/\pi)\log n_{\max}/n_{\max} + O(1/n_{\max})$ from L2. The propinquity is defined by

$$\Lambda(\mathcal{T}_1, \mathcal{T}_2) := \inf_\tau \mathrm{length}(\tau), \qquad \mathrm{length}(\tau) = \max(\mathrm{reach}(\tau), \mathrm{height}(\tau)),$$

with the infimum over **tunnels** $\tau$ connecting $\mathcal{T}_1$ and $\mathcal{T}_2$. A tunnel between two metric spectral triples is a pair of UCP (unital completely positive) maps relating their underlying operator systems, satisfying compatibility with their Dirac-operator-induced Lipschitz seminorms. The reach measures the approximate-identity error of the UCP maps; the height measures their Lipschitz-distortion envelope.

The propinquity completion's natural limit object is, per Paper 38 Proposition (limit identification), the **Wasserstein–Kantorovich state space** on $S^3$: $(\mathcal{P}(S^3), d_{\mathrm{Wass}})$. Kantorovich–Rubinstein duality identifies the Connes distance on $\mathcal{T}_{S^3}$ with $d_{\mathrm{Wass}}$ on $\mathcal{P}(S^3)$ in the dual sense.

### §1.2 The Lorentzian deliverable — first-pass formalisation

The L3 deliverable is a Lorentzian-signature analog of the above statement. The candidate Krein-space continuum object is the van den Dungen Krein spectral triple $(C_c^\infty(S^3 \times \mathbb{R}), L^2(\mathbb{S}), i \not{\!\!D}_{g_r}, J = \gamma^0)$ at $(s, t) = (3, 1)$, West-coast metric. Call this object $\mathcal{T}_{S^3 \times \mathbb{R}, \mathrm{Krein}}$.

The candidate truncated object is the Sprint L2 finite-cutoff Krein space $\mathcal{K}_{n_{\max}, N_t} = \mathcal{H}_{\mathrm{GV}}^{n_{\max}} \otimes \mathbb{C}^{N_t}$ equipped with the Lorentzian Dirac $D_L = i[\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t}]$ and fundamental symmetry $J = \gamma^0 \otimes I_{N_t}$, in the chiral basis at the BBB $(m, n) = (4, 6)$ classification. Call this object $\mathcal{T}_{n_{\max}, N_t}^{\mathrm{Krein}}$.

**L3 Claim (placeholder):** There exists a notion of **Lorentzian propinquity** $\Lambda^L(\cdot, \cdot)$ such that $\Lambda^L(\mathcal{T}_{n_{\max}, N_t}^{\mathrm{Krein}}, \mathcal{T}_{S^3 \times \mathbb{R}, \mathrm{Krein}}) \to 0$ as $(n_{\max}, N_t) \to (\infty, \infty)$ in a specified joint limit.

This is the **strong form** of L3 and is most likely **unreachable in the published literature's current state**.

### §1.3 What changes structurally in the Lorentzian setting

Three structural shifts make the Riemannian template not transferable verbatim:

**(a) The Hilbert space is replaced by a Krein space.** The natural inner product $\langle \psi, \phi \rangle_{\mathcal{K}} = \langle \psi, J\phi \rangle$ is **not positive-definite**. The Lorentzian Dirac is **Krein-self-adjoint** ($D_L^\times = J D_L^\dagger J = D_L$), not Hilbert-space-self-adjoint. The Connes distance formula

$$d(\varphi, \psi) = \sup\{|\varphi(a) - \psi(a)| : a \in \mathcal{A}, \|[D, a]\|_{\mathrm{op}} \le 1\}$$

uses an operator norm; on a Krein space, the natural norm replacement is either (i) the **Hilbert-space operator norm** of the Krein-positive completion (van den Dungen 2016 §2), or (ii) the **graph norm** on $D_L$-bounded operators. Either choice yields a non-trivial Lipschitz seminorm, but the relationship between (i), (ii), and the Krein-positive sub-cone is not standardly worked out.

**(b) The Wasserstein–Kantorovich state space requires a Lorentzian-distance ground.** Wasserstein–Kantorovich on $S^3$ uses the round-$S^3$ geodesic distance — a positive-definite metric. On Lorentzian $S^3 \times \mathbb{R}$, the natural causal-set structure replaces "distance" with **Lorentzian path-length** (a signed quantity that is imaginary along spacelike intervals and zero on null lines). There is no published Wasserstein-on-Lorentzian-spacetime construction that produces a propinquity-compatible state-space metric. (Eckstein–Miller 2017 "Causality for nonlocal phenomena" and Mainiero 2024 "Lorentzian optimal transport" are starting points but do not give a propinquity-completion.)

**(c) The propinquity limit object is no longer (clearly) "states on the algebra."** In the Riemannian case, $(\mathcal{P}(S^3), d_{\mathrm{Wass}})$ is the natural propinquity limit by Kantorovich–Rubinstein duality. In the Lorentzian case, the analog of "states on the algebra" must respect the Krein-product structure; a candidate is **Krein-positive states** (functionals $\varphi$ with $\varphi(a^* J a) \ge 0$), but the resulting state space is not a simplex in the usual sense and lacks a direct Wasserstein analog.

### §1.4 A weaker form that may be tractable

**L3-Krein-positive-restricted Claim:** Restrict the propinquity to the **Krein-positive cone** $\mathcal{K}^+_{n_{\max}, N_t}$ (the subspace where $\langle\psi, J\psi\rangle > 0$) on both sides. On $\mathcal{K}^+$, the Krein inner product reduces to a genuine Hilbert-space inner product. The Krein-positive completion of van den Dungen 2016 §2 gives a Type-I-factor Hilbert-Schmidt structure on this restricted cone. The L3-Krein-positive-restricted claim is then that the **Riemannian-propinquity machinery transfers to the Krein-positive cone with the Krein-positive Hilbert space replacing $\mathcal{H}_{\mathrm{GV}}$**. This is structurally the deliverable Sprint L2-E used at finite cutoff (the wedge KMS state $\rho_W^L$ is supported on the Krein-positive cone by construction); the L3 question is whether the limit $(n_{\max}, N_t) \to (\infty, \infty)$ exists in the same propinquity sense as Riemannian.

This is the **likely realistic shape** of the eventual L3 deliverable. It is partial convergence (positive-cone only) rather than full Krein-space convergence; but it lifts Sprint L2-E's finite-cutoff closure to the continuum, which is the load-bearing physics question (the Wick-rotation theorem at the operator-system level in the continuum limit).

---

## §2. Five-lemma transfer table (L1' / L2 / L3 / L4 / L5)

For each Riemannian lemma in Paper 38, the assessment is for the **L3-Krein-positive-restricted claim** (§1.4), since the strong claim (§1.2) is likely unreachable without 6–12 months of new Krein-Wasserstein optimal-transport mathematics.

### L1' — Chirality-doubled operator system

**Riemannian statement.** The truncated operator system $\mathcal{O}_{n_{\max}}$ inherits a chirality block-diagonal structure from $\mathcal{H}_{\mathrm{GV}} = \mathcal{H}^+ \oplus \mathcal{H}^-$; the offdiag CH extension is the SDP-bounding device for Connes-distance finiteness on cross-shell pure-state pairs.

**Riemannian proof technique.** Direct construction of $\mathcal{O}_{n_{\max}} = P_{n_{\max}} C^\infty(S^3) P_{n_{\max}}$ with $P$ a spectral projection; verification of multiplier block-diagonality; SDP-finiteness via numerical computation at $n_{\max} \in \{2, 3\}$.

**Lorentzian transfer assessment: TRANSFERS_FREELY (Krein-positive cone) / NEEDS_NEW_MATH (full Krein).**
- Krein-positive restriction: the chirality-grading $\gamma^5$ is signature-blind in the chiral basis. The chirality block decomposition $\mathcal{K}^+ = \mathcal{K}^{+,+} \oplus \mathcal{K}^{+,-}$ inherits from the spatial $\gamma^5$-grading; the temporal slot $\mathbb{C}^{N_t}$ adds no chirality structure.
- Full Krein: the SDP-bounding role of offdiag CH needs a Krein-self-adjoint replacement; the BBB universal axiom $\chi D = -D\chi$ structurally fails on truthful $D_{\mathrm{GV}}$ (L2-D finding), so the offdiag construction must be re-derived to be Krein-self-adjoint.

**Closest published analog.** Bizi–Brouder–Besnard 2018 §3 (operator system on indefinite spectral triples) gives the framework for L1' lift; no published treatment of cross-shell Connes-distance SDP in the Krein setting.

**Risk if NEEDS_NEW_MATH.** Low for Krein-positive restriction (mostly bookkeeping); medium for full Krein (requires reformulating the offdiag SDP in Krein-self-adjoint form).

**Confidence: HIGH.** L1' transfer is a re-statement, not a new theorem.

### L2 — Central spectral Fejér kernel with $4/\pi$ asymptote

**Riemannian statement.** The central spectral Fejér kernel $K_{n_{\max}} = (1/Z_{n_{\max}})|\sum_{j \le j_{\max}} \sqrt{2j+1}\chi_j(g)|^2$ on $\mathrm{SU}(2)$ has cb-norm $2/(n_{\max}+1)$, Plancherel symbol $(2j+1)/Z_{n_{\max}}$, and mass-concentration moment $\gamma_{n_{\max}} \to 0$ with asymptote $(4/\pi)\log n_{\max}/n_{\max}$.

**Riemannian proof technique.** Stein–Weiss IBP on the SU(2) Haar measure; Bożejko–Fendler cb-norm equality on central Schur multipliers on amenable compact groups; Peter–Weyl decomposition.

**Lorentzian transfer assessment: STRUCTURALLY_BLOCKED on full Krein / NEEDS_NEW_MATH on Krein-positive restriction.**

This is the **most consequential blocker** in the L3 transfer table. The L2 kernel lives on $\mathrm{SU}(2)$, a **compact** group. The Lorentzian continuum is $S^3 \times \mathbb{R}$, with $\mathbb{R}$ **non-compact**. The Stein–Weiss IBP on a non-compact Lie group does not produce a Cesàro-2-style $\log n/n$ rate; instead, on $\mathbb{R}$ the natural Fejér kernel is the **Fejér-on-the-line kernel** $(\sin(Tt/2))^2/(\pi t^2 T)$ with a continuous spectral parameter $T$, and the rate is $O(1/T)$ rather than $O(\log n/n)$.

Combining a compact-direction kernel (the SU(2) part) with a non-compact-direction kernel (the $\mathbb{R}_t$ part) requires constructing a **joint Schur multiplier on the SU(2) × $\mathbb{R}$ Peter–Weyl × Fourier structure**, with two different rate regimes. The product structure works at the level of cb-norm (multiplicative), but the joint $\gamma$ moment is a **product of two convergence rates** $\gamma_{n_{\max}} \cdot \gamma_{T_{\max}}$ that the L3 limit must take simultaneously. This is original work; Leimbach–van Suijlekom 2024 covers $T^d$ (purely compact); Hekkelman–McDonald 2024 covers $T^d$ + spinor sector but still purely compact.

**Closest published analog.** Leimbach–van Suijlekom 2024 (Adv. Math.) for the compact factor; no published analog for the non-compact factor in propinquity context. (Fejér-on-the-line is standard harmonic analysis, but its lift to a Lipschitz-compatible UCP map between propinquity-related operator systems is unpublished.)

**Risk.** This is the **dominant cost of L3**. Estimated 2–4 months of focused work to lift the joint kernel construction. The Krein-positive restriction does not simplify this much: the joint compact × non-compact kernel is needed in either case, because Krein-positivity is a spatial-chirality condition that does not project away the $\mathbb{R}_t$ direction.

**Mitigation.** **Temporal compactification.** If the temporal direction is replaced by a torus $S^1_T$ with periodic boundary conditions ($T \to \infty$), then both factors are compact and the joint Peter–Weyl × Fourier kernel is a straightforward tensor product. **This is the realistic L3a entry point** (see §6).

**Confidence: MEDIUM-HIGH** on the assessment; **LOW** on which specific path through the obstruction is best.

### L3 — Lipschitz comparison with $C_3 = 1$

**Riemannian statement.** For every $f \in C^\infty(S^3)$ and every $n_{\max} \ge 1$, $\|[D_{\mathrm{CH}}, M_f]\|_{\mathrm{op}} \le C_3 \cdot \|\nabla f\|_{L^\infty(S^3)}$ with $C_3 = 1$ (asymptotically sharp, $< 1$ at every finite cutoff via the $\sqrt{(N-1)/(N+1)}$ ratio).

**Riemannian proof technique.** Direct computation in the Camporesi–Higuchi basis using $D_{\mathrm{CH}}$ diagonality and the SO(4) selection rule $|n - n'| \le N - 1$; per-harmonic Lipschitz scaling $\|\nabla Y^{(3)}_{NLM}\|_{L^\infty} \sim \sqrt{N^2-1}$; combination via linear triangle inequality.

**Lorentzian transfer assessment: NEEDS_NEW_MATH (different obstruction class than L2).**
- Truthful $D_L = i[\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t}]$ is **no longer diagonal** in the natural product basis (it has off-diagonal entries through $\gamma^0$ acting on chirality). The per-harmonic shell-difference factorisation that gave Paper 38 the clean $(N-1)/\sqrt{N^2-1}$ ratio does not transfer.
- A Lorentzian gradient norm $\|\nabla f\|_{L^\infty}$ on $S^3 \times \mathbb{R}$ has **mixed signature**. The natural replacement is the **Euclidean-Wick-rotated gradient norm** $\|\nabla_{g_r} f\|_{L^\infty}$ (positive-definite, using the Wick-rotated $g_r$), which is signature-blind. Under this convention, the L3 transfer is closer to verbatim.
- The key technical question is whether **the per-harmonic Lipschitz bound transfers to the Lorentzian Dirac under the Wick-rotated gradient norm**. The answer is likely yes for the spatial part (Paper 38 verbatim) and yes for the temporal part (1D Fejér kernel on $\mathbb{R}_t$ or $S^1_T$), with a cross-term that must be bounded separately.

**Closest published analog.** Hekkelman–McDonald 2024 (arXiv:2412.00628) lifts L3-style Lipschitz bounds to spinor bundles on $T^d$ via the Lichnerowicz formula. The mechanism transports to the Lorentzian case via Wick rotation if the gradient norm is interpreted as the $g_r$-gradient.

**Risk.** Low if the Wick-rotated gradient norm convention is accepted; medium if the genuine Lorentzian gradient (with signature) is demanded. The latter is structurally distinct and would require **Lichnerowicz-on-Lorentzian-manifolds** machinery — not standard.

**Confidence: HIGH** for Wick-rotated convention; **LOW** for genuine Lorentzian.

### L4 — Berezin reconstruction

**Riemannian statement.** The Berezin map $B_{n_{\max}}: C^\infty(S^3) \to \mathcal{O}_{n_{\max}}$ defined by Peter–Weyl coefficients with Plancherel weight $\hat{K}_{n_{\max}}(N) = N/Z_{n_{\max}}$ is (a) positivity-preserving on $f \ge 0$, (b) contractive, (c) an approximate identity at rate $\gamma_{n_{\max}}$, and (d) L3-compatible (Lipschitz norm preserved under $B_{n_{\max}}$).

**Riemannian proof technique.** Direct verification via convolution form $B_{n_{\max}}(f) = P_{n_{\max}}(K_{n_{\max}} * f) P_{n_{\max}}$; positivity from kernel positivity + projection; contractivity from Young's inequality; approximate identity from $\gamma$ rate; L3 compatibility from $\hat{K} \le 1$.

**Lorentzian transfer assessment: NEEDS_NEW_MATH (joint kernel construction inherited from L2 + Krein-positive completion).**
- Construction of $B_{n_{\max}, N_t}: C_c^\infty(S^3 \times \mathbb{R}) \to \mathcal{O}_{n_{\max}, N_t}^{\mathrm{Krein}}$ as the joint Peter–Weyl × Fourier convolution by the joint kernel from L2.
- Positivity in the Krein-positive sense: $B(f) \ge 0$ on the Krein-positive cone requires the **Krein-positive completion** of van den Dungen 2016 §2. This is structurally available — vdD 2016 explicitly constructs the Krein-positive Hilbert-Schmidt completion for spectral triples — but the **propinquity-compatible version** (where positivity translates into reach-control on the UCP tunneling map) is not standard.
- Contractivity: the Krein-norm version of Young's inequality is a folklore result; needs verification.
- Approximate identity: inherits from L2's joint $\gamma$.
- L3 compatibility: inherits from L3's per-harmonic bound + temporal counterpart.

**Closest published analog.** Hawkins 2000 (Commun. Math. Phys.) for the Riemannian Berezin construction (on Kähler manifolds, which $S^3$ isn't); Paper 38 L4 is the non-Kähler Berezin construction the Riemannian-side proof needed to invent. Lorentzian Berezin is a further extension.

**Risk.** Medium. Berezin construction is the most mechanical of the five lemmas; the Lorentzian version is mostly a re-derivation with the Krein-positive completion replacing the Hilbert-space completion. **Conditional on L2 landing**, L4 is straightforward.

**Confidence: MEDIUM-HIGH.**

### L5 — Latrémolière propinquity assembly

**Riemannian statement.** The tunneling pair $(B_{n_{\max}}, P_{n_{\max}})$ between $\mathcal{T}_{S^3}$ and $\mathcal{T}_{n_{\max}}$ has reach $\le \gamma_{n_{\max}}$ and height $\le \gamma_{n_{\max}}$, giving $\Lambda \le \gamma_{n_{\max}}$.

**Riemannian proof technique.** Assembly via Latrémolière 2017 §4 propinquity-for-metric-spectral-triples definition; reach controlled by L4(c) + L3; height controlled by L4(b) + L3.

**Lorentzian transfer assessment: STRUCTURALLY_BLOCKED unless Lorentzian propinquity exists.**

**This is the hard structural block.** **No published Lorentzian propinquity exists.** Latrémolière 2017/2018/2023/2026, Hekkelman 2022, Hekkelman–McDonald 2024 a/b, Leimbach–van Suijlekom 2024, Toyota 2023, Farsi–Latrémolière 2024/2025: all strictly Riemannian / compact-Hilbert-space. The propinquity definition itself assumes Hilbert-space-based operator norms, positive-definite Lipschitz seminorms, and CP-map tunnels. None of these directly extends to Krein-space settings without re-construction.

**Three candidate workarounds.**

(i) **Wick-rotated propinquity.** Define $\Lambda^L(\mathcal{T}_1^{\mathrm{Krein}}, \mathcal{T}_2^{\mathrm{Krein}}) := \Lambda(\mathcal{T}_1^{\mathrm{Wick}}, \mathcal{T}_2^{\mathrm{Wick}})$ where the Wick-rotated triples have positive-definite metric $g_r$. This is the trivial "lift" and reduces L3 to the Riemannian Paper 38 result on $\mathcal{T}_{S^3 \times S^1_T}$ (with $T \to \infty$ as a separate limit). Pros: zero new mathematics. Cons: reproduces Paper 42's metric-functional-level Wick-rotation correspondence, not a genuine operator-system-level Lorentzian-propinquity statement.

(ii) **Krein-positive propinquity.** Define the propinquity only on the Krein-positive cone (per §1.4). The Krein-positive completion gives a genuine Hilbert space, and the Latrémolière definition applies verbatim. Pros: structurally clean. Cons: covers only positive states; the operator-system on $\mathcal{K}^+$ is a sub-system, not the full $\mathcal{O}_{n_{\max}}^{\mathrm{Krein}}$.

(iii) **Nieuviarts twisted morphism propinquity.** Nieuviarts 2025 (arXiv:2502.18105v3) constructs a morphism Riemannian → pseudo-Riemannian via twisted spectral triples. If the morphism is propinquity-compatible (an open question), then Riemannian propinquity convergence transfers via the morphism. Pros: clean shortcut. Cons: even-dim restriction in Nieuviarts (S³ is odd, requires verification); morphism's propinquity compatibility is unstudied.

**Closest published analog.** None. This is the named **multi-month NCG-math work** that the Sprint L0 audit named as blocker B1.

**Risk.** **HIGH.** If a Lorentzian propinquity exists at all, candidate (ii) is the most likely realistic path. Strong-form L3 (full Krein) requires original mathematics at the L5 level that is not present in the literature.

**Confidence: MEDIUM** that the Krein-positive workaround (ii) is achievable in 6–12 months of focused work; **LOW** that the full Krein L5 is achievable without new propinquity definitions.

### Summary table

| Lemma | Krein-positive verdict | Full Krein verdict | Confidence | Dominant cost |
|:------|:----------------------:|:------------------:|:----------:|:--------------|
| L1' | TRANSFERS_FREELY | NEEDS_NEW_MATH | HIGH | Low: SDP reformulation |
| L2 | NEEDS_NEW_MATH | NEEDS_NEW_MATH | MED-HIGH | **2–4 months: joint compact × non-compact kernel** |
| L3 | NEEDS_NEW_MATH | STRUCTURALLY_BLOCKED | HIGH (Wick) / LOW (genuine) | Low-Medium |
| L4 | NEEDS_NEW_MATH | NEEDS_NEW_MATH | MEDIUM-HIGH | Medium: Berezin + Krein-positive |
| L5 | NEEDS_NEW_MATH | STRUCTURALLY_BLOCKED | MEDIUM (Krein-pos) / LOW (full) | **6–12 months: original NCG-math** |

---

## §3. Krein-positivity restrictions — the realistic shape

Three convergent observations make the Krein-positive-restricted L3 the realistic deliverable:

**(a) L1-tighten (Sprint L1-tighten, 2026-05-16) operates on the Krein-positive cone naturally.** The BW-γ Tomita-Takesaki construction works on the GNS Hilbert-Schmidt space $H_{\mathrm{GNS}} = M_{\dim_W}(\mathbb{C})$ with cyclic vector $\Omega = \rho_W^{1/2}$. This is positive-definite by construction; the Tomita modular operator $\Delta$ has real positive eigenvalues. The Riemannian-side L1-tighten verification at $n_{\max} \in \{2, 3, 4, 5\}$ closed bit-exactly.

**(b) L2-E (Sprint L2-E, 2026-05-16) operates on the wedge KMS state, which is Krein-positive by construction.** The wedge KMS state $\rho_W^L = e^{-K_L^{\alpha, W}}/Z$ is a density matrix with positive trace; its support is the Krein-positive cone of the wedge-restricted Krein space.

**(c) BBB 2018 §5 (the Krein-positive construction).** Bizi–Brouder–Besnard explicitly construct the Krein-positive completion of an indefinite spectral triple at $(m, n) = (4, 6)$ (West-coast Lorentzian); the BBB Hilbert-Schmidt completion is a Type-I factor with tracial Gibbs structure. This is the **published infrastructure** for the Krein-positive restriction.

**What does propinquity look like on the Krein-positive cone?**

Restrict the algebra of operators to $\mathcal{O}_{n_{\max}, N_t}^{+} := \{O \in \mathcal{O}_{n_{\max}, N_t}^{\mathrm{Krein}} : O \cdot \mathcal{K}^+ \subset \mathcal{K}^+\}$. These are the operators preserving Krein-positivity. The Tomita modular Hamiltonian $K_L^{\mathrm{TT}}$ acts on $\mathcal{O}^+$ (Sprint L2-E verifies this bit-exactly). The Latrémolière propinquity definition transfers verbatim to $\mathcal{O}^+$ with the Krein-positive Hilbert-Schmidt operator norm replacing the Hilbert-space operator norm. The propinquity limit object is then the **Krein-positive state space** $\mathcal{P}^+(S^3 \times \mathbb{R}, \mathrm{Krein})$ with a Wasserstein-style distance derived from the Wick-rotated geodesic distance on $g_r$.

**Partial convergence is the realistic deliverable:** $\Lambda^{L,+}(\mathcal{T}_{n_{\max}, N_t}^{\mathrm{Krein}, +}, \mathcal{T}_{S^3 \times \mathbb{R}, \mathrm{Krein}, +}) \to 0$ as $(n_{\max}, N_t) \to (\infty, \infty)$ in a specified joint limit.

This is **not** "full Lorentzian propinquity convergence" but it is the **physics-relevant** statement: the wedge KMS state and the modular Hamiltonian both live on the Krein-positive cone; the L3-Krein-positive claim is exactly the statement that the Riemannian-side WH1 PROVEN result lifts to the Lorentzian setting on the physical (positive) sector.

---

## §4. Open mathematical questions (named explicitly)

### Q1. Joint compact × non-compact Fejér kernel with controllable rate

**Question:** Does there exist a UCP central Schur multiplier $K_{n_{\max}, T_{\max}}$ on $\mathrm{SU}(2) \times \mathbb{R}$ (or $\mathrm{SU}(2) \times S^1_T$) with cb-norm $\le 1$ and joint mass-concentration rate $\gamma_{n_{\max}} \cdot \gamma_{T_{\max}} \to 0$, with explicit asymptotic constants in each factor?

**Positive resolution:** Joint kernel constructed explicitly; cb-norm and Stein-Weiss-style rate verified at finite cutoff and asymptotic limit. **Difficulty: month (compact case with $S^1_T$); 2–4 months (non-compact $\mathbb{R}$).**

**Negative resolution:** No joint kernel with simultaneously bounded cb-norm and controlled rate exists. **Difficulty: unsolved-in-literature** (negative results in this area are typically not published; would require a structural obstruction theorem).

**Risk that this question is the structural block:** MEDIUM-HIGH.

### Q2. Lorentzian Lichnerowicz formula / per-harmonic Lipschitz bound

**Question:** Does the per-harmonic Lipschitz bound $\|[D_L, M_{Y_{NLM, k}}]\|_{\mathrm{op}} \le C \cdot \|\nabla_{g_r} Y_{NLM, k}\|_{L^\infty}$ hold for the Lorentzian Dirac $D_L = i[\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I]$ with a constant $C$ that is bounded above by 1 asymptotically?

**Positive resolution:** Direct computation extends Paper 38's L3 calculation to include the temporal slot; the additional cross-term $\gamma^0 \otimes \partial_t$ contributes a bounded correction. **Difficulty: 1–2 weeks** if Wick-rotated convention accepted.

**Negative resolution:** The Lorentzian Dirac's off-diagonal $\gamma^0$ couplings prevent a per-harmonic bound; only an integrated $L^2$ bound is available. **Difficulty: weeks** to verify; the implication is that L3 transfers only as an $L^2$-Lipschitz bound, which would weaken the propinquity rate.

**Risk:** LOW.

### Q3. Krein-positive Berezin reconstruction

**Question:** Does the joint Berezin map $B_{n_{\max}, N_t}$ preserve Krein-positivity (i.e., $B(f) \cdot \mathcal{K}^+ \subset \mathcal{K}^+$ for $f \ge 0$)?

**Positive resolution:** Direct construction via the joint kernel convolution; the kernel positivity transfers to Krein-positivity by the BBB 2018 Type-I-factor argument. **Difficulty: 1 week** conditional on L2.

**Negative resolution:** The Krein-positive cone is not preserved by convolution; an additional projection step is needed, which may break the approximate-identity rate. **Difficulty: month** to fix.

**Risk:** LOW-MEDIUM.

### Q4. Krein-positive propinquity definition

**Question:** Is there a well-defined Latrémolière-style propinquity on Krein-positive spectral triples that recovers the Riemannian propinquity in the positive-definite limit?

**Positive resolution:** Direct adaptation of Latrémolière 2017 §3 propinquity to the Krein-positive Hilbert-Schmidt setting; the tunnel/reach/height architecture transports verbatim. **Difficulty: 2–4 months** to verify rigorously and identify the correct limit object.

**Negative resolution:** The Krein-positive cone is not closed under the natural tunneling-pair operations; the propinquity definition does not transfer. **Difficulty: unsolved-in-literature** in negative-result form.

**Risk that this is the dominant block:** MEDIUM-HIGH. **This is the keystone question of L3.**

### Q5. Wasserstein–Kantorovich-on-Krein-positive-cone

**Question:** Does Kantorovich–Rubinstein duality on the Krein-positive state space $\mathcal{P}^+(\cdot)$ identify the Connes-distance-on-$\mathcal{O}^+$ with a Wasserstein distance derived from the Wick-rotated geodesic distance on $g_r$?

**Positive resolution:** Direct lift of the Riemannian Kantorovich–Rubinstein argument, using the Krein-positive completion's Hilbert-space structure. **Difficulty: 1–2 months.**

**Negative resolution:** The Krein-positive state space has a non-trivial action of the Krein conjugation $J$ that prevents direct Kantorovich-Rubinstein duality; only a partial identification is available. **Difficulty: months** to characterise the partial duality.

**Risk:** MEDIUM.

### Q6. Joint limit ordering

**Question:** In the joint limit $(n_{\max}, N_t) \to (\infty, \infty)$, does the limit depend on the order of limits, or is it independent (commutative diagram)?

**Positive resolution (commutative):** Both orders give the same propinquity limit; corresponds to the Riemannian-side observation that the spatial truncation and the temporal regulator commute. **Difficulty: 1 month** to verify, conditional on L2.

**Negative resolution (non-commutative):** Limit ordering matters; one ordering produces a propinquity-convergent sequence, the other doesn't. **Difficulty: months** to characterise the right ordering.

**Risk:** LOW-MEDIUM. Most likely commutative based on the L2-E temporal-cutoff insensitivity at finite $n_{\max}$.

### Q7. Nieuviarts morphism applicability to $S^3 = \mathrm{SU}(2)$

**Question:** Does the Nieuviarts 2025 twisted-spectral-triple morphism (arXiv:2502.18105v3) apply to $S^3$, an odd-dimensional manifold?

**Positive resolution:** A careful reading of Nieuviarts §3 confirms the morphism Def 2.2 generalizes from even-dim to odd-dim. **Difficulty: 1–2 weeks** for a scoping pass. If yes, the Sprint L0 audit's named L2-Nieuviarts-scoping shortcut activates and L3 may shorten significantly.

**Negative resolution:** The even-dim restriction is structural and odd-dim requires re-construction. **Difficulty: weeks** to confirm. In this case, L3 reverts to the BBB Krein-space construction with all the original-NCG-math cost.

**Risk:** UNCERTAIN. Paper 42 §3.2 footnote currently reads "NO-GO for odd-dim" but this is not a confirmed structural negative — the L2-Nieuviarts-scoping pass was never executed.

---

## §5. Risk assessment per sub-task

Realistic per-lemma estimates assuming focused work by a researcher with NCG/operator-algebra background:

| Sub-task | Best-case (no block) | Realistic | Worst-case (blocker resolution) | Partial deliverable if blocked |
|:---------|:---------------------|:----------|:--------------------------------|:-------------------------------|
| L1' transfer | 1 week | 2–4 weeks | 2 months | Krein-positive substrate without offdiag SDP |
| L2 joint kernel (compact × $S^1_T$) | 3 weeks | 1–2 months | 3 months | Single-factor kernels only; no joint rate |
| L2 joint kernel (compact × $\mathbb{R}$) | 1 month | 2–4 months | 6 months | $S^1_T$ result only; non-compact deferred |
| L3 per-harmonic bound (Wick gradient) | 1 week | 2 weeks | 1 month | Integrated $L^2$ bound instead of $L^\infty$ |
| L3 per-harmonic bound (Lorentzian gradient) | not started | 2–3 months | unsolved | None |
| L4 Berezin reconstruction (Krein-positive) | 2 weeks | 1 month | 2 months | Riemannian-side Berezin only |
| L5 propinquity assembly (Krein-positive) | 1 month | 3–6 months | unsolved | Wick-rotated propinquity statement |
| L5 propinquity assembly (full Krein) | not started | 6–12 months | unsolved | Krein-positive version only |
| Q7 Nieuviarts scoping | 1 week | 2 weeks | 1 month | Either positive (shortens L3) or negative (confirms BBB path) |

**Aggregated estimates:**

- **L3a (Sprint L3a first move, 1–3 weeks):** Q7 Nieuviarts scoping + L1' Krein-positive substrate write-up + L2 compact-$S^1_T$ kernel scoping. Low-risk, lays groundwork.
- **L3b (compact-temporal restriction, 2–4 months):** Full L1' + L2 + L3 + L4 + L5 on $S^3 \times S^1_T$ with $T$ as a separate limit. Krein-positive restriction. This is the realistic medium-term target.
- **L3c (full Krein, $\mathbb{R}_t$ direction, 6–12 months):** Open-ended NCG-math work. Likely original publishable mathematics; unclear whether it lands or terminates as a partial result.

---

## §6. Realistic timeline and path-of-least-resistance

### Best-case 2–3 month path (compact-temporal restriction)

**Sprint L3a (1–3 weeks, low-risk groundwork):**
- Q7 Nieuviarts scoping pass — does the twist morphism apply to S³? 1–2 weeks. Settles whether L2 needs original construction or just morphism-induced transfer.
- L1' Krein-positive substrate write-up — re-derive Paper 38 L1' on the Krein-positive cone using BBB §5 Type-I-factor completion. 1 week.
- L2 compact-$S^1_T$ kernel scoping — write down the joint Peter–Weyl × Fourier kernel on $\mathrm{SU}(2) \times S^1_T$, compute its cb-norm and joint $\gamma$ rate. 1 week.

**Sprint L3b (4–8 weeks, compact-temporal main course):**
- Full L1'–L5 on $S^3 \times S^1_T$ with periodic boundary conditions on the temporal direction, Krein-positive restriction. 4–8 weeks.
- The output is a theorem: $\Lambda^{L,+}(\mathcal{T}_{n_{\max}, N_t}^{\mathrm{Krein}, +}, \mathcal{T}_{S^3 \times S^1_T, \mathrm{Krein}, +}) \to 0$ at a controllable rate on the Krein-positive cone.

**Sprint L3c (separate sprint, 4–8 weeks):**
- Take the de-compactification limit $T \to \infty$ to recover $\mathcal{T}_{S^3 \times \mathbb{R}, \mathrm{Krein}}$. May require a separate convergence argument.
- Net: a 4–6 month total path to a Lorentzian-propinquity result on the Krein-positive cone of $S^3 \times \mathbb{R}$, at the cost of going through $S^3 \times S^1_T$ as an intermediate.

### Worst-case 12–18 month path

**If L2 joint kernel resists (Q1 negative):** Reformulate L2 as multiple separate kernels with different rates; the propinquity rate becomes a product of rates with explicit cross-terms. Estimated 6 months.

**If L5 Krein-positive propinquity definition resists (Q4 negative):** Construct a custom propinquity-style metric on Krein-positive spectral triples; verify it satisfies Latrémolière-style metric axioms (triangle inequality, separability, etc.). Estimated 6 months.

**If Lorentzian gradient (Q2 negative) demands non-Wick-rotated treatment:** Original Lichnerowicz-on-Lorentzian-manifolds work; estimated 4–6 months.

**Cumulative worst-case:** 12–18 months with significant probability of partial-result termination at month 6–9.

### Path-of-least-resistance: Sprint L3a

**The FIRST CONCRETE TRACK that could be dispatched after this scoping is Sprint L3a as named above (1–3 weeks, three parallel sub-tracks).** L3a does not commit to opening the full L3 stack; it lays groundwork and triggers the natural fork:

1. **If Q7 positive (Nieuviarts applies):** Sprint L3 reduces to morphism-induced transfer. 2–3 months total.
2. **If Q7 negative (Nieuviarts does not apply):** L3 reverts to BBB-construction path; L3b is the natural follow-up.
3. **If L2 scoping reveals a structural block:** Either Sprint L3 terminates as a partial deliverable, or original NCG-math work begins.

L3a is the right move because it is **tractable as a 3-week sprint with three independent sub-deliverables**, each of which produces standalone documentation value regardless of whether L3 itself is opened.

---

## §7. What L3 closes if it succeeds

The L3-Krein-positive-restricted form gives the framework operator-system-level Lorentzian-propinquity convergence on the physical (positive) sector. This is the **continuum-limit lift of Paper 43's finite-cutoff Krein-level closure**.

### §7.1 Direct closures

- **Paper 42 §10 open question O1** (Lorentzian extension at signature (3, 1)) **fully closed** at qualitative-rate level, parallel to WH1 PROVEN on the Riemannian side.
- **Paper 32 §VIII** master Mellin engine case-exhaustion theorem **extends to signature (3, 1)** at the GH-convergence-limit level: every $\pi$ source in any finite chain of Paper 34 projections evaluated on the Krein-positive cone of the continuum Lorentzian Camporesi–Higuchi triple sits in M1, M2, or M3.
- **Paper 38 §6.3 open question (ii)** (Lorentzian extension of the propinquity theorem) closed in the Krein-positive sense.

### §7.2 Multi-focal-composition wall implications

The multi-focal-composition wall (CLAUDE.md §3 entry, Sprint HF synthesis) names five independent observables hitting one wall: **the framework has no native composition theorem for multiple Fock-style projections at once.** L3 does **not** directly address this — the wall is about cross-register two-body coordinate operators (W1a/W1b/W1c), not about Lorentzian extension.

However, L3 closure **structurally extends** the architectural foundation on which composition theorems live. The W2a (multi-loop UV/IR composition) sub-wall has a natural Lorentzian-side reading: Z_2/δm renormalization counterterms in two-loop QED on Dirac-S³ × ℝ at Lorentzian signature. Currently, the LS-8a sprint verdict (Paper 35 Refined Prediction 1) is that bare iterated CC spectral action reproduces the UV-divergent integrand faithfully but cannot autonomously generate counterterms. **With L3 closure, the Lorentzian Krein-positive completion gives a natural setting for renormalisation in the spectral-action language**; whether this enables autonomous counterterm generation is the open follow-up.

### §7.3 What stays external after L3

- **Calibration data (W3):** Yukawas, gauge couplings, mass spectra. L3 is structurally additive at the spectral-triple machinery level; it does not generate calibration data.
- **Cross-manifold W2b** ($\mathcal{T}_{S^3} \otimes \mathcal{T}_{\mathrm{Hardy}(S^5)}$): blocked at NCG-framework level by Coulomb/HO category mismatch (Paper 24 §V). L3 does not touch the cross-manifold question.
- **Multi-loop renormalisation autonomously (LS-8a-renorm)**: the natural Lorentzian-side setting for renormalisation becomes available with L3 closure, but autonomous counterterm generation is a separate question.
- **Two-body Dirac / Breit content (16th projection)**: the operator-level verification of §III.16 at the Krein-positive Lorentzian level becomes possible with L3, but the bound-state Bethe-Salpeter expansion remains a separate technical extension.

### §7.4 Net physics-side impact of L3 closure

Modest at the immediate observable level, structural at the architectural level. L3 closure does **not** open new precision-physics observables (the multi-focal-composition wall observables are bound by W1a/W1b/W1c, not by L3 limitation). What L3 closes is the **structural-completeness** of the framework's NCG content at the operator-system level on the physical (Krein-positive, Lorentzian) sector.

**Compare:** Paper 38 WH1 PROVEN at the Riemannian level shifted the framework's structural-skeleton scope from "alignment with the spectral triple framework at finite cutoff" to "theorem at qualitative-rate GH-convergence level." L3 closure does the analog on the Lorentzian side. It is the structural completion of the operator-system NCG arc.

---

## §8. Net verdict

**GO-WITH-PREREQS** on opening Sprint L3a (the first concrete track) within the next 2 months. The recommended sequence is:

1. **Sprint L3a (1–3 weeks, low-risk):** Q7 Nieuviarts scoping + L1' Krein-positive substrate documentation + L2 compact-$S^1_T$ kernel scoping. Three parallel sub-tracks producing standalone deliverables.
2. **Decision gate at end of L3a:** Based on Q7 outcome and L2 scoping verdict, decide whether to open L3b (compact-temporal restriction) or terminate L3 with a partial deliverable.
3. **Sprint L3b (4–8 weeks, conditional):** Full L1'–L5 on $S^3 \times S^1_T$ Krein-positive cone if L3a verdict is GO.
4. **Sprint L3c (4–8 weeks, conditional on L3b):** De-compactification $T \to \infty$ if L3b lands.

**Recommended decision:** open Sprint L3a as the next concrete track. L3a is **tractable at sprint scale** (3 weeks, three sub-deliverables), produces documentation value regardless of whether full L3 is pursued, and provides the natural decision gate for whether to commit further.

**NO-GO on opening the full L3 stack (Sprints L3a/L3b/L3c) in a single multi-month commitment.** The L2 joint compact × non-compact kernel construction and the L5 Krein-positive propinquity definition are both substantial original-NCG-math work with significant risk of partial termination. The decision-gate architecture (open L3a only, decide L3b/L3c later) is the right risk management.

**Honest framing of the verdict:** The Sprint L2 finite-cutoff closure is the load-bearing physics deliverable for the Lorentzian-extension arc. The four-witness Wick-rotation theorem is closed at the operator-system level at finite cutoff. L3 closure would lift this to the continuum, completing the structural arc, but the physics impact is **modest at the immediate observable level** and **structural at the architectural level**. The 4–8 month investment of opening L3b/L3c should be weighed against alternative directions (state-side dictionary, multi-focal arc extensions, calibration-data W3 investigations) on the strategic question: does the framework benefit more from completing its NCG-math content or from extending its physics-side reach?

The PI's strategic judgment should decide. This scoping memo's role is to clarify what L3 actually is, what it would close, and what it costs — not to push for it. The L3a entry point is named explicitly so that opening it does not commit to the full multi-month stack; a decision gate at end-of-L3a preserves optionality.

---

End of memo.
