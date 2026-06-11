# Strong-Form Lorentzian Propinquity — Math Architecture Scoping Memo

**Date:** 2026-05-22
**Status:** Scoping only (no production code, no proofs, no paper edits).
**Predecessor:** `papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex` (K⁺-restricted weak-form closure, 2026-05-18).
**Goal:** Identify the math architecture for a Latrémolière-style propinquity living on the **full Krein space without K⁺ restriction**, and pick the opening direction for the multi-month sprint sequence (named Sprint L3b-2a onward).

---

## §1. Precise statement of the obstacle

### 1.1 What Latrémolière propinquity requires

The standard Latrémolière propinquity (Latrémolière 2018, Trans. AMS; metric-spectral-triple extension Latrémolière 2017/2023, Adv. Math.) is built on **metric spectral triples** $(\mathcal{A}, \mathcal{H}, D)$ where:

(i) $\mathcal{H}$ is a Hilbert space with a positive-definite inner product $\langle \cdot, \cdot \rangle$ inducing operator norm $\|A\|_{\mathrm{op}} = \sup_{\|\psi\|=1} \|A\psi\|$.

(ii) $D$ is a self-adjoint operator: $D^\dagger = D$.

(iii) The **Lipschitz seminorm** is
$$L(a) := \|[D, a]\|_{\mathrm{op}}$$
defined on a dense $*$-subalgebra of $\mathcal{A}$. This requires $\|\cdot\|_{\mathrm{op}}$ to be a *norm* (positive-definite, satisfying the triangle inequality, sub-multiplicative).

(iv) The propinquity assembly L5 uses tunneling pairs $(B, P)$ where each component is bounded and Lipschitz under $L(\cdot)$; reach and height bounds factor through $L(\cdot)$.

Concretely, Paper 38 §3 Lemma L3 establishes
$$\|[D_{\mathrm{CH}}, M_f]\|_{\mathrm{op}} \le C_3 \cdot \|\nabla f\|_{L^\infty(S^3)}, \qquad C_3 = 1,$$
and Paper 38 §3 Lemma L5 assembles a Latrémolière tunneling pair using this Lipschitz norm. The whole framework rides on $\|\cdot\|_{\mathrm{op}}$.

### 1.2 What breaks on a Krein space

A Krein space (per Paper 44 §2, citing van den Dungen 2016 [vdD16] Prop. 4.1) is $(\mathcal{K}, \langle \cdot, \cdot \rangle_J)$ where $\mathcal{K}$ is a Hilbert space carrying a fundamental symmetry $J$ with
$$J^* = J, \qquad J^2 = +I,$$
and the **Krein indefinite inner product** is
$$\langle \psi, \phi \rangle_J := \langle \psi, J \phi \rangle.$$
The Krein product is *not* positive-definite: it has signature splitting $\mathcal{K} = \mathcal{K}^+ \oplus \mathcal{K}^-$ where $\mathcal{K}^\pm = \ker(J \mp I)$, with $\langle \psi, \psi \rangle_J = +\|\psi\|^2$ on $\mathcal{K}^+$ and $-\|\psi\|^2$ on $\mathcal{K}^-$.

The Lorentzian Dirac is **Krein-self-adjoint** (Paper 45 §2.3, eq. (DL_def)):
$$D_L^\times := J D_L^\dagger J = D_L,$$
which is *not* the same as Hilbert-self-adjointness $D_L^\dagger = D_L$. Concretely, in the GeoVac setting at BBB $(m, n) = (4, 6)$ (Paper 44 §2.3, Bizi-Brouder-Besnard 2018), $D_L = i(\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I_{N_t})$ is Krein-self-adjoint but anti-Hermitian on the underlying auxiliary Hilbert structure (the factor of $i$ in front, plus the anti-Hermitian $\partial_t$).

The two natural objects on the Krein space are therefore:

(a) The **Hilbert-space operator norm** $\|A\|_{\mathrm{op}} := \sup_{\|\psi\|=1} \|A\psi\|$ inherited from the auxiliary Hilbert structure (the positive-definite $\langle \cdot, \cdot \rangle$ before $J$ is applied). This *is* a norm; but the natural adjoint that respects $\|[D_L, \cdot]\|$ is the *Krein adjoint*, not the Hilbert adjoint.

(b) The **Krein adjoint** $A^\times = J A^\dagger J$. The Krein-adjoint pairing $\langle \psi, A \phi \rangle_J = \langle A^\times \psi, \phi \rangle_J$ is well-defined but the resulting Krein-pairing seminorm
$$|A|_J := \sup_{\langle \psi, \psi \rangle_J = 1} |\langle \psi, A \psi \rangle_J|$$
fails to be a norm — it can vanish on non-zero operators (every operator that's symmetric across the K⁺/K⁻ block decomposition with opposite sign has $|A|_J = 0$), and worse, it can be negative on indefinite states.

The strong-form construction must therefore choose: (a) use the Hilbert operator norm and lose contact with the Krein-natural Krein-adjoint structure; (b) use the Krein-natural pairing and lose the norm property; or (c) find a hybrid object that respects both.

### 1.3 How Paper 45 dodged this

Paper 45 §2.6 (`sec:Kplus_definition`) constructs the orthogonal projection
$$P^+ := (I + J_L)/2 : \mathcal{K} \to \mathcal{K}^+$$
and defines the weak-form propinquity (Definition 2.3, `def:weak_form_propinquity`) as
$$\Lambda^L(\mathcal{T}^L_1, \mathcal{T}^L_2) := \Lambda(\mathcal{T}^+_1, \mathcal{T}^+_2),$$
where $\mathcal{T}^+_i = (P^+_i \mathcal{O}^L_i P^+_i, \mathcal{K}^+_i, P^+_i D^L_i P^+_i)$ is the K⁺-restricted **standard** metric spectral triple and $\Lambda$ is the standard Latrémolière construction.

On $\mathcal{K}^+$, the Krein product reduces to the standard Hilbert inner product (since $J|_{\mathcal{K}^+} = +I$). Hence $P^+ D^L P^+$ is genuinely self-adjoint as a Hilbert operator on $\mathcal{K}^+$, the operator norm restricted to $P^+ \mathcal{B}(\mathcal{K}) P^+$ is a genuine norm, and the entire Latrémolière apparatus applies verbatim.

This is honest but limited: it ignores the K⁻ sector entirely. The strong form must handle the full $\mathcal{K}$ natively.

### 1.4 Position in the literature

The audit handed to this scoping (cf. memo header) confirms: **no published Latrémolière-style propinquity on Krein-signature spectral triples exists** as of May 2026. The closest partial pieces are:

- **vdD 2016** (Krein-spectral-triple foundational framework; no metric).
- **Strohmaier 2006** (pseudo-Riemannian distance, temporal-commutative case only).
- **Franco-Eckstein 2014** (algebraic causality, no propinquity).
- **Bizi-Brouder-Besnard 2018** ((m, n) classification, no metric).
- **Devastato-Lizzi-Martinetti 2018** (twisted spectral triples → Krein space, no propinquity).
- **Nieuviarts 2025** (Riemannian-twisted → pseudo-Riemannian via twist morphism; restricted to even-dim manifolds and NOT applicable to GeoVac's S³ at KO-dim 3, per Paper 44 §3.2 footnote).
- **Latrémolière 2025** (most recent propinquity work; strictly Riemannian, analytic-path continuity).
- **Mondino-Sämann 2025 / Minguzzi-Suhr 2022** (synthetic Lorentzian GH on pre-length spaces with causal diamonds — categorically disjoint mathematical object; not operator-algebraic).
- **Krein-space operator theory** (arXiv:1401.3238, 1003.1908, 0904.4095): operator convexity and Lipschitz operator functions on indefinite spaces, relevant *machinery* for §3 below.

The phantom citation "Connes-vS 2021 on Lorentzian propinquity" does **not** exist; Connes-vS 2021 is the spectral-truncation paper (already in the bibliographies of Papers 38/44).

The strong-form construction is therefore a **genuinely open math.OA problem**. This memo's job is to pick the architecture.

---

## §2. Five-lemma audit: what changes from K⁺-weak-form to strong form?

Paper 45 closed the K⁺-weak-form via five lemmas (L1'/L2/L3/L4/L5). Each transports differently to the strong form.

### L1' — Operator-system substrate

**K⁺-weak-form (Paper 45 §3 / Paper 44):** The Lorentzian truncated operator system $\mathcal{O}^L_{n_{\max}, N_t, T}$ is the pure-tensor span of spatial Avery-Wen-Avery multipliers (chirality-doubled, blkdiag(W, W)) and temporal momentum-polynomial diagonal multipliers. It is $*$-closed (under Hilbert adjoint), contains identity, has propagation number prop = 2 on the Weyl-doubled achievable envelope, prop = ∞ on the full $\mathcal{B}(\mathcal{K})$ envelope. The Krein-positive substrate triviality (Paper 44 §5, Proposition 5.1, `prop:krein_positivity_trivial`) says every generator commutes with $J_L$, so $\mathcal{O}^{L,+} = \mathcal{O}^L$ as a vector subspace.

**Strong form:** *No change at the algebra level*. The same operator system $\mathcal{O}^L$ is the natural substrate. What changes is what we do with it: in K⁺-weak-form, we project via $P^+ \mathcal{O}^L P^+$ and lose half the matrix dimension; in strong form, we keep all of $\mathcal{O}^L$ acting on the full $\mathcal{K}$.

**Verdict on L1':** Transports verbatim. The substrate is *already* the right object. The question is the Lipschitz seminorm we put on it.

**Caveat:** If the chosen strong-form seminorm in §3 requires the operator system to be enlarged to include chirality-flipping operators (off-block-diagonal generators) or non-commutative temporal operators (e.g., $D_L$ itself), then L1' must be re-derived on the enlarged substrate. Paper 44 §5 Remark 5.1 (`rem:state_level`) flags both possibilities; neither has been worked out. Honest scope: I cannot rule out at the scoping level whether the natural strong-form seminorm forces such an enlargement.

### L2 — Central spectral Fejér kernel

**K⁺-weak-form (Paper 45 §4.1):** Joint kernel $K^{\mathrm{joint}}_{n_{\max}, N_t, T} = K^{\mathrm{SU}(2)}_{n_{\max}} \otimes K^{U(1)}_{N_t, T}$ with factorized Plancherel symbol. The Plancherel weights act on the operator-system via $\sqrt{2j+1}$-type multiplicative factors. The cb-norm bound $\|S_{K^{\mathrm{joint}}}\|_{\mathrm{cb}} = 2/(n_{\max} + 1)$ is $N_t$-independent (Paper 45 §4.1 Theorem 4.2). The Berezin reconstruction (L4) uses this kernel to convert smooth symbols $f \in C^\infty(S^3 \times S^1_T)$ to operator-system elements.

**Strong form:** The kernel itself is a *geometric / group-theoretic* object — it lives on $\mathrm{SU}(2) \times U(1)_T$ before any Krein structure is invoked. Hence it transports unchanged. What changes is the cb-norm computation: the cb-norm is taken with respect to whichever operator norm the chosen Lipschitz seminorm in §3 induces on the operator system. If the strong-form seminorm uses the *same* underlying Hilbert operator norm (candidates $L_{\mathrm{op}}, L_{\mathrm{twist}}, L_{\mathrm{block}}$ in §3 below), the cb-norm bound transports verbatim. If it uses a fundamentally different norm structure (candidates $L_J, L_{\mathrm{pair}}, L_{\mathrm{spec}}$), the cb-norm must be re-derived.

**Verdict on L2:** Transports verbatim *if* the strong-form seminorm respects the Hilbert operator norm. Otherwise needs partial rederivation. The kernel construction is robust either way.

### L3 — Lipschitz comparison

**K⁺-weak-form (Paper 45 §4.2):** Joint Lichnerowicz bound
$$\|[D_L, a_s \otimes a_t]\|_{\mathrm{op}} = i \|[D_{\mathrm{GV}}, a_s]\|_{\mathrm{op}} \cdot \|a_t\|_{\mathrm{op}} \le C_3^{\mathrm{joint}}(n_{\max}, N_t) \cdot \|\nabla^{\mathrm{joint}} a\|_\infty.$$
The vanishing time-chirality cross term $\{\gamma^0, D_{\mathrm{GV}}\} = 0$ on the chiral basis (Paper 45 §4.2 Lemma 4.3 proof, eq. `L3_struct_id`) makes the temporal direction contribute *zero* to the joint commutator. So $C_3^{\mathrm{joint}} = C_3^{\mathrm{SU}(2)}$ verbatim — the spatial Lichnerowicz constant (Paper 38 §3 L3).

**Strong form:** **This is the dominant change.** The Lipschitz seminorm definition is exactly what shifts; everything downstream depends on which candidate from §3 we pick.

The structural identity $[D_L, a_s \otimes a_t] = i[D_{\mathrm{GV}}, a_s] \otimes a_t$ (Paper 45 eq. `L3_struct_id`) is a property of the *operators*, not of the seminorm. It holds for any candidate seminorm we put on top. What changes is how the seminorm bounds the right-hand side.

For candidates that respect the Hilbert operator norm, the per-pure-tensor bound transports verbatim. For Krein-natural candidates ($L_J$, $L_{\mathrm{pair}}$), the upper-bound estimate must be rebuilt — and may produce a different constant (likely larger than $C_3 = 1$, possibly with a Krein-signature correction $\sim 2$ that survives in the limit).

**Verdict on L3:** Dominant risk. The strong form will need a *new* Lichnerowicz bound matching the chosen seminorm. If the chosen seminorm is "morally Hilbertian on the achievable envelope" (i.e., reduces to operator-norm Lipschitz on operators preserving the K⁺/K⁻ block decomposition — every scalar-multiplier generator does, by Paper 44 Prop. 5.1), then the bound transports. Otherwise it must be rederived.

### L4 — Berezin reconstruction

**K⁺-weak-form (Paper 45 §4.3):** Joint Berezin $B^{\mathrm{joint}} = B^{\mathrm{SU}(2)} \otimes B^{U(1)}$ as a pure-tensor decomposition. Four properties (positivity, contractivity, approximate-identity inheriting the L2 rate, L3-compatibility) verified factor-by-factor.

**Strong form:** The Berezin map itself ($f \mapsto P (K * f) P$ where $K$ is the Fejér kernel and $P$ the truncation projection) is constructed from geometric data and transports verbatim. Properties:
- *Positivity* ($B(f) \ge 0$ for $f \ge 0$): This is *Krein-positivity-sensitive*. In K⁺-weak-form, positivity is in the standard sense $\langle \psi, B(f) \psi \rangle \ge 0$ on $\mathcal{K}^+$. In strong form, the natural notion is *Krein-positive* $\langle \psi, B(f) \psi \rangle_J \ge 0$ for "physical" states (those in $\mathcal{K}^+$ or appropriate Krein-positive cone). Since $B(f)$ commutes with $J$ (Paper 44 Prop. 5.1 generators do, hence so do their images), Krein-positivity reduces to block-positivity on K⁺ ⊕ K⁻ separately. The K⁻ block is unobserved in K⁺-weak-form; in strong form it must be addressed. If $B(f)$ acts the same way on K⁺ and K⁻ (which it does, by the chirality-doubled structure $M = M_W \oplus M_W$), then Krein-positivity holds iff $B(f)|_{\mathcal{K}^+} \ge 0$, which is the K⁺-weak-form result.

- *Contractivity, approximate identity, L3 compatibility:* The contractivity bound $\|B(f)\|_{\mathrm{op}} \le \|f\|_\infty$ transports under the Hilbert operator norm. Under a Krein-natural seminorm the contractivity factor will gain a Krein-norm correction (possibly factor of $\|J\| = 1$, or a Cauchy-Schwarz-type bound — to be worked out).

**Verdict on L4:** The Berezin construction transports verbatim; the *properties* require rederivation with the chosen seminorm. The structural reason this is easier than L3 is that Berezin reconstruction is fundamentally about the *map* $C^\infty \to \mathcal{O}^L$, not about the norm on $\mathcal{O}^L$.

### L5 — Latrémolière propinquity assembly

**K⁺-weak-form (Paper 45 §5):** Standard Latrémolière 2017 tunneling pair $(B^{\mathrm{joint}}, P^{\mathrm{joint}})$. Reach bounded by $\gamma^{\mathrm{joint}}$ (L4 approximate-identity rate); height by L4-contractivity; the propinquity bound assembles.

**Strong form:** L5 is the *assembly* step — it takes the L1'/L2/L3/L4 ingredients and produces the propinquity bound. If those ingredients are in place with the chosen seminorm, the assembly transports.

**The key question is:** Does Latrémolière's 2017/2023 tunneling-pair machinery extend to seminorms that aren't operator-norm Lipschitz? The 2017 paper's Definition 3.1 of a metric spectral triple uses $L(a) = \|[D, a]\|_{\mathrm{op}}$ explicitly; the entire Adv. Math. 2023 framework rides on this. Extending requires either:

(a) Showing the chosen strong-form seminorm $L^{\mathrm{strong}}$ satisfies all properties of Latrémolière's $L(a)$ (lower semicontinuity, $*$-closure, density of the Lipschitz subalgebra, etc.) so the existing machinery applies verbatim with a different definition of $L$.

(b) Rewriting Latrémolière's propinquity construction with the indefinite-metric machinery (e.g., using Krein-adjoint instead of Hilbert-adjoint throughout).

Path (a) is much more tractable. The strong-form seminorm should be designed so that path (a) goes through.

**Verdict on L5:** Inherits from L3 + L4 in the standard case (path a). Genuinely new mathematics only if we go path (b). Honest assessment: path (a) probably works for the cleanest seminorm candidates (§3 below); but in pessimistic-realistic terms, there are probably bookkeeping issues at the reach/height definitions that will eat 4–6 weeks even on path (a).

### Audit summary

| Lemma | Transports verbatim? | New work needed? | Risk |
|:------|:---------------------|:-----------------|:-----|
| L1' (substrate) | ✓ Yes | None (unless §3 enlarges algebra) | LOW |
| L2 (Fejér kernel + cb-norm) | ✓ Yes (kernel); partial (cb-norm) | Cb-norm rederivation only if seminorm changes norm structure | LOW–MED |
| L3 (Lichnerowicz / Lipschitz) | ✗ No | **Dominant new work** — new bound matching chosen seminorm | **HIGH** |
| L4 (Berezin reconstruction) | ✓ (map); partial (properties) | Rederive 4 properties under chosen seminorm | MED |
| L5 (propinquity assembly) | ✓ via path (a); ✗ via path (b) | Latrémolière 2017/2023 framework extension if path (b) | MED–HIGH |

**Net:** L3 is the dominant risk and dominant new work. L5 is secondary risk. L1', L2, L4 are largely mechanical given the §3 choice.

---

## §3. Candidate indefinite-metric Lipschitz seminorm definitions

The strong form pivots on the choice of Lipschitz seminorm on the operator system $\mathcal{O}^L$ when $D = D_L$ is Krein-self-adjoint but not Hilbert-self-adjoint. Below I propose six candidates with pros/cons addressing six axiomatic checks:

(a) Is it a seminorm? (positive homogeneity, triangle inequality)
(b) Is it closed-graph / lower semicontinuous on the Lipschitz domain?
(c) Is it $*$-closed (real on Krein-self-adjoint, or Hilbert-self-adjoint, elements)?
(d) Lichnerowicz upper bound by spectral-gap data (the L3 input)?
(e) Behavior under Berezin reconstruction (L4 compatibility)?
(f) Does it reduce to the standard Latrémolière seminorm when restricted to K⁺, recovering Paper 45 as a corollary?

I flag uncertainty honestly where I can't answer at the scoping level.

### Candidate 1: $L_{\mathrm{op}}(a) = \|[D_L, a]\|_{\mathrm{op}}$ — operator norm in the auxiliary Hilbert structure

**Definition:** Use the Hilbert operator norm from the auxiliary Hilbert space underlying $\mathcal{K}$. $D_L$ is well-defined as a bounded operator (at finite cutoff); the commutator $[D_L, a]$ is a Hilbert operator and has a Hilbert operator norm.

**Pros:**
- Genuine norm; (a), (b) automatic.
- Latrémolière 2017/2023 framework applies essentially verbatim — path (a) of L5 above.
- (f) Reduces to Paper 45 on K⁺ trivially (operator norm restricted to a subspace is operator norm).
- Cb-norm (L2) computation transports.

**Cons:**
- Loses contact with Krein-natural structure. $D_L^\dagger \neq D_L$ in general (it's $D_L^\dagger = -D_L$ if we're naive, or related via $J$); so $[D_L, a]^\dagger = -[D_L^\dagger, a^\dagger] = +[D_L, a^\dagger]$ (treating $D_L$ as anti-Hermitian on the Hilbert structure). The result: $L_{\mathrm{op}}(a^\dagger) = L_{\mathrm{op}}(a)$ — (c) holds, but in a "Hermitian-by-accident" way, not via Krein-self-adjointness.
- (d) The L3 bound likely transports, *but the bound mixes Hilbert and Krein structures* — the natural Lipschitz norm on $C^\infty(S^3 \times S^1_T)$ should treat time and space symmetrically, but the operator-norm bound treats them via the Hilbert auxiliary structure rather than the Lorentzian metric. This is philosophically unsatisfying.
- The K⁻ sector contributes to $L_{\mathrm{op}}(a)$ symmetrically with K⁺, but the *physical* interpretation says K⁻ states are unphysical (Krein-negative norm). So Candidate 1 over-counts K⁻ contributions.

**Verdict:** Workable; closest to "vanilla strong form." Lowest mathematical risk but loses Lorentzian-structural information.

### Candidate 2: $L_J(a) = \|[D_L, a]\|_J$ — Krein-norm-induced seminorm

**Definition:** Define a Krein-norm $\|A\|_J := \sup_{\psi \neq 0} \frac{|\langle A\psi, A\psi \rangle_J|^{1/2}}{|\langle \psi, \psi \rangle_J|^{1/2}}$ or some variant. Then $L_J(a) := \|[D_L, a]\|_J$.

**Pros:**
- Krein-natural; respects the indefinite metric on $\mathcal{K}$.

**Cons:**
- **Not a norm.** $|\langle \psi, \psi \rangle_J|$ vanishes on light-like vectors (any $\psi \in \mathcal{K}^+ \oplus \mathcal{K}^-$ with equal K⁺ and K⁻ components has $\langle \psi, \psi \rangle_J = 0$). The supremum is then infinite or undefined.
- Salvaging by restricting to non-null vectors gives a Finsler-like structure but not a Banach norm. (a) likely fails.
- (b), (c), (d), (e), (f) all unclear or fail.

**Verdict:** Structurally elegant but mathematically broken at the norm level. Reject without further work.

### Candidate 3: $L_{\mathrm{twist}}(a) = \|J [D_L, a] J\|_{\mathrm{op}}$ — J-twisted operator norm

**Definition:** Conjugate the commutator by the fundamental symmetry and take the Hilbert operator norm.

**Pros:**
- (a), (b) automatic — same as Candidate 1.
- $J^2 = +I$ so $\|J A J\|_{\mathrm{op}} = \|A\|_{\mathrm{op}}$ for any $A$ (since $J$ is unitary — $J^* = J$ and $J^2 = I$ imply $J^\dagger J = I$). Hence $L_{\mathrm{twist}}(a) = L_{\mathrm{op}}(a)$ identically.

**Cons:**
- **Identical to Candidate 1.** $J$ is unitary on the auxiliary Hilbert space, so $J$-conjugation is an isometry and doesn't change the operator norm. The "twist" is invisible.

**Verdict:** Reduces to Candidate 1. Not a separate candidate. (This is a useful sanity-check observation: any norm that depends only on the Hilbert structure and ignores the indefinite structure will collapse to the operator norm under $J$-conjugation.)

### Candidate 4: $L_{\mathrm{spec}}(a) = \rho([D_L, a])$ — spectral radius under Krein-self-adjointness

**Definition:** The spectral radius $\rho(A) = \max\{|\lambda| : \lambda \in \mathrm{spec}(A)\}$.

**Pros:**
- Well-defined for any bounded operator.
- For *normal* operators, $\rho(A) = \|A\|_{\mathrm{op}}$. If $[D_L, a]$ is Krein-self-adjoint (which it is when $a$ is — since $D_L$ is Krein-self-adjoint and Krein-adjointness is an anti-automorphism), then we have access to Krein-spectral-theorem machinery (the operator has a Krein-self-adjoint decomposition, though its spectrum may be complex due to non-normality).
- Reduces to operator norm in the normal case.

**Cons:**
- **Spectral radius is not a norm.** $\rho(A) = 0$ for non-zero nilpotent $A$. Triangle inequality fails: $\rho(A + B) \le \rho(A) + \rho(B)$ is *false* in general.
- Krein-self-adjoint operators on Krein spaces can have *complex* spectrum (this is precisely the Krein-vs-Hilbert distinction); the spectral radius captures only the largest modulus, missing geometric information.
- (a) fails (no triangle inequality). (b)–(f) moot.

**Verdict:** Not a seminorm. Reject.

### Candidate 5: $L_{\mathrm{pair}}(a) = \sup_{\|\psi\|_{\mathcal{K}} = 1, \psi \in \mathcal{K}^+} |\langle \psi, [D_L, a] \psi \rangle_J|$ — Krein-pairing seminorm on K⁺

**Definition:** Take the Krein-pairing supremum, but **restrict the supremum to Krein-positive states** $\psi \in \mathcal{K}^+$.

**Pros:**
- Krein-natural.
- (a) Seminorm — yes, by standard sup-of-linear-functionals construction. Triangle inequality and positive homogeneity automatic.
- (b) Lower semicontinuous as a sup over a closed set of states.
- (c) $*$-closure: $\langle \psi, [D_L, a]^\dagger \psi \rangle_J = \overline{\langle \psi, [D_L, a] \psi \rangle_J}$ if $\psi \in \mathcal{K}^+$ (where Krein product = Hilbert product); so $L_{\mathrm{pair}}(a^\dagger) = L_{\mathrm{pair}}(a)$.
- (f) Reduces to the standard Latrémolière seminorm on K⁺ — possibly. The pairing $|\langle \psi, A \psi \rangle|$ is a *weak* seminorm; it lower-bounds $\|A\|_{\mathrm{op}}$ by $\sup |\langle \psi, A \psi \rangle|$ but is generally strictly weaker (numerical range vs operator norm).

**Cons:**
- **It's strictly weaker than the operator norm.** $L_{\mathrm{pair}}(a) \le \|[D_L, a]\|_{\mathrm{op}}$ with strict inequality in general. So this is a *weaker* Lipschitz seminorm.
- (f) Recovery of Paper 45: the standard Latrémolière seminorm on K⁺ uses operator norm, not numerical range. So $L_{\mathrm{pair}}$ restricted to K⁺ is *not* Paper 45's seminorm — it's smaller. The K⁺-weak-form would be *strictly contained in* the strong-form with this candidate; the propinquity bound would be at best as good (i.e., the upper bound would not get tighter).
- (d) Lichnerowicz bound: a numerical-range Lichnerowicz lemma exists but the constant likely differs from operator-norm Lichnerowicz.

**Verdict:** Workable as a *weaker* Lipschitz seminorm. Might be useful if the operator-norm version fails for L3 or L5 — provides a fallback. Honest scope: a seminorm strictly weaker than operator norm leaves room for *multiple* operators with the same Lipschitz norm to disagree on $\|\cdot\|_{\mathrm{op}}$, which may or may not be compatible with Latrémolière's reach/height bookkeeping.

### Candidate 6: $L_{\mathrm{block}}(a) = \max(\|[D_L, a]|_{\mathcal{K}^+}\|_{\mathrm{op}}, \|[D_L, a]|_{\mathcal{K}^-}\|_{\mathrm{op}})$ — K⁺/K⁻ block decomposition

**Definition:** Take the max of the K⁺-block operator norm and the K⁻-block operator norm.

**Pros:**
- **All six axioms work cleanly.** (a) max of two seminorms is a seminorm. (b) lower semicontinuous. (c) $*$-closure on each block separately. (d) Lichnerowicz bound applies on each block — and on each block, $[D_L, a]$ is a Hilbert operator on a Hilbert subspace, so Paper 45's L3 applies *within each block*. The constant $C_3^{\mathrm{joint}}$ transports verbatim because the block restrictions are projections of operators that already commute with $J$ (Paper 44 Prop. 5.1).
- (e) Berezin reconstruction: $B(f)$ commutes with $J$, so $B(f)|_{\mathcal{K}^\pm}$ are Hilbert operators; L4 properties transport block-by-block.
- (f) **Restriction to K⁺ recovers Paper 45 exactly.** On K⁺, the K⁻-block contribution is zero (no operators reach K⁻), so $L_{\mathrm{block}}(a)|_{\mathrm{K^+}} = \|[D_L, a]|_{\mathcal{K}^+}\|_{\mathrm{op}}$, which is Paper 45's seminorm.
- The chirality-doubled structure $M_W \oplus M_W$ (Paper 44 §3.2) means $L_{\mathrm{block}}(a) = \|[D_L, a]|_{\mathcal{K}^+}\|_{\mathrm{op}}$ for every scalar-multiplier generator (the two blocks have identical operator-norm content). So the strong form gives identical Lipschitz numerical values *on the natural operator system* as the K⁺-weak-form does.
- This is, in effect, a "free upgrade" from K⁺-weak to strong: the *numerical answer* is the same on the natural operator system, but the *theorem* is strong-form because the supremum is taken over the entire $\mathcal{K}$, not just K⁺.

**Cons:**
- The "free upgrade" reading is suspicious. If the strong-form propinquity numerically equals the K⁺-weak-form propinquity *on the chirality-doubled operator system*, then arguably the K⁺-weak-form theorem already *was* the strong-form theorem with a different presentation. This is honest content but partially deflates the "multi-month NCG-math problem" framing.
- **However**, $L_{\mathrm{block}}$ is *not* equal to $L_{\mathrm{op}}$. On chirality-flipping operators (off-block-diagonal generators not in the present operator system), $L_{\mathrm{op}}$ may be much smaller than $L_{\mathrm{block}}$ because cross-block content cancels in operator norm but enters block-norms separately. So $L_{\mathrm{block}}$ would yield a different theorem on an enlarged operator system that includes chirality-flipping operators. This is genuinely new content — the strong form is *non-trivial* once the operator system is enlarged.
- This means the strong-form theorem on the *natural* operator system is "free" (modulo the new L3 / L4 bookkeeping), but the *interesting* strong-form theorem requires enlarging the operator system as Paper 44 §5 Remark 5.1 suggests.
- Reach / height bookkeeping in L5 must be re-examined on the block structure; this is straightforward but tedious.

**Verdict:** **The clearest candidate.** Cleanly closed (a)–(f); L3/L4/L5 inherit by block decomposition; recovers Paper 45 as a corollary; and the "free-on-natural-substrate, non-trivial on enlarged substrate" structure is honest about scope.

### Candidate 7 (bonus): $L_{\mathrm{cb}}(a) = \|[D_L, a]\|_{\mathrm{cb}, J}$ — completely bounded Krein norm via Devastato-Lizzi-Martinetti twist

**Definition:** Use Devastato-Lizzi-Martinetti's twisted-spectral-triple framework: the twisted commutator is $[D_L, a]_\rho := D_L a - \rho(a) D_L$ where $\rho$ is a $*$-automorphism (here $\rho(a) = J a J$). Then $L_\mathrm{cb}(a) := \|[D_L, a]_\rho\|_{\mathrm{op}}$ where the operator norm is in the cb-sense (completely bounded).

**Pros:**
- Embeds GeoVac into a published framework (DLM 2018) — gives structural credibility and external machinery.
- The twist $\rho(a) = JaJ$ is identity on the chirality-doubled scalar-multiplier algebra (Paper 44 Prop. 5.1: every generator commutes with $J$). So $[D_L, a]_\rho = [D_L, a]$ for every $a$ in the natural operator system. **Reduces to Candidate 1 on the natural substrate** but with DLM-credibility framing.
- On an enlarged operator system with chirality-flipping generators, the twist is non-trivial and gives genuinely new content.

**Cons:**
- DLM machinery isn't a propinquity construction; it's an action principle. Integrating with Latrémolière 2017/2023 requires bridging work.
- The Krein-Lichnerowicz bound and Berezin properties under twisted commutators are partially worked out in DLM but not in propinquity context.

**Verdict:** Promising for a *future* sprint, especially when enlarging the operator system. Not the right opening choice — the bridge work to Latrémolière is itself a multi-month task.

---

## §4. Recommended candidate

**The recommended seminorm is Candidate 6 ($L_{\mathrm{block}}$).**

### Reasons

1. **Clean axiomatic closure.** All six checks (a)–(f) close cleanly using existing machinery; no genuinely new norm-theoretic work needed.

2. **Recovers Paper 45 exactly.** The K⁺-weak-form result is a strict corollary; no possibility of inconsistency with the closed theorem.

3. **Paper 44 substrate compatibility.** The propagation number prop = 2 result (Paper 44 Theorem 5.3, `thm:prop_envelope_dependent`) is on the *Weyl-doubled achievable envelope*. $L_{\mathrm{block}}$ respects this envelope structure perfectly — block-by-block restriction is exactly the natural way to read the propagation number invariant.

4. **L3 inherits.** The Lichnerowicz bound transports verbatim because the K⁺ and K⁻ blocks each carry a Hilbert spectral triple structure, and on each, the Paper 45 L3 / Paper 38 L3 argument is bit-exact. The constant $C_3^{\mathrm{joint}}$ is the maximum over the two blocks — and on chirality-doubled scalar multipliers, the two blocks give identical numerical values, so the constant transports verbatim.

5. **Strong form is genuinely new on the enlarged operator system.** On the natural operator system, $L_{\mathrm{block}}$ gives the K⁺-weak-form theorem (with a stronger interpretation). On any enlarged operator system including chirality-flipping or non-commutative-temporal generators, $L_{\mathrm{block}}$ gives genuinely new content — and this is where the multi-month NCG-math problem lives.

6. **L5 path (a) works.** Latrémolière 2017/2023 reach/height/propinquity bookkeeping transports because the block-restriction is a contractive completely-positive map. The propinquity bound assembles via standard machinery applied to each block, then maximized.

### What this candidate does NOT do (honest scope)

- It does not give a *categorically* new propinquity on Krein spaces in the sense Paper 45 §7.2 Q1 asked. The candidate is conservative: it leverages the Hilbert structure on each block separately and combines them via max. A more radical candidate would invoke the indefinite structure directly (Candidate 5 was the radical option but failed the seminorm test).

- It does not resolve the open question of what the "right" Lipschitz norm on a fully indefinite Krein space is. It uses a *coordinate* (the J-block decomposition) to reduce to two Hilbert problems.

- The K⁻ block, while now visible, may not carry physical interpretation. The Krein product is negative on K⁻; states there are "negative-norm" in the indefinite metric. The strong-form theorem with $L_{\mathrm{block}}$ treats K⁻ symmetrically with K⁺, which is reasonable mathematically but may not be the physically right answer.

### Fallback candidates

If $L_{\mathrm{block}}$ runs into a bookkeeping issue in L5 (Latrémolière 2017/2023 might not directly support max-of-two-seminorms; this is the named risk), the fallback hierarchy is:

(i) **Candidate 1** ($L_{\mathrm{op}}$). Slightly weaker — collapses the K⁺ and K⁻ blocks via a single operator-norm sup. Loses block structure but stays pure-Hilbert. Most likely to succeed if Candidate 6 has L5-assembly issues.

(ii) **Candidate 7** ($L_{\mathrm{cb}}$ via DLM twist). Higher risk but gives a structurally interesting answer with external published-framework support.

(iii) **Candidate 5** ($L_{\mathrm{pair}}$). Strictly weaker than operator norm; only useful if everything stronger fails.

---

## §5. Roadmap — 3–6 month sprint sequence

Numbering extends the Paper 45 internal taxonomy (L3b-2 was the closure sprint). Strong-form sprints are L3b-2a onward.

### Sprint L3b-2a — Candidate validation + axiom check (4–6 weeks)

**Goal:** Formal verification that Candidate 6 ($L_{\mathrm{block}}$) satisfies all six checks (a)–(f) on the natural operator system.

**Deliverables:**
- Memo formalizing the six checks.
- Symbolic computations at $n_{\max} = 2, 3$ verifying $L_{\mathrm{block}}(a) = \|[D_L, a]|_{\mathcal{K}^+}\|_{\mathrm{op}}$ for every scalar-multiplier generator. (Should be bit-exact by chirality-doubling.)
- Confirmation that Paper 45's K⁺-weak-form is the corollary at $L_{\mathrm{block}}$-Lipschitz on the natural substrate.
- Decision point: proceed with $L_{\mathrm{block}}$ or fall back to Candidate 1.

**Named falsifier:** If the symbolic computations show $L_{\mathrm{block}}(a) \neq \|[D_L, a]|_{\mathcal{K}^+}\|_{\mathrm{op}}$ at any $n_{\max}$ on the natural substrate, the chirality-doubling argument is wrong somewhere; falls back to Candidate 1.

### Sprint L3b-2b — L3 (Lichnerowicz) under $L_{\mathrm{block}}$ (4–8 weeks)

**Goal:** Rederive Paper 45 §4.2 Lemma 4.3 (joint Lichnerowicz) with $L_{\mathrm{block}}$ as the seminorm. Confirm $C_3^{\mathrm{joint}}$ transports verbatim or document the correction.

**Deliverables:**
- Block-by-block Lichnerowicz bound.
- Combined max-bound with explicit constant.
- Comparison with Paper 38 L3 / Paper 45 L3.

**Named falsifier:** Block-by-block bounds give different $C_3$ values for K⁺ and K⁻ blocks. If so, the max is the larger one, and the strong-form bound is strictly worse than the K⁺-weak-form bound. This would be a substantive new finding but still a valid bound.

**Dominant risk:** L3 is the named dominant-risk lemma. The honest probability assessment is ~70% that the bound transports verbatim under $L_{\mathrm{block}}$ on the natural substrate, ~30% that a Krein-structural correction enters.

### Sprint L3b-2c — L4 (Berezin) under $L_{\mathrm{block}}$ (3–5 weeks)

**Goal:** Rederive Paper 45 §4.3 Lemma 4.4 (joint Berezin reconstruction) four properties under $L_{\mathrm{block}}$.

**Deliverables:**
- Block-wise positivity, contractivity, approximate identity, L3 compatibility.
- The structural reason for "free upgrade": Berezin image commutes with $J$, so block restriction is functorial.

**Named falsifier:** Berezin positivity on K⁻ block fails (e.g., $B(f)|_{\mathcal{K}^-}$ negative for $f \ge 0$). This would force restriction to the K⁺ block, collapsing back toward Paper 45. Unlikely given the chirality-doubled structure (which makes K⁺ and K⁻ block-content identical for scalar-multiplier generators).

### Sprint L3b-2d — L5 (Latrémolière propinquity assembly) (4–8 weeks)

**Goal:** Apply Latrémolière 2017 tunneling-pair machinery to the strong-form $\mathcal{T}^L = (\mathcal{O}^L, \mathcal{K}, D_L, L_{\mathrm{block}})$.

**Deliverables:**
- Reach bound under $L_{\mathrm{block}}$.
- Height bound under $L_{\mathrm{block}}$.
- Propinquity assembly theorem.
- Verification that the numerical-panel values transport bit-exactly from Paper 45 §6 (since on the natural substrate, $L_{\mathrm{block}} = L_{\mathrm{op}}$ restricted to K⁺ block, the numerical values should match Paper 45's $\Lambda(2, 3) = 2.0746$, $\Lambda(3, 5) = 1.6101$, $\Lambda(4, 7) = 1.3223$).

**Named falsifier:** Latrémolière's reach/height definition (Latrémolière 2017 §3) uses the operator norm; max-of-two-operator-norms isn't formally a Latrémolière Lipschitz seminorm in the same sense. If this is a real obstruction, fall back to $L_{\mathrm{op}}$ (Candidate 1) and rederive.

**Secondary risk:** L5 is the secondary risk. ~80% chance of clean transport via path (a) of §2 above; ~20% chance of needing partial framework extension.

### Sprint L3b-2e — Strong-form theorem statement and writing (3–4 weeks)

**Goal:** Synthesize L3b-2a through L3b-2d into a coherent strong-form theorem; write Paper 46 (math.OA standalone, sibling of Paper 45).

**Deliverables:**
- Paper 46 draft.
- Bit-exact numerical-panel reproduction of Paper 45 values via the $L_{\mathrm{block}}$ reading.
- Honest scope: the strong-form theorem on the *natural* chirality-doubled scalar-multiplier operator system gives numerically the same propinquity as the K⁺-weak-form, with a strict-strong-form interpretation. The interesting strong form on an *enlarged* operator system (chirality-flipping generators) is the next-sprint target.

### Sprint L3b-2f — Enlargement (optional, 6–10 weeks, deferred)

**Goal:** Enlarge the operator system to include chirality-flipping generators ($\gamma$-matrix-valued multipliers) and re-derive L1' through L5.

**Deliverables (deferred, no commitment):**
- New propagation-number computation on the enlarged substrate.
- L1'/L2/L3/L4/L5 on the enlarged substrate.
- Net theorem: a genuinely-new propinquity convergence that wasn't a Paper 45 corollary.

**Status:** This is where the *interesting* strong-form work lives. Honest scope: deferring this to the second tranche keeps the initial roadmap tractable (4–6 months for the first tranche).

### Roadmap summary table

| Sprint | Duration | Risk | Builds on | Output |
|:-------|:--------:|:----:|:----------|:-------|
| L3b-2a | 4–6 wk | LOW | Paper 44 §5 | Candidate validation |
| L3b-2b | 4–8 wk | **HIGH** | L3b-2a, Paper 38 L3 | Block-L3 |
| L3b-2c | 3–5 wk | MED | L3b-2b, Paper 45 L4 | Block-L4 |
| L3b-2d | 4–8 wk | MED | L3b-2b/c, Latrémolière 2017 | Propinquity assembly |
| L3b-2e | 3–4 wk | LOW | L3b-2a-d | Paper 46 |
| L3b-2f | 6–10 wk | HIGH | Paper 46 | Enlarged operator system |

**Total first tranche (L3b-2a through L3b-2e):** 18–31 weeks ≈ 4.5–8 months.

### Dominant risk assessment

**Dominant risk: L3b-2b L3 transport.**

If the block-by-block Lichnerowicz bound fails to give a unified $C_3$ constant (e.g., the K⁻ block gives a strictly worse constant), the strong-form bound is strictly worse than the K⁺-weak-form bound. This is still a valid theorem but partially deflates the strong-form value proposition. Probability ~30%, mitigated by Candidate 1 fallback.

**Secondary risk: L3b-2d Latrémolière framework extension.**

If $L_{\mathrm{block}}$ as max-of-two doesn't fit Latrémolière 2017/2023 Lipschitz-seminorm axioms cleanly, partial framework extension is needed. Probability ~20%.

**Tertiary risk: scope creep.**

The L3b-2f enlargement is the *interesting* strong form. Resisting the temptation to combine first-tranche (L3b-2a–e) with second-tranche (L3b-2f) is a discipline issue. Recommendation: keep L3b-2f deferred and ship Paper 46 on the natural substrate.

---

## §6. Tractable sub-problems and toy cases

Several smaller units of new work could land before committing the full L3b-2a–e sequence.

### Toy case 1: $N_t = 1$ (Riemannian limit)

Paper 45 §6.2 documents bit-exact Riemannian-limit recovery at $N_t = 1$. The strong-form construction at $N_t = 1$ collapses to a *Riemannian Krein* (Krein structure on a Riemannian spectral triple via the chirality-grading), which is much simpler. Verifying $L_{\mathrm{block}}$ on the $N_t = 1$ slice gives a 2–3 week confidence check before committing to the full L3b-2a–e.

**Deliverable:** Memo + symbolic verification at $n_{\max} \in \{2, 3, 4\}, N_t = 1$. ~2 weeks.

### Toy case 2: SU(2)-only ($S^3$, no temporal slot)

The strong-form question makes sense even on the bare $S^3$ Camporesi-Higuchi spectral triple with chirality grading viewed Krein-style. This is essentially Paper 38 with Krein-block structure layered on. The $L_{\mathrm{block}}$ candidate should give bit-exact Paper 38 results on the natural Weyl-sector operator system. 2–3 week deliverable.

**Deliverable:** Memo + symbolic checks at $n_{\max} \in \{2, 3, 4\}$. ~3 weeks.

### Toy case 3: Abelian special case (replace SU(2) with U(1))

Paper 38 §4 explicitly notes that the SU(2) case is harder than the torus / abelian case. The strong-form Krein construction on a $U(1) \times U(1)_T$ toy case might be analyzable in closed form. This would test the architectural choice without the SU(2)-specific machinery (Avery-Wen-Avery integrals, Peter-Weyl basis).

**Deliverable:** Toy case write-up. 3–4 weeks. Could be done in parallel with L3b-2a.

### Toy case 4: Smallest non-trivial cutoff

At $n_{\max} = 2, N_t = 3$, the Krein space has dimension $16 \times 3 = 48$ (Paper 44 Table after Thm 5.3). $\mathcal{O}^L$ has dimension 42 (same table). All computations are tractable in exact arithmetic via sympy. This is the natural pilot computation for L3b-2a, L3b-2b, L3b-2c.

**Deliverable:** Exact-arithmetic computation of $L_{\mathrm{block}}$ for all 42 generators; comparison with $L_{\mathrm{op}}$. 1–2 weeks.

### Toy case 5: Block-diagonal sanity check

The Paper 44 §5 Proposition 5.1 says every scalar-multiplier generator commutes with $J$, so the K⁺ and K⁻ blocks of every generator are independently computable. Verify by direct computation that:
- $L_{\mathrm{block}}(a) = \|[D_L, a]|_{\mathcal{K}^+}\|_{\mathrm{op}}$ (chirality-doubled, K⁺ = K⁻ blocks equal).
- $L_{\mathrm{op}}(a) = L_{\mathrm{block}}(a)$ on the natural substrate (since both blocks contribute identically).

**Deliverable:** 1-week symbolic verification.

### Recommendation

The natural opening is **Toy case 4 + Toy case 5 combined, completed within Sprint L3b-2a**. This grounds the candidate validation in concrete computations before the four-month commitment.

---

## §7. Honest scope

### What is NOT in this scoping

- **No production code.** No modules in `geovac/`, no test files, no Python implementation of $L_{\mathrm{block}}$.
- **No proofs.** The §3 candidate analyses are heuristic axiom checks at the scoping level; honest "likely works" / "likely fails" reads, not theorems.
- **No paper edits.** No modifications to Paper 45, Paper 44, Paper 32, or any companion paper.
- **No bibliographic deep-dive.** The literature audit was handed in pre-completed (memo header); this scoping uses those findings but did not re-verify them.
- **No commitment to the L3b-2f enlargement.** First-tranche scope is L3b-2a through L3b-2e (4.5–8 months); the genuinely-new content on the enlarged operator system is deferred to a second tranche.

### What the PI should expect from the FIRST sub-sprint (Sprint L3b-2a, 4–8 weeks)

After 4–8 weeks of L3b-2a:

1. A scoping-grade-but-deeper memo formally verifying that Candidate 6 ($L_{\mathrm{block}}$) satisfies axiomatic checks (a)–(f) on the natural chirality-doubled scalar-multiplier operator system, with explicit symbolic computations at $n_{\max} \in \{2, 3\}, N_t \in \{1, 3, 5\}$.

2. A go/no-go decision point: proceed with $L_{\mathrm{block}}$ through L3b-2b/c/d, or fall back to Candidate 1 ($L_{\mathrm{op}}$).

3. Toy-case 4 / toy-case 5 computational verification.

4. Updated risk assessment for L3b-2b (the dominant-risk sprint).

5. An honest "we're not just rewriting Paper 45" check: confirmation that the strong-form theorem on the enlarged operator system (deferred L3b-2f) is genuinely new content, ensuring the multi-month investment is justified.

### What this scoping does *not* commit to

- **Specific theorem statement.** The strong-form propinquity theorem will be written in L3b-2e once L3b-2a–d are closed. Its precise statement depends on which constants survive and what the K⁻-block contribution turns out to be.
- **Specific timeline.** 4.5–8 months is a range, not a deadline. Sprint L3b-2b's dominant risk could extend that to 9–10 months in the bad case.
- **Specific paper venue.** Paper 46 if drafted; could go to J. Geom. Phys. (Paper 45 sibling) or to Adv. Math. (Latrémolière 2017/2023 home). Decision deferred to L3b-2e.

### Confidence assessment

**High confidence:**
- The architectural choice of $L_{\mathrm{block}}$ is the right opening move.
- The natural-substrate-strong-form will recover Paper 45 as a corollary.
- Toy cases 4 + 5 close in 2–3 weeks.

**Medium confidence:**
- L3b-2b L3 will transport verbatim (~70% probability).
- L3b-2d L5 will go path (a) of §2.5 (~80% probability).
- 4.5–8 month timeline for first tranche.

**Low confidence (honest uncertainty):**
- Whether $L_{\mathrm{block}}$ on an *enlarged* operator system (L3b-2f) gives a genuinely-new theorem or whether the enlargement structurally collapses.
- Whether Candidate 6 is *the* right opening or whether Candidate 7 (DLM-twist) is the more durable architecture.
- Whether the K⁻-sector treatment under $L_{\mathrm{block}}$ corresponds to anything physically meaningful in the Lorentzian-extension reading of Papers 43 / 45.

### What I cannot answer at scoping level

- Whether the L3b-2b Lichnerowicz bound gives a clean unified $C_3$ or splits between K⁺ and K⁻ blocks (probably unified, but ~30% probability of split).
- Whether L3b-2d L5 needs a partial Latrémolière 2017/2023 framework extension (~20% probability).
- Whether L3b-2f enlargement produces a categorically-new result vs a quantitative refinement.

These are exactly the questions L3b-2a–d are designed to answer.

---

## Closing note

The strong-form Lorentzian propinquity is genuinely open NCG mathematics. The architectural choice recommended here ($L_{\mathrm{block}}$, Candidate 6) is the conservative path that maximizes leverage on existing machinery (Paper 38, Paper 45, Latrémolière 2017/2023) while preserving honest strong-form interpretation. The interesting categorical content (the enlarged operator system) is deferred to a second tranche to keep the first-tranche commitment tractable.

The dominant scientific risk is the L3 Lichnerowicz bound on the K⁻ block; the dominant scope risk is conflating L3b-2a–e with L3b-2f. Both are manageable with the proposed sprint structure.

**Next action:** Sprint L3b-2a (4–8 weeks), opening with Toy cases 4 + 5 within the first 2 weeks.

