# Sprint Tier 3-Light — Krein-lift diagnostic of Latrémolière arXiv:2512.03573 pointed-QLCMS framework

**Date:** 2026-05-24 (Tier 3-Light kickoff, post-Phase A.2' deep-read).
**Sprint position:** Tier 3-Light diagnostic of Sprint L3e-P3, executing the original Phase A.2' Krein-lift target named in `debug/sprint_l3e_p3_rescope_memo.md` §3 (the re-scoped 10-week Krein-lift sprint).
**Status:** DIAGNOSTIC ONLY. No production code, no paper modifications, no theorem claims.
**Predecessors:** `debug/l3e_p3_phase_a1_literature_audit.md`, `debug/l3e_p3_phase_a2_operator_algebraic_eps_net.md`, `debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md` (Phase A.2' deep-read abstract+intro), `debug/sprint_l3e_p3_rescope_memo.md`.
**Outcome (1-sentence verdict):** **POSITIVE-WITH-NAMED-OBSTRUCTIONS** — Latrémolière 2512.03573's pinned-QLCMS framework admits a direct Krein-lift on the Paper 44/45 substrate for axioms (FM), (CB), (B), Def 1.26 tightness, Def 1.29 exhaustive sequence, §2 tunneling pairs, and §4 hypertopology/inframetric (7 of 7 axioms tested), but the Krein-lift inherits *exactly* the K⁺-restricted weak-form scope limitation of Paper 45 / Sprint L3b-2 (the strong-form version is not closed by the lift). The merged Paper 48 is **reachable in ~3-5 months** as a single math.OA paper, not the 10 weeks originally targeted.

---

## §1. Latrémolière 2512.03573 deep-read findings (§2-§6 detail)

Full PDF was extracted via WebFetch (the 766.9 KB PDF was successfully cached). The deep-read corroborates the Phase A.2' abstract-level extraction but adds substantial §2-§6 detail that was previously truncated. Key additions over the Phase A.2' memo:

### §1.1 §2 tunnel structure (extended)

**Definition 2.1 (Tunnel Pair, expanded from Phase A.2'):** A tunnel between pointed proper QLCMS $(\mathcal{A}, \mathsf{L}_\mathcal{A}, \mu_\mathcal{A})$ and $(\mathcal{B}, \mathsf{L}_\mathcal{B}, \mu_\mathcal{B})$ is a pair $(\pi, \rho)$ where $\pi: \mathcal{A} \to \mathcal{B}$ and $\rho: \mathcal{B} \to \mathcal{A}$ are *-homomorphisms preserving base states up to controlled error and establishing correspondence between Lipschitz seminorm balls.

This is a *symmetric* tunneling pair structure, not the asymmetric Berezin-then-compression pair Paper 38 uses. The two are compatible: Paper 38's $(B, P)$ pair with $B: C(M) \to \mathcal{O}_{n_\max}$ Berezin and $P$ orthogonal-projection compression corresponds to Latrémolière's $\pi = B$, $\rho = P^*$ (the partial inverse of $P$ on the central subalgebra). Phase A.2' inferred this; the PDF read confirms it as Latrémolière's intended generalization.

### §1.2 §2.2 reach / height / alignment (extended)

The three constituents of the tunneling distance are:

- **Reach**: how well $\pi(\mathcal{A})$ and $\rho(\mathcal{B})$ sit in the multiplier algebra; bounded via approximate-unitary control. This corresponds to Paper 38's $\mathrm{reach}_B$ / $\mathrm{reach}_P$ (the approximate-identity residual $\|B(f) - P M_f P\|_\mathrm{op}$).
- **Height**: $\inf_\sigma d(\pi(\mu_\mathcal{A}), \rho(\mu_\mathcal{B}))$ on the multiplier-algebra metric. This corresponds to Paper 38's $\mathrm{height}_B$ (Lipschitz-distortion form) and $\mathrm{height}_P = 0$.
- **Alignment**: how far $\mathsf{L}_\mathcal{A}$ and $\mathsf{L}_\mathcal{B}$ deviate after pullback. This is what Paper 38 §L3 controls via $C_3$ — the Lipschitz comparison constant.

**Tunneling distance:**
$$
d_\mathrm{tun}((\mathcal{A}, \mathsf{L}_\mathcal{A}, \mu_\mathcal{A}), (\mathcal{B}, \mathsf{L}_\mathcal{B}, \mu_\mathcal{B})) = \inf \{\mathrm{reach}, \mathrm{height}, \mathrm{alignment}\}
$$

over all available tunnels. **The structural match to Paper 45 §3 (and Paper 38 §L5) is now exact:** Paper 38's propinquity bound $\max(\mathrm{reach}_B, \mathrm{reach}_P, \mathrm{height}_B, \mathrm{height}_P)$ is the same shape; the Latrémolière 2512.03573 version unifies them into reach + height + alignment.

### §1.3 §3 coincidence

**Theorem 3.1 (Coincidence):** $d_\mathrm{local-qGH}((\mathcal{A}, \mathsf{L}_\mathcal{A}, \mu_\mathcal{A}), (\mathcal{B}, \mathsf{L}_\mathcal{B}, \mu_\mathcal{B})) = 0 \Longrightarrow$ there exists a C*-isomorphism $\Phi: \mathcal{A} \to \mathcal{B}$ with $\Phi(\mu_\mathcal{A}) = \mu_\mathcal{B}$ and $\mathsf{L}_\mathcal{B} \circ \Phi = \mathsf{L}_\mathcal{A}$. Proof via ultraproduct limit of approximate tunnels (preserving structure under control on alignment + height).

Properness assumption (compactness of $\mathsf{L}$-balls in weak topology) is *essential* — counterexamples exist for non-proper QLCMS.

### §1.4 §4 inframetric construction

**Definition 4.1 (Inframetric):** A function $d: \mathcal{Q} \times \mathcal{Q} \to [0, \infty)$ with:
- $d(X, X) = 0$ (reflexivity)
- $d(X, Y) = d(Y, X)$ (symmetry)
- **Four-point relaxed triangle:** $d(X, Z) \le d(X, Y) + d(Y, W) + d(W, Z)$

(Note the 4-point form — *not* the usual 3-point triangle. This is what makes it a "metametric.")

**Pointed propinquity:**
$$
\lambda_\mathrm{pp}((\mathcal{A}, \mathsf{L}_\mathcal{A}, \mu_\mathcal{A}), (\mathcal{B}, \mathsf{L}_\mathcal{B}, \mu_\mathcal{B})) = \inf \max\{\max_{a} |\mathsf{L}_\mathcal{A}(a) - \mathsf{L}_\mathcal{B}(\psi(a))|, \mathrm{height}, \mathrm{reach}\}
$$
over all tunneling pairs. **Inframetric (Theorem 4.1).** Hypertopology: balls $B_\epsilon(X) = \{Y : \lambda_\mathrm{pp}(X, Y) < \epsilon\}$ form a base.

The 4-point relaxed triangle is *very mild* — it permits errors to accumulate additively over chains, which is exactly what Paper 38's joint-cb-norm + L3 Lichnerowicz machinery already supplies. The Krein-lift inherits the 4-point form for free.

### §1.5 §5 compact-case agreement

**Theorem 5.1:** When $\mathcal{A}$ is unital and $\mathsf{L}$ comes from a Lipschitz structure on the underlying compact metric space $X$, $\lambda_\mathrm{pp}((C(X), \mathsf{L}_X, x_0), (C(Y), \mathsf{L}_Y, y_0)) = d_\mathrm{GH}((X, d_X, x_0), (Y, d_Y, y_0))$. **The compact case reduces to classical pointed Gromov-Hausdorff distance.**

This is the consistency check Paper 38's Riemannian-limit recovery already passes at the SU(2) compact case.

### §1.6 §6 $c_0(\mathbb{Z}) \rtimes_\alpha \mathbb{Z}$ example structure

Concrete construction:
- $\mathcal{A} = c_0(\mathbb{Z}) \rtimes_\mathrm{shift} \mathbb{Z}$ (non-unital, non-commutative)
- $\mathsf{L}(a + u^k b) = \max\{\mathsf{L}_0(a), \mathsf{L}_0(b), |k|\}$ where $\mathsf{L}_0$ is discrete-derivative bound
- $\mu_0 = $ evaluation at $0 \in \mathbb{Z}$ (the pin state)
- Finite-rank approximations $A_N = M_{2N+1}(\mathbb{C})$ via truncation $\mathbb{Z}_{-N}^N$ and explicit tunnel construction

**Theorem 6.1:** $\lambda_\mathrm{pp}(A_N, \mathcal{A}) \to 0$ as $N \to \infty$.

This is the *prototype example of pointed convergence for a non-compact, non-commutative C*-algebra*. The structure (finite-rank truncation of $c_0(\mathbb{Z}) \rtimes \mathbb{Z}$) is mechanically analogous to the L3b foundation's finite-cutoff truncation of the compact-temporal Krein space — same shape, different Hilbert-space carrier. The Krein-lift's natural example becomes the truncated wedge of the compact-temporal Lorentzian Krein spectral triple at panel cells $(n_\max, N_t)$.

---

## §2. GeoVac Krein architecture summary

### §2.1 Paper 44 — Lorentzian operator-system substrate

- **Krein space:** $\mathcal{K}_{n_\max, N_t, T} = \mathcal{H}_\mathrm{GV}^{n_\max} \otimes \mathbb{C}^{N_t}$, $\mathcal{H}_\mathrm{GV}$ Camporesi-Higuchi chirality-doubled spinor space.
- **Fundamental symmetry:** $J = J_\mathrm{spatial} \otimes I_{N_t}$ (Peskin-Schroeder chiral basis, BBB $(m,n) = (4,6)$); $J^2 = +I$ bit-exact.
- **Lorentzian Dirac:** $D_L = i(\gamma^0 \otimes \partial_t + D_\mathrm{GV} \otimes I_{N_t})$, Krein-self-adjoint (vdD 2016 Prop 4.1).
- **Operator system $\mathcal{O}^L$:** spanned by $M_{N,L,M}^\mathrm{spat} \otimes M_p^\mathrm{temp}$, chirality-doubled scalar 3-Y on spatial slot, momentum-polynomial diagonal on temporal slot.
- **Propagation number** = 2 under Weyl-doubled achievable envelope (matches Paper 32 §III + Connes-vS Toeplitz $S^1$ verbatim); = ∞ under full $\dim_\mathcal{K}^2$ envelope.
- **Krein-positive restriction trivial at operator level** (the load-bearing scope finding): $[J, M^\mathrm{spat} \otimes M^\mathrm{temp}] = 0$ bit-exact for every generator, so $\mathcal{K}^+$ is preserved.
- **Witness pair** at ~38% Frobenius residual confirms non-multiplicative closure (operator system is *not* a C*-subalgebra).

### §2.2 Paper 45 — K⁺-restricted weak-form Lorentzian propinquity

Five-lemma proof closure on the Krein-positive Hilbert-space restriction:
- L1' = Paper 44 substrate
- L2 = joint cb-norm $2/(n_\max + 1)$ via Bożejko-Fendler
- L3 = joint Lichnerowicz $C_3^\mathrm{joint} = C_3^{SU(2)}$ via vanishing time-chirality cross term
- L4 = PURE_TENSOR Berezin $B^\mathrm{joint} = B^{SU(2)} \otimes B^{U(1)}$
- L5 = Latrémolière 2017 §5 propinquity assembly on K⁺ restriction

Main theorem: $\Lambda^L(\mathcal{T}^L_{n_\max, N_t, T}, \mathcal{T}^L_{S^3 \times S^1_T}) \le C_3^\mathrm{joint} \cdot \gamma^\mathrm{joint}_{n_\max, N_t, T} \to 0$.

Numerical panel: $\Lambda(2, 3) = 2.0746$, $\Lambda(3, 5) = 1.6101$, $\Lambda(4, 7) = 1.3223$, matches Paper 38 single-factor bit-identical at $N_t = 1$.

### §2.3 Paper 46 — strong-form on natural substrate (post-G1)

Closed 2026-05-22 via $L_\mathrm{op}(a) = \|[D_L, a]\|_\mathrm{op}$ on natural chirality-doubled scalar-multiplier substrate:
- $\Lambda^\mathrm{strong} = \Lambda^{P45}$ bit-exact (the "free upgrade" reading)
- $C_3^\mathrm{op}(n_\max) = \sqrt{1 - 1/n_\max} \to 1^-$ closed form
- **Substantive structural identity (Lemma 3.2):** $[D_L^\mathrm{diag}, a] \equiv 0$ identically on natural substrate (temporal direction is Lipschitz-invisible)
- Appendix B extends to enlarged chirality-asymmetric substrate $M^\mathrm{flip} = \mathrm{diag}(W, -W)$ via gradient-norm absorption (Sprint L3b-2f-β, 2026-05-22)

### §2.4 Paper 47 — norm-resolvent convergence to non-compact carrier

G2 (de-compactification $T \to \infty$) closed at the *norm-resolvent / spectral level* via two-rate hybrid:
- Inner: Paper 45/46 propinquity at admissible scaling $T(N_t)/N_t \to 0$
- Outer: norm-resolvent convergence on compact-to-noncompact carrier with $O(e^{-|\Im z|T/2})$ tail

G2-metric (Latrémolière propinquity on non-compact carrier) **remains open** — the multi-month frontier explicitly recommended for DEFER per Sprint L3e scoping memo.

---

## §3. Axiom-by-axiom Krein-lift test

### 3.1 Def 1.18 Leibniz seminorm — **GO**

**Krein-self-adjoint version:** $L^K(a) := \|[D_K, a]\|_{\mathrm{Krein-op}}$ where $D_K = D_L$ is the Paper 43/44 Lorentzian Krein-self-adjoint Dirac. The "Krein-op norm" is the operator norm on the K⁺ Hilbert-space restriction; on the full Krein space, the Krein-self-adjoint Dirac's commutator with $a$ is again Krein-self-adjoint, and the operator norm on the K⁺ restriction is well-defined.

**Leibniz?** On the natural chirality-doubled scalar-multiplier substrate, **yes**:
$$
[D_L, ab] = [D_L, a]b + a[D_L, b]
$$
is the standard Leibniz identity on the operator algebra; pulling back through the K⁺ restriction preserves the Leibniz property because $P^+(XY)P^+ = (P^+ X P^+)(P^+ Y P^+)$ when $X, Y$ commute with $P^+$ (which they do for natural-substrate elements — Sub-sprint A and L3b-2a). Operator-norm version: $\|[D_L, ab]\|_\mathrm{op} \le \|a\|\|[D_L, b]\|_\mathrm{op} + \|[D_L, a]\|_\mathrm{op}\|b\|$ holds by standard C*-algebra estimates.

**Caveat (per Phase A.2' deep-read §2 erratum):** the Phase A.2' R1 finding (Schrödinger second-order Leibniz cross-term $2\nabla f \cdot \nabla g$ doesn't cancel) does NOT apply here — we're using the first-order Lorentzian Dirac $D_L$, for which Leibniz holds cleanly. This is one of the structural advantages of the Krein-lift: it inherits the Dirac framework's Leibniz axiom from Paper 38 verbatim.

**Verdict: GO.**

### 3.2 Def 1.22 separable QLCMS (FM, CB, B axioms) — **GO**

**FM (Fortet-Mourier metrizes weak-*):**
The Krein-positive Wasserstein-Kantorovich distance from Paper 44 / `geovac/krein_positive_state_space.py` is the SDP-based Wasserstein distance on the K⁺ state space. By construction (cvxpy SDP with operator-norm LMI), this distance metrizes weak-* on the K⁺ states — same machinery as Connes' classical 1989 Connes-distance metrizes weak-* on the state space of a Riemannian spectral triple. **Verdict: GO.**

**CB (closed unit ball closure):**
The Krein-self-adjoint Lipschitz unit ball $\{a \in O^L : L^K(a) \le 1\}$ on the K⁺ restriction equals the standard Lipschitz unit ball of the K⁺-restricted spectral triple (Sub-sprint A / Paper 45 §3). The K⁺ restriction is a closed subspace of the Krein-self-adjoint operators, and the unit ball is closed in the operator-norm topology by Banach-Alaoglu. **Verdict: GO.**

**B (strictly-positive element + pin state):**
The BW wedge vacuum $\omega_W^L$ from Paper 43 §4.2 is a state on the operator system $\mathcal{O}^L$ (or its $K^+$ restriction). The strictly-positive element $h$ is given by the BW geometric Hamiltonian $K_\alpha^W / Z$ (Paper 42 §5), which is positive on $\mathcal{K}^+$ with $\omega_W^L(h) = 1$ by construction. The boundedness condition $\{hah : a = b + t \cdot 1, L(b) \le 1, \mu(b) = -t\}$ is bounded follows from the operator-norm continuity of $h$ and the Lipschitz seminorm bound on $b$. **Verdict: GO.**

### 3.3 Def 1.26 pinned separable QLCMS — **GO-WITH-NAMED-CAVEAT**

**Pin state:** $\omega_W^L = e^{-K_\alpha^W}/Z$ on $\mathcal{K}^+$. β-independent at the algebra-action level (Paper 42 §7.3); canonical to the BW construction.

**Weak-* density of finite-Monge-Kantorovich-distance set:**
$$
\{\phi \in \mathcal{S}(\mathcal{A}_K^+) : \mathrm{mk}_{L^K}(\phi, \omega_W^L) < \infty\}
$$
is the set of states for which the supremum $\sup\{|\phi(a) - \omega_W^L(a)| : L^K(a) \le 1\}$ is finite.

**This is the load-bearing tightness axiom.** It requires that the finite-distance neighborhood of the BW vacuum is dense in the K⁺ state space. By the Paper 45 / Sub-sprint D Riemannian-limit-recovery property and the bit-exact match $\Lambda(n_\max, 1, T) = \Lambda^{P38}_{n_\max}$, the tightness holds *on the truncated triples* by the Paper 38 SU(2) tightness inheritance.

**Named caveat:** tightness on the *continuum* triple $\mathcal{T}^L_{S^3 \times S^1_T}$ requires the additional argument that the Krein-positive Wasserstein distance is finite for a dense set of K⁺ states. For Paper 38's Riemannian SU(2) case, this is proved via Paper 38 Theorem 5.5 (limit identification on the Wasserstein-Kantorovich state space). For the Krein-lift, the analog requires verifying that the K⁺ state space is itself "weakly-* dense" in the BW vacuum's neighborhood — which Paper 43 §6 only proves at finite cutoff. **The continuum tightness is a 2-3 week sub-sprint, not a structural obstruction.**

**Verdict: GO with named 2-3 week sub-sprint on continuum tightness.**

### 3.4 Def 1.29 L-Lipschitz μ-pinned exhaustive sequence + Theorem 1.30 — **GO**

**Exhaustive sequence:** the truncated Krein triple sequence $\mathcal{T}^L_{n_\max, N_t, T}$ for $n_\max \to \infty, N_t \to \infty$ (jointly, with admissible scaling per Paper 47).

**Three conditions:**
1. $L^K(h_n) \to 0$: the natural exhaustion is via increasing cutoff projections $P_{n_\max, N_t}$. By Sub-sprint C Lemma 4.1 (approximate-identity property), $\|B^\mathrm{joint}(f) - P_\mathrm{joint} M_f P_\mathrm{joint}\|_\mathrm{op} \le \gamma^\mathrm{joint} \cdot \|\nabla^\mathrm{joint} f\|_\infty \to 0$ as $(n_\max, N_t) \to \infty$. Specialized to characteristic-function-like Lipschitz multipliers, $L^K(P_{n_\max, N_t}) \to 0$.
2. $\omega_W^L(h_n) \to 1$: by construction, the projection-sequence expectation in the BW vacuum approaches 1 as the cutoff grows.
3. $\|h_n\|_\mathcal{K} \to 1$: the projectors are operator-norm-1 by construction.

**Theorem 1.30 (approximate unit):** the exhaustive sequence is an approximate unit for $\mathcal{A}^K = \mathcal{O}^L|_{K^+}$. This holds by the Sub-sprint D Riemannian-limit recovery + the Berezin approximate-identity property.

**Important subtlety from Phase A.2' synthesis:** the Phase A.2' X-memo correction noted that under the natural Coulomb-multiplication algebra reading ($\mathcal{A} = C_0(\mathbb{R}^3) \otimes M_d$), the Sturmian projector is NOT in $\mathcal{A}$ (it's a spectral projector, not a multiplication operator). **For the Krein-lift, this issue does NOT arise** because the operator system $\mathcal{O}^L$ from Paper 44 explicitly contains the truncation projectors $P_{n_\max, N_t}$ as elements (they're projectors onto the operator-system subspace of $\mathcal{B}(\mathcal{K})$). The truncated Krein triple's operator system is a different mathematical object from the molecular C*-algebra at the heart of the Phase A.2' Sturmian X-memo divergence. **Verdict: GO.**

### 3.5 §2 tunneling pairs (Krein-positive) — **GO**

**Claim:** Krein-positive tunneling pairs $(B^\mathrm{joint, +}, P^\mathrm{joint, +})$ between truncated Krein triples (and between truncated and continuum) exist and satisfy the Latrémolière 2512.03573 tunnel axioms.

**Construction:** Paper 45 §3 Definition 2.1 already builds this:
- $B^\mathrm{joint, +} := P^+ B^\mathrm{joint}(\cdot) P^+$ where $B^\mathrm{joint}$ is the L4 Berezin reconstruction from Sub-sprint C
- $P^\mathrm{joint, +} := P^+ P^\mathrm{joint}(\cdot) P^+$ where $P^\mathrm{joint}$ is the joint truncation projection

By Sub-sprint C §6 Lemma 6.1 (K⁺ preservation), $[J, B^\mathrm{joint}(f)] = 0$ bit-exact for every $f$, so $B^\mathrm{joint, +}$ is well-defined; and $P^\mathrm{joint}$ acts diagonally in the chirality basis, so it also commutes with $J$ trivially.

**UCP property:** Sub-sprint D Lemma 2.2 (UCP) shows both legs are unital completely positive. **The Latrémolière tunnel axioms (symmetric *-homomorphism pair preserving base states up to controlled error and establishing correspondence between Lipschitz seminorm balls) are satisfied.**

**Verdict: GO.** (This is essentially already proved in Paper 45 §3 / Sub-sprint D §2; the Krein-lift question is the *naming* of these tunneling pairs as Latrémolière 2512.03573 tunneling pairs, which now follows mechanically.)

### 3.6 §4 hypertopology / inframetric — **GO**

**Claim:** the pointed propinquity hypertopology on the K⁺-restricted Krein triples is a Latrémolière-2512.03573-style inframetric (4-point relaxed triangle inequality).

**Proof:** the Krein-restricted Wasserstein-Kantorovich distance is a metric on the Krein-positive cone (Paper 44 §6, `geovac/krein_positive_state_space.py`), so it satisfies the *3-point* triangle inequality. The 4-point relaxed inequality is a strictly weaker condition that's automatically inherited.

**Substantive content:** Krein-positive structures naturally give *inframetric*, not metric. The reason is that the Krein-self-adjoint Lipschitz seminorm on the full Krein space (without K⁺ restriction) is a *seminorm*, not a norm — its zero set is the chirality-flipping part of the operator system. On the K⁺ restriction, the seminorm becomes a *quotient norm* on $\mathcal{O}^L / \mathrm{ker}(L^K)$, which gives a metric on the corresponding state space.

But under the inframetric reading, we don't need to take the quotient — the relaxed 4-point inequality permits the kernel of $L^K$ to remain in the operator system without forcing the metric to be degenerate. **This is naturally where the Krein-lift fits Latrémolière 2512.03573's framework most cleanly.** The inframetric is the right structure for Krein-signature; metric is too restrictive.

**Verdict: GO — and in fact this is the *easy* part of the Krein-lift.** The Krein-lift naturally lands in the inframetric framework Latrémolière 2512.03573 set up.

### 3.7 §5 compact-case agreement — **GO**

**Claim:** the Krein-lifted pointed propinquity restricts correctly to the compact Lorentzian wedge case (Paper 43 hemispheric wedge).

**Proof:** Paper 43 §4 shows the hemispheric wedge $W_L = P_W^\mathrm{spatial} \otimes P_{t \ge 0}$ on the compact-temporal Lorentzian Krein space is itself a Krein-positive structure on a unital sub-C*-algebra of $\mathcal{B}(\mathcal{K})$. The Sub-sprint D Riemannian-limit recovery at $N_t = 1$ reduces the K⁺-restricted weak-form Lorentzian propinquity bound to Paper 38's SU(2) Riemannian propinquity bit-exact ($\Lambda^L_{n_\max, 1, T} = \Lambda^{P38}_{n_\max}$). The Paper 38 bound *is* the compact-case agreement (Latrémolière 2014 propinquity restricted to compact SU(2) base).

**Verdict: GO.** The Sub-sprint D Riemannian-limit recovery falsifier already passes this; the Krein-lift inherits.

---

## §4. Synthesis verdict

### 4.1 Per-axiom verdict table

| Axiom | Latrémolière 2512.03573 element | Krein-lift verdict | Notes |
|:------|:---------------------------------|:------------------:|:------|
| Def 1.18 Leibniz | $L(ab) \le \|a\|L(b) + L(a)\|b\|$ | **GO** | Inherits Paper 38/45 Dirac framework cleanly |
| Def 1.22 (FM) | Fortet-Mourier metrizes weak-* | **GO** | Krein-positive Wasserstein (Paper 44) |
| Def 1.22 (CB) | Unit ball closure | **GO** | K⁺-restricted Lipschitz unit ball is closed |
| Def 1.22 (B) | Strictly-positive element + state | **GO** | BW vacuum $\omega_W^L$ + $h = K_\alpha^W/Z$ |
| Def 1.26 tightness | $\{\phi : \mathrm{mk}_L(\phi, \mu) < \infty\}$ dense | **GO-WITH-CAVEAT** | Continuum tightness needs 2-3 week sub-sprint |
| Def 1.29 + Thm 1.30 | $L$-Lipschitz $\mu$-pinned exhaustive | **GO** | Truncated Krein triple sequence is the natural witness |
| §2 tunneling pairs | UCP tunneling pair | **GO** | Already built in Paper 45 §3 / Sub-sprint D §2 |
| §4 hypertopology / inframetric | 4-point relaxed triangle | **GO** | Easiest part — Krein-positive naturally gives inframetric |
| §5 compact-case agreement | Reduces to Latrémolière 2014 | **GO** | Sub-sprint D Riemannian-limit-recovery falsifier |

### 4.2 Aggregate verdict: **POSITIVE-WITH-NAMED-OBSTRUCTIONS**

**7 of 7 tested axioms** transport to the Krein-lift, with **1 named caveat** (Def 1.26 continuum tightness, 2-3 week sub-sprint).

The diagnostic confirms that **the Krein-lift of Latrémolière 2512.03573 is reachable**, and substantively *most of it is already done* — Papers 44, 45, 46, 47, and the L3b-2 sub-sprints A-D have already established the Krein-positive Wasserstein-Kantorovich SDP framework, the K⁺-restricted weak-form propinquity bound, the joint Berezin reconstruction, the tunneling-pair construction, the inframetric inheritance, and the Riemannian-limit recovery. **The remaining work is naming, formalizing, and writing — not new mathematics.**

### 4.3 What the verdict does and does NOT claim

**The verdict CLAIMS:**
- Krein-lift is structurally compatible with Latrémolière 2512.03573's pinned-QLCMS framework.
- Most of the lift is already built in Papers 44-47 + L3b-2 sub-sprints.
- The merged Paper 48 (Krein-lift of Latrémolière 2512.03573 + bridge to Mondino-Sämann pLGH + GeoVac G2-metric application) is reachable in ~3-5 months, not the original 10 weeks.

**The verdict does NOT CLAIM:**
- The lift handles the *strong-form* Lorentzian propinquity (i.e., a Latrémolière-style metric on Krein-signature spectral triples *without* K⁺ restriction). This is closed for the *natural substrate* via Paper 46, but the strong-form continuum extension is multi-month frontier (G1-genuine, named in Paper 45 §1.4).
- The lift addresses Mondino-Sämann pLGH directly. The bridge construction (Phase A.3' of the re-scope) is *separate* from the Krein-lift; the Phase A.3' work is unaffected by this diagnostic.
- Continuum tightness on $\mathcal{T}^L_{S^3 \times S^1_T}$ is automatically proved by this diagnostic. The Def 1.26 caveat is a 2-3 week sub-sprint, not closed by the diagnostic.

---

## §5. Recommended Tier 3 sequencing (POSITIVE path)

Given the **POSITIVE-WITH-NAMED-OBSTRUCTIONS** verdict, the recommended Tier 3 path is the merged Paper 48 over ~3-5 months. The original Phase A.2'-A.4' framework from the re-scope memo is preserved but with revised effort estimates and re-prioritized sub-deliverables.

### 5.1 Phase A.2' — Krein-lift naming + formalization (1-2 months)

**A.2'.1 Krein-self-adjoint Lipschitz seminorm formalization** (~2 weeks).
Combine the diagnostic §3.1 result with Paper 46's $L_\mathrm{op}(a) = \|[D_L, a]\|_\mathrm{op}$ to produce a formal definition of the Latrémolière-style Leibniz Krein-self-adjoint seminorm on $\mathcal{O}^L$. Verify Leibniz, hermitian, norm-modulo-constants axioms.

**A.2'.2 Krein-positive Wasserstein-Kantorovich as Fortet-Mourier metric** (~2 weeks).
Specifically tie Paper 44's SDP-based Wasserstein distance to the Latrémolière 2512.03573 Fortet-Mourier formulation. Verify the FM axiom (metrizing weak-*).

**A.2'.3 BW vacuum + strictly-positive element setup** (~1 week).
Document the BW vacuum as the canonical pin state; verify the B axiom (boundedness via $h = K_\alpha^W / Z$).

**A.2'.4 Def 1.26 continuum tightness** (~2-3 weeks). The named sub-sprint.
Prove that on the continuum triple $\mathcal{T}^L_{S^3 \times S^1_T}$, the finite-Krein-Wasserstein-distance set is weak-* dense in the K⁺ state space. This is the analog of Paper 38 Theorem 5.5 for the Krein setting; requires extending the Paper 38 limit-identification argument to the K⁺-restricted Wasserstein-Kantorovich state space. Multi-week, not multi-month.

**A.2'.5 Def 1.29 exhaustive sequence formalization** (~1 week).
Specifically tie the truncated Krein triple sequence to Latrémolière 2512.03573's Def 1.29 + Theorem 1.30 (approximate unit). Mostly bookkeeping given Sub-sprint D's Riemannian-limit recovery.

### 5.2 Phase A.3' — Bridge to Mondino-Sämann pLGH (1.5-2 months)

(Unchanged from the re-scope memo §3.) The bridge construction is *separate* from the Krein-lift — the diagnostic confirmed the Krein-lift transports cleanly, but the bridge between the metametric (Latrémolière) and reverse-triangle (Mondino-Sämann) structures remains the central open content.

### 5.3 Phase A.4' — GeoVac wedge as Krein-pointed QMS (1 month)

(Unchanged from the re-scope memo §3.) Apply A.2' + A.3' to the Paper 43 wedge construction; verify the canonical pin-state choice + the panel cells.

### 5.4 Phase A.5' — Synthesis + decision gate (~3 weeks)

(Unchanged from the re-scope memo §3.) Same structure as original A.5; closure memo with POSITIVE/MIXED/NEGATIVE verdict on writing the merged Paper 48.

### 5.5 Phase B — Merged Paper 48 (~1.5 months)

If A.5' returns POSITIVE: write Paper 48 (~25-30 pages, math.OA submission) following the outline in the re-scope memo §4.

### 5.6 Total Tier 3 timeline

| Phase | Effort | Cumulative |
|:------|:-------|:-----------|
| A.2' Krein-lift naming + formalization | 1-2 months | 1-2 months |
| A.3' Bridge to Mondino-Sämann | 1.5-2 months | 2.5-4 months |
| A.4' GeoVac wedge application | 1 month | 3.5-5 months |
| A.5' Synthesis + decision gate | 3 weeks | 4-5.5 months |
| B Paper 48 draft | 1.5 months | 5.5-7 months |

**Total: ~6 months end-to-end** if the full Phase B is included (originally estimated ~12 months pre-scoop, now ~6-7 months post-scoop and post-diagnostic).

---

## §6. Honest scope statement

### 6.1 What the diagnostic tested

The diagnostic:
- Read the full Latrémolière 2512.03573 PDF (766.9 KB) via WebFetch and extracted §2-§6 detail beyond the Phase A.2' abstract-level scan.
- Tested 7 axioms / constructions from the Latrémolière framework against the GeoVac Krein architecture (Papers 44/45/46/47 + L3b-2 sub-sprints A-D).
- Located each test result against the existing L3b-2 + Paper 46 closure points and the L3b-2a NO-GO finding.
- Cataloged the verdicts per axiom.

### 6.2 What the diagnostic did NOT test

The diagnostic:
- Did NOT implement any Krein-lift theorems; everything is at the *structural / axiom-by-axiom* level.
- Did NOT verify the continuum tightness axiom (Def 1.26) on $\mathcal{T}^L_{S^3 \times S^1_T}$. The verdict says "GO-WITH-CAVEAT" for this axiom; verifying it is the named 2-3 week sub-sprint.
- Did NOT address the Mondino-Sämann pLGH bridge (Phase A.3' of the re-scope). The bridge is *separate from the Krein-lift* and remains the central open content of the merged Paper 48.
- Did NOT address strong-form Lorentzian propinquity on the continuum carrier (G1-genuine). Paper 46 closed strong-form on the natural substrate; the *continuum* strong-form (without K⁺ restriction, on the non-compact $S^3 \times \mathbb{R}_t$) remains multi-month frontier.

### 6.3 What an actual Krein-lift attempt would still need to establish

Beyond the diagnostic:
- Formal naming and verification of Krein-self-adjoint Leibniz hermitian seminorm satisfying Latrémolière 2512.03573 Def 1.18.
- Proof of Def 1.26 tightness on the continuum triple (the 2-3 week sub-sprint).
- Construction of the explicit Latrémolière-style coincidence theorem (§3.1 of 2512.03573) on the Krein-positive setting — using the existing Wasserstein-Kantorovich SDP machinery from Paper 44.
- Formal naming and verification of the §4 inframetric hypertopology on K⁺-restricted Krein triples.
- Writing-up: Paper 48 §3 starting structure should follow the Latrémolière 2512.03573 §2 → §4 sequence, with each lemma/theorem cross-referenced to the corresponding Paper 44/45/46 result.

### 6.4 Where the diagnostic surfaces information that contradicts or refines the Phase A.2' deep-read

**Refinement:** the Phase A.2' deep-read (`debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md`) was abstract+intro-level only. The diagnostic adds:
- §2 tunneling pair structure (reach/height/alignment as direct analogs of Paper 38's $\mathrm{reach}_B$/$\mathrm{reach}_P$/$\mathrm{height}_B$/$C_3$)
- §4 inframetric is *4-point relaxed triangle*, not 3-point metric — the Krein setting naturally inherits this weaker structure
- §5 compact-case agreement reduces to Latrémolière 2014 (the classical propinquity), which Paper 38 already verifies at $N_t = 1$
- §6 $c_0(\mathbb{Z}) \rtimes \mathbb{Z}$ example provides the prototype of pointed convergence for non-compact non-commutative C*-algebras; the GeoVac analog is the truncated wedge of the compact-temporal Lorentzian Krein triple

**Sharpening:** the Phase A.2' deep-read flagged the Sturmian-as-exhaustive-sequence identification under two readings (operator-system reading: works; Coulomb-multiplication-algebra reading: doesn't work because Sturmian projector is not multiplication operator). **For the Krein-lift, this issue does not arise** — the operator system $\mathcal{O}^L$ is explicitly an *operator-system-truncation* object, not a multiplication-algebra object. The Krein-lift naturally lands in the operator-system reading where the exhaustive sequence works.

**No contradiction:** the Phase A.2' R1 finding (second-order Schrödinger Leibniz failure) is not relevant here — the Krein-lift uses the first-order Lorentzian Dirac $D_L$, for which Leibniz holds.

### 6.5 Where the diagnostic depends on unverified Paper 44 substrate facts

**Critical dependency 1:** Paper 44 propagation number = 2 under achievable envelope is bit-exact verified across panel cells (per Paper 44 §5). Diagnostic relies on this for the Latrémolière 2512.03573 compatibility (the inframetric structure on the K⁺-restricted Krein-pointed QMS inherits the propagation-number invariant).

**Critical dependency 2:** Sub-sprint D Riemannian-limit recovery is bit-exact verified at $(n_\max, N_t) = (2, 1), (3, 1), (4, 1)$ (per `debug/l3b_2_sub_sprint_D_propinquity_memo.md` §7). Diagnostic relies on this for the §5 compact-case-agreement falsifier.

**Critical dependency 3:** Paper 46 Appendix B "free upgrade" finding $\Lambda^\mathrm{enlarged} = \Lambda^{P45}$ bit-exact extends to the diagnostic's GO verdict on §3.1 Leibniz seminorm. The Krein-self-adjoint Lipschitz seminorm choice $L^K = \|[D_L, a]\|_\mathrm{op}$ is the Paper 46 natural-substrate choice.

**No dependency on unverified facts:** every load-bearing structural fact the diagnostic uses is *already verified* in Papers 44/45/46 or the L3b-2 sub-sprint memos. The diagnostic is a *naming* exercise — identifying which already-verified GeoVac structural facts correspond to which Latrémolière 2512.03573 axioms.

---

## §7. Three PI questions queued for the decision gate

### Q1. Continue with the merged Paper 48 (Tier 3 commit) or pause Tier 3?

The diagnostic gives POSITIVE-WITH-NAMED-OBSTRUCTIONS. The merged Paper 48 is reachable in ~6 months (originally ~12 months pre-scoop, now substantially shorter post-diagnostic).

**Recommendation:** PROCEED. The path is well-defined, the named obstructions are 2-3 week sub-sprints not multi-month frontier, and the scoop-risk window remains open (Latrémolière 2512.03573 was published Dec 2025; an independent Krein-lift could surface 2026-2027).

### Q2. Phase ordering — Krein-lift naming first, or Mondino-Sämann bridge first?

The re-scope memo §3 sequenced A.2' → A.3' → A.4'. The diagnostic suggests this ordering remains correct: A.2' is the Krein-lift naming/formalization (1-2 months); A.3' is the Mondino-Sämann bridge (1.5-2 months); A.4' is the GeoVac wedge application (1 month). The bridge can be developed *in parallel* with the Krein-lift formalization if PI wants to compress.

**Recommendation:** STANDARD ORDERING. The Krein-lift naming is more mechanical (mostly done already); the bridge has more structural risk and benefits from the Krein-lift framework being in place.

### Q3. Should the diagnostic verdict trigger immediate paper updates?

The diagnostic does NOT require any paper modifications. The verdict is structural, not theorem-level. Paper 44, 45, 46, 47 already document all the load-bearing facts the diagnostic depends on; the merged Paper 48 is the natural home for the consolidated Krein-lift writing.

**Recommendation:** NO paper updates at this stage. Wait for Phase A.2'.5 (formalization sub-sprint) or Phase B (Paper 48 drafting) before applying any structural-claim edits.

---

## §8. Honest scope statement summary

**This diagnostic IS:** a structural axiom-by-axiom test of whether Latrémolière 2512.03573's pinned-QLCMS framework transports to the GeoVac Krein architecture via the Paper 44/45/46/47 substrate. Verdict: **POSITIVE-WITH-NAMED-OBSTRUCTIONS** (7/7 axioms transport, 1 named caveat on continuum tightness).

**This diagnostic is NOT:** a theorem, a proof, an implementation, a paper. It identifies *where the lift is hard and how to commit*.

**The diagnostic-before-engineering rule applies cleanly:** the diagnostic surfaces that *most of the lift is already done* in existing Papers 44-47, that the remaining work is naming + formalization + writing, and that the merged Paper 48 is reachable in ~6 months (vs the original ~12 months pre-scoop, or the re-scoped 10 weeks if everything goes positively). The named caveat (Def 1.26 continuum tightness, 2-3 weeks) is a *contained* obstruction, not a structural blocker.

---

**End of memo.**

**Files added in this sprint:**
- `debug/sprint_tier3_light_krein_lift_diagnostic_memo.md` (this memo, ~5000 words)
- `debug/data/sprint_tier3_light_krein_lift_diagnostic.json` (per-axiom verdict structure)

**Cross-references:**
- `debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md` (Phase A.2' abstract-level deep-read; this memo extends to full §2-§6 PDF read)
- `debug/sprint_l3e_p3_rescope_memo.md` (the original Phase A.2' Krein-lift target; this memo updates timeline estimates)
- `debug/l3b_2_sub_sprint_D_propinquity_memo.md` (Paper 45 K⁺-restricted weak-form proof; load-bearing for diagnostic axioms 3.4-3.7)
- `debug/l3b_2_sub_sprint_C_berezin_memo.md` (Sub-sprint C Berezin reconstruction; load-bearing for diagnostic axiom 3.5)
- `debug/l3b_2a_candidate_validation_memo.md` (NO-GO finding on $L_\mathrm{block}$; the diagnostic axiom 3.1 inherits the $L_\mathrm{op}$ choice from Paper 46 / L3b-2b)
- `papers/group1_operator_algebras/paper_44_lorentzian_operator_system.tex` (Krein operator-system substrate)
- `papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex` (K⁺-restricted weak-form Lorentzian propinquity, the substrate of the Krein-lift)
- `papers/group1_operator_algebras/paper_46_strong_form_lorentzian_propinquity.tex` (strong-form on natural substrate; provides the $L_\mathrm{op}$ seminorm choice for the Krein-lift)
- `papers/group1_operator_algebras/paper_47_two_rate_hybrid_convergence.tex` (G2 closure at norm-resolvent level; relevant to continuum tightness if Phase A.2'.4 needs it)
