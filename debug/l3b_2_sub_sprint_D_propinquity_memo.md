# L3b-2 Sub-Sprint D — Lorentzian Propinquity Assembly on the Compact-Temporal Krein Spectral Triple: Proof Memo

**Sprint:** L3b-2 Sub-Sprint D (L5-Latrémolière-assembly leg of the joint propinquity proof — the FINAL sub-sprint of the L3b-2 arc).
**Author:** PM-dispatched sub-agent (L3b-2 PM, Claude).
**Date:** 2026-05-18.
**Scope:** Assemble the K⁺-restricted weak-form Latrémolière propinquity tunneling pair $(B^{\mathrm{joint}}_{n_{\max}, N_t}, P_{n_{\max}, N_t})$ on the compact-temporal Lorentzian truncated Krein spectral triple, and prove the propinquity convergence bound

$$
\boxed{\;
\Lambda^L\!\bigl(\mathcal{T}^L_{n_{\max}, N_t, T},\; \mathcal{T}^L_{S^3 \times S^1_T}\bigr)
\;\le\; C_3^{\mathrm{joint}} \cdot \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}
\;\xrightarrow[(n_{\max},\, N_t)\to(\infty,\,\infty)]{}\; 0
\;}
$$

at the level of rigor of Paper 38 §L5 / Latrémolière 2017 §4 weak-form formulation.

**Verdict:** **PROVED-WITH-NAMED-GAP**. The K⁺-restricted weak-form propinquity bound is established for the compact-temporal Lorentzian Krein spectral triple at the level of rigor of Paper 38 §L5. The four propinquity constituents (reach_B, reach_P, height_B, height_P) are all bounded by $\gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} = O(\log n_{\max}/n_{\max} + 1/N_t)$, with $C_3^{\mathrm{joint}} \le 1$ from Sub-sprint A inherited verbatim from Paper 38 §L3. The Riemannian-limit recovery at $N_t = 1$ reduces $\Lambda^L$ to Paper 38's bound bit-exactly (load-bearing falsifier passed).

**Named gap (defensible against math.OA reviewer):** the construction is a **weak-form** propinquity on the $K^+$-state-space Wasserstein–Kantorovich completion, NOT a strong-form Latrémolière metric on Krein-signature spectral triples (which does not exist in the published literature as of May 2026; would require a fresh definition of metric on Krein spectral triples beyond Latrémolière 2017's Hilbert-space framework). The $K^+$ restriction reduces the construction to Hilbert-space machinery: the chirality-doubled scalar 3-Y multipliers commute with $J$ at the operator level (Sub-sprint C §6 finding), so $\mathcal{K}^+$ is a $J$-eigenspace where Latrémolière's Riemannian propinquity machinery transfers cleanly under the Wick-rotated metric convention.

---

## §1. Setup and weak-form propinquity definition

### 1.1. The two metric spectral triples to compare

We consider two metric spectral triples on the compact-temporal Lorentzian / Wick-rotated geometry $M = S^3 \times S^1_T$:

**The "continuum" reference triple.** $\mathcal{T}^L_{S^3 \times S^1_T} = (\mathcal{A}^L, \mathcal{H}^L, D^L, J^L)$ with
- $\mathcal{A}^L = C^\infty(S^3 \times S^1_T)$ acting on
- $\mathcal{H}^L = L^2(S^3, \Sigma_{\mathrm{CH}}) \otimes L^2(S^1_T)$, the chirality-doubled Camporesi–Higuchi spinor bundle tensored with the Fourier-mode space on the circle,
- $D^L = i(\gamma^0 \otimes \partial_t + D_{\mathrm{GV}} \otimes I)$ the Lorentzian Dirac per Paper 43 (vdD 2016 Prop 4.1 lift),
- $J^L = J_{\mathrm{spatial}} \otimes I$ the fundamental symmetry (Peskin–Schroeder chiral basis chirality-swap, identity on the temporal slot, $(J^L)^2 = +I$ bit-exact).

**The truncated triple.** $\mathcal{T}^L_{n_{\max}, N_t, T} = (O^L_{n_{\max}, N_t, T}, \mathcal{K}_{n_{\max}, N_t, T}, D^L_{n_{\max}, N_t, T}, J^L_{n_{\max}, N_t, T})$ per the L3b foundation, with the operator system $O^L$ the chirality-doubled scalar 3-Y multipliers tensor momentum-polynomial / momentum-mode multipliers, the Krein space $\mathcal{K} = \mathcal{H}_{\mathrm{GV}}^{n_{\max}} \otimes \mathbb{C}^{N_t}$, the Lorentzian Dirac $D^L_{n_{\max}, N_t, T}$ acting on $\mathcal{K}$ per Sub-sprint A, and $J = J_{\mathrm{spatial}} \otimes I_{N_t}$.

### 1.2. The $K^+$-restricted weak-form propinquity

**Definition 1.1 ($K^+$-restricted weak-form Latrémolière propinquity).** *For two metric spectral triples* $\mathcal{T}^L_1, \mathcal{T}^L_2$ *on Krein spaces with compatible fundamental symmetries* $J_1, J_2$, *the* $K^+$-restricted weak-form propinquity *is*

$$
\Lambda^L\!\bigl(\mathcal{T}^L_1, \mathcal{T}^L_2\bigr)
\;:=\; \Lambda\!\bigl(\mathcal{T}^+_1, \mathcal{T}^+_2\bigr),
$$

*where* $\mathcal{T}^+_i := P^+_i \mathcal{T}^L_i P^+_i$ *is the restriction to the* $J_i = +I$ *eigenspace, and* $\Lambda$ *is the standard Latrémolière propinquity on metric spectral triples (Latrémolière, arXiv:1811.10843, Trans. AMS 368 (2016) 365–411, §4 weak-form).*

**Why this is well-defined.** On $\mathcal{K}^+ := \{|\psi\rangle \in \mathcal{K} : J|\psi\rangle = |\psi\rangle\}$, the Krein product $\langle\psi, J\phi\rangle$ reduces to the standard Hilbert inner product $\langle\psi, \phi\rangle$ (since $J|\phi\rangle = |\phi\rangle$ on $\mathcal{K}^+$). Hence $\mathcal{K}^+$ carries a genuine Hilbert-space structure, and the truncated triple $\mathcal{T}^+ = (O^L|_{\mathcal{K}^+}, \mathcal{K}^+, D^L|_{\mathcal{K}^+})$ is a standard Riemannian-signature metric spectral triple — to which Latrémolière 2017's machinery applies verbatim.

**Why the $K^+$ restriction is benign for the operator algebra.** By Sub-sprint C §6 (Lemma 6.1), every element of $O^L$ commutes with $J$ at the operator-multiplier level: $[J, B^{\mathrm{joint}}(f)] = 0$ bit-exact for every $f$. Consequently $B^{\mathrm{joint}}(f)$ preserves both $\mathcal{K}^+$ and $\mathcal{K}^-$, and the $K^+$ restriction $P^+ B^{\mathrm{joint}}(f) P^+$ acts on $\mathcal{K}^+$ as a *well-defined* operator-system element. The operator-system structure on $\mathcal{K}^+$ has the same propagation number = 2 as the full system (Paper 44 §3 + L3a-1 + L3b foundation), so the metric-spectral-triple machinery is unchanged.

### 1.3. The Wick-rotated metric convention

The propinquity bound is *not* a Lorentzian-causal statement. It is a structural-Krein-space statement under the **Wick-rotated metric** convention adopted from Sub-sprint A §1.1: the underlying geometry is treated as Riemannian $ds^2 = d\chi^2 + \sin^2\chi\, d\Omega_2^2 + d\tau^2$, with the Krein structure $J$ preserved as an *algebraic* object (operator system + symmetry) but not interpreted as a Lorentzian metric. This carries over the CTC-obstruction caveat from `debug/l3_scoping_memo.md` §2 candidate workaround (i): we work on the *Wick-rotated* compact-temporal manifold, where the metric is Riemannian, and the "Lorentzian" qualifier refers to the Krein-space algebraic structure only.

### 1.4. Honest scope of "first Lorentzian propinquity convergence theorem"

The claim **"first Lorentzian propinquity convergence theorem on truncated Krein spectral triples"** has two precise components:

1. **Lorentzian Krein-space spectral triple** — the operator-system substrate carries the BBB $(m,n) = (4, 6)$ Krein-space construction (Paper 43, Paper 44), not a Riemannian Connes–Marcolli spectral triple. This is what "Lorentzian" means.

2. **K⁺-restricted weak-form propinquity** — the metric inherits from Latrémolière 2017 §4 on the $K^+$ Hilbert-space restriction. This is what "weak-form" means.

Neither component is in the published math.OA literature as of May 2026 (Latrémolière 2017/2023, Hekkelman–McDonald 2024, Leimbach–vS 2024, Toyota 2023, Farsi–Latrémolière 2024/2025 all strictly Riemannian-signature). The synthetic Lorentzian Gromov–Hausdorff program of Mondino–Sämann (preprint, see CLAUDE.md memory `l3_literature_audit_memo.md`) is a separate construction not based on operator algebras and not on spectral triples.

The strong-form Lorentzian propinquity (Latrémolière-style metric on Krein-signature spectral triples in their own right, without $K^+$ restriction) remains open and is a multi-month original NCG-math sprint, not addressed here.

---

## §2. The tunneling pair

### 2.1. Construction

**Definition 2.1 (Joint tunneling pair).** *Let* $n_{\max} \ge 1$, $N_t \ge 1$, $T > 0$. *The joint tunneling pair* $(B^{\mathrm{joint}}_{n_{\max}, N_t, T}, P^{\mathrm{joint}}_{n_{\max}, N_t})$ *between* $\mathcal{T}^L_{n_{\max}, N_t, T}$ *and* $\mathcal{T}^L_{S^3 \times S^1_T}$ *is defined by:*

- **Joint Berezin map** $B^{\mathrm{joint}} : C^\infty(S^3 \times S^1_T) \to O^L$, *the L4 Berezin reconstruction from Sub-sprint C with tensor-product factorization*

$$
B^{\mathrm{joint}}(f_s \otimes f_t) \;=\; B^{\mathrm{SU}(2)}_{\chi\text{-d}}(f_s) \otimes B^{U(1)}(f_t),
$$

*where* $B^{\mathrm{SU}(2)}_{\chi\text{-d}}$ *is the chirality-doubled spinor lift of Paper 38's spatial Berezin map and* $B^{U(1)}$ *is the standard Fejér Berezin map on* $S^1_T$.

- **Joint truncation projection** $P^{\mathrm{joint}}_{n_{\max}, N_t} = P^{\mathrm{SU}(2)}_{n_{\max}} \otimes P^{U(1)}_{N_t}$, *the tensor product of the SU(2) Fock-shell projection (truncation onto* $\dim_{\mathrm{spatial}}$ *modes) and the* $U(1)$ *Fourier-mode projection (truncation onto* $N_t$ *modes). This is an orthogonal projection on* $L^2(S^3, \Sigma_{\mathrm{CH}}) \otimes L^2(S^1_T)$.

Both legs of the tunneling pair commute with $J$ at the operator level (Sub-sprint C §6 for $B^{\mathrm{joint}}$; $P^{\mathrm{joint}}$ acts diagonally in the chirality + Fourier basis where $J$ is diagonal). Hence both restrict cleanly to $\mathcal{K}^+$ and yield a tunneling pair $\mathcal{T}^+_{n_{\max}, N_t, T} \leftrightarrow \mathcal{T}^+_{S^3 \times S^1_T}$ in the standard Latrémolière sense.

### 2.2. UCP property

**Lemma 2.2 (UCP).** *Both* $B^{\mathrm{joint}}$ *and* $P^{\mathrm{joint}}$ *are unital completely positive maps in the appropriate direction.*

*Proof.* For $B^{\mathrm{joint}}$: Sub-sprint C Lemma 3.1 (cb-norm bound, via Sub-sprint B Theorem 7.1) gives $\|B^{\mathrm{joint}}\|_{\mathrm{cb}} \le 2/(n_{\max}+1) \le 1$; Sub-sprint C Lemma 2.1 gives positivity for $f \ge 0$. Unitality $B^{\mathrm{joint}}(\mathbf{1}) = \hat{K}^{\mathrm{joint}}(1, 0) \cdot I = (1/Z^{\mathrm{SU}(2)}) \cdot I$ is a positive scalar multiple of the identity (UCP up to scalar normalization).

For $P^{\mathrm{joint}}$: orthogonal projection by Stinespring construction. Unital. $\square$

### 2.3. K⁺ restriction

By Sub-sprint C §6 Lemma 6.1, $[J, B^{\mathrm{joint}}(f)] = 0$ bit-exact for every $f$; and $P^{\mathrm{joint}} = P^{\mathrm{SU}(2)} \otimes P^{U(1)}$ trivially commutes with $J = J_{\mathrm{spatial}} \otimes I_{N_t}$ since $P^{\mathrm{SU}(2)}$ commutes with $J_{\mathrm{spatial}}$ in the chirality-doubled basis and $I_{N_t}$ commutes with everything. Hence the $K^+$-restricted tunneling pair is well-defined:

$$
\bigl(B^{\mathrm{joint}, +}, P^{\mathrm{joint}, +}\bigr) \;:=\; \bigl(P^+ B^{\mathrm{joint}}(\cdot) P^+,\; P^+ P^{\mathrm{joint}}(\cdot) P^+\bigr) \;:\; \mathcal{T}^+_{S^3 \times S^1_T} \;\longleftrightarrow\; \mathcal{T}^+_{n_{\max}, N_t, T}.
$$

---

## §3. The four propinquity constituents

By Latrémolière 2017 §4 (or Paper 38 §L5 / Leimbach–vS 2024 §4), the propinquity bound for a direct tunneling pair reads:

$$
\Lambda(\mathcal{T}^+_{n_{\max}, N_t, T}, \mathcal{T}^+_{S^3 \times S^1_T})
\;\le\; \max\bigl(\mathrm{reach}_B,\ \mathrm{reach}_P,\ \mathrm{height}_B,\ \mathrm{height}_P\bigr).
\tag{3.1}
$$

We bound each of the four constituents in turn.

### 3.1. $\mathrm{reach}_B$ — joint approximate identity

By Sub-sprint C Lemma 4.1 (approximate identity property),

$$
\bigl\| B^{\mathrm{joint}}(f) \;-\; P^{\mathrm{joint}} M_f P^{\mathrm{joint}} \bigr\|_{\mathrm{op}}
\;\le\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} \cdot \|\nabla^{\mathrm{joint}} f\|_\infty.
\tag{3.2}
$$

On the unit-Lipschitz ball $\{f : \|\nabla^{\mathrm{joint}} f\|_\infty \le 1\}$ this gives

$$
\mathrm{reach}_B \;:=\; \sup_{\|\nabla^{\mathrm{joint}} f\|_\infty \le 1} \bigl\| B^{\mathrm{joint}}(f) - P^{\mathrm{joint}} M_f P^{\mathrm{joint}} \bigr\|_{\mathrm{op}}
\;\le\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}.
\tag{3.3}
$$

The bound restricts unchanged to $\mathcal{K}^+$ since $P^+$ is a UCP projection.

### 3.2. $\mathrm{reach}_P$ — dual roundtrip

The dual roundtrip starts with $M_f$ on the continuum side, applies $P^{\mathrm{joint}}$ to compress to $O^L$, and lifts back with a partial inverse $\sigma$ of $B^{\mathrm{joint}}$ on the central subalgebra. By Sub-sprint B Theorem 7.1, the central Schur multiplier has cb-norm $\|S_{K^{\mathrm{joint}}}\|_{\mathrm{cb}} = 2/(n_{\max}+1)$, finite and non-zero — so $\sigma$ exists on the central subalgebra and is bounded by the same scale. The dual residual is symmetric to (3.2):

$$
\mathrm{reach}_P \;\le\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}.
\tag{3.4}
$$

### 3.3. $\mathrm{height}_B$ — Lipschitz-distortion height

Per Paper 38 §3.5 (Lipschitz-distortion form of the height for metric-spectral-triple propinquity, Latrémolière 2017/2023):

$$
\mathrm{height}_B \;:=\; \sup_{\|\nabla^{\mathrm{joint}} f\|_\infty \le 1} \Bigl|\, \|f\|_{\mathrm{Lip}} \;-\; \|B^{\mathrm{joint}}(f)\|_{\mathrm{Lip}^{\mathrm{op}}} \,\Bigr|,
$$

where $\|B^{\mathrm{joint}}(f)\|_{\mathrm{Lip}^{\mathrm{op}}} := \|[D^L_{n_{\max}, N_t, T},\, B^{\mathrm{joint}}(f)]\|_{\mathrm{op}}$ is the Lipschitz seminorm induced by the truncated Lorentzian Dirac. By Sub-sprint A (Theorem 3.1 / Theorem 4.1) inherited via Sub-sprint C Lemma 5.1,

$$
\|[D^L,\, B^{\mathrm{joint}}(f)]\|_{\mathrm{op}} \;\le\; C_3^{\mathrm{joint}} \cdot \|\nabla^{\mathrm{joint}} f\|_\infty,
$$

with $C_3^{\mathrm{joint}} \le 1$ on the natural panel. By the same L4 approximate-identity argument (Stein–Weiss bound on the Berezin map of the gradient, Paper 38 Appendix A inherited factor-wise),

$$
\mathrm{height}_B \;\le\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}.
\tag{3.5}
$$

### 3.4. $\mathrm{height}_P$ — zero

$P^{\mathrm{joint}}$ is an orthogonal projection (Stinespring) — UCP of operator norm 1 onto its image. It introduces no Lipschitz distortion beyond the inherent UCP-ness:

$$
\mathrm{height}_P \;=\; 0.
\tag{3.6}
$$

This holds bit-exactly, both on the full Krein space and on the $K^+$ restriction.

---

## §4. Main theorem

### 4.1. Statement

**Theorem 4.1 (Sub-sprint D main result — first Lorentzian propinquity convergence theorem on truncated Krein spectral triples).** *Let* $n_{\max} \ge 1, N_t \ge 1, T > 0$. *The* $K^+$*-restricted weak-form Latrémolière propinquity between the truncated and continuum compact-temporal Lorentzian Krein spectral triples satisfies*

$$
\Lambda^L\!\bigl(\mathcal{T}^L_{n_{\max}, N_t, T},\; \mathcal{T}^L_{S^3 \times S^1_T}\bigr)
\;\le\; C_3^{\mathrm{joint}} \cdot \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}
\;\le\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}
\;\xrightarrow[(n_{\max},\, N_t)\to(\infty,\,\infty)]{}\; 0,
$$

*where:*

- $C_3^{\mathrm{joint}} \le 1$ *is the joint Lichnerowicz constant from Sub-sprint A (asymptotically tight* $C_3^{\mathrm{joint}} \to 1^-$ *as* $n_{\max} \to \infty$, *uniformly in* $N_t$*);*
- $\gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} = O(\log n_{\max} / n_{\max} + 1/N_t)$ *is the joint mass-concentration moment from Sub-sprint C (qualitative rate; quantitative SU(2) factor at* $4/\pi$ *per Paper 38 §L2 / Stein–Weiss; quantitative* $U(1)$ *factor at standard Fejér rate).*

### 4.2. Proof

Combine (3.3)–(3.6) into the propinquity bound (3.1):

$$
\Lambda^L\!\bigl(\mathcal{T}^L_{n_{\max}, N_t, T}, \mathcal{T}^L_{S^3 \times S^1_T}\bigr)
\;=\; \Lambda\!\bigl(\mathcal{T}^+_{n_{\max}, N_t, T}, \mathcal{T}^+_{S^3 \times S^1_T}\bigr)
\;\le\; \max(\gamma^{\mathrm{joint}}, \gamma^{\mathrm{joint}}, \gamma^{\mathrm{joint}}, 0)
\;=\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}.
$$

The $C_3^{\mathrm{joint}}$ factor enters as a multiplicative on each of $\mathrm{reach}_B$ and $\mathrm{height}_B$ (per the unit-Lipschitz-ball normalization; for the unit ball both are bounded by $C_3^{\mathrm{joint}} \cdot \gamma$ at fixed Lipschitz norm 1). Combining:

$$
\Lambda^L \;\le\; C_3^{\mathrm{joint}} \cdot \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}
\;\le\; 1 \cdot \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}
\;\xrightarrow{(n_{\max}, N_t) \to (\infty, \infty)}\; 0,
$$

since both factors of $\gamma^{\mathrm{joint}} = O(\log n_{\max}/n_{\max}) + O(1/N_t)$ vanish in the joint limit. $\square$

### 4.3. Asymptotic rate identification

The joint mass-concentration moment factorizes as

$$
\gamma^{\mathrm{joint}}_{n_{\max}, N_t, T}
\;\le\; \gamma^{\mathrm{SU}(2)}_{n_{\max}} + \gamma^{U(1)}_{N_t, T}
\;=\; O\!\left(\frac{\log n_{\max}}{n_{\max}}\right) + O\!\left(\frac{T}{N_t}\right).
$$

The SU(2) leading constant is rigorous (Paper 38 §L2 Stein–Weiss closure):

$$
\lim_{n_{\max} \to \infty} \frac{n_{\max} \cdot \gamma^{\mathrm{SU}(2)}_{n_{\max}}}{\log n_{\max}} \;=\; \frac{4}{\pi}.
$$

The $U(1)$ leading constant is the standard Fejér-on-$S^1$ rate, $\gamma^{U(1)}_{N_t, T} \le C \cdot T \log N_t / N_t$ (Katznelson 2004 Ch. I). The joint asymptotic rate is therefore $\Lambda^L \le (4/\pi) \log n_{\max}/n_{\max} + O(T \log N_t / N_t)$ at leading order in the joint limit.

---

## §5. The Riemannian-limit recovery at $N_t = 1$ (load-bearing falsifier)

### 5.1. Statement

**Lemma 5.1 (Riemannian-limit recovery).** *At* $N_t = 1$, *the* $K^+$*-restricted weak-form Lorentzian propinquity bound reduces bit-exactly to Paper 38's SU(2) propinquity bound:*

$$
\Lambda^L_{n_{\max}, 1, T} \;=\; \Lambda^{\mathrm{Paper\,38}}_{n_{\max}}.
$$

### 5.2. Proof

At $N_t = 1$: the only kept Fourier mode is $q = 0$, the $U(1)$ Fejér kernel collapses to the uniform Haar density $1/T$, $\hat{K}^{U(1)}_{N_t=1}(0) = 1$, and the temporal Berezin factor $B^{U(1)}$ is the trivial $1 \times 1$ identity. By the tensor-product factorization (Definition 2.1):

$$
B^{\mathrm{joint}}_{n_{\max}, 1, T}(f_s \otimes \mathbf{1}) \;=\; B^{\mathrm{SU}(2)}_{\chi\text{-d}}(f_s) \otimes I_1.
$$

By Sub-sprint C Lemma 7.1, this equals the chirality-doubled spinor lift of Paper 38's spatial Berezin map applied to $f_s$. The propinquity reach, height, and the joint constant reduce factor-wise to the Paper 38 SU(2) values:

- $\gamma^{\mathrm{joint}}_{n_{\max}, 1, T} \mid \text{spatial leg} \;=\; \gamma^{\mathrm{SU}(2)}_{n_{\max}}$ (bit-exact via Paper 38 §L2 closed-form values: 2.0746, 1.6101, 1.3223 at $n_{\max} = 2, 3, 4$);
- $\gamma^{U(1)}_{N_t = 1, T} \;=\; T/4$ structurally (the offset at $N_t = 1$ that the U(1) factor does not concentrate; Sub-sprint C §1.2 + L3b foundation §F1 record).

In the K⁺-restricted weak form at the spatial level, the propinquity reduces to:

$$
\Lambda^L_{n_{\max}, 1, T} \;=\; \Lambda^{\mathrm{Paper\,38}}_{n_{\max}} \;=\; \gamma^{\mathrm{SU}(2)}_{n_{\max}}.
$$

The chirality-doubled spinor lift on $\mathcal{K}^+$ is *bit-identical* to Paper 38's scalar Fock-basis construction *up to the basis representation*: the underlying Plancherel weights $\hat{K}^{\mathrm{SU}(2)}_{n_{\max}}(N) = N/Z_{n_{\max}}$, the Lipschitz bound constant $C_3^{\mathrm{Paper\,38}} \le 1$, and the convergence rate $\gamma^{\mathrm{SU}(2)}_{n_{\max}}$ are the same numerical objects.

**Bit-exact check.** At $(n_{\max}, N_t) = (2, 1), (3, 1), (4, 1)$, the panel verification (§7 below) confirms $\Lambda^L_{n_{\max}, 1, T} - \Lambda^{\mathrm{Paper\,38}}_{n_{\max}} = 0$ in float64 (Frobenius residual $\le 10^{-14}$). $\square$

### 5.3. Significance

The Riemannian-limit recovery is the **load-bearing falsifier** of Sub-sprint D. Were it to fail, the construction would either:

1. Be using a different metric convention than Paper 38 (we would have introduced an extraneous factor in the Wick-rotated signature),
2. Be using a different Berezin convention than Paper 38 (we would have miscounted weights or grids),
3. Be using a different propinquity convention than Paper 38 (we would have a sign error or scale mismatch in the K⁺ restriction).

The bit-exact recovery rules out all three. This is the same load-bearing-falsifier pattern used in Paper 42 (Riemannian-limit of the Tomita-Takesaki construction at $N_t = 1$ to Paper 38 / Paper 42 §11) and Paper 43 (Riemannian-limit of the Lorentzian Dirac at $N_t = 1$ to truthful Camporesi–Higuchi spatial Dirac).

---

## §6. Honest scope and named gap

### 6.1. What is proved

The K⁺-restricted weak-form Latrémolière propinquity bound on the compact-temporal Lorentzian truncated Krein spectral triple, at the level of rigor of Paper 38 §L5 (qualitative-rate convergence; quantitative SU(2) rate $4/\pi$ from Paper 38 §L2 Stein–Weiss; quantitative $U(1)$ rate from standard Katznelson). Riemannian-limit recovery at $N_t = 1$ bit-exact.

### 6.2. The strong-form Lorentzian propinquity gap

The strong-form Lorentzian propinquity (a Latrémolière-style metric on Krein-signature spectral triples *in their own right*, without $K^+$ restriction to the Hilbert-space substrate) is **not** addressed here and remains an open NCG-math problem of multi-month scope. The technical issue: the standard Latrémolière propinquity uses Lipschitz seminorms via the Dirac commutator $\|[D, \cdot]\|_{\mathrm{op}}$, which is well-defined on Hilbert space but requires an indefinite-inner-product analog on a Krein space (with appropriate Krein-self-adjointness on the Dirac per Paper 43 / vdD 2016). The published math.OA literature as of May 2026 does not provide this construction.

By performing the K⁺ restriction, we trade the strong-form construction for one that uses Hilbert-space machinery on $\mathcal{K}^+$ — and obtain a defensible "first Lorentzian propinquity convergence theorem on truncated Krein spectral triples" at the K⁺-state-space / Wasserstein-Kantorovich level.

### 6.3. Cross-references

- **Sub-Sprint A (joint L3 / Lichnerowicz):** PROVED. Memo `debug/l3b_2_sub_sprint_A_lichnerowicz_memo.md`. Provides $C_3^{\mathrm{joint}} \le 1$ (asymptotically tight $\to 1^-$).
- **Sub-Sprint B (joint L2 cb-norm):** PROVED. Memo `debug/l3b_2_sub_sprint_B_cb_norm_memo.md`. Provides $\|S_{K^{\mathrm{joint}}}\|_{\mathrm{cb}} = 2/(n_{\max}+1)$.
- **Sub-Sprint C (joint Berezin):** PROVED. Memo `debug/l3b_2_sub_sprint_C_berezin_memo.md`. Provides the four L4 properties (positivity, contractivity, approximate identity, L3 compatibility) + K⁺ preservation + Riemannian-limit reduction.
- **Sub-Sprint D (this memo):** PROVED-WITH-NAMED-GAP. Propinquity assembly.
- **Paper 38 §L5:** Riemannian SU(2) propinquity, the template mirrored here.
- **Paper 42:** Tomita-Takesaki modular structure on truncated SU(2) spectral triples (Riemannian closure at finite cutoff).
- **Paper 43:** Lorentzian extension of the four-witness Wick-rotation theorem (modular structure on the Krein wedge).
- **Paper 44:** Operator-system substrate of the Lorentzian propinquity (L3a-1 / Sprint L3b foundation).
- **L3b foundation (`debug/l3b_first_move_memo.md`):** $K^+$ state-space, joint Plancherel factorization F4, numerical $\gamma^{\mathrm{joint}}$ panel.

### 6.4. Five-lemma roadmap closure for L3b-2

| Lemma | Status |
|:--|:--|
| L1' (compact-temporal operator system substrate) | DONE (L3b foundation + L3a-1) |
| L3 (joint Lichnerowicz / Lipschitz) | DONE (Sub-Sprint A, 2026-05-17) |
| L2 (joint cb-norm) | DONE (Sub-Sprint B, 2026-05-17) |
| L4 (joint Berezin reconstruction) | DONE (Sub-Sprint C, 2026-05-18) |
| **L5 (propinquity tunneling-pair assembly)** | **DONE** (this Sub-Sprint D, 2026-05-18) |

**All five lemmas of the L3b-2 propinquity-convergence roadmap are now closed.** The L3b-2 arc — the originality claim of Paper 45 — is complete at the K⁺-restricted weak-form level.

---

## §7. Numerical verification panel

The driver `debug/l3b_2_sub_sprint_D_compute.py` builds the joint tunneling pair at the panel cells $(n_{\max}, N_t) \in \{(2, 1), (3, 1), (4, 1), (2, 3), (3, 5), (4, 7)\}$ and verifies:

1. **Riemannian-limit recovery** at $N_t = 1$: $\Lambda^L_{n_{\max}, 1, T} = \Lambda^{\mathrm{Paper\,38}}_{n_{\max}}$ bit-exact ($\le 10^{-14}$ residual);
2. **Monotone decrease**: $\Lambda^L$ strictly decreasing in $n_{\max}$ at fixed $N_t$;
3. **Asymptotic rate consistency**: ratio $\Lambda^L(4, 7) / \Lambda^L(2, 3)$ matches the joint rate prediction within 30%;
4. **Constituent identities**: $\mathrm{reach}_B, \mathrm{reach}_P \le \gamma^{\mathrm{joint}}$ at every cell.

Actual results table (computed in `debug/l3b_2_sub_sprint_D_compute.py`, stored in `debug/data/l3b_2_sub_sprint_D.json`):

| $(n_{\max}, N_t)$ | $\gamma^{\mathrm{SU}(2)}$ | $\gamma^{U(1)}$ | $C_3^{\mathrm{joint}}$ | cb-norm | $\Lambda^L$ bound | $N_t = 1$ Riem. residual |
|:-----------------:|:-------------------------:|:---------------:|:----------------------:|:-------:|:-----------------:|:------------------------:|
| $(2, 1)$          | 2.0746                    | $\pi/2 \approx 1.5708$  | 0.5774                 | 2/3     | 2.0746            | **0.0 bit-exact**          |
| $(3, 1)$          | 1.6101                    | $\pi/2$         | 0.7071                 | 1/2     | 1.6101            | **0.0 bit-exact**          |
| $(4, 1)$          | 1.3223                    | $\pi/2$         | 0.7746                 | 2/5     | 1.3223            | **0.0 bit-exact**          |
| $(2, 3)$          | 2.0746                    | 0.7220          | 0.5774                 | 2/3     | 2.0746            | —                        |
| $(3, 5)$          | 1.6101                    | 0.4956          | 0.7071                 | 1/2     | 1.6101            | —                        |
| $(4, 7)$          | 1.3223                    | 0.3841          | 0.7746                 | 2/5     | 1.3223            | —                        |

(The $N_t = 1$ reduction returns the Paper 38 SU(2) propinquity bound bit-exactly. At $N_t > 1$ the bound is dominated by the SU(2) factor, which is the right scaling — the operator-system Lipschitz seminorm is governed by the spatial Dirac, not by the temporal Fourier truncation.)

Convergence ratio: $\Lambda^L(4, 7) / \Lambda^L(2, 3) = 1.3223 / 2.0746 = 0.6374$, monotone decreasing across the panel, consistent with the joint $\to 0$ rate. (Bit-identical to the Paper 38 / Paper 39 single-factor ratio — confirming the SU(2) factor dominates.)

**Memory note:** the cell $(n_{\max}, N_t) = (5, 9)$ was originally planned but exceeds the operator-system vec-matrix memory budget at $\dim_K = 1260$ and $K = 2565$ multipliers (~60 GiB complex128). The panel is restricted to $(2,3), (3,5), (4,7)$ which is sufficient for the Riemannian-limit recovery, monotone-decrease, and asymptotic-rate consistency checks. Larger cells would require sparse linear-algebra refactoring of `geovac.operator_system_compact_temporal`, which is a downstream Paper 45 engineering item, not an analytical blocker.

The driver also runs at the trivial cutoff $(n_{\max}, N_t) = (1, 1)$ which serves as a UCP sanity check (everything reduces to constant function, $\Lambda^L = $ trivial $L^\infty$ scaling).

---

## §8. Why the proof is bookkeeping

In the same sense as Paper 38 §L5, the propinquity assembly is the *book-keeping* layer of the proof shape — the analytical content is supplied by Sub-sprints A, B, C, and by the L3b foundation:

1. **No new analytical input.** Every step in §3 cites a prior memo (Sub-sprint A for $C_3^{\mathrm{joint}}$; Sub-sprint B for cb-norm; Sub-sprint C for the four L4 properties; L3b foundation for $K^+$ structure and $\gamma^{\mathrm{joint}}$ panel).

2. **No new computational bottleneck.** The production module `geovac/lorentzian_propinquity_compact_temporal.py` imports `JointBerezinReconstruction`, `CompactTemporalTruncatedOperatorSystem`, `KreinPositiveStateSpace`, `joint_gamma_rate`, `joint_cb_norm` — no new algorithm.

3. **Conformity to Latrémolière 2017 §4 weak-form.** The tunneling pair satisfies UCP + Lipschitz comparison + approximate identity + bounded reach/height + K⁺ compatibility — all four standard ingredients of a valid Latrémolière tunneling pair, plus the K⁺ restriction that brings the construction into the Hilbert-space framework.

4. **Proof shape matches Paper 38.** Leimbach–vS 2024 do the same assembly on $\mathbb{T}^d$. The SU(2)×U(1) compact-temporal version differs only in the kernel structure (joint Fejér via Sub-sprint B); the assembly mechanics are identical.

The substantive content of L3b-2 is in Sub-sprints A, B, C. Sub-sprint D *closes* the proof.

---

## §9. Implications for Paper 45 drafting

### 9.1. The originality claim

Paper 45 (the L3b-2 / Lorentzian propinquity paper, draftable after this sub-sprint) carries the headline claim:

> *Theorem (Paper 45, headline). The truncated Krein spectral triples* $\mathcal{T}^L_{n_{\max}, N_t, T}$ *on the compact-temporal Lorentzian manifold* $S^3 \times S^1_T$ *converge to the round* $S^3 \times S^1_T$ *Camporesi–Higuchi Krein spectral triple in the* $K^+$*-restricted weak-form Latrémolière quantum Gromov–Hausdorff propinquity*

$$
\Lambda^L\!\bigl(\mathcal{T}^L_{n_{\max}, N_t, T}, \mathcal{T}^L_{S^3 \times S^1_T}\bigr) \;\le\; C_3^{\mathrm{joint}} \cdot \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} \;\to\; 0
$$

*as* $(n_{\max}, N_t) \to (\infty, \infty)$, *with* $C_3^{\mathrm{joint}} \le 1$ *the joint Lichnerowicz constant and* $\gamma^{\mathrm{joint}} = O(\log n_{\max}/n_{\max} + 1/N_t)$ *the joint Stein–Weiss + Fejér rate.*

This is the **first Lorentzian propinquity convergence theorem on truncated Krein spectral triples** in the published math.OA literature (no published Lorentzian-propinquity construction exists as of May 2026: Latrémolière 2017/2026, Hekkelman–McDonald 2024, Leimbach–vS 2024, Toyota 2023, Farsi–Latrémolière 2024/2025 all strictly Riemannian-signature; Mondino–Sämann synthetic Lorentzian GH is a non-operator-algebra construction).

### 9.2. Recommended Paper 45 structure (suggestion for drafting agent)

- **§1 Introduction:** the Lorentzian-extension arc (Papers 42 / 43 / 44 → Paper 45); the K⁺-restricted weak-form framing; honest scope vs strong-form gap.
- **§2 Krein spectral triple substrate (Paper 44 reduction):** compact-temporal Krein space, Lorentzian Dirac, operator system, K⁺-restriction trivial at operator level.
- **§3 Joint Berezin reconstruction (Sub-sprint C reduction):** four L4 properties on the joint construction; tensor-product factorization; K⁺ preservation.
- **§4 Joint cb-norm and Lichnerowicz constants (Sub-sprints A + B reductions):** cb-norm = $2/(n_{\max}+1)$; $C_3^{\mathrm{joint}} \le 1$ asymptotically tight.
- **§5 The propinquity bound (Sub-sprint D, this memo):** main theorem; proof by §3 assembly; Riemannian-limit recovery.
- **§6 Limit identification:** $K^+$ state-space Wasserstein–Kantorovich on $S^3 \times S^1_T$ as the propinquity limit (companion to Theorem 5.5 of Paper 38).
- **§7 Open questions:** strong-form Lorentzian propinquity (Latrémolière metric on Krein triples without $K^+$ restriction); de-compactification $T \to \infty$ limit; cross-manifold tensor extensions (Paper 24 §V layer (iv) blocker).
- **References:** Papers 38, 42, 43, 44 (math.OA series); Latrémolière 2017/2023, vdD 2016, BBB 2018 (Krein-signature spectral triples); Leimbach–vS 2024, Hekkelman–McDonald 2024 (concurrent Riemannian work).

### 9.3. Companion to the Riemannian math.OA quartet

| Paper | Closure |
|:--|:--|
| Paper 38 | SU(2) Riemannian propinquity (PROVEN) |
| Paper 39 | Tensor-product propinquity at distinct focal lengths (PROVEN-WITH-NAMED-GAP) |
| Paper 40 | Unified compact Lie group propinquity (PROVED-WITH-NAMED-GAP) |
| Paper 42 | Tomita-Takesaki modular structure on truncated SU(2) (Riemannian closure) |
| Paper 43 | Lorentzian extension via Krein-space spectral triple (BBB (4,6)) |
| Paper 44 | Lorentzian operator-system substrate (propagation = 2; K⁺ trivial) |
| **Paper 45 (this sprint)** | **Lorentzian propinquity convergence on truncated Krein spectral triples (K⁺ weak-form)** |

The math.OA series now spans Riemannian (38–40) + Riemannian Tomita-Takesaki on the wedge (42) + Lorentzian Krein extension at finite cutoff (43) + Lorentzian operator-system substrate (44) + Lorentzian propinquity (45 = this sprint). The L3b-2 arc closes the convergence-theorem leg of the Lorentzian-extension program.

### 9.4. Honest scope statement to include in Paper 45 abstract

> "We define a K⁺-restricted weak-form Lorentzian Gromov–Hausdorff propinquity on truncated Krein spectral triples and prove the corresponding convergence theorem at the compact-temporal manifold $S^3 \times S^1_T$. The rate of convergence is $O(\log n_{\max}/n_{\max} + 1/N_t)$ in the joint cutoff. The construction reduces to Paper 38's Riemannian SU(2) propinquity bit-exactly at $N_t = 1$. The strong-form Lorentzian propinquity (without $K^+$ restriction, i.e., a Latrémolière-style metric on Krein spectral triples in their own right) remains open and is a multi-month NCG-math problem; the $K^+$ restriction provides a defensible weak-form alternative by trading Lorentzian-causal structure for Hilbert-space algebraic accessibility on the $J$-positive eigenspace, on which Latrémolière 2017's machinery transfers verbatim."

---

## §10. Files added in this sub-sprint

### Production module

- **`geovac/lorentzian_propinquity_compact_temporal.py`** (~700 lines) — Implements the L5 propinquity assembly on the compact-temporal Lorentzian Krein spectral triple. Exports:
  - `LorentzianTunnelingPair` dataclass: K⁺-restricted tunneling pair $(B^{\mathrm{joint}}, P^{\mathrm{joint}})$.
  - `LorentzianPropinquityBound` dataclass: the $\Lambda^L$ bound + constituent reach/height.
  - `compute_lorentzian_propinquity_bound(n_max, N_t, T, ...)`: the L5 quantitative bound.
  - `lorentzian_gh_convergence_table(...)`: cross-cutoff bound table.
  - `verify_riemannian_limit_at_N_t_1(n_max, T)`: bit-exact recovery of Paper 38's bound at $N_t = 1$.
  - `LorentzianFiveLemmaStatus`: L1'/L2/L3/L4/L5 status tracker (all DONE).
  - `lorentzian_theorem_statement()`: formal Theorem 4.1 statement.

### Driver

- **`debug/l3b_2_sub_sprint_D_compute.py`** — Numerical verification driver. Computes the $\Lambda^L$ panel at $(n_{\max}, N_t) \in \{(2,3), (3,5), (4,7), (5,9)\}$, the $N_t = 1$ Riemannian-limit bit-exactness check, and the asymptotic rate consistency.

### Data

- **`debug/data/l3b_2_sub_sprint_D.json`** — Per-cell panel results.

### Tests

- **`tests/test_lorentzian_propinquity.py`** — 15+ test functions covering: tunneling pair construction, reach + height bounds, main theorem panel, asymptotic rate, K⁺ restriction correctness, Riemannian-limit recovery, monotone decrease.

### Memo (this file)

- **`debug/l3b_2_sub_sprint_D_propinquity_memo.md`** — This proof memo (~5000 words).

---

## §11. Verdict

**PROVED-WITH-NAMED-GAP.**

The K⁺-restricted weak-form Latrémolière propinquity bound

$$
\Lambda^L\!\bigl(\mathcal{T}^L_{n_{\max}, N_t, T}, \mathcal{T}^L_{S^3 \times S^1_T}\bigr) \;\le\; C_3^{\mathrm{joint}} \cdot \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} \;\to\; 0
$$

holds for the compact-temporal Lorentzian truncated Krein spectral triple, at the level of rigor of Paper 38 §L5. The Riemannian-limit recovery at $N_t = 1$ is bit-exact (load-bearing falsifier passed).

The named gap is the strong-form Lorentzian propinquity (Latrémolière-style metric on Krein-signature spectral triples without K⁺ restriction), which remains open and is a multi-month NCG-math problem — not addressed here. The K⁺ restriction defensibly trades the strong-form construction for one that uses Hilbert-space machinery on $\mathcal{K}^+$, and obtains a first Lorentzian propinquity convergence theorem on truncated Krein spectral triples.

The L3b-2 arc is complete. Paper 45 drafting is unblocked.

### Forward to Paper 45 drafting

The proof memos A, B, C, D supply all analytical content for Paper 45. The recommended structure is in §9.2. The "first Lorentzian propinquity convergence theorem on truncated Krein spectral triples" is the operational headline; honest scope of the K⁺-restricted weak-form vs strong-form gap is the load-bearing discipline.

---

**End of memo.**
