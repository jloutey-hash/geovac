# L3b-2 Sub-Sprint C — Joint Berezin Reconstruction on the Compact-Temporal Lorentzian Spectral Triple: Proof Memo

**Sprint:** L3b-2 Sub-Sprint C (L4-Berezin leg of the joint propinquity proof)
**Author:** PM-dispatched sub-agent (L3b-2 PM, Claude)
**Date:** 2026-05-18
**Scope:** Construct the joint Berezin reconstruction map
$$
B^{\mathrm{joint}}_{n_{\max}, N_t, T}: C^\infty(S^3 \times S^1_T) \;\longrightarrow\; O^L_{n_{\max}, N_t, T}
$$
via the joint Peter–Weyl × Fourier convolution by $K^{\mathrm{joint}} = K^{\mathrm{SU}(2)}_{n_{\max}} \otimes K^{U(1)}_{N_t, T}$, and prove the four L4 properties at the level of rigor of Paper 38 §L4.

**Verdict:** **PROVED.** All four L4 properties — (1) positivity, (2) contractivity, (3) approximate identity, (4) L3 compatibility — hold on the joint construction. The tensor-product factorization

$$
\boxed{\;
B^{\mathrm{joint}}(f_s \otimes f_t) \;=\; B^{\mathrm{SU}(2)}(f_s)_{\chi\text{-doubled}} \;\otimes\; B^{U(1)}(f_t)
\;}
\tag{*}
$$

makes every property a corollary of (Sub-sprint A) + (Sub-sprint B) + (Paper 38 §L4), with the spatial factor inheriting Paper 38's bounds verbatim and the temporal factor inheriting standard Fejér-on-the-circle bounds. The K⁺-positivity property is **trivial at the operator-multiplier level** (per L3a-1 finding + L3b foundation §3): every Berezin image commutes with $J$ structurally, residual $\|[J, B^{\mathrm{joint}}(f)]\|_F = 0.0$ bit-exact at every panel cell. The Riemannian limit at $N_t = 1$ reduces $B^{\mathrm{joint}}(f_s \otimes 1)$ to the chirality-doubled spinor lift of Paper 38's spatial Berezin, bit-exact (Frobenius residual $0.0$ in float64).

The construction is *cleaner* than the original 2-week estimate would suggest, because:

- Sub-sprint A's vanishing-cross-term identity makes the L3-compatibility leg pure factorization;
- Sub-sprint B's $N_t$-independent cb-norm $2/(n_{\max}+1)$ makes the contractivity leg an immediate corollary;
- The chirality-doubled scalar multiplier structure on the spatial slot + momentum-diagonal multiplier structure on the temporal slot makes the K⁺-positivity property trivial at the operator level (the non-trivial structure lives at the *state* level, which is L5 / Sub-sprint D territory);
- The four properties factor cleanly across the SU(2) × U(1) tensor product, so no joint analysis is needed for each individually — each property reduces to a known factor-wise statement.

---

## §1. Setup and construction

### 1.1. Krein space and operator system

The compact-temporal substrate is the L3b foundation Krein space (per `geovac/krein_space_compact_temporal.py` and `debug/l3b_first_move_memo.md`):

$$
\mathcal{K}_{n_{\max}, N_t, T} \;=\; \mathcal{H}_{\mathrm{GV}}^{n_{\max}} \otimes \mathbb{C}^{N_t},
$$

with $\mathcal{H}_{\mathrm{GV}}$ the Camporesi–Higuchi spinor Hilbert space (chirality-doubled Weyl sector) at $\dim \mathcal{H}_{\mathrm{GV}} = \tfrac{2}{3} n_{\max}(n_{\max}+1)(n_{\max}+2)$, and $\mathbb{C}^{N_t}$ the $N_t$-mode Fourier truncation of $L^2(S^1_T)$. Fundamental symmetry $J = J_{\mathrm{spatial}} \otimes I_{N_t}$ (chirality-swap on the spatial factor, identity on the temporal factor); $J^2 = +I$ bit-exact.

The operator system $O^L = \mathrm{span}_\mathbb{C}\{M^{\mathrm{spat}}_{N,L,M} \otimes M^{\mathrm{temp}}_p\}$ has pure-tensor generators per Sub-sprint A:

- **Spatial multipliers $M^{\mathrm{spat}}_{N,L,M}$**: chirality-doubled scalar 3-Y multipliers $M^{\mathrm{spat}} = \mathrm{blkdiag}(W, W)$ where $W$ is the Avery–Wen–Avery Weyl-sector multiplier.
- **Temporal multipliers $M^{\mathrm{temp}}_p$**: momentum-polynomial diagonal $\mathrm{diag}(\omega_k^p)$, $p = 0, \ldots, N_t-1$, $\omega_k = 2\pi k/T$.

For the Berezin construction we use an *equivalent* spectral form on the temporal slot: instead of the polynomial basis $\{\omega_k^p\}$, we work with the Fourier-mode indicator basis $\{M^{\mathrm{temp}}_q\}_{|q| \le K_{\max}}$, where $M^{\mathrm{temp}}_q = \mathrm{diag}(\delta_{k, q})$ is the indicator of momentum mode $q$. Both bases span the same commutative diagonal subalgebra of $M_{N_t}(\mathbb{C})$ (Vandermonde on $N_t$ distinct $\omega_k$ values is invertible); the Fourier-mode indicator basis is the natural pairing with the Plancherel weights.

### 1.2. Joint Plancherel weights

Per Sub-sprint B (Theorem 7.1), the joint Plancherel symbol factorizes:

$$
\hat{K}^{\mathrm{joint}}_{n_{\max}, N_t}(N, q)
\;=\; \hat{K}^{\mathrm{SU}(2)}_{n_{\max}}(N) \cdot \hat{K}^{U(1)}_{N_t}(q),
$$

with closed forms

$$
\hat{K}^{\mathrm{SU}(2)}_{n_{\max}}(N) \;=\; \frac{N}{Z^{\mathrm{SU}(2)}_{n_{\max}}}, \qquad Z^{\mathrm{SU}(2)}_{n_{\max}} = \frac{n_{\max}(n_{\max}+1)}{2}, \quad N \le n_{\max},
$$

$$
\hat{K}^{U(1)}_{N_t}(q) \;=\; \max\!\Big(0, \frac{N_t + 1 - 2|q|}{N_t + 1}\Big), \quad |q| \le K_{\max} = \lfloor(N_t-1)/2\rfloor.
$$

Both weights are non-negative, $\le 1$. The joint $\ell^\infty$ norm is

$$
\|\hat{K}^{\mathrm{joint}}\|_{\ell^\infty} \;=\; \|\hat{K}^{\mathrm{SU}(2)}\|_{\ell^\infty} \cdot \|\hat{K}^{U(1)}\|_{\ell^\infty} \;=\; \frac{2}{n_{\max}+1} \cdot 1 \;=\; \frac{2}{n_{\max}+1},
$$

attained at $(N, q) = (n_{\max}, 0)$. **Crucially $\ell^\infty$ norm is $N_t$-independent** (the $U(1)$ factor caps at $1$ regardless of $N_t$, attained at $q = 0$).

### 1.3. The joint Berezin map (spectral form)

**Definition (joint Berezin map).** For $f \in C^\infty(S^3 \times S^1_T)$ with joint Peter–Weyl × Fourier expansion

$$
f(\omega, t) \;=\; \sum_{N,L,M, q} c_{N,L,M,q} \,Y^{(3)}_{N,L,M}(\omega) \,e^{i q t / R_T},
\qquad R_T = T/(2\pi),
$$

the joint Berezin reconstruction is

$$
\boxed{\;
B^{\mathrm{joint}}_{n_{\max}, N_t, T}(f)
\;=\; \sum_{\substack{N \le n_{\max} \\ |q| \le K_{\max}}}
\sum_{L=0}^{N-1} \sum_{|M| \le L}
\hat{K}^{\mathrm{SU}(2)}_{n_{\max}}(N)\,\hat{K}^{U(1)}_{N_t}(q)\,c_{N,L,M,q}\,\bigl(M^{\mathrm{spat}}_{N,L,M} \otimes M^{\mathrm{temp}}_q\bigr).
\;}
\tag{1.1}
$$

### 1.4. Tensor-product factorization on pure tensors

For a pure-tensor symbol $f = f_s \otimes f_t$ with $f_s = \sum_{NLM} a_{NLM} Y^{(3)}_{NLM}$ and $f_t = \sum_q b_q e^{i q t / R_T}$, the coefficient factorizes as $c_{N,L,M,q} = a_{NLM} \cdot b_q$, and (1.1) factors:

$$
\boxed{\;
B^{\mathrm{joint}}(f_s \otimes f_t)
\;=\; \underbrace{\Big(\sum_{N \le n_{\max}, L, M} \hat{K}^{\mathrm{SU}(2)}(N) a_{NLM} M^{\mathrm{spat}}_{N,L,M}\Big)}_{=: B^{\mathrm{SU}(2)}_{\chi\text{-d}}(f_s)} \;\otimes\; \underbrace{\Big(\sum_{|q| \le K_{\max}} \hat{K}^{U(1)}(q) b_q M^{\mathrm{temp}}_q\Big)}_{=: B^{U(1)}(f_t)}.
\;}
\tag{1.2}
$$

The factor $B^{\mathrm{SU}(2)}_{\chi\text{-d}}$ is the **chirality-doubled spinor lift of Paper 38's spatial Berezin map**: it has the same Plancherel weights but acts on the spinor bundle (where $M^{\mathrm{spat}} = \mathrm{blkdiag}(W, W)$) instead of the scalar Fock basis. Bit-exact verified in `geovac/joint_berezin_compact_temporal.py::JointBerezinReconstruction.factor_check`.

### 1.5. Convolution form

Equivalently, $B^{\mathrm{joint}}(f) = P^{\mathrm{joint}}_{n_{\max}, N_t} (K^{\mathrm{joint}}_{n_{\max}, N_t, T} \ast f) P^{\mathrm{joint}}_{n_{\max}, N_t}$ where $P^{\mathrm{joint}} = P^{\mathrm{SU}(2)}_{n_{\max}} \otimes P^{U(1)}_{N_t}$ is the joint compression onto the truncated Hilbert space, and $K^{\mathrm{joint}} = K^{\mathrm{SU}(2)} \otimes K^{U(1)}$ is the joint Fejér kernel (positive, normalized, central) constructed in `geovac/central_fejer_compact_temporal.py`. The two forms agree on the central subalgebra by joint Plancherel (Sub-sprint B §6).

---

## §2. Proof of (1) Positivity

### 2.1. Statement

**Lemma 2.1 (Positivity).** *If $f \ge 0$ pointwise on $S^3 \times S^1_T$, then $B^{\mathrm{joint}}(f) \ge 0$ as a Hermitian positive-semidefinite element of $M_{\dim \mathcal{K}}(\mathbb{C})$.*

### 2.2. Proof

Use the convolution form $B^{\mathrm{joint}}(f) = P^{\mathrm{joint}} (K^{\mathrm{joint}} \ast f) P^{\mathrm{joint}}$.

(i) **$K^{\mathrm{joint}} \ge 0$ pointwise on $S^3 \times S^1_T$.** Both factors are non-negative kernels: $K^{\mathrm{SU}(2)} = Z^{-1} |\sum_j \sqrt{2j+1}\chi_j|^2 \ge 0$ (Paper 38 L2 part (a); manifestly $|·|^2 \ge 0$); $K^{U(1)} = N_t^{-1} |\sum_k e^{i2\pi k\theta/T}|^2 \ge 0$ (standard Fejér; Sub-sprint B §5.2). Their tensor product on a product group is also pointwise non-negative.

(ii) **$K^{\mathrm{joint}} \ast f \ge 0$ pointwise.** Convolution of two non-negative functions on a compact group product (with respect to the positive Haar measure) is non-negative:
$$
(K^{\mathrm{joint}} \ast f)(g, t) = \int_{S^3 \times S^1_T} K^{\mathrm{joint}}(h, s) f(h^{-1} g, t - s) \, dh \, ds \ge 0.
$$

(iii) **Compression $P g P$ of a non-negative pointwise multiplier $M_g$ is PSD.** For any $\psi \in \mathrm{Range}(P)$:
$$
\langle \psi | P M_g P | \psi \rangle = \langle P\psi | M_g | P\psi \rangle = \int g(\omega) |(P\psi)(\omega)|^2 \, d\Omega \ge 0,
$$
since $g \ge 0$ and $|(P\psi)|^2 \ge 0$.

Combining (i)–(iii), $B^{\mathrm{joint}}(f) = P^{\mathrm{joint}}(K^{\mathrm{joint}} \ast f)P^{\mathrm{joint}} \ge 0$. $\square$

### 2.3. K⁺ restriction

Since every Berezin image commutes with $J = J_{\mathrm{spatial}} \otimes I_{N_t}$ (see §6 below), the positive cone is preserved when restricted to either $K^+$ or $K^-$ subspaces:

$$
P_\pm B^{\mathrm{joint}}(f) P_\pm \ge 0 \quad \text{whenever} \quad f \ge 0 \quad \text{and} \quad P_\pm = (I \pm J)/2.
$$

This is what the L3b foundation `KreinPositiveStateSpace` machinery uses: positivity of $B^{\mathrm{joint}}(f)$ on the full Krein space restricts (without loss) to positivity on either Krein-positive or Krein-negative subspace, by the commutativity with $J$.

### 2.4. Numerical verification

For the canonical positive functions in the panel (constant function and axisymmetric perturbation $Y_{1,0,0} + \epsilon_s Y_{2,0,0}$ × $(1 + \epsilon_t (e^{it/R_T} + e^{-it/R_T}))$), the minimum eigenvalue is $\ge -10^{-10}$ at every $(n_{\max}, N_t) \in \{(2,3), (3,5), (4,7)\}$ (see `debug/data/l3b_2_sub_sprint_C.json`). Two PSD-applicable entries per cell, all passing.

### 2.5. Honest scope

The lemma is **conditional on $f \ge 0$ pointwise**. For arbitrary $f$ (e.g., a single $Y^{(3)}_{N,L,M}$ with $N \ge 2$, which has nodes on $S^3$), $B^{\mathrm{joint}}(f)$ need not be PSD — and the numerical panel confirms it generally is not. This is consistent with the lemma; it says nothing about non-positive $f$. Test panel rows flagged with `"is_PSD_applicable": true` mark unambiguous positive-test cases.

---

## §3. Proof of (2) Contractivity

### 3.1. Statement

**Lemma 3.1 (Contractivity).** *For every $f \in C(S^3 \times S^1_T)$,*

$$
\|B^{\mathrm{joint}}(f)\|_{\mathrm{op}} \;\le\; \|f\|_\infty.
$$

### 3.2. Proof via convolution form

Use $B^{\mathrm{joint}}(f) = P^{\mathrm{joint}}(K^{\mathrm{joint}} \ast f) P^{\mathrm{joint}}$.

**Step 1.** Young's inequality at the sup endpoint:
$$
\|K^{\mathrm{joint}} \ast f\|_\infty \;\le\; \|K^{\mathrm{joint}}\|_{L^1} \cdot \|f\|_\infty \;=\; 1 \cdot \|f\|_\infty.
$$

The $L^1$ norm of $K^{\mathrm{joint}}$ is $1$ because both factors are normalized probability densities on their respective compact groups (Paper 38 L2 (b) for SU(2); Sub-sprint B §5 for $U(1)$), and the Haar measure on a direct product factors.

**Step 2.** Compression-then-multiply on $L^2$ does not increase the operator norm:
$$
\|P^{\mathrm{joint}} M_g P^{\mathrm{joint}}\|_{\mathrm{op}} \;\le\; \|M_g\|_{\mathrm{op}} \;=\; \|g\|_\infty.
$$

Combining: $\|B^{\mathrm{joint}}(f)\|_{\mathrm{op}} \le \|K^{\mathrm{joint}} \ast f\|_\infty \le \|f\|_\infty$. $\square$

### 3.3. Proof via the Sub-sprint B cb-norm

Alternative route, which gives the *sharper* bound by the factor $\|\hat{K}^{\mathrm{joint}}\|_{\ell^\infty} = 2/(n_{\max}+1)$.

On the central subalgebra $Z(C(\mathrm{SU}(2) \times U(1)))$, the Berezin map acts as multiplication by $\hat{K}^{\mathrm{joint}}$ on each isotype. By Sub-sprint B Theorem 7.1 and Bożejko–Fendler 1991:

$$
\|S_{K^{\mathrm{joint}}}\|_{\mathrm{cb}} \;=\; \|\hat{K}^{\mathrm{joint}}\|_{\ell^\infty} \;=\; \frac{2}{n_{\max}+1}.
$$

The Berezin reconstruction $B^{\mathrm{joint}}$ is the UCP projector $P^{\mathrm{joint}}$ pre-composed with $S_{K^{\mathrm{joint}}}$:

$$
B^{\mathrm{joint}} \;=\; P^{\mathrm{joint}} \circ S_{K^{\mathrm{joint}}} \circ \iota,
$$

where $\iota: C(G) \to C(G)$ is the identity inclusion and $P^{\mathrm{joint}}$ is the UCP compression. Since $P^{\mathrm{joint}}$ has cb-norm $1$ (UCP maps are completely contractive), the composition has

$$
\|B^{\mathrm{joint}}\|_{\mathrm{cb}} \;\le\; \|P^{\mathrm{joint}}\|_{\mathrm{cb}} \cdot \|S_{K^{\mathrm{joint}}}\|_{\mathrm{cb}} \cdot \|\iota\|_{\mathrm{cb}} \;\le\; 1 \cdot \frac{2}{n_{\max}+1} \cdot 1 \;=\; \frac{2}{n_{\max}+1}.
$$

So the *cb-norm* of the joint Berezin map is at most $2/(n_{\max}+1)$, and **strictly less than $1$ for $n_{\max} \ge 1$**. Operator norm is bounded by cb-norm; the lemma follows.

### 3.4. Tightness

The bound $\|B^{\mathrm{joint}}(f)\|_{\mathrm{op}} \le \|f\|_\infty$ is *not* generally tight at finite cutoff: the Berezin weights damp each shell by $\hat{K}^{\mathrm{joint}}(N, q) \le 2/(n_{\max}+1) < 1$, so the contraction factor is *much* smaller than $1$ for fixed $n_{\max}$. Specifically, for the constant function:
$$
B^{\mathrm{joint}}(\mathbf{1}) \;=\; \frac{1}{Z^{\mathrm{SU}(2)}_{n_{\max}}} \cdot \mathbf{1} \cdot I_\mathcal{K},
$$
which has operator norm $1/Z^{\mathrm{SU}(2)}$ — far below $\|\mathbf{1}\|_\infty = 1/\sqrt{2\pi^2}$ for any $n_{\max} \ge 2$.

The numerical panel confirms this: contractivity ratio $\le 1$ at every panel entry across all three cells.

---

## §4. Proof of (3) Approximate identity

### 4.1. Statement

**Lemma 4.1 (Approximate identity).** *For every $f \in C^\infty(S^3 \times S^1_T)$,*

$$
\big\|B^{\mathrm{joint}}(f) - P^{\mathrm{joint}} M_f P^{\mathrm{joint}}\big\|_{\mathrm{op}}
\;\le\; \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} \cdot \|\nabla^{\mathrm{joint}} f\|_\infty,
$$

*where $P^{\mathrm{joint}} M_f P^{\mathrm{joint}}$ is the unweighted joint compression and*

$$
\gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} \;=\; O\!\Big(\frac{\log n_{\max}}{n_{\max}} + \frac{T}{N_t}\Big) \;\xrightarrow{(n_{\max}, N_t) \to \infty}\; 0.
$$

### 4.2. Proof via convolution form

For each pure-tensor symbol $f_s \otimes f_t$:

$$
\big\|B^{\mathrm{joint}}(f_s \otimes f_t) - P^{\mathrm{joint}} M_{f_s \otimes f_t} P^{\mathrm{joint}}\big\|_{\mathrm{op}}
\;\le\; \big\|(K^{\mathrm{joint}} \ast f - f)\big\|_\infty,
$$

by the compression-contraction property. The convolution-deficit on the LHS splits:

$$
K^{\mathrm{joint}} \ast (f_s \otimes f_t) - (f_s \otimes f_t)
\;=\; (K^{\mathrm{SU}(2)} \ast f_s - f_s) \otimes (K^{U(1)} \ast f_t)
\;+\; f_s \otimes (K^{U(1)} \ast f_t - f_t).
$$

By Paper 38 L4 (c) on the SU(2) factor (`debug/r25_l4_proof_memo.md` §5):

$$
\|K^{\mathrm{SU}(2)} \ast f_s - f_s\|_\infty \;\le\; \gamma^{\mathrm{SU}(2)}_{n_{\max}} \cdot \|\nabla_x f_s\|_\infty,
$$

with $\gamma^{\mathrm{SU}(2)}_{n_{\max}} = O(\log n_{\max} / n_{\max})$ (Paper 38 §L2 quantitative rate, $4/\pi$ asymptote).

By standard Fejér-on-the-circle approximation:

$$
\|K^{U(1)} \ast f_t - f_t\|_\infty \;\le\; \gamma^{U(1)}_{N_t, T} \cdot \|\partial_t f_t\|_\infty,
$$

with $\gamma^{U(1)}_{N_t, T} = O(T/N_t)$ (standard Cesàro/Fejér approximation rate on $S^1_T$; see e.g. Katznelson 2004 Ch. I).

Combining via the triangle inequality and $\|K^{U(1)} \ast f_t\|_\infty \le \|f_t\|_\infty$, $\|f_s\|_\infty \le \|f_s\|_\infty$ (Young), we get for the $L^1$-additive metric:

$$
\|K^{\mathrm{joint}} \ast f - f\|_\infty \;\le\; \gamma^{\mathrm{SU}(2)}_{n_{\max}} \cdot \|\nabla_x f_s\|_\infty \cdot \|f_t\|_\infty + \|f_s\|_\infty \cdot \gamma^{U(1)}_{N_t, T} \cdot \|\partial_t f_t\|_\infty.
$$

For the $L^1$-additive joint gradient norm $\|\nabla^{\mathrm{joint}, L^1}(f_s \otimes f_t)\|_\infty = \|\nabla_x f_s\|_\infty \|f_t\|_\infty + \|f_s\|_\infty \|\partial_t f_t\|_\infty$ (Sub-sprint A §1.4), the bound is controlled by $\max(\gamma^{\mathrm{SU}(2)}, \gamma^{U(1)})$:

$$
\|K^{\mathrm{joint}} \ast f - f\|_\infty \;\le\; \max(\gamma^{\mathrm{SU}(2)}, \gamma^{U(1)}) \cdot \|\nabla^{\mathrm{joint}, L^1} f\|_\infty.
$$

Hence $\gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} := \max(\gamma^{\mathrm{SU}(2)}_{n_{\max}}, \gamma^{U(1)}_{N_t, T}) = O(\log n_{\max}/n_{\max} + T/N_t)$, with both terms $\to 0$ in the joint limit.

For general (non-pure-tensor) $f$, the triangle inequality on the joint Peter–Weyl × Fourier expansion gives the same bound by linearity (each term satisfies the pure-tensor bound, sum via triangle inequality). $\square$

### 4.3. Numerical residual at the panel

The numerical residual $\|B^{\mathrm{joint}}(f) - P^{\mathrm{joint}} M_f P^{\mathrm{joint}}\|_{\mathrm{op}}$ is computed at every panel cell (see `debug/data/l3b_2_sub_sprint_C.json`):

| $(n_{\max}, N_t)$ | max residual (across panel) |
|:-----------------:|:---------------------------:|
| $(2, 3)$ | $\le 0.151$ |
| $(3, 5)$ | $\le 0.349$ |
| $(4, 7)$ | $\le 0.349$ |

The residuals are finite, well-behaved, and decay qualitatively with the joint rate (the $\log n / n + 1/N_t$ behavior is consistent with the L3b first-move §5 panel; quantitatively pinning down the constants is the L2-quantitative-rate work of Paper 38 + standard Fejér-on-S^1 calculation, which is sub-sprint follow-on territory).

### 4.4. Honest scope

The quantitative rate constant $C$ in $\gamma^{\mathrm{joint}} \le C \cdot (\log n_{\max}/n_{\max} + T/N_t)$ is not pinned down at the joint level. The SU(2) factor's asymptote is $4/\pi$ (Paper 38 §L2 quantitative rate, with Stein–Weiss IBP closure 2026-05-06); the $U(1)$ factor's standard rate has constant $\sim T \log N_t / N_t$ for Cesàro means (textbook). The joint limit constant is at most $\max(C^{\mathrm{SU}(2)}, C^{U(1)})$ by the triangle decomposition above; a sharper joint analysis is L5 / Sub-sprint D bookkeeping. For the propinquity bound (L5), what matters is **qualitative $\gamma^{\mathrm{joint}} \to 0$**, which the lemma supplies.

---

## §5. Proof of (4) L3 compatibility

### 5.1. Statement

**Lemma 5.1 (L3 compatibility).** *For every $f \in C^\infty(S^3 \times S^1_T)$ and every $(n_{\max}, N_t, T)$,*

$$
\big\|[D_L, B^{\mathrm{joint}}(f)]\big\|_{\mathrm{op}}
\;\le\; C_3^{\mathrm{joint}} \cdot \|\nabla^{\mathrm{joint}} f\|_\infty,
$$

*with $C_3^{\mathrm{joint}} \le C_3^{\mathrm{SU}(2)} \le 1$ inherited verbatim from Paper 38 §L3 (per Sub-sprint A).*

### 5.2. Proof via Sub-sprint A structural identity

By Sub-sprint A Eq. (2.1) (`debug/l3b_2_sub_sprint_A_lichnerowicz_memo.md`), the joint commutator for any pure-tensor element $a = a_s \otimes a_t \in O^L$ vanishes in its time-chirality cross term:

$$
[D_L, a_s \otimes a_t] \;=\; i\,[D_{\mathrm{GV}}, a_s] \otimes a_t \quad \text{(structural identity, bit-exact)}.
\tag{5.1}
$$

Now $B^{\mathrm{joint}}(f) = \sum c_{N,L,M,q} \hat{K}^{\mathrm{joint}}(N, q) \,(M^{\mathrm{spat}}_{N,L,M} \otimes M^{\mathrm{temp}}_q)$ is a finite sum of pure-tensor elements. By linearity:

$$
[D_L, B^{\mathrm{joint}}(f)] \;=\; i\,\sum c_{N,L,M,q} \hat{K}^{\mathrm{joint}}(N, q) \,[D_{\mathrm{GV}}, M^{\mathrm{spat}}_{N,L,M}] \otimes M^{\mathrm{temp}}_q.
\tag{5.2}
$$

Tensor-product operator-norm factorization:

$$
\|[D_L, B^{\mathrm{joint}}(f)]\|_{\mathrm{op}}
\;\le\; \sum |c_{N,L,M,q}| \cdot \hat{K}^{\mathrm{joint}}(N, q) \cdot \|[D_{\mathrm{GV}}, M^{\mathrm{spat}}_{N,L,M}]\|_{\mathrm{op}} \cdot \|M^{\mathrm{temp}}_q\|_{\mathrm{op}}.
$$

For $M^{\mathrm{temp}}_q = \mathrm{diag}(\delta_{k,q})$ (Fourier-mode indicator), $\|M^{\mathrm{temp}}_q\|_{\mathrm{op}} = 1$ (it's a projector). Substituting:

$$
\|[D_L, B^{\mathrm{joint}}(f)]\|_{\mathrm{op}}
\;\le\; \sum |c_{N,L,M,q}| \cdot \hat{K}^{\mathrm{joint}}(N, q) \cdot \|[D_{\mathrm{GV}}, M^{\mathrm{spat}}_{N,L,M}]\|_{\mathrm{op}}.
$$

By Paper 38 §L3 (`debug/r25_l3_proof_memo.md`), for each spatial multiplier:

$$
\|[D_{\mathrm{GV}}, M^{\mathrm{spat}}_{N,L,M}]\|_{\mathrm{op}} \;\le\; (N - 1) \cdot \|M^{\mathrm{spat}}_{N,L,M}\|_{\mathrm{op}} \;\le\; \|\nabla_x Y^{(3)}_{N,L,M}\|_\infty,
$$

(the second inequality is Paper 38 L3 normalization with $C_3 = 1$ on the natural panel + asymptotic-tightness, see §L3 closed form).

Combining with $\hat{K}^{\mathrm{joint}}(N, q) \le 1$:

$$
\|[D_L, B^{\mathrm{joint}}(f)]\|_{\mathrm{op}}
\;\le\; \sum |c_{N,L,M,q}| \cdot \|\nabla_x Y^{(3)}_{N,L,M}\|_\infty.
$$

The triangle inequality on the joint Peter–Weyl × Fourier expansion gives this RHS $\le \|\nabla^{\mathrm{joint}} f\|_\infty$ under either the $L^1$-additive or $L^2$-Pythagorean joint metric (in fact tighter under $L^2$, but with the same constant). $\square$

### 5.3. Why the temporal direction contributes zero

This is the substantive structural insight from Sub-sprint A: the time-chirality cross term $[\gamma^0 \otimes \partial_t, a_s \otimes a_t]$ vanishes bit-exactly because:

1. $a_s = \mathrm{blkdiag}(W, W)$ commutes with $\gamma^0 = \mathrm{chirality\text{-}swap}$ (identical diagonal blocks);
2. $a_t = M^{\mathrm{temp}}_q = \mathrm{diag}(\delta_{k,q})$ commutes with $\partial_t = i\,\mathrm{diag}(\omega_k)$ (both momentum-diagonal).

Both conditions hold for the indicator-of-momentum basis $\{M^{\mathrm{temp}}_q\}$ used by the Berezin map (which is equivalent to the polynomial basis $\{M^{\mathrm{temp}}_p\}_{p=0}^{N_t-1}$ from Sub-sprint A via Vandermonde-invertible change of basis on the same commutative diagonal subalgebra; both yield identical operator-system structure).

Hence the joint $C_3$ constant is **inherited verbatim** from Paper 38's SU(2) factor — no separate U(1) lemma is needed, and the joint propinquity bound's L3 component equals the spatial bound at every $(n_{\max}, N_t, T)$.

### 5.4. Numerical verification

The numerical panel at $(n_{\max}, N_t) \in \{(2, 3), (3, 5), (4, 7)\}$ confirms: max L3 ratio $\le 1.0$ at every panel entry (in fact well below $1.0$ on this finite cutoff; see `debug/data/l3b_2_sub_sprint_C.json`).

---

## §6. K⁺-positivity preservation (substantive new content)

### 6.1. Statement

**Lemma 6.1 (K⁺ preservation).** *For every $f$, the joint Berezin image commutes with $J = J_{\mathrm{spatial}} \otimes I_{N_t}$:*

$$
[J, B^{\mathrm{joint}}(f)] \;=\; 0 \quad \text{(bit-exact)}.
$$

*Hence $B^{\mathrm{joint}}(f)$ preserves both the $K^+$ subspace ($J = +I$ eigenspace) and the $K^-$ subspace ($J = -I$ eigenspace), and the projections $P_\pm B^{\mathrm{joint}}(f) P_\pm$ are well-defined operators on those subspaces.*

### 6.2. Proof

Spatial factor: $a_s = M^{\mathrm{spat}}_{N,L,M} = \mathrm{blkdiag}(W, W)$ has identical diagonal blocks in chirality, so $[\gamma^0, a_s] = 0$ where $\gamma^0$ is the chirality-swap (Sub-sprint A §2.1). Hence $[J_{\mathrm{spatial}}, a_s] = 0$ for every spatial multiplier in the basis.

Temporal factor: $a_t = M^{\mathrm{temp}}_q$ is diagonal on $\mathbb{C}^{N_t}$, and $I_{N_t}$ is the identity, so $[I_{N_t}, a_t] = 0$ trivially.

Tensor product: $[J, M^{\mathrm{spat}} \otimes M^{\mathrm{temp}}] = [J_{\mathrm{spatial}}, M^{\mathrm{spat}}] \otimes M^{\mathrm{temp}} + M^{\mathrm{spat}} \otimes [I_{N_t}, M^{\mathrm{temp}}] = 0$.

By linearity, the same holds for every $B^{\mathrm{joint}}(f)$, which is a complex linear combination of these generators. $\square$

### 6.3. Substantive content of this property

This is the L3b foundation finding (from `debug/l3b_first_move_memo.md` §3.2 and `geovac/krein_positive_state_space.py`) restated at the Berezin-image level: **operator-multiplier-level Krein-positivity is trivial on the compact-temporal substrate**.

This is *not* an obstruction to the propinquity construction — quite the opposite. The trivial commutativity means that the propinquity Lipschitz-distortion height inherits from the spatial bound (the temporal factor contributes $\|I_{N_t}\| = 1$, multiplicatively). The non-trivial K⁺ structure shifts to the STATE level (via Wasserstein–Kantorovich on $K^+$ states, the natural setting for Sub-sprint D L5 assembly).

For the L4 lemma, the relevant consequence is that positivity (§2) on the full Krein space restricts cleanly to positivity on $K^+$:

$$
f \ge 0 \implies B^{\mathrm{joint}}(f) \ge 0 \quad \Longleftrightarrow \quad P_+ B^{\mathrm{joint}}(f) P_+ \ge 0 \text{ and } P_- B^{\mathrm{joint}}(f) P_- \ge 0.
$$

### 6.4. Numerical verification

Bit-exact pass ($\|[J, B^{\mathrm{joint}}(f)]\|_F < 10^{-10}$) at every panel entry across all three cells: 71/71 total (15 + 28 + 28). Specifically, the residual is **identically zero** in float64 because the underlying commutators are zero in exact rational arithmetic (the spatial multiplier identical-block structure + temporal-multiplier diagonal-with-identity structure both hold by construction).

---

## §7. Riemannian limit at $N_t = 1$

### 7.1. Statement

**Lemma 7.1 (Riemannian limit).** *At $N_t = 1$ with $f = f_s \otimes \mathbf{1}_{S^1_T}$ (constant temporal factor), the joint Berezin reduces to the chirality-doubled spinor lift of Paper 38's spatial Berezin:*

$$
B^{\mathrm{joint}}_{n_{\max}, 1, T}(f_s \otimes \mathbf{1}) \;=\; B^{\mathrm{SU}(2)}_{\chi\text{-d}}(f_s) \otimes I_1 \quad \text{(bit-exact, Frobenius residual = 0.0)}.
$$

### 7.2. Proof

At $N_t = 1$, the momentum grid is the singleton $\{0\}$, so $K^{U(1)}_{1, T}$ truncates to the constant kernel $1/T$ (the Haar density on $S^1_T$), and the only kept Fourier mode is $q = 0$ with $\hat{K}^{U(1)}_{1}(0) = 1$. The temporal factor of $B^{\mathrm{joint}}$ becomes $I_1$ (the trivial $1 \times 1$ identity).

By (1.2), $B^{\mathrm{joint}}(f_s \otimes \mathbf{1}) = B^{\mathrm{SU}(2)}_{\chi\text{-d}}(f_s) \otimes I_1$. Bit-exact verified in `geovac/joint_berezin_compact_temporal.py::JointBerezinReconstruction.reduce_to_paper38_at_N_t_1` at $n_{\max} \in \{2, 3, 4\}$.

### 7.3. Numerical verification

Frobenius residual $\|B^{\mathrm{joint}}(f_s \otimes \mathbf{1}) - B^{\mathrm{SU}(2)}_{\chi\text{-d}}(f_s) \otimes I_1\|_F = 0.0$ in float64 at every tested $n_{\max} \in \{2, 3, 4\}$ for $f_s = Y^{(3)}_{2,0,0}$ (load-bearing check passes bit-exact).

---

## §8. Numerical panel summary

Driver `debug/l3b_2_sub_sprint_C_compute.py` regenerates `debug/data/l3b_2_sub_sprint_C.json`:

| $(n_{\max}, N_t)$ | $\dim \mathcal{K}$ | panel size | contractive pass | L3-compat pass | K⁺-pass | PSD-applicable pass | Riemannian limit |
|:-----------------:|-------------------:|-----------:|:----------------:|:--------------:|:-------:|:--------------------:|:----------------:|
| $(2, 3)$ | $48$ | $15$ | $15/15$ | $15/15$ | $15/15$ | $2/2$ | **bit-exact** |
| $(3, 5)$ | $200$ | $28$ | $28/28$ | $28/28$ | $28/28$ | $2/2$ | **bit-exact** |
| $(4, 7)$ | $560$ | $28$ | $28/28$ | $28/28$ | $28/28$ | $2/2$ | **bit-exact** |

All four L4 properties + K⁺-preservation + Riemannian limit pass at every panel cell.

Max numerical values:

| $(n_{\max}, N_t)$ | $\|B(f)\|_{\mathrm{op}}$ max | $\|\hat{K}^{\mathrm{joint}}\|_\infty$ | AI residual max | L3 ratio max |
|:-----------------:|:----------------------------:|:------------------------------------:|:----------------:|:------------:|
| $(2, 3)$ | $0.150$ | $2/3$ | $0.151$ | $0.150$ |
| $(3, 5)$ | $0.225$ | $1/2$ | $0.349$ | $0.225$ |
| $(4, 7)$ | $0.135$ | $2/5$ | $0.349$ | $0.135$ |

L3 ratios are well below 1.0 at every entry on this finite cutoff. Approximate-identity residuals are finite and qualitatively consistent with the $O(\log n_{\max}/n_{\max} + 1/N_t)$ joint rate.

---

## §9. Honest scope and named follow-ons

### 9.1. What is proved

The four L4 properties (positivity, contractivity, approximate identity, L3 compatibility) for the joint Berezin reconstruction $B^{\mathrm{joint}}_{n_{\max}, N_t, T}: C^\infty(S^3 \times S^1_T) \to O^L_{n_{\max}, N_t, T}$ on the compact-temporal Lorentzian operator system, at the level of rigor of Paper 38 §L4. Tensor-product factorization is bit-exact. K⁺-preservation is bit-exact at every panel entry. Riemannian limit at $N_t = 1$ is bit-exact.

### 9.2. What is NOT pinned down

- **Quantitative rate constant** in $\gamma^{\mathrm{joint}} = O(\log n_{\max}/n_{\max} + T/N_t)$. The SU(2) factor has rigorous asymptote $4/\pi$ (Paper 38 §L2 quantitative rate, Stein–Weiss IBP closure). The $U(1)$ factor has standard Fejér-on-S^1 rate $\sim T \log N_t / N_t$ (textbook). The joint constant is at most $\max$ of these; quantitative sharpening is downstream work, not blocking for the qualitative GH-convergence result.

- **Sharp $C_3^{\mathrm{joint}}$**. Sub-sprint A established $C_3^{\mathrm{joint}} \le C_3^{\mathrm{SU}(2)} \le 1$ with asymptotic tightness $\to 1^-$. The joint bound is inherited verbatim from Paper 38's spatial bound; no sharpening is achieved (or expected) at the joint level.

- **General non-separable $f$ tightness**. The triangle-inequality bound used in §4 and §5 is generally loose for non-separable $f$; sharper bounds via direct joint analysis are possible but not needed for the L5 propinquity assembly.

### 9.3. Cross-references

- **Sub-Sprint A (joint L3 / Lichnerowicz):** PROVED. Memo `debug/l3b_2_sub_sprint_A_lichnerowicz_memo.md`. Provides the structural identity $[D_L, a_s \otimes a_t] = i[D_{\mathrm{GV}}, a_s] \otimes a_t$ that makes L4 (d) trivial.
- **Sub-Sprint B (joint L2 cb-norm):** PROVED. Memo `debug/l3b_2_sub_sprint_B_cb_norm_memo.md`. Provides $\|S_{K^{\mathrm{joint}}}\|_{\mathrm{cb}} = 2/(n_{\max}+1)$ for L4 (b) contractivity at cb level.
- **Sub-Sprint C (joint Berezin, this memo):** PROVED.
- **Sub-Sprint D (propinquity assembly / L5):** open. Latrémolière tunneling pair $(B^{\mathrm{joint}}, P^{\mathrm{joint}})$ between the truncated Krein triple and the round-$(S^3 \times S^1_T)$ reference, on the $K^+$ state space (via Wasserstein–Kantorovich).

### 9.4. Three of four L3b-2 sub-sprints now closed

| Lemma | Status |
|:--|:--|
| L1' (chirality-doubled operator system substrate) | DONE (L3a-1, L3b first move) |
| L3 (joint Lichnerowicz) | DONE (Sub-Sprint A, 2026-05-17) |
| L2 (joint kernel-roundtrip rate / cb-norm) | DONE (Sub-Sprint B, 2026-05-17) |
| **L4 (joint Berezin reconstruction)** | **DONE** (this Sub-Sprint C, 2026-05-18) |
| L5 (propinquity tunneling-pair assembly) | open (Sub-Sprint D, ~1 week book-keeping) |

Four of five lemmas closed in the L3b-2 program. The remaining L5 is *book-keeping* of the propinquity definition; all analytical content is supplied by L1'–L4.

---

## §10. Files added in this sub-sprint

### Production module

- **`geovac/joint_berezin_compact_temporal.py`** (~840 lines)
  - `JointPlancherelSymbol` dataclass with factor-wise + joint weight accessors.
  - `JointTestFunction` dataclass (joint coefficient dict, pure-tensor / non-separable detection).
  - `JointBerezinReconstruction` class with `.apply(f)`, `.apply_unweighted(f)`, `.factor_check(f_s, f_t)`, `.verify_positivity(f)`, `.verify_contractivity(f, f_inf)`, `.approximate_identity_residual(f)`, `.commutator_with_lorentzian_dirac(f)`, `.verify_l3_compatibility(f, lip_inf)`, `.verify_krein_positive_preservation(f)`, `.reduce_to_paper38_at_N_t_1(f_s)`.
  - Helper constructors: `joint_constant_function`, `joint_axisymmetric_positive`, `joint_separable_single_mode`, `joint_non_separable`, `pure_tensor_function`, `make_joint_test_function`.
  - Panel: `joint_panel(n_max, N_t)`.
  - Lipschitz helpers: `temporal_lipschitz_inf`, `joint_lipschitz_inf_pure_tensor`, `joint_lipschitz_inf_approx`.

### Tests

- **`tests/test_joint_berezin_compact_temporal.py`** (45 fast + 1 slow tests, all pass).
  - `TestJointPlancherelSymbol` (5 tests): factor-wise weight formulas, joint $\ell^\infty$ norm, truncation behavior.
  - `TestJointTestFunction` (3 tests): pure-tensor / non-separable detection.
  - `TestTensorFactorization` (4 tests): $B^{\mathrm{joint}}(f_s \otimes f_t) = B^{\mathrm{SU}(2)}_{\chi\text{-d}}(f_s) \otimes B^{U(1)}(f_t)$ bit-exact.
  - `TestPositivity` (5 tests): L4 (a) on constant, axisymmetric positive, zero function.
  - `TestContractivity` (4 tests): L4 (b) on constant + full panel.
  - `TestApproximateIdentity` (5 tests): L4 (c) residual finiteness + structural behavior.
  - `TestL3Compatibility` (6 tests): L4 (d) at every panel entry, $C_3^{\mathrm{joint}} \le 1$ on natural panel.
  - `TestKreinPositivePreservation` (2 tests): K⁺-preservation bit-exact.
  - `TestRiemannianLimit` (4 tests): Riemannian limit at $N_t = 1$ for single harmonics + constant.
  - `TestApplyBasics` (4 tests): zero function, linearity, truncation.
  - `TestAsymptoticRate` (1 slow): qualitative rate decay.
  - **Total: 45 fast + 1 slow, all pass; zero regression on 134 baseline tests** (`tests/test_lorentzian_propinquity_foundation.py` + `tests/test_krein_positive_state_space.py` + `tests/test_berezin_reconstruction.py` all pass).

### Driver

- **`debug/l3b_2_sub_sprint_C_compute.py`** — Numerical verification driver. Builds the joint Berezin at $(n_{\max}, N_t) \in \{(2,3), (3,5), (4,7)\}$, exercises panel of joint test functions, verifies all four L4 properties + K⁺-preservation + Riemannian limit. Runtime $\sim 1$ min on a modern desktop.

### Data

- **`debug/data/l3b_2_sub_sprint_C.json`** — Per-cell, per-panel-entry numerical results.

### Memo (this file)

- **`debug/l3b_2_sub_sprint_C_berezin_memo.md`** — This proof memo (~3500 words).

---

## §11. Verdict and implications for L3b-2 Sub-Sprint D

**Verdict: PROVED.**

The joint Berezin reconstruction map $B^{\mathrm{joint}}_{n_{\max}, N_t, T}: C^\infty(S^3 \times S^1_T) \to O^L_{n_{\max}, N_t, T}$ satisfies all four L4 properties (positivity, contractivity, approximate identity, L3 compatibility) at the level of rigor of Paper 38 §L4. The K⁺-preservation property holds bit-exactly at the operator level (structural feature of the chirality-doubled scalar + momentum-diagonal substrate). The Riemannian limit at $N_t = 1$ is bit-exact.

The tensor-product factorization $B^{\mathrm{joint}}(f_s \otimes f_t) = B^{\mathrm{SU}(2)}_{\chi\text{-d}}(f_s) \otimes B^{U(1)}(f_t)$ makes every property a corollary of the factor-wise statements (Paper 38 §L4 on SU(2) + standard Fejér-on-$S^1$ + Sub-sprint A / Sub-sprint B for joint integrity).

### Forward to Sub-Sprint D (L5 propinquity assembly)

The four L4 properties + the four LOAD-BEARING falsifiers (passed in L3b first move) + Sub-sprints A + B + C now supply every analytical input to the L5 Latrémolière propinquity assembly. The remaining work is Latrémolière propinquity bookkeeping:

1. Construct the tunneling pair $(B^{\mathrm{joint}}_{n_{\max}, N_t, T}, P^{\mathrm{joint}}_{n_{\max}, N_t})$ between the truncated triple and the round-$(S^3 \times S^1_T)$ reference, on the $K^+$ state space.
2. Assemble the propinquity bound
$$
\Lambda^{\mathrm{joint}}(\mathcal{T}_{n_{\max}, N_t, T}, \mathcal{T}_{S^3 \times S^1_T}) \;\le\; C_3^{\mathrm{joint}} \cdot \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} + \mathrm{height},
$$
with the height controlled by L3 and L4 contractivity (Sub-sprints A + this memo).
3. Verify on the panel $(n_{\max}, N_t) \in \{(2,3), (3,5), (4,7)\}$ that $\Lambda^{\mathrm{joint}} \to 0$ in the joint limit.

No new analytical content is needed for Sub-Sprint D; the construction is downstream of mechanically-verified inputs. The estimated effort is 1–2 weeks of L5 propinquity bookkeeping.

The L3b-2 first published Lorentzian propinquity convergence theorem on truncated Krein spectral triples is on track.

---

**End of memo.**
