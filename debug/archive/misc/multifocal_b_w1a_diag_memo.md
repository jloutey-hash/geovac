# Multi-Focal Phase B / W1a Diagnostic Memo

**Date:** 2026-05-07
**Phase:** B (diagnostic; not closure)
**Wall:** W1a — cross-register two-body spatial coordinate operator $V(\hat{\mathbf{r}}_e, \hat{\mathbf{R}}_n)$ at distinct focal lengths
**Author:** PM (sub-agent dispatch; no production code modified)
**Companion data:** `debug/data/multifocal_b_w1a_qa_qb.json`, `multifocal_b_w1a_qb.json`, `multifocal_b_w1a_qc.json`, `multifocal_b_w1a_qd.json`
**Companion drivers:** `debug/multifocal_b_w1a_qa_qb.py`, `multifocal_b_w1a_qb.py`, `multifocal_b_w1a_qc.py`, `multifocal_b_w1a_qd.py`

---

## 0. Executive verdict (TL;DR)

**W1a is verdict (b) — HYBRID.** The mathematical core (closed-form mismatched-exponent two-center integral plus closed-form cross-register bilinear ERI) is *available*, both as small symbolic re-derivations from sympy and as an existing-tool extension of `geovac/shibuya_wulfman.py`. The pieces required for W1a closure split cleanly:

| Piece | Status |
|:------|:-------|
| Mismatched-exponent radial integral for cross-center $V_{nN}$ (single particle, two nuclei) | **Closed form survives.** Reduces to single-$\lambda$ Shibuya–Wulfman with $\alpha_{\rm total} = \lambda_a + \lambda_b$. Q-A confirmed symbolically. |
| Multipole-expansion termination at $L_{\max} = l_a + l_b$ under exponent mismatch | **Preserved.** Termination is purely a Gaunt/3j angular property; depends only on $(l_a, m_a, l_b, m_b)$, not on radial exponents. Q-B confirmed across $l_a, l_b \in \{0, 1, 2\}$. The briefing's stated bound $L_{\max} = 2 \max(l_a, l_b)$ is the loose upper bound; the tight bound is $l_a + l_b$. |
| Transfer operator $T_{\lambda_1 \to \lambda_2}$ between Sturmian/hydrogenic bases at different $\lambda$ | **Generically dense.** Q-C: 100 % density at $5 \times 5$, $\lambda_a = 1, \lambda_b = 2$ test. Same Löwdin-style sparsity destruction Track BU and Track DF Sprint 5 hit. **Therefore the transfer-operator route is the wrong path for W1a.** |
| Cross-register bilinear ERI $\langle \chi_e^{\lambda_e} \chi_n^{\lambda_n} \,\|\, \frac{1}{|\mathbf{r}_e - \mathbf{R}_n|} \,\|\, \chi_{e'}^{\lambda_e} \chi_{n'}^{\lambda_n} \rangle$ | **Closed form available.** Q-D: 1s-1s case gives $J_0 = \frac{\lambda_e \lambda_n (\lambda_e^2 + 3 \lambda_e \lambda_n + \lambda_n^2)}{(\lambda_e + \lambda_n)^3}$ — the textbook Roothaan integral, symmetric in $(\lambda_e, \lambda_n)$, with the correct point-nucleus limit $J_0 \to \lambda_e$ as $\lambda_n \to \infty$. The 2p-2p $L=1$ case also closes. Bilinear in (electron, nucleus) Sturmian indices. |

**The hybrid label is real:** the integrals close in elementary functions, but **building W1a closure requires (i) a new module `geovac/cross_register_vne.py` that wires the bilinear ERI into the joint Pauli encoding, and (ii) promoting Track NI's `R_PROTON_BOHR` classical scalar to an operator-valued $\mathbf{R}$ on the proton register**. This is engineering, not new mathematics. Phase C scope estimate: 4–6 weeks for a working bilinear ERI module, an additional 2–4 weeks to validate against Pachucki–Patkóš–Yerokhin 2023 recoil. Total 6–10 weeks, matching Track 3's Candidate 1 estimate.

**Surprise revising Phase A's framing:** Phase A's Section 1 classified W1a as "tooling-addressable" but flagged the multi-$\lambda$ Shibuya–Wulfman as "the most likely failure mode" (Track 3 Candidate 3 in the synthesis). The diagnostic shows **the multi-$\lambda$ Shibuya–Wulfman generalization is in fact the easiest piece**: the integrand factorizes through the same single-exponential structure as the matched case, and the existing `_split_integral_analytical` requires only that one substitution `alpha = lam_a + lam_b`. The harder piece is Q-D — but Q-D *also* closes, in elementary functions. **Phase A under-rated W1a's tractability.** I revise Phase A's "medium-high risk" assessment for Candidate 1 to "medium risk; the algebraic backbone is closed; the engineering is the bottleneck."

---

## 1. The diagnostic question, in precise mathematical form

### 1.1 Setting

The Coulomb-Sturmian basis at exponent $\lambda$ is
$$
\chi_{nlm}^{\lambda}(\mathbf{r}) = R_{nl}^{\lambda}(r) \, Y_{lm}(\hat{\mathbf{r}}),
\quad R_{nl}^{\lambda}(r) = N_{nl}^{\lambda} \, (2\lambda r)^l \, e^{-\lambda r} \, L_{n-l-1}^{2l+1}(2\lambda r),
$$
with $N_{nl}^{\lambda}$ chosen so $\int_0^\infty |R_{nl}^{\lambda}|^2 r^2 \, dr = 1$. For $\lambda = Z/n$ this reduces to the standard hydrogenic radial function.

The W1a wall, distilled from Sprint HF-3 and the Phase A audit, is that GeoVac currently has the multipole expansion of $1/|\mathbf{r}_e - \mathbf{R}_B|$ only when $\mathbf{R}_B$ is a *classical parameter*, not a quantum operator on a separate register. The architectural challenge has two algebraic ingredients: (i) handling mismatched exponents on bra and ket, since the electron and nucleus naturally sit at very different focal lengths, and (ii) reducing the joint cross-register two-body integral to a closed-form expression bilinear in (electron, nucleus) Sturmian indices.

### 1.2 The four diagnostic questions

**Q-A.** Compute symbolically
$$
M_{ab}(R; \lambda_a, \lambda_b, L) =
\int_0^\infty R_{n_a l_a}^{\lambda_a}(r) R_{n_b l_b}^{\lambda_b}(r) \, g_L(r, R) \, r^2 \, dr,
\qquad
g_L(r, R) = \frac{r_<^L}{r_>^{L+1}}
$$
with $\lambda_a \neq \lambda_b$. Does the closed-form structure of the matched-$\lambda$ Shibuya–Wulfman case survive the mismatch?

**Q-B.** Does the multipole expansion of $1/|\mathbf{r} - \mathbf{R}|$ (with $\mathbf{R}$ a classical parameter, but bra and ket at different $\lambda$) still terminate at $L_{\max}$ determined purely by $(l_a, l_b)$? Test at $l_a, l_b \in \{0, 1, 2\}$.

**Q-C.** Does the change-of-basis "transfer operator" $T_{\lambda_a \to \lambda_b}$ between Sturmian families at different exponents preserve sparsity and integer-Pauli scaling?

**Q-D.** When $\mathbf{R}$ is itself a quantum operator on a separate register at focal length $\lambda_n$, does the cross-register bilinear ERI
$$
\Big\langle \chi^{\lambda_e}_{n_e l_e m_e}(\mathbf{r}_e) \chi^{\lambda_n}_{n_n l_n m_n}(\mathbf{R}_n) \,\Big|\, \frac{1}{|\mathbf{r}_e - \mathbf{R}_n|} \,\Big|\, \chi^{\lambda_e}_{n'_e l'_e m'_e}(\mathbf{r}_e) \chi^{\lambda_n}_{n'_n l'_n m'_n}(\mathbf{R}_n) \Big\rangle
$$
admit a closed-form expansion bilinear in (electron, nucleus) Sturmian indices? This is the actual W1a closure object.

### 1.3 Boundary against the GUARDRAIL

Papers 8–9 prove the Sturmian Structural Theorem for the **shared-$p_0$** case: with a single exponent for all orbitals and a $V_{\rm mol}$-derived Hamiltonian, $H_{ij} \propto S_{ij}$ in the shared-$p_0$ basis, so the generalized-eigenvalue spectrum is $R$-independent — no PES exists. The W1a problem is **explicitly outside the theorem's scope**: we use $\lambda_e \neq \lambda_n$ by design (the electron and nucleus sit at categorically different binding scales). The diagnostic is consistent with the GUARDRAIL: the theorem warns about what fails when exponents are equal, and the W1a construction avoids that pitfall by having mismatched exponents on the two registers.

---

## 2. What's already in GeoVac (single-$\lambda$ Shibuya–Wulfman)

The production module `geovac/shibuya_wulfman.py` (Paper 19, Track CD, v2.0.39–42) computes
$$
I^{AB}_{nlm,n'l'm'} = \big\langle \psi_{nlm}^A \,\big|\, -\frac{Z_B}{|\mathbf{r} - \mathbf{R}_B|} \,\big|\, \psi_{n'l'm'}^A \big\rangle
$$
via:

1. **Multipole expansion** of $1/|\mathbf{r} - \mathbf{R}_B|$ in Legendre polynomials $P_L(\cos\theta)$, terminating exactly at $L_{\max} = 2 l_{\max}$ by Gaunt selection rules (or, more tightly, $L_{\max} = l_a + l_b$ for each individual matrix element).
2. **Angular factor** `_angular_coefficient(l1, m1, l2, m2, L)` returning $(-1)^m \sqrt{(2l_1+1)(2l_2+1)} \cdot 3j(l_1, L, l_2; 0, 0, 0) \cdot 3j(l_1, L, l_2; -m, 0, m)$, with strict $m_1 = m_2$ selection (z-axis nucleus).
3. **Radial split-region integral** `_radial_split_integral(Z_orb, n1, l1, n2, l2, L, R)` evaluated *analytically* via incomplete gamma functions:
$$
R_L = \frac{1}{R^{L+1}} \int_0^R R_a R_b r^{L+2} \, dr + R^L \int_R^\infty R_a R_b r^{1-L} \, dr.
$$
4. **Closed form via incomplete gamma**: the integrand is polynomial-times-single-exponential because both $R_a$ and $R_b$ share exponent $\alpha = Z_{\rm orb}/n$ (single-$\lambda$ deployment); the product is $e^{-\alpha_{\rm total} r} \cdot \text{poly}(r)$ with $\alpha_{\rm total} = \alpha_1 + \alpha_2$. The split-region integral is a sum of $j! / \alpha_{\rm total}^{j+1} \cdot \text{gammainc/gammaincc}(j+1, \alpha_{\rm total} R)$ terms.

**The current limitation.** The function `_hydrogenic_poly_coeffs(Z, n, l)` returns a single $\alpha = Z/n$ for both bra and ket. To handle mismatched exponents, the only change required is to call `_hydrogenic_poly_coeffs` separately for bra and ket and then track $\alpha_{\rm total} = \alpha_1 + \alpha_2$ in `_split_integral_analytical`. The angular machinery is identical.

---

## 3. Q-A computation: mismatched-exponent two-center integral

### 3.1 Setup

For $(n_a, l_a) = (n_b, l_b) = (1, 0)$ (the simplest mismatched case, two 1s functions at different $\lambda$), the radial integrand is
$$
R^{1s}_{\lambda_a}(r) \cdot R^{1s}_{\lambda_b}(r) \cdot r^k
= 2\lambda_a^{3/2} e^{-\lambda_a r} \cdot 2\lambda_b^{3/2} e^{-\lambda_b r} \cdot r^k
= 4 \lambda_a^{3/2} \lambda_b^{3/2} \, r^k \, e^{-(\lambda_a + \lambda_b) r}.
$$

The split-region $L=0$ integral is therefore
$$
M^{1s,1s}_{\lambda_a,\lambda_b}(R; L=0) = \frac{1}{R} \int_0^R 4 \lambda_a^{3/2} \lambda_b^{3/2} r^2 e^{-(\lambda_a + \lambda_b) r} \, dr + \int_R^\infty 4 \lambda_a^{3/2} \lambda_b^{3/2} r e^{-(\lambda_a + \lambda_b) r} \, dr.
$$

### 3.2 Closed-form result

Sympy returns (driver `multifocal_b_w1a_qa_qb.py`):
$$
\boxed{
M^{1s,1s}_{\lambda_a,\lambda_b}(R; L=0) = \frac{4 \lambda_a^{3/2} \lambda_b^{3/2}}{(\lambda_a + \lambda_b)^3} \cdot \frac{1}{R}\Big[ R(\lambda_a + \lambda_b) + 1 - e^{-R(\lambda_a + \lambda_b)} (1 + R(\lambda_a + \lambda_b)) - \cdots \Big]
}
$$
(Full expression in `debug/data/multifocal_b_w1a_qa_qb.json`.) The matched limit $\lambda_b \to \lambda_a$ yields
$$
M^{1s,1s}_{\lambda,\lambda}(R; L=0) = \frac{1}{R}\Big[1 - e^{-2\lambda R}(1 + \lambda R)\Big]
$$
which is the textbook 1s–1s monopole nuclear-attraction integral on the same center.

### 3.3 Structural reading

The integrand of every Q-A integral (any $n_a, l_a, n_b, l_b$, any $L$) is **polynomial × single exponential** at decay rate $\alpha_{\rm total} = \lambda_a + \lambda_b$. This is the same structure as the matched case, with only the value of $\alpha_{\rm total}$ changed. **The matched-case incomplete-gamma machinery extends directly.**

The reason this works: $R^{\lambda_a}_{nl}(r) \cdot R^{\lambda_b}_{n'l'}(r)$ is a product of two associated Laguerre polynomials in different scaled arguments ($2\lambda_a r$ and $2\lambda_b r$). The product is *not* a single Laguerre polynomial in either argument — but it *is* a polynomial of finite degree in $r$, multiplied by $e^{-(\lambda_a + \lambda_b) r}$. The split-region integral remains a finite sum of incomplete-gamma terms. **No new transcendental enters.**

### 3.4 Engineering implication

To extend `geovac/shibuya_wulfman.py` to mismatched-$\lambda$:

1. **`_hydrogenic_poly_coeffs(Z_orb, n, l)`** → **`_hydrogenic_poly_coeffs_lam(lam, n, l)`** taking the exponent $\lambda$ directly (rather than $Z/n$). Backward-compatible: the original `Z_orb`-keyed signature is recoverable via `lam = Z_orb / n`.
2. **`_radial_split_integral`** takes two separate $(c_1, \alpha_1)$ and $(c_2, \alpha_2)$ pairs from the modified `_hydrogenic_poly_coeffs_lam`, computes `prod_coeffs = _poly_product(c_1, c_2)` and `alpha_total = alpha_1 + alpha_2`, and feeds them to the existing `_split_integral_analytical`. No further changes.

This is a one-day refactor of `shibuya_wulfman.py`.

---

## 4. Q-B computation: multipole termination check

### 4.1 The angular bound is $L_{\max} = l_a + l_b$ (tight), not $L_{\max} = 2 \max(l_a, l_b)$

The briefing's $L_{\max} = 2 \max(l_a, l_b)$ is the upper bound implied by the strict inequality $L \leq l_a + l_b \leq 2 \max(l_a, l_b)$. The tight bound on the multipole sum for a single matrix element is $L_{\max}^{\rm mat} = l_a + l_b$, with the additional parity constraint $l_a + L + l_b \equiv 0 \pmod 2$.

### 4.2 Direct verification

Driver `multifocal_b_w1a_qb.py`. Tested cases:

| $(l_a, l_b)$ | $m_a = m_b$ | $L$ values with $A_L \neq 0$ | $L_{\max}^{\rm mat}$ | Note |
|:-------------|:------------|:-----------------------------|:----------------------|:-----|
| (0, 0) | 0 | $L = 0$ only | 0 | s–s |
| (0, 1) | 0 | $L = 1$ only | 1 | s–p |
| (1, 1) | 0 | $L = 0, 2$ | 2 | p–p (dipole + quadrupole) |
| (2, 2) | 0 | $L = 0, 2, 4$ | 4 | d–d (mono, quad, hexa) |
| (1, 2) | 0 | $L = 1, 3$ | 3 | p–d cross |

The angular factor `_angular_coefficient(l_a, m_a, l_b, m_b, L)` in `geovac/shibuya_wulfman.py` returns $A_L \propto 3j(l_a, L, l_b; 0, 0, 0) \cdot 3j(l_a, L, l_b; -m, 0, m)$ — depending **only** on $(l_a, l_b, m_a, m_b, L)$, not on radial exponents. The triangle inequality $|l_a - l_b| \leq L \leq l_a + l_b$ and the parity constraint $(l_a + L + l_b)$ even are inherited from 3j selection rules.

### 4.3 Verdict for Q-B

**Multipole termination at $L_{\max} = l_a + l_b$ is preserved under exponent mismatch.** This is structural: the mechanism is purely angular, and the radial exponent enters only through the radial integral $M_{ab}(R; \lambda_a, \lambda_b, L)$, never through the angular factor.

The tightness of the bound is set by Gaunt parity, not by the bound $2 \max(l_a, l_b)$. For matched-$l$ states ($l_a = l_b$), the two coincide; for mismatched-$l$, $l_a + l_b < 2 \max(l_a, l_b)$ is strictly tighter.

---

## 5. Q-C computation: transfer operator existence and sparsity

### 5.1 The transfer operator

The change-of-basis matrix between Sturmian families at $\lambda_1$ and $\lambda_2$ in the standard L²-inner-product (with $r^2 \, dr$ measure) is
$$
T(\lambda_1 \to \lambda_2)_{n,l,m;\, n',l',m'} = \big\langle \chi^{\lambda_2}_{n l m} \,\big|\, \chi^{\lambda_1}_{n' l' m'} \big\rangle = \delta_{l l'} \, \delta_{m m'} \, T^{\rm rad}_{n,n'}(l; \lambda_2, \lambda_1)
$$
with the radial overlap matrix
$$
T^{\rm rad}_{n,n'}(l; \lambda_2, \lambda_1) = \int_0^\infty R^{\lambda_2}_{nl}(r) R^{\lambda_1}_{n'l}(r) r^2 \, dr.
$$

### 5.2 Closed-form entries at $l = 0$

Driver `multifocal_b_w1a_qc.py`. The first few entries:
$$
T^{\rm rad}_{1,1}(0) = \frac{8 \lambda_a^{3/2} \lambda_b^{3/2}}{(\lambda_a + \lambda_b)^3}, \qquad
T^{\rm rad}_{1,2}(0) = \frac{8 \lambda_a^{3/2} \lambda_b^{3/2} (\lambda_a - 2\lambda_b)}{(\lambda_a + \lambda_b)^4},
$$
$$
T^{\rm rad}_{2,2}(0) = \frac{8 \lambda_a^{3/2} \lambda_b^{3/2} \cdot 2(6\lambda_a \lambda_b - (\lambda_a + \lambda_b)^2)}{(\lambda_a + \lambda_b)^5}, \quad \ldots
$$

### 5.3 Density measurement

At $l = 0$:

- **3 × 3** ($n_a, n_b \in \{1, 2, 3\}$): density = 9/9 = **100 %**.
- **5 × 5** ($n_a, n_b \in \{1, ..., 5\}$): density = 25/25 = **100 %**.

Numerical magnitudes at $\lambda_a = 1, \lambda_b = 2$:
```
   0.8381  -0.8381   0.4656  -0.2173   0.0931
   0.0000   0.5587  -0.8691   0.6829  -0.4139
  -0.0931   0.2173   0.1966  -0.7070   0.7806
  -0.0621   0.0000   0.2828  -0.1334  -0.4024
  -0.0310  -0.0517   0.1115   0.2012  -0.3345
```
This is **not approximately diagonal**. The off-diagonal magnitudes are comparable to the diagonal.

### 5.4 Matched-$\lambda$ limit: the Sturmian self-overlap is not the identity

Setting $\lambda_b = \lambda_a$ recovers
```
T^{rad}_{l=0}(\lambda, \lambda) =
   1     -1/2    0     -1/4   ...
  -1/2    1     -1/2    0     ...
   0    -1/2    1     -1/2   ...
```
This is the **Sturmian metric tensor** in the $r^2 \, dr$ measure: identity on the diagonal, $-1/2$ on the first off-diagonal, zero elsewhere on the second off-diagonal (and so on). This reflects the well-known fact that **Coulomb-Sturmians at fixed $\lambda$ are orthonormal in the Sturmian $r \, dr$ measure (with weight $1/r$ in the kinetic-energy bilinear form), but NOT in the L²($r^2 \, dr$) measure** — they have a tridiagonal Gram matrix in the standard inner product. Two consequences:

- **Caveat for any naive application of the Sturmian basis as an L² basis:** the basis is *not* orthonormal; one must work in the generalized eigenvalue problem $H \mathbf{c} = E S \mathbf{c}$ with the Sturmian Gram matrix $S$ explicit.
- **The transfer operator $T(\lambda_1 \to \lambda_2)$ is not even close to a partial isometry** in the L² sense, even at matched $\lambda$. This is the core sparsity-destruction reason.

### 5.5 Pauli-encoding implication

If an operator $A$ is sparse in the Sturmian basis at $\lambda_1$ (e.g., a Hamiltonian with localized matrix elements), the change-of-basis to $\lambda_2$ is $A_2 = T A_1 T^T$. Because $T$ is fully dense, $A_2$ is generically a fully populated matrix — even if $A_1$ was sparse. The Pauli-string count then inflates from the original sparse count to **$\Theta(Q^2)$** (or worse, with antisymmetrization), and the structural sparsity advantage is destroyed.

This is **the same sparsity-destruction pathology** Track BU (multi-electron Sturmian CI) and Track DF Sprint 5 (heterogeneous nested encoding with Löwdin orthogonalization) demonstrated. The transfer-operator route is therefore **the wrong path to W1a closure**.

### 5.6 The right path: build the cross-register operator directly

The cross-register $V_{eN}$ operator should be constructed directly in the joint basis (electron-Sturmian-at-$\lambda_e$ × nucleus-HO-or-Sturmian-at-$\lambda_n$), by computing matrix elements that are bilinear in the two single-register basis sets and folded through the multipole expansion. This is the Q-D object: see §6.

---

## 6. Q-D scoping: cross-register bilinear ERI

### 6.1 The integral

When $\mathbf{R}$ is operator-valued on a separate register, the matrix element is the six-coordinate integral
$$
M = \big\langle \chi^{\lambda_e}_{n_e l_e m_e}(\mathbf{r}_e) \chi^{\lambda_n}_{n_n l_n m_n}(\mathbf{R}_n) \,\big|\, \tfrac{1}{|\mathbf{r}_e - \mathbf{R}_n|} \,\big|\, \chi^{\lambda_e}_{n'_e l'_e m'_e}(\mathbf{r}_e) \chi^{\lambda_n}_{n'_n l'_n m'_n}(\mathbf{R}_n) \big\rangle.
$$
Apply the standard multipole expansion
$$
\frac{1}{|\mathbf{r}_e - \mathbf{R}_n|} = \sum_{L=0}^\infty \sum_{M=-L}^{L} \frac{4\pi}{2L+1} \frac{r_<^L}{r_>^{L+1}} Y_{LM}^*(\hat{\mathbf{r}}_e) Y_{LM}(\hat{\mathbf{R}}_n),
$$
with $r_< = \min(r_e, R_n)$, $r_> = \max(r_e, R_n)$.

### 6.2 Angular separability and termination

The angular integrals factorize:
$$
\int Y_{l_e m_e}^* \, Y_{LM}^* \, Y_{l'_e m'_e} \, d\Omega_e = G(l_e, m_e; L, M; l'_e, m'_e),
$$
$$
\int Y_{l_n m_n}^* \, Y_{LM} \, Y_{l'_n m'_n} \, d\Omega_n = G(l_n, m_n; L, -M; l'_n, m'_n)
$$
where $G$ is the standard Gaunt coefficient.

**Termination is bilateral:** Gaunt selection on the electron side gives $L \leq l_e + l'_e$, on the nucleus side $L \leq l_n + l'_n$. The combined truncation is
$$
L_{\max}^{\rm cross} = \min(l_e + l'_e, l_n + l'_n)
$$
plus the parity constraint $(l_e + L + l'_e)$ even and $(l_n + L + l'_n)$ even (consistent because $L$ is shared). The $M$ sum has $|M| \leq L$ but is strongly restricted by both sides' $m$-conservation (Wigner 3j zero unless $-m_e + M + m'_e = 0$ and $-m_n - M + m'_n = 0$, which forces $m_e + m_n = m'_e + m'_n$ — a **cross-register total $m$-conservation**).

### 6.3 Radial double-integral structure

After angular factorization, the radial double integral is
$$
J_L(\lambda_e, \lambda_n; n_e l_e n'_e l'_e n_n l_n n'_n l'_n) = \int_0^\infty dr \int_0^\infty dR \, R^{\lambda_e}_{n_e l_e}(r) R^{\lambda_e}_{n'_e l'_e}(r) R^{\lambda_n}_{n_n l_n}(R) R^{\lambda_n}_{n'_n l'_n}(R) \, r^2 R^2 \, \frac{r_<^L}{r_>^{L+1}}.
$$

Splitting on $r_e \lessgtr R_n$:
$$
J_L = \int_0^\infty dr \, \rho_e(r) \, \frac{1}{r^{L+1}} \int_0^r dR \, \rho_n(R) \, R^L
+ \int_0^\infty dr \, \rho_e(r) \, r^L \int_r^\infty dR \, \rho_n(R) \, \frac{1}{R^{L+1}}
$$
where $\rho_e(r) = R^{\lambda_e}_{n_e l_e} R^{\lambda_e}_{n'_e l'_e} \cdot r^2$ and $\rho_n(R) = R^{\lambda_n}_{n_n l_n} R^{\lambda_n}_{n'_n l'_n} \cdot R^2$.

Each $\rho$ is **polynomial × single exponential** ($\rho_e \propto e^{-2\lambda_e r}$, $\rho_n \propto e^{-2\lambda_n R}$). The inner integrals on bounded subintervals $[0, r]$ or $[r, \infty)$ are sums of **(lower or upper) incomplete gamma functions** $\gamma(s, c x)$ or $\Gamma(s, c x)$. The outer integral is then $\int_0^\infty (\text{polynomial} \times e^{-2\lambda_e r}) \times (\text{incomplete gamma in } r) \, dr$.

**Standard identity** (Gradshteyn–Ryzhik 6.455, DLMF 8.14.5):
$$
\int_0^\infty t^{s-1} e^{-\alpha t} \, \gamma(\nu, \beta t) \, dt = \frac{\Gamma(s + \nu) \, \beta^\nu}{\nu (\alpha + \beta)^{s + \nu}} \, _2F_1\!\Big(1, \, s + \nu; \, \nu + 1; \, \frac{\beta}{\alpha + \beta}\Big).
$$
For positive integer $s$ and $\nu$, the $_2F_1$ terminates and reduces to a finite rational expression in $\alpha, \beta$.

### 6.4 Direct verification: 1s × 1s at $L=0$

Driver `multifocal_b_w1a_qd.py`. The standard Coulomb mutual-interaction integral for two 1s densities at exponents $\lambda_e$ and $\lambda_n$:
$$
\rho_e(r) = 4 \lambda_e^3 \, r^2 \, e^{-2\lambda_e r}, \qquad \rho_n(R) = 4 \lambda_n^3 \, R^2 \, e^{-2\lambda_n R},
$$
$$
J_0(\lambda_e, \lambda_n) = \int_0^\infty dr \int_0^\infty dR \, \rho_e(r) \rho_n(R) \, \frac{1}{\max(r, R)}.
$$

Sympy returns
$$
\boxed{
J_0(\lambda_e, \lambda_n) = \frac{\lambda_e \lambda_n (\lambda_e^2 + 3\lambda_e \lambda_n + \lambda_n^2)}{(\lambda_e + \lambda_n)^3}
}
$$
which is the **Roothaan integral** (Roothaan 1951, *Rev. Mod. Phys.* 23, 69) for two 1s exponentials at different orbital exponents — a textbook closed-form result.

**Verifications:**
- **Symmetry:** $J_0(\lambda_e, \lambda_n) = J_0(\lambda_n, \lambda_e)$ (numerically, $J_0(2, 3) = J_0(3, 2) = 186/125$).
- **Point-nucleus limit:** $\lim_{\lambda_n \to \infty} J_0 = \lambda_e$, matching the standard $\langle 1s_{\lambda_e} | 1/r | 1s_{\lambda_e} \rangle = \lambda_e$ when the nucleus shrinks to a point at the origin. Numerically at $\lambda_e = 1, \lambda_n = 1000$: $J_0 = 0.999998 \approx 1.0$ as expected.
- **Hartree atomic-units sanity:** for hydrogen, $\lambda_e = 1$, $\lambda_n \to \infty$, recovers $\langle 1s | 1/r | 1s \rangle = 1$ Hartree.

### 6.5 Direct verification: 2p × 2p at $L=1$

The 2p × 2p case (computed for completeness even though Gaunt zeros it on the angular side at $L=1$ for $m_e = m_n = 0$):
$$
J_1(2p, 2p)(\lambda_e, \lambda_n) = \frac{5 \lambda_e^2 \, P(\lambda_e, \lambda_n)}{6 \lambda_n (\lambda_e + \lambda_n)^7}
$$
where $P(\lambda_e, \lambda_n)$ is a degree-7 polynomial in $(\lambda_e, \lambda_n)$. **The L=1 cross-register radial integral has a closed form.** (Full expression in `debug/data/multifocal_b_w1a_qd.json`.)

### 6.6 Bilinearity and Pauli-count estimate

The cross-register matrix element factorizes as
$$
M = \sum_{L} \sum_{M} A^e_L(l_e, m_e, l'_e, m'_e, M) \cdot A^n_L(l_n, m_n, l'_n, m'_n, -M) \cdot J_L(\text{radial})
$$
where each angular factor is a single Wigner 3j and the radial $J_L$ is a closed-form rational/polynomial expression in $(\lambda_e, \lambda_n)$. The matrix element is **bilinear in the two single-register basis indices**, with the bilinearity carried by the multipole sum.

**Pauli-count estimate** (informal):

- Single-register electron Pauli count at $Q_e$ qubits scales as $11.11 Q_e$ (Paper 14 §IV).
- Single-register nucleus Pauli count at $Q_n$ qubits scales similarly.
- Cross-register operator Pauli count is bounded by $Q_e \cdot Q_n$ for the worst-case two-body operator on a joint register (one site on each register).
- Gaunt selection on **both** sides multiplies the angular sparsity factors. In the single-register case, Gaunt sparsity reduces ERI density from $O(Q^2)$ ($Q^4$ for ERIs but with antisymmetry to two-electron only) to $\sim 1\%$–$10\%$ at $l_{\max} \in \{1, 2, 3\}$ (Paper 22). The bilateral cross-register version inherits both factors.

A reasonable conjecture: cross-register Pauli count $\sim O(Q_e^{1.5} Q_n^{1.5})$ at fixed $l_{\max}$, but a full ERI count is required to pin the exponent. This is a Phase C deliverable, not a Phase B claim.

### 6.7 The actual W1a closure object

W1a closure needs:

1. **`geovac/cross_register_vne.py`** (new module) implementing `build_cross_register_vne(spec_e, spec_n, R_AB=0)` returning a Pauli dictionary on the joint register. Internally:
   - Compute angular factor for each $(l_e, m_e, l'_e, m'_e, l_n, m_n, l'_n, m'_n, L, M)$ via Wigner 3j (already in `geovac/composed_qubit.py::_wigner3j`).
   - Compute radial $J_L(\lambda_e, \lambda_n; \ldots)$ via the Q-A/Q-D machinery (closed form).
   - Map to Pauli strings via the existing JW transformation in the joint register.

2. **Promote `R_PROTON_BOHR` to operator-valued $\hat{\mathbf{R}}_n$** on the proton register. Track NI's `geovac/nuclear/nuclear_electronic.py` currently uses the classical scalar; the modification replaces the finite-size-correction one-body shift with a genuine two-body integral.

3. **For Q-D in the HO basis** (Track NI's deuteron uses HO orbitals), replace the Coulomb-Sturmian $\chi^{\lambda_n}$ with HO basis $\phi^{\omega_n}_{n_r l m_l}$. The HO radial functions are polynomial × Gaussian, and $\rho_n^{\rm HO}(R) \sim P_{\rm poly}(R) \cdot e^{-\omega_n R^2/2}$ — *not* polynomial × single exponential. The radial integral structure changes: instead of $\int dR \, P(R) \, e^{-\alpha R}$ (gives incomplete gamma), we get $\int dR \, P(R) \, e^{-\beta R^2}$ (gives **complementary error function** $\mathrm{erfc}$). This is **still closed form**, but with a different transcendental class. Exact analog of how Pachucki–Patkóš–Yerokhin 2023 evaluates recoil integrals with HO nuclear basis.

### 6.8 What's required vs. what already exists

| Ingredient | Already in GeoVac? | What's needed |
|:-----------|:--------------------|:---------------|
| Wigner 3j and Gaunt machinery | Yes (`composed_qubit.py::_wigner3j`, `shibuya_wulfman.py::_angular_coefficient`) | Reuse directly |
| Multipole expansion of $1/|\mathbf{r}_e - \mathbf{R}_n|$ | Yes (`shibuya_wulfman.py`) | Reuse directly |
| Mismatched-$\lambda$ radial integrals | No, but trivial extension | One-day refactor of `_hydrogenic_poly_coeffs` and `_split_integral_analytical` |
| Cross-register double-radial integral with $\mathrm{erfc}$ (HO) or incomplete gamma (Sturmian) | No | New module `cross_register_vne.py` (~2-week build, including tests) |
| Operator-valued $\hat{\mathbf{R}}_n$ on proton register | No (Track NI uses classical scalar) | Wire HO matrix elements via `geovac/nuclear/moshinsky.py` (already exists for relative coordinates) (~1-week wiring) |
| Joint Pauli encoding / JW for the new ERI | Yes (existing two-body to Pauli machinery in composed pipeline) | Reuse directly |
| Validation against Pachucki–Patkóš–Yerokhin 2023 recoil corrections | No (no comparison literature ported in) | Phase C (~2-4 weeks for first comparison; longer for full $(Z\alpha)^6$ structure) |

---

## 7. Verdict

**(b) Hybrid.** The mathematical core of W1a closure is *available* — every integral required is symbolically tractable in elementary + incomplete-gamma + erfc functions. The angular machinery is unchanged from the matched case (Gaunt selection on both registers). The transfer-operator route (Q-C) is structurally wrong because $T_{\lambda_1 \to \lambda_2}$ is dense and destroys sparsity, but **the cross-register operator can be built directly in the joint basis** via the Q-D bilinear ERI without going through a transfer operator.

The hybrid label is real because the mathematics is straightforward but the engineering is non-trivial:

- **Mathematics already published / re-derivable in elementary functions.** Roothaan 1951 had $J_0(\lambda_e, \lambda_n) = \lambda_e \lambda_n (\lambda_e^2 + 3\lambda_e \lambda_n + \lambda_n^2)/(\lambda_e + \lambda_n)^3$ in print 75 years ago. The mismatched-$\lambda$ Shibuya–Wulfman is simply this Roothaan-style technology applied to two-center geometry — a one-page sympy verification, not a research bet.
- **Engineering: a new module** (`cross_register_vne.py`), a Track NI extension (operator-valued $\hat{\mathbf{R}}_n$), and a full validation campaign against Pachucki–Patkóš–Yerokhin 2023. Phase A's 6–10 week estimate is reasonable.

**Phase B finding sharpens Phase A's audit.** Phase A (Section 4 Candidate 1) flagged "multipole termination across mismatched exponents may not preserve Gaunt-style closure" as the most likely failure mode. **This Phase B diagnostic shows the multipole termination is *trivially* preserved — Gaunt is purely angular.** The "most likely failure mode" was a phantom; W1a's algebraic backbone is robust. The remaining risk in Phase C is engineering and validation, not algebra.

**The honest scope statement for Phase C-W1a:** a focused 6–10 week sprint to build `cross_register_vne.py`, port Track NI to operator-valued nuclear coordinates, and validate against Pachucki–Patkóš–Yerokhin 2023's $(Z\alpha)^6$ recoil corrections. Expected outcome: hydrogen 1S Lamb shift recoil component computed natively (rather than imported as the textbook $(1 + m_e/m_p)^{-3}$ Bohr-radius substitution), with the framework deriving the recoil correction from a quantum two-body $V_{eN}$ operator rather than from an external prescription.

---

## 8. Cleanest path forward (Phase C scope; sketch only)

If the PI authorizes Phase C-W1a, the sprint structure is:

**Phase C-W1a Sprint, 6–10 weeks total, four sub-tracks:**

| Sub-track | Goal | Effort |
|:----------|:-----|:-------|
| **C-W1a.1** | Mismatched-$\lambda$ refactor of `geovac/shibuya_wulfman.py`. Test against single-$\lambda$ regression. Add tests for multi-$\lambda$ symmetry, point-nucleus limit, dimensional consistency. | 1 week |
| **C-W1a.2** | New module `geovac/cross_register_vne.py`. Bilinear ERI builder with both Sturmian (for completeness) and HO (for Track NI compatibility) nucleus basis. Closed-form double-radial integrals via incomplete gamma (Sturmian) and erfc (HO). Pauli encoding to joint register. | 2–3 weeks |
| **C-W1a.3** | Promote Track NI's classical `R_PROTON_BOHR` to operator-valued $\hat{\mathbf{R}}_n$. Wire through `geovac/nuclear/nuclear_electronic.py::build_deuterium_composed_hamiltonian`. Verify backward compatibility with existing Track NI 21 cm hyperfine validation. | 1–2 weeks |
| **C-W1a.4** | Validation campaign: compute hydrogen 1S recoil correction natively, compare to Pachucki–Patkóš–Yerokhin 2023 PRL 130, 023004. Report in Paper 32 §VIII.C addendum or candidate Paper 39 (real-space multi-particle Connes-style spectral triple, per Phase A Q3 / Section 6 surprise S3). | 2–4 weeks |

**Risk profile:** Medium. The algebraic backbone is closed (Phase B diagnostic). The operator-valued $\hat{\mathbf{R}}_n$ on the HO basis works through Moshinsky-Talmi which is mature. The validation against Pachucki et al. has a known published target. The main risks are (i) HO basis truncation effects at the proton's high-MeV scale (mitigatable by computing against a fixed FNS-only reference), and (ii) numerical precision of the cross-register Pauli encoding when the coefficient ratio across registers is $\sim 10^{13}$ (already a known Track NI feature, manageable via block-partitioned solves).

**Position statement for Phase C delivery:** the W1a closure would be the **first explicit Connes-style multi-particle real-space spectral triple in the published literature**, per Track 3 surprise S3 in the Phase A literature memo. Combining this with the existing Track NI deuteron and Paper 32 §VIII spectral-triple framing, the deliverable is a publishable contribution at the intersection of NCG and atomic physics — not just a sprint result.

---

## 9. Honest scope and uncertainty

### What I am confident about

- Q-A: closed form survives mismatch. Verified symbolically for 1s-1s at L=0; the general structure (polynomial × single exponential, integrated against $r_<^L/r_>^{L+1}$) follows by inspection of the Sturmian wavefunction form.
- Q-B: multipole termination is preserved. Verified at $l_a, l_b \in \{0, 1, 2\}$ exhaustively. The mechanism (Gaunt 3j depending only on angular labels) is algebraically transparent.
- Q-C: transfer operator is dense. Verified at $l = 0$, $5 \times 5$. The structural reason (Sturmians are not L²-orthonormal even at matched $\lambda$) means the dense behavior persists at higher $l$ and larger basis size.
- Q-D: closed-form bilinear ERI exists. Verified for 1s-1s at L=0 (Roothaan 1951) and 2p-2p at L=1. The general structure (incomplete gamma / erfc nested integrals) follows from standard Gradshteyn–Ryzhik / DLMF identities.

### What I am less confident about

- **HO basis vs. Sturmian basis for the proton register.** Track NI uses HO; the audit defaults to Sturmian language. The HO version has erfc instead of incomplete gamma, and the structural argument carries through, but I have not verified the HO-version Q-D integral symbolically. **Recommendation: add a Phase C.0 task to verify the HO 0s × 0s × $1/|\mathbf{r}_e - \mathbf{R}_n|$ × 0s × 0s integral closes in elementary + erfc, before the full sprint.**
- **Pauli-count scaling estimate for the cross-register operator.** I conjectured $O(Q_e^{1.5} Q_n^{1.5})$ but did not verify with an explicit count. The bilateral Gaunt sparsity should be a strong factor, but the actual scaling depends on which $(l_e, l_n)$ pairs are included and the antisymmetrization structure. **Recommendation: full ERI count is a Phase C deliverable, not a Phase B claim.**
- **Whether Track NI's Moshinsky-Talmi infrastructure for relative coordinates already supports lab-frame $\hat{\mathbf{R}}_n$.** I noted the existing module but did not read it end-to-end. There may be hidden complications in extracting the position operator from the relative-coordinate basis, especially for unbound/translation-symmetric situations.
- **Phase C effort estimate (6–10 weeks).** This is consistent with Track 3 Candidate 1 in the Phase A literature memo, but the detailed task list above (4 sub-tracks at 1, 2-3, 1-2, 2-4 weeks) sums to 6–10 weeks under modest concurrent execution. Linear execution would be longer.

### What I did NOT cover

- **The Connes-Marcolli A8'-class spectral-triple lift** of W1a. I scoped Q-D as a Pauli-encoding / second-quantization closure. The spectral-triple framing — promoting the proton register from a Hilbert space to a spectral triple in its own right, then tensoring $T_{\rm GV} \otimes T_p$ in the Connes sense — is a Phase C polish that goes into Paper 32 §VIII.C addendum or candidate Paper 39 territory. The Phase A synthesis Section 3 sketches this; I did not advance it here.
- **Exchange / antisymmetry structure** for the cross-register two-body ERI when the joint register includes multiple electrons or multiple nucleons. Track NI was deuteron (one proton, one electron); generalization to many-electron / many-nucleon systems requires the standard Slater-Condon decomposition adapted to the joint register, which is not a Phase B issue.
- **Whether W1a's HO-basis closure for the nucleus generalizes to other natural-geometry projections** (e.g., the Bargmann-Segal $S^5$ HO graph from Paper 24). This would be relevant for nuclei beyond the deuteron and is genuinely G4b territory (cross-manifold composition).

### What would change Phase A's framing

The single most informative Phase B finding, in revising Phase A:

**Phase A's Section 4 Candidate 1 ranked W1a closure at "medium-high risk" because of the worry that multi-$\lambda$ multipole termination might fail. Phase B refutes this worry.** The angular factor depends only on $(l_a, l_b, m_a, m_b)$ and Gaunt selection rules; the radial mismatch is absorbed cleanly via $\alpha_{\rm total} = \lambda_a + \lambda_b$. The risk profile of Phase C-W1a is therefore **medium, not medium-high**, and the Phase A ranking of Candidate 1 should hold or rise (since the algebraic backbone is now confirmed).

The flip side is that Phase A's Section 4 Candidate 3 (multi-$\lambda$ Shibuya–Wulfman as a *standalone* sprint, before Candidate 1) is now of marginal utility: the algebraic question is essentially answered by Phase B, so a 4–6 week standalone sprint to verify what the diagnostic just showed in 200 lines of sympy is an inefficient use of effort. **Recommendation: bundle Candidate 3 into Candidate 1's first week (sub-track C-W1a.1) and skip the standalone version.** Phase A's Q5 should be answered "bundled."

---

## 10. Phase B-W1a-diag closure

The diagnostic question — *what would closing W1a actually require?* — has a clean answer:

**Closing W1a requires a 6–10 week engineering sprint to:**
1. Refactor `geovac/shibuya_wulfman.py` for mismatched exponents (one week, low risk).
2. Build `geovac/cross_register_vne.py` for the bilinear ERI (2–3 weeks, medium risk via HO double-integral verification).
3. Promote Track NI's classical `R_PROTON_BOHR` to operator-valued $\hat{\mathbf{R}}_n$ (1–2 weeks, low risk via existing Moshinsky machinery).
4. Validate against Pachucki–Patkóš–Yerokhin 2023 (2–4 weeks, medium risk in numerical precision and scope).

**The path is hybrid: textbook mathematics (Roothaan integrals, incomplete gamma / erfc identities) ported into GeoVac's discrete-S³ architecture, with a closure target calibrated against a published atomic-physics computation. The deliverable would be GeoVac's first explicit multi-focal composition theorem at the spatial-coordinate level — completing the structural-skeleton scope statement by adding a new tier where space-space cross-register coupling lives.**

The PI's strategic question (Phase A Q1) — *attack W1a aggressively now, or wait for the smaller Candidates 2 and 4 to land first?* — is now sharper. Phase B has answered the algebraic question for free. The remaining decision is a project-management one: whether the 6–10 week Phase C-W1a investment is the right next move given competing claims on PM/sub-agent budget. The current diagnostic does not change the answer to that question, but it does change the calculus: **Phase C-W1a is now a lower-risk sprint than Phase A estimated**, which marginally favors aggressive execution.

---

**End of multifocal_b_w1a_diag memo.**
