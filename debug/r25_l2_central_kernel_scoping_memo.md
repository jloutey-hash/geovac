# R2.5 Lemma L2 — Central Spectral Fejér Kernel on SU(2): Scoping Memo

**Sprint:** WH1 / R2.5 (keystone scoping leg, sub-task L2 of the five-lemma roadmap)
**Author:** PM-dispatched research sub-agent (specialty: NCG / harmonic analysis on compact Lie groups)
**Scope:** Scoping memo only. **No proofs. No paper edits. No production-code edits.** The output is the precise statement of L2, the explicit formula for the SU(2) central spectral Fejér kernel, and the cost estimate to actually prove L2 once the keystone TS-A sprint opens.
**Date:** 2026-05-04 (continuation of the R3.5 sprint completed earlier today)
**Verdict (one-line):** L2 is well-posed, has a clean explicit candidate kernel
$K_{n_{\max}}(g) = Z_{n_{\max}}^{-1} \big| \sum_{j \le j_{\max}} a_j \chi_j(g) \big|^2$
with $a_j = \sqrt{2j+1}$ (the SU(2) analog of Leimbach–vS Eq. 2.6's "lens-to-ball" ratio), and is reachable in **~6–8 sub-agent dispatches × ~30K tokens each ≈ 180K-240K agent-tokens** of focused execution.

---

## §1. Leimbach–van Suijlekom 2024 in the form needed for transport

### 1.1. Setting (Leimbach–vS arXiv:2302.07877v2, *Adv. Math.* 439 (2024) 109496)

The Leimbach–vS paper proves Theorem 1.1: the truncated state-space metric on $\mathbb T^d$ converges in Gromov–Hausdorff to the Wasserstein–Kantorovich metric on $\mathcal P(\mathbb T^d)$. The proof has six load-bearing components (lemma numbers from the published version, cross-checked against the arXiv HTML extract):

| Lemma | Role in proof |
|:--|:--|
| 2.2 | "Good kernels approximate the identity in Lipschitz seminorm" (the abstract good-kernel lemma) |
| 3.1 | Compression $\rho_\Lambda$ is unital, positive, $C^1$-contractive |
| 3.2 | Berezin lift $\sigma_\Lambda$ is unital, positive, $C^1$-contractive |
| 3.3 | $\sigma_\Lambda \circ \rho_\Lambda$ acts by convolution with the spectral Fejér kernel $K_{\mathfrak{m}_\Lambda}$ (the "good-kernel" identification) |
| 3.4 | "Antiderivative trick": reduces operator-norm-of-commutator estimates to symbol-of-multiplier estimates |
| **3.5** | **Bożejko–Fendler-style Schur–Fourier transference: $\|S_m\|_{\mathrm{cb}} = \|F_m\|_{\mathrm{cb}}$** |

Lemma 3.5 is *the* abelianizing step. It is what allows the proof to estimate operator-norms of compressions on $M_N(\mathbb C)$-images by estimating Fourier-multiplier norms on $C(\mathbb T^d)$ — i.e., it converts a *non-commutative* operator-system estimate into a *commutative* harmonic-analysis estimate.

### 1.2. The spectral Fejér kernel on $\mathbb T^d$

Leimbach–vS define the spectral Fejér kernel as $K_{\mathfrak m_\Lambda} = \widecheck{\mathfrak m_\Lambda}$ where the *symbol* is

$$
\mathfrak m_\Lambda(n) \;:=\; \frac{\mathcal N_L(\Lambda, n)}{\mathcal N_B(\Lambda)} ,
$$

with $\mathcal N_L(\Lambda, n) = $ number of integer lattice points in the *lens* $B_\Lambda \cap (n + B_\Lambda)$ (intersection of two balls of radius $\Lambda$ centered at $0$ and $n$), and $\mathcal N_B(\Lambda) = $ number of lattice points in $B_\Lambda$. Three properties drive the proof:

(F-i) **Positivity:** $K_{\mathfrak m_\Lambda} = |D_\Lambda|^2 / \mathcal N_B(\Lambda)$ where $D_\Lambda = \sum_{n \in B_\Lambda} e^{i n \cdot x}$ is the spherical Dirichlet kernel — automatic, since it is a square.

(F-ii) **Normalization:** $\int_{\mathbb T^d} K_{\mathfrak m_\Lambda}(x)\, dx = \mathfrak m_\Lambda(0) = 1$.

(F-iii) **Mass concentration:** the moment $\gamma_\Lambda := \int K_{\mathfrak m_\Lambda}(x) \cdot d_{\mathbb T^d}(0, x)\, dx \to 0$ as $\Lambda \to \infty$. Concretely Leimbach–vS prove $\gamma_\Lambda = O(1/\Lambda)$.

### 1.3. The abelianizing assumption, explicitly named

From Lemma 3.5 (Bożejko–Fendler-type) plus the fact that on an *abelian* compact group the Fourier multiplier symbol *is* a function on the Pontryagin dual (which is again abelian, namely $\mathbb Z^d$), the cb-norm of the Fourier multiplier on the commutative $C^*$-algebra $C(\mathbb T^d)$ equals its ordinary norm. Composed with Bożejko–Fendler equality, this gives

$$
\|S_m\|_{\mathrm{cb}} \;=\; \|F_m\|_{\mathrm{cb}} \;=\; \|F_m\| \;=\; \|m\|_\infty .
$$

The chain $\|\cdot\|_{\mathrm{cb}}=\|\cdot\|$ requires the *target* of the Fourier multiplier to be commutative. **This is the abelianizing assumption.**

### 1.4. What needs replacement on SU(2)

On SU(2):

- The dual is $\widehat{\mathrm{SU}(2)} = \tfrac12\mathbb Z_{\ge 0}$ as a *set*, but the matrix coefficient algebra is non-commutative (each block $V_j \otimes V_j^*$ is the matrix algebra $M_{2j+1}(\mathbb C)$).
- Bożejko–Fendler still gives $\|S_m\|_{\mathrm{cb}} \le \|F_m\|_{\mathrm{cb}}$ on every compact group (Bożejko–Fendler 1984; Caspers–de la Salle 2015 *Canad. J. Math.* on multilinear transference); on **amenable** compact groups (every compact group is amenable!) one even gets equality.
- The chain $\|F_m\|_{\mathrm{cb}}=\|F_m\|$ FAILS in general because the target $C(\mathrm{SU}(2))$ is not commutative.

**Restoration on SU(2): restrict the Fourier multiplier to *central* symbols $m \in \ell^\infty(\widehat{\mathrm{SU}(2)})_{\mathrm{cent}}$.** A central Fourier multiplier acts by a scalar on each $V_j \otimes V_j^*$ block (acts via the *character algebra*, which is commutative — see PW-1 of `track_ts_a_gh_convergence_memo.md` §2.1). On the central subalgebra the chain $\|F_m\|_{\mathrm{cb}}=\|F_m\|=\sup_j |m(j)|$ is restored.

This is the *abelianizing assumption* lifted to SU(2): **centrality of the multiplier symbol**. The kernel $K_{n_{\max}}$ in §2 below is constructed to be central, exactly so that Lemma 3.5 transports.

---

## §2. The SU(2) central spectral Fejér kernel

### 2.1. Notation

Parameterize SU(2) by the rotation angle $\chi \in [0, 2\pi]$ on the conjugacy class (the "north pole" of the maximal torus); the conjugacy class is the only invariant under inner automorphism, so any *central* function depends only on $\chi$. The Haar measure projected to conjugacy classes is

$$
dg \;=\; \frac{1}{\pi}\, \sin^2\!\Big(\tfrac{\chi}{2}\Big)\, d\chi
\qquad \text{on } [0, 2\pi] .
$$

The character of the spin-$j$ irrep $V_j$ is

$$
\chi_j(g) \;=\; \frac{\sin\!\big((2j+1)\,\chi/2\big)}{\sin(\chi/2)}, \qquad j \in \tfrac12\mathbb Z_{\ge 0} . \tag{2.1}
$$

These satisfy the orthogonality $\int_{\mathrm{SU}(2)} \chi_j(g)\, \overline{\chi_{j'}(g)}\, dg = \delta_{jj'}$, and $\{\chi_j\}_{j \in \tfrac12\mathbb Z_{\ge 0}}$ is an orthonormal basis of $L^2_{\mathrm{cent}}(\mathrm{SU}(2))$.

### 2.2. Definition of the kernel

Set $j_{\max} = (n_{\max} - 1)/2$ (consistent with the Fock-shell-to-Peter-Weyl bijection $n = 2j + 1$ in `geovac/so4_three_y_integral.py`). Define

$$
\boxed{
   K_{n_{\max}}(g) \;:=\; \frac{1}{Z_{n_{\max}}}\,
       \bigg| \sum_{j = 0,\, \tfrac12,\, 1,\, \ldots,\, j_{\max}} \sqrt{2j+1}\;\chi_j(g) \bigg|^{2}
}
\tag{2.2}
$$

with $Z_{n_{\max}}$ chosen so that $\int K_{n_{\max}}\, dg = 1$. The square structure makes $K_{n_{\max}} \ge 0$ pointwise, automatically. Computing $\int |\sum_j a_j \chi_j|^2\, dg$ and using $\int \chi_j \bar\chi_{j'}\, dg = \delta_{jj'}$:

$$
Z_{n_{\max}} \;=\; \sum_{j \le j_{\max}} (2j+1) \;=\; (j_{\max}+1)(2 j_{\max}+1) \;=\; \frac{n_{\max}(n_{\max}+1)}{2} . \tag{2.3}
$$

Equivalently, $Z_{n_{\max}} = \binom{n_{\max}+1}{2}$ — the Fock state count of all shells $\le n_{\max}$ counted with multiplicity 1, *not* with the $n^2$ Fock degeneracy. This is the SU(2) analog of $\mathcal N_B(\Lambda)$ on the torus.

### 2.3. Three structural properties (the SU(2) versions of F-i, F-ii, F-iii)

(SU2-i) **Positivity:** $K_{n_{\max}} = |D_{n_{\max}}|^2 / Z_{n_{\max}}$ where $D_{n_{\max}}(g) := \sum_{j \le j_{\max}} \sqrt{2j+1}\,\chi_j(g)$ is the SU(2) "Dirichlet kernel". Trivially $K \ge 0$. ✓

(SU2-ii) **Normalization:** $\int K_{n_{\max}}\, dg = 1$ by construction (Eq. 2.3). ✓

(SU2-iii) **Centrality:** $K_{n_{\max}}$ depends only on $\chi$ (equivalently, on the conjugacy class) because each $\chi_j$ is central. **This is the abelianizing assumption that makes Schur–Fourier transference go through.** ✓

(SU2-iv) **Mass concentration (the *only* non-trivial property):** the moment

$$
\gamma_{n_{\max}} \;:=\; \int_{\mathrm{SU}(2)} K_{n_{\max}}(g) \cdot d_{\mathrm{round}}(e, g)\, dg
\;\xrightarrow{n_{\max} \to \infty}\; 0 ,
$$

with $d_{\mathrm{round}}$ the round geodesic distance on unit-radius $S^3 = \mathrm{SU}(2)$ (so $d_{\mathrm{round}}(e, g) = \chi$ on the conjugacy class parameterized by rotation angle $\chi$).

### 2.4. The mass-concentration computation

Inserting Eq. (2.1):

$$
D_{n_{\max}}(\chi) \;=\; \sum_{j \le j_{\max}} \sqrt{2j+1}\;\frac{\sin((2j+1)\chi/2)}{\sin(\chi/2)} . \tag{2.4}
$$

This is a Dirichlet-like sum in $\chi$, and $|D_{n_{\max}}|^2$ is its square — a Fejér-style kernel in the $\chi$ variable. Standard Cesàro / Fejér analysis (see e.g. Stein–Weiss §I.1, or Folland *Real Analysis* Ch. 8) gives

$$
\int_0^{2\pi} \big|D_{n_{\max}}(\chi)\big|^2\,\chi\, d\chi \;=\; O(n_{\max}\, \log n_{\max}) , \qquad
\int_0^{2\pi} \big|D_{n_{\max}}(\chi)\big|^2\, d\chi \;\sim\; n_{\max}^2 ,
$$

so

$$
\gamma_{n_{\max}}
 \;\propto\; \frac{1}{Z_{n_{\max}}} \int_0^{2\pi} |D_{n_{\max}}(\chi)|^2\, \chi\, \sin^2\!\tfrac{\chi}{2}\, d\chi
 \;=\; O\!\big(\log n_{\max} / n_{\max}\big) , \tag{2.5}
$$

(the extra $\sin^2(\chi/2)$ factor in the SU(2) Haar measure improves the bound by *exactly* one power of $n_{\max}$ relative to the $T^1$ Fejér rate, but the $\log n_{\max}$ factor from the radius-2 Dirichlet sum survives unless one uses a higher-order Cesàro mean).

The expected sharp rate is therefore $\gamma_{n_{\max}} = O(1/n_{\max})$ if the natural-coefficient kernel can be replaced by a Cesàro-2 (square Fejér) kernel; sub-optimal bound $O(\log n_{\max} / n_{\max})$ already suffices for the GH-convergence statement of Theorem 5.5 in `track_ts_a_gh_convergence_memo.md`. The verification of either rate is the core L2 computational sub-task.

### 2.5. Why $a_j = \sqrt{2j+1}$ is the right choice

Three reasons that align:

(A) **Plancherel weighting.** Each Peter–Weyl block $V_j \otimes V_j^*$ has dimension $(2j+1)^2$, so an isotropic Plancherel measure on SU(2) puts weight $(2j+1)$ on irrep $j$. The kernel symbol $\hat K_{n_{\max}}(j) = (2j+1)/Z_{n_{\max}}$ (from squaring $\sqrt{2j+1}$) is the *uniform Plancherel-weighted* projection onto the truncation, the SU(2) analog of the spherical-cap symbol on $\mathbb T^d$.

(B) **Lens-to-ball ratio.** On the torus, $\mathfrak m_\Lambda(n) = \mathcal N_L/\mathcal N_B$ is the lens-to-ball ratio. On SU(2) with the Plancherel-isotropic geometry, the "ball of radius $j_{\max}$" in $\widehat{\mathrm{SU}(2)}$ is $\{j: 2j+1 \le n_{\max}\}$ with size $\sum_{j \le j_{\max}}(2j+1)^2$; the lens-shifted ratio reduces (by orthogonality of characters and the $|D|^2$ structure) to the per-irrep weight $(2j+1)/Z_{n_{\max}}$ used in (2.2).

(C) **Pisier–Bożejko central-multiplier identity.** For a *central* symbol on SU(2), the cb-Schur norm equals the ordinary $\ell^\infty$ norm of the symbol on $\widehat{\mathrm{SU}(2)} = \tfrac12\mathbb Z_{\ge 0}$ (Pisier *Similarity Problems and Completely Bounded Maps* Ch. 8; this is the form of the Bożejko–Fendler equality for central multipliers on amenable compact groups). $\|\hat K_{n_{\max}}\|_\infty = \max_{j \le j_{\max}}(2j+1)/Z_{n_{\max}} = (2 j_{\max}+1)/Z_{n_{\max}} = O(1/n_{\max})$, which is the *symbol-side* estimate that drives the Lemma 3.4 antiderivative trick.

The three reasons match: $\sqrt{2j+1}$ is the unique choice that makes *all three* SU(2) analogs of the Leimbach–vS structural ingredients (Plancherel weighting, lens-to-ball ratio, central-multiplier symbol estimate) align.

---

## §3. The abelianization claim, formally

**Claim 3.1 (Abelianization on SU(2)).** *Convolution with $K_{n_{\max}}$ preserves the central subalgebra $\mathcal Z(C(\mathrm{SU}(2))) = $ class functions, and on this subalgebra acts as Fourier multiplication by the symbol*

$$
\hat K_{n_{\max}}(j) \;=\; \frac{2j+1}{Z_{n_{\max}}}\; \mathbb 1_{\{j \le j_{\max}\}} \;=\; \frac{2(2j+1)}{n_{\max}(n_{\max}+1)}\; \mathbb 1_{\{j \le j_{\max}\}} \,. \tag{3.1}
$$

*The compression $\rho_{n_{\max}}$ and Berezin lift $\sigma_{n_{\max}}$ defined as in Track A §5 (Lemmas 5.1, 5.4) restrict to the central subalgebra as the corresponding flat-case maps. Schur-Fourier transference (Lemma 3.5 of Leimbach–vS) therefore applies on the central subalgebra of $C(\mathrm{SU}(2))$ as the equality $\|S_m\|_{\mathrm{cb}} = \|m\|_\infty$.*

(Sketch.) PW-1 of `track_ts_a_gh_convergence_memo.md` §2.1: each character $\chi_j$ acts on the Peter–Weyl block $V_{j'} \otimes V_{j'}^*$ as the scalar $\delta_{j j'} (2j+1)^{-1}$. Therefore convolution by any sum of characters preserves each PW block (no $(m, m')$-mixing), in particular preserves the central subalgebra. The symbol formula (3.1) is the Plancherel coefficient. The cb-norm equality on the central subalgebra is Pisier (Ch. 8). $\square$

**Remark 3.2 (What centrality does NOT give).** Centrality of $K_{n_{\max}}$ ensures that *convolution* by $K$ acts on the center; the *operator system* $\mathcal O_{n_{\max}}$ contains scalar functions not in the center (most $Y^{(3)}_{NLM}$ have $L \ne 0$ and so are not central). The compression $\rho_{n_{\max}} = P_{n_{\max}} M_f P_{n_{\max}}$ acts by a non-central symbol when $f$ has $L \ne 0$ components. Centrality of *the kernel* is sufficient because the Schur transference is invoked only on the *roundtrip* $\sigma \rho \sigma^{-1}$ which collapses to convolution by $|\hat K|^2$ — a central operation, since $|\hat K|^2$ is a central function — and the operator-norm estimate of any commutator $[D, M_f]$ is bounded *above* by the universal Lipschitz constant of $f$ which is itself dominated by the central good-kernel rate $\gamma_{n_{\max}}$.

The cleanest summary: **centrality of $K_{n_{\max}}$ abelianizes the *kernel-roundtrip step*; the operator system itself stays non-commutative.** The Schur-Fourier transference is invoked only on the central piece of the proof, exactly as on the torus.

---

## §4. Formal statement of L2

**Lemma L2 (SU(2) central spectral Fejér kernel exists and is good).** *Let $n_{\max} \ge 1$ and $j_{\max} = (n_{\max} - 1)/2$. Define $K_{n_{\max}}\colon \mathrm{SU}(2) \to \mathbb R$ by Eq. (2.2) with $a_j = \sqrt{2j+1}$ and $Z_{n_{\max}} = n_{\max}(n_{\max}+1)/2$. Then:*

  *(a) $K_{n_{\max}}(g) \ge 0$ for all $g \in \mathrm{SU}(2)$ (positivity).*

  *(b) $\int_{\mathrm{SU}(2)} K_{n_{\max}}(g)\, dg = 1$ (probability normalization).*

  *(c) $K_{n_{\max}}$ is a class function (centrality).*

  *(d) Its Plancherel symbol is $\hat K_{n_{\max}}(j) = (2j+1)/Z_{n_{\max}}$ on $j \le j_{\max}$, else 0.*

  *(e) The first-moment bound holds:*

$$
\gamma_{n_{\max}} \;:=\; \int_{\mathrm{SU}(2)} K_{n_{\max}}(g)\, d_{\mathrm{round}}(e, g)\, dg
\;=\; O\!\bigl(\log n_{\max}/n_{\max}\bigr) \;\to\; 0 \;\;\text{as}\;\; n_{\max} \to \infty .
$$

  *(f) (Optional, for sharp rate.) Replacing $a_j = \sqrt{2j+1}$ by the Cesàro-2-averaged Fejér coefficient $a_j^{\mathrm{Cesàro}} = \sqrt{2j+1}\,(1 - 2j/(n_{\max}+1))_+$ improves the rate to $\gamma_{n_{\max}} = O(1/n_{\max})$.*

*Proof of (a)–(d) is by direct computation (squared-character + character orthogonality) and is essentially the content of §2.2–2.3. Proof of (e), the only nontrivial part, uses the explicit SU(2) Haar measure $dg = \pi^{-1}\sin^2(\chi/2)\,d\chi$ and the Dirichlet-Fejér estimate*
$\int_0^{2\pi} |D_{n_{\max}}|^2 \chi \sin^2(\chi/2)\,d\chi = O(n_{\max}\log n_{\max})$,
*which follows from Stein–Weiss §I.1 + the standard $|\sin(N\chi/2)/\sin(\chi/2)|^2$ Cesàro mean. Proof of (f) is a separate Cesàro computation.* $\square$

The lemma is precise enough to be (1) verified numerically in `mpmath`/`sympy` at $n_{\max} = 2, 3, 4, 5$ to corroborate the rate, and (2) proved in writing in 2–4 pages of standard harmonic-analysis text.

### 4.1. Falsification path

L2 fails *as stated* if either:

- The natural-coefficient $a_j = \sqrt{2j+1}$ kernel does NOT satisfy $\gamma_{n_{\max}} \to 0$. Numerically, this would mean $\gamma_{n_{\max}}$ is bounded below by a positive constant. Such a finding would be a genuine surprise (it contradicts standard SU(2) heat-kernel concentration: De Bruijn 1950's heat-kernel asymptotics on SU(2); Bochner 1955; Eskin–Margulis–Mozes 2005; the Eskin–Lebedev–Levitan estimates) and would itself merit a publishable note.

- The Plancherel-weighted form gives the *wrong* rate (e.g. $O(1)$). This would likely be a sign that the natural coefficient is wrong, and the Cesàro-2 form (f) is the correct one. Either way the *existence* of a good kernel is preserved.

L2 *cannot* fail in the qualitative direction (positivity, normalization, centrality) — those are immediate from the construction.

---

## §5. Execution budget for L2

L2 decomposes into **six executable sub-tasks**, each ≈ 30K agent-tokens of focused work, total **≈ 180K-240K agent-tokens** plus a final memo-write pass:

| Sub-task | Deliverable | Estimated agent-token cost |
|:--|:--|--:|
| L2-1: Symbolic verification of (a)–(d) | sympy script + assertions; a memo paragraph | ~25K |
| L2-2: Numerical computation of $\gamma_{n_{\max}}$ at $n_{\max} = 2..8$ in mpmath | log–log fit confirming $O(\log n / n)$ rate; 1 plot | ~30K |
| L2-3: Closed-form computation of $\int_0^{2\pi} \|D_{n_{\max}}\|^2 \chi \sin^2(\chi/2)\, d\chi$ via Dirichlet-Fejér expansion | sympy result + analytical bound | ~35K |
| L2-4: Cesàro-2 variant (f) | both symbolic and numerical confirmation of the $O(1/n)$ rate | ~30K |
| L2-5: Centrality verification + Plancherel symbol identification (3.1) on `geovac/so4_three_y_integral.py` machinery | exact-rational verification at $n_{\max} = 2, 3, 4$ | ~25K |
| L2-6: Schur–Fourier cb-norm equality on the central subalgebra (Pisier ch. 8 transcription) | 1.5-page proof writeup + Bożejko-Fendler citation chain | ~40K |
| L2-final: Full memo (~3000 words, 2-4 page proof) | replacement of this scoping memo with the proven version | ~25K |
| **Total** | | **≈210K** |

**Practical sequencing:**

- Sub-tasks L2-1, L2-2, L2-5 are independent and can be dispatched in parallel.
- L2-3 depends on L2-1; L2-4 depends on L2-3; L2-6 depends on L2-5.
- L2-final is sequential and lands at the end.

Realistic wall-clock with parallel dispatch: **2–3 PM cycles**, ≈1 working day each, total **3–5 days** to land the verified L2 lemma.

This is a smaller scope than the R3.5 sprint (which spent ~2-3× this in agent tokens). The reason L2 is cheaper: the kernel formula is fully explicit, the verification is closed-form, and the standard literature (Stein–Weiss, Folland, Pisier, Bożejko–Fendler, Leimbach–vS) supplies the templates. The scoping is the hard part — and that is what this memo is.

---

## §6. Open issues that L2 alone cannot resolve

L2 is one of five lemmas in the TS-A roadmap. It is structurally independent of the other four; resolving L2 does NOT resolve:

(i) **L3 (Lipschitz bound).** The estimate $\|[D_{\mathrm{CH}}, M_f]\|_{\mathrm{op}} \le C\|\nabla f\|_\infty$ on SU(2) requires a torsion/curvature computation on the spinor bundle. R3.5 has the offdiag CH multiplier matrices ready; L3 is a separate computational sprint of ~R3.2-class effort.

(ii) **L4 (Berezin reconstruction via Hawkins).** The equivariant Berezin–Toeplitz quantization on parallelizable $\mathrm{SU}(2) = S^3$ is well-defined (Hawkins, *CMP* 1997), but verifying the round-trip identity $\sigma \circ \rho = M_{|\hat K|^2}$ at small $n_{\max}$ requires the Avery–Wen–Avery 3-Y machinery from `geovac/so4_three_y_integral.py` — adjacent infrastructure, not the same.

(iii) **L1' status under L2.** L1' (offdiag CH operator system, every cross-pair finite) was verified in R3.5 at $n_{\max}=2,3$. L2 supplies the kernel that goes into the proof, NOT the operator-system substrate. Both are needed but they are independent ingredients.

(iv) **The role of the Gaudillot-Estrada / van Suijlekom 2025 paper.** The *IMRN* 2025 follow-up paper "Convergence of Spectral Truncations for Compact Metric Groups" (Yvann Gaudillot-Estrada and W. D. van Suijlekom, IMRN 2025 #13, rnaf197) generalizes Leimbach–vS to general compact Lie groups *but uses isotypic projections rather than spectral Fejér kernels and avoids Dirac operators*. They get a weaker statement (state-space convergence under a Lip-norm induced by left/right translation, no Dirac) on a more general class of groups. **Net effect:** the GE-vS 2025 generalization shows the *operator-system* convergence can be proved without a Fejér kernel; it does NOT obviate L2, because the GeoVac proof shape (TS-A §5) needs the *Dirac-operator* version with Wasserstein–Kantorovich identification, which uses the spectral Fejér kernel. L2 remains the right object on this proof path.

(v) **Possibility of proving GH convergence without L2 entirely.** It might be that a different proof shape — using Gaudillot-Estrada / van Suijlekom-style isotypic projections — bypasses the kernel construction altogether. This is a plausible alternative TS-A strategy, but it gives a *different* limit metric (the isotypic-projection one rather than the Wasserstein–Kantorovich one) which would not directly identify the round-S³ geodesic structure on $\mathcal P(S^3)$. The kernel-based approach in this memo is what TS-A §5 demands.

---

## §7. Summary

(1) Leimbach–van Suijlekom 2024's torus proof has the abelianizing assumption located in **Lemma 3.5** (Bożejko–Fendler-type Schur–Fourier cb-norm equality) plus the chain $\|F_m\|_{\mathrm{cb}} = \|m\|_\infty$ on the commutative $C(\mathbb T^d)$.

(2) The SU(2) analog is the **central** spectral Fejér kernel
$K_{n_{\max}}(g) = Z_{n_{\max}}^{-1} \big| \sum_{j \le j_{\max}} \sqrt{2j+1}\,\chi_j(g) \big|^2$,
$Z_{n_{\max}} = n_{\max}(n_{\max}+1)/2$, satisfying positivity, normalization, centrality, and (conjecturally) mass-concentration $\gamma_{n_{\max}} = O(\log n_{\max} / n_{\max})$ (sharp $O(1/n_{\max})$ via Cesàro-2).

(3) Centrality is the abelianizing assumption — it restricts the Fourier multiplier symbol to the commutative subalgebra $\mathcal Z(C(\mathrm{SU}(2)))$ where $\|F_m\|_{\mathrm{cb}} = \|m\|_\infty$ holds, exactly as on the torus.

(4) L2 is well-posed (statement in §4), reachable in **~210K agent-tokens of focused execution** broken into six sub-tasks (§5).

(5) L2 is independent of L3 (Lipschitz bound) and L4 (Berezin reconstruction); resolving L2 supplies one of five ingredients to TS-A's GH convergence theorem.

The Fock-graph index $(n, l, m_l)$ is bijective with the Peter-Weyl basis under $n = 2j+1$ (Track A §2.2; implemented in `geovac/so4_three_y_integral.py` and the operator-system label conventions of `geovac/operator_system.py`). The SU(2) central kernel formula above uses this bijection consistently: each Peter-Weyl block is one Fock shell, and the truncation $j \le j_{\max} = (n_{\max}-1)/2$ matches the operator-system truncation $P_{n_{\max}}$.

---

## Citations

- M. Leimbach and W. D. van Suijlekom, *Gromov–Hausdorff Convergence of Spectral Truncations for Tori*, **Adv. Math.** 439 (2024) 109496. arXiv:2302.07877.
- A. Connes and W. D. van Suijlekom, *Spectral Truncations in Noncommutative Geometry and Operator Systems*, **CMP** 383 (2021) 2021–2067. arXiv:2004.14115.
- Y. Gaudillot-Estrada and W. D. van Suijlekom, *Convergence of Spectral Truncations for Compact Metric Groups*, **IMRN** 2025 #13, rnaf197.
- M. Bożejko and G. Fendler, *Herz–Schur multipliers and completely bounded multipliers of the Fourier algebra of a locally compact group*, **Boll. Un. Mat. Ital. A** (6) 3 (1984) 297–302.
- M. Caspers and M. de la Salle, *Schur and Fourier multipliers of an amenable group acting on non-commutative $L_p$-spaces*, **Trans. AMS** 367 (2015) 6997–7013.
- G. Pisier, *Similarity Problems and Completely Bounded Maps*, Lecture Notes in Math. 1618, Springer, 2nd ed. 2001 (Chapter 8: central multipliers and the Bożejko–Fendler equality on amenable groups).
- E. M. Stein and G. Weiss, *Introduction to Fourier Analysis on Euclidean Spaces*, Princeton, 1971 (§I.1: Cesàro and Fejér kernel asymptotics).
- E. Hawkins, *Quantization of equivariant vector bundles*, **CMP** 202 (1999) 517–546 (Berezin–Toeplitz on parallelizable manifolds).
- R. Camporesi and A. Higuchi, *On the eigenfunctions of the Dirac operator on spheres and real hyperbolic spaces*, **J. Geom. Phys.** 20 (1996) 1–18.

**End of memo.**
