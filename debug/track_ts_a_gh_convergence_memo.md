# Track TS-A — Sketch of Gromov–Hausdorff Convergence on S³ via Peter–Weyl on SU(2)

**Sprint:** WH1 / R2.5 (keystone scoping memo, *not* a sprint deliverable)
**Author:** PM-dispatched research sub-agent
**Scope:** Research memo only. No paper edits. No production-code edits. The task is to lay out the proof *shape* for GH convergence of the Connes–van Suijlekom (CvS) state-space metric on the GeoVac truncated operator system to the round-S³ Wasserstein–Kantorovich metric, identify the lemma chain, and stress-test it against the published proofs on flat structures (S¹: Hekkelman 2022; T^d: Leimbach–van Suijlekom 2024 [Adv. Math. 439, 109496]; discrete polynomial-growth: arXiv:2309.13469; UCP-maps: arXiv:2410.15454; quantum spheres: Indian J. Pure Appl. Math. 2024).
**Date:** 2026-05-04
**Verdict (one-line):** **(b) reachable, but harder than the flat case** — the proof shape transports cleanly through Peter–Weyl on SU(2), and the spectral-Fejér / good-kernel mechanism has a non-abelian analog via *positive central trigonometric polynomials in characters*; however three obstructions visible already in the round-3 / R3.1 / R3.2 numerics need to be threaded carefully. Time estimate: a focused two-track sprint (one mathematician-track, one computational-track) on the order of 4–8 weeks could get a Theorem-+-Proof manuscript to first-draft state. The output is a J. Geom. Phys. or Adv. Math. companion to Leimbach–vS 2024; a non-trivial NCG result worth the effort independently of WH1.

---

## §1. The Connes–vS framework recalled, in the form needed for S³

### 1.1. Operator-system spectral truncation

Following Connes–vS 2021 (arXiv:2004.14115v2, lines 1058–1101 of the local PDF
extract): given a spectral triple $(\mathcal A, \mathcal H, D)$ on a compact
Riemannian spin manifold $M$, fix a self-adjoint orthogonal projection
$P : \mathcal H \to \mathcal H$ that commutes with $D$. The spectral truncation
is then the *operator-system spectral triple*
$(P\mathcal A P,\, P\mathcal H,\, P D P)$, where $P\mathcal A P$ is a $*$-closed
(but in general *not* multiplicatively closed) self-adjoint subspace of
$\mathcal B(P\mathcal H)$ — an operator system in the sense of Choi–Effros. The
state space $\mathcal S(P\mathcal A P)$ is a compact convex set in the weak-$*$
topology; the operator-system Connes distance (Connes–vS Eq. (2)) is

$$
   d_D(\varphi, \psi) \;:=\; \sup\big\{\, |\varphi(x) - \psi(x)| \,:\, x = x^{*} \in P\mathcal A P,\ \|[D, x]\|_{\mathrm{op}} \le 1 \,\big\}. \tag{1}
$$

When $P = \mathbb 1$ this reduces to Connes' formula on the full $C^{*}$-algebra
(Connes 1989/1994).

### 1.2. The GeoVac / S³ instance

GeoVac fixes:

- $M = S^3$, equipped with the round metric of radius 1 and the canonical spin
  structure (Camporesi–Higuchi 1996, KO-dim $= 3$ mod $8$).
- $\mathcal A = C^\infty(S^3)$ acting by pointwise multiplication.
- $\mathcal H = L^2(S^3)$ for the scalar layer, lifted to
  $L^2(S^3, \Sigma) = L^2(S^3) \otimes \mathbb C^2$ for the Weyl spinor sector.
- $D = $ Camporesi–Higuchi Dirac operator with eigenvalues
  $\pm(n_{ch} + 3/2)$, equivalently $\pm(n_{\mathrm{Fock}} + 1/2)$, and
  multiplicity $g_n^{\mathrm D} = 2(n+1)(n+2)$ on each chirality.
- $P_{n_{\max}}$ = orthogonal projection onto the span of the first
  $n_{\max}$ Fock shells, i.e. eigenvectors with $1 \le n \le n_{\max}$. This is
  an *index truncation*, not a Borel functional-calculus cutoff
  $\chi_{[-\Lambda,\Lambda]}(D)$. Correction confirmed in `wh1_connes_vs_pdf_verification.md` against Connes–vS line 1113.

The truncated operator system at level $n_{\max}$ is therefore

$$
   \mathcal O_{n_{\max}} \;=\; P_{n_{\max}} \, C^\infty(S^3) \, P_{n_{\max}}
   \;\subset\; M_{N(n_{\max})}(\mathbb C),
$$

with $N(n_{\max}) = \sum_{n=1}^{n_{\max}} n^2 = 1, 5, 14, 30, 55, 91, \ldots$,
constructed explicitly in `geovac/operator_system.py` (R2.1) and lifted to
the Weyl sector in `geovac/spinor_operator_system.py` (R3.2).

### 1.3. The two ingredients that need to converge

The GH-convergence statement we want has two pieces:

(GH-a) The state space $(\mathcal S(\mathcal O_{n_{\max}}), d_{D_{n_{\max}}})$
       converges in Gromov–Hausdorff distance as $n_{\max} \to \infty$ to a
       limit metric space.

(GH-b) The limit metric space *is* $(\mathcal P(S^3), d_{\mathrm{Wass}})$, the
       Borel probability measures on the round $S^3$ with the
       Monge–Kantorovich (= Wasserstein-1) distance computed against the round
       geodesic distance on $S^3$.

This is the precise non-abelian analog of Theorem 1.1 of Leimbach–van Suijlekom
2024 (arXiv:2302.07877v2) for $\mathbb T^d$.

---

## §2. The Peter–Weyl basis on SU(2) and the Fock-graph bijection

### 2.1. Peter–Weyl, written out

Identify $S^3 = \mathrm{SU}(2)$ via the standard parameterization
$g(\chi, \theta, \phi) = \exp(i\,\chi\,\hat n(\theta, \phi)\!\cdot\!\vec\sigma)$.
The Peter–Weyl theorem (Peter & Weyl 1927; e.g. Sugiura, Knapp) gives

$$
   L^2(\mathrm{SU}(2)) \;=\; \widehat\bigoplus_{j \in \tfrac12 \mathbb Z_{\ge 0}} \pi_j \otimes \pi_j^{*}
   \;\cong\; \bigoplus_{j \in \tfrac12\mathbb Z_{\ge 0}} V_j \otimes V_j^{*}, \tag{2}
$$

where $V_j = \mathbb C^{2j+1}$ is the unique irreducible rep of $\mathrm{SU}(2)$
of dimension $2j+1$. An orthonormal basis of $L^2(\mathrm{SU}(2))$ is the
*Wigner $D$-functions*

$$
   \sqrt{2j+1}\; D^{j}_{m\, m'}(g), \qquad
     j = 0, \tfrac12, 1, \tfrac32, \ldots,\ \ -j \le m, m' \le j. \tag{3}
$$

Two structural facts are used below:

(PW-1) **Convolution acts diagonally.** Any *central* function (constant on
       conjugacy classes) acts on each $V_j \otimes V_j^*$ block by a scalar
       — this is the SU(2) generalization of "Fourier multipliers act
       diagonally" in the abelian case. Non-central functions mix
       $(m, m')$ within a block but preserve $j$ (this is what's NEW vs the
       abelian case).

(PW-2) **Convolution by characters annihilates non-isotypic components.** For
       any irrep $\sigma$ the character $\chi_\sigma(g) = \mathrm{tr}\,\pi_\sigma(g)$
       satisfies
       $\chi_\sigma * D^{j}_{m m'} = \delta_{j\sigma}\,(2j+1)^{-1}\,D^{j}_{m m'}$.

### 2.2. Bijection to the GeoVac Fock-graph index $(n, l, m_l)$

The Fock projection of Paper 7 sends the hydrogenic $(n, l, m_l)$ shell
labels onto $S^3$ harmonics. The standard identification (Paper 7 §VI; Avery
1989; Wen–Avery JMP 26, 396 (1985)) is

$$
   Y^{(3)}_{n,\, l,\, m_l}(\chi, \theta, \phi)
   \;=\; R_{n,l}(\chi)\, Y_{l m_l}(\theta, \phi)
   \quad\Leftrightarrow\quad
   D^{(n-1)/2}_{m_L, m_R}(g),
$$

with the precise correspondence between the SO(4) principal QN $n$ and the
SU(2)$_L \times$ SU(2)$_R$ irreps $((n-1)/2, (n-1)/2)$ — the *symmetric*
SU(2) bi-rep of total dimension $n^2$ — being

$$
   n \;=\; 2j + 1, \qquad
   (l, m_l) \;\leftrightarrow\; \text{diagonal Clebsch–Gordan basis of}\ V_j \otimes V_j^*. \tag{4}
$$

In words: the Fock $n$-shell is exactly one Peter–Weyl block. The Fock
$(l, m_l)$ labels diagonalize the *Casimir of the SO(3) diagonal subgroup*
inside SO(4) = SU(2)$_L \times$ SU(2)$_R / \mathbb Z_2$. The dimension count
checks: $V_j \otimes V_j^* = \bigoplus_{l=0}^{2j} V_l$, so $\sum_{l=0}^{n-1}
(2l+1) = n^2$.

The bijection $(n, l, m_l) \leftrightarrow (j; \text{CG reduction})$ is
*explicit and computable*, and is in fact already what
`geovac/so4_three_y_integral.py` (R3.1) implements — its Avery–Wen–Avery
3-Y integral is the SO(4) Clebsch–Gordan integral written in the radial-Gegenbauer
form rather than in the Wigner-$3j \otimes 3j$ form.

### 2.3. Index truncation in PW language

The projection $P_{n_{\max}}$ becomes

$$
   P_{n_{\max}} \;=\; \bigoplus_{2j+1 \le n_{\max}} \mathbb 1_{V_j \otimes V_j^*}, \tag{5}
$$

i.e. retain all Peter–Weyl blocks with $j \le (n_{\max}-1)/2$. This is
*exactly* the SU(2) analog of Connes–vS's $P_n = \sum_{|k| \le n-1}
|e_k\rangle\langle e_k|$ on $L^2(S^1)$.

---

## §3. Propagation number prop(O_{n_max}) = 2 in PW language

The R2.1 result $\mathrm{prop}(\mathcal O_{n_{\max}}) = 2$ for $n_{\max} \in
\{2, 3, 4\}$ has a structural reading on the PW side. The multiplier matrix
of $D^{N}_{M_L M_R} \in V_J \otimes V_J^*$ (with $2J+1 = N$) acting on
$L^2(\mathrm{SU}(2))$ by left-translation has support on triples $(j, j', J)$
satisfying the SU(2) triangle $|j - j'| \le J \le j + j'$. The map

$$
   V_j \otimes V_j^* \;\longrightarrow\; V_{j'} \otimes V_{j'}^*
$$

induced by such a multiplier is non-zero precisely when
$|n - n'| + 1 \le N \le n + n' - 1$ — which is the SO(4) selection rule
already in `geovac/operator_system.py`, here re-derived from the SU(2)
Clebsch–Gordan triangle.

**Proposition 3.1 (PW reading of prop = 2).** *Products of two
Peter–Weyl multipliers $(D^{N_1}_{*}, D^{N_2}_{*})$ on $L^2(\mathrm{SU}(2))$
restricted to the truncation $j \le (n_{\max}-1)/2$ span $M_{N(n_{\max})}(\mathbb C)$ for
every $n_{\max} \ge 2$.*

*Sketch.* The single-multiplier image already supplies all rank-1 matrix
units $|j, m_L, m_R\rangle\langle j', m_L', m_R'|$ with $j' \in \{j, j \pm
\tfrac12, j \pm 1\}$ — the off-diagonal ones via $D^{1/2}_{**}$ multipliers
(SO(4) raising/lowering). The diagonal-by-shell rank-1 units that come from
"raise-then-lower" within the same shell live in $\mathcal O^2$ but not
$\mathcal O$, and these are exactly the boundary-deficit entries identified in
the R2.1 witness pair (Memo §2.5: the $(2,0,0) \to *(2,1,0)\to (2,0,0)$
"corner" diagonal entry has $14.9\%$ Frobenius residual against $\mathcal O$).
Two products supply them via the spectral identity
$D^{1/2}_* \cdot (D^{1/2}_*)^* = $ (diagonal in shell, off-diagonal in $l$).

This is the same mechanism Connes–vS use for the Toeplitz $S^1$ system
(Proposition 4.2, lines 1626–1653). The structural alignment of WH1
operator-system invariants with the Toeplitz benchmark is therefore not
just a numerical coincidence at $n_{\max} \le 4$ but a direct
consequence of the SU(2) tensor-product / Clebsch–Gordan algebra.
$\square$

---

## §4. The Wasserstein–Kantorovich state-space metric

### 4.1. Pure-state distance vs full-state-space distance

The R2.3 / R3.1 / R3.2 sprints computed the Connes distance Eq. (1) on
*pure node-evaluation states* $\varphi_v(x) = \langle v|x|v\rangle$ for $v$ a
single basis vector. **This is not the right comparison object for GH
convergence.** The natural GH-convergence target is the *full state space*
$\mathcal S(\mathcal O_{n_{\max}})$ — a compact convex set whose extremal
points are determined by the algebraic geometry of the operator system (per
Connes–vS §3.2, Theorem 3.6 for $S^1$ at $n=3$: the extreme states form
a Möbius strip).

The right limit object is

$$
   \mathcal P(S^3) \;=\; \{\text{Borel probability measures on}\ S^3\}
$$

with the *Monge–Kantorovich distance*

$$
   d_{\mathrm{Wass}}(\mu, \nu) \;:=\; \sup_{f\colon \|\nabla f\|_\infty \le 1}
   \Big| \int f\, d\mu - \int f\, d\nu \Big| , \tag{6}
$$

(equivalently the optimal-transport cost with cost = round-$S^3$ geodesic
distance, by Kantorovich–Rubinstein duality, all $\mu, \nu$ probability
measures).

### 4.2. The R3.1 / R3.2 anti-correlation finding, properly framed

R3.1 and R3.2 found that on *pure node-evaluation states* the Connes
distance is anti-correlated with the Fock-graph distance (Pearson nz
$-0.36$ at $n_{\max}=3$ scalar; $-0.26$ spinor R3.2). This is **not a
falsification of GH convergence**. Pure node-evaluation states $\varphi_v$
are not the natural limit objects — they are *singular* probability
measures concentrated on conjugacy classes of $g \in \mathrm{SU}(2)$ (in
the limit) or on individual harmonics $Y^{(3)}_{n l m}$ (at finite cutoff).
What GH convergence asks is whether *every* state $\varphi \in \mathcal
S(\mathcal O_{n_{\max}})$ can be matched to a probability measure $\mu \in
\mathcal P(S^3)$ with metric distortion vanishing as $n_{\max} \to \infty$,
*after* allowing convex combinations and weak-$*$ limits.

The cleanest reading of the R3.1 / R3.2 anti-correlation is therefore:
**at finite $n_{\max}$ the pure node-evaluation states $\{\varphi_v\}_v$
are not a GH-meaningful approximation of the Dirac masses
$\{\delta_g\}_{g \in S^3}$; the natural finite approximation is the
*Berezin lift* (see §5.1 below).**

This squares with Connes–vS's own $S^1$ result (their Theorem 4.20):
on $S^1$ the truncated-system Connes distance is *larger* than the
Kantorovich distance — they prove an inequality, not an identity —
and they hand the GH-convergence statement off to a future paper. The
GeoVac analog should be expected to behave the same way, and the
finite-$n_{\max}$ pure-state non-monotonicity is a feature, not a bug.

---

## §5. The proof shape on SU(2)

### 5.1. The Berezin / spectral-Fejér map

The key technical tool in Leimbach–van Suijlekom 2024 is the *spectral
Fejér kernel*: a positive trigonometric polynomial $K_{n_{\max}} \in
C^\infty(\mathbb T^d)$ with $\int K_{n_{\max}} = 1$, supported in Fourier
modes $|k| \le n_{\max}$, whose $L^1$ mass concentrates near $g = e$ as
$n_{\max} \to \infty$. Convolution by $K_{n_{\max}}$ defines a
unital completely positive (UCP) map

$$
   \rho_{n_{\max}} \colon C^\infty(\mathbb T^d) \;\to\; \mathcal O_{n_{\max}}
$$

(by Fourier-coefficient truncation) with a UCP partial inverse
$\sigma_{n_{\max}} \colon \mathcal O_{n_{\max}} \to C^\infty(\mathbb T^d)$
(by Schur multiplier). The composite $\sigma \circ \rho$ acts as Fourier
multiplication by $|\hat K_{n_{\max}}|^2$, which converges to identity
because $K_{n_{\max}}$ is a *good kernel* (mass concentrates).

The "good-kernel" property is quantified by

$$
   \gamma_{n_{\max}} \;:=\; \int_{\mathbb T^d} |K_{n_{\max}}(g)|\,
   d_{\mathrm{round}}(e, g)\, dg \;\xrightarrow{n_{\max}\to\infty}\; 0,
$$

which gives the operator-Lipschitz bound

$$
   \big\| f - K_{n_{\max}} * f \big\|_\infty
   \;\le\; \gamma_{n_{\max}} \cdot \big\|[D, f]\big\|_{\mathrm{op}} . \tag{7}
$$

Eq. (7) is the heart of the Leimbach–vS proof. It is what makes
$\rho \circ \sigma$ and $\sigma \circ \rho$ *near-isometries* in the
appropriate compact-quantum-metric-space sense.

### 5.2. Non-abelian analog: positive central trigonometric polynomials

The flat-case proof uses the abelian Fourier basis $\{e^{ik\cdot x}\}$.
The SU(2) analog has two levels:

(a) **The character lattice.** Replace $\{e^{ik\cdot x}\}$ by
    $\{\chi_j(g)\}_{j \in \tfrac12\mathbb Z_{\ge 0}}$, the irreducible
    characters. These are central (constant on conjugacy classes), so the
    convolution algebra they generate is abelian — even though SU(2) itself
    is not. This is the **first-level abelianization**.

(b) **Schur multipliers.** Replace the abelian Schur transference
    (Leimbach–vS §3) by the SU(2) Schur transference, which for a
    *central* multiplier $m$ is exactly the same statement as in the
    abelian case: $\|S_m\|_{\mathrm{cb}} = \|m\|_{\mathrm{Fou}}$. For
    *non-central* multipliers the statement requires the
    Bożejko–Picardello / Haagerup non-abelian Schur–Fourier transference
    inequality, which is known to hold on compact connected Lie groups
    with a constant that depends on the group (for SU(2) it is bounded by
    a small explicit constant — De Cannière, Haagerup, *Multipliers of
    the Fourier algebras of some simple Lie groups and their discrete
    subgroups*, Amer. J. Math. 1985).

The candidate central good kernel on SU(2) is

$$
   K_{n_{\max}}(g) \;:=\; \frac{1}{Z_{n_{\max}}} \,
     \bigg|\sum_{j \le (n_{\max}-1)/2} a_j\, \chi_j(g)\bigg|^2, \tag{8}
$$

with $a_j$ chosen so that $K_{n_{\max}}$ is positive (automatic from the
$|\cdot|^2$ structure), $\int K = 1$ (set by the normalizing constant
$Z_{n_{\max}}$), and the *mass-concentration* property

$$
   \gamma_{n_{\max}} \;=\; \int_{\mathrm{SU}(2)} K_{n_{\max}}(g)\,
   d_{\mathrm{round}}(e, g)\, dg \;\to\; 0
$$

holds. The natural choice for $a_j$ is the *spherical-harmonic Dirichlet
kernel coefficients squared*, which on SU(2) is the $j$-th Chebyshev-of-the-second-kind
polynomial weight: $a_j = \sqrt{2j+1}$. This is the SU(2) analog of the
square of the spherical Dirichlet kernel that appears in Leimbach–vS
Eq. (2.6) ($\hat K = \mathcal N_L / \mathcal N_B$, the ratio of lens-count
to ball-count).

### 5.3. The proof shape, lemma-by-lemma

I now state the lemma chain that the SU(2) analog should follow. Each lemma
mirrors a known Leimbach–vS or Hekkelman lemma; SU(2)-specific obstructions
are flagged.

**Lemma 5.1 (Peter–Weyl truncation as UCP map).** Let
$P_{n_{\max}}$ be the index truncation Eq. (5). Then the *compression*
$\rho_{n_{\max}}\colon C^\infty(\mathrm{SU}(2)) \to \mathcal O_{n_{\max}}$
defined by $\rho_{n_{\max}}(f) := P_{n_{\max}} M_f P_{n_{\max}}$ (with
$M_f$ multiplication by $f$) is a UCP map.

*Proof sketch.* Standard: compression by a self-adjoint orthogonal projection
of a $C^*$-algebra-image into $\mathcal B(P\mathcal H)$ is UCP (Stinespring
dilation, Choi–Effros). The spectral structure of $D$ is irrelevant. $\square$

**Lemma 5.2 (Spectral Fejér kernel exists on SU(2)).** Define $K_{n_{\max}}$
by Eq. (8) with $a_j := \sqrt{2j+1}$ for $2j+1 \le n_{\max}$, $a_j := 0$
otherwise. Then:

  (i) $K_{n_{\max}} \ge 0$ (trivial from $|\cdot|^2$).

  (ii) $\int_{\mathrm{SU}(2)} K_{n_{\max}}(g)\,dg = 1$ (set by $Z_{n_{\max}}$).

  (iii) $K_{n_{\max}}$ is central (depends only on the conjugacy class, i.e.
        on the rotation angle $\chi(g)$).

  (iv) $\gamma_{n_{\max}} \to 0$ as $n_{\max} \to \infty$ at the rate
       $\gamma_{n_{\max}} = O(1/n_{\max})$.

*Proof sketch.* (i)–(iii) are direct from Peter–Weyl. (iv) is a Cesàro /
Fejér computation: with $\chi_j(g) = \sin((2j+1)\chi/2)/\sin(\chi/2)$
and $\chi = $ rotation angle, the kernel reduces to a square of a
Dirichlet sum in $\chi$, integrated against the SU(2) Haar measure
$dg = (1/\pi) \sin^2(\chi/2)\, d\chi$ on $[0, 2\pi]$. The standard Fejér
estimate $\int |K_n(\chi)|\, |\chi|\, d\chi = O(1/n)$ then transports
verbatim. $\square$

**Lemma 5.3 (Lipschitz bound).** For $f \in C^\infty(\mathrm{SU}(2))$,

$$
   \big\| f - K_{n_{\max}} * f \big\|_\infty
   \;\le\; \gamma_{n_{\max}} \cdot \big\| [D_{\mathrm{CH}}, M_f] \big\|_{\mathrm{op}} ,
$$

where $D_{\mathrm{CH}}$ is the Camporesi–Higuchi Dirac on $L^2(S^3,\Sigma)$
and $M_f$ is multiplication by $f$.

*Proof sketch.* This is where the *first non-flat obstruction* appears.
On $\mathbb T^d$ the operator $\|[D, M_f]\|$ equals the gradient norm
$\|\nabla f\|_\infty$ exactly. On $S^3$ with the Camporesi–Higuchi Dirac the
identity becomes

$$
   \|[D_{\mathrm{CH}}, M_f]\|_{\mathrm{op}} \;=\; \|\nabla f\|_\infty
   \;+\; (\text{torsion-and-curvature correction}),
$$

with the correction bounded uniformly by the supremum of the spinor-bundle
connection coefficients (i.e. by the SU(2) structure constants times
$\|\nabla f\|_\infty$). Both terms are bounded by the round-S³ Lipschitz
constant of $f$ up to a *fixed multiplicative constant* depending only on
the unit S³ geometry; this constant is benign and the inequality
$\|f - K * f\|_\infty \le C \gamma_{n_{\max}} \|f\|_{\mathrm{Lip}}$ goes
through with a different — but still vanishing — rate. **Verification path:**
Lemma 5.3 is the *one* lemma where direct computation in
`geovac/spinor_operator_system.py` would be load-bearing — compare
$\|[D_{\mathrm{CH}}, M_f]\|_{\mathrm{op}}$ to $\|\nabla f\|_\infty$ for
the multiplier matrices already built. R3.2 has the infrastructure ready.

**Lemma 5.4 (Berezin reconstruction).** Define the Berezin lift
$\sigma_{n_{\max}}\colon \mathcal O_{n_{\max}} \to C^\infty(\mathrm{SU}(2))$
by

$$
   \sigma_{n_{\max}}(T)(g) \;:=\; (2j_{\max}+1)\, \mathrm{tr}\big( T \cdot
   K_{n_{\max}}(g \cdot )\big), \qquad j_{\max} = (n_{\max}-1)/2,
$$

where $K_{n_{\max}}(g\cdot)$ acts on $\mathcal H_{n_{\max}}$ by
left-translation by $g$. Then $\sigma_{n_{\max}}$ is UCP, and
$\sigma_{n_{\max}} \circ \rho_{n_{\max}} = M_{|\hat K_{n_{\max}}|^2}$
(Fourier multiplication by the central function $|\hat K|^2 \in
\mathcal Z(C(\mathrm{SU}(2)))$).

*Proof sketch.* Berezin–Toeplitz quantization on $S^3 = \mathrm{SU}(2)$
is a standard construction (Hawkins, *Quantization of equivariant vector
bundles*, CMP 1997; Bordemann–Meinrenken–Schlichenmaier, Berezin–Toeplitz
on Kähler manifolds). The CG decomposition $V_j \otimes V_j^* \to
\bigoplus_l V_l$ is the SU(2)-equivariant analog of the abelian
Fourier-mode decomposition. UCP follows from the positivity of $K$. The
identity $\sigma\rho = M_{|\hat K|^2}$ is the convolution-vs-multiplication
duality, written in the SU(2) PW basis. $\square$

**Theorem 5.5 (GH convergence on S³).** *The state-space metric of
$\mathcal O_{n_{\max}}$ with the operator-system Connes distance Eq. (1)
converges in Gromov–Hausdorff distance to the round S³ Wasserstein–Kantorovich
metric:*

$$
   d_{\mathrm{GH}}\Big( (\mathcal S(\mathcal O_{n_{\max}}), d_{D_{n_{\max}}}),\,
   (\mathcal P(S^3), d_{\mathrm{Wass}}) \Big) \;\xrightarrow{n_{\max}\to\infty}\; 0.
$$

*Proof shape.* Use the Hekkelman criterion (LMP 2022, Hekkelman’s master’s
thesis Theorem 4.3.1) plus Leimbach–vS Theorem 1.1:

  (1) The maps $\rho_{n_{\max}}, \sigma_{n_{\max}}$ are UCP and unital.

  (2) The composites $\sigma_{n_{\max}} \rho_{n_{\max}} \to \mathrm{id}$
      and $\rho_{n_{\max}} \sigma_{n_{\max}} \to \mathrm{id}$ in
      operator-Lipschitz norm by Lemmas 5.3 and 5.4.

  (3) The Lipschitz norms agree up to a vanishing error:
      $\|[D_{n_{\max}}, \rho_{n_{\max}}(f)]\| \le \|[D, M_f]\| + o(1)$,
      directly from Lemma 5.3.

  (4) Apply Latrémolière's quantum GH propinquity (Trans. AMS 368, 2016)
      or Rieffel's quantum GH distance (Mem. AMS 168, 2004) to convert
      the UCP-near-isometry into a GH bound.

The dual statement (b) — limit identification with $(\mathcal P(S^3),
d_{\mathrm{Wass}})$ — uses the Kantorovich–Rubinstein duality and the
fact that *any* compact-quantum-metric-space limit of compressions by good
kernels of $C(M)$ is forced to be the Wasserstein–Kantorovich space of
$M$ (Rieffel 1999/2004; D'Andrea–Lizzi–Martinetti 2014, J. Geom. Phys. 82).
$\square$

---

## §6. Stress-test: three obstructions visible in R3.1 / R3.2

### 6.1. Obstruction A — Pure-state anti-correlation does not transport

R3.1 / R3.2 found Pearson nz $-0.36$ / $-0.26$ between the Connes distance
and the Fock-graph hop distance on pure node-evaluation states. The proof
shape in §5 does *not* claim monotonicity on pure states. It claims
GH convergence on the *full state space* with the Wasserstein metric,
which involves convex combinations and is sensitive only to *macroscopic*
distance distributions. **The obstruction is therefore not a real
obstruction to Theorem 5.5; but it is a real obstruction to a
*stronger* statement** — namely that the truncation is a discretization
of the round-$S^3$ geodesic on individual points. Connes–vS Theorem 4.20
on the circle is the precedent: at finite truncation the distance is
*larger* than Kantorovich (inequality, not identity), with equality only
in the limit. The GeoVac analog should be expected to be the same.

### 6.2. Obstruction B — Truthful-CH Dirac shell-degeneracy obstruction

R3.2 found that on the Weyl spinor sector with the *truthful*
Camporesi–Higuchi Dirac (eigenvalue $n_{\mathrm{Fock}} + 1/2$ depending
only on $n$), 24 of 28 cross-shell pure-state pairs at $n_{\max}=2$ have
*infinite* Connes distance. The kernel of $[\tilde D, \cdot]$ on the
truncated operator system contains every multiplier $M_{N, L=0, M=0}$ that
is shell-diagonal (large family).

**This obstructs the simple-minded Leimbach–vS proof but does not obstruct
the proof shape in §5.** The SU(2) Dirac on $L^2(S^3, \Sigma)$ has
*both* chiralities, with eigenvalues $\pm(n_{\mathrm{Fock}} + 1/2)$. The
*full Dirac sector* (R3.5) has eigenvalues that *do* depend on the
chirality grading, breaking the shell degeneracy at the operator level
because the $\gamma_5$-grading anti-commutes with $D$ but commutes with
multiplication operators. The full Dirac sector is necessary even for the
$S^1$ case in Hekkelman's proof — Hekkelman uses $D = -i\,d/d\theta$ which
is *unbounded on every shell* and breaks degeneracy automatically.

**Implication for the Track-TS-A sprint:** R3.5 (extending the spinor lift
from Weyl to full Dirac) is a *prerequisite* for the proof shape in §5,
not an optional follow-up. It is also a tractable computational sprint
on the order of the R3.2 effort.

### 6.3. Obstruction C — Non-abelian Schur multiplier inequality

The Leimbach–vS proof uses Schur–Fourier transference at *full strength*
on $\mathbb T^d$ (i.e., equality of cb-norms). On SU(2) the analog is
the Bożejko–Picardello–Haagerup inequality, which gives
$\|S_m\|_{\mathrm{cb}} \le C_G \|F_m\|_{\mathrm{cb}}$ with $C_G > 1$ a
group-dependent constant. For *central* multipliers (constant on
conjugacy classes) the inequality is an equality and the flat proof
transports directly.

The candidate spectral Fejér kernel Eq. (8) is *central* by construction
(a function of $\chi_j(g)$, which are central). So the relevant Schur
transference instance is the easy one. The non-central multipliers that
appear in $\mathcal O_{n_{\max}}$ from the SO(4) Clebsch–Gordan structure
do not enter the Berezin reconstruction $\sigma$ — they enter only the
*compression* $\rho$, which is automatically UCP.

**Net assessment:** Obstruction C is a non-obstruction *for this proof*,
provided one is careful about which multipliers enter which step.

---

## §7. Where the path can break

Three places where the §5 proof can break, in decreasing order of
likelihood:

(i) **Lemma 5.3 — $\|[D_{\mathrm{CH}}, M_f]\| \approx \|\nabla f\|_\infty$
    fails by a non-uniform amount.** The torsion-and-curvature correction on
    a non-flat manifold can in principle blow up at high frequencies. On
    constant-curvature compact manifolds this does *not* happen (the
    Lichnerowicz formula on a parallelizable manifold gives a uniform
    correction), and SU(2) is parallelizable, so this should hold. But
    "should" is doing real work here, and the verification cost is a
    multi-line spinor-bundle-curvature computation.

(ii) **Lemma 5.2(iv) — the SU(2) good-kernel rate $\gamma_{n_{\max}}$ may be
     sub-optimal.** The analog of the $T^d$ rate is $O(1/n_{\max})$. SU(2)
     compactness might give $O(\log n_{\max} / n_{\max})$ or
     $O(1/n_{\max})$ depending on how the lens-vs-ball ratio behaves on a
     non-abelian compact group. Sub-optimal rate doesn't kill the theorem,
     but it weakens the GH-distance estimate.

(iii) **Berezin–Toeplitz quantization on $S^3 \to S^2$, not on $S^3$ itself.**
      The standard Berezin–Toeplitz is on Kähler manifolds, and $S^3$ is
      not Kähler. The construction in Lemma 5.4 uses the parallelizable
      structure of $S^3 = \mathrm{SU}(2)$ as a substitute (Hawkins 1997).
      Whether Hawkins's setup gives a "good" Berezin quantization in the
      sense needed for Theorem 5.5 is not obvious from the abstract; it
      probably does, but a careful audit is needed.

If (i) or (iii) actually breaks, the proof modifies:

- (i) breaks → use the *scalar* Laplacian Dirac proxy from Lemma 5.2 in
      place of the spinor CH Dirac. The flat proof works for the scalar
      Laplacian on $T^d$; the same should be true on SU(2). This sidesteps
      the spinor bundle entirely, at the cost of a less physically
      meaningful Dirac. **This is the (d) modified statement in the
      verdict line: prove GH convergence with the scalar Laplacian Dirac
      proxy, leaving the spinor-Dirac case to a follow-up.**

- (iii) breaks → use the *equivariant Berezin* construction of Hawkins
      (CMP 1997) directly, with characters of SU(2) replacing characters
      of $T^d$. This is what R3.5 (full Dirac sector) is set up for.

---

## §8. Verdict and lemma roadmap

**(b) Reachable, but harder than the flat case. Time estimate: 4–8 weeks
of focused work.**

The proof shape in §5 transports cleanly through Peter–Weyl on SU(2). The
key new ingredients are Eq. (8)'s central spectral Fejér kernel and the
non-flat torsion-correction in Lemma 5.3. Obstruction A (pure-state
anti-correlation) is a feature, not a bug. Obstruction B (CH shell
degeneracy) requires R3.5 as a prerequisite. Obstruction C is benign for
the central kernel.

**Next 5 lemmas, in order of dependency:**

(L1) **R3.5 — full Dirac sector.** Lift `geovac/spinor_operator_system.py`
     from Weyl to both chiralities. Verify that on the full sector the
     Connes distance is finite on every cross-shell pair under the
     truthful CH Dirac. *Estimated effort: 1–2 weeks (R3.2-class sprint).*
     Produces the sectoral substrate Lemmas 5.3 and 5.4 need.

(L2) **Lemma 5.2 (SU(2) good kernel).** Construct $K_{n_{\max}}$ as in
     Eq. (8); verify positivity, normalization, centrality, and the
     mass-concentration rate $\gamma_{n_{\max}} = O(1/n_{\max})$ at
     $n_{\max} = 2, 3, 4, 5$ in `sympy`/`mpmath`. *Estimated effort:
     1 week.* Produces the kernel that does the work.

(L3) **Lemma 5.3 (Lipschitz bound).** Compute
     $\|[D_{\mathrm{CH}}, M_f]\|$ vs $\|\nabla f\|_\infty$ for the
     existing R3.2 multiplier matrices on SU(2). Identify the
     torsion-correction prefactor explicitly. *Estimated effort:
     1–2 weeks (mathematician + computational verification).*
     Produces the Lipschitz comparison constant; if it diverges,
     fall back to (d).

(L4) **Lemma 5.4 (Berezin reconstruction).** Define $\sigma_{n_{\max}}$
     concretely on SU(2) using the equivariant Hawkins quantization.
     Verify $\sigma\rho = M_{|\hat K|^2}$ at $n_{\max} = 2, 3$ in
     exact arithmetic (the `geovac/so4_three_y_integral.py` machinery
     supports this). *Estimated effort: 1 week.* Produces the round-trip
     near-identity.

(L5) **Theorem 5.5 (GH convergence).** Assemble L1–L4 with
     Latrémolière propinquity / Rieffel quantum-GH-distance
     bounds. *Estimated effort: 1–2 weeks of mathematical writeup.*
     Produces the manuscript.

The expected output is a manuscript "Gromov–Hausdorff Convergence of
Spectral Truncations on $S^3$" suitable for J. Geom. Phys. or Adv. Math.,
naturally companion-cited with Leimbach–vS (Adv. Math. 439, 2024) and
Hekkelman (LMP 2022).

The R3.5 prerequisite suggests an interim sprint plan: **R3.5 first
(2 weeks), then a focused TS-A sprint (4–8 weeks)** for the full
manuscript. In the meantime, the R3.5 result is itself a clean WH1
output: it would either close the n-degeneracy obstruction
*by construction*, validating the Connes-vS framework on the
Camporesi–Higuchi Dirac, or surface a new structural obstruction that
itself becomes an interesting NCG observation.

---

## §9. Honest scope-statement

This memo proves nothing. It identifies what the proof shape *should* be,
and what the obstacles are. Three points to flag:

(a) The *form* of Eq. (8) is dictated by the Peter–Weyl structure, but the
    *coefficients* $a_j$ that make it a "good" kernel in the
    Leimbach–vS sense are guessed from the abelian analog. Lemma 5.2(iv)
    is the place where the guess is verified or refuted. If the natural
    coefficients give a worse rate, the GH theorem still goes through with
    an explicit weaker rate.

(b) The Berezin–Toeplitz step Lemma 5.4 invokes Hawkins's
    equivariant quantization on parallelizable manifolds. The literature
    on this is solid but not completely automated; an explicit verification
    on SU(2) at small $n_{\max}$ in the GeoVac codebase would be the
    safest way to ground the final theorem.

(c) The argument is *not* a no-go theorem in the modified-statement
    direction. If the spinor-CH version fails because of obstruction (i)
    in §7, the scalar-Laplacian version still goes through, and that
    weaker theorem is also publishable as a non-flat extension of
    Leimbach–vS. The verdict is therefore (b) for the strong version
    and (d)-as-fallback for the scalar version, with R3.5 sitting on the
    critical path either way.

**No-go fallback:** the only way the *whole* proof shape collapses is if
SU(2) does not admit a positive central good kernel — i.e. if
$\gamma_{n_{\max}}$ in Lemma 5.2(iv) does not vanish. This would be a
genuine surprise (it would in particular contradict the Lebedev–Levitan /
Bożejko results on SU(2) heat-kernel concentration), and it would *itself*
be a publishable observation. So even the worst-case outcome of TS-A is
non-trivial.

---

## Executive summary (200 words)

The Connes–van Suijlekom Gromov–Hausdorff convergence theorem on the
torus (Leimbach–van Suijlekom 2024, Adv. Math. 439, 109496) transports to
the GeoVac S³ case via Peter–Weyl analysis on SU(2). The key technical
ingredient — the spectral Fejér kernel — has a clean non-abelian analog
as a positive *central* trigonometric polynomial in SU(2) characters
(Eq. 8), making the Schur–Fourier transference of the flat proof go
through automatically (centrality is the abelianizing assumption). Two
obstructions visible in the round-3 / R3.1 / R3.2 numerics — the
pure-state anti-correlation and the truthful-CH-Dirac shell degeneracy —
are not obstructions to GH convergence on the *full state space* with
the Wasserstein–Kantorovich metric, but the second one requires extending
the spinor lift from the Weyl sector to the full Dirac sector (R3.5) as
a prerequisite. With that prerequisite in place, the proof reduces to
five lemmas (full Dirac sector R3.5; central good kernel; Lipschitz
bound; Berezin reconstruction; assembly) that map onto a focused 4–8
week sprint culminating in a J. Geom. Phys. / Adv. Math. manuscript.
The verdict is **(b) reachable but harder than the flat case**, with a
named modified statement **(d)** as the scalar-Laplacian fallback if
the spinor-bundle Lipschitz bound (Lemma 5.3) fails.

---

**End of memo.**
