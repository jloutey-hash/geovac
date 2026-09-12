# Literature scan: is there a NO-GO THEOREM behind the two-centre overcompleteness wall?

**Date:** 2026-09-11 | **Target:** the two-centre Coulomb–Sturmian metric wall
(`papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex`, `sec:molecular`;
CLAUDE.md §2 v5.10.15 "the conditioning cost is BASIS overcompleteness, not the metric")
| **Scope:** frames / Riesz bases, density theorems, Balian–Low, shift-invariant systems of
translates, locality-vs-conditioning trade-offs.

**Read-only scan.** Nothing under `papers/`, `geovac/`, `tests/`, `CLAUDE.md`, `CHANGELOG.md`
was modified.

**Overall verdict: GO — but the theorem that closes it is elementary, and it is not the one the
task expected.** The fancy machinery (Balian–Low, Beurling density, BCHL localization) does
*not* transfer. What does transfer is (i) a two-line Gram/completeness argument that proves
`lambda_min -> 0` outright with *one* verifiable hypothesis, (ii) the multiplication-operator
(fiberization) form of the Ron–Shen criterion, which Paper 60 has already derived without
naming it, and (iii) a one-line block-preservation lemma that proves the l-sparsity loss
*independently of conditioning* — which means GeoVac is currently conflating two distinct walls.

---

## 0. The object, restated in the right language

Let `H` be the molecular Hilbert space with the working inner product (the Shibuya–Wulfman
metric `S`, Paper 60 Eq. for `S_{mu'mu}`, equivalently the `V_0`-weighted product; the choice
matters — see V1 below).

* `F_A = {chi^A_i}` — Coulomb–Sturmians at shared scale `k`, centred at A. Complete in `H`.
* `F_B = {chi^B_j} = {U_R chi^A_j}` — the same set translated by `R`. Also complete.
* The molecular basis is `F = F_A ∪ F_B`; `G_N` = Gram of the first `N` of each.
* Measured: `lambda_min(G_N) ~ N^-2`, `lambda_max < 2`, `cond ~ N^1.85` (fitted window;
  derived asymptote `N^2`), all severity in `1 - sigma_max`.

Frame-theoretic vocabulary (Christensen; Heil):

* `F` is a **Bessel sequence** with bound `B` if `sum_k |<f,f_k>|^2 <= B ||f||^2`.
* `F` is a **frame** if additionally `A ||f||^2 <= sum_k |<f,f_k>|^2`.
* `F` is a **Riesz sequence** (= Riesz basis for its closed span) if
  `A sum |c_k|^2 <= || sum c_k f_k ||^2 <= B sum |c_k|^2` for all finitely supported `c`.

**The identification that unlocks the literature:** "the truncated Gram is bounded below
uniformly in `N`" is *exactly* "`F` is a Riesz sequence". `lambda_min(G_N) >= A` for all `N`
⟺ Riesz sequence with lower bound `A`. So GeoVac's measured wall is literally the failure of
the Riesz-sequence property, and the question "is there a no-go theorem?" is the question
"can a complete family plus a translate of itself be a Riesz sequence?"

**The distinction that must be kept straight, because the task's phrasing blurs it:** the
*lower frame bound* of the union does **not** collapse. It is `>= A_A + A_B`, i.e. better than
either half. What collapses is the *lower Riesz bound*. Frame operator `S = TT*` and Gram
`G = T*T` share their non-zero spectrum; the union has a perfectly good `S` and a `G` with a
kernel. Any search for "the lower frame bound collapsing for redundant systems" will come back
empty, because that is not what happens. Redundancy kills injectivity of the synthesis
operator `T`, not surjectivity.

---

## Q1. FRAME vs RIESZ BASIS — **GO, with proof**

### Q1.1 The qualitative answer, from the standard equivalence

**Theorem (Christensen, Thm 5.4.7 / Ch. 7 in the 2nd ed.; Heil, *A Basis Theory Primer*,
Ch. 7).** For a frame `{f_k}` in a Hilbert space, the following are equivalent:
(i) `{f_k}` is a Riesz basis; (ii) `{f_k}` is `omega`-independent (no non-trivial
`c in l^2` with `sum c_k f_k = 0`); (iii) `{f_k}` is minimal; (iv) `{f_k}` is **exact**
(removing any single element destroys the frame property); (v) the synthesis (pre-frame)
operator `T : l^2 -> H`, `Tc = sum c_k f_k`, is injective.

Apply it. The union of two Bessel sequences is Bessel (bounds add), and the union of a frame
with a Bessel sequence is a frame. So `F = F_A ∪ F_B` is a frame for `H`. It is **inexact**:
delete any `chi^B_j` and `F_A` alone still spans `H`. Hence by (iv) ⇒ (i), `F` is **never** a
Riesz basis, and by (v) `ker T != 0` — there exist genuine `l^2` linear dependencies among the
two-centre Sturmians.

**Answer to the literal question: no. The union of a complete Riesz basis with a complete
translate of itself is never again a Riesz basis, and this is not a near-miss — it fails by
the single cleanest characterisation in the subject.**

### Q1.2 The quantitative answer — and this one needs no frame hypothesis at all

The frame route above needs `F_A` to be a frame and `F_B` to be Bessel *in the same inner
product*. That is a real hypothesis in our weighted-metric setting (V1/V2 below). It can be
dispensed with entirely:

> **Proposition A (elementary, exact).** Let `H` be a Hilbert space, `{f_i}_{i>=1} ⊂ H`, and
> let `g != 0` lie in the closed linear span of `{f_i}`. Write `V_N = span{f_1,...,f_N}` and
> let `G_N` be the Gram matrix of `{f_1,...,f_N} ∪ {g_1,...,g_N}` where `g_1 = g`. Then
>
>     lambda_min(G_N)  <=  dist(g, V_N)^2 / (1 + ||c_N||^2)  <=  dist(g, V_N)^2  ->  0,
>
> where `c_N` are the least-squares coefficients of `g` on `V_N`.
>
> *Proof.* For any `c in C^N`, the unit vector `v = (c, -e_1)/sqrt(||c||^2+1) in C^{2N}` gives
> `v* G_N v = || sum_i c_i f_i - g ||^2 / (||c||^2 + 1)`. Since `lambda_min(G_N) <= v* G_N v`
> for every unit `v`, minimise over `c`. Completeness gives `dist(g,V_N) -> 0`. ∎

No Bessel hypothesis, no frame hypothesis, no localization, no translate structure. The only
input is: **one** member of the second family lies in the closed span of the first.

> **Corollary (the GeoVac wall, proved).** If the one-centre Coulomb–Sturmian set at shared
> scale `k` is complete in the molecular Hilbert space `H` in the working metric, then for the
> two-centre union `lambda_min(G_N) -> 0`, hence (with `lambda_max >= 1`) `cond(G_N) -> inf`.
> The wall is a theorem, not a measurement.

> **Corollary (rate, and it is falsifiable).** `lambda_min(G_N) <= eps_N^2`, where `eps_N` is
> the error of the best `N`-term *one-centre* approximation to a single displaced Sturmian.
> The measured `lambda_min ~ N^-2` therefore *predicts* `eps_N = O(N^-1)` for the one-centre
> expansion of a displaced Sturmian — a quantity GeoVac can compute independently. If the
> measured `eps_N` decays faster than `N^-1`, Proposition A's bound is not tight and something
> else is also going on; if it matches, the two-centre expansion rate *is* the conditioning law.

**Generalisation, free of charge.** `{f_i}` complete and `U` any bounded operator with
`U f_i != 0` ⇒ `{f_i} ∪ {U f_i}` is not a Riesz sequence. "Translate" is incidental; it is
completeness that does the work. This is the correct level of generality for the corpus's
claim, and it is *stronger* than what the task asked for.

### Q1.3 The operator-theoretic form — which is Paper 60's own algebra

If each half is orthonormal in the working metric (Paper 60 measures exactly this for the SW
intra-centre block: identity to `~1e-16`), then

    G = [[I, C],[C*, I]],   spec(G) = {1 ± sigma_k(C)},   cond(G) = (1+sigma_max)/(1-sigma_max),

with `sigma_k` the cosines of the principal angles between the two centre subspaces (Jordan
1875; Björck–Golub 1973; in QC, Amos–Hall 1961, King et al. 1967 — Paper 60 already cites
these). Now: `C = T_A* T_B` with `T_A, T_B` isometries, so
`||C|| = sup{ |<u,v>| : u in ran T_A, v in ran T_B, ||u||=||v||=1 }`. An isometry has closed
range, so `ran T_A = ran T_B = H` when both halves are complete, and the supremum is attained
at `u = v`: `||C|| = 1` **exactly**. Hence `lambda_min = 1 - sigma_max = 0`.

Note the eigenvector: the small eigenvalue `1 - sigma` belongs to `(1,-1)` — the **ungerade**
combination. This is an independent confirmation of Paper 60's measured observation that the
ill-conditioning is confined to the ungerade sector. The mechanism and the measurement agree.

The canonical-form reference for two subspaces at angle zero is **Halmos, "Two subspaces",
Trans. AMS 144 (1969) 381–389** (generic position, the `(P,Q)` canonical form that *is* the
principal-angle decomposition). The union-of-two-Riesz-sequences criterion in the literature is
**Christensen–Kim–Kim–Lim, "Angle criteria for frame sequences and frames containing a Riesz
basis", J. Math. Anal. Appl. (2008)**, and the underlying subspace fact is
`M + N closed <=> Friedrichs angle c(M,N) < 1` (Deutsch's survey). Our case fails the criterion
in its *first* clause, not the angle clause: a union of Riesz sequences is a Riesz sequence only
if `span(F_A) ∩ span(F_B) = {0}`, and here the two spans are **equal**. That is maximal failure,
not marginal failure — which is why no amount of increasing `R` fixes it (Paper 60's measured
`cond(S) -> 2.2` at `R = 10` bohr is a *prefactor* improvement at fixed `N`, not a change of
asymptotics; `sigma_max -> 1` at every `R > 0`).

---

## Q2. DENSITY THEOREMS — **BORDERLINE: right shape, hypothesis fails, and the conclusion is
weaker than Q1 anyway**

### Q2.1 Ramanathan–Steger

**J. Ramanathan and T. Steger, "Incompleteness of sparse coherent states", Appl. Comput.
Harmon. Anal. 2 (1995) 148–153** (DOI 10.1006/acha.1995.1010). For a Gabor/coherent system
`G(g,Lambda) = {e^{2 pi i b x} g(x-a)}_{(a,b) in Lambda}` with `Lambda ⊂ R^{2d}` an arbitrary
discrete set: if `G(g,Lambda)` is a frame then `D^-(Lambda) >= 1`; if it is a Riesz sequence
then `D^+(Lambda) <= 1`; hence a Riesz basis forces critical Beurling density
`D^- = D^+ = 1`. The proof runs through the Homogeneous Approximation Property plus a
comparison argument.

**Does it transfer? No.** The theorem is about phase-space translates of a *single* generator,
and "density" is Beurling density of a set in `R^{2d}`. Our index set is `(n,l,m)` at two
centres; there is no `R^{2d}`, no Beurling density, and no single generator. The *moral* —
redundancy is measured by a density and a Riesz basis forces the critical value — is exactly
right, and Proposition A delivers that moral in redundancy-2 form without needing a density at
all.

### Q2.2 Balan–Casazza–Heil–Landau

**R. Balan, P. G. Casazza, C. Heil, Z. Landau, "Density, overcompleteness, and localization of
frames. I. Theory", J. Fourier Anal. Appl. 12 (2006) 105–143** (DOI 10.1007/s00041-006-6022-0;
arXiv:math/0510360) and **"…II. Gabor systems", JFAA 12 (2006) 307–344**
(DOI 10.1007/s00041-005-5035-4; arXiv:math/0510361).

The relevant machinery, quoting the structure of the paper:

* **Def. 2.7(a)** `F = {f_i}_{i in I}` is `l^p`-*localized* with respect to a reference
  `E = {e_j}_{j in G}` (`G` a discrete abelian group) via a map `a : I -> G` if there is
  `r in l^p(G)` with `|<f_i, e_j>| <= r_{a(i)-j}` for all `i, j`.
* **Def. 2.1** lower/upper density `D^-(I,a)`, `D^+(I,a)` by counting `|I_N(j)|/|S_N(j)|`.
* **Def. 2.17** the *measure* `M(F;p,c)` — an average of `<f_i, f~_i>` (frame element against
  its canonical dual).
* **Example 2.18(a)** if `F` is a Riesz sequence then `M(F) = M^+(F) = M^-(F) = 1`.
* **Theorem 3.4 (Density–Relative Measure)** `M^-(F) = 1/D^+(I,a)` and `M^+(F) = 1/D^-(I,a)`.

**The derived statement the task was hoping for is real:** combining 2.18(a) with 3.4, an
`l^1`-localized frame with `D^-(I,a) > 1` has `M^+(F) = 1/D^- < 1 != 1`, hence is **not** a
Riesz sequence. That is precisely "a localized frame with redundancy > 1 cannot be a Riesz
basis", with redundancy identified with density.

**Does it transfer? Shape yes, hypothesis no.** The natural set-up would be: reference
`E = F_A` (must be a Riesz basis), `G = Z` or `Z^3` indexing `(n,l,m)`, `a` the identity on
each half so that `a(I)` covers `G` twice and `D^± = 2`. The obstruction is the
**`l^1`-localization hypothesis**: it demands `|<chi^B_j, chi~^A_i>| <= r_{i-j}` with
`r in l^1`. Paper 60's own derivation says the cross block is the finite section of a
multiplication operator with symbol `j_0(kR cot(chi/2))` on `chi in (0,pi)`. That symbol is
smooth at `chi = pi` but *oscillates without bound* as `chi -> 0+` (argument `-> inf`,
`j_0(x) = sin x / x`). A symbol with an oscillatory singularity at an endpoint has slowly
decaying Fourier coefficients; `l^1` off-diagonal decay of `C` is therefore **doubtful and
must not be assumed**. GeoVac can test it directly: does `sum_d max_n |C_{n,n+d}|` converge?

Second objection, structural rather than technical: BCHL's density is a genuinely
group-theoretic quantity, and `(n,l,m)` carries no translation structure — the map `a` would be
an artifice imposed to make the theorem quotable. **Verdict: do not lean on BCHL.** It would
buy a weaker conclusion (`not a Riesz basis`) than Proposition A already gives
(`lambda_min <= eps_N^2`, with a rate), at the cost of an unverified and probably false
hypothesis.

### Q2.3 Gröchenig localization — the right formalism, but for the *repair*, not the *no-go*

**K. Gröchenig, "Localization of frames, Banach frames, and the invertibility of the frame
operator", J. Fourier Anal. Appl. 10 (2004) 105–132** (DOI 10.1007/s00041-004-8007-1).
Localization theory (introduced independently by Gröchenig and by BCHL) is a machine for
proving *good* news: if a frame is intrinsically localized (polynomial or exponential
off-diagonal decay of `<f_i,f_j>`) **and the frame operator is invertible**, then the canonical
dual is localized too, and frame expansions converge in all the associated Banach spaces. The
engine is Wiener's-lemma inverse-closedness (Jaffard; Baskakov; Sjöstrand; Gröchenig–Leinert;
Sun).

This is the correct formalism for the question "can the dual/orthogonalized basis stay local?"
— and its answer for GeoVac is *no by hypothesis failure*: the frame operator is precisely what
is **not** boundedly invertible here. That is a sharper and more honest reading of the failure
of the biorthogonal repair (the Artacho–Milans del Bosch dual-basis route) than "the dual is
dense": the dual is dense *because the Wiener-lemma hypothesis is violated*, not as a brute
fact. See Q5 for the quantitative version.

---

## Q3. BALIAN–LOW — **STOP. It does not transfer, and the reason is instructive.**

### Q3.1 Exact statement

**Balian–Low Theorem (BLT).** Let `g in L^2(R)`, `a,b > 0` with `ab = 1` (critical density),
and let

    G(g,a,b) = { e^{2 pi i b m t} g(t - n a) }_{m,n in Z}.

If `G(g,a,b)` is a Riesz basis (in particular, an orthonormal basis) for `L^2(R)`, then

    ( integral t^2 |g(t)|^2 dt ) * ( integral xi^2 |g^(xi)|^2 dxi )  =  infinity.

Equivalently `t g notin L^2` or `xi g^ notin L^2`: the generator cannot be simultaneously
well localized in time and in frequency.

**History.** Balian (1981) and Low (1985) independently, both for *orthonormal* bases and both
with a gap. Battle (1988) gave the uncertainty-principle proof. The extension to *Riesz* bases
is attributed in the reference literature to **Daubechies, Coifman and Semmes** (and the
Zak-transform machinery to **Daubechies–Janssen 1993**); this attribution is the one point in
this scan where secondary sources disagree with the version in the task prompt — "Coifman /
Daubechies / **Jaffard / Journé**" conflates the BLT extension with the *Wilson basis* paper
(Daubechies–Jaffard–Journé 1991). The survey of record is **Benedetto–Heil–Walnut, JFAA 1
(1995) 355–402**.

### Q3.2 Which hypotheses are essential

| Hypothesis | Essential? | Why |
|:---|:---|:---|
| `ab = 1` (critical density) | **ESSENTIAL — and fatal for us** | At `ab < 1` the Gaussian generates a perfectly localized Gabor *frame*. Overcompleteness **removes** the obstruction entirely. BLT lives only at critical density. |
| Lattice / group structure | **ESSENTIAL to the mechanism** | The Zak transform `Zg` exists because the index set is `aZ x bZ`. Riesz basis ⟺ `Zg` bounded above and below; but a continuous quasi-periodic `Zg` has non-zero winding on the fundamental torus and **must vanish**. The obstruction is topological (winding/Chern), not analytic. |
| Time–frequency plane / Heisenberg | Essential to Battle's proof | uses `[t, d/d xi]`. |
| Riesz basis (vs frame) | Essential | BLT constrains systems that *are* Riesz bases. |

### Q3.3 How far it generalises — and the wall

Verified extensions, all still inside the group-action world:

* Higher dimensions, symplectic lattices (Gröchenig, *Foundations of Time-Frequency Analysis*,
  Ch. 8 "The Zak Transform"; exact theorem number not retrieved — **UNVERIFIED numbering**).
* Amalgam / weak BLT: Benedetto–Heil–Walnut (1995).
* **Irregular (non-lattice) Gabor systems at critical density:** Ascensi, Feichtinger,
  Kaiblinger, "Dilation of the Weyl symbol and Balian–Low theorem", Trans. AMS 366 (2014)
  3865–3880.
* Quantitative: Nitzan–Olsen, "A quantitative Balian–Low theorem", JFAA (2013),
  arXiv:1205.0163; Gautam, "A critical-exponent Balian–Low theorem", Math. Res. Lett. 15
  (2008) 471–483.
* Homogeneous groups: *Balian–Low type theorems on homogeneous groups*, Anal. Math. (2020).
* Finite dimensions: arXiv:1707.06449.
* **Finitely generated shift-invariant spaces** (the closest thing to our setting):
  **D. P. Hardin, M. C. Northington V., A. M. Powell, "A sharp Balian–Low uncertainty principle
  for shift-invariant spaces", Appl. Comput. Harmon. Anal. (2016)**,
  DOI 10.1016/j.acha.2016.05.001, arXiv:1510.04855. Generators translated along a **lattice**
  form a frame or Riesz basis for `V`; in the Riesz-basis case `V` must have extra invariance
  under a finer lattice. Conclusion: some generator satisfies
  `integral |x| |f_k(x)|^2 dx = infinity` (Fourier transform not in `H^{1/2}`), and `H^{1/2}`
  is sharp.

**Does any version apply to {complete set} ∪ {its translate}? No — and the failure is not
marginal.** Every version requires (a) the system to be generated by a group action (lattice of
translations, or translations × modulations), and (b) the system to be a Riesz basis at
critical density. Our system has neither: a two-element translate set `{0,R}` is not a group,
and our system is at redundancy 2 and is provably *not* a Riesz basis (Q1). BLT constrains the
localization of families that **are** Riesz bases; it is silent about families that are not.
It is the wrong side of the boundary.

### Q3.4 The instructive part — Wilson bases, and a live lead

**Daubechies, Jaffard, Journé, "A simple Wilson orthonormal basis with exponential decay",
SIAM J. Math. Anal. 22 (1991) 554–572** (DOI 10.1137/0522035). Wilson bases are exponentially
localized *orthonormal* bases built out of Gabor systems at twice critical density, by taking
real `cos`/`sin` (i.e. `±`) combinations of pairs of atoms at `±` frequency. They evade the
BLT completely. **The BLT obstruction is fragile to re-indexing.**

Two consequences for GeoVac:

1. *Negative:* stop expecting a Balian–Low-style no-go here. Even in the setting where such a
   theorem exists, a modest change of index structure defeats it.
2. *Positive, and worth noting:* the Wilson move — pass to `±` combinations of the two copies —
   is **exactly** the gerade/ungerade recombination Paper 60 already performs, and Paper 60
   already measures that it confines the ill-conditioning to one sector. The Q1.3 algebra
   explains why (the near-null eigenvector is `(1,-1)`). This does *not* fix `lambda_min -> 0`
   — the singular sector is still there — but it identifies the parity sector as the exact
   carrier of the obstruction. Any repair must act on the ungerade sector and only there.

---

## Q4. SHIFT-INVARIANT SPACES / SYSTEMS OF TRANSLATES — **GO in substance, with the base
space changed; the lattice form does NOT apply**

### Q4.1 The standard theorem, and whose it is

For `Phi = {phi_1,...,phi_r} ⊂ L^2(R^d)` and `E(Phi) = {phi_i(· - k) : k in Z^d}`, define the
**Gramian** fiber

    G_Phi(w)  =  [ sum_{k in Z^d} phi^_i(w + 2 pi k) conj(phi^_j(w + 2 pi k)) ]_{i,j=1..r},
    w in T^d.

Then:

* `E(Phi)` is Bessel with bound `B` ⟺ `||G_Phi(w)|| <= B` a.e.;
* `E(Phi)` is a **Riesz sequence** (a "stable basis" for `S(Phi)`) with bounds `A,B` ⟺
  `A I <= G_Phi(w) <= B I` a.e. — i.e. **the matrix symbol is bounded away from 0**;
* `E(Phi)` is a frame for `S(Phi)` ⟺ the same bounds hold on `ran G_Phi(w)` a.e.

**Attribution.** The fiberization / Gramian and dual-Gramian technique is **A. Ron and Z. Shen,
"Frames and stable bases for shift-invariant subspaces of `L_2(R^d)`", Canad. J. Math. 47
(1995) 1051–1094**; the structural precursor is **C. de Boor, R. DeVore, A. Ron, "The structure
of finitely generated shift-invariant spaces in `L_2(R^d)`", J. Funct. Anal. 119 (1994) 37–78**;
the range-function reformulation is **M. Bownik, "The structure of shift-invariant subspaces of
`L^2(R^n)`", J. Funct. Anal. 177 (2000) 282–309**; the standard survey is **A. Aldroubi and
K. Gröchenig, "Nonuniform sampling and reconstruction in shift-invariant spaces", SIAM Review
43 (2001) 585–620**. (Theorem *numbers* not retrieved — PDFs were not machine-readable;
the statements above are the standard ones and are safe, the numbering is **UNVERIFIED**.)

### Q4.2 Does the two-centre problem fit?

**Literally, no.** There is no `Z^d` action, no torus of translations, no Fourier fiberization.
`{0, R}` is a two-element set, not a lattice; every attempt to force the lattice picture fails.

**In substance, yes — and GeoVac has already built the right fibration without naming it.**
The *abstract content* of Ron–Shen is: *the Gram is a matrix-valued multiplication operator
over some base, and the Riesz property is uniform positivity of the fibre matrices.* The base
does not have to be a torus of translations. Paper 60 `sec:molecular` derives, for the
two-centre `s` sector in the sine (Fock-angle) basis, that the cross block `C` is the finite
section of a **multiplication operator** with symbol

    sigma(chi)  =  j_0( k R cot(chi/2) ),    chi in (0, pi),

so the union's Gram is unitarily equivalent to the `2 x 2` matrix-valued multiplication operator

    [[ 1, sigma(chi) ], [ sigma(chi), 1 ]],   fibre spectrum { 1 ± |sigma(chi)| }.

By the spectral theorem for multiplication operators, the Riesz criterion reads

    Riesz sequence  ⟺  ess-inf_chi ( 1 - |sigma(chi)| ) > 0  ⟺  ||sigma||_inf < 1.

Here `||sigma||_inf = 1`, **attained** at `chi = pi`. This is exactly the Ron–Shen criterion
transplanted onto the Fock base, and it fails for the sharpest possible reason: the symbol
touches its critical value rather than approaching it.

### Q4.3 The mechanism, coordinate-free — state it this way in the paper

`cot(chi/2)` is the Fock momentum variable: `p = k cot(chi/2)`, and `chi = pi` is `p = 0`. The
symbol is the angular average of the translation phase, `j_0(pR) = <e^{i p·R}>_angles`, and

> **translation by `R` acts on the Fock sphere as multiplication by a phase that becomes
> trivial at zero momentum.**

At the bottom of the momentum spectrum the two centres are *indistinguishable*, so the two
halves of the basis coincide there and the Gram symbol reaches 1. The maximum is quadratic
(`j_0(x) = 1 - x^2/6 + ...`, `cot(chi/2) ≈ (pi-chi)/2` near `chi = pi`), and band-limited
concentration at a quadratic symbol maximum gives the finite-section law
`1 - sigma_max = (kR)^2 pi^2 / (24 n^2) + o(n^-2)` — Paper 60's Eq. `eq:sigma_law`, i.e. the
measured `N^-2`. **The exponent 2 is the order of the symbol's maximum**, nothing else.

*Cross-reference:* the finite-section asymptotics at a quadratic symbol extremum are classical
Toeplitz theory — the sibling scan `debug/lit_scan/toeplitz_finite_section_memo.md` (same day)
identifies `eq:sigma_law` as Kac–Murdock–Szegő, with the constant factoring as
`pi^2 * (kR)^2/24`. This memo supplies the *qualitative* statement (`||sigma||_inf = 1`, hence
not a Riesz sequence — true for every `l` sector and every `R`); that memo supplies the *rate*
and its provenance. They agree and do not overlap.

This is the most transferable statement in the whole scan, because it names the escape route as
well as the wall: *the obstruction is localised at the `p -> 0` end.* A basis whose Fock-momentum
content is bounded away from `p = 0` has `||sigma||_inf < 1` and **is** a Riesz sequence. (This
is a structural observation, not a recommendation — excising the diffuse/low-momentum end is
what canonical orthogonalization with a threshold already does empirically in quantum chemistry,
and it costs completeness. Whether it can be done `l`-block-coherently is open and is *not*
settled by anything in this scan.)

**Scope caveat:** Paper 60's symbol derivation is for the two-centre **`s`-orbital** basis.
The `l > 0` sectors have not been reduced to a symbol, and the generalisation is not automatic
(the cross block becomes matrix-valued in `l` at each `chi`). Proposition A of Q1 covers all
sectors unconditionally; the symbol picture is currently `s`-only.

---

## Q5. OVERCOMPLETENESS AS AN OBSTRUCTION TO SPARSITY — **the sharpest result of the scan, and
it splits the wall in two**

The task calls this "the sharpest form of the question and the most valuable hit". It is, and
the answer is that GeoVac is currently running two *independent* walls together.

### Q5.1 Wall 1 — the `l`-sparsity loss has **nothing to do with conditioning**

> **Proposition D (block-preservation lemma; elementary, exact, no hypotheses).** Let
> `H_N = ⊕_l H_l` be a block decomposition. Suppose `X` is invertible and block-diagonal with
> respect to it, and `X* S X = D` with `D` block-diagonal. Then
> `S = X^{-*} D X^{-1}` is block-diagonal.
>
> **Contrapositive:** if `S` is *not* block-diagonal, then **no** block-diagonal congruence
> orthogonalizes it — not Löwdin `S^{-1/2}`, not canonical orthogonalization, not Cholesky,
> not any other scheme, and not for any value of the condition number.

The two-centre metric couples `l` (the expansion of a displaced function carries all `l`) while
preserving `m` (cylindrical symmetry about the bond axis). Therefore:

* **`m`-block selection survives orthogonalization** (block-diagonalise within each `m`);
* **`l`-block selection cannot survive it, ever** — even for a perfectly conditioned `S`.

This is exactly true and needs no analysis. It means the standing corpus sentence "Löwdin
`S^{-1/2}` is dense and destroys `l` block-selection" is *understated as stated*: it is not a
side-effect of ill-conditioning, it is a structural consequence of the metric being `l`-coupled.
The consequence for the framework is a genuine dichotomy with no middle: **keep `l`-sparsity and
a non-orthogonal metric, or take an orthogonal basis and lose Gaunt selection.** No third option
exists. (And by the same lemma the biorthogonal route fails identically: `S^{-1}` is
`l`-coupled for the same reason `S` is — which is the structural content of the Artacho–Milans
del Bosch dual-basis observation, **Phys. Rev. A 43, 5770 (1991)**,
DOI 10.1103/PhysRevA.43.5770.)

### Q5.2 Wall 2 — the *distance*-locality loss IS conditioning, quantitatively

For off-diagonal decay (as opposed to block structure), there is a genuine quantitative
trade-off theorem, and it is the best literature hit in the scan:

**S. Demko, W. F. Moss, P. W. Smith, "Decay rates for inverses of band matrices", Math. Comp.
43 (1984) 491–499.** For `A` symmetric positive definite with bandwidth `m` and spectrum in
`[lambda_min, lambda_max]`, `kappa = lambda_max/lambda_min`:

    | (A^{-1})_{ij} |  <=  C q^{|i-j|},    q = ( (sqrt(kappa) - 1)/(sqrt(kappa) + 1) )^{2/m}

(constants in the paper's normalisation). **As `kappa -> inf`, `q -> 1`**: the guaranteed decay
length grows like `m sqrt(kappa)` and the bound becomes vacuous. The bounds are qualitatively
sharp — the paper exhibits matching examples.

Extensions to general matrix functions, including `A^{-1/2}` (Chebyshev/Bernstein-ellipse
argument, so the rate is governed by the ellipse for `[lambda_min, lambda_max]` and degenerates
as `lambda_min -> 0`):

* **M. Benzi and G. H. Golub, "Bounds for the entries of matrix functions with applications to
  preconditioning", BIT 39 (1999) 417–438**, DOI 10.1023/A:1022362401426.
* **M. Benzi and N. Razouk, "Decay bounds and O(n) algorithms for approximating functions of
  sparse matrices", ETNA 28 (2007) 16–39.**
* **M. Benzi, P. Boito, N. Razouk, "Decay properties of spectral projectors with applications to
  electronic structure", SIAM Review 55 (2013) 3–64**, arXiv:1203.3953 — the definitive survey,
  and it treats the **non-orthogonal (overlap-matrix) representation explicitly**. Abstract,
  verbatim: *"our theory leads to a rigorous proof of the exponential off-diagonal decay
  ('nearsightedness') for the density matrix of gapped systems at zero electronic temperature in
  both orthogonal and non-orthogonal representations."* The gap is the hypothesis; GeoVac's
  metric is gapless in the limit.

For non-banded off-diagonal decay the same content appears as Wiener's-lemma inverse-closedness:

* **S. Jaffard, "Propriétés des matrices « bien localisées » près de leur diagonale et quelques
  applications", Ann. Inst. H. Poincaré Anal. Non Linéaire 7 (1990) 461–476** — invertibility
  on `l^2` implies the inverse inherits the same polynomial/exponential decay.
* **K. Gröchenig and M. Leinert, "Wiener's lemma for twisted convolution and Gabor frames",
  J. Amer. Math. Soc. 17 (2004) 1–18**; **A. Baskakov (1990)**; **J. Sjöstrand (1994/95)**;
  **Q. Sun, "Wiener's lemma for infinite matrices", Trans. AMS 359 (2007) 3099–3123**.
* **Norm-controlled** versions, which make the degradation explicit: **K. Gröchenig and
  A. Klotz, "Norm-controlled inversion in smooth Banach algebras, I", J. London Math. Soc. 88
  (2013) 49–64; II, Math. Nachr. 287 (2014) 917–937**; Shin–Sun. These bound the decay constants
  of `A^{-1}` by an explicit function of `||A^{-1}||` — so as `lambda_min -> 0` the guaranteed
  localization of `S^{-1}` and `S^{-1/2}` degrades in a controlled, computable way.

**The honest logical status.** "Well-conditioned + localized ⇒ the orthogonalizer is localized"
is a **theorem**. Its contrapositive at the level of bounds — "as conditioning fails, no
locality is guaranteed" — is what GeoVac needs, and it is a *degradation of guarantees*, not a
proof of impossibility. Demko–Moss–Smith's sharpness examples make it close to one, but the
literal statement "a well-conditioned orthogonalization of a redundant localized family MUST be
non-local" was **not found as a theorem** and I do not believe it is one in that generality.
What *is* a theorem in exactly the needed generality is Proposition D (Q5.1), which is stronger
where it applies (block structure) because it is unconditional.

### Q5.3 The nearest true "localized AND orthonormal is impossible" theorem

For completeness: a genuine topological no-go of this shape does exist, in condensed matter.

* **D. J. Thouless, "Wannier functions for magnetic sub-bands", J. Phys. C 17 (1984) L325**,
  DOI 10.1088/0022-3719/17/12/003 — a complete set of well-localized Wannier functions exists
  iff the sub-band carries no Hall current (Chern number 0).
* **C. Brouder, G. Panati, M. Calandra, C. Mourougane, N. Marzari, "Exponential localization of
  Wannier functions in insulators", Phys. Rev. Lett. 98, 046402 (2007)**,
  DOI 10.1103/PhysRevLett.98.046402 — exponentially localized Wannier basis ⟺ vanishing Chern
  class of the Bloch bundle.
* **G. Panati, "Triviality of Bloch and Bloch–Dirac bundles", Ann. Henri Poincaré 8 (2007)
  995–1011**, DOI 10.1007/s00023-007-0326-8.
* **D. Monaco, G. Panati, A. Pisante, S. Teufel, "Optimal decay of Wannier functions in Chern
  and quantum Hall insulators", Comm. Math. Phys. 359 (2018) 61–100**,
  DOI 10.1007/s00220-017-3067-7 — localization dichotomy: either exponentially localized
  composite Wannier functions exist, or **every** composite Wannier function has divergent
  `<x^2>`.

These are the real thing: orthonormality + localization is *topologically* obstructed. The
mechanism is the same winding/Chern argument that drives the BLT (Q3.2). **It does not apply to
GeoVac** — there is no Bloch bundle, no Brillouin torus, no Chern class in the two-centre
Sturmian problem. It is listed because it is the correct answer to "what would a genuine no-go
of this shape look like", and because it says where one would have to look: at a bundle over the
Fock base with a non-trivial invariant. Nothing in the corpus currently suggests one exists.
Do not chase this without a new structural reason.

---

## 6. Verdict and what GeoVac must still verify

### 6.1 Verdict: **GO**

The two-centre overcompleteness wall upgrades from measured to **proved**, via
**Proposition A** (Q1.2) — completeness of the one-centre set ⇒ `lambda_min(G_N) -> 0` ⇒
`cond -> inf`. The frame-theoretic restatement (Q1.1, Christensen/Heil) supplies the standard
vocabulary and the fact that genuine `l^2` null relations exist; the multiplication-operator
form (Q4.2/Q4.3) supplies the mechanism and the exponent; the block-preservation lemma
(Q5.1) separately and unconditionally proves the `l`-sparsity dichotomy.

### 6.2 The hypotheses that remain to be verified

**V1 (load-bearing, and the only one Proposition A needs).** *Is the one-centre Coulomb–Sturmian
set at shared scale `k` complete in the molecular Hilbert space, in the working metric?*
This is **not** a formality. The Sturmians are orthonormal/complete in the `1/r_A`-weighted
product; the working metric is Shibuya–Wulfman (`V_0 = Z_A/r_A + Z_B/r_B`, or Paper 60's
Hermitian `(2k^2)^{-1}<grad,grad> + (1/2)<·,·>` form). The two weighted spaces
`L^2(r_A^{-1} d^3r)` and `L^2(r_B^{-1} d^3r)` are **not** the same space — neither contains the
other (a function may be square-integrable against one weight and not the other near the
opposite nucleus). What must be checked is completeness of `F_A` in the *SW* completion.
Note the half-Gram is automatically bounded **below** in that metric
(`<chi^A_i, chi^A_j>_{V_0} = Z_A delta_{ij} + Z_B <chi^A_i, r_B^{-1} chi^A_j>`, second term PSD),
so the lower bound is free; it is the upper (Bessel) bound and the completeness that need work.

**V2 (needed only for the frame route, Q1.1).** `F_A` a frame *and* `F_B` Bessel in the same
metric. The Bessel bound is the fragile one: `r_B^{-1}` is not form-bounded by `r_A^{-1}`, so
`lambda_max` of the half-Gram may be unbounded as `N -> inf`. Paper 60 measures `lambda_max < 2`
on its window — consistent with the idealised block form, but that is a window reading, not a
proof. **If V2 fails, drop the frame language and quote Proposition A only**; Proposition A does
not need it.

**V3 (turns the bound into the law).** Measure `eps_N` = best `N`-term one-centre approximation
error for a displaced Sturmian. Proposition A predicts `lambda_min <= eps_N^2`; the measured
`N^-2` predicts `eps_N = O(N^-1)`. Agreement would tie the conditioning law to a one-centre
convergence rate; disagreement would mean the bound is slack and the `N^-2` has a second source.

**V4 (only if BCHL is to be cited at all — expected to FAIL).** `l^1` off-diagonal decay of the
cross block: does `sum_d max_n |C_{n,n+d}|` converge? The symbol's oscillatory endpoint
singularity at `chi -> 0` argues no. **Recommendation: do not cite BCHL.**

**V5 (scope).** The symbol/`sigma(chi)` picture is derived for the `s` sector only. Either
extend it to general `l` (matrix-valued symbol) or state the scope explicitly.

### 6.3 What must NOT be claimed

* Do **not** say "the lower frame bound collapses". It does not; the lower *Riesz* bound does.
  (§0.)
* Do **not** invoke Balian–Low. It applies at critical density to group-generated Riesz bases;
  we are at redundancy 2 and are not a Riesz basis. (Q3.)
* Do **not** invoke Beurling density / Ramanathan–Steger. No phase-space, no single generator,
  no density. (Q2.1.)
* Do **not** claim a theorem that a well-conditioned orthogonalization of a redundant localized
  family must be non-local. That was **not found** and is probably false in that generality;
  what is true is the `kappa`-dependent degradation of guarantees (Q5.2) and the unconditional
  block-preservation lemma (Q5.1).
* Do **not** claim the `l`-sparsity loss is caused by ill-conditioning. It is not — it is caused
  by the metric being `l`-coupled, and it would happen at `kappa = 1`. (Q5.1.) These are two
  independent walls and the corpus currently runs them together.

---

## 7. Citation list

Verification key: **[V]** opened/corroborated at publisher, arXiv or AMS; **[P]** bibliographic
data corroborated by two or more independent secondary sources but primary text not opened;
**[U]** unverified (numbering or a specific claim I could not confirm).

### Frames and Riesz bases
1. O. Christensen, *An Introduction to Frames and Riesz Bases*, 2nd ed., Applied and Numerical
   Harmonic Analysis, Birkhäuser/Springer, 2016. DOI 10.1007/978-3-319-25613-9. ISBN
   978-3-319-25611-0. Ch. 7 "Frames Versus Riesz Bases"; the frame ⟺ Riesz-basis equivalence is
   Thm 5.4.7 in the 1st edition. **[V]** book; **[U]** exact 2nd-ed. theorem number.
2. C. Heil, *A Basis Theory Primer*, expanded ed., Birkhäuser, 2011. ISBN 978-0-8176-4686-8.
   Ch. 7 (frames; exact frames are Riesz bases). **[P]**
3. O. Christensen, H. O. Kim, R. Y. Kim, J. K. Lim, "Angle criteria for frame sequences and
   frames containing a Riesz basis", *J. Math. Anal. Appl.* (2008),
   article S0022-247X(08)00607-0. **[P]** (publisher page returned HTTP 403; title, content and
   venue corroborated.)
4. P. R. Halmos, "Two subspaces", *Trans. Amer. Math. Soc.* **144** (1969) 381–389.
   DOI 10.1090/S0002-9947-1969-0251519-5. **[V]**
5. F. Deutsch, "The angle between subspaces of a Hilbert space", in *Approximation Theory,
   Wavelets and Applications* (S. P. Singh, ed.), Kluwer, 1995, 107–130. — `M+N` closed ⟺
   Friedrichs angle `< 1`. **[P]**
6. Å. Björck and G. H. Golub, "Numerical methods for computing angles between linear subspaces",
   *Math. Comp.* **27** (1973) 579–594. **[P]**

### Density / localization of frames
7. J. Ramanathan and T. Steger, "Incompleteness of sparse coherent states", *Appl. Comput.
   Harmon. Anal.* **2** (1995) 148–153. DOI 10.1006/acha.1995.1010. **[V]**
8. R. Balan, P. G. Casazza, C. Heil, Z. Landau, "Density, overcompleteness, and localization of
   frames. I. Theory", *J. Fourier Anal. Appl.* **12** (2006) 105–143.
   DOI 10.1007/s00041-006-6022-0; arXiv:math/0510360. **[V]** (abstract + Defs 2.1/2.7/2.12/2.17,
   Example 2.18(a), Thm 3.4 read from the arXiv HTML.)
9. R. Balan, P. G. Casazza, C. Heil, Z. Landau, "Density, overcompleteness, and localization of
   frames. II. Gabor systems", *J. Fourier Anal. Appl.* **12** (2006) 307–344.
   DOI 10.1007/s00041-005-5035-4; arXiv:math/0510361. **[V]**
10. K. Gröchenig, "Localization of frames, Banach frames, and the invertibility of the frame
    operator", *J. Fourier Anal. Appl.* **10** (2004) 105–132. DOI 10.1007/s00041-004-8007-1.
    **[V]**
11. C. Heil, "History and evolution of the density theorem for Gabor frames", *J. Fourier Anal.
    Appl.* **13** (2007) 113–166. **[P]**

### Balian–Low
12. R. Balian, "Un principe d'incertitude fort en théorie du signal ou en mécanique quantique",
    *C. R. Acad. Sci. Paris* **292**, Sér. II (1981) 1357–1362. Orthonormal bases only. **[P]**
13. F. Low, "Complete sets of wave packets", in *A Passion for Physics — Essays in Honor of
    Geoffrey Chew* (C. DeTar et al., eds.), World Scientific, 1985, pp. 17–22. **[P]**
14. G. Battle, "Heisenberg proof of the Balian–Low theorem", *Lett. Math. Phys.* **15** (1988)
    175–177. DOI 10.1007/BF00397840. **[V]**
15. I. Daubechies and A. J. E. M. Janssen, "Two theorems on lattice expansions", *IEEE Trans.
    Inform. Theory* **39**(1) (1993) 3–6. **[V]** venue/pages; **[U]** its exact role in the
    Riesz-basis extension.
16. J. J. Benedetto, C. Heil, D. F. Walnut, "Differentiation and the Balian–Low theorem",
    *J. Fourier Anal. Appl.* **1**(4) (1995) 355–402. DOI 10.1007/s00041-001-4016-5. Survey of
    record. **[V]**
17. I. Daubechies, S. Jaffard, J.-L. Journé, "A simple Wilson orthonormal basis with exponential
    decay", *SIAM J. Math. Anal.* **22**(2) (1991) 554–572. DOI 10.1137/0522035. **[P]**
18. K. Gröchenig, *Foundations of Time-Frequency Analysis*, Birkhäuser, 2001. Ch. 8, The Zak
    Transform (BLT). **[U]** theorem number.
19. J. Ascensi, H. G. Feichtinger, N. Kaiblinger, "Dilation of the Weyl symbol and Balian–Low
    theorem", *Trans. Amer. Math. Soc.* **366** (2014) 3865–3880. Irregular (non-lattice) Gabor.
    **[P]**
20. S. Nitzan and J.-F. Olsen, "A quantitative Balian–Low theorem", *J. Fourier Anal. Appl.*
    (2013). DOI 10.1007/s00041-013-9289-y; arXiv:1205.0163. **[P]**
21. S. Z. Gautam, "A critical-exponent Balian–Low theorem", *Math. Res. Lett.* **15** (2008)
    471–483. **[P]**
22. D. P. Hardin, M. C. Northington V., A. M. Powell, "A sharp Balian–Low uncertainty principle
    for shift-invariant spaces", *Appl. Comput. Harmon. Anal.* (2016).
    DOI 10.1016/j.acha.2016.05.001; arXiv:1510.04855. **[V]** (abstract read).

*Attribution note.* The extension of BLT from orthonormal to Riesz bases is credited in the
reference literature to **Daubechies, Coifman and Semmes** (with the Zak-transform machinery to
Daubechies–Janssen). The pairing "Coifman/Daubechies/**Jaffard/Journé**" that circulates for
this extension appears to be a conflation with item 17, which is the Wilson-basis paper and a
different result. **[P]**

### Shift-invariant spaces
23. C. de Boor, R. DeVore, A. Ron, "The structure of finitely generated shift-invariant spaces in
    `L_2(R^d)`", *J. Funct. Anal.* **119** (1994) 37–78. **[P]**
24. A. Ron and Z. Shen, "Frames and stable bases for shift-invariant subspaces of `L_2(R^d)`",
    *Canad. J. Math.* **47** (1995) 1051–1094. Gramian / dual-Gramian fiberization. **[V]** venue;
    **[U]** theorem numbers.
25. M. Bownik, "The structure of shift-invariant subspaces of `L^2(R^n)`", *J. Funct. Anal.*
    **177** (2000) 282–309. Range-function reformulation. **[P]**
26. A. Aldroubi and K. Gröchenig, "Nonuniform sampling and reconstruction in shift-invariant
    spaces", *SIAM Review* **43** (2001) 585–620. **[P]**

### Locality vs conditioning
27. S. Demko, W. F. Moss, P. W. Smith, "Decay rates for inverses of band matrices", *Math. Comp.*
    **43**(168) (1984) 491–499. **[V]**
28. M. Benzi and G. H. Golub, "Bounds for the entries of matrix functions with applications to
    preconditioning", *BIT* **39**(3) (1999) 417–438. DOI 10.1023/A:1022362401426. **[P]**
29. M. Benzi and N. Razouk, "Decay bounds and O(n) algorithms for approximating functions of
    sparse matrices", *Electron. Trans. Numer. Anal.* **28** (2007) 16–39. **[P]**
30. M. Benzi, P. Boito, N. Razouk, "Decay properties of spectral projectors with applications to
    electronic structure", *SIAM Review* **55**(1) (2013) 3–64; arXiv:1203.3953. **[V]**
    (abstract read; explicitly covers the non-orthogonal / overlap-matrix representation.)
31. S. Jaffard, "Propriétés des matrices « bien localisées » près de leur diagonale et quelques
    applications", *Ann. Inst. H. Poincaré Anal. Non Linéaire* **7**(5) (1990) 461–476. **[V]**
32. K. Gröchenig and M. Leinert, "Wiener's lemma for twisted convolution and Gabor frames",
    *J. Amer. Math. Soc.* **17** (2004) 1–18. **[P]**
33. Q. Sun, "Wiener's lemma for infinite matrices", *Trans. Amer. Math. Soc.* **359** (2007)
    3099–3123. **[V]**
34. K. Gröchenig and A. Klotz, "Norm-controlled inversion in smooth Banach algebras, I",
    *J. London Math. Soc.* **88** (2013) 49–64; "…II", *Math. Nachr.* **287** (2014) 917–937;
    arXiv:1207.1269, arXiv:1211.2974. **[P]**

### Topological no-go (localized AND orthonormal)
35. D. J. Thouless, "Wannier functions for magnetic sub-bands", *J. Phys. C: Solid State Phys.*
    **17** (1984) L325. DOI 10.1088/0022-3719/17/12/003. **[V]**
36. C. Brouder, G. Panati, M. Calandra, C. Mourougane, N. Marzari, "Exponential localization of
    Wannier functions in insulators", *Phys. Rev. Lett.* **98**, 046402 (2007).
    DOI 10.1103/PhysRevLett.98.046402; arXiv:cond-mat/0606726. **[V]**
37. G. Panati, "Triviality of Bloch and Bloch–Dirac bundles", *Ann. Henri Poincaré* **8** (2007)
    995–1011. DOI 10.1007/s00023-007-0326-8; arXiv:math-ph/0601034. **[V]**
38. D. Monaco, G. Panati, A. Pisante, S. Teufel, "Optimal decay of Wannier functions in Chern and
    quantum Hall insulators", *Comm. Math. Phys.* **359** (2018) 61–100.
    DOI 10.1007/s00220-017-3067-7; arXiv:1612.09552. **[V]** abstract; **[U]** exact exponent.

### Quantum-chemistry side
39. E. Artacho and L. Miláns del Bosch, "Nonorthogonal basis sets in quantum mechanics:
    Representations and second quantization", *Phys. Rev. A* **43**, 5770 (1991).
    DOI 10.1103/PhysRevA.43.5770. **[V]** citation; **[U]** whether the delocalization of the
    dual basis is stated there as a named theorem (full text paywalled).
40. P.-O. Löwdin, "On the non-orthogonality problem connected with the use of atomic wave
    functions in the theory of molecules and crystals", *J. Chem. Phys.* **18**(3) (1950)
    365–375. DOI 10.1063/1.1747632. **[V]**
41. P.-O. Löwdin, "Quantum theory of cohesive properties of solids", *Adv. Phys.* **5** (1956)
    1–171 — canonical orthogonalization. **[P]**
42. B. Klahn and W. A. Bingel, "The convergence of the Rayleigh–Ritz method in quantum
    chemistry", *Theor. Chim. Acta* **44** (1977) 9–26 (I) and 27–43 (II).
    DOI 10.1007/BF00548027. Relevant because it establishes that `L^2`-completeness of a basis is
    **not** sufficient for energy convergence. **[P]** (one secondary listing gives vol. 47 for
    part I; unresolved.)

*Not found, despite direct search:* any named theorem in the quantum-chemistry literature
asserting that symmetric orthogonalization of an overcomplete or near-linearly-dependent basis
*necessarily* delocalizes. The phenomenon is well documented as a numerical-stability folklore
fact (small overlap eigenvalues amplify error; canonical orthogonalization / pivoted Cholesky as
the standard remedy) but is not, as far as this scan reaches, a citable theorem. Proposition D
above is the citable substitute, and it is stronger.

---

## 8. One-paragraph summary for the paper

> The two-centre Coulomb–Sturmian basis is a complete set adjoined to a translate of itself.
> Such a union is never a Riesz sequence: for any Hilbert space, if one member of the second
> family lies in the closed span of the first, the smallest eigenvalue of the truncated Gram is
> bounded by the squared `N`-term approximation error of that member by the first family, and
> therefore tends to zero. In frame language (Christensen; Heil) the union is a frame that is
> inexact, hence not a Riesz basis, hence its synthesis operator has a non-trivial `l^2` kernel;
> the vanishing `lambda_min` is the finite-section shadow of that kernel. The rate is set by a
> symbol: translation by `R` acts on the Fock sphere as multiplication by the angular-averaged
> phase `j_0(pR)`, which becomes trivial at zero momentum, so the `2 x 2` Gram symbol
> `[[1, sigma],[sigma, 1]]` attains `sigma = 1` at `p = 0`; the maximum is quadratic, giving
> `1 - sigma_max ~ (kR)^2 pi^2 / 24 n^2` and the measured `cond ~ N^2`. The near-null direction
> is the ungerade combination, matching the measured parity confinement. Separately and
> unconditionally: no congruence that preserves the angular-momentum block decomposition can
> orthogonalize a metric that couples `l`, so the loss of Gaunt selection under Löwdin (or any
> other) orthogonalization is not a conditioning artefact and would occur at `kappa = 1`.
> Balian–Low, Beurling density and the Balan–Casazza–Heil–Landau density theorems do not apply:
> all three require a group action and a system at critical density, and this system is at
> redundancy two.
