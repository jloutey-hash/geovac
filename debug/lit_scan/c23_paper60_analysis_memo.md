# C23 inverse-citation scan — Paper 60, harmonic-analysis / asymptotics half

**Date:** 2026-09-12. **Criterion:** `docs/qa/criteria.md` C23 ("is a result the paper
presents as its OWN derivation already a named theorem?"). **Scope:** claims B1–B4
(chirp Fourier asymptotics; the `j_0`-of-cotangent symbol; the Wiener-lemma /
off-diagonal-decay mechanism; the `M`-centre all-ones degeneracy). Read-only; nothing
under `papers/`, `geovac/`, `tests/`, `docs/` was touched.

**Verification discipline.** Every citation below is marked `[READ]` (I read the
statement in a source I can name and quote) or `[BIB-ONLY]` (bibliographic data
verified in a reference list I read, content NOT read). Per C23's hard rule, no
theorem is stated here that I have not seen stated somewhere I can name. Search-engine
summaries were not treated as verification. **Two independent re-derivations were done
locally** (B1's asymptotic, B1's exact Bessel evaluation) rather than taken on report.

---

## Consolidated verdict

| Claim | Verdict | One line |
|---|---|---|
| **B1** | **PRIOR ART on the derivation chain; ABSENT as a statement** | The `j^-5/4` law is not a named theorem for this function, but the model integral is *exactly* a modified Bessel function, so exponent, constant, `sqrt(j)` phase and the `pi/4` are all DLMF 10.32.10 + 10.40.2. The paper's "stationary phase" re-derives a closed form. |
| **B2** | **ABSENT** | No work composes `j_0` with a cotangent / Möbius reparametrisation and asks for coefficient asymptotics. One false friend worth knowing: "Toeplitz operators on the Fock space" in the OA literature means *Bargmann*–Fock. |
| **B3** | **PRIOR ART — and GeoVac's reading of the hypothesis is correct** | Jaffard (1990): polynomial decay `r > d` **plus invertibility on `l^2`**. Our `r = 5/4 > 1 = d` is comfortably inside the class; the obstruction to a localized `S^{-1/2}` is `0 in spec`, exactly as claimed. Three sharpenings owed. |
| **B4** | **PRIOR ART (different literature)** | This is the "flat limit" of kernel matrices: `K = K(0,0) 1 1^T + O(eps)`, characteristic root of multiplicity `n-1`, point set fixed and arbitrary. Barthelmé–Usevich state it verbatim. Not found in the basis-set linear-dependence literature. |

**Decision-gate reading.** B1's load-bearing outcome is *not* a clean "PRIOR ART" that
moves a headline result to "an instance of" — nobody has stated this law for this
function. But the finding is arguably more useful than a prior-art hit would have been:
the asymptotic did not need deriving at all, and the closed form hands the paper a
sharper claim than it currently makes. B3 is a clean PASS on the mechanism claim, with
the hypothesis correctly identified.

---

## B1 — the `j^-5/4` chirp law

### The claim, restated

`a(chi) = j_0(c cot(chi/2))` on `(0, pi)`; near `chi -> 0`, `cot(chi/2) ~ 2/chi`, so
`a(chi) ~ (chi/2c) sin(2c/chi)` — amplitude linear, phase `2c/chi`. Paper 60 asserts
(L838–839) only `|c_j| ~ j^-5/4`; the sharper form under test is

    |c_j| = (2 pi)^-1/2 * 2^-3/4 * c^-1/4 * j^-5/4 * |sin(2 sqrt(2 c j) + pi/4)| + o(j^-5/4).

### (i) Is the `j^-5/4` envelope a named result? **Not found.**

Searched, and not found stated for this class: Fourier asymptotics of chirps;
oscillatory integrals with a stationary point escaping to the boundary; stationary
phase applied to `sin(1/x)`-type singularities; Fourier coefficients of functions with
an oscillatory endpoint singularity. Two specific misses worth recording, because they
are where one would expect a hit:

- **Darboux-method coefficient asymptotics.** Olver, *Asymptotics and Special
  Functions*, has a section titled "Behavior of the Coefficients at a Singularity",
  pp. 386–388 in the 1997 A K Peters reprint `[BIB-ONLY — section title and page range
  verified via Crossref chapter DOIs 10.1201/9781439864548-123 and -124; text not
  read]`. Darboux's method keys off *algebraic and logarithmic* singularities. Ours is
  an essential, infinitely-oscillating singularity, which that machinery does not cover
  — so the natural textbook home for "coefficients from the singularity" is the wrong
  one here, which is plausibly why no statement exists.
- **Van der Corput / amplitude-singularity refinements.** The modern literature on the
  "interaction of amplitude singularities with stationary points" treats an *integrable
  amplitude* singularity at an endpoint, not an essential *phase* singularity with an
  escaping stationary point. Searched; no match.

### (ii) The real finding: the model integral is a Bessel function, exactly. `[READ]`

The stationary-phase computation is unnecessary. Split
`cos(j chi) sin(2c/chi) = 1/2 [ sin(2c/chi + j chi) + sin(2c/chi - j chi) ]`; only the
first term has a stationary point. The resulting model integral has a closed form:

    I(j,c) := int_0^inf x exp(i(2c/x + j x)) dx = 2 (P/Q) K_2( 2 sqrt(P Q) ),
    P = -2 i c,  Q = -i j,  sqrt(PQ) = -i sqrt(2 c j).

This is **DLMF 10.32.10** `[READ — quoted verbatim from dlmf.nist.gov/10.32]`:

> `K_nu(z) = (1/2) (z/2)^nu int_0^inf exp(-t - z^2/(4t)) dt / t^(nu+1)`,  `|ph z| < pi/4`

with `t = Q x` and `nu = -2`. (Algebra: substituting `t = q x` and `p = z^2/(4q)` turns
10.32.10 into `int_0^inf x^(-nu-1) exp(-p/x - q x) dx = 2 (q/p)^(nu/2) K_nu(2 sqrt(pq))`;
`-nu-1 = 1` gives `nu = -2`, `K_-2 = K_2`, and `(q/p)^-1 = p/q`.) The identity is stated
for `|ph z| < pi/4` while our `z` sits at `ph z = -pi/2`, so it is used by analytic
continuation — **verified numerically here** (table below), including the branch.

Then **DLMF 10.40.2** `[READ — quoted verbatim from dlmf.nist.gov/10.40]`:

> `K_nu(z) ~ (pi/2z)^(1/2) e^(-z) sum_k a_k(nu)/z^k`,  `|ph z| <= 3 pi/2 - delta`

covers `ph z = -pi/2` **directly** — no continuation needed on this leg. Feeding
`z = -2 i sqrt(2 c j)` through it, and carrying the `1/pi`, the `1/2` from the
sine-product split and the `1/(2c)` amplitude factor, reproduces
`(2 pi)^-1/2 2^-3/4 c^-1/4 j^-5/4` **and** the phase `2 sqrt(2 c j) + pi/4` — the `pi/4`
being the standard `sqrt(pi/2z)` branch phase, not a stationary-phase signature. I did
this algebra independently and it agrees with the claimed constant exactly.

**So every element of B1 — exponent, constant, sqrt-`j` phase, `pi/4` — is a corollary
of two DLMF entries.** The paper derived, correctly, something that has a closed form.

### (iii) The class IS named, though the general law is UNVERIFIABLE from here

Coefficient sequences of the shape `j^-beta e^(i c j^alpha)` are the classical
**oscillating-multiplier / "special trigonometric series"** family; ours is exactly
`alpha = 1/2, beta = 5/4`. Attributions read from the body text and reference list of
*Oscillating singular integral operators on compact Lie groups revisited*
(arXiv:2202.10531), whose introduction says the `L^p` properties of symbols
`(1+|xi|^2)^(-n theta/4) e^(i (1+|xi|^2)^(theta/2))` "were firstly studied by Hardy
[32], Hirschman [34] and Wainger [46]" `[READ for that sentence and for the three
bibliographic entries; the three originals NOT read]`:

- G. H. Hardy, *A theorem concerning Taylor's series*, Quart. J. Pure Appl. Math. **44**
  (1913) 147–160.
- I. I. Hirschman, *Multiplier transformations I*, Duke Math. J. (1956) 222–242.
  (No volume number was printed in the entry I read; it is commonly cited as vol. 26 —
  **do not print a volume number without checking it**.)
- S. Wainger, *Special trigonometric series in k-dimensions*, Mem. Amer. Math. Soc.
  **59** (1965).
- Later in the same lineage: Hörmander, C. Fefferman, Fefferman–Stein, Sjölin, Miyachi,
  Peral.

The standard duality in this family — phase `|x|^alpha` in one variable pairing with
phase `|xi|^(alpha/(alpha-1))` in the other — gives `alpha = -1 -> 1/2`, i.e. exactly
our `sqrt(j)`. **I could not reach Wainger's memoir**, so whether it contains the
general amplitude-exponent law `|c_j| ~ j^(-(2a+3)/4)` for amplitude `chi^a` is
**UNVERIFIABLE** from here. That general law (specialising to `5/4` at `a = 1`) is also
not something I found stated anywhere.

### Independent numerical re-derivation (done here, not taken on report)

`b1_fast2.py` — two quadrature resolutions (16-pt vs 24-pt Gauss–Legendre; 8 vs 13
panels per `cos(j chi)` period), panels on the union of the chirp zeros
`chi_k = 2c/(k pi)` (`k <= 4e5`) and a uniform grid:

| `c` | `j` | numeric `c_j` | asymptotic | ratio |
|---:|---:|---:|---:|---:|
| 1 | 64 | -1.3145e-3 | -1.2959e-3 | 1.0144 |
| 1 | 256 | 2.0212e-4 | 2.0471e-4 | 0.9873 |
| 1 | 1024 | -8.2295e-6 | -7.6889e-6 | 1.0703 |
| 1 | 4096 | -2.8258e-6 | -2.8699e-6 | 0.9846 |
| 1 | 16384 | -1.27930e-6 | -1.27914e-6 | 1.0001 |
| 2 | 64 | 1.0878e-3 | 1.0797e-3 | 1.0075 |
| 2 | 4096 | -4.4859e-6 | -4.4723e-6 | 1.0031 |
| 2 | 16384 | -6.9677e-7 | -6.9801e-7 | 0.9982 |

(The one outlier in the full run, `c=2, j=1024`, ratio 1.31, sits where
`sin(2 sqrt(2cj) + pi/4)` is near a zero and the ratio is meaningless.) **Sign and phase
reproduce, not just the envelope.**

`bessel_id.py` — the Bessel identity, quadrature vs `2(P/Q) K_2(2 sqrt(PQ))` on the
`-i` branch, cutoff-averaged over 8 half-periods:

| `c` | `j` | quad | `2(P/Q) K_2` | rel. diff |
|---:|---:|---|---|---:|
| 1 | 50 | -0.0105129+0.0198668i | -0.0099515+0.0201491i | 2.8e-2 |
| 1 | 200 | -0.0040031+0.0000221i | -0.0039655+0.0000335i | 9.8e-3 |
| 1 | 800 | 0.0004482-0.0005358i | 0.0004505-0.0005367i | 3.5e-3 |
| 2 | 800 | 0.0007828+0.0008798i | 0.0007850+0.0008788i | 2.1e-3 |

The residual falls like `1/j` — it is the quadrature's own finite-cutoff error on a
conditionally convergent integral, not a mismatch. The `+i` branch is rejected
(rel. diff ~1.5, i.e. the conjugate).

### What Paper 60 should do

1. **L838–839: cite, don't derive.** Replace "the symbol is a chirp … whose Fourier
   coefficients decay only as `j^-5/4`" with the closed-form chain: the chirp integral
   is a modified Bessel function (DLMF 10.32.10 at `nu = -2`, by continuation), and
   `|c_j| ~ (2 pi)^-1/2 2^-3/4 (kR)^-1/4 j^-5/4 |sin(2 sqrt(2 kR j) + pi/4)|` is its
   large-argument asymptotic (DLMF 10.40.2). This is an **upgrade**: the paper currently
   claims an order of magnitude and can claim the constant and the phase, cited.
2. **Add one orientation sentence** on class membership (Hardy / Hirschman / Wainger
   oscillating multipliers, `alpha = 1/2`, `beta = 5/4`) — with the honest cap that the
   general amplitude law was not located.
3. **Do not** print Hirschman's volume number without checking it.

---

## B2 — has anyone studied `j_0` composed with a cotangent?

**Verdict: ABSENT.**

Searched: momentum-space quantum chemistry (Shibuya–Wulfman lineage, Avery, the
Aquilanti–Cavalli–Coletti–Caligiana lineage), hyperspherical harmonics as momentum-space
orbitals, and the Toeplitz-symbol literature. Crossref bibliographic sweeps returned the
expected canon — Avery, *Many-center Coulomb Sturmians and Shibuya–Wulfman integrals*,
Int. J. Quantum Chem. **100** (2003) 121–130, DOI 10.1002/qua.10820 `[BIB-ONLY]`;
Avery, Hansen, Wang & Antonsen, *Sturmian basis sets in momentum space*, Int. J. Quantum
Chem. **57** (1996) 401–411 `[BIB-ONLY]`; Aquilanti, Cavalli & Coletti, *The
d-dimensional hydrogen atom…*, Chem. Phys. **214** (1997) 1–13 `[BIB-ONLY]` — and
**nothing** that treats these overlaps as a *multiplication operator with a symbol*, or
asks for coefficient asymptotics in the Fock chart. Independent of, and consistent with,
the 2026-09-11 scan's Q4 ("the symbol has never appeared; the *operator* is
Shibuya–Wulfman").

**One false friend worth putting on the record.** A search for "Toeplitz" + "Fock"
returns a substantial literature — *Toeplitz operators on the Fock space*,
hyponormality, commutative `C*`-algebras generated by them. That is the **Bargmann–Fock
space** of entire functions, unrelated to Fock's 1935 conformal projection. Anyone
re-running this scan will hit those and must not read them as coverage.

**Caveat (unchanged from the earlier scan):** the two Avery books were not reachable in
full text. The verdict is "not found in the reachable literature", not "proven absent".

---

## B3 — the Wiener-lemma / off-diagonal-decay mechanism

**Verdict: PRIOR ART, and GeoVac's reading of the invertibility hypothesis is correct
and standard.**

### The theorem, verbatim `[READ]`

From Gröchenig & Klotz, *Noncommutative approximation: inverse-closed subalgebras and
off-diagonal decay of matrices*, arXiv:0904.0386, Introduction (their Eq. (1.1)):

> "Jaffard's theorem [33]: If the matrix `A` with entries `A(k,l)`, `k,l in Z`, is
> **invertible on `l^2(Z)`** and if, for `r > 1`, `|A(k,l)| <= C(1+|k-l|)^-r` for all
> `k,l in Z`, then also `|(A^-1)(k,l)| <= C(1+|k-l|)^-r` for all `k,l in Z`."

and their Eq. (2.8) `[READ]`: "the **Jaffard algebra** `J_r`, `r > d`, is defined by the
norm `||A||_{J_r} = sup_{k,l in Z^d} |A(k,l)| v_r(k-l)`."

Original: S. Jaffard, *Propriétés des matrices « bien localisées » près de leur diagonale
et quelques applications*, Ann. Inst. H. Poincaré Anal. Non Linéaire **7**(5) (1990)
461–476 `[BIB-ONLY — entry read in Gröchenig–Klotz ref [33] and in the EMS Press
listing; the paper itself not read]`.

### Answering the question exactly as posed

> *for a matrix algebra with polynomial off-diagonal decay of order `> 1`, is
> inverse-closedness known, and what exactly is the invertibility hypothesis?*

**Yes, and the hypothesis is invertibility as an operator on `l^2`.** In 1D (`d = 1`)
the threshold is exactly `r > 1`. **`5/4 > 1`, so Paper 60's matrix is comfortably
inside the Jaffard class** — the margin is not the issue and should not be framed as
one.

A small observation the paper can use: the object is Toeplitz-*minus*-Hankel,
`M_{nm} = c_{n-m} - c_{n+m}`, and the Hankel part does **not** break the class, because
`n + m >= |n - m|` for `n, m >= 1` gives `|c_{n+m}| <= C(1+|n-m|)^-5/4`; hence
`|M_{nm}| <= 2C(1+|n-m|)^-5/4` and `M in J_{5/4}`.

### The `S^-1/2` leg `[READ]`

Inverse-closedness is exactly what licenses the functional calculus. Gröchenig,
*Wiener's Lemma: Theme and Variations* (Inzell Lectures, slide 23) states:

> "**Corollary.** If `A` is inverse-closed in `B`, then the Riesz functional calculi for
> `A` and `B` coincide." — with the listed examples "square roots, powers,
> pseudoinverse; Theorem of Wiener–Lévy".

`f(a) = (1/2 pi i) int_Gamma f(z) (z e - a)^-1 dz` requires `f` analytic on an open
neighbourhood of `sigma_B(a)`. For `f(z) = z^-1/2` that means **`0 not in spec`**.
Published form: K. Gröchenig, *Wiener's Lemma: Theme and Variations. An Introduction to
Spectral Invariance*, Applied and Numerical Harmonic Analysis, Birkhäuser, Boston, 2009
(Inzell Lectures on Harmonic Analysis) `[BIB-ONLY — entry as printed in Gröchenig–Klotz
ref [27]]`. The matrix-specific companion is Gröchenig & Leinert, *Symmetry and
inverse-closedness of matrix algebras and functional calculus for infinite matrices*,
Trans. Amer. Math. Soc. **358**(6) (2006) 2695–2711 `[BIB-ONLY — entry read in
Gröchenig–Klotz ref [28]]`.

**So GeoVac's mechanism claim is right:** the inverse square root escapes the
localized-inverse theorems because the symbol touches zero (`0 in spec`, so `z^-1/2` is
not holomorphic on a neighbourhood of the spectrum), **not** because `5/4` is too slow a
decay. The decay class is fine.

### Three sharpenings owed

1. **Citation correction.** The prompt's "Gröchenig–Leinert, JAMS 17 (2004) 1" is a real
   paper — *Wiener's lemma for twisted convolution and Gabor frames* `[BIB-ONLY — entry
   read in the Inzell slides reference list]` — but it is the twisted-convolution /
   Gabor result, **not** the matrix-algebra one. For matrices cite Jaffard (1990) and
   Gröchenig–Leinert TAMS **358** (2006).
2. **The quantitative form is the one that actually bites.** Every finite section is
   trivially invertible; what degrades is the *constant*. The relevant refinement is
   **norm-controlled inversion**, where the decay constant of `A^-1` is controlled by
   `||A^-1||_{B(l^2)}`: Gröchenig & Klotz, *Norm-controlled inversion in smooth Banach
   algebras I / II* (II = arXiv:1211.2974, Math. Nachr. 2014) and the Shin–Sun line
   `[BIB-ONLY — titles/venues seen in listings; papers NOT read, so do not quote a
   theorem from them without reading]`. If Paper 60 wants "escapes only because it is
   not boundedly invertible" to carry weight, this is the family that makes it
   quantitative. There is also a genuine gap to flag honestly: Jaffard's theorem is
   about an *infinite* matrix invertible on `l^2`, whereas Paper 60's object is a
   *sequence of finite sections*; the bridge is the finite-section / limit-operator
   theory, and this scan did not verify that bridge.
3. **Do not conflate two different borderline conditions.** `sum_j |c_j| < inf`
   (`5/4 > 1`) holds — the symbol **is** in the Wiener algebra. What fails is
   `sum_j j |c_j| < inf`, the Böttcher–Widom smoothness hypothesis, a different
   condition serving a different theorem. Paper 60's L836–843 paragraph already keeps
   these apart; keep it that way.

---

## B4 — the `M`-centre all-ones degeneracy

**Verdict: PRIOR ART, in the kernel / radial-basis-function literature, not in
quantum chemistry.**

The mathematical content — all pairwise couplings tend to a common value, so the matrix
degenerates to `phi(0) * J_M` of rank one, the null space has dimension `M-1`, and the
structure does not depend on where the centres are — is the **"flat limit"** of kernel
matrices, and it is stated explicitly. From Barthelmé & Usevich, *Spectral properties of
kernel matrices in the flat limit*, SIAM J. Matrix Anal. Appl. **42**(1) (2021) 17–57,
arXiv:1910.14067, Introduction `[READ — quoted from the arXiv PDF]`:

> "The difficulty comes from the fact that `K = K(0,0) 1 1^T + O(eps)`, where
> `1 = [1 … 1]^T`, i.e., we are dealing with a **singular perturbation problem**."

and their footnote 2 `[READ]`:

> "Seen from the point of view of the characteristic polynomial, the equation
> `det(K - lambda I) = 0` has a solution of **multiplicity `n-1`** at `eps = 0`, but
> these roots immediately separate when `eps > 0`."

Geometry-independence is explicit in their framing `[READ]`: "The point set `X` is
considered fixed, with arbitrary geometry (i.e., not lying in general on a regular
grid)". They also note that "it is impossible to retrieve the information about the
limiting projectors just from `K(0)` (which is **rank-one**)" — i.e. the rank-one limit
is their starting point, and the content is *which directions* the `M-1` eigenvalues
separate along, which is the part Paper 60 does not yet address.

Lineage `[BIB-ONLY — all four entries read in Barthelmé–Usevich's reference list; the
papers themselves not read]`:

- T. A. Driscoll & B. Fornberg, *Interpolation in the limit of increasingly flat radial
  basis functions*, Comput. Math. Appl. **43** (2002) 413–422,
  DOI 10.1016/S0898-1221(01)00295-4 — coined "flat limit".
- R. Schaback, *Multivariate interpolation by polynomials and radial basis functions*,
  Constr. Approx. **21** (2005) 293–317 — Theorem 6, orders of the eigenvalues.
- A. J. Wathen & S. Zhu, *On spectral distribution of kernel matrices related to radial
  basis functions*, Numer. Algorithms **70** (2015) 709–726.
- B. Fornberg, G. Wright & E. Larsson, *Some observations regarding interpolants in the
  limit of flat radial basis functions*, Comput. Math. Appl. **47** (2004) 37–55.

**Not found in the quantum-chemistry basis-set literature.** Crossref sweeps over
basis-set linear dependence, Löwdin canonical orthogonalization, quasi-linear dependence
and overlap-matrix spectra returned the expected canon (Naidu & Srivastava, IJQC **99**
(2004) 882–888; Aissing & Monkhorst, IJQC **43** (1992) 733–745; Jiao & Ho, IJQC **115**
(2015) 434–441; Neymeyr & Engel, IJQC **53** (1995) 537–540 — the last being the
*diatomic*, i.e. `M = 2`, overlap-spectrum case) `[all BIB-ONLY]`, and none states the
`M`-centre `J_M` / rank-`M-1` / fixed-direction observation.

**Reading for Paper 60.** Not new mathematics; the honest line is the one C23 prescribes
for `eq:sigma_law` — cite Barthelmé–Usevich (and Driscoll–Fornberg for the flat-limit
concept) and claim **the identification**: that the two-centre Sturmian block symbols
realise this structure at a *symbol point* (`chi = pi`, i.e. `p = 0`) rather than in a
shape-parameter limit, and that this is what licenses the geometry-independent
preconditioner. That identification is genuine and is where the value is. One structural
difference worth stating rather than hiding: in the RBF flat limit the degeneracy is a
global limit of the whole matrix, whereas here it is local to one end of the symbol with
the rest of the spectrum unaffected — which is precisely why a rank-`M-1` rotation fixes
it and a global rescaling does not.

---

## Files

- This memo. No files created or modified outside it and the session scratchpad.
- Scratchpad drivers (not in the repo): `b1_fast2.py` (B1 asymptotic, two quadrature
  resolutions), `bessel_id.py` (the `K_2` closed form, branch resolved).
