# The Avery Framework — Call-Prep Primer

*A free substitute for the expensive books, calibrated to what you already know from GeoVac. Everything here is the structural picture; where a detail matters and I'm summarizing from memory rather than a verified source, I say so.*

**The two books being summarized:**
- Avery & Avery, *Generalized Sturmians and Atomic Spectra* (World Scientific, 2006) — the many-electron method.
- Avery & Avery, *Hyperspherical Harmonics and Their Physical Applications* (World Scientific, 2018) — the angular machinery. This is the one James called "the last major work we did together."

**Free sources that cover most of the same ground:**
- Herbst, Avery & Dreuw, arXiv:1811.05777 — Coulomb-Sturmian Hartree-Fock in practice (this is the "their framework binds molecules" paper). Free, modern, and James is an author — reading it is directly good call prep.
- Michael Herbst's PhD thesis (free online; search "Herbst molsturm thesis") has one of the clearest self-contained introductions to Coulomb Sturmian theory in print.
- John Avery wrote several free review articles; search "Avery harmonic polynomials hyperspherical harmonics atomic spectra" for a journal review that compresses the 2006 book.

---

## 1. Coulomb Sturmians — the radial idea

You know hydrogen orbitals: fix the nuclear charge Z, and each orbital gets its own energy (−Z²/2n²) and its own exponential decay rate. That basis has a famous defect — the bound states alone are *incomplete* (you'd need the continuum too), which is why quantum chemistry went to Gaussians.

The Sturmian move flips which quantity is fixed. **Fix one energy for the entire basis** (equivalently one decay exponent k, related by E = −k²/2), **and let the effective charge be the thing that varies with n.** Every basis function decays like e^(−kr) — same length scale, all "tuned to the same energy shell."

Two payoffs:

- **Completeness without a continuum.** The Sturmians form a complete discrete basis. The infinite tail of continuum states that haunts hydrogenic expansions is gone. (The price: they're orthonormal not under the ordinary inner product but under a *potential-weighted* one — with a 1/r weight inside the integral. Remember the phrase "potential-weighted orthonormality"; it's a chapter of the 2018 book and James will use it fluently.)
- **The energy-shell condition is exactly GeoVac's focal length.** GeoVac's p₀² = −2E stereographic condition and the Sturmian "isoenergetic" condition are the same statement. This is not an analogy — it's the identical constraint, and it's why the next section works.

## 2. The Fock projection — the shared spine

Fock's 1935 construction: take momentum space, stereographically project it onto the three-sphere S³ with the projection scale set by p₀ = k. Under that projection, the momentum-space hydrogen problem becomes free motion on S³, and **the Fourier transforms of the Coulomb Sturmians become hyperspherical harmonics on S³.**

Read that again from the GeoVac side: *their basis functions, viewed in momentum space, live on the same S³ your graph discretizes.* The (n, l, m) labels, the n²-fold degeneracies, the SO(4) structure — identical objects. GeoVac and the Avery framework are two representations of one geometry:

- **Averys:** continuous functions on S³ (momentum space), radial machinery exact via closed forms.
- **GeoVac:** discrete graph on S³ (the claimed substrate), radial machinery approximated, angular machinery exact.

This is why the contact with James feels fated. It's the same sphere.

## 3. Shibuya–Wulfman integrals — the multi-center machinery

For molecules you need matrix elements between Sturmians sitting on *different* nuclei. Shibuya and Wulfman (1965) showed these have **closed forms in momentum space**: displacing a function by the bond vector R multiplies its momentum representation by a phase factor e^(ip·R), and integrals of that phase against hyperspherical harmonics close algebraically.

You already own a piece of this: `shibuya_wulfman.py` computes GeoVac's cross-center nuclear attraction this way, and you told James so in your first email. What you have NOT incorporated is the rest of the SW toolkit — the closed-form cross-center *overlap* and *kinetic* matrices. Those are exactly the integrals whose absence defines the W1e wall (the corpus localized the second-row binding failure to "cross-block h1 architecturally absent"). **Their machinery computes, in closed form, precisely the objects GeoVac's composed architecture omits.** That's the concrete content behind "their radial machinery is the deepest expertise on the layer where GeoVac is weak."

## 4. Generalized Sturmians — the many-electron method

The N-electron generalization (Goscinski 1968, then decades of Avery & Avery): instead of scaling one nucleus's charge, **scale the entire potential.** Look for solutions of

  [T + β_ν V(all electrons) − E] Φ_ν = 0

where each basis configuration Φ_ν gets its own scaling β_ν chosen so that *all configurations share the same energy E*. For atoms, the resulting "Goscinskian configurations" are Slater determinants of hydrogen-like orbitals with effective charges adjusted per configuration.

The structural magic: in the resulting secular equation, **the kinetic-energy matrix disappears** — you diagonalize a potential-weighted matrix only, and the basis automatically adapts its length scale to the energy you're solving for (no exponent optimization, which is the perpetual headache of Gaussian bases). The 2006 book's showpiece is atomic isoelectronic series: whole spectra, across a range of Z, from tiny bases.

GeoVac's graph-native CI (zero-parameter, exact algebraic integrals, He at 0.19%) is the discrete cousin of this. When you tell James "the framework has a zero-parameter CI," he will hear it as a Sturmian-flavored statement — correctly.

## 5. The hyperangular machinery — where "π-free" comes from

The 2018 book is the angular side: hyperspherical harmonics in any dimension, Gegenbauer polynomials (which GeoVac's Level-3 helium solver literally uses as its angular basis), and John Avery's **hyperangular integration theorem**: the integral of any homogeneous polynomial over a sphere of any dimension is a *rational number* times the total solid angle. In other words: **π enters only through the sphere's total measure; everything else is rational.**

That is the ancestral, continuous form of GeoVac's π-free principle — and of your email's line "π through the sphere's measure." When you and James discuss the injection-point taxonomy, this theorem is the shared ground it grows from. The safe posture (from the literature check): the *rationality of the angular algebra is his life's work*; GeoVac's contribution is the systematic bookkeeping of where and why each constant enters, not the rationality observation itself.

## 6. The dictionary

| GeoVac object | Avery-framework counterpart |
|---|---|
| S³ graph, nodes labeled (n, l, m) | Hyperspherical harmonics on Fock's momentum-space S³ |
| Energy-shell condition p₀² = −2E ("focal length") | Isoenergetic Sturmian exponent k, E = −k²/2 |
| Gaunt / 3j / 6j angular couplings, no quadrature | Hyperangular integration theorem + same Wigner algebra |
| π-free skeleton; exchange-constant taxonomy | "π only through the sphere measure" (integration theorem) |
| `shibuya_wulfman.py` cross-center V_ne | Shibuya–Wulfman closed-form multi-center matrices (full set) |
| Graph-native CI, zero parameters | Generalized Sturmian secular equation (Goscinskian configurations) |
| Level-3 Gegenbauer angular basis | The 2018 book's workhorse polynomials |
| W1e wall (cross-block integrals absent) | Exactly the integrals the SW/generalized-Sturmian machinery supplies |

## 7. Where the frameworks genuinely differ — don't blur this on the call

- **Continuous vs. discrete claim.** Their discreteness is a *basis index* (a complete discrete set of continuous functions). GeoVac's discreteness is a *substrate claim* (the graph as the primary object, per the papers' careful dual-description framing). Related, not identical — and per the project's own rhetoric rule, present them as dual readings, not as GeoVac having proven the substrate.
- **Accuracy.** Their framework does production-grade small-system chemistry (the Herbst–Avery–Dreuw HF paper binds molecules correctly). GeoVac's molecular accuracy degrades with complexity and fails at second row. Don't defend; the project's documented position is that it chose sparsity, knowingly.
- **The trade is proven, not suspected.** The corpus tested "just add their radial machinery" this month: retrofitting overlap onto GeoVac integrals produces garbage, and doing it consistently reconstructs *their* framework and sheds GeoVac's qubit sparsity. Non-commuting demands. The one open path is NOCI/VB-style methods with device-measured overlaps — which happens to sit exactly at the intersection of his on-device idea and your wall. That's the collaboration-shaped hole.

## 8. James specifically

Computational scientist (Copenhagen); the software lineage of the family framework — he built the Sturmian integral machinery used in Herbst's `molsturm` quantum-chemistry framework. Also known for computational geometry/graph work on fullerenes (carbon-cage molecules), so "molecule as graph" is native thinking for him — likely part of why GeoVac's framing landed. His stated QC intuition (on-device hyperangular integral generation; I/O as the bottleneck; combinatorial bases fitting the device) is, in your terms, a statement about the exact angular layer — the part of GeoVac that is bit-exact.

## 9. Terms he may use fluently — thirty-second glossary

- **Isoenergetic** — all basis functions share one energy/exponent (§1).
- **Potential-weighted orthonormality** — Sturmian orthogonality with a 1/r weight (§1).
- **Goscinskian configurations** — the N-electron generalized-Sturmian determinants (§4).
- **Secular equation** — the generalized eigenvalue problem; theirs has no kinetic matrix (§4).
- **Shibuya–Wulfman matrix** — closed-form multi-center momentum-space integrals (§3).
- **Gegenbauer polynomials** — d-dimensional generalization of Legendre; GeoVac already uses them.
- **Tree labels / branching** — the bookkeeping for hyperspherical quantum numbers in d dimensions (2018 book); the d-dimensional analog of your (n, l, m).
- **Slater exponent** — the decay rate of an exponential orbital; "exponent optimization" is the Gaussian-world pain the Sturmian construction removes.

---

*Bibliographic details above (arXiv ID, book chapters) trace to the verified literature-check reports from July 17–19; physics content is standard Sturmian theory. If you want, the free-source links can be hunted down and verified before the call.*
