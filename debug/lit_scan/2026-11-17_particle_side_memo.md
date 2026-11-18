# Prior-art scan, particle side: records AB..FXK (v6.7.26-v6.7.36)

Date: 2026-11-17. Scope: the bond-angle operators, contact exponents, constant classes and tree rung results of
CHANGELOG v6.7.26-v6.7.36, `docs/chemistry_period_inventory.md` rows 4-8, Paper 63 `prop:tree_bound`,
`prop:tree_rung`, `prop:tree_rung_l`, `sec:four_centre`. Read-only; no code, papers or git touched.
Already credited and not re-reported: Broadhurst (kites 2212.01962, walks), BBBG 2008, Vanhove 2014, Bloch-Vanhove,
Brown, Fresán-Sabbah-Yu (quadratic relations), Remiddi 1991, Fromm-Hill 1987, Nickel 1978, Davydychev-Delbourgo,
Almkvist-Zudilin, AESZ, Bönisch, Zhou (Wick rotations, Wronskian; Paper 59 only), Bouferguene-Jones, Sack, Darboux /
Bernstein ellipse (water records), Hannesdottir-McLeod-Schwartz-Vergu 2024 non-simple pinches (water_pinch memo).

Key: **[V]** = abstract page or full text opened this session; **[V-text]** = PDF text extracted and the passage read;
**[V-sec]** = confirmed only through a citing source or catalogue page; **[UNV]** = not confirmed.

---

## Q1. Exponents at the contact point x*

**What the singular point is.** x* = (D1²+D2²)/(2D1D2) is where the complexified squared distance |R_B − R_C|²
vanishes: a complex-null coincidence of two *external* points of a position-space graph. It is not a momentum-space
threshold in a Mandelstam variable, so the textbook momentum-space formula does not apply as is.

**The standard formula.** Hannesdottir-Mizera, *What is the iε for the S-matrix?*, arXiv:2204.02988 [V-text], Sec. 7.1:
for an isolated, non-degenerate saddle (simple pinch), I ∝ φ^γ (log at γ = 0) with Eq. (7.12)
γ = (E−1)/2 − d = (LD − E − 2d_N − 1)/2, where d = E − LD/2 + d_N (Eq. 7.11 context). Sec. 7.2 then bounds γ from below
in local QFT ("only simple poles, square-root singularities and logarithmic branch cuts"). This is the
Landau / Polkinghorne-Screaton / Eden et al. exponent in modern form. With general propagator powers ν_e the same
saddle count gives γ = (E−1)/2 − (Σν_e − LD/2 + d_N) [our generalisation; standard but not quoted from the source].
Note for the dictionary: a 1s density e^{−ar} is −∂_a of the Yukawa propagator e^{−ar}/r, i.e. ν = 2 in momentum space.

**Does it give 7/2?** Not directly: our pinch is not isolated. At R = B − C complex null, the points with
|r−B|² = |r−C|² = 0 and parallel normals form the whole complex null line through B and C (Morse-Bott, critical set of
dimension k = 1). That is the same non-simple-pinch mechanism the water_pinch memo already credits to Hannesdottir et
al. 2024 App. B. Our own evaluation (heuristic, [OUR ARGUMENT], not found in print): for an integrand Π Q_i^{α_i}
(Q_i squared distances; α = +1/2 per cusp |r−X|, −1/2 per Coulomb 1/|r−X|) in n = 3 dimensions, two pinching factors
and a k = 1 critical line, γ = Σα + (n − k + E − 1)/2 = Σα + 3/2. Cross-check by real-space Fourier counting
(Gel'fand-Shilov: FT|x|^λ ∝ Γ((λ+3)/2)/Γ(−λ/2) |k|^{−λ−3}), which gives the same powers of |R|:

| Case | Σα | γ (count) | Corpus |
|:--|:--|:--|:--|
| one-electron 3-centre, cusp × Coulomb | 0 | 3/2 (|R|³) | 3/2, 5/2 (row 4) ✓ |
| (AA\|BC), cusp × cusp, smooth A potential | 1 | 5/2 (|R|⁵) | **no 5/2**; 7/2 leads |
| next cusp orders u³v, uv³ | 2 | 7/2 | 7/2 ✓ |
| u³v³ / p,d numerators | 3 | 9/2 | 9/2 (p, d at equal arms) ✓ |

So the counting reproduces the *ladder* (half-integers in steps of one) and row 4 exactly, but predicts a 5/2 term in
(AA|BC) with coefficient ∝ ab × (A-potential averaged over the critical null segment with the Morse-Bott weight).
The corpus finds no 5/2 in the minimal operator (record MO). Either that segment average vanishes identically in the
angle (a cancellation worth a sentence in Paper 59) or the heuristic is wrong for complex-null critical lines (the
quadratic forms are degenerate along a null line, so the Morse-Bott hypothesis is not obviously met).

| Our result | Status | Reference | What it says | Credit to add | Shortcut / test |
|:--|:--|:--|:--|:--|:--|
| x* half-integer exponents read from the indicial equation | PARTIALLY KNOWN | Hannesdottir-Mizera 2204.02988 [V-text] Eq. (7.12)-(7.13) | φ^γ, γ = (LD−E−2d_N−1)/2 for simple pinches | cite as the momentum-space rule; say the contact point is a non-simple, position-space pinch | **Test:** compute the 5/2 coefficient (null-segment average); if it vanishes, the absence of 5/2 is explained |
| Position-space (coordinate) Landau analysis | PARTIALLY KNOWN | Erdoğan, PRD 89, 085016 (2014), arXiv:1312.0058 [V] | coordinate-space pinch surfaces and power counting (massless gauge theory) | optional context | none |
| Computing Landau loci algebraically | KNOWN (tool) | Fevola-Mizera-Telen, PRL 132, 101601 (2024), arXiv:2311.14669 [V-sec via arXiv listing]; companion 2311.16219 | principal Landau determinant; handles massless lines | none needed | **Shortcut:** PrincipalLandauDeterminants.jl could list the four-centre contact loci before an operator is computed |
| Legendre-coefficient decay set by the nearest singularity in cos θ | KNOWN | Lehmann, Nuovo Cim. 10, 579 (1958) [V-sec: JHEP ref lists and reviews only] | partial-wave series converges in an ellipse with foci ±1 set by the nearest singularity | add "Lehmann ellipse" beside Darboux / Bernstein | L^{−γ−1/2} rate already used |

## Q2. Minimal versus non-minimal operators, pieces, LCLM

| Our result | Status | Reference | What it says | Credit to add | Shortcut / test |
|:--|:--|:--|:--|:--|:--|
| AB's order-17, degree-48 operator is a left multiple of the minimal order-14, degree-111 operator; its non-integer "moving" points were spurious (record MO) | KNOWN phenomenon | Chen-Jaroschek-Kauers-Singer, arXiv:1301.0917 (2013) [V] | "Desingularization is the problem of finding a left multiple of a given Ore operator"; order-degree curves | cite for "a higher-order, lower-degree left multiple exists and carries extra singular points" | none; MO already did the right thing |
| Minimal-order operators for Feynman integrals via factorisation | KNOWN | Lairez-Vanhove, Lett. Math. Phys. 113 (2023), arXiv:2209.10962 [V] | extended Griffiths-Dwork reduction; "Using recent factorisation algorithms we give the minimal order differential operator" | cite next to the MO paragraph | their code is a cross-check for the order-14 claim |
| Minimal operators in general D | KNOWN | de la Cruz-Vanhove, arXiv:2401.09908 (2024) [V] | "an algorithm for determining the minimal order differential equations associated to a given Feynman integral" | optional | — |
| Operator factorisation | KNOWN | Müller-Stach-Weinzierl-Zayadeh, CMP (2014), DOI 10.1007/s00220-013-1838-3, arXiv:1212.4389 [V] | Fuchsian PF equations; "interesting factorisation properties ... when D is specialised to integer dimensions" | optional | — |
| Per-piece / per-constant operators, LCLM bounds the order | KNOWN as technique | textbook (LCLM of annihilators annihilates the sum); same logic as sub-topology inhomogeneities in the papers above [UNV for an explicit statement] | — | none; state it as standard | — |
| zeta_A-moving points logarithmic, integer exponents | PARTIALLY KNOWN | mass-dependent Landau loci (pseudo-thresholds) are standard in the papers above | — | none | — |

Our per-key search (FXK, a key operator carrying 7/2 at each contact) has no direct counterpart; the method is standard,
and the result is specific to this family.

## Q3. Exponential / E1 constants and depth-2 iterated integrals

| Our result | Status | Reference | What it says | Credit to add | Shortcut / test |
|:--|:--|:--|:--|:--|:--|
| Constant classes {exp, E1, ln} (and γ_E where present) | KNOWN as a class | Kontsevich-Zagier, *Periods* (2001) Sec. 4.3 [V-sec: Lagarias survey arXiv:1303.1856 only; IHES PDF did not extract] | exponential periods; γ is one | cite for the class name | — |
| — | KNOWN, rigorous frame | Commelin-Habegger-Huber, *Exponential periods and o-minimality*, arXiv:2007.08280 (to appear Ann. Inst. Fourier) [V-sec: listing] | compares naive, cohomological and Nori-motive exponential periods | optional | — |
| — | KNOWN | Fresán-Jossen, *Exponential motives* (manuscript, 2020) [UNV; seen only in citing sources] | Tannakian theory; γ appears in an exponential period matrix | cite if the class is named | — |
| Independence of the constants | PARTIALLY KNOWN, weak | Rivoal, Michigan Math. J. 61 (2012) 239-254 [V-text] Thm 2(ii) | for algebraic z, trdeg Q(γ + log z, e^z, G0(z)) ≥ 2, with G0(z) = ∫e^{−t}/(t+z)dt = e^z E1(z); so at least one of γ+log z, G0(z) is transcendental | add if the papers say anything about whether the constants are independent; nothing stronger is known | do not claim independence of our 21 constants; PSLQ "no relation" is the honest ceiling |
| Newer class | NEW (2026) | Snodgrass, *Periods of E-operators*, arXiv:2608.06005 [V] | "We define an enlargement of the class of exponential periods" | watch-list | — |
| Depth-2 iterated exponential integrals ∫E1(αt)e^{−βt}R(t)dt, E1×E1 | PARTIALLY KNOWN | Geller-Ng, J. Res. NBS 73B (1969) 191-210 [V-sec: NBS/NTRS catalogue; PDF too large to fetch]; Milgram, Math. Comp. 44 (1985) 443-458 [V-sec: AMS record] | tables of definite and indefinite integrals of E1 products; generalized integro-exponential function | cite Geller-Ng for the depth-2 reductions | **Test:** whether ∫e^{−pt}E1(at)E1(bt)dt brings in dilogarithms [UNV]. If it does, the four E1×E1 constants at M2r are the first place the "dilog-free" class could break |

No result found on relations among depth-2 exponential iterated integrals beyond integration by parts (searched under
exponential periods, E-functions and the Siegel-Shidlovskii setting).

## Q4. 3-D position-space graphs with Coulomb and Yukawa lines

| Our result | Status | Reference | What it says | Credit to add | Shortcut / test |
|:--|:--|:--|:--|:--|:--|
| 3-D massive graphs close in elementary functions plus few special functions | KNOWN, momentum space | Rajantie, Nucl. Phys. B480 (1996) 729, arXiv:hep-ph/9606216 [V] | three-loop 3-D vacuum diagrams: "For almost every integral an expression in terms of elementary and dilogarithm functions is obtained"; master and Mercedes (tetrahedron) only as 1-D integrals of elementary functions | yes, Paper 63 / inventory: the 3-D precedent for "elementary closes except at the tetrahedron" | compare our dilog-free (AA\|BC) class with his dilog entries |
| Banana in x-space = 1-D integral of products of propagators; closed form in odd D | KNOWN | Groote-Körner-Pivovarov, Nucl. Phys. B542 (1999) 515, arXiv:hep-ph/9806402 [V] | "Substantial simplifications occur for odd space-time dimensions"; explicit 3-D formulas | yes, next to `prop:tree_bound` (the tree's leaf-integration is their x-space product structure) | — |
| Legendre expansion in a bond angle of a position-space graph with fixed external points | NOT FOUND | — | — | — | — |
| Dimension-shift links | KNOWN tool | Tarasov, PRD 54 (1996) 6479, arXiv:hep-th/9606018 [V] | relates d and d−2 integrals by a differential operator | only if a 3-D ↔ 1-D/5-D link is used; it does not link 3-D to the 2-D K0 bananas (parity of d is kept) | — |

## Q5. Tree rung, Bessel-moment operator orders

| Our result | Status | Reference | What it says | Credit to add | Shortcut / test |
|:--|:--|:--|:--|:--|:--|
| Recurrences for ∫t^k K0^n | KNOWN | Borwein-Salvy, *A proof of a recursion for Bessel moments*, arXiv:0706.1409, Exp. Math. 2008 [V] | proves the BBBC 2007 recursion conjecture (orders not in the abstract) | yes, next to "annihilator of K0^m, order m+1" | — |
| Q-span dimension of Bessel moments from Vanhove operators | KNOWN | Zhou, Ramanujan J. 58 (2022) 723-746, arXiv:1911.04141 [V] | bound ⌊(a+b−1)/2⌋ on the Q-dimension, improving Borwein-Salvy's ⌊(a+b+1)/2⌋; exceptional relation at a=2, b=6 | Paper 63 does not cite Zhou (Paper 59 does): add | **Test:** the K0^n moments in the tree period obey these Q-relations; a free check of `prop:tree_rung` outputs |
| Motives of symmetric powers (attainment by irreducibility) | KNOWN frame | Fresán-Sabbah-Yu, *Hodge theory of Kloosterman connections*, Duke 171 (2022) 1649, arXiv:1810.06454 [V] | motives for symmetric power moments; every Hodge number 0 or 1 | possibly; irreducibility for all n was not seen in the abstract [UNV] | if it holds for all n, attainment n ≤ 8 could be lifted to all n. Read the paper before saying so |
| Tree to 1-D radial recursion via the shell theorem | PARTIALLY KNOWN | GKP 1999 (x-space product) + Newton's shell theorem | — | as above | — |

---

## Top 5 things the PI should know

1. **The exponent at x* is a position-space, non-simple pinch.** The standard momentum-space formula
   (Hannesdottir-Mizera Eq. 7.12) does not apply directly. A Morse-Bott / Fourier power count (ours, heuristic)
   reproduces row 4's 3/2 and the 7/2, 9/2 ladder. It also predicts a **5/2 term in (AA|BC) that the corpus does not
   see**. Either a computable null-segment average vanishes, or the heuristic fails. That check is cheap and would
   either explain the gap or flag it.
2. **MO's finding is textbook computer algebra.** AB's operator was a desingularised left multiple whose extra
   singular points are artefacts. Cite Chen-Jaroschek-Kauers-Singer and Lairez-Vanhove; Lairez-Vanhove's software is an
   independent check of order 14.
3. **No independence theorem covers our constants.** The best known (Rivoal 2012) is a transcendence-degree ≥ 2
   statement mixing γ+log z, e^z and e^z E1(z). Claims must stay at "no relation found" (PSLQ).
4. **Rajantie 1996 and Groote-Körner-Pivovarov 1999 are the 3-D precedents.** In 3-D, massive graphs close in
   elementary functions and dilogarithms, and only the tetrahedron-type masters stay as 1-D integrals. Both are uncited.
   Ask whether E1×E1 products (four at M2r) bring dilogarithms back once integrated (Geller-Ng tables).
5. **Paper 63 should cite Zhou 2022 and Borwein-Salvy.** Zhou's Q-linear relations among K0 moments give a free
   consistency test of the tree-rung outputs. Fresán-Sabbah-Yu's Duke paper may give irreducibility for all n, which
   would lift the n ≤ 8 limit; this is unverified and needs reading first.

Not found anywhere: a Legendre/bond-angle ODE for multicentre position-space graphs; the 21-constant fixed set; the
cancellation of the 3/2 branch. These still read as new, scoped to "not found in this search".
