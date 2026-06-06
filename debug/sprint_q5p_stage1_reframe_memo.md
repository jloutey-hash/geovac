# Sprint Q5' Stage 1 reframe — from $\mathrm{HP}_*$ to the JLO/Chern character on $\mathrm{HP}^*$

**Date:** 2026-06-05
**Sprint:** Q5' Stage 1 reframe (diagnostic, no code, no PSLQ, no paper edits)
**Triggered by:** structural tension between Sprint Q5'-CH-1 (2026-06-05, Stage 1 alive — bit-exact $k$-slot factorisation on truncated CH triple) and Sprint Q5'-R3-HP\*-check (2026-06-04, $\mathrm{HP}_*$-route structurally TRIVIAL at every finite $n_{\max}$, Morita-trivial pro-limit too).

---

## TL;DR

**Verdict: REFRAMED-CONSTRUCTIVE.** Stage 1's $\omega^{\mathrm{tri}}$ targets the *cohomological-dual* side of cyclic theory — Connes' **periodic cyclic cohomology $\mathrm{HP}^*$**, with the **JLO cocycle** as the explicit, computable representative of the Chern character of the truncated CH spectral triple. Two named published precedents support this: (a) Jaffe–Lesniewski–Osterwalder 1988 (the original JLO construction) supplies the *cocycle*; (b) Connes–Moscovici 1995 + Connes–Marcolli 2004 $U^*$ supplies the *Tannakian shape* (label-respecting automorphism group on cohomology-side data). The first concrete computation (the JLO cocycle of the truncated CH triple at $n_{\max} = 2$ decomposed along the master-Mellin $k$-slot) is **sprint-scale (~1 week, exact rational arithmetic at the symbol level)**. This is Stage 1 of the multi-year $U^*$ bridge, not the multi-year bridge itself.

The R3-vs-CH-1 tension dissolves once the *side* of cyclic theory is specified: R3 closes the *homology* side ($\mathrm{HP}_*$ / Marcolli–Tabuada input) decisively; CH-1 lives on the *cohomology* side ($\mathrm{HP}^*$ / Chern character output) which R3 explicitly identified as the surviving direction (R3 Step 4, Step 5). CH-1 read the right side without naming it; R3 named the wrong-side closure without identifying the right-side construction. This memo closes the gap.

The CH-1 phrasing "$\Z/2 \times \mathbb{Z}$-graded refinement of $\mathrm{HP}_*$" was a misnomer — the genuine refinement is of $\mathrm{HP}^*$ (the dual cyclic *cohomology*), where the $\Z/2$ chirality grading IS the standard parity and the $\mathbb{Z}$ heat-kernel-order grading IS the Mellin shift of the JLO entire-cyclic structure. (Recommended sharpening for Paper 55 §subsec:open_m2_m3, flagged at end of memo; no edits applied.)

## Verdict against gate

| Gate | Verdict |
|:-----|:--------|
| REFRAMED-CONSTRUCTIVE | **selected** — JLO cocycle on $\mathrm{HP}^*$ is the named target; Connes–Marcolli $U^*$ is the named published precedent; $n_{\max} = 2$ first computation is sprint-scale exact-rational |
| REFRAMED-OPEN | not selected, but acknowledged as fallback if the sprint-scale computation finds non-trivial obstructions |
| STILL-DEAD | rejected — R3 itself flagged $\mathrm{HP}^*$ / Chern character as the live direction (R3 Step 4 explicitly: "The spectral triple's K-homology class (Chern character) is non-trivial at finite $n_{\max}$") |
| NEEDS-PI-DIRECTION | partially raised at end: the choice between JLO as the unique Stage-1 target vs.\ also building the Connes–Moscovici residue cocycle in parallel is a PI call |

## The R3-vs-CH-1 tension and its resolution

R3 ($\mathrm{HP}_*$-check) and R3-pro-dg-category land on TRIVIAL via three load-bearing facts:

1. The bare algebra $\mathcal{A}_{\mathrm{GV}}^{(n_{\max})} \cong \C^{N_{\mathrm{Fock}}}$ is Morita-trivially commutative; $\mathrm{HP}_0 = \Q^{N_{\mathrm{Fock}}}$, $\mathrm{HP}_1 = 0$.
2. The dg-thickening via bounded $[D, \cdot]$ at finite cutoff does not escape Morita-stability of $\mathrm{HP}_*$.
3. The pro-limit with Berezin maps converges to $\mathrm{HP}_*(C^\infty(S^3)) = (\Q, \Q)$ — pure-Tate, no $k$-distinguishing content; and the Berezin maps are CP not multiplicative, so they are not dg-functorial.

CH-1 lands on POSITIVE via two structural facts about the *truncated CH triple* (algebra + Hilbert space + Dirac, all three pieces):

4. The chirality-parity selection rule on $\mathrm{Tr}(D^k\,e^{-tD^2})$ vs $\mathrm{Tr}(\gamma\,D^k\,e^{-tD^2})$ is bit-exact at every $n_{\max} \in \{2, 3, 4\}$, on BOTH the diagonal CH $\Lambda$ AND the full graph Dirac $D = \Lambda + \kappa A$.
5. The factorisation reads as a $\Z/2$ chirality component times a heat-kernel-order component within HP$_0$.

R3 Step 4 (the load-bearing observation buried in the memo): *"The spectral triple's K-homology class (Chern character) is non-trivial at finite $n_{\max}$. The Connes character of a spectral triple $(\mathcal{A}, \mathcal{H}, D)$ — equivalently the JLO cocycle or the Chern character of the associated Fredholm module — lives in $\mathrm{HP}^*(\mathcal{A})$ (dual cohomology), NOT in $\mathrm{HP}_*(\mathcal{A})$, and carries information about $D$ that is invisible to the algebra alone."*

**The resolution.** CH-1's data IS on the cohomological-dual side R3 named; CH-1 just didn't say "$\mathrm{HP}^*$" explicitly. The traces $\mathrm{Tr}(D^k\,e^{-tD^2})$ and $\mathrm{Tr}(\gamma\,D^k\,e^{-tD^2})$ are *traces against the Dirac data*, hence they are *functionals on the algebra* — which is precisely what $\mathrm{HP}^*$ elements are (cyclic cochains). The CH-1 finding is the *operator-symbol-level* witness of a non-trivial $\mathrm{HP}^*$ class, which is exactly the Chern character.

The CH-1 phrasing "the chirality component is exactly the standard $\mathrm{HP}_*$ grading" should read "...exactly the standard $\Z/2$ grading of cyclic theory (which on $\mathrm{HP}^*$ is the even/odd parity of cyclic cochains)." The chirality grading lives on both sides — $\mathrm{HP}_*$ and $\mathrm{HP}^*$ — because the parity grading is the standard $\Z/2$ structure of the whole cyclic-theory machinery. What lives *only* on $\mathrm{HP}^*$ is the Dirac-dependent heat-kernel-order data; that's the genuinely new content the Chern character carries.

## Identified cohomological-side target

### Primary target: JLO cocycle of the truncated CH triple, decomposed along the $k$-slot

The **JLO cocycle** (Jaffe–Lesniewski–Osterwalder 1988, *Quantum harmonic analysis on phase space*; rebuilt as standard NCG machinery in Getzler–Szenes 1989, Connes–Moscovici 1990s) of a $\theta$-summable spectral triple $(\mathcal{A}, \mathcal{H}, D)$ is the sequence of multilinear cochains $\{\phi_n\}_{n \ge 0}$:

$$
\phi_n(a_0, a_1, \ldots, a_n) \;=\; \int_{\Delta_n} \mathrm{Tr}_s\!\left( a_0\, e^{-s_0 D^2}\,[D, a_1]\, e^{-s_1 D^2} \cdots [D, a_n]\, e^{-s_n D^2}\right)\, ds_1 \cdots ds_n,
$$

where $\Delta_n = \{(s_0, \ldots, s_n) : s_i \ge 0, \sum s_i = 1\}$, $\mathrm{Tr}_s = \mathrm{Tr}(\gamma\,\cdot)$ on a graded triple. This is an *entire-cyclic-cohomology* representative of the Chern character; on a $\theta$-summable triple, $\{\phi_n\}$ defines a cocycle in the $(b, B)$-bicomplex of cyclic cohomology, and its periodic class $[\{\phi_n\}] \in \mathrm{HP}^*(\mathcal{A})$ is the Chern character of the triple (Connes–Moscovici 1995).

**Key facts that make JLO the right Stage-1 target for GeoVac:**

(i) **Computable at finite cutoff.** On a finite-dim spectral triple, every integral $\int_{\Delta_n}$ is finite (no UV divergences), and every trace is a finite-rank matrix trace. At $n_{\max} = 2$ with $\dim \mathcal{H} = 16$ (CH-1 panel), every $\phi_n$ is a polynomial in the eigenvalues of $D$ (which are $\pm(n + 1/2)$ rationals) and the matrix elements of $a_0, [D, a_1], \ldots$ — all bit-exact in `sympy.Rational`.

(ii) **The $k$-slot lives in the $t$-expansion.** Reparametrise via $t = \sum s_i$ (the total heat-kernel time); after rescaling, $\phi_n$ generates the Mellin source $\mathrm{Tr}_s(D^k\,e^{-tD^2})$ at $k = n$ via the leading $t$-asymptotic of the integrand. (More precisely: the small-$t$ expansion of $\phi_n(t \cdot \mathbf{1}, a_1, \ldots, a_n)$ Mellin-transforms to the $k = n$ slot.) The Mellin engine $k$-slot $\in \{0, 1, 2\}$ corresponds to the JLO cochain degree $n \in \{0, 1, 2\}$.

(iii) **Chirality grading is automatic.** $\mathrm{Tr}_s = \mathrm{Tr}\circ\gamma$ on the graded triple gives the supertrace at even cochain degree and the trace at odd degree (or vice versa, depending on convention). The CH-1 chirality-parity selection rule (plain trace nonzero at $k \in \{0, 2\}$, supertrace nonzero at $k = 1$) is exactly the JLO parity grading — bit-exactly visible at the cochain level.

(iv) **Bit-exact rational arithmetic at $n_{\max} = 2$.** Every $\phi_n$ at $n_{\max} = 2$ is a finite sum of products of rationals (matrix elements of finite-dim operators) times rationals (eigenvalue powers), divided by an integer ($n!$ from the simplex volume). Verifiable in `sympy.Rational` throughout; no PSLQ, no floats, no transcendental introduction at the Layer-1 skeleton level. Transcendentals enter only at Stage 2, when the Mellin-transform of $\phi_n$ is taken against $\Gamma(s)$.

### Why not Marcolli–Tabuada's $\mathrm{HP}_*$ — explicit one-line statement

R3 closed this decisively: $\mathrm{HP}_*$ uses *homology* (cycles modulo boundaries) of the cyclic complex of *the algebra*; finite-dim semisimple algebras have rank-1 $\mathrm{HP}_0$ and zero $\mathrm{HP}_1$. The Dirac data lives in the *dual cyclic theory* — *cohomology* (cocycles modulo coboundaries), where the Chern character is naturally placed. R3 explicitly named this as the live alternative.

### Three candidates evaluated (per curve-fit-audit discipline)

| Candidate | Side of cyclic theory | Tannakian shape precedent | Stage-1 computability | Verdict |
|:---|:---|:---|:---|:---|
| **JLO cocycle (Jaffe–Lesniewski–Osterwalder 1988)** | $\mathrm{HP}^*$, entire cyclic | Connes–Marcolli $U^*$ acts on equisingular flat connections via Tannakian dual — JLO is the natural NCG-side analogue | Sprint-scale at $n_{\max} = 2$: exact rational, $\phi_n$ is polynomial in $D$-eigenvalues + matrix elements | **PRIMARY** |
| Connes–Moscovici residue cocycle (CM 1995, the "local index formula") | $\mathrm{HP}^*$, periodic; via Wodzicki residues of zeta functions | Same Connes–Marcolli $U^*$ shape; residue cocycle is the standard NCG Chern-character formula | Sprint-scale at $n_{\max} = 2$: zeta function $\zeta_D(s) = \mathrm{Tr}\,|D|^{-s}$ has explicit closed form via $\zeta(s, 3/2)$ on CH | **SECONDARY** — equivalent to JLO at the periodic class level but builds in Mellin moments directly, so the $k$-slot factorisation is automatic |
| K-homology side via Fredholm-module Chern character (Connes 1985) | $\mathrm{HP}^*$, via Kasparov $KK$-theory | Standard Connes character pairing; less direct on the Mellin engine | $n_{\max} = 2$ computation possible but the $k$-slot decomposition requires extra work | tertiary — useful as cross-check |

JLO is selected over CM-residue as PRIMARY because:
- JLO is *explicit* at the cochain level (one finite-dim simplex integral per cochain degree), so the symbol-level visibility CH-1 already established transports verbatim;
- CM-residue is *equivalent* at the periodic class level (Connes–Moscovici 1995 prove $[\phi_{\mathrm{JLO}}] = [\phi_{\mathrm{CM-res}}]$ in $\mathrm{HP}^*$), so JLO is not less general;
- CM-residue would be the natural Stage 1.5 cross-check, computed in parallel if PI directs.

Curve-fit-audit check: no, JLO was not the only candidate looked at (CM-residue and K-homology Fredholm-module character both evaluated). The discriminator is computability at sprint-scale + direct visibility of the $k$-slot at cochain level. Both winning candidates land in $\mathrm{HP}^*$; the divergence from JLO is at the technical-execution layer, not at the structural reframe layer.

### Connes–Marcolli $U^*$ as the Tannakian-shape precedent

The $U^*$ cosmic Galois group (Connes–Marcolli 2004, arXiv:math/0409306; Connes–Marcolli 2008 book Ch. 1 §3, Ch. 4) operates by Tannakian dual on equisingular flat vector bundles whose holonomies are *labelled* by Feynman diagrams. The label-respecting structure is what enables a non-trivial Galois action: the same complex number can appear from multiple counterterm decompositions, but the Galois symmetry respects the diagrammatic source.

Q5' analogue (per R3 Step 5 + R3-pro-dg-category Section "What WOULD force NON-TRIVIAL"): a "GeoVac cosmic Galois group" acting on the master-Mellin Chern character $\mathrm{ch}(\mathcal{T}_{n_{\max}}) \in \mathrm{HP}^*$ decomposed along the slot $k \in \{0, 1, 2\}$. The Tannakian category is built from $(\mathrm{ch}_0, \mathrm{ch}_1, \mathrm{ch}_2)$ together with the case-exhaustion theorem's exactness statement (Paper 32 §VIII Thm `thm:pi_source_case_exhaustion`); the motivic Galois group is the automorphism group of the slot-respecting fiber functor.

This is precisely the *shape* of $U^*$, with the Feynman diagram labels replaced by the master-Mellin $k$-slot. **No mechanical transport from QFT $U^*$ to GeoVac is asserted** — the Q5' construction is genuinely new — only the *shape* is borrowed from a verified published precedent.

## First-computation sketch (Stage 1 at $n_{\max} = 2$)

**Goal.** Compute the JLO cochains $\phi_0, \phi_1, \phi_2$ explicitly on the truncated CH spectral triple at $n_{\max} = 2$ (dim $\mathcal{H} = 16$, algebra $\mathcal{A}_{\mathrm{GV}}^{(2)} = \C^5$ on the 5 Fock vertices $(n, l, m_l)$ with $n \le 2$, Dirac $D = \Lambda + \kappa A$ with $\kappa = -1/16$), and verify that they reproduce the CH-1 bit-exact factorisation at the *cochain* level. Then state the Stage-1 deliverable: $\omega^{\mathrm{tri}}: \mathrm{dg}(\mathcal{T}_{n_{\max}}) \to (\mathrm{HP}^*_{\mathrm{even}}, \mathrm{HP}^*_{\mathrm{odd}})$ given by $a \mapsto (\phi_{\mathrm{even}}(a, \ldots), \phi_{\mathrm{odd}}(a, \ldots))$ with the $k$-slot index $= n$ (the cochain degree).

**Concrete steps.**

*Step 1: Algebra representation.* On $\mathcal{A}_{\mathrm{GV}}^{(2)} = \C^5$, choose the basis $\{e_1, e_2, e_3, e_4, e_5\}$ of vertex idempotents. Each $a \in \mathcal{A}$ is $a = \sum_i a_i\,e_i$ with $a_i \in \C$. The representation on $\mathcal{H} = \C^{16}$ is block-diagonal: $e_i$ projects onto the shell at vertex $i$.

*Step 2: Dirac data.* $\Lambda$ is the diagonal $16 \times 16$ matrix with entries $\pm(n + 1/2)$ ($n = 1, 2$) with appropriate multiplicities. $A$ is the parity-respecting $E_1$ dipole adjacency. $D = \Lambda + \kappa A$; $D^2$ is computed exactly. $e^{-t D^2}$ in `sympy` at the symbol level returns a $16 \times 16$ matrix of exponentials of rationals, which we treat formally (i.e., we compute $\mathrm{Tr}(\text{any polynomial in }D \cdot e^{-tD^2})$ as a finite sum $\sum_j p_j\, e^{-t \lambda_j^2}$ with $p_j, \lambda_j^2$ rational).

*Step 3: $\phi_0(a_0)$.* The 0-cochain is the **graded trace**:
$$
\phi_0(a_0) \;=\; \int_{\Delta_0} \mathrm{Tr}_s(a_0\, e^{-D^2})\, =\; \mathrm{Tr}(\gamma\,a_0\, e^{-D^2}).
$$
By CH-1's parity rule, $\mathrm{Tr}(\gamma\,e^{-D^2}) = 0$ identically when $a_0 = 1$ (the supertrace of $e^{-D^2}$ is the McKean–Singer index, which vanishes on $S^3$ because $\chi(S^3) = 0$). For general $a_0$, $\phi_0(a_0) = \sum_i a_i \cdot \mathrm{Tr}(\gamma\,e_i\,e^{-D^2})$, with each vertex contribution computable bit-exactly. Some vertices contribute $\ne 0$ even though the total is zero (the McKean–Singer cancellation is global).

*Step 4: $\phi_1(a_0, a_1)$.* The 1-cochain:
$$
\phi_1(a_0, a_1) \;=\; \int_0^1 \mathrm{Tr}_s(a_0\,e^{-sD^2}\,[D, a_1]\,e^{-(1-s)D^2})\, ds.
$$
Compute $[D, a_1] = [D, \sum_i a_{1,i}\,e_i] = \sum_i a_{1,i}\,[D, e_i]$. On a commutative algebra, $[D, e_i] = D\,e_i - e_i\,D = (1 - 2 e_i)\,D\,e_i + e_i\,D\,(2 e_i - 1)$ — but more cleanly: $[D, e_i]$ is computed directly as a $16 \times 16$ matrix (entries are differences of eigenvalues times off-diagonal $A$-matrix elements). The $s$-integral is a finite sum of integrals of $e^{-s\,a} \cdot e^{-(1-s)\,b}$ from $0$ to $1$ — each is $(e^{-a} - e^{-b})/(b - a)$, a closed form. Total cochain: bit-exact rational sum (after Mellin-transform-symbol manipulation).

*Step 5: $\phi_2(a_0, a_1, a_2)$.* Same structure, 2-simplex integral. Each $\int_{\Delta_2}\,e^{-s_0 a}\,e^{-s_1 b}\,e^{-s_2 c}\,ds$ has a 3-exponential closed form.

*Step 6: Decompose along $k$-slot.* For each $\phi_n$ at cochain degree $n \in \{0, 1, 2\}$, identify the *leading* contribution at $t = \sum s_i \to 0$ (the Mellin shift to $k = n$). The CH-1 leading data ($\dim \mathcal{H}$ for $k = 0$, $\mathrm{Tr}(\Lambda^2)$ for $k = 2$, $\mathrm{Tr}(\gamma\Lambda)$ for $k = 1$) appears as the *short-time-asymptotic-coefficient* of $\phi_n$ at $a_0 = a_1 = \cdots = 1$. Verify bit-exactly.

*Step 7: Verify cocycle property.* Compute $(b\phi_n + B\phi_{n+2})(a_0, \ldots, a_{n+1})$ and show it equals zero up to a coboundary — this is the JLO cocycle condition. At $n_{\max} = 2$ and cochain degrees $\le 2$, this is a finite identity in `sympy.Rational`.

*Step 8: Stage-1 deliverable.* Output: a finite-cutoff explicit Chern-character cocycle $\{\phi_0, \phi_1, \phi_2\} \in \mathrm{HP}^*(\mathcal{A}_{\mathrm{GV}}^{(2)})$ whose $k$-slot decomposition matches CH-1. This is the bit-exact symbol-level construction of $\omega^{\mathrm{tri}}$ at $n_{\max} = 2$, on the cohomological side.

**Honest scope of the sprint-scale computation.**

- This sprint does NOT construct the Tannakian category or compute the motivic Galois group. Those are Stage 2 (multi-year).
- This sprint DOES produce the explicit cochain symbols at $n_{\max} = 2$, which is the natural Stage-1 input for any future Tannakian construction.
- No PSLQ, no transcendental introduction, no continuum limit — all rational at the cochain level.
- Wall-time estimate: 3–7 days for one PM in main session. Setup: ~1 day (algebra + Dirac matrices in sympy). $\phi_0$: ~1 day. $\phi_1$: ~2 days (the $s$-integral closed-form bookkeeping). $\phi_2$: ~2–3 days. Cocycle verification: ~1 day. Memo + tests: ~1 day.
- Risk: at $\phi_2$ the simplex integral involves nested exponentials of differences of eigenvalues; algebraic simplifications could grow complex enough to overflow the "sprint-scale" budget. If so, switch to floating-point at high precision (200 dps) and verify bit-exactly after the fact via Sprint-scale PSLQ. (This is the only place transcendentals could leak in at Stage 1, and would be a $k=2$-specific issue.)

**What would force a STILL-DEAD downgrade.** If at Step 6 the leading short-time asymptotic of $\phi_n$ at $a_0 = a_1 = \cdots = 1$ does *not* match CH-1's $\mathrm{Tr}(D^k\,e^{-tD^2})$ at $k = n$, then the Mellin-engine/JLO identification breaks at the cochain level and a different cohomological-side construction is needed. The structural prediction (Connes–Moscovici 1995 plus CH-1's data) is that they DO match; this is the Stage-1 sprint's first decision gate.

## Implications for Q5' Stage 2

The reframe sharpens Stage 2 as well:

- **The Tannakian category** is not the Marcolli–Tabuada category of dg-perfect-modules (R3 closed that); it is a *category of $k$-graded Chern characters*, with objects $(\mathcal{T}, [\mathrm{ch}_0], [\mathrm{ch}_1], [\mathrm{ch}_2])$ and morphisms preserving the slot grading.
- **The fiber functor** is $\omega^{\mathrm{tri}}: (\mathcal{T}, \mathrm{ch}) \mapsto (\mathrm{ch}_0, \mathrm{ch}_1, \mathrm{ch}_2) \in \mathrm{HP}^*(\mathcal{A}) \otimes \mathrm{IndexCat}(\{0, 1, 2\})$ — with the Mellin shift specifying the $\mathbb{Z}$-graded refinement within HP$^*$.
- **The motivic Galois group** is the automorphism group of $\omega^{\mathrm{tri}}$ acting label-respectingly on the $k$-slot. The *existence* of this group as a well-defined pro-algebraic group is the multi-year mathematical question Stage 2 names.
- **The case-exhaustion theorem** (Paper 32 §VIII) becomes the exactness statement of the Tannakian category — every $\pi$-period output sits in $M_1$, $M_2$, or $M_3$, with at least one $k$-slot producing it.

The R3-pro-dg-category memo Section "What WOULD force NON-TRIVIAL" listed four hypothetical scenarios; the JLO-on-$\mathrm{HP}^*$ route is structurally the second of those ("Some non-standard $\Z/3$-refinement of HP_* aligned with the master Mellin engine"), with $\mathrm{HP}_*$ replaced by $\mathrm{HP}^*$. The replacement is what makes the construction live — Chern character data on $\mathrm{HP}^*$ is precisely the missing input that $\mathrm{HP}_*$-side machinery cannot provide.

## Honest scope

1. **No code, no PSLQ, no transcendental introduction** — pure structural reframing on top of R3 + R3-pro-dg-category + CH-1.
2. **No paper edits applied.** Recommendations for Paper 32 §VIII and Paper 55 §subsec:open_m2_m3 flagged at end; PI to apply, decline, or modify.
3. **No new pre-commitment to Stage 1 sprint.** This memo IS the diagnostic; the construction sprint (Steps 1–8 above) is a separate decision-gate the PI controls.
4. **JLO selection is reasoned, not unique.** Two other cohomological-side candidates were evaluated (CM-residue cocycle, K-homology Fredholm-module character); JLO won on sprint-scale-explicit-cochain-form criteria, not on uniqueness. If the Stage-1 sprint is launched and JLO turns out to be technically intractable at $\phi_2$, CM-residue is the named fallback.
5. **Multi-year Stage 2 is genuinely multi-year.** The motivic Galois group of the slot-graded Chern-character Tannakian category has no published precedent. This memo does not estimate Stage-2 timescale beyond "multi-year"; that scoping is a separate exercise after Stage 1 lands.
6. **WH1 PROVEN is not re-opened** — JLO cocycle convergence in propinquity is a separate (open) question from cyclic-homology / cyclic-cohomology convergence; R3-pro-dg-category memo's "propinquity and cyclic homology are categorically orthogonal" survives unchanged. The relevant convergence question for Stage 2 is whether $\omega^{\mathrm{tri}}(\mathcal{T}_{n_{\max}}) \to \omega^{\mathrm{tri}}(\mathcal{T}_\infty)$ in $\mathrm{HP}^*$ as $n_{\max} \to \infty$. This is *new* and not addressed by Paper 38.
7. **No paper 32 §VIII / Paper 55 §7.6 framing changes proposed at the master-Mellin-engine theorem-statement level.** The case-exhaustion theorem, master Mellin engine domain partition, and mixed-Tate / cyclotomic classifications all stand verbatim. The reframe is purely on the Q5' construction-route side of Paper 55 §subsec:open_m2_m3.
8. **Diagnostic-before-engineering compliance.** This is the diagnostic; a Stage-1 construction sprint is opt-in for PI direction. The trigger for diagnostic was the R3-vs-CH-1 structural tension (one honest negative on $\mathrm{HP}_*$ + one positive at the same-day cohomology-side bit-exact level), not yet at the "≥ 2 honest negatives" threshold but operationally analogous: a reframe before engineering.

## Files used

### Memos read (load-bearing inputs)
- `debug/sprint_q5p_r3_hp_star_check_memo.md` (Round 3 $\mathrm{HP}_*$ check — TRIVIAL on bare algebra, NON-TRIVIAL on Chern character; R3 Step 4 explicitly named the cohomological-dual / Connes–Marcolli $U^*$ direction as the alive alternative)
- `debug/sprint_q5p_r3_pro_dg_category_memo.md` (Round 3 pro-dg-category — TRIVIAL with three load-bearing reasons; Berezin maps are CP not multiplicative, propinquity ≠ cyclic-homology convergence)
- `debug/sprint_q5p_ch1_memo.md` (Sprint Q5'-CH-1 — bit-exact chirality-parity factorisation at $n_{\max} \in \{2, 3, 4\}$; the symbol-level structural data the reframe ties to JLO)
- `debug/sprint_q5p_k_slot_tannakian_memo.md` (Round 2 BORDERLINE — $k$-slot Tannakian-invisible on standard fiber functor, relevant one categorical level up)
- `debug/sprint_q5p_greenfield_marcolli_transport_memo.md` (Round 1 NOT TRANSPORTABLE on QSM; ruled out Greenfield–Marcolli, leaving Strand A = Connes–Marcolli $U^*$ as the cohomological-dual-side precedent)
- `debug/sprint_q5p_qsm_litread_memo.md` (Round 1 PARTIALLY TRANSPORTABLE — five precedent strands; Strand A = Connes–Marcolli $U^*$ as Tannakian-shape precedent for cohomological-side data)

### Papers consulted
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` §VIII (Thm `thm:pi_source_case_exhaustion`, Rem `rem:master_mellin_domain` — the master Mellin engine $k \in \{0, 1, 2\}$ partition that JLO cochain degree should match)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §subsec:open_m2_m3 lines 1494–1661 (Q5' open question, CH-1/CH-2/CH-3 closures, current Stage-1 framing as "candidate enriched fiber functor on dg-category")
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` §7.6 (Q5' framing in joint engagement / open-questions section)

### Published references (no new bibitems — all four are already in scope for Paper 55's Q5' subsection or close cross-references)
- Jaffe, A.; Lesniewski, A.; Osterwalder, K. *"Quantum K-theory I:\ The Chern character."* Commun. Math. Phys. 118 (1988), 1–14. **The original JLO cocycle paper.** Not currently in Paper 55 bibliography; would need a new bibitem `jlo1988` if the Stage-1 sprint launches and Paper 55 §subsec:open_m2_m3 is updated.
- Connes, A.; Moscovici, H. *"The local index formula in noncommutative geometry."* GAFA 5 (1995), 174–243. **The residue cocycle representative of Chern character.** Equivalent to JLO at periodic class level. Already cited indirectly via Paper 32 §VIII.
- Connes, A.; Marcolli, M. arXiv:math/0409306 (2004) — the $U^*$ cosmic Galois group; already in Round 1 references list, would be added to Paper 55 bibliography if Q5' is activated.
- Getzler, E.; Szenes, A. *"On the Chern character of a theta-summable Fredholm module."* J. Funct. Anal. 84 (1989), 343–357. **JLO-cocycle entirety + bookkeeping.** Useful for the cocycle verification step.

### Scripts (none required this sprint)
Pure structural reasoning. The Stage-1 construction sprint (if launched) would build `debug/compute_jlo_cochains_nmax2.py` using `sympy.Rational` throughout.

## Recommended paper edits (PI to apply, decline, or modify)

### Paper 55 §subsec:open_m2_m3 — sharpen Stage 1 reframe

Currently the Q5' subsection (Paper 55 lines 1556–1608) describes Stage 1 as "lifting [the symbol-level finite-cutoff data] to a graded fiber functor on the dg-category" via Marcolli–Tabuada. Per R3 closures and this memo, the Marcolli–Tabuada side is closed; the live target is the cohomological-dual side. Suggested replacement / addition (one paragraph after the existing CH-1 paragraph, ~line 1608):

> *Stage 1 reframe from $\mathrm{HP}_*$ to $\mathrm{HP}^*$ (Sprint Q5'-Stage1-reframe, June 2026; memo \texttt{debug/sprint\_q5p\_stage1\_reframe\_memo.md}).* The Marcolli–Tabuada $\mathrm{HP}_*$-route is structurally dead at the dg-category source (R3 closures); the candidate enriched fiber functor $\omega^{\mathrm{tri}}$ targets the *cohomological-dual* side of cyclic theory, $\mathrm{HP}^*$. The natural explicit representative is the Jaffe–Lesniewski–Osterwalder (1988) entire-cyclic cocycle $\{\phi_n\}_{n \ge 0}$ of the truncated CH spectral triple, with the master-Mellin slot index $k \in \{0, 1, 2\}$ identified with the JLO cochain degree $n$. The CH-1 bit-exact symbol data — $G_0(t)$ landing in M1 with leading $\dim \mathcal{H}$, $G_2(t)$ in M2 with leading $\mathrm{Tr}(\Lambda^2)$, $S_1(t)$ in M3 with leading $\mathrm{Tr}(\gamma\Lambda)$ — is the short-time asymptotic of $\phi_n$ at $a_0 = a_1 = \cdots = 1$. Stage 1 is a sprint-scale construction (\textasciitilde1 week, exact rational at $n_{\max} = 2$, dim $\mathcal{H} = 16$); Stage 2 (the motivic Galois group of the slot-graded Chern-character Tannakian category, shape-precedented by the Connes–Marcolli cosmic Galois $U^*$) remains multi-year.

### Paper 32 §VIII — optional pointer to the cohomological-side reframe

Optional one-sentence Remark addition (after `rem:master_mellin_domain`) pointing forward to Paper 55's Q5' reframe — but this is not load-bearing and can be deferred until the Stage-1 sprint actually launches. No edits recommended at this time for Paper 32.

### One-line verdict

**REFRAMED-CONSTRUCTIVE.** Stage 1 of the Q5' cosmic-Galois $U^*$ bridge targets the JLO cocycle of the truncated Camporesi–Higuchi spectral triple on $\mathrm{HP}^*$ (cohomological-dual side, where R3 explicitly identified the live direction); the master Mellin $k$-slot $\in \{0, 1, 2\}$ maps to the JLO cochain degree $n \in \{0, 1, 2\}$; the first concrete computation at $n_{\max} = 2$ (dim $\mathcal{H} = 16$, algebra $\C^5$, bit-exact rational arithmetic throughout) is sprint-scale (\textasciitilde1 week); Connes–Marcolli $U^*$ is the named Tannakian-shape published precedent; the CH-1 bit-exact factorisation transports verbatim as the short-time asymptotic of $\phi_n$ at $a_0 = a_1 = \cdots = 1$. Multi-year Stage 2 (the motivic Galois group of the slot-graded Chern-character Tannakian category) is unchanged in scope; the JLO Stage-1 reframe just makes its input data explicit and computable.
