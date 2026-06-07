# Tapušković 2023 methodology — read for adoption potential against GeoVac

**Author / harness:** Methodology survey sub-agent, 2026-06-06 (evening, post-HB-PSLQ negative).
**Status:** Read-only survey. No paper, code, or memory modifications.
**Primary source:** Tapušković 2023, arXiv:2303.17534v1, "The cosmic Galois group, the sunrise Feynman integral, and the relative completion of $\Gamma_1(6)$" (51 pp, published Comm. Number Theory Phys. 18 (2024) no. 2). Verified by direct PDF read of arXiv preprint.

**Sprint context.** Today's Hain–Brown PSLQ Test A (`debug/sprint_hb_pslq_test_memo.md`) returned NEGATIVE at the period level for $\mathrm{Sym}^2$-tagged GeoVac periods at $\tau = 2i$: every PSLQ-identified period lands in $\mathrm{MZV}(\mathrm{MTM}(\mathbb{Z}))$ extended by level-4 Dirichlet content $\{G, \beta(4)\}$, with no engagement of $E_4, E_6, \Delta$ or Eichler periods of $\Delta$. The categorical shape match (relative completion $1 \to \mathcal{U} \to \mathcal{G}^{\mathrm{rel}} \to SL_2 \to 1$ with $\mathrm{Sym}^k V_{\mathrm{fund}}$ tensor category and pro-unipotent radical) remains structurally real. Tapušković is the *only* published physics-side input into the relative-completion program; reading his methodology with care matters precisely because the HB shape match needs to be either (a) sharpened by adopting his apparatus, or (b) cleanly distinguished from his case if it cannot transfer.

---

## §1. Plain-language summary of Tapušković's setup

Tapušković studies the equal-mass two-loop **sunrise Feynman integral** as a motivic period. The sunrise graph $G$ has 2 vertices, 3 internal edges (the three "rays"), and 2 external legs carrying momentum $\pm q_1$. In $d = 2$ space-time dimensions and with all three internal masses equal to $m$, the integral
$$I_G(m, q) = \int_\sigma \frac{\Omega_G}{\Xi_G(m, q)}$$
(where $\Xi_G$ is the second Symanzik polynomial, $\Omega_G$ the standard projective volume form, and $\sigma$ the positive simplex in $\mathbb{P}^2(\mathbb{R})$) is known to evaluate to a combination of elliptic dilogarithms and integrals of modular forms (Bloch–Vanhove 2015). The elliptic curve $\mathcal{E}$ defined by $\Xi_G = 0$ is, in the equal-mass case, the universal elliptic curve over the modular curve $X_1(6)$. **This is why $\Gamma_1(6)$ appears, not $SL_2(\mathbb{Z})$**: the modular group is forced by the specific Feynman diagram's graph hypersurface, not chosen by hand.

**The categorical machinery.** Tapušković sidesteps the conjectural mixed motives category by working in the **Tannakian category of realizations** $\mathcal{H}(S)$ over a base $S$ (Brown 2017 "Notes on motivic periods"): objects are triples $V = (V_{\mathrm{B}}, V_{\mathrm{dR}}, c)$ of Betti cohomology, algebraic de Rham cohomology, and comparison isomorphism, typically given by the cohomology of a family of algebraic varieties over $S$. Motivic periods are equivalence classes $[\mathcal{V}, [\gamma], [\omega]]^{\mathrm{m}}$ where $[\gamma]$ is a section of $\mathcal{V}_{\mathrm{B}}^\vee$ and $[\omega]$ a section of $\mathcal{V}_{\mathrm{dR}}$. They form a ring $\mathcal{P}^{\mathrm{m}}_{\mathcal{H}(S)}$ equipped with a period homomorphism $\mathrm{per}$ to multivalued meromorphic functions.

The **cosmic Galois group** is the subquotient of the Tannaka group $G^{\mathrm{dR}}_{\mathcal{H}(S)}$ (with respect to a chosen fiber functor) that acts on the subring of $\mathcal{P}^{\mathrm{m}}$ spanned by motivic periods of "graph motives" $\mathrm{mot}_G$. Its action dualizes to a **coaction**
$$\Delta: \mathcal{P}^{\mathrm{m}} \to \mathcal{P}^{\mathrm{m}} \otimes \mathcal{P}^{\mathrm{dR}}$$
on motivic periods, which decomposes a motivic Feynman amplitude into a sum of "Galois conjugates" labeled by a de Rham basis of the underlying graph motive (Brown 2017 "Feynman amplitudes, coaction principle, and cosmic Galois group", arXiv:1512.06409).

**The graph motive of the sunrise.** Defined as the cohomology of a pair
$$\mathrm{mot}_G = H^2(P^G \setminus \mathcal{E},\ D \setminus \mathcal{E} \cap D),$$
where $P^G \to \mathbb{P}^2$ is an iterated blow-up at the three coordinate points $[1:0:0], [0:1:0], [0:0:1]$ (corresponding to *motic subgraphs* — pairs of edges of $G$), $\mathcal{E}$ is the strict transform of $\{\Xi_G = 0\}$, and $D$ is the total transform of the coordinate simplex. Geometrically: a smooth elliptic curve $\mathcal{E}$ meets six lines $D_{\pm 1}, D_{\pm 2}, D_{\pm 3}$ at six points $P_1, \ldots, P_6$ (three coordinate hyperplanes + three exceptional divisors); the graph motive lives on this configuration.

**The main first-half claim.** Tapušković computes the cosmic-Galois coaction on $I_G^{\mathrm{m}}$ explicitly. The coaction has **six terms** in the generic-mass case and **three terms** in the equal-mass case (where motivic logarithms vanish). The "small graphs principle" (Brown's conjecture, arXiv:1512.06409 §8.4) predicts conjugates should be motivic Feynman amplitudes of *subquotient* graphs — but for the sunrise, only *one* conjugate (the term $I^{\mathrm{m}}_{G \setminus e_3}$, the "bubble" graph obtained by deleting an edge) admits this reading. The remaining conjugates are expressed in terms of **graphs obtained by subdividing edges** — a new tool Tapušković introduces.

The second half of the paper develops the relative completion of $\pi_1$ on modular curves (Brown–Hain machinery) and uses it to reprove that the equal-mass sunrise is a $\mathbb{Q}$-linear combination of iterated Eichler integrals, periods of the elliptic curve, and powers of $2\pi i$. This is structurally an alternative proof of Bloch–Vanhove 2015 / Adams–Weinzierl 2018.

---

## §2. The edge-subdivision construction (§3 of Tapušković)

### §2.1 The construction itself

**Definition (Tapušković Def. 3.1).** Given a connected Feynman graph $G$ and an internal edge $e \in E_G$ joining vertices $v_{e,1}, v_{e,2}$, the **subdivided graph** $G_{s(e)}$ is obtained by deleting $e$ and adding a new vertex $v$ together with two new edges $e_1$ (joining $v_{e,1}$ to $v$) and $e_2$ (joining $v$ to $v_{e,2}$), with the new edges inheriting the mass $m_e$ of the deleted edge. For a multi-index $I = \{e_1^{k_1}, \ldots, e_{N_G}^{k_{N_G}}\}$, $G_{s(I)}$ subdivides edge $e_i$ a total of $k_i$ times. Concretely for the sunrise: $G_{s(e_1)}$ has 4 edges (one ray subdivided once), $G_{s(e_1, e_2^2)}$ has 6 edges, etc. See Tapušković's Figure 2.

**The key morphism (Tapušković Eq. 9 / §3).** With edges ordered so that $\alpha_{N_G}$ corresponds to the subdivided edge $e$ in $G$ and $\alpha_{N_G}, \alpha_{N_G+1}$ correspond to $e_1, e_2$ in $G_{s(e)}$, define
$$\rho: \mathbb{A}^1 \times \mathbb{P}^{N_G - 1} \to \mathbb{P}^{N_G}, \qquad (y, \alpha_1, \ldots, \alpha_{N_G}) \mapsto (\alpha_1, \ldots, y\alpha_{N_G}, (1-y)\alpha_{N_G}).$$
**Lemma 3.1.** $\rho^*(X_{\Xi_{G_{s(e)}}}) = X_{\Xi_G}$ and $\rho^*(X_{\Psi_{G_{s(e)}}}) = X_{\Psi_G}$ — the graph hypersurfaces pull back to themselves (proof: subdivision replaces $\alpha_{N_G}$ by $\alpha_{N_G} + \alpha_{N_G+1}$ in the Symanzik polynomials, and $\rho^*$ inverts this by sending $y\alpha + (1-y)\alpha \mapsto \alpha$). Moreover the Feynman integrand pulls back as
$$\rho^*(\omega_{G_{s(e)}}) = -\alpha_{N_G} \Psi_G^{z_1} \Xi_G^{z_2}\, \omega_G \wedge dy,$$
with exponents $z_1, z_2$ depending on the space-time dimensions $d_G, d_{G_{s(e)}}$.

**Proposition 3.3 + Corollary 3.4 (the headline result).** The morphism $\rho$ lifts to a morphism $\tilde\rho: \mathbb{A}^1 \times P^G \to P^{G_{s(e)}}$ of the blow-ups, which **induces a morphism of graph motives**
$$\tilde\rho^*: \mathrm{mot}_{G_{s(e)}} \longrightarrow \mathrm{mot}_G \otimes \mathbb{Q}(0).$$
Consequently, the motivic Feynman amplitude of the subdivided graph is a motivic period of the *original* graph's motive:
$$I^{\mathrm{m}}_{G_{s(I)}}(m, q) = \big[\mathrm{mot}_G,\ [\sigma_G],\ \big[\pi_G^*\big( (-1)^K \alpha_1^{k_1} \cdots \alpha_{N_G}^{k_{N_G}} \Psi_G^{z_1} \Xi_G^{z_2} \omega_G(m, q) \big)\big]\big]^{\mathrm{m}},$$
where $K = \sum k_i$. **This confirms Brown's conjecture (Tapušković Eq. 10):** in the functor $G \mapsto \mathcal{FP}^{\mathrm{m}}(G)$ from graphs to motivic Feynman amplitudes, $G \sim G'$ when $G'$ is obtained from $G$ by subdividing an edge (when $m = q = 0$ — Tapušković extends to the kinematic case).

### §2.2 The concrete sunrise example

For the sunrise $G$, four "subdivided cousins" appear in the coaction (Tapušković Fig. 2):
- $G_{s(e_1)}$: one ray subdivided once (4 edges), used in $d = 2$.
- $G_{s(e_1, e_2^2)}$: ray 1 subdivided once, ray 2 subdivided twice (6 edges), $d = 4$.
- $G_{s(e_1^2, e_3)}$: ray 1 subdivided twice, ray 3 once (6 edges), $d = 4$.
- $G_{s(e_2, e_3^2)}$: ray 2 once, ray 3 twice (6 edges), $d = 4$.

The de Rham realisation $(\mathrm{mot}_G)_{\mathrm{dR}}$ is 6-dimensional and Tapušković completes a basis using the Feynman integrands of these graphs (pulled back via $\tilde\rho$). Specifically (Proposition 4.3):

| Basis element | Form | Source |
|:---|:---|:---|
| $\phi_{\mathrm{dR}}(\nu_0)$ | weight-0, from bubble $G \setminus e_3$ | face map (subgraph deletion) |
| $\pi_G^*(\omega_G)$ | weight-3, holomorphic on $\mathcal{E}$ | Feynman integrand in $d = 2$ |
| $\eta_G$ | differential of second kind | pullback of $G_{s(e_1)}$ in $d = 2$ |
| $\nu_1$ | weight-4 | pullback of $G_{s(e_1, e_2^2)}$ in $d = 4$ |
| $\nu_2$ | weight-4 | pullback of $G_{s(e_1^2, e_3)}$ in $d = 4$ |
| $\nu_3$ | weight-4 | pullback of $G_{s(e_2, e_3^2)}$ in $d = 4$ |

Without edge subdivision, only $\phi_{\mathrm{dR}}(\nu_0)$ and $\pi_G^*(\omega_G)$ would be available from the "naive" sunrise data; the remaining four basis vectors are *only accessible* via the edge-subdivision pullback. **The construction's load-bearing function is to populate the de Rham basis of $\mathrm{mot}_G$ when subquotient-graph cohomology is too small to span it.**

### §2.3 What the construction is, structurally

In one sentence: **edge subdivision is a way to map motivic Feynman amplitudes of *richer* graphs (more edges) into motivic *periods* of the *same* underlying motive, by exploiting the fact that subdivision replaces $\alpha_e$ with $\alpha_{e_1} + \alpha_{e_2}$ in the Symanzik polynomials without changing the graph hypersurface up to a linear change of variables.** It is a tool for **enriching de Rham cohomology of a graph motive without changing the motive** — necessary precisely because for some graphs (the sunrise being the smallest example) the small-graphs principle alone does not suffice to express all coaction conjugates as motivic Feynman amplitudes of subquotient graphs.

---

## §3. Mapping to GeoVac concepts

### §3.1 What GeoVac has that could play the role of "graph motive"

GeoVac's substrate at cutoff $n_{\max}$ is the discrete spectral triple $(\mathcal{A}_{\mathrm{GV}}, \mathcal{H}_{\mathrm{GV}}, D_{\mathrm{GV}})$ on the Fock-projected $S^3$ packing graph (Paper 32 §III). The **closed-form periods** at each cutoff sit in three rings (Paper 55):
- $M_1 \subset \mathbb{Q}[\pi, \pi^{-1}]$ — Hopf-base measure ring.
- $M_2 \subset \bigoplus_k \pi^{2k}\cdot\mathbb{Q}$ — pure-Tate sub-ring of Fathizadeh–Marcolli mixed-Tate (smaller, no MZVs at SD level).
- $M_3 \subset \mathrm{MT}(\mathbb{Z}[i, 1/2])$ — cyclotomic-mixed-Tate at level $\le 4$.

The case-exhaustion theorem (Paper 32 §VIII) classifies every $\pi$ in any GeoVac observable as coming from $\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$ at $k \in \{0, 1, 2\}$. The master Mellin engine domain partition (MR-A/B/C, 2026-05-06) sharpens this: $k$ indexes *both* mechanism (M1/M2/M3) *and* the natural class of observable from which that mechanism's signature can be extracted.

### §3.2 Do GeoVac's $k \in \{0, 1, 2\}$ correspond to graph operations?

**Honest answer: NO, not in the direct sense one might hope.** The dictionary M1 ↔ vertices, M2 ↔ edges, M3 ↔ subdivisions does **not** hold, for three structural reasons:

**(i) The graph object is fundamentally different.** Tapušković's "graph" is a Feynman diagram — a *single* combinatorial object whose graph hypersurface defines an algebraic variety (the elliptic curve, in the sunrise case). GeoVac's "graph" is the **Fock-projected $S^3$ packing graph** at cutoff $n_{\max}$ — a *family* of finite combinatorial objects forming a pro-system (Paper 56 §sec:ps1). The Tapušković graph carries kinematic data $(m, q)$ as continuous parameters; the GeoVac graph has only the discrete labels $(n, \ell, m_\ell)$ and the cutoff $n_{\max}$. There is no GeoVac analog of the Symanzik polynomials $\Psi_G, \Xi_G$ because there are no internal edge parameters $\alpha_e$ in the Feynman-integral sense.

**(ii) The operator-order grading $k$ is over a single spectral triple, not over graph modifications.** $\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$ at $k = 0, 1, 2$ varies the *operator weight* applied to a *fixed* Dirac operator $D$ on a *fixed* substrate. Tapušković's edge subdivision varies the *graph* and pulls back differential forms via a $\mathbb{A}^1$-family of morphisms. The two operations live on different sides of the spectral-triple/Feynman-graph divide.

**(iii) The Mellin engine's three slots are not produced by an analog of subdivision.** M1 (Hopf-base measure $\pi$), M2 (Seeley–DeWitt $\pi^{2k}$), M3 (vertex-parity Hurwitz $\beta(2k)/G$) come from three *categorically distinct* domains of observables (state-space propinquity rates / heat-kernel spectral-action rates / vertex-restricted parity-character sums respectively). Edge subdivision produces *de Rham basis enrichment for a single motive*, all sitting in one category of motivic periods.

### §3.3 The partial structural parallel that does exist

There is, however, a genuine partial parallel at the level of **basis enrichment of a too-small de Rham realization**. Note the following analogy:

| Tapušković (sunrise) | GeoVac (Paper 56) |
|:---|:---|
| Small-graphs principle gives only $\phi_{\mathrm{dR}}(\nu_0)$ and $\pi_G^*(\omega_G)$ — 2 of 6 needed basis vectors. | Per-cutoff abelian factor $\mathbb{G}_a^{3 N(n_{\max})}$ acts on the natural-substrate panel — but the panel is fixed, not enriched. |
| Edge subdivision provides the remaining 4 basis vectors $\eta_G, \nu_1, \nu_2, \nu_3$ via $\tilde\rho^*$ morphisms. | The pro-system $\{\mathcal{H}_{\mathrm{GV}}(n_{\max})\}$ provides transition maps $P_{m,k}$ between cutoffs (Paper 56 Thm 4.2), enriching the substrate as $n_{\max} \to \infty$. |
| Subdivision $G \to G_{s(I)}$ is parametrised by multi-indices $I$ of non-negative integers. | Cutoff refinement $n_{\max} \to n_{\max} + 1$ is parametrised by a single non-negative integer. |

**The parallel that *might* be made precise:** if there is a Tannakian functor from GeoVac's pro-system of substrates to a category of "enriched" cohomologies analogous to Tapušković's edge-subdivision-enriched de Rham realizations, then cutoff refinement plays the role of edge subdivision. **This is suggestive but not concrete.** The Tapušković morphism $\tilde\rho^*$ has a clean geometric origin (blow-ups commute with the linearization map); the GeoVac pro-system transitions $P_{m,k}$ are closed-form 0/1 matrices forced by the packing axiom but without an analogous "geometric covering" interpretation.

### §3.4 What about the proposed dictionary M1 ↔ vertices, M2 ↔ edges, M3 ↔ subdivisions?

Direct verification against Tapušković's text:

- **Vertices.** In Tapušković, vertices play no role in the motivic coaction itself — the graph hypersurfaces are built from spanning trees (edges) and spanning 2-trees (cuts), not from vertex incidence directly. The number of motic subgraphs equals the number of "topologically essential" edge subsets, not vertices. There is no observable in Tapušković that scales like "number of vertices" and corresponds to $\pi$ from a Hopf-base measure. **M1 ↔ vertices is not supported by the source.**

- **Edges.** Edges enter through the integrand $\Omega_G$ and the Symanzik polynomials, and the small-graphs principle relates graphs by edge contraction/deletion. The number of edges $N_G$ does control the dimension of the de Rham realization (loosely). But Tapušković's edge-subdivision-induced basis vectors $\nu_1, \nu_2, \nu_3$ are *weight-4* objects in a *single graph motive's cohomology*, not "M2 transcendentals". The connection to $\pi^{2k}$ is absent. **M2 ↔ edges is not supported.**

- **Subdivisions.** Edge subdivision is structurally what Tapušković introduces; its role is basis enrichment, not transcendental generation. M3 ($\beta(2k)$ / Catalan $G$) lives in Brown's coaction *period output* once Eichler integrals of cusp forms are evaluated — but this is a level-$\Gamma_1(6)$ phenomenon specific to the equal-mass sunrise. **M3 ↔ subdivisions is not supported by Tapušković's text.**

**Honest verdict on the dictionary.** Tapušković's edge subdivision is a graph-side enrichment of the de Rham realization of a single graph motive. GeoVac's M1/M2/M3 partition is an operator-order grading on a single spectral triple's heat-kernel Mellin transforms. The two structures are not naturally identified. The proposed dictionary in the prompt (M1 ↔ vertices, M2 ↔ edges, M3 ↔ subdivisions) does not transfer.

---

## §4. The cosmic-Galois coaction on Tapušković's side

### §4.1 The coaction as written

Tapušković's main coaction formulae (Theorems 1.4, 1.5, Theorem 4.4 of the paper) read in the equal-mass case:
$$\Delta(I_G^{\mathrm{m}}) = I_G^{\mathrm{m}} \otimes K_1^{\mathrm{dR}} \mathbb{L}^{\mathrm{dR}} + I_{G_{s(e_1)}}^{\mathrm{m}} \otimes K_{2,\eta}^{\mathrm{dR}} \mathbb{L}^{\mathrm{dR}} + I_{G \setminus e_3}^{\mathrm{m}} \otimes I_{G, G \setminus e_3}^{\mathrm{dR}},$$
where $K_1^{\mathrm{dR}}, K_{2,\eta}^{\mathrm{dR}}$ are de Rham complete elliptic integrals of the first and second kind associated to the elliptic curve $\mathcal{E}$, $\mathbb{L}^{\mathrm{dR}}$ is the de Rham version of $2\pi i$, and $I_{G, G \setminus e_3}^{\mathrm{dR}}$ is a de Rham Feynman period built from a chosen weight-0 basis vector (the bubble graph $G \setminus e_3$). In the general-mass case (Eq. 11), three additional terms appear, each carrying an $F^{\mathrm{dR}}_{P_i, P_j}$ incomplete elliptic integral of the first kind (or, in the equivalent alternative basis Eq. 12, three motivic logarithms of mass ratios).

### §4.2 What the coaction structurally does

The cosmic-Galois coaction expresses an "uncertainty" about the motivic Feynman amplitude: a Galois conjugate is a way the amplitude could have been if the comparison isomorphism (Betti ↔ algebraic de Rham) were applied with respect to a different element of the Tannaka group. Each conjugate term has the shape (motivic period) ⊗ (de Rham period). The motivic-side factor identifies *which* graph (sub-quotient or edge-subdivided) the conjugate corresponds to; the de Rham factor identifies *which de Rham basis element* of the original graph motive it pairs with.

The small-graphs principle (Brown 2017 §8.4, 9.3) predicts that for "well-behaved" graphs, motivic-side factors are motivic Feynman amplitudes of *smaller* graphs (sub-quotients), so easy results for small graphs constrain large-graph amplitudes. Tapušković's edge-subdivision lemma extends this: when subquotient amplitudes don't suffice, **edge-subdivided amplitudes plug the gap** — and crucially, edge-subdivided amplitudes are still motivic periods of the *original* graph's motive (Cor. 3.4), so they are computable without leaving the graph motive's cohomology.

### §4.3 Numerical fingerprint we could check on GeoVac

**Honest assessment: there is no direct numerical fingerprint of the Tapušković coaction that would test against GeoVac data, because the categorical structures are different.** What Tapušković computes is a coaction whose period output, when evaluated on the equal-mass sunrise integral at numerical kinematics, would produce a polynomial of bounded depth in:

- $K_1, K_{2,\eta}$: complete elliptic integrals of first/second kind at $\tau = $ modular point of $X_1(6)$.
- $\mathbb{L} = 2\pi i$.
- $I_{G \setminus e_3}$: a logarithm of a kinematic ratio.

For comparison with GeoVac, the cleanest plausible test would be:

**Test T-FP1.** Take GeoVac's M3 sector outputs (where $L$-values like $G = \beta(2)$ and $\beta(4)$ appear) at a *fixed* representative cutoff (e.g., $n_{\max} = 3$ or 4). Form depth-2 combinations (joint Mellin of two M3 outputs). PSLQ against the Hain–Brown / Tapušković modular ring including $\{E_4(\tau), E_6(\tau), \Delta(\tau), L(\Delta, k)\}$ at the modular point $\tau$ of $X_1(6)$, **for several candidate $\tau$ on the imaginary axis**. Today's HB PSLQ Test A used $\tau = 2i$ on $SL_2(\mathbb{Z})$. **A genuine test of the Tapušković parallel would instead use $X_1(6)$-specific modular points and the $\Gamma_1(6)$-specific period basis** (which Tapušković's Theorem 5.1 explicitly characterizes: the periods of $\mathcal{O}(\pi^{\mathrm{rel}}_1(X_1(6), x, y))$ are products of iterated Eichler integrals, powers of periods of $\mathcal{E}$, and powers of $2\pi i$, all evaluated on $X_1(6)$, not $\mathcal{M}_{1,1}$).

**Honest scope of this test.** Even with a perfect $X_1(6)$ modular basis, the test would still be a *blind PSLQ search* unless GeoVac has independent structural reason to expect cusp-form $L$-values to appear. Today's HB negative on $SL_2(\mathbb{Z})$ already drops the prior for any modular identification by a substantial amount. The Tapušković extension to $\Gamma_1(6)$ would be informative only if there were a structural argument that $\Gamma_1(6)$ is more natural for GeoVac than $SL_2(\mathbb{Z})$ — and **there is none in the GeoVac corpus today** ($SL_2$ in GeoVac is forced by Bertrand × Fock on $S^3 = SU(2)$; $\Gamma_1(6)$ would require an additional discrete structural input that has no GeoVac-internal motivation).

### §4.4 What the Tapušković paper does NOT have that would matter

The paper does not contain:
- A general theorem identifying which Feynman graphs produce which $\Gamma_0(N)$ or $\Gamma_1(N)$ modular content.
- An explicit numerical evaluation of the equal-mass sunrise's coaction conjugates at specific kinematic points.
- A PSLQ-ready basis of $\Gamma_1(6)$ periods to multi-digit precision.

These absences mean **even if we wanted to test GeoVac against the Tapušković coaction directly, we would first have to build the period-basis numerical infrastructure ourselves** — and that infrastructure does not exist in the literature in deployable form.

---

## §5. Verdict against the decision gate

**METHODOLOGY-VALUABLE-NO-DIRECT-PARALLEL** (with structural elaboration).

The Tapušković methodology is informative for GeoVac community framing in three ways:

1. **It is the published precedent for adopting relative-completion machinery on the physics side without claiming literal modular-form identification.** Tapušković works with motivic periods, the cosmic-Galois coaction, and the relative completion of $\Gamma_1(6)$ — and the equal-mass sunrise is *not* claimed to be a modular form itself; it is a period of a motive sitting in a Tannakian category that *also* contains modular periods. This is precisely the framing GeoVac would adopt: $\mathcal{H}_{\mathrm{GV}}$ as a Tannakian Hopf algebra, $U^*_{\mathrm{GV}}$ as a candidate matching the shape of a relative completion's pro-unipotent radical, with the understanding that the period content need not literally be modular.

2. **It demonstrates that physics-driven motives can require de Rham basis enrichment beyond the small-graphs principle.** GeoVac's pro-system structure (Paper 56) is plausibly playing an analogous role, though through cutoff refinement rather than edge subdivision. Naming this structural parallel honestly clarifies what GeoVac's pro-system is doing.

3. **It establishes that physics-side input into the relative-completion program is publishable and informative even without literal identification.** Tapušković's paper is in Communications in Number Theory and Physics — the joint math-physics venue. GeoVac's reciprocity contributions (per the HB adoption survey §3 R1–R5) would target the same audience.

**But the methodology does NOT provide a direct adoption candidate** because:

- The edge-subdivision construction is graph-hypersurface-specific (relies on Symanzik polynomial structure), and GeoVac has no graph-hypersurface analog.
- The coaction's period output is on $\Gamma_1(6)$, which has no GeoVac-internal motivation.
- The dictionary M1/M2/M3 ↔ graph operations (vertices/edges/subdivisions) is not supported by the source.

**Recommendation.** Do **not** propose a sprint-scale adoption test based on Tapušković. Instead, cite Tapušković in Paper 55 §6 and/or Paper 56 §sec:open as the published precedent for the relative-completion framing, and flag that GeoVac's path forward is the $\mathcal{G}_4$ comparison target (Paper 56 §sec:open_g4) rather than a modular identification — consistent with today's HB PSLQ negative.

---

## §6. Sprint outline — N/A under METHODOLOGY-VALUABLE-NO-DIRECT-PARALLEL verdict

No sprint outline is proposed. Per §5 above, the next sprint-scale move on this thread is **editorial citation and framing**, not a computational adoption test. Specifically (suggestions only — no edits applied in this read-only sprint):

- **Citation pass.** Add Tapušković 2023 (arXiv:2303.17534, Comm. Num. Theor. Phys. 18 (2024) no. 2) as a bibitem in Paper 55 § Open targets and/or Paper 56 §sec:open with framing: *"The first physics-side input into the relative-completion program (Hain–Brown lineage) is Tapušković's 2023 computation of the cosmic-Galois coaction on the equal-mass sunrise Feynman integral, expressing it as a period of the relative completion of $\Gamma_1(6)$. Our shape match $U^*_{\mathrm{GV}} = \mathbb{G}_a^\infty \rtimes SL_2$ to the relative completion's pro-unipotent structure is analogous in form but the period content differs structurally: GeoVac's master Mellin engine produces $M_1 \cup M_2 \cup M_3$ periods (Hopf-base, Seeley–DeWitt, level-4 cyclotomic Dirichlet) without modular-form input, while Tapušković's coaction outputs iterated Eichler integrals on $X_1(6)$."*

- **HB-attribution bug fix.** Memory file `memory/hain_brown_identification.md` and the v3.78.0 strategic synthesis memo currently miscite Tapušković as "Bouillon 2023". A sibling sprint is already addressing this per the prompt's read-only constraint note.

- **No code, no paper math, no memory file edits in this sprint.**

---

## §7. Honest limitations of this read

- The PDF text was extracted via a single end-to-end PDF parse; minor character-level OCR artifacts in displayed equations are possible. All structural claims (theorem statements, definitions, diagram structures, decomposition counts) were cross-checked across the introduction, §2, §3, §4, §5 of the paper.
- I did not independently verify Tapušković's bibliographic citations [1]–[30] beyond confirming that Bloch–Vanhove 2015 [6], Adams–Weinzierl 2018 [3], and Brown's three foundational papers [9, 10, 11] correspond to identifiable arXiv preprints. The Hain–Matsumoto 2015 paper, which the HB adoption survey relies on, is *not* in Tapušković's bibliography — Tapušković uses Hain 1998 [20] and Hain 2016 [21] instead.
- The methodology assessment in §3–§5 is based on reading the paper alone; it does not include a separate literature audit of post-2023 follow-ups to Tapušković's work, except the Kleinschmidt et al. 2026 paper noted in the HB adoption survey §1.
- Final-verdict gate ("METHODOLOGY-VALUABLE-NO-DIRECT-PARALLEL") is conservative; a more permissive read might argue that the edge-subdivision-as-de-Rham-basis-enrichment analogy with GeoVac's pro-system warrants its own structural sprint. I do not propose this because the analogy is suggestive rather than concrete and would require multi-month construction work without a clear falsifier.

---

**End of memo. ~3,900 words.**
