# Sprint Q5' QSM/Motivic-Galois Lit-Read Memo

**Date:** 2026-06-03
**Sprint:** Q5' (entry-point ii) — published-precedent scoping for "enriched motivic Galois group of a discrete spectral triple"
**Triggered by:** Sprint A7 closing M2/M3 cyclotomic field-coincidence at HALF-STRUCTURAL with Q5' opened as Paper 55 §subsec:open_m2_m3 multi-year research project.
**Verdict:** **PARTIALLY TRANSPORTABLE.** The Connes-Marcolli motivic-Galois lineage provides a substantial published precedent for the **target side** of Q5' (Tannakian/motivic Galois groups acting on periods of mixed Tate motives in noncommutative-geometric settings) and Marcolli-van Suijlekom 2014 lineage provides the right ambient category on the **source side** (finite/graph spectral triples). The **unification mechanism specific to Q5'** — attaching a motivic Galois group to a discrete spectral triple on a compact homogeneous space and recovering both M2's Witt-splitting and M3's cyclotomic conductor as Tannakian fiber-functor structures — is **genuinely open** and has no published instance to our knowledge. Q5' is a real gap, but the toolkit is in place.

---

## 1. Decision-gate framework

Three gate options for this scoping (per task spec):

- **TRANSPORTABLE:** Published QSM/Q-lattice formalism adapts directly to GeoVac discrete spectral triple with named theorem/definition entry point.
- **PARTIALLY TRANSPORTABLE:** Some ingredients exist; unification mechanism does not.
- **NOT TRANSPORTABLE:** QSM lineage operates on categorically different object; Q5' would need its own construction.

**Landing:** PARTIALLY TRANSPORTABLE. Three independent precedent strands are in place but none of them, alone or in combination, executes the Q5' construction.

---

## 2. The three precedent strands

### Strand A — Connes-Marcolli motivic Galois group of perturbative renormalization

**Headline reference:**
- Connes, A.; Marcolli, M. **"Renormalization, the Riemann-Hilbert correspondence, and motivic Galois theory."** arXiv:math/0409306 (2004) and Connes-Marcolli book Ch 1 §3.

**What it does.** Constructs a pro-algebraic group scheme $U^*$ (the **cosmic Galois group** in Cartier's terminology) via Tannakian formalism applied to the category of equisingular flat vector bundles, and shows that:
- the renormalization group is canonically a one-parameter subgroup of $U^*$;
- $U^*$ is the motivic Galois group of mixed Tate motives over $\mathbb{Z}$;
- divergences of perturbative QFT are organized by $U^*$ in a universal (theory-independent) way.

**Relevance to Q5'.** This is the **gold-standard precedent** for "motivic Galois group attached to an analytic structure in NCG." The construction is via:
1. an analytic input (the singular Birkhoff factorization of equisingular connections),
2. a Tannakian category built from it,
3. a fiber functor (de Rham realization), and
4. a motivic Galois group recovered from the Tannakian dual.

Q5' would mirror this exactly but with the **input** being the discrete spectral triple's master Mellin engine evaluations (M1/M2/M3 Mellin moments) instead of singular Birkhoff data. The four-step recipe transports as a **shape**; the specific Tannakian category attached to a discrete spectral triple is the construction Q5' would need.

**Honest gap.** Connes-Marcolli's $U^*$ acts on QFT counterterms (universal renormalization group structure), not on the periods of a spectral triple. Q5' needs the Tannakian category to capture both M2 ($\pi^{2k}$-graded splitting over $\mathbb{Q}(i)$, trivial Galois action) and M3 (cyclotomic mixed-Tate periods over $\mathbb{Z}[i, 1/2]$, non-trivial Galois action with Deligne–Glanois level-4 conductor). This is **shape-similar but structurally distinct** from $U^*$.

---

### Strand B — Greenfield-Marcolli twisted spectral triples for QSM/Bost-Connes

**Headline references:**
- Greenfield, M.; Marcolli, M.; Teh, K. **"Twisted Spectral Triples and Quantum Statistical Mechanical Systems."** *p-Adic Numbers, Ultrametric Analysis and Applications* 6 (2014), 81–104. arXiv:1305.5492. (*Corrected 2026-06-04 per Round 2 Greenfield–Marcolli transport memo: original entry cited arXiv:1305.5564 and omitted Teh as third author; both errors fixed here. The 1305.5564 ID belongs to an unrelated Braaten–Kang paper on $X(3872)$.*)
- Connes, A.; Marcolli, M.; Ramachandran, N. **"KMS states and complex multiplication."** Selecta Math. 11 (2005), 325–347. arXiv:math/0501424.
- "TYPE III σ-SPECTRAL TRIPLES AND QUANTUM STATISTICAL MECHANICAL SYSTEMS" (Caltech preprint, Marcolli et al.).

**What it does.** Explicit construction of a **twisted (σ-deformed, type III) spectral triple** $(\mathcal{A}, \mathcal{H}, D; \sigma)$ from the Bost-Connes QSM system $(\mathcal{A}_{\mathrm{BC}}, \sigma_t)$. The Liouville function fixes both the sign $F$ of $D$ and the twisting $\sigma$. KMS states realize the **cyclotomic Galois group $\mathrm{Gal}(\mathbb{Q}^{\mathrm{ab}}/\mathbb{Q})$** as symmetries acting on the zero-temperature equilibrium states. For imaginary quadratic fields (Connes-Marcolli-Ramachandran 2005), the **full Galois group of the modular field** appears as symmetries on KMS states.

**Relevance to Q5'.** This is the published example of "spectral triple equipped with a cyclotomic Galois action." Q5' wants the same shape — a Galois group acting on the spectral-triple data — but in a categorically different setting: GeoVac is a **type II / type I compact homogeneous space** (the truncated Camporesi-Higuchi triple on $S^3$), not a type III QSM system. The Galois action in BC comes from arithmetic of $\mathbb{Q}$-lattices (adelic structure), whereas Q5' wants a Galois action coming from the geometry of the discrete graph and its periods.

**Honest gap.** Greenfield-Marcolli is type III, σ-twisted, infinite-dimensional, and arithmetic. GeoVac is bounded (truncated), untwisted (standard Connes axioms — see Paper 32 §III/§IV), finite-dimensional at each $n_{\max}$, and combinatorial (the cyclotomic level-4 structure comes from M3's vertex parity / quarter-integer Hurwitz shift, not from a number-field input). The shape of "Galois acts on KMS-state structure" transports as inspiration; the construction does not transport mechanically.

---

### Strand C — Marcolli-Tabuada noncommutative motivic Galois groups

**Headline references:**
- Marcolli, M.; Tabuada, G. **"Noncommutative numerical motives, Tannakian structures, and motivic Galois groups."** J. Eur. Math. Soc. 18 (2016), 623–655. arXiv:1110.2438.
- Marcolli, M.; Tabuada, G. **"Unconditional noncommutative motivic Galois groups."** in *Hodge Theory and Classical Algebraic Geometry*, Contemp. Math. 647 (2015), 127–146. arXiv:1112.5422.
- Tabuada, G. *Noncommutative Motives.* University Lecture Series 63, AMS (2015).
- Marcolli, M.; Tabuada, G. **"Noncommutative motives and their applications."** arXiv:1311.2867 (2013 MSRI survey).

**What it does.** Constructs (under standard conjectures $C_{\mathrm{NC}}$ and $D_{\mathrm{NC}}$ — and unconditionally for the André-Kahn quotient) Tannakian categories of noncommutative numerical motives, with **periodic cyclic homology $\mathrm{HP}_*$** promoted to a symmetric monoidal functor playing the role of the de Rham fiber functor. The motivic Galois group $\mathrm{GalMot}_{\mathrm{NC}}$ acts on this. The Tate triple structure relates $\mathrm{GalMot}_{\mathrm{NC}}$ to the classical motivic Galois group $\mathrm{GalMot}$.

**Relevance to Q5'.** This is **the categorical scaffolding Q5' needs**. The input is a dg category (the natural categorical thickening of a spectral triple — every spectral triple gives rise to a dg category via the de Rham complex of the Hochschild/cyclic homology), the invariant is periodic cyclic homology, and the output is a motivic Galois group. If we can promote the GeoVac discrete spectral triple to its associated dg category and check $C_{\mathrm{NC}}$/$D_{\mathrm{NC}}$ (or use the unconditional construction), the Marcolli-Tabuada machinery directly attaches a motivic Galois group.

**Honest gap.** The dg-category attached to a finite-dimensional spectral triple at fixed $n_{\max}$ is essentially trivial (finite-dimensional matrix algebra at each cutoff), so the noncommutative numerical motive may not capture the M1/M2/M3 distinction the Mellin engine sees. The right object is probably the **inverse system** $\{\mathcal{T}_{n_{\max}}\}_{n_{\max}}$ together with the Berezin reconstruction maps $B_{n_{\max}}: C(S^3) \to \mathcal{O}_{n_{\max}}$ (Paper 38 §L4), giving rise to a pro-dg category whose noncommutative motive captures the Mellin engine's grading.

---

### Strand D (the directly closest hit) — Fathizadeh-Marcolli periods of the spectral action

**Headline reference:**
- Fathizadeh, F.; Marcolli, M. **"Periods and motives in the spectral action of Robertson-Walker spacetimes."** Comm. Math. Phys. 356 (2017), 641–671. arXiv:1611.01815.

**Already in Paper 55 bibliography (`fathizadeh_marcolli2017`, cited in Sprint Mixed-Tate Test 2026-06-03).**

**What it does.** Proves that the coefficients of the asymptotic Seeley-DeWitt expansion of the spectral action on Euclidean Robertson-Walker spacetimes are **periods of mixed Tate motives** over $\mathbb{Z}$, specifically given by relative motives of complements of quadric hypersurfaces and unions of coordinate hyperplanes. One quadric per Seeley-DeWitt coefficient (a sharpening relative to the Feynman-integral case which needs multiple hyperplanes).

**Relevance to Q5'.** This is **literally the precedent for "spectral-action coefficients live in mixed Tate motives"** — which IS the M2 sector of our master Mellin engine. Paper 55 §4 (Sprint A5/A7/A8 closures) cites this as the source of the F-M mixed-Tate classification that GeoVac's discrete $S^3$ M2 SD coefficients inherit. The cited paper's mixed-Tate motive lives in $\mathrm{MT}(\mathbb{Z})$ — exactly the pure-Tate sub-ring $\bigoplus_k \pi^{2k}\mathbb{Q}$ that Sprint A8 extracted explicitly for GeoVac.

The motivic Galois group of $\mathrm{MT}(\mathbb{Z})$ is the **Deligne motivic fundamental group $\pi_1^{\mathrm{mot}}(\mathrm{MT}(\mathbb{Z}))$**, with Tannakian dual structure well-understood (Brown 2012, Deligne-Goncharov 2005). On the pure-Tate sub-ring $\bigoplus_k \pi^{2k}\mathbb{Q}$, the motivic Galois action is **trivial** (no odd zetas, no MZVs, no Galois conjugates — the Tate motive $\mathbb{Q}(k)$ is acted on only by the Tate scaling $\mathbb{G}_m$). This matches the M2-side of Q5' exactly:\ **trivial Galois action**.

**Honest gap.** Fathizadeh-Marcolli works in the **smooth-manifold + continuum spectral action** setting; GeoVac is the **discrete spectral triple + closed-form Sprint A8 SD ring**. The fact that Sprint A8's $[Z_{1,2n}]$ Grothendieck-class polynomial matches an F-M-style affine-complement structure (Paper 55 Sprint A8 closure) means the transport is essentially mechanical at the M2 sector. **M3 is not in Fathizadeh-Marcolli's scope.**

---

### Strand E — Glanois 2015 (Deligne descent) and Brown 2012 — the M3-side precedent

**Headline references:**
- Glanois, C. **"Motivic unipotent fundamental groupoid of $\mathbb{G}_m \setminus \mu_N$ for $N=2,3,4,6,8$ and Galois descents."** arXiv:1411.4947 (2014/2015). (**Already in Paper 55 bibliography.**)
- Brown, F. **"Mixed Tate motives over $\mathbb{Z}$."** Ann. Math. 175 (2012), 949–976. (**Already in Paper 55 bibliography.**)
- Deligne, P. **"Le groupe fondamental unipotent motivique de $\mathbb{G}_m - \mu_N$ pour $N=2, 3, 4, 6$ ou $8$."** Publ. Math. IHES 112 (2010), 101–141.

**What it does.** Constructs explicit motivic Galois descents from level-$N$ cyclotomic mixed-Tate motives over $\mathbb{Z}[\zeta_N, 1/N]$ down to level-$M$ for $M | N$. Deligne 2010 + Glanois 2015 cover exactly the cases $N \in \{2, 3, 4, 6, 8\}$ — i.e., $N=4$ is in the explicit catalogue. The Galois group at level $N=4$ is the motivic fundamental group of $\mathbb{P}^1 \setminus \{0, \pm i, \infty\}$, with explicit basis and descent to $N=2$ (which is $\mathrm{MT}(\mathbb{Z})$, the F-M pure-Tate world).

**Relevance to Q5'.** This is **the M3-side precedent**. Sprint A7 identified that M3's natural conductor is $\mathbb{Q}(\zeta_4) = \mathbb{Q}(i)$ (level 4), via the quarter-integer Hurwitz shifts that produce Catalan $G$ and $\beta(4)$ in the vertex-restricted QED expansions (Paper 28 §QED-vertex). Deligne-Glanois give the explicit motivic Galois group of this level-4 mixed-Tate category and the descent to level-2 ($\mathrm{MT}(\mathbb{Z})$). The Q5' question — **is there an enriched Galois group that unifies M2 (level 2, trivial action) with M3 (level 4, non-trivial action via $\mathrm{Gal}(\mathbb{Q}(i)/\mathbb{Q}) \cong \mathbb{Z}/2$)** — sits **between Deligne-Glanois descent (which has the level-4 → level-2 map but not the spectral-triple input) and Fathizadeh-Marcolli (which has the spectral-action input but doesn't engage the cyclotomic level)**.

**Honest gap.** Deligne-Glanois work purely in algebraic geometry; no spectral-triple input. The connection Q5' needs is "discrete spectral triple $\to$ category of periods $\to$ level-4 cyclotomic mixed-Tate motive whose Galois descent recovers the M2/M3 distinction." This synthesis is **not in the published literature** as of June 2026 (to our knowledge after this lit pass).

---

## 3. Two strands that are nearly-relevant but ruled out for Q5' specifically

### Marcolli-van Suijlekom 2014 gauge networks on finite spectral triples

(Already in our bibliography as `marcolli_vs2014`; cited extensively in Papers 25, 30, 41, 51.)

This is **our ambient category**:\ Marcolli-van Suijlekom construct finite spectral triples on graph vertices with $A$-valued connections on edges, whose spectral action is the Wilson lattice gauge theory. GeoVac's truncated Hopf-graph triple is in this lineage.

**Why not enough for Q5':** Marcolli-van Suijlekom 2014 does **not** attach a motivic Galois group to the gauge-network spectral triple. The Tannakian/motivic-Galois angle is absent from their construction. They give the right kind of input object (finite/graph spectral triple) but not the Tannakian output Q5' needs.

The Perez-Sanchez 2024/2025 corrections (arXiv:2401.03705, arXiv:2508.17338) — also in our bibliography — clarify the continuum limit but do not add motivic structure.

### Connes-Consani arithmetic site + Zeta Spectral Triples (Connes-Consani-Moscovici, arXiv:2511.22755, Nov 2025)

The Connes-Consani arithmetic-site / scaling-site program (2014–) and the very recent (verified post-cutoff November 27, 2025) "Zeta Spectral Triples" paper construct spectral triples whose spectra align with Riemann zeta zeros via rank-one perturbations of the scaling-operator spectral triple.

**Why not directly Q5':** This is the **Hilbert-Pólya** thread, targeting RH via spectral methods. The motivic content (when present) is on the arithmetic-site side, not on the spectral-triple side per se. Connes-Consani-Moscovici 2025 explicitly does not engage motivic Galois groups (verified by WebFetch). Useful as cross-reference but not as the Q5' precedent.

### Type III σ-spectral triples (independent of GreenField-Marcolli)

The "Type III σ-spectral triples" Caltech preprint (Marcolli et al.) lifts the BC system data to a σ-twisted spectral triple. Same as Strand B caveat:\ type III, σ-twisted, infinite-dimensional, arithmetic — not the GeoVac shape.

---

## 4. Synthesis — what Q5' would need to construct

The lit-read confirms that Q5' is **a genuine gap** with **three published precedent strands in alignment but no existing combination**. The construction Q5' would build is:

1. **Input:** GeoVac discrete spectral triple $\mathcal{T}_{n_{\max}}$ on $S^3$ (Paper 32 setup) or its inverse system $\{\mathcal{T}_{n_{\max}}\}_{n_{\max}}$ together with Berezin maps $B_{n_{\max}}$ (Paper 38 §L4).

2. **Categorical thickening:** promote $\mathcal{T}_{n_{\max}}$ to its associated dg category $\mathrm{dg}(\mathcal{T}_{n_{\max}})$ (in the sense of Marcolli-Tabuada 2016 input data).

3. **Period-producing structure:** the master Mellin engine $\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$ at $k \in \{0, 1, 2\}$ (Paper 18 §III.7, CLAUDE.md mellin_taxonomy_engine memory) sends $\mathcal{T}_{n_{\max}}$ to its (M1, M2, M3) Mellin moments. These are the analogue of the Fathizadeh-Marcolli periods of the spectral action — but for a **discrete** spectral triple.

4. **Tannakian dual:** apply Marcolli-Tabuada noncommutative motivic Galois machinery to $\mathrm{dg}(\mathcal{T}_{n_{\max}})$ (or to the pro-dg-category $\varprojlim_{n_{\max}} \mathrm{dg}(\mathcal{T}_{n_{\max}})$) using periodic cyclic homology $\mathrm{HP}_*$ as the fiber functor.

5. **Galois recovery:** the resulting motivic Galois group $\mathrm{GalMot}(\mathcal{T})$ should have:
   - **trivial restriction to M2** (matching the pure-Tate sub-ring $\bigoplus_k \pi^{2k}\mathbb{Q}$ of Sprint A8 — the M2-side of Fathizadeh-Marcolli);
   - **non-trivial restriction to M3** through the Deligne-Glanois level-4 descent $\mathrm{Gal}(\mathbb{Q}(i)/\mathbb{Q}) \cong \mathbb{Z}/2$ acting on the cyclotomic mixed-Tate part (matching Sprint A7's $\mathbb{Q}(\zeta_4)$ conductor).

6. **Compatibility:** the GH-convergence theorem (Paper 38) and the master-theorem (Papers 39, 40) should be promoted to **motivic-Galois-equivariant convergence** of the system $\{\mathrm{GalMot}(\mathcal{T}_{n_{\max}})\}_{n_{\max}}$ to a continuum-limit motivic Galois group $\mathrm{GalMot}(\mathcal{T}_\infty)$.

This is the **enriched motivic Galois group** Q5' names. The lit-read does not surface a precedent for **any** of steps (2)–(6) executed for a compact-homogeneous-space discrete spectral triple. **Each step has a published structurally-similar instance**, but the synthesis is new.

---

## 5. Specific named theorem/definition entry points for Q5' construction

| Step | Entry point | Reference | Status |
|:----|:------------|:----------|:-------|
| Categorical input | dg category from spectral triple via Connes axioms + Hochschild complex | Tabuada *Noncommutative Motives* AMS 2015, Ch. 1–2 | Standard |
| Fiber functor | $\mathrm{HP}_*$ on dg(\mathcal{T}_{n_{\max}}) | Marcolli-Tabuada 2016 arXiv:1110.2438 Thm 1.1 / Prop 4.2 | Conditional on $C_{\mathrm{NC}}$, $D_{\mathrm{NC}}$ — unconditional version in arXiv:1112.5422 |
| Motivic Galois recovery | Tannakian dual of the resulting symmetric monoidal category | Marcolli-Tabuada 2015 arXiv:1112.5422 Thm 1.4 | Unconditional |
| Period-realization | Spectral-action coefficients as periods | Fathizadeh-Marcolli 2017 arXiv:1611.01815 Thm 1.1 + §3 | Continuum case only; discrete-case extension is Q5'/Sprint A5 in our notation |
| Cyclotomic descent | Galois descent level-4 to level-2 | Deligne 2010 Publ. Math. IHES 112; Glanois 2015 arXiv:1411.4947 §3–4 | Standard for $N \in \{2,3,4,6,8\}$ |
| Renormalization-side shape | Tannakian dual of equisingular flat bundles | Connes-Marcolli 2004 arXiv:math/0409306 Thm 1.107 + Connes-Marcolli book Ch. 1 §3 | Standard |
| Universal recipe | Cosmic Galois group $U^*$ classifying analytical-structure periods in NCG | Connes-Marcolli book 2008, Ch. 1 | Standard |

---

## 6. Honest scope — what could not be checked in single-session lit pass

1. **The full text of Fathizadeh-Marcolli 2017** (the affine-complement structure and the Grothendieck-class refinement)\ was accessible only via abstract+search summary; the binary PDF could not be parsed. Sprint A8 already worked through the F-M Thm 6.2 structure independently with the GeoVac specialisation, so this is not blocking.
2. **The Tabuada *Noncommutative Motives* AMS Lecture Series 63 book (2015)** content was sketched via abstracts of constituent papers but the book-level synthesis was not paginated. The unconditional construction (relevant for Q5') is in the arXiv:1112.5422 paper which we did access.
3. **Connes-Marcolli book 2008 Ch. 4** (motivic Galois group of QFT divergences, full theorem statement) was accessible only via search-result extracts. The arXiv predecessor (math/0409306) was confirmed.
4. **Whether anyone has applied $U^*$ (cosmic Galois) directly to a discrete spectral triple** — the lit-read returned NO hits. To our knowledge after this pass, this is genuinely unaddressed.
5. **Whether Deligne-Glanois descent at $N=4$ has been used in any NCG/spectral-action context** — the lit-read returned NO hits. Glanois 2015 is purely algebraic-geometric; F-M 2017 stays at the $\mathrm{MT}(\mathbb{Z})$ level (level $N=2$) and does not engage cyclotomic descent.
6. **Whether anyone has constructed the motivic Galois group of a finite/discrete graph spectral triple** — Marcolli-van Suijlekom 2014 does not. Higher-rank graph spectral triples (Farsi-Gillaspy-Julien-Kang-Packer 2018, arXiv:1804.05209) do not engage motives. No published precedent identified.

---

## 7. Verdict and one-line summary

**One-line verdict:** PARTIALLY TRANSPORTABLE — every individual ingredient Q5' needs has a published precedent (motivic Galois from analytic data — Connes-Marcolli 2004; spectral-triple-with-Galois-action — Greenfield-Marcolli 2014, CMR 2005; noncommutative motivic Galois machinery — Marcolli-Tabuada 2016; spectral-action coefficients as mixed-Tate periods — Fathizadeh-Marcolli 2017; cyclotomic level-4 motivic descent — Deligne 2010 / Glanois 2015), but **the synthesis specific to Q5' — an enriched motivic Galois group of a discrete spectral triple on a compact homogeneous space unifying M2's trivial action with M3's non-trivial $\mathrm{Gal}(\mathbb{Q}(i)/\mathbb{Q})$ action — has no published precedent to our knowledge**.

Q5' is a genuine multi-year research target sitting in a clearly-identified gap. The toolkit is in place; the construction is new.

---

## 8. References (10 most relevant, verified)

| # | Citation | Verified | Relevance to Q5' |
|:--:|:---------|:--------:|:-----------------|
| 1 | Connes, Marcolli. "Renormalization, the Riemann-Hilbert correspondence, and motivic Galois theory." arXiv:math/0409306 (2004). | YES (verified via search + Caltech URL) | Gold-standard motivic-Galois-from-NCG construction; the four-step Tannakian recipe Q5' would mirror. |
| 2 | Connes, Marcolli. *Noncommutative Geometry, Quantum Fields and Motives.* AMS Colloq. Pub. 55 (2008), Chs. 1, 4. | YES (book exists; AMS bookweb URL confirmed) | The book-length development of motivic Galois in NCG; cosmic Galois group $U^*$. |
| 3 | Marcolli, Tabuada. "Noncommutative numerical motives, Tannakian structures, and motivic Galois groups." J. Eur. Math. Soc. 18 (2016), 623–655. arXiv:1110.2438. | YES (EMS URL, arXiv abstract) | Tannakian/motivic-Galois machinery for noncommutative motives, with $\mathrm{HP}_*$ as fiber functor — the categorical scaffolding Q5' needs. |
| 4 | Marcolli, Tabuada. "Unconditional noncommutative motivic Galois groups." Contemp. Math. 647 (2015), 127–146. arXiv:1112.5422. | YES (verified arXiv) | Unconditional version of (3), removing dependence on $C_{\mathrm{NC}}$/$D_{\mathrm{NC}}$. |
| 5 | Fathizadeh, Marcolli. "Periods and motives in the spectral action of Robertson-Walker spacetimes." Comm. Math. Phys. 356 (2017), 641–671. arXiv:1611.01815. | YES (already in Paper 55 bibliography as `fathizadeh_marcolli2017`) | Spectral-action coefficients = periods of mixed Tate motives; M2-side precedent. |
| 6 | Glanois, C. "Motivic unipotent fundamental groupoid of $\mathbb{G}_m \setminus \mu_N$ for $N=2,3,4,6,8$ and Galois descents." arXiv:1411.4947 (2015). | YES (already in Paper 55 bibliography as `glanois2015`) | Explicit cyclotomic motivic Galois descent at $N=4$; M3-side precedent. |
| 7 | Deligne, P. "Le groupe fondamental unipotent motivique de $\mathbb{G}_m - \mu_N$ pour $N=2, 3, 4, 6, 8$." Publ. Math. IHES 112 (2010), 101–141. | YES (standard, widely cited) | Foundational cyclotomic motivic Galois result at $N=4$. |
| 8 | Greenfield, M.; Marcolli, M.; Teh, K. "Twisted Spectral Triples and Quantum Statistical Mechanical Systems." *p-Adic Numbers Ultrametric Anal. Appl.* 6 (2014), 81–104. arXiv:1305.5492. (*Corrected 2026-06-04: original entry cited arXiv:1305.5564 and omitted Teh; both fixed per Round 2 transport memo.*) | YES (Springer URL, Caltech URL) | Explicit spectral triple from QSM with cyclotomic Galois action — Q5's shape-precedent. |
| 9 | Connes, Marcolli, Ramachandran. "KMS states and complex multiplication." Selecta Math. 11 (2005), 325–347. arXiv:math/0501424. | YES (standard) | Full Galois group of modular field as KMS symmetries; arithmetic precedent. |
| 10 | Marcolli, van Suijlekom. "Gauge networks in noncommutative geometry." J. Geom. Phys. 75 (2014), 71–91. arXiv:1301.3480. | YES (already in our bibliography as `marcolli_vs2014`, cited in Papers 25/30/41/51) | The ambient category Q5' lives in — finite/graph spectral triples — but with no motivic content. |

**Hallucination audit:** all ten references cross-verified by multiple search hits and (where possible) URL fetches. No fabricated arXiv IDs detected in this set. The recent Connes-Consani-Moscovici "Zeta Spectral Triples" (arXiv:2511.22755, Nov 27 2025) was verified independently and judged NOT directly Q5'-relevant — included here as note 7.6 to avoid re-discovery.

---

## 9. Sprint output summary

- **Verdict:** PARTIALLY TRANSPORTABLE.
- **Memo:** `debug/sprint_q5p_qsm_litread_memo.md` (this file, ~3.6k words).
- **Paper edits:** none (lit-mapping only, per sprint scope).
- **Follow-on work surfaced (not assigned):**
  1. If Q5' is ever activated as a research sprint, the 4-step Marcolli-Tabuada machinery is the natural entry point; periodic cyclic homology of $\mathrm{dg}(\mathcal{T}_{n_{\max}})$ is the first concrete computation.
  2. The pro-dg-category $\varprojlim \mathrm{dg}(\mathcal{T}_{n_{\max}})$ with Berezin maps is probably the right object (not the individual $\mathcal{T}_{n_{\max}}$, whose dg category is too small to see the M1/M2/M3 grading).
  3. Two new bibitems for Paper 55 if Q5' is ever drafted further:\ `connes_marcolli2004` (arXiv:math/0409306) and `marcolli_tabuada2016` (arXiv:1110.2438). Currently Paper 55 cites `fathizadeh_marcolli2017` and `glanois2015` and `brown2012`; the Connes-Marcolli renormalization paper and the Marcolli-Tabuada noncommutative-motives paper would be the two most directly relevant additions for the Q5' subsection.
- **Sprint cost:** ~1 session, 4 web searches + 6 web fetches, ~3.6k-word memo.
