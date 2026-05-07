# SM Gauge Content Forcing: Is U(1) × SU(2) × SU(3) the Maximal Wilson Content of GeoVac Sub-Manifolds?

**Date:** 2026-05-06
**Status:** Exploratory / structural
**Verdict:** **STRUCTURALLY PARTIAL, NUMEROLOGICALLY ATTRACTIVE.** The "natural" Wilson gauge groups on S^(2n−1) follow the sequence U(1), SU(2), SU(3), SU(4), … indexed by S^(2n−1) = SU(n)/SU(n−1). GeoVac produces S³ and S⁵ (and only those) as natural sub-manifolds, **truncating the sequence at SU(3)**. This is not a forcing argument that the SM gauge groups must appear; it is a forcing argument that **no other classical gauge group can appear naturally** in the current GeoVac framework. The two statements are very different. The SM gauge content match is real but compatible with several non-structural readings; the most honest framing is "GeoVac produces the SM gauge groups as the maximal Wilson content of its natural sub-manifolds, and no more, but the *selection* of sub-manifolds is mediated by which Hamiltonians (Coulomb, HO) the universe gives us, not by GeoVac-internal axioms."

---

## §1. The question and what it is and is not

**The question.** Are U(1), SU(2), SU(3) the *only* compact Lie groups that admit "natural" Wilson lattice gauge constructions on GeoVac-produced sub-manifolds (S³ from Paper 7, S⁵ from Paper 24)?

**Two readings of "natural":**

1. *Strict reading.* G admits a transitive action on the host manifold M (so M = G/H for some closed subgroup H), and the Wilson construction inherits this action structurally — link variables are tied to the transitive action, not bolted on.

2. *Generic reading.* Wilson lattice gauge theory works on **any** finite graph for **any** compact Lie group G (link variables in G, plaquettes from primitive cycles, action S_W = β·Σ_P (1 − (1/d_R)·Re Tr U_P), Haar measure normalizable for any compact G). Under the generic reading, G_2, F_4, E_6, E_7, E_8 all admit Wilson constructions on the Bargmann graph or the Coulomb graph. Nothing is forced.

This memo works under the strict reading. Under the generic reading the question is vacuous — anything goes — and the SM-content match is purely numerological.

**What the memo is.** A structural test of whether the SM gauge content emerges from a forcing argument of any kind, plus an honest caveat audit naming all readings under which the apparent forcing is coincidence.

**What the memo is not.** A claim that GeoVac predicts the Standard Model. The Sprint H1 / G3 four-gap analysis (`debug/sm_unified_gauge_synthesis_memo.md` §3) names four structural gaps: no autonomous Yukawa (G2), no chirality co-location (G3), no generation tripling (G1), and the cross-manifold obstruction (G4). The SM gauge group catalogue is at most a starting point; the SM is not in hand.

---

## §2. Mathematical fact: the sphere → gauge group sequence

The classical homogeneous-space identifications of odd-dimensional spheres:

| Sphere | Group identification | Hopf base | Natural gauge group |
|:------:|:---------------------|:----------|:--------------------|
| S¹ | U(1) | S⁰ (real) | U(1) |
| S³ | SU(2) | S² (complex Hopf) | SU(2), with U(1) sub-bundle |
| S⁵ | SU(3)/SU(2) | CP² (complex Hopf) | SU(3), with U(1) sub-bundle |
| S⁷ | Sp(2)/Sp(1), or Spin(7)/G₂, or SU(4)/SU(3) | S⁴ (quaternionic Hopf) | SU(4), Sp(2), and Spin(7) all act |
| S⁹ | SU(5)/SU(4) | CP⁴ | SU(5) |
| S^(2n−1) | SU(n)/SU(n−1) | CP^(n−1) | SU(n) |

**Two ways to read this.**

The complex Hopf sequence (Hopf base = CP^(n−1)) gives:
- S^(2n−1) admits SU(n) transitively
- The natural gauge group for S^(2n−1) under the complex Hopf reading is SU(n)
- The U(1) sub-bundle (the Hopf circle) gives the abelian sub-construction

So the sequence of natural gauge groups indexed by complex-Hopf spheres is: **U(1), SU(2), SU(3), SU(4), SU(5), …** — and the SM gauge content is the **truncation at the first three terms**.

The quaternionic Hopf and octonionic Hopf give:
- S⁷ → S⁴ admits Sp(1) ≅ SU(2) bundle (the quaternionic Hopf)
- S^15 → S⁸ admits Spin(7) bundle (the octonionic Hopf)

These are the only Hopf fibrations (Adams 1960; only finite-dimensional normed division algebras over ℝ are ℝ, ℂ, ℍ, 𝕆, giving exactly four Hopf fibrations).

**G₂ on S⁶.** S⁶ = G₂/SU(3) is the unique sphere admitting a transitive G₂ action (Borel 1949). G₂ also acts transitively on S⁷ via Spin(7) ⊃ G₂. S⁶ is even-dimensional (not odd) and so is not produced by the complex Hopf sequence; G₂ is exceptional.

---

## §3. What GeoVac selects from this sequence

GeoVac produces sub-manifolds via two mechanisms:

| Mechanism | Sub-manifold | Source |
|:----------|:-------------|:-------|
| Fock projection of the −Z/r Hamiltonian (Paper 7) | S³ | Bertrand-Coulomb maximal symmetry SO(4) |
| Bargmann-Segal of the 3D HO Hamiltonian (Paper 24) | S⁵ | HO maximal symmetry SU(3) on (N,0) Hardy sector |

These are **the only sub-manifolds GeoVac produces from first-principles physics**. Other manifolds (S⁷, S⁹, …, ℝH³ via Wick rotation, CP^n, etc.) have been explicitly tested and either (a) are not natural images of any GeoVac packing (the higher complex Hopf S⁷/S⁹ would require physical 4D/5D harmonic-type Hamiltonians the framework does not produce), or (b) are produced as quotients but not as gauge-host structures (the m_l-quotient → 12-sector quasi-CP², see Sprint 5 Track S5 negative on CP² spectrum match), or (c) were closed by negative result (Wick-rotated H³ on hydrogen via Bander-Itzykson 1966 — Track RH-B, no discrete Γ ⊂ SO(3,1) selectable from GeoVac invariants).

**Why two and only two natural sub-manifolds?** The deeper fact is *Bertrand's theorem* (1873): the only central potentials with all bound orbits closed are the Coulomb 1/r and the harmonic r². Their maximal symmetry groups are SO(4) (Coulomb) and SU(3) (HO). GeoVac's natural sub-manifolds are the symmetry-orbit spaces of these two maximally symmetric Hamiltonians. The "two and only two" comes from a classical result about classical mechanics, not from a GeoVac-internal axiom.

**Mapping to the sphere-gauge sequence.**

- S³ is at position n=2 in the complex Hopf sequence: natural gauge group SU(2), sub-bundle U(1).
- S⁵ is at position n=3: natural gauge group SU(3), sub-bundle U(1).
- S¹ (n=1) is the Hopf fiber — not produced as a top-level sub-manifold but appears as the fiber of the S³ → S² Hopf bundle (Paper 25 uses this).

So GeoVac's natural Wilson gauge content is the truncation of the complex-Hopf gauge sequence at n ≤ 3. The SM is U(1) × SU(2) × SU(3) — exactly this truncation.

---

## §4. Test cases for "more groups can fit"

I worked through five candidate gauge groups not currently in GeoVac and checked whether any of them admit a *strict-reading* natural Wilson construction on S³ or S⁵.

### 4.1 SU(4) on S⁵ — NO

S⁵ = SU(3)/SU(2). SU(4) does not act transitively on S⁵. The natural SU(4) home is S⁷ = SU(4)/SU(3) — and GeoVac does not produce S⁷.

A Wilson SU(4) construction can mechanically be put on the S⁵ Bargmann graph (link variables in fundamental SU(4), plaquettes, Wilson action — same machinery as Sprint ST-SU3 for SU(3)), but there is no transitive SU(4) action on S⁵ to give the construction structural meaning. Under the strict reading: not natural.

### 4.2 Sp(2)/USp(4) on S³ or S⁵ — NO

Sp(2) acts transitively on S⁷ = Sp(2)/Sp(1) (the quaternionic Hopf). On S³ or S⁵, Sp(2) admits no transitive action. Under the strict reading: not natural for either GeoVac sub-manifold.

(Sp(1) ≅ SU(2) acts on S³ trivially through the SU(2) identification — but this is just SU(2) again, not a new gauge group.)

### 4.3 G₂ on S³ or S⁵ — NO

G₂ acts transitively on S⁶ = G₂/SU(3) and on S⁷ via Spin(7). It does not act transitively on S³ or S⁵. GeoVac does not produce S⁶. Under the strict reading: not natural.

### 4.4 SO(5) on S⁴ or S³ — NO

SO(5) acts transitively on S⁴ = SO(5)/SO(4). GeoVac does not produce S⁴ as a natural sub-manifold (it appears as the quaternionic Hopf base S⁷ → S⁴, but S⁷ is not GeoVac). On S³, SO(4) acts transitively, not SO(5). On S⁵, SO(6) acts transitively, not SO(5).

The SO(5) candidate is interesting because Sp(2) ≅ Spin(5) and Sp(2) does act on S⁷, but that's the quaternionic Hopf again, not a GeoVac sub-manifold.

### 4.5 F₄, E₆, E₇, E₈ — NO

The exceptional Lie groups have no transitive sphere actions in low dimensions other than G₂ on S⁶ (already discussed). They are not natural for any GeoVac sub-manifold.

### 4.6 Summary table

| Group | Strict-natural sub-manifold | Is it a GeoVac sub-manifold? |
|:-----:|:---------------------------|:----------------------------:|
| U(1) | S¹, also CP^(n−1) Hopf base of any S^(2n−1) | YES (Hopf base of S³ in Paper 25) |
| SU(2) | S³ = SU(2) | YES (Paper 30) |
| SU(3) | S⁵ = SU(3)/SU(2) | YES (Sprint ST-SU3) |
| SU(4) | S⁷ = SU(4)/SU(3) | NO |
| SU(5)+ | S⁹+ | NO |
| Sp(n) | S^(4n−1) (quaternionic Hopf) | NO (no S⁷ in GeoVac) |
| G₂ | S⁶, S⁷ (via Spin(7)) | NO |
| SO(n) | S^(n−1) standard | partial — SO(4) acts on S³, SO(6) on S⁵, but these are not Wilson-natural in the gauge-group sense; their Wilson constructions are the SU(2)×SU(2) and SU(4) (resp.) double covers, which are already in the sequence |
| F₄, E₆, E₇, E₈ | none in low dim | NO |

**Conclusion of §4:** Under the strict reading (transitive G-action on the host manifold), **the only natural gauge groups for GeoVac sub-manifolds are exactly U(1), SU(2), SU(3)**, plus their sub-bundles. Higher SU(n>3), Sp(n), G₂, exceptional groups are all ruled out by the absence of S⁷ and S^(2n−1>5) as GeoVac sub-manifolds.

---

## §5. The four-layer Coulomb/HO asymmetry as forcing argument

Sprint ST-SU3 §5.4 stated the asymmetry in four layers:

| Layer | S³ Coulomb | S⁵ Bargmann |
|:-----:|:-----------|:-------------|
| (i) Spectrum-computing role of L₀ | YES (κ·(D−A) = Rydberg) | NO (HO spectrum in diagonal) |
| (ii) Calibration π content | YES (Fock projection) | NO (linear-affine projection) |
| (iii) Wilson gauge with natural matter coupling | YES (SU(2) Wilson + Dirac) | NO (CG intertwiner needed) |
| (iv) Universal Wilson-Hodge vocabulary | YES | YES |

Layer (iii) is the gauge-content layer. The asymmetry there says:
- On S³, the natural gauge group SU(2) couples to natural matter (spin-1/2 fermions in the same fundamental rep — the "manifold IS the group" coincidence makes matter and gauge live in the same Hilbert space).
- On S⁵, the natural gauge group SU(3) couples through CG intertwiners to (N,0) shells of growing dimension — Wilson SU(3) is well-defined as a *pure-gauge* theory but matter coupling rediscovers Sprint 5 Track S5's CG obstruction.

**Reframing as a forcing argument.** What the four-layer asymmetry actually says is:

> The natural gauge content of S³ is "complete" — a fully matter-coupled non-abelian gauge theory exists. The natural gauge content of S⁵ is "partial" — the gauge sector exists but matter coupling requires non-Wilson structure.

This is *not* a forcing argument that the SM gauge groups are inevitable. It is a structural observation that S³ is special in admitting a complete gauge-with-matter theory while S⁵ is special in admitting a gauge-only theory. The fact that the SM gauge groups happen to match the natural gauge content of these two sub-manifolds is consistent with this asymmetry — but the asymmetry by itself does not predict that the SM gauge content must be exactly U(1) × SU(2) × SU(3), nor that the SM matter content must follow.

**What the asymmetry DOES force:**
- If GeoVac is the right framework, then the matter sector should look "complete" on S³ (Dirac fermions living in fundamental rep of SU(2)) and "incomplete" on S⁵ (no natural matter coupling — color confinement?).
- The strong coupling sector being on S⁵ where matter coupling is structurally awkward (matter reps grow with N, requiring CG intertwiners between distinct irreps) is consistent with the empirical observation that QCD matter (quarks) cannot be observed in isolation.
- This is a **suggestive analogy**, not a derivation.

---

## §6. Honest caveat audit

Six readings under which the apparent SM-gauge-content forcing is coincidence:

### 6.1 Selection bias from the small finite set of "natural" structures

The classical compact Lie groups in low dimensions form a small finite list: U(1), SU(2), SU(3), SU(4), SU(5), SO(3), SO(4), SO(5), SO(6), Sp(1), Sp(2), G₂. A small enumeration of "natural" structures has high prior probability of overlapping with the SM by chance. The match U(1) × SU(2) × SU(3) is suggestive but not statistically improbable given the small population of low-dimensional gauge groups.

### 6.2 Bertrand's theorem is classical mechanics, not GeoVac-internal

The reason GeoVac produces only S³ and S⁵ as natural sub-manifolds is *Bertrand's theorem* — the Coulomb 1/r and harmonic r² are the only central potentials with all bound orbits closed. This is a 19th-century theorem about classical mechanics, not a derivation from GeoVac axioms. GeoVac inherits the "two natural sub-manifolds" feature from Bertrand, not from packing.

### 6.3 The complex-Hopf-sequence framing is one of several

S^(2n−1) admits multiple natural group actions. S⁷ admits SU(4), Sp(2), Spin(7) all transitively. The complex-Hopf reading SU(n) on S^(2n−1) is one organizing principle; the quaternionic and octonionic Hopf give different sequences. The match between SM and complex-Hopf-truncated-at-n=3 is a choice of organizing principle, not the unique one.

### 6.4 The match doesn't extend to representation content

The SM has specific hypercharge assignments (Y = −1, +2/3, −1/3, …), specific chirality assignments (left-handed doublets, right-handed singlets), three generations, the CKM matrix, Yukawa couplings, the Higgs sector. GeoVac's three Wilson constructions have none of this representation content built in. Sprint 4H Track SM-B closed the naive "shell → generation" map in the negative. Sprint G3 closed the "GeoVac chirality = electroweak chirality" map in the negative. The SM gauge GROUPS match; the SM gauge ASSIGNMENTS do not.

### 6.5 The cross-manifold obstruction is the dominant gap

Sprint H1 §3.4 / Sprint G4 (just yesterday): the SM is one spectral triple with one almost-commutative algebra; GeoVac is at present a direct sum of three Wilson constructions on three sub-manifolds. The four-way S³ / S⁵ tensor product would require an NCG-framework extension that has no published prescription. Even if the gauge groups match, the way the SM puts them together (single triple, single Hilbert space, internal symmetry decomposition) is structurally different from how GeoVac has them (three separate triples on different sub-manifolds).

### 6.6 The Wilson construction is generic — strictly nothing is forced

Under the *generic* reading of "natural" (any compact G admits a Wilson lattice gauge theory on any finite graph), the SM gauge content match is selection bias on the strict reading. The strict reading itself is a choice. A reasonable physicist could argue that putting F₄ Wilson on the Bargmann graph is "just as natural" — it works mathematically, has the universal 1/(4·N_c) coefficient, etc. — and the strict-reading "naturalness" is being defined to match the SM after the fact.

This is the most cutting caveat. The strict-reading natural-gauge-group condition is genuinely satisfied (transitive action on host manifold) but the choice to use the strict reading is, itself, a choice that aligns the framework with the SM.

---

## §7. Net verdict

**STRUCTURAL component (real):**

1. The complex-Hopf sequence S^(2n−1) → SU(n) is a published mathematical fact.
2. GeoVac produces only S³ and S⁵ as natural sub-manifolds, by mechanisms (Fock projection, Bargmann-Segal) that are GeoVac-internal up to Bertrand's theorem.
3. Under a "strict-natural" reading (transitive G-action on host), the only candidate gauge groups are U(1), SU(2), SU(3), with no candidate for SU(4)+, Sp(n), G₂, F₄, etc.
4. The four-layer Coulomb/HO asymmetry is structurally specific and constrains the matter-coupling content of each sector.
5. **No higher gauge group can be added to GeoVac without first adding a new sub-manifold**, and no new sub-manifold is naturally produced by the framework's two existing mechanisms.

**NUMEROLOGICAL component (acknowledged):**

1. The match to SM gauge groups is real but compatible with selection bias from small population of low-dim Lie groups.
2. Bertrand's theorem is classical, not GeoVac-internal.
3. Three of four SM features (hypercharge, chirality, generations) are structurally absent from GeoVac (all confirmed by negative results in Sprint 4H, Sprint G3, Sprint H1).
4. The cross-manifold obstruction is unsolved (G4b open, blocked at NCG-framework level).
5. The "strict-natural" reading is a choice.

**The honest framing:**

> GeoVac produces the SM gauge groups U(1) × SU(2) × SU(3) as the maximal natural Wilson content of its sub-manifolds, **and no more**. Any higher-rank or exceptional gauge group would require a sub-manifold (S⁷, S⁶, etc.) that GeoVac does not produce by either of its construction mechanisms. This is a *truncation forcing argument*: the framework cannot host more than these three gauge groups under the strict-natural reading. It is **not** a forcing argument that these three must appear with their SM rep content, hypercharges, generations, or chirality assignments — those four structural features are documented gaps (G1, G2, G3, G4 in Sprint H1's analysis).

This is interesting — it answers a different question than "does GeoVac predict the SM" — it answers "*does GeoVac admit any extension beyond the SM gauge content under its current sub-manifold structure*." The answer is **no**, under the strict reading. The framework is gauge-content-saturated by SM groups.

The numerological match itself — that the truncation lands exactly at the SM — is not a derivation. But the truncation existing at all (no SU(4), no G₂, no F₄ admissible) is a structural feature, not a coincidence.

**For the chat with the PI:** The right thing to say is something like "the framework can't host more gauge content than the SM has, but it also can't justify why the SM has exactly this content rather than a subset (e.g., just SU(3))." The forcing is on the upper bound, not on the lower bound.

---

## §8. What would close this question structurally

Three things would upgrade this from "compatible with SM" to "forces SM gauge content":

1. **Generation tripling structure.** A GeoVac mechanism that selects three copies of H_GV would address G1. None has been found in nine months of investigation.

2. **Closed forcing argument that S⁷ cannot exist in GeoVac.** Currently the absence of S⁷ is empirical (no mechanism has produced it). A structural theorem ("any GeoVac packing produces only spheres of dim 3 and 5") would upgrade the truncation argument from "currently no candidate" to "no candidate exists." Bertrand's theorem is the natural starting point but doesn't directly give this.

3. **Co-located SM gauge content on a single triple.** Sprint G4a (Connes SM on T_S³ alone with 𝒜_F = ℂ ⊕ ℍ ⊕ M₃(ℂ)) would address whether GeoVac can host the SM gauge content as inner automorphisms of one triple, rather than as three separate Wilson constructions. This is the agreed next sprint after WH1 PROVEN closes; predicted positive-thin per Sprint H1 / G3.

None of these is in hand. The SM gauge content match is structurally suggestive and numerologically attractive but is not a derivation.

---

## §9. Files referenced

- `papers/observations/paper_30_su2_wilson.tex` (SU(2) on S³)
- `papers/synthesis/paper_25_hopf_gauge_structure.tex` (U(1) on S³)
- `papers/core/paper_24_bargmann_segal.tex` (S⁵ Bargmann lattice; HO rigidity; Coulomb/HO asymmetry §V)
- `geovac/su2_wilson_gauge.py`, `geovac/su3_wilson_s5.py`
- `debug/st_su3_wilson_memo.md` (Sprint ST-SU3 — SU(3) on S⁵)
- `debug/s5_gauge_structure_memo.md` (Sprint 5 Track S5 — SU(3) negative on (N,0) tower)
- `debug/sm_unified_gauge_synthesis_memo.md` (May 4 synthesis; four-gap analysis)

External references:
- Adams, J. F. (1960). "On the non-existence of elements of Hopf invariant one." *Annals of Mathematics* **72** — only four normed division algebras / four Hopf fibrations.
- Borel, A. (1949). "Some remarks about Lie groups transitive on spheres and tori." *Bulletin AMS* **55** — classification of compact Lie groups acting transitively on spheres.
- Bertrand, J. (1873). "Théorème relatif au mouvement d'un point attiré vers un centre fixe." *Comptes Rendus* **77** — only Coulomb and harmonic central potentials have all bound orbits closed.
- Marcolli, M. & van Suijlekom, W. (2014). "Gauge networks in noncommutative geometry." J. Geom. Phys. **75**, arXiv:1301.3480.
