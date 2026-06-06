# Sprint Q5'-U*-Action-Scoping (T5) — scoping memo for L1 follow-on (b): "verify cosmic-Galois U* acts on H_Levi"

**Date:** 2026-06-06 (T5 of the Q5'-HardParts-Round3 sprint)
**Sprint:** Q5' multi-year follow-on scoping for the L1 sub-sprint's named "verification that the cosmic-Galois $U^*$ of Connes–Marcolli 2007 acts on $\mathcal{H}_{\mathrm{Levi}}$ in the expected way" (Sprint Q5'-Levi-Synthesis §10.1 named follow-on b).
**Discipline:** structural scoping, no production code touched, no transcendentals introduced.

---

## 1. What is $U^*$? (Connes–Marcolli book 2008 Ch.~4, arXiv:math/0409306)

The **cosmic-Galois group** $U^*$ is a pro-affine algebraic group defined as follows. Let $\mathcal{M}_T$ denote the Tannakian category of mixed Tate motives over $\mathrm{Spec}(\mathbb{Z})$. The motivic fundamental group of this category (with respect to the de Rham fibre functor) is
$$
U^* := \mathrm{Aut}^{\otimes}(\omega_{\mathrm{dR}}|_{\mathcal{M}_T}).
$$
It is the **Tannakian dual** of $\mathcal{M}_T$. Concrete structure:
- $U^*$ is the semidirect product $\mathbb{G}_m \ltimes U$ where $U$ is the pro-unipotent radical (motivic Galois group of $\mathcal{M}_T(\mathbb{Z})_{\mathbb{Q}}$ relative to its half-Tate twist) and $\mathbb{G}_m$ is the Tate-weight grading.
- $U^*$ acts on the algebra $\mathcal{O}(U^*) = \mathcal{H}^{\mathrm{CK}}$ of Connes–Kreimer renormalization Hopf algebra, on the algebra of multi-zeta values (MZVs), and on every motivic period of $\mathcal{M}_T$ via the de Rham realisation.
- The published shape of $U^*$ is the **Levi-decomposition** $\mathbb{G}_m \ltimes (\mathrm{pro\text{-}unipotent}\;\mathrm{semi\text{-}direct})$ — exactly the Levi shape that GeoVac's L1 substrate
$$
U^*_{\mathrm{GeoVac},\;\mathrm{Levi}} = \mathbb{G}_a^{3 N(n_{\max})} \times SL_2
$$
matches structurally (Sprint Q5'-Levi-Synthesis v3.63.0).

## 2. What does "acts on $\mathcal{H}_{\mathrm{Levi}}$" mean? — three reasonable interpretations

The phrase "the cosmic-Galois $U^*$ acts on $\mathcal{H}_{\mathrm{Levi}}$" admits at least three structurally distinct readings:

### Interpretation A — Hopf-algebra coaction

$U^*$ coacts on $\mathcal{H}_{\mathrm{Levi}}$ via a Hopf-algebra map
$$
\rho: \mathcal{H}_{\mathrm{Levi}} \to \mathcal{H}_{\mathrm{Levi}} \otimes \mathcal{O}(U^*).
$$
This is the strongest interpretation. It would require $\mathcal{H}_{\mathrm{Levi}}$ to be an $\mathcal{O}(U^*)$-comodule algebra. The Tannakian theorem (Deligne 1990; Saavedra-Rivano 1972) would then identify the category of $\mathcal{O}(U^*)$-comodule algebras with a subcategory of $\mathcal{M}_T$ at finite cutoff.

### Interpretation B — Automorphism action

$U^*$ acts on $\mathcal{H}_{\mathrm{Levi}}$ by Hopf-algebra automorphisms:
$$
U^*(\mathbb{Q}) \to \mathrm{Aut}_{\mathrm{Hopf}}(\mathcal{H}_{\mathrm{Levi}}).
$$
This is weaker than Interpretation A but more directly verifiable. It requires identifying the Hopf-automorphism group of $\mathcal{H}_{\mathrm{Levi}}$ at each finite cutoff and embedding $U^*$-points into it.

### Interpretation C — Period-pairing action

GeoVac's spectral data (transcendental period values from M1/M2/M3 traces) pair with elements of $\mathcal{H}_{\mathrm{Levi}}$ via a "period map"
$$
\pi: \mathcal{H}_{\mathrm{Levi}} \to \mathbb{C},
$$
and $U^*$ acts on the image $\pi(\mathcal{H}_{\mathrm{Levi}}) \subset \mathbb{C}$ via motivic Galois conjugation. The Mellin engine's $\pi$ and $\zeta(\mathrm{odd})$ values would be the "GeoVac periods" on which $U^*$ acts via the standard $\mathcal{M}_T$ action on periods.

**Each of the three interpretations corresponds to a distinct sprint-scale verification target.** Identifying which is the load-bearing one for the Q5' program is part of the scoping problem.

## 3. Verification framework — what would each interpretation require?

### Interpretation A (coaction) verification framework

Required: explicit construction of the coaction map $\rho$ on the 18-dim Levi substrate $\mathcal{H}_{\mathrm{Levi}}^{(n_{\max}=2, j_{\max}=1/2)}$, with bit-exact verification of:
- $\rho$ is a Hopf-algebra map: $(\Delta \otimes \mathrm{id}) \circ \rho = (\mathrm{id} \otimes \mathrm{coproduct}_{U^*}) \circ \rho$.
- $\rho$ is counital: $(\varepsilon \otimes \mathrm{id}) \circ \rho = \varepsilon_{U^*} \cdot 1$.
- The image $\rho(\mathcal{H}_{\mathrm{Levi}})$ is a sub-Hopf-algebra of $\mathcal{H}_{\mathrm{Levi}} \otimes \mathcal{O}(U^*)$.

Sprint-scale prerequisite: explicit presentation of $\mathcal{O}(U^*)$ at finite truncation (analogous to Brown's $\mathcal{O}(U)$ truncation at weight $\le N$). This requires the Bar–Brown 2012 motivic-MZV machinery.

Estimated sprint-scale steps: (i) construct $\mathcal{O}(U^*)$ at weight $\le 6$ (~1 week); (ii) hypothesize coaction $\rho$ from M1/M2/M3 master-Mellin partition (~1 week); (iii) verify Hopf-axiom panel bit-exactly (~3 days).

### Interpretation B (automorphism) verification framework

Required: enumerate $\mathrm{Aut}_{\mathrm{Hopf}}(\mathcal{H}_{\mathrm{Levi}}^{(n_{\max})})$ explicitly. For the Levi-decomposition substrate $\mathbb{G}_a^{3 N(n_{\max})} \times SL_2$ at finite cutoff, the algebraic-group automorphism group is structurally
$$
\mathrm{Aut}_{\mathrm{Hopf}}(\mathcal{H}_{\mathrm{Levi}}) \supset \mathrm{Aut}(\mathbb{G}_a^{3 N(n_{\max})}) \times \mathrm{Aut}(SL_2) = GL_{3 N(n_{\max})}(\mathbb{Q}) \times \mathrm{PGL}_2(\mathbb{Q})
$$
at the connected component, plus a semi-direct factor encoding the Levi-trivial-action coupling.

The cosmic-Galois $U^*$ would embed as a sub-group preserving the M1/M2/M3 grading via Tate-weight grading on the abelian factor. Verification: identify which $GL_{3 N(n_{\max})}$ matrices preserve the $k$-slot grading and pair with the $SL_2$ factor through a Levi-decomposition-compatible action.

Sprint-scale prerequisite: $GL_{3 N(n_{\max})}$ automorphism enumeration at $n_{\max} \in \{2, 3\}$ (sprint-scale, ~1 week).

### Interpretation C (period-pairing) verification framework

Required: the period map $\pi: \mathcal{H}_{\mathrm{Levi}} \to \mathbb{C}$ defined by trace-evaluation against M1/M2/M3 spectral data. The image consists of GeoVac's empirically observed transcendentals (Paper 18 §III.7; M1 $\pi$-content, M2 Seeley-DeWitt $\zeta(2k)$, M3 vertex-parity Hurwitz $\zeta(\mathrm{odd})$, Catalan $G$).

Verification: $U^*$ acts on these transcendentals via the standard motivic-period action (Brown 2012, Fathizadeh–Marcolli 2016) and the period map is $U^*$-equivariant. Concrete test: compute the GeoVac period $\pi(\eta_{(n, l)}) \in \mathbb{Q}$ explicitly (Sprint Q5'-Prosystem v3.60.0; values $3, 3, 5, 15, 10, \dots$) and verify they sit in the **rational** sub-ring of MT periods — hence $U^*$ acts trivially on them at depth 0.

**This is the most sprint-scale-attackable interpretation.** The GeoVac periods at the JLO/CM cocycle-class level are already known bit-exactly (v3.61.0–v3.63.0). Most of them are rational (depth 0 in $\mathcal{M}_T$), with the master-Mellin partition's M1/M2/M3 carrying the depth-graded structure (Paper 55 §M3 main result: M3 is depth-graded cyclotomic mixed-Tate at level $\le 4$).

## 4. Sprint-scale sub-questions

Each of these is sprint-scale (1-2 weeks) and would substantively advance the L1 follow-on (b):

**SQ1 (Interpretation C, easiest):** Verify that the GeoVac periods $\chi_{(n, l)}, \eta_{(n, l)}$ from the v3.60.0 prosystem closure sit in the **rational** sub-ring of MT periods (depth 0). Bit-exact: $\chi_{(n, l)} \in \mathbb{Q}$ ✓ already verified; $\eta_{(n, l)} \in \mathbb{Q}$ ✓ already verified. This is a near-immediate sprint follow-on.

**SQ2 (Interpretation C, medium):** Verify that GeoVac's continuum Mellin lift $F(s)$ (Sprint T1, this sprint) has image in the MT period ring at integer $s$. From the T1 panel: $F(6) = -53\zeta(5)/3 + \pi^2/36 + 4\pi^6/945 + 11\zeta(3)/3 + 107\pi^4/540$ — every term is in MT level 0 (powers of $\pi^2 \cdot \mathbb{Q}$ and weight-3, weight-5 $\zeta$-values). Bit-exact verification of the period-ring containment: sprint-scale ~1 day.

**SQ3 (Interpretation B, medium):** Identify the automorphism group $\mathrm{Aut}_{\mathrm{Hopf}}(\mathcal{H}_{\mathrm{Levi}}^{(n_{\max}=2)})$ explicitly. Decompose as $GL_{15}(\mathbb{Q}) \times \mathrm{PGL}_2(\mathbb{Q})$ at the connected level, plus a semi-direct extension. Sprint-scale ~2 days. The cosmic-Galois $U^*$ embedding into this automorphism group is then a multi-year follow-on.

**SQ4 (Interpretation A, hardest):** Construct a candidate coaction $\rho$ from the M1/M2/M3 grading on $\mathcal{H}_{\mathrm{Levi}}$ to $\mathcal{H}_{\mathrm{Levi}} \otimes \mathcal{O}(U^*_{\le N})$ at truncated weight $\le 6$. Sprint-scale ~1-2 weeks. Verification of the Hopf-axiom panel: sprint-scale ~3 days.

**SQ5 (cross-interpretation):** Identify the period-action compatibility between the M1/M2/M3 partition (Paper 18 §III.7) and the standard MT depth grading. From the v3.62.0 M3 continuum sprint: M3 is depth-graded cyclotomic mixed-Tate at level $\le 4$ over $\mathbb{Z}[i, 1/2]$, with no $\pi$-power reduction at depth $\ge 1$. The M1/M3 ring asymmetry IS the operational content of $U^*$ acting differently on M1 (depth 0) vs M3 (depth $\ge 1$). Bit-exact compatibility verification at a 3-cell panel: sprint-scale ~3-5 days.

## 5. Multi-year obstructions

**MY1.** Explicit presentation of $U^*$ at all weights. The motivic-MZV literature (Brown 2012, Glanois 2015, Deligne 2010) gives $\mathcal{O}(U^*_{\le N})$ at truncated weight $\le N$ but the full Hopf-algebra structure at general weight requires the full motivic period machinery. Multi-year mathematical project.

**MY2.** Tannakian closure proof (L1 follow-on a). The "cosmic-Galois $U^*$ acts" question is the dual of the Tannakian closure question: $\mathcal{H}_{\mathrm{Levi}}$ as a comodule algebra is $\mathcal{M}_T$-realisation of GeoVac's substrate. Closure requires showing that finite-cutoff $\mathcal{H}_{\mathrm{Levi}}^{(n_{\max})}$ extends to a well-defined pro-Hopf-algebra in the inverse limit $n_{\max} \to \infty$, with $U^*$-coaction preserved. Multi-year.

**MY3.** Period conjecture compatibility. If the GeoVac periods (from the v3.61.0–v3.63.0 cocycle classes plus M1/M2/M3 spectral data) saturate the **conjectural** Brown period rationality (Brown 2012, Conjecture 3.7), then the $U^*$-action is uniquely determined by the motivic conjecture. If not, the action is multi-valued and the choice requires additional structural input from GeoVac's natural-geometry hierarchy. Multi-year.

## 6. Sprint-scale recommendation

**The cleanest sprint-scale next step is SQ2:** verify that the T1 continuum Mellin lift's integer-$s$ panel sits in the MT period ring bit-exactly. This is doable in ~1 day and would confirm Interpretation C (period-pairing) at the integer-$s$ panel level. Combined with SQ1's already-verified rationality of $\chi, \eta$, it would close Interpretation C at the cocycle-class level, leaving Interpretations A and B as the multi-year frontiers.

**SQ3 (automorphism enumeration) is the second priority** as it makes Interpretation B concrete. The result would be the cocoset $\mathrm{Aut}_{\mathrm{Hopf}}(\mathcal{H}_{\mathrm{Levi}}^{(n_{\max}=2)}) / U^*(\mathbb{Q})$ at finite cutoff — a finite-rank algebraic group that can be enumerated explicitly.

## 7. Hard prohibitions check (§13.5)

- No changes to natural geometry hierarchy.
- No fitted/empirical parameters.
- No deletion of negative results.
- Paper 2 combination-rule "conjectural" label unchanged.

## 8. Tag-transcendentals (`feedback_tag_transcendentals`)

All MT period content cited here lives in Paper 18 §III.7 master Mellin engine tiers:
- M1 ($\pi$-content, depth 0): GeoVac Hopf-base measure mechanism.
- M2 ($\pi^{2k}$ Seeley–DeWitt, depth 0): integer-shifted $\zeta(2k)$.
- M3 ($\zeta(\mathrm{odd})$ + Catalan + Hurwitz at $1/4$, $3/4$, $1/2$, $3/2$): cyclotomic mixed-Tate at level $\le 4$.

No new transcendentals introduced.

## 9. WH1 PROVEN unaffected

This scoping memo does not modify any propinquity result. It scopes a multi-year frontier of the Stage-2 Tannakian closure.

## 10. One-line verdict

**Three distinct interpretations of "cosmic-Galois $U^*$ acts on $\mathcal{H}_{\mathrm{Levi}}$" identified (coaction A / automorphism B / period-pairing C), with five sprint-scale sub-questions (SQ1–SQ5) named and three multi-year obstructions (MY1–MY3) named. SQ2 (verify T1's Mellin-lift panel sits in the MT period ring bit-exactly) is the cleanest sprint-scale next step.**
