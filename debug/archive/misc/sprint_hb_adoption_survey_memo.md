# Hain–Brown adoption survey — what GeoVac can take, what GeoVac can give

**Author / harness:** Survey sub-agent, 2026-06-06.
**Status:** Read-only survey. No paper modifications, no code modifications.
**Scope:** Following the v3.78.0 closure-of-day memo (`debug/strategic_synthesis_2026_06_06_memo.md` §Addendum) which named Hain–Brown's relative completion of $SL_2(\mathbb{Z})$ as the structural-shape match for GeoVac's $U^* = \mathbb{G}_a^\infty \rtimes SL_2$, this memo identifies (i) concrete tools / theorems / techniques from the Hain–Brown program that GeoVac could productively adopt; (ii) reciprocal contributions GeoVac brings; (iii) an honest comparison; (iv) ranked recommendations.

**Bibliographic note up front (audit catch).** The memory file [[hain_brown_identification]] and the v3.78.0 strategic synthesis memo both cite **"Bouillon 2023 arXiv:2303.17534"** as the recent physics-side application of relative-completion machinery. Cross-checking against arXiv directly, the author of that paper is **Matija Tapušković**, *not* Bouillon. The paper title and arXiv ID are correct; only the author attribution is wrong. This is a clean fix (one bibitem in any paper that cites it). I have not corrected the memory file or any paper in this sprint per the read-only constraint, but flag it here for next sprint's bib pass.

---

## §1. Plain-language summary of Hain–Brown machinery

**The setup.** $\mathcal{M}_{1,1}$ is the moduli space of elliptic curves; its orbifold fundamental group is $SL_2(\mathbb{Z})$. Classical theory associates to a group $\Gamma$ and a reductive group $R$ with a Zariski-dense map $\rho: \Gamma \to R(\mathbb{Q})$ its **relative completion** $\mathcal{G}^{\mathrm{rel}}_{\Gamma}$, fitting into an exact sequence

$$
1 \to \mathcal{U} \to \mathcal{G}^{\mathrm{rel}} \to R \to 1,
$$

with $\mathcal{U}$ pro-unipotent. For $\Gamma = SL_2(\mathbb{Z})$ and $R = SL_2/\mathbb{Q}$ (with $\rho$ the inclusion), this is the relative completion at the heart of Hain–Brown. The pro-unipotent radical $\mathcal{U}$ is the "modular cousin" of Drinfeld's pro-unipotent completion of $\pi_1(\mathbb{P}^1 \setminus \{0,1,\infty\})$ in the mixed-Tate setting.

**Hain's theorem (1403.6443, with antecedents in his earlier work).** The coordinate ring $\mathcal{O}(\mathcal{G}^{\mathrm{rel}})$ carries a canonical **mixed Hodge structure** (MHS) compatible with the Hopf-algebra structure (product, coproduct, antipode). The Lie algebra $\mathfrak{u}$ of the pro-unipotent radical is freely topologically generated (non-canonically) by

$$
\bigoplus_{m \ge 0} H^1(SL_2(\mathbb{Z}), S^m H)^* \otimes S^m H,
$$

where $H = \mathbb{Q}^2$ is the standard rep of $SL_2$. By Eichler–Shimura, $H^1(SL_2(\mathbb{Z}), S^m H)$ is the space of modular forms of weight $m+2$ — so the **generators of $\mathfrak{u}$ are labeled by modular forms**. Eisenstein series $E_4, E_6, \ldots$ give Eisenstein generators (often denoted $e_{2n}$ or $\epsilon_{2n}$); cusp forms give cuspidal generators ($\check{M}_f$, where $f$ ranges over a basis of cusp forms including $\Delta$ of weight $12$, etc.).

**Pollack's quadratic relations and the Eisenstein quotient.** The "Eisenstein quotient" $\mathfrak{u}^{\mathrm{eis}}$ is obtained by setting the cuspidal generators to zero. Pollack discovered explicit quadratic relations among Eisenstein generators in $\mathfrak{u}^{\mathrm{eis}}$ (later extended in depth by Baumard–Schneps via Écalle moulds, and pushed to depth 3 by Pollack — see arXiv:1504.04737). These relations are visible at low weight and tie back to period polynomials of cusp forms.

**Brown's machinery (1407.5167).** Brown's *multiple modular values* (MMVs) are the periods of $\mathcal{O}(\mathcal{G}^{\mathrm{rel}})$ computed as **regularised iterated Eisenstein integrals** on the upper half plane. They generalise Manin's iterated Shimura integrals. The periods are explicitly:
- Powers of $\pi$ (from the Tate factor / Hopf-base);
- Classical zetas $\zeta(2k+1)$ at odd arguments (entering through Eisenstein period polynomials);
- Special values of $L$-functions of cusp forms (e.g., $L(\Delta, k)$ at integer arguments inside $[1, 11]$);
- Numerical mixtures of the above via the iterated-integral coaction.

**Hain–Matsumoto (1512.03975).** Bundles the above into a $\mathbb{Q}$-linear Tannakian category $\mathrm{MEM}_1$ of **universal mixed elliptic motives** over $\mathcal{M}_{1,1}$. Objects are mixed Tate motives unramified over $\mathbb{Z}$ equipped with a compatible $SL_2(\mathbb{Z})$-action. $\mathrm{MEM}_1$ contains $\mathrm{MTM}(\mathbb{Z})$ (Brown's mixed Tate motives) as a sub-Tannakian category, and its Tannakian fundamental group is exactly the relative completion described above.

**Tapušković 2023 (2303.17534, often miscited as "Bouillon").** First physics-side input into the Hain–Brown program. The equal-mass sunrise Feynman integral is realised as a period of the relative completion of $\Gamma_1(6)$ (not $SL_2(\mathbb{Z})$ — a different modular group). The cosmic Galois group acts on the motivic lift of the sunrise; conjugates are expressed in terms of motivic lifts of Feynman integrals associated to related (edge-subdivided) graphs. The methodology is: lift the physics integral to a motive in MEM$_{\Gamma_1(6)}$-like category, compute the cosmic-Galois coaction.

**2024–2026 frontier.** (a) Kleinschmidt et al., *Towards Motivic Coactions at Genus One from Zeta Generators* (arXiv:2508.02800, JHEP 05 (2026) 105) — proposes explicit coaction formulae for iterated Eisenstein integrals at genus one via "zeta generators" $\sigma_w$ with arithmetic ($z_w$) and geometric (Tsunogai derivation) parts. Closes a depth-1 computational gap left open by Brown 2014. Honest scope: still a *proposal*; no first-principles de Rham–periods identification yet. (b) Baune–Broedel–Moeckli, *Single-valued elliptic polylogarithms* (arXiv:2511.15240) — extends Brown's single-valued construction to the once-punctured elliptic curve. (c) Hain, *Hecke Actions on Loops and Periods of Iterated Shimura Integrals* (arXiv:2303.00143) — extends MHS-compatibility to Hecke operators; Annales ENS, revisions through June 2025. (d) Hain, *Algebraic de Rham theory for relative completion of $SL_2(\mathbb{Z})$* (arXiv:1804.06977) — gives the algebraic-de-Rham realisation needed for explicit period computations.

**Aggregate.** The Hain–Brown program is the most concrete realisation in mathematics of a pro-algebraic group of the shape (reductive $SL_2$ ⋉ pro-unipotent) with a canonical mixed Hodge structure and a working iterated-integral period theory. Generators are labeled by modular forms, relations are partially understood (Pollack + Baumard–Schneps), and the periods are computable as iterated Eisenstein integrals.

---

## §2. Concrete adoption candidates, ranked

I rank by **(information gain) / (assembly cost)**. Each candidate cites the structural fit, the GeoVac home (paper + section + file), the cost class (sprint / multi-month / multi-year), and what would change in the corpus if adopted.

### A1 — Hodge / Hodge–de-Rham realisation as the natural target

**Cost:** sprint-scale (1–2 weeks, framing pass).
**What it adopts:** the explicit decision to position GeoVac's $U^*$ against the *Hodge realisation* of Hain–Brown, not étale or $\ell$-adic.
**GeoVac home:** Paper 55 §1 (framing of period rings); Paper 56 §sec:open_g4 + §sec:open_profinite; memory [[hain_brown_identification]].
**Information gain:** Currently the closure-of-day memo says "Hodge / Hodge-de-Rham over étale / $\ell$-adic" but Paper 56 §sec:open_g4 doesn't make this commitment in the writeup. Locking it in does two things: (i) matches the substrate Hain–Brown built; (ii) closes the realisation question without committing to a $p$-adic comparison that GeoVac has no input for. The four named Connes axioms verified at finite cutoff (Paper 32 §IV) are real-Hodge structure prerequisites, not étale ones.
**Risk:** none — this is editorial clarity, not new math.

### A2 — Iterated Eisenstein integrals as a coaction language for GeoVac periods

**Cost:** sprint-scale (2–3 weeks, computational test).
**What it adopts:** Brown's regularised iterated Eisenstein integrals as the candidate period-theoretic basis against which to PSLQ GeoVac's M3 vertex-parity outputs.
**GeoVac home:** Paper 55 §6 (joint engagement / period catalogue); a new computational driver in `debug/`.
**Information gain:** This is the operationalisation of **Test A from the strategic memo**. Hain–Brown's basis includes the classical odd zetas $\zeta(2k+1)$ and $L$-values of cusp forms (specifically $L(\Delta, k)$). GeoVac's M3 sector produces $\zeta(2k+1)$, $\beta(2k)$ ($\beta = L(\cdot, \chi_{-4})$), and the depth-$k$ tower starting from $S_{\min}$. The 2026 Kleinschmidt et al. paper (arXiv:2508.02800) gives explicit coaction formulae for iterated Eisenstein integrals — which is the missing computational tool we'd need to PSLQ on the modular side. **The cleanest test:** compute a few M3 outputs at depth 2 (joint Mellin of two Hurwitz-at-quarter-integer products) and PSLQ against $\{ \zeta(3), \zeta(5), \beta(2), \beta(4), L(\Delta, 3), L(\Delta, 5), \pi^k \}$. If any cusp-form-$L$-values appear, **the identification is empirically supported and GeoVac is structurally in MEM rather than MTM**. If only MZV-cyclotomic stuff appears, fall back to $\mathcal{G}_4$ as comparison target.
**Risk:** the depth-2 Mellin transform is structurally the same object as the NA-1 depth-2 Mellin test (Reading A vs B) already named in the strategic memo §6. The two tests can share infrastructure — basically free.

### A3 — Eichler–Shimura labeling of GeoVac's primitive generators

**Cost:** multi-month (2–4 months structural assembly).
**What it adopts:** the labeling principle that Lie algebra generators sit naturally in modular-form spaces $H^1(SL_2(\mathbb{Z}), S^m H)$, with Eisenstein and cuspidal classes giving structurally distinct kinds of generators.
**GeoVac home:** Paper 56 §sec:open_na1 — Reading A vs B (abelianization vs shuffle enrichment); the explicit generators in `geovac/tannakian.py`.
**Information gain:** Currently GeoVac's primitive generators are labeled by (sector, Mellin slot) pairs at each $n_{\max}$. Eichler–Shimura provides a way to **classify these by modular-form weight** if the M3 sector is identified with a Hodge realisation of a modular extension. **Specifically:** if the depth-2 PSLQ test (A2) reveals $L$-values of $\Delta$, then GeoVac has a Pollack-quadratic-relation-shaped structure on its $SL_2$ piece. This would partially close the converse-equality direction at the level of the cuspidal-vs-Eisenstein split.
**Risk:** The labeling principle requires (i) the depth-2 test passing positively (A2), (ii) building the Eichler–Shimura translation, (iii) checking whether GeoVac's panel generators sit in the right $H^1$ classes. Each step is a sub-sprint. Worth it if A2 passes.

### A4 — Pollack quadratic relations as candidate relations on GeoVac's $\mathfrak{u}^*$

**Cost:** sprint-scale (1–2 weeks, falsifiability test).
**What it adopts:** explicit quadratic relations from Pollack (arXiv:1504.04737) and Baumard–Schneps (arXiv:1510.05549) among Eisenstein generators in the Eisenstein quotient $\mathfrak{u}^{\mathrm{eis}}$.
**GeoVac home:** Paper 56 §sec:tc1e — TC-1e and TC-1f verified the abelian $\mathbb{G}_a$ generators commute trivially. A natural strong-form question: do Pollack's quadratic relations *fail* on GeoVac's substrate (consistent with current abelian reading, Reading A), or do they appear modulo higher-cutoff corrections (consistent with Reading B, where the substrate should be enriched)?
**Information gain:** This is a **falsifiability test for Reading A vs B** that is *cleaner than depth-2 Mellin* because it's a finite-cutoff algebraic identity. Take Pollack's lowest-weight quadratic relation, write it in terms of GeoVac's Lie-algebra generators (after the modular-form labeling of A3), and check whether it is bit-exactly zero, bit-exactly nonzero, or something in between (corrections at higher cutoff).
**Risk:** depends on A3 being done first to know the labeling. If A3 hasn't happened, this is multi-month.

### A5 — Algebraic de Rham apparatus (Hain 1804.06977)

**Cost:** multi-month (3–5 months to import).
**What it adopts:** Hain's explicit algebraic-de-Rham description of $\mathcal{O}(\mathcal{G}^{\mathrm{rel}})$ — coordinate-ring presentation amenable to comparison with concrete Hopf algebras.
**GeoVac home:** Paper 56 explicit construction of $\Phi: U^*_{\mathrm{Levi}} \to \mathrm{Aut}^\otimes(\omega)$.
**Information gain:** GeoVac's substrate $\mathcal{H}_{GV}(n_{\max}) = \mathrm{Sym}_\mathbb{Q}(V_{n_{\max}})$ is already a concrete Hopf algebra; Hain's algebraic-de-Rham construction gives a concrete Hopf algebra on the modular side. A morphism of Hopf algebras between the two would be a **structural identification at the algebra level**, sharper than the period-ring-level identification in §sec:open_g4. Multi-month because Hain's construction is technical and the morphism has to commute with both Hopf structures.
**Risk:** moderate. The "abelian vs free non-abelian" obstruction (NA-1) is real and a morphism in either direction has to either abelianize Hain or shuffle-enrich GeoVac.

### A6 — Hecke action compatibility (Hain 2303.00143)

**Cost:** multi-year (no sprint-scale entry).
**What it adopts:** Hain's Hecke-equivariant MHS on iterated Shimura integrals; the algebra generated by Hecke operators is non-commutative on the periods.
**GeoVac home:** None currently. Would open a new front.
**Information gain:** GeoVac has no Hecke-type structure at present. A natural place it *could* appear is in the inner-factor sector (Paper 18 §IV.6 inner-factor input data tier) where mass eigenvalues / mixing angles arise. But there's no concrete handle today.
**Risk:** speculative. Park for now.

### A7 — Tapušković 2023 edge-subdivision as a model for physics-side input

**Cost:** sprint-scale (2 weeks, methodology pass).
**What it adopts:** Tapušković's edge-subdivision construction: for a Feynman graph $G$ and a subdivided graph $G'$, the motivic lifts are periods of *the same* relative completion. This is a structural device showing **how a physics integral relates to multiple motivic lifts** under coaction.
**GeoVac home:** Paper 18 §III (the projection chain framework); could appear as a new §III sub-projection (edge-subdivision-as-projection).
**Information gain:** This is the *only* published precedent for a physics-side input into the relative-completion program. GeoVac's discrete-graph nature means *every* matrix element is a "graph integral" in some sense, with the master Mellin engine's three slots playing roles analogous to graph operations. Studying Tapušković's edge-subdivision construction explicitly may reveal whether GeoVac's $D^k$ powers in the master Mellin engine ($k = 0, 1, 2$) correspond to graph operations in the modular setting.
**Risk:** Tapušković works with $\Gamma_1(6)$, not $SL_2(\mathbb{Z})$. The right comparison may be at a different level. But the methodology is generic.

### A8 — Kleinschmidt 2026 zeta-generator coaction as computational engine

**Cost:** sprint-scale (1–2 weeks).
**What it adopts:** explicit coaction formulae (Eq. 64 of arXiv:2508.02800) for iterated Eisenstein integrals at genus one, with zeta generators $\sigma_w$ having concrete arithmetic + geometric (Tsunogai) decompositions.
**GeoVac home:** A new computational driver in `debug/`; supports A2 (PSLQ test) and A3 (Eichler–Shimura labeling).
**Information gain:** This is the **most concrete recent computational tool** in the Hain–Brown lineage. Brown's 2014 coaction was structural; Kleinschmidt et al. 2026 give explicit closed-form formulae for the coaction on iterated Eisenstein integrals of low depth. If GeoVac wants to do A2 (the depth-2 PSLQ test) and have any chance of identifying *which* coaction structure is operating, Kleinschmidt et al. is the engine.
**Risk:** Kleinschmidt et al. is still "a proposal"; the first-principles match to the abstract motivic coaction is open. But for sprint-scale PSLQ tests, formal validity is enough.

### A9 — Chabauty–Kim–Kantor relative-completion framework (Best–Dogra et al. 2411.18846)

**Cost:** multi-year (full transport).
**What it adopts:** the unipotent Chabauty–Kim approach for relative completions extends to "arithmetic schemes over $\mathbb{Z}$" and includes Selmer-stack-type apparatus.
**GeoVac home:** Speculative — this is the arithmetic-substrate direction the BC-RH probe (CLAUDE.md §2) explicitly tabled.
**Information gain:** Park as **off the menu** until GeoVac opens an arithmetic-substrate parallel arc.
**Risk:** explicitly out of current scope; named to keep the literature search complete.

### A10 — Cuspidal-vs-Eisenstein dichotomy as inner-factor diagnostic

**Cost:** multi-month (3–4 months structural).
**What it adopts:** the Hain–Brown observation that cuspidal generators (labeled by cusp forms like $\Delta$) and Eisenstein generators (labeled by $E_4, E_6, \ldots$) play structurally distinct roles in $\mathcal{O}(\mathcal{G}^{\mathrm{rel}})$.
**GeoVac home:** Paper 18 §IV.6 inner-factor data tier; H1 Yukawa non-selection theorem (Sprint H1).
**Information gain:** GeoVac's outer factor is in the M1/M2/M3 partition; the **inner factor** (Yukawa Dirichlet ring) sits separately and the H1 sprint proved that GeoVac does not autonomously select Yukawa values. *If* the inner factor's data tier turns out to be cuspidal in the Hain–Brown sense, this would name the structural reason the values are calibration-external — they're cuspidal periods on a different geometric object than the M1/M2/M3 outer carrier. Highly speculative; deserves to be in the catalogue but not high priority.
**Risk:** the inner-factor program is still calibration-tier; this would mainly buy us a structural reading.

---

## §3. What GeoVac contributes back to the Hain–Brown program

The reciprocity direction is potentially more important than the adoption direction, because GeoVac has structural content that the modular / motivic community lacks. The cleanest contributions:

**R1 — Bit-exact Tannakian closure at finite cutoff.** Paper 56 establishes Deligne–Milne reconstruction at $n_{\max} \in \{1,2,3,4\}$ with **2,960 sympy-rational zero residuals**. The Hain–Brown program traditionally works at the level of the pro-finite limit and constructs Tannakian categories with *conjectural* fundamental groups (e.g., MEM$_1$ in Hain–Matsumoto is "a $\mathbb{Q}$-linear Tannakian category … contains MTM"; the **fundamental group is determined only at the level of low-weight relations and lowest-order Galois action**). GeoVac's per-cutoff Tannakian closure is a finite-resolution model **that explicitly closes the converse direction at each cutoff**. This is structurally different from anything in Hain–Brown.

**R2 — Master Mellin engine M1/M2/M3 as period-ring stratification.** Paper 32 §VIII case-exhaustion theorem proves every $\pi$ in any GeoVac observable comes from $\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$ at $k \in \{0, 1, 2\}$. **This is a structural classification of periods by operator order** — there is no analog in the Hain–Brown program. The closest published thing is the depth filtration on MZVs, but it operates differently: depth is a property of the iterated-integral *expression*, while $k$ is a property of the spectral *operator*. The two filtrations are likely related (depth-2 Mellin test is the candidate diagnostic) but the operator-order grading is GeoVac-original and a candidate addition to the period-classification toolkit.

**R3 — A concrete geometric realisation with explicit closed-form periods.** GeoVac produces **closed-form period rings**: $M_1 \subset \mathbb{Q}[\pi, \pi^{-1}]$, $M_2 \subset \bigoplus_k \pi^{2k}\mathbb{Q}$ (pure-Tate, smaller than generic Fathizadeh–Marcolli mixed-Tate), $M_3 \subset \mathrm{MT}(\mathbb{Z}[i, 1/2])$ at level $\le 4$. Hain–Brown's periods are computed as iterated integrals with closed-form expressions only at low depth; GeoVac's are bit-exact in sympy rationals at every cutoff. **For testing conjectures in the modular setting, GeoVac is a sandbox.** Specifically: the Pollack quadratic relations could be tested numerically on GeoVac (A4 above) before doing the full algebraic transport.

**R4 — A non-modular instance of the $(SL_2, \mathbb{G}_a^\infty)$ shape.** Hain–Brown's $SL_2$ comes from $\pi_1(\mathcal{M}_{1,1}) = SL_2(\mathbb{Z})$ — modular curves. GeoVac's $SL_2$ comes from Bertrand × Fock projection on $S^3 = SU(2) = \mathrm{Spin}(3)$ — Coulomb physics. The two routes to $SL_2$ are forced *differently*, but they produce algebraic-$SL_2$ over $\mathbb{Q}$ acting on $\mathrm{Sym}^k V_{\mathrm{fund}}$ on both sides. **The convergence is striking and may indicate a deeper underlying mechanism that selects the $SL_2$ shape for "natural" reductive factors of cosmic-Galois-like structures.** GeoVac is currently the only known non-modular example.

**R5 — Spectral-triple / operator-algebraic framing.** The GeoVac corpus has a fully developed math.OA standalone series (Papers 38–50) using propinquity convergence, Krein-space modular Hamiltonians, Lorentzian extensions. Hain–Brown traditionally works in pure number theory / algebraic geometry. **Translating the per-cutoff Tannakian closure into propinquity / spectral-triple language could open a new bridge between the Hain–Brown program and the operator-algebraic community.** Concretely: the Connes axioms verified on GeoVac at finite cutoff (Paper 32 §IV) are *exactly* the prerequisites for the Hain canonical MHS in a NCG setting — but Hain–Brown has never been stated NCG-style.

---

## §4. Honest comparison

**What Hain–Brown has and GeoVac structurally lacks:**

| HB has | GeoVac lacks | Why it matters |
|:---|:---|:---|
| Modular-form input to the generator labeling (Eichler–Shimura) | A natural map from modular forms to GeoVac's generators | Cuspidal vs Eisenstein dichotomy is a real distinction on HB's side, absent in GeoVac. Could be added (A3 above), multi-month. |
| Canonical mixed Hodge structure (Hain 2014, Theorem 3.8) | Not yet stated MHS-style; only verified at the Connes-axiom + reconstruction level | Stating the GeoVac MHS would be an editorial/clarifying sprint (A1). The Connes axioms supply the ingredients but the MHS is not currently named. |
| Pollack quadratic relations (explicit, depth-2, weight ≤ 14) | No relations beyond commutators in $\mathbb{G}_a^\infty$ | This is the structural distinction Reading A vs B (Paper 56 §sec:open_na1) — GeoVac's abelianness vs HB's relations is potentially diagnostic. |
| Hecke action compatibility (Hain 2023) | No Hecke action | Multi-year add-on; speculative whether even relevant. |
| Modular curves $\mathcal{M}_{1,1}$, $\Gamma_1(6)$, $\Gamma_0(N)$ as carriers | Carriers are spheres $S^3, S^5$ | Two genuinely different geometric substrates that happen to produce the same algebraic group shape. Possibly a hint, possibly a coincidence. |
| Iterated Eisenstein integrals as computational engine (Brown 2014, Kleinschmidt 2026) | Closed-form Mellin moments at each cutoff | Both are computational; HB's are *integrals*, GeoVac's are *Mellin transforms*. The depth-2 test (A2) would put them in the same arena. |

**What GeoVac has and Hain–Brown structurally lacks:**

| GeoVac has | HB lacks | Why it matters |
|:---|:---|:---|
| Bit-exact Tannakian closure at finite cutoff | Per-cutoff finite-resolution closure | Per-cutoff verification is GeoVac-novel methodology. |
| Master Mellin engine partition (operator-order $k$) | No analog | $k$-grading is GeoVac-original; potential contribution to period classification. |
| Physics witnesses (chemistry, QED, gravity, nuclear) | Only modular / Feynman-genus-1 witnesses | GeoVac has a sandbox of >100 physical observables; HB has a focused arithmetic / scattering-amplitude program. Different test surfaces. |
| Spectral-triple / Connes-axiom verification | No NCG framing | R5 above; potential bridge construction. |
| A non-modular route to $SL_2$ (Bertrand × Fock on $S^3$) | Only modular route ($\pi_1(\mathcal{M}_{1,1})$) | Two independent forced-by-geometry derivations of the same algebraic group. R4 above. |

**Symmetry of the comparison.** GeoVac has substrate-side richness (Connes axioms, propinquity convergence, master Mellin grading); Hain–Brown has period-side richness (modular forms, Hecke equivariance, established iterated-integral coaction). **Neither side dominates.** Both directions of adoption-and-contribution are legitimate.

---

## §5. Recommendations

**Top 2 candidates for next sprint cycle:**

### Recommendation 1 — A2 (iterated Eisenstein integral PSLQ test), merged with the strategic memo's "Test A"

**Why it's first.** A2 is *the same computational object* as Test A in the strategic memo (depth-2 Mellin of two M3 outputs vs the modular ring $\{ \zeta(3), \zeta(5), \beta(2k), L(\Delta, k), \pi^k \}$). The strategic memo already named Test A as the **highest-priority sprint-scale move** post-v3.77.0. This survey doesn't change that prioritization; it sharpens the test by:
- Naming the computational engine (Kleinschmidt et al. 2026 coaction formulae) — without this, the test is a black-box PSLQ; with it, structural interpretation of any positive identification becomes immediate.
- Naming what success and failure look like in HB-categorical language: success → GeoVac is in $\mathrm{MEM}_1$-like category, structural identification with HB; failure → GeoVac is in strict-$\mathcal{G}_4$ MZV cyclotomic land, no modular enrichment, and the Hain–Brown reading retires in favor of the Deligne 2010 / Glanois 2015 target already in Paper 56 §sec:open_g4.

**Estimated cost:** 1.5–2.5 weeks (one week sprint on the GeoVac Mellin side, half a week to import the Kleinschmidt coaction code or recompute, half a week PSLQ + writeup).

**Decision rule (binary).** If any $L(\Delta, k)$ or any non-MTM modular period appears in the PSLQ basis with low coefficient (≤ 100), Hain–Brown identification is empirically supported. Otherwise, retire to $\mathcal{G}_4$ comparison.

### Recommendation 2 — A1 (Hodge realisation framing) bundled into the post-Test-A writeup

**Why it's second.** A1 is cheap, editorial, and sets the framing for whatever Test A's outcome is. If Test A passes, A1 commits GeoVac to the Hodge realisation target. If Test A fails, A1 still locks in Hodge over étale (the right substrate for GeoVac regardless of whether the HB identification holds). **A1 is risk-free preparation work.** Bundle into Paper 56 §sec:open in the same sprint as Test A's writeup.

**Estimated cost:** 2–3 days, parallel to Test A.

**What to defer:** A3 (Eichler–Shimura labeling), A4 (Pollack relations test), A5 (algebraic de Rham apparatus) all become high-priority *if Test A passes*. If Test A fails, A3–A5 become low priority (the structural identification dissolves and the $\mathcal{G}_4$ comparison takes over).

**What to keep off the menu:** A6 (Hecke action), A9 (Chabauty–Kim–Kantor), A10 (cuspidal-inner-factor) — speculative or out of scope.

**What to keep in the back pocket:** A7 (Tapušković edge-subdivision) and A8 (Kleinschmidt computational engine) — these are sprint-scale tools to deploy *during* Test A's actual execution; named here as resources.

**Reciprocity moves (R1–R5).** R1 (bit-exact per-cutoff reconstruction) and R2 (master Mellin engine $k$-grading) are the cleanest community-facing contributions. **If Test A passes, R1–R5 become the natural slogans for a community talk or a NCG-meets-modular preprint.** Paper 56 currently does not advertise R1; explicit framing of "first finite-cutoff Deligne–Milne reconstruction in the Hain–Brown lineage" would be earned, and would land receptively at venues where MEM is discussed (see §6 below).

---

## §6. Researchers / groups active in this area, 2024–2026

**Hain (Duke).** Hain himself remains active. The 2023 Hecke-action paper (arXiv:2303.00143, Annales ENS) had revisions through June 2025. He has a YouTube lecture series (4 parts) on "Universal mixed elliptic motives" available on YouTube. *Reception likelihood for GeoVac:* moderate-high. Hain is sympathetic to physics-flavored work (cites string-amplitude papers in 2014). The Tannakian / NCG framing should land cleanly if presented as "discrete-graph realisation of a (sub-quotient of?) the Hain–Brown relative completion."

**Brown (IHES / Oxford).** Multiple modular values paper (2014) is current standard. His 2024 work continues in single-valued and motivic-coaction direction. He runs the IHES motivic-periods program (Eskandari–Murty–Nemoto 2025, already cited by GeoVac, is in this circle). *Reception likelihood:* very high. Brown has explicitly stated interest in NCG cross-connections in past talks. The "GeoVac as a sandbox with closed-form periods" framing (R3) would land especially well.

**Hain–Matsumoto.** Matsumoto (Hiroshima) collaborated on the 1512 paper; Hain–Matsumoto is the canonical MEM reference. Activity 2024–2026 less visible but Matsumoto's program on relations in mixed-elliptic Lie algebras continues. Reception unknown.

**Kleinschmidt (AEI Potsdam) + Broedel (Berlin) + Mafra et al.** Active 2024–2026 on string amplitudes / iterated Eisenstein integrals. Kleinschmidt et al. 2026 (arXiv:2508.02800) is the most recent computational tool from this circle. *Reception likelihood:* very high — they have a *physics audience* and are looking for new physics inputs to their motivic-coaction programme. The Tapušković 2023 paper is the precedent.

**Tapušković (Oxford).** Wrote the cosmic-Galois / sunrise paper (arXiv:2303.17534). Single-paper publication record visible at the listed Oxford page; his interest is the physics-input direction. *Reception likelihood:* high. If GeoVac wanted a sprint-scale collaboration on framing the Test-A outcome for the periods community, Tapušković is the right contact.

**Pollack (UCSD / IAS).** Multiple papers on zeta elements in depth, relations in $\mathfrak{u}$. *Reception likelihood:* moderate. Pollack's work is technical-relations-focused; GeoVac's contribution would be most natural at the *generator-classification* level (A3), not at the relations level (A4) where Pollack is the expert.

**Schneps, Baumard.** Active in Lie-algebra relations and moulds (Écalle theory) for mixed-elliptic motives. Reception unknown; cultural distance is moderate.

**Best, Dogra, Kantor.** Chabauty–Kim direction (arXiv:2411.18846). Arithmetic-substrate direction; off menu for GeoVac as currently constructed.

**Connes, Marcolli, Consani, Moscovici.** GeoVac's heritage line. The Hain–Brown identification doesn't fit into Connes–Marcolli's existing framing (their $U^*_{CM}$ has free non-abelian unipotent; the relation to MEM is via Tate-degenerate fibre per the closure-of-day memo). Marcolli has previously written on motivic Galois and gravity. Reception of "GeoVac is in the Hain–Brown lineage rather than (or one storey above) the Connes–Marcolli lineage" — probably positive, as it sharpens the comparison rather than competes.

**Aggregate landscape verdict.** The HB/MEM circle currently has *three* active sub-circles: (i) the original Hain / Brown / Matsumoto pure-arithmetic lineage; (ii) the string-amplitude / Kleinschmidt / Mafra physics-input lineage; (iii) the Chabauty–Kim arithmetic-applications lineage. **GeoVac sits most naturally adjacent to (ii)** because it's a physics-rooted framework, but with a structural angle on (i) via the Tannakian-closure machinery (Paper 56). The cleanest single-talk pitch for community engagement would be: *"GeoVac is a discrete spectral-triple cousin of Hain–Brown MEM with bit-exact finite-cutoff Tannakian closure and a master-Mellin operator-order grading $k = 0, 1, 2$ on its periods."* That pitch lands at venues like the IHES periods program (Brown), the AEI string-amplitudes program (Kleinschmidt), and possibly the math.OA / NCG meetings already in the Connes–Consani orbit.

---

## §7. Closing observation

The most striking finding of this survey: **the cleanest test of the Hain–Brown identification (Test A in the strategic memo, A2 + A8 in this memo) and the cleanest test of Reading A vs B in NA-1 (joint Mellin transform of M2 × M3 at depth 2) are computationally the same object.** A single 1.5–2.5 week sprint executes both, with two cleanly-separable decision rules (modular content present? → HB identification supported; symmetric vs asymmetric factorisation? → Reading A vs B). This is unusually efficient — one sprint disambiguates two named structural questions.

**Per the auto-mode bias and the strategic memo's recommendation, A2/A8 + A1 framing pass is the proposed next sprint.** Cost: 2.5–3 weeks combined. Outcomes: structural identification or disidentification with Hain–Brown plus structural Reading-A-vs-B verdict — either way, the foundations arc gets sharper.

---

## Sources

- Hain, *The Hodge–de Rham Theory of Modular Groups*, arXiv:1403.6443 (2014).
- Brown, *Multiple Modular Values and the relative completion of the fundamental group of $\mathcal{M}_{1,1}$*, arXiv:1407.5167 (v4 2017).
- Hain–Matsumoto, *Universal Mixed Elliptic Motives*, arXiv:1512.03975 (J. Inst. Math. Jussieu 19 (2020) 663–766).
- Hain, *Algebraic de Rham theory for relative completion of $SL_2(\mathbb{Z})$*, arXiv:1804.06977.
- Hain, *Hecke Actions on Loops and Periods of Iterated Shimura Integrals*, arXiv:2303.00143 (Annales ENS, rev. June 2025).
- Pollack, *Zeta elements in depth 3 and the fundamental Lie algebra of a punctured elliptic curve*, arXiv:1504.04737 (Forum Math. Sigma 8 (2020) e31).
- Baumard–Schneps, *On the derivation representation of the fundamental Lie algebra of mixed elliptic motives*, arXiv:1510.05549 (Ann. Math. Québec).
- Tapušković, *The cosmic Galois group, the sunrise Feynman integral, and the relative completion of $\Gamma_1(6)$*, arXiv:2303.17534 (Comm. Number Theory Phys. 18 (2024) no. 2). **Note: GeoVac memory file [[hain_brown_identification]] and the v3.77.0 strategic synthesis memo misattribute this to "Bouillon" — flagged for next bib pass.**
- Kleinschmidt et al., *Towards Motivic Coactions at Genus One from Zeta Generators*, arXiv:2508.02800 (JHEP 05 (2026) 105).
- Baune–Broedel–Moeckli, *A construction of single-valued elliptic polylogarithms*, arXiv:2511.15240 (Nov 2025).
- Best–Dogra et al., *The Unipotent Chabauty-Kim-Kantor Method for Relative Completions*, arXiv:2411.18846 (Nov 2024).
- Sun, *Equivariant primitives of Eisenstein series for congruence subgroups*, arXiv:2502.04752 (Feb 2025).
- Eskandari–Murty–Nemoto 2025, arXiv:2510.20648 (forces level 4 for Catalan $G$ — already cited by GeoVac).
- Deligne, *Le groupe fondamental unipotent motivique de $\mathbb{G}_m - \mu_N$*, Publ. Math. IHÉS 112 (2010) 101–141.
- Glanois, PhD thesis, UPMC 2015.
- Brown, *Mixed Tate motives over $\mathbb{Z}$*, Ann. Math. 175 (2012) 949–976.
- GeoVac internal:
  - Paper 55 `papers/group3_foundations/paper_55_periods_of_geovac.tex` (period ring trinity).
  - Paper 56 `papers/group3_foundations/paper_56_tannakian_substrate.tex` (Tannakian dual + §sec:open).
  - Paper 32 §VIII (case-exhaustion theorem, master Mellin engine).
  - Paper 18 §III.7 (M1/M2/M3 partition).
  - `memory/hain_brown_identification.md` (memory file).
  - `debug/strategic_synthesis_2026_06_06_memo.md` (closure-of-day, §Addendum).
  - `debug/sprint_q5p_na1_non_abelian_probe_memo.md` (NA-1 abelian-by-CMM finding).
