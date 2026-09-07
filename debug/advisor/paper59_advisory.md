# Advisory — Paper 59: The Three-Center ERI Is an Elliptic Bessel Moment — 2026-09-06

*Charter: `.claude/commands/advisor.md`. This is the /advisor pass (committee-chair judgment),
NOT /qa (no PASS verdict, no line-level defect gate). Advice only — the PM executes or declines.*

## The contribution, in one sentence

The genuine three-center Slater electron-repulsion integral, read in momentum space (the Fock
projection where three centers become three phases, not three foci), reduces to a two-scale
Bessel moment on a genus-one (Legendre/Γ(2)) elliptic curve — one genus above the two-center
{E₁, ln, γ} class — and the natural sunrise/elliptic-dilogarithm route to its closed form is
*proven obstructed*, because the object is a period of an irreducible rank-four irregular
connection; the finite closed form is left honestly OPEN. The message is a relocation, not a
solution: the third center raises the transcendence one genus.

That thesis is real, non-trivial, this paper's own, and — critically — it *is* in the abstract,
stated as its last sentence. The spine exists and is well-chosen.

## Verdict

**DEFENSIBLE AFTER NAMED FIXES.** The central contribution would survive a committee: the genus
argument is sound, the irreducibility result has a clean in-principle proof (Katz Fourier–Laplace
auto-equivalence of holonomic D-modules) backed by an explicit integer-monodromy computation, the
tier discipline ([MEASURED]/[SYMBOLIC]/[OBSERVATION]/[OPEN]) is doing genuine load-bearing work,
and the negatives are honest. None of the issues below touch the validity of the contribution —
they are altitude, audience, and hedge problems. But two of them are real defensibility exposures
(one factual characterization stated as fact on paywalled evidence; one abstract phrase stronger
than the body supports), and the paper has grown a very large number-theory tower (Secs. 6–8) that
is where it is simultaneously most exposed and least certified — the /qa cert is OWED precisely on
that newest material. The paper is close to as-is; a short, cheap list of framing fixes gets it
there.

## What is genuinely strong

- **The reframing is a genuine bridge, and it is correctly positioned as an OBSERVATION, not a
  theorem.** No chemistry paper has read an ERI as an elliptic period/Bessel moment; no amplitudes
  paper has pointed the elliptic machinery at a molecular integral. The paper's own targeted-search
  hedge ("an absence from a targeted search, not a proof of absence; the residual risk is an obscure
  single source") is exactly the right altitude, and it is corroborated by the corpus's independent
  adversarial lit-scan (`debug/lit_scan/elliptic_eri_feynman_memo.md`: OPEN frontier, disjoint
  literatures).
- **The honesty on the closed form is exemplary.** The abstract does not claim the closed form; it
  marks it [OPEN] and then does the harder, more valuable thing — shows *why* the obvious route
  (unequal-mass sunrise / elliptic dilogarithm) fails, via the in-module computation (the source
  stays on the elliptic curve; the sunrise's tadpole boundary term vanishes identically). This is a
  negative result with a pinned mechanism — the FORCED/FREE/WALL vocabulary of the mission at its best.
- **The two-center COLLISION is handled cleanly.** The paper leads with "the two-center problem has
  been solved since Roothaan and Ruedenberg," cites Harris–Michels / Barnett–Coulson / the
  Averys correctly as classical baseline, and claims no novelty there. The one subtle distinction —
  *decidable* (the AA|BB block, Paper 58) vs merely *counted* (the cross classes) vs *weight-one* — is
  stated precisely and not overreached. This is the kind of restraint that earns a reader's trust.
- **The irreducibility argument is the mathematically deepest and most defensible result.** Two
  independent routes (the Katz FL "simple-goes-to-simple" argument and the explicit exact-integer
  monodromy M₀) converge, and the paper is scrupulous about tier: the M₀-dependent clauses carry a
  "[SYMBOLIC + MEASURED] … obtained by computation, not derived in closed form" scope note rather
  than being dressed as pure theorem.

## The hostile-examiner question

**"You show N(D) is the Laplace transform of the now-fully-solved unequal-mass sunrise
(Bogner–Müller-Stach–Weinzierl 2019), and you offer 'Laplace-transform that solution directly' as
'a concrete route … we have not carried out.' Until that route is either executed or proven not to
terminate in closed form, isn't your headline negative — 'past the reach of the sunrise machinery'
— really the weaker statement 'the *leading elliptic-dilogarithm mechanism* does not engage, and we
have not finished the harder dual route'? What is actually OPEN here: a genuinely new transcendent,
or a known object nobody has yet written down?"**

Honest status: **(b), a legitimate scope boundary that the body already names — but the abstract
overstates it.** The paper *does* make the correct careful distinction (the Fuchsian all-orders
sunrise vs its irregular Laplace dual; irregularity established symbolically via the Newton polygon
/ Poincaré rank one). The irreducibility result is a real theorem that survives this question
untouched — no factorization gives a closed form regardless of the dual. What does *not* fully
survive is the abstract's phrase "past the reach of the sunrise machinery," which reads as a proven
wall when the body's honest position is "the leading mechanism is obstructed; the irregular dual is
a named, uncarried-out calculation." The fix is to align the abstract to the body (see MUST-FIX 2),
after which the answer to the examiner is clean: the reducibility question is *settled* (irreducible),
the transcendental closed form is *open*, and the sunrise's leading mechanism is *proven* not to
supply it.

## Recommendations (ranked, must-fix first)

1. **[MUST-FIX] Give Avery 2013 its own flag-planting sentence, and hedge "stops at expansion and
   numerics" to match the paper's own novelty-hedge.** (Sec. I "Lineage and scope," ~L137–140; Sec.
   VI literature, ~L684–686.) Avery, *"Fast electron repulsion integrals for molecular Coulomb
   Sturmians"* (Adv. Quantum Chem. 67, 2013) is — per the corpus's own adversarial lit-scan — the
   single most dangerous prior-art paper: same Fock-projection momentum route, same objects, and the
   one place a reader could have noticed the elliptic structure. Right now it is folded into compound
   `\cite{avery_avery_4center,avery2013}` cites and never named as *the* closest stander. Worse, the
   characterization "which stops, as all prior work does, at expansion and numerics" is stated as
   established fact, while the lit-scan flags that the "no elliptic reading" reading rests on
   *paywalled abstracts* (ScienceDirect 403; titles confirmed, absence-of-genus inferred). A committee
   member from the Sturmian camp could produce the chapter. Two acceptable fixes: (a) fetch the
   chapter and confirm, or (b) name Avery 2013 explicitly ("the route on which the elliptic structure
   was closest to visible") and soften the characterization to the same [OBSERVATION]-tier search-hedge
   the abstract already uses for the broader novelty claim. This is fair *to Avery* and *protects the
   paper* — it costs one sentence.

2. **[MUST-FIX] Tighten the abstract's final clause to match the body's obstruction statement.**
   ("past the reach of the sunrise machinery," abstract L71–72.) The body is careful: the *leading*
   elliptic-dilogarithm mechanism is proven obstructed (in-module source), and the *irregular
   Laplace dual* of the fully-solved sunrise is a named-but-uncarried-out route. The abstract's phrase
   flattens that into a proven wall. Reword to something like "past the reach of the sunrise's
   *elliptic-dilogarithm* mechanism, whose inhomogeneity is absent here" (which is exactly what Sec. V
   proves) — leaving the dual route honestly open. This is a tier-consistency fix between abstract and
   body, not a retreat: the irreducibility theorem is unaffected.

3. **[STRENGTHEN] Resolve the audience/altitude tension: decide whether Secs. 6–8 are this paper or a
   companion.** The paper carries two full-altitude theses — the chemistry-facing relocation (Secs.
   1–5) and a number-theory-facing structural tower (Sec. 6 modular/CM/cosmic-Galois, Sec. 7
   Bessel-moment period algebra / intersection form B=πΩ, Sec. 8 F12). No single reader is fluent in
   all three of {STO molecular integrals, modular forms + D-modules, resurgence}. A chemistry reader
   is lost by Sec. 6; an amplitudes reader must wade through chemistry to reach Sec. 7. The spine
   exists but is overgrown. Two options to PROPOSE (PI's call — this touches the deliberate "atlas"
   completeness the mission values, so it is a recommendation, not a blocker): (i) compress Secs. 6–8
   to the load-bearing minimum the thesis needs (Legendre/Γ(2) home, CM fibers visited, irreducible ⇒
   not a repackaged genus-0 object, closed form open) and move the rest to a companion; or (ii) split
   the number-theory tower into a second paper aimed squarely at the periods/amplitudes community
   (Broadhurst/Vanhove/Weinzierl/Fresán–Sabbah–Yu; and per the corpus's own note, Brown/Kleinschmidt).
   Either way the chemistry paper lands harder and the number-theory paper gets a referee who can
   actually evaluate the intersection-form and resurgent-skeleton claims.

4. **[STRENGTHEN] Cut or externalize the "five-for-five resurgent skeleton" program paragraph.** (Sec.
   4/5, the long [OBSERVATION] on the exchange class, integer Stokes charges, γ-tagging across *five*
   objects.) Only one of those five objects — N(D) — belongs to this paper; the other four (the P18
   seed, the P58 hybrid class, the exchange class, the second-cusp triple) are corpus-program results
   a reader cannot verify from this paper. It reads as a program-status update embedded mid-argument
   and dilutes the thesis. Keep N(D)'s own resurgence (which is load-bearing for the irregularity
   claim); move the cross-object "5/5 pattern" to the synthesis/CHANGELOG, or to the companion.

5. **[STRENGTHEN] Disambiguate the overloaded ρ.** ρ means the elliptic modulus c₂/c₁ in Secs. 3–5,
   but in Sec. 6 the same glyph is reused for the corner radial variable ρ=s+t (with ρ=σ²). A reader
   tracking "the modulus ρ" through the modular section hits a different ρ without warning. Rename one
   of them (the corner variable is the cheaper rename).

6. **[CONSIDER] Add a concrete "so what for chemistry" to the scope section.** Sec. IX honestly
   concedes the reduction "does not by itself deliver a polyatomic energy" and that "the Gaussian fits
   are already exact to fit quality." That honesty is correct and should stay — but as written the
   positive payoff is only "a precise diagnosis of why it is hard." If there is an operational claim
   ("an accurate, fast, geometry-native evaluator"), make it concrete: one comparison to the chemistry
   state of the art (Özdoğan–Ruiz: 20 digits in 25–30 terms), or a timing/accuracy-vs-cost line. Absent
   that, state plainly that the deliverable is diagnostic + the bridge to elliptic-Feynman technology +
   the F12-native observation — so a chemistry reader is not left expecting a capability the paper does
   not claim. Also motivate the T2 = 0.3953… value for what it is (a number-theory PSLQ target on one
   fixed geometry), not as a chemistry result.

## Positioning check

- **Altitude on the OPEN frontier — right, with two phrase-level drifts.** The abstract does *not*
  overclaim a solved closed form (it marks it [OPEN] and shows the obstruction) and does *not*
  undersell the genuine first (the elliptic bridge is stated as a real, if hedged, novelty). This is
  the correct altitude on the load-bearing question. The two drifts are both toward slight overclaim
  and both cheap: "past the reach of the sunrise machinery" (MUST-FIX 2) and "a genuinely new elliptic
  transcendent" (L63) — the latter is a strong noun phrase at [OBSERVATION] tier; consider "an
  apparently unclaimed elliptic transcendent," tying it to the search-negative the way Sec. VI already
  does.
- **Avery flag-planting — fair in substance, not prominent, and one characterization is stated as fact
  on unverified evidence.** See MUST-FIX 1. The *substance* is honest (Avery is cited, called closest
  in spirit, and correctly said to stop at numerics *if* the abstracts are representative); the
  *prominence* and the *evidential hedge* need work.
- **Overclaim, elsewhere: minimal.** The tier tags are pervasive and mostly accurate. The F12 section
  is honestly scoped ("we do not treat the many-electron RI machinery, nor a resource comparison with
  tuned Gaussian F12"). The intersection-form B=πΩ and Sp₄(ℤ) containment carry an explicit tier-scope
  note that they inherit M₀'s [SYMBOLIC + MEASURED] status. Good.
- **Undersell: two crisp results are buried.** (a) "The transcendence class is set by the *number of
  two-center transition densities* — zero or one gives genus zero, two gives genus one" (Sec. VIII) is
  a clean, general, defensible statement that currently sits at the end inside the F12 section; it is
  arguably a better one-line summary of the whole paper than the abstract's current framing. (b) The
  in-module obstruction mechanism (the vanishing tadpole boundary term) is the paper's most original
  *positive* insight and is somewhat lost inside a dense Sec. V.
- **Corpus tier/rhetoric (§1.5): compliant.** No ontological-priority language; the momentum route is
  correctly credited as classical (Shibuya–Wulfman, Niukkanen, the Averys), not GeoVac-invented; the
  Paper-18/34 transcendence-taxonomy framing is internal and labeled as such. The "genus-graded
  embedding tier" claim (Sec. III) is presented as a structural reading, appropriately.

## What I don't know

- **The certification gap is real and current.** `check_cert_staleness.py` reports paper_59 cert =
  OWED (last cert 2026-08-21; the paper was edited after). The Secs. 6–8 material — the intersection
  form B=πΩ, the co-area/Γ₀(2) reduction, the resurgent-skeleton 5/5, the 66-digit T2 — post-dates the
  last /qa and has *not* been through the adversarial correctness gate. I am advising on framing, not
  certifying correctness; a fresh /qa (or at least a claims+code reviewer pass) on Secs. 6–8 is owed
  before those claims are load-bearing anywhere else. The CHANGELOG chronicle (v4.103.x) does back the
  headline numbers, so the paper is consistent with the record — but "chronicled" is not "certified."
- **The Avery 2013 / Avery–Avery 2015/2017 "no elliptic reading" rests on abstracts.** I could not
  confirm the full-text characterization independently (the corpus's own lit-scan hit ScienceDirect
  403). If the paper wants MUST-FIX 1 discharged by confirmation rather than hedge, someone needs
  library/authenticated access to the chapters.
- **Whether the Laplace-transform-of-the-solved-sunrise route actually terminates** is the load-bearing
  open problem and is beyond a reading pass — it is the collaboration frontier the paper itself names.
  My advice assumes it stays OPEN; if it were carried out (either direction) the whole framing of the
  central negative would change, and the abstract would need rewriting accordingly.
- **The deepest structural claims (irreducibility forcing, B=πΩ forced-not-fitted, Sp₄(ℤ))** lean on
  M₀, an integer matrix "obtained by computation." The paper is honest about this. Whether a
  differential-Galois specialist would accept the M₀-forcing arguments as proof, or want the in-principle
  D-module route to carry the full weight alone, is a domain-expert call I cannot make.
