# Email draft — outreach to Carlos Perez-Sanchez

**To:** Carlos Perez-Sanchez (Heidelberg, Faculty of Physics & Astronomy)
**From:** Josh Loutey (independent researcher, GeoVac project)
**Subject:** Infinite-dimensional vertex extension of Bratteli-network spectral actions

---

Dear Dr. Perez-Sanchez,

I am writing about your recent work on Bratteli-network spectral actions on quivers (arXiv:2401.03705, arXiv:2409.03705, arXiv:2508.17338). I would value your view on a specific extension question.

I lead an independent research project (GeoVac: github.com/jlout/<repo>) that constructs molecular quantum Hamiltonians by composing atomic spectral triples — one Camporesi-Higuchi (CH) spinor spectral triple on S³ per atom, with charge-Z-dependent scale, coupled via multipole-expansion bimodules on edges of the molecular bonding graph. The construction sits in the Marcolli-van Suijlekom 2014 lineage, and the framework's gauge-theoretic and operator-algebraic content is documented in a series of math.OA papers (Papers 38, 39, 42-49 in the project, on propinquity convergence of truncated SU(2) spectral triples, Lorentzian extensions, and related results).

Your Bratteli-network construction looks like a *very* natural home for our chemistry construction at the spectral-action level — vertices = atoms, edges = bonds, bond bimodules = multipole-expansion couplings. The combinatorial structure matches almost verbatim.

The obstruction I see, and the question I would like to ask you, is **infinite-dimensionality at vertices**. Our atomic CH spinor sectors are intrinsically infinite-dimensional; we truncate at a cutoff n_max for computational practice, but the underlying mathematical object is the full infinite-dim CH bundle on S³. Reading your 2024 papers, the prespectral-triple definition itself is permissive (just (A, H) with A a *-algebra faithfully *-represented on an inner-product space H, no completion required), but:

1. The spectral action Tr f(D/Λ) is computed as a literal finite-matrix trace.
2. The loop-equation framework of arXiv:2409.03705 uses Haar integrals over compact U(N_e); for infinite-dim H_v the unitary group is not locally compact.
3. The Bratteli combinatorics requires finite multiplicities.

**My specific question:** Does your Bratteli-network spectral-action framework admit an extension to vertex prespectral triples (A_v, H_v) with H_v infinite-dimensional, provided the global Dirac D on H = ⊕_v H_v has compact resolvent and the spectral action is replaced by a ζ-regularized or heat-kernel-regularized trace? If so, what are the obstructions to the loop equations of arXiv:2409.03705 transporting? And is there a natural "infinite-dim Bratteli" combinatorial parameterization that I should be looking at, perhaps via Brain-Mesland-van Suijlekom unbounded Kasparov products?

A secondary question, since our bond bimodule is Hermitian but not unitary (it is the standard multipole-expansion cross-V_ne — a fixed Hermitian operator, not norm-preserving): is the unitary-intertwiner axiom on edges (φ_e, U_e) in your 2024a definition strict, or does it admit a relaxation to Hermitian bimodules in the inner-fluctuation sense of Connes-Chamseddine?

I would be very happy to share the relevant project context (a brief technical memo and citation list) if that would help shape your response. The decision this informs is whether we invest several months in proving a "GeoVac chemistry composition is a Bratteli network spectral action" theorem; a positive direction from you would make this concrete sprint-scale work; a structural obstruction would save us the investment cleanly.

Thank you for your time. Looking forward to your perspective.

Best regards,
Josh Loutey
GeoVac project (independent research)
jloutey@gmail.com

---

## Drafting notes (not part of the email)

- The email frames GeoVac as a project in the Marcolli-vS lineage (which it is, via WH1) rather than a novice ask. This is honest and signals we've done the homework.
- The specific question is taken from the Bratteli deep-dive agent's recommendation, lightly polished. It's narrow enough to answer in a paragraph if the answer is structural, and detailed enough to scaffold a more substantive reply if Perez-Sanchez is interested.
- The "share project context" offer is real — we can ship Papers 32, 38, 45-49 if he asks; the PI's call on what to include.
- Affiliation: per the Bratteli agent, Perez-Sanchez is at Heidelberg. Verify his current email address via his homepage before sending.
- Tone: technically substantial, not deferential, ends with a clear ask.

## Action

PI to review, edit if desired, then send. Verify current email at Heidelberg before sending.
