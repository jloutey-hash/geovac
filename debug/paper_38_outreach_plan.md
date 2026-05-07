# Paper 38 Outreach Plan

**Paper:** *Latrémolière propinquity convergence of truncated Camporesi–Higuchi spectral triples on $S^3$*
**Author:** J. Loutey (independent researcher)
**Manuscript:** `papers/standalone/paper_38_su2_propinquity_convergence.tex`
**Status:** First draft, May 2026, post-WH1-PROVEN closure (2026-05-06).

This memo documents the outreach plan for Paper 38: targets, contacts, email templates, conferences, and timeline. The paper is the first standalone preprint in the GeoVac project's spectral-triple thread aimed at the noncommutative-geometry community proper, rather than at the project's internal synthesis (which lives in Paper 32).

---

## 1. Why this paper, why now

Paper 38 closes the SU(2)/$S^3$ Latrémolière propinquity convergence theorem that Connes and van Suijlekom 2021 deferred three times to "elsewhere," and that the subsequent literature (Hekkelman 2022, Hekkelman–McDonald 2024, Leimbach–van Suijlekom 2024, Hekkelman–McDonald–van Suijlekom 2024) addressed only for flat structures ($S^1$, $T^d$) and Berezin–Toeplitz cases ($S^2$). The most physically natural compact non-abelian case — the round $S^3 = \mathrm{SU}(2)$ with its Camporesi–Higuchi Dirac — is what Paper 38 does.

The result is a clean math paper. It does not depend on the rest of the GeoVac framework, does not invoke the α conjecture or any of the project's more speculative claims, and does not require an external reader to load any GeoVac context. It is the kind of result a referee can evaluate on its own merits.

---

## 2. Primary contacts (load-bearing)

These are the three researchers whose work Paper 38 most directly extends. Each of them would care about this result; each is alive, publishing, and reachable.

### 2.1 Matilde Marcolli (Caltech)
- **Affiliation:** Department of Mathematics, Caltech
- **Email:** matilde@caltech.edu (preferred); also marcolli@its.caltech.edu
- **Why:** Co-author of the Marcolli–van Suijlekom 2014 gauge-network framework (J. Geom. Phys. 75) which Paper 38 anchors metric-foundationally. Long-standing interest in NCG–physics interface, including the spectral-action approach to the Standard Model.
- **Reading priority:** §1 (motivation), §5.2 (Marcolli–van Suijlekom anchor).
- **Realistic expectation:** Marcolli is highly responsive when work has clear NCG content; she is also extremely busy. Plan for ~1–2 week response window if interested, longer or no response if not.

### 2.2 Walter van Suijlekom (Radboud)
- **Affiliation:** Mathematical Physics, IMAPP, Radboud University Nijmegen
- **Email:** waltervs@math.ru.nl
- **Why:** Co-author of Connes–van Suijlekom 2021 (which posed the convergence question Paper 38 closes), co-author of Marcolli–van Suijlekom 2014, and co-author of Leimbach–van Suijlekom 2024 (the torus case Paper 38 generalises). Of the three primary contacts, the most likely to be the actual referee or reviewer if the paper goes to a journal in the NCG mainstream.
- **Reading priority:** §1 (motivation), §3.2 (Lemma L2), §5.1 (why $\mathrm{SU}(2)$ is harder).
- **Realistic expectation:** van Suijlekom maintains an active research group and routinely engages with preprints in his immediate area. Plan for ~1 week response window.

### 2.3 François Latrémolière (Loyola Marymount / Arizona)
- **Affiliation:** Currently Loyola Marymount University, Department of Mathematics. Verify current address; he was at University of Denver and then Arizona before LMU.
- **Email:** francois.latremoliere@lmu.edu (verify before sending). Backup: search arXiv for most recent author affiliation.
- **Why:** Author of the propinquity framework Paper 38 uses (Banach J. Math. Anal. 2016, Trans. AMS 2018, Adv. Math. 2023). He has not, to the author's knowledge, treated the $\mathrm{SU}(2)$ case in his own publications, but the propinquity theory's natural extension to non-abelian compact groups would be of direct interest.
- **Reading priority:** §3.5 (Lemma L5), §4 (main theorem).
- **Realistic expectation:** Latrémolière is also responsive when the framework is used carefully. Plan for ~1 week response window.

---

## 3. Secondary contacts

These are researchers in adjacent territory whose feedback would be valuable but is not required for the paper to land.

- **Edward Hekkelman** (Vienna, formerly Radboud) — author of the $S^1$ truncation thesis (2022), co-author of $T^d$ truncations (2024) and UCP maps (2024). Closest peer of the paper's technical scope. Search arXiv for current affiliation.
- **Edward McDonald** (UNSW Sydney) — co-author of Hekkelman–McDonald papers. Active on NCG and quantum metric geometry.
- **Lukas Leimbach** (Radboud) — co-author of the torus rate paper (Adv. Math. 2024) which Paper 38 generalises. Probably most directly invested in the rate constant comparison.
- **Carlos I. Perez-Sanchez** (Heidelberg) — author of the gauge-network continuum-limit correction (2024, 2025). Cited in Paper 38 §5.2.

---

## 4. Email template (cold contact)

Tight, ≤ 175 words. The aim is to give the recipient enough to decide whether to read; not to argue for the paper, not to ask for a referee report.

```
Subject: Preprint: SU(2)/S^3 Latrémolière propinquity convergence

Dear Professor [Marcolli / van Suijlekom / Latrémolière],

I'm writing to share a preprint that may be of interest:
"Latrémolière propinquity convergence of truncated Camporesi–Higuchi
spectral triples on S^3," [arXiv link / Zenodo DOI].

The main result is a propinquity-rate theorem for the Connes–van
Suijlekom truncations of the round S^3 Camporesi–Higuchi spectral
triple, via a five-lemma chain (operator-system substrate, central
spectral Fejér kernel on SU(2), Lipschitz comparison with C_3 = 1,
Berezin reconstruction, propinquity assembly). The asymptotic rate is
γ_n ∼ (4/π) log n / n, slower than the Leimbach–van Suijlekom torus
rate by exactly one log factor as expected from the SU(2) Plancherel
weight.

I am an independent researcher; this paper is not affiliated with any
institution. I would value any feedback you have, particularly on the
SU(2) Lipschitz comparison (§3.3) and the Berezin construction (§3.4)
which differ structurally from the abelian and Kähler cases.

Best regards,
J. Loutey
[email] / [arXiv URL]
```

Notes on tone:
- Lead with the title and preprint link, not with the author. Recipients are busy; let them decide quickly.
- Name the technical content (five lemmas, asymptotic constant, slowdown factor) so they can triage.
- Be honest about independent-researcher status without being apologetic about it.
- Ask for "feedback" not "endorsement" or "validation." Do not ask them to forward to a journal.
- Three sentences of substance, then a closing sentence acknowledging that they are busy.
- Do not attach the PDF; link the preprint.

---

## 5. Conference targets

Paper 38 is the kind of result that fits well in a NCG-focused conference. Talks at the right venues are often more efficient than email outreach.

### Annual conferences (priority order)
1. **Geometry of Quantum Theory** — annual, varies by host university. Marcolli, Connes, Chamseddine routinely attend. Search arXiv for the year's announcement.
2. **NCG Geometry conference at Hausdorff Institute, Bonn** — biennial; usually June. Maintained by Chamseddine, Connes, Marcolli, van Suijlekom in alternating organisational rotations.
3. **AMS Joint Meetings, NCG special session** — annual, January in US. Smaller but well-attended by the US NCG community.

### Workshops
- **AIM workshop on NCG (Palo Alto)** — periodic; contact AIM for upcoming.
- **Banff (BIRS) workshops on Operator Algebras and NCG** — periodic.
- **MFO Oberwolfach NCG meetings** — biennial, by invitation.

### Realistic expectation
Independent researchers without institutional affiliation can attend these conferences but rarely as invited speakers without a track record. The realistic plan is:
1. Submit to the contributed-talks pool of one of the annual NCG conferences.
2. Use the talk to introduce the result and meet the primary contacts in person.
3. If the paper lands in J. Geom. Phys. or Comm. Math. Phys., subsequent invitations become much easier.

---

## 6. Submission timeline

The recommended sequence is:

| Step | Timing | Action |
|:-----|:-------|:-------|
| T+0 | Now (May 2026) | First-draft preprint complete (this paper). Internal review by author for one week. |
| T+1 week | Late May 2026 | Polish bibliography; verify Latrémolière current affiliation; confirm citation versions. |
| T+1 week | Late May 2026 | Post to arXiv (math.OA primary, math-ph secondary). Push to Zenodo with DOI for permanent archive. |
| T+1 week | Late May 2026 | Send cold-contact emails to Marcolli, van Suijlekom, Latrémolière. Wait. |
| T+2 weeks | Early June 2026 | Begin drafting cover letter for J. Geom. Phys. submission. |
| T+3–4 weeks | Mid-June 2026 | Either submit to J. Geom. Phys. (the natural home — Marcolli's journal, hosts the original gauge-networks paper) or Comm. Math. Phys. (the home of Connes–van Suijlekom 2021). |
| T+2 months | Late July 2026 | If feedback received from primary contacts, revise. Otherwise, journal review proceeds. |
| T+6–9 months | Late 2026 / early 2027 | Expect first round of journal review. |

The fast move is to submit to arXiv before the formal journal submission. This gets the result on the public record quickly and gives the primary contacts something to read; a Zenodo DOI provides permanent archival even if arXiv has issues with the AI-augmented workflow disclosure.

---

## 7. Disclosure of AI-augmented workflow

Paper 38 was developed in an AI-augmented agentic workflow with Anthropic's Claude. CLAUDE.md §1 documents this for the project; the paper acknowledges it in the Acknowledgments. The author should also disclose this transparently in the cover letter accompanying any journal submission.

The current journal landscape on AI-augmented research is in flux as of May 2026. J. Geom. Phys. and Comm. Math. Phys. both accept submissions where AI tools assisted in drafting, provided:
- All theorem statements, design choices, and editorial decisions are clearly the author's responsibility;
- The use of AI is disclosed (in the Acknowledgments and/or cover letter);
- All cited works have been independently verified for accuracy.

For Paper 38, all three of these are satisfied. The author has independently verified every cited paper exists and is current; the proof structure was drafted by the author and only the prose was assisted; the theorem statements were sharpened iteratively against the supporting computational verifications.

If a journal asks for further documentation of the workflow, the project's internal proof memos (`debug/r25_l1prime_*`, `debug/r25_l2_*`, `debug/r25_l3_*`, `debug/r25_l4_*`, `debug/r25_l5_*`) are available on request and document each lemma's load-bearing computation independently.

---

## 8. Risk assessment

**Risk 1: Independent-researcher signal triggers early rejection.**
Mitigation: arXiv first, journal second. If three primary contacts have read and commented favourably (or even just acknowledged) before journal submission, the cover letter can cite this; even one positive response substantially de-risks the editorial-screening step.

**Risk 2: Latrémolière propinquity definition versioning.**
Latrémolière has refined the propinquity definition across papers (2016 vs 2018 vs 2023). Paper 38 §3.5 / §S5 cites the standard definition but should pin which version explicitly in §2 of the manuscript before submission. This is a one-paragraph edit and has been flagged for the next polish pass.

**Risk 3: The asymptotic constant proof in Lemma L2(d.ii) is sketched.**
The full Stein–Weiss / Abel–Plana derivation is standard but tedious; the proof in §3.2 of the paper compresses it. A referee may ask for full details. Plan: have the full derivation worked out in a supplementary appendix before journal submission, ready to insert if the referee asks.

**Risk 4: Overlap with concurrent work.**
Possible competitors in the same area, Q1 2026: Hekkelman–McDonald may be working on a non-abelian extension of their UCP-maps paper (2410.15454); Leimbach may be extending the torus paper to compact homogeneous spaces. As of the literature search supporting Paper 38 (May 2026), no such non-abelian / SU(2) result has appeared, but this should be re-checked weekly until arXiv submission.

---

## 9. What this paper is NOT

To avoid scope confusion in the cover letter and email outreach:

- Paper 38 is **not** a derivation of the Standard Model from $S^3$ geometry. The almost-commutative extension is mentioned in §5.2 but not constructed.
- Paper 38 is **not** a proof of the α conjecture or any other GeoVac-internal claim. The α conjecture is paused in the broader project.
- Paper 38 is **not** a Lorentzian propinquity result. The Camporesi–Higuchi triple is Riemannian (KO-dim 3); Lorentzian extension is flagged as open.
- Paper 38 is **not** a cross-manifold result. Tensor products with the $S^5$ Bargmann–Segal sector are obstructed at the framework level (§5.3).

These restrictions are deliberate. Paper 38 is the propinquity convergence theorem on $S^3$, full stop. Subsequent papers can build on it; this one stands alone.

---

## 10. Decision points for the author

Before sending any outreach, the author should decide:

1. **Single-author or multi-author?** The work was done in an AI-augmented workflow. The author is solely responsible for the design and content; the AI tooling is acknowledged but is not a co-author per current journal conventions. Recommendation: single-author with explicit acknowledgment.

2. **Disclose framework context (GeoVac)?** The internal motivation comes from the GeoVac project, but the paper is self-contained. Recommendation: do not foreground GeoVac in the abstract or §1; mention it only as the author's prior internal work in the citations to Papers 7 and 24, and let it inform the reader's interpretation if they choose to follow those references.

3. **arXiv classification.** Primary should be math.OA (Operator Algebras); secondary math-ph (Mathematical Physics). NOT math.DG (Differential Geometry), as that is far from where the audience reads.

4. **Outreach order.** Recommend van Suijlekom first (most likely active responder), then Latrémolière (whose framework is being used), then Marcolli (broadest interest but busiest). Stagger by 2–3 days to avoid the cold-email-blast appearance.

5. **Self-imposed deadline.** Recommend a hard deadline of 2 weeks from now (mid-June 2026) for arXiv submission, to avoid indefinite polish and to make the result public before any concurrent work appears.

---

**End of memo.**
