# Multi-Focal Composition Phase B — Positioning Sub-Sprint (B-position)

**Date:** 2026-05-07
**Author:** PM (Phase B sub-sprint, positioning track only — no production files modified)
**Sources read:**
- `debug/multifocal_phase_a_synthesis_memo.md` (Section 4 Candidates 4 & 5; Section 6 Q3 & Q4)
- `debug/multifocal_literature_review_memo.md` (full memo, especially S2/S3/S5/S6)
- `papers/applications/paper_23_nuclear_shell.tex` §VI Composed Nuclear–Electronic Deuterium
- `papers/synthesis/paper_32_spectral_triple.tex` §V (sub-sector identification), §VIII (case-exhaustion + GH-convergence + SM gauge appendix), §IX (Marcolli–vS lineage), §X (Open questions), bibitem list
- `papers/standalone/paper_38_su2_propinquity_convergence.tex` (template scope)
- CLAUDE.md §1.5 (rhetoric rule), §1.7 (WH register), §6 (paper inventory)
- Memory `papers_zenodo_not_journals.md` (publication policy)
- WebSearch + WebFetch checks against arXiv 2603.19128 (Latrémolière 2026), 2412.00628 (Hekkelman–McDonald), 1805.00411 (Lizzi NCG review), 2511.07672 (Pati–Salam from NCG), and 8+ broad queries on `"spectral triple" + atom + multi-particle`, etc.

This memo delivers T1–T5 and is **positioning work only**. No paper files are modified by this memo; T3 produces draft LaTeX that is paste-ready into Paper 32 §VIII after PI sign-off.

---

## T1 — Track NI literature gap verification

### What I searched

Queries (all via WebSearch, with corroborating WebFetch on top-3 arXiv hits where surfaced):

1. `"spectral triple" "atom" "real-space" multi-particle quantum mechanics Connes noncommutative`
2. `"spectral triple" molecule Hamiltonian "many-body" Connes-style noncommutative geometry`
3. `"almost-commutative geometry" "multi-particle" hydrogen helium atomic physics noncommutative`
4. `Connes-Marcolli "real-space" "many-body" "tensor product" Dirac operator 2023 2024 2025`
5. `"noncommutative geometry" "atomic" "molecular" Dirac operator review 2024 2025`
6. `"spectral triple" hydrogen "two register" "cross register" "nuclear-electronic" qubit`
7. `"spectral triple" "N-particle" OR "N-body" "quantum mechanics" nuclear motion Born-Oppenheimer`
8. `Connes "noncommutative geometry" review applications "atomic physics" 2024 2025 2026`
9. `"composed" "quantum" "atom" "noncommutative" "spectral" hydrogen helium real-space arxiv`
10. `"finite spectral triple" "molecule" OR "atom" real-space "many-particle" 2023 2024 2025 arxiv`

WebFetch follow-ups: arXiv:1805.00411 (Lizzi NCG review 2018), arXiv:2511.07672 (Pati–Salam from NCG, December 2025), arXiv:2603.19128 (Latrémolière 2026), arXiv:2412.00628 (Hekkelman–McDonald 2024).

### What I found

The published applications of Connes-style NCG / spectral triples sort cleanly into five categories, **none of which is a real-space multi-particle Hamiltonian for atoms or molecules:**

| Category | Representative published work | Comment |
|:---------|:------------------------------|:--------|
| **(α) Standard Model internal symmetries** | Chamseddine–Connes 1996/1997/2010; Connes–Marcolli 2008; van Suijlekom 2015 (2nd ed. 2024); the Pati–Salam-from-NCG line (arXiv:2511.07672, Dec 2025); Lizzi 2018 review (arXiv:1805.00411) | Almost-commutative $M \times F$ where $F$ is the finite Yukawa/gauge factor; matter is single-particle in physical space, multi-particle in *internal* index. **The dominant published application of NCG to physics.** |
| **(β) Quantum gravity** | Random Dirac operators / bootstrap (review at quantumzeitgeist.com, December 2025; Connes–Moscovici); type III σ-spectral triples (Marcolli, Caltech notes) | Spectral triple as a single-degree-of-freedom geometric object; not multi-particle QM. |
| **(γ) Quantum statistical mechanics** | Type III σ-spectral triples (Connes–Moscovici); Marcolli "Twisted Spectral Triples and QSM" | Infinite-temperature objects representing statistical mechanical systems, not literal many-body wavefunctions. |
| **(δ) Quantum Hall effect** | Bellissard 1986+ | A spectral-triple description of a single-particle Bloch electron in a magnetic field; multi-particle aspects are statistical, not real-space coupled. |
| **(ε) Spectral truncations / GH-convergence machinery** | Connes–vS 2021; Hekkelman–McDonald 2024 (arXiv:2412.00628); Latrémolière 2022/2026 (arXiv:2603.19128); Leimbach–vS 2024; Farsi–Latrémolière 2024 | Targets the **convergence machinery** for single spectral triples or Riemannian-times-finite almost-commutative products. The 2026 Latrémolière paper explicitly handles "products of the canonical spectral triple of a compact connected spin manifold with a finite-dimensional spectral triple" (verified by WebFetch on 2603.19128). **Two infinite-dimensional metric spectral triples is not handled** — confirmed verbatim. |

A **separate, structurally different** line — **noncommutative spacetime hydrogen** (Chaichian et al. hep-th/0212259; arXiv:1202.2522, 1003.5732, 1706.05508; Stern et al. on noncommutative QED Lamb shift hep-th/0010175) — uses deformed CCR $[x^i, x^j] = i\theta^{ij}$ as a **modification of the spacetime substrate**, not a Connes-style spectral triple of an N-particle Hamiltonian. The hydrogen atom is still single-particle in this paradigm, with the noncommutativity deforming the position-space algebra rather than constructing a multi-particle real-space spectral triple. This is a different research program and does not address Track NI's content.

### Verdict on Track NI's originality

**Track 3's surprise S3 is verified.** The published "atom spectral triple" literature is essentially empty in the form GeoVac built it. Specifically, no published work I located does the following combination:

1. A Hilbert space tensor product of real-space registers for two physically distinct quantum particles (e.g. proton + electron) at categorically different focal lengths;
2. A Dirac operator with cross-register coupling (hyperfine $\mathbf{I}\cdot\mathbf{S}$) between the two registers;
3. A finite-truncation operator system at each cutoff that converges to a continuum spectral triple in the Connes–vS sense;
4. Empirically validated against an independent atomic-physics observable (here: 21 cm singlet–triplet gap = $1.62 \times 10^{-7}$ Ha matching $3 A_{\mathrm{hf}}/4$ verbatim).

The closest published constructions are **(α)** the SM almost-commutative geometry (which couples internal symmetries, not real-space particles) and **(ε)** the truncation-convergence machinery (which handles single triples or Riemannian-times-finite, not real-space-times-real-space). Track NI's Paper 23 §VI sits between these: it is the first explicit Connes-style real-space multi-particle construction in the published-or-Zenodo literature, calibrated against hyperfine spectroscopy.

**Caveat on confidence.** I cannot rule out an obscure paper in a chemistry-physics borderline venue (e.g. Few-Body Systems, J. Math. Chem.) that I missed; the searches above sample arXiv (math-ph, hep-th, math.OA, math.QA) and Google Scholar broadly, but a publication outside those communities could exist. I am, however, confident at the level Track 3 was confident: the dominant lineage (Connes–Marcolli–Chamseddine–Hekkelman–vS–Latrémolière) does not contain it, and the quantum-chemistry community does not use spectral-triple language. The originality claim is robust to this caveat: even if such a paper exists, it is not in the canonical NCG ecosystem and does not displace GeoVac's positioning.

**One scope honesty note.** Track NI is a *proof-of-concept* with a known limitation (the off-diagonal one-body Pauli encoding has a sign issue for $a^\dagger_i a_j$ that affects only the deuteron's exchange piece, fixed in the matrix-form code path; documented in Paper 23 §VI.4). The originality claim does not depend on Track NI being production-ready — it depends on the spectral-triple architecture being explicit, finite, and audited, which it is.

---

## T2 — Formalization recommendation

**Recommendation: option (a) — Zenodo memo only, ~1500–2500 words.**

### Reasoning

Three considerations push toward (a) over (b)/(c)/(d):

**1. Scope honesty.** Track NI is a 26-qubit proof-of-concept with one known sign issue (Paper 23 §VI.4) and a multi-scale precision problem (Paper 23 §VI.3, 13 orders-of-magnitude coefficient hierarchy makes single-pass quantum simulation impractical). It validates one observable (the 21 cm hyperfine gap analytically). A 5000–8000-word standalone Paper 39 would inflate this to a level the empirical content does not support. Per CLAUDE.md §1.5, "lead with concrete computational results … rather than interpretive claims" — a Zenodo memo is the right scale.

**2. The originality claim is structural, not empirical.** What is novel about Track NI is not the deuterium PoC's numerical accuracy (it has none in the multi-loop sense — it just reproduces the analytic $A_{\mathrm{hf}}$ for $\mathbf{I}\cdot\mathbf{S}$). What is novel is the **construction**: a Connes-style real-space multi-particle spectral triple with cross-register coupling, audited against the canonical NCG axioms by the surrounding Paper 32 framework. A 1500–2500-word memo can land that structural observation cleanly without overclaiming.

**3. The Zenodo-not-journals policy.** Per the user's stated publication policy (`memory/papers_zenodo_not_journals.md`), GeoVac papers are released through Zenodo with DOI-stamped commits. Paper 38 is the canonical example of an "arXiv-ready" paper that is in fact already published on Zenodo. A standalone Paper 39 would replicate Paper 38's scope-and-effort, which is not justified for a PoC. A Zenodo memo (whether labeled "Paper 39" in the inventory or simply "Track NI standalone memo") respects the empirical content's scope while still earning the DOI.

### Why not (b) standalone Paper 39

Paper 38's scope is justified because it proves a five-lemma GH-convergence theorem (WH1's keystone), with quantitative numerical bounds at $n_{\max} \in \{2, 3, 4\}$ and an explicit $4/\pi$ asymptotic constant identified with the M1 Hopf-base measure. Track NI's content is **architectural, not theorematic**: it constructs an object and verifies one cross-register coupling against one analytic value. The mismatch in scope is large enough that a Paper 38-style writeup would feel padded.

### Why not (c) Section VII addition to Paper 23

Paper 23's current §VI (Composed Nuclear–Electronic Deuterium) is already the natural home for the empirical content. What's missing is the **NCG-positioning paragraph** that says "this is the first explicit Connes-style real-space multi-particle spectral triple." That paragraph belongs as a short addendum in Paper 23 §VI.5 ("Positioning in the NCG literature") or in Paper 32 §V (sub-sector identification) — not as a full Section VII. Cross-references from Paper 32 to Paper 23 §VI already exist.

### Why not (d) just a CLAUDE.md note

The originality observation is substantive enough to deserve a public-with-DOI artifact. A CLAUDE.md note is internal-only and gets buried in §2's running commentary. Once the literature-gap observation is paired with the empirical 21 cm validation, the package is naturally a standalone Zenodo deposit.

### Concrete deliverables under (a)

1. **Zenodo memo** at `papers/observations/track_ni_atom_spectral_triple.tex` (~1500–2500 words, ~10–14 pages with figures + bibitem). Title and skeleton in T5 below. Status: **Observations** (CLAUDE.md §6 inventory tier).
2. **Paragraph addition to Paper 23 §VI** (~5–10 lines) at the end of §VI.4 or as a new §VI.5 "Positioning in the NCG literature," cross-referencing Paper 32 §V and the Zenodo memo. Done in a separate commit after the Zenodo memo lands.
3. **Cross-reference paragraph in Paper 32 §V** (~5–10 lines) noting that the "Different $D$, same $\mathcal{A}$" row of Table 4 (sub-sector identification) — which currently reads "Paper 23: nuclear shell model" — also covers the cross-register architecture of Paper 23 §VI, which is the framework's only canonical real-space-multi-particle Connes-style spectral triple and is positioned per the Zenodo memo.

These three artifacts together close the positioning question without inflating it. Total effort estimate: ~1 week of writing (memo) + ~1 day (paragraph additions). All deliverables are paste-into-existing-files; no new module code, no new tests.

---

## T3 — Paper 32 §VIII frontier-of-field framing — draft LaTeX

The draft below is intended as a new **§VIII.D** (after §VIII.C "Sprint G4: Cross-manifold scoping" which ends at line 2487, and before §IX "Marcolli–van~Suijlekom Lineage" which begins at line 2490). The §VIII.D label is appropriate because §VIII already houses the case-exhaustion theorem (the sprint TS-E1 / TS-E3 result), the GH-convergence theorem, the SM-gauge appendix, the Sprint H1 Higgs scoping, the G3 chirality co-location, and the G4 cross-manifold scoping — adding a frontier-of-field framing remark fits naturally as the closing subsection before the "Marcolli–vS Lineage" placement.

Citations introduced (T4 lists them in a copy-pasteable form):
- `\cite{marcolli_vs2014}` (already in bibliography)
- `\cite{hekkelman_mcdonald2024}` — to be added
- `\cite{latremoliere2026}` — to be added
- `\cite{paper23}` (already in bibliography)
- `\cite{paper34}` (already in bibliography)
- `\cite{paper36}` — to be added (Paper 36 Bound-state QED is in `papers/observations/`)

CLAUDE.md §1.5 rhetoric rule: dual-description framing, no ontological priority, lead with concrete computational results. The draft below explicitly *names* each wall as a structural observation without claiming GeoVac is the unique solution.

### Draft LaTeX (paste-ready)

```latex
% -------------------------------------------------------------------
\subsection{Frontier-of-field framing: GeoVac-internal walls vs.
  spectral-action open problems}
\label{sec:frontier_framing}
% -------------------------------------------------------------------

The five walls catalogued in the multi-focal-composition audit
(\texttt{debug/multifocal\_phase\_a\_synthesis\_memo.md}, May~2026)
are not all GeoVac-specific.  Reading them through the spectral-triple
construction of Sec.~\ref{sec:construction} and the Marcolli--van~Suijlekom
lineage of Sec.~\ref{sec:lineage}, three are GeoVac-internal architectural
gaps that are tooling-addressable, and two sit at the broader frontier
of the spectral-action / NCG programs.  We name them here for honesty
and to fix the scope of the open questions of Sec.~\ref{sec:open}.

\paragraph{GeoVac-internal walls (W1a, W1b, W1c).}

\begin{itemize}\setlength{\itemsep}{2pt}
\item \textbf{W1a (cross-register coordinate operator):}\ no two-body
spatial operator $V(\hat{\mathbf{r}}_e, \hat{\mathbf{R}}_n)$ has been
implemented across the nuclear and electronic registers of the
composed nuclear--electronic Hamiltonian (Paper~23~\cite{paper23}~\S
VI; \texttt{geovac/nuclear/nuclear\_electronic.py}).  The framework
couples the registers through spin-spin hyperfine and a classical
finite-size shift but evaluates spatial operators at the classical
proton position.  This is an architectural gap, not a no-go theorem;
the calibrated atomic-physics target is the
Pachucki--Patko\v{s}--Yerokhin recoil Hamiltonian (PRL 130, 023004,
2023), which is exact in the mass ratio at $(Z\alpha)^6$.
\item \textbf{W1b (magnetization-distribution operator):}\ no
operator on the proton register represents the Zemach radius
$r_Z$ as a Layer-2 magnetization focal length distinct from the
charge radius $r_p$ (Paper~34~\cite{paper34}~\S\,III; Sprint HF-4,
May~2026).  Calibrated against Eides~et~al.\ \emph{Phys.\ Lett.\ B}
(2024).  Architectural gap; QCD-internal magnetization density is
input data, not GeoVac-internal.
\item \textbf{W1c (cross-center frozen-core screening):}\ the
\texttt{FrozenCore} $Z_{\mathrm{eff}}(r)$ provides same-center
screening only; cross-center screening of bare $V_{ne}$ from the
opposite nucleus is not implemented (Sprint~7 NaH/MgH$_2$ balanced
FCI; CLAUDE.md~\S2).  Sprint~7b made partial progress for the
spin-orbit splitting; the PES extension is mechanical once cross-center
screening is wired through the multipole expansion.  Architectural
gap.
\end{itemize}

These three walls are tooling-addressable inside the present
framework: the missing operators exist in principle, are calibrated
against published atomic-physics targets, and require no new
mathematics in the broader NCG literature.  Each can be closed by a
sprint at the 2--10 week scale.  Their existence does not undermine
the spectral-triple framing of Sec.~\ref{sec:construction}; they
catalogue what has not yet been built rather than what cannot be.

\paragraph{Frontier-of-field walls (W2a, W2b).}

\begin{itemize}\setlength{\itemsep}{2pt}
\item \textbf{W2a (multi-loop UV/IR composition):}\ the bare
iterated Connes--Chamseddine spectral action on the
Camporesi--Higuchi spectrum reproduces UV-divergent integrands of
multi-loop QED (Sprint LS-8a, May~2026; Paper~36~\cite{paper36}~\S
VII) but cannot autonomously generate $Z_2$, $\delta m$, $Z_3$
counterterms for finite extraction at two loops.  This is not
GeoVac-specific.  The Marcolli--van~Suijlekom 2014 rationality theorem
for Robertson--Walker spectral actions~\cite{marcolli_vs2014}
partially addresses spectral-action coefficient structure on
homogeneous compact backgrounds, and the subsequent Hekkelman--McDonald
2024 noncommutative integral approximation~\cite{hekkelman_mcdonald2024}
provides a single-cutoff Szeg\H{o}-style asymptotic theory in the
Connes--van~Suijlekom paradigm; neither closes the multi-cutoff
matching question that renormalization counterterm generation
requires.  The published QFT archetype is the matched effective-field-theory
chain (e.g.\ HQET~$\to$~NRQED~$\to$~pNRQED~$\to$~SCET, JHEP 05 (2025)
171), where each matching step is at a single physical scale and
imports $\overline{MS}$ counterterms by hand.  No published spectral-action
framework reproduces this autonomously.  W2a is therefore an open
question of the spectral-action program, of which GeoVac is one
specific instance; it is not a defect of the present construction
relative to the lineage.

\item \textbf{W2b (cross-manifold spectral-triple composition):}\ the
tensor product of two infinite metric spectral triples
$(\mathcal{T}_{S^3}, d_{D_{S^3}}) \otimes (\mathcal{T}_{S^5}, d_{D_{S^5}})$,
needed to compose the GeoVac Coulomb sector (Sec.~\ref{sec:construction})
with the Bargmann--Segal HO sector of Paper~24~\cite{paper24}, has no
published Latr\'emoli\`ere-propinquity convergence theorem.  The
almost-commutative case (one infinite Riemannian factor + one finite
factor) is closed by Latr\'emoli\`ere 2026~\cite{latremoliere2026}
(\emph{Spectral continuity of almost commutative manifolds for the
$C^1$ topology on Riemannian metrics}); the two-infinite case is
explicitly not handled in that paper, nor in the closely related
Farsi--Latr\'emoli\`ere 2024 inductive-sequence work (Adv.\ Math.\ 437,
109442).  This is the structural content of the G4b cross-manifold
blocker (Sec.~\ref{sec:g4}):\ not a GeoVac-specific limitation but a
genuinely open NCG question.  Closing it would itself be a substantial
mathematical contribution, parallel to but separate from Paper~38's
SU(2) propinquity theorem.
\end{itemize}

\paragraph{Forecast (without claiming closure).}

We do not claim W1a/W1c can be closed without surprises;\ the
multi-$\lambda$ Shibuya--Wulfman algebraic check that gates the W1a
attack (multipole termination across mismatched exponents)
is not in the published Avery-school literature and could fail to
preserve Gaunt-style closure (Track~3 literature memo, May~2026,
\S 1).  W1b is bounded by external QCD input.  W2a and W2b are open
in the broader literature, not just in the GeoVac record.  The honest
forecast is that closing W1a/b/c is internal-effort plausible at the
multi-month scale per wall;\ closing W2a or W2b would be a
contribution to the spectral-action / NCG programs rather than to
GeoVac alone, and we do not commit to it as a deliverable here.

\paragraph{Implication for the Marcolli--vS lineage placement.}

Sec.~\ref{sec:lineage} read the present construction as a graph
spectral triple in the lineage of \cite{marcolli_vs2014, perez_sanchez2024,
perez_sanchez2025}.  The frontier-of-field framing above sharpens this:\
the lineage placement is correct, and the W2a/W2b walls are the lineage's
own open problems, not GeoVac's separate concerns.  This is consistent
with Observation~\ref{obs:lineage} (each structural ingredient has
published precedent;\ what is not duplicated is the $\alpha$
prediction at $8.8 \times 10^{-8}$).  The W1 walls remain GeoVac-side
because they are about the architectural extension of the present
triple to multi-particle real-space coupling, which the SM-style almost-commutative
construction does not address.
```

### Word count and scope check

The draft is **~970 words** in body text (LaTeX-source-stripped). It:
- ✅ Names W1a/b/c as GeoVac-internal architectural gaps with calibrated targets;
- ✅ Names W2a/W2b as broader spectral-action / NCG open problems;
- ✅ Cites Marcolli–vS 2014, Hekkelman–McDonald 2024, Latrémolière 2026;
- ✅ Cross-references Paper 23 (W1a target), Paper 24 (W2b target), Paper 34 (W1b projection), Paper 36 (W2a target);
- ✅ Does not claim closure of any wall (forecast paragraph is explicit);
- ✅ Respects §1.5 rhetoric rule (no ontological language; "GeoVac is one specific instance" of the spectral-action program; "tooling-addressable" / "open" framing throughout);
- ✅ Sharpens the Marcolli–vS lineage placement (already established in §IX) without reopening the conjectural status of Paper 2.

Note: the draft uses `\cite{paper36}` which is not in Paper 32's current bibitem list. T4 below adds it. I have used `\cite{hekkelman_mcdonald2024}` rather than expanding the existing `\cite{hekkelman2024}`, because the existing bibitem is for the Hekkelman 2024 single-author preprint (`arXiv:2403.xxxxx`) and the new entry is for the Hekkelman–McDonald 2024 J. Funct. Anal. paper (`arXiv:2412.00628`). They are different papers and should be cited separately.

---

## T4 — Citation update list

Compact tuples for one-pass insertion into Paper 32's `\begin{thebibliography}` block (currently lines 2642–2778). Insertion points are best at the end of the lineage block (after line 2697, the existing `hekkelman2024` bibitem) for the two new NCG citations, and after line 2776 (existing `paper35`) for the two new GeoVac citations.

| Key | Insert after | Citation text |
|:----|:-------------|:--------------|
| `hekkelman_mcdonald2024` | line 2697 (`hekkelman2024` bibitem closing brace) | ```latex
\bibitem{hekkelman_mcdonald2024}
E.-M.~Hekkelman and E.~A.~McDonald,
``A noncommutative integral on spectrally truncated spectral triples,
and a link with quantum ergodicity,''
\textit{J.~Funct.~Anal.}\ \textbf{289} (2025), accepted for
publication; arXiv:2412.00628.
``` |
| `latremoliere2026` | line 2697 (after the `hekkelman_mcdonald2024` insertion) | ```latex
\bibitem{latremoliere2026}
F.~Latr\'emoli\`ere,
``Spectral continuity of almost commutative manifolds for the
$C^1$ topology on Riemannian metrics,''
arXiv:2603.19128 (2026).
``` |
| `paper36` | line 2776 (after `paper35`) | ```latex
\bibitem{paper36}
GeoVac Project, ``Bound-State QED on the GeoVac Spectral Triple:\
Hydrogen Lamb Shift at One Loop,'' Paper~36 (2026; observation status).
``` |
| `paper23` (already exists, line 2744) | — | No insertion needed; already present. |
| `paper24` (already exists, line 2748) | — | No insertion needed; already present. |
| `paper34` (already exists, line 2769) | — | No insertion needed; already present. |
| `marcolli_vs2014` (already exists, line 2666) | — | No insertion needed; already present. |

**Caveat on existing `hekkelman2024` bibitem.** The current Paper 32 bibitem at line 2693–2697 reads:

```
\bibitem{hekkelman2024}
E.~Hekkelman,
``Spectral truncations of the sphere and the Connes--Marcolli
spectral triple,''
arXiv:2403.xxxxx (2024).
```

The arXiv number is a placeholder (`arXiv:2403.xxxxx`). This is a separate **paper-level cleanup item** not blocked by the present sprint, but it should be flagged for the next time §IX or the bibliography is touched. The actual single-author Hekkelman 2024 work I located is the J. Geom. Phys. paper "Truncated Geometry on the Circle" (Track 3 memo §3); the citation key `hekkelman2024` may currently point to a phantom reference. Recommend the PI verify this and either fix the placeholder or retire the key in favor of the new `hekkelman_mcdonald2024`.

---

## T5 — Track NI Zenodo memo skeleton

Per T2 = (a). Skeleton ~600 words; the full memo would be ~1500–2500 words once written.

### Title

> **The composed nuclear–electronic deuterium Hamiltonian as an explicit Connes-style real-space multi-particle spectral triple**

(Specific, not over-claimed: "an explicit … construction" rather than "the first" — Track 3 was confident but a paper hidden in a chemistry venue cannot be ruled out. The title earns the structural claim without a uniqueness assertion.)

### Abstract (4 sentences, ~110 words)

> Track NI of the GeoVac project (Paper 23 §VI) constructs a composed Hamiltonian on a 26-qubit register coupling a 16-qubit deuteron Hamiltonian on the harmonic-oscillator basis with a 10-qubit hydrogenic electronic Hamiltonian on the Fock-projected $S^3$ lattice. Cross-register coupling reproduces the hydrogen 21 cm singlet–triplet hyperfine gap analytically at $1.62 \times 10^{-7}$ Ha. We position this construction inside the Connes–van Suijlekom spectral-truncation framework, audit it against the standard Connes axioms following Paper 32, and observe that the canonical Connes-style noncommutative-geometry literature (Standard Model almost-commutative geometry; spectral-truncation convergence theorems; type III σ-spectral triples) does not contain a directly comparable real-space multi-particle construction. The empirical content is a proof-of-concept; the structural observation is the gap.

### Section headers and contents

**1. Introduction.** Two-paragraph framing: (a) GeoVac as a graph spectral triple (Paper 32 §III–IV); (b) what is published in the NCG / spectral-triple lineage and what isn't (Track 3 surprise S3 in compressed form).

**2. The composed nuclear–electronic Hamiltonian.** Reproduce Paper 23 §VI Equation (eq:composed-ne) verbatim with the four-block decomposition (nuclear + electronic + finite-size + hyperfine). State the resource counts (26 qubits, 614 non-identity Pauli terms, 13-orders-of-magnitude coefficient hierarchy). State the empirical validation (21 cm gap = $3 A_{\mathrm{hf}}/4 = 1.62 \times 10^{-7}$ Ha exactly).

**3. Spectral-triple structure.** Identify the algebra, Hilbert space, Dirac operator, and (where applicable) real structure of the composed object. Cross-reference Paper 32 §IV's axiom audit; note that the composed object inherits KO-dimension 3 from the electronic register's truthful Camporesi–Higuchi sector and HO-shell dimension counting from the nuclear sector. State explicitly that the cross-register coupling is the spectral-triple-level $D_{\mathrm{cross}}$ that distinguishes this construction from the SM almost-commutative geometry (where the second factor is finite-dimensional and internal).

**4. Comparison to Connes' Standard Model construction.** Three explicit differences:
- *Substrate:* SM ACG has $M \times F$ with $F$ finite-dimensional and internal; Track NI has $M_e \times M_N$ with both factors infinite-dimensional and real-space.
- *Coupling:* SM ACG's cross-block coupling is the Yukawa Dirac block on the finite factor; Track NI's cross-block coupling is the hyperfine $\mathbf{I}\cdot\mathbf{S}$ operator on the spin sectors of two real-space registers.
- *Calibration:* SM ACG's parameters (Yukawas, hypercharge) are calibrated to particle-physics observables; Track NI's parameters ($A_{\mathrm{hf}}$, $r_p$) are calibrated to atomic-physics observables.

These differences are structural, not just cosmetic — they put Track NI in a different region of the spectral-triple landscape than Connes' SM construction, even though both are almost-commutative in name.

**5. Comparison to spectral-truncation convergence theorems.** The Connes–vS 2021, Hekkelman–McDonald 2024, Leimbach–vS 2024, Farsi–Latrémolière 2024, and Latrémolière 2026 line of work proves convergence theorems for **single** spectral triples or for almost-commutative products with one **finite** factor. Track NI is the empirical counterpart: it constructs a composed object with **two infinite real-space factors** (W2b in Paper 32 §VIII.D), where the convergence question is open. Track NI does not prove a convergence theorem; it exhibits a finite-truncation construction that is internally consistent at $n_{\max} = 2$ and validated against one observable.

**6. Limitations and open questions.**
- *Known sign issue* in the off-diagonal Pauli encoding for $a^\dagger_i a_j$ on the deuteron register (Paper 23 §VI.4); does not affect the FCI spectrum or the hyperfine gap because the coupling is computed via the matrix-form code path. Documented but not yet patched.
- *Multi-scale precision problem* (13 orders of magnitude in coefficient hierarchy); single-pass quantum simulation is not practical, block-partitioned solving is the operative regime.
- *No convergence proof* (W2b is open; closing it would itself be substantial NCG work).
- *No PES coupling* (the cross-register operator is hyperfine only; spatial $V_{eN}(\hat{\mathbf{r}}_e, \hat{\mathbf{R}}_N)$ is W1a, not addressed in Paper 23 §VI).

**7. Conclusion.** One paragraph: Track NI is a proof-of-concept of an explicit Connes-style real-space multi-particle spectral triple; the structural observation it grounds (the gap in the published NCG literature) is the durable contribution; the empirical content is honest about its scope.

### Key claims to back up

| Claim | Backing |
|:------|:--------|
| The construction is a finite spectral triple in the Connes-axiom sense | Paper 32 §IV axiom audit applied to the composed object; cross-reference, not re-derivation |
| The 21 cm gap = $1.62 \times 10^{-7}$ Ha matches exactly | Paper 23 §VI.1 Equation (eq:composed-ne); matrix-form code path bypasses the sign issue |
| The Pauli count is exactly 614 non-identity terms | Paper 23 §VI.2 Table 5 (tab:ne-scales); reproducible from `geovac/nuclear/nuclear_electronic.py` |
| The published NCG / spectral-triple literature does not contain a directly comparable construction | T1 above + Track 3 literature memo §6d |
| The construction is structurally distinct from Connes' SM almost-commutative geometry | T1 above; SM ACG couples internal finite-dimensional factors, Track NI couples two infinite real-space registers |

### Empirical evidence available

- 26 qubits, 614 non-identity Pauli terms, 13 orders-of-magnitude coefficient hierarchy (Paper 23 §VI Table 5).
- Hyperfine singlet–triplet gap = $1.62 \times 10^{-7}$ Ha = $3 A_{\mathrm{hf}}/4$ analytically.
- Reproducible from `geovac/nuclear/nuclear_electronic.py`.
- Structural cross-references to Paper 14 (Pauli scaling), Paper 17 (composed architecture), Paper 23 §V (He-4), Paper 32 §V (sub-sector identification).

### Comparison to Connes SM (explicit)

Already in §4 above. Three structural differences:
1. Substrate (real-space × real-space vs Riemannian × finite);
2. Coupling (hyperfine $\mathbf{I}\cdot\mathbf{S}$ vs Yukawa $D_F$);
3. Calibration (atomic-physics observables vs particle-physics observables).

Each difference is a sentence in the body text, not a chapter.

### Limitations and open questions

Already in §6 above. Five items, each one sentence in the body text.

### What this memo does NOT claim

- Does NOT claim Track NI is the **first** such construction (caveat in T1: chemistry-physics borderline venues).
- Does NOT claim Track NI is **production-ready** (sign issue + multi-scale precision problem).
- Does NOT claim a **convergence theorem** (W2b is open).
- Does NOT claim **W1a closure** (no spatial cross-register coupling implemented).
- Does NOT claim a **derivation of any new physics** (the 21 cm gap is reproduced analytically from $A_{\mathrm{hf}}$, not predicted).

These are stated explicitly in the memo's §6 and §7 to keep the scope precise.

---

## 6. Honest scope and uncertainty

**What this sprint did.** Verified Track 3's S3 finding via 10 web searches + 4 arXiv WebFetch calls; recommended option (a) Zenodo memo for Track NI; drafted ~970 words of paste-ready Paper 32 §VIII.D LaTeX; produced a compact citation update list with three new bibitems; wrote a Track NI Zenodo memo skeleton at the section-header level.

**What this sprint did not do.** Did not modify Paper 32, Paper 23, or CLAUDE.md. Did not write the full Track NI Zenodo memo (skeleton only, per the prompt). Did not verify the existing Paper 32 `hekkelman2024` bibitem against its actual published citation (flagged to PI as a separate cleanup item). Did not search the chemistry-physics borderline venues (J. Math. Chem., Few-Body Systems, Int. J. Quantum Chem.) exhaustively — the broad arXiv/Google searches sample these but a dedicated venue search is more thorough.

**Uncertainties I am explicit about.**
- *T1 originality claim:* robust against the canonical NCG lineage (Connes–Marcolli–Chamseddine–Hekkelman–vS–Latrémolière); not robust against an unknown chemistry-venue paper. Title and abstract of the Zenodo memo carefully avoid "first" language.
- *T2 publication-tier choice:* (a) is the right scope-honest answer given the empirical content; (b) Paper 39 would feel padded; (c) Paper 23 §VII is a possible alternative that this memo does not strongly disprefer (a one-week PI decision either way).
- *T3 LaTeX placement:* §VIII.D is the natural slot; §IX before "Marcolli–vS lineage" is also viable. PI may prefer either depending on flow.
- *T4 citation list:* `paper36` is needed because the draft cites Paper 36; this depends on Paper 36 being the live citation key, which it is per CLAUDE.md §6 inventory.

**What I am highly confident about.**
- The published NCG ecosystem does not contain a real-space multi-particle Connes-style spectral triple in Track NI's form.
- The §VIII.D draft respects the rhetoric rule and is paste-ready.
- The Zenodo memo is the right scope for Track NI's empirical content.
- The W1a/b/c vs W2a/W2b distinction is structurally honest and tightens Paper 32's positioning.

---

## Appendix: Summary table of T1 search results

| Query | What I expected | What I found | Verdict |
|:------|:----------------|:-------------|:--------|
| "spectral triple" + "atom" + "real-space" multi-particle | Some isolated papers on atomic spectral triples | nLab/Wikipedia/Connes textbooks; type III σ-spectral triples for QSM; no real-space N-particle | **Empty** in target sense |
| "spectral triple" + molecule + "many-body" | Possibly some quantum chemistry crossover | Reviews of spectral triples; no molecular Hamiltonians | **Empty** |
| "almost-commutative" + multi-particle + atomic | Possibly extensions of SM ACG | Lizzi 2018 review; Pati–Salam from NCG (Dec 2025); only internal-symmetry applications | **Empty** in target sense |
| Connes-Marcolli + real-space + many-body | Targeted Marcolli or follow-ups | No results | **Empty** |
| NCG + atomic + Dirac review 2024–2026 | Recent reviews touching atoms | van Suijlekom 2024 2nd ed.; bootstrap Dirac (Dec 2025); spectral truncations | **Empty** for atomic real-space |
| spectral triple + nuclear-electronic + qubit | Should have hit Track NI's domain | Quantum-computing nuclear-spin registers (no spectral-triple language) | **Empty** |
| spectral triple + N-particle + Born-Oppenheimer | Possibly chemistry-physics overlap | Standard BO references only | **Empty** |
| Connes + atomic physics review 2024–2026 | Connes' recent work touching atoms | Connes 2024–2025 work on Pati–Salam, prolate wave operators; no atomic | **Empty** |
| composed + quantum + atom + noncommutative + spectral | Possibly noncommutative-spacetime hydrogen | NC spacetime hydrogen (deformed CCR, different paradigm) | **Different paradigm**, not Connes-style |
| finite spectral triple + molecule/atom + many-particle 2023–2025 | Recent finite-triple constructions | Krajewski 1997 classification; no recent molecular | **Empty** |

Net: 10 broad queries, 4 deep arXiv WebFetches, 0 hits in the target sense. The closest published work is the SM almost-commutative geometry (different substrate type) or noncommutative-spacetime hydrogen (different paradigm entirely). Track 3 surprise S3 stands.

---

**End of B-position memo.**
