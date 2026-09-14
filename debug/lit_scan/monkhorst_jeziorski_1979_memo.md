# Monkhorst & Jeziorski 1979 — the "no linear dependence" claim, read

**Date:** 2026-09-12
**Trigger:** `memory/avery_method_and_prior_art_gaps.md` §8 — points 2 and 3
SUSPENDED pending this read. Surfaced as by-catch by
`debug/lit_scan/primaries_sw1965_toeplitz_pencil_memo.md` §1g.

## VERDICT: **SECONDARY-QUOTED**

The **abstract is REACHED and verified on two independent indexes**
(Crossref JATS + OpenAlex inverted index — byte-identical reconstructions),
together with the complete bibliographic record and the reference list. The
**two-page body was NOT reached**: the article is closed access with *zero*
repository copies anywhere (Unpaywall `oa_status: closed`, `oa_locations: []`;
OpenAlex `any_repository_has_fulltext: false`; Semantic Scholar
`openAccessPdf.status: CLOSED`). AIP is behind Cloudflare and 403s every route.

The **mathematical content of the method** is therefore established from the
lineage rather than from M–J's own page: an OA paper in the same construction
(arXiv:quant-ph/0403195) writes the identical equations, and the direct 1982
response note attributes the method by name. Everything below is tagged
VERIFIED / LINEAGE / INFERRED so the legs can be separated.

---

## 0. Bibliographic data (VERIFIED)

| Field | Value |
|:--|:--|
| **DOI** | **`10.1063/1.438337`** |
| Title | "No linear dependence or many-center integral problems in momentum space quantum chemistry" |
| Authors | Hendrik J. Monkhorst; Bogumił Jeziorski |
| Affiliation | Department of Physics, University of Utah, Salt Lake City, Utah 84112 (both) |
| Journal | *J. Chem. Phys.* **71**(12), 5268–5269 |
| Published | 15 December 1979 |
| Length | **2 pages** — a Note/Communication, not a full paper |
| References | 11 |
| Cited by | 53 (Crossref) / 54 (S2) / 55 (OpenAlex) |
| Access | Closed; no OA copy exists |

### DOI CORRECTION — the task's assumed DOI was WRONG

`10.1063/1.438305` is **a different paper**: Brescansin, Leite & Ferreira,
"A study of the ground states and ionization energies of H2, C2, N2, F2, and CO
molecules by the variational cellular method," *J. Chem. Phys.* **71**,
4923–4930 (1979). Verified against Crossref. **Do not propagate 438305.**

### Two further bibliographic traps

1. **Crossref and Semantic Scholar both misspell the second author
   "Jeriorski"** (and S2 additionally splits "Bogumil" and "Jeriorski" into two
   separate author records). The correct name is **Jeziorski** — confirmed on
   Monkhorst's own CV. A citation-gate regex keyed on the Crossref string will
   produce a wrong bibliography entry.
2. **Monkhorst's own CV gives the title with "and", not "or"**: "No Linear
   Dependence *and* Multi-Center Integral Problems in Momentum Space Quantum
   Chemistry" (item 37, `people.clas.ufl.edu/henk/cv-long/`). The journal's
   registered title uses "or" and "many-center". **Cite the journal form.**

### Reference list (VERIFIED via Crossref; 6 of 11 resolvable)

| # | Resolved |
|:--|:--|
| 1 | Monkhorst & Harris, *Int. J. Quantum Chem.* **6**, 601 (1972) — Fourier transform of two-center Slater orbital products |
| 2 | **Löwdin, "On the Nonorthogonality Problem", *Adv. Quantum Chem.* **5**, 185 (1970)** |
| 3 | **Fock, "Zur Theorie des Wasserstoffatoms", *Z. Phys.* **98**, 145 (1935)** |
| 5 | **Shibuya & Wulfman, *Proc. R. Soc. London* A **286**, 376 (1965)** |
| 7 | Biedenharn, "Wigner Coefficients for the R4 Group and Some Applications", *J. Math. Phys.* **2**, 433 (1961) |
| 9 | Löwdin, "Studies in Perturbation Theory. X.", *Phys. Rev.* **139**, A357 (1965) |
| 10 | *Adv. Phys.* **5**, 111 (1956) — unresolved, probably Löwdin |
| 4, 6, 8, 11 | Not resolvable (Crossref holds no metadata) |

---

## 1. What do they actually claim, precisely?

### The operative sentences (VERIFIED — full abstract, verbatim)

> "It is shown that the position space LCAO-type wave functions for a
> many-center **one electron** system can be easily obtained via a suitable
> momentum space approach. These wave functions and corresponding
> **variational approximations to the eigenvalues are calculated via
> diagonalizations of simple overlap matrices**. In the proposed approach the
> problems of many-center integrals and **instabilities due to
> overcompleteness of basis sets do not appear at all**."

Reconstructed independently from Crossref JATS and from OpenAlex's
`abstract_inverted_index`; the two agree word for word.

### Answer

It is **not a theorem**. It is a **property of a particular posing**, asserted
for a **particular problem class**, in a two-page note. Precisely:

- **Hypothesis 1 — one electron.** "a many-center **one electron** system."
  The claim is made for the one-electron multi-centre Coulomb problem (H2+ and
  relatives). It says nothing about N-electron systems, and nothing about
  two-electron integrals.
- **Hypothesis 2 — isoenergetic / Sturmian, shared scale.** The basis is the
  Fock-projected momentum-space set: all orbitals carry the **common exponent
  `|p0|` fixed by the energy**, not per-orbital exponents. (LINEAGE: this is
  what Fock ref 3 + Shibuya–Wulfman ref 5 + Biedenharn's R4 Wigner coefficients
  ref 7 are doing in an 11-item list.)
- **Hypothesis 3 — the object diagonalised is an overlap matrix.**
  "diagonalizations of **simple overlap matrices**" — an *ordinary* eigenvalue
  problem of the Shibuya–Wulfman-type matrix, not a generalized one.
- **Scope of "no linear dependence."** It is a claim about **absence of
  numerical instability from an overcomplete basis**, i.e. the classical
  Löwdin nonorthogonality problem (ref 2 is *literally* Löwdin's "On the
  Nonorthogonality Problem"). The framing is: the disease that afflicts
  position-space LCAO does not afflict this posing.

**The title is an advertisement for a re-posing, not a spectral theorem.**

---

## 2. Is it prior art for a conditioning / linear-dependence analysis?

**It is prior art for the TOPIC. It is not prior art for the ANALYSIS.**
These must be kept apart, and the corpus has been conflating them.

**Prior art for the topic — YES, unambiguously.** A 1979 *J. Chem. Phys.*
paper whose title is "No linear dependence ... in momentum space quantum
chemistry", citing Löwdin's nonorthogonality paper, is direct, explicit,
published engagement with linear dependence in exactly the momentum-space
Sturmian setting the corpus works in. Any corpus sentence implying that the
momentum-space / Sturmian community never noticed or never addressed linear
dependence is **false and must be retracted**.

**Prior art for the analysis — NO, on the evidence available.**
There is no indication anywhere in the reachable record that they compute or
bound: a condition number, an overlap-matrix spectrum, a smallest eigenvalue,
a singular-value distribution, or any growth law with basis size. The abstract
asserts *absence* ("do not appear at all") — a structural claim, with no
quantity attached. Supporting (not conclusive) evidence:

- It is **two pages**. A conditioning scaling law with a fitted or derived
  exponent does not fit alongside a method statement and a demonstration.
- The reference list contains **no numerical-linear-algebra source at all** —
  no Wilkinson, no Golub, no Szegő/Toeplitz, nothing on matrix conditioning.
  Every one of the 6 resolvable refs is quantum chemistry or group theory.
- The 55-item citing set contains no paper that cites it *for* a conditioning
  result; the citations are for the method (Koga's "One-electron diatomics in
  momentum space" series I–V, Avery's Sturmian papers, Aquilanti's group).

**HONEST CAP:** this leg is an argument from the abstract, the reference list
and the page count, **not from a read of the body**. It is the strongest form
available given that no copy of the body exists in any open repository. If the
body is ever obtained and contains a spectral statement, point 2 of the memory
file falls with it.

**So: the corpus's "no prior art for secular-matrix NORM GROWTH" (point 2)
SURVIVES. The corpus's "conditioning is a BLIND SPOT in the whole Avery canon"
(point 3) DOES NOT survive as phrased and must be rewritten.**

---

## 3. THE RESOLUTION — how both can be true

**Both claims are true, because they are claims about different objects: M–J's
is about the *posing*, the corpus's is about the *matrix*.**

The two formulations are the *same pencil*. The arXiv lineage paper
(Elesin, Podlivaev & Openov, quant-ph/0403195, following Shibuya–Wulfman, Fock,
Koga and Avery–Hansen) writes the many-centre one-electron momentum-space
problem as the homogeneous system (its Eqs. 5–7)

```
  sum_j H_ij a_j = 0 ,
  H_ij = |p0| S_NN'(R_k - R_k')  -  sum_k'' Z_k'' sum_N'' (1/n'') S_NN''(R_k - R_k'') S_N''N'(R_k'' - R_k')
  S_NN'(R) = integral dOmega  Y_N(Omega) Y_N'(Omega) exp[i p.R]
```

That `S` is **exactly** GeoVac's Shibuya–Wulfman matrix — the matrix of the
translation phase `exp(i p.R)` on the Fock sphere, which is the corpus's own
description of it — and `H = |p0| S - W` is, up to sign convention, exactly
Paper 60's `eq:sw`, `[W - k S] C = 0`. **The near-degeneracy is in their
matrix too.** `sigma_max -> 1` as the basis grows is a property of the symbol
`j_0(kR cot(chi/2))` and it does not care who is looking at it.

The difference is entirely in **how the eigenvalue is extracted**:

- **M–J / Novosadov / Koga:** `S` appears **only multiplicatively**. One solves
  the scalar nonlinear root problem `det H(|p0|) = 0` (the lineage paper is
  explicit: "Eq. (5) is the nonvariational equation ... one should solve the
  nonlinear equation det(H_ij) = 0 for |p0|"), or diagonalises the overlap
  matrix directly as M–J's abstract says. **Nothing is ever inverted.** A
  near-null direction of `S` produces a near-null direction of `H`, which
  contributes a spurious branch far from the physical root — it does not
  amplify error in the physical root. Redundant basis functions are inert.
- **GeoVac:** needs `S^{-1/2} W S^{-1/2}` to hand a *standard* Hermitian
  operator to a quantum block-encoding. Paper 60 says so in one sentence:
  cond(S) "multiplies the block-encoding cost: inverting or square-rooting S by
  quantum singular-value transformation costs a polynomial degree growing with
  kappa." **That is the only place the conditioning becomes a cost.**

Two further reinforcements of M–J's side, both independently corroborated by
the corpus's own measurements:

- **The basis is exactly orthonormal in the metric the problem is actually
  posed in.** `S_NN'(0) = delta_NN'` (lineage paper, explicit). The corpus
  measured the same thing: "the intra-center block of S is *exactly* the
  identity ... diagonal {1,1,1} with off-diagonal ~1e-16." The L2 Gram matrix —
  whose near-singularity *is* the classical linear-dependence disaster, and
  which the corpus measures at cond 5.8 and growing — **never enters the
  equation at all**. Avery & Hansen (IJQC 60, 201, 1996) put it functionally:
  the Sturmian set "is a basis of a **Sobolev space rather than a Hilbert
  space**."
- **"No many-centre integrals" is a resolution of identity, not an integral
  table.** The three-centre one-electron integrals never have to be evaluated
  because the `sum_N''` in `H_ij` — a completeness sum over intermediate
  hyperspherical states on the Fock sphere, weighted `1/n''` — reproduces them
  as **sums of products of two-centre SW integrals**. Exact in an infinite
  basis; truncation error only.

**Is it a lever?** Partially, and honestly: **a real trade, not a free win.**
Re-posing GeoVac's molecular problem determinantally would remove the
conditioning multiplier, but it buys that by reintroducing the **outer
nonlinear energy search** — which is precisely the cost Paper 60 counts as
*avoided* in the atomic isoenergetic case ("The outer energy search that would
otherwise multiply the cost does not arise, because the eigenvalues are the
energies"). M–J do not make the degeneracy go away; **they keep it out of the
denominator.** Whether the trade is favourable on a quantum device is an
unasked, answerable resource question and is the one genuinely new lever this
read produces.

**Candidate explanations from the task prompt, adjudicated:** different basis —
NO (same Coulomb-Sturmian shared-scale set, same Fock projection). Avoids the
`p = 0` degeneracy — NO (the same symbol, the same `sigma_max -> 1`). Functions
vs Gram matrix — **PARTLY YES**, and this is the sharp half: their claim is
about L2 linear dependence, which is *structurally absent* from a problem posed
in the `V_0`/Sobolev metric. About many-centre integrals being evaluable rather
than about conditioning — **PARTLY YES**, that is half the title and it is a
completeness-sum statement. Not shared-scale — NO, it is shared-scale.

---

## 4. Shibuya–Wulfman, Coulomb Sturmians, Fock projection?

**YES to all three, explicitly, in an 11-item reference list.**

- **Fock 1935** — ref 3, `Z. Phys.` **98**, 145, "Zur Theorie des
  Wasserstoffatoms". The Fock projection, cited directly. (VERIFIED)
- **Shibuya–Wulfman 1965** — ref 5, `Proc. R. Soc. London` A **286**, 376.
  The SW paper, cited directly. (VERIFIED)
- **Coulomb Sturmians** — not nameable from the abstract, but the shared-`|p0|`
  construction is what "diagonalizations of simple overlap matrices" on a
  Fock+SW foundation means, and the 1982 response note (below) classifies the
  method as **Sturmian** in its title. (LINEAGE)
- **Biedenharn's R4 Wigner coefficients** (ref 7) is the SO(4) recoupling
  machinery on the Fock sphere — the same object as the corpus's hyperspherical
  angular couplings.

**The corpus's `eq:general_v0` identification is independently corroborated.**
The lineage paper's Eq. (8) writes the SW integral as
`S_NN'(R) = (n/|p0|) integral dr phi_N(r - R_k) phi_N'(r - R_k') / |r - R_k|`,
i.e. **as a Coulomb-weighted overlap**. That is Paper 60's "the molecular
metric ... *is* `V_0`". Note the provenance line in Paper 60 stays correct as
written: the *integrals* are the canon's, the *identification of the matrix as
the `V_0`-weighted overlap* is the corpus's formalization — the canon writes
the formula without drawing the conclusion.

### Attribution correction — the method is NOT originally M–J's

**Duchon, Dumont-Lepage & Gazeau,** "On two Sturmian alternatives to the LCAO
method for a many-center one-electron system," *J. Chem. Phys.* **76**, 445–447
(1982), DOI `10.1063/1.442741` — a direct response note in the same journal —
states that the technique "was originally used by **Novosadov** and recently
revisited by **Monkhorst and Jeziorski**." (VERIFIED from the published
abstract.) So M–J are a *revival*, and **B. K. Novosadov** is the origin. Note
also that the 1986 citing paper "Use of overcomplete basis sets in
quantum-chemical calculations" (*J. Mol. Struct. THEOCHEM* **136**, 387) is by
**Gribov & Novosadov** — the same Novosadov, still on overcompleteness seven
years later. **There is a Russian-language overcompleteness thread here the
corpus has never touched.**

### A related Monkhorst paper the corpus half-knows

**Aissing & Monkhorst, "Linear Dependence in Basis Set Calculations for
Extended Systems," *Int. J. Quantum Chem.* **43**, 733 (1992)** — already noted
in `debug/lit_scan/{contraction_seam_e3,toeplitz_finite_section}_memo.md` and
correctly judged not to transfer. Worth recording that it is **the same
Monkhorst**, thirteen years later, on the same disease in a different setting.

---

## 5. Routes tried, with failure modes

| Route | Result |
|:--|:--|
| Crossref DOI `10.1063/1.438305` | **Wrong paper** (Brescansin et al., pp. 4923–4930) |
| Crossref bibliographic title search | **HIT** — correct DOI `10.1063/1.438337`, full record + 11 refs |
| Crossref JATS abstract | **HIT** — full abstract |
| OpenAlex `abstract_inverted_index` | **HIT** — independent confirmation, byte-identical |
| AIP article-PDF (`5268_1_online.pdf`) via curl | **403**, Cloudflare `cf-mitigated: challenge` |
| AIP abstract page via curl | **403** |
| AIP abstract page via WebFetch | **403** |
| AIP landing + PDF via `r.jina.ai` text proxy | **200 but useless** — returns the Cloudflare interstitial |
| Unpaywall | `is_oa: false`, `oa_locations: []` |
| OpenAlex OA locations | `best_oa_location: null`, `any_repository_has_fulltext: false` |
| Semantic Scholar `openAccessPdf` | `status: CLOSED` |
| CORE API | HTTP 500 (malformed server-side query on phrase search) |
| scholar.archive.org | 200, no work/PDF links |
| Wiley `pdfdirect` (Avery 2003 SW paper) | **403** Cloudflare |
| IOP PDF (1982 Sturmian many-centre) | **Radware bot-manager CAPTCHA** |
| Windsor ETD, aau.dk VBN, academia.edu, ResearchGate | **403 / SPA shell, no file** |
| **arXiv quant-ph/0403195** | **HIT** — full PDF, `pdftotext` clean; supplied the equations |
| **muroran-it repo (Koga 1985)** | **HIT on retry** — 406 with a rich `Accept` header, **200 with bare `-A Mozilla/5.0`**; PDF is a CJK-encoded scan, `pdftotext` yielded only 4 KB, unusable |
| Monkhorst CV (UF) | **HIT** — confirmed citation, author spelling, title variant; no companion long paper exists |
| Web search for the body | Returned a confident paraphrase citing "their equation (1.108)" — **impossible in a 2-page note; discarded as unreliable.** Not used anywhere above. |

*Reusable operational note:* the `-A Mozilla/5.0` + minimal-header form beat the
richer header set on the NII repo (406 → 200). Over-specifying `Accept` is
itself a failure mode.

---

## 6. What the corpus must now change

**Concrete edits. None applied — this scan was read-only.**

### E1 — `memory/avery_method_and_prior_art_gaps.md` point 3: RETRACT AS PHRASED (required)

Current: *"Conditioning is a blind spot in the whole Avery canon. ... The
linear-algebraic consequences ... are simply not studied."*

This is now false as a statement about the field. Replace with a two-part
claim, e.g.: *"Linear dependence was addressed head-on in the momentum-space
Sturmian literature — Novosadov, revived by Monkhorst & Jeziorski (JCP 71,
5268, 1979) — but as a STRUCTURAL claim of absence, not a quantitative one:
the re-posing never inverts the metric, so overcompleteness is inert. What
remains unstudied is the metric's SPECTRUM and its growth with basis size,
which only becomes a cost when the problem is converted to a standard
eigenproblem for block-encoding."*

### E2 — `memory/avery_method_and_prior_art_gaps.md` point 2: KEEP, NARROW (required)

"No prior art for secular-matrix norm growth" survives. Narrow it to the
1-norm / norm-growth axis explicitly, and attach the honest cap: the M–J leg is
abstract+reference-list+page-count, not a read of the body.

### E3 — `memory/avery_method_and_prior_art_gaps.md` §8: REPLACE the suspension (required)

Points 2 and 3 are no longer "suspended pending a read." Point 2 is reinstated
narrowed; point 3 is retracted and rewritten. Per CLAUDE.md §13.11 rule 9,
replace the text and move the superseded version to its history home; per the
`cited_by` rule, whatever entry carries point 3 must list Paper 60 as a
dependent if Paper 60's framing rests on it.

### E4 — Paper 60: ADD THE CITATION (required; PI call on wording)

Paper 60 currently cites **neither** M–J 1979 nor Novosadov anywhere (45
bibitems, zero hits on "Monkhorst"/"Jeziorski"/"Novosadov"). It should:

- cite M–J 1979 in `sec:molecular` as **prior art for the structural claim**
  that the momentum-space posing has no linear-dependence problem, and
- state the trade in one sentence: the conditioning multiplier is an artifact
  of demanding a *standard* eigenproblem for block-encoding; the classical
  literature avoids it by never inverting `S`, at the price of the outer
  nonlinear root search that `eq:secular` is designed to eliminate.

This **strengthens** Paper 60 rather than weakening it: it names the exact
boundary between what is inherited and what is the corpus's own, and it
converts an unexamined asymmetry into a stated trade. **Paper 60's own prose
is already careful** — the load-bearing over-claim lives in the memory file,
not in the paper. The `sec:quantum` novelty claim (the quantum algorithm) is
untouched by this read.

### E5 — `debug/lit_scan/sturmian_conditioning_prior_art_memo.md`: SUPERSEDE THE BOTTOM LINE (required)

Its "Bottom line up front" reads *"no genuine prior-art hit was found on either
Q1 or Q2 ... The two repo claims survive this scan intact."* That verdict was
reached without this paper in scope. Add a dated superseding header pointing
here: **one of the two claims does not survive.**

### E6 — Bibliographic hygiene (mechanical)

Record in the bib: DOI `10.1063/1.438337`; author **Jeziorski** (not the
Crossref/S2 "Jeriorski"); journal title form with "or ... many-center" (not the
CV's "and ... Multi-Center"). Flag `10.1063/1.438305` as a wrong DOI in case it
propagated from the task framing.

### E7 — NEW LEVER, queued not claimed (PI call)

The three-centre one-electron integrals are obtained in this lineage **without
three-centre integral evaluation**, via the completeness sum over intermediate
Fock-sphere states (products of two-centre SW integrals). This bears directly
on the corpus's three-centre wall (`memory/polyatomic_state_of_play.md`; §3
rows "Dropping the three-centre ERI block", "Closing the *reducible* half").

**Do not record this as a breach.** Two caveats, both load-bearing:
1. It is the **one-electron nuclear-attraction** three-centre integral, **not
   the two-electron three-centre ERI**, which is the corpus's actual wall.
2. The lineage paper itself reports that for `Nion >= 3` the three-centre terms
   "cannot be neglected and have to be accounted for on an equal footing" —
   i.e. the "no many-centre integral problems" advertisement **weakens at three
   or more centres**, which is independent corroboration that the corpus's
   three-centre wall is real.

The honest status is: a **named, published mechanism** for the one-electron
half of the three-centre problem that the corpus has not tested, with a stated
reason it may not transfer. That is a diagnostic probe, not a result.

---

## 7. One-line answers

1. **Claim:** a property of a *posing*, asserted in a 2-page note for the
   many-centre **one-electron** problem in a shared-scale Fock/Sturmian
   momentum-space basis, where eigenvalues come from "diagonalizations of
   simple overlap matrices". Not a theorem.
2. **Prior art:** for the **topic**, yes, decisively — retract "blind spot".
   For the **analysis** (condition numbers, spectra, norm growth), no evidence,
   on an abstract-plus-reference-list argument, not a read of the body.
3. **Resolution:** same pencil, same degeneracy; they never invert the metric,
   GeoVac must, because a quantum block-encoding wants a standard eigenproblem.
4. **Fock and Shibuya–Wulfman:** both cited directly, refs 3 and 5 of 11.
