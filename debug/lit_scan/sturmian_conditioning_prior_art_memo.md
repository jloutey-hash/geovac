# Prior-art scan: conditioning / linear dependence of multi-center Sturmian bases

**Date:** 2026-09-11
**Purpose:** Test two load-bearing claims currently asserted in repo memory
(`avery_method_and_prior_art_gaps.md`, informing Paper 60):

> "There is NO prior art for secular-matrix norm growth."
> "Conditioning is a blind spot in the whole Avery canon."

against the WIDER Sturmian / hyperspherical-harmonic quantum-chemistry
literature (Aquilanti-Cavalli-Coletti-Calderini lineage), the Shibuya-Wulfman
original papers, the principal-angles lineage (Amos-Hall, King et al.,
West-Ruedenberg), and standard remedies for basis-set linear dependence.

**Method:** web search + direct fetch of arXiv PDFs where available (full text
searched for "condition", "conditioning", "ill-conditioned", "linear
dependence", "singular value", "overlap matrix eigenvalue"). Where only an
abstract/secondary description was reachable, this is flagged explicitly.
No paywalled book text (the two Avery books) was reachable in full; only
search-indexed snippets and reviews of them were available.

**Bottom line up front:** no genuine prior-art hit was found on either Q1 or
Q2. The decision-gate condition ("a single genuine prior-art hit... would
force Paper 60 to cite rather than claim") is NOT triggered. The two repo
claims survive this scan intact, with the caveats noted below (untested
book text; a plausible nuclear-physics-Sturmian claim of the opposite
sign that needs a closer look if it becomes load-bearing elsewhere).

---

## Q1. Aquilanti / Cavalli / Coletti / Calderini lineage — VERDICT: ABSENT

Searched extensively: "Sturmian orbitals and molecular structure" (Avery,
Ostrovsky, Coletti, Aquilanti — *J. Mol. Struct. THEOCHEM*, 2004, credits
Coletti + Avery + Tarantelli); "Hyperspherical harmonics as Sturmian orbitals
in momentum space: a systematic approach to the few-body Coulomb problem"
(Aquilanti, Cavalli, Coletti, Di Domenico, Grossi — *Int. Rev. Phys. Chem.*);
"Hydrogenic orbitals in momentum space and hyperspherical harmonics: Elliptic
Sturmian basis sets" (Aquilanti et al., *Int. J. Quantum Chem.* 93, 2003);
"d-Dimensional Kepler-Coulomb Sturmians and Hyperspherical Harmonics as
Complete Orthonormal Atomic and Molecular Orbitals" (Aquilanti group, book
chapter); "Molecular integrals for Slater type orbitals using Coulomb
Sturmians" (2013); the review "How specific exponential type orbitals
recently became a viable basis set choice" (arXiv:1010.5425 — full PDF
fetched and searched directly).

Every one of these sources is about: (a) classification/completeness of
Kepler-Coulomb Sturmian sets, (b) closed-form/algebraic evaluation of
two-center overlap and multicenter integrals via the momentum-space (Fock
projection) route, (c) convergence of the *energy* with the number of basis
functions or angular-momentum cutoff. None of the reachable text or
abstracts contains any treatment of: the overlap (or Shibuya-Wulfman) matrix's
condition number, its eigenvalue spectrum, near-linear-dependence as the
number of centers or basis functions grows, or any explicit large-N /
large-basis scaling law for the metric itself. The one recurring theme that
could be *mistaken* for a conditioning result is repeated claims that
Sturmian bases have "no linear dependency problem" compared to Gaussians
(see Q1a below) — this is a qualitative, single-sentence contrast claim, not
an analysis of how conditioning scales, and it is not attributed to the
Aquilanti group specifically in the sources that state it.

**Q1a — a claim worth flagging even though it doesn't trigger the gate.**
Independent searches on Coulomb-Sturmian bases in *nuclear* few-body physics
(not the Aquilanti chemistry lineage, but the same underlying orbital family)
turned up secondary-source language of the form: "a major caveat in the
Gaussian basis set is its linear dependence, growing with the number of
basis functions... [but] Sturmian functions are virtually free of this
problem... no linear dependency problem occurred during calculations with
Sturmian basis sets." This is a search-engine paraphrase, not a verified
quotation — I could not get readable full text out of the two PDF candidates
(arXiv:1208.4156, "Coulomb-Sturmian basis for the nuclear many-body problem",
returned as corrupted/binary to the fetch tool) to confirm wording or see
whether any quantitative scaling accompanies it. Flag: **UNCERTAIN**, and
notably this is a *single-center* nuclear basis (one shared origin, no
multi-center SW-matrix structure), so even if confirmed it would not
contradict the GeoVac claim, which is specifically about the *multi-center*
metric's conditioning. Recommend not citing this without reading the actual
paper text; it's listed here so a future pass doesn't waste time
re-discovering it.

**Conclusion: ABSENT.** No conditioning/linear-dependence/spectrum analysis
of the metric found anywhere in the Aquilanti-Cavalli-Coletti-Calderini
corpus reachable by this scan.

---

## Q2. The Avery canon — VERDICT: ABSENT (directly verified for the one
full-text-reachable source; consistent absence across everything else)

**Herbst, Avery, Dreuw, "Quantum chemistry with Coulomb Sturmians:
Construction and convergence of Coulomb Sturmian basis sets at Hartree-Fock
level," Phys. Rev. A 99, 012512 (2019); arXiv:1811.05777.** This is the one
Q2 source with full-text PDF access (fetched directly, both v1 abstract page
and the PDF itself, and re-fetched a second time with an explicit
find-the-phrase instruction). Direct result: **"No such mentions were
found... discusses energy convergence rates as the basis size increases, but
does not address the numerical conditioning of the basis set or the
behavior of overlap matrix eigenvalues... no references to linear dependence
problems, singular value decomposition, or ill-conditioning."** The paper's
own framing of "convergence" is exclusively energy-vs-basis-size (its title
literally promises "construction and convergence," and delivers only the
former in the numerical-stability sense). This is the paper most likely to
contain a conditioning result if one existed in the Avery lineage (it is the
most recent, most numerically-oriented Avery-coauthored paper on Coulomb
Sturmians), and it does not.

**Books — "Hyperspherical Harmonics and Generalized Sturmians" (Avery,
Kluwer/Springer, Progress in Theoretical Chemistry and Physics vol. 4) and
"Generalized Sturmians and Atomic Spectra" (J. Avery & J. Avery, World
Scientific, 2006).** Full text not reachable (no open PDF, Google
Books/Amazon/Goodreads listing pages only). Targeted searches for "condition
number" + book title, and "Avery" + "Journal of Mathematical Chemistry" +
"condition number"/"conditioning"/"singular," returned **zero hits** in
either case — not "hits that don't address it," but no matching indexed
content at all, which is itself informative (Google-indexed snippets would
likely surface the phrase if it appeared in the book's search-inside index,
as happens for other books). This is the weakest link in the chain: **I did
not read the book text**, so this is UNCERTAIN-leaning-ABSENT rather than a
verified negative. Flag clearly for the user: if Paper 60's claim needs to
survive a hostile citation-reviewer pass, the two books remain the one
unverified gap, since they are exactly the kind of source (`the papers papers
may not be online`, per the task's own framing) that a web scan cannot
fully rule out.

**Generalized Sturmian Method papers** (e.g., "The Generalized Sturmian
Method for Calculating Spectra of Atoms and Ions," *J. Math. Chem.* 2003;
"Many-particle Sturmians," *J. Math. Chem.*; "Many-center Coulomb Sturmians
and Shibuya-Wulfman integrals," Avery, *Int. J. Quantum Chem.* 98, 2004,
DOI 10.1002/qua.10820; "Generalized Sturmian Functions in prolate spheroidal
coordinates," arXiv:2006.06616 — full PDF fetched and searched directly).
All describe: isoenergetic configuration construction, the vanishing of the
kinetic-energy term when V₀ = nuclear attraction, completeness/convergence
of the configuration expansion, and (in the prolate-spheroidal paper)
explicit "energy convergence" discussion with **no** treatment of secular-
matrix conditioning, ill-conditioning, or linear dependence — confirmed by
direct fetch: "does not contain substantive discussion of condition numbers,
ill-conditioning, or linear dependence issues as the basis size increases...
emphasizes 'energy convergence' characteristics but does not address the
numerical conditioning challenges."

**Conclusion: ABSENT**, with the one honest caveat that the two books were
not read in full (search-index absence is suggestive, not dispositive). No
sentence anywhere in the Avery canon reachable by this scan states or
implies anything about how the secular matrix's norm or condition number
scales with basis size K. **This directly confirms the repo's existing
claim rather than refuting it** — the decision gate is not triggered.

---

## Q3. Shibuya & Wulfman (1965) and successors — VERDICT: ABSENT

**Original paper:** Shibuya, T. and Wulfman, C. E., "Molecular orbitals in
momentum space," *Proc. R. Soc. Lond. A* **286**, 376-389 (1965), DOI
10.1098/rspa.1965.0151 (20 July 1965). Confirmed citation via ADS/
Semantic Scholar cross-check. Content (from abstract + secondary
descriptions, full text not reachable — this is a 1965 Royal Society paper,
almost certainly not open-access): addresses the cusp problem in molecular
wavefunctions by moving to momentum space, uses R⁴ hyperspherical harmonics
as basis, derives algebraic integral formulas for many-center Coulomb
potentials. **No indication anywhere (title, abstract, or the many later
papers that cite and summarize it) that the paper analyzes the resulting
matrix's spectrum or conditioning** — its entire contribution is the
elimination of continuum-wavefunction/cusp difficulties via the Fock
projection, an *integral-evaluation* result, not a *linear-algebra-of-the-
metric* result. Flagged: reached only via abstract/secondary description,
not the primary text.

**Successors:** "Derivation of a general formula for the Shibuya-Wulfman
matrix" (~2004; author not fully resolved by this scan, described in
secondary sources as following Weniger's Fourier-transform methods) —
description available says it "presents an explicit and general formula for
the Shibuya-Wulfman matrix, which may easily be programmed." This is a
closed-form-expression result (algebraic evaluation), not a spectral/
conditioning analysis. "Many-center Coulomb Sturmians and Shibuya-Wulfman
integrals" (Avery, 2004) — same pattern, extends the integral evaluation to
general many-center potentials and to crystal/band-structure calculations;
no conditioning content found. No paper anywhere in this lineage was found
to analyze the SW matrix's eigenvalue spectrum, condition number, or
large-basis / large-N-center scaling.

**Conclusion: ABSENT.** The entire post-1965 SW-matrix literature reachable
by this scan is about evaluating the matrix elements efficiently and in
closed form (the "how do you compute it" question), never about what the
resulting matrix looks like spectrally as the basis grows (the "how well-
conditioned is it" question) — precisely the distinction the task asked me
to hold apart from convergence-rate results, and precisely GeoVac's claimed
gap.

---

## Q4. Principal angles / canonical correlations between center subspaces —
VERDICT: FOUND (the mathematical formalism), but ABSENT for conditioning/
R→∞ use

**Amos, A. T. and Hall, G. G., "Single determinant wave functions,"
*Proc. R. Soc. Lond. A* **263**, 483-493 (1961), DOI 10.1098/rspa.1961.0175.**
Confirmed citation. This is the origin of the "corresponding orbitals"
transformation: given two sets of (possibly non-orthogonal, different-
center) spin-orbitals, there is a unitary transform of each set making the
cross-overlap matrix diagonal — i.e., exactly a principal-angle
decomposition (the diagonal entries are the cosines of the principal
angles), used to write a tractable Slater-Condon-type rule for matrix
elements between two non-orthogonal determinants. Full text not reached
(1961 Royal Society paper); description from secondary sources only.

**King, H. F., Stanton, R. E., Kim, H., Wyatt, R. E., and Parr, R. G.,
"Corresponding orbitals and the nonorthogonality problem in molecular
quantum mechanics," *J. Chem. Phys.* **47**, 1936-1941 (1967).** Confirmed
citation via ADS. Generalizes Amos-Hall to arbitrary spin-orbital sets:
"given any two sets of spin orbitals, there exist equivalent sets such that
their overlap matrix is diagonal" — again the principal-angle construction,
again deployed purely to generalize the Slater-Condon rules for computing
⟨Ψ|Ô|Ψ'⟩ between non-orthogonal determinants. **Neither Amos-Hall nor King
et al. addresses conditioning, the large-basis limit, or R→∞ behavior** — the
principal angles are a computational device for one matrix-element
evaluation, evaluated once per pair of determinants, not tracked as a basis
grows or as centers separate.

**West, A. C., Schmidt, M. W., Gordon, M. S., and Ruedenberg, K., "A
comprehensive analysis of molecule-intrinsic quasi-atomic, bonding, and
correlating orbitals. I. Hartree-Fock wave functions," *J. Chem. Phys.*
**139**, 234107 (2013).** Full text blocked (HTTP 403 from the publisher);
only abstract-level and secondary descriptions reachable. What is confirmed:
the quasi-atomic orbital (QUAO) method **does** use an SVD between the
molecule's valence MO subspace and a free-atom minimal-basis subspace, and
the singular values of that intersection matrix **are** the cosines of
principal angles between the two subspaces — so the mathematical machinery
of Q4 is present and in active current use in mainstream quantum chemistry.
But its use here is orbital *localization/bonding analysis within one
converged molecular wavefunction* (finding the atom-like orbitals a
Hartree-Fock MO space contains), not conditioning analysis of a multi-
center Sturmian basis, and not an R→∞ or basis-size-scaling study. **Flag:
I could not reach the paper's full text**, so I cannot rule out a
conditioning remark appearing somewhere in its ~15 pages; the abstract and
all secondary sources found describe only the bonding-analysis application.
This is the single most important UNCERTAIN item in this scan given the
task's explicit naming of this reference.

**Conclusion for Q4: FOUND** (the principal-angle/canonical-correlation
formalism itself is well established in the quantum-chemistry lineage named
in the task — Amos-Hall → King et al. → West-Ruedenberg/QUAO), **but ABSENT**
for the specific question asked (conditioning, large-basis limit, R→∞) in
every source whose content could be verified. The one item flagged
UNCERTAIN (West-Ruedenberg full text) should be read directly before this
verdict is treated as fully closed.

---

## Q5. Standard remedies for near-linear-dependence — current best practice,
and does any of it preserve angular-momentum block structure?

**Canonical orthogonalization with an eigenvalue threshold (Löwdin).**
Confirmed as current best practice in mainstream (Gaussian-basis) quantum
chemistry: diagonalize the overlap matrix S, discard eigenvectors with
eigenvalue below a threshold (typical figures found: "true" linear
dependence below ~1e-7; numerical instabilities generally avoided once all
eigenvalues exceed ~1e-6 to 1e-7; more stringent 1e-8 to 1e-10 thresholds
used in double precision for near-linear-dependence removal), then rescale
the surviving eigenvectors by λ^(-1/2). Also found: "Application of Löwdin's
canonical orthogonalization method to the Slater-type orbital
configuration-interaction basis set" (Jiao, *Int. J. Quantum Chem.*, 2015)
and "Löwdin's Canonical Orthogonalization: Getting Round the Restriction of
Linear Independence" — both confirm this is the standard tool specifically
for STO/ETO-type (Sturmian-adjacent) bases too, not just Gaussians.

**Does canonical orthogonalization preserve angular-momentum (l-selection)
block structure? NO, not for a genuinely multi-center basis.** This is
confirmed both by direct search results and by the repo's own prior finding
(Track DF record, cited in CLAUDE.md §3.5: "Löwdin orthogonalization of
mixed-exponent bases destroys Gaunt sparsity (1711 vs 120 Pauli, 14×
inflation)"). The mechanism found in the literature: canonical
orthogonalization diagonalizes the *full* overlap matrix; this only
respects an l-block structure if the overlap matrix is *already*
block-diagonal in l before orthogonalizing. For a single center, Sturmians
of different l are exactly orthogonal (the angular parts are orthogonal
spherical harmonics), so the overlap matrix is block-diagonal in l and
canonical orthogonalization acts within each l-block harmlessly. For a
multi-center basis, cross-center overlaps mix different l's (there is no
shared angular frame), so the full S is NOT block-diagonal in l, and
diagonalizing it produces eigenvectors that are generic linear combinations
across l — i.e., canonical orthogonalization is confirmed (by the general
symmetry-breaking literature on Löwdin schemes, e.g. "A Note on the
Symmetry Properties of Löwdin's Orthogonalization Schemes," and directly by
the repo's own measurement) to be an l-block-destroying operation in the
multi-center case, not an l-preserving one.

**Dual/biorthogonal bases (Artacho & del Bosch, *Phys. Rev. A* 43, 5770,
1991, "Nonorthogonal basis sets in quantum mechanics: representations and
second quantization").** Confirmed citation and content: this is precisely
the alternative that *avoids* orthogonalizing at all — it builds a tensorial
(covariant/contravariant, i.e. biorthogonal) formalism for doing quantum
mechanics and second quantization directly in a non-orthogonal basis, using
the metric (overlap matrix) and its inverse rather than diagonalizing it
away. **This is the one Q5 candidate that structurally preserves whatever
block structure (including l-selection) the original non-orthogonal basis
had**, precisely because it never mixes basis functions to force
orthogonality — it keeps the original (e.g. per-center, per-l) basis
functions as the covariant frame and works with a dual (contravariant) frame
built from S⁻¹, which inherits the same index/block structure S has. This
is the practice-relevant answer to the "does any of it preserve l-selection"
half of Q5: canonical orthogonalization (the dominant mainstream practice)
does NOT; the biorthogonal/dual-basis route (much less commonly adopted in
production codes, but exactly the mathematically appropriate tool for a
symmetry-structured non-orthogonal basis) DOES, by construction, because it
sidesteps the diagonalization step entirely.

**Even-tempered/Gaussian coverage-vs-linear-dependence practice.**
Confirmed as current mainstream practice for diffuse-function-heavy
(Gaussian, not Sturmian) basis sets: rule-of-thumb removal of exponents
below ~0.1, even-tempered exponent sequences chosen for transferability,
and (more modern) pivoted Cholesky decomposition of the overlap matrix to
prune over-complete auxiliary basis sets (arXiv:1911.10372, "Communication:
Curing basis set overcompleteness with pivoted Cholesky decompositions") —
this is a genuine, actively-used current remedy, but it is a Gaussian-basis
technique, not part of the Sturmian/hyperspherical-harmonic lineage, and (per
its own framing) it operates on the raw overlap matrix without any special
provision for preserving angular-momentum block structure either.

**Q5 summary:** current best practice for multi-center non-orthogonal
molecular bases is canonical (Löwdin) orthogonalization with an eigenvalue
cutoff — well documented, standard, and confirmed NOT to preserve l-block
structure once centers differ. The dual/biorthogonal route (Artacho & del
Bosch) is the one approach found that would preserve l-selection, precisely
because it declines to orthogonalize. Even-tempered/Cholesky-pruning
practice is Gaussian-specific and, like canonical orthogonalization, makes
no attempt to protect angular structure.

---

## Full citation list

1. V. Aquilanti, S. Cavalli, C. Coletti, D. Di Domenico, G. Grossi,
   "Hyperspherical harmonics as Sturmian orbitals in momentum space: a
   systematic approach to the few-body Coulomb problem," *Int. Rev. Phys.
   Chem.* (secondary description only; full text not reached).
2. J. Avery, D. Z. Ostrovsky, C. Coletti, V. Aquilanti, "Sturmian orbitals
   and molecular structure," *J. Mol. Struct.: THEOCHEM* (2004) (abstract
   only).
3. V. Aquilanti et al., "Hydrogenic orbitals in momentum space and
   hyperspherical harmonics: Elliptic Sturmian basis sets," *Int. J.
   Quantum Chem.* **93** (2003), Wiley DOI 10.1002/qua.10508 (abstract
   only).
4. "d-Dimensional Kepler-Coulomb Sturmians and Hyperspherical Harmonics as
   Complete Orthonormal Atomic and Molecular Orbitals," book chapter,
   Aquilanti group (secondary description only).
5. "How specific exponential type orbitals recently became a viable basis
   set choice," arXiv:1010.5425 — **full PDF fetched**; no conditioning
   content found (only convergence-with-basis-size discussion).
6. F. E. Herbst, J. Avery, A. Dreuw (order as indexed; author list per
   arXiv/PRA record), "Quantum chemistry with Coulomb Sturmians:
   Construction and convergence of Coulomb Sturmian basis sets at
   Hartree-Fock level," *Phys. Rev. A* **99**, 012512 (2019);
   arXiv:1811.05777 — **full PDF fetched and searched directly**; result:
   no mention of condition number/conditioning/ill-conditioning/linear
   dependence/singular values/overlap-matrix eigenvalues; only
   energy-convergence-rate content.
7. J. Avery, "Many-center Coulomb Sturmians and Shibuya-Wulfman integrals,"
   *Int. J. Quantum Chem.* **98** (2004), DOI 10.1002/qua.10820 (abstract
   only).
8. "Generalized Sturmian Functions in prolate spheroidal coordinates,"
   arXiv:2006.06616 — **full PDF fetched**; energy-convergence content
   only, no conditioning.
9. J. S. Avery, *Hyperspherical Harmonics and Generalized Sturmians*,
   Progress in Theoretical Chemistry and Physics vol. 4 (Kluwer/Springer) —
   **book text not reached**; zero indexed hits for "condition number" in
   combination with the title.
10. J. Avery & J. Avery, *Generalized Sturmians and Atomic Spectra* (World
    Scientific, 2006) — **book text not reached**; same null-search result.
11. "The Generalized Sturmian Method for Calculating Spectra of Atoms and
    Ions," *J. Math. Chem.* (2003), DOI-indexed as 10.1023/A:1023204016217
    (abstract only).
12. T. Shibuya, C. E. Wulfman, "Molecular orbitals in momentum space,"
    *Proc. R. Soc. Lond. A* **286**, 376-389 (1965), DOI
    10.1098/rspa.1965.0151 (abstract/secondary description only; primary
    1965 text not reached).
13. "Derivation of a general formula for the Shibuya-Wulfman matrix" (~2004;
    author/venue not fully resolved by this scan; described as following
    Weniger's Fourier-transform methods) (secondary description only).
14. A. T. Amos, G. G. Hall, "Single determinant wave functions," *Proc. R.
    Soc. Lond. A* **263**, 483-493 (1961), DOI 10.1098/rspa.1961.0175
    (secondary description only).
15. H. F. King, R. E. Stanton, H. Kim, R. E. Wyatt, R. G. Parr,
    "Corresponding orbitals and the nonorthogonality problem in molecular
    quantum mechanics," *J. Chem. Phys.* **47**, 1936-1941 (1967) (abstract
    only).
16. A. C. West, M. W. Schmidt, M. S. Gordon, K. Ruedenberg, "A comprehensive
    analysis of molecule-intrinsic quasi-atomic, bonding, and correlating
    orbitals. I. Hartree-Fock wave functions," *J. Chem. Phys.* **139**,
    234107 (2013) — **full text blocked (HTTP 403)**; abstract + secondary
    description only. **Flagged as the least-verified source in this scan.**
17. E. Artacho, L. Milans del Bosch, "Nonorthogonal basis sets in quantum
    mechanics: Representations and second quantization," *Phys. Rev. A*
    **43**, 5770 (1991) (abstract-level content confirmed via search
    description).
18. "Communication: Curing basis set overcompleteness with pivoted Cholesky
    decompositions," arXiv:1911.10372 (abstract only; Gaussian-basis
    context, not Sturmian).
19. "Application of Löwdin's canonical orthogonalization method to the
    Slater-type orbital configuration-interaction basis set," Jiao, *Int.
    J. Quantum Chem.* (2015) (abstract only).
20. (For the Q1a nuclear-Sturmian flag, unconfirmed) "Coulomb-Sturmian basis
    for the nuclear many-body problem," arXiv:1208.4156 — fetch returned
    corrupted/binary content; not read.

---

## What would close the remaining uncertainty (for a future pass)

1. Read West, Schmidt, Gordon, Ruedenberg (2013) full text directly (via
   institutional access or a non-blocked mirror) — the single highest-value
   unresolved item, since the task named this reference specifically.
2. Get physical/library access to the two Avery books' index/bibliography
   sections (not just search-engine indexing) — the only way to make the
   Q2 ABSENT verdict fully dispositive rather than search-index-suggestive.
3. Read arXiv:1208.4156 properly (the fetch tool returned binary garbage) to
   confirm or refute the "no linear dependency problem" nuclear-Sturmian
   claim and see whether it comes with any quantitative backing — currently
   UNCERTAIN and flagged as not contradicting GeoVac's multi-center claim
   even if confirmed (different basis topology: single center, not SW-matrix
   multi-center).
