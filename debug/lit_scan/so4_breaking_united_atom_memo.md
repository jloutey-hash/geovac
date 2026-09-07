# Lit scan — SO(4)-continuous-breaking / decay-length front / geometric-mean tail law

**Date:** 2026-09-06 · **Type:** Adversarial prior-art scan (no code/paper changes)
**Scope:** CLAUDE.md v5.10.5, Paper 58 §"The continuous side: the decompactification front"
(`papers/group2_quantum_chemistry/paper_58_abelian_residue.tex`, lines ~744-812),
CHANGELOG v5.10.2, `debug/sprint_decompactification_R_sweep_memo.md`.

**Method:** WebSearch (multiple rounds) + WebFetch (OSTI, Crossref API) for primary-source
confirmation. Every citation below was independently verified via at least one authoritative
source (Crossref API, OSTI, or a journal abstract page/PubMed record); none is asserted from
search-summary text alone without a matching structured record. No fabricated arXiv IDs/DOIs.

---

## Ranked findings

### Claim 1 — SO(4) continuous breaking with a measurable front

1. **T. P. Grozdanov and E. A. Solov'ev, "Separated- and united-atom limits for dynamical
   adiabatic states," Phys. Rev. A 44, 5605 (1991).** DOI 10.1103/PhysRevA.44.5605.
   Confirmed via OSTI record (biblio 5924982) and PRA abstract page. Abstract (paraphrased
   from OSTI): for one-electron collision systems, the separated-atom limit of the dynamical
   adiabatic states is obtained via the **O(4) symmetry of Coulomb bound states** (mapped
   onto hydrogen in crossed E/B fields); the united-atom limit predicts overlapping manifolds
   of potential-energy curves. **MOST DANGEROUS single paper for claim 1** — it is literally
   "O(4) symmetry, continuously connecting united-atom and separated-atom limits of the
   same two-center Coulomb problem." Difference from our claim: framed for time-dependent
   collision (dynamical adiabatic states along a classical trajectory), not the static
   bound-molecule electronic structure; I found no evidence in the abstract/available
   secondary material of a measurable front defined via a projector-commutator / principal-
   angle criterion, nor of the n/Z-vs-n²/Z distinction, nor of a geometric-mean law for
   unequal centers.

2. **E. A. Solov'ev, "The advanced adiabatic approach and inelastic transitions via hidden
   crossings," J. Phys. B: At. Mol. Opt. Phys. 38, R153–R194 (2005).** DOI
   10.1088/0953-4075/38/12/r01. Confirmed via Crossref API (`api.crossref.org` bibliographic
   search returned title/authors/venue/volume/page/DOI as an exact match). Topical review of
   the "hidden crossings" program Solov'ev started in 1981: the two-Coulomb-center problem's
   adiabatic terms, continued into the complex R-plane, with the prolate-spheroidal
   separation constant (which reduces to l only at R=0 — the same object our Paper 11 cites)
   tracked continuously as a function of R. This is the closest **qualitative** ancestor of
   "the united-atom SO(4)/separation-constant degeneracy breaks continuously with R" — a
   40+-year research program built on exactly that object. It does not (as far as secondary
   sources let me confirm) name it "SO(4) breaking," does not use a decay-length front, and
   does not carry a geometric-mean tail law.

3. **T. P. Grozdanov and E. A. Solov'ev, "Dynamical adiabatic theory of atomic collisions:
   The structure of hidden, avoided, and L₃ crossings," Phys. Rev. A 90, 032706 (2014).**
   DOI 10.1103/PhysRevA.90.032706. Confirmed via Crossref API. Same lineage, modern
   installment — establishes the program is still active and the underlying object (adiabatic
   terms of the two-center Coulomb problem, continuous in R) is standard machinery in atomic
   collision theory, not something GeoVac discovered independently as an *object*.

4. **von Neumann & Wigner (1929) + Herzberg (1950)** — see Claim 3 below; these are the
   textbook correlation-diagram/non-crossing-rule backbone that claim 1's "m survives, l
   dissolves" framing sits directly on top of, and Paper 58 already cites both correctly.

**No paper found** stating the specific quantitative content of claim 1: the front
R*(n) ≈ 2.14 n/Z (median, n ≥ 2), the explicit rejection of the mean-radius criterion
n²/Z with a measured 4× drift, or the fitted exponent R*(n) ∝ n^0.98 landing on the
decay-length (not mean-radius) scaling.

### Claim 2 — decay-length linear-in-R scaling + geometric-mean tail law

5. **Y. S. Kim and R. G. Gordon, "Atomic Distortion and the Combining Rule for Repulsive
   Potentials," Phys. Rev. A 5, 1708 (1972).** Confirmed via WebSearch abstract summary
   (PRA abstract page). Refines the older **geometric-mean combining rule** for the
   hardness/decay parameter of Born–Mayer-type short-range repulsion between *unlike*
   closed-shell species — repulsion energy from overlap of exponentially decaying electron
   clouds of two different species combines via (in the simplest, pre-Kim-Gordon case) a
   geometric mean of the two atoms' decay/hardness parameters, which Kim-Gordon's atomic-
   distortion model reduces to as a special case. **Most relevant single mechanism-level
   precedent for claim 2's geometric-mean law** — the general idea "two different exponential
   tails combine via a geometric mean of their decay lengths" is standard, decades-old content
   in the interatomic short-range-repulsion/combining-rule literature, just applied there to
   the repulsion *energy scale*, not to a projector-commutator mixing-*front distance*.

6. **Gor'kov & Pitaevskii (1963, cited via secondary sources) / Herring–Flicker exchange
   asymptotics** for H₂: 2J ≈ −1.641 R^(5/2) e^(−2R). Standard textbook result (see also
   Landau & Lifshitz, *Quantum Mechanics*, the exchange-interaction problem) that the
   long-range asymptotics of an overlap-driven quantity between two identical 1s orbitals is
   governed by the orbital's **exponential decay rate** (κ = 1 a.u. for H, i.e. decay length
   n/Z = 1), not its mean radius — this is the single-atom (equal-center) special case of the
   "decay length is what matters" principle claim 1/2 rely on, and it is standard in atomic
   collision theory (not GeoVac-specific). I could not locate a stated *geometric-mean*
   combining law for the *front position* (as opposed to the interaction magnitude) between
   two unequal decay lengths in this literature, nor a stated critical exponent-ratio
   threshold analogous to t_c = 2.664.

7. **Mulliken, Rieke, Orloff & Orloff (1949)** exact two-center 1s–1s overlap formulas
   (already cited in Paper 58 as `mulliken1949`) — the computational substrate every quantity
   above (and our own memo) is built on; not itself a claim about fronts or scaling laws.

**No paper found** stating R_rel ≈ 1.58 √(ℓ_A ℓ_B) for the *mixing front* between two
hydrogenic 1s tails, nor the absolute-45°-front existence condition t < t_c ≈ 2.664 from
S(0) = (2√t/(1+t))³ = 1/√2.

### Claim 3 — non-crossing rule / m as sorting label

8. **J. von Neumann and E. Wigner, "Über das Verhalten von Eigenwerten bei adiabatischen
   Prozessen," Phys. Z. 30, 467–470 (1929).** Confirmed via WebSearch (multiple independent
   hits agree on title/journal/pages/year; English translation in World Scientific's *20th
   Century Chemistry* Vol. 8). This is exactly the citation already used in Paper 58
   (`neumann_wigner1929`) — no correction needed.
9. **G. Herzberg, *Molecular Spectra and Molecular Structure I: Spectra of Diatomic
   Molecules*, 2nd ed. (Van Nostrand, 1950), Ch. VI** — the canonical correlation-diagram /
   non-crossing-rule textbook chapter. Already cited in Paper 58 (`herzberg1950`) with the
   correct chapter pointer — no correction needed.

---

## Verdicts

**(1) SO(4)-continuous-breaking-with-measurable-front: ADJACENT.** The *qualitative* claim
— that the united-atom degeneracy (the two-center Coulomb problem's hidden/dynamical
symmetry, whose separation constant reduces to l only at R=0) breaks continuously as a
function of R — is a real, 40+-year research program (Solov'ev 1981 onward; Grozdanov &
Solov'ev PRA 1991/2014; review J. Phys. B 2005), and Grozdanov–Solov'ev (1991) explicitly
invokes O(4) symmetry to connect the two limits. This is close enough that Paper 58 should
cite it, not just von Neumann–Wigner/Herzberg, when it makes the "continuous breaking"
claim. What appears NOT already published: a measurable front defined via the
projector-commutator/principal-angle (1/√2) criterion, or the explicit "m stays exact / l
alone is lost" packaging tied to a numeric front location.

**(2) decay-length linear-in-R scaling + geometric-mean tail law: OPEN (mechanism-adjacent).**
The specific numeric laws (R* ≈ 2.14 n/Z, exponent ≈0.98–1.02, R_rel ≈ 1.58√(ℓ_Aℓ_B),
t_c ≈ 2.664) were not found anywhere searched. But the *underlying principles* — that
decay length (not mean radius) sets long-range overlap/exchange scales (standard atomic
collision theory, e.g. the H₂ exchange asymptotics), and that two different exponential
decay/hardness parameters combine via a geometric mean (Born–Mayer combining-rule
literature, refined by Kim–Gordon 1972) — are each separately well established. Recommend
citing both families in Paper 58 to pre-empt a reviewer flagging "decay length, not radius"
and "geometric mean of two lengths" as textbook ideas, while keeping the specific measured
constants/threshold as the paper's own contribution.

**(3) non-crossing-rule / m-as-sorting-label: COLLISION (confirmed, correctly cited).**
Both canonical references (von Neumann & Wigner 1929; Herzberg 1950 Ch. VI) are already
in Paper 58's bibliography with correct bibliographic detail. No action needed beyond
possibly adding the Solov'ev hidden-crossings lineage as a "closer" adjacent citation for
the *continuous* half of the picture (claim 1), since Herzberg/von Neumann–Wigner alone
cover only the discrete/topological (non-crossing, label-sorting) half.
