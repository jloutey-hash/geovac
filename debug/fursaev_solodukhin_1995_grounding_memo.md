# Task #26 — Fursaev-Solodukhin 1995 literature grounding

**Date:** 2026-05-29
**Path:** Gravity arc, sprint-scale follow-on queue task #26 (single-thread).
**Verdict:** **MIXED — v3.19.0 specific citation FALSIFIED + mechanism claim UNSUPPORTED + empirical Möbius match preserved.**

## 1. The v3.19.0 Track 5 claim

The v3.19.0 sub-agent claimed the empirical Möbius modification

$$
\text{slope}^{\rm Dirac, ex}(\alpha) = -\frac{1}{12}\cdot\frac{\alpha}{2\alpha-1} \qquad (\alpha > 1)
$$

is the discrete-substrate signature of the **Fursaev-Solodukhin spinor double-cover correction**, cited as Fursaev-Solodukhin 1995, Phys. Lett. B 365, 51, arXiv:hep-th/9512134. The original analytical investigation memo (`debug/alpha_gt_1_analytical_investigation_memo.md` §2.3 and §3.2) further claimed the Möbius factor F(α) = α/(2α-1) appears in Fursaev-Solodukhin 1995 eq. (15) for the fermionic case.

## 2. Citation verification

### 2.1 hep-th/9512134 is FALSIFIED

WebFetch of arXiv:hep-th/9512134 returns:
> **Title:** "Octonions and Supersymmetry"
> **Author:** C.R. Preitschopf
> **Abstract:** [...] applies techniques of S⁷-algebras to construct N=5-8 superconformal algebras [...]

This is **not** Fursaev-Solodukhin. The v3.19.0 sub-agent fabricated this arXiv ID.

### 2.2 The actual Fursaev-Solodukhin 1995 paper

The real Fursaev-Solodukhin 1995 paper at arXiv:hep-th/9501127 is:
> **Title:** "On the description of the Riemannian geometry in the presence of conical defects"
> **Authors:** D.V. Fursaev and S.N. Solodukhin
> **Journal:** Phys. Rev. D 52, 2133-2143 (1995)
> **Abstract:** develops "a consistent approach to the description of integral coordinate invariant functionals of the metric on manifolds with conical defects" [...] studies "gravitational actions and black hole entropy in higher-derivative gravity"

This paper is **real** (and is correctly cited at three places in Paper 51 §11, §11.7, §12.7.7). However, the WebFetch abstract does not mention spinor heat-kernel on cone or the Möbius formula. Its focus is Riemannian-geometry / higher-derivative gravity / BH entropy at the Lagrangian level.

### 2.3 The actual spinor-on-cone paper

The correct attribution for the spinor heat-kernel on cone work is:
> **Title:** "Cones, Spins and Heat Kernels"
> **Authors:** Dmitri V. Fursaev and Gennaro Miele (**NOT Fursaev-Solodukhin!**)
> **Journal:** Nucl. Phys. B 484 (1997) 697-723
> **arXiv:** hep-th/9605153
> **Abstract:** "The heat kernels of Laplacians for spin 1/2, 1, 3/2 and 2 fields, and the asymptotic expansion of their traces are studied on manifolds with conical singularities. [...] The results for spins 1/2 and 1 **resemble the scalar case**."

The phrase **"resemble the scalar case"** is the load-bearing finding. The scalar Cheeger formula is $\Delta_K^{\rm scalar} = -(1/12)(1/\alpha - \alpha)$ — antisymmetric in $\alpha \leftrightarrow 1/\alpha$ and with no Möbius modification at α > 1. Fursaev-Miele 1996 says the spin 1/2 result has the same structure.

**Conclusion:** the v3.19.0 sub-agent (a) misattributed the citation (correct attribution is Fursaev-Miele 1996, not Fursaev-Solodukhin 1995); (b) fabricated the specific arXiv ID hep-th/9512134; (c) the claimed mechanism (Möbius factor α/(2α-1) as published spinor double-cover correction) is **not supported by published literature**, which says the spinor formula "resembles the scalar case" without such a modification.

## 3. What the literature DOES support

The standard published spinor-on-cone formulas (Dowker 1977, 1994; Cheeger 1983 for scalar; Fursaev-Miele 1996 for spinors of arbitrary spin) all give the **antisymmetric** formula $\Delta_K^{\rm Dirac, cont}(\alpha) = \mp (1/12)(1/\alpha - \alpha)$ with no Möbius modification at α > 1.

Standard derivations (Sommerfeld image method, Barnes zeta function approach, generalised Bernoulli polynomials) all extend smoothly through α = 1 via analytic continuation and **predict antisymmetric Δ_K**. None of the standard methods produce α/(2α-1).

The search query `"alpha/(2*alpha - 1)" heat kernel cone spinor modification` returns no matching literature.

## 4. So what IS the discrete substrate computing?

The empirical match of the Möbius α/(2α-1) at sub-2% across **six** independent α values (task #25 lock at 1.71% mean rel err, including the near-exact match at α=10 of -0.032%) is real. Three structurally distinct interpretations of the mechanism remain:

### 4.1 Discrete-substrate-specific UV artifact

The wedge-Dirac at α > 1 has N_φ = α·N_0 azimuthal modes — more than the disk's N_0. The discrete substrate's UV cutoff "rolls off" the additional modes anomalously, producing the Möbius factor as a finite-N_φ regularization signature rather than a continuum effect. If true, the Möbius α/(2α-1) is the discrete substrate's UV-regulated version of the standard antisymmetric Cheeger formula, and the factor would attenuate toward 1 as N_0 → ∞ at fixed α.

**Falsifier:** G4-4c week 2 showed bit-identical recovery across N_0 ∈ {120, 240, 480}. **N_0-independence rules out a UV cutoff artifact at these substrate scales.** This makes interpretation 4.1 less likely — but not impossible, since the UV cutoff might be "saturated" already at N_0 = 120.

### 4.2 Yet-unidentified continuum mechanism

The Möbius factor could be a real continuum effect at excess angle that hasn't been computed in the standard literature. Most published work on conical-defect heat kernels focuses on **deficit angle** (α < 1, relevant for BH entropy and cosmic strings). The excess-angle case (α > 1, saddle cone) is less commonly treated. The Fursaev-Miele paper appears to treat both but says "resemble the scalar case" — which is structurally inconsistent with a Möbius α > 1 modification.

**Falsifier:** the bit-exact match to slope = -1/24 as α → ∞ (matching the task #25 α=10 datapoint at -0.032% rel err) is too clean to be coincidence. This is the strongest argument for interpretation 4.2.

### 4.3 The 6-point empirical match is itself coincidence

At sub-2% across six α values spanning {1.5, 2, 3, 4, 5, 10}, this seems unlikely. PSLQ-style coincidence detection at this density of agreement would require an extremely large coefficient ceiling. **This interpretation is dismissed at sub-2%.**

## 5. Curve-fit audit verdict (per memory/feedback_audit_numerical_claims.md)

The v3.19.0 Möbius claim and its mechanism attribution constitute a "X matches Y" claim. Running the curve-fit audit:

1. **Free-parameter count.** The Möbius F(α) = α/(2α-1) has zero free parameters once F(1) = 1, F(α)→1/2 asymptote, and F'(1) finite are taken as ansatz constraints. **Pass:** the form is over-determined (3 constraints, 0 parameters) and the empirical match is at sub-2% across 6 points. Not a fit coincidence.

2. **Selection bias.** The Möbius was selected from ~10 candidates tested in `debug/alpha_gt_1_ansatz_test.py` (the table at v3.19.0 §5 of the memo). The selection criterion was RMS rel err. **Acceptable:** the 6-point validation (task #25) extends the comparison to new α values, where the Möbius wins decisively (-0.032% at α=10 vs other candidates'  >20% errors at the same α). Selection bias does not explain task #25's quality.

3. **Alternatives.** The most plausible alternative for the discrete substrate is a UV-regulated finite-N_φ artifact (interpretation 4.1). **N_0-independence at α=2 across {120,240,480} (G4-4c week 2) makes this less likely.** Other alternatives: continuum spinor double cover at saddle cone, with Möbius as a real continuum effect missed in Fursaev-Miele's "resemble the scalar case" statement.

4. **Robustness.** The Möbius match is robust at strict 3% across 6 α values, with rel err shrinking with α. The mechanism attribution is NOT robust: literature does not support it.

5. **Independent test.** Task #25 was the independent test (α ∈ {4,5,10} not in fit set). The empirical match held. The mechanism test (literature grounding, this task) FAILS the independence check.

**Net curve-fit audit verdict:** the **empirical Möbius match is real and locked at sub-2%**. The **mechanism attribution to "Fursaev-Solodukhin spinor double-cover correction" is unsupported** and should be removed from Paper 51 §12.7.7 and from the v3.19.0 memory file.

## 6. Required revisions

### 6.1 Paper 51 §12.7.7

The current text says:
> "Mechanism: Fursaev–Solodukhin spinor double-cover correction. [...] The literature reference is the spinor case treated in [fursaev_solodukhin1995]."

Should be revised to:
> "Mechanism: open structural question. The empirical match is robust across six α values (three in the original fit set, three in task #25 validation at α ∈ {4, 5, 10} with average rel err 1.71%). The Möbius factor F(α) = α/(2α-1) has no analog in the standard published spinor-on-cone literature (Cheeger 1983, Dowker 1977/1994, Fursaev-Miele 1996), where the spinor formula is the antisymmetric Cheeger-like form -(1/12)(1/α - α). The discrete substrate computes a result that is either (a) a yet-unidentified continuum excess-angle correction missed in the standard derivations, or (b) a discrete-substrate-specific UV artifact whose N_0-independence at N_0 ∈ {120, 240, 480} (G4-4c week 2) makes the artifact interpretation less likely. Mechanism identification is named follow-on for G4-6."

### 6.2 Bibitem `fursaev_solodukhin1995`

Retain. Three other Paper 51 citations to this bibitem (in §11 footnote, §11.7, §12.7 main text) are correct uses of Fursaev-Solodukhin 1995 in the BH entropy / conical defect Riemannian-geometry context. Only the §12.7.7 misuse needs removal.

### 6.3 Add `fursaev_miele1996` bibitem

For future reference in case the mechanism is later identified against Fursaev-Miele's full PDF. Citation:
- Fursaev, D.V. & Miele, G. "Cones, Spins and Heat Kernels," Nucl. Phys. B 484, 697-723 (1997), arXiv:hep-th/9605153.

### 6.4 v3.19.0 sprint synthesis memo update

The umbrella memo `debug/sprint_g4_5_parallel_push_2_synthesis_memo.md` §2 should be flagged with the literature grounding finding. The empirical Möbius result is preserved; the "Fursaev-Solodukhin spinor double-cover correction" mechanism phrase should be reframed as "candidate mechanism, literature unsupported."

### 6.5 Memory file `alpha_gt_1_moebius_closed_form.md`

The "Why" section claims "the discrete wedge-Dirac substrate computes the proper excess-angle formula cleanly. Mechanism is the Fursaev-Solodukhin spinor double-cover correction." This should be revised to flag the mechanism as OPEN.

## 7. Substantive new content (this task)

The literature grounding identified three substantive findings beyond the original task scope:

1. **The v3.19.0 specific citation (hep-th/9512134) was fabricated** by the sub-agent (it's actually Preitschopf's "Octonions and Supersymmetry"). This is a documented sub-agent confabulation that the curve-fit audit memory rule was designed to catch. The discipline worked: task #26 is exactly the audit, even if dispatched as a "literature grounding" task.

2. **The correct spinor-on-cone paper is Fursaev-Miele 1996 (hep-th/9605153)**, not Fursaev-Solodukhin 1995. Distinct authorships, distinct papers.

3. **The published Fursaev-Miele 1996 says spin 1/2 "resembles the scalar case"** — antisymmetric Cheeger-like, with no Möbius modification. The discrete substrate Möbius α/(2α-1) is therefore either a novel continuum effect or a discrete artifact; mechanism identification remains OPEN.

## 8. Honest scope

- **Closed at sub-2%:** the empirical Möbius match (task #25) survives this audit. The numerical observation is unchanged.
- **Closed-with-revision:** the mechanism attribution. The "Fursaev-Solodukhin spinor double-cover correction" claim is REMOVED. Mechanism is open.
- **Not closed at theorem-grade:** none of the three interpretations (4.1, 4.2, 4.3) of what the discrete substrate computes is established. Interpretation 4.2 (novel continuum effect) is the most consistent with task #25's data; interpretation 4.1 (UV artifact) is partially ruled out by N_0-independence.
- **PDF-level verification of Fursaev-Miele 1996 not attempted.** The "resemble the scalar case" abstract quote is from search results, not from the paper's equations. A direct reading of Fursaev-Miele §III (the spin-1/2 calculation section) is the natural next-step if the mechanism question is reopened.

## 9. Files

- `debug/fursaev_solodukhin_1995_grounding_memo.md` (this, ~2500 words)

No driver, no JSON, no code modifications. Pure literature work.

## 10. What's next in the queue

Task #27 (geometric-mean azimuthal discretization) unblocks. Note that task #28 (wedge-spectral-density heat-kernel expansion) is now elevated in importance: it could provide the missing theoretical derivation that would resolve the mechanism question opened by this task. If task #28 finds the Möbius emerges naturally from the wedge spectral density at excess angle, interpretation 4.2 (novel continuum effect) would be supported and the v3.19.0 result would gain theoretical backing without the falsified Fursaev-Solodukhin citation.

## 11. Cross-references

- v3.19.0 Track 5 original (with the falsified citation): `debug/alpha_gt_1_analytical_investigation_memo.md` §2.3, §3.2, §8
- v3.19.0 ansatz sweep: `debug/alpha_gt_1_ansatz_test.py` (lines 112-124 reference "Fursaev-Solodukhin 1995 Phys.Lett. B365, 51, hep-th/9512134" — both wrong)
- Task #25 (empirical lock at sub-2%): `debug/alpha_gt_1_moebius_validation_4_5_10_memo.md`
- Paper 51 §12.7.7 (needs revision per §6.1): `papers/group5_qed_gauge/paper_51_gravity_arc.tex`
- Memory file (needs revision per §6.5): `memory/alpha_gt_1_moebius_closed_form.md`
- Curve-fit audit rule: `memory/feedback_audit_numerical_claims.md`
