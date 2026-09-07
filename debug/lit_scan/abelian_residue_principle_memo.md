# Adversarial external-literature scan: Paper 58 Theorem 1 ("the abelian residue")

**Date:** 2026-09-06 · **Scope:** is Paper 58's Thm 1 (`papers/group2_quantum_chemistry/paper_58_abelian_residue.tex`, L~594-660) a NAMED, generally-stated principle in the literature, or folklore restated ad hoc?
**Method:** ~14 targeted web searches + 6 direct fetches. Every identifier below was returned by a search or fetched. Anything not confirmed is marked **UNVERIFIED**.

---

## Decomposition of the target claim

Thm 1 is a conjunction of four separable pieces. Prior art differs sharply by piece.

| # | Piece | Prior-art status |
|:--|:------|:-----------------|
| P1 | Under G to H, the surviving exact labels are the irreps/characters of H | **STANDARD since 1929** (branching / subduction / "descent in symmetry") |
| P2 | The residual group of superposed constraints is the **intersection** of their symmetry groups | **NAMED PRINCIPLE since 1894** (Curie) |
| P3 | Abelian residue implies 1-dim irreps implies additive integer charge (Pontryagin) | Elementary textbook rep. theory; not a *principle* |
| P4 | **The asymmetry**: abelian selection rules survive, non-abelian CG selection does not, stated generally | **NOT FOUND as a general statement.** Restated repeatedly as an *engineering fact* per subfield |

---

## Ranked findings

### 1. Curie's principle (1894) — **COLLISION with P2**

IUCr Teaching Pamphlet No. 18, *An Introduction to Crystal Physics*, section 8 "Curie's Principle" — http://ww1.iucr.org/iucr-top/comm/cteach/pamphlets/18/node8.html (fetched verbatim; **author name beyond pamphlet number 18 UNVERIFIED**).

> "The symmetry group of a crystal under an external influence is given by the greatest common subgroup of the symmetry group of the crystal without the influence and of the symmetry group of the external influence." Formally **K-tilde = K ∩ G**.

The intersection mechanism in Thm 1(ii)/(iii) — "the intersection of the two centers' rotation groups is exactly the axial U(1)" — is Curie's principle applied to two Coulomb centers. This is a **named, general, 132-year-old principle**. The page also notes Curie's rule is a special case of the Shubnikov-Koptsik "principle of superposition of symmetry groups" (**Shubnikov & Koptsik, *Symmetry in Science and Art*, ISBN 9780306307591; page/statement UNVERIFIED**).

### 2. Bethe 1929, descent in symmetry — **COLLISION with P1**

H. Bethe, "Termaufspaltung in Kristallen," *Ann. Phys.* **395**, 133 (1929), **DOI 10.1002/andp.19293950202** (verified).

First systematic treatment of what happens to atomic terms when SO(3) is reduced to a crystal point group: the surviving labels are the subgroup's irreps, degeneracies split per the branching. Origin of the correlation-table / "descent in symmetry" apparatus. P1 is textbook, not new.

Modern restatement, verified: B. Klee, "Many Correlation Tables are Molien Sequences," **arXiv:1602.08774** (fetched) — correlation tables determine eigenstate splitting when a perturbation breaks Hamiltonian symmetry from a group to a subgroup.

### 3. Mulliken correlation diagrams — **ADJACENT (the diatomic case, primary source)**

R. S. Mulliken, "The Interpretation of Band Spectra Part III," *Rev. Mod. Phys.* **4**, 1 (1932), **DOI 10.1103/RevModPhys.4.1** (verified). Parts I/IIa/IIb: *Rev. Mod. Phys.* **2**, 60 (1930), **DOI 10.1103/RevModPhys.2.60** (verified).

The united-atom/separated-atom correlation diagram: exactly Paper 58's own "correlation-diagram reading," at its primary source. Paper 58 currently attributes this to Herzberg 1950 and cites `mulliken1949`; the 1930/1932 RMP series is the actual origin.

### 4. Textbook diatomic folklore — **ADJACENT (P1+P2+P3 for one case, unnamed)**

Every molecular-physics text states: in a diatomic the potential has only cylindrical symmetry, so l is not a good quantum number but the axial projection Lambda = |m_l| is. Sources surfaced were course notes / LibreTexts (no DOI — correctly, because it is folklore). This is Thm 1(i)+(ii) for the specific case, stated without any general principle behind it. **This is the single most important honest point: the diatomic residue is universally known, and universally stated ad hoc.**

### 5. Michel, symmetry breaking / isotropy subgroups — **BACKGROUND**

L. Michel, "Symmetry defects and broken symmetry. Configurations. Hidden symmetry," *Rev. Mod. Phys.* **52**, 617 (1980), **DOI 10.1103/RevModPhys.52.617** (verified).

General theory of residual (isotropy) subgroups under symmetry breaking. Gives the *lattice of residual groups*, not a statement about which selection rules survive, and carries no abelian/non-abelian asymmetry.

### 6. Contemporary QC statements of the asymmetry — **ADJACENT to P4, in Paper 58's own domain** (recency guard)

- L. D. da Silva & M. P. Santos, "Lie-algebraic incompleteness of symmetry-adapted VQE for non-Abelian molecular point groups," **arXiv:2603.21009** (2026) (arXiv URL + INSPIRE record 3133028 verified; **author names taken from a search summary, abstract not fetched**). Symmetry-adapted VQE "effectively reduce[s] the parameter count for Abelian molecular point groups, yet systematically fail[s] for non-Abelian groups"; the abelian-subgroup restriction spuriously splits multidimensional irreps and discards cross-component excitations.
- A. Gandon, A. Baiardi, M. Rossmannek, W. Dobrautz, I. Tavernelli, "Quantum computing in spin-adapted representations for efficient simulations of spin systems," **arXiv:2412.14797** (fetched). "[E]fficiently accounting for non-Abelian symmetries, such as the SU(2) total-spin symmetry, remains a major challenge."
- Quantum-chemistry codes work in D2h and its abelian subgroups precisely because abelian irreps are 1-dimensional — *J. Phys. Chem. A* **130**(24), 4683 (**URL verified, DOI UNVERIFIED**).

These say *abelian is the tractable case* and *non-abelian is where methods break*. None states Thm 1's structural reason as a principle; all are method papers.

### 7. Bartlett, Rudolph & Spekkens — **BACKGROUND, and an honest negative**

*Rev. Mod. Phys.* **79**, 555 (2007), **DOI 10.1103/RevModPhys.79.555**, **arXiv:quant-ph/0610030** (both verified; abstract fetched).

Prime candidate for a general abelian/non-abelian superselection statement (superselection sectors, twirling over a group, noiseless subsystems). **I could only fetch the abstract and did NOT confirm any such statement in the body.** Treat as an unresolved lead, not as support. If this scan is ever extended, this is the first place to look.

### 8. Pontryagin duality in physics — **ADJACENT (same shape, other setting)**

The Brillouin torus is the Pontryagin dual of the translation lattice (nLab "Brillouin torus", https://ncatlab.org/nlab/show/Brillouin+torus, URL verified). Continuous translations broken to a lattice; crystal momentum k is the surviving character label — structurally identical to Thm 1(iii), in a different domain, and again not stated as a general principle.

### 9. Residual symmetries in selection rules, 2026 — **BACKGROUND, mild counterpoint**

J. Dong, T. Kobayashi, S. Miyamoto, R. Nishida, H. Otsuka, "Residual group-like symmetries in selection rules without group actions," **arXiv:2603.14836** (fetched). Non-invertible selection rules in heterotic string compactifications: "residual group-like symmetries, **including both Abelian and non-Abelian ones**, remain exact" under groupification. Different setting, but note it explicitly declines the abelian/non-abelian split.

### 10. Checked and empty (negative results, reported plainly)

- **No named principle** of the form "abelian selection rules are structurally more robust under subgroup restriction." Searched: anomaly / 't Hooft matching, superselection sectors, abelian-vs-non-abelian anyon robustness, additive-vs-vector-coupling quantum numbers, Wigner-Eckart under subgroup restriction, "descent in symmetry" as a general rule. **Negative.**
- Paschen-Back (j dies, m_j = m_l + m_s survives) and the Stark effect are the same residue in atomic physics and are, again, textbook-ad-hoc. No general-principle source found; no citable identifier pursued, and none should be invented.

---

## OVERALL NOVELTY VERDICT

**The general statement is NOT novel; it is two standard principles composed.** Specifically:

- **P1 is Bethe / branching theory (1929).** Not ours.
- **P2 is Curie's principle (1894), a named general principle with the exact intersection form K ∩ G.** Not ours. This is the sharpest hit in the scan.
- **P3 is elementary representation theory** (abelian implies 1-dim irreps; Pontryagin dual of a compact abelian group is discrete). Not ours.
- **P4 — the explicit abelian/non-abelian asymmetry, stated as a general principle — was NOT found in the literature.** It is restated per-subfield as an engineering fact (QC tapering, D2h-only quantum chemistry, SU(2) as a "major challenge") and as case-specific folklore (Lambda in diatomics, m_j in Paschen-Back, k in solids), never as a principle with a stated mechanism.

**What is genuinely ours:** not the principle, but (a) the *composition* — routing Curie-intersection plus branching through Pontryagin duality to name the residue as an additive integer charge; and (b) the **application target**: Paper 58 uses the residue to predict the *sparsity structure of the two-center integral tensor* (which l-selection dies, which m-selection survives, and the measured 13.8x-15.0x permitted-density inflation), where the entire prior literature uses symmetry descent to predict *spectroscopic level splitting*. That re-targeting — symmetry descent as a statement about matrix-element density rather than about term splitting — is the defensible contribution.

**Hostile-referee exposure (flag to PI).** A referee who knows Curie's principle will say Thm 1(ii)/(iii) is Curie plus branching. Worse: Paper 58's own paragraph immediately after the theorem concedes that the largest group fixing both nuclei is the **non-abelian** C_inf,v and that it *does* yield usable label structure. So the "abelian residue" headline is narrower than its name, and the paper says so itself one paragraph later. Recommendation: keep the [INTERNAL THEOREM] tier, add the prior art below, and soften any implication that the abelian/non-abelian asymmetry is being asserted as new.

---

## Best citations to add to Paper 58

1. **Curie's principle** — the intersection rule, named and general. Cite via IUCr Teaching Pamphlet 18 section 8 (URL above) and/or Shubnikov & Koptsik. **Highest priority**; it names Thm 1's own mechanism.
2. **Bethe (1929), DOI 10.1002/andp.19293950202** — descent in symmetry / branching as the origin of "which labels survive."
3. **Mulliken, Rev. Mod. Phys. 4, 1 (1932), DOI 10.1103/RevModPhys.4.1** (and RMP 2, 60 (1930), DOI 10.1103/RevModPhys.2.60) — primary source for the correlation diagram the paper already invokes via Herzberg.
4. **arXiv:2603.21009** and **arXiv:2412.14797** — recency guard: contemporary statements, in Paper 58's own domain, that abelian symmetry is the exploitable case and non-abelian is where methods fail. Omitting these is the most likely "you did not cite the obvious" hit.

*Verification note: identifiers marked UNVERIFIED were surfaced by search but not independently fetched. Do not promote them into the paper bibliography without confirming.*
