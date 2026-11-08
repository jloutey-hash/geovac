# Citation & Novelty Check: Paper 7 — *The Dimensionless Vacuum: Recovering the Schrödinger Equation from Scale-Invariant Graph Topology*

**Date checked:** 2026-06-01
**Verdict:** YELLOW (one MEDIUM CITE-DOESNT-SUPPORT, one MEDIUM novelty softening, one MEDIUM CITE-WRONG-METADATA on Bander–Itzykson)

---

## Citation table

| `\cite` key | claimed as | verdict | what I found |
|---|---|---|---|
| `fock1935` | Fock, *Zur Theorie des Wasserstoffatoms*, Z. Phys. **98**, 145–154 (1935). Used for: stereographic projection of momentum-space hydrogen onto $S^3$ (load-bearing). | **CITE-OK** | DOI [10.1007/BF01336904](https://link.springer.com/article/10.1007/BF01336904); NASA/ADS [1935ZPhy...98..145F](https://ui.adsabs.harvard.edu/abs/1935ZPhy...98..145F). Title, vol, pages, year all correct. Content matches use. |
| `bargmann1936` | Bargmann, *Zur Theorie des Wasserstoffatoms: Bemerkungen…*, Z. Phys. **99**, 576–582 (1936). Used for: hidden SO(4) ↔ isometry group of $S^3$. | **CITE-OK** | DOI [10.1007/BF01338811](https://link.springer.com/article/10.1007/BF01338811). Title, vol, pages, year all correct. Bargmann did make the Pauli↔Fock / SO(4) connection in this paper. |
| `barut1967` | Barut & Kleinert, *Transition probabilities of the hydrogen atom from noncompact dynamical groups*, Phys. Rev. **156**, 1541–1545 (1967). Used (line 35) for: "extended the symmetry to the full conformal group SO(4,2)"; also (line 99) for: SU(2)⊗SU(1,1) factorization of SO(4,2). | **CITE-DOESNT-SUPPORT** (medium) | Paper exists and the metadata is correct ([APS link](https://link.aps.org/doi/10.1103/PhysRev.156.1541)). But this is the **SO(4,1)** noncompact-dynamical-group paper, not the SO(4,2) paper. The same authors extended the algebra to **SO(4,2)** in *Current Operators and Majorana Equation for the Hydrogen Atom from Dynamical Groups*, Phys. Rev. **157**, 1180 (1967). Cf. review by Gallup et al. ([Symmetry **12** (8), 1323 (2020)](https://www.mdpi.com/2073-8994/12/8/1323)): Bander–Itzykson 1966 introduced SO(4,1); the SO(4,2) extension came later. The SU(2)⊗SU(1,1) factorization (line 99) is also not the standard decomposition discussed in this specific paper — Barut–Kleinert 1967 use parabolic-coordinate ladders, not an explicit SU(2)⊗SU(1,1) maximal-compact×noncompact split. |
| `bander1966` | Bander & Itzykson, *Group theory and the hydrogen atom (I)*, Rev. Mod. Phys. **38**, 330–345 (1966). Used (line 35) as evidence for "dynamical symmetry groups." | **CITE-WRONG-METADATA** (low/medium) | Paper exists ([APS](https://link.aps.org/doi/10.1103/RevModPhys.38.330)). Authors, journal, vol, year, page all correct. However, Bander–Itzykson published a two-part series — Paper I (pp 330–345, SO(4) bound-state) and Paper II (pp 346–358, SO(3,1) scattering). The cite refers only to Paper I, which is correct for the bound-state / SO(4) claim. The companion review [arXiv:2305.18229](https://arxiv.org/abs/2305.18229) notes that Bander–Itzykson introduced **SO(4,1)** as the spectrum-generating algebra; if Paper 7 means to attribute SO(4,1) here it is fine, but the sentence on line 35 lumps `bander1966` with `biedenharn1981` as "thoroughly studied in the context of dynamical symmetry groups" — that is benign. The metadata issue is just that no part numeral appears in the bibitem; it should read "(I)" for clarity. |
| `biedenharn1981` | Biedenharn & Louck, *Angular Momentum in Quantum Physics*, Encyclopedia of Mathematics and its Applications **Vol. 8** (Addison-Wesley, 1981). Used as a general reference for the SO(4)/conformal-group hydrogen literature. | **CITE-OK** | Verified: ISBN 0521302285, Cambridge / Addison-Wesley, 1981, Vol. 8 of EMA series. Title, authors, year, publisher, volume number all correct. |
| `loutey_paper0`, `loutey_paper1`, `loutey_paper11`, `loutey_paper13`, `loutey_paper2` | Internal self-citations to GeoVac papers. | **CITE-OK** | Files exist at paths listed in CLAUDE.md §6 (e.g. `papers/group3_foundations/Paper_0_Geometric_Packing.tex`, `paper_1_spectrum.tex`, `papers/group2_quantum_chemistry/paper_11_prolate_spheroidal.tex`, `paper_13_hyperspherical.tex`, `papers/group5_qed_gauge/paper_2_alpha.tex`). Internal consistency check only; cannot externally verify. |

---

## Problems found

### MEDIUM — `barut1967` (CITE-DOESNT-SUPPORT for the SO(4,2) attribution)

**Location.** Two affected sentences:

1. **§I, line 35:** "…by Barut and Kleinert \cite{barut1967}, who extended the symmetry to the full conformal group $SO(4,2)$."
2. **§III, line 99 (Eq. 4 area):** "These operators form the algebraic skeleton of the $SU(2) \otimes SU(1,1)$ dynamical symmetry group, which is the maximal compact $\times$ non-compact factorization of the hydrogen conformal group $SO(4,2)$ \cite{barut1967}."

**Problem.** Both claims are true at the level of the *literature*, but the **specific 1967 Barut–Kleinert paper at Phys. Rev. 156, 1541** does not contain the SO(4,2) extension. That paper introduces and uses the **non-compact dynamical group SO(4,1)** (spectrum-generating algebra) to compute transition probabilities. The extension to **SO(4,2)** comes in the same authors' very-next paper Phys. Rev. **157**, 1180 (1967), *Current Operators and Majorana Equation for the Hydrogen Atom from Dynamical Groups*, where adding the dipole operator promotes SO(4,1) → SO(4,2). (Independent confirmation: Gallup et al., *Symmetry* **12**, 1323 (2020), arXiv:2305.18229.)

**Recommendation.**
- For line 35: either (a) replace the cite with Barut & Kleinert, Phys. Rev. **157**, 1180 (1967), or (b) cite both, or (c) reword to "Barut and collaborators (1966–67), who developed the noncompact dynamical groups SO(4,1) and SO(4,2)."
- For line 99: the SU(2)⊗SU(1,1) factorization is not the headline content of either Barut–Kleinert 1967 paper. Better citation: Wybourne, *Classical Groups for Physicists* (Wiley, 1974), or Adams, Čížek & Paldus, *Adv. Quantum Chem.* **19**, 1 (1988), both of which lay out the SU(2)⊗SU(1,1) ⊂ SO(4,2) decomposition explicitly. The Biedenharn–Louck reference would also support the SU(2)⊗SU(1,1) claim (Ch. 6 of vol. 8). At minimum, add `biedenharn1981` alongside `barut1967` at line 99.

### LOW — `bander1966` part-number cosmetic

The cite is for Paper I (Rev. Mod. Phys. **38**, 330). Bibitem text would benefit from the explicit "(I)" suffix to distinguish from Paper II (pp 346). Not load-bearing for any specific claim.

---

## Priority / novelty claims

| claim (verbatim) | location | searched | prior art found? | recommendation |
|---|---|---|---|---|
| "The mathematical structure has been thoroughly studied in the context of dynamical symmetry groups." | §I, line 35 | n/a (acknowledging prior art, not claiming novelty) | n/a | OK — this is a *concession* to prior work, not a novelty claim. |
| Sec. *What is New vs. What is Known* §II.B (4 enumerated contributions): (1) discrete graph convergence to $\Delta_{S^3}$; (2) 18 symbolic proofs in sympy; (3) "scale-invariance framing… not been stated explicitly in the literature"; (4) sparse $O(V)$ computational methodology. | §II.B, lines 76–84 | Searched: graph discretization of Laplace–Beltrami, hydrogen $S^3$ continuum-limit work, Coulomb Sturmians in momentum space, Shibuya–Wulfman tradition. | (1) and (2): no direct prior art for *this specific* discrete-graph→$\Delta_{S^3}$ convergence story with sympy verification. Generic graph→Laplace–Beltrami spectral convergence is well-known ([arXiv:1301.2222](https://arxiv.org/abs/1301.2222), and Belkin–Niyogi line of work) — paper does not claim novelty there, only for the hydrogen-specific instance. (3) is the load-bearing novelty claim — see below. (4) is computational, not a publication-priority claim. | Items (1), (2), (4) are appropriately scoped. Item (3) is the one to soften — see next row. |
| "While implicit in Fock's original treatment, this framing — and its implications for the GeoVac computational approach — has not been stated explicitly in the literature." | §II.B item 3, line 81 | Searched: "scale invariance hydrogen Fock projection unit $S^3$"; "dimensionless vacuum" framing; Avery's *Hyperspherical Harmonics and Generalized Sturmians* and the Shibuya–Wulfman / Coulomb-Sturmian tradition. | The *idea* that Fock's $S^3$ is intrinsically unit-radius regardless of $p_0$ is implicit in every textbook treatment (e.g. Biedenharn–Louck Ch. 6, Avery's hyperspherical-harmonics books). I did **not** find a prior paper that articulates "the conformal factor $\Omega$ always produces the unit $S^3$ regardless of $p_0$" as an *explicit organizing principle for a computational discretization*. So this is a defensible "to our knowledge, the first explicit statement" claim — but absence-of-evidence caveat applies (Avery's 1989/2000/2006 books would be the natural prior-art site to check). | **Soften to** "to our knowledge, the first explicit statement" rather than "not been stated explicitly in the literature." Already phrased close to this on line 730 ("has not, to our knowledge, been stated as an explicit organizing principle"), so line 81 should match that softer phrasing. |
| "We have verified Eq. ([eigenvalues]) symbolically for $n = 1, 2, 3$…" (and the broader 18-proof claim) | §IV.B, line 273; §I, line 46 | n/a — verification claim, not novelty | Internal-to-codebase claim; tests at `tests/test_fock_projection.py` and `tests/test_fock_laplacian.py` (cited in CLAUDE.md §4 as the 18 symbolic-proof suite). | OK — this is a reproducibility claim verifiable in-repo, not a novelty claim. |
| "The same combinatorial skeleton reappears in prolate spheroidal coordinates… hyperspherical coordinates… composed fiber bundles." (concluding sec.) | §VII, line 736 | n/a — internal cross-reference | Cites Papers 11, 13, 0 — consistent with CLAUDE.md §5. | OK. |
| Sec. V density-overlap formula $F^0(a,b) = (4Z/\pi)\int \Phi_a \Phi_b \,dt$ on $S^3$. No novelty claim is made, but the implication is novel-derivation. | §V, lines 391–397 (boxed eq.) | Searched: Coulomb Sturmians, Avery hyperspherical Slater integrals, Shibuya–Wulfman density. | The Sturmian/hyperspherical literature (Avery, Aquilanti, et al.) computes Slater integrals in momentum space via Fock projection; the *specific* single-integral density-overlap form may already exist in Avery's work. Paper 7 does not claim novelty for it. | Consider adding one sentence acknowledging Avery's Coulomb-Sturmian tradition as the natural prior-art context for Sec. V, with a cite to Avery, *Hyperspherical Harmonics: Applications in Quantum Theory* (Kluwer, 1989) or Avery & Avery, *Generalized Sturmians and Atomic Spectra* (World Scientific, 2006). This is courteous-citation hygiene, not a correctness issue. |

---

## Bottom line

**GREEN / YELLOW / RED → YELLOW.**

Nothing here would embarrass the paper in front of an expert, **but** the SO(4,2) attribution to Barut–Kleinert PR 156 (1967) is the kind of detail an expert in dynamical symmetry groups *will* notice and is straightforward to fix (the correct paper is by the same authors in the same year, PR 157, 1180). The novelty claim on line 81 should be softened to the "to our knowledge" phrasing already used on line 730 — the conclusion has the right hedge, but the §II.B contributions list does not.

### Required fixes (medium severity, do these)

1. **`barut1967` — split or replace.** At line 35, attribute the SO(4,2) extension to Barut–Kleinert, *Phys. Rev.* **157**, 1180 (1967), not PR 156. At line 99, the SU(2)⊗SU(1,1) factorization claim needs a different citation (Wybourne 1974 or Adams–Čížek–Paldus 1988 or Biedenharn–Louck 1981) — Barut–Kleinert PR 156 does not develop this decomposition.

2. **Line 81 softening.** Change "has not been stated explicitly in the literature" → "has not, to our knowledge, been stated explicitly in the literature" (mirroring line 730).

### Optional polish (low severity)

3. Add "(I)" to the `bander1966` bibitem title.
4. Add a one-sentence courtesy acknowledgment of the Avery Coulomb-Sturmian hyperspherical-harmonics tradition as the natural literature context for §V (the density-overlap formula).

### Strong points

- Fock 1935 and Bargmann 1936 — the two load-bearing historical citations for the entire paper — are correctly identified and used.
- Bander–Itzykson 1966 metadata is correct (modulo the part-numeral cosmetic).
- Biedenharn–Louck 1981 is correctly identified as Vol. 8 of the Encyclopedia of Mathematics and its Applications.
- All internal self-citations (`loutey_*`) point to extant files.
- The paper is appropriately careful about novelty in §II.B and the conclusion — most of the novelty work is in §V.C–V.D (density-overlap derivation + chordal-ansatz failure table), which is fully derived in-paper without external prior-art reliance.
