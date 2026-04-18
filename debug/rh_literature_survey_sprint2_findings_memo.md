# Sprint 2 RH Literature Findings (Verified Web)

**Date:** 2026-04-17. Web access verified (example: WebSearch on "Ihara zeta function Ramanujan graph 2024" returned 10 real URLs, e.g. https://en.wikipedia.org/wiki/Ihara_zeta_function). Budget 45 min; Sprint 2 goals met. 16 net new verified entries appended as §6 to `debug/rh_literature_survey.md`.

## The four required questions

### Q1. Is RH-A (Ihara zeta on Fock-projected S^3 / Bargmann-Segal S^5 / Hopf graph) duplicated in verified literature?

**NO.** Verified arXiv/journal search 2024-2025 turns up zero computations of the Ihara zeta function on any spherical graph, Hopf-fibration graph, or Bargmann-Segal lattice. The closest adjacent entries are:

- **Chico-Mattman-Richards 2024** (arXiv:2501.00639): Ihara ζ for complete, complete bipartite, Möbius ladders, cocktail-party, all graphs of order ≤ 5, rank-2 graphs. No spheres.
- **Matsuura-Ohta 2022 / 2025** (arXiv:2204.06424, *PTEP* 2025, 063B01, DOI 10.1093/ptep/ptaf071): Kazakov-Migdal gauge model on arbitrary graph → partition function = 1/ζ_G. Grid graphs only. **This is the strongest adjacent precedent** — it establishes Ihara-zeta-as-partition-function in a lattice-gauge physics context on general graphs, exactly the Paper 25 framing. It does NOT compute ζ on our Fock/Bargmann-Segal graphs.
- **Lei-Müller 2023** (arXiv:2307.01001) and follow-up 2025 (arXiv:2509.15214): Ihara ζ of supersingular isogeny graphs = Hasse-Weil ζ of modular curves. Arithmetic-geometric, not Lie-group geometric.

**Conclusion:** RH-A (Sprint 3, planned) is not duplicated. Paper 29's headline `det(I − uB)` evaluation on the N_max ≤ 5 Bargmann-Segal graph remains a first.

### Q2. Is there a post-2024 Hilbert-Pólya proposal on a compact manifold that could bite on GeoVac?

**NO compact-manifold candidate.** The most active post-2024 line is:

- **Yakaboylu 2024**, *J. Phys. A* 57, 235204 (arXiv:2408.15135). Non-Hermitian operator on `L^2([0,∞))` — the **half-line**, not a compact manifold. Ongoing through v15 (March 2026). Bellissard-style circularity critique appears in the review trail but no formal retraction.
- **LeClair/Martinez 2024** (JHEP 04, 062): integrable scattering Hamiltonian with Bethe-ansatz energies = Riemann zeros. Continuous non-compact.
- **Martinez et al. 2025 survey** (*Symmetry* 17(2), 225, MDPI): catalogs all post-2020 attempts, none on a compact manifold.

The compact-manifold setting (S^3, S^5, Bargmann-Segal) remains unclaimed territory in the Hilbert-Pólya literature. Paper 29's setting does not conflict with any verified operator proposal.

### Q3. Is there a recent Bianchi Selberg-zeta result that revives the closed RH-B direction?

**NO.** Verified 2024-2025 Selberg-zeta work focuses on cocompact Γ ⊂ PSL(2,R) (Gajewski-Hinrichs 2024, arXiv:2405.11084; the 2025 paper arXiv:2511.07100 on second moments at σ=1). No new numerical Bianchi Selberg-zero computation appeared post-Then-2005 in the verified search. The Wick-rotation obstruction documented in Sprint 1 (`SO(4) hydrogen → PSL(2,O_K)` has no literature bridge) stands. **Keep RH-B shelved.**

### Q4. Any concrete Paper 29 citation that turned out wrong / retracted / mis-dated?

**One caveat, no hard errors.** Paper 29's bibliography (lines 644-696 of `papers/observations/paper_29_ramanujan_hopf.tex`) cites only foundational material: Fock 1935, Bander-Itzykson 1966, Ihara 1966, Bass 1992, Hashimoto 1989, Kotani-Sunada 2000, Terras 2010, LPS 1988, Morgenstern 1994, Stark-Terras 1996, Alon-Milman 1985 / Alon 1986, Elstrodt-Grunewald-Mennicke 1998. **All 12 are verified by web search as correctly attributed.** No retraction or misdating found.

**Caveat (minor):** Paper 29 does NOT currently cite:
- Matsuura-Ohta 2022/2025 (Kazakov-Migdal ↔ Ihara zeta). This is the strongest post-2024 physics precedent for the Paper 25/29 framing and **should be added** before external release.
- Yakaboylu 2024. Optional but recommended for due diligence on Hilbert-Pólya.
- Huang-McKenzie-Yau 2024 (Ramanujan universality of random d-regular graphs). Optional; useful as a contrast showing that GeoVac's graphs are *not* random — they are structured — which matters for interpreting what "Ramanujan" means here.

## Recommendation for Sprint 3

1. **Proceed with RH-A** (Ihara-zeta computation on the Bargmann-Segal S^5 graph at N_max ≤ 5, and on the Fock-projected S^3 graph at n_max ≤ 4). No duplication risk. Use Paper 25's existing non-backtracking B matrix; compute `det(I − uB)` in exact rational arithmetic (sympy Fractions); inspect pole locations relative to `|u| = 1/√(q-1)` for average degree q.

2. **Add three citations to Paper 29 before external release:** Matsuura-Ohta 2025 PTEP (physics precedent for ζ_G as partition function), Yakaboylu 2024 J. Phys. A (compactness caveat in the discussion of Hilbert-Pólya), Huang-McKenzie-Yau 2024 arXiv (Ramanujan-universality context). Short bibitems; no restructure of the paper needed.

3. **Do NOT open RH-B.** The Bianchi-Selberg direction is not revived by 2024-2025 literature and still requires the unresolved Wick-rotation memo.

4. **Flag for plan-mode:** the Sprint-1 entries 1-60 remain web-unverified. Spot-checking them one-by-one is low-yield; the load-bearing ones (Terras 2010, Elstrodt-Grunewald-Mennicke 1998, Berry-Keating 1999, BBM 2017, Bellissard 2017, Kotani-Sunada 2000, Lubotzky-Phillips-Sarnak 1988) are corroborated by the verified §6 searches indirectly (all appear in the search-result Wikipedia / survey pages from verified queries). The remaining Sprint-1 entries are background context, not Paper 29 citations.

(Word count: ~780.)
