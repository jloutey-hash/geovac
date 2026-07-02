# Project Overview Survey — 2026-07-01 (PI-requested)

Three parallel surveys (sonnet, per PI direction 2026-07-01): code quality, papers/QA state,
untaken research directions. Condensed results; PM synthesis delivered in-session.
Not a sprint memo — a state-of-project snapshot.

---

## 1. Code quality survey

**Inventory:** 178 modules in `geovac/` (excl. `_archive`; incl. `nuclear/` 11, `gravity/` 2),
143,263 LOC. Largest: `composed_qubit.py` (3,076), `hylleraas_r12.py` (2,687),
`n_electron_solver.py` (2,648), `hamiltonian.py` (2,531), `tannakian.py` (2,426).
Tests: 314 files, 8,372 collected, 519 slow-marked. Type-hint coverage 5,571/5,653 functions (99%).

**Health:** all sampled physics suites PASS — 18 S³ symbolic proofs (6s), P38 keystone,
paper suites (7/8/11/14/20/22/1), trunk QA claims, Casimir CI, algebraic angular. C14 gate PASS.

**Five concrete risks:**
1. `tests/test_composed_qubit.py:330` — `lih_result_n4` session fixture (max_n=4 build, 334s)
   NOT slow-marked; adds ~5–10 min to every plain `pytest` run.
2. Constant precision drift, no central module: `C_LIGHT` = 137.036 (`dirac_hamiltonian.py`)
   vs 137.035999084 (`hamiltonian.py:1484`) vs 137.035999177 CODATA-2022 (`hopf_bundle.py:414`);
   `HARTREE_TO_EV` duplicated in two modules.
3. `docs/code_architecture.md` — 3 stale module refs (`prolate_spheroidal.py`,
   `hyperspherical.py`, `nested_hyperspherical.py` → actual names differ / archived).
4. pytest rootdir collects sub-projects: `ADSCFT/` (3 FAIL) + `geovac-hamiltonians/`
   (1 import ERROR) pollute the canonical fast gate. Fix: `--ignore` entries in `pytest.ini` addopts.
5. Wildcard imports in `geovac/hyperspherical_coupling.py` and `geovac/pauli_projector.py`
   (the former on the active Level 3–4 solver path).

**Grade:** solid research-grade codebase; hygiene metrics excellent (4 TODO/FIXME in 143k LOC);
the two urgent items are the unmarked slow fixture and the rootdir pollution.

---

## 2. Papers / QA survey

**Corpus:** ~52 active papers across 6 groups + 3 synthesis + 7 archived.
Status flags: P45/46 DESCOPED, P47–49 PARTIAL, P52/53 DRAFT, 2 GUARDRAIL, 3 OBSERVATION (group5).

**Certification scoreboard:** trunk (v4.16.2), group3 (v4.21.1), group1 (v4.49.0),
group2 (v4.51.0) CERTIFIED — 40 of ~52 active papers. group4: 5 FAIL→remediated runs,
6th (certifying) run HELD per PI. group5, group6, synthesis: NOT STARTED (no DoD files;
zero rows in claim_test_matrix).

**Claims register:** 21 headline claims — 3 SYMBOLIC PROOF, 6 INTERNAL THEOREM, 7 MEASURED,
1 PANEL-VERIFIED, 2 OBSERVATION, 2 RETRACTED (preserved for DOI resolution).
Claim #17 (Mellin case-exhaustion) covers 15/28 Paper 34 projections.

**Top debt:**
1. ~443 dangling `debug/` refs (v4.21.0 C14 sweep; papers 34/32/28) — deferred, no reduction since.
2. group5/group6 never QA'd, matrix empty.
3. group4 certifying run held (PI timing).
4. Missing PDFs: `paper_33_qed_selection_rules`, `paper_26_entanglement`,
   `paper_27_entropy_projection`; `paper_34` PDF one day staler than its .tex.
5. Known NO-TEST gaps: P18 α²-Ihara + κ-B; P54 selection rules; P55 mixed-Tate (M2 leg only);
   P56 PS-4; P57 ×2; group4 7/14 gaps open.

---

## 3. Untaken research directions (48 threads, 8 themes)

Status codes: (a) explicitly queued, (b) named follow-on never scheduled, (c) open problem in paper.

**Precision/atomic (8):**
1. §V.C.2 H 21cm four-component autopsy (a, CLAUDE.md §1.8)
2. §V.C.3 muonic-H Lamb autopsy (a, §1.8)
3. §V.C.4 He 2³P fine-structure autopsy (a, §1.8)
4. §V.C.5 He 2¹P→1¹S oscillator strength (a, §1.8)
5. §V.C.6 Cs 6S₁/₂ hyperfine + PNC prospective (a, §1.8)
6. LS-8a/8b/8c multi-loop QED (b; LS-8a WEAK at renormalization gap; Paper 36 §Outlook)
7. LS-8d recoil+nuclear-size, LS-8e hyperfine averaging (~+4.4 MHz residual) (b, Paper 36)
8. W1e cross-system robustness (LiH n_max=3, MgH₂+F5) (a, CHANGELOG v3.95.0)

**Math.OA / Lorentzian NCG (9):**
9. Strong-form Lorentzian propinquity, Sprint L3e Q1 (c) — **PM note: the 2026-06-19
   compact-boost structural closure resolved P45 §open Q1 on the convention branch;
   remaining honest legs are WH7 B2 pointed-proper bookkeeping + the state/modular boost
   candidate, NOT a metric-level reopen. Agent's #1 leverage rank overstates current standing.**
10. G2-metric de-compactification T→∞ / pointed-propinquity bridge (c, group1 synthesis)
11. G3 cross-manifold T_S³⊗T_Hardy(S⁵) (c, blocked by 4-layer asymmetry)
12. **G4a Connes SM / Higgs inner fluctuation (a)** — scoped POSITIVE-THIN 6–7 wk;
    prerequisite Paper 38 now PROVEN; venue Paper 32 §VIII.D addendum
13. Wasserstein–Kantorovich metric on K⁺ Krein states (c, Paper 44)
14. Cocycle entropy-production deficit closed form (b, Paper 49 §10 Q1, 2–4 wk)
15. Q2' non-commutative Mondino–Sämann pre-length space (b, 3–12 mo)
16. F₄/E₆/E₈ rank-uniform interior-leg formal proof (b, Paper 43 §11 O4)
17. WH6 Dirac spectral zeta D(s) (paused 2026-04-18, not closed)

**Cosmic-Galois / periods (6):**
18. S_min^(S⁵) PSLQ closed form (a, 2–3 wk)
19. PSLQ Sym² periods vs Hain–Brown modular periods (a, 1–2 wk; gates P55/56 framing)
20. Tannakian converse Aut⊗(ω)=U*_Levi (c, multi-year)
21. Kleinschmidt depth-2 coaction module closure, 3 paths (a)
22. Inner-factor Yukawa Dirichlet ring mixed-Tate check (b, Paper 55 §7 exclusion)
23. GH rate constant b≈4.1093 PSLQ (c, Paper 38 §Open i)

**Quantum simulation / hardware (6):**
24. Phase-3 program: resource estimation → circuit compilation → experimental collab (b, §1.6)
25. M-vS-3 two-body ERI Bratteli reading (a, 3–4 wk; gates M-vS-4/5)
26. Hopf-Z2 tapering wiring into `ecosystem_export.hamiltonian()` (a, infra)
27. Z2 taper spectrum verification at Q=30 (31/37 systems unverified) (a)
28. J²-adapted QPT for relativistic Tier 2 (b, ~3–5 days)
29. M-vS-4/5 theorem formalization + companion paper (b, multi-month)

**Chemistry wall (5):**
30. W1e structural proof at projection-step level (b)
31. Explicit-core CCSD/CCSD(T) closure attempt (b, conditional, 3–6 mo, paused at diag gate)
32. Frozen-core V_cross R-dependence in `balanced_coupled.py` (a, sprint-scale, since v3.56.0)
33. Algebraic-implicit form for adiabatic μ(R) (c, group2 synthesis)
34. Second-row/polyatomic balanced extension (c, Paper 17 boundaries)

**Nuclear / cross-domain (4):**
35. Gravity arc G4-4/5/6 → discrete-substrate S_BH=A/4 (b, 6–12 mo; G4-3 done, Paper 51)
36. G4-6c Möbius continuum theorem, Route A (b, ~1 wk alt-discretization check)
37. Two-body Coulomb long-range NCG route (b, multi-year, category-mismatch wall)
38. Bargmann–Segal HO graph-Laplacian convergence theorem (c, Paper 23 §9)

**Infrastructure / papers (5):**
39. **Paper 38 companion manuscript** (b, 4–8 wk pure writing; L1'–L5 memos complete)
40. Corpus citation sweep — 18 identified mechanical corrections (a, single sprint)
41. Paper 32 §VIII.D + Paper 35 §VIII.A BW paragraphs, Track B (a)
42. Outreach van Suijlekom/Marcolli + Perez-Sanchez co-citation note (b, PI-timed)
43. Higher-Q VQE scaling + sector-projected UCCSD on tapered LiH (b)

**Foundations open problems (5):**
44. α combination-rule structural origin (c, WH5; 12 mechanisms eliminated, no handle)
45. Second packing axiom / W3 calibration-data principle (c, multi-year, no handle)
46. Lorentz boost as 26th Paper 34 projection (c, §VIII)
47. Paper 40 rank r≥2 log-power at higher cutoff (c)
48. Higgs direction n̂∈S² ↔ Hopf-base S² identification (b, 1–2 wk)

**PM additions — directions never named anywhere in the repo:**
- **Real-time dynamics / spectra on the graph.** The corpus is static (eigenvalues,
  ground states, a few oscillator strengths). Chebyshev/Krylov propagation exploits O(V)
  sparsity directly — absorption spectra, autocorrelation, response functions from the same
  sparse Hamiltonians; conceptual hook to WH7 (the observer's compact time window).
- **Machine-checked formalization (Lean) of the 18 S³ proofs and/or the Paper 38 rate.**
  The viability case rests on "verified structural results"; machine-checked proofs are the
  strongest possible form for an independent, non-journal corpus.
- (Weaker candidates, flagged not endorsed: CTQW on the GeoVac lattice as a native
  quantum-simulation primitive; periodic/solid-state composed extension — likely re-hits
  the cross-block walls.)
