# Sprint: `/qa group1` Bite B sub-bite 1 (Papers 42, 43, 44 + synthesis) — 2026-06-18 (v4.23.0)

PI-invoked `/qa group1`, Bite B (Lorentzian cluster) as smaller sub-bites; this
is sub-bite 1 = **Papers 42 (four-witness modular), 43 (Lorentzian extension),
44 (BBB Krein operator-system)** + synthesis Lorentzian section.

## 1. Verdict: FAIL (calibrated) → remediated

10-agent panel (claims ×3, citation ×3, code ×3, synthesis ×1), all
path-pinned to worktree `../geovac-qa-seed-group1-B1`. **Calibrated:
sensitivity 5/5** (S-claims-C8 p42, S-claims-C14 p43, S-citation-C4 p44,
S-code-C2 p42, S-synthesis-C9 all caught), **specificity 5/5** (honest
open-scope controls + correct cites all affirmed SOUND). Seed key
`debug/qa/group1_B1_seed_key.json`. Worktree removed; no seed leaked.

## 2. THE HEADLINE: two reviewer LARGE findings OVERTURNED by PM verification

The §9 reconcile rule (PM verifies every MATERIAL finding vs primary code/text)
earned its keep twice — both on paper_42:

- **code-p42 "derived finding is FALSE" — WRONG; the finding is CORRECT.** The
  reviewer claimed the half-integer D_W also closes σ_2π, so the H_local=K_α/β
  "derived finding" is a false-positive. Verified against `geovac/
  modular_hamiltonian.py` (l.491–509): the modular generator is **β·H_local**,
  so H_local=D_W gives generator **2π·D_W** → σ_2π residual **3.58 (does NOT
  close)**; K_α closes (1.7e-15). The reviewer (and my first quick check) tested
  K=D_W *directly*, omitting the β=2π factor — the wrong object. **No reframe.**
- **claims-p42 "§10 re-asserts a descoped Lorentzian result" — OVER-FLAGGED.**
  §10 Theorem 10.2 is the *finite-cutoff Krein* four-witness **period closure**
  (σ_2π=id, backed by `test_modular_hamiltonian_lorentzian`, 137 pass) — the
  Lorentzian analog of the Riemannian §5, NOT the descoped Lorentzian-
  *propinquity convergence* (§10 itself names that the open target). Not a
  zombie. The real residual defect was narrower: a stale **intro disclaimer**
  (l.266) predating §10. Fix = reconcile the intro (continuum/propinquity open;
  finite-cutoff Krein closure done in §10), NOT descope §10.

Lesson reinforced: a fresh adversary finds plausible faults; PM verification
against primary code is what separates real from artifact. (Same shape as the
Bite-A Paper-40 resurrection.)

## 3. Genuine defects fixed

**Citations (verified vs arXiv + primary PDFs):**
- `hekkelman_mcdonald2024` ("Spectral truncations of T^d") **fabricated**
  (arXiv:2403.18619 = an OpenMP paper) — removed from **42/43/44**, tori
  re-attributed to Leimbach–vS. [same defect as Bite-A paper_40]
- `hekkelman2022` p42 wrong-ID (2206.13744 = Kerr-Melvin BH) → "Truncated
  geometry on the circle," LMP 112 (2022) 20, arXiv:2111.13865.
- `zhu_casini2020` p42 **fabricated authors** (Zhu/Casini/Hauke) → real
  Zhang–Calabrese–Dalmonte–Rajabpour; + 3 prose mentions fixed.
- `latremoliere2018` **42/43/44** wrong vol/year → "The *quantum* GH
  propinquity," Trans. AMS **368 (2016)** 365–411.
- `avery_wen_avery2002` p44 wrong title/initial → "Some properties of
  hyperspherical harmonics," Z.~Y.~Wen.
- `devastato_lizzi_martinetti2018` p44 wrong-ID → "Lorentz signature and
  twisted spectral triples," Devastato–**Farnsworth**–Lizzi–Martinetti, JHEP
  03 (2018) 089, arXiv:1710.04965.
- Connes-vS "Definition 2.39" → honest section ref (§2.7); the fetch found
  ~Def 2.29 but couldn't fully confirm, so referenced the named section.
- **Theorem-number verification (PI-requested):** van den Dungen Prop 4.1
  (iᵗD̸ Krein triple) + Nieuviarts Def 2.2 (twist-by-grading, even-dim — the
  odd-S³ NO-GO) both **GROUNDED** vs primary PDFs (no fix).
- p43 dangling `\cite{paper44}` → added the missing bibitem.

**Status / C7 labels:**
- Paper 38 "propinquity" → "state-space Gromov–Hausdorff" (p42 ×6; "strengthens"
  → "complements" at l.250); P39 zombie title "Tensor-product propinquity" →
  synced to its real state-space-GH title (p44).
- synthesis Paper-45 subsection title "…propinquity *convergence*" + "closes the
  convergence theorem" → degeneracy-theorem framing (Paper 45 is a NEGATIVE
  result).
- paper_42 intro/§10 reconcile (the genuine residual from the over-flagged §10).

All 4 files compile errors=0 / undefined=0; C11/C13/C14 gates PASS.

## 4. Honest scope
- **Sound (no change):** paper_42 derived-finding (verified correct); §10
  finite-cutoff Krein closure (backed). paper_43/44 finite-cutoff results
  (BACKED-SOUND per code reviewers; 311 + 177 tests pass).
- **Deferred:** UNVERIFIABLE-but-low-stakes theorem-numbers fully checked where
  load-bearing (van den Dungen, Nieuviarts grounded); the remaining advisory
  debug/ refs in 42/44 are part of the corpus-wide 443-ref sweep.

## 5. Next
Bite B sub-bite 2 (Papers 45–49, the descoped/partial core) then sub-bite 3
(39, 52, 53).
