# Group 1 (Operator algebras / NCG) — `/qa` profile

> **Inherits the shared criteria in [`docs/qa/criteria.md`](criteria.md).** This
> file supplies only group1-specific scope + deltas + the branch-specific
> criterion C14.

> **STATUS: DRAFTED 2026-06-16 — NOT YET FROZEN.** Sequenced *after* group3
> bite 2. Awaiting PI confirmation of scope + first bite before FREEZE.

**Scope (non-trunk group1):** Papers **29, 39, 40, 42, 43, 44, 45, 46, 47, 48,
49, 50, 52, 53** + the **group1 operator-algebras synthesis**. Trunk papers
**32, 38** taken as already-certified (`/qa trunk` PASS); in scope only where a
group1 paper restates them (C7).

**Deterministic `--gate`:** `group1`.

## Branch deltas (the only non-inherited content)

- **C14 — Descope / partial status accuracy (branch-specific; the defining
  criterion).** Every DESCOPED or PARTIAL paper presents withdrawn claims at
  current status:
  - Papers **45, 46 = DESCOPED** — the "first Lorentzian propinquity
    (convergence) theorem" is **retracted**; Paper 45's headline is the K⁺
    annihilation theorem; no body text presents weak/strong-form Lorentzian
    propinquity convergence as established.
  - Papers **47, 48, 49 = PARTIAL** — surviving arrows stated as surviving (47
    norm-resolvent + three-carrier; 48 bridge as conditional design; 49
    cocycle-deficit / OSLPLS algebra); descoped metric-level / Λ-inheritance
    claims flagged descoped.
  - The **product-carrier S³×S¹ convergence** is presented as
    **signature-agnostic (Euclidean), NOT a Lorentzian claim**.
  - Each carries an in-place Status note; no pre-descope claim cited as live.
  *(If a second branch later needs C14, promote it into `criteria.md`.)*
- **C4 high-priority for this branch.** The math.OA citation apparatus
  (Connes–vS, Latrémolière, van Suijlekom state-space GH, Marcolli,
  Camporesi–Higuchi, Bizi–Brouder–Besnard, Mondino–Sämann, Datta) is the
  corpus's highest fabrication surface and *already produced a real defect here*
  (the nonexistent "Latrémolière Thm 5.5 / Def 3.4" refs that drove the Paper 45
  descope; the Paper 38 `avery_wen_avery` v4.14.1 fix). Every cited theorem/def
  *number* is verified.
- **C7 (trunk-dependent status).** Paper 38/WH1 as PROVEN scoped to the van
  Suijlekom state-space GH distance; Paper 32 axioms / κ at current tiers.
- **C8 (headline honesty), per-paper.** Paper 29 = Ramanujan + integer-algebraic
  Ihara zeros; Paper 40 = 4/π rate **universal** (derived-numerics-pinned + PRV
  tightness, not full symbolic proof beyond the verified panel); Paper 50 =
  **bit-exact** CFT-on-sphere F-theorem match (arithmetic fact, not a CFT
  derivation); Paper 45 headline = the **annihilation theorem** (a NEGATIVE
  result), never a convergence theorem.

## Proposed first bite (PI to confirm at FREEZE)

- **(recommended) Papers 29, 40, 50 + synthesis** — the ACTIVE math.OA backbone
  (Ramanujan / universal 4/π rate / CFT F-theorem); calibrate the panel on the
  non-descoped core, C9 scoped to the synthesis's coverage of those three plus
  its Lorentzian-descope narrative (coherence-passed v4.19.0).
- **Bite 2 (proposed):** the **Lorentzian cluster 42–49** — the C14-heavy,
  highest-risk subset, after the panel is calibrated.
- *Alternative: attack the Lorentzian cluster first if you'd rather hit the risk
  head-on while the descope is fresh.*

## Change log
- 2026-06-16 — **DRAFTED** by PM for PI review (third pre-registered `/qa`
  target). Adds C14 (descope/partial status); flags C4 high-priority. Scope
  excludes trunk-certified 32/38. **Awaiting PI freeze** (sequenced after group3
  bite 2).
- 2026-06-17 — **Paper 29 (single-paper bite) = FAIL → FIXED** (v4.21.3).
  PI-scoped to one paper (token slow-roll); dimensions deterministic + code +
  claims + citation (synthesis excluded, single paper). Calibrated panel
  (sensitivity 4/4, specificity 4/4). Science sound; fixed a degree-arithmetic
  error (84→80), two C4 citation conflations (matsuura → JHEP 09(2022)178,
  yakaboylu → H-P-Hamiltonian arXiv:2309.00405), the McKenzie initial, + closed
  the S⁵-N3 closed-form coverage gap (new `test_s5_N3_zeta_matches_paper_closed_form`).
  Seed key `debug/qa/group1_seed_key.json`.
- 2026-06-17 — **Bite A (Papers 40, 50 + synthesis) = FAIL → REMEDIATED**
  (v4.22.0). Calibrated panel (sensitivity 4/4, specificity 5/5). Genuine
  defects fixed: 2 wrong-ID + 1 fabricated citation (beccaria→1406.3542,
  hartman→1807.11401, hekkelman-T^d removed→Leimbach-vS); paper_50 F-theorem
  false-positive backing → genuine `test_paper50_f_theorem.py` (ζ′(0) from
  framework spectrum → KPS 1e-82); 3 synthesis/paper-50 Lorentzian/propinquity
  zombies; predictive-CFT + over-rigor calibration. **Paper 40 rank≥2
  universality**: pruned `l2_universal_rate_memo` resurrected from git history
  → genuine computation but **fit-sensitive** extraction; calibrated to rank-1
  rigorous + rank≥2 A-over-B robust (named gap), backing
  `test_paper40_universal_rate.py`. Seed key `group1_A_seed_key.json`, memo
  `debug/sprint_qa_group1_biteA_memo.md`. **Remaining: Bite B (39, 42–49, 52,
  53) — Lorentzian cluster, smaller sub-bites.**
- 2026-06-18 — **Bite B sub-bite 1 (Papers 42, 43, 44 + synthesis) = FAIL →
  REMEDIATED** (v4.23.0). 10-agent panel calibrated (sensitivity 5/5,
  specificity 5/5). **Two reviewer LARGE findings OVERTURNED by PM verification**
  (§9 reconcile): paper_42's "derived finding is false" (reviewer omitted the
  β factor — D_W generator is 2π·D_W, does NOT close; finding CORRECT, no
  reframe) + §10 "descoped-zombie" (over-flag — genuine finite-cutoff Krein
  closure; fixed the stale intro disclaimer instead). Genuine fixes: 6 wrong-ID/
  fabricated citations (hekkelman_mcdonald2024 fabricated ×3, hekkelman2022,
  zhu_casini2020 authors, latremoliere2018 ×3, avery, devastato) + theorem-
  number verification (van den Dungen Prop 4.1 + Nieuviarts Def 2.2 GROUNDED) +
  Paper-38 propinquity→state-space-GH labels + P39 zombie title + synthesis
  Paper-45 degeneracy framing + paper44 bibitem. All compile clean; C11/C13/C14
  PASS. Seed key `group1_B1_seed_key.json`, memo
  `debug/sprint_qa_group1_biteB1_memo.md`. **Remaining: sub-bite 2 (45–49),
  sub-bite 3 (39, 52, 53).**
