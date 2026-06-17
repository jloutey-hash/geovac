# Paper 40 propinquity determination memo

**Date:** 2026-06-16
**Question:** Does Paper 40 (`paper_40_unified_propinquity_convergence.tex`) GENUINELY establish
Latrémolière propinquity convergence for all compact Lie groups, or only the weaker
van Suijlekom state-space Gromov–Hausdorff convergence — i.e., is its title + abstract
the same overclaim Paper 38 was corrected away from?

**Verdict: STATE-SPACE-GH-ONLY.** Paper 40's title, abstract, main theorem, and corollaries
claim full **Latrémolière propinquity** convergence, but the paper's own L5 proof concedes that
the load-bearing **dual-reach (`reach_P`) step is a NAMED GAP** that is "corroborated at the
state-space Gromov–Hausdorff level" — exactly the gap for which Paper 38 was descoped from
"Latrémolière propinquity" to "van Suijlekom state-space GH." The reason transfers in full.
Recommend descoping Paper 40's title / abstract / main theorem / corollaries to state-space GH,
**and** the general-G case is strictly weaker than even Paper 38's SU(2) repair (see §4).

---

## 1. What Paper 40 claims vs. what its proof delivers

### What it CLAIMS (propinquity, unconditional)
- **Title:** "Latrémolière propinquity convergence of spectral truncations of compact Lie groups…"
- **Abstract (L93–119):** "We prove that 𝒯_Λ converges to 𝒯_G in the **Latrémolière quantum
  Gromov–Hausdorff propinquity** for metric spectral triples… The result **upgrades** the
  state-space Gromov–Hausdorff convergence proved by Gaudillot-Estrada and van Suijlekom… **to
  Latrémolière propinquity** with explicit Dirac-controlled rate." (No conditional qualifier in
  the abstract.)
- **Main Theorem `thm:main` (L1708–1730):** "Then 𝒯_Λ converges to 𝒯_G in the **Latrémolière
  quantum Gromov–Hausdorff propinquity**: Λ_prop(𝒯_Λ, 𝒯_G) ≤ C₃(G)·γ_Λ(G)."
- **Corollaries (L1768+):** SU(N), Spin(n), Sp(n), exceptional groups "converge … in
  **Latrémolière propinquity** at the universal rate."

### What the PROOF delivers (state-space GH, with a hole at general G)
The propinquity number is built in **Lemma L5** (`lem:L5`, L1644–1702) from four tunnel
constituents: forward reach, **dual reach (`reach_P`)**, forward height, dual height. Three are
proved. The fourth is conceded:

> **L1666–1673:** "The dual-reach estimate reach_P ≤ γ_Λ(G) is a **named gap**
> (Remark `rem:correction_L2`): the smoothing map has cb-norm 1, and **a forward cb-norm controls
> no partial inverse.** A corrected proof requires an antiderivative-transference lemma in the
> style of Leimbach–van Suijlekom; for the scalar sector **the conclusion is corroborated at the
> state-space Gromov–Hausdorff level** for all compact metric groups by [Gaudillot-Estrada–vS]."

> **L1699–1701:** "The combination is **stated modulo the reach_P named gap** of
> Remark `rem:correction_L2`; all other constituents are proved above."

And `rem:correction_L2` (L711–746) is even more explicit about the framework actually in force:

> **L729–731:** "(ii) the tunneling-pair vocabulary of §L5 follows Paper 38; **the device actually
> used is van Suijlekom's state-space Gromov–Hausdorff framework.**"

So the object Paper 40 actually controls is the **state-space GH distance** (the same object
Paper 38 ends up proving), and the **propinquity number is bounded only modulo an unproved
dual-reach lemma.** The title/abstract/theorem assert propinquity; the proof delivers state-space
GH plus a named gap. This is a prose > delivered-content overclaim of exactly the Paper 38 type.

---

## 2. Why Paper 38 was descoped — the SPECIFIC reason

Source chain: CHANGELOG v3.106.0 (L808–818, L899–926), `debug/sprint_p45_hardening_phase1_memo.md`,
`debug/p45_descope_siblings_log.md`, Paper 38's in-paper `rem:history38` (L532–543), CLAUDE.md
§6 live-status flags, memory `wh1_proven`.

The 2026-06-09 adversarial-hardening sprint (run on Paper 45, cascaded to the whole propinquity
cluster) found **three compounding failures**, two of which are the Paper-38 descope reason:

1. **`reach_P` dual-reach gap (the load-bearing one).** The original propinquity proofs bounded the
   dual reach via the kernel's **cb-norm**. That is invalid: a forward cb-norm = 1 controls no
   partial inverse (the inverse is controlled by the symbol *minimum*, O(n²) on the full envelope).
   So the propinquity number was never actually bounded — only the state-space GH object was.

2. **Fabricated Latrémolière citations.** Paper 45's "closure-by-citation" sprint (2026-05-24) had
   cited "**Latrémolière Thm 5.5 / Def 3.4 / Def 3.5**" as the mechanism that closed the L5
   propinquity bookkeeping. Both candidate source PDFs were fetched: **arXiv:1811.10843 has no
   Theorem 5.5 (its §5 ends at Remark 5.3); items 3.4/3.5 are unrelated.** The framework actually
   used is **van Suijlekom's state-space GH** (arXiv:2005.08544; Leimbach–vS 2024;
   Gaudillot-Estrada–vS arXiv:2310.14733), dressed in Latrémolière propinquity vocabulary. The
   "Latrémolière propinquity" label across the arc was the overclaim.

3. **L2 cb-norm / symbol conflation** (the Plancherel *mass distribution* (2j+1)/Z was confused with
   the convolution *multiplier symbol*). This is what made the invalid cb-norm `reach_P` argument
   look plausible.

**Resolution (2026-06-10).** Paper 38's SU(2) theorem was rewritten and made **unconditional** —
but NOT as Latrémolière propinquity. It was reframed to **van Suijlekom's state-space GH distance**,
metrized by the **left-translation Lipschitz seminorm** (Def `def:translation_seminorm`,
`thm:main_unconditional`, L449–459). The `reach_P` gap was **dissolved, not bridged**: both
almost-inverse defects collapse to Fejér smoothings at the same moment γ_n, so **no inverse /
transference estimate is used at all** (`rem:no_inverse`, L498–516). Paper 38's title is now
"**State-space** Gromov–Hausdorff convergence …"; its `rem:history38` records the May-2026 draft's
"Latrémolière propinquity" statement as an overclaim resting on the L2 conflation + invalid
dual-reach argument.

Net: **the Paper 38 descope reason = "the Latrémolière-propinquity claim was bounded only via an
invalid cb-norm dual-reach step (backed by nonexistent Latrémolière theorem numbers); the object
genuinely proved is van Suijlekom state-space GH."**

---

## 3. Does the same reason apply to Paper 40? — YES, by the paper's own text

| Descope criterion (Paper 38) | Paper 40 status |
|---|---|
| Invalid cb-norm `reach_P` dual-reach argument | **PRESENT and CONCEDED.** L1666–1673 says reach_P is a "named gap"; "a forward cb-norm controls no partial inverse." Identical defect, openly labeled. |
| Object actually controlled = state-space GH, not propinquity | **CONCEDED.** L729–731: "the device actually used is van Suijlekom's state-space GH framework." Corroboration cited is GE–vS state-space GH (L1671–1673). |
| Title/abstract/theorem still say "Latrémolière propinquity" | **YES** — title, abstract (L93, L116–119), `thm:main_intro`, `thm:main`, every corollary. The descope was applied to L2 and to the L5 *proof body* (named-gap note), but **NOT to the title/abstract/theorem statements.** |
| Fabricated Latrémolière Thm/Def numbers | **NOT present in Paper 40.** Paper 40 cites the propinquity framework only generically (`latremoliere_metric_st_2017`, which is real — arXiv:1811.10843 / Adv. Math. 404 (2022), web-confirmed) and `latremoliere2016` / `latremoliere2018` (the standard lineage). It does **not** cite the nonexistent Thm 5.5 / Def 3.4 / 3.5. So Paper 40 is cleaner than Paper 45 on the citation axis — but the **substantive** gap (the invalid reach_P) is the same, and the title overclaim is the same. |

The 2026-06-09 sibling-descope (`debug/p45_descope_siblings_log.md` §"Paper 40") applied an
**L2-correction-only** treatment to Paper 40 (no full Lorentzian template). It fixed L2(b)/(c), the
Berezin definition, and the L5 dual-reach proof step (→ named gap), and added `rem:erratum_L2`. It
**did not touch the title, abstract, or main-theorem statements** — leaving Paper 40 in exactly the
inconsistent state Paper 38 was in *before* its 2026-06-10 in-place rewrite: corrected proof body,
but a title/abstract that still assert the stronger result the proof no longer delivers.

---

## 4. Extra finding: general-G is strictly WEAKER than Paper 38's SU(2) repair

Paper 38's unconditional state-space-GH result depends on a **group-specific construction**: the
exact spinor-window inclusion (`lem:lifted_state`(a), Paper 38 L393–427) where the SU(2)_L×SU(2)_R
decomposition V_j ⊗ V_{j±1/2} matches the Camporesi–Higuchi shells exactly, so the lifted state
sits in complete shells with no spillover, and both defects reduce to one Fejér moment.

Paper 40's own `rem:correction_L2` 2026-06-10 update (L733–745) says this repair route does **not
yet** exist for general G:

> "The same translation-seminorm route is the natural repair path for general G (the lifted-state
> construction requires the spin-module window decomposition of the relevant Dirac shell structure
> per group); **the general-G statement remains a named gap until that decomposition is checked
> case-by-case or uniformly.**"

So for general G, Paper 40 has **neither** a proven Latrémolière propinquity result **nor** a proven
unconditional state-space GH result — it has the state-space GH framework with a per-group
lifted-state construction that is explicitly unverified beyond SU(2). The honest current scope of
Paper 40's general-G theorem is: *conditional* state-space GH convergence (modulo the per-group
spinor-window decomposition), with the rate γ_Λ(G) and the universal 4/π constant being the genuine
new content (those are kernel moments, independent of the gap).

---

## 5. Recommendation (no edits made, per task)

**Descope Paper 40 to state-space GH, mirroring the Paper 38 repair**, and additionally flag the
general-G conditionality:

1. **Title:** "Latrémolière propinquity convergence …" → "State-space Gromov–Hausdorff convergence
   …" (parallel to Paper 38's retitle).
2. **Abstract:** drop "converges … in the Latrémolière quantum Gromov–Hausdorff propinquity" →
   "converges … in van Suijlekom's state-space Gromov–Hausdorff distance"; drop "**upgrades** …
   **to Latrémolière propinquity**" → "**supplies an explicit Dirac-controlled rate for** the
   state-space GH convergence." State the general-G lifted-state decomposition as a named gap.
3. **`thm:main_intro` / `thm:main` / all corollaries:** "Latrémolière propinquity" → "state-space
   GH distance"; add the per-group lifted-state caveat (or restrict the unconditional statement to
   SU(2)/rank-1 and state the rest as conditional).
4. **Keep unchanged** (genuine content): the closed-form rate moment γ_Λ(G), the universal 4/π
   theorem (`thm:universal_constant_analytical` + numerics), the L3 Lipschitz constant C₃(G) ≤ 1,
   the L1′/L4 substrate. These are not affected by the propinquity-vs-state-space distinction.
5. Consider whether Paper 40's general-G claim should be downgraded harder than Paper 38's SU(2)
   claim, since §4 shows the unconditional translation-seminorm repair is **not yet available** for
   general G. Paper 38 (SU(2)) is unconditional; Paper 40 (general G) is conditional even at the
   state-space-GH level. This is the opposite of "Paper 38 may be inconsistently under-scoped" — if
   anything Paper 40 is currently *over*-scoped relative to Paper 38.

**Process note:** This is a clean instance of the descope having been applied to a paper's *proof
body* but not propagated to its *title/abstract/theorem statements* — the same lag the branch-QA
sweep is catching elsewhere (cf. the 24 stale "conjectural" K-labels). The Paper 40 title/abstract
were left asserting "Latrémolière propinquity" after the 2026-06-09 sibling pass corrected only L2
and the L5 proof step. The fix is mechanical and parallels the already-completed Paper 38 retitle.

---

## Evidence index
- Paper 40: title L71–73; abstract L93–119; `thm:main_intro` L203–221; `rem:correction_L2`
  L711–746 (esp. L729–731 device, L733–745 general-G gap); `lem:L5` proof L1657–1702 (esp. L1666–1673
  named gap, L1699–1701 "modulo the reach_P named gap"); `thm:main` L1708–1730; corollaries L1768+;
  bibitems L2238–2253 (Latrémolière refs all real, generic).
- Paper 38: retitle L62–63; `thm:main_unconditional` L449–459; translation seminorm L313–332;
  `rem:no_inverse` L498–516 (gap dissolved); `rem:history38` L532–543 (overclaim record).
- CHANGELOG v3.106.0 L808–818 (P38 unconditional), L899–926 (descope sprint, fabricated-citation
  finding L902).
- `debug/sprint_p45_hardening_phase1_memo.md` §1–§4 (three failures; §4 "Paper 40 L2(b),(c) repeats
  it for all compact groups"; §5a execution: "Papers 40/46/47/48/49 … erratum blocks/remarks
  applied per template").
- `debug/p45_descope_siblings_log.md` §"Paper 40" (L2-correction-only; title/abstract/theorem NOT
  touched).
- Web: arXiv:1811.10843 = Latrémolière, "The Gromov–Hausdorff propinquity for metric Spectral
  Triples," Adv. Math. 404 (2022), Paper 108393 — EXISTS (the framework reference is sound; only the
  specific Thm 5.5/Def 3.4/3.5 numbers cited in Paper 45, not in Paper 40, were fabricated).
