# Sprint: /qa group4 4th cert (PASS-confirmation, FAIL) + remediation — canonical memo

**Date:** 2026-06-29 · **Version:** v4.57.0 · **PI direction:** released the v4.56.0 recert hold
(`/qa group4`, the certified-PASS confirmation run); then **"let's start remediation"**; the rel-λ
conflict precedent → I ran a **"diagnose first"** check on the Trenev provenance, which produced
the headline outcome below. Predecessors: v4.52–v4.56 (pre-work → 3 prior cert/remediation cycles).
Answer key (gitignored): `debug/qa/group4_seed_key.json`.

## 1. The 4th cert run (@ 4c17182/v4.56.0)

Whole-group, frozen DoD, 8 fresh seeds, fresh path-pinned panel (4 code 1/paper, 2 claims chunks,
1 citation, 1 synthesis; opus) + a Trenev-Table-5 web verification. **Verdict: FAIL — fully
calibrated.** Sensitivity **8/8** (all seeds caught), specificity **6/6** — the panel *confirmed*
every v4.55/v4.56 fix (the corrected λ table 40.59/143.96/218.78/413.90, R_eq=3.227, the Chawla
cite) rather than re-flagging it. Worktree torn down, leak-clean.

A thin converging layer: the headline numbers are all now correct + test-pinned; what FAILed was
**incomplete propagation of prior fixes into the synthesis + two zombies + a NIT cluster** — plus
the Trenev reversal below.

## 2. THE HEADLINE: the Trenev attribution was reversed — TWO prior cycles had it backwards

The cert flagged a code-vs-paper contradiction (Z1): `composed_qubit.py` hardcodes the Gaussian
baselines (LiH 276/5851/63519, H2O 551/8921/107382) with a **"Source: Trenev Table 5"** comment,
while the v4.55.0 papers claimed **"GeoVac's own OpenFermion recompute; Trenev is vibrational with
no electronic counts; cited for methodology only."** Per "diagnose first" I web-verified Trenev
(arXiv:2311.03719) **directly**:

> Trenev's **Appendix B / Table 5 ("Electronic structure vs Vibrational structure")** DOES tabulate
> electronic-structure Gaussian-basis JW Pauli counts for LiH/H2O — and **all six GeoVac values match
> Table 5 EXACTLY**, including qubit counts (190× = 63519/334 reconciles).

So the **code was right all along**, and the **v4.54.0 + v4.55.0 trenev remediation was a double
mistake** — it concluded "Trenev vibrational-only, no electronic counts → re-attribute to GeoVac's
recompute" because that earlier citation check saw the vibrational title/body but **missed
Appendix B**. The current papers carried a *false* "GeoVac's own recompute" provenance. **Reverted
corpus-wide** to the verified truth: *the Gaussian baseline counts are Trenev Table 5 (App. B); the
Q^3.9–4.3 scaling exponents are GeoVac's own log-log fit of those published counts.* (P14 11 loci +
bibitem, P20 caption + refs.bib, synthesis body + bibitem; code docstrings restored "vibrational" +
"Appendix B" — kept the correct "Trenev Table 5" source.) **No numeric value changed** (counts,
exponents, multipliers all stand). Grep-verified: 0 false-framing loci remain.

**Lesson:** "diagnose first" (the rel-λ precedent) prevented corrupting correct code to match a
wrong paper — and caught that a citation-reviewer's "vibrational title ⇒ no electronic counts"
inference had propagated through two remediation cycles. Read the cited *table*, not just the title.

## 3. The other findings (all remediated)

| # | Finding | Fix |
|---|---|---|
| **F2** | P20:957 "composed**/balanced** … does not bind under qubit FCI" — the "/balanced" is the retracted R3-A/B claim (balanced LiH *binds* at R_eq=3.227, the flagship result) | dropped "/balanced" → "LiH composed and NaH balanced" |
| **S1** | synthesis L238 "binds LiH at experimental R_eq=**3.015**" — the v4.56.0 M2 fix didn't propagate to the synthesis | reworded (computed R_eq=3.227, 7% high; 0.20% = n=3 single-point at the experimental geometry) |
| **F3/F4/F5** | P14:2967 "matched-accuracy" → matched-qubit; P20 17–32× bullet missing the RaH (Z=88) species-mismatch caveat (F4 "STO-3G-equivalent accuracy" phrase not found — covered by the P14:56 matched-qubit caveat) | matched-qubit wording + species caveat added |
| **NIT** | ScH "278 … Pauli/Q=9.23" convention mix; rocca 4th author; synthesis "integral values" | → 277 non-identity (333 for LiH); P.D.~Johnson→**P.J.~Ollitrault**; "integral values"→"published Pauli counts (Table 5, App. B)" |

**New tests (close recurring coverage flags):** `tests/test_paper20_library.py` (library = 37
distinct systems — was unasserted) + `tests/test_paper16_dirac_metric.py` (the §VI metric-vs-
topological claim: μ_free(137)=36720 + δ_1s table values + smoothness — was a multi-cert recurring
no-test flag; adjudicated NIT = correct *analytic* claim, now formula-verified). 6/6 green.

## 4. Verification

- New tests 6/6 green; deterministic gates C11/C13/C14/C15/C16 `--gate group4` PASS.
- All edited papers compile clean: P14 26 / P20 12 / synthesis 4 pp; 0 undefined. composed_qubit
  imports (docstring-only edit, no logic; 37 systems).
- Trenev reversal grep-verified: 0 "own recompute / methodology only / no electronic counts" loci;
  correct "Table 5 (Appendix B)" framing present.

## 5. Deferred NITs (noted, non-blocking per the secondary-number/provenance carve-out)

Pachucki prose "2023" vs bibitem-2018 + the FW-two-particle attribution (genuine source-check
needed); Navrátil bibitem title↔venue conflation; duplicate caesura bibitem (same arXiv:2501.06165,
two keys); 3 dangling P20 `.bib` entries (McArdle2020/Bauer2020/Cao2019); "M = n_max²" derivation
text (should be Σn²; tables correct — needs a careful exponent rework, empirical M^{-0.85} stands);
BeH2 32%-vs-11.7% (two methods, clarity); P14:981 "7–8.8%" mislabeled "composed" (balanced); H2O
exponent 2.22-vs-2.52 (2-pt vs 3-pt fit); rel λ/Pauli **n=3** + the 0.20% n=3 accuracy (heavy
84-qubit builds — coverage deferred; n=2 is pinned).

## 6. Honest scope

- **Theorem grade:** none. **Process verdict:** the 4th cert was fully calibrated (8/8, 6/6) →
  trustworthy. Its value was catching the **Trenev reversal** (a real provenance error that survived
  TWO prior remediation cycles — the fresh-adversary + "diagnose-first" discipline working as
  designed, a fourth time) and the synthesis propagation gaps.
- **No numbers changed** in the Trenev reversal — it is a provenance correction (counts ARE Trenev's
  published Table 5; exponents ARE GeoVac's fit). The market-test multipliers (2.7×/190×) stand on
  Trenev's verified electronic counts.
- **Re-cert HELD** per the gate (PI timing) — the certified-PASS confirmation is the next `/qa group4`.
