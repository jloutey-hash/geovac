# Sprint: /qa group4 5th cert (FAIL) + remediation — canonical memo

**Date:** 2026-06-29 · **Version:** v4.58.0 · **PI direction:** `/qa group4` (5th cert run);
then **"you can start remediation."** Predecessors: v4.52–v4.57 (pre-work → 4 prior cert/remediation cycles).
Answer key (gitignored): `debug/qa/group4_seed_key.json`.

## 1. The 5th cert run (@ 6e4d137/v4.57.0) — the closest to PASS yet

Whole-group, frozen DoD, 8 fresh seeds (new loci, incl. seeds on the two new v4.57.0 tests), fresh
path-pinned panel (4 code 1/paper, 2 claims chunks, 1 citation, 1 synthesis; opus). **Verdict: FAIL —
fully calibrated.** Sensitivity **8/8** (every seed caught), specificity **6/6** — every accumulated
fix (v4.55–4.57) was *confirmed accepted*, not re-flagged: **the Trenev reversal triple-confirmed**
(two code reviewers + the citation reviewer reading the PDF Appendix B), plus F2, S1, λ table, Chawla,
eri disclosure. **6 of 8 dimensions came back clean-except-their-seed.**

## 2. Genuine findings (all remediated unless flagged)

| # | Finding | Fix |
|---|---|---|
| **#2** (claims/P20, LARGE) | The Paper-38 GH-**metric** rate γ oversold as a chemistry **energy** bound — §VI "basis-truncation error bound" + "first quantitative … in any framework" + the metadata message "this energy sits within γ", contradicted by the paper's own footnote (energy lift overshoots 3–4 OoM). | §VI reworded: γ is a **proven state-space-GH convergence rate, NOT a direct energy-error bound** (heading "…rate"; "convergence-rate estimate"; metadata message → "operator within GH distance γ of the n_max→∞ limit"). |
| **M-A** (code/P20, LARGE) | The raw STO-3G **907** (the 2.7×/13× market-test denominator) was swept under "the Gaussian counts are Trenev Table 5 (2-qubit reduction)" — but Table 5's reduced STO-3G LiH is **276**, and 907 is the raw-JW value. **Caused by my own v4.57.0 Trenev reversal over-blanketing.** | Caption carves it out: the **2-qubit-reduced** counts (LiH 6-31G/cc-pVDZ + all H₂O) are Trenev Table 5 (App. B); the **raw 907** STO-3G LiH is the raw-JW matched-raw baseline, not a Trenev-reduced count. |
| **M-C** (code/P20, MATERIAL) | "Z=1–36" contradicted by **SrH (Z=38) + BaH (Z=56)** in the 37-registry (probe-confirmed). | "Z=1–36" → **"Z=1–56 (H through Ba)"** at the library-description loci (P20 ×2, synthesis ×1). The scaling-context "three rows" (11.10·Q main-group universality) is accurate and kept. |
| **M-B** (code/P20) | 0.20% @ n_max=3 abstract headline has no CI test (84-qubit FCI). | Logged as a deferred coverage gap (value correct, body-computed; CI test heavy). |
| **claims SMALL** | #3 P14:2782 "comparable accuracy", #4 P20:475 "comparable basis quality" (§1.5); #5 P14:1580 CF-1 lead missing pair-diagonal qualifier; #6 P14:980 "composed 7.0–8.8% [Paper 17]" mislabel (composed=5.3%; 7–8.8% is **balanced** [Paper 19]). | #4 fixed (me); #3/#5/#6 fixed (agent). |
| **NITs** | P20 §6.3 NaH lstlisting shows 239 (balanced) for composed `hamiltonian('NaH')`; atomic_classifier docstring 'A'–'E' missing 'F'/'d_block'; Pachucki prose "2023"→2018 (now resolvable — bibitem 2018 verified to support the claim); P14 "recomputed Gaussian"→"published Trenev counts". | NaH→composed value (me); docstring/Pachucki/wording (agent). |

## 3. Deferred / flagged (non-blocking carve-out or PI-level)

- **M-A market-direction decision (FLAG to PI):** the raw-907 is the matched-raw market-test baseline (2.7×); under Trenev's reduced 276 it's parity-to-worse. The caption brackets both (raw-vs-raw and reduced-vs-reduced). Whether the *headline* should lead raw-vs-raw or reduced is a PI call (it sets the market-test direction). The raw 907's own provenance/backing is also open.
- **"three rows" framing** — pervasive (~7 loci across P20/synthesis/P14) and now loose vs Z=1–56; the scaling-context ones are accurate. A consistency sweep, NIT.
- secondary-number reconciliations (H₂O 1-norm 175.6 vs 361; LiH 1-norm 33.3 vs 32.59; d-only 55 vs 56; BeH₂ 11.7% vs 32%) — need a which-is-right check; NIT carve-out.
- citation NITs: burkat title-paraphrase, duplicate caesura bibitem, cosmetic key-misnomers.
- synthesis market-test raw-vs-reduced completeness + dangling Paper 38/24 in-prose refs.

## 4. Honest scope

- **Theorem grade:** none. **Process verdict:** fully calibrated (8/8, 6/6) → trustworthy.
- **The branch is one thin layer from clean:** no headline NUMBER was wrong except the scope Z-range
  (M-C); #2 and M-A are framing/attribution-precision; the rest SMALL/NIT. The 4 prior rounds of
  remediation were all confirmed accepted (the reversal triple-verified).
- **M-A is mine to own** — the v4.57.0 reversal, correct for the reduced counts, over-blanketed the
  raw-907 baseline; carved out now.
- **Re-cert HELD** per the gate (PI timing).
