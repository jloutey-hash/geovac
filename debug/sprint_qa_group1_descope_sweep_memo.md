# Sprint memo — group1 Lorentzian-cluster inline descope-tagging sweep

**Date:** 2026-06-22
**Version:** v4.37.0
**Trigger:** PI direction (AskUserQuestion: "Descope-tagging sweep, then cert4") after cert3 (v4.36.0) identified the systemic issue.

## 1. Why

The cert re-runs (cert1=4 defects, cert2=1, cert3=3) were not converging because the Lorentzian cluster (Papers 45–49) carries many theorem statements whose descope lives only in a distant Status note, not inline on the statement. Each fresh adversarial read surfaced another un-tagged metric-level-convergence claim:
- cert1 (v4.35.0): §3 axiom lift (p48).
- cert3 (v4.36.0): Bridge-Theorem B4 convergence-transport + T2 (p48).

Same defect class. Blind cert4 would likely find the next one. This dedicated sweep drains the tail in one pass (analogous to the folded-in C4 citation sweeps that drained the citation tail).

## 2. Method

Enumerated EVERY `\begin{theorem}/proposition/corollary` environment in P45/P46/P47/P48/P49 (done directly via grep+read — the 3 enumeration agents hit a transient server classifier outage; read-only ops don't need it). For each: classified METRIC-CONV-DESCOPED vs SURVIVING vs OTHER, and checked whether the statement itself (title or adjacent clause, NOT a distant Status note) carries an inline descope tag.

## 3. Gaps tagged (6)

| Paper | Statement | What it asserts | Fix |
|---|---|---|---|
| p47 | `thm:inner` | inner-arrow propinquity convergence (Λ_prop ≤ C₃·γ→0) on the degenerate Λ_prop | title → "Inner arrow — descoped" + inline caveat |
| p48 | `thm:bridge` (B4 item) | convergence transport (BigDeth→0 ⇒ pLGH conv.) | "(descoped)" + degenerate-Krein-metric caveat (abstract summary was tagged in cert3, not the theorem) |
| p48 | `thm:relaxed_triangle` | 4-pt relaxed triangle of the Krein metametric | inline structural/categorical note (metametric built on degenerate seminorm) |
| p48 | `thm:coincidence` | d=0 ⇔ isometry | inline note: the "only if" leg is exactly what degeneracy breaks → conditional on repair |
| p48 | `thm:inframetric` | BigDeth is an inframetric (non-degeneracy) | inline structural/categorical note |
| p46 | `prop:reach_height_op` | the four Latrémolière propinquity constituents | inline note: rate-formula quantities feeding the descoped thm:main, not a metric convergence |

## 4. Deliberately NOT tagged (surviving — tagging would falsely descope)

p47 thm:outer / thm:three_carriers (norm-resolvent, proof-by-argument); p47 thm:main (norm-resolvent convergence; R_inner descoped in the adjacent Remark); p48 thm:compact_agreement / N_t=1 reduction / prop=2 / B1 / B3 / T1 (wedge-is-a-pLPS, structural) / T4 (four-witness thermal) / T5 (Pythagorean foliations); p49 B1'/B2'/B3' (B3' cardinality bound is combinatorial, not metric-dependent — the cert3 "NIT" is genuinely surviving, so left untagged) / TICI / D_max chain; p45 K⁺ annihilation (negative) / product-carrier (signature-agnostic). P45/P49 descoped items (P45 negative; P49 B4' DESCOPED) were already inline-tagged.

## 5. Verification

p46/p47/p48 compile errors=0/undef=0; C5/C11/C13/C15 group1 PASS. p45/p49 unchanged (already covered).

## 6. Honest scope

This is a corrective tagging sweep, not a /qa run (no seeds/calibration). It does NOT certify Batch 1 — it removes the systemic gap so the next calibrated cert re-run (cert4) can converge to a clean PASS. The claim is: every theorem/proposition in the cluster that asserts a metric-level Lorentzian/propinquity *convergence* (or a degenerate-metric property) now carries an inline descope tag, so a fresh adversarial reader checking statement-by-statement should find no new "Status-note-only descope" zombie. If cert4 still finds one, it is a genuinely-missed statement (the enumeration was line-by-line, so this should be exhaustive).

---

## Round 2 — prose/abstract descope sweep (2026-06-22, v4.39.0)

**Trigger:** cert4 (v4.38.0) found that the v4.37.0 theorem-statement sweep held perfectly, but descope zombies remained in **§1-Introduction prose** (p47) and **abstract enumerations** (p48 B4). PI direction: run the prose/abstract sweep (analogous to round 1's theorem sweep) to drain that layer before cert5.

**Method.** Read all five papers' abstract + §1 Introduction in full (P45 abstract 130–183 / §1 199–423; P46 148–215 / 268–472; P47 119–177 / 213–390; P48 147–231 / 268–630; P49 178–305 / 349–681), checking every metric-level Lorentzian/propinquity *convergence* (or degenerate-metric property) statement for an inline descope/structural tag. (PM-direct, same as round 1; cert5 is the fresh-adversary validation.)

**Findings — 2 genuine of 10 regions:**
| Region | Verdict | Fix |
|---|---|---|
| P45 abstract + §1 | CLEAN | negative paper; K⁺ annihilation + Remark[History and retraction] explicit |
| P46 abstract + §1 | CLEAN | bold "descoped" + "convergence not established" + rate-formula tags throughout |
| P47 abstract | CLEAN | inner/propinquity arrow tagged descoped inline; outer/norm-resolvent arrow (survives) correctly scoped; §1 fixed in cert4 |
| **P48 abstract B2** | **FIX** | "Main theorem...holds at theorem-grade rigor:\ (B2) reverse triangle inequality" listed B2 flatly — but its own §1 (l.385–411) honestly calls B2 "structurally vacuous off-orbit." Abstract was *less honest than its own §1*. Restructured: B1/B3 structural (theorem-grade); B2 structural/vacuous; B4 descoped (B4 already tagged in cert4). |
| P48 §1 | CLEAN | §1.2 overview already partitions structural (B1,B3) vs descoped (B4); T6/wedge-G2 descoped inline |
| **P49 §1 l.510** | **FIX** | honest-scope bullet "Bridge Theorem 6.4′-Q1′ holds at theorem-grade rigor on the enlarged substrate" lacked the per-leg qualifier its abstract carries. Added "for its structural and thermal-time-stack legs (B1′/B2′/B3′); the metric-level convergence leg B4′ is descoped." |
| P49 abstract | CLEAN | B4′ descoped; B2′ explained as surviving state-level thermal-time-stack (D_max chain); "generic" strictness + illustrative-Umegaki labels held |

**Verification.** p48/p49 compile errors=0, zero undefined cites/refs; C5/C11/C13/C14/C15 group1 PASS.

**Honest scope.** Corrective tagging sweep, not a /qa run (no seeds/calibration). Claim: every abstract + §1 statement in the cluster asserting a metric-level Lorentzian/propinquity convergence (or degenerate-metric property) now carries an inline descope/structural tag, so a fresh statement-by-statement read of the front matter should find no new prose/abstract zombie. Combined with round 1 (theorem statements), both the theorem-statement layer and the prose/abstract layer are now swept. If cert5 still finds a descope zombie, it is in the paper *body* (post-§1 prose), the remaining untouched layer.
