# Audit: Sprint W1e period-class diagnostic

**Auditor:** adversarial review per `memory/feedback_audit_numerical_claims.md`.
**Subject:** `debug/sprint_w1e_period_class_memo.md` (2026-06-04).
**Sprint headline:** *STOP because 0/11 W1e correction terms identify with low-coefficient outer-factor periods (M1/M2/M3); W1e lives in the inner-factor input-data tier (Paper 18 §IV.6 chemistry-side analog).*
**Audit verdict (one line at end): CONFIRMED PARTIAL.**

---

## §1. Free-parameter count

The sprint headline is a NEGATIVE result — "zero M1/M2/M3 hits across 11 terms." A negative result has very few free parameters in the headline itself:

- 4 basis choices (M1 / M2 / M3 / INNER).
- 1 PSLQ ceiling (audit C=100; permissive C=10⁶).
- 3 audit filter classes (rational-only, basis-internal, CFA-failure).
- 11 target terms.
- mpmath dps=100.

The decision to **exclude the constant `1` from M1/M2/M3 bases** is a non-trivial parameter, defended in §1.2 of the memo. I accept this defense:\ a basis-with-`1`-and-no-transcendental-content would always find a low-denominator rational fit at 100 dps against a float64 input (3-4 significant digits), which would not be a period identification. The exclusion makes the test *stronger* against the framework, not weaker.

The INNER basis has **6 generators built from only 3 independent CR exponents** ($\zeta^2/2n^2$ relations). The sprint correctly identifies this as basis-overcompleteness and uses it as a filter, not as evidence. This is a legitimate methodological choice, not a free-parameter inflation.

**Free-parameter audit verdict:** Acceptable. The negative-result headline is robust to reasonable choices of basis ceiling and filter thresholds.

---

## §2. Selection bias

Four bases were tested in parallel with no a-priori preference, and the sprint reports the result for *all four* including INNER, with a defended filter classification. This is methodologically clean — no cherry-picking.

**However**, there is a hidden selection in the *targets*:

- 11 terms were extracted from sprints F4, F5, F6 (the most recent ones).
- Sprints F1, F2, F3 (earlier W1c/W1d closure attempts) had analog terms (multi-zeta differential, kernel-shape corrections) that were *not* extracted.
- The hypothesis "W1e lives in M3" could have been tested on a larger set, but the sprint chose the F4-F6 subset because they are the latest and best-characterized.

I do not flag this as a problem — it is reasonable scoping — but the negative result is "0/11 from the F4-F6 subset" not "0/many from the whole W1c-W1e arc." If a single F1 or F2 term *had* identified in M3, the conclusion would change. The memo §5 honest-scope correctly names cross-system robustness as a follow-on but does NOT name "cross-sprint robustness within NaH."

**Selection-bias verdict:** Mostly clean; minor scope-creep that the memo does not flag.

---

## §3. Independent test

**Held-out check that was not part of the fitting process:** None in the sprint itself. The decision-gate ("STOP if zero M1/M2/M3 hits") was set up as a null hypothesis test, which is good practice, but every term was part of the fit.

**An honest independent test would look like:**

1. **Random-rational null hypothesis on each basis** — I ran this: 50 random rationals in [0.001, 1000] against each basis at ceiling 100:
   - **M1, M2, M3: 0/50 false positives each** (the "0 hits" result is genuinely informative — the bases are sparse enough that random rationals do not fit).
   - **INNER: 16/50 false positives** (32% noise rate).

This null-hypothesis test was *not* included in the sprint but it strongly supports the headline: zero hits in M2/M3 against the 11 terms is not noise (would expect 0/11 even on random data). Zero hits in M1 likewise. The INNER "hits" are noise-consistent.

2. **The memo's claim that 0.194 = 97/500 is a 'rational-only M2 artifact'** is *wrong as stated*. Running PSLQ on `[0.194, π², π⁴, π⁶]` at maxcoeff=100 returns **None**. The memo §2.1 table claims three "rational-only M2 artifacts" (0.194, 4.374, 0.0746) but the JSON shows **zero M2 hits at any ceiling** and `spurious_rational_hits_count: 0`. This is a memo-vs-JSON inconsistency. The memo's auditor-friendly story about "0.194 = 97/500 is an input-precision artifact" is *true in principle* but the M2 hits the memo describes do not actually exist. I flag this as a **memo writing error** — the structural conclusion is unchanged (no M2 hits), but the cited evidence is fictional.

**Independent-test verdict:** Sprint did not include null-hypothesis test; I ran one and it confirms the headline for M1/M2/M3. **Memo §2.1 contains a factual error**: the "3 rational-only M2 artifacts" table entries do not exist in the JSON; the actual M2 hit count is 0, not 3.

---

## §4. Robustness

**Robustness check I ran (1 minute):** Shift `zeta_1s_CR` by +0.1% (10.6259 → 10.63654) and re-run PSLQ on `F4_PK_barrier_Ha = 0.194` against the INNER basis at ceiling 100.

- Original relation: `[11, 33, 20, -25, -9, -8, 41]` (target coefficient 11)
- Perturbed relation: `[-38, -35, -14, 2, 0, 72, 12]` (target coefficient -38)

**Completely different relation** under a 0.1% input shift. This is the *signature of a tautological hit*: when the relation depends sensitively on input value rather than encoding a structural identity, the hit is overfitting noise. The memo's CFA-failure-mode flag was correct; my robustness check confirms it independently.

**The memo's CFA arithmetic** (`p/(n-1) < log10(C)`) is the standard PSLQ noise-floor heuristic. For INNER with n=6, p=3-4 digit precision, C=100: threshold is `4/5 = 0.8 < log10(100) = 2`. Violated, as the memo states.

**Cross-system robustness:** The sprint explicitly declines to test LiH / MgH₂ / H₂O analogs. This is named as a follow-on. Without these, the structural conclusion is supported by *one* system. The memo §5 honest-scope is appropriately conservative.

**Cross-precision robustness:** All 11 inputs are float64 from prior JSONs. The memo §5 notes that 100-dps recomputation would not change the *structural* reading (FCI eigenvalues are algebraic-implicit even at infinite precision), but it would tighten the M1/M2/M3 ceiling argument. **The memo does not run this check.**

**Robustness verdict:** The negative result is robust to basis-ceiling choice (audit C=100 vs permissive C=10⁶ give the same M1/M2/M3 verdict). The INNER hits are demonstrably tautological. Cross-system and cross-precision robustness are explicitly deferred.

---

## §5. Honest scope

The memo's verdict line says:

> *STOP because 0/11 W1e correction terms identify with low-coefficient outer-factor periods (M1, M2, or M3) after audit filtering; period-class framing was the wrong axis — W1e corrections are FCI eigenvalues and Hartree integrals built from external chemistry-side calibration data (Clementi-Raimondi exponents), structurally in the Class 1 / inner-factor input-data tier (Paper 18 §IV.6 chemistry-side analog).*

Two parts to this verdict:

**(a) Numerical claim:** "0/11 W1e correction terms identify with low-coefficient outer-factor periods at audit ceiling 100."
This is CONFIRMED by my re-check. The bases are clean (M2/M3 give 0/50 false positives on random rationals); the inputs at 100 dps mpmath are exact rationals (no float noise to suppress); the basis-without-constant choice is defended.

**(b) Structural claim:** "Period-class framing was the wrong axis — W1e lives in the inner-factor input-data tier (Paper 18 §IV.6 chemistry-side analog), categorically disjoint from M1/M2/M3."
This is **PARTIALLY supported, mildly overclaimed**. The negative result (0/11 outer-factor hits) is *consistent with* the inner-factor-tier classification but does not *prove* it. Logically:
- "No M1/M2/M3 hit" rules out *one* hypothesis (W1e is an outer-factor period).
- It does NOT rule out *other* hypotheses (W1e has higher-weight cyclotomic content, W1e involves periods of *different* algebraic varieties not yet built into the master Mellin engine, W1e involves M2/M3 *combinations* at higher coefficient ceiling than 100).

The structural argument the memo *does* give — "M3 lives at k=1 vertex-parity, W1e has zero vertex-parity content" — is *independent of the PSLQ negative* and is the actual load-bearing argument. The PSLQ run is a *consistency check* on this structural argument, not a derivation of it.

**Recommended tighter scope:** "The PSLQ result is consistent with the structural argument (M1/M2/M3 are *a priori* the wrong sector for FCI eigenvalues over fitted hydrogenic orbitals). It does not constitute an independent proof; the structural argument carries the verdict, and the PSLQ test rules out the loophole 'W1e might *happen* to land in M1/M2/M3 by coincidence.'"

---

## §6. Track-specific check (PSLQ period identification)

Per task prompt: *did the PSLQ pass actually identify periods (low coefficient ceiling, large basis, ratio of basis-element magnitudes), or did it find tautological hits (e.g., the term IS one of the basis elements)?*

| Check | Status |
|:------|:------:|
| Low coefficient ceiling (≤100 for audit run) | **PASS** — primary verdict at C=100 |
| Large transcendental basis (≥3) | M1=4, M2=3, M3=5 — PASS for all three |
| Ratio of basis-element magnitudes within ~10⁴ | M1: π and 1/π² differ by ~31 — OK. M2: π² and π⁶ differ by ~100 — OK. M3: log2 and ζ(5) differ by ~1.5 — tight, near-degenerate but OK. INNER: $\zeta_{1s}=10.6$ and $E_{1s}=56$ differ by 5× — basis has small dynamic range, *bad* for PSLQ discrimination |
| Tautological hit check (target IS a basis element) | **PASS** — none of the 11 targets are in any basis |
| Tautological hit check (target IS Q-linear in basis) | **FAIL for INNER**: 5 of 8 INNER "hits" are basis-internal Q-linear dependencies (target coefficient = 0). Caught by sprint's audit filter, not pre-empted. Filter correctness verified. |
| Robustness under 0.1% basis perturbation | **FAIL for INNER**: relation flips completely (independent verification §4 above). PASS for M2/M3 (no hits to perturb). |

**Track-specific verdict:** The PSLQ pass is **legitimate as a null hypothesis test on M1/M2/M3** and **inconclusive for INNER**. The sprint correctly identifies INNER hits as noise. The headline "0/11 outer-factor hits" survives all checks.

---

## §7. Pass/fail itemization

| Item | Verdict |
|:-----|:-------:|
| Free-parameter count reasonable for a negative result | PASS |
| Selection bias minimal (4 bases tested in parallel) | PASS |
| Within-system target selection (F4-F6 only) not flagged | minor MISS (not in honest-scope §5) |
| Independent random-rational null hypothesis | PASS (I ran it; memo did not) — supports headline |
| Memo §2.1 "3 rational-only M2 artifacts" matches JSON | **FAIL** — JSON has 0 M2 hits, not 3. Memo writing error. |
| Memo internal consistency on filter counts | mixed — §0 table says M2=3 raw hits / 3 spurious_rational; JSON `spurious_rational_hits_count = 0`. Counting convention drift. |
| Robustness under 0.1% INNER basis perturbation | PASS (tautological signature confirmed) |
| CFA failure-mode threshold arithmetic | PASS — `p/(n-1)=0.8 < log10(100)=2` |
| Cross-system robustness deferred | NOTED (named follow-on) |
| Cross-precision robustness deferred | NOTED (named follow-on) |
| Structural argument independent of PSLQ | PASS — M3-at-k=1-vertex-parity vs W1e-no-vertex-parity is the load-bearing arg |
| Verdict scope (numerical claim) | CONFIRMED |
| Verdict scope (structural claim) | mildly OVERCLAIMED — PSLQ rules out one loophole, structural arg is the close |
| Hard-prohibition compliance (no production code or papers modified) | PASS — diagnostic-only, defers .tex edits to PROPOSED §4 |

---

## §8. Recommended modifications to the sprint's verdict

1. **Fix the memo factual error.** §0 and §2.1 cite "3 rational-only M2 artifacts (0.194=97/500, 4.374=2187/500, 0.0746=373/5000)" but the JSON shows zero M2 hits and `spurious_rational_hits_count=0`. The structural conclusion is unchanged, but the cited evidence is fictional. Either: (a) re-run with constant-included M2 basis to *generate* these artifacts and report honestly, or (b) delete the §2.1 claim and replace with "M2 returned zero hits at any ceiling — sparser than M1, consistent with the structural argument that pure π-power identifications would require an explicit graph-spectral mechanism on the chemistry side that does not exist."

2. **Add the random-rational null-hypothesis check.** "PSLQ at ceiling 100 against M2 (3 transcendental generators, no constant) returns zero hits on 50 random rationals in [0.001, 1000]. Against M3 (5 generators): zero hits. Against INNER (6 generators): 16 hits = 32% false-positive rate. The zero-hits-across-11-W1e-terms result is informative for M1/M2/M3 and noise-consistent for INNER." This is a 1-minute computation that strengthens the negative result substantially.

3. **Tighten the structural-claim scope.** Change "period-class framing was the wrong axis" to "*outer-factor* period-class framing was the wrong axis, on the structural argument that M3's vertex-parity sector has no chemistry-side analog. The PSLQ run is a consistency check, not an independent derivation; the structural argument is the close."

4. **Add the 0.1% perturbation robustness sentence.** "INNER 'hits' under 0.1% perturbation of zeta_1s_CR flip from `[11, 33, 20, -25, -9, -8, 41]` to `[-38, -35, -14, 2, 0, 72, 12]` — tautological signature confirmed."

5. **(Optional) Add cross-sprint robustness check within NaH.** Pull F1/F2/F3 analog terms from earlier sprint JSONs and run PSLQ on them too. If still 0 outer-factor hits, the negative result strengthens; if any hit appears, the structural conclusion needs revisiting.

---

## §9. Final audit verdict

The headline numerical claim ("0/11 outer-factor period identifications at audit ceiling 100") is **CONFIRMED** by my independent re-check including a 50-random-rational null test that returns 0/50 false positives on M2 and M3.

The headline structural claim ("W1e lives in the inner-factor input-data tier") is **PARTIALLY supported**: it relies on the *a priori* structural argument (M3-at-k=1-vertex-parity has no chemistry analog) rather than on the PSLQ negative. The PSLQ negative is a clean consistency check, not a derivation.

The memo contains **one factual error** in §2.1 ("3 rational-only M2 artifacts") that contradicts the JSON. The error is in the memo's narrative, not in the underlying computation; the headline survives.

**Audit verdict line: CONFIRMED PARTIAL because the numerical negative is solid (independent null-hypothesis 50-random-rational test gives 0/50 on M2/M3), but the structural claim is carried by the a-priori argument (not by PSLQ), and the memo §2.1 cites fictional "rational-only M2 artifacts" that the JSON contradicts (`spurious_rational_hits_count=0` vs memo's "3").**
