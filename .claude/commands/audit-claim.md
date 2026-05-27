---
description: Curve-fit-audit pattern — before letting a numerical match become a structural claim
---

Audit the claim just made (typically a numerical coincidence, a PSLQ relation, a structural match, or a "this looks like that" pattern). Use the protocol that caught W3 and c₂ retrospectively. Reference:\ `docs/curve_fit_audit_memo.md`.

**Step 1. Restate the claim cleanly.** One sentence. What is being matched? At what precision? Against what reference?

**Step 2. Free-parameter count.** How many real numbers does the claim explain? How many small-integer free parameters did the claim need (PSLQ coefficient ceiling, basis size, anything tuned)? **If the free parameter count is comparable to or exceeds the explained-real-number count, the claim is curve-fit territory.** State the ratio explicitly.

**Step 3. Selection-bias check.** Was the basis frozen before the test, or was it expanded until the relation closed? Was a single-data-point hit treated as evidence, or was the relation tested on independent data the basis was not chosen for? If the relation came from a brute-force basis sweep on one data point, mark it as **selection-attributable until tested on a frozen-before-data second point**.

**Step 4. Alternative explanations.** List at least three. What else could produce the observed numerical match? (Accidental rational approximant, related-but-distinct identity, generic property of the ring being searched, numerical artifact at the precision tested, etc.) For each alternative:\ how would you tell it apart from the structural claim?

**Step 5. Robustness perturbation.** What happens if you nudge a parameter slightly (precision, basis ceiling, cutoff)? If the relation evaporates under perturbations smaller than the test precision, the relation is fragile and the claim is at risk.

**Step 6. Independent test.** What is the cheapest independent computation that would either reinforce or kill the claim? Schedule it. **If you cannot name an independent test, the claim is not yet ready to become structural.**

**Step 7. Verdict.** GO (claim survives) / BORDERLINE (claim survives audit but does not yet warrant paper-level structural identification) / STOP (claim is selection-attributable or fragile; demote to numerical observation).

**Discipline.** The framework has four irreducible-but-natural constants outside the master Mellin engine (K, Wolfenstein, S_full(GS), L2 next-order c). Not every PSLQ-relatable value is a structural identity. The cost of the audit is small; the cost of building a paper on a false-positive is enormous.
