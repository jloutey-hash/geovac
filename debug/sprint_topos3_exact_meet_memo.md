# Sprint Topos-3 — the two-center frame meet (exact redo)

**Date:** 2026-07-05 · **PI-directed** ("let's do it now" after the quadrature attempt was parked) · **Verdict: GO — the pre-registered prediction CONFIRMED exactly.** The two-center frame meet is the m-grading algebra; the meet-dimension ladder 14 → 9 → 5 completes.

## 1. The two attempts (the process story first)

**Attempt 1 (parked, same day):** nested 2D mpmath quadrature in prolate coordinates. The machinery was *correct* (reproduced ⟨1s|1s⟩(R) = (1+R+R²/3)e⁻ᴿ to 12 digits in the smoke test) but the *cost* was disqualifying — hours with no partial output — and quadrature is categorically unable to prove the exact zeros a support question needs. Killed at PI call; §3 row recorded (the algebraic-first violation).

**Attempt 2 (this memo, ~minutes):** the Mulliken/Ruedenberg route. In prolate spheroidal coordinates every overlap integrand clears to a **polynomial × e^(−pξ)e^(−qη)** (the classical cancellation — *asserted* in code as a built-in correctness check), and the exact auxiliary integrals A_i(p), B_j(q) give

  S = e^(−p) · (U·e^(q) + V·e^(−q)),  U, V ∈ ℚ exactly,

so **S = 0 ⇔ U = V = 0** (linear independence of e^(±q) over ℚ for rational q ≠ 0; single-exponential rational case at q = 0). Support is exactly decidable. Same result the quadrature couldn't finish in hours: under a minute, exact.

One implementation bug caught by the cross-check before any conclusion was drawn: the exponent substitution was missing its /2 (p = (a+b)R instead of (a+b)R/2), producing an R-dependent mismatch against the classical closed form — fixed, after which the cross-check became an **exact rational identity**: normalized ⟨1s|1s⟩(R) = 7/3, 13/3, 103/12 at R = 1, 2, 7/2, matching (1+R+R²/3) exactly.

## 2. Results (all exact; panel Z′ ∈ {1,3} × R ∈ {1,2} at n_max = 3)

- **m-conservation** is structural (azimuthal integral). Within every m-block, the exact support is **FULL** — not merely connected — in every cell: the same-center rate-coincidence zeros (Topos-2's lemma) do **not** survive displacement.
- **Meet = the m-grading algebra, exactly.** The pre-registered prediction (Paper 57 §sec:open_bohr honest-scope, written before this computation) confirms.
- **The meet-dimension ladder** (n_max = 3): identical frames **14** → same-center different-focal **9** (the (l,m)-grading, Topos-2) → two-center **5** (the m-grading). The meet dimension is a measure of the composition's symmetry.
- **Join obstruction is strictly worse than same-center**: KS-obstructed block dimensions {6, 3} (m = 0 has d = 6; m = ±1 have d = 3) vs same-center's single d = 3 block — the site restatement of the corpus's empirical hierarchy (cross-center composition = the hardest wall; the chemistry record).

## 3. Honest scope

- Panel: Z′ ∈ {1, 3}, R ∈ {1, 2}, n_max = 3, homonuclear and heteronuclear. Full support at every tested cell; the claim is panel-scoped genericity, not an all-(Z′,R) theorem (specific (Z′,R) could in principle produce accidental zeros — but connectivity, which is all the meet theorem needs, is robust to isolated zeros).
- The join = M_d step is the same connectivity + Burnside argument as Topos-2.
- Remaining Family-1 follow-ons: the 17-entry census; the general-parameter vanishing lemma. (The two corpus-faithful *instances* — mixed-exponent and two-center — are now both done.)
- §13.5 untouched.

## 4. Artifacts

Driver `debug/compute_topos3_exact_meet.py` (exact; supersedes the parked `compute_topos3_two_center_meet.py`, kept with its status header as the incident artifact) · data `debug/data/sprint_topos3_exact_meet.json` · pins `tests/test_topos3_exact_meet.py` (3/3, self-contained, exact; the m=0 d=6 block runs in ~2 s) · Paper 57 §sec:open_bohr item 6 + honest-scope update.
