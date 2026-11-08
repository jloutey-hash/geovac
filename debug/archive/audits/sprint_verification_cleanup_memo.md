# Sprint Verification Cleanup — 2026-06-04

## 1. Provenance

- **Trigger.** PI question:\ "we identified low hanging fruit that we didn't follow up on. is there anything like that?"
- **Workflow ID.** `wf_39fc232e-906` (task id `w2lysgmxc`).
- **Wall time.** ~25 min for 12 parallel verification agents; ~10 min of sequential file edits + cleanup.
- **Sub-agent cost.** 2.3M tokens, 313 tool uses across 12 agents.
- **Authorization.** PI explicit:\ "yep, let's just launch 6 at a time tho, so two waves" + "I want to address all the follow ons before I commit."

## 2. Setup

Three Explore-agent scans (math.OA arc, gravity arc, equation-verification §13.4a protocol) identified twelve named theorems with no `tests/` coverage across Papers 34, 45, 46, 51, 55. Cross-checked against `debug/followon_register.md` — confirmed none were tracked. Plus drift in two additional places:\ Paper 51 missing from CLAUDE.md §6 inventory entirely (despite published v3.14.0), and a convention split between Paper 55 §4 closed forms and production code that future readers would hit.

PI authorized two waves of six parallel agents.

## 3. Wave 1 — Paper 51 gravity arc

Six theorems from `papers/group5_qed_gauge/paper_51_gravity_arc.tex`:

| Theorem | Test file | Tests | Runtime | Verdict |
|:--------|:----------|:-----:|:-------:|:-------:|
| `thm:zeta_unit_neg_k` (Bernoulli identity ζ_{S³,R}(-k)=0) | `tests/test_paper51_zeta_unit_neg_k.py` | 40 | 0.56 s | PASS |
| `thm:scalar_ak` (a_k^Δ = 2π²/k!) | `tests/test_paper51_scalar_ak.py` | 10 | 1.09 s | PASS |
| `thm:dS_vacuum` (s_min = -Λ/(12√6)) | `tests/test_paper51_dS_vacuum.py` | 4 | 0.67 s | PASS |
| `thm:j_blindness` (Schur S^(2) on (1,1)) | `tests/test_paper51_j_blindness.py` | 6 | 0.68 s | PASS |
| `thm:cutoff_dep` + `thm:sector_mellin` | `tests/test_paper51_cutoff_mellin.py` | 32 | 0.64 s | PASS |
| `thm:L6` F1–F7 | `tests/test_paper51_L6_full.py` | 22 | 0.63 s | PASS |

Key findings:

- **`thm:zeta_unit_neg_k` agent self-corrected a Mellin-residue Jacobian error during its first draft.** Initial test would have shown 4 failures from an incorrect prefactor (missing the $d(2s-2)/ds = 2$ Jacobian); agent caught, fixed, re-ran. Real falsifier preserved, not a tautology.
- **`thm:L6` gap was smaller than initially claimed.** Existing `tests/test_warped_dirac.py` already covered F1, F2, F3, F4, F6, F7 across 75 tests. Only F5 (cone-Dirac saturation) and explicit L6 theorem content (Lemma L6.3 Gaussian envelope, alpha-uniformity, $\lim_n$ / $d/d\alpha$ commute) was genuinely missing. New file fills exactly that gap (22 tests).
- **`thm:cutoff_dep` paper-table values reproduced bit-exactly** for all three paper cutoffs (Gaussian, polynomial $e^{-x^2}$, sharp); `tab:g8_cutoffs` $G_{\rm eff}\cdot\Lambda^2$ and $\Lambda_{\rm cc}/\Lambda^2$ columns match column-by-column.
- **`thm:j_blindness` Schur eigenvalue equality** verified at $n_{\max} \in \{1, 2\}$ across $J = 0, 1, 2$ sectors at $10^{-6}$ tolerance. Multiplicity ratios 1:1:3 (within-sector after Hermitisation) and 1:3:5 (cross-sector) match Paper 51 §6.1 Table.

## 4. Wave 2 — Periods + Lorentzian + Paper 34

| Theorem | Test file | Tests | Runtime | Verdict |
|:--------|:----------|:-----:|:-------:|:-------:|
| Paper 55 `thm:m1_pure_tate` | `tests/test_paper55_m1_pure_tate.py` | 9 | 0.87 s | PASS |
| Paper 55 `thm:m2_mixed_tate` | `tests/test_paper55_M2_mixed_tate.py` | 18 | 2.84 s | PASS |
| Paper 55 `thm:m3_cyclotomic_mixed_tate` | `tests/test_paper55_m3_cyclotomic_mixed_tate.py` | 19 | 0.69 s | PASS |
| Paper 45 `thm:main` asymptotic $4/\pi$ | `tests/test_paper45_asymptotic_rate.py` | 9 | 215 s (slow) | PASS |
| Paper 46 $C_3^{\rm op} = \sqrt{1 - 1/n_{\max}}$ | `tests/test_paper46_C3op_closed_form.py` | 10 | 8.94 s (slow) | PASS |
| Paper 34 6-of-28 projection spot-checks | `tests/test_paper34_projection_spot_checks.py` | 45 | 1.14 s | PASS |

Key findings:

- **Paper 55 M2 PSLQ negative test** (no $\zeta(3)$, no MZV, no Catalan G at the SD level) verified at 120 dps with `maxsteps=3000`. Two-term exactness of Dirac SD verified independently via Jacobi $\vartheta_2$ modular tail (residual $1.5 \times 10^{-39}$ at $t = 0.1$).
- **Paper 55 M3** verified that Catalan G and $\beta(4)$ are genuine level-4 Deligne–Glanois generators — PSLQ-disprovable in the lower-level basis $\{1, \pi, \pi^2, \ldots, \pi^6, \zeta(3), \zeta(5), \zeta(7)\}$ at 100 dps with `maxcoeff=1000`. Sprint 3 RH-J closed form $D_{\rm even}(s) - D_{\rm odd}(s) = 2^{s-1}(\beta(s) - \beta(s-2))$ verified at $s = 4, 5, 6$.
- **Paper 45 asymptotic memory ceiling.** Full `LorentzianTunnelingPair.build` at $(n_{\max}=5, N_t=9)$ hits ~60 GB; infeasible at (6, 11). Agent substituted the cheap $\gamma_{\rm rate}$ (Paper 38 L2 SU(2) Fejér rate); equivalence verified because `Lambda_bound == gamma_joint_su2` bit-exactly at every buildable panel cell (height_P = 0). $4/\pi$ asymptote verified at $n_{\max}$ up to 500 via the cheap path.
- **Paper 46 $C_3^{\rm op}$** bit-exact at $n_{\max} \in \{2, \ldots, 10\}$ and $\{20, 100, 1000\}$ via sympy symbolic enumeration of the per-harmonic supremum $C_3^{(N)} = (N-1)/\sqrt{N^2-1}$ over the envelope range $2 \le N \le 2n_{\max} - 1$. Boundary values $\sqrt{1/2}, \sqrt{2/3}, \sqrt{3/4}, \sqrt{4/5}$ symbolically exact. Asymptote $1 - 1/(2n_{\max}) + O(1/n_{\max}^2)$ verified via `sp.series`.
- **Paper 34 spot-checks** covered §III.2 (Hopf bundle $\pi$, $2\pi^2$, $2/\pi$), §III.6 (Seeley–DeWitt $a_0 = a_1 = \sqrt\pi$, $a_2 = \sqrt\pi/8$ — cross-consistent with `thm:scalar_ak` after $(4\pi)^{3/2}$ rescaling), §III.7 (CH spinor $|\lambda_n| = n + 3/2$, $g_n = 2(n+1)(n+2)$), §III.8 (Wigner 3j Gaunt identity on 5 triples), §III.14 (rest-mass KG-S³×ℝ rationality + S³ conformal-scalar Casimir = $1/240$ from $\zeta_R(-3)$ Bernoulli), §III.27 (Wick $\sigma_{2\pi}(O) = O$ STRONG_IDENTIFICATION at BW canonical $\beta = 2\pi$, $\kappa_g = 1$).

## 5. Side-fixes (PM edits)

- **`CLAUDE.md` §6 Group 5 table.** Paper 51 row added (was missing entirely from the §6 inventory tables, despite Paper 51 being published v3.14.0 / 2026-05-28). The grep for `paper_51` in CLAUDE.md returned no matches before this edit — drift had silently grown over ~3 months.
- **`papers/group3_foundations/paper_55_periods_of_geovac.tex` §4.** New `rem:m2_convention` inserted after `rem:m2_specialisation`, documenting the volume-normalised-asymptotic vs Vassilevich–Branson–Gilkey curvature-polynomial SD convention split. Both readings produce the same M2 pure-Tate ring; rescaling factor is $(4\pi)^{d/2}$. Production code `geovac/qed_vacuum_polarization.py::seeley_dewitt_coefficients_s3()` returns $a_0 = \sqrt\pi$ (raw); $\dim_S \cdot \Vol(S^3) = 4 \cdot 2\pi^2$ recovers the bare $8\pi^2 \in \pi^2 \cdot \Q$ of Paper 51 line 513. Future readers reproducing M2 verification with production code can now land in the right ring without confusion.
- **`debug/followon_register.md`.** Six consolidated A-entries A9–A14 added at top of Section A:
  - **A9** Paper 34 — 22 of 28 projections still uncovered (named per §III.N).
  - **A10** Paper 40 §sec:open trio (log-power ansatz, rank-3 SU(4), G/K symmetric spaces).
  - **A11** Paper 42 §10 O3 + O4 (Pythagorean formal proof, higher-rank modular Hamiltonian).
  - **A12** Paper 50 §8 (Maxwell partition function, squashed-S³, quantitative propinquity rate).
  - **A13** Paper 53 Q1–Q4 cluster (sharp propinquity constant, interior Berezin rate, boundary-adapted operator system, higher-dim warped carriers).
  - **A14** Gravity arc (G4-6a Q1–Q3 dichotomy, Möbius v3.19.0 sign-convention audit, L6 standalone L1'–L5 write-out, L6 sharp lattice rate + B1 constant).

## 6. Honest scope

**Theorem-grade closed.** All 12 named theorems have numerical verification at framework grade — symbolic bit-exact where possible (Bernoulli identities, $C_3^{\rm op}$ closed form, $a_k^\Delta$, $s_{\min}$, $D_{\rm even} - D_{\rm odd}$ identity), numerical cross-check at $\le 10^{-6}$ tolerance otherwise, PSLQ at 80–120 dps for ring classifications. No FAILs, no ANOMALIES, no NOT_VERIFIABLE.

**Important caveat:\ the tests do NOT prove the theorems from scratch.** They verify that the paper's stated closed forms agree with independent numerical computation, at the precision the paper claims. This is exactly what §13.4a asks for:\ a specification test that pins the equation to a computational ground truth.

**Structural sketch — not closed by this sprint, named in register:**

- Paper 42 O3 Pythagorean diagonal/off-diagonal formal proof at general $n_{\max}$ (A11.O3).
- Paper 45 asymptotic $4/\pi$ verified at the cheap-$\gamma$-rate level only; the full `LorentzianTunnelingPair` construction is memory-bound at $(5, 9)$. Substitution is mathematically equivalent at every buildable cell but a sweep over the full construction at higher $n_{\max}$ is computationally infeasible. Honest reading:\ asymptote verified via the equivalent cheap path, not via direct full-pair sweep.

**Numerical observation only:**

- Paper 51 `thm:cutoff_dep` Mellin-master identity verified to $\sim 10^{-7}$ across 4 cutoffs $\times$ 5 $\Lambda$ values (limited by quadrature, not by the closed form).
- Paper 55 M1 PSLQ at 80 dps over $\pi^n$, $n \in [-4, 4]$ — sufficient given $1/\pi^2$ is the deepest M1 witness in the corpus.

**Named open follow-ons (post-sprint register snapshot):** A9–A14 added this sprint; pre-existing A6-followon ($S_{\min}^{(S^5)}$ PSLQ irreducibility), A.M2-M3-unification (multi-year), B1–B4 (collab outreach), C1–C2 (mechanical). Total open:\ ~30 items across Sections A (24), B (4), C (2).

**What this sprint did NOT verify:**

- The 22 remaining Paper 34 projections (now A9 in register).
- The asymptotic rate for Paper 45 beyond the cheap-$\gamma$ substitution.
- Paper 38, 39, 40 propinquity bounds — already covered by `tests/test_gh_convergence.py` per the third Explore scan, so not re-tested here.
- Inner-factor calibration (orthogonal to all checked theorems; multi-year per A.M2-M3-unification).

**What changed in the corpus's epistemic state.** Twelve published theorems went from "stated, no test" to "stated, verified numerically at framework grade." That is the precise content of the §13.4a closure. It is bookkeeping — the theorems' truth is unchanged. What is changed is that the corpus now has a falsifier installed at each of those equations:\ if the paper's claim drifts in a future edit, the test fails and the drift is caught.

## 7. What's ready for commit

- 12 new test files under `tests/test_paper{34,45,46,51,55}_*.py` (224 tests total, all PASS in the workflow run).
- 3 modified files:\ `CLAUDE.md`, `papers/group3_foundations/paper_55_periods_of_geovac.tex`, `debug/followon_register.md`.
- 1 new memo:\ this file.
- Suggested version:\ v3.50.0 (minor — completed verification-cleanup arc).

## 8. Cross-references

- Workflow run transcript:\ `C:\Users\jlout\.claude\projects\C--Users-jlout-Desktop-Project-Geometric\93eb2141-6f3f-4e33-8315-32f20aeb4b9d\subagents\workflows\wf_39fc232e-906\`
- Workflow script:\ `.../workflows/scripts/verification-cleanup-2026-06-wf_39fc232e-906.js`
- Updated register:\ `debug/followon_register.md` (Last updated:\ 2026-06-04).
- New Paper 55 Remark:\ `rem:m2_convention` after `rem:m2_specialisation` (~ line 605 + 18 new lines).
- New CLAUDE.md row:\ §6 Group 5, after Paper 41 row at line 615.
