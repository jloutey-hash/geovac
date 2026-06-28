# /qa group2 — re-cert run #2 (2026-06-27): FAIL → remediated

Canonical record for the second `/qa group2` invocation (confirmation re-cert of the
v4.50.0 remediated corpus). Run #1 (v4.50.0) = FAIL→remediated; this run #2 peeled a
**new layer** of genuine defects (the expected fresh-adversary "one scope over" pattern;
group1/group3 each took several cycles).

## Verdict: FAIL (trustworthy)

All five dimensions exercised + calibrated; verified material defects in every
calibrated LLM dimension → FAIL. Deterministic C10–C16 PASS.

### Calibration scorecard
- **Sensitivity 8/9.** Caught rs1 (P11 "derives exact from first-principles"),
  rs2 (P17 "cc-pVTZ-competitive"), rs3 (FCI-A "matching exact"), rs4 (P15 kolos1968
  vol/year), rs5 (P19 bunge1993 venue-swap), rs6 (test_neumann_vee `lhs-lhs`),
  rs7 (test_composed_diatomic `errs>=0`), rs9 (synth "matching the exact D_e").
  **Missed rs8** (test_fcim `E<-7.5`→`E<0`) — code-FCI-M didn't flag the vacuous guard
  ⇒ FCI-M code strictly INCONCLUSIVE this run (but the dimension found genuine material
  defects elsewhere). Next re-cert: sharper code-FCI-M prompt (enumerate each
  assertion's discriminating power).
- **Specificity 6/6.** No control (M1–M6) flagged material.

### Per-dimension scorecard
| Dimension | Exercised | Calibrated | Clean | Verdict |
|:--|:--|:--|:--|:--|
| Code/test (C1/C2) | yes | partial (rs8 missed) | no | FAIL |
| Claims (C3/5/6/8) | yes | 3/3 | no | FAIL |
| Citations (C4) | yes | 2/2 | no | FAIL |
| Synthesis (C9) | yes | 1/1 | no | FAIL |
| Deterministic (C10–C16) | yes | n/a | yes | PASS |

## Verified material findings + disposition

**Investigate-before-downgrade overturned TWO apparent downgrades (again):**
- **H2-rovib (P13 §IX):** code-P13 claimed 4435/4157 not reproducible (only-runnable
  artifact gives +10.5%). FALSE — the +11.7%/+10.5% is a Morse **fit-range artifact**:
  the SAME Neumann PES (D_e=0.161, 92%) gives ω_e=4396 (−0.1%) on the near-minimum fit
  [1.0,2.0] (§IX recipe) vs 4918 (+11.7%) on the wide fit [0.8,6.0] (bad dissociation
  tail). §IX is reproducible. My initial "honesty downgrade" edits were **reverted**;
  the genuine gap (no pytest) is closed by `tests/test_h2_rovib_morse.py` (3✓) +
  corrected `results.md` note.
- **P13 cusp "properly variational":** code-P13 (computed) vs claims-B (reasoned)
  conflict. Recompute confirms code-P13: E_corrected=−2.90383 < exact −2.90372
  (non-variational, 0.0035%). Applied: restrict "variational" to the raw 0.022%; label
  cusp 0.004% a non-variational extrapolation (3 loci) + test now asserts
  `E_corrected < E_EXACT`.

**LARGE — applied:**
1. **TC qubit re-assertion (P15 + P17)** [C5/§3]: "TC validated, eliminates 5.3%→8.2%
   divergence" re-asserted the Track BX-3 qubit-space-diag **false positive**
   (corrected in TC-V; standard projected FCI converges 5.3%→2.0%, better than TC's
   3.4% plateau). Reframed both papers to the §3 honest negative. **Added C16 registry
   entry `tc-qubit-validated-success` (group2, fail-severity).**
2. **P15 HeH⁺ "parity"** [C8(d)]: 93.1% is the *adiabatic* (disavowed) solver. **2D
   variational recompute (same grid as H2 94.1%): E_total=−2.9077 (variational) → D_e
   0.00396 → only 5.3% of D_e** (abs error 0.071 Ha ≈ the entire 0.075 Ha binding).
   Adiabatic 93.1% is a non-variational over-binding coincidence; parity withdrawn;
   HeH⁺ converges far slower than H2. Driver `debug/heh_2d_recompute.py`.
3. **FCI-A LiH "equilibrium minimum near R=2.5"** [C8(b)/C5]: contradicted the FCI-M
   no-equilibrium guardrail. Reframed to monotonically-attractive / no interior minimum
   (R-independent kinetic, not basis incompleteness).

**LARGE — flagged to PI, NOT auto-patched:**
4. **P8 Sturmian structural theorem** [C1]: displayed eq (1/β_j−1), proof (1−β_j⁻¹,
   opposite sign), code (1−β_j) are three inconsistent coefficients; numerics show
   neither the pure proportionality nor (fixed-β) R-independence hold for the
   `compute_h1_matrix` general two-center construction — but that is the binding-capable
   object (Corollary 1), NOT the degenerate single-center bond-sphere object the theorem
   isolates (the v0.9.34 diagnostic, likely archived). Correct form needs careful
   re-derivation; patching risks making a guardrail theorem more wrong. Headline
   (single-center Sturmian can't bind without R-dependent β) robust via Track DF + §3.

**§13.5 judgment (PI):** FCI-molecules fitted λ (Fock-weighted correction) — KEEP as a
labeled qualitative diagnostic (it demonstrates the negative; it still fails to
reproduce both R_eq and D_e). PI-directed.

**SMALL — applied:** P11 abstract −0.6025→−0.6026; FCI-A graph-native 0.22%→0.25% (n_max=5);
P17 conclusion "2D reduces drift 4×" deleted (v2.0.32 retraction); P13 Abdouraman2016
PRA 94,023403→J.Phys.B 49,235005; FCI-M eq:cp composition note (free-atom rows → D_e^raw;
ghost rows → D_e^CP); synthesis Li 1.07% n_max=4→n_max=5; P12 "25 tests"→32.

## Test backfills (all green)
- `tests/test_h2_rovib_morse.py` (NEW, 3✓) — near-min vs wide-range Morse fit; pins ω_e=4396 (−0.1%).
- `tests/test_level3_variational.py::test_cusp_correction_sub_01pct` — now asserts raw>exact (variational) AND corrected<exact (non-variational).
- `tests/test_l_dependent_pk.py::test_r_eq_comparison` — pins l_dependent err < channel_blind err and <7% (was NO-TEST).
- `tests/test_fcim_lih_molecular.py::test_lih_pes_monotone_no_interior_minimum` (NEW, 1✓ @slow) — discriminating no-equilibrium R-scan.

## Verification
- Deterministic C10 (compile): all 8 touched papers content-err=0 modulo missing-figure
  draft-mode (environmental); no new undefined refs. C11–C16 PASS (C16 includes the new
  group2 TC entry, clean).
- Worktree removed; leak-scan clean.

## Next
- Release the remediation (patch); re-run `/qa group2` (run #3) for the certified PASS.
- PI decision needed on P8 theorem-formula (raise-to-PI item).
