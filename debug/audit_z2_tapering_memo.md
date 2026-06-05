# Audit — Sprint Z2 tapering (Hopf-U(1) m -> -m systematic application)

**Date:** 2026-06-04
**Auditor role:** curve-fit-audit per `feedback_audit_numerical_claims.md`
**Sprint memo audited:** `debug/sprint_z2_tapering_memo.md`
**Driver:** `debug/sprint_z2_tapering.py`
**Data:** `debug/data/sprint_z2_tapering.json` (37 entries)
**Audit verdict:** **CONFIRMED with caveats** — the headline "uniform ΔQ = −3 across 37/37 molecules with machine-precision spectrum preservation on every system small enough to diagonalise directly (Q ≤ 20)" is supported by the artefacts, but several pieces of small print should be reflected in the sprint's framing before any paper claim is made.

---

## 1. Free-parameter count

How many tuneable knobs were dialled to reach ΔQ = −3 with machine-precision δE?

| Knob | Set how / where | Free? |
|:-----|:---------------|:------|
| Choice of symmetry to taper (αβ parity + Hopf-U(1) m→−m) | Three Z₂'s identified up-front, all three applied uniformly | **No** — set before any numeric run |
| Basis rotation U (sym/antisym combinations) | Single canonical choice $\phi_\pm = (\phi_{+m} \pm \phi_{-m})/\sqrt2$ | **No** — orthogonal, dictated by the reflection generator |
| Pairing key for orbitals (block label vs block index) | Implementation bug found and fixed mid-sprint | Quasi-free; the fix is structurally correct (across-atom pairing is geometrically wrong for homonuclear cases) |
| Sector sign $(s_1, s_2, s_3) \in \{\pm 1\}^3$ for tapering | Swept all 8 sectors, minimum-energy chosen as GS | **8-way per molecule**; but the choice is justified post-hoc by ground-state minimisation, not a fit |
| Tolerance for "spectrum preserved" (1e-8) | Hard-coded in driver | Free; loose by design (actual δE ≤ 1.2e-11) |
| max_n = 2, R = default | Production defaults | Reasonable; one robustness shift run by this audit (max_n=3) |

**Verdict:** ΔQ = −3 is a structural consequence of the construction once you accept three commuting Z₂'s — not a fit. The Z₂'s themselves are forced by (a) electron-number conservation, (b) S_z conservation, and (c) the Paper 29 §5.3 Observation 5.1 m_l → −m_l reflection (which is the only sub-action of Hopf U(1) that commutes with a real-integer adjacency). The free-parameter count for the headline ΔQ result is effectively zero. The 8-sector sweep is correct ground-state determination, not a fit.

The **only genuinely empirical** decision was the block-index pairing key — which was caught by a spectrum-preservation gate (N₂/F₂ failed P_commutes when label-only was used). This is a feature, not a fit.

## 2. Selection bias

Did the sprint cherry-pick the molecules that match?

- **No.** The sprint ran every molecule in `ecosystem_export._SYSTEM_REGISTRY` (37 systems), reported 37/37 = pass, and the audit confirmed registry size = 37 (matches sprint claim).
- Three CLAUDE.md §1.5-listed multi-center molecules (CH₂O, C₂H₂, C₂H₆) are excluded; the audit confirmed that no `ch2o_spec` / `c2h2_spec` / `c2h6_spec` exists in `geovac.molecular_spec` or `geovac.composed_qubit`. The exclusion is forced by the codebase, not by selection. **However, the CLAUDE.md §1.5 entry claiming these are in the library is stale** — surfacing this is a side-finding the sprint correctly flagged.
- **One legitimate ambiguity:** the sprint says "37/37 = 132%" against a stated 20/28 decision gate. The 28 referred to Paper 14's composed table cohort; the actual run is 37 because it includes everything the registry holds (TM hydrides + alkaline-earth + H2/He). The 132% is a cosmetic overshoot of a slightly mismatched denominator, not a substantive issue, but the "decision gate" framing is loose.

**Verdict:** No selection bias. The 28 vs 37 denominator mismatch is cosmetic.

## 3. Independent test

Is there a held-out check that was not part of the design?

- **Spectrum preservation** is the natural held-out check, and it is genuinely held out: the construction (orthogonal rotation U + stabiliser projection) is built to commute with H by symmetry; the δE comparison to direct diagonalisation is the independent verification.
- **Six of 37 molecules were spectrum-checked** (Q ≤ 20): H2, He, NaH, KH, SrH, BaH. All six passed (δE between 5e-16 and 1.2e-11 Ha — machine precision through Lanczos tolerance).
- **31 molecules were unverified.** The sprint argues this is OK because the construction is a similarity transform. **This argument is rigorous** — orthogonal basis rotation preserves spectrum exactly; tapering onto the +1 eigenspace of a true H-commuting stabiliser gives the spectrum restricted to that sector. The numerical confirmation P_commutes=True at every n_max ∈ {2,3} for 37 molecules (audit's max_n=3 spot-check + sprint's max_n=2 run) gates the construction.
- **A stronger held-out check would be:** Q=22-26 molecules (HCl, HBr, NaCl, NH3? — at slightly smaller n_basis). The driver's gs_cutoff defaults to 22, and Q_naive=20 systems verify in ~75 s each. The audit attempted Q=28 spectrum checks at max_n=3 (H2, He) but they were too expensive to complete in the audit window (2^28 ≈ 268M dim sparse matrix).

**Honest gap to name:** no spectrum verification at Q in {30..50} (LiH, BeH2, TM hydrides). The structural argument is sound but reads stronger if even one Q=30 case is verified. Q=30 sparse-Lanczos is feasible (~10-30 min); follow-on recommendation.

## 4. Robustness

Audit ran the driver at one shifted parameter (max_n=2 → max_n=3) for H2, He, NaH:

| Mol | max_n | M | Q_naive | Q_tapered | ΔQ | P_comm |
|:----|:-----:|:--|:-------:|:---------:|:--:|:------:|
| H2  | 3 | 14 | 28 | 25 | 3 | True |
| He  | 3 | 14 | 28 | 25 | 3 | True |
| NaH | 3 | 28 | 56 | 53 | 3 | True |

ΔQ = −3 and P_commutes = True persist at max_n=3 for all three. Pauli term count is identical pre/post (Pauli_naive = Pauli_tap), 1-norm identical — consistent with the sprint's observation that "for Q ≥ 20 the Pauli term count and 1-norm are unchanged by tapering" (the merging happens only at Q=10).

Spectrum verification at max_n=3 (Q=28) was attempted by the audit but did not complete in the audit window. The structural similarity-transform argument carries.

**Verdict:** Robust at max_n=3 on three test molecules. The construction is basis-intrinsic (Gaunt selection rule), so it should hold at every max_n.

## 5. Honest scope

Does the verdict overclaim?

**Two scope tightenings recommended:**

1. **"Machine-precision spectrum preservation" is honest only for Q ≤ 20.** The sprint says "≤ 1.2 × 10⁻¹¹ Ha for every system small enough to verify directly". This is correct. But the abstract framing "37/37 with machine-precision spectrum preservation" is slightly tight; honest reading: 6/37 verified to ≤ 1.2e-11 Ha, 31/37 trusted by structural similarity-transform argument. The memo §3 footnote captures this honestly; the abstract/headline should mirror it.

2. **The ΔQ count is the metric; Pauli/λ reduction is small or zero at Q ≥ 20.** The sprint is honest about this in §3 ("Cross-family pattern") — for Q ≥ 20 the Pauli count and 1-norm are unchanged, only qubit count drops by 3. A reader skimming the headline might infer Pauli reduction; the body correctly says no. Paper 14 update text should foreground "qubit count drops by 3, Pauli count and 1-norm essentially unchanged for production-size molecules".

3. **The "GeoVac-specific advantage over Gaussian encodings" framing in §4 needs care.** The claim is "Gaussian basis on the bond axis has only axial rotation symmetry; the m_l-reflection is a property of the natural-geometry (n,l,m) basis". This is correct for atom-centred Cartesian Gaussians but a Gaussian basis built on spherical harmonics also has m → −m reflection symmetry. The honest difference is that GeoVac's full (n,l,m) basis is **the** angular-momentum eigenbasis (no truncation, no AO orthogonalisation across centres in standard composed builds), so the symmetry survives bit-exact, while general Gaussian multi-center bases break it via Löwdin orthogonalisation (which mixes m-channels). This nuance should land in the Paper 14 text rather than the blanket "Gaussian only has axial rotation".

## Track-specific check (Z2 tapering)

The audit prompt asks specifically:
**Did the sprint actually verify spectrum preservation to ≤ 1e-10 Ha relative, or just count qubits?**

| Mol | Q_naive | |ΔE| (Ha) | |ΔE|/|E| (relative) | Verdict |
|:----|:-------:|:----------:|:------------------:|:-------:|
| H2 | 10 | 5.0e-16 | ~2e-15 | OK |
| He | 10 | 4.0e-15 | ~1.4e-15 | OK |
| NaH | 20 | 1.1e-13 | ~7e-16 | OK |
| KH | 20 | 3.3e-12 | ~5.5e-15 | OK |
| SrH | 20 | 1.2e-11 | ~3.8e-15 | OK |
| BaH | 20 | 1.0e-11 | ~1.3e-15 | OK |

Relative δE/|E| is ≤ 6e-15 for all six verified cases — significantly tighter than the 1e-10 threshold. The driver's tolerance is `1e-8` absolute, which is much looser than what the construction actually delivers.

**Are there molecules where it silently failed?**

No silent failures. The block-index fix was caught by the spectrum-preservation gate (memo §1.5 — N2/F2 first pass failed P_commutes; fix re-ran clean). The sprint records 37/37 P_commutes=True after the fix.

**One pattern that needs naming:** the sprint reports `delta_E = None, spectrum_preserved = None` for 31 of 37 systems because they are above the gs_cutoff. This is honestly reported as 'skip' in the table. The structural argument fills the gap.

## Pass / fail list

| Item | Status |
|:-----|:------:|
| Free-parameter count = 0 for ΔQ | PASS |
| All molecules in registry tested | PASS |
| No cherry-picking | PASS |
| Spectrum preserved to < 1e-10 Ha rel. on every checked system | PASS (actual: < 6e-15) |
| Held-out check exists (spectrum diag) | PASS (limited to Q ≤ 20) |
| Robust at max_n=3 on H2/He/NaH | PASS (P_commutes carries) |
| Honest scope on "machine-precision spectrum preservation across 37/37" | NEEDS_TIGHTEN (verified on 6/37, others trusted by structural argument) |
| Pauli / 1-norm reduction claims | PASS (sprint is honest in body; headline should mirror) |
| Gaussian-comparison framing | NEEDS_TIGHTEN (spherical-Gaussian basis also has m→−m; honest difference is AO mixing) |
| Block-index pairing fix | PASS (correct structural fix, caught by gate) |

## Recommended modifications to verdict

The sprint's "GO" verdict stands. Two small additions:

1. **Reword headline** from "machine-precision spectrum preservation for every system" → **"machine-precision spectrum preservation (δE/|E| ≤ 6 × 10⁻¹⁵) verified directly on 6/37 systems (Q ≤ 20); the construction is a similarity transform onto the +1 eigenspace of an H-commuting stabiliser, so the spectrum is preserved by construction on the remaining 31/37."** This is what the body of the memo already says; the headline just needs to mirror it.

2. **Reword Paper 14 insertion** to replace "Gaussian encodings do not [admit this Z₂]" with "Gaussian encodings break this Z₂ structurally via the multi-center AO orthogonalisation that mixes $m$-channels; only an angular-momentum eigenbasis that does not cross-center-mix (such as GeoVac's natural composed build) preserves it bit-exact." The Paper-29 / Hopf framing is unaffected.

3. **Run one Q=30 spectrum verification** (e.g. LiH, Q_naive=30, single n=27 sector). ~10-30 min wall, would extend the held-out class from "Q ≤ 20" to "Q ≤ 30" and remove the last residual scope tightening. Optional; the structural argument is sufficient.

## Audit verdict line

**CONFIRMED with caveats because** the construction is rigorous (orthogonal basis rotation + projection onto +1 eigenspace of three commuting Z₂'s), the spectrum-preservation gate runs cleanly at machine precision (δE/|E| ≤ 6 × 10⁻¹⁵) on all 6 verified systems, P_commutes = True for all 37/37 molecules at max_n=2 and audit-spot-checked at max_n=3, and the ΔQ = −3 headline carries zero fitted parameters; the only caveats are that 31 of 37 cases rely on the (correct) similarity-transform argument rather than direct numerical verification, the Pauli/λ reduction is small or zero for production-size molecules (the qubit count is the load-bearing metric), and the "Gaussian encodings cannot do this" framing should be sharpened to "Gaussian encodings break it via AO orthogonalisation, GeoVac preserves it because it does not cross-center-mix the angular eigenbasis".
