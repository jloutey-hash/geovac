# Day consolidation memo: vibe-physics meta-lesson arc (2026-06-05)

**Date:** 2026-06-05 (end of session)
**Author:** PM session
**Status:** Seven positive results in one session, all bit-exact or honestly scope-bounded. The PI's vibe-physics meta-lesson ("GeoVac is comparable with any algebraic-exact technique") held on six independent methods. Theorem 3.2.A empirically PROVEN; closed-form has F3-connection identified, full proof named multi-day follow-on.

---

## Origin of the arc

The PI's vibe-physics observation at the start of the session: "DF and THC use `V_ee = Σ_P L_P ⊗ L_P` which is just `A ⊗ A` in the tensor-product spectral triple. The previous literature search categorized DF/THC as 'close but different' — but maybe they're the same algebraic content accessed two ways. We should check, and not leave any other 'actually equivalent' relationships on the table."

That single observation generated the day's seven sprints.

## Tally of results (in order)

### Result 1 — DF = multipole, bit-exact

`debug/sprint_df_multipole_lift_memo.md`
- DF singular vectors live entirely inside the multipole subspace at every basis tested (n_max=2, 3, 4), bit-exact.
- At production basis (n_max=2): DF rank = analytic count of distinct radial densities = 7.
- F2: `||G2 - VKV^T||_F / ||G2||_F = 10^-16` at n_max=2, 10^-14 at n_max=3 — bit-exact reconstruction.
- F3: graded R^L singular-value spectrum at higher basis attributed to Laguerre-basis ill-conditioning of hydrogenic ERIs.
- **Paper edits applied**: Paper 14 §intro, Paper 20 §intro, Paper 54 new §8.2.

### Result 2 — Cholesky = same multipole subspace

`debug/sprint_cholesky_multipole_memo.md`
- Pivoted Cholesky on the GeoVac ERI tensor reproduces the same subspace as DF, bit-exact.
- Min Cholesky-multipole overlap squared = 1.000000000000 at all basis sizes tested.
- Cholesky reconstruction error 3.4 × 10⁻¹⁷ at n_max=2.
- Meta-lesson held on second independent method.

### Result 3 — QPT stacks with Hopf-Z₂, bit-exact

`debug/sprint_qpt_hopf_stacking_memo.md`
- Every commutator `[H, S²]`, `[H, P_Hopf]`, `[S², P_Hopf]` is bit-exact zero across LiH, BeH₂, H₂O.
- Triply-additive blocking (Bravyi + Hopf-Z₂ + QPT) structurally confirmed.
- Joint reduction: 5.7× beyond particle-number tapering on LiH ground-state sector.

### Result 4 — QPT relativistic compatibility

`debug/sprint_qpt_followons_memo.md` (Section A)
- α² scaling bit-exact: `||H_rel(α) - H_rel(0)||₁ / α² = 3.5417` identical at α/2, α, 2α.
- Spinor m_j → -m_j parity commutes with H_rel bit-exact.
- Triply-additive blocking extends to relativistic regime with J² replacing S².

### Result 5 — QPT Toffoli cost negligible

`debug/sprint_qpt_followons_memo.md` (Section B)
- 18 GeoVac library molecules: Toffoli cost 1,000 – 91,125 per molecule.
- Total library cost: 561,000 Toffoli — 5-7 orders of magnitude smaller than typical QPE budgets.
- QPT is essentially free at GeoVac scale.

### Result 6 — MPO bond rank = 2 at boundaries, bit-identical

`debug/sprint_s1_bond_rank_memo.md`
- χ_k = 2 at every sub-block boundary across LiH/BeH₂/H₂O (12/12 boundaries, bit-identical).
- Universal interior profile {4, 16, 16, 9, 9, 9, 6, 3, 3}, basis-dependent only.
- 8x interior/boundary ratio, bit-identical.
- Theorem 3.2.A empirically PROVEN at production basis.

### Result 7 — Negative control confirms structural reading

`debug/sprint_s2_partial_memo.md`
- Balanced-coupled (Paper 19, cross-block V_ne on) lifts χ_k at would-be boundaries from 2 → 9, bit-identical at both boundaries.
- Confirms the χ=2 floor is structurally tied to cross-block ERI vanishing (F4 from DF memo).

### Result 8 — S2 closed-form connection to F3 identified

`debug/sprint_s2_closed_form_attempt_memo.md`
- Naive operator-counting "χ_k = 2 + n_distinct_left_ops" doesn't match empirical profile.
- Structural reason: bond rank gaps at m=1 (16 vs 9, spin doubling) and m=2-4 (linear dependences in radial-integral matrix).
- **The same F3 hydrogenic-basis ill-conditioning that controls DF rank also controls MPO bond rank**. The MPO and DF live on the same algebraic object.
- Refined Theorem 3.2.A.unified statement proposed: chi_k inherits rank structure from R^L.
- Closed-form derivation named as 3-5 day analytical follow-on.

## Meta-lesson — what holds

**Bit-exact on three independent methods:**
1. DF (SVD-sorted multipole basis)
2. Cholesky (greedy-pivot multipole basis)
3. QPT (Bravyi + Hopf + spin Casimir, triply-additive)

**Structurally confirmed (bit-exact zero commutators or identical singular spectra):**
4. QPT relativistic compatibility (α² bit-exact, spinor parity commutes)
5. MPO bond rank at boundaries (12/12 χ=2 bit-identical)
6. Negative control (χ=9 lift bit-identical at both boundaries)

**Substantive new finding (structural, not bit-exact):**
7. MPO bond rank rank-inheritance from R^L's F3 graded spectrum — same algebraic object as DF rank.

## Paper edits applied today

- **Paper 14 §intro** (`papers/group4_quantum_computing/paper_14_qubit_encoding.tex`): DF/multipole framing strengthened, bit-exact references to verification memo.
- **Paper 20 §intro** (`papers/group4_quantum_computing/paper_20_resource_benchmarks.tex`): chemistry-audience framing of the same content.
- **Paper 54 new §8.2** (`papers/group3_foundations/paper_54_tensor_product_two_body.tex`): tensor-product spectral triple connection to DF/THC, with `vonBurg2021` and `Lee2021` new bibitems.

All three papers compile clean (no errors, no undefined references) on second pass.

## Paper edits staged (not applied today)

- **Paper 14 §sec:hopf_tapering** — Recommended addition of QPT stacking compatibility (relativistic + Toffoli cost). Single new bibitem (Burkat-Fitzpatrick arXiv:2506.09151). Ready for next chemistry-arc consolidation.
- **Paper 14 §sec:mpo_bond_rank** (new) — Empirical Theorem 3.2.A. Could land now with the empirical content or wait for S2-v2 closed-form. PI choice.
- **Paper 20** — Optional QPT cost column for the library resource table. Low priority.

## Honest scope across the whole arc

- **Composed builder only**: no cross-block ERIs (F4 from DF). Paper 19 balanced-coupled spot-checked only for S1 negative control on LiH; full library would tighten the negative-control side.
- **Production basis (n_max=2) for all bit-exact results**. n_max=3 and n_max=4 tested for F1 (multipole containment) only; F2 saturation is n_max=2 specific.
- **Non-relativistic for QPT main test**. Relativistic Tier 2 tested separately and structurally compatible.
- **JW fermion-to-qubit only**. BK / parity / ternary tree encodings not tested.

## Multi-day follow-ons named

| sprint | scope | wall-clock | grade |
|---|---|---|---|
| S2-v2 closed-form | Derive χ_k = f(R^L rank, spin) explicitly | 3-5 days | analytical |
| S3 Berezin transport | Paper 38 L4 → math.OA companion paper | ~2 weeks | math-grade |
| Paper 19 full negative control | Bond rank on balanced-coupled across library | 1 day | empirical |
| n_max=3 bond rank | Confirm χ=2 boundary preserved at larger basis | 1 day | empirical |
| BK / parity bond rank | Effect of fermion mapping on χ profile | 1 day | empirical |

## Where the project stands at end-of-day

The day started with PI vibing the question "what else is leaving DF/multipole-on-the-table?" The lit-scan flagged Cholesky as the easy test, QPT as the next, MPO bond-rank as the deepest. All three landed positive.

The Connes-vS / MPO bridge was scoped as multi-month → 3 sprints → S1 done in main session. S2 partial converges to a structural F3 connection. S3 (math.OA companion paper) is the remaining sprint of the original 3.

**Net: the PI's bet — "scoping turns multimonth into a couple sprints" — held twice in one day.**

## Files produced

Sprint memos:
- `debug/sprint_df_multipole_lift_memo.md` (revised with Cholesky cross-ref)
- `debug/sprint_cholesky_multipole_memo.md`
- `debug/sprint_qpt_hopf_stacking_memo.md`
- `debug/sprint_qpt_followons_memo.md`
- `debug/sprint_s1_bond_rank_memo.md`
- `debug/sprint_s2_partial_memo.md`
- `debug/sprint_s2_closed_form_attempt_memo.md`
- `debug/sprint_day_consolidation_memo.md` (this file)

Scoping / lit-scan agent outputs:
- `debug/lit_search_algebraic_methods_compatibility_2026_06_05.md`
- `debug/sprint_connes_vs_mpo_scoping_memo.md`

Test scripts:
- `debug/df_vs_multipole_rank_test.py`, `debug/df_radial_density_count.py`, `debug/df_factor_structure.py`, `debug/df_per_L_rank_decomposition.py`, `debug/df_verify_decomposition.py`, `debug/df_cross_L_dependences.py`
- `debug/cholesky_vs_multipole_test.py`, `debug/cholesky_atomic_nmax_sweep.py`
- `debug/qpt_hopf_stacking_test.py`, `debug/qpt_joint_block_dim.py`, `debug/qpt_relativistic_test.py`, `debug/qpt_cost_table.py`
- `debug/sprint_s1_bond_rank.py`, `debug/sprint_s1_negative_control.py`
- `debug/sprint_s2_decode_profile.py`, `debug/sprint_s2_closed_form.py`

Plus ~15 JSON data files under `debug/data/`.

## Recommended next-session start

If picking back up tomorrow, the cleanest entry point is **S2-v2 closed-form derivation**: with F3 connection identified, the analytical path to chi_k = f(R^L rank, spin combinatorics) is now navigable. 3-5 focused days → Theorem 3.2.A.unified → Paper 14 §sec:mpo_bond_rank ready to ship.

Alternative entry points:
- **Apply Paper 14 §sec:hopf_tapering QPT addition** (~30 min, quick win).
- **Apply Paper 14 §sec:mpo_bond_rank empirical version** (~30 min, defer closed-form).
- **Run Paper 19 full library negative control** (~1 day, strengthens S1 result).

The vibe-physics arc is at a clean pause point.
